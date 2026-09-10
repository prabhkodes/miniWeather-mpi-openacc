!> @brief Physics driver for the 2-D compressible Euler solver.
!! @details Owns prognostic/reference state, initializes fields, advances the
!!          solution via dimensionally split Runge-Kutta, and reports conserved
!!          quantities. Routines are instrumented with timers for MPI-wide
!!          profiling and support CPU or OpenACC offload paths.
module module_physics
  use mpi
#ifdef USE_OPENACC
  use openacc
  use module_nvtx
#endif
  use calculation_types, only : wp
  use physical_constants
  use physical_parameters
  use parallel_parameters
  use indexing
  use legendre_quadrature
  use dimensions
  use iodir
  use module_types
  use parallel_timer

  implicit none

  private

  public :: init
  public :: finalize
  public :: rungekutta
  public :: total_mass_energy

  real(wp), public :: dt       !< Stable time step selected from CFL condition.
  real(wp) :: dx, dz           !< Grid spacing in x and z (local but consistent with global nx).

  type(atmospheric_state), public :: oldstat !< Current prognostic state (density, momentum, rho*theta).
  type(atmospheric_state), public :: newstat !< Next-step prognostic state.
  type(atmospheric_tendency), public :: tend !< Tendency workspace for flux divergences and sources.
  type(atmospheric_flux), public :: flux     !< Flux workspace reused across sweeps.
  type(reference_state), public :: ref       !< Hydrostatic background state.

  contains

  !> @brief Initialize grid spacing, states, reference profiles, and dt.
  !! @param[out] etime Simulation time (reset to 0).
  !! @param[out] output_counter Accumulator for output cadence (reset to 0).
  !! @param[out] dt Chosen stable time step based on min(dx,dz), `max_speed`, and `cfl`.
  !! @details Allocates state/flux/tendency arrays, computes quadrature-weighted
  !!          initial perturbation fields using `thermal`, builds hydrostatic
  !!          reference columns, and logs configuration on rank 0.
  subroutine init(etime,output_counter,dt)
    implicit none

    type(CTimer) :: pll_timer 
    real(wp), intent(out) :: etime, output_counter, dt
    integer :: i, k, ii, kk, ierr
    real(wp) :: x, z, r, u, w, t, hr, ht
    integer :: i_global

    call pll_timer%start("INIT")

    ! Use global nx for dx calculation to maintain consistent grid spacing
    dx = xlen / nx_global
    dz = zlen / nz

    call oldstat%new_state( )
    call newstat%new_state( )
    call flux%new_flux( )
    call tend%new_tendency( )
    call ref%new_ref( )

    call init_halo_buffers()

    dt = min(dx,dz) / max_speed * cfl
    etime = 0.0_wp
    output_counter = 0.0_wp

    if (myrank == 0) then
      write(stdout,*) 'INITIALIZING MODEL STATUS.'
      write(stdout,*) 'nx_global  : ', nx_global
      write(stdout,*) 'nx_local   : ', nx
      write(stdout,*) 'nz         : ', nz
      write(stdout,*) 'dx         : ', dx
      write(stdout,*) 'dz         : ', dz
      write(stdout,*) 'dt         : ', dt
      write(stdout,*) 'final time : ', sim_time
    end if

    call oldstat%set_state(0.0_wp)

    do k = 1-hs, nz+hs
      do i = 1-hs, nx+hs
        ! Convert local index to global index
        i_global = i_start_global + i - 1
        do kk = 1, nqpoints
          do ii = 1, nqpoints
            ! Use global index for x coordinate calculation
            x = (i_global - 0.5_wp) * dx + (qpoints(ii)-0.5_wp)*dx
            z = (k_beg-1 + k-0.5_wp) * dz + (qpoints(kk)-0.5_wp)*dz
            call thermal(x,z,r,u,w,t,hr,ht)
            oldstat%dens(i,k) = oldstat%dens(i,k) + &
                       r * qweights(ii)*qweights(kk)
            oldstat%umom(i,k) = oldstat%umom(i,k) + &
                      (r+hr)*u * qweights(ii)*qweights(kk)
            oldstat%wmom(i,k) = oldstat%wmom(i,k) + &
                      (r+hr)*w * qweights(ii)*qweights(kk)
            oldstat%rhot(i,k) = oldstat%rhot(i,k) + &
                    ( (r+hr)*(t+ht) - hr*ht ) * qweights(ii)*qweights(kk)
          end do
        end do
      end do
    end do
    newstat = oldstat
    ref%density(:) = 0.0_wp
    ref%denstheta(:) = 0.0_wp
    do k = 1-hs, nz+hs
      do kk = 1, nqpoints
        z = (k_beg-1 + k-0.5_wp) * dz + (qpoints(kk)-0.5_wp)*dz
        call thermal(0.0_wp,z,r,u,w,t,hr,ht)
        ref%density(k) = ref%density(k) + hr * qweights(kk)
        ref%denstheta(k) = ref%denstheta(k) + hr*ht * qweights(kk)
      end do
    end do
    do k = 1, nz+1
      z = (k_beg-1 + k-1) * dz
      call thermal(0.0_wp,z,r,u,w,t,hr,ht)
      ref%idens(k) = hr
      ref%idenstheta(k) = hr*ht
      ref%pressure(k) = c0*(hr*ht)**cdocv
    end do
    if (myrank == 0) then
      write(stdout,*) 'MODEL STATUS INITIALIZED.'
    end if
  end subroutine init

  !> @brief Third-order Runge-Kutta driver with dimensional splitting.
  !! @param[inout] s0 Baseline state used for stages (often `oldstat`).
  !! @param[inout] s1 Scratch/target state for intermediate stages.
  !! @param[inout] fl Flux workspace shared across stages.
  !! @param[inout] tend Tendency workspace shared across stages.
  !! @param[in] dt Time step size for this full RK update.
  !! @details Alternates the order of x- then z-sweeps each call to achieve
  !!          overall second-order accuracy in time when splitting.
  subroutine rungekutta(s0,s1,fl,tend,dt)
    implicit none

    type(CTimer) :: pll_timer 
    type(atmospheric_state), intent(inout) :: s0
    type(atmospheric_state), intent(inout) :: s1
    type(atmospheric_flux), intent(inout) :: fl
    type(atmospheric_tendency), intent(inout) :: tend
    real(wp), intent(in) :: dt
    real(wp) :: dt1, dt2, dt3
    logical, save :: dimswitch = .true.

    call pll_timer%start("Computation: rungekutta")

#ifdef USE_OPENACC
    call nvtx_push("Runge-Kutta Integration")
#endif

    dt1 = dt/1.0_wp
    dt2 = dt/2.0_wp
    dt3 = dt/3.0_wp
    if ( dimswitch ) then
      call step(s0, s0, s1, dt3, DIR_X, fl, tend)
      call step(s0, s1, s1, dt2, DIR_X, fl, tend)
      call step(s0, s1, s0, dt1, DIR_X, fl, tend)
      call step(s0, s0, s1, dt3, DIR_Z, fl, tend)
      call step(s0, s1, s1, dt2, DIR_Z, fl, tend)
      call step(s0, s1, s0, dt1, DIR_Z, fl, tend)
    else
      call step(s0, s0, s1, dt3, DIR_Z, fl, tend)
      call step(s0, s1, s1, dt2, DIR_Z, fl, tend)
      call step(s0, s1, s0, dt1, DIR_Z, fl, tend)
      call step(s0, s0, s1, dt3, DIR_X, fl, tend)
      call step(s0, s1, s1, dt2, DIR_X, fl, tend)
      call step(s0, s1, s0, dt1, DIR_X, fl, tend)
    end if
    dimswitch = .not. dimswitch

#ifdef USE_OPENACC
    call nvtx_pop()
#endif
  end subroutine rungekutta

  ! Semi-discretized step in time:
  ! s2 = s0 + dt * rhs(s1)
  !> @brief Perform one directional sweep update.
  !! @param[in] s0 Baseline state (input).
  !! @param[inout] s1 State to differentiate (also receives halo updates).
  !! @param[inout] s2 Destination state after applying `tend`.
  !! @param[in] dt Time step for this sub-stage.
  !! @param[in] dir Direction flag (`DIR_X` or `DIR_Z`) selecting flux routine.
  !! @param[inout] fl Flux workspace.
  !! @param[inout] tend Tendency workspace.
  subroutine step(s0, s1, s2, dt, dir, fl, tend)
    implicit none

    type(CTimer) :: pll_timer 
    type(atmospheric_state), intent(in) :: s0
    type(atmospheric_state), intent(inout) :: s1
    type(atmospheric_state), intent(inout) :: s2
    type(atmospheric_flux), intent(inout) :: fl
    type(atmospheric_tendency), intent(inout) :: tend
    real(wp), intent(in) :: dt
    integer, intent(in) :: dir

    call pll_timer%start("Computation: step")

    if (dir == DIR_X) then
      call tend%xtend(fl,ref,s1,dx,dt)
    else if (dir == DIR_Z) then
      call tend%ztend(fl,ref,s1,dz,dt)
    end if
    call s2%update(s0,tend,dt)
  end subroutine step

  !> @brief Compute initial perturbation and hydrostatic background at a point.
  !! @param[in] x Horizontal coordinate [m].
  !! @param[in] z Vertical coordinate [m].
  !! @param[out] r Density perturbation (rho').
  !! @param[out] u x-velocity perturbation.
  !! @param[out] w z-velocity perturbation.
  !! @param[out] t Potential temperature perturbation.
  !! @param[out] hr Hydrostatic background density.
  !! @param[out] ht Hydrostatic background potential temperature.
  subroutine thermal(x,z,r,u,w,t,hr,ht)
    implicit none
    type(CTimer) :: pll_timer 
    real(wp), intent(in) :: x, z
    real(wp), intent(out) :: r, u, w, t
    real(wp), intent(out) :: hr, ht
    call pll_timer%start("Computation: thermal")
    call hydrostatic_const_theta(z,hr,ht)
    r = 0.0_wp
    t = 0.0_wp
    u = 0.0_wp
    w = 0.0_wp
    t = t + ellipse(x,z,3.0_wp,hxlen,p1,p1,p1)
  end subroutine thermal

  !> @brief Hydrostatic background with constant potential temperature.
  !! @param[in] z Vertical coordinate [m].
  !! @param[out] r Hydrostatic density at height z.
  !! @param[out] t Hydrostatic potential temperature (constant here).
  !! @details Integrates hydrostatic balance analytically using Exner function
  !!          to obtain density profile consistent with constant theta.
  subroutine hydrostatic_const_theta(z,r,t)
    implicit none

    type(CTimer) :: pll_timer 
    real(wp), intent(in) :: z
    real(wp), intent(out) :: r, t
    real(wp) :: p,exner,rt

    call pll_timer%start("Computation: hydrostatic_const_theta")

    t = theta0
    exner = exner0 - grav * z / (cp * theta0)
    p = p0 * exner**(cp/rd)
    rt = (p / c0)**cvocd
    r = rt / t
  end subroutine hydrostatic_const_theta

  !> @brief Smooth cosine-squared elliptical perturbation.
  !! @param[in] x Horizontal coordinate [m].
  !! @param[in] z Vertical coordinate [m].
  !! @param[in] amp Amplitude of the perturbation.
  !! @param[in] x0,z0 Center of ellipse.
  !! @param[in] x1,z1 Horizontal/vertical radii.
  !! @return val Perturbation value at (x,z).
  elemental function ellipse(x,z,amp,x0,z0,x1,z1) result(val)
    implicit none
    
    real(wp), intent(in) :: x, z
    real(wp), intent(in) :: amp
    real(wp), intent(in) :: x0, z0
    real(wp), intent(in) :: x1, z1
    real(wp) :: val
    real(wp) :: dist

    dist = sqrt( ((x-x0)/x1)**2 + ((z-z0)/z1)**2 ) * hpi
    if (dist <= hpi) then
      val = amp * cos(dist)**2
    else
      val = 0.0_wp
    end if
  end function ellipse

  !> @brief Release allocated state/flux/tendency/reference storage.
  subroutine finalize()
    implicit none
    call oldstat%del_state( )
    call newstat%del_state( )
    call flux%del_flux( )
    call tend%del_tendency( )
    call ref%del_ref( )
    call free_halo_buffers()
  end subroutine finalize

  !> @brief Compute global total mass and total energy.
  !! @param[out] mass Integrated mass over the domain.
  !! @param[out] te Integrated total energy (kinetic + internal).
  !! @details Accumulates local contributions (optionally on GPU), then performs
  !!          MPI all-reduce to obtain global totals. Uses reference state to
  !!          reconstruct full density/potential temperature.
  subroutine total_mass_energy(mass,te)
    implicit none
    
    type(CTimer) :: pll_timer 
    type(CTimer) :: mpi_timer
    real(wp), intent(out) :: mass, te
    integer :: i, k, ierr
    real(wp) :: r, u, w, th, p, t, ke, ie
    real(wp) :: local_mass, local_te

    call pll_timer%start("Computation: total_mass_energy")

#ifdef USE_OPENACC
    call nvtx_push("Conservation Check")
#endif

    local_mass = 0.0_wp
    local_te = 0.0_wp
#ifdef USE_OPENACC
    if (acc_is_present(oldstat%mem)) then
      !$acc parallel loop collapse(2) present(oldstat%mem, ref%density, ref%denstheta) reduction(+:local_mass, local_te)
      do k = 1, nz
        do i = 1, nx
          r = oldstat%dens(i,k) + ref%density(k)
          u = oldstat%umom(i,k) / r
          w = oldstat%wmom(i,k) / r
          th = (oldstat%rhot(i,k) + ref%denstheta(k) ) / r
          p = c0*(r*th)**cdocv
          t = th / (p0/p)**rdocp
          ke = r*(u*u+w*w)
          ie = r*cv*t
          local_mass = local_mass + r * dx * dz
          local_te = local_te + (ke + r*cv*t) * dx * dz
        end do
      end do
    else
#endif
      do k = 1, nz
        do i = 1, nx
          r = oldstat%dens(i,k) + ref%density(k)
          u = oldstat%umom(i,k) / r
          w = oldstat%wmom(i,k) / r
          th = (oldstat%rhot(i,k) + ref%denstheta(k) ) / r
          p = c0*(r*th)**cdocv
          t = th / (p0/p)**rdocp
          ke = r*(u*u+w*w)
          ie = r*cv*t
          local_mass = local_mass + r * dx * dz
          local_te = local_te + (ke + r*cv*t) * dx * dz
        end do
      end do
#ifdef USE_OPENACC
    end if
#endif

    call mpi_timer%start("MPI: Communication")
    ! Reduce to global values across all MPI processes
    call MPI_Allreduce(local_mass, mass, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_te, te, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call mpi_timer%stop()

#ifdef USE_OPENACC
    call nvtx_pop()
#endif

  end subroutine total_mass_energy

end module module_physics
