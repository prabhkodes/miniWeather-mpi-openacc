!> @brief Core data structures and halo exchange for prognostic fields.
!! @details Defines reference/prognostic/flux/tendency containers, allocation
!!          helpers, halo exchange in x/z, and stencil-based flux divergence
!!          routines. Supports MPI domain decomposition with optional OpenMP or
!!          OpenACC offload.
module module_types
  use mpi
#ifdef USE_OPENACC
  use openacc
  use module_nvtx
#endif
  use calculation_types
  use physical_constants
  use physical_parameters
  use parallel_parameters
  use indexing
  use legendre_quadrature
  use dimensions
  use iodir
  use parallel_timer

  implicit none

  private

  public :: reference_state
  public :: atmospheric_state
  public :: atmospheric_flux
  public :: atmospheric_tendency

  public :: assignment(=)
  public :: init_halo_buffers
  public :: free_halo_buffers

  real(wp), allocatable, save :: send_left(:,:,:), send_right(:,:,:)
  real(wp), allocatable, save :: recv_left(:,:,:), recv_right(:,:,:)
  integer, save :: halo_size = 0
  logical, save :: halos_allocated = .false.

  !> @brief Hydrostatic reference profiles and interface values.
  type reference_state
	real(wp), allocatable, dimension(:) :: density     !< Cell-centered hydrostatic density.
	real(wp), allocatable, dimension(:) :: denstheta   !< Cell-centered rho*theta for hydrostatic state.
	real(wp), allocatable, dimension(:) :: idens       !< Interface hydrostatic density.
	real(wp), allocatable, dimension(:) :: idenstheta  !< Interface hydrostatic rho*theta.
	real(wp), allocatable, dimension(:) :: pressure    !< Interface hydrostatic pressure.
	contains
	procedure, public :: new_ref
	procedure, public :: del_ref
  end type reference_state

  !> @brief Prognostic state (density, momentum, rho*theta) with halo padding.
  type atmospheric_state
	real(wp), pointer, dimension(:,:,:) :: mem => null( ) !< Backing contiguous storage.
	real(wp), pointer, dimension(:,:) :: dens             !< Density perturbation (rho').
	real(wp), pointer, dimension(:,:) :: umom             !< x-momentum.
	real(wp), pointer, dimension(:,:) :: wmom             !< z-momentum.
	real(wp), pointer, dimension(:,:) :: rhot             !< rho*theta perturbation.
	contains
	procedure, public :: new_state
	procedure, public :: set_state
	procedure, public :: del_state
	procedure, public :: update
	procedure, public :: exchange_halo_x
	procedure, public :: exchange_halo_z
  end type atmospheric_state

  !> @brief Flux container used during directional sweeps.
  type atmospheric_flux
	real(wp), pointer, dimension(:,:,:) :: mem => null( ) !< Backing contiguous storage.
	real(wp), pointer, dimension(:,:) :: dens             !< Flux of density.
	real(wp), pointer, dimension(:,:) :: umom             !< Flux of x-momentum.
	real(wp), pointer, dimension(:,:) :: wmom             !< Flux of z-momentum.
	real(wp), pointer, dimension(:,:) :: rhot             !< Flux of rho*theta.
	contains
	procedure, public :: new_flux
	procedure, public :: set_flux
	procedure, public :: del_flux
  end type atmospheric_flux

  !> @brief Right-hand-side tendencies from flux divergences and sources.
  type atmospheric_tendency
	real(wp), pointer, dimension(:,:,:) :: mem => null( ) !< Backing contiguous storage.
	real(wp), pointer, dimension(:,:) :: dens             !< d(rho)/dt tendency.
	real(wp), pointer, dimension(:,:) :: umom             !< d(ru)/dt tendency.
	real(wp), pointer, dimension(:,:) :: wmom             !< d(rw)/dt tendency.
	real(wp), pointer, dimension(:,:) :: rhot             !< d(rho*theta)/dt tendency.
	contains
	procedure, public :: new_tendency
	procedure, public :: set_tendency
	procedure, public :: del_tendency
	procedure, public :: xtend
	procedure, public :: ztend
  end type atmospheric_tendency

  interface assignment(=)
	module procedure state_equal_to_state
  end interface assignment(=)

  contains

  subroutine init_halo_buffers()
	implicit none

	if (halos_allocated) return

	halo_size = hs * nz * NVARS
	allocate(send_left(hs, nz, NVARS), send_right(hs, nz, NVARS))
	allocate(recv_left(hs, nz, NVARS), recv_right(hs, nz, NVARS))
#ifdef USE_OPENACC
	! Keep halo buffers resident on the device to avoid per-call staging
	!$acc enter data create(send_left, send_right, recv_left, recv_right)
#endif
	halos_allocated = .true.
  end subroutine init_halo_buffers

  subroutine free_halo_buffers()
	implicit none
	if (.not. halos_allocated) return
#ifdef USE_OPENACC
	if (acc_is_present(send_left)) then
	  !$acc exit data delete(send_left, send_right, recv_left, recv_right)
	end if
#endif
	if (allocated(send_left)) deallocate(send_left)
	if (allocated(send_right)) deallocate(send_right)
	if (allocated(recv_left)) deallocate(recv_left)
	if (allocated(recv_right)) deallocate(recv_right)
	halo_size = 0
	halos_allocated = .false.
  end subroutine free_halo_buffers

  !> @brief Allocate state storage with halos and set component pointers.
  subroutine new_state(atmo)
	implicit none
	class(atmospheric_state), intent(inout) :: atmo
	if ( associated(atmo%mem) ) deallocate(atmo%mem)
	allocate(atmo%mem(1-hs:nx+hs, 1-hs:nz+hs, NVARS))
	atmo%dens(1-hs:,1-hs:) => atmo%mem(:,:,I_DENS)
	atmo%umom(1-hs:,1-hs:) => atmo%mem(:,:,I_UMOM)
	atmo%wmom(1-hs:,1-hs:) => atmo%mem(:,:,I_WMOM)
	atmo%rhot(1-hs:,1-hs:) => atmo%mem(:,:,I_RHOT)
  end subroutine new_state

  !> @brief Fill entire state (including halos) with a constant value.
  subroutine set_state(atmo, xval)
	implicit none
	class(atmospheric_state), intent(inout) :: atmo
	real(wp), intent(in) :: xval
	if ( .not. associated(atmo%mem) ) then
	  write(stderr,*) 'NOT ALLOCATED STATE ERROR AT LINE ', __LINE__
	  stop
	end if
	atmo%mem(:,:,:) = xval
  end subroutine set_state

  !> @brief Deallocate state storage and nullify component pointers.
  subroutine del_state(atmo)
	implicit none
	class(atmospheric_state), intent(inout) :: atmo
	if ( associated(atmo%mem) ) deallocate(atmo%mem)
	nullify(atmo%dens)
	nullify(atmo%umom)
	nullify(atmo%wmom)
	nullify(atmo%rhot)
  end subroutine del_state

  !> @brief Explicit Euler update: s2 = s0 + dt * tendency.
  !! @param[inout] s2 Destination state updated in place.
  !! @param[in] s0 Baseline state.
  !! @param[in] tend Right-hand-side tendencies.
  !! @param[in] dt Time step for this update.
  subroutine update(s2,s0,tend,dt)
	implicit none
	class(atmospheric_state), intent(inout) :: s2
	class(atmospheric_state), intent(in) :: s0
	class(atmospheric_tendency), intent(in) :: tend
	real(wp), intent(in) :: dt
	integer :: ll, k, i
	
#ifdef USE_OPENACC
	  !$acc parallel loop collapse(3) present(s0%mem, tend%mem, s2%mem)
	  do ll = 1, NVARS
		do k = 1, nz
		  do i = 1, nx
			s2%mem(i,k,ll) = s0%mem(i,k,ll) + dt * tend%mem(i,k,ll)
		  end do
		end do
	  end do
	  !$acc end parallel loop
#else
	  do ll = 1, NVARS
		do k = 1, nz
		  do i = 1, nx
			s2%mem(i,k,ll) = s0%mem(i,k,ll) + dt * tend%mem(i,k,ll)
		  end do
		end do
	  end do
#endif
	
  end subroutine update


  !> @brief Compute x-direction fluxes and tendencies (with hyper-viscosity).
  !! @param[inout] tendency Accumulates d/dt contributions for all variables.
  !! @param[inout] flux Workspace to store face fluxes.
  !! @param[in] ref Hydrostatic reference used to recover full state.
  !! @param[inout] atmostat State being differentiated (halos updated inside).
  !! @param[in] dx Cell width in x.
  !! @param[in] dt Time step (used for hyper-viscosity coefficient).
  subroutine xtend(tendency,flux,ref,atmostat,dx,dt)
	implicit none
	class(atmospheric_tendency), intent(inout) :: tendency
	class(atmospheric_flux), intent(inout) :: flux
	class(reference_state), intent(in) :: ref
	class(atmospheric_state), intent(inout) :: atmostat
	real(wp), intent(in) :: dx, dt
	integer :: i, k, ll, s
	real(wp) :: r, u, w, t, p, hv_coef

	real(wp), dimension(STEN_SIZE) :: stencil
	real(wp), dimension(NVARS) :: d3_vals, vals

#ifdef USE_OPENACC
	call nvtx_push("X-Direction Tendency")
#endif

	call atmostat%exchange_halo_x( )

	hv_coef = -hv_beta * dx / (16.0_wp*dt)

#ifdef USE_OPENACC
	  !$acc parallel loop collapse(2) present(atmostat%mem, ref%density, ref%denstheta, flux%mem) &
	  !$acc   private(stencil, vals, d3_vals, r, u, w, t, p)
	  do k = 1, nz
		do i = 1, nx+1
		  do ll = 1, NVARS
			do s = 1, STEN_SIZE
			  stencil(s) = atmostat%mem(i-hs-1+s,k,ll)
			end do
			vals(ll) = - 1.0_wp * stencil(1)/12.0_wp &
						+ 7.0_wp * stencil(2)/12.0_wp &
						+ 7.0_wp * stencil(3)/12.0_wp &
						- 1.0_wp * stencil(4)/12.0_wp
			d3_vals(ll) = - 1.0_wp * stencil(1) &
						  + 3.0_wp * stencil(2) &
						  - 3.0_wp * stencil(3) &
						  + 1.0_wp * stencil(4)
		  end do
		  r = vals(I_DENS) + ref%density(k)
		  u = vals(I_UMOM) / r
		  w = vals(I_WMOM) / r
		  t = ( vals(I_RHOT) + ref%denstheta(k) ) / r
		  p = c0*(r*t)**cdocv
		  flux%mem(i,k,I_DENS) = r*u - hv_coef*d3_vals(I_DENS)
		  flux%mem(i,k,I_UMOM) = r*u*u+p - hv_coef*d3_vals(I_UMOM)
		  flux%mem(i,k,I_WMOM) = r*u*w - hv_coef*d3_vals(I_WMOM)
		  flux%mem(i,k,I_RHOT) = r*u*t - hv_coef*d3_vals(I_RHOT)
		end do
	  end do
	  !$acc end parallel loop

	  !$acc parallel loop collapse(3) present(tendency%mem, flux%mem)
	  do ll = 1, NVARS
		do k = 1, nz
		  do i = 1, nx
			tendency%mem(i,k,ll) = &
				-( flux%mem(i+1,k,ll) - flux%mem(i,k,ll) ) / dx
		  end do
		end do
	  end do
	  !$acc end parallel loop
#else
	  !$omp parallel do default(shared) &
	  !$omp private(i, k, ll, s, stencil, vals, d3_vals, r, u, w, t, p)
	  do k = 1, nz
		do i = 1, nx+1
		  do ll = 1, NVARS
			do s = 1, STEN_SIZE
			  stencil(s) = atmostat%mem(i-hs-1+s,k,ll)
			end do
			vals(ll) = - 1.0_wp * stencil(1)/12.0_wp &
						+ 7.0_wp * stencil(2)/12.0_wp &
						+ 7.0_wp * stencil(3)/12.0_wp &
						- 1.0_wp * stencil(4)/12.0_wp
			d3_vals(ll) = - 1.0_wp * stencil(1) &
						  + 3.0_wp * stencil(2) &
						  - 3.0_wp * stencil(3) &
						  + 1.0_wp * stencil(4)
		  end do
		  r = vals(I_DENS) + ref%density(k)
		  u = vals(I_UMOM) / r
		  w = vals(I_WMOM) / r
		  t = ( vals(I_RHOT) + ref%denstheta(k) ) / r
		  p = c0*(r*t)**cdocv
		  flux%dens(i,k) = r*u - hv_coef*d3_vals(I_DENS)
		  flux%umom(i,k) = r*u*u+p - hv_coef*d3_vals(I_UMOM)
		  flux%wmom(i,k) = r*u*w - hv_coef*d3_vals(I_WMOM)
		  flux%rhot(i,k) = r*u*t - hv_coef*d3_vals(I_RHOT)
		end do
	  end do
	  !$omp end parallel do

	  !$omp parallel do default(shared) private(i, k, ll) collapse(2)
	  do ll = 1, NVARS
		do k = 1, nz
		  do i = 1, nx
			tendency%mem(i,k,ll) = &
				-( flux%mem(i+1,k,ll) - flux%mem(i,k,ll) ) / dx
		  end do
		end do
	  end do
	  !$omp end parallel do
#endif

#ifdef USE_OPENACC
	call nvtx_pop()
#endif

  end subroutine xtend


  !> @brief Compute z-direction fluxes/tendencies and add gravity source.
  !! @param[inout] tendency Accumulates d/dt contributions for all variables.
  !! @param[inout] flux Workspace to store face fluxes.
  !! @param[in] ref Hydrostatic reference used to recover full state and pressure perturbations.
  !! @param[inout] atmostat State being differentiated (halos updated inside).
  !! @param[in] dz Cell height in z.
  !! @param[in] dt Time step (used for hyper-viscosity coefficient).
  subroutine ztend(tendency,flux,ref,atmostat,dz,dt)
	implicit none
	class(atmospheric_tendency), intent(inout) :: tendency
	class(atmospheric_flux), intent(inout) :: flux
	class(reference_state), intent(in) :: ref
	class(atmospheric_state), intent(inout) :: atmostat
	real(wp), intent(in) :: dz, dt
	integer :: i, k, ll, s
	real(wp) :: r, u, w, t, p, hv_coef
	real(wp), dimension(STEN_SIZE) :: stencil
	real(wp), dimension(NVARS) :: d3_vals, vals

#ifdef USE_OPENACC
	call nvtx_push("Z-Direction Tendency")
#endif

	call atmostat%exchange_halo_z(ref)

	hv_coef = -hv_beta * dz / (16.0_wp*dt)

#ifdef USE_OPENACC
	  !$acc parallel loop collapse(2) present(atmostat%mem, ref%idens, ref%idenstheta, ref%pressure, flux%mem) &
	  !$acc   private(stencil, vals, d3_vals, r, u, w, t, p)
	  do k = 1, nz+1
		do i = 1, nx
		  do ll = 1, NVARS
			do s = 1, STEN_SIZE
			  stencil(s) = atmostat%mem(i,k-hs-1+s,ll)
			end do
			vals(ll) = - 1.0_wp * stencil(1)/12.0_wp &
					  + 7.0_wp * stencil(2)/12.0_wp &
					  + 7.0_wp * stencil(3)/12.0_wp &
					  - 1.0_wp * stencil(4)/12.0_wp
			d3_vals(ll) = - 1.0_wp * stencil(1) &
						  + 3.0_wp * stencil(2) &
						  - 3.0_wp * stencil(3) &
						  + 1.0_wp * stencil(4)
		  end do
		  r = vals(I_DENS) + ref%idens(k)
		  u = vals(I_UMOM) / r
		  w = vals(I_WMOM) / r
		  t = ( vals(I_RHOT) + ref%idenstheta(k) ) / r
		  p = c0*(r*t)**cdocv - ref%pressure(k)
		  
		  if (k == 1 .or. k == nz+1) then
			w = 0.0_wp
			d3_vals(I_DENS) = 0.0_wp
		  end if
		  
		  flux%mem(i,k,I_DENS) = r*w - hv_coef*d3_vals(I_DENS)
		  flux%mem(i,k,I_UMOM) = r*w*u - hv_coef*d3_vals(I_UMOM)
		  flux%mem(i,k,I_WMOM) = r*w*w+p - hv_coef*d3_vals(I_WMOM)
		  flux%mem(i,k,I_RHOT) = r*w*t - hv_coef*d3_vals(I_RHOT)
		end do
	  end do
	  !$acc end parallel loop

	  !$acc parallel loop collapse(3) present(tendency%mem, flux%mem)
	  do ll = 1, NVARS
		do k = 1, nz
		  do i = 1, nx
			tendency%mem(i,k,ll) = &
				-( flux%mem(i,k+1,ll) - flux%mem(i,k,ll) ) / dz
		  end do
		end do
	  end do
	  !$acc end parallel loop
	  
	  !$acc parallel loop collapse(2) present(tendency%mem, atmostat%mem)
	  do k = 1, nz
		do i = 1, nx
		  tendency%mem(i,k,I_WMOM) = tendency%mem(i,k,I_WMOM) - atmostat%mem(i,k,I_DENS)*grav
		end do
	  end do
	  !$acc end parallel loop

#else
	!$omp parallel do default(shared) &
	!$omp private(i, k, ll, s, stencil, vals, d3_vals, r, u, w, t, p)
	do k = 1, nz+1
	  do i = 1, nx
		do ll = 1, NVARS
		  do s = 1, STEN_SIZE
			stencil(s) = atmostat%mem(i,k-hs-1+s,ll)
		  end do
		  vals(ll) = - 1.0_wp * stencil(1)/12.0_wp &
					 + 7.0_wp * stencil(2)/12.0_wp &
					 + 7.0_wp * stencil(3)/12.0_wp &
					 - 1.0_wp * stencil(4)/12.0_wp
		  d3_vals(ll) = - 1.0_wp * stencil(1) &
						+ 3.0_wp * stencil(2) &
						- 3.0_wp * stencil(3) &
						+ 1.0_wp * stencil(4)
		end do
		r = vals(I_DENS) + ref%idens(k)
		u = vals(I_UMOM) / r
		w = vals(I_WMOM) / r
		t = ( vals(I_RHOT) + ref%idenstheta(k) ) / r
		p = c0*(r*t)**cdocv - ref%pressure(k)
		
		! This IF statement is thread-safe because it relies only on 'k'
		if (k == 1 .or. k == nz+1) then
		  w = 0.0_wp
		  d3_vals(I_DENS) = 0.0_wp
		end if
		
		flux%dens(i,k) = r*w - hv_coef*d3_vals(I_DENS)
		flux%umom(i,k) = r*w*u - hv_coef*d3_vals(I_UMOM)
		flux%wmom(i,k) = r*w*w+p - hv_coef*d3_vals(I_WMOM)
		flux%rhot(i,k) = r*w*t - hv_coef*d3_vals(I_RHOT)
	  end do
	end do
	!$omp end parallel do

	!$omp parallel do default(shared) private(i, k, ll) collapse(2)
	do ll = 1, NVARS
	  do k = 1, nz
		do i = 1, nx
		  tendency%mem(i,k,ll) = &
			  -( flux%mem(i,k+1,ll) - flux%mem(i,k,ll) ) / dz
		  
		  ! Adding gravity source term
		  if (ll == I_WMOM) then
			tendency%wmom(i,k) = tendency%wmom(i,k) - atmostat%dens(i,k)*grav
		  end if
		end do
	  end do
	  end do
	  !$omp end parallel do
#endif

#ifdef USE_OPENACC
	call nvtx_pop()
#endif

  end subroutine ztend
  !> @brief Exchange halo cells in x-direction with periodic boundaries.
  !! @param[inout] s State whose halos are updated.
  !! @details For serial runs, halos wrap locally. For MPI runs, exchanges
  !!          boundary layers with left/right ranks using non-blocking sends.
  subroutine exchange_halo_x(s)
    use mpi
    implicit none
    class(atmospheric_state), intent(inout) :: s
    integer :: k, ll, ierr
    integer :: send_size
    integer :: req(4), status(MPI_STATUS_SIZE, 4)
    type(CTimer) :: mpi_timer

    if (.not. halos_allocated) call init_halo_buffers()
    send_size = halo_size

#ifdef USE_OPENACC
    call nvtx_push("Halo Exchange X")
#endif

    if (nprocs == 1) then
      ! Serial case: local periodic boundary conditions
#ifdef USE_OPENACC
      !$acc parallel loop collapse(2) present(s%mem)
      do ll = 1, NVARS
        do k = 1, nz
          s%mem(-1,k,ll)   = s%mem(nx-1,k,ll)
          s%mem(0,k,ll)    = s%mem(nx,k,ll)
          s%mem(nx+1,k,ll) = s%mem(1,k,ll)
          s%mem(nx+2,k,ll) = s%mem(2,k,ll)
        end do
      end do
#else
      do ll = 1, NVARS
        do k = 1, nz
          s%mem(-1,k,ll)   = s%mem(nx-1,k,ll)
          s%mem(0,k,ll)    = s%mem(nx,k,ll)
          s%mem(nx+1,k,ll) = s%mem(1,k,ll)
          s%mem(nx+2,k,ll) = s%mem(2,k,ll)
        end do
      end do
#endif
    else
      ! Parallel case: MPI communication

#ifdef USE_OPENACC
      if (acc_is_present(s%mem)) then
        ! Pack halos on device to avoid host staging
        !$acc parallel loop collapse(2) present(send_left, send_right, s%mem)
        do ll = 1, NVARS
          do k = 1, nz
            send_left(1,k,ll)  = s%mem(1,k,ll)
            send_left(2,k,ll)  = s%mem(2,k,ll)
            send_right(1,k,ll) = s%mem(nx-1,k,ll)
            send_right(2,k,ll) = s%mem(nx,k,ll)
          end do
        end do
        !$acc end parallel loop

        call mpi_timer%start("MPI: Communication")
        !$acc host_data use_device(send_left, send_right, recv_left, recv_right)
        call MPI_Irecv(recv_left, send_size, MPI_DOUBLE_PRECISION, &
                       left_rank, 1, MPI_COMM_WORLD, req(1), ierr)
        call MPI_Irecv(recv_right, send_size, MPI_DOUBLE_PRECISION, &
                       right_rank, 2, MPI_COMM_WORLD, req(2), ierr)
        call MPI_Isend(send_right, send_size, MPI_DOUBLE_PRECISION, &
                       right_rank, 1, MPI_COMM_WORLD, req(3), ierr)
        call MPI_Isend(send_left, send_size, MPI_DOUBLE_PRECISION, &
                       left_rank, 2, MPI_COMM_WORLD, req(4), ierr)
        call MPI_Waitall(4, req, status, ierr)
        !$acc end host_data
        call mpi_timer%stop()

        ! Unpack halos back on device
        !$acc parallel loop collapse(2) present(s%mem, recv_left, recv_right)
        do ll = 1, NVARS
          do k = 1, nz
            s%mem(-1,k,ll)   = recv_left(1,k,ll)
            s%mem(0,k,ll)    = recv_left(2,k,ll)
            s%mem(nx+1,k,ll) = recv_right(1,k,ll)
            s%mem(nx+2,k,ll) = recv_right(2,k,ll)
          end do
        end do
        !$acc end parallel loop

      else
#endif
        ! Prepare data to send on host
        do ll = 1, NVARS
          do k = 1, nz
            send_left(1,k,ll)  = s%mem(1,k,ll)
            send_left(2,k,ll)  = s%mem(2,k,ll)
            send_right(1,k,ll) = s%mem(nx-1,k,ll)
            send_right(2,k,ll) = s%mem(nx,k,ll)
          end do
        end do

        call mpi_timer%start("MPI: Communication")

        ! Non-blocking send/receive with neighbors
        call MPI_Irecv(recv_left, send_size, MPI_DOUBLE_PRECISION, &
                       left_rank, 1, MPI_COMM_WORLD, req(1), ierr)
        call MPI_Irecv(recv_right, send_size, MPI_DOUBLE_PRECISION, &
                       right_rank, 2, MPI_COMM_WORLD, req(2), ierr)
        call MPI_Isend(send_right, send_size, MPI_DOUBLE_PRECISION, &
                       right_rank, 1, MPI_COMM_WORLD, req(3), ierr)
        call MPI_Isend(send_left, send_size, MPI_DOUBLE_PRECISION, &
                       left_rank, 2, MPI_COMM_WORLD, req(4), ierr)
        call MPI_Waitall(4, req, status, ierr)

        call mpi_timer%stop()

        ! Copy received data to halos
        do ll = 1, NVARS
          do k = 1, nz
            s%mem(-1,k,ll)   = recv_left(1,k,ll)
            s%mem(0,k,ll)    = recv_left(2,k,ll)
            s%mem(nx+1,k,ll) = recv_right(1,k,ll)
            s%mem(nx+2,k,ll) = recv_right(2,k,ll)
          end do
        end do
#ifdef USE_OPENACC
      end if
#endif

    end if

#ifdef USE_OPENACC
    call nvtx_pop()
#endif

  end subroutine exchange_halo_x

  !> @brief Apply vertical boundary conditions into z-direction halos.
  !! @param[inout] s State whose vertical halos are updated.
  !! @param[in] ref Reference density used for momentum scaling at vertical edges.
  !! @details Enforces w=0 at top/bottom, scales u to maintain mass consistency,
  !!          and mirrors scalar fields.
  subroutine exchange_halo_z(s,ref)
	implicit none
	class(atmospheric_state), intent(inout) :: s
	class(reference_state), intent(in) :: ref
	integer :: i, ll

#ifdef USE_OPENACC
	call nvtx_push("Halo Exchange Z")
#endif

#ifdef USE_OPENACC
		  !$acc parallel loop collapse(2) present(s%mem, ref%density)
		  do ll = 1, NVARS
			do i = 1-hs,nx+hs
			  if (ll == I_WMOM) then
				s%mem(i,-1,ll) = 0.0_wp
			s%mem(i,0,ll) = 0.0_wp
			s%mem(i,nz+1,ll) = 0.0_wp
			s%mem(i,nz+2,ll) = 0.0_wp
		  else if (ll == I_UMOM) then
			s%mem(i,-1,ll)   = s%mem(i,1,ll) /  &
				ref%density(1) * ref%density(-1)
			s%mem(i,0,ll)    = s%mem(i,1,ll) /  &
				ref%density(1) * ref%density(0)
			s%mem(i,nz+1,ll) = s%mem(i,nz,ll) / &
				ref%density(nz) * ref%density(nz+1)
			s%mem(i,nz+2,ll) = s%mem(i,nz,ll) / &
				ref%density(nz) * ref%density(nz+2)
		  else
			s%mem(i,-1,ll) = s%mem(i,1,ll)
			s%mem(i,0,ll) = s%mem(i,1,ll)
			s%mem(i,nz+1,ll) = s%mem(i,nz,ll)
				s%mem(i,nz+2,ll) = s%mem(i,nz,ll)
			  end if
			end do
		  end do
		  !$acc end parallel loop
#else
		  do ll = 1, NVARS
			do i = 1-hs,nx+hs
			  if (ll == I_WMOM) then
				s%mem(i,-1,ll) = 0.0_wp
			s%mem(i,0,ll) = 0.0_wp
			s%mem(i,nz+1,ll) = 0.0_wp
			s%mem(i,nz+2,ll) = 0.0_wp
		  else if (ll == I_UMOM) then
			s%mem(i,-1,ll)   = s%mem(i,1,ll) /  &
				ref%density(1) * ref%density(-1)
			s%mem(i,0,ll)    = s%mem(i,1,ll) /  &
				ref%density(1) * ref%density(0)
			s%mem(i,nz+1,ll) = s%mem(i,nz,ll) / &
				ref%density(nz) * ref%density(nz+1)
			s%mem(i,nz+2,ll) = s%mem(i,nz,ll) / &
				ref%density(nz) * ref%density(nz+2)
		  else
			s%mem(i,-1,ll) = s%mem(i,1,ll)
			s%mem(i,0,ll) = s%mem(i,1,ll)
			s%mem(i,nz+1,ll) = s%mem(i,nz,ll)
			s%mem(i,nz+2,ll) = s%mem(i,nz,ll)
		end if
		end do
	  end do
#endif

#ifdef USE_OPENACC
	call nvtx_pop()
#endif

  end subroutine exchange_halo_z

  !> @brief Allocate reference profiles.
  subroutine new_ref(ref)
	implicit none
	class(reference_state), intent(inout) :: ref
	allocate(ref%density(1-hs:nz+hs))
	allocate(ref%denstheta(1-hs:nz+hs))
	allocate(ref%idens(nz+1))
	allocate(ref%idenstheta(nz+1))
	allocate(ref%pressure(nz+1))
  end subroutine new_ref

  !> @brief Deallocate reference profiles.
  subroutine del_ref(ref)
	implicit none
	class(reference_state), intent(inout) :: ref
	deallocate(ref%density)
	deallocate(ref%denstheta)
	deallocate(ref%idens)
	deallocate(ref%idenstheta)
	deallocate(ref%pressure)
  end subroutine del_ref

  !> @brief Allocate flux storage and component views.
  subroutine new_flux(flux)
	implicit none
	class(atmospheric_flux), intent(inout) :: flux
	if ( associated(flux%mem) ) deallocate(flux%mem)
	allocate(flux%mem(1:nx+1, 1:nz+1,NVARS))
	flux%dens => flux%mem(:,:,I_DENS)
	flux%umom => flux%mem(:,:,I_UMOM)
	flux%wmom => flux%mem(:,:,I_WMOM)
	flux%rhot => flux%mem(:,:,I_RHOT)
  end subroutine new_flux

  !> @brief Fill entire flux array with a constant value.
  subroutine set_flux(flux, xval)
	implicit none
	class(atmospheric_flux), intent(inout) :: flux
	real(wp), intent(in) :: xval
	if ( .not. associated(flux%mem) ) then
	  write(stderr,*) 'NOT ALLOCATED FLUX ERROR AT LINE ', __LINE__
	  stop
	end if
	flux%mem(:,:,:) = xval
  end subroutine set_flux

  !> @brief Deallocate flux storage and nullify component pointers.
  subroutine del_flux(flux)
	implicit none
	class(atmospheric_flux), intent(inout) :: flux
	if ( associated(flux%mem) ) deallocate(flux%mem)
	nullify(flux%dens)
	nullify(flux%umom)
	nullify(flux%wmom)
	nullify(flux%rhot)
  end subroutine del_flux

  !> @brief Allocate tendency storage and component views.
  subroutine new_tendency(tend)
	implicit none
	class(atmospheric_tendency), intent(inout) :: tend
	if ( associated(tend%mem) ) deallocate(tend%mem)
	allocate(tend%mem(nx, nz,NVARS))
	tend%dens => tend%mem(:,:,I_DENS)
	tend%umom => tend%mem(:,:,I_UMOM)
	tend%wmom => tend%mem(:,:,I_WMOM)
	tend%rhot => tend%mem(:,:,I_RHOT)
  end subroutine new_tendency

  !> @brief Fill entire tendency array with a constant value.
  subroutine set_tendency(tend, xval)
	implicit none
	class(atmospheric_tendency), intent(inout) :: tend
	real(wp), intent(in) :: xval
	if ( .not. associated(tend%mem) ) then
	  write(stderr,*) 'NOT ALLOCATED FLUX ERROR AT LINE ', __LINE__
	  stop
	end if
	tend%mem(:,:,:) = xval
  end subroutine set_tendency

  !> @brief Deallocate tendency storage and nullify component pointers.
  subroutine del_tendency(tend)
	implicit none
	class(atmospheric_tendency), intent(inout) :: tend
	if ( associated(tend%mem) ) deallocate(tend%mem)
	nullify(tend%dens)
	nullify(tend%umom)
	nullify(tend%wmom)
	nullify(tend%rhot)
  end subroutine del_tendency

  !> @brief Assignment overload to copy entire atmospheric_state contents.
  subroutine state_equal_to_state(x,y)
	implicit none
	type(atmospheric_state), intent(inout) :: x
	type(atmospheric_state), intent(in) :: y
	x%mem(:,:,:) = y%mem(:,:,:)
  end subroutine state_equal_to_state

end module module_types
