!> @brief Numeric kinds shared across the solver.
!! @details Centralizes the choice of working precision for floating point,
!!          integers, and NetCDF I/O kinds so the rest of the code can stay
!!          type-agnostic.
module calculation_types
  use netcdf, only : nf90_real, nf90_double
  use iso_fortran_env
  implicit none
  public
  integer, parameter :: wp = real64       !< Working precision for floating point math.
  integer, parameter :: iowp = nf90_double !< NetCDF type corresponding to `wp`.
  integer, parameter :: ip = int32        !< Default integer precision for indices.
end module calculation_types

!> @brief Physical constants used by the compressible Euler equations.
!! @details Provides thermodynamic constants and derived coefficients such as
!!          Exner function factors. Values are expressed in SI units.
module physical_constants
  use calculation_types, only : wp
  implicit none
  public
  real(wp), parameter :: pi = 3.14159265358979323846264338327_wp !< Ratio of circumference to diameter.
  real(wp), parameter :: hpi = 0.5_wp * pi                        !< Half-pi convenience constant.
  real(wp), parameter :: grav = 9.80665_wp                        !< Acceleration due to gravity [m s^-2].
  real(wp), parameter :: boltzk = 1.3806490e-23_wp                !< Boltzmann constant [J K^-1].
  real(wp), parameter :: navgdr = 6.02214076e23_wp                !< Avogadro number [mol^-1].
  real(wp), parameter :: rgasmol = navgdr*boltzk                  !< Universal gas constant [J mol^-1 K^-1].
  real(wp), parameter :: amd = 28.96454_wp                        !< Mean molecular weight of dry air [g mol^-1].
  real(wp), parameter :: amw = 18.01528_wp                        !< Molecular weight of water vapor [g mol^-1].
  real(wp), parameter :: rgas = (rgasmol/amd)*1000.0_wp           !< Specific gas constant for dry air [J kg^-1 K^-1].
  real(wp), parameter :: cp = 3.5_wp*rgas                         !< Specific heat at constant pressure [J kg^-1 K^-1].
  real(wp), parameter :: cv = 2.5_wp*rgas                         !< Specific heat at constant volume [J kg^-1 K^-1].
  real(wp), parameter :: rd = rgas                                !< Alias for dry gas constant to match notation.
  real(wp), parameter :: p0 = 1.e5_wp                             !< Reference surface pressure [Pa].
  real(wp), parameter :: t0 = 298.0_wp                            !< Reference temperature [K].
  real(wp), parameter :: cdocv = 3.5_wp/2.5_wp                    !< cp/cv ratio for dry air.
  real(wp), parameter :: cvocd = 2.5_wp/3.5_wp                    !< cv/cp ratio (inverse of cdocv).
  real(wp), parameter :: rdocp = rd/cp                            !< Rd/Cp used in Exner function.
  real(wp), parameter :: rdocv = rd/cv                            !< Rd/Cv exponent in equation of state.
  real(wp), parameter :: c0 = (rd**cdocv)*p0**(-rdocv)            !< Ideal gas constant re-cast for pressure relation.
  real(wp), parameter :: theta0 = 300.0_wp                        !< Background potential temperature [K].
  real(wp), parameter :: exner0 = 1.0_wp                          !< Background Exner function (non-dimensional).
end module physical_constants

!> @brief Domain-scale physical and numerical parameters.
!! @details Contains domain geometry, feature placement, hyper-viscosity
!!          strength, CFL targets, and perturbation amplitudes used by the
!!          solver and initialization routines.
module physical_parameters
  use calculation_types, only : wp
  implicit none
  public
  real(wp), parameter :: xlen = 20000.0_wp     !< Domain length in x-direction [m].
  real(wp), parameter :: zlen = 10000.0_wp     !< Domain height in z-direction [m].
  real(wp), parameter :: hxlen = 0.5_wp*xlen   !< Half-domain length used in analytic profiles [m].
  real(wp), parameter :: p1 = xlen/10.0_wp     !< Characteristic wavelength for perturbations [m].
  real(wp), parameter :: p2 = 4.0_wp*p1        !< Secondary wavelength scale [m].
  real(wp), parameter :: p3 = 2.0_wp*p1        !< Tertiary wavelength scale [m].
  real(wp), parameter :: mnt_width = xlen/8.0_wp !< Mountain half-width for topography [m].
  real(wp), parameter :: hv_beta = 0.05_wp     !< Hyper-viscosity strength parameter (0-1).
  real(wp), parameter :: cfl = 1.5_wp          !< CFL target for time stepping stability.
  real(wp), parameter :: max_speed = 450.0_wp  !< Assumed max signal speed (wind + sound) [m s^-1].
  real(wp), parameter :: x0 = mnt_width        !< Perturbation center in x [m].
  real(wp), parameter :: z0 = 1000.0_wp        !< Perturbation center in z [m].
  real(wp), parameter :: xrad = 500.0_wp       !< Perturbation e-folding radius in x [m].
  real(wp), parameter :: zrad = 500.0_wp       !< Perturbation e-folding radius in z [m].
  real(wp), parameter :: amp = 0.01_wp         !< Potential temperature perturbation amplitude [non-dimensional].
end module physical_parameters

!> @brief Indexing helper constants for state vectors and directions.
!! @details Encapsulates common integer codes used throughout the solver to
!!          reference variables in packed arrays and to distinguish spatial
!!          directions.
module indexing
  implicit none
  public
  integer, parameter :: NVARS = 4        !< Number of prognostic variables.
  integer, parameter :: STEN_SIZE = 4    !< Stencil width (cells) for reconstruction.
  integer, parameter :: I_DENS  = 1      !< Index of density in packed vectors.
  integer, parameter :: I_UMOM  = 2      !< Index of x-momentum in packed vectors.
  integer, parameter :: I_WMOM  = 3      !< Index of z-momentum in packed vectors.
  integer, parameter :: I_RHOT  = 4      !< Index of rho*theta in packed vectors.
  integer, parameter :: DIR_X = 1        !< Code for x-direction sweeps.
  integer, parameter :: DIR_Z = 2        !< Code for z-direction sweeps.
end module indexing

!> @brief Convenience access to standard I/O units.
!! @details Wraps compiler-provided unit numbers to keep I/O statements readable.
module iodir
  use iso_fortran_env
  public
  integer, public, parameter :: stdin = input_unit   !< Standard input unit.
  integer, public, parameter :: stdout = output_unit !< Standard output unit.
  integer, public, parameter :: stderr = error_unit  !< Standard error unit.
end module iodir

!> @brief Halo and indexing offsets for domain decomposition.
!! @details Defines the start indices for interior computation and the halo
!!          size (`hs`) assumed by halo exchange routines.
module parallel_parameters
  implicit none
  public
  integer, parameter :: i_beg = 1 !< First i-index of owned interior cells.
  integer, parameter :: k_beg = 1 !< First k-index of owned interior cells.
  integer, parameter :: hs = 2    !< Halo size (ghost cell thickness) in each direction.
end module parallel_parameters

!> @brief Gauss-Legendre quadrature constants.
!! @details Stores three-point quadrature abscissas and weights for integrating
!!          vertical source terms and flux divergences.
module legendre_quadrature
  use calculation_types, only : wp
  implicit none
  public
  integer , parameter :: nqpoints = 3 !< Number of quadrature points.
  real(wp), parameter :: qpoints(nqpoints) =    &
      [ 0.112701665379258311482073460022E0_wp , & !< Quadrature abscissa #1 (0-1 interval).
        0.500000000000000000000000000000E0_wp , & !< Quadrature abscissa #2 (0-1 interval).
        0.887298334620741688517926539980E0_wp ]  !< Quadrature abscissa #3 (0-1 interval).
  real(wp), parameter :: qweights(nqpoints) =   &
      [ 0.277777777777777777777777777779E0_wp , & !< Weight for abscissa #1.
        0.444444444444444444444444444444E0_wp , & !< Weight for abscissa #2.
        0.277777777777777777777777777779E0_wp ]   !< Weight for abscissa #3.
end module legendre_quadrature

!> @brief Domain sizing and MPI decomposition utilities.
!! @details Tracks global/local grid sizes, simulation length, output cadence,
!!          and neighbor ranks. Provides helpers to set dimensions from CLI and
!!          to compute the Cartesian 1-D partitioning in x.
module dimensions
  use calculation_types, only : wp
  use physical_parameters, only : zlen, xlen
  use indexing
  implicit none
  public
  integer, parameter :: default_nx = 100            !< Default number of global grid cells in x.
  real(wp), parameter :: default_sim_time = 1000.0_wp !< Default simulation duration [s].
  real(wp), parameter :: default_output_freq = 10.0_wp !< Default output cadence [s].

  integer :: nx = default_nx                        !< Local number of x cells (overwritten by decomposition).
  integer :: nz = int(default_nx * zlen/xlen)       !< Number of vertical cells derived from aspect ratio.
  real(wp) :: sim_time = default_sim_time           !< Total simulation time requested [s].
  real(wp) :: output_freq = default_output_freq     !< Interval between diagnostics [s].

  ! MPI domain decomposition variables
  integer :: nprocs = 1           !< Total number of MPI processes.
  integer :: myrank = 0           !< Rank of this process.
  integer :: nx_local             !< Local nx after partitioning.
  integer :: nx_global            !< Global nx (pre-partition).
  integer :: i_start_global       !< Global starting i-index owned by this rank (1-based).
  integer :: i_end_global         !< Global ending i-index owned by this rank (1-based).
  integer :: left_rank            !< MPI rank to the immediate left (periodic).
  integer :: right_rank           !< MPI rank to the immediate right (periodic).

contains

  !> @brief Set global grid size, total simulation time, and output cadence.
  !! @param[in] nx_in Global number of x cells requested.
  !! @param[in] sim_time_in Total simulation time [s].
  !! @param[in] output_freq_in Output interval [s].
  !! @details Updates `nx`, `nx_global`, `sim_time`, `output_freq`, and recomputes
  !!          the vertical resolution assuming a fixed aspect ratio between x
  !!          and z extents.
  subroutine set_dimensions(nx_in, sim_time_in, output_freq_in)
    integer, intent(in) :: nx_in
    real(wp), intent(in) :: sim_time_in
    real(wp), intent(in) :: output_freq_in

    nx = nx_in
    nx_global = nx_in
    sim_time = sim_time_in
    output_freq = output_freq_in
    nz = int(nx * zlen/xlen)
  end subroutine set_dimensions

  !> @brief Partition the x-dimension among MPI ranks.
  !! @param[in] rank MPI rank of the caller.
  !! @param[in] num_procs Total number of ranks in the communicator.
  !! @details Splits `nx_global` into approximately equal contiguous chunks,
  !!          distributing any remainder to the lowest ranks. Records local
  !!          extents and periodic neighbors, and overwrites `nx` with the local
  !!          cell count for downstream routines.
  subroutine setup_domain_decomposition(rank, num_procs)
    implicit none
    integer, intent(in) :: rank, num_procs
    integer :: base_nx, remainder

    myrank = rank
    nprocs = num_procs
    nx_global = nx

    ! Divide nx among processes
    base_nx = nx_global / nprocs
    remainder = mod(nx_global, nprocs)

    ! Distribute remainder among first processes
    if (myrank < remainder) then
      nx_local = base_nx + 1
      i_start_global = myrank * (base_nx + 1) + 1
    else
      nx_local = base_nx
      i_start_global = remainder * (base_nx + 1) + (myrank - remainder) * base_nx + 1
    end if
    i_end_global = i_start_global + nx_local - 1

    ! Update nx to local value for this process
    nx = nx_local

    ! Neighbors (periodic in x)
    left_rank = mod(myrank - 1 + nprocs, nprocs)
    right_rank = mod(myrank + 1, nprocs)

  end subroutine setup_domain_decomposition

end module dimensions
