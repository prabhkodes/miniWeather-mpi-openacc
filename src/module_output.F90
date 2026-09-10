!> @brief Parallel NetCDF output helpers.
!! @details Collects state slices from each MPI rank and writes them serially to
!!          a single NetCDF file. Rank 0 owns the file handle and receives data
!!          from workers; local arrays are used as staging buffers to avoid
!!          exposing solver internals.
module module_output
  use mpi
  use calculation_types, only : wp, iowp
  use parallel_parameters, only : i_beg, k_beg
  use dimensions, only : nx, nz, nx_global, nprocs, myrank, i_start_global
  use module_types, only : atmospheric_state, reference_state
  use iodir, only : stderr, stdout
  use netcdf

  implicit none

  private

  public :: create_output
  public :: write_record
  public :: close_output

  ! Local arrays (Small slices)
  real(wp), allocatable, dimension(:,:) :: dens   !< Local density slice [kg m^-3].
  real(wp), allocatable, dimension(:,:) :: uwnd   !< Local u-wind (x-momentum / rho) slice [m s^-1].
  real(wp), allocatable, dimension(:,:) :: wwnd   !< Local w-wind (z-momentum / rho) slice [m s^-1].
  real(wp), allocatable, dimension(:,:) :: theta  !< Local potential temperature slice [K].

  ! NetCDF IDs (Only valid on Rank 0)
  integer :: ncid                                  !< NetCDF file identifier (rank 0).
  integer :: dens_varid, uwnd_varid, wwnd_varid, theta_varid, t_varid !< Variable identifiers in the NetCDF file.
  integer :: rec_out                               !< Current output record index (1-based NetCDF unlimited dimension).

  contains

  !> @brief Create the NetCDF file and define variables.
  !! @details Allocates local staging arrays on every rank. Rank 0 creates the
  !!          `output.nc` file, defines dimensions/variables using the global
  !!          grid sizes, and ends define mode. The record counter is reset.
  subroutine create_output
    implicit none
    integer :: t_dimid, x_dimid, z_dimid

    ! 1. Allocate LOCAL arrays (Size: local nx, nz)
    allocate(dens(nx,nz))
    allocate(uwnd(nx,nz))
    allocate(wwnd(nx,nz))
    allocate(theta(nx,nz))

    ! 2. Setup File (RANK 0 ONLY)
    if (myrank == 0) then
      ! Standard Serial Create (clobber overwrites existing files)
      call ncwrap(nf90_create('output.nc', nf90_clobber, ncid), __LINE__)

      ! Define dimensions using GLOBAL size
      call ncwrap(nf90_def_dim(ncid, 'time', nf90_unlimited, t_dimid), __LINE__)
      call ncwrap(nf90_def_dim(ncid, 'x', nx_global, x_dimid), __LINE__) 
      call ncwrap(nf90_def_dim(ncid, 'z', nz, z_dimid), __LINE__)

      call ncwrap(nf90_def_var(ncid,'time',iowp,[t_dimid],t_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'rho',iowp, [x_dimid,z_dimid,t_dimid],dens_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'u',iowp,   [x_dimid,z_dimid,t_dimid],uwnd_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'w',iowp,   [x_dimid,z_dimid,t_dimid],wwnd_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'theta',iowp,[x_dimid,z_dimid,t_dimid],theta_varid), __LINE__)
      
      call ncwrap(nf90_enddef(ncid), __LINE__)
    end if

    rec_out = 1
  end subroutine create_output

  !> @brief Gather local state from all ranks and append one time record.
  !! @param[in] atmostat Prognostic state (density, momentum, rho*theta) local to the rank.
  !! @param[in] ref Reference hydrostatic profile needed to convert to velocities/potential temperature.
  !! @param[in] etime Simulation time stamp to associate with this record [s].
  !! @details Rank 0 writes its own chunk then receives each worker's chunk in
  !!          order (dimension and data messages) before writing them to their
  !!          appropriate offsets. Workers send derived velocity and theta
  !!          fields as dense buffers.
  subroutine write_record(atmostat,ref,etime)
    implicit none
    type(atmospheric_state), intent(in) :: atmostat
    type(reference_state), intent(in) :: ref
    real(wp), intent(in) :: etime
    integer :: i, k, src, ierr
    integer, dimension(3) :: st3, ct3
    integer, dimension(1) :: st1, ct1
    real(wp), dimension(1) :: etime_arr
    integer :: status(MPI_STATUS_SIZE)
    
    ! Temp variables for Rank 0 receiving
    integer :: r_nx, r_istart
    real(wp), allocatable, dimension(:,:) :: temp_buf

    ! --- 1. Compute Local Fields ---
    do k = 1, nz
      do i = 1, nx
        dens(i,k) = atmostat%dens(i,k)
        uwnd(i,k) = atmostat%umom(i,k)/(ref%density(k)+dens(i,k))
        wwnd(i,k) = atmostat%wmom(i,k)/(ref%density(k)+dens(i,k))
        theta(i,k) = (atmostat%rhot(i,k) + ref%denstheta(k)) / &
            (ref%density(k) + dens(i,k)) - ref%denstheta(k)/ref%density(k)
      end do
    end do

    ! --- 2. Serial I/O Strategy ---
    if (myrank == 0) then
       
       ! A. Write Time
       st1 = [ rec_out ]
       ct1 = [ 1 ]
       etime_arr(1) = etime
       call ncwrap(nf90_put_var(ncid, t_varid, etime_arr, start=st1, count=ct1), __LINE__)

       ! B. Write MY OWN data (Rank 0)
       st3 = [ i_start_global, 1, rec_out ]
       ct3 = [ nx, nz, 1 ]
       call ncwrap(nf90_put_var(ncid, dens_varid, dens, start=st3, count=ct3), __LINE__)
       call ncwrap(nf90_put_var(ncid, uwnd_varid, uwnd, start=st3, count=ct3), __LINE__)
       call ncwrap(nf90_put_var(ncid, wwnd_varid, wwnd, start=st3, count=ct3), __LINE__)
       call ncwrap(nf90_put_var(ncid, theta_varid,theta,start=st3, count=ct3), __LINE__)

       ! C. Receive and Write Workers data ONE BY ONE
       do src = 1, nprocs - 1
          ! 1. Receive Dimensions
          call MPI_Recv(r_nx, 1, MPI_INTEGER, src, 100, MPI_COMM_WORLD, status, ierr)
          call MPI_Recv(r_istart, 1, MPI_INTEGER, src, 101, MPI_COMM_WORLD, status, ierr)
          
          ! 2. Allocate Buffer just for this chunk
          allocate(temp_buf(r_nx, nz))
          
          ! 3. Define where to write in file
          st3 = [ r_istart, 1, rec_out ]
          ct3 = [ r_nx, nz, 1 ]

          ! 4. Recv & Write DENS
          call MPI_Recv(temp_buf, r_nx*nz, MPI_DOUBLE_PRECISION, src, 200, MPI_COMM_WORLD, status, ierr)
          call ncwrap(nf90_put_var(ncid, dens_varid, temp_buf, start=st3, count=ct3), __LINE__)

          ! 5. Recv & Write UWND
          call MPI_Recv(temp_buf, r_nx*nz, MPI_DOUBLE_PRECISION, src, 201, MPI_COMM_WORLD, status, ierr)
          call ncwrap(nf90_put_var(ncid, uwnd_varid, temp_buf, start=st3, count=ct3), __LINE__)

          ! 6. Recv & Write WWND
          call MPI_Recv(temp_buf, r_nx*nz, MPI_DOUBLE_PRECISION, src, 202, MPI_COMM_WORLD, status, ierr)
          call ncwrap(nf90_put_var(ncid, wwnd_varid, temp_buf, start=st3, count=ct3), __LINE__)

          ! 7. Recv & Write THETA
          call MPI_Recv(temp_buf, r_nx*nz, MPI_DOUBLE_PRECISION, src, 203, MPI_COMM_WORLD, status, ierr)
          call ncwrap(nf90_put_var(ncid, theta_varid, temp_buf, start=st3, count=ct3), __LINE__)

          deallocate(temp_buf)
       end do

    else
       ! --- WORKER: Send Data to Rank 0 ---
       call MPI_Send(nx, 1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, ierr)
       call MPI_Send(i_start_global, 1, MPI_INTEGER, 0, 101, MPI_COMM_WORLD, ierr)
       
       call MPI_Send(dens, nx*nz, MPI_DOUBLE_PRECISION, 0, 200, MPI_COMM_WORLD, ierr)
       call MPI_Send(uwnd, nx*nz, MPI_DOUBLE_PRECISION, 0, 201, MPI_COMM_WORLD, ierr)
       call MPI_Send(wwnd, nx*nz, MPI_DOUBLE_PRECISION, 0, 202, MPI_COMM_WORLD, ierr)
       call MPI_Send(theta, nx*nz, MPI_DOUBLE_PRECISION, 0, 203, MPI_COMM_WORLD, ierr)
    end if

    rec_out = rec_out + 1
  end subroutine write_record

  !> @brief Close the NetCDF output and release buffers.
  !! @details Deallocates local staging arrays on every rank. Rank 0 closes the
  !!          NetCDF file handle.
  subroutine close_output
    implicit none
    
    if ( allocated(dens) ) deallocate(dens, uwnd, wwnd, theta)
    
    if (myrank == 0) then
      call ncwrap(nf90_close(ncid), __LINE__)
    end if
    
  end subroutine close_output

  !> @brief Abort on NetCDF errors with line context.
  !! @param[in] ierr Status code returned by a NetCDF call.
  !! @param[in] line Source code line number supplied by caller.
  !! @details If `ierr` is non-zero, prints an informative message to stderr and
  !!          aborts the MPI communicator.
  subroutine ncwrap(ierr,line)
    implicit none
    integer, intent(in) :: ierr
    integer, intent(in) :: line
    integer :: mpi_err
    if (ierr /= nf90_noerr) then
      write(stderr,*) 'NetCDF Error (Rank ', myrank, ') at line: ', line
      write(stderr,*) nf90_strerror(ierr)
      call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
    end if
  end subroutine ncwrap

end module module_output
