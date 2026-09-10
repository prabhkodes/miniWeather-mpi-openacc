!> @brief Main driver for the 2-D compressible Euler solver.
!! @details Initializes domain/physics, runs the Runge-Kutta time stepping with
!!          dimensional splitting, optionally outputs diagnostics, and reports
!!          mass/energy conservation plus wall-clock timing. Command line
!!          arguments allow overriding grid size, simulation length, and output
!!          frequency.
program atmosphere_model
  use mpi
  use calculation_types, only : wp
  use module_physics, only : dt, oldstat, newstat, flux, tend, ref
  use module_physics, only : init, finalize
  use module_physics, only : rungekutta, total_mass_energy
  use module_output, only : create_output, write_record, close_output
  use dimensions , only : sim_time, output_freq, set_dimensions
  use dimensions , only : default_nx, default_sim_time, default_output_freq
  use dimensions , only : setup_domain_decomposition, nx_global, nx, nprocs
  use iodir, only : stdout
  use parallel_timer
#ifdef USE_OPENACC
  use openacc
  use module_nvtx
#endif
#ifdef USE_OPENMP
  use omp_lib
#endif
  implicit none

  real(wp) :: etime            !< Current simulation time [s].
  real(wp) :: ptime            !< Interval used for periodic progress reporting.
  real(wp) :: output_counter   !< Accumulator that triggers output writes.
  real(wp) :: pctime           !< Percentage of simulation completed.
  real(wp) :: mass0, te0       !< Initial total mass and energy.
  real(wp) :: mass1, te1       !< Final total mass and energy.
  integer(8) :: t1, t2, rate   !< Raw tick counters and tick rate from `system_clock`.

  integer :: ierr              !< MPI error status.
  integer :: my_rank           !< Rank of current MPI task.
  integer :: dt_values(8)      !< Values returned by `date_and_time`.
  integer :: arg_count         !< Number of CLI arguments detected.
  integer :: arg_status        !< I/O status for CLI parsing.
  integer :: nx_cli            !< Grid cells in x supplied on CLI (global).
  real(wp) :: sim_time_cli     !< Total simulation time supplied on CLI [s].
  real(wp) :: output_freq_cli  !< Output cadence supplied on CLI [s].
  character(len=64) :: arg     !< Temporary buffer used when reading CLI values.
  real(8) :: final_duration    !< Wall-clock duration in seconds computed from ticks.
  integer :: num_procs         !< Number of MPI ranks participating.


  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

#ifdef USE_OPENACC
  ! Assign each MPI rank to a different GPU (round-robin)
  call nvtx_push("GPU Initialization")
  block
    integer :: num_gpus, my_gpu
    num_gpus = acc_get_num_devices(acc_device_nvidia)
    if (num_gpus > 0) then
      my_gpu = mod(my_rank, num_gpus)
      call acc_set_device_num(my_gpu, acc_device_nvidia)
    end if
  end block
  call nvtx_pop()
#endif

  if (my_rank == 0) then
    write(stdout,*) "================= Execution Info =================="
    write(stdout,'(a,i4)') "  Number of MPI tasks: ", num_procs
#ifdef USE_OPENMP
    write(stdout,'(a,i4)') "  Number of OpenMP threads: ", omp_get_max_threads()
    write(stdout,*) "  OpenMP: ENABLED"
#else
    write(stdout,*) "  OpenMP: DISABLED"
#endif
#ifdef USE_OPENACC
    write(stdout,'(a,i4)') "  Number of GPUs available: ", acc_get_num_devices(acc_device_nvidia)
    write(stdout,*) "  OpenACC: ENABLED"
#else
    write(stdout,*) "  OpenACC: DISABLED"
#endif
    write(stdout,*) "==================================================="
  end if

#ifdef USE_OPENACC
  ! Print GPU assignment for each MPI rank
  write(stdout,'(a,i4,a,i2)') "  MPI Rank ", my_rank, " -> GPU ", acc_get_device_num(acc_device_nvidia)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

  nx_cli = default_nx
  sim_time_cli = default_sim_time
  output_freq_cli = default_output_freq

  arg_count = command_argument_count()

  if (arg_count >= 1) then
    call get_command_argument(1, arg)
    read(arg, *, iostat=arg_status) nx_cli
    if (arg_status /= 0) nx_cli = default_nx
  end if

  if (arg_count >= 2) then
    call get_command_argument(2, arg)
    read(arg, *, iostat=arg_status) sim_time_cli
    if (arg_status /= 0) sim_time_cli = default_sim_time
  end if

  if (arg_count >= 3) then
    call get_command_argument(3, arg)
    read(arg, *, iostat=arg_status) output_freq_cli
    if (arg_status /= 0) output_freq_cli = default_output_freq
  end if

  call set_dimensions(nx_cli, sim_time_cli, output_freq_cli)

  ! Setup MPI domain decomposition (divides nx among processes)
  call setup_domain_decomposition(my_rank, num_procs)

  if (my_rank == 0) then
    write(stdout,*) "--------------- Domain Decomposition --------------"
    write(stdout,'(a,i6)') "  Global nx: ", nx_global
    write(stdout,'(a,i6)') "  Local nx per process: ~", nx_global/num_procs
    write(stdout,*) "---------------------------------------------------"
  end if

  ! --- PRINT START TIME (Rank 0 Only) ---
  if (my_rank == 0) then
      call date_and_time(values=dt_values)
      write(stdout, *) 'SIMPLE ATMOSPHERIC MODEL STARTING.'
      ! Format: HH:MM:SS.Milliseconds
      write(stdout, '(a,I2.2,a,I2.2,a,I2.2,a,I3.3)') &
        "Wall Clock Start: ", dt_values(5), ":", dt_values(6), ":", dt_values(7), ".", dt_values(8)
  end if

  call init(etime,output_counter,dt)
  call total_mass_energy(mass0,te0)

  ! --- For benchmark, removing file i/o
  ! call create_output( )
  ! call write_record(oldstat,ref,etime)

#ifdef USE_OPENACC
  call nvtx_push("Data Transfer to GPU")
  !$acc enter data create(oldstat, newstat, flux, tend, ref)
  !$acc enter data copyin(oldstat%mem, newstat%mem, flux%mem, tend%mem)
  !$acc enter data copyin(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
  !$acc enter data attach(oldstat%mem, newstat%mem, flux%mem, tend%mem)
  !$acc enter data attach(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
  call nvtx_pop()
#endif

  ! Get initial tick count and the clock rate (ticks per second)
  call system_clock(t1, rate) 

#ifdef USE_OPENACC
  call nvtx_push("Main Time Loop")
#endif
  ptime = int(sim_time/10.0)
  do while (etime < sim_time)

    if (etime + dt > sim_time) dt = sim_time - etime

    call rungekutta(oldstat,newstat,flux,tend,dt)

    if ( mod(etime,ptime) < dt ) then
      pctime = (etime/sim_time)*100.0_wp
      if (my_rank == 0) then
         write(stdout,'(1x,a,i2,a)') 'TIME PERCENT : ', int(pctime), '%'
      end if
    end if

    etime = etime + dt
    output_counter = output_counter + dt

    if (output_counter >= output_freq) then
      output_counter = output_counter - output_freq
#ifdef USE_OPENACC
      call nvtx_push("Data Update from GPU")
      !$acc update self(oldstat%mem)
      call nvtx_pop()
#endif
        ! --- For benchmark, removing file i/o
      ! call write_record(oldstat,ref,etime)
    end if

  end do
#ifdef USE_OPENACC
  call nvtx_pop()  ! End Main Time Loop
#endif

  call total_mass_energy(mass1,te1)
  
  ! --- For benchmark, removing file i/o
  ! call close_output( )

#ifdef USE_OPENACC
  call nvtx_push("Data Cleanup from GPU")
  !$acc exit data delete(oldstat%mem, newstat%mem, flux%mem, tend%mem)
  !$acc exit data delete(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
  !$acc exit data delete(oldstat, newstat, flux, tend, ref)
  call nvtx_pop()
#endif

  if (my_rank == 0) then
      write(stdout,*) "----------------- Atmosphere check ----------------"
      write(stdout,*) "Fractional Delta Mass  : ", (mass1-mass0)/mass0
      write(stdout,*) "Fractional Delta Energy: ", (te1-te0)/te0
      write(stdout,*) "---------------------------------------------------"
  end if

  call finalize()
  call print_timing_results()
  
  ! Get final tick count
  call system_clock(t2)

  ! --- PRINT END TIME (Rank 0 Only) ---
  if (my_rank == 0) then
      call date_and_time(values=dt_values)
      write(stdout,*) "SIMPLE ATMOSPHERIC MODEL RUN COMPLETED."
      
      ! 1. Print Wall Clock with Milliseconds
      write(stdout, '(a,I2.2,a,I2.2,a,I2.2,a,I3.3)') &
        "Wall Clock End:   ", dt_values(5), ":", dt_values(6), ":", dt_values(7), ".", dt_values(8)
      
      ! 2. Print Duration with Microsecond Precision
      final_duration = dble(t2-t1)/dble(rate)
      write(stdout, '(a,F18.6,a)') "USED CPU TIME:    ", final_duration, " seconds"
  end if

  call MPI_FINALIZE(ierr)

end program atmosphere_model
