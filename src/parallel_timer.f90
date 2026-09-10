!> @brief Hierarchical MPI wall-clock timer utilities.
!! @details Provides scoped timers that can be started/stopped manually or
!!          automatically (via finalizers) to gather inclusive/exclusive
!!          microsecond timings across ranks. A stack tracks parent-child
!!          relationships so exclusive time discounts the time spent in nested
!!          regions. `print_timing_results` reduces the data across MPI ranks.
module parallel_timer
    use mpi
    implicit none
    
    ! --- Internal Storage ---
    integer, parameter :: MAX_TIMERS = 100 !< Maximum number of distinct labels that can be tracked.
    integer, parameter :: MAX_DEPTH = 50  !< Maximum nesting depth of simultaneous timers.
    integer, parameter :: LABEL_LEN = 64  !< Maximum length of each timer label.
    
    !> @brief Internal storage for one timer entry.
    !! @details Tracks aggregated timing statistics per unique label. `total_time`
    !!          stores inclusive wall time in microseconds, while `children_time`
    !!          captures time spent in nested timers to allow computation of
    !!          exclusive time.
    type :: timer_store_t
        character(len=LABEL_LEN) :: label = "" !< User supplied timer name.
        real(8) :: total_time = 0.0d0          !< Inclusive accumulated wall time [microseconds].
        real(8) :: children_time = 0.0d0       !< Time spent inside child timers [microseconds].
        integer :: calls = 0                   !< Number of times this timer has been closed.
        logical :: active = .false.            !< Flag indicating that the entry has been initialized.
    end type timer_store_t
    
    type(timer_store_t), save :: database(MAX_TIMERS) !< Persistent table of timers indexed by creation order.
    integer, save :: num_timers = 0                   !< Number of timers currently in the table.
    
    integer, save :: timer_stack(MAX_DEPTH) !< Stack of active timer indices for computing exclusivity.
    integer, save :: stack_ptr = 0          !< Pointer to the top of the stack (0 means empty).

    ! --- The Scoped Timer Object ---
    !> @brief RAII-style timer that records elapsed wall time between construction and destruction.
    !! @details `CTimer` instances look up or create a timer entry in `database`.
    !!          The `start` method pushes the timer index on a stack; `stop`
    !!          pops the stack, updates inclusive counts, and attributes time to
    !!          the parent timer's children total.
    type :: CTimer
        integer :: db_index = 0  !< Index into `database` for this timer instance.
        real(8) :: start_time = 0.0d0 !< Time returned by `MPI_Wtime` when the timer was started.
    contains
        procedure :: start => start_timer_sub
        ! 1. Use the new CLASS version for manual calls
        procedure :: stop => stop_timer 
        ! 2. Use the new TYPE version for auto cleanup
        final :: stop_timer_final 
    end type CTimer

    interface CTimer
        module procedure new_timer
    end interface

contains

    !> @brief Start or resume a named timer.
    !! @param[inout] this The scoped timer object being started.
    !! @param[in] label Human-readable label that keys the timer entry.
    !! @details If a label has been seen before, the existing entry is reused so
    !!          callers can repeatedly time the same routine. The timer index is
    !!          pushed to `timer_stack` to capture nesting for exclusive timing.
    subroutine start_timer_sub(this, label)
        class(CTimer), intent(inout) :: this
        character(len=*), intent(in) :: label
        integer :: i
        
        this%start_time = MPI_Wtime()
        this%db_index = -1
        
        do i = 1, num_timers
            if (trim(database(i)%label) == trim(label)) then
                this%db_index = i
                exit
            end if
        end do
        
        if (this%db_index == -1) then
            if (num_timers < MAX_TIMERS) then
                num_timers = num_timers + 1
                this%db_index = num_timers
                database(num_timers)%label = label
                database(num_timers)%active = .true.
            else
                print *, "Timer Error: Max timers exceeded"
            end if
        end if

        if (this%db_index > 0) then
            if (stack_ptr < MAX_DEPTH) then
                stack_ptr = stack_ptr + 1
                timer_stack(stack_ptr) = this%db_index
            end if
        end if
    end subroutine start_timer_sub

    !> @brief Factory for a scoped timer that starts immediately.
    !! @param[in] label Name of the timer entry.
    !! @return t Initialized timer with an active start time.
    function new_timer(label) result(t)
        character(len=*), intent(in) :: label
        type(CTimer) :: t
        call t%start(label)
    end function new_timer

    ! --- MAIN LOGIC (Changed to CLASS) ---
    !> @brief Stop a scoped timer and attribute elapsed wall time.
    !! @param[inout] this Timer to stop; no-op if it was not started.
    !! @details Computes elapsed microseconds using `MPI_Wtime`, increments call
    !!          counters, and propagates time spent to the parent timer via the
    !!          stack so that exclusive time can be derived.
    impure subroutine stop_timer(this)
        class(CTimer), intent(inout) :: this ! <--- Changed from TYPE to CLASS
        real(8) :: end_time, duration
        integer :: parent_idx
        
        if (this%db_index > 0) then
            end_time = MPI_Wtime()
            duration = (end_time - this%start_time) * 1.0d6 
            
            database(this%db_index)%total_time = database(this%db_index)%total_time + duration
            database(this%db_index)%calls = database(this%db_index)%calls + 1
            
            if (stack_ptr > 0) then
                if (timer_stack(stack_ptr) == this%db_index) then
                    stack_ptr = stack_ptr - 1
                    if (stack_ptr > 0) then
                        parent_idx = timer_stack(stack_ptr)
                        database(parent_idx)%children_time = database(parent_idx)%children_time + duration
                    end if
                end if
            end if
            
            this%db_index = 0 
        end if
    end subroutine stop_timer

    ! --- WRAPPER FOR FINALIZER (Takes TYPE) ---
    !> @brief Finalizer wrapper that stops timers that fall out of scope.
    !! @param[inout] this Timer being finalized.
    !! @details Enables RAII-style usage where `CTimer` local variables
    !!          automatically record elapsed time when the block exits.
    impure subroutine stop_timer_final(this)
        type(CTimer), intent(inout) :: this ! <--- Strictly TYPE for Final
        ! Just forward the call to the main logic
        call stop_timer(this) 
    end subroutine stop_timer_final

    ! --- Parallel Reporting ---
    !> @brief Reduce and print timing statistics across all MPI ranks.
    !! @details For each timer label, computes the maximum inclusive time,
    !!          maximum exclusive time (discounting nested timers), average
    !!          inclusive time, and total number of calls aggregated over the
    !!          communicator. Only rank 0 prints the table.
    subroutine print_timing_results()
        integer :: i, rank, ierr, world_size
        real(8) :: local_total, local_excl
        real(8) :: max_total, avg_total
        real(8) :: max_excl
        integer :: sum_calls
        
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)

        if (rank == 0) then
            print *, repeat("-", 100)
            print *, "PARALLEL TIMING STATISTICS (Microseconds)"
            print "(A30, A15, A15, A15, A10)", "Function", "Max Total", "Max Excl", "Avg Total", "Calls"
            print *, repeat("-", 100)
        end if

        do i = 1, num_timers
            local_total = database(i)%total_time
            local_excl  = database(i)%total_time - database(i)%children_time
            if (local_excl < 0.0d0) local_excl = 0.0d0
            
            call MPI_Reduce(local_total, max_total, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(local_total, avg_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            avg_total = avg_total / world_size

            call MPI_Reduce(local_excl, max_excl, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(database(i)%calls, sum_calls, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

            if (rank == 0) then
                print "(A30, F15.2, F15.2, F15.2, I10)", &
                      trim(database(i)%label), max_total, max_excl, avg_total, sum_calls
            end if
        end do
        
        if (rank == 0) then
            print *, repeat("-", 100)
            print *, "* Max Total/Excl: The slowest rank for that function."
            print *, "* Avg Total: Average wall time across all ranks."
        end if
    end subroutine print_timing_results

end module parallel_timer
