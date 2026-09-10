!> @brief NVTX range annotations for Nsight Systems.
!! @details When built with -DUSE_NVTX (set automatically by the OpenACC build,
!!          which links NVIDIA's libnvToolsExt) the routines push and pop named
!!          ranges that appear as labelled spans in an Nsight Systems timeline.
!!          Otherwise they compile to no-ops, so callers need no preprocessor
!!          guards and the model still builds on hosts without the CUDA toolkit.
module module_nvtx
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
  implicit none

#ifdef USE_NVTX
  interface
     function nvtxRangePushA(name) bind(C, name="nvtxRangePushA")
       import :: c_int, c_char
       integer(c_int) :: nvtxRangePushA
       character(kind=c_char), dimension(*) :: name
     end function nvtxRangePushA

     function nvtxRangePop() bind(C, name="nvtxRangePop")
       import :: c_int
       integer(c_int) :: nvtxRangePop
     end function nvtxRangePop
  end interface
#endif

contains

  !> Open a named NVTX range.
  subroutine nvtx_push(label)
    character(len=*), intent(in) :: label
#ifdef USE_NVTX
    character(kind=c_char,len=:), allocatable :: c_label
    integer(c_int) :: istat

    ! Convert Fortran string to C string (null-terminated)
    c_label = trim(label)//c_null_char
    istat = nvtxRangePushA(c_label)
#else
    associate(unused => label); end associate
#endif
  end subroutine nvtx_push

  !> Close the most recently opened NVTX range.
  subroutine nvtx_pop()
#ifdef USE_NVTX
    integer(c_int) :: istat
    istat = nvtxRangePop()
#endif
  end subroutine nvtx_pop

end module module_nvtx
