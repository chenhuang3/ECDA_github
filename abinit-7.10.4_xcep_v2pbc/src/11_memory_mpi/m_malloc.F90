!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_malloc
!! NAME
!!  m_malloc
!!
!! FUNCTION
!!   This module is used for tracing memory allocations/deallocations when we compile the code in DEBUG_MODE
!!
!! COPYRIGHT
!!  Copyright (C) 2014 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!!  *  This module uses several global variables to simplify the API. 
!!     Please do *not* follow this approach in other parts of the code.
!!
!!  *  This module is not thread-safe. However, this should not represent a significant limitation
!!     since memory-tracing is only enabled in debug mode.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_malloc
    
 use defs_basis, only : sp, dp, dpc, spc, std_out, i8b

#ifdef HAVE_MPI2
 use mpi
#endif

 implicit none

 private

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

 public :: malloc_init         ! Initialize the library
 public :: malloc_shutdown     ! Shutdown the library.
 !public :: malloc_enable      
 !public :: malloc_disable      
 !public :: malloc_reset
 !public :: malloc_change_level
 !public :: malloc_report
 !public :: malloc_dump  
 public :: abi_malloc          ! Allocate an allocatable entity and register the event.
 public :: abi_malloc_ptr      ! Allocate a pointer and register the event.
 public :: abi_free            ! Deallocate an allocatable entity and register the event.
 public :: abi_free_ptr        ! Dellocate a pointer and register the event.
!!***

 ! Generic interface for the allocation of allocatable arrays.
 interface abi_malloc
   ! Single precision real
   module procedure malloc_allocatable_real_sp1
   module procedure malloc_allocatable_real_sp2
   module procedure malloc_allocatable_real_sp3
   module procedure malloc_allocatable_real_sp4
   module procedure malloc_allocatable_real_sp5
   module procedure malloc_allocatable_real_sp6
   module procedure malloc_allocatable_real_sp7
   !
 end interface abi_malloc

 ! Generic interface for the allocation of pointers.
 interface abi_malloc_ptr
   ! Single precision real
   module procedure malloc_pointer_real_sp1
   module procedure malloc_pointer_real_sp2
   module procedure malloc_pointer_real_sp3
   module procedure malloc_pointer_real_sp4
   module procedure malloc_pointer_real_sp5
   module procedure malloc_pointer_real_sp6
   module procedure malloc_pointer_real_sp7
 end interface abi_malloc_ptr

 ! Generic interface for the deallocation of allocatable arrays.
 interface abi_free
   ! Single precision real
    module procedure free_allocatable_real_sp1
    module procedure free_allocatable_real_sp2
    module procedure free_allocatable_real_sp3
    module procedure free_allocatable_real_sp4
    module procedure free_allocatable_real_sp5
    module procedure free_allocatable_real_sp6
    module procedure free_allocatable_real_sp7
 end interface abi_free

 ! Generic interface for the deallocation of pointer.
 interface abi_free_ptr
   ! Single precision real
    module procedure free_pointer_real_sp1
    module procedure free_pointer_real_sp2
    module procedure free_pointer_real_sp3
    module procedure free_pointer_real_sp4
    module procedure free_pointer_real_sp5
    module procedure free_pointer_real_sp6
    module procedure free_pointer_real_sp7
 end interface abi_free_ptr

!=============
!private stuff
!=============
  !integer,public,save :: ABI_ALLOC_STAT_ABI
  integer,save :: MEM_STAT  
    ! allocation/deallocation status returned by the last call to allocate/deallocate

  integer,parameter :: MSG_LEN=500


!!****t* m_malloc/minfo_t
!! NAME
!!  minfo_t
!! 
!! FUNCTION
!!  Store information on the memory allocated at run-time
!! 
!! SOURCE

type :: minfo_t

  integer(kind=i8b) :: tot_memory=0 
    ! Total memory allocated in bytes

  integer(kind=i8b) :: peak=0   
    ! Maximum memory allocated in bytes

  integer :: line=-1
    ! Line number in filename where the allocation of peak was performed.

  character(len=MSG_LEN) :: func_name="Unknown" 
    ! Name of the procedure where the allocation was performed.

  character(len=MSG_LEN) :: file_name="Unknown"
    ! Filename containing the procedure 

  character(len=MSG_LEN) :: array="_Unknown"
    ! Name of the array

end type minfo_t
!***

! Global structure used to store the total amount of memory, and the position of the peak
 type(minfo_t),private,save :: MINFO  

!integer :: LEVEL=2
!real,save :: memorylimit_abi = 0.e0
!logical,save :: meminit_abi = .false.
!integer, parameter :: mallocFile = 98
!type(memstat),save :: memloc_abi,memtot_abi
!integer,save :: memalloc_abi,memdealloc_abi,memproc_abi = 0
!!Debug option for memocc_abi, set in the input file
!!logical :: memdebug=.true.
!integer,save :: malloc_level_abi=2

 integer,save :: MALLOC_LEVEL=0

 ! Output file.
 ! Use same unit as the one employed in m_profiling_abi
 ! This should become a reserved unit to avoid problems with the rest of the code!
 character(len=MSG_LEN),save :: MALLOC_FILE="__malloc__.txt"
 integer,parameter :: MALLOC_FUNIT = 98 

 ! Counters: number of allocations, deallocations
 integer,save :: NUM_ALLOC=0, NUM_DEALLOC=0

 ! Selective memory tracing
#define NONE_STRING "__NONE_STRING__"
 character(MSG_LEN),save :: FILTER_FILENAME = NONE_STRING 
 character(MSG_LEN),save :: FILTER_FUNCNAME = NONE_STRING

 integer,save :: MYRANK   ! Rank of processor in MPI_COMM_WORLD
 integer,save :: NPROCS   ! Number of processors in MPI_COMM_WORLD

contains 
!!***

!!****f* m_malloc/malloc_init
!! NAME
!!  malloc_init
!!
!! FUNCTION
!!  Initialize the module
!!
!! INPUT
!!  level:
!!      Integer selecting the operation mode:
!!         0 no file malloc.prc is created, only memory allocation counters running
!!         1 file malloc.prc is created in a light version (only information on the memory peak is written)
!!         2 file malloc.prc is created with full information inside (default state if not specified)
!!      The status can only be downgraded. A stop signal is produced if status is increased
!!
!!  mpimode="master", "personal"
!!     Used only if level>0. 
!!     If "master" only master node i.e. rank==0 in comm_world will write data to malloc.prc
!!     If "personal" each node will write data to its own malloc_PROC#.prc file.
!!
!! PARENTS
!!
!! CHILDREN
!!      malloc_abort
!!
!! SOURCE

subroutine malloc_init(level, mpimode, file_name, func_name)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'malloc_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: level
 character(len=*),intent(in) :: mpimode
 character(len=*),optional,intent(in) :: file_name, func_name

!Local variables-------------------------------
 integer :: iostat,ierr
 
! *************************************************************************

 MALLOC_LEVEL = level

 MYRANK=0; NPROCS=1; ierr=0
#ifdef HAVE_MPI
 call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,ierr)
#endif

 if (present(file_name)) FILTER_FILENAME = file_name
 if (present(func_name)) FILTER_FUNCNAME = func_name

 if (level> 0) then ! Open file.
   
   if (mpimode == "master") then
     if (MYRANK==0) then
       ! Only master writes malloc info to file.
       open(unit=MALLOC_FUNIT, file=MALLOC_FILE, form="formatted", iostat=iostat)
     else
       MALLOC_FILE = NONE_STRING
     end if
   else
     ! Each node writes malloc info to its own file.
     write(MALLOC_FUNIT,"(a,i0,a)")"__malloc__P",MYRANK,".txt"
     open(unit=MALLOC_FUNIT, file=MALLOC_FILE, form="formatted", iostat=iostat)
   end if

   if (iostat /= 0) then 
!    MT feb 2014: split in several lines because __FILE__ can be expanded as a long line
     call malloc_abort(iostat,"Opening MALLOC_FILE",&
&     __FILE__,&
&     "malloc_init",__LINE__)
   end if

 end if

end subroutine malloc_init
!!***

!----------------------------------------------------------------------

!!****f* m_malloc/malloc_shutdown
!! NAME
!!  malloc_init
!!
!! FUNCTION
!!  Finalize the memory tracing.
!!
!! INPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine malloc_shutdown()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'malloc_shutdown'
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
! integer :: ierr                                     
 
! *************************************************************************

 if (NPROCS > 1) then
   ! Collect important info from the other MPI nodes.
   ! TODO
 end if

 ! Write final report to std_out
 !call malloc_report(unit=std_out)
 !if (MYRANK == 0) then
 !  call malloc_report(unit=ab_out)
 !end if

 ! Close file.
 if (MALLOC_FILE /= NONE_STRING) then
    close(unit=MALLOC_FUNIT)
 end if 

end subroutine malloc_shutdown
!!***

!----------------------------------------------------------------------

!!****f* m_malloc/malloc_report
!! NAME
!!  malloc_report
!!
!! FUNCTION
!!   Write memory consumption report on the specified unit.
!!
!! INPUT
!!   Fortran logical unit (already open in formatted mode)
!!
!! PARENTS
!!
!! CHILDREN
!!      malloc_abort
!!
!! SOURCE

subroutine malloc_report(unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'malloc_report'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unit

! *************************************************************************

 ! Write report
 write(unit,'(1x,a)')'-------------------------MEMORY CONSUMPTION REPORT-----------------------------'
 write(unit,'(1x,2(i0,a,1x),i0)')NUM_ALLOC,' allocations and ',NUM_DEALLOC,' deallocations, remaining memory(B): ',MINFO%tot_memory
 write(unit,'(1x,a,i0,a)') 'memory occupation peak: ',MINFO%peak/int(1024**2,kind=8),' Mb'
 write(unit,'(4(1x,a))')   'for the array ',trim(MINFO%array),' in procedure ',trim(MINFO%func_name)

end subroutine malloc_report
!!***

!----------------------------------------------------------------------

!!****f* m_malloc/malloc_abort
!! NAME
!!  malloc_abort
!!
!! FUNCTION
!!  Stop the code if an error occurs.
!!
!! INPUT
!!  istat=Status error of the parent.
!!  msg=Error message 
!!  file_name=File name
!!  func_name=Function name.
!!  line=Line number
!!
!! PARENTS
!!      m_malloc
!!
!! CHILDREN
!!      malloc_abort
!!
!! SOURCE

subroutine malloc_abort(istat,msg,file_name,func_name,line)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'malloc_abort'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: istat,line
 character(len=*),intent(in) :: msg,file_name,func_name

!Local variables-------------------------------
 integer :: ierr                                     
 
! *************************************************************************

 write(std_out,*)istat,msg,file_name,func_name,line

 ierr = 0
#ifdef HAVE_MPI
 call MPI_ABORT(MPI_COMM_WORLD,MPI_ERR_UNKNOWN,ierr)
#endif
 stop

end subroutine malloc_abort
!!***

!----------------------------------------------------------------------

!!****f* m_malloc/malloc_register
!! NAME
!!  malloc_register
!!
!! FUNCTION
!!  Register an allocation/deallocation event
!!
!! INPUT
!!  istat=Status error of the parent who has performed an (allocation/deallocation)
!!  mem_size=Memory size in bytes (>=0 for allocation, < 0 for deallocations)
!!  arr_name=Array name
!!  file_name=File name
!!  func_name=Function name.
!!  line=Line number
!!
!! PARENTS
!!
!! CHILDREN
!!      malloc_abort
!!
!! SOURCE

subroutine malloc_register(istat,mem_size,arr_name,file_name,func_name,line)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'malloc_register'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: istat,line
 integer(i8b),intent(in) :: mem_size
 character(len=*),intent(in) :: arr_name,file_name,func_name

!Local variables-------------------------------
 !integer :: ierr                                     
 logical :: has_peak,do_log
 character(len=MSG_LEN) :: msg

! *************************************************************************

 MEM_STAT = istat

 ! Handle allocate/deallocate failures
 if (istat /= 0) then
   write(msg,('(5a,i0,/,3a,i0,a)'))&
&    'Procedure: ',trim(func_name),"@",trim(file_name),":",line,&
&    'problem of allocation of array: ',trim(arr_name),'. Error code = ',istat,' Aborting now...'
   call malloc_abort(istat,msg,file_name,func_name,line)
  end if 

 ! Selective memory tracing
 do_log = .True.
 if (FILTER_FILENAME /= NONE_STRING) then
   do_log = (FILTER_FILENAME == file_name)
 end if

 if (FILTER_FUNCNAME /= NONE_STRING) then
   do_log = do_log .and. (FILTER_FUNCNAME == func_name)
 end if

 ! Register allocation/deallocation
 if (mem_size<0) then
   ! Deallocation
   NUM_DEALLOC = NUM_DEALLOC + 1
 else
   ! Allocation
   NUM_ALLOC = NUM_ALLOC + 1

   MINFO%tot_memory = MINFO%tot_memory + mem_size

   has_peak = .False.
   if (mem_size > MINFO%peak) then
      ! Update info.
      has_peak = .True.
      MINFO%peak = mem_size
      MINFO%func_name=func_name
      MINFO%file_name=file_name
      MINFO%line=line
   end if

   if (MALLOC_LEVEL == 1 .and. has_peak .and. do_log) then
     ! Write only if we have a peak in memory allocation.
     if (MALLOC_FILE /= NONE_STRING) then
       rewind(unit=MALLOC_FUNIT)
       !write(MALLOC_UNIT,*)"..."
     end if
   end if

   if (MALLOC_LEVEL == 2 .and. do_log) then
     ! Append mode. Register everything.
     if (MALLOC_FILE /= NONE_STRING) then
       !write(MALLOC_FUNIT,*)"..."
     end if
   end if

 end if 

end subroutine malloc_register
!!***

!----------------------------------------------------------------------

! Include the routines generated by genmalloc.py
#include "malloc.finc"

end module m_malloc
!!***
