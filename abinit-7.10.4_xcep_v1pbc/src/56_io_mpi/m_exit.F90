!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_exit
!! NAME
!! m_exit
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2013 ABINIT group (MG, DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_exit

 use defs_basis
 use m_xmpi
 use m_errors

 use m_time,      only : abi_wtime
 use m_io_tools,  only : get_unit

 implicit none

 private 
!!***

 !public :: exit_setup  ! Initialize the global variables of the module
 public :: exit_check     ! Test if we should try to stop the code gracefully 
                          ! and to create a restart point.

! Global variables.
 integer,private,save :: MODULE_INIT=0
! 0 if the module has not been initialized, 1 otherwise

 real(dp),private,save :: WALL0          
! Origin of time. 

 real(dp),private,save :: WALLTIME_LIMIT 
! Wall time limit in seconds.

 logical,private,save :: DOTEST_ABIEXIT
! True if we have to test for the presence of STOPFILE in the working directory

 character(len=fnlen),parameter :: STOPFILE = "abinit.exit"
! Name of the file that will trigger the exit if DOTEST_ABIEXIT

!----------------------------------------------------------------------

CONTAINS  
!!***

!!****f* m_exit/exit_setup
!! NAME
!! exit_setup
!!
!! FUNCTION
!!  Initialize the global variables of the modules.
!!  This is a collective function that should be called by all the 
!!  nodes in COMM_WORLD
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      inupper,timein,wrtout,xmpi_bcast
!!
!! SOURCE

subroutine exit_setup(wtime_limit,test_abiexit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exit_setup'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 real(dp),intent(in) :: wtime_limit
 logical,intent(in) :: test_abiexit

! *************************************************************************
 
 MODULE_INIT = MODULE_INIT + 1

 DOTEST_ABIEXIT = test_abiexit
 WALLTIME_LIMIT = wtime_limit
 WALL0 = abi_wtime()

end subroutine exit_setup
!!***

!----------------------------------------------------------------------

!!****f* m_exit/exit_check
!! NAME
!! exit_check
!!
!! FUNCTION
!! This routine checks whether the CPU time limit is exceeded or not.
!! If openexit is non-zero, it also checks the "filename" file
!! for the "exit" character string in its first line and returns the location
!! of the string on the line (0 if not found).  Maps both strings to upper case
!! before attempting to match them. Also checks for the existence
!! of the "abinit.exit" file in the directory where the job was started.
!! Finally, checks whether the CPU time limit was not exceeded.
!! If one of these conditions occurs, will induce graceful exit of iterations.
!!
!! INPUTS
!!  cpus = CPU time limit
!!  filename = character string giving name of file to be opened
!!  iout = unit number to print output to
!!  openexit = if 1, open the "filename" and "abinit.exit" files
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  iexit = index of "exit" on first line of file (0 if not found),
!!      or -1 if the exit was ordered through the existence of the "exit" file
!!      or -2 if the exit was ordered through the CPU time limit.
!!
!! PARENTS
!!      driver,gstate,loper3,respfn,scprqt
!!
!! CHILDREN
!!      inupper,timein,wrtout,xmpi_bcast
!!
!! SOURCE

subroutine exit_check(cpus,filename,iexit,iout,comm,openexit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exit_check'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: comm
 real(dp),intent(in) :: cpus
 character(len=*),intent(in) :: filename
 integer,intent(in) :: openexit,iout
 integer,intent(out) :: iexit

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer,save :: iexit_save=0
 integer :: ierr,temp_unit,ierrmpi
 logical :: ex
 real(dp),save :: tcpu_last=zero
 character(len=500) :: message
 character(len=fnlen) :: line
 character(len=4), parameter :: string='EXIT'
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

 if (iexit_save==0) then   
   ! ABINIT will pass again in this routine even after exit call has been detected

   if (xcomm_rank(comm)==master) then
     ! Master tests and broadcast the result to others
     iexit=0

     ! Is it worth to test the cpu time ?
     if (abs(cpus)>1.0d-5 .or. openexit==1) then
       call timein(tsec(1),tsec(2))
     end if

     ! A first way of exiting: the cpu time limit
     if (abs(cpus)>1.0d-5) then
       if(cpus<tsec(1))iexit=-2
     end if

     ! Test the content of files only when sufficient time (2 sec) has elapsed from last time it was tested.
     if (openexit==1 .and. iexit==0 .and. tsec(1)-tcpu_last>two ) then
       ! TODO Remove this approach. Use abinit.exit!
       tcpu_last=tsec(1)
       ! Open file and read first line as character string
       temp_unit = get_unit()
       open(unit=temp_unit,file=filename,form='formatted',status='old')
       rewind (unit=temp_unit)
       read (unit=temp_unit,fmt='(a)',iostat=ierr) line
       if(ierr/=0)then
         write(message, '(a,a,a,i5,a,a)' )&
&         'Problem when reading file=',TRIM(filename),'iostat =',ierr,ch10,&
&         'Action : check whether this file is OK.'
         MSG_ERROR(message)
       end if
       ! Make a local copy of matching string of length equal to nonblank length of input string
       ! Map to upper case
       call inupper(line)
       iexit=index(line,string)
       close (unit=temp_unit)

       ! This is another way of exiting : the previous one does not work
       ! on some machines, may be because they keep a copy of the initial input file.
       if(iexit==0)then
         inquire(file='abinit.exit',exist=ex)
         if(ex)iexit=-1
       end if

     end if
   end if

   call xmpi_bcast(iexit,master,comm,ierrmpi)

 else 
   ! In case the exit mechanism has already been activated
   iexit=iexit_save
 end if

 if (iexit/=0) then
   if (iexit>0) write(message, '(a,a,a,a,a,a,a)' ) ch10,&
&   ' chkexi: WARNING -',ch10,&
&   '  Exit has been requested from file ',trim(filename),'.',ch10
   if (iexit==-1) write(message, '(a,a,a,a,a)' ) ch10,&
&   ' chkexi: WARNING -',ch10,&
&   '  Exit has been requested from file "abinit.exit".',ch10
   if (iexit==-2) write(message, '(a,a,a,a,a)' ) ch10,&
&   ' chkexi: WARNING -',ch10,&
&   '  Exit due to cpu time limit exceeded.',ch10
   if (iout/=std_out) then
     call wrtout(iout,message,'COLL')
   end if
   call wrtout(std_out,  message,'COLL')
 end if

 iexit_save=iexit

end subroutine exit_check
!!***

END MODULE m_exit
!!***
