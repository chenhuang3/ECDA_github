!{\src2tex{textfont=tt}}
!!****f* ABINIT/iofn1
!! NAME
!! iofn1
!!
!! FUNCTION
!! Begin by eventual redefinition of unit std_in and std_out
!! Then, print greetings for interactive user.
!! Next, Read filenames from unit std_in, AND check that new
!! output file does not already exist.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  character(len=fnlen) :: filnam(5)=character strings giving file names
!!  character(len=fnlen) :: filstat=character strings giving name of status file
!!
!! NOTES
!! If it does exist, isfile will create a new name to avoid overwriting the output file.
!! Also create name of status file
!!
!! File names refer to following files, in order:
!!  (1) Formatted input file  (std_in)
!!  (2) Formatted output file (std_out)
!!  (3) Root name for generic input files (wavefunctions, potential, density ...)
!!  (4) Root name for generic output files (wavefunctions, potential, density,
!!                                          DOS, hessian ...)
!!  (5) Root name for generic temporary files (wftmp1,wftmp2,kgunit,status ...)
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      abi_log_status_state,int2char4,isfile,xmpi_barrier,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine iofn1(filnam,filstat,comm)

 use defs_basis
 use m_profiling_abi
 use m_xmpi
 use m_errors

 use m_fstrings, only : int2char4
 use m_io_tools, only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iofn1'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: comm
 character(len=fnlen), intent(out) :: filstat
 character(len=fnlen), intent(out) :: filnam(5)

!Local variables-------------------------------
 character(len=1) :: blank
 integer,parameter :: master=0
 integer :: me,ios,nproc,ierr
 logical :: ex
 character(len=fnlen) :: fillog,tmpfil
 character(len=10) :: tag
 character(len=500) :: message

!*************************************************************************

!Initialise (beautification MS)
 blank = ' '; tmpfil = ''

!Determine who I am in comm
 me = xcomm_rank(comm)
 nproc = xcomm_size(comm)

!Define values of do_write_log and do_write_status parameters
!if a _NOLOG file exists no LOG file and no STATUS file are created for each cpu core
!if a _LOG file exists, a LOG file and a STATUS file are created for each cpu core
!if the #_of_cpu_core>NPROC_NO_EXTRA_LOG, LOG file is only created for master proc
!if the #_of_cpu_core>NPROC_NO_EXTRA_STATUS, STATUS file is only created for master proc
 inquire(file=ABI_NO_LOG_FILE,iostat=ios,exist=ex)
 if (ios/=0) ex=.false.
 if (ex) then
   call abi_log_status_state(new_do_write_log=.false.,new_do_write_status=.false.)
 else
   inquire(file=ABI_ENFORCE_LOG_FILE,iostat=ios,exist=ex)
   if (ios/=0) ex=.false.
   if (ex) then
     call abi_log_status_state(new_do_write_log=.true.,new_do_write_status=.true.)
   else
     if (do_write_log.and.me/=0) then
       call abi_log_status_state(new_do_write_log=(nproc<NPROC_NO_EXTRA_LOG))
     end if
     if (do_write_status.and.me/=0) then
       call abi_log_status_state(new_do_write_status=(nproc<NPROC_NO_EXTRA_STATUS))
     end if
   end if
 end if

 if (me==master) then

!  Eventually redefine standard input and standard output

   if (do_write_log) then
#if defined READ_FROM_FILE
!    Take care of the output file
     tmpfil(1:fnlen)=blank
     tmpfil(1:3)='log'
     call isfile(tmpfil,'new')
     close(std_out)
     if (open_file(tmpfil,message,unit=std_out,form='formatted',status='new') /= 0) then
       MSG_ERROR(message)
     end if 
#endif
   else
!    Redirect standard output to null
     close(std_out)
     if (open_file(NULL_FILE,message,unit=std_out) /= 0) then
       MSG_ERROR(message)
     end if
   end if

#if defined READ_FROM_FILE
!  Now take care of the "files" file
   tmpfil(1:fnlen)=blank
   tmpfil(1:9)='ab.files'
   write(message, '(4a)' )&
&   '  Because of cpp option READ_FROM_FILE,',ch10,&
&   '  read file "ab.files" instead of standard input ' ,ch10
   MSG_COMMENT(message)
   call isfile(tmpfil,'old')
   close(std_in)
   if (open_file(tmpfil,message,unit=std_in,form='formatted',status='old') /= 0) then
     MSG_ERROR(message)
   end if
#endif

!  Print greetings for interactive user
   write(std_out,*)' ABINIT '
   write(std_out,*)' '

!  Read name of input file (std_in):
   write(std_out,*)' Give name for formatted input file: '
   read(std_in, '(a)' ) filnam(1)
   write(std_out, '(a)' ) trim(filnam(1))
   write(std_out,*)' Give name for formatted output file:'
   read (std_in, '(a)' ) filnam(2)
   write (std_out, '(a)' ) trim(filnam(2))
   write(std_out,*)' Give root name for generic input files:'
   read (std_in, '(a)' ) filnam(3)
   write (std_out, '(a)' ) trim(filnam(3))
   write(std_out,*)' Give root name for generic output files:'
   read (std_in, '(a)' ) filnam(4)
   write (std_out, '(a)' ) trim(filnam(4))
   write(std_out,*)' Give root name for generic temporary files:'
   read (std_in, '(a)' ) filnam(5)
   write (std_out, '(a)' ) trim(filnam(5))

!  Check that old input file exists
   call isfile(filnam(1),'old')

!  Check that new output file does NOT exist
   call isfile(filnam(2),'new')

!  Check that root name for generic input and output differ
   if ( trim(filnam(3))==trim(filnam(4)) ) then
     write(message, '(a,a,a)' )&
&     'Root name for generic input and output files must differ ',ch10,&
&     'Action : correct your "file" file.'
     MSG_ERROR(message)
   end if

!  Check that root names are at least 20 characters less than fnlen
   if ( len_trim(filnam(3)) >= (fnlen-20) ) then
     write(message, '(a,a,a,a,a,i4,a,i4,a,a)' )&
&     'Root name for generic input files is too long. ',ch10,&
&     'It must be 20 characters less than the maximal allowed ',ch10,&
&     'length of names, that is ',fnlen,', while it is ',len_trim(filnam(3)),ch10,&
&     'Action : correct your "file" file.'
     MSG_ERROR(message)
   end if
   if ( len_trim(filnam(4)) >= (fnlen-20) ) then
     write(message, '(a,a,a,a,a,i4,a,i4,a,a)' )&
&     'Root name for generic output files is too long. ',ch10,&
&     'It must be 20 characters less than the maximal allowed ',ch10,&
&     'length of names, that is ',fnlen,', while it is ',len_trim(filnam(4)),ch10,&
&     'Action : correct your "file" file.'
     MSG_ERROR(message)
   end if
   if ( len_trim(filnam(5)) >= (fnlen-20) ) then
     write(message, '(a,a,a,a,a,i4,a,i4,a,a)' )&
&     'Root name for generic temporary files is too long. ',ch10,&
&     'It must be 20 characters less than the maximal allowed ',ch10,&
&     'length of names, that is ',fnlen,', while it is ',len_trim(filnam(5)),ch10,&
&     'Action : correct your "file" file.'
     MSG_ERROR(message)
   end if

 end if ! master only

!Communicate filenames to all processors
 call xmpi_bcast(filnam,master,comm,ierr)

!Check
!Create a name for the status file, based on filnam(5)
 filstat=trim(filnam(5))//'_STATUS'

!Redefine the log unit if not the master
 if (me/=master) then
   call int2char4(me,tag)
   ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
   filstat=trim(filstat)//'_P-'//trim(tag)
   if (do_write_log) then
     fillog=trim(filnam(5))//'_LOG_'//trim(tag)
     close(std_out)
     if (open_file(fillog,message,unit=std_out,status='unknown') /= 0) then
       MSG_ERROR(message)
     end if
   else
     close(std_out)
     if (open_file(NULL_FILE,message,unit=std_out) /= 0) then
       MSG_ERROR(message)
     end if
   end if
 end if

 call xmpi_barrier(comm)

end subroutine iofn1
!!***
