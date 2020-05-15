!{\src2tex{textfont=tt}}
!!****f* ABINIT/isfile
!! NAME
!! isfile
!!
!! FUNCTION
!! Inquire Status of FILE
!! Checks that for status =
!! 'old': file already exists
!! 'new': file does not exist; if file exists,
!! filnam is modified to filnam.A or filnam.B,....
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! filnam=character string to specify filename
!! status='old' or 'new'
!!
!! OUTPUT
!! stops processing if old file does not exist; changes name
!! and returns new name in redefined filnam if new file already exists.
!!
!! PARENTS
!!      anaddb,iofn1,m_vcoul,mrgscr,ujdet
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine isfile(filnam,status)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isfile'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=3),intent(in) :: status
 character(len=fnlen),intent(inout) :: filnam

!Local variables-------------------------------
!scalars
 integer :: ii,ios,jj
 logical :: ex,found
 character(len=500) :: message
 character(len=fnlen) :: filnam_tmp
 character(len=fnlen) :: trialnam
!arrays
 character(len=1) :: alpha(27)

! *************************************************************************

 alpha(1:27)=(/' ','A','B','C','D','E','F','G','H','I','J','K','L','M','N',&
& 'O','P','Q','R','S','T','U','V','W','X','Y','Z'/)

 filnam_tmp=filnam

 if (status=='old') then !  Check that old file exists
   inquire(file=filnam,iostat=ios,exist=ex)

   if (ios/=0) then
     write(message,'(a,a,a,a,i8,a,a)')&
&     'Checks for existence of file  ',trim(filnam),ch10,&
&     'but INQUIRE statement returns error code',ios,ch10,&
&     'Action: identify which problem appears with this file.'
     MSG_ERROR(message)
   else if (.not.ex) then
     write(message, '(a,a,a,a,a)' )&
&     'Checks for existence of file  ',trim(filnam),ch10,&
&     'but INQUIRE finds file does not exist.',&
&     'Action: check file name and re-run.'
     MSG_ERROR(message)
   end if

 else if (status=='new') then

   found=.false.
   do jj=1,27
     filnam_tmp=trim(filnam)//alpha(jj)
     trialnam=trim(filnam_tmp)

     do ii=1,28
       if (jj>1 .and. ii==1) cycle

!      Check that new output file does NOT exist
       inquire(file=trim(trialnam),iostat=ios,exist=ex)

       if (ios/=0) then

!        There is a problem => stop
         write(message, '(a,a,a,a,i8,a,a)' )&
&         'Checks for existence of file  ',trim(trialnam),ch10,&
&         'but INQUIRE statement returns error code',ios,ch10,&
&         'Action: identify which problem appears with this file.'
         MSG_ERROR(message)

       else if (ex) then

         write(message,'(4a)')'Finds that output file ',trim(trialnam),ch10,' already exists.'
         MSG_WARNING(message)
!        'New' file already exists; define a new file name
         if (jj==27 .and. ii==28) then
           write(message,'(a,a,a)')&
&           'Have used up all names of the form filename.[A-Z]{,2}',ch10,&
&           'Action: clean up your directory and start over.'
           MSG_ERROR(message)
         else if ( jj <27 .and. ii==28) then
           cycle
         else
           trialnam=trim(filnam_tmp)//alpha(ii)
           write(message, '(a,a,a)' ) ' new name assigned:',trim(trialnam),ch10
           call wrtout(std_out,message,'PERS')
           cycle
         end if

       else ! The name (or the new name) is correct
         found=.true.
         exit
       end if

!      End loop on ii : scan the alphabet.
!      There is a "cycle" and an "exit" in the loop.
     end do
     if (found) exit
   end do

   filnam=trim(trialnam)

 else ! status not recognized
   write(message,'(3a)')'  Input status= ',status,' not recognized.'
   MSG_BUG(message)
 end if

end subroutine isfile
!!***
