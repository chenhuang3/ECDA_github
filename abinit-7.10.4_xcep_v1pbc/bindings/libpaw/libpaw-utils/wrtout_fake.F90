!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrtout_fake
!! NAME
!!  wrtout_fake
!!
!! FUNCTION
!!  Fake wrtout function.
!!  Print out a message.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2014 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  msg=(character(len=*)) message to be written
!!  unit=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!  [mode_paral]= --optional argument-- ignored in this fake version
!!
!! OUTPUT
!!  (only writing)
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wrtout(unit,msg,mode_paral)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=*),optional,intent(in) :: mode_paral
 character(len=*),intent(in) :: msg

!******************************************************************

 write(unit,fmt='(a)') msg

end subroutine wrtout
!!***
