!{\src2tex{textfont=tt}}
!!****f* ABINIT/rrho
!! NAME
!! rrho
!!
!! FUNCTION
!! Reads in the charge in mkdens3D format
!! The file was opened in the calling program, unit number 19.
!! The header was already read in the case of the unformatted file
!!
!! COPYRIGHT
!! Copyright (C) 2000-2014 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! densfileformat=integer flag :
!! 0=formatted (ASCII) or 1=unformatted (BINARY) or 2=unformatted NETCDF-ETSF file 
!! nr1=grid_full size along x
!! nr2=grid_full size along y
!! nr3=grid_full size along z
!! nspden=number of spin polartized densities (1 for non-spin polarized, 2 for spin-polarized)
!! wff=structure type containing the density information
!!
!! OUTPUT
!! grid_full(nr1,nr2,nr3)=grid_full matrix
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!      etsf_io_main_get,wffclose
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rrho(densfileformat,grid_full,nr1,nr2,nr3,nspden,wff)

 use defs_basis
 use m_errors
 use m_profiling_abi
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rrho'
!End of the abilint section

 implicit none

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: densfileformat,nr1,nr2,nr3,nspden
 type(wffile_type),intent(inout) :: wff
!arrays
 real(dp),intent(out),target :: grid_full(nr1,nr2,nr3,nspden)

!Local variables--------------------------------------------------------
!scalars
 character(len=500) :: message
 integer :: ierr,ir1,ir2,ir3,ispden
#if defined HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(etsf_main), target :: main_folder
 type(etsf_io_low_error) :: error
#endif

! *************************************************************************

 select case (densfileformat)

!  Formatted (only one spin component is allowed)
   case (0)
     do ir3=1,nr3
       do ir2=1,nr2
         do ir1=1,nr1
           read(unit=wff%unwff,fmt=*) grid_full(ir1,ir2,ir3,1)
         end do
       end do
     end do

!    Unformatted, on one record
   case (1)
     do ispden=1,nspden
       read(unit=wff%unwff) grid_full(1:nr1,1:nr2,1:nr3,ispden)
     end do
!    ETSF case
   case (2)
#ifdef HAVE_TRIO_ETSF_IO
     main_folder%density%data4D => grid_full
     call etsf_io_main_get(wff%unwff, main_folder, lstat, error)
     ETSF_CHECK_ERROR(lstat,error)
#else
     message = 'ETSF_IO support is not compiled. Reconfigure with --enable-etsf-io.'
     MSG_ERROR(message)
#endif

   case default
     write(message,'(a,i0)')"value for 3D function file format is invalid: ",densfileformat
     MSG_BUG(message)
 end select

 call wffclose(wff, ierr)

end subroutine rrho
!!***
