!{\src2tex{textfont=tt}}
!!****f* ABINIT/polyn_coeff
!! NAME
!!  polyn_coeff
!!
!! FUNCTION
!!  For N function values Y(X) compute coefficients of
!!  N-1 degree interpolating polynomial. 
!!  Due to G. Rybicki
!!  tested for linear and parabolic interpolation
!!
!! COPYRIGHT
!!  Copyright (C) 2014 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! n = number of points 
!! x = array of abcissa values
!! y = array of ordinate values
!!
!! OUTPUT
!! coeff(n) = array of polynomial coefficients
!!
!! PARENTS
!!      get_susd_null,geteexc_uc,prtsusd
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine polyn_coeff(n,x,y,coeff)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polyn_coeff'
!End of the abilint section

 implicit none

!arguments
 integer,intent(in) :: n
 real(dp),intent(in) ::  x(n),y(n)
 real(dp),intent(out) ::  coeff(n)

!local variables
 integer :: ii,jj,kk
 real(dp) :: acc,ff,den
 real(dp), allocatable :: s(:)

 ABI_ALLOCATE(s,(n))
 s(:)=0.0d0
 coeff(:)=0.0d0
 s(n)=-x(1)

 do ii=2,n
   do jj=n+1-ii,n-1
     s(jj)=s(jj)-x(ii)*s(jj+1)
   end do
   s(n)=s(n)-x(ii)
 end do

 do jj=1,n
   den=n
   do kk=n-1,1,-1
     den=kk*s(kk+1)+x(jj)*den
   end do
   ff=y(jj)/den
   acc=1.0d0
   do kk=n,1,-1
     coeff(kk)=coeff(kk)+acc*ff
     acc=s(kk)+x(jj)*acc
   enddo
 end do

 ABI_DEALLOCATE(s)

end subroutine polyn_coeff
!!***
