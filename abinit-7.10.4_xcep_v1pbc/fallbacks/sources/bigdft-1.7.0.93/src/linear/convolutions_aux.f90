!> @file
!! Convolutions for linear version (with quartic potentials)
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calculates the effective filter for the operator [kineticEnergy + (x-x0)^4].
subroutine getEffectiveFilterQuartic(parabPrefac,hgrid, x0, eff, filterCode)
use filterModule
implicit none

! Calling arguments
real(kind=8),intent(in) :: parabPrefac
real(kind=8),intent(in) :: hgrid                  !< Grid spacing
real(kind=8),intent(in) :: x0                     !< The center of the parabolic potential (x-x0)^2
real(kind=8),dimension(lb:ub),intent(out) :: eff  !< The effective filter
character(len=*), intent(in) :: filterCode

! Local variables
integer :: i
real(kind=8) :: fac, fac2, prefac1, hgrid2, hgrid3, x02, x03
prefac1=-.5d0/hgrid**2
fac=parabPrefac
fac2=parabPrefac*hgrid
hgrid2=hgrid**2
hgrid3=hgrid**3
x02=x0**2
x03=x0**3
! Determine which filter we have to calculate

if(filterCode=='a')then
    do i=lb,ub
        eff(i)=prefac1*a(i) + fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
    end do
    eff(0)=eff(0)+fac*x0**4
elseif(filterCode=='b') then
    do i=lb,ub
        eff(i)=prefac1*b(i) + fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
    end do
elseif(filterCode=='c') then
    do i=lb,ub
        eff(i)=prefac1*c(i) + fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
    end do
elseif(filterCode=='e') then
    do i=lb,ub
        eff(i)=prefac1*e(i) + fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
    end do
    eff(0)=eff(0)+fac*x0**4
else
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end if

end subroutine getEffectiveFilterQuartic


!> Calculates the effective filter for the operator (x-x0)^4
subroutine getFilterQuartic(parabPrefac,hgrid, x0, eff, filterCode)

use filterModule
implicit none

! Calling arguments
! Calling arguments
real(kind=8),intent(in) :: parabPrefac
real(kind=8),intent(in) :: hgrid                  !< Grid spacing
real(kind=8),intent(in) :: x0                     !< The center of the quartic potential (x-x0)^4
real(kind=8),dimension(lb:ub),intent(out) :: eff  !< The effective filter
character(len=*), intent(in) :: filterCode

! Local variables
integer :: i
real(kind=8) :: fac, fac2,  hgrid2, hgrid3, x02, x03
real(kind=8) :: scale
scale=1.d0
!scale=1.d-1
!scale=0.d-1
!scale=5.d-2
!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
fac=parabPrefac*scale
fac2=parabPrefac*hgrid*scale
hgrid2=hgrid**2
hgrid3=hgrid**3
x02=x0**2
x03=x0**3
! Determine which filter we have to calculate

if(filterCode=='a') then
    do i=lb,ub
        !eff(i)=prefac1*a(i) + fac2*(hgrid*a2(i)+2*x0*a1(i))
        eff(i) = fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
    end do
    !eff(0)=eff(0)+fac*x0**2
    eff(0)=eff(0)+fac*x0**4
elseif(filterCode=='b') then
    do i=lb,ub
        !eff(i)=prefac1*b(i) + fac2*(hgrid*b2(i)+2*x0*b1(i))
        eff(i) = fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
    end do
elseif(filterCode=='c') then
    do i=lb,ub
        !eff(i)=prefac1*c(i) + fac2*(hgrid*c2(i)+2*x0*c1(i))
        eff(i) = fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
    end do
elseif(filterCode=='e') then
    do i=lb,ub
        !eff(i)=prefac1*e(i) + fac2*(hgrid*e2(i)+2*x0*e1(i))
        eff(i) = fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
    end do
    !eff(0)=eff(0)+fac*x0**2
    eff(0)=eff(0)+fac*x0**4
else
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end if

end subroutine getFilterQuartic


!> Calculates the effective filter for the operator (x-x0)^2
subroutine getFilterQuadratic(parabPrefac,hgrid, x0, eff, filterCode)
use filterModule
implicit none

! Calling arguments
real(kind=8),intent(in) :: parabPrefac
real(kind=8),intent(in) :: hgrid                  !< Grid spacing
real(kind=8),intent(in) :: x0                     !< The center of the parabolic potential (x-x0)^2
real(kind=8),dimension(lb:ub),intent(out) :: eff  !< The effective filter
character(len=*), intent(in) :: filterCode


! Local variables
integer :: i
real(kind=8) :: fac, fac2, hgrid2, hgrid3, x02, x03
real(kind=8) :: scale
scale=1.d0
!scale=1.d-1
!scale=0.d-1
!scale=5.d-2
!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
fac=parabPrefac*scale
fac2=parabPrefac*hgrid*scale
hgrid2=hgrid**2
hgrid3=hgrid**3
x02=x0**2
x03=x0**3
! Determine which filter we have to calculate

if(filterCode=='a') then
    do i=lb,ub
        eff(i) = fac2*( hgrid*a2(i) + 2.d0*x0*a1(i) )
        !eff(i) = fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
    end do
    eff(0)=eff(0)+fac*x0**2
    !eff(0)=eff(0)+fac*x0**4
elseif(filterCode=='b') then
    do i=lb,ub
        eff(i) = fac2*( hgrid*b2(i) + 2.d0*x0*b1(i) )
        !eff(i) = fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
    end do
elseif(filterCode=='c') then
    do i=lb,ub
        eff(i) = fac2*( hgrid*c2(i) + 2.d0*x0*c1(i) )
        !eff(i) = fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
    end do
elseif(filterCode=='e') then
    do i=lb,ub
        eff(i) = fac2*( hgrid*e2(i) + 2.d0*x0*e1(i) )
        !eff(i) = fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
    end do
    eff(0)=eff(0)+fac*x0**2
    !eff(0)=eff(0)+fac*x0**4
else
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end if

end subroutine getFilterQuadratic


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_for_quartic_convolutions(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  & 
     mseg_c, mvctr_c, keyg_c, keyv_c,  & 
     mseg_f, mvctr_f, keyg_f, keyv_f,  & 
     scal, psi_c, psi_f, &
     work)
  use module_base
  use module_types
  implicit none
  integer,intent(in) :: n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, mseg_c, mvctr_c, mseg_f, mvctr_f
  integer,dimension(mseg_c),intent(in) :: keyv_c
  integer,dimension(mseg_f),intent(in) :: keyv_f
  integer,dimension(2,mseg_c),intent(in) :: keyg_c
  integer,dimension(2,mseg_f),intent(in) :: keyg_f
  real(wp),dimension(0:3),intent(in) :: scal
  real(wp),dimension(mvctr_c),intent(in) :: psi_c
  real(wp),dimension(7,mvctr_f),intent(in) :: psi_f
  type(workarrays_quartic_convolutions),intent(inout) :: work
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(private) &
  !$omp shared(scal) &
  !$omp shared(psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,mseg_c,mseg_f,work)
  !!! coarse part
  !$omp do
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
        work%xy_c(i2,i,i3)=psi_c(i-i0+jj)*scal(0)
        work%xz_c(i3,i,i2)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo
  !$omp enddo
  !!! fine part
  !$omp do
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_f1(i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xx_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xy_f(1,i2,i,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xz_f(1,i3,i,i2)=psi_f(1,i-i0+jj)*scal(1)

        work%xy_f2(i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xx_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xy_f(2,i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xz_f(2,i3,i,i2)=psi_f(2,i-i0+jj)*scal(1)

        work%xx_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xy_f(3,i2,i,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xz_f(3,i3,i,i2)=psi_f(3,i-i0+jj)*scal(2)

        work%xz_f4(i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)
        work%xx_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xy_f(4,i2,i,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xz_f(4,i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)

        work%xx_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xy_f(5,i2,i,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xz_f(5,i3,i,i2)=psi_f(5,i-i0+jj)*scal(2)

        work%xx_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xy_f(6,i2,i,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xz_f(6,i3,i,i2)=psi_f(6,i-i0+jj)*scal(2)

        work%xx_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xy_f(7,i2,i,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xz_f(7,i3,i,i2)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo
 !$omp enddo
 !$omp end parallel
 

END SUBROUTINE uncompress_for_quartic_convolutions
