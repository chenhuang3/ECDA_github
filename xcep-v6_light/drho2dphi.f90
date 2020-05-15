!
!   Solve delta_rho for given delta_phi using conjugate gradiemnt method 
!   Rho is a function of phi
!
!   rho = phi**2*Ne/<phi,phi>
! 
!   we then have 
!
!     delta_rho(r) = \int \frac{\delta rho(r)}{\delta phi(r1)} delta_phi(r1) dr_1
!   
!   INPUT:   drho -- desired change of rho 
!            phi  -- current phi
!
!   OUTPUT:  dphi -- for making the desired delta_rho
!
!   create by Chen Huang (12/15/2019)
!
subroutine   drho2dphi(file_unit,nfft,dvol,q,phi,drho,dphi)
  implicit none 

  integer :: nfft, file_unit
  real(8) :: drho(nfft), dphi(nfft), & 
             dvol, q, phi(nfft)

  ! local vars =============
  integer :: iter
  real(8) :: arr_tmp(nfft), new_dphi(nfft), norm


  ! >>>>>>>>>>>>>> function <<<<<<<<<<<<< !
  write(file_unit,'(a)')NEW_LINE('a')//'enter drho2dphi()'

  !=============================
  ! iterative solver for dphi 
  !=============================
  norm = sum(dvol*phi**2)
  dphi = 0.d0 
  do iter=1,100
!    arr_tmp  = drho + 2.d0*phi**2*q/norm**2*sum(dvol*phi*dphi)
!    new_dphi = arr_tmp / (2.d0*q/norm) / phi
    new_dphi = drho / 2.d0 / phi
    write(file_unit,*)'resid: ',sum((new_dphi-dphi)**2*dvol),'  dphi:',minval(dphi),maxval(dphi)
    dphi = (new_dphi + dphi)/2.d0
  enddo

end subroutine
