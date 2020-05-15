!======================================
! compute total energy and dE/dphi
!
! electron denisty is regularized as 
!
!  rho_reg = (1-exp(-sqrt(rho)/c))*rho
!======================================
subroutine compute_etotal_vwreg3(nfft,n1,n2,n3,qvec,dvol,phi,q,vext,etotal,dE_dphi,ef)
  use comm_data
  implicit none 
  integer :: nfft, n1,n2,n3
  real(8) :: phi(nfft),etotal,dvol,vext(nfft), & 
             dd,rho(nfft), q, f(nfft), & 
             qvec(3,(n1/2+1)*n2*n3), & 
             expo(nfft), & 
             tf_pot(nfft), tf_energy,  & 
             rho_reg(nfft), vwreg_phi, vw_pot(nfft), vw_energy,  & 
             vOF(nfft), ef, lap(nfft), dE_dphi(nfft), & 
             y(nfft)


  ! threshold for phi
  ! this is for reg_type = 3
  vwreg_phi = vw_reg**0.5

  ! compute electron density 
  dd = sum(phi**2*dvol)  
  f = phi/sqrt(dd)
  rho = f**2*q

  call tf(nfft,rho,dvol,tf_pot,tf_energy)

  ! compute vw energy 
  expo    = sqrt(rho)/vwreg_phi
  rho_reg = (1.d0-exp(-expo))*rho
  call laplacian(nfft,n1,n2,n3,sqrt(rho_reg),qvec,lap)
  vw_energy = sum(-0.5d0*sqrt(rho_reg)*lap)*dvol

  ! total energy 
  etotal = vw_lam*vw_energy + tf_energy + dvol*sum(vext*rho)


  !======================================
  ! compute dE/dphi
  !======================================

  ! d(sqrt(rho_reg))/d(rho)
  y = (1.d0 - exp(-expo) + exp(-expo)/2.d0/vwreg_phi*sqrt(rho)) / & 
       2.d0/sqrt(rho_reg + 1e-12) ! prevent rho_reg becomes zero

  ! dE/drho
  vOF = tf_pot + vext + vw_lam*(-lap)*y
  ef  = dvol*sum(vOF*rho)/q

  !
  ! since rho depends on the norm of phi. 
  ! d rho / d phi is a matrix and is a little tricky. 
  ! See paper for the derivations 
  !
  dE_dphi = vOF*2.d0*phi*q/dd - 2.d0*phi*q/dd**2*sum(dvol*vOF*phi**2)

end subroutine compute_etotal_vwreg3


