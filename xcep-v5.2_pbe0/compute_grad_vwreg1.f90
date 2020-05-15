!====================================================
!
! gradient of total energy with respect to the phi, 
!
! note that rho depends on the phi as
!
! rho = phi**2/<phi,phi>*N
!
!====================================================
subroutine compute_grad_vwreg1(nfft,n1,n2,n3,qvec,dvol,phi,q,vext,etotal,dE_dphi,ef)

  use comm_data 

  implicit none 

  integer :: nfft,n1,n2,n3
  real(8) :: phi(nfft),dE_dphi(nfft),dvol,vext(nfft), & 
             dd,rho(nfft),q,f(nfft), & 
             qvec(3,(n1/2+1)*n2*n3), ef, & 
             vOF(nfft), tf_pot(nfft), etotal, & 
             tf_energy, vw_pot(nfft), vw_energy

  ! compute electron density 
  dd = sum(phi**2*dvol)
  f = phi/sqrt(dd)
  rho = f**2*q

  call tf(nfft,rho,dvol,tf_pot,tf_energy)
  call vw(nfft,n1,n2,n3,rho,qvec,dvol,vw_reg,vw_pot,vw_energy)

  vOF = tf_pot + vext + vw_pot*vw_lam 
  ef  = dvol*sum(vOF*rho)/q

  ! total energy 
  etotal = vw_lam*vw_energy + tf_energy + dvol*sum(vext*rho)

  ! gradient 
  dE_dphi = 2.d0*q*vOF*f/sqrt(dd) - 2.d0*q*sum(dvol*vOF*f*phi)*phi/dd**(1.5d0)


end subroutine compute_grad_vwreg1

