!====================================================
! compute the total energy. 
! The electron density depends on the phi as
! rho = phi**2/<phi,phi>*N
!====================================================
subroutine compute_etotal_vwreg1(nfft,n1,n2,n3,qvec,dvol,phi,q,vext,etotal,ef)
  use comm_data 
  implicit none 
  integer :: nfft,n1,n2,n3
  real(8) :: phi(nfft),etotal,dvol,vext(nfft), & 
             qvec(3,(n1/2+1)*n2*n3), & 
             dd,rho(nfft), q, f(nfft), & 
             vOF(nfft), tf_pot(nfft), tf_energy,  & 
             vw_pot(nfft), vw_energy, ef

  ! compute electron density 
  dd = sum(phi**2*dvol)
  f = phi/sqrt(dd)
  rho = f**2*q

  call tf(nfft,rho,dvol,tf_pot,tf_energy)
  call vw(nfft,n1,n2,n3,rho,qvec,dvol,vw_reg,vw_pot,vw_energy)

  etotal = vw_lam*vw_energy + tf_energy + dvol*sum(vext*rho)

  vOF = tf_pot + vw_lam*vw_pot + vext
  ef  = dvol*sum(vOF*rho)/q

end subroutine 

