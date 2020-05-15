
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! compute dE/d\sqrt(rho)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine  compute_OF_g(nfft,n1,n2,n3,nspin,vw_lam,qvec,ucvol,phi,v,dvol,g)
  use comm_data, only: vw_reg
  implicit none 
  integer :: nfft,n1,n2,n3,nspin
  real(8) :: g(nfft), vw_lam, & 
             cTF = 3.d0/10.d0*(3.d0*3.1415926d0**2)**(2.d0/3.d0), & 
             phi(nfft), & 
             phi_reg(nfft), & 
             ucvol, & 
             qvec(3,(n1/2+1)*n2*n3), & 
             lap(nfft), & 
             v(nfft), eh, vh(nfft), & 
             dvol,exc,vw_pot(nfft),vw_energy

  !call calculate_xc_pw92lsda(n1,n2,n3,nspin,phi**2,dvol,vxc,Exc)
  !call hartree(n1,n2,n3,nspin,ucvol,qvec,phi**2,vh,eh)

  if (nspin==2) then 
    print *,'OF_cgmin() only works with non-spin polarized case, now stop!'
    stop
  endif 

  phi_reg = sqrt(phi**2 + vw_reg)

  call laplacian(nfft,n1,n2,n3,phi_reg,qvec,lap)

  g =   10.d0/3.d0*cTF*(phi**2)**(2.d0/3.d0)*phi &  ! TF term 
      - vw_lam*lap*phi/phi_reg & 
      + 2.d0*phi*v                   ! external term

end subroutine compute_OF_g
