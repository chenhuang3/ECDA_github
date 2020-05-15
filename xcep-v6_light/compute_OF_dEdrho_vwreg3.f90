subroutine compute_OF_dEdrho_vwreg3(nfft,n1,n2,n3,qvec,ucvol,rho,vext,dE_drho) 

  use comm_data 

  implicit none 
  integer :: nfft, n1, n2, n3
  real(8) :: rho(nfft), dE_drho(nfft), vext(nfft), & 
             tf_pot(nfft), vw_pot(nfft), & 
             vwreg_phi, arr_tmp(nfft),  &
             rho_tmp(nfft), vtmp(nfft), ucvol, & 
             dvol, dtmp, lap(nfft),qvec(3,(n1/2+1)*n2*n3)


  dvol = ucvol/dble(nfft)
  dE_drho = vext

  ! ========= VW KEDF =============
  ! regularized rho  for vW KEDF 
  vwreg_phi = sqrt(vw_reg)
  arr_tmp = rho**0.5/vwreg_phi
  rho_tmp = (1.d0-exp(-arr_tmp))*rho

  call laplacian(nfft,n1,n2,n3,sqrt(rho_tmp),qvec,lap)

  ! compute d sqrt(rho')/d rho
  vtmp = (1.d0 - exp(-arr_tmp) + exp(-arr_tmp)*sqrt(rho)/2.d0/vwreg_phi)/&
          2.d0/sqrt(rho_tmp)

!  vtmp = sqrt(rho_tmp)/2.d0/rho & 
!       + sqrt(rho)/4.d0/vwreg_phi/sqrt(rho_tmp) & 
!       - sqrt(rho_tmp)/4.d0/vwreg_phi/sqrt(rho)

  dE_drho = dE_drho - vw_lam*lap*vtmp  ! vw potential 

  ! ========= TF KEDF ============
  call tf(nfft,rho,dvol,tf_pot,dtmp)
  dE_drho = dE_drho + tf_pot
end subroutine 

