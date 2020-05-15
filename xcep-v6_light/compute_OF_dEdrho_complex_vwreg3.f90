!
! dE/drho is of complex type for given complex density 
! 
subroutine compute_OF_dEdrho_complex_vwreg3(nfft,n1,n2,n3,qvec,ucvol,rho,vext,dE_drho) 

  use comm_data 

  implicit none 
  integer    :: nfft, n1, n2, n3
  complex(8) :: rho(nfft), dE_drho(nfft), & 
                tf_pot(nfft), vw_pot(nfft), & 
                arr_tmp(nfft), rho_tmp(nfft), & 
                vtmp(nfft), lap(nfft)

  real(8) :: vext(nfft), & 
             vwreg_phi, &
             ucvol, & 
             dvol, dtmp, & 
             lap_re(nfft), lap_im(nfft), & 
             qvec(3,(n1/2+1)*n2*n3)
             
  real(8) :: ctf = 2.87123400018819d0


  dvol = ucvol / dble(nfft)
  dE_drho = vext

  ! ========= VW KEDF =============

  ! regularized rho  for vW KEDF 
  vwreg_phi = sqrt(vw_reg)
  arr_tmp = sqrt(rho)/vwreg_phi
  rho_tmp = (1.d0-exp(-arr_tmp))*rho

  call laplacian(nfft,n1,n2,n3,real(sqrt(rho_tmp)), qvec,lap_re)
  call laplacian(nfft,n1,n2,n3,aimag(sqrt(rho_tmp)),qvec,lap_im)

  lap = cmplx(lap_re,lap_im)

  ! compute d sqrt(rho')/d rho
  vtmp = (1.d0 - exp(-arr_tmp) + exp(-arr_tmp)*sqrt(rho)/2.d0/vwreg_phi)/&
          2.d0/sqrt(rho_tmp)

  dE_drho = dE_drho - vw_lam*lap*vtmp  ! vw potential 

  tf_pot = ctf*rho**(2.0d0/3.d0)*(5.d0/3.d0) ! TF pot

  dE_drho = dE_drho + tf_pot

end subroutine 

