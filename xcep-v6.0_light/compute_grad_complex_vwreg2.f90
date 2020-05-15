!======================================
! compute total energy and dE/dphi
!
! electron denisty is regularized as 
!
!  rho_reg = (1-exp(-sqrt(rho)/c))*rho
!======================================
subroutine compute_grad_complex_vwreg2(nfft,n1,n2,n3,qvec,dvol,phi,q,vext,dE_dphi)

  use comm_data

  implicit none 

  integer :: nfft, n1,n2,n3
  complex(8) :: dd, phi(nfft), expo(nfft),  & 
                f(nfft), rho(nfft), tf_pot(nfft),  & 
                vOF(nfft), lap(nfft), & 
                rho_reg(nfft), dE_dphi(nfft), y(nfft)

  integer :: ii
  real(8) :: ctf = 2.87123400018819d0
  real(8) :: lap_re(nfft), lap_im(nfft), & 
             dvol, vext(nfft), & 
             q, qvec(3,(n1/2+1)*n2*n3)


  ! compute electron density 
  dd  = sum(phi**2*dvol) 
  f   = phi/sqrt(dd)
  rho = f**2*q

  ! TF potential 
  tf_pot = ctf*rho**(2.0d0/3.d0)*(5.d0/3.d0)

  !========= vW potential ==========
  expo    = rho/vw_reg
  rho_reg = (1.d0-exp(-expo))*rho

  ! real part 
  call laplacian(nfft,n1,n2,n3, real(sqrt(rho_reg)),qvec,lap_re)
  call laplacian(nfft,n1,n2,n3,aimag(sqrt(rho_reg)),qvec,lap_im)
  lap = cmplx(lap_re,lap_im)


  !======================================
  ! compute dE/dphi
  !======================================

  ! d(sqrt(rho_reg))/d(rho)
  y = (1.d0 - exp(-expo) + exp(-expo)/vw_reg*rho)/2.d0/sqrt(rho_reg+1e-20)

  ! dE/drho
  vOF = tf_pot + vext + vw_lam*(-lap)*y

  ! Since rho depends on the norm of phi. 
  ! d rho / d phi is a matrix and is a little tricky. 
  ! See paper for the derivations 
  dE_dphi = vOF*2.d0*phi*q/dd - 2.d0*phi*q/dd**2*sum(dvol*vOF*phi**2)

end subroutine


