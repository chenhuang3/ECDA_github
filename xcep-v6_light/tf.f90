!=========================================
! compute TF kinetic energy and potential
! with incoming electron density
!
!  Chen Huang  Jan/2011 
!  
!=========================================
subroutine TF(nfft,rhor,dvol,tf_pot,tf_energy)
 implicit none

 integer,intent(in) :: nfft
 real(kind=8),intent(in ) :: rhor(nfft),dvol
 real(kind=8),intent(out) :: tf_pot(nfft),tf_energy

!! local vars
 real(kind=8) :: ctf = 2.87123400018819d0
 
! print *,'get in TF ...'

 tf_pot    = ctf * rhor**(2.0d0/3.d0) * (5.d0/3.d0)
 tf_energy = sum( ctf * rhor**(5.d0/3.d0) ) * dvol
 
! print *,'leave TF.'
 return 
end subroutine TF
