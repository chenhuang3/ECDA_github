
subroutine compute_OF_E(nfft,n1,n2,n3,nspin,vw_lam,qvec,ucvol,phi,v,dvol,etotal)
  use comm_data, only: vw_reg
  implicit none 
  integer :: nfft,n1,n2,n3,nspin
  real(8) :: etotal, ucvol, vw_lam, & 
             qvec(3,(n1/2+1)*n2*n3), & 
             cTF = 3.d0/10.d0*(3.d0*3.1415926d0**2)**(2.d0/3.d0), & 
             phi_reg(nfft), phi(nfft),  & 
             v(nfft), vh(nfft), eh, & 
             dvol, Exc, vw_pot(nfft), vw_energy

  
  call vw(nfft,n1,n2,n3,phi**2,qvec,dvol,vw_reg,vw_pot,vw_energy)

  etotal =   vw_energy*vw_lam  &    ! vW energy 
           + cTF*sum((phi**2)**(5.d0/3.d0))*dvol   &  ! TF energy 
           + sum(v*phi**2)*dvol                       ! external energy 

!  print *,'TF:  ',cTF*sum((phi**2)**(5.d0/3.d0))*dvol
!  print *,'vW:  ',-0.5d0*sum(phi*lap_phi)*dvol
!  print *,'sum(v*phi**2)*dvol: ',sum(v*phi**2)*dvol
!  stop
end subroutine compute_OF_E
