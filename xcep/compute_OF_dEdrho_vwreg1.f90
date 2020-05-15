  subroutine compute_OF_dEdrho_vwreg1(nfft,n1,n2,n3,qvec,ucvol,rho,vext,dEdrho)
     use comm_data 
     implicit none 
     integer :: nfft,n1,n2,n3
     real(8) :: rho(nfft), dEdrho(nfft), & 
                ucvol,  & 
                qvec(3,(n1/2+1)*n2*n3), & 
                dtmp, tfpot(nfft), vw_pot(nfft), & 
                vext(nfft), dvol

     dvol = ucvol / dble(nfft) 

     call vw(nfft,n1,n2,n3,rho,qvec,dvol,vw_reg,vw_pot,dtmp)
     call tf(nfft,rho,dvol,tfpot,dtmp)
     dEdrho = tfpot + vw_lam*vw_pot + vext
  end subroutine 
