!===========================================
! solve sqrt(rho) as a fixed point problem 
!===========================================
subroutine OF_fixpoint_solver(nfft,n1,n2,n3,ucvol,qvec,vext,rho,q)
   
  use comm_data 

  implicit none 

  integer,parameter :: npulay=5
  integer :: m, nfft, n1,n2,n3
  real(8) :: rho(nfft),  ucvol, dvol, & 
             qvec(3,(n1/2+1)*n2*n3), & 
             wf(nfft),  rho_new(nfft), & 
             wf_new(nfft), wf_reg(nfft), & 
             mu,mu_reg, & 
             q,vext(nfft),  & 
             tf_pot(nfft), tf_energy, & 
             etotal,lap(nfft), resid, & 
             pulay_beta = 0.5d0, & 
             pulay_work(nfft,npulay,2)


  write(logOF,*)
  write(logOF,'(a)')'>>> OF_fixpoint_solver() << '
  m = 0
  dvol = ucvol / dble(nfft)

  do while(.true.) 
    m = m + 1

    wf     = sqrt(rho)
    wf_reg = sqrt(rho + vw_reg)

    call compute_etotal_vwreg1(nfft,n1,n2,n3,qvec,dvol,wf,q,vext,etotal,mu)
    call laplacian(nfft,n1,n2,n3,wf_reg,qvec,lap)
    call tf(nfft,rho,dvol,tf_pot,tf_energy)

    wf_new = (vw_lam*(-0.5d0*lap)*wf/wf_reg + (tf_pot+vext)*wf)/mu
    wf_new = wf_new / sqrt(sum(dvol*wf_new**2)) * sqrt(q) ! renomalize 
    rho_new = wf_new**2


    ! check convergence 
    resid = sqrt(sum(dvol*(rho_new-rho)**2))
    write(logOF,'(a,i3,a,es12.4,a,2es12.4,a,f12.6,a,f14.8)') & 
      'fp_iter: ',m,' resid:',resid, & 
      ' rho:',minval(rho_new),maxval(rho_new),' ef:',mu, ' etotal:',etotal
    call flush(logOF)
    if (resid<1e-6) exit 

    rho = (rho+rho_new)/2.d0
    !call pulay_mix(nfft,m,rho,rho_new,pulay_work,npulay,pulay_beta)
    !rho = rho_new
  enddo 

  write(logOF,*)'OF_fixpoint_solver() converged.'
  write(logOF,*)
end subroutine 

