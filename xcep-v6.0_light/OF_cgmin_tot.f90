!
! Use conjugate gradient method to direct minimize total energy w.r.t. electron density 
! The basic variable is phi, the denisty is 
!
!  rho = phi**2/<phi,phi>*Ne
!
! therefore we perform unconstraied minimization
!
subroutine OF_cgmin_tot(n1,n2,n3,nfft,nspin,qvec,ucvol,vext,rho,chempot)
 
 use comm_data 

 implicit none 

 logical :: user_init

 integer :: n1,n2,n3, & 
            nfft, nspin, isp 

 real(8) :: vext(nfft,nspin), & 
            rho(nfft,nspin), & 
            rho_tmp(nfft), & 
            qvec(3,(n1/2+1)*n2*n3), & 
            resid, ucvol, dvol, chempot, & 
            q(nspin), q_tmp, & 
            coeff, phi(nfft), & 
            phi_tmp(nfft), & 
            lam, dconv, & 
            theta, & 
            dE_dphi(nfft), & 
            dE_dphi_tmp(nfft), & 
            g(nfft), & 
            g_old(nfft), & 
            d(nfft), & 
            cg_beta, & 
            cg_alpha

 logical :: print_info = .false., & 
            finished, & 
            do_precond 

 integer :: ii, opt_method, & 
            iter, mix_case, & 
            ipt,  array(1), cell_nfft(3)

 real(8) :: ctf = 2.87123400018819d0, & 
            vtf(nfft), vscf(nfft), mu, etotal, & 
            pre_ke, pi=3.1415926d0, & 
            etotal1, etotal2, etotal3, & 
            etotal_tmp, & 
            etotal_tmp2, & 
            dE_dtheta, theta0,  theta1, der1, der2, & 
            dtmp, pg(nfft),pg_old(nfft),lap(nfft), & 
            dd, f(nfft), wf(nfft), wf_reg(nfft), & 
            delta_rho(nfft),  & 
            delta_mu, & 
            wf_regN(nfft), arr_tmp(nfft)

 ! >>>>>>>>>>>>>> function <<<<<<<<<<<<!

 write(logOF,'(a)')NEW_LINE('a')//'OF_cgmin_tot() [revised on 4/2/2020] ...'

 dvol = ucvol / dble(nfft)
 cell_nfft(1) = n1
 cell_nfft(2) = n2
 cell_nfft(3) = n3

 do isp=1,nspin 
   q(isp) = dvol*sum(rho(:,isp))
   write(logOF,'(a,2es12.4,a,i3)')'vext: ',& 
     minval(vext(:,isp)),maxval(vext(:,isp)),' for spin:',isp
 enddo
 write(logOF,'(a,f12.5)')'q:        ',q(1:nspin)
 write(logOF,'(a,f12.4)')'vw_lambda:',vw_lam
 write(logOF,'(a,es12.4)')'vw_reg:   ',vw_reg
 write(logOF,'(a,i3)')'vwreg_type:',vwreg_type
 call flush(logOF)



 !===============================
 ! loop over spins 
 !===============================
 do isp=1,nspin 

   phi = sqrt(rho(:,isp))
   do_precond = .true. 
   pre_ke = 1.0  ! pre_ke for precondition
   theta = 0.0d0

!   ! randomize phi 
!   do ipt=1,nfft 
!     phi(ipt) = rand()
!   enddo 
!   phi = phi/sqrt(sum(phi**2*dvol)) * sqrt(q(isp))


   !=======================================
   ! cg method for minimize system energy 
   !=======================================
   iter = 0
   finished = .false. 
   write(logOF,'(a,f12.4,a)')'*** precond is on with pre_ke: ',pre_ke,' ***'

   do while(.true.)

     iter = iter + 1
     call flush(logOF)

     ! update electron density 
     dd = sum(phi**2*dvol)
     f  = phi/sqrt(dd)
     rho(:,isp) = f**2*q(isp)

     ! compute the gradient dE/dphi
     select case(vwreg_type) 
     case (1)  
       call compute_grad_vwreg1(nfft,n1,n2,n3,qvec,dvol,phi,q(isp), & 
                                vext(:,isp),etotal,dE_dphi,chempot)
     case (2)
       call compute_etotal_vwreg2(nfft,n1,n2,n3,qvec,dvol,phi,q(isp), & 
                                  vext(:,isp),etotal,dE_dphi,chempot)
     case (3)
       call compute_etotal_vwreg3(nfft,n1,n2,n3,qvec,dvol,phi,q(isp), & 
                                  vext(:,isp),etotal,dE_dphi,chempot)
     endselect

     g = dE_dphi

     
     ! print information 
     write(logOF,'(a,i4,a,f12.7,a,f9.5,a,es9.2,a,es8.2,a,es10.2)') & 
       'it: ',iter,'  etot:',etotal,'   Ef: ',chempot, & 
       '  theta: ',theta,'  |g|: ',sqrt(sum(g**2)*dvol),'  gmax:',maxval(abs(g))
     call flush(logOF)


     ! job fininshed? 
     if ( maxval(abs(g))<1e-3 .AND. iter>=20) then 
       write(logOF,'(a)')'gmax<1e-3 and iter>20. OF-DFT minimized for this spin'
       finished = .true.
     endif 
     if (iter>=100) then 
       write(logOF,*)'maxiter=100 is reached. exit OFDFT'
       finished = .true. 
     endif 
     if (finished) then 
       ! update electron density 
       dd = sum(phi**2*dvol)
       f = phi/sqrt(dd)
       rho(:,isp) = f**2*q(isp)
       call flush(logOF)
       exit 
     endif 


     ! precondition g => pg
     !================================
     if (do_precond) then 
       !if (vwreg_type==1) then  
       !  pg = g * sqrt(phi**2+vw_reg)*sum(dvol*phi**2)/q(isp)/(phi+vw_reg**0.5/10.d0)
       !endif 
       call precond_teter(n1,n2,n3,nfft,qvec,pre_ke,g,pg)
       !if (vwreg_type==1) then  
       !  pg = pg * sqrt(phi**2+vw_reg)*sum(dvol*phi**2)/q(isp)/(phi+vw_reg**0.5/10.d0)
       !endif 
     endif 


777 continue 

     ! compute CG direction 
     !=====================
     if (iter==1 .or. mod(iter,50)==1) then
       if (iter>1) write(logOF,*)'restart cg'
       if (do_precond) d=-pg
       if (.not. do_precond) d = -g
     else
       if (do_precond) then 
         !cg_beta = sum(pg*g)/sum(pg_old*g_old)
         !Polak–Ribière formula 
         cg_beta = max(0.0d0, sum(pg*(g-g_old))/sum(pg_old*g_old))
         d = -pg + cg_beta*d
       else 
         !cg_beta = sum(g*g)/sum(g_old*g_old)
         !Polak–Ribière formula 
         cg_beta = max(0.d0, sum(g*(g-g_old))/sum(g_old*g_old))
         d = -g + cg_beta*d
       endif 
     endif 

     g_old  = g
     pg_old = pg


     ! test dE/dtheta 
     !================
     dE_dtheta = dvol*sum(dE_dphi*d)
     if (dE_dtheta>0.d0) then 
       write(logOF,'(a)')' ---- dE/dtheta>0, switch to steepest decent ---'
       iter = 1
       goto 777
     endif 


     ! line search (compute the optimal theta)
     !=========================================
     der1 = dvol*sum(dE_dphi*d)
     theta1 = 1e-3
     phi_tmp = phi + theta1*d
     select case (vwreg_type) 
     case (1) 
       call compute_grad_vwreg1(nfft,n1,n2,n3,qvec,dvol,phi_tmp, &
                                q(isp),vext(:,isp),dtmp,dE_dphi_tmp,chempot)
     case (2) 
       call compute_etotal_vwreg2(nfft,n1,n2,n3,qvec,dvol,phi_tmp,q(isp),vext(:,isp), & 
                                  etotal_tmp,dE_dphi_tmp,chempot)
     case (3) 
       call compute_etotal_vwreg3(nfft,n1,n2,n3,qvec,dvol,phi_tmp,q(isp),vext(:,isp), & 
                                  etotal_tmp,dE_dphi_tmp,chempot)
     end select
     der2 = dvol*sum(dE_dphi_tmp*d)

     ! compute optimal theta (just to find at what theta der=0)
     theta = theta1/(der1-der2)*der1
     if (theta<0.d0) then 
       write(logOF,'(a)')'Warning theta<0, restart with steepest decent direction'
       iter = 1
       goto 777
     endif 

     ! update phi
     phi = phi + d*theta 

   enddo ! cg minimization 

!   ! DEBUG ==================================
!   call OF_fixpoint_solver(nfft,n1,n2,n3,ucvol,qvec,vext(:,isp),rho(:,isp),q(isp))
!   stop
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 enddo ! end of spin  





 write(logOF,'(a)')'done OF_cgmin_tot(). '
 call flush(logOF)



contains 



  subroutine vwreg3_energy(rho,etotal)
     real(8) :: rho(nfft), etotal, vwreg_phi,tf_energy,  & 
                lap(nfft), & 
                expo(nfft), rho_reg(nfft), tf_pot(nfft), vw_energy

     vwreg_phi=sqrt(vw_reg)
     call tf(nfft,rho,dvol,tf_pot,tf_energy)
     ! compute vw energy 
     expo    = sqrt(rho)/vwreg_phi
     rho_reg = (1.d0-exp(-expo))*rho
     call laplacian(nfft,n1,n2,n3,sqrt(rho_reg),qvec,lap)
     vw_energy = sum(-0.5d0*sqrt(rho_reg)*lap)*dvol
     etotal = vw_lam*vw_energy + tf_energy + dvol*sum(vext(:,1)*rho)
  end subroutine 


  !==================================================
  ! apply d\rho/d\phi to a vector 
  ! note that d\rho / dphi is a matrix, because
  ! rho is defined as 
  ! 
  !  rho = phi*phi/<phi,phi>*Ne
  !==================================================
  subroutine apply_drho_dphi(nfft,phi,dvol,q,x,Ax)

    implicit none 
    integer :: nfft
    real(8) :: phi(nfft),dvol, & 
               x(nfft),Ax(nfft),q,norm

    norm = sum(dvol*phi**2)

    Ax = 2.d0*phi*x*q/norm - phi**2*2.d0*q/norm**2*sum(dvol*phi*x)

  end subroutine apply_drho_dphi



  subroutine apply_drho_dphi_tr(nfft,phi,dvol,q,x,Ax)

    implicit none 
    integer :: nfft
    real(8) :: phi(nfft),dvol, & 
               x(nfft),Ax(nfft),q,norm

    norm = sum(dvol*phi**2)
    Ax = 2.d0*phi*x*q/norm - phi*2.d0*q/norm**2*sum(dvol*phi**2*x)

  end subroutine apply_drho_dphi_tr

end subroutine  OF_cgmin_tot






