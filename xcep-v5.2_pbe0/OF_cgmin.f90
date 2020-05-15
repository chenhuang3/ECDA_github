!
! Use L-BFGS or conjugate gradient method to direct minimize total energy w.r.t. electron density 
! The square root of electron density is used as the basic variable 
! the new density is 
! 
!    phi_new = phi_old * cos(theta) + lbfgs_dir * sin(theta)
!
! on input: 
!   1) rho is the initial guess 
!   2) vext is the external  potential 
!
! on output: 
!   1) rho is the solution to TF+vW OFDFT

subroutine OF_cgmin(n1,n2,n3,nfft,nspin,qvec,ucvol,vext,rho,chempot)
 
 use comm_data 

 implicit none 

 logical :: user_init

 integer :: n1,n2,n3, & 
            nfft, nspin, isp 

 real(8) :: vext(nfft,nspin), & 
            rho(nfft,nspin), & 
            new_sq_rho(nfft), & 
            qvec(3,(n1/2+1)*n2*n3), & 
            resid, & 
            ucvol, & 
            dvol, &  
            chempot, & 
            q(nspin), q_tmp, & 
            coeff,phi0(nfft), & 
            phi(nfft), & 
            phi_tmp(nfft), & 
            lam, dconv, & 
            theta, & 
            dE_dphi(nfft), & 
            g(nfft), & 
            g_old(nfft), & 
            d(nfft), d_proj(nfft), & 
            cg_beta, & 
            cg_alpha

 logical :: print_info = .false., & 
            do_precond 

 integer :: ii, opt_method, & 
            iter, mix_case, & 
            ipt

 real(8) :: etotal, ke, pi=3.1415926d0, & 
            etotal_tmp, & 
            etotal_tmp2, & 
            dE_dtheta, theta0,  & 
            pg(nfft),pg_old(nfft),lap(nfft), & 
            AA, BB, orth 

 ! >>>>>>>>>>>>>> function <<<<<<<<<<<<!

 write(logOF,'(a)')NEW_LINE('a')//'OF_cgmin() ...'

 dvol = ucvol / dble(nfft)

 do isp=1,nspin 
   q(isp) = dvol*sum(rho(:,isp))
   write(logOF,'(a,2es12.4,a,i3)')'vext: ',& 
     minval(vext(:,isp)),maxval(vext(:,isp)),' for spin:',isp
 enddo
 write(logOF,'(a,f12.5,a,f8.4)')'q: ',q(1:nspin),' vw_lambda:',vw_lam
 call flush(logOF)


 do isp=1,nspin 

   phi = sqrt(rho(:,isp))
   do_precond = .false.

!   ! randomize phi 
!   do ipt=1,nfft 
!     phi(ipt) = rand()
!   enddo 
!   phi = phi/sqrt(sum(phi**2*dvol)) * sqrt(q(isp))

   
   theta = 0.0d0

   ! cg method for minimize system energy 
   !======================================
   iter = 0
   do while(.true.)
     iter = iter + 1

     call compute_OF_E(nfft,n1,n2,n3,nspin,vw_lam,qvec,ucvol,phi,vext(:,isp),dvol,etotal) 
     call compute_OF_g(nfft,n1,n2,n3,nspin,vw_lam,qvec,ucvol,phi,vext(:,isp),dvol,dE_dphi) 

     chempot = dvol*sum(dE_dphi*phi)/2.d0/q(isp)  ! chemical potential 

     ! the true gradient d(Lagragian) / d(phi)
     g = dE_dphi - chempot*2.d0*phi


     ! switch on precondition
     !========================
     if (.not. do_precond) then 
       if (maxval(abs(g))<1e-1 ) then 
         ! compute kinetic energy 
         phi_tmp = phi/sqrt(sum(dvol*phi**2))
         call laplacian(nfft,n1,n2,n3,phi_tmp,qvec,lap)
         ke = -sum(dvol*phi_tmp*lap)/2.d0
         write(logOF,'(a)') '  **** switch on precondition ****'
         do_precond = .true.
         iter = 1 ! restart cg
       endif 
     endif 

     ! precondition g => pg
     !================================
     if (do_precond) then 
       call precond_teter(n1,n2,n3,nfft,qvec,ke,g,pg)
       ! make pg orthorgnal to phi
       pg = pg - phi*sum(pg*phi)*dvol/sum(dvol*phi**2)
     endif 

     !! restart CG every 50 steps
     !!===========================
     !if (iter==50) then 
     !  iter = 1
     !  do_precond = .false. 
     !  write(logOF,'(a)') '  *** restart cg *** '
     !endif 

777 continue 

     ! compute CG direction 
     !=====================
     if (iter==1) then
       if (do_precond) d=-pg
       if (.not. do_precond) d = -g
     else
       if (do_precond) then 
         cg_beta = sum(pg*g)/sum(pg_old*g_old)
         d = -pg + cg_beta*d
       else 
         cg_beta = sum(g*g)/sum(g_old*g_old)
         d = -g + cg_beta*d
       endif 
     endif 

     g_old  = g
     pg_old = pg


     ! project and normalize 
     d_proj = d - phi*sum(d*phi)/sum(phi*phi)
     d_proj = d_proj/sqrt(sum(d_proj**2*dvol))*sqrt(q(isp))


    ! if ( mod(iter,10)==1) then 
       write(logOF,'(a,i4,a,f12.7,a,f9.5,a,es10.2,a,es8.2,a,es10.2)') & 
       'it: ',iter,'  etot: ',etotal,'   Ef: ',chempot, & 
       '  theta: ',theta,'  |g|: ',sqrt(sum(g**2)*dvol),' gmax:',maxval(abs(g))
       call flush(logOF)
    ! endif 


     ! computing theta 
     !================
     ! first derivative dE/d(theta)
     !
     ! Eq. 5.29 (without the 0.5 coefficient) 
     !  in "Iterative minimization techniques for ab initio total-energy
     !      calculations: molecular dynamics and conjugate gradients"
     !
     !=============================================
     dE_dtheta = dvol*sum(dE_dphi*d_proj)
     if (dE_dtheta>0.d0) then 
       write(logOF,*)'dE_dtheta>0, swtich to steepst decent'
       iter = 1
       do_precond = .false. 
       goto 777
     endif 

     ! second derivative d^2 E / d(theta^2)
     !
     ! Eq. 5.28 in "Iterative minimization techniques for ab initio total-energy 
     !              calculations: molecular dynamics and conjugate gradients"
     !
     ! compute the curvature (central finite difference)
     ! 
     theta0  = 1e-2
     phi_tmp = cos(theta0)*phi + sin(theta0)*d_proj
     call compute_OF_E(nfft,n1,n2,n3,nspin,vw_lam,qvec,ucvol,phi_tmp,vext(:,isp),dvol,etotal_tmp)

     phi_tmp = cos(-theta0)*phi + sin(-theta0)*d_proj
     call compute_OF_E(nfft,n1,n2,n3,nspin,vw_lam,qvec,ucvol,phi_tmp,vext(:,isp),dvol,etotal_tmp2)

     ! second derivative: d^2 E/d(theta^2) at theta=0
     AA = (etotal_tmp + etotal_tmp2 - 2.d0*etotal) / theta0**2

     theta = 0.5d0*atan(-2.d0*dE_dtheta/AA)   ! Eq. 5.37 in Above paper
     if (theta<0) then 
        write(logOF,*)'theta: ',theta
        do ii=1,100
           theta = (ii-1)*0.002d0
           phi_tmp = cos(theta)*phi + sin(theta)*d_proj
           call compute_OF_E(nfft,n1,n2,n3,nspin,vw_lam,qvec,ucvol,phi_tmp,vext(:,isp),dvol,etotal_tmp)
           write(logOF,*)theta,etotal_tmp
        enddo
        stop
     endif 

     ! job fininshed? 
     if (sqrt(sum(g**2)*dvol)<1e-3 .AND. iter>=20) then 
       write(logOF,'(a,i4,a,f12.7,a,f9.5,a,es9.2,a,es8.2)') & 
       'it: ',iter,'  etot: ',etotal,' ef: ',chempot,' theta: ',theta, & 
       ' |g|: ',sqrt(sum(g**2)*dvol)
       call flush(logOF)
       write(logOF,'(a)')'OF-DFT minimized for this spin'
       exit 
     endif 

     phi =  phi*cos(theta) + sin(theta)*d_proj

   enddo ! minimization 
   
   rho(:,isp) = phi**2

 enddo ! end of spin  


 write(logOF,'(a)')'done OF_cgmin(). '
 call flush(logOF)


end subroutine  OF_cgmin
