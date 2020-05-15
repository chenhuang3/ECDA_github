!
! INPUT: 
!  rho    -- current density 
!  vperb  -- perturbing potential 
! 
! OUTPUT: 
!  delta_rho -- first order change of density 
!  delta_mu  -- first order change of chemical potential
!
subroutine OF_dfpt(nfft,cell_nfft,ucvol,qvec,rho,vperb,delta_rho,delta_mu)

 use comm_data 

 implicit none 

 ! external vars 
 integer :: nfft,cell_nfft(3)

 real(8),parameter :: pi=3.1415926d0, & 
                      cTF = 3.d0/10.d0*(3.d0*pi**2)**(2.d0/3.d0)

 real(8) :: delta_rho(nfft), &     
            delta_mu, &        
            chempot, rho(nfft), & 
            vperb(nfft), mach_eps, & 
            qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), & 
            ucvol 

 ! local vars
 logical :: do_precond

 integer :: k, n1, n2, n3, j
 real(8) :: vwreg_phi, pre_ke,  & 
            q, fstep,  etotal, phi(nfft), & 
            arr_tmp(nfft), rho_tmp(nfft), & 
            mu,ef, vext(nfft), drho_sum, shift = 1.d0, dnorm, & 
            cg_z(nfft), cg_zold(nfft), & 
            cg_r(nfft), & 
            cg_x(nfft), cg_Ax(nfft), & 
            cg_p(nfft),  & 
            cg_Ap(nfft), cg_b(nfft), cg_rold(nfft), & 
            cg_alpha, cg_beta, dtmp, norm, & 
            delta_vtf(nfft), resid, & 
            dvol, tf_energy, vw_energy, & 
            wf_reg(nfft), wf_regN(nfft), & 
            tf_kernel(nfft), & ! TF kernel 
            dE_dphi1(nfft), & 
            dE_dphi2(nfft), & 
            lap(nfft), &     
            vscf(nfft), & 
            vw_pot(nfft), & 
            vtf(nfft), & 
            nelectr, & 
            fake_chempot

 ! for complex step diff
 complex(8) :: rho_complex(nfft), & 
               dEdrho_complex(nfft), & 
               phi_complex(nfft),  & 
               dE_dphi_complex(nfft) 

 ! >>>>>>>>>> function begins <<<<<<<<<<< !

 write(logOF,'(a)') NEW_LINE('a')//'enter OF_dfpt()'

 n1=cell_nfft(1)
 n2=cell_nfft(2)
 n3=cell_nfft(3)

 dvol  = ucvol/dble(nfft)
 q = sum(rho*dvol)

 
 !====================================
 ! compute external potential for rho
 !====================================
 call tf(nfft,rho,dvol,vtf,tf_energy)

 ! compute VW potential -- dVW/drho
 select case(vwreg_type)
 case (1)
   call vw(nfft,n1,n2,n3,rho,qvec,dvol,vw_reg,vw_pot,vw_energy)
 case (2)
   arr_tmp = rho/vw_reg
   rho_tmp = (1.d0-exp(-arr_tmp))*rho
   call laplacian(nfft,n1,n2,n3,sqrt(rho_tmp),qvec,lap)
   ! compute d sqrt(rho')/d rho
   vw_pot = (1.d0-exp(-arr_tmp) & 
           + exp(-arr_tmp) * rho/vw_reg)/2.d0/sqrt(rho_tmp)
   vw_pot = -lap * vw_pot
 case (3) 
   vwreg_phi = sqrt(vw_reg)
   arr_tmp = rho**0.5/vwreg_phi
   rho_tmp = (1.d0-exp(-arr_tmp))*rho
   call laplacian(nfft,n1,n2,n3,sqrt(rho_tmp),qvec,lap)
   ! compute d sqrt(rho')/d rho
   vw_pot = (1.d0 - exp(-arr_tmp) + exp(-arr_tmp)*sqrt(rho)/2.d0/vwreg_phi)/2.d0/sqrt(rho_tmp)
   vw_pot = -lap*vw_pot
 endselect

 fake_chempot = 0.d0
 vext = fake_chempot - vw_lam*vw_pot - vtf


 ! compute vscf and Ef for |wf_reg>.
 ! wf_reg is the solution of the equation
 !
 ! -1/2*nabla^2 |wf_reg> + vscf|wf_reg> = mu|wf_reg>
 !
 vscf = (vext + vtf)/vw_lam   ! for wf_reg, vscf is scaled by vw_lam
 mu   = fake_chempot/vw_lam   ! for wf_reg, mu is scaled by vw_lam 

 tf_kernel = 5.d0/3.d0*cTF*2.d0/3.d0*rho**(-1.d0/3.d0) 

 wf_reg  = sqrt(rho + vw_reg)
 wf_regN = wf_reg/sqrt(sum(wf_reg**2*dvol))


 write(logOF,'(a,f12.4)')  'total_Q:  ',sum(rho)*dvol
 write(logOF,'(a,f12.4)')  'vw_lam:   ',vw_lam
 write(logOF,'(a,es12.4)') 'vw_reg:   ',vw_reg
 write(logOF,'(a,2es12.4)')'min/max(vks):   ',minval(vext),maxval(vext)
 write(logOF,'(a,2es12.4)')'min/max(vperb): ',minval(vperb),maxval(vperb)


 !======================================
 ! vwreg_type = 1
 !
 if (vwreg_type==1) then 
   !======================================================================
   ! Conjugate gradient method for solving modified Sternheimer equation 
   !
   ! [ H_scf - mu + shift*Pv + 
   !   2*Pc*wf_reg*tf_kernel*wf_reg ] |delta_wf_reg> = - Pc*wf_reg*vperb
   !
   ! with  H_scf = -1/2*nabla + vscf 
   !       vscf  = (vtf + vext)/vw_lam
   !
   ! NOTE: This is a modified Sternheimer equation which consider
   ! the change of vscf for a variation of wf_reg. This is only 
   ! possible for OF-DFT, and is not possible KS-DFT.
   !======================================================================
   write(logOF,'(a,i3)')'vwreg_type:',vwreg_type

   do_precond = .false. 
   cg_x = 0.d0


!=============================
! restart point for the cg
999 continue 

   call compute_Ax_vwreg1(vscf,mu,wf_regN,wf_reg,tf_kernel,shift,cg_x,cg_Ax)

   cg_b = wf_reg*vperb/vw_lam 
   cg_b = -(cg_b-wf_reg*sum(cg_b*wf_reg)/sum(wf_reg**2))

   cg_r = cg_b - cg_Ax

   ! precondition for the laplacian operator
   ! damp out q^2 terms 
   if (do_precond) then 
     call precond_teter(n1,n2,n3,nfft,qvec,pre_ke,cg_r,cg_z)
     cg_p = cg_z
   else 
     cg_p = cg_r
   endif 
  
  
   !==============================
   ! conjugate gradient loop
   !
   k = 0
   do while (.true.) 
     k = k + 1
  
     call compute_Ax_vwreg1(vscf,mu,wf_regN,wf_reg,tf_kernel,shift,cg_p,cg_Ap)

     if (.not. do_precond) then 
       cg_alpha = sum(cg_r**2)/sum(cg_p*cg_Ap)
     else 
       cg_alpha = sum(cg_r*cg_z)/sum(cg_p*cg_Ap)
     endif 

     cg_x     = cg_x + cg_alpha*cg_p
     cg_rold  = cg_r
     cg_zold  = cg_z
     cg_r     = cg_r - cg_alpha*cg_Ap


     ! precondition for the laplacian operator
     ! damp out q^2 terms 
     if (do_precond) then 
       call precond_teter(n1,n2,n3,nfft,qvec,pre_ke,cg_r,cg_z)
       cg_beta = sum(cg_z*cg_r)/sum(cg_zold*cg_rold)
       cg_p    = cg_z + cg_beta*cg_p
     else 
       cg_beta = sum(cg_r**2)/sum(cg_rold**2)
       cg_p    = cg_r + cg_beta*cg_p
     endif 
  
  
     ! information 
     resid = dvol*sum(cg_r**2)
     delta_rho = 2.d0*wf_reg*cg_x
     delta_vtf = tf_kernel*delta_rho
     delta_mu  = sum(wf_reg*(delta_vtf+vperb)*wf_reg)*dvol/sum(wf_reg**2*dvol)
  
  
     if (mod(k,10)==1) then 
       write(logOF,'(a,i4,a,es8.2,a,f8.5,a,2es12.3)') & 
         '(OF_dfpt) it:',k,'  |r|: ',dvol*sum(cg_r**2),'   delta_ef: ',delta_mu, & 
         '  delta_rho:',minval(delta_rho),maxval(delta_rho) 
       call flush(logOF)
     endif
     if (k>1000) then 
       print *,'error! cg in OF_dfpt() does not converge after 10000 iterations, stop'
       call flush(6)
       stop
     endif 
     if (k>=20 .and. resid<1e-6) then 
       write(logOF,'(a,i4,a,es8.2,a,f8.5,a,2es12.3)') & 
       '(OF_dfpt) it:',k,'  |r|: ',dvol*sum(cg_r**2),'   delta_ef: ',delta_mu, & 
       '  delta_rho:',minval(delta_rho),maxval(delta_rho) 
       call flush(logOF)
       exit 
     endif


     ! switch on preconditioner?
     !==========================
     if ((.not. do_precond) .and. resid<0.1d0) then 
       ! Note: Teter's preconditioner does not work.
       !       The convergence becomes really worse, I do not know why
       !       So I switched preconditioner off.
       do_precond = .true.
       call laplacian(nfft,n1,n2,n3,cg_x,qvec,lap)
       pre_ke = -sum(dvol*cg_x*lap)/2.d0
       if (do_precond) then 
         write(logOF,'(a,es12.4,a)')'*** precondition is ON, with pre_ke: ',pre_ke,' ***'
       endif 
       goto 999
     endif 

   enddo
 endif 






 !======================================
 ! vwreg_type = 2 or 3
 !
 if ( (vwreg_type==2) .or. (vwreg_type == 3))  then 
   write(logOF,'(a,i3)')'vwreg_type:',vwreg_type
   write(logOF,'(a,f14.6)')'q: ',q

   do_precond = .true.
   if (do_precond) write(logOF,'(a)')'*** precond is ON ***'

   pre_ke = 1.0d0

   ! This is the simplest way to set phi
   ! however, for this case, phi is general since rho is automatically 
   ! normalized as rho = phi**2/<phi,phi>*Ne
   phi = sqrt(rho)

   ! make cg_b 
   !===================
   norm = sqrt(dvol*sum(phi**2))
   cg_b = -vperb*2.d0*phi*q/norm**2 + 2.d0*q*phi/norm**4*sum(dvol*vperb*phi**2)
   cg_r = cg_b 

   if (do_precond) then 
     ! damp out q^2 term in laplacian operator to improve 
     ! condition number of the hessian matrix 
     call precond_teter(n1,n2,n3,nfft,qvec,pre_ke,cg_r,cg_z)
     cg_p = cg_z
   else 
     cg_p = cg_r
   endif 

   cg_x = 0.d0 
   k = 0

   !=============================
   ! conjugate gradient loop
   !=============================
   do while (.true.) 
     k = k + 1
  
     ! compute A x (finite difference) 
     !fstep = 1e-4
     !!mach_eps = 1e-16
     !!fstep = 2.d0*sqrt(mach_eps)*(1.d0+sqrt(dvol*sum(phi**2)))/sqrt(sum(dvol*cg_p**2))
     !call compute_etotal_vwreg3(nfft,n1,n2,n3,qvec,dvol,phi+fstep*cg_p,q,vext,etotal,dE_dphi1,ef)
     !call compute_etotal_vwreg3(nfft,n1,n2,n3,qvec,dvol,phi-fstep*cg_p,q,vext,etotal,dE_dphi2,ef)
     !cg_Ap = (dE_dphi1-dE_dphi2)/2.d0/fstep

     ! complex step differentiation 
     fstep = 1e-12
     phi_complex = cmplx(phi,fstep*cg_p)
     select case (vwreg_type)
     case(1)
       fstep = 1e-5
       !mach_eps = 1e-16
       !fstep = 2.d0*sqrt(mach_eps)*(1.d0+sqrt(dvol*sum(phi**2)))/sqrt(sum(dvol*cg_p**2))
       call compute_grad_vwreg1(nfft,n1,n2,n3,qvec,dvol,phi+fstep*cg_p,q,vext,dtmp,dE_dphi1,ef)
       call compute_grad_vwreg1(nfft,n1,n2,n3,qvec,dvol,phi-fstep*cg_p,q,vext,dtmp,dE_dphi2,ef)
       cg_Ap = (dE_dphi1-dE_dphi2)/2.d0/fstep
       !call compute_grad_complex_vwreg1(nfft,n1,n2,n3,qvec,dvol,phi_complex,q,vext,dE_dphi_complex)
       !cg_Ap = aimag(dE_dphi_complex)/fstep
     case(2)
       call compute_grad_complex_vwreg2(nfft,n1,n2,n3,qvec,dvol,phi_complex,q,vext,dE_dphi_complex)
       cg_Ap = aimag(dE_dphi_complex)/fstep
     case(3)
       call compute_grad_complex_vwreg3(nfft,n1,n2,n3,qvec,dvol,phi_complex,q,vext,dE_dphi_complex)
       cg_Ap = aimag(dE_dphi_complex)/fstep
     endselect

   
     ! cg parameters 
     if (do_precond) then 
       cg_alpha = sum(cg_r*cg_z)/sum(cg_p*cg_Ap)
     else 
       cg_alpha = sum(cg_r**2)/sum(cg_p*cg_Ap)
     endif 

     cg_x     = cg_x + cg_alpha*cg_p

     cg_rold  = cg_r
     cg_zold  = cg_z
     cg_r     = cg_r - cg_alpha*cg_Ap

     ! update beta and cg_p
     if (do_precond) then 
       ! damp out q^2 term in laplacian operator to improve 
       ! condition number of the hessian matrix 
       call precond_teter(n1,n2,n3,nfft,qvec,pre_ke,cg_r,cg_z)
       cg_beta  = sum(cg_z*cg_r)/sum(cg_zold*cg_rold)
       cg_p     = cg_z + cg_beta*cg_p
     else 
       cg_beta  = sum(cg_r**2)/sum(cg_rold**2)
       cg_p     = cg_r + cg_beta*cg_p
     endif 
  
  
     ! information 
     resid = dvol*sum(cg_r**2)
     norm = sqrt(sum(dvol*phi**2))


     ! Delta_rho = \int delta(rho) / delta(phi) * Delta_phi
     delta_rho = 2.d0*phi*cg_x/norm**2*q & 
               - phi**2/norm**4*2.d0*q*sum(dvol*phi*cg_x)

  
     write(logOF,'(a,i4,a,es8.2,a,2es12.3)') & 
       '(OF_dfpt) it:',k,'  |r|: ',dvol*sum(cg_r**2),& 
       '  delta_rho:',minval(delta_rho),maxval(delta_rho) 
     call flush(logOF)
     if (k>10000) then 
       print *,'error! cg in OF_dfpt() does not converge after 10000 iterations, stop'
       call flush(6)
       stop
     endif 

     ! check convergence 
     if ((k>=20) .and. (resid<1e-8)) then 
       write(logOF,'(a,i4,a,es8.2,a,f8.5,a,2es12.3)') & 
       '(OF_dfpt) it:',k,'  |r|: ',dvol*sum(cg_r**2),'   delta_ef: ',delta_mu, & 
       '  delta_rho:',minval(delta_rho),maxval(delta_rho) 
       write(logOF,'(a,2es14.6)')'delta_phi: ',minval(cg_x),maxval(cg_x)
       call flush(logOF)

       ! compute first order change of chemical potential 
       !
       ! mu = <dE/drho,rho>/Ne 
       !
       ! \Delta mu =  <\partial mu/\partial v, \Delta v> 
       !             +<\partial mu/\partial rho, \Delta rho>
       !
       ! these two parts are computed below 
       !
       ! part 1 due to change of vext 
       delta_mu  = sum(rho*vperb*dvol)/q
       
       ! part 2 due to change of density 
       fstep = 1e-4
       select case(vwreg_type)
       case (1)
         call compute_OF_dEdrho_vwreg1(nfft,n1,n2,n3,qvec,ucvol,rho+fstep*rho,vext,dE_dphi1)
         call compute_OF_dEdrho_vwreg1(nfft,n1,n2,n3,qvec,ucvol,rho-fstep*rho,vext,dE_dphi2)
       case (2)
         call compute_OF_dEdrho_vwreg2(nfft,n1,n2,n3,qvec,ucvol,rho+fstep*rho,vext,dE_dphi1)
         call compute_OF_dEdrho_vwreg2(nfft,n1,n2,n3,qvec,ucvol,rho-fstep*rho,vext,dE_dphi2)
       case (3)
         call compute_OF_dEdrho_vwreg3(nfft,n1,n2,n3,qvec,ucvol,rho+fstep*rho,vext,dE_dphi1)
         call compute_OF_dEdrho_vwreg3(nfft,n1,n2,n3,qvec,ucvol,rho-fstep*rho,vext,dE_dphi2)
       endselect
       arr_tmp = (dE_dphi1-dE_dphi2)/2.d0/fstep
       delta_mu = delta_mu + sum(dvol*arr_tmp*delta_rho)/q

       exit 

     endif
   enddo  ! end of CG
 endif  ! vwreg_type = 3



 !
 ! check the sum of delta_rho
 !
 drho_sum = abs(sum(delta_rho)*dvol)
 write(logOF,'(a,es12.4)')'[OF_dfpt] sum(delta_rho)*dvol: ',sum(delta_rho)*dvol
 if ( drho_sum>1e-2 ) then 
   write(logOF,*)'OF_dfpt(): sum of delta_rho is not zero! stop'
   stop
 endif 
 
 

 !
 ! computer first order change of total energy for the vperb 
 !
 call tf(nfft,rho,dvol,vtf,dtmp)
 call vw(nfft,n1,n2,n3,rho,qvec,dvol,vw_reg,vw_pot,dtmp)
 write(logOF,'(a,f12.4)')'first order change in etotal <rho,vperb>: ',dvol*sum(rho*vperb)
 write(logOF,'(a,f12.8)')'first order change in Ef (dE/drho):       ',delta_mu
 write(logOF,'(a)')'OF_dfpt finished.'
 call flush(logOF)



contains 

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


  
  subroutine compute_Ax_vwreg1(vscf,mu,wf_regN,wf_reg,tf_kernel,shift,x,Ax)
     implicit none 
     real(8) :: vscf(nfft), mu, &
                shift, & 
                tf_kernel(nfft),  & 
                wf_reg(nfft),  &   ! unnormalized 
                wf_regN(nfft), &   ! normalized to one 
                x(nfft), Ax(nfft),  & 
                lap(nfft), arr_tmp(nfft)

     call laplacian(nfft,n1,n2,n3,x,qvec,lap)

     Ax = - lap/2.d0 + vscf*x - mu*x &       ! H_scf*|cg_p>
          + shift*wf_regN*dvol*sum(wf_regN*x)   ! shift*Pv*|cg_p> 
  
     arr_tmp = x - wf_regN*sum(dvol*x*wf_regN)       ! project out wf_reg
     arr_tmp = wf_reg*tf_kernel*wf_reg*arr_tmp/vw_lam
     arr_tmp = arr_tmp - wf_regN*sum(dvol*arr_tmp*wf_regN) ! project out wf_reg

     Ax = Ax + 2.d0*arr_tmp
  end subroutine compute_Ax_vwreg1


end subroutine OF_dfpt

