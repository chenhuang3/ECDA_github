!
!  Change of electron denisty for a change electron number with 
!  the v_ext fixed. 
!
!  different from the OF_dfpt.f90 code, 
!  here we consider the wf_reg which is normalized to one. 
!  only implemented for vwreg_type = 1
!
!  Created by Chen Huang 12/28/2019
!
subroutine OF_drho_dN_vfix(nfft,cell_nfft,ucvol,qvec,rho,drho_dN)

 use comm_data 

 implicit none 

 ! external vars 
 integer :: nfft,cell_nfft(3)

 real(8),parameter :: pi=3.1415926d0, & 
                      cTF = 3.d0/10.d0*(3.d0*pi**2)**(2.d0/3.d0)

 real(8) :: drho_dN(nfft), &    
            chempot, rho(nfft), & 
            vperb(nfft), mach_eps, & 
            qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), & 
            ucvol 

 ! local vars
 logical :: do_precond

 integer :: k, n1, n2, n3, j
 real(8) :: vwreg_phi, pre_ke,  & 
            q, fstep, ef, etotal, phi(nfft), & 
            arr_tmp(nfft), rho_tmp(nfft), & 
            mu, vext(nfft), drho_sum, shift = 1.d0, dnorm, & 
            cg_z(nfft), cg_zold(nfft), & 
            cg_r(nfft), cg_x(nfft),  &
            cg_Ax(nfft), cg_p(nfft),  & 
            cg_Ap(nfft), cg_b(nfft), cg_rold(nfft), & 
            cg_alpha, cg_beta, dtmp, norm, & 
            resid, & 
            dvol, tf_energy, vw_energy, & 
            wfn(nfft), wf_reg(nfft), wf_regN(nfft), & 
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

 write(logOF,'(a)') NEW_LINE('a')//'enter OF_drho_dN_vfix()'

 if (vwreg_type==3 .or. vwreg_type==2) then 
   write(logOF,*)'OF_drho_dN_vfix is not coded for vwreg_type=2 and 3 stop'
   stop
 endif 

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

 ! normalized wfn 
 wfn = sqrt(rho)
 wfn = wfn/sqrt(sum(wfn**2*dvol))

 wf_reg  = sqrt(rho+vw_reg)
 wf_regN = wf_reg/sqrt(sum(wf_reg**2*dvol))

 write(logOF,'(a,f12.4)')  'total_Q:  ',sum(rho)*dvol
 write(logOF,'(a,f12.4)')  'vw_lam:   ',vw_lam
 write(logOF,'(a,es12.4)') 'vw_reg:   ',vw_reg


 !======================================
 ! vwreg_type = 1
 !
 if (vwreg_type==1) then 

   write(logOF,'(a,i3)')'vwreg_type:',vwreg_type
   cg_x = 0.d0 
   do_precond = .false.   ! precondition is turned on automatically in CG


! restart point for CG 
999 continue 

   call compute_Ax_vwreg1(shift,vscf,mu,q,tf_kernel,wf_regN,cg_x,cg_Ax)

   ! Note: different from OF_dfpt code, 
   ! the formalism is based on normalized wf_regN here
   cg_b = wf_regN*tf_kernel*wf_regN**2/vw_lam
   cg_b = -(cg_b-wf_regN*sum(dvol*cg_b*wf_regN))  !project out wf_reg

   cg_r = cg_b - cg_Ax

   if (do_precond) then 
     call precond_teter(n1,n2,n3,nfft,qvec,pre_ke,cg_r,cg_z)
     cg_p = cg_z
   else 
     cg_p = cg_r
   endif 
  
  
   !=======================================================================
   ! conjugate gradient loop
   !
   ! Note that, different from OF_dfpt code,
   ! the code is based on wf_regN which is normalized to 1,
   !
   !-----------------------------------------------------------------------
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
   !
   !-----------------------------------------------------------------------
   k = 0
   do while (.true.) 
     k = k + 1
     
     call compute_Ax_vwreg1(shift,vscf,mu,q,tf_kernel,wf_regN,cg_p,cg_Ap)

     ! cg parameters 
     if (.not. do_precond) then 
       cg_alpha = sum(cg_r**2)/sum(cg_p*cg_Ap)
     else 
       cg_alpha = sum(cg_r*cg_z)/sum(cg_p*cg_Ap)
     endif 

     cg_x     = cg_x + cg_alpha*cg_p
     cg_rold  = cg_r
     cg_zold  = cg_z
     cg_r     = cg_r - cg_alpha*cg_Ap

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
     if (mod(k,10)==1) then 
       write(logOF,'(a,i4,a,es8.2,a,2es12.4)') & 
         'cg_iter:',k,'  |r|: ',dvol*sum(cg_r**2),' cg_x:',minval(cg_x),maxval(cg_x)
       call flush(logOF)
     endif
     if (k>1000) then 
       print *,'error! cg in OF_dfpt() does not converge after 10000 iterations, stop'
       call flush(6)
       stop
     endif 
     if (resid<1e-6) then 
       write(logOF,*)'converged. resid < 1e-6'
       write(logOF,'(a,i4,a,es8.2,a,2es12.4)') & 
         'cg_iter:',k,'  |r|: ',dvol*sum(cg_r**2),' cg_x:',minval(cg_x),maxval(cg_x)
       call flush(logOF)
       exit 
     endif


     ! switch on preconditioner? 
     !==========================
     !!if ( (.not. do_precond) .and. resid<0.1d0 ) then 
     if ( .not. do_precond ) then 
       ! Note: Teter's preconditioner does not work.
       !       The convergence becomes really worse, I do not know why
       !       So I switched preconditioner off.
       do_precond = .true.
     !  ! compute the kinetic energy of cg_x
     !  call laplacian(nfft,n1,n2,n3,cg_x,qvec,lap) 
     !  pre_ke = -sum(dvol*cg_x*lap)/2.d0
       pre_ke = 0.1d0
       if (do_precond) then 
         write(logOF,'(a,es12.4,a)')'*** precondition is ON, pre_ke: ',pre_ke,' ***'
       endif 
       goto 999
     endif 

   enddo ! end of cg 

    
   ! compute drho/dN, with rho defined as
   ! rho = wf_regN**2 * q - c
   !
   drho_dN = wf_regN**2 + 2.d0*wf_regN*q*cg_x

   write(logOF,'(a,2es14.6)')'drho_dN: ',minval(drho_dN),maxval(drho_dN)
   write(logOF,'(a,f12.4)')'sum[drho_dN]',sum(drho_dN)*dvol

 endif 


 write(logOF,'(a)')'OF_drho_dN_vfix() finished.'
 call flush(logOF)



contains  



   subroutine compute_Ax_vwreg1(shift,vscf,mu,q,tf_kernel,wf_regN,x,Ax)
     implicit none 
     real(8) :: shift, vscf(nfft), mu, q, & 
                tf_kernel(nfft), wf_regN(nfft), & 
                x(nfft), Ax(nfft)

     ! compute H_scf |cg_p>
     call laplacian(nfft,n1,n2,n3,x,qvec,lap)

     Ax = - lap/2.d0 + vscf*x - mu*x &          ! (H_scf-eigen)*|cg_p>
          + shift*wf_regN*dvol*sum(wf_regN*x)   ! shift*Pv*|cg_p> 
  

     ! Note: different OF-dfpt code, here the code 
     ! is based on normalized wf_regN. 
     !
     arr_tmp = x - wf_regN*sum(dvol*wf_regN*x)       ! project out wf_reg
     arr_tmp = wf_regN*tf_kernel*wf_regN*arr_tmp/vw_lam*q
     arr_tmp = arr_tmp - wf_regN*sum(dvol*wf_regN*arr_tmp) ! project out wf_reg
     Ax = Ax + 2.d0*arr_tmp

   endsubroutine 

end subroutine OF_drho_dN_vfix

