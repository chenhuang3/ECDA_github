!
!  inverse of the OF perturbation theory
!  for given delta_rho, solve for vperb
!
subroutine OF_dfpt_inv(nfft,cell_nfft,ucvol,qvec,rho,delta_rho,vperb)

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
 integer :: k, n1, n2, n3, j
 real(8) :: vwreg_phi, pre_ke,  & 
            q, fstep, ef, etotal, phi(nfft), & 
            arr_tmp(nfft), rho_tmp(nfft), & 
            mu, vext(nfft), drho_sum, shift = 1.d0, dnorm, & 
            norm, resid, & 
            dvol, tf_energy, vw_energy, & 
            fvec(nfft), dmu2, dmu1, & 
            vperb_new(nfft), & 
            dE_dphi1(nfft), & 
            dE_dphi2(nfft), & 
            lap(nfft), &     
            vw_pot(nfft), & 
            vtf(nfft), & 
            nelectr, & 
            fake_chempot

 ! for complex step diff
 complex(8) :: rho_complex(nfft),  & 
               dEdrho_complex(nfft) 

 ! >>>>>>>>>> function begins <<<<<<<<<<< !

 write(logOF,'(a)') NEW_LINE('a')//'enter OF_dfpt_inv()'

 n1=cell_nfft(1)
 n2=cell_nfft(2)
 n3=cell_nfft(3)
 dvol  = ucvol/dble(nfft)
 q = sum(rho*dvol)


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

 write(logOF,'(a,f12.4)')  'total_Q:  ',sum(rho)*dvol
 write(logOF,'(a,f12.4)')  'vw_lam:   ',vw_lam
 write(logOF,'(a,es12.4)') 'vw_reg:   ',vw_reg
 write(logOF,'(a,2es12.4)')'min/max(vext):  ',minval(vext),maxval(vext)
 write(logOF,'(a,2es12.4)')'min/max(d_rho): ',minval(delta_rho),maxval(delta_rho)


 ! vwreg_type = 2
 if (vwreg_type==2) then 
   write(logOF,*)'OF_dfpt_inv() for vwreg_type=2 has not been coded yet. stop!'
   call flush(logOF)
   stop
 endif 



 !======================================
 ! vwreg_type = 1 or 3
 !
 if ((vwreg_type == 1) .or. (vwreg_type==3)) then 

   select case(vwreg_type)
   case (1)
     write(logOF,'(a)')'vwreg_type=1'
     !! complex step differentiation. 
     !!
     !! NOTE: We have to use complex step differentiation instead of 
     !!       traditional finite differentiation, because this avoids 
     !!       the problem that density cannot be negative, 
     !fstep = 1e-10
     !rho_complex = cmplx(rho,fstep*delta_rho)
     !call compute_OF_dEdrho_complex_vwreg1(nfft,n1,n2,n3,ucvol,qvec,rho_complex,dEdrho_complex)
     !fvec = aimag(dEdrho_complex)/fstep
     !write(logOF,*)'fvec: ',minval(fvec),maxval(fvec)
     
     ! Analytical way to compute fvec
     call apply_hess_vwreg1(nfft,n1,n2,n3,qvec,rho,delta_rho,fvec)
     !!write(logOF,*)'fvec: ',minval(fvec),maxval(fvec)
   case (3)
     write(logOF,'(a)')'vwreg_type=3'
     ! complex step differentiation. 
     !
     ! NOTE: We have to use complex step differentiation instead of 
     !       traditional finite differentiation, because this avoids 
     !       the problem that density cannot be negative, 
     fstep = 1e-10
     rho_complex = cmplx(rho,fstep*delta_rho)
     call compute_OF_dEdrho_complex_vwreg3(nfft,n1,n2,n3,qvec,ucvol,rho_complex,vext,dEdrho_complex)
     fvec = aimag(dEdrho_complex)/fstep
     !fstep = 1e-4
     !call compute_OF_dEdrho_vwreg3(nfft,n1,n2,n3,qvec,ucvol,rho+fstep*delta_rho,vext,dE_dphi1)
     !call compute_OF_dEdrho_vwreg3(nfft,n1,n2,n3,qvec,ucvol,rho-fstep*delta_rho,vext,dE_dphi2)
     !fvec = (dE_dphi1-dE_dphi2)/2.d0/fstep 
   end select

   ! Note: Any shift of v_ext is a solution, which is also confirmed 
   ! in my derivation (see ECDA/OFDFT paper for details). 
   ! Therefore, here we just need to set 
   !
   !   Delta v_{ext} = d^2 E/(d rho(r) d rho(r')) Delta rho
   !
   vperb = -fvec 

 endif  ! vwreg_type = 3

 write(logOF,'(a,2es14.6)')'vperb: ',minval(vperb),maxval(vperb)
 write(logOF,'(a)')'leave OF_dfpt_inv().'
 call flush(logOF)


end subroutine OF_dfpt_inv

