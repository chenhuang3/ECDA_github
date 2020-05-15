!
! convert recip_input from recip space to real space
! or from real space to recip space
!
! Often, this code just transform real space array to 
! reciprocal space, and back and forth. 
!
! If you set check_norm = .True., the code will check the norm of wf to make sure the norm  is 1 
! the real space vector's imag part is always set to zero after the transformation 
!
! fft_dir = 1: recip to real
!              input: recip_input
!              output: res_real
!
! fft_dir = -1: real to recip
!              input: res_real
!              output: recip_real
!
subroutine recip_to_real_NEW(fft_dir,mpi_enreg,paral_kgb,nkpt,mpw,npw,ikpt,istwf_k, & 
                         nspinor,mgfft,nfft,ucvol,ngfft,&
                         kg_k,npw_k,recip_input,res_real,check_norm)  
 use defs_datatypes
 use defs_abitypes  
 use interfaces_53_ffts  ! In FORTRAN, if optinal vars are used, we must use interface!!
                         ! otherwise the optional vars will be screwed up
                         ! here, we must include the interface for fourwf() subroutine
                         ! this cost me four hours!

 use m_cgtools, only: dotprod_g

 implicit none 
 
 integer,intent(in)  :: fft_dir,npw,nspinor,nfft,nkpt,ngfft(18), & 
                        istwf_k,mgfft,ikpt,kg_k(3,mpw), & 
                        npw_k,paral_kgb,mpw

 logical :: check_norm
                 
 type(MPI_type):: mpi_enreg                 
                 
 real(kind=8)  :: recip_input(2,npw*nspinor),ucvol
 real(kind=8)  :: res_real(nfft,1)
 real(kind=8)  :: fftg_out(2,npw)
   
 ! local vars 
 integer       :: gbound(2*mgfft+8,2), nsppol, q1, q2, q3
 real(kind=8)  :: tmp_real(ngfft(4),ngfft(5),ngfft(6)), & 
                  wfraug(2,ngfft(4),ngfft(5),ngfft(6)), & 
                  weight,dummy(2,1), thetaF, dtmp,  & 
                  theta , sum1,sum2, re_norm, im_norm


!! debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real(8)          :: gvnl1(2,npw_k*nspinor)                  
 double precision :: dotr, doti
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 ! get gbound
 call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)       
 weight = 1.d0/ucvol
 
 if (fft_dir == 1) then 

   ! ----------------------
   ! recip -> real space
   ! ----------------------

 !  call dotprod_g(dotr,doti,istwf_k,npw_k*nspinor, & 
 !      2,recip_input,recip_input,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
 !   print *,'dotr, doti: ',dotr, doti

   dtmp = sum(recip_input**2)*ucvol/dble(nfft)
   if ( dtmp < 1e-14 ) then 
      print *,'(recip_to_real_NEW): norm of recip_input is < 1e-14, just retrun res_real=zero'
      res_real = 0.d0
      return 
   endif 

   wfraug = 0.d0 ! This is a must, since not all elements in wfraug will be filled 
   call fourwf(1,tmp_real,recip_input,dummy,wfraug,gbound,gbound,istwf_k,&
               kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,1, & 
               ngfft(4),ngfft(5),ngfft(6),0,paral_kgb,0,weight,weight)

   wfraug = wfraug / sqrt(ucvol)

   !
   ! remove the phase of wave function => enfore wave function to be real
   !
   ! the trick is that the wfr might be complex, we need to remove the phase
   !    wfr = phi_1 + phi_2*I 
   ! apply an phase to it 
   !    wfr=(phi_1+phi_2*I)*(cos(theta)+I*sin(theta))
   ! the imag part is 
   !    IM= phi_1*sin(theta)+phi_2*cos(theta)
   !
   ! minimize/maximize the norm of IM with respect to theta
   ! compute d IM/ d theta = 0, we have 
   !    tan(*2theta) = -2*sum(phi_1j*phi2j)/ sum(phi_1j^2 - phi_2j^2) = 0
   !
   sum1 = sum(wfraug(1,:,:,:)*wfraug(2,:,:,:))
   sum2 = sum(wfraug(1,:,:,:)**2 - wfraug(2,:,:,:)**2)
   theta = atan(-2.0*sum1/sum2)/2.d0

   !
   ! remove the phase: new_wf = wfraug*(cos(-theta)+I*sin(-theta))
   !
   tmp_real        = wfraug(1,:,:,:)*cos(theta) - wfraug(2,:,:,:)*sin(theta)
   wfraug(2,:,:,:) = wfraug(1,:,:,:)*sin(theta) + wfraug(2,:,:,:)*cos(theta)
   wfraug(1,:,:,:) = tmp_real

   !
   ! Note, this might MAXIMIZW the imaginary part as well.
   ! therefore, we check which part (real or imag) is larger later.
   ! if the imaginary part is large, we flip the real and imag parts
   !
   ! if the real part is much smaller than the imaginary part
   ! then the real and imaginary part should be flipped
   !
   re_norm = sum(wfraug(1,:,:,:)**2)*ucvol/dble(nfft)
   im_norm = sum(wfraug(2,:,:,:)**2)*ucvol/dble(nfft)
   
   if ( re_norm < im_norm ) then 
      tmp_real        = wfraug(1,:,:,:)
      wfraug(1,:,:,:) = wfraug(2,:,:,:) 
      wfraug(2,:,:,:) = tmp_real
   endif 
   
   ! check the norm of the imginary part 
   if ( sum(wfraug(2,:,:,:)**2)*ucvol/dble(nfft) > 1e-4 ) then 
     print *,'maxval(abs(wfraug(2,:,:,:))): ',maxval(abs(wfraug(2,:,:,:)))
     print *,'maxval(abs(wfraug(1,:,:,:))): ',maxval(abs(wfraug(1,:,:,:)))
     print *,'sum(wfraug(2,:,:,:)**2)*ucvol/dble(nfft): ',sum(wfraug(2,:,:,:)**2)*ucvol/dble(nfft)
     print *,'the imagniary part is nonzero stop,'
     print *,'Solution: I guess that you set istwfk = 1, so you can increase nbdbuf.'
     stop 
   endif 

   !!print *,maxval(wfraug(1,:,:,:)),maxval(wfraug(2,:,:,:))

   !  print *,'re+im: ',sum(wfraug**2)*ucvol/dble(nfft)
   !  print *,'=>>', maxval(abs(wfraug(2,:,:,:)))

   ! print *,'after mod, max(abs(real)): ',maxval(abs(wfraug(1,:,:,:)))
   ! print *,'after mod, max(abs(imag)): ',maxval(abs(wfraug(2,:,:,:)))
   ! print *,wfraug(1,1,1,1),wfraug(2,1,1,1)
   ! print *,wfraug(1,2,1,1),wfraug(2,2,1,1)
   ! print *,wfraug(1,3,1,1),wfraug(2,3,1,1)
   ! print *,wfraug(1,34,23,14),wfraug(2,34,23,14)
   ! print *,wfraug(1,34,2,14), wfraug(2,34,2,14)
   ! print *,wfraug(1,34,45,14),wfraug(2,34,45,14)
   ! print *,'theta: ',theta

   ! transfer wfraug(1,n4,n5,n6) array ==> tmp_phir(nfft,1) array
   nsppol = 1
   res_real = 0.d0
   call fftpac(1,mpi_enreg,nsppol, & 
               ngfft(1),ngfft(2),ngfft(3), & 
               ngfft(4),ngfft(5),ngfft(6), & 
               ngfft,res_real,wfraug(1,:,:,:),1)

   ! check norm of wf
   if (check_norm) then 
     if (  abs(1.d0-sum(res_real*res_real*ucvol/dble(nfft)))>1e-3 ) then 
       print *,'in 94_scfcv/recip_to_real.F90: norm of wf is not one'
       print *,'norm is :',sum(res_real*res_real*ucvol/dble(nfft)),' theta: ',theta
       stop 
     endif
   endif

   !write(*,*)'[recip_to_real] min/max(phir)=',minval(res_real),maxval(res_real), &
   !' norm(phir^2)=',sum(res_real**2) * ucvol/dble(nfft)

 else
   ! ----------------------------
   ! real space -> recip space
   ! ----------------------------

!   write(6,'(a)')'(recip_to_real) real -> recip transform '
   nsppol = 1
   wfraug = 0.d0
   call fftpac(1,mpi_enreg,nsppol,& 
               ngfft(1),ngfft(2),ngfft(3), & 
               ngfft(4),ngfft(5),ngfft(6), & 
               ngfft,res_real,wfraug(1,:,:,:),2)

   wfraug = wfraug * sqrt(ucvol)
   fftg_out = 0.d0

   call fourwf(1,tmp_real,recip_input,fftg_out,wfraug,gbound,gbound,istwf_k,&
               kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,npw,ngfft(4),ngfft(5),ngfft(6),3,&
               paral_kgb,0,weight,weight)       

   recip_input = fftg_out

 endif
     
end subroutine recip_to_real_NEW

