!===========================================================================================
!
! Solve global system's v_xc, with the regularization term 
! 
! 1/2*lambda<vks-vh-vpsp,vks-vh-vpsp> + 1/2*lambda*<grad(vks-vh-vpsp),grad(vks-vh-vpsp))>
!
! The linear system to solve here is 
! 
!  [-\chi + \lambda*(I-\grad^2)] * x = b
!
! \chi is the system's KS linear response function. Instead of solving x, 
! we make the variable change 
!   
!  x = P y
!
! where P = (I-\grad^2)^(-1/2) to remove the laplacian in the linear equation 
!
!===========================================================================================

subroutine solve_total_oep_regvxc_trans(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
                           logf,ucvol,global_vks,vhart,vpsp,reg_lambda,init,dExc_dvks,vxc)

  use mpi

  implicit none

  integer :: iter_scf, reg_type, & 
             nfft, nspin, logf, init, cell_nfft(3)
  real(8) :: dExc_dvks(nfft,nspin)
  real(8) :: vxc(nfft,nspin),ucvol, & 
             qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), &
             global_vks(nfft,nspin), & 
             vpsp(nfft), & 
             vhart(nfft), & 
             reg_lambda


  ! local vars ====================
  logical :: conv,do_precond, & 
             direct_resid 

  integer :: n1,n2,n3, i, isp, iter_cg
  real(8) :: lapvxc(nfft,nspin), &  
             intvxc,intvxc_old1,intvxc_old2, & 
             cg_x(nfft,nspin),cg_b(nfft,nspin),cg_resid(nfft,nspin), & 
             cg_Ap(nfft,nspin), & 
             cg_p(nfft,nspin),  & 
             cg_z(nfft,nspin), & 
             cg_zold(nfft,nspin), & 
             px(nfft,nspin), & 
             cg_resid_old(nfft,nspin), & 
             vin(nfft,nspin), &
             vout(nfft,nspin), &
             start_time, end_time,  & 
             cg_alpha, cg_beta, dtmp, dvol


  !>>>>>>>>> function begins <<<<<<<<<!

  dvol = ucvol / dble(nfft)

  n1 = cell_nfft(1)
  n2 = cell_nfft(2)
  n3 = cell_nfft(3)

  do_precond = .false.
  direct_resid = .false.

  write(logf,'(a)') NEW_LINE('a')//'    >>> enter solve_total_oep_regvxc_trans() <<<'
  write(logf,'(a)')       'variable is transformed to remove laplaican in the linear system.'
  write(logf,'(a,es12.2)')'regularization parameter: ',reg_lambda
  write(logf,'(a,i3)')    'reg_type: ',reg_type
  write(logf,'(a)')'regularize both |x| and grad(x)'
  if (do_precond) write(logf,*)'*** precond is ON ***'
  if (direct_resid) then 
    write(logf,'(a)')'direct_resid = T, resid is calculated exactly for new x'
  else 
    write(logf,'(a)')'direct_resid = F, cg recursion is used for new resid'
  endif 


  if (nspin==1) then 
    write(logf,'(a,2e12.4,a,es12.4)')'dExc/dvks: ', & 
    minval(dExc_dvks(:,1)),maxval(dExc_dvks(:,1)),' avg: ',sum(dExc_dvks(:,1))/nfft
  endif 
  if (nspin==2) then 
    write(logf,'(a,2e12.4,a,es12.4)')'dExc/dvks(spin-up): ', & 
    minval(dExc_dvks(:,1)),maxval(dExc_dvks(:,1)),' avg: ',sum(dExc_dvks(:,1))/nfft
    write(logf,'(a,2e12.4,a,es12.4)')'dExc/dvks(spin-dn): ', & 
    minval(dExc_dvks(:,2)),maxval(dExc_dvks(:,2)),' avg: ',sum(dExc_dvks(:,2))/nfft
  endif 
  call flush(logf)

  ! test 
  if (reg_type/=3) then 
    write(logf,*)'solve_total_oep_regvxc_trans() is not coded for reg_type=1,2 stop '
    stop
  endif 
  
  ! make b =====================
  do isp=1,nspin 
    vin(:,isp) = global_vks(:,isp)-vpsp-vhart
  enddo
  !
  ! compute Ax with A = -\chi * (lambda*f_H - lambda*f_H*\lap)
  !
  call apply_A(vin,vout) 
  cg_b = - dExc_dvks - vout

  ! apply P = inv[sqrt(I-\laplacian)] to cg_b 
  do isp=1,nspin 
    call apply_trans(1,n1,n2,n3,nfft,qvec,cg_b(:,isp),px(:,isp))
  enddo 
  cg_b = px

 
  !==========================
  ! set up cg_resid and cg_p
  !
  if (init<0) then
    init = 1
    cg_x = 0.d0
    cg_resid = cg_b

    ! precondition
    if (do_precond) then 
      do isp=1,nspin
        call precond_teter_inv(n1,n2,n3,nfft,qvec,1.d0,cg_resid(:,isp),cg_z(:,isp))
      enddo
    endif 
  else
    ! transform vxc to cg_x as cg_x = P^(-1) * vxc
    do isp=1,nspin 
      call apply_trans(-1,n1,n2,n3,nfft,qvec,vxc(:,isp),cg_x(:,isp))
    enddo
    call apply_chi_reg(cg_x,cg_Ap)
    cg_resid = cg_b - cg_Ap

    ! precondition 
    if (do_precond) then 
      do isp=1,nspin
        call precond_teter_inv(n1,n2,n3,nfft,qvec,1.d0,cg_resid(:,isp),cg_z(:,isp))
      enddo
    endif 
  endif 

  if (do_precond) then 
    cg_p = cg_z 
  else 
    cg_p = cg_resid
  endif 


  iter_cg = 0
  cg_alpha = 0.0d0
  cg_beta  = 0.0d0 
  start_time = MPI_Wtime()
  write(logf,'(a,2e12.4)')'initial cg_x: ',minval(cg_x(:,1)),maxval(cg_x(:,1))


  !============================
  ! cg loop 
  !============================
  do while (.true.) 
    end_time = MPI_Wtime()
    call iter_info() 
    start_time = MPI_Wtime()

    call conv_test(conv)
    if (conv) exit 
    iter_cg = iter_cg +  1

    call apply_chi_reg(cg_p,cg_Ap)

    ! precondition? 
    if (do_precond) then 
      cg_alpha = sum(cg_resid*cg_z)/sum(cg_p*cg_Ap)
    else 
      cg_alpha = sum(cg_resid**2)/sum(cg_p*cg_Ap)
    endif 

    cg_x = cg_x + cg_alpha*cg_p


!!!!!!!!!!!!!
!  do isp=1,nspin 
!    call apply_trans(1,n1,n2,n3,nfft,qvec,cg_x(:,isp),px(:,isp))
!  enddo 
!  write(logf,*)'px: ',minval(px(:,1)),maxval(px(:,1))

    cg_resid_old = cg_resid 

    !==============================
    ! compute resid exactly?
    !==============================
    if (direct_resid) then 
      ! direct compute new resid 
      call apply_chi_reg(cg_x,cg_resid)
      cg_resid = cg_b - cg_resid
    else 
      cg_resid = cg_resid - cg_alpha*cg_Ap
    endif 


    ! precondition?
    cg_zold = cg_z 
    if (do_precond) then 
      do isp=1,nspin
        call precond_teter_inv(n1,n2,n3,nfft,qvec,1.d0,cg_resid(:,isp),cg_z(:,isp))
      enddo
      cg_beta = sum(cg_resid*cg_z)/sum(cg_resid_old*cg_zold) !! Fletcher–Reeves
      cg_p    = cg_z + cg_beta*cg_p
    else 
      cg_beta = sum(cg_resid**2)/sum(cg_resid_old**2) !! Fletcher–Reeves
      cg_p    = cg_resid + cg_beta*cg_p
    endif 

  enddo 
  ! end of cg =====================
 


  ! write info 
  if (nspin==1) then 
    write(logf,'(a,2es14.6)')'cg_x: ',minval(cg_x(:,1)),maxval(cg_x(:,1))
  endif 
  if (nspin==2) then 
    write(logf,'(a,2es14.6)')'cg_x (up): ',minval(cg_x(:,1)),maxval(cg_x(:,1))
    write(logf,'(a,2es14.6)')'cg_x (dn): ',minval(cg_x(:,2)),maxval(cg_x(:,2))
  endif 

  ! transform cg_x to get the final result 
  write(logf,'(a)') 'transform cg_x back to get vxc.'
  do isp=1,nspin 
    call apply_trans(1,n1,n2,n3,nfft,qvec,cg_x(:,isp),vxc(:,isp))
  enddo 

  ! write info 
  if (nspin==1) then 
    write(logf,'(a,2es14.6)')'vxc: ',minval(vxc(:,1)),maxval(vxc(:,1))
  endif 
  if (nspin==2) then 
    write(logf,'(a,2es14.6)')'vxc (up): ',minval(vxc(:,1)),maxval(vxc(:,1))
    write(logf,'(a,2es14.6)')'vxc (dn): ',minval(vxc(:,2)),maxval(vxc(:,2))
  endif 
  call flush(logf)




!======================================
! supporting subroutines
!======================================
contains 

  subroutine iter_info 
    dtmp = sum(cg_resid**2)*ucvol/nfft
    intvxc_old2 = intvxc_old1
    intvxc_old1 = intvxc
    write(logf,'(a,i3,a,es8.2,a,2es11.2,a,f6.1,a)') & 
     "",iter_cg, & 
     "  res: ",dtmp,'   cg_x(up):',minval(cg_x(:,1)),maxval(cg_x(:,1)), & 
     '   t: ',end_time-start_time,' s'
    call flush(logf)
  end subroutine 


  subroutine conv_test(conv)
    logical :: conv
    real(8) :: dtmp
    
    conv = .false. 
    dtmp = sum(cg_resid**2)*ucvol/nfft

    !if (iter_cg==20 .or. (iter_cg>=10 .and. dtmp<1e-11)) then 
    if (iter_cg==20) then 
      !!write(logf,'(a)')'iter_cg=20 or resid<1e-11. done'
      write(logf,'(a)')'iter_cg=20. done'
      conv = .true. 
      call flush(logf)
    endif 
  end subroutine conv_test


  !============================================
  ! compute A|vin> = |vout> with A defined as
  ! 
  ! A = -\chi * (lambda*f_H - lambda*f_H*\lap)
  !
  subroutine apply_A(vin, vout)
    implicit none 
    integer :: system_type,isp
    real(8) :: vin(nfft,nspin), & 
               vout(nfft,nspin), & 
               avec(nfft,nspin), & 
               vh(nfft), etmp, & 
               lap(nfft,nspin)

    ! compute laplacian of vin
    do isp=1,nspin 
      call laplacian(nfft,n1,n2,n3,vin(:,isp),qvec,lap(:,isp))
    enddo 

    do isp=1,nspin
      call hartree(n1,n2,n3,1,ucvol,qvec,-vin(:,isp)+lap(:,isp),vh,etmp) 
      avec(:,isp) = reg_lambda*vh
    enddo

    system_type = 0 ! global system 
    call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,avec,vout)
    call remove_mean(nspin,nfft,vout)
  end subroutine



  
  !==========================================
  ! compute [P\chiP + lamdba*I] x
  ! where P is the transformation matrix 
  !
  subroutine apply_chi_reg(vin, vout) 
    implicit none 
    real(8),intent(in)  :: vin(nfft,nspin)
    real(8),intent(out) :: vout(nfft,nspin)
    ! local vars 
    integer :: system_type,isp
    real(8) :: lap(nfft,nspin),px(nfft,nspin)

    ! apply P to x 
    do isp=1,nspin 
      call apply_trans(1,n1,n2,n3,nfft,qvec,vin(:,isp),px(:,isp))
    enddo

    system_type = 0 ! global system 
    call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,px,vout)
    call remove_mean(nspin,nfft,vout)
    vout = -vout  ! KS response is negative definite  

    ! apply P to x 
    do isp=1,nspin 
      call apply_trans(1,n1,n2,n3,nfft,qvec,vout(:,isp),px(:,isp))
    enddo

    vout = px + reg_lambda*vin
  endsubroutine 



  !===================================================
  ! compute Px where P is the transformation matrix 
  !
  ! P = inv[sqrt[I-\laplacian]]
  ! 
  subroutine apply_trans(dir,n1,n2,n3,nfft,qvec,x,px)
    implicit none 

    integer :: n1,n2,n3,nfft,dir
    real(8) :: x(nfft), px(nfft), & 
               qvec(3,(n1/2+1)*n2*n3)
    
    ! local vars 
    integer :: dim1, ix, iy, iz
    real(8) :: q3d(3,n1/2+1,n2,n3), & 
               rtmp(n1,n2,n3),qq

    complex(kind=8) :: fft1(n1/2+1,n2,n3), & 
                       fft2(n1/2+1,n2,n3)

    dim1=(n1/2+1)
    call FFT(n1,n2,n3,reshape(x,(/n1,n2,n3/)),fft1,1)

    ! convert qvec to 3-dimension 
    q3d(1,:,:,:) = reshape(qvec(1,:),(/dim1,n2,n3/))
    q3d(2,:,:,:) = reshape(qvec(2,:),(/dim1,n2,n3/))
    q3d(3,:,:,:) = reshape(qvec(3,:),(/dim1,n2,n3/))

    do ix=1,dim1
      do iy=1,n2
        do iz=1,n3
          qq = q3d(1,ix,iy,iz)**2 + & 
               q3d(2,ix,iy,iz)**2 + & 
               q3d(3,ix,iy,iz)**2
          
          if (dir == 1)  fft2(ix,iy,iz) = fft1(ix,iy,iz)*1.d0/sqrt(1.d0+qq)   ! P*x
          if (dir == -1) fft2(ix,iy,iz) = fft1(ix,iy,iz)/(1.d0/sqrt(1.d0+qq)) ! P^(-1)*x
        enddo
      enddo
    enddo

    call FFT(n1,n2,n3,rtmp,fft2,-1)
    px = reshape(rtmp,(/n1*n2*n3/))

  endsubroutine apply_trans

end subroutine 

