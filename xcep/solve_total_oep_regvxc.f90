!===========================================================================================
!
! Solve global system's v_xc, with the regularization term 
! 
! 1/2*lambda<vks-vh-vpsp,vks-vh-vpsp> + 1/2*lambda*<grad(vks-vh-vpsp),grad(vks-vh-vpsp))>
!
!===========================================================================================

subroutine solve_total_oep_regvxc(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
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
  logical :: conv, & 
             do_precond 

  integer :: n1,n2,n3, i, isp, system_type, iter_cg
  real(8) :: lapvxc(nfft,nspin), &  
             intvxc,intvxc_old1,intvxc_old2
  real(8) :: cg_x(nfft,nspin),cg_b(nfft,nspin),cg_resid(nfft,nspin), & 
             cg_Ap(nfft,nspin), & 
             cg_p(nfft,nspin),  & 
             cg_z(nfft,nspin), & 
             cg_zold(nfft,nspin), & 
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

  write(logf,'(a)') NEW_LINE('a')//'    >>> enter solve_total_oep_regvxc() <<<'
  write(logf,'(a,es12.2)')'regularization parameter: ',reg_lambda
  write(logf,'(a,i3)')'reg_type: ',reg_type
 
  system_type = 0 ! global system 

  select case(reg_type) 
    case(1); write(logf,'(a)')'regularize |x|'
    case(2); write(logf,'(a)')'regularize grad(x)'
    case(3); write(logf,'(a)')'regularize both |x| and grad(x)'
  end select

  if (do_precond) write(logf,*)'*** precond is ON ***'


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

  
  ! make b =====================
  do isp=1,nspin 
    vin(:,isp) = global_vks(:,isp)-vpsp-vhart
  enddo
  call apply_A(vin,vout) 
  cg_b = - dExc_dvks - vout


  !==========================
  ! set up cg_resid and cg_p
  !==========================
  if (init<0) then
    cg_x = 0.d0
    cg_resid = cg_b

    ! precondition
    if (do_precond) then 
      do isp=1,nspin
        call precond_teter_inv(n1,n2,n3,nfft,qvec,1.d0,cg_resid(:,isp),cg_z(:,isp))
      enddo
    endif 

    init = 1
  else 
    cg_x = vxc
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

  !~~~~~~~~~~~~~~~~~~~~~~~~
  ! cg loop 
  !~~~~~~~~~~~~~~~~~~~~~~~~
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

    cg_resid_old = cg_resid 
    cg_resid     = cg_resid - cg_alpha*cg_Ap

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

  enddo ! end of cg 
 
  vxc = cg_x 





!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SUPPORTING SUBROUTINES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    if (iter_cg==20 .or. (iter_cg>=10 .and. dtmp<1e-9)) then 
      write(logf,'(a)')'iter_cg=20 or resid<1e-10. done solve_total_oep().'
      if (nspin==1) then 
        write(logf,'(a,2e12.4)')'cg_x: ', & 
        minval(cg_x(:,1)),maxval(cg_x(:,1))
      endif 
      if (nspin==2) then 
        write(logf,'(a,2e12.4)')'cg_x (up): ', & 
        minval(cg_x(:,1)),maxval(cg_x(:,1))
        write(logf,'(a,2e12.4)')'cg_x (dn): ', & 
        minval(cg_x(:,2)),maxval(cg_x(:,2))
      endif 
      conv = .true. 
      call flush(logf)
    endif 
  end subroutine conv_test


  !============================================
  ! compute A|vin> = |vout> with A defined as
  ! 
  ! A = -\chi * (lambda*f_H - lambda*f_H*\lap)
  !============================================
  subroutine apply_A(vin, vout)
    implicit none 
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
      if (reg_type==1) then
        call hartree(n1,n2,n3,1,ucvol,qvec,-vin(:,isp),vh,etmp) 
        avec(:,isp) = reg_lambda*vh
      elseif (reg_type==2) then
        call hartree(n1,n2,n3,1,ucvol,qvec,lap(:,isp),vh,etmp)
        avec(:,isp) = reg_lambda*vh
      elseif (reg_type==3) then 
        call hartree(n1,n2,n3,1,ucvol,qvec,-vin(:,isp)+lap(:,isp),vh,etmp) 
        avec(:,isp) = reg_lambda*vh
      endif 
    enddo

    system_type = 0 ! global system 
    call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,avec,vout)
    call remove_mean(nspin,nfft,vout)
  end subroutine




  subroutine apply_chi_reg(vin, vout) 
    implicit none 
    real(8), intent(in) :: vin(nfft,nspin)
    real(8), intent(out) :: vout(nfft,nspin)
    ! local vars 
    real(8) :: lap(nfft,nspin)

    system_type = 0 ! global system 
    call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,vin,vout)
    call remove_mean(nspin,nfft,vout)
    vout = -vout  ! KS response is negative definite  

    ! compute laplacian of vin
    do isp=1,nspin 
      call laplacian(nfft,n1,n2,n3,vin(:,isp),qvec,lap(:,isp))
    enddo 

    if (reg_type==1) vout = vout + reg_lambda*vin
    if (reg_type==2) vout = vout - reg_lambda*lap
    if (reg_type==3) vout = vout - reg_lambda*lap + reg_lambda*vin
  endsubroutine 

end subroutine solve_total_oep_regvxc

