!=========================================
!
! solve for global system's v_xc
!
!=========================================

subroutine solve_total_oep(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
                           logf,ucvol,global_vks,reg_lambda,init,dExc_dvks,vxc)

  use mpi

  implicit none

  integer :: iter_scf, reg_type, & 
             nfft, nspin, logf, init, cell_nfft(3)
  real(8) :: dExc_dvks(nfft,nspin)
  real(8) :: vxc(nfft,nspin),ucvol, & 
             qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), &
             global_vks(nfft,nspin), & 
             reg_lambda

  ! local vars 
  logical :: conv
  integer :: i, isp, system_type, iter_cg
  real(8) :: intvxc,intvxc_old1,intvxc_old2
  real(8) :: cg_x(nfft,nspin),cg_b(nfft,nspin),cg_resid(nfft,nspin), & 
             cg_Ap(nfft,nspin), & 
             cg_p(nfft,nspin), &
             cg_resid_old(nfft,nspin), & 
             start_time, end_time,  & 
             cg_alpha, cg_beta, dtmp, dvol


  !>>>>>>>>> function begins <<<<<<<<<!

  dvol = ucvol / dble(nfft)

  write(logf,'(a)') NEW_LINE('a')//'    >>> enter solve_total_oep() <<<'
  write(logf,'(a,es12.2)')'regularization parameter: ',reg_lambda
 
  system_type = 0 ! global system 

  select case(reg_type) 
    case(1); write(logf,'(a)')'regularize |x|'
    case(2); write(logf,'(a)')'regularize grad(x)'
    case(3); write(logf,'(a)')'regularize both |x| and grad(x)'
  end select


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

  
  cg_b = -dExc_dvks
  ! clean up cg_b, otherwise the CG cannot be solved
  !call remove_mean(nspin,nfft,cg_b)

  if (init<0) then
    cg_x = 0.d0
    cg_resid = cg_b
    init = 1
  else 
    cg_x = vxc
    call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,cg_x,cg_Ap)
    call remove_mean(nspin,nfft,cg_Ap)
    cg_Ap = -cg_Ap  ! KS response is negative definite 
    call regularization(reg_type,reg_lambda,nspin,nfft,cell_nfft,qvec,cg_Ap,cg_x)
    cg_resid = cg_b - cg_Ap
  endif 
  cg_p = cg_resid


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

    call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,cg_p,cg_Ap)
    call remove_mean(nspin,nfft,cg_Ap)
    cg_Ap = -cg_Ap  ! since KS response function is negative definite 

    !! >>> checking  <<<
    !if (abs(sum(cg_Ap)*dvol)>1e-3) then 
    !   write(logf,*) &
    !   'abs(sum(cg_Ap)*dvol)>1e-3, stop! abs(sum(cg_Ap)*dvol): ',abs(sum(cg_Ap)*dvol)
    !   stop
    !endif 
    if (sum(cg_p*cg_Ap)*dvol<0.d0) then 
       write(logf,*)'WARNING! <p,Ap>:',sum(cg_p*cg_Ap)*dvol,' <0 '
       call flush(logf)
       if ( abs(sum(cg_p*cg_Ap)*dvol)>1e-5 ) then 
         write(logf,'(a)')'ERROR: abs(sum(cg_p*cg_Ap)*dvol) > 1e-5, STOP'
         call flush(logf)
         stop
       endif  
    endif 

    call regularization(reg_type,reg_lambda,nspin,nfft,cell_nfft,qvec,cg_Ap,cg_p)

    cg_alpha = sum(cg_resid**2)/sum(cg_p*cg_Ap)
    cg_x     = cg_x + cg_alpha*cg_p

    cg_resid_old = cg_resid 
    cg_resid     = cg_resid - cg_alpha*cg_Ap

    cg_beta      = sum(cg_resid**2)/sum(cg_resid_old**2) !! Fletcher–Reeves
    cg_p         = cg_resid + cg_beta*cg_p

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
    write(logf,'(a,i3,a,es8.2,a,2es13.4,a,f8.1)') & 
     "it: ",iter_cg, & 
     "  resid: ",dtmp,'   cg_x(up):',minval(cg_x(:,1)),maxval(cg_x(:,1)), & 
     '   t(s): ',end_time-start_time
    call flush(logf)
  end subroutine 


  subroutine conv_test(conv)
    logical :: conv
    real(8) :: dtmp
    
    conv = .false. 
    dtmp = sum(cg_resid**2)*ucvol/nfft

    if (iter_cg==40 .or. dtmp<1e-9) then 
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


end subroutine solve_total_oep

