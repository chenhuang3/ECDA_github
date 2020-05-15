!=========================================
! solve for global system's v_xc
!=========================================

subroutine solve_total_oep_lbfgs(myrank,iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
                                 logf,ucvol,global_vks,reg_lambda,init,dExc_dvks,vxc)

  use mpi

  implicit none

  integer :: iter_scf, reg_type, myrank, & 
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
  real(8) :: cg_x(nfft,nspin), & 
             cg_b(nfft,nspin), & 
             cg_Ap(nfft,nspin), & 
             start_time, end_time,  & 
             grad(nfft,nspin), & 
             cg_beta, dtmp, dvol, ee

 ! L-BFGS from Nocedal and Zhu
 integer,parameter   :: mhist = 10
 logical             :: lsave(4)
 character*60        :: task, csave
 integer             :: isave(44)
 real(8)             :: dsave(29), gmax, gnorm
 integer :: lbfgs_nbd(nfft*nspin), & 
            lbfgs_iwa(3*nfft*nspin)
 real    :: lbfgs_time 
 real(8) :: lbfgs_l(nfft*nspin), & 
            lbfgs_u(nfft*nspin), & 
            lbfgs_wa(2*mhist*nfft*nspin+4*nfft*nspin+12*mhist*mhist+12*mhist)

  !>>>>>>>>> function begins <<<<<<<<<!

  dvol = ucvol / dble(nfft)

  write(logf,'(a)') NEW_LINE('a')//'    >>> enter solve_total_oep_lbfgs() <<<'
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
    init = 1
    cg_Ap = 0.d0
  else 
    cg_x = vxc
    call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,cg_x,cg_Ap)
    cg_Ap = -cg_Ap  ! KS response is negative definite 
    call remove_mean(nspin,nfft,cg_Ap)
    call regularization(reg_type,reg_lambda,nspin,nfft,cell_nfft,qvec,cg_Ap,cg_x)
  endif 

  grad = -cg_b+cg_Ap
  ee = 0.5d0*sum(cg_x*cg_Ap)*dvol - sum(cg_b*cg_x)*dvol
  task='START'

!
! ee = 1/2<x,\chi,x> - <cg_b,x> + 1/2*reg_vxc*<x,x>
! grad = \chi x - cg_b + reg_vxc*x
!

9999 continue

  lbfgs_l = 0.d0 
  lbfgs_u = 0.d0 
  lbfgs_nbd = 0    ! unbounded 

  call setulb(nfft*nspin, mhist, cg_x, lbfgs_l, lbfgs_u, lbfgs_nbd,  & 
              ee, grad, 1.d1, 0.d0, lbfgs_wa, lbfgs_iwa, & 
              task, 99, csave, lsave, isave, dsave, lbfgs_time, myrank)
   

  if (task(1:2) == 'FG') then 
    write(logf,'(a)')'L-BFGS returned FG'
    call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,cg_x,cg_Ap)
    cg_Ap = -cg_Ap  ! KS response is negative definite 
    call remove_mean(nspin,nfft,cg_Ap)
    call regularization(reg_type,reg_lambda,nspin,nfft,cell_nfft,qvec,cg_Ap,cg_x)
    grad = -cg_b+cg_Ap
    ee = 0.5d0*sum(cg_x*cg_Ap)*dvol - sum(cg_b*cg_x)*dvol

    ! display
    dtmp = sum(grad**2)*ucvol/nfft
    intvxc_old2 = intvxc_old1
    intvxc_old1 = intvxc
    write(logf,'(a,f14.6,a,es8.2,a,2es13.4,a,f8.1)') & 
     "ee: ",ee,"  resid: ",dtmp,'   cg_x(up):',minval(cg_x(:,1)),maxval(cg_x(:,1)), & 
     '   t(s): ',end_time-start_time
    call flush(logf)
    goto 9999
  elseif (task(1:5) == 'NEW_X') then 
    write(logf,'(a)')'L-BFGS returned NEW_X'
    goto 9999
  else
    write(logf,*)'code stopped due to unknow message from L-BFGS'
    write(logf,*)'task: ',trim(task)
  endif 

 
  vxc = cg_x 




end subroutine solve_total_oep_lbfgs

