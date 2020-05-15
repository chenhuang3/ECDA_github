!
! For each cluster, compute its z vector which is defined as 
!
! z = y*(\chi_clu + \chi_env)^(-1)
!
! solve z by L-BFGS optimization, line search is exact since the systme is quadratic 
!
subroutine solve_z_lbfgs_exlsr(natom,natom_clu,natom_env,nfft,nspin,& 
                   maxiter,j_atom,myrank,logf,cell_nfft,ucvol, &
                   qvec,yvec,vks_clu,vks_env,sub_rhor, & 
                   zp_coeff,zp_alpha,zvec_init,zvec)

  use comm_data
  use mpi 
  use interface_funcs, only: calc_subsystem_dfpt

  implicit none 

  integer :: nfft, nspin, j_atom, maxiter, & 
             natom, natom_clu, natom_env, & 
             myrank, logf, cell_nfft(3), zvec_init

  real(8) :: zp_alpha, zp_coeff, & 
             sub_rhor(nfft,nspin,2), & 
             yvec(nfft,nspin),  & 
             vks_clu(nfft,nspin), & 
             vks_env(nfft,nspin), & 
             qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), & 
             zvec(nfft,nspin), & 
             ucvol 


  ! local vars 
  integer :: iter_newx, & 
             isp, system_type, iter_cg, s, n1, n2, n3
  real(8) :: fourpi=4.d0*3.14159265359d0, & 
             lap(nfft), accu(nspin), pen(nspin), coeff_lap, & 
             rho_env0(nfft,nspin), &   ! electron density of the enviroment 
             coeff_norm, & 
             Ad(nfft,nspin), & ! fake array 
             vec1(nfft), &      ! fake array for env_weight 
             vec2(nfft), &      ! fake array for cluster_weight
             vec3(nfft,nspin,2) ! fake array for dfermi_dvks

  ! cg ------------------------
  logical :: conv, do_precond
  real(8) :: cg_b(nfft,nspin), & 
             cg_Az(nfft,nspin), & 
             dir_tn(nfft,nspin), & 
             g(nfft,nspin), & 
             start_time, end_time, & 
             objfun,  dstep, & 
             cg_alpha, cg_beta, dtmp, d_resid, dvol

  ! LBFGS 
  integer,parameter   :: mhist = 5
  integer             :: iter, j
  real(8)             :: ghist(nspin*nfft,mhist+1), & 
                         xhist(nspin*nfft,mhist+1)

             
  ! >>>>>>>>>> function begins <<<<<<<<<<<

  write(715,'(a)')'enter solve_z_lbfgs_exlsr()'
  write(715,*)''
  call flush(715)

  n1 = cell_nfft(1)
  n2 = cell_nfft(2)
  n3 = cell_nfft(3)

  call print_information()

  do_precond = .false.

  ! regularization coeff due to Zhao-Parr method 
  coeff_norm = zp_alpha**2/fourpi/zp_coeff
  coeff_lap  = 1.d0/fourpi/zp_coeff

  dvol = ucvol/dble(nfft)

  cg_b = -yvec  ! since KS linear response is negative definite, 
                ! we in fact solve -\chi x = -yvec

  !========================
  ! My L-BFGS   
  !========================
  do iter=1,maxiter

    !====================================================================
    !
    ! NOTE: We must re-compute Az here, even if the computational cost 
    ! will be doubled by doing so. In principle, we do not need to do this;
    ! however, for very corrlated system, the solutions for clusters'  
    ! Sternheimer equations are not very accurate and 
    ! this error will accumulate and cause significant
    ! error in the solution (zvec). 
    !
    ! It caused me weeks to figure this out. Therefore, we need to 
    ! compute A*zvec here and also in the conjugate gradinent based 
    ! solve_z() code. 
    !
    call compute_Ax(nspin,nfft,zvec,cg_Az)

    ! objective funciton 
    objfun = sum(zvec*cg_Az)*0.5d0*dvol - sum(cg_b*zvec)*dvol

    ! gradient 
    g = cg_Az - cg_b

    call print_info(iter,objfun,g,zvec,dstep)

    ! L-BFGS direction 
    !====================
    call lbfgs_dir(nspin,nfft,mhist,iter,g,zvec,ghist,xhist,dir_tn)
    
    ! compute exact step size 
    !=========================
    call compute_Ax(nspin,nfft,dir_tn,Ad)

    ! Since it is a quadratic function, the step size 
    ! can be calculated exactly. This is just the step size (alpha) 
    ! used in the CG method, which is
    !
    !  alpha = <p,b-Ax>/<p,Ap>
    !
    dstep = sum(dir_tn*(-g))/sum(dir_tn*Ad)
    zvec = zvec + dir_tn*dstep

    ! make use of the Ad to reduce one evaluation of Ax
    !!cg_Az = cg_Az + Ad*dstep
  enddo 

  if (j_atom==1) then 
    write(logf,'(a)')'left solve_z_lbfgs_exlsr()'//NEW_LINE('a')
    call flush(logf)
  endif

  write(715,'(a)')'left solve_z_lbfgs_exlsr()'//NEW_LINE('a')//NEW_LINE('a')
  write(715,'(a)') '==== done solve_z_lbfgs_exlsr() ===='
  call flush(715)





!~~~~~~~~~~~~~~~~~~~~~~~~~~
! SUPPORTING FUNCTIONS
!~~~~~~~~~~~~~~~~~~~~~~~~~~
contains

  subroutine compute_Ax(nspin,nfft,x,Ax)
    implicit none
    integer :: nspin, nfft
    real(8) :: x(nfft,nspin), & 
               Ax(nfft,nspin) 
               

    ! compute object function and gradient 
    !=====================================
    call calc_subsystem_dfpt(natom,natom_clu,natom_env,j_atom,nfft,nspin,ucvol, & 
                             vec1,vec2,qvec,sub_rhor,cell_nfft, & 
                             vec3,vks_clu,vks_env,x,Ax)

    Ax = -Ax  ! since KS response function is negative definite 
    call remove_mean(nspin,nfft,Ax)

    ! include regularization due to Zhao-Parr method 
    if (dfet_pen_type==1) then 
      call regularization(1,coeff_norm,nspin,nfft,cell_nfft,qvec,Ax,x)
      call regularization(2,coeff_lap, nspin,nfft,cell_nfft,qvec,Ax,x)
    endif 
    if (dfet_pen_type==2) then
      call regularization(1,1.d0/zp_coeff,nspin,nfft,cell_nfft,qvec,Ax,x)
    endif 

  end subroutine 




  subroutine print_information()
    if (myrank==0) then 
      write(logf,*)''
      write(logf,'(a)')'   >>> enter solve_z_lbfgs_exlsr() <<<'
      write(logf,'(a)')'L-BFGS for solving zvec, Line search is exact.'
      write(logf,'(a,i3)')'dfet_pen_type: ',dfet_pen_type
      write(logf,'(a,L2)')'do_envOF: ',do_envOF
      write(logf,'(a,es8.1,a,f10.3)')& 
        'zhao-parr reg. is considered, zp_coeff: ',zp_coeff, '  zp_alpha: ',zp_alpha
      write(logf,'(a,i3)')'maxiter: ',maxiter
      if (nspin==1) write(logf,'(a,2es12.4,a,es10.2)') & 
        'yvec: ',minval(yvec(:,1)),maxval(yvec(:,1))," avg(y):",sum(yvec)/dble(nfft)
      if (nspin==2) then 
        write(logf,'(a,2es14.4,a,es10.2)') & 
          'yvec(up): ',minval(yvec(:,1)),maxval(yvec(:,1)),' avg(y): ',sum(yvec(:,1))/dble(nfft)
        write(logf,'(a,2es14.4,a,es10.2)') & 
        'yvec(dn): ',minval(yvec(:,2)),maxval(yvec(:,2)),' avg(y): ',sum(yvec(:,2))/dble(nfft)
      endif 
      call flush(logf)
    endif 
  end subroutine 



  subroutine print_info(iter,objfun,g,z,dstep) 
    integer :: iter
    real(8) :: objfun, dstep, & 
               g(nfft,nspin), & 
               z(nfft,nspin)

    ! print information 
    if (myrank==0) then 
      write(logf,'(a,i3,a,es14.6,a,es10.2,a,2es11.3,a,es8.1)') & 
       'it:',iter,'  F:',objfun,'  gmax:',maxval(abs(g)), & 
        '  z(up):',minval(z(:,1)),maxval(z(:,1)),'  stp:',dstep
      call flush(logf)


      write(715,'(a,i3,a,es14.6,a,es10.2,a,2es11.3,a,es8.1)') & 
        'it:',iter,'  F:',objfun,'  gmax:',maxval(abs(g)), & 
        '  z(up):',minval(z(:,1)),maxval(z(:,1)),'  stp:',dstep
      call flush(715)

    endif 
  end subroutine 


  subroutine print_cg_header()
    ! print header 
    if (j_atom==1) then 
      if (do_precond) write(logf,'(a)')'*** precond is ON ***'
      if (.not. do_precond) write(logf,'(a)')'*** precond is OFF ***'
      if (nspin==2) & 
        write(logf,'(a)')' iter  |resid|  cg_alpha  cg_beta  min/max(cg_x)(spin up)  avg(cg_x)  time(s)'
      if (nspin==1) &
        write(logf,'(a)')' iter  |resid|  cg_alpha  cg_beta      min/max(cg_x)       avg(cg_x)  time(s)'
      call flush(logf)
    endif 
    if (nspin==2) & 
      write(715,'(a)')' iter  |resid|  cg_alpha  cg_beta  min/max(cg_x)(spin up)  avg(cg_x)  time(s)'
    if (nspin==1) &
      write(715,'(a)')' iter  |resid|  cg_alpha  cg_beta      min/max(cg_x)       avg(cg_x)  time(s)'
    call flush(715)
  endsubroutine 



end subroutine 







