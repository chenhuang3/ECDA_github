!
! For each cluster, compute its z vector which is defined as 
!
! z = y*(\chi_clu + \chi_env)^(-1)
!
subroutine solve_z_lbfgs(natom,natom_clu,natom_env,nfft,nspin,& 
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
             vec1(nfft), &      ! fake array for env_weight 
             vec2(nfft), &      ! fake array for cluster_weight
             regv(nfft,nspin), &      ! fake array for cluster_weight
             vec3(nfft,nspin,2) ! fake array for dfermi_dvks

  ! cg ------------------------
  logical :: conv, do_precond
  real(8) :: cg_b(nfft,nspin), & 
             cg_z(nfft,nspin), & 
             cg_Az(nfft,nspin), & 
             g(nfft,nspin), & 
             start_time, end_time, & 
             objfun, & 
             cg_alpha, cg_beta, dtmp, d_resid, dvol

  ! LBFGS 
  integer,parameter   :: mhist = 5
  logical             :: lsave(4)
  character*60        :: task, csave
  integer             :: isave(44), nfg
  real(8)             :: dsave(29), gmax, gnorm
  real                :: lbfgs_time 
  integer   :: lbfgs_nbd(nfft*nspin), lbfgs_iwa(3*nfft*nspin)
  real(8)   :: lbfgs_l(nfft*nspin),  & 
               lbfgs_u(nfft*nspin),  &  
               lbfgs_wa(2*mhist*nfft*nspin+4*nfft*nspin+12*mhist*mhist+12*mhist)

             
  ! >>>>>>>>>> function begins <<<<<<<<<<<

  write(715,'(a)')'enter solve_z_lbfgs()'
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

  cg_z =  zvec
  cg_b = -yvec  ! since KS linear response is negative definite, 
                ! we in fact solve -\chi x = -yvec

  !========================
  ! L-BFGS   
  !========================
  task = 'START'
  lbfgs_l = 0.d0 
  lbfgs_u = 0.d0 
  lbfgs_nbd = 0    ! unbounded 
  iter_newx = 0
  nfg = 0


9999 continue 

  call setulb(nfft*nspin, mhist, cg_z, lbfgs_l, lbfgs_u, lbfgs_nbd,  & 
              objfun, g, 1.d1, 0.d0, lbfgs_wa, lbfgs_iwa, & 
              task, 99, csave, lsave, isave, dsave, lbfgs_time, myrank)

  if (task(1:2) == 'FG') then 
    nfg = nfg + 1
    call calc_subsystem_dfpt(natom,natom_clu,natom_env,j_atom,nfft,nspin,ucvol, & 
                             vec1,vec2,qvec,sub_rhor,cell_nfft, & 
                             vec3,vks_clu,vks_env,cg_z,cg_Az)

    cg_Az = -cg_Az  ! since KS response function is negative definite 
    call remove_mean(nspin,nfft,cg_Az)

    ! include regularization due to Zhao-Parr method 
    if (dfet_pen_type==1) then 
      call regularization(1,coeff_norm,nspin,nfft,cell_nfft,qvec,cg_Az,cg_z)
      call regularization(2,coeff_lap, nspin,nfft,cell_nfft,qvec,cg_Az,cg_z)
    endif 
    if (dfet_pen_type==2) then
      call regularization(1,1.d0/zp_coeff,nspin,nfft,cell_nfft,qvec,cg_Az,cg_z)
    endif 

    ! objective funciton 
    objfun = sum(cg_z*cg_Az)*0.5d0*dvol - sum(cg_b*cg_z)*dvol

    ! gradient 
    g = cg_Az - cg_b

    call print_info()

    goto 9999

  elseif (task(1:5) == 'NEW_X') then 
    nfg = 0
    iter_newx = iter_newx + 1
    if (iter_newx==20) then 
      write(logf,'(a)')'iter_new_x = 20, exit L-BFGS.'
      call flush(logf)
      goto 10000
    endif 
    goto 9999

  else
    write(logf,*)'code stopped due to unknow message from L-BFGS'
    write(logf,*)'task: ',trim(task)

  endif 
  !============================= 
  ! end of L-BFGS 
  !=============================

10000 continue 

  zvec = cg_z

  if (j_atom==1) then 
    write(logf,'(a)')'left solve_z_lbfgs()'//NEW_LINE('a')
    call flush(logf)
  endif

  write(715,'(a)')'left solve_z_lbfgs()'//NEW_LINE('a')//NEW_LINE('a')
  write(715,'(a)') '==== done solve_z_lbfgs() ===='
  call flush(715)





!~~~~~~~~~~~~~~~~~~~~~~~~~~
! SUPPORTING FUNCTIONS
!~~~~~~~~~~~~~~~~~~~~~~~~~~
contains


  subroutine print_information()
    if (myrank==0) then 
      write(logf,*)''
      write(logf,'(a)')'   >>> enter solve_z_lbfgs() <<<'
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

  subroutine print_info
    ! print information 
    if (myrank==0) then 
      if (nfg==1) then 
        write(logf,'(a,es14.6,a,es10.2,a,2es11.3,a,es9.1)') & 
        'FG=> F:',objfun,'  gmax:',maxval(abs(g)), & 
        '  z:',minval(cg_z),maxval(cg_z), ' dstep: ',dsave(14)
      else
        write(logf,'(a,es14.6,a,es10.2,a,2es11.3,a,es9.1,a)') & 
        'FG=> F:',objfun,'  gmax:',maxval(abs(g)), & 
        '  z:',minval(cg_z),maxval(cg_z), ' dstep: ',dsave(14),' (LS)'
      endif 
      call flush(logf)


      if (nfg==1) then 
        write(715,'(a,es14.6,a,es10.2,a,2es11.3,a,es9.1)') & 
        'FG=> F:',objfun,'  gmax:',maxval(abs(g)), & 
        '  z:',minval(cg_z),maxval(cg_z), ' dstep: ',dsave(14)
      else
        write(715,'(a,es14.6,a,es10.2,a,2es11.3,a,es9.1,a)') & 
        'FG=> F:',objfun,'  gmax:',maxval(abs(g)), & 
        '  z:',minval(cg_z),maxval(cg_z), ' dstep: ',dsave(14),' (LS)'
      endif 
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







