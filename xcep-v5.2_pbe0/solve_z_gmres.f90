!
! For each cluster, compute its z vector which is defined as 
!
! z = y*(\chi_clu + \chi_env)^(-1)
!
subroutine solve_z_gmres(natom,natom_clu,natom_env,nfft,nspin,& 
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
  integer :: isp, system_type, iter_cg, s, n1, n2, n3
  real(8) :: fourpi=4.d0*3.14159265359d0, & 
             dvol, lap(nfft), accu(nspin), pen(nspin), coeff_lap, & 
             rho_env0(nfft,nspin), &   ! electron density of the enviroment 
             coeff_norm, & 
             vec1(nfft), &      ! fake array for env_weight 
             vec2(nfft), &      ! fake array for cluster_weight
             regv(nfft,nspin), &      ! fake array for cluster_weight
             vec3(nfft,nspin,2) ! fake array for dfermi_dvks

  ! gmres 
  integer             :: nbscal, revcom, & 
                         colx, coly, colz, & 
                         i,n,nloc,m,lwork, & 
                         irc(5),icntl(8),info(3)
  real(8)             :: cntl(5), rinfo(2), & 
                         start_time, end_time
  real(8),allocatable :: work(:)
  complex*16          :: zdotc
  external               zdotc




             
  ! >>>>>>>>>> function begins <<<<<<<<<<<

  write(715,'(a)')'enter solve_z_gmres()'
  write(715,*)''
  call flush(715)

  n1 = cell_nfft(1)
  n2 = cell_nfft(2)
  n3 = cell_nfft(3)

  call print_information()

  ! regularization coeff due to Zhao-Parr method 
  coeff_norm = zp_alpha**2/fourpi/zp_coeff
  coeff_lap  = 1.d0/fourpi/zp_coeff

  dvol = ucvol/dble(nfft)

  start_time = MPI_Wtime()


  !==================================================================
  ! GMRES 
  ! User's guide: 
  ! http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.414.5967&rep=rep1&type=pdf
  !==================================================================
  m = 10
  n = nfft*nspin
  nloc = n
  lwork = m*m + m*(n+5) + 5*n + 1
  allocate(work(lwork))
  icntl(1) = 6 ! stdout for error messages
  icntl(2) = 6 ! stdout for warnings
  icntl(3) = logf ! stdout for convergence 
  icntl(4) = 0 ! 0 - no preconditioning
               ! 1 - left preconditioning
               ! 2 - right preconditioning
               ! 3 - double side preconditioning
               ! 4 - error, default set in Init
  icntl(5) = 0 ! modified Gram-Schmidt
  icntl(6) = 1 ! user supplied initial guess
  icntl(7) = -1  ! max iteration 
  icntl(8) = 1   ! compute the true residual at each restart

  cntl(1) = 1e-10  ! convergence tolerance for backward erorr
  cntl(2:5) = 0.d0 ! normoalizing factors
 
  work(1:n)     =  reshape(zvec,(/n/))  ! initial x vector    
  work(n+1:2*n) = -reshape(yvec,(/n/))  ! set b vector 


! GMRES loop 
!===========
999 continue 

  call drive_dgmres(n,nloc,m,lwork,work,irc,icntl,cntl,info,rinfo)

  revcom = irc(1)
  colx   = irc(2)
  coly   = irc(3)
  colz   = irc(4)
  nbscal = irc(5)

  select case(revcom)
  case (0) 
    write(logf,*)'normal exit.'
  case (1) 
    write(logf,*)'compute Ax'
    ! work(colz)<--A * work(colx)
    zvec = reshape(work(colx:colx+n-1),(/nfft,nspin/))
    call compute_Ax(nspin,nfft,zvec,work(colz:colz+n-1))
    write(logf,*)'zvec: ',minval(zvec),maxval(zvec)
    call flush(logf)
    goto 999
  case (2) 
    write(logf,*)'compute left preconditioning, ERROR'
    stop
  case (3) 
    write(logf,*)'compute right preconditioning, ERROR'
    stop
  case (4) 
    ! scalar product
    ! work(colz) <- work(colx) work(coly)
    write(logf,*)'compute <x,y>'
    do i=0,nbscal-1 
      work(colz+i)= zdotc(n,work(colx+i*n),1,work(coly),1) 
    enddo
    call flush(logf)
    goto 999
  endselect
  ! GMRES end ==============================



  zvec = reshape(work(1:2*n),(/nfft,nspin/))
  deallocate(work)

  if (j_atom==1) then 
    write(logf,'(a)')'left solve_z_gmres()'//NEW_LINE('a')
    call flush(logf)
  endif

  write(715,'(a)')'left solve_z_gmres()'//NEW_LINE('a')//NEW_LINE('a')
  write(715,'(a)') '==== done solve_z_gmres() ===='
  call flush(715)





!~~~~~~~~~~~~~~~~~~~~~~~~~~
! SUPPORTING FUNCTIONS
!~~~~~~~~~~~~~~~~~~~~~~~~~~
contains


  subroutine print_information()
    if (myrank==0) then 
      write(logf,*)''
      write(logf,'(a)')'   >>> enter solve_z_gmres() <<<'
      write(logf,'(a,i3)')'maxiter: ',maxiter
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




  subroutine compute_Ax(nspin,nfft,x,Ax)
    implicit none
    integer :: nspin, nfft
    real(8) :: x(nfft,nspin), & 
               Ax(nfft,nspin) 
               
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
  end subroutine compute_Ax










end subroutine solve_z_gmres







