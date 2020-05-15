!
! For each cluster, compute its z vector which is defined as 
!
! z = y*(\chi_clu + \chi_env)^(-1)
!
subroutine solve_z_trans(natom,natom_clu,natom_env,nfft,nspin,& 
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
  real(8) :: pi=3.1415926d0, & 
             fourpi=4.d0*3.14159265359d0, & 
             lap(nfft), accu(nspin), pen(nspin), coeff_lap, & 
             rho_env0(nfft,nspin), &   ! electron density of the enviroment 
             coeff_norm, px(nfft,nspin),  & 
             vec1(nfft), &      ! fake array for env_weight 
             vec2(nfft), &      ! fake array for cluster_weight
             regv(nfft,nspin), &      ! fake array for cluster_weight
             vec3(nfft,nspin,2) ! fake array for dfermi_dvks

  ! cg ------------------------
  logical :: conv, do_precond
  real(8) :: cg_x(nfft,nspin), & 
             cg_b(nfft,nspin), & 
             cg_Ax(nfft,nspin), & 
             cg_z(nfft,nspin), & 
             cg_zold(nfft,nspin), & 
             cg_resid(nfft,nspin), & 
             cg_Ap(nfft,nspin), & 
             cg_p(nfft,nspin), &
             cg_resid_old(nfft,nspin), & 
             start_time, end_time, & 
             cg_alpha, cg_beta, dtmp, d_resid, dvol


             
  ! >>>>>>>>>> function begins <<<<<<<<<<<

  write(715,'(a)')'enter solve_z_trans()'
  write(715,*)''
  call flush(715)

  n1 = cell_nfft(1)
  n2 = cell_nfft(2)
  n3 = cell_nfft(3)

  call print_information()

!  if (do_envOF) then 
!    do_precond = .true. 
!  else 
    do_precond = .false.
!  endif 

  ! regularization coeff due to Zhao-Parr method 
  coeff_norm = zp_alpha**2/fourpi/zp_coeff
  coeff_lap  = 1.d0/fourpi/zp_coeff

  dvol = ucvol/dble(nfft)

  !=====================
  ! compute cg_b 
  !
  cg_b = -yvec  ! since KS linear response is negative definite, 
                ! we in fact solve -\chi x = -yvec
  ! apply transform matrix P 
  do isp=1,nspin 
    call apply_trans(1,n1,n2,n3,nfft,qvec,cg_b(:,isp),px(:,isp))
  enddo 
  cg_b = px



  if ((nspin==2) .and. (do_precond)) then
    write(logf,*)'precond_teter() doe snto with spin=2 yet, code stop'
    stop
  endif 

  if (zvec_init<0) then 
    zvec = 0.d0 
    cg_x = 0.d0
    cg_resid = cg_b
    if (do_precond) then 
      do isp=1,nspin
        call precond_teter_inv(n1,n2,n3,nfft,qvec,1.d0,cg_resid(:,isp),cg_z(:,isp))
      enddo
    endif 
    zvec_init = 1
  else 
    !=====================================
    ! use previous solution
    !
    ! transform vxc to cg_x as cg_x = P^(-1) * vxc
    do isp=1,nspin 
      call apply_trans(-1,n1,n2,n3,nfft,qvec,zvec(:,isp),cg_x(:,isp))
    enddo
    ! compute Ax
    call calc_subsystem_dfpt(natom,natom_clu,natom_env,j_atom,nfft,nspin,ucvol, & 
                             vec1,vec2,qvec,sub_rhor,cell_nfft, & 
                             vec3,vks_clu,vks_env,zvec,px)
    px = -px ! since KS response function is negative definite
    call remove_mean(nspin,nfft,cg_Ap)
    ! transform vxc to cg_x as cg_x = P^(-1) * vxc
    do isp=1,nspin 
      call apply_trans(1,n1,n2,n3,nfft,qvec,px(:,isp),cg_Ap(:,isp))
    enddo
    ! regularization term 
    cg_Ap = cg_Ap + 1.d0/4.d0/pi/zp_coeff*cg_x

    cg_resid = cg_b - cg_Ap

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
  call print_cg_header() 


  !===========================
  !       cg loop 
  !===========================
  do while (.true.) 
    iter_cg = iter_cg +  1

    end_time = MPI_Wtime()
    call iter_info()
    start_time = MPI_Wtime()


    ! check convergence 
    d_resid = sum(cg_resid**2)*ucvol/nfft
    if (iter_cg==maxiter .or. (iter_cg>=10 .and. d_resid<1e-10)) then 
      if (nspin==1) then 
        write(logf,'(a,2es12.4)')'z: ',minval(cg_x(:,1)),maxval(cg_x(:,1))
      endif 
      if (nspin==2) then 
        write(logf,'(a,2es12.4,a)')'z:   ',minval(cg_x(:,1)),maxval(cg_x(:,1)),' spin up'
        write(logf,'(a,2es12.4,a)')'     ',minval(cg_x(:,2)),maxval(cg_x(:,2)),' spin down'
      endif 
      exit 
    endif 


    ! compute P*cg_p
    do isp=1,nspin 
      call apply_trans(1,n1,n2,n3,nfft,qvec,cg_p(:,isp),px(:,isp))
    enddo
    call calc_subsystem_dfpt(natom,natom_clu,natom_env,j_atom,nfft,nspin,ucvol, & 
                             vec1,vec2,qvec,sub_rhor,cell_nfft, & 
                             vec3,vks_clu,vks_env,px,cg_Ap)
    cg_Ap = -cg_Ap  ! since KS response function is negative definite 
    call remove_mean(nspin,nfft,cg_Ap)

    px = cg_Ap
    do isp=1,nspin 
      call apply_trans(1,n1,n2,n3,nfft,qvec,px(:,isp),cg_Ap(:,isp))
    enddo
    ! regularization term 
    cg_Ap = cg_Ap + 1.d0/4.d0/pi/zp_coeff*cg_p


    if (sum(cg_p*cg_Ap)*dvol<0.d0) then 
      write(715,*)'WARNINIG: pAp: ',sum(cg_p*cg_Ap)*dvol,' <0 '
      write(715,*)'d_resid: ',d_resid
      write(715,*)'WARNING: abs(pAp) < 1e-4 or d_resid<1e-5, move on'
    endif 
  

    if (do_precond) then 
      cg_alpha = sum(cg_resid*cg_z)/sum(cg_p*cg_Ap)
    else   
      cg_alpha = sum(cg_resid**2)/sum(cg_p*cg_Ap)
    endif 

    cg_x = cg_x + cg_alpha*cg_p

    ! back up old resids 
    cg_resid_old = cg_resid 
    cg_zold = cg_z

    ! new resid
    cg_resid  = cg_resid - cg_alpha*cg_Ap

    if (do_precond) then 
      do isp=1,nspin
        call precond_teter_inv(n1,n2,n3,nfft,qvec,1.d0,cg_resid(:,isp),cg_z(:,isp))
      enddo 
      cg_beta = sum(cg_z*cg_resid)/sum(cg_zold*cg_resid_old)
      cg_p = cg_z + cg_beta*cg_p
    else 
      cg_beta = sum(cg_resid**2)/sum(cg_resid_old**2) !! Fletcher–Reeves
      cg_p = cg_resid + cg_beta*cg_p
    endif 

  enddo 
  ! end of cg ============================


  ! apply transform matrix P 
  write(logf,*)'transform cg_x back to zvec.'
  do isp=1,nspin 
    call apply_trans(1,n1,n2,n3,nfft,qvec,cg_x(:,isp),zvec(:,isp))
  enddo 
  write(logf,'(a,2es12.4)')'zvec: ',minval(zvec,1),maxval(zvec,1)
  if (nspin==2)write(logf,'(a,2es12.4,a)')'zvec: ',minval(zvec,2),maxval(zvec,2),' spin down '
 

 ! ! penalty term 
 ! do s=1,nspin
 !   call laplacian(nfft,n1,n2,n3,zvec(:,s),qvec,lap)
 !   pen(s) = zp_alpha**2/2.d0/fourpi/zp_coeff*sum(zvec(:,s)**2)*dvol  & 
 !           -sum(zvec(:,s)*lap)*dvol/2.d0/fourpi/zp_coeff
 ! enddo
 ! ! compute |Ax-b| without regularization 
 ! ! to see how accurate the solution is 
 ! call calc_subsystem_dfpt(natom,natom_clu,natom_env,j_atom,nfft,nspin,ucvol, & 
 !                          vec1,vec2,qvec,sub_rhor,cell_nfft, & 
 !                          vec3,vks_clu,vks_env,zvec,cg_Ap)
 ! cg_Ap = -cg_Ap  ! since KS response function is negative definite 
 ! call remove_mean(nspin,nfft,cg_Ap)
 ! do s=1,nspin
 !   accu(s) = sum((cg_Ap(:,s)-cg_b(:,s))**2)*dvol
 ! enddo
 !
 !
  if (j_atom==1) then 
 !   write(logf,*)'solve_z => penalty term: ', pen
 !   write(logf,*)'solve-z => |Ax-b| (without reg. term): ',accu
    write(logf,'(a)')'left solve_z()'//NEW_LINE('a')
    call flush(logf)
  endif

 ! write(715,*)'solve_z => penalty: ',pen
 ! write(715,*)'solve_z => |Ax-b| (without reg. term): ',accu
  write(715,'(a)')'left solve_z()'//NEW_LINE('a')//NEW_LINE('a')
  write(715,'(a)') '==== done solve_z() ===='
  call flush(715)





!~~~~~~~~~~~~~~~~~~~~~~~~~~
! SUPPORTING FUNCTIONS
!~~~~~~~~~~~~~~~~~~~~~~~~~~
contains


  subroutine print_information()
    if (myrank==0) then 
      write(logf,*)''
      write(logf,'(a)')'   >>> enter solve_z_trans() <<<'
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



  subroutine iter_info 
   d_resid = sum(cg_resid**2)*ucvol/nfft
   write(715,'(i3,a,es8.2,es10.2,es10.2,2es10.2,es12.2,f8.1)') & 
     iter_cg-1,'   ',d_resid,cg_alpha,cg_beta, & 
     minval(cg_x(:,1)),maxval(cg_x(:,1)),sum(cg_x(:,1))/nfft,end_time-start_time
   call flush(715)
   if (myrank==0) then 
      write(logf,'(i3,a,es8.2,es10.2,es10.2,2es10.2,es12.2,f8.1)') & 
        iter_cg-1,'   ',d_resid,cg_alpha,cg_beta, & 
        minval(cg_x(:,1)),maxval(cg_x(:,1)),sum(cg_x(:,1))/nfft,end_time-start_time
      call flush(logf)
   endif 
  endsubroutine iter_info 


  !===================================================
  ! compute Px where P is the transformation matrix 
  !
  ! P = inv[sqrt( alpha**2 - laplacian)]
  ! 
  subroutine apply_trans(dir,n1,n2,n3,nfft,qvec,x,px)
    implicit none 

    integer :: n1,n2,n3,nfft,dir
    real(8) :: x(nfft), px(nfft), & 
               qvec(3,(n1/2+1)*n2*n3)
    
    ! local vars 
    integer :: dim1, ix, iy, iz
    real(8) :: q3d(3,n1/2+1,n2,n3), & 
               rtmp(n1,n2,n3),qq, factor 

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
          
          factor = sqrt(zp_alpha**2 + qq)

          if (dir == 1)  fft2(ix,iy,iz) = fft1(ix,iy,iz)/factor    ! P*x
          if (dir == -1) fft2(ix,iy,iz) = fft1(ix,iy,iz)*factor    ! P^(-1)*x
        enddo
      enddo
    enddo

    call FFT(n1,n2,n3,rtmp,fft2,-1)
    px = reshape(rtmp,(/n1*n2*n3/))

  endsubroutine apply_trans









end subroutine solve_z_trans







