!
!  Embedded clustetr density approximation (ECDA) 
!
!  Department of Scientific Computing 
!  Florida State University 
!  
!  Chen Huang  (Jan/2011)        
!
PROGRAM ECDA

 use mpi
 use interface_funcs
 use comm_data

 implicit none

 character(len=500) :: message,str,scmd
 character(len=100) :: filename,fname,ss,fname2,dtime
 character(len=pspfile_len),dimension(max_psp_file) :: psps_file
 integer,parameter  :: dp=8, &
                       mxband = 200, & ! the maximum number of band in subsystem 
                       logf = 9944, & 
                       mlist=100, & 
                       melement=200, & 
                       job_type=300   ! job_type=100: generate embedding potential
                                      ! job_type=200: self-consistent dfet
                                      ! job_type=300: xc_patch program

 logical            :: in_list, only_lda, & 
                       do_force = .false., & 
                       direct_min = .false. , & 
                       solve_vxc = .false., & 
                       file_ext, & 
                       plot_sub_rho


 integer, dimension(3) :: tarray
 integer            :: reg_type, & 
                       stat, & 
                       scf_mix_scheme = 1, & ! 1: just mix vxc
                                             ! 2: mix uemb as well
                                             ! 3: mix subsystem electron number as well
                       dft_flavor = 2, & 
                       npulay_gvks, &       ! historical global vks to mix
                       npulay_uemb = 2, &   ! historical global vks to mix
                       geo_method = 1, &  ! 1: ABINIT's BFGS
                                          ! 2: my L-BFGS 
                       exx_method = -1, & ! -1: real space code (poisson solver)
                                          !  0: set v(q=0)=0
                                          !  1: HSE (error function 
                                          !  2: Spencer Avali method 
                       iter_geo, & 
                       restart = -1, & 
                       relax_nscf, &        
                       relax_geo = -1


 real(8)            :: de, fd_step, step=0.d0, tsmear, & 
                       hse_range, & 
                       cov_radius( melement ), &
                       band_energy, & 
                       ehart_tmp, & 
                       zp_alpha, zp_conv, zp_coeff, & 
                       reg_vxc, patch_err, &
                       vac, & ! total system's fermi + vac defines the vaccum level
                       cgW_tol, dtmp, dtmp1, dtmp2

 integer            :: k, kk, rss, nspin, ia, & 
                       natom, npsps, natom_clu, natom_env, & 
                       do_xcep, do_exx=1, do_global_exx, & 
                       file_len,nfft,jj,nsys,i,ii,is, & 
                       date_time(8), n1,n2,n3, & 
                       zvec_init, gvxc_init, & 
                       global_fermi, system_type, & 
                       nopt_uemb, mix_gvks, & 
                       j,s,isp,start_u, & 
                       iter_scf,time_start, &
                       cell_nfft(3)

 ! Buffer atoms 
 integer,allocatable :: nlist(:),  &   ! number of buffer atoms 
                        atom_class(:), & ! 1: cluster atom, 2: env atom 
                        ilist(:,:), &  ! indeces of buffer atoms
                        nshare(:), &   ! the number of boundary atoms to be shared 
                        ishare(:,:)    ! the boundary atoms to be shared between cluster and env.

 ! mpi related 
 integer :: ierr
 integer :: p_num
 integer :: myrank 

 real(kind=dp),allocatable :: & 
   watom(:,:), &      ! atom weight function. sum(watom)=1.0 all over the space
   atom_info(:,:), & 
   lap(:), & 
   grad(:,:), & 
   exc_eps(:), & 
   vhart(:), &
   vtmp(:,:), &
   zvec(:,:), & 
   zvec_old(:,:), & 
   pvec(:,:), & 
   points(:,:), &      ! r points in the cell 
   global_vks(:,:), & 
   global_vOF(:,:), & 
   global_vks0(:,:), & 
   global_vks_new(:,:), & 
   global_vpsp(:), &
   global_vxc(:,:), & 
   global_vxc_old(:,:), & 
   dExc_dvks(:,:), & 
   dE_dvks(:,:), & 
   dE_dvks0(:,:), & 
   rho_perb(:,:), & 
   patch_rho(:,:), & 
   metric(:), & 
   dfermi_dvks(:,:), & 
   gvec(:,:), &  
   ref_rho(:,:),  & 
   rho_old(:,:),  & 
   rho_tmp(:,:),  & 
   arr_tmp(:),  & 
   rho_tmp2(:,:),  & 
   vperb(:,:),  & 
   ke_pot(:),  & 
   total_dExc_dvks(:,:), & 
   total_fermi(:), & 
   qvec(:,:), &        ! q-vector in q-space
   u(:,:), den_diff_old(:,:), & 
   den_diff(:,:), & 
   u_old(:,:), &           ! u(r) potential to be find
   pulay_work_gvks(:,:,:), & 
   pulay_work_uemb(:,:,:), & 
   xcart(:,:), & 
   xcart_tmp(:,:), & 
   xcart_clu(:,:), & 
   xcart_env(:,:), & 
   znucl(:), & 
   zion(:), & 
   force_pen(:,:), & 
   force_abinit(:,:), &    ! force due to ion-ion and ion-electron 
   force_total(:,:), & 
   f_noavg(:,:)


 real(kind=dp) :: & 
   dipole(3), dmu, vwreg_phi, & 
   xpbe, cpbe, & 
   test_fft(3,2,2), & 
   resid_gvks_dtmp, & 
   x1,x2, & 
   pulay_gvks_beta, & 
   cgW, int_exc, & 
   exc_exx, exc_exx_tmp,  & 
   total_ehart, & 
   total_ehart_tmp, & 
   total_ke,  & 
   total_ke_tmp,  & 
   total_ewald, total_corepsp, & 
   total_nlpsp, total_nlpsp_tmp, & 
   total_TS, total_TS_tmp, & 
   total_elocal, total_elocal_tmp, &
   etotal,& 
   etotal_old_geo, & 
   etotal_old=0.d0, & 
   etotal_old2=0.d0, & 
   reg_energy, & 
   cell_acell(3,3), &  ! the a b c of cell
   gmet(3,3),  &       ! 
   abinit_gmet(3,3),  &   ! differ by gmet by a factor of (2*pi)**2
   rmet(3,3),  &     ! 
   gprim(3,3) , &    ! prime vectors in q-space
   rprimd(3,3), &    ! lattice vector in bohr
   ucvol , &         ! cell volume
   dvol              ! volume / nfft


 ! debug my L-BFGS dir code 
 integer,parameter :: db_mhist=3
 real(8) :: db_ghist(5,db_mhist+1), & 
            db_xhist(5,db_mhist+1), & 
            db_lbfgs_z(5,1), &  
            Amat_tmp(5,5), & 
            bvec_tmp(5,1), & 
            xvec_tmp(5,1) , & 
            gvec_tmp(5,1) 


 ! timing 
 real(8) :: start_time, end_time 
 real(8) :: dotprod_four


 ! My L-BFGS 
 integer :: lbfgs_mhist = 5
 real(8),allocatable :: lbfgs_ghist(:,:), & 
                        lbfgs_xhist(:,:), & 
                        lbfgs_z(:)

 ! L-BFGS for geometry relaxation 
 integer             :: geo_mhist = 5
 real(8),allocatable :: geo_ghist(:,:), & 
                        geo_xhist(:,:), & 
                        geo_z(:)

 ! L-BFGS from Nocedal and Zhu
 integer             :: mhist = 5
 logical             :: lsave(4)
 character*60        :: task, csave
 integer             :: isave(44)
 real(8)             :: dsave(29), gmax, gnorm
 integer,allocatable :: lbfgs_nbd(:), lbfgs_iwa(:)
 real(8),allocatable :: lbfgs_l(:), lbfgs_u(:), lbfgs_wa(:)
 real                :: lbfgs_time 

 !=========================================================

 call MPI_Init ( ierr )
 call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, ierr )
 call MPI_Comm_size ( MPI_COMM_WORLD, p_num, ierr ) 
 call clean_folder()

 if (myrank==0) then 
   call system("rm *.xsf")
   call system("rm xcep.out log_cluster_*.txt iterate_*.dat  cluster_*.vasp")
   open(file='xcep.out',action='write',unit=logf,form='formatted')
   call welcome_message(logf)
 endif
 write(logf,'(a,i4)') 'number of processors: ',p_num; call flush(logf)
 call mpi_barrier(mpi_comm_world, ierr)


 ! OF-DFT log files
 write(fname,*)myrank+1
 fname = 'log_OF_solver_'//trim(adjustl(fname))//'.txt'
 open(file=fname,action='write',unit=logOF,form='formatted')

 zp_alpha = 1/5.d0  ! 5 bohr as screening length
 call load_parameters(restart,logf,natom,nspin,npsps,cell_nfft,nsys,global_fermi,mix_gvks,& 
                      start_u,dft_flavor,nopt_uemb,npulay_gvks,pulay_gvks_beta,cgW_tol, & 
                      cell_acell,ucvol,rprimd,rmet,tsmear,relax_geo,only_lda,& 
                      reg_type,zp_alpha,reg_vxc,zp_coeff,zp_conv,vac,plot_sub_rho,relax_nscf, & 
                      exx_method,hse_range)

 if (p_num /= nsys) then 
   print *,"ERROR: number of processors is not set to the number of atoms"
   stop
 endif

 n1 = cell_nfft(1)
 n2 = cell_nfft(2)
 n3 = cell_nfft(3)
 nfft = n1*n2*n3
 dvol = ucvol/dble(nfft)

 allocate(xcart(3,natom), & 
          xcart_tmp(3,natom), & 
          xcart_clu(3,natom), & 
          xcart_env(3,natom), & 
          znucl(natom),zion(natom), & 
          u(nfft,nspin), &         ! embedding potential 
          u_old(nfft,nspin), &      ! embedding potential 
          den_diff(nfft,nspin), &      ! embedding potential 
          den_diff_old(nfft,nspin), &      ! embedding potential 
          atom_info(5,natom), &    ! (subsystem_index, type, x,y,z) x,y,z are in angstrom 
          points(3,nfft), & 
          global_vks(nfft,nspin), & 
          global_vOF(nfft,nspin), & 
          global_vks0(nfft,nspin), & 
          global_vpsp(nfft), &
          metric(nfft), & 
          force_pen(3,natom), & 
          force_abinit(3,natom), & 
          force_total(3,natom), & 
          f_noavg(3,natom), & 
          global_vxc(nfft,nspin), & 
          global_vxc_old(nfft,nspin), & 
          dExc_dvks(nfft,nspin), & 
          dE_dvks(nfft,nspin), & 
          dE_dvks0(nfft,nspin), & 
          vperb(nfft,nspin), & 
          rho_perb(nfft,nspin), & 
          patch_rho(nfft,nspin), & 
          dfermi_dvks(nfft,nspin), & 
          gvec(nfft,nspin), & 
          zvec(nfft,nspin), & 
          zvec_old(nfft,nspin), & 
          pvec(nfft,nspin), & 
          exc_eps(nfft), & 
          lap(nfft), & 
          grad(3,nfft), & 
          total_dExc_dvks(nfft,nspin),& 
          total_fermi(nspin),& 
          global_vks_new(nfft,nspin), & 
          pulay_work_gvks(nfft*nspin,npulay_gvks,2), & ! pulay array for global vks 
          pulay_work_uemb(nfft*nspin,npulay_gvks,2), & ! pulay array for global vks 
          qvec(3,(n1/2+1)*n2*n3), & 
          vhart(nfft), & 
          arr_tmp(nfft), & 
          vtmp(nfft,nspin), & 
          nlist(natom), & 
          atom_class(natom), & 
          ilist(mlist,natom),& 
          nshare(natom), & 
          watom(nfft,nsys), & 
          ishare(mlist,natom), & 
          ref_rho(nfft,nspin), & 
          rho_old(nfft,nspin), & 
          rho_tmp(nfft,nspin), & 
          rho_tmp2(nfft,nspin), & 
          ke_pot(nfft), & 
          lbfgs_ghist(nfft*nspin,lbfgs_mhist+1),  &   ! L-BFGS 
          lbfgs_xhist(nfft*nspin,lbfgs_mhist+1), &    ! L-BFGS 
          lbfgs_z(nfft*nspin), & ! my l-bfgs
          geo_ghist(3*natom,geo_mhist+1),  &   ! L-BFGS 
          geo_xhist(3*natom,geo_mhist+1), &    ! L-BFGS 
          geo_z(3*natom), &                      ! my l-bfgs
          lbfgs_nbd(nfft*nspin), & 
          lbfgs_u(nfft*nspin), & 
          lbfgs_l(nfft*nspin), & 
          lbfgs_iwa(3*nfft*nspin), & 
          lbfgs_wa(2*mhist*nfft*nspin+4*nfft*nspin+12*mhist*mhist+12*mhist))
 
 
 call load_atom_info(natom,xcart,znucl,atom_info,mlist,nlist,ilist,nshare,ishare)
 call load_psp_info(npsps,psps_file)
 call make_gprim(cell_acell,gprim,gmet,abinit_gmet)
 call make_q_vector(cell_nfft(1),cell_nfft(2),cell_nfft(3),gprim,qvec)
 call make_points(transpose(rprimd),nfft,n1,n2,n3,points)
 call init_conv_radius(melement, cov_radius)

 natom_clu = nlist(myrank+1)+1
 natom_env = natom-natom_clu

 zvec_init = -1
 gvxc_init = -1

 if (global_fermi<=0) then 
   print *,"it is highly recommended to set to global_fermi==1"
   stop
 endif 

 select case (dft_flavor)
 case (-101,-102) 
   do_xcep = 1
   do_exx  = 1
   do_global_exx = -1
 case (11,12)
   do_xcep = -1
   do_exx  = 1
   do_global_exx = 1
 case (1,2) 
   do_exx = -1
   do_xcep = -1 
   do_global_exx = -1
   reg_vxc = 0.d0  ! no penalty functional, just LDA and PBE 
 case default 
   print *,'undefined dft_flavor'
   stop 
 endselect 

 ! each cluster has its log file 
 write(message,'(i4)')myrank+1
 write(message,'(a)')'log_cluster_'//trim(adjustl(message))//'.txt'
 open(file=adjustl(trim(message)),unit=715,action='write',form='formatted')


 if (myrank==0) write(logf,'(a)')'making atom weights...'; call flush(logf)
 if (do_pbc) then 
   call make_weights_pbc(logf,natom,nspin,nfft,nsys,melement,atom_info,& 
                     cell_nfft,xcart,cell_acell,znucl,rprimd,cov_radius,watom)
 else 
   call make_weights(logf,natom,nspin,nfft,nsys,melement,atom_info,& 
                     cell_nfft,xcart,cell_acell,znucl,rprimd,cov_radius,watom)
 endif 

 call init_uemb(nspin,nfft,start_u,u)

 ! load initial guess of gloabl vks 
 if (restart<=0) then 
   if (myrank==0) then 
     call check_file_length(nspin,n1,n2,n3)
   endif 
   call mpi_barrier(mpi_comm_world, ierr)
   call load_init_guess_vks(nfft,nspin,logf,global_vks)
   write(logf,'(a)')NEW_LINE('a')//'Load initial vks from init_guess_total_vks.dat...'
   write(logf,'(a)')'removed the mean of vks'
 else 
   ! restat, loading vks
   if (myrank==0) then 
     write(logf,*)
     write(logf,'(a)')'*** restart XCEP based on previous global_vks (in vks_restart.dat) ***'
     write(logf,'(a)')'removed the mean of vks'
     write(logf,*)
   endif 
   open(file='vks_restart.dat', action='read',form='unformatted',unit=111)
   read(111)global_vks
   close(111)
 endif 
 if (myrank==0) then 
   write(logf,'(a,2es12.4,a)')'min/max init_vks: ', & 
     minval(global_vks(:,1)),maxval(global_vks(:,1)),' spin up'
   if (nspin==2) & 
     write(logf,'(a,2es12.4,a)')'min/max init_vks: ', & 
       minval(global_vks(:,2)),maxval(global_vks(:,2)),' spin down'
   write(logf,'(a)')''
 endif 


 ! display time and memory needs
 call mpi_barrier(mpi_comm_world, ierr)
 call itime(tarray)
 time_start = tarray(1)*60.d0*60.d0 + tarray(2)*60.d0 + tarray(3)
 call system_mem_usage(rss)
 write(715,'(a,f12.4,a,f12.4,a)')'memory used (RSS): ',rss/1000.0, ' MByte', rss/1.0e6,' GB'

 if (myrank==0)  & 
   call write_cluster_env_structure_XSF(natom,rprimd,xcart,ilist,nlist,mlist,nshare,ishare)

 ! load cluster and env atom coordinates from ABINIT subsystems 
 if (do_global_exx<0)  & 
   call load_xcart_cluster_env(myrank,natom,natom_clu,natom_env,xcart_clu,xcart_env)

 ! load vpsp from global system 
 do while (.true.)
   open(file='./global_system/vpsp.dat',unit=111,action='read',form='unformatted',iostat=stat)
   if (stat/=0) then 
     call sleep(1)
     cycle 
   endif 
   read(111,iostat=stat) global_vpsp
   if (stat/=0) then  
     close(111); call sleep(1)
     cycle 
   else
     close(111)
     exit
   endif 
 enddo 


 ! Reset the average of global_vks from LDA calculations. 
 ! By using the regularized total energy functional, it can be shown that 
 ! the integration of vxc in space is zero. Therefore, the integration of vks 
 ! in space is equal to that of vpsp(q=0), due to the fact that vhart(q=0)=0
 ! (note that vhart is computed with FFT by setting q=0 to zero)
 do isp=1,nspin
   global_vks(:,isp) = global_vks(:,isp) & 
     - sum(global_vks(:,isp))/nfft + sum(global_vpsp)/nfft
 enddo


 ! write system information file 
 open(file='system_info.dat',action='write',unit=111,form='unformatted')
 write(111) exx_method 
 write(111) hse_range
 close(111)
 call flush(6)



 !============================================================
 ! SCF loop for optimizing system's gloable KS potential 
 !============================================================
 iter_scf = 0
 iter_geo = 1

 do while (.true.)

   iter_scf = iter_scf + 1
   if (iter_scf==1) task='START'
   if (myrank==0 .and. relax_geo==1) & 
     write(logf,'(a,i4,a)')'   ==> geometry relaxation step: ',iter_geo,' <=='
   if (myrank==0) call print_header_optimize_global_system(logf,iter_scf)


   !============================
   ! run global KS with new vks
   !============================
   if (myrank==0) then 
     write(logf,'(a,2f12.4)')'vks(spin up): ',minval(global_vks(:,1)),maxval(global_vks(:,1))
     if (nspin==2)write(logf,'(a,2f12.4)')'vks(spin dn): ',minval(global_vks(:,2)),maxval(global_vks(:,2))
     write(logf,'(a)')'calling run_global_ks() for new global density ...'; call flush(logf)

     call run_global_ks(n1,n2,n3,nfft,nspin,natom,qvec,ucvol,do_global_exx,do_xcep, & 
                        global_vks,ref_rho,total_ehart,total_ke,total_nlpsp,total_TS, &
                        total_elocal,total_ewald,total_corepsp,total_fermi,dfermi_dvks, & 
                        band_energy,exc_eps,exc_exx,total_dExc_dvks,force_abinit,vhart)

     call read_system_zion(natom,zion)
     call dipole_moment(logf,natom,nspin,n1,n2,n3,dvol,points,ref_rho,xcart,zion,dipole)
     call print_info()

     ! PBE0? 
     if (dft_flavor==12) then 
       write(logf,*)'PBE0 ==> modify total_dExc_dvks and PBE0 energy'
       ! compute dExc/dvks for the PBE part
       call calculate_x_pbe(n1,n2,n3,nspin,ref_rho,qvec,dvol,vtmp,xpbe)  ! PBE exchange 
       vperb = 0.75d0*vtmp
       call calculate_c_pbe(n1,n2,n3,nspin,ref_rho,qvec,dvol,vtmp,cpbe)  ! PBE correlation 
       vperb = vperb + vtmp
       ! apply chi_total 
       system_type = 0 ! 0: global system
       call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,vperb,rho_perb)
       call remove_mean(nspin,nfft,rho_perb)
       ! get dExc/dvks for PBE0
       total_dExc_dvks = total_dExc_dvks*0.25d0 + rho_perb
       int_exc = exc_exx*0.25d0 + xpbe*0.75d0 + cpbe  ! PBE0 energy
     endif 
   endif

   ! send results from global_system to all processors for dfet
   call send_data_master_to_all(nspin,myrank,total_fermi)
   call send_data_master_to_all(1,myrank,total_ehart)
   call send_data_master_to_all(3*natom,myrank,force_abinit)
   call send_data_master_to_all(nspin*nfft,myrank,dfermi_dvks)
   call send_data_master_to_all(nspin*nfft,myrank,ref_rho)
   call send_data_master_to_all(nfft,      myrank,exc_eps)
   call send_data_master_to_all(nspin*nfft,myrank,total_dExc_dvks)
   call send_data_master_to_all(nfft,myrank,vhart)


   !==========================================
   ! compute XC potential for global EXX case
   !==========================================
   if (do_global_exx>0) then 
     if (myrank==0) then 
       ! set avg of total_dExc_dvks to zero 
       call print_dExc_dvs(nspin,nfft,total_dExc_dvks)
       write(logf,'(a)')'set mean of dExc_dvxc to zero'
       do isp=1,nspin 
         total_dExc_dvks(:,isp) = total_dExc_dvks(:,isp) - sum(total_dExc_dvks(:,isp))/dble(nfft)
       enddo

       !!call solve_total_oep(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
       !!  logf,ucvol,global_vks,reg_vxc,gvxc_init,total_dExc_dvks,global_vxc)
       call solve_total_oep_regvxc(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
           logf,ucvol,global_vks,vhart,global_vpsp,reg_vxc,gvxc_init,total_dExc_dvks,global_vxc)
       !!The trans method does not work well in practice, vxc convergence is slower
       !call solve_total_oep_regvxc_trans(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
       !    logf,ucvol,global_vks,vhart,global_vpsp,reg_vxc,gvxc_init,total_dExc_dvks,global_vxc)

       write(logf,'(a,2f12.4)')'global_vxc (up): ',minval(global_vxc(:,1)),maxval(global_vxc(:,1))
       if (nspin==2) & 
       write(logf,'(a,2f12.4)')'global_vxc (dn): ',minval(global_vxc(:,2)),maxval(global_vxc(:,2))
     endif
     call send_data_master_to_all(nspin*nfft,myrank,global_vxc)
   endif 


   if (myrank==0) then 
     call force_penalty(nspin,nfft,n1,n2,n3,natom,qvec,vhart, & 
                        global_vks,global_vpsp,reg_vxc,force_pen)
   endif  
   call send_data_master_to_all(natom*3,myrank,force_pen)


   !=========================
   ! load atom coordinates 
   !=========================
   if (iter_scf==1) then 
     ! load new atom coordinates from global system (change during geometry relaxation)
     open(file='./global_system/atom_coords.dat',unit=111,action='read',form='unformatted')
     read(111) xcart
     close(111)
     ! get the atom index for cluster and env atoms 
     ! only do this at the very first iteration 
     if (iter_scf==1 .and. iter_geo==1) then 
       do i=1,natom 
          atom_class(i)=2  ! env atom 
          ! check if the atom is in cluster?
          do ia=1,natom_clu 
            if (sqrt(sum(xcart(:,i)-xcart_clu(:,ia))**2)<1e-3) then 
              atom_class(i)=1  ! cluster atom 
            endif 
          enddo 
       enddo 
     endif 

     call make_weights(logf,natom,nspin,nfft,nsys,melement,atom_info,& 
                       cell_nfft,xcart,cell_acell,znucl,rprimd,cov_radius,watom)

     if (myrank==0) then 
       write(logf,'(a)')NEW_LINE('a')//'Atoms from global system are loaded and atom weights are updated'
       write(logf,'(a)')'new coordinates (angstrom):'
       do ia=1,natom
         write(logf,'(a,i4,a,3f14.4)')"atom ",ia," ",xcart(:,ia)*0.52917721067d0
       enddo
       write(logf,*)''
       call flush(logf)
       filename = 'geo.xsf'
       call wrt_xsf_geo(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename)
     endif 

     if (do_global_exx<0) then 
       ! write new cluste and env coordiates to ABINIT subsystems
       ! xcart_clu and xcart_env will be updated as well 
       call update_coords_subsystem(logf,myrank,natom,natom_clu,natom_env, & 
                                    atom_class,xcart,xcart_clu,xcart_env)
     endif 
   endif


   !=======================================================
   ! compute vks_OF due to OFDFT, KEDF = vw_lam * vW + TF
   !=======================================================
   if (do_envOF) then 
     do s=1,nspin 
       !
       ! different regularization type for KEDF
       !
       select case(vwreg_type)
       case (1) 
         call vw(nfft,n1,n2,n3,ref_rho(:,s),qvec,dvol,vw_reg,ke_pot,dtmp)
       case (2)
         ! regularized rho 
         arr_tmp = ref_rho(:,s)/vw_reg
         rho_tmp(:,s) = (1.d0-exp(-arr_tmp))*ref_rho(:,s)
         call laplacian(nfft,n1,n2,n3,sqrt(rho_tmp(:,s)),qvec,lap)

         ! compute d sqrt(rho')/d rho
         vtmp(:,s) = (1.d0-exp(-arr_tmp) & 
                      + exp(-arr_tmp) * ref_rho(:,s)/vw_reg)/2.d0/sqrt(rho_tmp(:,s))

         ke_pot = -lap * vtmp(:,s)
       case (3)
         ! regularized rho 
         vwreg_phi = sqrt(vw_reg)
         arr_tmp = ref_rho(:,s)**0.5/vwreg_phi
         rho_tmp(:,s) = (1.d0-exp(-arr_tmp))*ref_rho(:,s)
         call laplacian(nfft,n1,n2,n3,sqrt(rho_tmp(:,s)),qvec,lap)

         ! compute d sqrt(rho')/d rho
         vtmp(:,s) = (1.d0 - exp(-arr_tmp) + exp(-arr_tmp)*sqrt(ref_rho(:,s))/2.d0/vwreg_phi)/&
                      2.d0/sqrt(rho_tmp(:,s))

         ke_pot = -lap * vtmp(:,s)
       endselect 

       global_vOF(:,s) = total_fermi(s) - ke_pot*vw_lam
       call tf(nfft,ref_rho(:,s),dvol,ke_pot,dtmp)
       global_vOF(:,s) = global_vOF(:,s) - ke_pot
       write(logf,'(a,2es12.4,a,i3,a,i2)') & 
         'global_vOF: ',minval(global_vOF(:,s)),maxval(global_vOF(:,s)), & 
         ' for spin:',s,'   vwreg_type: ',vwreg_type
     enddo 

     call send_data_master_to_all(nspin*nfft,myrank,global_vOF)

     if (myrank==0) then
       filename = 'global_vks_OF_spin_up.xsf'
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,global_vOF(:,1))
       filename = 'global_vks_OF_spin_dn.xsf'
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,global_vOF(:,2))
     endif
   endif 


   !!!!!!!!!!!!!!!! DEBUG OF-DFPT !!!!!!!!!!!!!!!
   !call debug_OF_dfpt(nspin)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
   !================================================================
   ! compute XC potential based on different XC functionals or ECDA
   !================================================================
   if (dft_flavor==11) then 
     ! exx
     write(logf,'(a)')'dft_flavor=100 => percent EXX (full OEP)'
     int_exc     = exc_exx
     force_total = force_abinit + force_pen
     if (myrank==0) then 
       write(logf,'(a)')NEW_LINE('a')//'total forces (hartree/bohr): '
       do ia=1,natom 
         write(logf,'(a,i4,3f14.6)')'  atom: ',ia,force_total(1:3,ia)
       enddo
     endif 

   else if (dft_flavor==2) then 
     write(logf,'(a)')'dft_flavor => PBE'
     call calculate_xc_pbe(n1,n2,n3,nspin,ref_rho,qvec,dvol,global_vxc,int_exc)
     force_total = force_abinit

   else if (dft_flavor==1) then 
     write(logf,'(a)')'dft_flavor => PW92LDA'
     call calculate_xc_pw92lsda(n1,n2,n3,nspin,ref_rho,dvol,global_vxc,int_exc)
     force_total = force_abinit

   else if (dft_flavor<0) then 
     ! patch_xc()
     !  input: ref_rho and global_vks 
     !  output: dExc_dvks and int_exc
     write(logf,'(a)')'dft_flavor => patched EXX, call patch_xc()'
     call patch_xc(natom,natom_clu,natom_env,nfft,nspin,mxband,myrank,cell_nfft, & 
                   atom_info,cell_acell,plot_sub_rho,logf,mlist,ilist,nlist,nshare,ishare,qvec, & 
                   solve_vxc,nopt_uemb,reg_vxc,zp_alpha,zp_coeff,zp_conv,relax_geo,rprimd,tsmear,ucvol, & 
                   dft_flavor,vac,xcart,znucl,global_fermi,cov_radius,melement,do_force,do_xcep, & 
                   dvol,iter_scf,total_fermi,watom,dfermi_dvks,gvxc_init,zvec_init,cgW_tol, & 
                   ref_rho,u,den_diff,do_exx,zvec,global_vks,global_vOF,patch_rho,int_exc, & 
                   force_pen,force_abinit,force_total,dExc_dvks,only_lda,exx_method,hse_range)
     write(logf,'(a,f18.8)')'[main] exc from xcpatch: ',int_exc
   endif 


   !============================================
   ! regularization energy (note the 1/2 factor) 
   !============================================
   if (.not. direct_min) then 
     do isp=1,nspin
       pvec(:,isp) = global_vks(:,isp) - global_vpsp - vhart
     enddo 
     select case(reg_type)
     case (1)
       reg_energy = 0.5d0*reg_vxc*sum(pvec**2)*dvol
     case (2)
       reg_energy = 0.d0 
       !<grad vks, grad vks>
       do isp=1,nspin 
         call laplacian(nfft,n1,n2,n3,pvec(:,isp),qvec,lap)
         reg_energy = reg_energy - 0.5d0*reg_vxc*sum(pvec(:,isp)*lap)*dvol
       enddo 
     case (3)
       reg_energy = 0.5d0*reg_vxc*sum(pvec**2)*dvol
       !<grad vks, grad vks>
       do isp=1,nspin 
         call laplacian(nfft,n1,n2,n3,pvec(:,isp),qvec,lap)
         reg_energy = reg_energy - 0.5d0*reg_vxc*sum(pvec(:,isp)*lap)*dvol
       enddo 
     end select
   else 
     reg_energy = 0.d0
   endif 


   !==================================
   ! solve oep equation to get new vxc
   !==================================
   if (.not. direct_min  .and. do_global_exx<0 .and. do_xcep>0) then 
     if (myrank==0) then 
       write(logf,'(a)')'solving for system KS potential ...'
       global_vxc_old = global_vxc

       ! set avg of dExc_dvks to zero 
       call print_dExc_dvs(nspin,nfft,dExc_dvks)
       write(logf,'(a)')'set mean of dExc_dvxc to zero'
       do isp=1,nspin 
         dExc_dvks(:,isp) = dExc_dvks(:,isp) - sum(dExc_dvks(:,isp))/dble(nfft)
       enddo

       call solve_total_oep_regvxc(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
         logf,ucvol,global_vks,vhart,global_vpsp,reg_vxc,gvxc_init,dExc_dvks,global_vxc)
       !!The trans method does not work well in practice, vxc convergence is slower
       !!   call solve_total_oep_regvxc_trans(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
       !!     logf,ucvol,global_vks,vhart,global_vpsp,reg_vxc,gvxc_init,dExc_dvks,global_vxc)
       !!
       !!call solve_total_oep(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
       !!  logf,ucvol,global_vks,reg_vxc,gvxc_init,dExc_dvks,global_vxc)

       write(715,'(a,2f14.6)')'new global vxc: ',minval(global_vxc),maxval(global_vxc)
     endif 
     call send_data_master_to_all(nspin*nfft,myrank,global_vxc)
   endif 


   ! write vxc files for VESTA
   if (.not. direct_min .and. myrank==0) then 
     if (nspin==2) then 
       if (iter_scf==1) filename = "global_vxc_spin_up_iter0.xsf"
       if (iter_scf>1)  filename = "global_vxc_spin_up.xsf"
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,global_vxc(:,1))
       if (iter_scf==1) filename = "global_vxc_spin_dn_iter0.xsf"
       if (iter_scf>1)  filename = "global_vxc_spin_dn.xsf"
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,global_vxc(:,2))
     else
       if (iter_scf==1) filename = "global_vxc_iter0.xsf"
       if (iter_scf>1)  filename = "global_vxc.xsf"
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,global_vxc(:,1))
     endif 
   endif

   ! correction to XC since patch rho is not exactly ref_rho
   patch_err = sum((ref_rho-patch_rho)*global_vxc)*dvol

   ! get total energy (direct sum)
   etotal =  total_ehart  + total_ke      + total_nlpsp & 
           + total_ewald  + total_corepsp + total_TS & 
           + total_elocal + int_exc       + reg_energy 


   if (iter_scf==1) etotal_old = etotal 
   de = etotal - etotal_old
   etotal_old2 = etotal_old
   etotal_old  = etotal


   !=========================
   ! pulay mixing of the vks
   !=========================
   if (.not. direct_min) then 
     ! new global vks 
     do isp=1,nspin
       global_vks_new(:,isp) = global_vxc(:,isp) + vhart + global_vpsp
     enddo 
     if (iter_scf==1) resid_gvks_dtmp = 0.d0
     if (iter_scf>=2) resid_gvks_dtmp = sqrt(sum((global_vks_new - global_vks)**2))

     ! mixing 
     select case (scf_mix_scheme) 
     case (1) 
       if (precond_vks) then 
          write(logf,*)'pulay mixing of system vks (kerker).'  
          call pulay_mix_kerker(n1,n2,n3,nfft*nspin,iter_scf,global_vks,global_vks_new, & 
                               pulay_work_gvks,npulay_gvks,pulay_gvks_beta,qvec,ucvol)
       else 
          write(logf,*)'pulay mixing of system vks.'  
          call pulay_mix(nfft*nspin,iter_scf,global_vks,global_vks_new, & 
                         pulay_work_gvks,npulay_gvks,pulay_gvks_beta)
       endif 
       global_vks = global_vks_new ! update system vks 
     case (2) 
       ! mix vks and uemb together 
       ! the pulay mixing has to be run in parallel to save memory
       write(logf,'(a)')'pulay mixing of [vks,uemb].'
       call pulay_mix_vks_uemb(nfft*nspin,iter_scf,global_vks,global_vks_new, & 
                               u_old,u,pulay_work_gvks,pulay_work_uemb, & 
                               npulay_gvks,pulay_gvks_beta)
       global_vks = global_vks_new ! update system vks 
     end select
   endif 

   write(715,'(a,2f12.6)')'next global_vks (spin_up): ',minval(global_vks(:,1)),maxval(global_vks(:,1))
   if (nspin==2) & 
   write(715,'(a,2f12.6)')'next global_vks (spin_dn): ',minval(global_vks(:,2)),maxval(global_vks(:,2))

   u_old = u ! backup embedding potential 


   !=====================
   ! direct minimization 
   !=====================
   if (direct_min .and. do_global_exx<0) then 
     ! make dEtot/dvks 
     if (myrank==0) then 
       write (logf,'(a)')NEW_LINE('a')//'direct update global_vks using L-BFGS'
       do isp=1,nspin
         vtmp(:,isp) = vhart + global_vpsp - global_vks(:,isp)
       enddo
       system_type = 0 ! 0: global system
       write(logf,'(a)')'computing d(J+ext+Ts)/dvks ...'
       call flush(logf)
       call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,vtmp,rho_perb)
     endif 
     call send_data_master_to_all(nspin*nfft,myrank,rho_perb)

     dE_dVks = dExc_dVks + rho_perb
     call display_dE_dVks()

     gmax  = maxval(abs(dE_dVks))
     gnorm = sqrt(sum(dE_dVks**2)*dvol)

     call lbfgs_dir(nspin,nfft,lbfgs_mhist,iter_scf,dE_dvks,global_vks, & 
                    lbfgs_ghist,lbfgs_xhist,lbfgs_z)

     global_vks = global_vks + reshape(lbfgs_z,(/nfft,nspin/))*0.5d0
     call send_data_master_to_all(nspin*nfft,myrank,global_vks)

!!!!!!!!!!!!!!!!!! DEBUG !!!!!!!!!!!!!!!
!     if (iter_scf==1) then 
!       dE_dvks0    = dE_dvks
!       global_vks0 = global_vks
!     endif
!     if (myrank==0) then 
!       write(logf,*)'etotal/dE_dstep/step: ',etotal,sum(-dE_dvks0*dE_dVks)*dvol,step
!     endif 
!     step = step + 0.02d0
!     global_vks = global_vks0 - dE_dVks0*step
!!!!!!!!!!!!!!!!!! DEBUG !!!!!!!!!!!!!!!

!9999 continue 
!     ! L-BFGS on master processor 
!     if (myrank==0) then 
!       lbfgs_l = 0.d0 
!       lbfgs_u = 0.d0 
!       lbfgs_nbd = 0    ! unbounded 
!       call setulb(nfft*nspin, mhist, global_vks, lbfgs_l, lbfgs_u, lbfgs_nbd,  & 
!                   etotal, dE_dvks, 1.d1, 0.d0, lbfgs_wa, lbfgs_iwa, & 
!                   task, 99, csave, lsave, isave, dsave, lbfgs_time, myrank)
!
!       if (task(1:2) == 'FG') then 
!         write(logf,'(a)')'L-BFGS returned FG'
!       elseif (task(1:5) == 'NEW_X') then 
!         write(logf,'(a)')'L-BFGS returned NEW_X'
!         goto 9999
!       else
!         write(logf,*)'code stopped due to unknow message from L-BFGS'
!         write(logf,*)'task: ',trim(task)
!       endif 
!     endif 
!     call send_data_master_to_all(nspin*nfft,myrank,global_vks)
   endif 
   ! end of direct minimizaiton ==================


   if (myrank==0) call print_energy_info()
   

   !===================================================
   ! geometry relaxation, send forces to global system 
   ! and global system will move the atoms
   !===================================================
   if ( relax_geo==1 .and. iter_scf==relax_nscf )  then 
     if (myrank==0) then 
       if (iter_geo>1) then 
         write(logf,'(a,i4,a,f14.8,a,a,es10.2,a,es10.2)')"@ (geo. relaxation) iter_geo: ", & 
           iter_geo,"  etotal= ",etotal,'  Ha',' dE: ',etotal - etotal_old_geo, &
           ' fmax: ',maxval(abs(force_total))
       else
         write(logf,'(a,i4,a,f14.8,a,a,es10.2)')"@ (geo. relaxation) iter_geo: ", & 
           iter_geo,"  etotal= ",etotal,'  Ha', ' fmax: ',maxval(abs(force_total))
       endif  
       write(logf,'(a)')"finished ECDA for this geometry"
       write(logf,*)""

       if (geo_method==1) then ! ABINIT's BFGS 
         write(logf,'(a)')NEW_LINE('a')//'Total forces: (Ha/bohr)'
         do ia=1,natom 
           write(logf,'(a,i4,3f14.6)')'atom: ',ia,force_total(1:3,ia)
         enddo 
         write(logf,'(a)')NEW_LINE('a')
         ! send force to global_KS
         open(file='global_system/xcep_force.dat',unit=111,action='write',form='unformatted')
         write(111) force_total
         close(111)
         write(logf,'(a)')'Use ABINIT BFGS optimizer. Sending forces to global ABINIT ...'
         ! tell global system XCEP is finished
         open(file='./global_system/control.dat',unit=111,action='write')
         write(111,*) 2        ! type of system, global system KS
         write(111,*) 10,1     ! calc_mode: we are doing geometry relaxation
         write(111,*) -1       ! ask ABINIT for EXX energy density and potential
         close(111)

       elseif (geo_method==2) then ! my L-BFGS 
         write(logf,'(a)') 'use my L-BFGS for geometry minimization.'
         ! remove avg(force) 
         f_noavg(1,:) = force_total(1,:) - sum(force_total(1,:))/natom 
         f_noavg(2,:) = force_total(2,:) - sum(force_total(2,:))/natom 
         f_noavg(3,:) = force_total(3,:) - sum(force_total(3,:))/natom 
         call lbfgs_dir(3,natom,geo_mhist,iter_geo,f_noavg,xcart,geo_ghist,geo_xhist,geo_z)
         xcart_tmp = xcart + reshape(geo_z,(/3,natom/))*0.2d0 
         ! send xcart to global system  
         open(file='./global_system/new_coords.dat',unit=111,action='write',form='unformatted')
         write(111) xcart_tmp
         close(111)
         open(file='./global_system/control.dat',unit=111,action='write')
         write(111,*) 2        ! type of system, global system KS
         write(111,*) 11,1     ! calc_mode=> exit scfcv(). In gstate() load new coords and re-run scfcv()
         write(111,*) -1       ! ask ABINIT for EXX energy density and potential
         close(111)
       endif 

       ! remove vpsp.dat
       open(file='./global_system/vpsp.dat',unit=111,action='read',form='unformatted')
       close(111,status='delete')

       ! remove done_end.dat
       open(file='./global_system/done_end.dat',unit=111,iostat=stat,status='old')
       close(111,status='delete')
       !
       ! Activte global system
       ! message_from_dfet and message_from_dfet2
       !
       open(file='./global_system/message_from_dfet2',unit=111,action='write',form='formatted')
       close(111,status='delete')
       open(file='./global_system/message_from_dfet',unit=111,action='write',form='formatted')
       write(111,'(a)')'RUN';close(111)
       open(file='./global_system/message_from_dfet2',unit=111,action='write',form='formatted')
       write(111,'(a)')'RUN';close(111)
       write(logf,'(a)')'global_system is activated, waiting it to finish ...'
       call flush(logf)
       !
       ! wait for global system to finish 
       !
       file_ext = .false.
       do while (.not. file_ext)
         inquire(file='./global_system/done_end.dat',exist=file_ext)
         call sleep(1)
       enddo
       file_ext = .false.
       do while (.not. file_ext)
         inquire(file='./global_system/done2_end.dat',exist=file_ext)
         call sleep(1)
       enddo
       write(logf,'(a)')'atom positions in global_system are updated.'
       write(logf,*)""
       write(logf,*)""
     endif  

     etotal_old_geo = etotal


     ! wait master to finish running global system 
     ! then we load vpsp on all the processors 
     call mpi_barrier(mpi_comm_world, ierr)

     ! remove old vpsp from system's KS potential
     do isp=1,nspin
       global_vks(:,isp) = global_vks(:,isp) - global_vpsp
     enddo 
     ! load new vpsp (all processors)
     do while(.true.)
       open(file='./global_system/vpsp.dat',unit=111,action='read',form='unformatted',IOSTAT=stat)  
       if (stat/=0) then 
         call sleep(1)
         cycle 
       endif 
       read(111,iostat=stat) global_vpsp
       if (stat/=0) then
         close(111)
         call sleep(1)
         cycle 
       else 
         close(111)
         exit
       endif 
     enddo   
     ! update vks using new vpsp (all processors) 
     do isp=1,nspin
       global_vks(:,isp) = global_vks(:,isp) + global_vpsp
     enddo 

     iter_geo = iter_geo + 1
     iter_scf = 0  ! reset iter_scf
   endif

 enddo ! loop for optimizing total system's KS potential

! if (myrank==0) then 
!    write(logf,*)'>>> job finished <<<'
!    call flush(logf)
! endif 
! call mpi_barrier(mpi_comm_world, ierr)
! call mpi_finalize(ierr);




contains 
 

!------------------------------------------------
!--   SUPPORTING SUBROUTINES BELOW   ------------
!------------------------------------------------

  subroutine clean_folder()
    if (myrank==0) then 
      call system("rm  *.xsf")
      call system("rm  log_OF_solver*.txt")
    endif
  end subroutine 


  subroutine print_energy_info

     integer, dimension(3) :: tarray
     integer :: time_end 

     call itime(tarray)
     time_end  = tarray(1)*60.d0*60.d0 + &
                 tarray(2)*60.d0 + tarray(3)

     write(logf,*)""
     write(logf,'(a,f14.8,a)')' total_ehartree:  ',total_ehart,' [abinit]'
     write(logf,'(a,f14.8)')' total_ke:        ',total_ke
     write(logf,'(a,f14.8)')' total_nlpsp:     ',total_nlpsp
     write(logf,'(a,f14.8)')' total_TS:        ',total_TS
     write(logf,'(a,f14.8)')' total_elocal:    ',total_elocal
     write(logf,'(a,f14.8)')' total_ewald:     ',total_ewald
     write(logf,'(a,f14.8)')' total_corepsp:   ',total_corepsp
     write(logf,'(a,f14.8)')' total_xc:        ',int_exc
     write(logf,'(a,f14.8)')' reg_energy:      ',reg_energy
     write(logf,'(a)')      '--------------------------------'
     write(logf,'(a,f14.8)')' etotal(sum):     ',etotal
     write(logf,'(a,f14.8)')' etotal(w/o reg): ',etotal-reg_energy
!!     write(logf,'(a,f14.8,a,es12.4)')' etotal+patch_err:',etotal+patch_err,'  patch_err: ',patch_err
     write(logf,*)''

     if (.not. direct_min) then 
       if (nspin==2) then 
         write(logf,'(a,i3,a,f14.9,a,es10.2,a,es10.2,a,f6.3,a,f6.3,a,i6)') & 
         "@it: ",iter_scf, & 
         "  etot: ",etotal, & 
         "   dE: ",de, &  
         '  |dvks|:',resid_gvks_dtmp, & 
         '  mag:  ',sum(ref_rho(:,1)-ref_rho(:,2))*dvol, & 
         '  dipole: ',sqrt(sum(dipole**2)), & 
         '  t(s):',time_end - time_start
       else 
         write(logf,'(a,i3,a,f14.9,a,es10.2,a,es10.2,a,f6.3,a,i6)') & 
         "@it: ",iter_scf, & 
         "  etot: ",etotal, &
         "   dE: ",de, &  
         '  |dvks|:',resid_gvks_dtmp, & 
         '  dipole: ',sqrt(sum(dipole**2)), & 
         '  t(s): ',time_end - time_start
       endif 
     else 
       !========== direct minimization =========
       if (nspin==2) then 
       write(logf,'(a,i4,a,f14.9,a,es10.2,a,es8.2,a,es8.2,a,f6.3,a,i7,a)') & 
         "@it: ",iter_scf, & 
         "  etot: ",etotal, &
         "   dE: ",de, &  
         '  gmax: ',gmax,' gnorm: ',gnorm, & 
         '  mag:  ',sum(ref_rho(:,1)-ref_rho(:,2))*dvol, & 
         '  t(s):',time_end - time_start,'  task: '//trim(task(1:2))
       else 
         write(logf,'(a,i4,a,f14.9,a,es10.2,a,es8.2,a,es8.2,a,i7,a)') & 
         "@it: ",iter_scf, & 
         "  etot: ",etotal, &
         "   dE: ",de, &  
         '  gmax: ',gmax,' gnorm: ',gnorm, & 
         '  t(s):',time_end - time_start,'  task: '//trim(task(1:2))
       endif 
     endif 
     write(logf,*)''
     write(logf,*)''
     write(logf,*)''
     write(logf,*)''

     ! -------- track time ---------
     call itime(tarray)
     time_start = tarray(1)*60.d0*60.d0 + tarray(2)*60.d0 + tarray(3)
     print *,''
     print *,''
     call flush(6)
     call flush(logf)
  end subroutine print_energy_info


  subroutine print_info()
     write(logf,'(a,2f12.4)')'Fermi:  ',total_fermi(1:nspin)
     write(logf,'(a,f16.6)') 'total q:',sum(ref_rho)*dvol
     if (do_global_exx>0) write(logf,'(a,f16.6)')'exx: ',exc_exx
     if (nspin==2) write(logf,'(a,f16.6)') 'mag mom:',sum(ref_rho(:,1)-ref_rho(:,2))*dvol
     filename = "rhor_up.xsf"
     call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,ref_rho(:,1))
     if (nspin==2) then 
       filename = "rhor_dn.xsf"
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,ref_rho(:,2))
       filename = "rhor_sppol.xsf"
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,ref_rho(:,1)-ref_rho(:,2))
     endif 
     if (do_global_exx>0) then   
       ! EXX energy on each atom 
       write(logf,'(a)')NEW_LINE('a')//'EXX energy on each atom:'
       do ia=1,natom
         write(logf,'(a,i4,a,f14.6)') & 
           'atom ',ia,'  EXX_atom: ',sum(exc_eps*watom(:,ia))*dvol
       enddo 
       if (iter_scf==1) filename = "dExc_dvks_up_iter0.xsf"
       if (iter_scf/=1) filename = "dExc_dvks_up.xsf"
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,total_dExc_dvks(:,1))
     endif 
     call flush(logf)
  end subroutine print_info


  subroutine display_dE_dVks
  if (nspin==1) then 
     write(logf,'(a,2es12.4,a,es12.4)') & 
    'dE_dvks: ',minval(dE_dVks),maxval(dE_dVks),'  norm: ',sqrt(sum(dE_dVks**2)*dvol)
  else 
     write(logf,'(a,2es12.4,a,es12.4)') & 
    'dE_dvks(spin up): ',minval(dE_dVks(:,1)),maxval(dE_dVks(:,1)), & 
    '  norm: ',sqrt(sum(dE_dVks(:,1)**2)*dvol)
     write(logf,'(a,2es12.4,a,es12.4)') & 
    'dE_dvks(spin dn): ',minval(dE_dVks(:,2)),maxval(dE_dVks(:,2)), & 
    '  norm: ',sqrt(sum(dE_dVks(:,2)**2)*dvol)
  endif 
  end subroutine


  subroutine print_dExc_dvs(nspin,nfft,dExc_dvks)
    implicit none 
    integer :: nspin, nfft
    real(8) :: dExc_dvks(nfft,nspin)
    if (nspin==1) then 
      write(logf,'(a,2e12.4,a,es12.4)')'(main): dExc/dvks: ', & 
      minval(dExc_dvks(:,1)),maxval(dExc_dvks(:,1)),' avg: ',sum(dExc_dvks(:,1))/nfft
    endif 
    if (nspin==2) then 
      write(logf,'(a,2e12.4,a,es12.4)')'(main): dExc/dvks(spin-up): ', & 
      minval(dExc_dvks(:,1)),maxval(dExc_dvks(:,1)),' avg: ',sum(dExc_dvks(:,1))/nfft
      write(logf,'(a,2e12.4,a,es12.4)')'(main): dExc/dvks(spin-dn): ', & 
      minval(dExc_dvks(:,2)),maxval(dExc_dvks(:,2)),' avg: ',sum(dExc_dvks(:,2))/nfft
    endif 
    call flush(logf)
  end subroutine 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  DEBUG OF-DFPT !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine debug_OF_dfpt(nspin)
   implicit none 
   integer :: nspin
   real(8) :: OF_chempot(nspin)

   write(logf,*)'enter debug_OF_dfpt().'
   if (myrank==0) then 
     rho_tmp = ref_rho
     filename ='rho_OF.xsf'
     call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,rho_tmp(:,1))
     call OF_cgmin_tot(n1,n2,n3,nfft,nspin,qvec,ucvol,global_vOF(:,1),rho_tmp,dtmp1)

     ! finite difference =====================
     vperb(:,1) = ref_rho(:,1)
     fd_step = 0.01d0
     rho_tmp = ref_rho
     !call OF_cgmin(n1,n2,n3,nfft,nspin,qvec,ucvol,global_vOF(:,1)+vperb(:,1)*fd_step,rho_tmp,dtmp1)
     call OF_cgmin_tot(n1,n2,n3,nfft,nspin,qvec,ucvol,global_vOF(:,1)+vperb(:,1)*fd_step,rho_tmp,dtmp1)

     rho_tmp2 = ref_rho
     !call OF_cgmin(n1,n2,n3,nfft,nspin,qvec,ucvol,global_vOF(:,1)-vperb(:,1)*fd_step,rho_tmp2,dtmp2)
     call OF_cgmin_tot(n1,n2,n3,nfft,nspin,qvec,ucvol,global_vOF(:,1)-vperb(:,1)*fd_step,rho_tmp2,dtmp2)

     rho_perb = (rho_tmp - rho_tmp2)/2.d0/fd_step
     write(logOF,'(a,2es12.4)')'rho_perb [finite diff]: ',minval(rho_perb),maxval(rho_perb)
     write(logOF,'(a,f12.8)')'ef_perb [finite diff]:  ',(dtmp1-dtmp2)/2.d0/fd_step
     filename = "rho_perb_finite_diff.xsf"
     call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,rho_perb(:,1))

     ! analytical 
     call OF_dfpt(nfft,cell_nfft,ucvol,qvec,ref_rho(:,1),vperb(:,1),rho_perb(:,1),dmu)
     filename = "rho_perb_dfpt.xsf"
     call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,rho_perb(:,1))

     call OF_dfpt_inv(nfft,cell_nfft,ucvol,qvec,ref_rho(:,1),rho_perb(:,1),vperb(:,1))
     call OF_dfpt(nfft,cell_nfft,ucvol,qvec,ref_rho(:,1),vperb(:,1),rho_perb(:,1),dmu)


     ! debug drho/dN code 
     !========================
     !finite difference 
     rho_tmp  = ref_rho/sum(ref_rho*dvol)*(sum(ref_rho*dvol)+0.02d0)
     call OF_cgmin_tot(n1,n2,n3,nfft,nspin,qvec,ucvol,global_vOF(:,1),rho_tmp,dtmp2)
     rho_tmp2 = ref_rho/sum(ref_rho*dvol)*(sum(ref_rho*dvol)-0.02d0)
     call OF_cgmin_tot(n1,n2,n3,nfft,nspin,qvec,ucvol,global_vOF(:,1),rho_tmp2,dtmp2)
     rho_perb = (rho_tmp-rho_tmp2)/0.02/2.d0
     write(logOF,'(a,2es12.4)')'drho/dN [finite diff]: ',minval(rho_perb),maxval(rho_perb)
     filename = "drho_dN_fd.xsf"
     call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,rho_perb(:,1))

     ! analytical 
     call OF_drho_dN_vfix(nfft,cell_nfft,ucvol,qvec,ref_rho(:,1),rho_perb(:,1))
     filename = "drho_dN.xsf"
     call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,rho_perb(:,1))
   endif 
   call mpi_barrier(mpi_comm_world,ierr)
   write(logf,*)'stop in debug_OF_dfpt()'
   stop
  end subroutine  debug_OF_dfpt




  subroutine debug_my_lbfgs_dir()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! DEBUG lbfgs_chen 
 Amat_tmp(1,:) = (/1, 2, 3, 4, 5/)
 Amat_tmp(2,:) = (/6, 7, 8, 9, 10/)
 Amat_tmp(3,:) = (/5, 4, 3, 2, 1/)
 Amat_tmp(4,:) = (/10,9, 8, 7, 6/)
 Amat_tmp(5,:) = (/5, 4, 3, 2, 1/)

 Amat_tmp = matmul(Amat_tmp,transpose(Amat_tmp))
 bvec_tmp(:,1) = (/1,2,3,4,5/)
 xvec_tmp(:,1) = (/2,4,6,8,10/)

 do iter_scf=1,50
   
   gvec_tmp = matmul(Amat_tmp,xvec_tmp) - bvec_tmp
   dtmp = 0.5d0*sum(xvec_tmp*matmul(Amat_tmp,xvec_tmp)) - sum(bvec_tmp*xvec_tmp)

   write(logf,*)'dtmp: ',dtmp
   write(logf,*)'gvec: ',gvec_tmp

   call lbfgs_dir(1,5,db_mhist,iter_scf,gvec_tmp,xvec_tmp, & 
                  db_ghist,db_xhist,db_lbfgs_z)
   
   write(logf,*)'z: ',db_lbfgs_z

   xvec_tmp = xvec_tmp + db_lbfgs_z*1.d0
 enddo 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  endsubroutine 



!include "debug_line_search.f90" 
!include "make_global_vks_debug.f90"

end program ecda
