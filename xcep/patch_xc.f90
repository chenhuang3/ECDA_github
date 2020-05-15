
subroutine patch_xc(natom,natom_clu,natom_env,nfft,nspin,mxband,myrank,cell_nfft,& 
                    atom_info,cell_acell,plot_sub_rho,logf,mlist,ilist,nlist,nshare,ishare,qvec, & 
                    solve_vxc,nopt_uemb,reg_vxc,zp_alpha,zp_coeff,zp_conv,relax_geo,rprimd,tsmear,ucvol, & 
                    dft_flavor,vac,xcart,znucl,global_fermi,cov_radius,melement,do_force,do_xcep, & 
                    dvol,iter_scf,total_fermi,watom,dfermi_dvks,gvxc_init,zvec_init,cgW_tol, & 
                    ref_rho,u,den_diff,do_exx,zvec,global_vks,global_vOF,patch_rho,int_exc, & 
                    force_pen,force_abinit,force_total,dExc_dvks,only_lda)

   use mpi
   use interface_funcs
   use comm_data

   implicit none 
   
   logical :: do_force,plot_sub_rho, & 
              solve_vxc, only_lda

   integer :: natom, do_exx, nfft, nspin, & 
              cell_nfft(3), myrank, iter_scf, & 
              melement, do_xcep, relax_geo, & 
              global_fermi, mxband, & 
              nopt_uemb,gvxc_init,zvec_init, & 
              logf, natom_clu, natom_env, & 
              nshare(natom), dft_flavor, & 
              ishare(mlist,natom), & 
              mlist, nlist(natom), ilist(mlist,natom)

   real(8) :: cgW_tol,ref_rho(nfft,nspin), &
              qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), & 
              vhart(nfft), reg_vxc, zp_alpha,  & 
              rprimd(3,3), &  
              xcart(3,natom), znucl(natom), zp_coeff,zp_conv, & 
              cov_radius( melement ), dvol, tsmear, ucvol, vac,&
              prj_vks(nfft,nspin,2), &   ! one for env, another for cluster 
              vks_clu(nfft,nspin), & 
              vks_env(nfft,nspin), & 
              dExc_dVks(nfft,nspin), & 
              zvec(nfft,nspin), & 
              yvec(nfft,nspin), & 
              sub_rhor(nfft,nspin,2), & 
              watom(nfft,natom), & 
              atom_info(5,natom), &    ! (subsystem_index, type, x,y,z) x,y,z are in angstrom 
              cell_acell(3,3), &  ! the a b c of cell
              total_fermi(nspin), & 
              force_abinit(3,natom), & 
              force_total(3,natom), & 
              force_pen(3,natom),  &
              global_vks(nfft,nspin), & 
              global_vOF(nfft,nspin), & 
              patch_rho(nfft,nspin), & 
              u(nfft,nspin), &
              dfermi_dvks(nfft,nspin), & 
              int_exc, &
              total_exc(nfft), & 
              total_vxc(nfft,nspin), & 
              exx_mix

   ! local vars 
   logical :: trans, & 
              exact_regularization
   character(len=100) :: filename,fname,ss,fname2,dtime
   integer,parameter :: npulay = 5
   integer :: n1,n2,n3,iter_z, system_type, & 
              iter_total_vxc,ierr,ii,kk,isp,k,j_atom, & 
              maxiter_solvez, & 
              i, shared, & 
              status(MPI_STATUS_SIZE)
   real(8) :: dtmp,g1,g2,gN, & 
              exc_exx,  & 
              exc_envcorr,  & 
              xpbe_atom, cpbe_atom, & 
              lda_eps_tot(nfft),  & 
              lda_eps_clu(nfft), & 
              xpbe_eps_clu(nfft), & 
              cpbe_eps_clu(nfft), & 
              vxpbe_clu(nfft,nspin), & 
              vcpbe_clu(nfft,nspin), & 
              clu_weight(nfft), & 
              env_weight(nfft), & 
              force_xc1(3,natom), & 
              force_xc2(3,natom), & 
              force_xcN(3,natom), & 
              force_xc(3,natom), & 
              den_diff(nfft,nspin),  & 
              arr_tmp(nfft,nspin), & 
              rho_perb(nfft,nspin), & 
              zsum(nfft,nspin), & 
              depsLDA_dN_clu(nfft), &
              deps_dN_clu(nfft,nspin), &
              drho_dN_clu(nfft,nspin), & 
              drho_dN_env(nfft,nspin), & 
              nele(natom,2), & 
              nele_mag(natom,2), & 
              q_clu(nspin), & 
              q_env(nspin),& 
              fermi_sub(nspin,natom,2), & ! (:,:,1) => cluster, (:,:,2) => env
              vxclda_clu(nfft,nspin), & 
              vxclda_tot(nfft,nspin), & 
              tmp_vxc(nfft,nspin), & 
              tmpvec(nfft), & 
              tmpvec2(nfft,nspin), & 
              yvec2(nfft,nspin), & 
              z_chi_clu(nfft,nspin), & 
              z_chi_env(nfft,nspin), & 
              wenv_chi(nfft,nspin), & 
              gvec_old(nfft,nspin), & 
              gvec_change(nspin), & 
              sigma_g(nfft,nspin), & 
              sigma_x(nfft,nspin), & 
              z_change(nspin), & 
              shift, vd(nfft), &
              dNclu_dR, dNenv_dR, & 
              vxc_change(nspin), & 
              pulay_array(nfft*nspin,npulay,2), & 
              pulay_beta

   ! forces 
   integer :: id, ipt, iR, p1, p2
   real(8) :: hvec(nfft,nspin), & 
              dAtomW_dR(3,nfft,natom), & 
              dAtomW_dR2(3,nfft,natom), & 
              dClusterW_dR(3,nfft,natom), & 
              dClusterW_dR2(3,nfft,natom), & 
              cluster_exc_eps(nfft,nspin)

   integer :: local_do_exx, & 
              local_do_xcep

   logical :: last_iter


   ! PBE0
   if (dft_flavor==-102) exx_mix = 0.25d0

   total_exc = 0.0d0 
   total_vxc = 0.0d0
   fermi_sub = 0.d0 
   n1 = cell_nfft(1)
   n2 = cell_nfft(2)
   n3 = cell_nfft(3)


   ! write ref_density 
   if (myrank==0) then 
     call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,'ref_rho_up.xsf',ref_rho(:,1))
   endif 

   !===============================
   !  Big loop over all the atoms
   !===============================
   do j_atom = 1,natom     

     if (j_atom /= myrank+1) cycle

     if (myrank==0) then 
       write(logf,'(a)')'enter patch_xc() ...'
       write(logf,'(a,i3)')'xc_patch => working on atom: ',j_atom
       call flush(logf)
     endif


     ! make cluster weight, loop over all atoms
     clu_weight = watom(:,j_atom)
     do kk=1,nlist(j_atom)
       clu_weight = clu_weight + watom(:,ilist(kk,j_atom))
     enddo
     env_weight = 1.d0 - clu_weight


     ! new electron numbers for cluster and env
     ! send them to cluster and environment 
     do isp=1,nspin
       q_clu(isp) = dvol*sum(ref_rho(:,isp)*clu_weight)
       q_env(isp) = dvol*sum(ref_rho(:,isp)*env_weight)
     enddo
     if (global_fermi>0) then 
       write(logf,'(a)')'set charge in cluster and env.'
       call set_clu_env_q(nspin,j_atom,q_clu,q_env,natom_clu,natom_env)
     endif 


     ! projecting global vks to subsystem vks
     ! For xc_patching, we have only two subsystems for the case of xc_path 
     if (myrank==0) write(logf,'(a,f12.4)')'vac_shift: ',vac; call flush(logf)
     call proj_dump_cluster_env_vks(j_atom,nspin,nfft,natom,total_fermi, & 
                                    prj_vks,vac,clu_weight,env_weight,global_vks)

     ! For the case of OF-DFT for env
     if (do_envOF) then  
       do isp=1,nspin
         prj_vks(:,isp,2) = (global_vOF(:,isp)-total_fermi(isp)-vac)*env_weight
       enddo
     endif 

     ! prepare embedding potential 
     if (iter_scf==1 .and. relax_geo/=1) then 
       den_diff = 0.0d0
       u = 0.d0
     endif 

     if (myrank==0) write(logf,'(a)')'check log_cluster_* files for details of other clusters ...'
     call flush(logf)


     !===========================
     ! density-based embedding 
     !===========================
     if (.not. only_lda) then 
       call dfet_zp(mxband,nfft,nspin,natom,j_atom,ucvol,logf,myrank,do_exx,do_xcep, & 
              zp_alpha,nopt_uemb,global_fermi,cell_nfft,cgW_tol,qvec,zp_coeff,zp_conv, & 
              ref_rho,global_vks,global_vOF,prj_vks,q_clu,q_env,env_weight,tsmear,sub_rhor,fermi_sub,& 
              natom_clu,natom_env,nele_mag,nele,exc_exx,cluster_exc_eps,u,den_diff,yvec, & 
              drho_dN_clu,drho_dN_env,deps_dN_clu)

       ! PBE0? 
       if (dft_flavor==-102) then 
         yvec = exx_mix*yvec                     ! 1/4 for PBE0
         deps_dN_clu = exx_mix*deps_dN_clu       ! 1/4 for PBE0
         call calculate_x_pbe_eps(n1,n2,n3,nspin,sub_rhor(:,:,1),qvec,dvol,vxpbe_clu,xpbe_eps_clu)
         call calculate_c_pbe_eps(n1,n2,n3,nspin,sub_rhor(:,:,1),qvec,dvol,vcpbe_clu,cpbe_eps_clu)
         xpbe_atom = dvol*sum(watom(:,j_atom)*xpbe_eps_clu)
         cpbe_atom = dvol*sum(watom(:,j_atom)*cpbe_eps_clu)
         exc_exx   = exx_mix*exc_exx + (1.d0-exx_mix)*xpbe_atom + cpbe_atom
       endif 

       ! LDA correction for cluster XC energy
       exc_envcorr = 0.d0 
       if (do_env_correction) then 
         call calculate_xc_pw92LSDA_eps(n1,n2,n3,nspin,sub_rhor(:,:,1),dvol,vxclda_clu,dtmp,lda_eps_clu)
         call calculate_xc_pw92LSDA_eps(n1,n2,n3,nspin,ref_rho,dvol,vxclda_tot,dtmp,lda_eps_tot)
         ! correction from env (LDA)
         exc_envcorr = exc_envcorr + dvol*sum(watom(:,j_atom)*(lda_eps_tot-lda_eps_clu)) 
         exc_exx     = exc_exx + exc_envcorr
       else 
         vxclda_clu = 0.d0 
         vxclda_tot = 0.d0 
         lda_eps_tot = 0.d0 
         lda_eps_clu = 0.d0 
       endif 

       ! get patched XC energy 
       write(715,'(a,f14.6)')'atomic EXX energy: ',exc_exx
       write(715,'(a,f14.6)')'exc_envcorr      : ',exc_envcorr
       write(715,'(a)') 'waiting others to finish ...'; call flush(715)
       call mpi_allreduce(mpi_in_place,exc_exx,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       call mpi_allreduce(mpi_in_place,exc_envcorr,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       int_exc = exc_exx

       if (myrank==0) write(logf,'(a,f14.6)')'patched EXX energy: ',exc_exx
       if (myrank==0) write(logf,'(a,f14.6,a,f14.6,a)') & 
         'exc_envcorr:        ',exc_envcorr,', ',exc_envcorr/dble(natom),'  per atom'
     endif 

     ! cluster and environment's KS potential 
     vks_clu = prj_vks(:,:,1) + u 
     vks_env = prj_vks(:,:,2) + u 


     !==================================================
     ! System is treated by LDA (PW92), for debug only
     !==================================================
     if (only_lda) then 
       call calculate_xc_pw92LSDA_eps(n1,n2,n3,nspin,ref_rho,dvol,vxclda_tot,dtmp,lda_eps_tot)
       exc_exx = 0.d0
       exc_exx = exc_exx + dvol*sum(watom(:,j_atom)*lda_eps_tot)  ! correction with PW92LDA 
       call mpi_allreduce(mpi_in_place,exc_exx,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       int_exc = exc_exx
       if (myrank==0) write(logf,'(a,f14.6)')'patched LDA (PW92) energy: ',exc_exx

       ! get dExc/dvks
       do isp=1,nspin
         tmp_vxc(:,isp) = watom(:,j_atom)*vxclda_tot(:,isp)  ! due to total system's LDA xc
       enddo
       call mpi_allreduce(MPI_IN_PLACE,tmp_vxc,nfft*nspin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

       ! apply chi_total 
       if (myrank==0) then 
         system_type = 0 ! 0: global system
         call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,tmp_vxc,rho_perb)
         call remove_mean(nspin,nfft,rho_perb)
       else 
         rho_perb = 0.d0
       endif

       ! send to all processors
       call mpi_allreduce(rho_perb,dExc_dvks,nfft*nspin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

       ! compute force (xc force is zero because rho_total is fixed for fixed vks)
       force_total = force_abinit + force_pen
       return 
     endif 


     !===================================================
     ! Following codes are for only_lda = false
     !===================================================
     if (myrank==0) then 
       if (nspin==2) then 
         write(logf,'(a,2es12.4,a,es12.4,a)')'min/max(yvec_alpha): ',minval(yvec(:,1)),maxval(yvec(:,1)), & 
           ' sum(yvec):',sum(yvec(:,1))*dvol,' (process 0)'
         write(logf,'(a,2es12.4,a,es12.4,a)')'min/max(yvec_beta):  ',minval(yvec(:,2)),maxval(yvec(:,2)), & 
           ' sum(yvec):',sum(yvec(:,2))*dvol,' (process 0)'
       else 
         write(logf,'(a,2es12.4,a,es12.4,a)')'min/max(yvec): ',minval(yvec(:,1)),maxval(yvec(:,1)), & 
           ' sum(yvec):',sum(yvec(:,1))*dvol,' (process 0)'
       endif 
       call flush(logf)
     endif 
     ! check yvec 
     if (abs(sum(yvec)*dvol)>1e-4) then 
       write(logf,*)'abs(sum(yvec)*dvol)>1e-4, stop! abs(sum(yvec)*dvol)>1e-4: ',abs(sum(yvec)*dvol)>1e-4
       stop
     endif 


     ! plot cluster densigty 
     if (plot_sub_rho) then 
       ! write embedding potential 
       write(ss,*)j_atom
       filename = "uemb_"//trim(adjustl(ss))//".xsf"
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,u(:,1))
       ! write cluster density 
       filename = "rho_cluster_"//trim(adjustl(ss))//".xsf"
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,sub_rhor(:,1,1))
       filename = "rho_env_"//trim(adjustl(ss))//".xsf"
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,sub_rhor(:,1,2))
     endif 


     ! system's density obtained from patching, 
     ! note that there are errors associated with clusters' densities 
     do isp =1,nspin 
       patch_rho(:,isp) = watom(:,j_atom)*sub_rhor(:,isp,1)
     enddo 
     call mpi_allreduce(mpi_in_place,patch_rho,nfft*nspin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
     write(logf,*)'Q(patch_rho): ',sum(patch_rho)*dvol


     ! get cluster's electron numbers 
     call mpi_allreduce(MPI_IN_PLACE,nele, natom*2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(MPI_IN_PLACE,nele_mag, natom*2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(MPI_IN_PLACE,fermi_sub,nspin*natom*2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
     call print_cluster_info()


     !========================
     ! Compute dExc/dVks
     !========================

     !======= compute yvec2 due to env's LDA XC correction ===========
     if (do_env_correction) then 
       write(logf,'(a)')NEW_LINE('a')//'computing yvec2 due to LDA XC energy of cluster ...';call flush(logf)

       select case(dft_flavor)
       case (-101) ! EXX
         write(logf,'(a)')'cluster is treated by EXX.';call flush(logf)
         do isp=1,nspin
           tmpvec2(:,isp) = watom(:,j_atom)*vxclda_clu(:,isp)
         enddo
       case (-102) ! PBE0
         write(logf,'(a)')'cluster is treated by PBE0.';call flush(logf)
         do isp=1,nspin
           tmpvec2(:,isp) = watom(:,j_atom)*(-(1.d0-exx_mix)*vxpbe_clu(:,isp)-vcpbe_clu(:,isp)+vxclda_clu(:,isp))
         enddo
       endselect 

       system_type = 1 ! cluster
       call apply_chi(-1,system_type,j_atom,nfft,nspin,vks_clu,tmpvec2,yvec2)
       call remove_mean(nspin,nfft,yvec2)
       write(logf,'(a,2e20.8)')'yvec2: ',minval(yvec2(:,1)),maxval(yvec2(:,1))
       call flush(logf)
     else 
       write(logf,'(a)')NEW_LINE('a')//'no env correction, set yvec2 to 0.'
       yvec2 = 0.d0 
       call flush(logf)
     endif 
     
     !=======================================================
     ! Now we have yvec and yvec2, let's compute dExc/dVks
     !=======================================================

     !======= Part I of dExc/dvks ========================================
     ! This is due to the derivative of Exc with respect to vks_cluster
     ! due to the change of system's vks with v_emb fixed 
     !
     do isp=1,nspin
       dExc_dvks(:,isp) = (yvec(:,isp)-yvec2(:,isp))*clu_weight
     enddo
     write(715,'(a,2es12.4)')'contribution to dExc_dvks (p1):',minval(dExc_dvks),maxval(dExc_dvks)
     call mpi_allreduce(MPI_IN_PLACE,dExc_dvks,nfft*nspin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)


     !======= Part II of dExc/dvks ================================
     ! This is the derivative of Exc with respect to vks_cluster 
     ! however vks is fixed but the v_emb is change, since v_emb is a functional of vks.
     !
     if (myrank==0) write(logf,'(a)')"computing z = y*(\chi_cluster + \chi_env)^(-1) ... "
     maxiter_solvez = 20
     call solve_z(natom,natom_clu,natom_env,nfft,nspin, & 
                  maxiter_solvez,j_atom,myrank,logf,cell_nfft,ucvol,qvec,yvec-yvec2, & 
                  vks_clu,vks_env,sub_rhor,zp_coeff,zp_alpha,zvec_init,zvec)

     write(715,'(a)') 'main(): waiting others to finish ...';call flush(715)

     ! sum over z = y*(chi_clu+chi_env)^{-1} from all processors 
     call mpi_allreduce(zvec,zsum,nfft*nspin, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
     write(715,'(a)') 'main(): done mpi_allreduce() to sum up zvec from all processes.';call flush(715)


     !==========================
     ! add zvec*\chi_total 
     !==========================
     if (myrank==0) then 
       write(logf,'(a)')'computing z*\chi_total (process 0) ...'; call flush(logf)
       write(715,'(a)') 'computing z*\chi_total (process 0) ...'; call flush(715)
       system_type = 0 ! 0: global system
       call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,zsum,rho_perb)
       call remove_mean(nspin,nfft,rho_perb)
     else 
       rho_perb = 0.d0
     endif
     write(715,'(a,2es12.4,a)')'contribution to dExc_dvks: ', & 
       minval(rho_perb),maxval(rho_perb),' [zero for other nodes]'
     ! broadcast rho_perb to all processes
     call mpi_allreduce(rho_perb, tmp_vxc, nfft*nspin, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
     dExc_dvks = dExc_dvks + tmp_vxc


     !===================================================
     ! add -zvec*\chi_cluster*weight_cluster
     !
     ! Note that: the ECDA+force paper is wrong about 
     !            the order of Weight and chi_cluster,
     !            however, the original code is correct.
     !===================================================
     write(logf,'(a)')'computing z*\chi_cluster (all processes) ...'; call flush(logf)
     write(715,'(a)') 'computing z*\chi_cluster (all processes) ...'; call flush(715)
     system_type = 1 ! cluster
     call apply_chi(-1,system_type,j_atom,nfft,nspin,vks_clu,zvec,rho_perb)
     call remove_mean(nspin,nfft,rho_perb)
     z_chi_clu = rho_perb

     ! apply weight_cluster
     do isp=1,nspin
       rho_perb(:,isp) = rho_perb(:,isp)*clu_weight
     enddo 
     write(715,'(a,2es12.4)')'contribution to dExc_dvks:',minval(rho_perb),maxval(rho_perb)
     call mpi_allreduce(rho_perb,tmp_vxc,nfft*nspin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
     dExc_dvks = dExc_dvks - tmp_vxc


     !=================================================
     ! add -zvec*\chi_env*weight_env
     !
     ! Note that: the ECDA/force paper is wrong about the order of Weight and chi_cluster
     !=================================================
     if (natom_env>0) then 
       write(logf,'(a)')'computing z*\chi_env (all processes) ...'; call flush(logf)
       write(715 ,'(a)')'computing z*\chi_env (all processes) ...'; call flush(715)
       if (.NOT. do_envOF) then 
         ! env is treated by KSDFT
         system_type = 2 ! env
         call apply_chi(-1,system_type,j_atom,nfft,nspin,vks_env,zvec,rho_perb)
         call remove_mean(nspin,nfft,rho_perb)
       else 
         ! env is treated by OF-DFT
         do isp=1,nspin 
           call OF_dfpt(nfft,cell_nfft,ucvol,qvec,sub_rhor(:,isp,2),& 
                        zvec(:,isp),rho_perb(:,isp),dtmp)
         enddo
       endif 
       z_chi_env = rho_perb

       ! apply weight_env
       do isp=1,nspin
         rho_perb(:,isp) = rho_perb(:,isp)*env_weight
       enddo 
     else
       write(715 ,'(a)')'natom_env=0, skip computing z*\chi_env '; call flush(715)
       rho_perb  = 0.d0 
       z_chi_env = 0.0d0
     endif 

     if (.not. do_envOF) then 
       write(715,'(a,2es12.4)')'contribution to dExc_dvks:',minval(rho_perb),maxval(rho_perb)
       call mpi_allreduce(mpi_in_place,rho_perb,nfft*nspin,mpi_double_precision, & 
                          mpi_sum,mpi_comm_world,ierr)
     endif


     !=================================
     ! extra steps for OFDFT case 
     if (do_envOF) then 
       tmp_vxc = rho_perb

       call mpi_allreduce(mpi_in_place,rho_perb,nfft*nspin, & 
         mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

       ! only apply system's chi_OF and chi_KS
       ! on the master node and then send to all other nodes 
       !
       if (myrank==0) then 
         ! apply inverse of system's OF-DFT linear response 
         do isp=1,nspin
           arr_tmp(:,isp) = rho_perb(:,isp)
           call OF_dfpt_inv(nfft,cell_nfft,ucvol,qvec,ref_rho(:,isp), & 
                            arr_tmp(:,isp),rho_perb(:,isp))
         enddo

         ! apply system's KS-DFT linear response 
         !
         system_type = 0 ! 0: global system
         arr_tmp = rho_perb
         call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,arr_tmp,rho_perb)
         call remove_mean(nspin,nfft,rho_perb)
       endif 
       call send_data_master_to_all(nfft*nspin,myrank,rho_perb)
     endif

     dExc_dvks = dExc_dvks - rho_perb
 

     !===============================================
     ! Terms due to change of system's fermi level 
     !
     ! + zvec*\chi_clu*|w_clu><dfermi_dvks| 
     ! + zvec*\chi_env*|w_env><dfermi_dvks|
     !===============================================
     do isp=1,nspin
        tmp_vxc(:,isp) = sum(z_chi_clu(:,isp)*clu_weight)*dvol*dfermi_dvks(:,isp)   &
                        +sum(z_chi_env(:,isp)*env_weight)*dvol*dfermi_dvks(:,isp)
     enddo 
     call mpi_allreduce(MPI_IN_PLACE,tmp_vxc,nfft*nspin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
     dExc_dvks = dExc_dvks + tmp_vxc


     !=========================================================
     ! Terms due to change of cluster and env electron numbers
     !
     ! add - z*(drho_env/dN_env)*<w_clu|*\chi
     !     - z*(drho_clu/dN_clu)*<w_env|*\chi
     !
     ! also include 
     !  1) the part that is due to change of N_cluster
     !  2) the part due to cluster's LDA XC (change N_cluster)
     !  3) the part due to system's LDA XC (change rho_total)
     !=========================================================

     ! compute drho/dN of environment (OF-DFT)
     if (do_envOF .and. natom_env>0) then 
       do isp=1,nspin
         call OF_drho_dN_vfix(nfft,cell_nfft,ucvol,qvec,sub_rhor(:,isp,2),drho_dN_env(:,isp))
       enddo
     else 
       drho_dN_env = 0.d0
     endif 

     do isp=1,nspin
       tmp_vxc(:,isp) = & 
         - sum(zvec(:,isp)*drho_dN_env(:,isp))*dvol*env_weight &
         - sum(zvec(:,isp)*drho_dN_clu(:,isp))*dvol*clu_weight & 
         + sum(watom(:,j_atom)*deps_dN_clu(:,isp))*dvol*clu_weight & 
         + watom(:,j_atom)*vxclda_tot(:,isp)  ! due to total system's LDA xc
       
       ! due to cluster's LDA XC energy
       select case (dft_flavor)
       case (-101) ! EXX
         tmp_vxc(:,isp) = tmp_vxc(:,isp) & 
           - sum(vxclda_clu(:,isp)*drho_dN_clu(:,isp)*watom(:,j_atom))*dvol*clu_weight 
       case (-102) ! PBE0
         tmp_vxc(:,isp) = tmp_vxc(:,isp) & 
           - sum((-(1.d0-exx_mix)*vxpbe_clu(:,isp)-vcpbe_clu(:,isp)+vxclda_clu(:,isp)) & 
                  *drho_dN_clu(:,isp)*watom(:,j_atom))*dvol*clu_weight
       endselect 
     enddo

     call mpi_allreduce(MPI_IN_PLACE,tmp_vxc,nfft*nspin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

     ! apply chi_total 
     if (myrank==0) then 
       system_type = 0 ! 0: global system
       call apply_chi(-1,system_type,-1,nfft,nspin,global_vks,tmp_vxc,rho_perb)
       call remove_mean(nspin,nfft,rho_perb)
     else 
       rho_perb = 0.d0
     endif

     ! send to all processors
     call mpi_allreduce(rho_perb,tmp_vxc,nfft*nspin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
     dExc_dvks = dExc_dvks + tmp_vxc


     !========================================================
     ! PART III
     ! change of v_cluster due to change of chemical potential 
     !========================================================
     do isp=1,nspin
       tmp_vxc(:,isp) = -sum((yvec(:,isp)-yvec2(:,isp))*clu_weight)*dvol*dfermi_dvks(:,isp)
     enddo 
     call mpi_allreduce(MPI_IN_PLACE,tmp_vxc,nfft*nspin,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
     dExc_dvks = dExc_dvks + tmp_vxc


     ! writing dExc_dvks XSF files 
     if (myrank==0) then 
       if (iter_scf==1) then 
          filename = 'dExc_dvks_up_patched_iter0.xsf'
       else 
          filename = "dExc_dvks_up_patched.xsf"
       endif
       call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,dExc_dvks(:,1))
     endif 



     !====================================
     ! compute Non-Hellmann-Feynman force 
     !====================================
     if (do_force) then 
       if (nspin==2) then 
         write(logf,*)'force has not been coded for spin polarized case, mainly due to the following '
         write(logf,*)'variables are not spin-polarized yet: lda_eps_clu, lda_eps_tot, xpbe_eps_clu, cpbe_eps_clu'
         write(logf,*)'also forces has not been coded for PBE0 case yet.'
         stop
       endif 

       ! commpute the change of atom and cluster weights due to moving atom j_atom
       call calc_dwatom_dR(logf,natom,nspin,nfft,natom,melement,atom_info, &
                           cell_nfft,xcart,cell_acell,znucl,rprimd,cov_radius, &
                           mlist,nlist,ilist,j_atom,dAtomW_dR2,dClusterW_dR2)
       
       ! note that, for processor j_atom, we only have zvec, yvec, etc for that cluster j_atom. 
       ! Thus, we create a new array dClusterW_dR which is d(cluter_weight[j_atom])/dR for R = 1,2.,,, natom 
       ! Using MPI commands, we need to exchange the data between processors.
       do p2 = 1,natom
         do p1 = 1,natom 
           ! send data process #p1 --> process #p2
           if (p1/=p2) then 
             if (myrank == p1-1) then
               call mpi_send(dClusterW_dR2(:,:,p2),3*nfft,MPI_DOUBLE_PRECISION,p2-1,1,MPI_COMM_WORLD,ierr)
               call mpi_send(dAtomW_dR2(:,:,p2),3*nfft,MPI_DOUBLE_PRECISION,p2-1,1,MPI_COMM_WORLD,ierr)
             elseif (myrank == p2-1)  then 
               call mpi_recv(dClusterW_dR(:,:,p1),3*nfft,MPI_DOUBLE_PRECISION,p1-1,1,MPI_COMM_WORLD,status,ierr)
               call mpi_recv(dAtomW_dR(:,:,p1),3*nfft,MPI_DOUBLE_PRECISION,p1-1,1,MPI_COMM_WORLD,status,ierr)
             endif 
             call mpi_barrier(mpi_comm_world,ierr)
           else
             if (myrank==p1-1) then
               dClusterW_dR(:,:,p2) = dClusterW_dR2(:,:,p2)
               dAtomW_dR(:,:,p2) = dAtomW_dR2(:,:,p2)
             endif 
           endif 
         enddo
       enddo 

       ! loop the atoms on which forces to be calculated 
       ! each node is in charge of its own XC energy density 
       do iR = 1,natom 
         ! x, y, and z directions
         do id = 1,3 
           g1 = sum((sum(cluster_exc_eps,2)-lda_eps_clu+lda_eps_tot)*dAtomW_dR(id,:,iR))*dvol ! change of atom weight
           g2 = 0.d0 
           gN = 0.d0  

           do isp=1,nspin
             dNclu_dR =  sum(ref_rho(:,isp)*dClusterW_dR(id,:,iR))*dvol 
             dNenv_dR = -sum(ref_rho(:,isp)*dClusterW_dR(id,:,iR))*dvol
             vd = global_vks(:,isp)-total_fermi(isp)-vac
 
             ! force due to the change of cluster's KS potential
             g2 = g2 & 
                + sum((yvec(:,isp)-yvec2(:,isp))*vd*dClusterW_dR(id,:,iR))*dvol  &  ! F_A
                - sum(z_chi_clu(:,isp)*vd*  dClusterW_dR(id,:,iR) )*dvol         &  ! F_c
                - sum(z_chi_env(:,isp)*vd*(-dClusterW_dR(id,:,iR)))*dvol         &  ! F_c
                - sum(zvec(:,isp)*drho_dN_clu(:,isp))*dvol*dNclu_dR              &  ! F_c
                - sum(zvec(:,isp)*drho_dN_env(:,isp))*dvol*dNenv_dR 

             ! force due to change of cluster's electron number 
             depsLDA_dN_clu = vxclda_clu(:,isp)*drho_dN_clu(:,isp)
             gN = gN + sum(watom(:,j_atom)*(deps_dN_clu(:,isp)-depsLDA_dN_clu))*dvol*dNclu_dR ! F_N
           enddo ! spin

           ! sum up contributions from all clustsers
           call mpi_allreduce(-g1,force_xc1(id,iR),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
           call mpi_allreduce(-g2,force_xc2(id,iR),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
           call mpi_allreduce(-gN,force_xcN(id,iR),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
         enddo ! x,y,z
       enddo ! iR (moving atom)

       force_xc = force_xc1 + force_xc2 + force_xcN 
       force_total = force_xc + force_abinit + force_pen
       
       if (myrank==0) then 
         call display_force(logf,natom,force_xc,force_xc1,force_xc2,force_xcN,force_pen,force_abinit,force_total)
       endif  
     endif 
     ! end of force calculation
     !=========================

     exit 

   enddo ! j_atom, loop for all the atoms (for job_type=300)

   call mpi_barrier(mpi_comm_world,ierr)
   write(logf,'(a)')'done patch_xc()'
   call flush(logf)







contains 


subroutine print_cluster_info 
  ! display information
  if (myrank==0) then 
    write(logf,'(a)')" "
    if (nspin==2) then 
      write(logf,'(a)') & 
       "atom_id   q_clu    q_env      q_tot      mag_clu    mag_env    mag_tot  clu_Ef_up  clu_Ef_dn  env_Ef_up  env_Ef_dn" 
    else 
      write(logf,'(a)') & 
       "atom_id   q_clu    q_env      q_tot    clu_fermi  env_fermi" 
    endif 
    do ii=1,natom
      if (nspin==2) then 
        write(logf,'(i3,a,f10.4,a,f10.4,a,f10.4,a,f10.4,a,f10.4,a,f10.4,a,4f10.4)') & 
        ii," ",nele(ii,1)," ",nele(ii,2)," ",sum(nele(ii,:))," ",& 
        nele_mag(ii,1)," ", nele_mag(ii,2)," ",sum(nele_mag(ii,:)), " ", fermi_sub(:,ii,1),fermi_sub(:,ii,2)
      else 
        write(logf,'(i3,a,f10.4,a,f10.4,a,f10.4,2f10.4)') & 
        ii," ",nele(ii,1)," ",nele(ii,2)," ",sum(nele(ii,:)),fermi_sub(:,ii,1),fermi_sub(:,ii,2)
      endif 
    enddo
    write(logf,'(a)')'partition of total system is done for all cluster-env pairs'
    call flush(logf)
    call flush(6)
  endif
end subroutine print_cluster_info
 
end subroutine patch_xc
