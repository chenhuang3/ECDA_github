!
! partition the cluster and env according to the 
! density functional embedding theory 
! 
!  outputs: 
!     1) u - embedding potential 
!     2) yvec
!     3) nele_cluster_mag and nele_cluster electron numbers in cluster and env 
!     4) fermi - fermi levels 
!     5) exc_exx - XC energy for j_atom 
!
subroutine dfet_zp_envOF(logOF,vw_lam,vw_reg,mxband,nfft,nspin,natom,j_atom, & 
                   ucvol,logf,myrank,do_exx,do_xcep, & 
                   zp_alpha,nopt_uemb,global_fermi,cell_nfft,cgW_tol,qvec,zp_coeff,zp_conv,& 
                   ref_rho,global_vks,global_vOF,prj_vks,q_clu,q_env, & 
                   env_weight,tsmear,sub_rhor,fermi_sub,& 
                   natom_clu,natom_env,nele_cluster_mag,nele_cluster, & 
                   exc_exx,cluster_exc_eps,u,den_diff,yvec,& 
                   drho_dN_clu,drho_dN_env,deps_dN_clu)
   
   implicit none 

   integer :: logOF, nfft, nspin, logf, natom, j_atom, & 
              do_exx, do_xcep, myrank, & 
              global_fermi, nopt_uemb,&  
              mxband, cell_nfft(3), & 
              natom_clu, & 
              natom_env

   real(8) :: vw_reg,vw_lam, & 
              ucvol, cgW_tol, zp_alpha, & 
              global_vks(nfft,nspin), & 
              global_vOF(nfft,nspin), & 
              prj_vks(nfft,nspin,2), & 
              q_clu(nspin), & 
              q_env(nspin),& 
              env_weight(nfft), & 
              deps_dN_clu(nfft,nspin), &  
              ref_rho(nfft,nspin), & 
              sum_rho(nfft,nspin), & 
              zp_coeff, zp_conv, tsmear, & 
              drho_dN_clu(nfft,nspin), & 
              drho_dN_env(nfft,nspin), & 
              fermi_sub(nspin,natom,2), & 
              qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3))

   real(8), intent(out) :: u(nfft,nspin), & 
                           yvec(nfft,nspin), & 
                           sub_rhor(nfft,nspin,2), &    ! the last index is for cluster and env 
                           nele_cluster_mag(natom,2), & 
                           nele_cluster(natom,2), & 
                           exc_exx, cluster_exc_eps(nfft,nspin) 


   ! local vars ================
   logical :: last_iter, do_precond
   real(8) :: start_time, end_time, & 
              start_OFtime, & 
              end_OFtime
   integer :: isp, s, nnew,  & 
              n1,n2,n3, ii,& 
              iter, precond_type,  &
              local_do_exx, & 
              load_dEdVks, & 
              load_new_occ, &
              local_do_xcep, & 
              nband(2), & 
              dim1,ix,iy,iz

   real(8) :: nelectr(nspin), & 
              factor, & 
              vOF_env(nfft,nspin), & 
              qq,coeff, OF_chempot, & 
              comm_fermi,  & 
              ne_clu(nspin), ne_env(nspin), & 
              ee_clu_old, & 
              ee_clu_old2, & 
              ee_env_old, & 
              ee_env_old2, & 
              max_coeff,& 
              u_new(nfft,nspin), &  
              env_vOF(nfft), & 
              dtmpv1, dtmpv2, & 
              gnorm, dvol,& 
              total_entropy, & 
              tmp_entropy, & 
              band_energy, & 
              eigenvalues(mxband,nspin,2), & 
              occ(mxband,nspin,2), & 
              vh_tot(nfft), & 
              rho_diff(nfft), & 
              den_diff_old(nfft,nspin), & 
              den_diff(nfft,nspin), & 
              tmp_vec(nfft,nspin), & 
              tmp_vec2(nfft,nspin), & 
              tmp_vec3(nfft), & 
              cc = 1e-4, ee_clu, ee_env, e_coul, & 
              d_rho,dtmp, &  
              u_best(nfft,nspin), & 
              sub_ke(natom), & 
              sub_ts(natom),sub_nlps(natom), & 
              sub_etotal(natom), & 
              vec1(nfft), vec2(nfft), & 
              accu(nspin), pen(nspin), lap(nfft), & 
              q3d(3,cell_nfft(1)/2+1,cell_nfft(2),cell_nfft(3))

   complex(kind=8) :: fft1((cell_nfft(1)/2+1),cell_nfft(2),cell_nfft(3))
   

   ! pulay 
   integer, parameter :: npulay = 5
   real(8) :: pulay_beta = 0.01d0
   integer :: pulay_iwork(5,nspin)
   real(8) :: pulay_work(nfft,npulay,2,nspin)
   real(8) :: resid(nfft,npulay,nspin)


   ! >>>>>>>>>> function begins <<<<<<<<<<<<<

   do_precond = .false.
   precond_type = 3  ! 1: precond by dielectric matrix
                     ! 2: density based precondition (very important)
                     ! 3: kerker method (focus on long wavelength)
      

   write(715,'(a)')NEW_LINE('a')//'>>> enter dfet_zp() <<<'
   if (natom_env==0) write(715, '(a)')'NOTE: natom_env=0'
   call flush(715)

   pulay_beta = 0.1d0


   !========================
   ! information 
   !========================
   if(myrank==0) then 
     write(logf,'(a)')NEW_LINE('a')//'   >>> enter dfet_zp() <<< '
     write(logf,*)'*** potential mixing ***'
     write(logf,'(a,L2)')'do_precond: ',do_precond
     if (do_precond) then 
       if (precond_type==1) write(logf,'(a)')'** precondition is on. residual is preconditioned **'
       if (precond_type==2) write(logf,'(a)')'** precondition is on. density-based metric **'
       if (precond_type==3) write(logf,'(a)')'** precondition is on. kerker mixing  **'
     endif 
     write(logf,'(a,f12.4)') 'zp_coeff:   ',zp_coeff
     write(logf,'(a,es12.4,a,f10.4,a)')'zp_alpha: ',zp_alpha, & 
       '   1/zp_alpha=',1/zp_alpha,' bohr (screen length)' 
     write(logf,'(a,es12.4)')'zp_conv:    ',zp_conv
     write(logf,'(a,i4)')    'nopt_uemb:  ',nopt_uemb
     write(logf,'(a,f8.4)')  'pulay_beta: ',pulay_beta
     write(logf,'(a,i3)')    'npulay:     ',npulay

     if (global_fermi==2) then 
       write(logf,'(a)')'not fully implemented'
       stop
     else if (global_fermi==1) then 
       write(logf,'(a)')'charge in cluster and env are determined via partitioning'
     endif 
     call flush(logf)
   endif 


   !=================
   ! variables 
   !=================
   n1 = cell_nfft(1)
   n2 = cell_nfft(2)
   n3 = cell_nfft(3)
   
   dvol = ucvol / dble(nfft) 
   nelectr = sum(ref_rho,1)*dvol

   if (myrank==0) & 
     call print_monitor_electron_numbers(logf,nspin,nelectr)

   last_iter = .false.
   nele_cluster = 0.d0 
   nele_cluster_mag = 0.d0
   iter = 0
   call print_header()

   call cpu_time(start_time) 


   ! for OF-DFT, initialize the environment's density 
   write(logOF,'(a,f12.6)')'get in dfet_zp.f90, q_env:',q_env
   do isp=1,nspin
     sub_rhor(:,isp,2) = max(1e-8,ref_rho(:,isp)*env_weight)  ! make sure sub_rhor is positive
     sub_rhor(:,isp,2) = sub_rhor(:,isp,2)/sum(dvol*sub_rhor(:,isp,2))*q_env(isp) ! normalize
   enddo


   !=====================================
   !  Zhao-Parr method solving for vemb
   !=====================================
   do while (.true.) 

     iter = iter + 1

     ! restore the best u to compute the EXX related quantities 
     if (last_iter) u = u_best

     call dump_uemb_xcpatch(j_atom,u,nfft,nspin)

     local_do_exx = -1
     load_new_occ = -1
     load_dEdVks  = -1
     local_do_xcep = -1

     if (last_iter .and. do_exx>0)  then 
       local_do_exx = 1
     endif 

     if (last_iter .and. do_xcep>0) then 
       load_dEdVks = 1
       local_do_xcep = 1
     endif

     call run_subsystem_xcpatch(j_atom,local_do_exx,load_new_occ,local_do_xcep,last_iter)

     !=========================================================
     ! Obtain energies, eigenvalues, densities for subsystems
     !=========================================================
     do s=1,2

       ! skip environment if no atom is in env
       ! set env density to zero 
       if (natom_env==0 .and. s==2) then 
         sub_rhor(:,:,2) = 0.0d0
         fermi_sub(:,j_atom,2)=0.d0
         exit
       endif 

       ! treat env with OF-DFT
       if (s==2) then 
         do isp=1,nspin 
           env_vOF = u(:,isp)+prj_vks(:,isp,2)
           call cpu_time(start_OFtime)
           call OF_cgmin_tot(n1,n2,n3,nfft,nspin,qvec,ucvol,env_vOF, & 
                             sub_rhor(:,isp,2),fermi_sub(isp,j_atom,2))
           call cpu_time(end_OFtime)
         enddo
         exit 
       endif 

       ! for s=1, doing subsys folder
       ! for s=2, doing subsys_env folder
       call wait_for_system_end_xcpatch(j_atom,s)

       local_do_exx = -1
       if (s==1 .and. last_iter .and. do_exx>0) local_do_exx = 1

       if (s==1) then 
         ! cluster 
         call get_subsystem_data_xcpatch(j_atom,s,nfft,dvol,nspin,sub_etotal(s), & 
          sub_nlps(s),sub_TS(s),sub_ke(s),sub_rhor(:,:,s),fermi_sub(:,j_atom,1), & 
          mxband,eigenvalues(:,:,s),occ(:,:,s),nband(s), & 
          local_do_exx,exc_exx,cluster_exc_eps,yvec,load_dEdVks,drho_dN_clu,deps_dN_clu)        
       else
         ! env
         call get_subsystem_data_xcpatch(j_atom,s,nfft,dvol,nspin,sub_etotal(s), & 
          sub_nlps(s),sub_TS(s),sub_ke(s),sub_rhor(:,:,s),fermi_sub(:,j_atom,2), & 
          mxband,eigenvalues(:,:,s),occ(:,:,s),nband(s),& 
          local_do_exx,dtmp,tmp_vec,tmp_vec2,load_dEdVks,drho_dN_env,& 
          deps_dN_clu) ! deps_dN_clu is not updated for env
       endif
     enddo

     sum_rho = sum(sub_rhor,3)
     d_rho = maxval(abs(ref_rho-sum(sub_rhor,3)))

     
     ! last iteration 
     !=========================
     if (last_iter) then 
       write(715,*) 'exc_exx: ',exc_exx
       if (myrank==0) then 
         write(logf,'(a)')'optimization of v_emb on process 0 is done.'
         call flush(logf)
       endif 
       exit 
     endif

     nele_cluster(j_atom,1) = dvol*sum(sub_rhor(:,:,1)) ! cluster's electron 
     nele_cluster(j_atom,2) = dvol*sum(sub_rhor(:,:,2)) ! env's electron 

     if (nspin==2) then 
       nele_cluster_mag(j_atom,1) = dvol*(sum(sub_rhor(:,1,1))-sum(sub_rhor(:,2,1)))  ! magnetic mom
       nele_cluster_mag(j_atom,2) = dvol*(sum(sub_rhor(:,1,2))-sum(sub_rhor(:,2,2)))  ! magnetic mom
     endif 


     ! measure the mismatch 
     !==========================
     e_coul = 0.d0 
     ee_clu = 0.0d0
     ee_env = 0.d0
     do s=1,nspin
       call hartree(n1,n2,n3,1,ucvol,qvec,sub_rhor(:,s,1),tmp_vec3,dtmp); ee_clu = ee_clu + dtmp
       call hartree(n1,n2,n3,1,ucvol,qvec,sub_rhor(:,s,2),tmp_vec3,dtmp); ee_env = ee_env + dtmp
       call hartree(n1,n2,n3,1,ucvol,qvec,sub_rhor(:,s,1)+sub_rhor(:,s,2)-ref_rho(:,s),vh_tot,dtmp)
       e_coul = e_coul + dtmp
       accu(s) = dtmp
     enddo 



     ! potential mixing 
     !================================
     ! make new embedding potential (Zhao-Parr method)
     ! Yukawa potential for penalty
     do s=1,nspin
       rho_diff = sum_rho(:,s)-ref_rho(:,s)
       call yukawa(n1,n2,n3,1,ucvol,qvec,zp_alpha,rho_diff,vh_tot)
       !!call hartree(n1,n2,n3,1,ucvol,qvec,rho_diff,vh_tot,dtmp)
       if (.not. do_precond) then 
         ! new formula to update embedding potential 
         ! better convergence rate, safe to use with large pulay_beta
         factor = 1.0d0
         !factor = 1.d0
         u_new(:,s) = factor*vh_tot - factor*u(:,s)/zp_coeff + u(:,s)
       else
         if (precond_type==1) u_new(:,s) = zp_coeff*vh_tot
         if (precond_type==2) u_new(:,s) = vh_tot - u(:,s)/zp_coeff + u(:,s)
         if (precond_type==3) u_new(:,s) = vh_tot - u(:,s)/zp_coeff + u(:,s)
       endif 
     enddo 

     ! compute penalty term
     do s=1,nspin
         call laplacian(nfft,n1,n2,n3,u_new(:,s),qvec,lap)
         pen(s) = - sum(lap*u_new(:,s))*dvol/8.d0/3.1415926d0/zp_coeff  & 
                  + zp_alpha**2/8.d0/3.1415926d0/zp_coeff*sum(u_new(:,s)**2)*dvol
     enddo

     ! some information (potential mixing)
     dtmpv1 = sum((u(:,1)-u_new(:,1))**2)*ucvol/nfft
     if (nspin==2) then 
       dtmpv2 = sum((u(:,2)-u_new(:,2))**2)*ucvol/nfft
     else
       dtmpv2 = 0.d0 
     endif 
     call cpu_time(end_time)
     call output_iter_info()
     call cpu_time(start_time)


     ! check convergence 
     !======================
     if (nopt_uemb<0 .and. iter>=20 .and. & 
       abs(ee_clu_old -ee_clu)<zp_conv .and. & 
       abs(ee_clu_old2-ee_clu)<zp_conv .and. & 
       abs(ee_env_old -ee_env)<zp_conv .and. & 
       abs(ee_env_old2-ee_env)<zp_conv ) then 
       write(logf,'(a,es8.2,a)')'cluster and env hartree energies converged to ',zp_conv,' Ha.'
       write(715, '(a,es8.2,a)')'cluster and env hartree energies converged to ',zp_conv,' Ha.'
       u_best = u
       last_iter = .true. 
       cycle 
     endif 
     ee_clu_old2 = ee_clu_old
     ee_clu_old  = ee_clu
     ee_env_old2 = ee_env_old 
     ee_env_old  = ee_env 


     if (nopt_uemb>0 .and. iter==nopt_uemb) then 
       write(logf,'(a,i4)')'reach max iteration number => ',nopt_uemb
       u_best = u
       last_iter = .true.
       cycle  
     endif 

     !if (natom_env==0) then 
     !  write(715,'(a)')'natom_env=0, finish dfet.'
     !  write(logf,'(a)')'natom_env=0, finish dfet.'
     !  u_best = u
     !  last_iter = .true.
     !endif 

     ! potential mixing 
     !=======================
     do isp=1,nspin
       if (.not. do_precond) then 
         call pulay_mix(nfft,iter,u(:,isp),u_new(:,isp),pulay_work(:,:,:,isp),npulay,pulay_beta)
       else
         ! precondition 
         select case (precond_type)
         case(1) 
           call pulay_mix_dfet_metric(n1,n2,n3,nfft,iter,u(:,isp),u_new(:,isp), & 
             pulay_work(:,:,:,isp),npulay,pulay_beta,zp_alpha,zp_coeff,qvec,ucvol)
         case(2) 
           call pulay_mix_den_metric(nfft,iter,ref_rho(:,isp),u(:,isp),u_new(:,isp), & 
             pulay_work(:,:,:,isp),npulay,pulay_beta)
         case(3)
           call pulay_mix_kerker_dfet(n1,n2,n3,nfft,iter,u(:,isp),u_new(:,isp), & 
             pulay_work(:,:,:,isp),npulay,pulay_beta,qvec,ucvol)
         endselect
       endif
     enddo
     u = u_new  


   enddo
   !============ done ZMP method ===========


   if (nspin==1) then 
     write(715,'(a,2es12.4)')'yvec: ',minval(yvec),maxval(yvec)
   else
     write(715,'(a,2es12.4)')'yvec (spin-up): ',minval(yvec(:,1)),maxval(yvec(:,1))
     write(715,'(a,2es12.4)')'     (spin-dn): ',minval(yvec(:,2)),maxval(yvec(:,2))
   endif 

   if(myrank==0) then 
     write(logf,'(a)')'done dfet_zp()'
     call flush(logf)
   endif 
   write(715,'(a)')'done dfet_zp()'
   call flush(715)






contains 

subroutine  print_header 
   write(715,'(a)')'*** env is treated by OF-DFT ***'
   write(logf,*)'vw_lam: ',vw_lam
   write(logf,'(a)')'*** env is treated by OF-DFT ***'


   ! potential mixing
   !=======================
     if (nspin==1) then 
         ! env is treated by OF 
         write(715,'(2a)')'iter max(d_rho)  |d_u|     min(u)    max(u)    e_coul      J_clu      J_env   time_tot time_OF'
     else
       ! spin polarized 
       write(715,'(a,a)')'iter   max(d_rho)  |d_u(up)|   |d_u(dn)| ', & 
       '  min/max(u(up))          min/max(u(dn))          e_coul     time(s) '
     endif 
     call flush(715)
     if (myrank==0) then 
      if (nspin==1) then 
          ! env is treated by OF 
          write(logf,'(a)')'iter max(d_rho)  |d_u|     min(u)    max(u)    e_coul      J_clu      J_env   time_tot time_OF'
      else
        ! spin polarized 
        write(logf,'(a,a)')'iter   max(d_rho)  |d_u(up)|   |d_u(dn)| ', & 
        '  min/max(u(up))          min/max(u(dn))          e_coul        J_clu       J_env   time(s)'
      endif 
      call flush(logf)
     endif
end subroutine print_header



subroutine output_iter_info_denmix
     ! output iteration info
     !
     ! xcep.out file
     !
     if (myrank==0) then
       if (nspin==2) then 
         write(logf,'(i3,es12.2,2es12.2,5es12.2,f9.1)')iter,d_rho,dtmpv1,dtmpv2,& 
          minval(u(:,1)),maxval(u(:,1)),minval(u(:,2)),maxval(u(:,2)),e_coul,end_time-start_time
       else 
         write(logf,'(i3,es12.2,a,es12.2,a,3es12.2,2f12.6,f9.1,f9.1)') & 
          iter,d_rho,' ',dtmpv1,' ',minval(u),maxval(u),e_coul,ee_clu,ee_env, & 
          end_time-start_time,end_OFtime-start_OFtime
       endif 
     endif 
     !
     ! unit=715
     !
     if (nspin==2) then 
       write(715,'(i3,es12.2,2es12.2,5es12.2,f9.1)')iter,d_rho,dtmpv1,dtmpv2,& 
         minval(u(:,1)),maxval(u(:,1)),minval(u(:,2)),maxval(u(:,2)),e_coul,end_time-start_time
     else 
       write(715,'(i3,es12.2,a,es12.2,a,3es12.2,2f12.6,f9.1,f9.1)') & 
        iter,d_rho,' ',dtmpv1,' ',minval(u),maxval(u),e_coul,ee_clu,ee_env, & 
        end_time-start_time,end_OFtime-start_OFtime
     endif 
     call flush(logf)
     call flush(715)
end subroutine output_iter_info_denmix




subroutine output_iter_info
     ! output iteration info
     if (myrank==0) then
       if (nspin==2) then 
         write(logf,'(i3,es12.2,2es11.2,5es11.2,f9.1)')iter,d_rho,dtmpv1,dtmpv2,& 
          minval(u(:,1)),maxval(u(:,1)),minval(u(:,2)),maxval(u(:,2)),e_coul,end_time-start_time
       else 
           write(logf,'(i3,es10.2,a,es9.2,a,3es10.2,2f11.5,f8.1,f8.1)') & 
            iter,d_rho,' ',dtmpv1,' ',minval(u),maxval(u),e_coul,ee_clu,ee_env, & 
            end_time-start_time,end_OFtime-start_OFtime
       endif 
     endif 
     !
     ! unit=715
     !
     if (nspin==2) then 
       write(715,'(i3,es12.2,2es12.2,5es12.2,f9.1)')iter,d_rho,dtmpv1,dtmpv2,& 
         minval(u(:,1)),maxval(u(:,1)),minval(u(:,2)),maxval(u(:,2)),e_coul,end_time-start_time
     else 
       write(715,'(i3,es12.2,a,es12.2,a,3es12.2,2f12.6,f9.1,f9.1,2f8.3)') & 
        iter,d_rho,' ',dtmpv1,' ',minval(u),maxval(u),e_coul,ee_clu,ee_env,& 
        end_time-start_time,end_OFtime-start_OFtime,sum(ne_clu),sum(ne_env)
     endif 
     call flush(logf)
     call flush(715)
end subroutine output_iter_info



!
! Yukawa potential defined as exp(-alpha*r)/r
!
subroutine yukawa(n1,n2,n3,nspin,ucvol,qvec,alpha,rhor,vhart)
  implicit none 

  integer :: n1,n2,n3,ix,iy,iz,dim1,kk,isp,nspin

  real(kind=8) :: rhor(n1*n2*n3,nspin), &   ! spin density
                  vhart(n1*n2*n3), &        ! hartree potential due to rho_up + rho_dn
                  qvec(3,(n1/2+1)*n2*n3), &
                  alpha, & 
                  rtmp(n1,n2,n3), &         ! temporary array in real space
                  ucvol,        &           ! volume of cell
                  q3d(3,n1/2+1,n2,n3), &
                  pi=3.14159265359d0,  & 
                  qq, &                     ! q^2
                  ehart                     ! hartree energy
  
  complex(kind=8) :: fft1(n1/2+1,n2,n3), & 
                     fft2(n1/2+1,n2,n3)

  dim1=(n1/2+1)

  call FFT(n1,n2,n3,reshape(sum(rhor,2),(/n1,n2,n3/)),fft1,1)

  ! convert qvec to 3-dimension 
  q3d(1,:,:,:) = reshape(qvec(1,:),(/dim1,n2,n3/))
  q3d(2,:,:,:) = reshape(qvec(2,:),(/dim1,n2,n3/))
  q3d(3,:,:,:) = reshape(qvec(3,:),(/dim1,n2,n3/))

  vhart = 0.0d0

  do ix=1,dim1
    do iy=1,n2
      do iz=1,n3
        qq = q3d(1,ix,iy,iz)**2 + & 
             q3d(2,ix,iy,iz)**2 + & 
             q3d(3,ix,iy,iz)**2

        fft2(ix,iy,iz) = fft1(ix,iy,iz)*4.0d0*pi/(qq+alpha**2)
      enddo
    enddo
  enddo

  call FFT(n1,n2,n3,rtmp,fft2,-1)
  vhart = reshape(rtmp,(/n1*n2*n3/))

end subroutine yukawa

end subroutine dfet_zp_envOF



