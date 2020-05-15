module interface_funcs

 implicit none 

 integer, parameter  :: max_psp_file = 20,  & 
                        pspfile_len=100

 interface 


   subroutine solve_z(natom,natom_clu,natom_env,nfft,nspin, & 
                      maxiter,j_atom,myrank,logf,cell_nfft,ucvol,qvec,yvec, & 
                      vks_clu,vks_env,sub_rhor,zp_coeff,reg_z,zvec_init,zvec)
   implicit none 
   integer :: nfft, nspin, j_atom, maxiter, & 
              myrank, logf, cell_nfft(3), zvec_init, & 
              natom, natom_clu, natom_env

   real(8) :: yvec(nfft,nspin),  & 
              sub_rhor(nfft,nspin,2), & 
              vks_clu(nfft,nspin), & 
              qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), & 
              vks_env(nfft,nspin), & 
              zvec(nfft,nspin), & 
              ucvol, reg_z, & 
              zp_coeff
   end subroutine 



   subroutine OF_dfpt(nfft,cell_nfft,ucvol,qvec, & 
                   rho,v_perb,delta_rho,delta_mu)
    integer :: nfft, cell_nfft(3)
    real(8) :: delta_rho(nfft), &     
               delta_mu, &        
               chempot, rho(nfft), & 
               v_perb(nfft), & 
               qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), & 
               ucvol 
   end subroutine OF_dfpt 



   subroutine get_subsystem_data_xcpatch(j_atom,sys,nfft,dvol, & 
                                         nspin,etotal,nlpsp,TS,ke,rhor,fermi, & 
                                         mxband,eigenvalues,occ,nband, & 
                                         do_exx,exc_exx,exc_eps,yvec,load_dEdVks,drho_dN,deps_dN)
    implicit none
    integer      :: sys,nfft,nspin,j_atom,do_exx
    integer      :: mxband,nband,load_dEdVks
    real(kind=8) :: eigenvalues(mxband,nspin), &   ! eigenvalues from the ABINIT program 
                    occ(mxband,nspin)              ! occupation number of KS oribtals 
    real(kind=8) :: etotal, dvol, &      ! etotal: total energies for each subsystem 
                    TS, nlpsp,ke, exc, ehart, elocal, & !
                    yvec(nfft,nspin), & 
                    drho_dN(nfft,nspin), &
                    deps_dN(nfft,nspin), & 
                    fermi(nspin), & 
                    exc_exx, exc_eps(nfft),  &
                    rhor(nfft,nspin)     ! subsystem electron density 
   end subroutine 

   subroutine solve_total_vxc(iter_scf,nfft,cell_nfft,qvec,nspin,reg_type, & 
       logf,ucvol,global_vks,reg_lambda,init,dExc_dvks,vxc)
    use MPI
    implicit none
    integer :: iter_scf,nfft,reg_type,nspin, logf, init, cell_nfft(3)
    real(8) :: dExc_dvks(nfft,nspin), & 
               vxc(nfft,nspin),ucvol, & 
               qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), &
               global_vks(nfft,nspin), & 
               reg_lambda
   end subroutine 

   subroutine solve_total_vxc_cgls(nfft,cell_nfft,nspin,logf,ucvol,qvec, & 
                         global_vks,rho,reg_lambda,init,dExc_dvks,vxc)

    use mpi
    implicit none
    integer :: nfft, cell_nfft(3), nspin, logf, init
    real(8) :: dExc_dvks(nfft,nspin), & 
               reg_lambda, & 
               vxc(nfft,nspin),ucvol, & 
               qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), &
               global_vks(nfft,nspin), & 
               rho(nfft,nspin)
   endsubroutine

   subroutine apply_chi(info,system_type,j_atom,nfft,nspin,vks,vperb,rho_perb)
    implicit none 
    integer :: system_type, &   ! 0: total system 
                              ! 1: cluster
                              ! 2: env 
               j_atom, nfft, nspin, info, & 
               type_of_system
    real(8) :: vks(nfft,nspin), & 
               vperb(nfft,nspin) , &   ! perturbing potnetial 
               rho_perb(nfft,nspin)
   end subroutine 


   subroutine calc_subsystem_dfpt(natom,natom_clu,natom_env,j_atom,nfft,nspin,ucvol, & 
                                  cluster_weight,env_weight,qvec,sub_rhor,cell_nfft, & 
                                  dfermi_dvks,vks_cluster,vks_env,vperb,rho_perb)
    implicit none 
    integer :: j_atom, nfft, nspin, natom, natom_clu, natom_env, cell_nfft(3)
    real(8) :: ucvol, & 
               sub_rhor(nfft,nspin,2) , & 
               qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), & 
               dfermi_dvks(nfft,nspin,2), & 
               vks_cluster(nfft,nspin), & 
               vks_env(nfft,nspin), & 
               vperb(nfft,nspin) , &   ! perturbing potnetial 
               rho_clu(nfft,nspin), &  ! perturbed density of cluster 
               rho_perb(nfft,nspin), & 
               rho_env(nfft,nspin), & 
               cluster_weight(nfft), & 
               env_weight(nfft) 
   endsubroutine 

   subroutine calc_subsystem_rho_xcep(j_atom,nfft,nspin,q_cluster,q_env, & 
                     mag_cluster,mag_env,vks_cluster,vks_env,rho_tot)
    implicit none 
    integer :: j_atom, nfft, nspin
    real(8) :: vks_cluster(nfft,nspin), & 
               vks_env(nfft,nspin), & 
               rho_clu(nfft,nspin), & 
               rho_tot(nfft,nspin), & 
               rho_env(nfft,nspin), & 
               q_cluster, q_env, & 
               mag_cluster, mag_env 
   end subroutine 



   subroutine grand_chempot(nsys,nspin,q_alpha,q_beta,eigen,nband, & 
                      nconfig_max,nconfig,config,grand_temp,sub_etotal,weight,fermi)
    implicit none
    integer :: nsys, nspin, nband(nsys), nconfig_max, nconfig(nsys)
    real(8) :: q_alpha, q_beta, eigen(nconfig_max,2,100,nsys), & 
               config(nspin,nconfig_max,nsys), grand_temp, &
               sub_etotal(nconfig_max,nsys), & 
               weight(nconfig_max,nsys), & 
               fermi(nspin) 
   endsubroutine

   subroutine make_weights(file_dfet_out,natom,nspin,nfft,nsys,melement,atom_info, & 
                           cell_nfft,xcart,cell_acell,znucl,rprimd,cov_radius,watom)
     implicit none 
     integer :: nfft, nspin, nsys, natom, & 
                melement, file_dfet_out, cell_nfft(3)
     real(8) :: cov_radius(melement), & 
                watom(nfft,nsys), & 
                xcart(3,natom), & 
                cell_acell(3,3), & 
                znucl(natom) , & 
                rprimd(3,3), & 
                atom_info(5,natom)
   endsubroutine make_weights

   subroutine load_parameters(restart,file_dfet_out,natom,nspin,npsps,cell_nfft,nsys,global_fermi,mix_gvks, & 
                           start_u,dft_type,nopt_uemb,npulay_gvks,pulay_gvks_beta,cgW_tol, & 
                           cell_acell,ucvol,rprimd,rmet,tsmear,relax_geo,only_lda, & 
                           reg_type,reg_x,reg_vxc,zp_coeff,zp_conv,vac,plot_sub_rho,relax_nscf, & 
                           exx_method,hse_range)

    implicit none 
 
     ! external vars
     integer,parameter    :: dp=8
     logical              :: plot_sub_rho,only_lda
     integer              :: reg_type,restart,nopt_uemb,global_fermi,mix_gvks,relax_geo
     integer              :: file_dfet_out,dft_type,relax_nscf,exx_method
     integer,intent(out)  :: natom,nspin,cell_nfft(3),nsys,start_u,npsps,npulay_gvks
     real(dp),intent(out) :: ucvol,cell_acell(3,3),rprimd(3,3),rmet(3,3), cgW_tol, & 
                             pulay_gvks_beta,tsmear,hse_range, & 
                             reg_x,reg_vxc,zp_coeff,zp_conv,vac
   end subroutine 


subroutine run_global_ks(n1,n2,n3,nfft,nspin,natom,qvec,ucvol, & 
                         do_exx,do_xcep,global_vks,ref_rho,ehart, & 
                         ke,nlpsp,TS,elocal,ewald,corepsp,total_fermi,dfermi_dvks,band_energy, & 
                         exc_eps,exc_exx,dExc_dvks,force_hf,vhart)

    implicit none 

    ! external vars 
    integer :: nfft,nspin,n1,n2,n3,do_exx,do_xcep,natom
    real(8) :: ref_rho(nfft,nspin),etotal,ucvol, &
               dExc_dvks(nfft,nspin), total_fermi(nspin), & 
               dfermi_dvks(nfft,nspin), & 
               exc_eps(nfft), &  
               band_energy, &  
               global_vks(nfft,nspin),qvec(3,(n1/2+1)*n2*n3), & 
               ehart,ke,nlpsp,TS,elocal,ewald,corepsp, & 
               exc_exx, & 
               force_hf(3,natom), & 
               vhart(nfft) 

end subroutine 


subroutine load_atom_info(natom,xcart,znucl,atom_info,mlist,nlist,ilist,nshare,ishare)

  implicit none 

  character(len=100) :: string
  integer :: natom,ii,file_unit=200,iost,int_znucl(natom), & 
             jj,mlist,nlist(natom),ilist(mlist,natom), & 
             nshare(natom), ishare(natom,mlist)
  real(8) :: xcart(3,natom), znucl(natom), & 
             atom_info(5,natom)
end subroutine             

 
subroutine init_conv_radius( melement, cov_radius)
    implicit none 
    integer :: melement
    real(8) :: cov_radius(melement)
end subroutine     



subroutine init_uemb(nspin,nfft,start_u,u)

 implicit none 

 integer :: nspin, nfft, start_u, isp, s
 real(8) :: u(nfft,nspin)
end subroutine  



   subroutine atom2cell(s,celldata,cellnfft,cellacell,ucvol,atbox_acell,atbox_npt,xatom,atomdata)
     implicit none
     integer :: s,atbox_npt,cellnfft(3)
     real(8) :: celldata(cellnfft(1),cellnfft(2),cellnfft(3))
     real(8) :: cellacell(3),ucvol,atbox_acell,xatom(3),atomdata(atbox_npt,atbox_npt,atbox_npt)
   endsubroutine 

   subroutine cell2atom(s,cell_data,cell_nfft,cell_acell, & 
                        atbox_acell,atbox_npt,xatom,atom_data) 
     implicit none 
     integer :: s,cell_nfft(3),atbox_npt
     real(8) :: cell_data(cell_nfft(1),cell_nfft(2),cell_nfft(3)), & 
                cell_acell(3), atbox_acell, xatom(3), & 
                atom_data(atbox_npt,atbox_npt,atbox_npt)
   end subroutine 

   subroutine check_bfgs_flag(task, new_x, nfft, nspin, u, iter_umin, etot)
     character(len=60) :: task
     logical           :: new_x
     integer           :: iter_umin, nfft, nspin
     real(kind=8)      :: etot, u(nfft,nspin)
   end subroutine 


   subroutine make_gprim(cell_acell,gprim,gmet,abinit_gmet)
     implicit none
     real(kind=8) :: cell_acell(3)
     real(kind=8) :: gprim(3,3),gmet(3,3),abinit_gmet(3,3)
   end subroutine

   subroutine make_q_vector(n1,n2,n3,gprim,qvec)
     implicit none
     integer,intent(in) :: n1,n2,n3
     real(kind=8),intent(out) :: qvec(3,(n1/2+1)*n2*n3)
     real(kind=8),intent(in)  :: gprim(3,3)
   end subroutine make_q_vector

   subroutine optimize_charge(iter_enum,nsys,nfft,u,sub_nele)
     implicit none
     integer :: iter_enum, nsys, nfft
     real(8) :: u(nfft), sub_nele(nsys), nele0(nsys)
   end subroutine optimize_charge

   subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, & 
                     task, iprint,  csave, lsave, isave, dsave, bfgs_time, myrank)
     implicit none 
     integer, intent(in) :: n, m, myrank
     character*60     :: task, csave
     logical          :: lsave(4)
     integer          :: iprint, nbd(n), iwa(3*n), isave(44)
     double precision :: f, factr, pgtol, x(n), l(n), u(n), g(n), wa(2*m*n+4*n+12*m*m+12*m), dsave(29)
     real :: bfgs_time
   end subroutine setulb

   subroutine wait_for_system(s)
     implicit none
     integer :: s
   end subroutine 


   subroutine hartree(n1,n2,n3,nspin,ucvol,qvec,rhor,vhart,ehart)
     implicit none
     integer :: n1,n2,n3,nspin
     real(kind=8) :: ucvol,qvec(3,(n1/2+1)*n2*n3),rhor(n1*n2*n3,nspin), & 
                     vhart(n1*n2*n3),ehart
   end subroutine 

   subroutine obj_gradient(iter,ike,ixc,nsys,n1,n2,n3,nfft,nspin,nconfig_max,nconfig,config, & 
                           grand_temp,dvol,xatom,atbox_acell,atbox_npt,cell_acell,epsatm, & 
                           lambda,int_ewald,int_pspcore,rhor,etotal,weight,qvec, & 
                           sub_vloc,sub_ke,u,int_energy,tot_nlpsp,tot_entropy,grad)
     implicit none
     integer,intent(in)       :: iter,ixc,ike,nsys,n1,n2,n3,nfft,nspin, & 
                                 nconfig_max,nconfig(nsys),atbox_npt
     real(kind=8),intent(in)  :: dvol, weight(nconfig_max,nsys), & 
                                 epsatm(nsys), sub_ke(nsys), &
                                 rhor(nfft,nspin,nsys), & 
                                 config(2,nconfig_max,nsys), & 
                                 xatom(3,nsys), atbox_acell, cell_acell(3), & 
                                 etotal(nconfig_max,nsys), &
                                 qvec(3,(n1/2+1)*n2*n3), & 
                                 u(nfft,nspin),lambda(nsys), &
                                 grand_temp,int_ewald,int_pspcore, & 
                                 sub_vloc(nfft,nsys)

     real(kind=8),intent(out) :: grad(nfft,nspin), & 
                                 int_energy,  & 
                                 tot_nlpsp, & 
                                 tot_entropy
   end subroutine obj_gradient

   subroutine total_energy(nsys,nspin,ike,sub_etotal,sub_ewald,sub_pspcore,& 
                           sub_nlpsp,sub_entropy,int_energy, &
                           tot_ewald,tot_pspcore,tot_nlpsp,tot_entropy,etotal)
     implicit none 
     integer :: nsys, nspin, ike
     real(kind=8) :: sub_etotal(nsys), & 
                     sub_ewald(nsys),  &
                     sub_pspcore(nsys), &
                     sub_nlpsp(nspin,nsys),   &
                     sub_entropy(nsys), & 
                     int_energy,  & 
                     tot_ewald, & 
                     tot_pspcore,  & 
                     tot_nlpsp,  & 
                     tot_entropy, & 
                     etotal
   end subroutine total_energy

   subroutine calculate_kedf(isys,n1,n2,n3,nspin,ike,rhor,qvec,dvol,ts_pot,ts_energy)
     implicit none
     integer      :: n1,n2,n3,nspin,ike,isys
     real(kind=8) :: rhor(n1*n2*n3,nspin), &
                     qvec(3,(n1/2+1)*n2*n3), &
                     dvol, &
                     ts_pot(n1*n2*n3,nspin), &
                     ts_energy
   end subroutine calculate_kedf

   subroutine calculate_xc_pbe(n1,n2,n3,nspin,rho,qvec,dvol,xc_pot,xc_energy)
     implicit none
     integer      :: n1,n2,n3,nspin
     real(kind=8) :: dvol
     real(kind=8) :: rho(n1,n2,n3,nspin)
     real(kind=8) :: xc_pot(n1,n2,n3,nspin)
     real(kind=8) :: qvec(3,(n1/2+1)*n2*n3)
     real(kind=8) :: xc_energy
   end subroutine calculate_xc_pbe 

   subroutine  calculate_xc_pw92LSDA(n1,n2,n3,nspin,rho,dvol,xc_pot,xc_energy)
     implicit none 
     integer                 :: n1,n2,n3,ixc,nspin
     real(kind=8)            :: rho(n1*n2*n3,nspin), & 
                                dvol
     real(kind=8),intent(out) :: xc_pot(n1*n2*n3,nspin), & 
                                 xc_energy
   end subroutine

 end interface 

end module interface_funcs
