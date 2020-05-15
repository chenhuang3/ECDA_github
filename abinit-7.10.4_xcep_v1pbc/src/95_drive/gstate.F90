!{\src2tex{textfont=tt}}
!!****f* ABINIT/gstate
!! NAME
!! gstate
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations by CG minimization.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, JYR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial CPU time
!!
!! OUTPUT
!!  npwtot(nkpt) = total number of plane waves at each k point
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!
!! SIDE EFFECTS
!!  acell(3)=unit cell length scales (bohr)
!!  amu(ntypat)=mass of each atom type
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | mband =maximum number of bands (IN)
!!   | mgfft =maximum single fft dimension (IN)
!!   | mkmem =number of k points treated by this processor (IN)
!!   | mpw   =maximum number of planewaves in basis sphere (large number) (IN)
!!   | natom =number of atoms in unit cell (IN)
!!   | nfft  =(effective) number of FFT grid points (for this processor) (IN)
!!   | nkpt  =number of k points (IN)
!!   | nspden=number of spin-density components (IN)
!!   | nsppol=number of channels for spin-polarization (1 or 2) (IN)
!!   | nsym  =number of symmetry elements in space group
!!  iexit= exit flag
!!  initialized= 0 for the first GS calculation (not initialized), else 1
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!   some others to be initialized here)
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in gstate, a significant part of
!!   psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,
!!     ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!   All the remaining components of psps are to be initialized in the call
!!   to pspini .
!!   The next time the code enters gstate, psps might be identical to the
!!   one of the previous dtset, in which case, no reinitialisation is scheduled
!!   in pspini.f .
!!  rprim(3,3)=dimensionless real space primitive translations
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  vel(3,natom)=value of velocity
!!  vel_cell(3,3)=value of cell parameters velocity
!!  xred(3,natom) = reduced atomic coordinates
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!! In case of wavelets:
!! --------------------
!!    - Only the usual FFT grid (defined by wvl_crmult) is used.
!!      It is defined by nfft, ngfft, mgfft, ... This is strictly not
!!      an FFT grid since its dimensions are not suited for FFTs. They are
!!      defined by wvl_setngfft().
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! TODO
!! Not yet possible to use restartxf in parallel when localrdwf==0
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!      alloc_hamilt_gpu,cleanrec,clnup1,clnup2,ddb_free,ddb_malloc
!!      dealloc_hamilt_gpu,bandfft_kpt_destroy_array,destroy_bfield,destroy_efield
!!      destroy_electronpositron,destroy_gemm_nonlop,destroy_invovl
!!      destroy_sc_dmft,ebands_free,ebands_init,energies_init,exit_check,fixsym
!!      fourdp,getph,hdr_free,hdr_init,hdr_update,bandfft_kpt_init1
!!      init_eb_field_vars,init_electronpositron,init_gemm_nonlop,init_invovl
!!      init_sc_dmft,init_scf_history,initrec,initrhoij,initro,initylmg,inwffil
!!      ioarr,ioddb8_out,jellium,kpgio,local_potential_dimensions
!!      lotfparam_init,matr3inv,mkradim,mkrho,mover,newocc
!!      nullify_gaussian_basis,nullify_wvl_data,outqmc,outwf,outxfhist,paw2wvl
!!      paw_gencond,pawfgr_destroy,pawfgr_init,pawinit,pawpuxinit,pawrhoij_copy
!!      pawrhoij_destroy,pawuj_drive,print_sc_dmft,prtefield,prtene,psddb8
!!      psolver_kernel,pspini,readocc_dmft,reportgap,scfcv_init,scfcv_init2
!!      scfcv_new2,setsym,setsymrhoij,setup1,setup2,status,timab,transgrid
!!      wffclose,wffreadskiprec,write_blok8,wrtout,wvl_denspot_free
!!      wvl_denspot_set,wvl_descr_atoms_set,wvl_descr_atoms_set_sym
!!      wvl_descr_free,wvl_descr_psp_set,wvl_initro,wvl_mkrho,wvl_occ_abi2big
!!      wvl_paw_free,wvl_projectors_free,wvl_projectors_set,wvl_setboxgeometry
!!      wvl_setngfft,wvl_timing,wvl_wfs_free,wvl_wfs_lr_copy,wvl_wfs_set
!!      xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gstate(acell,amu,codvsn,cpui,dtfil,dtset,iexit,initialized,&
&                 mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,&
&                 psps,results_gs,rprim,scf_history,vel,vel_cell,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 !use defs_scftypes
 use defs_wvltypes
 !use defs_scfcvargs
 use defs_parameters
 use defs_rectypes
 use m_scf_history
 use m_abimover
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_exit
 use m_wffile
 use m_rec
 use m_bfield
 use m_efield
 use m_ddb
 use m_bandfft_kpt
 use m_invovl
 use m_gemm_nonlop

 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_pawfgr,           only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_abi2big,          only : wvl_occ_abi2big
 use m_energies,         only : energies_type, energies_init
 use m_results_gs ,      only : results_gs_type
 use m_pawrhoij,         only : pawrhoij_type, pawrhoij_copy, pawrhoij_destroy
 use m_paw_dmft,         only : init_sc_dmft,destroy_sc_dmft,print_sc_dmft,paw_dmft_type,readocc_dmft
 use m_data4entropyDMFT, only : data4entropyDMFT_t, data4entropyDMFT_init, data4entropyDMFT_destroy
 use m_electronpositron, only : electronpositron_type,init_electronpositron,destroy_electronpositron, &
&                               electronpositron_calctype
 use m_header,           only : hdr_init, hdr_free, hdr_update
 use m_ebands,           only : ebands_init, ebands_free, ReportGap
 use m_scfcv,            only : scfcv_t,scfcv_init, scfcv_destroy, scfcv_run

#if defined HAVE_DFT_BIGDFT
 use BigDFT_API,         only : wvl_timing => timing, local_potential_dimensions,nullify_gaussian_basis, copy_coulomb_operator,&
& deallocate_coulomb_operator
#endif
#if defined HAVE_LOTF
 use defs_param_lotf,    only : lotfparam_init
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gstate'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_43_wvl_wrappers
#if defined HAVE_GPU_CUDA
 use interfaces_52_manage_cuda
#endif
 use interfaces_53_ffts
 use interfaces_56_io_mpi
 use interfaces_56_recipspace
 use interfaces_57_iovars
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_62_poisson
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_72_response
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => gstate
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: iexit,initialized
 real(dp),intent(in) :: cpui
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 type(results_gs_type),intent(inout) :: results_gs
 type(scf_history_type),target,intent(inout) :: scf_history
!arrays
 integer,intent(out) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: acell(3),amu(psps%ntypat),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: rprim(3,3),vel(3,dtset%natom),vel_cell(3,3),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might change soon ...
!2   for wavefunction file, new format (version 2.0 and after)    (fform)   NOT USED
!52  for density rho(r)       (fformr)
!102 for potential V(r) file. (fformv)  NOT USED
!scalars
 integer,parameter :: formeig=0,level=101,response=0 
 integer :: ndtpawuj=0  ! Cannot use parameter because scfargs points to this! Have to get rid of pointers to scalars!
#if defined HAVE_DFT_BIGDFT
 integer :: icoulomb
#endif
 integer :: accessfil,ask_accurate,bantot,choice,comm_psp,fformr=52,fullinit
 integer :: gnt_option,gscase,iatom,idir,ierr,ii,indx,jj,kk,ios,itypat
 integer :: ixfh,izero,master,mcg,me,mgfftf,mpert,msize,mu,my_natom,my_nspinor
 integer :: nblok,nfftf,nfftot
 integer :: openexit,option,optorth,psp_gencond,conv_retcode
 integer :: pwind_alloc,rdwr,rdwrpaw,spaceworld,tim_mkrho,use_sc_dmft
 real(dp) :: cpus,ecore,ecut_eff,ecutdg_eff,etot,fermie
 real(dp) :: gsqcut_eff,gsqcut_shp,gsqcutc_eff,residm,tolwfr,ucvol
 logical :: read_wf_or_den,has_to_init,call_pawinit,write_wfk
 logical :: wvlbigdft=.false.
 character(len=500) :: message
 character(len=fnlen) :: ddbnm,dscrpt,filnam
 real(dp) :: fatvshift
 type(ebands_t) :: bstruct,tmp_Bst
 type(bfield_type) :: dtbfield
 type(efield_type) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(hdr_type) :: hdr
 !type(macro_uj_type) :: dtpawuj(0:ndtpawuj)
 type(macro_uj_type) :: dtpawuj(0)
 type(paw_dmft_type) :: paw_dmft
 type(pawfgr_type) :: pawfgr
 type(recursion_type) ::rec_set
 type(wffile_type) :: wff1,wffnew,wffnow
 type(wvl_data) :: wvl
 !type(ab_scfcv_args_in) :: ab_scfcv_in
 !type(ab_scfcv_args_inout) :: ab_scfcv_inout
 type(ab_xfh_type) :: ab_xfh
 !type(ab_scfcvargs) :: scfcv_args
 type(ddb_type) :: ddb
 type(scfcv_t) :: scfcv_args
!arrays
 integer :: ngfft(18),ngfftf(18)
 integer,allocatable :: atindx(:),atindx1(:),indsym(:,:,:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),nattyp(:),symrec(:,:,:)
 integer,allocatable,target :: npwarr(:)
 integer,pointer :: npwarr_(:),pwind(:,:,:)
!real(dp) :: xcart(3,dtset%natom)
 real(dp) :: efield_band(3),gmet(3,3),gprimd(3,3)
 real(dp) :: rmet(3,3),rprimd(3,3),rprimd_orig(3,3),tsec(2)
 real(dp),allocatable :: amass(:),cg(:,:),doccde(:)
 real(dp),allocatable :: eigen(:),ph1df(:,:),phnons(:,:,:),resid(:),rhowfg(:,:)
 real(dp),allocatable :: rhowfr(:,:),spinat_dum(:,:),start(:,:),work(:)
 real(dp),allocatable :: ylm(:,:),ylmgr(:,:,:)
 real(dp),pointer :: pwnsfac(:,:),rhog(:,:),rhor(:,:)
 real(dp),pointer :: taug(:,:),taur(:,:),xred_old(:,:)
 type(pawrhoij_type),pointer :: pawrhoij(:)
 type(coulomb_operator) :: kernel_dummy
! ***********************************************************************


! chen 
 real(dp),allocatable :: chen_xcart(:,:)




 DBG_ENTER("COLL")

 call timab(32,1,tsec)
 call timab(33,3,tsec)

!###########################################################
!### 01. Initializations XML, MPI, WVL, etc

!Init MPI data
 master =0
 spaceworld=mpi_enreg%comm_cell
 me=xcomm_rank(spaceworld)
!Set up MPI informations from the dataset
 my_natom=mpi_enreg%my_natom

!Nullify wvl_data. It is important to do so irregardless of the value of usewvl:
 call nullify_wvl_data(wvl)

!Set up informations when wavelets are in use
 if (dtset%usewvl == 1) then

!  If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
   if(dtset%wvl_bigdft_comp==1) wvlbigdft=.true.

!  Default value, to be set-up elsewhere.
   wvl%descr%h(:)                 = dtset%wvl_hgrid

#if defined HAVE_DFT_BIGDFT
   wvl%descr%paw%usepaw=psps%usepaw
   wvl%descr%paw%natom=dtset%natom
#endif

!  We set the atom-related internal wvl variables.
   call wvl_descr_atoms_set(acell, dtset%icoulomb, dtset%natom, &
&   dtset%ntypat, dtset%typat, wvl%descr)
   if(dtset%usepaw==0) then
!    nullify PAW proj_G in NC case:
#if defined HAVE_DFT_BIGDFT
     ABI_DATATYPE_ALLOCATE(wvl%projectors%G,(dtset%ntypat))
     do itypat=1,dtset%ntypat
       call nullify_gaussian_basis(wvl%projectors%G(itypat))
     end do
#endif
   end if

   wvl%descr%exctxpar = "OP2P"
 end if

 if (me == 0 .and. dtset%prtxml == 1) then
!  gstate() will handle a dataset, so we output the dataSet markup.
   write(ab_xml_out, "(A)") '  <dataSet>'
!  We output the variables of the dataset given in argument.
!  call outvarsXML()
 end if

!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)

!Set flag for IO fileformat.
 accessfil = 0
 if (dtset%accesswff == IO_MODE_MPI ) accessfil = 4
 if (dtset%accesswff == IO_MODE_ETSF) accessfil = 3

!Structured debugging if prtvol==-level
 if(dtset%prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,' gstate : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

!###########################################################
!### 02. Calls setup1, kpgio, initylmg

 ecore=zero
 results_gs%pel(1:3)   =zero
 results_gs%grewtn(:,:)=zero
!MT Feb 2012: I dont know why but grvdw has to be allocated
!when using BigDFT to ensure success on inca_gcc44_sdebug
 if (dtset%vdw_xc==5.or.dtset%usewvl==1) then
   results_gs%ngrvdw=dtset%natom
   if (allocated(results_gs%grvdw)) then
     ABI_DEALLOCATE(results_gs%grvdw)
   end if
   ABI_ALLOCATE(results_gs%grvdw,(3,dtset%natom))
   results_gs%grvdw(:,:)=zero
 end if
 call energies_init(results_gs%energies)

!Set up for iterations
 ABI_ALLOCATE(amass,(dtset%natom))
 call setup1(acell,amass,amu,bantot,dtset,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
& dtset%natom,ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
& response,rmet,rprim,rprimd,ucvol,psps%usepaw)

 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 if (dtset%usewvl == 0 .and. dtset%tfkinfunc /= 2) then
!  Set up the basis sphere of planewaves
   ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
   call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,dtfil%fnametmp_kg, &
&   dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,&
&   dtset%mpw,npwarr,npwtot,dtset%nsppol,dtfil%unkg)
   call bandfft_kpt_init1(bandfft_kpt,dtset%istwfk,kg,dtset%mgfft,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,dtset%nkpt,npwarr,dtset%nsppol)
 else
   ABI_ALLOCATE(kg,(0,0))
   npwarr(:) = 0
   npwtot(:) = 0
 end if

 if(dtset%wfoptalg == 1 .and. psps%usepaw == 1) then
   call init_invovl(dtset%nkpt)
 end if

 if(dtset%use_gemm_nonlop == 1 .and. dtset%use_gpu_cuda/=1) then
   ! set global variable
   gemm_nonlop_use_gemm = .true.
   call init_gemm_nonlop(dtset%nkpt)
 else
   gemm_nonlop_use_gemm = .false.
 end if
 
 ! PW basis set: test if the problem is ill-defined.
 if (dtset%usewvl == 0) then
   if (dtset%mband > dtset%mpw) then 
     ! No way we can solve the problem. Abort now!
     write(message,"(2(a,i0),4a)")&
&     "Number of bands mband ",dtset%mband," > number of planewaves mpw = ",dtset%mpw,ch10,&
&     "The number of eigenvectors cannot be greater that the size of the Hamiltonian!",ch10,&
&     "Action: decrease nband or, alternatively, increase ecut"
     ! FIXME: This should be an Error but LOFT in paral stops. Xavier will fix it.
     MSG_WARNING(message)
     !MSG_ERROR(message)

   else if (dtset%mband >= 0.9 * dtset%mpw) then 
     ! Warn the user
     write(message,"(a,i0,a,f6.1,4a)")&
&     "Number of bands mband ",dtset%mband," >= 0.9 * maximum number of planewaves = ",0.9 * dtset%mpw,ch10,&
&     "The problem is ill-defined and the GS algorithm will show numerical instabilities!",ch10,&
&     "Assume experienced user. Execution will continue."
     MSG_WARNING(message)
   end if 
 end if
 
!Set up the Ylm for each k point
 if ( dtset%tfkinfunc /= 2) then
   ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_ALLOCATE(ylmgr,(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm))
   if (psps%useylm==1) then
     option=0
     if (dtset%prtstm==0.and.dtset%iscf>0.and.dtset%positron/=1) option=1 ! compute gradients of YLM
     if (dtset%orbmag > 0) option=1 ! compute gradients of YLM
     if (dtset%berryopt==4 .and. dtset%optstress /= 0 .and. psps%usepaw==1) option = 1 ! compute gradients of YLM
     call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,&
&     psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
   end if
 else
   ABI_ALLOCATE(ylm,(0,0))
   ABI_ALLOCATE(ylmgr,(0,0,0))
 end if

!SCF history management (allocate it at first call)
 has_to_init=(initialized==0.or.scf_history%history_size<0)
 if (initialized==0) then
!  This call has to be done before any use of SCF history
   call init_scf_history(dtset,mpi_enreg,scf_history)
 end if

 call timab(33,2,tsec)
 call timab(701,3,tsec)

!###########################################################
!### 03. Calls pspini

!Open and read pseudopotential files
 comm_psp=mpi_enreg%comm_cell;if (dtset%usewvl==1) comm_psp=mpi_enreg%comm_wvl
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,level,&
& pawrad,pawtab,psps,rprimd,comm_mpi=comm_psp)

 call timab(701,2,tsec)
 call timab(33,3,tsec)

!In case of isolated computations, ecore must set to zero
!because its contribution is counted in the ewald energy as the ion-ion interaction.
 if (dtset%icoulomb == 1) ecore = zero

!WVL - Now that psp data are available, we compute rprimd, acell... from the atomic positions.
 if (dtset%usewvl == 1) then
   call wvl_descr_psp_set(trim(dtfil%filnam_ds(3)) // "_OCCUP", dtset%nsppol, psps, wvl%descr)
   call wvl_setBoxGeometry(mpi_enreg%me_wvl, dtset%prtvol, psps%gth_params%radii_cf, rprimd, xred, &
&   wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult)
   call mkradim(acell,rprim,rprimd)
   call wvl_denspot_set(wvl%den, psps%gth_params, dtset%ixc, &
&   me, dtset%natom, mpi_enreg%nproc_wvl, dtset%nsppol, rprimd, wvl%descr, &
&   dtset%wvl_crmult, dtset%wvl_frmult, xred)
!  TODO: to be moved in a routine.
#if defined HAVE_DFT_BIGDFT
   if (wvl%descr%atoms%astruct%geocode == "F") then
     icoulomb = 1
   else if (wvl%descr%atoms%astruct%geocode == "S") then
     icoulomb = 2
   else
     icoulomb = 0
   end if
!  calculation of the Poisson kernel anticipated to reduce memory peak for small systems
   call psolver_kernel( wvl%den%denspot%dpbox%hgrids, 1, icoulomb, mpi_enreg%me_wvl, wvl%den%denspot%pkernel , &
&   mpi_enreg%comm_wvl, wvl%den%denspot%dpbox%ndims, mpi_enreg%nproc_wvl, dtset%nscforder)
   nullify(wvl%den%denspot%pkernelseq%kernel)
   !call copy_coulomb_operator(wvl%den%denspot%pkernel,wvl%den%denspot%pkernelseq,ABI_FUNC)
!  Associate the denspot distribution into mpi_enreg.
   mpi_enreg%nscatterarr  => wvl%den%denspot%dpbox%nscatterarr
   mpi_enreg%ngatherarr   => wvl%den%denspot%dpbox%ngatherarr
   mpi_enreg%ngfft3_ionic =  wvl%den%denspot%dpbox%n3pi
   call wvl_setngfft(mpi_enreg%me_wvl, dtset%mgfft, dtset%nfft, &
&   dtset%ngfft, mpi_enreg%nproc_wvl, wvl%den%denspot%dpbox%ndims(1), &
&   wvl%den%denspot%dpbox%ndims(2),wvl%den%denspot%dpbox%ndims(3),&
&   wvl%den%denspot%dpbox%nscatterarr(mpi_enreg%me_wvl, 2))
#endif
   nfftf     = dtset%nfft
   mgfftf    = dtset%mgfft
   ngfftf(:) = dtset%ngfft(:)
   ngfft     = dtset%nfft
!  Recalculate gprimd
   call matr3inv(rprimd,gprimd)
!  PAW section
   if(psps%usepaw==1) then
!    Reinitialize Pawfgr with new values of ngfft
     call pawfgr_destroy(pawfgr)
     call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)
!    fill wvl objects from paw objects
!    wvl%descr%paw%usepaw=dtset%usepaw
     call paw2wvl(pawtab,wvl%projectors,wvl%descr)
   end if
 end if

!Initialize band structure datatype
 ABI_ALLOCATE(doccde,(bantot))
 ABI_ALLOCATE(eigen,(bantot))
 doccde(:)=zero ; eigen(:)=zero
 if (dtset%paral_kgb/=0) then     !  We decide to store total npw in bstruct,
   ABI_ALLOCATE(npwarr_,(dtset%nkpt))
   npwarr_(:)=npwarr(:)
   call xmpi_sum(npwarr_,mpi_enreg%comm_bandfft,ierr)
 else
   npwarr_ => npwarr
 end if

 call ebands_init(bantot,bstruct,dtset%nelect,doccde,eigen,dtset%istwfk,dtset%kptns,&
& dtset%nband,dtset%nkpt,npwarr_,dtset%nsppol,dtset%nspinor,dtset%tphysel,&
& dtset%tsmear,dtset%occopt,occ,dtset%wtk)

 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen)
 if (dtset%paral_kgb/=0)  then
   ABI_DEALLOCATE(npwarr_)
 end if
 nullify(npwarr_)

!Initialize PAW atomic occupancies
 if (scf_history%history_size>=0) then
   pawrhoij => scf_history%pawrhoij_last
 else
   ABI_DATATYPE_ALLOCATE(pawrhoij,(my_natom*psps%usepaw))
 end if
 if (psps%usepaw==1.and.has_to_init) then
   call initrhoij(dtset%pawcpxocc,dtset%lexexch,&
&   dtset%lpawu,my_natom,dtset%natom,dtset%nspden,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,pawrhoij,dtset%pawspnorb,pawtab,dtset%spinat,dtset%typat,&
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps,wvl%descr,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
 call hdr_update(bantot,etot,fermie,hdr,dtset%natom,residm,rprimd,occ,pawrhoij,psps%usepaw,xred,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!Clean band structure datatype (should use it more in the future !)
 call ebands_free(bstruct)

!###########################################################
!### 04. Symmetry operations when nsym>1

!Do symmetry stuff only for nsym>1
 if (dtset%usewvl == 0) then
   nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 else
#if defined HAVE_DFT_BIGDFT
   nfftot=product(wvl%den%denspot%dpbox%ndims)
#else
   BIGDFT_NOTENABLED_ERROR()
#endif
 end if
 ABI_ALLOCATE(irrzon,(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(indsym,(4,dtset%nsym,dtset%natom))
 ABI_ALLOCATE(symrec,(3,3,dtset%nsym))
 irrzon(:,:,:)=0
 phnons(:,:,:)=zero
 indsym(:,:,:)=0
 symrec(:,:,:)=0

 if (dtset%nsym>1) then
   call setsym(indsym,irrzon,dtset%iscf,dtset%natom,&
&   nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
&   phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!  Make sure dtset%iatfix does not break symmetry
   call fixsym(dtset%iatfix,indsym,dtset%natom,dtset%nsym)
 else
!  The symrec array is used by initberry even in case nsym = 1
   symrec(:,:,1) = 0
   symrec(1,1,1) = 1 ; symrec(2,2,1) = 1 ; symrec(3,3,1) = 1
 end if
 if (dtset%usewvl == 1) then
   call wvl_descr_atoms_set_sym(wvl%descr, dtset%efield, irrzon, &
&   dtset%nsym, phnons, dtset%symafm, dtset%symrel, dtset%tnons, dtset%tolsym)
#if defined HAVE_DFT_BIGDFT
   wvl%den%symObj = wvl%descr%atoms%astruct%sym%symObj
#endif
 end if

!###########################################################
!### 05. Calls inwffil

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 mcg=dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol

 ABI_ALLOCATE(cg,(2,mcg))
 ABI_CHECK_ALLOC("out of memory in cg")

 ABI_ALLOCATE(eigen,(dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(resid,(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen(:)=zero ; resid(:)=zero
!mpi_enreg%paralbd=0 ; ask_accurate=0
 ask_accurate=0

!WVL - Branching, allocating wavefunctions as wavelets.
 if (dtset%usewvl == 1) then
   call wvl_wfs_lr_copy(wvl%wfs, wvl%descr)
!  Create access arrays for wavefunctions and allocate wvl%wfs%psi (other arrays are left unallocated).
   call wvl_wfs_set(dtset%strprecon,dtset%spinmagntarget, dtset%kpt, mpi_enreg%me_wvl,&
&   dtset%natom, sum(dtset%nband), &
&   dtset%nkpt, mpi_enreg%nproc_wvl, dtset%nspinor, dtset%nsppol, dtset%nwfshist, dtset%occ_orig, &
&   psps, rprimd, wvl%wfs, dtset%wtk, wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult, &
&   xred)
!  We transfer wavelets informations to the hdr structure.
#if defined HAVE_DFT_BIGDFT
   call local_potential_dimensions(wvl%wfs%ks%lzd,wvl%wfs%ks%orbs,wvl%den%denspot%dpbox%ngatherarr(0,1))
   hdr%nwvlarr(1) = wvl%wfs%ks%lzd%Glr%wfd%nvctr_c
   hdr%nwvlarr(2) = 7 * wvl%wfs%ks%lzd%Glr%wfd%nvctr_f
#endif
!  Create access arrays for projectors and allocate them.
!  Compute projectors from each atom.
   call wvl_projectors_set(mpi_enreg%me_wvl, dtset%natom, wvl%projectors, psps, rprimd, &
&   wvl%wfs, wvl%descr, dtset%wvl_frmult, xred)
 end if

 read_wf_or_den=(dtset%iscf<=0.or.dtfil%ireadwf/=0.or.(dtfil%ireadden/=0.and.dtset%positron<=0))
 read_wf_or_den=(read_wf_or_den.and.has_to_init)

!RECURSION -  initialization
 if(has_to_init .and. dtset%userec==1) then
   call InitRec(dtset,mpi_enreg,rec_set,rmet,maxval(psps%indlmn(3,:,:)))
 end if

!LOTF - initialization
#if defined HAVE_LOTF
 if(has_to_init .and. dtset%ionmov==23) then
   call lotfparam_init(dtset%natom,dtset%lotf_version,1,&
&   dtset%lotf_nitex,dtset%lotf_nneigx,&
&   dtset%lotf_classic,1,1)
 end if
#endif

!Initialize wavefunctions.
 if(dtset%tfkinfunc /=2) then
   wff1%unwff=dtfil%unwff1
   optorth=1   !if (psps%usepaw==1) optorth=0
   if(psps%usepaw==1 .and. dtfil%ireadwf==1)optorth=0
   call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen,&
&   dtset%exchn2n3d,formeig,gmet,hdr,dtfil%ireadwf,dtset%istwfk,kg,&
&   dtset%kptns,dtset%localrdwf,dtset%mband,mcg,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,ngfft,dtset%nkpt,npwarr,&
&   dtset%nsppol,dtset%nsym,occ,optorth,rprimd,dtset%symafm,&
&   dtset%symrel,dtset%tnons,dtfil%unkg,wff1,wffnow,dtfil%unwff1,&
&   dtfil%unwft1,dtfil%fnamewffk,dtfil%fnametmp_wf1,wvl)
 end if

 if (psps%usepaw==1.and.dtfil%ireadwf==1)then
   call pawrhoij_copy(hdr%pawrhoij,pawrhoij,mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
!  Has to update header again (because pawrhoij has changed)  -  MT 2007-10-22: Why ?
!  call hdr_update(bantot,etot,fermie,hdr,residm,rprimd,occ,pawrhoij,psps%usepaw,xred, &
!  &                  mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

!###########################################################
!### 06. Operations related to restartxf (Old version)

!Initialize xf history (should be put in inwffil)
 ab_xfh%nxfh=0
 if(dtset%restartxf>=1 .and. dtfil%ireadwf==1)then

!  Should exchange the data about history in parallel localrdwf==0
   if(xmpi_paral==1 .and. dtset%localrdwf==0)then
     write(message, '(a,a,a)' )&
&     'It is not yet possible to use non-zero restartxf,',ch10,&
&     'in parallel, when localrdwf=0. Sorry for this ...'
     MSG_BUG(message)
   end if

   ABI_ALLOCATE(ab_xfh%xfhist,(3,dtset%natom+4,2,0))
   call outxfhist(ab_xfh,dtset%natom,2,wff1,ios)
   ABI_DEALLOCATE(ab_xfh%xfhist)

   if(ios>0)then
     write(message,'(a,a,a)')&
&     'An error occurred reading the input wavefunction file,',ch10,&
&     'with restartxf=1.'
     MSG_ERROR(message)
   else if(ios==0)then
     write(message, '(a,a,i4,a)' )ch10,&
&     ' gstate : reading',ab_xfh%nxfh,' (x,f) history pairs from input wf file.'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
   end if
!  WARNING : should check that restartxf is not negative
!  WARNING : should check that restartxf /= only when dtfil%ireadwf is activated
 end if

!Allocate the xf history array : takes into account the existing
!pairs, minus those that will be discarded, then those that will
!be computed, governed by dtset%ntime, and some additional pairs
!(needed when it will be possible to use xfhist for move.f)
 ab_xfh%mxfh=(ab_xfh%nxfh-dtset%restartxf+1)+dtset%ntime+5
 ABI_ALLOCATE(ab_xfh%xfhist,(3,dtset%natom+4,2,ab_xfh%mxfh))
 ab_xfh%xfhist(:,:,:,:) = zero
!WARNING : should check that the number of atoms in the wf file and natom are the same

!Initialize the xf history array
 if(ab_xfh%nxfh>=dtset%restartxf .and. ab_xfh%nxfh>0)then
!  Eventually skip some of the previous history
   if(dtset%restartxf>=2)then
     do ixfh=1,dtset%restartxf-1
       call WffReadSkipRec(ios,1,wff1)
     end do
   end if

!  Read and store the relevant history
   ab_xfh%nxfhr=ab_xfh%nxfh-dtset%restartxf+1
   call outxfhist(ab_xfh,dtset%natom,3,wff1,ios)
 end if

!Close wff1, if it was ever opened (in inwffil)
 if (dtfil%ireadwf==1) then
   call WffClose(wff1,ierr)
 end if

!###########################################################
!### 07. Calls setup2

!Further setup
 ABI_ALLOCATE(start,(3,dtset%natom))
 call setup2(dtset,npwtot,start,wvl%wfs,xred)

!Allocation of previous atomic positions
 if (scf_history%history_size>=0) then
   xred_old => scf_history%xred_last
 else
   ABI_ALLOCATE(xred_old,(3,dtset%natom))
 end if
 if (has_to_init) xred_old=xred

!Timing for initialisation period
 call timab(33,2,tsec)
 call timab(34,3,tsec)

!###########################################################
!### 08. Compute new occupation numbers

!Compute new occupation numbers, in case wavefunctions and eigenenergies
!were read from disk, occupation scheme is metallic (this excludes iscf=-1),
!and occupation numbers are required by iscf
 if( dtfil%ireadwf==1 .and. &
& (dtset%occopt>=3.and.dtset%occopt<=8) .and. &
& (dtset%iscf>0 .or. dtset%iscf==-3) .and. dtset%positron/=1 ) then

   ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
!  Warning : ideally, results_gs%entropy should not be set up here XG 20011007
!  Do not take into account the possible STM bias
   call newocc(doccde,eigen,results_gs%energies%entropy,&
&   results_gs%energies%e_fermie,&
&   dtset%spinmagntarget,dtset%mband,dtset%nband,&
&   dtset%nelect,dtset%nkpt,dtset%nspinor,dtset%nsppol,occ,&
&   dtset%occopt,dtset%prtvol,zero,dtset%tphysel,dtset%tsmear,dtset%wtk)
   ABI_DEALLOCATE(doccde)

!  Transfer occupations to bigdft object:
   if(dtset%usewvl==1 .and. .not. wvlbigdft) then
     call wvl_occ_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,occ,1,wvl%wfs)
!    call wvl_energies_abi2big(results_gs%energies,wvl%wfs,2)
   end if

 else
!  Warning : ideally, results_gs%entropy should not be set up here XG 20011007
   results_gs%energies%entropy=zero
 end if

!###########################################################
!### 09. Generate an index table of atoms

!Definition of atindx array
!Generate an index table of atoms, in order for them to be used type after type.
 ABI_ALLOCATE(atindx,(dtset%natom))
 ABI_ALLOCATE(atindx1,(dtset%natom))
 ABI_ALLOCATE(nattyp,(psps%ntypat))
 indx=1
 do itypat=1,psps%ntypat
   nattyp(itypat)=0
   do iatom=1,dtset%natom
     if(dtset%typat(iatom)==itypat)then
       atindx(iatom)=indx
       atindx1(indx)=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

!Compute structure factor phases for current atomic pos:
 if ((.not.read_wf_or_den).or.(scf_history%history_size>0.and.has_to_init)) then
   ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*dtset%natom))
   call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 end if

!Here allocation of GPU for vtorho calculations
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call alloc_hamilt_gpu(atindx1,dtset,gprimd,mpi_enreg,nattyp,npwarr,2,psps,dtset%use_gpu_cuda)
 end if
#endif

!###########################################################
!### 10. PAW related operations

!Initialize paw_dmft, even if neither dmft not paw are used
!write(std_out,*) "dtset%usedmft",dtset%usedmft
 use_sc_dmft=dtset%usedmft
 if(dtset%paral_kgb>0) use_sc_dmft=0
 call init_sc_dmft(dtset%nbandkss,dtset%dmftbandi,dtset%dmftbandf,dtset%dmft_read_occnd,dtset%mband,&
& dtset%nband,dtset%nkpt,dtset%nspden, &
& dtset%nspinor,dtset%nsppol,occ,dtset%usedmft,paw_dmft,use_sc_dmft,dtset%dmft_solv,mpi_enreg)
 if (paw_dmft%use_dmft==1.and.me==0) then
   call readocc_dmft(paw_dmft,dtfil%filnam_ds(3),dtfil%filnam_ds(4))
 end if
 !Should be done inside init_sc_dmft
 if ( dtset%usedmft /= 0 ) then
   call data4entropyDMFT_init(paw_dmft%forentropyDMFT,&
                             dtset%natom,&
                             dtset%typat,&
                             dtset%lpawu,&
                             dtset%dmft_t2g==1, &
                             dtset%upawu,&
                             dtset%jpawu)
 end if
!write(std_out,*) "paw_dmft%use_dmft",paw_dmft%use_dmft

!PAW: 1- Initialize values for several arrays unchanged during iterations
!2- Initialize data for LDA+U
!3- Eventually open temporary storage file
 if(psps%usepaw==1) then
!  1-
   gnt_option=1;if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option=2

   ! Test if we have to call pawinit
   call paw_gencond(Dtset,gnt_option,"test",call_pawinit) 

   if (psp_gencond==1.or. call_pawinit) then
!    some gen-cond have to be added...
     call timab(553,1,tsec)
     gsqcut_shp=two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
     call pawinit(gnt_option,gsqcut_shp,dtset%pawlcutd,dtset%pawlmix,&
&     psps%mpsang,dtset%pawnphi,dtset%nsym,dtset%pawntheta,&
&     pawang,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%xclevel,dtset%usepotzero)

     ! Update internal values
     call paw_gencond(Dtset,gnt_option,"save",call_pawinit) 
     call timab(553,2,tsec)
#if defined HAVE_DFT_BIGDFT
!    In the PAW+WVL case, copy sij:
     if(dtset%usewvl==1) then
       do itypat=1,dtset%ntypat
         wvl%descr%paw%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
       end do
     end if
#endif
   end if
   psps%n1xccc=maxval(pawtab(1:psps%ntypat)%usetcore)
   call setsymrhoij(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,rprimd,symrec,pawang%zarot)

!  2-Initialize and compute data for LDA+U, EXX, or LDA+DMFT
   pawtab(:)%usepawu=0
   pawtab(:)%useexexch=0
   pawtab(:)%exchmix=zero
   if(paw_dmft%use_dmft==1) then
     call print_sc_dmft(paw_dmft,dtset%pawprtvol)
   end if
   if (dtset%usepawu>0.or.dtset%useexexch>0.or.paw_dmft%use_dmft>0) then
     call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%f4of2_sla,dtset%f6of2_sla,&
&     dtset%jpawu,dtset%lexexch,dtset%lpawu,dtset%ntypat,pawang,dtset%pawprtvol,&
&     pawrad,pawtab,dtset%upawu,dtset%usedmft,dtset%useexexch,dtset%usepawu)
   end if
 end if

!###########################################################
!### 11. Initialize (eventually) electron-positron data and 
!###     electric and magnetic field data

!Initialize (eventually) electron-positron data
 nullify (electronpositron)
 if (dtset%positron/=0) then
   call init_electronpositron(dtfil%ireadwf,dtset,electronpositron,mpi_enreg,nfftf,pawrhoij,pawtab)
 end if

!!Electric and magnetic field: initialization stage
!!further initialization and updates happen in scfcv.F90
!call init_eb_field_vars(dtbfield,dtefield,dtset,gmet,gprimd,kg,&
!& mpi_enreg,npwarr,occ,pawang,pawrad,pawtab,psps,&
!& pwind,pwind_alloc,pwnsfac,rprimd,symrec,xred)


!###########################################################
!### 12. Operations dependent of iscf value

!Get starting charge density : rhor as well as rhog
!Also initialize the kinetic energy density
 if (scf_history%history_size>=0) then
   rhor => scf_history%rhor_last
   taur => scf_history%taur_last
 else
   ABI_ALLOCATE(rhor,(nfftf,dtset%nspden))
   ABI_ALLOCATE(taur,(nfftf,dtset%nspden*dtset%usekden))
 end if
 ABI_ALLOCATE(rhog,(2,nfftf))
 ABI_ALLOCATE(taug,(2,nfftf*dtset%usekden))

 if (has_to_init) then
   if (dtset%iscf>0 .or. (dtset%iscf==0 .and. dtset%usewvl==1 )) then ! .and. dtset%usepaw==1)) then
     if(dtfil%ireadden/=0.and.dtset%positron<=0)then

       rdwr=1;rdwrpaw=psps%usepaw;if(dtfil%ireadwf/=0) rdwrpaw=0
       call ioarr(accessfil,rhor,dtset,results_gs%etotal,fformr,dtfil%fildensin,hdr,&
&       mpi_enreg,nfftf,pawrhoij,rdwr,rdwrpaw,wvl%den)
       if(dtfil%ireadkden/=0 .and. dtset%usekden==1 )then
         call ioarr(accessfil,taur,dtset,results_gs%etotal,fformr,dtfil%filkdensin,hdr,&
&         mpi_enreg,nfftf,pawrhoij,rdwr,rdwrpaw,wvl%den)
       end if
       if (rdwrpaw/=0) then
         call hdr_update(bantot,etot,fermie,hdr,dtset%natom,residm,rprimd,occ,pawrhoij,psps%usepaw,xred,&
&         mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
       end if
!      Compute up+down rho(G) by fft
       ABI_ALLOCATE(work,(nfftf))
       work(:)=rhor(:,1)
       call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
       if(dtset%usekden==1)then
         work(:)=taur(:,1)
         call fourdp(1,taug,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
       end if
       ABI_DEALLOCATE(work)

     else if(dtfil%ireadwf/=0)then
       izero=0
!      Obtain the charge density from wfs that were read previously
!      Be careful: in PAW, rho does not include the compensation
!      density (to be added in scfcv.F90) !
!      tim_mkrho=1 ; mpi_enreg%paralbd=0
       tim_mkrho=1
       if (psps%usepaw==1) then
         ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
         ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
!        write(std_out,*) "mkrhogstate"
         call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&         mpi_enreg,npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,wvl%wfs,wvl%den)
         call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
         ABI_DEALLOCATE(rhowfg)
         ABI_DEALLOCATE(rhowfr)
       else
         call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&         mpi_enreg,npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,wvl%wfs,wvl%den)
         if(dtset%usekden==1)then
           call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&           mpi_enreg,npwarr,occ,paw_dmft,phnons,taug,taur,rprimd,tim_mkrho,ucvol,wvl%wfs,wvl%den,option=1)
         end if

       end if

     else if(dtfil%ireadwf==0.and.dtset%positron/=1)then

!      Crude, but realistic initialisation of the density
!      There is not point to compute it from random wavefunctions
!      except with wavelets.
       if (dtset%usewvl == 0) then
         call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,&
&         mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,nattyp,nfftf,&
&         ngfftf,dtset%nspden,psps%ntypat,dtset%paral_kgb,pawtab,ph1df,&
&         psps%qgrid_vl,rhog,rhor,dtset%spinat,ucvol,psps%usepaw,&
&         dtset%ziontypat,dtset%znucl)
!        Update initialized density taking into account jellium slab
         if(dtset%jellslab/=0) then
           option=2
           ABI_ALLOCATE(work,(nfftf))
           call jellium(gmet,gsqcut_eff,mpi_enreg,nfftf,ngfftf,dtset%nspden,&
&           option,dtset%paral_kgb,dtset%slabwsrad,rhog,rhor,rprimd,work,dtset%slabzbeg,dtset%slabzend)
           ABI_DEALLOCATE(work)
         end if ! of usejell
!        Kinetic energy density initialized to zero (used only in metaGGAs ... )
         if(dtset%usekden==1)then
           taur=zero ; taug=zero
         end if
       else if (dtset%usewvl/=0) then
!        skip for the moment for wvl+paw, not yet coded
         if(dtset%usepaw==0) then
           call wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wvl%wfs, wvl%den)
         else !usepaw
#if defined HAVE_DFT_BIGDFT
           call wvl_initro(atindx1,wvl%descr%atoms%astruct%geocode,wvl%descr%h,mpi_enreg%me_wvl,&
&           dtset%natom,nattyp,nfftf,dtset%nspden,psps%ntypat,&
&           wvl%descr%Glr%d%n1,wvl%descr%Glr%d%n1i,&
&           wvl%descr%Glr%d%n2,wvl%descr%Glr%d%n2i,&
&           wvl%descr%Glr%d%n3,&
&           pawrad,pawtab,psps%gth_params%psppar,rhor,rprimd,&
&           dtset%spinat,wvl%den,dtset%xc_denpos,xred,dtset%ziontypat)

#endif
         end if
       end if

     end if

   else if ((dtset%iscf==-1.or.dtset%iscf==-2.or.dtset%iscf==-3).and.dtset%positron<=0) then

     call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
!    Read rho(r) from a disk file
     rdwr=1;rdwrpaw=psps%usepaw
!    Note : results_gs%etotal is read here,
!    and might serve in the tddft routine, but it is contrary to the intended use of results_gs ...
!    Warning : should check the use of results_gs%e_fermie
!    Warning : should check the use of results_gs%residm
!    One might make them separate variables.

     call ioarr(accessfil,rhor,dtset,results_gs%etotal,fformr,dtfil%fildensin,hdr,&
&     mpi_enreg,nfftf,pawrhoij,rdwr,rdwrpaw,wvl%den)
     if(dtfil%ireadkden/=0 .and. dtset%usekden==1)then
       call ioarr(accessfil,taur,dtset,results_gs%etotal,fformr,dtfil%filkdensin,hdr,&
&       mpi_enreg,nfftf,pawrhoij,rdwr,rdwrpaw,wvl%den)
     end if

!    Compute up+down rho(G) by fft
     ABI_ALLOCATE(work,(nfftf))
     work(:)=rhor(:,1)
     call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
     if(dtset%usekden==1)then
       work(:)=taur(:,1)
       call fourdp(1,taug,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
     end if
     ABI_DEALLOCATE(work)

   end if
 end if ! has_to_init

!###########################################################
!### 13. If needed, initialize SCF history variables

!If needed, initialize atomic density in SCF history
 if (scf_history%history_size>0.and.has_to_init) then
!  If rhor is an atomic density, just store it in history
   if (.not.read_wf_or_den) then
     scf_history%atmrho_last(:)=rhor(:,1)
   else
!    If rhor is not an atomic density, has to compute rho_at(r)
     ABI_ALLOCATE(rhowfg,(2,nfftf))
     ABI_ALLOCATE(rhowfr,(nfftf,1))
     ABI_ALLOCATE(spinat_dum,(3,dtset%natom))
     spinat_dum=zero
     call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,mgfftf,mpi_enreg,&
&     psps%mqgrid_vl,dtset%natom,nattyp,nfftf,ngfftf,1,psps%ntypat,dtset%paral_kgb,pawtab,&
&     ph1df,psps%qgrid_vl,rhowfg,rhowfr,spinat_dum,ucvol,&
&     psps%usepaw,dtset%ziontypat,dtset%znucl)
     scf_history%atmrho_last(:)=rhowfr(:,1)
     ABI_DEALLOCATE(rhowfg)
     ABI_DEALLOCATE(rhowfr)
     ABI_DEALLOCATE(spinat_dum)
   end if
 end if

 if ((.not.read_wf_or_den).or.(scf_history%history_size>0.and.has_to_init))  then
   ABI_DEALLOCATE(ph1df)
 end if

!!Electric and magnetic field: initialization stage
!!further initialization and updates happen in scfcv.F90
 call init_eb_field_vars(dtbfield,dtefield,dtset,gmet,gprimd,kg,&
& mpi_enreg,npwarr,occ,pawang,pawrad,pawtab,psps,&
& pwind,pwind_alloc,pwnsfac,rprimd,symrec,xred)

 fatvshift=one
 rprimd_orig(:,:)=rprimd

#if defined HAVE_DFT_BIGDFT
 if (dtset%usewvl == 1) then
   call wvl_timing(me,'INIT','PR')
 end if
#endif

!Check whether exiting was required by the user. If found then do not start minimization steps
!At this first call to chkexi, initialize cpus, if it
!is non-zero (which would mean that no action has to be taken)
!Should do this in driver ...
 cpus=dtset%cpus
 if(abs(cpus)>1.0d-5)cpus=cpus+cpui
 openexit=1 ; if(dtset%chkexit==0) openexit=0
 call exit_check(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg%comm_cell,openexit)

!If immediate exit, and wavefunctions were not read, must zero eigenvalues
 if (iexit/=0) eigen(:)=zero

 call timab(34,2,tsec)

 conv_retcode = -1

 if (iexit==0) then

!  ###########################################################
!  ### 14. Move atoms and acell acording to ionmov value


!  Eventually symmetrize atomic coordinates over space group elements:
!  call symzat(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

!  If atoms are not being moved and U should not be determined,
!  use scfcv directly; else
!  call move, pawuj_drive or brdmin which in turn calls scfcv.

   call timab(35,3,tsec)

!   call scfcv_init(ab_scfcv_in,ab_scfcv_inout,atindx,atindx1,cg,cpus,&
!&   dtbfield,dtefield,dtfil,dtpawuj,dtset,ecore,eigen,hdr,iapp,&
!&   indsym,initialized,irrzon,kg,mcg,mpi_enreg,my_natom,nattyp,ndtpawuj,&
!&   nfftf,npwarr,occ,pawang,pawfgr,pawrad,&
!&   pawrhoij,pawtab,phnons,psps,pwind,pwind_alloc,pwnsfac,&
!&   rec_set,resid,results_gs,scf_history,&
!&   fatvshift,symrec,taug,taur,wvl,ylm,ylmgr)
!
!   call scfcv_init2(scfcv_args,ab_scfcv_in,ab_scfcv_inout,dtset,paw_dmft,wffnew,wffnow)
   call scfcv_init(scfcv_args,atindx,atindx1,cg,cpus,&
&  dtbfield,dtefield,dtfil,dtpawuj,dtset,ecore,eigen,hdr,&
&  indsym,initialized,irrzon,kg,mcg,mpi_enreg,my_natom,nattyp,ndtpawuj,&
&  nfftf,npwarr,occ,pawang,pawfgr,pawrad,pawrhoij,&
&  pawtab,phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,&
&  resid,results_gs,scf_history,fatvshift,&
&  symrec,taug,taur,wvl,ylm,ylmgr,paw_dmft,wffnew,wffnow)
   call dtfil_init_time(dtfil,0,mpi_enreg)

   write(message,'(a,80a)')ch10,('=',mu=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   if (dtset%ionmov==0) then

!    Should merge this call with the call for dtset%ionmov==4 and 5

     if (dtset%macro_uj==0) then


! chen 
       if ( .NOT. allocated(chen_xcart)) then 
         allocate(chen_xcart(3,scfcv_args%dtset%natom))
       endif 
800 continue  
! end of hack 


       !call scfcv_new2(scfcv_args,electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)
       call scfcv_run(scfcv_args,electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)

!=================================================================================
! For cluster and enviroment, we read xred from ECDA and re-run the scfcv_run 
!=================================================================================
       if (scfcv_args%dtset%calc_mode==11) then 
         ! ECDA asks env and cluster to update its atom coordiates 
         call c_wrtlog('[gstate] loading in new atom cooridates')
         open(file='new_coords.dat',unit=111,action='read',form='unformatted')
         if (scfcv_args%dtset%type_of_system == 1) then 
           do kk=1,scfcv_args%dtset%natom 
             read(111) chen_xcart(:,kk)
           enddo 
         endif 
         if (scfcv_args%dtset%type_of_system == 2) then 
           read(111) chen_xcart
         endif 
         close(111)
         call xcart2xred(scfcv_args%dtset%natom,rprimd,chen_xcart,xred)
         ! the first lock 
         open(file='done_end.dat',unit=111,action='write')
         write(111,*)11
         close(111)
         ! the 2nd lock
         open(file='done2_end.dat',unit=111,action='write')
         write(111,*)11
         close(111)
         call c_wrtlog('[gstate] loaded new_coords.dat and updated xred.')
         goto 800
       endif 
       deallocate(chen_xcart)
! end of hack 
!==============================       


     else
!      Conduct determination of U

       call pawuj_drive(scfcv_args,dtset,electronpositron,rhog,rhor,rprimd,xred,xred_old)
     end if

!    ========================================
!    New structure for geometry optimization
!    ========================================
   else if (dtset%ionmov>50.or.dtset%ionmov<=23) then

     ! TODO: return conv_retcode
     call mover(scfcv_args,ab_xfh,acell,amass,dtfil,&
&     electronpositron,rhog,rhor,rprimd,vel,vel_cell,xred,xred_old)

!    Compute rprim from rprimd and acell
     do kk=1,3
       do jj=1,3
         rprim(jj,kk)=rprimd(jj,kk)/acell(kk)
       end do
     end do

!    =========================================
!    New structure for geometry optimization
!    =========================================

   else ! Not an allowed option
     write(message, '(a,i12,a,a)' )&
&     'Disallowed value for ionmov=',dtset%ionmov,ch10,&
&     'Allowed values are:1,2,3,4,5,6,7,8,9,10,11,12,13,14,20,21 and 30'
     MSG_BUG(message)
   end if

   call scfcv_destroy(scfcv_args)

   call timab(35,2,tsec)

!  ###########################################################
!  ### 15. Final operations and output for gstate

 end if !  End of the check of hasty exit

 call timab(36,3,tsec)

 write(message, '(80a,a,a,a,a)' ) ('=',mu=1,80),ch10,ch10,&
& ' ----iterations are completed or convergence reached----',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Mark this GS computation as done
 initialized=1
#if defined HAVE_DFT_BIGDFT
 if (dtset%usewvl == 1) then
   call wvl_timing(me,'WFN_OPT','PR')
 end if
#endif

!Will be put here later.
!! ! WVL - maybe compute the tail corrections to energy
!! if (dtset%tl_radius > zero) then
!!    ! Store xcart for each atom
!!    allocate(xcart(3, dtset%natom))
!!    call xred2xcart(dtset%natom, rprimd, xcart, xred)
!!    ! Use the tails to improve energy precision.
!!    call wvl_tail_corrections(dtset, results_gs%energies, results_gs%etotal, &
!!         & mpi_enreg, occ, psps, vtrial, wvl, xcart)
!!    deallocate(xcart)
!! end if

!Update the header, before using it
!WARNING : There is a problem now (time of writing, 6.0.4, but was in ABINITv5 and had ever been there) to update
!the header with change of rprim. Might be due to the planewave basis set definition.
!Put the original rprimd .
 call hdr_update(bantot,results_gs%etotal,results_gs%energies%e_fermie,hdr,dtset%natom,&
& results_gs%residm,rprimd_orig,occ,pawrhoij,psps%usepaw,xred,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

 if(dtset%nqpt==0)filnam=dtfil%fnameabo_wfk
 if(dtset%nqpt==1)filnam=dtfil%fnameabo_wfq

 write_wfk = .True.
 ! Write wavefunctions file only if convergence was not achieved.
 if (dtset%prtwf==-1 .and. conv_retcode == 0) then
   write_wfk = .False.
   call wrtout(ab_out,"- gstate: GS cycle converged and prtwf=-1  --> Skip output of WFK file.","COLL")
 end if

 if (write_wfk) then
   call outwf(cg,dtset,eigen,filnam,hdr,kg,dtset%kptns,&
&   dtset%mband,mcg,dtset%mkmem,mpi_enreg,dtset%mpw,ab_xfh%mxfh,dtset%natom,&
&   dtset%nband,dtset%nkpt,npwarr,dtset%nsppol,dtset%nstep,&
&   ab_xfh%nxfh,occ,resid,response,dtfil%unwff2,wvl%wfs,wvl%descr,ab_xfh%xfhist)
 end if

 if (dtset%prtwf==2)then
   call outqmc(cg,dtset,eigen,gprimd,hdr,kg,mcg,mpi_enreg,npwarr,occ,psps,results_gs)
 end if

 call clnup1(acell,dtset,eigen,results_gs%energies%e_fermie,&
& dtfil%fnameabo_dos,dtfil%fnameabo_eig,results_gs%fred,&
& mpi_enreg,nfftf,ngfftf,occ,dtset%optforces,&
& resid,rhor,rprimd,results_gs%vxcavg,xred)

 if ( (dtset%iscf>=0 .or. dtset%iscf==-3) .and. dtset%prtstm==0) then
   call prtene(dtset,results_gs%energies,ab_out,psps%usepaw)
 end if

!write final electric field components HONG

 if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt ==7 .or.  &
& dtset%berryopt == 14 .or. dtset%berryopt == 16 .or. dtset%berryopt ==17 ) then ! output final elctric field data    !!HONG
   if (dtset%berryopt == 4) then
     write(message,'(a,a)')   ch10, 'Constant unreduced E calculation  - final values:'
   end if
   if (dtset%berryopt == 6 ) then
     write(message,'(a,a)')   ch10, 'Constant unreduced D calculation  - final values:'
   end if
   if (dtset%berryopt == 14) then
     write(message,'(a,a)')   ch10, 'Constant reduced ebar calculation  - final values:'
   end if
   if (dtset%berryopt == 16 ) then
     write(message,'(a,a)')   ch10, 'Constant reduced d calculation  - final values:'
   end if

   if (dtset%berryopt == 17) then
     write(message,'(a,a)')   ch10, 'Constant reduced ebar and d calculation  - final values:'
   end if

   call wrtout(ab_out,message,'COLL')
   call prtefield(dtset,dtefield,ab_out,rprimd)

   call wrtout(std_out,message,'COLL')
   call prtefield(dtset,dtefield,std_out,rprimd)

!  To check if the final electric field is below the critical field
   do kk = 1, 3
     efield_band(kk) = abs(dtset%red_efieldbar(kk))*dtefield%nkstr(kk)
   end do
!  eg = maxval(eg_dir)
!  eg_ev = eg*Ha_eV
   write(message,'(a,a,a,a,a,a,a,a,f7.2,a,a)')ch10,&
&   ' Please check: COMMENT - ',ch10,&
&   '  As a rough estimate,',ch10,&
&   '  to be below the critical field, the bandgap of your system',ch10,&
&   '  should be larger than ',maxval(efield_band)*Ha_eV,' eV.',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   write(message,'(a)')  '--------------------------------------------------------------------------------'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

 end if

 ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
 doccde=zero

 call ebands_init(bantot,tmp_Bst,dtset%nelect,doccde,eigen,hdr%istwfk,hdr%kptns,hdr%nband,&
& hdr%nkpt,hdr%npwarr,hdr%nsppol,hdr%nspinor,hdr%tphysel,hdr%tsmear,hdr%occopt,hdr%occ,hdr%wtk)

 tmp_Bst%fermie = results_gs%energies%e_fermie

 ABI_DEALLOCATE(doccde)

!Compute and print the gaps. Store them in results_gs
 call ReportGap(tmp_BSt,header="Gap info",unit=std_out,mode_paral="COLL",gaps=results_gs%gaps)

 call ebands_free(tmp_BSt)

!Open the formatted derivative database file, and write the preliminary information
!In the // case, only one processor writes the energy and the gradients to the DDB

 if ((me==0).and.(dtset%nimage==1).and.((dtset%iscf > 0).or.&
& (dtset%berryopt == -1).or.(dtset%berryopt) == -3)) then

   dscrpt=' Note : temporary (transfer) database '
   ddbnm=trim(dtfil%filnam_ds(4))//'_DDB'
!  tolwfr must be initialized here, but it is a dummy value
   tolwfr=1.0_dp
   call ioddb8_out (dscrpt,ddbnm,dtset%natom,dtset%mband,&
&   dtset%nkpt,dtset%nsym,psps%ntypat,dtfil%unddb,DDB_VERSION,&
&   acell,amu,dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&   dtset%intxc,dtset%iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&   dtset%natom,dtset%nband,ngfft,dtset%nkpt,dtset%nspden,dtset%nspinor,&
&   dtset%nsppol,dtset%nsym,psps%ntypat,occ,dtset%occopt,dtset%pawecutdg,&
&   rprim,dtset%sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&   dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&   dtset%typat,dtset%usepaw,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

   if (dtset%iscf > 0) then
     nblok = 2          ! 1st blok = energy, 2nd blok = gradients
   else
     nblok = 1
   end if
   fullinit = 0 ; choice=2
   call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&   psps%lmnmax,nblok,psps%ntypat,dtfil%unddb,pawtab,&
&   psps%pspso,psps%usepaw,psps%useylm,DDB_VERSION)

   mpert = dtset%natom + 6   ; msize = 3*mpert

!  create a ddb structure with just one blok
   call ddb_malloc(ddb,msize,1,dtset%natom,dtset%ntypat)

   ddb%flg = 0
   ddb%qpt = zero
   ddb%nrm = one
   ddb%val = zero

!  Write total energy to the DDB
   if (dtset%iscf > 0) then
     ddb%typ(1) = 0
     ddb%val(1,1,1) = results_gs%etotal
     ddb%flg(1,1) = 1
     call write_blok8(ddb,1,choice,dtset%mband,&
&     mpert,msize,dtset%nkpt,dtfil%unddb)
   end if

!  Write gradients to the DDB
   ddb%typ = 4
   ddb%flg = 0
   ddb%val = zero
   indx = 0
   if (dtset%iscf > 0) then
     do iatom = 1, dtset%natom
       do idir = 1, 3
         indx = indx + 1
         ddb%flg(indx,1) = 1
         ddb%val(1,indx,1) = results_gs%fred(idir,iatom)
       end do
     end do
   end if

   indx = 3*dtset%natom + 3
   if ((abs(dtset%berryopt) == 1).or.(abs(dtset%berryopt) == 3)) then
     do idir = 1, 3
       indx = indx + 1
       if (dtset%rfdir(idir) == 1) then
         ddb%flg(indx,1) = 1
         ddb%val(1,indx,1) = results_gs%pel(idir)
       end if
     end do
   end if

   indx = 3*dtset%natom + 6
   if (dtset%iscf > 0) then
     ddb%flg(indx+1:indx+6,1) = 1
     ddb%val(1,indx+1:indx+6,1) = results_gs%strten(1:6)
   end if

   call write_blok8(ddb,1,choice,dtset%mband,mpert,msize,dtset%nkpt,dtfil%unddb)

   call ddb_free(ddb)

!  Close DDB
   close(dtfil%unddb)
 end if

 if (dtset%nstep>0 .and. dtset%prtstm==0 .and. dtset%positron/=1) then
   call clnup2(psps%n1xccc,results_gs%fred,results_gs%gresid,&
&   results_gs%grewtn,results_gs%grvdw,results_gs%grxc,dtset%iscf,dtset%natom,&
&   results_gs%ngrvdw,dtset%optforces,dtset%optstress,dtset%prtvol,start,&
&   results_gs%strten,results_gs%synlgr,xred)
 end if

!Deallocate arrays
 ABI_DEALLOCATE(amass)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(cg)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(resid)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(start)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(taug)
 ABI_DEALLOCATE(ab_xfh%xfhist)
 call pawfgr_destroy(pawfgr)

 if (dtset%usewvl == 0 .or. dtset%nsym <= 1) then
!  In wavelet case, irrzon and phnons are deallocated by wavelet object.
   ABI_DEALLOCATE(irrzon)
   ABI_DEALLOCATE(phnons)
 end if

 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylmgr)

 if (scf_history%history_size<0) then
   if (psps%usepaw==1) then
     call pawrhoij_destroy(pawrhoij)
   end if
   ABI_DEALLOCATE(rhor)
   ABI_DEALLOCATE(taur)
   ABI_DATATYPE_DEALLOCATE(pawrhoij)
   ABI_DEALLOCATE(xred_old)
 else
   nullify(rhor,taur,pawrhoij,xred_old)
 end if

!PAW+DMFT
 call destroy_sc_dmft(paw_dmft)
 ! This call should be done inside destroy_sc_dmft
 if ( dtset%usedmft /= 0 ) then
   call data4entropyDMFT_destroy(paw_dmft%forentropyDMFT)
 end if

!Destroy electronpositron datastructure
 if (dtset%positron/=0) then
   call destroy_electronpositron(electronpositron)
 end if

!Deallocating the basis set.
 if (dtset%usewvl == 1) then
   call wvl_projectors_free(wvl%projectors)
   call wvl_wfs_free(wvl%wfs)
   call wvl_descr_free(wvl%descr)
   call wvl_denspot_free(wvl%den)
   if(dtset%usepaw == 1) then
     call wvl_paw_free(psps%ntypat,wvl%descr,wvl%projectors)
   end if
 end if

 ABI_DEALLOCATE(kg)

 if (dtset%icoulomb /= 0) then
   call psolver_kernel((/ 0._dp, 0._dp, 0._dp /), 0, dtset%icoulomb, 0, kernel_dummy, &
&   0, dtset%ngfft, 1, dtset%nscforder)
 end if

 if ((dtset%berryopt<0).or.&
& (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
& dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17)) then
   ABI_DEALLOCATE(pwind)
   ABI_DEALLOCATE(pwnsfac)
   if (xmpi_paral == 1) then
     ABI_DEALLOCATE(mpi_enreg%kptdstrb)
     if (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
&     dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17) then
       ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
     end if
   end if
   if (allocated(mpi_enreg%kpt_loc2ibz_sp))  then
     ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
   end if
   if (allocated(mpi_enreg%kpt_loc2fbz_sp)) then
     ABI_DEALLOCATE(mpi_enreg%kpt_loc2fbz_sp)
   end if
   if (allocated(mpi_enreg%mkmem)) then
     ABI_DEALLOCATE(mpi_enreg%mkmem)
   end if
 else
   ABI_DEALLOCATE(pwind)
   ABI_DEALLOCATE(pwnsfac)
 end if

 call destroy_efield(dtefield)
 call destroy_bfield(dtbfield)

!deallocate Recursion
 if (dtset%userec == 1) then
   call CleanRec(rec_set)
 end if

!Clean the header
 call hdr_free(hdr)

 if (me == 0 .and. dtset%prtxml == 1) then
!  The dataset given in argument has been treated, then we output its variables.
!  call outvarsXML()
!  gstate() will handle a dataset, so we output the dataSet markup.
   write(ab_xml_out, "(A)") '  </dataSet>'
 end if

!Clean the MPI informations
 if (dtset%usewvl == 0.and. dtset%tfkinfunc /= 2) then
!  Plane-wave case
   call bandfft_kpt_destroy_array(bandfft_kpt,mpi_enreg)
 end if

 if(dtset%wfoptalg == 1 .and. psps%usepaw == 1) then
   call destroy_invovl(dtset%nkpt)
 end if
 
 if(gemm_nonlop_use_gemm) then
   call destroy_gemm_nonlop(dtset%nkpt)
   gemm_nonlop_use_gemm = .false.
 end if

!Eventually clean cuda runtime
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call dealloc_hamilt_gpu(2,dtset%use_gpu_cuda)
 end if
#endif

 call timab(36,2,tsec)
 call timab(32,2,tsec)

 DBG_EXIT("COLL")

end subroutine gstate
!!***
