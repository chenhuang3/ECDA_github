!{\src2tex{textfont=tt}}
!!****f* ABINIT/indefo
!! NAME
!! indefo
!!
!! FUNCTION
!! Initialisation phase : defaults values for most input variables
!! (some are initialized earlier, see indefo1 routine)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2014 ABINIT group (XG,MM,FF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!!  nprocs=Number of MPI processors available.
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are given a default value here.
!!   The dataset with number 0 should be the reference default value
!!   in the remaining of the code.
!!
!! NOTES
!! The outputs of this routine are the defaults values of input
!! variables, stored at the index 0 of the last dimension of their
!! multi-dataset representation.
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine indefo(dtsets,ndtset_alloc,nprocs)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_gwdefs
#if defined DEV_YP_VDWXC
 use m_xc_vdw
#endif

 use defs_abitypes,  only : dataset_type
 use m_fftcore,      only : get_cache_kb, fftalg_for_npfft

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'indefo'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset_alloc,nprocs
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc) !vz_i

!Local variables -------------------------------
!scalars
 integer :: idtset,ii,jdtset,paral_atom_default
 logical :: wvl_bigdft
#if defined DEV_YP_VDWXC
 type(xc_vdw_type) :: vdw_defaults
#endif

!******************************************************************

 DBG_ENTER("COLL")

!Set up default values. All variables to be output in outvars.f
!should have a default, even if a nonsensible one can be
!chosen to garantee print in that routine.

!These variables have already been initialized, for idtset/=0
 dtsets(0)%istatr=49
 dtsets(0)%istatshft=1
 dtsets(0)%kptrlatt(1:3,1:3)=0
 dtsets(0)%qptrlatt(1:3,1:3)=0
 dtsets(0)%ptgroupma=0
 dtsets(0)%spgroup=0
 dtsets(0)%shiftk(:,:)=half
 dtsets(0)%tolsym=tol8
 dtsets(0)%znucl(:)=zero
 dtsets(0)%usedmft=0

 paral_atom_default=0
 if (nprocs>1.and.maxval(dtsets(:)%usepaw)>0) paral_atom_default=1

!WARNING : set default in all datasets, including idtset=0 !!!
!Use alphabetic order

 do idtset=0,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset

   wvl_bigdft=.false.
   if(dtsets(idtset)%usewvl==1 .and. dtsets(idtset)%wvl_bigdft_comp==1) then
     wvl_bigdft=.true.
   end if
!  Special case of use_gpu_cuda (can be undertermined at this point)
!  use_gpu_cuda=-1 means undetermined ; here impose its value due to some restrictions
   if (dtsets(idtset)%use_gpu_cuda==-1) then
     if (dtsets(idtset)%optdriver/=0.or.&
&     dtsets(idtset)%tfkinfunc/=0.or.&
&     dtsets(idtset)%nspinor/=1) then
       dtsets(idtset)%use_gpu_cuda=0
     else
       dtsets(idtset)%use_gpu_cuda=1
     end if
   end if

!  A
!  Note that this default value might be overriden for specific datasets later, in case of MPI_IO
   dtsets(idtset)%accesswff=IO_MODE_FORTRAN
   dtsets(idtset)%adpimd=0
   dtsets(idtset)%adpimd_gamma=one
   dtsets(idtset)%accuracy=0
   dtsets(idtset)%atvshift(:,:,:)=zero
   dtsets(idtset)%awtr=1
!  B
   dtsets(idtset)%bdberry(1:4)=0
   dtsets(idtset)%bdeigrf=-1
   dtsets(idtset)%bdgw=0
   dtsets(idtset)%berrystep=1
   dtsets(idtset)%bmass=ten
   dtsets(idtset)%boxcenter(1:3)=half
   dtsets(idtset)%boxcutmin=two
   dtsets(idtset)%brvltt=0
   dtsets(idtset)%bs_nstates=0
!  dtsets(idtset)%bs_hayd_term=0
   dtsets(idtset)%bs_hayd_term=1
   dtsets(idtset)%builtintest=0
   dtsets(idtset)%bxctmindg=two
!  C
   dtsets(idtset)%cd_halfway_freq=3.674930883_dp !(100 eV)
   dtsets(idtset)%cd_imfrqs(:) = zero
   dtsets(idtset)%cd_max_freq=36.74930883_dp     !(1000 eV)
   dtsets(idtset)%cd_subset_freq(1:2)=0
   dtsets(idtset)%cd_frqim_method=1
   dtsets(idtset)%cd_full_grid=0
   dtsets(idtset)%charge=zero
   dtsets(idtset)%chkexit=0
   dtsets(idtset)%chkgwcomp=0
   dtsets(idtset)%chksymbreak=1
   dtsets(idtset)%cineb_start=7
   dtsets(idtset)%corecs(:) = zero
   dtsets(idtset)%cgtyphf = 0
!  D
   dtsets(idtset)%ddamp=0.1_dp
   dtsets(idtset)%delayperm=0
   dtsets(idtset)%diecut=2.2_dp
   dtsets(idtset)%dielng=1.0774841_dp
   dtsets(idtset)%diemac=1.0d6
   if (dtsets(idtset)%usepaw==0) then
     dtsets(idtset)%diemix=one
   else
     dtsets(idtset)%diemix=0.7_dp
   end if
   dtsets(idtset)%diemixmag=dtsets(idtset)%diemix
   dtsets(idtset)%diegap=0.1_dp
   dtsets(idtset)%dielam=half
   dtsets(idtset)%diismemory=8
   dtsets(idtset)%dilatmx=one
!  dtsets(idtset)%dilatmx=1.5 !vz_i
   dtsets(idtset)%dmatpuopt=2
   if (size(dtsets(idtset)%dmatpawu,4)>0) dtsets(idtset)%dmatpawu=-10._dp
   dtsets(idtset)%dmatudiag=0
   dtsets(idtset)%dmft_entropy=0
   dtsets(idtset)%dmft_dc  =1
   dtsets(idtset)%dmft_iter=0
   dtsets(idtset)%dmft_nlambda=6
   dtsets(idtset)%dmft_nwli=0
   dtsets(idtset)%dmft_nwlo=0
   dtsets(idtset)%dmft_mxsf=0.3_dp
   dtsets(idtset)%dmft_read_occnd=0
   dtsets(idtset)%dmft_rslf=0
   dtsets(idtset)%dmft_solv=5
   dtsets(idtset)%dmft_t2g=0
   dtsets(idtset)%dmft_tollc=tol5
   dtsets(idtset)%dmftbandi=0
   dtsets(idtset)%dmftbandf=0
   dtsets(idtset)%dmftcheck=0
   dtsets(idtset)%dmftctqmc_check =0
   dtsets(idtset)%dmftctqmc_correl=0
   dtsets(idtset)%dmftctqmc_grnns =0
   dtsets(idtset)%dmftctqmc_meas  =1
   dtsets(idtset)%dmftctqmc_mrka  =0
   dtsets(idtset)%dmftctqmc_mov   =0
   dtsets(idtset)%dmftctqmc_order =0
   dtsets(idtset)%dmftqmc_l=0
   dtsets(idtset)%dmftqmc_n=0.0_dp
   dtsets(idtset)%dmftqmc_seed=jdtset
   dtsets(idtset)%dmftqmc_therm=1000
   dtsets(idtset)%dmftctqmc_gmove = dtsets(idtset)%dmftqmc_therm / 10
   dtsets(idtset)%dosdeltae=zero
   dtsets(idtset)%dtion=100.0_dp
!  E
   dtsets(idtset)%ecut=-one
   dtsets(idtset)%ecuteps=zero
   dtsets(idtset)%ecutsigx=zero ! The true default value is ecut . This is defined in invars2.F90
   dtsets(idtset)%ecutsm=zero
   dtsets(idtset)%ecutwfn=zero ! The true default value is ecut . This is defined in invars2.F90
   dtsets(idtset)%effmass=one
   dtsets(idtset)%elph2_imagden=zero
   dtsets(idtset)%enunit=0
   dtsets(idtset)%eshift=zero
   dtsets(idtset)%esmear=0.04_dp
   dtsets(idtset)%exchn2n3d=0
   dtsets(idtset)%exchmix=quarter
!  F
   dtsets(idtset)%fermie_nest=zero
   dtsets(idtset)%fftgw=21
   dtsets(idtset)%fft_opt_lob=1;if(dtsets(idtset)%paral_kgb>0)dtsets(idtset)%fft_opt_lob=2
   dtsets(idtset)%freqsusin=one
   dtsets(idtset)%freqsuslo=one
   dtsets(idtset)%optfreqsus=2
   dtsets(idtset)%freqim_alpha=five
   dtsets(idtset)%freqremin=zero
   dtsets(idtset)%freqremax=zero
   dtsets(idtset)%freqspmin=zero
   dtsets(idtset)%freqspmax=zero
   dtsets(idtset)%friction=0.001_dp
   dtsets(idtset)%frzfermi=0
   dtsets(idtset)%fxcartfactor=one ! Should be adjusted to the H2 conversion factor
!  G
   dtsets(idtset)%ga_algor =1
   dtsets(idtset)%ga_fitness =1
   dtsets(idtset)%ga_opt_percent =0.2_dp
   dtsets(idtset)%ga_rules(:) =1
   dtsets(idtset)%goprecon =0
   dtsets(idtset)%goprecprm(:)=0
   dtsets(idtset)%gpu_linalg_limit=2000000
   if (dtsets(idtset)%gw_customnfreqsp/=0) dtsets(idtset)%gw_freqsp(:) = zero
   dtsets(idtset)%gw_nstep =30
   dtsets(idtset)%gwgamma =0
   if ( dtsets(idtset)%gw_nqlwl > 0 ) then
     dtsets(idtset)%gw_qlwl(:,:)=zero
     dtsets(idtset)%gw_qlwl(1,1)=0.00001_dp
     dtsets(idtset)%gw_qlwl(2,1)=0.00002_dp
     dtsets(idtset)%gw_qlwl(3,1)=0.00003_dp
   end if
   dtsets(idtset)%gw_eet=-1
   dtsets(idtset)%gw_eet_nband=-1
   dtsets(idtset)%gw_eet_inclvkb=0
   dtsets(idtset)%gw_eet_scale=0.0_dp
   dtsets(idtset)%gw_frqim_inzgrid=0
   dtsets(idtset)%gw_frqre_inzgrid=0
   dtsets(idtset)%gw_frqre_tangrid=0
   dtsets(idtset)%gw_qprange=0
   dtsets(idtset)%gw_npoles=0
   dtsets(idtset)%gw_reconst_scr=0
   dtsets(idtset)%gw_sigxcore=0
   dtsets(idtset)%gw_sctype =  GWSC_one_shot
   dtsets(idtset)%gw_toldfeig=0.1/Ha_eV
   dtsets(idtset)%gw_use_pole_scr=0
   dtsets(idtset)%getbseig=0
   dtsets(idtset)%getbsreso=0
   dtsets(idtset)%getbscoup=0
   dtsets(idtset)%getcell =0
   dtsets(idtset)%getddk  =0
   dtsets(idtset)%getden  =0
   dtsets(idtset)%getgam_eig2nkq  =0
   dtsets(idtset)%gethaydock=0
   dtsets(idtset)%getkss  =0
   dtsets(idtset)%getocc  =0
   dtsets(idtset)%getpawden=0
   dtsets(idtset)%getqps  =0
   dtsets(idtset)%getscr  =0
   dtsets(idtset)%getsuscep=0
   dtsets(idtset)%getvel  =0
   dtsets(idtset)%getwfk  =0
   dtsets(idtset)%getwfkfine = 0
   dtsets(idtset)%getwfq  =0
   dtsets(idtset)%getxcart=0
   dtsets(idtset)%getxred =0
   dtsets(idtset)%get1den =0
   dtsets(idtset)%get1wf  =0
   dtsets(idtset)%gwcalctyp=0
   dtsets(idtset)%gwcomp=0
   dtsets(idtset)%gwencomp=2.0_dp
   dtsets(idtset)%gwmem=11
   dtsets(idtset)%gwpara=1
   dtsets(idtset)%gwrpacorr=0

!  I
   if(dtsets(idtset)%natsph/=0) then
!    do not use iatsph(:) but explicit boundaries
!    to avoid to read to far away in the built array (/ ... /)
     dtsets(idtset)%iatsph(1:dtsets(idtset)%natsph)=(/ (ii,ii=1,dtsets(idtset)%natsph) /)
   else
     dtsets(idtset)%iatsph(:)=0
   end if
   dtsets(idtset)%iboxcut=0
   dtsets(idtset)%idyson=1
   dtsets(idtset)%icutcoul=6
   dtsets(idtset)%ieig2rf=0
   dtsets(idtset)%ikhxc=0
   dtsets(idtset)%inclvkb=1
   dtsets(idtset)%intexact=0
   dtsets(idtset)%intxc=0
!  if (dtsets(idtset)%paral_kgb>0.and.idtset>0) dtsets(idtset)%intxc=0
   dtsets(idtset)%ionmov=0
   dtsets(idtset)%iprcch=2
   if (dtsets(idtset)%paral_kgb>0.and.idtset>0) dtsets(idtset)%iprcch=6 ! Recommended for band-FFT parallelism
#ifdef HAVE_LOTF
   if (dtsets(idtset)%ionmov==23) dtsets(idtset)%iprcch=2 ! Recommended for LOTF
#endif
   dtsets(idtset)%iprcel=0
   dtsets(idtset)%iprctfvw=0
   dtsets(idtset)%iprcfc=0
   dtsets(idtset)%irandom=3
   dtsets(idtset)%irdbseig=0
   dtsets(idtset)%irdbsreso=0
   dtsets(idtset)%irdbscoup=0
   dtsets(idtset)%irdddk=0
   dtsets(idtset)%irdden=0
   dtsets(idtset)%irdhaydock=0
   dtsets(idtset)%irdkss=0
   dtsets(idtset)%irdpawden=0
   dtsets(idtset)%irdqps=0
   dtsets(idtset)%irdscr=0
   dtsets(idtset)%irdsuscep=0
   dtsets(idtset)%irdvdw=0
   dtsets(idtset)%irdwfk=0
   dtsets(idtset)%irdwfkfine=0
   dtsets(idtset)%irdwfq=0
   dtsets(idtset)%ird1den=0
   dtsets(idtset)%ird1wf=0
!iscf
   if(wvl_bigdft) then
     dtsets(idtset)%iscf=0
   else
     if(dtsets(idtset)%usepaw==0) then
       dtsets(idtset)%iscf=7
     else
       dtsets(idtset)%iscf=17
     end if
   end if
   dtsets(idtset)%isecur=0
   dtsets(idtset)%istatimg = 1
   dtsets(idtset)%istwfk(:)=0
   dtsets(idtset)%ixc=1
   dtsets(idtset)%ixcpositron=1
!  J
   dtsets(idtset)%f4of2_sla(:)=-one
   dtsets(idtset)%f6of2_sla(:)=-one
   dtsets(idtset)%jpawu(:)=zero
!  K
!  We should not do that, kberry can be smaller than 1:20
!  dtsets(idtset)%kberry(1:3,1:20)=0
   dtsets(idtset)%kberry(1:3,:)=0
   dtsets(idtset)%kpt(:,:)=zero
   dtsets(idtset)%kptgw(:,:)=zero
   dtsets(idtset)%kptnrm=one
   dtsets(idtset)%kptopt=1
   dtsets(idtset)%kptrlen=30.0_dp
   dtsets(idtset)%kssform=1
!  dtsets(idtset)%kssform=3
!  L
   dtsets(idtset)%ldgapp=0
   dtsets(idtset)%localrdwf=1
#if defined HAVE_LOTF
   dtsets(idtset)%lotf_classic=5
   dtsets(idtset)%lotf_nitex=10  
   dtsets(idtset)%lotf_nneigx=40 
   dtsets(idtset)%lotf_version=2
#endif
!  M
   dtsets(idtset)%magconon = 0
   dtsets(idtset)%magcon_lambda = 10.0_dp
   dtsets(idtset)%max_ncpus = 0
   dtsets(idtset)%mband = -1
   dtsets(idtset)%mdf_epsinf = zero
   dtsets(idtset)%mdtemp(:)=300.0_dp
   dtsets(idtset)%mdwall=10000_dp
   dtsets(idtset)%mep_mxstep=100._dp
   dtsets(idtset)%mep_solver=0
   dtsets(idtset)%mffmem=1
   dtsets(idtset)%mgfft = -1
   dtsets(idtset)%mgfftdg = -1
   dtsets(idtset)%mpw = -1
   dtsets(idtset)%mqgrid=0
   dtsets(idtset)%mqgriddg=0
!  N
   dtsets(idtset)%natrd = -1
   dtsets(idtset)%nband(:)=0
   dtsets(idtset)%nbandhf=0
   dtsets(idtset)%nbandsus=-1
   dtsets(idtset)%nbdblock=1
   dtsets(idtset)%nbdbuf=0
   dtsets(idtset)%nberry=1
   dtsets(idtset)%nbandkss=0
   dtsets(idtset)%nctime=0
   dtsets(idtset)%ndtset = -1
   dtsets(idtset)%neb_algo=1
   dtsets(idtset)%neb_spring(1:2)=(/0.05_dp,0.05_dp/)
   dtsets(idtset)%npwkss=0
   dtsets(idtset)%ndyson=-1
   dtsets(idtset)%nfft = -1
   dtsets(idtset)%nfftdg = -1

   dtsets(idtset)%nfreqim=-1
   dtsets(idtset)%nfreqre=-1
   dtsets(idtset)%nfreqsp=0
   dtsets(idtset)%nfreqsus=0

   dtsets(idtset)%npulayit=7

!  ngfft is a special case
   dtsets(idtset)%ngfft(1:8)=0
   dtsets(idtset)%ngfft(7) = fftalg_for_npfft(1)
!  fftcache=ngfft(8) is machine-dependent.
   dtsets(idtset)%ngfft(8) = get_cache_kb()

   dtsets(idtset)%ngfftdg(:)=dtsets(idtset)%ngfft(:)
!
   !nline
   dtsets(idtset)%nline=4
   if(dtsets(idtset)%usewvl==1 .and. .not. wvl_bigdft) then
     if(dtsets(idtset)%usepaw==1) then
       dtsets(idtset)%nline=4
     else
       dtsets(idtset)%nline=2
     end if
   end if

!  nloalg is also a special case
   dtsets(idtset)%nloalg(1)=4
   dtsets(idtset)%nloalg(2)=4                      ! governs, in fine, the use of registers in the most cpu
   dtsets(idtset)%nloalg(3)=199                    ! sets the size of blocks of planewaves
   dtsets(idtset)%nloalg(4)=10                     ! mincat
   dtsets(idtset)%nloalg(5)=dtsets(idtset)%usepaw
#if defined i386
!  nloalg(2) is machine-dependent, the default is 4, but 8 is better for the P6 and hp .
   dtsets(idtset)%nloalg(2)=8
#elif defined HAVE_FFT_ASL
   dtsets(idtset)%nloalg(1)=2    ! Important to use opernl2 for this vector machine
#endif
   dtsets(idtset)%ngkpt=0
   dtsets(idtset)%nnsclo=0
   dtsets(idtset)%nnsclohf=0
   dtsets(idtset)%nomegasf=100
   dtsets(idtset)%nomegasrd=9
   dtsets(idtset)%nomegasi=12
   dtsets(idtset)%noseinert=1.0d5
   dtsets(idtset)%npvel=0
   dtsets(idtset)%npweps=0
   dtsets(idtset)%npwsigx=0
   dtsets(idtset)%npwwfn=0
   dtsets(idtset)%nqpt=0
   dtsets(idtset)%nscforder=16
   dtsets(idtset)%nsheps=0
   dtsets(idtset)%nshiftk=1
   dtsets(idtset)%nshiftk_orig=0
   dtsets(idtset)%nshsigx=0
   dtsets(idtset)%nshwfn=0
   dtsets(idtset)%nstep=30
   dtsets(idtset)%ntime=1
   dtsets(idtset)%nwfshist=0
   if(dtsets(idtset)%usewvl==1 .and. .not. wvl_bigdft) then
     if(dtsets(idtset)%usepaw==1) then
       dtsets(idtset)%nwfshist=4
     else
       dtsets(idtset)%nwfshist=2
     end if
   end if
!  O
   dtsets(idtset)%occopt=1
   dtsets(idtset)%occ_orig(:)=zero
   dtsets(idtset)%omegasrdmax=1.0_dp/Ha_eV  ! = 1eV
   dtsets(idtset)%omegasimax=50/Ha_eV
   dtsets(idtset)%orbmag=0
   dtsets(idtset)%optcell=0
   dtsets(idtset)%optforces=2
   dtsets(idtset)%optstress=1
   dtsets(idtset)%optnlxccc=1
   if (dtsets(idtset)%usepaw==0) then
     dtsets(idtset)%ortalg=2
!    dtsets(idtset)%ortalg=999
   else
     dtsets(idtset)%ortalg=-2
!    dtsets(idtset)%ortalg=999
   end if
!  P
   dtsets(idtset)%paral_atom=paral_atom_default
   dtsets(idtset)%pawcpxocc=1
   dtsets(idtset)%pawcross=0
   dtsets(idtset)%pawecutdg=-one
   dtsets(idtset)%pawfatbnd=0
   dtsets(idtset)%pawlcutd=10
   dtsets(idtset)%pawlmix=10
   dtsets(idtset)%pawmixdg=0 ! Will be set to 1 when npfft>1
   dtsets(idtset)%pawnhatxc=1
   dtsets(idtset)%pawntheta=12
   dtsets(idtset)%pawnphi=13
   dtsets(idtset)%pawnzlm=1
   dtsets(idtset)%pawoptmix=0
   dtsets(idtset)%pawoptosc=0
   dtsets(idtset)%pawovlp=5._dp
   dtsets(idtset)%pawprtdos=0
   dtsets(idtset)%pawprtvol=0
   dtsets(idtset)%pawprtwf=0
   dtsets(idtset)%pawprt_k=0
   dtsets(idtset)%pawprt_b=0
   dtsets(idtset)%pawstgylm=1
   dtsets(idtset)%pawsushat=0
   dtsets(idtset)%pawujat=1
   dtsets(idtset)%pawujrad=20.0_dp
   dtsets(idtset)%pawujv=0.1_dp/Ha_eV
   dtsets(idtset)%pawusecp=1
   dtsets(idtset)%pawxcdev=1
   dtsets(idtset)%pitransform=0
   dtsets(idtset)%ptcharge(:) = zero
   dtsets(idtset)%plowan_bandi=0
   dtsets(idtset)%plowan_bandf=0
   if(dtsets(idtset)%plowan_compute>0) then
     dtsets(idtset)%plowan_it(:)=0
     dtsets(idtset)%plowan_iatom(:)=0
     dtsets(idtset)%plowan_lcalc(:)=-1
     dtsets(idtset)%plowan_projcalc(:)=0
     dtsets(idtset)%plowan_nbl(:)=0
   end if
   dtsets(idtset)%plowan_natom=0
   dtsets(idtset)%plowan_nt=0
   dtsets(idtset)%plowan_realspace=0
   dtsets(idtset)%pol(:)=zero
   dtsets(idtset)%polcen(:)=zero
   dtsets(idtset)%posdoppler=0
   dtsets(idtset)%positron=0
   dtsets(idtset)%posnstep=50
   dtsets(idtset)%posocc=one
   dtsets(idtset)%postoldfe=0.000001_dp
   dtsets(idtset)%postoldff=zero
   dtsets(idtset)%ppmodel=1
   dtsets(idtset)%ppmfrq=zero
   dtsets(idtset)%prepanl=0
   dtsets(idtset)%prepgkk=0
   dtsets(idtset)%prtbbb=0
   dtsets(idtset)%prtbltztrp=0
   dtsets(idtset)%prtcif=0
   dtsets(idtset)%prtcml=0
   dtsets(idtset)%prtcs=0
   dtsets(idtset)%prtden=1;if (dtsets(idtset)%nimage>1) dtsets(idtset)%prtden=0
   dtsets(idtset)%prtdensph=0
   dtsets(idtset)%prtdipole=0
   dtsets(idtset)%prtdos=0
   dtsets(idtset)%prtdosm=0
   dtsets(idtset)%prtefg=0
   dtsets(idtset)%prteig=1;if (dtsets(idtset)%nimage>1) dtsets(idtset)%prteig=0
   dtsets(idtset)%prtelf=0
   dtsets(idtset)%prtfc=0
   dtsets(idtset)%prtfsurf=0
   dtsets(idtset)%prtgden=0
   dtsets(idtset)%prtgeo=0
   dtsets(idtset)%prtgkk=0
   dtsets(idtset)%prtkden=0
   dtsets(idtset)%prtkpt = -1
   dtsets(idtset)%prtlden=0
   dtsets(idtset)%prtnabla=0
   dtsets(idtset)%prtnest=0
   dtsets(idtset)%prtposcar=0
   dtsets(idtset)%prtpot=0
   dtsets(idtset)%prtspcur=0
   dtsets(idtset)%prtsuscep=1
   dtsets(idtset)%prtstm=0
   dtsets(idtset)%prtvdw=0
   dtsets(idtset)%prtvha=0
   dtsets(idtset)%prtvhxc=0
   dtsets(idtset)%prtvxc=0
   dtsets(idtset)%prtvol=0
   dtsets(idtset)%prtvolimg=0
   dtsets(idtset)%prtwant=0
   dtsets(idtset)%prtwf=1;;if (dtsets(idtset)%nimage>1) dtsets(idtset)%prtwf=0
   dtsets(idtset)%prtxml = 0
   do ii=1,dtsets(idtset)%natom,1
     dtsets(idtset)%prtatlist(ii)=ii
   end do
   dtsets(idtset)%prt1dm=0
   dtsets(idtset)%pvelmax(:)=one 
!  Q
   dtsets(idtset)%qmass(:)=ten
   dtsets(idtset)%qprtrb(1:3)=0
   dtsets(idtset)%qptdm(:,:)=zero
   dtsets(idtset)%quadmom(:) = zero
!  R
   dtsets(idtset)%random_atpos=zero
   dtsets(idtset)%ratsph_extra=two
   dtsets(idtset)%recefermi=zero
   dtsets(idtset)%recgratio=1
   dtsets(idtset)%recnpath=500
   dtsets(idtset)%recnrec=10
   dtsets(idtset)%recrcut=zero
   dtsets(idtset)%recptrott=0
   dtsets(idtset)%rectesteg=0
   dtsets(idtset)%rectolden=zero
   dtsets(idtset)%rcut=zero
   dtsets(idtset)%restartxf=0
   dtsets(idtset)%rfasr=0
   dtsets(idtset)%rfatpol(1:2)=1
   dtsets(idtset)%rfddk=0
   dtsets(idtset)%rfdir(1:3)=0
   dtsets(idtset)%rfelfd=0
   dtsets(idtset)%rfmeth=1
   dtsets(idtset)%rfphon=0
   dtsets(idtset)%rfstrs=0
   dtsets(idtset)%rfuser=0

   dtsets(idtset)%rf1atpol(1:2)=1
   dtsets(idtset)%rf1dir(1:3)=0
   dtsets(idtset)%rf1elfd=0
   dtsets(idtset)%rf1phon=0

   dtsets(idtset)%rf2atpol(1:2)=1
   dtsets(idtset)%rf2dir(1:3)=0
   dtsets(idtset)%rf2elfd=0
   dtsets(idtset)%rf2phon=0

   dtsets(idtset)%rf3atpol(1:2)=1
   dtsets(idtset)%rf3dir(1:3)=0
   dtsets(idtset)%rf3elfd=0
   dtsets(idtset)%rf3phon=0

   dtsets(idtset)%rhoqpmix=one

!  S
   dtsets(idtset)%sciss=zero
   dtsets(idtset)%shiftk_orig(:,:)=one
   dtsets(idtset)%signperm=1
   dtsets(idtset)%slabwsrad=zero
   dtsets(idtset)%smdelta=0
   dtsets(idtset)%soenergy=zero
   dtsets(idtset)%spbroad=0.1
   dtsets(idtset)%spgaxor = -1
   dtsets(idtset)%spgorig = -1
   dtsets(idtset)%spinmagntarget=-99.99_dp
   dtsets(idtset)%spmeth=0
   dtsets(idtset)%spnorbscl=one
   dtsets(idtset)%stmbias=zero
   dtsets(idtset)%strfact=100.0_dp
   dtsets(idtset)%string_algo=1
   dtsets(idtset)%strprecon=one
   dtsets(idtset)%strtarget(1:6)=zero
   dtsets(idtset)%supercell(:)=1
   dtsets(idtset)%suskxcrs=0
   dtsets(idtset)%symchi=0
   dtsets(idtset)%symsigma=0
!  T
   dtsets(idtset)%td_maxene=zero
   dtsets(idtset)%td_mexcit=0
   dtsets(idtset)%tl_nprccg = 30
   dtsets(idtset)%tl_radius = zero
   dtsets(idtset)%tphysel=zero
   dtsets(idtset)%toldfe=zero
   dtsets(idtset)%toldff=zero
   dtsets(idtset)%tolimg=5.0d-5
   dtsets(idtset)%tolrde=0.005_dp
   dtsets(idtset)%tolrff=zero
   dtsets(idtset)%tolmxf=5.0d-5
   dtsets(idtset)%tolvrs=zero
   dtsets(idtset)%tolwfr=zero
   dtsets(idtset)%tsmear=0.04_dp
!  U
   dtsets(idtset)%ucrpa=0
   dtsets(idtset)%ucrpa_bands(:)=-1
   dtsets(idtset)%ucrpa_window(:)=-1.0_dp
   dtsets(idtset)%upawu(:)=zero
   dtsets(idtset)%usefock=0
   dtsets(idtset)%usekden=0
   dtsets(idtset)%use_gemm_nonlop=0
   dtsets(idtset)%use_nonscf_gkk=0 !1 ! deactivate by default, for now 6 Oct 2013
   dtsets(idtset)%userec=0
   dtsets(idtset)%usexcnhat_orig=-1
   dtsets(idtset)%useylm=0
!  V
   dtsets(idtset)%vacnum = -1
   dtsets(idtset)%vcutgeo(:)=zero
   dtsets(idtset)%vdw_nfrag = 1
#if defined DEV_YP_VDWXC
   dtsets(idtset)%vdw_df_acutmin = vdw_defaults%acutmin
   dtsets(idtset)%vdw_df_aratio = vdw_defaults%aratio
   dtsets(idtset)%vdw_df_damax = vdw_defaults%damax
   dtsets(idtset)%vdw_df_damin = vdw_defaults%damin
   dtsets(idtset)%vdw_df_dcut = vdw_defaults%dcut
   dtsets(idtset)%vdw_df_dratio = vdw_defaults%dratio
   dtsets(idtset)%vdw_df_dsoft = vdw_defaults%dsoft
   dtsets(idtset)%vdw_df_gcut = vdw_defaults%gcut
   dtsets(idtset)%vdw_df_ndpts = vdw_defaults%ndpts
   dtsets(idtset)%vdw_df_ngpts = vdw_defaults%ngpts
   dtsets(idtset)%vdw_df_nqpts = vdw_defaults%nqpts
   dtsets(idtset)%vdw_df_nrpts = vdw_defaults%nrpts
   dtsets(idtset)%vdw_df_nsmooth = vdw_defaults%nsmooth
   dtsets(idtset)%vdw_df_phisoft = vdw_defaults%phisoft
   dtsets(idtset)%vdw_df_qcut = vdw_defaults%qcut
   dtsets(idtset)%vdw_df_qratio = vdw_defaults%qratio
   dtsets(idtset)%vdw_df_rcut = vdw_defaults%rcut
   dtsets(idtset)%vdw_df_rsoft = vdw_defaults%rsoft
   dtsets(idtset)%vdw_df_tolerance = vdw_defaults%tolerance
   dtsets(idtset)%vdw_df_tweaks = vdw_defaults%tweaks
   dtsets(idtset)%vdw_df_zab = vdw_defaults%zab
   dtsets(idtset)%vdw_df_threshold = 1.0d-2
#endif
   dtsets(idtset)%vdw_supercell(:) = 0
   dtsets(idtset)%vdw_tol = tol10
   dtsets(idtset)%vdw_typfrag(:) = 1
   dtsets(idtset)%vdw_xc = 0
   dtsets(idtset)%vis=100.0_dp
   dtsets(idtset)%vprtrb(1:2)=zero
!  W
   dtsets(idtset)%wtatcon(:,:,:)=zero
   dtsets(idtset)%wtk=one
   dtsets(idtset)%wvl_crmult  = 6._dp
   dtsets(idtset)%wvl_frmult  = 10._dp
   dtsets(idtset)%wvl_hgrid   = 0.5_dp
   dtsets(idtset)%wvl_ngauss  =(1,100)
   dtsets(idtset)%wvl_nprccg  = 10
   dtsets(idtset)%w90iniprj   = 1
   dtsets(idtset)%w90prtunk   = 0

!  X
   dtsets(idtset)%xclevel  = 0
   dtsets(idtset)%xc_denpos = tol14
   dtsets(idtset)%xc_tb09_c = 99.99_dp
!  Y
!  Z
   dtsets(idtset)%zcut=3.67493260d-03  ! = 0.1eV
   dtsets(idtset)%ziontypat(:)=zero

!  BEGIN VARIABLES FOR @Bethe-Salpeter
   dtsets(idtset)%bs_algorithm    =1
   dtsets(idtset)%bs_haydock_niter=100
   dtsets(idtset)%bs_exchange_term=1
   dtsets(idtset)%bs_coulomb_term=11
   dtsets(idtset)%bs_calctype=1
   dtsets(idtset)%bs_coupling=1

   dtsets(idtset)%bs_haydock_tol=(0.02_dp,zero)

   dtsets(idtset)%bs_loband=0
!  Take big absolute value numbers, but the the biggest ones, otherwise overflow can happen
   dtsets(idtset)%bs_eh_cutoff = [smallest_real*tol6,greatest_real*tol6]
   dtsets(idtset)%bs_freq_mesh = [zero,zero,0.01_dp/Ha_eV]

!  Interpolation
   dtsets(idtset)%bs_interp_mode=0 ! No interpolation
   dtsets(idtset)%bs_interp_prep=0 ! Do not prepare interp
   dtsets(idtset)%bs_interp_kmult=(/zero,zero,zero/)

!  END VARIABLES FOR @Bethe-Salpeter.

! JB:UNINITIALIZED VALUES (not found in this file neither indefo1)
! They might be initialized somewhereelse, I don't know.
! That might cause unitialized error with valgrind depending on the compilo
! chkprim
! maxnsym
! nsym
! macro_uj
! prtpmp
! timopt
! useria
! userib
! useric
! userid
! userie
! bravais
! symafm
! symrel
! fband
! nelect
! userra
! userrb
! userrc
! userrd
! userre
! vacwidth
! genafm
! kptns
! rprimd_orig
! tnons
! xredsph_extra

 end do

 DBG_EXIT("COLL")

end subroutine indefo
!!***
