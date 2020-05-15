! {\src2tex{textfont=tt}}
!!****f* ABINIT/outscfcv
!! NAME
!! outscfcv
!!
!! FUNCTION
!! Output routine for the scfcv.F90 routine
!!
!! COPYRIGHT
!! Copyright (C) 2005-2014 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cg(2,mcg)=planewave coefficients of wavefunctions (see also side effects)
!!  compch_fft=compensation charge, from FFT grid
!!  compch_sph=compensation charge, from sphere
!!  cprj(natom,mcprj*usecprj)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!          and each |p_lmn> non-local projector. See also side effects
!!  dimcprj(natom*usecprj)=array of dimensions of array cprj (not ordered)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  elfr(nfft,nspden(+1))=electron localization function, real space.
!!   (+1) if spin-polarized in order to get total, spin up and spin down elf
!!  etotal=total energy
!!  fermie= Fermi energy
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  grhor(nfft,nspden,3)= gradient of electron density in electrons/bohr**4, real space
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  lrhor(nfft,nspden)= Laplacian of electron density in electrons/bohr**5, real space
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw**mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfftc=maximum size of 1D FFTs for the PAW coarse grid
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nhat(nfft,nspden*usepaw)= compensation charge density  (PAW)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms in unit cell.
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  occ(mband*nkpt*nsppol)=occupation number for each band (usually 2) for each k.
!!  paw_an(my_natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(natom) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)> tables on PAW fine grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  paw_ij(my_natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!            note:structure factors are given on the coarse grid for PAW
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!     forces and its components, the stress tensor) of a ground-state computation
!!  rhor(nfft,nspden)=total electron density in electrons/bohr**3, real space.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  taur(nfft,nspden)=total kinetic energy density in bohr**(-5), real space.
!!  ucvol=unit cell volume (bohr**3)
!!  usecprj=1 if cprj datastructure has been allocated
!!  vhartr(nfft)=Hartree potential
!!  vxc(nfft,nspden)=xc potential
!!  vtrial(nfft,nspden)=the trial potential
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  If prtwant==3 the following quantitities are updated using the unitary transformation
!!  defining the QP amplitudes in terms of the KS basis set:
!!   cg(2,mcg)=planewave coefficients of wavefunctions.
!!   cprj(natom,mcprj*usecpyj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      abi_etsf_electrons_put,abi_etsf_geo_put,bonds_lgth_angles,calc_cs
!!      calc_efg,calc_fc,calcdensph,crystal_free,crystal_init,crystal_ncwrite
!!      datafordmft,denfgr,destroy_dmft,destroy_oper,dos_degeneratewfs
!!      ebands_free,ebands_init,ebands_ncwrite,gaus_dos,init_dmft,init_oper
!!      ioarr,mag_constr_e,mknesting,mlwfovlp,mlwfovlp_qp,multipoles_out
!!      optics_paw,optics_paw_core,optics_vloc,out1dm,outkss,outwant
!!      partial_dos_fractions,partial_dos_fractions_paw,pawmkaewf,pawprt
!!      pawrhoij_copy,pawrhoij_nullify,posdoppler,poslifetime,print_dmft
!!      prt_cif,prt_cml2,prtbltztrp_out,prtfatbands,read_atomden
!!      results_gs_ncwrite,tetrahedron,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dtefield,dtfil,dtset,&
& ecut,eigen,electronpositron,elfr,etotal,fermie,gmet,gprimd,grhor,hdr,kg,&
& lrhor,mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpsang,mpw,my_natom,natom,&
& nattyp,nfft,ngfft,nhat,nkpt,npwarr,nspden,nsppol,nsym,ntypat,n3xccc,occ,&
& paw_dmft,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,paw_an,paw_ij,prtvol,psps,results_gs,rhor,rprimd,&
& taur,ucvol,usecprj,vhartr,vtrial,vxc,wvl_den,xccc3d,xred)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
 use defs_abitypes
 use defs_parameters
 use m_profiling_abi
 use m_wffile
 use m_efield
 use m_errors
 use m_xmpi
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif
 use m_abi_etsf
 use m_ncfile
 use m_header

 use m_io_tools,         only : get_unit, open_file
 use m_fstrings,         only : strcat
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use m_oper,             only : oper_type,init_oper,destroy_oper
 use m_crystal,          only : crystal_init, crystal_free, crystal_t, isalchemical
 use m_crystal_io,       only : crystal_ncwrite
 use m_ebands,           only : ebands_init, ebands_ncwrite, ebands_free, ebands_3dprint, get_nelect_per_spin
 use m_results_gs,       only : results_gs_type, results_gs_ncwrite
 use m_nesting,          only : mknesting
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type, simp_gen
 use m_pawtab,           only : pawtab_type
 use m_paw_an,           only : paw_an_type
 use m_paw_ij,           only : paw_ij_type
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_pawrhoij,         only : pawrhoij_type, pawrhoij_nullify, pawrhoij_copy
 use m_pawcprj,          only : pawcprj_type
 use m_pawfgr,           only : pawfgr_type
 use m_paw_dmft,         only : paw_dmft_type,init_dmft,destroy_dmft,print_dmft
 use m_plowannier

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outscfcv'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_47_xml
 use interfaces_53_spacepar
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_68_dmft
 use interfaces_69_wfdesc
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mgfftc,mkmem,mpsang,mpw,n3xccc,my_natom,natom,nfft
 integer,intent(in) :: nkpt,nspden,nsppol,nsym,ntypat,prtvol,usecprj
 real(dp),intent(in) :: compch_fft,compch_sph,ecut,fermie,ucvol
 real(dp),intent(inout) :: etotal
 type(efield_type),intent(in) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(results_gs_type),intent(in) :: results_gs
 type(wvl_denspot_type), intent(in) :: wvl_den
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom*usecprj)
 integer,intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat),ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: rprimd(3,3),vhartr(nfft),xccc3d(n3xccc)
 real(dp),intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: nhat(nfft,nspden*psps%usepaw)
 real(dp),intent(inout) :: rhor(nfft,nspden),vtrial(nfft,nspden)
 real(dp),intent(inout) :: vxc(nfft,nspden),xred(3,natom)
 real(dp),pointer :: elfr(:,:),grhor(:,:,:),lrhor(:,:),taur(:,:)
 type(pawcprj_type),intent(inout) :: cprj(natom,mcprj*usecprj)
 type(paw_an_type),intent(inout) :: paw_an(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(my_natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),target,intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: bantot !ncerr,
 integer :: accessfil,coordn,fatbands_flag,fformr,fformv
 integer :: ii,ierr,ifft,ikpt,ispden,isppol
 integer :: prtdosm,mbesslang,me,my_natom_tmp,ndosfraction
 integer :: occopt,partial_dos_flag,paw_dos_flag,pawfatbnd,prt1dm
 integer :: prtcml,prtcs,prtden,prtdos,prtefg,prtelf,prtfc,prtgden,prtgeo,prtkden,prtlden,prtnabla
 integer :: prtpot,prtstm,prtvha,prtvhxc,prtvxc,rdwr,rdwrpaw,pawprtden
 integer :: nqpath, skipnest,iband,nocc,spacecomm,spacecomm_fft,tmp_unt
 real(dp) :: gammastate_core,norm,occ_norm,unocc_norm
 real(dp) :: invgauwidth,prefact,rate_dum
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 integer :: unitmatrix(9)
 real(dp) :: tsec(2),nt_ntone_norm(nspden),nelect_per_spin(nsppol)
 real(dp),allocatable :: dos_fractions(:,:,:,:),dos_fractions_m(:,:,:,:),dos_fractions_average_m(:,:,:,:)
 real(dp),allocatable :: dos_fractions_paw1(:,:,:,:)
 real(dp),allocatable :: dos_fractions_pawt1(:,:,:,:),eigen2(:)
 real(dp),allocatable :: elfr_down(:,:),elfr_up(:,:)
 real(dp),allocatable :: rhor_paw(:,:),rhor_paw_core(:,:),rhor_paw_val(:,:),vwork(:,:)
 real(dp),allocatable :: rhor_n_one(:,:),rhor_nt_one(:,:),ps_norms(:,:,:)
 real(dp), allocatable :: qpath_vertices(:,:)
 real(dp), allocatable :: fs_weights(:,:,:),doccde(:)
 type(pawrhoij_type) :: pawrhoij_dum(0)
 type(pawrhoij_type),pointer :: pawrhoij_all(:)
 logical :: remove_inv
 real(dp) :: e_fermie
 type(oper_type) :: lda_occup
 type(crystal_t) :: Crystal
 type(ebands_t) :: Bands
 type(ncfile_t) :: ncf
 type(plowannier_type) :: wan
! *************************************************************************

 DBG_ENTER("COLL")
 call timab(950,1,tsec) ! outscfcv




 
 if ((usecprj==0.or.mcprj==0).and.psps%usepaw==1.and. &
& (dtset%prtwant==2.or.dtset%prtwant==3.or.dtset%prtnabla>0.or.dtset%prtdos==3 &
& .or.dtset%kssform==3.or.dtset%pawfatbnd>0.or.dtset%pawprtwf>0)) then
   write (message,'(5a)')&
&   'cprj datastructure must be allocated',ch10,&
&   'with options prtwant=2,3, prtnabla>0, prtdos>3, kssform==3, pawfatbnd>0, pawprtwf>0',ch10,&
&   'Action: change pawusecp input keyword.'
   MSG_ERROR(message)
 end if

!Initialize two objects to facilitate the propagation of info.
!These objects should used more frequently, actually they should
!become basic objects used in abinit.

!Crystalline structure.
 remove_inv=.false.
 if(dtset%nspden==4 .and. dtset%usedmft==1) remove_inv=.true. ! MG: why this?

 call crystal_init(Crystal,dtset%spgroup,natom,dtset%npsp,ntypat, &
& dtset%nsym,rprimd,dtset%typat,xred,dtset%ziontypat,dtset%znucl,1,&
& dtset%nspden==2.and.dtset%nsppol==1,remove_inv,hdr%title,&
& dtset%symrel,dtset%tnons,dtset%symafm)

!Electron band energies.
 bantot= dtset%mband*dtset%nkpt*dtset%nsppol
 ABI_MALLOC(doccde,(bantot))
 doccde=zero

 call ebands_init(bantot,Bands,dtset%nelect,doccde,eigen,hdr%istwfk,hdr%kptns,hdr%nband,&
& hdr%nkpt,hdr%npwarr,hdr%nsppol,hdr%nspinor,hdr%tphysel,hdr%tsmear,hdr%occopt,hdr%occ,hdr%wtk)

 ABI_FREE(doccde)

 Bands%fermie  = results_gs%energies%e_fermie
 e_fermie=results_gs%energies%e_fermie
 Bands%entropy = results_gs%energies%entropy

!wannier interface
 call timab(951,1,tsec)
 if (dtset%prtwant==2) then

   call mlwfovlp(atindx1,cg,cprj,dtset,dtfil,eigen,gprimd,hdr,kg,&
&   mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpw,natom,&
&   nattyp,nfft,ngfft,nkpt,npwarr,nsppol,ntypat,occ,&
&   pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)

 else if (dtset%prtwant==3) then

!  Convert cg and eigen to GW quasiparticle wave functions and eigenvalues in mlwfovlp_qp
   ABI_ALLOCATE(eigen2,(mband*nkpt*nsppol))
   eigen2=eigen

   call mlwfovlp_qp(cg,cprj,dtset,dtfil,eigen2,mband,mcg,mcprj,mkmem,mpw,natom,&
&   nkpt,npwarr,nspden,nsppol,ntypat,Hdr,pawtab,rprimd,MPI_enreg)

!  Call Wannier90
   call mlwfovlp(atindx1,cg,cprj,dtset,dtfil,eigen2,gprimd,hdr,kg,&
&   mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpw,natom,&
&   nattyp,nfft,ngfft,nkpt,npwarr,nsppol,ntypat,occ,&
&   pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)

!  this is the old implementation, risky due to unpredictable size effects
!  now eigen is not overwritten, one should use other ways to print the GW corrections
!  eigen=eigen2
   ABI_DEALLOCATE(eigen2)
 end if !prtwant
 call timab(951,2,tsec)

!
!if accesswff == 2 then set all outputs to netcdf format
!if accesswff == 3 then set all outputs to ETSF format
!
 accessfil = 0
 if (dtset%accesswff == IO_MODE_NETCDF) accessfil = 1
 if (dtset%accesswff == IO_MODE_ETSF) accessfil = 3
 if (dtset%accesswff == IO_MODE_MPI) accessfil = 4

 occopt=dtset%occopt;

 pawfatbnd=dtset%pawfatbnd
 prtden=dtset%prtden ; prtpot=dtset%prtpot ; prtgeo=dtset%prtgeo
 prtcml=dtset%prtcml ; prtdos=dtset%prtdos ; prtstm=dtset%prtstm
 prt1dm=dtset%prt1dm ; prtvha=dtset%prtvha ; prtvhxc=dtset%prtvhxc
 prtvxc=dtset%prtvxc ; prtnabla=dtset%prtnabla; prtefg=dtset%prtefg
 prtcs=dtset%prtcs   ; prtfc=dtset%prtfc ; prtkden=dtset%prtkden
 prtelf=dtset%prtelf ; prtgden=dtset%prtgden; prtlden=dtset%prtlden
 pawprtden=dtset%prtden-1

 call timab(952,1,tsec)

 spacecomm=mpi_enreg%comm_cell
 me=xcomm_rank(spacecomm)

!Warnings :
!- core charge is excluded from the charge density;
!- the potential is the INPUT vtrial.
 if( xmpi_paral==0 .or. &
& (mpi_enreg%paral_kgb==0 .and. me==0 ) .or. &
& (mpi_enreg%paral_kgb==1 .and. mpi_enreg%me_band==0 .and. &
& mpi_enreg%me_kpt==0 .and. mpi_enreg%me_spinor==0)) then

!  We output the density.
   if (prtden/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
!    We create the file name.
     call ioarr(accessfil,rhor,dtset,etotal,fformr,dtfil%fnameabo_app_den, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if (dtset%positron/=0) then
       rdwr=2; fformr=52; rdwrpaw=0
       fname=trim(dtfil%fnameabo_app_den)//'_POSITRON'
       call ioarr(accessfil,electronpositron%rhor_ep, dtset, etotal, fformr,&
&       fname,hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     end if
     if (accessfil==3) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_den, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_den)
     end if
   end if

 end if ! if master

!! MS - Printing of PAWDEN parallellised and several possible options included
!We output the total electron density in the PAW case
!this requires removing nhat from rhor and making PAW on-site corrections
 if (pawprtden>0 .and. psps%usepaw==1) then
!  pawprtden 1 --> output PAW valence density
!  "     2 --> output PAW valence+core density
!  "     3 --> output core, valence and full atomic protodensity
!  "     4 --> options 1+3
!  "     5 --> options 2+3
!  "     6 --> output all individual PAW density contributions
   spacecomm_fft=spacecomm;if (mpi_enreg%paral_kgb==1) spacecomm_fft=mpi_enreg%comm_fft
   if (pawprtden/=3) then ! calc PAW valence density
     ABI_ALLOCATE(rhor_paw,(pawfgr%nfft,nspden))
     ABI_ALLOCATE(rhor_n_one,(pawfgr%nfft,nspden))
     ABI_ALLOCATE(rhor_nt_one,(pawfgr%nfft,nspden))
!    If the communicator used for denfgr is kpt_comm, it is not compatible with paral_atom
     if (mpi_enreg%paral_kgb==0.and.my_natom/=natom) then
       my_natom_tmp=natom
       ABI_DATATYPE_ALLOCATE(pawrhoij_all,(natom))
       call pawrhoij_nullify(pawrhoij_all)
       call pawrhoij_copy(pawrhoij,pawrhoij_all,mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       keep_cplex=.false.,keep_itypat=.false.,keep_nspden=.false.)
     else
       my_natom_tmp=my_natom
       pawrhoij_all => pawrhoij
     end if
     if (pawprtden/=6) then
       call denfgr(atindx1,gmet,spacecomm_fft,my_natom_tmp,natom,nattyp,ngfft,nhat,dtset%nspinor,nsppol,nspden,&
&       ntypat,pawfgr,pawrad,pawrhoij_all,pawtab,prtvol,rhor,rhor_paw,rhor_n_one,&
&       rhor_nt_one,rprimd,dtset%typat,ucvol,xred,&
&       mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     else
       call denfgr(atindx1,gmet,spacecomm_fft,my_natom_tmp,natom,nattyp,ngfft,nhat,dtset%nspinor,nsppol,nspden,&
&       ntypat,pawfgr,pawrad,pawrhoij_all,pawtab,prtvol,rhor,rhor_paw,rhor_n_one,&
&       rhor_nt_one,rprimd,dtset%typat,ucvol,xred,&
&       abs_n_tilde_nt_diff=nt_ntone_norm,znucl=dtset%znucl,&
&       mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if
     if (mpi_enreg%paral_kgb==1.and.my_natom/=natom) then
       ABI_DATATYPE_DEALLOCATE(pawrhoij_all)
     end if

!    Check normalisation
     if (prtvol>9) then
       norm = SUM(rhor_paw(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
       write(message,'(a,F8.4)') '  PAWDEN - NORM OF DENSITY: ',norm
       call wrtout(std_out,message,'COLL')
     end if
   end if

   if (pawprtden>1.AND.pawprtden<6) then ! We will need the core density
     ABI_ALLOCATE(rhor_paw_core,(pawfgr%nfft,nspden))
     call read_atomden(mpi_enreg,natom,nspden,ntypat,pawfgr,rhor_paw_core,&
&     dtset%typat,rprimd,xred,prtvol,file_prefix='core   ')
!    Check normalisation
     if (prtvol>9) then
       norm = SUM(rhor_paw_core(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
       write(message,'(a,F8.4)') '  ATMDEN - NORM OF CORE DENSITY: ', norm
       call wrtout(std_out,message,'COLL')
     end if
   end if

   if (pawprtden>2.AND.pawprtden<6) then ! We will need the valence protodensity
     ABI_ALLOCATE(rhor_paw_val,(pawfgr%nfft,nspden))
     call read_atomden(mpi_enreg,natom,nspden,ntypat,pawfgr,rhor_paw_val,&
&     dtset%typat,rprimd,xred,prtvol,file_prefix='valence')
!    Check normalisation
     if (prtvol>9) then
       norm = SUM(rhor_paw_val(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
       write(message,'(a,F8.4)') '  ATMDEN - NORM OF VALENCE PROTODENSITY: ', norm
       call wrtout(std_out,message,'COLL')
     end if
   end if

   if( xmpi_paral==0 .or. &
&   (mpi_enreg%paral_kgb==0 .and. me==0 ) .or. &
&   (mpi_enreg%paral_kgb==1 .and. mpi_enreg%me_band==0 .and. mpi_enreg%me_kpt==0 .and. mpi_enreg%me_spinor == 0)) then ! if master
     if (pawprtden/=3) then
       if (pawprtden==2.or.pawprtden==5) rhor_paw = rhor_paw + rhor_paw_core
!      PAWDEN
       rdwr=2 ; fformr=52 ; rdwrpaw=0
       call ioarr(accessfil,rhor_paw,dtset,etotal,fformr,dtfil%fnameabo_app_pawden, &
&       hdr,mpi_enreg,pawfgr%nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
       if ( accessfil == 3 ) then
!        Complete the geometry informations with missing values from hdr_io().
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_pawden, psps)
!        Complete the electrons definition with missing values from hdr_io().
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_pawden)
       end if
     end if

     if (pawprtden>2.AND.pawprtden<6) then
!      ATMDEN_CORE
       rdwr=2 ; fformr=52 ; rdwrpaw=0
       call ioarr(accessfil,rhor_paw_core,dtset,etotal,fformr,dtfil%fnameabo_app_atmden_core, &
&       hdr,mpi_enreg,pawfgr%nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
       if ( accessfil == 3 ) then
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_atmden_core, psps)
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_atmden_core)
       end if
!      We create the file name for valence protodensity. ATMDEN_VAL
       call ioarr(accessfil,rhor_paw_val,dtset,etotal,fformr,dtfil%fnameabo_app_atmden_val, &
&       hdr,mpi_enreg,pawfgr%nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
       if ( accessfil == 3 ) then
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_atmden_val, psps)
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_atmden_val)
       end if
!      We create the file name for full protodensity. ATMDEN_FULL
       rhor_paw_val = rhor_paw_val + rhor_paw_core
       call ioarr(accessfil,rhor_paw_val,dtset,etotal,fformr,dtfil%fnameabo_app_atmden_full, &
&       hdr,mpi_enreg,pawfgr%nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
       if ( accessfil == 3 ) then
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_atmden_full, psps)
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_atmden_full)
       end if
     end if

     if (pawprtden==6) then ! Print all individual contributions to the density
!      N_TILDE - N_HAT
!      Use rhor_paw_val as temporary array
       if (.not.allocated(rhor_paw_val))  then
         ABI_ALLOCATE(rhor_paw_val,(pawfgr%nfft,nspden))
       end if
       rhor_paw_val = rhor - nhat
       rdwr=2 ; fformr=52 ; rdwrpaw=0
       call ioarr(accessfil,rhor_paw_val,dtset,etotal,fformr,dtfil%fnameabo_app_n_tilde, &
&       hdr,mpi_enreg,pawfgr%nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
       if ( accessfil == 3 ) then
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_n_tilde, psps)
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_n_tilde)
       end if
!      N_ONSITE
       call ioarr(accessfil,rhor_n_one,dtset,etotal,fformr,dtfil%fnameabo_app_n_one, &
&       hdr,mpi_enreg,pawfgr%nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
       if ( accessfil == 3 ) then
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_n_one, psps)
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_n_one)
       end if
!      N_TILDE_ONSITE
       call ioarr(accessfil,rhor_nt_one,dtset,etotal,fformr,dtfil%fnameabo_app_nt_one, &
&       hdr,mpi_enreg,pawfgr%nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
       if ( accessfil == 3 ) then
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_nt_one, psps)
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_nt_one)
       end if

     end if ! All indivdual density cont.
   end if ! if master

   if (allocated(rhor_paw))  then
     ABI_DEALLOCATE(rhor_paw)
   end if
   if (allocated(rhor_paw_core))  then
     ABI_DEALLOCATE(rhor_paw_core)
   end if
   if (allocated(rhor_paw_val))  then
     ABI_DEALLOCATE(rhor_paw_val)
   end if
   if (allocated(rhor_n_one))  then
     ABI_DEALLOCATE(rhor_n_one)
   end if
   if (allocated(rhor_nt_one))  then
     ABI_DEALLOCATE(rhor_nt_one)
   end if

 end if ! if paw+pawprtden

 call timab(952,2,tsec)
 call timab(953,1,tsec)

 if( xmpi_paral==0                         .or. &
& (mpi_enreg%me==0 .and. mpi_enreg%paral_kgb==0 ) .or. &
& (mpi_enreg%paral_kgb==1 .and. mpi_enreg%me_band==0 .and. &
& mpi_enreg%me_kpt==0 .and. mpi_enreg%me_spinor==0)) then ! if master

!  We output the electron localization function ELF
   if (prtelf/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call ioarr(accessfil,elfr,dtset,etotal,fformr,dtfil%fnameabo_app_elf, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_elf, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_elf)
     end if
     if (nspden==2)then
       ABI_ALLOCATE(elfr_up,(nfft,nspden))
       elfr_up(:,:) = zero
       do ifft=1,nfft
         elfr_up(ifft,1) = elfr(ifft,2)
       end do
!      ELF_UP
       call ioarr(accessfil,elfr_up,dtset,etotal,fformr,dtfil%fnameabo_app_elf_up, &
&       hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
       if ( accessfil == 3 ) then
!        Complete the geometry informations with missing values from hdr_io().
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_elf_up, psps)
!        Complete the electrons definition with missing values from hdr_io().
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_elf_up)
       end if
       ABI_ALLOCATE(elfr_down,(nfft,nspden))
       elfr_down(:,:) = zero
       do ifft=1,nfft
         elfr_down(ifft,1) = elfr(ifft,3)
       end do
!      ELF_DOWN'
       call ioarr(accessfil,elfr_down,dtset,etotal,fformr,dtfil%fnameabo_app_elf_down, &
&       hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
       if ( accessfil == 3 ) then
!        Complete the geometry informations with missing values from hdr_io().
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_elf_down, psps)
!        Complete the electrons definition with missing values from hdr_io().
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_elf_down)
       end if
       ABI_DEALLOCATE(elfr_up)
       ABI_DEALLOCATE(elfr_down)
     end if
   end if

   call timab(953,2,tsec)
   call timab(954,1,tsec)

!  We output the gradient of density
   if (prtgden/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call ioarr(accessfil,grhor(:,:,1),dtset,etotal,fformr,dtfil%fnameabo_app_gden1, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_gden1, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_gden1)
     end if
     call ioarr(accessfil,grhor(:,:,2),dtset,etotal,fformr,dtfil%fnameabo_app_gden2, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_gden2, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_gden2)
     end if
     call ioarr(accessfil,grhor(:,:,3),dtset,etotal,fformr,dtfil%fnameabo_app_gden3, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_gden3, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_gden3)
     end if
   end if

!  We output the total kinetic energy density KDEN
   if (prtkden/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call ioarr(accessfil,taur,dtset,etotal,fformr,dtfil%fnameabo_app_kden, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_kden, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_kden)
     end if
   end if

!  We output the Laplacian of density
   if (prtlden/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call ioarr(accessfil,lrhor,dtset,etotal,fformr,dtfil%fnameabo_app_lden, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_lden, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_lden)
     end if
   end if

   call timab(954,2,tsec)
   call timab(955,1,tsec)

!  We handle the output of wavefunctions. WFK
   if (dtset%prtwf == 1) then
!    In ETSF, some geometric informations are required for wave functions files.
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_wfk, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_wfk)
     end if
   end if

   call timab(955,2,tsec)
   call timab(956,1,tsec)

!  POT
   if (prtpot>0) then
     rdwr=2 ; fformv=102 ; rdwrpaw=0
!    MJV note: why is accessfil forced to 0???? This disables the writing of ETSF
!    format potentials!
!    
!    set to 1 for netcdf output
     accessfil = 0
     call ioarr(accessfil,vtrial,dtset,etotal,fformv,dtfil%fnameabo_app_pot, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_pot, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_pot)
     end if
   end if

   call timab(956,2,tsec)
   call timab(957,1,tsec)

   if (prtgeo>0) then
     coordn=prtgeo
     if ( accessfil == 3 ) then
       call abi_etsf_geo_put(dtset,dtfil%fnameabo_app, psps)
     else
       call bonds_lgth_angles(coordn,dtfil%fnameabo_app_geo,natom,psps%ntypat,&
&       rprimd,dtset%typat,xred,dtset%znucl)
     end if
   end if

   if (any(prtcml==[1,2]))  then
     call prt_cml2(prtcml,dtfil%fnameabo_app_cml_xml,natom,dtset%nsym,psps%ntypat,&
&     rprimd,dtset%spgroup,dtset%symrel,dtset%tnons,dtset%typat,xred,dtset%znucl)
   end if

   if (dtset%prtcif > 0) then
     call prt_cif(dtset%brvltt, dtfil%fnameabo_app_cif, natom, dtset%nsym, dtset%ntypat, rprimd, &
&     dtset%spgaxor, dtset%spgroup, dtset%spgorig, dtset%symrel, dtset%tnons, dtset%typat, xred, dtset%znucl)
   end if

   call timab(957,2,tsec)
   call timab(958,1,tsec)

!  STM
   if (prtstm>0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
!    set to 1 for netcdf output
     call ioarr(accessfil,rhor,dtset,etotal,fformr,dtfil%fnameabo_app_stm, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if ( accessfil == 3 ) then
!      Complete the geometry information with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_stm, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_stm)
     end if
   end if

   if (prt1dm>0) then
     call out1dm(dtfil%fnameabo_app_1dm,mpi_enreg,natom,nfft,ngfft,nspden,psps%ntypat,&
&     rhor,rprimd,dtset%typat,ucvol,vtrial,xred,dtset%znucl)
   end if

!  VHA
   if (prtvha>0) then
     rdwr=2 ; fformv=102 ; rdwrpaw=0
!    set to 1 for netcdf output
     ABI_ALLOCATE(vwork,(nfft,nspden))
     do ispden=1,nspden
       vwork(:,ispden)=vhartr(:)
     end do
     call ioarr(accessfil,vwork,dtset,etotal,fformv,dtfil%fnameabo_app_vha, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     ABI_DEALLOCATE(vwork)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_vha, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_vha)
     end if
   end if

!  VHXC
   if (prtvhxc>0) then
     rdwr=2 ; fformv=102 ; rdwrpaw=0
!    set to 1 for netcdf output
     ABI_ALLOCATE(vwork,(nfft,nspden))
     do ispden=1,nspden
       vwork(:,ispden)=vhartr(:)+vxc(:,ispden)
     end do
     call ioarr(accessfil,vwork,dtset,etotal,fformv,dtfil%fnameabo_app_vhxc, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     ABI_DEALLOCATE(vwork)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_vhxc, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_vhxc)
     end if
   end if

!  VXC
   if (prtvxc>0) then
     rdwr=2 ; fformv=102 ; rdwrpaw=0
!    set to 1 for netcdf output
     call ioarr(accessfil,vxc,dtset,etotal,fformv,dtfil%fnameabo_app_vxc, &
&     hdr,mpi_enreg,nfft,pawrhoij_dum,rdwr,rdwrpaw,wvl_den)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_vxc, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_vxc)
     end if
   end if

   call timab(958,2,tsec)

 end if ! if master

 call timab(959,1,tsec)

!Generate DOS using the tetrahedron method or using Gaussians
!FIXME: Should centralize all calculations of DOS here in outscfcv
 partial_dos_flag = 0
 if (prtdos>=2.or.pawfatbnd>0) then

   if(prtdos==2) partial_dos_flag = 0
   if(prtdos==3) partial_dos_flag = 1
   if(prtdos==4) partial_dos_flag = 1
   prtdosm=0
   if (partial_dos_flag==1) prtdosm=dtset%prtdosm
   paw_dos_flag=0
   if (psps%usepaw==1.and.partial_dos_flag==1.and.dtset%pawprtdos==1) paw_dos_flag=1
!  paw_dos_flag= 1 if both PAW contributions are evaluated AND stored
   fatbands_flag=0
   if(pawfatbnd>0.and.prtdosm==0) fatbands_flag=1
   if(prtdosm==1.and.pawfatbnd>0)then
     message = 'pawfatbnd>0  and prtdosm=1 are not compatible '
!    because they compute quantities in real and complex harmonics respectively
     MSG_ERROR(message)
   end if

!  mjv : initialization is needed as mbesslang is used for allocation below
!  NOTE: 10/5/2010 the whole of this could be looped over ndosfraction,
!  to store much less in memory. The DOS is accumulated in an array
!  and then printed to file at the end.
   mbesslang = 1
   if(partial_dos_flag==1.or.fatbands_flag==1)then
     mbesslang = 5
     ndosfraction=(dtset%natsph + dtset%natsph_extra)*mbesslang
   else
     ndosfraction = 1
     mbesslang = 0
   end if

   ABI_ALLOCATE(dos_fractions,(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction))
   if (prtdosm>=1.or.fatbands_flag==1) then
     ABI_ALLOCATE(dos_fractions_m,(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang))
     ABI_ALLOCATE(dos_fractions_average_m,(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang))
   else
     ABI_ALLOCATE(dos_fractions_m,(0,0,0,0))
     ABI_ALLOCATE(dos_fractions_average_m,(0,0,0,0))
   end if
   if (psps%usepaw==1.and.(partial_dos_flag==1)) then
     ABI_ALLOCATE(dos_fractions_paw1,(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction))
     ABI_ALLOCATE(dos_fractions_pawt1,(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction))
   else
     ABI_ALLOCATE(dos_fractions_paw1,(0,0,0,0))
     ABI_ALLOCATE(dos_fractions_pawt1,(0,0,0,0))
   end if
   if( partial_dos_flag==1.or.fatbands_flag==1)then
!    Generate fractions for partial DOSs if needed
!    partial_dos 1,2,3,4  give different decompositions
     if ((psps%usepaw==0.or.dtset%pawprtdos/=2).and.partial_dos_flag==1) then
       call partial_dos_fractions(cg,dos_fractions,dos_fractions_m,dtfil,dtset,hdr,mbesslang,mcg,mpi_enreg, &
&       prtdosm,ndosfraction,partial_dos_flag)
     else
       dos_fractions=zero;if (prtdosm>=1.or.fatbands_flag==1) dos_fractions_m=zero
     end if
     if (psps%usepaw==1) then
!      TODO: update partial_dos_fractions_paw for extra atoms - no PAW contribution normally, but check bounds and so on.
       call partial_dos_fractions_paw(atindx1,cprj,dimcprj,dos_fractions,dos_fractions_m,&
&       dos_fractions_paw1,dos_fractions_pawt1,dtfil,dtset,fatbands_flag,&
&       mbesslang,mcprj,mkmem,mpi_enreg,prtdosm,ndosfraction,&
&       paw_dos_flag,pawrad,pawtab)
     end if
     if(prtdosm>=1)then
       call dos_degeneratewfs(dos_fractions_m,dos_fractions_average_m,&
&       eigen,mband,dtset%nband,ndosfraction*mbesslang,dtset%nkpt,dtset%nsppol)
     end if
   else
     dos_fractions(:,:,:,1)=one
   end if

!  Here, computation of fatbands for the k-point given. _FATBANDS
   if(pawfatbnd>0.and.fatbands_flag==1) then
     call prtfatbands(dos_fractions_m,dtset,dtfil%fnameabo_app_fatbands,fermie,eigen,&
&     mbesslang,prtdosm,ndosfraction,pawfatbnd,pawtab)
   end if

!  Here, computation and output of DOS and partial DOS  _DOS
   if(fatbands_flag==0) then
     if(prtdos/=4) then
       call tetrahedron (dos_fractions,dos_fractions_average_m,dos_fractions_paw1,dos_fractions_pawt1,&
&       dtset,fermie,eigen,dtfil%fnameabo_app_dos,mbesslang,prtdosm,ndosfraction,paw_dos_flag,rprimd)
     else
!      this option is not documented in input variables: is it working?
       call gaus_dos(dos_fractions,dos_fractions_paw1,dos_fractions_pawt1,dtset,&
&       fermie,eigen,dtfil%fnameabo_app_dos,mbesslang,prtdosm,ndosfraction,paw_dos_flag)
     end if
   end if

   ABI_DEALLOCATE(dos_fractions)
   ABI_DEALLOCATE(dos_fractions_m)
   !if (prtdosm==1.or.fatbands_flag==1)  then
   ABI_DEALLOCATE(dos_fractions_average_m)
   !end if
   !if (psps%usepaw==1.and.(partial_dos_flag==1)) then
   ABI_DEALLOCATE(dos_fractions_paw1)
   ABI_DEALLOCATE(dos_fractions_pawt1)
   !end if

 end if ! prtdos > 1
 call timab(959,2,tsec)
 call timab(960,1,tsec)

!Output of integrated density inside atomic spheres
 if (dtset%prtdensph==1)then
   call calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,&
&   ntypat,ab_out,dtset%ratsph,rhor,rprimd,dtset%typat,ucvol,xred)
 end if

 call timab(960,2,tsec)
 
 if (dtset%magconon /= 0) then
!  calculate final value of terms for magnetic constraint: "energy" term, lagrange multiplier term, and atomic contributions
   call mag_constr_e(dtset%magconon,dtset%magcon_lambda,mpi_enreg,&
&   natom,nfft,ngfft,nspden,ntypat,dtset%ratsph,rhor,rprimd,dtset%spinat,dtset%typat,xred)
 end if

 call timab(961,1,tsec)

!If PAW, provide additional outputs
 if (psps%usepaw==1) then
!  Output of compensation charge
   if (dtset%nstep>0.or.dtfil%ireadwf/=0) then
     write(message, '(4a)' )ch10,' PAW TEST:',ch10,&
&     ' ==== Compensation charge inside spheres ============'
     if (compch_sph>-1.d4.and.compch_fft>-1.d4) &
&     write(message, '(3a)' ) trim(message),ch10,' The following values must be close to each other ...'
     if (compch_sph>-1.d4) write(message, '(3a,f22.15)' ) trim(message),ch10,&
&     ' Compensation charge over spherical meshes = ',compch_sph
     if (compch_fft>-1.d4) then
       if (pawfgr%usefinegrid==1) then
         write(message, '(3a,f22.15)' ) trim(message),ch10,&
&         ' Compensation charge over fine fft grid    = ',compch_fft
       else
         write(message, '(3a,f22.15)' ) trim(message),ch10,&
&         ' Compensation charge over fft grid         = ',compch_fft
       end if
     end if
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
!  Output of pseudopotential strength Dij and augmentation occupancies Rhoij
   call pawprt(dtset,my_natom,paw_ij,pawrhoij,pawtab,&
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&   electronpositron=electronpositron)
 end if

 call timab(961,2,tsec)
 call timab(962,1,tsec)


!PAW + output for optical conductivity   _OPT and _OPT2
 if (psps%usepaw==1.and.prtnabla>0) then
   call optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,eigen,gprimd,hdr,kg,&
&   mband,mcg,mcprj,mkmem,mpi_enreg,mpsang,mpw,natom,nkpt,npwarr,nsppol,pawrad,pawtab)
   if (prtnabla>1) then
     call optics_paw_core(atindx1,cprj,dimcprj,dtfil,dtset,eigen,hdr,&
&     mband,mcprj,mkmem,mpi_enreg,mpsang,natom,nkpt,nsppol,pawrad,pawtab)
   end if
 end if
 if (prtnabla<0) then
   call optics_vloc(cg,dtfil,dtset,eigen,gprimd,hdr,kg,&
&   mband,mcg,mkmem,mpi_enreg,mpw,nkpt,npwarr,nsppol)
 end if

 call timab(962,2,tsec)
 call timab(963,1,tsec)



!Optionally provide output for AE wavefunctions (only for PAW)
 if (psps%usepaw==1 .and. dtset%pawprtwf==1) then
   ABI_ALLOCATE(ps_norms,(nsppol,nkpt,mband))

   call pawmkaewf(Dtset,my_natom,natom,mpw,mband,mcg,mcprj,nkpt,mkmem,nsppol,ntypat,Dtset%nband,&
&   Dtset%istwfk,npwarr,Dtset%kptns,Dtset%paral_kgb,Dtset%ngfftdg,kg,dimcprj,pawfgrtab,&
&   Pawrad,Pawtab,gmet,rprimd,ucvol,Psps,Hdr,Dtfil,eigen,occ,cg,Cprj,&
&   MPI_enreg,ierr,pseudo_norms=ps_norms,set_k=dtset%pawprt_k,set_band=dtset%pawprt_b,&
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   if (dtset%pawprt_b==0) then
     if (open_file(Dtfil%fnameabo_pawstat,message,newunit=tmp_unt,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(message)
     end if
     write(tmp_unt,'(5a)') '# This file contains the statistics on the cancellation of',ch10,&
&     '# the onsite pseudo component of the all-electron wavefunction',ch10,&
&     '# with the plane wave part'
     do isppol=1,nsppol
       write(tmp_unt,'(a,i0)') '# isppol = ',isppol
       do ikpt=1,nkpt
         write(tmp_unt,'(a,i0)') '# ikpt = ',ikpt
         write(tmp_unt,'(a)') '#    band      norm'
         occ_norm = zero; unocc_norm = zero; nocc = 0
         do iband=1,mband
           write(tmp_unt,'(i8,ES16.6)') iband,ps_norms(isppol,ikpt,iband)
! WARNING : this usage of occ looks wrong to me... mjv. The dimensionality is ok, but the offset is bizarre.
           if (occ(iband*ikpt*isppol)==zero) then
             unocc_norm = unocc_norm + ps_norms(isppol,ikpt,iband)
           else
             occ_norm = occ_norm + ps_norms(isppol,ikpt,iband)
             nocc = nocc + 1
           end if
         end do
         if(mband/=nocc)then
           write(tmp_unt,'(2(a,ES16.6))') '# occ average: ',occ_norm/real(nocc),&
&           ' unocc average: ',unocc_norm/real(mband-nocc)
         else
           write(tmp_unt,'(2(a,ES16.6))') '# occ average: ',occ_norm/real(nocc)
         end if
       end do
     end do
     close(tmp_unt)
   end if
   ABI_DEALLOCATE(ps_norms)
 end if

 call timab(963,2,tsec)
 if(dtset%plowan_compute>0) then
   write(message,'(2a,i3)') ch10,&
&   ' ====================================================================================== '
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   write(message,'(2a,i3)') ch10,&
&   ' == Start computation of Projected Local Orbitals Wannier functions == ',dtset%nbandkss
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   
!  ==  compute psichi

   call init_plowannier(dtset,wan)
   call compute_coeff_plowannier(Crystal,cprj,dimcprj,dtset,eigen,e_fermie,&
&   mpi_enreg,occ,wan,pawtab,psps,&
&   usecprj,dtfil%unpaw,pawrad,dtfil)
   call destroy_plowannier(wan)
   
 end if

!Optionally provide output for the GW part of ABINIT
 if (dtset%nbandkss/=0) then
   ! Use DMFT to compute wannier function for cRPA calculation.
   if(dtset%usedmft==1) then
     write(message,'(2a,i3)') ch10,&
&     '  Warning: Psichi are renormalized in datafordmft because nbandkss is used',dtset%nbandkss
     call wrtout(std_out,message,'COLL')
     call init_dmft(dtset,e_fermie,dtfil%fnameabo_app,dtset%nspinor,paw_dmft,pawtab,psps,dtset%typat)
     call print_dmft(paw_dmft,dtset%pawprtvol)
     
!    ==  compute psichi
     call init_oper(paw_dmft,lda_occup)

     call datafordmft(Crystal,cprj,dimcprj,dtset,eigen,e_fermie,&
&     lda_occup,dtset%mband,dtset%mband,dtset%mkmem,mpi_enreg,&
&     dtset%nkpt,dtset%nspinor,dtset%nsppol,occ,&
&     paw_dmft,paw_ij,pawang,pawtab,psps,usecprj,dtfil%unpaw,dtset%nbandkss) 

     call destroy_dmft(paw_dmft)
     call destroy_oper(lda_occup)
   end if

   call timab(964,1,tsec) ! outscfcv(outkss)

   call outkss(dtfil,dtset,ecut,gmet,gprimd,hdr,&
&   dtset%kssform,mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,natom,&
&   nfft,nkpt,npwarr,nspden,nsppol,nsym,psps%ntypat,occ,pawtab,pawfgr,paw_ij,&
&   prtvol,psps,rprimd,vtrial,xred,cg,usecprj,cprj,eigen,ierr)
   call timab(964,2,tsec) ! outscfcv(outkss)
   if (ierr/=0) then
     MSG_WARNING("outkss returned a non zero status error, check log")
   end if
 end if
!
!Output of the GSR file.
 if (me==master) then
   fname = strcat(dtfil%filnam_ds(4),"_GSR.nc")
   if (isalchemical(Crystal)) then
     MSG_WARNING("Alchemical pseudos are not supported by ETSF-IO, GSR file won't be produced")
   else
#ifdef HAVE_TRIO_ETSF_IO
     NCF_CHECK(ncfile_create(ncf,fname,NF90_CLOBBER),"Creating GSR file")
     call crystal_ncwrite(Crystal,ncf%ncid)
     call ebands_ncwrite(Bands,dtset%nshiftk_orig,dtset%shiftk_orig,dtset%nshiftk,dtset%shiftk,dtset%ngkpt,dtset%kptrlatt,ncf%ncid)
     call results_gs_ncwrite(results_gs,dtset%ecut,ncf%ncid)
     NCF_CHECK(ncfile_close(ncf),"Closing GSR file")
#else
     ABI_UNUSED(ncf%ncid)
#endif
   end if
 end if


 if (electronpositron_calctype(electronpositron)/=0) then

!  Optionally provide output for  positron life time calculation
   call timab(965,1,tsec)
   gammastate_core=zero
   call poslifetime(dtset,electronpositron,gammastate_core,gprimd,my_natom,&
&   mpi_enreg,n3xccc,nfft,ngfft,nhat,1,pawang,&
&   pawrad,pawrhoij,pawtab,rate_dum,rate_dum,&
&   rhor,ucvol,xccc3d)
   call timab(965,2,tsec)

!  Optionally provide output for momentum distribution of annihilation radiation
   if (dtset%posdoppler>0) then
     call posdoppler(cg,cprj,Crystal,dimcprj,dtfil,dtset,electronpositron,psps%filpsp,&
&     gammastate_core,kg,mcg,mcprj,mpi_enreg,my_natom,n3xccc,nfft,ngfft,nhat,npwarr,&
&     occ,pawang,pawrad,pawrhoij,pawtab,rhor,xccc3d)
   end if
 end if

!Optionally provide output for WanT
 if (dtset%prtwant==1) then
   call timab(966,1,tsec)
   ! WARNING: mpi_enreg not used --> MPI is not supported
   call outwant(dtset,eigen,cg,kg,npwarr,mband,mcg,mpi_enreg,nkpt,nsppol,mkmem,mpw,dtset%prtwant)
   call timab(966,2,tsec)
 end if

!Optionally provide output for chemical shielding calculation
 if (prtcs > 0) then
   call timab(967,1,tsec)
   call calc_cs(cg,dtefield,gprimd,kg,dtset%kptns,mcg,mkmem,mpi_enreg,mpw,natom,nkpt,npwarr,occopt,psps%usepaw)
   call timab(967,2,tsec)
 end if

!Optionally provide output for electric field gradient calculation
 if (prtefg > 0) then
   call timab(967,1,tsec)
   call calc_efg(gprimd,my_natom,natom,nfft,ngfft,nspden,dtset%nsym,ntypat,dtset%paral_kgb,&
&   paw_an,pawang,pawrad,pawrhoij,pawtab,&
&   dtset%ptcharge,prtefg,dtset%quadmom,rhor,rprimd,dtset%symrel,&
&   dtset%tnons,dtset%typat,ucvol,psps%usepaw,xred,psps%zionpsp,&
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call timab(967,2,tsec)
 end if

!Optionally provide output for Fermi-contact term at nuclear positions
 if (prtfc > 0) then
   call timab(967,1,tsec)
   call calc_fc(my_natom,natom,nspden,ntypat,pawrad,pawrhoij,pawtab,dtset%typat,psps%usepaw,&
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call timab(967,2,tsec)
 end if

!Optionally provide Xcrysden output for the Fermi surface (Only master writes)
 if (dtset%prtfsurf==1.and.me==0) then
   if (ebands_3dprint(Bands,Crystal,dtset%kptrlatt,dtset%nshiftk,dtset%shiftk,dtfil%fnameabo_app_bxsf) /= 0) then
     MSG_WARNING("Cannot produce file for Fermi surface, check log file for more info")
   end if
 end if ! prtfsurf==1

!output nesting factor for Fermi surface
 if (dtset%prtnest>0 .and. me==0) then
   call timab(968,1,tsec)
!  
!  FIXME: these should become input variables for abinit as well
   nqpath = 5
   ABI_ALLOCATE(qpath_vertices,(3,nqpath))
   qpath_vertices(:,1) = (/zero, zero, zero/)
   qpath_vertices(:,2) = (/half, half, zero/)
   qpath_vertices(:,3) = (/half, zero, zero/)
   qpath_vertices(:,4) = (/zero, zero, zero/)
   qpath_vertices(:,5) = (/half, half, half/)

   skipnest = 0
   do ikpt = 1, dtset%nkpt
     if (dtset%nband(ikpt) /= dtset%nband(1)) then
       message = 'Error: mknesting can not handle variable nband(1:nkpt). Skipped.'//&
&       ch10//' Correct input file to get nesting output'
       call wrtout(ab_out,message,'COLL')
       skipnest = 1
     end if
   end do
   if (skipnest == 0) then
!    FIXME: needs to be generalized to complete the k grid for one of the arguments to mknesting

!    generate weights for FS
     ABI_ALLOCATE(fs_weights,(dtset%nband(1),dtset%nkpt,dtset%nsppol))
     invgauwidth = one / (0.1_dp * Ha_eV) ! default is sigma = 0.1 eV
     if (dtset%tsmear > tol10) invgauwidth = one / dtset%tsmear
     prefact = one / sqrt(pi) * invgauwidth
     ii=1
     do isppol = 1, dtset%nsppol
       do ikpt = 1, dtset%nkpt
         do iband = 1, dtset%nband(1)
           fs_weights(iband, ikpt, isppol) = prefact*exp(-(invgauwidth*(eigen(ii)-(fermie+dtset%fermie_nest)))**2)
           ii = ii+1
         end do
       end do
     end do
     unitmatrix = (/1,0,0,0,1,0,0,0,1/)
     if (dtset%kptopt == 3) then
       call mknesting(dtset%nkpt,hdr%kptns,dtset%kptrlatt,dtset%nband(1),fs_weights,nqpath,&
&       qpath_vertices,1,(/zero, zero, zero/),dtfil%fnameabo_app_nesting,gprimd,gmet,dtset%prtnest, unitmatrix)
     else
       call mknesting(dtset%nkpt,hdr%kptns,dtset%kptrlatt,dtset%nband(1),fs_weights,nqpath,&
&       qpath_vertices,1, (/zero, zero, zero/), dtfil%fnameabo_app_nesting,gprimd,gmet,dtset%prtnest, unitmatrix, &
&       nsym, Crystal%symrec)
     end if
     ABI_DEALLOCATE(fs_weights)
   end if
   ABI_DEALLOCATE(qpath_vertices)
   call timab(968,2,tsec)
 end if ! prtnest=1

 call timab(969,1,tsec)
!
 if (dtset%prtdipole == 1) then
!  FIXME: need to add ionic part of multipoles
   call multipoles_out(rhor,mpi_enreg,natom,nfft,ngfft,dtset%nspden,dtset%ntypat,rprimd,&
&   dtset%typat,ucvol,xred,dtset%ziontypat)
 end if ! prtmultipoles

!BoltzTraP output files in GENEric format
 if (dtset%prtbltztrp == 1) then

   nelect_per_spin = get_nelect_per_spin(Bands) 

   call prtbltztrp_out (eigen, fermie, dtfil%filnam_ds(4), hdr%kptns, natom, dtset%nband(1), &
&   nelect_per_spin, dtset%nkpt, dtset%nspinor, dtset%nsppol, nsym, rprimd, dtset%symrel)
 end if !prtbltztrp

!Free memory.
 call crystal_free(Crystal)
 call ebands_free(Bands)

 call timab(969,2,tsec)
 call timab(950,2,tsec) ! outscfcv

 DBG_EXIT("COLL")

end subroutine outscfcv
!!***
