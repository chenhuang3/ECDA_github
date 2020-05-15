!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtorho3
!! NAME
!! vtorho3
!!
!! FUNCTION
!! This routine compute the new 1-density from a fixed 1-potential (vtrial1)
!! but might also simply compute eigenvectors and eigenvalues.
!! The main part of it is a wf update over all k points
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, AR, DRH, MB, XW, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)=pw coefficients of GS wavefunctions at k+q.
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL; if 2, COMPLEX
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!              projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  dbl_nnsclo=if 1, will double the value of dtset%nnsclo
!!  dimcprj(natom*usepaw)=array of dimensions of arrays cprj, cprjq (ordered by atom-type)
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  dtefield = variables related to response Berry-phase calculation
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  fermie1=derivative of fermi energy wrt (strain) perturbation
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir=direction of the perturbation
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points.
!!  mband=maximum number of bands
!!  mkmem =number of k points treated by this node (GS data).
!!  mkqmem =number of k+q points treated by this node (GS data)
!!  mk1mem =number of k points treated by this node (RF data)
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs - see comment in respfn.F90)
!!  nkpt_rbz=number of k points in the IBZ for this perturbation
!!  mpi_enreg=informations about MPI parallelization
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with
!!   perturbation
!!  ntypat=number of types of atoms in unit cell.
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k+q point of the reduced Brillouin zone.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band and k
!!   (usually 2)
!!  optres=0: the new value of the density is computed in place of the input value
!!         1: only the density residual is computed ; the input density is kept
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  paw_ij1(natom) <type(paw_ij_type)>=1st-order paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawang1 <type(pawang_type)>=pawang datastructure containing only the symmetries preserving the perturbation
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  ph1d(2,3*(2*dtset%mgfft+1)*natom)=one-dimensional structure factor information
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!!  qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3) =
!!  inverse of the overlap matrix
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional real space primitive translations
!!  symaf1(nsym1)=(anti)ferromagnetic part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  symrl1(3,3,nsym1)=symmetry operations in real space
!!  ucvol=unit cell volume in bohr**3.
!!  usecprj= 1 if cprj, cprjq, cprj1 arrays are stored in memory
!!  useylmgr1= 1 if ylmgr1 array is allocated
!!  wffddk=struct info for wf dot (ddk) file.
!!  vtrial(nfftf,nspden)=GS Vtrial(r).
!!  vtrial1(cplex*nfftf,nspden)=INPUT RF Vtrial(r).
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+g point
!!  ylmgr1(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics for each G and k+g point
!!
!!
!! OUTPUT
!!  cg1(2,mpw*nspinor*mband*mk1mem*nsppol)=updated wavefunctions, orthogonalized to the occupied states
!!  cg1_active(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q. They are orthogonalized to the active.
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues
!!    (hartree)
!!  edocc=correction to 2nd-order total energy coming from changes of occupation
!!  eeig0=0th-order eigenenergies part of 2nd-order total energy
!!  ek0=0th-order kinetic energy part of 2nd-order total energy.
!!  ek1=1st-order kinetic energy part of 2nd-order total energy
!!    (not for phonons)
!!  eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!!  enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!!  enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!!  gh1c_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(1)}|nK>
!!  gh0c1_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(0)}|\Psi^{(1)}>
!!      The wavefunction is orthogonal to the active space (for metals). It is not
!!      coherent with cg1.
!!  resid(mband*nkpt_rbz*nsppol)=residuals for each band over all k points.
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!  rhog1(2,nfftf)=RF electron density in reciprocal space
!!  ==== if optres==1
!!    nres2=square of the norm of the residual
!!    nvresid1(cplex*nfftf,nspden)=1st-order density residual
!!  ==== if psps%usepaw==1
!!    cprj1(natom,nspinor*mband*mk1mem*nsppol*usecprj)=
!!              1st-order wave functions at k,q projected with non-local projectors:
!!                       cprj1=<p_i|C1nk,q> where p_i is a non-local projector
!!    nhat1(cplex*nfftf,nspden*psps%usepaw)=1st-order compensation charge density
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!  pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!  rhor1(cplex*nfftf,nspden)=RF electron density in electrons/bohr**3.
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      clsopn,destroy_hamiltonian,destroy_rf_hamiltonian,fftpac
!!      finalize_hamiltonian,gbefd3,gradberry3,hdr_skip,init_hamiltonian
!!      init_rf_hamiltonian,kpg3,kpgstr,load_paw_hamiltonian
!!      load_paw_rf_hamiltonian,mkffnl,mkkin,mkkpg,occeig,pawmkrho
!!      pawrhoij_alloc,pawrhoij_destroy,pawrhoij_init_unpacked
!!      pawrhoij_mpisum_unpacked,ph1d3d,sphereboundary,sqnorm_v,symrhg,timab
!!      transgrid,vtowfk3,wffreadskipk,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine vtorho3(atindx,atindx1,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,dbl_nnsclo,&
& dimcprj,dim_eig2rf,doccde_rbz,docckqde,dtefield,dtfil,dtset,&
& edocc,eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,enl0,enl1,&
& fermie1,gh0c1_set,gh1c_set,gmet,gprimd,idir,indsy1,&
& ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,mband,&
& mkmem,mkqmem,mk1mem,mpi_enreg,mpw,mpw1,my_natom,&
& natom,nband_rbz,ncpgr,nfftf,nhat1,nkpt_rbz,npwarr,npwar1,nres2,nspden,&
& nsppol,nsym1,ntypat,nvresid1,occkq,occ_rbz,optres,&
& paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrhoij,pawrhoij1,pawtab,&
& phnons1,ph1d,prtvol,psps,pwindall,qmat,resid,residm,rhog1,rhor1,rmet,rprimd,symaf1,symrc1,symrl1,ucvol,&
& usecprj,useylmgr1,wffddk,vtrial,vtrial1,wtk_rbz,xred,ylm,ylm1,ylmgr1)


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_wffile
 use m_efield
 use m_hamiltonian
 use m_wfk
 use m_cgtools

 use m_header,   only : hdr_skip, hdr_io
 use m_pawang,   only : pawang_type
 use m_pawtab,   only : pawtab_type
 use m_paw_ij,   only : paw_ij_type
 use m_pawfgrtab,only : pawfgrtab_type
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_alloc, pawrhoij_destroy, &
&                       pawrhoij_init_unpacked, pawrhoij_destroy_unpacked, &
&                       pawrhoij_mpisum_unpacked
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_destroy, pawcprj_get
 use m_pawfgr,   only : pawfgr_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vtorho3'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_72_response, except_this_one => vtorho3
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,dbl_nnsclo,dim_eig2rf,idir,ipert,mband,mk1mem,mkmem
 integer,intent(in) :: mkqmem,mpw,mpw1,my_natom,natom,ncpgr,nfftf,nkpt_rbz,nspden
 integer,intent(in) :: nsppol,nsym1,ntypat,optres,prtvol,usecprj,useylmgr1
 real(dp),intent(in) :: fermie1,ucvol
 real(dp),intent(out) :: edocc,eeig0,ek0,ek1,eloc0,enl0,enl1,nres2,residm
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(efield_type),intent(in) :: dtefield
 type(pawang_type),intent(in) :: pawang,pawang1
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffddk
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),dimcprj(natom*psps%usepaw)
 integer,intent(in) :: indsy1(4,nsym1,natom)
!                      nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),npwar1(nkpt_rbz,2)
 integer,intent(in) :: npwarr(nkpt_rbz,2),symaf1(nsym1),symrc1(3,3,nsym1),symrl1(3,3,nsym1)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)
 real(dp),intent(inout) :: cg1(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol)
 real(dp),intent(out):: cg1_active(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh1c_set(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh0c1_set(2,mpw1*dtset%nspinor*mband*mk1mem*nsppol*dim_eig2rf)
 real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*mband*mkqmem*nsppol)
 real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
 real(dp),intent(out) :: eigen1(2*mband*mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigenq(mband*nkpt_rbz*nsppol),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz),occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: occkq(mband*nkpt_rbz*nsppol),ph1d(2,3*(2*dtset%mgfft+1)*natom)
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
 real(dp), intent(out) :: nhat1(cplex*nfftf,dtset%nspden*psps%usepaw)
 real(dp),intent(out) :: resid(mband*nkpt_rbz*nsppol),rhog1(2,nfftf)
 real(dp),intent(inout) :: nvresid1(cplex*nfftf,nspden),rhor1(cplex*nfftf,nspden)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3)
 real(dp),intent(in),target :: vtrial(nfftf,nspden)
 real(dp),intent(inout) :: vtrial1(cplex*nfftf,nspden)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt_rbz,2,3)
 type(pawcprj_type),intent(in) :: cprj (natom,dtset%nspinor*mband*mkmem *nsppol*usecprj)
 type(pawcprj_type),intent(in) :: cprjq(natom,dtset%nspinor*mband*mkqmem*nsppol*usecprj)
 type(pawcprj_type),intent(inout) :: cprj1(natom,dtset%nspinor*mband*mk1mem*nsppol*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw),paw_ij1(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*psps%usepaw)
 type(pawrhoij_type),target,intent(inout) :: pawrhoij1(my_natom*psps%usepaw)
 type(pawtab_type), intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=13
 integer :: bd2tot_index,bdtot_index,buffer_size,counter,cplex_rhoij,dimffnl1
 integer :: dimffnlk,dimphkxred,ia,iatom,iband,nlines_done
 integer :: ibdkpt,ibg,ibg1,ibgq,icg,icg1,icgq,ider,idir0,ierr
 integer :: ii,ikg,ikg1,ikpt,ilm,ilmn,index1
 integer :: ipw,iscf_mod,isppol,istr,istwf_k
 integer :: itypat,mbd2kpsp,mbdkpsp,mcgq,mcgq_disk,mcprjq
 integer :: mcprjq_disk,me,n1,n2,n3,n4,n5,n6,nband_k,nband_kq,nkpg,nkpg1
 integer :: nnsclo_now,npw1_k,npw_k,nspden_rhoij !nspinor_,
 integer :: spaceworld,test_dot,tim_rwwf
 logical :: paral_atom
 real(dp) :: arg,wtk_k
 !character(len=500) :: msg
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
!arrays
 integer,allocatable :: gbound(:,:),kg1_k(:,:)
 integer,allocatable :: kg_k(:,:)
 integer, pointer :: my_atmtab(:)
 real(dp) :: kpoint(3),kpq(3),rhodum(1)
 real(dp) :: tsec(2),ylmgr_dum(1,1,1)
 real(dp),allocatable :: buffer1(:),cgrvtrial(:,:),dkinpw(:)
 real(dp),allocatable :: doccde_k(:),doccde_kq(:)
 real(dp),allocatable :: edocc_k(:),eeig0_k(:),eig0_k(:),eig0_kq(:),eig1_k(:)
 real(dp),allocatable :: ek0_k(:),ek1_k(:)
 real(dp),allocatable :: eloc0_k(:),enl0_k(:),enl1_k(:),ffnl1(:,:,:,:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnlkq(:,:,:,:),ffspl_typ(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:,:),kinpw1(:),kpg1_k(:,:)
 real(dp),allocatable :: kpg_k(:,:),occ_k(:),occ_kq(:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),resid_k(:)
 real(dp),allocatable :: rho1wfg(:,:),rho1wfr(:,:),rhoaug1(:,:,:),rocceig(:,:)
 real(dp),allocatable :: vlocal(:,:,:),vlocal1(:,:,:)
 real(dp),allocatable :: ylm1_k(:,:),ylm_k(:,:),ylmgr1_k(:,:,:)
 real(dp),pointer :: vtrial_ptr(:,:)
 type(pawrhoij_type),pointer :: pawrhoij1_unsym(:) 

! *********************************************************************

 DBG_ENTER('COLL')

 ABI_UNUSED(atindx1(1))

!Keep track of total time spent in this routine
 call timab(121,1,tsec)
 call timab(124,1,tsec)

!Retrieve parallelism data
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
 paral_atom=(my_natom/=natom)
 my_atmtab=>mpi_enreg%my_atmtab
 
 if (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
& dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17) then
   ABI_ALLOCATE(grad_berry,(2,mpw1,dtefield%nband_occ))
 else
   ABI_ALLOCATE(grad_berry,(0,0,0))
 end if

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
   MSG_BUG('wrong values for nfft, nfftf!')
 end if

!The value of iscf must be modified if ddk perturbation, see loper3.f
 iscf_mod=dtset%iscf;if(ipert==natom+1 ) iscf_mod=-3

 edocc=zero ; eeig0=zero ; ek0=zero  ; ek1=zero
 eloc0=zero ; enl0=zero ; enl1=zero
 bdtot_index=0
 bd2tot_index=0
 ibg=0;icg=0
 ibgq=0;icgq=0
 ibg1=0;icg1=0
 mbdkpsp=mband*nkpt_rbz*nsppol
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol

 dimphkxred=0
 if (ipert>=1.and.ipert<=natom) then
   if (psps%usepaw==0) dimphkxred=1
   if (psps%usepaw==1) dimphkxred=natom
 end if

 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(kg1_k,(3,mpw1))

!Initialize PW 1st-order density if needed
!Also store old rho1 in case of density mixing
 if (iscf_mod>0) then
   if (optres==1) nvresid1=rhor1
   if (psps%usepaw==0) then
     rhor1(:,:)=zero
   else
     ABI_ALLOCATE(rho1wfr,(cplex*dtset%nfft,dtset%nspden))
     ABI_ALLOCATE(rho1wfg,(2,dtset%nfft))
     rho1wfr(:,:)=zero
   end if
 end if

!Set max number of non-self-consistent loops nnsclo_now for use in vtowfk3
 if(iscf_mod<=0 .and. iscf_mod/=-3)then
   nnsclo_now=dtset%nstep
 else
   if(dtset%nnsclo>0)then
     nnsclo_now=dtset%nnsclo
   else
     nnsclo_now=1
   end if
   if(dbl_nnsclo==1) nnsclo_now=nnsclo_now*2
 end if

 n1=dtset%ngfft(1); n2=dtset%ngfft(2); n3=dtset%ngfft(3)
 n4=dtset%ngfft(4); n5=dtset%ngfft(5); n6=dtset%ngfft(6)
 ABI_ALLOCATE(rhoaug1,(cplex*n4,n5,n6))
 ABI_ALLOCATE(vlocal,(n4,n5,n6))
 ABI_ALLOCATE(vlocal1,(cplex*n4,n5,n6))

!Prepare GS k+q wf 
 mcgq=mpw1*dtset%nspinor*mband*mkqmem*nsppol;mcgq_disk=0

!Prepare RF PAW files 
 if (psps%usepaw==1) then
   mcprjq=dtset%nspinor*mband*mkqmem*nsppol*usecprj;mcprjq_disk=0
 else
   mcprjq=0;mcprjq_disk=0
 end if

!Initialisation of the wfdot file in case of electric field
 test_dot=0
 if ((ipert==natom+2.or.ipert==natom+5).and.sum((dtset%qptn(1:3))**2 )<=tol7.and.&
& (dtset%berryopt/= 4.and.dtset%berryopt/= 6.and.dtset%berryopt/= 7.and.&
& dtset%berryopt/=14.and.dtset%berryopt/=16.and.dtset%berryopt/=17)) then
   test_dot=1
   call clsopn(wffddk)
   call hdr_skip(wffddk,ierr)
 end if

!==== Initialize most of the Hamiltonian (and derivative) ====
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
!* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,nspden,natom,&
& ntypat,dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
& ph1d=ph1d,use_gpu_cuda=dtset%use_gpu_cuda)

 call init_rf_hamiltonian(cplex,gs_hamkq,ipert,dtset%nspinor,rf_hamkq,dtset%typat,has_e1kbsc=1)

!PAW:allocate memory for non-symetrized 1st-order occupancies matrix (pawrhoij1)
 pawrhoij1_unsym => pawrhoij1
 if (psps%usepaw==1.and.iscf_mod>0) then
   if (paral_atom) then
     ABI_DATATYPE_ALLOCATE(pawrhoij1_unsym,(natom))
     cplex_rhoij=max(cplex,dtset%pawcpxocc);nspden_rhoij=dtset%nspden
     call pawrhoij_alloc(pawrhoij1_unsym,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&     dtset%nsppol,dtset%typat,pawtab=pawtab,use_rhoijp=0,use_rhoij_=1)
   else
     pawrhoij1_unsym => pawrhoij1
     call pawrhoij_init_unpacked(pawrhoij1_unsym)
   end if
 end if

 nlines_done = 0

!LOOP OVER SPINS
 do isppol=1,nsppol

!  Rewind kpgsph data file if needed:
   ikg=0;ikg1=0

!  Set up local potential vlocal1 with proper dimensioning, from vtrial1
!  Same thing for vlocal from vtrial Also take into account the spin.
   vtrial_ptr=>vtrial
   if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
     call fftpac(isppol,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial_ptr,vlocal,2)
     call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,vtrial1,vlocal1,2)
   else
     ABI_ALLOCATE(cgrvtrial,(dtset%nfft,nspden))
     call transgrid(1,mpi_enreg,nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial_ptr)
     call fftpac(isppol,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
     ABI_DEALLOCATE(cgrvtrial)
     ABI_ALLOCATE(cgrvtrial,(cplex*dtset%nfft,nspden))
     call transgrid(cplex,mpi_enreg,nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial1)
     call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal1,2)
     ABI_DEALLOCATE(cgrvtrial)
   end if
   nullify(vtrial_ptr)

!  Continue to initialize the Hamiltonian
   call load_paw_hamiltonian(gs_hamkq,isppol,paw_ij, &
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call load_paw_rf_hamiltonian(gs_hamkq,ipert,isppol,paw_ij1,rf_hamkq,dtset%typat, &
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!  Nullify contribution to 1st-order density from this k-point
   rhoaug1(:,:,:)=zero

   call timab(125,1,tsec)

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz
     counter=100*ikpt+isppol

     nband_k = nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k = istwfk_rbz(ikpt)
     npw_k   = npwarr(ikpt,1)
     npw1_k  = npwar1(ikpt,1)
     wtk_k   = wtk_rbz(ikpt)

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       eigen1(1+bd2tot_index : 2*nband_k**2+bd2tot_index) = zero
       resid(1+bdtot_index : nband_k+bdtot_index) = zero
       bdtot_index=bdtot_index+nband_k
       bd2tot_index=bd2tot_index+2*nband_k**2

       if (test_dot==1) then
!        Skip the wavefunction block for ikpt,isppol
         tim_rwwf=0
         call WffReadSkipK(1,0,ikpt,isppol,mpi_enreg,wffddk)
!        End the treatment of the dot file
       end if
       cycle ! Skip the rest of the k-point loop
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylm1_k,(npw1_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

!    Finalize to initialize the Hamiltonian at k+q
     kpoint(:)=kpt_rbz(:,ikpt)
     kpq(:)=kpoint(:);if (ipert<natom+3) kpq(:)=kpq(:)+dtset%qptn(1:3)
     call finalize_hamiltonian(gs_hamkq,npw1_k,istwf_k,kpq)

     ABI_ALLOCATE(doccde_k,(nband_k))
     ABI_ALLOCATE(doccde_kq,(nband_k))
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(eig0_kq,(nband_k))
     ABI_ALLOCATE(eig1_k,(2*nband_k**2))
     ABI_ALLOCATE(edocc_k,(nband_k))
     ABI_ALLOCATE(eeig0_k,(nband_k))
     ABI_ALLOCATE(ek0_k,(nband_k))
     ABI_ALLOCATE(ek1_k,(nband_k))
     ABI_ALLOCATE(eloc0_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(occ_kq,(nband_k))
     ABI_ALLOCATE(resid_k,(nband_k))
     ABI_ALLOCATE(gbound,(2*dtset%mgfft+8,2))
     ABI_ALLOCATE(rocceig,(nband_k,nband_k))
     ABI_ALLOCATE(enl0_k,(nband_k))
     ABI_ALLOCATE(enl1_k,(nband_k))

     eig1_k(:)=zero
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
     eig0_kq(:)=eigenq(1+bdtot_index:nband_k+bdtot_index)
     edocc_k(:)=zero
     eeig0_k(:)=zero ; ek0_k(:)=zero  ; ek1_k(:)=zero
     eloc0_k(:)=zero ; enl0_k(:)=zero ; enl1_k(:)=zero
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
     occ_kq(:)=occkq(1+bdtot_index:nband_k+bdtot_index)
     doccde_k(:)=doccde_rbz(1+bdtot_index:nband_k+bdtot_index)
     doccde_kq(:)=docckqde(1+bdtot_index:nband_k+bdtot_index)
     resid_k(:)=zero

!    For each pair of active bands (m,n), generates the ratios
!    rocceig(m,n)=(occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n))
!    and decide to which band to attribute it.
     call occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,dtset%occopt,occ_k,occ_kq,rocceig)

     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Get (k+q+G) wave vectors and associated spherical harmonics
     kg1_k(:,:) = 0
     kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
     call sphereboundary(gs_hamkq%gbound,istwf_k,kg1_k,dtset%mgfft,npw1_k)

     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
       end do
       if (useylmgr1==1) then
         do ilm=1,psps%mpsang*psps%mpsang
           do ii=1,3
             ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
           end do
         end do
       end if
     end if

!    Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian

!    Note that not all these arrays should be allocated in the general case
!    when wtk_k vanishes
     dimffnlk=1
     ABI_ALLOCATE(dkinpw,(npw_k))
     dkinpw(:)=zero
     ABI_ALLOCATE(kinpw1,(npw1_k))
     kinpw1(:)=zero
     ABI_ALLOCATE(ffnlk,(npw_k,dimffnlk,psps%lmnmax,1+psps%usepaw*(ntypat-1)))

!    Compute (k+G) vectors (only if useylm=1)
     nkpg=0;if(ipert>=1.and.ipert<=natom) nkpg=3*dtset%nloalg(5)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

!    Preparation of the kinetic and non-local contributions

!    === Atomic displacement perturbation ====
     if( ipert>=1 .and. ipert<=natom )then
       ider=0;idir0=0
!      Compute nonlocal form factors ffnlk at (k+G)
!      If PAW, need all ffnl(k), else only need ffnl for the displaced atom
       if (psps%usepaw==1) then
         call mkffnl(psps%dimekb,dimffnlk,psps%ekb,ffnlk,psps%ffspl,&
&         gmet,gprimd,ider,idir0,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,ntypat,&
&         psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
       else
         itypat=dtset%typat(ipert)
         ABI_ALLOCATE(ffspl_typ,(psps%mqgrid_ff,2,psps%lnmax,1))
         ffspl_typ(:,:,:,1)=psps%ffspl(:,:,:,itypat)
         call mkffnl(psps%dimekb,dimffnlk,rf_hamkq%ekb_typ(:,1,1),ffnlk,ffspl_typ,&
&         gmet,gprimd,ider,idir0,rf_hamkq%indlmn_typ,kg_k,kpg_k,kpoint,psps%lmnmax,&
&         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,1,&
&         rf_hamkq%pspso_typ,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
         ABI_DEALLOCATE(ffspl_typ)
       end if

!      === k wave vector perturbation ====
     else if(ipert==natom+1)then
       ider=1;idir0=idir
!      Compute the derivative of the kinetic operator vs k in kinpw
       call kpg3(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,idir,kg_k,kpoint,npw_k)
!      No non-local form factor is to be calculated here in this case

!      === Electric field perturbation ====
     else if(ipert==natom+2 )then
       if (psps%usepaw==1) then
         ider=1; idir0 = idir ! PAW electric field needs to apply a derivative in getgh1c
       else
         ider=0;idir0=0 ! NCPP electric field doesn't need ffnlk at all (already done in ddk)
       end if

!      === Strain perturbation ===
     else if(ipert==natom+3 .or. ipert==natom+4)then
!      istr is 1,2,...,6 and indicates cartesian strain component
       if(ipert==natom+3) then
         istr=idir
       else
         istr=idir+3
       end if
       ider=1;idir0=-istr
!      Compute the derivative of the kinetic operator vs strain in dkinpw
       call kpgstr(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,gprimd,istr,kg_k,kpoint,npw_k)

!      === Magnetic field perturbation ===
     else if(ipert==natom+5)then
!      (set equal to NCPP electric field for testing at the moment)
       ider=0;idir0=0
!      dkinpw and ffnlk are not needed ...
     end if ! ipert

!    Compute (1/2) (2 Pi)**2 (k+q+G)**2:
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg1_k,kinpw1,kpq,npw1_k)

!    Compute (k+q+G) vectors (only if useylm=1)
     nkpg1=0;if(ipert>=1.and.ipert<=natom ) nkpg1=3*dtset%nloalg(5)
     ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
     if (nkpg1>0) then
       call mkkpg(kg1_k,kpg1_k,gs_hamkq%kpoint(:),nkpg1,npw1_k)
     end if

!    Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
     dimffnl1=1+ider;if (ider==1.and.idir0==0) dimffnl1=dimffnl1+2*psps%useylm
     ABI_ALLOCATE(ffnlkq,(npw1_k,dimffnl1,psps%lmnmax,1))
     ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,ntypat))

     call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gmet,gprimd,ider,idir0,&
&     psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
&     npw1_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)

!    Atomic displ. perturbation: copy the part needed for the displaced atom, in ffnlkq.
     if( ipert>=1 .and. ipert<=natom )then
       itypat=dtset%typat(ipert)
!$OMP PARALLEL DO COLLAPSE(2)
       do ilmn=1,psps%lmnmax
         do ii=1,dimffnl1
           do ipw=1,npw1_k
             ffnlkq(ipw,ii,ilmn,1)=ffnl1(ipw,ii,ilmn,itypat)
           end do
         end do
       end do
     end if
     if (psps%usepaw==1.and.ipert==natom+2) then
!$OMP PARALLEL DO COLLAPSE(2)
       do itypat=1,ntypat
         do ilmn=1,psps%lmnmax
           do ii=1,min(dimffnlk,dimffnl1)
             do ipw=1,npw1_k
               ffnlk(ipw,ii,ilmn,itypat)=ffnl1(ipw,ii,ilmn,itypat)
             end do
           end do
         end do
       end do
     end if

!    Compute phkxred (at k and k+q) and eventually ph3d
     ABI_ALLOCATE(phkxred,(2,dimphkxred))
     if (ipert>=1.and.ipert<=natom) then
       if (psps%usepaw==0) then
         arg=two_pi*DOT_PRODUCT(kpoint,xred(:,ipert))
         phkxred(1,1)=cos(arg) ; phkxred(2,1)=sin(arg)
       else
!        If PAW and atom displ., also need all phkxred for k
         do ia=1,natom
           iatom=atindx(ia)
           arg=two_pi*DOT_PRODUCT(kpoint,xred(:,ia))
           phkxred(1,iatom)=cos(arg) ; phkxred(2,iatom)=sin(arg)
         end do
       end if
     end if

     ABI_ALLOCATE(ph3d,(2,npw1_k,gs_hamkq%matblk))
     if(dtset%nloalg(1)>0) then
!      Here, precomputation of ph3d
       call ph1d3d(1,natom,kg1_k,gs_hamkq%matblk,natom,npw1_k,n1,n2,n3,gs_hamkq%phkxred,ph1d,ph3d)
     end if

!    Compute the gradient of the Berry-phase term
     if (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
&     dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17) then
       if (ipert<=natom) then
!        phonon perturbation
         call gradberry3(cg,cg1,dtefield,grad_berry,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt_rbz,&
&         npwarr,npwar1,dtset%nspinor,nsppol,qmat,pwindall)
       else
!        electric field perturbation
         call gbefd3(cg,cg1,dtefield,grad_berry,idir,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt_rbz,&
&         npwarr,npwar1,dtset%nspinor,&
&         nsppol,qmat,pwindall,rprimd)
       end if
     end if

!    Compute the eigenvalues, wavefunction, residuals,
!    contributions to kinetic energy, nonlocal energy, forces,
!    and update of 1st-order density to this k-point and this spin polarization.
     nband_kq = nband_k  !Note that the calculation only works for same number of bandes on all K points.
!    Note that vtowfk3 is called with kpoint, while kpt is used inside vtowfk3
     call vtowfk3(cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,dimcprj,&
&     dim_eig2rf,dimffnlk,dimffnl1,dimphkxred,dkinpw,dtfil,dtset,&
&     edocc_k,eeig0_k,eig0_k,eig0_kq,eig1_k,&
&     ek0_k,ek1_k,eloc0_k,enl0_k,enl1_k,fermie1,ffnlk,ffnlkq,ffnl1,&
&     gbound,gh0c1_set,gh1c_set,grad_berry,gs_hamkq,&
&     ibg,ibgq,ibg1,icg,icgq,icg1,idir,ikpt,ipert,&
&     isppol,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpoint,&
&     mband,mcgq,mcprjq,mkmem,mk1mem,mpi_enreg,&
&     mpw,mpw1,natom,nband_k,ncpgr,nkpg,nkpg1,&
&     nnsclo_now,npw_k,npw1_k,dtset%nspinor,nsppol,ntypat,n4,n5,n6,&
&     occ_k,pawrhoij1_unsym,ph3d,phkxred,prtvol,psps,resid_k,rf_hamkq,rhoaug1,rocceig,&
&     usecprj,wffddk,vlocal,vlocal1,wtk_k,nlines_done)

!    Free temporary storage
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(kpg1_k)
     ABI_DEALLOCATE(dkinpw)
     ABI_DEALLOCATE(ffnlk)
     ABI_DEALLOCATE(ffnlkq)
     ABI_DEALLOCATE(ffnl1)
     ABI_DEALLOCATE(phkxred)
     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kinpw1)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(doccde_k)
     ABI_DEALLOCATE(doccde_kq)
     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(eig0_kq)
     ABI_DEALLOCATE(occ_kq)
     ABI_DEALLOCATE(rocceig)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)
     ABI_DEALLOCATE(ylmgr1_k)

!    Save eigenvalues (hartree), residuals (hartree**2)
     eigen1 (1+bd2tot_index : 2*nband_k**2+bd2tot_index) = eig1_k(:)
     resid  (1+bdtot_index : nband_k+bdtot_index) = resid_k(:)

!    Accumulate sum over k points for nonlocal and kinetic energies,
!    also accumulate gradients of Enonlocal:
     if (iscf_mod>0 .or. iscf_mod==-3 .or. iscf_mod==-2)then
       do iband=1,nband_k
         edocc=edocc+wtk_k*occ_k(iband)*edocc_k(iband)
         eeig0=eeig0+wtk_k*occ_k(iband)*eeig0_k(iband)
         ek0=ek0+wtk_k*occ_k(iband)*ek0_k(iband)
         ek1=ek1+wtk_k*occ_k(iband)*ek1_k(iband)
         eloc0=eloc0+wtk_k*occ_k(iband)*eloc0_k(iband)
         enl0=enl0+wtk_k*occ_k(iband)*enl0_k(iband)
         enl1=enl1+wtk_k*occ_k(iband)*enl1_k(iband)
       end do
     end if

     ABI_DEALLOCATE(eig1_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(resid_k)
     ABI_DEALLOCATE(edocc_k)
     ABI_DEALLOCATE(eeig0_k)
     ABI_DEALLOCATE(ek0_k)
     ABI_DEALLOCATE(ek1_k)
     ABI_DEALLOCATE(eloc0_k)
     ABI_DEALLOCATE(enl0_k)
     ABI_DEALLOCATE(enl1_k)

!    Keep track of total number of bands (all k points so far, even for k points not treated by me)
     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2

!    Shift array memory
     if (mkmem/=0) then
       ibg=ibg+nband_k
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if
     if (mkqmem/=0) then
       ibgq=ibgq+dtset%nspinor*nband_k
       icgq=icgq+npw1_k*dtset%nspinor*nband_k
     end if
     if (mk1mem/=0) then
       ibg1=ibg1+nband_k
       icg1=icg1+npw1_k*dtset%nspinor*nband_k
       ikg1=ikg1+npw1_k
     end if

   end do ! End big k point loop

   call timab(125,2,tsec)

!  Transfer density on augmented fft grid to normal fft grid in real space. Also take into account the spin.
   if(iscf_mod>0) then
     if (psps%usepaw==0) then
       call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rhor1  ,rhoaug1,1)
     else
       call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rho1wfr,rhoaug1,1)
     end if

   end if

 end do !  End loop over spins

!More memory cleaning
 call destroy_hamiltonian(gs_hamkq)
 call destroy_rf_hamiltonian(rf_hamkq)
 ABI_DEALLOCATE(rhoaug1)
 ABI_DEALLOCATE(vlocal)
 ABI_DEALLOCATE(vlocal1)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kg1_k)

 call timab(124,2,tsec)

!=== MPI communications ==================
 if(xmpi_paral==1)then
   call timab(129,1,tsec)

!  Compute buffer size
   buffer_size=7+mbd2kpsp+mbdkpsp
   if (iscf_mod>0) then
     buffer_size=buffer_size+cplex*dtset%nfft*nspden
   end if
   ABI_ALLOCATE(buffer1,(buffer_size))

!  Pack rhor1,edocc,eeig0,ek0,ek1,eloc0,enl0,enl1,eigen1,resid
   if (iscf_mod>0) then
     index1=cplex*dtset%nfft*nspden
     if (psps%usepaw==0) then
       buffer1(1:index1)=reshape(rhor1  ,(/index1/))
     else
       buffer1(1:index1)=reshape(rho1wfr,(/index1/))
     end if
   else
     index1=0
   end if
   buffer1(index1+1)=edocc;buffer1(index1+2)=eeig0
   buffer1(index1+3)=ek0  ;buffer1(index1+4)=ek1
   buffer1(index1+5)=eloc0;buffer1(index1+6)=enl0
   buffer1(index1+7)=enl1
   index1=index1+7
   bdtot_index=0;bd2tot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       buffer1(index1+1:index1+2*nband_k**2) = eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2)
       buffer1(index1+2*nband_k**2+1:index1+2*nband_k**2+nband_k)= resid(bdtot_index+1:bdtot_index+nband_k)
       bdtot_index=bdtot_index+nband_k
       bd2tot_index=bd2tot_index+2*nband_k**2
       index1=index1+2*nband_k**2+nband_k
     end do
   end do
   if(index1<buffer_size)buffer1(index1+1:buffer_size)=zero

!  Build sum of everything
   call timab(48,1,tsec)
   call xmpi_sum(buffer1,buffer_size,spaceworld,ierr)
   call timab(48,2,tsec)

!  Unpack the final result
   if(iscf_mod>0) then
     index1=cplex*dtset%nfft*nspden
     if (psps%usepaw==0) then
       rhor1(:,:)  =reshape(buffer1(1:index1),(/cplex*dtset%nfft,nspden/))
     else
       rho1wfr(:,:)=reshape(buffer1(1:index1),(/cplex*dtset%nfft,nspden/))
     end if
   else
     index1=0
   end if

   edocc=buffer1(index1+1);eeig0=buffer1(index1+2)
   ek0=buffer1(index1+3)  ;ek1=buffer1(index1+4)
   eloc0=buffer1(index1+5);enl0=buffer1(index1+6)
   enl1=buffer1(index1+7)
   index1=index1+7
   bdtot_index=0;bd2tot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2) = buffer1(index1+1:index1+2*nband_k**2)
       resid(bdtot_index+1:bdtot_index+nband_k)= buffer1(index1+2*nband_k**2+1:index1+2*nband_k**2+nband_k)
       bdtot_index=bdtot_index+nband_k
       bd2tot_index=bd2tot_index+2*nband_k**2
       index1=index1+2*nband_k**2+nband_k
     end do
   end do
   ABI_DEALLOCATE(buffer1)

!  Accumulate PAW occupancies
   if (psps%usepaw==1.and.iscf_mod>0) then
     call pawrhoij_mpisum_unpacked(pawrhoij1_unsym,spaceworld)
   end if

   call timab(129,2,tsec)
 end if ! if kpt parallel

 call timab(127,1,tsec)

!If needed, compute rhog1, and symmetrizes the density
 if (iscf_mod > 0) then

!  In order to have the symrhg working in parallel on FFT coefficients, the size
!  of irzzon1 and phnons1 should be set to nfftot. Therefore, nsym\=1 does not work.
   if (psps%usepaw==0) then
     call symrhg(cplex,gprimd,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,&
&     nspden,nsppol,nsym1,dtset%paral_kgb,phnons1,rhog1  ,rhor1  ,rprimd,symaf1,symrl1)
   else
     call symrhg(cplex,gprimd,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,&
&     nspden,nsppol,nsym1,dtset%paral_kgb,phnons1,rho1wfg,rho1wfr,rprimd,symaf1,symrl1)
   end if
!  We now have both rho(r) and rho(G), symmetrized, and if nsppol=2
!  we also have the spin-up density, symmetrized, in rhor1(:,2).
 end if

 ABI_DEALLOCATE(grad_berry)

!Find largest residual over bands, k points, and spins except for nbdbuf highest bands
 ibdkpt=1
 residm=zero
 do isppol=1,nsppol
   do ikpt=1,nkpt_rbz
     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     nband_k=max(1,nband_k-dtset%nbdbuf)
     residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_k-1)))
     ibdkpt=ibdkpt+nband_k
   end do
 end do

 call timab(127,2,tsec)

 if (iscf_mod>0) then

!  PAW: Build new 1st-order rhoij quantities then symetrize them
!  Compute and add the 1st-order compensation density to rho1wfr
!  to get the total 1st-order density
   if (psps%usepaw==1) then
     call pawmkrho(arg,cplex,gprimd,idir,indsy1,ipert,mpi_enreg,&
&     my_natom,natom,nspden,nsym1,ntypat,dtset%paral_kgb,pawang,pawfgr,pawfgrtab,&
&     dtset%pawprtvol,pawrhoij1,pawrhoij1_unsym,pawtab,dtset%qptn,rho1wfg,rho1wfr,&
&     rhor1,rprimd,symaf1,symrc1,dtset%typat,ucvol,dtset%usewvl,xred,&
&     pawang_sym=pawang1,pawnhat=nhat1,pawrhoij0=pawrhoij,rhog=rhog1)
     ABI_DEALLOCATE(rho1wfr)
     ABI_DEALLOCATE(rho1wfg)
     if (paral_atom) then
       call pawrhoij_destroy(pawrhoij1_unsym)
       ABI_DATATYPE_DEALLOCATE(pawrhoij1_unsym)
     end if
   end if

!  Compute density residual (if required) and its squared norm
   if (optres==1) then
     nvresid1=rhor1-nvresid1
     call sqnorm_v(1,nfftf,nres2,dtset%nspden,optres,nvresid1)
   end if
 end if ! iscf>0

 call timab(121,2,tsec)

 DBG_EXIT('COLL')

end subroutine vtorho3
!!***
