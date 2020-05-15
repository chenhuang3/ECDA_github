!{\src2tex{textfont=tt}}
!!****f* ABINIT/d2frnl
!! NAME
!! d2frnl
!!
!! FUNCTION
!! Compute the frozen-wavefunction non-local contribution for reponse functions
!! (strain and/or phonon) 
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GM, AR, MB, MT, AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnuC.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of WF
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)=<p_lmn|C> coefficients for WF |C> (and 1st derivatives)
!!  dimcprj(natom*usepaw)=array of dimensions of array cprj (ordered by atom-type)
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  dyfr_nondiag=1 if dyfrnl is non diagonal with respect to atoms; 0 otherwise
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gsqcut_eff=Fourier cutoff on G^2 for "large sphere" of radius double that of the basis sphere
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  kptns(3,nkpt)=coordinates of k points in terms of reciprocal space
!!   primitive translations
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mgfftf=maximum size of 1D FFTs for the fine FFT grid (PAW)
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimension for number of planewaves
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nband(nkpt*nsppol)=number of bands being considered per k point
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  nkpt=number of k points
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  nsym=number of symmetry elements in space group (at least 1)
!!  ntypat=integer specification of atom type (1, 2, ...)
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each k point
!!  rfphon=1   if non local contribution of dynamical matrix have to be computed
!!  rfstrs!=0  if non local contribution of elastic tensor have to be computed
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawbec= flag for the computation of Born Effective Charge within PAW ; set to 1 if yes
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  ph1df(2,3*(2*mgfftf+1)*natom)=phase information related to structure factor on the fine FFT grid (PAW)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=wavevector of the phonon
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries in reciprocal space (dimensionless)
!!  typat(natom)=type integer for each atom in cell
!!  unpaw=unit number for cprj PAW data (if used)
!!  usecprj=1 if cprj coefficients are already in memory (PAW only)
!!  vtrial(nfftf,nspden)=total potential (Hartree+XC+loc)
!!  vxc(nfftf,nspden)=XC potential
!!  wtk(nkpt)=k point weights
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,9,mpsang*mpsang*useylm)= gradients of real spherical harmonics for each G and k point
!! OUTPUT
!!  becfrnl(3,natom,3*pawbec)=NL frozen contribution to Born Effective Charges (PAW only)
!!                            (3,natom) = derivative wr to the displ. of one atom in one direction
!!                            (3)       = derivative wr to electric field in one direction
!!  dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=
!!         non-symmetrized non-local contribution to the dynamical matrix
!!         If NCPP, it depends on one atom
!!         If PAW,  it depends on two atoms
!!  eltfrnl(6+3*natom,6)=non-symmetrized non-local contribution to the
!!                    elastic tensor
!!
!! SIDE EFFECTS
!!  ===== if psps%usepaw==1
!!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!                          pawfgrtab(:)%gylmgr2 are deallocated here
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      metric,mkffnl,mkkpg,nonlop,paw_ij_destroy,paw_ij_init,paw_ij_nullify
!!      pawaccrhoij,pawcprj_alloc,pawcprj_copy,pawcprj_destroy
!!      pawcprj_diskinit_r,pawdijfr,pawfgrtab_destroy,pawfgrtab_init,pawgrnl
!!      pawrhoij_destroy,pawrhoij_gather,pawrhoij_nullify,ph1d3d,strconv
!!      symrhoij,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine d2frnl(atindx,atindx1,becfrnl,cg,cprj,dimcprj,dyfrnl,dyfr_cplex,dyfr_nondiag,eigen,eltfrnl,gsqcut,indsym,&
&          istwfk,kg,kptns,kptopt,mband,mgfft,mgfftf,mkmem,mpi_enreg,mpsang,mpw,my_natom,natom,nattyp,nband,&
&          nfftf,ngfft,ngfftf,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,nsym,ntypat,occ,&
&          paw_ij,pawang,pawbec,pawprtvol,pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,ph1df,psps,&
&          qphon,rprimd,rfphon,rfstrs,symafm,symrec,typat,unpaw,usecprj,vtrial,&
&          vxc,wtk,xred,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_wffile
 use m_header, only : hdr_skip

 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type
 use m_pawtab,   only : pawtab_type
 use m_pawfgrtab,only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_destroy
 use m_paw_ij,   only : paw_ij_type, paw_ij_init, paw_ij_destroy, paw_ij_nullify
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_copy, pawrhoij_destroy, pawrhoij_gather,&
&                       pawrhoij_nullify, symrhoij
 use m_pawcprj,  only : pawcprj_type, pawcprj_diskinit_r, pawcprj_alloc, pawcprj_get,&
&                       pawcprj_copy, pawcprj_destroy
 use m_pawdij,   only : pawdijfr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2frnl'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_65_nonlocal
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dyfr_cplex,dyfr_nondiag,kptopt,mband,mgfft,mgfftf,mkmem
 integer,intent(in) :: mpsang,mpw,my_natom,natom,nfftf,nkpt,nspden,nspinor
 integer,intent(in) :: nsppol,nsym,ntypat,pawbec,pawprtvol,rfphon,rfstrs,unpaw
 integer,intent(in) :: usecprj
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),dimcprj(natom*psps%usepaw)
 integer,intent(in) :: indsym(4,nsym,natom),istwfk(nkpt),kg(3,mpw*mkmem)
 integer,intent(in) :: nattyp(ntypat),nband(nkpt*nsppol),ngfft(18),ngfftf(18)
 integer,intent(in) :: nloalg(5),npwarr(nkpt),symafm(nsym),symrec(3,3,nsym)
 integer,intent(in) :: typat(ntypat)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*natom),qphon(3),rprimd(3,3)
 real(dp),intent(in) :: vxc(nfftf,nspden),wtk(nkpt),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,9,mpsang*mpsang*psps%useylm)
 real(dp),intent(in),target :: vtrial(nfftf,nspden)
 real(dp),intent(out) :: becfrnl(3,natom,3*pawbec)
 real(dp),intent(out) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
 real(dp),intent(out) :: eltfrnl(6+3*natom,6)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom)
 type(pawcprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(inout),target :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,bufdim,choice_bec,choice_phon,choice_strs,cplex
 integer :: cplex_dij,cplx,cpopt,cpopt_bec,dimdij,dimekb1,dimekb2,dimffnl
 integer :: dimnhat,ia,iatom,iatom_tot,iband,ibg,ibsp,icg,ider,idir,ielt,ieltx
 integer :: ierr,ii,ikg,ikpt,ilm,iorder_cprj,ipw,isp,ispden,isppol,istwf_k
 integer :: itypat,jj,klmn,master,matblk,me,mu,my_comm_atom,n1,n2,n3,nband_k
 integer :: ncpgr,ngrhoij,nkpg,nnlout_bec,nnlout_phon,nnlout_strs,npw_k,nsploop
 integer :: nu,optgr,optgr2,option,option_rhoij,optstr,optstr2,paw_opt
 integer :: paw_opt_bec,shift_rhoij,signs,spaceworld,sz2,sz3,tim_nonlop
 real(dp) :: arg,eig_k,enl,enlk,occ_k,ucvol,wtk_k
 logical :: paral_atom,usetimerev
!arrays
 integer,allocatable :: dimlmn(:),kg_k(:,:),l_size_atm(:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: dummy(0),gmet(3,3),gprimd(3,3),grhoij(3),kpoint(3),nonlop_dum(1,1)
 real(dp) :: rmet(3,3),tsec(2)
 real(dp),allocatable :: becfrnl_tmp(:,:,:),becfrnlk(:,:,:),becij(:,:,:,:)
 real(dp),allocatable :: cwavef(:,:),dyfrnlk(:,:),ekb(:,:,:),elt_work(:,:)
 real(dp),allocatable :: eltfrnlk(:,:),enlout_bec(:),enlout_phon(:)
 real(dp),allocatable :: enlout_strs(:),ffnl(:,:,:,:),kpg_k(:,:),mpibuf(:)
 real(dp),allocatable :: nhat_dum(:,:),ph3d(:,:,:),phkxred(:,:),sij(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 type(paw_ij_type),allocatable :: paw_ij_tmp(:)
 type(pawcprj_type),allocatable :: cprj_disk(:,:),cwaveprj(:,:)
 type(pawfgrtab_type),allocatable :: pawfgrtab_tmp(:)
 type(pawrhoij_type),pointer :: pawrhoij_tot(:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(159,1,tsec)

!Set up parallelism
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
 master=0
 paral_atom=(my_natom/=natom)
 my_comm_atom=mpi_enreg%comm_atom
 my_atmtab=>mpi_enreg%my_atmtab

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!PAW file
 iorder_cprj=0
 if (usecprj==1) then
   call pawcprj_diskinit_r(atindx1,natom,iorder_cprj,mkmem,natom,3,dimcprj,nspinor,unpaw)
 end if

!Initialization of frozen non local array
 if(rfphon==1) then
   dyfrnl(:,:,:,:,:)=zero
   ABI_ALLOCATE(dyfrnlk,(6,natom))
 end if

 if(rfstrs/=0)then
   eltfrnl(:,:)=zero;enl=zero
   ABI_ALLOCATE(eltfrnlk,(6+3*natom,6))
 end if
 if (pawbec==1) then
   becfrnl(:,:,:)=zero
   ABI_ALLOCATE(becfrnlk,(3,natom,3))
 end if

!Common initialization
 bdtot_index=0;ibg=0;icg=0
 nsploop=nsppol;if (nspden==4) nsploop=4
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(cwavef,(2,mpw*nspinor))
 ABI_ALLOCATE(phkxred,(2,natom))

!Common data for "nonlop" routine
 signs=1 ; idir=0 ; eig_k=zero
 choice_phon=0;choice_strs=0
 if(rfphon==1)then
   shift_rhoij=0
   choice_phon=4
   nnlout_phon=max(1,6*natom)
   tim_nonlop=4
   ABI_ALLOCATE(enlout_phon,(nnlout_phon))
 end if
 if(rfstrs/=0)then
   shift_rhoij=6
   choice_strs=6
   nnlout_strs=6*(3*natom+6)
   tim_nonlop=6
   ABI_ALLOCATE(enlout_strs,(nnlout_strs))
 end if
 if (psps%usepaw==0) then
   paw_opt=0 ; cpopt=-1
 else
   paw_opt=2 ; cpopt=1+2*usecprj
 end if
 if (pawbec==1) then
   choice_bec=2 ; nnlout_bec=max(1,3*natom) ; paw_opt_bec=1 ; cpopt_bec=4
   ABI_ALLOCATE(enlout_bec,(nnlout_bec))
 else
   choice_bec=0 ; nnlout_bec=0 ; paw_opt_bec=0 ; cpopt_bec=0
 end if

!===== Norm-conserving PSPS
 if (psps%usepaw==0) then   !NCPP

!  Non-local factors: kleimann-Bylander energies
   dimekb1=psps%dimekb;dimekb2=ntypat
   ABI_ALLOCATE(ekb,(psps%dimekb,ntypat,nspinor**2))
   ekb(:,:,1)=psps%ekb(:,:)
   if (nspinor==2) then
     ekb(:,:,2)=psps%ekb(:,:)
     ekb(:,:,3:4)=zero
   end if
   ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
   ABI_ALLOCATE(sij,(0,0))
 else

!===== PAW

!  Define several sizes & flags
   ncpgr=0;ngrhoij=0
   if(rfstrs/=0)then
     ncpgr=9;ngrhoij=9
   end if
   if(rfphon==1)then
     ncpgr=3;ngrhoij=3
   end if
   if(rfstrs/=0.and.rfphon==1)then
     ncpgr=9;ngrhoij=9
   end if

!  Non-local factors: Dij coefficients and overlap coefficients
   cplex_dij=max(1,nspinor)
   dimekb1=psps%dimekb*cplex_dij;dimekb2=natom
   ABI_ALLOCATE(ekb,(dimekb1,dimekb2,nspinor**2))
   ABI_ALLOCATE(sij,(dimekb1,ntypat))
   ekb(:,:,:)=zero;sij(:,:)=zero
   do itypat=1,ntypat
     if (cplex_dij==1) then
       sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
     else
       do klmn=1,pawtab(itypat)%lmn2_size
         sij(2*klmn-1,itypat)=pawtab(itypat)%sij(klmn)
       end do
     end if
   end do

!  If PAW and Born Eff., Charges has to compute some additional data:
!  For each atom and for electric field direction k:
!  becij(k)=<Phi_i|r_k-R_k|Phi_j>-<tPhi_i|r_k-R_k|tPhi_j> + sij.R_k
   if (pawbec==1) then
     ABI_ALLOCATE(becij,(psps%dimekb,dimekb2,nspinor**2,3))
     becij=zero
     ABI_DATATYPE_ALLOCATE(paw_ij_tmp,(my_natom))
     ABI_DATATYPE_ALLOCATE(pawfgrtab_tmp,(my_natom))
     call paw_ij_nullify(paw_ij_tmp)
     call paw_ij_init(paw_ij_tmp,1,1,1,1,0,natom,ntypat,typat,pawtab,has_dijfr=1,&
&     mpi_comm_atom=my_comm_atom,mpi_atmtab=my_atmtab )
     ABI_ALLOCATE(l_size_atm,(my_natom))
     do iatom=1,my_natom
       iatom_tot=iatom; if(paral_atom) iatom_tot=my_atmtab(iatom)
       itypat=typat(iatom_tot)
       l_size_atm(iatom)=pawtab(itypat)%lcut_size
     end do
     call pawfgrtab_init(pawfgrtab_tmp,1,l_size_atm,nspden,typat,&
&     mpi_atmtab=my_atmtab,mpi_comm_atom=my_comm_atom)
     ABI_DEALLOCATE(l_size_atm)
     do ii=1,3 ! Loop over direction of electric field
       call pawdijfr(1,gprimd,ii,natom+2,my_natom,natom,nfftf,ngfftf,nspden,ntypat,&
&       0,paw_ij_tmp,pawang,pawfgrtab_tmp,pawrad,pawtab,&
&       (/zero,zero,zero/),rprimd,ucvol,vtrial,vtrial,vxc,xred,&
&       mpi_comm_atom=my_comm_atom, mpi_atmtab=my_atmtab ) ! vtrial not used here
       do iatom=1,my_natom
         iatom_tot=iatom; if(paral_atom) iatom_tot=my_atmtab(iatom)
         itypat=typat(iatom_tot);dimdij=pawtab(itypat)%lmn2_size
!        Add contribution from Phi_i/Phi_j and tPhi_i/tPhi_j to becij
         do klmn=1,dimdij
           becij(klmn,iatom_tot,1,ii)=paw_ij_tmp(iatom)%dijfr(klmn,1)
         end do
!        Add contribution from sij to becij
!        xred are atomic positions in reduced cooridinates of real space
!        need them in reduced coordinates of reciprocal space
         arg=gmet(ii,1)*xred(1,iatom_tot)+gmet(ii,2)*xred(2,iatom_tot)+gmet(ii,3)*xred(3,iatom_tot)
         if (cplex_dij==1) then
           becij(1:dimdij,iatom_tot,1,ii)=becij(1:dimdij,iatom_tot,1,ii)+arg*sij(1:dimdij,itypat)
         else
           do klmn=1,dimdij
             becij(klmn,iatom_tot,1,ii)=becij(dimdij,iatom_tot,1,ii)+arg*sij(2*klmn-1,itypat)
           end do
         end if
       end do
     end do
     if (nspinor==2) becij(:,:,2,:)=becij(:,:,1,:)
     if (paral_atom) then
       call xmpi_sum(becij,my_comm_atom,ierr)
     end if
     call paw_ij_destroy(paw_ij_tmp)
     call pawfgrtab_destroy(pawfgrtab_tmp)
     ABI_DATATYPE_DEALLOCATE(paw_ij_tmp)
     ABI_DATATYPE_DEALLOCATE(pawfgrtab_tmp)
   end if

!  PAW occupancies: need to communicate when paral atom is activated
   if (paral_atom) then
     ABI_DATATYPE_ALLOCATE(pawrhoij_tot,(natom))
     call pawrhoij_nullify(pawrhoij_tot)
     call pawrhoij_gather(pawrhoij,pawrhoij_tot,-1,my_comm_atom)
   else
     pawrhoij_tot => pawrhoij
   end if

!  Projected WF (cprj) and PAW occupancies (& gradients)
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,nspinor))
   call pawcprj_alloc(cwaveprj,ncpgr,dimcprj)
   do iatom=1,natom
     sz2=pawrhoij_tot(iatom)%cplex*pawrhoij_tot(iatom)%lmn2_size
     sz3=pawrhoij_tot(iatom)%nspden
     ABI_ALLOCATE(pawrhoij_tot(iatom)%grhoij,(ngrhoij,sz2,sz3))
     pawrhoij_tot(iatom)%ngrhoij=ngrhoij
     pawrhoij_tot(iatom)%grhoij=zero
   end do
   usetimerev=(kptopt>0.and.kptopt<3)

 end if !PAW

!LOOP OVER SPINS
 do isppol=1,nsppol

   ikg=0

!  PAW: retrieve Dij coefficients
   if (psps%usepaw==1) then
     ekb(:,:,:)=zero
     if (my_natom>0) then
       do ispden=1,nspinor**2
         isp=isppol;if (nspinor==2) isp=ispden
         do iatom=1,my_natom
           iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
           dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
           do klmn=1,dimdij
             ekb(klmn,iatom_tot,ispden)=paw_ij(iatom)%dij(klmn,isp)
           end do
         end do
       end do
     end if
     if (paral_atom) then
       call xmpi_sum(ekb,my_comm_atom,ierr)
     end if
   end if

!  Loop over k points
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)
     wtk_k=wtk(ikpt)

!    Skip this k-point if not the proper processor
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     if(rfstrs/=0)then
       ABI_ALLOCATE(ylmgr_k,(npw_k,9,mpsang*mpsang*psps%useylm))
     else
       ABI_ALLOCATE(ylmgr_k,(0,0,0))
     end if

     kpoint(:)=kptns(:,ikpt)
     kg_k(:,:) = 0

!$OMP PARALLEL DO
     do ipw=1,npw_k
       kg_k(1,ipw)=kg(1,ipw+ikg)
       kg_k(2,ipw)=kg(2,ipw+ikg)
       kg_k(3,ipw)=kg(3,ipw+ikg)
     end do
     if (psps%useylm==1) then
!SOMP PARALLEL DO COLLAPSE(2)
       do ilm=1,mpsang*mpsang
         do ipw=1,npw_k
           ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
         end do
       end do
       if(rfstrs/=0)then
!SOMP PARALLEL DO COLLAPSE(3)
         do ilm=1,mpsang*mpsang
           do ii=1,9
             do ipw=1,npw_k
               ylmgr_k(ipw,ii,ilm)=ylmgr(ipw+ikg,ii,ilm)
             end do
           end do
         end do
       end if
     end if

     cplex=2;if (istwf_k>1) cplex=1

!    Compute nonlocal psp energy

!    Compute (k+G) vectors (only if useylm=1)     
     if(rfstrs/=0) nkpg=3*nloalg(5)
     if (rfphon==1) nkpg=9*nloalg(5)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

!    Compute nonlocal form factors ffnl at all (k+G):
     ider=0;dimffnl=1
     if(rfstrs/=0)then
       ider=2;dimffnl=3+7*psps%useylm 
     end if
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,&
&     ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

!    Compute phkxred and eventually ph3d.
     do iatom=1,natom
       ia=atindx1(iatom)
       arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
       phkxred(1,iatom)=cos(arg)
       phkxred(2,iatom)=sin(arg)
     end do
     if(nloalg(1)<=0)then
!      Only the allocation, not the precomputation.
       matblk=nloalg(4)
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
     else
!      Here, allocation as well as precomputation
       matblk=natom
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
       call ph1d3d(1,natom,kg_k,matblk,natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)
     end if

!    Initialize contributions from current k point
     if(rfphon==1) dyfrnlk(:,:)=zero
     if(rfstrs/=0)then
       enlk=zero;eltfrnlk(:,:)=zero
     end if
     if (pawbec==1) becfrnlk(:,:,:)=zero

!    Loop over bands
     do iband=1,nband_k

       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,isppol,me)) cycle

       occ_k=occ(iband+bdtot_index)

       cwavef(:,1:npw_k*nspinor) = cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)

       if (psps%usepaw==1.and.usecprj==1) then
         ibsp=(iband-1)*nspinor+ibg
         call pawcprj_copy(cprj     (:,ibsp+1:ibsp+nspinor),cwaveprj)
       end if

!      Compute non-local contributions from n,k
       if (psps%usepaw==1) eig_k=eigen(iband+bdtot_index)

!      === Dynamical matrix
       if(rfphon==1) then
         call nonlop(atindx1,choice_phon,cpopt,cwaveprj,dimekb1,dimekb2,dimffnl,dimffnl,ekb,&
&         enlout_phon,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&         kpoint,(/eig_k/),psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,nattyp,1,&
&         ngfft,nkpg,nkpg,nloalg,nnlout_phon,npw_k,npw_k,nspinor,nspinor,ntypat,0,paw_opt,phkxred,&
&         phkxred,ph1d,ph3d,ph3d,signs,sij,nonlop_dum,4,ucvol,&
&         psps%useylm,cwavef,cwavef)
!        Accumulate non-local contributions from n,k
         dyfrnlk(:,:)=dyfrnlk(:,:)+occ_k*reshape(enlout_phon(:),(/6,natom/))
       end if

!      === Elastic tensor
       if(rfstrs/=0) then
         call nonlop(atindx1,choice_strs,cpopt,cwaveprj,dimekb1,dimekb2,dimffnl,dimffnl,ekb,&
&         enlout_strs,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&         kpoint,(/eig_k/),psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,nattyp,1,&
&         ngfft,nkpg,nkpg,nloalg,nnlout_strs,npw_k,npw_k,nspinor,nspinor,ntypat,0,paw_opt,phkxred,&
&         phkxred,ph1d,ph3d,ph3d,signs,sij,nonlop_dum,6,ucvol,&
&         psps%useylm,cwavef,cwavef)
!        Accumulate non-local contributions from n,k
         eltfrnlk(:,:)=eltfrnlk(:,:)+occ_k*reshape(enlout_strs(:),(/3*natom+6,6/))
!        PAW: accumulate gradients of rhoij
       end if !endxo if strs

       if (psps%usepaw==1) then 
         call pawaccrhoij(atindx,cplex,cwaveprj,cwaveprj,0,isppol,natom,&
&         natom,nspinor,occ_k,3,pawrhoij_tot,usetimerev,wtk_k)
       end if
       

!      PAW: Compute frozen contribution to Born Effective Charges
       if (pawbec==1) then
         do ii=1,3 ! Loop over elect. field directions
           call nonlop(atindx1,choice_bec,cpopt_bec,cwaveprj,psps%dimekb,dimekb2,dimffnl,dimffnl,&
&           becij(:,:,:,ii),enlout_bec,ffnl,ffnl,gmet,gprimd,ii,psps%indlmn,istwf_k,kg_k,kg_k,&
&           kpg_k,kpg_k,kpoint,kpoint,(/eig_k/),psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,&
&           natom,nattyp,1,ngfft,nkpg,nkpg,nloalg,nnlout_bec,npw_k,npw_k,nspinor,nspinor,ntypat,0,paw_opt_bec,phkxred,&
&           phkxred,ph1d,ph3d,ph3d,signs,sij,nonlop_dum,tim_nonlop,ucvol,&
&           psps%useylm,cwavef,cwavef)
           becfrnlk(:,:,ii)=becfrnlk(:,:,ii)+occ_k*reshape(enlout_bec(:),(/3,natom/))
         end do
       end if
     end do ! End of loop on bands

     if(rfphon==1) then
       do iatom=1,natom
         ia=iatom;if (dyfr_nondiag==0) ia=1
         dyfrnl(1,1,1,iatom,ia)=dyfrnl(1,1,1,iatom,ia)+wtk_k*dyfrnlk(1,iatom)
         dyfrnl(1,2,2,iatom,ia)=dyfrnl(1,2,2,iatom,ia)+wtk_k*dyfrnlk(2,iatom)
         dyfrnl(1,3,3,iatom,ia)=dyfrnl(1,3,3,iatom,ia)+wtk_k*dyfrnlk(3,iatom)
         dyfrnl(1,2,3,iatom,ia)=dyfrnl(1,2,3,iatom,ia)+wtk_k*dyfrnlk(4,iatom)
         dyfrnl(1,1,3,iatom,ia)=dyfrnl(1,1,3,iatom,ia)+wtk_k*dyfrnlk(5,iatom)
         dyfrnl(1,1,2,iatom,ia)=dyfrnl(1,1,2,iatom,ia)+wtk_k*dyfrnlk(6,iatom)
       end do
     end if ! end if rfphon
     if(rfstrs/=0)then
       eltfrnl(:,:)=eltfrnl(:,:)+wtk(ikpt)*eltfrnlk(:,:)
     end if
     if(pawbec==1)then
       becfrnl(:,:,:)=becfrnl(:,:,:)+wtk(ikpt)*becfrnlk(:,:,:)
     end if

!    Increment indexes
     bdtot_index=bdtot_index+nband_k
     if (mkmem/=0) then
       ibg=ibg+nband_k*nspinor
       icg=icg+npw_k*nspinor*nband_k
       ikg=ikg+npw_k
     else if (usecprj==1) then
       call pawcprj_destroy(cprj_disk)
       ABI_DATATYPE_DEALLOCATE(cprj_disk)
     end if

     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(ylm_k)
     !if(rfstrs/=0)then
     ABI_DEALLOCATE(ylmgr_k)
     !end if

   end do ! End loops on isppol and ikpt
 end do

 ABI_DEALLOCATE(phkxred)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(ekb)
 if(rfphon==1) then
   ABI_DEALLOCATE(dyfrnlk)
   ABI_DEALLOCATE(enlout_phon)
 end if

 if(rfstrs/=0) then
   ABI_DEALLOCATE(eltfrnlk)
   ABI_DEALLOCATE(enlout_strs)
 end if
 if (pawbec==1)  then
   ABI_DEALLOCATE(becfrnlk)
   ABI_DEALLOCATE(enlout_bec)
   ABI_DEALLOCATE(becij)
 end if
 if (psps%usepaw==1) then
   call pawcprj_destroy(cwaveprj)
 end if
 ABI_DEALLOCATE(sij)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

!Fill in lower triangle of matrixes
 if(rfphon==1)then
   do iatom=1,natom
     ia=iatom;if (dyfr_nondiag==0) ia=1
     dyfrnl(1,3,2,iatom,ia)=dyfrnl(1,2,3,iatom,ia)
     dyfrnl(1,3,1,iatom,ia)=dyfrnl(1,1,3,iatom,ia)
     dyfrnl(1,2,1,iatom,ia)=dyfrnl(1,1,2,iatom,ia)
   end do
 end if
 if(rfstrs/=0)then
   do jj=2,6
     do ii=1,jj-1
       eltfrnl(jj,ii)=eltfrnl(ii,jj)
     end do
   end do
 end if


!Parallel case: accumulate (n,k) contributions
 if (xmpi_paral==1) then
   call timab(48,1,tsec)
!  Accumulate dyfrnl
   if(rfphon==1)then
     call xmpi_sum(dyfrnl,spaceworld,ierr)
   end if
!  Accumulate eltfrnl.
   if(rfstrs/=0)then
     call xmpi_sum(eltfrnl,spaceworld,ierr)
   end if
!  Accumulate becfrnl
   if (pawbec==1) then
     call xmpi_sum(becfrnl,spaceworld,ierr)
   end if
!  PAW: accumulate gradients of rhoij
   if (psps%usepaw==1) then
     ABI_ALLOCATE(dimlmn,(natom))
     dimlmn(1:natom)=pawrhoij_tot(1:natom)%cplex*pawrhoij_tot(1:natom)%lmn2_size
     bufdim=ncpgr*sum(dimlmn)*nsploop
     ABI_ALLOCATE(mpibuf,(bufdim))
     ii=0;mpibuf=zero
     do iatom=1,natom
       do isppol=1,nsploop
         do mu=1,ncpgr
           mpibuf(ii+1:ii+dimlmn(iatom))=pawrhoij_tot(iatom)%grhoij(mu,1:dimlmn(iatom),isppol)
           ii=ii+dimlmn(iatom)
         end do
       end do
     end do
     call xmpi_sum(mpibuf,spaceworld,ierr)
     ii=0
     do iatom=1,natom
       do isppol=1,nsploop
         do mu=1,ncpgr
           pawrhoij_tot(iatom)%grhoij(mu,1:dimlmn(iatom),isppol)=mpibuf(ii+1:ii+dimlmn(iatom))
           ii=ii+dimlmn(iatom)
         end do
       end do
     end do
     ABI_DEALLOCATE(mpibuf)
     ABI_DEALLOCATE(dimlmn)
   end if
   call timab(48,2,tsec)
 end if

!The indexing array atindx is used to reestablish the correct order of atoms
 if(rfstrs/=0)then
   ABI_ALLOCATE(elt_work,(6+3*natom,6))
   elt_work(1:6,1:6)=eltfrnl(1:6,1:6)
   do ia=1,natom
     ielt=7+3*(ia-1)
     ieltx=7+3*(atindx(ia)-1)
     elt_work(ielt:ielt+2,1:6)=eltfrnl(ieltx:ieltx+2,1:6)
   end do
   eltfrnl(:,:)=elt_work(:,:)
   ABI_DEALLOCATE(elt_work)
 end if

!====== PAW: Addiditonal steps

 if (psps%usepaw==1) then
!  Symmetrize rhoij gradients and transfer to cartesian (reciprocal space) coord.
!  This symetrization is necessary in the antiferromagnetic case...
   if (rfphon==1.and.rfstrs==0) then
     option_rhoij=2;option=0
     call symrhoij(pawrhoij_tot,pawrhoij_tot,option_rhoij,gprimd,indsym,0,natom,nsym,&
&     ntypat,option,pawang,pawprtvol,pawtab,rprimd,symafm,symrec,typat,&
&     mpi_comm_atom=my_comm_atom, mpi_atmtab=my_atmtab)
   else if (rfphon==1.and.rfstrs==1) then
     option_rhoij=23;option=0
     call symrhoij(pawrhoij_tot,pawrhoij_tot,option_rhoij,gprimd,indsym,0,natom,nsym,&
&     ntypat,option,pawang,pawprtvol,pawtab,rprimd,symafm,symrec,typat,&
&     mpi_comm_atom=my_comm_atom, mpi_atmtab=my_atmtab)
   end if

!  Coordinates translation
   do iatom=1,natom
     cplx=pawrhoij_tot(iatom)%cplex
     do isppol=1,nsploop
       do klmn=1,pawrhoij_tot(iatom)%lmn2_size
         do ii=1,cplx
           if(rfphon==1.or.rfstrs/=0)then
             grhoij(1:3)=pawrhoij_tot(iatom)%grhoij(shift_rhoij+1:shift_rhoij+3,cplx*(klmn-1)+ii,isppol)
             do mu=1,3
               pawrhoij_tot(iatom)%grhoij(shift_rhoij+mu,cplx*(klmn-1)+ii,isppol)=gprimd(mu,1)*grhoij(1)&
&               +gprimd(mu,2)*grhoij(2)+gprimd(mu,3)*grhoij(3)
             end do
           end if
           if(rfstrs/=0)then
             call strconv(pawrhoij_tot(iatom)%grhoij(1:6,cplx*(klmn-1)+ii,isppol),gprimd,&
&             pawrhoij_tot(iatom)%grhoij(1:6,cplx*(klmn-1)+ii,isppol))
           end if
         end do
       end do
     end do
   end do

!  In case of elastic tensor computation, add diagonal contribution:
!     -delta_{alphabeta} rhoi_{ij} to drhoij/d_eps
   if(rfstrs/=0)then
     do iatom=1,natom
       cplx=pawrhoij_tot(iatom)%cplex
       do isppol=1,nsploop
         do nu=1,pawrhoij_tot(iatom)%nrhoijsel
           klmn=pawrhoij_tot(iatom)%rhoijselect(nu)
           do ii=1,cplx
             pawrhoij_tot(iatom)%grhoij(1:3,cplx*(klmn-1)+ii,isppol)= &
&             pawrhoij_tot(iatom)%grhoij(1:3,cplx*(klmn-1)+ii,isppol)&
&             -pawrhoij_tot(iatom)%rhoijp(cplx*(nu-1)+ii,isppol)
           end do
         end do
       end do
     end do
   end if

!  Add gradients due to Dij derivatives to dynamical matrix/stress tensor
   dimnhat=0;optgr=0;optgr2=0;optstr=0;optstr2=0
   if (rfphon==1) optgr2=1
   if (rfstrs/=0) optstr2=1
   ABI_ALLOCATE(nhat_dum,(1,0))
   call pawgrnl(atindx1,dimnhat,dyfrnl,dyfr_cplex,eltfrnl,dummy,gsqcut,mgfftf,my_natom,natom,&
&   nattyp,nfftf,ngfftf,nhat_dum,dummy,nspden,nsym,ntypat,optgr,optgr2,optstr,optstr2,&
&   pawang,pawfgrtab,pawrhoij_tot,pawtab,ph1df,psps,qphon,rprimd,symrec,typat,vtrial,vxc,xred,&
&   mpi_atmtab=my_atmtab,mpi_comm_atom=my_comm_atom)
   ABI_DEALLOCATE(nhat_dum)

 end if !PAW

!Born Effective Charges and PAW:
!1-Re-order atoms -- 2-Add contribution from rhoij
 if (pawbec==1) then
   ABI_ALLOCATE(becfrnl_tmp,(3,natom,3))
   becfrnl_tmp=becfrnl
   do ia=1,natom         ! Atom (sorted by type)
     iatom=atindx1(ia)   ! Atom (not sorted)
     itypat=typat(iatom)
     arg=zero            ! Computation of Sum [Rhoij.Sij]
!    do isp=1,nsppol
!      do mu=1,pawrhoij_tot(iatom)%nrhoijsel
!        klmn=pawrhoij_tot(iatom)%rhoijselect(mu)
!        arg=arg+pawrhoij_tot(iatom)%rhoijp(mu,isp)*pawtab(itypat)%sij(klmn)
!      end do
!    end do
     do ii=1,3           ! Direction of electric field
       do jj=1,3         ! Direction of atom
         becfrnl(jj,iatom,ii)=becfrnl_tmp(jj,ia,ii) !+arg*gmet(ii,jj)
       end do
     end do
   end do
   ABI_DEALLOCATE(becfrnl_tmp)
!  For testing purpose
   becfrnl=zero
 end if

!PAW: release now useless memory
 if (psps%usepaw==1) then
   if (paral_atom) then
     call pawrhoij_destroy(pawrhoij_tot)
     ABI_DATATYPE_DEALLOCATE(pawrhoij_tot)
   else
     do iatom=1,my_natom
       ABI_DEALLOCATE(pawrhoij(iatom)%grhoij)
       pawrhoij(iatom)%ngrhoij=0
     end do
   end if
 end if

 call timab(159,2,tsec)

 DBG_EXIT("COLL")

end subroutine d2frnl
!!***
