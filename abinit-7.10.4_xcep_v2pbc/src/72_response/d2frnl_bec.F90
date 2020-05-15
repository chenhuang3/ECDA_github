!{\src2tex{textfont=tt}}
!!****f* ABINIT/d2frnl_bec
!! NAME
!! d2frnl_bec
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
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  dyfr_nondiag=1 if dyfrnl is non diagonal with respect to atoms; 0 otherwise
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gsqcut_eff=Fourier cutoff on G^2 for "large sphere" of radius double that of the basis sphere
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!   primitive translations
!!  mgfftf=maximum size of 1D FFTs for the fine FFT grid (PAW)
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngs_rbzfftf=ngfft for norm-conserving potential runs)
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  ntypat=integer specification of atom type (1, 2, ...)
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each k point
!!  rfphon=1   if non local contribution of dynamical matrix have to be computed
!!  rfstrs!=0  if non local contribution of elastic tensor have to be computed
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawbec= flag for the computation of Born Effective Charge within PAW ; set to 1 if yes
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  ph1df(2,3*(2*mgfftf+1)*natom)=phase information related to structure factor on the fine FFT grid (PAW)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space (dimensionless)
!!  usecprj=1 if cprj coefficients are already in memory (PAW only)
!!  vtrial(nfftf,nspden)=total potential (Hartree+XC+loc)
!!  vxc(nfftf,nspden)=XC potential
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
!!      pawdijfr,pawfgrtab_destroy,pawfgrtab_init,pawgrnl
!!      pawrhoij_destroy,pawrhoij_gather,pawrhoij_nullify,ph1d3d,strconv
!!      symrhoij,timab,xmpi_sum
!!
!! SOURCE
#include "abi_common.h"

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine d2frnl_bec(atindx,atindx1,becfrnl,cg,cprj,dimcprj,dtfil,dtset,dyfrnl,dyfr_cplex,dyfr_nondiag,eigen,&
&          eltfrnl,gsqcut,indsym,kg,mgfftf,mpi_enreg,mpsang,my_natom,natom,nattyp,nfftf,ngfft,ngfftf,npwarr,&
&          occ,paw_ij,pawang,pawbec,pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,ph1df,psps,&
&          rprimd,rfphon,rfstrs,symrec,usecprj,vtrial,vxc,xred,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_cgtools
 use m_wffile
 use m_header,   only : hdr_skip
 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type
 use m_pawtab,   only : pawtab_type
 use m_pawfgrtab,only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_destroy
 use m_paw_ij,   only : paw_ij_type, paw_ij_init, paw_ij_destroy, paw_ij_nullify,&
&                       paw_ij_reset_flags                       
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_copy, pawrhoij_destroy, pawrhoij_gather,&
&                       pawrhoij_nullify, symrhoij
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_get,&
&                       pawcprj_copy, pawcprj_destroy
 use m_pawdij,   only : pawdijfr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2frnl_bec'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_62_iowfdenpot
 use interfaces_65_nonlocal
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dyfr_cplex,dyfr_nondiag,mgfftf,mpsang,my_natom,natom
 integer,intent(in) :: nfftf,pawbec,rfphon,rfstrs,usecprj
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),dimcprj(natom*psps%usepaw)
 integer,intent(in) :: indsym(4,dtset%nsym,natom),kg(3,dtset%mpw*dtset%mkmem)
 integer,intent(in) :: nattyp(psps%ntypat),ngfft(18),ngfftf(18),npwarr(dtset%nkpt)
 integer,intent(in) :: symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*natom),rprimd(3,3)
 real(dp),intent(in) :: vxc(nfftf,dtset%nspden),xred(3,natom)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,9,mpsang*mpsang*psps%useylm)
 real(dp),intent(in),target :: vtrial(nfftf,dtset%nspden)
 real(dp),intent(out) :: becfrnl(3,natom,3*pawbec)
 real(dp),intent(out) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
 real(dp),intent(out) :: eltfrnl(6+3*natom,6)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom)
 type(pawcprj_type) :: cprj(natom,dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*usecprj)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat)
 type(pawrhoij_type),intent(inout),target :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,bufdim,choice_bec2,choice_bec54,choice_phon,choice_strs
 integer :: cplex,cplex_dij,cplx,cpopt,cpopt_bec,ddkcase,dimdij,dimekb1,dimekb2,dimekb2_der
 integer :: dimffnl,dimffnl_bec,dimnhat,ia,iatom,iashift,iatom_tot,iband,ibg,ibsp,icg,ider,idir
 integer :: ider_bec,idir_bec,ielt,ieltx,ierr,ii,ikg,ikpt,ikpt_,ilm,ipw,isp,ispden
 integer :: isppol,ispinor,istwf_k,itypat,jj,klmn,master,matblk,matblk_der,me,mu
 integer :: my_comm_atom,n1,n2,n3,natom_der,nband_,nband_k,ncpgr,nfftot,nskip,ngrhoij,nkpg,nnlout_bec1,nnlout_bec2
 integer :: nnlout_phon,nnlout_strs,npw_,npw_k,ntypat_der,nspinor_,nsploop,nu
 integer :: optgr,optgr2,option,option_rhoij,optstr,optstr2,paw_opt,paw_opt_bec1,paw_opt_bec3
 integer :: shift_rhoij,shift1,shift2,shift3,signs,signs_bec1,signs_bec2,spaceworld,sz2,sz3,t_iostat,tim_nonlop
 real(dp) :: arg,eig_k,enl,enlk,occ_k,ucvol,wtk_k
 logical :: has_ddk_file,need_becfr,paral_atom,t_test,usetimerev
 character(len=500) :: msg
!arrays
 integer :: atindx1_der(1),dimcprj_der(1)
 integer :: ddkfil(3),ikpt_fbz(3),ikpt_fbz_previous(3),nattyp_der(1),nloalg_der(5),skipddk(3)
 integer,allocatable :: dimlmn(:),indlmn_der(:,:,:),kg_k(:,:),l_size_atm(:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: dummy(0),gmet(3,3),gprimd(3,3),grhoij(3),kpoint(3),phkxredin(2,1),phkxredout(2,1),nonlop_dum(1,1)
 real(dp) :: dotprod(2),rmet(3,3),tsec(2),xred_der(3)
 real(dp),allocatable :: becfrnl_tmp(:,:,:),becfrnlk(:,:,:),becij(:,:,:,:)
 real(dp),allocatable :: cwavef(:,:),ddk(:,:),dyfrnlk(:,:),ekb(:,:,:),ekb_der(:,:,:)
 real(dp),allocatable :: elt_work(:,:),eltfrnlk(:,:),enlout_bec1(:),enlout_bec2(:)
 real(dp),allocatable :: enlout_phon(:),enlout_strs(:),ffnl(:,:,:,:),ffnl_bec(:,:,:,:),ffnlk_der(:,:,:,:),kpg_k(:,:)
 real(dp),allocatable :: mpibuf(:),nhat_dum(:,:),ph1d_der(:,:),ph3d(:,:,:),ph3din(:,:,:),ph3dout(:,:,:),phkxred(:,:)
 real(dp),allocatable :: sij(:,:),sij_der(:,:),svectout(:,:),ylm_k(:,:),ylmgr_k(:,:,:)
 character(len=fnlen) :: fiwfddk(3)
 type(paw_ij_type),allocatable :: paw_ij_tmp(:)
 type(pawcprj_type),allocatable :: cprj_der(:,:)
 type(pawcprj_type),allocatable,target :: cwaveprj(:,:)
 type(pawfgrtab_type),allocatable :: pawfgrtab_tmp(:)
 type(pawrhoij_type),pointer :: pawrhoij_tot(:)
 type(wffile_type) :: wffddk(3)

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

!If needed, check for ddk files (used for effective charges)
 if (pawbec==1) then
   ddkfil(:)=0
   do idir=1,3
     ddkcase=idir+natom*3
     call appdig(ddkcase,dtfil%fnamewffddk,fiwfddk(idir))
     inquire(file=fiwfddk(idir),iostat=t_iostat,exist=t_test)
     if (t_iostat/=0) then
       write(msg,fmt='(5a,i8)') 'Check for existence of file ',trim(fiwfddk(idir)),',',ch10,&
&                               'but INQUIRE statement returns error code',t_iostat
       MSG_ERROR(msg)
     else if (t_test) then
       ddkfil(idir)=20+idir ! Note the use of unit numbers 21, 22 and 23
     end if
   end do
   has_ddk_file=(any(ddkfil(:)>0))
 else
   has_ddk_file=.FALSE.
 end if
 need_becfr=(pawbec==1.and.has_ddk_file)

!Initialization of frozen non local array
 if(rfphon==1) then
   dyfrnl(:,:,:,:,:)=zero
   ABI_ALLOCATE(dyfrnlk,(6,natom))
 end if
 if(rfstrs/=0)then
   eltfrnl(:,:)=zero;enl=zero
   ABI_ALLOCATE(eltfrnlk,(6+3*natom,6))
 end if
 if (need_becfr) then
   becfrnl(:,:,:)=zero
   ABI_ALLOCATE(becfrnlk,(3,natom,3))
 end if

!Common initialization
 bdtot_index=0;ibg=0;icg=0
 nsploop=dtset%nsppol;if (dtset%nspden==4) nsploop=4
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 nfftot=ngfftf(1)*ngfftf(2)*ngfftf(3)
 ABI_ALLOCATE(kg_k,(3,dtset%mpw))
 ABI_ALLOCATE(cwavef,(2,dtset%mpw*dtset%nspinor))
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
 if (need_becfr) then
   choice_bec2=2 ; choice_bec54=54
   signs_bec1=1 ; signs_bec2=2
   nnlout_bec1=max(1,3*natom) ; nnlout_bec2=max(1,18*natom); 
   paw_opt_bec1=1 ; paw_opt_bec3=3 ; cpopt_bec=-1
   ABI_ALLOCATE(enlout_bec1,(nnlout_bec1))
   ABI_ALLOCATE(enlout_bec2,(nnlout_bec2))
 else
   choice_bec2=0 ; choice_bec54=0; signs_bec1=0; signs_bec2=0
   nnlout_bec1=0 ; paw_opt_bec1=0 ; paw_opt_bec3=0 ;cpopt_bec=0
 end if

!===== Norm-conserving PSPS
 if (psps%usepaw==0) then

!  Non-local factors: kleimann-Bylander energies
   dimekb1=psps%dimekb;dimekb2=psps%ntypat
   ABI_ALLOCATE(ekb,(psps%dimekb,psps%ntypat,dtset%nspinor**2))
   ekb(:,:,1)=psps%ekb(:,:)
   if (dtset%nspinor==2) then
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
   cplex_dij=max(1,dtset%nspinor)
   dimekb1=psps%dimekb*cplex_dij;dimekb2=natom
   ABI_ALLOCATE(ekb,(dimekb1,dimekb2,dtset%nspinor**2))
   ABI_ALLOCATE(sij,(dimekb1,psps%ntypat))
   ekb(:,:,:)=zero;sij(:,:)=zero
   do itypat=1,psps%ntypat
     if (cplex_dij==1) then
       sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
     else
       do klmn=1,pawtab(itypat)%lmn2_size
         sij(2*klmn-1,itypat)=pawtab(itypat)%sij(klmn)
       end do
     end if
   end do

!  If PAW and Born Eff. Charges, one has to compute some additional data:
!  For each atom and for electric field direction k:
!  becij(k)=<Phi_i|r_k-R_k|Phi_j>-<tPhi_i|r_k-R_k|tPhi_j> + sij.R_k
   if (need_becfr) then
     ABI_ALLOCATE(becij,(psps%dimekb,dimekb2,dtset%nspinor**2,3))
     becij=zero
     ABI_DATATYPE_ALLOCATE(paw_ij_tmp,(my_natom))
     ABI_DATATYPE_ALLOCATE(pawfgrtab_tmp,(my_natom))
     call paw_ij_nullify(paw_ij_tmp)
     call paw_ij_init(paw_ij_tmp,1,1,1,1,0,natom,psps%ntypat,dtset%typat,pawtab,has_dijfr=1,&
&     mpi_comm_atom=my_comm_atom,mpi_atmtab=my_atmtab )
     ABI_ALLOCATE(l_size_atm,(my_natom))
     do iatom=1,my_natom
       iatom_tot=iatom; if(paral_atom) iatom_tot=my_atmtab(iatom)
       itypat=dtset%typat(iatom_tot)
       l_size_atm(iatom)=pawtab(itypat)%lcut_size
     end do
     call pawfgrtab_init(pawfgrtab_tmp,1,l_size_atm,dtset%nspden,dtset%typat,&
&     mpi_atmtab=my_atmtab,mpi_comm_atom=my_comm_atom)
     ABI_DEALLOCATE(l_size_atm)
     do ii=1,3 ! Loop over direction of electric field
       call paw_ij_reset_flags(paw_ij_tmp,all=.True.)
       call pawdijfr(1,gprimd,ii,natom+2,my_natom,natom,nfftf,ngfftf,dtset%nspden,psps%ntypat,&
&       0,paw_ij_tmp,pawang,pawfgrtab_tmp,pawrad,pawtab,&
&       (/zero,zero,zero/),rprimd,ucvol,vtrial,vtrial,vxc,xred,&
&       mpi_comm_atom=my_comm_atom, mpi_atmtab=my_atmtab ) ! vtrial not used here
       do iatom=1,my_natom
         iatom_tot=iatom; if(paral_atom) iatom_tot=my_atmtab(iatom)
         itypat=dtset%typat(iatom_tot);dimdij=pawtab(itypat)%lmn2_size
!        Add contribution from Phi_i/Phi_j and tPhi_i/tPhi_j to becij
          do klmn=1,dimdij
           becij(klmn,iatom_tot,1,ii)=paw_ij_tmp(iatom)%dijfr(klmn,1)
         end do
       end do
     end do
     if (dtset%nspinor==2) becij(:,:,2,:)=becij(:,:,1,:)
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
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,dtset%nspinor))
   call pawcprj_alloc(cwaveprj,ncpgr,dimcprj)
   do iatom=1,natom
     sz2=pawrhoij_tot(iatom)%cplex*pawrhoij_tot(iatom)%lmn2_size
     sz3=pawrhoij_tot(iatom)%nspden
     ABI_ALLOCATE(pawrhoij_tot(iatom)%grhoij,(ngrhoij,sz2,sz3))
     pawrhoij_tot(iatom)%ngrhoij=ngrhoij
     pawrhoij_tot(iatom)%grhoij=zero
   end do
   usetimerev=(dtset%kptopt>0.and.dtset%kptopt<3)

 end if !PAW
 
!If needed, manage ddk files
 if (need_becfr) then
   do ii=1,3 ! Loop over elect. field directions
     if (ddkfil(ii)/=0) then
!      Open ddk WF file(s)
       write(msg, '(a,a)') '-open ddk wf file :',fiwfddk(ii)
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       call WffOpen(dtset%accesswff,spaceworld,fiwfddk(ii),ierr,wffddk(ii),master,me,ddkfil(ii))
!      Prepare DDK files for reading
       skipddk(:)=0;ikpt_fbz(1:3)=0;ikpt_fbz_previous(1:3)=0
       call clsopn(wffddk(ii))
       call hdr_skip(wffddk(ii),ierr)
     end if
   end do
 end if

!LOOP OVER SPINS
 do isppol=1,dtset%nsppol
   ikg=0

!  PAW: retrieve Dij coefficients
   if (psps%usepaw==1) then
     ekb(:,:,:)=zero
     if (my_natom>0) then
       do ispden=1,dtset%nspinor**2
         isp=isppol;if (dtset%nspinor==2) isp=ispden
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

!  If needed, manage ddk files
   if (need_becfr) then
     do ii=1,3 ! Loop over elect. field directions
       if (ddkfil(ii)/=0) then
!        ddk files: skip the remaining isppol=1 records
         if ((isppol==2).and.(skipddk(ii)<dtset%nkpt)) then
           do ikpt=1,(dtset%nkpt-skipddk(ii))
             call WffReadNpwRec(ierr,ikpt,isppol,nband_,npw_,nspinor_,wffddk(ii))
             call WffReadSkipRec(ierr,1+2*nband_,wffddk(ii))
           end do
         end if
       end if
     end do
   end if
   
!  Loop over k points
   do ikpt=1,dtset%nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     istwf_k=dtset%istwfk(ikpt)
     npw_k=npwarr(ikpt)
     wtk_k=dtset%wtk(ikpt)

!    Skip this k-point if not the proper processor
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

!    If needed, manage ddk files
     if (need_becfr) then
       do ii=1,3 ! Loop over elect. field directions
         if (ddkfil(ii)/=0)then
!        Skip records in ddk file
           ikpt_fbz_previous(ii)=ikpt_fbz(ii)
           ikpt_fbz(ii)=ikpt
!        Number of k points to skip in the full set of k pointsp
           nskip=ikpt_fbz(ii)-ikpt_fbz_previous(ii)-1
           skipddk(ii)=skipddk(ii)+nskip+1
           if (nskip/=0) then
             do ikpt_=1+ikpt_fbz_previous(ii),ikpt_fbz(ii)-1
               call WffReadSkipK(1,0,ikpt_,isppol,mpi_enreg,wffddk(ii))
             end do
           end if
!          Begin to read current record (k+G)
           call WffReadNpwRec(ierr,ikpt,isppol,nband_,npw_,nspinor_,wffddk(ii)) 
           if (npw_/=npw_k) then
            write(unit=msg,fmt='(a,i3,a,i5,a,i3,a,a,i5,a,a,i5)')&
&             'For isppol = ',isppol,', ikpt = ',ikpt,' and idir = ',ii,ch10,&
&             'the number of plane waves in the ddk file is equal to', npw_,ch10,&
&             'while it should be ',npw_k
               MSG_ERROR(msg)
             end if
             call WffReadSkipRec(ierr,1,wffddk(ii))
         end if
       end do
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     if(rfstrs/=0.or.need_becfr)then
       ABI_ALLOCATE(ylmgr_k,(npw_k,9,mpsang*mpsang*psps%useylm))
     else
       ABI_ALLOCATE(ylmgr_k,(0,0,0))
     end if

     kpoint(:)=dtset%kptns(:,ikpt)
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
       if(rfstrs/=0.or.need_becfr)then
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

!    Compute (k+G) vectors (only if useylm=1)     
     if (rfstrs/=0) nkpg=3*dtset%nloalg(5)
     if (rfphon==1) nkpg=9*dtset%nloalg(5)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

!    Compute nonlocal form factors ffnl at all (k+G):
     ider=0;dimffnl=1;
     dimffnl_bec=0;
     if(rfstrs/=0)then
       ider=2;dimffnl=3+7*psps%useylm
     end if

     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,&
&     psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

     if(need_becfr)then
       ider_bec=1; idir_bec=0; dimffnl_bec=4;
       ABI_ALLOCATE(ffnl_bec,(npw_k,dimffnl_bec,psps%lmnmax,psps%ntypat))
       call mkffnl(psps%dimekb,dimffnl_bec,psps%ekb,ffnl_bec,psps%ffspl,&
&       gmet,gprimd,ider_bec,idir_bec,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,&
&       psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
     end if

!    Compute phkxred and eventually ph3d.
     do iatom=1,natom
       ia=atindx1(iatom)
       arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
       phkxred(1,iatom)=cos(arg)
       phkxred(2,iatom)=sin(arg)
     end do
     if(dtset%nloalg(1)<=0)then
!      Only the allocation, not the precomputation.
       matblk=dtset%nloalg(4)
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
     if (need_becfr) becfrnlk(:,:,:)=zero

!    Loop over bands
     do iband=1,nband_k

       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,isppol,me)) cycle

       occ_k=occ(iband+bdtot_index)
       cwavef(:,1:npw_k*dtset%nspinor) = cg(:,1+(iband-1)*npw_k*dtset%nspinor+icg:iband*npw_k*dtset%nspinor+icg)

       if (psps%usepaw==1.and.usecprj==1) then
         ibsp=(iband-1)*dtset%nspinor+ibg
         call pawcprj_copy(cprj     (:,ibsp+1:ibsp+dtset%nspinor),cwaveprj)
       end if

!      Compute non-local contributions from n,k
       if (psps%usepaw==1) eig_k=eigen(iband+bdtot_index)

!      === Dynamical matrix
       if(rfphon==1) then
         call nonlop(atindx1,choice_phon,cpopt,cwaveprj,dimekb1,dimekb2,dimffnl,dimffnl,ekb,&
&         enlout_phon,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&         kpoint,(/eig_k/),psps%lmnmax,matblk,dtset%mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,nattyp,1,&
&         ngfft,nkpg,nkpg,dtset%nloalg,nnlout_phon,npw_k,npw_k,dtset%nspinor,dtset%nspinor,psps%ntypat,0,paw_opt,phkxred,&
&         phkxred,ph1d,ph3d,ph3d,signs,sij,nonlop_dum,4,ucvol,&
&         psps%useylm,cwavef,cwavef)
!        Accumulate non-local contributions from n,k
         dyfrnlk(:,:)=dyfrnlk(:,:)+occ_k*reshape(enlout_phon(:),(/6,natom/))
       end if

!      === Elastic tensor
       if(rfstrs/=0) then
         call nonlop(atindx1,choice_strs,cpopt,cwaveprj,dimekb1,dimekb2,dimffnl,dimffnl,ekb,&
&         enlout_strs,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&         kpoint,(/eig_k/),psps%lmnmax,matblk,dtset%mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,nattyp,1,&
&         ngfft,nkpg,nkpg,dtset%nloalg,nnlout_strs,npw_k,npw_k,dtset%nspinor,dtset%nspinor,psps%ntypat,0,paw_opt,phkxred,&
&         phkxred,ph1d,ph3d,ph3d,signs,sij,nonlop_dum,6,ucvol,&
&         psps%useylm,cwavef,cwavef)
!        Accumulate non-local contribut ions from n,k
         eltfrnlk(:,:)=eltfrnlk(:,:)+occ_k*reshape(enlout_strs(:),(/3*natom+6,6/))
       end if !endxo if strs

!      PAW: accumulate gradients of rhoij
       if (psps%usepaw==1) then 
         call pawaccrhoij(atindx,cplex,cwaveprj,cwaveprj,0,isppol,natom,&
&         natom,dtset%nspinor,occ_k,3,pawrhoij_tot,usetimerev,wtk_k)
       end if
       
!      PAW: Compute frozen contribution to Born Effective Charges
       if (need_becfr) then
         do ii=1,3 ! Loop over elect. field directions
           call nonlop(atindx1,choice_bec2,cpopt_bec,cwaveprj,psps%dimekb,dimekb2,dimffnl_bec,dimffnl_bec,&
&            becij(:,:,:,ii),enlout_bec1,ffnl_bec,ffnl_bec,gmet,gprimd,0,psps%indlmn,istwf_k,kg_k,kg_k,&
&            kpg_k,kpg_k,kpoint,kpoint,(/zero/),psps%lmnmax,matblk,dtset%mgfft,mpi_enreg,mpsang,psps%mpssoang,&
&            natom,nattyp,1,ngfft,nkpg,nkpg,dtset%nloalg,nnlout_bec1,npw_k,npw_k,dtset%nspinor,dtset%nspinor,psps%ntypat,&
&            0,paw_opt_bec1,phkxred,phkxred,ph1d,ph3d,ph3d,signs_bec1,sij,nonlop_dum,tim_nonlop,ucvol,&
&            psps%useylm,cwavef,cwavef)
           becfrnlk(:,:,ii)=becfrnlk(:,:,ii)+occ_k*reshape(enlout_bec1(:),(/3,natom/))
         end do !end do ii

         ABI_ALLOCATE(svectout,(2,npw_k*dtset%nspinor))
         do ii=1,3 ! Loop over elect. field directions
!          Not able to compute if ipert=(Elect. field) and no ddk WF file
           if (ddkfil(ii)==0) cycle
!            Read ddk wave function 
           ABI_ALLOCATE(ddk,(2,npw_k*dtset%nspinor))
           if (ddkfil(ii)/=0) then
             call WffReadSkipRec(ierr,1,wffddk(ii))
             call WffReadDataRec(ddk,ierr,2,npw_k*dtset%nspinor,wffddk(ii))
!            Multiply ddk by +i
             do jj=1,npw_k*dtset%nspinor
               arg=ddk(1,jj)
               ddk(1,jj)=-ddk(2,jj);ddk(2,jj)=arg
             end do
           else
             ddk=zero
           end if

           do iatom=1,natom !Loop over atom
             natom_der=1 ; nattyp_der(1)=1 ; ntypat_der=1
             dimekb2_der=1
             matblk_der=1
             xred_der(:)=xred(:,iatom)
             atindx1_der(1)=1
             ia=atindx(iatom)
             ABI_ALLOCATE(indlmn_der,(6,psps%lmnmax,1))
             indlmn_der(:,:,1)=psps%indlmn(:,:,dtset%typat(iatom))
!            Store at the right place the 1d phases
             ABI_ALLOCATE(ph1d_der,(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
             shift1=(ia-1)*(2*n1+1)
             ph1d_der(:,1:2*n1+1)=ph1d(:,1+shift1:2*n1+1+shift1)
             shift2=(ia-1)*(2*n2+1)+natom*(2*n1+1)
             ph1d_der(:,1+2*n1+1:2*n2+1+2*n1+1)=ph1d(:,1+shift2:2*n2+1+shift2)
             shift3=(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
             ph1d_der(:,1+2*n1+1+2*n2+1:2*n3+1+2*n2+1+2*n1+1)=ph1d(:,1+shift3:2*n3+1+shift3)
!            Will compute the 3D phase factors inside nonlop
             ABI_ALLOCATE(ph3din,(2,npw_k,1))
             ABI_ALLOCATE(ph3dout,(2,npw_k,1))
             nloalg_der(:)=dtset%nloalg(:)
             nloalg_der(1)=-abs(dtset%nloalg(1))
             nloalg_der(4)=1
!            Retrieve here phkxred for kpt and kpq
             phkxredin(:,1) =phkxred(:,ia)
             phkxredout(:,1)=phkxred(:,ia)
!            Retrieve Sij for this atom type
             ABI_ALLOCATE(ekb_der,(dimekb1,1,dtset%nspinor**2))
             ABI_ALLOCATE(sij_der,(dimekb1,1))
             sij_der(:,1) = sij(:,dtset%typat(iatom))
            
!            PAW: retrieve ffnlk and cwaveprj for the displaced atom
             if (psps%usepaw==1) then
               ABI_ALLOCATE(ffnlk_der,(npw_k,dimffnl_bec,psps%lmnmax,1))
               ffnlk_der(:,:,:,1)=ffnl_bec(:,:,:,dtset%typat(iatom))
               if (usecprj==1) then
                 dimcprj_der(1)=dimcprj(ia)
                 ABI_DATATYPE_ALLOCATE(cprj_der,(1,dtset%nspinor))
                 call pawcprj_alloc(cprj_der,1,dimcprj_der)
                 do ispinor=1,dtset%nspinor
                   call pawcprj_copy(cwaveprj(ia:ia,ispinor:ispinor),cprj_der(1:1,ispinor:ispinor))
                 end do
               end if
             end if
             do mu=1,3 !loop over atom direction 
               svectout = zero
               call nonlop(atindx1_der,choice_bec2,cpopt_bec,cprj_der,psps%dimekb,dimekb2_der,dimffnl_bec,dimffnl_bec,&
&                ekb_der,enlout_bec1,ffnlk_der,ffnlk_der,gmet,gprimd,mu,indlmn_der,istwf_k,kg_k,kg_k,&
&                kpg_k,kpg_k,kpoint,kpoint,(/zero/),psps%lmnmax,matblk_der,dtset%mgfft,mpi_enreg,mpsang,psps%mpssoang,&
&                natom_der,nattyp_der,1,ngfft,nkpg,nkpg,nloalg_der,nnlout_bec1,npw_k,npw_k,dtset%nspinor,dtset%nspinor,ntypat_der,&
&                0,paw_opt_bec3,phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,signs_bec2,sij_der,svectout,tim_nonlop,ucvol,&
&                psps%useylm,cwavef,cwavef)
               
               call dotprod_g(dotprod(1),dotprod(2),istwf_k,npw_k*dtset%nspinor,2,&
&                             svectout,ddk,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
               becfrnlk(mu,iatom,ii)=becfrnlk(mu,iatom,ii)+occ_k*dotprod(1)
             end do ! End loop atom direction

             ABI_DEALLOCATE(ph1d_der)
             ABI_DEALLOCATE(ph3din)
             ABI_DEALLOCATE(ph3dout)
             ABI_DEALLOCATE(sij_der)
             ABI_DEALLOCATE(ekb_der)
             ABI_DEALLOCATE(indlmn_der)
             if (psps%usepaw==1) then
               ABI_DEALLOCATE(ffnlk_der)
               if (usecprj==1) then
                 call pawcprj_destroy(cprj_der)
                 ABI_DATATYPE_DEALLOCATE(cprj_der)
               end if
             end if
           end do ! End loop atom
           ABI_DEALLOCATE(ddk)
         end do ! End loop ddk file
         ABI_DEALLOCATE(svectout)
           
         call nonlop(atindx1,choice_bec54,cpopt_bec,cwaveprj,psps%dimekb,dimekb2,dimffnl_bec,dimffnl_bec,&
&          ekb,enlout_bec2,ffnl_bec,ffnl_bec,gmet,gprimd,0,psps%indlmn,istwf_k,kg_k,kg_k,&
&          kpg_k,kpg_k,kpoint,kpoint,(/zero/),psps%lmnmax,matblk,dtset%mgfft,mpi_enreg,mpsang,psps%mpssoang,&
&         natom,nattyp,1,ngfft,nkpg,nkpg,dtset%nloalg,nnlout_bec2,npw_k,npw_k,dtset%nspinor,dtset%nspinor,psps%ntypat,&
&         0,paw_opt_bec3,phkxred,phkxred,ph1d,ph3d,ph3d,signs_bec1,sij,nonlop_dum,tim_nonlop,ucvol,&
&         psps%useylm,cwavef,cwavef)
         
!        Multiply enlout by +i
         iashift = 1
         do iatom=1,natom ! atm
           do mu=1,3 ! atm pos.
             do nu=1,3 ! k
               becfrnlk(mu,iatom,nu)=becfrnlk(mu,iatom,nu)-occ_k*(enlout_bec2(iashift+1))! Real part
!              becfrnlk(mu,iatom,nu)=becfrnlk(mu,iatom,nu)+occ_k*(enlout_bec2(iashift  )) ! Imaginary part 
               iashift = iashift + 2
             end do
           end do
         end do
       end if
     end do! End of loop on bands

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
     if(rfstrs/=0) then
       eltfrnl(:,:)=eltfrnl(:,:)+dtset%wtk(ikpt)*eltfrnlk(:,:)
     end if
    if(need_becfr)then
       becfrnl(:,:,:)=becfrnl(:,:,:)+dtset%wtk(ikpt)*becfrnlk(:,:,:)
     end if

!    Increment indexes
     bdtot_index=bdtot_index+nband_k
     if (dtset%mkmem/=0) then
       ibg=ibg+nband_k*dtset%nspinor
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_DEALLOCATE(ffnl)
     if (need_becfr) then
       ABI_DEALLOCATE(ffnl_bec)
     end if
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)

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
 if (need_becfr)  then
   ABI_DEALLOCATE(becfrnlk)
   ABI_DEALLOCATE(enlout_bec1)
   ABI_DEALLOCATE(enlout_bec2)
 end if
 if (psps%usepaw==1) then
   if (need_becfr)  then
     ABI_DEALLOCATE(becij)
   end if
   call pawcprj_destroy(cwaveprj)
 end if
 ABI_DEALLOCATE(sij)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

!Fill in lower triangle of matrixes
 if (rfphon==1) then
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
   if (need_becfr) then
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
 if (rfstrs/=0)then
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

!====== PAW: Additional steps
 if (psps%usepaw==1) then

!  Symmetrize rhoij gradients and transfer to cartesian (reciprocal space) coord.
!  This symetrization is necessary in the antiferromagnetic case...
   if (rfphon==1.and.rfstrs==0) then
     option_rhoij=2;option=0
     call symrhoij(pawrhoij_tot,pawrhoij_tot,option_rhoij,gprimd,indsym,0,natom,dtset%nsym,&
&     psps%ntypat,option,pawang,dtset%pawprtvol,pawtab,rprimd,dtset%symafm,symrec,dtset%typat,&
&     mpi_comm_atom=my_comm_atom, mpi_atmtab=my_atmtab)
   else if (rfphon==1.and.rfstrs==1) then
     option_rhoij=23;option=0
     call symrhoij(pawrhoij_tot,pawrhoij_tot,option_rhoij,gprimd,indsym,0,natom,dtset%nsym,&
&     psps%ntypat,option,pawang,dtset%pawprtvol,pawtab,rprimd,dtset%symafm,symrec,dtset%typat,&
&     mpi_comm_atom=my_comm_atom, mpi_atmtab=my_atmtab)
   end if

!  Translate coordinates
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
&   nattyp,nfftf,ngfftf,nhat_dum,dummy,dtset%nspden,dtset%nsym,psps%ntypat,optgr,optgr2,optstr,optstr2,&
&   pawang,pawfgrtab,pawrhoij_tot,pawtab,ph1df,psps,dtset%qptn,rprimd,symrec,dtset%typat,vtrial,vxc,xred,&
&   mpi_atmtab=my_atmtab,mpi_comm_atom=my_comm_atom)
   ABI_DEALLOCATE(nhat_dum)
 end if !PAW

!Born Effective Charges and PAW:
!1-Re-order atoms -- 2-Add diagonal contribution from rhoij
!3-Multiply by -1 because that the effective charges
!  are minus the second derivatives of the energy
 if (need_becfr) then
   ABI_ALLOCATE(becfrnl_tmp,(3,natom,3))
   becfrnl_tmp=-becfrnl
   do ia=1,natom         ! Atom (sorted by type)
     iatom=atindx1(ia)   ! Atom (not sorted)
     itypat=dtset%typat(iatom)
     do ii=1,3           ! Direction of electric field
       do jj=1,3         ! Direction of atom
         becfrnl(jj,iatom,ii)=becfrnl_tmp(jj,ia,ii)
       end do
     end do
   end do
   ABI_DEALLOCATE(becfrnl_tmp)
 end if
 
!Close the ddk files
 if (need_becfr) then
   do ii=1,3
     if (ddkfil(ii)/=0)then
       call WffClose(wffddk(ii),ierr)
     end if
   end do
 end if

!PAW: release now useless memory
 if (psps%usepaw==1) then
   do iatom=1,natom
     ABI_DEALLOCATE(pawrhoij_tot(iatom)%grhoij)
     pawrhoij_tot(iatom)%ngrhoij=0
   end do
   if (paral_atom) then
     call pawrhoij_destroy(pawrhoij_tot)
     ABI_DATATYPE_DEALLOCATE(pawrhoij_tot)
   end if
 end if

 call timab(159,2,tsec)

 DBG_EXIT("COLL")

end subroutine d2frnl_bec
!!***
