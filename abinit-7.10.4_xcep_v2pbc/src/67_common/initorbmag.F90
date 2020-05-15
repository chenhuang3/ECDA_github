!{\src2tex{textfont=tt}}
!!****f* ABINIT/initorbmag
!! NAME
!! initorbmag
!!
!! FUNCTION
!! Initialization of orbital magnetism calculation,  
!! and the response of an insulator to a homogenous magnetic field.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2014 ABINIT group (JWZ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3) = primitive translations in recip space
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  mband = maximum number of bands
!!  mkmem = maximum number of k-points in core memory
!!  mpw = maximum number of plane waves
!!  natom = number of atoms in unit cell
!!  nkpt = number of k points
!!  npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  nsym = number of symmetry operations
!!  ntypat = number of types of atoms in unit cell
!!  occ(mband*nkpt*nsppol) = occup number for each band at each k point
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive vectors
!!  symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!    reciprocal space primitive translations
!!  typat = typat(natom) list of atom types
!!  usepaw = flag for PAW (1 PAW, 0 NCPP)
!!  xred(3,natom) = location of atoms in reduced units
!!
!! OUTPUT
!!  dtbfield <type(bfield_type)> = variables related to orbital magnetization calculations
!!  pwind(pwind_alloc,2,3) = array used to compute the overlap matrix smat
!!                         between k-points k and k +- dk where dk is
!!                         parallel to the direction idir
!!    jpw = pwind(ipw,ifor,idir)
!!      * ipw = index of plane wave vector G for a given k-point k
!!      * ifor = 1: k + dk
!!               2: k - dk
!!      * idir = direction of the polarization/ddk calculation [dk(idir)
!!               is the only non-zero element of dk(:)]
!!      * jpw = index of plane wave vector G (+dG) at k +- dk
!!              where dG is a shift of one reciprocal lattice vector
!!              (required to close the strings of k-points using the
!!               periodic gauge condition)
!!    In case a G-vector of the basis sphere of plane waves at k
!!    does not belong to the basis sphere of plane waves at k+dk, jpw = 0.
!!   pwind_alloc = first dimension of pwind and pwnsfac
!!   pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!
!! SIDE EFFECTS
!!  mpi_enreg = information about MPI parallelization
!!    kptdstrb(nproc,nneighbour,fmkmem_max*nsppol) : Array required
!!      for MPI // over k-points. Defined
!!      for k-points in the fBZ
!!      but for k-points in the iBZ. Used by vtorho.f
!!           nproc = number of cpus
!!           nneighbour = number of neighbours for each k-point (= 6)
!!
!! TO DO
!!
!! NOTES
!! This code is derived from initberry.F90
!!
!! PARENTS
!!      init_eb_field_vars
!!
!! CHILDREN
!!      expibi,kpgsph,lij,listkk,pawcprj_alloc,pawcprj_getdim,pawtwdij_1
!!      pawtwdij_2a,pawtwdij_2b,pawtwdij_2c,pawtwdij_2d,pawtwdij_2e,pawtwdij_2f
!!      qijb_kk,set_twind,setsymrhoij,smpbz,symatm,timab,twexpibi,twqijb_kk
!!      wrtout,xmpi_max,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initorbmag(dtbfield,dtset,gmet,gprimd,kg,mband,&
&              mkmem,mpi_enreg,mpw,natom,nkpt,npwarr,nsppol,&
&              nsym,ntypat,occ,pawang,pawrad,pawtab,psps,&
&              pwind,pwind_alloc,pwnsfac,&
&              rprimd,symrec,typat,usepaw,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_bfield

 use m_fftcore, only : kpgsph
 use m_pawang,  only : pawang_type
 use m_pawrad,  only : pawrad_type
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_alloc, pawcprj_getdim

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initorbmag'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_66_paw
 use interfaces_67_common, except_this_one => initorbmag
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,natom,nkpt,nsppol,nsym,ntypat,usepaw
 integer,intent(out) :: pwind_alloc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
 type(bfield_type),intent(inout) :: dtbfield
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 integer,intent(in) :: symrec(3,3,nsym),typat(natom)
 integer,pointer :: pwind(:,:,:)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),pointer :: pwnsfac(:,:)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: exchn2n3d,flag,flag_kpt,fnkpt_computed,iband,icg,icprj
 integer :: idir,idum,idum1,ierr,ifor,ikg,ikg1,ikpt,ikpt1,ikpt1f
 integer :: ikpt1i,ikpt_loc,ikptf,ikpti,index,ineigh,ipw,ipwnsfac
 integer :: isppol,istwf_k,isym,isym1,itrs,itypat,jpw,lmax,lmn2_size_max
 integer :: me,me_g0,mkmem_,my_nspinor,nband_k,nband_occ_k,ncpgr,nproc,npw_k,npw_k1,spaceComm
 integer :: option, brav, mkpt, nkptlatt, isign
 real(dp) :: ecut_eff,rdum,diffk1,diffk2,diffk3
 real(dp) :: kpt_shifted1,kpt_shifted2,kpt_shifted3

 character(len=500) :: message
 logical :: use_symrec
!arrays
 integer :: dg(3),iadum(3),iadum1(3),neigh(6)
 integer,allocatable :: dimlmn(:),kg1_k(:,:),nattyp_dum(:)
 real(dp) :: diffk(3),dk(3),dum33(3,3)
 real(dp) :: kpt1(3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: spkpt(:,:)

! *************************************************************************

 call timab(1001,1,tsec)
 call timab(1002,1,tsec)

!Compatibility tests
 if (usepaw /= 1) then
   message = ' usepaw /= 1 but orbital magnetization requires PAW '
   MSG_ERROR(message)
 end if

!save the current value of berryopt
 dtbfield%orbmag = dtset%orbmag
!save the current value of nspinor
 dtbfield%nspinor = dtset%nspinor
 dtbfield%ndij = (dtbfield%nspinor)**2

!----------------------------------------------------------------------------
!-------------------- Obtain k-point grid in the full BZ --------------------
!----------------------------------------------------------------------------

 if(dtset%kptopt==1 .or. dtset%kptopt==2 .or. dtset%kptopt==4)then
!  Compute the number of k points in the G-space unit cell
   nkptlatt=dtset%kptrlatt(1,1)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,3) &
&   +dtset%kptrlatt(1,2)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,1) &
&   +dtset%kptrlatt(1,3)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,2) &
&   -dtset%kptrlatt(1,2)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,3) &
&   -dtset%kptrlatt(1,3)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,1) &
&   -dtset%kptrlatt(1,1)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,2)

!  Call smpbz to obtain the list of k-point in the full BZ - without symmetry reduction
   option = 0
   brav = 1
   mkpt=nkptlatt*dtset%nshiftk
   ABI_ALLOCATE(spkpt,(3,mkpt))
   call smpbz(1,ab_out,dtset%kptrlatt,mkpt,fnkpt_computed,dtset%nshiftk,option,dtset%shiftk,spkpt)
   dtbfield%fnkpt = fnkpt_computed
   ABI_ALLOCATE(dtbfield%fkptns,(3,dtbfield%fnkpt))
   dtbfield%fkptns(:,:)=spkpt(:,1:dtbfield%fnkpt)
   ABI_DEALLOCATE(spkpt)
 else if(dtset%kptopt==3.or.dtset%kptopt==0)then
   dtbfield%fnkpt=nkpt
   ABI_ALLOCATE(dtbfield%fkptns,(3,dtbfield%fnkpt))
   dtbfield%fkptns(1:3,1:dtbfield%fnkpt)=dtset%kpt(1:3,1:dtbfield%fnkpt)
   if(dtset%kptopt==0)then
     write(message,'(10a)') ch10,&
&     ' initorbmag : WARNING -',ch10,&
&     '  you have defined manually the k-point grid with kptopt = 0',ch10,&
&     '  the orbital magnetization calculation works only with a regular k-points grid,',ch10,&
&     '  abinit doesn''t check if your grid is regular...'
     call wrtout(std_out,message,'PERS')
   end if
 end if

!call listkk to get mapping from FBZ to IBZ
 rdum=1.0d-5  ! cutoff distance to decide when two k points match
 ABI_ALLOCATE(dtbfield%indkk_f2ibz,(dtbfield%fnkpt,6))

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!ji: The following may need modification in the future
!**** no spin-polarization doubling ; allow use of time reversal symmetry ****

!Here is original call
!
!call listkk(rdum,gmet,dtefield%indkk_f2ibz,dtset%kptns,dtefield%fkptns,nkpt,&
!& dtefield%fnkpt,dtset%nsym,1,dtset%symafm,dtset%symrel,1)

 call timab(1002,2,tsec)
 call timab(1003,1,tsec)

 use_symrec = .TRUE.
 call listkk(rdum,gmet,dtbfield%indkk_f2ibz,dtset%kptns,dtbfield%fkptns,nkpt,&
& dtbfield%fnkpt,dtset%nsym,1,dtset%symafm,symrec,1,use_symrec)

 call timab(1003,2,tsec)
 call timab(1004,1,tsec)

!Construct i2fbz and f2ibz
 ABI_ALLOCATE(dtbfield%i2fbz,(nkpt))
 idum=0
 do ikpt=1,dtbfield%fnkpt
   if (dtbfield%indkk_f2ibz(ikpt,2)==1 .and. &
&   dtbfield%indkk_f2ibz(ikpt,6) == 0 .and. &
&   maxval(abs(dtbfield%indkk_f2ibz(ikpt,3:5))) == 0 ) then
     dtbfield%i2fbz(dtbfield%indkk_f2ibz(ikpt,1))=ikpt
     idum=idum+1
   end if
 end do
 if (idum/=nkpt)then
   message = ' Found wrong number of k-points in IBZ'
   MSG_ERROR(message)
 end if

!----------------------------------------------------------------------------
!------------- Allocate PAW space -------------------------------------------
!----------------------------------------------------------------------------

 dtbfield%usepaw   = usepaw
 dtbfield%natom    = natom
 dtbfield%my_natom = mpi_enreg%my_natom
 lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
 dtbfield%lmn2max = lmn2_size_max

 ABI_ALLOCATE(dtbfield%lmn_size,(ntypat))
 ABI_ALLOCATE(dtbfield%lmn2_size,(ntypat))
 do itypat = 1, ntypat
   dtbfield%lmn_size(itypat) = pawtab(itypat)%lmn_size
   dtbfield%lmn2_size(itypat) = pawtab(itypat)%lmn2_size
 end do

 if( dtbfield%has_qijb == 0) then 
   ABI_ALLOCATE(dtbfield%qijb_kk,(2,lmn2_size_max,natom,3))
   ABI_ALLOCATE(dtbfield%twqijb_kk,(2,lmn2_size_max,natom,6))
   dtbfield%has_qijb = 1
 end if

 if( dtbfield%has_expibi == 0) then 
   ABI_ALLOCATE(dtbfield%expibi,(2,dtbfield%my_natom,3))
   ABI_ALLOCATE(dtbfield%twexpibi,(2,dtbfield%my_natom,6))
   dtbfield%has_expibi = 1
 end if

 if(dtbfield%has_Lij==0) then
   ABI_ALLOCATE(dtbfield%Lij,(2,lmn2_size_max,ntypat,3))
   dtbfield%has_Lij = 1
 end if

 if(dtbfield%has_Lijr3==0) then
   ABI_ALLOCATE(dtbfield%Lijr3,(2,lmn2_size_max,ntypat,3))
   dtbfield%has_Lijr3 = 1
 end if

 if (dtbfield%usecprj == 0) then
   ABI_ALLOCATE(dimlmn,(natom))
   call pawcprj_getdim(dimlmn,natom,nattyp_dum,ntypat,typat,pawtab,'R')
!  allocate space for cprj at kpts in BZ (IBZ or FBZ)
   ABI_DATATYPE_ALLOCATE(dtbfield%cprj,(natom, mband*dtset%nspinor*dtset%nkpt*nsppol))
   ncpgr = 3 ! gradients for forces (forces not yet coded though)
   call pawcprj_alloc(dtbfield%cprj,ncpgr,dimlmn)
   dtbfield%usecprj = 1
   ABI_DEALLOCATE(dimlmn)
 end if

 ABI_ALLOCATE(dtbfield%cprjindex,(nkpt,nsppol))
 dtbfield%cprjindex(:,:) = 0

 if (dtset%kptopt /= 3) then
   ABI_ALLOCATE(dtbfield%atom_indsym,(4,nsym,natom))
   call symatm(dtbfield%atom_indsym,natom,nsym,symrec,dtset%tnons,tol8,typat,xred)
   lmax = psps%mpsang - 1
   ABI_ALLOCATE(dtbfield%zarot,(2*lmax+1,2*lmax+1,lmax+1,nsym))
   call setsymrhoij(gprimd,lmax,nsym,1,rprimd,symrec,dtbfield%zarot)
   dtbfield%nsym = nsym
   dtbfield%lmax = lmax
   dtbfield%lmnmax = psps%lmnmax
 end if

 ABI_ALLOCATE(dtbfield%mag_local_k,(3,nkpt*nsppol))
 ABI_ALLOCATE(dtbfield%mag_k,(2,nkpt*nsppol,3))

 ABI_ALLOCATE(dtbfield%twdij0,(2,psps%lmnmax,psps%lmnmax,dtbfield%my_natom,24))
 dtbfield%twdij0(:,:,:,:,:) = zero
 dtbfield%has_twdij0 = 1

 ABI_ALLOCATE(dtbfield%tweijkl,(2,lmn2_size_max,lmn2_size_max,dtbfield%my_natom,6))
 dtbfield%tweijkl(:,:,:,:,:) = zero
 dtbfield%has_tweijkl = 1

 ABI_ALLOCATE(dtbfield%twdij,(2,psps%lmnmax,psps%lmnmax,natom,24,dtbfield%ndij))
 dtbfield%twdij(:,:,:,:,:,:) = zero

!------------------------------------------------------------------------------
!------------------- Compute variables related to MPI // ----------------------
!------------------------------------------------------------------------------

 spaceComm=mpi_enreg%comm_cell
 nproc=xcomm_size(spaceComm)
 me=xcomm_rank(spaceComm)

 if (xmpi_paral== 0) then  ! no MPI //

   dtbfield%fmkmem = dtbfield%fnkpt
   dtbfield%fmkmem_max = dtbfield%fnkpt
   dtbfield%mkmem_max = nkpt

!  we allocate also those
   ABI_ALLOCATE(mpi_enreg%kpt_loc2fbz_sp,(0:nproc-1,1:dtbfield%fmkmem_max*nsppol, 1:2))
   ABI_ALLOCATE(mpi_enreg%mkmem,(0:nproc-1))
   mpi_enreg%kpt_loc2fbz_sp(:,:,:) = 0
   mpi_enreg%mkmem(:) = 0

   ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:nproc-1,1:dtbfield%mkmem_max*nsppol, 1:2))
   mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0

 else    ! MPI //

!  Number of k-points in the FBZ for each cpu

   dtbfield%fmkmem = 0
   do ikpt = 1, dtbfield%fnkpt
     ikpti = dtbfield%indkk_f2ibz(ikpt,1)
     nband_k = dtset%nband(ikpti)
     if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,-1,me))) &
&     dtbfield%fmkmem = dtbfield%fmkmem + 1
   end do

!  Maximum value of mkmem and fmkmem
   call xmpi_max(dtbfield%fmkmem,dtbfield%fmkmem_max,spaceComm,ierr)

   mkmem_ = mkmem   ! I have to use the dummy variable mkmem_ because
!  mkmem is declared as intent(in) while the first
!  argument of xmpi_max must be intent(inout)

   call xmpi_max(mkmem_,dtbfield%mkmem_max,spaceComm,ierr)

   ABI_ALLOCATE(mpi_enreg%kptdstrb,(nproc,6,dtbfield%fmkmem_max*nsppol*2))
   mpi_enreg%kptdstrb(:,:,:) = 0

   ABI_ALLOCATE(mpi_enreg%kpt_loc2fbz_sp,(0:nproc-1,1:dtbfield%fmkmem_max*nsppol, 1:2))
   ABI_ALLOCATE(mpi_enreg%mkmem,(0:nproc-1))
   mpi_enreg%kpt_loc2fbz_sp(:,:,:) = 0
   mpi_enreg%mkmem(:) = 0

   ABI_ALLOCATE(dtbfield%cgqindex,(3,6,nkpt*nsppol))
   ABI_ALLOCATE(dtbfield%nneigh,(nkpt))
   dtbfield%cgqindex(:,:,:) = 0 ; dtbfield%nneigh(:) = 0

   ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:nproc-1,1:dtbfield%mkmem_max*nsppol, 1:2))
   mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0

 end if

 pwind_alloc = mpw*dtbfield%fmkmem_max
 ABI_ALLOCATE(pwind,(pwind_alloc,2,3))
 ABI_ALLOCATE(pwnsfac,(2,pwind_alloc))

!------------------------------------------------------------------------------
!---------------------- Compute bfield_type variables -------------------------
!------------------------------------------------------------------------------

!Initialization of bfield_type variables
 dtbfield%dkvecs(:,:) = zero
 ABI_ALLOCATE(dtbfield%ikpt_dk,(dtbfield%fnkpt,2,3))
 ABI_ALLOCATE(dtbfield%cgindex,(nkpt,nsppol))
 ABI_ALLOCATE(dtbfield%kgindex,(nkpt))
 ABI_ALLOCATE(dtbfield%fkgindex,(dtbfield%fnkpt))
 dtbfield%ikpt_dk(:,:,:) = 0
 dtbfield%cgindex(:,:) = 0
 dtbfield%nband_occ = 0
 dtbfield%kgindex(:) = 0
 dtbfield%fkgindex(:) = 0

 dtset%rfdir(1:3) = 1

!Compute spin degeneracy
 if (nsppol == 1 .and. dtset%nspinor == 1) then
   dtbfield%sdeg = two
 else if (nsppol == 2 .or. my_nspinor == 2) then
   dtbfield%sdeg = one
 end if

!Compute the number of occupied bands and check that
!it is the same for each k-point

 index = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt

     nband_occ_k = 0
     nband_k = dtset%nband(ikpt + (isppol - 1)*nkpt)

     do iband = 1, nband_k
       index = index + 1
       if (abs(occ(index) - dtbfield%sdeg) < tol8) nband_occ_k = nband_occ_k + 1
     end do

     if (nband_k /= nband_occ_k) then
       write(message,'(a,a,a)')&
&       '  For magnetization calculation,  nband must be equal ',ch10,&
&       '  to the number of valence bands.'
       MSG_ERROR(message)
     end if

     if ((ikpt > 1).or.(isppol > 1)) then
       if (dtbfield%nband_occ /= nband_occ_k) then
         message = "The number of valence bands is not the same for every k-point"
         MSG_ERROR(message)
       end if
     else
       dtbfield%nband_occ = nband_occ_k
     end if

   end do                ! close loop over ikpt
 end do                ! close loop over isppol

 ABI_ALLOCATE(dtbfield%smat,(2,dtbfield%nband_occ,dtbfield%nband_occ,nkpt*nsppol,2,3))
 dtbfield%smat(:,:,:,:,:,:) = zero

 ABI_ALLOCATE(dtbfield%chern_k,(2,nkpt*nsppol,3))
 dtbfield%chern_k(:,:,:) = zero
 
 ABI_ALLOCATE(dtbfield%emat,(2,dtbfield%nband_occ,nkpt*nsppol))
 dtbfield%emat(:,:,:) = zero

!if (abs(dtset%berryopt)==5) then
!ABI_ALLOCATE(dtefield%twh,(2,dtefield%nband_occ,dtefield%nband_occ,nkpt*nsppol,24))
!dtefield%twh(:,:,:,:,:) = zero
!ABI_ALLOCATE(dtefield%chern_k,(2,nkpt*nsppol,3))
!dtefield%chern_k(:,:,:) = zero
!
!ABI_ALLOCATE(dtefield%emat,(2,dtefield%nband_occ,nkpt*nsppol))
!
!end if 

!if (dtset%berryopt == 5) then
!ihcg = 0
!do ikpt = 1, nkpt
!ihcg = ihcg + npwarr(ikpt)*dtefield%nband_occ
!end do
!dtefield%mhcg = ihcg
!ABI_ALLOCATE(dtefield%hcg,(2,dtefield%mhcg,24))
!dtefield%hcg(:,:,:) = zero
!
!I think that dtefield%cgindex (defined below) will also give location
!in hcg

!end if

 ABI_ALLOCATE(dtbfield%sflag,(dtbfield%nband_occ,nkpt*nsppol,2,3))
 dtbfield%sflag(:,:,:,:) = 0

!Compute the location of each wavefunction

 icg = 0
 icprj = 0
!ikg = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt

     nband_k = dtset%nband(ikpt + (isppol-1)*nkpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle
     
     dtbfield%cgindex(ikpt,isppol) = icg
     npw_k = npwarr(ikpt)
     icg = icg + npw_k*dtbfield%nspinor*nband_k

     dtbfield%cprjindex(ikpt,isppol) = icprj
     icprj = icprj + dtbfield%nspinor*nband_k

   end do
 end do

 ikg = 0
 do ikpt = 1, nkpt
   if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,1,me)).and.&
&   (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,nsppol,me))) cycle
   
   npw_k = npwarr(ikpt)
   dtbfield%kgindex(ikpt) = ikg
   ikg = ikg + npw_k
 end do

!Store magnetic field
!if (dtset%berryopt == 5) dtefield%bfield(:)=dtset%bfield(:)

 call timab(1004,2,tsec)

!------------------------------------------------------------------------------
!---------------------- Compute dk --------------------------------------------
!------------------------------------------------------------------------------

 call timab(1005,1,tsec)

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

!    Compute dk(:), the vector between a k-point and its nearest
!    neighbour along the direction idir

     dk(:) = zero
     dk(idir) = 1._dp   ! 1 mean there is no other k-point un the direction idir
     do ikpt = 2, dtbfield%fnkpt
       diffk(:) = abs(dtbfield%fkptns(:,ikpt) - dtbfield%fkptns(:,1))
       if ((diffk(1) < dk(1)+tol8).and.(diffk(2) < dk(2)+tol8).and.&
&       (diffk(3) < dk(3)+tol8)) dk(:) = diffk(:)
     end do
     dtbfield%dkvecs(:,idir) = dk(:)
!    DEBUG
!    write(std_out,*)' initberry : idir, dk', idir, dk
!    ENDDEBUG

!    For each k point, find k_prim such that k_prim= k + dk mod(G)
!    where G is a vector of the reciprocal lattice
     do ikpt = 1, dtbfield%fnkpt

!      First k+dk, then k-dk
       do isign=-1,1,2
         kpt_shifted1=dtbfield%fkptns(1,ikpt)- isign*dk(1)
         kpt_shifted2=dtbfield%fkptns(2,ikpt)- isign*dk(2)
         kpt_shifted3=dtbfield%fkptns(3,ikpt)- isign*dk(3)
!        Note that this is still a order fnkpt**2 algorithm.
!        It is possible to implement a order fnkpt algorithm, see listkk.F90.
         do ikpt1 = 1, dtbfield%fnkpt
           diffk1=dtbfield%fkptns(1,ikpt1) - kpt_shifted1
           if(abs(diffk1-nint(diffk1))>tol8)cycle
           diffk2=dtbfield%fkptns(2,ikpt1) - kpt_shifted2
           if(abs(diffk2-nint(diffk2))>tol8)cycle
           diffk3=dtbfield%fkptns(3,ikpt1) - kpt_shifted3
           if(abs(diffk3-nint(diffk3))>tol8)cycle
           dtbfield%ikpt_dk(ikpt,(isign+3)/2,idir) = ikpt1
           exit
         end do   ! ikpt1
       end do     ! isign

!      OLD CODING
!      First: k + dk 
!      do ikpt1 = 1, dtefield%fnkpt
!      diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
!      &         dtefield%fkptns(:,ikpt) - dk(:))
!      if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
!      dtefield%ikpt_dk(ikpt,1,idir) = ikpt1
!      exit
!      end if
!      end do
       
!      Second: k - dk
!      do ikpt1 = 1, dtefield%fnkpt
!      diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
!      &         dtefield%fkptns(:,ikpt) + dk(:))
!      if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
!      dtefield%ikpt_dk(ikpt,2,idir) = ikpt1
!      exit
!      end if
!      end do 

     end do     ! ikpt

   end if ! end check on rfdir == 1

 end do     ! close loop over idir

 call timab(1005,2,tsec)
 call timab(1006,1,tsec)

!------------------------------------------------------------------------------
!------------ Compute PAW on-site terms if necessary --------------------------
!------------------------------------------------------------------------------
 
!write(std_out,'(a)')' JWZ Debug: dkvecs set to zero for testing!! '
!dtbfield%dkvecs = zero

 if (dtbfield%has_expibi == 1) then

   call expibi(dtbfield%expibi,dtbfield%dkvecs,gprimd,natom,rprimd,xred)
   call twexpibi(dtbfield%twexpibi,dtbfield%dkvecs,gprimd,&
&   dtbfield%my_natom,natom,rprimd,xred,&
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   dtbfield%has_expibi = 2

 end if

 if (dtbfield%has_qijb == 1) then

   call qijb_kk(dtbfield%qijb_kk,dtbfield%dkvecs,dtbfield%expibi,&
&   gprimd,dtbfield%lmn2max,natom,ntypat,pawang,pawrad,pawtab,typat)
   call twqijb_kk(dtbfield%twqijb_kk,dtbfield%dkvecs,dtbfield%twexpibi,&
&   gprimd,dtbfield%lmn2max,dtbfield%my_natom,natom,ntypat,pawang,pawrad,pawtab,typat,&
&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

 end if

 call set_twind(dtbfield)

 if (dtbfield%has_Lij == 1 .and. dtbfield%has_Lijr3 == 1 ) then

   call Lij(dtbfield,ntypat,pawrad,pawtab)

 end if

!
 write(std_out,'(a)')' making pawtwdij0 terms for magnetic field ... '
 write(std_out,'(a)')' term 1 for magnetic field ... '
!pawtwdij_1 has been checked at dk(:,:) = 0 against output of atompaw, gives same result
 call pawtwdij_1(dtbfield,gprimd,psps%mpsang,natom,ntypat,pawrad,pawtab,typat,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 write(std_out,'(a)')' term 2b for magnetic field ... '
!pawtwdij_2b has been checked at dk(:,:) = 0 against output of atompaw, gives same result
 call pawtwdij_2b(dtbfield,gprimd,psps%mpsang,natom,ntypat,pawrad,pawtab,typat,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 write(std_out,'(a)')' term 2e for magnetic field ... '
!pawtwdij_2e has been checked at dk(:,:) = 0 against output of atompaw, gives same result
 call pawtwdij_2e(dtbfield,gprimd,natom,ntypat,pawrad,pawtab,typat,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 dtbfield%has_twdij0 = 2
 write(std_out,'(a)')' done making pawtwdij0 terms for magnetic field ... '
!
!!  all four terms checked perfect at dk = 0 
 write(std_out,'(a)')' making pawtweijkl terms for magnetic field ... '
!
 write(std_out,'(a)')' term 2a* for magnetic field ... '
 call pawtwdij_2a(dtbfield,gprimd,natom,ntypat,pawrad,pawtab,typat,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 write(std_out,'(a)')' term 2c* for magnetic field ... '
 call pawtwdij_2c(dtbfield,gprimd,natom,ntypat,pawrad,pawtab,typat,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 write(std_out,'(a)')' term 2d* for magnetic field ... '
 call pawtwdij_2d(dtbfield,gprimd,natom,ntypat,pawrad,pawtab,typat,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 write(std_out,'(a)')' term 2f* for magnetic field ... '
 call pawtwdij_2f(dtbfield,gprimd,natom,ntypat,pawrad,pawtab,typat,&
& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
!!
 dtbfield%has_tweijkl = 2
!
 write(std_out,'(a)')' done making pawtweijkl terms for magnetic field ... '

 call timab(1007,2,tsec)
 call timab(1008,1,tsec)

!------------------------------------------------------------------------------
!------------ Build the array pwind that is needed to compute the -------------
!------------ overlap matrices at k +- dk                         -------------
!------------------------------------------------------------------------------

 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0
 pwind(:,:,:) = 0
 pwnsfac(1,:) = 1.0_dp
 pwnsfac(2,:) = 0.0_dp
 ABI_ALLOCATE(kg1_k,(3,mpw))

 ipwnsfac = 0

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     dk(:) = dtbfield%dkvecs(:,idir)

     do ifor = 1, 2

       if (ifor == 2) dk(:) = -1._dp*dk(:)

!      Build pwind and kgindex
!      NOTE: The array kgindex is important for parallel execution.
!      In case nsppol = 2, it may happent that a particular processor
!      treats k-points at different spin polarizations.
!      In this case, it is not possible to address the elements of
!      pwind correctly without making use of the kgindex array.

       ikg = 0 ; ikpt_loc = 0 ; isppol = 1
       do ikpt = 1, dtbfield%fnkpt

         ikpti = dtbfield%indkk_f2ibz(ikpt,1)
         nband_k = dtset%nband(ikpti)
         ikpt1f = dtbfield%ikpt_dk(ikpt,ifor,idir)
         ikpt1i = dtbfield%indkk_f2ibz(ikpt1f,1)

         if (xmpi_paral== 1) then
           if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,1,me)).and.&
&           (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,nsppol,me))) cycle

           ikpt_loc = ikpt_loc + 1

         end if

!        Build basis sphere of plane waves for the nearest neighbour of
!        the k-point (important for MPI //)

         kg1_k(:,:) = 0
         kpt1(:) = dtset%kptns(:,ikpt1i)
         call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg1_k,kpt1,&
&         1,mpi_enreg,mpw,npw_k1)
         me_g0=mpi_enreg%me_g0


!        ji: fkgindex is defined here !
         dtbfield%fkgindex(ikpt) = ikg

!        
!        Deal with symmetry transformations
!        

!        bra k-point k(b) and IBZ k-point kIBZ(b) related by
!        k(b) = alpha(b) S(b)^t kIBZ(b) + G(b)
!        where alpha(b), S(b) and G(b) are given by indkk_f2ibz
!        
!        For the ket k-point:
!        k(k) = alpha(k) S(k)^t kIBZ(k) + G(k) - GBZ(k)
!        where GBZ(k) takes k(k) to the BZ
!        

         isym  = dtbfield%indkk_f2ibz(ikpt,2)
         isym1 = dtbfield%indkk_f2ibz(ikpt1f,2)

!        Construct transformed G vector that enters the matching condition:
!        alpha(k) S(k)^{t,-1} ( -G(b) - GBZ(k) + G(k) )

         dg(:) = -dtbfield%indkk_f2ibz(ikpt,3:5) &
&         -nint(-dtbfield%fkptns(:,ikpt) - dk(:) - tol10 + &
&         dtbfield%fkptns(:,ikpt1f)) &
&         +dtbfield%indkk_f2ibz(ikpt1f,3:5)

!        old code
!        iadum(:)=0
!        do idum=1,3
!        iadum(:)=iadum(:)+ symrec(:,idum,isym1)*dg(idum)
!        end do

!        new code
         iadum(:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),dg(:))

         dg(:) = iadum(:)

         if ( dtbfield%indkk_f2ibz(ikpt1f,6) == 1 ) dg(:) = -dg(:)

!        Construct S(k)^{t,-1} S(b)^{t}

         dum33(:,:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),symrec(:,:,isym))

!        Construct alpha(k) alpha(b)

         if (dtbfield%indkk_f2ibz(ikpt,6) == dtbfield%indkk_f2ibz(ikpt1f,6)) then
           itrs=0
         else
           itrs=1
         end if


         npw_k  = npwarr(ikpti)
!        npw_k1 = npwarr(ikpt1i)

!        loop over bra G vectors
         do ipw = 1, npw_k

!          NOTE: the bra G vector is taken for the sym-related IBZ k point,
!          not for the FBZ k point
           iadum(:) = kg(:,dtbfield%kgindex(ikpti) + ipw)

!          Store non-symmorphic operation phase factor exp[i2\pi \alpha G \cdot t]

           if ( ipwnsfac == 0 ) then
!            old code
             rdum=0.0_dp
             do idum=1,3
               rdum=rdum+dble(iadum(idum))*dtset%tnons(idum,isym)
             end do
             rdum=two_pi*rdum
             if ( dtbfield%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
             pwnsfac(1,ikg+ipw) = cos(rdum)
             pwnsfac(2,ikg+ipw) = sin(rdum)
!            
!            new code
!            rdum = DOT_PRODUCT(dble(iadum(:)),dtset%tnons(:,isym))
!            rdum= two_pi*rdum
!            if ( dtefield%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
!            pwnsfac(1,ikg+ipw) = cos(rdum)
!            pwnsfac(2,ikg+ipw) = sin(rdum)

           end if

!          to determine r.l.v. matchings, we transformed the bra vector
!          Rotation
           iadum1(:)=0
           do idum1=1,3
             iadum1(:)=iadum1(:)+dum33(:,idum1)*iadum(idum1)
           end do
           iadum(:)=iadum1(:)
!          Time reversal
           if (itrs==1) iadum(:)=-iadum(:)
!          Translation
           iadum(:) = iadum(:) + dg(:)

           do jpw = 1, npw_k1
             iadum1(1:3) = kg1_k(1:3,jpw)
             if ( (iadum(1) == iadum1(1)).and. &
&             (iadum(2) == iadum1(2)).and. &
&             (iadum(3) == iadum1(3)) ) then
               pwind(ikg + ipw,ifor,idir) = jpw
!              write(std_out,'(a,2x,3i4,2x,i4)') 'Found !:',iadum1(:),jpw
               exit
             end if
           end do
         end do

         ikg  = ikg + npw_k

       end do    ! close loop over ikpt

       ipwnsfac = 1

     end do    ! close loop over ifor

   end if      ! rfdir(idir) == 1

 end do        ! close loop over idir


 call timab(1008,2,tsec)
 call timab(1009,1,tsec)

!Build mpi_enreg%kptdstrb
!array required to communicat the WFs between cpus in berryphase_new.f
!(MPI // over k-points)

 if (xmpi_paral== 1) then
   do idir = 1, 3
     if (dtset%rfdir(idir) == 1) then
       do ifor = 1, 2

         ikpt_loc = 0
         do isppol = 1, nsppol

           do ikpt = 1, dtbfield%fnkpt

             ikpti = dtbfield%indkk_f2ibz(ikpt,1)
             nband_k = dtset%nband(ikpti)
             ikpt1f = dtbfield%ikpt_dk(ikpt,ifor,idir)
             ikpt1i = dtbfield%indkk_f2ibz(ikpt1f,1)

             if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle
             
             ikpt_loc = ikpt_loc + 1
             mpi_enreg%kptdstrb(me + 1,ifor+2*(idir-1),ikpt_loc) = &
&             ikpt1i + (isppol - 1)*nkpt

             mpi_enreg%kptdstrb(me+1,ifor+2*(idir-1),&
&             ikpt_loc+dtbfield%fmkmem_max*nsppol) = &
&             ikpt1f + (isppol - 1)*dtbfield%fnkpt

           end do   ! ikpt
         end do     ! isppol
       end do       ! ifor
     end if         ! dtset%rfdir(idir) == 1
   end do           ! idir
 end if             ! xmpi_paral == 1

!build mpi_enreg%kpt_loc2fbz_sp 
 ikpt_loc = 0
 do isppol = 1, nsppol
   do ikpt = 1, dtbfield%fnkpt

     ikpti = dtbfield%indkk_f2ibz(ikpt,1)
     nband_k = dtset%nband(ikpti)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle
     
     ikpt_loc = ikpt_loc + 1

     mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 1) = ikpt
     mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 2) = isppol

   end do
 end do


!parallel case only :
!build mpi_enreg%kpt_loc2ibz_sp, dtefield%cgqindex and dtefield%nneigh
 if (xmpi_paral== 1) then
   ikpt_loc = 0
   do isppol = 1, nsppol
     do ikpt = 1, nkpt

       ikptf = dtbfield%i2fbz(ikpt)
       nband_k = dtset%nband(ikpti)

       neigh(:) = 0 ; icg = 0 ; ikg = 0 ; flag_kpt = 0; icprj = 0
       do idir=1, 3

         do ifor = 1, 2

           flag = 0

           ikpt1f = dtbfield%ikpt_dk(ikptf,ifor,idir)
           ikpt1i = dtbfield%indkk_f2ibz(ikpt1f,1)

           dtbfield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ikg
           ikg = ikg + npwarr(ikpt1i)

!          check if this neighbour is also a previous neighbour
           do ineigh = 1, (ifor+2*(idir-1))
             if (neigh(ineigh) == ikpt1i) then
               flag = 1
               dtbfield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ineigh
               dtbfield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
&               dtbfield%cgqindex(2,ineigh,ikpt+(isppol-1)*nkpt)
               exit
             end if
           end do
!          create the cgqindex of the neighbour if necessary
           if (flag == 0) then
             neigh(ifor+2*(idir-1)) = ikpt1i
             dtbfield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
&             ifor+2*(idir-1)
             dtbfield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = icg
             if (isppol == 1) dtbfield%nneigh(ikpt) = dtbfield%nneigh(ikpt) + 1
             icg = icg + npwarr(ikpt1i)*dtbfield%nspinor*nband_k
           end if

!          if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
!          !ikpt is not one of my kpt_loc
!          cycle
!          end if

         end do !ifor
       end do !idir

       if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me))) then
!        ikpt is one of my kpt_loc
         ikpt_loc = ikpt_loc + 1

         mpi_enreg%kpt_loc2ibz_sp(me, ikpt_loc, 1) = ikpt
         mpi_enreg%kpt_loc2ibz_sp(me, ikpt_loc, 2) = isppol
       end if

     end do !ikpt

   end do !isppol

 end if !xmpi_paral== 1

!should be temporary
!unassigned mpi_enreg%kpt_loc2fbz_sp are empty ; inform other cpu (there are better ways...)
 mpi_enreg%mkmem(me) = mkmem
!do ii=ikpt_loc+1,dtefield%fmkmem_max
!mpi_enreg%kpt_loc2fbz_sp(me, ii, 1) = -1
!end do


!(same as mpi_enreg%kptdstrb but for k-points in the iBZ),
!dtefield%cgqindex and dtefield%nneigh

 if (xmpi_paral== 1) then

   ikpt_loc = 1
   do isppol = 1, nsppol
     do ikpt = 1, nkpt

       nband_k = dtset%nband(ikpt)
       ikptf = dtbfield%i2fbz(ikpt)

       neigh(:) = 0 ; icg = 0 ; ikg = 0 ; flag_kpt = 0; icprj = 0
       do idir = 1, 3

         do ifor = 1, 2

           ikpt1f = dtbfield%ikpt_dk(ikptf,ifor,idir)
           ikpt1i = dtbfield%indkk_f2ibz(ikpt1f,1)

!          dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ikg
           ikg = ikg + npwarr(ikpt1i)

           flag = 0
           do ineigh = 1, (ifor+2*(idir-1))
             if (neigh(ineigh) == ikpt1i) then
               flag = 1
!              dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ineigh
!              dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
!              &               dtefield%cgqindex(2,ineigh,ikpt+(isppol-1)*nkpt)
               exit
             end if
           end do
           if (flag == 0) then
!            neigh(ifor+2*(idir-1)) = ikpt1i
!            dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
!            &             ifor+2*(idir-1)
!            dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = icg
!            if (isppol == 1) dtefield%nneigh(ikpt) = dtefield%nneigh(ikpt) + 1
!            icg = icg + npwarr(ikpt1i)*dtset%nspinor*nband_k
           end if

           if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle

           flag_kpt = 1

!          MVeithen: the if condition allows to avoid that the same wavefunction
!          is send several times to a particular cpu

         end do    ! ifor
       end do    ! idir

       if (flag_kpt == 1) ikpt_loc = ikpt_loc + 1

     end do    ! ikpt
   end do    ! isppol

 end if   ! fieldflag and xmpi-paral

 call xmpi_sum(mpi_enreg%kptdstrb,spaceComm,ierr)
 call xmpi_sum(mpi_enreg%kpt_loc2fbz_sp,spaceComm,ierr)
 call xmpi_sum(mpi_enreg%kpt_loc2ibz_sp,spaceComm,ierr)
 call xmpi_sum(mpi_enreg%mkmem,spaceComm,ierr)

 ABI_DEALLOCATE(kg1_k)

!DEBUG
!write(std_out,*)'initorbmag: exit'
!ENDDEBUG
 call timab(1009,2,tsec)
 call timab(1001,2,tsec)

end subroutine initorbmag
!!***
