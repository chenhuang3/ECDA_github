!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ebands
!! NAME
!!  m_ebands
!!
!! FUNCTION
!!  This module contains utilities to analyze and retrieve information
!!  from the ebands_t.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! TODO
!! 1) Remove npwarr, istwfk.
!! 2) Use 3d arrays for Bands%nband
!! 3) Solve issue with Hdr dependency
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_ebands

 use defs_basis
 use m_errors
 use m_profiling_abi
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 use defs_datatypes,   only : ebands_t
 use defs_abitypes,    only : hdr_type
 use m_copy,           only : alloc_copy
 use m_fstrings,       only : tolower, itoa, sjoin
 use m_numeric_tools,  only : arth, imin_loc, imax_loc, bisect, stats_t, stats_eval, simpson_int
 use m_special_funcs,  only : gaussian
 use m_tetrahedron,    only : get_tetra_weight
 use m_geometry,       only : normv
 use m_crystal,        only : crystal_t
 use m_header,         only : hdr_free, hdr_get_nelect_byocc, hdr_io_etsf
 use m_bz_mesh,        only : kmesh_t, isamek

 implicit none

 private

 public :: ebands_init             ! Main creation method.
 public :: ebands_from_hdr         ! Init from abinit header.
 public :: ebands_free             ! Destruction method.
 public :: ebands_copy             ! Deep copy of the ebands_t.
 public :: ebands_print            ! Printout basic info on the data type.
 public :: unpack_eneocc           ! Helper function for reshaping (energies|occupancies|derivate of occupancies).
 public :: pack_eneocc             ! Helper function for reshaping (energies|occupancies|derivate of occupancies).
 public :: get_eneocc_vect         ! Reshape (ene|occ|docdde) returning a matrix instead of a vector.
 public :: put_eneocc_vect         ! Put (ene|occ|doccde) in vectorial form into the data type doing a reshape.
 public :: get_bandenergy          ! Returns the band energy of the system.
 public :: get_valence_idx         ! Gives the index of the (valence|bands at E_f).
 public :: apply_scissor           ! Apply a scissor operator (no k-dependency)
 public :: get_occupied            ! Returns band indeces after wich occupations are less than an input value.
 public :: enclose_degbands        ! Adjust band indeces such that all degenerate states are treated.
 public :: get_nelect_per_spin     ! Returns number of electrons per spin channel
 public :: get_minmax              ! Returns min and Max value of (eig|occ|doccde).
 public :: bands_edstats           ! Compute statistical parameters of the energy differences e_ks[b+1] - e_ks[b]
 public :: get_FS                  ! Returns the Fermi surface (k-points and bands at E_f).
 public :: bst_metallic_scheme     ! .TRUE. if metallic occupation scheme is used.
 public :: ebands_3dprint          ! Write 3D energies for Fermi surface visualization (XSF format)
 public :: update_occ              ! Update the occupation numbers.
 public :: ReportGap               ! Print info on the fundamental and optical gap.
 public :: ExpandBands             ! Returns a new object defined in the BZ starting from the IBZ.
 public :: ebands_ncwrite          ! Dump the object into NETCDF file.
 public :: ebands_ncread           ! Initialize the object from a NETCDF file.
!!***

!----------------------------------------------------------------------

!!****t* m_ebands/edos_t
!! NAME
!! edos_t
!! 
!! FUNCTION
!! Store the electron DOS 
!! 
!! SOURCE

 type,public :: edos_t

   integer :: nsppol
    ! Number of spins.

   integer :: nw
   ! Number of points in the frequency mesh.

   integer :: ief=0
   ! Rightmost Index of the energy mesh such as IDOS[mesh[ief]] < nelect.
   ! 0 if Fermi level could not be computed
   ! Note the value of gef stored in edos_t is computed by performing 
   ! a linear interpolation between ief and ief+1

   real(dp) :: broad=zero
   ! Gaussian broadening

   real(dp) :: step
   ! Step of the mesh

   character(len=500) :: method
   ! gaussian or tetra

   real(dp),allocatable :: mesh(:)
   ! mesh(nw)

   real(dp),allocatable :: dos(:,:)
   ! dos(nw,0:nsppol)
   ! Total DOS, spin up and spin down component.

   real(dp),allocatable :: idos(:,:)
   ! idos(nw,0:nsppol)
   ! Integrated DOS, spin up and spin down component.

   !real(dp),allocatable gef(:)
   ! gef(0:nsppol)
   ! DOS at the Fermi level. Total, spin up, spin down

   type(ebands_t),pointer :: Bands => null()
   ! Reference to the bandstructure.

   type(kmesh_t),pointer :: Kmesh => null()
   ! Reference to the Kpoint mesh
 end type edos_t

 public :: edos_init     ! Compute the dos from the band structure.
 public :: edos_free     ! Free memory
 public :: edos_print    ! Print results to file (formatted mode)
!!***

!----------------------------------------------------------------------

!!****t* m_ebands/gaps_t
!! NAME
!! gaps_t
!! 
!! FUNCTION
!! Structure with information on the fundamental and optical gaps returned by ReportGap.
!! 
!! SOURCE

 type,public :: gaps_t

   integer :: nsppol
    ! Number of spins.

   integer,allocatable :: fo_kpos(:,:)
    ! fo_kpos(3,nsppol)
    ! fo_kpos(1:2,spin) ==> Indices of the k-points where the homo, lumo states are located (for each spin).
    ! fo_kpos(3,spin)   ==> the index of k-point where the optical gap is located (for each spin).

   integer,allocatable :: ierr(:)
     ! The third index corresponds to a "status" : 0.0dp if gaps were not computed
     ! (because there are only valence bands) ; -1.0dp if the system (or spin-channel) is metallic ; 1.0dp if the gap was computed

   real(dp),allocatable :: fo_values(:,:)
     ! fo_values(2,nsppol)]
     ! Fundamental and optical gaps (in Hartree) for each spin.

   real(dp),pointer :: kpoints(:,:) => null()
     ! Reference to the k-points of the band structure used to compute the gaps.

   character(len=500),allocatable :: errmsg_spin(:)
     ! errmsg_spin(nsppol)
     ! String with human-readable error messages if ierr(spin) != 0.

 end type gaps_t

 public :: get_gaps      ! Build the object from a bandstructure.
 public :: gaps_free     ! Free the structure.
 public :: gaps_print    ! Print info on the gaps 
!!***


CONTAINS  !=====================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_gaps
!! NAME
!! get_gaps
!!
!! FUNCTION
!!  Returns a structure with info on the fundamental and optical gap.
!!
!! INPUTS
!!  Bands<ebands_t>=Info on the band structure, the smearing technique and the physical temperature used.
!!  [kmask]=Logical mask used to exclude k-points.
!!
!! OUTPUT
!!  retcode=Return code (!=0 signals failure)
!!  gaps<gaps_t>=object with info on the gaps (parent is responsible for freeing the object).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_gaps(Bands,gaps,kmask) result(retcode)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_gaps'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),target,intent(in)  :: Bands
 type(gaps_t),intent(out) :: gaps
!arrays
 logical,optional,intent(in) :: kmask(Bands%nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikibz,nband_k,spin,nsppol,ikopt,ivk,ick,ivb,icb,retcode
 real(dp),parameter :: tol_fermi=tol6
 real(dp) :: fun_gap,opt_gap
 logical :: ismetal
!arrays
 integer :: val_idx(Bands%nkpt,Bands%nsppol)
 real(dp) :: top_valence(Bands%nkpt),bot_conduct(Bands%nkpt) 
 logical :: my_kmask(Bands%nkpt)

! *********************************************************************

 nsppol = Bands%nsppol

 ! Initialize gaps_t
 gaps%nsppol = nsppol
 ABI_MALLOC(gaps%fo_kpos, (3,nsppol))
 ABI_MALLOC(gaps%ierr, (nsppol))
 ABI_MALLOC(gaps%fo_values, (2, nsppol))
 ABI_MALLOC(gaps%errmsg_spin, (nsppol))
 gaps%kpoints => Bands%kptns

 gaps%fo_kpos = 0
 gaps%ierr = 0
 gaps%fo_values = zero
 gaps%errmsg_spin(:) = ""

 my_kmask=.TRUE.; if (PRESENT(kmask)) my_kmask=kmask

 val_idx(:,:) = get_valence_idx(Bands,tol_fermi)

 spin_loop: &
&  do spin=1,nsppol

   ! No output if system i metallic
   ismetal=ANY(val_idx(:,spin)/=val_idx(1,spin)) 
   if (ismetal) then
     gaps%ierr(spin) = 1
     write(gaps%errmsg_spin(spin), "(a,i0)")"Metallic system for spin channel ",spin
     CYCLE 
   endif

   ivb=val_idx(1,spin)
   icb=ivb+1

   do ikibz=1,Bands%nkpt
     if (.not.my_kmask(ikibz)) CYCLE
     nband_k=Bands%nband(ikibz+(spin-1)*Bands%nkpt)
     top_valence(ikibz)=Bands%eig(ivb,ikibz,spin)
     if (icb>nband_k) then
       gaps%ierr(spin) = 2
       gaps%errmsg_spin(spin) = "Not enough states to calculate the band gap."
       CYCLE spin_loop
     endif
     bot_conduct(ikibz)=Bands%eig(icb,ikibz,spin)
   end do

   ! Minimum of the optical Gaps
   ikopt= imin_loc(bot_conduct-top_valence,MASK=my_kmask)
   opt_gap=bot_conduct(ikopt)-top_valence(ikopt)

   ! Fundamental Gap ===
   ick = imin_loc(bot_conduct,MASK=my_kmask)
   ivk = imax_loc(top_valence,MASK=my_kmask)
   fun_gap = Bands%eig(icb,ick,spin)-Bands%eig(ivb,ivk,spin)

   gaps%fo_values(:,spin) = [fun_gap, opt_gap]
   gaps%fo_kpos(:,spin) = [ivk, ick, ikopt]
 end do spin_loop

 retcode = MAXVAL(gaps%ierr)

end function get_gaps
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/gaps_free
!! NAME
!!  gaps_free 
!!
!! FUNCTION
!!  Free the memory allocated in gaps_t
!!
!! PARENTS
!!      setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine gaps_free(gaps)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaps_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(gaps_t),intent(inout) :: gaps

! *********************************************************************

 !@gaps_t

!integer
 if (allocated(gaps%fo_kpos)) then
   ABI_FREE(gaps%fo_kpos)
 end if

 if (allocated(gaps%ierr)) then
   ABI_FREE(gaps%ierr)
 end if

!real
 if (allocated(gaps%fo_values)) then
   ABI_FREE(gaps%fo_values)
 end if

!chars
 if (allocated(gaps%errmsg_spin)) then
   ABI_FREE(gaps%errmsg_spin)
 end if

! nullify pointers
 nullify(gaps%kpoints)

end subroutine gaps_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/gaps_print
!! NAME
!! gaps_print
!!
!! FUNCTION
!!  Print info on the fundamental and optical gap.
!!
!! INPUTS
!!  gaps<gaps_t>=Object with info on the gaps.
!!  [header]=Optional title.
!!  [unit]=Optional unit for output (std_out if not specified)
!!  [mode_paral]=Either "COLL" or "PERS", former is default.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine gaps_print(gaps,header,unit,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaps_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 type(gaps_t),intent(in)  :: gaps

!Local variables-------------------------------
!scalars
 integer :: spin,ikopt,ivk,ick,my_unt
 real(dp) :: fun_gap,opt_gap
 character(len=4) :: my_mode
 character(len=500) :: msg

! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

 do spin=1,gaps%nsppol

   if (spin==1) then
     msg=ch10
     if (PRESENT(header)) msg=ch10//' === '//TRIM(ADJUSTL(header))//' === '
     call wrtout(my_unt,msg,my_mode) 
   end if

   if (gaps%ierr(spin) /= 0) then
     call wrtout(my_unt,gaps%errmsg_spin(spin), my_mode)
     continue
   end if

   ! Get minimum of the optical Gap.
   fun_gap = gaps%fo_values(1,spin)
   opt_gap = gaps%fo_values(2,spin)

   if (any(gaps%fo_kpos(:,spin) == 0)) then
     call wrtout(my_unt,sjoin("Cannot detect gap for spin: ",itoa(spin)),"COLL")
     cycle
   end if

   ivk = gaps%fo_kpos(1,spin)
   ick = gaps%fo_kpos(2,spin)
   ikopt = gaps%fo_kpos(3,spin)

   write(msg,'(a,i2,a,2(a,f8.4,a,3f8.4,a),33x,a,3f8.4)')&
&    '  >>>> For spin ',spin,ch10,&
&    '   Minimum optical gap = ',opt_gap*Ha_eV,' [eV], located at k-point      : ',gaps%kpoints(:,ikopt),ch10,&
&    '   Fundamental gap     = ',fun_gap*Ha_eV,' [eV], Top of valence bands at : ',gaps%kpoints(:,ivk),ch10,  &
&                                              '       Bottom of conduction at : ',gaps%kpoints(:,ick)
   call wrtout(my_unt,msg,my_mode) 
 end do !spin

end subroutine gaps_print
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_init
!! NAME
!! ebands_init
!!
!! FUNCTION
!! This subroutine initializes the ebands_t structured datatype
!!
!! INPUTS
!! bantot=total number of bands (=sum(nband(:))
!! doccde(bantot)=derivative of the occupation numbers with respect to the energy (Ha)
!! eig(bantot)=eigenvalues (hartree)
!! istwfk(nkpt)=parameter that describes the storage of wfs.
!! kptns(3,nkpt)=k points in terms of recip primitive translations
!! nband(nkpt*nsppol)=number of bands
!! nelect=Number of electrons.
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves at each k point
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nspinor=Number of spinor components
!! occopt=Occupation options (see input variable)
!! occ(bantot)=occupation numbers
!! tphysel=Physical temperature (input variable)
!! tsmear=Temperature of smearing.
!! wtk(nkpt)=weight assigned to each k point
!!
!! OUTPUT
!! Bands<ebands_t>=the ebands_t datatype
!!
!! SIDE EFFECTS
!!  %entropy and %fermie initialized to zero.
!!
!! PARENTS
!!      eig2tot,gstate,loper3,m_ebands,m_shirley,mlwfovlp_qp,nonlinear,optic
!!      outscfcv,respfn,setup_bse,setup_bse_interp,setup_screening,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_init(bantot,Bands,nelect,doccde,eig,istwfk,kptns,&
& nband,nkpt,npwarr,nsppol,nspinor,tphysel,tsmear,occopt,occ,wtk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot,nkpt,nsppol,nspinor,occopt
 real(dp),intent(in) :: nelect,tphysel,tsmear
 type(ebands_t),intent(out) :: Bands
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt)
 real(dp),intent(in) :: doccde(bantot),eig(bantot),kptns(3,nkpt),occ(bantot)
 real(dp),intent(in) :: wtk(nkpt)
! *************************************************************************

 ! Copy the scalars
 ! MG TODO here there is a inconsistency in the way occ are treated in the header
 ! (only the states used, bantot. are saved, and the way occ. and energies
 ! are passed to routines (mband,nkpt,nsppol). It might happen that bantot<mband*nktp*nsppol
 ! this should not lead to problems since arrays are passed by reference
 ! anyway the treatment of these arrays have to be rationalized
 Bands%bantot =bantot
 Bands%mband  =MAXVAL(nband(1:nkpt*nsppol))
 Bands%nkpt   =nkpt
 Bands%nspinor=nspinor
 Bands%nsppol =nsppol
 Bands%occopt =occopt
 
 Bands%entropy=zero  ! Initialize results
 Bands%fermie =zero  ! Initialize results
 Bands%nelect =nelect
 Bands%tphysel=tphysel
 Bands%tsmear =tsmear

 ! Allocate the components
 ABI_MALLOC(Bands%nband,(nkpt*nsppol))
 ABI_MALLOC(Bands%istwfk,(nkpt))
 ABI_MALLOC(Bands%npwarr,(nkpt))
 ABI_MALLOC(Bands%kptns,(3,nkpt))

 ! Copy the arrays
 Bands%nband(1:nkpt*nsppol)=nband(1:nkpt*nsppol)
 Bands%istwfk(1:nkpt)      =istwfk(1:nkpt)
 Bands%npwarr(1:nkpt)      =npwarr(1:nkpt)
 Bands%kptns(1:3,1:nkpt)   =kptns(1:3,1:nkpt)

 ! In Bands, energies and occupations are stored in a matrix (mband,nkpt,nsppol).
 ! put_eneocc_vect is used to reshape the values stored in vectorial form.
 ABI_MALLOC(Bands%eig   ,(Bands%mband,nkpt,nsppol))
 ABI_MALLOC(Bands%occ   ,(Bands%mband,nkpt,nsppol))
 ABI_MALLOC(Bands%doccde,(Bands%mband,nkpt,nsppol))
 Bands%eig=HUGE(one); Bands%occ=zero; Bands%doccde=zero

 call put_eneocc_vect(Bands,'eig',   eig   ) 
 call put_eneocc_vect(Bands,'occ',   occ   ) 
 call put_eneocc_vect(Bands,'doccde',doccde) 

 ABI_MALLOC(Bands%wtk,(nkpt))
 Bands%wtk(1:nkpt)=wtk(1:nkpt)

end subroutine ebands_init
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_from_hdr
!! NAME
!! ebands_from_hdr
!!
!! FUNCTION
!! This subroutine initializes the ebands_t datatype from the abinit header by
!! calling the main creation method.
!!
!! INPUTS
!!  Hdr<hdr_type>=Abinit header.
!!  mband=Maximun number of bands.
!!  ene3d(mband,Hdr%nkpt,Hdr%nsppol)=Energies.
!!  [nelect]=Number of electrons per unit cell.
!!    Optional argument that can be used for performing a ridid shift of the fermi level. 
!!    in the case of metallic occupancies.
!!    If not specified, nelect will be initialized from Hdr.
!!
!! OUTPUT
!!  Bands<ebands_t>=The ebands_t datatype completely initialized.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_from_hdr(Bands,Hdr,mband,ene3d,nelect)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_from_hdr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband
 type(hdr_type),intent(in) :: Hdr
 type(ebands_t),intent(out) :: Bands
 real(dp),optional,intent(in) :: nelect
!arrays
 real(dp),intent(in) :: ene3d(mband,Hdr%nkpt,Hdr%nsppol)

!Local variables-------------------------------
!scalars
 real(dp) :: my_nelect
!arrays
 real(dp),allocatable :: ugly_doccde(:),ugly_ene(:)
! *************************************************************************

 if (PRESENT(nelect)) then
   my_nelect = nelect
 else 
   ! TODO Have to add nelect to the header
   my_nelect = hdr_get_nelect_byocc(Hdr) 
 end if
 !
 ! Have to use ugly 1d vectors to call ebands_init
 ABI_MALLOC(ugly_doccde,(Hdr%bantot))
 ugly_doccde=zero

 ABI_MALLOC(ugly_ene,(Hdr%bantot))
 call pack_eneocc(Hdr%nkpt,Hdr%nsppol,mband,Hdr%nband,Hdr%bantot,ene3d,ugly_ene)

 call ebands_init(Hdr%bantot,Bands,my_nelect,ugly_doccde,ugly_ene,Hdr%istwfk,Hdr%kptns,Hdr%nband,Hdr%nkpt,&
&  Hdr%npwarr,Hdr%nsppol,Hdr%nspinor,Hdr%tphysel,Hdr%tsmear,Hdr%occopt,Hdr%occ,Hdr%wtk)

 ABI_FREE(ugly_doccde)
 ABI_FREE(ugly_ene)

end subroutine ebands_from_hdr
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_free
!! NAME
!! ebands_free
!!
!! FUNCTION
!! Deallocates the components of the ebands_t structured datatype
!!
!! INPUTS
!!  Bands<ebands_t>=The data type to be deallocated.
!!
!! OUTPUT
!!  Deallocate the dynamic arrays in the ebands_t type.
!!  (only deallocate)
!!
!! PARENTS
!!      bethe_salpeter,eig2tot,elphon,gstate,loper3,m_shirley,mlwfovlp_qp
!!      nonlinear,optic,outscfcv,respfn,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_free(Bands)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(inout) :: Bands
! *************************************************************************

 DBG_ENTER("COLL")

!Deallocate all components of bstruct
 if (allocated(Bands%istwfk)) then
   ABI_FREE(Bands%istwfk)
 end if

 if (allocated(Bands%nband)) then
   ABI_FREE(Bands%nband)
 end if

 if (allocated(Bands%npwarr)) then
   ABI_FREE(Bands%npwarr)
 end if

 if (allocated(Bands%kptns)) then
   ABI_FREE(Bands%kptns)
 end if

 if (allocated(Bands%eig)) then
   ABI_FREE(Bands%eig)
 end if

 if (allocated(Bands%occ)) then
   ABI_FREE(Bands%occ)
 end if

 if (allocated(Bands%doccde)) then
   ABI_FREE(Bands%doccde)
 end if

 if (allocated(Bands%wtk)) then
   ABI_FREE(Bands%wtk)
 end if

 DBG_EXIT("COLL")

end subroutine ebands_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_copy
!! NAME
!!  ebands_copy
!!
!! FUNCTION
!! This subroutine performs a deep copy of a ebands_t datatype.
!! All the associated pointers in the input object will be copied preserving the shape.
!! If a pointer in Bands_in happens to be not associated, the corresponding
!! pointer in the copied object will be nullified.
!!
!! INPUTS
!!  Bands_in<ebands_t>=The data type to be copied.
!!
!! OUTPUT
!!  Bands_cp<ebands_t>=The copy.
!!
!! TODO 
!!  To be on the safe side one should nullify all pointers in the ebands_t 
!!  in the creation method. We have to follow F90 specifications and the initial status 
!!  of a pointer is not defined. This might lead to problem in deep_copy.
!!
!! PARENTS
!!      screening,setup_bse,setup_bse_interp,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_copy(Bands_in,Bands_cp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_copy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in)  :: Bands_in
 type(ebands_t),intent(out) :: Bands_cp

! *********************************************************************

 ! Copy scalars
 Bands_cp%bantot  = Bands_in%bantot  
 Bands_cp%mband   = Bands_in%mband
 Bands_cp%nkpt    = Bands_in%nkpt    
 Bands_cp%nspinor = Bands_in%nspinor 
 Bands_cp%nsppol  = Bands_in%nsppol  
 Bands_cp%occopt  = Bands_in%occopt  

 Bands_cp%entropy = Bands_in%entropy 
 Bands_cp%fermie  = Bands_in%fermie  
 Bands_cp%nelect  = Bands_in%nelect  
 Bands_cp%tphysel = Bands_in%tphysel 
 Bands_cp%tsmear  = Bands_in%tsmear  

 ! Copy allocatable arrays
 call alloc_copy(Bands_in%istwfk, Bands_cp%istwfk)
 call alloc_copy(Bands_in%nband , Bands_cp%nband )     
 call alloc_copy(Bands_in%npwarr, Bands_cp%npwarr)    

 call alloc_copy(Bands_in%kptns , Bands_cp%kptns ) 
 call alloc_copy(Bands_in%eig   , Bands_cp%eig   )  
 call alloc_copy(Bands_in%occ   , Bands_cp%occ   )   
 call alloc_copy(Bands_in%doccde, Bands_cp%doccde)   
 call alloc_copy(Bands_in%wtk   , Bands_cp%wtk   )

end subroutine ebands_copy
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_print
!! NAME
!! ebands_print
!!
!! FUNCTION
!! Print the content of the object.
!!
!! INPUTS
!!  Bands<ebands_t>The type containing the data.
!!  [unit]=Unit number (std_out if None)
!!  [header]=title for info
!!  [prtvol]=Verbosity level (0 if None)
!!  [mode_paral]=Either 'COLL' or 'PERS' ('COLL' if None).
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      setup_bse,setup_bse_interp
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_print(Bands,header,unit,prtvol,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=*),optional,intent(in) :: header
 character(len=4),optional,intent(in) :: mode_paral
 type(ebands_t),intent(in) :: Bands

!Local variables-------------------------------
 integer :: spin,ikpt,my_unt,my_prtvol,ii
 character(len=4) :: my_mode
 character(len=500) :: msg
! *************************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the ebands_t ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(5(a,i5,a))')&
&  '  Number of spinorial components ...... ',Bands%nspinor,ch10,&
&  '  Number of spin polarizations ........ ',Bands%nsppol,ch10,&
&  '  Number of k-points in the IBZ ....... ',Bands%nkpt,ch10,&
&  '  Maximum number of bands ............. ',Bands%mband,ch10,&
&  '  Occupation option ................... ',Bands%occopt,ch10
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a,f14.2,a,4(a,f14.6,a))')&
&  '  Number of valence electrons ......... ',Bands%nelect,ch10,&
&  '  Fermi level  ........................ ',Bands%fermie,ch10,&
&  '  Entropy ............................. ',Bands%entropy,ch10,&
&  '  Tsmear value ........................ ',Bands%tsmear,ch10,&
&  '  Tphysel value ....................... ',Bands%tphysel,ch10
 call wrtout(my_unt,msg,my_mode)

 if (my_prtvol>0) then

   if (Bands%nsppol==1)then
     write(msg,'(a,i0,a)')' New occ. numbers for occopt= ',Bands%occopt,' , spin-unpolarized case. '
     call wrtout(my_unt,msg,my_mode)
   end if

   do spin=1,Bands%nsppol
     if (Bands%nsppol==2) then
       write(msg,'(a,i4,a,i2)')' New occ. numbers for occopt= ',Bands%occopt,' spin ',spin
       call wrtout(my_unt,msg,my_mode)
     end if

     do ikpt=1,Bands%nkpt
       write(msg,'(2a,i4,a,3f12.6,a,f6.3)')ch10,&
&        ' k-point number ',ikpt,') ',Bands%kptns(:,ikpt),'; weight: ',Bands%wtk(ikpt)
       call wrtout(my_unt,msg,my_mode)
       do ii=1,Bands%nband(ikpt+(spin-1)*Bands%nkpt)
         write(msg,'(3(f7.3,1x))')Bands%eig(ii,ikpt,spin)*Ha_eV,Bands%occ(ii,ikpt,spin),Bands%doccde(ii,ikpt,spin)
         call wrtout(my_unt,msg,my_mode)
       end do
     end do !ikpt

   end do !spin

   !TODO add additional info useful for debugging)
 end if !my_prtvol

end subroutine ebands_print
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/unpack_eneocc
!! NAME
!! unpack_eneocc
!!
!! FUNCTION
!!  Helper function to do a reshape of (energies|occupancies|derivate of occupancies)
!!  initially stored in a vector. Return a 3D array index by (band,ikpt,spin) 
!!
!! INPUTS
!!  nkpt=number of k-points
!!  nsppol=number of spin polarizations
!!  mband=Max number of bands over k-points (just to dimension the output)
!!  nbands(nkpt*nsppol)=Number of bands at eack k and spin
!!  vect(:)=The input values to reshape
!!
!! OUTPUT
!!  array3d(mband,nkpt,nsppol)=Arrays containing the values of vect. 
!!   Note that the first dimension is usually larger than the 
!!   number of bands really used for a particular k-point and spin.
!!
!! PARENTS
!!      cchi0q0_intraband,dfpt_write_cg,kss2wfk,m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine unpack_eneocc(nkpt,nsppol,mband,nband,vect,array3d)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unpack_eneocc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsppol,mband
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: vect(:)
 real(dp),intent(out) :: array3d(mband,nkpt,nsppol)

!Local variables-------------------------------
 integer :: spin,ikpt,band,idx
! *************************************************************************

 array3d=HUGE(zero)

 idx=0
 ! elements in vect are packed in the first positions.
 do spin=1,nsppol
   do ikpt=1,nkpt
     do band=1,nband(ikpt+(spin-1)*nkpt)
      idx=idx+1
      array3d(band,ikpt,spin)=vect(idx)
     end do
   end do
 end do

end subroutine unpack_eneocc
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/pack_eneocc
!! NAME
!! pack_eneocc
!!
!! FUNCTION
!!  Helper function to do a reshape of (energies|occupancies|derivate of occupancies)
!!  initially stored in a 3D arrays returning a vector. 
!!
!! INPUTS
!!  nkpt=number of k-points
!!  nsppol=number of spin polarizations
!!  mband=Max number of bands over k-points (just to dimension the output)
!!  nbands(nkpt*nsppol)=Number of bands at eack k and spin
!!  bantot=Total number of bands
!!  array3d(mband,nkpt,nsppol)=Arrays containing the values to reshape.
!!
!! OUTPUT
!!  vect(bantot)=The input values stored in vector mode. Only the values really
!!   considered at each k-point and spin are copied.
!!
!! PARENTS
!!      cchi0q0_intraband,m_ebands,m_shirley
!!
!! CHILDREN
!!
!! SOURCE

subroutine pack_eneocc(nkpt,nsppol,mband,nband,bantot,array3d,vect)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pack_eneocc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsppol,mband,bantot
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: array3d(mband,nkpt,nsppol)
 real(dp),intent(out) :: vect(bantot)

!Local variables-------------------------------
 integer :: spin,ikpt,band,idx

! *************************************************************************

 vect(:)=zero
 idx=0
 do spin=1,nsppol
   do ikpt=1,nkpt
     do band=1,nband(ikpt+(spin-1)*nkpt)
       idx=idx+1
       vect(idx)=array3d(band,ikpt,spin)
     end do
   end do
 end do

end subroutine pack_eneocc 
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_eneocc_vect
!! NAME
!! get_eneocc_vect
!!
!! FUNCTION
!!  Retrieve energies or occupations from a ebands_t structure accessing by name. 
!!  Results are reported in a vector to facilitate the interface with other abinit routines.
!!
!! INPUTS
!!  Bands<ebands_t>The type containing the data.
!!  arr_name=The name of the quantity to retrieve. Allowed values are
!!   == "eig"    == For the eigenvalues.
!!   == "occ"    == For the occupation numbers.
!!   == "doccde" == For the derivative of the occupancies wrt the energy.
!!
!! OUTPUT
!!  vect(Bands%bantot)=The values required.
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_eneocc_vect(Bands,arr_name,vect)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_eneocc_vect'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arr_name
 type(ebands_t),intent(in) :: Bands
 real(dp),intent(out) :: vect(Bands%bantot)

!Local variables-------------------------------
 integer :: nkpt,nsppol,mband,bantot
 character(len=500) :: msg
! *************************************************************************

 mband =Bands%mband; bantot=Bands%bantot; nkpt=Bands%nkpt; nsppol=Bands%nsppol

 SELECT CASE (arr_name)
 CASE ('occ')
   call pack_eneocc(nkpt,nsppol,mband,Bands%nband,bantot,Bands%occ,vect)
 CASE ('eig')
   call pack_eneocc(nkpt,nsppol,mband,Bands%nband,bantot,Bands%eig,vect)
 CASE ('doccde')
   call pack_eneocc(nkpt,nsppol,mband,Bands%nband,bantot,Bands%doccde,vect)
 CASE DEFAULT
   write(msg,'(2a)')' Wrong value of arr_name= ',TRIM(arr_name)
   MSG_BUG(msg)
 END SELECT

end subroutine get_eneocc_vect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/put_eneocc_vect
!! NAME
!! put_eneocc_vect
!!
!! FUNCTION
!!  Update the energies or the occupations stored in a ebands_t structure. 
!!  The input values are stored in a vector according to the abinit convention
!!  In the data type, on the contrary,  we use 3D arrays (mband,nkpt,nsspol) 
!!  which are much easier to use inside loops.
!!
!! INPUTS
!!  vect(Bands%bantot)=The new values to be stored in the structure.
!!  arr_name=The name of the quantity to be saved (CASE insensitive). 
!!  Allowed values are
!!   == "eig"    == For the eigenvalues.
!!   == "occ"    == For the occupation numbers.
!!   == "doccde" == For the derivative of the occupancies wrt the energy.
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  Bands<ebands_t>=The object with updated values depending on the value of arr_name
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine put_eneocc_vect(Bands,arr_name,vect) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'put_eneocc_vect'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arr_name
 type(ebands_t),intent(inout) :: Bands
 real(dp),intent(in) :: vect(Bands%bantot)

!Local variables-------------------------------
 integer :: nkpt,nsppol,mband,bantot
 character(len=500) :: msg
! *************************************************************************

 mband =Bands%mband 
 bantot=Bands%bantot
 nkpt  =Bands%nkpt
 nsppol=Bands%nsppol

 SELECT CASE (tolower(arr_name))
 CASE ('occ')
   call unpack_eneocc(nkpt,nsppol,mband,Bands%nband,vect,Bands%occ)
 CASE ('eig')
   call unpack_eneocc(nkpt,nsppol,mband,Bands%nband,vect,Bands%eig)
 CASE ('doccde')
   call unpack_eneocc(nkpt,nsppol,mband,Bands%nband,vect,Bands%doccde)
 CASE DEFAULT 
   write(msg,'(2a)')' Wrong value of arr_name= ',TRIM(arr_name)
   MSG_BUG(msg)
 END SELECT

end subroutine put_eneocc_vect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_bandenergy
!! NAME
!! get_bandenergy
!!
!! FUNCTION
!!  Return the band energy (weighted sum of occupied eigenvalues)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!! TODO Likely this expression is not accurate since it is not variatonal
!!  One should use 
!!   band_energy = \int e N(e) de   for e<Ef , where N(e) is the e-DOS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_bandenergy(Bands) result(band_energy)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_bandenergy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: Bands
 real(dp) :: band_energy

!Local variables-------------------------------
 integer :: spin,ikibz,nband_k
 real(dp) :: wtk
! *********************************************************************

 band_energy=zero
 do spin=1,Bands%nsppol
   do ikibz=1,Bands%nkpt
     wtk=Bands%wtk(ikibz)
     nband_k=Bands%nband(ikibz+(spin-1)*Bands%nkpt) 
     band_energy = band_energy + wtk*SUM( Bands%eig(1:nband_k,ikibz,spin)*Bands%occ(1:nband_k,ikibz,spin) )
   end do
 end do

end function get_bandenergy
!!***

!!****f* m_ebands/get_valence_idx
!! NAME
!!  get_valence_idx
!!
!! FUNCTION
!!  For each k-point and spin polarisation, report: 
!!   The index of the valence in case of Semiconductors.
!!   The index of the band at the Fermi energy+toldfe
!!
!! INPUTS
!!  Bands<ebands_t>=The object describing the band structure.
!!  tol_fermi[optional]
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_valence_idx(Bands,tol_fermi) result(val_idx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_valence_idx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: tol_fermi
 type(ebands_t),intent(in) :: Bands
!arrays
 integer :: val_idx(Bands%nkpt,Bands%nsppol)

!Local variables-------------------------------
 integer :: band,ikpt,spin,idx,nband_k 
 real(dp) :: tol_

! *************************************************************************

 tol_=tol6; if (PRESENT(tol_fermi)) tol_=tol_fermi

 do spin=1,Bands%nsppol
   do ikpt=1,Bands%nkpt
     nband_k=Bands%nband(ikpt+(spin-1)*Bands%nkpt)

     idx=0
     do band=1,nband_k
       if (Bands%eig(band,ikpt,spin) > Bands%fermie+ABS(tol_)) then
         idx=band; EXIT
       end if
     end do
     val_idx(ikpt,spin)=idx-1
     if (idx==1) val_idx(ikpt,spin)=idx
     if (idx==0) val_idx(ikpt,spin)=nband_k

   end do
 end do

end function get_valence_idx
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/apply_scissor
!! NAME
!!  apply_scissor
!!
!! FUNCTION
!!  Apply a scissor operator of amplitude scissor_energy.
!!
!! INPUTS
!!  scissor_energy=The energy shift
!!
!! OUTPUT
!!
!! SIDE EFFECT
!!  Bands<ebands_t>=The following quantities are modified:
!!   %eig(mband,nkpt,nsppol)=The band structure after the application of the scissor operator
!!   %fermi_energy
!!
!! PARENTS
!!      screening,setup_bse,setup_bse_interp
!!
!! CHILDREN
!!
!! SOURCE

subroutine apply_scissor(Bands,scissor_energy)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'apply_scissor'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: scissor_energy
 type(ebands_t),intent(inout) :: Bands

!Local variables-------------------------------
 integer :: ikpt,spin,ival,nband_k 
 real(dp) :: spinmagntarget_
 character(len=500) :: msg
!arrays
 integer :: val_idx(Bands%nkpt,Bands%nsppol)
! *************************************************************************

 ! === Get the valence band index for each k and spin ===
 val_idx(:,:) = get_valence_idx(Bands)

 do spin=1,Bands%nsppol
   if (ANY(val_idx(:,spin)/=val_idx(1,spin))) then
     write(msg,'(a,i2,a)')&
&     ' Trying to apply a scissor operator on a metallic band structure for spin = ',spin,&
&     ' Assuming you know what you are doing, continuing anyway! '
     MSG_COMMENT(msg)
     !Likely newocc will stop, unless the system is semimetallic ?
   end if
 end do

 ! === Apply the scissor ===
 do spin=1,Bands%nsppol
   do ikpt=1,Bands%nkpt
     nband_k=Bands%nband(ikpt+(spin-1)*Bands%nkpt)
     ival=val_idx(ikpt,spin)

     if (nband_k>=ival+1) then
       Bands%eig(ival+1:,ikpt,spin) = Bands%eig(ival+1:,ikpt,spin)+scissor_energy
     else 
       write(msg,'(2a,4(a,i4))')&
&       ' Not enough bands to apply the scissor operator. ',ch10,&
&       ' spin = ',spin,' ikpt = ',ikpt,' nband_k = ',nband_k,' but valence index = ',ival
       MSG_COMMENT(msg)
     end if

   end do
 end do

 ! === Recalculate the fermi level and occ. factors ===
 ! * For Semiconductors only the Fermi level is changed (in the middle of the new gap) 
 spinmagntarget_=-99.99_dp !?; if (PRESENT(spinmagntarget)) spinmagntarget_=spinmagntarget
 call update_occ(Bands,spinmagntarget_)

end subroutine apply_scissor
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_occupied
!! NAME
!!  get_occupied
!!
!! FUNCTION
!!  For each k-point and spin polarisation, report the band index
!!  after which the occupation numbers are less than tol_occ.
!!
!! INPUTS
!!  Bands<ebands_t>=The object describing the band structure.
!!  tol_occ[Optional]=Tollerance on the occupation factors.
!!
!! OUTPUT
!!
!! NOTES
!!  We assume that the occupation factors monotonically decrease as a function of energy.
!!  This is not always true for eavery smearing technique implemented in Abinit.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_occupied(Bands,tol_occ) result(occ_idx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_occupied'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: tol_occ
 type(ebands_t),intent(in) :: Bands
!arrays
 integer :: occ_idx(Bands%nkpt,Bands%nsppol)

!Local variables-------------------------------
 integer :: band,ikpt,spin,idx,nband_k 
 real(dp) :: tol_

! *************************************************************************

 tol_=tol8 ; if (PRESENT(tol_occ)) tol_=tol_occ

 do spin=1,Bands%nsppol
   do ikpt=1,Bands%nkpt
     nband_k=Bands%nband(ikpt+(spin-1)*Bands%nkpt)

     idx=0
     do band=1,nband_k
       if (Bands%occ(band,ikpt,spin)<ABS(tol_)) then
         idx=band; EXIT
       end if
     end do
     occ_idx(ikpt,spin)=idx-1
     if (idx==1) occ_idx(ikpt,spin)=idx
     if (idx==0) occ_idx(ikpt,spin)=nband_k

   end do
 end do

end function get_occupied
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/enclose_degbands
!! NAME
!!  enclose_degbands
!!
!! FUNCTION
!!  Adjust ibmin and ibmax such that all the degenerate states are enclosed
!!  between ibmin and ibmax. The routine works for a given k-point a spin.
!!
!! INPUTS
!!  Bands<ebands_t>=The object describing the band structure.
!!  ikibz=Index of the k-point.
!!  spin=Spin index.
!!  tol_enedif=Tolerance on the energy difference.
!!
!! OUTPUT 
!!  changed=.TRUE. if ibmin or ibmax has been changed.
!!
!! SIDE EFFECTS
!!  ibmin,ibmax=
!!    Input: initial guess for the indeces
!!    Output: All the denerate states are between ibmin and ibmax 
!!
!! PARENTS
!!      setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine enclose_degbands(Bands,ikibz,spin,ibmin,ibmax,changed,tol_enedif)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'enclose_degbands'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikibz,spin
 integer,intent(inout) :: ibmin,ibmax
 real(dp),intent(in) :: tol_enedif
 logical,intent(out) :: changed
 type(ebands_t),intent(in) :: Bands

!Local variables-------------------------------
!scalars
 integer :: ib,ibmin_bkp,ibmax_bkp
 real(dp) :: emin,emax

! *************************************************************************

 ibmin_bkp = ibmin
 ibmax_bkp = ibmax

 emin =  Bands%eig(ibmin,ikibz,spin)
 do ib=ibmin-1,1,-1
   if ( ABS(Bands%eig(ib,ikibz,spin) - emin) > tol_enedif) then
     ibmin = ib +1 
     EXIT 
   else 
     ibmin = ib
   end if
 end do

 emax =  Bands%eig(ibmax,ikibz,spin)
 do ib=ibmax+1,Bands%nband(ikibz+(spin-1)*Bands%nkpt)
   if ( ABS(Bands%eig(ib,ikibz,spin) - emax) > tol_enedif) then
     ibmax = ib - 1 
     EXIT 
   else 
     ibmax = ib
   end if
 end do

 changed = (ibmin /= ibmin_bkp) .or. (ibmax /= ibmax_bkp)

end subroutine enclose_degbands
!!***


!----------------------------------------------------------------------

!!****f* m_ebands/get_nelect_per_spin
!! NAME
!!  get_nelect_per_spin
!!
!! FUNCTION
!!   return number of electrons in each spin channel
!!
!! INPUTS
!!  BSt<ebands_t>=The object describing the band structure.
!!
!! OUTPUT
!!  nelect_per_spin(BSt%nsppol)=For each spin the number of electrons (eventually fractional)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_nelect_per_spin(BSt) result(nelect_per_spin)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_nelect_per_spin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: BSt
!arrays
 real(dp) :: nelect_per_spin(BSt%nsppol)

!Local variables-------------------------------
!scalars
 integer :: iband,ikpt,isppol, ii
!arrays

! *************************************************************************

 nelect_per_spin = BSt%nelect
 if (BSt%nsppol > 1) then
   ii = 0
   nelect_per_spin = zero
   do isppol =1, BSt%nsppol
     do ikpt = 1, BSt%nkpt
       do iband = 1, BSt%nband(ikpt+BSt%nkpt*(isppol-1))
         ii=ii+1
         nelect_per_spin(isppol) = nelect_per_spin(isppol) + BSt%wtk(ikpt)*BSt%occ(iband, ikpt, isppol)
       end do
     end do
   end do
 end if

end function get_nelect_per_spin
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_minmax
!! NAME
!!  get_minmax
!!
!! FUNCTION
!!  Report the min and max value over k-points and bands of (eig|occ|doccde) for each 
!!  spin. Cannot use F90 array syntax due to the internal storage used in abinit.
!!
!! INPUTS
!!  Bands<ebands_t>=The object describing the band structure.
!!  arr_name=The name of the array whose min and Max value has to be calculated.
!!   Possible values: 'occ', 'eig' 'doccde'
!!
!! OUTPUT
!! minmax(2,Bands%nsppol)=For each spin the min and max value of the quantity specified by "arr_name"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_minmax(Bands,arr_name) result(minmax)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_minmax'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),target,intent(in) :: Bands
 character(len=*),intent(in) :: arr_name
!arrays
 real(dp) :: minmax(2,Bands%nsppol)

!Local variables-------------------------------
!scalars
 integer :: band,ikpt,spin,nband_k 
 real(dp) :: datum
 character(len=500) :: msg
!arrays
 real(dp),pointer :: rdata(:,:,:)

! *************************************************************************

 SELECT CASE (tolower(arr_name))
 CASE ('occ')
   rdata => Bands%occ
 CASE ('eig')
   rdata => Bands%eig
 CASE ('doccde')
   rdata => Bands%doccde
 CASE DEFAULT
   write(msg,'(2a)')' Wrong value of arr_name = ',TRIM(arr_name)
   MSG_BUG(msg)
 END SELECT

 minmax(1,:)=greatest_real
 minmax(2,:)=smallest_real
 
 do spin=1,Bands%nsppol
   do ikpt=1,Bands%nkpt
     nband_k=Bands%nband(ikpt+(spin-1)*Bands%nkpt)
     do band=1,nband_k
       datum=rdata(band,ikpt,spin)
       minmax(1,spin)=MIN(minmax(1,spin),datum)
       minmax(2,spin)=MAX(minmax(2,spin),datum)
     end do
   end do
 end do

end function get_minmax
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/bands_edstats
!! NAME
!! bands_edstats
!!
!! FUNCTION
!!  Compute statistical parameters of the energy differences e_ks[b+1] - e_ks[b]
!!  Returns stats_t record with the results (mean, stdev, min, max)
!!
!! INPUTS
!!  Bands<ebands_t>=Band energies.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function bands_edstats(Bands) result(Stats)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bands_edstats'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: Bands
 type(stats_t) :: Stats

!Local variables-------------------------------
!scalars
 integer :: ikibz,nband_k,spin
!arrays
 real(dp),allocatable :: ediffs(:,:,:),edvec(:)

! *************************************************************************

! Compute energy difference between b+1 and b.
 ABI_CALLOC(ediffs, (Bands%mband-1,Bands%nkpt,Bands%nsppol))

 do spin=1,Bands%nsppol
   do ikibz=1,Bands%nkpt
     nband_k=Bands%nband(ikibz+(spin-1)*Bands%nkpt) 
     if (nband_k > 1) then
       ediffs(1:nband_k-1,ikibz,spin) = Bands%eig(2:nband_k,ikibz,spin) - Bands%eig(1:nband_k-1,ikibz,spin)
     end if
   end do
 end do

 ! Calculate the statistical parameters 
 ! Not completely correct if nband_k varies but ...
 ABI_MALLOC(edvec, ((Bands%mband-1)*Bands%nkpt*Bands%nsppol))
 edvec = reshape(ediffs, [(Bands%mband-1)*Bands%nkpt*Bands%nsppol])

 Stats = stats_eval(edvec)

 ABI_FREE(ediffs)
 ABI_FREE(edvec)

end function bands_edstats
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_FS
!! NAME
!!  get_FS
!!
!! FUNCTION
!!  Returns the indeces of the k-points belonging to the Fermi surface as well
!!  as the minimum and Maximum index of the bands crossing the Fermi level.
!!
!! INPUTS
!!  Bands<ebands_t>=The object describing the band structure.
!!  tolweight=Tolerange on the weighs for FS integrations. A k-point belongs
!!    to the Fermi surface if for some band its weight is > tolweight
!!
!! OUTPUT
!!  nFSkpt=Number of points on the Fermi surface.
!!  fs2ibz(1:nFSkpt)=Index of the FS points in the input IBZ array.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_FS(Bands,spin,tolweight,fs_weight,bmin,bmax,nFSkpt,fs2ibz)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_FS'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin
 integer,intent(out) :: bmin,bmax
 real(dp),intent(in) :: tolweight
 type(ebands_t),intent(in) :: Bands
 integer,intent(out) ::  nFSkpt
!arrays
 integer,intent(out) :: fs2ibz(Bands%nkpt)
 real(dp),intent(in) :: fs_weight(Bands%mband,Bands%nkpt)

!Local variables-------------------------------
 integer :: band,ikpt,nband_k,i1,i2
 logical :: seen

! *************************************************************************
 
 nFSkpt=0
 bmin = HUGE(1)
 bmax = 0
 fs2ibz(:)=0

 do ikpt=1,Bands%nkpt
   nband_k=Bands%nband(ikpt+(spin-1)*Bands%nkpt)

   seen=.FALSE.
   do band=1,nband_k
    if (fs_weight(band,ikpt) > tolweight) then 
      if (.not.seen) then
        seen=.TRUE.
        i1=band
        i2=band
      end if
    else 
      i2=band
    end if
   end do

   if (seen) then 
     nFSkpt = nFSkpt+1 
     fs2ibz(nFSkpt) = ikpt
     bmin = MIN(bmin,i1)
     bmax = MAX(bmax,i2)
   end if
 end do !ikpt

end subroutine get_FS
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/bst_metallic_scheme
!! NAME
!! bst_metallic_scheme
!!
!! FUNCTION
!! Returns .TRUE. if metallic occupation scheme is used.
!!
!! INPUTS
!! Bands<ebands_t>=The ebands_t datatype
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function bst_metallic_scheme(Bands)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bst_metallic_scheme'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical :: bst_metallic_scheme
 type(ebands_t),intent(in) :: Bands

! *************************************************************************

 bst_metallic_scheme = ( ANY(Bands%occopt == (/3,4,5,6,7,8/)) )

end function bst_metallic_scheme
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_3dprint
!! NAME
!!  ebands_3dprint
!!
!! FUNCTION
!!   Write 3D energies for Fermi surface visualization (XSF format)
!!
!! INPUTS
!!  Bands<ebands_t>=The object describing the band structure.
!!  Crystal<crystal_t>=Info on unit cell and symmetries.
!!  fname=File name for output.
!!  kptrlatt(3,3)=Matrix partly defining the mesh. See related input variable.
!!  nshiftk=Number of shifts in K-mesh (usually 1)
!!  shiftk(3,nshiftk)=The shifts of the mesh. Xcrysden requires Gamma-centered k-meshes.
!!  [format]=String specifying the format to be used:
!!     1 => BXSF Xcrysden file format [DEFAULT]
!!
!! OUTPUT
!!  ierr=Status error.
!!
!! SIDE EFFECTS
!!  Produce BXSF file.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ebands_3dprint(Bands,Crystal,kptrlatt,nshiftk,shiftk,fname,format) result(ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_3dprint'
 use interfaces_62_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nshiftk
 integer :: ierr
 character(len=*),optional,intent(in) :: format
 character(len=*),intent(in) :: fname
 type(ebands_t),intent(in) :: Bands
 type(crystal_t),intent(in) :: Crystal
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)

!Local variables-------------------------------
 character(len=500) :: myform
 logical :: use_timrev
 character(len=500) :: msg

! *************************************************************************

 ierr=0
 myform="bxsf"; if (PRESENT(format)) myform=format

 use_timrev = (Crystal%timrev==2)

 select case (myform)
 case ("xsf", "bxsf")
   call printbxsf(Bands%eig,zero,Bands%fermie,Crystal%gprimd,kptrlatt,Bands%mband,Bands%nkpt,Bands%kptns,&
&    Crystal%nsym,Crystal%use_antiferro,Crystal%symrec,Crystal%symafm,use_timrev,Bands%nsppol,shiftk,nshiftk,fname,ierr)

 case default
   ierr = ierr+1
   msg = "Unsupported value for format: "//trim(myform)//". Fermi surface file won't be produced."
   MSG_WARNING(msg)
 end select

end function ebands_3dprint
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/update_occ
!! NAME
!! update_occ
!!
!! FUNCTION
!! Calculate new occupation numbers, the Fermi level and the Max occupied band index 
!! for each spin channel starting from the the knowledge of eigenvalues.
!!
!! INPUTS
!!  spinmagntarget=if differ from -99.99d0, fix the spin polarization (in Bohr magneton)
!!  [stmbias]=
!!  [prtvol]=Verbosity level (0 for lowest level)
!!  Bands<ebands_t>=Info on the band structure, the smearing technique and the physical temperature used.
!!
!! OUTPUT
!!  see also SIDE EFFECTS.
!!
!! SIDE EFFECTS
!!  === For metallic occupation the following quantites are recalculated ===
!!   %fermie=the new Fermi energy
!!   %entropy=the new entropy associated with the smearing.
!!   %occ(mband,nkpt,nsppol)=occupation numbers
!!   %doccde(mband,nkpt,nsppol)=derivative of occupancies wrt the energy for each band and k point
!!  === In case of semiconductors ===
!!   All the quantitities in Bands are left unchanged with the exception of:
!!   %fermie=Redefined so that it is in the middle of the gap
!!   %entropy=Set to zero
!!
!! PARENTS
!!      bethe_salpeter,elphon,get_nv_fs_temp,get_tau_k,m_ebands,optic,screening
!!      setup_bse,setup_bse_interp,setup_sigma,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine update_occ(Bands,spinmagntarget,stmbias,prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'update_occ'
 use interfaces_14_hidewrite
 use interfaces_62_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(inout) :: Bands
 integer,optional,intent(in) :: prtvol
 real(dp),intent(in) :: spinmagntarget
 real(dp),optional,intent(in) :: stmbias
!arrays

!Local variables-------------------------------
!scalars
 integer :: band,mband,ikibz,nkpt,spin,nsppol,my_prtvol,nband_k
 real(dp) :: entropy,fermie,stmbias_local,ndiff,cbot,vtop,maxocc
 character(len=500) :: msg
!arrays
 real(dp) :: nelect_spin(Bands%nsppol),condbottom(Bands%nsppol),valencetop(Bands%nsppol)
 real(dp),allocatable :: doccdet(:),occt(:),eigent(:)

! *************************************************************************

 my_prtvol    =0   ; if (PRESENT(prtvol )) my_prtvol    =prtvol
 stmbias_local=zero; if (PRESENT(stmbias)) stmbias_local=stmbias

 if (Bands%occopt>=3.and.Bands%occopt<=8) then 
   !  If occupation is metallic have to compute new occupation numbers.
   write(msg,'(a,f9.5)')' metallic scheme, calling newocc with spinmagntarget = ',spinmagntarget
   call wrtout(std_out,msg,'COLL')

   mband  = Bands%mband ! to make the interface of newocc happy.
   nkpt   = Bands%nkpt
   nsppol = Bands%nsppol

   ABI_MALLOC(eigent,(mband*nkpt*nsppol))
   call get_eneocc_vect(Bands,'eig',eigent)

   ABI_MALLOC(occt,(mband*nkpt*nsppol))
   ABI_MALLOC(doccdet,(mband*nkpt*nsppol))

   call newocc(doccdet,eigent,entropy,fermie,spinmagntarget,mband,Bands%nband,&
&    Bands%nelect,Bands%nkpt,Bands%nspinor,Bands%nsppol,occt,Bands%occopt,&
&    my_prtvol,stmbias_local,Bands%tphysel,Bands%tsmear,Bands%wtk)
   !
   ! Save output in Bands%. 
   Bands%entropy = entropy
   Bands%fermie  = fermie
   call put_eneocc_vect(Bands,'occ'   ,occt   ) 
   call put_eneocc_vect(Bands,'doccde',doccdet) 
   ABI_FREE(eigent)
   ABI_FREE(occt)
   ABI_FREE(doccdet)

 else  
   !  Semiconductor or Insulator.
   ! 
   ! FIXME here there is an inconsistency btw GW and Abinit 
   ! In abinit Fermi is set to HOMO while in GW fermi is in the middle
   ! of Gap. In case of crystal systems, the later convention should be preferable.
   ! Anyway we have to decide and follow a unique convention to avoid problems.
   !
   ! occupation factors MUST be initialized
   if (ALL(ABS(Bands%occ) < tol6)) then
     msg = "occupation factors are not initialized, likely due to the use of iscf=-2"
     MSG_ERROR(msg)
   end if

   maxocc=two/(Bands%nsppol*Bands%nspinor)

   ! * Calculate the valence index for each spin channel.
   do spin=1,Bands%nsppol
     valencetop(spin)= smallest_real
     condbottom(spin)= greatest_real

     do ikibz=1,Bands%nkpt
       nband_k=Bands%nband(ikibz+(spin-1)*Bands%nkpt) 
       do band=1,nband_k
         if (Bands%occ(band,ikibz,spin)/maxocc>one-tol6 .and. valencetop(spin)<Bands%eig(band,ikibz,spin)) then 
           valencetop(spin)=Bands%eig(band,ikibz,spin)
         end if
         if (Bands%occ(band,ikibz,spin)/maxocc<tol6 .and. condbottom(spin)>Bands%eig(band,ikibz,spin)) then 
           condbottom(spin)=Bands%eig(band,ikibz,spin)
         end if
       end do
     end do 

   end do 

   vtop=MAXVAL(valencetop)
   cbot=MINVAL(condbottom)
   write(msg,'(a,f6.2,2a,f6.2)')&
&    ' top of valence       [eV] ',vtop*Ha_eV,ch10,&
&    ' bottom of conduction [eV] ',cbot*Ha_eV
   call wrtout(std_out,msg,'COLL')
   if (Bands%nsppol==2) then 
     if (ABS(vtop-MINVAL(valencetop))>tol6) then 
       write(msg,'(a,i2)')' top of valence is spin ',MAXLOC(valencetop)
       call wrtout(std_out,msg,'COLL')
     end if
     if (ABS(cbot-MAXVAL(condbottom))>tol6) then 
       write(msg,'(a,i2)')' bottom of conduction is spin ',MINLOC(condbottom)
       call wrtout(std_out,msg,'COLL')
     end if
   end if

   ! === Save output === 
   ! Here I dont know if it is better to be consistent with the abinit convention i.e fermi=vtop
   Bands%entropy=zero
   Bands%fermie=(vtop+cbot)/2 
   if (ABS(cbot-vtop)<1.d-4) Bands%fermie=vtop ! To avoid error on the last digit FIXME is it really needed
 end if

 write(msg,'(a,f6.2,a)')' Fermi energy         [eV] ',Bands%fermie*Ha_eV,ch10
 call wrtout(std_out,msg,'COLL')
 !
 ! === Compute number of electrons for each spin channel ===
 nelect_spin(:)=zero 
 do spin=1,Bands%nsppol
   do ikibz=1,Bands%nkpt
     nband_k=Bands%nband(ikibz+(spin-1)*Bands%nkpt)
     nelect_spin(spin)= nelect_spin(spin) + Bands%wtk(ikibz)*SUM(Bands%occ(1:nband_k,ikibz,spin))
   end do
 end do

 ndiff=Bands%nelect-SUM(nelect_spin)
 if (my_prtvol>0) then
   write(msg,'(2a,f6.2,2a,f7.4)')ch10,&
&    ' total number of electrons = ',SUM(nelect_spin),ch10,&
&    ' input and calculated no. of electrons differ by ',ndiff 
   call wrtout(std_out,msg,'COLL')
 end if

 if (ABS(ndiff)>5.d-2*Bands%nelect) then
   write(msg,'(2a,2(a,f6.2))')&
&    'Too large difference in no. of electrons:,',ch10,&
&    'Expected= ',Bands%nelect,' Calculated= ',SUM(nelect_spin)
   MSG_ERROR(msg)
 end if

end subroutine update_occ
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ReportGap
!! NAME
!! ReportGap
!!
!! FUNCTION
!!  Print info on the fundamental and optical gap.
!!
!! INPUTS
!!  Bands<ebands_t>=Info on the band structure, the smearing technique and the physical temperature used.
!!  [header]=Optional title.
!!  [kmask]=Logical mask used to exclude k-points.
!!  [unit]=Optional unit for output (std_out if not specified)
!!  [mode_paral]=Either "COLL" or "PERS", former is default.
!!
!! OUTPUT
!!  writing.
!!  [gaps(3,nsppol)]=Fundamental and optical gaps. The third index corresponds to a "status" : 0.0dp if gaps were not computed
!!     (because there are only valence bands) ; -1.0dp if the system (or spin-channel) is metallic ; 1.0dp if the gap was computed
!!
!! PARENTS
!!      gstate,m_exc_diago,setup_bse,setup_bse_interp,setup_sigma,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ReportGap(Bands,header,kmask,unit,mode_paral,gaps)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ReportGap'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 type(ebands_t),intent(in)  :: Bands
!arrays
 real(dp),optional,intent(out) :: gaps(3,Bands%nsppol)
 logical,optional,intent(in) ::  kmask(Bands%nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikibz,nband_k,spin,nsppol,ikopt,ivk,ick,ivb,icb,my_unt,first
 real(dp),parameter :: tol_fermi=tol6
 real(dp) :: fun_gap,opt_gap
 logical :: ismetal
 character(len=4) :: my_mode
 character(len=500) :: msg
!arrays
 integer :: val_idx(Bands%nkpt,Bands%nsppol)
 real(dp) :: top_valence(Bands%nkpt),bot_conduct(Bands%nkpt) 
 logical :: my_kmask(Bands%nkpt)

! *********************************************************************

 nsppol = Bands%nsppol

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral
 my_kmask=.TRUE.; if (PRESENT(kmask     )) my_kmask=kmask

 if (PRESENT(gaps)) gaps=zero

 val_idx(:,:) = get_valence_idx(Bands,tol_fermi)
 first=0

!Initialize the return status for the gaps
 if (PRESENT(gaps)) gaps(1:3,1:nsppol)=zero

 do spin=1,nsppol

   ! No output if system i metallic
   ismetal=ANY(val_idx(:,spin)/=val_idx(1,spin)) 
   if (ismetal) then
     if (PRESENT(gaps)) gaps(3,nsppol)=-one
     CYCLE 
   endif

   first=first+1
   if (first==1) then
     msg=ch10
     if (PRESENT(header)) msg=ch10//' === '//TRIM(ADJUSTL(header))//' === '
     call wrtout(my_unt,msg,my_mode) 
   end if

   ivb=val_idx(1,spin)
   icb=ivb+1

   do ikibz=1,Bands%nkpt
     if (.not.my_kmask(ikibz)) CYCLE
     nband_k=Bands%nband(ikibz+(spin-1)*Bands%nkpt)
     top_valence(ikibz)=Bands%eig(ivb,ikibz,spin)
     if (icb>nband_k) then
       GOTO 10 ! Only occupied states are present, no output!
     endif
     bot_conduct(ikibz)=Bands%eig(icb,ikibz,spin)
   end do

   ! === Get minimum of the optical Gap ===
   ikopt= imin_loc(bot_conduct-top_valence,MASK=my_kmask)
   opt_gap=bot_conduct(ikopt)-top_valence(ikopt)

   ! === Get fundamental Gap ===
   ick = imin_loc(bot_conduct,MASK=my_kmask)
   ivk = imax_loc(top_valence,MASK=my_kmask)
   fun_gap = Bands%eig(icb,ick,spin)-Bands%eig(ivb,ivk,spin)

   write(msg,'(a,i2,a,2(a,f8.4,a,3f8.4,a),33x,a,3f8.4)')&
&    '  >>>> For spin ',spin,ch10,&
&    '   Minimum optical gap = ',opt_gap*Ha_eV,' [eV], located at k-point      : ',Bands%kptns(:,ikopt),ch10,&
&    '   Fundamental gap     = ',fun_gap*Ha_eV,' [eV], Top of valence bands at : ',Bands%kptns(:,ivk),ch10,  &
&                                              '       Bottom of conduction at : ',Bands%kptns(:,ick)
   call wrtout(my_unt,msg,my_mode) 

   if (PRESENT(gaps)) then
     gaps(:,spin) = (/fun_gap,opt_gap,one/)
   end if

 end do !spin

 return

 10 continue
 MSG_WARNING("Not enough states to calculate the band gap.")

end subroutine ReportGap
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ExpandBands
!! NAME
!! ExpandBands
!!
!! FUNCTION
!!  Return a new object of type ebands_t corresponding to a list of k-points 
!!  specified in input. Symmetry properties of the eigenvectors are used to 
!!  symmetrize energies and occupation numbers.
!!
!! INPUTS
!!  Bands<ebands_t>=Info on the band structure 
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  In info/=0, some of the k-points in klist have no corresponding image in Bands_in%kptns
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ExpandBands(Bands_in,nklist,klist,use_tr,use_afm,nsym,symrec,symafm,info) result(Bands_out)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ExpandBands'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym,nklist
 integer,intent(out) :: info
 logical,intent(in) :: use_tr,use_afm
 type(ebands_t),intent(in) :: Bands_in 
 type(ebands_t) :: Bands_out
!arrays
 real(dp),intent(in) :: klist(3,nklist)
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(in) :: symafm(nsym)

!Local variables-------------------------------
!scalars
 integer :: iklist,itim,timrev,isym,ikibz,nkfound,ii,jj,nband_k,spin,ikdx
 character(len=500) :: msg
 logical :: found
!arrays
 integer :: G0(3)
 real(dp) :: krot(3)
 integer,allocatable :: klist2ibz(:)
! *********************************************************************

 do ii=1,nklist-1
   do jj=ii+1,nklist
     if (isamek(klist(:,ii),klist(:,jj),G0)) then
       write(msg,'(a,i4,a,i4,5a)')&
&       ' Points ',ii,' and ',jj,' in klist are equivalent ',ch10,&
&       ' This is not allowed in the present implementation.',ch10,&
&       ' Change the input values to avoid duplicated k-points. '
       MSG_ERROR(msg)
     end if
   end do
 end do

 info=0
 timrev=1 ; if (use_tr) timrev=2
 ABI_MALLOC(klist2ibz,(nklist))
 klist2ibz=0

 do iklist=1,nklist
   found=.FALSE.

ibzl: do ikibz=1,Bands_in%nkpt
        do itim=1,timrev
          do isym=1,nsym
            if (use_afm.and.symafm(isym)==-1) CYCLE
            ! * Form IS k ===
            krot(:)=(3-2*itim)*MATMUL(symrec(:,:,isym),Bands_in%kptns(:,ikibz))

            ! * Check whether it is equal to klist(:,ilist) within a RL vector.
            ! FIXME see notes below related to this function
            if (isamek(krot,klist(:,iklist),G0)) then 
              found=.TRUE.
              klist2ibz(iklist)=ikibz
              EXIT ibzl
            end if

          end do !isym 
        end do !itim
      end do ibzl

      if (.not.found) info=info+1
 end do !iklist

 nkfound=COUNT(klist2ibz/=0)

 Bands_out%nkpt    = nkfound
 Bands_out%nspinor = Bands_in%nspinor      
 Bands_out%nsppol  = Bands_in%nsppol
 Bands_out%occopt  = Bands_in%occopt        

 Bands_out%entropy = Bands_in%entropy 
 Bands_out%fermie  = Bands_in%fermie  
 Bands_out%nelect  = Bands_in%nelect  
 Bands_out%tphysel = Bands_in%tphysel 
 Bands_out%tsmear  = Bands_in%tsmear  

 ABI_MALLOC(Bands_out%nband,(Bands_out%nkpt*Bands_out%nsppol))
 ikdx=0
 do spin=1,Bands_out%nsppol
   do iklist=1,nklist
     ikibz=klist2ibz(iklist)
     if (ikibz==0) CYCLE
     ikdx=ikdx+1
     nband_k=Bands_in%nband(ikibz+(spin-1)*Bands_in%nkpt)
     Bands_out%nband(ikdx)=nband_k
   end do
 end do

 Bands_out%mband =MAXVAL(Bands_out%nband)
 Bands_out%bantot=SUM   (Bands_out%nband)

 ABI_MALLOC(Bands_out%istwfk,(Bands_out%nkpt))
 ABI_MALLOC(Bands_out%nband,(Bands_out%nkpt*Bands_out%nsppol))
 ABI_MALLOC(Bands_out%npwarr,(Bands_out%nkpt))
 ABI_MALLOC(Bands_out%kptns,(3,Bands_out%nkpt))
 ABI_MALLOC(Bands_out%eig   ,(Bands_out%mband,Bands_out%nkpt,Bands_out%nsppol))
 ABI_MALLOC(Bands_out%occ   ,(Bands_out%mband,Bands_out%nkpt,Bands_out%nsppol))
 ABI_MALLOC(Bands_out%doccde,(Bands_out%mband,Bands_out%nkpt,Bands_out%nsppol))
 ABI_MALLOC(Bands_out%wtk,(Bands_out%nkpt))
 !
 ! === Copy arrays that depend only on k ===
 ikdx=0
 do iklist=1,nklist
   ikibz=klist2ibz(iklist)
   if (ikibz==0) CYCLE
   ikdx=ikdx+1
   Bands_out%istwfk(ikdx) =Bands_in%istwfk(ikibz)      
   Bands_out%npwarr(ikdx) =Bands_in%npwarr(ikibz)      
   Bands_out%kptns(:,ikdx)=klist(:,iklist) !Use klist, not Bands_in%kptns
   Bands_out%wtk(ikdx)    =Bands_in%wtk(ikibz)        
 end do
 !
 ! * Renormalize weights
 Bands_out%wtk=Bands_out%wtk/SUM(Bands_out%wtk)

 ! === Copy arrays that depend on b-k-s ===
 do spin=1,Bands_in%nsppol
   ikdx=0
   do iklist=1,nklist
     ikibz=klist2ibz(iklist)
     if (ikibz==0) CYCLE
     ikdx=ikdx+1
     Bands_out%eig   (:,ikdx,spin)=Bands_in%eig   (:,ikibz,spin)        
     Bands_out%occ   (:,ikdx,spin)=Bands_in%occ   (:,ikibz,spin)        
     Bands_out%doccde(:,ikdx,spin)=Bands_in%doccde(:,ikibz,spin)     
   end do !ikpt
 end do !spin

 ABI_FREE(klist2ibz)

end function ExpandBands
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_ncwrite
!! NAME
!! ebands_ncwrite
!!
!! FUNCTION
!!  Writes the content of a ebands_t object to a NETCDF file 
!!  according to the ETSF-IO specifications.
!!
!! INPUTS
!!  nshiftk_orig=number of original shifts (0 means that the mesh was specified with kptrlatt)
!!  shiftk_orig(3,8)=Original shifts given in input.
!!  nshiftk=Number of shiftks computed in inkpts.
!!  shiftk=Shifts computed in inkpts.F90
!!  ngkpt(3)=Number of divisions for MP sampling (0 means that the mesh was specified with kptrlatt).
!!  kptrlatt(3,3)=kptrlatt input variable
!!  ncid =NC file handle
!!
!! OUTPUT
!!
!! PARENTS
!!      eig2tot,loper3,m_shirley,outscfcv,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_ncwrite(Bands,nshiftk_orig,shiftk_orig,nshiftk,shiftk,ngkpt,kptrlatt,ncid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid 
 integer,target,intent(in) :: nshiftk,nshiftk_orig
 type(ebands_t),target,intent(in) :: Bands
!arrays
 integer,target,intent(in) :: ngkpt(3),kptrlatt(3,3)
 real(dp),target,intent(in) :: shiftk(3,nshiftk),shiftk_orig(3,nshiftk_orig)

!Local variables-------------------------------
!scalars
 integer,target :: nelect_int
 logical :: lstat,write_kptrlatt,write_ngkpt
#ifdef HAVE_TRIO_ETSF_IO
 type(ETSF_io_low_error) :: Error_data
 type(ETSF_electrons),target :: Electrons
 type(ETSF_kpoints),target :: Kpoints
#endif

! *************************************************************************

#ifdef HAVE_TRIO_ETSF_IO
 ! ==============================================
 ! === Write the dimensions specified by ETSF ===
 ! ==============================================
 !FIXME: do not handle k_dependent = 1
 call etsf_io_low_set_define_mode(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'max_number_of_states',Bands%mband,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_spinor_components',Bands%nspinor,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_spins',Bands%nsppol,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_kpoints',Bands%nkpt,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! Unofficial variables. Notes:
 ! 1) ETSF-IO does not support nshifts > 1 FIXME
 ! 2) shiftk_orig, nshiftk_orig refers to the values specified in the input (most useful ones).
 ! 3) shiftk, kptrlatt refers to the values computed in inkpts.
 write_kptrlatt = (SUM(ABS(kptrlatt))/=0)
 write_ngkpt = (ANY(ngkpt/=0))

 if (write_ngkpt) then
   call etsf_io_low_write_dim(ncid,'ngkpt_nshiftk',nshiftk_orig,lstat,Error_data=Error_data)  
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if

 if (write_kptrlatt) then
   ABI_CHECK(nshiftk==1,"nshiftk!=1")
   call etsf_io_low_write_dim(ncid,'nshiftk',nshiftk,lstat,Error_data=Error_data)  
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if

 call etsf_io_kpoints_def(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_electrons_def(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 Electrons%fermi_energy            => Bands%fermie
 nelect_int=NINT(Bands%nelect)
 Electrons%number_of_electrons     => nelect_int  !FIXME in ETSF this is integer
 Electrons%smearing_width          => Bands%tsmear 
 Electrons%number_of_states%data1D => Bands%nband 

 Electrons%eigenvalues%data3D => Bands%eig
 Electrons%occupations%data3D => Bands%occ

 Kpoints%reduced_coordinates_of_kpoints => Bands%kptns
 Kpoints%kpoint_weights                 => Bands%wtk

 if (write_ngkpt) then
   Kpoints%monkhorst_pack_folding => ngkpt
 end if

 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_electrons_put(ncid, Electrons, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_kpoints_put(ncid, Kpoints, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! ===========================================================
 ! === Write abinit-related stuff (not covered by ETSF-IO) ===
 ! ===========================================================
 ! Define variables.
 call etsf_io_low_set_define_mode(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'tphysel',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'occopt',etsf_io_low_integer,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'istwfk',etsf_io_low_integer,(/'number_of_kpoints'/),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! Abinit variables defining the K-point sampling.
 if (write_kptrlatt) then
   call etsf_io_low_def_var(ncid,'shiftk',etsf_io_low_double,(/pad('number_of_reduced_dimensions'), pad('nshiftk')/),&
&    lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)

   call etsf_io_low_def_var(ncid,'kptrlatt',etsf_io_low_integer,&
&   (/pad('number_of_reduced_dimensions'), pad('number_of_reduced_dimensions')/), lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if

 if (write_ngkpt) then
   call etsf_io_low_def_var(ncid,'ngkpt_shiftk',etsf_io_low_double,&
&    (/pad('number_of_reduced_dimensions'), pad('ngkpt_nshiftk')/),lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if

 ! Write variables
 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'tphysel',Bands%tphysel,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'occopt',Bands%occopt,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'istwfk',Bands%istwfk,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 if (write_kptrlatt) then
   call etsf_io_low_write_var(ncid,'shiftk',shiftk,lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
   call etsf_io_low_write_var(ncid,'kptrlatt',kptrlatt,lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if

 if (write_ngkpt) then
   !write(std_out,*)"nshiftk_orig",nshiftk_orig
   !write(std_out,*)"shiftk_orig",shiftk_orig
   call etsf_io_low_write_var(ncid,'ngkpt_shiftk',shiftk_orig,lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if

#else 
 MSG_ERROR("ETSF-IO support is not activated. ")
#endif

end subroutine ebands_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_ncread
!! NAME
!! ebands_ncread
!!
!! FUNCTION
!!  Initialize an instance of the ebands_t from a NETCDF file
!!  written according to the ETSF-IO specifications.
!!
!! INPUTS
!!  fname=Name of the netcdf file.
!!
!! OUTPUT
!!  Bands<type(ebands_t)>=Band energies and occupation factors.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_ncread(Bands,fname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_ncread'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: fname
 type(ebands_t),target,intent(out) :: Bands
!arrays

!Local variables-------------------------------
!scalars
 integer :: ncid,fform,usewvl,formeig,idx,spin,nkpt,band,ikpt
 integer,target :: nelect_int
 real(dp) :: spinmagntarget_
 logical :: lstat
#ifdef HAVE_TRIO_ETSF_IO
 type(ETSF_dims) :: Dims
 type(ETSF_io_low_error) :: Error_data
 type(ETSF_basisdata) :: Basisdata
 type(ETSF_groups) :: GroupFolder
 type(ETSF_electrons),target :: Electrons
 type(ETSF_kpoints),target :: Kpoints
#endif
 type(Hdr_type) :: Hdr
!arrays
 real(dp),allocatable,target :: eig_vec(:),occ_vec(:)

! *************************************************************************

 DBG_ENTER("COLL")

#ifdef HAVE_TRIO_ETSF_IO
 ! === Open the file ===
 call wrtout(std_out,' ebands_ncread : about to read file '//TRIM(fname),'COLL')

 call etsf_io_low_open_read(ncid,fname,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! === Read dimensions handled by ETSF ===
 call etsf_io_dims_get(ncid,Dims,lstat,Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! TODO : do not handle k_dependent = 1
 Bands%bantot   = Dims%max_number_of_states * Dims%number_of_kpoints * Dims%number_of_spins
 Bands%mband    = Dims%max_number_of_states
 Bands%nkpt     = Dims%number_of_kpoints
 Bands%nspinor  = Dims%number_of_spinor_components
 Bands%nsppol   = Dims%number_of_spins
 !? Bands%nspden   = Dims%number_of_components

 ABI_MALLOC(Bands%istwfk,(Bands%nkpt))
 ABI_MALLOC(Bands%nband,(Bands%nkpt*Bands%nsppol))
 ABI_MALLOC(Bands%npwarr,(Bands%nkpt))
 ABI_MALLOC(Bands%kptns,(3,Bands%nkpt))
 ABI_MALLOC(Bands%eig,(Bands%mband,Bands%nkpt,Bands%nsppol))
 ABI_MALLOC(Bands%occ,(Bands%mband,Bands%nkpt,Bands%nsppol))
 ABI_MALLOC(Bands%doccde,(Bands%mband,Bands%nkpt,Bands%nsppol))
 ABI_MALLOC(Bands%wtk,(Bands%nkpt))

 usewvl=0
 call etsf_io_low_read_var(ncid,'usewvl',usewvl,lstat,error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! Get all variables included in ETSF
 if (usewvl==0) then
   BasisData%number_of_coefficients => Bands%npwarr
   call etsf_io_basisdata_get(ncid,Basisdata,lstat,Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if

 Electrons%fermi_energy            => Bands%fermie
 !FIXME this is integer, therefore alchemy or charged cells won"t work!
 Electrons%number_of_electrons     => nelect_int 
 Electrons%smearing_width          => Bands%tsmear 
 Electrons%number_of_states%data1D => Bands%nband 

 !TODO DFPT not treated, we read correctly but the object should be modified a bit. 
 formeig=0
 ABI_MALLOC(eig_vec,((2*Bands%mband)**formeig*Bands%mband*Bands%nkpt*Bands%nsppol))
 ABI_MALLOC(occ_vec,(Bands%mband*Bands%nkpt*Bands%nsppol))
 Electrons%eigenvalues%data1D      => eig_vec  !then we have to unpack
 Electrons%occupations%data1D      => occ_vec

 Kpoints%reduced_coordinates_of_kpoints => Bands%kptns
 Kpoints%kpoint_weights                 => Bands%wtk

 GroupFolder%Electrons => Electrons
 GroupFolder%Kpoints   => Kpoints

 call etsf_io_data_read(fname,GroupFolder,lstat,Error_data)

 idx=0 ! TODO call my helper function once modules will be supported.
 do spin=1,Bands%nsppol
   do ikpt=1,Bands%nkpt
     do band=1,Bands%nband(ikpt+(spin-1)*nkpt)
       ! MG020313 this write is needed to avoid sigfaults with open64, psc and ibm
       write(std_out,*)spin, ikpt, band
       idx=idx+1
       Bands%eig(band,ikpt,spin)=eig_vec(idx)
       Bands%occ(band,ikpt,spin)=occ_vec(idx)
     end do
   end do
 end do
 ABI_FREE(eig_vec)
 ABI_FREE(occ_vec)

 ! === Read the abinit header ===
 ! * Fill missing quantities in Bands using Hdr.
 call hdr_io_etsf(fform,Hdr,1,ncid)

 Bands%tphysel = Hdr%tphysel
 Bands%istwfk  = Hdr%istwfk
 Bands%occopt  = Hdr%occopt
 call hdr_free(Hdr)

 ! === Close the file ===
 call etsf_io_low_close(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! === Finalize the object ===
 !FIXME At the moment this has to be done outside the routine
 ! Moreover I have to solve the problem with the convention for the Fermi level
 spinmagntarget_=99.99_dp 
 !call update_occ(Bands,spinmagntarget_)

#else 
 MSG_ERROR("ETSF-IO support is not activated.")
#endif
 
 DBG_EXIT("COLL")

end subroutine ebands_ncread
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/edos_init
!! NAME
!!  edos_init
!!
!! FUNCTION
!!  Calculate the electronic density of states.
!!
!! INPUTS
!!  Bands<ebands_t>=The object describing the band structure.
!!  Kmesh<kmesh_t>=Strcuture defining the BZ sampling.
!!  method="gaussian" or "tetra"
!!  step=Step on the linear mesh in Ha. If <0, the routine will use the 
!!    mean of the energy level spacing
!!  broad=Gaussian broadening, If <0, the routine will use a default
!!    value for the broadening computed from the mean of the energy level spacing. 
!!    No meaning if method == "tetra"
!!
!! OUTPUT 
!!  ierr=Exit status. The routine can return ierr/=0 if method=="tetra"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine edos_init(Edos,Bands,Kmesh,method,step,broad,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'edos_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ierr
 real(dp),intent(in) :: step,broad
 character(len=*),intent(in) :: method
 type(ebands_t),target,intent(in)  :: Bands
 type(kmesh_t),target,intent(in)  :: Kmesh
 type(edos_t),intent(out) :: Edos

!Local variables-------------------------------
!scalars
 integer :: iw,nw,spin,band,ikpt,ief
 real(dp) :: max_ene,min_ene,wtk,max_occ,rcvol
 logical :: ltest
 character(len=500) :: msg
 type(stats_t) :: Ediffs
!arrays
 integer :: g0(3)
 real(dp) :: eminmax_spin(2,Bands%nsppol),gprimd(3,3)
 real(dp),allocatable :: wme0(:)
 real(dp),allocatable :: dtweightde(:,:),tmp_eigen(:),tweight(:,:)

! *********************************************************************

 ierr = 0

 ! Keep a reference to Bands and Kmesh
 Edos%Bands => Bands
 Edos%Kmesh => Kmesh

 Edos%nsppol = Bands%nsppol
 Edos%method = method

 ! Compute the mean value of the energy spacing.
 Ediffs = bands_edstats(Bands)
 Edos%broad = broad; if (broad <= tol16) Edos%broad = three * Ediffs%mean
 Edos%step = step; if (step <= tol16) Edos%step =  Ediffs%mean

 ! Compute the linear mesh so that it encloses all bands.
 eminmax_spin = get_minmax(Bands, "eig")
 min_ene = minval(eminmax_spin(1,:)); min_ene = min_ene - 0.1_dp * abs(min_ene)
 max_ene = maxval(eminmax_spin(2,:)); max_ene = max_ene + 0.1_dp * abs(max_ene)

 nw = nint((max_ene - min_ene)/Edos%step) + 1
 Edos%nw = nw

 ABI_MALLOC(Edos%mesh, (nw))
 Edos%mesh = arth(min_ene, Edos%step, nw)

 ABI_CALLOC(Edos%dos,  (nw, 0:Edos%nsppol))
 ABI_CALLOC(Edos%idos, (nw, 0:Edos%nsppol))
 max_occ=two/(Bands%nspinor*Bands%nsppol)  

 select case (method)
 case ("gaussian")
   ABI_MALLOC(wme0, (nw))

   do spin=1,Edos%nsppol
     do ikpt=1,Bands%nkpt
       wtk = Kmesh%wt(ikpt)
       do band=1,Bands%nband(ikpt+(spin-1)*Bands%nkpt)
          wme0 = Edos%mesh - Bands%eig(band, ikpt, spin)
          Edos%dos(:, spin) = Edos%dos(:, spin) + wtk * gaussian(wme0, Edos%broad)
       end do
     end do
   end do

   ABI_FREE(wme0)

 case ("tetra")
   ! Consisteny test: return if something goes wrong.
   if (Bands%nkpt<2) then
     MSG_WARNING('At least 2 points are needed for tetrahedrons')
     ierr = ierr + 1
   end if

   if ( ANY(Bands%nband/=Bands%nband(1)) ) then
     MSG_WARNING('For tetrahedrons, nband(:) must be constant')
     ierr = ierr + 1
   end if

   ! kpoints must be the same and with the same ordering.
   if (Bands%nkpt/=Kmesh%nibz) then
     ierr = ierr + 1
     MSG_WARNING('Mismatch in number of k-points')
   end if

   ltest=.TRUE.
   do ikpt=1,Bands%nkpt
     ltest = (ltest.and.isamek(Bands%kptns(:,ikpt),Kmesh%ibz(:,ikpt),g0))
   end do
   if (.not.ltest) then
     MSG_WARNING('k-points in Kmesh and Bands are not equivalent')
     ierr = ierr + 1
   end if

   if (Kmesh%nshift>1) then 
     write(msg,'(a,i0)')" For tetrahedrons, nshift must be (0,1) but found: ",Kmesh%nshift
     MSG_WARNING(msg)
     ierr = ierr + 1
   end if

   if (Kmesh%has_tetra/=1) then
     MSG_WARNING("Tetrahedrons are not available. call Kmesh%init_tetra_tabs first")
     ierr = ierr + 1
   end if

   if (ierr/=0) then 
     return
   end if

   gprimd = Kmesh%gprimd
   rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
&   -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
&   +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

   ABI_MALLOC(tweight,(Bands%nkpt,nw))
   ABI_MALLOC(dtweightde,(Bands%nkpt,nw))

   ! For each spin and band, interpolate over kpoints, 
   ! calculate integration weights and DOS contribution.
   ABI_MALLOC(tmp_eigen,(Bands%nkpt))
   do spin=1,Bands%nsppol
     do band=1,Bands%nband(1)
       ! For each band get its contribution
       tmp_eigen(:) = Bands%eig(band,:,spin)

       ! Calculate integration weights at each irred kpoint (Blochl et al PRB 49 16223)
       call get_tetra_weight(tmp_eigen,min_ene,max_ene,one,nw,Bands%nkpt,Kmesh%Tetra,tweight,dtweightde)
   
       do iw=1,nw
         do ikpt=1,Bands%nkpt
           Edos%dos(iw,spin) = Edos%dos(iw,spin) + dtweightde(ikpt,iw) 
           ! IDOS is computed afterwards with simpson
           !Edos%idos(iw,spin) = Edos%idos(iw,spin) + tweight(ikpt,iw) 
         end do
       end do
     end do ! band
   end do !spin

   ! Free memory
   ABI_FREE(tmp_eigen)
   ABI_FREE(tweight)
   ABI_FREE(dtweightde)

   ! Filter so that dos[i] is always >= 0 and idos is monotonic
   ! IDOS is computed afterwards with simpson
   where (Edos%dos(:,1:) <= zero) Edos%dos(:,1:) = zero

 case default
   MSG_ERROR("Wrong method: "//trim(method))
 end select

 ! Compute total DOS and IDOS
 Edos%dos(:, 0) = max_occ * sum(Edos%dos(:,1:), dim=2)

 do spin=1,Edos%nsppol
   call simpson_int(nw,Edos%step,Edos%dos(:,spin),Edos%idos(:,spin))
 end do
 Edos%idos(:, 0) = max_occ * sum(Edos%idos(:,1:), dim=2)

 ! Use bisection to find fermi level.
 ! Warning: this code assumes idos[i+1] >= idos[i]. This condition may not be
 ! fullfilled if we use tetra and this is the reason why we have filtered the DOS.
 ief = bisect(Edos%idos(:,0), Bands%nelect) 

 ! Handle out of range condition.
 if (ief ==0 .or. ief == nw) then 
   write(msg,"(3a)")&
&    "Bisection could not find an initial guess for the Fermi level!",ch10,&
&    "Possible reasons: not enough bands or wrong number of electrons"
   MSG_WARNING(msg)
 end if

 ! TODO: Use Linear interpolation to find an improved estimate of the Fermi level?
 Edos%ief = ief

end subroutine edos_init
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/edos_free
!! NAME
!!  edos_free 
!!
!! FUNCTION
!!  Free the memory allocated in edos_t
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine edos_free(edos)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'edos_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(edos_t),intent(inout) :: edos

! *********************************************************************

 !@edos_t
!real
 if (allocated(edos%mesh)) then
   ABI_FREE(edos%mesh)
 end if

 if (allocated(edos%dos)) then
   ABI_FREE(edos%dos)
 end if

 if (allocated(edos%idos)) then
   ABI_FREE(edos%idos)
 end if

! Nullify pointers.
 nullify(Edos%Bands)
 nullify(Edos%Kmesh)

end subroutine edos_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/edos_print
!! NAME
!! edos_print
!!
!! FUNCTION
!! Print out electron DOS to a Fortran unit file open in formatted mode.
!!
!! INPUTS
!!  Edos<edos_t>=DOS container
!!  unit=Fortran logical unit open in formatted mode
!!  [write_data]=.False. if only generic info (DOS(ef) are wanted). Default: .True.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine edos_print(Edos, unit, write_data)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'edos_print'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unit
 logical,optional,intent(in) :: write_data
 type(edos_t),intent(in) :: Edos

!Local variables-------------------------------
 integer :: iw,spin
 real(dp) :: cfact

! *************************************************************************

 ! Convert everything into eV 
 ! I know that Abinit should use Ha but Hartree are not readable.
 ! Plaas don't change this code, in case add an optional argument to specify different units.
 cfact=Ha_eV

 ! Write header.
 write(unit,'(3a)')'# Electron density of states: Energy in eV, DOS in states/eV'

 select case (Edos%method)
 case ("gaussian")
   write(unit,'(a,es16.8,a,i0)')'# Gaussian method with smearing = ',Edos%broad*cfact,'[eV], nkibz = ',Edos%Kmesh%nibz

 case ("tetra") 
   write(unit,'(a,i0)')'# Tetrahedron method, nkibz = ',Edos%Kmesh%nibz

 case default
   MSG_ERROR("Wrong method: "//trim(Edos%method))
 end select

 if (Edos%ief==0) then
   write(unit,'(a)')'# Fermi level: None'
 else
   write(unit,'(a,es16.8)')'# Fermi level: ',Edos%mesh(Edos%ief)*cfact
 end if

 if (present(write_data)) then
   if (.not.write_data) return
 end if

 ! Write data.
 write(unit,"(a)")"# energy           DOS_TOT          IDOS_TOT         DOS[spin=UP]     IDOS[spin=UP] ..."
 do iw=1,Edos%nw
   write(unit,'(es17.8)',advance='no')Edos%mesh(iw)*cfact
   do spin=0,Edos%nsppol
     write(unit,'(2es17.8)',advance='no')Edos%dos(iw,spin)/cfact,Edos%idos(iw,spin)
   end do
   write(unit,*)
 end do

end subroutine edos_print
!!***

END MODULE m_ebands
!!***
