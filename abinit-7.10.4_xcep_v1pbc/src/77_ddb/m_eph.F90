!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_eph
!! NAME
!!  m_eph
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods 
!!  used to perform electron-phonon calculations. 
!! ***** Work in progress, do not rely on this implementation. type declarations might be changed *****
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (MG,MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! Most important requirements of the new e-ph implementation:
!!
!! 1) Reduce the memory requirements:
!!    * Memory distribution via MPI.
!!    * Only states in an energy window close to E_F should be stored
!!    * Only the irreducible k-points for a given q and (ipert,idir) should be stored.
!!      Specialized routines will be used to retrieve the symmetrized matrix elements in the full BZ.
!!    * Only the irreducible pertubations at given (q, idir, ipert) should be stored and 
!!      explicitly calculated at the DFPT level.
!!
!! 2) Methods hiding the ugly details of the implementation.
!!
!! 3) out-of-core and in-core solutions (ETSF-IO support should facilitate the implementation but
!!    final specifications are still missing)
!!
!! 4) The objects should be designed taking int account a possible use of the Wannier interpolation.
!!    Since non-homogeneous k-meshes lead to a dramatic increase in the number of k-points, all the
!!    arrays dimensioned with nkibz and nkbz should be carefully designed to reduce memory.
!!
!! Issues:
!!
!! 1) Several quantities depend on the spin : Fermi surface, E_f, eigen, occ, gkk, set of bands around E_f
!!
!! 2) Spinorial case: each gkk is a two-by-two matrix.
!!
!! 3) Only those symmetries which preserve q, and (iatom, idir) can be used to reconstruct the gkk's in the full BZ.
!!    ==> for given q and (idir,ipert) we have IBZ_q^{iatom,idir}, 
!!
!! 4) To facilitate the implementation, the new E-PH part will make use of several data types already introduced in 
!!    the GW part:
!!     crystal_t 
!!     band_structure
!!     kmesh_t
!!
!!    TODO: Better integration between GW data types and the main abinit code is needed
!!    e.g treatment of k-points, time-reversal, Fermi level for semiconductors...
!!
!! 5) To speed up the k-point search either we use the previous approach (rank and invrank)
!!    or we order the points according to their length. The later approach has the advantage
!!    of being applicable also to non-homogeneous k-meshes. Moreover it is less memory demanding 
!!    as the storage of invrank is not needed anymre. On the other hand, looping over shells  
!!    is expected to be a bit slower than the algorithm based on the k-point rank.
!!
!! 6) Which representation for the gkk? The one presently used i.e. (q,idir,ipert) 
!!    or the phonon representation (q,nu)?
!!    One should check whether (q,nu) leads to some simplification when symmetries are used 
!!    to complete the matrix elements.
!!
!! Questions:
!! 1) At Gamma and border zone, \Delta V_{SCF} is Hermitian thus we can store the upper triangle 
!!    of the (b1,b2) matrix. Is it worth taking into account this possibility?
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_eph

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

 private 
!!***

!!****t* m_eph/fsurf_t
!! NAME
!! fsurf_t
!!
!! FUNCTION
!!  Stores data and table related to the Fermi surface for a single spin polarization.
!!  Only the set of k-points where at least one band falls within an energy 
!!  window centered at Ef are stored in memory.
!!
!! SOURCE

 type,public :: fsurf_t

  integer :: nkbz
  ! Number of k-points on the Fermi surface in the full Brillouin zone (FS BZ).

  integer :: nkibz
  ! Number of k-points on the Fermi surface in the irreducible Brillouin zone (FS IBZ).

  integer :: nband
  ! Number of bands around the Fermi level falling within the energy window [ef-ewidth, ef+ewidth].

  integer :: bstart
  ! Initial band index.

  integer :: nsym
  ! Number of symmetries in the full point group.

  integer :: timrev
  ! 1 or 2 depending whether time-reversal cannot or can be used (respectively). 
  ! TODO use abinit conventions.

  integer :: occopt
  ! Option for metallic occupation.

  real(dp) :: efermi
  ! Fermi level.

  real(dp) :: ewidth
  ! Width of the energy window around Ef

  real(dp) :: tsmear
  ! Value of the broadening for occupation factors.

  real(dp) :: dosef 
  ! Dos at the Fermi level for this spin channel.

  real(dp) :: tol_weight
  ! Tolerance on the weight for FS integration. 
  ! A (k,b) state whose weight is smaller than tolweight won"t be stored.

  integer :: mesh_type
  ! 1 for homogeneous mesh
  ! 2 for random mesh (e.g. for Wannier interpolation)

  integer :: int_opt
  ! Flag defining the technique used to perform the integration of the FS.
  !  1 for gaussian.
  !  2 for tetrahedrons.

  real(dp) :: elphsmear
  ! Standard deviation of the gaussian function used to integrate over the FS. Only used if int_opt=1

  integer :: kptrlatt(3,3)=RESHAPE((/0,0,0,0,0,0,0,0,0/),(/3,3/))
  ! Coordinates of three vectors in real space, expressed in reduced coordinates.
  ! They defines a super-lattice in real space. The k point lattice is the reciprocal of
  ! this super-lattice, eventually shifted by shift. ONLY used if mesh_type == 1

  ! real(dp),pointer :: shift(3) !shift(:,:)
   !  shift(3,nshift)
   !  shift for k-points, usually nshift=1

  integer,allocatable :: bz2ibz(:)  
  ! bz2ibz(nkpt) 
  ! For each point in the FS BZ, it gives the index of the symmetrical image in the IBZ.

  !integer,allocatable :: symopbz2ibz(:)
  ! symopbz2ibz(nkbz)
  ! Index of the symmetry operation such that IS k_ibz = k_bz
  ! where S is one of the symrec ops and I is the identity or the inversion, depending on symtrbz2ibz.

  !integer,allocatable :: symtrbz2ibz(:)
  ! symtrbz2ibz(nkbz)
  ! Index of the symmetry operation such that IS k_ibz = k_bz

  !integer, allocatable :: kbz2sh(:)
  ! kbz2sh(nkbz)
  ! For each k-point in the FS BZ, it gives the index of the shell to which it belongs.
  
  !integer,allocatable :: shlim(:)
  ! shlim(nksh+1)
  ! Index of the first k-point in each shell, =nkbz+1 for nksh+1

  !integer :: nksh                             ! Number of shells of k-points.

  !real(dp),allocatable :: shlen(:)
  ! shlen(nksh)
  ! Radius of each shell.

  real(dp) :: gprimd(3,3)
   ! Dimensional reciprocal space primitive translations (Bohr^-1)

  real(dp),allocatable :: kibz(:,:) 
  ! kibz(3,nkibz)
  ! Reduced coordinates of the points in the FS IBZ. Ordered by increasing module.

  real(dp),allocatable :: kbz(:,:) 
  ! kbz(3,nkbz)
  ! Reduced coordinates of the points in the FS. Ordered by increasing module.
  ! TODO is it really needed? It might be large. Might be recalculated on-the-fly
  
  real(dp),allocatable :: eigen(:,:)
  ! eigen(nband,nkibz)
  ! Energies in the FS IBZ.

  real(dp),allocatable :: occ(:,:)
  ! occ(nband,nkibz)
  ! Occupation number in the FS IBZ.

  real(dp),allocatable :: fs_weight(:,:)
  ! fs_weight(nband,nkibz)
  ! Weights due to the delta function centered at the Fermi level.
  ! Two methods are available, standard gaussian and tetrahedron methods, depending on int_opt.

  real(dp),allocatable :: wtk(:)
  ! wtk(nkibz)
  ! Weights for each point on the FS IBZ. Normalized to one.

 end type fsurf_t
!!***
 
 ! Bound methods:
 !public :: init_fermi_surface
 public :: fsurf_free
 !public :: get_fs_ibz

 ! example:
 !type(fsurf_t),allocatable :: Fsurf(:)
 !
 !allocate(Fsurf(nsppol))
 !do spin=1,nsppol
 ! call init_Fermi_surface(Fsurf(spin),Cryst,Kmesh,Bands,tolweight)
 !end do

 interface fsurf_free
   module procedure fsurf_free_0D
   module procedure fsurf_free_1D
 end interface fsurf_free

!----------------------------------------------------------------------

!!****t* m_eph/gkk_t
!! NAME
!! gkk_t
!!
!! FUNCTION
!! Structure used to store the matrix elements at fixed (q, idir, ipert)
!!
!! SOURCE

 type,public :: gkk_t

  integer :: nkibz_gkk
  ! Number of points in the IBZ_q^{idir,ipert}

  integer :: fs_nkbz
  ! Numbe of points on the FS BZ.

  integer :: fs_nband
  ! Number of states around E_f. Equivalent to the value stored in the fsurf_t.

  integer :: bstart
  ! First band index. Equivalent to the value stored in the fsurf_t.

  integer :: nsym_gkk
  ! Number of symmetries preserving (q,idir,ipert)

  integer :: timrev
  ! 1 or 2 depending whether time-reversal cannot or can be used (respectively). 
  ! TODO use abinit conventions.

  integer,allocatable :: symgkk2symrec(:)
  ! symgkk2symrec(nsym_gkk)
  ! Index in the full array symrec of the symmetries preserving (qpt,idir,iper)

  integer,allocatable :: ibzgkk2fsbz(:)
  ! ibzgkk2fsbz(nkibz_gkk)
  ! Index in the full FS BZ of each point in the IBZ_q^{idir,ipert}
  
  integer,allocatable :: fsbz2ibzgkk(:) 
  ! fsbz2ibzgkk(fs_nkbz)
  ! For each point in the full BZ gives the index in my IBZ of the symmetrical point. 

  integer,allocatable :: symtrbz2ibz(:)
  ! symtrbz2ibz(fs_nkbz)
  ! 1 or 2, depending wheter time-reversal has to be used to obtain the point.

  integer,allocatable :: symopbz2ibz(:) 
  ! ksymopbz2ibz(fs_nkbz)
  ! Index of the symmetry operation in the set of my symmetries such that I S kibz = kbz
  ! TODO Do we need unklapp vectors

  real(dp),allocatable :: gkk(:,:,:,:)  
  ! gkk(2,fs_nband,fs_nband,nkibz_gkk)

  !integer,allocatable :: symrec_gkk(:,:,:)
  ! symrec_gkk(3,3,nsym_gkk)
  ! symmetry operations in reciprocal space preserving the external q and the perturbation (idir,ipert)
                                                                                                  
  !integer,allocatable :: symafm_gkk(:)
  ! symafm_gkk(nsym_gkk)

 end type gkk_t

 public :: gkk_free
 !gkk_init
 !gkk_read
 !gkk_full_fsbz        ! complete gkk on the full FS BZ.
!!***

 interface gkk_free
   module procedure gkk_free_0D
   module procedure gkk_free_1D
 end interface gkk_free

!----------------------------------------------------------------------

!!****t* m_eph/gkk_handler_type
!! NAME
!! gkk_handler_type
!!
!! FUNCTION
!! Structure used to store the matrix elements for a given q.
!!
!! SOURCE

 type,public :: gkk_handler_type

  integer :: natom
  ! Number of atoms in the unit cell.

  integer :: nirr_perts
  ! Number of irreducible perturbations (idir,ipert) <= 3*natom.

  character(len=fnlen) :: fname
  ! Name of the file storing the gkk matrix elements (either complete database of single q-point and spin index).
  ! Used in the case of out-of-core solution.

  integer,allocatable :: pert_list(:,:)
  ! pert_list(2,nirr_perts)
  ! gives (idir,ipert) for each irreducible perturbation.

  real(dp) :: qpt(3)
  ! The q-point in reduced coordinates.

  type(gkk_t),allocatable :: Gkk(:)   
  ! gkk(nirr_perts)
  ! gkk matrix elements for this q-point 

 end type gkk_handler_type 
!!***

! Bound Methods: 
 !public :: gkk_free_handler
 !% init_gkk_handler(Gkk,FSurf,Cryst,Cryst,qpt,fname)
 !% get_gammaq
 !% symmetrize_gkk_over_perts

 !interface gkk_free_handler
 !  module procedure gkk_free_handler_0D
 !  module procedure gkk_free_handler_1D
 !end interface gkk_free_handler

 !type(gkk_handler_type),allocatable :: Gkk(:,:)
 !allocate(Gkk(nqibz,nsppol))

 !do spin=1,nsppol
 ! do iqibz=1,nqibz
 !  if (I_treat(iqibz,spin)) then
 !   call init_gkk_handler(Gkk(iqibz,spin),Fsurf(spin),Cryst,qpt,fname)
 !  end if
 ! end do
 !end do

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_eph/fsurf_free_0D
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the fsurf_t structure
!!  0 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine fsurf_free_0D(FSurf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fsurf_free_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(fsurf_t),intent(inout) :: FSurf

! ************************************************************************

 !@fsurf_t

 ! integer 
 if (allocated(Fsurf%bz2ibz)) then
   ABI_FREE(Fsurf%bz2ibz)
 end if
 !if (allocated(Fsurf%symopbz2ibz)) !deallocate(Fsurf%symopbz2ibz)
 !if (allocated(Fsurf%symtrbz2ibz)) !deallocate(Fsurf%symtrbz2ibz)
 !if (allocated(Fsurf%kbz2sh       )) !deallocate(Fsurf%kbz2sh       )
 !if (allocated(Fsurf%shlim      )) !deallocate(Fsurf%shlim      )
 !if (allocated(Fsurf%shlen      )) !deallocate(Fsurf%shlen      )

 ! real
 if (allocated(Fsurf%kibz)) then
   ABI_FREE(Fsurf%kibz)
 end if
 if (allocated(Fsurf%kbz)) then
   ABI_FREE(Fsurf%kbz)
 end if
 if (allocated(Fsurf%eigen)) then
   ABI_FREE(Fsurf%eigen)
 end if
 if (allocated(Fsurf%occ)) then
   ABI_FREE(Fsurf%occ)
 end if
 if (allocated(Fsurf%fs_weight)) then
   ABI_FREE(Fsurf%fs_weight)
 end if
 if (allocated(Fsurf%wtk)) then
   ABI_FREE(Fsurf%wtk)
 end if

end subroutine fsurf_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/fsurf_free_1D 
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the fsurf_t structure
!!  1 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine fsurf_free_1D(FSurf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fsurf_free_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(fsurf_t),intent(inout) :: FSurf(:)

!Local variables-------------------------------
!scalars
 integer :: spin

! ************************************************************************

 do spin=1,SIZE(FSurf)
   call fsurf_free_0D(FSurf(spin))
 end do

end subroutine fsurf_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/gkk_free_0D
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the gkk_t structure
!!  0 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gkk_free_0D(Gkk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkk_free_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_t),intent(inout) :: Gkk

! ************************************************************************

 !@gkk_t
                            
 ! integer 
 if (allocated(Gkk%symgkk2symrec)) then
   ABI_FREE(Gkk%symgkk2symrec)
 end if
 if (allocated(Gkk%ibzgkk2fsbz)) then
   ABI_FREE(Gkk%ibzgkk2fsbz)
 end if
 if (allocated(Gkk%fsbz2ibzgkk)) then
   ABI_FREE(Gkk%fsbz2ibzgkk)
 end if
 if (allocated(Gkk%symtrbz2ibz)) then
   ABI_FREE(Gkk%symtrbz2ibz)
 end if
 if (allocated(Gkk%symopbz2ibz)) then
   ABI_FREE(Gkk%symopbz2ibz)
 end if
                            
 ! real
 if (allocated(Gkk%gkk)) then
   ABI_FREE(Gkk%gkk)
 end if
 !if (allocated(Gkk%symrec_gkk)  deallocate(Gkk%symrec_gkk)
 !if (allocated(Gkk%symafm_gkk)  deallocate(Gkk%symafm_gkk)

end subroutine gkk_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/gkk_free_1D 
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the gkk_t structure
!!  1 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gkk_free_1D(Gkk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkk_free_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_t),intent(inout) :: Gkk(:)

!Local variables-------------------------------
!scalars
 integer :: ii

! ************************************************************************

 do ii=1,SIZE(Gkk)
   call gkk_free_0D(Gkk(ii))
 end do

end subroutine gkk_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/gkk_free_handler_0D
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the gkk_handler_type structure
!!  0 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_gkk_handler_0D(Gkk_hdl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_gkk_handler_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_handler_type),intent(inout) :: Gkk_hdl

! ************************************************************************

 !@gkk_handler_type
                            
 ! integer 
 if (allocated(Gkk_hdl%pert_list)) then
   ABI_FREE(Gkk_hdl%pert_list)
 end if

 ! substructures
 call destroy_gkk_1D(Gkk_hdl%Gkk)

end subroutine destroy_gkk_handler_0D
!!***

!----------------------------------------------------------------------

!!****f* m_eph/destroy_gkk_handler_1D
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the gkk_handler_type structure
!!  1 dimensional Fermi surface case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_gkk_handler_1D(Gkk_hdl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_gkk_handler_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkk_handler_type),intent(inout) :: Gkk_hdl(:)

!Local variables-------------------------------
!scalars
 integer :: ii

! ************************************************************************

 do ii=1,SIZE(Gkk_hdl)
  call destroy_gkk_handler_0D(Gkk_hdl(ii))
 end do

end subroutine destroy_gkk_handler_1D
!!***

!----------------------------------------------------------------------

END MODULE m_eph
!!***
