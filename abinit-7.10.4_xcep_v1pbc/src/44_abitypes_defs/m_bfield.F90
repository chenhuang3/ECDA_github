!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_bfield
!! NAME
!!  m_bfield
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle magnetic fields
!!
!! COPYRIGHT
!! Copyright (C) 2011-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_bfield

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_pawcprj, only : pawcprj_type, pawcprj_destroy

 implicit none

 private
!!***


!!****t* defs_datatypes/bfield_type
!! NAME
!! bfield_type
!!
!! FUNCTION
!! First-principles calculations in a finite magnetic field
!!
!! SOURCE

 type, public :: bfield_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer variables
  integer :: orbmag              ! value of orbmag in use
  integer :: fmkmem              ! number of k-points in the FBZ per cpu
  integer :: fmkmem_max          ! max of fmkmem
  integer :: fnkpt               ! number of k-points in the FBZ
  integer :: has_expibi          ! 2 if expibi computed, 1 if only allocated, zero else
  integer :: has_twexpibr         ! 2 if twexpibr computed, 1 if only allocated, zero else
  integer :: has_Lij             ! 2 if paw_Lij computed, 1 if only allocated, zero else
  integer :: has_Lijr3           ! 2 if Lijr3 computed, 1 if only allocated, zero else
  integer :: has_twdij0          ! 2 if twdij0 computed, 1 if only allocated, zero else
  integer :: has_tweijkl         ! 2 if tweijkl computed, 1 if only allocated, zero else
  integer :: has_qijb            ! 2 if paw_qijb computed, 1 if only allocated, zero else
  integer :: lmax
  integer :: lmnmax
  integer :: lmn2max
  integer :: mhcg                ! 2nd dimension of hcg array
  integer :: mkmem_max           ! max of mkmem
  integer :: natom               ! number of atoms in unit cell
  integer :: my_natom            ! number of atoms treated by current proc
  integer :: nband_occ           ! number of occupied bands
                                 ! this number must be the same for every k
  integer :: ndij                ! size of dij due to spinors (nspinor**2)
  integer :: nspinor             ! nspinor input from data set
  integer :: nsym
  integer :: usecprj             ! 1 if bfield%cprj allocated (see below), 0 else
  integer :: usepaw              ! 1 if a PAW calculation, 0 else

! Integer arrays
  integer :: indhk(6,6)          ! index of phase twist terms for <u_k1|H_k2|u_k3> type structures.
                                 ! in this case there are 24 distinct terms. See set_twind.F90

  integer :: twind(6,6)          ! index of phase twist terms. See set_twind.F90

! Real(dp) scalars
  real(dp) :: sdeg               ! spin degeneracy: sdeg = 2 if nsppol = 1
                                 !                         1 if nsppol = 2

! Real(dp) arrays
  real(dp) :: bfield(3)          ! bfield vector in real space, atomic units
  real(dp) :: dkvecs(3,3)        ! dkvec(:,idir) = vector between a k-poinit
                                 ! and its nearest neighbour along idir
  real(dp) :: mag_cart(3)        ! magnetization in cartesian coordinates

! Integer arrays
  integer, allocatable :: atom_indsym(:,:,:) ! atom_indsym(4,natom,nsym)
                                         ! this is data on how the symmetries map the atoms in the cell
                                         ! see symatm.F90 for full description
  integer, allocatable :: cgindex(:,:)    ! cgindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cg array
  integer, allocatable :: cgqindex(:,:,:) ! cgqindex(3,6,nkpt*nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cgq and pwnsfacq
                                      ! arrays
                                      ! (see vtorho.f and initberry.f)
  integer, allocatable :: cprjindex(:,:)  ! cprjindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the cprj in the cprj array (used only
                                      ! for PAW calculations)
  integer, allocatable :: fkgindex(:)     ! same as kgindex, but defined
                                      ! for the FBZ and intended to use
                                      ! with pwindf
  integer, allocatable :: ikpt_dk(:,:,:)  ! ikpt_dk(nkpt,2,3)
                                      ! ikpt_dp(ikpt,ii,idir) = index of the
                                      ! k-point at k+dk (ii=1) and k-dk (ii=2)
  integer, allocatable :: indkk_f2ibz(:,:)   ! indkk_f2ibz(1:dtbfield%fnkpt,1:6)
                                         ! information needed to fold a
                                         ! k-point in the FBZ into the IBZ;
                                         ! the second index (1:6)
                                         ! is as described in listkk
  integer, allocatable :: i2fbz(:)           ! i2fbz(1:nkpt) gives index of IBZ
                                         ! k-points in the FBZ k-point list

  integer, allocatable :: kgindex(:)      ! kgind(nkpt)
                                      ! kgind(ikpt) = ikg

  integer, allocatable :: lmn_size(:)        ! lmn_size(ntypat)
  integer, allocatable :: lmn2_size(:)       ! lmn2_size(ntypat)

  integer, allocatable :: nneigh(:)          ! nneigh(nkpt)
                                         ! for each k-point, nneigh stores
                                         ! the number of its nearest neighbours
                                         ! that are not related by symmetry
  integer, allocatable :: sflag(:,:,:,:)  ! sflag(nband_occ,nkpt*nsppol,2,3)
                                      ! sflag = 0 : compute the whole row of
                                      !             smat
                                      ! sflag = 1 : the row is up to date

! Real(dp) arrays

  real(dp), allocatable :: chern_k(:,:,:)
! chern_k(2,nkpt,dir)
! Chern-form at each k-point. Used for magnetization.

  real(dp),allocatable :: emat(:,:,:)
! emat(2,nband_occ,nkpt*nsppol)
! Diagonal H matrix for every k-point. Used for magnetization.
! stores <u_nk|H_k|u_nk>

  real(dp), allocatable :: expibi(:,:,:)
! expibi(2,my_natom,3)
! used for PAW field calculations (distributed over atomic sites)
! stores the on-site phase factors arising from
! $\langle\phi_{i,k}|\phi_{j,k+\sigma_k k}\rangle$
! where $\sigma = \pm 1$. The on-site phase factor is 
! $\exp[i( - k_k)\cdot I]$ where
! $I$ is the nuclear position. 

  real(dp), allocatable :: twexpibr(:,:,:,:)
! twexpibr(2,my_natom,nfgd,6)  (distributed over atomic sites)
! stores the on-site phase factors arising from
! $\langle\phi_{i,k+\sigma_b k_b}|\phi_{j,k+\sigma_k k_k}\rangle$
! where $\sigma = \pm 1$. The on-site phase factor is 
! $\exp[i( \sigma_b k_b - \sigma_k k_k)\cdot r]$ where
! $r$ is the position on the fine grid. Only 6 values are saved:
! -b1 - (-k2)
! -b1 - (+k2)
! -b2 - (-k3)
! -b2 - (+k3)
! -b3 - (-k1)
! -b3 - (+k1)

  real(dp), allocatable :: fkptns(:,:)       ! fkptns(3,1:dtbfield%fnkpt)
                                         ! k-points in FBZ

  real(dp), allocatable :: hcg(:,:,:)
! hcg(2,mhcg,3)
! Script H wavefunctions at each k point, used only for finite field
! magnetization. Stored as
! hcg(2,npw_k*nband*nkpt,idir) like cg

! pointer to on-site angular momentum
  real(dp),allocatable :: Lij(:,:,:,:) ! Lij(2,lmn2_size_max,ntypat,3)
! gives <(r-R) \wedge p> at each atom type in each of 3 directions
! these are used only in the PAW case with magnetic field

! pointer to on-site angular momentum
  real(dp),allocatable :: Lijr3(:,:,:,:) ! Lijr3(2,lmn2_size_max,ntypat,3)
! gives <(r-R) \wedge p/|r-R|^3> at each atom type in each of 3 directions
! these are used only in the PAW case with magnetic field to get shielding

! pointer to magnetization
! this list gives the magnetization at each k point
! in each direction. 
  real(dp),allocatable :: mag_k(:,:,:) ! mag_k(2,idir,ikpt)

! pointer to onsite local magnetization
! this list gives the on-site <L> part of the magnetization at each k point
! in each direction. Used for berryopt = 5 and -5.
  real(dp),allocatable :: mag_local_k(:,:) ! mag_local_k(idir,ikpt)

  real(dp), allocatable :: qijb_kk(:,:,:,:)
! qijb_kk(2,lmnmax,natom,3)
! on-site part of <u_nk|u_mk+b> matrix elements, 
!
! others are obtained by appropriate complex conjugation

  real(dp), allocatable :: smat(:,:,:,:,:,:)
! smat(2,nband_occ,nband_occ,nkpt*nsppol,2,3)
! Overlap matrix for every k-point. In an electric field calculation,
! smat is updated at every iteration.

  real(dp), allocatable :: twexpibi(:,:,:)
! twexpibi(2,my_natom,6)
! used for PAW field calculations (distributed over atomic sites)
! stores the on-site phase factors arising from
! $\langle\phi_{i,k+\sigma_b k_b}|\phi_{j,k+\sigma_k k_k}\rangle$
! where $\sigma = \pm 1$. The on-site phase factor is 
! $\exp[i( \sigma_b k_b - \sigma_k k_k)\cdot I]$ where
! $I$ is the nuclear position. Only 6 values are saved:
! -b1 - (-k2)
! -b1 - (+k2)
! -b2 - (-k3)
! -b2 - (+k3)
! -b3 - (-k1)
! -b3 - (+k1)
! 

  real(dp), allocatable :: twh(:,:,:,:,:)
! twh(2,nband_occ,nband_occ,nkpt*nsppol,tind)
! Overlap H matrix for every k-point. Used for magnetization.
! stores <u_nk_b|H_k|u_mk_k>
! tind is an index combining the bra (bdir, bfor) and ket
! (kdir, kfor) directions. Recall that bdir = 1,2,3 are the
! three directions in the cell in recip space, and for each one
! we have bfor = 1,2 for forward, backward. Likewise on the
! ket side for kdir and kfor. Then the ket elements are
! indexed as 2*kdir-kfor+1, and the bra elements as
! 2*bdir-bfor+1. These are combined as
! tind = 6*( (2*kdir-kfor+1) - 1) + 2*bdir-bfor+1

  real(dp), allocatable :: twdij0(:,:,:,:,:)
! twdij0(2,lmnmax,lmnmax,my_natom,24)  (distributed over atomic sites)
! for each atom, on-site Dij0 terms including phase shifts
! these include (Torrent CMS 42, 337 (2008) appendix E)
! (1) kinetic <phi_i|exp(I*b_b.r)(-del^2/2)exp(-I*b_k.r)|phi_j>
! (2b) Hartree n_ZC
! (2e) Hartree n_ZC charge compensation
! needed for magnetic fields

  real(dp), allocatable :: twdij(:,:,:,:,:,:)
! twdij(2,lmnmax,lmnmax,natom,24,ndij)
! for each atom, on-site Dij terms including phase shifts
! needed for magnetic fields

  real(dp), allocatable :: tweijkl(:,:,:,:,:)
! tweijkl(2,lmn2max,lmn2max,my_natom,6)  (distributed over atomic sites)
! for each atom, on-site e_ijkl terms including phase shifts
! these include (Torrent CMS 42, 337 (2008) appendix E)
! (2a*) on-site Hartree
! (2c*) \hat{n} Hartree
! (2d*) \tilde{n}^1 Hartree and charge compensation
! (2f*) \hat{n} Hartree and charge compensation
! needed for magnetic fields
! only six values are saved as in dtbfield%twexpibi

  real(dp), allocatable :: twqijb_kk(:,:,:,:)
! qijb_kk(2,lmnmax,natom,6)
! on-site part of <u_nk+sig_ib_i|u_mk+sig_jbj> matrix elements, 
! only 6 values are saved:
! -b1 - (-k2)
! -b1 - (+k2)
! -b2 - (-k3)
! -b2 - (+k3)
! -b3 - (-k1)
! -b3 - (+k1)
! 

  real(dp), allocatable :: zarot(:,:,:,:)
   !  zarot(l_size_max,l_size_max,l_max,nsym)
   !  Coeffs of the transformation of real spherical
   !  harmonics under the symmetry operations. These are needed when the
   ! cprj's need to be computed in the full BZ, that is,
   ! in the PAW case with kptopt /= 3.

! pointer to cprj
   type(pawcprj_type),allocatable :: cprj(:,:)
! used with finite bfield and PAW

 end type bfield_type

 ! Bound methods:
 public :: destroy_bfield
!!***

contains

!!****f* m_bfield/destroy_bfield
!! NAME
!!
!! FUNCTION
!!   deallocate fields in bfield structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_bfield(dtbfield)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_bfield'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(bfield_type),intent(out) :: dtbfield

! ************************************************************************

! Integer pointers
  if(allocated(dtbfield%atom_indsym))  then
    ABI_DEALLOCATE(dtbfield%atom_indsym)
  end if
  if(allocated(dtbfield%cgindex))  then
    ABI_DEALLOCATE(dtbfield%cgindex)
  end if
  if(allocated(dtbfield%cgqindex))  then
    ABI_DEALLOCATE(dtbfield%cgqindex)
  end if
  if(allocated(dtbfield%cprjindex))  then
    ABI_DEALLOCATE(dtbfield%cprjindex)
  end if
  if(allocated(dtbfield%fkgindex))  then
    ABI_DEALLOCATE(dtbfield%fkgindex)
  end if
  if(allocated(dtbfield%ikpt_dk))  then
    ABI_DEALLOCATE(dtbfield%ikpt_dk)
  end if
  if(allocated(dtbfield%indkk_f2ibz))  then
    ABI_DEALLOCATE(dtbfield%indkk_f2ibz)
  end if
  if(allocated(dtbfield%i2fbz))  then
    ABI_DEALLOCATE(dtbfield%i2fbz)
  end if
  if(allocated(dtbfield%kgindex))  then
    ABI_DEALLOCATE(dtbfield%kgindex)
  end if
  if(allocated(dtbfield%lmn_size))  then
    ABI_DEALLOCATE(dtbfield%lmn_size)
  end if
  if(allocated(dtbfield%lmn2_size))  then
    ABI_DEALLOCATE(dtbfield%lmn2_size)
  end if
  if(allocated(dtbfield%nneigh))  then
    ABI_DEALLOCATE(dtbfield%nneigh)
  end if
  if(allocated(dtbfield%sflag))  then
    ABI_DEALLOCATE(dtbfield%sflag)
  end if

! Real(dp) pointers

  if(allocated(dtbfield%chern_k))  then
    ABI_DEALLOCATE(dtbfield%chern_k)
  end if
  if(allocated(dtbfield%emat))  then
    ABI_DEALLOCATE(dtbfield%emat)
  end if
  if(allocated(dtbfield%expibi))  then
    ABI_DEALLOCATE(dtbfield%expibi)
  end if
  if(allocated(dtbfield%twexpibr))  then
    ABI_DEALLOCATE(dtbfield%twexpibr)
  end if
  if(allocated(dtbfield%fkptns))  then
    ABI_DEALLOCATE(dtbfield%fkptns)
  end if
  if(allocated(dtbfield%hcg))  then
    ABI_DEALLOCATE(dtbfield%hcg)
  end if
  if(allocated(dtbfield%Lij))  then
    ABI_DEALLOCATE(dtbfield%Lij)
  end if
  if(allocated(dtbfield%Lijr3))  then
    ABI_DEALLOCATE(dtbfield%Lijr3)
  end if
  if(allocated(dtbfield%mag_k))  then
    ABI_DEALLOCATE(dtbfield%mag_k)
  end if
  if(allocated(dtbfield%mag_local_k))  then
    ABI_DEALLOCATE(dtbfield%mag_local_k)
  end if
  if(allocated(dtbfield%qijb_kk))  then
    ABI_DEALLOCATE(dtbfield%qijb_kk)
  end if
  if(allocated(dtbfield%smat))  then
    ABI_DEALLOCATE(dtbfield%smat)
  end if
  if(allocated(dtbfield%twexpibi))  then
    ABI_DEALLOCATE(dtbfield%twexpibi)
  end if
  if(allocated(dtbfield%twh))  then
    ABI_DEALLOCATE(dtbfield%twh)
  end if
  if(allocated(dtbfield%twdij0))  then
    ABI_DEALLOCATE(dtbfield%twdij0)
  end if
  if(allocated(dtbfield%twdij))  then
    ABI_DEALLOCATE(dtbfield%twdij)
  end if
  if(allocated(dtbfield%tweijkl))  then
    ABI_DEALLOCATE(dtbfield%tweijkl)
  end if
  if(allocated(dtbfield%twqijb_kk))  then
    ABI_DEALLOCATE(dtbfield%twqijb_kk)
  end if
  if(allocated(dtbfield%zarot))  then
    ABI_DEALLOCATE(dtbfield%zarot)
  end if

! pointer to cprj
  if(allocated(dtbfield%cprj)) then
    call pawcprj_destroy(dtbfield%cprj)
    ABI_DATATYPE_DEALLOCATE(dtbfield%cprj)
  end if

end subroutine destroy_bfield
!!***

end module m_bfield
!!***
