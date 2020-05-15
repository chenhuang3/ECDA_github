!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ifc
!! NAME
!!  m_ifc
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods 
!!  used to handle interatomic force constant sets 
!!
!! COPYRIGHT
!! Copyright (C) 2011-2014 ABINIT group (XG,MJV,EB,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_ifc

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_io_tools,    only : get_unit
 use m_copy,        only : alloc_copy
 use m_anaddb_dataset, only : anaddb_dataset_type
 use m_ewald,       only : ewald9
 use m_crystal,     only : crystal_t
 use m_dynmat,      only : canct9, dist9 , ifclo9, axial9, q0dy3_apply, q0dy3_calc, asrif9, &
&                          make_bigbox, canat9, chkrp9, ftifc_q2r, wght9, symdm9, nanal9, gtdyn9, dymfz9
 use m_ddb,         only : ddb_type

 implicit none

 private 
!!***

!!****t* m_ifc/ifc_type
!! NAME
!! ifc_type
!!
!! FUNCTION
!!  contains the necessary data to interpolate the
!!  phonon bandstructure and eigenvectors in reciprocal space (ie.
!!  interatomic force constants and corresponding real space grid info).
!!
!! SOURCE

 type,public :: ifc_type

   integer :: natom
     ! Number of atoms in the unit cell.

   integer :: mpert 
     ! Maximum number of ipert.

   !integer :: ifcflag
   ! 0 if Fourier interpolation should not be used, 1 otherwise

   integer :: asr
     ! Option for the treatment of the Acoustic Sum Rule.

   integer :: brav
     ! Option for the sampling of the BZ (anaddb input variable)

   integer :: dipdip
     ! dipole dipole interaction flag.
                                                                                              
   integer :: symdynmat
     ! If equal to 1, the dynamical matrix is symmetrized in phfrq3 before the diagonalization.

   integer :: nqshft
     ! Number of shifts in the q-mesh (usually 1 since the mesh is gamma-centered!)

   integer :: nqbz
     ! Number of points in the full BZ

   integer :: nrpt
     ! Number of real space points used to integrate IFC (for interpolation of dynamical matrices)

   integer :: ngqpt(3)
    ! Number of division in the Q mesh.

   real(dp) :: rprim(3,3),gprim(3,3),acell(3)
     ! These values are used to call anaddb routines the don't use rprimd, gprimd

   ! These values will be stored in ddb but then we have to pass ddb to ifc_fourq
    real(dp) :: dielt(3,3)
     ! Dielectric tensor

   real(dp),allocatable :: amu(:)
     ! amu(ntypat)
     ! mass of the atoms (atomic mass unit)

   real(dp),allocatable :: atmfrc(:,:,:,:,:,:)
     ! atmfrc(2,3,natom,3,natom,nrpt)
     ! Inter atomic forces in real space

   real(dp),allocatable :: qshft(:,:)
    ! qshft(3,nqshft)
    ! The shifts of the q-mesh

   real(dp), allocatable :: rpt(:,:)
     ! rpt(3,nrpt)
     ! Real space points in canonical type coordinates.

   real(dp),allocatable :: wghatm(:,:,:)
     ! wghatm(natom,natom,nrpt)
     ! Weights for each point and atom in the Wigner Seitz supercell in real space.

   real(dp),allocatable :: rcan(:,:)
     ! rcan(3,natom) 
     ! Atomic position in canonical coordinates.

   real(dp),allocatable :: trans(:,:)
     ! trans(3,natom)
     ! Atomic translations: xred = rcan + trans

   real(dp),allocatable :: dyewq0(:,:,:)
     ! dyewq0(3,3,natom)
     ! Atomic self-interaction correction to the dynamical matrix (only when dipdip=1).

   real(dp),allocatable :: zeff(:,:,:)
     ! zeff(3,3,natom)
     ! Effective charge on each atom, versus electric field and atomic displacement.

   real(dp),allocatable :: qbz(:,:)
     ! qbz(3,nqpt))
     ! List of q-points in the full BZ 

   real(dp),allocatable :: dynmat(:,:,:,:,:,:)
     ! dynmat(2,3,natom,3,natom,nqpt))
     ! dynamical matrices relative to the q points of the B.Z. sampling
     ! Note that the long-range dip-dip part has been removed if dipdip=1
     ! Moreover the array is multiplied by a phase shift in mkifc9.

   !real(dp),allocatable :: dynmat_lr(:,:,:,:,:,:) 
    ! dynmat_lr(2,3,natom,3,natom,nqpt))
    ! Long-range part of dynmat in q-space
 end type ifc_type
 
 public :: ifc_init     ! Constructor
 public :: ifc_free     ! Release memory
 public :: ifc_fourq    ! Use Fourier interpolation to compute interpolated frequencies w(q) and eigenvectors e(q)
 !public :: ifc_diagoq  ! Compute phonon frequencies via direct diagonalization (mainly for debugging purposes)
!!***

!----------------------------------------------------------------------

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_free
!! NAME
!!
!! FUNCTION
!!  Deallocate memory for the ifc_type structure
!!
!! SOURCE

subroutine ifc_free(ifc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(ifc_type),intent(inout) :: ifc

! ************************************************************************
 ! real
 if (allocated(ifc%amu)) then
   ABI_FREE(ifc%amu)
 end if

 if (allocated(ifc%atmfrc)) then
   ABI_FREE(ifc%atmfrc)
 end if

 if (allocated(ifc%qshft)) then
   ABI_FREE(ifc%qshft)
 end if

 if (allocated(ifc%rpt)) then
   ABI_FREE(ifc%rpt)
 end if

 if (allocated(ifc%wghatm)) then
   ABI_FREE(ifc%wghatm)
 end if

 if (allocated(ifc%rcan)) then
   ABI_FREE(ifc%rcan)
 end if

 if (allocated(ifc%trans)) then
   ABI_FREE(ifc%trans)
 end if

 if (allocated(ifc%dyewq0)) then
   ABI_FREE(ifc%dyewq0)
 end if

 if (allocated(ifc%qbz)) then
   ABI_FREE(ifc%qbz)
 end if

 if (allocated(ifc%zeff))then
   ABI_FREE(ifc%zeff)
 end if

 if (allocated(ifc%dynmat)) then
   ABI_FREE(ifc%dynmat)
 end if

 !if (allocated(ifc%dynmat_lr)) then
 !  ABI_FREE(ifc%dynmat_lr)
 !end if

end subroutine ifc_free
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_init
!! NAME
!!  ifc_init
!!
!! FUNCTION
!!  Initialize the dynamical matrix as well as the IFCs.
!!  taking into account the dipole-dipole interaction.
!!
!! INPUTS
!! Crystal<type(crystal_t)> = Information on the crystalline structure.
!! ddb<type(ddb_type)> = Database with derivatives.
!! anaddb_dtset<type(anaddb_dataset_type)> = TODO: TO BE REMOVED
!! dielt(3,3)=dielectric tensor TODO: Should be computed from DDB
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!! [Ifc_coarse]=Optional. 

!! anaddb_dtset= (derived datatype) contains all the input variables
!! ddb = storage object for ddb information read in from DDB file
!! dielt(3,3)=dielectric tensor
!! iout=unit number for output of formatted data
!! mpert =maximum number of ipert
!! natom=number of atoms in cell
!! ngqpt_in = input values of ngqpt
!! ntypat=number of atom types
!! tcpui,twalli=initial values of cpu and wall clocktime
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!!
!! OUTPUT
!! Ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!!   dyewq0(3,3,natom)=atomic self-interaction correction to the dynamical matrix. (only when dipdip=1)
!!   rcan(3,natom) = Atomic position in canonical coordinates
!!   trans(3,natom) = Atomic translations : xred = rcan + trans
!!   Ifc%nrpt = number of vectors of the lattice in real space
!!   Ifc%atmfrc(2,3,natom,3,natom,nrpt)= Interatomic Forces in real space
!!   Ifc%rpt(3,nprt) = Canonical coordinates of the R points in the unit cell
!!   Ifc%wghatm(natom,natom,nrpt)= Weight associated to the couple of atoms and the R vector
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine ifc_init(Ifc,Crystal,ddb,anaddb_dtset,ngqpt_in,dielt,zeff,Ifc_coarse)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_init'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_56_recipspace
 use interfaces_72_response
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ifc_type),intent(inout) :: Ifc
 type(crystal_t),intent(in) :: Crystal
 type(ddb_type),intent(in) :: ddb
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
 type(ifc_type),optional,intent(in) :: Ifc_coarse
!arrays
 integer,intent(in) :: ngqpt_in(3)
 real(dp),intent(in) :: dielt(3,3),zeff(3,3,Crystal%natom)

!Local variables -------------------------
!scalars
 integer :: mpert,iout,iqpt,mqpt,nsym,ntypat,iq_bz,ii,natom
 integer :: nqbz,nqshft,option,plus,sumg0,irpt,irpt_new,rfmeth
 real(dp),parameter :: qphnrm=one
 real(dp) :: tcpui,twalli
 character(len=500) :: message
 type(ifc_type) :: ifc_tmp
!arrays
 integer :: ngqpt(9),qptrlatt(3,3)
 integer,allocatable :: atifc(:),qmissing(:)
 real(dp) :: gprim(3,3),rprim(3,3),q1shft(3,4),qpt(3),rprimd(3,3)
 real(dp):: rcan(3,Crystal%natom),trans(3,Crystal%natom),dyewq0(3,3,Crystal%natom) 
 real(dp) :: displ(2*3*Crystal%natom*3*Crystal%natom)
 real(dp) :: phfrq(3*Crystal%natom) !eigval(3,Crystal%natom),
 real(dp) :: eigvec(2,3,Crystal%natom,3,Crystal%natom)
 real(dp),allocatable :: dyew(:,:,:,:,:),spqpt(:,:)
 real(dp),allocatable :: dynmatfull(:,:,:,:,:,:),dynmat_sr(:,:,:,:,:,:),dynmat_lr(:,:,:,:,:,:) ! for OmegaSRLR
 real(dp),allocatable :: out_d2cart(:,:,:,:,:)

!******************************************************************

 ! This dimension should be encapsulated somewhere. We don't want to 
 ! change the entire code if someone adds a new kind of perturbation.
 mpert = Crystal%natom + 6; iout = ab_out

 rprim = ddb%rprim; gprim = ddb%gprim

 nsym = Crystal%nsym
 natom = Crystal%natom
 ntypat = Crystal%ntypat
 rprimd = Crystal%rprimd

 ngqpt=0; ngqpt(1:3)=ngqpt_in(1:3)

 nqshft = anaddb_dtset%nqshft
 rfmeth = anaddb_dtset%rfmeth
 q1shft = anaddb_dtset%q1shft

 ABI_MALLOC(atifc,(natom))
 atifc(:)=anaddb_dtset%atifc(:)

! Copy important parameters in Ifc
 Ifc%natom = natom
 Ifc%mpert = mpert
 Ifc%asr = anaddb_dtset%asr
 Ifc%brav = anaddb_dtset%brav
 Ifc%dipdip = anaddb_dtset%dipdip
 Ifc%symdynmat = anaddb_dtset%symdynmat
 Ifc%nqshft = anaddb_dtset%nqshft
 call alloc_copy(anaddb_dtset%q1shft(:,1:Ifc%nqshft),Ifc%qshft)
 Ifc%ngqpt = ngqpt_in(1:3)
 Ifc%rprim = ddb%rprim
 Ifc%gprim = ddb%gprim
 Ifc%acell = ddb%acell

 ! Check if the rprim are coherent with the choice used in the interatomic forces generation
 call chkrp9(Ifc%brav,rprim)

 dyewq0 = zero
 if (Ifc%dipdip==1 .and. (Ifc%asr==1.or.Ifc%asr==2)) then
   ! Calculation of the non-analytical part for q=0
   sumg0=0
   qpt(:)=zero
   ABI_MALLOC(dyew,(2,3,natom,3,natom))
   call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,Crystal%xred,zeff)
   option=Ifc%asr
   call q0dy3_calc(natom,dyewq0,dyew,option)
   ABI_FREE(dyew)
 end if

 ! Sample the Brillouin zone
 option=1
 qptrlatt(:,:)=0
 qptrlatt(1,1)=ngqpt(1)
 qptrlatt(2,2)=ngqpt(2)
 qptrlatt(3,3)=ngqpt(3)
 mqpt=ngqpt(1)*ngqpt(2)*ngqpt(3)*nqshft
 if (Ifc%brav==2) mqpt=mqpt/2
 if (Ifc%brav==3) mqpt=mqpt/4

 ABI_MALLOC(spqpt,(3,mqpt))
 call smpbz(Ifc%brav,ab_out,qptrlatt,mqpt,nqbz,nqshft,option,q1shft,spqpt)

 ABI_MALLOC(Ifc%dynmat,(2,3,natom,3,natom,nqbz))

! Find symmetrical dynamical matrices
 if (.not.present(Ifc_coarse)) then

   ! Each q-point in the BZ mush be the symmetrical of one of the qpts in the ddb file.
   call symdm9(ddb%flg,ddb%nrm,ddb%qpt,ddb%typ,ddb%val,&
&    Ifc%dynmat,gprim,Crystal%indsym,mpert,natom,ddb%nblok,nqbz,nsym,rfmeth,rprim,spqpt,&
&    Crystal%symrec,Crystal%symrel)

 else
   ! Symmetrize the qpts in the BZ using the q-points in the ddb.
   ! Then use Ifc_coarse to fill the missing entries with Fourier interpolated matrices.
   !
   ! TODO: The previous version of refineblk was hacking the DDB database to add the q-points in the **IBZ**
   ! Then D(q) was symmetrized in symdm9. This version avoids the symmetrization: the q-points 
   ! in the BZ that are not in the coarse q-mesh are obtained by an explicit FT.
   ! This means that the final D(q) may break some symmetry in q-space if the FT does not preserve it.
   ! The most elegant approach would be to get D(q_ibz) via FT if q_ibz is not in the coarse mesh and then
   ! call symdm9 to get D(q) for each q point in the star of q_ibz.
   call wrtout(std_out,"Will fill missing qpoints in the full BZ using the coarse q-mesh","COLL")

   call symdm9(ddb%flg,ddb%nrm,ddb%qpt,ddb%typ,ddb%val,&
&    Ifc%dynmat,gprim,Crystal%indsym,mpert,natom,ddb%nblok,nqbz,nsym,rfmeth,rprim,spqpt,&
&    Crystal%symrec,Crystal%symrel,qmissing=qmissing)

   ! Compute dynamical matrix with Fourier interpolation on the coarse q-mesh.
   write(message,"(a,i0,a)")"Will use Fourier interpolation to construct D(q) for ",size(qmissing)," q-points"
   call wrtout(std_out,message,"COLL")

   ABI_MALLOC(out_d2cart, (2,3,natom,3,natom))
   do ii=1,size(qmissing)
     iq_bz = qmissing(ii)
     qpt = spqpt(:,iq_bz)
     ! TODO: check dipdip option and phase, but I think this is correct!
     call ifc_fourq(Ifc_coarse,Crystal,qpt,phfrq,displ,out_d2cart=out_d2cart)
     Ifc%dynmat(:,:,:,:,:,iq_bz) = out_d2cart
   end do

   ABI_FREE(qmissing)
   ABI_FREE(out_d2cart)
 end if

!OmegaSRLR: Store full dynamical matrix for decomposition into short- and long-range parts
 ABI_MALLOC(dynmatfull,(2,3,natom,3,natom,nqbz))
 dynmatfull=Ifc%dynmat

 if (Ifc%dipdip==1) then
   ! Take off the dipole-dipole part of the dynamical matrix
   write(std_out, '(a)' )' mkifc9 : will extract the dipole-dipole part,'
   write(std_out, '(a)' )' using ewald9, q0dy3 and nanal9 for every wavevector.'
   ABI_MALLOC(dyew,(2,3,natom,3,natom))

   do iqpt=1,nqbz
     qpt(:)=spqpt(:,iqpt)
     sumg0=0
     call ewald9(ddb%acell,dielt,dyew,Crystal%gmet,gprim,natom,qpt,Crystal%rmet,rprim,sumg0,Crystal%ucvol,Crystal%xred,zeff)
     call q0dy3_apply(natom,dyewq0,dyew)

     plus=0
     call nanal9(dyew,Ifc%dynmat,iqpt,natom,nqbz,plus)
   end do

   ABI_FREE(dyew)
 end if

! OmegaSRLR: Store the short-range dynmat and compute long-range as difference
 ABI_MALLOC(dynmat_sr,(2,3,natom,3,natom,nqbz))
 ABI_MALLOC(dynmat_lr,(2,3,natom,3,natom,nqbz))
 dynmat_sr=Ifc%dynmat
 dynmat_lr=dynmatfull-dynmat_sr

! Now, take care of the remaining part of the dynamical matrix
! Move to canonical normalized coordinates
 call canat9(Ifc%brav,natom,rcan,rprim,trans,Crystal%xred)

! Multiply the dynamical matrix by a phase shift
 option=1
 call dymfz9(Ifc%dynmat,natom,nqbz,gprim,option,spqpt,trans)

! Create the Big Box of R vectors in real space and compute the number of points (cells) in real space
 call make_bigbox(Ifc%brav,ngqpt,nqshft,rprim,ifc_tmp%nrpt,ifc_tmp%rpt)

! Weights associated to these R points and to atomic pairs
! MG FIXME: Why ngqpt is intent(inout)?
 ABI_MALLOC(ifc_tmp%wghatm,(natom,natom,ifc_tmp%nrpt))
 call wght9(Ifc%brav,gprim,natom,ngqpt,nqbz,nqshft,ifc_tmp%nrpt,q1shft,rcan,ifc_tmp%rpt,ifc_tmp%wghatm)

! Fourier transformation of the dynamical matrices
 ABI_MALLOC(ifc_tmp%atmfrc,(2,3,natom,3,natom,ifc_tmp%nrpt))

 call ftifc_q2r(ifc_tmp%atmfrc,Ifc%dynmat,gprim,natom,nqbz,ifc_tmp%nrpt,ifc_tmp%rpt,spqpt)

! Eventually impose Acoustic Sum Rule to the interatomic forces
 if (Ifc%asr>0) then
   call asrif9(Ifc%asr,ifc_tmp%atmfrc,natom,ifc_tmp%nrpt,ifc_tmp%rpt,ifc_tmp%wghatm)
 end if

!*** The interatomic forces have been calculated ! ***
 write(message, '(2a)')ch10,' The interatomic forces have been obtained '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!Analysis of the real-space interatomic force constants
! MG: This is a post-processing tool. It should not be called here!
 if (anaddb_dtset%nsphere/=0 .or. anaddb_dtset%rifcsph>tol10 .or. anaddb_dtset%ifcout/=0) then
   call wrtout(std_out, 'mkifc9: analysis of real-space IFCs', "COLL")

   call timein(tcpui,twalli)

   call rsiaf9(ddb%acell,atifc,ifc_tmp%atmfrc,dielt,Ifc%dipdip,dyewq0,&
&   gprim,anaddb_dtset%ifcana,anaddb_dtset%ifcout,iout,&
&   natom,ifc_tmp%nrpt,anaddb_dtset%nsphere,anaddb_dtset%prt_ifc,&
&   rcan,anaddb_dtset%rifcsph,rprim,ifc_tmp%rpt,&
&   tcpui,twalli,ifc_tmp%wghatm,zeff)
 end if

 ABI_FREE(atifc)
 !ABI_FREE(dynmat)

!Eventually impose Acoustic Sum Rule to the interatomic forces
!(Note : here, after the analysis, in which the range
!of the short-range IFCs may have been changed
!That is why it is asked that asr be negative)
!FIXME: asr < 0 is not tested in abinit suite and does not appear to work
!(I get frequencies of 10^105 Hartree...) Modifying this 12/6/2011
!if(asr<0)then
!asr=-asr

 ! Be careful here: if I move this call to anaddb, several tests will fail
 ! due to the different number of points in the big box (see code below.)
 if (Ifc%asr > 0 .and. (anaddb_dtset%nsphere/=0 .or. anaddb_dtset%rifcsph>tol10)) then
   call asrif9(Ifc%asr,ifc_tmp%atmfrc,natom,ifc_tmp%nrpt,ifc_tmp%rpt,ifc_tmp%wghatm)
 end if

 ! Only conserve the necessary points in rpt: in the FT algorithm the order of the points is unimportant
 Ifc%nrpt = 0
 do irpt=1,ifc_tmp%nrpt
   if (sum(ifc_tmp%wghatm(:,:,irpt)) /= 0) then
     Ifc%nrpt = Ifc%nrpt+1
   end if
 end do

 ABI_MALLOC(Ifc%atmfrc,(2,3,natom,3,natom,Ifc%nrpt))
 ABI_MALLOC(Ifc%rpt,(3,Ifc%nrpt))
 ABI_MALLOC(Ifc%wghatm,(natom,natom,Ifc%nrpt))

 irpt_new = 1
 do irpt = 1, ifc_tmp%nrpt
   if (sum(ifc_tmp%wghatm(:,:,irpt)) /= 0) then
     Ifc%atmfrc(:,:,:,:,:,irpt_new) = ifc_tmp%atmfrc(:,:,:,:,:,irpt)
     Ifc%rpt(:,irpt_new) = ifc_tmp%rpt(:,irpt)
     Ifc%wghatm(:,:,irpt_new) = ifc_tmp%wghatm(:,:,irpt)
     irpt_new = irpt_new + 1
   end if
 end do

 ! Copy other useful arrays.
 Ifc%dielt = dielt
 Ifc%nqbz = nqbz

 call alloc_copy(rcan, Ifc%rcan)
 call alloc_copy(trans, Ifc%trans)
 call alloc_copy(dyewq0, Ifc%dyewq0)
 call alloc_copy(spqpt(:,1:nqbz), Ifc%qbz)
 call alloc_copy(zeff, Ifc%zeff)
 call alloc_copy(ddb%amu, Ifc%amu)

 call ifc_free(ifc_tmp)

 ! Check that the starting values are well reproduced.
 ! (This is to be suppressed in a future version)
 write(std_out, '(a,a)' )' mkifc9 : now check that the starting values ',&
& ' are reproduced after the use of interatomic forces '

 do iqpt=1,nqbz
   qpt(:)=Ifc%qbz(:,iqpt)

   call ifc_fourq(Ifc,Crystal,qpt,phfrq,displ,out_eigvec=eigvec)

   ! OmegaSRLR: Perform decomposition of dynamical matrix
   ! MG: FIXME I don't think the implementation is correct when q !=0
   if (anaddb_dtset%prtsrlr==1) then
     call omega_decomp(ddb%amu,natom,ntypat,Crystal%typat,dynmatfull,dynmat_sr,dynmat_lr,iqpt,nqbz,eigvec)
   end if

   ! Write the phonon frequencies (this is for checking purposes).
   ! Note: these phonon frequencies are not written on unit iout, only on unit std_out.
   call prtph3(displ,0,anaddb_dtset%enunit,-1,natom,phfrq,qphnrm,qpt)
 end do

 ABI_FREE(spqpt)

!OmegaSRLR: deallocate memory used by dynmat decomposition
 ABI_FREE(dynmatfull)
 ABI_FREE(dynmat_sr)
 ABI_FREE(dynmat_lr)

end subroutine ifc_init
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_fourq
!! NAME
!!  ifc_fourq
!!
!! FUNCTION
!!  Compute the phonon frequencies at the specified q-point by performing
!!  a Fourier transform on the IFCs matrix in real space.
!!
!! INPUTS
!!  Ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!!  Crystal<type(crystal_t)> = Information on the crystalline structure.
!!  qpt(3)=q-point in reduced coordinates (unless nanaqdir is specified)
!!  [nanaqdir]=If present, the qpt will be treated as a vector specifying the 
!!    direction in q-space along which the non-analytic behaviour of the dynamical
!!    matrix will be treated. Possible values:
!!       "cart" if qpt defines a direction in Cartesian coordinates
!!       "reduced" if qpt defines a direction in reduced coordinates
!!
!! OUTPUT
!!  phfrq(3*Crystal%natom)=Phonon frequencies in Hartree
!!  displ(2,3,Crystal%natom,3*Crystal%natom)=Phonon displacement in Cartesian coordinates
!!  [out_d2cart(2,3,3*natom,3,3*natom)]=Optional. The (interpolated) dynamical matrix for this q-point
!!  [out_eigvec(2*3*natom*3*natom)= Optional. The (interpolated) eigenvectors of the dynamical matrix.
!!
!! PARENTS
!!      get_nv_fs_en,get_tau_k,harmonic_thermo,interpolate_gkk,m_gamma,m_ifc
!!      m_phonons,mka2f,mka2f_tr,mka2f_tr_lova,mkph_linwid,outphbtrap,read_gkk
!!
!! CHILDREN
!!
!! SOURCE


subroutine ifc_fourq(Ifc,Crystal,qpt,phfrq,displ,nanaqdir,out_d2cart,out_eigvec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_fourq'
 use interfaces_72_response
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),optional,intent(in) :: nanaqdir
 type(ifc_type),intent(in) :: Ifc
 type(crystal_t),intent(in) :: Crystal
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: displ(2,3,Crystal%natom,3*Crystal%natom)
 real(dp),intent(out) :: phfrq(3*Crystal%natom)
 real(dp),optional,intent(out) :: out_d2cart(2,3,Crystal%natom,3,Crystal%natom)
 real(dp),optional,intent(out) :: out_eigvec(2,3,Crystal%natom,3*Crystal%natom)

!Local variables-------------------------------
!scalars
 integer :: natom
 real(dp) :: qphnrm
!arrays 
 real(dp) :: my_qpt(3),eigvec(2,3,Crystal%natom,3*Crystal%natom),eigval(3*Crystal%natom)
 real(dp) :: d2cart(2,3,Ifc%mpert,3,Ifc%mpert)

! ************************************************************************

 natom = Crystal%natom

 ! Use my_qpt because phfrq3 can change the q-point (very bad design)
 qphnrm = one; my_qpt = qpt

 if (present(nanaqdir)) then
   ! This will break backward compatibility because qpt is **always** in reduced coordinates.
   ! while phfrq3 assume cartesian coordinates !!!!!!!!!!!
   ! It does not make sense to change API just to treat this particular case
   ! We should **alwayse use q-points in reduced coordinates.
   qphnrm = zero
   select case (nanaqdir)
   case ("reduced")
     ! Convert to Cartesian.
     my_qpt = matmul(Crystal%gprimd, qpt) 
   case ("cart")
     continue
   case default
     MSG_ERROR("Wrong value for nanaqdir: "//trim(nanaqdir))
   end select
 end if

 ! The dynamical matrix d2cart is calculated here:
 call gtdyn9(Ifc%acell,Ifc%atmfrc,Ifc%dielt,Ifc%dipdip,Ifc%dyewq0,d2cart,Crystal%gmet,Ifc%gprim,Ifc%mpert,natom,&
&  Ifc%nrpt,qphnrm,my_qpt,Crystal%rmet,Ifc%rprim,Ifc%rpt,Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,Ifc%zeff)

 ! Calculate the eigenvectors and eigenvalues of the dynamical matrix
 call phfrq3(Ifc%amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
&  Ifc%mpert,Crystal%nsym,natom,Crystal%nsym,Crystal%ntypat,phfrq,qphnrm,my_qpt,&
&  Crystal%rprimd,Ifc%symdynmat,Crystal%symrel,Crystal%typat,Crystal%ucvol)

 ! OmegaSRLR: Perform decomposition of dynamical matrix
 !if (srlr==1)then
 !  call omega_decomp(amu,natom,ntypat,typat,dynmatfull,dynmatsr,dynmatlr,iqpt,nqpt,eigvec)
 !end if

 ! Returns the interpolated dynamical matrix and the eigenvector for this q-point
 if (present(out_d2cart)) out_d2cart = d2cart(:,:,:natom,:,:natom)
 if (present(out_eigvec)) out_eigvec = eigvec

 ! Option to get vectors in reduced coordinates?
 !call phdispl_cart2red(natom,Crystal%gprimd,displ_cart,displ_red)
 !call phdispl_cart2red(natom,Crystal%gprimd,out_eigvec,displ_red)

end subroutine ifc_fourq
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/ifc_diagoq
!! NAME
!!  ifc_diagoq
!!
!! FUNCTION
!!  Compute the phonon frequencies at the specified q-point by performing
!!  a direct diagonalizatin of the dynamical matrix. The q-point **MUST** be
!!  one the points used for the DFPT calculation or one of its symmetrical image.
!!
!! INPUTS
!!  Ifc<type(ifc_type)>=Object containing the dynamical matrix and the IFCs.
!!  Crystal<type(crystal_t)> = Information on the crystalline structure.
!!  qpt(3)=q-point in reduced coordinates (unless nanaqdir is specified)
!!  [nanaqdir]=If present, the qpt will be treated as a vector specifying the 
!!    direction in q-space along which the non-analytic behaviour of the dynamical
!!    matrix will be treated. Possible values:
!!       "cart" if qpt defines a direction in Cartesian coordinates
!!       "reduced" if qpt defines a direction in reduced coordinates
!!
!! OUTPUT
!!  phfrq(3*Crystal%natom)=Phonon frequencies in Hartree
!!  displ(2,3*Crystal%natom,3*Crystal%natom)=Phonon displacement in Cartesian coordinates
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if 0

subroutine ifc_diagoq(Ifc,Crystal,qpt,phfrq,displ,nanaqdir)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_diagoq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),optional,intent(in) :: nanaqdir
 type(ifc_type),intent(in) :: Ifc
 type(crystal_t),intent(in) :: Crystal
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: phfrq(3*Crystal%natom)
 real(dp),intent(out) :: displ(2,3,Crystal%natom,3,Crystal%natom)

!Local variables-------------------------------
!scalars
 real(dp) :: qphnrm
!arrays 
 real(dp) :: my_qpt(3) !,eigvec(2,3*Crystal%natom,3*Crystal%natom),eigval(3*Crystal%natom)
 !real(dp) :: d2cart(2,3,Ifc%mpert,3,Ifc%mpert)

! ************************************************************************
 MSG_ERROR("Not implemented error")

 ! Use my_qpt because phfrq3 can change the q-point (very bad design)
 qphnrm = one; my_qpt = qpt
                                                                                              
 if (present(nanaqdir)) then
   ! This will break backward compatibility because qpt is **always** in reduced coordinates.
   ! while phfrq3 assume cartesian coordinates !!!!!!!!!!!
   ! It does not make sense to change API just to treat this particular case
   ! We should **alwayse use q-points in reduced coordinates.
   qphnrm = zero
   select case (nanaqdir)
   case ("reduced")
     ! Convert to Cartesian.
     my_qpt = matmul(Crystal%gprimd, qpt) 
   case ("cart")
     continue
   case default
     MSG_ERROR("Wrong value for nanaqdir: "//trim(nanaqdir))
   end select
 end if

 ! See mkphbs
 ! Copy the dynamical matrix in d2cart
 !d2cart(:,1:msize)=ddb%val(:,:,iblok)

 ! Eventually impose the acoustic sum rule based on previously calculated d2asr
 !if (anaddb_dtset%asr==1 .or. anaddb_dtset%asr==2 .or. anaddb_dtset%asr==5) then
 !  call asria_corr(anaddb_dtset%asr,d2asr,d2cart,mpert,natom)
 !end if

 ! Impose acoustic sum rule plus rotational symmetry for 0D and 1D systems
 !if (anaddb_dtset%asr==3 .or. anaddb_dtset%asr==4) then
 !  call asrprs(anaddb_dtset%asr,2,3,uinvers,vtinvers,singular,d2cart,mpert,natom,xcart)
 !end if

 ! Calculate the eigenvectors and eigenvalues of the dynamical matrix
! call phfrq3(amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
!&  mpert,Cryst%nsym,Crystal%natom,Crystal%nsym,Crystal%ntypat,phfrq,qphnrm,qpt,&
!&  Crystal%rprimd,symdynmat,Crystal%symrel,Crystal%typat,Crystal%ucvol)

end subroutine ifc_diagoq
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_ifc/rsiaf9
!!
!! NAME
!! rsiaf9
!!
!! FUNCTION
!! Compute the real-space interatomic force constants, including
!! both analytical (short-range) and non-analytical (long-range contribution)
!!
!! INPUTS
!! acell(3)=length scales by which rprim is to be multiplied
!! atifc(natom) =  atifc(ia) equals 1 if the analysis of ifc
!!  has to be done for atom ia; otherwise 0.
!! atmfrc(2,3,natom,3,natom,nrpt)
!!  = Analytical part of the Interatomic Forces in real space.
!!  We used the imaginary part just for debugging
!! dielt(3,3)=dielectric tensor
!! dipdip= if 0, no dipole-dipole interaction was subtracted in atmfrc
!!  if 1, atmfrc has been build without dipole-dipole part
!! dyewq0(3,3,natom)=contraction of the Ewald dynamical matrix at q=0
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! ifcana= 0 => no analysis of ifc ; 1 => full analysis
!! ifcout= Number of interatomic force constants written in the output file
!! iout=unit number for nice output
!! natom=number of atoms in unit cell
!! nrpt= Number of R points in the Big Box
!! nsphere=number of atoms to be included in the cut-off sphere for interatomic
!!  force constant; if = 0 : maximum extent allowed by the grid.
!! prt_ifc = flag to print out ifc information for dynamical matrix (AI2PS) 
!! rcan(3,natom)=canonical coordinates of atoms
!! rifcsph=radius for cutoff of IFC
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nrpt)=canonical coordinates of the points in the BigBox.
!! tcpui,twalli=initial values of cpu and wall clocktime
!! wghatm(natom,natom,nrpt) = Weights associated to a pair of atoms and to a R vector
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!!
!! OUTPUT
!!   written in the output file.
!!
!! NOTES
!! This routine should be executed by one processor only
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine rsiaf9(acell,atifc,atmfrc,dielt,dipdip,dyewq0,&
& gprim,ifcana,ifcout,iout,natom,nrpt,nsphere,prt_ifc,rcan,&
& rifcsph,rprim,rpt,tcpui,twalli,wghatm,zeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rsiaf9'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: dipdip,ifcana,ifcout,iout,natom,nrpt,nsphere
 real(dp),intent(in) :: rifcsph,tcpui,twalli
 integer, intent(in) :: prt_ifc
!arrays
 integer,intent(in) :: atifc(natom)
 real(dp),intent(in) :: acell(3),atmfrc(2,3,natom,3,natom,nrpt),dielt(3,3)
 real(dp),intent(in) :: dyewq0(3,3,natom),gprim(3,3),rcan(3,natom)
 real(dp),intent(in) :: rprim(3,3),rpt(3,nrpt),zeff(3,3,natom)
 real(dp),intent(inout) :: wghatm(natom,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: flag,ia,ib,ii,index,irpt,jj,kk,mu,nu
! unit number to print out ifc information for dynamical matrix (AI2PS) 
 integer :: unit_ifc
 real(dp) :: detdlt,dist1,dist2,rsq,scprod,tcpu,trace1,trace2,trace3
 real(dp) :: twall,yy
 character(len=500) :: message
!arrays
 integer,allocatable :: list(:),sorted(:)
 real(dp) :: ewiaf0(3,3),ewiaf1(3,3),ewloc(3,3),ifcloc(3,3),invdlt(3,3),ra(3)
 real(dp) :: rbpos(3),rcart(3),rdiff(3),rsiaf(3,3),rsloc(3,3),sriaf(3,3)
 real(dp) :: srloc(3,3),vect1(3),vect2(3),vect3(3),work(3),xred(3),xx(3)
 real(dp),allocatable :: dist(:,:,:),wkdist(:)

! *********************************************************************

 ! Calculating the inverse (transpose) of the dielectric tensor
 call matr3inv(dielt,invdlt)

 ! Calculating the determinant of the dielectric tensor
 detdlt=dielt(1,1)*dielt(2,2)*dielt(3,3)+dielt(1,3)*dielt(2,1)*&
& dielt(3,2)+dielt(1,2)*dielt(2,3)*dielt(3,1)-dielt(1,3)*&
& dielt(2,2)*dielt(3,1)-dielt(1,1)*dielt(2,3)*dielt(3,2)-&
& dielt(1,2)*dielt(2,1)*dielt(3,3)

 ! Compute the distances between atoms
 ABI_MALLOC(dist,(natom,natom,nrpt))
 call dist9(acell,dist,gprim,natom,nrpt,rcan,rprim,rpt)
 ! Now dist(ia,ib,irpt) contains the distance from atom ia to atom ib in unit cell irpt.

 write(std_out,'(a)' )' rsiaf9 : analysis of interatomic force constants '
 if (iout > 0) then
   write(iout, '(/,a,/)' )' Analysis of interatomic force constants '
   if(dipdip==1)then
     write(iout, '(a)' )' Are given : column(1-3), the total force constant'
     write(iout, '(a)' )'       then  column(4-6), the Ewald part'
     write(iout, '(a)' )'       then  column(7-9), the short-range part'
     write(iout, '(a)' )' Column 1, 4 and 7 are related to the displacement'
     write(iout, '(a)' )'       of the generic atom along x,               '
     write(iout, '(a)' )' column 2, 5 and 8 are related to the displacement'
     write(iout, '(a)' )'       of the generic atom along y,               '
     write(iout, '(a)' )' column 3, 6 and 9 are related to the displacement'
     write(iout, '(a)')'       of the generic atom along z.               '
   else if(dipdip==0)then
     write(iout, '(a)' )' column 1 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along x,    '
     write(iout, '(a)' )' column 2 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along y,    '
     write(iout, '(a)' )' column 3 is related to the displacement'
     write(iout, '(a)' )'        of the generic atom along z,    '
   end if
 end if

 ABI_MALLOC(list,(natom*nrpt))
 ABI_MALLOC(sorted,(natom*nrpt))

 ! set up file for real space ifc output, if required
 if (prt_ifc > 0) then
   unit_ifc = get_unit()
   open(unit_ifc, file='ifcinfo.out', status="replace")
   write(iout, '(a,a)' )ch10,&
&   '  NOTE : Open file ifcinfo.out, for the output of interatomic force constants. This is because prt_ifc==1. '
 end if

 ! BIG loop on all generic atoms
 do ia=1,natom

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-before big ia loop at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')

   ! First transform canonical coordinates to reduced coordinates
   do ii=1,3
     xred(ii)=gprim(1,ii)*rcan(1,ia)+gprim(2,ii)*rcan(2,ia)+gprim(3,ii)*rcan(3,ia)
   end do

   ! Then to cartesian coordinates
   ra(:)=xred(1)*acell(1)*rprim(:,1)+ xred(2)*acell(2)*rprim(:,2)+ xred(3)*acell(3)*rprim(:,3)

   ! write(std_out,*) ' nsphere, rifcsph, atifc(ia) =', nsphere, rifcsph, atifc(ia)

   if( (nsphere/=0.and.nsphere<natom*nrpt) .or. rifcsph > tol10 .or. atifc(ia)==1)then
     ! This sorting algorithm is slow ...
     ABI_MALLOC(wkdist,(natom*nrpt))
     wkdist(:)=reshape(dist(ia,:,:),(/natom*nrpt/))
     do ii=1,natom*nrpt
       list(ii)=ii
     end do
     call sort_dp(natom*nrpt,wkdist,list,tol14)
   end if

   ! In case some cut-off has to be imposed on the IFCs,
   ! zero the outside IFCs now : act on wghatm
   if(nsphere/=0.and.nsphere<natom*nrpt)then
     do ii=nsphere+1,natom*nrpt
       index=list(ii)
       irpt=(index-1)/natom+1
       ib=index-natom*(irpt-1)
       wghatm(ia,ib,irpt)=0._dp
     end do
   end if

   if(rifcsph>tol10)then
     ! write(std_out,*) ' rsiaf9 : wkdist = ',wkdist
     do ii=nsphere+1,natom*nrpt
       index=list(ii)
       ! preserve weights for atoms inside sphere of radius rifcsph
       if (wkdist(ii) < rifcsph) cycle
       irpt=(index-1)/natom+1
       ib=index-natom*(irpt-1)
       wghatm(ia,ib,irpt)=0._dp
     end do
   end if

   ! deallocate wkdist if used
   if(allocated(wkdist))  then
     ABI_FREE(wkdist)
   end if

   if(atifc(ia)==1)then

     if (iout > 0) then
       write(iout, '(a)' )
       write(std_out,'(a,i4)' )' generic atom number',ia
       write(iout, '(a,i4)' )' generic atom number',ia
       write(std_out,'(a,3es16.8)' ) ' with cartesian coordinates',ra(1:3)
       write(iout,'(a,3es16.8)' ) ' with cartesian coordinates',ra(1:3)
       write(iout, '(a)' )
     end if

     if(ifcana==1)then
       ! Generate the local coordinate system for the atom ia
       index=list(2)
       write(std_out,*)index
       call canct9(acell,gprim,ib,index,irpt,natom,nrpt,rcan,rcart,rprim,rpt)
       dist1=dist(ia,ib,irpt)
       vect2(1)=rcart(1)-ra(1)
       vect2(2)=rcart(2)-ra(2)
       vect2(3)=rcart(3)-ra(3)
       flag=0
       do ii=3,natom*nrpt
         index=list(ii)
         call canct9(acell,gprim,ib,index,irpt,natom,nrpt,rcan,rcart,rprim,rpt)
         dist2=dist(ia,ib,irpt)
         vect1(1)=(rcart(1)-ra(1))-vect2(1)
         vect1(2)=(rcart(2)-ra(2))-vect2(2)
         vect1(3)=(rcart(3)-ra(3))-vect2(3)
         scprod=0.0_dp
         do jj=1,3
           scprod=scprod+vect1(jj)**2
         end do
         do jj=1,3
           vect1(jj)=vect1(jj)/scprod**0.5
         end do
         scprod=0.0_dp
         do jj=1,3
           scprod=scprod+vect2(jj)*vect1(jj)
         end do
         do jj=1,3
           work(jj)=vect2(jj)-vect1(jj)*scprod
         end do
         scprod=0.0_dp
         do jj=1,3
           scprod=scprod+work(jj)**2
         end do
         if(scprod>1.0d-10)then
           flag=1
         end if
         if(flag==1)exit
       end do
       if(flag==0)then
         write(message, '(3a)' )&
&         'Unable to find a third atom not aligned',ch10,&
&         'with the two selected ones.'
         MSG_BUG(message)
       end if
       vect2(1)=work(1)/scprod**0.5
       vect2(2)=work(2)/scprod**0.5
       vect2(3)=work(3)/scprod**0.5
       vect3(1)=vect1(2)*vect2(3)-vect1(3)*vect2(2)
       vect3(2)=vect1(3)*vect2(1)-vect1(1)*vect2(3)
       vect3(3)=vect1(1)*vect2(2)-vect1(2)*vect2(1)
       if (iout > 0) then
         write(iout, '(a)' )' Third atom defining local coordinates : '
         write(iout, '(a,i4,a,i4)' )'     ib = ',ib,'   irpt = ',irpt
       end if
     end if

     ! Analysis and output of force constants, ordered with respect to the distance from atom ia
     do ii=1,ifcout
       index=list(ii)
       call canct9(acell,gprim,ib,index,irpt,natom,nrpt,rcan,rbpos,rprim,rpt)
       if (iout > 0) then
         write(iout, '(a)' )
         write(iout, '(i4,a,i6,a,i8)' )ii,' interaction with atom',ib,' cell',irpt
         write(iout, '(a,3es16.6)' )' with coordinates ',rbpos(1:3)*(one+tol8)
         write(iout, '(a,es16.6)' )' and distance ',dist(ia,ib,irpt)
       end if

       if(ifcana==1.and.ii/=1)then
         dist1=dist(ia,ib,irpt)
         vect1(1)=(rbpos(1)-ra(1))/dist1
         vect1(2)=(rbpos(2)-ra(2))/dist1
         vect1(3)=(rbpos(3)-ra(3))/dist1
       end if
       
       if (prt_ifc == 1) then
         write(unit_ifc,'(i6,i6)') ia,ii
         write(unit_ifc,'(3es28.16)') rbpos(1:3)
       end if

       if(dipdip==0)then

         ! Get the "total" force constants (=real space FC)
         ! without taking into account the dipole-dipole interaction
         do mu=1,3
           do nu=1,3
             rsiaf(mu,nu)=atmfrc(1,mu,ia,nu,ib,irpt) * wghatm(ia,ib,irpt)
           end do
         end do

         ! Output of the ifcs in cartesian coordinates
         do nu=1,3
           if (iout > 0) then
             write(iout, '(1x,3f9.5)' )(rsiaf(mu,nu)+tol10,mu=1,3)
           end if
           if (prt_ifc == 1) then
             write(unit_ifc,   '(3f28.16)'  )(rsiaf(nu,mu),mu=1,3)
           end if
         end do

         if(ifcana==1)then
           ! Further analysis
           trace1=rsiaf(1,1)+rsiaf(2,2)+rsiaf(3,3)
           if (iout > 0) then
             write(iout, '(a,f9.5)' ) '  Trace         ',trace1+tol10
           end if
           if(ii/=1)then
             call axial9(rsiaf,vect1,vect2,vect3)
           end if
           if (iout > 0) then
             write(iout, '(a)' )' Transformation to local coordinates '
             write(iout, '(a,3f16.6)' ) ' First  local vector :',vect1
             write(iout, '(a,3f16.6)' ) ' Second local vector :',vect2
             write(iout, '(a,3f16.6)' ) ' Third  local vector :',vect3
           end if
           call ifclo9(rsiaf,ifcloc,vect1,vect2,vect3)
           if (iout > 0) then
             do nu=1,3
               write(iout, '(1x,3f9.5)' )(ifcloc(mu,nu)+tol10,mu=1,3)
             end do
           end if
         end if ! Further analysis finished

       else if(dipdip==1)then

         ! Get the Coulomb part
         do jj=1,3
           rdiff(jj)=ra(jj)-rbpos(jj)
         end do
         rsq=0.0_dp
         xx(1:3)=0.0_dp
         do jj=1,3
           do kk=1,3
             ewiaf0(jj,kk)=0.0_dp
             rsq=rsq+rdiff(jj)*invdlt(kk,jj)*rdiff(kk)
             xx(kk)=xx(kk)+invdlt(kk,jj)*rdiff(jj)
           end do
         end do
         yy=sqrt(rsq)
         !  Avoid zero denominators in term:
         if (sqrt(rsq)>=tol12) then
           do mu=1,3
             do nu=1,3
               ewiaf0(mu,nu)=(-3*xx(nu)*xx(mu)+invdlt(nu,mu)*yy**2)/yy**5/sqrt(detdlt)
             end do
           end do
         else
           if (ia/=ib)then
             write(message, '(a,a,a,a,a,i5,a,i5,a)' )&
&             'The distance between two atoms vanishes.',ch10,&
&             'This is not allowed.',ch10,&
&             'Action: check the input for the atoms number',ia,' and',ib,'.'
             MSG_ERROR(message)
           end if
         end if

         ! Take into account the effective charge tensor
         do mu=1,3
           do nu=1,3
             ewiaf1(mu,nu)=0.0_dp
             if(ii==1)then
               ewiaf1(mu,nu)=-dyewq0(mu,nu,ia)
             end if
             do jj=1,3
               do kk=1,3
                 ewiaf1(mu,nu)=ewiaf1(mu,nu) +zeff(jj,mu,ia)*zeff(kk,nu,ib)* ewiaf0(jj,kk)
               end do
             end do
           end do
         end do

         ! Get the short-range force constants and the
         ! "total" force constants (=real space FC)
         do mu=1,3
           do nu=1,3
             sriaf(mu,nu)=atmfrc(1,mu,ia,nu,ib,irpt) * wghatm(ia,ib,irpt)
             rsiaf(mu,nu)=ewiaf1(mu,nu)+sriaf(mu,nu)
           end do
         end do

         ! Output of the results
         do nu=1,3
           if (iout > 0) then
             write(iout, '(1x,3(3f9.5,1x))' )&
&             (rsiaf(mu,nu) +tol10,mu=1,3),&
&             (ewiaf1(mu,nu)+tol10,mu=1,3),&
&             (sriaf(mu,nu) +tol10,mu=1,3)
           end if
           if (prt_ifc == 1) then
             write(unit_ifc,  '(3f28.16)'  )(rsiaf(nu,mu),mu=1,3)
           end if
         end do

         if(ifcana==1)then
           ! Further analysis
           if (iout > 0) then
             write(iout, '(a)' )' Traces (and ratios) :'
           end if
           trace1=rsiaf(1,1)+rsiaf(2,2)+rsiaf(3,3)
           trace2=ewiaf1(1,1)+ewiaf1(2,2)+ewiaf1(3,3)
           trace3=sriaf(1,1)+sriaf(2,2)+sriaf(3,3)
           if (iout > 0) then
             write(iout,'(3(f9.5,17x))')trace1+tol10,trace2+tol10,trace3+tol10
             write(iout,'(3(f9.5,17x))')1.0,trace2/trace1+tol10,trace3/trace1+tol10
           end if

           if(ii/=1)then
             call axial9(rsiaf,vect1,vect2,vect3)
           end if
           if (iout > 0) then
             write(iout, '(a)' )' Transformation to local coordinates '
             write(iout, '(a,3f16.6)' )' First  local vector :',vect1
             write(iout, '(a,3f16.6)' )' Second local vector :',vect2
             write(iout, '(a,3f16.6)' )' Third  local vector :',vect3
           end if
           call ifclo9(rsiaf,rsloc,vect1,vect2,vect3)
           call ifclo9(ewiaf1,ewloc,vect1,vect2,vect3)
           call ifclo9(sriaf,srloc,vect1,vect2,vect3)
           if (iout > 0) then
             do nu=1,3
               write(iout, '(1x,3(3f9.5,1x))' )&
&               (rsloc(mu,nu)+tol10,mu=1,3),&
&               (ewloc(mu,nu)+tol10,mu=1,3),&
&               (srloc(mu,nu)+tol10,mu=1,3)
             end do
             if(ii/=1)then
               write(iout, '(a)' )' Ratio with respect to the longitudinal ifc'
             else
               write(iout, '(a)' )' Ratio with respect to the (1,1) element'
             end if
             do nu=1,3
               write(iout, '(1x,3(3f9.5,1x))' )&
&               (rsloc(mu,nu)/rsloc(1,1)+tol10,mu=1,3),&
&               (ewloc(mu,nu)/rsloc(1,1)+tol10,mu=1,3),&
&               (srloc(mu,nu)/rsloc(1,1)+tol10,mu=1,3)
             end do
           end if

         end if ! Further analysis finished
       end if ! End the condition on dipdip
     end do ! End loop over all atoms in BigBox:

   end if
 end do ! End Big loop on atoms in the unit cell, and corresponding test .

 ABI_FREE(dist)
 ABI_FREE(list)
 ABI_FREE(sorted)

 if (prt_ifc > 0) close(unit_ifc)

end subroutine rsiaf9
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/omega_decomp
!! NAME
!!  omega_decomp
!!
!! FUNCTION
!! Compute and return the eigenvalues (frequencies) of the short-range and 
!! long-range part of the dynamical matrix  See Europhys. Lett. 33 p.713 (1996) for details.
!! (included by U. Aschauer and EB)
!!  
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine omega_decomp(amu,natom,ntypat,typat,dynmatfl,dynmatsr,dynmatlr,iqpt,nqpt,eigenvec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'omega_decomp'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 integer,intent(in) :: iqpt,nqpt
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat)
 real(dp),intent(inout) :: dynmatfl(2,3,natom,3,natom,nqpt)
 real(dp),intent(inout) :: dynmatsr(2,3,natom,3,natom,nqpt)
 real(dp),intent(inout) :: dynmatlr(2,3,natom,3,natom,nqpt)
 real(dp),intent(in)    :: eigenvec(2*3*natom*3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,i2,idir1,idir2,imode,ipert1,ipert2,index1,index2
 real(dp),parameter :: break_symm=1.0d-12
 real(dp) :: fac
!arrays
 real(dp) :: nearidentity(3,3)
 real(dp) :: omegafl, omegasr, omegalr
 real(dp) :: sumfl,sumlr,sumsr,asr
! *********************************************************************

!write(ab_out,*)''
!write(std_out,*) 'SR/LR decomposition: enter for wavevector number :',iqpt
 
!apply asr (note the lr part is already asred by construction in mkifc9)
 do ipert1=1,natom
   do idir1=1,3
     do idir2=1,3
       asr=0.0d0
       do ipert2=1,natom
         asr=asr+dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)
       end do
       dynmatfl(1,idir1,ipert1,idir2,ipert1,iqpt)=dynmatfl(1,idir1,ipert1,idir2,ipert1,iqpt)-asr
       dynmatsr(1,idir1,ipert1,idir2,ipert1,iqpt)=dynmatsr(1,idir1,ipert1,idir2,ipert1,iqpt)-asr
     end do
   end do
 end do


!This slight breaking of the symmetry allows the
!results to be more portable between machines
 nearidentity(:,:)=1.0
 nearidentity(1,1)=1.0+break_symm
 nearidentity(3,3)=1.0-break_symm


!Include Mass
 do ipert1=1,natom
   do ipert2=1,natom
     
     fac=1.0d0/sqrt(amu(typat(ipert1))*amu(typat(ipert2)))/amu_emass
     
     do idir1=1,3
       do idir2=1,3
         
         dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&         dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)*&
&         fac*nearidentity(idir1,idir2)
         
         dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&         dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)*&
&         fac*nearidentity(idir1,idir2)

         dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&         dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)*&
&         fac*nearidentity(idir1,idir2)

         ! This is to break slightly the translation invariance, and make
         ! the automatic tests more portable
         if(ipert1==ipert2 .and. idir1==idir2)then
           dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&           dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)+&
&           break_symm*natom/amu_emass/idir1*0.01d0

           dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&           dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)+&
&           break_symm*natom/amu_emass/idir1*0.01d0

           dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&           dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)+&
&           break_symm*natom/amu_emass/idir1*0.01d0
         end if

       end do
     end do
   end do
 end do


!Calculation of <eigvec|Dyn_tot,Dyn_SR,Dyn_LR|eigenvec>=omega**2

!write(ab_out,*)''
!write(ab_out,*)'==============================================================================='
 write(ab_out,*)''
 write(ab_out,*) 'Long-Range/Short-Range decomposed phonon freq. (cm-1)**2'
 write(ab_out,*) 'at wavevector number:',iqpt
 write(ab_out,*)''
 write(ab_out,'(a13,1x,a16,2x,a16,2x,a16)') ' Mode number.','tot**2','SR**2','LR**2'
 write(std_out,'(a13,1x,a16,2x,a16,2x,a16)') ' Mode number.','tot**2','SR**2','LR**2'
!write(ab_out,'(a12,2x,a10,2x,a10,2x,a10,2x,a16,2x,a16,2x,a16)') 'Mode number.','tot','SR','LR','tot**2','SR**2','LR**2'
!write(std_out,'(a12,2x,a10,2x,a10,2x,a10,2x,a16,2x,a16,2x,a16)') 'Mode number.','tot','SR','LR','tot**2','SR**2','LR**2'

 do imode=1,3*natom
   
   sumfl=0.0d0
   sumlr=0.0d0
   sumsr=0.0d0
   
   do ipert1=1,natom
     do ipert2=1,natom
       do i1=1,3
         do i2=1,3
           
           index1=i1+(ipert1-1)*3+3*natom*(imode-1)
           index2=i2+(ipert2-1)*3+3*natom*(imode-1)
           ! MG FIXME: I don't think these expressions are correct when q != 0
           ! We should also include the imaginary part
           
           sumfl = sumfl + eigenvec(2*index1-1) * dynmatfl(1,i1,ipert1,i2,ipert2,iqpt) * eigenvec(2*index2-1)
           sumlr = sumlr + eigenvec(2*index1-1) * dynmatlr(1,i1,ipert1,i2,ipert2,iqpt) * eigenvec(2*index2-1)
           sumsr = sumsr + eigenvec(2*index1-1) * dynmatsr(1,i1,ipert1,i2,ipert2,iqpt) * eigenvec(2*index2-1)
         end do
       end do
     end do
   end do
   
   sumfl = sumfl * Ha_cmm1 * Ha_cmm1
   sumsr = sumsr * Ha_cmm1 * Ha_cmm1
   sumlr = sumlr * Ha_cmm1 * Ha_cmm1
   
!  Compute omega=sqrt(omega**2)
   if(sumfl>=1.0d-16)then
     omegafl=sqrt(sumfl)
   else if(sumfl>=-1.0d-16)then
     omegafl=zero
   else
     omegafl=-sqrt(-sumfl)
   end if

   if(sumsr>=1.0d-16)then
     omegasr=sqrt(sumsr)
   else if(sumsr>=-1.0d-16)then
     omegasr=zero
   else
     omegasr=-sqrt(-sumsr)
   end if

   if(sumlr>=1.0d-16)then
     omegalr=sqrt(sumlr)
   else if(sumlr>=-1.0d-16)then
     omegalr=zero
   else
     omegalr=-sqrt(-sumlr)
   end if

!  Output
   write(ab_out,'(i4,10x,s,f16.4,2x,f16.4,2x,f16.4)') imode,sumfl,sumsr,sumlr  !vz_d
   write(std_out,'(i4,10x,s,f16.4,2x,f16.4,2x,f16.4)') imode,sumfl,sumsr,sumlr  !vz_d
!  write(ab_out,'(i4,8x,f10.4,2x,f10.4,2x,f10.4,2x,s,f16.6,2x,f16.6,2x,f16.6)') imode,omegafl,omegasr,omegalr,sumfl,sumsr,sumlr
!  write(std_out,'(i4,8x,f10.4,2x,f10.4,2x,f10.4,2x,s,f16.6,2x,f16.6,2x,f16.6)') imode,omegafl,omegasr,omegalr,sumfl,sumsr,sumlr
 end do

end subroutine omega_decomp
!!***

!----------------------------------------------------------------------

END MODULE m_ifc
