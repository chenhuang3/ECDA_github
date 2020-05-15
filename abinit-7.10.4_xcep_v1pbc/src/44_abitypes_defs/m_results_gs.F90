!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_results_gs
!! NAME
!!  m_results_gs
!!
!! FUNCTION
!!  This module provides the definition of the results_gs_type
!!  used to store results from GS calculations.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2014 ABINIT group (MT)
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

MODULE m_results_gs

 use defs_basis
 use m_profiling_abi
 use m_energies
 use m_errors
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 implicit none

 private
!!***

!!****t* m_results_gs/results_gs_type
!! NAME
!! results_gs_type
!!
!! FUNCTION
!! This structured datatype contains the results of a GS calculation :
!! energy and its decomposition, forces and their decompositions, stresses
!! and their decompositions
!!
!! SOURCE

 type, public :: results_gs_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: natom
   ! The number of atoms for this dataset

  integer :: nsppol
   ! The number of spin channels for this dataset

  integer :: ngrvdw
   ! Size of grvdw array
   ! Can be 0 (not allocated) or natom

! Real (real(dp)) scalars

  real(dp) :: deltae
   ! change in energy (Hartree)

  real(dp) :: diffor
   ! maximal absolute value of changes in the components of force

! All the energies are in Hartree, obtained "per unit cell".
  type(energies_type) :: energies
!!!  real(dp) :: eei      ! local pseudopotential energy (Hartree)
!!!  real(dp) :: eeig     ! sum of eigenvalue energy (Hartree)
!!!  real(dp) :: eew      ! Ewald energy (Hartree)
!!!  real(dp) :: ehart    ! Hartree part of total energy (Hartree)
!!!  real(dp) :: eii      ! pseudopotential core-core energy
!!!  real(dp) :: ek       ! kinetic energy (Hartree)
!!!  real(dp) :: enefield ! the term of the energy functional that depends
!!!                       ! explicitely on the electric field
!!!                       ! enefield = -ucvol*E*P
!!!  real(dp) :: enl      ! nonlocal pseudopotential energy (Hartree)
  real(dp) :: entropy  ! entropy (Hartree)
!!!  real(dp) :: enxc     ! exchange-correlation energy (Hartree)
!!!  real(dp) :: enxcdc   ! exchange-correlation double-counting energy (Hartree)
!!!  real(dp) :: epaw     ! PAW spherical energy (Hartree)
!!!  real(dp) :: epawdc   ! PAW spherical double-counting energy (Hartree)
  real(dp) :: etotal   ! total energy (Hartree)
                       ! for fixed occupation numbers (occopt==0,1,or 2):
                       !   etotal=ek+ehart+enxc+eei+eew+eii+enl+PAW_spherical_part
                       ! for varying occupation numbers (occopt>=3):
                       !   etotal=ek+ehart+enxc+eei+eew+eii+enl - tsmear*entropy +PAW_spherical_part
  real(dp) :: fermie   ! Fermi energy (Hartree)
  real(dp) :: residm   ! maximum value for the residual over all bands, all k points,
                       !   and all spins (Hartree or Hartree**2, to be checked !)
  real(dp) :: res2     ! density/potential residual (squared)
  real(dp) :: vxcavg   ! Average of the exchange-correlation energy. The average
                       ! of the local psp pot and the Hartree pot is set to zero (due
                       ! to the usual problem at G=0 for Coulombic system, so vxcavg
                       ! is also the average of the local part of the Hamiltonian

! Real (real(dp)) arrays

  real(dp), allocatable :: fcart(:,:)
   ! fcart(3,natom)
   ! Cartesian forces (Hartree/Bohr)
   ! Note : unlike fred, this array has been corrected by enforcing
   ! the translational symmetry, namely that the sum of force
   ! on all atoms is zero.

  real(dp), allocatable :: fred(:,:)
   ! fred(3,natom)
   ! Forces in reduced coordinates (Hartree)
   ! Actually, gradient of the total energy with respect
   ! to change of reduced coordinates

  real(dp),allocatable :: gaps(:,:)
   ! gaps(3,nsppol)
   ! gaps(1,:) : fundamental gap
   ! gaps(2,:) : optical gap
   ! gaps(3,:) : "status" for each channel : 0.0dp if the gap was not computed
   !   (because there are only valence bands) ; -1.0dp if the system (or spin-channel) is metallic ; 1.0dp if the
   !   gap was computed

  real(dp), allocatable :: gresid(:,:)
   ! gresid(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the residual
   ! of the potential

  real(dp), allocatable :: grewtn(:,:)
   ! grewtn(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the Ewald energy

  real(dp), allocatable :: grvdw(:,:)
   ! grvdw(3,ngrvdw)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from
   ! Van der Waals DFT-D2 dispersion (hartree)
   ! ngrvdw can be 0 or natom

  real(dp), allocatable :: grxc(:,:)
   ! grxc(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the XC energy

  real(dp) :: pel(3)
   ! ucvol times the electronic polarization in reduced coordinates

  real(dp) :: strten(6)
   ! Stress tensor in cartesian coordinates (Hartree/Bohr^3)
   ! 6 unique components of this symmetric 3x3 tensor:
   ! Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).

  real(dp), allocatable :: synlgr(:,:)
   ! synlgr(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the non-local energy
   ! The "sy" prefix refer to the fact that this gradient has been
   ! symmetrized.

 end type results_gs_type

!public procedures.
 public :: init_results_gs
 public :: init_results_gs_array
 public :: destroy_results_gs
 public :: destroy_results_gs_array
 public :: copy_results_gs
 public :: results_gs_ncwrite
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_results_gs/init_results_gs
!! NAME
!!  init_results_gs
!!
!! FUNCTION
!!  Init all (or part of) scalars and allocatables in a results_gs datastructure
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  nsppol=number of spin channels for this dataset
!!  only_part= --optional, default=false--
!!            if this flag is activated only the following parts of results_gs
!!            are initalized: all scalars, fcart,fred,strten
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs=<type(results_gs_type)>=results_gs datastructure
!!
!! PARENTS
!!      m_results_img
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine init_results_gs(natom,nsppol,results_gs,only_part)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_results_gs'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsppol
 logical,optional,intent(in) :: only_part
!arrays
 type(results_gs_type),intent(inout) :: results_gs
!Local variables-------------------------------
!scalars
 logical :: full_init
!arrays

!************************************************************************

 !@results_gs_type

 full_init=.true.;if (present(only_part)) full_init=(.not.only_part)

 results_gs%natom  =natom
 results_gs%ngrvdw =0
 results_gs%nsppol =nsppol
 results_gs%deltae =zero
 results_gs%diffor =zero
 results_gs%entropy=zero
 results_gs%etotal =zero
 results_gs%fermie =zero
 results_gs%residm =zero
 results_gs%res2   =zero
 results_gs%vxcavg =zero

 call energies_init(results_gs%energies)

 results_gs%strten=zero
 ABI_ALLOCATE(results_gs%fcart,(3,natom))
 results_gs%fcart=zero
 ABI_ALLOCATE(results_gs%fred,(3,natom))
 results_gs%fred =zero
 ABI_ALLOCATE(results_gs%gaps,(3,nsppol))
 results_gs%gaps =zero
 if (full_init) then
   results_gs%pel=zero
   ABI_ALLOCATE(results_gs%gresid,(3,natom))
   results_gs%gresid=zero
   ABI_ALLOCATE(results_gs%grewtn,(3,natom))
   results_gs%grewtn=zero
   ABI_ALLOCATE(results_gs%grxc,(3,natom))
   results_gs%grxc  =zero
   ABI_ALLOCATE(results_gs%synlgr,(3,natom))
   results_gs%synlgr=zero
   ABI_ALLOCATE(results_gs%grvdw,(3,results_gs%ngrvdw))
 end if

end subroutine init_results_gs
!!***

!----------------------------------------------------------------------

!!****f* m_results_gs/init_results_gs_array
!! NAME
!!  init_results_gs_array
!!
!! FUNCTION
!!  Init all (or part of) scalars and allocatables in a 2D-array of results_gs datastructures
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  nsppol=number of spin channels for this dataset
!!  only_part= --optional, default=false--
!!            if this flag is activated only the following parts of results_gs
!!            are initalized: all scalars, fcart,fred,strten
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs(:)=<type(results_gs_type)>=results_gs datastructure 2Darray
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine init_results_gs_array(natom,nsppol,results_gs,only_part)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_results_gs_array'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsppol
 logical,optional,intent(in) :: only_part
!arrays
 type(results_gs_type),intent(inout) :: results_gs(:,:)
!Local variables-------------------------------
!scalars
 integer :: ii,jj,results_gs_size1,results_gs_size2
 logical :: full_init
!arrays

!************************************************************************

 !@results_gs_type

 results_gs_size1=size(results_gs,1)
 results_gs_size2=size(results_gs,2)
 full_init=.true.;if (present(only_part)) full_init=(.not.only_part)

 if (results_gs_size1>0.and.results_gs_size2>0) then

   do ii=1,results_gs_size2
     do jj=1,results_gs_size1

       results_gs(jj,ii)%natom  =natom
       results_gs(jj,ii)%ngrvdw =0
       results_gs(jj,ii)%nsppol =nsppol
       results_gs(jj,ii)%deltae =zero
       results_gs(jj,ii)%diffor =zero
       results_gs(jj,ii)%entropy=zero
       results_gs(jj,ii)%etotal =zero
       results_gs(jj,ii)%fermie =zero
       results_gs(jj,ii)%residm =zero
       results_gs(jj,ii)%res2   =zero
       results_gs(jj,ii)%vxcavg =zero

       call energies_init(results_gs(jj,ii)%energies)

       results_gs(jj,ii)%strten=zero
       ABI_ALLOCATE(results_gs(jj,ii)%fcart,(3,natom))
       results_gs(jj,ii)%fcart=zero
       ABI_ALLOCATE(results_gs(jj,ii)%fred,(3,natom))
       results_gs(jj,ii)%fred =zero
       ABI_ALLOCATE(results_gs(jj,ii)%gaps,(3,nsppol))
       results_gs(jj,ii)%gaps =zero
       if (full_init) then
         results_gs(jj,ii)%pel=zero
         ABI_ALLOCATE(results_gs(jj,ii)%gresid,(3,natom))
         results_gs(jj,ii)%gresid=zero
         ABI_ALLOCATE(results_gs(jj,ii)%grewtn,(3,natom))
         results_gs(jj,ii)%grewtn=zero
         ABI_ALLOCATE(results_gs(jj,ii)%grxc,(3,natom))
         results_gs(jj,ii)%grxc  =zero
         ABI_ALLOCATE(results_gs(jj,ii)%synlgr,(3,natom))
         results_gs(jj,ii)%synlgr=zero
         ABI_ALLOCATE(results_gs(jj,ii)%grvdw,(3,results_gs(jj,ii)%ngrvdw))
       end if

     end do
   end do
 end if

end subroutine init_results_gs_array
!!***

!----------------------------------------------------------------------

!!****f* m_results_gs/destroy_results_gs
!! NAME
!!  destroy_results_gs
!!
!! FUNCTION
!!  Clean and destroy a results_gs datastructure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs(:)=<type(results_gs_type)>=results_gs datastructure
!!
!! PARENTS
!!      m_results_img
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine destroy_results_gs(results_gs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_results_gs'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(results_gs_type),intent(inout) :: results_gs
!Local variables-------------------------------

!************************************************************************

 !@results_gs_type

 results_gs%natom =0
 results_gs%ngrvdw=0
 results_gs%nsppol=0
 if (allocated(results_gs%fcart))   then
   ABI_DEALLOCATE(results_gs%fcart)
 end if
 if (allocated(results_gs%fred))    then
   ABI_DEALLOCATE(results_gs%fred)
 end if
 if (allocated(results_gs%gaps))    then
   ABI_DEALLOCATE(results_gs%gaps)
 end if
 if (allocated(results_gs%gresid))  then
   ABI_DEALLOCATE(results_gs%gresid)
 end if
 if (allocated(results_gs%grewtn))  then
   ABI_DEALLOCATE(results_gs%grewtn)
 end if
 if (allocated(results_gs%grvdw))  then
   ABI_DEALLOCATE(results_gs%grvdw)
 end if
 if (allocated(results_gs%grxc))    then
   ABI_DEALLOCATE(results_gs%grxc)
 end if
 if (allocated(results_gs%synlgr))  then
   ABI_DEALLOCATE(results_gs%synlgr)
 end if

end subroutine destroy_results_gs
!!***

!----------------------------------------------------------------------

!!****f* m_results_gs/destroy_results_gs_array
!! NAME
!!  destroy_results_gs_array
!!
!! FUNCTION
!!  Clean and destroy a 2D-array of results_gs datastructures
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  results_gs(:)=<type(results_gs_type)>=results_gs datastructure 2D-array
!!
!! PARENTS
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine destroy_results_gs_array(results_gs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_results_gs_array'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(results_gs_type),intent(inout) :: results_gs(:,:)
!Local variables-------------------------------
!scalars
 integer :: ii,jj,results_gs_size1,results_gs_size2

!************************************************************************

 !@results_gs_type

 results_gs_size1=size(results_gs,1)
 results_gs_size2=size(results_gs,2)
 if (results_gs_size1>0.and.results_gs_size2>0) then

   do ii=1,results_gs_size2
     do jj=1,results_gs_size1
       results_gs(jj,ii)%natom =0
       results_gs(jj,ii)%ngrvdw=0
       results_gs(jj,ii)%nsppol=0
       if (allocated(results_gs(jj,ii)%fcart))   then
         ABI_DEALLOCATE(results_gs(jj,ii)%fcart)
       end if
       if (allocated(results_gs(jj,ii)%fred))    then
         ABI_DEALLOCATE(results_gs(jj,ii)%fred)
       end if
       if (allocated(results_gs(jj,ii)%gaps))   then
         ABI_DEALLOCATE(results_gs(jj,ii)%gaps)
       end if
       if (allocated(results_gs(jj,ii)%gresid))  then
         ABI_DEALLOCATE(results_gs(jj,ii)%gresid)
       end if
       if (allocated(results_gs(jj,ii)%grewtn))  then
         ABI_DEALLOCATE(results_gs(jj,ii)%grewtn)
       end if
       if (allocated(results_gs(jj,ii)%grvdw))  then
         ABI_DEALLOCATE(results_gs(jj,ii)%grvdw)
       end if
       if (allocated(results_gs(jj,ii)%grxc))    then
         ABI_DEALLOCATE(results_gs(jj,ii)%grxc)
       end if
       if (allocated(results_gs(jj,ii)%synlgr))  then
         ABI_DEALLOCATE(results_gs(jj,ii)%synlgr)
       end if
     end do
   end do

 end if

end subroutine destroy_results_gs_array
!!***
!----------------------------------------------------------------------

!!****f* m_results_gs/copy_results_gs
!! NAME
!!  copy_results_gs
!!
!! FUNCTION
!!  Copy a results_gs datastructure into another
!!
!! INPUTS
!!  results_gs_in=<type(results_gs_type)>=input results_gs datastructure
!!
!! OUTPUT
!!  results_gs_out=<type(results_gs_type)>=output results_gs datastructure
!!
!! PARENTS
!!      m_results_img
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine copy_results_gs(results_gs_in,results_gs_out)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_results_gs'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(results_gs_type),intent(in) :: results_gs_in
 type(results_gs_type),intent(inout) :: results_gs_out !vz_i
!Local variables-------------------------------
!scalars
 integer :: natom_in,natom_out,ngrvdw_in,nsppol_in,nsppol_out

!************************************************************************

 !@results_gs_type

 natom_in =results_gs_in%natom
 natom_out=results_gs_out%natom
 ngrvdw_in =results_gs_in%ngrvdw
 nsppol_in =results_gs_in%nsppol
 nsppol_out =results_gs_out%nsppol

 if (natom_in>natom_out) then
   if (allocated(results_gs_out%fcart))   then
     ABI_DEALLOCATE(results_gs_out%fcart)
   end if
   if (allocated(results_gs_out%fred))    then
     ABI_DEALLOCATE(results_gs_out%fred)
   end if
   if (allocated(results_gs_out%gresid))  then
     ABI_DEALLOCATE(results_gs_out%gresid)
   end if
   if (allocated(results_gs_out%grewtn))  then
     ABI_DEALLOCATE(results_gs_out%grewtn)
   end if
   if (allocated(results_gs_out%grvdw))  then
     ABI_DEALLOCATE(results_gs_out%grvdw)
   end if
   if (allocated(results_gs_out%grxc))    then
     ABI_DEALLOCATE(results_gs_out%grxc)
   end if
   if (allocated(results_gs_out%synlgr))  then
     ABI_DEALLOCATE(results_gs_out%synlgr)
   end if

   if (allocated(results_gs_in%fcart))   then
     ABI_ALLOCATE(results_gs_out%fcart,(3,natom_in))
   end if
   if (allocated(results_gs_in%fred))    then
     ABI_ALLOCATE(results_gs_out%fred,(3,natom_in))
   end if
   if (allocated(results_gs_in%gresid))  then
     ABI_ALLOCATE(results_gs_out%gresid,(3,natom_in))
   end if
   if (allocated(results_gs_in%grewtn))  then
     ABI_ALLOCATE(results_gs_out%grewtn,(3,natom_in))
   end if
   if (allocated(results_gs_in%grvdw))  then
     ABI_ALLOCATE(results_gs_out%grvdw,(3,ngrvdw_in))
   end if
   if (allocated(results_gs_in%grxc))    then
     ABI_ALLOCATE(results_gs_out%grxc,(3,natom_in))
   end if
   if (allocated(results_gs_in%synlgr))  then
     ABI_ALLOCATE(results_gs_out%synlgr,(3,natom_in))
   end if
 end if

 if (nsppol_in>nsppol_out) then
   if (allocated(results_gs_out%gaps))   then
     ABI_DEALLOCATE(results_gs_out%gaps)  
   end if
   if (allocated(results_gs_in%gaps))    then
     ABI_ALLOCATE(results_gs_out%gaps,(3,nsppol_in))
   end if 
 endif


 results_gs_out%natom  =results_gs_in%natom
 results_gs_out%ngrvdw =results_gs_in%ngrvdw
 results_gs_out%nsppol =results_gs_in%nsppol
 results_gs_out%deltae =results_gs_in%deltae
 results_gs_out%diffor =results_gs_in%diffor
 results_gs_out%entropy=results_gs_in%entropy
 results_gs_out%etotal =results_gs_in%etotal
 results_gs_out%fermie =results_gs_in%fermie
 results_gs_out%residm =results_gs_in%residm
 results_gs_out%res2   =results_gs_in%res2
 results_gs_out%vxcavg =results_gs_in%vxcavg

 call energies_copy(results_gs_in%energies,results_gs_out%energies)

 results_gs_out%pel(:)=results_gs_in%pel(:)
 results_gs_out%strten(:)=results_gs_in%strten(:)

 if (allocated(results_gs_in%fcart))  results_gs_out%fcart(:,1:natom_in) =results_gs_in%fcart(:,1:natom_in)
 if (allocated(results_gs_in%fred))   results_gs_out%fred(:,1:natom_in)  =results_gs_in%fred(:,1:natom_in)
 if (allocated(results_gs_in%gaps))   results_gs_out%gaps(:,1:nsppol_in) =results_gs_in%gaps(:,1:nsppol_in)
 if (allocated(results_gs_in%gresid)) results_gs_out%gresid(:,1:natom_in)=results_gs_in%gresid(:,1:natom_in)
 if (allocated(results_gs_in%grewtn)) results_gs_out%grewtn(:,1:natom_in)=results_gs_in%grewtn(:,1:natom_in)
 if (allocated(results_gs_in%grxc))   results_gs_out%grxc(:,1:natom_in)  =results_gs_in%grxc(:,1:natom_in)
 if (allocated(results_gs_in%synlgr)) results_gs_out%synlgr(:,1:natom_in)=results_gs_in%synlgr(:,1:natom_in)
 if (allocated(results_gs_in%grvdw).and.ngrvdw_in>0) then
   results_gs_out%grvdw(:,1:ngrvdw_in)=results_gs_in%grvdw(:,1:ngrvdw_in)
 end if

end subroutine copy_results_gs
!!***

!----------------------------------------------------------------------


!!****f* m_results_gs/results_gs_ncwrite
!! NAME
!! results_gs_ncwrite
!!
!! FUNCTION
!!
!! INPUTS
!!  ncid=NC file handle
!!  ecut=Cutoff enery in Ha.
!!
!! OUTPUT
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine results_gs_ncwrite(results_gs,ecut,ncid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'results_gs_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 real(dp),intent(in) :: ecut
 type(results_gs_type),intent(in) :: results_gs
!arrays

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(ETSF_io_low_error) :: Error_data

! *************************************************************************
 ! ==============================================
 ! === Write the dimensions specified by ETSF ===
 ! ==============================================
 !FIXME: do not handle k_dependent = 1
 call etsf_io_low_set_define_mode(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! TODO check whether this name is already used
 call etsf_io_low_write_dim(ncid,"dimension_six",6,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,"number_of_cartesian_dimensions",3,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_atoms',results_gs%natom,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_spins',results_gs%nsppol,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

! Define variables.
! scalars passed in input (not belonging to results_gs)
 call etsf_io_low_def_var(ncid,'ecut',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

! scalars defined in results_gs
 call etsf_io_low_def_var(ncid,'deltae',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'diffor',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'entropy',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'etotal',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'fermie',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'residm',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'res2',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! arrays
 !
 ! Note: unlike fred, this array has been corrected by enforcing
 ! the translational symmetry, namely that the sum of force on all atoms is zero.
 call etsf_io_low_def_var(ncid,'cartesian_forces',etsf_io_low_double,&
&  (/pad("number_of_cartesian_dimensions"),pad("number_of_atoms")/),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'cartesian_stress_tensor',etsf_io_low_double,&
&  (/'dimension_six'/),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

! Write data.
 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_atoms',results_gs%natom,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_spins',results_gs%nsppol,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)


! Write variables passed in input
 call etsf_io_low_write_var(ncid,'ecut',ecut,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

! Write results_gs variables
 call etsf_io_low_write_var(ncid,'deltae',results_gs%deltae,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'diffor',results_gs%diffor,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'entropy',results_gs%entropy,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'etotal',results_gs%etotal,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'fermie',results_gs%fermie,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'residm',results_gs%residm,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'res2',results_gs%res2,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'cartesian_forces',results_gs%fcart,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'cartesian_stress_tensor',results_gs%strten,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

#else 
 MSG_ERROR("ETSF-IO support is not activated.")
#endif

end subroutine results_gs_ncwrite
!!***

!----------------------------------------------------------------------

END MODULE m_results_gs
!!***
