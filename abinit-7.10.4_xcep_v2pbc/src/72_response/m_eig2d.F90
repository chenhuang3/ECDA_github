!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_eig2d
!! NAME
!!  m_eig2d
!!
!! FUNCTION
!!  This module contains utilities to analyze and retrieve information
!!  from the second order derivative of the eigen-energies wrt
!!  displacements.
!!
!! COPYRIGHT
!! Copyright (C) 2014-2014 ABINIT group (SP)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!  loper3
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_eig2d

 use defs_basis
 use m_errors
 use m_profiling_abi
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 implicit none

 private

 public :: eigr2d_init              ! Main creation method of EIG2D.nc files.
 public :: eigr2d_ncwrite           ! Dump the object into NETCDF file.
 public :: eigr2d_free              ! Destruction method.
 public :: fan_init                 ! Main creation method of Fan.nc files.
 public :: fan_ncwrite              ! Dump the object into NETCDF file.
 public :: fan_free                 ! Destruction method.

!!***

!!****t* m_eig2d/eigr2d_t
!! NAME
!! eig2d_t
!!
!! FUNCTION
!! It contains informations about the second-order derivative of the
!! eigenenergies wrt atomic displacement
!!
!! SOURCE

 type,public :: eigr2d_t

! WARNING : if you modify this datatype, please check whether there might be
! creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your
! modification.

  integer :: mband                 ! Max number of bands i.e MAXVAL(nband) (to dimension arrays)
  integer :: nsppol                ! number of spin-polarization
  integer :: nkpt                  ! number of k points
  integer :: natom                 ! number of atoms

  real(dp),allocatable :: eigr2d(:,:,:,:,:,:,:)
  ! eigr2d(2,mband*nsppol,nkpt,3,natom,3,natom)
  ! Second-order derivative of eigenergies (real,im) at each
  ! band,k-point,dir1,dir2,natom1,natom2 .


 end type eigr2d_t
!!***

!!****t* m_eig2d/fan_t
!! NAME
!! fan_t
!!
!! FUNCTION
!! It contains informations about the second-order derivative of the 
!! eigenenergies wrt atomic displacement
!!
!! SOURCE

 type,public :: fan_t

! WARNING : if you modify this datatype, please check whether there might be
! creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your
! modification.

  integer :: mband                 ! Max number of bands i.e MAXVAL(nband) (to dimension arrays)
  integer :: nsppol                ! number of spin-polarization
  integer :: nkpt                  ! number of k points
  integer :: natom                 ! number of atoms

  real(dp),allocatable :: fan2d(:,:,:,:,:,:,:)
  ! fan2d(2*mband*nsppol,nkpt,3,natom,3,natom,mband*nsppol)
  ! Second-order derivative of the eigenergies (real,im) at each
  ! iband(real,im),k-point,dir1,dir2,natom1,natom2,jband

 end type fan_t
!!***




CONTAINS
!=====================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/eigr2d_init
!! NAME
!! eigr2d_init
!!
!! FUNCTION
!! This subroutine initializes the eigr2d_t structured datatype
!!
!! INPUTS
!! mbands=maximum number of bands
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! natom=number of atoms
!! eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom)=second-order derivative of the
!!     eigenenergies wrt phononic displacements
!!
!! OUTPUT
!! eigr2d<eigr2d_t>=the eigr2d_t datatype
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      eig2tot,loper3
!!
!! CHILDREN
!!
!! SOURCE

subroutine eigr2d_init(eig2nkq,eigr2d,mband,nsppol,nkpt,natom)

 use defs_basis
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eigr2d_init'
!End of the abilint section

 integer,intent(in) ::mband,nsppol,nkpt,natom
 type(eigr2d_t),intent(out) :: eigr2d
!arrays
 real(dp), intent(in) :: eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom)

! *************************************************************************

 eigr2d%mband = mband
 eigr2d%nsppol = nsppol
 eigr2d%nkpt = nkpt
 eigr2d%natom = natom
 
 ABI_MALLOC(eigr2d%eigr2d  ,(2,mband*nsppol,nkpt,3,natom,3,natom))
 eigr2d%eigr2d=eig2nkq

end subroutine eigr2d_init
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/eigr2d_ncwrite
!! NAME
!! eigr2d_ncwrite
!!
!! FUNCTION
!!  Writes the content of a eigr2d_t object to a NETCDF file 
!!  according to the ETSF-IO specifications.
!!    
!! INPUTS
!!  ncid =NC file handle
!! 
!! OUTPUT
!!  
!! PARENTS
!!      eig2tot,loper3
!!
!! CHILDREN
!!
!! SOURCE

subroutine eigr2d_ncwrite(eigr2d,iqpt,wtq,ncid)
 
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eigr2d_ncwrite'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) ::ncid
 real(dp), intent(in) :: iqpt(3),wtq
 type(eigr2d_t),target,intent(in) :: eigr2d

!Local variables-------------------------------
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 integer :: cplex,cart_dir,one_dim
 type(ETSF_io_low_error) :: Error_data
#endif

! *************************************************************************

#ifdef HAVE_TRIO_ETSF_IO
 ! ==============================================
 ! === Write the dimensions specified by ETSF ===
 ! ==============================================
 one_dim=1
 cplex=2
 cart_dir=3

 call etsf_io_low_set_define_mode(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'max_number_of_states',eigr2d%mband,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_spins',eigr2d%nsppol,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_kpoints',eigr2d%nkpt,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_atoms',eigr2d%natom,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_cartesian_directions',cart_dir,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'current_one_dim',one_dim,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'cplex',cplex,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'product_mband_nsppol',eigr2d%mband*eigr2d%nsppol,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'current_q_point',etsf_io_low_double,(/pad('number_of_cartesian_directions')/),&
&                        lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'current_q_point_weight',etsf_io_low_double,(/pad('current_one_dim')/),&
&                        lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'second_derivative_eigenenergies',etsf_io_low_double,(/pad('cplex'),&
&  pad('product_mband_nsppol'),pad('number_of_kpoints'),pad('number_of_cartesian_directions'),&
&  pad('number_of_atoms'),pad('number_of_cartesian_directions'),pad('number_of_atoms')/), lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'current_q_point',iqpt,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'current_q_point_weight',wtq,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)
 
 call etsf_io_low_write_var(ncid,'second_derivative_eigenenergies',eigr2d%eigr2d,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

#else 
 MSG_ERROR("ETSF-IO support is not activated. ")
#endif

end subroutine eigr2d_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/eigr2d_free
!! NAME
!! eigr2d_free
!!
!! FUNCTION
!! Deallocates the components of the eigr2d_t structured datatype
!!
!! INPUTS
!!  eigr2d<eigr2d_t>=The data type to be deallocated.
!!
!! OUTPUT
!!  Deallocate the dynamic arrays in the ebands_t type.
!!  (only deallocate)
!!
!! PARENTS
!!      eig2tot,loper3
!!
!! CHILDREN
!!
!! SOURCE

subroutine eigr2d_free(eigr2d)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eigr2d_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(eigr2d_t),intent(inout) :: eigr2d
! *************************************************************************
DBG_ENTER("COLL")

!Deallocate all components of bstruct

 if (allocated(eigr2d%eigr2d)) then
   ABI_FREE(eigr2d%eigr2d)
 end if

 DBG_EXIT("COLL")

end subroutine eigr2d_free
!!***

!!****f* m_eig2d/fan_init
!! NAME
!! fan_init
!!
!! FUNCTION
!! This subroutine initializes the fan_t structured datatype
!!
!! INPUTS
!! mbands=maximum number of bands
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! natom=number of atoms
!! fan2d(2*mband*nsppol,nkpt,3,natom,3,natom,mband*nsppol)=second-order derivative of the
!!     eigenenergies wrt phononic displacements
!!
!! OUTPUT
!! fan2d<fan_t>=the fan_t datatype
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine fan_init(fan,fan2d,mband,nsppol,nkpt,natom)

 use defs_basis
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fan_init'
!End of the abilint section

 integer,intent(in) ::mband,nsppol,nkpt,natom
 type(fan_t),intent(out) :: fan2d
!arrays
 real(dp), intent(in) :: fan(2*mband*nsppol,nkpt,3,natom,3,natom,mband*nsppol)
! *************************************************************************

 fan2d%mband = mband
 fan2d%nsppol = nsppol
 fan2d%nkpt = nkpt
 fan2d%natom = natom

 ABI_MALLOC(fan2d%fan2d,(2*mband*nsppol,nkpt,3,natom,3,natom,mband*nsppol))
 fan2d%fan2d=fan

end subroutine fan_init
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/fan_ncwrite
!! NAME
!! fan_ncwrite
!!
!! FUNCTION
!!  Writes the content of a fan_t object to a NETCDF file 
!!  according to the ETSF-IO specifications.
!!    
!! INPUTS
!!  ncid =NC file handle
!! 
!! OUTPUT
!!  
!! PARENTS
!!      eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine fan_ncwrite(fan2d,iqpt,wtq,ncid)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fan_ncwrite'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) ::ncid
 real(dp), intent(in) :: iqpt(3),wtq
 type(fan_t),target,intent(in) :: fan2d

!Local variables-------------------------------
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 integer :: cplex,cart_dir,one_dim
 type(ETSF_io_low_error) :: Error_data
#endif

! *************************************************************************

#ifdef HAVE_TRIO_ETSF_IO
 ! ==============================================
 ! === Write the dimensions specified by ETSF ===
 ! ==============================================
 one_dim=1
 cplex=2
 cart_dir=3

 call etsf_io_low_set_define_mode(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'max_number_of_states',fan2d%mband,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_spins',fan2d%nsppol,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_kpoints',fan2d%nkpt,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_atoms',fan2d%natom,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'3_number_of_atoms',3*(fan2d%natom),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_cartesian_directions',cart_dir,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'current_one_dim',one_dim,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'cplex',cplex,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'product_mband_nsppol',fan2d%mband*fan2d%nsppol,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'product_mband_nsppol2',fan2d%mband*fan2d%nsppol*2,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'current_q_point',etsf_io_low_double,(/pad('number_of_cartesian_directions')/),&
&                        lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'current_q_point_weight',etsf_io_low_double,(/pad('current_one_dim')/),&
&                        lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'second_derivative_eigenenergies_actif',etsf_io_low_double,(/&
&  pad('product_mband_nsppol2'),pad('number_of_kpoints'),pad('number_of_cartesian_directions'),&
&  pad('number_of_atoms'),pad('number_of_cartesian_directions'),pad('number_of_atoms'),&
&  pad('product_mband_nsppol')/), lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'current_q_point',iqpt,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'current_q_point_weight',wtq,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'second_derivative_eigenenergies_actif',fan2d%fan2d,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

#else 
 MSG_ERROR("ETSF-IO support is not activated. ")
#endif

end subroutine fan_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_eig2d/fan_free
!! NAME
!! fan_free
!!
!! FUNCTION
!! Deallocates the components of the fan_t structured datatype
!!
!! INPUTS
!!  fan2d<fan_t>=The data type to be deallocated.
!!
!! OUTPUT
!!  Deallocate the dynamic arrays in the fan_t type.
!!  (only deallocate)
!!
!! PARENTS
!!      eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine fan_free(fan2d)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fan_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(fan_t),intent(inout) :: fan2d
! *************************************************************************
DBG_ENTER("COLL")

!Deallocate all components of bstruct

 if (allocated(fan2d%fan2d)) then
   ABI_FREE(fan2d%fan2d)
 end if

 DBG_EXIT("COLL")

end subroutine fan_free



END MODULE m_eig2d
!!***

