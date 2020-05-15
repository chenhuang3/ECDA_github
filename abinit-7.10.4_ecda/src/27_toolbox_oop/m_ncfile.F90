!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ncfile
!! NAME
!! m_ncfile
!!
!! FUNCTION
!!  This module provides tools for the IO of NETCDF files.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2014 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_ncfile

 use defs_basis
 use m_profiling_abi
 use m_errors
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif

 use m_fstrings,  only : itoa

 implicit none

 private 
!!***

 ! Constants mapping the name used in ABINIT to ETSF-IO name (just to improve the readability of the code)
 !character(len=*),parameter :: natom_vname = 'number_of_atoms'
 !character(len=*),parameter :: nsppol_vname = 'number_of_spins'
 !character(len=*),parameter :: nspinor_vname = 'number_of_spinor_components',
 !character(len=*),parameter :: nspden_vname = 'number_of_spin_density_components'

!!****t* m_ncfile/ncfile_t
!! NAME
!! ncfile_t
!! 
!! FUNCTION
!!  File descriptor used for Netcdf files.
!! 
!! SOURCE

 type,public :: ncfile_t
   integer :: ncid                 ! File handle
   integer :: cmode                ! cmode passed to ncfile_create
   integer :: comm                 ! MPI communicator.
   character(len=fnlen) :: fname   ! File name
 end type ncfile_t
!!***

#if defined(HAVE_TRIO_ETSF_IO) || defined(HAVE_TRIO_NETCDF)
 public :: ncfile_create       ! Create a new netcdf file.
 !public :: ncfile_open        ! Close the netcdf file.
 public :: ncfile_close        ! Define the basic dimensions used in ETSF-IO files. 
 public :: ab_define_var       ! Helper function used to declara a netcdf variable.

#ifdef HAVE_TRIO_ETSF_IO
 ! Helper functions
 public :: ncid_define_basedims
#endif

CONTAINS

!!****f* m_ncfile/ncfile_create
!! NAME
!!  ncfile_create
!!
!! FUNCTION
!!  Create a new netcdf file.
!!
!! INPUTS
!!  path=Filepath.
!!  cmode=Creation mode
!!
!! OUTPUT
!!  ncf<ncfile_t>=netcdf file descriptor.
!!  ncerr=Status error.
!!
!! PARENTS
!!
!! SOURCE

function ncfile_create(ncf,path,cmode,comm) result(ncerr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ncfile_create'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ncfile_t),intent(out) :: ncf
 character(len=*),intent(in) :: path
 integer,intent(in) :: cmode
 integer,optional,intent(in) :: comm
 integer :: ncerr

! *********************************************************************

 ncf%fname = path
 ncf%cmode = cmode

 ncerr = nf90_create(path, cmode, ncf%ncid)

 RETURN
 ABI_UNUSED(comm)

end function ncfile_create
!!***

!----------------------------------------------------------------------

!!****f* m_ncfile/ncfile_close
!! NAME
!!  ncfile_close
!!
!! FUNCTION
!!  Close the netcdf file descriptor.
!!
!! INPUTS
!!  ncf<ncfile_t>=netcdf file descriptor.
!!
!! OUTPUT
!!  ncerr=Status error.
!!
!! PARENTS
!!
!! SOURCE

function ncfile_close(ncf) result(ncerr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ncfile_close'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ncfile_t),intent(inout) :: ncf
 integer :: ncerr

! *********************************************************************

 ncerr = nf90_close(ncf%ncid)

end function ncfile_close
!!***

!----------------------------------------------------------------------

#ifdef HAVE_TRIO_ETSF_IO
!!****f* m_ncfile/ncid_define_basedims
!! NAME
!!  ncid_define_basedims
!!
!! FUNCTION
!!  Define the basic dimensions used in ETSF-IO files. 
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dfptdb,m_header,m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine ncid_define_basedims(ncid) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ncid_define_basedims'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ncid

!Local variables-------------------------------
!scalars
! integer :: ii
 logical :: lstat
 type(ETSF_io_low_error) :: Error_data
 !character(len=20) :: numbers(10) 

! *********************************************************************

 call etsf_io_low_set_define_mode(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Basic ETSF-IO dimensions that should be always present in the file.
 call etsf_io_low_write_dim(ncid,'complex',2,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,"symbol_length",2,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'character_string_length',80,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,"number_of_cartesian_directions",3,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,"number_of_reduced_dimensions",3,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,"number_of_vectors",3,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! Useful integers. XLF@IBM6 does not like it!
 ! numbers = &
 !&  [pad("one"), pad("two"), pad("three"), pad("four"), pad("five"), pad("six"), pad("seven"), pad("eight"), pad("nine"), pad("ten")]
 !do ii=1,size(numbers)
 !  call etsf_io_low_write_dim(ncid,itoa(index),index,lstat,Error_data=Error_data)  
 !  call etsf_io_low_write_dim(ncid,numbers(ii),ii,lstat,Error_data=Error_data)  
 !  ETSF_CHECK_ERROR(lstat,Error_data)
 !end do

end subroutine ncid_define_basedims
!!***
#endif

!!****f* ABINIT/ab_define_var
!!
!! NAME
!! ab_define_var
!!
!! FUNCTION
!! Write the definition of a variable, including units and mnemonics
!!
!! INPUTS
!! ncid = Identifier of the netcdf dataset
!! var_dim_id = Identifier of the Dimensions
!! var_id     = Identifier of the variable
!! var_mnemo  = String of mnemonics
!! var_name   = String with the name of the variable
!! var_type   = NetCDF type of variable (NF90_DOUBLE, etc)
!! var_units  = String of units
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      m_abihist,m_bse_io,write_eig
!!
!! CHILDREN
!!
!! SOURCE

subroutine ab_define_var(ncid, var_dim_id, var_id, var_type, var_name, var_mnemo, var_units)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ab_define_var'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ncid 
 integer, intent(out) :: var_id
 character(len=*), intent(in) :: var_mnemo,var_units
 character(len=*), intent(in) :: var_name
 integer,intent(in) :: var_type
!arrays
 integer,intent(in) :: var_dim_id(:)

!Local variables-------------------------------
!scalars
 integer :: ncerr

! *************************************************************************

#if defined HAVE_TRIO_NETCDF
 ncerr = nf90_def_var(ncid, trim(var_name), var_type, var_dim_id, var_id)
 NCF_CHECK(ncerr," define variable "//trim(var_name))
 
 ncerr = nf90_put_att(ncid, var_id,  "units",trim(var_units))
 NCF_CHECK(ncerr," define attribute for "//trim(var_name))

 ncerr = nf90_put_att(ncid, var_id,  "mnemonics", trim(var_mnemo))
 NCF_CHECK(ncerr," define attribute for "//trim(var_name))
#endif

end subroutine ab_define_var
!!***

#endif

END MODULE m_ncfile
!!***

