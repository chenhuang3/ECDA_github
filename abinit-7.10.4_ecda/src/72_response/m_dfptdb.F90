!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfptdb
!! NAME
!!  m_dfptdb
!!
!! FUNCTION
!!  This module provides an object, dfptdb_t, that allows one to read/write/merge 
!!  the first order change of the density and the first order change of the potential
!!  obtained in the DFPT code. Results are in netcdf format.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2013 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! NOTES
!!  1) DFPT results are written in netcdf format (etsf-io names are used whenever possible)
!!
!!  2) This version does not support MPI-IO. When writing, we assume the data is fully available
!!     on the master node.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_dfptdb

 use defs_basis
 use m_profiling_abi
 use m_errors
!!#ifdef HAVE_TRIO_ETSF_IO
!! use netcdf
!! use etsf_io
!! use m_ncfile,     only : ncid_define_basedims
!!#endif

 use m_fstrings,   only : strcat, sjoin, itoa
 use m_io_tools,   only : delete_file, file_exists
 use m_copy,       only : alloc_copy
 use m_pawrhoij,   only : pawrhoij_type, pawrhoij_destroy

 implicit none

 private
!!***

 ! FIXME
 real(dp),public,parameter :: DDB_QTOL=2.0d-8
 ! Tolerance for the identification of two wavevectors

!----------------------------------------------------------------------

!!****t* m_dfptdb/dfptdb_t
!! NAME
!!  dfptdb_t
!!
!! FUNCTION
!!  Database of DFPT retults. The database contains npert perturbations
!!  and the corresponding first order change of the density and of the local potential
!!  in real space on the FFT mesh. Note that one can have different FFT meshes 
!!  for the different perturbations
!!
!! SOURCE

 type,public :: dfptdb_t

  integer :: ncid
   ! Netcdf file handler

  integer :: npert=0
  ! Number of perturbations present in file.

  ! natom, nspden, nspinor, and usepaw are global
  ! in the sense that it's not possible to add
  ! new entries in the database with different values.

  integer :: natom
   ! Number of atoms

  integer :: nspden
   ! Number of spin density components

  integer :: nspinor
   ! Number of spinor components.

  integer :: usepaw
   ! 1 if PAW calculation, 0 otherwise

  character(len=fnlen) :: filename = ABI_NOFILE
   ! File name 

  integer,allocatable :: idir_pert(:), ipert_pert(:)
   ! idir_pert(npert), ipert_pert(npert)
   ! Value of idir and ipert associated to the perturbation.

  real(dp),allocatable :: qpt_pert(:,:)
   ! qpt_pert(3,npert)
   ! q-points of the perturbation in reduced coordinates

  integer,allocatable :: cplex_pert(:)
   ! 2 if density and potential are complex, 1 if real
   ! cplex_pert(npert)

  integer,allocatable :: ngfft3_pert(:,:)
   ! ngfft3_pert(3,npert)
   ! FFT mesh for the different perturbations.

  ! Stuff needed for pawrhoij1
  !integer :: ntypat

  !integer,allocatable :: nlmn_type(:)
  !  nlmn_type(ntypat)= Number of (l,m,n) elements for the paw basis for each type of atom. Only used for reading.

  !integer,allocatable :: typat(:)
  ! typat(natom) =Type of each atom.

 end type dfptdb_t

#if 0
 public :: dfptdb_free                 ! Release memory
 public :: dfptdb_close                ! Close the file and release the memory allocated.
 public :: dfptdb_open_read            ! Open the file in read mode.
 public :: dfptdb_open_write           ! Open the file in write mode.
 public :: dfptdb_write                ! Write a new perturbation in the dabase.
 public :: dfptdb_read                 ! Read data.
 public :: dfptdb_merge                ! Merge multiple databases
 !!public :: dfptdb_prepare_fourier
!!***

 interface add_record
   module procedure add_record_int0d
   module procedure add_record_int1d
 end interface add_record

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/dfptdb_has_pert
!! NAME
!!  dfptdb_has_pert
!!
!! FUNCTION
!!  Returns the index of the given perturbation in the file. 0 if not found.
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

integer function dfptdb_find_pert(db,qpt,ipert,idir) result(ip)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptdb_find_pert'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ipert,idir
 type(dfptdb_t),intent(in) :: db
!arrays
 real(dp),intent(in) :: qpt(3)

!************************************************************************

 do ip=1,db%npert
   if (db%idir_pert(ip) == idir .and.   &
&      db%ipert_pert(ip) == ipert .and. &  
&      all(abs(db%qpt_pert(:,ip) - qpt) < DDB_QTOL)) return
 end do
 ip = 0

end function dfptdb_find_pert
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/dfptdb_free
!! NAME
!!  dfptdb_free
!!
!! FUNCTION
!!  Release memory
!!
!! PARENTS
!!      m_dfptdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfptdb_free(db)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptdb_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dfptdb_t),intent(inout) :: db

!************************************************************************

 if (allocated(db%idir_pert)) then
   ABI_FREE(db%idir_pert)
 end if

 if (allocated(db%ipert_pert)) then
   ABI_FREE(db%ipert_pert)
 end if

 if (allocated(db%qpt_pert)) then
   ABI_FREE(db%qpt_pert)
 end if

 if (allocated(db%cplex_pert)) then
   ABI_FREE(db%cplex_pert)
 end if

 if (allocated(db%ngfft3_pert)) then
   ABI_FREE(db%ngfft3_pert)
 end if

end subroutine dfptdb_free
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/dfptdb_close
!! NAME
!!  dfptdb_close
!!
!! FUNCTION
!!  Close the file and release the memory allocated
!!
!! INPUT
!!  [delete]=True if file should be deleted. Default: False
!!
!! PARENTS
!!      m_dfptdb,scfcv3
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfptdb_close(db,delete)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptdb_close'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dfptdb_t),intent(inout) :: db
 logical,optional,intent(in) :: delete

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 integer :: ierr
 logical :: my_delete,lstat
 type(etsf_io_low_error) :: Error

! *************************************************************************

 call etsf_io_low_close(db%ncid, lstat, Error_data=Error)
 ETSF_CHECK_ERROR(lstat, Error)

 if (present(delete)) then 
   if (delete) call delete_file(db%filename,ierr)
 end if

 ! Free memory
 call dfptdb_free(db)
#endif

end subroutine dfptdb_close
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/dfptdb_open_read
!! NAME
!!  dfptdb_open_read
!!
!! FUNCTION
!!  Open the file in read mode.
!!
!! INPUTS
!!  filename = Name of the file
!!
!! OUTPUT
!!  db<type(dfptdb_t)> = File handler initialized and set in read mode
!!
!! PARENTS
!!      m_dfptdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfptdb_open_read(db,filename)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptdb_open_read'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: filename
 type(dfptdb_t),intent(out) :: db

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 integer :: ip,fform
 character(len=500) :: msg
 logical :: lstat
 type(etsf_io_low_error) :: Error

!************************************************************************

 db%filename = filename

 call etsf_io_low_open_read(db%ncid,db%filename,lstat,Error_data=Error,with_etsf_header=.FALSE.)
 ETSF_CHECK_ERROR(lstat, Error)

 ! Read the list of perturbations stored in the file.
 call etsf_io_low_read_var(db%ncid,'number_of_perturbations',db%npert,lstat,error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 ! Read global variables.
 call etsf_io_low_read_dim(db%ncid,'number_of_atoms',db%natom,lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 !?
 call etsf_io_low_read_dim(db%ncid,'number_of_spinor_components',db%nspinor,lstat,error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_read_dim(db%ncid,'number_of_spin_density_components',db%nspden,lstat,error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_read_dim(db%ncid,'usepaw',db%usepaw,lstat,error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)
 ABI_CHECK(db%usepaw==0,"PAW not yet coded")
 
 ! Allocate and read list of perturbations.
 ABI_MALLOC(db%idir_pert, (db%npert))
 ABI_MALLOC(db%ipert_pert, (db%npert))
 ABI_MALLOC(db%qpt_pert, (3,db%npert))
 ABI_MALLOC(db%cplex_pert, (db%npert))
 ABI_MALLOC(db%ngfft3_pert, (3,db%npert))

 ! Read up to npert dimensions (last dimension is unlimited)
 call etsf_io_low_read_var(db%ncid,'idir_pert',db%idir_pert,lstat,start=[1],count=[db%npert],error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)
                                                                                                              
 call etsf_io_low_read_var(db%ncid,'ipert_pert',db%ipert_pert,lstat,start=[1],count=[db%npert],error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_read_var(db%ncid,'qpt_pert',db%qpt_pert,lstat,start=[1,1],count=[3,db%npert],error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_read_var(db%ncid,'cplex_pert',db%cplex_pert,lstat,start=[1],count=[db%npert],error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_read_var(db%ncid,'ngfft3_pert',db%ngfft3_pert,lstat,start=[1,1],count=[3,db%npert],error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)
#endif

end subroutine dfptdb_open_read
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/dfptdb_open_write
!! NAME
!!  dfptdb_open_write
!!
!! FUNCTION
!!  Open the file in write mode.
!!
!! INPUTS
!!  filename = Name of the file
!!  If.false., an IO error is raised if a file already exists.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dfptdb,scfcv3
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfptdb_open_write(db,filename,natom,nspinor,nspden,usepaw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptdb_open_write'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nspinor,nspden,usepaw 
 character(len=*),intent(in) :: filename
 type(dfptdb_t),intent(out) :: db

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat,first_write
 character(len=500) :: msg
!arrays
 type(etsf_io_low_error) :: Error

!************************************************************************

 ! Handle case in which file is present! Do we overwrite or append?
 ! append option?
 db%filename = filename
 db%npert = 0
 db%natom = natom
 db%nspinor = nspinor
 db%nspden = nspden
 db%usepaw = usepaw
 ABI_CHECK(usepaw==0,"PAW not yet coded")

 ABI_MALLOC(db%ipert_pert, (db%npert))
 ABI_MALLOC(db%idir_pert, (db%npert))
 ABI_MALLOC(db%qpt_pert, (3,db%npert))
 ABI_MALLOC(db%cplex_pert, (db%npert))
 ABI_MALLOC(db%ngfft3_pert, (3,db%npert))

 first_write = (.not.file_exists(filename))

 if (first_write) then
   call etsf_io_low_open_create(db%ncid, db%filename, etsf_file_format_version, lstat, & !& title, history, 
&  error_data=Error, with_etsf_header=.False.)
 else
   call etsf_io_low_open_modify(db%ncid, db%filename, lstat, Error_data=Error, with_etsf_header=.FALSE.)
 end if
 ETSF_CHECK_ERROR(lstat, Error)

 if (.not.first_write) return

 call etsf_io_low_set_define_mode(db%ncid, lstat, Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 ! Define the basic dimensions used in ETSF-IO files. 
 call ncid_define_basedims(db%ncid)

 ! 0 --> Unlimited dimension
 !call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )
 !call etsf_io_low_write_dim(db%ncid,'unlimited_npert',NF90_UNLIMITED,lstat,Error_data=Error)  
 call etsf_io_low_write_dim(db%ncid,'unlimited_npert',0,lstat,Error_data=Error)  
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_write_dim(db%ncid,'number_of_atoms',db%natom,lstat,Error_data=Error)  
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_write_dim(db%ncid,'number_of_spinor_components',db%nspinor,lstat,Error_data=Error)  
 ETSF_CHECK_ERROR(lstat,Error)

 !?
 call etsf_io_low_write_dim(db%ncid,'number_of_spin_density_components',db%nspden,lstat,Error_data=Error)  
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_def_var(db%ncid,"usepaw",etsf_io_low_integer,lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_def_var(db%ncid,'number_of_perturbations',etsf_io_low_integer,lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_def_var(db%ncid,"qpt_pert",etsf_io_low_double,&
&  [pad('number_of_reduced_dimensions'), pad('unlimited_npert')],lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_def_var(db%ncid,"cplex_pert",etsf_io_low_integer,['unlimited_npert'],lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_def_var(db%ncid,"idir_pert",etsf_io_low_integer,['unlimited_npert'],lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)
                                                                                                                              
 call etsf_io_low_def_var(db%ncid,"ipert_pert",etsf_io_low_integer,['unlimited_npert'],lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)
                                                                                                                              
 call etsf_io_low_def_var(db%ncid,"ngfft3_pert",etsf_io_low_integer,[pad('3'),pad('unlimited_npert')],lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_set_write_mode(db%ncid,lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 ! Write variables
 call etsf_io_low_set_write_mode(db%ncid,lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_write_var(db%ncid,'usepaw',db%usepaw,lstat,Error_data=Error)  
 ETSF_CHECK_ERROR(lstat,Error)

 !call etsf_io_low_write_var(db%ncid,"idir_pert",db%idir_pert,lstat,start=[db%npert],Error_data=Error)
 !ETSF_CHECK_ERROR(lstat,Error)
 !                                                                                                        
 !call etsf_io_low_write_var(db%ncid,"ipert_pert",db%ipert_pert,lstat,count=[db%npert],Error_data=Error)
 !ETSF_CHECK_ERROR(lstat,Error)

 !call etsf_io_low_write_var(db%ncid,"qpt_pert",db%qpt_pert,lstat,count=[3,db%npert],Error_data=Error)
 !ETSF_CHECK_ERROR(lstat,Error)

 !call etsf_io_low_write_var(db%ncid,"cplex_pert",db%cplex_pert,lstat,count=[db%npert],Error_data=Error)
 !ETSF_CHECK_ERROR(lstat,Error)

 !call etsf_io_low_write_var(db%ncid,"ngfft3_pert",db%ngfft3_pert,lstat,count=[db%npert],Error_data=Error)
 !ETSF_CHECK_ERROR(lstat,Error)
#endif

end subroutine dfptdb_open_write
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/dfptdb_write
!! NAME
!!  dfptdb_write
!!
!! FUNCTION
!!
!! INPUTS
!!  db<>
!!  rhor1(cplex*nfft,nspden)=array for RF electron density in real space in electrons/bohr**3.
!!  vtrial1(cplex*nfft,nspden)= 1st-order potential
!!  === if usepaw==1 
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!
!! PARENTS
!!      m_dfptdb,scfcv3
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfptdb_write(db,qpt,ipert,idir,cplex,nfft,nspden,ngfft3,rhor1,vtrial1,pawrhoij1,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptdb_write'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,cplex,nfft,nspden
 integer,intent(out) :: ierr
 type(dfptdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft3(3)
 real(dp),intent(in) :: qpt(3),rhor1(cplex*nfft,nspden)
 real(dp),intent(in) :: vtrial1(cplex*nfft,nspden)
 type(pawrhoij_type),intent(in) :: pawrhoij1(db%natom*db%usepaw)

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 !integer :: ncerr,qpt_id
 logical :: lstat
 character(len=500) :: msg,cplex_dimname,nfft_dimname,rhor1_varname,vtrial1_varname 
 character(len=500) :: cplexnfft_dimname
 type(etsf_io_low_error) :: Error
!arrays
 integer,allocatable :: i1_work(:),i2_work(:,:)
 real(dp),allocatable :: r2_work(:,:)

!************************************************************************

 ierr = 0

 if (dfptdb_find_pert(db,qpt,ipert,idir) /= 0) then
   write(msg,"(a,3f8.5,2(a,i0))")&
&    "Database already contains perturbation qpt: ",qpt,", ipert: ",ipert,", idir:",idir
   MSG_ERROR(msg)
 end if

 ! Add info on the perturbation in the internal tables.
 !call dfptdb_add_pert(qtp,idir,ipert)

 ! Increment the number of perturbations.
 db%npert = db%npert + 1

 call alloc_copy(db%ipert_pert, i1_work)
 ABI_FREE(db%ipert_pert)
 ABI_MALLOC(db%ipert_pert, (db%npert))
 db%ipert_pert(:db%npert-1) = i1_work
 db%ipert_pert(db%npert) = ipert
 ABI_FREE(i1_work)

 call alloc_copy(db%idir_pert, i1_work)
 ABI_FREE(db%idir_pert)
 ABI_MALLOC(db%idir_pert, (db%npert))
 db%idir_pert(:db%npert-1) = i1_work
 db%idir_pert(db%npert) = idir
 ABI_FREE(i1_work)

 call alloc_copy(db%qpt_pert, r2_work)
 ABI_FREE(db%qpt_pert)
 ABI_MALLOC(db%qpt_pert, (3,db%npert))
 db%qpt_pert(:,:db%npert-1) = r2_work
 db%qpt_pert(:,db%npert) = qpt
 ABI_FREE(r2_work)

 call alloc_copy(db%cplex_pert, i1_work)
 ABI_FREE(db%cplex_pert)
 ABI_MALLOC(db%cplex_pert, (db%npert))
 db%cplex_pert(:db%npert-1) = i1_work
 db%cplex_pert(db%npert) = cplex
 ABI_FREE(i1_work)

 call alloc_copy(db%ngfft3_pert, i2_work)
 ABI_FREE(db%ngfft3_pert)
 ABI_MALLOC(db%ngfft3_pert, (3,db%npert))
 db%ngfft3_pert(:,:db%npert-1) = i2_work
 db%ngfft3_pert(:,db%npert) = ngfft3
 ABI_FREE(i2_work)

 ! Create dimensions for this perturbation.
 call etsf_io_low_set_define_mode(db%ncid, lstat, Error)
 ETSF_CHECK_ERROR(lstat,Error)

 cplex_dimname = strcat("cplex_pertnum_",itoa(db%npert))
 call etsf_io_low_write_dim(db%ncid,cplex_dimname,cplex,lstat,Error_data=Error)  
 ETSF_CHECK_ERROR(lstat,Error)

 nfft_dimname = strcat("nfft_pertnum_",itoa(db%npert))
 call etsf_io_low_write_dim(db%ncid,nfft_dimname,nfft,lstat,Error_data=Error)  
 ETSF_CHECK_ERROR(lstat,Error)

 cplexnfft_dimname = strcat("cplexnfft_pertnum_",itoa(db%npert))
 call etsf_io_low_write_dim(db%ncid,cplexnfft_dimname,cplex*nfft,lstat,Error_data=Error)  
 ETSF_CHECK_ERROR(lstat,Error)

 ! Define new variables.
 rhor1_varname = strcat("rhor1_pertnum_",itoa(db%npert))
 call etsf_io_low_def_var(db%ncid,rhor1_varname,etsf_io_low_double,&
&  [pad(cplexnfft_dimname), pad('number_of_spin_density_components')],lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 vtrial1_varname = strcat("vtrial1_pertnum_",itoa(db%npert))
 call etsf_io_low_def_var(db%ncid,vtrial1_varname,etsf_io_low_double,&
&  [pad(cplexnfft_dimname), pad('number_of_spin_density_components')],lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 ! Write variables
 call etsf_io_low_set_write_mode(db%ncid,lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 ! Update the variable with the total number of perturbations
 call etsf_io_low_write_var(db%ncid,'number_of_perturbations',db%npert,lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 ! This one is needed to avoid problems with ETSF-IO
 call add_record(db%ncid,"ipert_pert",ipert,start=[db%npert],count=[1]) 

 call etsf_io_low_write_var(db%ncid,"ipert_pert",ipert,lstat,start=[db%npert],count=[1],Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_write_var(db%ncid,"idir_pert",idir,lstat,start=[db%npert],count=[1],Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_write_var(db%ncid,"qpt_pert",qpt,lstat,start=[1,db%npert],count=[3,1],Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_write_var(db%ncid,"cplex_pert",cplex,lstat,start=[db%npert],count=[1],Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_write_var(db%ncid,"ngfft3_pert",ngfft3(1:3),lstat,start=[1,db%npert],count=[3,1],Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 ! Write rhor1 and vtrial1
 call etsf_io_low_write_var(db%ncid,rhor1_varname,rhor1,lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)
                                                                                    
 call etsf_io_low_write_var(db%ncid,vtrial1_varname,vtrial1,lstat,Error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 ! Write pawrhoij1
 if (db%usepaw==1) then
!   call pawrhoij_io(pawrhoij,unitfi,nsppol_in,nspinor_in,nspden_in,nlmn_type,typat,&
!&                   headform,rdwr_mode,form,natinc,mpi_atmtab)
 end if
#endif

end subroutine dfptdb_write
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/dfptdb_read
!! NAME
!!  dfptdb_read
!!
!! FUNCTION
!!
!! INPUTS
!!  db<>
!!  rhor1(cplex*nfft,nspden)=array for RF electron density in real space in electrons/bohr**3.
!!  vtrial1(cplex*nfft,nspden)= 1st-order potential
!!  === if usepaw==1 
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!
!! PARENTS
!!      m_dfptdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfptdb_read(db,qpt,ipert,idir,cplex,nfft,nspden,ngfft3,rhor1,vtrial1,pawrhoij1,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptdb_read'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,cplex,nfft,nspden
 integer,intent(out) :: ierr
 type(dfptdb_t),intent(in) :: db
!arrays
 integer,intent(in) :: ngfft3(3)
 real(dp),intent(in) :: qpt(3)
 real(dp),target,intent(out) ::rhor1(cplex*nfft,nspden),vtrial1(cplex*nfft,nspden)
 type(pawrhoij_type),allocatable,intent(out) :: pawrhoij1(:) ! natom*usepaw

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 integer :: ip,file_cplex,file_nfft
 logical :: lstat,use_fourier
 character(len=500) :: msg !,cplex_dimname,nfft_dimname,rhor1_varname,vtrial1_varname
 type(etsf_io_low_error) :: Error
!arrays
 integer :: file_ngfft3(3)
 real(dp), ABI_CONTIGUOUS pointer :: work_rhor1(:,:),work_vtrial1(:,:)

!************************************************************************

 ip = dfptdb_find_pert(db,qpt,ipert,idir)
 
 ierr = 0
 if (ip==0) then
   write(msg,"(a,3f8.5,2(a,i0))")&
&    "Database does not contain perturbation qpt: ",qpt,", ipert: ",ipert,", idir:",idir
   MSG_WARNING(msg)
   ierr = 1
   return 
 end if

 file_cplex = db%cplex_pert(ip)
 file_ngfft3 = db%ngfft3_pert(1:3,ip)
 file_nfft = product(file_ngfft3(1:3))

 work_rhor1 => rhor1
 work_vtrial1 => vtrial1
 use_fourier = (any(file_ngfft3 /= ngfft3(1:3)))

 if (use_fourier) then
   ! allocate data with correct shape, then use Fourier interpolation
   ! to have rhor1 and vtrial1 on the input ngfft(1:3) FFT mesh.
   MSG_ERROR("Will perform FFT interpolation of vtrial1 and rhor1")
   ABI_MALLOC(work_rhor1, (file_cplex*file_nfft,nspden))
   ABI_MALLOC(work_vtrial1, (file_cplex*file_nfft,nspden))
 end if

 ! Read rhor1 and vtrial1
 call etsf_io_low_read_var(db%ncid,strcat("rhor1_pertnum_",itoa(db%npert)),work_rhor1,lstat,error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 call etsf_io_low_read_var(db%ncid,strcat("vtrial1_pertnum_",itoa(db%npert)),work_vtrial1,lstat,error_data=Error)
 ETSF_CHECK_ERROR(lstat,Error)

 if (use_fourier) then
   !call fourier_interpol(cplex,nspden,optin,optout,nfft_in,ngfft_in,nfft_out,ngfft_out,&
   !& paral_kgb,MPI_enreg,rhor_in,rhor_out,rhog_in,rhog_out)
   ABI_FREE(work_rhor1)
   ABI_FREE(work_vtrial1)
 end if

 if (db%usepaw==1) then
    ABI_DT_MALLOC(pawrhoij1, (db%natom*db%usepaw))
! TODO need netcdf fmethod for pawrhoij

!call pawrhoij_alloc(pawrhoij,cplex,nspden,nspinor,nsppol,typat,&              ! Mandatory arguments
!&                      lmnsize,ngrhoij,nlmnmix,pawtab,use_rhoij_,use_rhoijp,& ! Optional arguments
!&                      use_rhoijres,mpi_comm_atom,mpi_atmtab)                 ! Optional arguments

!   call pawrhoij_io(pawrhoij,unitfi,nsppol_in,nspinor_in,nspden_in,nlmn_type,typat,&
!&                   headform,rdwr_mode,form,natinc,mpi_atmtab)
 end if
#endif

end subroutine dfptdb_read
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/dfptdb_merge
!! NAME
!!  dfptdb_merge
!!
!! FUNCTION
!!  merge multiple databases
!!
!! INPUT 
!!  delete=True if files should be deleted after the merge
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfptdb_merge(odb,out_filename,filenames,delete)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptdb_merge'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: delete
 character(len=*),intent(in) :: out_filename
 type(dfptdb_t),intent(out) :: odb
!arrays
 character(len=*),intent(in) :: filenames(:)

!Local variables-------------------------------
#ifdef HAVE_TRIO_ETSF_IO
!scalars
 integer :: ii,ierr,ip,nfiles,cplex,nfft,nspinor,nspden,natom,usepaw,idir,ipert
 logical :: lstat
 type(etsf_io_low_error) :: Error
!arrays
 integer :: ngfft3(3)
 real(dp) :: qpt(3)
 real(dp),allocatable :: rhor1(:,:),vtrial1(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij1(:)
 type(dfptdb_t),allocatable :: idbs(:)

! *************************************************************************
 
 nfiles = size(filenames)

 ABI_DT_MALLOC(idbs, (nfiles))
 do ii=1,nfiles
   call dfptdb_open_read(idbs(ii),filenames(ii))
 end do

 ! Consistency check
 !call dfptdb_validate(idbs)

 ! Get dimensions and write new database.
 natom = idbs(1)%natom
 usepaw = idbs(1)%usepaw
 nspinor = idbs(1)%nspinor 
 nspden = idbs(1)%nspden

 call dfptdb_open_write(odb,out_filename,natom,nspinor,nspden,usepaw)

 do ii=1,nfiles
   do ip=1,idbs(ii)%npert
     ! Get info on the perturbation.
     qpt = idbs(ii)%qpt_pert(:,ip)
     idir = idbs(ii)%idir_pert(ip)
     ipert = idbs(ii)%ipert_pert(ip)
     cplex = idbs(ii)%cplex_pert(ip)
     ngfft3 = idbs(ii)%ngfft3_pert(1:3,ip)
     nfft = product(ngfft3(1:3))

     ABI_MALLOC(rhor1, (cplex*nfft,nspden))
     ABI_MALLOC(vtrial1, (cplex*nfft,nspden))

     call dfptdb_read(idbs(ii),qpt,ipert,idir,cplex,nfft,nspden,ngfft3,rhor1,vtrial1,pawrhoij1,ierr)

     call dfptdb_write(odb,qpt,ipert,idir,cplex,nfft,nspden,ngfft3,rhor1,vtrial1,pawrhoij1,ierr)

     ABI_FREE(rhor1)
     ABI_FREE(vtrial1)

     if (usepaw == 1) then
       call pawrhoij_destroy(pawrhoij1)
       ABI_DT_FREE(pawrhoij1)
     end if
   end do
 end do

 ! Close the temporary databases.
 do ii=1,nfiles
   call dfptdb_close(idbs(ii),delete=delete)
 end do
 ABI_DT_FREE(idbs)
#endif

end subroutine dfptdb_merge
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/add_record_int0d
!! NAME
!!  add_record_int0d
!!
!! FUNCTION
!!
!! INPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine add_record_int0d(ncid,varname,value,start,count)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'add_record_int0d'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,value
 character(len=*),intent(in) :: varname
!arrays
 integer,intent(in) :: start(:),count(:)

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 integer :: work(1)

! *************************************************************************

 work(1) = value
 call add_record_int1d(ncid,varname,work,start,count)
#endif

end subroutine add_record_int0d
!!***

!----------------------------------------------------------------------

!!****f* m_dfptdb/add_record_int1d
!! NAME
!!  add_record_int1d
!!
!! FUNCTION
!!
!! INPUT
!!
!! PARENTS
!!      m_dfptdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine add_record_int1d(ncid,varname,value,start,count)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'add_record_int1d'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 character(len=*),intent(in) :: varname
!arrays
 integer,intent(in) :: value(:),start(:),count(:)

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 integer :: ncerr,var_id

! *************************************************************************

 ncerr = nf90_inq_varid(ncid, "qpt_pert", var_id)
 ABI_CHECK(ncerr == NF90_NOERR, sjoin("inquring variable: ",varname))

 ncerr = nf90_put_var(ncid,var_id,value,start=start, count=count)
 ABI_CHECK(ncerr == NF90_NOERR, sjoin("putting variable: ",varname))
#endif

end subroutine add_record_int1d
!!***

#endif

!----------------------------------------------------------------------

END MODULE m_dfptdb
!!***
