!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_wfk
!! NAME
!!  m_wfk
!!
!! FUNCTION
!!  This module defines the file handler wfk_t used to perform IO operations on the 
!!  wavefunction files produced by ABINIT.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! NOTES
!!  1) The object supports random access also when plain Fortran IO is used.
!!     One can easily *read* the block of wavefunctions with a given (kpt,spin) 
!!     by simply passing the appropriate indices (ik_ibz,spin) to the wfk_read_ routines.
!!     Note however that the same feature is not available in write mode when Fortran IO
!!     is used. In this case indeed one should access the block of wavefunctions 
!!     according to their (kpt,spin) indices in order to write the correct record markers 
!!     MPI-IO and NETCDF do not have such limitation.
!!
!!  2) MPI-IO reads are done with file views even for contiguous data of the same type.
!!     I found, indeed, that mixing views with explicit offset calls causes
!!     wrong results unless the file is closed and re-open! Very strange since, according to
!!     the documentation, the two APIs do not interfere and can be mixed. Calls to MPI_FILE_SEEK
!!     to reset the pointer to the start of the file do not solve the problem. Don't know if it's
!!     a feature or a bug (the problem showed up with MPICH2, I haven't tested other MPI libraries)
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

MODULE m_wfk

 use defs_basis
 use m_profiling_abi
 use m_errors
#ifdef HAVE_MPI2 
 use mpi
#endif
 use m_xmpi
 use m_mpiotk
 use m_header
 use m_wffile
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 use defs_abitypes,  only : hdr_type, MPI_type 
 use m_time,         only : cwtime
 use m_fstrings,     only : strcat
 use m_io_tools,     only : get_unit, mvrecord, iomode_from_fname, open_file, close_unit
 use m_fftcore,      only : get_kg
 use m_distribfft,   only : init_distribfft_seq
 use m_mpinfo,       only : destroy_mpi_enreg 

 implicit none

#ifdef HAVE_MPI1 
 include 'mpif.h'
#endif

 private

 integer,private,parameter :: WFK_NOMODE    = 0
 integer,private,parameter :: WFK_READMODE  = 1
 integer,private,parameter :: WFK_WRITEMODE = 2
!!***

!----------------------------------------------------------------------

!!****t* m_wfk/wfk_t
!! NAME
!!  wfk_t
!!
!! FUNCTION
!!  File handler for the WFK file.
!!
!! SOURCE

 type,public :: wfk_t

  integer :: fh
   !  unit number if IO_MODE_FORTRAN
   !  MPI file handler if IO_MODE_MPI
   !  Netcdf file handler if IO_MODE_NETCDF 

  integer :: iomode
   ! Method used to access the WFK file:
   !   IO_MODE_FORTRAN for usual Fortran IO routines
   !   IO_MODE_MPI if MPI/IO routines.
   !   IO_MODE_ETSF, NetCDF format read/written via etsf-io.

  integer :: mband
  ! Max number of bands stored on file (MAX(Hdr%nband))

  integer :: nkpt
  ! Number of k-points.

  integer :: nsppol
  ! Number of spins

  integer :: nspinor
  ! Number of spinor components.

  integer :: formeig
   ! formeig=format of the eigenvalues
   !    0 => vector of eigenvalues
   !    1 => hermitian matrix of eigenvalues
   ! TODO: this should be reported somewhere in the WFK file, at present is passed to wfk_open

  integer :: fform
   ! File type format of the header

  integer :: rw_mode = WFK_NOMODE
   ! (Read|Write) mode 

  character(len=fnlen) :: fname = ABI_NOFILE
   ! File name 

  integer :: master
   ! master node of the IO procedure 

  integer :: my_rank
   ! index of my processor in the MPI communicator comm

  integer :: nproc
   ! number of processors in comm

  integer :: comm
   ! MPI communicator 

  integer :: recn_eof
   ! EOF record number (used for Fortran IO)

  integer(XMPI_OFFSET_KIND) :: offset_eof
  ! EOF offset (used for MPI-IO access)

  logical :: debug=.FALSE.
  !logical :: debug=.TRUE.

  type(hdr_type) :: Hdr
   ! Abinit header.

  integer,allocatable :: nband(:,:) 
  ! nband(nkpt,nsppol) = Number of bands at each (k,s)

  integer :: f90_fptr(3) = (/0,0,0/)
  ! The position of the file pointer used for sequential access with Fortran-IO.
  !  f90_fprt(1) = Index of the k-point associated to the block.
  !  f90_fprt(2) = the spin associated to the block.
  !  f90_fprt(3) = Record Type (see REC_* variables).
  !  (/0,0,0/) corresponds to the beginning of the file.
  !  FPTR_EOF signals the end of file

  integer,allocatable :: recn_ks(:,:,:)
   ! recn_ks(k,s,1) : record number of  (npw, nspinor, nband_disk)
   ! recn_ks(k,s,2) : record number of the (k+G) vectors.
   ! recn_ks(k,s,3) : record number of the eigenvalues.
   ! recn_ks(k,s,4) : record number of the first wavefunction in  the wf coefficients block.

  integer(XMPI_OFFSET_KIND),allocatable :: offset_ks(:,:,:) 
   ! offset_ks(k,s,1) : offset of the record: npw, nspinor, nband_disk.
   ! offset_ks(k,s,2) : offset of the Second record: (k+G) vectors.
   ! offset_ks(k,s,3) : offset of the third record eigenvalues.
   ! offset_ks(k,s,4) : offset of the fourth record (wavefunction coefficients).
   ! NB: The offset point to the Fortran record marker and not to the data.

  integer(XMPI_OFFSET_KIND) :: hdr_offset
   ! offset of the header
   ! TODO this should be the output of a hdr method!

  integer(XMPI_OFFSET_KIND) :: chunk_bsize
   ! IO is performed in chunks of max size chunk_bsize [bytes] 

 end type wfk_t

!public procedures.
 public :: wfk_open_read           ! Open the WFK file in read mode.
 public :: wfk_open_write          ! Open the WFK file in write mode.
 public :: wfk_close               ! Close the WFK file and release the memory allocated in wfk_t.

 public :: wfk_read_band_block     ! Read a contiguous block of bands for a given (kpoint, spin)
 public :: wfk_write_band_block    ! Write a contiguous block of bands for a given (kpoint, spin)
 public :: wfk_read_bmask          ! Read a scattered set of bands for a given (kpoint, spin).
 !public :: wfk_write_bmask        ! Write a scattered set of bands for a given (kpoint, spin).

 public :: wfk_read_eigk           ! Read the eigenvalues at a given (kpoint,spin).
 public :: wfk_read_eigenvalues    ! Read all the GS eigenvalues stored in the WFK file fname.
 !public :: wfk_read_gkk           ! Write all the GS eigenvalues stored in the WFK file fname.
 !public :: wfk_write_eigk         ! Write the eigenvalues at a given (kpoint,spin).
 !public :: wfk_write_eigenvalues  ! Write all the GS eigenvalues.
 public :: wfk_write_allgkk        ! Write all the GKK matrix elements.
 public :: wfk_cp

 ! Unitary tests and profiling tools
 public :: wfk_prof                ! Profiling tool.
 public :: wfk_diff                ! Compare two WFK file for binary equality.
 public :: wfk_create_wfkfile      ! Create a FAKE WFK file.
 public :: wfk_check_wfkfile       ! Read a FAKE WFK file and perform basic tests.
!!***

! Indices associated to the start of the different records of the WFK file.
 integer,private,parameter :: REC_HDR=0
 integer,private,parameter :: REC_NPW=1
 integer,private,parameter :: REC_KG =2
 integer,private,parameter :: REC_EIG=3
 integer,private,parameter :: REC_CG =4
 integer,private,parameter :: REC_NUM=REC_CG

 !integer(XMPI_OFFSET_KIND),private,parameter :: WFK_ORIGIN = 0

 integer,private,parameter :: FPTR_EOF(3) = (/-1,-1,-1/)

 integer(XMPI_OFFSET_KIND),private,parameter :: WFK_CHUNK_BSIZE = 2000 * (1024.0_dp**2)
   ! Maximum size (in bytes) of the block of wavefunctions that are (read|written) 
   ! in a single MPI-IO call. (Some MPI-IO implementation crashes if we try to 
   ! (read|write) a big chunk of data with a single call.

!----------------------------------------------------------------------

!!****t* m_wfk/kvars_t
!! NAME
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: kvars_t
   integer,allocatable  :: kg_k(:,:)
   real(dp),pointer :: occ_k(:)   => null()
   real(dp),pointer :: eig_k(:)   => null()
 end type kvars_t

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_open_read
!! NAME
!!  wfk_open_read
!!
!! FUNCTION
!!  Open the WFK file in read mode.
!!
!! INPUTS
!!  fname = Name of the file
!!  formeig = 0 for GS wavefunctions, 1 for RF wavefunctions.
!!  iomode = access mode 
!!  funt = Fortran unit numer for  Only used if iomode == IO_MODE_FORTRAN
!!  comm = MPI communicator (used for MPI-IO)
!!
!! OUTPUT
!!  Wfk<type(wfk_t)> = WFK handler initialized and set in read mode
!!  [Hdr_out]=Copy of the abinit header
!!
!! PARENTS
!!      initwf,m_wfk,m_wfs,mrggkk,wffile
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_open_read(Wfk,fname,formeig,iomode,funt,comm,Hdr_out)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_open_read'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,comm,formeig,funt
 character(len=*),intent(in) :: fname
 type(wfk_t),intent(inout) :: Wfk !vz_i
 type(hdr_type),optional,intent(inout) :: Hdr_out  ! should be intent(out), but psc miscompiles the call!

!Local variables-------------------------------
!scalars
 integer :: ierr,mpierr
 character(len=500) :: msg
#ifdef HAVE_MPI_IO
 integer :: fform 
#endif
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(etsf_io_low_error) :: Error
#endif

!************************************************************************

 DBG_ENTER("COLL")

 !Initialize the mandatory data of the Wfk datastructure
 !@wfk_t
 Wfk%rw_mode     = WFK_READMODE
 Wfk%chunk_bsize = WFK_CHUNK_BSIZE

 Wfk%fname     = fname
 Wfk%formeig   = formeig
 Wfk%iomode    = iomode
 Wfk%comm      = comm

 Wfk%master    = 0
 Wfk%my_rank   = xcomm_rank(comm)
 Wfk%nproc     = xcomm_size(comm)
 !
 ! Reads fform and the Header.
 call hdr_read_from_fname(Wfk%Hdr,fname,Wfk%fform,comm) 
 ABI_CHECK(Wfk%fform/=0,"fform ==0")

 if (Wfk%debug) then
   call hdr_echo(Wfk%Hdr,Wfk%fform,4,unit=std_out)
 end if
 !
 ! Copy the header if required.
 if (PRESENT(Hdr_out)) then
   call hdr_copy(Wfk%Hdr,Hdr_out)
 end if
 !
 ! Useful dimensions
 Wfk%mband   = MAXVAL(Wfk%Hdr%nband)
 Wfk%nkpt    = Wfk%Hdr%nkpt
 Wfk%nsppol  = Wfk%Hdr%nsppol
 Wfk%nspinor = Wfk%Hdr%nspinor

 ABI_MALLOC(Wfk%nband, (Wfk%nkpt,Wfk%nsppol))
 Wfk%nband = RESHAPE(Wfk%Hdr%nband, (/Wfk%nkpt,Wfk%nsppol/))

 ierr=0
 SELECT CASE (iomode)

 CASE (IO_MODE_FORTRAN) 
   ! All processors see a local Fortran binary file.
   ! Each node opens the file, skip the header and set f90_fptr.
   Wfk%fh = funt
   if (open_file(Wfk%fname,msg,unit=Wfk%fh,form="unformatted", status="old", action="read") /= 0) then
     MSG_ERROR(msg)
   end if

   ! Precompute number of records for Fortran IO.
   call wfk_compute_offsets(Wfk)

   call hdr_skip(Wfk%fh,ierr)
   Wfk%f90_fptr = (/1,1,REC_NPW/)

 CASE (IO_MODE_MPI) 
#ifdef HAVE_MPI_IO
   call MPI_FILE_OPEN(Wfk%comm, Wfk%fname, MPI_MODE_RDONLY, MPI_INFO_NULL, Wfk%fh, mpierr)
   ABI_CHECK_MPI(mpierr,"MPI_FILE_OPEN")

   !call MPI_FILE_SET_VIEW(Wfk%fh,origin,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,mpierr)

   call hdr_mpio_skip(Wfk%fh,fform,Wfk%hdr_offset)

   ! Precompute offsets for MPI-IO access 
   if (Wfk%hdr_offset > 0) then
     call wfk_compute_offsets(Wfk)
   else
     MSG_ERROR("hdr_offset <=0")
   end if

#else
   MSG_ERROR("MPI-IO not enabled")
#endif

 CASE (IO_MODE_ETSF)
#ifdef HAVE_TRIO_ETSF_IO
   call etsf_io_low_open_read(Wfk%fh, Wfk%fname, lstat, Error_data= Error, with_etsf_header=.FALSE.)
   ETSF_CHECK_ERROR(lstat, Error)

  ! TODO: 1) should become a method.
  !       2) ETSF_IO does not export it! 
  ! CHECK whether the file is ETSF compliant
  !if  (Wfk%iomode == IO_MODE_ETSF)
  !call etsf_io_file_check_wavefunctions_data(Wfk%fh, lstat, Error)
  !ETSF_CHECK_ERROR(lstat, Error)
  !endif

#else
   MSG_ERROR("ETSF-IO not enabled")
#endif

 CASE DEFAULT
   write(msg,'(a,i0)')' Wrong value of iomode = ',iomode
   MSG_ERROR(msg)
 END SELECT

 DBG_EXIT("COLL")

end subroutine wfk_open_read
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_open_write
!! NAME
!!  wfk_open_write
!!
!! FUNCTION
!!  Open the WFK file in write mode.
!!
!! INPUTS
!!  fname = Name of the file
!!  formeig = 0 for GS wavefunctions, 1 for RF wavefunctions.
!!  iomode = access mode 
!!  funt = Fortran unit numer for  Only used if iomode == IO_MODE_FORTRAN
!!  comm = MPI communicator (used for MPI-IO)
!!
!! OUTPUT
!!  Wfk<type(wfk_t)> = WFK handler initialized and set in read mode
!!
!! PARENTS
!!      dfpt_write_cg,m_wfk,m_wfs
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_open_write(Wfk,Hdr,fname,formeig,iomode,funt,comm,write_hdr,write_frm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_open_write'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,comm,formeig,funt
 character(len=*),intent(in) :: fname
 logical,optional,intent(in) :: write_hdr,write_frm
 type(hdr_type),intent(in) :: Hdr
 type(wfk_t),intent(out) :: Wfk

!Local variables-------------------------------
!scalars
 integer :: mpierr,ierr
 real(dp) :: cpu,wall,gflops
 logical :: do_write_frm,do_write_hdr
 character(len=500) :: msg
#ifdef HAVE_MPI_IO
 integer :: fform,nfrec,sc_mode
 integer(XMPI_OFFSET_KIND) :: offset
 integer(XMPI_OFFSET_KIND),pointer :: bsize_frecords(:)
#endif
#ifdef HAVE_TRIO_ETSF_IO
 integer :: ivar,mpw,flags_electrons
 logical :: lstat
 !character(len=80) :: file_title
 !type(etsf_groups_flags) :: flags
 type(etsf_io_low_error) :: Error
#endif
!arrays

!************************************************************************

 DBG_ENTER("COLL")

 do_write_hdr = .TRUE.; if (PRESENT(write_hdr)) do_write_hdr = write_hdr
 do_write_frm = .TRUE.; if (PRESENT(write_frm)) do_write_frm = write_frm

 !Initialize the mandatory data of the Wfk datastructure
 !@wfk_t
 Wfk%rw_mode     = WFK_WRITEMODE
 Wfk%chunk_bsize = WFK_CHUNK_BSIZE

 Wfk%fname     = fname
 Wfk%formeig   = formeig
 Wfk%iomode    = iomode
 Wfk%comm      = comm

 Wfk%master    = 0
 Wfk%my_rank   = xcomm_rank(comm)
 Wfk%nproc     = xcomm_size(comm)
 Wfk%fform     = hdr_ftype2fform(HDR_WF_PW) ! latest file format for the WFK file (planewave case)
 ABI_CHECK(Wfk%fform/=0,"hdr_ftype2fform")  ! TODO wavelets
 !write(std_out,*)"fform = ",Wfk%fform

 ! Copy the header
 call hdr_copy(Hdr,Wfk%Hdr)
 !
 ! Master writes fform and the Header with Fortran IO
 if (Wfk%my_rank==Wfk%master .and. do_write_hdr) then
   call hdr_write_to_fname(Wfk%Hdr,Wfk%fname,Wfk%fform)
   !if (Wfk%debug) then
   !  call hdr_echo(Wfk%Hdr,Wfk%fform,4,unit=std_out)
   !end if
 end if
 !
 ! Useful dimensions
 Wfk%mband   = MAXVAL(Wfk%Hdr%nband)
 Wfk%nkpt    = Wfk%Hdr%nkpt
 Wfk%nsppol  = Wfk%Hdr%nsppol
 Wfk%nspinor = Wfk%Hdr%nspinor

 ABI_MALLOC(Wfk%nband, (Wfk%nkpt,Wfk%nsppol))
 Wfk%nband = RESHAPE(Wfk%Hdr%nband, (/Wfk%nkpt,Wfk%nsppol/))

 call xmpi_barrier(Wfk%comm)

 ierr=0
 SELECT CASE (iomode)

 CASE (IO_MODE_FORTRAN) 
   ! All processors see a local Fortran binary file.
   ! Each node opens the file, skip the header and set f90_fptr.
   Wfk%fh = funt
   if (open_file(Wfk%fname,msg,unit=Wfk%fh,form="unformatted", status="unknown", action="readwrite") /= 0) then
     MSG_ERROR(msg)
   end if

   ! Precompute number of records for Fortran IO.
   call wfk_compute_offsets(Wfk)

   call hdr_skip(Wfk%fh,ierr)
   Wfk%f90_fptr = (/1,1,REC_NPW/)

 CASE (IO_MODE_MPI) 

#ifdef HAVE_MPI_IO
   ! FIXME: mode flags should be rationalized
   !call MPI_FILE_OPEN(Wfk%comm, Wfk%fname, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, Wfk%fh, mpierr)
   call cwtime(cpu,wall,gflops,"start")

   call MPI_FILE_OPEN(Wfk%comm, Wfk%fname, MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, Wfk%fh, mpierr)
   ABI_CHECK_MPI(mpierr,"MPI_FILE_OPEN")

   !call MPI_FILE_SET_VIEW(Wfk%fh,origin,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,mpierr)

   ! TODO
   !%% call MPI_File_set_size(Wfk%fh, MPI_Offset size, mpierr)
   !ABI_CHECK_MPI(mpierr,"MPI_FILE_SET_SIZE")

   call hdr_mpio_skip(Wfk%fh,fform,Wfk%hdr_offset)
   ABI_CHECK(fform == Wfk%fform,"fform != Wfk%fform")
   !call hdr_io(Wfk%fform,Wfk%Hdr,4,std_out)

   ! Precompute offsets for MPI-IO access 
   if (Wfk%hdr_offset > 0) then
     call wfk_compute_offsets(Wfk)
   else
     MSG_ERROR("hdr_offset <=0")
   end if

   call cwtime(cpu,wall,gflops,"stop")
   write(msg,"(2(a,f8.2))")" FILE_OPEN cpu: ",cpu,", wall:",wall
   call wrtout(std_out,msg,"PERS")

   ! Write Fortran record markers.
   !if (.FALSE.) then
   if (do_write_frm) then
     call cwtime(cpu,wall,gflops,"start")
     call hdr_bsize_frecords(Wfk%Hdr,Wfk%formeig,nfrec,bsize_frecords)

     sc_mode = xmpio_collective
     offset = Wfk%hdr_offset 

     if (sc_mode == xmpio_collective) then
       call xmpio_write_frmarkers(Wfk%fh,offset,sc_mode,nfrec,bsize_frecords,ierr)
     else 
       ierr = 0
       if (Wfk%my_rank == Wfk%master) then
         call xmpio_write_frmarkers(Wfk%fh,offset,xmpio_single,nfrec,bsize_frecords,ierr)
       end if
     end if
     ABI_CHECK(ierr==0,"xmpio_write_frmarkers returned ierr!=0")

     !call MPI_FILE_SYNC(Wfk%fh,mpierr)
     !ABI_CHECK_MPI(mpierr,"FILE_SYNC")

     if (Wfk%debug) then
       call xmpio_check_frmarkers(Wfk%fh,offset,sc_mode,nfrec,bsize_frecords,ierr)
       ABI_CHECK(ierr==0,"xmpio_check_frmarkers returned ierr!=0")
     end if

     ABI_FREE(bsize_frecords)

     call cwtime(cpu,wall,gflops,"stop")
     write(msg,"(2(a,f8.2))")" write_frmarkers cpu: ",cpu,", wall:",wall
     call wrtout(std_out,msg,"PERS")
   end if
#else
   MSG_ERROR("MPI-IO not enabled")
#endif

 CASE (IO_MODE_ETSF)
#ifdef HAVE_TRIO_ETSF_IO
   call etsf_io_low_open_modify(Wfk%fh, Wfk%fname, lstat, Error_data=Error, with_etsf_header=.FALSE.)
   ETSF_CHECK_ERROR(lstat, Error)

   !call etsf_io_low_write_att(ncid, ivar, "k_dependent", "yes", lstat, error_data = error_data)
   !file_title = "Wavefunctions file"
   !flags%geometry  = etsf_geometry_all
   !flags%kpoints   = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights
   !flags%electrons = etsf_electrons_all - etsf_electrons_x_functional - etsf_electrons_c_functional
   !flags%basisdata = etsf_basisdata_basis_set
   !flags%basisdata = flags%basisdata + etsf_basisdata_kin_cutoff + etsf_basisdata_n_coeff
   !flags%basisdata = etsf_basisdata_red_coord_pw
   !call etsf_io_basisdata_def(Wfk%fh, lstat, error, k_dependent=.TRUE., flags=flags%basisdata)
   !ETSF_CHECK_ERROR(lstat, error)

   call etsf_io_low_write_dim(Wfk%fh,"real_or_complex_coefficients",2,lstat,error_data=Error)
   ETSF_CHECK_ERROR(lstat,Error)

   mpw = MAXVAL(Hdr%npwarr)
   call etsf_io_low_write_dim(Wfk%fh,"max_number_of_coefficients",mpw,lstat,error_data=Error)
   ETSF_CHECK_ERROR(lstat,Error)

   ! Define kg_k
   call etsf_io_low_def_var(Wfk%fh, "reduced_coordinates_of_plane_waves",etsf_io_low_integer,&
&    (/ pad("number_of_reduced_dimensions"), pad("max_number_of_coefficients"), pad("number_of_kpoints") /),&
&    lstat, ncvarid=ivar, error_data=Error)
   ETSF_CHECK_ERROR(lstat,Error)

   call etsf_io_low_write_att(Wfk%fh, ivar, "k_dependent", "yes", lstat, error_data=error)
   ETSF_CHECK_ERROR(lstat,Error)

   ! this comes from
   !subroutine etsf_io_main_def(ncid, lstat, error_data, k_dependent, flags, split)

   flags_electrons = etsf_electrons_eigenvalues 

   call etsf_io_electrons_def(Wfk%fh, lstat, Error, k_dependent=.TRUE., flags=flags_electrons)
   ETSF_CHECK_ERROR(lstat,Error)

#if 1
   call etsf_io_low_def_var(Wfk%fh, "coefficients_of_wavefunctions",etsf_io_low_double, &
&    (/pad("real_or_complex_coefficients"), pad("max_number_of_coefficients"), pad("number_of_spinor_components"),&
&      pad("max_number_of_states"),pad("number_of_kpoints"), pad("number_of_spins")/), lstat, ncvarid=ivar, error_data=error)
   ETSF_CHECK_ERROR(lstat,Error)

#else
   flags%main = etsf_main_wfs_coeff
   call etsf_io_main_def(Wfk%fh, lstat, error, k_dependent = .TRUE., flags = flags%main)
   ETSF_CHECK_ERROR(lstat, error)
#endif

   ! Switch to write mode.
   call etsf_io_low_set_write_mode(Wfk%fh,lstat,Error_data=Error)
   ETSF_CHECK_ERROR(lstat,Error)

#else
   MSG_ERROR("ETSF-IO not enabled")
#endif

 CASE DEFAULT
   write(msg,'(a,i0)')' Wrong value of iomode = ',iomode
   MSG_ERROR(msg)
 END SELECT

 DBG_EXIT("COLL")

end subroutine wfk_open_write
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_close
!! NAME
!!  wfk_close
!!
!! FUNCTION
!!  Close the wavefunction file handler and release the memory allocated
!!
!! PARENTS
!!      dfpt_write_cg,initwf,m_wfk,m_wfs,mrggkk,wffile
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_close(Wfk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_close'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfk_t),intent(inout) :: Wfk

!Local variables-------------------------------
!scalars
 integer :: ierr
 character(len=500) :: msg
#ifdef HAVE_MPI_IO
 integer :: mpierr,nfrec
 integer(XMPI_OFFSET_KIND),pointer :: bsize_frecords(:)
#endif
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(etsf_io_low_error) :: Error
#endif

! *************************************************************************

 DBG_ENTER("COLL")

 !@wfk_t
 Wfk%rw_mode = WFK_NOMODE

 SELECT CASE (Wfk%iomode)

 CASE (IO_MODE_FORTRAN) 
    ABI_FCLOSE(Wfk%fh,msg)

 CASE (IO_MODE_MPI)
#ifdef HAVE_MPI_IO
   call MPI_FILE_CLOSE(Wfk%fh,mpierr)
   ABI_CHECK_MPI(mpierr,"FILE_CLOSE!")

   if (Wfk%debug .and. Wfk%my_rank == Wfk%master) then
     ! Check the fortran records.
     call MPI_FILE_OPEN(xmpi_self, Wfk%fname, MPI_MODE_RDONLY, MPI_INFO_NULL, Wfk%fh, mpierr)
     ABI_CHECK_MPI(mpierr,"MPI_FILE_OPEN")
     call hdr_bsize_frecords(Wfk%Hdr,Wfk%formeig,nfrec,bsize_frecords)
     call xmpio_check_frmarkers(Wfk%fh,Wfk%hdr_offset,xmpio_single,nfrec,bsize_frecords,ierr)
     ABI_CHECK(ierr==0,"xmpio_check_frmarkers returned ierr!=0")
     ABI_FREE(bsize_frecords)
     call MPI_FILE_CLOSE(Wfk%fh,mpierr)
     ABI_CHECK_MPI(mpierr,"FILE_CLOSE!")
   end if
#endif

 CASE (IO_MODE_ETSF)
#ifdef HAVE_TRIO_ETSF_IO
   call etsf_io_low_close(Wfk%fh, lstat, Error_data=Error)
   ETSF_CHECK_ERROR(lstat, Error)
#endif

 CASE DEFAULT
   write(msg,'(a,i0)')' Wrong value of iomode = ',Wfk%iomode
   MSG_ERROR(msg)
 END SELECT
 !
 ! Free memory.
 call hdr_free(Wfk%Hdr)

 if (allocated(Wfk%nband)) then
   ABI_FREE(Wfk%nband)
 end if

 if (allocated(Wfk%recn_ks)) then
   ABI_FREE(Wfk%recn_ks)
 end if

 if (allocated(Wfk%offset_ks)) then
   ABI_FREE(Wfk%offset_ks)
 end if

 DBG_EXIT("COLL")

end subroutine wfk_close
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_band_block
!! NAME
!!  wfk_read_band_block
!!
!! FUNCTION
!!  Read a block of contigous bands.
!!
!! INPUTS
!!  Wfk<type(wfk_t)>=
!!  band_block(2)=Initial and final band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  [kg_k=(:,:)] = G-vectors
!!  [eig_k(:)] = Eigenvectors
!!  [cg_k(:,:)]  = Fourier coefficients
!!
!! PARENTS
!!      initwf,m_wfk,m_wfs,wffile
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_read_band_block(Wfk,band_block,ik_ibz,spin,sc_mode,kg_k,cg_k,eig_k,occ_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_read_band_block'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,sc_mode 
 type(wfk_t),intent(inout) :: Wfk
!arrays
 integer,intent(in) :: band_block(2)
 integer,intent(out), DEV_CONTARRD  optional,target :: kg_k(:,:)  !(3,npw_k)
 real(dp),intent(out), DEV_CONTARRD optional,target :: cg_k(:,:) !(2,npw_k*nspinor*nband)
 real(dp),intent(out),optional,target :: eig_k((2*Wfk%mband)**Wfk%formeig*Wfk%mband)
 real(dp),intent(out),optional,target :: occ_k(Wfk%mband)  

!Local variables-------------------------------
!scalars
 integer :: ierr,npw_disk,nspinor_disk,nband_disk,band 
 integer :: ipw,my_bcount,npwso,npw_tot,nb_block,base
 integer :: npw_read,nspinor_read,nband_read
 character(len=500) :: msg
!arrays
 real(dp),ABI_CONTIGUOUS pointer :: tmp_eigk(:),tmp_occk(:)
#ifdef HAVE_MPI_IO
 integer :: mpierr,bufsz,gkk_type,cgblock_type
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad
 integer :: sizes(2),subsizes(2),starts(2),types(2)
#endif
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(etsf_main) :: main_folder
 type(etsf_basisdata) :: wave_folder
 type(etsf_electrons) :: electrons_folder
 type(etsf_io_low_error) :: Error
#endif

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfk%rw_mode==WFK_READMODE, "Wfk must be in READMODE")

 if (Wfk%Hdr%headform<40) then 
   write(msg,'(a,i0)')"Too old headform : ",Wfk%Hdr%headform
   MSG_ERROR(msg)
 end if
 !
 ! Look before you leap.
 npw_disk     = Wfk%Hdr%npwarr(ik_ibz)
 nspinor_disk = Wfk%nspinor 
 nband_disk   = Wfk%nband(ik_ibz,spin)
 nb_block     = (band_block(2) - band_block(1) + 1)
 ABI_CHECK(nb_block>0,"nband <=0")
 npw_tot      = npw_disk * nspinor_disk * nb_block

 if (PRESENT(kg_k)) then
   ABI_CHECK(SIZE(kg_k,DIM=2) >= npw_disk,"kg_k too small")
 end if

 if (PRESENT(cg_k)) then
   ABI_CHECK(SIZE(cg_k, DIM=2) >= npw_tot,"cg_k too small")
 end if

 if (PRESENT(eig_k)) then
   if (Wfk%formeig==0) then 
      ABI_CHECK(SIZE(eig_k) >= nband_disk, "GS eig_k too small")
   else if (Wfk%formeig==1) then
      ABI_CHECK(SIZE(eig_k) >= 2*nband_disk**2, "DFPT eig_k too small")
   else
     MSG_ERROR("formeig != [0,1]")
   end if
 end if

 if (PRESENT(occ_k)) then
   ABI_CHECK(SIZE(occ_k) >= nband_disk, "GS eig_k too small")
   !if (Wfk%formeig==0) then
   !  MSG_WARNING("occ_k with formeig != 0")
   !end if
 end if

 SELECT CASE (Wfk%iomode)

 CASE (IO_MODE_FORTRAN) 

   ! Rewind the file to have the correct (k,s) block (if needed)
   call wfk_seek(Wfk,ik_ibz,spin)
   !
   ! Read the first record: npw, nspinor, nband_disk
   read(Wfk%fh, iostat=ierr) npw_read, nspinor_read, nband_read

   if (ierr/=0)then
     msg = "Reading the (npw,nspinor,nband) record of file: "//TRIM(Wfk%fname)
     MSG_ERROR(msg)
   end if

   if ( ANY( (/npw_read, nspinor_read, nband_read/) /= (/npw_disk, nspinor_disk, nband_disk/) )) then
     write(msg,"(a,6(i0,2x))")"Mismatch between (npw, nspinor, nband) read from WFK and those found in HDR ",&
&      npw_read, nspinor_read, nband_read, npw_disk, nspinor_disk, nband_disk
     MSG_ERROR(msg)
   end if
   !
   ! The second record: (k+G) vectors
   if (PRESENT(kg_k)) then
     read(Wfk%fh, iostat=ierr) kg_k(1:3,1:npw_disk)
   else
     read(Wfk%fh, iostat=ierr) ! kg_k(1:3,1:npw_disk)
   end if

   if (ierr/=0)then
     MSG_ERROR("Reading the kg_k record of file: "//TRIM(Wfk%fname))
   end if
   !
   SELECT CASE (Wfk%formeig)
   CASE (0)
     !
     ! The third record: eigenvalues and occupation factors.
     ! write(unitwf) (eigen(iband),iband=1,nband_disk),(occ(iband),iband=1,nband_disk)
     if (PRESENT(eig_k) .or. PRESENT(occ_k)) then

       ABI_MALLOC(tmp_eigk, (nband_disk))
       ABI_MALLOC(tmp_occk, (nband_disk))

       read(Wfk%fh, iostat=ierr) tmp_eigk, tmp_occk

       if (PRESENT(eig_k)) eig_k = tmp_eigk
       if (PRESENT(occ_k)) occ_k = tmp_occk

       ABI_FREE(tmp_eigk)
       ABI_FREE(tmp_occk)

     else
       read(Wfk%fh, iostat=ierr) ! eig_k(1:nband_disk), occ_k(1:nband_k)
     end if
    
     if (ierr/=0)then
       msg = " Reading the GS eigenvalue record of file: "//TRIM(Wfk%fname)
       MSG_ERROR(msg)
     end if
     !
     ! The wave-functions.
     if (PRESENT(cg_k)) then
       !
       npwso = npw_disk*nspinor_disk
       my_bcount = 0
       do band=1,nband_disk
         if (band >= band_block(1) .and. band <= band_block(2)) then
           ipw = my_bcount * npwso
           my_bcount = my_bcount + 1
           read(Wfk%fh, iostat=ierr) cg_k(1:2,ipw+1:ipw+npwso)
         else
           read(Wfk%fh, iostat=ierr) ! cg_k(1:2,ipw+1:ipw+npwso)
         end if
         if (ierr/=0) EXIT
       end do
       !
     else
       !
       do band=1,nband_disk
         read(Wfk%fh, iostat=ierr) ! cg_k(1:2,ipw+1:ipw+npwso)
         if (ierr/=0) EXIT
       end do
       !
     end if
                                                                        
     if (ierr/=0)then
       msg = "Reading the cg record of file: "//TRIM(Wfk%fname)
       MSG_ERROR(msg)
     end if

   CASE (1)
     !
     npwso = npw_disk*nspinor_disk
     my_bcount = 0
     do band=1,nband_disk
       !
       if (PRESENT(eig_k)) then ! Read column matrix of size (2*nband_k**2)
         base = 2*(band-1)*nband_disk
         read(Wfk%fh, iostat=ierr) eig_k(base+1:base+2*nband_disk)
       else
         read(Wfk%fh, iostat=ierr) ! eig_k(2*nband_disk)
       end if

       if (ierr/=0)then
         msg = "Reading the DFPT eigenvalue record of file: "//TRIM(Wfk%fname)
         MSG_ERROR(msg)
       end if

       if (PRESENT(cg_k) .and. (band >= band_block(1) .and. band <= band_block(2)) ) then
         ipw = my_bcount * npwso
         my_bcount = my_bcount + 1
         read(Wfk%fh, iostat=ierr) cg_k(1:2,ipw+1:ipw+npwso)
       else
         read(Wfk%fh, iostat=ierr) ! cg_k(1:2,ipw+1:ipw+npwso)
       end if

       if (ierr/=0)then
         msg = "Reading the DFPT cg_k record of file: "//TRIM(Wfk%fname)
         MSG_ERROR(msg)
       end if
       !
     end do

   CASE DEFAULT
     MSG_ERROR("formeig != [0,1]")
   END SELECT
   ! 
   ! Reached the end of the (k,s) block. Update f90_fptr
   if (ik_ibz < Wfk%nkpt) then
     Wfk%f90_fptr = (/ik_ibz+1,spin,REC_NPW/)
   else
     if (spin==Wfk%nsppol) then
       Wfk%f90_fptr = FPTR_EOF ! EOF condition 
     else
       Wfk%f90_fptr = (/1,spin+1,REC_NPW/)
     end if
   end if

 CASE (IO_MODE_MPI) 

#ifdef HAVE_MPI_IO
   if (PRESENT(kg_k)) then
     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_KG) + xmpio_bsize_frm

     call mpio_read_kg_k(Wfk%fh,my_offset,npw_disk,sc_mode,kg_k,mpierr)
     ABI_CHECK_MPI(mpierr,"reading kg")
   end if
   !
   ! formeig=0 =>  Read both eig and occ in tmp_eigk
   ! formeig=1 =>  Read (nband_k,nband_k) matrix of complex numbers.
   !
   SELECT CASE (Wfk%formeig)

   CASE (0)

     if (PRESENT(eig_k) .or. PRESENT(occ_k)) then
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + xmpio_bsize_frm
       
       call mpio_read_eigocc_k(Wfk%fh,my_offset,nband_disk,Wfk%formeig,sc_mode,tmp_eigk,mpierr)
       ABI_CHECK_MPI(mpierr,"reading eigocc")

       if (PRESENT(eig_k)) eig_k(1:nband_disk) = tmp_eigk(1:nband_disk)
       if (PRESENT(occ_k)) occ_k(1:nband_disk) = tmp_eigk(nband_disk+1:)

       ABI_FREE(tmp_eigk)
     end if

     if (PRESENT(cg_k)) then
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG)
       sizes     = (/npw_disk*nspinor_disk, nband_disk/)
       subsizes  = (/npw_disk*nspinor_disk, band_block(2)-band_block(1)+1/)
       bufsz     = 2 * npw_disk * nspinor_disk * nb_block 
       starts    = (/1, band_block(1)/)

       call mpiotk_read_fsuba_dp2D(Wfk%fh,my_offset,sizes,subsizes,starts,bufsz,cg_k,Wfk%chunk_bsize,sc_mode,Wfk%comm,ierr) 
       ABI_CHECK(ierr==0,"Fortran record too big")
     end if

   CASE (1)

     if (PRESENT(eig_k)) then
       types = (/MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX/)
       sizes = (/nband_disk,npw_disk*nspinor_disk/)
                                                                                                 
       call xmpio_create_fstripes(nband_disk,sizes,types,gkk_type,my_offpad,mpierr)
       ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")
                                                                                                 
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + my_offpad 
                                                                                                 
       call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,gkk_type,'native',MPI_INFO_NULL,mpierr)
       ABI_CHECK_MPI(mpierr,"")
                                                                                                 
       call MPI_TYPE_FREE(gkk_type,mpierr)
       ABI_CHECK_MPI(mpierr,"")
                                                                                                 
       bufsz = (nband_disk**2)

       if (sc_mode==xmpio_collective) then
         call MPI_FILE_READ_ALL(Wfk%fh,eig_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else if (sc_mode==xmpio_single) then
         call MPI_FILE_READ(Wfk%fh,eig_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else 
         MSG_ERROR("Wrong sc_mode")
       end if

       ABI_CHECK_MPI(mpierr,"FILE_READ")
     end if

     if (PRESENT(cg_k)) then
       ABI_CHECK(band_block(1)==1,"band_block(1) !=1 not coded")

       types = (/MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX/)
       sizes = (/npw_disk*nspinor_disk, nband_disk/)
                                                                                                 
       call xmpio_create_fstripes(nb_block,sizes,types,cgblock_type,my_offpad,mpierr)
       ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")
                                                                                                 
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG) + my_offpad 
                                                                                                 
       call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,cgblock_type,'native',MPI_INFO_NULL,mpierr)
       ABI_CHECK_MPI(mpierr,"SET_VIEW")
                                                                                                 
       call MPI_TYPE_FREE(cgblock_type,mpierr)
       ABI_CHECK_MPI(mpierr,"")
                                                                                                 
       bufsz = npw_disk * nspinor_disk * nb_block 

       if (sc_mode==xmpio_collective) then
         call MPI_FILE_READ_ALL(Wfk%fh,cg_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else if (sc_mode==xmpio_single) then
         call MPI_FILE_READ(Wfk%fh,cg_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else 
         MSG_ERROR("Wrong sc_mode")
       end if
       ABI_CHECK_MPI(mpierr,"FILE_READ")
     end if

   CASE DEFAULT
     MSG_ERROR("formeig != [0,1]")
   END SELECT
#endif

 CASE (IO_MODE_ETSF)

#ifdef HAVE_TRIO_ETSF_IO
   if (PRESENT(kg_k)) then
     ! Read the reduced_coordinates_of_plane_waves for this k point.
     wave_folder%red_coord_pw__kpoint_access               = ik_ibz
     wave_folder%red_coord_pw__number_of_coefficients      = npw_disk
     wave_folder%reduced_coordinates_of_plane_waves%data2D => kg_k(:,1:npw_disk)

     call etsf_io_basisdata_get(Wfk%fh, wave_folder, lstat, error)
     ETSF_CHECK_ERROR(lstat,error)
   end if

   ! Read eigenvalues and occupations.
   if (Wfk%formeig==0) then
     if (PRESENT(eig_k) .or. PRESENT(occ_k)) then
       !
       if (PRESENT(eig_k)) then
         electrons_folder%eigenvalues__kpoint_access = ik_ibz
         electrons_folder%eigenvalues__spin_access   = spin
         !electrons_folder%eigenvalues__number_of_states = mband
         electrons_folder%eigenvalues%data1D => eig_k(1:nband_disk)
       end if

       if (PRESENT(occ_k)) then
         electrons_folder%occupations__kpoint_access = ik_ibz
         electrons_folder%occupations__spin_access = spin
         !electrons_folder%occupations__number_of_states = mband
         electrons_folder%occupations%data1D => occ_k(1:nband_disk)
       end if

       call etsf_io_electrons_get(Wfk%fh, electrons_folder, lstat, error)
       ETSF_CHECK_ERROR(lstat,error)
     end if
   else
     MSG_ERROR("formeig !=0 not compatible with ETSF-IO")
   end if

   if (PRESENT(cg_k)) then
     ABI_CHECK(band_block(1)==1,"band_block(1) != 1 not coded")
     ! Read the coefficients_of_wavefunctions
     main_folder%wfs_coeff__kpoint_access             = ik_ibz
     main_folder%wfs_coeff__spin_access               = spin
     main_folder%wfs_coeff__number_of_states          = nb_block ! Number of bands to read
     main_folder%wfs_coeff__number_of_coefficients    = npw_disk
     main_folder%coefficients_of_wavefunctions%data2D => cg_k

     ! With g95, the association done above sometime leads to segfaults.
     ! So we allocate a temporary array to store the wfs of our kpt.
     !ABI_MALLOC(main_folder%coefficients_of_wavefunctions%data2D,(2,npw*nspinor*nband))
     ! Now we copy our values and deallocate the temporary array.
     ! cg(:,icg+1:icg+npw*nspinor*nband)=main_folder%coefficients_of_wavefunctions%data2D
     ! this is better than the previous instruction to optimize virtual memory.
     !do iband=1,nband
     !  ipw=(iband-1)*npwso
     !  cg(:,icg+ipw+1:icg+ipw+npwso)=main_folder%coefficients_of_wavefunctions%data2D(:,ipw+1:ipw+npwso)
     !end do
     !ABI_FREE(main_folder%coefficients_of_wavefunctions%data2D)
     call etsf_io_main_get(Wfk%fh, main_folder, lstat, error)
     ETSF_CHECK_ERROR(lstat,error)
  end if
#endif

 CASE DEFAULT
   write(msg,'(a,i0)')' Wrong value of iomode = ',Wfk%iomode
   MSG_ERROR(msg)
 END SELECT

 DBG_EXIT("COLL")

end subroutine wfk_read_band_block
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_write_band_block
!! NAME
!!  wfk_write_band_block
!!
!! FUNCTION
!!  Write a block of contigous bands.
!!
!! INPUTS
!!  Wfk<type(wfk_t)>=
!!  band_block(2)=Initial and final band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  [kg_k=(:,:)] = G-vectors
!!  [cg_k(:,:)]  = Fourier coefficients
!!  [eig_k(:)] = Eigenvectors
!!  [occ_k(:)] = Eigenvectors
!!
!! PARENTS
!!      dfpt_write_cg,m_wfk,m_wfs
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_write_band_block(Wfk,band_block,ik_ibz,spin,sc_mode,kg_k,cg_k,eig_k,occ_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_write_band_block'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,sc_mode !,mband,rdcg,rdeig,npw_k,nband_k
 type(wfk_t),intent(inout) :: Wfk
!arrays
 integer,intent(in) :: band_block(2)
 integer,intent(in),optional,target :: kg_k(:,:)  !(3,npw_k)
 real(dp),intent(in),optional,target :: cg_k(:,:) ! cg_k(2,rdcg*cgsize2) !(2,npw_k*nspinor*nband)
 real(dp),intent(in),optional,target :: eig_k((2*Wfk%mband)**Wfk%formeig*Wfk%mband)
 real(dp),intent(in),optional,target :: occ_k(Wfk%mband)  

!Local variables-------------------------------
!scalars
 integer :: ierr,npw_disk,nspinor_disk,nband_disk,band 
 integer :: ipw,my_bcount,npwso,npw_tot,nb_block,base
 character(len=500) :: msg
!arrays
 real(dp),ABI_CONTIGUOUS pointer :: tmp_eigk(:) 
#ifdef HAVE_MPI_IO
 integer :: mpierr,bufsz,recnpw_type,gkk_type,cgblock_type
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad 
 integer :: sizes(2),subsizes(2),starts(2),dims(3),types(2) !,statux(MPI_STATUS_SIZE)
 integer(XMPI_OFFSET_KIND) :: bsize_rec(1)
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecords(:)
#endif
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(etsf_main) :: main_folder
 type(etsf_basisdata) :: wave_folder
 type(etsf_electrons) :: electrons_folder
 type(etsf_io_low_error) :: Error
#endif

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfk%rw_mode==WFK_WRITEMODE, "Wfk must be in WRITEMODE")

 if (Wfk%Hdr%headform<40) then 
   write(msg,'(a,i0)')"Too old headform : ",Wfk%Hdr%headform
   MSG_ERROR(msg)
 end if
 !
 ! Look before you leap.
 npw_disk     = Wfk%Hdr%npwarr(ik_ibz)
 nspinor_disk = Wfk%nspinor 
 nband_disk   = Wfk%nband(ik_ibz,spin)
 nb_block     = (band_block(2) - band_block(1) + 1)
 npw_tot      = npw_disk * nspinor_disk * nb_block

 if (PRESENT(kg_k)) then
   ABI_CHECK(SIZE(kg_k,DIM=2) >= npw_disk,"kg_k too small")
 end if

 if (PRESENT(cg_k)) then
   ABI_CHECK(SIZE(cg_k, DIM=2) >= npw_tot,"cg_k too small")
 end if

 if (PRESENT(eig_k)) then
   if (Wfk%formeig==0) then 
      ABI_CHECK(SIZE(eig_k) >= nband_disk, "GS eig_k too small")
      ABI_CHECK(PRESENT(occ_k),"both eig_k and occ_k must be present")
   else if (Wfk%formeig==1) then
      ABI_CHECK(SIZE(eig_k) >= 2*nband_disk**2, "DFPT eig_k too small")
   else
     MSG_ERROR("formeig != [0,1]")
   end if
 end if

 if (PRESENT(occ_k)) then
   ABI_CHECK(SIZE(occ_k) >= nband_disk, "GS eig_k too small")
   ABI_CHECK(PRESENT(eig_k),"both eig_k and occ_k must be present")
   !if (Wfk%formeig==0) then
   !  MSG_WARNING("occ_k with formeig==0")
   !end if
 end if

 SELECT CASE (Wfk%iomode)

 CASE (IO_MODE_FORTRAN) 

   ! Rewind the file to have the correct (k,s) block (if needed)
   call wfk_seek(Wfk,ik_ibz,spin)
   !
   ! Write the first record: npw, nspinor, nband_disk
   write(Wfk%fh, iostat=ierr) npw_disk, nspinor_disk, nband_disk

   if (ierr/=0)then
     msg = "Writing the (npw,nspinor,nband) record of the wfk file:"//TRIM(Wfk%fname)
     MSG_ERROR(msg)
   end if
   !
   ! The second record: (k+G) vectors
   if (PRESENT(kg_k)) then
     write(Wfk%fh, iostat=ierr) kg_k(1:3,1:npw_disk)
   else
     read(Wfk%fh, iostat=ierr) ! kg_k(1:3,1:npw_disk)
   end if

   if (ierr/=0)then
     msg = "Reading the k+g record of the WFK file: "//TRIM(Wfk%fname)
     MSG_ERROR(msg)
   end if
   !
   ! The third record: eigenvalues and occupation factors.
   SELECT CASE (Wfk%formeig)
   CASE (0)
     !write(unitwf) (eigen(iband),iband=1,nband_disk),(occ(iband),iband=1,nband_disk)

     if (PRESENT(eig_k) .and. PRESENT(occ_k)) then
       write(Wfk%fh, iostat=ierr) eig_k, occ_k
     else
       MSG_ERROR("Not coded")
       read(Wfk%fh, iostat=ierr) ! eig_k(1:nband_disk), occ_k(1:nband_k)
     end if
    
     if (ierr/=0)then
       msg = "Writing the GS eigenvalue record of the WFK file: "//TRIM(Wfk%fname)
       MSG_ERROR(msg)
     end if
     !
     ! The wave-functions.
     if (PRESENT(cg_k)) then
       !
       npwso = npw_disk*nspinor_disk
       my_bcount = 0
       do band=1,nband_disk
         if (band >= band_block(1) .and. band <= band_block(2)) then
           ipw = my_bcount * npwso
           my_bcount = my_bcount + 1
           write(Wfk%fh, iostat=ierr) cg_k(1:2,ipw+1:ipw+npwso)
         else
           MSG_ERROR("Not coded")
           read(Wfk%fh, iostat=ierr) ! cg_k(1:2,ipw+1:ipw+npwso)
         end if
         if (ierr/=0) EXIT
       end do
       !
     else
       !
       MSG_ERROR("Not coded")
       do band=1,nband_disk
         read(Wfk%fh, iostat=ierr) ! cg_k(1:2,ipw+1:ipw+npwso)
         if (ierr/=0) EXIT
       end do
       !
     end if
                                                                        
     if (ierr/=0)then
       msg = "Reading the cg record of the WFK file: "//TRIM(Wfk%fname)
       MSG_ERROR(msg)
     end if

   CASE (1)

     ! Write matrix of size (2*nband_k**2)
     ! The wave-functions.
     npwso = npw_disk*nspinor_disk
     my_bcount = 0

     do band=1,nband_disk
       base = 2*(band-1)*nband_disk
       write(Wfk%fh, iostat=ierr) eig_k(base+1:base+2*nband_disk)
       if (ierr/=0)then
         msg = "Writing the DFPT eigenvalue record of the WFK file: "//TRIM(Wfk%fname)
         MSG_ERROR(msg)
       end if
       if (band >= band_block(1) .and. band <= band_block(2)) then
         ipw = my_bcount * npwso
         my_bcount = my_bcount + 1
         write(Wfk%fh, iostat=ierr) cg_k(1:2,ipw+1:ipw+npwso)
       else
         MSG_ERROR("Not coded")
         read(Wfk%fh, iostat=ierr) ! cg_k(1:2,ipw+1:ipw+npwso)
       end if

       if (ierr/=0)then
         msg = "Reading the cg record of the WFK file: "//TRIM(Wfk%fname)
         MSG_ERROR(msg)
       end if
     end do

   CASE DEFAULT
     MSG_ERROR("formeig != [0,1]")
   END SELECT
   ! 
   ! Reached the end of the (k,s) block. Update f90_fptr
   if (ik_ibz < Wfk%nkpt) then
     Wfk%f90_fptr = (/ik_ibz+1,spin,REC_NPW/)
   else
     if (spin==Wfk%nsppol) then
       Wfk%f90_fptr = FPTR_EOF ! EOF condition 
     else
       Wfk%f90_fptr = (/1,spin+1,REC_NPW/)
     end if
   end if

 CASE (IO_MODE_MPI) 

#ifdef HAVE_MPI_IO
   my_offset = Wfk%offset_ks(ik_ibz,spin,REC_NPW) 

   bsize_rec(1) = 3 * xmpi_bsize_int
   call xmpio_write_frmarkers(Wfk%fh,my_offset,sc_mode,1,bsize_rec,ierr)
   ABI_CHECK(ierr==0,"ierr!=0")

   my_offset = Wfk%offset_ks(ik_ibz,spin,REC_NPW) + xmpio_bsize_frm

   call MPI_TYPE_CONTIGUOUS(3, MPI_INTEGER, recnpw_type, mpierr)
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   call MPI_TYPE_COMMIT(recnpw_type,mpierr)
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,recnpw_type,'native',MPI_INFO_NULL,mpierr)
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   call MPI_TYPE_FREE(recnpw_type,mpierr)
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   dims = (/npw_disk, nspinor_disk, nband_disk/)

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(Wfk%fh,dims,SIZE(dims),MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE(Wfk%fh,dims,SIZE(dims),MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
   else 
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_CHECK_MPI(mpierr,"writing REC_NPW")

   if (PRESENT(kg_k)) then
     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_KG) 

     bsize_rec(1) = 3 * npw_disk * xmpi_bsize_int
     call xmpio_write_frmarkers(Wfk%fh,my_offset,sc_mode,1,bsize_rec,ierr)

     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_KG) + xmpio_bsize_frm

     call mpio_write_kg_k(Wfk%fh,my_offset,npw_disk,sc_mode,kg_k,mpierr)
     ABI_CHECK_MPI(mpierr,"mpio_write_kg_k")
   end if

   if (Wfk%formeig==0) then

     if (PRESENT(eig_k) .and. PRESENT(occ_k)) then

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG)

       bsize_rec(1) = 2 * nband_disk * xmpi_bsize_dp
       call xmpio_write_frmarkers(Wfk%fh,my_offset,sc_mode,1,bsize_rec,ierr)

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + xmpio_bsize_frm
       !
       ! Write both eig and occ in tmp_eigk
       bufsz = 2*nband_disk 
       ABI_MALLOC(tmp_eigk, (bufsz))

       tmp_eigk(1:nband_disk)  = eig_k(1:nband_disk)
       tmp_eigk(nband_disk+1:) = occ_k(1:nband_disk)

       call mpio_write_eigocc_k(Wfk%fh,my_offset,nband_disk,Wfk%formeig,sc_mode,tmp_eigk,mpierr)
       ABI_CHECK_MPI(mpierr,"mpio_write_eigocc_k")

       ABI_FREE(tmp_eigk)
     end if

     if (PRESENT(cg_k)) then
       ABI_MALLOC(bsize_frecords, (nb_block))
       bsize_frecords = 2 * npw_disk * nspinor_disk * xmpi_bsize_dp 
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG) + (band_block(1)-1) * (bsize_frecords(1) + 2*xmpio_bsize_frm)
       call xmpio_write_frmarkers(Wfk%fh,my_offset,sc_mode,nb_block,bsize_frecords,ierr)
       ABI_CHECK(ierr==0,"ierr!=0")
       ABI_FREE(bsize_frecords)

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG)
       sizes    = (/npw_disk*nspinor_disk, nband_disk/)
       subsizes = (/npw_disk*nspinor_disk, band_block(2)-band_block(1)+1/)
       bufsz = 2 * npw_disk * nspinor_disk * nb_block 
       starts = (/1, band_block(1)/)

       call mpiotk_write_fsuba_dp2D(Wfk%fh,my_offset,sizes,subsizes,starts,bufsz,cg_k,Wfk%chunk_bsize,sc_mode,Wfk%comm,ierr) 
       ABI_CHECK(ierr==0,"ierr!=0")
     end if

   else if (Wfk%formeig==1) then

     if (PRESENT(eig_k)) then
       types = (/MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX/)
       sizes = (/nband_disk,npw_disk*nspinor_disk/)

       call xmpio_create_fstripes(nband_disk,sizes,types,gkk_type,my_offpad,mpierr)
       ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + my_offpad 

       call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,gkk_type,'native',MPI_INFO_NULL,mpierr)
       ABI_CHECK_MPI(mpierr,"SET_VIEW")

       call MPI_TYPE_FREE(gkk_type,mpierr)
       ABI_CHECK_MPI(mpierr,"")

       bufsz = (nband_disk**2)

       if (sc_mode==xmpio_collective) then
         call MPI_FILE_WRITE_ALL(Wfk%fh,eig_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else if (sc_mode==xmpio_single) then
         call MPI_FILE_WRITE(Wfk%fh,eig_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else 
         MSG_ERROR("Wrong sc_mode")
       end if

       ABI_CHECK_MPI(mpierr,"FILE_WRITE")
     end if

     if (PRESENT(cg_k)) then
       ABI_CHECK(band_block(1)==1,"band_block(1) !=1 not coded")

       types = (/MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX/)
       sizes = (/npw_disk*nspinor_disk, nband_disk/)
                                                                                                 
       call xmpio_create_fstripes(nb_block,sizes,types,cgblock_type,my_offpad,mpierr)
       ABI_CHECK_MPI(mpierr,"xmpio_create_fstripes")
                                                                                                 
       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG) + my_offpad 
                                                                                                 
       call MPI_FILE_SET_VIEW(Wfk%fh,my_offset,MPI_BYTE,cgblock_type,'native',MPI_INFO_NULL,mpierr)
       ABI_CHECK_MPI(mpierr,"SET_VIEW")
                                                                                                 
       call MPI_TYPE_FREE(cgblock_type,mpierr)
       ABI_CHECK_MPI(mpierr,"")
                                                                                                 
       bufsz = npw_disk * nspinor_disk * nb_block 
       if (sc_mode==xmpio_collective) then
         call MPI_FILE_WRITE_ALL(Wfk%fh,cg_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else if (sc_mode==xmpio_single) then
         call MPI_FILE_WRITE(Wfk%fh,cg_k,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
       else 
         MSG_ERROR("Wrong sc_mode")
       end if
       ABI_CHECK_MPI(mpierr,"FILE_WRITE")
     end if

   else 
     MSG_ERROR("formeig not in [0,1]")
   end if
#endif

 CASE (IO_MODE_ETSF)
#ifdef HAVE_TRIO_ETSF_IO
   if (PRESENT(kg_k)) then
     ! Write the reduced_coordinates_of_plane_waves for this k point.
     wave_folder%red_coord_pw__kpoint_access               = ik_ibz
     wave_folder%red_coord_pw__number_of_coefficients      = npw_disk
     wave_folder%reduced_coordinates_of_plane_waves%data2D => kg_k(:,1:npw_disk)

     call etsf_io_basisdata_put(Wfk%fh,wave_folder,lstat,error)
     ETSF_CHECK_ERROR(lstat,error)
   end if

   ! Write eigenvalues and occupation factors.
   if (Wfk%formeig==0) then
     !
     if (PRESENT(eig_k) .or. PRESENT(occ_k)) then
       if (PRESENT(eig_k)) then
         electrons_folder%eigenvalues__kpoint_access    = ik_ibz
         electrons_folder%eigenvalues__spin_access      = spin
         electrons_folder%eigenvalues__number_of_states = Wfk%mband
         electrons_folder%eigenvalues%data1D            => eig_k
       end if

       if (PRESENT(occ_k)) then
         electrons_folder%occupations__kpoint_access    = ik_ibz
         electrons_folder%occupations__spin_access      = spin
         electrons_folder%occupations__number_of_states = Wfk%mband
         electrons_folder%occupations%data1D            => occ_k
       end if

       call etsf_io_electrons_put(Wfk%fh, electrons_folder, lstat, error)
       ETSF_CHECK_ERROR(lstat,error)
    end if

   else if (Wfk%formeig==1) then
     MSG_ERROR("formeig==1 not compatible with ETSF-IO")

   else
     MSG_ERROR("formeig != [0,1]")
   end if

   if (PRESENT(cg_k)) then
     ! Write the wavefunctions.
     ABI_CHECK(band_block(1)==1,"band_block(1) != 1 not coded")
     main_folder%wfs_coeff__kpoint_access             = ik_ibz
     main_folder%wfs_coeff__spin_access               = spin
     main_folder%wfs_coeff__number_of_states          = nb_block ! Number of bands to write
     main_folder%wfs_coeff__number_of_coefficients    = npw_disk
     main_folder%coefficients_of_wavefunctions%data2D => cg_k
     !
     ! With g95, the association done above sometime leads to segfaults.
     ! So we allocate a temporary array to store the wfs of our kpt.
     !ABI_MALLOC(main_folder%coefficients_of_wavefunctions%data2D,(2,npw*nspinor*nband))
     !The following instruction leads to segmentation faults sometimes
     !I changed it for a longer expression !Tonatiuh Rangel 2009
     !main_folder%coefficients_of_wavefunctions%data2D = cg(:, icg + 1:icg + npw * nspinor * nband)
     !jj=0
     !do ii=icg+1,icg+npw*nspinor*nband
     !  jj=jj+1
     !  main_folder%coefficients_of_wavefunctions%data2D(:,jj)=cg(:,ii)
     !end do
     call etsf_io_main_put(Wfk%fh, main_folder, lstat, error)
     ETSF_CHECK_ERROR(lstat,error)
     !ABI_FREE(main_folder%coefficients_of_wavefunctions%data2D)
  end if
#endif

 CASE DEFAULT
   write(msg,'(a,i0)')' Wrong value of iomode = ',Wfk%iomode
   MSG_ERROR(msg)
 END SELECT

 DBG_EXIT("COLL")

end subroutine wfk_write_band_block
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_bmask
!! NAME
!!  wfk_read_bmask
!!
!! FUNCTION
!!  Read a set of bands specified by a mask
!!
!! INPUTS
!!  Wfk<type(wfk_t)>=
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  [kg_k=(:,:)] = G-vectors
!!  [cg_k(:,:)]  = Fourier coefficients
!!  [eig_k(:)] = Eigenvectors
!!  [occ_k(:)] = Occupation 
!!
!! PARENTS
!!      m_wfk,m_wfs
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_read_bmask(Wfk,bmask,ik_ibz,spin,sc_mode,kg_k,cg_k,eig_k,occ_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_read_bmask'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,sc_mode
 type(wfk_t),intent(inout) :: Wfk
!arrays
 logical,intent(in) :: bmask(Wfk%mband)
 integer,intent(out), DEV_CONTARRD optional,target :: kg_k(:,:)  !(3,npw_k)
 real(dp),intent(out), DEV_CONTARRD optional,target :: cg_k(:,:) !(2,npw_k*nspinor*nband)
 real(dp),intent(out),optional,target :: eig_k((2*Wfk%mband)**Wfk%formeig*Wfk%mband)
 real(dp),intent(out),optional,target :: occ_k(Wfk%mband)  

!Local variables-------------------------------
!scalars
 integer :: ierr,npw_disk,nspinor_disk,nband_disk,ipw,my_bcount,cnt,npwso,npw_tot,pt1,pt2,band
 integer :: npw_read,nspinor_read,nband_read,nb_tot,ncount,my_bcnt,my_maxb,base
 character(len=500) :: msg
!arrays 
 real(dp),ABI_CONTIGUOUS pointer :: tmp_eigk(:),tmp_occk(:)
#ifdef HAVE_MPI_IO
 integer :: mpierr,cgscatter_type,cg_type,method,block,nblocks,nbxblock
 integer :: bstart,bstop,bufsz,ugsz,brest,max_nband
 integer(XMPI_OFFSET_KIND) :: my_offset,base_ofs 
 integer :: band_block(2),sizes(2),subsizes(2),starts(2) !,statux(MPI_STATUS_SIZE)
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)
 real(dp),allocatable :: buffer(:,:)
#endif

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfk%rw_mode==WFK_READMODE, "Wfk must be in READMODE")

 if (Wfk%Hdr%headform<40) then 
   write(msg,'(a,i0)')"Too old headform : ",Wfk%Hdr%headform
   MSG_ERROR(msg)
 end if
 !
 ! Look before you leap.
 npw_disk = Wfk%Hdr%npwarr(ik_ibz)
 nspinor_disk = Wfk%nspinor
 nband_disk = Wfk%nband(ik_ibz,spin)
 nb_tot = COUNT(bmask)
 npw_tot = npw_disk * nspinor_disk * nb_tot

 if (PRESENT(kg_k)) then
   ABI_CHECK((SIZE(kg_k,DIM=2) >= npw_disk),"kg_k too small")
 end if

 if (PRESENT(cg_k)) then
  ABI_CHECK(SIZE(cg_k, DIM=2) >= npw_tot, "Too small cg_k")
 end if

 if (PRESENT(eig_k)) then
   if (Wfk%formeig==0) then 
      ABI_CHECK(SIZE(eig_k) >= nband_disk, "GS eig_k too small")
   else if (Wfk%formeig==1) then
      ABI_CHECK(SIZE(eig_k) >= 2*nband_disk**2, "DFPT eig_k too small")
   else
     MSG_ERROR("formeig != [0,1]")
   end if
 end if

 if (PRESENT(occ_k)) then
   ABI_CHECK(Wfk%formeig==0,"occ_k with formeig != 0")
   ABI_CHECK(SIZE(occ_k) >= nband_disk, "GS eig_k too small")
 end if

 SELECT CASE (Wfk%iomode)

 CASE (IO_MODE_FORTRAN) 

   ! Rewind the file to have the correct (k,s) block (if needed)
   call wfk_seek(Wfk,ik_ibz,spin)
   !
   ! Read the first record: npw, nspinor, nband_disk
   read(Wfk%fh, iostat=ierr) npw_read, nspinor_read, nband_read

   if (ierr/=0)then
     msg = " Reading the (npw,nspinor,nband) record of the WFK file: "//TRIM(Wfk%fname)
     MSG_ERROR(msg)
   end if

   if ( ANY( (/npw_read, nspinor_read, nband_read/) /= (/npw_disk, nspinor_disk, nband_disk/) )) then
     write(msg,"(a,6(i0,2x))")"Mismatch between (npw, nspinor, nband) read from WFK and those found in HDR ",&
&      npw_read, nspinor_read, nband_read, npw_disk, nspinor_disk, nband_disk
     MSG_ERROR(msg)
   end if
   !
   ! The second record: (k+G) vectors
   if (PRESENT(kg_k)) then
     read(Wfk%fh, iostat=ierr) kg_k(1:3,1:npw_disk)
   else
     read(Wfk%fh, iostat=ierr) ! kg_k(1:3,1:npw_disk)
   end if

   if (ierr/=0)then
     msg = "Reading the k+g record of the WFK file: "//TRIM(Wfk%fname)
     MSG_ERROR(msg)
   end if
   !
   ! The third record: eigenvalues and occupation factors.
   if (Wfk%formeig==0) then

     if (PRESENT(eig_k) .or. PRESENT(occ_k)) then
       ABI_MALLOC(tmp_eigk, (nband_disk))
       ABI_MALLOC(tmp_occk, (nband_disk))
                                                    
       read(Wfk%fh, iostat=ierr) tmp_eigk, tmp_occk
                                                    
       if (PRESENT(eig_k)) eig_k = tmp_eigk
       if (PRESENT(occ_k)) occ_k = tmp_occk
                                                    
       ABI_FREE(tmp_eigk)
       ABI_FREE(tmp_occk)

     else
       read(Wfk%fh, iostat=ierr) ! eig_k(1:nband_disk)
     end if
    
     if (ierr/=0)then
       msg = "Reading the eigenvalue record of the WFK file: "//TRIM(Wfk%fname)
       MSG_ERROR(msg)
     end if
     !
     ! The wave-functions.
     if (PRESENT(cg_k)) then
       npwso = npw_disk*nspinor_disk
       my_bcount = 0
       do band=1,nband_disk
         if (bmask(band)) then
           ipw = my_bcount * npwso
           my_bcount = my_bcount + 1
           read(Wfk%fh, iostat=ierr) cg_k(1:2,ipw+1:ipw+npwso)
         else
           read(Wfk%fh, iostat=ierr) ! cg_k(1:2,ipw+1:ipw+npwso)
         end if
         if (ierr/=0) EXIT
       end do
       !
     else
       do band=1,nband_disk
         read(Wfk%fh, iostat=ierr) ! cg_k(1:2,ipw+1:ipw+npwso)
         if (ierr/=0) EXIT
       end do
     end if
                                                                        
     if (ierr/=0)then
       msg = "Reading the cg record of the WFK file: "//TRIM(Wfk%fname)
       MSG_ERROR(msg)
     end if

   else if (Wfk%formeig==1) then
     ! Read matrix of size (2*nband_k**2)
     npwso = npw_disk*nspinor_disk
     my_bcount = 0

     do band=1,nband_disk
       !
       base = 2*(band-1)*nband_disk
       if (PRESENT(eig_k)) then
         read(Wfk%fh, iostat=ierr) eig_k(base+1:base+2*nband_disk)
       else
         read(Wfk%fh, iostat=ierr) ! eig_k(base+1:base+2*nband_disk)
       end if
       if (ierr/=0)then
         msg = "Reading the DFPT eigenvalue record of the WFK file: "//TRIM(Wfk%fname)
         MSG_ERROR(msg)
       end if

       if (bmask(band).and.PRESENT(cg_k)) then
         ipw = my_bcount * npwso
         my_bcount = my_bcount + 1
         read(Wfk%fh, iostat=ierr) cg_k(1:2,ipw+1:ipw+npwso)
       else
         read(Wfk%fh, iostat=ierr) ! cg_k(1:2,ipw+1:ipw+npwso)
       end if

       if (ierr/=0)then
         msg = "Reading the cg record of the WFK file: "//TRIM(Wfk%fname)
         MSG_ERROR(msg)
       end if
     end do

   else 
     MSG_ERROR("formeig != [0,1]")
   end if
   ! 
   ! Reached the end of the (k,s) block. Update f90_fptr
   if (ik_ibz < Wfk%nkpt) then
     Wfk%f90_fptr = (/ik_ibz+1,spin,REC_NPW/)
   else
     if (spin==Wfk%nsppol) then
       Wfk%f90_fptr = FPTR_EOF ! EOF condition
     else
       Wfk%f90_fptr = (/1,spin+1,REC_NPW/)
     end if
   end if

 CASE (IO_MODE_MPI) 
#ifdef HAVE_MPI_IO

   if (PRESENT(kg_k)) then
     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_KG) + xmpio_bsize_frm

     call mpio_read_kg_k(Wfk%fh,my_offset,npw_disk,sc_mode,kg_k,mpierr)
     ABI_CHECK_MPI(mpierr,"mpio_read_kg_k")
   end if

   ! The third record: eigenvalues and occupation factors.
   if (PRESENT(eig_k) .or. PRESENT(occ_k)) then
     my_offset = Wfk%offset_ks(ik_ibz,spin,REC_EIG) + xmpio_bsize_frm
     !
     ! formeig=0 =>  Read both eig and occ in tmp_eigk.
     ! formeig=1 =>  Read (nband_k,nband_k) matrix of complex numbers.
     !
     call mpio_read_eigocc_k(Wfk%fh,my_offset,nband_disk,Wfk%formeig,sc_mode,tmp_eigk,mpierr)
     ABI_CHECK_MPI(mpierr,"mpio_read_eigocc")

     if (Wfk%formeig==0) then
       if (PRESENT(eig_k)) eig_k(1:nband_disk) = tmp_eigk(1:nband_disk)
       if (PRESENT(occ_k)) occ_k(1:nband_disk) = tmp_eigk(nband_disk+1:)
     else if (Wfk%formeig==1) then
       if (PRESENT(eig_k)) eig_k(1:2*nband_disk**2) = tmp_eigk(1:2*nband_disk**2)
     else 
       MSG_ERROR("formeig not in [0,1]")
     end if
                                                                        
     ABI_FREE(tmp_eigk)
   end if

   if (PRESENT(cg_k)) then
     method = 0

     SELECT CASE (method)
     CASE (0)
       ! DATA SIEVING: 
       !   read max_nband states in chuncks of nbxblock, then extract my states according to bmask.
       !
       ! MAX number of bands read by the procs in the communicator 
       my_maxb = nband_disk
       do band=nband_disk,1,-1
         if (bmask(band)) then
           my_maxb = band
           EXIT
         end if
       end do
       call xmpi_max(my_maxb,max_nband,Wfk%comm,mpierr)
       !max_nband = nband_disk
       !
       ! MPI-IO crashes if we try to read a large number of bands in a single call.
       nbxblock = max_nband
       if ((2*npw_disk*nspinor_disk*nbxblock*xmpi_bsize_dp) > Wfk%chunk_bsize) then
         nbxblock = Wfk%chunk_bsize / (2*npw_disk*nspinor_disk*xmpi_bsize_dp)
         if (nbxblock == 0) nbxblock = 50
       end if
       !nbxblock = 2

       nblocks = max_nband / nbxblock
       brest   = MOD(max_nband, nbxblock)
       if (brest /= 0) nblocks = nblocks + 1

       !write(std_out,*)"in buffered bmask"
       !write(std_out,*)"nb_tot",nb_tot,"nblocks",nblocks,"nbxblock",nbxblock

       base_ofs = Wfk%offset_ks(ik_ibz,spin,REC_CG)
       sizes = (/npw_disk*nspinor_disk, nband_disk/)

       my_bcnt = 0  ! index of my band in cg_k
       do block=1,nblocks
         bstart = 1 + (block-1) * nbxblock
         bstop  = bstart + nbxblock - 1 
         if (bstop > max_nband) bstop = max_nband
         !
         ! Allocate and read the buffer
         band_block = (/bstart, bstop/)
         ugsz = npw_disk*nspinor_disk
         bufsz = 2*ugsz*(bstop-bstart+1)
         ABI_MALLOC(buffer, (2,bufsz))

         !write(std_out,*)"  bstart,bstop, ",band_block
         subsizes = (/npw_disk*nspinor_disk, band_block(2)-band_block(1)+1/)
         starts = (/1, bstart/)

         call mpiotk_read_fsuba_dp2D(Wfk%fh,base_ofs,sizes,subsizes,starts,&
&          bufsz,buffer,Wfk%chunk_bsize,sc_mode,Wfk%comm,ierr) 
         ABI_CHECK(ierr==0,"Fortran record too big")

         ! Extract my bands from buffer.
         do band=bstart,bstop
           if (bmask(band)) then
             my_bcnt = my_bcnt + 1
             pt1 = 1 + (my_bcnt - 1) * ugsz
             pt2 = 1 + (band - bstart) * ugsz
             cg_k(:,pt1:pt1+ugsz-1) = buffer(:,pt2:pt2+ugsz-1)
           end if
         end do

         ABI_FREE(buffer)
       end do

     CASE (1,2)
       call MPI_TYPE_CONTIGUOUS(npw_disk*nspinor_disk,MPI_DOUBLE_COMPLEX,cg_type,mpierr)
       ABI_CHECK_MPI(mpierr,"type_contigous")

       if (method==1) then
         ncount = nb_tot
         ABI_MALLOC(block_length,(ncount+2))
         ABI_MALLOC(block_type, (ncount+2))
         ABI_MALLOC(block_displ,(ncount+2))

         block_length(1)=1
         block_displ (1)=0
         block_type  (1)=MPI_LB

         my_bcount = 1
         do band=1,Wfk%mband
           if (bmask(band)) then
             my_bcount = my_bcount + 1
             block_length(my_bcount) = 1
             block_type(my_bcount) = cg_type
             block_displ(my_bcount) = xmpio_bsize_frm + &
&              (band-1) * (2*npw_disk*nspinor_disk*xmpi_bsize_dp + 2*xmpio_bsize_frm)
           end if
         end do

         block_length(ncount+2) = 1
         block_displ (ncount+2) = block_displ(my_bcount)
         block_type  (ncount+2) = MPI_UB

       else if (method==2) then
         ! this file view is not efficient but it's similar to the 
         ! one used in wff_readwrite. Let's see if MPI-IO likes it!
         ncount = nb_tot* nspinor_disk * npw_disk

         ABI_MALLOC(block_length,(ncount+2))
         ABI_MALLOC(block_type, (ncount+2))
         ABI_MALLOC(block_displ,(ncount+2))
                                             
         block_length(1)=1
         block_displ (1)=0
         block_type  (1)=MPI_LB
         !
         ! The view starts at REC_CG
         cnt = 1
         do band=1,Wfk%mband
           if (bmask(band)) then
             base_ofs =  xmpio_bsize_frm + &
&              (band-1) * (2*npw_disk*nspinor_disk*xmpi_bsize_dp + 2*xmpio_bsize_frm)
             do ipw=1,npw_disk*nspinor_disk
               cnt = cnt + 1
               block_length(cnt) = 1
               block_type(cnt)   = MPI_DOUBLE_COMPLEX
               block_displ(cnt)  = base_ofs + 2*(ipw-1)*xmpi_bsize_dp 
             end do
           end if
         end do

         block_length(ncount+2) = 1
         block_displ (ncount+2) = block_displ(cnt)
         block_type  (ncount+2) = MPI_UB
       end if

       call xmpio_type_struct(ncount+2,block_length,block_displ,block_type,cgscatter_type,mpierr)
       ABI_CHECK_MPI(mpierr,"type_struct")

       ABI_FREE(block_length)
       ABI_FREE(block_type)
       ABI_FREE(block_displ)

       call MPI_TYPE_FREE(cg_type, mpierr)
       ABI_CHECK_MPI(mpierr,"MPI_TYPE_FREE")

       my_offset = Wfk%offset_ks(ik_ibz,spin,REC_CG)

       call MPI_FILE_SET_VIEW(Wfk%fh, my_offset, MPI_BYTE, cgscatter_type, 'native', MPI_INFO_NULL, mpierr)
       ABI_CHECK_MPI(mpierr,"SET_VIEW")
                                                                                                                 
       call MPI_TYPE_FREE(cgscatter_type, mpierr)
       ABI_CHECK_MPI(mpierr,"MPI_TYPE_FREE")

       call MPI_FILE_READ_ALL(Wfk%fh, cg_k, npw_tot, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
       ABI_CHECK_MPI(mpierr,"FILE_READ_ALL")

     CASE DEFAULT
       MSG_ERROR("Wrong method")
     END SELECT
   end if
#endif

 CASE (IO_MODE_ETSF)
   MSG_ERROR("Not Coded")

 CASE DEFAULT
   write(msg,'(a,i0)')' Wrong value of iomode = ',Wfk%iomode
   MSG_ERROR(msg)
 END SELECT

 DBG_EXIT("COLL")

end subroutine wfk_read_bmask
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_eigk
!! NAME
!!  wfk_read_eigk
!!
!! FUNCTION
!!  Helper function to read all the eigenvalues for a given (k-point,spin)
!!
!! INPUTS
!!  Wfk<type(wfk_t)>= WFK file handler
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=spin index
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  eig_k(1:nband_k) = GS Eigenvalues for the given (k,s)
!!  occ_k(1:nband_k) = Occupation factors for the given (k,s)
!!
!! PARENTS
!!      m_wfk,mrggkk
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_read_eigk(Wfk,ik_ibz,spin,sc_mode,eig_k,occ_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_read_eigk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,sc_mode
 type(wfk_t),intent(inout) :: Wfk
!arrays
 real(dp),intent(out) :: eig_k((2*Wfk%mband)**Wfk%formeig*Wfk%mband)
 real(dp),optional,intent(out) :: occ_k(Wfk%mband)

!Local variables-------------------------------
!scalars
 integer,parameter :: band_block00(2)=(/0,0/)

!************************************************************************

 if (PRESENT(occ_k)) then
   ABI_CHECK(Wfk%formeig==0,"formeig !=0")
   call wfk_read_band_block(Wfk,band_block00,ik_ibz,spin,sc_mode,eig_k=eig_k,occ_k=occ_k)
 else
   call wfk_read_band_block(Wfk,band_block00,ik_ibz,spin,sc_mode,eig_k=eig_k)
 end if

end subroutine wfk_read_eigk
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_read_eigenvalues
!! NAME
!!  wfk_read_eigenvalues
!!
!! FUNCTION
!!  Read all the GS eigenvalues stored in the WFK file fname.
!!
!! INPUTS
!!  fname=Name of the file
!!  comm=MPI communicator.
!!
!! OUTPUTS
!!  eigen = In input: nullified pointer
!!          In output: eigen(mband,nkpt,nsppol) contains the GS eigevalues.
!!  Hdr_out<hdr_type>=The header of the file
!!
!! PARENTS
!!      loper3,setup_bse,setup_bse_interp,setup_screening,setup_sigma
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_read_eigenvalues(fname,eigen,Hdr_out,comm,occ)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_read_eigenvalues'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=*),intent(in) :: fname
 type(hdr_type),intent(out) :: Hdr_out
!arrays
 real(dp),pointer :: eigen(:,:,:)
 real(dp),pointer,optional :: occ(:,:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,formeig0=0
 integer :: ik_ibz,spin,my_rank,nprocs,ierr,iomode,funt,sc_mode,mband
 type(wfk_t) :: Wfk

!************************************************************************

 my_rank = xcomm_rank(comm)
 nprocs  = xcomm_size(comm)

 iomode = iomode_from_fname(fname) 

 if (my_rank==master) then
   ! Open the file and read the GS eigenvalues.
   sc_mode = xmpio_single
   funt = get_unit()
   call wfk_open_read(Wfk,fname,formeig0,iomode,funt,xmpi_self,Hdr_out=Hdr_out)

   ABI_MALLOC(eigen, (Wfk%mband,Wfk%nkpt,Wfk%nsppol))
   eigen = HUGE(zero)

   if (PRESENT(occ)) then 
     ABI_MALLOC(occ, (Wfk%mband,Wfk%nkpt,Wfk%nsppol))
     occ = HUGE(zero)

     do spin=1,Wfk%nsppol
       do ik_ibz=1,Wfk%nkpt
         call wfk_read_eigk(Wfk,ik_ibz,spin,sc_mode,eigen(:,ik_ibz,spin),occ_k=occ(:,ik_ibz,spin))
       end do
     end do

   else
     do spin=1,Wfk%nsppol
       do ik_ibz=1,Wfk%nkpt
         call wfk_read_eigk(Wfk,ik_ibz,spin,sc_mode,eigen(:,ik_ibz,spin))
       end do
     end do
   end if
   !
   ! Close the file.
   call wfk_close(Wfk)
 end if

 ! Broadcast data
 if (nprocs>1) then
   call hdr_comm(Hdr_out,master,my_rank,comm)
   mband = MAXVAL(Hdr_out%nband)
   if (my_rank/=master) then
     ABI_MALLOC(eigen, (mband,Hdr_out%nkpt,Hdr_out%nsppol))
     eigen = HUGE(zero)
     if (PRESENT(occ)) then 
       ABI_MALLOC(occ, (mband,Hdr_out%nkpt,Hdr_out%nsppol))
       occ = HUGE(zero)
     end if
   end if
   call xmpi_bcast(eigen,master,comm,ierr)
   if (PRESENT(occ)) then
     call xmpi_bcast(occ,master,comm,ierr)
   end if
 end if

end subroutine wfk_read_eigenvalues
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_write_allgkk
!! NAME
!!  wfk_write_allgkk
!!
!! FUNCTION
!!  Write all GKK matrix elements in the WFK file fname.
!!
!! INPUTS
!!
!! OUTPUTS
!!
!! PARENTS
!!      dfpt_write_cg
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_write_allgkk(Wfk,sc_mode,eigen)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_write_allgkk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sc_mode
 type(wfk_t),intent(inout) :: Wfk
!arrays
 real(dp),intent(in) :: eigen(2*Wfk%mband**2*Wfk%nkpt*Wfk%nsppol)

!Local variables-------------------------------
!scalars
 integer :: spin,ik_ibz,nband_k,ptr
!arrays
 integer,parameter :: band_block00(2)=(/0,0/)

!************************************************************************

 ptr=1
 do spin=1,Wfk%nsppol
   do ik_ibz=1,Wfk%nkpt
     nband_k = Wfk%nband(ik_ibz,spin)
     call wfk_write_band_block(Wfk,band_block00,ik_ibz,spin,sc_mode,eig_k=eigen(ptr:))
     ptr = ptr + 2*nband_k**2
   end do
 end do

end subroutine wfk_write_allgkk
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_seek
!! NAME
!!  wfk_seek
!!
!! FUNCTION
!!   Move the internal file pointer so that it points to the
!!   block (ik_ibz, spin). Needed only if iomode==IO_MODE_FORTRAN 
!!
!! INPUTS
!!   ik_ibz,spin = (k-point,spin) indices
!! 
!! SIDE EFFECTS
!!   Wfk<type(Wfk_t)> : modifies Wfk%f90_fptr and the internal F90 file pointer.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_seek(Wfk,ik_ibz,spin)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_seek'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: ik_ibz,spin
 type(Wfk_t),intent(inout) :: Wfk

!Local variables-------------------------------
 integer :: ierr,ik_fpt,spin_fpt,recn_wanted,recn_fpt,rec_type
 character(len=500) :: msg

! *************************************************************************

 SELECT CASE (Wfk%iomode)
 CASE (IO_MODE_FORTRAN)
   !
   ! Find the position inside the file.
   !
   if (ALL(Wfk%f90_fptr==FPTR_EOF)) then ! handle the EOF condition
     if (Wfk%debug) then
       call wrtout(std_out,"EOF condition","PERS")
     end if
     recn_fpt = Wfk%recn_eof
   else
     ik_fpt   = Wfk%f90_fptr(1)
     spin_fpt = Wfk%f90_fptr(2)
     rec_type = Wfk%f90_fptr(3)
     recn_fpt = Wfk%recn_ks(ik_fpt,spin_fpt, rec_type)
   end if

   recn_wanted = Wfk%recn_ks(ik_ibz,spin, REC_NPW)

   if (Wfk%debug) then
     write(msg,'(a,3(i0,2x))')"seeking ik_ibz, spin, recn_wanted-recn_fpt: ",ik_ibz,spin,recn_wanted - recn_fpt
     call wrtout(std_out,msg,"PERS")
   end if

   call mvrecord(Wfk%fh, (recn_wanted - recn_fpt) ,ierr)
   ABI_CHECK(ierr==0,"error in mvrecord")

   Wfk%f90_fptr = (/ik_ibz, spin, REC_NPW/)

 CASE DEFAULT
   msg = ABI_FUNC//" should not be called when Wfk%iomode /= IO_MODE_FORTRAN" 
   MSG_ERROR(msg)
 END SELECT 

end subroutine wfk_seek
!!***

!----------------------------------------------------------------------

!!****f* m_wfkfile/wfk_compute_offsets
!! NAME
!!  wfk_compute_offsets
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_compute_offsets(Wfk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_compute_offsets'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wfk_t),intent(inout) :: Wfk

!Local variables-------------------------------
!scalars
 integer :: spin,ik_ibz,npw_k,nband_k,bsize_frm,mpi_type_frm,base !,band
 integer(XMPI_OFFSET_KIND) :: offset

! *************************************************************************

 SELECT CASE (Wfk%iomode)

 CASE (IO_MODE_FORTRAN) 
   !
   ! Compute record number for Fortran IO
   ABI_MALLOC(Wfk%recn_ks,(Wfk%nkpt,Wfk%nsppol,REC_NUM))

   ! TODO this should point to the end of the Header!
   ! if we want to have meaningful absolute recn.
   base = 0 
   do spin=1,Wfk%nsppol
     do ik_ibz=1,Wfk%nkpt
       nband_k = Wfk%nband(ik_ibz,spin)
       Wfk%recn_ks(ik_ibz,spin, REC_NPW) = base + 1
       Wfk%recn_ks(ik_ibz,spin, REC_KG)  = base + 2
       Wfk%recn_ks(ik_ibz,spin, REC_EIG) = base + 3
       Wfk%recn_ks(ik_ibz,spin, REC_CG)  = base + 4
       base = Wfk%recn_ks(ik_ibz,spin,REC_CG)
       if (Wfk%formeig==0) then
         base = base + (nband_k-1)
       else if (Wfk%formeig==1) then
         base = base + 2*(nband_k-1)
       else
         MSG_ERROR("formeig != [0,1]")
       end if
     end do
   end do
   !
   ! Save EOF position
   Wfk%recn_eof = base + 1

 CASE (IO_MODE_MPI) 
   !
   ! Compute offsets for MPI-IO.
   ABI_MALLOC(Wfk%offset_ks,(Wfk%nkpt,Wfk%nsppol,REC_NUM))

   bsize_frm    = xmpio_bsize_frm    ! Byte length of the Fortran record marker.
   mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

   ! fform must be 0 ??
   ! The offset of the Header. TODO
   offset = Wfk%hdr_offset ! hdr_offset(Hdr)

   do spin=1,Wfk%nsppol
     do ik_ibz=1,Wfk%nkpt
       npw_k   = Wfk%Hdr%npwarr(ik_ibz)
       nband_k = Wfk%nband(ik_ibz,spin)
       !
       !---------------------------------------------------------------------------
       ! First record: npw, nspinor, nband_disk
       !---------------------------------------------------------------------------
       Wfk%offset_ks(ik_ibz,spin,REC_NPW) = offset

       if (Wfk%Hdr%headform>=40) then 
         ! npw, nspinor, nband_disk
         offset = offset +  3*xmpi_bsize_int + 2*bsize_frm
       else 
         MSG_ERROR("Old headforms < 40 are not supported")
       end if
       Wfk%offset_ks(ik_ibz,spin,REC_KG) = offset
       !
       !---------------------------------------------------------------------------
       ! Second record: (k+G) vectors
       ! kg_k(1:3,1:npw_k)
       !---------------------------------------------------------------------------
       offset = offset + 3*npw_k*xmpi_bsize_int + 2*bsize_frm
       Wfk%offset_ks(ik_ibz,spin,REC_EIG) = offset
       !
       !---------------------------------------------------------------------------
       ! Third record: eigenvalues
       !---------------------------------------------------------------------------
       if (Wfk%formeig==0) then
         ! eigen(1:nband_k), occ(1:nband_k)
         offset = offset + 2*nband_k*xmpi_bsize_dp + 2*bsize_frm 
         Wfk%offset_ks(ik_ibz,spin,REC_CG) = offset
         !
         ! Wavefunction coefficients 
         ! do band=1,nband_k; write(unitwf) cg_k(1:2,npw_k*nspinor); end do
         offset = offset + nband_k * (2*npw_k*Wfk%nspinor*xmpi_bsize_dp + 2*bsize_frm)

       else if (Wfk%formeig==1) then
         ! read(unitwf) eigen(2*nband_k)
         Wfk%offset_ks(ik_ibz,spin,REC_CG) = offset + (2*nband_k*xmpi_bsize_dp + 2*bsize_frm)

         offset = offset + & 
&          nband_k * (2*npw_k*Wfk%nspinor*xmpi_bsize_dp + 2*bsize_frm) + &
&          nband_k * (2*nband_k*xmpi_bsize_dp + 2*bsize_frm)

       else 
         MSG_ERROR("Wrong formeig")
       end if
       
     end do ! ik_ibz
   end do ! spin
   !
   ! Save EOF offset
   Wfk%offset_eof = offset

   ! Check for possible wraparound errors.
   if (ANY(Wfk%offset_ks <= 0) .or. Wfk%offset_eof < 0) then
     MSG_ERROR("Found negative offset. File too large for MPI-IO!!!")
   end if

 END SELECT

 if (Wfk%debug) then
   call wfk_show_offsets(Wfk)
 end if

end subroutine wfk_compute_offsets
!!***

!----------------------------------------------------------------------

!!****f* m_wfkfile/wfk_show_offsets
!! NAME
!!  wfk_show_offsets
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_show_offsets(Wfk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_show_offsets'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wfk_t),intent(inout) :: Wfk

!Local variables-------------------------------
!scalars
 integer :: spin,ik_ibz

! *************************************************************************

 SELECT CASE (Wfk%iomode)

 CASE (IO_MODE_FORTRAN) 
   ! TODO this should point to the end of the Header!
   ! if we want to have meaningful absolute recn.
   do spin=1,Wfk%nsppol
     do ik_ibz=1,Wfk%nkpt
       write(std_out,"(a,2(i0,2x),a,4(a,i0,a))")                   &
&        "(ik_ibz, spin) ",ik_ibz,spin,ch10,                       &
&        "  recn(REC_NPW): ",Wfk%recn_ks(ik_ibz,spin,REC_NPW),ch10,&
&        "  recn(REC_KG) : ",Wfk%recn_ks(ik_ibz,spin,REC_KG), ch10,&
&        "  recn(REC_EIG): ",Wfk%recn_ks(ik_ibz,spin,REC_EIG),ch10,& 
&        "  recn(REC_CG) : ",Wfk%recn_ks(ik_ibz,spin,REC_CG),ch10
     end do
   end do
   !
   ! Write EOF position
   write(std_out,"(a,i0)")"recn_eof ",Wfk%recn_eof

 CASE (IO_MODE_MPI) 
   write(std_out,"(a,i0)")"hdr_offset ",Wfk%hdr_offset

   do spin=1,Wfk%nsppol
     do ik_ibz=1,Wfk%nkpt
       write(std_out,"(a,2(i0,2x),a,4(a,i0,a))")                       &
&        "(ik_ibz, spin) ",ik_ibz,spin,ch10,                           &
&        "  offset(REC_NPW): ",Wfk%offset_ks(ik_ibz,spin,REC_NPW),ch10,&
&        "  offset(REC_KG) : ",Wfk%offset_ks(ik_ibz,spin,REC_KG), ch10,&
&        "  offset(REC_EIG): ",Wfk%offset_ks(ik_ibz,spin,REC_EIG),ch10,& 
&        "  offset(REC_CG) : ",Wfk%offset_ks(ik_ibz,spin,REC_CG),ch10
     end do ! ik_ibz
   end do ! spin
   !
   ! Write EOF position
   write(std_out,"(a,i0)")"offset_eof ",Wfk%offset_eof
 END SELECT

end subroutine wfk_show_offsets
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/mpio_read_kg_k
!! NAME
!!  mpio_read_kg_k
!!
!! FUNCTION
!!
!! INPUTS
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  kg_k=(3,npw_disk) = G-vectors
!!  mpierr=MPI error status (error check is delegated to the caller)
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine mpio_read_kg_k(fh,offset,npw_disk,sc_mode,kg_k,mpierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mpio_read_kg_k'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,npw_disk,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: mpierr
!arrays
 integer,intent(out) :: kg_k(3,npw_disk)

!Local variables-------------------------------
!scalars
 integer :: kg_k_type,ncount,myfh
 integer(XMPI_OFFSET_KIND) :: my_offset

!************************************************************************

 ! Workarounds for XLF
 myfh      = fh
 ncount    = 3*npw_disk 
 my_offset = offset

 call MPI_TYPE_CONTIGUOUS(ncount, MPI_INTEGER, kg_k_type, mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_COMMIT(kg_k_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,kg_k_type,'native',MPI_INFO_NULL,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 if (sc_mode==xmpio_collective) then
   call MPI_FILE_READ_ALL(myfh,kg_k,ncount,MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
 else if (sc_mode==xmpio_single) then
   !call MPI_File_seek(myfh, 0, MPI_SEEK_SET,mpierr)
   call MPI_FILE_READ(myfh,kg_k,ncount,MPI_INTEGER, MPI_STATUS_IGNORE,mpierr)
 else 
   MSG_ERROR("Wrong sc_mode")
 end if

 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_FREE(kg_k_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

end subroutine mpio_read_kg_k
#endif
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/mpio_write_kg_k
!! NAME
!!  mpio_write_kg_k
!!
!! FUNCTION
!!
!! INPUTS
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for writing by current proc.
!!    xmpio_collective ==> for collective write.
!!  kg_k=(3,npw_disk) = G-vectors
!!
!! OUTPUTS
!!  mpierr=MPI error status (error check is delegated to the caller)
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine mpio_write_kg_k(fh,offset,npw_disk,sc_mode,kg_k,mpierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mpio_write_kg_k'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,npw_disk,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: mpierr
!arrays
 integer,intent(in) :: kg_k(3,npw_disk)

!Local variables-------------------------------
!scalars
 integer :: myfh,kg_k_type,ncount
 integer(XMPI_OFFSET_KIND) :: my_offset

!************************************************************************

 DBG_ENTER("COLL")

 ! Workarounds for XLF
 myfh      = fh
 ncount    = 3*npw_disk 
 my_offset = offset

 call MPI_TYPE_CONTIGUOUS(ncount, MPI_INTEGER, kg_k_type, mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_COMMIT(kg_k_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,kg_k_type,'native',MPI_INFO_NULL,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_FREE(kg_k_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 if (sc_mode==xmpio_collective) then
   call MPI_FILE_WRITE_ALL(myfh,kg_k,ncount,MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
 else if (sc_mode==xmpio_single) then
   call MPI_FILE_WRITE(myfh,kg_k,ncount,MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
 else 
   MSG_ERROR("Wrong sc_mode")
 end if

 ABI_HANDLE_MPIERR(mpierr)

 DBG_EXIT("COLL")

end subroutine mpio_write_kg_k
#endif
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/mpio_read_eigocc_k
!! NAME
!!  mpio_read_eigocc_k
!!
!! FUNCTION
!!
!! INPUTS
!!  fh
!!  offset
!!  nband_disk
!!  formeig
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!
!! OUTPUTS
!!  buffer(:)
!!  mpierr=MPI error status.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine mpio_read_eigocc_k(fh,offset,nband_disk,formeig,sc_mode,buffer,mpierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mpio_read_eigocc_k'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,nband_disk,formeig,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: mpierr
!arrays
 real(dp),pointer :: buffer(:) 

!Local variables-------------------------------
!scalars
 integer :: myfh,bufsz,gkk_type,eneocc_type
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad !,fmarker
!arrays
 integer :: sizes(2),subsizes(2),starts(2) 

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 SELECT CASE (formeig)
 CASE (0)
   !
   ! Read both eig and occ in buffer
   bufsz = 2*nband_disk 
   my_offset = offset
   ABI_MALLOC(buffer, (bufsz))

   call MPI_TYPE_CONTIGUOUS(bufsz, MPI_DOUBLE_PRECISION, eneocc_type, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_COMMIT(eneocc_type,mpierr)
   ABI_HANDLE_MPIERR(mpierr)
                                                                                                   
   call MPI_FILE_SET_VIEW(myfh, my_offset, MPI_BYTE, eneocc_type, 'native', MPI_INFO_NULL, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_FREE(eneocc_type,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_READ_ALL(myfh,buffer,bufsz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_READ(myfh,buffer,bufsz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,mpierr)
   else 
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_HANDLE_MPIERR(mpierr)

 CASE (1)
   ! Read the (nband_k,nband_k) matrix with the (complex) GKK matrix elements.
   bufsz    = (nband_disk**2)
   sizes    = (/nband_disk, nband_disk/)
   subsizes = (/nband_disk, nband_disk/)
   starts   = (/1, 1/)

   ABI_MALLOC(buffer, (2*bufsz))

   !my_offset = offset - xmpio_bsize_frm
   !call xmpio_read_dp(myfh,my_offset,sc_mode,2*nband_disk,buffer,fmarker,mpierr)
   !write(std_out,*)buffer(1:2*nband_disk)
   !MSG_ERROR("Done")

   call xmpio_create_fsubarray_2D(sizes,subsizes,starts,MPI_DOUBLE_COMPLEX,gkk_type,my_offpad,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   ! TODO: Rationalize the offsets
   my_offset = offset + my_offpad - xmpio_bsize_frm

   call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,gkk_type,'native',MPI_INFO_NULL,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_FREE(gkk_type, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_READ_ALL(myfh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_READ(myfh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
   else
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_HANDLE_MPIERR(mpierr)

 CASE DEFAULT
   MSG_ERROR("formeig not in [0,1]")
 END SELECT

end subroutine mpio_read_eigocc_k
#endif
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/mpio_write_eigocc_k
!! NAME
!!  mpio_write_eigocc_k
!!
!! FUNCTION
!!
!! INPUTS
!!  fh
!!  offset
!!  nband_disk
!!  formeig
!!  sc_mode= MPI-IO option
!!    xmpio_single     ==> for writing  by current proc.
!!    xmpio_collective ==> for collective write.
!!
!! OUTPUTS
!!  buffer(:)
!!  mpierr=MPI error status.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine mpio_write_eigocc_k(fh,offset,nband_disk,formeig,sc_mode,buffer,mpierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mpio_write_eigocc_k'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,nband_disk,formeig,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: mpierr
!arrays
 real(dp),intent(in) :: buffer(:) 

!Local variables-------------------------------
!scalars
 integer :: bufsz,gkk_type,eneocc_type,myfh
 integer(XMPI_OFFSET_KIND) :: my_offset,my_offpad 
!arrays
 integer :: sizes(2),subsizes(2),starts(2) 

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 SELECT CASE (formeig)
 CASE (0)
   !
   ! write both eig and occ in buffer
   my_offset = offset 

   bufsz = 2*nband_disk 
   ABI_CHECK(SIZE(buffer) >= bufsz, "buffer too small")

   call MPI_TYPE_CONTIGUOUS(bufsz, MPI_DOUBLE_PRECISION, eneocc_type, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_COMMIT(eneocc_type,mpierr)
   ABI_HANDLE_MPIERR(mpierr)
                                                                                                   
   call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,eneocc_type,'native',MPI_INFO_NULL,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_FREE(eneocc_type,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(myfh,buffer,bufsz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE(myfh,buffer,bufsz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,mpierr)
   else 
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_HANDLE_MPIERR(mpierr)

 CASE (1)
   !MSG_ERROR("formeig ==1 with MPI-IO not tested")
   ! write the (nband_k,nband_k) matrix with the (complex) GKK matrix elements.
   bufsz    = (nband_disk**2)
   sizes    = (/nband_disk, nband_disk/)
   subsizes = (/nband_disk, nband_disk/)
   starts   = (/1, 1/)

   ABI_CHECK(SIZE(buffer) >= bufsz, "buffer too small")

   call xmpio_create_fsubarray_2D(sizes,subsizes,starts,MPI_DOUBLE_COMPLEX,gkk_type,my_offpad,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   ! TODO: Rationalize the offsets
   my_offset = offset + my_offpad - xmpio_bsize_frm

   call MPI_FILE_SET_VIEW(myfh,my_offset,MPI_BYTE,gkk_type,'native',MPI_INFO_NULL,mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   call MPI_TYPE_FREE(gkk_type, mpierr)
   ABI_HANDLE_MPIERR(mpierr)

   if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(myfh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE(myfh,buffer,bufsz,MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,mpierr)
   else
     MSG_ERROR("Wrong sc_mode")
   end if
   ABI_HANDLE_MPIERR(mpierr)

 CASE DEFAULT
   MSG_ERROR("formeig not in [0,1]")
 END SELECT

end subroutine mpio_write_eigocc_k
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_cp
!! NAME
!!  wfk_cp
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_cp(source,dest,formeig)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_cp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig
 character(len=*),intent(in) :: source,dest

!Local variables-------------------------------
!scalars
 integer :: iomode_source,iomode_dest,mband,nspinor
 integer :: ikpt,spin,nband_k,mpw
 character(len=len(source)+len(dest)+500) :: command
#ifdef HAVE_FC_COMMAND_LINE
 character(len=len(source)+len(dest)+500) :: msg
#endif
 type(wfk_t) :: Wfk_source,Wfk_dest
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp),allocatable :: cg_k(:,:),eig_k(:),occ_k(:)

! *************************************************************************

 !command = "mv "//TRIM(source)//" "//TRIM(dest)
 command = "cp "//TRIM(source)//" "//TRIM(dest)

 !call EXECUTE_COMMAND_LINE(command, wait=.TRUE., exitstat=ierr, cmdmsg=msg) 

 !if (ierr/=0) then
 !  MSG_WARNING(TRIM(command)//ch10//TRIM(msg))
 !end if

 iomode_source = iomode_from_fname(source)
 iomode_dest   = iomode_from_fname(dest)

 call wfk_open_read(Wfk_source,source,formeig,iomode_source,get_unit(),xmpi_self)

 call wfk_open_write(Wfk_dest,Wfk_source%Hdr,dest,formeig,iomode_dest,get_unit(),xmpi_self)

 mband   = Wfk_source%mband
 mpw     = MAXVAL(Wfk_source%Hdr%npwarr)
 nspinor = Wfk_source%nspinor

 ABI_MALLOC(kg_k, (3,mpw))
 ABI_MALLOC(cg_k, (2,mpw*nspinor*mband))
 ABI_MALLOC(eig_k, ((2*mband)**formeig*mband) )
 ABI_MALLOC(occ_k, (mband))

 do spin=1,Wfk_source%nsppol
   do ikpt=1,Wfk_source%nkpt
      nband_k = Wfk_source%nband(ikpt,spin)

      call wfk_read_band_block(Wfk_source,(/1,nband_k/),ikpt,spin,xmpio_single,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)

      call wfk_write_band_block(Wfk_dest,(/1,nband_k/),ikpt,spin,xmpio_single, kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
   end do
 end do

 ABI_FREE(kg_k)
 ABI_FREE(cg_k)
 ABI_FREE(eig_k)
 ABI_FREE(occ_k)

 call wfk_close(Wfk_source)
 call wfk_close(Wfk_dest)

end subroutine wfk_cp
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_prof
!! NAME
!!  wfk_prof
!!
!! FUNCTION
!!
!! INPUTS
!! 
!! PARENTS
!!      ioprof
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_prof(wfk_fname,formeig,nband,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_prof'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
 use interfaces_56_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nband,formeig,comm
 character(len=*),intent(in) :: wfk_fname

!Local variables-------------------------------
!scalars
 integer,parameter :: rdwr1=1,master=0,optkg1=1,option1=1,tim_rwwf0=0,icg0=0,headform0=0
 integer :: iomode,wfk_unt,ik_ibz,spin,ierr,ii,option,mband
 integer :: npw_disk,nband_disk,mcg,fform,nband_read,sc_mode,my_rank,nproc
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
 type(hdr_type) :: Hdr
 type(Wfk_t) :: Wfk
 type(wffile_type) :: wff
 type(MPI_type) :: MPI_enreg_seq
!arrays
 !integer,parameter :: io_modes(2) = (/IO_MODE_FORTRAN, IO_MODE_MPI/)
 integer,parameter :: io_modes(1) = (/IO_MODE_MPI/)
 integer :: ngfft(18)
 logical,allocatable :: my_bmask(:)
 integer,allocatable :: kg_k(:,:)
 real(dp),allocatable :: eig_k(:),cg_k(:,:),occ_k(:)

! *************************************************************************

 my_rank = xcomm_rank(comm)
 nproc   = xcomm_size(comm)
 sc_mode = xmpio_collective

 if (wfk_fname==ABI_NOFILE) then
   MSG_ERROR("Not coded!")
   !call wfk_create_wfkfile(new_fname,Hdr,iomode,formeig,Kvars,cwtimes,xmpi_self)
   !call kvars_free(Kvars)
 else
   call wrtout(std_out,ABI_FUNC//": about to read "//TRIM(wfk_fname),"COLL")
 end if

 call hdr_read_from_fname(Hdr,wfk_fname,fform,comm)

 wfk_unt = get_unit()

 do ii=1,SIZE(io_modes)
   iomode = io_modes(ii)
   !do option=1,3
   do option=1,3,2
     write(std_out,*)"iomode, option",iomode,option
     call cwtime(cpu,wall,gflops,"start")

     SELECT CASE (option)

     CASE (1)
       !
       call wfk_open_read(Wfk,wfk_fname,formeig,iomode,wfk_unt,comm)

       do spin=1,Hdr%nsppol
         do ik_ibz=1,Hdr%nkpt
           npw_disk   = Hdr%npwarr(ik_ibz)
           nband_disk = Wfk%nband(ik_ibz,spin)

           nband_read = nband
           if (nband_read <=0) nband_read = nband_disk
           if (nband_read > nband_disk) then
             write(msg,'(a,2(i0,1x))')" nband_read cannot be greater than nband_disk while: ",nband_read,nband_disk
             MSG_ERROR(msg)
           end if

           mcg = npw_disk*Hdr%nspinor*nband_read

           ABI_MALLOC(eig_k,((2*Wfk%mband)**formeig*Wfk%mband))
           ABI_MALLOC(occ_k,(Wfk%mband))

           ABI_MALLOC(kg_k,(3,npw_disk))
           ABI_MALLOC(cg_k,(2,mcg))
           ABI_CHECK_ALLOC("out of memory in cg_k") 

           ! Read the block of bands for this (k,s).
           call wfk_read_band_block(Wfk,(/1,nband_read/),ik_ibz,spin,xmpio_collective,&
&            kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
           !
           ABI_FREE(eig_k)
           ABI_FREE(occ_k)
           ABI_FREE(kg_k)
           ABI_FREE(cg_k)
         end do !ik_ibz
       end do !spin

       call wfk_close(Wfk)

     CASE (2)

       call wfk_open_read(Wfk,wfk_fname,formeig,iomode,wfk_unt,comm)

       do spin=1,Hdr%nsppol
         do ik_ibz=1,Hdr%nkpt
           npw_disk   = Hdr%npwarr(ik_ibz)
           nband_disk = Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)

           nband_read = nband
           if (nband_read <=0) nband_read = nband_disk
           if (nband_read > nband_disk) then
             write(msg,'(a,2(i0,1x))')"nband_read cannot be greater than nband_disk while: ",nband_read,nband_disk
             MSG_ERROR(msg)
           end if

           ABI_MALLOC(my_bmask,(MAXVAL(Hdr%nband)))
           my_bmask=.FALSE.
           my_bmask(1:nband_read) = .TRUE.

           ABI_MALLOC(eig_k,((2*nband_disk)**formeig*nband_disk))
           ABI_MALLOC(kg_k,(3,npw_disk))
           ABI_MALLOC(occ_k,(nband_disk))

           mcg = npw_disk*Hdr%nspinor*COUNT(my_bmask)
           ABI_MALLOC(cg_k,(2,mcg))
           ABI_CHECK_ALLOC("out of memory in cg_k") 

           call wfk_read_bmask(Wfk,my_bmask,ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
           !call wfk_read_band_block(Wfk,(/1,nband_read/),ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)

           ABI_FREE(my_bmask)
           ABI_FREE(eig_k)
           ABI_FREE(occ_k)
           ABI_FREE(kg_k)
           ABI_FREE(cg_k)
         end do !ik_ibz
       end do !spin

       call wfk_close(Wfk)
       !
     CASE (3)
       !Fake MPI_type for the sequential part.
       ngfft(1:6) = (/12,12,12,13,13,13/)
       call initmpi_seq(MPI_enreg_seq)
       call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft(2),ngfft(3),'all')
       call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfft(2),ngfft(3),'all')

       call WffOpen(iomode,comm,wfk_fname,ierr,wff,master,my_rank,wfk_unt) !,spaceComm_mpiio) ! optional argument
       ABI_CHECK(ierr==0,"ierr!=0")

       call hdr_free(Hdr)
       call hdr_io(fform,Hdr,1,wff)
       call WffKg(wff,optkg1)

       do spin=1,Hdr%nsppol
         do ik_ibz=1,Hdr%nkpt

           npw_disk   = Hdr%npwarr(ik_ibz)
           nband_disk = Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)
                                                                                                                               
           nband_read = nband
           if (nband_read <=0) nband_read = nband_disk
           if (nband_read > nband_disk) then
             write(msg,'(a,2(i0,1x))')" nband_read cannot be greater than nband_disk while: ",nband_read,nband_disk
             MSG_ERROR(msg)
           end if

           mband = MAXVAL(Hdr%nband)
           mcg = npw_disk*Hdr%nspinor*nband_read
                                                                                                                               
           ABI_MALLOC(eig_k,((2*mband)**formeig*mband))
           ABI_MALLOC(occ_k,(mband))
                                                                                                                               
           ABI_MALLOC(kg_k,(3,optkg1*npw_disk))
           ABI_MALLOC(cg_k,(2,mcg))
           ABI_CHECK_ALLOC("out of memory in cg_k") 
           !
           ! Read the block of bands for this (k,s).
           call rwwf(cg_k,eig_k,formeig,headform0,icg0,ik_ibz,spin,kg_k,mband,mcg,MPI_enreg_seq,nband_read,&
&           nband_disk,npw_disk,Hdr%nspinor,occ_k,option1,optkg1,tim_rwwf0,Wff)

           ABI_FREE(eig_k)
           ABI_FREE(occ_k)
           ABI_FREE(kg_k)
           ABI_FREE(cg_k)

         end do !ik_ibz
       end do !spin

       call WffClose(wff,ierr)
       call destroy_mpi_enreg(MPI_enreg_seq)

     CASE DEFAULT
       MSG_ERROR("Wrong method")
     END SELECT

     call cwtime(cpu,wall,gflops,"stop")
     write(msg,'(3(a,i2),2(a,f8.2))')&
&      " iomode: ",iomode,", nproc: ",nproc,", option: ",option,", cpu: ",cpu,", wall:",wall
     call wrtout(std_out,msg,"COLL")
   end do
 end do

 call hdr_free(Hdr)

end subroutine wfk_prof
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_create_wfkfile
!! NAME
!!  wfk_create_wfkfile
!!
!! FUNCTION
!!
!! INPUTS
!! 
!! PARENTS
!!      ioprof
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_create_wfkfile(wfk_fname,Hdr,iomode,formeig,Kvars,cwtimes,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_create_wfkfile'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,formeig,comm
 character(len=*),intent(in) :: wfk_fname
!arrays
 real(dp),intent(out) :: cwtimes(2)
 type(hdr_type),intent(in) :: Hdr
 type(kvars_t),target,intent(out) :: Kvars(Hdr%nkpt)

!Local variables-------------------------------
!scalars
 integer :: nkpt,nsppol,nspinor,ierr,sc_mode
 integer :: ik_ibz,spin,funt,nband_k,npw_k,istwfk_k
 real(dp) :: cpu,wall,gflops,ucvol
 type(wfk_t) :: Wfk
!arrays
 integer :: nband(Hdr%nkpt,Hdr%nsppol)
 integer,pointer :: kg_k(:,:)
 real(dp) :: kpoint(3),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: cg_k(:,:),eig_k(:),occ_k(:)

!************************************************************************

 cwtimes = zero

 nband   = RESHAPE(Hdr%nband, (/Hdr%nkpt,Hdr%nsppol/) )
 nkpt    = Hdr%nkpt
 nsppol  = Hdr%nsppol
 nspinor = Hdr%nspinor

 call metric(gmet,gprimd,dev_null,rmet,Hdr%rprimd,ucvol)

 ! Generate the G-vectors from input Hdr%ecut.
 do ik_ibz=1,nkpt
   kpoint   = Hdr%kptns(:,ik_ibz)
   istwfk_k = Hdr%istwfk(ik_ibz) 
   call get_kg(kpoint,istwfk_k,Hdr%ecut,gmet,npw_k,Kvars(ik_ibz)%kg_k)
   ABI_CHECK(npw_k == Hdr%npwarr(ik_ibz),"npw_k != Hdr%npwarr(ik)")
 end do
 !
 ! Open the file for writing.
 sc_mode = xmpio_collective

 call cwtime(cpu,wall,gflops,"start")
 funt = get_unit()

 call wfk_open_write(Wfk,Hdr,wfk_fname,formeig,iomode,funt,comm,write_frm=.TRUE.)

 call cwtime(cpu,wall,gflops,"stop")
 cwtimes = cwtimes + (/cpu,wall/)

 do spin=1,nsppol
   do ik_ibz=1,nkpt

     nband_k = nband(ik_ibz,spin)
     npw_k   = Hdr%npwarr(ik_ibz)
     ABI_MALLOC(cg_k, (2,npw_k*nspinor*nband_k))
     ABI_MALLOC(eig_k, ((2*Wfk%mband)**formeig*Wfk%mband) )
     ABI_MALLOC(occ_k, (Wfk%mband))

     kg_k => Kvars(ik_ibz)%kg_k
     !
     ! Fill cg_k, eig_k, occ_k using a deterministic algorithm so that 
     ! we can check the correctness of the reading.
     call fill_or_check("fill",Hdr,Kvars(ik_ibz),ik_ibz,spin,formeig,kg_k,cg_k,eig_k,occ_k,ierr)
     ABI_CHECK(ierr==0,"filling")
     
     call cwtime(cpu,wall,gflops,"start")

     call wfk_write_band_block(Wfk,(/1,nband_k/),ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)

     call cwtime(cpu,wall,gflops,"stop")
     cwtimes = cwtimes + (/cpu,wall/)

     ABI_FREE(cg_k)
     ABI_FREE(eig_k)
     ABI_FREE(occ_k)
   end do
 end do

 ! Close the file
 call wfk_close(Wfk)

end subroutine wfk_create_wfkfile
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_check_wfkfile
!! NAME
!!  wfk_check_wfkfile
!!
!! FUNCTION
!!
!! INPUTS
!! 
!! PARENTS
!!      ioprof
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_check_wfkfile(wfk_fname,Hdr,iomode,method,formeig,Kvars,cwtimes,comm,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_check_wfkfile'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,formeig,comm,method
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: wfk_fname
!arrays
 real(dp),intent(out) :: cwtimes(2)
 type(hdr_type),intent(in) :: Hdr
 type(kvars_t),intent(in) :: Kvars(Hdr%nkpt)

!Local variables-------------------------------
!scalars
 integer :: nkpt,nsppol,nspinor,ik_ibz,spin,funt,nband_k,npw_k,sc_mode 
 integer :: my_ierr,restart,restartpaw,is,ik,ntests,test,mband
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
 type(wfk_t) :: Wfk
!arrays
 integer :: nband(Hdr%nkpt,Hdr%nsppol),spins(Hdr%nsppol),kindices(Hdr%nkpt)
 integer,allocatable :: kg_k(:,:)
 real(dp),allocatable :: cg_k(:,:),eig_k(:),occ_k(:)
 logical,allocatable :: bmask(:)

!************************************************************************

 !write(msg,"(3a,i2)")"Checking file: ",TRIM(wfk_fname),", with iomode = ",iomode
 !call wrtout(std_out,msg,"COLL")

 ierr = 0
 cwtimes = zero

 nband   = RESHAPE(Hdr%nband, (/Hdr%nkpt,Hdr%nsppol/) )
 nkpt    = Hdr%nkpt
 nsppol  = Hdr%nsppol
 nspinor = Hdr%nspinor

 ! Open the file for writing.
 call cwtime(cpu,wall,gflops,"start")
 funt = get_unit()

 call wfk_open_read(Wfk,wfk_fname,formeig,iomode,funt,comm)
 mband = Wfk%mband

 call cwtime(cpu,wall,gflops,"stop")
 cwtimes = cwtimes + (/cpu,wall/)

 ntests = 2

 do test=1,ntests
   spins    = (/(spin, spin=1,Hdr%nsppol)/)
   kindices = (/(ik_ibz, ik_ibz=1,Hdr%nkpt)/)

   if (test==2) then ! Reverse the indices
     spins    = (/(spin, spin=Hdr%nsppol,1,-1)/)
     kindices = (/(ik_ibz, ik_ibz=Hdr%nkpt,1,-1)/)
   end if
   !
   do is=1,SIZE(spins)
     spin = spins(is)
     do ik=1,SIZE(kindices)
       ik_ibz = kindices(ik)

       if (Wfk%debug) then 
         call hdr_check(Wfk%fform,Wfk%fform,Hdr,Wfk%Hdr,"COLL",restart,restartpaw)
       end if

       nband_k = nband(ik_ibz,spin)
       npw_k   = Hdr%npwarr(ik_ibz)

       ABI_MALLOC(kg_k, (3,npw_k))
       ABI_MALLOC(cg_k, (2,npw_k*nspinor*nband_k))
       ABI_MALLOC(eig_k, ((2*mband)**Wfk%formeig*mband) )
       ABI_MALLOC(occ_k, (mband))
       !
       !sc_mode = xmpio_collective
       sc_mode = xmpio_single
       ABI_MALLOC(bmask, (mband))
       bmask = .FALSE.
       bmask(1:nband_k) = .TRUE.

       call cwtime(cpu,wall,gflops,"start")

       if (method==0) then
         call wfk_read_band_block(Wfk,(/1,nband_k/),ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
       else if (method==1) then
         call wfk_read_bmask(Wfk,bmask,ik_ibz,spin,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
       else 
         MSG_ERROR("Wrong method")
       end if

       !call wfk_read_eigk(Wfk,ik_ibz,spin,sc_mode,eig_k)
       !write(std_out,*)"eig_k",eig_k

       call cwtime(cpu,wall,gflops,"stop")
       cwtimes = cwtimes + (/cpu,wall/)

       ! Check the correctness of the reading.
       call fill_or_check("check",Hdr,Kvars(ik_ibz),ik_ibz,spin,formeig,kg_k,cg_k,eig_k,occ_k,my_ierr)

       if (my_ierr/=0) then
         write(msg,"(a,i0)")"fill_or_check returned my_ierr: ",my_ierr
         ierr = my_ierr
         MSG_WARNING(msg)
       end if

       ABI_FREE(kg_k)
       ABI_FREE(cg_k)
       ABI_FREE(eig_k)
       ABI_FREE(occ_k)

       ABI_FREE(bmask)
     end do
   end do
   !
 end do ! test

 ! Close the file
 call wfk_close(Wfk)

end subroutine wfk_check_wfkfile
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/fill_or_check
!! NAME
!!  fill_or_check
!!
!! FUNCTION
!!
!! INPUTS
!! 
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine fill_or_check(task,Hdr,Kvars,ik_ibz,spin,formeig,kg_k,cg_k,eig_k,occ_k,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fill_or_check'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,formeig
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: task
!arrays
 integer,intent(in) :: kg_k(:,:)
 real(dp),intent(inout) :: cg_k(:,:),eig_k(:),occ_k(:)
 type(hdr_type),intent(in) :: Hdr
 type(kvars_t),intent(in) :: Kvars

!Local variables-------------------------------
!scalars
 integer :: nkpt,nsppol,nspinor,nband_k,npw_k,band,ipw,kspad,ii,base,idx,mpw,eigsz
 character(len=500) :: msg
!arrays
 integer,allocatable :: ref_kg_k(:,:)
 real(dp),allocatable :: ref_eig_k(:),ref_occ_k(:),ref_cg_k(:,:)
!************************************************************************

 ierr = 0

 nkpt    = Hdr%nkpt
 nsppol  = Hdr%nsppol
 nspinor = Hdr%nspinor

 nband_k = Hdr%nband(ik_ibz + (spin-1)*nkpt)
 npw_k   = Hdr%npwarr(ik_ibz)

 ABI_MALLOC(ref_kg_k,(3,npw_k))
 ABI_MALLOC(ref_eig_k,((2*nband_k)**formeig*nband_k))
 ABI_MALLOC(ref_occ_k,(nband_k))
 ABI_MALLOC(ref_cg_k,(2,npw_k*nspinor*nband_k))

 ref_kg_k = Kvars%kg_k

 ! Pad values according to (k,s).
 kspad = (spin-1)*nkpt + (ik_ibz-1) * nband_k

 if (formeig==0) then
   eigsz = nband_k
   do band=1,nband_k
     ref_eig_k(band) = half * (kspad + band)
     ref_occ_k(band) = two  * (kspad + band)
   end do
 else if (formeig==1) then
   eigsz = 2*nband_k**2
   base=0
   do band=1,nband_k
     do ii=1,2*nband_k
       idx = base + ii
       ref_eig_k(idx) = ii*(kspad + band)
     end do
     base = base + 2*nband_k
   end do
 end if

 mpw = npw_k*nspinor*nband_k
 do ipw=1,mpw
   ref_cg_k(1,ipw) =  ipw + kspad
   ref_cg_k(2,ipw) = -ipw + kspad
 end do

 SELECT CASE (task)
 CASE ("fill")
   cg_k(:,1:mpw) = ref_cg_k(:,1:mpw)
   if (formeig==0) then 
     eig_k(1:nband_k) = ref_eig_k
     occ_k(1:nband_k) = ref_occ_k
   else 
     eig_k(1:2*nband_k**2) = ref_eig_k
   end if

 CASE ("check")

   if (ANY( ABS(cg_k(:,1:mpw) - ref_cg_k) > zero)) then
     ierr = ierr + 1
     MSG_WARNING("Difference in cg_k")
   end if

   if (ANY( ABS(kg_k - ref_kg_k) > zero)) then
     ierr = ierr + 2
     MSG_WARNING("Difference in kg_k")
     !write(std_out,*)"ref_kg_k",ref_kg_k
     !write(std_out,*)"kg_k",kg_k
   end if

   if (ANY( ABS(eig_k(1:eigsz) - ref_eig_k) > zero)) then 
     ierr = ierr + 4
     MSG_WARNING("Difference in eig_k")
     !write(std_out,*)"ref_eig_k",ref_eig_k
     !write(std_out,*)"eig_k",eig_k
   end if

   if (formeig==0) then
     if (ANY( ABS(occ_k(1:nband_k) - ref_occ_k) > zero)) then
       ierr = ierr + 8
       MSG_WARNING("occ_k")
       !write(std_out,*)"ref_occ_k",ref_occ_k
       !write(std_out,*)"occ_k",occ_k
     end if
   end if

   write(msg,"(a,3(i0,2x))")" (ik_ibz, spin, ierr) ",ik_ibz,spin,ierr
   if (ierr/=0) then
     MSG_WARNING(TRIM(msg)//": FAILED")
   else
     call wrtout(std_out,TRIM(msg)//": OK","COLL")
   end if

 CASE DEFAULT
   MSG_ERROR("Wrong task")
 END SELECT

 ABI_FREE(ref_kg_k)
 ABI_FREE(ref_eig_k)
 ABI_FREE(ref_occ_k)
 ABI_FREE(ref_cg_k)

end subroutine fill_or_check
!!***

!----------------------------------------------------------------------

!!****f* m_wfk/wfk_diff
!! NAME
!!  wfk_diff
!!
!! FUNCTION
!!  Compare two WFK file for binary equality
!!
!! INPUTS
!! 
!! PARENTS
!!      outwf
!!
!! CHILDREN
!!      hdr_free,hdr_read_from_fname,wfk_close,wfk_open_read
!!      wfk_read_band_block,wrtout
!!
!! SOURCE

subroutine wfk_diff(fname1,fname2,formeig,comm,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfk_diff'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: formeig,comm
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: fname1,fname2

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iomode1,iomode2,ik_ibz,spin,mband,nband_k
 integer :: npw_k,mcg,fform1,fform2,sc_mode,my_rank,nproc
 character(len=500) :: msg
 type(hdr_type) :: Hdr1,Hdr2
 type(Wfk_t) :: Wfk1,Wfk2
!arrays
 integer,allocatable :: kg1_k(:,:),kg2_k(:,:)
 real(dp),allocatable :: eig1_k(:),cg1_k(:,:),occ1_k(:)
 real(dp),allocatable :: eig2_k(:),cg2_k(:,:),occ2_k(:)

! *************************************************************************

 call wrtout(std_out,ABI_FUNC//": comparing "//TRIM(fname1)//" "//TRIM(fname2),"COLL")

 my_rank = xcomm_rank(comm)
 nproc   = xcomm_size(comm)
 sc_mode = xmpio_collective

 call hdr_read_from_fname(Hdr1,fname1,fform1,comm)
 call hdr_read_from_fname(Hdr2,fname2,fform2,comm)

 ABI_CHECK(fform1==fform2,"fform1 != fform2")
 ABI_CHECK(Hdr1%nsppol==Hdr2%nsppol,"nsppol1 != nsppol2")
 ABI_CHECK(Hdr1%nspinor==Hdr2%nspinor,"nspinor1 != nspinor2")
 ABI_CHECK(Hdr1%nkpt==Hdr2%nkpt,"nkpt1 != nkpt2")

 !call hdr_check(fform,fform0,hdr1,hdr2,"COLL",restart,restartpaw)

 iomode1 = iomode_from_fname(fname1)
 iomode2 = iomode_from_fname(fname1)

 ABI_CHECK(iomode1==iomode2,"iomode1 != iomode2")

 call wfk_open_read(Wfk1,fname1,formeig,iomode1,get_unit(),comm)
 call wfk_open_read(Wfk2,fname2,formeig,iomode2,get_unit(),comm)

 mband = Wfk1%mband
 ABI_CHECK(mband==Wfk2%mband,"different mband")

 ierr = 0
 do spin=1,Wfk1%nsppol
   do ik_ibz=1,Wfk1%nkpt
     npw_k    = Wfk1%Hdr%npwarr(ik_ibz)
     nband_k  = Wfk1%nband(ik_ibz,spin)
     ABI_CHECK(npw_k  ==Wfk2%Hdr%npwarr(ik_ibz),"different npw_k")
     ABI_CHECK(nband_k==Wfk2%nband(ik_ibz,spin),"different nband_k")

     mcg = npw_k*Hdr1%nspinor*nband_k

     ABI_MALLOC(eig1_k,((2*mband)**formeig*mband))
     ABI_MALLOC(occ1_k,(mband))
     ABI_MALLOC(kg1_k,(3,npw_k))
     ABI_MALLOC(cg1_k,(2,mcg))
     ABI_CHECK_ALLOC("out of memory in cg1_k") 

     ABI_MALLOC(eig2_k,((2*mband)**formeig*mband))
     ABI_MALLOC(occ2_k,(mband))
     ABI_MALLOC(kg2_k,(3,npw_k))
     ABI_MALLOC(cg2_k,(2,mcg))
     ABI_CHECK_ALLOC("out of memory in cg2_k") 

     ! Read the block of bands for this (k,s).
     call wfk_read_band_block(Wfk1,(/1,nband_k/),ik_ibz,spin,sc_mode,kg_k=kg1_k,cg_k=cg1_k,eig_k=eig1_k,occ_k=occ1_k)

     call wfk_read_band_block(Wfk2,(/1,nband_k/),ik_ibz,spin,sc_mode,kg_k=kg2_k,cg_k=cg2_k,eig_k=eig2_k,occ_k=occ2_k)

     if (ANY( ABS(cg1_k - cg2_k) > zero)) then
       ierr = ierr + 1
       MSG_WARNING("Difference in cg_k")
     end if
                                                                        
     if (ANY( ABS(kg1_k - kg2_k) > zero)) then
       ierr = ierr + 2
       MSG_WARNING("Difference in kg_k")
       !write(std_out,*)"kg1_k",kg1_k
       !write(std_out,*)"kg2_k",kg2_k
     end if
                                                                        
     if (ANY( ABS(eig1_k - eig2_k) > zero)) then 
       ierr = ierr + 4
       MSG_WARNING("Difference in eig_k")
       !write(std_out,*)"eig1_k",eig1_k
       !write(std_out,*)"eig2_k",eig2_k
     end if
                                                                        
     if (formeig==0) then
       if (ANY( ABS(occ1_k - occ2_k) > zero)) then
         ierr = ierr + 8
         MSG_WARNING("occ_k")
         write(std_out,*)"occ1_k",occ1_k
         write(std_out,*)"occ2_k",occ2_k
       end if
     end if
                                                                        
     write(msg,"(a,3(i0,2x))")" (ik_ibz, spin, ierr) ",ik_ibz,spin,ierr
     if (ierr/=0) then
       MSG_WARNING(TRIM(msg)//": FAILED")
     else
       call wrtout(std_out,TRIM(msg)//": OK","COLL")
     end if

     ABI_FREE(eig1_k)
     ABI_FREE(occ1_k)
     ABI_FREE(kg1_k)
     ABI_FREE(cg1_k)

     ABI_FREE(eig2_k)
     ABI_FREE(occ2_k)
     ABI_FREE(kg2_k)
     ABI_FREE(cg2_k)
   end do !ik_ibz
 end do !spin

 call wfk_close(Wfk1)
 call wfk_close(Wfk2)

 call hdr_free(Hdr1)
 call hdr_free(Hdr2)

end subroutine wfk_diff
!!***

!----------------------------------------------------------------------

END MODULE m_wfk
!!***
