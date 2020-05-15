!{\src2tex{textfont=tt}}
!!****f* ABINIT/ioarr
!!
!! NAME
!! ioarr
!!
!! FUNCTION
!! Read or write rho(r) or v(r), either ground-state or response-functions.
!! If ground-state, these arrays are real, if response-functions, these arrays are complex.
!! (in general, an array stored in unformatted form on a real space fft grid).
!! rdwr=1 to read, 2 to write
!!
!! This subroutine should be called only by one processor in the writing mode
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, MVer, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! (some may be output)
!! accessfil=
!!    0 for FORTRAN_IO
!!    3 for ETSF_IO
!!    4 for MPI_IO
!! dtset <type(dataset_type)>=all input variables for this dataset
!! fform=integer specification for data type:
!!   2 for wf; 52 for density; 102 for potential
!!   old format (prior to ABINITv2.0): 1, 51 and 101.
!! fildata=file name
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!!  if rdwr=1 , used to compare with the hdr of the read disk file
!!  if rdwr=2 , used as the header of the written disk file
!! mpi_enreg=information about MPI parallelization
!! rdwr=choice parameter, see above
!! rdwrpaw=1 only if rhoij PAW quantities have to be read (if rdwr=1)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! arr(ncplxfft,dtset%nspden)=array on real space grid, returned for rdwr=1, input for rdwr=2
!! etotal=total energy (Ha), returned for rdwr=1
!! === if rdwrpaw/=0 ===
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! PARENTS
!!      gstate,gw_tools,loop3dte,loper3,nonlinear,outscfcv,respfn,scfcv,scfcv3
!!      setup_positron,sigma,suscep
!!
!! CHILDREN
!!      etsf_io_low_close,etsf_io_low_open_modify,etsf_io_low_open_read
!!      etsf_io_main_get,etsf_io_main_put,hdr_check,hdr_free,hdr_io,hdr_io_etsf
!!      pawrhoij_copy,wffclose,wffopen,wrtout,xderiveread,xderiverrecend
!!      xderiverrecinit,xderivewrecend,xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ioarr(accessfil,arr,dtset,etotal,fform,fildata,hdr,mpi_enreg, &
&                ncplxfft,pawrhoij,rdwr,rdwrpaw,wvl_den)
 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_xmpi
 use m_wffile
 use m_errors
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 use m_header,   only : hdr_free, hdr_io_etsf, hdr_io, hdr_check
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_copy

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ioarr'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accessfil,ncplxfft,rdwr,rdwrpaw
 integer,intent(inout) :: fform
 real(dp),intent(inout) :: etotal
 character(len=fnlen),intent(in) :: fildata
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(wvl_denspot_type), intent(in) :: wvl_den
!arrays
 real(dp),intent(inout),target :: arr(ncplxfft,dtset%nspden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:) 

!Local variables-------------------------------
#ifdef HAVE_TRIO_ETSF_IO
 integer :: ncid
 logical :: lstat
 character(len=fnlen) :: file_etsf
 type(etsf_main), target :: main_folder
 type(etsf_io_low_error) :: error
#endif
!scalars
 integer :: accesswff,fform_dum,i,i1,i2,i3,ia,iarr,ierr,ind,ispden,me,me_fft
 integer :: restart,restartpaw,spaceComm,spaceComm_io
 integer :: zindex,zstart,zstop,n1,n2,n3
 character(len=500) :: message
 type(hdr_type) :: hdr0
 type(wffile_type) :: wff
!arrays
 real(dp),pointer :: my_density(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 restartpaw=0

!Check validity of arguments--only rho(r) (51,52) and V(r) (101,102) are presently supported

 if ( (fform-1)/2 /=25 .and. (fform-1)/2 /=50 ) then
   write(message,'(a,i10,a)')' Input fform=',fform,' not allowed.'
   MSG_BUG(message)
 end if

!Print input fform
 if ( (fform-1)/2==25 .and. rdwr==1) then
   message = ' ioarr: reading density data '
 else if ( (fform-1)/2==25 .and. rdwr==2) then
   message = ' ioarr: writing density data'
 else if ( (fform-1)/2==50 .and. rdwr==1) then
   message = ' ioarr: reading potential data'
 else if ( (fform-1)/2==50 .and. rdwr==2) then
   message = ' ioarr: writing potential data'
 end if
 call wrtout(std_out,message,'COLL')

 call wrtout(std_out,ABI_FUNC//': file name is '//TRIM(fildata),'COLL')

#ifdef HAVE_TRIO_ETSF_IO
 if (accessfil == 3) then ! Initialize filename in case of ETSF file.
   file_etsf = TRIM(fildata) // '-etsf.nc'
   call wrtout(std_out,'created file name for ETSF access '//TRIM(file_etsf),'COLL')
 end if
#endif

!Some definitions for MPI-IO access
 if (accessfil == 4) then
   accesswff=IO_MODE_MPI
   if (rdwr==1) then
     spaceComm=mpi_enreg%comm_cell
   else
     spaceComm=mpi_enreg%comm_fft
   end if
   me=xcomm_rank(spaceComm)
   if (mpi_enreg%nproc_fft>1) then
     me_fft=mpi_enreg%me_fft
     spaceComm_io=mpi_enreg%comm_fft
   else
     me_fft=0
     spaceComm_io=xmpi_self
   end if
 end if
 if (dtset%usewvl==1) then
   spaceComm=mpi_enreg%comm_cell
   me=xcomm_rank(spaceComm)
 end if

!=======================================
!Handle input from disk file
!=======================================

 if (rdwr==1) then
   if (accessfil == 0 .or. accessfil == 4) then
     if(accessfil == 4) then
       call WffOpen(accesswff,spaceComm,fildata,ierr,wff,0,me,tmp_unit,spaceComm_io)
       call hdr_io(fform_dum,hdr0,rdwr,wff)
!      Compare the internal header and the header from the file
       call hdr_check(fform,fform_dum,hdr,hdr0,'COLL',restart,restartpaw)

     else
       open (unit=tmp_unit,file=fildata,form='unformatted',status='old')
!      Initialize hdr0, thanks to reading of unwff1
       call hdr_io(fform_dum,hdr0,rdwr,tmp_unit)
!      Compare the internal header and the header from the file
       call hdr_check(fform,fform_dum,hdr,hdr0,'COLL',restart,restartpaw)
     end if
     etotal=hdr0%etot

!    NOTE : should check that restart is possible !!

!    Read data
     do ispden=1,dtset%nspden
       if(accessfil == 4) then
         call xderiveRRecInit(wff,ierr)
         call xderiveRead(wff,arr(1:ncplxfft,ispden),ncplxfft,spaceComm_io,ierr)
         call xderiveRRecEnd(wff,ierr)
       else
         read (tmp_unit) (arr(iarr,ispden),iarr=1,ncplxfft)
       end if
     end do

     if(accessfil == 4) then
       call wffclose(wff,ierr)
     else
       close (unit=tmp_unit)
     end if

#ifdef HAVE_TRIO_ETSF_IO
   else if ( accessfil == 3 ) then

!    Open the file
     call etsf_io_low_open_read(ncid, file_etsf, lstat, error_data = error)
     ETSF_CHECK_ERROR(lstat,error) 

!    Read the header
     call hdr_io_etsf(fform_dum, hdr0, rdwr, ncid)

!    Compare the internal header and the header from the file
     call hdr_check(fform, fform_dum, hdr, hdr0, 'COLL', restart, restartpaw)

!    Read the array
     if (fform==52) then ! density
       main_folder%density%data2D => arr
     else if (fform==102) then ! all potential forms!!!!
       main_folder%exchange_correlation_potential%data2D => arr
     end if

     call etsf_io_main_get(ncid, main_folder, lstat, error)
     ETSF_CHECK_ERROR(lstat,error) 

!    Close the file
     call etsf_io_low_close(ncid, lstat, error_data = error)
     ETSF_CHECK_ERROR(lstat,error) 
#endif

   else
     write(message,'(a,i0,a)')'Bad value for accessfil', accessfil, ' on read '
     MSG_BUG(message)
   end if

   call wrtout(std_out,ABI_FUNC//': data read from disk file '//TRIM(fildata),'COLL')

   etotal=hdr0%etot
!  Eventually copy (or distribute) PAW data
   if (rdwrpaw==1.and.restartpaw/=0) then
     if (size(hdr0%pawrhoij) /= size(pawrhoij)) then
       call pawrhoij_copy(hdr0%pawrhoij,pawrhoij,&
&       mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     else
       call pawrhoij_copy(hdr0%pawrhoij,pawrhoij)
     end if
   end if

   if (accessfil == 0 .or. accessfil == 3 .or. accessfil == 4) then
     call hdr_free(hdr0)
   end if

!  =======================================
!  Set up for writing data
!  =======================================
 else if (rdwr==2) then

!  In the wavelet case (isolated boundary counditions), the
!  arr array has a buffer that we need to remove.
   if (dtset%usewvl == 1) then
#ifdef HAVE_DFT_BIGDFT
     zindex = wvl_den%denspot%dpbox%nscatterarr(me, 3)
     if (wvl_den%denspot%rhod%geocode == 'F') then
       n1 = (wvl_den%denspot%dpbox%ndims(1) - 31) / 2
       n2 = (wvl_den%denspot%dpbox%ndims(2) - 31) / 2
       n3 = (wvl_den%denspot%dpbox%ndims(3) - 31) / 2
       zstart = max(15 - zindex, 0)
       zstop  = wvl_den%denspot%dpbox%nscatterarr(me, 2) + &
&       wvl_den%denspot%dpbox%nscatterarr(me, 4) - &
&       max(zindex + wvl_den%denspot%dpbox%nscatterarr(me, 2) &
&       - 2 * n3 - 15, 0)
     else
       MSG_ERROR('ioarr: WVL not implemented yet.')
     end if
     if (zstop - zstart + 1 > 0) then
!      Our slab contains (zstop - zstart + 1) elements
       ABI_ALLOCATE(my_density,((n1*2)*(n2*2)*(zstop-zstart),dtset%nspden))
!      We copy the data except the buffer to my_density
       ind = 0

       do i3 = zstart, zstop - 1, 1
         ia = (i3 - 1) * dtset%ngfft(1) * dtset%ngfft(2)
         do i2 = 0, 2 * n2 - 1, 1
           i = ia + (i2 + 14) * dtset%ngfft(1) + 14
           do i1 = 0, 2 * n1 - 1, 1
             i   = i + 1
             ind = ind + 1
             my_density(ind, :) = arr(i, :)
           end do
         end do
       end do
     else
       nullify(my_density)
     end if
#else
     BIGDFT_NOTENABLED_ERROR()
#endif
   end if

   if (accessfil == 0 .or. accessfil == 4) then
     if(accessfil == 4) then
       call WffOpen(accesswff,spaceComm,fildata,ierr,wff,0,me,tmp_unit)
       call hdr_io(fform,hdr,rdwr,wff)
     else
       open(unit=tmp_unit,file=fildata,form='unformatted',status='unknown')
!      Write header
       call hdr_io(fform,hdr,rdwr,tmp_unit)
     end if

!    Write actual data
     do ispden=1,dtset%nspden
       if(accessfil == 4) then
         call xderiveWRecInit(wff,ierr,me_fft)
         call xderiveWrite(wff,arr(1:ncplxfft,ispden),ncplxfft,spaceComm_io,ierr)
         call xderiveWRecEnd(wff,ierr,me_fft)
       else
         if (dtset%usewvl == 0) then
           write(tmp_unit) (arr(iarr,ispden),iarr=1,ncplxfft)
         else
           write(tmp_unit) (my_density(iarr,ispden),iarr=1,size(my_density, 1))
         end if
       end if
     end do

     if(accessfil == 4) then
       call WffClose(wff,ierr)
     else
       close (tmp_unit)
     end if

#ifdef HAVE_TRIO_ETSF_IO
   else if ( accessfil == 3 ) then
!    Open the file
     call etsf_io_low_open_modify(ncid, trim(file_etsf), lstat, error_data = error)
     ETSF_CHECK_ERROR(lstat,error) 

!    Write the header
     call hdr_io_etsf(fform, hdr, rdwr, ncid)

!    Write the array
     if (fform==52) then ! density
       if (dtset%usewvl == 0) then
         main_folder%density%data2D => arr
       else
         main_folder%density%data2D => my_density
       end if
     else if (fform==102) then ! all potential forms!!!!
       main_folder%exchange_correlation_potential%data2D => arr
     end if

     call etsf_io_main_put(ncid, main_folder, lstat, error)
     ETSF_CHECK_ERROR(lstat,error) 

!    Close the file
     call etsf_io_low_close(ncid, lstat, error_data = error)
     ETSF_CHECK_ERROR(lstat,error) 
#endif

   else
     write(message,'(a,i0,a)')'Bad value for accessfil', accessfil, ' on write '
     MSG_ERROR(message)
   end if

   if (dtset%usewvl == 1) then
     if (associated(my_density))  then
       ABI_DEALLOCATE(my_density)
     end if
   end if

   call wrtout(std_out,ABI_FUNC//': data written to disk file '//TRIM(fildata),'COLL')

 else
   write(message,'(a,i0,a)')'Called with rdwr = ',rdwr,' not allowed.'
   MSG_BUG(message)
 end if

 DBG_EXIT("COLL")

end subroutine ioarr
!!***
