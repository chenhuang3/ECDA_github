!{\src2tex{textfont=tt}}
!!****p* ABINIT/ioprof
!! NAME
!! ioprof
!!
!! FUNCTION
!! Tool for frofiling and and testing the IO routines used in abinit
!!
!! COPYRIGHT
!! Copyright (C) 2004-2014 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,delete_file,hdr_echo,hdr_free,hdr_read_from_fname
!!      herald,lower,m_header_init,wfk_check_wfkfile,wfk_create_wfkfile
!!      wfk_prof,wrtout,xmpi_bcast,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program ioprof

 use defs_basis
 use m_build_info
 use m_errors
 use m_xmpi
 use m_wfk
 use m_profiling_abi
 use m_header
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif

 use defs_abitypes,    only : hdr_type
 use m_fstrings,       only : lower
 use m_io_tools,       only : delete_file, file_exists

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ioprof'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,MAX_NFILES=50
 integer :: comm,my_rank,nprocs,iomode,formeig,ierr,fform
 integer :: ii,io,check_iomode,method,feg,ount,nband2read
 logical :: verbose=.FALSE.
 character(len=24) :: codename
 character(len=500) :: msg
 character(len=fnlen) :: new_fname
 type(hdr_type) :: Hdr
!arrays
 integer,parameter :: formeigs(2) = (/0,1/)
 integer,parameter :: io_modes(1) = (/IO_MODE_FORTRAN/)
 !integer,parameter :: io_modes(1) = (/IO_MODE_MPI/)
 !integer,parameter :: io_modes(1) = (/IO_MODE_ETSF/)
 !integer,parameter :: io_modes(2) = (/IO_MODE_FORTRAN, IO_MODE_MPI/)
 !integer,parameter :: io_modes(3) = (/IO_MODE_FORTRAN, IO_MODE_MPI, IO_MODE_ETSF/)
 real(dp) :: cwtimes(2)
 type(kvars_t),allocatable :: Kvars(:)
! ==========  INPUT FILE ==============
 character(len=fnlen) :: wfk_fname = ABI_NOFILE
 character(len=fnlen) :: hdr_fnames(MAX_NFILES) = ABI_NOFILE
 character(len=500) :: tasks=""
 namelist /CONTROL/ tasks, hdr_fnames, wfk_fname

! *************************************************************************

!call clib_mtrace(ierr)

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

 call xmpi_init()

!TODO perhaps I can avoid this boring initialization.
 call m_header_init(ierr)
 ABI_CHECK(ierr==0,"m_header_init returned ierr != 0")

 codename='IOPROF'//REPEAT(' ',17)
 call herald(codename,abinit_version,std_out)

 comm    = xmpi_world 
 my_rank = xcomm_rank(comm)
 nprocs  = xcomm_size(comm) 

!verbose = .TRUE.
 ount = dev_null; if (verbose) ount = std_out

!call abi_io_redirect(new_ab_out=ount,new_std_out=ount)

 if (nprocs>1) then
   MSG_ERROR("not programmed for parallel execution, Run it in sequential")
 end if

!Read input file
 read(std_in, NML=CONTROL)
!write(std_out, NML=CONTROL)

 call lower(tasks)

!replace "+" with white spaces
 do ii=1,LEN_TRIM(tasks)
   if (tasks(ii:ii) == "+") tasks(ii:ii) = " "
 end do
!%call str_replace_chars(tasks,"+"," ")
 tasks = " "//TRIM(tasks)//""

 if (nprocs > 1) then
   call xmpi_bcast(tasks,master,comm,ierr)
   call xmpi_bcast(wfk_fname,master,comm,ierr)
   call xmpi_bcast(hdr_fnames,master,comm,ierr)
 end if

!Benchmark
 if (INDEX(tasks," bench")>0) then
   formeig=0; nband2read=1500
   call wfk_prof(wfk_fname,formeig,nband2read,comm)
 end if

!Section with unitary tests.
 if (INDEX(tasks," utests")>0) then
!  
   do ii=1,COUNT(hdr_fnames/=ABI_NOFILE)
!    
!    Read the header from an external netcdf file
!    This trick is needed because the initialization of 
!    the header is *IMPOSSIBLE* if we don't have a full ABINIT input file!
!    
     call hdr_read_from_fname(Hdr,hdr_fnames(ii),fform,comm)
     ABI_CHECK(fform/=0,"fform==0")

     call hdr_echo(hdr,fform,4,unit=std_out)

     do feg=1,SIZE(formeigs)
       formeig = formeigs(feg)
       do io=1,SIZE(io_modes)
         iomode = io_modes(io) 

         new_fname = "NEW_WFK"
         if (iomode==IO_MODE_ETSF) new_fname = "NEW_WFK.nc"

         if (file_exists(new_fname)) then
           call delete_file(new_fname,ierr)
         end if

!        TODO
         if (formeig==1 .and. iomode==IO_MODE_ETSF) then
           MSG_WARNING("iomode==1 with ETSF_IO not coded yet")
           CYCLE
         end if

         ABI_DT_MALLOC(Kvars, (Hdr%nkpt))

         if (my_rank == master) then
           call wrtout(ount,"Calling wfk_create_wfkfile","COLL")
           call wfk_create_wfkfile(new_fname,Hdr,iomode,formeig,Kvars,cwtimes,xmpi_self)
!          call wfk_create_wfkfile(new_fname,Hdr,IO_MODE_FORTRAN,formeig,Kvars,cwtimes,xmpi_self)
           call wrtout(ount,"Done wfk_create_wfkfile","COLL")
         end if

         if (nprocs > 1) then
           MSG_ERROR("Not coded")
!          call kvars_mpicast(Kvars,master,comm)
         end if

         method = 0

         write(msg,"(3a,i2)")" Checking file: ",TRIM(new_fname),", with iomode = ",iomode
         call wrtout(std_out,msg,"COLL")

         call wfk_check_wfkfile(new_fname,Hdr,iomode,method,formeig,Kvars,cwtimes,comm,ierr)

         write(msg,'(2(a,i2),2(a,f8.2))')&
&         " Read with iomode: ",iomode,", nproc: ",nprocs,", cpu: ",cwtimes(1),", wall:",cwtimes(2)
         call wrtout(ount,msg,"COLL")

         if (ierr/=0) then
           write(msg,"(a,i0)")" wfk_check_wfkfile returned ierr ",ierr
           MSG_ERROR(msg)
         end if

!        If not netcdf file, try to read the file with the other mode that is compatible with it.
         check_iomode = -100
         if (iomode == IO_MODE_FORTRAN) check_iomode = IO_MODE_MPI
         if (iomode == IO_MODE_MPI)     check_iomode = IO_MODE_FORTRAN

         if (check_iomode /= -100) then
           write(msg,"(3a,i2)")" Trying to read file: ",TRIM(new_fname),", with check_iomode = ",check_iomode
           call wrtout(std_out,msg,"COLL")

           call wfk_check_wfkfile(new_fname,Hdr,check_iomode,method,formeig,Kvars,cwtimes,comm,ierr)

           write(msg,'(2(a,i2),2(a,f8.2))')&
&           " Read with iomode: ",iomode,", nproc: ",nprocs,", cpu: ",cwtimes(1),", wall:",cwtimes(2)
           call wrtout(ount,msg,"COLL")

           if (ierr/=0) then
             write(msg,"(a,i0)")" wfk_check_wfkfile returned ierr ",ierr
             MSG_ERROR(msg)
           end if
         end if

         ABI_DT_FREE(Kvars) 
       end do ! iomode
     end do ! formeig

     call hdr_free(Hdr)
   end do
 end if

 call wrtout(std_out,ch10//" Analysis completed.","COLL")

 call xmpi_end()

 end program ioprof
!!***
