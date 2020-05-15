!{\src2tex{textfont=tt}}
!!****p* ABINIT/mrggkk
!! NAME
!! mrggkk
!!
!! FUNCTION
!! This program merges a GS file and several 1WF or GKK files for
!! different q-vectors and perturbations.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2014 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!! GKK file structure is composed of header records and eigenvalue arrays,
!! in binary or ascii:
!!   GS header = hdr
!!   GS eigenvalues = eigen
!!   number of perturbations = ntot
!!   for each perturbation
!!      1WF header = hdr1
!!      1st order eigenvalues = eigen1
!!
!! CHILDREN
!!      abi_io_redirect,destroy_mpi_enreg,flush_unit,hdr_free,hdr_io,herald
!!      initmpi_seq,rwwf,wrtout,xmpi_end,xmpi_init
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,destroy_mpi_enreg,flush_unit,hdr_free,hdr_io,herald
!!      initmpi_seq,wfk_close,wfk_open_read,wfk_read_eigk,wrtout,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program mrggkk

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wfk
 use m_build_info
 use m_argparse

 use m_io_tools,        only : flush_unit, open_file
 use m_header,          only : hdr_free, hdr_io
 use m_mpinfo,          only : destroy_mpi_enreg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mrggkk'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
!scalars
 integer,parameter :: unit1wf=22,unitgkk=24,unitgs=21,unitout=23,formeig0=0,formeig1=1
 integer :: binascii,fform,headform,i1wf,igkk,ik_ibz,ios,spin,mband
 integer :: n1wf,ngkk,ntot,ntotgkk,rdwr,comm,iomode
 integer :: base,iband,jband,nband_k,rdwrout,ierr,ipos
 real(dp) :: tolgkk=tol6
 character(len=1),parameter :: comment="#"
 character(len=24) :: codename
 character(len=500) :: message
 character(len=fnlen) :: file1wf,filegkk,filegs,outfile
 type(args_t) :: args
 type(MPI_type) :: mpi_enreg
 type(hdr_type) :: hdr,hdr1
 type(wfk_t) :: GS_wfk,PH_wfk
!arrays
 real(dp),allocatable :: eig_k(:)

! *************************************************************************

!Change communicator for I/O (mandatory!)
 comm  = xmpi_world
 call abi_io_redirect(new_io_comm=comm)

!Initialize MPI
 call xmpi_init()

 ! parse command line arguments
 args = args_parser(); if (args%exit /= 0) goto 100

 iomode = IO_MODE_FORTRAN
#ifdef HAVE_MPI_IO
 iomode = IO_MODE_MPI
#endif

!Default for sequential use
 call initmpi_seq(mpi_enreg)

 codename='MRGGKK'//repeat(' ',18)

!write greating,read the file names, etc.
 call herald(codename,abinit_version,std_out)

 write(message,'(17a)')&
& ' Files file format: ',ch10,ch10,&
& '  Name of the output file',ch10,&
& '  Integer flag: 0 --> binary output,   1 --> ascii formatted output',ch10,&
& '  Name of the groud state wavefunction file WF',ch10,&
& '  Number of 1WF, of GKK files, and number of 1WF files in all the GKK files',ch10,&
& '  Names of the 1WF files...',ch10,&
& '  Names of the GKK files...',ch10,ch10,&
& ' Enter name of output file: '
 call wrtout(std_out,message,'COLL')

!get file with filenames and number of 1wf files
 read(*,'(a)') outfile
 ipos=INDEX(outfile,comment)
 if (ipos/=0) outfile=outfile(:ipos-1)

 read(*,*) binascii

 read(*,'(a)') filegs
 ipos=INDEX(filegs,comment)
 if (ipos/=0) filegs=filegs(:ipos-1)

 read(*,*) n1wf,ngkk,ntotgkk

 write(message,'(7a,i4,2a,i4,2a,i4,a)')&
& ' Output                     = ',trim(outfile),ch10,&
& ' Ground State file          = ',trim(filegs),ch10,&
& ' Number of 1WF files        = ',n1wf,ch10,&
& ' Number of GKK files        = ',ngkk,ch10,&
& ' Total Number of 1WF in GKK = ',ntotgkk,ch10
 call wrtout(std_out,message,'COLL')


!output without rewinding the file
 if (binascii == 0) then
   ! open output file
   ios = open_file(outfile,message,unit=unitout,form='unformatted')
   rdwrout = 6
 else if (binascii == 1) then
!  rdwrout=4 ! use for screen output and change writes of eigen to (*,*)
!  MJV 27/5/2008 removed 'new' constraint on gkk files: presume competent user!
   ios = open_file(outfile,message,unit=unitout,form='formatted')
   rdwrout = 4
 else if (binascii == 2) then
!  this is for simple "short" output of the matrices, without headers or imaginary part
   ios = open_file(outfile,message,unit=unitout,form='formatted')
   rdwrout = 4
 else
   MSG_ERROR(' binascii must be 0 or 1')
 end if

 ABI_CHECK(ios==0,message)
 rewind (unitout)

!-------------------------------------------------------
!now read and write information for GS file
!-------------------------------------------------------

!open GS wf file
 call wrtout(std_out,' normal input for GS file',"COLL")

 call wfk_open_read(GS_wfk,filegs,formeig0,iomode,unitgs,comm)

!Copy header of GS file to output.
 if (binascii /= 2) then 
   call hdr_io(GS_wfk%fform,GS_wfk%Hdr,rdwrout,unitout)
 end if

 call wrtout(std_out,' header echoed to output file',"COLL")

 ABI_MALLOC(eig_k,(GS_wfk%mband))

!Retrieve GS eigenvalues from GS wf file and echo to output
 do spin=1,GS_wfk%nsppol
   do ik_ibz=1,GS_wfk%nkpt
     nband_k = GS_wfk%nband(ik_ibz,spin)

     call wfk_read_eigk(GS_wfk,ik_ibz,spin,xmpio_single,eig_k)
!    
     if (binascii==0) then
       write(unitout) eig_k(1:nband_k)
     else
       write(unitout,*) eig_k(1:nband_k)
     end if
!    
   end do
 end do

 ABI_FREE(eig_k)

!Close GS wf file
 call wfk_close(GS_wfk)

 ntot = n1wf + ntotgkk
 if (binascii==0) then
   write (unitout) ntot
 else
   write (unitout,*) ntot
 end if

!-------------------------------------------------------
!now read and write information for 1WF files
!-------------------------------------------------------
 do i1wf=1,n1wf
!  for each 1wf file, get name...
   read(*,'(a)') file1wf
   ipos=INDEX(file1wf,comment)
   if (ipos/=0) file1wf=file1wf(:ipos-1)

!  open 1wf file
   call wrtout(std_out,' normal input for 1WF file ',"COLL")

   call wfk_open_read(PH_wfk,file1wf,formeig1,iomode,unit1wf,comm,Hdr_out=hdr1)

!  copy header of 1WF file to output
!  WARNING: cant use normal hdr_io because it rewinds the output file and destroys previous content.
   if (binascii /= 2) then
     call hdr_io(PH_wfk%fform,hdr1,rdwrout,unitout)
   else
     write (unitout,'(a,3E20.10)') "qpt ", hdr1%qptn
     write (unitout,'(a,I6)') "pertnum ", hdr1%pertcase
   end if

!  retrieve 1WF <psi_k+q | H | psi_k> from 1wf file and echo to output
   mband = maxval(hdr1%nband)
   headform=hdr1%headform
   ABI_MALLOC(eig_k,(2*mband*mband))

   ABI_CHECK(ALL(PH_wfk%nband == PH_wfk%nband(1,1)),"nband must be constant")

   do spin=1,hdr1%nsppol
     do ik_ibz=1,hdr1%nkpt
!      write(std_out,*) 'spin,ik_ibz = ', spin,ik_ibz
       nband_k = PH_wfk%nband(ik_ibz,spin)

       call wfk_read_eigk(PH_wfk,ik_ibz,spin,xmpio_single,eig_k)

       base = 0
       do jband=1,nband_k
         base = 2*(jband-1)*nband_k
         do iband=1,2*nband_k
           write(777,*) iband,jband,eig_k(base+iband)
         end do
       end do

       if (binascii==0) then
         write(unitout) eig_k(1:2*nband_k**2)
       else if (binascii==1) then
         write(unitout,*) eig_k(1:2*nband_k**2)
       else if (binascii==2) then
         do iband=1,nband_k
           do jband=1,nband_k
             if (abs(eig_k(2*nband_k*(iband-1)+2*(jband-1)+1))>tolgkk) then
               write(unitout,'(E18.7, 2x)', ADVANCE='NO') eig_k(2*nband_k*(iband-1)+2*(jband-1)+1)
             else
               write(unitout,'(I18, 2x)', ADVANCE='NO') 0
             end if 
             if (abs(eig_k(2*nband_k*(iband-1)+2*(jband-1)+2))>tolgkk) then
               write(unitout,'(E18.7, 2x)', ADVANCE='NO') eig_k(2*nband_k*(iband-1)+2*(jband-1)+2)
             else 
               write(unitout,'(I18, 2x)', ADVANCE='NO') 0
             end if 
           end do
           write(unitout,*)
         end do
         write(unitout,*)
       end if
!      
     end do
     if (binascii==2) write(unitout,'(2a)') ch10, ch10
   end do

   ABI_FREE(eig_k)

!  clean header to deallocate everything
   call hdr_free(hdr1)

   call wfk_close(PH_wfk)
 end do

!-------------------------------------------------------
!now read and write information for small GKK files
!-------------------------------------------------------
 do igkk=1,ngkk
!  
!  for each gkk file, get name...
   read(*,'(a)') filegkk
   ipos=INDEX(filegkk,comment)
   if (ipos/=0) filegkk=filegkk(:ipos-1)

!  open gkk file
   call wrtout(std_out,' normal input for GKK file',"COLL")

   if (open_file(filegkk,message,unit=unitgkk,form='unformatted',status='old') /= 0) then
     MSG_ERROR(message)
   end if
   rewind (unitgkk)

!  read in header of GS file and eigenvalues
!  could force a comparison of header with global header above for consistency
   call hdr_io(fform,hdr,5,unitgkk)

   mband = maxval(hdr%nband)
   ABI_MALLOC(eig_k,(mband))
   call wrtout(std_out,'mrggkk : try to reread GS eigenvalues','COLL')

   do spin=1,hdr%nsppol
     do ik_ibz=1,hdr%nkpt
       nband_k = hdr%nband(ik_ibz + (spin-1)* hdr%nkpt)
       read (unitgkk,IOSTAT=ierr) eig_k(1:nband_k)
       ABI_CHECK(ierr==0,'error reading eigen from gkk file')
     end do
   end do
   ABI_FREE(eig_k)

   read(unitgkk,IOSTAT=ierr) n1wf
   ABI_CHECK(ierr==0,'error reading n1wf record')

   ABI_MALLOC(eig_k,(2*mband*mband))
   do i1wf=1,n1wf
!    read in header of 1WF file
     rdwr = 5
     call hdr_io(fform,hdr1,rdwr,unitgkk)
     if (fform == 0) then
       write(message,'(a,i0,a)')' 1WF header number ',i1wf,' was mis-read. fform == 0'
       MSG_ERROR(message)
     end if

!    copy header of 1WF file to output
     if (binascii /= 2) then
       call hdr_io(fform,hdr1,rdwrout,unitout)
     else
       write (unitout,'(a,3E20.10)') "qpt ", hdr1%qptn
       write (unitout,'(a,I6)') "pertnum ", hdr1%pertcase
     end if

!    retrieve 1WF <psi_k+q | H | psi_k> from gkk file and echo to output
     do spin=1,hdr1%nsppol
       do ik_ibz=1,hdr1%nkpt
         nband_k = hdr%nband(ik_ibz + (spin-1)* hdr1%nkpt)
         read (unitgkk,IOSTAT=ierr) eig_k(1:2*nband_k**2)
         if (ierr /= 0) write (std_out,*) 'error reading eigen2 from gkk file',spin,ik_ibz

         if (binascii==0) then
           write (unitout) eig_k(1:2*nband_k**2)
         else if (binascii==1) then
           write (unitout,*) eig_k(1:2*nband_k**2)
         else if (binascii==2) then
           do iband=1,nband_k
             do jband=1,nband_k
               if (abs(eig_k(2*nband_k*(iband-1)+2*(jband-1)+1))>tolgkk) then
                 write(unitout,'(E18.7, 2x)', ADVANCE='NO') eig_k(2*nband_k*(iband-1)+2*(jband-1)+1)
               else
                 write(unitout,'(I18, 2x)', ADVANCE='NO') 0
               end if 
               if (abs(eig_k(2*nband_k*(iband-1)+2*(jband-1)+2))>tolgkk) then
                 write(unitout,'(E18.7, 2x)', ADVANCE='NO') eig_k(2*nband_k*(iband-1)+2*(jband-1)+2)
               else 
                 write(unitout,'(I18, 2x)', ADVANCE='NO') 0
               end if 
             end do
             write(unitout,*)
           end do
           write(unitout,*)
         end if
!        
       end do
       if (binascii==2) write(unitout,'(2a)') ch10, ch10
     end do
     call hdr_free(hdr1)
   end do !  end loop over 1wf segments in small gkk file

   ABI_FREE(eig_k)

   close (unitgkk)
   call hdr_free(hdr)
 end do !end loop over small gkk files

 close (unitout)

 write(message,'(2a)')ch10,' Done'
 call wrtout(std_out,message,'COLL')

 call flush_unit(std_out)

 call destroy_mpi_enreg(mpi_enreg)
 100 call xmpi_end()

 end program mrggkk
!!***
