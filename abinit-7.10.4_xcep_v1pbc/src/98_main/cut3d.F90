!{\src2tex{textfont=tt}}
!!****p* ABINIT/cut3d
!! NAME
!! cut3d
!!
!! FUNCTION
!! Main routine for the analysis of the density and potential files,
!! as well as other files with the ABINIT header.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2014 ABINIT group (GMR, RC, LSI, XG, NCJ, JFB, MCote, LPizzagalli)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!  (main program)
!!
!! NOTES
!! natom = number of atoms in the unit cell
!! nr1,nr2,nr3 = grid size (nr1 x nr2 x nr3 = filrho dimension)
!! ntypat = number of atom types
!! ucvol = unit cell volume (> 0)
!! densfileformat = flag for the format of the density file:
!!         0 = ASCII
!!         1 = binary
!!         2 = binary (ETSF)
!! denval = density value exported by interpol, to be wrote in the output file
!! filrho = name of the density file (ASCII or binary)
!! filtau = name of the atomic position file (Xmol format)
!! acell = unit cell parameters
!! rprim = orientation of the unit cell axes
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,destroy_mpi_enreg,flush_unit,hdr_free,hdr_io
!!      hdr_io_etsf,herald,hirsh,initmpi_seq,lineint,localorb_s,metric,planeint
!!      pointint,rrho,rtau,timein,volumeint,wffile,wffopen,xmpi_end,xmpi_init
!!      xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program cut3d

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_wffile
 use m_build_info
 use m_xmpi
#if defined FC_NAG
 use f90_unix_proc
#endif

 use m_header,          only : hdr_free, hdr_io_etsf, hdr_io
 use m_mpinfo,          only : destroy_mpi_enreg
 use m_io_tools,        only : flush_unit, file_exists, open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cut3d'
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
 use interfaces_83_cut3d
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
 character(len=*),parameter :: INPUTfile='cut.in'
 character(len=1) :: outputchar,blank=' '
!scalars
 integer,parameter :: paral_kgb0=0
 integer :: densfileformat,exchn2n3d,fform0,gridshift1,gridshift2,gridshift3,i1,i2,i3
 integer :: iatom,ierr,ifiles,ii,ii1,ii2,ii3,index,iprompt,ir1,ir2,ir3,ispden
 integer :: itask,jfiles,mfiles,natom,nfiles,nr1,nr2
 integer :: nr3,nr1_stored,nr2_stored,nr3_stored,nrws,nspden,nspden_stored,ntypat,rdwr,unitfi
 real(dp) :: dotdenpot,maxmz,normz,sumdenpot,ucvol,xm,xnow,xp,ym,ynow,yp,zm,znow,zp,tcpui,twalli
 character(len=24) :: codename
 character(len=50) :: chr_inputfname
 character(len=fnlen) :: filetsf,filnam,filrho,filrho_tmp,filtau
 type(hdr_type) :: hdr
 type(MPI_type) :: mpi_enreg
 type(wffile_type) :: wff
!arrays
 integer, allocatable :: isdenpot(:)
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3),shift_tau(3),tsec(2)
 real(dp) :: xcart2(3),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: grid(:,:,:),grid_full(:,:,:,:),grid_full_stored(:,:,:,:,:),gridtt(:,:,:)
 real(dp),allocatable :: tau2(:,:),xcart(:,:),xred(:,:),rhomacu(:,:),gridmz(:,:,:),gridmy(:,:,:),gridmx(:,:,:)
 character(len=fnlen),allocatable :: filrho_stored(:)
 character(len=500) :: message

!******************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()

 call timein(tcpui,twalli)

!Default for sequential use
 call initmpi_seq(mpi_enreg)

!Signal MPI I/O compilation has been activated
#if defined HAVE_MPI_IO
 if(xmpi_paral==0)then
   write(message,'(3a)')&
&   ' In order to use MPI_IO, you must compile with the MPI flag ',ch10,&
&   ' Action : recompile your code with different CPP flags.'
   MSG_ERROR(message)
 end if
#endif

!Other values of mpi_enreg are dataset dependent, and should NOT be initialized inside cut3d.F90.

 codename='CUT3D '//repeat(' ',18)
 call herald(codename,abinit_version,std_out)

!BIG LOOP on files
 mfiles=10
 ABI_MALLOC(isdenpot,(mfiles))
 isdenpot=0
 ABI_MALLOC(filrho_stored,(mfiles))

 do ifiles=1,mfiles

!  Get name of density file
   write(std_out,*)
   write(std_out,*) ' What is the name of the 3D function (density, potential or wavef) file ?'
   read(5,'(a)')filrho
   filrho_tmp=adjustl(filrho)
   do ii=1,len_trim(filrho_tmp)
     if(filrho_tmp(ii:ii)==blank)then
       filrho=trim(filrho_tmp(1:ii-1))
       exit
     end if
   end do
   write(std_out,*) ' => Your 3D function file is : ',trim(filrho)
   write(std_out,*)

!  Checking the existence of data file
   if (.not.file_exists(filrho)) then
     MSG_ERROR('Missing data file: '//TRIM(filrho))
   end if

!  Get its type
   write(std_out,*) ' Does this file contain formatted 3D ASCII data (=0)  '
   write(std_out,*) '  or unformatted binary header + 3D data        (=1) ?'
   write(std_out,*) '  or ETSF binary                                (=2) ?'
   read(5,*)densfileformat

!  Treat the different cases : formatted or unformatted
   if (densfileformat==1 .or. densfileformat == 2) then

     if (densfileformat == 1) then
       write(std_out,*) ' 1 => Your file contains unformatted binary header + 3D data'
     else
       write(std_out,*) ' 2 => Your file contains ETSF data'
     end if
     write(std_out,*) ' The information it contains should be sufficient.'

     if (densfileformat == 1) then
       if (open_file(filrho,message,unit=19,form='unformatted',status='old') /=0) then
         MSG_ERROR(message)
       end if 
       write(std_out,'(a,a,a,i4)')'  cut3d : read file ',trim(filrho),' from unit number 19.'
       write(std_out,*)
!      Read the header
       rdwr=1 ; unitfi=19
       call hdr_io(fform0,hdr,rdwr,unitfi)
       wff%unwff=19
       wff%fname=filrho
       wff%accesswff=IO_MODE_FORTRAN
     else
!      We remove -etsf.nc from the file name.
       write(filetsf,"(A)") filrho(1:len(trim(filrho)) - 8)
!      Note that the MPI information are dummy, to avoid errors when consistency checks inside wffopen
       call WffOpen(3, mpi_enreg%comm_world, filetsf, ierr, wff, mpi_enreg%me, mpi_enreg%me, 19)
       write(std_out,'(a,a,a,i4)' )'  cut3d : read file ',trim(filrho),'.'
       write(std_out,*)
!      Read the header
       rdwr=1; unitfi=wff%unwff
       call hdr_io_etsf(fform0,hdr,rdwr,unitfi)
     end if

!    Echo part of the header
     rdwr=4; unitfi=6
     call hdr_io(fform0,hdr,rdwr,unitfi)

     natom=hdr%natom
     nr1=hdr%ngfft(1)
     nr2=hdr%ngfft(2)
     nr3=hdr%ngfft(3)
     nspden=hdr%nspden
     ntypat=hdr%ntypat
     rprimd(:,:)=hdr%rprimd(:,:)

!    Need to know natom in order to allocate xcart
     ABI_MALLOC(xcart,(3,natom))
     ABI_MALLOC(xred,(3,natom))
     xred(:,:)=hdr%xred(:,:)
     call xred2xcart(natom,rprimd,xcart,xred)

     if(fform0>50)then
       ispden=0
       if(nspden/=1)then
         write(std_out,'(a)' )' '
         write(std_out,'(a)' )' * This file contains more than one spin component,'
         write(std_out,'(a,i3,a)' )'  (indeed, nspden=',nspden,' )'
         write(std_out,'(a)' )'  Some of the tasks that you will define later will concern all spin components.'
         write(std_out,'(a)' )'  Others tasks might require you to have chosen among the following :'
       end if
       if(nspden==2)then
         write(std_out,'(a)' )'   ispden= 0 ==> Total density'
         write(std_out,'(a)' )'   ispden= 1 ==> spin-up density'
         write(std_out,'(a)' )'   ispden= 2 ==> spin-down density'
         write(std_out,'(a)' )'   ispden= 3 ==> spin-polarization (or magnetization) density'
         write(std_out,'(a)' )'                 spin up - spin down difference.'
       end if
       if(nspden==4)then
         write(std_out,'(a)' )'   ispden= 0 ==> Total density'
         write(std_out,'(a)' )'   ispden= 1 ==> magnetization in the x direction'
         write(std_out,'(a)' )'   ispden= 2 ==> magnetization in the y direction'
         write(std_out,'(a)' )'   ispden= 3 ==> magnetization in the z direction'
         write(std_out,'(a)' )'   ispden= 4 might be used to plot the magnetization (3D) in the XCrysDen format,'
       end if
       if(nspden/=1)then
         write(std_out,*)'  Please define ispden :'
         read(5,*)ispden
         write(std_out,'(a,i3)' )' You entered ispden=',ispden
       end if
     end if

   else if(densfileformat==0)then

     write(std_out,*) ' 0 => Your file contains formatted 3D ASCII data'
     write(std_out,*) ' The complementary information is taken from ',trim(INPUTfile)

!    Checking the existence of input file
     if (.not.file_exists(INPUTfile)) then
       MSG_ERROR('Missing input file: '//TRIM(INPUTfile))
     end if

!    Read in the input file INPUTfile
     write(std_out,*)
     write(std_out,*) 'READING FROM FILE ',INPUTfile

     if (open_file(INPUTfile,message,unit=32,status='old') /= 0) then
       MSG_ERROR(message)
     end if
     read(32,'(a)') filtau

!    Checking the existence of atomic position file
     if (.not.file_exists(filtau)) then
       MSG_ERROR('Missing atomic position file: '//TRIM(filtau))
     end if
     write(std_out,*) 'atomic position file (Xmol format):',trim(filtau)

!    Read cell parameters
     read(32,*) acell(1),acell(2),acell(3)
     do ii=1,3
       if (acell(ii) <= zero) then
         write(message,'(a,i0,a,f8.4)')' Invalid value for acell(',ii,'): ',acell(ii)
         MSG_ERROR(message)
       end if
     end do
     read(32,*) rprim

     do ii=1,3
       rprimd(:,ii)=rprim(:,ii)*acell(ii)
     end do

!    FFT grid, number of atoms, number of atom types
     read(32,*) nr1,nr2,nr3
     read(32,*) natom
     read(32,*) ntypat

     close(32)

!    Need to know natom in order to allocate xcart
     ABI_MALLOC(xcart,(3,natom))

     call rtau(filtau,xcart,natom,ntypat)

     write(std_out,*)

!    By default there is only one spin component, and one works with a total density
     nspden=1; ispden=0; fform0=52
     
     wff%unwff=19
     wff%fname=filrho
     if (open_file(filrho,message,unit=wff%unwff,form='formatted',status='old') /= 0) then
       MSG_ERROR(message)
     end if

   else
     write(message,'(a,i0)')' Value for density file format is invalid: ',densfileformat
     MSG_ERROR(message)
   end if

   write(std_out,*)
   write(std_out,*) '==========================================================='
   write(std_out,*)

!  Echo the value of different input parameters
   write(std_out,*)'ECHO important input variables ...'
   write(std_out,*)
   write(std_out,*) ' Dimensional primitive vectors (ABINIT equivalent : rprimd):'
   write(std_out,'(3es16.6)' ) rprimd(1:3,1)
   write(std_out,'(3es16.6)' ) rprimd(1:3,2)
   write(std_out,'(3es16.6)' ) rprimd(1:3,3)

!  Compute ucvol and test the non-collinearity of rprimd vectors.
   call metric(gmet,gprimd,dev_null,rmet,rprimd,ucvol)

   write(std_out,'(a,3i5)' ) '  Grid density (ABINIT equivalent : ngfft): ',nr1,nr2,nr3
   write(std_out,*) ' Number of atoms       :',natom
   write(std_out,*) ' Number of atomic types:',ntypat

   write(std_out,*)
   write(std_out,*) '  #    Atomic positions (cartesian coordinates - Bohr)'
   do iatom=1,natom
     write(std_out,'(i4,3es16.6)' )iatom,xcart(1:3,iatom)
   end do
   write(std_out,*)

!  ------------------------------------------------------------------------
!  Branching : either WF file, or DEN/POT file.

   if((densfileformat==1 .or. densfileformat == 2) .and. fform0<50)then
     write(std_out,*)' This file is a WF file. '
     isdenpot(ifiles)=0
     iprompt = 0 ! this needs to be initialized, as it is used after the loop on files...

     exchn2n3d=0

     write(std_out,*)" If you want to analyze one wavefunction,                   type  0 "
     write(std_out,*)" If you want to construct Wannier-type Localized Orbitals,  type  2 "
     read(*,*)ii1
     write(std_out,*)" You typed ",ii1

     if(ii1==0)then
!      MG: Close wff here because we are gonna reopen fname in wffile with the new routines wfk_open_read.
!      A bit dirty but the problem is wffile_type that is a programming sin!
       close(wff%unwff)
       call wffile(wff%fname,hdr%ecut_eff,exchn2n3d,hdr%istwfk,hdr%kptns,natom,hdr%nband,hdr%nkpt,hdr%npwarr,&
&       nr1,nr2,nr3,hdr%nspinor,hdr%nsppol,ntypat,paral_kgb0,rprimd,xcart,hdr%typat,hdr%znucltypat)

     else if(ii1==2)then
!      Read the name of the input file name :
       read(5,'(a)')chr_inputfname
       call localorb_S(chr_inputfname,hdr%ecut_eff,exchn2n3d,hdr%headform,hdr%istwfk,hdr%kptns,&
&       natom,hdr%nband,hdr%nkpt,hdr%npwarr,&
&       nr1,nr2,nr3,hdr%nspinor,hdr%nsppol,ntypat,paral_kgb0,rprimd,xcart,hdr%typat,hdr%znucltypat)

       write(std_out,*)" "
       write(std_out,*)" ###################################################### "
       write(std_out,*)" "
       write(std_out,*)" Localized orbital files fort.1**1 for spin up "
       write(std_out,*)"                     and fort.1**2 for spin dn written."
       write(std_out,*)" "
       write(std_out,*)" ###################################################### "
     else
       write(std_out,*)" Option ",ii1," is not allowed => stop "
     end if

     call hdr_free(hdr)

!    -------------------------------------------------------------------------
   else ! This is a DEN/POT file

!    This should become a subroutine
     write(std_out,*)' This file is a Density or Potential file '
     isdenpot(ifiles)=1
!    FJ: Pb on ibm. Anyway, there is no chance to have a header for this formatted file.
!    In the wavelet case (with isolated boundary conditions), ngfft is buffered.
!    if (hdr%usewvl == 1) then
!    nr1 = nr1 - 31
!    nr2 = nr2 - 31
!    nr3 = nr3 - 31
!    end if

!    Read the function on the 3D grid
     ABI_MALLOC(grid,(nr1,nr2,nr3))
     ABI_MALLOC(grid_full,(nr1,nr2,nr3,nspden))
     ABI_MALLOC(gridtt,(nr1,nr2,nr3))
     ABI_MALLOC(gridmx,(nr1,nr2,nr3))
     ABI_MALLOC(gridmy,(nr1,nr2,nr3))
     ABI_MALLOC(gridmz,(nr1,nr2,nr3))

     call rrho(densfileformat,grid_full,nr1,nr2,nr3,nspden,wff)

!    Do not forget that the first sub-array of a density file is the total density,
!    while the first sub-array of a potential file is the spin-up potential
     if(fform0==51 .or. fform0==52)then   ! Density case

!      gridtt= grid --> Total density or potential.
!      gridmx= grid --> spin-Up density, or magnetization density in X direction.
!      gridmy= grid --> spin-Down density, or magnetization density in Y direction.
!      gridmz= grid --> spin-polarization density (Magnetization),
!      or magnetization density in Z direction.
       gridtt(:,:,:)=grid_full(:,:,:,1)
       if(nspden==2)then
         gridmx(:,:,:)=grid_full(:,:,:,2)
         gridmy(:,:,:)=grid_full(:,:,:,1)-grid_full(:,:,:,2)
         gridmz(:,:,:)=-grid_full(:,:,:,1)+two*grid_full(:,:,:,2)
       else if(nspden==4)then
         gridmx(:,:,:)=grid_full(:,:,:,2)
         gridmy(:,:,:)=grid_full(:,:,:,3)
         gridmz(:,:,:)=grid_full(:,:,:,4)
       end if

       if(nspden==1)then
         grid(:,:,:)=grid_full(:,:,:,1)
       else
         if(ispden==0)then
           grid(:,:,:)=gridtt(:,:,:)
         else if(ispden==1)then
           grid(:,:,:)=gridmx(:,:,:)
         else if(ispden==2)then
           grid(:,:,:)=gridmy(:,:,:)
         else if(ispden==3)then
           grid(:,:,:)=gridmz(:,:,:)
!          if(ispden==0)then
!          grid(:,:,:)=grid_full(:,:,:,1)
!          else if(ispden==1)then
!          grid(:,:,:)=grid_full(:,:,:,2)
!          else if(ispden==2)then
!          grid(:,:,:)=grid_full(:,:,:,1)-grid_full(:,:,:,2)
!          else if(ispden==-1)then
!          grid(:,:,:)=-grid_full(:,:,:,1)+two*grid_full(:,:,:,2)
         else if(ispden==4)then
           write(std_out,*) ' '
         else
           write(message,'(a,i0)')' bad ispden value = ',ispden
           MSG_ERROR(message)
         end if
       end if

     else if(fform0==101 .or. fform0==102)then    ! Potential case
       if(ispden==0)then
         grid(:,:,:)=grid_full(:,:,:,1)
       else if(ispden==1 .or. ispden==2)then
         grid(:,:,:)=grid_full(:,:,:,ispden)
       else
         write(message,'(a,i0)')' bad ispden value = ',ispden
         MSG_ERROR(message)
       end if
       gridtt = grid
     end if

     write(std_out,*)
     write(std_out,*) ' 3D function was read. Ready for further treatment.'
     write(std_out,*)
     write(std_out,*) '==========================================================='
     write(std_out,*)

!    ------------------------------------------------------------------------

!    At this moment all the input is done
!    The code knows the geometry of the system,
!    and the data file (electron density, potential, etc).
!    It will further calculate the electron density by interpolation in
!    a point, along a line or in a plane.

     do
       do
         write(std_out,*) ' What is your choice ? Type:'
         write(std_out,*) '  0 => exit'
         write(std_out,*) '  1 => point  (interpolation of data for a single point)'
         write(std_out,*) '  2 => line   (interpolation of data along a line)'
         write(std_out,*) '  3 => plane  (interpolation of data in a plane)'
         write(std_out,*) '  4 => volume (interpolation of data in a volume)'
         write(std_out,*) '  5 => 3D formatted data (output the bare 3D data - one column)'
         write(std_out,*) '  6 => 3D indexed data (bare 3D data, preceeded by 3D index)'
         write(std_out,*) '  7 => 3D Molekel formatted data '
         write(std_out,*) '  8 => 3D data with coordinates (tecplot ASCII format)'
         write(std_out,*) '  9 => output .xsf file for XCrysDen'
         write(std_out,*) ' 10 => output .dx file for OpenDx'
         write(std_out,*) ' 11 => compute atomic charge using the Hirshfeld method'
         write(std_out,*) ' 12 => NetCDF file'
         write(std_out,*) ' 14 => Gaussian/cube wavefunction module'
         read(*,*) itask
         write(std_out,'(a,a,i2,a)' ) ch10,' Your choice is ',itask,ch10

         if( 5<=itask .and. itask<=10 .or. itask==12 .or. itask==14 )then
           write(std_out,*) ch10,'  Enter the name of an output file:'
           read(*,*) filnam
           write(std_out,*) '  The name of your file is : ',trim(filnam)
         end if

         select case(itask)

           case(1) ! point calculation
             call pointint(gridtt,gridmx,gridmy,gridmz,nr1,nr2,nr3,nspden,rprimd)
             exit

           case(2) ! line calculation
             call lineint(gridtt,gridmx,gridmy,gridmz,nr1,nr2,nr3,nspden,rprimd)
             exit

           case(3) ! plane calculation
             call planeint(gridtt,gridmx,gridmy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,xcart)
             exit

           case(4) ! volume calculation
             write(std_out,*) ' Enter volume calculation'
             call volumeint(gridtt,gridmx,gridmy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,xcart)
             exit

           case(5)
!            Rewrite the data on a formatted file, just in one (or four) column(s)
             if (open_file(filnam,message,unit=31,status='unknown') /= 0) then
               MSG_ERROR(message)
             end if

             if(nspden==1)then
               do i3=1,nr3
                 do i2=1,nr2
                   do i1=1,nr1
                     write(31,'(4(es22.12))') grid(i1,i2,i3)
                   end do
                 end do
               end do
             else
               do i3=1,nr3
                 do i2=1,nr2
                   do i1=1,nr1
                     write(31,'(4(es22.12))') gridtt(i1,i2,i3), gridmx(i1,i2,i3), gridmy(i1,i2,i3), gridmz(i1,i2,i3)
                   end do
                 end do
               end do
             end if
             close(31)
             exit

           case(6)
!            Rewrite the data on a formatted file, 3D index + density
             if (open_file(filnam,message,unit=31,status='unknown') /= 0) then
               MSG_ERROR(message)
             end if

             if(nspden==1)then
               write(31,*)'   i1    i2    i3      data '
               do i3=1,nr3
                 do i2=1,nr2
                   do i1=1,nr1
                     write(31,'(3i6,4(es24.14))') i1,i2,i3,grid(i1,i2,i3)
                   end do
                 end do
               end do
             else
               if(nspden==2)then
                 write(31,*)'   i1    i2    i3     non-spin-polarized spin up  spin down  difference  '
               else if(nspden==4)then
                 write(31,*)'   i1    i2    i3     non-spin-polarized   x       y      z   '
               end if
               do i3=1,nr3
                 do i2=1,nr2
                   do i1=1,nr1
                     write(31,'(3i6,4(es24.14))') i1,i2,i3,gridtt(i1,i2,i3),gridmx(i1,i2,i3),gridmy(i1,i2,i3),gridmz(i1,i2,i3)
                   end do
                 end do
               end do
             end if ! nspden
             close(31)
             exit

           case(7)
             if (open_file(filnam,message,unit=31,form='unformatted') /= 0) then
               MSG_ERROR(message)
             end if

             xm=0 ; xp=rprimd(1,1)*Bohr_Ang
             ym=0 ; yp=rprimd(2,2)*Bohr_Ang
             zm=0 ; zp=rprimd(3,3)*Bohr_Ang
             write(std_out,'(/,a,/)' )' Extremas (x,y,z) of the cube in which the molecule is placed, in Angstroms'
             write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
             write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1,nr2,nr3
             write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', nr1*nr2*nr3
             write(31) xm,xp,ym,yp,zm,zp,nr1,nr2,nr3
             ABI_MALLOC(rhomacu,(nr1,nr2))
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   rhomacu(i1,i2)=grid(i1,i2,i3)
                 end do
               end do
               write(31) rhomacu(:,:)
             end do
             close(31)
             exit

           case (8)
             if (open_file(filnam,message,unit=31,form='formatted') /= 0) then
               MSG_ERROR(message)
             end if

             write(std_out,'(/,a,/)' )' Extremas (x,y,z) of the cube in which the molecule is placed, in Angstroms'
             write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
             write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1,nr2,nr3
             write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', nr1*nr2*nr3
             write(31,'(a)') 'TITLE = "  " '
             write(31,'(a)') 'VARIABLES = "X"  "Y"  "Z" (all three in Angstrom)  "DENSITY or POTENTIAL" (atomic units) '
             write(31,'(3(a,i6),a)') 'ZONE I=',nr1, ' J=', nr2, ' K=', nr3, ' F=POINT'
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   xnow = rprimd(1,1)*(i1-1)/nr1 + rprimd(1,2)*(i2-1)/nr2 + rprimd(1,3)*(i3-1)/nr3
                   ynow = rprimd(2,1)*(i1-1)/nr1 + rprimd(2,2)*(i2-1)/nr2 + rprimd(2,3)*(i3-1)/nr3
                   znow = rprimd(3,1)*(i1-1)/nr1 + rprimd(3,2)*(i2-1)/nr2 + rprimd(3,3)*(i3-1)/nr3
                   write(31,'(4es22.15)') Bohr_Ang*xnow, Bohr_Ang*ynow, Bohr_Ang*znow, grid (i1,i2,i3)
                 end do
               end do
             end do
             close(31)
             exit

           case (9)
             if (open_file(filnam,message,unit=31,form='formatted') /= 0) then
               MSG_ERROR(message)
             end if
             xm=0 ; xp=rprimd(1,1)*Bohr_Ang
             ym=0 ; yp=rprimd(2,2)*Bohr_Ang
             zm=0 ; zp=rprimd(3,3)*Bohr_Ang
             write(std_out,'(/,a,/)' )' Extremas (x,y,z) of the cube in which the molecule is placed, in Angstroms'
             write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
             write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1+1,nr2+1,nr3+1
             write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', (nr1+1)*(nr2+1)*(nr3+1)
             write(std_out,*) '  znucl = ', hdr%znucltypat, ' type = ', hdr%typat, ' ntypat = ', ntypat

             gridshift1 = 0
             gridshift2 = 0
             gridshift3 = 0
             write(std_out,*) 'Do you want to shift the grid along the x,y or z axis (y/n)?'
             write(std_out,*)
             shift_tau(:) = zero
             read (*,*) outputchar
             if (outputchar == 'y' .or. outputchar == 'Y') then
               write(std_out,*) 'Give the three shifts (x,y,z < ',nr1,nr2,nr3,') :'
               write(std_out,*)
               read (*,*) gridshift1, gridshift2, gridshift3
               shift_tau(:) = gridshift1*rprimd(:,1)/(nr1+1) + gridshift2*rprimd(:,2)/(nr2+1) + gridshift3*rprimd(:,3)/(nr3+1)
             end if
!            
!            Generate translated coordinates to match density shift
!            
             ABI_MALLOC(tau2,(3,natom))
             do iatom = 1,natom
               tau2(:,iatom) = xcart(:,iatom) - shift_tau(:)
             end do
!            ################################################################### (LD)
!            Option only available for "xcrysden" format as documented at the beginning
             if (ispden==4) then
!              It is necessary to know previously how many atoms will be used.
!              in order to plot the necessary magnetization arrows only.
               write(std_out,*)'Is it possible to decrease the number of arrows in order to improve the'
               write(std_out,*)'visualization in the screen, and decrease the size of the xcrysden output file.'
               write(std_out,*)'How many arrows would you like to skip? 0 = take all. 1 = skip every other point...'
               read (*,*) nrws
               nrws=nrws+1
               index=natom
               maxmz=0.0
               do i1=1,nr1,nrws
                 do i2=1,nr2,nrws
                   do i3=1,nr3,nrws
                     normz=gridmx(i1,i2,i3)**2+gridmy(i1,i2,i3)**2+gridmz(i1,i2,i3)**2
                     if(normz > maxmz) maxmz=normz
                   end do
                 end do
               end do
               if(abs(maxmz)<tol10)then
                 write(std_out,*)' At least, one of the components must differ from zero.'
                 MSG_ERROR("stopping now")
               end if
               do i1=1,nr1,nrws
                 do i2=1,nr2,nrws
                   do i3=1,nr3,nrws
                     normz=gridmx(i1,i2,i3)**2+gridmy(i1,i2,i3)**2+gridmz(i1,i2,i3)**2
                     if(0.1*maxmz <= normz) index=index+1
                   end do
                 end do
               end do

               write(31,'(1X,A)') 'CRYSTAL'
               write(31,'(1X,A)') 'PRIMVEC'
               do i1 = 1,3
                 write(31,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(i2,i1), i2=1,3)
               end do
               write(31,'(1X,A)') 'PRIMCOORD'
               write(31,*) index, '1'
!              
!              write out atom types and positions
!              
               do iatom = 1,natom
                 write(31,'(i9,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),Bohr_Ang*tau2(1:3,iatom)
               end do
!              
!              write out magnetization vectors.
!              xcrysden consider these as X (dummy) atoms.
!              
               do i1=1,nr1,nrws
                 do i2=1,nr2,nrws
                   do i3=1,nr3,nrws
                     normz=gridmx(i1,i2,i3)**2+gridmy(i1,i2,i3)**2+gridmz(i1,i2,i3)**2
                     if(0.1*maxmz <= normz) then
                       xcart2 = matmul (rprimd, (/(i1-one)/nr1, (i2-one)/nr2, (i3-one)/nr3/))
                       write(31,'(A,1X,6(ES17.10,2X))')'X',&
                       Bohr_Ang*(xcart2(1)-shift_tau(1)),&
                       Bohr_Ang*(xcart2(2)-shift_tau(2)),&
                       Bohr_Ang*(xcart2(3)-shift_tau(3)),&
                       gridmx(i1,i2,i3),&
                       gridmy(i1,i2,i3),&
                       gridmz(i1,i2,i3)
                     end if
                   end do
                 end do
               end do
             else
!              ################################################################### (LD)
!              
!              normal case: output density or potential (scalar field)
!              
               write(31,'(1X,A)')  'DIM-GROUP'
               write(31,*) '3  1'
               write(31,'(1X,A)') 'PRIMVEC'
               do i1 = 1,3
                 write(31,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(i2,i1), i2=1,3)
               end do
               write(31,'(1X,A)') 'PRIMCOORD'
               write(31,*) natom, ' 1'
               do iatom = 1,natom
                 write(31,'(i9,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),Bohr_Ang*tau2(1:3,iatom)
               end do
               write(31,'(1X,A)') 'ATOMS'
               do iatom = 1,natom
                 write(31,'(i9,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),Bohr_Ang*tau2(1:3,iatom)
               end do
!              write(31,'(1X,A)') 'FRAMES'
               write(31,'(1X,A)') 'BEGIN_BLOCK_DATAGRID3D'
               write(31,*) 'datagrids'
               write(31,'(1X,A)') 'DATAGRID_3D_DENSITY'
               write(31,*) nr1+1,nr2+1,nr3+1
               write(31,*) '0.0 0.0 0.0 '
               do i1 = 1,3
                 write(31,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(i2,i1), i2=1,3)
               end do

               index = 0
               do ir3=gridshift3+1,nr3+1
                 ii3=mod(ir3-1,nr3) + 1
                 do ir2=gridshift2+1,nr2+1
                   ii2=mod(ir2-1,nr2) + 1
                   do ir1=gridshift1+1,nr1+1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                   do ir1=1,gridshift1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                 end do
                 do ir2=1,gridshift2
                   ii2=mod(ir2-1,nr2) + 1
                   do ir1=gridshift1+1,nr1+1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                   do ir1=1,gridshift1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                 end do
               end do
               do ir3=1,gridshift3
                 ii3=mod(ir3-1,nr3) + 1
                 do ir2=gridshift2+1,nr2+1
                   ii2=mod(ir2-1,nr2) + 1
                   do ir1=gridshift1+1,nr1+1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                   do ir1=1,gridshift1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                 end do
                 do ir2=1,gridshift2
                   ii2=mod(ir2-1,nr2) + 1
                   do ir1=gridshift1+1,nr1+1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                   do ir1=1,gridshift1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                 end do
               end do
               write (31,*)
               write(31,'(1X,A)') 'END_DATAGRID_3D'
               write(31,'(1X,A)') 'END_BLOCK_DATAGRID3D'

             end if

             close(31)
             exit

           case (10) ! formatted for OpenDX
             if (open_file(trim(filnam)//'.dx', message,unit=31,form='formatted') /= 0) then
               MSG_ERROR(message)
             end if
             if (open_file(trim(filnam)//'.xyz',message,unit=32,form='formatted') /= 0) then
               MSG_ERROR(message)
             end if

             write(31, '(a,2x,3i5)' )'object 1 class gridpositions counts',nr3,nr2,nr1
             write(31, '(a)' )'origin 0 0 0'
             write(31,'(a,3(es17.10,2x))')'delta ',(Bohr_Ang*rprimd(i1,3)/nr3, i1=1,3)
             write(31,'(a,3(es17.10,2x))')'delta ',(Bohr_Ang*rprimd(i1,2)/nr2, i1=1,3)
             write(31,'(a,3(es17.10,2x))')'delta ',(Bohr_Ang*rprimd(i1,1)/nr1, i1=1,3)
             write(31, '(a,2x,3i5)' )'object 2 class gridconnections counts',nr3,nr2,nr1
             write(31, '(a)' )'attribute "element type" string "cubes"'
             write(31, '(a)' )'attribute "ref" string "positions"'
             write(31, '(a,1x,i10,1x,a)' )'object 3 class array type float rank 0 items',nr1*nr2*nr3,' data follows'
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   write(31,'(es24.14)') grid(i1,i2,i3)
                 end do
               end do
             end do
             write(31, '(a)' )'attribute "dep" string "positions"'
             write(32, '(i6,/)' ) natom
             do iatom=1,natom
               do ii=1,3
                 xcart2(ii)=xcart(ii,iatom)
                 if (xred(ii,iatom)<-1e-4) then
                   xcart2(ii)=xcart2(ii)+rprimd(ii,ii)
                 else if (xred(ii,iatom)>rprimd(ii,ii)+1e-4) then
                   xcart2(ii)=xcart2(ii)-rprimd(ii,ii)
                 end if
               end do
               write(32,'(i8,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),Bohr_Ang*xcart2(1:3)
             end do
             write(31, '(a)' )'object "density" class field'
             write(31, '(a)' )'component "positions" value 1'
             write(31, '(a)' )'component "connections" value 2'
             write(31, '(a)' )'component "data" value 3'
             close(31)
             close(32)
             exit

           case(11)

             call hirsh(grid,natom,nr1,nr2,nr3,ntypat,rprimd,xcart,hdr%typat,hdr%zionpsp,hdr%znucltypat)
             exit

           case(12)
             write(std_out,*) 'NetCDF is not defined. You must choose another option'
             exit

           case(14) ! CUBE file format from GAUSSIAN

             write(std_out,*)
             write(std_out,*) 'Output a cube file of 3D volumetric data'
             write(std_out,*)

!            EXAMPLE FROM THE WEB
!            CPMD CUBE FILE.
!            OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
!            3    0.000000    0.000000    0.000000
!            40    0.283459    0.000000    0.000000
!            40    0.000000    0.283459    0.000000
!            40    0.000000    0.000000    0.283459
!            8    0.000000    5.570575    5.669178    5.593517
!            1    0.000000    5.562867    5.669178    7.428055
!            1    0.000000    7.340606    5.669178    5.111259
!            -0.25568E-04  0.59213E-05  0.81068E-05  0.10868E-04  0.11313E-04  0.35999E-05

             if (open_file(filnam,message,unit=31,status='unknown',form='formatted') /= 0) then
               MSG_ERROR(message)
             end if

!            %% call print_fofr_cube(nr1,nr2,n3,nr1,nr2,nr3,fofr,rprimd,natom,znucl_atom,xcart,unit=31)
             write(31,'(a)') 'ABINIT generated cube file'
             write(31,'(a)') 'from cut3d tool'

             write(31,'(i9,3(1x,f12.6))') natom,0.,0.,0.
             write(31,'(i9,3(1x,f12.6))') nr1,(rprimd(ir2,1)/nr1, ir2=1,3)
             write(31,'(i9,3(1x,f12.6))') nr2,(rprimd(ir2,2)/nr2, ir2=1,3)
             write(31,'(i9,3(1x,f12.6))') nr3,(rprimd(ir2,3)/nr3, ir2=1,3)

             do iatom=1,natom
               write(31,'(i9,4(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),0.d0, &
&               xcart(1,iatom),xcart(2,iatom),xcart(3,iatom)
             end do

!            C ordering of the indexes
             do i1=1,nr1
               do i2=1,nr2
                 do i3=1,nr3
                   write(31,'(6(f12.6,2x))') grid(i1,i2,i3)
                 end do
               end do
             end do

             close(31)

             exit

           case(0)
             write(std_out,*)' Exit requested by user'
             exit

           case default
             cycle

         end select

       end do

       write(std_out,*) ' Task ',itask,' has been done !'
       write(std_out,*)
       write(std_out,'(a)') ' More analysis of the 3D file ? ( 0=no ; 1=default=yes ; 2= treat another file - restricted usage)'
       read(*,*) iprompt
       if(iprompt/=1) then
         if(densfileformat==1 .or. densfileformat == 2)then
           call hdr_free(hdr)
         end if
         exit
       else
         cycle
       end if
     end do

   end if ! WF file or DEN/POT file

!  A maximum number of files had been previously specified, but set the actual number of files to 1 if one does not read at least one other.
   if(ifiles==1)then
     nfiles=1
     if(iprompt==2)nfiles=mfiles

!    A data structure for storing the important information should be created ...
!    Here, one supposes that the files are compatible ...
     if(isdenpot(ifiles)==1)then
       ABI_MALLOC(grid_full_stored,(nr1,nr2,nr3,nspden,nfiles))
       nr1_stored=nr1
       nr2_stored=nr2
       nr3_stored=nr3
       nspden_stored=nspden
     else if(isdenpot(ifiles)/=1 .and. iprompt==2)then
       message = " in case of storage mode, the first file must be a density/potential file."
       MSG_ERROR(message)
     end if
   end if

   if(isdenpot(ifiles)==1) grid_full_stored(:,:,:,:,ifiles)=grid_full(:,:,:,:)
   if(isdenpot(ifiles)==1) filrho_stored(ifiles)=filrho

   if(allocated(xcart)) then
     ABI_FREE(xcart)
   end if
   if(allocated(xred)) then
     ABI_FREE(xred)
   end if
   if(allocated(grid)) then
     ABI_FREE(grid)
   end if
   if(allocated(grid_full)) then
     ABI_FREE(grid_full)
   end if
   if(allocated(gridtt)) then
     ABI_FREE(gridtt)
   end if
   if(allocated(gridmx)) then
     ABI_FREE(gridmx)
   end if
   if(allocated(gridmy)) then
     ABI_FREE(gridmy)
   end if
   if(allocated(gridmz)) then
     ABI_FREE(gridmz)
   end if
   if(allocated(rhomacu)) then
     ABI_FREE(rhomacu)
   end if
   if(allocated(tau2)) then
     ABI_FREE(tau2)
   end if

   if(iprompt/=2) then
     exit
   end if

 end do ! End big loop on files

!Will provide different information on the density and potential files
 do ifiles=1,nfiles
   if(isdenpot(ifiles)==1)then
     write(std_out,*)
     write(std_out,*) ' Provide some global information about the density and/or potential file(s)'
     exit
   end if
 end do
 do ifiles=1,nfiles
   if(isdenpot(ifiles)==1)then
     write(std_out,*)
     write(std_out, '(a,i5,3a)' ) '  File number ',ifiles,', with name "',trim(filrho_stored(ifiles)),'"'
     write(std_out, '(a,i12,a,es14.6)' ) '  Number of grid points =',nr1*nr2*nr3,' ; Volume of real space cell (Bohr^3)=',ucvol
     do ispden=1,nspden
       sumdenpot=sum(grid_full_stored(:,:,:,ispden,ifiles))
       write(std_out, '(a,i5,3a)' ) '   Spin-component number ',ispden
       write(std_out, '(a,3es16.6)' ) '      Sum of values, mean, mean times cell volume=',&
&       sumdenpot,sumdenpot/real(nr1*nr2*nr3),sumdenpot*ucvol/real(nr1*nr2*nr3)
     end do
   end if
 end do

 if(nspden==1)then
!  At present, only nspden=1 is correctly implemented, due to specificities of the treatment of the spin-density
   do ifiles=1,nfiles
     if(isdenpot(ifiles)==1)then
       write(std_out,*)
       write(std_out,'(a)') ' Provide some global joint information about the stored density and potential file(s)'
       exit
     end if
   end do
   do ifiles=1,nfiles
     if(isdenpot(ifiles)==1)then
       do jfiles=ifiles,nfiles
         if(isdenpot(jfiles)==1)then
           write(std_out,*)
           write(std_out, '(a,2i5)' )'  File numbers :',ifiles,jfiles
           do ispden=1,nspden
             dotdenpot=zero
             do ir1=1,nr1
               do ir2=1,nr2
                 do ir3=1,nr3
                   dotdenpot=dotdenpot+grid_full_stored(ir1,ir2,ir3,ispden,ifiles)*grid_full_stored(ir1,ir2,ir3,ispden,jfiles)
                 end do
               end do
             end do
             write(std_out, '(a,i5,3a)' ) '   Spin-component number ',ispden
             write(std_out, '(a,3es16.6)' ) '      Dot product of values, mean, mean times cell volume=',&
!            write(std_out, '(a,3es20.10)' ) '      Dot product of values, mean, mean times cell volume=',&
&             dotdenpot,dotdenpot/real(nr1*nr2*nr3),dotdenpot*ucvol/real(nr1*nr2*nr3)
           end do
         end if
       end do
     end if
   end do
 end if

 if(allocated(grid_full_stored)) then
   ABI_FREE(grid_full_stored)
 end if

 if(allocated(isdenpot)) then
   ABI_FREE(isdenpot)
 end if

 call timein(tsec(1),tsec(2))
 tsec(1)=tsec(1)-tcpui
 tsec(2)=tsec(2)-twalli

 write(std_out, '(3a,f13.1,a,f13.1)' )'-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

 write(std_out,*)
 write(std_out,*) ' Thank you for using me'
 write(std_out,*)

 call flush_unit(std_out)

 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program cut3d
!!***
