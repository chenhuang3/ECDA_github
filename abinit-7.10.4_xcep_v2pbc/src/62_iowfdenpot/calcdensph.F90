!{\src2tex{textfont=tt}}
!!****f* ABINIT/calcdensph
!! NAME
!! calcdensph
!!
!! FUNCTION
!! Compute and print integral of total density inside spheres around atoms.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (MT,ILuk,MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of atom types
!!  nunit=number of the unit for writing
!!  ratsph(ntypat)=radius of spheres around atoms
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type of each atom
!!  ucvol=unit cell volume in bohr**3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  intgden(nspden, natom)=intgrated density (magnetization...) for each atom in a sphere of radius ratsph. Optional arg
!!  Rest is printing
!!
!! PARENTS
!!      mag_constr,mag_constr_e,outscfcv,vtorho
!!
!! CHILDREN
!!      timab,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,nunit,ratsph,rhor,rprimd,typat,ucvol,xred, &
&    intgden)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_xmpi, only : xmpi_sum

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calcdensph'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,nunit
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: gmet(3,3),ratsph(ntypat),rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out),optional :: intgden(nspden,natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: ishift=5
 integer :: i1,i2,i3,iatom,ierr,ifft_local,ix,iy,iz,izloc,n1,n1a,n1b,n2
 integer :: n2a,n2b,n3,n3a,n3b,nd3,nfftot
 real(dp),parameter :: delta=0.99_dp
 real(dp) :: difx,dify,difz,r2,r2atsph,rr1,rr2,rr3,rx,ry,rz
 real(dp) :: fsm, ratsm, ratsm2
 logical   :: grid_found
 character(len=500) :: message
!arrays
 real(dp) :: intgden_(nspden,natom)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: tsec(2)

! *************************************************************************

 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3);nd3=n3/mpi_enreg%nproc_fft
 nfftot=n1*n2*n3
 intgden_=zero

 ratsm = zero
 if (present(intgden)) then
   ratsm = 0.05_dp ! default value for the smearing region radius - may become input variable later
 end if

!Get the distrib associated with this fft_grid
 grid_found=.false.
 if (n2 == mpi_enreg%distribfft%n2_coarse ) then
   if (n3== size(mpi_enreg%distribfft%tab_fftdp3_distrib)) then
     fftn3_distrib => mpi_enreg%distribfft%tab_fftdp3_distrib
     ffti3_local => mpi_enreg%distribfft%tab_fftdp3_local
     grid_found=.true.
   end if
 end if
 if (n2 == mpi_enreg%distribfft%n2_fine ) then
   if (n3 == size(mpi_enreg%distribfft%tab_fftdp3dg_distrib)) then
     fftn3_distrib => mpi_enreg%distribfft%tab_fftdp3dg_distrib
     ffti3_local => mpi_enreg%distribfft%tab_fftdp3dg_local
     grid_found = .true.
   end if
 end if
 if (.not.(grid_found)) then
   MSG_BUG("Unable to find an allocated distrib for this fft grid")
 end if

!Loop over atoms
!-------------------------------------------
 do iatom=1,natom

!  Define a "box" around the atom
   r2atsph=1.0000001_dp*ratsph(typat(iatom))**2
   rr1=sqrt(r2atsph*gmet(1,1))
   rr2=sqrt(r2atsph*gmet(2,2))
   rr3=sqrt(r2atsph*gmet(3,3))

   n1a=int((xred(1,iatom)-rr1+ishift)*n1+delta)-ishift*n1
   n1b=int((xred(1,iatom)+rr1+ishift)*n1      )-ishift*n1
   n2a=int((xred(2,iatom)-rr2+ishift)*n2+delta)-ishift*n2
   n2b=int((xred(2,iatom)+rr2+ishift)*n2      )-ishift*n2
   n3a=int((xred(3,iatom)-rr3+ishift)*n3+delta)-ishift*n3
   n3b=int((xred(3,iatom)+rr3+ishift)*n3      )-ishift*n3

   ratsm2 = (2*ratsph(typat(iatom))-ratsm)*ratsm 

   do i3=n3a,n3b
     iz=mod(i3+ishift*n3,n3)
     if(fftn3_distrib(iz+1)==mpi_enreg%me_fft) then
       izloc = ffti3_local(iz+1) - 1
       difz=dble(i3)/dble(n3)-xred(3,iatom)
       do i2=n2a,n2b
         iy=mod(i2+ishift*n2,n2)
         dify=dble(i2)/dble(n2)-xred(2,iatom)
         do i1=n1a,n1b
           ix=mod(i1+ishift*n1,n1)
           difx=dble(i1)/dble(n1)-xred(1,iatom)
           rx=difx*rprimd(1,1)+dify*rprimd(1,2)+difz*rprimd(1,3)
           ry=difx*rprimd(2,1)+dify*rprimd(2,2)+difz*rprimd(2,3)
           rz=difx*rprimd(3,1)+dify*rprimd(3,2)+difz*rprimd(3,3)
           r2=rx**2+ry**2+rz**2


!          Identify the fft indexes of the rectangular grid around the atom
           if (r2 > r2atsph) cycle

           fsm = radsmear(r2, r2atsph, ratsm2)

           ifft_local=1+ix+n1*(iy+n2*izloc)
           if (nspden==1) then
!            intgden_(1,iatom)= integral of total density
             intgden_(1,iatom)=intgden_(1,iatom)+fsm*rhor(ifft_local,1)
           else if (nspden==2) then
!            intgden_(1,iatom)= integral of up density
!            intgden_(1,iatom)= integral of dn density
             intgden_(1,iatom)=intgden_(1,iatom)+fsm*rhor(ifft_local,2)
             intgden_(2,iatom)=intgden_(2,iatom)+fsm*rhor(ifft_local,1)-rhor(ifft_local,2)
           else
!            intgden_(1,iatom)= integral of total density
!            intgden_(2,iatom)= integral of magnetization, x-component
!            intgden_(3,iatom)= integral of magnetization, y-component
!            intgden_(4,iatom)= integral of magnetization, z-component
             intgden_(1,iatom)=intgden_(1,iatom)+fsm*rhor(ifft_local,1)
             intgden_(2,iatom)=intgden_(2,iatom)+fsm*rhor(ifft_local,2)
             intgden_(3,iatom)=intgden_(3,iatom)+fsm*rhor(ifft_local,3)
             intgden_(4,iatom)=intgden_(4,iatom)+fsm*rhor(ifft_local,4)
           end if

         end do
       end do
     end if
   end do

   intgden_(:,iatom)=intgden_(:,iatom)*ucvol/dble(nfftot)

!  End loop over atoms
!  -------------------------------------------
 end do

!MPI parallelization
 if(mpi_enreg%nproc_fft>1)then
   call timab(48,1,tsec)
   call xmpi_sum(intgden_,mpi_enreg%comm_fft,ierr)
   call timab(48,2,tsec)
 end if

!Printing
 write(message, '(4a)' ) ch10,&
& ' Integrated total density in atomic spheres:',ch10,&
& ' -------------------------------------------'
 call wrtout(nunit,message,'COLL')
 if (nspden==1) then
   write(message, '(a)' ) ' Atom  Sphere_radius  Integrated_density'
   call wrtout(nunit,message,'COLL')
   do iatom=1,natom
     write(message, '(i5,f15.5,f20.8)' ) iatom,ratsph(typat(iatom)),intgden_(1,iatom)
     call wrtout(nunit,message,'COLL')
   end do
 else if(nspden==2) then
   write(message, '(a)' ) ' Atom  Sphere radius  Integrated_up_density  Integrated_dn_density  Total(up+dn)   Diff(up-dn)'
   call wrtout(nunit,message,'COLL')
   do iatom=1,natom
     write(message, '(i5,f15.5,2f23.8,2f14.8)' ) iatom,ratsph(typat(iatom)),intgden_(1,iatom),intgden_(2,iatom),&
&     (intgden_(1,iatom)+intgden_(2,iatom)),(intgden_(1,iatom)-intgden_(2,iatom))
     call wrtout(nunit,message,'COLL')
   end do
   write(message, '(3a)' ) ' Note: Diff(up-dn) can be considered as a rough ',&
&   'approximation of a local magnetic moment.'
   call wrtout(nunit,message,'COLL')
 else if(nspden==4) then
   write(message, '(a)' ) ' Atom  Sphere radius  Total_density Integrated_x_magnetiz Integrated_y_magnetiz Integrated_z_magnetiz'
   call wrtout(nunit,message,'COLL')
   do iatom=1,natom
     write(message, '(i5,f14.5,f14.8,3f20.8)' ) iatom,ratsph(typat(iatom)),(intgden_(ix,iatom),ix=1,4)
     call wrtout(nunit,message,'COLL')
   end do
 end if
 write(message, '(a)' ) ch10
 call wrtout(nunit,message,'COLL')

 if (present(intgden)) then
   intgden = intgden_
 end if

end subroutine calcdensph
!!***
