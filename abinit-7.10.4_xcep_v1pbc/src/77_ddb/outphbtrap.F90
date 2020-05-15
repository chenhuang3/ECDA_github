!{\src2tex{textfont=tt}}
!!****f* ABINIT/outphbtrap
!! NAME
!! outphbtrap
!!
!! FUNCTION
!!  Print out phonon frequencies on regular grid for BoltzTrap
!!  Flag in input file is outboltztrap=1
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2014 ABINIT group (MVer)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Ifc<ifc_type>=Stores data related to interatomic force constants.
!!  Crystal<crystal_t>=Info on the crystal structure
!!  anaddb_dtset= (derived datatype) contains all the input variables
!!  basename = file name for output to disk
!!
!! OUTPUT
!!  only write to file
!!
!! TODO
!!  anaddb_dtset can be removed, params are now passed via Ifc
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      appdig,ifc_fourq,smpbz,symkpt,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outphbtrap(Ifc,Crystal,anaddb_dtset,basename)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_io_tools,  only : open_file
 use m_dynmat,    only : gtdyn9
 use m_crystal,   only : crystal_t
 use m_ifc,       only : ifc_type, ifc_fourq
 use m_anaddb_dataset, only : anaddb_dataset_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outphbtrap'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 character(len=*),intent(in) :: basename
 type(ifc_type),intent(in) :: Ifc
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
 type(crystal_t),intent(in) :: Crystal

!Local variables -------------------------
!scalars
 integer,parameter :: brav1=1,chksymbreak0=0
 integer :: natom,nsym
 integer :: facbrv,imode,iq_ibz,msym
 integer :: nqbz,nqpt_max,nqshft,option,timrev
 integer :: nqibz, nreals,unit_btrap
 integer :: iatom, idir
 real(dp) :: bzvol
 character(len=500) :: message
 character(len=500) :: format_nreals,format_line_btrap
 character(len=fnlen) :: outfile
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: ibz2bz(:)
 real(dp) :: d2cart(2,3,Crystal%natom,3,Crystal%natom),displ(2*3*Crystal%natom*3*Crystal%natom)
 real(dp) :: gprimd(3,3),phfrq(3*Crystal%natom)
 real(dp) :: qphon(3)
 real(dp),allocatable :: qbz(:,:),qibz(:,:),qshft(:,:)
 real(dp),allocatable :: wtq(:),wtq_folded(:),wtqibz(:)

! *********************************************************************

 DBG_ENTER("COLL")

 natom = Crystal%natom
 nsym  = Crystal%nsym
 msym  = nsym

 outfile = trim(basename) // '_BTRAP'
 write(message, '(3a)')ch10,&
& ' Will write phonon FREQS in BoltzTrap format to file ',trim(outfile)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if (open_file(outfile,message,newunit=unit_btrap,status="replace") /= 0) then
   MSG_ERROR(message)
 end if

 write (unit_btrap,'(a)') '#'
 write (unit_btrap,'(a)') '# ABINIT package : Boltztrap phonon file. Remove this header before feeding to BT'
 write (unit_btrap,'(a)') '#    for compatibility with PHON output the freq are in Ry (before the square)'
 write (unit_btrap,'(a)') '#'
 write (unit_btrap,'(a)') '#    nq, nband  '
 write (unit_btrap,'(a)') '#  qx, qy, qz   '
 write (unit_btrap,'(a)') '#  qpt weight   '
 write (unit_btrap,'(a)') '#  freq_1^2, dynmat column for mode 1 '
 write (unit_btrap,'(a)') '#  etc for mode 2,3,4... qpt 2,3,4... '

 gprimd = Crystal%gprimd
 bzvol=ABS ( gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
& -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
& +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))
!
!Save memory during the generation of the q-mesh in the full BZ  
!Take into account the type of Bravais lattice
 facbrv=1
 if (brav1==2) facbrv=2
 if (brav1==3) facbrv=4

 nqshft=1 !always 1 
 ABI_ALLOCATE(qshft,(3,nqshft))
 qshft(:,1)=anaddb_dtset%q2shft(:)  ! FIXME small inconsistency in the dimension of q1shft

 nqpt_max=(anaddb_dtset%ng2qpt(1)*anaddb_dtset%ng2qpt(2)*anaddb_dtset%ng2qpt(3)*nqshft)/facbrv
 ABI_ALLOCATE(qibz,(3,nqpt_max))
 ABI_ALLOCATE(qbz,(3,nqpt_max))

 qptrlatt(:,:)=0
 qptrlatt(1,1)=anaddb_dtset%ng2qpt(1)
 qptrlatt(2,2)=anaddb_dtset%ng2qpt(2)
 qptrlatt(3,3)=anaddb_dtset%ng2qpt(3)
 option=1 
!
!here I noticed a problem in the declaration of q1shft in the anaddb datatype 
!FIXME we write on unit std_out just to avoid problem with automatic tests
 call smpbz(brav1,std_out,qptrlatt,nqpt_max,nqbz,nqshft,option,qshft,qbz)
!
!Reduce the number of such points by symmetrization.
 ABI_ALLOCATE(ibz2bz,(nqbz))
 ABI_ALLOCATE(wtq,(nqbz))
 ABI_ALLOCATE(wtq_folded,(nqbz))
 wtq(:)=one/nqbz         ! Weights sum up to one
 timrev=1; option=1     ! TODO timrev should be input 
!
 call symkpt(chksymbreak0,Crystal%gmet,ibz2bz,std_out,qbz,nqbz,nqibz,nsym,Crystal%symrec,timrev,wtq,wtq_folded)
 write(std_out,*) 'nqibz = ', nqibz

 ABI_ALLOCATE(wtqibz,(nqibz))
 do iq_ibz=1,nqibz
   wtqibz(iq_ibz)=wtq_folded(ibz2bz(iq_ibz))
   qibz(:,iq_ibz)=qbz(:,ibz2bz(iq_ibz))
 end do
 ABI_DEALLOCATE(wtq_folded)
 ABI_DEALLOCATE(qshft)

 write (unit_btrap,'(2I6)') nqibz, 3*natom 

! Loop over irreducible q-points
 do iq_ibz=1,nqibz
   qphon(:)=qibz(:,iq_ibz);

   call ifc_fourq(Ifc,Crystal,qphon,phfrq,displ,out_d2cart=d2cart)

   write (unit_btrap,'(3E20.10)') qphon
   write (unit_btrap,'(E20.10)') wtqibz(iq_ibz)
   nreals=1+2*3*natom
   call appdig(nreals,'(',format_nreals)
   format_line_btrap=trim(format_nreals)//'E20.10)'
   do iatom = 1, natom
     do idir = 1, 3
       imode = idir + 3*(iatom-1)
!      factor two for Ry output - this may change in definitive BT and abinit formats 
       write (unit_btrap,trim(format_line_btrap))phfrq(imode)*two,d2cart(1:2,1:3,1:natom,idir,iatom)
!      XG130409 : This was the old coding for the previous line, but was not portable...
!      write (unit_btrap,'(E20.10)', ADVANCE='NO') phfrq(imode)*two
!      do jatom = 1, natom
!      do jdir = 1, 3
!      write (unit_btrap,'(2E20.10)', ADVANCE='NO') d2cart(:,jdir,jatom,idir, iatom)
!      end do
!      end do
!      write (unit_btrap,*)
     end do
   end do
   
 end do !irred q-points

 ABI_DEALLOCATE(ibz2bz)
 ABI_DEALLOCATE(qibz)
 ABI_DEALLOCATE(qbz)
 ABI_DEALLOCATE(wtq)
 ABI_DEALLOCATE(wtqibz)

 close (unit=unit_btrap)

 DBG_EXIT("COLL")

end subroutine outphbtrap
!!***
