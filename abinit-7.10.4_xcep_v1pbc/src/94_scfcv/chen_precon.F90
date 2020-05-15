!{\src2tex{textfont=tt}}
!!****f* ABINIT/precon
!!
!! NAME
!! precon
!!
!! FUNCTION
!! precondition $<G|(H-e_{n,k})|C_{n,k}>$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2011 ABINIT group (DCA, XG, GMR, MT))
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  $cg(2,npw)=<G|C_{n,k}>$.
!!  $eval=current band eigenvalue=<C_{n,k}|H|C_{n,k}>$.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  nspinor=number of spinorial components of the wavefunctions
!!  $vect(2,npw)=<G|H|C_{n,k}>$.
!!  npw=number of planewaves at this k point.
!!  optekin= 1 if the kinetic energy used in preconditionning is modified
!!             according to Kresse, Furthmuller, PRB 54, 11169 (1996)
!!           0 otherwise
!!
!!
!! OUTPUT
!!  vect(2,npw*nspinor)=<G|(H-eval)|C_{n,k}>*(polynomial ratio)
!!
!! PARENTS
!!      cgwf,cgwf3
!!
!! CHILDREN
!!      timab,wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine chen_precon(max_ke,kinpw,npw,nspinor,optekin,vect)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: npw,nspinor,optekin
!arrays
 real(dp),intent(in) :: kinpw(npw), max_ke
 real(dp),intent(inout) :: vect(2,npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ierr,ig,igs,ipw1,ispinor,old_paral_level,spaceComm
 real(dp) :: ek0,ek0_inv,fac,poly,xx
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

!DEBUG
!write(6,*)' precon: debug, enter.'
!ENDDEBUG

! get the max of the kinetic energy 
 ek0 = max_ke 

 if(ek0<1.0d-10)then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' precon : WARNING -',ch10,&
&   '  The mean kinetic energy of a wavefunction vanishes.',ch10,&
&   '  It is reset to 0.1Ha.'
   call wrtout(std_out,message,'PERS')
   ek0=0.1_dp
 end if

 if (optekin==1) then
   ek0_inv=2.0_dp/(3._dp*ek0)
 else
   ek0_inv=1.0_dp/ek0
 end if
!
!Carry out preconditioning
 do ispinor=1,nspinor
   igs=(ispinor-1)*npw
   do ig=1+igs,npw+igs
     if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
       xx=kinpw(ig-igs)*ek0_inv
!      Teter polynomial ratio
       poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
       fac=poly/(poly+16._dp*xx**4)
       if (optekin==1) fac=two*fac
       vect(1,ig)=( vect(1,ig) )*fac
       vect(2,ig)=( vect(2,ig) )*fac
     else
       vect(1,ig)=zero
       vect(2,ig)=zero
     end if
   end do
!  $OMP END PARALLEL DO
 end do

!DEBUG
!write(6,*)' precon: debug, exit.'
!ENDDEBUG

end subroutine chen_precon
!!***
