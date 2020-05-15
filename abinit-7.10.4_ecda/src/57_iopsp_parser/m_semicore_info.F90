!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_semicore_info
!! NAME
!!   m_semicore_info
!!
!! FUNCTION
!!   module for handling semicore psp information in XML format shared with SIESTA
!!
!! COPYRIGHT
!! Copyright (C) 2014 ABINIT group (AGarcia)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_semicore_info

! Routines to determine whether semicore states
! have been implicitly considered in the generation
! of a pseudopotential.
! This is relevant for the generation of the appropriate
! number of KB projectors in the non-local form of the
! pseudopotential.
!
! This version uses the new PSML library syntax, which is
! still in flux. To use the library, define HAVE_PSML_SUPPORT.

! Wrappers are provided if the new library is not available, and one has
! to resort to the old interface with explicit handling of the derived
! type "pseudo_t".
!
! Alberto Garcia, June 4, 2014
!
!#ifdef HAVE_PSML_SUPPORT
#if 0
use m_psml, only: psml_t
use m_psml, only: psxmlAtomicNumber
use m_psml, only: psxmlPotentialsDown
use m_psml, only: psxmlPrincipalN
use m_psml, only: psxmlPotAngMomentum
#else
use m_xml_pseudo_types, only: psml_t => pseudo_t
#endif /* HAVE_PSML_SUPPORT */

use defs_basis

public :: get_n_semicore_shells

CONTAINS
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_n_generation
!! NAME
!!   get_n_generation
!!
!! FUNCTION
!!   
!!
!! COPYRIGHT
!! Copyright (C) 2014 ABINIT group (AGarcia).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_semicore_info
!!
!! CHILDREN
!!      cnfig,get_n_generation
!!
!! SOURCE
  subroutine get_n_generation(p,nval)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_n_generation'
!End of the abilint section

    type(psml_t), intent(in) :: p
    integer, intent(out) :: nval(0:3)
  !
  ! Returns an array with the n-quantum numbers of the states used
  ! to generate the pseudopotential on channel l (0..3)
  ! If no pseudo is available for a given l, a -1 is returned.

    integer :: i, l, npots

    nval(:) = -1
    npots = psxmlPotentialsDown(p)
    do i = 1, npots
       l = psxmlPotAngMomentum(p,"d",i)
       nval(l) = psxmlPrincipalN(p,"d",i)
    enddo
  end subroutine get_n_generation
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_n_semicore_shells
!! NAME
!!   get_n_semicore_shells
!!
!! FUNCTION
!!   
!!
!! COPYRIGHT
!! Copyright (C) 2014 ABINIT group (AGarcia).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      inpspheads,psp9in
!!
!! CHILDREN
!!      cnfig,get_n_generation
!!
!! SOURCE
  subroutine get_n_semicore_shells(p,nsemic)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_n_semicore_shells'
!End of the abilint section

    type(psml_t), intent(in) :: p
    integer, intent(out)       :: nsemic(0:3)
  !
  ! Returns an array with the number of semicore shells
  ! on channel l (0..3)
  ! If no pseudo is available for a given l, a 0 is returned.
  ! If for some reason the generation step did not pseudize
  ! a proper valence state (e.g. Cu 3d), the number of semicore
  ! shells is returned as -1. This should be checked by the caller.

    integer :: l, z
    integer :: nval(0:3)
    integer :: nval_gs(0:3)

    z = int(psxmlAtomicNumber(p))
    call cnfig(z,nval_gs)
    call get_n_generation(p,nval)
    nsemic(:) = 0
    do l = 0, 3
!!     DEBUG
!       write(std_out,*)' l, nval_gs(l), nval(l) = ', l, nval_gs(l), nval(l)
!!     END DEBUG
       if (nval(l) > 0) then
          nsemic(l) = nval_gs(l) - nval(l)
       endif
    enddo
  end subroutine get_n_semicore_shells
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/CNFIG
!! NAME
!!   CNFIG
!!
!! FUNCTION
!!   
!!
!! COPYRIGHT
!! Copyright (C) 2014 ABINIT group (AGarcia).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!
! This routine encodes some choices regarding the core-valence split,
! which might not be universal.
! 
      SUBROUTINE CNFIG( Z, NVAL ) 
! Returns the valence configuration for atomic ground state, i.e.
! the principal quantum number NVAL of the valence orbilas for each L
! Originally written by A.R.Williams. Modified by J.M.Soler

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'CNFIG'
!End of the abilint section

      integer,intent(in) :: Z        ! Atomic number
      integer,intent(out):: NVAL(0:3) ! Valence electrons for each L

      integer, parameter :: LMAX=3, NCHNG=15
      integer :: ICHNG, L, LCHNG(NCHNG), ZCHNG(NCHNG)

      ! Originally: s valence orbital switched for p occupation = 4
      !           Li, F,Na,Cl, K,Ga,Br,Rb,In, I,Cs,Hf,Tl,At,Fr
!!     DATA ZCHNG / 3, 9,11,17,19,31,35,37,49,53,55,72,81,85,87/

      ! Changed to: s valence orbital switched for full p occupation
      !           Li,Na,Na, K, K,Ga,Rb,Rb,In,Cs,Cs,Hf,Tl,Fr,Fr
      DATA ZCHNG / 3,11,11,19,19,31,37,37,49,55,55,72,81,87,87/
      DATA LCHNG / 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/
      DO L=0,LMAX
         NVAL(L)=L+1
      END DO
      DO ICHNG=1,NCHNG
         IF (ZCHNG(ICHNG).GT.Z) EXIT
         L=LCHNG(ICHNG)
         NVAL(L)=NVAL(L)+1
      END DO

      END subroutine cnfig
!!***

!#ifndef HAVE_PSML_SUPPORT
#if 1


!
!  Provide emulators for the psml functions
!
!{\src2tex{textfont=tt}}
!!****f* ABINIT/psxmlAtomicNumber
!! NAME
!!   psxmlAtomicNumber
!!
!! FUNCTION
!!   
!!
!! COPYRIGHT
!! Copyright (C) 2014 ABINIT group (AGarcia).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
function psxmlAtomicNumber(psxml) result(z)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psxmlAtomicNumber'
!End of the abilint section

type(psml_t), intent(in) :: psxml
real(dp) :: z
 z = psxml%header%atomicnumber
end function psxmlAtomicNumber
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/psxmlPotentialsUp
!! NAME
!!   psxmlPotentialsUp
!!
!! FUNCTION
!!   
!!
!! COPYRIGHT
!! Copyright (C) 2014 ABINIT group (AGarcia).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
function psxmlPotentialsUp(psxml) result(n)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psxmlPotentialsUp'
!End of the abilint section

type(psml_t), intent(in) :: psxml
integer                    :: n
n = psxml%npots_up
end function psxmlPotentialsUp
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/psxmlPotentialsDown
!! NAME
!!   psxmlPotentialsDown
!!
!! FUNCTION
!!   
!!
!! COPYRIGHT
!! Copyright (C) 2014 ABINIT group (AGarcia).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
function psxmlPotentialsDown(psxml) result(n)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psxmlPotentialsDown'
!End of the abilint section

type(psml_t), intent(in) :: psxml
integer                    :: n
n = psxml%npots_down
end function psxmlPotentialsDown
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/psxmlPotAngMomentum
!! NAME
!!   psxmlPotAngMomentum
!!
!! FUNCTION
!!   
!!
!! COPYRIGHT
!! Copyright (C) 2014 ABINIT group (AGarcia).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
function psxmlPotAngMomentum(psxml,ud,i) result(l)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psxmlPotAngMomentum'
!End of the abilint section

type(psml_t), intent(in) :: psxml
character, intent(in)      :: ud
integer,   intent(in)      :: i
integer                    :: l

character(len=1), dimension(0:4) :: sym = (/ "s", "p", "d", "f", "g" /)
integer :: ndown
character(len=1) :: str

ndown = psxmlPotentialsDown(psxml)

select case(ud)
   case ( "u", "U")
      if (i > psxmlPotentialsUp(psxml)) then
!         call die("attempt to get l from non-existing Up potential")
      endif
      str = psxml%pot(ndown+i)%l
   case ( "d", "D")
      if (i > ndown) then
!         call die("attempt to get l from non-existing Down potential")
      endif
      str = psxml%pot(i)%l
end select
!
do l = 0,4
   if (str == sym(l)) RETURN
enddo
!call die("Wrong l symbol in potential")

end function psxmlPotAngMomentum
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/psxmlPrincipalN
!! NAME
!!   psxmlPrincipalN
!!
!! FUNCTION
!!   
!!
!! COPYRIGHT
!! Copyright (C) 2014 ABINIT group (AGarcia).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
function psxmlPrincipalN(psxml,ud,i) result(n)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psxmlPrincipalN'
!End of the abilint section

type(psml_t), intent(in) :: psxml
character, intent(in)      :: ud
integer,   intent(in)      :: i
integer                    :: n

integer :: ndown

ndown = psxmlPotentialsDown(psxml)

select case(ud)
   case ( "u", "U")
      if (i > psxmlPotentialsUp(psxml)) then
!         call die("attempt to get n from non-existing Up potential")
      endif
      n = psxml%pot(ndown+i)%n
   case ( "d", "D")
      if (i > ndown) then
!         call die("attempt to get n from non-existing Down potential")
      endif
      n = psxml%pot(i)%n
end select
end function psxmlPrincipalN
!!***
#endif       /*HAVE_PSML_SUPPORT*/

end module m_semicore_info
!!***
