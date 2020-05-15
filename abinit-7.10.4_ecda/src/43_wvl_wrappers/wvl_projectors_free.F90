
!!****f* ABINIT/wvl_projectors_free
!!
!! NAME
!! wvl_projectors_free
!!
!! FUNCTION
!! Freeing routine.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  proj <type(wvl_projectors_type)>=projectors informations in a wavelet basis.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      deallocate_proj_descr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_projectors_free(proj)

 use defs_basis
 use defs_wvltypes
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_projectors_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wvl_projectors_type),intent(inout) :: proj
!Local variables -------------------------
#ifndef HAVE_DFT_BIGDFT
 character(len=500) :: message
#endif
  ! *********************************************************************

#if defined HAVE_DFT_BIGDFT
 call deallocate_proj_descr(proj%nlpspd, ABI_FUNC)

#else
 BIGDFT_NOTENABLED_ERROR()
#endif

 ABI_DEALLOCATE(proj%proj)

end subroutine wvl_projectors_free
!!***
