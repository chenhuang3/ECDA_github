!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_paw_free
!! NAME
!!  wvl_paw_free
!!
!! FUNCTION
!!  Frees memory for WVL+PAW implementation
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ntypat = number of atom types
!!  wvl= wvl type
!!  wvl_proj= wvl projector type
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_paw_free(ntypat,wvl,wvl_proj)
    
 use defs_basis
 use defs_wvltypes
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_paw_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)::ntypat
 type(wvl_internal_type),intent(inout)::wvl
 type(wvl_projectors_type),intent(inout)::wvl_proj
  
!Local variables-------------------------------
 integer:: itypat
!character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
#if defined HAVE_DFT_BIGDFT
 if( associated(wvl%paw%spsi)) then
   ABI_DEALLOCATE(wvl%paw%spsi)
 end if
 if( associated(wvl%paw%indlmn)) then
   ABI_DEALLOCATE(wvl%paw%indlmn)
 end if
 if( associated(wvl%paw%sij)) then
   ABI_DEALLOCATE(wvl%paw%sij)
 end if
 if( associated(wvl%paw%rpaw)) then
   ABI_DEALLOCATE(wvl%paw%rpaw)
 end if


!proj_G
 do itypat=1,ntypat
   if( associated(wvl_proj%G(itypat)%ndoc)) then
     ABI_DEALLOCATE(wvl_proj%G(itypat)%ndoc)
   end if
   if( associated(wvl_proj%G(itypat)%nam)) then
     ABI_DEALLOCATE(wvl_proj%G(itypat)%nam)
   end if
   if( associated(wvl_proj%G(itypat)%xp)) then
     ABI_DEALLOCATE(wvl_proj%G(itypat)%xp)
   end if
   if( associated(wvl_proj%G(itypat)%psiat)) then
     ABI_DEALLOCATE(wvl_proj%G(itypat)%psiat)
   end if
 end do
 if( allocated(wvl_proj%G)) then
   ABI_DATATYPE_DEALLOCATE(wvl_proj%G)
 end if

!rholoc
 if( associated(wvl%rholoc%msz )) then
   ABI_DEALLOCATE(wvl%rholoc%msz)
 end if
 if( associated(wvl%rholoc%d )) then
   ABI_DEALLOCATE(wvl%rholoc%d)
 end if
 if( associated(wvl%rholoc%rad)) then
   ABI_DEALLOCATE(wvl%rholoc%rad)
 end if
 if( associated(wvl%rholoc%radius)) then
   ABI_DEALLOCATE(wvl%rholoc%radius)
 end if

#endif
!
!paw%paw_ij and paw%cprj are allocated and deallocated inside vtorho

end subroutine wvl_paw_free
!!***
