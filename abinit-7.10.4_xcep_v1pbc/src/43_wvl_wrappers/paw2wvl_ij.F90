!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw2wvl_ij
!! NAME
!!  paw2wvl_ij
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop,vtorho
!!
!! CHILDREN
!!      nullify_paw_ij_objects
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw2wvl_ij(my_natom,option,paw_ij,wvl)
    
 use defs_basis
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_paw_ij, only : paw_ij_type

#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only : nullify_paw_ij_objects
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw2wvl_ij'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)::my_natom,option
 type(wvl_internal_type), intent(inout)::wvl
 type(paw_ij_type),intent(in) :: paw_ij(my_natom)
!Local variables-------------------------------
 integer :: iatom,iaux                                     ! to be filled, if needed
 character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
 DBG_ENTER("COLL")

#if defined HAVE_DFT_BIGDFT
!Option==1: allocate and copy
 if(option==1) then
   ABI_DATATYPE_ALLOCATE(wvl%paw%paw_ij,(my_natom))
   do iatom=1,my_natom
     call nullify_paw_ij_objects(wvl%paw%paw_ij(iatom))
   end do
   do iatom=1,my_natom
     iaux=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
     ABI_ALLOCATE(wvl%paw%paw_ij(iatom)%dij,(iaux,paw_ij(iatom)%ndij))
!    
     wvl%paw%paw_ij(iatom)%cplex          =paw_ij(iatom)%cplex
     wvl%paw%paw_ij(iatom)%cplex_dij      =paw_ij(iatom)%cplex_dij
     wvl%paw%paw_ij(iatom)%has_dij        =paw_ij(iatom)%has_dij
     wvl%paw%paw_ij(iatom)%has_dijfr      =paw_ij(iatom)%has_dijfr
     wvl%paw%paw_ij(iatom)%has_dijhartree =paw_ij(iatom)%has_dijhartree
     wvl%paw%paw_ij(iatom)%has_dijhat     =paw_ij(iatom)%has_dijhat
     wvl%paw%paw_ij(iatom)%has_dijso      =paw_ij(iatom)%has_dijso
     wvl%paw%paw_ij(iatom)%has_dijU       =paw_ij(iatom)%has_dijU
     wvl%paw%paw_ij(iatom)%has_dijxc      =paw_ij(iatom)%has_dijxc
     wvl%paw%paw_ij(iatom)%has_dijxc_val  =paw_ij(iatom)%has_dijxc_val
     wvl%paw%paw_ij(iatom)%has_exexch_pot =paw_ij(iatom)%has_exexch_pot
     wvl%paw%paw_ij(iatom)%has_pawu_occ   =paw_ij(iatom)%has_pawu_occ
     wvl%paw%paw_ij(iatom)%lmn_size       =paw_ij(iatom)%lmn_size
     wvl%paw%paw_ij(iatom)%lmn2_size      =paw_ij(iatom)%lmn2_size
     wvl%paw%paw_ij(iatom)%ndij           =paw_ij(iatom)%ndij
     wvl%paw%paw_ij(iatom)%nspden         =paw_ij(iatom)%nspden
     wvl%paw%paw_ij(iatom)%nsppol         =paw_ij(iatom)%nsppol
     wvl%paw%paw_ij(iatom)%dij(:,:)             =paw_ij(iatom)%dij(:,:)

!    debug,
!    write(*,*)'paw2wvl_ij, erase me, set dij=0'
!    wvl%paw%paw_ij(iatom)%dij(:,:)             =zero

!    write(*,*)'paw2wvl_ij, erase me, set dij=hij'
!    !l2,m2, l1,m1 (l and m are diagonal in hgh)
!    !             !for this case i=1
!    wvl%paw%paw_ij(iatom)%dij(1 ,1)= 1.858811d0 !1,1, 1,1
!    wvl%paw%paw_ij(iatom)%dij(2 ,1)= zero       !2,1, 1,1
!    wvl%paw%paw_ij(iatom)%dij(3 ,1)=-0.005895d0 !2,1, 2,1
!    wvl%paw%paw_ij(iatom)%dij(4 ,1)= zero       !2,2, 1,1
!    wvl%paw%paw_ij(iatom)%dij(5 ,1)= zero       !2,2, 2,1
!    wvl%paw%paw_ij(iatom)%dij(6 ,1)=-0.005895d0 !2,2, 2,2
!    wvl%paw%paw_ij(iatom)%dij(7 ,1)= zero       !2,3, 1,1
!    wvl%paw%paw_ij(iatom)%dij(8 ,1)= zero       !2,3, 2,1
!    wvl%paw%paw_ij(iatom)%dij(9 ,1)= zero       !2,3, 2,2
!    wvl%paw%paw_ij(iatom)%dij(10,1)=-0.005895d0 !2,3, 2,3
!    check only 10 terms exists

   end do
!  Option==2: deallocate
 elseif(option==2) then
   do iatom=1,my_natom
     ABI_DEALLOCATE(wvl%paw%paw_ij(iatom)%dij)
   end do
   ABI_DATATYPE_DEALLOCATE(wvl%paw%paw_ij)
 else 
   message = 'paw2wvl_ij: option should be equal to 1 or 2'
   MSG_ERROR(message)
 end if


#endif

 DBG_EXIT("COLL")

end subroutine paw2wvl_ij
!!***
