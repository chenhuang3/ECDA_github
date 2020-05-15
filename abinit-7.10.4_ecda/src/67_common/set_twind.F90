!{\src2tex{textfont=tt}}
!!****f* ABINIT/set_twind
!! NAME
!! set_twind
!!
!! FUNCTION
!! set index tables for mappings used in magnetization
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group
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
!!  dtbfield <type(bfield_type)> = variables related to Berry phase
!! 
!! NOTES
!!
!! PARENTS
!!      initorbmag
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine set_twind(dtbfield)

 use m_profiling_abi

 use defs_basis
 use m_bfield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'set_twind'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 type(bfield_type), intent(inout) :: dtbfield

!Local
 integer :: bdir,bdx,bsig,kdir,kdx,ksig,kstep,twdx

! *********************************************************************
!fill in indexing information
!for some quantities only six combinations of bra and ket shift vectors are stored
!see m_bfield%twexpibi for documentation
!
!1) -b1  - -b2
!2) -b1  - +b2
!3) -b2  - -b3
!4) -b2  - +b3
!5) -b3  - -b1
!6) -b3  - +b1
 dtbfield%twind(:,:) = 0
 dtbfield%indhk(:,:) = 0

!the following macro maps direction (vdir) and +/- (vsig) to
!the numbers 1-6 as follows:
!-k_1 -> 1
!+k_1 -> 2
!-k_2 -> 3
!+k_2 -> 4
!-k_3 -> 5
!+k_3 -> 6
#define PIND(vdir,vsig) 2*(vdir-1)+(vsig+3)/2

!here is indexing for structures like expibi( +i (sig_b*k_b - sig_k*k_k) . R)
!only 6 values are saved, namely -k_b - (+/- k_k), with 1,2 , 2,3 and 3,1 pairs
!of directions. the rest are obtained by complex conjugation of the proper
!term as indicatd by a negative value of the index number.

 dtbfield%twind(PIND(1,-1),PIND(2,-1)) = +1 ! -b1 - (-b2)
 dtbfield%twind(PIND(1,-1),PIND(2,+1)) = +2 ! -b1 - (+b2)
 dtbfield%twind(PIND(1,+1),PIND(2,-1)) = -2 ! +b1 - (-b2)
 dtbfield%twind(PIND(1,+1),PIND(2,+1)) = -1 ! +b1 - (+b2)

 dtbfield%twind(PIND(2,-1),PIND(1,-1)) = -1 ! -b2 - (-b1)
 dtbfield%twind(PIND(2,-1),PIND(1,+1)) = +2 ! -b2 - (+b1)
 dtbfield%twind(PIND(2,+1),PIND(1,-1)) = -2 ! +b2 - (-b1)
 dtbfield%twind(PIND(2,+1),PIND(1,+1)) = +1 ! +b2 - (+b1)

 dtbfield%twind(PIND(2,-1),PIND(3,-1)) = +3
 dtbfield%twind(PIND(2,-1),PIND(3,+1)) = +4
 dtbfield%twind(PIND(2,+1),PIND(3,-1)) = -4
 dtbfield%twind(PIND(2,+1),PIND(3,+1)) = -3

 dtbfield%twind(PIND(3,-1),PIND(2,-1)) = -3
 dtbfield%twind(PIND(3,-1),PIND(2,+1)) = +4
 dtbfield%twind(PIND(3,+1),PIND(2,-1)) = -4
 dtbfield%twind(PIND(3,+1),PIND(2,+1)) = +3

 dtbfield%twind(PIND(3,-1),PIND(1,-1)) = +5
 dtbfield%twind(PIND(3,-1),PIND(1,+1)) = +6
 dtbfield%twind(PIND(3,+1),PIND(1,-1)) = -6
 dtbfield%twind(PIND(3,+1),PIND(1,+1)) = -5
 
 dtbfield%twind(PIND(1,-1),PIND(3,-1)) = -5
 dtbfield%twind(PIND(1,-1),PIND(3,+1)) = +6
 dtbfield%twind(PIND(1,+1),PIND(3,-1)) = -6
 dtbfield%twind(PIND(1,+1),PIND(3,+1)) = +5

 do bdir = 1, 3
   do bsig = -1, 1, 2
     bdx = PIND(bdir,bsig)
     do kstep = 1, 2
       kdir = modulo(bdir-1+kstep,3)+1
       do ksig = -1, 1, 2
         kdx = 2*(kstep-1)+(ksig+3)/2
         twdx = 4*(bdx-1)+kdx
         dtbfield%indhk(PIND(bdir,bsig),PIND(kdir,ksig)) = twdx
       end do ! end loop on ksig
     end do ! end loop on kstep
   end do ! end loop on bsig
 end do ! end loop on bdir

end subroutine set_twind

!!***
