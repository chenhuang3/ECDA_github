!{\src2tex{textfont=tt}}
!!****f* ABINIT/smatrix_k_paw_bfield
!! NAME
!! smatrix_k_paw_bfield
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  bdir :: integer giving direction along which overlap is computed for bra
!!  bfor :: integer indicating whether to compute forward (1) or backward (2) or
!!          neither (0) along kpt string
!!  cprj_k (pawcprj_type) :: cprj for occupied bands at point k
!!  cprj_kb :: cprj for occupied bands at point k+b
!!  dtbfield :: structure referring to all bfield variables
!!  kdir :: integer giving direction along which overlap is computed for ket
!!  kfor :: integer indicating whether to compute forward (1) or backward (2)
!!    along kpt string
!!  natom :: number of atoms in cell
!!  typat :: typat(natom) type of each atom
!!
!! OUTPUT
!! smat_k_paw :: array of the on-site PAW parts of the overlaps between Bloch states at points
!!   k and k+b, for the various pairs of bands, that is, the on-site part of 
!!   <u_nk|u_mk+b>
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine assumes that the cprj are not explicitly ordered by 
!! atom type.
!!
!! PARENTS
!!      update_orbmag
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

!the following macro maps direction (vdir) and +/- (vsig) to
! the numbers 1-6 as follows:
! -k_1 -> 1
! +k_1 -> 2
! -k_2 -> 3
! +k_2 -> 4
! -k_3 -> 5
! +k_3 -> 6
#define PIND(vdir,vsig) 2*(vdir-1)+(vsig+3)/2

 subroutine smatrix_k_paw_bfield(bdir,bfor,cprj_k,cprj_kb,dtbfield,kdir,kfor,mband,natom,smat_k_paw,typat)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_bfield
 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'smatrix_k_paw_bfield'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: bdir,bfor,kdir,kfor,mband,natom
 type(bfield_type),intent(in) :: dtbfield
 type(pawcprj_type),intent(in) :: cprj_k(natom,dtbfield%nspinor*mband)
 type(pawcprj_type),intent(in) :: cprj_kb(natom,dtbfield%nspinor*mband)

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(out) :: smat_k_paw(2,dtbfield%nband_occ,dtbfield%nband_occ)

!Local variables---------------------------
!scalars
 integer :: bsig,iatom,iband,ibs,ilmn,ispinor,itypat
 integer :: jband,jbs,jlmn,klmn,ksig,nspinor,twdx
 logical :: need_conjg
 complex(dpc) :: cpk,cpkb,cterm,paw_onsite

! *************************************************************************

!initialize smat_k_paw
 smat_k_paw(:,:,:) = zero

 nspinor = dtbfield%nspinor

 bsig = -2*bfor+3; ksig=-2*kfor+3
 if (bfor == 0) bsig = 0
 need_conjg = .false.
 if (bfor == 0) then
   twdx = kdir
   if (ksig == -1) need_conjg = .true.
 else
   twdx = dtbfield%twind(PIND(bdir,bsig),PIND(kdir,ksig))
   if(twdx < 0) then
     twdx = -twdx
     need_conjg = .true.
   end if
 end if

 do iatom = 1, natom
   itypat = typat(iatom)

   do ilmn=1,dtbfield%lmn_size(itypat)
     do jlmn=1,dtbfield%lmn_size(itypat)
       klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
       if (bfor == 0) then
         paw_onsite = cmplx(dtbfield%qijb_kk(1,klmn,iatom,twdx),&
&         dtbfield%qijb_kk(2,klmn,iatom,twdx))
       else
         paw_onsite = cmplx(dtbfield%twqijb_kk(1,klmn,iatom,twdx),&
&         dtbfield%twqijb_kk(2,klmn,iatom,twdx))
       end if
       if (need_conjg) paw_onsite = conjg(paw_onsite)
       do iband = 1, dtbfield%nband_occ
         do jband = 1, dtbfield%nband_occ
           do ispinor = 1, nspinor
             ibs = nspinor*(iband-1) + ispinor
             jbs = nspinor*(jband-1) + ispinor
             cpk=cmplx(cprj_k(iatom,ibs)%cp(1,ilmn),cprj_k(iatom,ibs)%cp(2,ilmn))
             cpkb=cmplx(cprj_kb(iatom,jbs)%cp(1,jlmn),cprj_kb(iatom,jbs)%cp(2,jlmn))
             cterm = conjg(cpk)*paw_onsite*cpkb
             smat_k_paw(1,iband,jband) = smat_k_paw(1,iband,jband)+real(cterm)
             smat_k_paw(2,iband,jband) = smat_k_paw(2,iband,jband)+aimag(cterm)
           end do ! end loop over ispinor
         end do ! end loop over jband
       end do ! end loop over iband
     end do ! end loop over ilmn
   end do ! end loop over jlmn

 end do ! end loop over atoms

 end subroutine    smatrix_k_paw_bfield
!!***

