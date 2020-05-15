!{\src2tex{textfont=tt}}
!!****f* abinit/sym_pawcprj_kn
!! NAME
!! sym_pawcprj_kn
!!
!! FUNCTION
!! compute cprj for a given band and k point based on cprj at a symmetry-related
!! k point.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprj_ikn (pawcprj_type) :: cprj for a single band and k point, typically a k point in the IBZ
!!  cprj_sym(4,nsym,natom) :: 1:3 shift, and 4 final atom, of symmetry isym operating on iatom 
!!                            (S^{-1}(R - t) = r0 + L, see symatm.F90
!!  dimlmn(natom) :: ln dimension of each atom
!!  iband :: number of bands to treat, use -1 to treat all nband bands
!!  indlmn(6,lmnmax,ntypat) :: n,l,m dimensions for each atom type (see psps type)
!!  isym :: symmetry element used in current application
!!  itim :: 1 if time reversal also used, 0 else
!!  kpt(3) :: kpt vector used
!!  lmax :: max l value 
!!  lmnmax :: max lmn value
!!  mband :: maximum number of bands
!!  natom :: number of atoms in cell
!!  nband :: number of bands in cprj_ikn
!!  nspinor :: number of spinors
!!  nsym :: total number of symmetry elements
!!  ntypat :: number of types of atoms
!!  typat(natom) :: type of each atom
!!  zarot(2*lmax+1,2*lmax+1,lmax+1,nsym) :: elements of rotation matrix for angular momentum states 
!!                                          and symmetry operations. See setsymrhoij.F90
!!
!! OUTPUT
!!  cprj_fkn (pawcprj_type) :: cprj for a single band and k point where the k point is related to 
!!    the input k point by a symmetry operation
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  This routine is based on M. Giantomassi's doctoral dissertation, formula 7.77. It is not clear
!!  whether it is implemented correctly for nonsymmorphic symmetries.
!!
!! PARENTS
!!      berryphase_new,cgwf,make_grad_berry,update_orbmag
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine sym_pawcprj_kn(cprj_fkn,cprj_ikn,cprj_sym,dimlmn,iband,indlmn,&
&                       isym,itim,kpt,lmax,lmnmax,mband,natom,nband,nspinor,nsym,ntypat,&
&                       typat,zarot)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sym_pawcprj_kn'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: iband,isym,itim,lmax,lmnmax,mband
 integer,intent(in) :: natom,nband,nspinor,nsym,ntypat

!arrays
 integer,intent(in) :: cprj_sym(4,nsym,natom),dimlmn(natom)
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),typat(natom)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(in) :: zarot(2*lmax+1,2*lmax+1,lmax+1,nsym) 
 type(pawcprj_type),intent(in) :: cprj_ikn(natom,mband*nspinor)
 type(pawcprj_type),intent(inout) :: cprj_fkn(natom,mband*nspinor) !vz_i

!Local variables---------------------------
!scalars
 integer :: iatom, ibct, ibnd, ibsp, ibst, icpgr, iin, il, il0, im
 integer :: ilmn, iln, iln0, ilpm, indexi, ispinor, itypat, jatom, mm, nlmn
 real(dp) :: kdotL, phr, phi

!arrays
 real(dp) :: rl(3), t1(2), t2(2)

! *************************************************************************

!write(std_out,*)' jwz debug: enter sym_pawcprj_kn '
 
 if (iband == -1) then
   ibst = 1
   ibnd = nband
 else
   ibst = iband
   ibnd = iband
 end if

 do iatom = 1, natom
   
   itypat = typat(iatom)
   nlmn = dimlmn(iatom)
   jatom = cprj_sym(4,isym,iatom)
   rl(:) = cprj_sym(1:3,isym,iatom)
   kdotL = dot_product(rl,kpt)
   phr = cos(two_pi*kdotL)
   phi = sin(two_pi*kdotL)

   il0 = -1; iln0 = -1; indexi = 1 
   do ilmn = 1, nlmn

     il = indlmn(1,ilmn,itypat)
     im = indlmn(2,ilmn,itypat)
     iin = indlmn(3,ilmn,itypat)
     iln = indlmn(5,ilmn,itypat)

     ilpm = 1 + il + im
     if (iln /= iln0) indexi = indexi + 2*il0 + 1

     do ibct = ibst, ibnd

       do ispinor = 1, nspinor

         ibsp = nspinor*(ibct-1) + ispinor

         t1(:) = zero
         do mm = 1, 2*il+1
           t1(1) = t1(1) + zarot(mm,ilpm,il+1,isym)*cprj_ikn(jatom,ibsp)%cp(1,indexi+mm)
           t1(2) = t1(2) + zarot(mm,ilpm,il+1,isym)*cprj_ikn(jatom,ibsp)%cp(2,indexi+mm)
         end do

         t2(1) = t1(1)*phr - t1(2)*phi
         t2(2) = t1(2)*phr + t1(1)*phi

         if (itim == 1) t2(2) = -t2(2)

         cprj_fkn(iatom,ibsp)%cp(1,ilmn) = t2(1)
         cprj_fkn(iatom,ibsp)%cp(2,ilmn) = t2(2)

! do same transformations for gradients of cprj_ikn
! note that ncpgr = 0 if no gradients present so this loop will not be executed
! in this case
         do icpgr = 1, cprj_ikn(jatom,ibsp)%ncpgr
           t1(:) = zero
           do mm = 1, 2*il+1
             t1(1) = t1(1) + zarot(mm,ilpm,il+1,isym)*cprj_ikn(jatom,ibsp)%dcp(1,icpgr,indexi+mm)
             t1(2) = t1(2) + zarot(mm,ilpm,il+1,isym)*cprj_ikn(jatom,ibsp)%dcp(2,icpgr,indexi+mm)
           end do
           
           t2(1) = t1(1)*phr - t1(2)*phi
           t2(2) = t1(2)*phr + t1(1)*phi
           
           if (itim == 1) t2(2) = -t2(2)
           
           cprj_fkn(iatom,ibsp)%dcp(1,icpgr,ilmn) = t2(1)
           cprj_fkn(iatom,ibsp)%dcp(2,icpgr,ilmn) = t2(2)

         end do ! end loop over ncpgr

       end do ! end loop over nspinor

     end do ! end loop over bands

     il0 = il; iln0 = iln

   end do ! end loop over ilmn
 end do ! end loop over atoms

!write(std_out,*)' jwz debug: leave sym_pawcprj_kn '

 end subroutine sym_pawcprj_kn 
!!***
