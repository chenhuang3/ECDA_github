!{\src2tex{textfont=tt}}
!!****f* ABINIT/expibi
!! NAME
!! expibi
!!
!! FUNCTION
!! Routine that computes exp(i (-b_ket).R) at each site.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dkvecs(3,3) :: $\Delta k$ increments in three directions for current k-point grid
!!  gprimd(3,3) :: dimensioned primitive translations of reciprocal lattice
!!  natom :: number of atoms in unit cell
!!  rprimd(3,3) :: dimensioned primitive translations of real space lattice
!!  xred(natom,3) :: reduced coordinates of atoms in unit cell
!!
!! OUTPUT
!!  calc_expibi(2,natom,3) :: phase factors at each atom for different vector shifts
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      initberry,initorbmag
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine expibi(calc_expibi,dkvecs,gprimd,natom,rprimd,xred)

 use m_profiling_abi
 use m_errors
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'expibi'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: natom
 real(dp),intent(out) :: calc_expibi(2,natom,3)
!arrays
 real(dp),intent(in) :: dkvecs(3,3),gprimd(3,3),rprimd(3,3),xred(3,natom)

!Local variables---------------------------
!scalars
 integer :: iatom,kdir,mu
 real(dp) :: bdotr
!arrays
 real(dp) :: bb(3),bcart(3),xcart(3)

! *************************************************************************

 calc_expibi(:,:,:) = zero

!calc_expibi(2,my_natom,3)
!used for PAW field calculations (distributed over atomic sites)
!stores the on-site phase factors arising from
!$\langle\phi_{i,k}|\phi_{j,k+\sigma_k k_k}\rangle$
!where $\sigma = \pm 1$. These overlaps arise in various Berry
!phase calculations of electric and magnetic polarization. The on-site
!phase factor is $\exp[-i\sigma_k k_k)\cdot I]$ where
!$I$ is the nuclear position. 

 do iatom = 1, natom

   do kdir = 1, 3    

!    note the definition used for the k-dependence of the PAW basis functions:
!$|\phi_{i,k}\rangle = exp(-i k\cdot r)|\phi_i\rangle
!    see Umari, Gonze, and Pasquarello, PRB 69,235102 Eq. 23. 
     bb(:) = -dkvecs(:,kdir)

!    get cartesian positions of atom in cell
     do mu=1,3
       bcart(mu)=dot_product(bb(:),gprimd(mu,:))
       xcart(mu)=dot_product(rprimd(mu,:),xred(:,iatom))
     end do
     bdotr = dot_product(xcart,bcart)
!    here is exp(i b.R) for the given site
     calc_expibi(1,iatom,kdir) = cos(two_pi*bdotr)
     calc_expibi(2,iatom,kdir) = sin(two_pi*bdotr)

   end do ! end loop over kdir

 end do ! end loop on my_natom

 end subroutine expibi
!!***
