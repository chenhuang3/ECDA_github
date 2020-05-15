!{\src2tex{textfont=tt}}
!!****f* ABINIT/twexpibi
!! NAME
!! twexpibi
!!
!! FUNCTION
!! Routine that computes exp(i (b_bra-b_ket).R) at each site.
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
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  my_natom :: controls atom loop in parallelism
!!  natom :: number of atoms in unit cell
!!  rprimd(3,3) :: dimensioned primitive translations of real space lattice
!!  xred(natom,3) :: reduced coordinates of atoms in unit cell
!!
!! OUTPUT
!!  calc_expibi(2,natom,6) :: phase factors at each atom for different vector shifts
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      initorbmag
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine twexpibi(calc_expibi,dkvecs,gprimd,my_natom,natom,rprimd,xred, &
&                  mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_self

 use m_paral_atom, only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'twexpibi'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: my_natom,natom
 integer,optional,intent(in) :: mpi_comm_atom
 real(dp),intent(out) :: calc_expibi(2,my_natom,6)
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: dkvecs(3,3),gprimd(3,3),rprimd(3,3),xred(3,natom)

!Local variables---------------------------
!scalars
 integer :: bdir,bsig,calc_indx,iatom,iatom_tot,kdir,kdx,ksig,mu,my_comm_atom
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: bdotr
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: bb(3),bcart(3),xcart(3)

! *************************************************************************

 calc_expibi(:,:,:) = zero

!calc_expibi(2,my_natom,12)
!used for PAW field calculations (distributed over atomic sites)
!stores the on-site phase factors arising from
!$\langle\phi_{i,k}|\phi_{j,k+\sigma_k k_k}\rangle$
!where $\sigma = \pm 1$. These overlaps arise in various Berry
!phase calculations of electric and magnetic polarization. The on-site
!phase factor is $\exp[i(\sigma_b k_b - \sigma_k k_k)\cdot I]$ where
!$I$ is the nuclear position. Only the following
!are computed and saved, in the given order:
!

!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.my_natom/=natom)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_self;if (present(mpi_comm_atom)) my_comm_atom=mpi_comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!the magnetic field case is much more complicated, as we need full 
!sig_b*b_k - sig_k*k_k. However, we never
!need b_k // k_k. 

 do iatom = 1, my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   do bdir = 1, 3
     bsig = -1 ! obtain other cases from complex conjugation
!    
     kdir = modulo(bdir,3)+1
!    
     do ksig = -1, 1, 2
       kdx = (ksig+3)/2 ! 
       calc_indx = 2*(bdir-1)+kdx

!      note the definition used for the k-dependence of the PAW basis functions:
!$|\phi_{i,k}\rangle = exp(-i k\cdot r)|\phi_i\rangle
!      see Umari, Gonze, and Pasquarello, PRB 69,235102 Eq. 23. Thus the k-vector on the
!      bra side enters as k, while on the ket side it enters as -k.
       bb(:) = bsig*dkvecs(:,bdir) - ksig*dkvecs(:,kdir)

!      get cartesian positions of atom in cell
       do mu=1,3
         bcart(mu)=dot_product(bb(:),gprimd(mu,:))
         xcart(mu)=dot_product(rprimd(mu,:),xred(:,iatom_tot))
       end do
       bdotr = dot_product(xcart,bcart)
!      here is exp(i b.R) for the given site
       calc_expibi(1,iatom,calc_indx) = cos(two_pi*bdotr)
       calc_expibi(2,iatom,calc_indx) = sin(two_pi*bdotr)

     end do ! end loop over ksig
   end do ! end loop over bdir

 end do ! end loop on my_natom

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 end subroutine twexpibi
!!***
