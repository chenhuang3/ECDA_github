!{\src2tex{textfont=tt}}
!!****f* ABINIT/expibr
!! NAME
!! expibr
!!
!! FUNCTION
!! Routine that computes exp(i (b_bra-b_ket).r) on the fine grid around
!! each PAW sphere. These quantities are stored in an bfield structure 
!! and used in PAW calculations of \hat{D}_ij in magnetization.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3) :: dimensioned primitive translations of reciprocal lattice
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom :: number of atoms in unit cell
!!  pawfgrtab <type(pawfgrtab_type)>= atomic data given on fine rectangular grid
!!  xred(natom,3) :: reduced coordinates of atoms in unit cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtbfield :: bfield structure 
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine expibr(dtbfield,gprimd,my_natom,natom,pawfgrtab,xred,&
&                  mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_xmpi, only : xmpi_self

 use m_bfield

 use m_pawfgrtab, only : pawfgrtab_type
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'expibr'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: my_natom,natom
 integer,optional,intent(in) :: mpi_comm_atom
 type(bfield_type),intent(inout) :: dtbfield
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gprimd(3,3),xred(3,natom)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(my_natom)

!Local variables---------------------------
!scalars
  integer :: calc_indx,iatom,iatom_tot,ic,bdir,bsig,kdir,kdx,ksig,my_comm_atom,mu,nfgd_max
  logical :: my_atmtab_allocated,paral_atom
  real(dp) :: phase,phase_xred
!arrays
  integer,pointer :: my_atmtab(:)
  real(dp) :: bb(3),bcart(3)

! *************************************************************************
!real(dp), pointer :: twexpibr(:,:,:,:)
!twexpibr(2,my_natom,nfgd,6)  (distributed over atomic sites)
!stores the on-site phase factors arising from
!$\langle\phi_{i,k+\sigma_b k_b}|\phi_{j,k+\sigma_k k_k}\rangle$
!where $\sigma = \pm 1$. The on-site phase factor is 
!$\exp[i( \sigma_b k_b - \sigma_k k_k)\cdot r]$ where
!$r$ is the position on the fine grid. Only 6 values are saved:
!-b1 - (-k2)
!-b1 - (+k2)
!-b2 - (-k3)
!-b2 - (+k3)
!-b3 - (-k1)
!-b3 - (+k1)

!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.my_natom/=natom)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_self;if (present(mpi_comm_atom)) my_comm_atom=mpi_comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 if (dtbfield%has_twexpibr == 0) then
   nfgd_max=0
   do iatom=1,my_natom
     if(nfgd_max<pawfgrtab(iatom)%nfgd) nfgd_max=pawfgrtab(iatom)%nfgd
   end do
   ABI_ALLOCATE(dtbfield%twexpibr,(2,my_natom,nfgd_max,6))
   dtbfield%has_twexpibr = 1
 end if

 if (dtbfield%has_twexpibr>0.and.my_natom>0) then
   dtbfield%twexpibr(:,:,:,:) = zero

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

         bb(:) = bsig*dtbfield%dkvecs(:,bdir) - ksig*dtbfield%dkvecs(:,kdir)

         do mu=1,3
           bcart(mu)=dot_product(bb(:),gprimd(mu,:))
         end do
         phase_xred = two_pi*dot_product(bb(:),xred(:,iatom_tot))

         do ic = 1, pawfgrtab(iatom)%nfgd
           phase = two_pi*dot_product(bcart(:),pawfgrtab(iatom)%rfgd(:,ic))+phase_xred
           dtbfield%twexpibr(1,iatom,ic,calc_indx)=cos(phase)
           dtbfield%twexpibr(2,iatom,ic,calc_indx)=sin(phase)
         end do ! end loop over nfgd for this atom
       end do ! end loop over ksig
     end do ! end loop over bdir
   end do ! end loop over natom

   dtbfield%has_twexpibr = 2

 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 end subroutine expibr
!!***
