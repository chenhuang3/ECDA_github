!{\src2tex{textfont=tt}}
!!****f* ABINIT/twqijb_kk
!! NAME
!! twqijb_kk
!!
!! FUNCTION
!! Routine which computes PAW onsite part of wavefunction overlap for Bloch
!! functions at two k-points k and k+b. These
!! quantities are used in PAW-based computations of polarization and magnetization.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dkvecs(3,3) :: $\Delta k$ increments in three directions for current k-point grid
!!  twexpibi(2,my_natom,6) :: phase factors at each atomic site for given k offset
!!  gprimd(3,3)=dimensioned primitive translations of reciprocal lattice
!!  lmn2max :: lmnmax*(lmnmax+1)/2
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  my_natom :: controls atom loop in parallelism
!!  natom=number of atoms in unit cell
!!  ntypat=number of types of atoms in unit cell
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat=typat(natom) list of atom types
!!
!! OUTPUT
!!  calc_qijb(2,lmn2max,natom,6) :: PAW on-site overlaps of wavefunctions at neighboring
!!                                   k points 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      initorbmag
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,initylmr,sbf8,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine twqijb_kk(calc_qijb,dkvecs,twexpibi,gprimd,lmn2max,my_natom,natom,ntypat,&
&                   pawang,pawrad,pawtab,typat, &
&                   mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_xmpi, only : xmpi_self,xmpi_sum
 use m_pawang, only : pawang_type
 use m_pawrad, only : pawrad_type, simp_gen
 use m_pawtab, only : pawtab_type
 use m_sphharm, only : initylmr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'twqijb_kk'
 use interfaces_65_psp
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: lmn2max,my_natom,natom,ntypat
 integer,optional,intent(in) :: mpi_comm_atom
 type(pawang_type),intent(in) :: pawang
 real(dp),intent(out) :: calc_qijb(2,lmn2max,natom,6)
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: dkvecs(3,3),twexpibi(2,my_natom,6),gprimd(3,3)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: bdir,bsig,calc_indx
!integer :: bdirend,bdirstart,bdx,bfor,bsigend,bsigstart 
 integer :: iatom,iatom_tot,ir,isel,itypat,kdir
!integer :: ierr
 integer :: kdx,ksig
!integer :: kfor,ksigend,ksigstart,kstep,kstepstart,kstepen
 integer :: klm,kln,klmn,lbess,lbesslm,lmin,lmax,mbess,mesh_size,mu,my_comm_atom
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 logical :: my_atmtab_allocated,paral_atom
!character(len=500) :: message
 real(dp) :: arg,bessg,bnorm,intg,rterm
 complex(dpc) :: cterm,etb,ifac
!arrays 
 integer,pointer :: my_atmtab(:)
 real(dp) :: bb(3),bbn(3),bcart(3),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),j_bessel(:,:),ylmb(:)
! the following is (i)^L mod 4.
 complex(dpc),dimension(0:3) :: il(0:3)=(/(1.0,0.0),(0.0,1.0),(-1.0,0.0),(0.0,-1.0)/)

! *************************************************************************
 
!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_self;if (present(mpi_comm_atom)) my_comm_atom=mpi_comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 calc_qijb(:,:,:,:) = zero

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr

!If field_type is 3 (electric field) then calc_qijb returns
!b1   b2
!(1)   0  +k1
!(2)   0  +k2
!(3)   0  +k3
!
!If field_type is 9 (magnetic field) then calc_qijb returns
!b1   b2

 do iatom = 1, my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   itypat = typat(iatom_tot)
   mesh_size = pawrad(itypat)%mesh_size

   ABI_ALLOCATE(j_bessel,(mesh_size,pawang%l_size_max))
   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(ylmb,(pawang%l_size_max*pawang%l_size_max))

   do bdir = 1, 3
     bsig = -1 ! obtain other cases from complex conjugation
!    
     kdir = modulo(bdir,3)+1
!    
     do ksig = -1, 1, 2
       kdx = (ksig+3)/2 ! 
       calc_indx = 2*(bdir-1)+kdx

!      here is exp(-i b.R) for current atom: recall storage in expibi
       etb = cmplx(twexpibi(1,iatom,calc_indx),twexpibi(2,iatom,calc_indx))

!      note the definition used for the k-dependence of the PAW basis functions:
!$|\phi_{i,k}\rangle = exp(-i k\cdot r)|\phi_i\rangle
!      see Umari, Gonze, and Pasquarello, PRB 69,235102 Eq. 23. Thus the k-vector on the
!      bra side enters as k, while on the ket side it enters as -k.
       bb(:) = bsig*dkvecs(:,bdir) - ksig*dkvecs(:,kdir)

!      get cartesian positions of atom in cell
       do mu=1,3
         bcart(mu)=dot_product(bb(:),gprimd(mu,:))
       end do

!      bbn is b-hat (the unit vector in the b direction) 
       bnorm=dsqrt(dot_product(bcart,bcart))
       bbn(:) = bcart(:)/bnorm

!      as an argument to the bessel function, need 2pi*b*r = 1 so b is re-normed to two_pi
       bnorm = two_pi*bnorm
       do ir=1,mesh_size
         arg=bnorm*pawrad(itypat)%rad(ir)
         call sbf8(pawang%l_size_max,arg,j_bessel(ir,:)) ! spherical bessel functions at each mesh point
       end do ! end loop over mesh
!      compute Y_LM(b) here
       call initylmr(pawang%l_size_max,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,bbn,ylmb(:),ylmgr)
       
       do klmn = 1, pawtab(itypat)%lmn2_size
         klm =pawtab(itypat)%indklmn(1,klmn)
         kln =pawtab(itypat)%indklmn(2,klmn)
         lmin=pawtab(itypat)%indklmn(3,klmn)
         lmax=pawtab(itypat)%indklmn(4,klmn)
         do lbess = lmin, lmax, 2    ! only possible choices for L s.t. Gaunt integrals
!          will be non-zero
           ifac = il(mod(lbess,4))
           do mbess = -lbess, lbess
             lbesslm = lbess*lbess+lbess+mbess+1
             isel=pawang%gntselect(lbesslm,klm)
             if (isel > 0) then
               bessg = pawang%realgnt(isel)
               ff(1:mesh_size)=(pawtab(itypat)%phiphj(1:mesh_size,kln)&
&               -pawtab(itypat)%tphitphj(1:mesh_size,kln))&
&               *j_bessel(1:mesh_size,lbess+1)
               call simp_gen(intg,ff,pawrad(itypat))
               rterm = four_pi*bessg*intg*ylmb(lbesslm)
               cterm = etb*ifac*rterm
               calc_qijb(1,klmn,iatom,calc_indx) = &
&               calc_qijb(1,klmn,iatom,calc_indx) + real(cterm)
               calc_qijb(2,klmn,iatom,calc_indx) = &
&               calc_qijb(2,klmn,iatom,calc_indx) + aimag(cterm)
               
             end if ! end selection on non-zero Gaunt factors
           end do ! end loop on mbess = -lbess, lbess
         end do ! end loop on lmin-lmax bessel l values
       end do ! end loop on lmn2_size klmn basis pairs

     end do ! end loop over ksig
   end do ! end loop over bdir

   ABI_DEALLOCATE(j_bessel)
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(ylmb)

 end do ! end loop over atoms

!if (paral_atom) then
!  call xmpi_sum(calc_qijb,my_comm_atom,ierr)
!end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 end subroutine twqijb_kk
!!***

