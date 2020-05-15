!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawtwdij_2e
!! NAME
!! pawtwdij_2e
!!
!! FUNCTION
!! compute phase-twisted contribution to ${D}^1_{ij}-\tilde{D}^1_{ij}$ due to 
!! Hartree potential of core charge and compensation charge moments.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3) = primitive translations in recip space
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  natom = number of atoms in unit cell
!!  ntypat = number of types of atoms in unit cell
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat = typat(natom) list of atom types
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtbfield <type(bfield_type)> = dtbfield%twdij0 is updated
!!
!! NOTES
!! This term corresponds to term (2e) of Eq. 49 in Torrent et al.,
!! CMS 42, 337 (2008), including phase shifts as needed for orbital
!! magnetization:
!! $\sum_{LM}\int_{\Omega_R}e^{i\mathbf{b.r}}v_H[\tilde{n}_{Zc}]\hat{Q}^{LM}_{ij}(\mathbf{r})e^{-i\mathbf{k.r}}$
!! where $\mathbf{b}$ is the bra shift vector and $\mathbf{k}$ is the ket shift vector, both are determined in
!! initorbmag.F90.
!!
!! PARENTS
!!      initorbmag
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,initylmr,jbessel,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! macro to go from row-column indexing to combined indexing
#define RCC(glmn,hlmn) max(glmn,hlmn)*(max(glmn,hlmn)-1)/2+min(glmn,hlmn)

!macro to go from l,m angular momentum indexing to combined indexing
#define LMC(lval,mval) lval*lval+lval+mval+1

!the following macro maps direction (vdir) and +/- (vsig) to
! the numbers 1-6 as follows:
! -k_1 -> 1
! +k_1 -> 2
! -k_2 -> 3
! +k_2 -> 4
! -k_3 -> 5
! +k_3 -> 6
#define PIND(vdir,vsig) 2*(vdir-1)+(vsig+3)/2

#include "abi_common.h"

 subroutine pawtwdij_2e(dtbfield,gprimd,natom,ntypat,pawrad,pawtab,typat,&
&                       mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_xmpi, only : xmpi_self

 use m_bfield
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_pawrad, only : pawrad_type, simp_gen
 use m_pawtab, only : pawtab_type
 use m_paw_numeric, only: jbessel
 use m_sphharm, only : initylmr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtwdij_2e'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 integer,optional,intent(in) :: mpi_comm_atom
 type(bfield_type),intent(inout) :: dtbfield
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gprimd(3,3)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: bdir,bl,bm,blm,bln,blmn,bsig,clmn
 integer :: expibi_idx,iatom,iatom_tot,itypat,ilm2
 integer :: ir,kdir,kl,km,klm,kln,klmn,ksig
 integer :: lcmax,ll,mesh_size,mm,mu,my_comm_atom,twdx
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 logical :: my_atmtab_allocated,need_conjg,paral_atom
 real(dp) :: bessarg,intgrl,jlx,jldx,jldx2,knorm,xx
 complex(dpc) :: cexpibi,vijsum
!arrays 
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: bb(3),bcart(3),kb(3),kcart(3),kbn(3),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),ylmk(:)
! the following is (i)^L mod 4.
 complex(dpc),dimension(0:3) :: il(0:3)=(/(1.0,0.0),(0.0,1.0),(-1.0,0.0),(0.0,-1.0)/)

! *************************************************************************

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr

!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.(dtbfield%my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_self;if (present(mpi_comm_atom)) my_comm_atom=mpi_comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=dtbfield%my_natom)

 do itypat = 1, ntypat

   lcmax = pawtab(itypat)%l_size  ! lcmax - 1 is highest angular momentum state used in expansion
   ABI_ALLOCATE(ylmk,(lcmax*lcmax))

   indlmn => pawtab(itypat)%indlmn
   mesh_size = pawrad(itypat)%mesh_size
   ABI_ALLOCATE(ff,(mesh_size))

   do bdir = 1, 3
     do kdir = 1, 3
       if (kdir == bdir) cycle ! never need the kdir // bdir terms
       do bsig = -1, 1, 2

         bb(:) = bsig*dtbfield%dkvecs(:,bdir) ! bra vector
         do mu=1,3
           bcart(mu)=dot_product(bb(:),gprimd(mu,:))
         end do

         do ksig = -1, 1, 2

           twdx = dtbfield%indhk(PIND(bdir,bsig),PIND(kdir,ksig))
           expibi_idx = dtbfield%twind(PIND(bdir,bsig),PIND(kdir,ksig))
           need_conjg=.false.
           if(expibi_idx < 0) then
             expibi_idx = -expibi_idx
             need_conjg = .true.
           end if
           
           kb(:) = ksig*dtbfield%dkvecs(:,kdir)  ! ket vector
           do mu=1,3
             kcart(mu)=dot_product(kb(:),gprimd(mu,:))
           end do

!          form delta_k = b_b - b_k vector
           kcart(1:3) = bcart(1:3) - kcart(1:3)
           knorm=dsqrt(dot_product(kcart,kcart))
           if (knorm < tol12) then
             kbn(:) = zero
             ylmk(:) = zero; ylmk(1) = one/sqrt(four_pi)
           else
             kbn(:) = kcart(:)/knorm ! unit vector in kb direction
             call initylmr(lcmax,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,kbn,ylmk(:),ylmgr)
           end if

           knorm = two_pi*knorm ! re-normed kb for calls to bessel fnc

           do blmn = 1, pawtab(itypat)%lmn_size
             bl = indlmn(1,blmn)
             bm = indlmn(2,blmn)
             blm = indlmn(4,blmn)
             bln = indlmn(5,blmn)
             
             do klmn = 1, pawtab(itypat)%lmn_size
               kl = indlmn(1,klmn)
               km = indlmn(2,klmn)
               klm = indlmn(4,klmn)
               kln = indlmn(5,klmn)
               
               clmn = RCC(blmn,klmn)
               
               vijsum = cmplx(zero,zero)
               
               do ll = abs(kl-bl),kl+bl,2
                 do ir = 1, mesh_size
                   xx = pawrad(itypat)%rad(ir)
                   bessarg = xx*knorm
                   call jbessel(jlx,jldx,jldx2,ll,0,bessarg)
                   ff(ir) = xx*xx*jlx*pawtab(itypat)%VHtnZC(ir)*pawtab(itypat)%shapefunc(ir,ll+1)
                 end do
                 call simp_gen(intgrl,ff,pawrad(itypat))
                 
                 do mm = -ll, ll
                   ilm2 = LMC(ll,mm)
                   vijsum = vijsum + il(mod(ll,4))*pawtab(itypat)%qijl(ilm2,clmn)*ylmk(ilm2)*intgrl 
                 end do ! end loop over mm
               end do ! end loop over ll
               
               vijsum = 4.d0*pi*vijsum
               do iatom = 1, dtbfield%my_natom ! store result
                 iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
                 if(typat(iatom_tot) == itypat) then
                   cexpibi = cmplx(dtbfield%twexpibi(1,iatom,expibi_idx),&
&                   dtbfield%twexpibi(2,iatom,expibi_idx))
                   if(need_conjg) cexpibi = conjg(cexpibi)
                   vijsum = cexpibi*vijsum
                   dtbfield%twdij0(1,blmn,klmn,iatom,twdx) = dtbfield%twdij0(1,blmn,klmn,iatom,twdx) - &
&                   real(vijsum)
                   dtbfield%twdij0(2,blmn,klmn,iatom,twdx) = dtbfield%twdij0(2,blmn,klmn,iatom,twdx) - &
&                   aimag(vijsum)
                 end if
               end do
               
             end do ! end loop over ket states
           end do ! end loop over bra states

         end do ! end loop over ksig
       end do ! end loop over bsig

     end do ! end loop over kdir
   end do ! end loop over bdir

   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(ylmk)

 end do ! end loop on ntypat

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 end subroutine pawtwdij_2e
!!***
