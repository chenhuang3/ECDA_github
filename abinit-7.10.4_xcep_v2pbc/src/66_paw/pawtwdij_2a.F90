!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawtwdij_2a
!! NAME
!! pawtwdij_2a
!!
!! FUNCTION
!! compute phase-twisted contribution to ${D}^1_{ij}-\tilde{D}^1_{ij}$ due to 
!! Hartree potential of $\tilde{n}^1$ interacting with onsite densities.
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
!!  dtbfield <type(bfield_type)> = dtbfield%tweijkl is updated
!!
!! NOTES
!! This term corresponds to term (2a) of Eq. 49 in Torrent et al.,
!! CMS 42, 337 (2008), including phase shifts as needed for orbital
!! magnetization:
!! $\langle\phi_i|e^{-i\mathbf{(k-b).r}}v_H[n^1]|\phi_j\rangle - 
!!  \langle\tilde{phi}_i|e^{-i\mathbf{(k-b).r}}v_H[\tilde{n}^1]|\tilde{phi}_j\rangle$
!! where $\mathbf{b}$ is the bra shift vector and $\mathbf{k}$ is the ket shift vector, both are determined in
!! initorbmag.F90.
!!
!! PARENTS
!!      initorbmag
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,initylmr,pawrad_deducer0,poisson,realgaunt
!!      sbf8,simp_gen
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

 subroutine pawtwdij_2a(dtbfield,gprimd,natom,ntypat,pawrad,pawtab,typat,&
&                       mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_xmpi, only : xmpi_self

 use m_bfield
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_pawang, only : realgaunt
 use m_pawrad, only : pawrad_type, pawrad_deducer0, simp_gen, poisson
 use m_pawtab, only : pawtab_type
 use m_sphharm, only : initylmr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtwdij_2a'
 use interfaces_65_psp
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
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: bdir,bigl,bigm,biglmin,biglmax,bsig,biglm,clm
 integer :: eijkl_ind,gij,gkl,glm
 integer :: ijlmn,ijlm,ijln,ir,iatom,iatom_tot,itypat
 integer :: kdir,kllmn,kllm,klln,ksig,lcmax
 integer :: litl,litm,litlm,litlmin,litlmax,lt
 integer :: mesh_size,meshsz,mt,mu,my_comm_atom,ngnt
 integer :: tlm,ylmr_normchoice,ylmr_npts,ylmr_option
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: bessarg,knorm,tuujv,uujv,qq
 complex(dpc) :: afac,tijkl
!arrays
 integer,allocatable :: gntselect(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: bb(3),bcart(3),kb(3),kcart(3),kbn(3),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),j_bessel(:,:),realgnt(:),rvl(:),trvl(:),ylmk(:)
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

   lcmax = pawtab(itypat)%l_size
   ABI_ALLOCATE(ylmk,(lcmax*lcmax))
   ABI_ALLOCATE(gntselect,((2*lcmax-1)**2,lcmax**2*(lcmax**2+1)/2))
   ABI_ALLOCATE(realgnt,((2*lcmax-1)**2*(lcmax)**4))
   call realgaunt(lcmax,ngnt,gntselect,realgnt)

   mesh_size = pawrad(itypat)%mesh_size

   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(rvl,(mesh_size))
   ABI_ALLOCATE(trvl,(mesh_size))
   ABI_ALLOCATE(j_bessel,(mesh_size,lcmax))
   meshsz=pawrad(itypat)%int_meshsz;if (meshsz>mesh_size) ff(meshsz+1:mesh_size)=zero

   do bdir = 1, 3
     bsig = -1

     bb(:) = bsig*dtbfield%dkvecs(:,bdir) ! bra vector
     do mu=1,3
       bcart(mu)=dot_product(bb(:),gprimd(mu,:))
     end do
     
     do ksig = -1, 1, 2

       kdir = mod(bdir,3)+1
       eijkl_ind = dtbfield%twind(PIND(bdir,bsig),PIND(kdir,ksig))
       kb(:) = ksig*dtbfield%dkvecs(:,kdir)  ! ket vector
       do mu=1,3
         kcart(mu)=dot_product(kb(:),gprimd(mu,:))
       end do
!      form b_b - b_k vector
       kcart(1:3) = bcart(1:3) - kcart(1:3)
       knorm=dsqrt(dot_product(kcart,kcart))
       if (knorm < tol12) then
         kbn(:) = zero
         ylmk(:) = zero
         ylmk(1) = 1.d0/sqrt(four_pi)
       else
         kbn(:) = kcart(:)/knorm ! unit vector in kb direction
         call initylmr(lcmax,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,kbn,ylmk(:),ylmgr)
       end if
       knorm = two_pi*knorm ! re-normed kb for calls to bessel fnc

!      compute bessel functions here
       do ir = 1, mesh_size
         bessarg = knorm*pawrad(itypat)%rad(ir)
         call sbf8(lcmax,bessarg,j_bessel(ir,:))
       end do

       do ijlmn = 1, pawtab(itypat)%lmn2_size
         ijlm = pawtab(itypat)%indklmn(1,ijlmn)
         ijln = pawtab(itypat)%indklmn(2,ijlmn)
         biglmin = pawtab(itypat)%indklmn(3,ijlmn)
         biglmax = pawtab(itypat)%indklmn(4,ijlmn)

         do kllmn = 1, pawtab(itypat)%lmn2_size
           kllm = pawtab(itypat)%indklmn(1,kllmn)
           klln = pawtab(itypat)%indklmn(2,kllmn)
           litlmin = pawtab(itypat)%indklmn(3,kllmn)
           litlmax = pawtab(itypat)%indklmn(4,kllmn)

           tijkl = cmplx(zero,zero)

           do litl = litlmin, litlmax, 2
             
!            get v_l, remember that rvl = r*v_l
             do ir = 1, meshsz
               ff(ir) = four_pi*pawtab(itypat)%phiphj(ir,klln)
             end do
             call poisson(ff,litl,qq,pawrad(itypat),rvl) 
             do ir = 1, meshsz
               ff(ir) = four_pi*pawtab(itypat)%tphitphj(ir,klln)
             end do 
             call poisson(ff,litl,qq,pawrad(itypat),trvl) 
             
             do bigl = biglmin, biglmax, 2
               
               do lt = 0, lcmax-1

                 if( lt+bigl-litl >= 0 .and. &
&                 lt-bigl+litl >= 0 .and. &
&                 -lt+bigl+litl >= 0 ) then

!                  get radial integral (note division by r, because rv output from
!                  poisson gave r*vl, not vl)
                   do ir = 2, meshsz
                     ff(ir)=rvl(ir)*j_bessel(ir,lt+1)*pawtab(itypat)%phiphj(ir,ijln)/&
&                     pawrad(itypat)%rad(ir)
                   end do
                   call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
                   call simp_gen(uujv,ff,pawrad(itypat))
                   do ir = 2, meshsz
                     ff(ir)=trvl(ir)*j_bessel(ir,lt+1)*pawtab(itypat)%tphitphj(ir,ijln)/&
&                     pawrad(itypat)%rad(ir)
                   end do
                   call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
                   call simp_gen(tuujv,ff,pawrad(itypat))
                   
                   do bigm = -bigl, bigl
                     biglm = LMC(bigl,bigm)
                     gij = gntselect(biglm,ijlm)
                     if(gij == 0) cycle
                     
                     do litm = -litl, litl
                       litlm = LMC(litl,litm)
                       gkl = gntselect(litlm,kllm)
                       if(gkl == 0) cycle
                       
                       do mt = -lt, lt
                         tlm = LMC(lt,mt)
                         clm = RCC(tlm,litlm)
                         glm = gntselect(biglm,clm)
                         if(glm == 0) cycle
                         
                         tijkl = tijkl + four_pi*il(mod(lt,4))*ylmk(tlm)*&
&                         realgnt(gij)*realgnt(gkl)*realgnt(glm)*(uujv-tuujv)
                         
                       end do ! end loop over mt
                     end do ! end loop over litm
                   end do ! end loop over bigm
                 end if ! end check on triangle conditions
               end do ! end loop over lt
             end do ! end loop over big l
           end do ! end loop over little l

           do iatom = 1, dtbfield%my_natom
             iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
             if(typat(iatom_tot) == itypat) then
               afac = cmplx(dtbfield%twexpibi(1,iatom,eijkl_ind),&
&               dtbfield%twexpibi(2,iatom,eijkl_ind))
               dtbfield%tweijkl(1,ijlmn,kllmn,iatom,eijkl_ind) = &
&               dtbfield%tweijkl(1,ijlmn,kllmn,iatom,eijkl_ind) + &
&               real(tijkl*afac)
               dtbfield%tweijkl(2,ijlmn,kllmn,iatom,eijkl_ind) = &
&               dtbfield%tweijkl(2,ijlmn,kllmn,iatom,eijkl_ind) + &
&               aimag(tijkl*afac)             
             end if
           end do

         end do ! end loop over kl states

       end do ! end loop over ij states

     end do ! end loop over ksig
   end do ! end loop over bdir

   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(rvl)
   ABI_DEALLOCATE(trvl)
   ABI_DEALLOCATE(j_bessel)
   ABI_DEALLOCATE(ylmk)
   ABI_DEALLOCATE(gntselect)
   ABI_DEALLOCATE(realgnt)

 end do ! end loop on ntypat

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 end subroutine pawtwdij_2a
!!***
