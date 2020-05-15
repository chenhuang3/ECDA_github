!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawtwdij
!! NAME
!! pawtwdij
!!
!! FUNCTION
!! Compute the pseudopotential strengths Dij of the PAW non local operator,
!! including phase twists due to different k points as needed in orbital 
!! magnetization.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  vtrial(cplex*nfft,nspden)=GS potential
!!  vxc(cplex*nfft,nspden)=XC potential (Hartree) on the fine FFT mesh
!!
!! OUTPUT
!!        dtbfield%twdij(...,1) contains Dij^up-up
!!        dtbfield%twdij(...,2) contains Dij^dn-dn
!!        dtbfield%twdij(...,3) contains Dij^up-dn (only if nspinor=2)
!!        dtbfield%twdij(...,4) contains Dij^dn-up (only if nspinor=2)
!!
!! NOTES
!!      Parallelisation over atoms (MT, april 2012):
!!      As the MPI communicator over atoms can be the same as the communicator
!!      over k-points, dtfield%twdij cannot be distributed over atomic sites
!!      (to do this, it is necessary to modify update_mmat.F90)
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,initylmr,nderiv_gen,pawgylm
!!      pawrad_deducer0,realgaunt,sbf8,simp_gen,xmpi_sum
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

subroutine pawtwdij(cplex,dtbfield,gprimd,my_natom,natom,nfft,nspden,ntypat,&
&                   paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawspnorb,pawtab,vtrial,vxc,&
&                   mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_self,xmpi_sum

 use m_bfield

 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_pawang, only       : pawang_type, realgaunt
 use m_pawrad, only       : pawrad_type, pawrad_deducer0, simp_gen, nderiv_gen
 use m_pawtab, only       : pawtab_type
 use m_paw_an, only       : paw_an_type
 use m_pawfgrtab, only    : pawfgrtab_type
 use m_pawrhoij, only     : pawrhoij_type
 use m_paw_finegrid, only : pawgylm
 use m_sphharm, only      : initylmr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtwdij'
 use interfaces_65_psp
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,my_natom,natom,nfft,nspden,ntypat,pawspnorb
 integer,optional,intent(in) :: mpi_comm_atom
 type(bfield_type),intent(inout) :: dtbfield
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gprimd(3,3),vtrial(cplex*nfft,nspden)
 real(dp),intent(in) :: vxc(cplex*nfft,nspden)
 type(paw_an_type),intent(in) :: paw_an(my_natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrad_type), intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: bigl,biglt,bigltmin,bigltmax,bigm,bigmt,bdir,bsig,cbiglit,cbiglt,cbigl
 integer :: clitl,cljm,gbiglt,gbiglit,glm,iatom,iatom_tot,ic,ierr,idij,ir,itypat,indhk
 integer :: iil,iim,iilm,iiln,iilmn,ijlm,ijlmn,ilslm,jc,jjl,jjm,jjlm,jjln,jjlmn
 integer :: irhoij,jrhoij,kdir,ksig,klmn,lcmax,lilj,ljlj,litl,litm,lm_size,lmn_size,lmn2_size
 integer :: mesh_size,meshsz,mu,my_comm_atom,ngnt,nfgd,nspdiag,nsploop,twind
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 real(dp) :: bessarg,fact,knorm,rhoij,uujv,vr
 real(dp), parameter :: HalfFineStruct2=half/InvFineStruct**2
 complex(dpc) :: afac,dijhartree,dijxc,eijkl,lsme,tij
 logical :: my_atmtab_allocated,need_conjg,paral_atom
!arrays
 integer,allocatable :: gntselect(:,:)
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:) 
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: have_radint_3a(:,:)
 real(dp) :: bb(3),bcart(3),dij0(2),kb(3),kcart(3),kbn(3),gylmgr2(1),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: dv1dr(:),ff(:),j_bessel(:,:),radint_3a(:,:),realgnt(:),ylmk(:)
 complex(dpc) :: tijso(4)
! the following is (i)^L mod 4.
 complex(dpc),dimension(0:3) :: il(0:3)=(/(1.0,0.0),(0.0,1.0),(-1.0,0.0),(0.0,-1.0)/)

! *************************************************************************

!DBG_ENTER("COLL")

 dtbfield%twdij = zero

 nsploop = 1
 if (dtbfield%ndij == 4) nsploop = 4

!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.(dtbfield%my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_self;if (present(mpi_comm_atom)) my_comm_atom=mpi_comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=dtbfield%my_natom)

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr

!------------------------------------------------------------------------
!----- Big loop over atoms
!------------------------------------------------------------------------

 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   lm_size=paw_an(iatom)%lm_size
   lmn_size = pawrhoij(iatom)%lmn_size
   lmn2_size = pawrhoij(iatom)%lmn2_size
   nfgd=pawfgrtab(iatom)%nfgd
   nspdiag=1;if (pawrhoij(iatom)%nspden==2) nspdiag=2

   itypat=pawrhoij(iatom)%itypat
   indlmn => pawtab(itypat)%indlmn

   lcmax = pawtab(itypat)%l_size
   ABI_ALLOCATE(ylmk,(lcmax*lcmax))
   ABI_ALLOCATE(gntselect,((2*lcmax-1)**2,lcmax**2*(lcmax**2+1)/2))
   ABI_ALLOCATE(realgnt,((2*lcmax-1)**2*(lcmax)**4))
   call realgaunt(lcmax,ngnt,gntselect,realgnt)

   mesh_size = pawrad(itypat)%mesh_size
   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(j_bessel,(mesh_size,lcmax))
   meshsz=pawrad(itypat)%int_meshsz;if (meshsz>mesh_size) ff(meshsz+1:mesh_size)=zero

!  Compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (pawfgrtab(iatom)%gylm_allocated==0) then
     if (allocated(pawfgrtab(iatom)%gylm))  then
       ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,lm_size))
     pawfgrtab(iatom)%gylm_allocated=2
     call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,gylmgr2,&
&     lm_size,nfgd,1,0,0,pawtab(itypat),pawfgrtab(iatom)%rfgd)
   end if

   do bdir = 1, 3
     do kdir = 1, 3
       if (kdir == bdir) cycle ! never need the kdir // bdir terms
       
       do bsig = -1, 1, 2
         
         bb(:) = bsig*dtbfield%dkvecs(:,bdir) ! bra vector
         do mu=1,3
           bcart(mu)=dot_product(bb(:),gprimd(mu,:))
         end do
         
         do ksig = -1, 1, 2
           
           indhk = dtbfield%indhk(PIND(bdir,bsig),PIND(kdir,ksig))
           twind = dtbfield%twind(PIND(bdir,bsig),PIND(kdir,ksig))
           need_conjg=.false.
           if(twind < 0) then
             twind = -twind
             need_conjg = .true.
           end if
!          here is the atomic site twist factor
           afac = cmplx(dtbfield%twexpibi(1,iatom,twind),&
&           dtbfield%twexpibi(2,iatom,twind))
           if(need_conjg) afac = conjg(afac)
           
           kb(:) = dtbfield%dkvecs(:,kdir)  ! ket vector
           do mu=1,3
             kcart(mu)=dot_product(kb(:),gprimd(mu,:))
           end do
!          form b_b - b_k vector
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
!          compute bessel functions here
           do ir = 1, mesh_size
             bessarg = knorm*pawrad(itypat)%rad(ir)
             call sbf8(lcmax,bessarg,j_bessel(ir,:))
           end do

           
!          ------------------------------------------------------------------------
!          ----------- Load atomic Dij0 into Dij
!          -- these are the phase-twisted versions of terms 1, 2b, and 2e from
!          Torrent et al, CMS 42 337-351 (2008)
!          ------------------------------------------------------------------------

           do iilmn = 1, lmn_size
             do jjlmn = 1, lmn_size

               dij0(1) = dtbfield%twdij0(1,iilmn,jjlmn,iatom_tot,indhk)
               dij0(2) = dtbfield%twdij0(2,iilmn,jjlmn,iatom_tot,indhk)

               do idij=1,MIN(2,nsploop) ! if nsploop = 1, just one term. if 4, then two terms: up-up, down-down
                 dtbfield%twdij(1,iilmn,jjlmn,iatom_tot,indhk,idij) = dij0(1)
                 dtbfield%twdij(2,iilmn,jjlmn,iatom_tot,indhk,idij) = dij0(2)
               end do
               
             end do
           end do

!          ------------------------------------------------------------------------
!          ----------- Add Dij_Hartree to Dij
!          -- these are the phase-twisted versions of terms 2a, 2c, 2d, and 2f from
!          Torrent et al, CMS 42 337-351 (2008)
!          ------------------------------------------------------------------------

           do iilmn = 1, lmn_size
             do jjlmn = 1, lmn_size
               
               ijlmn = RCC(iilmn,jjlmn)
               jrhoij=1
               do irhoij=1,pawrhoij(iatom)%nrhoijsel
                 klmn=pawrhoij(iatom)%rhoijselect(irhoij)
!                in pawrhoij nspden is assumed == 1
                 rhoij=pawrhoij(iatom)%rhoijp(jrhoij,1)*pawtab(itypat)%dltij(klmn)
                 eijkl = cmplx(dtbfield%tweijkl(1,ijlmn,klmn,iatom,twind),&
&                 dtbfield%tweijkl(2,ijlmn,klmn,iatom,twind))
                 if(need_conjg) eijkl = conjg(eijkl)

                 dijhartree = rhoij*eijkl
                 do idij=1,MIN(2,nsploop) ! if nsploop = 1, just one term. if 4, then two terms: up-up, down-down
                   dtbfield%twdij(1,iilmn,jjlmn,iatom_tot,indhk,idij) = &
&                   dtbfield%twdij(1,iilmn,jjlmn,iatom_tot,indhk,idij) + real(dijhartree)
                   dtbfield%twdij(2,iilmn,jjlmn,iatom_tot,indhk,idij) = &
&                   dtbfield%twdij(2,iilmn,jjlmn,iatom_tot,indhk,idij) + aimag(dijhartree)
                 end do

               end do ! end loop over nrhoijsel

             end do ! end loop over jjlmn
           end do ! end loop over iilmn

!          ------------------------------------------------------------------------
!          ----------- Add Dij_xc to Dij
!          -- these are the phase-twisted versions of terms 3a and 3b from
!          Torrent et al, CMS 42 337-351 (2008)
!          note: we require (see chkinp.F90) that:
!          usexcnhat = 0 (thus \hat{n} DOES NOT appear in any v_xc, see also
!          Torrent et al., Comp. Phys. Commun. 181 1862-1867 (2010)
!          pawxcdev = 1 so v_xc is delivered as moments, not as values on a grid
!          because usexcnhat = 0, term 3b is also zero. This is because this term
!          arises from a functional derivative of E_xc[\tilde{n}^1 + \hat{n} + \tilde{n}_c],
!          but in the usexcnhat = 0 case, E_xc DOES NOT depend on \hat{n} (the compensation
!          charge \hat{n} only appears in the Hartree terms, to cancel Coulomb interactions
!          between neighboring PAW spheres). 
!          ------------------------------------------------------------------------  

           do iilmn = 1, lmn_size
             iil = indlmn(1,iilmn)
             iim = indlmn(2,iilmn)
             iilm = indlmn(4,iilmn)
             iiln = indlmn(5,iilmn)

             do jjlmn = 1, lmn_size
               jjl = indlmn(1,jjlmn)
               jjm = indlmn(2,jjlmn)
               jjlm = indlmn(4,jjlmn)
               jjln = indlmn(5,jjlmn)

               ijlm = RCC(iilm,jjlm)
               ijlmn = RCC(iilmn,jjlmn)

               bigltmin = abs(iil-jjl)
               bigltmax = iil+jjl

!              radint_3a and have_radint_3a hold the radial integrals to avoid recomputation
!              radial integral is int u_i u_j bessel(l;b*r) v_xc^LM(r) 
!              thus the dimension is for little l (0:bigltmax) and the coupled LM states 
               ABI_ALLOCATE(radint_3a,(0:bigltmax,(bigltmax+1)**2))
               ABI_ALLOCATE(have_radint_3a,(0:bigltmax,(bigltmax+1)**2))
               radint_3a(:,:) = zero
               have_radint_3a(:,:) = .FALSE.

               tij = cmplx(zero,zero)

               do biglt = bigltmin, bigltmax, 2
                 do bigmt = -biglt, biglt
                   cbiglt = LMC(biglt,bigmt)
                   gbiglt = gntselect(cbiglt,ijlm)
                   if(gbiglt == 0) cycle

                   do bigl = 0, biglt
                     do bigm = -bigl, bigl
                       cbigl = LMC(bigl,bigm)
                       do litl = 0, biglt
                         if (.NOT.have_radint_3a(litl,cbigl)) then
                           do ir = 1, mesh_size ! we are assuming here cplex = 1 and nspden = 1
                             ff(ir) = paw_an(iatom)%vxc1(ir,cbigl,1)*&
&                             pawtab(itypat)%phi(ir,iiln)*pawtab(itypat)%phi(ir,jjln)*j_bessel(ir,litl+1) - &
&                             paw_an(iatom)%vxct1(ir,cbigl,1)*&
&                             pawtab(itypat)%tphi(ir,iiln)*pawtab(itypat)%tphi(ir,jjln)*j_bessel(ir,litl+1)
                           end do
                           call simp_gen(uujv,ff,pawrad(itypat))
                           radint_3a(litl,cbigl) = uujv
                           have_radint_3a(litl,cbigl) = .TRUE.
                         end if
                         do litm = -litl, litl
                           clitl = LMC(litl,litm)
                           cbiglit = RCC(cbigl,clitl)
                           gbiglit = gntselect(cbiglt,cbiglit)

                           if(gbiglit == 0) cycle

                           tij = tij + four_pi*il(mod(litl,4))*ylmk(clitl)*&
&                           realgnt(gbiglt)*realgnt(gbiglit)*radint_3a(litl,cbigl)

                         end do ! end loop over litm
                       end do ! end loop over litl
                     end do ! end loop over bigm
                   end do ! end loop over bigl
                 end do ! end loop over bigmt
               end do ! end loop over biglt

               dijxc = afac*tij

               do idij=1,MIN(2,nsploop) ! if nsploop = 1, just one term. if 4, then two terms: up-up, down-down

                 dtbfield%twdij(1,iilmn,jjlmn,iatom,indhk,idij) = &
&                 dtbfield%twdij(1,iilmn,jjlmn,iatom,indhk,idij) + real(dijxc)
                 dtbfield%twdij(2,iilmn,jjlmn,iatom,indhk,idij) = &
&                 dtbfield%twdij(2,iilmn,jjlmn,iatom,indhk,idij) + aimag(dijxc)

               end do

               ABI_DEALLOCATE(radint_3a)
               ABI_DEALLOCATE(have_radint_3a)

             end do ! end loop over jjlmn
           end do ! end loop over iilmn

!          ------------------------------------------------------------------------
!          ----------- Add Dij_hat to Dij
!          ------------------------------------------------------------------------

           do iilmn = 1, lmn_size
             do jjlmn = 1, lmn_size

               ijlmn = RCC(iilmn,jjlmn)
               tij = cmplx(zero,zero) 
               
               do ilslm=1,lm_size
                 do ic=1,nfgd
                   jc=pawfgrtab(iatom)%ifftsph(ic)
!                  note that vxc is subtracted from vtrial, this is because usexcnhat = 0
                   vr=vtrial(jc,1)-vxc(jc,1) ! nspden = 1 implicitly
                   afac=cmplx(dtbfield%twexpibr(1,iatom,ic,twind),&
&                   dtbfield%twexpibr(2,iatom,ic,twind))
                   if(need_conjg) afac=conjg(afac)
                   tij=tij+afac*vr*pawtab(itypat)%qijl(ilslm,ijlmn)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end do
               end do

               do idij=1,MIN(2,nsploop) ! if nsploop = 1, just one term. if 4, then two terms: up-up, down-down

                 dtbfield%twdij(1,iilmn,jjlmn,iatom,indhk,idij) = &
&                 dtbfield%twdij(1,iilmn,jjlmn,iatom,indhk,idij) + real(tij)
                 dtbfield%twdij(2,iilmn,jjlmn,iatom,indhk,idij) = &
&                 dtbfield%twdij(2,iilmn,jjlmn,iatom,indhk,idij) + aimag(tij)
               end do 

             end do ! end loop over jjlmn
           end do ! end loop over iilmn

!          ------------------------------------------------------------------------
!          ----------- Add Dij_SO to Dij
!          ------------------------------------------------------------------------

           if (pawspnorb == 1) then 

             do iilmn = 1, lmn_size
               iil = indlmn(1,iilmn)
               iim = indlmn(2,iilmn)
               iilm = indlmn(4,iilmn)
               iiln = indlmn(5,iilmn)

               do jjlmn = 1, lmn_size
                 jjl = indlmn(1,jjlmn)
                 jjm = indlmn(2,jjlmn)
                 jjlm = indlmn(4,jjlmn)
                 jjln = indlmn(5,jjlmn)

                 tij = cmplx(zero,zero)

                 do litl = abs(iil-jjl),iil+jjl,2
                   ABI_ALLOCATE(dv1dr,(mesh_size))
                   fact=one/sqrt(four_pi) ! Y_00
!                  L=0 moment of vxc1 and vh1, nspden = 1 assumed
                   ff(1:mesh_size)=paw_an(iatom)%vxc1(1:mesh_size,1,1)                  
                   ff(1:mesh_size)=fact*(ff(1:mesh_size)+paw_an(iatom)%vh1(1:mesh_size,1,1)) 
                   call nderiv_gen(dv1dr,ff,1,pawrad(itypat))
                   do ir = 2, mesh_size
                     dv1dr(ir)=HalfFineStruct2*(one/(one-ff(ir)/InvFineStruct**2)) &
&                     *dv1dr(ir)/pawrad(itypat)%rad(ir)
                   end do
                   call pawrad_deducer0(dv1dr,mesh_size,pawrad(itypat))
                   do ir = 1, mesh_size
                     ff(ir)= dv1dr(ir)*j_bessel(ir,litl+1)*&
&                     pawtab(itypat)%phi(ir,iiln)*pawtab(itypat)%phi(ir,jjln)
                   end do
                   call simp_gen(uujv,ff,pawrad(itypat))
                   ABI_DEALLOCATE(dv1dr)
                   do litm = -litl, litl
                     clitl = LMC(litl,litm) ! joint {lm}
                     do bigmt = -jjl,jjl
                       cljm = LMC(jjl,bigmt) ! joint {l_j,\tilde{M}}
                       lilj = RCC(iilm,cljm) ! joint {l_im_i,l_j\tilde{M}}
                       glm = gntselect(clitl,lilj) 
                       if(glm == 0) cycle
                       do idij = 1, nsploop ! nsploop should always be 4 if we are doing SO calc
                         ljlj = RCC(cljm,jjlm) ! joint {l_j\tilde{M},l_jm_j} 
                         fact = uujv; if (cljm > jjlm) fact = -fact
                         select case (idij)
                           case (1) ! up-up
                             lsme = cmplx(pawang%ls_ylm(1,ljlj,1),&
&                             pawang%ls_ylm(2,ljlj,1))
                             tijso(1) = tijso(1) + four_pi*il(mod(litl,4))*ylmk(clitl)*&
&                             realgnt(glm)*fact*lsme
                           case (2) ! down-down
                             lsme = cmplx(-pawang%ls_ylm(1,ljlj,1),&
&                             -pawang%ls_ylm(2,ljlj,1))
                             tijso(2) = tijso(2) + four_pi*il(mod(litl,4))*ylmk(clitl)*&
&                             realgnt(glm)*fact*lsme
                           case(3) ! up-down 
                             lsme = cmplx(pawang%ls_ylm(1,ljlj,2),&
&                             pawang%ls_ylm(2,ljlj,2))
                             tijso(3) = tijso(3) + four_pi*il(mod(litl,4))*ylmk(clitl)*&
&                             realgnt(glm)*fact*lsme
                           case(4) ! down-up
                             lsme = cmplx(-pawang%ls_ylm(1,ljlj,2),&
&                             pawang%ls_ylm(2,ljlj,2))
                             tijso(4) = tijso(4) + four_pi*il(mod(litl,4))*ylmk(clitl)*&
&                             realgnt(glm)*fact*lsme
                         end select
                       end do ! end loop over nsploop
                     end do ! end loop over bigmt
                   end do ! end loop over litm
                 end do ! end loop over litl

                 do idij=1,nsploop ! nsploop = 4 in SO case
                   dtbfield%twdij(1,iilmn,jjlmn,iatom,indhk,idij) = &
&                   dtbfield%twdij(1,iilmn,jjlmn,iatom,indhk,idij) + real(afac*tijso(idij))
                   dtbfield%twdij(2,iilmn,jjlmn,iatom,indhk,idij) = &
&                   dtbfield%twdij(2,iilmn,jjlmn,iatom,indhk,idij) + aimag(afac*tijso(idij))
                 end do 

               end do ! end loop over jjlmn
             end do ! end loop over iilmn

           end if ! end check on spin-orbit case

         end do ! end loop over ksig
       end do ! end loop over bsig
     end do ! end loop over kdir
   end do ! end loop over bdir

   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(j_bessel)
   ABI_DEALLOCATE(ylmk)
   ABI_DEALLOCATE(gntselect)
   ABI_DEALLOCATE(realgnt)

 end do ! end loop over my_natom

 if (paral_atom) then
   call xmpi_sum(dtbfield%twdij,my_comm_atom,ierr)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

!DBG_EXIT("COLL")

end subroutine pawtwdij
!!***
