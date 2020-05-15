!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawtwdij_1
!! NAME
!! pawtwdij_1
!!
!! FUNCTION
!! compute phase-twisted contribution to ${D}^1_{ij}-\tilde{D}^1_{ij}$ due to 
!! kinetic energy.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt.
!!
!! INPUTS
!!  gprimd(3,3) = primitive translations in recip space
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  mpsang=1+max. angular momentum
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
!! This term corresponds to term (1) of Eq. 49 in Torrent et al.,
!! CMS 42, 337 (2008), including phase shifts as needed for orbital
!! magnetization:
!! $\langle\Phi_i|e^{i\mathbf{b.r}}(-\frac{1}{2}\nabla^2)e^{-i\mathbf{k.r}}|\Phi_j\rangle -
!!  \langle\tilde{\Phi}_i|e^{i\mathbf{b.r}}(-\frac{1}{2}\nabla^2)e^{-i\mathbf{k.r}}|\tilde{\Phi}_j\rangle $
!! where $\mathbf{b}$ is the bra shift vector and $\mathbf{k}$ is the ket shift vector, both are determined in
!! initorbmag.F90
!!
!! PARENTS
!!      initorbmag
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,initylmr,jbessel,nderiv_gen,realgaunt
!!      simp_gen
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

 subroutine pawtwdij_1(dtbfield,gprimd,mpsang,natom,ntypat,pawrad,pawtab,typat,&
&                      mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_xmpi, only : xmpi_self

 use m_bfield
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_pawang, only : realgaunt
 use m_pawrad, only : pawrad_type, simp_gen, nderiv_gen
 use m_pawtab, only : pawtab_type
 use m_sphharm, only : initylmr
 use m_paw_numeric, only: jbessel

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtwdij_1'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: mpsang,natom,ntypat
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
 integer :: bdir,bsig,expibi_idx,iatom,iatom_tot,itypat,kdir,ksig,mu
 integer :: lcmax, lcfac, mesh_size,ir
 integer :: ilm1, ilm1_k, ilm2, ilmp1, ilmp1_k, cbl, cbm
 integer :: bl, bm, blm, bln, blmn, bgnt
 integer :: kl, km, klm, kln, klmn, kgnt
 integer :: klm1, klm1_k
 integer :: lb, mb, lk, mk, my_comm_atom, ngnt, twdx
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 logical :: my_atmtab_allocated,need_conjg,paral_atom
 real(dp) :: bnorm,knorm,intgrl
 real(dp) :: bessarg,jlx,jldx,jldx2

! real(dp) :: rpaw, testb, testc

!arrays
 integer,allocatable :: gntselect(:,:)
 logical,allocatable :: have_radint(:,:,:)
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: bb(3),bbn(3),bcart(3)
 real(dp) :: kb(3),kbn(3),kcart(3)
 real(dp) :: ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: realgnt(:)
 real(dp),allocatable :: ylmb(:),ylmk(:)
 real(dp),allocatable :: buj(:),btuj(:),tuk2(:,:),uk2(:,:)
 real(dp),allocatable :: kuj(:),ktuj(:),intgrnd(:),radint(:,:,:)
 complex(dpc) :: afac,kijsum
! the following is (i)^L mod 4.
 complex(dpc),dimension(0:3) :: ipl(0:3)=(/(1.0,0.0),(0.0,1.0),(-1.0,0.0),(0.0,-1.0)/)
! the following is (-i)^L mod 4.
 complex(dpc),dimension(0:3) :: iml(0:3)=(/(1.0,0.0),(0.0,-1.0),(-1.0,0.0),(0.0,1.0)/)

! *************************************************************************

!this should probably be a user input convergence parameter
 lcfac = 2 ! extra number of angular momentum states to use in sums beyond max in PAW sets

 lcmax = mpsang+lcfac ! lcmax - 1 is highest angular momentum state used in expansion

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr

!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.(dtbfield%my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_self;if (present(mpi_comm_atom)) my_comm_atom=mpi_comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=dtbfield%my_natom)

 ABI_ALLOCATE(ylmb,(lcmax*lcmax))
 ABI_ALLOCATE(ylmk,(lcmax*lcmax))
 ABI_ALLOCATE(gntselect,((2*lcmax-1)**2,lcmax**2*(lcmax**2+1)/2))
 ABI_ALLOCATE(realgnt,((2*lcmax-1)**2*(lcmax)**4))
 call realgaunt(lcmax,ngnt,gntselect,realgnt)

 ABI_ALLOCATE(radint,(0:lcmax-1,0:lcmax-1,0:2*lcmax))
 ABI_ALLOCATE(have_radint,(0:lcmax-1,0:lcmax-1,0:2*lcmax))

 do itypat = 1, ntypat

   mesh_size = pawrad(itypat)%mesh_size
   indlmn => pawtab(itypat)%indlmn
!  uncomment the following if you need the exact radius of the integrations
!  write(std_out,'(f12.8)')pawrad(itypat)%rad(pawrad(itypat)%int_meshsz)
!  rpaw = pawrad(itypat)%rad(pawrad(itypat)%int_meshsz)
!  testb = -2.0*(rpaw*rpaw-2.0)/rpaw**3
!  testc = (rpaw*rpaw - 2.0)/rpaw**4

   ABI_ALLOCATE(buj,(mesh_size))
   ABI_ALLOCATE(btuj,(mesh_size))
   ABI_ALLOCATE(tuk2,(mesh_size,2))
   ABI_ALLOCATE(uk2,(mesh_size,2))
   ABI_ALLOCATE(kuj,(mesh_size))
   ABI_ALLOCATE(ktuj,(mesh_size))
   ABI_ALLOCATE(intgrnd,(mesh_size))

   intgrnd(:) = zero

   do bdir = 1, 3
     do kdir = 1, 3
       if (kdir == bdir) cycle ! never need the kdir // bdir terms

       do bsig = -1, 1, 2
         bb(:) = bsig*dtbfield%dkvecs(:,bdir) ! bra vector

         do mu=1,3
           bcart(mu)=dot_product(bb(:),gprimd(mu,:))
         end do

!        JWZ override for testing
!        bcart(1) = 0.25d0; bcart(2) = 0.0d0; bcart(3) = 0.0d0
!        end override

         bnorm=dsqrt(dot_product(bcart,bcart))

         if (bnorm < tol12) then
           bbn(:) = zero
           ylmb(:) = zero; ylmb(1) = one/sqrt(four_pi)
         else
           bbn(:) = bcart(:)/bnorm ! unit vector in bb direction
           call initylmr(lcmax,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,bbn,ylmb(:),ylmgr)
         end if

         bnorm = two_pi*bnorm ! re-normed bb for calls to bessel fnc

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

!          JWZ override for testing
!          kcart(1)=0.0d0;kcart(2)=0.25d0;kcart(3)=0.00d0
!          end override

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

!              for testing pure k = 0 case
!              if ( (kl /= bl) .OR. (bm /= km)) cycle 
!              end test addition

!              get 1st and 2nd derivatives of ket-side u_nl functions
               call nderiv_gen(uk2,pawtab(itypat)%phi(:,kln),2,pawrad(itypat))
               call nderiv_gen(tuk2,pawtab(itypat)%tphi(:,kln),2,pawrad(itypat))

               kijsum = cmplx(zero,zero)
               radint(:,:,:) = zero
               have_radint(:,:,:) = .FALSE.

               do lb = 0, max(bl,kl)+lcfac ! loop over bra-side l expansion

!                make bra-side j_nl(b*r)*u_nili(r)
                 do ir = 1, mesh_size

                   bessarg = pawrad(itypat)%rad(ir)*bnorm
                   call jbessel(jlx,jldx,jldx2,lb,0,bessarg)
                   buj(ir) = jlx*pawtab(itypat)%phi(ir,bln)
                   btuj(ir) = jlx*pawtab(itypat)%tphi(ir,bln)

!                  JWZ override for testing
!                  buj(ir) = jlx*sin(pawrad(itypat)%rad(ir))
!                  radval = pawrad(itypat)%rad(ir)
!                  buj(ir) = jlx*(radval*(1.0+radval*(testb+testc*radval)))
!                  end override

                 end do

                 do lk = 0, max(bl,kl)+lcfac ! loop over ket side l expansion

                   do cbl = abs(lb-bl),lb+bl,2 ! only possible L values
!                    do cbl = 0, max((lb+bl),(lk+kl))

                     do ir = 2, mesh_size
                       bessarg = pawrad(itypat)%rad(ir)*knorm
                       call jbessel(jlx,jldx,jldx2,lk,2,bessarg)

                       kuj(ir)=knorm*knorm*jldx2*pawtab(itypat)%phi(ir,kln)+2.0D0*knorm*jldx*uk2(ir,1)+jlx*uk2(ir,2)
                       ktuj(ir)=knorm*knorm*jldx2*pawtab(itypat)%tphi(ir,kln)+2.0D0*knorm*jldx*tuk2(ir,1)+jlx*tuk2(ir,2)

!                      JWZ override for testing
!                      kuj(ir)=knorm*knorm*jldx2*sin(pawrad(itypat)%rad(ir))+&
!                      &  2.0d0*knorm*jldx*cos(pawrad(itypat)%rad(ir))-&
!                      &  jlx*sin(pawrad(itypat)%rad(ir))

!                      radval = pawrad(itypat)%rad(ir)
!                      kuj(ir) = jldx2*(radval*(1.0+radval*(testb+testc*radval)))+&
!                      &                        2.0d0*jldx*(1.0+radval*(2.0*testb+3.0*testc*radval))+&
!                      &                        jlx*(2.0*testb+6.0*testc*radval)
!                      end override
                       
                       if (cbl > 0) then
                         kuj(ir) = kuj(ir)-cbl*(cbl+1)*pawtab(itypat)%phi(ir,kln)*jlx/&
&                         (pawrad(itypat)%rad(ir)**2)
                         ktuj(ir) = ktuj(ir)-cbl*(cbl+1)*pawtab(itypat)%tphi(ir,kln)*jlx/&
&                         (pawrad(itypat)%rad(ir)**2)
!                        
!                        JWZ override for testing
!                        kuj(ir)=kuj(ir)-cbl*(cbl+1.0)*jlx*sin(pawrad(itypat)%rad(ir))/&
!                        &  (pawrad(itypat)%rad(ir)**2)
!                        kuj(ir)=kuj(ir)-cbl*(cbl+1.0)*jlx*(radval*(1.0+radval*(testb+testc*radval)))/&
!                        &                          radval**2
!                        end override
                       end if

                       intgrnd(ir) = buj(ir)*kuj(ir) - btuj(ir)*ktuj(ir)
!                      JWZ override for testing
!                      intgrnd(ir) = buj(ir)*kuj(ir)
!                      end override

                     end do ! end loop over mesh to make ket-side u function
!                    call simp_gen(intgrl,intgrnd,pawrad(itypat))
!                    radint(lb,lk,cbl) = intgrl

                     do mb = -lb, lb ! loop over bra expansion m values

                       ilmp1 = LMC(lb,mb)
                       ilm1  = LMC(bl,bm)
                       klm1 = RCC(ilmp1,ilm1)

                       do cbm = -cbl, cbl ! loop over M values
                         ilm2 = LMC(cbl,cbm)
                         bgnt = gntselect(ilm2,klm1)

                         if(bgnt == 0) cycle

                         do mk = -lk, lk ! loop over ket expansion m values

                           ilmp1_k = LMC(lk,mk)
                           ilm1_k = LMC(kl,km)
                           klm1_k = RCC(ilm1_k,ilmp1_k)

!                          note use of ilm2 in the following, it forces the LM state on the ket side to be identical
!                          to the LM state on the bra side, which is the result of the <LM_bra|LM_ket> integration
                           kgnt = gntselect(ilm2,klm1_k)

                           if (kgnt == 0) cycle

                           if(have_radint(lb,lk,cbl)) then
                             intgrl = radint(lb,lk,cbl)
                           else
                             call simp_gen(intgrl,intgrnd,pawrad(itypat))
                             radint(lb,lk,cbl) = intgrl
                             have_radint(lb,lk,cbl) = .TRUE.
                           end if

                           kijsum = kijsum + ipl(mod(lb,4))*iml(mod(lk,4))*&
&                           ylmb(ilmp1)*ylmk(ilmp1_k)*realgnt(bgnt)*realgnt(kgnt)*radint(lb,lk,cbl)

                         end do ! end loop over mk
                       end do ! end loop over cbm
                     end do ! end loop over mb
                   end do ! end loop over lk
                 end do ! end loop over cbl
               end do ! end loop over lb

               kijsum = -8.0*pi*pi*kijsum

               do iatom = 1,dtbfield%my_natom ! store result multiplied by exp(i*(k_b-k_k).I)
                 iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
                 if(typat(iatom_tot) == itypat) then
                   afac = cmplx(dtbfield%twexpibi(1,iatom,expibi_idx),dtbfield%twexpibi(2,iatom,expibi_idx))
                   if(need_conjg) afac = conjg(afac)
                   dtbfield%twdij0(1,blmn,klmn,iatom,twdx) = real(kijsum*afac)
                   dtbfield%twdij0(2,blmn,klmn,iatom,twdx) = aimag(kijsum*afac)
                 end if
               end do

!              write(std_out,'(a,24,2es16.8)')'JWZ Debug: blmn,klmn,kijsum',blmn,klmn,real(kijsum),aimag(kijsum)

             end do ! end loop over ket states
           end do ! end loop over bra states

         end do ! end loop over ksig
       end do ! end loop over kdir

     end do ! end loop over bsig
   end do ! end loop over bdir

   ABI_DEALLOCATE(buj)
   ABI_DEALLOCATE(btuj)
   ABI_DEALLOCATE(tuk2)
   ABI_DEALLOCATE(uk2)
   ABI_DEALLOCATE(kuj)
   ABI_DEALLOCATE(ktuj)
   ABI_DEALLOCATE(intgrnd)

 end do ! end loop on ntypat

 ABI_DEALLOCATE(gntselect)
 ABI_DEALLOCATE(realgnt)
 ABI_DEALLOCATE(radint)
 ABI_DEALLOCATE(have_radint)
 ABI_DEALLOCATE(ylmb)
 ABI_DEALLOCATE(ylmk)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 end subroutine pawtwdij_1
!!***
