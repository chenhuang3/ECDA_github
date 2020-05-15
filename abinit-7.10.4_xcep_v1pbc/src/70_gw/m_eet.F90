!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_eet
!! NAME
!! m_eet
!!
!! FUNCTION
!! Module containing the EET subroutines
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2014 ABINIT group (AB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


MODULE m_eet

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_gwdefs
 use m_errors 
 use m_profiling_abi
 use m_xomp
 use m_fft            
 use m_splines
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_ppmodel

 use m_io_tools,      only : get_unit
 use m_blas,          only : xgerc, xgemv, xgemm, xdotc
 use m_geometry,      only : normv, vdotw
 use m_vcoul,         only : vcoul_t
 use m_sigma,         only : sigma_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_screening,     only : epsilonm1_results

 implicit none

 private 
!!***

 public :: gw_eet_sigma         ! Wrapper routine for the calculation of the matrix elements of Sigma using the EET
 public :: gw_eet_chi0          ! Wrapper routine for the calculation of the polarizability using the EET
 public :: gw_eet_sigma_cd      ! Wrapper routine for the calculation of the matrix elements of Sigma 
                                ! using the EET and contour deformation

CONTAINS  !=========================================================================================================================
!!***

!!****f* m_eet/gw_eet_sigma
!! NAME
!! gw_eet_sigma
!!
!! FUNCTION
!! Wrapper routine for the calculation of the matrix elements of Sigma using the EET
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine gw_eet_sigma(Sigp,Sr,Dtset,Cryst,Wfs,Kmesh,Qmesh,Gsph_Max,Gsph_c,Psps,Vcp,QP_BSt,PPm, &
&                       isppol,iq_bz,ik_bz,jk_bz,ik_ibz,jk_ibz,itim_q,isym_q,iq_ibz,tabr_ki, &
&                       tabr_kj,spinrot_ki,spinrot_kj,ph_mkit,ph_mkjt,nfftot_gw,ngfft_gw, &
&                       use_padfft,igfftcg0,gw_gbound,gw_mgfft,ib1,ib2,nomega_tot,nomega_sigc, &
&                       fact_sp,nspinor,botsq,otq,sigcme_tmp,sigc,nbhomo,tim_fourdp,wtqp,wtqm, &
&                       extrapolar_distrb,can_symmetrize)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_eet_sigma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: isppol,itim_q,isym_q,iq_ibz,nspinor
 integer,intent(in) :: nomega_tot,nomega_sigc
 integer,intent(in) :: wtqp,wtqm
 integer,intent(in) :: iq_bz,ik_bz,jk_bz,ik_ibz,jk_ibz,ib1,ib2
 integer,intent(out) :: nbhomo
 real(dp),intent(in) :: fact_sp
 complex(dpc),intent(in) ::  ph_mkit,ph_mkjt
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_Max
 type(gsphere_t),intent(in) :: Gsph_c
 type(Dataset_type),intent(in) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Sigma_parameters),intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr
 type(wfd_t),target,intent(inout) :: Wfs
 type(ebands_t),target,intent(in) :: QP_BSt
 type(ppmodel_t),intent(in) :: PPm
!arrays
 integer,intent(in) :: igfftcg0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: tabr_ki(nfftot_gw),tabr_kj(nfftot_gw)
 integer,intent(in) :: extrapolar_distrb(ib1:ib2,ib1:ib2,Kmesh%nbz,Wfs%nsppol)
 real(dp),intent(in) :: spinrot_ki(4),spinrot_kj(4)
 complex(dpc),intent(inout) :: sigcme_tmp(nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nsppol*Sigp%nsig_ab)
 complex(dpc),intent(inout) :: sigc(2,nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nsppol*Sigp%nsig_ab)
 logical,intent(in) :: can_symmetrize(Wfs%nsppol)
 complex(gwpc),intent(in) :: otq(Sigp%npwc,PPm%dm2_otq)
 complex(gwpc),intent(in) :: botsq(Sigp%npwc,PPm%dm2_botsq)

!Local variables ------------------------------
!scalars
 integer :: io,kb,jb
 integer :: niter,nptwg,iter,nbmax
 integer :: ig,igp,ib,ibv
 integer :: isym_kgw,isym_ki,iik,jik
!arrays
 integer :: fnlloc(Cryst%ntypat,2),fnlmax(Cryst%ntypat)
 real(dp),pointer :: qp_ene(:,:,:),qp_occ(:,:,:)
 real(dp),allocatable :: qplg(:,:)
 real(dp),allocatable :: kplqg(:)
 real(dp),allocatable :: omegame0k(:),omegame0lumo(:)
 real(dp),allocatable :: qbzpg(:) 
 complex(gwpc),allocatable :: vc_sqrt_qbz(:)
 complex(gwpc),allocatable :: fnlkr(:,:,:)
 complex(gwpc),allocatable :: fnlkpr(:,:,:)
 complex(gwpc),allocatable :: wfr1(:,:)
 complex(gwpc),allocatable :: mtwk(:,:)
 complex(gwpc),allocatable :: mtwkp(:,:)
 complex(gwpc),allocatable :: ptwsq(:,:,:)
 complex(dpc) :: sigctmp(nomega_sigc)

!************************************************************************

 qp_ene => QP_BSt%eig(:,:,:)
 qp_occ => QP_BSt%occ(:,:,:)

 !@Arjan:
 ! Gsph_c gamma-centered sphere for W and thus Sigma_x
 ! Gsph_Max gamma-centered sphere with gvec(3,npwvec) where Sigp%npwvec=MAX(Sigp%npwwfn,Sigp%npwx)
 ! it is used for the wavefunctions but it will be REMOVED when we switch to k-centered
 ! G-spheres for the wavefunctions

 !call gsph_init(Gsph_Max,Cryst,Sigp%npwvec,gvec=gvec_kss)

 ! min and Max band indeces for GW corrections (for this k-point)

 nbhomo=1
 do ib = 2, Sigp%nbnds 
   if (fact_sp*qp_occ(ib,jk_ibz,isppol)<GW_TOL_DOCC) exit
   nbhomo=nbhomo+1
 enddo
 nbmax=max(nbhomo,Dtset%gw_eet_nband)

 niter = Dtset%gw_eet
 nptwg=1
 do iter = 1, niter
   nptwg=nptwg+mod(iter,2)
 enddo

 ABI_ALLOCATE(vc_sqrt_qbz,(Sigp%npwc))
 do ig=1,Sigp%npwc
   vc_sqrt_qbz(Gsph_c%rottb(ig,itim_q,isym_q))=Vcp%vc_sqrt(ig,iq_ibz)
 end do

 ABI_ALLOCATE(omegame0k,(nomega_tot))
 ABI_ALLOCATE(omegame0lumo,(nomega_tot))
 ABI_ALLOCATE(qplg,(3,Sigp%npwc))
 ABI_ALLOCATE(kplqg,(Sigp%npwc))
 ABI_ALLOCATE(qbzpg,(Sigp%npwc))

 isym_kgw = Kmesh%tabo(jk_bz)
 jik = (3-Kmesh%tabi(jk_bz))/2

 isym_ki = Kmesh%tabo(ik_bz)
 iik = (3-Kmesh%tabi(ik_bz))/2
!
!  Calculate some auxiliary quantities
!
 do ig=1,Sigp%npwc
   qplg(:,ig) = Qmesh%bz(:,iq_bz) + Gsph_c%gvec(:,ig)
   kplqg(ig) = -vdotw(Kmesh%bz(:,jk_bz),qplg(:,ig),Cryst%gmet,"G")
   qbzpg(ig) = normv(qplg(:,ig),Cryst%gmet,"G")
 end do

 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then
!
! Calculate FFTs of f_{nlk}(G) and f_{nlkp}(G) and the auxiliary functions 
! Mtw_k(r) = \sum_{nl} f_{nlk}(r) * \sum_{G} u_{vk}(G) f_{nlk}(G) and 
! Mtw_kp(r) = \sum_{nl} f_{nlkp}(r) * \sum_{G} u_{vkp}(G) f_{nlkp}(G)
!
   ABI_ALLOCATE(fnlkr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlkpr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(mtwk,(Wfs%nfftot*nspinor,nbmax))
   ABI_ALLOCATE(mtwkp,(Wfs%nfftot*nspinor,ib1:ib2))

   call gw_eet_sigma_vkb(Sigp,Cryst,Wfs,Kmesh,Psps,isppol,ik_ibz,jk_ibz,ib1,ib2,nspinor, &
&                        tim_fourdp,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)
 endif

 ABI_ALLOCATE(wfr1,(Wfs%nfftot*nspinor,nbmax))
!
! Get FFTs of the wave functions
!
 do ibv = 1, nbmax
   call wfd_get_ur(Wfs,ibv,ik_ibz,isppol,wfr1(:,ibv))
 enddo

 do jb = ib1,ib2
   do kb = ib1,ib2
     if (kb/=jb) CYCLE
     if (extrapolar_distrb(jb,kb,ik_bz,isppol)/=Wfs%my_rank) CYCLE
     ABI_ALLOCATE(ptwsq,(Sigp%npwc,Sigp%npwc,niter+1))
     do io=1,Sr%nomega4sd
       omegame0k(io)  = real(Sr%omega4sd(kb,jk_ibz,io,isppol)) - qp_ene(kb,jk_ibz,isppol)
       omegame0lumo(io)= real(Sr%omega4sd(kb,jk_ibz,io,isppol)) - qp_ene(nbmax+1,ik_ibz,isppol)
     end do
     sigctmp=czero_gw
     if (niter==0.or.Dtset%gw_eet_inclvkb==0) then
!
! Calculate f^{pp}, f^{pj}, and f{jj}. See Eqs. (25),(26),(27) of PRB 85, 085126, for their expressions.
! They are stored in ptwsq(:,:,1),ptwsq(:,:,2),ptwsq(:,:,3)
! Case: without vkb
!
       call fft4eet_sig(Sigp,Cryst,Wfs,Kmesh,Gsph_c,Sr,nbhomo,nbmax,nomega_tot,isppol,nfftot_gw, &
&                   ngfft_gw,use_padfft,igfftcg0,gw_gbound,gw_mgfft,iik,tabr_ki,ph_mkit,spinrot_ki, &
&                   ik_ibz,jk_ibz,isym_kgw,jik,tabr_kj,ph_mkjt,spinrot_kj,Gsph_Max, &
&                   nspinor,tim_fourdp,wfr1,vc_sqrt_qbz,Vcp%i_sz,kb,qplg,kplqg,niter, &
&                   ptwsq,ik_bz,jk_bz,PPm%dm2_botsq,PPm%dm2_otq,botsq,otq,sigctmp)

     else
!
! Case: with vkb
!
       call fft4eet_sig_kb(Sigp,Cryst,Wfs,Kmesh,Gsph_c,Psps,Sr,nbhomo,nbmax,nomega_tot,isppol,nfftot_gw, &
&                   ngfft_gw,use_padfft,igfftcg0,gw_gbound,gw_mgfft,iik,tabr_ki,ph_mkit,spinrot_ki, &
&                   ik_ibz,jk_ibz,isym_kgw,jik,tabr_kj,ph_mkjt,spinrot_kj,Gsph_Max, &
&                   nspinor,tim_fourdp,fnlloc,fnlmax,fnlkr,mtwk,mtwkp(:,kb),wfr1, &
&                   vc_sqrt_qbz,Vcp%i_sz,kb,qplg,kplqg,niter,ptwsq,ik_bz,jk_bz,PPm%dm2_botsq,PPm%dm2_otq,botsq,otq, &
&                   sigctmp)

     endif
!
! Calculation of the ptwsq(:,:,2)/ptwsq(:,:,1) and ptwsq(:,:,3)/ptwsq(:,:,2) and store in ptwsq(:,:,2) and ptwsq(:,:,3)
!
     call calc_delta_ppm(Sigp,ptwsq,niter)
!
! Divide ptwsq(:,:,1) by symmetrized Coulomb potential, i.e, |q+G||q+G'|
! Treat separately the case (q,G,Gp)=(0,0,0)
!
     do ig = 1, Sigp%npwc
       do igp = 1, Sigp%npwc
         if (ik_bz==jk_bz) then
           if (ig/=1.and.igp/=1) then
             ptwsq(ig,igp,1)= ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)*vc_sqrt_qbz(igp)
           else
             if (jb/=kb.or.kb<=nbmax.or.jb<=nbmax) then
               ptwsq(ig,igp,1)=(0.0,0.0)
             else
               if (ig==1.and.igp==1) then
                 ptwsq(ig,igp,1) = cmplx(Vcp%i_sz,0.0_gwp)
               elseif (ig==1.and.igp/=1) then
                 ptwsq(ig,igp,1) = cmplx(sqrt(Vcp%i_sz),0.0_gwp)*ptwsq(ig,igp,1)*vc_sqrt_qbz(igp)
               elseif (igp==1.and.ig/=1) then
                 ptwsq(ig,igp,1) = cmplx(sqrt(Vcp%i_sz),0.0_gwp)*ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)
               else
                 ptwsq(ig,igp,1)= ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)*vc_sqrt_qbz(igp)
               endif
             endif
           endif
         else
           ptwsq(ig,igp,1)= ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)*vc_sqrt_qbz(igp)
         endif
       enddo
     enddo
!
! Calculate the matrix elements of the self-energy (sigctmp) using the closure relation
!
     call calc_sig_ppm_delta_clos(Sigp%npwc,nomega_tot,ik_bz,jk_bz,qbzpg,botsq,otq,omegame0k,omegame0lumo,Sigp%zcut, &
&                                 ptwsq,sigctmp,PPm%dm2_botsq,PPm%dm2_otq,Dtset%gw_eet_scale,niter)
!
! Add sigctmp to sigcme_tmp with or without using symmetry
!
     if (can_symmetrize(isppol)) then
       sigcme_tmp(:,jb,kb,isppol)=sigcme_tmp(:,jb,kb,isppol) + &
         (wtqp+wtqm)*DBLE(sigctmp(:)) + (wtqp-wtqm)*j_gw*AIMAG(sigctmp(:))
       sigc(1,:,jb,kb,isppol)=sigc(1,:,jb,kb,isppol) + wtqp*      sigctmp(:)
       sigc(2,:,jb,kb,isppol)=sigc(2,:,jb,kb,isppol) + wtqm*CONJG(sigctmp(:))
     else
       sigcme_tmp(:,jb,kb,isppol)=sigcme_tmp(:,jb,kb,isppol)+sigctmp(:)
     endif

     ABI_DEALLOCATE(ptwsq)

   enddo
 enddo

 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then
   ABI_DEALLOCATE(mtwk)
   ABI_DEALLOCATE(mtwkp)
   ABI_DEALLOCATE(fnlkr)
   ABI_DEALLOCATE(fnlkpr)
 endif

 ABI_DEALLOCATE(wfr1)
 ABI_DEALLOCATE(vc_sqrt_qbz)
 ABI_DEALLOCATE(omegame0k)
 ABI_DEALLOCATE(omegame0lumo)
 ABI_DEALLOCATE(qplg)
 ABI_DEALLOCATE(kplqg)
 ABI_DEALLOCATE(qbzpg)

end subroutine gw_eet_sigma
!!***

!!****f* m_eet/gw_eet_sigma_vkb
!! NAME
!! gw_eet_sigma_vkb
!!
!! FUNCTION
!!      Calculate FFTs of f_{nlk}(G) and f_{nlkp}(G) and the auxiliary functions 
!!      Mtw_k(r) = \sum_{nl} f_{nlk}(r) * \sum_{G} u_{vk}(G) f_{nlk}(G) and 
!!      Mtw_kp(r) = \sum_{nl} f_{nlkp}(r) * \sum_{G} u_{vkp}(G) f_{nlkp}(G)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine gw_eet_sigma_vkb(Sigp,Cryst,Wfs,Kmesh,Psps,isppol,ik_ibz,jk_ibz,ib1,ib2,nspinor, &
&                           tim_fourdp,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_eet_sigma_vkb'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: tim_fourdp
 integer,intent(in) :: isppol
 integer,intent(in) :: nspinor
 integer,intent(in) :: ik_ibz,jk_ibz,ib1,ib2
 integer,intent(in) :: nbmax
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Sigma_parameters),intent(in) :: Sigp
 type(wfd_t),target,intent(inout) :: Wfs
 integer,intent(out) :: fnlloc(Cryst%ntypat,2),fnlmax(Cryst%ntypat)
 complex(gwpc),intent(out) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(out) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,ib1:ib2)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: istwf_ki,istwf_kj,npw_ki,npw_kj,iat,ilm,i
 character(len=fnlen) :: title
 integer :: lloc,lmax,mmax,pspcod,pspdat,pspxc,temp_unit
 integer :: ityp
 real(dp) :: r2well,zion,znucl
 type(kb_potential) :: KBff_ki,KBff_kj
!arrays
 integer,pointer :: gbound_ki(:,:),gbound_kj(:,:)
 integer,pointer :: kg_ki(:,:),kg_kj(:,:)
 complex(gwpc),allocatable :: sfnl_tmp(:) 

!************************************************************************

 ABI_UNUSED(tim_fourdp)

 fnlloc(:,:)=0
 fnlmax(:)=0
 do ityp = 1, Cryst%ntypat
   temp_unit = get_unit()
   open (unit=temp_unit,file=psps%filpsp(ityp),form='formatted',status='old')
   rewind (unit=temp_unit)
   read (temp_unit,'(a)') title
   read (temp_unit,*) znucl,zion,pspdat
   read (temp_unit,*) pspcod,pspxc,lmax,lloc,mmax,r2well
   close(temp_unit)
   do i = 1, lloc
     fnlloc(ityp,1) = fnlloc(ityp,1) + 2*(lloc-1)+1
   enddo
   fnlloc(ityp,1)=fnlloc(ityp,1)+1
   do i = 1, lloc+1
     fnlloc(ityp,2) = fnlloc(ityp,2) + 2*(lloc-1)+1
   enddo
   do i = 1, lmax+1
     fnlmax(ityp) = fnlmax(ityp) + 2*(lmax-1)+1
   enddo
 enddo

 istwf_ki  =  Wfs%istwfk(ik_ibz)
 npw_ki    =  Wfs%npwarr(ik_ibz)
 kg_ki     => Wfs%Kdata(ik_ibz)%kg_k
 gbound_ki => Wfs%Kdata(ik_ibz)%gbound
 
 istwf_kj  =  Wfs%istwfk(jk_ibz)
 npw_kj    =  Wfs%npwarr(jk_ibz)
 kg_kj     => Wfs%Kdata(jk_ibz)%kg_k
 gbound_kj => Wfs%Kdata(jk_ibz)%gbound

 !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
 call init_kb_potential(KBff_ki,Cryst,Psps,2,istwf_ki,npw_ki,Kmesh%ibz(:,ik_ibz),kg_ki)
 call init_kb_potential(KBff_kj,Cryst,Psps,2,istwf_kj,npw_kj,Kmesh%ibz(:,jk_ibz),kg_kj)
 ABI_DEALLOCATE(KBff_ki%fnld)
 ABI_DEALLOCATE(KBff_kj%fnld)

 !MG: note that it is not necessary to rotate fnl in the full BZ to
 !calculate
 !<SK|V_nl|Sk>.
 !MG check this part in parallel. Keep in mind the difference between
 !Wfs_braket (bdgw states= and Wfs (full set distributed
 !across the nod

 ! MG: TODO Be careful here as the symmetry properties of fnl are not the same as the 
 !     ones of the wavefunctions. One should write a separate methods for the FFT of fnl.
 if (istwf_kj>1.or.istwf_ki>1) then
   MSG_ERROR("istwfk /= 1 not coded")
 end if

 ABI_MALLOC(sfnl_tmp,(MAX(npw_ki,npw_kj)) )
!
! Calculate FFTs of fnl_k(G) and fnl_kp(G)
!
 do iat = 1, Cryst%natom
   do ilm = 1, Psps%mpsang*Psps%mpsang
     if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
     if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE

     sfnl_tmp(1:npw_ki) = conjg(KBff_ki%fnl(:,ilm,iat))

     call fft_ug(npw_ki,Wfs%nfftot,nspinor,ndat1,Wfs%mgfft,Wfs%ngfft,istwf_ki,kg_ki,gbound_ki,sfnl_tmp,fnlkr(:,ilm,iat))

     sfnl_tmp(1:npw_kj) = conjg(KBff_kj%fnl(:,ilm,iat))

     call fft_ug(npw_kj,Wfs%nfftot,nspinor,ndat1,Wfs%mgfft,Wfs%ngfft,istwf_kj,kg_kj,gbound_kj,sfnl_tmp,fnlkpr(:,ilm,iat))

   enddo
 enddo

 ABI_FREE(sfnl_tmp)
 call destroy_kb_potential(KBff_ki)
 call destroy_kb_potential(KBff_kj)
!
! Calculate auxiliary functions Mtw_k(r) = \sum_{nl} fnl_k(r) * \sum_{G} u_{vk}(G) fnl_k(G)
! and Mtw_kp(r) = \sum_{nl} fnl_kp(r) * \sum_{G} u_{vkp}(G) fnl_kp(G)
!
 call calc_eet_sig_prep(Sigp,Cryst,Wfs,Kmesh,Psps,isppol,nbmax,ib1,ib2,ik_ibz,jk_ibz, &
&                       nspinor,fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)

end subroutine gw_eet_sigma_vkb
!!***

!!****f* m_eet/calc_eet_sig_prep
!! NAME                  
!! calc_eet_sig_prep
!!
!! FUNCTION
!!      Calculate auxiliary functions Mtw_k(r) = \sum_{nl} f_{nlk}(r) * \sum_{G} u_{vk}(G) f_{nlk}(G)
!!      and Mtw_kp(r) = \sum_{nl} f_{nlkp}(r) * \sum_{G} u_{vkp}(G) f_{nlkp}(G)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_eet_sig_prep(Sigp,Cryst,Wfs,Kmesh,Psps,is,nbmax,ib1,ib2,ik_ibz, &
&                            jk_ibz,nspinor,fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_eet_sig_prep'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),target,intent(in) :: Wfs
 type(Sigma_parameters),intent(in) :: Sigp
 type(kb_potential) :: KBff_ki,KBff_kj

 integer,intent(in) :: is,nbmax,ib1,ib2
 integer,intent(in) :: nspinor
 integer,intent(in) :: ik_ibz,jk_ibz
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)

 complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(in) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)

 complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,ib1:ib2)

!Local variables-------------------------------
!scalars
 integer :: ibv,kb,ilm,iat,ig
 integer :: istwf_ki,istwf_kj,npw_ki,npw_kj
!arrays 
 integer,pointer :: kg_ki(:,:),kg_kj(:,:)
 complex(gwpc),allocatable :: maux(:,:)

!************************************************************************

 !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
 istwf_ki = Wfs%istwfk(ik_ibz)
 istwf_kj = Wfs%istwfk(jk_ibz)

 npw_ki = Wfs%npwarr(ik_ibz)
 npw_kj = Wfs%npwarr(jk_ibz)

 kg_ki => Wfs%Kdata(ik_ibz)%kg_k
 kg_kj => Wfs%Kdata(jk_ibz)%kg_k

 call init_kb_potential(KBff_ki,Cryst,Psps,2,istwf_ki,npw_ki,Kmesh%ibz(:,ik_ibz),kg_ki)
 call init_kb_potential(KBff_kj,Cryst,Psps,2,istwf_kj,npw_kj,Kmesh%ibz(:,jk_ibz),kg_kj)
 ABI_DEALLOCATE(KBff_ki%fnld)
 ABI_DEALLOCATE(KBff_kj%fnld)

 ABI_ALLOCATE(maux,(Psps%mpsang*Psps%mpsang,Cryst%natom))
!
! Calculate auxiliary function Mtw_k(r) = \sum_{nl} fnl_k(r) * \sum_{G} u_{vk}(G) fnl_k(G)
!
 mtwk(:,:)=(0.0,0.0)
 do ibv = 1, nbmax
   maux(:,:)=(0.0,0)
   do ig = 1, npw_ki
     maux(:,:) = maux(:,:) + Wfs%Wave(ibv,ik_ibz,is)%ug(ig)*KBff_ki%fnl(ig,:,:)
   enddo
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       mtwk(:,ibv)=mtwk(:,ibv)+maux(ilm,iat)*fnlkr(:,ilm,iat)
     enddo
   enddo
 enddo
!
! Calculate auxiliary function Mtw_kp(r) = \sum_{nl} fnl_kp(r) * \sum_{G} u_{vkp}(G) fnl_kp(G)
!
 mtwkp(:,:)=(0.0,0.0)
 do kb = ib1, ib2
   maux(:,:)=(0.0,0)
   do ig = 1, npw_kj
     maux(:,:) = maux(:,:) + Wfs%Wave(kb,jk_ibz,is)%ug(ig)*KBff_kj%fnl(ig,:,:)
   enddo
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       mtwkp(:,kb)=mtwkp(:,kb)+maux(ilm,iat)*fnlkpr(:,ilm,iat)
     enddo
   enddo
 enddo

 ABI_DEALLOCATE(maux)

 call destroy_kb_potential(KBff_ki)
 call destroy_kb_potential(KBff_kj)

 RETURN 
 ABI_UNUSED(Sigp%npwc)

end subroutine calc_eet_sig_prep
!!***

!!****f* m_eet/fft4eet_sig
!! NAME                  
!! fft4eet_sig
!!
!! FUNCTION
!! Calculate f^{pp}, f^{pj}, and f{jj}. See Eqs. (25),(26),(27) of PRB 85, 085126, for their expressions.
!! They are stored in ptwsq(:,:,1),ptwsq(:,:,2),ptwsq(:,:,3). 
!! Terms involving the nonlocal part of the pseudopotential are not included.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine fft4eet_sig(Sigp,Cryst,Wfs,Kmesh,Gsph_c,Sr,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw, &
&                  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,Gsph_max, &
&                  nspinor,tim_fourdp,wfr1,vc_sqrt_qbz,i_sz,kb,qplg,kplqg,niter, &
&                  ptwsq,ik_bz,ikmq_bz,npwc1,npwc2,botsq,otq,sigmac)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_sig'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(Sigma_parameters),intent(in) :: Sigp
 type(wfd_t),target,intent(inout) :: Wfs
 type(sigma_t),intent(in) :: Sr
 type(gsphere_t),intent(in) :: Gsph_c,Gsph_max
 type(kmesh_t),intent(in) :: Kmesh
 integer,intent(in) :: nbhomo,nbmax,is,kb,niter,npwc1,npwc2,nomega
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_kmq,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qplg(3,Sigp%npwc),kplqg(Sigp%npwc)
 real(dp),intent(in) :: i_sz
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt
!arrays
 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)
 complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
 complex(dpc),intent(inout) :: sigmac(nomega)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: istwf_kmq,nbz,npw_kmq,ibv,ig,igp,i,j
 integer :: ig4,ig4x,ig4y,ig4z,outofbox
 integer,save :: enough=0
 complex(gwpc) :: minusone
 character(len=500) :: msg
!arrays
 integer :: gmgp(3)
 integer,pointer :: gbound_kmq(:,:),kg_kmq(:,:),igfft_kmq(:)
 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux(3)
 complex(gwpc),allocatable :: wfr2(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)
 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)

!************************************************************************

 ABI_UNUSED((/isym_kmq,tim_fourdp/))

!Dummy statement, to keep Kmesh arg
 nbz=Kmesh%nbz

 gbound_kmq => Wfs%Kdata(ikmq_ibz)%gbound
 igfft_kmq  => Wfs%Kdata(ikmq_ibz)%igfft0
 kg_kmq     => Wfs%Kdata(ikmq_ibz)%kg_k
 istwf_kmq  = Wfs%istwfk(ikmq_ibz)
 npw_kmq    = Wfs%npwarr(ikmq_ibz)

 ABI_ALLOCATE(wfr2,(Wfs%nfftot*nspinor))

 ABI_ALLOCATE(rhotwg,(Sigp%npwc*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Sigp%npwc*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(npw_kmq))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
 endif

 if (niter>1) then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
 endif

 call wfd_get_ur(Wfs,kb,ikmq_ibz,is,wfr2)
!
! Calculate matrix element <wfn2|exp[-iG.r]|wfn2>
!
 call calc_wfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,wfr2,wfr2,wfwfg(:))

 do ibv = 1, nbmax
!
! Calculate matrix element <wfn1|exp[-i(q+G).r]|wfn2>
!
   call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                spinrot_kmq,nspinor,rhotwg(:,ibv))

   if (ibv>nbhomo) then
!
! Add sum-over-states correction to the matrix element of the GW self-energy obtained with the EET
!
     call calc_corr_sig(Sigp,Sr,nomega,nspinor,npwc1,npwc2,botsq,otq,rhotwg(:,ibv), &
&                       is,ibv,kb,ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,i_sz,vc_sqrt_qbz,sigmac)
   endif
 enddo

 if (niter>0) then
!
! Calculate FFT of i nabla u(r)
!
  ! MG: FIXME: this won't work if k-centered G-spheres are used. We need ISKg
  call wfd_dur_isk(Wfs,Cryst,Kmesh,kb,ikmq_bz,is,npw_kmq,Gsph_max%rottbm1,dwfr)  !,ISkg,ik_ibz
!
! Calculate matrix element <wfn1|exp[-i(q+G).r]|nabla wfn2>
!
   do ibv = 1, nbmax
     do i = 1, 3
       call drho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibv),itim_k,tabr_k,ph_mkt,dwfr(:,i), &
&                     nspinor,drhotwg(:,ibv,i))
     enddo
   enddo
!
! Calculate matrix element <wfn2|exp[-iG.r]|nabla wfn2> and <nabla wfn2|exp[-iG.r]|nabla wfn2>
!
   do i = 1, 3
     call calc_dwfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt,wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(itim_kmq,nfftot_gw,ngfft_gw,dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo

   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Sigp%npwc,nbmax))
!
! Calculate (q+G).nabla u(G)
!
   cauxg(:,:)=(0.0,0.0)
   do ibv = 1, nbmax
     do ig = 1, Sigp%npwc
       drhaux(:)=cmplx(real(drhotwg(ig,ibv,:)),aimag(drhotwg(ig,ibv,:)))
       cauxg(ig,ibv)=vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
     enddo
   enddo

 endif
 ABI_DEALLOCATE(wfr2)
!
! Calculate <wfn2|exp[i(G-Gp).r]|wfn2>, <wfn2|exp[i(G-Gp).r]|nabla wfn2>.(q+Gp), and
! (q+G).<nabla wfn2|exp[i(G-Gp).r]|nabla wfn2>.(q+Gp)
! They are stored in ptwsq(:,:,1),ptwsq(:,:,2), and ptwsq(:,:,3)
!
 ptwsq(:,:,:)=(0.0,0.0)
 outofbox=0
 do igp=1,Sigp%npwc
   do ig=1,Sigp%npwc
     gmgp(:)=Gsph_c%gvec(:,igp)-Gsph_c%gvec(:,ig)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)
     if (igp>=ig) then
       ptwsq(ig,igp,1)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2) + vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")
       if (niter>1.and.igp>=ig) then
         paux(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux(i)=paux(i) + vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
         enddo
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) + vdotw(qplg(:,igp),paux,Cryst%gmet,"G")
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if
!
! Calculate sum_v <wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|wfn2> and subtract from ptwsq(:,:,1)
!
 do ibv = 1, nbmax
   do igp=1,Sigp%npwc
     do ig=1,igp
       ptwsq(ig,igp,1)=ptwsq(ig,igp,1)-conjg(rhotwg(ig,ibv))*rhotwg(igp,ibv)
     enddo
   end do
 end do

 do igp = 1, Sigp%npwc
   do ig = igp+1, Sigp%npwc
     ptwsq(ig,igp,1)=conjg(ptwsq(igp,ig,1))
   enddo
 enddo
!
! Calculate sum_v <wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|nabla wfn2>.(q+Gp) 
! and subtract from ptwsq(:,:,2)
!
 if (niter>0) then
   minusone=(-1.,0.)
   do ibv = 1, nbmax
     call XGERC(Sigp%npwc,Sigp%npwc,minusone,conjg(rhotwg(:,ibv)),1,conjg(cauxg(:,ibv)),1,ptwsq(:,:,2),Sigp%npwc)
   enddo
   do igp=1,Sigp%npwc
     do ig=1,Sigp%npwc
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig
 endif
!
! Calculate sum_v (q+G).<nabla wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|nabla wfn2>.(q+Gp) 
! and subtract from ptwsq(:,:,3)
!
 if (niter>1) then
   do ibv = 1, nbmax
     do igp=1,Sigp%npwc
       do ig=1,igp
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3)-conjg(cauxg(ig,ibv))*cauxg(igp,ibv)
       enddo
     enddo
   enddo
   do igp=1,Sigp%npwc
     do ig=1,igp
       ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+kplqg(igp)*ptwsq(ig,igp,2)+kplqg(ig)*conjg(ptwsq(igp,ig,2))&
&                                     -kplqg(ig)*ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig
   do igp = 1, Sigp%npwc
     do ig = igp+1, Sigp%npwc
       ptwsq(ig,igp,3)=conjg(ptwsq(igp,ig,3))
     enddo
   enddo
 endif

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
 endif
 if (niter>1) then
   ABI_DEALLOCATE(ddwfwfg)
 endif

end subroutine fft4eet_sig

!!***
!----------------------------------------------------------------------

!!****f* m_eet/fft4eet_sig_kb
!! NAME                  
!! fft4eet_sig_kb
!!
!! FUNCTION
!! Calculate f^{pp}, f^{pj}, and f{jj}. See Eqs. (25),(26),(27) of PRB 85, 085126, for their expressions.
!! They are stored in ptwsq(:,:,1),ptwsq(:,:,2),ptwsq(:,:,3). 
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine fft4eet_sig_kb(Sigp,Cryst,Wfs,Kmesh,Gsph_c,Psps,Sr,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw, &
&                  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,Gsph_max, &
&                  nspinor,tim_fourdp,fnlloc,fnlmax,fnlkr,mtwk,mtwkp,wfr1, &
&                  vc_sqrt_qbz,i_sz,kb,qplg,kplqg,niter,ptwsq,ik_bz,ikmq_bz, &
&                  npwc1,npwc2,botsq,otq,sigmac)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_sig_kb'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Sigma_parameters),intent(in) :: Sigp
 type(wfd_t),target,intent(inout) :: Wfs
 type(sigma_t),intent(in) :: Sr
 type(gsphere_t),intent(in) :: Gsph_c, Gsph_max

 integer,intent(in) :: nbhomo,nbmax,is,kb,niter,npwc1,npwc2,nomega
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_kmq,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qplg(3,Sigp%npwc),kplqg(Sigp%npwc)
 real(dp),intent(in) :: i_sz
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt

 complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(in) :: mtwk(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: mtwkp(Wfs%nfftot*nspinor)

 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)

 complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)

 complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
 complex(dpc),intent(inout) :: sigmac(nomega)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: istwf_kmq,npw_kmq,i,j,ibv,ilm,iat,ig,igp 
 integer :: ig4,ig4x,ig4y,ig4z
 integer :: ig5,ig5x,ig5y,ig5z,nlx,outofbox
 integer,save :: enough=0
 complex(gwpc) :: minusone
 character(len=500) :: msg
!arrays
 integer,pointer :: gbound_kmq(:,:),kg_kmq(:,:),igfft_kmq(:)
 integer :: gmgp(3)
 complex(gwpc),allocatable :: wfr2(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg2(:,:)
 complex(gwpc),allocatable :: fnltwg3(:,:)
 complex(gwpc),allocatable :: kns(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)
 complex(gwpc),allocatable :: fnlwfg(:)
 complex(gwpc),allocatable :: fkdwfg(:,:)
 complex(gwpc),allocatable :: fdrhotwg(:,:,:,:)
 complex(gwpc),allocatable :: lnkp(:)
 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)
 complex(gwpc),allocatable :: ff(:,:,:,:)
 complex(gwpc),allocatable :: vzn(:,:,:)
 complex(gwpc),allocatable :: paux(:,:)
 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux2(3)

!************************************************************************

 ABI_UNUSED((/isym_kmq,tim_fourdp/))

 gbound_kmq => Wfs%Kdata(ikmq_ibz)%gbound
 igfft_kmq  => Wfs%Kdata(ikmq_ibz)%igfft0
 kg_kmq     => Wfs%Kdata(ikmq_ibz)%kg_k
 istwf_kmq  = Wfs%istwfk(ikmq_ibz)
 npw_kmq    = Wfs%npwarr(ikmq_ibz)

 nlx = min(Psps%mpsang,4)

 ABI_ALLOCATE(wfr2,(Wfs%nfftot*nspinor))
 ABI_ALLOCATE(rhotwg,(Sigp%npwc*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Sigp%npwc*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(npw_kmq))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
   ABI_ALLOCATE(fnltwg,(Sigp%npwc*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlwfg,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(fnltwg2,(Sigp%npwc*nspinor**2,nbmax))
   ABI_ALLOCATE(fnltwg3,(Sigp%npwc*nspinor**2,nbmax))
 endif

 if (niter>1) then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
   ABI_ALLOCATE(lnkp,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(kns,(Sigp%npwc*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fdrhotwg,(Sigp%npwc*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom,3))
   ABI_ALLOCATE(fkdwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(ff,(nlx*nlx,Cryst%natom,nlx*nlx,Cryst%natom))
   ABI_ALLOCATE(vzn,(Sigp%npwc*nspinor**2,nlx*nlx,Cryst%natom))
 endif

 call wfd_get_ur(Wfs,kb,ikmq_ibz,is,wfr2)
!
! Calculate matrix element <wfn2|exp[-iG.r]|wfn2>
!
 call calc_wfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,wfr2,wfr2,wfwfg)

 do ibv = 1, nbmax
!
! Calculate matrix element <wfn1|exp[-i(q+G).r]|wfn2>
!
   call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                spinrot_kmq,nspinor,rhotwg(:,ibv))

   if (ibv>nbhomo) then
!
! Add sum-over-states correction to the matrix element of the GW self-energy obtained with the EET
!
     call calc_corr_sig(Sigp,Sr,nomega,nspinor,npwc1,npwc2,botsq,otq,rhotwg(:,ibv), &
&                       is,ibv,kb,ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,i_sz,vc_sqrt_qbz,sigmac)
   endif
 enddo

 if (niter>0) then
!
! Calculate matrix element <fnl|exp[-i(q+G).r]|wfn2>
!
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                    fnlkr(:,ilm,iat),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                    spinrot_kmq,nspinor,fnltwg(:,ilm,iat))
     enddo
   enddo
!
! Calculate matrix element <wfn2|exp[-iG.r]|Mtw>
!
   call calc_wfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,wfr2,mtwkp,fnlwfg)
!
! Calculate FFT of i nabla u(r)
!
  ! FIXME: this won't work if k-centered G-spheres are used, we need ISkg
   call wfd_dur_isk(Wfs,Cryst,Kmesh,kb,ikmq_bz,is,npw_kmq,Gsph_max%rottbm1,dwfr)  !,ISkg,ik_ibz
!
! Calculate matrix elements <wfn1|exp[-i(q+G).r]|nabla wfn2>, <wfn1|exp[-i(q+G).r]|Mtw_kp>, and
! <Mtw_k|exp[-i(q+G).r]|wfn2>
!
   do ibv = 1, nbmax
     do i = 1, 3
       call drho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibv),itim_k,tabr_k,ph_mkt,dwfr(:,i), &
&                     nspinor,drhotwg(:,ibv,i))
     enddo
     call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  wfr1(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,mtwkp,itim_kmq,tabr_kmq,ph_mkmqt, &
&                  spinrot_kmq,nspinor,fnltwg2(:,ibv))

     call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  mtwk(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                  spinrot_kmq,nspinor,fnltwg3(:,ibv))
   enddo
!
! Calculate matrix element <wfn2|exp[-i(q+G).r]|nabla wfn2> and <nabla wfn2|exp[-iG.r]|nabla wfn2>
!
   do i = 1, 3
     call calc_dwfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt,wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(itim_kmq,nfftot_gw,ngfft_gw,dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo
!
! Calculate matrix elements <fnl|exp[-i(q+G).r]|nabla wfn2> ,<fnl|exp[-i(q+G).r]|Mtw>, 
! <fnl|exp[-iG.r]|Mtw> and <Mtw|exp[-iG.r]|Mtw>
!
   if (niter>1) then
     do iat = 1, Cryst%natom
       do ilm = 1, Psps%mpsang*Psps%mpsang
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         do i = 1, 3
           call drho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                         fnlkr(:,ilm,iat),itim_k,tabr_k,ph_mkt,dwfr(:,i), &
&                         nspinor,fdrhotwg(:,ilm,iat,i))
         enddo
         call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                      fnlkr(:,ilm,iat),itim_k,tabr_k,ph_mkt,spinrot_k,mtwkp,itim_kmq,tabr_kmq,ph_mkmqt, &
&                      spinrot_kmq,nspinor,kns(:,ilm,iat))
       enddo
     enddo
     do i = 1, 3
       call calc_dwfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt,mtwkp,dwfr(:,i),fkdwfg(:,i))
     enddo
     call calc_wfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,mtwkp,mtwkp,lnkp)
   endif

   ABI_DEALLOCATE(wfr2)
   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Sigp%npwc,nbmax))
!
! Calculate (q+G).nabla u(G)-fnltwg2+fnltwg3
!
   cauxg(:,:)=(0.0,0.0)
   do ibv = 1, nbmax
     do ig=1,Sigp%npwc
       drhaux(:)=cmplx(real(drhotwg(ig,ibv,:)),aimag(drhotwg(ig,ibv,:)))
       cauxg(ig,ibv)=vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
     enddo
   enddo
   do ig=1,Sigp%npwc
     cauxg(ig,:)=cauxg(ig,:)-fnltwg2(ig,:)+fnltwg3(ig,:)
   enddo

 endif
!
! Calculate <wfn2|exp[i(G-Gp).r]|wfn2>, <wfn2|exp[i(G-Gp).r]|nabla wfn2>.(q+Gp)-fnlwfg, and
! (q+G).<nabla wfn2|exp[i(G-Gp).r]|nabla wfn2>.(q+Gp)-(q+G).fkdwfg+lnkp
! They are stored in ptwsq(:,:,1),ptwsq(:,:,2), and ptwsq(:,:,3)
!
 ptwsq(:,:,:)=(0.0,0.0)
 outofbox=0
 do igp=1,Sigp%npwc
   do ig=1,Sigp%npwc
     gmgp(:)=Gsph_c%gvec(:,igp)-Gsph_c%gvec(:,ig)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)

     ig5x= modulo(-gmgp(1),ngfft_gw(1))
     ig5y= modulo(-gmgp(2),ngfft_gw(2))
     ig5z= modulo(-gmgp(3),ngfft_gw(3))
     ig5= 1+ig5x+ig5y*ngfft_gw(1)+ig5z*ngfft_gw(1)*ngfft_gw(2)

     if (igp>=ig) then
       ptwsq(ig,igp,1)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2) + vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")-fnlwfg(ig4)
       if (niter>1.and.igp>=ig) then
         paux2(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux2(i)=paux2(i) + vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
         enddo
         drhaux(:)=paux2(:)-cmplx(real(fkdwfg(ig4,:)),aimag(fkdwfg(ig4,:)))
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) + vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")
         drhaux(:)=cmplx(real(fkdwfg(ig5,:)),-aimag(fkdwfg(ig5,:)))
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) - vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")+lnkp(ig4)
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if
!
! Calculate sum_v <wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|wfn2> and subtract from ptwsq(:,:,1)
!
 do ibv = 1, nbmax
   do igp=1,Sigp%npwc
     do ig=1,igp
       ptwsq(ig,igp,1)=ptwsq(ig,igp,1)-conjg(rhotwg(ig,ibv))*rhotwg(igp,ibv)
     enddo
   end do !igp
 end do !ig

 do igp = 1, Sigp%npwc
   do ig = igp+1, Sigp%npwc
     ptwsq(ig,igp,1)=conjg(ptwsq(igp,ig,1))
   enddo
 enddo
!
! Calculate sum_v <wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|nabla wfn2>.(q+Gp) + fnltw(G)^* fnltw(Gp)
! and subtract from ptwsq(:,:,2)
!
 if (niter>0) then
   ABI_ALLOCATE(paux,(Sigp%npwc,Sigp%npwc))
   paux(:,:)=(0.0,0.0)
   minusone=(-1.,0.)
   do ibv = 1, nbmax
     call XGERC(Sigp%npwc,Sigp%npwc,minusone,conjg(rhotwg(:,ibv)),1,conjg(cauxg(:,ibv)),1,ptwsq(:,:,2),Sigp%npwc)
   enddo
   do igp=1,Sigp%npwc
     do ig=1,Sigp%npwc
       if (ig>=igp) then
         do iat = 1, Cryst%natom
           do ilm = 1, nlx*nlx
             if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
             if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
             paux(ig,igp)=paux(ig,igp)+conjg(fnltwg(ig,ilm,iat))*fnltwg(igp,ilm,iat)
           enddo
         enddo
       endif
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig
   do ig=1,Sigp%npwc
     do igp=1,Sigp%npwc
       if (ig>=igp) then
         ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+paux(ig,igp)
       else
         ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+conjg(paux(igp,ig))
       endif
     end do !igp
   end do !ig
   ABI_DEALLOCATE(paux)
 endif
!
! Calculate sum_v (q+G).<nabla wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|nabla wfn2>.(q+Gp) + 
! ff fnltw(G)^* fnltw(Gp) + fnltw(G)^* fdrhotwg(Gp).(q+Gp) - kns(G) fnltw(Gp)
! and subtract from ptwsq(:,:,3)
!
 if (niter>1) then
!
! Calculate sum_G fnl_kp(G) fnl^*_kp(G) and store in ff
! FIXME: this won't work if k-centered G-spheres are used
!
   call compute_ff(Cryst,Kmesh,Psps,ikmq_bz,nlx,istwf_kmq,npw_kmq,kg_kmq,Gsph_max%rottbm1,fnlmax,fnlloc,ff)

   vzn(:,:,:)=(0.0,0.0)
   do igp=1,Sigp%npwc
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         vzn(igp,:,:) = vzn(igp,:,:) + half*ff(ilm,iat,:,:)*fnltwg(igp,ilm,iat)
       enddo
     enddo
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         drhaux(:)=cmplx(real(fdrhotwg(igp,ilm,iat,:)),aimag(fdrhotwg(igp,ilm,iat,:)))
         vzn(igp,ilm,iat)=vzn(igp,ilm,iat)+vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")-kns(igp,ilm,iat)
       enddo
     enddo
     do ig=1,igp
       do ibv = 1, nbmax
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3)-conjg(cauxg(ig,ibv))*cauxg(igp,ibv)
       enddo
       ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+kplqg(igp)*ptwsq(ig,igp,2)+kplqg(ig)*conjg(ptwsq(igp,ig,2))&
&                                     -kplqg(ig)*ptwsq(ig,igp,1)*kplqg(igp)
       do iat = 1, Cryst%natom
         do ilm = 1, nlx*nlx
           if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
           if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
           ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+conjg(fnltwg(ig,ilm,iat))*vzn(igp,ilm,iat)+ &
&                                                fnltwg(igp,ilm,iat)*conjg(vzn(ig,ilm,iat))
         enddo
       enddo
     end do !ig
   end do !igp

   do igp = 1, Sigp%npwc
     do ig = igp+1, Sigp%npwc
       ptwsq(ig,igp,3)=conjg(ptwsq(igp,ig,3))
     enddo
   enddo

 endif

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
   ABI_DEALLOCATE(fnltwg)
   ABI_DEALLOCATE(fnltwg2)
   ABI_DEALLOCATE(fnltwg3)
   ABI_DEALLOCATE(fnlwfg)
 endif
 if (niter>1) then
   ABI_DEALLOCATE(ddwfwfg)
   ABI_DEALLOCATE(lnkp)
   ABI_DEALLOCATE(kns)
   ABI_DEALLOCATE(fdrhotwg)
   ABI_DEALLOCATE(fkdwfg)
   ABI_DEALLOCATE(vzn)
   ABI_DEALLOCATE(ff)
 endif

end subroutine fft4eet_sig_kb
!!***

!!****f* m_eet/calc_corr_sig
!!
!! NAME
!! calc_corr_sig
!!
!! FUNCTION
!! Calculate sum-over-states correction to the matrix element of the GW self-energy obtained with the EET
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_corr_sig(Sigp,Sr,nomega,nspinor,npwc1,npwc2,botsq,otq,rhotwg,is,ibv,kb, &
&                        ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,i_sz,vc_sqrt_qbz,sigmac)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_corr_sig'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_parameters),intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr

 integer,intent(in) :: nomega,nspinor,npwc1,npwc2
 integer,intent(in) :: is,ibv,kb
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz
 real(dp),intent(in) :: i_sz
 complex(gwpc),intent(in) :: rhotwg(Sigp%npwc)
 complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)
 complex(dpc),intent(inout) :: sigmac(nomega)

 complex(gwpc),allocatable :: rtaux(:)

 integer :: ig,igp,ios
 real(dp) :: otw,twofm1_zcut
 real(dp) :: den,omegame0i
 complex(gwpc) :: num

!************************************************************************

 ABI_ALLOCATE(rtaux,(Sigp%npwc*nspinor**2))

 twofm1_zcut=-Sigp%zcut

 do ig = 1,Sigp%npwc
   rtaux(ig)=rhotwg(ig)*vc_sqrt_qbz(ig)
 enddo
 if (ik_bz==ikmq_bz) then
   rtaux(1)=czero_gw
   if (ibv==kb) then
     rtaux(1)=cmplx(sqrt(i_sz),0.0_gwp)
   endif
 endif
 do ios=1,nomega
   omegame0i = real(Sr%omega4sd(kb,ikmq_ibz,ios,is)) - Sr%e0(ibv,ik_ibz,is)
   do ig=1,Sigp%npwc
     do igp=1,Sigp%npwc
       otw = DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
       num = botsq(ig,igp)*conjg(rtaux(ig))*rtaux(igp)
       den = omegame0i-otw
       if (real(den*den)>Sigp%zcut**2) then
         sigmac(ios) = sigmac(ios) + 0.5*num/(den*otw)
       else
         sigmac(ios) = sigmac(ios) + 0.5*num*cmplx(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)
       end if
     end do !igp
   end do !ig
 end do !ios

 ABI_DEALLOCATE(rtaux)

end subroutine calc_corr_sig
!!***

!!****f* m_eet/drho_tw_g
!! NAME                  
!! drho_tw_g
!!
!! FUNCTION
!!      Calculate matrix element <wfn1|exp[-i(q+G).r]|nabla wfn1>
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine drho_tw_g(nspinor,npwvec,nr,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
&                    wfn1,i1,ktabr1,ktabp1,wfn2,dim_rtwg,rhotwg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'drho_tw_g'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: i1,npwvec,nr,nspinor,dim_rtwg,map2sphere,use_padfft
 complex(dpc),intent(in) :: ktabp1
!arrays
 !integer,intent(in) :: gbound(2*mgfft+8,2)
 integer,intent(in) :: gbound(:,:)
 integer,intent(in) :: igfftg0(npwvec*map2sphere),ngfft(18)
 integer,intent(in) :: ktabr1(nr)
 complex(gwpc),intent(in) :: wfn1(nr*nspinor),wfn2(nr*nspinor)
 complex(gwpc),intent(out) :: rhotwg(npwvec*dim_rtwg)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: ig,igfft,nx,ny,nz,ldx,ldy,ldz,mgfft
 type(fftbox_plan3_t) :: plan
!arrays
 complex(dpc),allocatable :: usk(:),uu(:)

! *************************************************************************

 SELECT CASE (nspinor)

 CASE (1) ! Collinear case.
  !
  ! Form rho-twiddle(r)=u_1^*(r,b1,kbz1) u_2(r,b2,kbz2), to account for
  ! symmetries:
  ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz)
  !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
  !
  ABI_ALLOCATE(uu,(nr))
  ABI_ALLOCATE(usk,(nr))

  uu  = wfn1(ktabr1)*ktabp1; if (i1==1) uu  = CONJG(uu)
  usk = wfn2
  uu  = uu * usk

  SELECT CASE (map2sphere)

  CASE (0) ! Need results on the full FFT box thus cannot use zero-padded FFT.

    call fftbox_plan3(plan,ngfft(1:3),ngfft(7),-1) 
    call fftbox_execute(plan,uu)
    rhotwg=uu

  CASE (1) ! Need results on the G-sphere. Call zero-padded FFT routines if required.

    if (use_padfft==1) then
      nx =ngfft(1); ny =ngfft(2); nz =ngfft(3); mgfft = MAXVAL(ngfft(1:3))
      ldx=nx      ; ldy=ny      ; ldz=nz
      call fftpad(uu,ngfft,nx,ny,nz,ldx,ldy,ldz,ndat1,mgfft,-1,gbound)
    else
      call fftbox_plan3(plan,ngfft(1:3),ngfft(7),-1) 
      call fftbox_execute(plan,uu)
    end if

    do ig=1,npwvec       ! Have to map FFT to G-sphere.
      igfft=igfftg0(ig)
      if (igfft/=0) then ! G-G0 belong to the FFT mesh.
        rhotwg(ig)=uu(igfft)
      else               ! Set this component to zero.
        rhotwg(ig)=czero_gw
      end if
    end do

  CASE DEFAULT
    MSG_BUG("Wrong map2sphere")
  END SELECT

  ABI_DEALLOCATE(uu)
  ABI_DEALLOCATE(usk)

  RETURN

 CASE DEFAULT
   MSG_BUG('Wrong nspinor')
 END SELECT

end subroutine drho_tw_g
!!***
 
!----------------------------------------------------------------------

!!****f* m_eet/calc_dwfwfg
!! NAME                  
!! calc_dwfwfg
!!
!! FUNCTION
!!      Calculate matrix element <nabla wfn2|exp[-i(q+G).r]|nabla wfn2>
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_dwfwfg(ktabr_k,ktabi_k,nfftot,ngfft_gw,ph_mkt,wfr_jb,wfr_kb,wfg2_jk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_dwfwfg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ktabi_k,nfftot
 complex(dpc),intent(in) :: ph_mkt
!arrays
 integer,intent(in) :: ktabr_k(nfftot),ngfft_gw(18)
 complex(gwpc),intent(in) :: wfr_jb(nfftot),wfr_kb(nfftot)
 complex(gwpc),intent(out) :: wfg2_jk(nfftot)

!Local variables-------------------------------
!arrays
 integer,parameter :: ndat1=1
 complex(dpc),allocatable :: wfr2_dpcplx(:)
#ifndef HAVE_GW_DPC
 complex(dpc),allocatable :: wfg2_dpcplx(:)
#endif
 type(fftbox_plan3_t) :: plan

! *************************************************************************

 ABI_ALLOCATE(wfr2_dpcplx,(nfftot))

 SELECT CASE (ktabi_k)

 CASE (1)
   wfr2_dpcplx = CONJG(ph_mkt*wfr_jb(ktabr_k)) * wfr_kb

 CASE (2) ! Conjugate the product if time-reversal is used to reconstruct this k-point
   wfr2_dpcplx = ph_mkt*wfr_jb(ktabr_k) * wfr_kb

 CASE DEFAULT
   MSG_ERROR("Wrong ktabi_k")
 END SELECT

 ! Transform to Fourier space (result in wfg2_jk)
 call fftbox_plan3(plan,ngfft_gw(1:3),ngfft_gw(7),-1) 
#ifdef HAVE_GW_DPC
 call fftbox_execute(plan,wfr2_dpcplx,wfg2_jk)
#else
 ABI_ALLOCATE(wfg2_dpcplx,(nfftot))
 call fftbox_execute(plan,wfr2_dpcplx,wfg2_dpcplx)
 wfg2_jk=wfg2_dpcplx
 ABI_DEALLOCATE(wfg2_dpcplx)
#endif

 ABI_DEALLOCATE(wfr2_dpcplx)

end subroutine calc_dwfwfg
!!***

!!****f* m_eet/calc_ddwfwfg
!! NAME                  
!! calc_ddwfwfg
!!
!! FUNCTION
!! Calculate matrix element <nabla wfn2|exp[-iG.r]|nabla wfn2>
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_ddwfwfg(ktabi_k,nfftot,ngfft_gw,wfr_jb,wfr_kb,wfg2_jk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_ddwfwfg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ktabi_k,nfftot
!arrays
 integer,intent(in) :: ngfft_gw(18)
 complex(gwpc),intent(in) :: wfr_jb(nfftot),wfr_kb(nfftot)
 complex(gwpc),intent(out) :: wfg2_jk(nfftot)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 type(fftbox_plan3_t) :: plan
!arrays
 complex(dpc),allocatable :: wfr2_dpcplx(:)
#ifndef HAVE_GW_DPC
 complex(dpc),allocatable :: wfg2_dpcplx(:)
#endif

! *************************************************************************

 ! There is no need to take into account phases arising from non-symmorphic
 ! operations since the wavefunctions are evaluated at the same k-point.
 ABI_ALLOCATE(wfr2_dpcplx,(nfftot))

 SELECT CASE (ktabi_k)

 CASE (1)
   wfr2_dpcplx = CONJG(wfr_jb) * wfr_kb

 CASE (2) ! Conjugate the product if time-reversal is used to reconstruct this k-point
   wfr2_dpcplx = CONJG(wfr_jb) * wfr_kb

 CASE DEFAULT
   MSG_ERROR("Wrong ktabi_k")
 END SELECT

 ! Transform to Fourier space (result in wfg2_jk)
 call fftbox_plan3(plan,ngfft_gw(1:3),ngfft_gw(7),-1) 
#ifdef HAVE_GW_DPC
 call fftbox_execute(plan,wfr2_dpcplx,wfg2_jk)
#else
 ABI_ALLOCATE(wfg2_dpcplx,(nfftot))
 call fftbox_execute(plan,wfr2_dpcplx,wfg2_dpcplx)
 wfg2_jk=wfg2_dpcplx
 ABI_DEALLOCATE(wfg2_dpcplx)
#endif

 ABI_DEALLOCATE(wfr2_dpcplx)

end subroutine calc_ddwfwfg
!!***

!!****f* m_eet/calc_delta_ppm
!! NAME
!! calc_delta_ppm
!!
!! FUNCTION
!!      Calculation of the ptwsq(:,:,2)/ptwsq(:,:,1) and ptwsq(:,:,3)/ptwsq(:,:,2) 
!!      necessary for the calculation of the effective energy.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_delta_ppm(Sigp,ptwsq,niter)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_delta_ppm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter
 type(Sigma_parameters),intent(in) :: Sigp
!arrays
 complex(gwpc),intent(inout) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)

!Local variables-------------------------------
!scalars
 integer :: ig,igp
!arrays
 real(dp) :: delta_huge,denchk
 complex(gwpc) :: num,den

!*************************************************************************

 delta_huge = 1.0d8

 if(niter==1) then
!
! Calculation of the ptwsq(:,:,2)/ptwsq(:,:,1)
!
   do igp=1,Sigp%npwc
     do ig=igp,Sigp%npwc
       num = ptwsq(ig,igp,2)+conjg(ptwsq(igp,ig,2))
       den = ptwsq(ig,igp,1)+conjg(ptwsq(igp,ig,1))
       denchk=abs(real(den))+abs(aimag(den))
       if (denchk>0.0) then
         ptwsq(ig,igp,2)=num/den
       else
         ptwsq(ig,igp,2)=cmplx(delta_huge,0.0)
       endif
     end do !ig
   end do !igp

   do igp=1,Sigp%npwc
     do ig=1,igp-1
       ptwsq(ig,igp,2)=conjg(ptwsq(igp,ig,2))
     enddo
   enddo

 elseif(niter==2) then
!
! Calculation of the ptwsq(:,:,3)/ptwsq(:,:,2) and ptwsq(:,:,2)/ptwsq(:,:,1)
!
   do igp=1,Sigp%npwc
     do ig=igp,Sigp%npwc
       num = ptwsq(ig,igp,3)+conjg(ptwsq(igp,ig,3))
       den = ptwsq(ig,igp,2)+conjg(ptwsq(igp,ig,2))
       denchk=abs(real(den))+abs(aimag(den))
       if (denchk>0.0) then
         ptwsq(ig,igp,3)=num/den
       else
         ptwsq(ig,igp,3)=cmplx(delta_huge,0.0)
       endif
       num = ptwsq(ig,igp,2)+conjg(ptwsq(igp,ig,2))
       den = ptwsq(ig,igp,1)+conjg(ptwsq(igp,ig,1))
       denchk=abs(real(den))+abs(aimag(den))
       if (denchk>0.0) then
         ptwsq(ig,igp,2)=num/den
       else
         ptwsq(ig,igp,2)=cmplx(delta_huge,0.0)
       endif
     end do !ig
   end do !igp

   do igp=1,Sigp%npwc
     do ig=1,igp-1
       ptwsq(ig,igp,2)=conjg(ptwsq(igp,ig,2))
       ptwsq(ig,igp,3)=conjg(ptwsq(igp,ig,3))
     enddo
   enddo

 endif

end subroutine calc_delta_ppm
!!***

!!****f* m_eet/calc_sig_ppm_delta
!! NAME
!! calc_sig_ppm_delta
!!
!! FUNCTION
!!      Calculation of the part of the matrix elements of the self-energy coming from the sum
!!      over the valence states in OPTIMAL GW
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_sig_ppm_delta(npwc,nomega,rhotwgp,botsq,otq,omegame0i,zcut,theta_mu_minus_e0i, &
&                             ket,npwx,npwc1,npwc2,omega4sd,e0,delta)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sig_ppm_delta'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwc1,npwc2,npwx
 real(dp),intent(in) :: theta_mu_minus_e0i,zcut
!arrays
 real(dp),intent(in) :: omegame0i(nomega),e0
 complex(dpc),intent(in) :: omega4sd(nomega)
 complex(gwpc),intent(in) :: otq(npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(npwc,npwc1)
 complex(gwpc),intent(in) :: rhotwgp(npwx)
 complex(gwpc),intent(in) :: delta(npwc,npwc,nomega)
 complex(gwpc),intent(inout) :: ket(npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ios
 real(dp) :: den,omegame0i_io,otw,twofm1_zcut
 complex(gwpc) :: num,den2,rhotwgdp_igp
 logical :: fully_occupied,totally_empty

!*************************************************************************

  fully_occupied=(abs(theta_mu_minus_e0i-1.)<0.001)
  totally_empty=(abs(theta_mu_minus_e0i)<0.001)

  if (.not.(totally_empty)) then
   twofm1_zcut=zcut

   do ios=1,nomega
    omegame0i_io=omegame0i(ios)
    do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
      otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
      num = botsq(ig,igp)*rhotwgdp_igp
      den = omegame0i_io+otw

      if (den**2>zcut**2)then
       ket(ig,ios) = ket(ig,ios) + num/(den*otw) * theta_mu_minus_e0i
      else
       ket(ig,ios) = ket(ig,ios) + num*cmplx(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)&
&                            *theta_mu_minus_e0i
      end if
     end do !ig
    end do !igp
   end do !ios

  end if !not totally empty


  if (.not.(totally_empty)) then
   twofm1_zcut=-zcut

   do ios=1,nomega
    omegame0i_io=DBLE(omega4sd(ios)) - e0
    do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
      otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
      num = botsq(ig,igp)*rhotwgdp_igp
      den2 = omegame0i_io-otw-delta(ig,igp,ios)

      if (real(conjg(den2)*den2)>zcut**2) then
       ket(ig,ios) = ket(ig,ios) - num/(den2*otw)*theta_mu_minus_e0i
      else
       ket(ig,ios) = ket(ig,ios) - num*(den2+cmplx(0.0,twofm1_zcut))/((den2**2+twofm1_zcut**2)*otw) &
&                           *theta_mu_minus_e0i
      end if
     end do !ig
    end do !igp
   end do !ios

  end if

  if (.not.(fully_occupied)) then
   twofm1_zcut=-zcut

   do ios=1,nomega
    omegame0i_io=DBLE(omega4sd(ios)) - e0
    do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
      otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
      num = botsq(ig,igp)*rhotwgdp_igp
      den2 = omegame0i_io-otw-delta(ig,igp,ios)

      if (real(conjg(den2)*den2)>zcut**2) then
       ket(ig,ios) = ket(ig,ios) - num/(den2*otw)*(1.-theta_mu_minus_e0i)
      else
       ket(ig,ios) = ket(ig,ios) - num*(den2+cmplx(0.0,twofm1_zcut))/((den2**2+twofm1_zcut**2)*otw) &
&                           *(1.-theta_mu_minus_e0i)
      end if
     end do !ig
    end do !igp
   end do !ios

   do ios=1,nomega
    omegame0i_io=omegame0i(ios)
    do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
      otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
      num = botsq(ig,igp)*rhotwgdp_igp

      den = omegame0i_io-otw
      if (den**2>zcut**2) then
       ket(ig,ios) = ket(ig,ios) + num/(den*otw)*(1.-theta_mu_minus_e0i)
      else
       ket(ig,ios) = ket(ig,ios) + num*cmplx(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw) &
&                           *(1.-theta_mu_minus_e0i)
      end if
     end do !ig
    end do !igp
   end do !ios

  end if

  ket(:,:)=ket(:,:)*0.5

end subroutine calc_sig_ppm_delta
!!***

!----------------------------------------------------------------------

!!****f* m_eet/calc_sig_ppm_delta_clos
!! NAME
!! calc_sig_ppm_delta_clos
!!
!! FUNCTION
!! Calculation of the part of the GW matrix elements of the self-energy coming from
!! the sum over all the states using the EET
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE
subroutine calc_sig_ppm_delta_clos(npwc,nomega,ikbz,jkbz,qbzpg,botsq,otq,omegame0k,omegame0lumo,zcut,ptwsq, &
&                                   ket,npwc1,npwc2,gw_eet_scale,niter)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sig_ppm_delta_clos'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwc1,npwc2,niter
 integer,intent(in) :: ikbz,jkbz
 real(dp),intent(in) :: zcut,gw_eet_scale
!arrays
 real(dp),intent(in) :: omegame0k(nomega),omegame0lumo(nomega)
 real(dp),intent(in) :: qbzpg(npwc)
 complex(gwpc),intent(in) :: otq(npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(npwc,npwc1)
 complex(gwpc),intent(in) :: ptwsq(npwc,npwc,niter+1)
 complex(dpc),intent(inout) :: ket(nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ios
 real(dp) :: omegame0i_io,otw,omegakg,twofm1_zcut,bq
 complex(gwpc) :: den,num,delta,deltaux
 logical :: ggpnonzero
!arrays
 real(dp), allocatable :: qpgsq(:)

!*************************************************************************

   twofm1_zcut=-zcut

   ABI_ALLOCATE(qpgsq,(npwc))

   do ig = 1, npwc
     qpgsq(ig)=half*qbzpg(ig)*qbzpg(ig)
   enddo

   do igp=1,npwc
     do ig=1,npwc
       otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
       ggpnonzero=(ig/=1.and.igp/=1)
       if (ikbz/=jkbz.or.ggpnonzero) then
         bq=half*(qpgsq(ig)+qpgsq(igp))
         if (niter==0) then
           delta = bq
         elseif(niter==1) then
           delta = bq+ptwsq(ig,igp,2)
         endif
         do ios=1,nomega
           if(niter==2) then
             omegakg=omegame0k(ios)-otw
             delta = bq + ptwsq(ig,igp,2)*(omegakg-bq-ptwsq(ig,igp,2))/(omegakg-bq-ptwsq(ig,igp,3))
           endif
           if (niter<2) deltaux=delta
           call check_delta_sigma(qpgsq(ig),qpgsq(igp),delta,omegame0k(ios),omegame0lumo(ios),ig,igp,gw_eet_scale,niter)
           omegame0i_io=omegame0k(ios)
           num = botsq(ig,igp)*ptwsq(ig,igp,1)
           den = omegame0i_io-otw-delta
           if (real(conjg(den)*den)>zcut**2) then
             ket(ios) = ket(ios) + 0.5*num/(den*otw)
           else
             ket(ios) = ket(ios) + 0.5*num*(den+cmplx(0.0,twofm1_zcut))/((den**2+twofm1_zcut**2)*otw)
           end if
           if (niter<2) delta=deltaux
         end do !ios
       else
         do ios=1,nomega
           omegame0i_io=omegame0k(ios)
           num = botsq(ig,igp)*ptwsq(ig,igp,1)
           den = omegame0i_io-otw
           if (real(conjg(den)*den)>zcut**2) then
             ket(ios) = ket(ios) + 0.5*num/(den*otw)
           else
             ket(ios) = ket(ios) + 0.5*num*(den+cmplx(0.0,twofm1_zcut))/((den**2+twofm1_zcut**2)*otw)
           end if
         end do !ios
       endif
     end do !ig
   end do !igp

 ABI_DEALLOCATE(qpgsq)

end subroutine calc_sig_ppm_delta_clos
!!***

!!****f* m_eet/check_delta_sigma
!!
!! NAME
!! check_delta_sigma
!!
!! FUNCTION
!! Check if the effective energy (delta) is in the correct range
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine check_delta_sigma(qpgsq,qpgpsq,delta,omegame0k,omegame0lumo,ig,igp,gw_eet_scale,niter)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_delta_sigma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: ig,igp,niter
 real(dp),intent(in) :: qpgsq,qpgpsq
 real(dp),intent(in) :: omegame0k,omegame0lumo
 real(dp),intent(in) :: gw_eet_scale
 complex(gwpc),intent(inout) :: delta

!Local variables-------------------------------
!scalars
 real(dp) :: test

!*************************************************************************

 if (gw_eet_scale>0.01) then
   delta=gw_eet_scale*delta
 endif

 if (ig==igp) then
   delta=real(delta)
   test=omegame0k-real(delta)
   if (test>omegame0lumo.and.niter>0) then
     delta=half*(qpgsq+qpgpsq)
     test=omegame0k-real(delta)
   endif
   if (test>omegame0lumo) then
     delta=omegame0k-omegame0lumo
   endif
 else
   test=omegame0k-real(delta)
   if (test>omegame0lumo) then
     delta=omegame0k-omegame0lumo
   endif
 endif

end subroutine check_delta_sigma
!!***

!!****f* calc_sig_ppm_eet/compute_ff
!! NAME
!!  compute_ff
!!
!! FUNCTION
!!      Calculate sum_G fnl_kp(G) fnl^*_kp(G)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine compute_ff(Cryst,Kmesh,Psps,ik_bz,nlx,istwf_k,npw_k,kg_k,grottbm1,fnlmax,fnlloc,ff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compute_ff'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,istwf_k,nlx,ik_bz
 type(crystal_t),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 type(kmesh_t),intent(in) :: Kmesh
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 integer,intent(in) :: grottbm1(npw_k,2,Cryst%nsym)
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)
 complex(gwpc),intent(out) :: ff(nlx*nlx,Cryst%natom,nlx*nlx,Cryst%natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: inclvkb2=2
 integer :: iat,iat2,ilm,ilm2,ig,igbz,itim_k,isym_k,ik_ibz,istwf_kbz
 type(kb_potential) :: KBff_k_ibz
!arrays
 real(dp) :: kbz(3)

!************************************************************************

 ! @Arjan:
 ! Could you check if it's possible to use the version below that uses the Gvectors
 ! centered on the point in the full BZ and thus avoids the symmetrization (groottbm1)?
 ! When we enable the reading of the WFK file, one should call the routine with the 
 ! kg_k array centered on the point in the full brillouin zone.

 call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k)

 !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz

 if (.TRUE.) then
   ! old version with symmetrization
                                                                                              
   call init_kb_potential(KBff_k_ibz,Cryst,Psps,2,istwf_k,npw_k,Kmesh%ibz(:,ik_ibz),kg_k)
   ABI_DEALLOCATE(KBff_k_ibz%fnld)
                                                                                              
   ff(:,:,:,:)=czero_gw
   do iat = 1, Cryst%natom
     do ilm = 1, nlx*nlx
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       do iat2 = 1, Cryst%natom
         do ilm2 = 1, nlx*nlx
           if (ilm2>fnlmax(Cryst%typat(iat2))) CYCLE
           if (ilm2>=fnlloc(Cryst%typat(iat2),1).and.ilm2<=fnlloc(Cryst%typat(iat2),2)) CYCLE
           do ig = 1, npw_k
             ! FIXME: this won't work if k-centered G-spheres are used
!
! Here grottbm1 should be replaced by iskg_tabs_t
!
             igbz = grottbm1(ig,itim_k,isym_k)
             if (itim_k==1) then
               ff(ilm,iat,ilm2,iat2)=ff(ilm,iat,ilm2,iat2)+ &
&                conjg(KBff_k_ibz%fnl(igbz,ilm,iat))*KBff_k_ibz%fnl(igbz,ilm2,iat2)
             else
               ff(ilm,iat,ilm2,iat2)=ff(ilm,iat,ilm2,iat2)+ &
&                KBff_k_ibz%fnl(igbz,ilm,iat)*conjg(KBff_k_ibz%fnl(igbz,ilm2,iat2))
             endif
           enddo
         enddo
       enddo
     enddo
   enddo

 else
   ! This one uses the point in the BZ thus avoiding the symmetrization
   istwf_kbz = set_istwfk(kbz)
   ABI_CHECK(istwf_kbz==1, "istwf_kbz!=1 not coded")

   call init_kb_potential(KBff_k_ibz,Cryst,Psps,inclvkb2,istwf_kbz,npw_k,kbz,kg_k)
   ABI_DEALLOCATE(KBff_k_ibz%fnld)
                                                                                              
   ff(:,:,:,:)=czero_gw
   do iat = 1, Cryst%natom
     do ilm = 1, nlx*nlx
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       do iat2 = 1, Cryst%natom
         do ilm2 = 1, nlx*nlx
           if (ilm2>fnlmax(Cryst%typat(iat2))) CYCLE
           if (ilm2>=fnlloc(Cryst%typat(iat2),1).and.ilm2<=fnlloc(Cryst%typat(iat2),2)) CYCLE
           do ig = 1, npw_k
             ! FIXME: this won't work if k-centered G-spheres are used
             !igbz = grottbm1(ig,itim_k,isym_k)
             igbz = ig
             !if (itim_k==1) then
               ff(ilm,iat,ilm2,iat2)=ff(ilm,iat,ilm2,iat2)+ &
&                conjg(KBff_k_ibz%fnl(igbz,ilm,iat))*KBff_k_ibz%fnl(igbz,ilm2,iat2)
             !else
             !  ff(ilm,iat,ilm2,iat2)=ff(ilm,iat,ilm2,iat2)+ &
             !   KBff_k_ibz%fnl(igbz,ilm,iat)*conjg(KBff_k_ibz%fnl(igbz,ilm2,iat2))
             !endif
           enddo
         enddo
       enddo
     enddo
   enddo

 end if
                                                                                              
 call destroy_kb_potential(KBff_k_ibz)

end subroutine compute_ff
!!***

!!****f* m_eet/gw_eet_chi0
!! NAME
!! gw_eet_chi0
!!
!! FUNCTION
!! Wrapper routine for the calculation of the polarizability using the EET
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine gw_eet_chi0(Ep,Dtset,Cryst,Wfs,Kmesh,Gsph_epsG0,Gsph_wfn,Psps,Ltg_q,nbvw,qpoint, &
&                      nfftot_gw,ngfft_gw,use_padfft,igfftepsG0,gw_gbound,gw_mgfft,is, &
&                      ik_bz,ik_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                      ikmq_ibz,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg, &
&                      qp_energy,chi0,spin_fact,qp_occ,nspinor,tim_fourdp,bbp_ks_distrb,nbmax)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_eet_chi0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,ik_ibz,ikmq_ibz
 integer,intent(in) :: nbvw,nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: is,isym_k,itim_k,itim_kmq
 integer,intent(in) :: dim_rtwg
 integer,intent(out) :: nbmax
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(gsphere_t),intent(in) :: Gsph_epsG0,Gsph_wfn
 type(Dataset_type),intent(in) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfd_t),target,intent(inout) :: Wfs
 type(littlegroup_t),intent(in) :: Ltg_q
!arrays
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: bbp_ks_distrb(nbvw,Kmesh%nbz,Wfs%nsppol)
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Ep%npwepG0)
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4),spin_fact
 real(dp),intent(in) :: qpoint(3)
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: ib,ibv
 integer :: niter,nptwg,iter
 integer :: ig
!arrays
 integer :: nbhomo(2)
 integer :: fnlloc(Cryst%ntypat,2),fnlmax(Cryst%ntypat)
 real(dp),allocatable :: qpgsq(:)
 real(dp),allocatable :: qplg(:,:)
 real(dp),allocatable :: kplqg(:)
 complex(gwpc),allocatable :: fnlkr(:,:,:)
 complex(gwpc),allocatable :: fnlkpr(:,:,:)
 complex(gwpc),allocatable :: wfr1(:,:)
 complex(gwpc),allocatable :: wfr2(:)
 complex(gwpc),allocatable :: mtwk(:,:)
 complex(gwpc),allocatable :: mtwkp(:,:)
 complex(gwpc),allocatable :: frhorho(:)
 complex(gwpc),allocatable :: frhoj(:,:)
 complex(gwpc),allocatable :: fjj(:)
!************************************************************************

 nbhomo(:)=1
 do ib = 2, Ep%nbnds
   if (spin_fact*qp_occ(ib,ikmq_ibz,is)<GW_TOL_DOCC) exit
   nbhomo(1)=nbhomo(1)+1
 enddo
 do ib = 2, Ep%nbnds
   if (spin_fact*qp_occ(ib,ikmq_ibz,is)<(one-GW_TOL_DOCC)) exit
   nbhomo(2)=nbhomo(2)+1
 enddo
 nbmax=max(nbhomo(1),Dtset%gw_eet_nband)

 niter = Dtset%gw_eet
 nptwg=1
 do iter = 1, niter
   nptwg=nptwg+mod(iter,2)
 enddo

 ABI_ALLOCATE(qpgsq,(Ep%npwe))
 ABI_ALLOCATE(qplg,(3,Ep%npwe))
 ABI_ALLOCATE(kplqg,(Ep%npwe))
!
!  Calculate some auxiliary quantities
!
 do ig=1,Ep%npwe
   qplg(:,ig) = qpoint(:)+ Gsph_epsG0%gvec(:,ig)
   kplqg(ig)=-vdotw(Kmesh%bz(:,ik_bz),qplg(:,ig),Cryst%gmet,"G")
   qpgsq(ig) = normv(qplg(:,ig),Cryst%gmet,"G")
 enddo
 qpgsq(:)=half*qpgsq(:)**2

 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then
   ABI_ALLOCATE(fnlkpr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlkr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(mtwk,(Wfs%nfftot*nspinor,nbhomo(1)))
   ABI_ALLOCATE(mtwkp,(Wfs%nfftot*nspinor,nbmax))
!
! Calculate FFTs of f_{nlk}(G) and f_{nlkp}(G) and the auxiliary functions 
! Mtw_k(r) = \sum_{nl} f_{nlk}(r) * \sum_{G} u_{vk}(G) f_{nlk}(G) and 
! Mtw_kp(r) = \sum_{nl} f_{nlkp}(r) * \sum_{G} u_{vkp}(G) f_{nlkp}(G)
!
   call gw_eet_chi0_vkb(Ep,Cryst,Wfs,Kmesh,Psps,is,ik_ibz,ikmq_ibz,nspinor,tim_fourdp, &
&                       nbhomo,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)
 endif

 ABI_ALLOCATE(wfr1,(Wfs%nfftot*nspinor,nbmax))

 do ibv = 1, nbmax
   call wfd_get_ur(Wfs,ibv,ikmq_ibz,is,wfr1(:,ibv))
 enddo

 do ibv = 1, nbhomo(1)

   if ((bbp_ks_distrb(ibv,ik_bz,is) /= Wfs%my_rank)) CYCLE

   ABI_ALLOCATE(frhorho,(Ep%npwe*(Ep%npwe+1)/2))
   ABI_CHECK_ALLOC('out of memory in frhorho')
   if (niter>0) then
     ABI_ALLOCATE(frhoj,(Ep%npwe,Ep%npwe))
     ABI_CHECK_ALLOC('out of memory in frhoj')
     ABI_ALLOCATE(fjj,(Ep%npwe*(Ep%npwe+1)/2*(niter-1)))
     ABI_CHECK_ALLOC('out of memory fij')
   endif

   ABI_ALLOCATE(wfr2,(Wfs%nfftot*nspinor))
!
! Get FFTs of the wave functions
!
   call wfd_get_ur(Wfs,ibv,ik_ibz,is,wfr2)
!
! Calculate f^{pp}, f^{pj}, and f{jj}. See Eqs. (25),(26),(27) of PRB 85, 085126, for their expressions.
! They are stored in frhorho, frhoj, and fjj.
!
   if (niter==0) then
! Case: zero order
     call fft4eet_0(ik_bz,Ep,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax,is,nfftot_gw,ngfft_gw, &
&                   use_padfft,igfftepsG0,gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,itim_k,tabr_k, &
&                   ph_mkt,spinrot_k,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg, &
&                   nspinor,tim_fourdp,wfr1,wfr2,ibv,frhorho,spin_fact,qp_occ,qp_energy,chi0,1)
   else
! Case: with vkb
     if (Dtset%gw_eet_inclvkb==1) then
       call fft4eet_kb(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,Psps,nbhomo,nbmax, &
&                      is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0, &
&                      gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                      itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,Gsph_wfn%rottbm1, &
&                      nspinor,tim_fourdp,fnlloc,fnlmax,fnlkpr,mtwk,mtwkp,wfr1,wfr2,ibv,qplg,kplqg, &
&                      niter,frhorho,frhoj,fjj,spin_fact,qp_occ,qp_energy,chi0,1)
     else
! Case: without vkb
       call fft4eet(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax, &
&                   is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0, &
&                   gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                   itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,Gsph_wfn%rottbm1, &
&                   nspinor,tim_fourdp,wfr1,wfr2,ibv,qplg,kplqg,niter,frhorho,frhoj,fjj, &
&                   spin_fact,qp_occ,qp_energy,chi0,1)
     endif
   endif

   ABI_DEALLOCATE(wfr2)
   if (niter==0) then
     frhorho=spin_fact*qp_occ(ibv,ik_ibz,is)*frhorho
!
!    Calculate the polarizability using the closure relation
!
     call calc_chi0_delta0(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,qp_energy(ibv,ik_ibz,is), &
&                          qp_energy(nbmax+1,ikmq_ibz,is),qpgsq,frhorho,chi0)
   else
!
!    Calculation of the frhoj/frhorho and fjj/frhoj and store in frhoj and fjj
!
     call calc_delta_chi0(Ep,frhorho,frhoj,fjj,niter)
     frhorho=spin_fact*qp_occ(ibv,ik_ibz,is)*frhorho
!
!    Calculate the polarizability using the closure relation
!
     call calc_chi0_delta_clos(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,niter,qp_energy(ibv,ik_ibz,is), &
&                              qp_energy(nbmax+1,ikmq_ibz,is),qpgsq,frhorho,frhoj,fjj,chi0)
     ABI_DEALLOCATE(frhoj)
     ABI_DEALLOCATE(fjj)
   endif

   ABI_DEALLOCATE(frhorho)

 enddo

 ABI_DEALLOCATE(wfr1)
 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then
   ABI_DEALLOCATE(mtwk)
   ABI_DEALLOCATE(mtwkp)
   ABI_DEALLOCATE(fnlkr)
   ABI_DEALLOCATE(fnlkpr)
 endif

 ABI_DEALLOCATE(qpgsq)
 ABI_DEALLOCATE(qplg)
 ABI_DEALLOCATE(kplqg)

 nbmax=0

end subroutine gw_eet_chi0
!!***

!!****f* m_eet/gw_eet_chi0_vkb
!! NAME
!! gw_eet_chi0_vkb
!!
!! FUNCTION
!!      Calculate FFTs of f_{nlk}(G) and f_{nlkp}(G) and the auxiliary functions 
!!      Mtw_k(r) = \sum_{nl} f_{nlk}(r) * \sum_{G} u_{vk}(G) f_{nlk}(G) and 
!!      Mtw_kp(r) = \sum_{nl} f_{nlkp}(r) * \sum_{G} u_{vkp}(G) f_{nlkp}(G)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine gw_eet_chi0_vkb(Ep,Cryst,Wfs,Kmesh,Psps,is,ik_ibz,ikmq_ibz,nspinor,tim_fourdp, &
&                          nbhomo,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_eet_chi0_vkb'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfd_t),target,intent(inout) :: Wfs
 type(kb_potential) :: KBff_k_ibz,KBff_kmq_ibz

 integer,intent(in) :: nspinor,tim_fourdp
 integer,intent(in) :: is
 integer,intent(in) :: nbhomo(2), nbmax
 integer,intent(in) :: ik_ibz,ikmq_ibz

 integer,intent(out) :: fnlloc(Cryst%ntypat,2),fnlmax(Cryst%ntypat)

 complex(gwpc),intent(out) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(out) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbhomo(1))
 complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,nbmax)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: npw_k,npw_kmq,i,ilm,iat,istwf_k,istwf_kmq,temp_unit
 integer :: lloc,lmax,mmax,pspcod,pspdat,pspxc,ityp
 real(dp) :: r2well,zion,znucl
 character(len=fnlen) :: title
!arrays
 integer,pointer :: gbound_k(:,:),gbound_kmq(:,:),kg_k(:,:),kg_kmq(:,:)
 complex(gwpc),allocatable :: sfnl_tmp(:) 

!************************************************************************

 ABI_UNUSED(tim_fourdp)

   fnlloc(:,:)=0
   fnlmax(:)=0
   do ityp = 1, Cryst%ntypat
     temp_unit = get_unit()
     open (unit=temp_unit,file=psps%filpsp(ityp),form='formatted',status='old')
     rewind (unit=temp_unit)
     read (temp_unit,'(a)') title
     read (temp_unit,*) znucl,zion,pspdat
     read (temp_unit,*) pspcod,pspxc,lmax,lloc,mmax,r2well
     close(temp_unit)
     do i = 1, lloc
       fnlloc(ityp,1) = fnlloc(ityp,1) + 2*(lloc-1)+1
     enddo
     fnlloc(ityp,1)=fnlloc(ityp,1)+1
     do i = 1, lloc+1
       fnlloc(ityp,2) = fnlloc(ityp,2) + 2*(lloc-1)+1
     enddo
     do i = 1, lmax+1
       fnlmax(ityp) = fnlmax(ityp) + 2*(lmax-1)+1
     enddo
   enddo

   !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
   npw_k    = Wfs%npwarr(ik_ibz)
   istwf_k  = Wfs%istwfk(ik_ibz)
   kg_k     => Wfs%Kdata(ik_ibz)%kg_k
   gbound_k => Wfs%Kdata(ik_ibz)%gbound

   npw_kmq   = Wfs%npwarr(ikmq_ibz)
   istwf_kmq = Wfs%istwfk(ikmq_ibz)
   kg_kmq     => Wfs%Kdata(ikmq_ibz)%kg_k
   gbound_kmq => Wfs%Kdata(ikmq_ibz)%gbound

   call init_kb_potential(KBff_k_ibz  ,Cryst,Psps,2,istwf_k  ,npw_k,  Kmesh%ibz(:,ik_ibz),kg_k)
   call init_kb_potential(KBff_kmq_ibz,Cryst,Psps,2,istwf_kmq,npw_kmq,Kmesh%ibz(:,ikmq_ibz),kg_kmq)

   ABI_DEALLOCATE(KBff_k_ibz%fnld)
   ABI_DEALLOCATE(KBff_kmq_ibz%fnld)

   ! MG: TODO Be careful here as the symmetry properties of fnl are not the same as the 
   !     ones of the wavefunctions. One should write a separate methods for the FFT of fnl.
   if (istwf_k>1.or.istwf_kmq>1) then
     MSG_ERROR("istwfk /= 1 not coded")
   end if

   ABI_MALLOC(sfnl_tmp, (MAX(npw_k,npw_kmq)) )
!
! Calculate FFTs of fnl_k(G) and fnl_kp(G)
!
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE

       sfnl_tmp(1:npw_k) = conjg(KBff_k_ibz%fnl(:,ilm,iat))

       call fft_ug(npw_k,Wfs%nfftot,nspinor,ndat1,Wfs%mgfft,Wfs%ngfft,istwf_k,kg_k,gbound_k,sfnl_tmp,fnlkr(:,ilm,iat))

       sfnl_tmp(1:npw_kmq) = conjg(KBff_kmq_ibz%fnl(:,ilm,iat))

       call fft_ug(npw_kmq,Wfs%nfftot,nspinor,ndat1,Wfs%mgfft,Wfs%ngfft,istwf_kmq,kg_kmq,gbound_kmq,sfnl_tmp,fnlkpr(:,ilm,iat))

     enddo
   enddo

   ABI_FREE(sfnl_tmp)

   call destroy_kb_potential(KBff_k_ibz  )
   call destroy_kb_potential(KBff_kmq_ibz)
!
! Calculate auxiliary functions Mtw_k(r) = \sum_{nl} fnl_k(r) * \sum_{G} u_{vk}(G) fnl_k(G)
! and Mtw_kp(r) = \sum_{nl} fnl_kp(r) * \sum_{G} u_{vkp}(G) fnl_kp(G)
!
   call calc_eet_prep(Ep,Cryst,Wfs,Kmesh,Psps,is,nbhomo(1),nbmax,ik_ibz,ikmq_ibz,nspinor, &
&                     fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)

end subroutine gw_eet_chi0_vkb
!!***

!!****f* m_eet/calc_eet_prep
!! NAME                  
!! calc_eet_prep
!!
!! FUNCTION
!!      Calculate auxiliary functions Mtw_k(r) = \sum_{nl} f_{nlk}(r) * \sum_{G} u_{vk}(G) f_{nlk}(G)
!!      and Mtw_kp(r) = \sum_{nl} f_{nlkp}(r) * \sum_{G} u_{vkp}(G) f_{nlkp}(G)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_eet_prep(Ep,Cryst,Wfs,Kmesh,Psps,is,nbhomo,nbmax,ik_ibz,ikmq_ibz,nspinor, &
&                        fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_eet_prep'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: is,nbhomo,nbmax,nspinor,ik_ibz,ikmq_ibz
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfd_t),target,intent(in) :: Wfs
!arrays
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)
 complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(in) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbhomo)
 complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,nbmax)

!Local variables-------------------------------
!scalars
 integer :: ibv,ilm,iat,ig,istwf_k,istwf_kmq,npw_k,npw_kmq
 type(kb_potential) :: KBff_k_ibz,KBff_kmq_ibz
!arrays
 integer,pointer :: kg_k(:,:),kg_kmq(:,:)
 complex(gwpc),allocatable :: maux(:,:)

!************************************************************************

 ABI_UNUSED(Ep%npwe)

 !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
 npw_k   = Wfs%npwarr(ik_ibz)
 istwf_k = Wfs%istwfk(ik_ibz)
 kg_k    => Wfs%Kdata(ik_ibz)%kg_k

 npw_kmq   = Wfs%npwarr(ikmq_ibz)
 istwf_kmq = Wfs%istwfk(ikmq_ibz)
 kg_kmq    => Wfs%Kdata(ikmq_ibz)%kg_k

 ! MG: TODO Be careful here as the symmetry properties of fnl are not the same as the 
 !     ones of the wavefunctions. One should write a separate methods for the FFT of fnl.
 if (istwf_k>1.or.istwf_k>1) then
   MSG_ERROR("istwfk /= 1 not coded")
 end if

 call init_kb_potential(KBff_k_ibz  ,Cryst,Psps,2,istwf_k,  npw_k,  Kmesh%ibz(:,  ik_ibz),kg_k)
 call init_kb_potential(KBff_kmq_ibz,Cryst,Psps,2,istwf_kmq,npw_kmq,Kmesh%ibz(:,ikmq_ibz),kg_kmq)

 ABI_DEALLOCATE(KBff_k_ibz%fnld)
 ABI_DEALLOCATE(KBff_kmq_ibz%fnld)

 ABI_ALLOCATE(maux,(Psps%mpsang*Psps%mpsang,Cryst%natom))
!
! Calculate auxiliary function Mtw_k(r) = \sum_{nl} fnl_k(r) * \sum_{G} u_{vk}(G) fnl_k(G)
!
 mtwk(:,:)=(0.0,0.0)
 mtwkp(:,:)=(0.0,0.0)
 do ibv = 1, nbhomo
   maux(:,:)=(0.0,0)
   do ig=1,npw_k
     maux(:,:) = maux(:,:) + Wfs%Wave(ibv,ik_ibz,is)%ug(ig) * KBff_k_ibz%fnl(ig,:,:)
   enddo
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       mtwk(:,ibv)=mtwk(:,ibv)+maux(ilm,iat)*fnlkr(:,ilm,iat)
     enddo
   enddo
 enddo
!
! Calculate auxiliary function Mtw_kp(r) = \sum_{nl} fnl_kp(r) * \sum_{G} u_{vkp}(G) fnl_kp(G)
!
 do ibv = 1, nbmax
   maux(:,:)=(0.0,0)
   do ig=1,npw_kmq
     maux(:,:) = maux(:,:) + Wfs%Wave(ibv,ikmq_ibz,is)%ug(ig)* KBff_kmq_ibz%fnl(ig,:,:)
   enddo
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       mtwkp(:,ibv)=mtwkp(:,ibv)+maux(ilm,iat)*fnlkpr(:,ilm,iat)
     enddo
   enddo
 enddo

 ABI_DEALLOCATE(maux)

 call destroy_kb_potential(KBff_k_ibz)
 call destroy_kb_potential(KBff_kmq_ibz)

end subroutine calc_eet_prep
!!***

!!****f* m_eet/fft4eet_0
!! NAME                  
!! fft4eet_0
!!
!! FUNCTION
!! Calculate f^{pp}. See Eqs. (25) of PRB 85, 085126, for their expressions.
!! It is stored in frhorho
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine fft4eet_0(ik_bz,Ep,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax,is,nfftot_gw,ngfft_gw, &
&                    use_padfft,igfftepsG0,gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,itim_k,tabr_k,ph_mkt, &
&                    spinrot_k,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,nspinor,tim_fourdp, &
&                    wfr1,wfr2,ibv,frhorho,spin_fact,qp_occ,qp_energy,chi0,igstart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_0'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfd_t),target,intent(inout) :: Wfs
 type(kmesh_t),intent(in) :: Kmesh
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(littlegroup_t),intent(in) :: Ltg_q
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt
 real(dp),intent(in) :: spin_fact
 integer,intent(in) :: dim_rtwg
 integer,intent(in) :: ik_ibz,ikmq_ibz
 integer,intent(in) :: itim_k,itim_kmq
 integer,intent(in) :: nbmax,is,ibv,igstart,ik_bz
 integer,intent(in) :: nfftot_gw,nspinor,tim_fourdp,use_padfft,gw_mgfft
!arrays
 integer,intent(in) :: nbhomo(2),ngfft_gw(18)
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Ep%npwepG0)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: wfr2(Wfs%nfftot*nspinor)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
 complex(gwpc),intent(out) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: ibvp,ig,igp,istwf_k,npw_k,ig4,ig4x,ig4y,ig4z,igaux
 integer :: outofbox
 integer,save :: enough=0
 character(len=500) :: msg
!arrays
 integer :: gmgp(3)
 integer,pointer :: gbound_k(:,:),kg_k(:,:),igfft_k(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: wfwfg(:)

!************************************************************************

 ABI_UNUSED(tim_fourdp)

 ABI_ALLOCATE(rhotwg,(Ep%npwepG0*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 istwf_k  =  Wfs%istwfk(ik_ibz)
 npw_k    =  Wfs%npwarr(ik_ibz)
 gbound_k => Wfs%Kdata(ik_ibz)%gbound
 kg_k     => Wfs%Kdata(ik_ibz)%kg_k    
 igfft_k  => Wfs%Kdata(ik_ibz)%igfft0
!
! Calculate matrix element <wfn2|exp[-iG.r]|wfn2>
!
 call calc_wfwfg(tabr_k,itim_k,nfftot_gw,ngfft_gw, wfr2,wfr2,wfwfg(:))

 do ibvp = 1, nbmax
!
! Calculate matrix element <wfn1|exp[-i(q+G).r]|wfn2>
!
   call rho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                spinrot_k,dim_rtwg,rhotwg(:,ibvp))
!
! Add sum-over-states correction to the polarizability obtained with the EET
!
   if (ibvp>nbhomo(2)) then
     call calc_corr_chi0(ik_bz,Ep,Kmesh,Gsph_epsG0,Ltg_q,rhotwg(:,ibvp),spin_fact,qp_occ,qp_energy, &
&                        ibv,ibvp,ik_ibz,ikmq_ibz,is,chi0)
   endif

 enddo
!
! Calculate <wfn2|exp[i(G-Gp).r]|wfn2>. It is stored in frhorho
!
 frhorho(:)=czero_gw
 outofbox=0
 do ig=igstart,Ep%npwe
   do igp=ig,Ep%npwe
     gmgp(:)=kg_k(:,ig)-kg_k(:,igp)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)
     igaux=igp*(igp-1)/2+ig
     frhorho(igaux)=wfwfg(ig4)
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if
!
! Calculate sum_v <wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|wfn2> and subtract from frhorho
!
 do ibvp = 1, nbmax
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=frhorho(igaux)-conjg(rhotwg(igp,ibvp))*rhotwg(ig,ibvp)
     enddo
   end do
 end do

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)

end subroutine fft4eet_0
!!***

!!****f* m_eet/fft4eet
!! NAME                  
!! fft4eet
!!
!! FUNCTION
!! Calculate f^{pp}, f^{pj}, and f{jj}. See Eqs. (25),(26),(27) of PRB 85, 085126, for their expressions.
!! They are stored in frhorho, frhoj, and fjj.
!! Terms involving the nonlocal part of the pseudopotential are not included.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine fft4eet(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax, &
&                  is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0, &
&                  gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                  itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,grottbm1, &
&                  nspinor,tim_fourdp,wfr1,wfr2,ibv,qplg,kplqg, &
&                  niter,frhorho,frhoj,fjj,spin_fact,qp_occ,qp_energy,chi0,igstart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfd_t),target,intent(inout) :: Wfs
 type(kmesh_t),intent(in) :: Kmesh
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(littlegroup_t),intent(in) :: Ltg_q
 integer,intent(in) :: ik_bz,nbhomo(2),nbmax,is,ibv,niter,igstart
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_k,itim_k,itim_kmq
 integer,intent(in) :: dim_rtwg
 real(dp),intent(in) :: spin_fact
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt
!arrays
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Ep%npwepG0)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: grottbm1(Ep%npwvec,2,Cryst%nsym)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qplg(3,Ep%npwe),kplqg(Ep%npwe)
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: wfr2(Wfs%nfftot*nspinor)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
 complex(gwpc),intent(out) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(out) :: frhoj(Ep%npwe,Ep%npwe)
 complex(gwpc),intent(out) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: ibvp,ig,igp,igaux,istwf_k,npw_k,i,j,ig4,ig4x,ig4y,ig4z,outofbox
 integer,save :: enough=0
 character(len=500) :: msg
 complex(gwpc) :: faux,minusone
!arrays
 integer :: gmgp(3)
 integer,pointer :: gbound_k(:,:),kg_k(:,:),igfft_k(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)
 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)
 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux(3)

!************************************************************************

 ABI_UNUSED((/isym_k,tim_fourdp/))

 istwf_k  =  Wfs%istwfk(ik_ibz)
 npw_k    =  Wfs%npwarr(ik_ibz)
 gbound_k => Wfs%Kdata(ik_ibz)%gbound
 kg_k     => Wfs%Kdata(ik_ibz)%kg_k    
 igfft_k  => Wfs%Kdata(ik_ibz)%igfft0

 ABI_ALLOCATE(rhotwg,(Ep%npwepG0*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Ep%npwepG0*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(npw_k))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
 endif

 if (niter>1)  then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
 end if
!
! Calculate matrix element <wfn2|exp[-iG.r]|wfn2>
!
 call calc_wfwfg(tabr_k,itim_k,nfftot_gw,ngfft_gw,wfr2,wfr2,wfwfg(:))

 do ibvp = 1, nbmax
!
! Calculate matrix element <wfn1|exp[-i(q+G).r]|wfn2>
!
   call rho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                spinrot_k,dim_rtwg,rhotwg(:,ibvp))
   if (ibvp>nbhomo(2)) then
!
! Add sum-over-states correction to the polarizability obtained with the EET
!
     call calc_corr_chi0(ik_bz,Ep,Kmesh,Gsph_epsG0,Ltg_q,rhotwg(:,ibvp),spin_fact,qp_occ,qp_energy, &
&                        ibv,ibvp,ik_ibz,ikmq_ibz,is,chi0)
   endif

 enddo

 if (niter>0) then
!
! Calculate FFT of i nabla u(r)
!
   ! MG: FIXME this won't work if I switch on k-centered G-spheres
! Here grottbm1 should be replaced by iskg_tabs_t
!
   call wfd_dur_isk(Wfs,Cryst,Kmesh,ibv,ik_bz,is,npw_k,grottbm1,dwfr)  !,ISkg,ik_ibz
!
! Calculate matrix element <wfn1|exp[-i(q+G).r]|nabla wfn2>
!
   do ibvp = 1, nbmax
     do i = 1, 3
       call drho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,dwfr(:,i), &
&                     dim_rtwg,drhotwg(:,ibvp,i))
     enddo
   enddo
!
! Calculate matrix element <wfn2|exp[-iG.r]|nabla wfn2> and <nabla wfn2|exp[-iG.r]|nabla wfn2>
!
   do i = 1, 3
     call calc_dwfwfg(tabr_k,itim_k,nfftot_gw,ngfft_gw,ph_mkt,wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(itim_k,nfftot_gw,ngfft_gw,dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo

   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Ep%npwe,nbmax))
!
! Calculate (q+G).nabla u(G)
!
   do ibvp = 1, nbmax
     do ig=1,Ep%npwe
       drhaux(:)=cmplx(real(drhotwg(ig,ibvp,:)),aimag(drhotwg(ig,ibvp,:)))
       cauxg(ig,ibvp)=vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
     enddo
   enddo
 endif
!
! Calculate <wfn2|exp[i(G-Gp).r]|wfn2>, <wfn2|exp[i(G-Gp).r]|nabla wfn2>.(q+Gp), and
! (q+G).<nabla wfn2|exp[i(G-Gp).r]|nabla wfn2>.(q+Gp)
! They are stored in frhorho, frhoj, and fjj.
!
 frhorho=czero_gw
 frhoj=czero_gw
 if (niter>1) fjj=czero_gw
 outofbox=0
 do ig=igstart,Ep%npwe
   do igp=igstart,Ep%npwe
     gmgp(:)=kg_k(:,ig)-kg_k(:,igp)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)
     if (ig<=igp) then
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       frhoj(ig,igp)=frhoj(ig,igp) + vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
       if (niter>1.and.ig<=igp) then
         igaux=igp*(igp-1)/2+ig
         paux(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux(i)=paux(i) + vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")
         enddo
         fjj(igaux)=fjj(igaux) + vdotw(qplg(:,ig),paux,Cryst%gmet,"G")
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if
!
! Calculate sum_v <wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|wfn2> and subtract from frhorho
!
 do ibvp = 1, nbmax
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=frhorho(igaux)-conjg(rhotwg(igp,ibvp))*rhotwg(ig,ibvp)
     enddo
   end do
 end do
!
! Calculate sum_v <wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|nabla wfn2>.(q+Gp) 
! and subtract from frhoj
!
 if (niter>0) then
   minusone=(-1.,0.)
   do ibvp = 1, nbmax
     call XGERC(Ep%npwe,Ep%npwe,minusone,cauxg(:,ibvp),1,rhotwg(:,ibvp),1,frhoj,Ep%npwe)
   end do
   do ig=igstart,Ep%npwe
     do igp=igstart,Ep%npwe
       if (ig<=igp) then
         igaux=igp*(igp-1)/2+ig
         faux=frhorho(igaux)
       else
         igaux=ig*(ig-1)/2+igp
         faux=conjg(frhorho(igaux))
       endif
       frhoj(ig,igp)=frhoj(ig,igp)+faux*kplqg(ig)
     end do !igp
   end do !ig
 end if
!
! Calculate sum_v (q+G).<nabla wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|nabla wfn2>.(q+Gp) 
! and subtract from fjj
!
 if (niter>1) then
   do ibvp = 1, nbmax
     do ig=igstart,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         fjj(igaux)=fjj(igaux)-conjg(cauxg(igp,ibvp))*cauxg(ig,ibvp)
       enddo
     end do
   end do
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       faux=frhorho(igaux)
       fjj(igaux)=fjj(igaux)+kplqg(igp)*frhoj(ig,igp)+ kplqg(ig)*conjg(frhoj(igp,ig)) - &
&                            kplqg(igp)*faux*kplqg(ig)
     end do !igp
   end do !ig
 endif

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
 endif
 if (niter>1)  then
   ABI_DEALLOCATE(ddwfwfg)
 end if

end subroutine fft4eet
!!***

!!****f* m_eet/fft4eet_kb
!! NAME                  
!! fft4eet_kb
!!
!! FUNCTION
!! Calculate f^{pp}, f^{pj}, and f{jj}. See Eqs. (25),(26),(27) of PRB 85, 085126, for their expressions.
!! They are stored in frhorho, frhoj, and fjj.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine fft4eet_kb(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,Psps,nbhomo,nbmax, &
&                     is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0, &
&                     gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                     itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,grottbm1, &
&                     nspinor,tim_fourdp,fnlloc,fnlmax,fnlkpr,mtwk,mtwkp,wfr1,wfr2,ibv,qplg,kplqg, &
&                     niter,frhorho,frhoj,fjj,spin_fact,qp_occ,qp_energy,chi0,igstart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_kb'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfd_t),target,intent(inout) :: Wfs
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(littlegroup_t),intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz,nbhomo(2),nbmax,is,ibv,niter,igstart
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_k,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: dim_rtwg
 integer,intent(in) :: igfftepsG0(Ep%npwepG0)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: grottbm1(Ep%npwvec,2,Cryst%nsym)
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: spin_fact
 real(dp),intent(in) :: qplg(3,Ep%npwe),kplqg(Ep%npwe)
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt

 complex(gwpc),intent(in) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(in) :: mtwk(Wfs%nfftot*nspinor,nbhomo(1))
 complex(gwpc),intent(in) :: mtwkp(Wfs%nfftot*nspinor,nbmax)

 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: wfr2(Wfs%nfftot*nspinor)

 complex(gwpc),intent(out) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(out) :: frhoj(Ep%npwe,Ep%npwe)
 complex(gwpc),intent(out) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: istwf_k,npw_k,nlx,outofbox
 integer :: ilm,iat,ibvp,ig,igp,igaux,i,j
 integer :: ig4,ig4x,ig4y,ig4z
 integer :: ig5,ig5x,ig5y,ig5z
 integer,save :: enough=0
 character(len=500) :: msg
 complex(gwpc) :: faux,minusone
!arrays
 integer :: gmgp(3),ngfft(3)
 integer,pointer :: gbound_k(:,:),kg_k(:,:),igfft_k(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg2(:,:)
 complex(gwpc),allocatable :: fnltwg3(:,:)
 complex(gwpc),allocatable :: kns(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)
 complex(gwpc),allocatable :: fnlwfg(:)
 complex(gwpc),allocatable :: fkdwfg(:,:)
 complex(gwpc),allocatable :: fdrhotwg(:,:,:,:)
 complex(gwpc),allocatable :: lnkp(:)
 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)
 complex(gwpc),allocatable :: ff(:,:,:,:)
 complex(gwpc),allocatable :: vzn(:,:,:)
 complex(gwpc),allocatable :: paux(:,:)
 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux2(3)

!************************************************************************

 ABI_UNUSED((/isym_k,tim_fourdp/))

 istwf_k  =  Wfs%istwfk(ik_ibz)
 npw_k    =  Wfs%npwarr(ik_ibz)
 gbound_k => Wfs%Kdata(ik_ibz)%gbound
 kg_k     => Wfs%Kdata(ik_ibz)%kg_k    
 igfft_k  => Wfs%Kdata(ik_ibz)%igfft0

 nlx = min(Psps%mpsang,4)

 ABI_ALLOCATE(rhotwg,(Ep%npwepG0*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Ep%npwepG0*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(npw_k))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
   ABI_ALLOCATE(fnltwg,(Ep%npwepG0*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlwfg,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(fnltwg2,(Ep%npwepG0*nspinor**2,nbmax))
   ABI_ALLOCATE(fnltwg3,(Ep%npwepG0*nspinor**2,nbmax))
 endif

 if (niter>1) then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
   ABI_ALLOCATE(lnkp,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(kns,(Ep%npwepG0*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fdrhotwg,(Ep%npwepG0*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom,3))
   ABI_ALLOCATE(fkdwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(ff,(nlx*nlx,Cryst%natom,nlx*nlx,Cryst%natom))
   ABI_ALLOCATE(vzn,(Ep%npwepG0*nspinor**2,nlx*nlx,Cryst%natom))
 endif
!
! Calculate matrix element <wfn2|exp[-iG.r]|wfn2>
!
 call calc_wfwfg(tabr_k,itim_k,nfftot_gw,ngfft_gw,wfr2,wfr2,wfwfg)

 do ibvp = 1, nbmax
!
! Calculate matrix element <wfn1|exp[-i(q+G).r]|wfn2>
!
   call rho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                spinrot_k,dim_rtwg,rhotwg(:,ibvp))

   if (ibvp>nbhomo(2)) then
!
! Add sum-over-states correction to the polarizability obtained with the EET
!
     call calc_corr_chi0(ik_bz,Ep,Kmesh,Gsph_epsG0,Ltg_q,rhotwg(:,ibvp),spin_fact,qp_occ,qp_energy, &
&                        ibv,ibvp,ik_ibz,ikmq_ibz,is,chi0)
   endif
 enddo

 if (niter>0) then
!
! Calculate FFT of i nabla u(r)
!
! MG: FIXME this won't work if I switch on k-centered G-spheres
! Here grottbm1 should be replaced by iskg_tabs_t
!
   call wfd_dur_isk(Wfs,Cryst,Kmesh,ibv,ik_bz,is,npw_k,grottbm1,dwfr)  !,ISkg,ik_ibz
!
! Calculate matrix elements <wfn1|exp[-i(q+G).r]|nabla wfn2>
!
   do ibvp = 1, nbmax
     do i = 1, 3
       call drho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,dwfr(:,i), &
&                     dim_rtwg,drhotwg(:,ibvp,i))
     enddo
   enddo
!
! Calculate matrix element <wfn2|exp[-i(q+G).r]|nabla wfn2> and <nabla wfn2|exp[-iG.r]|nabla wfn2>
!
   do i = 1, 3
     call calc_dwfwfg(tabr_k,itim_k,nfftot_gw,ngfft_gw,ph_mkt,wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(itim_k,nfftot_gw,ngfft_gw,dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo
!
! Calculate matrix element <fnl|exp[-i(q+G).r]|wfn2>
!
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       call rho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                    fnlkpr(:,ilm,iat),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                    spinrot_k,dim_rtwg,fnltwg(:,ilm,iat))
     enddo
   enddo
!
! Calculate matrix element <wfn2|exp[-iG.r]|Mtw>
!
   call calc_wfwfg(tabr_k,itim_k,nfftot_gw,ngfft_gw,wfr2,mtwk(:,ibv),fnlwfg)
!
! Calculate matrix elements <wfn1|exp[-i(q+G).r]|Mtw_kp> and <Mtw_k|exp[-i(q+G).r]|wfn2>
!
   do ibvp = 1, nbmax
     call rho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,mtwk(:,ibv),itim_k,tabr_k,ph_mkt, &
&                  spinrot_k,dim_rtwg,fnltwg2(:,ibvp))

     call rho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  mtwkp(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                  spinrot_k,dim_rtwg,fnltwg3(:,ibvp))
  
   enddo
!
! Calculate matrix elements <fnl|exp[-i(q+G).r]|nabla wfn2> ,<fnl|exp[-i(q+G).r]|Mtw>, 
! <fnl|exp[-iG.r]|Mtw> and <Mtw|exp[-iG.r]|Mtw>
!
   if (niter>1) then
     do iat = 1, Cryst%natom
       do ilm = 1, Psps%mpsang*Psps%mpsang
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         do i = 1, 3
           call drho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                         fnlkpr(:,ilm,iat),itim_kmq,tabr_kmq,ph_mkmqt, &
&                         dwfr(:,i),dim_rtwg,fdrhotwg(:,ilm,iat,i))
         enddo
         call rho_tw_g(nspinor,Ep%npwepG0,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                       fnlkpr(:,ilm,iat),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,mtwk(:,ibv),itim_k,tabr_k,ph_mkt, &
&                       spinrot_k,dim_rtwg,kns(:,ilm,iat))
       enddo
     enddo
     do i = 1, 3
       call calc_dwfwfg(tabr_k,itim_k,nfftot_gw,ngfft_gw,ph_mkt,mtwk(:,ibv),dwfr(:,i),fkdwfg(:,i))
     enddo
     call calc_wfwfg(tabr_k,itim_k,nfftot_gw,ngfft_gw,mtwk(:,ibv),mtwk(:,ibv),lnkp)
   endif

   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Ep%npwe,nbmax))
!
! Calculate (q+G).nabla u(G)-fnltwg2+fnltwg3
!
   do ibvp = 1, nbmax
     do ig=1,Ep%npwe
       drhaux(:)=cmplx(real(drhotwg(ig,ibvp,:)),aimag(drhotwg(ig,ibvp,:)))
       cauxg(ig,ibvp)=vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
     enddo
   enddo
   do ig=1,Ep%npwe
     cauxg(ig,:)=cauxg(ig,:)-fnltwg2(ig,:)+fnltwg3(ig,:)
   enddo
 endif
!
! Calculate <wfn2|exp[i(G-Gp).r]|wfn2>, <wfn2|exp[i(G-Gp).r]|nabla wfn2>.(q+Gp)-fnlwfg, and
! (q+G).<nabla wfn2|exp[i(G-Gp).r]|nabla wfn2>.(q+Gp)-(q+G).fkdwfg+lnkp
! They are stored in frhorho, frhoj, and fjj.
!
 frhorho=czero_gw
 frhoj=czero_gw
 if (niter>1) fjj=czero_gw
 outofbox=0
 do ig=igstart,Ep%npwe
   do igp=igstart,Ep%npwe
     gmgp(:)=Gsph_epsG0%gvec(:,ig)-Gsph_epsG0%gvec(:,igp)
     ngfft(1)=ngfft_gw(1)
     ngfft(2)=ngfft_gw(2)
     ngfft(3)=ngfft_gw(3)
     if (ANY(gmgp(:)>ngfft(1:3)/2) .or. ANY(gmgp(:)<-(ngfft(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft(1))
     ig4y= modulo(gmgp(2),ngfft(2))
     ig4z= modulo(gmgp(3),ngfft(3))
     ig4= 1+ig4x+ig4y*ngfft(1)+ig4z*ngfft(1)*ngfft(2)
     ig5x= modulo(-gmgp(1),ngfft(1))
     ig5y= modulo(-gmgp(2),ngfft(2))
     ig5z= modulo(-gmgp(3),ngfft(3))
     ig5= 1+ig5x+ig5y*ngfft(1)+ig5z*ngfft(1)*ngfft(2)
     if (ig<=igp) then
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       frhoj(ig,igp)=frhoj(ig,igp) + vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")-fnlwfg(ig4)
       if (niter>1.and.ig<=igp) then
         igaux=igp*(igp-1)/2+ig
         paux2(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux2(i)=paux2(i) + vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")
         enddo
         drhaux(:)=paux2(:)-cmplx(real(fkdwfg(ig4,:)),aimag(fkdwfg(ig4,:)))
         fjj(igaux)=fjj(igaux) + vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
         drhaux(:)=cmplx(real(fkdwfg(ig5,:)),-aimag(fkdwfg(ig5,:)))
         fjj(igaux)=fjj(igaux) - vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")+lnkp(ig4)
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if
!
! Calculate sum_v <wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|wfn2> and subtract from frhorho
!
 do ibvp = 1, nbmax
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=frhorho(igaux)-conjg(rhotwg(igp,ibvp))*rhotwg(ig,ibvp)
     enddo
   end do
 end do
!
! Calculate sum_v <wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|nabla wfn2>.(q+Gp) + fnltw(G)^* fnltw(Gp)
! and subtract from frhoj
!
 if (niter>0) then
   ABI_ALLOCATE(paux,(Ep%npwe,Ep%npwe))
   paux(:,:)=(0.0,0.0)
   minusone=(-1.,0.)
   do ibvp = 1, nbmax
     call XGERC(Ep%npwe,Ep%npwe,minusone,cauxg(:,ibvp),1,rhotwg(:,ibvp),1,frhoj,Ep%npwe)
   enddo
   do ig=igstart,Ep%npwe
     do igp=igstart,Ep%npwe
       if (ig>=igp) then
         do iat = 1, Cryst%natom
           do ilm = 1, nlx*nlx
             if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
             if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
             paux(ig,igp)=paux(ig,igp)+conjg(fnltwg(igp,ilm,iat))*fnltwg(ig,ilm,iat)
           enddo
         enddo
       endif
       if (ig<=igp) then
         igaux=igp*(igp-1)/2+ig
         faux=frhorho(igaux)
       else
         igaux=ig*(ig-1)/2+igp
         faux=conjg(frhorho(igaux))
       endif
       frhoj(ig,igp)=frhoj(ig,igp)+faux*kplqg(ig)
     end do !igp
   end do !ig
   do ig=igstart,Ep%npwe
     do igp=igstart,Ep%npwe
       if (ig>=igp) then
         frhoj(ig,igp)=frhoj(ig,igp)+paux(ig,igp)
       else
         frhoj(ig,igp)=frhoj(ig,igp)+conjg(paux(igp,ig))
       endif
     end do !igp
   end do !ig
   ABI_DEALLOCATE(paux)
 end if
!
! Calculate sum_v (q+G).<nabla wfn2|exp[iG.r]|wfn1><wfn1|exp[-iGp.r]|nabla wfn2>.(q+Gp) + 
! ff fnltw(G)^* fnltw(Gp) + fnltw(G)^* fdrhotwg(Gp).(q+Gp) - kns(G) fnltw(Gp)
! and subtract from fjj
!
 if (niter>1) then
!
! Calculate sum_G fnl_kp(G) fnl^*_kp(G) and store in ff
!
!
! FIXME: Here grottbm1 should be replaced by iskg_tabs_t
   call compute_ff(Cryst,Kmesh,Psps,ik_bz,nlx,istwf_k,npw_k,kg_k,grottbm1,fnlmax,fnlloc,ff)

   vzn(:,:,:)=czero_gw
   do ig=igstart,Ep%npwe
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         vzn(ig,:,:) = vzn(ig,:,:) + half*ff(ilm,iat,:,:)*fnltwg(ig,ilm,iat)
       enddo
     enddo
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         drhaux(:)=cmplx(real(fdrhotwg(ig,ilm,iat,:)),aimag(fdrhotwg(ig,ilm,iat,:)))
         vzn(ig,ilm,iat)=vzn(ig,ilm,iat)+vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")-kns(ig,ilm,iat)
       enddo
     enddo
   enddo
   do ibvp = 1, nbmax
     do ig=igstart,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         fjj(igaux)=fjj(igaux)-conjg(cauxg(igp,ibvp))*cauxg(ig,ibvp)
       enddo
     enddo
   enddo
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       faux=frhorho(igaux)
       fjj(igaux)=fjj(igaux)+kplqg(igp)*frhoj(ig,igp) + kplqg(ig)*conjg(frhoj(igp,ig))&
&                           -kplqg(igp)*faux*kplqg(ig)
       do iat = 1, Cryst%natom
         do ilm = 1, nlx*nlx
           if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
           if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
           fjj(igaux)=fjj(igaux)+conjg(fnltwg(igp,ilm,iat))*vzn(ig,ilm,iat)+ &
&                                      fnltwg(ig,ilm,iat)*conjg(vzn(igp,ilm,iat))
         enddo
       enddo
     end do !igp
   end do !ig

 end if

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
   ABI_DEALLOCATE(fnltwg)
   ABI_DEALLOCATE(fnltwg2)
   ABI_DEALLOCATE(fnltwg3)
   ABI_DEALLOCATE(fnlwfg)
 endif
 if (niter>1) then
   ABI_DEALLOCATE(ddwfwfg)
   ABI_DEALLOCATE(lnkp)
   ABI_DEALLOCATE(kns)
   ABI_DEALLOCATE(fdrhotwg)
   ABI_DEALLOCATE(fkdwfg)
   ABI_DEALLOCATE(vzn)
   ABI_DEALLOCATE(ff)
 endif

end subroutine fft4eet_kb
!!***

!!****f* m_eet/calc_delta_chi0
!!
!! NAME
!! calc_delta_chi0
!!
!! FUNCTION
!!      Calculation of frhoj/frhorho and fjj/frhoj which are necessary 
!!      for the calculation of the effective energy.
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_delta_chi0(Ep,frhorho,frhoj,fjj,niter)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_delta_chi0'
!End of the abilint section

 implicit none

 type(Epsilonm1_parameters),intent(in) :: Ep

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter
!arrays
 complex(gwpc),intent(in) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(inout) :: frhoj(Ep%npwe,Ep%npwe)
 complex(gwpc),intent(inout) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))

!Local variables-------------------------------
!scalars

 integer :: ig,igp,igaux!,ios
! integer :: iter

! real(dp) :: bq
 real(dp) :: delta_huge,denchk

 complex(gwpc) :: num,den

!*************************************************************************

   delta_huge = 1.0d8

   if (niter==1) then
!
! Calculation of frhoj/frhorho
!
     do ig=1,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         num = frhoj(ig,igp)+conjg(frhoj(igp,ig))
         den = 2*frhorho(igaux)
         denchk=abs(real(den))+abs(aimag(den))
         if (denchk>0.0) then
           frhoj(ig,igp)=num/den
         else
           frhoj(ig,igp)=cmplx(delta_huge,0.0)
         endif
       end do !igp
     end do !ig

     do ig=1,Ep%npwe
       do igp=1,ig-1
         frhoj(ig,igp)=conjg(frhoj(igp,ig))
       end do !igp
     end do !ig

   elseif (niter==2) then
!
! Calculation of frhoj/frhorho and frhoj/fjj
!
     do ig=1,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         num = 2*fjj(igaux)
         den = frhoj(ig,igp)+conjg(frhoj(igp,ig))
         denchk=abs(real(den))+abs(aimag(den))
         if (denchk>0.0) then
           fjj(igaux)=num/den
         else
           fjj(igaux)=cmplx(delta_huge,0.0)
         endif
       end do !igp
     end do !ig

     do ig=1,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         num = frhoj(ig,igp)+conjg(frhoj(igp,ig))
         den = 2*frhorho(igaux)
         denchk=abs(real(den))+abs(aimag(den))
         if (denchk>0.0) then
           frhoj(ig,igp)=num/den
         else
           frhoj(ig,igp)=cmplx(delta_huge,0.0)
         endif
       end do !igp
     end do !ig

     do ig=1,Ep%npwe
       do igp=1,ig-1
         frhoj(ig,igp)=conjg(frhoj(igp,ig))
       end do !igp
     end do !ig

   endif

end subroutine calc_delta_chi0
!!***

!!****f* m_eet/calc_chi0_delta_clos
!!
!! NAME
!! calc_chi0_delta_clos
!!
!! FUNCTION
!! Calculation of the part of the polarizability coming from
!! the sum over all the states using the EET
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_chi0_delta_clos(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,niter,epsv,epslumo,qpgsq,frhorho,frhoj,fjj,chi0)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_chi0_delta_clos'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars

 type(Dataset_type),intent(in) :: Dtset
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(gsphere_t),target,intent(in) :: Gsph_epsG0
 type(littlegroup_t),target,intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz
 integer,intent(in) :: niter
!arrays
 real(dp),intent(in) :: qpgsq(Ep%npwe)
 real(dp),intent(in) :: epsv,epslumo
 complex(gwpc),intent(in) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(in) :: frhoj(Ep%npwe,Ep%npwe)
 complex(gwpc),intent(in) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,igaux,igaux2,ios,isym,itim
 real(dp) :: bq, delta_huge
 complex(gwpc) :: delta1,delta2,den1,den2,faux
 complex(gwpc) :: frhorho_sym,frhoj_sym,fjj_sym
 character(len=500) :: msg
!arrays
 integer,pointer :: gmG0(:)
 integer,allocatable :: Sm1_gmG0(:)
 complex(gwpc),pointer :: phmGt(:)

!*************************************************************************

 delta_huge=1.0d8

 SELECT CASE (Ep%symchi)

 CASE (0)

   do ig=1,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       bq=half*(qpgsq(ig)+qpgsq(igp))
       faux=frhorho(igaux)

       SELECT CASE (niter)
       CASE (1)
         delta1=bq+frhoj(ig,igp)
         if (Dtset%gw_eet_scale>0.01) delta1=Dtset%gw_eet_scale*delta1
         delta2=delta1
         do ios=1,Ep%nomega
           if (abs(aimag(Ep%omega(ios)))<0.001) then
             call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
           else
             delta1=delta2
           endif
           chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta1) &
                                               + faux/(Ep%omega(ios)-delta1)
         end do !ios
       CASE (2)
         do ios=1,Ep%nomega
           delta1=delta_huge
           delta2=delta_huge
           den1=Ep%omega(ios)+bq+fjj(igaux)
           den2=Ep%omega(ios)-bq-fjj(igaux)
           if(abs(den1)>tol12) delta1=bq+frhoj(ig,igp)*(Ep%omega(ios)+bq+frhoj(ig,igp))/den1
           if(abs(den2)>tol12) delta2=bq+frhoj(ig,igp)*(Ep%omega(ios)-bq-frhoj(ig,igp))/den2
           if (Dtset%gw_eet_scale>0.01) then
             delta1=Dtset%gw_eet_scale*delta1
             delta2=Dtset%gw_eet_scale*delta2
           endif
           if (abs(aimag(Ep%omega(ios)))<0.001) then
             call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
             call check_delta(Ep,qpgsq,delta2,epsv,epslumo,ig,igp)
           endif
           chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta1) &
                                               + faux/(Ep%omega(ios)-delta2)
         end do !ios
       END SELECT

     end do !igp
   end do !ig

 CASE (1)

   ABI_ALLOCATE(Sm1_gmG0,(Ep%npwe))

   do isym=1,Ltg_q%nsym_sg
     do itim = 1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then

         phmGt => Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         gmG0  => Ltg_q%igmG0(1:Ep%npwe,itim,isym)
         Sm1_gmG0(1:Ep%npwe)=Gsph_epsG0%rottbm1(gmG0(1:Ep%npwe),itim,isym)

         do ig = 1, Ep%npwe
           do igp = ig, Ep%npwe
             igaux=igp*(igp-1)/2+ig
             bq=half*(qpgsq(ig)+qpgsq(igp))
             if (Sm1_gmG0(ig)<=Sm1_gmG0(igp)) then
               igaux2=Sm1_gmG0(igp)*(Sm1_gmG0(igp)-1)/2+Sm1_gmG0(ig)
               frhorho_sym=frhorho(igaux2)
               if (niter==2) fjj_sym=fjj(igaux2)
             else
               igaux2=Sm1_gmG0(ig)*(Sm1_gmG0(ig)-1)/2+Sm1_gmG0(igp)
               frhorho_sym=conjg(frhorho(igaux2))
               if (niter==2) fjj_sym=conjg(fjj(igaux2))
             endif
             frhoj_sym=frhoj(Sm1_gmG0(ig),Sm1_gmG0(igp))
             if (itim==2) then
               frhorho_sym=conjg(frhorho_sym)
               frhoj_sym=conjg(frhoj_sym)
               if (niter==2) fjj_sym=conjg(fjj_sym)
             endif
             frhorho_sym=frhorho_sym*conjg(phmGt(igp))*phmGt(ig)

             SELECT CASE (niter)
             CASE (1)
               delta1=bq+frhoj_sym
               if (Dtset%gw_eet_scale>0.01) delta1=Dtset%gw_eet_scale*delta1
               delta2=delta1
               do ios=1,Ep%nomega
                 if (abs(aimag(Ep%omega(ios)))<0.001) then
                   call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
                 else
                   delta1=delta2
                 endif
                 chi0(ig,igp,ios) = chi0(ig,igp,ios) - frhorho_sym/(Ep%omega(ios)+delta1) &
                                                     + frhorho_sym/(Ep%omega(ios)-delta1)
               enddo
             CASE (2)
               do ios=1,Ep%nomega
                 delta1=bq+frhoj_sym*(Ep%omega(ios)+bq+frhoj_sym)/(Ep%omega(ios)+bq+fjj_sym)
                 delta2=bq+frhoj_sym*(Ep%omega(ios)-bq-frhoj_sym)/(Ep%omega(ios)-bq-fjj_sym)
                 if (Dtset%gw_eet_scale>0.01) then
                   delta1=Dtset%gw_eet_scale*delta1
                   delta2=Dtset%gw_eet_scale*delta2
                 endif
                 if (abs(aimag(Ep%omega(ios)))<0.001) then
                   call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
                   call check_delta(Ep,qpgsq,delta2,epsv,epslumo,ig,igp)
                 endif
                 chi0(ig,igp,ios) = chi0(ig,igp,ios) - frhorho_sym/(Ep%omega(ios)+delta1) & 
                                                     + frhorho_sym/(Ep%omega(ios)-delta2)
               enddo
             END SELECT

           end do !igp
         end do !ig

       endif
     enddo
   enddo

   ABI_DEALLOCATE(Sm1_gmG0)

 CASE DEFAULT
  write(msg,'(a,i3)')'Wrong symchi= ',Ep%symchi
  MSG_BUG(msg)
 END SELECT

end subroutine calc_chi0_delta_clos
!!***

!!****f* m_eet/calc_chi0_delta0
!!
!! NAME
!! calc_chi0_delta0
!!
!! FUNCTION
!! Calculation of the part of the polarizability coming from
!! the sum over all the states using the EET in its zero-order approximation
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_chi0_delta0(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,epsv,epslumo,qpgsq,frhorho,chi0)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_chi0_delta0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz
 type(Dataset_type),intent(in) :: Dtset
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(gsphere_t),target,intent(in) :: Gsph_epsG0
 type(littlegroup_t),target,intent(in) :: Ltg_q

!arrays
 real(dp),intent(in) :: qpgsq(Ep%npwe)
 real(dp),intent(in) :: epsv,epslumo
 complex(gwpc),intent(in) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ios,isym,itim,igaux,igaux2
 complex(gwpc) :: faux,delta1,delta2,frhorho_sym
 character(len=500) :: msg
!arrays
 integer,pointer :: gmG0(:)
 integer,allocatable :: Sm1_gmG0(:)
 complex(gwpc),pointer :: phmGt(:)

!*************************************************************************

 SELECT CASE (Ep%symchi)

 CASE (0)

   do ig=1,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       faux=frhorho(igaux)
       delta1=half*(qpgsq(ig)+qpgsq(igp))
       if (Dtset%gw_eet_scale>0.01) delta1=Dtset%gw_eet_scale*delta1
       delta2=delta1
       do ios=1,Ep%nomega
         if (abs(aimag(Ep%omega(ios)))<0.001) then
           call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
         else
           delta1=delta2
         endif
         chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta1) &
                                             + faux/(Ep%omega(ios)-delta1)
       end do !ios
     end do !igp
   end do !ig

 CASE (1)

   ABI_ALLOCATE(Sm1_gmG0,(Ep%npwe))

   do isym=1,Ltg_q%nsym_sg
     do itim = 1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then

         phmGt => Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         gmG0  => Ltg_q%igmG0(1:Ep%npwe,itim,isym)
         Sm1_gmG0(1:Ep%npwe)=Gsph_epsG0%rottbm1(gmG0(1:Ep%npwe),itim,isym)

         do ig = 1, Ep%npwe
           do igp = ig, Ep%npwe
             igaux=igp*(igp-1)/2+ig
             if (Sm1_gmG0(ig)<=Sm1_gmG0(igp)) then
               igaux2=Sm1_gmG0(igp)*(Sm1_gmG0(igp)-1)/2+Sm1_gmG0(ig)
               faux=frhorho(igaux2)
             else
               igaux2=Sm1_gmG0(ig)*(Sm1_gmG0(ig)-1)/2+Sm1_gmG0(igp)
               faux=conjg(frhorho(igaux2))
             endif
             frhorho_sym=faux*conjg(phmGt(igp))*phmGt(ig)
             if (itim==2) frhorho_sym=conjg(frhorho_sym)
             delta1=half*(qpgsq(ig)+qpgsq(igp))
             if (Dtset%gw_eet_scale>0.01) delta1=Dtset%gw_eet_scale*delta1
             delta2=delta1
             do ios=1,Ep%nomega
               if (abs(aimag(Ep%omega(ios)))<0.001) then
                 call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
               else
                 delta1=delta2
               endif
               chi0(ig,igp,ios) = chi0(ig,igp,ios) - frhorho_sym/(Ep%omega(ios)+delta1) &
                                                   + frhorho_sym/(Ep%omega(ios)-delta1)
             end do !ios
           end do !igp
         end do !ig

       endif
     enddo
   enddo

   ABI_DEALLOCATE(Sm1_gmG0)

 CASE DEFAULT
  write(msg,'(a,i3)')'Wrong symchi= ',Ep%symchi
  MSG_BUG(msg)
 END SELECT

end subroutine calc_chi0_delta0
!!***

!!****f* m_eet/calc_corr_chi0
!!
!! NAME
!! calc_corr_chi0
!!
!! FUNCTION
!! Calculate sum-over-states correction to the polarizability obtained with the EET
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_corr_chi0(ik_bz,Ep,Kmesh,Gsph_epsG0,Ltg_q,rhotwg,spin_fact,qp_occ,qp_energy, &
&                         ibv,ibc,ik_ibz,ikmq_ibz,is,chi0)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_corr_chi0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(kmesh_t),intent(in) :: Kmesh
 type(gsphere_t),target,intent(in) :: Gsph_epsG0
 type(littlegroup_t),target,intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz,ibv,ibc,ik_ibz,ikmq_ibz,is
 real(dp),intent(in) :: spin_fact
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 complex(gwpc),intent(in) :: rhotwg(Ep%npwepG0)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

 integer :: ig,igp,ios,isym,itim
 real(dp) :: deltae,deltaf
 complex(dpc) :: green_w
 complex(gwpc) :: dd

 character(len=500) :: msg
!arrays
 integer,pointer :: gmG0(:)
 integer,allocatable :: Sm1_gmG0(:)
 complex(gwpc),allocatable :: rhotwg_sym(:)
 complex(gwpc),pointer :: phmGt(:)

!************************************************************************

 deltaf=spin_fact*(qp_occ(ibc,ikmq_ibz,is)-qp_occ(ibv,ik_ibz,is))
 deltae=qp_energy(ibc,ikmq_ibz,is)-qp_energy(ibv,ik_ibz,is)

 SELECT CASE (Ep%symchi)

 CASE (0) ! Do not use symmetries

   do ios = 1, Ep%nomega
     if (ibc==ibv) then
       green_w = g0g0w(Ep%omega(ios),deltaf,deltae,Ep%zcut,GW_TOL_W0,1)
     else
       green_w = g0g0w(Ep%omega(ios),deltaf,deltae,Ep%zcut,GW_TOL_W0,2)
     endif
     dd = green_w
     if (ABS(REAL(Ep%omega(ios)))<0.00001) then
       do ig=1,Ep%npwe
         do igp=ig,Ep%npwe
           chi0(ig,igp,ios) = chi0(ig,igp,ios)+rhotwg(ig)*conjg(rhotwg(igp))*dd
         end do !igp
       end do !ig
     else
       call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,ios),Ep%npwe)
     endif
   end do !ios

 CASE (1) ! Use symmetries to reconstruct the integrand in the BZ.

   ABI_ALLOCATE(rhotwg_sym,(Ep%npwe))
   ABI_ALLOCATE(Sm1_gmG0  ,(Ep%npwe))

   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then

        phmGt => Gsph_epsG0%phmGt(1:Ep%npwe,isym)
        gmG0  => Ltg_q%igmG0     (1:Ep%npwe,itim,isym)
        Sm1_gmG0(1:Ep%npwe)=Gsph_epsG0%rottbm1(gmG0(1:Ep%npwe),itim,isym)

        SELECT CASE (itim)
        CASE (1)
          rhotwg_sym(1:Ep%npwe)=rhotwg(Sm1_gmG0)*phmGt(1:Ep%npwe)
        CASE (2)
          rhotwg_sym(1:Ep%npwe)=CONJG(rhotwg(Sm1_gmG0))*phmGt(1:Ep%npwe)
        CASE DEFAULT
          write(msg,'(a,i3)')'Wrong itim= ',itim
          MSG_BUG(msg)
        END SELECT

        do ios=1, Ep%nomega
          if (ibc==ibv) then
            green_w = g0g0w(Ep%omega(ios),deltaf,deltae,Ep%zcut,GW_TOL_W0,1)
          else
            green_w = g0g0w(Ep%omega(ios),deltaf,deltae,Ep%zcut,GW_TOL_W0,2)
          endif
          dd = green_w
          if (ABS(REAL(Ep%omega(ios)))<0.00001) then
            do ig=1,Ep%npwe
              do igp=ig,Ep%npwe
                chi0(ig,igp,ios) = chi0(ig,igp,ios)+rhotwg_sym(ig)*conjg(rhotwg_sym(igp))*dd
              end do !igp
            end do !ig
          else
            call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,ios),Ep%npwe)
          endif
        end do

       end if
     end do
   end do

   ABI_DEALLOCATE(rhotwg_sym)
   ABI_DEALLOCATE(Sm1_gmG0)

 CASE DEFAULT
   write(msg,'(a,i3)')'Wrong symchi= ',Ep%symchi
   MSG_BUG(msg)
 END SELECT

end subroutine calc_corr_chi0
!!***

!!****f* m_eet/check_delta
!!
!! NAME
!! check_delta
!!
!! FUNCTION
!! Check if the effective energy (delta) is in the correct range
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine check_delta(Ep,qpgsq,delta,epsv,epslumo,ig,igp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_delta'
!End of the abilint section

 implicit none

 type(Epsilonm1_parameters),intent(in) :: Ep

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: ig,igp
 real(dp),intent(in) :: qpgsq(Ep%npwe)
 real(dp),intent(in) :: epsv,epslumo
 complex(gwpc),intent(inout) :: delta

!Local variables-------------------------------
!scalars

 real(dp) :: test

!*************************************************************************

 if (ig==igp) then
   delta=real(delta)
   test=delta+epsv
   if (test<epslumo) then
     delta=half*(qpgsq(ig)+qpgsq(igp))
     test=delta+epsv
   endif
   if (test<epslumo) then
     delta=epslumo-epsv
   endif
 else
   test=delta+epsv
   if (test<epslumo) then
     delta=epslumo-epsv
   endif
 endif

end subroutine check_delta
!!***

!!****f* m_eet/calc_sig_ppm_delta_clos_cd
!! NAME
!! calc_sig_ppm_delta_clos_cd
!!
!! FUNCTION
!! Calculation of the part of the matrix elements of the self-energy coming from
!! the sum over all the states in OPTIMAL GW
!! 
!! INPUTS
!! 
!! OUTPUT
!! 
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_sig_ppm_delta_clos_cd(Er,npwc,nomega,ikbz,jkbz,qbzpg,omegame0k,omegame0lumo,ptwsq,epsm1_qbz, &
&                                     ket,gw_eet_scale,niter)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sig_ppm_delta_clos_cd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_results),intent(in) :: Er
 integer,intent(in) :: nomega,npwc,niter
 integer,intent(in) :: ikbz,jkbz
 real(dp),intent(in) :: gw_eet_scale
!arrays
 real(dp),intent(in) :: omegame0k(nomega),omegame0lumo(nomega)
 real(dp),intent(in) :: qbzpg(npwc)
 complex(gwpc),intent(in) :: ptwsq(npwc,npwc,niter+1)
 complex(gwpc),intent(in) :: epsm1_qbz(npwc,npwc,Er%nomega)
 complex(dpc),intent(inout) :: ket(nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ios,io,iaux
 real(dp) :: bq
 complex(gwpc) :: num,delta
 logical :: ggpzero
!arrays
 real(dp), allocatable :: qpgsq(:)
 real(dp) :: omegame0k_tmp(nomega)
 real(dp) :: left(nomega),right(nomega)
 complex(gwpc) :: weight(Er%nomega_i+1,nomega)
 complex(dpc) :: omega_imag(Er%nomega_i+1)
 complex(dpc) :: domegaleft,domegaright

!*************************************************************************

 ABI_ALLOCATE(qpgsq,(npwc))

 do ig = 1, npwc
   qpgsq(ig)=half*qbzpg(ig)*qbzpg(ig)
 enddo

 ! Frequency mesh for integral along the imaginary axis.    
 omega_imag(1)=Er%omega(1)
 omega_imag(2:Er%nomega_i+1)=Er%omega(Er%nomega_r+1:Er%nomega)

 do ig = 1, npwc
   do igp = 1, npwc

     bq=half*(qpgsq(ig)+qpgsq(igp))
     if (niter==0) then
       delta = bq
     elseif(niter==1) then
       delta = bq+ptwsq(ig,igp,2)
     endif
     ggpzero=(ig==1.and.igp==1)
     if (ikbz==jkbz.and.ggpzero) delta = 0.0_gwp

     do ios = 1, nomega
       call check_delta_sigma(qpgsq(ig),qpgsq(igp),delta,omegame0k(ios),omegame0lumo(ios), &
&                             ig,igp,gw_eet_scale,niter)
       omegame0k_tmp(ios)=omegame0k(ios) - delta
       ! Avoid divergences in $\omega - \omega_s$.
       if (ABS(omegame0k_tmp(ios))<tol6) omegame0k_tmp(ios)=sign(tol6,omegame0k_tmp(ios))
     end do

     weight(1,:) = ATAN(-half*AIMAG(omega_imag(2))/REAL(omegame0k_tmp(:)))
     domegaleft  = (three*omega_imag(Er%nomega_i+1)-omega_imag(Er%nomega_i))
     domegaright = (omega_imag(Er%nomega_i+1)+omega_imag(Er%nomega_i))
     right(:)    = -AIMAG(omega_imag(Er%nomega_i+1)-omega_imag(Er%nomega_i))*REAL(omegame0k_tmp(:))
     left(:)     = quarter*AIMAG(domegaleft)*AIMAG(domegaright) & 
&                        +REAL(omegame0k_tmp(:))*REAL(omegame0k_tmp(:))
     weight(Er%nomega_i+1,:) = ATAN(right(:)/left(:))
     do io=2,Er%nomega_i
       domegaleft  = (omega_imag(io  )+omega_imag(io-1))
       domegaright = (omega_imag(io+1)+omega_imag(io  ))
       right(:)    = -half*AIMAG(omega_imag(io+1)-omega_imag(io-1))*REAL(omegame0k_tmp(:))
       left(:)     = REAL(omegame0k_tmp(:))*REAL(omegame0k_tmp(:)) & 
&        +quarter*AIMAG(domegaleft)*AIMAG(domegaright)
       weight(io,:) = ATAN(right(:)/left(:))
     end do

     do io = 1, Er%nomega_i+1
       if (io==1) then
         iaux = 1
       else
         iaux = Er%nomega_r+1
       endif
       num = piinv*epsm1_qbz(ig,igp,iaux)*ptwsq(ig,igp,1)
       do ios=1,nomega
         ket(ios) = ket(ios) + num*weight(io,ios)
       enddo
     end do

   enddo
 enddo

 ABI_DEALLOCATE(qpgsq)

end subroutine calc_sig_ppm_delta_clos_cd
!!***

!-----------------------------------------------------

!!****f* m_eet/gw_eet_sigma_cd
!! NAME
!! gw_eet_sigma_cd
!!
!! FUNCTION
!! Wrapper routine for the calculation of the matrix elements of Sigma 
!! using the EET and contour deformation
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine gw_eet_sigma_cd(Sigp,Sr,Dtset,Cryst,Wfs,Kmesh,Qmesh,Gsph_Max,Gsph_c,Psps,Vcp,QP_BSt,Er, &
&                       isppol,iq_bz,ik_bz,jk_bz,ik_ibz,jk_ibz,itim_q,isym_q,iq_ibz,tabr_ki, &
&                       tabr_kj,spinrot_ki,spinrot_kj,ph_mkit,ph_mkjt,nfftot_gw,ngfft_gw, &
&                       use_padfft,igfftcg0,gw_gbound,gw_mgfft,ib1,ib2,nomega_tot,nomega_sigc, &
&                       fact_sp,nspinor,sigcme_tmp,sigc,nbhomo,tim_fourdp,wtqp,wtqm, &
&                       extrapolar_distrb,can_symmetrize,epsm1_qbz,npoles_missing)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_eet_sigma_cd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: isppol,itim_q,isym_q,iq_ibz,nspinor
 integer,intent(in) :: nomega_tot,nomega_sigc
 integer,intent(in) :: wtqp,wtqm
 integer,intent(in) :: iq_bz,ik_bz,jk_bz,ik_ibz,jk_ibz,ib1,ib2
 integer,intent(inout) :: npoles_missing(ib1:ib2)
 integer,intent(out) :: nbhomo
 real(dp),intent(in) :: fact_sp
 complex(dpc),intent(in) ::  ph_mkit,ph_mkjt
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_Max
 type(gsphere_t),intent(in) :: Gsph_c
 type(Dataset_type),intent(in) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Sigma_parameters),intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr
 type(wfd_t),target,intent(inout) :: Wfs
 type(ebands_t),target,intent(in) :: QP_BSt
 type(Epsilonm1_results),intent(in) :: Er
!arrays
 integer,intent(in) :: igfftcg0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: tabr_ki(nfftot_gw),tabr_kj(nfftot_gw)
 integer,intent(in) :: extrapolar_distrb(ib1:ib2,ib1:ib2,Kmesh%nbz,Wfs%nsppol)
 real(dp),intent(in) :: spinrot_ki(4),spinrot_kj(4)
 complex(gwpc),intent(in) :: epsm1_qbz(Sigp%npwc,Sigp%npwc,Er%nomega)
 complex(dpc),intent(inout) :: sigcme_tmp(nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nsppol*Sigp%nsig_ab)
 complex(dpc),intent(inout) :: sigc(2,nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nsppol*Sigp%nsig_ab)
 logical,intent(in) :: can_symmetrize(Wfs%nsppol)

!Local variables ------------------------------
!scalars
 integer :: io,kb,jb
 integer :: niter,nptwg,iter,nbmax
 integer :: ig,igp,ib,ibv
 integer :: isym_kgw,isym_ki,iik,jik
!arrays
 integer :: fnlloc(Cryst%ntypat,2),fnlmax(Cryst%ntypat)
 real(dp),pointer :: qp_ene(:,:,:),qp_occ(:,:,:)
 real(dp),allocatable :: qplg(:,:)
 real(dp),allocatable :: kplqg(:)
 real(dp),allocatable :: omegame0k(:),omegame0lumo(:)
 real(dp),allocatable :: qbzpg(:) 
 complex(gwpc),allocatable :: vc_sqrt_qbz(:)
 complex(gwpc),allocatable :: fnlkr(:,:,:)
 complex(gwpc),allocatable :: fnlkpr(:,:,:)
 complex(gwpc),allocatable :: wfr1(:,:)
 complex(gwpc),allocatable :: mtwk(:,:)
 complex(gwpc),allocatable :: mtwkp(:,:)
 complex(gwpc),allocatable :: ptwsq(:,:,:)
 complex(dpc) :: sigctmp(nomega_sigc)

!************************************************************************

 qp_ene => QP_BSt%eig(:,:,:)
 qp_occ => QP_BSt%occ(:,:,:)

 !@Arjan:
 ! Gsph_c gamma-centered sphere for W and thus Sigma_x
 ! Gsph_Max gamma-centered sphere with gvec(3,npwvec) where Sigp%npwvec=MAX(Sigp%npwwfn,Sigp%npwx)
 ! it is used for the wavefunctions but it will be REMOVED when we switch to k-centered
 ! G-spheres for the wavefunctions

 !call gsph_init(Gsph_Max,Cryst,Sigp%npwvec,gvec=gvec_kss)

 ! min and Max band indeces for GW corrections (for this k-point)

 nbhomo=1
 do ib = 2, Sigp%nbnds 
   if (fact_sp*qp_occ(ib,jk_ibz,isppol)<GW_TOL_DOCC) exit
   nbhomo=nbhomo+1
 enddo
 nbmax=max(nbhomo,Dtset%gw_eet_nband)

 !allocate(val_idx(QP_BSt%nkpt,QP_BSt%nsppol))
 !val_idx = get_valence_idx(QP_BSt,tol3)
 !do isppol=1,nsppol
 ! nbhomo(isppol) = val_idx(1,isppol)
 ! ltest = ALL(val_idx(:,isppol))==val_idx(1,isppol),
 ! ABI_CHECK(ltest,"Optimal GW + metals not coded")
 !end do
 !deallocate(val_idx)

 niter = Dtset%gw_eet
 nptwg=1
 do iter = 1, niter
   nptwg=nptwg+mod(iter,2)
 enddo

 !call init_iskg_tabs(ISkg,Cryst,k_tim,k_sym,Wfd%ecut,istwf_k,npw_k,k_ibz,kg_k,Cryst%gmet,Wfd%ngfft,ierr)
 !ABI_CHECK(ierr==0," init_iskg_tabs returned ierr/=0")

 ABI_ALLOCATE(vc_sqrt_qbz,(Sigp%npwc))
 do ig=1,Sigp%npwc
   vc_sqrt_qbz(Gsph_c%rottb(ig,itim_q,isym_q))=Vcp%vc_sqrt(ig,iq_ibz)
 end do

 ABI_ALLOCATE(omegame0k,(nomega_tot))
 ABI_ALLOCATE(omegame0lumo,(nomega_tot))
 ABI_ALLOCATE(qplg,(3,Sigp%npwc))
 ABI_ALLOCATE(kplqg,(Sigp%npwc))
 ABI_ALLOCATE(qbzpg,(Sigp%npwc))

 isym_kgw = Kmesh%tabo(jk_bz)
 jik = (3-Kmesh%tabi(jk_bz))/2

 isym_ki = Kmesh%tabo(ik_bz)
 iik = (3-Kmesh%tabi(ik_bz))/2

 do ig=1,Sigp%npwc
   qplg(:,ig) = Qmesh%bz(:,iq_bz) + Gsph_c%gvec(:,ig)
   kplqg(ig) = -vdotw(Kmesh%bz(:,jk_bz),qplg(:,ig),Cryst%gmet,"G")
   qbzpg(ig) = normv(qplg(:,ig),Cryst%gmet,"G")
 end do

 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then

   ABI_ALLOCATE(fnlkr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlkpr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(mtwk,(Wfs%nfftot*nspinor,nbmax))
   ABI_ALLOCATE(mtwkp,(Wfs%nfftot*nspinor,ib1:ib2))

   call gw_eet_sigma_vkb(Sigp,Cryst,Wfs,Kmesh,Psps,isppol,ik_ibz,jk_ibz,ib1,ib2,nspinor, &
&                        tim_fourdp,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)
 endif

 ABI_ALLOCATE(wfr1,(Wfs%nfftot*nspinor,nbmax))

 do ibv = 1, nbmax
   call wfd_get_ur(Wfs,ibv,ik_ibz,isppol,wfr1(:,ibv))
 enddo

 do jb = ib1,ib2

   do kb = ib1,ib2

     if (Sigp%gwcalctyp/=28.and.kb/=jb) CYCLE

     if (extrapolar_distrb(jb,kb,ik_bz,isppol)/=Wfs%my_rank) CYCLE

     ABI_ALLOCATE(ptwsq,(Sigp%npwc,Sigp%npwc,niter+1))

     do io=1,Sr%nomega4sd
       omegame0k(io)  = real(Sr%omega4sd(kb,jk_ibz,io,isppol)) - qp_ene(kb,jk_ibz,isppol)
       omegame0lumo(io)= real(Sr%omega4sd(kb,jk_ibz,io,isppol)) - qp_ene(nbmax+1,ik_ibz,isppol)
     end do

     sigctmp=czero_gw

     if (niter==0.or.Dtset%gw_eet_inclvkb==0) then

       call fft4eet_sig_cd(Sigp,Cryst,Wfs,Kmesh,Gsph_c,Sr,Er,nbhomo,nbmax,nomega_tot,isppol,nfftot_gw, &
&                   ngfft_gw,use_padfft,igfftcg0,gw_gbound,gw_mgfft,iik,tabr_ki,ph_mkit,spinrot_ki, &
&                   ik_ibz,jk_ibz,isym_kgw,jik,tabr_kj,ph_mkjt,spinrot_kj,Gsph_Max, &
&                   nspinor,tim_fourdp,wfr1,vc_sqrt_qbz,Vcp%i_sz,kb,qplg,kplqg,niter, &
&                   ptwsq,ik_bz,jk_bz,npoles_missing(kb),epsm1_qbz,sigctmp)

     else

       call fft4eet_sig_kb_cd(Sigp,Cryst,Wfs,Kmesh,Gsph_c,Psps,Sr,Er,nbhomo,nbmax,nomega_tot,isppol,nfftot_gw, &
&                   ngfft_gw,use_padfft,igfftcg0,gw_gbound,gw_mgfft,iik,tabr_ki,ph_mkit,spinrot_ki, &
&                   ik_ibz,jk_ibz,isym_kgw,jik,tabr_kj,ph_mkjt,spinrot_kj,Gsph_Max, &
&                   nspinor,tim_fourdp,fnlloc,fnlmax,fnlkr,mtwk,mtwkp(:,kb),wfr1,vc_sqrt_qbz,Vcp%i_sz, &
&                   kb,qplg,kplqg,niter,ptwsq,ik_bz,jk_bz,npoles_missing(kb),epsm1_qbz,sigctmp)

     endif

     call calc_delta_ppm(Sigp,ptwsq,niter)

     do ig = 1, Sigp%npwc
       do igp = 1, Sigp%npwc
         if (ik_bz==jk_bz) then
           if (ig/=1.and.igp/=1) then
             ptwsq(ig,igp,1)= ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)*vc_sqrt_qbz(igp)
           else
             if (jb/=kb.or.kb<=nbmax.or.jb<=nbmax) then
               ptwsq(ig,igp,1)=(0.0,0.0)
             else
               if (ig==1.and.igp==1) then
                 ptwsq(ig,igp,1) = cmplx(Vcp%i_sz,0.0_gwp)
               elseif (ig==1.and.igp/=1) then
                 ptwsq(ig,igp,1) = cmplx(sqrt(Vcp%i_sz),0.0_gwp)*ptwsq(ig,igp,1)*vc_sqrt_qbz(igp)
               elseif (igp==1.and.ig/=1) then
                 ptwsq(ig,igp,1) = cmplx(sqrt(Vcp%i_sz),0.0_gwp)*ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)
               else
                 ptwsq(ig,igp,1)= ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)*vc_sqrt_qbz(igp)
               endif
             endif
           endif
         else
           ptwsq(ig,igp,1)= ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)*vc_sqrt_qbz(igp)
         endif
       enddo
     enddo

     call calc_sig_ppm_delta_clos_cd(Er,Sigp%npwc,nomega_tot,ik_bz,jk_bz,qbzpg,omegame0k,omegame0lumo, &
&                                    ptwsq,epsm1_qbz,sigctmp,Dtset%gw_eet_scale,niter)

     if (can_symmetrize(isppol)) then
       sigcme_tmp(:,jb,kb,isppol)=sigcme_tmp(:,jb,kb,isppol) + &
&         (wtqp+wtqm)*DBLE(sigctmp(:)) + (wtqp-wtqm)*j_gw*AIMAG(sigctmp(:))
       sigc(1,:,jb,kb,isppol)=sigc(1,:,jb,kb,isppol) + wtqp*      sigctmp(:)
       sigc(2,:,jb,kb,isppol)=sigc(2,:,jb,kb,isppol) + wtqm*CONJG(sigctmp(:))
     else
       sigcme_tmp(:,jb,kb,isppol)=sigcme_tmp(:,jb,kb,isppol)+sigctmp(:)
     endif

     ABI_DEALLOCATE(ptwsq)

   enddo
 enddo

 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then
   ABI_DEALLOCATE(mtwk)
   ABI_DEALLOCATE(mtwkp)
   ABI_DEALLOCATE(fnlkr)
   ABI_DEALLOCATE(fnlkpr)
 endif

 ABI_DEALLOCATE(wfr1)
 ABI_DEALLOCATE(vc_sqrt_qbz)
 ABI_DEALLOCATE(omegame0k)
 ABI_DEALLOCATE(omegame0lumo)
 ABI_DEALLOCATE(qplg)
 ABI_DEALLOCATE(kplqg)
 ABI_DEALLOCATE(qbzpg)

end subroutine gw_eet_sigma_cd
!!***
 
!----------------------------------------------------------------------

!!****f* m_eet/fft4eet_sig_cd
!! NAME                  
!! fft4eet_sig_cd
!!
!! FUNCTION
!! Calculate f^{pp}, f^{pj}, and f{jj}. See Eqs. (25),(26),(27) of PRB 85, 085126, for their expressions.
!! They are stored in ptwsq(:,:,1),ptwsq(:,:,2),ptwsq(:,:,3). 
!! Terms involving the nonlocal part of the pseudopotential are not included.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine fft4eet_sig_cd(Sigp,Cryst,Wfs,Kmesh,Gsph_c,Sr,Er,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw, &
&                  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,Gsph_max, &
&                  nspinor,tim_fourdp,wfr1,vc_sqrt_qbz,i_sz,kb,qplg,kplqg,niter, &
&                  ptwsq,ik_bz,ikmq_bz,npoles_missing,epsm1q,sigmac)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_sig_cd'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(Sigma_parameters),intent(in) :: Sigp
 type(wfd_t),target,intent(inout) :: Wfs
 type(sigma_t),intent(in) :: Sr
 type(gsphere_t),intent(in) :: Gsph_c,Gsph_max
 type(kmesh_t),intent(in) :: Kmesh
 type(Epsilonm1_results),intent(in) :: Er
 integer,intent(in) :: nbhomo,nbmax,is,kb,niter,nomega
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_kmq,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(inout) :: npoles_missing
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qplg(3,Sigp%npwc),kplqg(Sigp%npwc)
 real(dp),intent(in) :: i_sz
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt
!arrays
 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 complex(gwpc),intent(in) :: epsm1q(Sigp%npwc,Sigp%npwc,Er%nomega)
 complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
 complex(dpc),intent(inout) :: sigmac(nomega)

!Local variables-------------------------------
!scalars
 integer :: istwf_kmq,nbz,npw_kmq,ibv,ig,igp,i,j
 integer :: ig4,ig4x,ig4y,ig4z,outofbox
 integer,save :: enough=0
 complex(gwpc) :: minusone
 character(len=500) :: msg
!arrays
 integer :: gmgp(3)
 integer,pointer :: gbound_kmq(:,:),kg_kmq(:,:),igfft_kmq(:)
 integer,parameter :: ndat1=1
 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux(3)
 complex(gwpc),allocatable :: wfr2(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)
 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)

!************************************************************************

 ABI_UNUSED((/isym_kmq,tim_fourdp/))

!Dummy statement, to keep Kmesh arg
 nbz=Kmesh%nbz

 gbound_kmq => Wfs%Kdata(ikmq_ibz)%gbound
 igfft_kmq  => Wfs%Kdata(ikmq_ibz)%igfft0
 kg_kmq     => Wfs%Kdata(ikmq_ibz)%kg_k
 istwf_kmq  = Wfs%istwfk(ikmq_ibz)
 npw_kmq    = Wfs%npwarr(ikmq_ibz)

 ABI_ALLOCATE(wfr2,(Wfs%nfftot*nspinor))

 ABI_ALLOCATE(rhotwg,(Sigp%npwc*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Sigp%npwc*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(npw_kmq))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
 endif

 if (niter>1) then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
 endif

 call wfd_get_ur(Wfs,kb,ikmq_ibz,is,wfr2)
 call calc_wfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,wfr2,wfr2,wfwfg(:))

 do ibv = 1, nbmax
   call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                spinrot_kmq,nspinor,rhotwg(:,ibv))

   if (ibv>nbhomo) then
      call calc_corr_sig_cd(Sr,Sigp%npwc,Sigp%npwx,ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,ibv,kb,is,nomega, &
&                           Er%nomega,Er%nomega_r,Er%nomega_i,rhotwg(:,ibv), &
&                           Er%omega,epsm1q,npoles_missing,i_sz,vc_sqrt_qbz,sigmac)
   endif
 enddo

 if (niter>0) then

   call wfd_dur_isk(Wfs,Cryst,Kmesh,kb,ikmq_bz,is,npw_kmq,Gsph_max%rottbm1,dwfr)  !,ISkg,ik_ibz

   do ibv = 1, nbmax
     do i = 1, 3
       call drho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibv),itim_k,tabr_k,ph_mkt,dwfr(:,i), &
&                     nspinor,drhotwg(:,ibv,i))
     enddo
   enddo

   do i = 1, 3
     call calc_dwfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt,wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(itim_kmq,nfftot_gw,ngfft_gw,dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo

   !ABI_DEALLOCATE(wfr2)
   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Sigp%npwc,nbmax))

   cauxg(:,:)=(0.0,0.0)
   do ibv = 1, nbmax
     do ig = 1, Sigp%npwc
       drhaux(:)=cmplx(real(drhotwg(ig,ibv,:)),aimag(drhotwg(ig,ibv,:)))
       cauxg(ig,ibv)=vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
     enddo
   enddo

 endif
 ABI_DEALLOCATE(wfr2)

 ptwsq(:,:,:)=(0.0,0.0)
 outofbox=0
 do igp=1,Sigp%npwc
   do ig=1,Sigp%npwc
     gmgp(:)=Gsph_c%gvec(:,igp)-Gsph_c%gvec(:,ig)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)
     if (igp>=ig) then
       ptwsq(ig,igp,1)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2) + vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")
       if (niter>1.and.igp>=ig) then
         paux(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux(i)=paux(i) + vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
         enddo
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) + vdotw(qplg(:,igp),paux,Cryst%gmet,"G")
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 do ibv = 1, nbmax
   do igp=1,Sigp%npwc
     do ig=1,igp
       ptwsq(ig,igp,1)=ptwsq(ig,igp,1)-conjg(rhotwg(ig,ibv))*rhotwg(igp,ibv)
     enddo
   end do
 end do

 do igp = 1, Sigp%npwc
   do ig = igp+1, Sigp%npwc
     ptwsq(ig,igp,1)=conjg(ptwsq(igp,ig,1))
   enddo
 enddo

 if (niter>0) then
   minusone=(-1.,0.)
   do ibv = 1, nbmax
     call XGERC(Sigp%npwc,Sigp%npwc,minusone,conjg(rhotwg(:,ibv)),1,conjg(cauxg(:,ibv)),1,ptwsq(:,:,2),Sigp%npwc)
   enddo
   do igp=1,Sigp%npwc
     do ig=1,Sigp%npwc
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig
 endif

 if (niter>1) then

   do ibv = 1, nbmax
     do igp=1,Sigp%npwc
       do ig=1,igp
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3)-conjg(cauxg(ig,ibv))*cauxg(igp,ibv)
       enddo
     enddo
   enddo
   do igp=1,Sigp%npwc
     do ig=1,igp
       ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+kplqg(igp)*ptwsq(ig,igp,2)+kplqg(ig)*conjg(ptwsq(igp,ig,2))&
&                                     -kplqg(ig)*ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig

   do igp = 1, Sigp%npwc
     do ig = igp+1, Sigp%npwc
       ptwsq(ig,igp,3)=conjg(ptwsq(igp,ig,3))
     enddo
   enddo

 endif

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
 endif
 if (niter>1) then
   ABI_DEALLOCATE(ddwfwfg)
 endif

end subroutine fft4eet_sig_cd
!!***

!----------------------------------------------------------------------

!!****f* m_eet/fft4eet_sig_kb_cd
!! NAME                  
!! fft4eet_sig_kb_cd
!!
!! FUNCTION
!! Calculate f^{pp}, f^{pj}, and f{jj}. See Eqs. (25),(26),(27) of PRB 85, 085126, for their expressions.
!! They are stored in ptwsq(:,:,1),ptwsq(:,:,2),ptwsq(:,:,3). 
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine fft4eet_sig_kb_cd(Sigp,Cryst,Wfs,Kmesh,Gsph_c,Psps,Sr,Er,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw, &
&                  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,Gsph_max, &
&                  nspinor,tim_fourdp,fnlloc,fnlmax,fnlkr,mtwk,mtwkp,wfr1, &
&                  vc_sqrt_qbz,i_sz,kb,qplg,kplqg,niter,ptwsq,ik_bz,ikmq_bz,npoles_missing,epsm1q,sigmac)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_sig_kb_cd'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Sigma_parameters),intent(in) :: Sigp
 type(wfd_t),target,intent(inout) :: Wfs
 !type(kb_potential) :: KBff_kmq_ibz
 type(sigma_t),intent(in) :: Sr
 type(gsphere_t),intent(in) :: Gsph_c, Gsph_max
 type(Epsilonm1_results),intent(in) :: Er

 integer,intent(in) :: nbhomo,nbmax,is,kb,niter,nomega
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_kmq,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)
 integer,intent(inout) :: npoles_missing
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qplg(3,Sigp%npwc),kplqg(Sigp%npwc)
 real(dp),intent(in) :: i_sz
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt

 complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(in) :: mtwk(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: mtwkp(Wfs%nfftot*nspinor)

 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)

 complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 complex(gwpc),intent(in) :: epsm1q(Sigp%npwc,Sigp%npwc,Er%nomega)

 complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
 complex(dpc),intent(inout) :: sigmac(nomega)

!Local variables-------------------------------
!scalars
 integer :: istwf_kmq,npw_kmq,i,j,ibv,ilm,iat,ig,igp 
 integer :: ig4,ig4x,ig4y,ig4z
 integer :: ig5,ig5x,ig5y,ig5z,nlx,outofbox
 integer,save :: enough=0
 complex(gwpc) :: minusone
 character(len=500) :: msg
!arrays
 integer,pointer :: gbound_kmq(:,:),kg_kmq(:,:),igfft_kmq(:)
 integer :: gmgp(3)
 integer,parameter :: ndat1=1
 complex(gwpc),allocatable :: wfr2(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg2(:,:)
 complex(gwpc),allocatable :: fnltwg3(:,:)
 complex(gwpc),allocatable :: kns(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)
 complex(gwpc),allocatable :: fnlwfg(:)
 complex(gwpc),allocatable :: fkdwfg(:,:)
 complex(gwpc),allocatable :: fdrhotwg(:,:,:,:)
 complex(gwpc),allocatable :: lnkp(:)
 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)
 complex(gwpc),allocatable :: ff(:,:,:,:)
 complex(gwpc),allocatable :: vzn(:,:,:)
 complex(gwpc),allocatable :: paux(:,:)
 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux2(3)

!************************************************************************

 ABI_UNUSED((/isym_kmq,tim_fourdp/))

 gbound_kmq => Wfs%Kdata(ikmq_ibz)%gbound
 igfft_kmq  => Wfs%Kdata(ikmq_ibz)%igfft0
 kg_kmq     => Wfs%Kdata(ikmq_ibz)%kg_k
 istwf_kmq  = Wfs%istwfk(ikmq_ibz)
 npw_kmq    = Wfs%npwarr(ikmq_ibz)

 nlx = min(Psps%mpsang,4)

 ABI_ALLOCATE(wfr2,(Wfs%nfftot*nspinor))
 ABI_ALLOCATE(rhotwg,(Sigp%npwc*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Sigp%npwc*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(npw_kmq))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
   ABI_ALLOCATE(fnltwg,(Sigp%npwc*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlwfg,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(fnltwg2,(Sigp%npwc*nspinor**2,nbmax))
   ABI_ALLOCATE(fnltwg3,(Sigp%npwc*nspinor**2,nbmax))
 endif

 if (niter>1) then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
   ABI_ALLOCATE(lnkp,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(kns,(Sigp%npwc*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fdrhotwg,(Sigp%npwc*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom,3))
   ABI_ALLOCATE(fkdwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(ff,(nlx*nlx,Cryst%natom,nlx*nlx,Cryst%natom))
   ABI_ALLOCATE(vzn,(Sigp%npwc*nspinor**2,nlx*nlx,Cryst%natom))
 endif

 call wfd_get_ur(Wfs,kb,ikmq_ibz,is,wfr2)
 call calc_wfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,wfr2,wfr2,wfwfg)

 do ibv = 1, nbmax
   call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                spinrot_kmq,nspinor,rhotwg(:,ibv))

   if (ibv>nbhomo) then
      call calc_corr_sig_cd(Sr,Sigp%npwc,Sigp%npwx,ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,ibv,kb,is,nomega, &
&                           Er%nomega,Er%nomega_r,Er%nomega_i,rhotwg(:,ibv), &
&                           Er%omega,epsm1q,npoles_missing,i_sz,vc_sqrt_qbz,sigmac)
   endif
 enddo

 if (niter>0) then

   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                    fnlkr(:,ilm,iat),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                    spinrot_kmq,nspinor,fnltwg(:,ilm,iat))
     enddo
   enddo

   call calc_wfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,wfr2,mtwkp,fnlwfg)

   call wfd_dur_isk(Wfs,Cryst,Kmesh,kb,ikmq_bz,is,npw_kmq,Gsph_max%rottbm1,dwfr)  !,ISkg,ik_ibz

   do ibv = 1, nbmax
     do i = 1, 3
       call drho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibv),itim_k,tabr_k,ph_mkt,dwfr(:,i), &
&                     nspinor,drhotwg(:,ibv,i))
     enddo
     call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  wfr1(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,mtwkp,itim_kmq,tabr_kmq,ph_mkmqt, &
&                  spinrot_kmq,nspinor,fnltwg2(:,ibv))

     call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  mtwk(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                  spinrot_kmq,nspinor,fnltwg3(:,ibv))
   enddo

   do i = 1, 3
     call calc_dwfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt,wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(itim_kmq,nfftot_gw,ngfft_gw,dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo

   if (niter>1) then

     do iat = 1, Cryst%natom
       do ilm = 1, Psps%mpsang*Psps%mpsang
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         do i = 1, 3
           call drho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                         fnlkr(:,ilm,iat),itim_k,tabr_k,ph_mkt,dwfr(:,i), &
&                         nspinor,fdrhotwg(:,ilm,iat,i))
         enddo
         call rho_tw_g(nspinor,Sigp%npwc,nfftot_gw,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                      fnlkr(:,ilm,iat),itim_k,tabr_k,ph_mkt,spinrot_k,mtwkp,itim_kmq,tabr_kmq,ph_mkmqt, &
&                      spinrot_kmq,nspinor,kns(:,ilm,iat))
       enddo
     enddo

     do i = 1, 3
       call calc_dwfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt,mtwkp,dwfr(:,i),fkdwfg(:,i))
     enddo
     call calc_wfwfg(tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,mtwkp,mtwkp,lnkp)

   endif

   ABI_DEALLOCATE(wfr2)
   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Sigp%npwc,nbmax))

   cauxg(:,:)=(0.0,0.0)
   do ibv = 1, nbmax
     do ig=1,Sigp%npwc
       drhaux(:)=cmplx(real(drhotwg(ig,ibv,:)),aimag(drhotwg(ig,ibv,:)))
       cauxg(ig,ibv)=vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
     enddo
   enddo
   do ig=1,Sigp%npwc
     cauxg(ig,:)=cauxg(ig,:)-fnltwg2(ig,:)+fnltwg3(ig,:)
   enddo

 endif

 ptwsq(:,:,:)=(0.0,0.0)
 outofbox=0
 do igp=1,Sigp%npwc
   do ig=1,Sigp%npwc
     gmgp(:)=Gsph_c%gvec(:,igp)-Gsph_c%gvec(:,ig)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)

     ig5x= modulo(-gmgp(1),ngfft_gw(1))
     ig5y= modulo(-gmgp(2),ngfft_gw(2))
     ig5z= modulo(-gmgp(3),ngfft_gw(3))
     ig5= 1+ig5x+ig5y*ngfft_gw(1)+ig5z*ngfft_gw(1)*ngfft_gw(2)

     if (igp>=ig) then
       ptwsq(ig,igp,1)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2) + vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")-fnlwfg(ig4)
       if (niter>1.and.igp>=ig) then
         paux2(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux2(i)=paux2(i) + vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")
         enddo
         drhaux(:)=paux2(:)-cmplx(real(fkdwfg(ig4,:)),aimag(fkdwfg(ig4,:)))
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) + vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")
         drhaux(:)=cmplx(real(fkdwfg(ig5,:)),-aimag(fkdwfg(ig5,:)))
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) - vdotw(qplg(:,ig),drhaux,Cryst%gmet,"G")+lnkp(ig4)
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 do ibv = 1, nbmax
   do igp=1,Sigp%npwc
     do ig=1,igp
       ptwsq(ig,igp,1)=ptwsq(ig,igp,1)-conjg(rhotwg(ig,ibv))*rhotwg(igp,ibv)
     enddo
   end do !igp
 end do !ig

 do igp = 1, Sigp%npwc
   do ig = igp+1, Sigp%npwc
     ptwsq(ig,igp,1)=conjg(ptwsq(igp,ig,1))
   enddo
 enddo

 if (niter>0) then
   ABI_ALLOCATE(paux,(Sigp%npwc,Sigp%npwc))
   paux(:,:)=(0.0,0.0)
   minusone=(-1.,0.)
   do ibv = 1, nbmax
     call XGERC(Sigp%npwc,Sigp%npwc,minusone,conjg(rhotwg(:,ibv)),1,conjg(cauxg(:,ibv)),1,ptwsq(:,:,2),Sigp%npwc)
   enddo
   do igp=1,Sigp%npwc
     do ig=1,Sigp%npwc
       if (ig>=igp) then
         do iat = 1, Cryst%natom
           do ilm = 1, nlx*nlx
             if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
             if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
             paux(ig,igp)=paux(ig,igp)+conjg(fnltwg(ig,ilm,iat))*fnltwg(igp,ilm,iat)
           enddo
         enddo
       endif
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig
   do ig=1,Sigp%npwc
     do igp=1,Sigp%npwc
       if (ig>=igp) then
         ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+paux(ig,igp)
       else
         ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+conjg(paux(igp,ig))
       endif
     end do !igp
   end do !ig
   ABI_DEALLOCATE(paux)
 endif

 if (niter>1) then

   call compute_ff(Cryst,Kmesh,Psps,ikmq_bz,nlx,istwf_kmq,npw_kmq,kg_kmq,Gsph_max%rottbm1,fnlmax,fnlloc,ff)

   vzn(:,:,:)=(0.0,0.0)
   do igp=1,Sigp%npwc
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         vzn(igp,:,:) = vzn(igp,:,:) + half*ff(ilm,iat,:,:)*fnltwg(igp,ilm,iat)
       enddo
     enddo
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         drhaux(:)=cmplx(real(fdrhotwg(igp,ilm,iat,:)),aimag(fdrhotwg(igp,ilm,iat,:)))
         vzn(igp,ilm,iat)=vzn(igp,ilm,iat)+vdotw(qplg(:,igp),drhaux,Cryst%gmet,"G")-kns(igp,ilm,iat)
       enddo
     enddo
     do ig=1,igp
       do ibv = 1, nbmax
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3)-conjg(cauxg(ig,ibv))*cauxg(igp,ibv)
       enddo
       ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+kplqg(igp)*ptwsq(ig,igp,2)+kplqg(ig)*conjg(ptwsq(igp,ig,2))&
&                                     -kplqg(ig)*ptwsq(ig,igp,1)*kplqg(igp)
       do iat = 1, Cryst%natom
         do ilm = 1, nlx*nlx
           if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
           if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
           ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+conjg(fnltwg(ig,ilm,iat))*vzn(igp,ilm,iat)+ &
&                                                fnltwg(igp,ilm,iat)*conjg(vzn(ig,ilm,iat))
         enddo
       enddo
     end do !igp
   end do !ig

   do igp = 1, Sigp%npwc
     do ig = igp+1, Sigp%npwc
       ptwsq(ig,igp,3)=conjg(ptwsq(igp,ig,3))
     enddo
   enddo

 endif

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
   ABI_DEALLOCATE(fnltwg)
   ABI_DEALLOCATE(fnltwg2)
   ABI_DEALLOCATE(fnltwg3)
   ABI_DEALLOCATE(fnlwfg)
 endif
 if (niter>1) then
   ABI_DEALLOCATE(ddwfwfg)
   ABI_DEALLOCATE(lnkp)
   ABI_DEALLOCATE(kns)
   ABI_DEALLOCATE(fdrhotwg)
   ABI_DEALLOCATE(fkdwfg)
   ABI_DEALLOCATE(vzn)
   ABI_DEALLOCATE(ff)
 endif

end subroutine fft4eet_sig_kb_cd
!!***

!----------------------------------------------------------------------

!!****f* m_eet/calc_corr_sig_cd
!! NAME
!! calc_corr_sig_cd
!!
!! FUNCTION
!! Calculate sum-over-states correction to the matrix element of the GW self-energy obtained with the EET
!! using the contour deformation method
!!
!! INPUTS
!!  nomega=Total number of frequencies where $\Sigma_c$ matrix elements are evaluated.
!!  nomegae=Number of frequencies where $\epsilon^{-1}$ has been evaluated.
!!  nomegaei=Number of imaginary frequencies for $\epsilon^{-1}$ (non zero).
!!  nomegaer=Number of real frequencies for $\epsilon^{-1}$
!!  npwc=Number of G vectors for the correlation part.
!!  npwx=Number of G vectors in rhotwgp for each spinorial component.
!!  nspinor=Number of spinorial components.
!!  theta_mu_minus_e0i=1 if e0i is occupied, 0 otherwise. Fractional occupancy in case of metals. 
!!  omegame0i(nomega)=Contains $\omega-\epsilon_{k-q,b1,\sigma}$
!!  epsm1q(npwc,npwc,nomegae)=Symmetrized inverse dielectric matrix (exchange part is subtracted).
!!  omega(nomegae)=Set of frequencies for $\epsilon^{-1}$.
!!  rhotwgp(npwx*nspinor)=Matrix elements: $<k-q,b1,\sigma|e^{-i(q+G)r} |k,b2,\sigma>*vc_sqrt$
!!
!! OUTPUT
!! ket(npwc,nomega)=Contains \Sigma_c(\omega)|\phi> in reciprocal space. 
!!
!! SIDE EFFECTS
!! npoles_missing=Incremented with the number of poles whose contribution has not been taken into account due to
!!  limited frequency mesh used for W.
!!
!! PARENTS
!!      m_eet
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_corr_sig_cd(Sr,npwc,npwx,ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,ibv,kb,is,nomega,nomegae, &
&                           nomegaer,nomegaei,rhotwg,omega,epsm1q,npoles_missing,i_sz,vc_sqrt_qbz,sigmac)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_corr_sig_cd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(sigma_t),intent(in) :: Sr

 integer,intent(in) :: nomega,nomegae,nomegaei,nomegaer,npwc,npwx
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,ibv,kb,is
 integer,intent(inout) :: npoles_missing
!arrays
 real(dp),intent(in) :: i_sz
 complex(dpc),intent(in) :: omega(nomegae)
 complex(gwpc) :: epsm1q(npwc,npwc,nomegae) 
 complex(gwpc),intent(in) :: rhotwg(npwx)
 complex(gwpc),intent(in) :: vc_sqrt_qbz(npwc)
 complex(dpc),intent(inout) :: sigmac(nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,io,ios,my_err,ierr
 real(dp) :: rt_imag,rt_real,local_one,local_zero
 complex(dpc) :: ct,domegaleft,domegaright
 complex(gwpc) :: fact
!arrays
 real(dp) :: omegame0i(nomega),omegame0i_tmp(nomega)
 real(dp) :: tmp_x(2),tmp_y(2)
 real(dp) :: left(nomega),right(nomega)
 real(dp) :: rtmp_r(nomegaer),rtmp_i(nomegaer)
 complex(dpc) :: omega_imag(nomegaei+1)
 complex(gwpc) :: epsrho(npwc,nomegae),epsrho_imag(npwc,nomegaei+1)
 complex(gwpc) :: weight(nomegaei+1,nomega)
 logical :: my_calc_poles(nomega)
!
 complex(gwpc),allocatable :: rtaux(:)
 complex(gwpc),allocatable :: ket(:,:)
!*************************************************************************

 ABI_ALLOCATE(rtaux,(npwc))
 ABI_ALLOCATE(ket,(npwc,nomega))

 my_calc_poles=.TRUE.
 my_err=0

 ! Avoid divergences in $\omega - \omega_s$.
 omegame0i(:) = real(Sr%omega4sd(kb,ikmq_ibz,:,is)) - Sr%e0(ibv,ik_ibz,is)
 omegame0i_tmp(:)=omegame0i(:)
 do ios=1,nomega
   if (ABS(omegame0i_tmp(ios))<tol6) omegame0i_tmp(ios)=sign(tol6,omegame0i_tmp(ios))
 end do

 do ig = 1,npwc
   rtaux(ig)=rhotwg(ig)*vc_sqrt_qbz(ig)
 enddo
 if (ik_bz==ikmq_bz) then
   rtaux(1)=czero_gw
   if (ibv==kb) then
     rtaux(1)=cmplx(sqrt(i_sz),0.0_gwp)
   endif
 endif

 !
 ! Calculate $ \sum_{Gp} (\epsilon^{-1}_{G Gp}(\omega)-\delta_{G Gp}) \rhotwgp(Gp) $
!$omp parallel do
 do io=1,nomegae
   call XGEMV('N',npwc,npwc,cone_gw,epsm1q(:,:,io),npwc,rtaux(:),1,czero_gw,epsrho(:,io),1)
 end do

 ! Integrand along the imaginary axis.
 epsrho_imag(:,1)=epsrho(:,1)
 epsrho_imag(:,2:nomegaei+1)=epsrho(:,nomegaer+1:nomegae)

 ! Frequency mesh for integral along the imaginary axis.    
 omega_imag(1)=omega(1)
 omega_imag(2:nomegaei+1)=omega(nomegaer+1:nomegae)

 weight(1,:) = ATAN(-half*AIMAG(omega_imag(2))/REAL(omegame0i_tmp(:)))
 domegaleft  = (three*omega_imag(nomegaei+1)-omega_imag(nomegaei))
 domegaright = (omega_imag(nomegaei+1)+omega_imag(nomegaei))
 right(:)    = -AIMAG(omega_imag(nomegaei+1)-omega_imag(nomegaei))*REAL(omegame0i_tmp(:))
 left(:)     = quarter*AIMAG(domegaleft)*AIMAG(domegaright) &
&                +REAL(omegame0i_tmp(:))*REAL(omegame0i_tmp(:))
 weight(nomegaei+1,:) = ATAN(right(:)/left(:))
 ! Calculate the rest of the weights
 do io=2,nomegaei
   domegaleft  = (omega_imag(io  )+omega_imag(io-1))
   domegaright = (omega_imag(io+1)+omega_imag(io  ))
   right(:)    = -half*AIMAG(omega_imag(io+1)-omega_imag(io-1))*REAL(omegame0i_tmp(:))
   left(:)     = REAL(omegame0i_tmp(:))*REAL(omegame0i_tmp(:)) &
    +quarter*AIMAG(domegaleft)*AIMAG(domegaright)
   weight(io,:) = ATAN(right(:)/left(:))
 end do

 ! Use BLAS call to perform matrix-matrix multiplication and accumulation
 fact = CMPLX(piinv,zero) 

 ket = 0.0_gwp

 call xgemm('N','N',npwc,nomega,nomegaei+1,fact,epsrho_imag,npwc,&
&  weight,nomegaei+1,cone_gw,ket(1:npwc,:),npwc)

 local_one = one
 local_zero = zero

 ! ============================================
 ! ==== Add contribution coming from poles ====
 ! ============================================
 ! First see if the contribution has been checked before the routine is entered
 do ios=1,nomega
   if (omegame0i_tmp(ios)>tol12) then
     if (local_one<tol12) my_calc_poles(ios) = .FALSE.
   end if
   if (omegame0i_tmp(ios)<-tol12) my_calc_poles(ios) = .FALSE.
 end do !ios

 if (ANY(my_calc_poles(:))) then ! Make sure we only enter if necessary

! *** OPENMP SECTION *** Added by MS
!!OMP !write(std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',xomp_get_num_threads()
!$OMP PARALLEL SHARED(npwc,nomega,nomegaer,local_one,local_zero, &
!$OMP                    omega,epsrho,omegame0i_tmp,ket,my_calc_poles) &
!$OMP PRIVATE(ig,ios,rtmp_r,rtmp_i,tmp_x,tmp_y,rt_real,rt_imag,ct,ierr) REDUCTION(+:my_err) 
!!OMP $ write(std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',xomp_get_num_threads()
!$OMP DO 
   do ig=1,npwc
     !
     ! * Prepare the spline interpolation by filling at once the arrays rtmp_r, rtmp_i
     call spline(DBLE(omega(1:nomegaer)),DBLE(epsrho(ig,1:nomegaer)),nomegaer,local_zero,local_zero,rtmp_r)
     call spline(DBLE(omega(1:nomegaer)),DBLE(AIMAG(epsrho(ig,1:nomegaer))),nomegaer,local_zero,local_zero,rtmp_i)

     !! call spline_complex( DBLE(omega(1:nomegaer)), epsrho(ig,1:nomegaer), nomegaer, zero, zero, rtmp )

     do ios=1,nomega

       if (.NOT.my_calc_poles(ios)) CYCLE

       ! * Interpolate real and imaginary part of epsrho at |omegame0i_tmp|.
       tmp_x(1) = ABS(omegame0i_tmp(ios))
       call splint(nomegaer,DBLE(omega(1:nomegaer)),DBLE(epsrho(ig,1:nomegaer)),rtmp_r,1,tmp_x,tmp_y,ierr=ierr)
       if (ig==1) my_err = my_err + ierr
       rt_real = tmp_y(1)

       tmp_x(1) = ABS(omegame0i_tmp(ios))
       call splint(nomegaer,DBLE(omega(1:nomegaer)),DBLE(AIMAG(epsrho(ig,1:nomegaer))),rtmp_i,1,tmp_x,tmp_y)
       rt_imag = tmp_y(1)

       !!call splint_complex(nomegaer,DBLE(omega(1:nomegaer)),epsrho(ig,1:nomegaer),rtmp,1,tmp_x,yfit)

       ct=DCMPLX(rt_real,rt_imag)

       ket(ig,ios)=ket(ig,ios)+ct*local_one

     end do !ios
   end do !ig
!$OMP END DO
!$OMP END PARALLEL
 end if ! ANY(my_calc_poles)

 do ios=1,nomega
   sigmac(ios)=sigmac(ios) + XDOTC(npwc,rtaux(:),1,ket(:,ios),1)
 end do

 npoles_missing = npoles_missing + my_err

 ABI_DEALLOCATE(rtaux)
 ABI_DEALLOCATE(ket)

end subroutine calc_corr_sig_cd
!!***

!----------------------------------------------------------------------

END MODULE m_eet
!!***
