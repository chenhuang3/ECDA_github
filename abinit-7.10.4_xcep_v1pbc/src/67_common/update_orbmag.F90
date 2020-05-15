!{\src2tex{textfont=tt}}
!!****f* ABINIT/update_orbmag
!! NAME
!! update_orbmag
!!
!! FUNCTION
!! This routine updates the orbital magnetization
!!
!! COPYRIGHT
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cprj(natom,mcprj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!  dimffnl = 2nd dimension of ffnl (1 + number of derivatives)
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat) = nonlocal form factors
!!  filstat=name of the status file
!!  gs_hamk <type(gs_hamiltonian_type> = current ground state hamiltonian at this k point
!!  ikpt=index of kpt currently treated
!!  isppol=index of spin polarization currently treated
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kinpw(npw_k)=(modified) kinetic energy for each plane wave (Hartree)
!!  mband=maximum number of bands
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkmem=number of k points treated by this node.
!!  mpw=maximum dimensioned size of npw
!!  natom=number of atoms in cell
!!  nkpt=number of k-points
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  npw_k = number of plane waves for this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  paral_kgb = flag defining parallelism in getghc
!!  ph3d(2,npw_k,matblk) = 3D structure factors for each atom and planewave
!!  prtvol = flag controlling verbosity of output
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtbfield <type(bfield_type)> = variables related to orbital magnetization
!!  mpi_enreg=informations about MPI parallelization
!!
!! NOTES
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      getghc,pawcprj_alloc,pawcprj_copy,pawcprj_destroy,pawcprj_get,smatrix
!!      smatrix_k_paw_bfield,store_bfield_cprj,sym_pawcprj_kn
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

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

subroutine update_orbmag(cg,cgq,cprj,dimffnl,dtbfield,ffnl,filstat,gs_hamk,&
&            icg,ikpt,isppol,kg,kinpw,mband,mcg,mcgq,mcprj,mkgq,mkmem,&
&            mpi_enreg,mpw,natom,nkpt,npw_k,npwarr,nsppol,ntypat,paral_kgb,&
&            pawtab,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,pwnsfacq,vlocal)


 use defs_basis
 use defs_abitypes
 use m_bfield
 use m_profiling_abi
 use m_xmpi

 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawtab,      only : pawtab_type
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc,  pawcprj_copy, pawcprj_get, pawcprj_destroy
 use m_fock,        only : fock_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'update_orbmag'
 use interfaces_66_paw
 use interfaces_66_wfs
 use interfaces_67_common, except_this_one => update_orbmag
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: dimffnl,icg,ikpt,isppol,mband,mcg,mcgq,mcprj,mkgq,mkmem,mpw
 integer, intent(in) :: natom,nkpt,npw_k,nsppol,ntypat,paral_kgb,prtvol,pwind_alloc
 character(len=*),intent(in) :: filstat
 type(gs_hamiltonian_type),intent(in) :: gs_hamk
 type(MPI_type), intent(inout) :: mpi_enreg
 type(bfield_type), intent(inout) :: dtbfield
!arrays
 real(dp),intent(in) :: ffnl(npw_k,dimffnl,gs_hamk%lmnmax,ntypat)
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp),intent(in) :: cg(2,mcg),cgq(2,mcgq)
 real(dp),intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq)
 real(dp),intent(inout) :: kinpw(npw_k),ph3d(2,npw_k,gs_hamk%matblk)
 real(dp), intent(inout) :: vlocal(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc) ! this variable is inout in getghc
 type(pawcprj_type),intent(in) :: cprj(natom,mcprj)

!Local variables-------------------------------
!scalars
 integer :: bband,bbs,bdir,bfor,blmn,bsig,cpopt,ddkflag,iat,iatom,icg1,icp,icp2
 integer :: idij,idir,idum1,ifor,ikg,ikgf,ikptf,ikpt2,ikpt2f,indhk
 integer :: istep,itrs,itypat,job,kband,kbs,kdir,kfor,klmn,ksig
 integer :: nband_k,ndat,nn,nnp,nnpp,nnppp,npw_k2,nspinor,mcg1_k,mcg_q
 integer :: shiftbd,sij_opt,tim_getghc,type_calc
 real(dp) :: lambda,sfac
 complex(dpc) :: cpb,cpk,dterm,IA,IA1,IA2,IA3,IB,IB1,IB2,IB3,IB4
 complex(dpc) :: IIA,IIA2,IIA3,IIA4,IIIA,IIIA1,IIIA2,IIIA3,IIIA4,twdij
 type(fock_type),pointer :: fock => null()
!arrays
 integer,allocatable :: dimlmn(:),dimlmn_srt(:),kg_k(:,:),pwind_k(:),sflag_k(:)
 real(dp) :: dk(3),dotri(2),dtm_k(2)
 real(dp),allocatable :: bwave(:,:),cg1_k(:,:),cgq_k(:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
 real(dp),allocatable :: kwave(:,:),pwnsfac_k(:,:),smat_inv(:,:,:),smat_k(:,:,:),smat_k_paw(:,:,:)
 real(dp),allocatable :: tcg(:,:,:,:)
 complex(dpc),allocatable :: omat(:,:,:,:,:,:)
 type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:),cprj_k(:,:),cprj_kb(:,:)
 type(pawcprj_type),allocatable :: kcprj(:,:),tcprj(:,:,:,:)


!!debug scalars
! integer :: bra_start,bra_end,ipw,ket_start,ket_end
! integer :: ilmn,ispinor,jlmn,spnshft,spnipw
! real(dp) :: avg_err_ovlp,err_ovlp,mag_ovlp,max_err_ovlp,ovlp_i,ovlp_r,paw_i,paw_r,tot_r,tot_i
! complex(dpc) :: cterm
!!debug arrays
! real(dp),allocatable :: bra(:,:),ket(:,:)
 type(pawtab_type),intent(in) :: pawtab(ntypat) ! has to be passed in as argument

! *********************************************************************

!DEBUG
!write(std_out,'(a)')' entering update_orbmag '
!END DEBUG

!XG130427
!This is to fool the check on unused dummy arguments. 
!pawtab is indeed used only in case of debugging.
!Should be removed when debugging phase is finished ...
 if(.false.)then
   write(std_out,*)pawtab(1)%sij(1)
 end if

!==============================================
!store new cprj into bfield structure
!==============================================

 call store_bfield_cprj(gs_hamk%atindx1,cprj,dtbfield,ikpt,isppol,mband,&
& mcprj,mkmem,mpi_enreg,natom,nkpt,nsppol)

!=============================================================================
!update smatrix for neighboring k points. Note that S^{-1} is not used in this
!implementation, unlike what is done in electric field case.
!=============================================================================

 nband_k = dtbfield%nband_occ
 nspinor = dtbfield%nspinor
 mcg1_k = mpw*nband_k*nspinor

 ABI_ALLOCATE(tcg,(2,mcg1_k,3,2))
 ABI_ALLOCATE(cg1_k,(2,mcg1_k))
 tcg(:,:,:,:) = zero

 mcg_q = mpw*mband*nspinor
 ABI_ALLOCATE(cgq_k,(2,mcg_q))
 ABI_ALLOCATE(sflag_k,(nband_k))
 ABI_ALLOCATE(pwind_k,(mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,mpw))

 ABI_DATATYPE_ALLOCATE(cprj_k,(natom,nband_k*nspinor))
 ABI_DATATYPE_ALLOCATE(cprj_kb,(natom,nband_k*nspinor))
 ABI_ALLOCATE(dimlmn_srt,(natom))
 iatom = 0
 do itypat = 1, gs_hamk%ntypat
   do iat = 1, gs_hamk%nattyp(itypat)
     iatom = iatom + 1
     dimlmn_srt(iatom)=dtbfield%lmn_size(itypat)
   end do
 end do
 ABI_ALLOCATE(dimlmn,(natom))
 do iatom = 1, natom
   itypat = gs_hamk%typat(iatom)
   dimlmn(iatom)=dtbfield%lmn_size(itypat)
 end do
 call pawcprj_alloc(cprj_k,0,dimlmn)
 call pawcprj_alloc(cprj_kb,0,dimlmn)
 if (nkpt /= dtbfield%fnkpt) then
   ABI_DATATYPE_ALLOCATE(cprj_fkn,(natom,nband_k*nspinor))
   ABI_DATATYPE_ALLOCATE(cprj_ikn,(natom,nband_k*nspinor))
   call pawcprj_alloc(cprj_fkn,0,dimlmn)
   call pawcprj_alloc(cprj_ikn,0,dimlmn)
 end if

 ABI_DATATYPE_ALLOCATE(tcprj,(natom,nband_k*nspinor,3,2))

 do idir = 1, 3
   do ifor = 1, 2
     call pawcprj_alloc(tcprj(:,:,idir,ifor),0,dimlmn)
   end do
 end do

 ikptf = dtbfield%i2fbz(ikpt)
 ikgf = dtbfield%fkgindex(ikptf)  ! this is the shift for pwind

 icp=dtbfield%cprjindex(ikpt,isppol)
 call pawcprj_get(gs_hamk%atindx1,cprj_k,dtbfield%cprj,natom,1,icp,ikpt,0,isppol,&
& nband_k,nkpt,natom,nband_k,nband_k,nspinor,nsppol,0,&
& mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

!!=======================================
!! code to test orthonormality of cg_k
!!=======================================
!
!ABI_ALLOCATE(bra,(2,npw_k*nspinor))
!ABI_ALLOCATE(ket,(2,npw_k*nspinor))
!max_err_ovlp=0.0
!do bband = 1, nband_k
!bra_start = dtbfield%cgindex(ikpt,nsppol)+1+(bband-1)*npw_k*nspinor
!bra_end = bra_start + npw_k*nspinor - 1
!bra(1:2,1:npw_k*nspinor) = cg(1:2,bra_start:bra_end)
!do kband = 1, nband_k
!ket_start = dtbfield%cgindex(ikpt,nsppol)+1+(kband-1)*npw_k*nspinor
!ket_end = ket_start + npw_k*nspinor - 1
!ket(1:2,1:npw_k*nspinor) = cg(1:2,ket_start:ket_end)
!
!tot_r = 0.0; tot_i = 0.0     
!do ispinor = 1, nspinor
!ovlp_r = 0.0; ovlp_i = 0.0
!spnshft = (ispinor-1)*npw_k
!do ipw = 1, npw_k
!spnipw = ipw + spnshft
!ovlp_r = ovlp_r + bra(1,spnipw)*ket(1,spnipw)+bra(2,spnipw)*ket(2,spnipw)
!ovlp_i = ovlp_i - bra(2,spnipw)*ket(1,spnipw)+bra(1,spnipw)*ket(2,spnipw)
!end do ! end loop over ipw
!paw_r = 0.0
!paw_i = 0.0
!do iatom = 1, natom
!itypat = gs_hamk%typat(iatom)
!do ilmn = 1, dtbfield%lmn_size(itypat)
!do jlmn = 1, dtbfield%lmn_size(itypat)
!klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
!bbs = nspinor*(bband-1)+ispinor
!kbs = nspinor*(kband-1)+ispinor
!cpb=cmplx(cprj_k(iatom,bbs)%cp(1,ilmn),cprj_k(iatom,bbs)%cp(2,ilmn))
!cpk=cmplx(cprj_k(iatom,kbs)%cp(1,jlmn),cprj_k(iatom,kbs)%cp(2,jlmn))
!cterm = conjg(cpb)*pawtab(itypat)%sij(klmn)*cpk
!paw_r = paw_r + real(cterm)
!paw_i = paw_i + aimag(cterm)
!end do ! end loop over jlmn
!end do ! end loop over ilmn
!end do ! end loop over iatom
!tot_r = tot_r + ovlp_r + paw_r
!tot_i = tot_i + ovlp_i + paw_i
!end do ! end loop over ispinor
!!     write(std_out,'(a,3i4,2es16.8)')' JWZ Debug: update_orbmag ikpt bband kband ovlp : ',&
!!&           ikpt,bband,kband,tot_r,tot_i
!mag_ovlp = tot_r*tot_r + tot_i*tot_i
!if(bband==kband) then
!err_ovlp=abs(mag_ovlp-1.0)
!else
!err_ovlp=abs(mag_ovlp)
!end if
!max_err_ovlp=MAX(max_err_ovlp,err_ovlp)
!end do ! end loop over kband
!end do ! end loop over bband
!write(std_out,'(a,i4,es16.8)')' JWZ Debug: update_orbmag ikpt ovlp err : ',&
!&           ikpt,max_err_ovlp
!ABI_DEALLOCATE(bra)
!ABI_DEALLOCATE(ket)
!
!!=========================================
!! end code to test orthonormality of cg_k
!!=========================================


 ABI_ALLOCATE(smat_k,(2,nband_k,nband_k))
 ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
 ABI_ALLOCATE(smat_k_paw,(2,nband_k,nband_k))

 job = 20 ! update overlap matrix and transfer cgq to cg1_k without S^{-1}
 shiftbd = 1 ! cg contains all bands
 ddkflag = 0 ! do not multiply wavefunctions at k by S^{-1}

 do idir = 1, 3

   dk(:) = dtbfield%dkvecs(:,idir)

   do ifor = 1, 2

     ikpt2f = dtbfield%ikpt_dk(ikptf,ifor,idir)
     if (dtbfield%indkk_f2ibz(ikpt2f,6) == 1) then
       itrs = 10
     else
       itrs = 0
     end if
     ikpt2 = dtbfield%indkk_f2ibz(ikpt2f,1)
     npw_k2 = npwarr(ikpt2)
     pwind_k(1:npw_k) = pwind(ikgf+1:ikgf+npw_k,ifor,idir)
     pwnsfac_k(1:2,1:npw_k) = pwnsfac(1:2,ikgf+1:ikgf+npw_k)
     sflag_k(:) = dtbfield%sflag(:,ikpt,ifor,idir)
     smat_k(:,:,:) = dtbfield%smat(:,:,:,ikpt,ifor,idir)

     if (xmpi_paral== 1) then
       icg1 = dtbfield%cgqindex(2,ifor+2*(idir-1),ikpt) ! nsppol implicitly = 1
       cgq_k(:,1:nband_k*nspinor*npw_k2) = &
&       cgq(:,icg1+1:icg1+nband_k*nspinor*npw_k2)
       idum1 = dtbfield%cgqindex(3,ifor+2*(idir-1),ikpt) ! nsppol implicitly = 1
       pwnsfac_k(3:4,1:npw_k2) = pwnsfacq(1:2,idum1+1:idum1+npw_k2)
     else
       icg1 = dtbfield%cgindex(ikpt2,1) ! nsppol implicitly = 1
       cgq_k(:,1:nband_k*nspinor*npw_k2) = &
&       cg(:,icg1+1:icg1+nband_k*nspinor*npw_k2)
       idum1 = dtbfield%fkgindex(ikpt2f)
       pwnsfac_k(3:4,1:npw_k2) = pwnsfac(1:2,idum1+1:idum1+npw_k2)
     end if

!    NOTE: kptopt 3 case seems robust, kptopt 4 may need more testing

     icp2=dtbfield%cprjindex(ikpt2,isppol)
     call pawcprj_get(gs_hamk%atindx1,cprj_kb,dtbfield%cprj,natom,1,icp2,ikpt,0,isppol,&
&     nband_k,nkpt,natom,nband_k,nband_k,nspinor,nsppol,0,&
&     mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

     if (ikpt2 /= ikpt2f) then ! construct cprj_kb by symmetry
       call pawcprj_copy(cprj_kb,cprj_ikn)
       call sym_pawcprj_kn(cprj_fkn,cprj_ikn,dtbfield%atom_indsym,dimlmn,-1,gs_hamk%indlmn,&
&       dtbfield%indkk_f2ibz(ikpt2f,2),dtbfield%indkk_f2ibz(ikpt2f,6),&
&       dtbfield%fkptns(:,dtbfield%i2fbz(ikpt2)),&
&       dtbfield%lmax,dtbfield%lmnmax,mband,natom,nband_k,nspinor,&
&       dtbfield%nsym,gs_hamk%ntypat,gs_hamk%typat,dtbfield%zarot)
       call pawcprj_copy(cprj_fkn,cprj_kb)
     end if

     call smatrix_k_paw_bfield(mod(idir,3)+1,0,cprj_k,cprj_kb,dtbfield,idir,ifor,mband,natom,&
&     smat_k_paw,gs_hamk%typat)

     icg1 = 0
     call smatrix(cg,cgq_k,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,&
&     job,nband_k,mcg,mcg_q,mcg1_k,1,mpw,nband_k,&
&     npw_k,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&     shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
     dtbfield%sflag(:,ikpt,ifor,idir) = sflag_k(:)
     dtbfield%smat(:,:,:,ikpt,ifor,idir) = smat_k(:,:,:)

!    now cg1_k contains |u}_k,k+b>, and these have also been shifted by pwind so
!    that they are expanded in exactly the same G vectors as |u_k>

!    !=======================================
!    ! code to test orthonormality of cg1_k
!    !=======================================
!    
!    ABI_ALLOCATE(bra,(2,npw_k*nspinor))
!    ABI_ALLOCATE(ket,(2,npw_k*nspinor))
!    avg_err_ovlp=0.0
!    max_err_ovlp=0.0
!    do bband = 1, nband_k
!    bra_start = 1+(bband-1)*npw_k*nspinor
!    bra_end = bra_start + npw_k*nspinor - 1
!    bra(1:2,1:npw_k*nspinor) = cg1_k(1:2,bra_start:bra_end)
!    do kband = 1, nband_k
!    ket_start = 1+(kband-1)*npw_k*nspinor
!    ket_end = ket_start + npw_k*nspinor - 1
!    ket(1:2,1:npw_k*nspinor) = cg1_k(1:2,ket_start:ket_end)
!    
!    tot_r = 0.0; tot_i = 0.0     
!    do ispinor = 1, nspinor
!    ovlp_r = 0.0; ovlp_i = 0.0
!    spnshft = (ispinor-1)*npw_k
!    do ipw = 1, npw_k
!    spnipw = ipw + spnshft
!    ovlp_r = ovlp_r + bra(1,spnipw)*ket(1,spnipw)+bra(2,spnipw)*ket(2,spnipw)
!    ovlp_i = ovlp_i - bra(2,spnipw)*ket(1,spnipw)+bra(1,spnipw)*ket(2,spnipw)
!    end do ! end loop over ipw
!    paw_r = 0.0
!    paw_i = 0.0
!    do iatom = 1, natom
!    itypat = gs_hamk%typat(iatom)
!    do ilmn = 1, dtbfield%lmn_size(itypat)
!    do jlmn = 1, dtbfield%lmn_size(itypat)
!    klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
!    bbs = nspinor*(bband-1)+ispinor
!    kbs = nspinor*(kband-1)+ispinor
!    cpb=cmplx(cprj_kb(iatom,bbs)%cp(1,ilmn),cprj_kb(iatom,bbs)%cp(2,ilmn))
!    cpk=cmplx(cprj_kb(iatom,kbs)%cp(1,jlmn),cprj_kb(iatom,kbs)%cp(2,jlmn))
!    cterm = conjg(cpb)*pawtab(itypat)%sij(klmn)*cpk
!    paw_r = paw_r + real(cterm)
!    paw_i = paw_i + aimag(cterm)
!    end do ! end loop over jlmn
!    end do ! end loop over ilmn
!    end do ! end loop over iatom
!    tot_r = tot_r + ovlp_r + paw_r
!    tot_i = tot_i + ovlp_i + paw_i
!    end do ! end loop over ispinor
!    !     write(std_out,'(a,5i4,2es16.8)')' JWZ Debug: update_orbmag ikpt idir ifor bband kband ovlp : ',&
!    !&           ikpt,idir,ifor,bband,kband,tot_r,tot_i
!    mag_ovlp = tot_r*tot_r + tot_i*tot_i
!    if(bband==kband) then
!    err_ovlp=abs(mag_ovlp-1.0)
!    else
!    err_ovlp=abs(mag_ovlp)
!    end if
!    max_err_ovlp=MAX(max_err_ovlp,err_ovlp)
!    avg_err_ovlp=avg_err_ovlp+err_ovlp
!    end do ! end loop over kband
!    end do ! end loop over bband
!    write(std_out,'(a,3i4,2es16.8)')' JWZ Debug: update_orbmag ikpt idir ifor avg_err max_err : ',&
!    &           ikpt,idir,ifor,avg_err_ovlp/(1.0*nband_k*nband_k),max_err_ovlp
!    ABI_DEALLOCATE(bra)
!    ABI_DEALLOCATE(ket)
!    
!    !=========================================
!    ! end code to test orthonormality of cg1_k
!    !=========================================

     tcg(:,:,idir,ifor) = cg1_k(:,:)
     call pawcprj_copy(cprj_kb,tcprj(:,:,idir,ifor))

   end do ! end loop over ifor

 end do ! end loop over idir

!compute emat(2,kband,ikpt) = <u_k|H_k|u_k>

 sij_opt = 0 ! compute <G|H|u> only, not gsc
 cpopt = 2 ! cpopt in memory
 lambda = 0.0 ! no shift to H used here
 ndat = 1 ! number of FFTs to do in parallel
 tim_getghc=8 
 
 type_calc = 0 ! use entire Hamiltonian
!
 ABI_ALLOCATE(kwave,(2,npw_k*nspinor))
 ABI_ALLOCATE(ghc,(2,npw_k*nspinor))
 ABI_ALLOCATE(gvnlc,(2,npw_k*nspinor))
 ABI_ALLOCATE(gsc,(2,npw_k*nspinor*ndat*(sij_opt+1)/2))

 ghc(:,:) = zero
 gvnlc(:,:) = zero

 ABI_DATATYPE_ALLOCATE(kcprj,(natom,nspinor*((cpopt+5)/5)))
 call pawcprj_alloc(kcprj,0,dimlmn_srt)

 ABI_ALLOCATE(kg_k,(3,mpw))
 ikg = dtbfield%kgindex(ikpt)
 kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

 do kband = 1, nband_k

   kwave(:,:) = cg(:,icg+(kband-1)*npw_k*nspinor+1:icg+kband*npw_k*nspinor)

!  copy cprj for kband into kcprj structure, change order from input atom order
!  to atom type sort order as needed by getghc
   call pawcprj_get(gs_hamk%atindx,kcprj,cprj_k,natom,kband,0,ikpt,1,1,&
&   nband_k,1,natom,1,nband_k,nspinor,nsppol,0,&
&   mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

   call  getghc(cpopt,kwave,kcprj,dimffnl,ffnl,filstat,ghc,gsc,&
&   gs_hamk,gvnlc,kg_k,kinpw,lambda,mpi_enreg,natom,ndat,npw_k,nspinor,&
&   paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,type_calc,vlocal,fock)

   dtbfield%emat(1,kband,ikpt) = dot_product(kwave(1,:),ghc(1,:)) + dot_product(kwave(2,:),ghc(2,:))
   dtbfield%emat(2,kband,ikpt) = dot_product(kwave(1,:),ghc(2,:)) - dot_product(kwave(2,:),ghc(1,:))

!  write(std_out,'(a,2i4,2es16.8)')' JWZ Debug: ikpt kband emat',ikpt,kband,&
!  &    dtbfield%emat(1,kband,ikpt),dtbfield%emat(2,kband,ikpt)

 end do

!compute omat(bband,kband,bdir,bfor,kdir,kfor) = <_bband,k+bdir*bfor|u_kband,k+kdir*kfor>
 ABI_ALLOCATE(omat,(nband_k,nband_k,3,2,3,2))
 ABI_ALLOCATE(bwave,(2,npw_k*nspinor))
 omat(:,:,:,:,:,:) = czero
 do bdir = 1, 3
   do kdir = 1, 3
     if (bdir == kdir) cycle ! never need omat for bdir // kdir
     do bfor = 1, 2
       do kfor = 1, 2

         call smatrix_k_paw_bfield(bdir,bfor,tcprj(:,:,bdir,bfor),tcprj(:,:,kdir,kfor),dtbfield,&
&         kdir,kfor,mband,natom,smat_k_paw,gs_hamk%typat)

         do bband = 1, nband_k
           bwave(:,:) = tcg(:,(bband-1)*nspinor*npw_k+1:bband*nspinor*npw_k,bdir,bfor)

           do kband = 1, nband_k

             kwave(:,:) = tcg(:,(kband-1)*nspinor*npw_k+1:kband*nspinor*npw_k,kdir,kfor)

             dotri(1) = dot_product(bwave(1,:),kwave(1,:))+dot_product(bwave(2,:),kwave(2,:))
             dotri(2) = -dot_product(bwave(2,:),kwave(1,:))+dot_product(bwave(1,:),kwave(2,:))

             omat(bband,kband,bdir,bfor,kdir,kfor) = cmplx(dotri(1)+smat_k_paw(1,bband,kband),&
&             dotri(2)+smat_k_paw(2,bband,kband))

           end do ! end loop over kband
         end do ! end loop over bband
       end do ! end loop over kfor
     end do ! end loop over bfor
   end do ! end loop over kdir
 end do ! end loop over bdir

!!compute chern_k(2,ikpt,idir)
!
 dtbfield%chern_k(:,ikpt,:) = zero

 do idir = 1, 3
   bdir = mod(idir,3)+1
   kdir = mod(bdir,3)+1
   sfac = 0.25
   do istep = 1, 2
     do bsig = -1, 1, 2
       do ksig= -1, 1, 2
         bfor = (-bsig+3)/2; kfor = (-ksig+3)/2

         IA = cmplx(zero,zero)
         IB = cmplx(zero,zero)
         do nn = 1, nband_k
           do nnp = 1, nband_k
             IA1 = cmplx(dtbfield%smat(1,nn,nnp,ikpt,bfor,bdir),&
&             dtbfield%smat(2,nn,nnp,ikpt,bfor,bdir))
             IB1 = IA1
             do nnpp = 1, nband_k
               IA2 = omat(nnp,nnpp,bdir,bfor,kdir,kfor)
               IA3 = cmplx(dtbfield%smat(1,nn,nnpp,ikpt,kfor,kdir),&
&               dtbfield%smat(2,nn,nnpp,ikpt,kfor,kdir))
               IA = IA + IA1*IA2*conjg(IA3)

               IB2 = cmplx(dtbfield%smat(1,nnpp,nnp,ikpt,bfor,bdir),&
&               dtbfield%smat(2,nnpp,nnp,ikpt,bfor,bdir))

               do nnppp = 1, nband_k
                 IB3 = cmplx(dtbfield%smat(1,nnpp,nnppp,ikpt,kfor,kdir),&
&                 dtbfield%smat(2,nnpp,nnppp,ikpt,kfor,kdir))
                 IB4 = cmplx(dtbfield%smat(1,nn,nnppp,ikpt,kfor,kdir),&
&                 dtbfield%smat(2,nn,nnppp,ikpt,kfor,kdir))
                 IB = IB + IB1*conjg(IB2)*IB3*conjg(IB4)
               end do ! end loop over nnppp
             end do ! end loop over nnpp
           end do ! end loop over nnp
         end do ! end loop over nn

         dtbfield%chern_k(1,ikpt,idir) = &
&         dtbfield%chern_k(1,ikpt,idir) + sfac*bsig*ksig*real(IA-IB)
         dtbfield%chern_k(2,ikpt,idir) = &
&         dtbfield%chern_k(2,ikpt,idir) + sfac*bsig*ksig*aimag(IA-IB)

       end do ! end loop over ksig
     end do ! end loop over bsig

     sfac = -sfac
     idum1 = bdir; bdir = kdir; kdir = idum1

   end do ! loop over istep

 end do ! loop over idir

!!compute mag_k(2,ikpt,idir)
!

 sij_opt = 0 ! compute <G|H|u> only, not gsc
 cpopt = 2 ! cpopt in memory
 lambda = 0.0 ! no shift to H used here
 ndat = 1 ! number of FFTs to do in parallel
 tim_getghc=8 
 
 type_calc = 3 ! use only local + kinetic (phase twisted Dij will be applied following)

 dtbfield%mag_k(:,ikpt,:) = zero

 do idir = 1, 3
   bdir = mod(idir,3)+1
   kdir = mod(bdir,3)+1
   sfac = 0.25
   do istep = 1, 2
     do bsig = -1, 1, 2
       do ksig= -1, 1, 2

         indhk = dtbfield%indhk(PIND(bdir,bsig),PIND(kdir,ksig))
         bfor = (-bsig+3)/2; kfor = (-ksig+3)/2

         do nnpp = 1, nband_k

           kwave(:,:) = tcg(:,(nnpp-1)*nspinor*npw_k+1:nnpp*nspinor*npw_k,bdir,bfor)

!          copy cprj for nnpp into kcprj structure, change order from input atom order
!          to atom type sort order as needed by getghc
           call pawcprj_get(gs_hamk%atindx,kcprj,cprj_k,natom,nnpp,0,ikpt,1,1,&
&           nband_k,1,natom,1,nband_k,nspinor,nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

           call  getghc(cpopt,kwave,kcprj,dimffnl,ffnl,filstat,ghc,gsc,&
&           gs_hamk,gvnlc,kg_k,kinpw,lambda,mpi_enreg,natom,ndat,npw_k,nspinor,&
&           paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,type_calc,vlocal,fock)

           do nnp = 1, nband_k
             IIIA3 = omat(nnp,nnpp,bdir,bfor,kdir,kfor)

             bwave(:,:) = tcg(:,(nnp-1)*nspinor*npw_k+1:nnp*nspinor*npw_k,kdir,kfor)
             dotri(1)=dot_product(bwave(1,:),ghc(1,:))+dot_product(bwave(2,:),ghc(2,:))
             dotri(2)=dot_product(bwave(1,:),ghc(2,:))-dot_product(bwave(2,:),ghc(1,:))
             do iatom = 1, natom
               do idij = 1, dtbfield%ndij
                 select case (idij)
                   case (1) 
                     bbs = nspinor*(nnp-1)+1
                     kbs = nspinor*(nnpp-1)+1
                   case (2) 
                     bbs = nspinor*(nnp-1)+2
                     kbs = nspinor*(nnpp-1)+2
                   case (3) 
                     bbs = nspinor*(nnp-1)+1
                     kbs = nspinor*(nnpp-1)+2
                   case (4) 
                     bbs = nspinor*(nnp-1)+2
                     kbs = nspinor*(nnpp-1)+1
                 end select
                 do blmn = 1, dimlmn(iatom)
                   cpb = cmplx(tcprj(iatom,bbs,kdir,kfor)%cp(1,blmn),&
&                   tcprj(iatom,bbs,kdir,kfor)%cp(2,blmn))
                   do klmn = 1, dimlmn(iatom)
                     cpk = cmplx(tcprj(iatom,kbs,bdir,bfor)%cp(1,klmn),&
&                     tcprj(iatom,kbs,bdir,bfor)%cp(2,klmn))
                     twdij = cmplx(dtbfield%twdij(1,blmn,klmn,iatom,indhk,idij),&
&                     dtbfield%twdij(2,blmn,klmn,iatom,indhk,idij))
                     dterm = conjg(cpb)*twdij*cpk
                     dotri(1) = dotri(1) + real(dterm)
                     dotri(2) = dotri(2) + aimag(dterm)
                   end do ! end loop over ket lmn
                 end do ! end loop over bra lmn
               end do ! end loop over idij
             end do ! end loop over atoms

             IIA3=cmplx(dotri(1),dotri(2))

             do nn = 1, nband_k

               IIIA1 = cmplx(dtbfield%emat(1,nn,ikpt),&
&               dtbfield%emat(2,nn,ikpt))

               IIIA2 = cmplx(dtbfield%smat(1,nn,nnp,ikpt,bfor,bdir),&
&               dtbfield%smat(2,nn,nnp,ikpt,bfor,bdir))

               IIIA4 = cmplx(dtbfield%smat(1,nn,nnpp,ikpt,kfor,kdir),&
&               dtbfield%smat(2,nn,nnpp,ikpt,kfor,kdir))

               IIIA = -IIIA1*IIIA2*IIIA3*conjg(IIIA4)

               IIA2 = cmplx(dtbfield%smat(1,nn,nnp,ikpt,kfor,kdir),&
&               dtbfield%smat(2,nn,nnp,ikpt,kfor,kdir))

               IIA4 = cmplx(dtbfield%smat(1,nn,nnpp,ikpt,bfor,bdir),&
&               dtbfield%smat(2,nn,nnpp,ikpt,bfor,bdir))

               IIA = IIA2*IIA3*conjg(IIA4)

             end do ! end loop over nn
           end do ! end loop over nnp
         end do ! end loop over nnpp

         dtbfield%mag_k(1,ikpt,idir) = &
&         dtbfield%mag_k(1,ikpt,idir) + sfac*bsig*ksig*real(IIA+IIIA)
         dtbfield%mag_k(2,ikpt,idir) = &
&         dtbfield%mag_k(2,ikpt,idir) + sfac*bsig*ksig*aimag(IIA+IIIA)

       end do ! end loop over ksig
     end do ! end loop over bsig

     sfac = -sfac
     idum1 = bdir; bdir = kdir; kdir = idum1

   end do ! loop over istep

 end do ! loop over idir

 
 ABI_DEALLOCATE(kwave)
 ABI_DEALLOCATE(bwave)
 ABI_DEALLOCATE(omat)
 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(gvnlc)
 ABI_DEALLOCATE(gsc)
 call pawcprj_destroy(kcprj)
 ABI_DATATYPE_DEALLOCATE(kcprj)

 ABI_DEALLOCATE(tcg)
 ABI_DEALLOCATE(cg1_k)
 ABI_DEALLOCATE(cgq_k)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(pwind_k)
 ABI_DEALLOCATE(pwnsfac_k)
 ABI_DEALLOCATE(dimlmn)
 ABI_DEALLOCATE(dimlmn_srt)
 ABI_DEALLOCATE(smat_k)
 ABI_DEALLOCATE(smat_inv)
 ABI_DEALLOCATE(smat_k_paw)

 call pawcprj_destroy(cprj_k)
 call pawcprj_destroy(cprj_kb)
 ABI_DATATYPE_DEALLOCATE(cprj_k)
 ABI_DATATYPE_DEALLOCATE(cprj_kb)

 do idir = 1, 3
   do ifor = 1, 2
     call pawcprj_destroy(tcprj(:,:,idir,ifor))
   end do
 end do
 ABI_DATATYPE_DEALLOCATE(tcprj)

 if (nkpt /= dtbfield%fnkpt) then
   call pawcprj_destroy(cprj_fkn)
   call pawcprj_destroy(cprj_ikn)
   ABI_DATATYPE_DEALLOCATE(cprj_fkn)
   ABI_DATATYPE_DEALLOCATE(cprj_ikn)
 end if



!DEBUG
!write(std_out,'(a)')' leaving update_orbmag '
!END DEBUG

end subroutine update_orbmag
!!***
