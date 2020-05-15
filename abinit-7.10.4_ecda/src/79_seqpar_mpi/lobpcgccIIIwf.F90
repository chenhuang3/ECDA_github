!{\src2tex{textfont=tt}}
!!****f* ABINIT/lobpcgccIIIwf
!! NAME
!! lobpcgccIIIwf
!!
!! FUNCTION
!! TO BE DESCRIBED 090830
!!
!! PARENTS
!!      lobpcgccIIwf
!!
!! CHILDREN
!!      abi_xorthonormalize,cg_zprecon_block,getghc,nonlop,zgemm,zhegv,ztrsm
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine lobpcgccIIIwf(dimffnl,dtfil,dtset,ffnl,gs_hamk,iterationnumber,&
     & kg_k,kinpw,mpi_enreg,natom,npw_k,&
     & pcon,ph3d,prtvol,vlocal,&
     & blocksize,bblocksize,vectsize,pflag,&
     & blockvectorx,blockvectorbx,blockvectorax,blockvectory,blockvectorby,lambda,&
     & blockvectorp,blockvectorbp,blockvectorap,&
     & vxctaulocal) ! optional argument

 use m_profiling_abi
 use defs_basis
 use defs_abitypes
 use m_abi_linalg
 use m_cgtools

 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj, only : pawcprj_type
 use m_fock,    only : fock_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lobpcgccIIIwf'
 use interfaces_65_nonlocal
 use interfaces_66_wfs
!End of the abilint section

  implicit none

!Arguments -------------------------------

  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  integer,intent(in) :: dimffnl,natom,npw_k,prtvol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: kg_k(3,npw_k)
  integer,intent(in) :: bblocksize
  real(dp),intent(in) :: ffnl(npw_k,dimffnl,gs_hamk%lmnmax,gs_hamk%ntypat)
  real(dp),intent(in) :: kinpw(npw_k)
  real(dp),intent(inout) :: ph3d(2,npw_k,gs_hamk%matblk)
  real(dp),intent(inout) :: vlocal(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc)
  integer,intent(in) :: blocksize,vectsize
  integer,intent(in) :: iterationnumber
  real(dp),intent(inout) :: lambda
  real(dp),intent(inout) :: pcon(npw_k,blocksize)
  complex(dp),intent(out) :: blockvectorx(vectsize,blocksize),blockvectorax(vectsize,blocksize),&
&                             blockvectorbx(vectsize,blocksize),blockvectory(vectsize,bblocksize),&
&                             blockvectorby(vectsize,bblocksize)!,lambda(blocksize,blocksize)
  complex(dp),intent(inout) :: blockvectorp(vectsize,blocksize),blockvectorap(vectsize,blocksize),&
&                               blockvectorbp(vectsize,blocksize)
  logical,intent(inout) :: pflag(blocksize)
  real(dp), intent(inout), optional :: vxctaulocal(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc,4)

!Local variables -------------------------
  real(dp) :: sq2
  integer:: activersize,bigorder,choice,cond_try,cpopt,iblocksize
! integer :: activepsize,cgindex ! appears in comment below
  integer:: i1,i2,i3,i4,idir,info
  integer:: istwf_k,lwork
  integer:: maxiterations,my_nspinor,nnlout,nkpg,optekin,optpcon,paw_opt
  integer:: restart,signs,sij_opt,tim_getghc,tim_nonlop
  logical::gen_eigenpb
  real(dp) :: dum
  real(dp) :: dummy_lambda(1)
  type(fock_type),pointer  :: fock => null()
  complex(dp) ::  blockvectorz(vectsize,blocksize),blockvectoraz(vectsize,blocksize),&
       &                 blockvectorbz(vectsize,blocksize),blockvectorr1(vectsize,blocksize)
  real(dp) :: residualnorms1(blocksize)
  real(dp), allocatable :: dummy1(:),dummy2(:,:),dummy3(:,:,:)
  complex(dp) :: zlambda(1,1)
  type(pawcprj_type) :: cprj_dum(1,1)
  complex(dp), allocatable ::  blockvectorr(:,:),blockvectorar(:,:),blockvectorbr(:,:),&
       & blockvectordumm(:,:),&
       & gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),&
       & grampap(:,:),&
       & gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),&
       & grampbp(:,:),&
       & coordx(:,:),&
       & grama(:,:),gramb(:,:),gramyx(:,:),&
       & work(:)
  real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlc(:,:)
  real(dp), allocatable :: residualnorms(:),eigen(:),rwork(:),kpg_dum(:,:)
  real(dp), parameter :: tolerance1 = 1.e-13, tolerance2 = 1.e2

! *********************************************************************

!correspondence with abinit. here for real wf but in complex mode
!this is the index of a given band
!cgindex(iblocksize)=npw_k*my_nspinor*(iblocksize-1)+icg+1
 gen_eigenpb=(gs_hamk%usepaw==1)
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 optekin=0;if (dtset%wfoptalg>10) optekin=0
 optpcon=1;if (dtset%wfoptalg>10) optpcon=0
 sq2=sqrt(two)
!vectsize=npw_k*my_nspinor
!blocksize=(nband_k-1)/nbdblock+1
!bblocksize=(iblock-1)*blocksize
 istwf_k=gs_hamk%istwf_k
 maxiterations=dtset%nline



!passing x into z
 blockvectorz = blockvectorx
 blockvectoraz = blockvectorax
 blockvectorbz = blockvectorbx

!allocations
 ABI_ALLOCATE(blockvectorr,(vectsize,blocksize))
 ABI_ALLOCATE(blockvectorar,(vectsize,blocksize))
 ABI_ALLOCATE(blockvectorbr,(vectsize,blocksize))
 ABI_ALLOCATE(blockvectordumm,(vectsize,blocksize))
 ABI_ALLOCATE(gramyx,(bblocksize,blocksize))
 ABI_ALLOCATE(gramxax,(blocksize,blocksize))
 ABI_ALLOCATE(gramxar,(blocksize,blocksize))
 ABI_ALLOCATE(gramxap,(blocksize,blocksize))
 ABI_ALLOCATE(gramrar,(blocksize,blocksize))
 ABI_ALLOCATE(gramrap,(blocksize,blocksize))
 ABI_ALLOCATE(grampap,(blocksize,blocksize))
 ABI_ALLOCATE(gramxbx,(blocksize,blocksize))
 ABI_ALLOCATE(gramxbr,(blocksize,blocksize))
 ABI_ALLOCATE(gramxbp,(blocksize,blocksize))
 ABI_ALLOCATE(gramrbr,(blocksize,blocksize))
 ABI_ALLOCATE(gramrbp,(blocksize,blocksize))
 ABI_ALLOCATE(grampbp,(blocksize,blocksize))
 ABI_ALLOCATE(residualnorms,(blocksize))
!write(std_out,*)'iterationnumber',iterationnumber
!construct residual
!blockvectorr=blockvectorax-matmul(blockvectorx,lambda)
 do iblocksize=1,blocksize
   zlambda(1,1)=lambda
   call cg_zprecon_block(blockvectorbx(:,iblocksize),zlambda,1,& !,iblocksize),&
&   iterationnumber,kinpw,npw_k,my_nspinor,&
&   optekin,optpcon,pcon,blockvectorax(:,iblocksize),blockvectorr(:,iblocksize),vectsize,mpi_enreg%comm_bandspinorfft)

!  blockvectorr(:,iblocksize)=blockvectorax(:,iblocksize)-lambda(iblocksize,iblocksize)*blockvectorbx(:,iblocksize)
 end do
 residualnorms=sqrt(sum(abs(blockvectorr)**2,dim=1))
 write(std_out,*)'residualnorm',residualnorms
!if(abs(sum(residualnorms)) < 1.d-10) exit
!not yet masking

!DEBUG le if suivant ne marche que si blocksize = 1
 if (residualnorms(1)>tolerance1) then !sinon on masque

   if(bblocksize>0) then !residuals orthogonal to blockvectorby
!    blockvectorr=blockvectorr-&
!    &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorr))
     call zgemm('c','n',bblocksize,blocksize,vectsize,cone,blockvectorby,&
&     vectsize,blockvectorr,vectsize,czero,gramyx,bblocksize)
     call zgemm('n','n',vectsize,blocksize,bblocksize,cone,blockvectory,&
&     vectsize,gramyx,bblocksize,czero,blockvectordumm,vectsize)
     blockvectorr=blockvectorr-blockvectordumm
   end if
!  residuals orthogonal to blockvectorx
!  blockvectorr=blockvectorr-&
!  &matmul(blockvectorx,matmul(transpose(blockvectorbx),blockvectorr))
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&   vectsize,blockvectorr,vectsize,czero,gramxax,blocksize)
   call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&   vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
   blockvectorr=blockvectorr-blockvectordumm
!  and now (b)orthornormalize r
!  if(iterationnumber >1) stop('xx2')
!  call operators(blockvectorr,blockvectorbr)
   ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
   ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
   do iblocksize=1,blocksize
     cwavef(1,:)=real (blockvectorr(1:npw_k*my_nspinor,iblocksize))
     cwavef(2,:)=aimag(blockvectorr(1:npw_k*my_nspinor,iblocksize))
     if(gen_eigenpb) then
!      Call to nonlop: compute <g|S|c>
       choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
       call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
&       dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
&       istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dummy_lambda,gs_hamk%lmnmax,gs_hamk%matblk,&
&       gs_hamk%mgfft,mpi_enreg,gs_hamk%mpsang,gs_hamk%mpssoang,natom,gs_hamk%nattyp,1,gs_hamk%ngfft,nkpg,nkpg,&
&       gs_hamk%nloalg,nnlout,npw_k,npw_k,my_nspinor,dtset%nspinor,gs_hamk%ntypat,0,paw_opt,gs_hamk%phkxred,&
&       gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,signs,gs_hamk%sij,&
&       gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,use_gpu_cuda=dtset%use_gpu_cuda)
     else
       gwavef(:,:)=cwavef(:,:)
     end if
     blockvectorbr(1:npw_k*my_nspinor,iblocksize)=dcmplx(gwavef(1,:),gwavef(2,:))
   end do
   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(gwavef)
!  call zorthonormalize(blockvectorr,blockvectorbr)
   call abi_xorthonormalize(blockvectorr,blockvectorbr,blocksize,mpi_enreg%comm_bandspinorfft,gramrbr,vectsize)
   call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramrbr,blocksize,&
&   blockvectorbr,vectsize)
!  compute ar
!  blockvectorar=matmul(operatora,blockvectorr)
!  if(iterationnumber >1) stop('xx3')
!  call operatorh(blockvectorr,blockvectorar)
   ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
   ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor))
   ABI_ALLOCATE(gvnlc,(2,npw_k*my_nspinor))
   do iblocksize=1,blocksize
     cwavef(1,:)=real(blockvectorr(1:npw_k*my_nspinor,iblocksize))
     cwavef(2,:)=aimag(blockvectorr(1:npw_k*my_nspinor,iblocksize))
     tim_getghc=7; sij_opt=0
     if(present(vxctaulocal))then
       call getghc(-1,cwavef,cprj_dum,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
&       kinpw,dum,mpi_enreg,natom,blocksize,npw_k,my_nspinor,dtset%paral_kgb,ph3d,prtvol,&
&       sij_opt,tim_getghc,0,vlocal,fock,vxctaulocal=vxctaulocal)
     else
       call getghc(-1,cwavef,cprj_dum,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
&       kinpw,dum,mpi_enreg,natom,blocksize,npw_k,my_nspinor,dtset%paral_kgb,ph3d,prtvol,&
&       sij_opt,tim_getghc,0,vlocal,fock)
     end if
     blockvectorar(1:npw_k*my_nspinor,iblocksize)=dcmplx(gwavef(1,:),gwavef(2,:))
   end do
   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(gwavef)
   ABI_DEALLOCATE(gvnlc)
   if(pflag(1)) then  !ce n'est pas la bonne condition si blocksize /= 1
!    write(std_out,*)'blockvectorp,blockvectorbp,blockvectorap'
!    write(std_out,*)blockvectorp
!    write(std_out,*)blockvectorbp
!    write(std_out,*)blockvectorap
!    if(iterationnumber >1) stop('xx4')
!    call zorthonormalize(blockvectorp,blockvectorbp,blockvectorap)
     call abi_xorthonormalize(blockvectorp,blockvectorbp,blocksize,mpi_enreg%comm_bandspinorfft,grampbp,vectsize)
     call ztrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,&
&     blockvectorbp,vectsize)
     call ztrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,&
&     blockvectorap,vectsize)

!    if(iterationnumber >1) stop('xx4')
!    blockvectorap=matmul(blockvectorap,grampbp)
   end if
   activersize=blocksize
   if (.not.pflag(1)) then   !ce n'est pas la bonne condition si blocksize /= 1
!    activepsize=0
     restart=1
   else
!    activepsize=blocksize
     restart=0
   end if
!  gramxar=matmul(transpose(blockvectorax),blockvectorr)
!  gramrar=matmul(transpose(blockvectorar),blockvectorr)
!  gramxax=matmul(transpose(blockvectorax),blockvectorx)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorax,&
&   vectsize,blockvectorr,vectsize,czero,gramxar,blocksize)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorar,&
&   vectsize,blockvectorr,vectsize,czero,gramrar,blocksize)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorax,&
&   vectsize,blockvectorx,vectsize,czero,gramxax,blocksize)
!  gramxbx=matmul(transpose(blockvectorbx),blockvectorx)
!  gramrbr=matmul(transpose(blockvectorbr),blockvectorr)
!  gramxbr=matmul(transpose(blockvectorbx),blockvectorr)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&   vectsize,blockvectorx,vectsize,czero,gramxbx,blocksize)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbr,&
&   vectsize,blockvectorr,vectsize,czero,gramrbr,blocksize)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&   vectsize,blockvectorr,vectsize,czero,gramxbr,blocksize)
!  if(iterationnumber >1) stop('xx1')
   i1=0;i2=blocksize;i3=2*blocksize;i4=3*blocksize
   cond: do cond_try=1,2 !2 when restart, but not implemented
     if (restart==0) then
!      gramxap=matmul(transpose(blockvectorax),blockvectorp)
!      gramrap=matmul(transpose(blockvectorar),blockvectorp)
!      grampap=matmul(transpose(blockvectorap),blockvectorp)
       call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorax,&
&       vectsize,blockvectorp,vectsize,czero,gramxap,blocksize)
       call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorar,&
&       vectsize,blockvectorp,vectsize,czero,gramrap,blocksize)
       call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorap,&
&       vectsize,blockvectorp,vectsize,czero,grampap,blocksize)
       bigorder=i4
       ABI_ALLOCATE(grama,(i4,i4))
       ABI_ALLOCATE(gramb,(i4,i4))
       ABI_ALLOCATE(eigen,(i4))
       ABI_ALLOCATE(coordx,(i4,blocksize))
       grama(i1+1:i2,i1+1:i2)=gramxax
       grama(i1+1:i2,i2+1:i3)=gramxar
       grama(i1+1:i2,i3+1:i4)=gramxap
!      grama(i2+1:i3,i1+1:i2)=transpos(gramxar)
       grama(i2+1:i3,i2+1:i3)=gramrar
       grama(i2+1:i3,i3+1:i4)=gramrap
!      grama(i3+1:i4,i1+1:i2)=transpos(gramxap)
!      grama(i3+1:i4,i2+1:i3)=transpos(gramrap)
       grama(i3+1:i4,i3+1:i4)=grampap

!      gramxbp=matmul(transpose(blockvectorbx),blockvectorp)
!      gramrbp=matmul(transpose(blockvectorbr),blockvectorp)
!      grampbp=matmul(transpose(blockvectorbp),blockvectorp)
       call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&       vectsize,blockvectorp,vectsize,czero,gramxbp,blocksize)
       call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbr,&
&       vectsize,blockvectorp,vectsize,czero,gramrbp,blocksize)
       call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbp,&
&       vectsize,blockvectorp,vectsize,czero,grampbp,blocksize)

       gramb(i1+1:i2,i1+1:i2)=gramxbx
       gramb(i1+1:i2,i2+1:i3)=gramxbr
       gramb(i1+1:i2,i3+1:i4)=gramxbp
!      gramb(i2+1:i3,i1+1:i2)=transpos(gramxbr)
       gramb(i2+1:i3,i2+1:i3)=gramrbr
       gramb(i2+1:i3,i3+1:i4)=gramrbp
!      gramb(i3+1:i4,i1+1:i2)=transpos(gramxbp)
!      gramb(i3+1:i4,i2+1:i3)=transpos(gramrbp)
       gramb(i3+1:i4,i3+1:i4)=grampbp
     else
       bigorder=i3
       ABI_ALLOCATE(grama,(i3,i3))
       ABI_ALLOCATE(gramb,(i3,i3))
       ABI_ALLOCATE(eigen,(i3))
       ABI_ALLOCATE(coordx,(i3,blocksize))
       grama(i1+1:i2,i1+1:i2)=gramxax
       grama(i1+1:i2,i2+1:i3)=gramxar
!      grama(i2+1:i3,i1+1:i2)=transpos(gramxar)
       grama(i2+1:i3,i2+1:i3)=gramrar
       gramb(i1+1:i2,i1+1:i2)=gramxbx
       gramb(i1+1:i2,i2+1:i3)=gramxbr
!      gramb(i2+1:i3,i1+1:i2)=transpos(gramxbr)
       gramb(i2+1:i3,i2+1:i3)=gramrbr
     end if
!    !!     end do cond
!    call la_sygv(grama,gramb,eigen,itype=1,jobz='v')
!    if(iterationnumber >1) stop('xx')
     lwork=3*bigorder-2
     ABI_ALLOCATE(work,(lwork))
     ABI_ALLOCATE(rwork,(lwork))
     call zhegv(1,'v','u',bigorder,grama,bigorder,gramb,bigorder,eigen,&
&     work,lwork,rwork,info)
     ABI_DEALLOCATE(work)
     ABI_DEALLOCATE(rwork)
     do iblocksize=1,blocksize
       lambda=eigen(iblocksize)
     end do
     write(std_out,*)'eigen',eigen(1:blocksize)
     coordx=grama(:,1:blocksize)
     ABI_DEALLOCATE(grama)
     ABI_DEALLOCATE(gramb)
     ABI_DEALLOCATE(eigen)
     if (restart==0) then
!      blockvectorp=matmul(blockvectorr,coordx(i2+1:i3,:))+&
!      &matmul(blockvectorp,coordx(i3+1:i4,:))
       call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorr,&
&       vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
       call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorp,&
&       vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
       blockvectorp=blockvectordumm
!      blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))+&
!      &matmul(blockvectorap,coordx(i3+1:i4,:))
       call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorar,&
&       vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
       call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorap,&
&       vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
       blockvectorap=blockvectordumm
!      blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))+&
!      &matmul(blockvectorbp,coordx(i3+1:i4,:))
       call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbr,&
&       vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
       call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbp,&
&       vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
       blockvectorbp=blockvectordumm
     else
!      blockvectorp =matmul(blockvectorr,coordx(i2+1:i3,:))
       call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorr,&
&       vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorp,vectsize)
!      blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))
       call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorar,&
&       vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorap,vectsize)
!      blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))
       call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbr,&
&       vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorbp,vectsize)
     end if

!    blockvectorx = matmul(blockvectorx,coordx(i1+1:i2,:))+blockvectorp
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&     vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
     blockvectorx = blockvectordumm+blockvectorp
!    blockvectorax= matmul(blockvectorax,coordx(i1+1:i2,:))+blockvectorap
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
&     vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
     blockvectorax = blockvectordumm+blockvectorap
!    blockvectorbx= matmul(blockvectorbx,coordx(i1+1:i2,:))+blockvectorbp
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
&     vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
     blockvectorbx = blockvectordumm+blockvectorbp
     ABI_DEALLOCATE(coordx)

     do iblocksize=1,blocksize
       zlambda(1,1) = lambda
       call cg_zprecon_block(blockvectorbx(:,iblocksize),zlambda,1,& !iblocksize),&
&       iterationnumber,kinpw,npw_k,my_nspinor,&
&       optekin,optpcon,pcon,blockvectorax(:,iblocksize),blockvectorr1(:,iblocksize),vectsize,mpi_enreg%comm_bandspinorfft)
     end do
     residualnorms1=sqrt(sum(abs(blockvectorr1)**2,dim=1))
     write(std_out,*)'residualnorm after lobpcgcciii',residualnorms1

     do iblocksize=1,blocksize
       if (residualnorms1(iblocksize) > tolerance2*residualnorms(iblocksize)) then
         blockvectorx = blockvectorz
         blockvectorax = blockvectoraz
         blockvectorbx = blockvectorbz
         if (cond_try == 1) then
           write(std_out,*)'restart here for this eig'
           restart = 1
         else
           blockvectorp = czero
           blockvectorap = czero
           blockvectorbp = czero
           pflag = .false.
           write(std_out,*)'restart was unuseful'
         end if
       else
         pflag = .true.
         exit cond   !DEBUG exact que si blocksize = 1 (sinon il faut tester les autres )
       end if
     end do
   end do cond

 else
   blockvectorp = czero
   blockvectorap = czero
   blockvectorbp = czero
   pflag = .false.
 end if

!write(std_out,*)'blockvectorr',blockvectorr
!write(std_out,*)'blockvectorx',blockvectorx
!write(std_out,*)'blockvectorax',blockvectorax
 ABI_DEALLOCATE(blockvectorr)
 ABI_DEALLOCATE(blockvectorar)
 ABI_DEALLOCATE(blockvectorbr)
 ABI_DEALLOCATE(gramyx)
 ABI_DEALLOCATE(blockvectordumm)
 ABI_DEALLOCATE(gramxax)
 ABI_DEALLOCATE(gramxar)
 ABI_DEALLOCATE(gramxap)
 ABI_DEALLOCATE(gramrar)
 ABI_DEALLOCATE(gramrap)
 ABI_DEALLOCATE(grampap)
 ABI_DEALLOCATE(gramxbx)
 ABI_DEALLOCATE(gramxbr)
 ABI_DEALLOCATE(gramxbp)
 ABI_DEALLOCATE(gramrbr)
 ABI_DEALLOCATE(gramrbp)
 ABI_DEALLOCATE(grampbp)
 ABI_DEALLOCATE(residualnorms)
end subroutine lobpcgccIIIwf
!!***
