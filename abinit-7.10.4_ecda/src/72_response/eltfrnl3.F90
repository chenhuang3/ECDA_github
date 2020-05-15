!{\src2tex{textfont=tt}}
!!****f* ABINIT/eltfrnl3
!! NAME
!! eltfrnl3
!!
!! FUNCTION
!! Compute the frozen-wavefunction non-local contribution to the elastic tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DRH, DCA, XG, GM, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>
!!              =Fourier coefficients of wavefunction
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  kptns(3,nkpt)=coordinates of k points in terms of reciprocal space
!!   primitive translations
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimension for number of planewaves
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nband(nkpt*nsppol)=number of bands being considered per k point
!!  nkpt=number of k points
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  ntypat=integer specification of atom type (1, 2, ...)
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2)
!!    at each k point
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  wtk(nkpt)=k point weights
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,9,mpsang*mpsang*useylm)= gradients of real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  eltfrnl(6+3*natom,6)=non-symmetrized non-local contribution to the
!!                    elastic tensor
!!
!! PARENTS
!!
!! CHILDREN
!!      metric,mkffnl,mkkpg,nonlop,ph1d3d,sphereboundary,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eltfrnl3(atindx,atindx1,cg,eltfrnl,istwfk,&
&  kg,kptns,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nband,&
&  nkpt,ngfft,nloalg,npwarr,nspinor,nsppol,ntypat,occ,ph1d,&
&  psps,rprimd,useylm,wtk,xred,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_wffile

 use m_pawcprj, only : pawcprj_type
 use m_header,  only : hdr_skip

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eltfrnl3'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_52_fft_mpi_noabirule
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mkmem,mpsang,mpw,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: ntypat,useylm
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),istwfk(nkpt)
 integer,intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),nloalg(5),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),wtk(nkpt),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,9,mpsang*mpsang*useylm)
 real(dp),intent(out) :: eltfrnl(6+3*natom,6)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,choice,cpopt,dimekb1,dimekb2,dimffnl,ia,iatom
 integer :: iband,icg,ider,idir,ielt,ieltx,ierr,ii,ikg,ikpt,ilm
 integer :: index,ipw,isppol,istwf_k,jj,master,matblk,me,n1,n2,n3
 integer :: nband_k,nkpg,nnlout,npw_k,paw_opt,signs,spaceComm
 integer :: tim_nonlop
 real(dp) :: arg,enl,enlk,ucvol
!arrays
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 real(dp) :: dum(1)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),dum_sij(1,1),dum_svectout(1,1)
 real(dp) :: rmet(3,3),tsec(2)
 real(dp),allocatable :: cwavef(:,:),ekb(:,:,:)
 real(dp),allocatable :: elt_work(:,:),eltfrnlk(:,:),enlout(:),ffnl(:,:,:,:)
 real(dp),allocatable :: kpg_k(:,:),ph3d(:,:,:),phkxred(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

!Default for sequential use
 master=0
 spaceComm=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt

 nnlout=6*(3*natom+6)

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 enl=0.0_dp
 eltfrnl(:,:)=0.0_dp
 bdtot_index=0
 icg=0

!Non-local factors:
!Norm-conserving: kleimann-Bylander energies
 if (psps%usepaw==0) then
   dimekb1=psps%dimekb;dimekb2=ntypat
   ABI_ALLOCATE(ekb,(psps%dimekb,ntypat,nspinor**2))
   ekb(:,:,1)=psps%ekb(:,:)
   if (nspinor==2) then
     ekb(:,:,2)=psps%ekb(:,:)
     ekb(:,:,3:4)=zero
   end if
 else
!  Not available within PAW
   ABI_ALLOCATE(ekb,(psps%dimekb,natom,nspinor**2))
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(cwavef,(2,mpw*nspinor))
 ABI_ALLOCATE(eltfrnlk,(6+3*natom,6))

!Define k-points distribution

!LOOP OVER SPINS
 do isppol=1,nsppol

   ikg=0

!  Loop over k points
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)

!    Skip this k-point if not the proper processor
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*useylm))
     ABI_ALLOCATE(ylmgr_k,(npw_k,9,mpsang*mpsang*psps%useylm))
     kpoint(:)=kptns(:,ikpt)

     kg_k(:,:) = 0

!$OMP PARALLEL DO PRIVATE(ipw) SHARED(ikg,kg,kg_k,npw_k)
     do ipw=1,npw_k
       kg_k(1,ipw)=kg(1,ipw+ikg)
       kg_k(2,ipw)=kg(2,ipw+ikg)
       kg_k(3,ipw)=kg(3,ipw+ikg)
     end do

     call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

     if (psps%useylm==1) then
!$OMP PARALLEL DO 
       do ilm=1,mpsang*mpsang
         do ipw=1,npw_k
           ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
         end do
!$OMP PARALLEL DO 
         do ii=1,9
           do ipw=1,npw_k
             ylmgr_k(ipw,ii,ilm)=ylmgr(ipw+ikg,ii,ilm)
           end do
         end do
       end do
     end if

     index=1+icg

!    Compute nonlocal psp energy

!    Compute (k+G) vectors (only if useylm=1)
     nkpg=3*nloalg(5)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

!    Compute nonlocal form factors ffnl at all (k+G):
     ider=2;idir=0;dimffnl=3+7*psps%useylm
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,&
&     ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

     enlk=0.0_dp
     eltfrnlk(:,:)=0.0_dp

!    Allocate the arrays phkxred and ph3d, compute phkxred and eventually ph3d.
     ABI_ALLOCATE(phkxred,(2,natom))
     do ia=1,natom
       iatom=atindx(ia)
       arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
       phkxred(1,iatom)=cos(arg)
       phkxred(2,iatom)=sin(arg)
     end do
     if(nloalg(1)<=0)then
!      Only the allocation, not the precomputation.
       matblk=nloalg(4)
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
     else
!      Here, allocation as well as precomputation
       matblk=natom
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
       call ph1d3d(1,natom,kg_k,matblk,natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)
     end if

     nnlout=6*(3*natom+6)
     ABI_ALLOCATE(enlout,(nnlout))

     do iband=1,nband_k

       if(mpi_enreg%proc_distrb(ikpt, iband,isppol) /= me) then
         cycle
       end if

       cwavef(:,1:npw_k*nspinor) = cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)

       signs=1 ; choice=6 ; idir=0 ; tim_nonlop=6 ; paw_opt=0 ; cpopt=-1
       call nonlop(atindx1,choice,cpopt,cprj_dum,dimekb1,dimekb2,&
&       dimffnl,dimffnl,ekb,enlout,&
&       ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&       kpoint,dum,psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,nattyp,1,&
&       ngfft,nkpg,nkpg,nloalg,nnlout,npw_k,npw_k,nspinor,nspinor,ntypat,0,paw_opt,phkxred,&
&       phkxred,ph1d,ph3d,ph3d,signs,&
&       dum_sij,dum_svectout,tim_nonlop,ucvol,psps%useylm,cwavef,cwavef)

       eltfrnlk(:,:)=eltfrnlk(:,:)+ &
&       occ(iband+bdtot_index)* reshape(enlout(:), (/6+3*natom,6/) )
     end do

     ABI_DEALLOCATE(enlout)

     eltfrnl(:,:)=eltfrnl(:,:)+wtk(ikpt)*eltfrnlk(:,:)

     ABI_DEALLOCATE(gbound)

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
!      Handle case in which kg, cg, are kept in core
       icg=icg+npw_k*nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(phkxred)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)

   end do
 end do ! End loops on isppol and ikpt

!Fill in lower triangle
 do jj=2,6
   do ii=1,jj-1
     eltfrnl(jj,ii)=eltfrnl(ii,jj)
   end do
 end do

!Accumulate eltfrnl on all proc.
 call timab(48,1,tsec)
 call xmpi_sum(eltfrnl,spaceComm,ierr)
 call timab(48,2,tsec)

!The indexing array atindx is used to reestablish the correct order of atoms
 ABI_ALLOCATE(elt_work,(6+3*natom,6))
 elt_work(1:6,1:6)=eltfrnl(1:6,1:6)
 do ia=1,natom
   ielt=7+3*(ia-1)
   ieltx=7+3*(atindx(ia)-1)
   elt_work(ielt:ielt+2,1:6)=eltfrnl(ieltx:ieltx+2,1:6)
 end do
 eltfrnl(:,:)=elt_work(:,:)

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(eltfrnlk)
 ABI_DEALLOCATE(elt_work)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(ekb)

!DEBUG
!write(std_out,*)' dyfnl3 : exit '
!stop
!ENDDEBUG

end subroutine eltfrnl3
!!***
