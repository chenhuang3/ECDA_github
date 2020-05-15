!{\src2tex{textfont=tt}}
!!****f* abinit/prep_nonlop
!! NAME
!! prep_nonlop
!!
!! FUNCTION
!! this routine prepares the data to the call of nonlop.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (FB,MT,GZ,MD,FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  choice: chooses possible output:
!!    choice=1 => a non-local energy contribution
!!          =2 => a gradient with respect to atomic position(s)
!!          =3 => a gradient with respect to strain(s)
!!          =23=> a gradient with respect to atm. pos. and strain(s)
!!          =4 => a 2nd derivative with respect to atomic pos.
!!          =24=> a gradient and 2nd derivative with respect to atomic pos.
!!          =5 => a gradient with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!          =7 => no operator, just projections
!!  blocksize= size of block for FFT
!!  cpopt=flag defining the status of cwaveprj=<Proj_i|Cnk> scalars (see below, side effects)
!!  cwavef(2,npw*nspinortot*blocksize)=planewave coefficients of wavefunction.
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  enl(dimenl1,dimenl2,nspinortot**2)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!  ->PAW :             ==== when paw_opt=1 or 4 ====
!!                      (Real, symmetric) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!                      These are complex numbers if cplex_enl=2
!!                        enl(:,:,1) contains Dij^up-up
!!                        enl(:,:,2) contains Dij^dn-dn
!!                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gvnlc=matrix elements <G|Vnonlocal|C>
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kpg(npw,nkpg)= if nkpg==3 (k+G) components
!!                 if nkpg==9 [(k+G)_a].[(k+G)_b] quantities
!!  kpt(3)=k point in terms of recip. translations
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                        - strain component (1:6) in the case (choice=2,signs=2) or (choice=6,signs=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                                  or i=lmn (if useylm=1)
!!  istwf_k=option parameter that describes the storage of wfs
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mgfft=maximum size of 1d ffts
!!  mpi_enreg=informations about mpi parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)=number of atoms of each type
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpg=second dimension of array kpg
!!  nloalg(5)=governs the choice of the algorithm for nonlocal operator
!!  nnlout=dimension of enlout (when signs=1):
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nspinortot=total number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell
!!  paw_opt= define the nonlocal operator concerned with
!!  phkxred(2,natom)=phase factors exp(2 pi kpt.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  sij(dimenl1,ntypat*(paw_opt/3))=overlap matrix components (only if paw_opt=2, 3 or 4)
!!  tim_nonlop=timing code of the calling routine (can be set to 0 if not attributed)
!!  ucvol=unit cell volume (bohr^3)
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!
!! OUTPUT
!!  ==== if (signs==1) ====
!!  enlout_block(nnlout)=
!!    if paw_opt==0, 1 or 2: contribution of this block of states to the nl part of various properties
!!    if paw_opt==3:         contribution of this block of states to <c|S|c>  (where S=overlap when PAW)
!!  ==== if (signs==2) ====
!!    if paw_opt==0, 1, 2 or 4:
!!       gvnlc(2,nspinor*npw)=result of the aplication of the nl operator
!!                        or one of its derivative to the input vect.
!!    if paw_opt==3 or 4:
!!       gsc(2,nspinor*npw*(paw_opt/3))=result of the aplication of (I+S)
!!                        to the input vect. (where S=overlap when PAW)
!!
!! SIDE EFFECTS
!!  ==== ONLY IF useylm=1
!!  cwaveprj(natom,nspinortot) <type(pawcprj_type)>=projected input wave function |c> on non-local projector
!!                                  =<p_lmn|c> and derivatives
!!                                  Treatment depends on cpopt parameter:
!!                     if cpopt=-1, <p_lmn|in> (and derivatives)
!!                                  have to be computed (and not saved)
!!                     if cpopt= 0, <p_lmn|in> have to be computed and saved
!!                                  derivatives are eventually computed but not saved
!!                     if cpopt= 1, <p_lmn|in> and first derivatives have to be computed and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 2  <p_lmn|in> are already in memory;
!!                                  only derivatives are computed here and not saved
!! (if useylm=0, should have cpopt=-1)
!!
!! PARENTS
!!      forstrnps,m_invovl,vtowfk
!!
!! CHILDREN
!!      dcopy,nonlop,prep_index_wavef_bandpp,prep_sort_wavef_spin,timab
!!      xmpi_allgather,xmpi_alltoallv
!!
!! NOTES
!!  cprj (as well as cg) is distributed over band processors.
!!  Only the mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) projected WFs are stored on each proc.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_nonlop(atindx1,choice,cpopt,cwaveprj,dimenl1,dimenl2,dimffnl,enl,enlout_block,&
&                       gmet,gprimd,idir,ikpt,indlmn,istwf_k,&
&                       kpt,lambdablock,lmnmax,matblk,&
&                       blocksize,mgfft,mpi_enreg,mpsang,mpssoang,&
&                       natom,nattyp,ngfft,nkpg,nloalg,nnlout,npw_k,&
&                       nspinor,nspinortot,ntypat,paw_opt,phkxred,ph1d,signs,sij,gsc,&
&                       tim_nonlop,ucvol,useylm,cwavef,gvnlc,use_gpu_cuda, already_transposed)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_bandfft_kpt

 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_nonlop'
 use interfaces_18_timing
 use interfaces_65_nonlocal
 use interfaces_66_wfs, except_this_one => prep_nonlop
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: blocksize,choice,cpopt,dimenl1,dimenl2,dimffnl,idir,ikpt,istwf_k
 integer,intent(in) :: lmnmax,matblk,mgfft,mpsang,mpssoang,signs
 integer,intent(in) :: natom,nkpg,nnlout,npw_k,nspinor,nspinortot,ntypat,paw_opt,useylm
 integer,intent(in),optional :: use_gpu_cuda
 real(dp),intent(in) :: ucvol
 type(mpi_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),nloalg(5)
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),lambdablock(blocksize)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),phkxred(2,natom)
 real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
 real(dp),intent(out) :: enlout_block(nnlout*blocksize)
 real(dp),intent(out) :: gvnlc(2,npw_k*nspinor*blocksize)
 real(dp),intent(out) :: gsc  (2,npw_k*nspinor*blocksize*(paw_opt/3))
 real(dp),intent(inout) :: cwavef(2,npw_k*nspinortot*blocksize)
 type(pawcprj_type),intent(inout) :: cwaveprj(natom,nspinortot*mpi_enreg%bandpp*((cpopt+5)/5))
 logical, optional, intent(in) :: already_transposed

!Local variables-------------------------------
 integer :: spaceComm=0
 integer :: ier,ikpt_this_proc,old_me_g0,tim_nonlop,usegpu
 real(dp) :: lambda,lambda_nonlop(mpi_enreg%bandpp)
 real(dp) :: tsec(2)
 real(dp), allocatable :: enlout(:)
 logical :: do_transpose

!local variables for mpialltoallv
 integer, pointer :: ndatarecv
 integer,ABI_CONTIGUOUS  pointer :: kg_k_gather(:,:)
 integer,ABI_CONTIGUOUS  pointer :: recvcounts(:),sendcounts(:),sdispls(:),rdispls(:)
 integer, allocatable :: sendcountsloc(:),sdisplsloc(:),recvcountsloc(:),rdisplsloc(:)
 real(dp), ABI_CONTIGUOUS pointer :: kpg_k_gather(:,:),ffnl_gather(:,:,:,:),ph3d_gather(:,:,:)
 real(dp), allocatable :: cwavef_alltoall(:,:),gvnlc_alltoall(:,:),gsc_alltoall(:,:)

!local variables for bandpp
 integer :: nproc_band,bandpp
 integer,allocatable :: index_wavef_band(:)

! *************************************************************************

 DBG_ENTER('COLL')

 call timab(570,1,tsec)

 do_transpose = .true.
 if(present(already_transposed)) then
   if(already_transposed) do_transpose = .false.
 end if

 nproc_band = mpi_enreg%nproc_band
 bandpp     = mpi_enreg%bandpp
 spaceComm=mpi_enreg%comm_fft
 if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_band
 ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)

 usegpu=0;if(present(use_gpu_cuda)) usegpu=use_gpu_cuda

 ABI_ALLOCATE(sendcountsloc,(nproc_band))
 ABI_ALLOCATE(sdisplsloc   ,(nproc_band))
 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc   ,(nproc_band))

 recvcounts   => bandfft_kpt(ikpt_this_proc)%recvcounts(:)
 sendcounts   => bandfft_kpt(ikpt_this_proc)%sendcounts(:)
 rdispls      => bandfft_kpt(ikpt_this_proc)%rdispls   (:)
 sdispls      => bandfft_kpt(ikpt_this_proc)%sdispls   (:)
 ndatarecv    => bandfft_kpt(ikpt_this_proc)%ndatarecv

 kg_k_gather  => bandfft_kpt(ikpt_this_proc)%kg_k_gather(:,:)
 kpg_k_gather => bandfft_kpt(ikpt_this_proc)%kpg_k_gather(:,:)
 ffnl_gather  => bandfft_kpt(ikpt_this_proc)%ffnl_gather(:,:,:,:)
 ph3d_gather  => bandfft_kpt(ikpt_this_proc)%ph3d_gather(:,:,:)

 ABI_ALLOCATE(cwavef_alltoall,(2,ndatarecv*nspinor*bandpp))
 ABI_ALLOCATE(gsc_alltoall,(2,ndatarecv*nspinor*(paw_opt/3)*bandpp))
 ABI_ALLOCATE(gvnlc_alltoall,(2,ndatarecv*nspinor*bandpp))

 ABI_ALLOCATE(enlout,(nnlout*bandpp))
 enlout = zero

 recvcountsloc(:)=recvcounts(:)*2*nspinor*bandpp
 rdisplsloc(:)=rdispls(:)*2*nspinor*bandpp
 sendcountsloc(:)=sendcounts(:)*2*nspinor
 sdisplsloc(:)=sdispls(:)*2*nspinor

 if(do_transpose) then
   call timab(581,1,tsec)
   call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall,&
&   recvcountsloc,rdisplsloc,spaceComm,ier)
   call timab(581,2,tsec)
 else
   ! Here, we cheat, and use DCOPY to bypass some compiler's overzealous bound-checking
   ! (ndatarecv*nspinor*bandpp might be greater than the declared size of cwavef)
   call DCOPY(2*ndatarecv*nspinor*bandpp, cwavef, 1, cwavef_alltoall, 1)
 end if

 if(istwf_k==2) then
   old_me_g0=mpi_enreg%me_g0
   if (mpi_enreg%me_fft==0) then
     mpi_enreg%me_g0=1
   else
     mpi_enreg%me_g0=0
   end if
 end if

 if (paw_opt==2) lambda=lambdablock(mpi_enreg%me_band+1)

!=====================================================================
 if (bandpp==1) then

   if (do_transpose .and. mpi_enreg%paral_spinor==0.and.nspinortot==2) then !Sort WF by spin
     call prep_sort_wavef_spin(nproc_band,nspinor,ndatarecv,recvcounts,rdispls,index_wavef_band)
     cwavef_alltoall(:, :)=cwavef_alltoall(:,index_wavef_band)
   end if

   call nonlop(atindx1,choice,cpopt,cwaveprj,dimenl1,dimenl2,dimffnl,dimffnl,&
&   enl,enlout,ffnl_gather,ffnl_gather,gmet,gprimd,idir,indlmn,&
&   istwf_k,kg_k_gather,kg_k_gather,kpg_k_gather,kpg_k_gather,kpt,kpt,(/lambda/),lmnmax,matblk,mgfft,&
&   mpi_enreg,mpsang,mpssoang,natom,nattyp,1,ngfft,nkpg,nkpg,nloalg,&
&   nnlout,ndatarecv,ndatarecv,nspinor,nspinortot,ntypat,0,paw_opt,phkxred,&
&   phkxred,ph1d,ph3d_gather,ph3d_gather,signs,sij,gsc_alltoall,&
&   tim_nonlop,ucvol,useylm,cwavef_alltoall,gvnlc_alltoall,use_gpu_cuda=usegpu)

   if (do_transpose .and. mpi_enreg%paral_spinor == 0 .and. nspinortot==2.and.signs==2) then
     if (paw_opt/=3) gvnlc_alltoall(:,index_wavef_band)=gvnlc_alltoall(:,:)
     if (paw_opt==3.or.paw_opt==4) gsc_alltoall(:,index_wavef_band)=gsc_alltoall(:,:)
   end if

 else

!  -------------------------------------------------------------
!  Computation of the index used to sort the waves functions below bandpp
!  -------------------------------------------------------------
   if(do_transpose) then
     call prep_index_wavef_bandpp(nproc_band,bandpp,&
&     nspinor,ndatarecv,recvcounts,rdispls,index_wavef_band)

!  -------------------------------------------------------
!  Sorting of the waves functions below bandpp
!  -------------------------------------------------------
     cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)
     if (signs==2) then
       if (paw_opt/=3) gvnlc_alltoall(:,:)=gvnlc_alltoall(:,index_wavef_band)
       if (paw_opt==3.or.paw_opt==4) gsc_alltoall(:,:)=gsc_alltoall(:,index_wavef_band)
     end if
   end if

!  -------------------------------------------------------
!  Call nonlop
!  -------------------------------------------------------
   if(paw_opt == 2) lambda_nonlop = lambdablock((mpi_enreg%me_band*bandpp)+1:((mpi_enreg%me_band+1)*bandpp))
   call nonlop(atindx1,choice,cpopt,cwaveprj,dimenl1,dimenl2,dimffnl,dimffnl,&
&   enl,enlout,ffnl_gather,ffnl_gather,gmet,gprimd,idir,indlmn,&
&   istwf_k,kg_k_gather,kg_k_gather,kpg_k_gather,kpg_k_gather,kpt,kpt,&
&   lambda_nonlop,lmnmax,matblk,mgfft,&
&   mpi_enreg,mpsang,mpssoang,natom,nattyp,bandpp,ngfft,nkpg,nkpg,nloalg,&
&   nnlout,ndatarecv,ndatarecv,nspinor,nspinortot,ntypat,0,paw_opt,phkxred,&
&   phkxred,ph1d,ph3d_gather,ph3d_gather,signs,sij,gsc_alltoall,&
&   tim_nonlop,ucvol,useylm,cwavef_alltoall,gvnlc_alltoall,use_gpu_cuda=usegpu)

!  -----------------------------------------------------
!  Sorting of waves functions below the processors
!  -----------------------------------------------------
   if(do_transpose) then
     cwavef_alltoall(:,index_wavef_band)=cwavef_alltoall(:,:)
     if (signs==2) then
       if (paw_opt/=3) gvnlc_alltoall(:,index_wavef_band)=gvnlc_alltoall(:,:)
       if (paw_opt==3.or.paw_opt==4) gsc_alltoall(:,index_wavef_band)=gsc_alltoall(:,:)
     end if
   end if

 end if

!=====================================================================
!  -------------------------------------------------------
!  Deallocation
!  -------------------------------------------------------
 if (allocated(index_wavef_band)) then
   ABI_DEALLOCATE(index_wavef_band)
 end if

!Transpose the gsc_alltoall or gvlnc_alltoall tabs
!according to the paw_opt and signs values
 if(do_transpose) then
   if (signs==2) then
     call timab(581,1,tsec)
     if (paw_opt/=3) then
       call xmpi_alltoallv(gvnlc_alltoall,recvcountsloc,rdisplsloc,gvnlc,&
&       sendcountsloc,sdisplsloc,spaceComm,ier)
     end if
     if (paw_opt==3.or.paw_opt==4) then
       call xmpi_alltoallv(gsc_alltoall,recvcountsloc,rdisplsloc,gsc,&
&       sendcountsloc,sdisplsloc,spaceComm,ier)
     end if
     call timab(581,2,tsec)
   end if
 else
   ! TODO check other usages, maybe
   call DCOPY(2*ndatarecv*nspinor*bandpp, gsc_alltoall, 1, gsc, 1)
 end if
 if (istwf_k==2) mpi_enreg%me_g0=old_me_g0

 if (nnlout>0) then
   call xmpi_allgather(enlout,nnlout*bandpp,enlout_block,spaceComm,ier)
 end if
 ABI_DEALLOCATE(enlout)
 ABI_DEALLOCATE(sendcountsloc)
 ABI_DEALLOCATE(sdisplsloc)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)
 ABI_DEALLOCATE(cwavef_alltoall)
 ABI_DEALLOCATE(gvnlc_alltoall)
 ABI_DEALLOCATE(gsc_alltoall)

 call timab(570,2,tsec)

 DBG_EXIT('COLL')

end subroutine prep_nonlop
!!***
