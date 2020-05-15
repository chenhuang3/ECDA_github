!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_bandfft_tabs
!! NAME
!! prep_bandfft_tabs
!!
!! FUNCTION
!! This routine transpose various tabs needed in bandfft parallelization
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  kpoint(3)=k point in terms of recip. translations
!!  ikpt=number of the k-point
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
!!  npw_k=number of plane waves at this k point
!!  ntypat=number of types of atoms in unit cell.
!!  option=1 if gbound have to be updated
!!         0 otherwise
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!
!! OUTPUT
!!  if option=1
!!     gbound(2*mgfft+8,2)=G sphere boundary
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! PARENTS
!!      forstrnps,vtorho
!!
!! CHILDREN
!!      bandfft_kpt_init2,mkkpg,timab,xmpi_allgatherv
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine prep_bandfft_tabs(dimffnl,ffnl,gbound,ikpt,kinpw,kpoint,lmnmax,&
& matblk,mgfft,mkmem,mpi_enreg,nkpg,npw_k,ntypat,option,ph3d)


 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_bandfft_kpt

 use m_hamiltonian, only : gs_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_bandfft_tabs'
 use interfaces_18_timing
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments -------------------------------
 integer, intent(in)    :: dimffnl,ikpt,lmnmax,matblk,mgfft,mkmem,nkpg,npw_k,ntypat,option
 integer, intent(inout) :: gbound(2*mgfft+8,2*option)
 real(dp), intent(in)   :: kpoint(3)
 real(dp), intent(in)   :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
 real(dp), intent(in)   :: kinpw(npw_k*option)
 real(dp), intent(in)   :: ph3d(2,npw_k,matblk)
 type(MPI_type), intent(inout) :: mpi_enreg
!Local variables-------------------------------
 integer  :: ierr,ikpt_this_proc,ipw,ndatarecv,spaceComm
 real(dp) :: tsec(2)
 character(len=500)   :: message
 integer, allocatable :: recvcounts(:),rdispls(:)
 integer, allocatable :: recvcountsloc(:),rdisplsloc(:)
 real(dp),allocatable :: ffnl_gather(:,:,:,:),ffnl_little(:,:,:,:),ffnl_little_gather(:,:,:,:)
 real(dp),allocatable :: kinpw_gather(:),kpg_k_gather(:,:)
 real(dp),allocatable :: ph3d_gather(:,:,:),ph3d_little(:,:,:),ph3d_little_gather(:,:,:)

! *********************************************************************

 call timab(575,1,tsec)

 ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
 if ((bandfft_kpt(ikpt_this_proc)%flag1_is_allocated==1).and.&
& (ikpt_this_proc <= mkmem).and.(ikpt_this_proc/=0)) then
   spaceComm          =mpi_enreg%comm_band
   ndatarecv          =bandfft_kpt(ikpt_this_proc)%ndatarecv
   ABI_ALLOCATE(rdispls       ,(mpi_enreg%nproc_band))
   ABI_ALLOCATE(recvcounts    ,(mpi_enreg%nproc_band))
   if (option==1) then
     gbound(:,:)       =bandfft_kpt(ikpt_this_proc)%gbound(:,:)
   end if
   rdispls(:)         =bandfft_kpt(ikpt_this_proc)%rdispls(:)
   recvcounts(:)      =bandfft_kpt(ikpt_this_proc)%recvcounts(:)
 else
   message = ' the bandfft tabs are not allocated !'
   MSG_BUG(message)
 end if

 ABI_ALLOCATE(ffnl_gather,(ndatarecv,dimffnl,lmnmax,ntypat))
 if (option==1)  then
   ABI_ALLOCATE(kinpw_gather,(ndatarecv))
 else
   ABI_ALLOCATE(kinpw_gather,(0))
 end if
 ABI_ALLOCATE(kpg_k_gather,(ndatarecv,nkpg))
 ABI_ALLOCATE(ph3d_gather,(2,ndatarecv,matblk))
 ABI_ALLOCATE(rdisplsloc    ,(mpi_enreg%nproc_band))
 ABI_ALLOCATE(recvcountsloc ,(mpi_enreg%nproc_band))

 ABI_ALLOCATE(ffnl_little,(dimffnl,lmnmax,ntypat,npw_k))
 ABI_ALLOCATE(ffnl_little_gather,(dimffnl,lmnmax,ntypat,ndatarecv))
 do ipw=1,npw_k
   ffnl_little(:,:,:,ipw)=ffnl(ipw,:,:,:)
 end do
 recvcountsloc(:)=recvcounts(:)*dimffnl*lmnmax*ntypat
 rdisplsloc(:)=rdispls(:)*dimffnl*lmnmax*ntypat
 call xmpi_allgatherv(ffnl_little,npw_k*dimffnl*lmnmax*ntypat,ffnl_little_gather,&
& recvcountsloc(:),rdisplsloc,spaceComm,ierr)
 do ipw=1,ndatarecv
   ffnl_gather(ipw,:,:,:)=ffnl_little_gather(:,:,:,ipw)
 end do
 ABI_DEALLOCATE(ffnl_little)
 ABI_DEALLOCATE(ffnl_little_gather)

 if (option==1) then
   recvcountsloc(:)=recvcounts(:)
   rdisplsloc(:)=rdispls(:)
   call xmpi_allgatherv(kinpw,npw_k,kinpw_gather,recvcountsloc(:),rdisplsloc,spaceComm,ierr)
 end if

 ABI_ALLOCATE(ph3d_little,(2,matblk,npw_k))
 ABI_ALLOCATE(ph3d_little_gather,(2,matblk,ndatarecv))
 recvcountsloc(:)=recvcounts(:)*2*matblk
 rdisplsloc(:)=rdispls(:)*2*matblk
 do ipw=1,npw_k
   ph3d_little(:,:,ipw)=ph3d(:,ipw,:)
 end do
 call xmpi_allgatherv(ph3d_little,npw_k*2*matblk,ph3d_little_gather,recvcountsloc(:),rdisplsloc,spaceComm,ierr)
 do ipw=1,ndatarecv
   ph3d_gather(:,ipw,:)=ph3d_little_gather(:,:,ipw)
 end do
 ABI_DEALLOCATE(ph3d_little_gather)
 ABI_DEALLOCATE(ph3d_little)

 if (nkpg>0) then
   call mkkpg(bandfft_kpt(ikpt_this_proc)%kg_k_gather,kpg_k_gather,kpoint,nkpg,ndatarecv)
!  recvcountsloc(:)=recvcounts(:)*nkpg
!  rdisplsloc(:)=rdispls(:)*nkpg
!  call xmpi_allgatherv(kpg_k,npw_k*nkpg,kpg_k_gather,recvcountsloc(:),rdisplsloc,spaceComm,ierr)
 end if

 ABI_DEALLOCATE(recvcounts)
 ABI_DEALLOCATE(rdispls)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)

 call bandfft_kpt_init2(bandfft_kpt,dimffnl,ffnl_gather,ikpt_this_proc,kinpw_gather,kpg_k_gather,lmnmax,matblk,&
& mkmem,ndatarecv,nkpg,ntypat,option,ph3d_gather)

 ABI_DEALLOCATE(ffnl_gather)
 ABI_DEALLOCATE(ph3d_gather)
 ABI_DEALLOCATE(kinpw_gather)
 ABI_DEALLOCATE(kpg_k_gather)

 call timab(575,2,tsec)

end subroutine prep_bandfft_tabs
!!***
