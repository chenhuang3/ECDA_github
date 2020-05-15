!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkresi
!! NAME
!! mkresi
!!
!! FUNCTION
!! Make residuals from knowledge of wf in G space and application of Hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mcg)=<G|Cnk>=Fourier coefficients of wavefunction
!!  dimffnl=second dimension of ffnl (1+number of derivatives)g
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  filstat=name for the status file
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  kg_k(3,npw)=planewave reduced coordinates in basis sphere.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  mcg=second dimension of the cg array
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in unit cell.
!!  nband=number of bands involved in subspace matrix.
!!  npw=number of planewaves in basis sphere at this k point.
!!  nspinor=number of spinors (on current proc)
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vlocal(n4,n5,n6,nvloc)=local potential in real space, on the augmented fft grid
!!
!! OUTPUT
!!  eig_k(nband)$= \langle C_n \mid H \mid C_n \rangle $ for each band.
!!  resid_k(nband)=residual for each band
!!   $= \langle C_n \mid H H \mid C_n \rangle- \langle C_n \mid H \mid C_n \rangle^2 $.
!!
!! PARENTS
!!      energy
!!
!! CHILDREN
!!      dotprod_g,getghc,sqnorm_g,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkresi(cg,dimffnl,eig_k,ffnl,filstat,gs_hamk,icg,kg_k,kinpw,&
&                 mcg,mpi_enreg,natom,nband,npw,nspinor,paral_kgb,ph3d,prtvol,&
&                 resid_k,usepaw,vlocal)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_cgtools

 use m_pawcprj,     only : pawcprj_type
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_fock,        only : fock_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkresi'
 use interfaces_18_timing
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimffnl,icg,mcg,natom,nband,npw,nspinor
 integer,intent(in) :: paral_kgb,prtvol,usepaw
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_hamk
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(in) :: cg(2,mcg),ffnl(npw,dimffnl,gs_hamk%lmnmax,gs_hamk%ntypat),kinpw(npw)
 real(dp),intent(inout) :: ph3d(2,npw,gs_hamk%matblk),vlocal(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc)
 real(dp),intent(out) :: eig_k(nband),resid_k(nband)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_getghc=3
 integer :: cpopt,iband,ipw,istwf_k
 real(dp) :: doti,dotr
 type(fock_type),pointer :: fock => null()
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
 type(pawcprj_type) :: cwaveprj(1,1)

! *************************************************************************

!Keep track of total time spent in mkresi
 call timab(13,1,tsec)

 istwf_k=gs_hamk%istwf_k

 do iband=1,nband
   ABI_ALLOCATE(cwavef,(2,npw*nspinor))
   ABI_ALLOCATE(ghc,(2,npw*nspinor))
   ABI_ALLOCATE(gvnlc,(2,npw*nspinor))
   if (usepaw==1)  then
     ABI_ALLOCATE(gsc,(2,npw*nspinor))
   else
     ABI_ALLOCATE(gsc,(0,0))
   end if

!$OMP PARALLEL DO 
   do ipw=1,npw*nspinor
     cwavef(1,ipw)=cg(1,ipw+(iband-1)*npw*nspinor+icg)
     cwavef(2,ipw)=cg(2,ipw+(iband-1)*npw*nspinor+icg)
   end do

   cpopt=-1
   call getghc(cpopt,cwavef,cwaveprj,dimffnl,ffnl,filstat,ghc,gsc,&
&   gs_hamk,gvnlc,kg_k,kinpw,dotr,mpi_enreg,natom,1,npw,nspinor,&
&   paral_kgb,ph3d,prtvol,usepaw,tim_getghc,0,vlocal,fock)
   ABI_DEALLOCATE(gvnlc)

!  Compute the residual, <Cn|(H-<Cn|H|Cn>)**2|Cn>:
!  First get eigenvalue <Cn|H|Cn>:
   call dotprod_g(dotr,doti,istwf_k,npw*nspinor,1,cwavef,ghc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   eig_k(iband)=dotr

!  Next need <G|(H-<Cn|H|Cn>)|Cn> (in ghc):
!  ghc(:,:)=ghc(:,:)-eig_k(iband)*cwavef(:,:)
   if (usepaw==0) then
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(cwavef,eig_k,ghc,iband,npw,nspinor)
     do ipw=1,npw*nspinor
       ghc(1,ipw)=ghc(1,ipw)-eig_k(iband)*cwavef(1,ipw)
       ghc(2,ipw)=ghc(2,ipw)-eig_k(iband)*cwavef(2,ipw)
     end do
   else
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(gsc,eig_k,ghc,iband,npw,nspinor)
     do ipw=1,npw*nspinor
       ghc(1,ipw)=ghc(1,ipw)-eig_k(iband)*gsc(1,ipw)
       ghc(2,ipw)=ghc(2,ipw)-eig_k(iband)*gsc(2,ipw)
     end do
   end if

!  Then simply square the result:
   call sqnorm_g(dotr,istwf_k,npw*nspinor,ghc,mpi_enreg%me_g0,mpi_enreg%comm_fft)
   resid_k(iband)=dotr

   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(ghc)
   ABI_DEALLOCATE(gsc)

 end do

 call timab(13,2,tsec)

end subroutine mkresi
!!***
