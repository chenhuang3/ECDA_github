!{\src2tex{textfont=tt}}
!!****f* ABINIT/getgsc
!!
!! NAME
!! getgsc
!!
!! FUNCTION
!! Compute <G|S|C> for all input vectors |Cnk> at a given k-point,
!!              OR for one input vector |Cnk>.
!! |Cnk> are expressed in reciprocal space.
!! S is the overlap operator between |Cnk> (used for PAW).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  cprj(natom,mcprj)= wave functions projected with non-local projectors: cprj=<p_i|Cnk>
!!  dimcprj(natom)=array of dimensions of arrays cprj, cprjq (ordered by atom-type)
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  ibg=shift to be applied on the location of data in the array cprj (beginning of current k-point)
!!  icg=shift to be applied on the location of data in the array cg (beginning of current k-point)
!!  igsc=shift to be applied on the location of data in the array gsc (beginning of current k-point)
!!  ikpt,isppol=indexes of current (spin.kpoint)
!!  kg_k(3,npw_k)=G vec coordinates wrt recip lattice transl.
!!  mcg=second dimension of the cg array
!!  mcprj=second dimension of the cprj array
!!  mgsc=second dimension of the gsc array
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in unit cell.
!!  nband= if positive: number of bands at this k point for that spin polarization
!!         if negative: abs(nband) is the index of the only band to be computed
!!  npw_k=number of planewaves in basis for given k point.
!!  nspinor=number of spinorial components of the wavefunctions
!!  ph3d(2,npw_k,matblk)=3-dim structure factors, for each atom and plane wave.
!!
!! OUTPUT
!!  gsc(2,mgsc)= <g|S|Cnk> or <g|S^(1)|Cnk> (S=overlap)
!!
!! PARENTS
!!      vtowfk3
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_destroy,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getgsc(cg,cprj,dimcprj,dimffnl,ffnl,gs_ham,gsc,ibg,icg,igsc,ikpt,isppol,&
&          kg_k,mcg,mcprj,mgsc,mpi_enreg,natom,nband,npw_k,nspinor,ph3d)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors

 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_destroy, pawcprj_copy

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getgsc'
 use interfaces_18_timing
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: dimffnl,ibg,icg,igsc,ikpt,isppol,mcg,mcprj
 integer,intent(in) :: mgsc,natom,nband,npw_k,nspinor
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_ham
!arrays
 integer,intent(in) :: dimcprj(natom),kg_k(3,npw_k)
 real(dp),intent(in) :: cg(2,mcg),ffnl(npw_k,dimffnl,gs_ham%lmnmax,gs_ham%ntypat)
 real(dp),intent(inout) :: ph3d(2,npw_k,gs_ham%matblk)
 real(dp),intent(out) :: gsc(2,mgsc)
 type(pawcprj_type),intent(in) :: cprj(natom,mcprj)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimenl1,dimenl2,iband,iband1,iband2,ierr,index_cg,index_cprj
 integer :: index_gsc,me,my_nspinor,nkpg,paw_opt,signs,tim_nonlop,usecprj,useylm
 character(len=500) :: msg
!arrays
 real(dp) :: rdum(1)
 real(dp) :: dummy1(1),tsec(2)
 real(dp),allocatable :: cwavef(:,:),enl_dum(:,:,:),kpg_dum(:,:),scwavef(:,:)
 real(dp) :: vect_dum(0,0)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
 if(gs_ham%usepaw==0) then
   msg='Only compatible with PAW (usepaw=1) !'
   MSG_BUG(msg)
 end if
 if(nband<0.and.(mcg<npw_k*my_nspinor.or.mgsc<npw_k*my_nspinor.or.mcprj<my_nspinor)) then
   msg='Invalid value for mcg, mgsc or mcprj !'
   MSG_BUG(msg)
 end if

!Keep track of total time spent in getgsc:
 call timab(565,1,tsec)

!Prepare some data
 ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(scwavef,(2,npw_k*my_nspinor))
 if (mcprj>0) then
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,my_nspinor))
   call pawcprj_alloc(cwaveprj,0,dimcprj)
 end if
 usecprj=0;if (mcprj>0) usecprj=1
 dimenl1=gs_ham%dimekb1;dimenl2=natom;tim_nonlop=0;nkpg=0
 choice=1;signs=2;cpopt=-1+3*usecprj;paw_opt=3;useylm=1
 ABI_ALLOCATE(enl_dum,(dimenl1,natom,nspinor**2))
 enl_dum=zero
 ABI_ALLOCATE(kpg_dum,(3,nkpg))
 me=mpi_enreg%me_kpt

!Loop over bands
 index_cprj=ibg;index_cg=icg;index_gsc=igsc
 if (nband>0) then
   iband1=1;iband2=nband
 else if (nband<0) then
   iband1=abs(nband);iband2=iband1
   index_cprj=index_cprj+(iband1-1)*my_nspinor
   index_cg  =index_cg  +(iband1-1)*npw_k*my_nspinor
   index_gsc =index_gsc +(iband1-1)*npw_k*my_nspinor
 end if

 do iband=iband1,iband2

   if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me.and.nband>0) then
     gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=zero
     index_cprj=index_cprj+my_nspinor
     index_cg=index_cg+npw_k*my_nspinor
     index_gsc=index_gsc+npw_k*my_nspinor
     cycle
   end if
   
!  Retrieve WF at (n,k)
   cwavef(:,1:npw_k*my_nspinor)=cg(:,1+index_cg:npw_k*my_nspinor+index_cg)
   if (usecprj==1) then
     call pawcprj_copy(cprj(:,1+index_cprj:my_nspinor+index_cprj),cwaveprj)
   end if

!  Compute <g|S|Cnk>
   call nonlop(gs_ham%atindx1,choice,cpopt,cwaveprj,dimenl1,dimenl2,dimffnl,dimffnl,&
&   enl_dum,dummy1,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,0,gs_ham%indlmn,&
&   gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&   rdum,gs_ham%lmnmax,gs_ham%matblk,gs_ham%mgfft,mpi_enreg,gs_ham%mpsang,gs_ham%mpssoang,&
&   natom,gs_ham%nattyp,1,gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,1,npw_k,npw_k,&
&   my_nspinor,nspinor,gs_ham%ntypat,0,paw_opt,gs_ham%phkxred,gs_ham%phkxred,&
&   gs_ham%ph1d,ph3d,ph3d,signs,gs_ham%sij,scwavef,tim_nonlop,gs_ham%ucvol,&
&   useylm,cwavef,vect_dum,use_gpu_cuda=gs_ham%use_gpu_cuda)

   gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=scwavef(:,1:npw_k*my_nspinor)

!  End of loop over bands
   index_cprj=index_cprj+my_nspinor
   index_cg=index_cg+npw_k*my_nspinor
   index_gsc=index_gsc+npw_k*my_nspinor
 end do

!Reduction in case of parallelization
 if ((xmpi_paral==1)) then
   call timab(48,1,tsec)
   call xmpi_sum(gsc,mpi_enreg%comm_band,ierr)
   call timab(48,2,tsec)
 end if

!Memory deallocation
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(scwavef)
 ABI_DEALLOCATE(enl_dum)
 ABI_DEALLOCATE(kpg_dum)
 if (usecprj==1) then
   call pawcprj_destroy(cwaveprj)
   ABI_DATATYPE_DEALLOCATE(cwaveprj)
 end if

 call timab(565,2,tsec)

 DBG_EXIT("COLL")

end subroutine getgsc
!!***
