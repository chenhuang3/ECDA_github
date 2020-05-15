!{\src2tex{textfont=tt}}
!!****f* ABINIT/store_bfield_cprj
!! NAME
!! store_bfield_cprj
!!
!! FUNCTION
!! This routine stores cprj in bfield structure
!!
!! COPYRIGHT
!! Copyright (C) 2003-2014 ABINIT  group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! cprj(natom,mcprj*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!! ikpt=index of ikpt to store (0 for all)
!! isppol_index=index of spin polarization to store (0 for all)
!! mband=maximum number of bands
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mkmem=number of k points treated by this node.
!! natom=number of atoms in cell
!! nkpt=number of k-points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! Input/Output
!! dtbfield <type(bfield_type)> = variables related to magnetization
!! mpi_enreg=informations about MPI parallelization
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      update_eb_field_vars,update_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_destroy,pawcprj_get,pawcprj_mpi_allgather
!!      pawcprj_put,xmpi_allgather
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine store_bfield_cprj(atindx1,cprj,dtbfield,ikpt,isppol_index,&
&                            mband,mcprj,mkmem,mpi_enreg,natom,nkpt,nsppol)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_bfield
 use m_profiling_abi

 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_mpi_allgather, pawcprj_put, pawcprj_destroy

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'store_bfield_cprj'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: ikpt,isppol_index,mband,mcprj,mkmem,natom,nkpt,nsppol
 type(MPI_type), intent(inout) :: mpi_enreg
 type(bfield_type), intent(inout) :: dtbfield
!arrays
 integer, intent(in) :: atindx1(natom)
 type(pawcprj_type),intent(in) :: cprj(natom,mcprj)

!Local variables -------------------------
! scalars
 integer :: icp1,icp2,ierr,ikpt_loc,ikpt1,iproc,isppol,me,my_nspinor
 integer :: n2dim,nband_k,ncpgr,nproc,ntotcp,spaceComm
 character(len=500) :: message
! arrays
 integer,allocatable :: dimlmn(:),ikpt1_recv(:)
 type(pawcprj_type),allocatable :: cprj_gat(:,:),cprj_k(:,:)

! ***********************************************************************

!Init MPI
 spaceComm=mpi_enreg%comm_cell
 nproc=xcomm_size(spaceComm)
 me=mpi_enreg%me_kpt

 my_nspinor=max(1,dtbfield%nspinor/mpi_enreg%nproc_spinor)

 if (dtbfield%usecprj /= 1) then
   message = ' store_bfield_cprj: cprj datastructure has not been allocated !'
   MSG_BUG(message)
 end if

 ABI_ALLOCATE(dimlmn,(natom))
 dimlmn(:)=cprj(1:natom,1)%nlmn
 ncpgr = cprj(1,1)%ncpgr

 ABI_DATATYPE_ALLOCATE(cprj_k,(natom,dtbfield%nspinor*mband))
 ABI_DATATYPE_ALLOCATE(cprj_gat,(natom,nproc*dtbfield%nspinor*mband))
 call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
 call pawcprj_alloc(cprj_gat,ncpgr,dimlmn)

 n2dim = dtbfield%nspinor*mband
 ntotcp = n2dim*SUM(dimlmn(1:natom))

 do isppol = 1, nsppol

   if(isppol_index > 0 .and. (isppol /= isppol_index) ) cycle

   ikpt_loc = 0
   ikpt1 = 0
   do while (ikpt_loc < mkmem)

     if (ikpt_loc < mkmem) ikpt1 = ikpt1 + 1
     if ((ikpt1 > nkpt).and.(ikpt_loc < mkmem)) exit
     nband_k = dtbfield%nband_occ

     if ( (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt1,1,nband_k,isppol,me)) .and. &
&     (ikpt_loc <= mkmem) ) cycle

     ikpt_loc = ikpt_loc + 1
     
     if(ikpt>0 .and. (ikpt1 /= ikpt)) cycle ! if a specific kpt was requested but we are not at it yet, move on

     ABI_ALLOCATE(ikpt1_recv,(nproc))
     call xmpi_allgather(ikpt1,ikpt1_recv,spaceComm,ierr)
     call pawcprj_get(atindx1,cprj_k,cprj,natom,1,(ikpt_loc-1)*nband_k*my_nspinor,ikpt1,0,isppol,mband,&
&     mkmem,natom,nband_k,nband_k,my_nspinor,nsppol,0,&
&     mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
     call pawcprj_mpi_allgather(cprj_k,cprj_gat,natom,n2dim,dimlmn,ncpgr,nproc,spaceComm,ierr,rank_ordered=.true.)
     do iproc = 1, nproc
       icp2=nband_k*(iproc-1)*my_nspinor
       call pawcprj_get(atindx1,cprj_k,cprj_gat,natom,1,icp2,ikpt1,0,isppol,mband,&
&       nproc,natom,nband_k,nband_k,my_nspinor,1,0,&
&       mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       icp1 = nband_k*(ikpt1_recv(iproc)-1)*my_nspinor
       call pawcprj_put(atindx1,cprj_k,dtbfield%cprj,natom,1,icp1,ikpt1,0,isppol,&
&       mband,dtbfield%fnkpt,natom,nband_k,nband_k,dimlmn,my_nspinor,nsppol,0,&
&       mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
     end do
     ABI_DEALLOCATE(ikpt1_recv)

   end do ! close loop over k-points

 end do ! end loop over nsppol

 ABI_DEALLOCATE(dimlmn)

 call pawcprj_destroy(cprj_k)
 call pawcprj_destroy(cprj_gat)
 ABI_DATATYPE_DEALLOCATE(cprj_k)
 ABI_DATATYPE_DEALLOCATE(cprj_gat)


!DEBUG
!write(std_out,*)'store_bfield_cprj EXIT'
!END_DEBUG
end subroutine store_bfield_cprj
!!***
