!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_write_cg
!! NAME
!! dfpt_write_cg
!!
!! FUNCTION
!! Conduct output of a "wave-functions" file.
!!  - Compute the maximal residual
!!  - Then open a permanent file wff2 for final output of wf data
!!  - Create a new header for the file.
!!  - Write wave-functions (and energies)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mcg)=wavefunction array (storage if nkpt>1)
!!  eigen((2*mband)**response *mband*nkpt*nsppol)= eigenvalues (hartree) for all bands at each k point
!!  fname= character string giving the root to form the name of the
!!   output WFK or WFQ file if response==0, otherwise it is the filename.
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem=number of k-points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum number of plane waves
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outwf
!!
!! CHILDREN
!!      cg_zcopy,cwtime,flush_unit,m_header_init,mask2blocks,timab
!!      unpack_eneocc,wfk_close,wfk_open_write,wfk_write_band_block,wrtout
!!      xmpi_barrier,xmpi_exch
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_write_cg(fname,Hdr,mpw,mband,nband,nkpt,nsppol,nspinor,mcg,mkmem,eigen,cg,kg,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_wfk
 use m_header

 use m_cgtools,        only : cg_zcopy
 use m_io_tools,       only : iomode_from_fname, get_unit
 use m_numeric_tools,  only : mask2blocks
 use m_cgtools,        only : cg_zcopy
 use m_time,           only : cwtime

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_write_cg'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
!End of the abilint section

 implicit none


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mkmem,mpw,nkpt,nsppol,nspinor
 character(len=*),intent(in) :: fname
 type(MPI_type),intent(in) :: mpi_enreg
 type(hdr_type),intent(in) :: hdr
!arrays
 integer, intent(in) :: nband(nkpt*nsppol),kg(3,mpw*mkmem)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(2*mband*mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig1=1,tim_rwwf12=12,master=0
 integer :: iomode,icg,ikg,ikpt,spin,my_rank,my_nspinor,nband_k,npw_k,comm,blk 
 integer :: action,nblocks,ierr,nprocs,source !,comm_band
 logical :: ihave_data,iam_master,single_writer
 character(len=500) :: msg
 real(dp) :: cpu,wall,gflops
 type(wfk_t) :: Wfk
!arrays
 integer :: band_block(2)
 integer,allocatable :: kg_k(:,:)
 integer,pointer :: blocks(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cg_k(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 ABI_UNUSED(nband(1))

!FIXME
 call m_header_init(ierr)
 ABI_CHECK(ierr==0,"header_init returned ierr!=0")

 call timab(270+tim_rwwf12,1,tsec)

 if (mpi_enreg%paralbd==1) then
 end if

!Init comm and my_rank
 comm       = mpi_enreg%comm_cell
 nprocs     = xcomm_size(comm)
 my_rank    = mpi_enreg%me_kpt
 iam_master = (master==my_rank)
!comm_band = mpi_enreg%comm_band

 my_nspinor=MAX(1,nspinor/mpi_enreg%nproc_spinor)

 iomode = iomode_from_fname(fname)

 single_writer = (ANY(iomode == (/IO_MODE_FORTRAN, IO_MODE_ETSF/) ))

 write(msg,'(3a,i0)')ABI_FUNC//': writing DFPT WFK file ',trim(fname),", with iomode ",iomode
 call wrtout(std_out,msg,'PERS')
 write(std_out,*)"single_writer ",single_writer

 call cwtime(cpu,wall,gflops,"start")

 if (.not.single_writer) then

   if (my_rank==master) then
     call wfk_open_write(Wfk,Hdr,fname,formeig1,iomode,get_unit(),xmpi_self,write_hdr=.TRUE.,write_frm=.TRUE.)
     call wfk_write_allgkk(Wfk,xmpio_single,eigen)
   end if
   call xmpi_barrier(comm)

   if (my_rank/=master) then
     call wfk_open_write(Wfk,Hdr,fname,formeig1,iomode,get_unit(),xmpi_self,write_hdr=.FALSE.,write_frm=.FALSE.)
   end if

   call cwtime(cpu,wall,gflops,"stop")
   write(msg,'(2(a,f8.2))')" open_write, cpu: ",cpu,", wall: ",wall
   call wrtout(std_out,msg,"PERS")

   call cwtime(cpu,wall,gflops,"start")

   call cwtime(cpu,wall,gflops,"stop")
   write(msg,'(2(a,f8.2))')" write_allgkk cpu: ",cpu,", wall: ",wall
   call wrtout(std_out,msg,"PERS")

   call cwtime(cpu,wall,gflops,"start")

   icg=0
   do spin=1,nsppol
     ikg=0
     do ikpt=1,nkpt
       nband_k = Wfk%nband(ikpt,spin)
       npw_k   = Hdr%npwarr(ikpt)

!      Compute my block of bands for this k-point and spin.
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,my_rank)) CYCLE

       call mask2blocks(mpi_enreg%proc_distrb(ikpt,:,spin)==my_rank, nblocks,blocks)
       ABI_CHECK(nblocks==1,"nblocks !=1")

       ! write(msg,"(a,3(i0,2x))")"Will write (ikpt, spin, nblocks)",ikpt,spin,nblocks
       ! call wrtout(std_out,msg,"PERS")

       do blk=1,nblocks
         band_block = blocks(:,blk)

         if (band_block(1)==1) then
           ! This processor writes kg_k
           call wfk_write_band_block(Wfk,band_block,ikpt,spin,xmpio_single,kg_k=kg(:,1+ikg:),cg_k=cg(:,1+icg:))
         else 
           call wfk_write_band_block(Wfk,band_block,ikpt,spin,xmpio_single,cg_k=cg(:,1+icg:))
         end if
       end do

       ABI_FREE(blocks)

       icg = icg+npw_k*my_nspinor*nband_k
       ikg = ikg+npw_k
     end do 
   end do 

   call wfk_close(Wfk)
   call xmpi_barrier(comm)

 else
   MSG_ERROR("Not coded yet")

   ABI_MALLOC(kg_k,(3,mpw))
   ABI_MALLOC(cg_k,(2,mpw*my_nspinor*mband))
   ABI_CHECK_ALLOC("out of memory in cg_k")

   if (iam_master) then
     call wfk_open_write(Wfk,Hdr,fname,formeig1,iomode,get_unit(),xmpi_self,write_hdr=.TRUE., write_frm=.FALSE.)
   end if

   icg=0
   do spin=1,nsppol
     ikg=0
     do ikpt=1,nkpt

       nband_k = nband(ikpt+(spin-1)*nkpt)
       npw_k   = Hdr%npwarr(ikpt)

       call xmpi_barrier(comm)

       ! Transfer the wavefunctions and the g-vectors to the master processor
       source = MINVAL(mpi_enreg%proc_distrb(ikpt,1:nband_k,spin))
       ihave_data = (source==my_rank)


       action=0
       if (iam_master .and. ihave_data)    action=1 ! I am the master node, and I have the data in cg
       if (.not.iam_master.and.ihave_data) action=2 ! I am not the master, and I have the data => send to master
       if (iam_master.and..not.ihave_data) action=3 ! I am the master, and I receive the data

       if (action==1) then ! Copy from kg and cg
         kg_k(:,1:npw_k) = kg(:,ikg+1:ikg+npw_k)
         call cg_zcopy(npw_k*my_nspinor*nband_k, cg(1,icg+1), cg_k)
       end if

!      Exchange data
       if (action==2.or.action==3) then
         call timab(48,1,tsec)
         if (action==2) then
           call xmpi_exch(kg(:,1+ikg:npw_k+ikg),3*npw_k,source,kg_k,master,comm,ierr)
           call xmpi_exch(cg(:,icg+1:icg+nband_k*npw_k*my_nspinor),2*nband_k*npw_k*my_nspinor,&
&           source,cg_k,master,comm,ierr)
         else
           call xmpi_exch(kg_k,3*npw_k,source,kg_k,master,comm,ierr)
           call xmpi_exch(cg_k,2*nband_k*npw_k*my_nspinor,source,cg_k,master,comm,ierr)
         end if
         call timab(48,2,tsec)
       end if


!      Master writes this block of bands.
       if (iam_master) then
         band_block = (/1,nband_k/)
!        call wfk_write_band_block(Wfk,band_block,ikpt,spin,xmpio_single,kg_k=kg_k,cg_k=cg_k,&
!        &          eig_k=eigen3d(:,ikpt,spin))
       end if

       if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,my_rank))) then
         icg = icg+npw_k*my_nspinor*nband_k
         ikg = ikg+npw_k
       end if

     end do !ikpt
   end do !spin

   ABI_FREE(kg_k)
   ABI_FREE(cg_k)

   if (iam_master) then
     call wfk_close(Wfk)
   end if
   call xmpi_barrier(comm)

 end if

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))')" write all cg cpu: ",cpu,", wall: ",wall
 call wrtout(std_out,msg,"PERS")

 call timab(270+tim_rwwf12,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_write_cg
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/gs_write_cg
!! NAME
!! gs_write_cg
!!
!! FUNCTION
!! Conduct output of a "wave-functions" file.
!!  - Compute the maximal residual
!!  - Then open a permanent file wff2 for final output of wf data
!!  - Create a new header for the file.
!!  - Write wave-functions (and energies)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mcg)=wavefunction array (storage if nkpt>1)
!!  eigen((2*mband)**response *mband*nkpt*nsppol)= eigenvalues (hartree) for all bands at each k point
!!  fname= character string giving the root to form the name of the
!!   output WFK or WFQ file if response==0, otherwise it is the filename.
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem=maximum number of k-points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum number of plane waves
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outwf
!!
!! CHILDREN
!!      cg_zcopy,cwtime,flush_unit,m_header_init,mask2blocks,timab
!!      unpack_eneocc,wfk_close,wfk_open_write,wfk_write_band_block,wrtout
!!      xmpi_barrier,xmpi_exch
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gs_write_cg(fname,Hdr,mpw,mband,nband,nkpt,nsppol,nspinor,mcg,mkmem,eigen,occ,cg,kg,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_wfk
 use m_header

 use m_io_tools,       only : iomode_from_fname, get_unit, flush_unit
 use m_numeric_tools,  only : mask2blocks
 use m_ebands,         only : unpack_eneocc
 use m_cgtools,        only : cg_zcopy
 use m_time,           only : cwtime

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gs_write_cg'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mkmem,mpw,nkpt,nsppol,nspinor
 character(len=*),intent(in) :: fname
 type(MPI_type),intent(in) :: mpi_enreg
 type(hdr_type),intent(in) :: hdr
!arrays
 integer, intent(in) :: nband(nkpt*nsppol),kg(3,mpw*mkmem)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),occ(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig0=0,tim_rwwf6=6,master=0
 integer :: iomode,icg,ikg,ikpt,spin,my_rank,my_nspinor,nband_k,npw_k,comm,blk,nblocks 
 integer :: ierr,nprocs,action,source
 character(len=500) :: msg
 logical :: ihave_data,iam_master,single_writer
 real(dp) :: cpu,wall,gflops
 type(wfk_t) :: Wfk
!arrays
 integer :: band_block(2)
 integer,pointer :: blocks(:,:)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: eigen3d(:,:,:),occ3d(:,:,:),cg_k(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!FIXME
 call m_header_init(ierr)
 ABI_CHECK(ierr==0,"header_init returned ierr!=0")

 call timab(270+tim_rwwf6,1,tsec)

!Init comm and my_rank
 comm       = mpi_enreg%comm_cell
 nprocs     = xcomm_size(comm)
 my_rank    = mpi_enreg%me_kpt
 iam_master = (master==my_rank)

 my_nspinor=MAX(1,nspinor/mpi_enreg%nproc_spinor)

 ABI_MALLOC(eigen3d, (mband,nkpt,nsppol))
 ABI_MALLOC(occ3d,   (mband,nkpt,nsppol))

 call unpack_eneocc(nkpt,nsppol,mband,nband,eigen,eigen3d)
 call unpack_eneocc(nkpt,nsppol,mband,nband,occ,occ3d)

 iomode = iomode_from_fname(fname)
 iomode = IO_MODE_FORTRAN

 single_writer = (ANY(iomode == (/IO_MODE_FORTRAN, IO_MODE_ETSF/) ))

!if (single_writer .and. paral_kgb) then
!MSG_ERROR("")
!end if

 write(msg,'(3a,i0)')ABI_FUNC//': writing GS WFK file ',TRIM(fname),", with iomode ",iomode
 call wrtout(std_out,msg,'PERS')
 write(std_out,*)"single_writer ",single_writer
 call flush_unit(std_out)

 if (.not.single_writer) then

   call cwtime(cpu,wall,gflops,"start")

   if (iam_master) then
     call wfk_open_write(Wfk,Hdr,fname,formeig0,iomode,get_unit(),xmpi_self,write_hdr=.TRUE.,write_frm=.FALSE.)
   end if
   call xmpi_barrier(comm)
   
   if (.not.iam_master) then
     call wfk_open_write(Wfk,Hdr,fname,formeig0,iomode,get_unit(),xmpi_self,write_hdr=.FALSE.,write_frm=.FALSE.)
   end if

   ABI_CHECK(ALL(Hdr%nband == nband),"nband")

   call cwtime(cpu,wall,gflops,"stop")
   write(msg,'(2(a,f8.2))')" open_write, cpu: ",cpu,", wall: ",wall
   call wrtout(std_out,msg,"PERS")

   call cwtime(cpu,wall,gflops,"start")

   icg=0
   do spin=1,nsppol
     ikg=0
     do ikpt=1,nkpt
       nband_k = Wfk%nband(ikpt,spin)
       npw_k   = Hdr%npwarr(ikpt)

!      write(std_out,*) (mpi_enreg%proc_distrb(ikpt,:,spin)==my_rank)
!      call flush_unit(std_out)

!      Compute my block of bands for this k-point and spin.
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,my_rank)) CYCLE

       call mask2blocks(mpi_enreg%proc_distrb(ikpt,:,spin)==my_rank, nblocks,blocks)
       ABI_CHECK(nblocks==1,"nblocks !=1")

       write(msg,"(a,3(i0,2x))")"Will write (ikpt, spin, nblocks)",ikpt,spin,nblocks
       call wrtout(std_out,msg,"PERS")
       call flush_unit(std_out)

       do blk=1,nblocks
         band_block = blocks(:,blk)

         !write(std_out,*)1+ikg,1+ikg+npw_k,SIZE(kg,DIM=2)
         !write(std_out,*)1+icg,1+icg+npw_k*nband_k, SIZE(cg,DIM=2)

         if (band_block(1)==1) then
           ! This processor writes also kg_k, eig_k and occ_k
           call wfk_write_band_block(Wfk,band_block,ikpt,spin,xmpio_single,kg_k=kg(:,1+ikg:),cg_k=cg(:,1+icg:),&
&           eig_k=eigen3d(:,ikpt,spin),occ_k=occ3d(:,ikpt,spin))
         else 
           call wfk_write_band_block(Wfk,band_block,ikpt,spin,xmpio_single,cg_k=cg(:,1+icg:))
         end if
       end do

       ABI_FREE(blocks)

       !call wrtout(std_out,"Done","PERS")

       icg = icg+npw_k*my_nspinor*nband_k
       ikg = ikg+npw_k
     end do 
   end do 

   ABI_FREE(eigen3d)
   ABI_FREE(occ3d)

   call wfk_close(Wfk)
   call xmpi_barrier(comm)

   call cwtime(cpu,wall,gflops,"stop")
   write(msg,'(2(a,f8.2))')" write all cg cpu: ",cpu,", wall: ",wall
   call wrtout(std_out,msg,"PERS")

 else

   ABI_MALLOC(kg_k,(3,mpw))
   ABI_MALLOC(cg_k,(2,mpw*my_nspinor*mband))
   ABI_CHECK_ALLOC("out of memory in cg_k")

   if (iam_master) then
     call wfk_open_write(Wfk,Hdr,fname,formeig0,iomode,get_unit(),xmpi_self,write_hdr=.TRUE., write_frm=.FALSE.)
   end if

   icg=0
   do spin=1,nsppol
     ikg=0
     do ikpt=1,nkpt
       nband_k = nband(ikpt+(spin-1)*nkpt)
       npw_k   = Hdr%npwarr(ikpt)

       call xmpi_barrier(comm)

!      Transfer the wavefunctions and the g-vectors to the master processor
       source = MINVAL(mpi_enreg%proc_distrb(ikpt,1:nband_k,spin))
       ihave_data = (source==my_rank)

       action=0
       if (iam_master .and. ihave_data)    action=1 ! I am the master node, and I have the data in cg
       if (.not.iam_master.and.ihave_data) action=2 ! I am not the master, and I have the data => send to master
       if (iam_master.and..not.ihave_data) action=3 ! I am the master, and I receive the data

       if (action==1) then ! Copy from kg and cg
         kg_k(:,1:npw_k) = kg(:,ikg+1:ikg+npw_k)
         call cg_zcopy(npw_k*my_nspinor*nband_k, cg(1,icg+1), cg_k)
       end if

       ! Exchange data
       if (action==2.or.action==3) then
         call timab(48,1,tsec)
         if (action==2) then
           call xmpi_exch(kg(:,1+ikg:npw_k+ikg),3*npw_k,source,kg_k,master,comm,ierr)
           call xmpi_exch(cg(:,icg+1:icg+nband_k*npw_k*my_nspinor),2*nband_k*npw_k*my_nspinor,&
&           source,cg_k,master,comm,ierr)
         else
           call xmpi_exch(kg_k,3*npw_k,source,kg_k,master,comm,ierr)
           call xmpi_exch(cg_k,2*nband_k*npw_k*my_nspinor,source,cg_k,master,comm,ierr)
         end if
         call timab(48,2,tsec)
       end if

       ! Master writes this block of bands.
       if (iam_master) then
         band_block = (/1,nband_k/)
         call wfk_write_band_block(Wfk,band_block,ikpt,spin,xmpio_single,kg_k=kg_k,cg_k=cg_k,&
&         eig_k=eigen3d(:,ikpt,spin),occ_k=occ3d(:,ikpt,spin))
       end if

       if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,my_rank))) then
         icg = icg+npw_k*my_nspinor*nband_k
         ikg = ikg+npw_k
       end if

     end do !ikpt
   end do !spin

   ABI_FREE(kg_k)
   ABI_FREE(cg_k)

   if (iam_master) then
     call wfk_close(Wfk)
   end if
   call xmpi_barrier(comm)
 end if

 call timab(270+tim_rwwf6,2,tsec)

 DBG_EXIT("COLL")

end subroutine gs_write_cg
!!***
