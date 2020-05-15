!{\src2tex{textfont=tt}}
!!****f* ABINIT/outwf
!! NAME
!! outwf
!!
!! FUNCTION
!! Conduct output of a "wave-functions" file.
!!  - Compute the maximal residual
!!  - Then open a permanent file wff2 for final output of wf data
!!  - Create a new header for the file.
!!  - Write wave-functions (and energies)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, AR, MB, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mcg)=wavefunction array (storage if nkpt>1)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen( (2*mband)**response *mband*nkpt*nsppol)= eigenvalues (hartree) for all bands at each k point
!!  filnam= character string giving the root to form the name of the
!!   output WFK or WFQ file if response==0, otherwise it is the filename.
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kptns(3,nkpt)=k points in terms of recip primitive translations
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem=Number of k-points treated by this node.
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of plane waves
!!  mxfh=last dimension of the xfhist array
!!  natom=number of atoms in unit cell
!!  nband=number of bands
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nstep=desired number of electron iteration steps
!!  nxfh=actual number of (x,f) history pairs, see xfhist array.
!!  occ(mband*nkpt*nsppol)=occupations for all bands at each k point
!!  resid(mband*nkpt*nsppol)=squared residuals for each band and k point
!!   where resid(n,k)=|<C(n,k)|(H-e(n,k))|C(n,k)>|^2 for the ground state
!!  response: if == 0, GS wavefunctions , if == 1, RF wavefunctions
!!  unwff2=unit for output of wavefunction
!!  xfhist(3,natom+4,2,mxfh)=(x,f) history array, also includes rprim and stress
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! * The name of the file wff2 might be the same as that of the file wff1.
!!
!! PARENTS
!!      berryphase_new,gstate,loper3
!!
!! CHILDREN
!!      cwtime,dfpt_write_cg,gs_write_cg,hdr_io,hdr_io_etsf,rwwf,timab,wffclose
!!      wffkg,wffoffset,wffopen,wfk_diff,wrtout,wvl_write,xmpi_barrier
!!      xmpi_exch
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outwf(cg,dtset,eigen,filnam,hdr,kg,kptns,mband,mcg,mkmem,&
&                mpi_enreg,mpw,mxfh,natom,nband,nkpt,npwarr,&
&                nsppol,nstep,nxfh,occ,resid,response,unwff2,&
&                wfs,wvl,xfhist)

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_wffile
 use m_wfk

 use m_time,         only : cwtime
 use m_header,       only : hdr_skip, hdr_io_etsf, hdr_io

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outwf'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_56_io_mpi
 use interfaces_62_wvl_wfs
 use interfaces_67_common, except_this_one => outwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mkmem,mpw,mxfh,natom,nkpt,nsppol,nstep,nxfh,response,unwff2
 character(len=fnlen),intent(in) :: filnam
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(hdr_type), intent(inout) :: hdr
 type(wvl_wf_type),intent(in) :: wfs
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer, intent(in) :: kg(3,mpw*mkmem),nband(nkpt*nsppol),npwarr(nkpt)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp), intent(in) :: eigen((2*mband)**response*mband*nkpt*nsppol),kptns(3,nkpt)
 real(dp), intent(in) :: occ(mband*nkpt*nsppol),resid(mband*nkpt*nsppol)
 real(dp), intent(in) :: xfhist(3,natom+4,2,mxfh)

!Local variables-------------------------------
 integer,parameter :: nkpt_max=50
 integer :: accesswff,action,band_index,fform,formeig,iband,ibdkpt,icg
 integer :: ierr,ii,ikg,ikpt,spin,ixfh,master,mcg_disk,me,me0,my_nspinor
 integer :: nband_k,nkpt_eff,nmaster,npw_k,option,rdwr,sender,source
 integer :: spaceComm,spaceComm_io,spacecomsender,spaceWorld,sread,sskip,tim_rwwf,xfdim2
#ifdef HAVE_MPI
 integer :: ipwnbd
#endif
 real(dp) :: residk,residm,resims,cpu,wall,gflops
 logical :: ihave_data,iwrite,iam_master
 character(len=500) :: msg
 type(wffile_type) :: wff2
 integer,allocatable :: kg_disk(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cg_disk(:,:),eig_k(:),occ_k(:)

! *************************************************************************
!For readability of the source file, define a "me" variable also in the sequential case

 DBG_ENTER("COLL")

 xfdim2 = natom+4
!Init mpi_comm
 spaceWorld= mpi_enreg%comm_cell
 spaceComm=spaceWorld
 spaceComm_io=xmpi_self

 if (mpi_enreg%paral_kgb==1 ) spaceComm_io= mpi_enreg%comm_bandspinorfft
 if (mpi_enreg%paral_kgb==1 ) spaceComm= mpi_enreg%comm_cell

!Paral_kgb=1 and Fortran-I/O is not supported (only for testing purpose)
 if (mpi_enreg%paral_kgb==1.and.dtset%accesswff==IO_MODE_FORTRAN) then
   spaceWorld=mpi_enreg%comm_kpt
   write(msg,'(7a)') &
&   'WF file is written using standard Fortran I/O',ch10,&
&   'and Kpt-band-FFT parallelization is active !',ch10,&
&   'This is only allowed for testing purposes.',ch10,&
&   'The produced WF file will be incomplete and not useable.'
   MSG_WARNING(msg)
 end if

!If parallel HF calculation
 if (mpi_enreg%paral_hf==1 ) spaceComm_io= mpi_enreg%comm_hf
 if (mpi_enreg%paral_hf==1 ) spaceComm= mpi_enreg%comm_cell

!Paral_hf=1 and Fortran-I/O is not supported (copy from paral_kgb... not tested)
 if (mpi_enreg%paral_hf==1.and.dtset%accesswff==IO_MODE_FORTRAN) then
   spaceWorld=mpi_enreg%comm_kpt
   write(msg,'(7a)') &
&   'WF file is written using standard Fortran I/O',ch10,&
&   'and HF parallelization is active !',ch10,&
&   'This is only allowed for testing purposes.',ch10,&
&   'The produced WF file will be incomplete and not useable.'
   MSG_WARNING(msg)
 end if


!Init me
 me=mpi_enreg%me_kpt
 me0=me
!Define master
 master=0

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 tim_rwwf =0
 source = master
 sread = master
 iam_master=(master==me)
 iwrite=iam_master
 sender=-1

!Compute mean square and maximum residual over all bands and k points and spins
!(disregard k point weights and occupation numbers here)
 band_index=sum(nband(1:nkpt*nsppol))
 resims=sum(resid(1:band_index))/dble(band_index)

!Find largest residual over bands, k points, and spins, except for nbdbuf highest bands
!Already AVAILABLE in hdr ?!
 ibdkpt=1
 residm=zero
 do spin=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(spin-1)*nkpt)
     nband_k=max(1,nband_k-dtset%nbdbuf)
     residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_k-1)))
     ibdkpt=ibdkpt+nband_k
   end do
 end do

 write(msg,'(a,1p,e12.4,a,e12.4)')' Mean square residual over all n,k,spin= ',resims,'; max=',residm
 call wrtout(ab_out,msg,'COLL')

 band_index=0
 nkpt_eff=nkpt
 if( (dtset%prtvol==0 .or. dtset%prtvol==1) .and. nkpt_eff>nkpt_max ) nkpt_eff=nkpt_max

!Loop over spin again
 do spin=1,nsppol
!  Give (squared) residuals for all bands at each k
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(spin-1)*nkpt)
!    Will not print all residuals when prtvol=0 or 1
     if(ikpt<=nkpt_eff)then
!      Find largest residual over all bands for given k point
       residk=maxval(resid(1+band_index:nband_k+band_index))
       write(msg,'(1x,3f8.4,3x,i2,1p,e13.5,a)')kptns(1:3,ikpt),spin,residk,' kpt; spin; max resid(k); each band:'
       call wrtout(ab_out,msg,'COLL')
       do ii=0,(nband_k-1)/8
         write(msg,'(1x,1p,8e9.2)')(resid(iband+band_index),iband=1+ii*8,min(nband_k,8+ii*8))
         call wrtout(ab_out,msg,'COLL')
       end do
     else if(ikpt==nkpt_eff+1)then
       write(msg,'(2a)')' outwf : prtvol=0 or 1, do not print more k-points.',ch10
       call wrtout(ab_out,msg,'COLL')
     end if
     band_index=band_index+nband_k
   end do
 end do

!Will write the wavefunction file only when nstep>0
 if (nstep>0 .and. dtset%prtwf/=0) then

   call cwtime(cpu,wall,gflops,"start")

!  Only the master write the file, except if MPI I/O, but the
!  full wff dataset should be provided to WffOpen in this case
   accesswff=IO_MODE_FORTRAN_MASTER
   if (dtset%accesswff==IO_MODE_MPI)  accesswff = IO_MODE_MPI
   if (dtset%accesswff==IO_MODE_ETSF) accesswff = IO_MODE_ETSF
!  accesswff=IO_MODE_MPI

   write(msg,'(4a,i0)')ch10,' outwf: write wavefunction to file ',trim(filnam),", with accesswff ",accesswff
   call wrtout(std_out,msg,'COLL')

   call WffOpen(accesswff,spaceComm,filnam,ierr,wff2,master,me0,unwff2,spaceComm_io)
!  Conduct wavefunction output to wff2

   ABI_ALLOCATE(kg_disk,(3,mpw))

   mcg_disk=mpw*my_nspinor*mband
   formeig=0; if (response==1) formeig=1

   ABI_ALLOCATE(eig_k,( (2*mband)**formeig * mband))
   ABI_ALLOCATE(occ_k,(mband))

#ifdef HAVE_MPI
   call xmpi_barrier(spaceComm)
!  Compute mband and mpw
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))
   ABI_CHECK_ALLOC("out of memory in cg_disk")
#endif

   band_index=0
   icg=0
   if(mpi_enreg%paralbd==0) tim_rwwf=6
   if(mpi_enreg%paralbd==1) tim_rwwf=12

!  Write header info for new wf file
   rdwr=2
   if (dtset%usewvl==0) then
     fform=2
   else
     fform = 200 ! Use 200 as radical for naming file format used by wavelets.
   end if

   if (wff2%accesswff < 2) then
     call hdr_io(fform,hdr,rdwr,wff2)
     call WffKg(wff2,1)
   else if (wff2%accesswff==IO_MODE_ETSF .and. iam_master) then
     call hdr_io_etsf(fform, hdr, rdwr, wff2%unwff)
   end if

   do spin=1,nsppol
     ikg=0

     do ikpt=1,nkpt
       nband_k=nband(ikpt+(spin-1)*nkpt)
       npw_k=npwarr(ikpt)

#ifdef HAVE_MPI
       if (dtset%usewvl == 0) then
         call xmpi_barrier(spaceWorld)

!        Must transfer the wavefunctions to the master processor
!        Separate sections for paralbd=1 or other values ; might be merged
         if(mpi_enreg%paralbd==0)then
           nmaster=0
           source=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,spin))
           ihave_data=.false.
           if(source==me)ihave_data=.true.
           action=0
!          I am the master node, and I have the data in cg or cg_disk
           if((iam_master).and.(ihave_data))action=1
!          I am not the master, and I have the data => send to master
           if((.not.iam_master).and.(ihave_data))action=2
!          I am the master, and I receive the data
           if((iam_master).and.(.not.ihave_data))action=3

!          I have the data in cg or cg_disk ( MPI_IO case)
           if (accesswff==IO_MODE_MPI) then
             action = 0
             sender=-1
             iwrite=.false.
             if (ihave_data)then
               action=1
               iwrite=.true.
               sender=me
             end if
           end if

!          I am the master node, and I have the data in cg or cg_disk
!          I have the data in cg or cg_disk ( MPI_IO case)
           if(action==1)then
!            Copy from kg to kg_disk
             kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
!            Copy from cg to cg_disk
             do ipwnbd=1,nband_k*npw_k*my_nspinor
               cg_disk(1,ipwnbd)=cg(1,ipwnbd+icg)
               cg_disk(2,ipwnbd)=cg(2,ipwnbd+icg)
             end do
           end if

!          I am not the master, and I have the data => send to master
!          I am the master, and I receive the data
           if ( action==2.or.action==3) then
             !write(std_out,*)npw_k,nband_k
             call timab(48,1,tsec)
             if(action==2)then
               call xmpi_exch(kg(:,1+ikg:npw_k+ikg),3*npw_k,source,kg_disk,nmaster,spaceWorld,ierr)
               call xmpi_exch(cg(:,icg+1:icg+nband_k*npw_k*my_nspinor),2*nband_k*npw_k*my_nspinor, &
&               source,cg_disk,nmaster,spaceWorld,ierr)
             else
               call xmpi_exch(kg_disk,3*npw_k,source,kg_disk,nmaster,spaceWorld,ierr)
               call xmpi_exch(cg_disk,2*nband_k*npw_k*my_nspinor,source,cg_disk,nmaster,spaceWorld,ierr)
             end if
             call timab(48,2,tsec)
           end if


         else if(mpi_enreg%paralbd==1)then
           nmaster=0
#ifdef HAVE_MPI_IO
           sender=IO_MODE_FORTRAN_MASTER
           if( accesswff==IO_MODE_MPI) then
             nmaster=mpi_enreg%proc_distrb(ikpt,1,spin)
             sender=nmaster
           end if
#endif

!          Note the loop over bands
           do iband=1,nband_k

!            The message passing related to kg is counted as one band
             action=0

!            I am the master node, and I have the data in cg or cg_disk
             if( mpi_enreg%proc_distrb(ikpt,iband,spin)==nmaster .and. me==nmaster) then
               action=1
!              I am not the master, and I have the data => send to master
             elseif( mpi_enreg%proc_distrb(ikpt,iband,spin)==me .and. me/=nmaster ) then
               action = 2
!              I am the master, and I receive the data
             elseif( mpi_enreg%proc_distrb(ikpt,iband,spin)/=me .and. me==nmaster ) then
               action=3
             end if

             if(action==1) then
!              I am the master node, and I have the data in cg or cg_disk
!              Copy from kg to kg_disk
               if(iband==1)kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
!              Copy from cg to cg_disk
               do ipwnbd=1,npw_k*my_nspinor
                 cg_disk(1,(iband-1)*npw_k*my_nspinor+ipwnbd) = cg(1,(iband-1)*npw_k*my_nspinor+ipwnbd+icg)
                 cg_disk(2,(iband-1)*npw_k*my_nspinor+ipwnbd) = cg(2,(iband-1)*npw_k*my_nspinor+ipwnbd+icg)
               end do
             end if  ! action=1

             if ( action==2.or.action==3) then
!              action=2 :  I am not the master, and I have the data => send to master
!              action=3 :  I am the master, and I receive the data
               call timab(48,1,tsec)
               if ( iband == 1 ) then
                 if (action==2) then
                   call xmpi_exch(kg(:,1+ikg:npw_k+ikg),3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,spin), &
&                   kg_disk,nmaster,spaceWorld,ierr)
                 else
                   call xmpi_exch(kg_disk,3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,spin),  &
&                   kg_disk,nmaster,spaceWorld,ierr)
                 end if
               end if       ! iband =1
               ipwnbd=(iband-1)*npw_k*my_nspinor
               if (action==2) then
                 call xmpi_exch( cg(:,ipwnbd+icg+1:ipwnbd+icg+npw_k*my_nspinor),2*npw_k*my_nspinor &
&                 ,mpi_enreg%proc_distrb(ikpt,iband,spin)                    &
&                 ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*my_nspinor),nmaster,spaceWorld,ierr)
               else
                 call xmpi_exch( cg_disk(:,ipwnbd+1:ipwnbd+npw_k*my_nspinor),2*npw_k*my_nspinor    &
&                 ,mpi_enreg%proc_distrb(ikpt,iband,spin)                    &
&                 ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*my_nspinor),nmaster,spaceWorld,ierr)
               end if

               call timab(48,2,tsec)
             end if        ! action=2 or action=3

             if(accesswff==IO_MODE_MPI) then
!              I have the data in cg or cg_disk
               iwrite=.false.
               if (nmaster == me) iwrite=.true.
             end if

           end do ! End of loop over bands
         end if ! End of paralbd=1
       end if
#endif

!      Only the master will write to disk the final output wf file.
!      in MPI_IO case only iwrite will write to disk the final output wf file.
       if(iwrite) then
!        DEBUG
!        write(std_out,*) 'outwf : I am master and will write wf file'
!        ENDDEBUG
         if(formeig==0)then
           eig_k(1:nband_k)=eigen(1+band_index:nband_k+band_index)
           occ_k(1:nband_k)=occ(1+band_index:nband_k+band_index)
         else
           eig_k(1:2*nband_k*nband_k)=eigen(1+band_index:2*nband_k*nband_k+band_index)
         end if
         option=2
         if(dtset%prtwf==3)option=5
!        if (dtset%prtwf == 2 .and. mkmem/=0) option=4

         if (dtset%usewvl == 0) then
#ifdef HAVE_MPI
           call rwwf(cg_disk,eig_k,formeig,0,0,ikpt,spin,kg_disk,mband,mcg_disk,mpi_enreg, &
&           nband_k, nband_k,npw_k,my_nspinor,occ_k,option,1,tim_rwwf,wff2)

#else
           kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
           call rwwf(cg,eig_k,formeig,0,icg,ikpt,spin,kg_disk,mband,mcg,mpi_enreg,nband_k, &
&           nband_k, npw_k,my_nspinor,occ_k,option,1,tim_rwwf,wff2)
#endif
         else
           call wvl_write(dtset,eigen,mpi_enreg,option,hdr%rprimd,wff2,wfs,wvl,hdr%xred)
         end if
       end if

!      The wavefunctions for the present k point and spin are written
       if(response==0)band_index=band_index+nband_k
       if(response==1)band_index=band_index+2*nband_k*nband_k

       sskip=1
#ifdef HAVE_MPI
       if (dtset%usewvl == 0) then
         sskip=0
         if(.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,me)))sskip=1
       end if
#endif
       if(sskip==1)then
         icg=icg+npw_k*my_nspinor*nband_k
         ikg=ikg+npw_k
       end if


#ifdef HAVE_MPI_IO
       spacecomsender=spaceComm
       if (mpi_enreg%paral_kgb==1) spacecomsender =mpi_enreg%comm_kpt
       if (mpi_enreg%paral_hf==1) spacecomsender =mpi_enreg%comm_kpt
       call WffOffset(wff2,sender,spacecomsender,ierr)
#endif

     end do ! ikpt
   end do ! spin

   ABI_DEALLOCATE(kg_disk)
#ifdef HAVE_MPI
   ABI_DEALLOCATE(cg_disk)
#endif

   ABI_DEALLOCATE(eig_k)
   ABI_DEALLOCATE(occ_k)

!  Write the (x,f) history
   if(me0==0 .and. nxfh>0 .and. response==0)then
     if (wff2%accesswff /= 2) then
#ifdef HAVE_MPI_IO
       if(wff2%accesswff==IO_MODE_MPI) then
         close(unit=wff2%unwff)
!        the file is to be positioned at the terminal point
         open(unit=wff2%unwff,file=wff2%fname,form='unformatted',POSITION="APPEND")
       end if
#endif
       write(unit=wff2%unwff)nxfh
       do ixfh=1,nxfh
         write(unit=wff2%unwff)xfhist(:,:,:,ixfh)
       end do
     end if
   end if

!  Close the wavefunction file (and do NOT delete it !)
   if (wff2%accesswff /= IO_MODE_NETCDF) then
     call WffClose(wff2,ierr)
   end if

   call cwtime(cpu,wall,gflops,"stop")
   write(msg,'(a,i0,2(a,f8.2))')" outwf with iomode: ",accesswff,", cpu: ",cpu,", wall: ",wall
   call wrtout(std_out,msg,"PERS")
 end if ! End condition of nstep>0

 if (.FALSE.) then
!  if (dtset%usewvl==0 .and. mpi_enreg%paral_kgb==0 .and. nstep>0 .and. dtset%prtwf/=0) then
   if (response==0) then
     formeig=0
     call gs_write_cg(TRIM(filnam)//".TEST",Hdr,mpw,mband,nband,nkpt,nsppol,dtset%nspinor,mcg,mkmem,eigen,occ,cg,kg,mpi_enreg)
   else if (response==1) then
     formeig=1
     call dfpt_write_cg(TRIM(filnam)//".TEST",Hdr,mpw,mband,nband,nkpt,nsppol,dtset%nspinor,mcg,mkmem,eigen,cg,kg,mpi_enreg)
   else 
     MSG_ERROR("Wrong value for response")
   end if
   if (me==master) then
     call wfk_diff(filnam,TRIM(filnam)//".TEST",formeig,xmpi_self,ierr)
     ABI_CHECK(ierr==0,"wfk_diff returned ierr != 0")
   end if
 end if

 DBG_EXIT("COLL")

end subroutine outwf
!!***
