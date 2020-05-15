!{\src2tex{textfont=tt}}
!!****f* ABINIT/getghc
!!
!! NAME
!! getghc
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space
!! Result is put in array ghc.
!! <G|Vnonlocal|C> is also returned in gvnlc.
!! if required, <G|S|C> is returned in gsc (S=overlap - PAW only)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, LSI, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cpopt=flag defining the status of cwaveprj%cp(:)=<Proj_i|Cnk> scalars (PAW only)
!!       (same meaning as in nonlop.F90 routine)
!!       if cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
!!       if cpopt= 0, <p_lmn|in> are computed here and saved
!!       if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!       if cpopt= 2  <p_lmn|in> are already in memory;
!!       if cpopt= 3  <p_lmn|in> are already in memory; first derivatives are computed here and saved
!!       if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!! cwavef(2,npw*nspinor*ndat)=planewave coefficients of wavefunction.
!! dimffnl=second dimension of ffnl (1+number of derivatives)
!! ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!! filstat=name of the status file
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!! kg_k(3,npw)=G vec coordinates wrt recip lattice transl.
!! kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!! lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!        Typically lambda is the eigenvalue (or its guess)
!! mpi_enreg=informations about MPI parallelization
!! natom=number of atoms in unit cell.
!! ndat=number of FFT to do in //
!! npw=number of planewaves in basis for given k point.
!! nspinor=number of spinorial components of the wavefunctions (on current proc)
!! ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!! prtvol=control print volume and debugging output
!! sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!    (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                      if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!! tim_getghc=timing code of the calling subroutine(can be set to 0 if not attributed)
!! type_calc= option governing which part of Hamitonian is to be applied:
!             0: whole Hamiltonian
!!            1: local part only
!!            2: non-local+kinetic only (added to the exixting Hamiltonian)
!!            3: local + kinetic only (added to the existing Hamiltonian)
!! vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!! fock <type(fock_type)>= quantities to calculate Fock exact exchange
!! vxctaulocal(n4,n5,n6,nvloc,4)= local potential corresponding to the derivative of XC energy with respect to
!!  kinetic energy density, in real space, on the augmented fft grid. (optional argument)
!!  This array contains also the gradient of vxctaulocal (gvxctaulocal) in vxctaulocal(:,:,:,:,2:4).
!!
!! OUTPUT
!!  ghc(2,npw*nspinor*ndat)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                       or <G|H-lambda.S|C> (if sij_opt=-1)
!!  gvnlc(2,npw*nspinor*ndat)=matrix elements <G|Vnonlocal|C> (if sij_opt>=0)
!!                                         or <G|Vnonlocal-lambda.S|C> (if sij_opt=-1)
!!    (sometimes desired for computing nonlocal part of total energy, but can be ignored).
!!  if (sij_opt=1)
!!    gsc(2,npw*nspinor*ndat)=matrix elements <G|S|C> (S=overlap).
!!
!! SIDE EFFECTS
!!  ====== if gs_ham%usepaw==1 and
!!  cwaveprj(natom,nspinor*(1+cpopt)*ndat)= wave functions at k projected with nl projectors
!!
!! PARENTS
!!      cgwf,cgwf3,chebfi,ks_ddiago,lobpcgIIwf,lobpcgccIIIwf,lobpcgccIIwf
!!      lobpcgwf,m_lobpcgIIIwf,mkresi,prep_getghc,update_orbmag
!!
!! CHILDREN
!!      fock_getghc,fourwf,nonlop,status,timab,wrtout,xmpi_alltoallv,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getghc(cpopt,cwavef,cwaveprj,dimffnl,ffnl,filstat,ghc,gsc,gs_ham,&
&  gvnlc,kg_k,kinpw,lambda,mpi_enreg,natom,ndat,npw,nspinor,&
&  paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,type_calc,vlocal,fock,&
&  ikpt_this_proc,vxctaulocal) ! optional argument

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_pawcprj,     only : pawcprj_type
 use m_bandfft_kpt, only : bandfft_kpt
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_fock,        only : fock_type, fock_get_getghc_call

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getghc'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_53_ffts
 use interfaces_65_nonlocal
 use interfaces_66_fock
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpopt,dimffnl,natom,ndat,npw,nspinor,paral_kgb, prtvol
 integer,intent(in) :: sij_opt,tim_getghc,type_calc
 integer,intent(in),optional :: ikpt_this_proc
 real(dp),intent(in) :: lambda !! TODO this should be (ndat)
 character(len=*),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_ham
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(in) :: ffnl(npw,dimffnl,gs_ham%lmnmax,gs_ham%ntypat),kinpw(npw)
 real(dp),intent(inout) :: cwavef(2,npw*nspinor*ndat)
 real(dp),intent(inout) :: ghc   (2,npw*nspinor*ndat)
 real(dp),intent(inout) :: gvnlc (2,npw*nspinor*ndat)
 real(dp),intent(inout) :: ph3d(2,npw,gs_ham%matblk),vlocal(gs_ham%n4,gs_ham%n5,gs_ham%n6,gs_ham%nvloc)
 real(dp),intent(out) :: gsc(2,npw*nspinor*ndat*((sij_opt+1)/2))
 real(dp),intent(inout), optional :: vxctaulocal(gs_ham%n4,gs_ham%n5,gs_ham%n6,gs_ham%nvloc,4)
 type(fock_type),pointer, intent(inout) :: fock
 type(pawcprj_type),intent(inout) :: cwaveprj(natom,nspinor*((cpopt+5)/5)*gs_ham%usepaw*ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=114,re=1,im=2,tim_fourwf=1
 integer :: choice,cplex,cpopt_here,i1,i2,i3,idat,idir,ierr,iexit
 integer :: ig,igspinor,ii,iispinor,ipw,ispinor,nkpg,nnlout,npw_fft,nspinortot
 integer :: paw_opt,shift,signs,tim_nonlop
 logical :: have_to_reequilibrate,nspinor1TreatedByThisProc,nspinor2TreatedByThisProc 
 logical :: compute_fock
 real(dp) :: ghcim,ghcre,gp2pi1,gp2pi2,gp2pi3,kpt_cart,kg_k_cart,weight
 character(len=500) :: message
!arrays
 integer,ABI_CONTIGUOUS pointer :: indices_pw_fft(:),kg_k_fft(:,:)
 integer,ABI_CONTIGUOUS pointer :: recvcount_fft(:),recvdisp_fft(:),sendcount_fft(:),senddisp_fft(:)
 real(dp) :: enlout(ndat),nonlop_dum(1,1),nonlop_dum2(1,1),tsec(2)
 real(dp),allocatable :: buff_wf(:,:),cwavef1(:,:),cwavef2(:,:)
 real(dp),allocatable :: cwavef_fft(:,:),cwavef_fft_tr(:,:)
 real(dp),allocatable :: gcwavef(:,:,:),gcwavef1(:,:,:),gcwavef2(:,:,:)
 real(dp),allocatable :: ghc1(:,:),ghc2(:,:),ghc3(:,:),ghc4(:,:)
 real(dp) :: kpg_dum(0,0)
 real(dp),allocatable :: lcwavef(:,:),lcwavef1(:,:),lcwavef2(:,:)
 real(dp),allocatable :: vlocal_tmp(:,:,:),work(:,:,:,:)
 real(dp) :: lambda_ndat(ndat)

! *********************************************************************

!Keep track of total time spent in getghc:
 call timab(200+tim_getghc,1,tsec)

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' getghc : enter, debugging '
   call wrtout(std_out,message,'PERS')
 end if

!Parallelization over spinors management
 nspinortot=min(2,(1+mpi_enreg%paral_spinor)*nspinor)
 if (mpi_enreg%paral_spinor==0) then
   shift=npw
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(nspinortot==2)
 else
   shift=0
   nspinor1TreatedByThisProc=(mpi_enreg%me_spinor==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spinor==1)
 end if

!Determine if a FFT load balancing is OK or not
 have_to_reequilibrate = .false.
 if(mpi_enreg%paral_kgb == 1) then
   if(.not.(present(ikpt_this_proc))) then
     MSG_BUG("Reequilibration step is impossible without passing ikpt_this_proc as argument")
   end if
   have_to_reequilibrate = bandfft_kpt(ikpt_this_proc)%have_to_reequilibrate
 end if

 if ((type_calc==0).or.(type_calc==1).or.(type_calc==3)) then

!  Eventually adjust load balancing for FFT (by changing FFT distrib)
   if (have_to_reequilibrate) then
     npw_fft =  bandfft_kpt(ikpt_this_proc)%npw_fft
     sendcount_fft  => bandfft_kpt(ikpt_this_proc)%sendcount_fft(:)
     recvcount_fft  => bandfft_kpt(ikpt_this_proc)%recvcount_fft(:)
     senddisp_fft   => bandfft_kpt(ikpt_this_proc)%senddisp_fft(:)
     recvdisp_fft   => bandfft_kpt(ikpt_this_proc)%recvdisp_fft(:)
     indices_pw_fft => bandfft_kpt(ikpt_this_proc)%indices_pw_fft(:)
     kg_k_fft       => bandfft_kpt(ikpt_this_proc)%kg_k_fft(:,:)
     ABI_ALLOCATE(buff_wf,(2,npw*ndat) ) ! for cgwavef sorting
     ABI_ALLOCATE(cwavef_fft,(2,npw_fft*ndat) )
     if(ndat>1) then
       ABI_ALLOCATE(cwavef_fft_tr, (2,npw_fft*ndat))
     end if
!    filling of sorted send buffers before exchange
     do idat=1, ndat
       do ipw = 1 ,npw
         buff_wf(1:2, idat + ndat*(indices_pw_fft(ipw)-1) ) = cwavef(1:2,ipw + npw*(idat-1))
       end do
     end do
     if(ndat > 1) then
       call xmpi_alltoallv(buff_wf,2*ndat*sendcount_fft,2*ndat*senddisp_fft,  &
&       cwavef_fft_tr,2*ndat*recvcount_fft, 2*ndat*recvdisp_fft, mpi_enreg%comm_fft,ierr)
!      We need to transpose data
       do idat=1,ndat
         do ipw = 1 ,npw_fft
           cwavef_fft(1:2,  ipw + npw_fft*(idat-1)) = cwavef_fft_tr(1:2,  idat + ndat*(ipw-1))
         end do
       end do
     else
       call xmpi_alltoallv(buff_wf,2*sendcount_fft,2*senddisp_fft,  &
&       cwavef_fft,2*recvcount_fft, 2*recvdisp_fft, mpi_enreg%comm_fft,ierr)
     end if
   end if

!  Apply the local potential to the wavefunction
!  Start from wavefunction in reciprocal space cwavef
!  End with function ghc in reciprocal space also.
   ABI_ALLOCATE(work,(2,gs_ham%n4,gs_ham%n5,gs_ham%n6*ndat))
   weight=one

!  Application of the local potential
   if (nspinortot==2) then
     ABI_ALLOCATE(cwavef1,(2,npw*ndat))
     ABI_ALLOCATE(cwavef2,(2,npw*ndat))
     do idat=1,ndat
       do ipw=1,npw
         cwavef1(1:2,ipw+(idat-1)*npw)=cwavef(1:2,ipw+(idat-1)*nspinor*npw)
         cwavef2(1:2,ipw+(idat-1)*npw)=cwavef(1:2,ipw+(idat-1)*nspinor*npw+shift)
       end do
     end do
!    call cg_zcopy(npw*ndat,cwavef(1,1),cwavef1)
!    call cg_zcopy(npw*ndat,cwavef(1,1+shift),cwavef2)
   end if

!  Treat scalar local potentials
   if (gs_ham%nvloc==1) then

     if (nspinortot==1) then

       if (have_to_reequilibrate) then
         call fourwf(1,vlocal,cwavef_fft,cwavef_fft,work,gs_ham%gbound,gs_ham%gbound,&
&         gs_ham%istwf_k,kg_k_fft,kg_k_fft,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw_fft,npw_fft,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
       else
         call fourwf(1,vlocal,cwavef,ghc,work,gs_ham%gbound,gs_ham%gbound,&
&         gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
       end if

       if (present(vxctaulocal)) then  !metaGGA
!        Do it in 3 STEPs:
!        STEP1: Compute grad of cwavef and Laplacian of cwavef
         ABI_ALLOCATE(gcwavef,(2,npw*ndat,3))
         ABI_ALLOCATE(lcwavef,(2,npw*ndat))
!$OMP PARALLEL DO
         do idat=1,ndat
           do ipw=1,npw
             gcwavef(:,ipw+(idat-1)*npw,1:3)=zero
             lcwavef(:,ipw+(idat-1)*npw)  =zero
           end do
         end do
         do idir=1,3
           gp2pi1=gs_ham%gprimd(idir,1)*two_pi
           gp2pi2=gs_ham%gprimd(idir,2)*two_pi
           gp2pi3=gs_ham%gprimd(idir,3)*two_pi
           kpt_cart=gp2pi1*gs_ham%kpoint(1)+gp2pi2*gs_ham%kpoint(2)+gp2pi3*gs_ham%kpoint(3)
!          Multiplication by 2pi i (G+k)_idir for gradient
!          Multiplication by -(2pi (G+k)_idir )**2 for Laplacian
           do idat=1,ndat
             do ipw=1,npw
               kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
               gcwavef(1,ipw+(idat-1)*npw,idir)= cwavef(2,ipw+(idat-1)*npw)*kg_k_cart
               gcwavef(2,ipw+(idat-1)*npw,idir)=-cwavef(1,ipw+(idat-1)*npw)*kg_k_cart
               lcwavef(1,ipw+(idat-1)*npw)=lcwavef(1,ipw+(idat-1)*npw)-cwavef(1,ipw+(idat-1)*npw)*kg_k_cart**2
               lcwavef(2,ipw+(idat-1)*npw)=lcwavef(2,ipw+(idat-1)*npw)-cwavef(2,ipw+(idat-1)*npw)*kg_k_cart**2
             end do
           end do
         end do ! idir
!        STEP2: Compute (vxctaulocal)*(Laplacian of cwavef) and add it to ghc
         ABI_ALLOCATE(ghc1,(2,npw*ndat))
         call fourwf(1,vxctaulocal(:,:,:,:,1),lcwavef,ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&         gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
!$OMP PARALLEL DO
         do idat=1,ndat
           do ipw=1,npw
             ghc(:,ipw+(idat-1)*npw)=ghc(:,ipw+(idat-1)*npw)-half*ghc1(:,ipw+(idat-1)*npw)
           end do
         end do
         ABI_DEALLOCATE(ghc1)
         ABI_DEALLOCATE(lcwavef)
!        STEP3: Compute sum of (grad components of vxctaulocal)*(grad components of cwavef)
         ABI_ALLOCATE(ghc1,(2,npw*ndat))
         do idir=1,3
           call fourwf(1,vxctaulocal(:,:,:,:,1+idir),gcwavef(:,:,idir),ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&           gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&           npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&           use_gpu_cuda=gs_ham%use_gpu_cuda)
!$OMP PARALLEL DO
           do idat=1,ndat
             do ipw=1,npw
               ghc(:,ipw+(idat-1)*npw)=ghc(:,ipw+(idat-1)*npw)-half*ghc1(:,ipw+(idat-1)*npw)
             end do
           end do
         end do ! idir
         ABI_DEALLOCATE(ghc1)
         ABI_DEALLOCATE(gcwavef)
       end if ! if present(vxctaulocal) i.e. metaGGA

     else ! nspinortot==2

       if (nspinor1TreatedByThisProc) then

         ABI_ALLOCATE(ghc1,(2,npw*ndat))
         ghc1(:,:)=zero
         call fourwf(1,vlocal,cwavef1,ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&         gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)
         do idat=1,ndat
           do ipw =1, npw
             ghc(1:2,ipw+(idat-1)*nspinor*npw)    =ghc1(1:2,ipw+(idat-1)*npw)
           end do
         end do
         ABI_DEALLOCATE(ghc1)

         if(present(vxctaulocal))then  !metaGGA
!          Do it in 3 STEPs:
!          STEP1: Compute grad of cwavef and Laplacian of cwavef
           ABI_ALLOCATE(gcwavef1,(2,npw*ndat,3))
           ABI_ALLOCATE(lcwavef1,(2,npw*ndat))
!$OMP PARALLEL DO
           do idat=1,ndat
             do ipw=1,npw
               gcwavef1(:,ipw+(idat-1)*npw,1:3)=zero
               lcwavef1(:,ipw+(idat-1)*npw)  =zero
             end do
           end do
           do idir=1,3
             gp2pi1=gs_ham%gprimd(idir,1)*two_pi
             gp2pi2=gs_ham%gprimd(idir,2)*two_pi
             gp2pi3=gs_ham%gprimd(idir,3)*two_pi
             kpt_cart=gp2pi1*gs_ham%kpoint(1)+gp2pi2*gs_ham%kpoint(2)+gp2pi3*gs_ham%kpoint(3)
!            Multiplication by 2pi i (G+k)_idir for gradient
!            Multiplication by -(2pi (G+k)_idir )**2 for Laplacian
             do idat=1,ndat
               do ipw=1,npw
                 kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
                 gcwavef1(1,ipw+(idat-1)*npw,idir)= cwavef1(2,ipw+(idat-1)*npw)*kg_k_cart
                 gcwavef1(2,ipw+(idat-1)*npw,idir)=-cwavef1(1,ipw+(idat-1)*npw)*kg_k_cart
                 lcwavef1(1,ipw+(idat-1)*npw)=lcwavef1(1,ipw+(idat-1)*npw)-cwavef1(1,ipw+(idat-1)*npw)*kg_k_cart**2
                 lcwavef1(2,ipw+(idat-1)*npw)=lcwavef1(2,ipw+(idat-1)*npw)-cwavef1(2,ipw+(idat-1)*npw)*kg_k_cart**2
               end do
             end do
           end do ! idir
!          STEP2: Compute (vxctaulocal)*(Laplacian of cwavef) and add it to ghc
           ABI_ALLOCATE(ghc1,(2,npw*ndat))
           call fourwf(1,vxctaulocal(:,:,:,:,1),lcwavef1,ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&           gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&           npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&           use_gpu_cuda=gs_ham%use_gpu_cuda)
!$OMP PARALLEL DO
           do idat=1,ndat
             do ipw=1,npw
               ghc(:,ipw+(idat-1)*nspinor*npw) = ghc(:,ipw+(idat-1)*nspinor*npw) -half*ghc1(:,ipw+(idat-1)*npw)
             end do
           end do
           ABI_DEALLOCATE(ghc1)
           ABI_DEALLOCATE(lcwavef1)
!          STEP3: Compute (grad components of vxctaulocal)*(grad components of cwavef)
           ABI_ALLOCATE(ghc1,(2,npw*ndat))
           do idir=1,3
             call fourwf(1,vxctaulocal(:,:,:,:,1+idir),gcwavef1(:,:,idir),ghc1,work,gs_ham%gbound,&
&             gs_ham%gbound,gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&             npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&             use_gpu_cuda=gs_ham%use_gpu_cuda)
!$OMP PARALLEL DO
             do idat=1,ndat
               do ipw=1,npw
                 ghc(:,ipw+(idat-1)*nspinor*npw) = ghc(:,ipw+(idat-1)*nspinor*npw) - half*ghc1(:,ipw+(idat-1)*npw)
               end do
             end do
           end do ! idir
           ABI_DEALLOCATE(ghc1)
           ABI_DEALLOCATE(gcwavef1)
         end if ! if present(vxctaulocal) i.e. metaGGA

       end if ! spin 1 treated by this proc

       if (nspinor2TreatedByThisProc) then

         ABI_ALLOCATE(ghc2,(2,npw*ndat))
         ghc2(:,:)=zero

         call fourwf(1,vlocal,cwavef2,ghc2,work,gs_ham%gbound,gs_ham%gbound,&
&         gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&         npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&         use_gpu_cuda=gs_ham%use_gpu_cuda)

         do idat=1,ndat
           do ipw=1,npw
             ghc(1:2,ipw+(idat-1)*nspinor*npw+shift)=ghc2(1:2,ipw+(idat-1)*npw)
           end do
         end do
         ABI_DEALLOCATE(ghc2)

         if(present(vxctaulocal))then  !metaGGA
!          Do it in 3 STEPs:
!          STEP1: Compute grad of cwavef and Laplacian of cwavef
           ABI_ALLOCATE(gcwavef1,(2,npw*ndat,3))
           ABI_ALLOCATE(gcwavef2,(2,npw*ndat,3))
           ABI_ALLOCATE(lcwavef1,(2,npw*ndat))
           ABI_ALLOCATE(lcwavef2,(2,npw*ndat))
!$OMP PARALLEL DO
           do idat=1,ndat
             do ipw=1,npw
               gcwavef2(:,ipw+(idat-1)*npw,1:3)=zero
               lcwavef2(:,ipw+(idat-1)*npw)  =zero
             end do
           end do
           do idir=1,3
             gp2pi1=gs_ham%gprimd(idir,1)*two_pi
             gp2pi2=gs_ham%gprimd(idir,2)*two_pi
             gp2pi3=gs_ham%gprimd(idir,3)*two_pi
             kpt_cart=gp2pi1*gs_ham%kpoint(1)+gp2pi2*gs_ham%kpoint(2)+gp2pi3*gs_ham%kpoint(3)
!            Multiplication by 2pi i (G+k)_idir for gradient
!            Multiplication by -(2pi (G+k)_idir )**2 for Laplacian
             do idat=1,ndat
               do ipw=1,npw
                 kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
                 gcwavef2(1,ipw+(idat-1)*npw,idir)= cwavef2(2,ipw+(idat-1)*npw)*kg_k_cart
                 gcwavef2(2,ipw+(idat-1)*npw,idir)=-cwavef2(1,ipw+(idat-1)*npw)*kg_k_cart
                 lcwavef2(1,ipw+(idat-1)*npw)=lcwavef2(1,ipw+(idat-1)*npw)-cwavef2(1,ipw+(idat-1)*npw)*kg_k_cart**2
                 lcwavef2(2,ipw+(idat-1)*npw)=lcwavef2(2,ipw+(idat-1)*npw)-cwavef2(2,ipw+(idat-1)*npw)*kg_k_cart**2
               end do
             end do
           end do ! idir
!          STEP2: Compute (vxctaulocal)*(Laplacian of cwavef) and add it to ghc
           ABI_ALLOCATE(ghc2,(2,npw*ndat))
           call fourwf(1,vxctaulocal(:,:,:,:,1),lcwavef2,ghc2,work,gs_ham%gbound,gs_ham%gbound,&
&           gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&           npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&           use_gpu_cuda=gs_ham%use_gpu_cuda)
!$OMP PARALLEL DO
           do idat=1,ndat
             do ipw=1,npw
               ghc(:,ipw+(idat-1)*nspinor*npw)=ghc(:,ipw+(idat-1)*nspinor*npw) - half*ghc2(:,ipw+(idat-1)*npw)
             end do
           end do
           ABI_DEALLOCATE(ghc2)
           ABI_DEALLOCATE(lcwavef2)
!          STEP3: Compute sum of (grad components of vxctaulocal)*(grad components of cwavef)
           ABI_ALLOCATE(ghc2,(2,npw*ndat))
           do idir=1,3
             call fourwf(1,vxctaulocal(:,:,:,:,1+idir),gcwavef2(:,:,idir),ghc2,work,gs_ham%gbound,gs_ham%gbound,&
&             gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&             npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&             use_gpu_cuda=gs_ham%use_gpu_cuda)
!$OMP PARALLEL DO
             do idat=1,ndat
               do ipw=1,npw
                 ghc(:,ipw+(idat-1)*nspinor*npw) = ghc(:,ipw+(idat-1)*nspinor*npw) - half*ghc2(:,ipw+(idat-1)*npw)
               end do
             end do
           end do ! idir
           ABI_DEALLOCATE(ghc2)
           ABI_DEALLOCATE(gcwavef2)
         end if ! if present(vxctaulocal) i.e. metaGGA

       end if ! spin 2 treated by this proc
     end if ! npsinortot

!    Treat non-collinear local potentials
   else if (gs_ham%nvloc==4) then
     ABI_ALLOCATE(ghc1,(2,npw*ndat))
     ABI_ALLOCATE(ghc2,(2,npw*ndat))
     ABI_ALLOCATE(ghc3,(2,npw*ndat))
     ABI_ALLOCATE(ghc4,(2,npw*ndat))
     ghc1(:,:)=zero; ghc2(:,:)=zero; ghc3(:,:)=zero ;  ghc4(:,:)=zero
     ABI_ALLOCATE(vlocal_tmp,(gs_ham%n4,gs_ham%n5,gs_ham%n6))
!    ghc1=v11*phi1
     vlocal_tmp(:,:,:)=vlocal(:,:,:,1)
     if (nspinor1TreatedByThisProc) then
       call fourwf(1,vlocal_tmp,cwavef1,ghc1,work,gs_ham%gbound,gs_ham%gbound,&
&       gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
!    ghc2=v22*phi2
     vlocal_tmp(:,:,:)=vlocal(:,:,:,2)
     if (nspinor2TreatedByThisProc) then
       call fourwf(1,vlocal_tmp,cwavef2,ghc2,work,gs_ham%gbound,gs_ham%gbound,&
&       gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
     ABI_DEALLOCATE(vlocal_tmp)
     cplex=2
     ABI_ALLOCATE(vlocal_tmp,(cplex*gs_ham%n4,gs_ham%n5,gs_ham%n6))
!    ghc3=(re(v12)-im(v12))*phi1
     do i3=1,gs_ham%n6
       do i2=1,gs_ham%n5
         do i1=1,gs_ham%n4
           vlocal_tmp(2*i1-1,i2,i3)= vlocal(i1,i2,i3,3)
           vlocal_tmp(2*i1  ,i2,i3)=-vlocal(i1,i2,i3,4)
         end do
       end do
     end do
     if (nspinor1TreatedByThisProc) then
       call fourwf(cplex,vlocal_tmp,cwavef1,ghc3,work,gs_ham%gbound,gs_ham%gbound,&
&       gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
!    ghc4=(re(v12)+im(v12))*phi2
     if (nspinor2TreatedByThisProc) then
       do i3=1,gs_ham%n6
         do i2=1,gs_ham%n5
           do i1=1,gs_ham%n4
             vlocal_tmp(2*i1,i2,i3)=-vlocal_tmp(2*i1,i2,i3)
           end do
         end do
       end do
       call fourwf(cplex,vlocal_tmp,cwavef2,ghc4,work,gs_ham%gbound,gs_ham%gbound,&
&       gs_ham%istwf_k,kg_k,kg_k,gs_ham%mgfft,mpi_enreg,ndat,gs_ham%ngfft,&
&       npw,npw,gs_ham%n4,gs_ham%n5,gs_ham%n6,2,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
     ABI_DEALLOCATE(vlocal_tmp)
!    Build ghc from pieces
!    (v11,v22,Re(v12)+iIm(v12);Re(v12)-iIm(v12))(psi1;psi2): matrix product
     if (mpi_enreg%paral_spinor==0) then
       do idat=1,ndat
         do ipw=1,npw
           ghc(1:2,ipw+(idat-1)*nspinor*npw)    =ghc1(1:2,ipw+(idat-1)*npw)+ghc4(1:2,ipw+(idat-1)*npw)
           ghc(1:2,ipw+(idat-1)*nspinor*npw+npw)=ghc3(1:2,ipw+(idat-1)*npw)+ghc2(1:2,ipw+(idat-1)*npw)
         end do
       end do
     else
       call xmpi_sum(ghc4,mpi_enreg%comm_spinor,ierr)
       call xmpi_sum(ghc3,mpi_enreg%comm_spinor,ierr)
       if (nspinor1TreatedByThisProc) then
         do idat=1,ndat
           do ipw=1,npw
             ghc(1:2,ipw+(idat-1)*nspinor*npw)    =ghc1(1:2,ipw+(idat-1)*npw)+ghc4(1:2,ipw+(idat-1)*npw)
           end do
         end do
       else if (nspinor2TreatedByThisProc) then
         do idat=1,ndat
           do ipw=1,npw
             ghc(1:2,ipw+(idat-1)*nspinor*npw)=ghc3(1:2,ipw+(idat-1)*npw)+ghc2(1:2,ipw+(idat-1)*npw)
           end do
         end do
       end if
     end if
     ABI_DEALLOCATE(ghc1)
     ABI_DEALLOCATE(ghc2)
     ABI_DEALLOCATE(ghc3)
     ABI_DEALLOCATE(ghc4)
   end if ! nvloc

   if (nspinortot==2)  then
     ABI_DEALLOCATE(cwavef1)
     ABI_DEALLOCATE(cwavef2)
   end if
   ABI_DEALLOCATE(work)
   if(prtvol<0)then
     call status(0,filstat,iexit,level,'call nonlop   ')
   end if

!  Retrieve eventually original FFT distrib
   if(have_to_reequilibrate) then
     if(ndat > 1 ) then
       do idat=1,ndat
         do ipw = 1 ,npw_fft
           cwavef_fft_tr(1:2,  idat + ndat*(ipw-1)) = cwavef_fft(1:2,  ipw + npw_fft*(idat-1))
         end do
       end do
       call xmpi_alltoallv(cwavef_fft_tr,2*ndat*recvcount_fft, 2*ndat*recvdisp_fft, &
&       buff_wf,2*ndat*sendcount_fft,2*ndat*senddisp_fft, mpi_enreg%comm_fft,ierr)
     else
       call xmpi_alltoallv(cwavef_fft,2*recvcount_fft, 2*recvdisp_fft, &
&       buff_wf,2*sendcount_fft,2*senddisp_fft, mpi_enreg%comm_fft,ierr)
     end if
     do idat=1,ndat
       do ipw = 1 ,npw
         ghc(1:2,ipw + npw*(idat-1)) = buff_wf(1:2, idat + ndat*(indices_pw_fft(ipw)-1))
       end do
     end do
     ABI_DEALLOCATE(buff_wf)
     ABI_DEALLOCATE(cwavef_fft)
     if(ndat > 1) then
       ABI_DEALLOCATE(cwavef_fft_tr)
     end if
   end if

 end if ! type_calc

 if ((type_calc==0).or.(type_calc==2).or.(type_calc==3)) then

   if ((type_calc==0).or.(type_calc==2)) then
     signs=2 ; choice=1 ; nnlout=1 ; idir=0 ; tim_nonlop=1 ; nkpg=0
     cpopt_here=-1;if (gs_ham%usepaw==1) cpopt_here=cpopt
     paw_opt=gs_ham%usepaw ; if (sij_opt/=0) paw_opt=sij_opt+3

     lambda_ndat = lambda

     if(gs_ham%usepaw==0)then
       call nonlop(gs_ham%atindx1,choice,cpopt_here,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&       gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&       gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&       lambda_ndat,gs_ham%lmnmax,gs_ham%matblk,gs_ham%mgfft,mpi_enreg,gs_ham%mpsang,gs_ham%mpssoang,&
&       natom,gs_ham%nattyp,ndat,gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,nspinortot,&
&       gs_ham%ntypat,0,paw_opt,gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,&
&       signs,nonlop_dum,nonlop_dum2,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef,gvnlc,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     else
       call nonlop(gs_ham%atindx1,choice,cpopt_here,cwaveprj,gs_ham%dimekb1,gs_ham%dimekb2,dimffnl,dimffnl,&
&       gs_ham%ekb,enlout,ffnl,ffnl,gs_ham%gmet,gs_ham%gprimd,idir,&
&       gs_ham%indlmn,gs_ham%istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_ham%kpoint,gs_ham%kpoint,&
&       lambda_ndat,gs_ham%lmnmax,gs_ham%matblk,gs_ham%mgfft,mpi_enreg,gs_ham%mpsang,gs_ham%mpssoang,&
&       natom,gs_ham%nattyp,ndat,gs_ham%ngfft,nkpg,nkpg,gs_ham%nloalg,nnlout,npw,npw,nspinor,nspinortot,&
&       gs_ham%ntypat,0,paw_opt,gs_ham%phkxred,gs_ham%phkxred,gs_ham%ph1d,ph3d,ph3d,&
&       signs,gs_ham%sij,gsc,tim_nonlop,gs_ham%ucvol,gs_ham%useylm,cwavef,gvnlc,&
&       use_gpu_cuda=gs_ham%use_gpu_cuda)
     end if
   end if ! end type_calc 0 or 2 for nonlop application

!  Assemble modified kinetic, local and nonlocal contributions
!  to <G|H|C(n,k)>. Take also into account build-in debugging.
   if(prtvol/=-level)then
     do idat=1,ndat
!      MG: ifort11 miscompiles collapse(3), this one seems to work.
!      !$OMP PARALLEL DO PRIVATE(igspinor) COLLAPSE(2)
       do ispinor=1,nspinor
         do ig=1,npw
           igspinor=ig+npw*(ispinor-1)+npw*nspinor*(idat-1)
           if(kinpw(ig)<huge(zero)*1.d-11)then
             ghc(re,igspinor)= ghc(re,igspinor) + kinpw(ig)*cwavef(re,igspinor) + gvnlc(re,igspinor)
             ghc(im,igspinor)= ghc(im,igspinor) + kinpw(ig)*cwavef(im,igspinor) + gvnlc(im,igspinor)
           else
             ghc(re,igspinor)=zero
             ghc(im,igspinor)=zero
             if (sij_opt==1) then
               gsc(re,igspinor)=zero
               gsc(im,igspinor)=zero
             end if
           end if
         end do ! ig
       end do ! ispinor
     end do ! idat
   else
!    Here, debugging section
     call wrtout(std_out,' getghc : components of ghc ','PERS')
     write(message,'(a)')&
&     'icp ig ispinor igspinor re/im     ghc        kinpw         cwavef      glocc        gvnlc  gsc'
     call wrtout(std_out,message,'PERS')
     do idat=1,ndat
       do ispinor=1,nspinor
         do ig=1,npw
           igspinor=ig+npw*(ispinor-1)+npw*nspinor*(idat-1)
           if(kinpw(ig)<huge(zero)*1.d-11)then
             ghcre=kinpw(ig)*cwavef(re,igspinor)+ghc(re,igspinor)+gvnlc(re,igspinor)
             ghcim=kinpw(ig)*cwavef(im,igspinor)+ghc(im,igspinor)+gvnlc(im,igspinor)
           else
             ghcre=zero
             ghcim=zero
             if (sij_opt==1) then
               gsc(re,igspinor)=zero
               gsc(im,igspinor)=zero
             end if
           end if
           iispinor=ispinor;if (mpi_enreg%paral_spinor==1) iispinor=mpi_enreg%me_spinor+1
           if (sij_opt == 1) then
             write(message,'(a,3(1x,i5),6(1x,es13.6))') '  1 ', ig, iispinor, igspinor,ghcre,&
&             kinpw(ig),cwavef(re,igspinor),ghc(re,igspinor),gvnlc(re,igspinor), gsc(re,igspinor)
             call wrtout(std_out,message,'PERS')
             write(message,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
&             kinpw(ig),cwavef(im,igspinor),ghc(im,igspinor),gvnlc(im,igspinor), gsc(im,igspinor)
             call wrtout(std_out,message,'PERS')
           else
             write(message,'(a,3(1x,i5),6(1x,es13.6))') '  1 ', ig, iispinor, igspinor,ghcre,&
&             kinpw(ig),cwavef(re,igspinor),ghc(re,igspinor),gvnlc(re,igspinor)
             call wrtout(std_out,message,'PERS')
             write(message,'(a,3(1x,i5),6(1x,es13.6))') '  2 ', ig, iispinor, igspinor,ghcim,&
&             kinpw(ig),cwavef(im,igspinor),ghc(im,igspinor),gvnlc(im,igspinor)
             call wrtout(std_out,message,'PERS')
           end if
           ghc(re,igspinor)=ghcre
           ghc(im,igspinor)=ghcim
         end do ! ig
       end do ! ispinor
     end do ! idat
   end if

   ! Calculation of the Fock exact exchange term.
   compute_fock = .False.
   if (associated(fock)) compute_fock = (fock_get_getghc_call(fock) == 1)

   if (compute_fock) then  
     if (nspinortot==1) then
!       ABI_ALLOCATE(ghc1,(2,npw))
       do idat=1,ndat
!         ghc1(:,:)=zero
         call fock_getghc(cwavef(:,1+(idat-1)*npw:idat*npw),fock,gs_ham%gbound,ghc(:,1+(idat-1)*npw:idat*npw),&
&         gs_ham%gmet,gs_ham%istwf_k,gs_ham%kpoint,kg_k,gs_ham%mgfft,mpi_enreg,gs_ham%n4,gs_ham%n5,gs_ham%n6,&
&         gs_ham%nfft,gs_ham%ngfft,npw,paral_kgb,use_gpu_cuda=gs_ham%use_gpu_cuda)
!         do ipw=1,npw
!           ghc(:,ipw+(idat-1)*npw)=ghc(:,ipw+(idat-1)*npw)+fock%alpha*ghc1(:,ipw)
!         end do ! ipw
       end do ! idat
!       ABI_DEALLOCATE(ghc1)

     else ! nspinortot==2 
!* [BEGIN This section has not been tested. 
       if (nspinor1TreatedByThisProc) then
!         ABI_ALLOCATE(ghc1,(2,npw))
         do idat=1,ndat
!           ghc1(:,:)=zero
           call fock_getghc(cwavef1(:,1+(idat-1)*npw:idat*npw),fock,gs_ham%gbound,ghc(:,1+(idat-1)*npw:idat*npw),&
&           gs_ham%gmet,gs_ham%istwf_k,gs_ham%kpoint,kg_k,gs_ham%mgfft,mpi_enreg,gs_ham%n4,gs_ham%n5,gs_ham%n6,&
&           gs_ham%nfft,gs_ham%ngfft,npw,paral_kgb,use_gpu_cuda=gs_ham%use_gpu_cuda)
!           do ipw=1,npw
!             ghc(:,ipw+(idat-1)*nspinor*npw)=ghc(:,ipw+(idat-1)*nspinor*npw) +fock%alpha*ghc1(:,ipw)
!           end do ! ipw
         end do ! idat
!         ABI_DEALLOCATE(ghc1)
       end if ! spin 1 treated by this proc

       if (nspinor2TreatedByThisProc) then
!         ABI_ALLOCATE(ghc2,(2,npw))
         do idat=1,ndat
!           ghc2(:,:)=zero
           call fock_getghc(cwavef2(:,1+(idat-1)*npw:idat*npw),fock,gs_ham%gbound,ghc(:,1+(idat-1)*npw+shift:idat*npw+shift),&
&           gs_ham%gmet,gs_ham%istwf_k,gs_ham%kpoint,kg_k,gs_ham%mgfft,mpi_enreg,gs_ham%n4,gs_ham%n5,gs_ham%n6,&
&           gs_ham%nfft,gs_ham%ngfft,npw,paral_kgb,use_gpu_cuda=gs_ham%use_gpu_cuda)
!           do ipw=1,npw
!             ghc(:,ipw+(idat-1)*nspinor*npw+shift) = ghc(:,ipw+(idat-1)*nspinor*npw+shift) +fock%alpha*ghc2(:,ipw)
!           end do ! ipw
         end do ! idat
!         ABI_DEALLOCATE(ghc2)
       end if ! spin 2 treated by this proc
!* This section has not been tested. END]
     end if ! npsinortot
!     ghc=ghc/mpi_enreg%nproc_hf
!     call xmpi_sum(ghc,mpi_enreg%comm_hf,ierr)
   end if ! compute_fock

!  Structured debugging : if prtvol=-level, stop here.
   if(prtvol==-level)then
     write(message,'(a,i0,a)')' getghc : exit prtvol=-',level,', debugging mode => stop '
     MSG_ERROR(message)
   end if

 end if ! type_calc

 call timab(200+tim_getghc,2,tsec)

end subroutine getghc
!!***
