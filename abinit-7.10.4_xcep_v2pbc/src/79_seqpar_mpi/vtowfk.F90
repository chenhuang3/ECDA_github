!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtowfk
!! NAME
!! vtowfk
!!
!! FUNCTION
!! This routine compute the partial density at a given k-point,
!! for a given spin-polarization, from a fixed Hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cgq = array that holds the WF of the nearest neighbours of
!!        the current k-point (electric field, MPI //)
!!  cpus= cpu time limit in seconds
!!  dimcprj(natom*usepaw)=array of dimensions of array cprj (ordered by atom-type)
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  fixed_occ=true if electronic occupations are fixed (occopt<3)
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  fock <type(fock_type)>= [optional] quantities to calculate Fock exact exchange
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cprj
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point
!!  iscf=(<= 0  =>non-SCF), >0 => SCF
!!  isppol isppol=1 for unpolarized, 2 for spin-polarized
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  kpg_k(npw,nkpg)= (k+G) components (only if useylm=1)
!!  mcg=second dimension of the cg array
!!  mcgq=second dimension of the cgq array (electric field, MPI //)
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkgq = second dimension of pwnsfacq
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
!!  nkpt=number of k points.
!!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!!             (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!  npw_k=number of plane waves at this k point
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  ntypat=number of types of atoms in unit cell.
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  optforces=option for the computation of forces
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  pwind(pwind_alloc,2,3)= array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc= first dimension of pwind
!!  pwnsfac(2,pwind_alloc)= phase factors for non-symmorphic translations
!!                          (see initberry.f)
!!  pwnsfacq(2,mkgq)= phase factors for the nearest neighbours of the
!!                    current k-point (electric field, MPI //)
!!  usebanfft=flag for band-fft parallelism
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  usecprj= 1 if cprj array is stored in memory
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!  vxctaulocal(n4,n5,n6,nvloc,4)= local potential corresponding to the derivative of XC energy with respect to
!!   inetic energy density, in real space, on the augmented fft grid. (optional argument).
!!   This array contains also the gradient of vxctaulocal (gvxctaulocal) in vxctaulocal(:,:,:,:,2:4).
!!  wtk=weight assigned to the k point.
!!  zshift(nband_k)=energy shifts for the squared shifted hamiltonian algorithm
!!
!! OUTPUT
!!  dphase_k(3)=change in Zak phase for the current k-point
!!  eig_k(nband_k)=array for holding eigenvalues (hartree)
!!  ek_k(nband_k)=contribution from each band to kinetic energy, at this k-point
!!  ek_k_nd(2,nband_k,nband_k*use_dmft)=contribution to kinetic energy, including non-diagonal terms, at this k-point (usefull if use_dmft)
!!  resid_k(nband_k)=residuals for each band over all k points,
!!                   BEFORE the band rotation.
!!  ==== if optforces>0 ====
!!    grnl_k(3*natom,nband_k)=nonlocal gradients, at this k-point
!!  ==== if (gs_hamk%usepaw==0) ====
!!    enl_k(nband_k)=contribution from each band to nonlocal pseudopotential part of total energy, at this k-point
!!  ==== if (gs_hamk%usepaw==1) ====
!!    cprj(natom,mcprj*usecprj)= wave functions projected with non-local projectors:
!!                               cprj(n,k,i)=<p_i|Cnk> where p_i is a non-local projector.
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=updated wavefunctions
!!  rhoaug(n4,n5,n6,nvloc)= density in electrons/bohr**3, on the augmented fft grid.
!!                    (cumulative, so input as well as output). Update only
!!                    for occopt<3 (fixed occupation numbers)
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      cgwf,chebfi,dsymm,fourwf,fxphas,lobpcgcciiwf,lobpcgiiwf,lobpcgwf
!!      meanvalue_g,nonlop,pawcprj_alloc,pawcprj_destroy,pawcprj_put
!!      prep_fourwf,prep_nonlop,pw_orthon,subdiago,timab,wrtout,xmpi_sum,zhemm
!!
!! NOTES
!!  The cprj are distributed over band and spinors processors.
!!  One processor doesn't know all the cprj.
!!  Only the mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) projectors
!!  are stored on each proc.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine vtowfk(cg,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,dtset,&
& eig_k,ek_k,ek_k_nd,enl_k,fixed_occ,ffnl,grnl_k,gs_hamk,fock,&
& ibg,icg,ikpt,iscf,isppol,kg_k,kinpw,kpg_k,mband_cprj,mcg,mcgq,mcprj,mkgq,mpi_enreg,&
& mpw,natom,nband_k,nkpg,nkpt,nnsclo_now,npw_k,npwarr,ntypat,occ_k,optforces,ph3d,prtvol,&
& pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,rhoaug,paw_dmft,usecprj,vlocal,wtk,zshift,&
& vxctaulocal) ! optional argument

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_efield
 use m_linalg_interfaces
 use m_cgtools

 use m_hamiltonian, only : gs_hamiltonian_type
 use m_paw_dmft,    only : paw_dmft_type
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_destroy, pawcprj_put
 use m_paw_dmft,    only : paw_dmft_type
 use m_fock,        only : fock_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vtowfk'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_65_nonlocal
 use interfaces_66_wfs
 use interfaces_67_common
 use interfaces_79_seqpar_mpi, except_this_one => vtowfk
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(gs_hamiltonian_type), intent(inout) :: gs_hamk
 type(fock_type),pointer, intent(inout) :: fock
 integer, intent(in) :: dimffnl,ibg,icg,ikpt,iscf,isppol,mband_cprj,mcg,mcgq,mcprj,mkgq,mpw
 integer, intent(in) :: natom,nband_k,nkpg,nkpt,nnsclo_now,npw_k,ntypat,optforces
 integer, intent(in) :: prtvol,pwind_alloc,usecprj
 logical,intent(in) :: fixed_occ
 real(dp), intent(in) :: cpus,wtk
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(MPI_type), intent(inout) :: mpi_enreg
 type(efield_type), intent(inout) :: dtefield
 integer, intent(in) :: dimcprj(natom*gs_hamk%usepaw),kg_k(3,npw_k)
 integer, intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp), intent(in) :: cgq(2,mcgq),ffnl(npw_k,dimffnl,gs_hamk%lmnmax,ntypat)
 real(dp), intent(in) :: kinpw(npw_k),kpg_k(npw_k,nkpg),occ_k(nband_k)
 real(dp), intent(inout) :: ph3d(2,npw_k,gs_hamk%matblk)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq)
 real(dp), intent(in) :: zshift(nband_k)
 real(dp), intent(inout) :: vlocal(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc)
 real(dp), intent(out) :: eig_k(nband_k),ek_k(nband_k),dphase_k(3),ek_k_nd(2,nband_k,nband_k*paw_dmft%use_dmft)
 real(dp), intent(out) :: enl_k(nband_k*(1-gs_hamk%usepaw))
 real(dp), intent(out) :: grnl_k(3*natom,nband_k*optforces)
 real(dp), intent(out) :: resid_k(nband_k)
 real(dp), intent(inout) :: cg(2,mcg),rhoaug(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc)
 real(dp), intent(inout), optional :: vxctaulocal(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc,4)
 type(pawcprj_type),intent(inout) :: cprj(natom,mcprj*usecprj)

!Local variables-------------------------------
 integer,parameter :: level=112,tim_fourwf=2,tim_nonlop_prep=11
 integer,save :: nskip=0
!     Flag use_subovl: 1 if "subovl" array is computed (see below)
!     subovl should be Identity (in that case we should use use_subovl=0)
!     But this is true only if conjugate gradient algo. converges
 integer :: use_subovl=0
 integer :: bandpp_cprj,blocksize,choice,cpopt,iband,iband1
 integer :: iblock,iblocksize,ibs,idir,ierr,igs,igsc,ii,pidx,inonsc
 integer :: iorder_cprj,ipw,ispinor,ispinor_index,istwf_k,iwavef,jj,mgsc,my_nspinor,n1,n2,n3 !kk
 integer :: nband_k_cprj,nblockbd,ndat,nkpt_max,nnlout,ortalgo
 integer :: paw_opt,quit,signs,spaceComm,tim_nonlop,wfoptalg,wfopta10
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc,usefock
 real(dp) :: ar,ar_im,eshift,lambda_k,occblock
 real(dp) :: res,residk,weight
 character(len=500) :: message
 real(dp) :: dummy(2,1),nonlop_dum(1,1),tsec(2)
 real(dp),allocatable :: cwavef(:,:),cwavef1(:,:),cwavef_x(:,:),cwavef_y(:,:),cwavefb(:,:,:)
 real(dp),allocatable :: eig_save(:),enlout(:),evec(:,:),evec_loc(:,:),gsc(:,:)
 real(dp),allocatable :: lambda_loc(:),mat_loc(:,:),mat1(:,:,:),matvnl(:,:,:)
 real(dp),allocatable :: subham(:),subovl(:),subvnl(:),totvnl(:,:),wfraug(:,:,:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

 ! chen: exit criterion 
 integer :: pp
 real    :: finish_time, start_time, eigtol, & 
            chen_residk 
 real(8) :: tmp_eig(nband_k*dtset%nsppol), & 
            tmp_doccde(dtset%mband*dtset%nkpt*dtset%nsppol), & 
            tmp_ek, old_ek1, old_ek2, & 
            tmp_spinmagntarget, & 
            tmp_nelec_spin, & 
            tmp_entropy, & 
            tmp_fermie, & 
            tmp_occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp) :: old_eig_k(nband_k), sum_eig
 real(dp) :: old_old_eig_k(nband_k)

!! chen hack 
!! we only can use wfoptalg=14 now
! if (dtset%wfoptalg/=14 .and. dtset%wfoptalg/=4) then 
!   print *,'only works with wfoptalg==14!'
!   print *,'see how to use it in detail: http://www.abinit.org/doc/helpfiles/for-v7.6/input_variables/vardev.html#wfoptalg'
!   stop
! endif 


! **********************************************************************

 DBG_ENTER("COLL")

 call timab(28,1,tsec) ! Keep track of total time spent in "vtowfk"

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,'vtowfk : enter'
   call wrtout(std_out,message,'PERS')
 end if


!=========================================================================
!============= INITIALIZATIONS AND ALLOCATIONS ===========================
!=========================================================================

 nkpt_max=50; if(xmpi_paral==1)nkpt_max=-1

 wfoptalg=dtset%wfoptalg; wfopta10=mod(wfoptalg,10)
 istwf_k=gs_hamk%istwf_k
 quit=0

!Parallelization over spinors management
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 if (mpi_enreg%paral_spinor==0) then
   ispinor_index=1
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(dtset%nspinor==2)
 else
   ispinor_index=mpi_enreg%me_spinor+1
   nspinor1TreatedByThisProc=(mpi_enreg%me_spinor==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spinor==1)
 end if

!Define several block values for band parallelization
 if (wfopta10/=4.and.wfopta10/=5.and.wfopta10 /= 1) then
   nblockbd=nband_k
   blocksize=1
 elseif (wfopta10==5) then
   nblockbd=1
   blocksize=nband_k
 else
   nblockbd=nband_k/mpi_enreg%nproc_fft
   if (nband_k/=nblockbd*mpi_enreg%nproc_fft) nblockbd=nblockbd+1
   if(mpi_enreg%paral_kgb==1) then
     nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
!    DEBUG
!    write(std_out,*) 'starting lobpcg, with nblockbd,mpi_enreg%nproc_band',nblockbd,mpi_enreg%nproc_band
!    ENDDEBUG
   end if
   blocksize=nband_k/nblockbd
 end if

!Save eshift
 if(wfoptalg==3)then
   eshift=zshift(1)
   ABI_ALLOCATE(eig_save,(nband_k))
   eig_save(:)=eshift
 end if

 n1=gs_hamk%ngfft(1); n2=gs_hamk%ngfft(2); n3=gs_hamk%ngfft(3)

 igsc=0
 mgsc=nband_k*npw_k*my_nspinor*gs_hamk%usepaw

 ABI_ALLOCATE(gsc,(2,mgsc))
 ABI_CHECK_ALLOC("out of memory in gsc")
 gsc=zero

 if(wfopta10 /= 1) then !chebfi already does this stuff inside
   ABI_ALLOCATE(evec,(2*nband_k,nband_k))
   ABI_ALLOCATE(subham,(nband_k*(nband_k+1)))

   ABI_ALLOCATE(subvnl,(0))
   ABI_ALLOCATE(totvnl,(0,0))
   if (gs_hamk%usepaw==0) then
     if (wfopta10==4) then
       ABI_DEALLOCATE(totvnl)
       if (istwf_k==1) then
         ABI_ALLOCATE(totvnl,(2*nband_k,nband_k))
       else if (istwf_k==2) then
         ABI_ALLOCATE(totvnl,(nband_k,nband_k))
       end if
     else
       ABI_DEALLOCATE(subvnl)
       ABI_ALLOCATE(subvnl,(nband_k*(nband_k+1)))
     end if
   end if

   if (use_subovl==1) then
     ABI_ALLOCATE(subovl,(nband_k*(nband_k+1)))
   else
     ABI_ALLOCATE(subovl,(0))
   end if
 end if

!Carry out UP TO dtset%nline steps, or until resid for every band is < dtset%tolwfr

 if(prtvol>2 .or. ikpt<=nkpt_max)then
   write(message,'(a,i5,2x,a,3f9.5,2x,a)')' Non-SCF iterations; kpt # ',ikpt,', k= (',gs_hamk%kpoint,'), band residuals:'
   call wrtout(std_out,message,'PERS')
 end if

!Electric field: initialize dphase_k
 dphase_k(:) = zero

!Check that fock is present if want to use fock option
 usefock = .false.; if (dtset%usefock==1 .and. associated(fock)) usefock = .true.

!=========================================================================
!==================== NON-SELF-CONSISTENT LOOP ===========================
!=========================================================================

!nnsclo_now=number of non-self-consistent loops for the current vtrial
!(often 1 for SCF calculation, =nstep for non-SCF calculations)
 call timab(39,1,tsec) ! "vtowfk (loop)"



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! chen print header .............
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! compute the electron number for this spin channel 
 if (dtset%nsppol==2) then 
   if (isppol==1) tmp_nelec_spin = (dtset%nelect + dtset%spinmagntarget)/2.d0
   if (isppol==2) tmp_nelec_spin = (dtset%nelect - dtset%spinmagntarget)/2.d0
 else 
   tmp_nelec_spin = dtset%nelect
 endif 
 if (mpi_enreg%me==0) then
    open(unit=2222,action='write',file='invKS_log',form='formatted',access='append');
    write(2222,'(a,i4)')   "dtset%nbdbuf: ",dtset%nbdbuf
    write(2222,'(a,i4)')   "dtset%conv_unocc_orb: ",dtset%conv_unocc_orb
    write(2222,'(a,f12.6)')"tmp_nelec_spin: ",tmp_nelec_spin
    write(2222,'(a,i4)') "nline: ",dtset%nline
    if (nband_k>=4) then 
      write(2222,'(a)') & 
      'cgwf_call chen_residwf      eig(1)          eig(2)          eig(3)          eig(4)         ek_band      time(sec)'
    else
      if  (nband_k==3) write(2222,'(a)')' cg_iter  residwf         eig(1)           eig(2)          eig(3)   time(sec)'
      if  (nband_k==2) write(2222,'(a)')' cg_iter  residwf         eig(1)           eig(2)          time(sec)'
      if  (nband_k==1) write(2222,'(a)')' cg_iter  residwf         eig(1)           time(sec)'
    endif
    close(2222)
 endif







! do inonsc=1,nnsclo_now
 do inonsc=1,10000

   call cpu_time(start_time)

!  This initialisation is needed for the MPI-parallelisation (gathering using sum)
   if(wfopta10 /= 1) then
     subham(:)=zero
     if (gs_hamk%usepaw==0) then
       if (wfopta10==4) then
         totvnl(:,:)=zero
       else
         subvnl(:)=zero
       end if
     end if
     if (use_subovl==1)subovl(:)=zero
   end if
   resid_k(:)=zero

!  Filter the WFs when modified kinetic energy is too large (see routine mkkin.f)
!  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(igs,iwavef)
   do ispinor=1,my_nspinor
     do iband=1,nband_k
       igs=(ispinor-1)*npw_k
       iwavef=(iband-1)*npw_k*my_nspinor+icg
       do ipw=1+igs,npw_k+igs
         if(kinpw(ipw-igs)>huge(zero)*1.d-11)then
           cg(1,ipw+iwavef)=zero
           cg(2,ipw+iwavef)=zero
         end if
       end do
     end do
   end do

!  =========================================================================
!  ============ MINIMIZATION OF BANDS: LOBPCG ==============================
!  =========================================================================
   if(wfopta10==4.or.wfopta10==5.or.wfopta10==1) then
     if (istwf_k/=1.and.istwf_k/=2) then !no way to use lobpcg
       write(message,'(3a)')&
&       'Only istwfk=1 or 2 are allowed with wfoptalg=4/14 !',ch10,&
&       'Action: put istwfk to 1 or remove k points with half integer coordinates.'
       MSG_ERROR(message)
     end if
!      if(usefock)then
!        write(message,'(3a)')&
! &     ' vtowfk : Hartree-Fock option can not be used with lobpcg.'
!      MSG_ERROR(message)
!      end if
     if (wfopta10==4) then ! This is lobpcg I
       if(present(vxctaulocal))then
         call lobpcgwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,fock,gsc,icg,igsc,ikpt,&
&         kg_k,kinpw,mcg,mgsc,mpi_enreg,natom,nband_k,nblockbd,npw_k,ph3d,prtvol,&
&         resid_k,subham,totvnl,vlocal,vxctaulocal=vxctaulocal)
       else
         call lobpcgwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,fock,gsc,icg,igsc,ikpt,&
&         kg_k,kinpw,mcg,mgsc,mpi_enreg,natom,nband_k,nblockbd,npw_k,ph3d,prtvol,&
&         resid_k,subham,totvnl,vlocal)
       end if

     else if (wfopta10==5) then ! This is lobpcg II

       if(istwf_k==2) then ! we use the symmetric version
         if(present(vxctaulocal))then
           call lobpcgIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
&           kg_k,kinpw,mcg,mgsc,mpi_enreg,natom,nband_k,nblockbd,npw_k,ph3d,prtvol,&
&           resid_k,subham,subovl,subvnl,use_subovl,vlocal,vxctaulocal=vxctaulocal)
         else
           call lobpcgIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
&           kg_k,kinpw,mcg,mgsc,mpi_enreg,natom,nband_k,nblockbd,npw_k,ph3d,prtvol,&
&           resid_k,subham,subovl,subvnl,use_subovl,vlocal)
         end if
       else if(istwf_k==1) then !we use the hermitian version
         if(present(vxctaulocal))then
           call lobpcgccIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
&           kg_k,kinpw,mcg,mgsc,mpi_enreg,natom,nband_k,nblockbd,npw_k,ph3d,prtvol,&
&           resid_k,subham,subovl,subvnl,use_subovl,vlocal,vxctaulocal=vxctaulocal)
         else
           call lobpcgccIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
&           kg_k,kinpw,mcg,mgsc,mpi_enreg,natom,nband_k,nblockbd,npw_k,ph3d,prtvol,&
&           resid_k,subham,subovl,subvnl,use_subovl,vlocal)
         end if
       end if

     else if (wfopta10 == 1) then
       if(present(vxctaulocal))then
         call chebfi(cg(:, icg+1:),dimcprj,dimffnl,dtfil,dtset,eig_k,enl_k,&
&         ffnl,gs_hamk,gsc(:, igsc+1:),ikpt,kg_k,kinpw,mpi_enreg,natom,nband_k,npw_k,my_nspinor,&
&         ph3d,prtvol,resid_k,vlocal,vxctaulocal=vxctaulocal)
       else
         call chebfi(cg(:, icg+1:),dimcprj,dimffnl,dtfil,dtset,eig_k,enl_k,&
&         ffnl,gs_hamk,gsc(:, igsc+1:),ikpt,kg_k,kinpw,mpi_enreg,natom,nband_k,npw_k,my_nspinor,&
&         ph3d,prtvol,resid_k,vlocal)
       end if
     end if

!    In case of FFT parallelism, exchange subspace arrays
     spaceComm=mpi_enreg%comm_bandspinorfft
     if(wfopta10 /= 1) then
       call xmpi_sum(subham,spaceComm,ierr)

       if (gs_hamk%usepaw==0) then
         if (wfopta10==4) then
           call xmpi_sum(totvnl,spaceComm,ierr)
         else
           call xmpi_sum(subvnl,spaceComm,ierr)
         end if
       end if
       if (use_subovl==1) call xmpi_sum(subovl,spaceComm,ierr)
     end if

!    =========================================================================
!    ======== MINIMIZATION OF BANDS: CONJUGATE GRADIENT (Teter et al.) =======
!    =========================================================================
   else
     if(present(vxctaulocal))then
       call cgwf(dtset%berryopt,cg,cgq,dtset%chkexit,cpus,dimffnl,dphase_k,dtefield,&
&       ffnl,dtfil%filnam_ds(1),dtfil%filstat,gsc,gs_hamk,fock,icg,igsc,ikpt,inonsc,&
&       isppol,kg_k,kinpw,dtset%mband,mcg,mcgq,mgsc,mkgq,mpi_enreg,&
&       mpw,natom,nband_k,dtset%nbdblock,nkpt,dtset%nline,dtset%nloalg,&
&       npw_k,npwarr,my_nspinor,dtset%nsppol,dtset%ortalg,&
&       dtset%paral_kgb,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,&
&       pwnsfacq,quit,resid_k,subham,subovl,subvnl,dtset%tolrde,dtset%tolwfr,&
&       use_subovl,vlocal,wfoptalg,zshift,vxctaulocal=vxctaulocal)
     else
       call cgwf(dtset%berryopt,cg,cgq,dtset%chkexit,cpus,dimffnl,dphase_k,dtefield,&
&       ffnl,dtfil%filnam_ds(1),dtfil%filstat,gsc,gs_hamk,fock,icg,igsc,ikpt,inonsc,&
&       isppol,kg_k,kinpw,dtset%mband,mcg,mcgq,mgsc,mkgq,mpi_enreg,&
&       mpw,natom,nband_k,dtset%nbdblock,nkpt,dtset%nline,dtset%nloalg,&
&       npw_k,npwarr,my_nspinor,dtset%nsppol,dtset%ortalg,&
&       dtset%paral_kgb,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,&
&       pwnsfacq,quit,resid_k,subham,subovl,subvnl,dtset%tolrde,dtset%tolwfr,&
&       use_subovl,vlocal,wfoptalg,zshift)
     end if
   end if

!  =========================================================================
!  ===================== FIND LARGEST RESIDUAL =============================
!  =========================================================================

!  Find largest resid over bands at this k point
!  Note that this operation is done BEFORE rotation of bands :
!  it would be time-consuming to recompute the residuals after.
   residk=maxval(resid_k(1:max(1,nband_k-dtset%nbdbuf)))
 
!  Print residuals
   if(prtvol>2 .or. ikpt<=nkpt_max)then
     do ii=0,(nband_k-1)/8
       write(message,'(a,8es10.2)')' res:',(resid_k(iband),iband=1+ii*8,min(nband_k,8+ii*8))
       call wrtout(std_out,message,'PERS')
     end do
   end if

!  =========================================================================
!  ========== DIAGONALIZATION OF HAMILTONIAN IN WFs SUBSPACE ===============
!  =========================================================================

   if(.not. wfopta10 == 1) then
     call timab(585,1,tsec) !"vtowfk(subdiago)"
     call subdiago(cg,eig_k,evec,gsc,icg,igsc,istwf_k,&
&     mcg,mgsc,nband_k,npw_k,my_nspinor,dtset%paral_kgb,&
&     subham,subovl,use_subovl,gs_hamk%usepaw,mpi_enreg%me_g0)
     call timab(585,2,tsec)
   end if

!  Print energies
   if(prtvol>2 .or. ikpt<=nkpt_max)then
     do ii=0,(nband_k-1)/8
       write(message, '(a,8es10.2)' )' ene:',(eig_k(iband),iband=1+ii*8,min(nband_k,8+ii*8))
       call wrtout(std_out,message,'PERS')
     end do
   end if
   

!  THIS CHANGE OF SHIFT DOES NOT WORK WELL
!  Update zshift in the case of wfoptalg==3
!  if(wfoptalg==3 .and. inonsc/=1)then
!  do iband=1,nband_k
!  if(eig_k(iband)<eshift .and. eig_save(iband)<eshift)then
!  zshift(iband)=max(eig_k(iband),eig_save(iband))
!  end if
!  if(eig_k(iband)>eshift .and. eig_save(iband)>eshift)then
!  zshift(iband)=min(eig_k(iband),eig_save(iband))
!  end if
!  end do
!  eig_save(:)=eig_k(:)
!  end if

!  =========================================================================
!  =============== ORTHOGONALIZATION OF WFs (if needed) ====================
!  =========================================================================

!  Re-orthonormalize the wavefunctions at this k point--
!  this is redundant but is performed to combat rounding error in wavefunction orthogonality

   call timab(583,1,tsec) ! "vtowfk(pw_orthon)"
   ortalgo=mpi_enreg%paral_kgb
   if ((wfoptalg/=14 .and. wfoptalg /= 1).or.dtset%ortalg>0) then
     call pw_orthon(icg,igsc,istwf_k,mcg,mgsc,npw_k*my_nspinor,nband_k,ortalgo,gsc,gs_hamk%usepaw,cg,&
&     mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)
   end if
   call timab(583,2,tsec)

!  DEBUG seq==par comment next block
!  Fix phases of all bands
   if ((xmpi_paral/=1).or.(mpi_enreg%paral_kgb/=1)) then
     call fxphas(cg,gsc,icg,igsc,istwf_k,mcg,mgsc,mpi_enreg,nband_k,npw_k*my_nspinor,gs_hamk%usepaw)
   end if

!!   if (residk<dtset%tolwfr) exit  !  Exit loop over inonsc if converged



! ====================================================
! chen:  check convergence 
! ====================================================
  call fxphas(cg,gsc,icg,igsc,istwf_k,mcg,mgsc,mpi_enreg,nband_k,npw_k*my_nspinor,gs_hamk%usepaw)
  !
  ! get new occ, vtowfk only compute wfr for a given spin, 
  ! so we trick newocc by doubling the eigen array
  !
  if (dtset%nspden==2) then 
    tmp_spinmagntarget = 0.d0 
    tmp_eig(1:nband_k) = eig_k 
    tmp_eig(nband_k+1:2*nband_k) = eig_k 
    call newocc(tmp_doccde,tmp_eig,tmp_entropy,tmp_fermie,tmp_spinmagntarget,dtset%mband,dtset%nband(1),&
              tmp_nelec_spin*2.d0,dtset%nkpt,dtset%nspinor,dtset%nsppol,tmp_occ,dtset%occopt,prtvol, & 
              dtset%stmbias,dtset%tphysel,dtset%tsmear,dtset%wtk)
  else
    call newocc(tmp_doccde,eig_k,tmp_entropy,tmp_fermie,dtset%spinmagntarget,dtset%mband,dtset%nband(1),&
              dtset%nelect,dtset%nkpt,dtset%nspinor,dtset%nsppol,tmp_occ,dtset%occopt,prtvol, & 
              dtset%stmbias,dtset%tphysel,dtset%tsmear,dtset%wtk)
  endif 

  ! get kinetic energy 
  tmp_ek = 0.0d0
  do iblock=1,nblockbd
   ! Compute kinetic energy of each band
   do iblocksize=1,blocksize
     iband=(iblock-1)*blocksize+iblocksize
     call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),0)
!     ek_k(iband)=ar
     tmp_ek = tmp_ek + ar*tmp_occ(iband)
   end do
  end do

  ! get the highest occupied band 
  chen_residk = maxval(resid_k(1:nband_k-dtset%nbdbuf))

  ! === output residk ===
  call cpu_time(finish_time)
  if (mpi_enreg%me==0) then
    open(unit=2222,action='write',file='invKS_log',form='formatted',access='append');
    if (nband_k>=4) then 
      write(2222,'(i4,es16.4,5f16.8,f10.2)')inonsc,chen_residk,eig_k(1:4),tmp_ek,finish_time-start_time
    endif
    if (nband_k==3) write(2222,'(i4,es16.4,4f16.8,f10.2)')inonsc,residk,eig_k(1:3),tmp_ek,finish_time - start_time
    if (nband_k==2) write(2222,'(i4,es16.4,3f16.8,f10.2)')inonsc,residk,eig_k(1:2),tmp_ek,finish_time - start_time
    if (nband_k==1) write(2222,'(i4,es16.4,2f16.8,f10.2)')inonsc,residk,eig_k(1:1),tmp_ek,finish_time - start_time
    close(2222)
  endif

  ! === test residk ===
  if (dtset%conv_unocc_orb<0 .and. inonsc>=5) then 
    if (chen_residk<1e-8) then 
      if (mpi_enreg%me==0) then
        open(unit=2222,action='write',file='invKS_log',form='formatted',access='append');
        write(2222,'(a,es10.2,a)') & 
          '  chen_residk:',chen_residk,' < 1e-8, converged!'
        close(2222)
      endif
      exit
    endif
  endif
  !
  ! test kinetic energy 
  ! for calc_mode=7, we do the last run of cluster and env, and after that
  ! we will do solve_z(), we need to solve the unoccupied orbitals accurately
  !
 ! if (dtset%conv_unocc_orb<0) then  
 !   if (inonsc>=3) then 
 !     if ( abs(old_ek1-tmp_ek)<1e-4 .and. abs(old_ek2-tmp_ek)<1e-4)  then 
 !       open(unit=2222,action='write',file='invKS_log',form='formatted',access='append');
 !       write(2222,'(a)')'ek_band converged to 1e-4 in last two iterations, vtowfk() done.'
 !       close(2222)
 !       exit 
 !     endif 
 !   endif 
 ! endif 
  old_ek2 = old_ek1
  old_ek1 = tmp_ek 
  !
  ! convergence test 
  !
  if (dtset%conv_unocc_orb>0 .and. inonsc >= 5) then 
    eigtol = 1e-5   ! 0.27 meV tolerance
    if (maxval(abs(old_old_eig_k-eig_k))<eigtol .and. & 
        maxval(abs(old_eig_k-eig_k))<eigtol ) then 
      if (mpi_enreg%me==0) then
        open(unit=2222,action='write',file='invKS_log',form='formatted',access='append');
        write(2222,'(a,es10.4,a)')&  
         'max eigenvalue change smaller than ',eigtol,' in last steps. done!'
        close(2222)
      endif
      exit
    endif
  endif
  old_old_eig_k = old_eig_k
  old_eig_k = eig_k 

 end do !  End loop over inonsc (NON SELF-CONSISTENT LOOP)

















 call timab(39,2,tsec)
 call timab(30,1,tsec) ! "vtowfk  (afterloop)"

!###################################################################

!Compute kinetic energy and non-local energy for each band, and in the SCF
!case, contribution to forces, and eventually accumulate rhoaug

 ndat=1;if (mpi_enreg%paral_kgb==1) ndat=mpi_enreg%bandpp
 if(iscf>0 .and. fixed_occ)  then
   ABI_ALLOCATE(wfraug,(2,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6*ndat))
 end if

!"nonlop" routine input parameters
 nnlout=3*natom*optforces
 signs=1;idir=0
 if (gs_hamk%usepaw==0) then
   choice=1+optforces
   paw_opt=0;cpopt=-1;tim_nonlop=2
 else
   choice=2*optforces
   paw_opt=2;cpopt=0;tim_nonlop=10-8*optforces
 end if

 ABI_ALLOCATE(enlout,(nnlout*blocksize))

 if (wfopta10==5) then
   nblockbd=nband_k; blocksize=1
 end if

!Allocation of memory space for one WF
 ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
 if (gs_hamk%usepaw==1.and.iscf>0) then
   iorder_cprj=0
   nband_k_cprj=nband_k*(mband_cprj/dtset%mband)
   bandpp_cprj=1;if (mpi_enreg%paral_kgb==1) bandpp_cprj=mpi_enreg%bandpp
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,my_nspinor*bandpp_cprj))
   call pawcprj_alloc(cwaveprj,0,dimcprj)
 else
   ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
 end if

#undef DEV_NEW_CODE
!#define DEV_NEW_CODE

!The code below is more efficient if paral_kgb==1 (less MPI communications)
!however OMP is not compatible with paral_kgb since we should define
!which threads performs the call to MPI_ALL_REDUCE.
!This problem can be easily solved by removing MPI_enreg from meanvalue_g so that
!the MPI call is done only once outside the OMP parallel region.

#ifdef DEV_NEW_CODE
!Loop over bands or blocks of bands. Note that in sequential mode iblock=iband, nblockbd=nband_k and blocksize=1
!!$OMP PARALLEL DO PRIVATE(iband,ar)
 do iblock=1,nblockbd

!  Compute kinetic energy of each band
   do iblocksize=1,blocksize
     iband=(iblock-1)*blocksize+iblocksize

     call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),0)

     ek_k(iband)=ar

     if(paw_dmft%use_dmft==1) then
       do iband1=1,nband_k
         call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&         cg(:,1+(iband -1)*npw_k*my_nspinor+icg:iband *npw_k*my_nspinor+icg),&
&         cg(:,1+(iband1-1)*npw_k*my_nspinor+icg:iband1*npw_k*my_nspinor+icg),paw_dmft%use_dmft,ar_im=ar_im)
         ek_k_nd(1,iband,iband1)=ar
         ek_k_nd(2,iband,iband1)=ar_im
       end do
     end if
!    if(use_dmft==1) then
!    do iband1=1,nband_k
!    call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
!    &         cg(:,1+(iband -1)*npw_k*my_nspinor+icg:iband *npw_k*my_nspinor+icg),&
!    &         cg(:,1+(iband1-1)*npw_k*my_nspinor+icg:iband1*npw_k*my_nspinor+icg),use_dmft)
!    ek_k_nd(iband,iband1)=ar
!    end do
!    end if

   end do
 end do
!TODO: xmpi_sum is missing but I have to understand the logic used to deal with the different
!MPI options and communicators.
#endif


!Loop over bands or blocks of bands. Note that in sequential mode iblock=iband, nblockbd=nband_k and blocksize=1
 do iblock=1,nblockbd
   occblock=maxval(occ_k(1+(iblock-1)*blocksize:iblock*blocksize))
   cwavef(:,:)=cg(:,1+(iblock-1)*npw_k*my_nspinor*blocksize+icg:iblock*npw_k*my_nspinor*blocksize+icg)

#ifndef DEV_NEW_CODE
!  Compute kinetic energy of each band
   do iblocksize=1,blocksize
     iband=(iblock-1)*blocksize+iblocksize

     call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),0)

     ek_k(iband)=ar

     if(paw_dmft%use_dmft==1) then
       do iband1=1,nband_k
         call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&         cg(:,1+(iband -1)*npw_k*my_nspinor+icg:iband *npw_k*my_nspinor+icg),&
&         cg(:,1+(iband1-1)*npw_k*my_nspinor+icg:iband1*npw_k*my_nspinor+icg),paw_dmft%use_dmft,ar_im=ar_im)
         ek_k_nd(1,iband,iband1)=ar
         ek_k_nd(2,iband,iband1)=ar_im
       end do
     end if
   end do
#endif

   if(iscf>0)then ! In case of fixed occupation numbers, accumulates the partial density
     if (fixed_occ .and. mpi_enreg%paral_kgb/=1) then
       if (abs(occ_k(iblock))>=tol8) then
         weight=occ_k(iblock)*wtk/gs_hamk%ucvol
!        Accumulate charge density in real space in array rhoaug

!        The same section of code is also found in mkrho.F90 : should be rationalized !
         call fourwf(1,rhoaug(:,:,:,1),cwavef,dummy,wfraug,&
&         gs_hamk%gbound,gs_hamk%gbound,&
&         istwf_k,kg_k,kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&         gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,&
&         dtset%paral_kgb,tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

         if(dtset%nspinor==2)then
           ABI_ALLOCATE(cwavef1,(2,npw_k))
           cwavef1(:,:)=cwavef(:,1+npw_k:2*npw_k)

           if(dtset%nspden==1) then

             call fourwf(1,rhoaug(:,:,:,1),cwavef1,dummy,wfraug,&
&             gs_hamk%gbound,gs_hamk%gbound,&
&             istwf_k,kg_k,kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&             gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,&
&             dtset%paral_kgb,tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

           else if(dtset%nspden==4) then

!            Build the four components of rho. We use only norm quantities and, so fourwf.
!$\sum_{n} f_n \Psi^{* \alpha}_n \Psi^{\alpha}_n =\rho^{\alpha \alpha}$
!$\sum_{n} f_n (\Psi^{1}+\Psi^{2})^*_n (\Psi^{1}+\Psi^{2})_n=rho+m_x$
!$\sum_{n} f_n (\Psi^{1}-i \Psi^{2})^*_n (\Psi^{1}-i \Psi^{2})_n=rho+m_y$
             ABI_ALLOCATE(cwavef_x,(2,npw_k))
             ABI_ALLOCATE(cwavef_y,(2,npw_k))
!$(\Psi^{1}+\Psi^{2})$
             cwavef_x(:,:)=cwavef(:,1:npw_k)+cwavef1(:,1:npw_k)
!$(\Psi^{1}-i \Psi^{2})$
             cwavef_y(1,:)=cwavef(1,1:npw_k)+cwavef1(2,1:npw_k)
             cwavef_y(2,:)=cwavef(2,1:npw_k)-cwavef1(1,1:npw_k)

             call fourwf(1,rhoaug(:,:,:,4),cwavef1,dummy,wfraug,gs_hamk%gbound,gs_hamk%gbound,&
&             istwf_k,kg_k,kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&             gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,dtset%paral_kgb,&
&             tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

             call fourwf(1,rhoaug(:,:,:,2),cwavef_x,dummy,wfraug,gs_hamk%gbound,gs_hamk%gbound,&
&             istwf_k,kg_k,kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&             gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,dtset%paral_kgb,&
&             tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

             call fourwf(1,rhoaug(:,:,:,3),cwavef_y,dummy,wfraug,gs_hamk%gbound,gs_hamk%gbound,&
&             istwf_k,kg_k,kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&             gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,dtset%paral_kgb,&
&             tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

             ABI_DEALLOCATE(cwavef_x)
             ABI_DEALLOCATE(cwavef_y)

           end if ! dtset%nspden/=4
           ABI_DEALLOCATE(cwavef1)
         end if
       else
         nskip=nskip+1
       end if

!      In case of fixed occupation numbers,in bandFFT mode accumulates the partial density
     else if (fixed_occ .and. mpi_enreg%paral_kgb==1) then

       if (dtset%nspinor==1) then
         call timab(537,1,tsec) ! "prep_fourwf%vtow"
         call prep_fourwf(rhoaug(:,:,:,1),blocksize,cwavef,wfraug,iblock,ikpt,istwf_k,&
&         gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,npw_k,&
&         gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,&
&         1,gs_hamk%ucvol,wtk,use_gpu_cuda=dtset%use_gpu_cuda)
         call timab(537,2,tsec)
       else if (dtset%nspinor==2) then
         ABI_ALLOCATE(cwavefb,(2,npw_k*blocksize,2))
         ibs=(iblock-1)*npw_k*my_nspinor*blocksize+icg
!        --- No parallelization over spinors ---
         if (mpi_enreg%paral_spinor==0) then
           do iband=1,blocksize
             cwavefb(:,(iband-1)*npw_k+1:iband*npw_k,1)=cg(:,1+(2*iband-2)*npw_k+ibs:(iband*2-1)*npw_k+ibs)
             cwavefb(:,(iband-1)*npw_k+1:iband*npw_k,2)=cg(:,1+(2*iband-1)*npw_k+ibs:iband*2*npw_k+ibs)
           end do
         else
!          --- Parallelization over spinors ---
!          (split the work between 2 procs)
           cwavefb(:,:,3-ispinor_index)=zero
           do iband=1,blocksize
             cwavefb(:,(iband-1)*npw_k+1:iband*npw_k,ispinor_index) = cg(:,1+(iband-1)*npw_k+ibs:iband*npw_k+ibs)
           end do
           call xmpi_sum(cwavefb,mpi_enreg%comm_spinor,ierr)
         end if

         call timab(537,1,tsec) !"prep_fourwf%vtow"
         if (nspinor1TreatedByThisProc) then
           call prep_fourwf(rhoaug(:,:,:,1),blocksize,cwavefb(:,:,1),wfraug,iblock,&
&           ikpt,istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,npw_k,&
&           gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,gs_hamk%ucvol,wtk,&
&           use_gpu_cuda=dtset%use_gpu_cuda)
         end if
         if(dtset%nspden==1) then
           if (nspinor2TreatedByThisProc) then
             call prep_fourwf(rhoaug(:,:,:,1),blocksize,cwavefb(:,:,2),wfraug,&
&             iblock,ikpt,istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,&
&             gs_hamk%ngfft,npw_k,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,&
&             gs_hamk%ucvol,wtk,use_gpu_cuda=dtset%use_gpu_cuda)
           end if
         else if(dtset%nspden==4) then
           ABI_ALLOCATE(cwavef_x,(2,npw_k*blocksize))
           ABI_ALLOCATE(cwavef_y,(2,npw_k*blocksize))
           cwavef_x(:,:)=cwavefb(:,1:npw_k*blocksize,1)+cwavefb(:,:,2)
           cwavef_y(1,:)=cwavefb(1,1:npw_k*blocksize,1)+cwavefb(2,:,2)
           cwavef_y(2,:)=cwavefb(2,:,1)-cwavefb(1,:,2)
           if (nspinor1TreatedByThisProc) then
             call prep_fourwf(rhoaug(:,:,:,4),blocksize,cwavefb(:,:,2),wfraug,&
&             iblock,ikpt,istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,&
&             npw_k,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,gs_hamk%ucvol,wtk,use_gpu_cuda=dtset%use_gpu_cuda)
           end if
           if (nspinor2TreatedByThisProc) then
             call prep_fourwf(rhoaug(:,:,:,2),blocksize,cwavef_x,wfraug,&
&             iblock,ikpt,istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,&
&             npw_k,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,gs_hamk%ucvol,wtk,&
&             use_gpu_cuda=dtset%use_gpu_cuda)
             call prep_fourwf(rhoaug(:,:,:,3),blocksize,cwavef_y,wfraug,&
&             iblock,ikpt,istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,&
&             npw_k,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,gs_hamk%ucvol,wtk,&
&             use_gpu_cuda=dtset%use_gpu_cuda)
           end if
           ABI_DEALLOCATE(cwavef_x)
           ABI_DEALLOCATE(cwavef_y)
         end if
         call timab(537,2,tsec)
         ABI_DEALLOCATE(cwavefb)
       end if
     end if

!    Call to nonlocal operator:
!    - Compute nonlocal forces from most recent wfs
!    - PAW: compute contribution to augmentation occ. (rhoij)
     if (gs_hamk%usepaw==1.or.optforces/=0) then
!      Treat all wavefunctions in case of varying occupation numbers or PAW
!      Only treat occupied bands in case of fixed occupation numbers and NCPP
       if(fixed_occ.and.abs(occblock)<=tol8.and.gs_hamk%usepaw==0) then
         if (optforces>0) grnl_k(:,(iblock-1)*blocksize+1:iblock*blocksize)=zero
       else
         if(gs_hamk%usepaw==1) then
           call timab(554,1,tsec)  ! "vtowfk:rhoij"
         end if

         if (mpi_enreg%paral_kgb==1) then
           ABI_ALLOCATE(lambda_loc,(blocksize))
           call timab(572,1,tsec) ! 'prep_nonlop%vtowfk'
           lambda_loc(1:blocksize)=eig_k(1+(iblock-1)*blocksize:iblock*blocksize)

           call prep_nonlop(gs_hamk%atindx1,choice,cpopt,cwaveprj,gs_hamk%dimekb1,&
&           gs_hamk%dimekb2,dimffnl,gs_hamk%ekb,enlout,&
&           gs_hamk%gmet,gs_hamk%gprimd,idir,ikpt,gs_hamk%indlmn,istwf_k,&
&           gs_hamk%kpoint,lambda_loc,gs_hamk%lmnmax,gs_hamk%matblk,blocksize,gs_hamk%mgfft,&
&           mpi_enreg,gs_hamk%mpsang,gs_hamk%mpssoang,&
&           natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,gs_hamk%nloalg,nnlout,npw_k,&
&           my_nspinor,dtset%nspinor,ntypat,paw_opt,gs_hamk%phkxred,gs_hamk%ph1d,signs,&
&           gs_hamk%sij,nonlop_dum,tim_nonlop_prep,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,&
&           use_gpu_cuda=dtset%use_gpu_cuda)

           call timab(572,2,tsec)
           ABI_DEALLOCATE(lambda_loc)
         else
           lambda_k=eig_k(iblock)

           call nonlop(gs_hamk%atindx1,choice,cpopt,cwaveprj,&
&           gs_hamk%dimekb1,gs_hamk%dimekb2,dimffnl,dimffnl,gs_hamk%ekb,&
&           enlout,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,&
&           gs_hamk%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,gs_hamk%kpoint,gs_hamk%kpoint,&
&           (/lambda_k/),gs_hamk%lmnmax,gs_hamk%matblk,gs_hamk%mgfft,mpi_enreg,gs_hamk%mpsang,&
&           gs_hamk%mpssoang,natom,gs_hamk%nattyp,1,&
&           gs_hamk%ngfft,nkpg,nkpg,gs_hamk%nloalg,nnlout,npw_k,npw_k,my_nspinor,dtset%nspinor, &
&           ntypat,0,paw_opt,&
&           gs_hamk%phkxred,gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,&
&           signs,gs_hamk%sij,nonlop_dum,&
&           tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef,&
&           use_gpu_cuda=dtset%use_gpu_cuda)
         end if

         if(gs_hamk%usepaw==1) then
           call timab(554,2,tsec)
         end if

         if (optforces>0) then
           iband=(iblock-1)*blocksize
           do iblocksize=1,blocksize
             iband=iband+1;ibs=nnlout*(iblocksize-1)
             grnl_k(1:nnlout,iband)=enlout(ibs+1:ibs+nnlout)
           end do
         end if
         if (gs_hamk%usepaw==1.and.usecprj==1) then
           iband=1+(iblock-1)*bandpp_cprj
           call pawcprj_put(gs_hamk%atindx,cwaveprj,cprj,natom,iband,ibg,ikpt,iorder_cprj,isppol,&
&           mband_cprj,dtset%mkmem,natom,bandpp_cprj,nband_k_cprj,dimcprj,my_nspinor,dtset%nsppol,dtfil%unpaw,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         end if
       end if
     end if

   end if ! End of SCF calculation
 end do !  End of loop on blocks

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(enlout)

 if (gs_hamk%usepaw==1.and.iscf>0) then
   call pawcprj_destroy(cwaveprj)
 end if
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

 if (fixed_occ.and.iscf>0) then
   ABI_DEALLOCATE(wfraug)
 end if

!Write the number of one-way 3D ffts skipped until now (in case of fixed occupation numbers
 if(iscf>0 .and. fixed_occ .and. (prtvol>2 .or. ikpt<=nkpt_max) )then
   write(message,'(a,i0)')' vtowfk : number of one-way 3D ffts skipped in vtowfk until now =',nskip
   call wrtout(std_out,message,'PERS')
 end if

!Norm-conserving only: Compute nonlocal part of total energy : rotate subvnl
 if (gs_hamk%usepaw==0 .and. wfopta10 /= 1) then
   call timab(586,1,tsec)   ! 'vtowfk(nonlocalpart)'
   ABI_ALLOCATE(matvnl,(2,nband_k,nband_k))
   ABI_ALLOCATE(mat1,(2,nband_k,nband_k))
   mat1=zero

   if (wfopta10==4) then
     enl_k(1:nband_k)=zero

     if (istwf_k==1) then
       call zhemm('l','l',nband_k,nband_k,cone,totvnl,nband_k,evec,nband_k,czero,mat1,nband_k)
       do iband=1,nband_k
         res = cg_real_zdotc(nband_k,evec(:,iband),mat1(:,:,iband))
         enl_k(iband)= res
       end do
     else if (istwf_k==2) then
       ABI_ALLOCATE(evec_loc,(nband_k,nband_k))
       ABI_ALLOCATE(mat_loc,(nband_k,nband_k))
       do iband=1,nband_k
         do jj=1,nband_k
           evec_loc(iband,jj)=evec(2*iband-1,jj)
         end do
       end do
       call dsymm('l','l',nband_k,nband_k,one,totvnl,nband_k,evec_loc,nband_k,zero,mat_loc,nband_k)
       do iband=1,nband_k
         enl_k(iband)=ddot(nband_k,evec_loc(:,iband),1,mat_loc(:,iband),1)
       end do
       ABI_DEALLOCATE(evec_loc)
       ABI_DEALLOCATE(mat_loc)
     end if

   else
!    MG: This version is much faster with good OMP scalability.
!    Construct upper triangle of matvnl from subvnl using full storage mode.
     pidx=0
     do jj=1,nband_k
       do ii=1,jj
         pidx=pidx+1
         matvnl(1,ii,jj)=subvnl(2*pidx-1)
         matvnl(2,ii,jj)=subvnl(2*pidx  )
       end do
     end do

     call zhemm('L','U',nband_k,nband_k,cone,matvnl,nband_k,evec,nband_k,czero,mat1,nband_k)

!$OMP PARALLEL DO PRIVATE(res)
     do iband=1,nband_k
       res = cg_real_zdotc(nband_k,evec(:,iband),mat1(:,:,iband))
       enl_k(iband) = res
     end do
   end if

   ABI_DEALLOCATE(matvnl)
   ABI_DEALLOCATE(mat1)
   call timab(586,2,tsec)
 end if

!###################################################################

 if (iscf<=0 .and. residk>dtset%tolwfr) then
   write(message,'(a,2i5,a,es13.5)')&
&   'Wavefunctions not converged for nnsclo,ikpt=',nnsclo_now,ikpt,' max resid=',residk
   MSG_WARNING(message)
 end if

!Print out eigenvalues (hartree)
 if (prtvol>2 .or. ikpt<=nkpt_max) then
   write(message, '(5x,a,i5,2x,a,a,a,i4,a,i4,a)' ) &
&   'eigenvalues (hartree) for',nband_k,'bands',ch10,&
&   '              after ',inonsc,' non-SCF iterations with ',dtset%nline,' CG line minimizations'
   call wrtout(std_out,message,'PERS')
   do ii=0,(nband_k-1)/6
     write(message, '(1p,6e12.4)' ) (eig_k(iband),iband=1+6*ii,min(6+6*ii,nband_k))
     call wrtout(std_out,message,'PERS')
   end do
 else if(ikpt==nkpt_max+1)then
   call wrtout(std_out,' vtowfk : prtvol=0 or 1, do not print more k-points.','PERS')
 end if

!Print out decomposition of eigenvalues in the non-selfconsistent case or if prtvol>=10
 if( (iscf<0 .and. (prtvol>2 .or. ikpt<=nkpt_max)) .or. prtvol>=10)then
   write(message, '(5x,a,i5,2x,a,a,a,i4,a,i4,a)' ) &
&   ' mean kinetic energy (hartree) for ',nband_k,' bands',ch10,&
&   '              after ',inonsc,' non-SCF iterations with ',dtset%nline,' CG line minimizations'
   call wrtout(std_out,message,'PERS')

   do ii=0,(nband_k-1)/6
     write(message, '(1p,6e12.4)' ) (ek_k(iband),iband=1+6*ii,min(6+6*ii,nband_k))
     call wrtout(std_out,message,'PERS')
   end do

   if (gs_hamk%usepaw==0) then
     write(message, '(5x,a,i5,2x,a,a,a,i4,a,i4,a)' ) &
&     ' mean non-local energy (hartree) for ',nband_k,' bands',ch10,&
&     '              after ',inonsc,' non-SCF iterations with ',dtset%nline,' CG line minimizations'
     call wrtout(std_out,message,'PERS')

     do ii=0,(nband_k-1)/6
       write(message,'(1p,6e12.4)') (enl_k(iband),iband=1+6*ii,min(6+6*ii,nband_k))
       call wrtout(std_out,message,'PERS')
     end do
   end if
 end if

 if(wfopta10 /= 1) then
   ABI_DEALLOCATE(evec)
   ABI_DEALLOCATE(subham)
   !if (gs_hamk%usepaw==0) then
   !if (wfopta10==4) then
   ABI_DEALLOCATE(totvnl)
   !else
   ABI_DEALLOCATE(subvnl)
   !end if
   !end if
   ABI_DEALLOCATE(subovl)
 end if
 ABI_DEALLOCATE(gsc)

 if(wfoptalg==3) then
   ABI_DEALLOCATE(eig_save)
 end if

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
   write(message,'(a,a,a,i0,a)')' vtowfk : exit ',ch10,'  prtvol=-',level,', debugging mode => stop '
   MSG_ERROR(message)
 end if

 call timab(30,2,tsec)
 call timab(28,2,tsec)

 DBG_EXIT("COLL")

end subroutine vtowfk
!!***
