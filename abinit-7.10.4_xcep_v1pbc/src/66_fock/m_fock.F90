!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fock
!! NAME
!!  m_fock
!!
!! FUNCTION
!!  This module provides the definition of 
!!  the fock_type used to store data for the calculation of Fock exact exchange term
!!  and the procedures to perform this calculation.
!!
!! COPYRIGHT
!!  Copyright (C) 2012 ABINIT group (CMartins)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_fock

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_fstrings,   only : itoa, ftoa, sjoin

 !use m_pawang
 !use m_pawtab
 !use m_pawfgr
 !use m_pawfgrtab
 !use m_pawcprj

 implicit none

 private
!!***

!!****t* m_fock/fock_type
!! NAME
!!  fock_type
!!
!! FUNCTION
!!   This object stores the occupied wavefunctions and other quantities 
!!   needed to calculate Fock exact exchange
!!
!! NOTES
!!
!! SOURCE

 type, public :: fock_type

! Integer scalars
  !integer :: mcgocc_bz,mkg_bz,mocc
  !integer :: natom,ntypat

  integer :: usepaw=0
    ! 0 if norm-conserving psps, 1 for PAW (not implemented)

  integer :: nkpt_bz
    ! Number of k-points in the BZ for Fock operator

  integer :: ikpt,isppol,ieigen
    ! data relative to the current states.

  integer :: mkpt
    ! maximum number of k-points for Fock treated by this node

  integer :: mkptband
    ! size of occupied states stored by this node.

  integer :: mband
    ! maximum number of bands

  integer :: my_nsppol
   ! my_nsppol=1 when nsppol=1 or nsppol=2 and only one spin is treated by the processor.
   ! my_nsppol=2 when nsppol=2 and no parallelization over kpt (both spins are treated by the processor).

  integer :: nsppol 
   ! Number of indipendent spin polarizations, input variable
   ! Note that this value does not take into account the MPI distribution of the wavefunctions.

  integer :: cg_typ 
    ! Option to control the application of Vx in cgwf.F90

  integer :: nnsclo_hf
    ! Number of iterations with fixed occupied states when calculating the exact exchange contribution.

  integer :: ixc
    ! XC option (abinit input variable)

  integer ABI_PRIVATE :: getghc_call_ = 1
  ! 1 if fock_getghc should be called in getghc, 0 otherwise

! Real(dp) scalars
  real(dp) :: divgq0
    ! description of divergence in |q+G|=0

  real(dp) :: gsqcut
    !  cutoff value on G**2 for sphere inside fft box.
    !   (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2)). Used in hartre

  real(dp) :: alpha
    ! hybrid mixing coefficient

! Integer arrays
  !integer :: ngfft(18)
  !  FFT mesh used for the computation of the Fock operator
  !  Note that fock%ngfft may differ from the ngfft used for the application 
  !  of the local part of the KS Hamiltonian

  integer, allocatable :: kg_bz(:,:)
    ! kg_bz(3,mpw*mkpt)
    ! G-vectors for each k-point in the BZ treate by this node

  integer, allocatable :: nbandocc_bz(:,:) 
    ! nbandocc_bz,(mkpt,my_nsppol))
    ! nb of bands at each k point

  integer, allocatable :: istwfk_bz(:) 
    ! istwfk_bz,(mkpt))
    ! storage mode of the wavefunction at each k-point

  !integer, allocatable :: npwarr_bz(:)
    ! npwarr_bz(mkpt)
    ! Number of plane-waves used for the wavefunctions at each k-point

  integer, allocatable :: calc_phase(:)
    ! calc_phase,(mkpt))
    ! 1 if a phase factor must be considered (0 otherwise) at each k point

  integer, allocatable :: timerev(:)
    ! timerev,(mkpt))
    ! 1 if time reversal symmetry must be used (0 otherwise) at each k point 

  integer, allocatable :: tab_ibg(:,:)
    ! tab_ibg,(mkpt,my_nsppol))
    ! indices of cprj(ikpt)/occ(ikpt) in the arrays cprj/occ for each k-point jkpt

  integer, allocatable :: tab_icg(:,:)
    ! tab_icg,(mkpt,my_nsppol))
    ! indices of cg(ikpt) in the arrays cg for each k-point jkpt

  integer, allocatable :: tab_ikpt(:)
    ! tab_ikpt,(mkpt))
    ! indices of k-point ikpt in IBZ which corresponds to each k-point jkpt in full BZ

  !integer, allocatable :: indkk(:,:)

  integer, allocatable   :: gbound_bz(:,:,:)
    ! gbound_bz(2*mgfft+8,2,mkpt)
    ! Tables for zero-padded FFT of wavefunctions.

  !integer :: ham_nkpt

  !real(dp),allocatable :: ham_kptns(:,:)
  ! ham_kptns(3,ham_nkpt)

  !integer, allocatable :: ham_gbound_k(:,:,:)
    ! ham_gbound_k(2*mgfft+8,2,ham_nkpt)
    ! Table used for the zero-padded FFT of the input u(g) on which Vx will be applied

! Real(dp) arrays
  real(dp), allocatable :: cwaveocc_bz(:,:,:,:,:,:)
    ! (2,n4,n5,n6,mkptband,my_nsppol))
    ! occupied states of each bands at each k point (used to construct Fock operator)

  real(dp), allocatable :: occ_bz(:,:)
    ! occ_bz(mkptband,my_nsppol))
    ! occupancy of each bands at each k point 

  real(dp), allocatable :: wtk_bz(:)
    ! wtk_bz,(mkpt))
    ! weights assigned to each k point in the BZ
    ! Caution, the definition takes into account "ucvol" !

  real(dp), allocatable :: kptns_bz(:,:)
    ! kptns_bz(3,mkpt)
    ! k-points in full BZ

  real(dp), allocatable :: phase(:,:)
    ! phase(2,mpw*mkpt))
    ! phase factor the cg array will be multiplied with at each k point

  real(dp), allocatable :: eigen_ikpt(:)
    ! eigen_ikpt,(nband))
    !  Will contain the band index of the current state
    !  if the value is 0, the Fock contribution to the eigenvalue is not calculated.

   !integer,allocatable :: fockbz_kham_tovg(:,:)
     !fockbz_ksibz2_vg(mkpt, ks_nkpt)
     !Table giving the correspondence hf_kbz - ks_ham --> entry in fock%vg_box array

   !real(dp),allocatable :: vg_box(:,:)
     ! vg_box(nfft, vgsize)
     ! Coulomb interaction in G-space on the FFT box

!* [intermediate variables for the calculation]
!  real(dp), allocatable :: cwavef_r(:,:,:,:),vlocpsi_r(:,:,:,:)
!  real(dp), allocatable :: rhog_munu(:,:)
!  real(dp), allocatable :: dummytab3(:,:,:),dummytab2(:,:)
!  real(dp), allocatable :: work_tmp3(:) 

! Pointers to PAW-types (associated only if usepaw==1)
! Note that these are references to already existing objects.
!  type(pawang_type),pointer :: pawang => null()
!  type(pawtab_type), pointer :: pawtab(:) => null()
!  type(pawfgr_type),pointer :: pawfgr => null()
!  type(pawfgrtab_type),pointer :: pawfgrtab(:) => null()

 end type fock_type

 public :: fock_init                  ! Initialize the object.
 public :: fock_set_ieigen            ! Set the value of ieigen to the value given in argument.
 public :: fock_updateikpt            ! Update the value of energies%e_xc and energies%e_xcdc with Fock contribution.
 public :: fock_destroy               ! Free memory.
 public :: fock_calc_ene              ! Calculate the Fock contribution to the total energy.
 public :: fock_update_exc            ! Update the value of energies%e_xc and energies%e_xcdc with Fock contribution.
 public :: fock_set_getghc_call       ! Enable/disable the call to fock_getghc in getghc.
 public :: fock_get_getghc_call       ! Return the value of the flag used to enable/disable the call to fock_getghc in getghc.
 public :: fock_print                 ! Print info on the object.
!!***

 ! Help functions
 public :: bare_vqg

contains 
!!***

!!****f* ABINIT/m_fock/fock_create
!! NAME
!!  fock_create
!!
!! FUNCTION
!!  Create a fock_type structure.
!!
!! INPUTS
!!  cg(2,mcg)= wavefunctions
!!  dtset <type(dataset_type)>= all input variables for this dataset
!!  gsqcut= Fourier cutoff on G^2 used to calculate charge density
!!  kg(3,mpw*mkmem)= reduced planewave coordinates.
!!  mcg= size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  npwarr_bz(nkpt)= number of planewaves in basis at this k point
!!  occ(mband*nkpt*nsppol)= occupation number for each band (often 2) at each k point
!!  ucvol= unit cell volume ($\textrm{bohr}^{3}$)
!!
!! OUTPUT
!!  none
!!
!! SIDE EFFECTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange are allocated
!!
!! NOTES
!!
!!  ############################
!!  ### Not fully tested yet ###
!!  ############################
!!
!!  The current version is restricted to the case nsym=1, nspinor=1 and mkmem/=0.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine fock_create(fock,mgfft,mpw,mkpt,mkptband,my_nsppol,n4,n5,n6,nband)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_create'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: mgfft,mpw,mkpt,mkptband,my_nsppol,n4,n5,n6,nband
 type(fock_type),pointer :: fock

!Local variables-------------------------------
!scalars
!arrays
! character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
 !write (std_out,*) ' fock_create : enter'

!* Create the array %kptns_bz = the k points in full BZ
 ABI_ALLOCATE(fock%kptns_bz,(3,mkpt))
 fock%kptns_bz=zero
!* Create the array %jstwfk = how is stored the wavefunction at each k-point
!* By default, the table is initialized to 1 (do NOT take advantage of the time-reversal symmetry)
 ABI_ALLOCATE(fock%istwfk_bz,(mkpt))
 fock%istwfk_bz=1
!* Create the array %wtk_bz = weight assigned to each k point.
 ABI_ALLOCATE(fock%wtk_bz,(mkpt))
 fock%wtk_bz=zero
!* Create the array %npwarr_bz = number of planewaves in basis at each k-point
!    ABI_ALLOCATE(fock%npwarr_bz,(mkpt))
!    fock%npwarr_bz=0
!* Create the array %kg_bz = reduced planewave coordinates at each k-point
 ABI_ALLOCATE(fock%kg_bz,(3,mpw*mkpt))
 fock%kg_bz=0
!* Create the array %gbound_bz = boundary of the basis sphere of G vectors at each k-point
 ABI_ALLOCATE(fock%gbound_bz,(2*mgfft+8,2,mkpt))
 fock%gbound_bz=0

!* Create the array %tab_ikpt = indices of k-point ikpt in IBZ which corresponds to each k-point jkpt in full BZ
 ABI_ALLOCATE(fock%tab_ikpt,(mkpt))
 fock%tab_ikpt=0
!* Create the array %tab_ibg = indices of cprj(ikpt)/occ(ikpt) in the arrays cprj/occ for each k-point jkpt
 ABI_ALLOCATE(fock%tab_ibg,(mkpt,my_nsppol))
 fock%tab_ibg=0
!* Create the array %tab_icg = indices of cg(ikpt) in the arrays cg for each k-point jkpt
 ABI_ALLOCATE(fock%tab_icg,(mkpt,my_nsppol))
 fock%tab_icg=0

!* Create the array %calc_phase = 1 if a phase factor must be considered (0 otherwise) at each k point
 ABI_ALLOCATE(fock%calc_phase,(mkpt))
 fock%calc_phase=0
!* Create the array %phase = phase factor the cg array will be multiplied with at each k point
 ABI_ALLOCATE(fock%phase,(2,mpw*mkpt))
 fock%phase=zero

!* Create the array %timerev i= 1 if time reversal symmetry must be used (0 otherwise) at each k point 
 ABI_ALLOCATE(fock%timerev,(mkpt))
 fock%timerev=0

!* Create the array %cwaveocc_bz = wavefunctions of each bands at each k point
 ABI_ALLOCATE(fock%cwaveocc_bz,(2,n4,n5,n6,mkptband,my_nsppol))
 fock%cwaveocc_bz=zero
!* Create the array %occ_bz = occupancy of each bands at each k point => will be limited to only the occupied states
 ABI_ALLOCATE(fock%occ_bz,(mkptband,my_nsppol))
 fock%occ_bz=zero
!* Create the array %nbandocc_bz = nb of bands at each k point
 ABI_ALLOCATE(fock%nbandocc_bz,(mkpt,my_nsppol))
 fock%nbandocc_bz=0

! ========================================================
! === Set all the other state-dependent fields to zero ===
! ========================================================
 fock%ikpt= 0
!* Will contain the k-point ikpt of the current state
 fock%isppol= 0
!* Will contain the spin isppol of the current state
 fock%ieigen=0
!* Will contain the band index of the current state
!* if the value is 0, the Fock contribution to the eigenvalue is not calculated.
 ABI_ALLOCATE(fock%eigen_ikpt,(nband))
 fock%eigen_ikpt=0.d0
!* Will contain the Fock contributions to the eigenvalue of the current state

! ==============================================================
! === Allocate the memory workspace (intermediate variables) ===
! ===              to perform HF calculation                 ===
! ==============================================================
!!* Initialize the intermediate variable cwavef_r
!   ABI_ALLOCATE(fock%cwavef_r,(2,n4,n5,n6))
!!* Initialize the intermediate variable rhog_munu
!   ABI_ALLOCATE(fock%rhog_munu,(2,dtset%nfft))
!!* Initialize the intermediate variable vlocpsi_r
!   ABI_ALLOCATE(fock%vlocpsi_r,(2,n4,n5,n6))

!!* Initialize the dummy variables for fourwf
!   ABI_ALLOCATE(fock%dummytab3,(n4,n5,n6))
!   ABI_ALLOCATE(fock%dummytab2,(2,1)) ! max(ndat,ndat_occ)

!!* Initialize the variables for size-change (fftpac)
!   ABI_ALLOCATE(fock%work_tmp3,(2*nfft))


 !write (std_out,*) ' fock_create : exit'

end subroutine fock_create
!!***

!!****f* ABINIT/m_fock/fock_init
!! NAME
!!  fock_init
!!
!! FUNCTION
!!  Init all scalars and pointers in the structure.
!!
!! INPUTS
!!  cg(2,mcg)= wavefunctions
!!  dtset <type(dataset_type)>= all input variables for this dataset
!!  gsqcut= Fourier cutoff on G^2 used to calculate charge density
!!  kg(3,mpw*mkmem)= reduced planewave coordinates.
!!  mcg= size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  npwarr_bz(nkpt)= number of planewaves in basis at this k point
!!  occ(mband*nkpt*nsppol)= occupation number for each band (often 2) at each k point
!!  ucvol= unit cell volume ($\textrm{bohr}^{3}$)
!!
!! OUTPUT
!!  none
!!
!! SIDE EFFECTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange are initialized
!!
!! NOTES
!!
!!  ############################
!!  ### Not fully tested yet ###
!!  ############################
!!
!!  The current version is restricted to the case nsym=1, nspinor=1 and mkmem/=0.
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

subroutine fock_init(dtset,fock,gsqcut,kg,mpi_enreg,npwarr,rprimd,ucvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_init'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_52_fft_mpi_noabirule
 use interfaces_56_recipspace
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: gsqcut,ucvol
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
 type(fock_type),intent(inout),pointer :: fock
!arrays
 integer, intent(in) :: npwarr(dtset%nkpt)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 real(dp), intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ibg,icg,ier,ik,ikg,ikpt,isppol,jkpt,jpw,jsym,mband,mkpt,mkptband
 integer :: n1,n2,n3,n4,n5,n6,nband,nkpt,nkpt_bz,nproc_hf,npwj,timrev,v1,v2,v3
 integer :: my_jkpt,jkg_this_proc,my_nsppol
 real(dp) :: dksqmax,arg
 !character(len=500) :: msg
!arrays
 integer :: indx(1),shiftg(3),symm(3,3)
 real(dp) :: gmet(3,3),gprimd(3,3),tau_nons(3),phktnons(2,1),tsec(2)
 integer,allocatable :: indkk(:,:),kg_tmp(:),my_ikgtab(:),my_ibgtab(:,:),my_icgtab(:,:)
 real(dp),allocatable :: kptns_hf(:,:), phase1d(:,:)
 
! *************************************************************************
 
 !write (std_out,*) ' fock_init : enter'

 call timab(1500,1,tsec)

 if (dtset%nspinor/=1) then
   MSG_ERROR('Hartree-Fock option can be used only with option nspinor=1.')
 end if 

! =====================================
! === Define useful local variables ===
! =====================================           

 nkpt_bz=dtset%nkpthf
 nproc_hf=mpi_enreg%nproc_hf
 mband=dtset%nbandhf
 nband=dtset%mband
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

!* Allocations
 ABI_ALLOCATE(kptns_hf,(3,nkpt_bz))
 kptns_hf=zero
 ABI_ALLOCATE(indkk,(nkpt_bz,6))
 indkk=0
 ABI_ALLOCATE(phase1d,(2,(2*n1+1)*(2*n2+1)*(2*n3+1)))
 phase1d=zero
 ABI_ALLOCATE(kg_tmp,(3*dtset%mpw))
 
!* Initialize the array tab_indikpt = indices of kg(ikpt)/cprj(ikpt)-occ(ikpt) and cg(ikpt) associated to ikpt
! ABI_ALLOCATE(tab_indikpt,(1+2*dtset%nsppol,dtset%nkpt))
! tab_indikpt=0
! ikg=0; ibg=0; icg=0

! do ikpt=1,dtset%nkpt
!    tab_indikpt(1,ikpt)=ikg
!   tab_indikpt(2,ikpt)=ibg
!   tab_indikpt(2+dtset%nsppol,ikpt)=icg
!   ikg=ikg+npwarr(ikpt)
!   ibg=ibg+dtset%nband(ikpt)
!   icg=icg+npwarr(ikpt)*dtset%nband(ikpt)
! end do
! if (dtset%nsppol==2) then
!   do ikpt=1,dtset%nkpt
!     tab_indikpt(3,ikpt)=ibg
!     tab_indikpt(3+dtset%nsppol,ikpt)=icg
!     ibg=ibg+dtset%nband(ikpt)
!     icg=icg+npwarr(ikpt)*dtset%nband(ikpt)
!   end do
! end if

!* Initialize the array my_ikgtab = shifts in arrays kg(ikg) associated to ikpt
 ABI_ALLOCATE(my_ikgtab,(dtset%nkpt))
 ikg=0
 do ikpt=1,dtset%nkpt
   nband=dtset%nband(ikpt)
   if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband,-1,mpi_enreg%me_kpt))) then
!* The point ikpt is treated on this processor.
     my_ikgtab(ikpt)=ikg
!* The array kg is distributed, the shift ikg is incremented only on this proc.
     ikg=ikg+npwarr(ikpt)
   else
     my_ikgtab(ikpt)=-1
!* Default value is -1.
   end if
 end do
 
!* Initialize the array my_ibgtab = shifts in arrays occ(ibg) associated to ikpt
!* Initialize the array my_icgtab = shifts in arrays cg(icg) associated to ikpt
 ABI_ALLOCATE(my_ibgtab,(dtset%nkpt,dtset%nsppol))
 ABI_ALLOCATE(my_icgtab,(dtset%nkpt,dtset%nsppol))
 ibg=0; icg=0
 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     nband=dtset%nband(ikpt)
     if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband,isppol,mpi_enreg%me_kpt))) then
!* The states with (ikpt,isppol) are stored on this processor.
       my_icgtab(ikpt,isppol)=icg
!* The array cg is distributed, the shift icg is incremented only on this proc.
       icg=icg+npwarr(ikpt)*nband
     else
       my_icgtab(ikpt,isppol)=-1
!* Otherwise, the states with (ikpt,isspol) are not stored on this processor and default value is -1.
     end if
!* The array occ is shared among the proc, the shift ibg is always incremented.
     my_ibgtab(ikpt,isppol)=ibg
     ibg=ibg+nband
   end do
 end do
 
 if (.not.(associated(fock))) then

! =================================
! === Create the fock structure ===
! =================================
   ABI_DATATYPE_ALLOCATE(fock,)

!* Compute the dimension of arrays in "spin" w.r.t parallelism
   my_nsppol=dtset%nsppol
   if (mpi_enreg%nproc_kpt>1) my_nsppol=1
!* my_nsppol=1 when nsppol=1 or nsppol=2 and only one spin is treated by the processor.
!* my_nsppol=2 when nsppol=2 and no parallelization over kpt (both spins are treated by the processor).

!* Compute mkpt the size of arrays/pointers for k poitns w.r.t. parallelism
!* Compute mkptband the size of arrays/pointers for occupied states w.r.t. parallelism
   if (nproc_hf<nkpt_bz) then
!* Parallelization over kpts only 
     mkpt=nkpt_bz/nproc_hf
     if (mod(nkpt_bz,nproc_hf) /=0) mkpt=mkpt+1
     mkptband=mkpt*mband
   else
!* Parallelization over occupied states
     if (nproc_hf<nkpt_bz*mband) then 
       mkptband=(nkpt_bz*mband)/nproc_hf 
       if (mod((nkpt_bz*mband),nproc_hf) /=0) mkptband=mkptband+1
       mkpt=1
       if (mod(nproc_hf,nkpt_bz) /=0) mkpt=2
     else
       mkptband=1 
       mkpt=1
     end if
   end if 
   
   call fock_create(fock,dtset%mgfft,dtset%mpw,mkpt,mkptband,my_nsppol,n4,n5,n6,nband)

!* Initialize %mabnd, %mkpt, %mkptband = size of arrays
   fock%mband=mband
   fock%mkpt=mkpt
   fock%mkptband=mkptband
   fock%my_nsppol = my_nsppol
   fock%nsppol = dtset%nsppol

! ==========================================
! === Initialize the convergence options ===
! ==========================================
!* Type of conjugate gradient algorithm used for exact exchange calculation
 if ((dtset%cgtyphf<0).or.(dtset%cgtyphf>2)) then 
    MSG_ERROR('The parameter cgtyp must have an integer value between 0 and 2.')
  end if
  if (dtset%cgtyphf==0) then 
    fock%cg_typ=2
    write (std_out,*) ' fock_init : The parameter cgtyphf is set to its default value 2.'
!* Default value is set to 2 (calculation of exact exchange each time the function getghc is called in cgwf)
!* May be useful to put default to 1 (calculation of exact exchange only for the first call to getghc in cgwf)
  else 
    fock%cg_typ=dtset%cgtyphf
    write (std_out,*) ' fock_init : The parameter cgtyphf is set to the value:', dtset%cgtyphf
!* value chosen by the user : 1 or 2.
  end if
  
!* Number of iterations with fixed occupied states when calculating the exact exchange contribution.
  if (dtset%nnsclohf<0) then 
    MSG_ERROR('The parameter nnsclohf must be a non-negative integer.')
  end if
  if (dtset%nnsclohf==0) then 
    fock%nnsclo_hf=1
    write (std_out,*) ' fock_init : The parameter nnsclohf is set to its default value 1.'
!* Default value is set to 1 (updating cgocc at each step)
!* May be useful to put default to 3 
  else 
    fock%nnsclo_hf=dtset%nnsclohf
    write (std_out,*) ' fock_init : The parameter nnsclohf is set to the value:', dtset%nnsclohf
!* value chosen by the user 
  end if

! =========================================
! === Initialize the hybrid coefficient ===
! =========================================
  fock%ixc = dtset%ixc

  if (dtset%ixc==40) then
    fock%alpha=1.d0
    write (std_out,*) ' fock_init : This is an Hartree-Fock calculation. The mixing coefficient alpha is set to 1.'
  end if
  if (dtset%ixc==41) then
    fock%alpha=1/(4.d0)
    write (std_out,*) ' fock_init : This is a standard PBE0 calculation. The mixing coefficient alpha is set to 0.25.'
  end if
  if (dtset%ixc==42) then
    fock%alpha=1/(3.d0)
    write (std_out,*) ' fock_init : This is a modified PBE0 calculation. The mixing coefficient alpha is set to 0.33.'
  end if
 
 
! ======================================================
! === Initialize the data relative to Poisson solver ===
! ======================================================
!* In the case |q+G| = 0, the density will be considered constant in the small volume around q=0 
!* and the Coulomb term is integrated out, giving the constant divq0.
   MSG_WARNING('The divergence in q=0 is treated, assuming a spherical BZ.')
!* Compute unit cell volume (useless since already calculated in scfcv)
!    ucvol=rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3))+&
! &    rprimd(2,1)*(rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3))+&
! &    rprimd(3,1)*(rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3))
!* Evaluate the divergence term at q=0)
   fock%divgq0= 7.79*(nkpt_bz*ucvol)**(two_thirds)
   write(std_out,*)'fock%divgq0:  ',fock%divgq0
!  fock%divgq0= zero
!* divgq0 = constant to put in vhartre if q=0 
   fock%gsqcut= gsqcut
!* gsqcut = cutoff value on G^2 for sphere inside the fft box (input for vhartre).

! =======================================================
! === Initialize the properties of the k-points in BZ ===
! =======================================================
!* Initialize %nkpt_bz = nb of k point in BZ for the calculation of exchange
   fock%nkpt_bz=nkpt_bz
!* Initialize the array %wtk_bz = weight assigned to each k point.
   fock%wtk_bz=1.0_dp/(dble(nkpt_bz)*ucvol)
!* Caution, the definition takes into account "ucvol" !

   if (dtset%kptopt>=1 .and. dtset%kptopt<=4) then
     if (dtset%kptopt/=3) then
! ============================================
! === Initialize the set of k-points in BZ ===
! ============================================
!* Generate all the k-points in BZ (Monkhorst-Pack grid)
!* brav=1 to treat all Bravais lattices ; iout=0 since we do not want any output ; option=0 since we consider k-points 
       call smpbz(1,0,dtset%kptrlatt,nkpt_bz,nkpt,dtset%nshiftk,0,dtset%shiftk,kptns_hf)
!* kptns_hf contains the special k points obtained by the Monkhorst & Pack method, in reduced coordinates. (output)
       if (nkpt_bz/=nkpt) then
         MSG_ERROR('The value of nkpt_bz and the result of smpbz should be equal!')
       end if 

! =======================================================
! === Compute the transformation to go from IBZ to BZ ===
! =======================================================
!* Compute the reciprocal space metric.
       call matr3inv(rprimd,gprimd)
       gmet = MATMUL(TRANSPOSE(gprimd),gprimd)

!* Calculate the array indkk which describes how to get IBZ from BZ
!* dksqmax=maximal value of the norm**2 of the difference between a kpt2 vector and the closest k-point found from the kptns1 set, using symmetries. (output)
!* sppoldbl=1, no spin-polarisation doubling is required.
       timrev=1 ; if (dtset%kptopt==4) timrev=0
!* timrev=1 if the use of time-reversal is allowed ; 0 otherwise
       if (dtset%kptopt==2) then
!* Only time reversal symmetry is used.
         symm=0 ; symm(1,1)=1 ; symm(2,2)=1 ; symm(3,3)=1
         call listkk(dksqmax,gmet,indkk(1:nkpt_bz,:),dtset%kptns,kptns_hf,dtset%nkpt, & 
&            nkpt_bz,1,1,indx,symm,timrev)
       else
!* As in getkgrid, no use of antiferromagnetic symmetries thans to the option sppoldbl=1
         call listkk(dksqmax,gmet,indkk(1:nkpt_bz,:),dtset%kptns,kptns_hf,dtset%nkpt, & 
&            nkpt_bz,dtset%nsym,1,dtset%symafm,dtset%symrel,timrev)
       end if
!* indkk(nkpt_bz,6) describes the k point of IBZ that generates each k point of BZ
!*      indkk(:,1)   = k point of IBZ, kpt_ibz
!*      indkk(:,2)   = symmetry operation to apply to kpt_ibz to give the k point of BZ
!*                     (if 0, means no symmetry operation, equivalent to identity )
!*      indkk(:,3:5) = Umklapp vectors to apply to remain in BZ
!*      indkk(:,6)   = 1 if time-reversal was used to generate the k point of BZ, 0 otherwise
!* No use of symafm to generate spin down wfs from spin up wfs for the moment

     else ! In this case, dtset%kptopt=3, one deals with BZ directly.
! ============================================
! === Initialize the set of k-points in BZ ===
! ============================================
       if (nkpt_bz/=dtset%nkpt) then
         MSG_ERROR('In this version, the value of nkpt_bz and nkpt should be equal!')
       end if 

       kptns_hf=dtset%kptns

! ==========================================================
! === Initialize the transformation to go from IBZ to BZ ===
! ==========================================================
!* indkk(nkpt_bz,6) describes the k point of IBZ that generates each k point of BZ
       do ikpt=1,nkpt_bz
         indkk(ikpt,1)=ikpt
       end do
!* indkk(:,1)   = k point of BZ and kpt=kpt_hf
!* all the other field are zero because the Identity is the only symmetry operation.

!* In the most general case, use of listkk is certainly possible.
     end if

   else 
     if (dtset%kptopt==0) then
!* kptopt =0 : read directly nkpt, kpt, kptnrm and wtk in the input file 
!*              => this case is not allowed for the moment
       MSG_ERROR('Hartree-Fock option can not be used with option kptopt=0.')
     else 
!* kptopt <0 : rely on kptbounds, and ndivk to set up a band structure calculation 
!*              => a band structure calculation is not yet allowed. 
       MSG_ERROR('Hartree-Fock option can not be used with option kptopt<0.')
     end if
   end if
       
!! =======================================================
!! === Initialize the properties of the k-points in BZ ===
!! =======================================================
!       jkg=0
!!* Initialize the arrays %npwarr_bz, %kg_j, %phase_j, %gbound_j
!       do jkpt=1,nkpt_bz
!         ikpt=indkk(jkpt,1)
!!* ikpt = the point of IBZ that jkpt is an image of in BZ
!         npwj=npwarr(ikpt)
!!* npwj = number of planewaves in basis at point jkpt = at point ikpt
!         jsym=indkk(jkpt,2)
!!* jsym = symmetry operation to apply to get jkpt from ikpt
!         shiftg(:)=indkk(jkpt,3:5)
!!* shiftg = Bravais vector G0 to add to remain in BZ
!         if (jsym/=0) then
!           symm(:,:)=dtset%symrel(:,:,jsym) 
!           tau_nons(:)=dtset%tnons(:,jsym)
!!* The symmetry operation in k-space (symm) and the non-symorphic translation (tau_nons) are now defined. 
!           if(sum(tau_nons(:)**2)>tol8) then
!!* Initialize %calc_phase(jkpt) to 1
!             fock%calc_phase(jkpt)=1
!!* Compute the phase factor exp(i*2*pi*G.tau) for all G. 
!             indx(1)=1
!             phase1d=zero
!             call getph(indx,1,n1,n2,n3,phase1d,tau_nons)
!!* Although the routine getph is orignally written for atomic phase factors, it does precisely what we want
!             arg=two_pi*(dtset%kptns(1,ikpt)*tau_nons(1) + dtset%kptns(2,ikpt)*tau_nons(2) & 
!&                + dtset%kptns(3,ikpt)*tau_nons(3))
!             phktnons(1,1)=cos(arg)
!             phktnons(2,1)=sin(arg)
!!              phktnons(1,1)=one
!!              phktnons(2,1)=zero
!!* Convert 1D phase factors to 3D phase factors exp(i*2*pi*(k+G).tau) and store it in %phase_j
!             call ph1d3d(1,1,kg(:,1+tab_indikpt(1,ikpt):npwj+tab_indikpt(1,ikpt)),1,1,npwj,n1, & 
!&              n2,n3,phktnons,phase1d,fock%phase(:,1+jkg:npwj+jkg))
!           end if
!         else
!           symm=0 ; symm(1,1)=1 ; symm(2,2)=1 ; symm(3,3)=1
!           tau_nons(:)=zero
!           shiftg(:)=0
!         end if
!!* Apply time-reversal symmetry if required
!         if(indkk(jkpt,6)/=0) then
!!* Initialize %timerev(jkpt) to 1
!           fock%timerev(jkpt)=1
!           symm(:,:)=-symm(:,:)
!         end if

!!* Initialize %istwfk_bz(jkpt) to 
!         fock%istwfk_bz(jkpt)=dtset%istwfk(ikpt)

!!* Initialize %tab_ikpt and %tab_ibgcg
!         fock%tab_ikpt(jkpt)=ikpt
!         fock%tab_ibgcg(1:dtset%nsppol,jkpt)=tab_indikpt(2:1+dtset%nsppol,ikpt)
!         fock%tab_ibgcg(1+dtset%nsppol:2*dtset%nsppol,jkpt)= & 
!&          tab_indikpt(2+dtset%nsppol:2*dtset%nsppol+1,ikpt)

!!* Initialize %npwarr_bz
!         fock%npwarr_bz(jkpt)=npwj

!!* Initialize %kg_bz
!         do jpw=1,npwj
!           v1=kg(1,jpw+tab_indikpt(1,ikpt)) ; v2=kg(2,jpw+tab_indikpt(1,ikpt)) ; v3=kg(3,jpw+tab_indikpt(1,ikpt)) 
!           fock%kg_bz(1,jpw+jkg)=-shiftg(1)+symm(1,1)*v1+symm(2,1)*v2+symm(3,1)*v3
!           fock%kg_bz(2,jpw+jkg)=-shiftg(2)+symm(1,2)*v1+symm(2,2)*v2+symm(3,2)*v3
!           fock%kg_bz(3,jpw+jkg)=-shiftg(3)+symm(1,3)*v1+symm(2,3)*v2+symm(3,3)*v3
!!* The symmetry operation symm must be transposed when used. (cf. docs about wfconv)
!         end do

!!* Initialize %gbound_bz
!         call sphereboundary(fock%gbound_bz(:,:,jkpt),fock%istwfk_bz(jkpt), & 
!&          fock%kg_bz(:,1+jkg:npwj+jkg),dtset%mgfft,npwj)

!!* Update of the shift to be applied
!         jkg=jkg+npwj 
!       end do

! ==========================================================
! === Initialize the k-points in BZ and their properties ===
! ==========================================================
!   jkg=0;
   jkg_this_proc=0;my_jkpt=0
   do jkpt=1,nkpt_bz

!* If this processor does not calculate exchange with the k point jkpt, skip the rest of the k-point loop.
     if (proc_distrb_cycle(mpi_enreg%distrb_hf,jkpt,1,mband,1,mpi_enreg%me_hf)) cycle 
!       if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,jkpt,1,dtset%nbandhf,1,mpi_enreg%me_kpt))) then
!* The processor does own a copy of the array kg of ikpt ; increment the shift.
!         jkg=jkg+npwj
!        end if
! Skip the rest of the k-point loop
!       cycle
!     end if  
     my_jkpt=my_jkpt+1

     ikpt=indkk(jkpt,1)
!* ikpt = the point of IBZ that jkpt is an image of in BZ
     npwj=npwarr(ikpt)
!* npwj = number of planewaves in basis at point jkpt = at point ikpt
     jsym=indkk(jkpt,2)
!* jsym = symmetry operation to apply to get jkpt from ikpt
     shiftg(:)=indkk(jkpt,3:5)
!* shiftg = Bravais vector G0 to add to remain in BZ

!* Initialize the array %kptns_bz = the k points in full BZ
     fock%kptns_bz(:,my_jkpt)=kptns_hf(:,jkpt)

!* Initialize the array %jstwfk = how is stored the wavefunction at each k point
     if (dtset%istwfk(ikpt)/=1) then
       fock%istwfk_bz(my_jkpt)=set_istwfk(kptns_hf(:,jkpt))
     end if
!* One can take advantage of the time-reversal symmetry in this case.
!* Initialize the array %wtk_bz = weight assigned to each k point.
!     fock%wtk_bz(my_jkpt)=dtset%wtk(jkpt)/ucvol
!* Caution, the definition takes into account "ucvol" !

!* Initialize the array %npwarr_bz = number of planewaves in basis at each k point
!     fock%npwarr_bz(my_jkpt)=npwj

!!* Initialize the array %tab_ikpt = indices of k-point in IBZ ikpt for each k point jkpt in BZ (here,ikpt=jkpt)
     fock%tab_ikpt(my_jkpt)=ikpt

!!* Initialize the array %tab_ibgcg = indices of cprj(ikpt)/occ(ikpt) and cg(ikpt) for each k point jkpt
!     if (my_nsppol==2) then
!!* In this case, my_nsppol=dtset%nsppol=2 
!       fock%tab_ibgcg(1:2,my_jkpt)=tab_indikpt(2:3,ikpt)
!       fock%tab_ibgcg(3:4,my_jkpt)=tab_indikpt(4:5,ikpt)
!     else 
!       if(mpi_enreg%my_isppoltab(1)==1) then
!!* In this case, my_nsppol=1 and the up spin is treated (dtset%nsppol= 1 or 2)
!         fock%tab_ibgcg(1,my_jkpt)=tab_indikpt(2,ikpt)
!         fock%tab_ibgcg(2,my_jkpt)=tab_indikpt(2+dtset%nsppol,ikpt)
!       else
!!* In this case, my_nsppol=1 and the dn spin is treated (so dtset%nsppol=2)
!         fock%tab_ibgcg(1,my_jkpt)=tab_indikpt(3,ikpt)
!         fock%tab_ibgcg(2,my_jkpt)=tab_indikpt(5,ikpt)
!       end if
!     end if

!* Initialize the array %kg_bz = reduced planewave coordinates at each k point
     if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,-1,mpi_enreg%me_kpt))) then
!* We perform the test with isppol=-1 (both spins) and the occupied band (dtset%nbandhf). 
!* We assume that paral_kgb==0 (a k-point may not be present on several proc.)
!* The array kg for ikpt is stored on this processor and copied in kg_tmp.
       ikg=my_ikgtab(ikpt)
!* ikg = the shift in kg to get the G-vectors associated to ikpt
       do ik=1,3
!         kg_tmp(1+(ik-1)*npwj:ik*npwj)=kg(ik,1+tab_indikpt(1,ikpt):npwj+tab_indikpt(1,ikpt))
         kg_tmp(1+(ik-1)*npwj:ik*npwj)=kg(ik,1+ikg:npwj+ikg)
       end do
!       jkg=jkg+npwj
     end if
!* Broadcast the array kg_tmp to all the processors of comm_kpt.
!* Since paral_kgb==0, all the bands of a k-point are treated on the same proc. 
     call xmpi_bcast(kg_tmp,mpi_enreg%proc_distrb(ikpt,1,1),mpi_enreg%comm_kpt,ier)
     do ik=1,3
       fock%kg_bz(ik,1+jkg_this_proc:npwj+jkg_this_proc)=kg_tmp(1+(ik-1)*npwj:ik*npwj)
     end do

!* Apply a symmetry operation on kg_bz if necessary
     if (jsym/=0) then
       symm(:,:)=dtset%symrel(:,:,jsym) 
       tau_nons(:)=dtset%tnons(:,jsym)
!* The symmetry operation in k-space (symm) and the non-symorphic translation (tau_nons) are now defined. 
       if(sum(tau_nons(:)**2)>tol8) then
!* Initialize %calc_phase(jkpt) to 1
         fock%calc_phase(my_jkpt)=1
!* Compute the phase factor exp(i*2*pi*G.tau) for all G. 
         indx(1)=1
         phase1d=zero
         call getph(indx,1,n1,n2,n3,phase1d,tau_nons)
!* Although the routine getph is orignally written for atomic phase factors, it does precisely what we want
         arg=two_pi*(dtset%kptns(1,ikpt)*tau_nons(1) + dtset%kptns(2,ikpt)*tau_nons(2) & 
&            + dtset%kptns(3,ikpt)*tau_nons(3))
         phktnons(1,1)=cos(arg)
         phktnons(2,1)=sin(arg)
!          phktnons(1,1)=one
!          phktnons(2,1)=zero
!* Convert 1D phase factors to 3D phase factors exp(i*2*pi*(k+G).tau) and store it in %phase_j
         call ph1d3d(1,1,fock%kg_bz(:,1+jkg_this_proc:npwj+jkg_this_proc),1,1,npwj,n1,n2,n3, & 
&          phktnons,phase1d,fock%phase(:,1+jkg_this_proc:npwj+jkg_this_proc))
       end if
!* Apply time-reversal symmetry if required
       if(indkk(jkpt,6)/=0) then
!* Initialize %timerev(jkpt) to 1
         fock%timerev(my_jkpt)=1
         symm(:,:)=-symm(:,:)
       end if
!* Initialize %kg_bz
       do jpw=1,npwj
         v1=fock%kg_bz(1,jpw+jkg_this_proc) ; v2=fock%kg_bz(2,jpw+jkg_this_proc) ; v3=fock%kg_bz(3,jpw+jkg_this_proc) 
         fock%kg_bz(1,jpw+jkg_this_proc)=-shiftg(1)+symm(1,1)*v1+symm(2,1)*v2+symm(3,1)*v3
         fock%kg_bz(2,jpw+jkg_this_proc)=-shiftg(2)+symm(1,2)*v1+symm(2,2)*v2+symm(3,2)*v3
         fock%kg_bz(3,jpw+jkg_this_proc)=-shiftg(3)+symm(1,3)*v1+symm(2,3)*v2+symm(3,3)*v3
!* The symmetry operation symm must be transposed when used. (cf. docs about wfconv)
       end do
     else
!* Ths symmetry operation is the identity.         
!* Apply time-reversal symmetry if required
       if(indkk(jkpt,6)/=0) then
!* Initialize %timerev(jkpt) to 1
         fock%timerev(my_jkpt)=1
         fock%kg_bz(ik,1+jkg_this_proc:npwj+jkg_this_proc)=-fock%kg_bz(ik,1+jkg_this_proc:npwj+jkg_this_proc)
       end if
     end if

!* Initialize the array %gbound_bz = boundary of the basis sphere of G vectors at each k point
     call sphereboundary(fock%gbound_bz(:,:,my_jkpt),fock%istwfk_bz(my_jkpt),& 
&      fock%kg_bz(:,1+jkg_this_proc:npwj+jkg_this_proc),dtset%mgfft,npwj)

     jkg_this_proc=jkg_this_proc+npwj
     
!* Initialize the arrays %tab_ibg = shifts in arrays cprj and occ (ibg) for each k point jkpt
!* Initialize the arrays %tab_icg = shifts in arrays cg(icg) for each k point jkpt
     if (my_nsppol==1) then
         fock%tab_ibg(my_jkpt,1)=my_ibgtab(ikpt,1+mpi_enreg%my_isppoltab(2))
         fock%tab_icg(my_jkpt,1)=my_icgtab(ikpt,1+mpi_enreg%my_isppoltab(2))
!* if mpy_isppoltab(2)=0, the up spin is treated (dtset%nsppol= 1 or 2)
!* if mpy_isppoltab(2)=1, the dn spin is treated (so dtset%nsppol=2)
     
!       if(mpi_enreg%my_isppoltab(2)==1) then
!* In this case, my_nsppol=1 and the up spin is treated (dtset%nsppol= 1 or 2)
!         fock%tab_ibg(my_jkpt,1)=my_ibgtab(ikpt,1)
!         fock%tab_icg(my_jkpt,1)=my_icgtab(ikpt,1)
!       else
!* In this case, my_nsppol=1 and the dn spin is treated (so dtset%nsppol=2)
!         fock%tab_ibg(my_jkpt,1)=my_ibgtab(ikpt,2)
!         fock%tab_icg(my_jkpt,1)=my_icgtab(ikpt,2)
!       end if
     else
!* In this case, my_nsppol=dtset%nsppol=2 
       fock%tab_ibg(my_jkpt,:)=my_ibgtab(ikpt,:)
       fock%tab_icg(my_jkpt,:)=my_icgtab(ikpt,:)
     end if

   enddo

!* Deallocation
   ABI_DEALLOCATE(indkk)
   ABI_DEALLOCATE(kg_tmp)
   ABI_DEALLOCATE(kptns_hf)
   ABI_DEALLOCATE(my_ibgtab)
   ABI_DEALLOCATE(my_icgtab)
   ABI_DEALLOCATE(my_ikgtab)
   ABI_DEALLOCATE(phase1d)
 end if

 call fock_print(fock,unit=std_out)

 call timab(1500,2,tsec)
 
 !write (std_out,*) ' fock_init : exit'

end subroutine fock_init
!!***

!!****f* ABINIT/m_fock/fock_updateikpt
!! NAME
!!  fock_updateikpt
!!
!! FUNCTION
!!  Update the value of ikpt,isppol for the next exact exchange calculation.
!!
!! INPUTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!  gs_ham <type(gs_hamiltonian_type)>= all data for the Hamiltonian to be applied
!!  ikpt= reduced planewave coordinates.
!!  isppol= number of planewaves in basis at this k point
!!
!! OUTPUT
!!  none
!!
!! SIDE EFFECTS
!!   The field fock%eigen_ikpt is also set to 0.d0.
!!
!! NOTES
!!  May be improved to calculate the star of ikpt. => I think NO finally
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!
!! SOURCE

subroutine fock_updateikpt(fock,ikpt,isppol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_updateikpt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: ikpt,isppol
 type(fock_type),pointer :: fock

! *************************************************************************
 
 !write (std_out,*) ' fock_updateikpt : enter'

! ======================================================
! === Update the data relative to the current states ===
! ======================================================
!* Copy of the value ikpt in the field ikpt
   fock%ikpt=ikpt
!* Copy of the value isppol in the field isppol
   fock%isppol=isppol
!* Set all the Fock contributions to the eigenvalues to 0.d0.
   fock%eigen_ikpt=zero

end subroutine fock_updateikpt
!!***

!!****f* ABINIT/m_fock/fock_set_ieigen
!! NAME
!!  fock_set_ieigen
!!
!! FUNCTION
!!  Set the value of ieigen to the value given in argument.
!!
!! INPUTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!  iband= index of the band iband
!!
!! OUTPUT
!!  none
!!
!! PARENTS
!!      cgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine fock_set_ieigen(fock,iband)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_set_ieigen'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: iband
 type(fock_type),pointer :: fock
 
! *************************************************************************
 
! ======================================================
! === Update the data relative to the current states ===
! ======================================================
!* Copy of the value iband in the field ieigen
 fock%ieigen=iband

end subroutine fock_set_ieigen
!!***

!!****f* ABINIT/m_fock/fock_destroy
!! NAME
!!  fock_destroy
!!
!! FUNCTION
!!  Clean and destroy fock datastructure.
!!
!! INPUTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine fock_destroy(fock)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_destroy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(fock_type),pointer :: fock

! *************************************************************************

 DBG_ENTER("COLL")
 
 ! real arrays 
 if (allocated(fock%cwaveocc_bz)) then
   ABI_DEALLOCATE(fock%cwaveocc_bz)
 endif
 if (allocated(fock%occ_bz)) then
   ABI_DEALLOCATE(fock%occ_bz)
 endif

 ! Deallocate integer arrays
 if (allocated(fock%kg_bz)) then
   ABI_DEALLOCATE(fock%kg_bz)
 endif
 if (allocated(fock%nbandocc_bz)) then
   ABI_DEALLOCATE(fock%nbandocc_bz)
 endif
 if (allocated(fock%istwfk_bz)) then
   ABI_DEALLOCATE(fock%istwfk_bz)
 endif
 !if (allocated(fock%npwarr_bz)) then
 !   ABI_DEALLOCATE(fock%npwarr_bz)
 !endif
 if (allocated(fock%calc_phase)) then
    ABI_DEALLOCATE(fock%calc_phase)
 endif
 if (allocated(fock%timerev)) then
    ABI_DEALLOCATE(fock%timerev)
 endif
 if (allocated(fock%tab_ibg)) then
    ABI_DEALLOCATE(fock%tab_ibg)
 endif
 if (allocated(fock%tab_icg)) then
    ABI_DEALLOCATE(fock%tab_icg)
 endif
 if (allocated(fock%tab_ikpt)) then
    ABI_DEALLOCATE(fock%tab_ikpt)
 endif

!* [description of IBZ and BZ]
!* Deallocate real arrays 
 if (allocated(fock%wtk_bz)) then
   ABI_DEALLOCATE(fock%wtk_bz)
 endif
 if (allocated(fock%kptns_bz)) then
    ABI_DEALLOCATE(fock%kptns_bz)
 endif
 if (allocated(fock%phase)) then
    ABI_DEALLOCATE(fock%phase)
 endif
!* Put the integer to 0
   fock%nkpt_bz=0

!* Deallocate real arrays 
 if (allocated(fock%eigen_ikpt)) then
    ABI_DEALLOCATE(fock%eigen_ikpt)
 endif

 ! Put the integer to 0
 fock%ieigen=0
 fock%ikpt=0
 fock%isppol=0

!* [intermediate variables for the calculation]
!* Deallocate real arrays 
!   if (allocated(fock%cwavef_r)) then
!      ABI_DEALLOCATE(fock%cwavef_r)
!   endif
!   if (allocated(fock%vlocpsi_r)) then
!      ABI_DEALLOCATE(fock%vlocpsi_r)
!   endif
!   if (allocated(fock%rhog_munu)) then
!      ABI_DEALLOCATE(fock%rhog_munu)
!   endif
!   if (allocated(fock%dummytab3)) then
!      ABI_DEALLOCATE(fock%dummytab3)
!   endif
!   if (allocated(fock%dummytab2)) then
!      ABI_DEALLOCATE(fock%dummytab2)
!   endif
!   if (allocated(fock%work_tmp3)) then
!      ABI_DEALLOCATE(fock%work_tmp3)
!   endif
!* Deallocate integer arrays
 if (allocated(fock%gbound_bz)) then
    ABI_DEALLOCATE(fock%gbound_bz)
 endif

!* [description of divergence in |q+G|=0]
!* Put the real (dp) to 0
 fock%gsqcut=0
 fock%divgq0=0
 fock%alpha=0

!* [description of size of arrays/pointers]
!* Put the integer to 0
 fock%mkpt=0
 fock%mkptband=0

 ABI_DATATYPE_DEALLOCATE(fock)

 DBG_EXIT("COLL")

end subroutine fock_destroy
!!***

!!****f* ABINIT/m_fock/fock_calc_ene
!! NAME
!!  fock_calc_ene
!!
!! FUNCTION
!!  Calculate the Fock contribution to the total energy
!!
!! INPUTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!  ikpt= reduced planewave coordinates.
!!
!! OUTPUT
!!  none
!!
!! SIDE EFFECTS
!!
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_exactX = Fock contribution to the total energy (Hartree)
!! 
!! NOTES
!!
!! If the cgocc_bz are not updated at each iteration, be careful to calculate Fock energy at the same frequency.
!! TO CHECK == CHANGE IN SOME DEFINTIONS 
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!
!! SOURCE

subroutine fock_calc_ene(dtset,fock,fock_energy,ikpt,nband,occ)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_calc_ene'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,nband
 real(dp),intent(inout) :: fock_energy
 type(dataset_type),intent(in) :: dtset
 type(fock_type),pointer :: fock
!arrays
 real(dp),intent(in) :: occ(nband)

!Local variables-------------------------------
 integer :: iband
 
! *************************************************************************
 
 do iband=1,nband
   ! Select only the occupied states (such that fock%occ_bz > 10^-8)
   if (abs(occ(iband))>tol8) then
     fock_energy=fock_energy + half*fock%eigen_ikpt(iband)*occ(iband)*dtset%wtk(ikpt)
     !* Sum the contribution of each occupied states at point k_i
     !* No need to multiply %wtk by ucvol since there is no factor 1/ucvol in the definition of %wtk
   end if
 end do

end subroutine fock_calc_ene
!!***

!!****f* ABINIT/m_fock/fock_update_exc
!! NAME
!!  fock_update_exc
!!
!! FUNCTION
!!  Update the value of energies%e_xc and energies%e_xcdc with Fock contribution
!!
!! INPUTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!
!! OUTPUT
!!  none
!!
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_fock= Fock contribution to the total energy (Hartree)
!! 
!! NOTES
!!   If the cgocc_bz are not updated at each iteration, be careful to calculate Fock energy at the same frequency.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine fock_update_exc(fock,fock_energy,xc_energy,xcdc_energy)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_update_exc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(fock_type),pointer,intent(in) :: fock
 real(dp),intent(in) :: fock_energy
 real(dp),intent(inout) :: xc_energy,xcdc_energy

! *************************************************************************
 
  xc_energy = xc_energy + fock%alpha*fock_energy
  xcdc_energy = xcdc_energy + two*fock%alpha*fock_energy
! CMartins : For an atom, ewald should be set to zero (at the beginning of the loop) and 
! the contribution in !|q+G|=0 should be an approximation to the missing component of Vloc in G=0
! energies%e_ewald=energies%e_ewald-half*fock%divgq0*fock%wtk_bz(1)*piinv

end subroutine fock_update_exc
!!***

!----------------------------------------------------------------------

!!****f* m_fock/fock_set_getghc_call
!! NAME
!!  fock_set_getghc_call
!!
!! FUNCTION
!!  Set the value of fock%getghc_call, Returns the old value
!!
!! PARENTS
!!
!! SOURCE

integer function fock_set_getghc_call(fock, new) result(old)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_set_getghc_call'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: new
 type(fock_type),intent(inout) :: fock

! *************************************************************************

 old = fock%getghc_call_
 fock%getghc_call_ = new

end function fock_set_getghc_call
!!***

!----------------------------------------------------------------------

!!****f* m_fock/fock_get_getghc_call
!! NAME
!!  fock_get_getghc_call
!!
!! FUNCTION
!!  Returns the value of fock%getghc_call_
!!
!! PARENTS
!!
!! SOURCE

pure integer function fock_get_getghc_call(fock)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_get_getghc_call'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(fock_type),intent(in) :: fock

! *************************************************************************

 fock_get_getghc_call = fock%getghc_call_

end function fock_get_getghc_call
!!***

!----------------------------------------------------------------------

!!****f* m_fock/fock_print
!! NAME
!!  fock_print 
!!
!! FUNCTION
!!  Print info on the fock_type data type
!!
!! INPUTS
!!  fock<crystal_t>=The object
!!  [unit]=Unit number for output
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS"
!!  [header]=String to be printed as header for additional info.
!!
!! OUTPUT
!!  Only printing 
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine fock_print(fock,header,unit,mode_paral,prtvol) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
 type(fock_type),intent(in) :: fock

!Local variables-------------------------------
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg      

! ********************************************************************* 

 my_unt=std_out; if (PRESENT(unit)) my_unt=unit
 my_prtvol=0 ; if (PRESENT(prtvol)) my_prtvol=prtvol 
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

 msg=' ==== Info on fock_type ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 ! Important dimensions
 call wrtout(my_unt,sjoin(" my_nsppol ...",itoa(fock%my_nsppol)),my_mode)
 call wrtout(my_unt,sjoin(" nkpt_bz .....",itoa(fock%nkpt_bz)),my_mode)

 ! Options
 call wrtout(my_unt,sjoin(" cg_typ ......",itoa(fock%cg_typ)),my_mode)
 call wrtout(my_unt,sjoin(" nnsclo_hf ...",itoa(fock%nnsclo_hf)),my_mode)
 call wrtout(my_unt,sjoin(" ixc .........",itoa(fock%ixc)),my_mode)
 call wrtout(my_unt,sjoin(" alpha .......",ftoa(fock%alpha)),my_mode)

 write(msg,"(a,f12.1,a)")" Memory required for HF u(r) states: ",product(shape(fock%cwaveocc_bz)) * dp * b2Mb, " [Mb]"
 call wrtout(my_unt,msg,my_mode)

 ! Extra info.
 if (my_prtvol > 0) then
   call wrtout(my_unt,"Extra info not available",my_mode)
 end if

end subroutine fock_print
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/bare_vqg
!! NAME
!! bare_vqg
!!
!! FUNCTION
!! Compute bare coulomb term in G-space on the FFT mesh i.e. 4pi/(G+q)**2
!!
!! NOTES
!!
!! INPUTS
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  gsqcut=cutoff value on G**2 for sphere inside fft box. (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  divgq0= value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq. Used if q = Gamma
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  izero=if 1, unbalanced components of V(q,g) are set to zero
!!  nfft=Total number of FFT grid points.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  vqg(nfft)=4pi/(G+q)**2, G=0 component is set to divgq0/pi if q = Gamma. 
!!
!! NOTES
!!  This routine operates on the full FFT mesh. DO NOT PASS MPI_TYPE
!!  One can easily implemente MPI-FFT by just calling this routine and then
!!  extracting the G-vectors treated by the node.
!!
!! PARENTS
!!
!! SOURCE

subroutine bare_vqg(qphon,gsqcut,divgq0,gmet,izero,nfft,ngfft,vqg) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bare_vqg'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,nfft
 real(dp),intent(in) :: gsqcut,divgq0
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),qphon(3)
 real(dp),intent(out) ::  vqg(nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex1=1
 integer :: i1,i2,i23,i3,id1,id2,id3
 integer :: ig,ig1min,ig1,ig1max,ig2,ig2min,ig2max,ig3,ig3min,ig3max
 integer :: ii,ii1,ing,n1,n2,n3,qeq0,qeq05
 real(dp),parameter :: tolfix=1.000000001e0_dp ! Same value as the one used in hartre
 real(dp) :: cutoff,den,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3
!arrays
 integer :: id(3)
 real(dp),allocatable :: gq(:,:)

! *************************************************************************

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)

 ! Initialize a few quantities
 cutoff=gsqcut*tolfix

 qeq0=0; if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1
 qeq05=0
 if (qeq0==0) then
   if (abs(abs(qphon(1))-half)<tol12.or.abs(abs(qphon(2))-half)<tol12.or. &
&   abs(abs(qphon(3))-half)<tol12) qeq05=1
 end if

 ! In order to speed the routine, precompute the components of g+q
 ! Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qphon(ii)
   end do
 end do
 ig1max=-1;ig2max=-1;ig3max=-1
 ig1min=n1;ig2min=n2;ig3min=n3

 id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

 ! Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
     gqgm12=gq(2,i2)*gmet(1,2)*2
     gqg2p3=gqgm13+gqgm12

     i23=n1*(i2-1 +(n2)*(i3-1))
     ! Do the test that eliminates the Gamma point outside of the inner loop
     ii1=1
     if (i23==0 .and. qeq0==1  .and. ig2==0 .and. ig3==0) then
       ii1=2
       ! value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq
       vqg(1+i23)=divgq0*piinv
     end if

     ! Final inner loop on the first dimension (note the lower limit)
     do i1=ii1,n1
       gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
       ii=i1+i23

       if(gs<=cutoff)then
         ! Identify min/max indexes (to cancel unbalanced contributions later)
         ! Count (q+g)-vectors with similar norm
         if ((qeq05==1).and.(izero==1)) then
           ig1=i1-(i1/id1)*n1-1
           ig1max=max(ig1max,ig1); ig1min=min(ig1min,ig1)
           ig2max=max(ig2max,ig2); ig2min=min(ig2min,ig2)
           ig3max=max(ig3max,ig3); ig3min=min(ig3min,ig3)
         end if

         den=piinv/gs
         vqg(ii)=den
       else 
         ! gs>cutoff
         vqg(ii)=zero
       end if
     end do ! End loop on i1
   end do ! End loop on i2
 end do ! End loop on i3

 if (izero==1) then
   ! Set contribution of unbalanced components to zero

   if (qeq0==1) then !q=0
     call zerosym(vqg,cplex1,n1,n2,n3)

   else if (qeq05==1) then 
     !q=1/2; this doesn't work in parallel
     ig1=-1;if (mod(n1,2)==0) ig1=1+n1/2
     ig2=-1;if (mod(n2,2)==0) ig2=1+n2/2
     ig3=-1;if (mod(n3,2)==0) ig3=1+n3/2
     if (abs(abs(qphon(1))-half)<tol12) then
       if (abs(ig1min)<abs(ig1max)) ig1=abs(ig1max)
       if (abs(ig1min)>abs(ig1max)) ig1=n1-abs(ig1min)
     end if
     if (abs(abs(qphon(2))-half)<tol12) then
       if (abs(ig2min)<abs(ig2max)) ig2=abs(ig2max)
       if (abs(ig2min)>abs(ig2max)) ig2=n2-abs(ig2min)
     end if
     if (abs(abs(qphon(3))-half)<tol12) then
       if (abs(ig3min)<abs(ig3max)) ig3=abs(ig3max)
       if (abs(ig3min)>abs(ig3max)) ig3=n3-abs(ig3min)
     end if
     call zerosym(vqg,cplex1,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3)
   end if
 end if

 ABI_DEALLOCATE(gq)

end subroutine bare_vqg
!!***

end module m_fock
!!***
