!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hamiltonian
!! NAME
!! m_hamiltonian
!!
!! FUNCTION
!!  This module provides the definition of the gs_hamiltonian_type and rf_hamiltonian_type
!!  datastructures used in the "getghc" and "getgh1c" routines to apply the Hamiltonian (or
!!  its derivative) on a wavefunction.
!!  Methods to initialize or destroy the objects are defined here.
!!  It also defines the ddiago_ctl_type structures datatype used to control the algorithm
!!  used in ks_ddiago for performing the direct diagonalization of the KS Hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (MG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

module m_hamiltonian

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_fstrings,          only : toupper
 use m_pawtab,            only : pawtab_type
 use m_fftcore,           only : kpgsph
 use m_paw_ij,            only : paw_ij_type
 use m_paral_atom,        only : get_my_atmtab, free_my_atmtab
 use m_electronpositron,  only : electronpositron_type, electronpositron_calctype
 use m_mpinfo,            only : destroy_mpi_enreg

 implicit none

 private
 
 public ::  pawdij2ekb
 public ::  pawdij2e1kb
!!***

!----------------------------------------------------------------------

!!****t* m_hamiltonian/gs_hamiltonian_type
!! NAME
!! gs_hamiltonian_type
!!
!! FUNCTION
!! This datastructure contains the information about one Hamiltonian,
!! needed in the "getghc" routine, that apply the Hamiltonian
!! on a wavefunction.
!! About the non-local part of the Hamiltonian
!! The operator Onl has the following general form:
!! $Onl=sum_{R,lmn,l''m''n''} {|P_{Rlmn}> Enl^{R}_{lmn,l''m''n''} <P_{Rl''m''n''}|}$
!! Operator Onl is -- in the typical case -- the nonlocal potential.
!! - In a classical plane-wave calculation, $Enl^{R}_{lmn,l''m''n''}$ is the
!!   Kleinmann-Bylander energy $Ekb^{R}_{ln}$.
!! - In a PAW calculation, $Enl^{R}_{lmn,l''m''n''}$ can either be the nonlocal
!!   contribution to the Hamiltonian $D_{ij}$ or the overlap matrix $S_{ij}$.
!! - The |P_{Rlmn}> are the projector functions.
!!
!! SOURCE

 type,public :: gs_hamiltonian_type

! Integer scalars

  integer :: dimekb1
   ! First dimension of Ekb (see ekb in defs_datatypes.F90)
   ! Same as psps%dimekb
   ! ->Norm conserving : Max. number of Kleinman-Bylander energies
   !                     for each atom type
   !                     dimekb1=lnmax
   ! ->PAW : Max. number of Dij coefficients connecting projectors
   !                     for each atom
   !                     dimekb1=lmnmax*(lmnmax+1)/2

  integer :: dimekb2
   ! Second dimension of Ekb (see ekb in defs_datatypes.F90)
   ! ->Norm conserving psps: dimekb2=ntypat
   ! ->PAW                 : dimekb2=natom

  integer :: istwf_k
   ! option parameter that describes the storage of wfs (k-dependent)

  integer :: lmnmax
   ! Maximum number of different l,m,n components over all types of psps.
   ! same as dtset%lmnmax

  integer :: matblk
   ! dimension of the array ph3d

  integer :: mgfft
   ! maximum size for 1D FFTs (same as dtset%mgfft)

  integer :: mproj  ! TO BE SUPPRESSED LATER
   ! Maximum number of non-local projectors over all angular momenta
   !  and type of psps
   ! 0 only if all psps are local
   ! same as psps%mproj

  integer :: mpsang
   ! Highest angular momentum of non-local projectors over all type of psps.
   ! shifted by 1 : for all local psps, mpsang=0; for largest s, mpsang=1,
   ! for largest p, mpsang=2; for largest d, mpsang=3; for largest f, mpsang=4
   ! This gives also the number of non-local "channels"
   ! same as psps%mpsang

  integer :: mpssoang
   ! Maximum number of channels, including those for treating the spin-orbit coupling
   ! when mpspso=1, mpssoang=mpsang
   ! when mpspso=2, mpssoang=2*mpsang-1
   ! same as psps%mpssoang

  integer :: natom
   ! The number of atoms for this dataset same as dtset%natom

  integer :: nfft
   ! number of FFT grid points same as dtset%nfft

  integer :: npw
   ! number of plane waves (k-dependent)

  integer:: nspinor
   ! Number of spinorial components

  integer :: ntypat
   ! Number of types of pseudopotentials same as dtset%ntypat

  integer :: nvloc
   ! Number of components of vloc
   ! usually, nvloc=1, except in the non-collinear magnetism case, where nvloc=4

  integer :: n4,n5,n6
   ! same as ngfft(4:6)

  integer :: use_gpu_cuda
  ! governs wheter we do the hamiltonian calculation on gpu (1) or not

  integer :: usepaw
   ! if usepaw=0 , use norm-conserving psps part of the code
   ! is usepaw=1 , use paw part of the code

  integer :: useylm
   ! governs the way the nonlocal operator is to be applied:
   !   1=using Ylm, 0=using Legendre polynomials

! Integer arrays

  integer, allocatable :: atindx(:)
   ! atindx(natom)
   ! index table for atoms (see gstate.f)

  integer, allocatable :: atindx1(:)
   ! atindx1(natom)
   ! index table for atoms, inverse of atindx (see gstate.f)

  integer, allocatable :: gbound(:,:)
   ! gbound(2*mgfft+8,2)
   ! G sphere boundary

  integer, allocatable :: indlmn(:,:,:)
   ! indlmn(6,lmnmax,ntypat)
   ! For each type of psp,
   ! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
   !                                or i=lmn (if useylm=1)

! integer, allocatable :: indpw_k(:,:)
   ! indpw_k(4,npw)
   ! array which gives fft box index for given basis sphere
   ! This component was taken away : CPU time problem !

! integer, allocatable :: kg_k(:,:)
   ! kg_k(3,npw)
   ! G vector coordinates with respect to reciprocal lattice translations
   ! This component was taken away : CPU time problem !

  integer, allocatable :: nattyp(:) 
   ! nattyp(ntypat)
   ! # of atoms of each type

  integer :: ngfft(18)
   ! ngfft(1:3)=integer fft box dimensions
   ! ngfft(4:6)=integer fft box dimensions, might be augmented for CPU speed
   ! ngfft(7)=fftalg
   ! ngfft(8)=fftalg
   ! same as dtset%ngfft

  integer :: nloalg(5)
   ! governs the choice of the algorithm for non-local operator same as dtset%nloalg

  integer, allocatable :: pspso(:) 
   ! pspso(ntypat)
   ! For each type of psp, 1 if no spin-orbit component is taken
   ! into account, 2 if a spin-orbit component is used

  integer, allocatable :: typat(:)  
   ! typat(natom)
   ! type of each atom

! Real (real(dp)) scalar

  real(dp) :: ucvol
   ! unit cell volume (Bohr**3)

! Real (real(dp)) arrays

  real(dp), allocatable :: ekb(:,:,:)
   ! ekb(dimekb1,dimekb2,nspinor**2)
   !  ->Norm conserving : (Real) Kleinman-Bylander energies (hartree)
   !          for number of basis functions (l,n) (lnmax)
   !          and number of atom types (ntypat)
   !          dimekb1=lnmax ; dimekb2=ntypat
   !  ->PAW : (Real, symmetric) Frozen part of Dij coefficients
   !                            to connect projectors
   !          for number of basis functions (l,m,n) (lmnmax)
   !          and number of atom (natom)
   !          dimekb1=lmnmax*(lmnmax+1)/2 ; dimekb2=natom
   ! NOTE (MT) : ekb (norm-conserving) is now diagonal (one dimension
   !             lnmax); it would be easy to give it a second
   !             (symmetric) dimension by putting
   !             dimekb1=lnmax*(lnmax+1)/2
   !             in the place of dimekb1=lnmax.
   ! %ekb is spin dependent in the case of PAW calculations.

  real(dp), allocatable :: sij(:,:)  
   ! sij(dimekb1,ntypat*usepaw) = overlap matrix for paw calculation

! real(dp), allocatable :: ffnl(:,:,:,:)
   ! ffnl(npw,2,lmnmax,ntypat)
   ! nonlocal form factors
   ! This component was taken away : CPU time problem !

  real(dp) :: gmet(3,3)
   ! reciprocal space metric tensor in Bohr**-2

  real(dp) :: gprimd(3,3)
   ! dimensional reciprocal space primitive translations (Bohr^-1)

! real(dp), allocatable :: kinpw(:)
   ! kinpw(npw)
   ! (modified) kinetic energy for each plane wave (Hartree)
   ! This component was taken away : CPU time problem !

  real(dp) :: kpoint(3)
   ! dimensionless k point coordinates wrt reciprocal lattice vectors. (k-dependent).

  real(dp), allocatable :: phkxred(:,:)  
   ! phkxred(2,natom)
   ! phase factors exp(2 pi k.xred) (k-dependent)

  real(dp), allocatable :: ph1d(:,:) 
   ! ph1d(2,3*(2*mgfft+1)*natom)
   ! 1-dim phase arrays for structure factor (see getph.f).

! real(dp), allocatable :: ph3d(:,:,:)
   ! ph3d(2,npw,matblk)
   ! 3-dim structure factors, for each atom and plane wave
   ! This component was taken away : CPU time problem !

! real(dp), allocatable :: vlocal(:,:,:,:)
   ! vlocal(n4,n5,n6,nvloc)
   ! local potential in real space, on the augmented fft grid
   ! This component was taken away : CPU time problem !

  real(dp),allocatable :: xred(:,:) 
   ! xred(3,natom)
   ! reduced coordinates of atoms (dimensionless)

! real(dp),allocatable :: ylm(:,:)
   ! ylm(npw,mpsang*mpsang*useylm)
   ! Real spherical harmonics for each k+G
   ! This component was taken away : CPU time problem !
  
 end type gs_hamiltonian_type

 public :: init_hamiltonian         ! Initialize the GS Hamiltonian
 public :: destroy_hamiltonian      ! Free the memory in the GS Hamiltonian
 public :: load_paw_hamiltonian     ! Setup of the k-dependent part of the PAW Hamiltonian.
 public :: finalize_hamiltonian     ! Setup of the k-dependent part of the Hamiltonian.
!!***

!----------------------------------------------------------------------

!!****t* m_hamiltonian/rf_hamiltonian_type
!! NAME
!! rf_hamiltonian_type
!!
!! FUNCTION
!! This datastructure contains few data about one 1st-order Hamiltonian,
!! needed in the "getgh1c" routine, that apply the 1st-order Hamiltonian
!! on a wavefunction.
!!
!! SOURCE

 type,public :: rf_hamiltonian_type

! Integer scalars

  integer :: dime1kb
   ! First dimension of E1kb, derivative of Ekb (see ekb in defs_datatypes.F90)
   ! with respect to a perturbation

  integer :: dimekb1
   ! First dimension of Ekb (see ekb in defs_datatypes.F90)
   ! Same as psps%dimekb

  integer :: dimekb2
   ! Second dimension of Ekb and E1kb (see ekb in defs_datatypes.F90)
   ! ->Norm conserving psps: dimekb2=ntypat
   ! ->PAW                 : dimekb2=natom

!  integer :: has_e1kbfr
   ! If 1, e1kbfr (in this datastructure) is allocated

!  integer :: has_e1kbsc
   ! If 1, e1kbsc (in this datastructure) is allocated

  integer :: lmnmax
   ! Maximum number of different l,m,n components over all types of psps.
   ! same as dtset%lmnmax

  integer:: nspinor
   ! Number of spinorial components

  integer :: usepaw
   ! if usepaw=0 , use norm-conserving psps part of the code
   ! is usepaw=1 , use paw part of the code

! Integer arrays

  integer,allocatable :: indlmn_typ(:,:,:)
   ! Atomic displacement perturbation only:
   ! Array giving l,m,n,lm,ln,spin for i=lmn (for the displaced atom)
   ! indlmn_typ(6,lmnmax,1)

  integer,allocatable :: pspso_typ(:)
   ! Atomic displacement perturbation only:
   ! Spin-orbit characteristics for the displaced atom (1, 2, or 3)
   ! pspso_typ(1)

! Real arrays

  real(dp),allocatable :: e1kbfr(:,:,:)
   ! Frozen part of 1st derivative of ekb
   ! for the considered pertubation (not depending on VHxc^(1))
   ! e1kbfr(dime1kb,dimekb2,nspinor**2)

  real(dp),allocatable :: e1kbsc(:,:,:)
   ! Self-consistent part of 1st derivative of ekb
   ! for the considered pertubation (depending on VHxc^(1))
   ! e1kbsc(dime1kb,dimekb2,nspinor**2)

  real(dp),allocatable :: ekb_typ(:,:,:)
   ! Atomic displacement perturbation only:
   ! Kleinman-Bylander energies (norm-conserving psp) or Dij coefficients (PAW)
   ! for the displaced atom
   ! ekb_typ(dimekb1,1,nspinor**2)

  real(dp),allocatable :: sij_typ(:,:)
   ! Atomic displacement perturbation + PAW only:
   ! Overlap matrix components
   ! sij_typ(dimekb1,1)

 end type rf_hamiltonian_type

 public ::  init_rf_hamiltonian      ! Initialize the RF Hamiltonian
 public ::  destroy_rf_hamiltonian   ! Free the memory in the RF Hamiltonian
 public ::  load_paw_rf_hamiltonian  ! Setup of the k-dependent part of the PAW Hamiltonian.
!!***

!----------------------------------------------------------------------

!!****t* m_hamiltonian/ddiago_ctl_type
!! NAME
!!  ddiago_ctl_type
!!
!! FUNCTION
!!  Structure storing the variables controlling the direct diagonalization of the Kohn-Sham Hamiltonian.
!!  Mainly used for debugging (and in the KSS code!)
!!
!! SOURCE

 type, public :: ddiago_ctl_type

  integer :: isppol
   ! The spin component of the Hamiltonian (1 if nspinor==1 or nsppol==1).

  integer :: istwf_k
   ! Option defining whether time-reversal symmetry is used at particular k-points
   ! If 0, the code will automatically use TR symmetry if possible (depending on the k-point)

  integer :: nband_k
   ! Number of bands to be calculated.

  integer :: npw_k
  ! The number of planes waves for the wavefunctions taking into account time-reversal symmetry.

  integer :: npwtot
  ! The number of planes waves in the Hamiltonian without taking into account istwf_k

  integer :: nspinor
  ! Number of spinorial components.

  integer :: prtvol
   ! Flag controlling the verbosity.

  integer :: use_scalapack
  ! 0 if diagonalization is done in sequential on each node.
  ! 1 to use scalapack

  real(dp) :: abstol
   ! used fro RANGE="V","I", and "A" when do_full_diago=.FALSE.
   ! The absolute error tolerance for the eigenvalues. An approximate eigenvalue is accepted
   ! as converged when it is determined to lie in an interval [a,b] of width less than or equal to
   !
   !         ABSTOL + EPS *   max( |a|,|b| ) ,
   !
   ! where EPS is the machine precision.  If ABSTOL is less than or equal to zero, then  EPS*|T|  will be used in its place,
   ! where |T| is the 1-norm of the tridiagonal matrix obtained by reducing A to tridiagonal form.
   !
   ! Eigenvalues will be computed most accurately when ABSTOL is
   ! set to twice the underflow threshold 2*DLAMCH('S'), not zero.
   ! If this routine returns with INFO>0, indicating that some
   ! eigenvectors did not converge, try setting ABSTOL to 2*DLAMCH('S').

  real(dp) :: ecut
   ! The cutoff energy for the plane wave basis set.

  real(dp) :: ecutsm
   ! Smearing energy for plane wave kinetic energy (Ha)

  real(dp) :: effmass
   ! Effective mass for electrons (usually one).

  logical :: do_full_diago
  ! Specifies whether direct or partial diagonalization will be performed.
  ! Meaningful only if RANGE='A'.

  integer :: ilu(2)
   ! If RANGE='I', the indices (in ascending order) of the smallest and largest eigenvalues to be returned.
   ! il=ilu(1), iu=ilu(2) where
   ! 1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. NOT used if RANGE = 'A' or 'V'.

  integer :: nloalg(5)

  real(dp) :: kpoint(3)
   ! The k-point in reduced coordinates at which the Hamiltonian is diagonalized.

  real(dp) :: vlu(2)
   ! If RANGE='V', the lower and upper bounds of the interval to
   ! be searched for eigenvalues. vl=vlu(1) and vu=vlu(2) with VL < VU.
   ! Not referenced if RANGE = 'A' or 'I'.

  character(len=1) :: jobz
   ! character defining whether wavefunctions are required (lapack option).
   ! "N":  Compute eigenvalues only;
   ! "V":  Compute eigenvalues and eigenvectors.

  character(len=1) :: range
   ! character defining the subset of eigenstates that will be calculated (lapack option).
   ! "A": all eigenvalues will be found.
   ! "V": all eigenvalues in the half-open interval (VL,VU] will be found.
   ! "I": the IL-th through IU-th eigenvalues will be found.

  !$character(len=fnlen) :: fname
  ! The name of the file storing the eigenvectors and eigenvalues (only if jobz="V")

 end type ddiago_ctl_type

 public ::  init_ddiago_ctl
!!***

CONTAINS  !===========================================================

!----------------------------------------------------------------------

!!****f* m_hamiltonian/destroy_hamiltonian
!! NAME
!!  destroy_hamiltonian
!!
!! FUNCTION
!!  Clean and destroy gs_hamiltonian_type datastructure
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=All dynamic memory defined in the structure is deallocated.
!!
!! PARENTS
!!      energy,ks_ddiago,m_shirley,nstdy3,nstpaw3,rhofermi3,vtorho,vtorho3
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine destroy_hamiltonian(Ham)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_hamiltonian'
!End of the abilint section

 type(gs_hamiltonian_type),intent(inout) :: Ham

! *************************************************************************

 !@gs_hamiltonian_type

! Integer arrays.
 if (allocated(Ham%atindx))  then
   ABI_DEALLOCATE(Ham%atindx)
 end if
 if (allocated(Ham%atindx1))  then
   ABI_DEALLOCATE(Ham%atindx1)
 end if
 if (allocated(Ham%gbound))  then
   ABI_DEALLOCATE(Ham%gbound)
 end if
 if (allocated(Ham%indlmn ))  then
   ABI_DEALLOCATE(Ham%indlmn)
 end if
 !if (allocated(indpw_k)) deallocate(indpw_k)
 !if (allocated(kg_k)) deallocate(kg_k)
 if (allocated(Ham%nattyp))  then
   ABI_DEALLOCATE(Ham%nattyp)
 end if
 if (allocated(Ham%pspso))  then
   ABI_DEALLOCATE(Ham%pspso)
 end if
 if (allocated(Ham%typat))  then
   ABI_DEALLOCATE(Ham%typat)
 end if

! Real arrays.
 if (allocated(Ham%ekb))   then
   ABI_DEALLOCATE(Ham%ekb)
 end if
 if (allocated(Ham%sij))   then
   ABI_DEALLOCATE(Ham%sij)
 end if
 !if (allocated(Ham%ffnl))  deallocated(Ham%ffnl)
 !if (allocated(Ham%kinpw))  deallocated(Ham%kinpw)
 if (allocated(Ham%phkxred))   then
   ABI_DEALLOCATE(Ham%phkxred)
 end if
 if (allocated(Ham%ph1d))   then
   ABI_DEALLOCATE(Ham%ph1d)
 end if
 !if (allocated(Ham%ph3d))  deallocated(Ham%ph3d)
 !if (allocated(Ham%vlocal))  deallocated(Ham%vlocal)
 if (allocated(Ham%xred))   then
   ABI_DEALLOCATE(Ham%xred)
 end if
 !if (allocated(Ham%ylm))  deallocated(Ham%ylm)

#if defined HAVE_GPU_CUDA
 if(Ham%use_gpu_cuda==1) then
   call gpu_finalize_ham_data()
 end if
#endif

end subroutine destroy_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/init_hamiltonian
!! NAME
!!  init_hamiltonian
!!
!! FUNCTION
!!  Creation method for the gs_hamiltonian_type structure.
!!  It allocates memory and initializes all quantities that do not depend on the k-point or spin.
!!
!! INPUTS
!!  natom=Number of atoms in the unit cell.
!!  nfft=Number of FFT grid points (for this processors).
!!  ntypat=Number of type of atoms.
!!  nspinor=Number of spinorial components
!!  nspden=Number of spin density components.
!!  mgfft=Maximum size for 1D FFTs i.e., MAXVAL(ngfft(1:3))
!!  psps<pseudopotential_type>=structure datatype gathering data on the pseudopotentials.
!!  [electronpositron<electronpositron_type>]=Structured datatype storing data for the
!!    electron-positron two-component DFT (optional).
!!  ngfft(18)=integer array with FFT box dimensions and other information on FFTs, for the FINE rectangular grid.
!!  nloalg(5)=governs the choice of the algorithm for non-local operator
!!  typat(natom)=Type of each atom.
!!  [ph1d(2,3*(2*mgfft+1)*natom)]= 1-dimensions phase arrays for structure factor (see getph.f).
!!    Recalculated inside the routine if not present in input.
!!  rprimd(3,3)=Direct lattice vectors in Bohr.
!!  xred(3,natom)=Reduced coordinates of the atoms.
!!  pawtab(ntypat*psps%usepaw)<pawtab_type>=PAW TABulated data initialized at start.
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=Structured datatype almost completely initialized:
!!   * Basic variables and dimensions are transfered to the structure.
!!   * All pointers are allocated with correct dimensions.
!!   * Quantities that do not depend on the k-point or spin are initialized.
!!
!! PARENTS
!!      energy,ks_ddiago,m_shirley,nstdy3,nstpaw3,rhofermi3,vtorho,vtorho3
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine init_hamiltonian(gs_hamk,Psps,pawtab,nspinor,nspden,natom,ntypat,typat,&
&                           xred,nfft,mgfft,ngfft,rprimd,nloalg,ph1d,electronpositron,use_gpu_cuda)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_hamiltonian'
 use interfaces_41_geometry
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,natom,ntypat,nspinor,nspden,mgfft
 integer,optional,intent(in) :: use_gpu_cuda
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(pseudopotential_type),intent(in) :: psps
 type(electronpositron_type),optional,pointer :: electronpositron
!arrays
 integer,intent(in) :: ngfft(18),nloalg(5)
 integer,intent(in) :: typat(natom)
 real(dp),optional,intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 type(pawtab_type),intent(in)  :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: itypat,iat,indx,ilmn,cplex,cplex_dij
 real(dp) :: ucvol
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)

! *************************************************************************

 !@gs_hamiltonian_type
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ABI_CHECK(mgfft==MAXVAL(ngfft(1:3)),"Wrong mgfft")

!Allocate the arrays of the Hamiltonian whose dimensions do not depend on k
 ABI_ALLOCATE(gs_hamk%atindx,(natom))
 ABI_ALLOCATE(gs_hamk%atindx1,(natom))
 ABI_ALLOCATE(gs_hamk%typat,(natom))
 gs_hamk%typat=typat(1:natom)
 ABI_ALLOCATE(gs_hamk%gbound,(2*mgfft+8,2))
 gs_hamk%gbound(:,:)=0
 ABI_ALLOCATE(gs_hamk%indlmn,(6,psps%lmnmax,ntypat))
 ABI_ALLOCATE(gs_hamk%nattyp,(ntypat))
 ABI_ALLOCATE(gs_hamk%phkxred,(2,natom))
 ABI_ALLOCATE(gs_hamk%ph1d,(2,3*(2*mgfft+1)*natom))
 ABI_ALLOCATE(gs_hamk%pspso,(ntypat))
 ABI_ALLOCATE(gs_hamk%xred,(3,natom))

!Initialize most of the Hamiltonian
 indx=1
 do itypat=1,ntypat
   gs_hamk%nattyp(itypat)=0
   do iat=1,natom
     if (typat(iat)==itypat) then
       gs_hamk%atindx (iat )=indx
       gs_hamk%atindx1(indx)=iat
       indx=indx+1
       gs_hamk%nattyp(itypat)=gs_hamk%nattyp(itypat)+1
     end if
   end do
 end do

 gs_hamk%gmet(:,:)  =gmet(:,:)
 gs_hamk%gprimd(:,:)=gprimd(:,:)
 gs_hamk%indlmn(:,:,:)=psps%indlmn(:,:,:)
 gs_hamk%lmnmax     =psps%lmnmax
 gs_hamk%mgfft      =mgfft
 gs_hamk%mproj      =psps%mproj
 gs_hamk%mpsang     =psps%mpsang
 gs_hamk%mpssoang   =psps%mpssoang
 gs_hamk%natom      =natom
 gs_hamk%nfft       =nfft
 gs_hamk%ngfft(:)   =ngfft(:)
 gs_hamk%nloalg(:)  =nloalg(:)
 gs_hamk%matblk=nloalg(4); if (nloalg(1)>0) gs_hamk%matblk=natom
 gs_hamk%nspinor    =nspinor
 gs_hamk%ntypat     =ntypat

 gs_hamk%nvloc=1; if(nspden==4)gs_hamk%nvloc=4
 gs_hamk%n4         =ngfft(4)
 gs_hamk%n5         =ngfft(5)
 gs_hamk%n6         =ngfft(6)
 gs_hamk%usepaw     =psps%usepaw
 if (PRESENT(ph1d)) then
   gs_hamk%ph1d(:,:)  =ph1d(:,:)
 else ! Recalculate structure factor phases
   call getph(gs_hamk%atindx,natom,ngfft(1),ngfft(2),ngfft(3),gs_hamk%ph1d,xred)
 end if
 gs_hamk%pspso(:)   =psps%pspso(1:ntypat)
 gs_hamk%ucvol      =ucvol
 gs_hamk%useylm     =psps%useylm
 gs_hamk%use_gpu_cuda=0
 if(PRESENT(use_gpu_cuda)) gs_hamk%use_gpu_cuda=use_gpu_cuda
 gs_hamk%xred(:,:)  =xred(:,:)

! ===========================
! ==== Non-local factors ====
! ===========================

 if (gs_hamk%usepaw==0) then ! Norm-conserving: use constant Kleimann-Bylander energies.
   gs_hamk%dimekb1=psps%dimekb
   gs_hamk%dimekb2=ntypat
   ABI_ALLOCATE(gs_hamk%ekb,(psps%dimekb,ntypat,nspinor**2))
   ABI_ALLOCATE(gs_hamk%sij,(0,0))
   gs_hamk%ekb(:,:,1)=psps%ekb(:,:)
   if (nspinor==2) then
     gs_hamk%ekb(:,:,2)=psps%ekb(:,:)
     gs_hamk%ekb(:,:,3:4)=zero
   end if
   if (PRESENT(electronpositron)) then
     if (electronpositron_calctype(electronpositron)==1) gs_hamk%ekb(:,:,:)=-gs_hamk%ekb(:,:,:)
   end if

!  Update enl on GPU (will do it later for PAW)
#if defined HAVE_GPU_CUDA
   if (gs_hamk%use_gpu_cuda==1) then
     call gpu_update_ham_data(gs_hamk%ekb,gs_hamk%sij,gs_hamk%gprimd,gs_hamk%dimekb1,&
&                    gs_hamk%dimekb2,gs_hamk%nspinor,gs_hamk%ntypat,gs_hamk%usepaw)
   end if
#endif

 else ! PAW: store overlap coefficients and allocate memory for Dij coefficients (spin dependent)
   cplex=1;cplex_dij=max(cplex,nspinor)
   gs_hamk%dimekb1=psps%dimekb*cplex_dij
   gs_hamk%dimekb2=natom
   ABI_ALLOCATE(gs_hamk%ekb,(gs_hamk%dimekb1,gs_hamk%dimekb2,nspinor**2))
   gs_hamk%ekb(:,:,:)=zero
   ABI_ALLOCATE(gs_hamk%sij,(gs_hamk%dimekb1,ntypat))
   do itypat=1,ntypat
     if (cplex_dij==1) then
       gs_hamk%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
     else
       do ilmn=1,pawtab(itypat)%lmn2_size
         gs_hamk%sij(2*ilmn-1,itypat)=pawtab(itypat)%sij(ilmn)
         gs_hamk%sij(2*ilmn  ,itypat)=zero
       end do
     end if
     if (cplex_dij*pawtab(itypat)%lmn2_size<gs_hamk%dimekb1) then
       gs_hamk%sij(cplex_dij*pawtab(itypat)%lmn2_size+1:gs_hamk%dimekb1,itypat)=zero
     end if
   end do
 end if

end subroutine init_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/finalize_hamiltonian
!! NAME
!!  finalize_hamiltonian
!!
!! FUNCTION
!!  Setup of the k-dependent part of the Hamiltonian.
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=
!!
!! PARENTS
!!      energy,m_shirley,nstdy3,nstpaw3,rhofermi3,vtorho,vtorho3
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine finalize_hamiltonian(gs_hamk,npw_k,istwfk,kpoint)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'finalize_hamiltonian'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,istwfk
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
!arrays
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer :: iat,iatom
 real(dp) :: arg

! *************************************************************************

 !@gs_hamiltonian_type

!Setup of the k-dependent part of the Hamiltonian.
 gs_hamk%npw      = npw_k
 gs_hamk%istwf_k  = istwfk
 gs_hamk%kpoint(:)= kpoint(:)

! Allocate the arrays phkxred and ph3d, compute phkxred and eventually ph3d.
 do iat=1,gs_hamk%natom
   iatom=gs_hamk%atindx(iat)
   arg=two_pi*DOT_PRODUCT(kpoint,gs_hamk%xred(:,iat))
   gs_hamk%phkxred(1,iatom)=DCOS(arg)
   gs_hamk%phkxred(2,iatom)=DSIN(arg)
 end do

end subroutine finalize_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/load_paw_hamiltonian
!! NAME
!!  load_paw_hamiltonian
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!
!! FUNCTION
!!  Setup of the k-dependent part of the Hamiltonian.
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=
!!
!! PARENTS
!!      energy,m_shirley,nstdy3,nstpaw3,rhofermi3,vtorho,vtorho3
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine load_paw_hamiltonian(gs_hamk,isppol,paw_ij, &
&          mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'load_paw_hamiltonian'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol
 integer,optional,intent(in) :: mpi_comm_atom
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
!arrays
 integer,optional,target,intent(in)  :: mpi_atmtab(:)
 type(paw_ij_type),intent(in) :: paw_ij(:)

!Local variables-------------------------------
!scalars
 integer :: my_comm_atom,my_natom
 logical :: my_atmtab_allocated,paral_atom
!arrays
 integer,pointer :: my_atmtab(:)

! *************************************************************************

 !@gs_hamiltonian_type

 if (gs_hamk%usepaw==0) return

!Set up parallelism over atoms
 my_natom=size(paw_ij);my_comm_atom=xmpi_self
 paral_atom=(present(mpi_comm_atom).and.my_natom/=gs_hamk%natom)
 my_comm_atom=xmpi_self;if (paral_atom) my_comm_atom=mpi_comm_atom
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,gs_hamk%natom,my_natom_ref=my_natom)

!Retrieve PAW Dij coefficients for this spin component
 gs_hamk%ekb=zero ! Need otherwise valgrind complains ... (even if I knwow it is an output argument...)
 call pawdij2ekb(gs_hamk%ekb,paw_ij,gs_hamk%dimekb1,isppol,my_atmtab,my_comm_atom,gs_hamk%nspinor)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

!Update enl and sij on GPU
#if defined HAVE_GPU_CUDA
 if (gs_hamk%use_gpu_cuda==1) then
   call gpu_update_ham_data(gs_hamk%ekb,gs_hamk%sij,gs_hamk%gprimd,gs_hamk%dimekb1,&
&                  gs_hamk%dimekb2,gs_hamk%nspinor,gs_hamk%ntypat,gs_hamk%usepaw)
 end if
#endif

end subroutine load_paw_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/destroy_rf_hamiltonian
!! NAME
!!  destroy_rf_hamiltonian
!!
!! FUNCTION
!!  Clean and destroy rf_hamiltonian_type datastructure
!!
!! SIDE EFFECTS
!!  rf_Ham<rf_hamiltonian_type>=All dynamic memory defined in the structure is deallocated.
!!
!! PARENTS
!!      nstpaw3,nstwf3,rhofermi3,vtorho3
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine destroy_rf_hamiltonian(rf_Ham)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_rf_hamiltonian'
!End of the abilint section

 type(rf_hamiltonian_type),intent(inout) :: rf_Ham

! *************************************************************************

 !@rf_hamiltonian_type

! Integer pointers.
 if (allocated(rf_Ham%indlmn_typ))  then
   ABI_DEALLOCATE(rf_Ham%indlmn_typ)
 end if

! Real pointers
 if (allocated(rf_Ham%e1kbfr))   then
   ABI_DEALLOCATE(rf_Ham%e1kbfr)
 end if
 if (allocated(rf_Ham%e1kbsc))   then
   ABI_DEALLOCATE(rf_Ham%e1kbsc)
 end if
 if (allocated(rf_Ham%ekb_typ))   then
   ABI_DEALLOCATE(rf_Ham%ekb_typ)
 end if
 if (allocated(rf_Ham%sij_typ))   then
   ABI_DEALLOCATE(rf_Ham%sij_typ)
 end if
 if (allocated(rf_Ham%pspso_typ))   then
   ABI_DEALLOCATE(rf_Ham%pspso_typ)
 end if

end subroutine destroy_rf_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/init_rf_hamiltonian
!! NAME
!!  init_rf_hamiltonian
!!
!! FUNCTION
!!  Creation method for the rf_hamiltonian_type structure.
!!  It allocates memory and initializes all quantities that do not depend on the k-point or spin.
!!
!! INPUTS
!!  cplex_paw=1 if all on-site PAW quantities are real (GS), 2 if they are complex (RF)
!!  gs_Ham<gs_hamiltonian_type>=Structured datatype containing data for ground-state Hamiltonian at (k+q)
!!  has_e1kbsc= -optional- if 1, rf_Ham%e1kbsc is allocated (if necessary); otherwise not allocated.
!!             e1kbsc contains the self-consistent part of 1st-order PAW Dij coefficients.
!!  ipert=index of perturbation
!!  nspinor=Number of spinorial components
!!  typat(natom)=Type of each atom
!!
!! SIDE EFFECTS
!!  rf_Ham<rf_hamiltonian_type>=Structured datatype almost completely initialized:
!!   * Basic variables and dimensions are transfered to the structure.
!!   * All pointers are allocated with correct dimensions.
!!   * Quantities that do not depend on the k-point or spin are initialized.
!!
!! PARENTS
!!      nstpaw3,nstwf3,rhofermi3,vtorho3
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine init_rf_hamiltonian(cplex,gs_Ham,ipert,nspinor,rf_Ham,typat,&
&                              has_e1kbsc) ! optional argument

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_rf_hamiltonian'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ipert,nspinor
 integer,intent(in),optional :: has_e1kbsc
 type(gs_hamiltonian_type),intent(in) :: gs_Ham
 type(rf_hamiltonian_type),intent(inout) :: rf_Ham
!arrays
 integer,intent(in) :: typat(gs_Ham%natom)

!Local variables-------------------------------
!scalars
 integer :: cplex_dij1,has_e1kbsc_,ispden,itypat

! *************************************************************************

 !@rf_hamiltonian_type

 has_e1kbsc_=0;if (present(has_e1kbsc)) has_e1kbsc_=has_e1kbsc

 rf_Ham%usepaw  =gs_Ham%usepaw
 rf_Ham%lmnmax  =gs_Ham%lmnmax
 rf_Ham%nspinor =gs_Ham%nspinor
 rf_Ham%dimekb1 =gs_Ham%dimekb1
 rf_Ham%dimekb2 =gs_Ham%dimekb2

!Additional dimension in case of PAW
 if (rf_Ham%usepaw==0) then
   rf_Ham%dime1kb=0
 else
   if (ipert/=gs_Ham%natom+1.and.ipert/=gs_Ham%natom+5) then
     cplex_dij1=max(cplex,nspinor)
     rf_Ham%dime1kb=cplex_dij1*(rf_Ham%lmnmax*(rf_Ham%lmnmax+1))/2
   else
     rf_Ham%dime1kb=0
   end if
 end if

!Allocate the arrays of the 1st-order Hamiltonian
 if (ipert>=1 .and. ipert<=gs_Ham%natom) then
   ABI_ALLOCATE(rf_Ham%pspso_typ,(1))
   ABI_ALLOCATE(rf_Ham%indlmn_typ,(6,rf_Ham%lmnmax,1))
   ABI_ALLOCATE(rf_Ham%ekb_typ,(rf_Ham%dimekb1,1,rf_Ham%nspinor**2))
   if (rf_Ham%usepaw==1) then
     ABI_ALLOCATE(rf_Ham%sij_typ,(rf_Ham%dimekb1,1))
   end if
 end if

 if ((ipert>=1.and.ipert<=gs_Ham%natom).or.(ipert==gs_Ham%natom+2)&
&  .or.(ipert==gs_Ham%natom+3).or.(ipert==gs_Ham%natom+4)) then
   if (rf_Ham%usepaw==1.and.rf_Ham%dime1kb>0) then
     ABI_ALLOCATE(rf_Ham%e1kbfr,(rf_Ham%dime1kb,rf_Ham%dimekb2,rf_Ham%nspinor**2))
     if (has_e1kbsc_==1) then
       ABI_ALLOCATE(rf_Ham%e1kbsc,(rf_Ham%dime1kb,rf_Ham%dimekb2,rf_Ham%nspinor**2))
     end if
   end if
 end if

!Initialize arrays not depending on k point
 if (ipert>=1 .and. ipert<=gs_Ham%natom) then
   itypat=typat(ipert)
   rf_Ham%pspso_typ(1)=gs_Ham%pspso(itypat)
   rf_ham%indlmn_typ(:,:,1)=gs_Ham%indlmn(:,:,itypat)
   if (rf_Ham%usepaw==0) then
     do ispden=1,rf_Ham%nspinor**2
       rf_Ham%ekb_typ(:,1,ispden)=gs_Ham%ekb(:,itypat,ispden)
     end do
   else
     rf_Ham%sij_typ(:,1)=gs_Ham%sij(:,itypat)
   end if
 end if

end subroutine init_rf_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/load_paw_rf_hamiltonian
!! NAME
!!  load_paw_rf_hamiltonian
!!
!! FUNCTION
!!  Setup of the k-dependent part of the 1st-order Hamiltonian.
!!
!! INPUTS
!!  gs_Ham<gs_hamiltonian_type>=Structured datatype containing data for ground-state Hamiltonian at (k+q)
!!  ipert=index of perturbation
!!  isppol=index of current spin
!!  paw_ij1(natom*usepaw)<paw_ij_type>=Various 1st-order arrays given on (i,j) (partial waves)
!!                                     channels (paw_ij1%dij and paw_ij1%difr only used here).
!!  typat(natom)=Type of each atom
!!
!! SIDE EFFECTS
!!  rf_Ham<rf_hamiltonian_type>=data for first-order hamiltonian at (k,q)
!!   * Quantities that do depend on the k-point or spin are initialized.
!!
!! PARENTS
!!      rhofermi3,vtorho3
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine load_paw_rf_hamiltonian(gs_Ham,ipert,isppol,paw_ij1,rf_Ham,typat, &
&          mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'load_paw_rf_hamiltonian'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ipert,isppol
 integer,optional,intent(in) :: mpi_comm_atom
 type(gs_hamiltonian_type),intent(in) :: gs_Ham
 type(rf_hamiltonian_type),intent(inout) :: rf_Ham
!arrays
 integer,intent(in) :: typat(gs_Ham%natom)
 integer,optional,target,intent(in)  :: mpi_atmtab(:)
 type(paw_ij_type),intent(in) :: paw_ij1(:)

!Local variables-------------------------------
!scalars
 integer :: iatom,ispden,my_comm_atom,my_natom
 logical :: load_e1kbfr,load_e1kbsc,my_atmtab_allocated,paral_atom
!scalars
 integer,pointer :: my_atmtab(:)

! *************************************************************************

 !@rf_hamiltonian_type
 if (rf_Ham%usepaw==0) return

 if (allocated(rf_Ham%ekb_typ)) then
   iatom=typat(ipert);if (rf_Ham%usepaw==1) iatom=ipert
   do ispden=1,rf_Ham%nspinor**2
     rf_Ham%ekb_typ(:,1,ispden)=gs_Ham%ekb(:,iatom,ispden)
   end do
 end if

 load_e1kbfr=allocated(rf_Ham%e1kbfr)
 load_e1kbsc=allocated(rf_Ham%e1kbsc)
 if ((.not.load_e1kbfr).and.(.not.load_e1kbsc)) return

!Set up parallelism over atoms
 my_natom=size(paw_ij1);my_comm_atom=xmpi_self
 paral_atom=(present(mpi_comm_atom).and.my_natom/=gs_Ham%natom)
 my_comm_atom=xmpi_self;if (paral_atom) my_comm_atom=mpi_comm_atom
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,gs_ham%natom,my_natom_ref=my_natom)

!Retrieve PAW Dij(1) coefficients for this spin component
 if (load_e1kbfr.and.load_e1kbsc) then
   call pawdij2e1kb(paw_ij1,rf_Ham%dime1kb,isppol,my_atmtab,my_comm_atom,rf_Ham%nspinor,&
&                   e1kbfr=rf_Ham%e1kbfr,e1kbsc=rf_Ham%e1kbsc)
 else if (load_e1kbfr) then
   call pawdij2e1kb(paw_ij1,rf_Ham%dime1kb,isppol,my_atmtab,my_comm_atom,rf_Ham%nspinor,&
&                   e1kbfr=rf_Ham%e1kbfr)
 else if (load_e1kbsc) then
   call pawdij2e1kb(paw_ij1,rf_Ham%dime1kb,isppol,my_atmtab,my_comm_atom,rf_Ham%nspinor,&
&                   e1kbsc=rf_Ham%e1kbsc)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine load_paw_rf_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/pawdij2ekb
!! NAME
!!  pawdij2ekb
!!
!! FUNCTION
!!  Transfer PAW Dij (on-site GS Hamiltonian) values
!!  from paw_ij datastructure to ekb array
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hamiltonian
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine pawdij2ekb(ekb,paw_ij,dimekb1,isppol,mpi_atmtab,mpi_comm_atom,nspinor)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdij2ekb'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimekb1,isppol,mpi_comm_atom,nspinor
!arrays
 integer,pointer,intent(in)  :: mpi_atmtab(:)
 real(dp),intent(out) :: ekb(:,:,:)
 type(paw_ij_type),intent(in) :: paw_ij(:)

!Local variables-------------------------------
!scalars
 integer :: dimdij,iatom,iatom_tot,ierr,isp,ispden,my_natom
 logical :: paral_atom

! *************************************************************************

 ekb=zero
 paral_atom=(xcomm_size(mpi_comm_atom)>1)
 my_natom=size(paw_ij)

 if (my_natom>0) then
   if (allocated(paw_ij(1)%dij)) then
     do ispden=1,nspinor**2
       isp=isppol; if (nspinor==2) isp=ispden
       do iatom=1,my_natom
         iatom_tot=iatom;if (paral_atom) iatom_tot=mpi_atmtab(iatom)
         dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
         if (dimdij>dimekb1) then
           MSG_BUG(' size of paw_ij%dij>dimekb1 !')
         end if
         ekb(1:dimdij,iatom_tot,ispden)=paw_ij(iatom)%dij(1:dimdij,isp)
       end do
     end do
   end if
 end if

!Communication in case of distribution over atomic sites
 if (paral_atom) then
   call xmpi_sum(ekb,mpi_comm_atom,ierr)
 end if

end subroutine pawdij2ekb
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/pawdij2e1kb
!! NAME
!!  pawdij2e1kb
!!
!! FUNCTION
!!  Transfer PAW Dij (on-site RF Hamiltonian) values
!!  from paw_ij datastructure to e1kb array
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hamiltonian,nstpaw3
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine pawdij2e1kb(paw_ij1,dime1kb,isppol,mpi_atmtab,mpi_comm_atom,nspinor,e1kbfr,e1kbsc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdij2e1kb'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dime1kb,isppol,mpi_comm_atom,nspinor
!arrays
 integer,intent(in)  :: mpi_atmtab(:)
 real(dp),optional,intent(out) :: e1kbfr(:,:,:),e1kbsc(:,:,:)
 type(paw_ij_type),intent(in) :: paw_ij1(:)

!Local variables-------------------------------
!scalars
 integer :: dimdij1,iatom,iatom_tot,ierr,isp,ispden,my_natom
 logical :: paral_atom

! *************************************************************************

 if (present(e1kbfr)) e1kbfr=zero
 if (present(e1kbsc)) e1kbsc=zero
 paral_atom=(xcomm_size(mpi_comm_atom)>1)
 my_natom=size(paw_ij1)

 if (my_natom>0.and.present(e1kbfr)) then
   if (allocated(paw_ij1(1)%dijfr)) then
     do ispden=1,nspinor**2
       isp=isppol;if (nspinor==2) isp=ispden
       do iatom=1,my_natom
         iatom_tot=iatom;if (paral_atom) iatom_tot=mpi_atmtab(iatom)
         dimdij1=paw_ij1(iatom)%cplex_dij*paw_ij1(iatom)%lmn2_size
         if (dimdij1>dime1kb) then
           MSG_BUG(' size of paw_ij1%dij>dime1kb !')
         end if
         e1kbfr(1:dimdij1,iatom_tot,ispden)=paw_ij1(iatom)%dijfr(1:dimdij1,isp)
       end do
     end do
   end if
 end if

 if (my_natom>0.and.present(e1kbsc)) then
   if (allocated(paw_ij1(1)%dijfr).and.allocated(paw_ij1(1)%dij)) then
     do ispden=1,nspinor**2
       isp=isppol;if (nspinor==2) isp=ispden
       do iatom=1,my_natom   
         iatom_tot=iatom;if (paral_atom) iatom_tot=mpi_atmtab(iatom)
         dimdij1=paw_ij1(iatom)%cplex_dij*paw_ij1(iatom)%lmn2_size
         if (dimdij1>dime1kb) then
           MSG_BUG(' size of paw_ij1%dij>dime1kb !')
         end if
         e1kbsc(1:dimdij1,iatom_tot,ispden)=paw_ij1(iatom)%dij  (1:dimdij1,isp) &
&                                          -paw_ij1(iatom)%dijfr(1:dimdij1,isp)
       end do
     end do
   end if
 end if
 
!Communication in case of distribution over atomic sites
 if (paral_atom) then
   if (present(e1kbfr)) then
     call xmpi_sum(e1kbfr,mpi_comm_atom,ierr)
   end if
   if (present(e1kbsc)) then
     call xmpi_sum(e1kbsc,mpi_comm_atom,ierr)
   end if
 end if

end subroutine pawdij2e1kb
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/init_ddiago_ctl
!! NAME
!!  init_ddiago_ctl
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_shirley,outkss
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine init_ddiago_ctl(Dctl,jobz,isppol,nspinor,ecut,kpoint,nloalg,gmet,&
& nband_k,istwf_k,ecutsm,effmass,abstol,range,ilu,vlu,use_scalapack,prtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_ddiago_ctl'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol,nspinor
 integer,optional,intent(in) :: istwf_k,prtvol,use_scalapack,nband_k
 real(dp),intent(in) :: ecut
 real(dp),optional,intent(in) :: ecutsm,effmass
 real(dp),optional,intent(in) :: abstol
 character(len=*),intent(in) :: jobz
 character(len=*),optional,intent(in) :: range
 type(ddiago_ctl_type),intent(out) :: Dctl
!arrays
 integer,intent(in) :: nloalg(5)
 integer,optional,intent(in) :: ilu(2)
 real(dp),intent(in) :: kpoint(3)
 real(dp),optional,intent(in) :: vlu(2)
 real(dp),intent(in) :: gmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: npw_k
 logical :: ltest
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,allocatable :: kg_k(:,:)

! *************************************************************************

 call initmpi_seq(MPI_enreg_seq) ! Fake MPI_type.

 Dctl%isppol  = isppol
 Dctl%nspinor = nspinor
 Dctl%kpoint  = kpoint

 if (PRESENT(istwf_k)) then
  Dctl%istwf_k = istwf_k
 else
  Dctl%istwf_k = set_istwfk(kpoint)
 end if

 ABI_CHECK(Dctl%istwf_k==1,"istwf_k/=1 not coded")

 Dctl%jobz   = toupper(jobz(1:1))
 Dctl%range  = "A"
 if (PRESENT(range)) then
  Dctl%range = toupper(range)
 end if

 Dctl%ecut = ecut
 Dctl%ecutsm = zero; if (PRESENT(ecutsm)) Dctl%ecutsm = ecutsm
 Dctl%effmass = one; if (PRESENT(effmass)) Dctl%effmass = effmass
 Dctl%nloalg  = nloalg
 Dctl%prtvol = 0; if (PRESENT(prtvol)) Dctl%prtvol = prtvol
 Dctl%abstol = -tol8; if (PRESENT(abstol)) Dctl%abstol = abstol

 ABI_ALLOCATE(kg_k,(3,0))

! * Total number of G-vectors for this k-point with istwf_k=1.
 call kpgsph(ecut,0,gmet,0,0,1,kg_k,kpoint,0,MPI_enreg_seq,0,Dctl%npwtot)

! * G-vectors taking into account time-reversal symmetry.
 call kpgsph(ecut,0,gmet,0,0,istwf_k,kg_k,kpoint,0,MPI_enreg_seq,0,npw_k)

 Dctl%npw_k = npw_k
 ABI_DEALLOCATE(kg_k)

 Dctl%do_full_diago = .FALSE.

 SELECT CASE (Dctl%range)
  CASE ("A")

  ! Check on the number of stored bands.
  Dctl%nband_k=-1
  if (PRESENT(nband_k)) Dctl%nband_k=nband_k

  if (Dctl%nband_k==-1.or.Dctl%nband_k>=npw_k*nspinor) then
    Dctl%nband_k=npw_k*nspinor
    write(msg,'(4a)')ch10,&
&    ' Since the number of bands to be computed was (-1) or',ch10,&
&    ' too large, it has been set to the max. value npw_k*nspinor. '
    if (Dctl%prtvol>0) then
      call wrtout(std_out,msg,'COLL')
    end if
  else
    Dctl%nband_k=nband_k
  end if

  Dctl%do_full_diago = (Dctl%nband_k==npw_k*nspinor)

  if (Dctl%do_full_diago) then
    write(msg,'(6a)')ch10,&
&    ' Since the number of bands to be computed',ch10,&
&    ' is equal to the number of G-vectors found for this k-point,',ch10,&
&    ' the program will perform complete diagonalization.'
  else
    write(msg,'(6a)')ch10,&
&     ' Since the number of bands to be computed',ch10,&
&     ' is less than the number of G-vectors found,',ch10,&
&     ' the program will perform partial diagonalization.'
  end if
  if (Dctl%prtvol>0) then
    call wrtout(std_out,msg,'COLL')
  end if

 CASE ("I")
  if (.not.PRESENT(ilu)) then
    MSG_ERROR(" ilu must be specified when range=I ")
  end if
  Dctl%ilu = ilu

  ltest = ( ( ilu(2)>=ilu(1) ) .and. ilu(1)>=1 .and. ilu(2)<=Dctl%npwtot )
  write(msg,'(a,2i0)')" Illegal value for ilu: ",ilu
  ABI_CHECK(ltest,msg)
  Dctl%nband_k= ilu(2)-ilu(1)+1

 CASE ("V")
  if (.not.PRESENT(vlu)) then
    MSG_ERROR(" vlu must be specified when range=V ")
  end if
  Dctl%vlu = vlu

  Dctl%nband_k=-1 !??

  ltest = (vlu(2)>vlu(1))
  write(msg,'(a,2f0.3)')" Illegal value for vlu: ",vlu
  ABI_CHECK(ltest,msg)

 CASE DEFAULT
   msg = " Unknown value for range: "//TRIM(Dctl%range)
   MSG_ERROR(msg)
 END SELECT

 ! Consider the case in which we asked for the entire set of eigenvectors
 ! but the number of bands is less that npw_k. Therefore have to prepare the call to ZHEEVX.
 ! TODO this has to be done in a cleaner way.
 if (Dctl%range=="A".and. (.not.Dctl%do_full_diago)) then
   Dctl%range="I"
   Dctl%ilu(1) = 1
   Dctl%ilu(2) = npw_k*nspinor
   Dctl%nband_k= npw_k*nspinor
 end if

 Dctl%use_scalapack=0
 if (PRESENT(use_scalapack)) then
   Dctl%use_scalapack=use_scalapack
 end if
 ABI_CHECK(Dctl%use_scalapack==0," scalapack mode not coded")

 call destroy_mpi_enreg(MPI_enreg_seq)

end subroutine init_ddiago_ctl
!!***

!----------------------------------------------------------------------

END MODULE m_hamiltonian
!!***
