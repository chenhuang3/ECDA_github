!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_seq
!! NAME
!!  initmpi_seq
!!
!! FUNCTION
!!  Initializes the MPI information for a sequential use of other routines.
!!
!! COPYRIGHT
!!  Copyright (C) 2004-2014 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      atm2fft,atm2fft3,bethe_salpeter,bsepostproc,calc_vhxc_me,cut3d
!!      debug_tools,dieltcel,fftprof,ks_ddiago,kss2wfk,linear_optics_paw,m_fft
!!      m_fft_prof,m_fftcore,m_gsphere,m_hamiltonian,m_io_kss,m_kxc,m_mpinfo
!!      m_paw_pwij,m_ppmodel,m_screening,m_wfk,m_wfs,make_efg_el,mlwfovlp_qp
!!      mrgddb,mrggkk,mrgscr,optic,overlap_wf,pawmknhat,pawmknhat_psipsi
!!      pawsushat,scfcv,screening,sigma,suscep_stat,susk,suskmm,ujdet
!!      vdw_kernelgen,wffile,wfread
!!
!! CHILDREN
!!      nullify_mpi_enreg
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_seq(mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_mpinfo, only : nullify_mpi_enreg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_seq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(MPI_type),intent(out) :: mpi_enreg

! ***********************************************************************

 DBG_ENTER("COLL")

!Set default seq values for scalars
 mpi_enreg%bandpp=1
 mpi_enreg%me=0
 mpi_enreg%me_band=0
 mpi_enreg%me_cell=0
 mpi_enreg%me_fft=0
 mpi_enreg%me_g0=0
 mpi_enreg%me_img=0
 mpi_enreg%me_hf=0
 mpi_enreg%me_kpt=0
 mpi_enreg%me_pert=0
 mpi_enreg%me_spinor=0
 mpi_enreg%me_wvl=0
 mpi_enreg%my_natom=0       ! Should be natom
 mpi_enreg%my_isppoltab=0   ! Should be (1,0) if nsppol=1 or (1,1) if nsppol=2
 mpi_enreg%ngfft3_ionic=1
 mpi_enreg%my_nimage=1
 mpi_enreg%nproc=1
 mpi_enreg%nproc_atom=1
 mpi_enreg%nproc_band=1
 mpi_enreg%nproc_cell=1
 mpi_enreg%nproc_fft=1
 mpi_enreg%nproc_img=1
 mpi_enreg%nproc_hf=1
 mpi_enreg%nproc_kpt=1
 mpi_enreg%nproc_pert=1
 mpi_enreg%nproc_spinor=1
 mpi_enreg%nproc_wvl=1
 mpi_enreg%paralbd=0
 mpi_enreg%paral_img=0
 mpi_enreg%paral_hf=0
 mpi_enreg%paral_kgb=0
 mpi_enreg%paral_pert=0
 mpi_enreg%paral_spinor=0
 
!Set default seq values for communicators
 mpi_enreg%comm_world          = xmpi_world
 mpi_enreg%comm_atom           = xmpi_self
 mpi_enreg%comm_band           = xmpi_self
 mpi_enreg%comm_bandspinor     = xmpi_self
 mpi_enreg%comm_bandfft        = xmpi_self
 mpi_enreg%comm_bandspinorfft  = xmpi_self
 mpi_enreg%comm_cell           = xmpi_self
 mpi_enreg%comm_cell_pert      = xmpi_self
 mpi_enreg%comm_fft            = xmpi_self
 mpi_enreg%comm_hf             = xmpi_self 
 mpi_enreg%comm_img            = xmpi_self 
 mpi_enreg%comm_kpt            = xmpi_self
 mpi_enreg%comm_kptband        = xmpi_self
 mpi_enreg%comm_pert           = xmpi_self
 mpi_enreg%comm_spinor         = xmpi_self
 mpi_enreg%comm_spinorfft      = xmpi_self
 mpi_enreg%comm_wvl            = xmpi_self

!Nullify all pointers
 call nullify_mpi_enreg(mpi_enreg)

!Allocate and nullify distribfft datastructure
! This is not good since distribfft is not initialized here (even with 0s). 
! It can be dangerous if use with no care (Valgrind might complain)
 ABI_DATATYPE_ALLOCATE(mpi_enreg%distribfft,)

 DBG_EXIT("COLL")

end subroutine initmpi_seq
!!***
