!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_denspot_set
!! NAME
!!  wvl_denspot_set
!!
!! FUNCTION
!!  Fill in denspot datatype with information related
!!  to density and potential data.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (DC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      allocaterhopot,density_descriptors,dpbox_set
!!      initialize_dft_local_fields,wrtout,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_denspot_set(den,gth_params,ixc,&
& me,natom,nproc,nsppol,rprimd,wvl,&
& wvl_crmult,wvl_frmult,xred)
    
 use defs_basis
 use defs_datatypes
! use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_xmpi

#if defined HAVE_DFT_BIGDFT
  use BigDFT_API,only: initialize_DFT_local_fields,&
allocateRhoPot,input_variables
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_denspot_set'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
  integer,intent(in):: ixc,me,natom,nproc,nsppol
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(in) :: wvl_frmult,wvl_crmult
  real(dp), intent(inout)  :: xred(3,natom)
  type(wvl_denspot_type), intent(out) :: den
  type(wvl_internal_type),intent(in)  :: wvl
  type(pseudopotential_gth_type),intent(in)::gth_params

!Local variables-------------------------------
  real(dp), allocatable :: xcart(:,:)

  character(len=3),parameter :: rho_commun='DBL'
  character(len=500) :: message

#if defined HAVE_DFT_BIGDFT
  ! To be removed, waiting for BigDFT upgrade.
  type(input_variables) :: in
  type(local_zone_descriptors) :: Lzd
#endif

  ! *************************************************************************
 
!DEBUG
!write (std_out,*) ' wvl_denspot_set : enter'
!ENDDEBUG
 
 write(message, '(a,a)' ) ch10,&
& ' wvl_denspot_set: Create wavelet type denspot.'
 call wrtout(std_out,message,'COLL')

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, natom))
 call xred2xcart(natom, rprimd, xcart, xred)

#if defined HAVE_DFT_BIGDFT
 call initialize_DFT_local_fields(den%denspot)

!number of planes for the density
!dpbox%nscatterarr(jproc, 1) = ngfft3_density
!number of planes for the potential
!dpbox%nscatterarr(jproc, 2) = ngfft3_potential
!starting offset for the potential
!dpbox%nscatterarr(jproc, 3) = density_start + potential_shift - 1
!GGA XC shift between density and potential
!dpbox%nscatterarr(jproc, 4) = potential_shift

!DEBUG
 write(std_out,*) 'wvl_denspot_set: TODO, update BigDFT dpbox_set()'
!ENDDEBUG
 ! To be removed !!!!!
 Lzd%hgrids = wvl%h
 Lzd%Glr%d = wvl%Glr%d
 in%PSolver_groupsize = 0
 in%ixc = ixc
 in%nspin = nsppol
 in%SIC%approach = "NONE"
 ! To be removed !!!!!
 call dpbox_set(den%denspot%dpbox,Lzd,me,nproc,xmpi_world,in,wvl%atoms%astruct%geocode)

!here dpbox can be put as input
 call density_descriptors(me,nproc,nsppol,wvl_crmult,wvl_frmult,wvl%atoms,&
 den%denspot%dpbox,rho_commun,xcart,gth_params%radii_cf,den%denspot%rhod)

!allocate the arrays.
!"in" objects required: in%spin
!Note: change allocateRhoPot
 call allocateRhoPot(me,wvl%Glr,nsppol,wvl%atoms,xcart,den%denspot)

!Aditional informations.
 den%symObj = wvl%atoms%astruct%sym%symObj
#endif

 ABI_DEALLOCATE(xcart)

!DEBUG
!write (std_out,*) ' wvl_denspot_set : exit'
!stop
!ENDDEBUG

end subroutine wvl_denspot_set
!!***
