!{\src2tex{textfont=tt}}
!!****f* ABINIT/nullify_wvl_data
!! NAME
!!  nullify_wvl_data
!!
!! FUNCTION
!!  Nullify all wvl pointers
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2014 ABINIT group (T. Rangel)
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
!!      gstate
!!
!! CHILDREN
!!      nullify_atoms_data,nullify_collective_comms
!!      nullify_communications_arrays,nullify_dft_local_fields
!!      nullify_diis_objects,nullify_gaussian_basis,nullify_gpu_pointers
!!      nullify_local_zone_descriptors,nullify_locreg_descriptors
!!      nullify_orbitals_data,nullify_p2pcomms,nullify_paw_objects
!!      nullify_rholoc_objects
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nullify_wvl_data(wvl)
    
 use m_profiling_abi
 use m_errors
 use defs_basis
 use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : &
& nullify_gaussian_basis, nullify_local_zone_descriptors,&
& nullify_orbitals_data, nullify_communications_arrays,&
& nullify_diis_objects, nullify_p2pcomms,&
& nullify_GPU_pointers, nullify_locreg_descriptors,&
& nullify_atoms_data, nullify_rholoc_objects, nullify_paw_objects,&
& nullify_DFT_local_fields
!& nullify_foe_data, &
!& nullify_hamiltonian_descriptors,&
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_wvl_data'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wvl_data) , intent(inout)  :: wvl

!Local variables-------------------------------
!character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
 DBG_ENTER("COLL")

!1)   wvl_projectors_type: projectors
 nullify(wvl%projectors%proj)

!1.1) nonlocal_psp_descriptors: nlpspd
#if defined HAVE_DFT_BIGDFT
 nullify(wvl%projectors%nlpspd%plr)
#endif

!2)   wvl_wf_type: wfs
#if defined HAVE_DFT_BIGDFT
!2.1) DFT_wavefunction: ks
 nullify(wvl%wfs%ks%psi)
 nullify(wvl%wfs%ks%hpsi)
 nullify(wvl%wfs%ks%psit)
 nullify(wvl%wfs%ks%spsi)
 nullify(wvl%wfs%ks%gaucoeffs)
 nullify(wvl%wfs%ks%confdatarr)
 nullify(wvl%wfs%ks%oldpsis)
 nullify(wvl%wfs%ks%coeff)
!2.1.1) gaussian_basis: gbd 
 call nullify_gaussian_basis(wvl%wfs%ks%gbd)
!2.1.2) local_zone_descriptors: Lzd
 call nullify_local_zone_descriptors(wvl%wfs%ks%Lzd)
!2.1.3) orbitals_data: orbs
 call nullify_orbitals_data(wvl%wfs%ks%orbs)
!2.1.4) communications_arrays: comms
 call nullify_communications_arrays(wvl%wfs%ks%comms)
!2.1.5) diis_objects: diis
 call nullify_diis_objects(wvl%wfs%ks%diis)
!2.1.6) p2pComms: comgp
 call nullify_p2pcomms(wvl%wfs%ks%comgp)
!2.1.7) collective_comms: collcom
 call nullify_collective_comms(wvl%wfs%ks%collcom)
!2.1.8) collective_comms: collcom_sr
 call nullify_collective_comms(wvl%wfs%ks%collcom_sr)
!2.1.9) foe_data: foe_obj
! call nullify_foe_data(wvl%wfs%ks%foe_obj)
!2.1.10) linear matrices: 
!2.1.11) hamiltonian_descriptors: ham_descr
! call nullify_hamiltonian_descriptors(wvl%wfs%ks%ham_descr)

!2.2) GPU_pointers: GPU
 call nullify_GPU_pointers(wvl%wfs%GPU)
#endif

!3) wvl_internal_type: descr
#if defined HAVE_DFT_BIGDFT
!3.1) locreg_descriptors: Glr
 call nullify_locreg_descriptors(wvl%descr%Glr)
!3.2) atoms_data: atoms
 call nullify_atoms_data(wvl%descr%atoms)
!3.3) rholoc_objects: rholoc
 call nullify_rholoc_objects(wvl%descr%rholoc)
!3.4) paw_objects: paw 
 call nullify_paw_objects(wvl%descr%paw)
#endif

!4) wvl_denspot_type: den
#if defined HAVE_DFT_BIGDFT
!4.1) DFT_local_fields: denspot
 call nullify_DFT_local_fields(wvl%den%denspot)
#endif

!5) wvl_energy_terms: e
!there are no pointers here

 DBG_EXIT("COLL")

end subroutine nullify_wvl_data
!!***
