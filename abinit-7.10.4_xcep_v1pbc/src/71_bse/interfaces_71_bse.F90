!!****m* ABINIT/interfaces_71_bse
!! NAME
!! interfaces_71_bse
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/71_bse
!!
!! COPYRIGHT
!! Copyright (C) 2010-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_71_bse

 implicit none

interface
 subroutine exc_build_block(BSp,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,Wfd,W,Hdr_bse,&  
  &  nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff,rhxtwg_q0,is_resonant,fname)
  use m_vcoul
  use m_pawtab
  use m_wfs
  use defs_basis
  use m_bz_mesh
  use m_pawang
  use defs_datatypes
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  use m_gsphere
  use m_screen
  use m_paw_pwij
  implicit none
  integer,intent(in) :: nfftot_osc
  type(excparam),intent(in) :: BSp
  type(crystal_t),intent(in) :: Cryst
  type(gsphere_t),intent(in) :: Gsph_c
  type(gsphere_t),intent(in) :: Gsph_x
  type(hdr_type),intent(inout) :: Hdr_bse
  type(kmesh_t),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(kmesh_t),intent(in) :: Qmesh
  type(vcoul_t),intent(in) :: Vcp
  type(screen_t),intent(inout) :: W
  type(wfd_t),target,intent(inout) :: Wfd
  character(len=*),intent(in) :: fname
  logical,intent(in) :: is_resonant
  integer,intent(in) :: ngfft_osc(18)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Wfd%usepaw)
  integer,intent(in) :: ktabr(nfftot_osc,Kmesh%nbz)
  complex(gwpc),intent(in) :: rhxtwg_q0(BSp%npweps,BSp%lomo_min:BSp%humo_max, &
  &         BSp%lomo_min:BSp%humo_max,Wfd%nkibz,Wfd%nsppol)
 end subroutine exc_build_block
end interface

interface
 subroutine exc_build_ham(BSp,BS_files,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,&  
  &  Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff)
  use m_vcoul
  use m_pawtab
  use m_wfs
  use m_bz_mesh
  use m_pawang
  use defs_datatypes
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  use m_gsphere
  use m_screen
  use m_paw_pwij
  implicit none
  integer,intent(in) :: nfftot_osc
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(crystal_t),intent(in) :: Cryst
  type(gsphere_t),intent(in) :: Gsph_c
  type(gsphere_t),intent(in) :: Gsph_x
  type(hdr_type),intent(inout) :: Hdr_bse
  type(kmesh_t),intent(in) :: Kmesh
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(kmesh_t),intent(in) :: Qmesh
  type(vcoul_t),intent(in) :: Vcp
  type(screen_t),intent(inout) :: W
  type(wfd_t),target,intent(inout) :: Wfd
  integer,intent(in) :: ngfft_osc(18)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Wfd%usepaw)
  integer,intent(in) :: ktabr(nfftot_osc,Kmesh%nbz)
 end subroutine exc_build_ham
end interface

interface
 subroutine wfd_all_mgq0(Wfd,Cryst,Qmesh,Gsph_x,Vcp,&  
  &  Psps,Pawtab,Paw_pwff,lomo_spin,homo_spin,humo_spin,nfftot_osc,ngfft_osc,npweps,mgq0)
  use m_vcoul
  use m_pawtab
  use m_wfs
  use m_bz_mesh
  use defs_basis
  use m_crystal
  use m_gsphere
  use defs_datatypes
  use m_paw_pwij
  implicit none
  integer,intent(in) :: nfftot_osc
  integer,intent(in) :: npweps
  type(crystal_t),intent(in) :: Cryst
  type(gsphere_t),intent(in) :: Gsph_x
  type(pseudopotential_type),intent(in) :: Psps
  type(kmesh_t),intent(in) :: Qmesh
  type(vcoul_t),intent(in) :: Vcp
  type(wfd_t),target,intent(inout) :: Wfd
  integer,intent(in) :: ngfft_osc(18)
  type(paw_pwff_type),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  integer,intent(in) :: homo_spin(Wfd%nsppol)
  integer,intent(in) :: humo_spin(Wfd%nsppol)
  integer,intent(in) :: lomo_spin(Wfd%nsppol)
  complex(gwpc),allocatable,intent(out) :: mgq0(:,:,:,:,:)
 end subroutine wfd_all_mgq0
end interface

interface
 subroutine exc_den(BSp,BS_files,ngfft,nfftot,Kmesh,ktabr,Wfd)
  use m_bz_mesh
  use m_bs_defs
  use m_wfs
  implicit none
  integer,intent(in) :: nfftot
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(kmesh_t),intent(in) :: Kmesh
  type(wfd_t),intent(inout) :: Wfd
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ktabr(nfftot,BSp%nkbz)
 end subroutine exc_den
end interface

interface
 subroutine exc_plot(Bsp,Bs_files,Wfd,Kmesh,Cryst,Psps,Pawtab,Pawrad,paw_add_onsite,spin_opt,which_fixed,eh_rcoord,nrcell,ngfftf)
  use m_pawrad
  use m_wfs
  use m_bz_mesh
  use defs_basis
  use m_bs_defs
  use m_crystal
  use m_pawtab
  use defs_datatypes
  implicit none
  integer,intent(in) :: spin_opt
  integer,intent(in) :: which_fixed
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: Bsp
  type(crystal_t),intent(in) :: Cryst
  type(kmesh_t),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(wfd_t),intent(inout) :: Wfd
  logical,intent(in) :: paw_add_onsite
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: nrcell(3)
  type(pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
  real(dp),intent(in) :: eh_rcoord(3)
 end subroutine exc_plot
end interface

interface
 subroutine build_spectra(BSp,BS_files,Cryst,Kmesh,KS_BSt,QP_BSt,Psps,Pawtab,Wfd,Hur,drude_plsmf,comm)
  use defs_basis
  use m_wfs
  use m_bz_mesh
  use m_bs_defs
  use m_crystal
  use m_paw_commutator
  use m_pawtab
  use defs_datatypes
  implicit none
  integer,intent(in) :: comm
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  type(crystal_t),intent(in) :: Cryst
  type(ebands_t),intent(in) :: KS_BSt
  type(kmesh_t),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(ebands_t),intent(in) :: QP_BSt
  type(wfd_t),intent(inout) :: Wfd
  real(dp),intent(in) :: drude_plsmf
  type(hur_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 end subroutine build_spectra
end interface

interface
 subroutine exc_write_data(BSp,BS_files,what,eps,dos)
  use defs_basis
  use m_bs_defs
  implicit none
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  character(len=*),intent(in) :: what
  real(dp),optional,intent(in) :: dos(BSp%nomega)
  complex(dpc),intent(in) :: eps(BSp%nomega,BSp%nq)
 end subroutine exc_write_data
end interface

interface
 subroutine exc_eps_rpa(nbnds,lomo_spin,lomo_min,homo_spin,Kmesh,Bst,nq,nsppol,opt_cvk,ucvol,broad,nomega,omega,eps_rpa,dos)
  use m_bz_mesh
  use defs_datatypes
  use defs_basis
  implicit none
  integer,intent(in) :: lomo_min
  integer,intent(in) :: nbnds
  integer,intent(in) :: nomega
  integer,intent(in) :: nq
  integer,intent(in) :: nsppol
  type(ebands_t),intent(in) :: BSt
  type(kmesh_t),intent(in) :: Kmesh
  real(dp),intent(in) :: broad
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: dos(nomega)
  complex(dpc),intent(out) :: eps_rpa(nomega,nq)
  integer,intent(in) :: homo_spin(nsppol)
  integer,intent(in) :: lomo_spin(nsppol)
  complex(dpc),intent(in) :: omega(nomega)
  complex(dpc),intent(in) :: opt_cvk(lomo_min:nbnds,lomo_min:nbnds,Kmesh%nbz,nsppol,nq)
 end subroutine exc_eps_rpa
end interface

interface
 subroutine exc_eps_resonant(Bsp,BS_files,lomo_min,max_band,nkbz,nsppol,opt_cvk,ucvol,nomega,omega,eps_exc,dos_exc)
  use defs_basis
  use m_bs_defs
  implicit none
  integer,intent(in) :: lomo_min
  integer,intent(in) :: max_band
  integer,intent(in) :: nkbz
  integer,intent(in) :: nomega
  integer,intent(in) :: nsppol
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: dos_exc(nomega)
  complex(dpc),intent(out) :: eps_exc(nomega,BSp%nq)
  complex(dpc),intent(in) :: omega(nomega)
  complex(dpc),intent(in) :: opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,nsppol,BSp%nq)
 end subroutine exc_eps_resonant
end interface

interface
 subroutine exc_eps_coupling(Bsp,BS_files,lomo_min,max_band,nkbz,nsppol,opt_cvk,ucvol,nomega,omega,eps_exc,dos_exc)
  use defs_basis
  use m_bs_defs
  implicit none
  integer,intent(in) :: lomo_min
  integer,intent(in) :: max_band
  integer,intent(in) :: nkbz
  integer,intent(in) :: nomega
  integer,intent(in) :: nsppol
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: dos_exc(nomega)
  complex(dpc),intent(out) :: eps_exc(nomega,BSp%nq)
  complex(dpc),intent(in) :: omega(nomega)
  complex(dpc),intent(in) :: opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,nsppol,BSp%nq)
 end subroutine exc_eps_coupling
end interface

interface
 subroutine exc_write_tensor(BSp,BS_files,what,tensor)
  use defs_basis
  use m_bs_defs
  implicit none
  type(excfiles),intent(in) :: BS_files
  type(excparam),intent(in) :: BSp
  character(len=*),intent(in) :: what
  complex(dpc),intent(in) :: tensor(BSp%nomega,6)
 end subroutine exc_write_tensor
end interface

interface
 subroutine mdfs_ncwrite(ncid,Bsp,eps_exc,eps_rpanlf,eps_gwnlf)
  use defs_basis
  use m_bs_defs
  implicit none
  integer,intent(in) :: ncid
  type(excparam),intent(in) :: BSp
  complex(dpc),target,intent(in) :: eps_exc(BSp%nomega,BSp%nq)
  complex(dpc),target,intent(in) :: eps_gwnlf(BSp%nomega,BSp%nq)
  complex(dpc),target,intent(in) :: eps_rpanlf(BSp%nomega,BSp%nq)
 end subroutine mdfs_ncwrite
end interface

interface
 subroutine setup_bse(codvsn,acell,rprim,ngfftf,ngfft_osc,Dtset,Dtfil,BS_files,Psps,Pawtab,BSp,&  
  &  Cryst,Kmesh,Qmesh,KS_BSt,QP_bst,Hdr_wfk,Gsph_x,Gsph_c,Vcp,Hdr_bse,w_fname,comm,Wvl)
  use m_vcoul
  use m_pawtab
  use m_bz_mesh
  use defs_abitypes
  use m_bs_defs
  use m_gsphere
  use m_crystal
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: comm
  type(excfiles),intent(out) :: BS_files
  type(excparam),intent(inout) :: Bsp
  type(crystal_t),intent(out) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(gsphere_t),intent(out) :: Gsph_c
  type(gsphere_t),intent(out) :: Gsph_x
  type(hdr_type),intent(out) :: Hdr_bse
  type(hdr_type),intent(out) :: Hdr_wfk
  type(ebands_t),intent(out) :: KS_BSt
  type(kmesh_t),intent(out) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(ebands_t),intent(out) :: QP_Bst
  type(kmesh_t),intent(out) :: Qmesh
  type(vcoul_t),intent(out) :: Vcp
  type(wvl_internal_type), intent(in) :: Wvl
  character(len=6),intent(in) :: codvsn
  character(len=fnlen),intent(out) :: w_fname
  integer,intent(out) :: ngfft_osc(18)
  integer,intent(in) :: ngfftf(18)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine setup_bse
end interface

interface
 subroutine setup_bse_interp(Dtset,Dtfil,BSp,Cryst,Kmesh,&  
  &  Kmesh_dense,Qmesh_dense,KS_BSt_dense,QP_bst_dense,Gsph_x,Gsph_c,Vcp_dense,Hdr_wfk_dense,ngfftf,grid,comm)
  use m_vcoul
  use m_bz_mesh
  use m_double_grid
  use defs_abitypes
  use m_bs_defs
  use m_crystal
  use m_gsphere
  use defs_datatypes
  implicit none
  integer,intent(in) :: comm
  type(excparam),intent(inout) :: Bsp
  type(crystal_t),intent(in) :: Cryst
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(gsphere_t),intent(out) :: Gsph_c
  type(gsphere_t),intent(out) :: Gsph_x
  type(hdr_type),intent(out) :: Hdr_wfk_dense
  type(ebands_t),intent(out) :: KS_BSt_dense
  type(kmesh_t),intent(in) :: Kmesh
  type(kmesh_t),intent(out) :: Kmesh_dense
  type(ebands_t),intent(out) :: QP_Bst_dense
  type(kmesh_t),intent(out) :: Qmesh_dense
  type(vcoul_t),intent(out) :: Vcp_dense
  type(double_grid_t),intent(out) :: grid
  integer,intent(in) :: ngfftf(18)
 end subroutine setup_bse_interp
end interface

interface
 subroutine check_kramerskronig(n,o,eps)
  use defs_basis
  implicit none
  integer,intent(in) :: n
  complex(dpc),intent(in) :: eps(n)
  real(dp),intent(in) :: o(n)
 end subroutine check_kramerskronig
end interface

interface
 subroutine check_fsumrule(n,o,e2,omegaplasma)
  use defs_basis
  implicit none
  integer,intent(in) :: n
  real(dp),intent(in) :: omegaplasma
  real(dp),intent(in) :: e2(n)
  real(dp),intent(in) :: o(n)
 end subroutine check_fsumrule
end interface

end module interfaces_71_bse
!!***
