!!****m* ABINIT/interfaces_66_fock
!! NAME
!! interfaces_66_fock
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/66_fock
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

module interfaces_66_fock

 implicit none

interface
 subroutine fock_getghc(cwavef,fock,gbound,ghc,gmet,istwf_k,kpoint_i,kg_k,&  
  &  mgfft,mpi_enreg,n4,n5,n6,nfft,ngfft,npw,paral_kgb,use_gpu_cuda)
  use defs_basis
  use defs_abitypes
  use m_fock
  implicit none
  integer, intent(in) :: istwf_k
  integer, intent(in) :: mgfft
  integer, intent(in) :: n4
  integer, intent(in) :: n5
  integer, intent(in) :: n6
  integer, intent(in) :: nfft
  integer, intent(in) :: npw
  integer, intent(in) :: paral_kgb
  integer, intent(in),optional :: use_gpu_cuda
  type(fock_type),pointer,intent(inout) :: fock
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: cwavef(2,npw)
  integer,intent(in) :: gbound(2*mgfft+8,2)
  real(dp),intent(inout) :: ghc(2,npw)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kpoint_i(3)
 end subroutine fock_getghc
end interface

interface
 subroutine fock_updatecwaveocc(cg,dtset,fock,fock_energy,istep,mcg,mpi_enreg,npwarr,occ)
  use defs_basis
  use defs_abitypes
  use m_fock
  implicit none
  integer, intent(in) :: istep
  integer, intent(in) :: mcg
  type(dataset_type),intent(in) :: dtset
  type(fock_type),intent(inout),pointer :: fock
  real(dp) :: fock_energy
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 end subroutine fock_updatecwaveocc
end interface

end module interfaces_66_fock
!!***
