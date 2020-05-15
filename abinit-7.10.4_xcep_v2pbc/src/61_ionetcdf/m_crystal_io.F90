!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_crystal_io
!! NAME
!! m_crystal_io
!!
!! FUNCTION
!! Module containing the methods used to do IO on Crystal objects. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2014 ABINIT group (MG, YP, DC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_crystal_io

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_crystal
 use m_atomdata
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 use defs_abitypes,    only : Hdr_type

 implicit none

 private 
!!***

 public :: crystal_from_hdr   ! Initialize the object from the abinit header.
 public :: crystal_ncwrite    ! Dump the object in a netcdf file.

 ! TODO corresponding reading method is missing, at present use InitCrystalFromHdr

CONTAINS

!!****f* m_crystal_io/crystal_from_hdr 
!! NAME
!!  crystal_from_hdr
!!
!! FUNCTION
!!  Initializes a crystal_t data type starting from the abinit header.
!!
!! INPUTS
!!  Hdr<Hdr_type>=the abinit header
!!  timrev ==2 => take advantage of time-reversal symmetry
!!         ==1 ==> do not use time-reversal symmetry 
!!  remove_inv [optional]= if .TRUE. the inversion symmetry is removed from the set of operations
!!  even though it is present in the header
!!
!! OUTPUT
!!  Cryst<crystal_t>= the data type filled with data reported in the abinit header 
!!
!! TODO
!!  Add information on the use of time-reversal in the Abinit header.
!!
!! PARENTS
!!      m_wfs,mlwfovlp_qp,mrgscr,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!      atomdata_from_znucl,etsf_io_geometry_def,etsf_io_geometry_put
!!      etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine crystal_from_hdr(Cryst,Hdr,timrev,remove_inv)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_from_hdr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(Hdr_type),intent(in) :: Hdr
 type(crystal_t),intent(out) :: Cryst 
 integer,intent(in) :: timrev
 logical,optional,intent(in) :: remove_inv

!Local variables-------------------------------
 integer :: space_group
 logical :: rinv,ltest,use_antiferro
! *********************************************************************

 rinv=.FALSE.; if (PRESENT(remove_inv)) rinv=remove_inv
 use_antiferro=(Hdr%nspden==2.and.Hdr%nsppol==1)

 ! Consistency check
 ltest = (timrev==1.or.timrev==2)
 ABI_CHECK(ltest,"Wrong value for timrev (1|2)")
 if (use_antiferro) then
   ABI_CHECK(ANY(Hdr%symafm==-1),"Wrong nspden, nsppol, symafm.")
 end if

 space_group=0 !FIXME not known

 call crystal_init(Cryst,space_group,Hdr%natom,Hdr%npsp,Hdr%ntypat,Hdr%nsym,Hdr%rprimd,Hdr%typat,Hdr%xred,&
& Hdr%zionpsp,Hdr%znuclpsp,timrev,use_antiferro,rinv,Hdr%title,&
& symrel=Hdr%symrel,tnons=Hdr%tnons,symafm=Hdr%symafm) ! Optional

end subroutine crystal_from_hdr
!!***

!----------------------------------------------------------------------

!!****f* m_crystal_io/crystal_ncwrite
!! NAME
!! crystal_ncwrite
!!
!! FUNCTION
!! Output system geometry to a file, using the NETCDF file format and ETSF I/O.
!! Data are taken from the crystal_t object.
!!
!! INPUTS
!!  Cryst<crystal_t>=Object defining the unit cell and its symmetries.
!!  ncid=NC file handle.
!!
!! OUTPUT
!!  Only writing
!!
!! NOTES
!!  Alchemy not treated, since Crystal should be initialized at the beginning of the run.
!!
!! PARENTS
!!      anaddb,eig2tot,exc_spectra,loper3,m_haydock,m_phonons,m_shirley
!!      outscfcv,sigma
!!
!! CHILDREN
!!      atomdata_from_znucl,etsf_io_geometry_def,etsf_io_geometry_put
!!      etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine crystal_ncwrite(Cryst,ncid)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: ncid
 type(crystal_t),target,intent(in) :: Cryst
!arrays

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 integer :: itypat
 integer,target :: space_group
 logical :: lstat
 !character(len=500) :: msg
 character(len=etsf_charlen) :: symmorphic
 type(atomdata_t) :: atom
 !type(ETSF_groups) :: Group_folder
 type(ETSF_geometry),target :: Geo_folder
 type(ETSF_io_low_error) :: Error_data
#endif
!arrays
 character(len=2),allocatable,target :: symbols(:)
 character(len=80),allocatable,target :: psp_desc(:),symbols_long(:)

! *************************************************************************

#ifdef HAVE_TRIO_ETSF_IO
 !@crystal_t

! call etsf_io_low_open_create(ncid, "hello2_crystal", etsf_file_format_version, lstat,&
!&  Error_data = Error_data, with_etsf_header=.TRUE., overwrite=.FALSE.)
! ETSF_CHECK_ERROR(lstat,Error_data)

 ! === Define subset of ETSF-Dimensions and write them on file ===
 ! * Use low level procedures to write dimensions.
 ! * If dimensions already exist, check that definitions are coherent.

 !call ncfile_write_basedims(ncid)

 call etsf_io_low_set_define_mode(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Basic dimensions that should be always written
! TODO should write an helper function.
 call etsf_io_low_write_dim(ncid,'complex',2,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,"symbol_length",2,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'character_string_length',80,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,"number_of_cartesian_directions",3,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,"number_of_reduced_dimensions",3,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,"number_of_vectors",3,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)
!end obvious stuff

 call etsf_io_low_write_dim(ncid,'number_of_atoms',Cryst%natom,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_atom_species',Cryst%ntypat,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_symmetry_operations',Cryst%nsym,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)
 
 call etsf_io_geometry_def(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! === Fill-in ETSF geometry folder ===
 ! FIXME alchemy not treated since Cryst should be initialized in invars2
 ABI_CHECK(Cryst%npsp==Cryst%ntypat,"alchemy not supported")

 ! Set-up atomic symbols.
 ABI_MALLOC(symbols,(Cryst%ntypat))
 ABI_MALLOC(symbols_long,(Cryst%ntypat))
 ABI_MALLOC(psp_desc,(Cryst%ntypat))

 do itypat=1,Cryst%ntypat  
   call atomdata_from_znucl(atom,Cryst%znucl(itypat)) 
   symbols(itypat) = atom%symbol
   write(symbols_long(itypat),'(a2,a78)') symbols(itypat),REPEAT(CHAR(0),78)
   write(psp_desc(itypat),'(2a)') &
&    Cryst%title(itypat)(1:MIN(80,LEN_TRIM(Cryst%title(itypat)))),REPEAT(CHAR(0),MAX(0,80-LEN_TRIM(Cryst%title(itypat))))
 end do

 space_group=0; if (Cryst%space_group>0) space_group=Cryst%space_group
 Geo_folder%space_group                   => space_group
 Geo_folder%primitive_vectors             => Cryst%rprimd
 Geo_folder%reduced_symmetry_matrices     => Cryst%symrel
 Geo_folder%reduced_symmetry_translations => Cryst%tnons
 Geo_folder%atom_species                  => Cryst%typat
 Geo_folder%reduced_atom_positions        => Cryst%xred
 if (Cryst%npsp==Cryst%ntypat) then
  Geo_folder%valence_charges              => Cryst%zion
 end if
 Geo_folder%atomic_numbers                => Cryst%znucl
 Geo_folder%atom_species_names            => symbols_long
 Geo_folder%chemical_symbols              => symbols
 Geo_folder%pseudopotential_types         => psp_desc

 symmorphic='no'; if (isymmorphic(Cryst)) symmorphic='yes'

 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_geometry_put(ncid, Geo_folder, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ABI_FREE(symbols)
 ABI_FREE(symbols_long)
 ABI_FREE(psp_desc)
 !
 ! === At this point we have an ETSF-compliant file ===
 ! * Add additional stuff for internal use in abinit.

 ! TODO add spinat.

 call etsf_io_low_set_define_mode(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'symafm',etsf_io_low_integer,&
& (/'number_of_symmetry_operations'/),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'symafm',Cryst%symafm,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

#else
 MSG_ERROR('ETSF-IO support is not activated.')
#endif

end subroutine crystal_ncwrite
!!***

end MODULE m_crystal_io
!!***
