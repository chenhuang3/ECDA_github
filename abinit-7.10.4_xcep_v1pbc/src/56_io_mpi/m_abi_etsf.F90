!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abi_etsf
!! NAME
!! m_abi_etsf
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2014 ABINIT group (DCA,YP,MJV,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_abi_etsf

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_atomdata
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 implicit none

 private

#ifdef HAVE_TRIO_ETSF_IO
 public :: etsf_dims_nullify
 public :: abi_etsf_dims_init_low
 public :: abi_etsf_dims_init   
#endif
 public :: abi_etsf_init
 public :: abi_etsf_geo_put
 public :: abi_etsf_electrons_put
 public :: ini_wf_etsf

CONTAINS  !===========================================================
!!***

!!****f* m_abi_etsf/etsf_dims_nullify
!! NAME
!! etsf_dims_nullify
!!
!! FUNCTION
!!  Initialize the structure defining ETSF dimensions.
!!
!! INPUTS
!!
!! OUTPUT
!!  dims=structure with ETSF dimensions.
!!
!! PARENTS
!!      m_abi_etsf
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_write_dim
!!
!! SOURCE

#ifdef HAVE_TRIO_ETSF_IO

subroutine etsf_dims_nullify(Dims,dimlen)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'etsf_dims_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: dimlen
 type(etsf_dims),intent(inout) :: Dims

!Local variables-------------------------------
 integer :: my_dimlen

! *************************************************************************

  my_dimlen = 3; if (PRESENT(dimlen)) my_dimlen = dimlen

  Dims%max_number_of_angular_momenta                 = etsf_no_dimension
  Dims%max_number_of_basis_grid_points               = etsf_no_dimension 
  Dims%max_number_of_coefficients                    = etsf_no_dimension
  Dims%max_number_of_projectors                      = etsf_no_dimension
  Dims%max_number_of_states                          = etsf_no_dimension
  Dims%number_of_atoms                               = etsf_no_dimension
  Dims%number_of_atom_species                        = etsf_no_dimension
  Dims%number_of_cartesian_directions                = my_dimlen
  Dims%number_of_coefficients_dielectric_function    = etsf_no_dimension
  Dims%number_of_components                          = etsf_no_dimension
  Dims%number_of_frequencies_dielectric_function     = etsf_no_dimension
  Dims%number_of_grid_points_vector1                 = etsf_no_dimension
  Dims%number_of_grid_points_vector2                 = etsf_no_dimension
  Dims%number_of_grid_points_vector3                 = etsf_no_dimension
  Dims%number_of_kpoints                             = etsf_no_dimension 
  Dims%number_of_localization_regions                = etsf_no_dimension
  Dims%number_of_qpoints_dielectric_function         = etsf_no_dimension
  Dims%number_of_qpoints_gamma_limit                 = etsf_no_dimension
  Dims%number_of_reduced_dimensions                  = my_dimlen
  Dims%number_of_spinor_components                   = etsf_no_dimension
  Dims%number_of_spins                               = etsf_no_dimension
  Dims%number_of_symmetry_operations                 = etsf_no_dimension
  Dims%number_of_vectors                             = my_dimlen
  Dims%real_or_complex_coefficients                  = etsf_no_dimension
  Dims%real_or_complex_density                       = etsf_no_dimension
  Dims%real_or_complex_gw_corrections                = etsf_no_dimension
  Dims%real_or_complex_potential                     = etsf_no_dimension
  Dims%real_or_complex_wavefunctions                 = etsf_no_dimension
  Dims%symbol_length                                 = etsf_chemlen

  !Dimensions for variables that can be splitted.
  Dims%my_max_number_of_coefficients  = etsf_no_dimension
  Dims%my_max_number_of_states        = etsf_no_dimension
  Dims%my_number_of_components        = etsf_no_dimension
  Dims%my_number_of_grid_points_vect1 = etsf_no_dimension
  Dims%my_number_of_grid_points_vect2 = etsf_no_dimension
  Dims%my_number_of_grid_points_vect3 = etsf_no_dimension
  Dims%my_number_of_kpoints           = etsf_no_dimension
  Dims%my_number_of_spins             = etsf_no_dimension

end subroutine etsf_dims_nullify
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_abi_etsf/abi_etsf_dims_init_low
!! NAME
!! abi_etsf_dims_init_low
!!
!! FUNCTION
!!  Initialize the structure defining ETSF dimensions.
!!  starting from values stored in the dataset_type, the pseudopotential_type.
!!  and the wave function handler for the BIGDFT part.
!!
!! INPUTS
!!
!! OUTPUT
!!  Dims=structure with ETSF dimensions.
!!
!! PARENTS
!!      m_header
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_write_dim
!!
!! SOURCE

#ifdef HAVE_TRIO_ETSF_IO

subroutine abi_etsf_dims_init_low(Dims,natom,ntypat,mband,nkpt,nsppol,nspinor,nspden,nsym,&
&  mpsang,mproj,initialize) ! Optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_etsf_dims_init_low'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,natom,ntypat,nspden,nkpt,nspinor,nsppol,nsym
 integer,optional,intent(in) :: mpsang,mproj
 logical,optional,intent(in) :: initialize
 type(etsf_dims),intent(inout) :: Dims

!Local variables-------------------------------
 logical :: do_init

! *************************************************************************

! Init dimensions with default values.
 do_init = .TRUE.; if (PRESENT(initialize)) do_init = initialize

 if (do_init) then
   call etsf_dims_nullify(Dims)
 end if

!Set-up the dimensions
!=====================
 if (PRESENT(mpsang)) dims%max_number_of_angular_momenta = mpsang
 if (PRESENT(mproj))  dims%max_number_of_projectors      = mproj

 dims%max_number_of_states           = mband
 dims%number_of_atoms                = natom
 dims%number_of_atom_species         = ntypat
 dims%number_of_components           = nspden
 dims%number_of_kpoints              = nkpt
 dims%number_of_spinor_components    = nspinor
 dims%number_of_spins                = nsppol
 dims%number_of_symmetry_operations  = nsym

!!! !In the case of BigDFT, the number of coefficients are the number of wavelets.
!!!  if (dtset%usewvl==0) then
!!!    dims%max_number_of_coefficients      = dtset%mpw
!!!    dims%max_number_of_basis_grid_points = etsf_no_dimension
!!!  else
!!! #ifdef HAVE_DFT_BIGDFT
!!!    dims%max_number_of_coefficients      = wfs%ks%lzd%Glr%wfd%nvctr_c + 7 * wfs%ks%lzd%Glr%wfd%nvctr_f
!!!    dims%max_number_of_basis_grid_points = wfs%ks%lzd%Glr%wfd%nvctr_c
!!! #else
!!!    MSG_ERROR("BIGDFT support is missing")
!!! #endif
!!!  end if
!!! 
!!!  if (dtset%usepaw==1) then
!!!    dims%number_of_grid_points_vector1  = dtset%ngfftdg(1)
!!!    dims%number_of_grid_points_vector2  = dtset%ngfftdg(2)
!!!    dims%number_of_grid_points_vector3  = dtset%ngfftdg(3)
!!!  else if (dtset%usewvl==1) then
!!! #ifdef HAVE_DFT_BIGDFT
!!! !In the case of BigDFT, the grid size is not defined by ngfft.
!!!    dims%number_of_grid_points_vector1  = wfs%ks%lzd%Glr%d%n1 * 2
!!!    dims%number_of_grid_points_vector2  = wfs%ks%lzd%Glr%d%n2 * 2
!!!    dims%number_of_grid_points_vector3  = wfs%ks%lzd%Glr%d%n3 * 2
!!! #endif
!!!  else
!!!    dims%number_of_grid_points_vector1  = dtset%ngfft(1)
!!!    dims%number_of_grid_points_vector2  = dtset%ngfft(2)
!!!    dims%number_of_grid_points_vector3  = dtset%ngfft(3)
!!!  end if

end subroutine abi_etsf_dims_init_low
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_abi_etsf/abi_etsf_dims_init
!! NAME
!! abi_etsf_dims_init
!!
!! FUNCTION
!!  Initialize the structure defining ETSF dimensions.
!!  starting from values stored in the dataset_type, the pseudopotential_type.
!!  and the wave function handler for the BIGDFT part.
!!
!! INPUTS
!!  dtset<type(dataset_type)>=all input variables for this dataset
!!  psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!  wfs<wvl_wf_type>=Object to handle wave functions for the BIGDFT part
!!    Presently not used, likely will be needed when ETSF-IO will be generalized to deal with PAW.
!!  itype = an integer to define what to put in the output file. This can
!!          be one of the following values (maybe a sum latter):
!!          1 for a density file,
!!          2 for a wavefunction file,
!!          4 for a KSS file,
!!          8 for the exchange potential,
!!         16 for the correlation potential.
!!
!! OUTPUT
!!  dims=structure with ETSF dimensions.
!!
!! PARENTS
!!      m_abi_etsf,pawmkaewf
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_write_dim
!!
!! SOURCE

#ifdef HAVE_TRIO_ETSF_IO

subroutine abi_etsf_dims_init(dims, dtset, itype, psps, wfs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_etsf_dims_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itype
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),intent(in) :: wfs
 type(etsf_dims),intent(inout) :: dims

! *************************************************************************

!Set-up the dimensions
!=====================
 dims%max_number_of_angular_momenta = psps%mpsang
 dims%max_number_of_projectors      = 1

 dims%max_number_of_states   = dtset%mband
 dims%number_of_atoms        = dtset%natom
 dims%number_of_atom_species = dtset%ntypat
 dims%number_of_components   = dtset%nspden
                                                     
 dims%number_of_kpoints              = dtset%nkpt
 dims%number_of_spinor_components    = dtset%nspinor
 dims%number_of_spins                = dtset%nsppol
 dims%number_of_symmetry_operations  = dtset%nsym

!In the case of BigDFT, the number of coefficients are the number of wavelets.
 if (dtset%usewvl==0) then
   dims%max_number_of_coefficients      = dtset%mpw
   dims%max_number_of_basis_grid_points = etsf_no_dimension
 else
#ifdef HAVE_DFT_BIGDFT
   dims%max_number_of_coefficients      = wfs%ks%lzd%Glr%wfd%nvctr_c + 7 * wfs%ks%lzd%Glr%wfd%nvctr_f
   dims%max_number_of_basis_grid_points = wfs%ks%lzd%Glr%wfd%nvctr_c
#else
   MSG_ERROR("BIGDFT support is missing")
#endif
 end if

 if (dtset%usepaw==1) then
   dims%number_of_grid_points_vector1  = dtset%ngfftdg(1)
   dims%number_of_grid_points_vector2  = dtset%ngfftdg(2)
   dims%number_of_grid_points_vector3  = dtset%ngfftdg(3)
 else if (dtset%usewvl==1) then
#ifdef HAVE_DFT_BIGDFT
!In the case of BigDFT, the grid size is not defined by ngfft.
   dims%number_of_grid_points_vector1  = wfs%ks%lzd%Glr%d%n1 * 2
   dims%number_of_grid_points_vector2  = wfs%ks%lzd%Glr%d%n2 * 2
   dims%number_of_grid_points_vector3  = wfs%ks%lzd%Glr%d%n3 * 2
#endif
 else
   dims%number_of_grid_points_vector1  = dtset%ngfft(1)
   dims%number_of_grid_points_vector2  = dtset%ngfft(2)
   dims%number_of_grid_points_vector3  = dtset%ngfft(3)
 end if

!The density real_or_complex.
 dims%real_or_complex_density = etsf_no_dimension
 if (iand(itype, 1) /= 0) dims%real_or_complex_density = 1

!The coefficient of wavefunctions real_or_complex.
 dims%real_or_complex_coefficients   = etsf_no_dimension
 if (iand(itype, 2) /= 0 .or. iand(itype, 4) /= 0) then
   if (dtset%usewvl == 0) then
     dims%real_or_complex_coefficients = 2 ! used in plane waves
   else
     dims%real_or_complex_coefficients = 1 ! used in wavelets
   end if
 end if

!The gw corrections real_or_complex.
!Todo: Currently not exported.
!if (.false. .and. iand(itype, 4) /= 0) then
 dims%real_or_complex_gw_corrections = etsf_no_dimension
 if (iand(itype, 4) /= 0) then
   dims%real_or_complex_gw_corrections = 2 ! used in plane waves
 ! dims%real_or_complex_gw_corrections = 1 ! used in plane waves
 end if

!The potential real_or_complex.
 dims%real_or_complex_potential = etsf_no_dimension
 if (iand(itype, 8) /= 0 .or. iand(itype, 16) /= 0) dims%real_or_complex_potential = 1

 dims%real_or_complex_wavefunctions = etsf_no_dimension

end subroutine abi_etsf_dims_init
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_abi_etsf/abi_etsf_init
!! NAME
!! abi_etsf_init
!!
!! FUNCTION
!!  Create a NetCDF file following the ETSF file format specifications.
!!  It declares the dimensions and set-up the different variables, according
!!  to the itype argument.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  filapp = character string giving the root to form the name of the GEO file
!!  itype = an integer to define what to put in the output file. This can
!!          be one of the following values (maybe a sum latter):
!!          1 for a density file,
!!          2 for a wavefunction file,
!!          4 for a KSS file,
!!          8 for the exchange potential,
!!         16 for the correlation potential.
!!  kdep= .true. if the data for the array sizes are dependant on k points.
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!  Data written in file whose name is filapp//'-etsf.nc'
!!
!! PARENTS
!!      m_io_kss,scfcv
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_write_dim
!!
!! SOURCE

subroutine abi_etsf_init(dtset,filapp,itype,kdep,lmn_size,psps,wfs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_etsf_init'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itype
 logical,intent(in) :: kdep
 character(len=fnlen),intent(in) :: filapp
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),intent(in) :: wfs
!arrays
 integer,intent(in) :: lmn_size(psps%npsp)

!Local variables-------------------------------
#ifdef HAVE_TRIO_ETSF_IO
!scalars
 integer :: ncid,var_main,usewvl
 logical :: lstat
 character(len=80) :: file_title
 character(len=fnlen) :: filetsf
 type(etsf_dims) :: dims
 type(etsf_groups_flags) :: flags
 type(etsf_io_low_error) :: error

! *************************************************************************

!Initialize the filename
 filetsf = TRIM(filapp)//'-etsf.nc'
 call wrtout(std_out,ABI_FUNC//': about to create file '//TRIM(filetsf),'COLL')

 usewvl = dtset%usewvl

!Set-up the dimensions
!=====================
 call abi_etsf_dims_init(dims,dtset,itype,psps,wfs)

!Set-up the variables
!====================
!These mandatory values are always written by the hdr_io_etsf() routine.
 flags%geometry  = etsf_geometry_all
 flags%kpoints   = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights
 flags%electrons = etsf_electrons_all - etsf_electrons_x_functional - etsf_electrons_c_functional
 flags%basisdata = etsf_basisdata_basis_set

 if (usewvl==0) then
   flags%basisdata = flags%basisdata + etsf_basisdata_kin_cutoff + etsf_basisdata_n_coeff
 end if

!These variables may be written depending on prt<something> input variables.
 if (itype==1) then
   flags%main = etsf_main_density
   file_title = "Density file"
 else if (itype == 2) then
   if (usewvl==0) then
     flags%basisdata = flags%basisdata + etsf_basisdata_red_coord_pw
   else
     flags%basisdata = flags%basisdata + etsf_basisdata_coord_grid + etsf_basisdata_n_coeff_grid
   end if
   flags%main = etsf_main_wfs_coeff
   file_title = "Wavefunctions file"
 else if (itype==4) then
   if (usewvl==0) then
     flags%basisdata = flags%basisdata + etsf_basisdata_red_coord_pw
   else
     flags%basisdata = flags%basisdata + etsf_basisdata_coord_grid + etsf_basisdata_n_coeff_grid
   end if
   flags%main   = etsf_main_wfs_coeff
   flags%gwdata = etsf_gwdata_all
   file_title = "KSS file"
 else if (itype==8) then
   flags%main = etsf_main_pot_x_only
   file_title = "Exchange potential file"
 else if (itype==16) then
   flags%main = etsf_main_pot_c_only
   file_title = "Correlation potential file"
 else if (itype==24) then
   flags%main = etsf_main_pot_xc
   file_title = "Exchange-correlation potential file"
 end if

!Actually create the file
!========================
!If the group contains main, we remove it for a while to be sure to
!add it at the end, after ABINIT private variables.
 var_main   = flags%main
 flags%main = etsf_main_none

 call etsf_io_data_init(filetsf, flags, dims, file_title, &
& 'File generated by ABINIT with ETSF_IO', lstat, error, overwrite = .true., k_dependent = kdep)
 ETSF_CHECK_ERROR(lstat, error)

!Aadd the private ABINIT information when required.
 call etsf_io_low_open_modify(ncid, filetsf, lstat, error_data = error)
 ETSF_CHECK_ERROR(lstat, error)

!Add the private data
 call ini_wf_etsf(ncid,usewvl,lmn_size,psps%npsp,psps%ntypat)

! Add the main part as last variables in the ETSF file.
 call etsf_io_main_def(ncid, lstat, error, flags = var_main)
 ETSF_CHECK_ERROR(lstat, error)
 
! Close the file.
 call etsf_io_low_close(ncid, lstat, error_data = error)
 ETSF_CHECK_ERROR(lstat, error)
#endif

end subroutine abi_etsf_init
!!***

!----------------------------------------------------------------------

!!****f* m_abi_etsf/abi_etsf_geo_put
!! NAME
!! abi_etsf_geo_put
!!
!! FUNCTION
!!  Output system geometry to a file, using the ETSF I/O file format.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | natom  = number of atoms in unit cell
!!   | ntypat = number of types of atoms in unit cell.
!!   | typat(natom) = type integer for each atom in cell
!!   | znucl(ntypat)= real(dp), atomic number of atom type
!!  filapp = character string giving the root to form the name of the GEO file
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | title = a description for the pseudo.
!!
!! OUTPUT
!!  Data written in file whose name is filapp//'-etsf.nc'
!!
!! PARENTS
!!      m_io_kss,outscfcv,pawmkaewf,sigma
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_write_dim
!!
!! SOURCE

subroutine abi_etsf_geo_put(dtset,filapp,psps)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_etsf_geo_put'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filapp
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps

!Local variables-------------------------------
#ifdef HAVE_TRIO_ETSF_IO
!scalars
 integer :: i
 integer,target :: spgroup
 logical :: lstat
 character(len=fnlen) :: filgeom
 type(atomdata_t) :: atom
 type(etsf_groups) :: group_folder
 type(etsf_geometry),target :: geo_folder
 type(etsf_io_low_error) :: error
!arrays
 character(len=2),allocatable,target :: symbols(:)
 character(len=80),allocatable,target :: psp_desc(:),symbols_long(:)

! *************************************************************************

!Initialize filename
 filgeom=TRIM(filapp)//'-etsf.nc'

 call wrtout(std_out,ABI_FUNC//': about to open file '//TRIM(filgeom),'COLL')

!Set-up atomic symbols
 ABI_MALLOC(symbols,(dtset%ntypat))
 ABI_MALLOC(symbols_long,(dtset%ntypat))
 ABI_MALLOC(psp_desc,(dtset%ntypat))

 do i=1,dtset%ntypat
   call atomdata_from_znucl(atom,dtset%znucl(i))
   symbols(i) = atom%symbol
   write(symbols_long(i), "(A2,A78)") symbols(i), repeat(char(0), 78)
   write(psp_desc(i),"(A,A)") psps%title(i)(1:min(80, len_trim(psps%title(i)))), &
&   repeat(char(0), max(0, 80 - len_trim(psps%title(i))))
 end do

!Fill-in geometry folder
 if (dtset%spgroup > 0) then
   spgroup = dtset%spgroup
   geo_folder%space_group => spgroup
 end if
!To use rprimd or xred, add it as an argument.
!!! geo_folder%primitive_vectors => rprimd
!!! geo_folder%reduced_symmetry_matrices => dtset%symrel
!!! geo_folder%reduced_symmetry_translations => dtset%tnons
!!! geo_folder%atom_species => dtset%typat
!!! geo_folder%reduced_atom_positions => xred
!!! if (psps%npsp == psps%ntypat) then
!!!   geo_folder%valence_charges => psps%zionpsp
!!! end if
!!! geo_folder%atomic_numbers => dtset%znucl
 geo_folder%atom_species_names    => symbols_long
 geo_folder%chemical_symbols      => symbols
 geo_folder%pseudopotential_types => psp_desc

 group_folder%geometry => geo_folder

 call etsf_io_data_write(filgeom, group_folder, lstat, error)
 ETSF_CHECK_ERROR(lstat, error)

!Free memory
 ABI_FREE(symbols)
 ABI_FREE(symbols_long)
 ABI_FREE(psp_desc)
#endif

end subroutine abi_etsf_geo_put
!!***

!----------------------------------------------------------------------

!!****f* m_abi_etsf/abi_etsf_electrons_put
!! NAME
!! abi_etsf_electrons_put
!!
!! FUNCTION
!!  Output system of electrons to a file, using the ETSF I/O file format.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | natom  = number of atoms in unit cell
!!   | ntypat = number of types of atoms in unit cell.
!!   | typat(natom) = type integer for each atom in cell
!!   | znucl(ntypat)= real(dp), atomic number of atom type
!!  filapp = character string giving the root to form the name of the ETSF file
!!
!! OUTPUT
!!  Data written in file whose name is filapp//'-etsf.nc'
!!
!! PARENTS
!!      outscfcv,pawmkaewf,sigma
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_write_dim
!!
!! SOURCE

subroutine abi_etsf_electrons_put(dtset, filapp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_etsf_electrons_put'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filapp
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
#ifdef HAVE_TRIO_ETSF_IO
!scalars
 integer,target :: nelect
 real(dp),target :: tsmear
 logical :: lstat
 character(len = etsf_charlen), target :: smearing
 character(len=fnlen) :: filgeom
 type(etsf_groups) :: group_folder
 type(etsf_electrons),target :: electrons_folder
 type(etsf_io_low_error) :: error

! *************************************************************************

!Initialize filename
 filgeom=trim(filapp)//'-etsf.nc'
 call wrtout(std_out,ABI_FUNC//': about to open file '//TRIM(filgeom),'COLL')

!Fill-in electrons folder
!FIXME ! define XC and smearing scheme
!MG WARNING, in abinit, unlike ETSF, nelect is real to allow for charging and alchemy!!
 nelect = dtset%nelect
 electrons_folder%number_of_electrons => nelect

 if (dtset%occopt == 3) then
   smearing = "Fermi-Dirac"
 else if (dtset%occopt == 4) then
   smearing = "cold smearing of N. Marzari with minimization of the bump"
 else if (dtset%occopt == 5) then
   smearing = "cold smearing of N. Marzari with monotonic function in the tail"
 else if (dtset%occopt == 6) then
   smearing = "Methfessel and Paxton"
 else if (dtset%occopt == 7) then
   smearing = "gaussian"
 else if (dtset%occopt == 8) then
   smearing = "uniform"
 else
   smearing = "none"
 end if

 electrons_folder%smearing_scheme => smearing
 tsmear = dtset%tsmear
 electrons_folder%smearing_width => tsmear

 group_folder%electrons => electrons_folder

 call etsf_io_data_write(filgeom, group_folder, lstat, error)
 ETSF_CHECK_ERROR(lstat, error)
#endif

end subroutine abi_etsf_electrons_put
!!***

!----------------------------------------------------------------------

!!****f* m_abi_etsf/ini_wf_etsf
!! NAME
!! ini_wf_etsf
!!
!! FUNCTION
!! Do initialization of additional dimensions and variables in wavefunction files in ETSF format.
!!
!! INPUTS
!!  usewvl=1 if wavelets are used, 0 otherwise
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  npsp=number of pseudopotentials.
!!  ntypat=number of types of atoms in cell.
!!  ncid=the unit of the open NetCDF file.
!!
!! SIDE EFFECTS 
!!  New dimensions and variables are added to the initial NetCDF file.
!!
!! PARENTS
!!      m_abi_etsf,m_header,pawmkaewf
!!
!! CHILDREN
!!      etsf_io_low_def_var,etsf_io_low_write_dim
!!
!! SOURCE

subroutine ini_wf_etsf(ncid,usewvl,lmn_size,npsp,ntypat)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ini_wf_etsf'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,usewvl,npsp,ntypat
!arrays
 integer,intent(in) :: lmn_size(npsp)

!Local variables-------------------------------
#ifdef HAVE_TRIO_ETSF_IO
 integer :: rhoijdim1
 logical :: lstat
 type(etsf_io_low_error) :: error_data

! *************************************************************************

!Add the none-ETSF dimensions and variables.
!If dimensions already exist, it will check that definitions are coherent.

!Define dimensions.
 call etsf_io_low_write_dim(ncid, "npsp", npsp, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid, "codvsnlen", 6, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!version 44 add first dimension for rhoij = max(lmn_size)*(max(lmn_size)+1)/2
 rhoijdim1 = maxval(lmn_size)
 rhoijdim1 = rhoijdim1 * (rhoijdim1 + 1) / 2
!impose rhoijdim1 >= 1 : if 0, it defaults to NF90_UNLIMITED
 rhoijdim1 = max(rhoijdim1, 1)
 call etsf_io_low_write_dim(ncid, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid, "psptitlen", 132, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 if (usewvl==1) then ! Add the BigDFT private dimensions.
   call etsf_io_low_write_dim(ncid, "number_of_wavelet_resolutions", 2, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if
 
!If variables already exist, it will check that definitions are coherent.

!Define variables.
 call etsf_io_low_def_var(ncid, "date", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "codvsn", etsf_io_low_character, (/ "codvsnlen" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "ecut_eff", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "ecutsm", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "etot", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)
 
 call etsf_io_low_def_var(ncid, "headform", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "fform", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "intxc", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "istwfk", etsf_io_low_integer, (/"number_of_kpoints"/), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "ixc", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "occopt", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "pertcase", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "residm", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "stmbias", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "tphysel", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "tsmear", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Version 44 add ecutdg and usepaw
 call etsf_io_low_def_var(ncid, "ecutdg", etsf_io_low_double, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "usepaw", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Multi-dimensional variables.
 call etsf_io_low_def_var(ncid, "pspcod", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "pspdat", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "pspso", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "pspxc", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "qptn", etsf_io_low_double, (/ "number_of_reduced_dimensions" /), &
& lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "so_psp", etsf_io_low_integer, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "symafm", etsf_io_low_integer, (/ "number_of_symmetry_operations" /), &
& lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "title", etsf_io_low_character, (/ pad("psptitlen"), pad("npsp") /), &
& lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 !if (npsp /= ntypat) then
   call etsf_io_low_def_var(ncid, "zionpsp", etsf_io_low_double, (/ "npsp" /), lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 !end if

 call etsf_io_low_def_var(ncid, "znuclpsp", etsf_io_low_double, (/ "npsp" /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Version 44 add lmn_size and rhoij
 call etsf_io_low_def_var(ncid, "lmn_size", etsf_io_low_integer, (/ "number_of_atom_species" /), &
& lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid, "rhoij", etsf_io_low_double, &
& (/ pad("rhoijdim1"), pad("number_of_components"), pad("number_of_atoms") /), lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Add BigDFT variables.
 call etsf_io_low_def_var(ncid, "usewvl", etsf_io_low_integer, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

!Add the BigDFT private variables.
 if (usewvl == 1) then
   call etsf_io_low_def_var(ncid, "number_of_wavelets", etsf_io_low_integer, (/ "number_of_wavelet_resolutions" /), &
&   lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end if
#endif

!if ETSF_IO is undefined, do nothing
 ABI_UNUSED(ntypat)

end subroutine ini_wf_etsf
!!***

!----------------------------------------------------------------------

END MODULE m_abi_etsf
!!***
