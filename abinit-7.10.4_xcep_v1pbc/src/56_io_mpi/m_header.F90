!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_header
!! NAME
!! m_header
!!
!! FUNCTION
!! This module contains the definition of the abinit header (TODO) 
!! and methods acting on the data type.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (XG, MB, MT, DC, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

MODULE m_header

 use defs_basis
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_wffile
#ifdef HAVE_MPI2
 use mpi
#endif
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io_low_level
 use etsf_io
 use m_ncfile
 use m_abi_etsf
#endif

 use m_copy,          only : alloc_copy
 use m_io_tools,      only : flush_unit, isncfile, file_exists, open_file
 use defs_wvltypes,   only : wvl_internal_type
 use defs_datatypes,  only : ebands_t, pseudopotential_type
 use defs_abitypes,   only : hdr_type, dataset_type
 use m_pawtab,        only : pawtab_type
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_destroy, pawrhoij_io

 implicit none

 private
!!**

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 public :: m_header_init           ! Initialize the global variable HDR_fforms stored in this module.
 public :: hdr_fform2ftype         ! Return the filetype flag, a string with the name of the file and the patch level from fform.
 public :: hdr_ftype2fform         ! Returns the last fform associated to the filetype ftype.
 public :: hdr_init                ! Initialize the header structured datatype and most of its content from dtset and psps.
 public :: hdr_init_lowlvl         ! Low level initialization method for Hdr (no dtset).
 public :: hdr_free               ! Deallocates the components of the header.
 public :: hdr_copy                ! Deep copy of the Header.
 public :: hdr_get_nelect_byocc    ! Returns the number of electrons calculated from Hdr%occ
 public :: isknown_headform        ! Returns .TRUE. is headform is in HDR_KNOWN_HEADFORMS.
 public :: hdr_mpio_skip           ! Skip the abinit header using MPI-IO routines. Return the offset of the first Fortran record after the header.
 public :: hdr_bsize_frecords      ! Compute the size of the Fortran records from the header and formeig.
 public :: hdr_comm                ! Transmits the header datatype with MPI.
 public :: hdr_check               ! Compare two headers.
 public :: hdr_io_etsf             ! I/O of the hdr_type with ETSF-IO.
 public :: hdr_update              ! Update the header structured datatype.
 public :: hdr_read_from_fname     ! Read the header (requires a string with the file name).
 public :: hdr_write_to_fname      ! Write the header (requires a string with the file name).
 public :: hdr_skip                ! Skip the header.
 public :: hdr_io                  ! IO of the header.
 public :: hdr_echo                ! Echo the header.
 !public :: hdr_read
 public :: hdr_fort_read           ! Reads the header from a logical unit associated to a unformatted file.
 public :: hdr_ncread              ! Reads the header from a Netcdf file.
 !public :: hdr_write
 public :: hdr_fort_write          ! Writes the header and fform to unformatted file
 public :: hdr_ncwrite             ! Writes the header and fform to a Netcdf file.

! Generic interface of the routines hdr_skip
 interface hdr_skip  
   module procedure hdr_skip_int
   module procedure hdr_skip_wfftype
 end interface hdr_skip

! Generic interface of the routines hdr_io
 interface hdr_io  
   module procedure hdr_io_int
   module procedure hdr_io_wfftype
 end interface hdr_io

 integer,public,parameter :: HDR_NHEADFORMS=9
 ! Number of abinit header formats used so far.

 integer,public,parameter :: HDR_KNOWN_HEADFORMS( HDR_NHEADFORMS ) = (/23,34,40,41,42,44,53,56,57/)
 ! The list of headforms used so far.

 integer,public,parameter :: HDR_LATEST_HEADFORM = HDR_KNOWN_HEADFORMS( HDR_NHEADFORMS )
 ! The latest headform to be used for writing.

 ! Each filetype is denoted by the filetype flag HDR_* and a patch level.
 ! To add a new filetype do the following:
 !   1) add a new integer filetype flag
 !   2) increase HDR_NFTYPES
 !   3) Modify m_header_init to add a new entry in the array HDR_fforms.

 integer,public,parameter ::   &
&  HDR_WF_PW      = 1,         &
&  HDR_DENSITY    = 2,         &
&  HDR_POTENTIAL  = 3,         &
&  HDR_WF_WVL     = 4,         &
&  HDR_SCREENING  = 5,         &
&  HDR_BS_HAM     = 6,         &
&  HDR_BS_EIG     = 7,         &
&  HDR_BS_HAYDOCK = 8,         &
&  HDR_NFTYPES    = HDR_BS_HAYDOCK

!----------------------------------------------------------------------

 type,private :: fform_t
   !integer :: ftype
   character(len=500) :: fname
   integer,pointer :: fforms(:)
 end type fform_t

 type(fform_t),private,save,allocatable :: HDR_fforms(:)
 ! HDR_fforms(HDR_NFTYPES)
 ! The internal databases with the correspondence ftype ==> (fform, patch_level)

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_header/m_header_init
!! NAME
!! m_header_init
!!
!! FUNCTION
!!  Initialize the global variable HDR_fforms stored in this module.
!!
!! OUTPUT
!!  ierr=A nonzero values signals an inconsistency in the internal tables.
!!
!! SIDE EFFECTS
!!   HDR_fforms(HDR_NFTYPES): Database completely initialized.
!!
!! PARENTS
!!      dfpt_write_cg,ioprof
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine m_header_init(ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'm_header_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ierr
!Local variables-------------------------------
!scalars
 integer :: ft,nft

! *************************************************************************

 !@fform_t
 nft = HDR_NFTYPES
 ABI_DT_MALLOC(HDR_fforms,(nft))

 ierr = 0
 do ft=1,nft
   select case(ft)

   case (HDR_WF_PW)
     HDR_fforms(ft)%fname = "wf_planewave"
     ABI_MALLOC(HDR_fforms(ft)%fforms,(2))
     HDR_fforms(ft)%fforms = (/1,2/)

   case (HDR_DENSITY)
     HDR_fforms(ft)%fname = "density"
     ABI_MALLOC(HDR_fforms(ft)%fforms,(2))
     HDR_fforms(ft)%fforms = (/51,52/)

   case (HDR_POTENTIAL)
     HDR_fforms(ft)%fname = "potential"
     ABI_MALLOC(HDR_fforms(ft)%fforms,(2))
     HDR_fforms(ft)%fforms = (/101,102/)

   case (HDR_WF_WVL)
     HDR_fforms(ft)%fname = "wf_wavelet"
     ABI_MALLOC(HDR_fforms(ft)%fforms,(1))
     HDR_fforms(ft)%fforms = (/202/)

   ! FIXME Screening part
   case (HDR_SCREENING)
     HDR_fforms(ft)%fname = "screening"
     ABI_MALLOC(HDR_fforms(ft)%fforms,(2))
     HDR_fforms(ft)%fforms = (/1002,1102/)

   case (HDR_BS_HAM)
     HDR_fforms(ft)%fname = "bs_hamiltonian"
     ABI_MALLOC(HDR_fforms(ft)%fforms,(1))
     HDR_fforms(ft)%fforms = (/301/)

   case (HDR_BS_EIG)
     HDR_fforms(ft)%fname = "bs_eigenstates"
     ABI_MALLOC(HDR_fforms(ft)%fforms,(1))
     HDR_fforms(ft)%fforms = (/401/)

   case (HDR_BS_HAYDOCK)
     HDR_fforms(ft)%fname = "bs_haydock"
     ABI_MALLOC(HDR_fforms(ft)%fforms,(1))
     HDR_fforms(ft)%fforms = (/501/)

   case default
     ierr=ierr+1
   end select

 end do

end subroutine m_header_init
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_fform2ftype
!! NAME
!! hdr_fform2ftype
!!
!! FUNCTION
!!  Return the filetype flag, a string with the name of the file and the patch level from the
!!  input fform.
!!
!! INPUTS
!!  fform=The value of fform read from the header.
!!
!! OUTPUT
!!  fname=The name of the file.
!!  ftype=The integer flag giving the filetype.
!!  patch_level=The
!!  ierr=A nonzero value signals that fform is not allowed
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_fform2ftype(fform,fname,ftype,patch_level,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_fform2ftype'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fform
 integer,intent(out) :: ftype,patch_level,ierr
 character(len=500),intent(out):: fname

!Local variables-------------------------------
!scalars
 integer :: ii,patch

! *************************************************************************

 !@fform_t
 ierr=0
 do ii=1,SIZE(HDR_fforms)
   if ( ANY( fform == HDR_fforms(ii)%fforms) ) then
     ierr = ierr+1
     ftype = ii
     fname = HDR_fforms(ii)%fname
     do patch=1,SIZE(HDR_fforms(ii)%fforms)
       if (fform == HDR_fforms(ii)%fforms(patch)) then
         patch_level=patch-1
         EXIT
       end if
     end do
   end if
 end do

 if (ierr==0) then ! Invalid fform, raise the error.
   ierr=1
 else if (ierr>1) then
   MSG_ERROR("Internal database is inconsistent")
 else
  ierr=0
 end if

end subroutine hdr_fform2ftype
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_ftype2fform
!! NAME
!! hdr_ftype2fform
!!
!! FUNCTION
!!  Returns the last fform associated to the filetype ftype.
!!
!! INPUTS
!!  ftype
!!
!! OUTPUT
!!  fform = Negative value signal a wrong value of ftype.
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_barrier,mpi_bcast
!!
!! SOURCE

function hdr_ftype2fform(ftype) result(fform)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_ftype2fform'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ftype
 integer :: fform

!Local variables-------------------------------
!scalars
 integer :: last

! *************************************************************************

 !@fform_t
 fform=0
 if ( ftype<1 .or. ftype>SIZE(HDR_fforms) ) RETURN ! Wrong ftype.
 !
 ! Return the last fform associated to this ftype.
 last = SIZE(HDR_fforms(ftype)%fforms)
 fform = HDR_fforms(ftype)%fforms(last)

end function hdr_ftype2fform
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_init
!! NAME
!! hdr_init
!!
!! FUNCTION
!! This subroutine initializes the header structured datatype
!! and most of its content from dtset and psps, and put default values for
!! evolving variables.
!!
!! INPUTS
!! Bands <type(ebands_t)>=band structure information including Brillouin zone description
!! codvsn=code version
!! dtset <type(dataset_type)>=all input variables for this dataset
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! mpi_comm_atom=--optional-- MPI communicator over atoms
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! pertcase=index of the perturbation, or 0 if GS calculation
!! psps <type(pseudopotential_type)>=all the information about psps
!! my_atomtab(:)=Index of the atoms (in global numbering ) treated by current proc (Optional)
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      gstate,loper3,nonlinear,respfn,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_init(Bands,codvsn,dtset,hdr,pawtab,pertcase,psps,wvl, &
&                   mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: pertcase
 integer,intent(in),optional :: mpi_comm_atom
 character(len=6),intent(in) :: codvsn
 type(ebands_t),intent(in) :: Bands
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr !vz_i
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type),intent(in) :: wvl
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

 !@hdr_type

! More checking would be needed ...
 if (dtset%ntypat/=psps%ntypat) then
   write(msg,'(a,2(i0,2x))')' dtset%ntypat and psps%ntypat differs. They are :',dtset%ntypat,psps%ntypat
   MSG_ERROR(msg)
 end if

 if (dtset%npsp/=psps%npsp) then
   write(msg,'(a,2(i0,2x))')' dtset%npsp and psps%npsp differs. They are :',dtset%npsp,psps%npsp
   MSG_ERROR(msg)
 end if

 if (present(mpi_comm_atom)) then
   if (present(mpi_atmtab)) then
     call hdr_init_lowlvl(hdr,Bands,psps,pawtab,wvl,codvsn,&
&     pertcase,dtset%natom,dtset%nsym,dtset%nspden,dtset%ecut,dtset%pawecutdg,dtset%ecutsm,dtset%dilatmx,&
&     dtset%intxc,dtset%ixc,dtset%stmbias,dtset%usewvl,dtset%pawcpxocc,dtset%ngfft,dtset%ngfftdg,dtset%so_psp,&
&     dtset%qptn, dtset%rprimd_orig(:,:,1),dtset%xred_orig(:,:,1),dtset%symrel,dtset%tnons,dtset%symafm,dtset%typat,&
&     mpi_comm_atom=mpi_comm_atom,mpi_atmtab=mpi_atmtab)
   else
     call hdr_init_lowlvl(hdr,Bands,psps,pawtab,wvl,codvsn,&
&     pertcase,dtset%natom,dtset%nsym,dtset%nspden,dtset%ecut,dtset%pawecutdg,dtset%ecutsm,dtset%dilatmx,&
&     dtset%intxc,dtset%ixc,dtset%stmbias,dtset%usewvl,dtset%pawcpxocc,dtset%ngfft,dtset%ngfftdg,dtset%so_psp,&
&     dtset%qptn, dtset%rprimd_orig(:,:,1),dtset%xred_orig(:,:,1),dtset%symrel,dtset%tnons,dtset%symafm,dtset%typat,&
&     mpi_comm_atom=mpi_comm_atom)
   end if
 else
   call hdr_init_lowlvl(hdr,Bands,psps,pawtab,wvl,codvsn,&
&   pertcase,dtset%natom,dtset%nsym,dtset%nspden,dtset%ecut,dtset%pawecutdg,dtset%ecutsm,dtset%dilatmx,&
&   dtset%intxc,dtset%ixc,dtset%stmbias,dtset%usewvl,dtset%pawcpxocc,dtset%ngfft,dtset%ngfftdg,dtset%so_psp,&
&   dtset%qptn, dtset%rprimd_orig(:,:,1),dtset%xred_orig(:,:,1),dtset%symrel,dtset%tnons,dtset%symafm,dtset%typat)
 end if

end subroutine hdr_init
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_free
!! NAME
!! hdr_free
!!
!! FUNCTION
!! This subroutine deallocates the components of the header structured datatype
!!
!! INPUTS
!! hdr <type(hdr_type)>=the header
!!
!! OUTPUT
!!  (only deallocate)
!!
!! PARENTS
!!      bethe_salpeter,conducti_nc,conducti_paw,conducti_paw_core,cut3d,elphon
!!      emispec_paw,finddistrproc,gstate,gw_tools,initaim,inpgkk,inwffil,ioarr
!!      ioprof,kss2wfk,linear_optics_paw,loper3,m_bse_io,m_ebands,m_gamma
!!      m_io_gkk,m_io_kss,m_io_screening,m_wfk,m_wfs,macroave,mrggkk,nonlinear
!!      optic,read_el_veloc,read_gkk,respfn,screening,setup_screening
!!      setup_sigma,sigma,suscep
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_free(hdr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(hdr_type),intent(inout) :: hdr

! *************************************************************************

 DBG_ENTER("COLL")

 !@hdr_type

 !integer
 if (allocated(hdr%istwfk)) then
   ABI_FREE(hdr%istwfk)
 end if
 if (allocated(hdr%lmn_size)) then
   ABI_FREE(hdr%lmn_size)
 end if
 if (allocated(hdr%nband)) then
   ABI_FREE(hdr%nband)
 end if
 if (allocated(hdr%npwarr)) then
   ABI_FREE(hdr%npwarr)
 end if

 if (allocated(hdr%pspcod)) then
   ABI_FREE(hdr%pspcod)
 end if
 if (allocated(hdr%pspdat)) then
   ABI_FREE(hdr%pspdat)
 end if
 if (allocated(hdr%pspso)) then
   ABI_FREE(hdr%pspso)
 end if
 if (allocated(hdr%pspxc)) then
   ABI_FREE(hdr%pspxc)
 end if
 if (allocated(hdr%so_psp)) then
   ABI_FREE(hdr%so_psp)
 end if
 if (allocated(hdr%symafm)) then
   ABI_FREE(hdr%symafm)
 end if
 if (allocated(hdr%symrel)) then
   ABI_FREE(hdr%symrel)
 end if
 if (allocated(hdr%typat)) then
   ABI_FREE(hdr%typat)
 end if

 !real
 if (allocated(hdr%kptns)) then
   ABI_FREE(hdr%kptns)
 end if
 if (allocated(hdr%occ)) then
   ABI_FREE(hdr%occ)
 end if
 if (allocated(hdr%tnons)) then
   ABI_FREE(hdr%tnons)
 end if
 if (allocated(hdr%wtk)) then
   ABI_FREE(hdr%wtk)
 end if
 if (allocated(hdr%xred)) then
   ABI_FREE(hdr%xred)
 end if
 if (allocated(hdr%zionpsp)) then
   ABI_FREE(hdr%zionpsp)
 end if
 if (allocated(hdr%znuclpsp)) then
   ABI_FREE(hdr%znuclpsp)
 end if
 if (allocated(hdr%znucltypat)) then
   ABI_FREE(hdr%znucltypat)
 end if

 !string arrays
 if(allocated(hdr%title)) then
   ABI_FREE(hdr%title)
 end if

 if (hdr%usepaw==1 .and. allocated(hdr%pawrhoij) ) then
   call pawrhoij_destroy(hdr%pawrhoij)
   ABI_DT_FREE(hdr%pawrhoij)
 end if

 DBG_EXIT("COLL")

end subroutine hdr_free
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_copy
!! NAME
!! hdr_copy
!!
!! FUNCTION
!! Make a deep copy of the abinit header.
!!
!! INPUTS
!!  Hdr_in=The header to be copied.
!!
!! OUTPUT
!!  Hdr_cp=The deep copy of Hdr_in.
!!
!! NOTES
!!  The present version deals with versions of the header up to 56.
!!
!! PARENTS
!!      m_io_kss,m_io_screening,m_wfk
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_copy(Hdr_in,Hdr_cp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_copy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(hdr_type),intent(in) :: Hdr_in
 type(hdr_type),intent(inout) :: Hdr_cp

!Local variables-------------------------------
!scalars
 integer :: cplex
 character(len=500) :: msg
! *************************************************************************

 !@hdr_type
 if (Hdr_in%headform>57) then
   write(msg,'(3a,i3,2a)')&
&   'hdr_copy deals with versions of the header only up to 57. ',ch10,&
&   'However headform = ',Hdr_in%headform,ch10,&
&   'Change the source to add the changes done in the new version. '
   MSG_ERROR(msg)
 end if

!=== Integer values ===
 Hdr_cp%bantot   = Hdr_in%bantot
 Hdr_cp%date     = Hdr_in%date
 Hdr_cp%headform = Hdr_in%headform
 Hdr_cp%intxc    = Hdr_in%intxc
 Hdr_cp%ixc      = Hdr_in%ixc
 Hdr_cp%natom    = Hdr_in%natom
 Hdr_cp%nkpt     = Hdr_in%nkpt
 Hdr_cp%npsp     = Hdr_in%npsp
 Hdr_cp%nspden   = Hdr_in%nspden
 Hdr_cp%nspinor  = Hdr_in%nspinor
 Hdr_cp%nsppol   = Hdr_in%nsppol
 Hdr_cp%nsym     = Hdr_in%nsym
 Hdr_cp%ntypat   = Hdr_in%ntypat
 Hdr_cp%occopt   = Hdr_in%occopt
 Hdr_cp%pertcase = Hdr_in%pertcase
 Hdr_cp%usepaw   = Hdr_in%usepaw
 Hdr_cp%usewvl   = Hdr_in%usewvl

 ! === Integer arrays ===
 Hdr_cp%ngfft(:)   = Hdr_in%ngfft(:)
 Hdr_cp%nwvlarr(:) = Hdr_in%nwvlarr(:)

!=== Integer pointers ====
 call alloc_copy( Hdr_in%istwfk,  Hdr_cp%istwfk   )
 call alloc_copy( Hdr_in%lmn_size,Hdr_cp%lmn_size )
 call alloc_copy( Hdr_in%nband,   Hdr_cp%nband    )
 call alloc_copy( Hdr_in%npwarr,  Hdr_cp%npwarr   )
 call alloc_copy( Hdr_in%pspcod,  Hdr_cp%pspcod )
 call alloc_copy( Hdr_in%pspdat,  Hdr_cp%pspdat )
 call alloc_copy( Hdr_in%pspso ,  Hdr_cp%pspso  )
 call alloc_copy( Hdr_in%pspxc ,  Hdr_cp%pspxc  )
 call alloc_copy( Hdr_in%so_psp,  Hdr_cp%so_psp )
 call alloc_copy( Hdr_in%symafm,  Hdr_cp%symafm )
 call alloc_copy( Hdr_in%symrel,  Hdr_cp%symrel )
 call alloc_copy( Hdr_in%typat ,  Hdr_cp%typat  )

!=== Real variables ====
 Hdr_cp%ecut        = Hdr_in%ecut
 Hdr_cp%ecutdg      = Hdr_in%ecutdg
 Hdr_cp%ecutsm      = Hdr_in%ecutsm
 Hdr_cp%ecut_eff    = Hdr_in%ecut_eff
 Hdr_cp%etot        = Hdr_in%etot
 Hdr_cp%fermie      = Hdr_in%fermie
 Hdr_cp%residm      = Hdr_in%residm
 Hdr_cp%stmbias     = Hdr_in%stmbias
 Hdr_cp%tphysel     = Hdr_in%tphysel
 Hdr_cp%tsmear      = Hdr_in%tsmear

 Hdr_cp%qptn(:)     = Hdr_in%qptn(:)
 Hdr_cp%rprimd(:,:) = Hdr_in%rprimd(:,:)

!=== Real pointers ===
 call alloc_copy( Hdr_in%kptns     ,Hdr_cp%kptns     )
 call alloc_copy( Hdr_in%occ       ,Hdr_cp%occ       )
 call alloc_copy( Hdr_in%tnons     ,Hdr_cp%tnons     )
 call alloc_copy( Hdr_in%wtk       ,Hdr_cp%wtk       )
 call alloc_copy( Hdr_in%xred      ,Hdr_cp%xred      )
 call alloc_copy( Hdr_in%zionpsp   ,Hdr_cp%zionpsp   )
 call alloc_copy( Hdr_in%znuclpsp  ,Hdr_cp%znuclpsp  )
 call alloc_copy( Hdr_in%znucltypat,Hdr_cp%znucltypat)

!=== Character pointers ===
 Hdr_cp%codvsn = Hdr_in%codvsn
! THIS DOES NOT WORK ON XLF: Hdr_cp%title string length becomes huge and segfaults
! call alloc_copy( Hdr_in%title,Hdr_cp%title )
 ABI_MALLOC(Hdr_cp%title,(Hdr_cp%npsp))
 Hdr_cp%title = Hdr_in%title

!=== For PAW have to copy Pawrhoij ====
!* TODO alchemy requires a different treatment but for the moment it is not available within PAW.
 if (Hdr_in%usepaw==1) then
   cplex = Hdr_in%Pawrhoij(1)%cplex
   ABI_DT_MALLOC(Hdr_cp%Pawrhoij,(Hdr_in%natom))
   call pawrhoij_alloc(Hdr_cp%Pawrhoij,cplex,Hdr_in%nspden,Hdr_in%nspinor,Hdr_in%nsppol,Hdr_in%typat,&
&    lmnsize=Hdr_in%lmn_size(1:Hdr_in%ntypat))
   call pawrhoij_copy(Hdr_in%Pawrhoij,Hdr_cp%Pawrhoij)
 end if

end subroutine hdr_copy
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_get_nelect_byocc
!! NAME
!! hdr_get_nelect_byocc
!!
!! FUNCTION
!!  Return the number of electrons from the occupation numbers
!!  thus taking into account a possible additional charge or alchemy.
!!
!! INPUTS
!!  Hdr<hdr_type>
!!
!! OUTPUT
!!  nelect=Number of electrons in the unit cell.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function hdr_get_nelect_byocc(Hdr) result(nelect)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_get_nelect_byocc'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 type(hdr_type),intent(in) :: Hdr
 real(dp) :: nelect

!Local variables ---------------------------------------
!scalars
 integer :: idx,isppol,ikibz,nband_k

! *************************************************************************

!* Cannot use znucl because we might have additional charge or alchemy.
 nelect=zero ; idx=0
 do isppol=1,Hdr%nsppol
   do ikibz=1,Hdr%nkpt
     nband_k=Hdr%nband(ikibz+(isppol-1)*Hdr%nkpt)
     nelect = nelect + Hdr%wtk(ikibz)*SUM(Hdr%occ(idx+1:idx+nband_k))
     idx=idx+nband_k
   end do
 end do

!Might also check also Hdr%znuclpsp(:) to avoid round off errors

end function hdr_get_nelect_byocc
!!***

!----------------------------------------------------------------------

!!****f* m_header/isknown_headform
!! NAME
!!  isknown_headform
!!
!! FUNCTION
!!  Returns .TRUE. if headform is one of the allowed values.
!!
!! INPUTS
!!  headform
!!
!! SOURCE

function isknown_headform(headform) result(ans)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isknown_headform'
!End of the abilint section

 integer,intent(in) :: headform
 logical :: ans

! *************************************************************************

 ans = ANY(headform == HDR_KNOWN_HEADFORMS)

end function isknown_headform
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_init_lowlvl
!! NAME
!! hdr_init_lowlvl
!!
!! FUNCTION
!! This subroutine initializes the header structured datatype
!! and most of its content from psps and other input variables that
!! are passed explicitly. It also use default values for evolving variables.
!! Note that Dtset is not required thus rendering the initialization of the header
!! much easier.
!!
!! INPUTS
!! Bands <type(ebands_t)>=band structure information including Brillouin zone description
!! codvsn=code version
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! mpi_comm_atom=--optional-- MPI communicator over atoms
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! pertcase=index of the perturbation, or 0 if GS calculation
!! psps <type(pseudopotential_type)>=all the information about psps
!! For the meaning of the other varialble see the definition of dataset_type.
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      m_header,m_wfs
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_init_lowlvl(hdr,Bands,psps,pawtab,wvl,&
&  codvsn,pertcase,natom,nsym,nspden,ecut,pawecutdg,ecutsm,dilatmx,&
&  intxc,ixc,stmbias,usewvl,pawcpxocc,ngfft,ngfftdg,so_psp,qptn,&
&  rprimd,xred,symrel,tnons,symafm,typat, &
&  mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_init_lowlvl'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym,nspden,intxc,ixc,usewvl,pawcpxocc,pertcase
 integer, intent(in),optional :: mpi_comm_atom
 real(dp),intent(in) :: ecut,ecutsm,dilatmx,stmbias,pawecutdg
 character(len=6),intent(in) :: codvsn
 type(ebands_t),intent(in) :: Bands
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type),intent(in) :: wvl
 type(hdr_type),intent(inout) :: hdr !vz_i
!arrays
 integer,intent(in) :: typat(natom)
 integer,intent(in) ::  so_psp(psps%npsp)
 integer,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 integer,intent(in) :: ngfft(18),ngfftdg(18)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(in) :: qptn(3) ! the wavevector, in case of a perturbation
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: bantot,date,nkpt,npsp,ntypat,nsppol,nspinor
 integer :: idx,isppol,ikpt,iband,ipsp
#ifndef HAVE_DFT_BIGDFT
 character(len=500) :: msg
#endif
 character(len=8) :: date_time
!arrays

! *************************************************************************

 !@hdr_type
 call date_and_time(date_time)
 read(date_time,'(i8)')date

 npsp   = psps%npsp
 ntypat = psps%ntypat
 nkpt   = Bands%nkpt
 nsppol = Bands%nsppol
 nspinor= Bands%nspinor
 bantot = Bands%bantot

!Transfer dimensions and other scalars to hdr.
 hdr%intxc    =intxc
 hdr%ixc      =ixc
 hdr%natom    =natom
 hdr%npsp     =npsp
 hdr%nspden   =nspden
 hdr%nspinor  =nspinor
 hdr%nsym     =nsym
 hdr%ntypat   =ntypat
 hdr%bantot   =bantot
 hdr%nkpt     =nkpt
 hdr%nsppol   =nsppol
 hdr%usepaw   =psps%usepaw
 hdr%usewvl   =usewvl !hdr%nwvlarr will be set later since the number !of wavelets have not yet been computed.
 hdr%occopt   =Bands%occopt
 hdr%codvsn   =codvsn
 hdr%date     =date
 hdr%headform =HDR_LATEST_HEADFORM ! Initialize with the latest headform
 hdr%pertcase =pertcase
 hdr%ecut     =ecut
 hdr%ecutsm   =ecutsm
 hdr%ecut_eff =ecut * (dilatmx)**2
 hdr%stmbias  =stmbias
 hdr%tphysel  =Bands%tphysel
 hdr%tsmear   =Bands%tsmear
 hdr%qptn     =qptn
 hdr%rprimd   =rprimd      ! Evolving data

!Default for other data  (all evolving data)
 hdr%etot     =1.0d20
 hdr%fermie   =1.0d20
 hdr%residm   =1.0d20

!Allocate all components of hdr

!Transfer data from Bands
 ABI_MALLOC(hdr%istwfk,(nkpt))
 hdr%istwfk(1:nkpt)      =Bands%istwfk(1:nkpt)
 ABI_MALLOC(hdr%kptns,(3,nkpt))
 hdr%kptns(:,:)          =Bands%kptns(:,:)
 ABI_MALLOC(hdr%nband,(nkpt*nsppol))
 hdr%nband(1:nkpt*nsppol)=Bands%nband(1:nkpt*nsppol)
 ABI_MALLOC(hdr%npwarr,(nkpt))
 hdr%npwarr(:)           =Bands%npwarr(:)
 ABI_MALLOC(hdr%wtk,(nkpt))
 hdr%wtk(:)=Bands%wtk(:)

!Transfer data from psps
 ABI_MALLOC(hdr%pspcod,(npsp))
 hdr%pspcod    =psps%pspcod
 ABI_MALLOC(hdr%pspdat,(npsp))
 hdr%pspdat    =psps%pspdat
 ABI_MALLOC(hdr%pspso,(npsp))
 hdr%pspso     =psps%pspso
 ABI_MALLOC(hdr%pspxc,(npsp))
 hdr%pspxc     =psps%pspxc
 ABI_MALLOC(hdr%znuclpsp,(npsp))
 hdr%znuclpsp  =psps%znuclpsp
 ABI_MALLOC(hdr%znucltypat,(ntypat))
 hdr%znucltypat=psps%znucltypat
 ABI_MALLOC(hdr%zionpsp,(npsp))
 hdr%zionpsp   =psps%zionpsp
 ABI_MALLOC(hdr%title,(npsp))
 do ipsp=1,psps%npsp
   write(hdr%title(ipsp), "(A)") psps%title(ipsp)(1:132)
 end do

 ABI_MALLOC(hdr%so_psp,(npsp))
 hdr%so_psp=so_psp

 ABI_MALLOC(hdr%symafm,(nsym))
 hdr%symafm(1:min(size(symafm),size(hdr%symafm)))=symafm(1:min(size(symafm),size(hdr%symafm)))

 ABI_MALLOC(hdr%symrel,(3,3,nsym))
 hdr%symrel(:,:,1:min(size(symrel,3),size(hdr%symrel,3))) =symrel(:,:,1:min(size(symrel,3),size(hdr%symrel,3)))

 ABI_MALLOC(hdr%tnons,(3,nsym))
 hdr%tnons(:,1:min(size(tnons,2),size(hdr%tnons,2)))=tnons(:,1:min(size(tnons,2),size(hdr%tnons,2)))

 ABI_MALLOC(hdr%typat,(natom))
 hdr%typat(1:natom) =typat(1:natom)  ! PMA : in tests/v2/t11 size(dtset%typat) is bigger dtset%natom
 ABI_MALLOC(hdr%xred,(3,natom))
 hdr%xred(:,1:natom)=xred(:,1:natom) ! Evolving data

 if (psps%usepaw==1)then
   ABI_DT_MALLOC(hdr%pawrhoij,(natom))
   !Values of nspden/nspinor/nsppol are dummy ones; they are overwritten later (by hdr_update)
   if (present(mpi_comm_atom)) then
     if (present(mpi_atmtab)) then
       call pawrhoij_alloc(hdr%pawrhoij,pawcpxocc,nspden,nspinor,nsppol,typat, &
&                       pawtab=pawtab,mpi_comm_atom=mpi_comm_atom,mpi_atmtab=mpi_atmtab)
     else
       call pawrhoij_alloc(hdr%pawrhoij,pawcpxocc,nspden,nspinor,nsppol,typat, &
&                       pawtab=pawtab,mpi_comm_atom=mpi_comm_atom)
     end if
   else
     call pawrhoij_alloc(hdr%pawrhoij,pawcpxocc,nspden,nspinor,nsppol,typat,pawtab=pawtab)
   end if
 end if

 if (psps%usepaw==1 .and. usewvl ==0 ) then
   hdr%ngfft(:) =ngfftdg(1:3)
 else if (usewvl==1) then
#if defined HAVE_DFT_BIGDFT
   hdr%ngfft(:) = (/ wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i /)
#else
 BIGDFT_NOTENABLED_ERROR()
#endif
 else
   hdr%ngfft(:) =ngfft(1:3)
 end if

!Transfer paw data
 ABI_MALLOC(hdr%lmn_size,(npsp))
 if(psps%usepaw==1) then
   hdr%ecutdg   =pawecutdg
   hdr%lmn_size(1:npsp)=pawtab(1:npsp)%lmn_size
 else
   hdr%ecutdg=hdr%ecut
   hdr%lmn_size(:)=psps%lmnmax
 end if

 ABI_MALLOC(hdr%occ,(bantot))
 hdr%occ(:)=zero; idx=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband=1,hdr%nband(ikpt+(isppol-1)*nkpt)
       idx=idx+1
       hdr%occ(idx)=Bands%occ(iband,ikpt,isppol)
     end do
   end do
 end do

end subroutine hdr_init_lowlvl
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_read_from_fname
!! NAME
!! hdr_read_from_fname
!!
!! FUNCTION
!!  Read the header from file fname. 
!!  Use Fortran IO or Netcdf depending on the extension of the file
!!  Only rank0 process reads the header and then broadcast data to the other
!!  processes inside comm.
!!
!! INPUTS
!!  fname=String with the name of the file.
!!  comm = MPI communicator.
!!
!! OUTPUT
!!  Hdr<hdr_type>=The abinit header.
!!  fform=Kind of the array in the file (0 signals an error)
!!
!! PARENTS
!!      finddistrproc,ioprof,m_wfk,m_wfs
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_read_from_fname(Hdr,fname,fform,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_read_from_fname'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: comm
 integer,intent(out) :: fform
 character(len=*),intent(in) :: fname
 type(hdr_type),intent(inout) :: Hdr !vz_i

!Local variables-------------------------------
 integer,parameter :: rdwr1=1,master=0
 integer :: fh,my_rank,mpierr
 character(len=500) :: msg
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(etsf_io_low_error) :: Error
#endif

! *************************************************************************

 my_rank = xcomm_rank(comm)

 if (my_rank==master) then
   if (.not.isncfile(fname)) then
     ! Use Fortran IO to open the file and read the header.
     if (open_file(fname,msg,newunit=fh,form="unformatted", status="old") /= 0) then
       MSG_ERROR(msg)
     end if

     call hdr_fort_read(Hdr,fh,fform,rewind=(rdwr1==1))
     close(fh)

   else
     ! Use Netcdf to open the file and read the header.
#ifdef HAVE_TRIO_ETSF_IO
     call etsf_io_low_open_read(fh, fname, lstat, Error_data= Error, with_etsf_header=.FALSE.)
     ETSF_CHECK_ERROR(lstat, Error)

     call hdr_ncread(Hdr,fh, fform)

     call etsf_io_low_close(fh, lstat, Error_data=Error)
     ETSF_CHECK_ERROR(lstat, Error)
#else
     MSG_ERROR("ETSF-IO not enabled")
#endif
   end if
 end if
 !
 ! Broadcast fform and the header.
 if (xcomm_size(comm) > 1) then
   call hdr_comm(Hdr,master,my_rank,comm)
   call xmpi_bcast(fform,master,comm,mpierr)
 end if

end subroutine hdr_read_from_fname
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_write_to_fname
!! NAME
!! hdr_write_to_fname
!!
!! FUNCTION
!!  Write the header and fform to file fname. 
!!  Use Fortran IO or Netcdf depending on the extension of the file
!!
!! INPUTS
!!  fname=String with the name of the file.
!!  fform=Kind of the array in the file
!!  Hdr<hdr_type>=The abinit header.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_write_to_fname(Hdr,fname,fform) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_write_to_fname'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: fform
 character(len=*),intent(in) :: fname
 type(hdr_type),intent(inout) :: Hdr

!Local variables-------------------------------
 integer,parameter :: rdwr2=2
 integer :: fh,ierr
 character(len=500) :: msg
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(etsf_io_low_error) :: Error
#endif

! *************************************************************************

 ! TODO: split the (read|write) cases in the hdr_io_* routines.
 !fform_copy = fform ! So that we can use intent(in)

 if (.not.isncfile(fname)) then
   !
   ! Use Fortran IO to write the header.
   if (open_file(fname,msg,newunit=fh,form="unformatted", status="unknown") /= 0) then
     MSG_ERROR(msg)
   end if
   !call hdr_io(fform_copy,Hdr,rdwr2,fh)
   call hdr_fort_write(Hdr,fh,fform,ierr)
   ABI_CHECK(ierr==0," Error while writing Abinit header.")
   close(fh)

 else
   ! Use Netcdf to open the file and read the header.
#ifdef HAVE_TRIO_ETSF_IO
   if (file_exists(fname)) then
     call etsf_io_low_open_modify(fh,fname,lstat,Error_data=Error,with_etsf_header=.FALSE.)
   else
     call etsf_io_low_open_create(fh,fname,etsf_file_format_version,lstat,Error_data=Error,with_etsf_header=.TRUE.)
   end if
   ETSF_CHECK_ERROR(lstat, Error)

   call hdr_ncwrite(Hdr,fh,fform,nc_define=.True.)

   call etsf_io_low_close(fh, lstat, Error_data=Error)
   ETSF_CHECK_ERROR(lstat, Error)
#else
   MSG_ERROR("ETSF-IO not enabled")
#endif
 end if

end subroutine hdr_write_to_fname
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_mpio_skip
!! NAME
!!  hdr_mio_skip
!!
!! FUNCTION
!!   Skip the abinit header in MPI-IO mode. This routine uses local MPI-IO calls hence
!!   it can be safely called by master node only. Note however that in this case the
!!   offset has to be communicated to the other nodes.
!!
!! INPUTS
!!  mpio_fh=MPI-IO file handler
!!
!! OUTPUT
!!  fform=kind of the array in the file
!!  offset=The offset of the Fortran record located immediately below the Abinit header.
!!
!! SOURCE

subroutine hdr_mpio_skip(mpio_fh,fform,offset)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_mpio_skip'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: mpio_fh
 integer,intent(out) :: fform
 integer(kind=XMPI_OFFSET_KIND),intent(out) :: offset

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,mpi_type_frm
#ifdef HAVE_MPI_IO
 integer :: headform,ierr,mu,usepaw,npsp
!arrays
 integer(kind=MPI_OFFSET_KIND) :: fmarker,positloc
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *************************************************************************

 offset = 0; fform  = 0

 bsize_frm    = xmpio_bsize_frm    ! bsize_frm= Byte length of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

#ifdef HAVE_MPI_IO
!Reading the first record of the file -------------------------------------
!read (unitfi)   codvsn,headform,..............
 positloc = bsize_frm + 6*xmpi_bsize_ch
 call MPI_FILE_READ_AT(mpio_fh,positloc,fform,1,MPI_INTEGER,statux,ierr)

 if (ANY(fform == (/1,2,51,52,101,102/) )) then
   headform=22  ! This is the old format !read (unitfi) codvsn,fform
 else
   !read (unitfi)codvsn,headform,fform
   call MPI_FILE_READ_AT(mpio_fh,positloc,headform,1,MPI_INTEGER,statux,ierr)
   positloc = positloc + xmpi_bsize_int
   call MPI_FILE_READ_AT(mpio_fh,positloc,fform,1,MPI_INTEGER,statux,ierr)
 end if

 ! Skip first record.
 call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)

!Reading the second record of the file: read(unitfi) bantot, hdr%date, hdr%intxc.................
!Read npsp and usepaw.
 positloc  = offset + bsize_frm + xmpi_bsize_int*13
 call MPI_FILE_READ_AT(mpio_fh,positloc,npsp,1,MPI_INTEGER,statux,ierr)

 usepaw=0
 if (headform >= 44) then
   positloc = positloc +  xmpi_bsize_int*4
   call MPI_FILE_READ_AT(mpio_fh,positloc,usepaw,1,MPI_INTEGER,statux,ierr)
 end if

 ! Skip second record.
 call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)

 ! Skip the rest of the file ---------------------------------------------
 do mu=1,2+npsp
   call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)
 end do

 if (headform>=44.and.usepaw==1) then ! skip rhoij records.
   call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)
   call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)
 end if

#else
 MSG_ERROR("hdr_mpio_skip cannot be used when MPI-IO is not enabled")
 ABI_UNUSED(mpio_fh)
#endif

end subroutine hdr_mpio_skip
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_bsize_frecords
!! NAME
!!  hdr_bsize_frecords
!!
!! FUNCTION
!!  Compute the size of the Fortran records of the WFK file from the header and formeig.
!!
!! INPUTS
!!  Hdr<hdr_type>=The abinit header.
!!  formeig = 0 for GS WFK, 1 for response functio WFK.
!!
!! OUTPUTS
!!  nfrec = Number fof Fortran records
!!  bsize_frecords(nfrec) = Byte size of each records 
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_bsize_frecords(Hdr,formeig,nfrec,bsize_frecords)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_bsize_frecords'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig
 integer,intent(out) :: nfrec
 type(hdr_type),intent(in) :: Hdr
!arrays
 integer(XMPI_OFFSET_KIND),pointer :: bsize_frecords(:)

!Local variables-------------------------------
!scalars
 integer :: max_nfrec,ik_ibz,spin,mband,nband_k,npw_k,band
!arrays
 integer(XMPI_OFFSET_KIND),allocatable :: bsz_frec(:)

!************************************************************************

!@hdr_type
 mband = MAXVAL(Hdr%nband)
 max_nfrec = Hdr%nkpt*Hdr%nsppol * (3 + mband)

 if (formeig==1) max_nfrec = max_nfrec + Hdr%nkpt*Hdr%nsppol*mband
 ABI_MALLOC(bsz_frec, (max_nfrec))

 nfrec = 0
 do spin=1,Hdr%nsppol
   do ik_ibz=1,Hdr%nkpt
     nband_k = Hdr%nband(ik_ibz + (spin-1)*Hdr%nkpt)
     npw_k   = Hdr%npwarr(ik_ibz)

     ! First record: npw, nspinor, nband_disk
     nfrec = nfrec + 1
     bsz_frec(nfrec) = 3*xmpi_bsize_int
     !
     ! Record with kg_k(3,npw_k) vectors
     nfrec = nfrec + 1
     bsz_frec(nfrec) = 3*npw_k*xmpi_bsize_int
     !
     if (formeig==0) then 
       ! Record with the eigenvalues
       ! eig_k(nband_k), occ_k(nband_k)
       nfrec = nfrec + 1
       bsz_frec(nfrec) = 2*nband_k*xmpi_bsize_dp
       !
       ! cg_k record
       do band=1,nband_k
         nfrec = nfrec + 1
         bsz_frec(nfrec) = 2*npw_k*Hdr%nspinor*xmpi_bsize_dp
       end do

     else if (formeig==1) then
       do band=1,nband_k
         ! Record with the eigenvalues
         nfrec = nfrec + 1
         bsz_frec(nfrec) = 2*nband_k*xmpi_bsize_dp

         ! cg_k record
         nfrec = nfrec + 1
         bsz_frec(nfrec) = 2*npw_k*Hdr%nspinor*xmpi_bsize_dp
       end do
     else 
       MSG_ERROR("Wrong formeig")
     end if
     
   end do
 end do

 ABI_MALLOC(bsize_frecords, (nfrec))
 bsize_frecords = bsz_frec(1:nfrec)

 ABI_FREE(bsz_frec)

end subroutine hdr_bsize_frecords
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_io_wfftype
!! NAME
!! hdr_io_wfftype
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hdr_type
!! structured variables (read/write/echo).
!! According to the value of rdwr, it reads the header
!! of a file, writes it, or echo the value of the structured
!! variable to a file.
!! Note that, when reading, different records of hdr
!! are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated
!! correctly by a call to hdr_free when hdr is not used anymore.
!! Two instances of the hdr_io routines are defined :
!!  hdr_io_int to which only the unit number is given
!!  hdr_io_wfftype to which a wffil datatype is given
!!
!! INPUTS
!!  rdwr= if 1, read the hdr structured variable from the header of the file,
!!        if 2, write the header to unformatted file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the hdr without rewinding (unformatted)
!!        if 6, write the hdr without rewinding (unformatted)
!!  unitfi=unit number of the file (unformatted if rdwr=1, 2, 5 or 6 formatted if rdwr=3,4)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  hdr <type(hdr_type)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!! In all cases, the file is supposed to be open already
!! When reading (rdwr=1) or writing (rdwr=2), rewind the file
!! When echoing (rdwr=3) does not rewind the file.
!! When reading (rdwr=5) or writing (rdwr=6), DOES NOT rewind the file
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_io_wfftype(fform,hdr,rdwr,wff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_io_wfftype'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr
 type(hdr_type),intent(inout) :: hdr
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: ierr
#endif

! *************************************************************************

 DBG_ENTER("COLL")

 if ( wff%accesswff==IO_MODE_FORTRAN .or. &
& (wff%accesswff==IO_MODE_FORTRAN_MASTER .and.wff%master==wff%me).or. &
& (wff%accesswff==IO_MODE_MPI  .and.wff%master==wff%me)    ) then
   call hdr_io_int(fform,hdr,rdwr,wff%unwff)
   ! Master node **MUST** flush the output buffer so that the 
   ! other nodes can read headform and therefore the Fortran marker length when MPI-IO is used
   if (rdwr == 2) then
     call flush_unit(wff%unwff)
   end if
 end if

#if defined HAVE_MPI
!In the parallel case, if the files were not local, need to bcast the data
 if(rdwr==1)then
   if (wff%accesswff==IO_MODE_FORTRAN_MASTER .or. wff%accesswff==IO_MODE_MPI) then
     if (wff%spaceComm/=MPI_COMM_SELF) then
       call MPI_BCAST(fform,1,MPI_INTEGER,wff%master,wff%spaceComm,ierr)
       call hdr_comm(hdr,wff%master,wff%me,wff%spaceComm)
     end if
     wff%headform=hdr%headform
     if(wff%accesswff==IO_MODE_MPI)then
       call hdr_skip_wfftype(wff,ierr)
     end if
   end if
 end if
#if defined HAVE_MPI_IO
 if (rdwr == 2 .and. wff%accesswff==IO_MODE_MPI) then
   if (wff%spaceComm/=MPI_COMM_SELF) then
     call xmpi_barrier(wff%spaceComm)
   end if
   wff%headform=hdr%headform
   call hdr_skip_wfftype(wff,ierr)
 end if
#endif
 if (rdwr==5) wff%headform=hdr%headform
#else
 if (rdwr==1.or.rdwr==5) wff%headform=hdr%headform
#endif

 DBG_EXIT("COLL")

end subroutine hdr_io_wfftype
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_io_int
!! NAME
!! hdr_io_int
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hdr_type structured variables (read/write/echo).
!! According to the value of rdwr, it reads the header of a file, writes it, or echo the value of the structured
!! variable to a file. Note that, when reading, different records of hdr are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated correctly by a call to hdr_free when hdr is not used anymore.
!! Two instances of the hdr_io routines are defined :
!!   hdr_io_int to which only the unit number is given
!!   hdr_io_wfftype to which a wffil datatype is given
!!
!! INPUTS
!!  rdwr= if 1, read the hdr structured variable from the header of the file,
!!        if 2, write the header to unformatted file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the hdr without rewinding (unformatted)
!!        if 6, write the hdr without rewinding (unformatted)
!!  unitfi=unit number of the file (unformatted if rdwr=1, 2, 5 or 6 formatted if rdwr=3,4)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  hdr <type(hdr_type)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!! In all cases, the file is supposed to be open already
!! When reading (rdwr=1) or writing (rdwr=2), rewind the file
!! When echoing (rdwr=3) does not rewind the file.
!! When reading (rdwr=5) or writing (rdwr=6), DOES NOT rewind the file
!!
!! PARENTS
!!      m_header
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_io_int(fform,hdr,rdwr,unitfi)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_io_int'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr,unitfi
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
 integer :: ierr

!*************************************************************************

 DBG_ENTER("COLL")

 if(rdwr==1 .or. rdwr==5)then
   ! Reading the header of an unformatted file
    call hdr_fort_read(Hdr,unitfi,fform,rewind=(rdwr==1))

 else if (rdwr==2 .or. rdwr==6)then
   ! Writing the header of an unformatted file
   call hdr_fort_write(Hdr,unitfi,fform,ierr,rewind=(rdwr==2))

 else if (rdwr==3 .or. rdwr==4) then
   !  Writing the header of a formatted file
   call hdr_echo(Hdr,fform,rdwr,unit=unitfi)
 end if 

 DBG_EXIT("COLL")

end subroutine hdr_io_int
!!***
!----------------------------------------------------------------------

!!****f* m_header/hdr_echo
!! NAME
!! hdr_echo
!!
!! FUNCTION
!! Echo the header 
!!
!! INPUTS
!!  hdr <type(hdr_type)>=the header structured variable
!!  rdwr= if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!  fform=kind of the array in the file
!!  [unit]=unit number of the formatted file [DEFAULT: std_out]
!!
!! OUTPUT
!!  Only writing 
!!
!! PARENTS
!!      ioprof,m_header,m_wfk
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_echo(Hdr,fform,rdwr,unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_echo'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr 
 integer,optional,intent(in) :: unit
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
 integer :: iatom,ii,ikpt,ipsp,isym,ount

!*************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit

 if (ount == dev_null) RETURN

 write(ount,'(a)')' ==============================================================================='
 if (rdwr==3) write(ount, '(a)' ) ' ECHO of part of the ABINIT file header '
 if (rdwr==4) write(ount, '(a)' ) ' ECHO of the ABINIT file header '
 write(ount, '(a)' ) ' '
 write(ount, '(a)' ) ' First record :'
 write(ount, '(a,a6,2i5)' )  '.codvsn,headform,fform = ',hdr%codvsn, hdr%headform, fform   ! Do not worry about 22 format

 write(ount, '(a)' ) ' '
 write(ount, '(a)' ) ' Second record :'
 write(ount, '(a,4i6)') ' bantot,intxc,ixc,natom  =',hdr%bantot, hdr%intxc, hdr%ixc, hdr%natom
 write(ount, '(a,4i6)') ' ngfft(1:3),nkpt         =',hdr%ngfft(1:3), hdr%nkpt

 if (hdr%headform>=23) then
   write(ount, '(a,2i6)') ' nspden,nspinor          =',hdr%nspden, hdr%nspinor
 end if

 if (hdr%headform<=23) then
   write(ount, '(a,4i6)' ) ' nsppol,nsym,ntypat,occopt=',hdr%nsppol,hdr%nsym,hdr%ntypat,hdr%occopt
 else if(hdr%headform<=40)then
   write(ount, '(a,5i6)' ) ' nsppol,nsym,npsp,ntypat,occopt=',hdr%nsppol,hdr%nsym,hdr%npsp,hdr%ntypat,hdr%occopt

 else if (hdr%headform==41 .or. hdr%headform==42) then
   write(ount, '(a,6i6)' )&
&   ' nsppol,nsym,npsp,ntypat,occopt,pertcase=',hdr%nsppol,hdr%nsym,hdr%npsp,hdr%ntypat,hdr%occopt,hdr%pertcase
 else if (hdr%headform>=44) then
   write(ount, '(a,4i6)' ) ' nsppol,nsym,npsp,ntypat =',hdr%nsppol,hdr%nsym,hdr%npsp,hdr%ntypat
   write(ount, '(a,3i6)' ) ' occopt,pertcase,usepaw  =',hdr%occopt,hdr%pertcase,hdr%usepaw
 end if

 if(hdr%headform==40 .or. hdr%headform==41 .or. hdr%headform==42)then
   write(ount, '(a,2es18.10)') ' ecut,ecutsm             =',hdr%ecut, hdr%ecutsm
 else if(hdr%headform>=44)then
   write(ount, '(a,3es18.10)') ' ecut,ecutdg,ecutsm      =',hdr%ecut, hdr%ecutdg, hdr%ecutsm
 end if

 write(ount, '(a, es18.10)' ) ' ecut_eff                =',hdr%ecut_eff

 if(hdr%headform>=41)then
   write(ount, '(a,3es18.10)') ' qptn(1:3)               =',hdr%qptn(1:3)
 end if

 write(ount, '(a,3es18.10)' ) ' rprimd(1:3,1)           =',hdr%rprimd(1:3,1)
 write(ount, '(a,3es18.10)' ) ' rprimd(1:3,2)           =',hdr%rprimd(1:3,2)
 write(ount, '(a,3es18.10)' ) ' rprimd(1:3,3)           =',hdr%rprimd(1:3,3)

 if(hdr%headform==40.or.hdr%headform==41)then
   write(ount, '(a,2es18.10)') ' tphysel,tsmear          =',hdr%tphysel, hdr%tsmear
 else if(hdr%headform>=42)then
   write(ount, '(a,3es18.10)') ' stmbias,tphysel,tsmear  =',hdr%stmbias,hdr%tphysel, hdr%tsmear
 end if

 write(ount, '(a)' )
 if (rdwr==3)then
   write(ount, '(a,i3,a)' ) ' The header contain ',hdr%npsp+2,' additional records.'
 else
   write(ount, '(a)' ) ' Third record :'
   write(ount, '(a,(12i4,8x))') ' istwfk=',hdr%istwfk(:)
   write(ount, '(a,(12i4,8x))') ' nband =',hdr%nband(:)
   write(ount, '(a,(10i5,8x))') ' npwarr=',hdr%npwarr(:)

   if(hdr%headform>=40)then
     write(ount, '(a,(12i4,8x))') ' so_psp=',hdr%so_psp(:)
   end if

   if(hdr%headform>=40)then
     write(ount, '(a)') ' symafm='
     write(ount, '(8x,24i3,8x)') hdr%symafm(:)
   end if

   write(ount, '(a)' ) ' symrel='
   do isym=1,hdr%nsym/2
     write(ount, '(a,9i4,a,9i4)' ) '        ',hdr%symrel(:,:,2*isym-1),'  ',hdr%symrel(:,:,2*isym)
   end do
   if(2*(hdr%nsym/2)/=hdr%nsym)write(ount, '(a,9i4)' ) '        ',hdr%symrel(:,:,hdr%nsym)

   write(ount, '(a,(12i4,8x))') ' type  =',hdr%typat(:)
   write(ount, '(a)' ) ' kptns =                 (max 50 k-points will be written)'
   do ikpt=1,min(hdr%nkpt,50)
     write(ount, '(a,3es16.6)' ) '        ',hdr%kptns(:,ikpt)
   end do
   write(ount, '(a)' ) ' wtk ='
   do ikpt=1,hdr%nkpt,10
     write(ount, '(a,10f6.2)' ) '        ',hdr%wtk(ikpt:min(hdr%nkpt,ikpt + 10 - 1))
   end do
   write(ount, '(a)' ) '   occ ='
   do ii=1,hdr%bantot,10
     write(ount, '(a,10f6.2)') '        ',hdr%occ(ii:min(hdr%bantot,ii+10-1))
   end do
   write(ount, '(a)' ) ' tnons ='
   do isym=1,hdr%nsym/2
     write(ount, '(a,3f10.6,a,3f10.6)' ) '        ',hdr%tnons(:,2*isym-1),'  ',hdr%tnons(:,2*isym)
   end do
   if(2*(hdr%nsym/2)/=hdr%nsym)write(ount, '(a,3f10.6)' ) '        ',hdr%tnons(:,hdr%nsym)
   write(ount, '(a,(10f6.2,8x))') '  znucl=',hdr%znucltypat(:)
   write(ount,'(a)')

   write(ount, '(a)' ) ' Pseudopotential info :'
   do ipsp=1,hdr%npsp
     write(ount,'(a,a)' ) ' title=',trim(hdr%title(ipsp))
     write(ount,'(a,f6.2,a,f6.2,a,i3,a,i6,a,i3,a,i3)' ) &
&     '  znuclpsp=',hdr%znuclpsp(ipsp),    ', zionpsp=',  hdr%zionpsp(ipsp),&
&     ', pspso=' , hdr%pspso(ipsp),  ', pspdat=',hdr%pspdat(ipsp),          &
&     ', pspcod=', hdr%pspcod(ipsp), ', pspxc=', hdr%pspxc(ipsp)
     if(hdr%headform>=44)then
       if(hdr%usepaw==1)then
         write(ount,'(a,i3)' ) '  lmn_size=', hdr%lmn_size(ipsp)
       else
         write(ount,'(a,i3)' ) '  lmnmax  =', hdr%lmn_size(ipsp)
       end if
     end if

   end do

   write(ount, '(a)' ) ' '
   write(ount, '(a)' ) ' Last record :'
   write(ount, '(a,es16.6,es22.12,es16.6)' )' residm,etot,fermie=',hdr%residm, hdr%etot, hdr%fermie
   write(ount, '(a)' ) ' xred ='
   do iatom=1,hdr%natom
     write(ount, '(a,3es16.6)' ) '        ',hdr%xred(:,iatom)
   end do

   if (hdr%usepaw==1)then
     call pawrhoij_io(hdr%pawrhoij,ount,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,hdr%headform,"Echo")
   end if

   if (rdwr==3)write(ount, '(a)' ) ' End the ECHO of part of the ABINIT file header '
   if (rdwr==4)write(ount, '(a)' ) ' End the ECHO of the ABINIT file header '
   write(ount,'(a)')' ==============================================================================='
 end if ! rdwr is 3 or 4

 call flush_unit(ount)

end subroutine hdr_echo
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_io_etsf
!! NAME
!! hdr_io_etsf
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hdr_type structured variables (read/write/echo).
!! It handles variables according to the ETSF format, whenever
!! possible and uses new variables when not available in the ETSF format.
!! According to the value of rdwr, it reads the header
!! of a file, writes it, or echo the value of the structured
!! variable to a file.
!! Note that, when reading, different records of hdr
!! are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated
!! correctly by a call to hdr_free when hdr is not used anymore.
!!
!! INPUTS
!!  rdwr= if 1, read the hdr structured variable from the header of the netCDF file,
!!        if 2, write the header to unformatted netCDF file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the hdr without rewinding (unformatted), identical to 1 for netCDF
!!        if 6, read the hdr without rewinding (unformatted), identical to 2 for netCDF
!!  unitwff=the unit of the open NetCDF file.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  hdr <type(hdr_type)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!!
!! PARENTS
!!      cut3d,inwffil,ioarr,kss2wfk,m_ebands,m_io_kss,m_io_screening,outwf
!!      pawmkaewf
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_io_etsf(fform,hdr,rdwr,unitwff,nc_define)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_io_etsf'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rdwr,unitwff
 integer,intent(inout) :: fform
 logical,optional,intent(in) :: nc_define
 type(hdr_type),target,intent(inout) :: hdr

#ifdef HAVE_TRIO_ETSF_IO
!Local variables-------------------------------
!scalars
 logical :: my_define

! *************************************************************************

 if (rdwr==1 .or. rdwr==5) then
   ! Reading the header from a Netcdf  file
   call hdr_ncread(Hdr,unitwff,fform)

 else if (rdwr==2 .or. rdwr==6) then
   !  Writing the header to a Netcdf file

   my_define = .FALSE.; if (PRESENT(nc_define)) my_define = nc_define
   call hdr_ncwrite(hdr,unitwff,fform,nc_define=my_define)

 else if (rdwr==3 .or. rdwr==4) then
   ! Echo the header 
   call hdr_io_int(fform, hdr, rdwr, unitwff)
 end if 

#else
 MSG_ERROR("ETSF-IO support not activated")
#endif

end subroutine hdr_io_etsf
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_skip_int
!! NAME
!! hdr_skip_int
!!
!! FUNCTION
!! Skip wavefunction or density file header, after having rewound the file.
!! Two instances of the hdr_skip routines are defined :
!!  hdr_skip_int to which only the unit number is given
!!  hdr_skip_wfftype to which a wffil datatype is given
!!
!! INPUTS
!!  unit = number of unit to be read
!!
!! OUTPUT
!!  ierr = error code returned by the MPI calls
!!
!! SIDE EFFECTS
!!
!! NOTES
!! No checking performed, since hdr_skip is assumed to be used only
!! on temporary wavefunction files.
!! This initialize further reading and checking by rwwf
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_skip_int(unitfi,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_skip_int'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unitfi
 integer,intent(out) :: ierr

!Local variables-------------------------------
 type(wffile_type) :: wff

! *************************************************************************

!Use default values for wff
 wff%unwff=unitfi
 wff%accesswff=IO_MODE_FORTRAN
 wff%me=0
 wff%master=0
!Then, transmit to hdr_skip_wfftype
 call hdr_skip_wfftype(wff,ierr)

end subroutine hdr_skip_int
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_skip_wfftype
!! NAME
!! hdr_skip_wfftype
!!
!! FUNCTION
!! Skip wavefunction or density file header, after having rewound the file.
!! Two instances of the hdr_skip routines are defined :
!!  hdr_skip_int to which only the unit number is given
!!  hdr_skip_wfftype to which a wffil datatype is given
!!
!! INPUTS
!!  unit = number of unit to be read
!!
!! OUTPUT
!!  ierr = error code returned by the MPI calls
!!
!! NOTES
!! No checking performed, since hdr_skip is assumed to be used only
!! on temporary wavefunction files.
!! This initialize further reading and checking by rwwf
!!
!! PARENTS
!!      m_header
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_skip_wfftype(wff,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_skip_wfftype'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer, intent(out) :: ierr

!Local variables-------------------------------
 integer :: headform,mu,npsp,unit,usepaw
 integer :: integers(17)
 character(len=6) :: codvsn
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit,positloc
 integer :: statux(MPI_STATUS_SIZE)
#endif

!*************************************************************************

 unit=wff%unwff
 ierr=0

 if( wff%accesswff==IO_MODE_FORTRAN .or. (wff%accesswff==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me) ) then

   rewind (unit)

!  Pick off headform from WF file
   read(unit) codvsn,headform                  ! XG040806 This does not work, but I do not understand why ?!
!  read(unit) integers(1),headform             ! This works ...
!  ! MT012408 Because codvsn is char*6 and not an integer !
   if(headform==1   .or. headform==2   .or. &
&   headform==51  .or. headform==52  .or.   &
&   headform==101 .or. headform==102       ) headform=22

   if (headform<44) then
     read (unit) integers(1:13),npsp
   else
     read (unit) integers(1:13),npsp,integers(15:17),usepaw
   end if

!  Skip rest of header records
   do mu=1,2+npsp
     read (unit)
   end do
   if ((headform>=44).and.(usepaw==1)) then
     read (unit)
     read (unit)
   end if

#if defined HAVE_MPI_IO
 else if(wff%accesswff==IO_MODE_MPI)then

   headform=wff%headform
   if(headform==1   .or. headform==2   .or. &
&   headform==51  .or. headform==52  .or. &
&   headform==101 .or. headform==102) headform=22

!  Causes all previous writes to be transferred to the storage device
   call flush_unit(wff%unwff)
   call MPI_FILE_SYNC(wff%fhwff,ierr)

!  Check FORTRAN record marker length (only at first call)
   if (wff%nbOct_recMarker<=0) then
     call getRecordMarkerLength_wffile(wff)
   end if

   if (wff%master==wff%me) then

!    Reading the first record of the file -------------------------------------
!    read (unitfi)   codvsn,headform,..............
     posit = 0
     call rwRecordMarker(1,posit,delim_record,wff,ierr)

!    Reading the second record of the file ------------------------------------
!    read(unitfi) bantot, hdr%date, hdr%intxc.................
!    Pick off npsp and usepaw from WF file
     positloc  = posit + wff%nbOct_recMarker + wff%nbOct_int*13
     call MPI_FILE_READ_AT(wff%fhwff,positloc,npsp,1,MPI_INTEGER,statux,ierr)
!    call MPI_FILE_READ_AT_ALL(wff%fhwff,positloc,npsp,1,MPI_INTEGER,statux,ierr)
     if (headform >= 44) then
       positloc = positloc +  wff%nbOct_int*4
       call MPI_FILE_READ_AT(wff%fhwff,positloc,usepaw,1,MPI_INTEGER,statux,ierr)
!      call MPI_FILE_READ_AT_ALL(wff%fhwff,positloc,usepaw,1,MPI_INTEGER,statux,ierr)
     end if
     call rwRecordMarker(1,posit,delim_record,wff,ierr)

!    Reading the rest of the file ---------------------------------------------
     do mu=1,2+npsp
       call rwRecordMarker(1,posit,delim_record,wff,ierr)
     end do
     if ((headform>=44).and.(usepaw==1)) then
       call rwRecordMarker(1,posit,delim_record,wff,ierr)
       call rwRecordMarker(1,posit,delim_record,wff,ierr)
     end if

     wff%offwff=posit
   end if

   if (wff%spaceComm/=MPI_COMM_SELF) then
     call MPI_BCAST(wff%offwff,1,wff%offset_mpi_type,wff%master,wff%spaceComm,ierr)
   end if
#endif
 end if

end subroutine hdr_skip_wfftype
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_update
!! NAME
!! hdr_update
!!
!! FUNCTION
!! This subroutine update the header structured datatype.
!! Most of its records had been initialized correctly, but some corresponds
!! to evolving variables, or change with the context (like fform),
!! This routine is to be called before writing the header
!! to a file, in order to have up-to-date information.
!!
!! INPUTS
!! bantot=total number of bands
!! etot=total energy (Hartree)
!! fermie=Fermi energy (Hartree)
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! mpi_comm_atom=--optional-- MPI communicator over atoms
!! natom=number of atoms
!! residm=maximal residual
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! occ(bantot)=occupancies for each band and k point
!! pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      afterscfloop,gstate,loper3,nonlinear,respfn,scfcv,setup_bse
!!      setup_screening,setup_sigma
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_update(bantot,etot,fermie,hdr,natom,residm,rprimd,occ,pawrhoij,usepaw,xred, &
&                     mpi_comm_atom,mpi_atmtab) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_update'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot,natom,usepaw
 integer,optional,intent(in) :: mpi_comm_atom
 real(dp),intent(in) :: etot,fermie,residm
 type(hdr_type),intent(inout) :: hdr !vz_i
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: occ(bantot),rprimd(3,3),xred(3,natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------

! *************************************************************************

!Update of the "evolving" data
 hdr%etot     =etot
 hdr%fermie   =fermie
 hdr%residm   =residm
 hdr%rprimd(:,:)=rprimd(:,:)
 hdr%occ(:)   =occ(:)
 hdr%xred(:,:)=xred(:,:)

 if (usepaw==1) then
   if (present(mpi_comm_atom)) then
     if (present(mpi_atmtab)) then
       call pawrhoij_copy(pawrhoij,hdr%pawrhoij,mpi_comm_atom=mpi_comm_atom,mpi_atmtab=mpi_atmtab)
     else
       call pawrhoij_copy(pawrhoij,hdr%pawrhoij,mpi_comm_atom=mpi_comm_atom)
     end if
   else
     call pawrhoij_copy(pawrhoij,hdr%pawrhoij)
   end if
 end if

end subroutine hdr_update
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_comm
!! NAME
!! hdr_comm
!!
!! FUNCTION
!! This subroutine transmits the header structured datatype
!! initialized on one processor (or a group of processor),
!! to the other processors. It also allocate the needed
!! part of the header.
!!
!! INPUTS
!!  master = id of the master process
!!  me = id of the current process
!!  comm = id of the space communicator handler
!!
!! OUTPUT
!!  (no output)
!!
!! SIDE EFFECTS
!!  hdr <type(hdr_type)>=the header. For the master, it is already
!!   initialized entirely, while for the other procs, everything has
!!   to be transmitted.
!!
!! NOTES
!! This routine is called only in the case of MPI version of the code.
!!
!! PARENTS
!!      elphon,initaim,m_header,m_io_kss,m_io_screening,m_wfk,read_gkk
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_comm(hdr,master,me,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_comm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: master,me,comm
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
 integer :: bantot,cplex,iatom,ierr,index,index2,ipsp,ispden,list_size,list_size2,natom,nkpt
 integer :: npsp,nsel,nspden,nsppol,nsym,nrhoij,ntypat
 integer,allocatable :: list_int(:)
 real(dp),allocatable :: list_dpr(:)
 character(len=fnlen) :: list_tmp
 character(len=fnlen),allocatable :: list_char(:)

! *************************************************************************

 if (xcomm_size(comm) == 1) return ! Nothing to do

 DBG_ENTER("COLL")

!Transmit the integer scalars
 list_size=20
 ABI_MALLOC(list_int,(list_size))
 if (master==me)then
   list_int(1)=hdr%bantot
   list_int(2)=hdr%date
   list_int(3)=hdr%headform
   list_int(4)=hdr%intxc
   list_int(5)=hdr%ixc
   list_int(6)=hdr%natom
   list_int(7)=hdr%nkpt
   list_int(8)=hdr%npsp
   list_int(9)=hdr%nspden
   list_int(10)=hdr%nspinor
   list_int(11)=hdr%nsppol
   list_int(12)=hdr%nsym
   list_int(13)=hdr%ntypat
   list_int(14)=hdr%occopt
   list_int(15)=hdr%pertcase
   list_int(16)=hdr%usepaw
   list_int(17:19)=hdr%ngfft(1:3)
   list_int(20)=hdr%usewvl
 end if

 call xmpi_bcast(list_int,master,comm,ierr)

 if(master/=me)then
   hdr%bantot  =list_int(1)
   hdr%date    =list_int(2)
   hdr%headform=list_int(3)
   hdr%intxc   =list_int(4)
   hdr%ixc     =list_int(5)
   hdr%natom   =list_int(6)
   hdr%nkpt    =list_int(7)
   hdr%npsp    =list_int(8)
   hdr%nspden  =list_int(9)
   hdr%nspinor =list_int(10)
   hdr%nsppol  =list_int(11)
   hdr%nsym    =list_int(12)
   hdr%ntypat  =list_int(13)
   hdr%occopt  =list_int(14)
   hdr%pertcase=list_int(15)
   hdr%usepaw  =list_int(16)
   hdr%ngfft(1:3)=list_int(17:19)
   hdr%usewvl  =list_int(20)
 end if
 ABI_FREE(list_int)

 bantot=hdr%bantot
 natom =hdr%natom
 nkpt  =hdr%nkpt
 npsp  =hdr%npsp
 nspden=hdr%nspden
 nsppol=hdr%nsppol
 nsym  =hdr%nsym
 ntypat=hdr%ntypat

 if(master/=me)then
!  Allocate all components of hdr
   ABI_MALLOC(hdr%istwfk,(nkpt))
   ABI_MALLOC(hdr%nband,(nkpt*nsppol))
   ABI_MALLOC(hdr%npwarr,(nkpt))
   ABI_MALLOC(hdr%pspcod,(npsp))
   ABI_MALLOC(hdr%pspdat,(npsp))
   ABI_MALLOC(hdr%pspso,(npsp))
   ABI_MALLOC(hdr%pspxc,(npsp))
   ABI_MALLOC(hdr%lmn_size,(npsp))
   ABI_MALLOC(hdr%so_psp,(npsp))
   ABI_MALLOC(hdr%symafm,(nsym))
   ABI_MALLOC(hdr%symrel,(3,3,nsym))
   ABI_MALLOC(hdr%typat,(natom))
   ABI_MALLOC(hdr%kptns,(3,nkpt))
   ABI_MALLOC(hdr%occ,(bantot))
   ABI_MALLOC(hdr%tnons,(3,nsym))
   ABI_MALLOC(hdr%wtk,(nkpt))
   ABI_MALLOC(hdr%xred,(3,natom))
   ABI_MALLOC(hdr%zionpsp,(npsp))
   ABI_MALLOC(hdr%znuclpsp,(npsp))
   ABI_MALLOC(hdr%znucltypat,(ntypat))
   ABI_MALLOC(hdr%title,(npsp))
 end if

!Transmit the integer arrays
 list_size=nkpt*(2+nsppol)+6*npsp+10*nsym+natom
 ABI_MALLOC(list_int,(list_size))
 if (master==me)then
   list_int(1      :nkpt             )=hdr%istwfk ; index=nkpt
   list_int(1+index:nkpt*nsppol+index)=hdr%nband  ; index=index+nkpt*nsppol
   list_int(1+index:nkpt       +index)=hdr%npwarr ; index=index+nkpt
   list_int(1+index:npsp       +index)=hdr%pspcod ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%pspdat ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%pspso  ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%pspxc  ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%lmn_size ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%so_psp ; index=index+npsp
   list_int(1+index:nsym       +index)=hdr%symafm ; index=index+nsym
   list_int(1+index:nsym*3*3   +index)=reshape(hdr%symrel,(/3*3*nsym/))
   index=index+nsym*3*3
   list_int(1+index:natom      +index)=hdr%typat   ; index=index+natom
 end if

 call xmpi_bcast(list_int,master,comm,ierr)

 if(master/=me)then
   hdr%istwfk=list_int(1      :nkpt             ) ; index=nkpt
   hdr%nband =list_int(1+index:nkpt*nsppol+index) ; index=index+nkpt*nsppol
   hdr%npwarr=list_int(1+index:nkpt       +index) ; index=index+nkpt
   hdr%pspcod=list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%pspdat=list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%pspso =list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%pspxc =list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%lmn_size=list_int(1+index:npsp     +index) ; index=index+npsp
   hdr%so_psp =list_int(1+index:npsp   +index) ; index=index+npsp
   hdr%symafm=list_int(1+index:nsym       +index) ; index=index+nsym
   hdr%symrel=reshape(list_int(1+index:nsym*3*3   +index),(/3,3,nsym/))
   index=index+nsym*3*3
   hdr%typat  =list_int(1+index:natom      +index) ; index=index+natom
 end if
 ABI_FREE(list_int)

!Transmit the double precision scalars and arrays
 list_size=21+3*nkpt+nkpt+bantot+3*nsym+3*natom+2*npsp+ntypat
 ABI_MALLOC(list_dpr,(list_size))
 if (master==me)then
   list_dpr(1)=hdr%ecut_eff
   list_dpr(2)=hdr%etot
   list_dpr(3)=hdr%fermie
   list_dpr(4)=hdr%residm
   list_dpr(5:13)=reshape(hdr%rprimd(1:3,1:3),(/9/))
   list_dpr(14)=hdr%ecut
   list_dpr(15)=hdr%ecutdg
   list_dpr(16)=hdr%ecutsm
   list_dpr(17)=hdr%tphysel
   list_dpr(18)=hdr%tsmear
   list_dpr(19:21)=hdr%qptn(1:3)                                 ; index=21
   list_dpr(1+index:3*nkpt +index)=reshape(hdr%kptns,(/3*nkpt/)) ; index=index+3*nkpt
   list_dpr(1+index:nkpt   +index)=hdr%wtk                       ; index=index+nkpt
   list_dpr(1+index:bantot +index)=hdr%occ                       ; index=index+bantot
   list_dpr(1+index:3*nsym +index)=reshape(hdr%tnons,(/3*nsym/)) ; index=index+3*nsym
   list_dpr(1+index:3*natom+index)=reshape(hdr%xred,(/3*natom/)) ; index=index+3*natom
   list_dpr(1+index:npsp   +index)=hdr%zionpsp                   ; index=index+npsp
   list_dpr(1+index:npsp   +index)=hdr%znuclpsp                  ; index=index+npsp
   list_dpr(1+index:ntypat  +index)=hdr%znucltypat               ; index=index+ntypat
 end if

 call xmpi_bcast(list_dpr,master,comm,ierr)

 if(master/=me)then
   hdr%ecut_eff=list_dpr(1)
   hdr%etot    =list_dpr(2)
   hdr%fermie  =list_dpr(3)
   hdr%residm  =list_dpr(4)
   hdr%rprimd  =reshape(list_dpr(5:13),(/3,3/))
   hdr%ecut    =list_dpr(14)
   hdr%ecutdg  =list_dpr(15)
   hdr%ecutsm  =list_dpr(16)
   hdr%tphysel =list_dpr(17)
   hdr%tsmear  =list_dpr(18)
   hdr%qptn(1:3)=list_dpr(19:21)                                    ; index=21
   hdr%kptns   =reshape(list_dpr(1+index:3*nkpt +index),(/3,nkpt/)) ; index=index+3*nkpt
   hdr%wtk     =list_dpr(1+index:nkpt   +index)                     ; index=index+nkpt
   hdr%occ     =list_dpr(1+index:bantot +index)                     ; index=index+bantot
   hdr%tnons   =reshape(list_dpr(1+index:3*nsym +index),(/3,nsym/)) ; index=index+3*nsym
   hdr%xred    =reshape(list_dpr(1+index:3*natom+index),(/3,natom/)); index=index+3*natom
   hdr%zionpsp =list_dpr(1+index:npsp   +index)                     ; index=index+npsp
   hdr%znuclpsp=list_dpr(1+index:npsp   +index)                     ; index=index+npsp
   hdr%znucltypat=list_dpr(1+index:ntypat  +index)                  ; index=index+ntypat
 end if
 ABI_FREE(list_dpr)

!Transmit the characters
 list_size=npsp+1
 ABI_MALLOC(list_char,(list_size))
 if (master==me)then
   list_char(1)       =hdr%codvsn  ! Only 6 characters are stored in list_char(1)
   list_char(2:npsp+1)=hdr%title
 end if

 call xmpi_bcast(list_char,master,comm,ierr)

 if(master/=me)then
   list_tmp=list_char(1)
   hdr%codvsn=list_tmp(1:6)
   do ipsp=2,npsp+1
     list_tmp =list_char(ipsp)
     hdr%title =list_tmp(1:fnlen)
   end do
 end if
 ABI_FREE(list_char)

!Transmit the structured variables in case of PAW
 if (hdr%usepaw==1) then

   nrhoij=0
   if (master==me)then
     cplex=hdr%pawrhoij(1)%cplex
     nspden=hdr%pawrhoij(1)%nspden
     do iatom=1,natom
       nrhoij=nrhoij+hdr%pawrhoij(iatom)%nrhoijsel
     end do
   end if

   call xmpi_bcast(nrhoij,master,comm,ierr)
   call xmpi_bcast(cplex ,master,comm,ierr)
   call xmpi_bcast(nspden,master,comm,ierr)

   list_size=natom+nrhoij;list_size2=nspden*nrhoij*cplex
   ABI_MALLOC(list_int,(list_size))
   ABI_MALLOC(list_dpr,(list_size2))
   if (master==me)then
     index=0;index2=0
     do iatom=1,natom
       nsel=hdr%pawrhoij(iatom)%nrhoijsel
       list_int(1+index)=nsel
       list_int(2+index:1+nsel+index)=hdr%pawrhoij(iatom)%rhoijselect(1:nsel)
       index=index+1+nsel
       do ispden=1,nspden
         list_dpr(1+index2:nsel*cplex+index2)=hdr%pawrhoij(iatom)%rhoijp(1:nsel*cplex,ispden)
         index2=index2+nsel*cplex
       end do
     end do
   end if

   call xmpi_bcast(list_int,master,comm,ierr)
   call xmpi_bcast(list_dpr,master,comm,ierr)

   if(master/=me)then
     index=0;index2=0
     ABI_DT_MALLOC(hdr%pawrhoij,(natom))
     call pawrhoij_alloc(hdr%pawrhoij,cplex,nspden,hdr%nspinor,hdr%nsppol,hdr%typat,lmnsize=hdr%lmn_size)
     do iatom=1,natom
       nsel=list_int(1+index)
       hdr%pawrhoij(iatom)%nrhoijsel=nsel
       hdr%pawrhoij(iatom)%rhoijselect(1:nsel)=list_int(2+index:1+nsel+index)
       index=index+1+nsel
       do ispden=1,nspden
         hdr%pawrhoij(iatom)%rhoijp(1:nsel*cplex,ispden)=list_dpr(1+index2:nsel*cplex+index2)
         index2=index2+nsel*cplex
       end do
     end do
   end if
   ABI_FREE(list_int)
   ABI_FREE(list_dpr)
 end if

 DBG_EXIT("COLL")

end subroutine hdr_comm
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_check
!! NAME
!! hdr_check
!!
!! FUNCTION
!! This subroutine compare the header structured variable (hdr)
!! from input data (mostly dtset and psps) with the one (hdr0) of
!! an input data file (e.g. wf, density, potential).
!! Various values are checked for agreement or near agreement in the
!! case of floating point numbers.  The program will exit or produce
!! warning messages when unexpected values are found.
!! A record of the comparison of the headers is written to stdout.
!!
!! Decisions have been taken about whether a restart is allowed.
!! In the self-consistent case, a restart will always be allowed, but
!! one has to distinguish between a direct restart and a restart with
!! translation of wavefunction.
!! In the non-self-consistent case, the conditions below
!! must be fulfilled to allow a restart.
!!
!! INPUTS
!!  fform=integer specification of data type (expected)
!!  fform0=integer specification of data type (from disk file)
!!  mode_paral : COLL or PERS, for all leave_new and wrtout calls
!!  hdr <type(hdr_type)>=the header structured variable from dtset and psps
!!  hdr0<type(hdr_type)>=the header structured variable from the disk file
!!
!! OUTPUT
!!  restart=1 if direct restart, =2 if translation is needed, =0 if no
!!              restart is possible.
!!  restartpaw= deals with the additional informations in the PAW method
!!              =1 if direct restart, =0 if no restart from spherical data is possible.
!!              also 0 if no PAW
!!
!! NOTES
!! In the current version of the user interface restarts are allowed from
!! wavefunction files for self-consistent runs and from densities for
!! non-self-consistent runs. The precise conditions under which we will
!! allow a restart in this release are as follows.
!!
!!           self-consistent case : direct restarts
!!           ======================================
!!
!! A direct restart will be allowed provided the following quantities in
!! old and new calculations are the same:
!!
!!   (A) the primitive vectors                             (tprim)
!!   (B) the plane-wave cutoff                             (tecut)
!!   (C) nkpt, kpt(3,nkpt), wtk(nkpt)                      (tkpt)
!!   (D) istwfk(nkpt), the format of wavefunctions         (twfk)
!!   (E) nspinor, the scalar or spinor wf characteristics  (tspinor)
!! For PAW calculations:
!!   (F) the use of PAW method                             (tpaw)
!!   (G) the number of lmn elements for the paw basis      (tlmn)
!!   (H) the energy cutoff for the double (fine) grid      (tdg)
!! For WVL calculations:
!!   (I) the number of wavelets differs                    (twvl)
!!   (J) the space-grid size differs                       (tgrid)
!!
!!            non-self-consistent restarts
!!            ============================
!!
!! A restart will be allowed provided the following quantities in
!! old and new calculation are the same
!!
!!   (A) the primitive vectors                            (tprim)
!!   (B) the number of atoms of each type                 (tatty)
!!   (C) xred(3,natom)                                    (txred)
!!   (D) pseudopotentials (not just pseudocharges)        (tpseu)
!!   (E) the plane-wave cutoff                            (tecut)
!!   (F) ngfft(1:3)                                       (tng)
!! For PAW calculations:
!!   (G) the use of PAW method                            (tpaw)
!!   (H) the number of lmn elements for the paw basis     (tlmn)
!!   (I) the energy cutoff for the double (fine) grid     (tdg)
!!
!! PARENTS
!!      inwffil,ioarr,m_io_screening,m_wfk,setup_bse,setup_screening
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_check(fform,fform0,hdr,hdr0,mode_paral,restart,restartpaw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_check'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fform,fform0
 integer,intent(out) :: restart,restartpaw
 character(len=4),intent(in) :: mode_paral
 type(hdr_type),intent(in) :: hdr,hdr0

!Local variables-------------------------------
 character(len=1), parameter :: number(0:10)=(/'0','1','2','3','4','5','6','7','8','9',' '/)
 character(len=13), parameter :: filtypes(5)=(/'wf_planewave ','density      ','potential    ','screening    ', 'wf_wavelet   '/)
 character(len=24), save :: bndfmt='(2x, i4,t41,   a,2x, i4)'
 character(len=28), save :: occfmt='(2x, f4.1,t41,   a,2x, f4.1)'
 character(len=28), save :: wtkfmt='(2x, f7.3,t41,   a,2x, f7.3)'
 character(len=28), save :: zatfmt='(2x, f6.2,t41,   a,2x, f6.2)'
!scalars
 integer,parameter :: mwarning=5,nkpt_max=5
 integer :: bantot,bantot_eff,ii,ipsp,isppol,istart,istop,isym,itest,iwarning
 integer :: jj,mu,natom,nelm,nkpt,npsp,nsppol,nsym,ntypat,tatty,tband,tdg
 integer :: tecut,tgrid,tkpt,tlmn,tng,tpaw,tprim,tpsch,tpseu,tspinor,tsym,twfk
 integer :: twvl,txred
 real(dp) :: rms
 logical :: tfform2,tfform52
 character(len=13) :: filtyp,filtyp0
 character(len=26) :: typfmt
 character(len=500) :: msg

! *************************************************************************

 DBG_ENTER("COLL")

!We will adopt convention that if things agree between restart
!and current calculation then the tflag is 0. Begin by assuming
!that there is complete agreement between the files

 tatty = 0; tband = 0; tdg = 0 ; tecut = 0; tkpt = 0;
 tlmn = 0; tng = 0; tpaw = 0; tprim = 0; tpsch = 0; tpseu = 0;
 tspinor=0; tsym = 0; twfk = 0 ; txred = 0 ; twvl = 0 ; tgrid = 0

!Write out a header
 write(msg,'(a1,80a,2a1,10x,a,3a1,8x,a,25x,a,a1,8x,19a,25x,12a,a1)' )&
& ch10,('=',ii=1,80),ch10,ch10,&
& '- hdr_check: checking restart file header for consistency -',&
& (ch10,ii=1,3),'current calculation','restart file',ch10,('-',ii=1,19),('-',ii=1,12),ch10
 call wrtout(std_out,msg,mode_paral)

!Check validity of fform, and find filetype
 ii=1+fform /50
 if(ii>5)then
   if(fform==1002)then
     ii=4
   else
     write(msg, '(a,i0)' )'  Incorrect file format, fform=',fform
     MSG_BUG(msg)
   end if
 end if
 filtyp=filtypes(ii)

!Check validity of fform0, and find filetype
 ii=1+fform0/50
 if(ii>5)then
   if(fform0==1002)then
     ii=4
   else
     write(msg,'(a,i0,a)')&
&     'Incorrect file format, fform0=',fform0,'Action: it seems that the file you try to read is not an appropriate file.'
     MSG_ERROR(msg)
   end if
 end if
 filtyp0=filtypes(ii)

 write(msg,'(a,a,3x,3a)') &
& '  calculation expects a ',filtyp,'|','  input file contains a ',filtyp0
 call wrtout(std_out,msg,mode_paral)

 write(msg,'(a,a,11x,a,a,a)')&
& '. ABINIT  code version ',hdr%codvsn,'|','  ABINIT  code version ',hdr0%codvsn
 call wrtout(std_out,msg,mode_paral)

!Check fform from input, not from header file
 if ( (fform+1)/2 /= (fform0+1)/2  ) then
   write(msg,'(a,i0,a,i0,a)')&
&   'input fform=',fform,' differs from disk file fform=',fform0,'.'
   MSG_ERROR(msg)
 end if

 write(msg, '(a,i8,a,i4,a,i4,2x,a,a,i8,a,i4,a,i4)' ) &
& '. date ',hdr %date,' bantot ',hdr %bantot,' natom ',hdr %natom,'|',&
& '  date ',hdr0%date,' bantot ',hdr0%bantot,' natom ',hdr0%natom
 call wrtout(std_out,msg,mode_paral)

 write(msg, '(a,i4,a,i3,3(a,i4),2x,a,a,i4,a,i3,3(a,i4))' )&
& '  nkpt',hdr %nkpt,' nsym',hdr %nsym,' ngfft',hdr %ngfft(1),',',hdr %ngfft(2),',',hdr %ngfft(3),'|',&
& '  nkpt',hdr0%nkpt,' nsym',hdr0%nsym,' ngfft',hdr0%ngfft(1),',',hdr0%ngfft(2),',',hdr0%ngfft(3)
 call wrtout(std_out,msg,mode_paral)

 if (hdr%usewvl == 0) then
!  Note that the header actually contains ecut_eff=ecut*dilatmx**2
   write(msg,'(a,i3,a,f12.7,8x,a,a,i3,a,f12.7)')&
&   '  ntypat',hdr %ntypat,' ecut_eff',hdr %ecut_eff,'|',&
&   '  ntypat',hdr0%ntypat,' ecut_eff',hdr0%ecut_eff
   call wrtout(std_out,msg,mode_paral)
 else
   write(msg,'(a,i3,a,f12.7,8x,a,a,i3,a,f12.7)')&
&   '  ntypat',hdr %ntypat,' hgrid   ', 2. * hdr %rprimd(1,1) / (hdr %ngfft(1) - 31),'|',&
&   '  ntypat',hdr0%ntypat,' hgrid   ', 2. * hdr0%rprimd(1,1) / (hdr0%ngfft(1) - 31)
   call wrtout(std_out,msg,mode_paral)
!  Check hgrid and rprimd values.
   if (hdr0%rprimd(1,2) /= zero .or. hdr0%rprimd(1,3) /= zero .or. &
&   hdr0%rprimd(2,1) /= zero .or. hdr0%rprimd(2,3) /= zero .or. &
&   hdr0%rprimd(3,1) /= zero .or. hdr0%rprimd(3,2) /= zero) then
     msg = 'disk file rprimd is not parallelepipedic.'
     MSG_ERROR(msg)
   end if
   if (abs(hdr0%rprimd(1,1) / hdr0%ngfft(1) - hdr %rprimd(1,1) / hdr %ngfft(1)) > tol8) then
     write(msg,'(a,F7.4,a,F7.4)')&
&     'input wvl_hgrid=', 2. * hdr%rprimd(1,1) / hdr%ngfft(1), &
&     'not equal disk file wvl_hgrid=', 2. * hdr0%rprimd(1,1) / hdr0%ngfft(1)
     MSG_WARNING(msg)
     tgrid = 1
   end if
 end if

 write(msg, '(a,i3,29x,a,a,i3)' )&
& '  usepaw',hdr %usepaw,'|','  usepaw',hdr0%usepaw
 call wrtout(std_out,msg,mode_paral)

 write(msg, '(a,i3,29x,a,a,i3)' )&
& '  usewvl',hdr %usewvl,'|','  usewvl',hdr0%usewvl
 call wrtout(std_out,msg,mode_paral)

 write(msg,'(a,31x,a,a,3(a1,2x,3f12.7,2x,a,2x,3f12.7))')&
& '  rprimd:','|','  rprimd:',ch10,&
& hdr%rprimd(:,1),'|',hdr0%rprimd(:,1),ch10,&
& hdr%rprimd(:,2),'|',hdr0%rprimd(:,2),ch10,&
& hdr%rprimd(:,3),'|',hdr0%rprimd(:,3)
 call wrtout(std_out,msg,mode_paral)

 if (hdr%bantot/=hdr0%bantot) tband=1

 if (hdr%intxc/=hdr0%intxc) then
   write(msg,'(a,i0,a,i0)')'input intxc=',hdr%intxc,' not equal disk file intxc=',hdr0%intxc
   MSG_WARNING(msg)
 end if

 if (hdr%ixc/=hdr0%ixc) then
   write(msg,'(a,i0,a,i0)')'input ixc=',hdr%ixc,' not equal disk file ixc=',hdr0%ixc
   MSG_WARNING(msg)
 end if

 if (hdr%natom/=hdr0%natom) then
   write(msg,'(a,i0,a,i0)')'input natom=',hdr%natom,' not equal disk file natom=',hdr0%natom
   MSG_WARNING(msg)
   tatty=1
 end if

 if ( ANY(hdr%ngfft/=hdr0%ngfft) ) then
!  For sensible rho(r) or V(r) data, fft grid must be identical
!  MG TODO one should perform an FFT interpolation when the two ngfft differ!
   if (fform==52.or.fform==102) then
     write(msg, '(a,a,a,a,a)' )&
&     'fft grids must be the same for restart from a ',trim(filtyp),' file.',ch10,&
&     'Action: change your fft grid or your restart file.'
     MSG_ERROR(msg)
   end if
   tng=1
 end if

 if (hdr%nkpt/=hdr0%nkpt) then
   if (fform==2) then
     write(msg,'(a,i0,a,i0)' )'input nkpt=',hdr%nkpt,' not equal disk file nkpt=',hdr0%nkpt
     MSG_WARNING(msg)
   end if
   tkpt=1
   twfk=1
 end if

 if (hdr%nspinor/=hdr0%nspinor) then
   if (fform==2) then
     write(msg,'(a,i0,a,i0)')'input nspinor=',hdr%nspinor,' not equal disk file nspinor=',hdr0%nspinor
     MSG_WARNING(msg)
   end if
   tspinor=1
 end if

!No check is present for nspden

 if (hdr%nsppol/=hdr0%nsppol) then
   write(msg,'(a,i0,a,i0)')'input nsppol=',hdr%nsppol,'not equal disk file nsppol=',hdr0%nsppol
   MSG_WARNING(msg)
 end if

 if (hdr%nsym/=hdr0%nsym) then
   write(msg, '(a,i0,a,i0)' )'input nsym=',hdr%nsym,' not equal disk file nsym=',hdr0%nsym
   MSG_WARNING(msg)
   tsym=1
 end if

 if (hdr%ntypat/=hdr0%ntypat) then
   write(msg,'(a,i0,a,i0)')'input ntypat=',hdr%ntypat,' not equal disk file ntypat=',hdr0%ntypat
   call wrtout(std_out,msg,mode_paral)
   MSG_WARNING(msg)
   tatty=1
 end if

 if (hdr%usepaw/=hdr0%usepaw) then
   write(msg,'(a,i0,a,i0)')'input usepaw=',hdr%usepaw,' not equal disk file usepaw=',hdr0%usepaw
   MSG_WARNING(msg)
   tpaw=1
 end if

 if (hdr%usewvl/=hdr0%usewvl) then
   write(msg, '(a,i6,a,i6,a,a)' )&
&   'input usewvl=',hdr%usewvl,' not equal disk file usewvl=',hdr0%usewvl, ch10, &
&   'Action: change usewvl input variable or your restart file.'
   MSG_ERROR(msg)
 end if

!Also examine agreement of floating point data

 if (hdr%usewvl == 0 .and. abs(hdr%ecut_eff-hdr0%ecut_eff)>tol8) then
   write(msg,'(a,f12.6,a,f12.6,a)')'input ecut_eff=',hdr%ecut_eff,' /= disk file ecut_eff=',hdr0%ecut_eff,'.'
   MSG_WARNING(msg)
   tecut=1
 end if

 do ii=1,3
   do jj=1,3
     if (abs(hdr%rprimd(ii,jj)-hdr0%rprimd(ii,jj))>tol6) then
       write(msg, '(a,i1,a,i1,a,1p,e17.9,a,i1,a,i1,a,e17.9)' )&
&       'input rprimd(',ii,',',jj,')=',hdr%rprimd(ii,jj),' /= disk file rprimd(',ii,',',jj,')=',hdr0%rprimd(ii,jj)
       MSG_WARNING(msg)
       tprim=1
     end if
   end do
 end do

!Below this point many comparisons only make sense if
!certain things agree, e.g. nkpt, natom.  Also have to
!accomodate different amounts of data in general.

 if (hdr%usepaw==1 .and. hdr0%usepaw==1) then

!  Compare ecutdg (PAW)
   write(msg, '(a,f12.6,15x,a,a,f12.6)' )'  PAW: ecutdg',hdr %ecutdg,'|','  PAW: ecutdg',hdr0%ecutdg
   call wrtout(std_out,msg,mode_paral)
   if (hdr%ecutdg/=hdr0%ecutdg) then
     write(msg, '(a,f12.6,a,f12.6)' )'input ecutdg=',hdr%ecutdg,'not equal disk file ecutdg=',hdr0%ecutdg
     MSG_WARNING(msg)
     tdg=1
   end if
 end if

!Compare nband(nkpt*nsppol) (cannot compare if nkpt and nsppol not same)
 if (hdr%nkpt==hdr0%nkpt .and. hdr%nsppol==hdr0%nsppol) then
   nkpt=hdr%nkpt ; nsppol=hdr%nsppol
   write(msg,'(a,32x,a,a)') '  nband:','|','  nband:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,nsppol*nkpt,9
     istop = min(istart + 8,nsppol*nkpt)
     mu = istop - istart + 1
!    generate a format specifier
     bndfmt(5:5) = number(mu)
     bndfmt(21:21) = number(mu)
     write(msg,fmt=bndfmt) hdr%nband(istart:istop),'|',hdr0%nband(istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do

   do isppol=1,nsppol
     do ii=1,nkpt
       if (hdr%nband(ii)/=hdr0%nband(ii)) then
         tband=1
         if (fform == 2) then
           write(msg,'(a,i0,a,i0,a,i0)' )&
&           'kpt num',ii,' input nband=',hdr%nband(ii),'not equal disk file nband=',hdr0%nband(ii)
           MSG_WARNING(msg)
         end if
       end if
     end do
   end do
 end if

!Compare the number of wavelets in each resolution.
 if (hdr%usewvl == 1) then
   if (size(hdr%nwvlarr) /= size(hdr0%nwvlarr) .or. size(hdr%nwvlarr) /= 2) then
     write(msg, '(a,i6,a,i6,a,a)' )&
&     'input nwvlres=',size(hdr%nwvlarr),' not equal disk file nwvlres=',size(hdr0%nwvlarr),' or 2',&
&     ' ABINIT is not implemented for wavelet resolutions different from 2.'
     MSG_ERROR(msg)
   end if
!  This part is commented out since nwvlarr is not read from the
!  write(msg,'(a,30x,a,a)') '  nwvlres:','|', '  nwvlres:'
!  call wrtout(std_out,msg,mode_paral)
!  write(msg,'(a,2I6,24x,a,a,2I6)') '    ',hdr%nwvlarr, '|','    ', hdr0%nwvlarr
!  call wrtout(std_out,msg,mode_paral)
!  do ii=1,2
!  if (hdr%nwvlarr(ii)/=hdr0%nwvlarr(ii)) then
!  twvl=1
!  write(msg,'(a,a,a,a,i5,a,i8,a,i8)' ) ch10,&
!  &       ' hdr_check: WARNING -',ch10,&
!  &       '  nwvl resolution',ii,' input =',hdr%nwvlarr(ii),&
!  &       ' not equal disk file value=',hdr0%nwvlarr(ii)
!  call wrtout(std_out,msg,mode_paral)
!  end if
!  end do
 end if

!Compare symmetry arrays (integers) symafm(nsym)
!-- only for same number of symmetries nsym
 itest=0
 if (hdr%nsym==hdr0%nsym) then
   nsym=hdr%nsym
   write(msg,'(a,31x,a,a)') '  symafm:','|','  symafm:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,nsym,12
     istop=min(istart+11,nsym)
     nelm = istop - istart + 1
     call mk_hdr_check_fmt(nelm,typfmt)
     write(msg,fmt=typfmt) hdr%symafm(istart:istop),'|',hdr0%symafm(istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do
 end if

 if (itest/=0) then
   write(msg,'(a,i0,a)' )'For symmetry number',itest,' input symafm not equal disk file symafm'
   MSG_WARNING(msg)
   tsym=1
 end if

!Compare symmetry arrays (integers) symrel(3,3,nsym)
!-- only for same number of symmetries nsym
 itest=0
 if (hdr%nsym==hdr0%nsym) then
   nsym=hdr%nsym
   write(msg,'(a,31x,a,a)') '  symrel:','|','  symrel:'
   call wrtout(std_out,msg,mode_paral)
   do isym=1,nsym
     write(msg,'(2x,9i3,11x,a,2x,9i3)')hdr%symrel(:,:,isym),'|',hdr0%symrel(:,:,isym)
     call wrtout(std_out,msg,mode_paral)
     if(sum(abs(hdr%symrel(:,:,isym)-hdr0%symrel(:,:,isym)))/=0)then
       itest=isym
       exit
     end if
   end do
 end if

 if (itest/=0) then
   write(msg,'(a,i0,a)')'For symmetry number',itest,' input symrel not equal disk file symrel'
   MSG_WARNING(msg)
   tsym=1
 end if

!Compare typat(natom)
 if (hdr%natom==hdr0%natom) then
   natom=hdr%natom
   write(msg,'(a,32x,a,a)') '  typat:','|','  typat:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,natom,12
     istop=min(istart+11,natom)
     nelm = istop - istart + 1
     call mk_hdr_check_fmt(nelm,typfmt)
     write(msg,fmt=typfmt) hdr%typat(istart:istop),'|',hdr0%typat(istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do
   do ii=1,natom
     if (hdr%typat(ii)/=hdr0%typat(ii)) then
       write(msg, '(a,i0,a,i0,a,i0)' )&
&       'For atom number',ii,' input typat=',hdr%typat(ii),'not equal disk file typat=',hdr0%typat(ii)
       MSG_WARNING(msg)
       tatty=1
     end if
   end do
 end if


!Compare so_psp(npsp)
 if (hdr%npsp==hdr0%npsp) then
   npsp=hdr%npsp
   write(msg,'(a,29x,a,a)') '  so_psp  :','|','  so_psp  :'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,npsp  ,12
     istop=min(istart+11,npsp  )
     nelm = istop - istart + 1
     call mk_hdr_check_fmt(nelm,typfmt)
     write(msg,fmt=typfmt) hdr%so_psp  (istart:istop),'|',hdr0%so_psp  (istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do
   do ii=1,npsp
     if (hdr%so_psp  (ii)/=hdr0%so_psp  (ii)) then
       write(msg,'(a,i0,a,i0,a,i0)')&
&       'For pseudopotential number',ii,' input so_psp  =',hdr%so_psp(ii),'not equal disk file so_psp=',hdr0%so_psp(ii)
       MSG_WARNING(msg)
     end if
   end do
 end if


!Compare istwfk(nkpt)
 if (hdr%nkpt==hdr0%nkpt) then
   nkpt=hdr%nkpt
   write(msg,'(a,31x,a,a)') '  istwfk:','|','  istwfk:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,nkpt,12
     istop=min(istart+11,nkpt)
     nelm = istop - istart + 1
     call mk_hdr_check_fmt(nelm,typfmt)
     write(msg,fmt=typfmt) hdr%istwfk(istart:istop),'|',hdr0%istwfk(istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do
   do ii=1,nkpt
     if (hdr%istwfk(ii)/=hdr0%istwfk(ii)) then
       write(msg, '(a,i0,a,i0,a,i0)' )&
&       'For k point number',ii,' input istwfk=',hdr%istwfk(ii),'not equal disk file istwfk=',hdr0%istwfk(ii)
       MSG_WARNING(msg)
       twfk=1
     end if
   end do
 end if

!Compare kpt(3,nkpt)
 if (hdr%nkpt==hdr0%nkpt) then
   nkpt=hdr%nkpt
   write(msg,'(a,34x,a,a)') '  kpt:','|','  kpt:'
   call wrtout(std_out,msg,mode_paral)
   do ii = 1,min(nkpt,nkpt_max)
     write(msg,'(2x,3f12.7,2x,a,2x,3f12.7)')&
&     hdr%kptns(:,ii),'|',hdr0%kptns(:,ii)
     call wrtout(std_out,msg,mode_paral)
     if(ii>nkpt_max)then
       call wrtout(std_out,'The number of printed k points is sufficient... stop writing them.',mode_paral)
       exit
     end if
   end do
   iwarning=0
   do ii=1,nkpt
     itest=0
     do mu=1,3
       if(abs( hdr%kptns(mu,ii)-hdr0%kptns(mu,ii) )>tol6)itest=1
     end do
     if (itest==1) then
       write(msg, '(a,i5,a,3es17.7,a,a,3es17.7)' )&
&       'kpt num',ii,', input kpt=',hdr%kptns(:,ii),ch10,&
&       'not equal  disk file kpt=',hdr0%kptns(:,ii)
       MSG_WARNING(msg)
       tkpt=1 ; iwarning=iwarning+1
       if(iwarning>=mwarning)then
         call wrtout(std_out,'The number of warning messages is sufficient ... stop writing them.',mode_paral)
         exit
       end if
     end if
   end do
 end if

!Compare wtk(nkpt)
 if (hdr%nkpt==hdr0%nkpt) then
   nkpt=hdr%nkpt

   write(msg,'(a,34x,a,a)') '  wtk:','|','  wtk:'
   call wrtout(std_out,msg,mode_paral)
   istop = min(nkpt,nkpt_max)
   do ii = 1, istop, 5
     mu = min(5, istop - ii + 1)
     wtkfmt(5:5) = number(mu)
     wtkfmt(23:23) = number(mu)
     write(msg, wtkfmt)hdr%wtk(ii:min(istop, ii + 5 - 1)),'|',hdr0%wtk(ii:min(istop, ii + 5 - 1))
     call wrtout(std_out,msg,mode_paral)
   end do
   iwarning=0
   do ii=1,nkpt
     itest=0
     if (abs( hdr%wtk(ii)-hdr0%wtk(ii) )>tol6) then
       write(msg,'(a,i5,a,es17.7,a,a,es17.7)')&
&       'kpt num',ii,', input weight=',hdr%wtk(ii),ch10,&
&       'not equal  disk file weight=',hdr0%wtk(ii)
       MSG_WARNING(msg)

       tkpt=1 ; iwarning=iwarning+1
       if(iwarning>=mwarning)then
         call wrtout(std_out,'The number of warning messages is sufficient ... stop writing them.',mode_paral)
         exit
       end if
     end if
   end do
 end if

!Compare occ(bantot)
 if (hdr%nkpt==hdr0%nkpt.and. hdr%bantot==hdr0%bantot) then
   nkpt=hdr%nkpt
   bantot=hdr%bantot

   write(msg,'(a,34x,a,a)') '  occ:','|','  occ:'
   call wrtout(std_out,msg,mode_paral)
   bantot_eff=min(bantot,9*nkpt_max)
   do istart = 1,bantot_eff,9
     istop = min(istart+8,bantot_eff)
     mu = istop - istart + 1
     occfmt(5:5) = number(mu)
     occfmt(23:23) = number(mu)
     write(msg,fmt=occfmt)hdr%occ(istart:istop),'|', hdr0%occ(istart:istop)
     call wrtout(std_out,msg,mode_paral)
     if(istart>9*nkpt_max)then
       call wrtout(std_out,'The number of printed occupation numbers is sufficient ... stop writing them.',mode_paral)
       exit
     end if
   end do
   iwarning=0
   do ii=1,bantot
     if (abs( hdr%occ(ii)-hdr0%occ(ii) )>tol6) then
       write(msg,'(a,i10,a,1p,e15.7,a,e15.7)')'band,k',ii,', input occ=',hdr%occ(ii),' disk occ=',hdr0%occ(ii)
       MSG_WARNING(msg)
       tband=1 ; iwarning=iwarning+1
       if(iwarning>=mwarning)then
         call wrtout(std_out,'The number of warning msgs is sufficient ... stop writing them.',mode_paral)
         exit
       end if
     end if
   end do
 end if

!Compare tnons(3,nsym)
 if (hdr%nsym==hdr0%nsym) then
   nsym=hdr%nsym
   itest=0
   write(msg,'(a,32x,a,a)') '  tnons:','|','  tnons:'
   call wrtout(std_out,msg,mode_paral)
   do isym=1,nsym
     write(msg,'(2x,3f12.7,2x,a,2x,3f12.7)') hdr%tnons(:,isym),'|',hdr0%tnons(:,isym)
     call wrtout(std_out,msg,mode_paral)
   end do

   do isym=1,nsym
     if( sum(abs(  hdr%tnons(:,isym)-hdr0%tnons(:,isym) )) > tol6) then
       itest=isym
       exit
     end if
   end do
   if (itest/=0) then
     write(msg, '(a,i0,a)' )'For symmetry number',itest,' input tnons not equal disk file tnons'
     MSG_WARNING(msg)
   end if
 end if

!Compare znucltypat(ntypat)
 if (hdr%ntypat==hdr0%ntypat) then
   ntypat=hdr%ntypat

   write(msg,'(a,31x,a,a)') '   znucl:','|','   znucl:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,ntypat,6
     istop = min(istart+5,ntypat)
     mu = istop-istart+1
     zatfmt(5:5) = number(mu)
     zatfmt(23:23) = number(mu)
     write(msg,fmt=zatfmt) hdr%znucltypat(istart:istop),'|',hdr0%znucltypat(istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do

   do ii=1,ntypat
     if (abs(hdr%znucltypat(ii)-hdr0%znucltypat(ii))>tol6) then
       write(msg, '(a,i5,a,f12.6,a,f12.6)' )&
&       ' For atom number',ii,' input znucl=',hdr%znucltypat(ii),' not equal disk file znucl=',hdr0%znucltypat(ii)
       MSG_WARNING(msg)
     end if
   end do
 end if

!Should perform some checks related to pertcase and qptn,
!that have been introduced in the header in v4.1
!Warning : a GS file might be read, while the hdr corresponds
!to a RF file (to initialize k+q), and vice-versa (in nonlinear).

!Now check agreement of psp headers too
 if (hdr%npsp==hdr0%npsp) then
   npsp=hdr%npsp
   itest=0

   do ipsp=1,npsp

     write(msg,'(a,i3,a,9x,a,a,i3,a)')&
&     '  pseudopotential atom type',ipsp,':','|','  pseudopotential atom type',ipsp,':'
     call wrtout(std_out,msg,mode_paral)

     if (hdr%usepaw==1 .and. hdr0%usepaw==1) then
       write(msg,'(a,i3,a,i3,a,i3,5x,a,a,i3,a,i3,a,i3)')&
&       '  pspso ',hdr %pspso(ipsp),' pspxc ',hdr %pspxc(ipsp),&
&       '  lmn_size ',hdr%lmn_size(ipsp),'|',&
&       '  pspso ',hdr0%pspso(ipsp),' pspxc ',hdr0%pspxc(ipsp),&
&       '  lmn_size ',hdr0%lmn_size(ipsp)
       call wrtout(std_out,msg,mode_paral)
       if (hdr%lmn_size(ipsp)/=hdr0%lmn_size(ipsp)) then
         write(msg, '(a,i3,a,i3,a,i3)' )&
&         '  For atom type ',ipsp,' input lmn_size=',hdr%lmn_size(ipsp),&
&         ' not equal disk file lmn_size=',hdr0%lmn_size(ipsp)
         MSG_WARNING(msg)
         tlmn=1
       end if
     else
       write(msg,'(a,i3,a,i3,19x,a,a,i3,a,i3)')&
&       '  pspso ',hdr %pspso(ipsp),' pspxc ',hdr %pspxc(ipsp),'|',&
&       '  pspso ',hdr0%pspso(ipsp),' pspxc ',hdr0%pspxc(ipsp)
       call wrtout(std_out,msg,mode_paral)
     end if
     write(msg,'(a,i6,a,i4,a,f5.1,2x,a,a,i6,a,i4,a,f5.1)')&
&     '  pspdat ',hdr %pspdat(ipsp),' pspcod ',hdr %pspcod(ipsp),&
&     ' zion ',hdr %zionpsp(ipsp),'|',&
&     '  pspdat ',hdr0%pspdat(ipsp),' pspcod ',hdr0%pspcod(ipsp),&
&     ' zion ',hdr0%zionpsp(ipsp)
     call wrtout(std_out,msg,mode_paral)

!    Second, test
!    NOTE, XG 000719 : should do something about pspso
!    NOTE, XG 020716 : znucl and zion are not written
     if (abs(hdr%znuclpsp(ipsp)-hdr0%znuclpsp(ipsp))>tol6) itest=1
     if (abs(hdr%zionpsp(ipsp)-hdr0%zionpsp(ipsp))>tol6) then
       itest=1
       tpsch=1
     end if
     if (hdr%pspdat(ipsp)/= hdr0%pspdat(ipsp)) itest=1
     if (hdr%pspcod(ipsp)/= hdr0%pspcod(ipsp)) itest=1
     if (hdr%pspxc(ipsp) /= hdr0%pspxc(ipsp) )  itest=1
   end do

   if (itest==1) then
     msg = 'input psp header does not agree perfectly with disk file psp header.'
     MSG_WARNING(msg)
     tpseu=1
   end if
 end if

!Finally, read residm and etotal ("current value" not known), and check xred.

 if (hdr%natom==hdr0%natom) then

   natom=hdr%natom
   write(msg,'(a,33x,a,a)') '  xred:','|','  xred:'
   call wrtout(std_out,msg,mode_paral)
   do ii=1,natom
     write(msg,'(2x,3f12.7,2x,a,2x,3f12.7)') hdr%xred(:,ii),'|',hdr0%xred(:,ii)
     call wrtout(std_out,msg,mode_paral)
   end do

!  check atom positions one atom at a time and allow possibility
!  that there is a harmless translation of atoms by a cell vector.
   do ii=1,natom
     rms=0.0_dp
     do jj=1,3
       rms=rms+(hdr%xred(jj,ii)-hdr0%xred(jj,ii) - dble(nint((hdr%xred(jj,ii)-hdr0%xred(jj,ii)))) )**2
     end do
     rms=sqrt(rms/3.0_dp)
     if (rms>tol6) txred=1
   end do
 end if

!Run tests here to establish whether this is a valid restart

!tfform2 will be true if there is a problem for the wavefunctions
 tfform2 = (hdr%usewvl == 0 .and. &
& (tprim /= 0 .or. tecut /= 0 .or. tkpt /= 0 .or. &
& twfk /=0 .or. tspinor /= 0)) .or. &
& (hdr%usepaw == 1 .and. &
& (tpaw /= 0 .or. tlmn /= 0 .or. tdg /= 0)) .or. &
& (hdr%usewvl == 1 .and. &
& (tatty /= 0 .or. tband /= 0))
!tfform52 will be true if there is a problem for the format 52
 tfform52=tprim /= 0 .or. tatty /= 0 .or. txred /= 0 .or.&
& tpseu /= 0 .or. tecut /= 0 .or. tng /= 0 .or. &
& (hdr%usepaw == 1 .and. &
& (tpaw /= 0 .or. tlmn /= 0 .or. tdg /= 0))

 restart=1
 restartpaw=hdr%usepaw

!If there is a problem somewhere
 if ( (fform == 2  .and. tfform2  ) .or.  &
& (fform == 52 .and. tfform52 ) .or.  &
& (fform == 200 .and. tfform2 ) ) then

   if(fform==2)then
     restart=2
     msg = 'Restart of self-consistent calculation need translated wavefunctions.'
   else if(fform==52)then
     restart=0
     msg = 'Illegal restart of non-self-consistent calculation'
   end if
   MSG_WARNING(msg)

   write(msg,'(a,a1,a)') &
&   '  Indeed, critical differences between current calculation and',ch10,&
&   '  restart file have been detected in:'
   call wrtout(std_out,msg,mode_paral)

   if ( (fform==52 .or. fform == 200) .and. tatty /= 0 ) then
     write(msg, '(8x,a)' ) '* the number of atoms of each type'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( fform /= 200 .and. tecut /= 0 ) then
     write(msg, '(8x,a)' ) '* the plane-wave cutoff'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( fform == 200 .and. tband /= 0 ) then
     write(msg, '(8x,a)' ) '* the band and their occupation'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( fform==2 .and. tkpt /= 0 ) then
     write(msg, '(8x,a)' ) '* the number, position, or weight of k-points'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( fform==2 .and. twfk /= 0 ) then
     write(msg, '(8x,a)' ) '* the format of wavefunctions (istwfk)'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( fform==2 .and. tspinor /= 0 ) then
     write(msg, '(8x,a)' ) '* the scalar/spinor character of the wf (nspinor)'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( fform==52 .and. tng /= 0 ) then
     write(msg, '(8x,a)' ) '* the Fourier transform box dimensions'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( tprim /= 0 ) then
     write(msg, '(8x,a)' )'* the vectors defining the unit cell (obtained from rprim and acell)'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( fform==52 .and. tpseu /= 0 ) then
     write(msg, '(8x,a)' )'* the pseudopotential files'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( fform==52 .and. txred /= 0 ) then
     write(msg, '(8x,a)' ) '* the positions of the ions in the basis'
     call wrtout(std_out,msg,mode_paral)
   end if

!  Tests for a restart in the framework of the PAW method
   if (hdr%usepaw/=0 .or. hdr0%usepaw/=0) then
     if (tpaw /= 0 .or. tlmn /= 0) restartpaw=0
     if (restartpaw == 0) then
       write(msg,'(8x,a)') 'Critical differences for a restart within PAW method:'
       call wrtout(std_out,msg,mode_paral)
       if ( tpaw /= 0 ) then
         write(msg, '(8x,a)' ) '* the use of the PAW method'
         call wrtout(std_out,msg,mode_paral)
       else
         if(tlmn/=0)then
           write(msg, '(8x,a)' ) '* the number of lmn elements for the paw basis'
           call wrtout(std_out,msg,mode_paral)
         end if
       end if
     else if (tdg/=0) then
       write(msg,'(a,a,a,a,a,a)') ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  Restart of calculation within PAW may be inconsistent because of:"'
       call wrtout(std_out,msg,mode_paral)
       if(tdg/=0)then
         write(msg, '(8x,a)' )'* the cutoff energy of the paw double (fine) grid'
         call wrtout(std_out,msg,mode_paral)
       end if
     end if
   end if

 else

   if(fform==2 .or. fform == 200)then
     write(msg,'(a,a)') ' hdr_check: ',' Wavefunction file is OK for direct restart of calculation'
     call wrtout(std_out,msg,mode_paral)
   else if(fform==52)then
     write(msg,'(a,a)') ' hdr_check: ',' Density/Potential file is OK for restart of calculation'
     call wrtout(std_out,msg,mode_paral)
   end if
!  MG TODO add screening case! 
 end if

 write(msg,'(80a)') ('=',ii=1,80)
 call wrtout(std_out,msg,mode_paral)

 CONTAINS
!!***

!!****f* hdr_check/mk_hdr_check_fmt
!! NAME
!! mk_hdr_check_fmt
!!
!! FUNCTION
!! make a format needed in hdr_check, for arrays of nint integers each of format i3
!!
!! INPUTS
!!  nelm=number of elements to be printed
!!
!! OUTPUT
!!  character(len=26), typfmt= format needed
!!
!! PARENTS
!!      m_header
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine mk_hdr_check_fmt(nelm,typfmt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mk_hdr_check_fmt'
!End of the abilint section

   implicit none

!  Arguments ------------------------------------
!  scalars
   integer,intent(in) :: nelm
   character(len=26),intent(out) :: typfmt

!  Local variables-------------------------------
!  scalars
   integer :: ii
   character(len=1), parameter :: number(0:10)=(/'0','1','2','3','4','5','6','7','8','9',' '/)
   character(len=26), parameter :: templatefmt='(2x,  i3,t41   ,a,2x,  i3)'
!  *************************************************************************

!  Initialize the format
   typfmt=templatefmt

!  Generate the type format specifier
   ii=nelm/10
   if ( ii /= 0 ) then
     typfmt(5:5) = number(ii)
     typfmt(22:22) = number(ii)
   else
     typfmt(5:5) = ' '
     typfmt(22:22) = ' '
   end if
   ii = nelm - 10 * (nelm/10)
   typfmt(6:6) = number(ii)
   typfmt(23:23) = number(ii)
 end subroutine mk_hdr_check_fmt

end subroutine hdr_check
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_fort_read
!! NAME
!! hdr_fort_read
!!
!! FUNCTION
!! Reads the header from a logical unit associated to a unformatted file.
!! Note that, when reading, different records of hdr are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated correctly by a call to hdr_free when hdr is not used anymore.
!!
!! INPUTS
!!  unit=unit number of the unformatted file 
!!  [rewind]=True to rewind the file. Default: False
!!
!! OUTPUT
!!  Hdr<hdr_type>=The header of the file fully initialized (if fform /=0)
!!  fform=kind of the array in the file.  if the reading fail, return fform=0
!!
!! NOTES
!! The file is supposed to be open already
!!
!! PARENTS
!!      m_header
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_fort_read(Hdr,unit,fform,rewind)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_fort_read'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(out)  :: fform
 integer,intent(in) :: unit
 logical,optional,intent(in) :: rewind
 type(hdr_type),intent(out) :: hdr

!Local variables-------------------------------
 integer :: bantot,headform,ierr,ipsp,lloc,lmax,mmax,natom,nkpt,npsp,nsppol,nsym,ntypat
 character(len=500) :: msg
 character(len=6) :: codvsn
 real(dp) :: acell(3)

!*************************************************************************

 DBG_ENTER("COLL")

!-------------------------------------------------------------------------
!Reading the header of an unformatted file
!-------------------------------------------------------------------------
 if (present(rewind)) then 
   if (rewind) rewind(unit)
 end if

 ! Reading the first record of the file ------------------------------------
 read(unit,iostat=ierr)codvsn,fform
 if (ierr /=0) then
   fform=0
   return ! This is to allow treatment of old epsm1 format
 end if

 if(fform==1   .or. &
& fform==2   .or. &
& fform==51  .or. &
& fform==52  .or. &
& fform==101 .or. &
& fform==102       )then
!  This is the old format
   headform=22

 else
   ! Format beyond 22 have a different first line, so need reading again the first line
   backspace (unit)
   read (unit) codvsn,headform,fform

   if(headform/=23 .and. &
&   headform/=34 .and. &
&   headform/=40 .and. &
&   headform/=41 .and. &
&   headform/=42 .and. &
&   headform/=44 .and. &
&   headform/=53 .and. &
&   headform/=56 .and. &
&   headform/=57         )then
     write(msg,'(a,i0,3a,i0,3a)')&
&     'The first line of the (WF, DEN or POT) file read in unit ',unit,' is erroneous.',ch10,&
&     'headform is ',headform,', while it should be 23, 34, 40, 41, 42, 44, 53 or 56 or 57.',ch10,&
&     'Action : check the correctness of your file.'
     MSG_ERROR(msg)
   end if
 end if

 hdr%codvsn=codvsn
 hdr%headform=headform
!fform is not a record of hdr_type

!Reading the second record of the file ------------------------------------

!Initialize the values that are not present for all versions (exception : npsp)
 hdr%nspden=1
 hdr%nspinor=1
 hdr%occopt=1
 hdr%pertcase=1
 hdr%usepaw=0
 hdr%usewvl=0
 hdr%ecut=zero
 hdr%ecutdg=zero
 hdr%ecutsm=zero
 hdr%qptn(1:3)=zero
 hdr%stmbias=zero
 hdr%tphysel=zero
 hdr%tsmear=zero

 if(headform==22)then

   read(unit,err=10) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&   nkpt, nsppol, nsym, ntypat,acell, hdr%ecut_eff, hdr%rprimd
   npsp=ntypat

 else if(headform==23)then

!  Compared to v2.2, add nspden, nspinor, occopt
   read(unit,err=10) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&   nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, ntypat, hdr%occopt,&
&   acell, hdr%ecut_eff, hdr%rprimd
   npsp=ntypat

 else if(headform==34)then

!  Compared to v2.3, subtract acell, and add npsp
   read(unit,err=10) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&   nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt,&
&   hdr%ecut_eff, hdr%rprimd

 else if(headform==40)then

!  Compared to v3.4, add ecut, ecutsm, tphysel, tsmear
   read(unit,err=10) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&   nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt,&
&   hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%rprimd, hdr%tphysel, hdr%tsmear

 else if(headform==41)then

!  Compared to v4.0, add pertcase and qptn(3)
   read(unit,err=10) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&   nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
&   hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd, hdr%tphysel, hdr%tsmear

 else if(headform==42)then

!  Compared to v4.1, add stmbias
   read(unit,err=10) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&   nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
&   hdr%ecut, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
&   hdr%stmbias, hdr%tphysel, hdr%tsmear

 else if(headform>=44 .and. headform<57)then

!  Compared to v4.2, add usepaw and ecutdg
   read(unit,err=10) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&   nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
&   hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
&   hdr%stmbias, hdr%tphysel, hdr%tsmear

 else if(headform>=57)then

!  Compared to v4.4, add usewvl
   read(unit,err=10) bantot, hdr%date, hdr%intxc, hdr%ixc, natom, hdr%ngfft(1:3),&
&   nkpt, hdr%nspden, hdr%nspinor, nsppol, nsym, npsp, ntypat, hdr%occopt, hdr%pertcase,&
&   hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
&   hdr%stmbias, hdr%tphysel, hdr%tsmear, hdr%usewvl
 end if

 hdr%bantot=bantot
 hdr%natom =natom
 hdr%nkpt  =nkpt
 hdr%npsp  =npsp
 hdr%nsppol=nsppol
 hdr%nsym  =nsym
 hdr%ntypat =ntypat

 if(hdr%ecutsm>tol6 .and. headform<44 .and. .not.(fform==51.or.fform==52.or.fform==101.or.fform==102))then
   write(msg,'(a,es16.6,9a)')&
&   'The value of ecutsm is',hdr%ecutsm,', while the file has been produced prior to v4.4 .',ch10,&
&   'The definition of the smearing function has changed, so that you are not allowed',ch10,&
&   'to restart from a old wavefunction file. By contrast, you can restart from an old',ch10,&
&   'potential or density file, and perform a self-consistent cycle with a new ABINIT version.',ch10,&
&   'Action: produce a density or potential file using the old version of ABINIT, and restart from it.'
   MSG_ERROR(msg)
 end if

!Allocate all parts of hdr that need to be --------------------------------
 ABI_MALLOC(hdr%istwfk,(nkpt))
 ABI_MALLOC(hdr%kptns,(3,nkpt))
 ABI_MALLOC(hdr%lmn_size,(npsp))
 ABI_MALLOC(hdr%nband,(nkpt*nsppol))
 ABI_MALLOC(hdr%npwarr,(nkpt))
 ABI_MALLOC(hdr%occ,(bantot))
 ABI_MALLOC(hdr%pspcod,(npsp))
 ABI_MALLOC(hdr%pspdat,(npsp))
 ABI_MALLOC(hdr%pspso,(npsp))
 ABI_MALLOC(hdr%pspxc,(npsp))
 ABI_MALLOC(hdr%so_psp,(npsp))
 ABI_MALLOC(hdr%symafm,(nsym))
 ABI_MALLOC(hdr%symrel,(3,3,nsym))
 ABI_MALLOC(hdr%title,(npsp))
 ABI_MALLOC(hdr%tnons,(3,nsym))
 ABI_MALLOC(hdr%typat,(natom))
 ABI_MALLOC(hdr%wtk,(nkpt))
 ABI_MALLOC(hdr%xred,(3,natom))
 ABI_MALLOC(hdr%zionpsp,(npsp))
 ABI_MALLOC(hdr%znuclpsp,(npsp))
 ABI_MALLOC(hdr%znucltypat,(ntypat))

 if(hdr%usepaw==1)  then
   ABI_DT_MALLOC(hdr%pawrhoij,(natom))
 end if

!Reading the third record of the file ------------------------------------

!Initialize the values that are not present for all versions
 hdr%istwfk(:)=1
 hdr%so_psp(:)=1
 hdr%symafm(:)=1

 if(headform==22 .and. (fform==1 .or. fform==51 .or. fform==101))then

!  This is very old (pre-2.0) format !
   read(unit,err=10) hdr%nband(:), hdr%npwarr(:), hdr%symrel(:,:,:), &
&   hdr%typat(:), hdr%kptns(:,:), hdr%occ(:),hdr%tnons(:,:), hdr%znucltypat(:)

 else if(headform==22 .or. headform==23 .or. headform==34)then

!  Compared to pre v2.0, add istwfk
   read(unit,err=10) hdr%nband(:), hdr%npwarr(:), hdr%symrel(:,:,:), &
&   hdr%typat(:), hdr%istwfk(:), hdr%kptns(:,:), hdr%occ(:), &
&   hdr%tnons(:,:), hdr%znucltypat(:)

 else if(headform>=40 .and. headform < 50)then

!  Compared to pre v4.0, add so_psp and symafm, and switch istwfk

   read(unit,err=10)  hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
&   hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
&   hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
&   hdr%tnons(:,:), hdr%znucltypat(:)

 else if(headform>=50)then

!  Compared to pre v5.0, add wtk
   read(unit,err=10)  hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
&   hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
&   hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
&   hdr%tnons(:,:), hdr%znucltypat(:), hdr%wtk(:)

 end if

!Reading the records with psp information ---------------------------------

!Initialize the values that are not present for all versions
 hdr%pspso(:)=1
 hdr%lmn_size(:)=0

 if(headform==22)then

   do ipsp=1,npsp
     read(unit,err=10) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
&     hdr%zionpsp(ipsp), hdr%pspdat(ipsp), hdr%pspcod(ipsp), &
&     hdr%pspxc(ipsp), lmax, lloc, mmax
   end do

 else if(headform==23)then

!  Compared to 2.2, add pspso
   do ipsp=1,npsp
     read(unit,err=10) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
&     hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
&     hdr%pspcod(ipsp), hdr%pspxc(ipsp), lmax, lloc, mmax
   end do

 else if(headform==34 .or. headform==40 .or. headform==41 .or. headform==42)then

!  Compared to 2.3, suppress lmax, lloc, mmax
   do ipsp=1,npsp
     read(unit,err=10) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
&     hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
&     hdr%pspcod(ipsp), hdr%pspxc(ipsp)
   end do

 else if(headform>=44)then

!  Compared to 4.2, add lmn_size
   do ipsp=1,npsp
     read(unit,err=10) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
&     hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
&     hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp)
   end do

 end if

!Reading the final record of the header  ---------------------------------

!Initialize the values that are not present for all versions
 hdr%fermie=zero

 if(headform==22)then
   read(unit,err=10) hdr%residm, hdr%xred(:,:), hdr%etot
 else if(headform==23 .or. headform==34 .or. headform>=40)then
   read(unit,err=10) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie
 end if

 if (hdr%usepaw==1) then ! Reading the Rhoij tab if the PAW method was used.
   call pawrhoij_io(hdr%pawrhoij,unit,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,headform,"Read")
 end if

 DBG_EXIT("COLL")
 return 

! Handle IO-error.
10 fform=0

end subroutine hdr_fort_read
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_ncread
!! NAME
!! hdr_ncread
!!
!! FUNCTION
!! This subroutine deals with the reading of the hdr_type structured variables 
!! It handles variables according to the ETSF format, whenever
!! possible and uses new variables when not available in the ETSF format.
!! Note that, when reading, different records of hdr are allocated here, 
!! Records of hdr should be deallocated
!! correctly by a call to hdr_free when hdr is not used anymore.
!!
!! INPUTS
!!  ncid=the unit of the open NetCDF file.
!!
!! OUTPUT
!!  fform=kind of the array in the file. if the reading fails, return fform=0
!!
!! PARENTS
!!      m_header
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_ncread(Hdr,ncid,fform)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_ncread'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 integer,intent(out) :: fform
 type(hdr_type),target,intent(out) :: hdr

#ifdef HAVE_TRIO_ETSF_IO
!Local variables-------------------------------
!scalars
 integer :: rhoijdim1,rhoijdim2,nresolution,headform,iatom,itypat
 !integer :: 
 real(dp),target :: ecut, fermie
 logical :: lstat
 character(len=500) :: msg
 type(etsf_dims) :: dims
 type(etsf_kpoints) :: kpoints
 type(etsf_basisdata) :: basisdata
 type(etsf_geometry) :: geometry
 type(etsf_electrons) :: electrons
 type(etsf_io_low_error) :: error_data
!arrays
 real(dp), target :: rprimd(3, 3)
! integer :: cplex, ilmn, irhoij, ispden, lmn2_size, nselect
! real(dp), allocatable :: rhoij(:,:,:)

! *************************************************************************

 !Switch on write mode mode.
 call etsf_io_low_set_write_mode(ncid, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

!In case the file is just ETSF valid, we ignore the missing variables and we use default values.
 hdr%codvsn   = "ETSF  "
 fform        = 1
 headform     = 57

!First, we read the declaration of code, fform ...
!We ignore errors, assuming that the file is at least ETSF valid.
 call etsf_io_low_read_var(ncid, "codvsn", hdr%codvsn, 6, lstat, error_data = error_data)
 if (lstat) then ! We pad the returned string with " " instead of "\0"
   call strip(hdr%codvsn)
 end if

 call etsf_io_low_read_var(ncid, "fform", fform, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_read_var(ncid, "headform", hdr%headform, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 if (headform <= 42) then
   write(msg,'(a,i0,3a)')&
&   'headform is ',headform,', while it should be > 42.',ch10,&
&   'Action: check the correctness of your file.'
   MSG_ERROR(msg)
 end if

!Then, read dimensions handled by ETSF
 call etsf_io_dims_get(ncid, dims, lstat, error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

!Copy dimensions to hdr structure
!FIXME: don't handle k_dependent = 1
 hdr%bantot   = dims%max_number_of_states * dims%number_of_kpoints * dims%number_of_spins
 hdr%natom    = dims%number_of_atoms
 hdr%nkpt     = dims%number_of_kpoints
 hdr%nspden   = dims%number_of_components
 hdr%nspinor  = dims%number_of_spinor_components
 hdr%nsppol   = dims%number_of_spins
 hdr%nsym     = dims%number_of_symmetry_operations
 hdr%ntypat   = dims%number_of_atom_species
 hdr%ngfft(1) = dims%number_of_grid_points_vector1
 hdr%ngfft(2) = dims%number_of_grid_points_vector2
 hdr%ngfft(3) = dims%number_of_grid_points_vector3

!We read other dimensions, not handled by ETSF format.
!In case the file is just ETSF valid, we ignore the missing dimensions and we use default values.

 hdr%npsp    = hdr%ntypat
 call etsf_io_low_read_dim(ncid, "npsp", hdr%npsp, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 rhoijdim1   = 1
 call etsf_io_low_read_dim(ncid, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 rhoijdim2   = hdr%nspden
 call etsf_io_low_read_dim(ncid, "rhoijdim2", rhoijdim2, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%usepaw  = 0
 call etsf_io_low_read_var(ncid, "usepaw", hdr%usepaw, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%usewvl  = 0
 call etsf_io_low_read_var(ncid, "usewvl", hdr%usewvl, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 nresolution=0
 if (hdr%usewvl == 1) then
!  This value must be 2...
   call etsf_io_low_read_dim(ncid, "number_of_wavelet_resolutions",nresolution, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

!  We set the right ngfft, adding the padding space for wavelets.
   hdr%ngfft = hdr%ngfft + 31
 end if

!Allocate all parts of hdr that need to be
 ABI_MALLOC(hdr%istwfk,(hdr%nkpt))
 ABI_MALLOC(hdr%lmn_size,(hdr%npsp))
 ABI_MALLOC(hdr%nband,(hdr%nkpt*hdr%nsppol))
 ABI_MALLOC(hdr%npwarr,(hdr%nkpt))
 ABI_MALLOC(hdr%pspcod,(hdr%npsp))
 ABI_MALLOC(hdr%pspdat,(hdr%npsp))
 ABI_MALLOC(hdr%pspso,(hdr%npsp))
 ABI_MALLOC(hdr%pspxc,(hdr%npsp))
 ABI_MALLOC(hdr%so_psp,(hdr%npsp))
 ABI_MALLOC(hdr%symafm,(hdr%nsym))
 ABI_MALLOC(hdr%symrel,(3,3,hdr%nsym))
 ABI_MALLOC(hdr%typat,(hdr%natom))
 ABI_MALLOC(hdr%kptns,(3,hdr%nkpt))
 ABI_MALLOC(hdr%occ,(hdr%bantot))
 ABI_MALLOC(hdr%tnons,(3,hdr%nsym))
 ABI_MALLOC(hdr%wtk,(hdr%nkpt))
 ABI_MALLOC(hdr%xred,(3,hdr%natom))
 ABI_MALLOC(hdr%znuclpsp,(hdr%npsp))
 ABI_MALLOC(hdr%znucltypat,(hdr%ntypat))
 ABI_MALLOC(hdr%zionpsp,(hdr%npsp))
 ABI_MALLOC(hdr%title,(hdr%npsp))
 if(hdr%usepaw==1)  then
   ABI_DT_MALLOC(hdr%pawrhoij,(hdr%natom))
 end if

!We get then all variables included in ETSF
 if (hdr%usewvl==0) then
   basisdata%kinetic_energy_cutoff => ecut
   basisdata%number_of_coefficients => hdr%npwarr
   call etsf_io_basisdata_get(ncid, basisdata, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)
 else
   call etsf_io_low_read_var(ncid, "number_of_wavelets", hdr%nwvlarr, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)
 end if

 electrons%fermi_energy => fermie
 electrons%number_of_states%data1D => hdr%nband
 electrons%occupations%data1D => hdr%occ

 call etsf_io_electrons_get(ncid, electrons, lstat, error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 geometry%primitive_vectors => rprimd
 geometry%reduced_symmetry_matrices => hdr%symrel
 geometry%atom_species => hdr%typat
 geometry%reduced_symmetry_translations => hdr%tnons
 geometry%reduced_atom_positions => hdr%xred
 geometry%atomic_numbers => hdr%znucltypat

 if (hdr%npsp == hdr%ntypat) then
   geometry%valence_charges => hdr%zionpsp
 end if

 call etsf_io_geometry_get(ncid, geometry, lstat, error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 kpoints%reduced_coordinates_of_kpoints => hdr%kptns
 kpoints%kpoint_weights => hdr%wtk

 call etsf_io_kpoints_get(ncid, kpoints, lstat, error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 hdr%fermie = fermie
 hdr%ecut   = ecut
 hdr%rprimd = rprimd
 hdr%znuclpsp(1:hdr%npsp) = hdr%znucltypat(1:hdr%npsp)

!We get all other variables
!In case the file is just ETSF valid, we ignore the missing variables and we use default values.

 hdr%date = 0
 call etsf_io_low_read_var(ncid, "date", hdr%date, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%ecut_eff = hdr%ecut
 call etsf_io_low_read_var(ncid, "ecut_eff", hdr%ecut_eff, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%ecutsm = zero
 call etsf_io_low_read_var(ncid, "ecutsm", hdr%ecutsm, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%etot = zero
 call etsf_io_low_read_var(ncid, "etot", hdr%etot, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%intxc = 0
 call etsf_io_low_read_var(ncid, "intxc", hdr%intxc, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%ixc = 1
 call etsf_io_low_read_var(ncid, "ixc", hdr%ixc, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%occopt = 1
 call etsf_io_low_read_var(ncid, "occopt", hdr%occopt, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%pertcase = 0
 call etsf_io_low_read_var(ncid, "pertcase", hdr%pertcase, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%qptn(:) = 0
 call etsf_io_low_read_var(ncid, "qptn", hdr%qptn, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%residm = zero
 call etsf_io_low_read_var(ncid, "residm", hdr%residm, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%stmbias = zero
 call etsf_io_low_read_var(ncid, "stmbias", hdr%stmbias, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%tphysel  = zero
 call etsf_io_low_read_var(ncid, "tphysel", hdr%tphysel, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%tsmear = zero
 call etsf_io_low_read_var(ncid, "tsmear", hdr%tsmear, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%ecutdg = hdr%ecut
 call etsf_io_low_read_var(ncid, "ecutdg", hdr%ecutdg, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

!test for old wavefunction style
 if(hdr%ecutsm>tol6 .and. headform<44 .and. .not.(fform==51.or.fform==52.or.fform==101.or.fform==102)) then
   write(msg,'(a,es16.6,13a)' )&
&   'The value of ecutsm is',hdr%ecutsm,', while the file has been produced prior to v4.4 .',ch10,&
&   'The definition of the smearing function has changed,',' so that you are not allowed',ch10,&
&   'to restart from a old wavefunction file. By contrast,',' you can restart from an old',ch10,&
&   'potential or density file, and perform a self-consistent',' cycle with a new ABINIT version.',ch10,&
&   'Action: produce a density or potential file using the old',' version of ABINIT, and restart from it.'
   MSG_ERROR(msg)
 end if

!Multidimensional variables.
!The case of istwfk is always 1, since ETSF don't use the time reversal symetry.
 hdr%istwfk(:) = 1
 call etsf_io_low_read_var(ncid, "istwfk", hdr%istwfk, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%pspcod(:) = 0
 call etsf_io_low_read_var(ncid, "pspcod", hdr%pspcod, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%pspdat(:) = 0
 call etsf_io_low_read_var(ncid, "pspdat", hdr%pspdat, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%pspso(:) = 0
 call etsf_io_low_read_var(ncid, "pspso", hdr%pspso, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%pspxc(:) = 0
 call etsf_io_low_read_var(ncid, "pspxc", hdr%pspxc, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%so_psp(:) = 1
 call etsf_io_low_read_var(ncid, "so_psp", hdr%so_psp, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%symafm(:) = 1
 call etsf_io_low_read_var(ncid, "symafm", hdr%symafm, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%title(:) = ""
 call etsf_io_low_read_var(ncid, "title", hdr%title, 132, lstat, error_data = error_data)
 if (lstat) then ! Pad the returned string with " " instead of "\0"
   do itypat = 1, size(hdr%title), 1
     call strip(hdr%title(itypat))
   end do
 end if

 call etsf_io_low_read_var(ncid, "zionpsp", hdr%zionpsp, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 call etsf_io_low_read_var(ncid, "znuclpsp", hdr%znuclpsp, lstat, error_data = error_data)
 ETSF_WARN(lstat,error_data)

 hdr%lmn_size = 1
 if (headform>=44) then ! Compared to 4.2, add lmn_size and
   call etsf_io_low_read_var(ncid, "lmn_size", hdr%lmn_size, lstat, error_data = error_data)
!  DC: remove the PAW specific variables, they are not in accordance
!  with 56_io_mpi/hdr_io.F90. This part should be rewritten totally.
   if (hdr%usepaw==1) then
     write(msg, '(9a)' )&
&     'The support for the internal variables of PAW are not yet',ch10,&
&     'available with ETSF output. Restarting calculation from this',ch10,&
&     'will not be possible.',ch10,&
&     'Action: produce a density or potential file using the old',ch10,&
&     'binary format of ABINIT, and restart from it.'
     MSG_ERROR(msg)
     do iatom=1,hdr%natom
       hdr%pawrhoij(iatom)%ngrhoij = 0
       hdr%pawrhoij(iatom)%lmnmix_sz = 0
       hdr%pawrhoij(iatom)%use_rhoij_ = 0
       hdr%pawrhoij(iatom)%use_rhoijres = 0
     end do
!    !!    call etsf_io_low_read_var(ncid, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
!    !!    allocate (rhoij(rhoijdim1,hdr%nspden,hdr%natom))
!    !!    call etsf_io_low_read_var(ncid, "rhoij", rhoij, lstat, error_data = error_data)
!    !!    if (.not.lstat) goto 1000
!    !!
!    !!    cplex=1;if (rhoijdim1/=hdr%natom) cplex=2
!    !!    call pawrhoij_alloc(hdr%pawrhoij,cplex,hdr%nspden,hdr%nspinor,hdr%nsppol,hdr%typat,lmnsize=hdr%lmn_size)
!    !!    do iatom=1,hdr%natom
!    !!     itypat=hdr%typat(iatom)
!    !!     lmn2_size=hdr%lmn_size(itypat)*(hdr%lmn_size(itypat)+1)/2
!    !!     nselect=0
!    !!     if (cplex==1) then
!    !!      do ilmn=1,lmn2_size
!    !!       if (any(abs(rhoij(ilmn,:,iatom))>tol10)) then
!    !!        nselect=nselect+1
!    !!        hdr%pawrhoij(iatom)%rhoijselect(nselect)=ilmn
!    !!        do ispden=1,hdr%nspden
!    !!         hdr%pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoij(ilmn,ispden,iatom)
!    !!        end do
!    !!       end if
!    !!      end do
!    !!     else
!    !!      do ilmn=1,lmn2_size
!    !!       if (any(abs(rhoij(2*ilmn-1:2*ilmn,:,iatom))>tol10)) then
!    !!        nselect=nselect+1
!    !!        hdr%pawrhoij(iatom)%rhoijselect(nselect)=ilmn
!    !!        hdr%pawrhoij(iatom)%rhoijp(2*nselect-1,ispden)=rhoij(2*ilmn-1,ispden,iatom)
!    !!        hdr%pawrhoij(iatom)%rhoijp(2*nselect  ,ispden)=rhoij(2*ilmn  ,ispden,iatom)
!    !!       end if
!    !!      end do
!    !!     end if
!    !!     if (nselect<lmn2_size) then
!    !!      hdr%pawrhoij(iatom)%rhoijselect(nselect+1:lmn2_size)=0
!    !!      do ispden=1,hdr%nspden
!    !!       hdr%pawrhoij(iatom)%rhoijp(cplex*nselect+1:cplex*lmn2_size,ispden)=zero
!    !!      end do
!    !!     end if
!    !!     hdr%pawrhoij(iatom)%nrhoijsel=nselect
!    !!    end do
!    !!    deallocate(rhoij)
   end if
 end if

 !  BigDFT private variables.
 !  First implementation, we assume that the number of wavelet resolutions
 !  is 2. Latter, we may add this value to hdr.
 lstat = .true.

#else
 MSG_ERROR("ETSF-IO support not activated")
#endif

end subroutine hdr_ncread
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_read
!! NAME
!! hdr_read
!!
!! FUNCTION
!! This subroutine deals with the reading of the hdr_type structured variables 
!! It support both Fortran format and Netcdf.
!! Note that, when reading, different records of hdr are allocated here, 
!! Records of hdr should be deallocated correctly by a call to hdr_free when hdr is not used anymore.
!!
!! INPUTS
!!  unit=Fortran unit number or the unit of the open NetCDF file.
!!  io_mode= IO_MODE_FORTRAN for Fortran IO, IO_MODE_ETSF for netcdf file.
!!  [rewind]=True to rewind the file (valid only if IO_MODE_FORTRAN). Default: False
!!
!! OUTPUT
!!  fform=kind of the array in the file. if the reading fails, return fform=0
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_read(Hdr,fform,unit,io_mode,rewind)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_read'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit,io_mode
 integer,intent(out) :: fform
 type(hdr_type),target,intent(out) :: Hdr
 logical,optional,intent(in) :: rewind

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

 select case (io_mode)
 case (IO_MODE_FORTRAN, IO_MODE_MPI)
   if (present(rewind)) then
     call hdr_fort_read(Hdr,unit,fform,rewind=rewind)
   else
     call hdr_fort_read(Hdr,unit,fform)
   end if

 case (IO_MODE_ETSF)
   call hdr_ncread(Hdr,unit,fform)

 case default
   write(msg,'(a,i0)')"Wrong value for io_mode: ",io_mode
   MSG_ERROR(msg)
 end select

end subroutine hdr_read
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_fort_write
!! NAME
!! hdr_fort_write
!!
!! FUNCTION
!!  Writes the header and fform to unformatted file
!!
!! INPUTS
!!  Hdr<hdr_type>=The header of the file.
!!  fform=kind of the array in the file
!!  unit=unit number of the unformatted file 
!!  [rewind]=True to rewind the file. Default: False
!!
!! OUTPUT
!!  ierr=Exit status
!!
!! NOTES
!! The file is supposed to be open already
!!
!! PARENTS
!!      m_header
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_fort_write(Hdr,unit,fform,ierr,rewind)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_fort_write'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(out) :: ierr
 integer,intent(in) :: unit,fform
 logical,optional,intent(in) :: rewind
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
 integer :: headform,ipsp

!*************************************************************************

 ! Change intent to in. Change pawrhoij_io first!
 ierr = 0

! natom,nkpt,npsp,ntypat... are not defined in this section: always address them from hdr
 if (present(rewind)) then
   if (rewind) rewind(unit)
 end if

!Writing always use last format version
 headform=57
!headform= HDR_LATEST_HEADFORM TODO  rationalize hdr methods and deps. so that we can use the vars. in m_header.
 write(unit,err=10) hdr%codvsn, headform, fform

 write(unit,err=10) hdr%bantot, hdr%date, hdr%intxc, hdr%ixc, &
& hdr%natom, hdr%ngfft(1:3), hdr%nkpt, &
& hdr%nspden, hdr%nspinor, &
& hdr%nsppol, hdr%nsym, hdr%npsp, hdr%ntypat, hdr%occopt, hdr%pertcase,&
& hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, &
& hdr%qptn, hdr%rprimd, hdr%stmbias, hdr%tphysel, hdr%tsmear, &
& hdr%usewvl

 write(unit,err=10) hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:),&
& hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
& hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
& hdr%tnons(:,:), hdr%znucltypat(:), hdr%wtk(:)

 do ipsp=1,hdr%npsp
   write(unit,err=10) hdr%title(ipsp), hdr%znuclpsp(ipsp), &
&   hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
&   hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp)
 end do

 write(unit,err=10) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie

 if (hdr%usepaw==1) then
   call pawrhoij_io(hdr%pawrhoij,unit,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,headform,"Write")
 end if

 return

 ! Handle IO-error.
10 ierr=1

end subroutine hdr_fort_write
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_ncwrite
!! NAME
!! hdr_ncwrite
!!
!! FUNCTION
!! This subroutine deals with the output of the hdr_type structured variables in ETSF+NETCDF fornat.
!! It handles variables according to the ETSF format, whenever possible and uses new variables 
!!  when not available in the ETSF format. 
!!
!! INPUTS
!!  fform=kind of the array in the file
!!  ncid=the unit of the open NetCDF file.
!!  [nc_define]=Optional flag. If True, the basic dimensions required by the ETSF specification
!!    are written. Default: False.
!!
!! OUTPUT
!!  Only writing 
!!
!! PARENTS
!!      m_header
!!
!! CHILDREN
!!      abi_etsf_dims_init_low,etsf_io_basisdata_def,etsf_io_basisdata_put
!!      etsf_io_dims_def,etsf_io_electrons_def,etsf_io_electrons_put
!!      etsf_io_geometry_def,etsf_io_geometry_put,etsf_io_kpoints_def
!!      etsf_io_kpoints_put,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_var,ini_wf_etsf
!!      ncid_define_basedims
!!
!! SOURCE

subroutine hdr_ncwrite(hdr,ncid,fform,nc_define)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,fform
 logical,optional,intent(in) :: nc_define
 type(hdr_type),target,intent(in) :: hdr

#ifdef HAVE_TRIO_ETSF_IO
!Local variables-------------------------------
!scalars
 integer,parameter :: headform=44
 integer :: mband
 integer :: flags_electrons,flags_kpoints,flags_basisdata,flags_geometry
 real(dp),target :: ecut, fermie
 logical :: lstat,k_dep,my_define
 character(len=etsf_charlen),target :: basis_set
 character(len=500) :: msg
 type(etsf_dims) :: dims
 type(etsf_kpoints) :: kpoints
 type(etsf_basisdata) :: basisdata
 type(etsf_geometry) :: geometry
 type(etsf_electrons) :: electrons
 type(etsf_io_low_error) :: error_data
!arrays
 real(dp), target :: rprimd(3, 3)
! integer :: cplex, ilmn, irhoij, ispden, lmn2_size, nselect,rhoijdim1,iatom,itypat,
! real(dp), allocatable :: rhoij(:,:,:)

! *************************************************************************

 !  Writing the header of an unformatted file
 my_define = .FALSE.; if (PRESENT(nc_define)) my_define = nc_define

 if (my_define) then
   ! TODO: use ab_etsf_init but rewrite the interface so that it 
   ! can be easily called without passing dtset
   k_dep = .TRUE.
   call etsf_io_low_set_define_mode(ncid, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call ncid_define_basedims(ncid) 

   mband = MAXVAL(Hdr%nband)
   call abi_etsf_dims_init_low(Dims,Hdr%natom,Hdr%ntypat,mband,Hdr%nkpt,Hdr%nsppol,Hdr%nspinor,Hdr%nspden,Hdr%nsym)

   call etsf_io_dims_def(ncid, Dims,lstat,error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   flags_geometry = etsf_geometry_all - etsf_geometry_space_group - &
&     etsf_geometry_chemical_symbols - etsf_geometry_pseudo_types 

   call etsf_io_geometry_def(ncid, lstat, error_data, k_dependent=k_dep, flags=flags_geometry)
   ETSF_CHECK_ERROR(lstat,error_data)

   ! TODO some info on the k-points is missing in the header
   flags_kpoints = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights
                                                                                                
   call etsf_io_kpoints_def(ncid, lstat, error_data, k_dependent=k_dep, flags=flags_kpoints)
   ETSF_CHECK_ERROR(lstat,error_data)

   flags_electrons = etsf_electrons_all - etsf_electrons_x_functional - etsf_electrons_c_functional - &
&    etsf_electrons_eigenvalues 

   call etsf_io_electrons_def(ncid, lstat, error_data, k_dependent=k_dep, flags=flags_electrons)
   ETSF_CHECK_ERROR(lstat,error_data)

   call ini_wf_etsf(ncid,Hdr%usewvl,Hdr%lmn_size,Hdr%npsp,Hdr%ntypat)
   !ETSF_CHECK_ERROR(lstat,error_data)

   flags_basisdata = etsf_basisdata_basis_set
   if (Hdr%usewvl==0) then
     flags_basisdata = flags_basisdata + etsf_basisdata_kin_cutoff + etsf_basisdata_n_coeff
   end if

   call etsf_io_basisdata_def(ncid, lstat, error_data, k_dependent=k_dep, flags=flags_basisdata)
   ETSF_CHECK_ERROR(lstat,error_data)
 end if

 ! Switch to write mode.
 call etsf_io_low_set_write_mode(ncid, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 ! Associate and write values to ETSF groups.
 if (hdr%usewvl == 0) then  
   ! Plane wave case.
   ecut = hdr%ecut
   basis_set =  "plane_waves"
   basisdata%basis_set              => basis_set
   basisdata%kinetic_energy_cutoff  => ecut
   basisdata%number_of_coefficients => hdr%npwarr
 else  
   ! Wavelet case.
   basis_set = "daubechies_wavelets"
   basisdata%basis_set => basis_set
   ! Required variable than should enter the standard.
   call etsf_io_low_write_var(ncid, "number_of_wavelets", hdr%nwvlarr, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)
 end if

 call etsf_io_basisdata_put(ncid, basisdata, lstat, error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 fermie = hdr%fermie
 electrons%fermi_energy            => fermie
 electrons%number_of_states%data1D => hdr%nband
 electrons%occupations%data1D      => hdr%occ

 call etsf_io_electrons_put(ncid, electrons, lstat, error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 rprimd = hdr%rprimd
 geometry%primitive_vectors             => rprimd
 geometry%reduced_symmetry_matrices     => hdr%symrel
 geometry%atom_species                  => hdr%typat
 geometry%reduced_symmetry_translations => hdr%tnons
 geometry%reduced_atom_positions        => hdr%xred
 geometry%atomic_numbers                => hdr%znucltypat
 if (hdr%npsp == hdr%ntypat) then
   geometry%valence_charges => hdr%zionpsp
 end if

 call etsf_io_geometry_put(ncid, geometry, lstat, error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 kpoints%reduced_coordinates_of_kpoints => hdr%kptns
 kpoints%kpoint_weights                 => hdr%wtk

 call etsf_io_kpoints_put(ncid, kpoints, lstat, error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

!Write non-ETSF variables.
!TODO istwfk
 call etsf_io_low_write_var(ncid, "date", hdr%date, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "codvsn", hdr%codvsn, 6, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "ecut_eff", hdr%ecut_eff, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "ecutsm", hdr%ecutsm, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "etot", hdr%etot, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 !MG Why 44
 call etsf_io_low_write_var(ncid, "headform", headform, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "fform", fform, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "intxc", hdr%intxc, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "ixc", hdr%ixc, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "occopt", hdr%occopt, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "pertcase", hdr%pertcase, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "residm", hdr%residm, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "stmbias", hdr%stmbias, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "tphysel", hdr%tphysel, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "tsmear", hdr%tsmear, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

!Version 44 add usepaw ecutdg
 call etsf_io_low_write_var(ncid, "ecutdg", hdr%ecutdg, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "usepaw", hdr%usepaw, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

!Array variables.
 call etsf_io_low_write_var(ncid, "istwfk", hdr%istwfk, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "pspcod", hdr%pspcod, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "pspdat", hdr%pspdat, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "pspso", hdr%pspso, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "pspxc", hdr%pspxc, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "qptn", hdr%qptn, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "so_psp", hdr%so_psp, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "symafm", hdr%symafm, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "title", hdr%title, 132, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 call etsf_io_low_write_var(ncid, "znuclpsp", hdr%znuclpsp, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

 !if (hdr%npsp /= hdr%ntypat) then
 call etsf_io_low_write_var(ncid, "zionpsp", hdr%zionpsp, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)
 !end if

!Version 44 add lmn_size and rhoij
 call etsf_io_low_write_var(ncid, "lmn_size", hdr%lmn_size, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

!DC: remove the PAW specific variables, they are not in accordance
!with 56_io_mpi/hdr_io.F90. This part should be rewritten totally.
 if (hdr%usepaw == 1) then
   write(msg,'(9a)')&
&   'The support for the internal variables of PAW are not yet',ch10,&
&   'available with ETSF output. Restarting calculation from this',ch10,&
&   'will not be possible.',ch10,&
&   'Action: produce a density or potential file using the old',ch10,&
&   'binary format of ABINIT, and restart from it.'
   MSG_WARNING(msg)
!  !!   rhoijdim1 = maxval(hdr%lmn_size)
!  !!   rhoijdim1 = hdr%pawrhoij(1)%cplex*rhoijdim1*(rhoijdim1+1)/2
!  !!   call etsf_io_low_write_var(ncid, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
!  !!   allocate (rhoij(rhoijdim1,hdr%nspden,hdr%natom))
!  !!   do iatom=1,hdr%natom
!  !!    itypat=hdr%typat(iatom)
!  !!    lmn2_size = hdr%lmn_size(itypat)*(hdr%lmn_size(itypat)+1)/2
!  !!    cplex=hdr%pawrhoij(iatom)%cplex
!  !!    do ispden=1,hdr%nspden
!  !!     rhoij(1:cplex*lmn2_size,ispden,iatom)=zero
!  !!     if (cplex==1) then
!  !!      do irhoij=1,hdr%pawrhoij(iatom)%nrhoijsel
!  !!       ilmn=hdr%pawrhoij(iatom)%rhoijselect(irhoij)
!  !!       rhoij(ilmn,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(irhoij,ispden)
!  !!      end do
!  !!     else
!  !!      do irhoij=1,hdr%pawrhoij(iatom)%nrhoijsel
!  !!       ilmn=hdr%pawrhoij(iatom)%rhoijselect(irhoij)
!  !!       rhoij(2*ilmn-1,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(2*irhoij-1,ispden)
!  !!       rhoij(2*ilmn  ,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(2*irhoij  ,ispden)
!  !!      end do
!  !!     end if
!  !!    end do
!  !!   end do
!  !!   call etsf_io_low_write_var(ncid, "rhoij", rhoij, lstat, error_data = error_data)
!  !!   if (.not.lstat) goto 1000
!  !!   deallocate (rhoij)
 end if

 ! BigDFT variables.
 call etsf_io_low_write_var(ncid, "usewvl", hdr%usewvl, lstat, error_data = error_data)
 ETSF_CHECK_ERROR(lstat,error_data)

#else
 MSG_ERROR("ETSF-IO support not activated")
#endif

end subroutine hdr_ncwrite
!!***

END MODULE m_header
!!***
