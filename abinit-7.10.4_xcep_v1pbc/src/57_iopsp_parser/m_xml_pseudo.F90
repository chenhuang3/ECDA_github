!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xml_pseudo
!! NAME
!! m_xml_pseudo
!!
!! FUNCTION
!! This module reads a pseudopotential file written in XML.
!! A full example of the building up of a data structure using
!! the SAX paradigm.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2013 ABINIT group (JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! 
!! OUTPUT
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


module m_xml_pseudo

 use m_profiling_abi
 use m_errors
#if defined HAVE_TRIO_FOX
 use fox_sax
 use m_xml_pseudo_types
#endif

implicit none

private

#if defined HAVE_TRIO_FOX
!
! It defines the routines that are called from xml_parser in response
! to particular events.
!
public  :: begin_element, end_element, pcdata_chunk

logical, private  :: in_vps = .false. , in_radfunc = .false.
logical, private  :: in_semilocal = .false. , in_header = .false.
logical, private  :: in_coreCharge = .false. , in_data = .false.
logical, private  :: in_valenceCharge = .false.
logical, private  :: in_pseudowavefun = .false. , in_pswf = .false.

integer, private, save  :: ndata

type(pseudo_t), public, target, save :: pseudo
type(grid_t), private, save        :: grid
type(grid_t), private, save        :: global_grid
!
! Pointers to make it easier to manage the data
!
type(header_t), private, pointer   :: hp
type(vps_t), private, pointer      :: pp
type(pswf_t), private, pointer     :: pw
type(radfunc_t), private, pointer  :: rp

CONTAINS  !===========================================================
!!***

!!****f* m_xml_pseudo/begin_element
!! NAME
!! begin_element
!!
!! FUNCTION
!!  Read an XML tag with a given name.
!!  Fills the present module private data.
!!
!! INPUTS
!!  namespaceURI = universal resource indicator for XML namespace??? Not used.
!!  localName = local equivalent of tag name?? Not used.
!!  name = name of XML tag which has been read in
!!  attributes = attributes of XML tag
!!
!! OUTPUT
!!  Fills private data in present module.
!!
!! PARENTS
!!
!! CHILDREN
!!      build_data_array
!!
!! SOURCE
subroutine begin_element(namespaceURI,localName,name,attributes)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'begin_element'
!End of the abilint section

character(len=*),intent(in)   :: namespaceURI,localName,name
type(dictionary_t),intent(in) :: attributes

character(len=100)  :: value

!Just to fool abirules
 value=localName
 value=namespaceURI

select case(name)

      case ("pseudo")
         pseudo%npots  = 0
         pseudo%npswfs = 0
         global_grid%npts = 0
!         value = getValue(attributes,"version")
!         if (value == "0.5") then
!            write(std_out,*) "Processing a PSEUDO version 0.5 XML file"
!            pseudo%npots  = 0
!            pseudo%npswfs = 0
!            global_grid%npts = 0
!         else
!            write(std_out,*) "Can only work with PSEUDO version 0.5 XML files"
!            STOP
!         end if

      case ("header")
         in_header = .true.
         hp => pseudo%header

         value = getValue(attributes,"symbol")
         if (value == "" ) then
           MSG_ERROR("Cannot determine atomic symbol")
         end if
         hp%symbol = value

         value = getValue(attributes,"atomic-number")
         if (value == "" ) then
           MSG_ERROR("Cannot determine atomic-number")
         end if
         read(unit=value,fmt=*) hp%atomicnumber

         value = getValue(attributes,"zval")
         if (value == "" ) then
           MSG_ERROR("Cannot determine zval")
         end if
         read(unit=value,fmt=*) hp%zval

         value = getValue(attributes,"xc-functional-type")
         if (value == "" ) then
           MSG_ERROR("Cannot determine xc-functional-type")
         end if
         hp%xcfunctionaltype = value

         value = getValue(attributes,"xc-functional-parametrization")
         if (value == "" ) then
           MSG_ERROR("Cannot determine xc-functional-parametrization ")
         end if
         hp%xcfunctionalparametrization = value

         value = getValue(attributes,"creator")
         if (value == "" ) value = "unknown"
         hp%creator = value

         value = getValue(attributes,"date")
         if (value == "" ) value = "unknown"
         hp%date = value

         value = getValue(attributes,"flavor")
         if (value == "" ) value = "unknown"
         hp%flavor = value

         value = getValue(attributes,"relativistic")
         if (value == "" ) value = "no"
         hp%relativistic = (value == "yes")

         value = getValue(attributes,"polarized")
         if (value == "" ) value = "no"
         hp%polarized = (value == "yes")

         value = getValue(attributes,"core-corrections")
         if (value == "" ) value = "nc"
         hp%core_corrections = value

      case ("vps")
         in_vps = .true.

         pseudo%npots = pseudo%npots + 1
         pp => pseudo%pot(pseudo%npots)
         rp => pp%V                       ! Pointer to radial function

         value = getValue(attributes,"l")
         if (value == "" ) then
           MSG_ERROR("Cannot determine l for Vps")
         end if
         pp%l = value

         value = getValue(attributes,"principal-n")
         if (value == "" ) then
           MSG_ERROR("Cannot determine n for Vps")
         end if
         read(unit=value,fmt=*) pp%n

         value = getValue(attributes,"cutoff")
         if (value == "" ) then
           MSG_ERROR("Cannot determine cutoff for Vps")
         end if
         read(unit=value,fmt=*) pp%cutoff

         value = getValue(attributes,"occupation")
         if (value == "" ) then
           MSG_ERROR("Cannot determine occupation for Vps")
         end if
         read(unit=value,fmt=*) pp%occupation

         value = getValue(attributes,"spin")
         if (value == "" ) then
           MSG_ERROR("Cannot determine spin for Vps")
         end if
         read(unit=value,fmt=*) pp%spin

      case ("grid")

         value = getValue(attributes,"type")
         if (value == "" ) then
           MSG_ERROR("Cannot determine grid type")
         end if
         grid%type = value

         value = getValue(attributes,"npts")
         if (value == "" ) then
           MSG_ERROR("Cannot determine grid npts")
         end if
         read(unit=value,fmt=*) grid%npts

         value = getValue(attributes,"scale")
         if (value == "" ) then
           MSG_ERROR("Cannot determine grid scale")
         end if
         read(unit=value,fmt=*) grid%scale

         value = getValue(attributes,"step")
         if (value == "" ) then
           MSG_ERROR("Cannot determine grid step")
         end if
         read(unit=value,fmt=*) grid%step

         !
         ! In this way we allow for a private grid for each radfunc,
         ! or for a global grid specification
         !
         if (in_radfunc) then
            rp%grid = grid
         else
            global_grid = grid
         end if

      case ("data")
         in_data = .true.
         if (rp%grid%npts == 0) then 
           MSG_ERROR("Grid not specified correctly")
         end if
         ABI_ALLOCATE(rp%data,(rp%grid%npts))
         ndata = 0             ! To start the build up

      case ("radfunc")
         in_radfunc = .true.
         rp%grid = global_grid     ! Might be empty
                                   ! There should then be a local grid element
                                   ! read later

      case ("pseudocore-charge")
         in_coreCharge = .true.
         rp => pseudo%core_charge

      case ("valence-charge")
         in_valenceCharge = .true.
         rp => pseudo%valence_charge

      case ("semilocal")
         in_semilocal = .true.

         value = getValue(attributes,"npots-down")
         if (value == "" ) then
           MSG_ERROR("Cannot determine npots-down")
         end if
         read(unit=value,fmt=*) pseudo%npots_down

         value = getValue(attributes,"npots-up")
         if (value == "" ) then
           MSG_ERROR("Cannot determine npots-up")
         end if
         read(unit=value,fmt=*) pseudo%npots_up

      case ("pseudowave-functions")
         in_pseudowavefun = .true.

      case ("pswf")
         in_pswf = .true.

         pseudo%npswfs = pseudo%npswfs + 1

         pw => pseudo%pswf(pseudo%npswfs)
         rp => pw%V                       ! Pointer to radial function

         value = getValue(attributes,"l")
         if (value == "" ) then
           MSG_ERROR("Cannot determine l for Vps")
         end if
         pw%l = value

         value = getValue(attributes,"principal-n")
         if (value == "" ) then
           MSG_ERROR("Cannot determine n for Vps")
         end if
         read(unit=value,fmt=*) pw%n

         value = getValue(attributes,"spin")
         if (value == "" ) then
           MSG_ERROR("Cannot determine spin for Vps")
         end if
         read(unit=value,fmt=*) pw%spin

end select

end subroutine begin_element
!!***

!!****f* m_xml_pseudo/end_element
!! NAME
!! end_element
!!
!! FUNCTION
!!  End XML tag effect: switches flags in private data of this module
!!
!! INPUTS
!!  namespaceURI = universal resource indicator for XML namespace??? Not used.
!!  localName = local equivalent of tag name?? Not used.
!!  name = name of XML tag which has been read in
!!
!! OUTPUT
!!  side effect: private data flags in present module are turned to .false.
!!
!! PARENTS
!!
!! CHILDREN
!!      build_data_array
!!
!! SOURCE
subroutine end_element(namespaceURI,localName,name)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'end_element'
!End of the abilint section

character(len=*),intent(in) :: namespaceURI,localName,name

character(len=100)  :: value

!Just to fool abirules 
 value=localName
 value=namespaceURI

select case(name)

      case ("vps")
         in_vps = .false.

      case ("radfunc")
         in_radfunc = .false.

      case ("data")
      !
      ! We are done filling up the radfunc data
      ! Check that we got the advertised number of items
      !
         in_data = .false.
         if (ndata /= size(rp%data)) then 
           MSG_ERROR("npts mismatch")
         end if

      case ("pseudocore-charge")
         in_coreCharge = .false.

      case ("valence-charge")
         in_valenceCharge = .false.

      case ("semilocal")
         in_semilocal = .false.

      case ("pseudowave-functions")
         in_pseudowavefun = .false.

      case ("pswf")
         in_pswf = .false.

      case ("pseudo")
!         call dump_pseudo(pseudo)

end select

end subroutine end_element
!!***

!!****f* m_xml_pseudo/pcdata_chunk
!! NAME
!! pcdata_chunk
!!
!! FUNCTION
!!   Read in data from XML structure, if we are in a valid data field
!!   for the present XML structures
!!
!! INPUTS
!!   chunk = raw data for chunk of XML data
!!
!! OUTPUT
!!   copied and translated into module data (side effect)
!!
!! PARENTS
!!
!! CHILDREN
!!      build_data_array
!!
!! SOURCE
subroutine pcdata_chunk(chunk)

use m_xml_converters

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pcdata_chunk'
!End of the abilint section

character(len=*), intent(in) :: chunk


if (len_trim(chunk) == 0) RETURN     ! skip empty chunk

if (in_data) then
!
! Note that we know where we need to put it through the pointer rp...
!

     call build_data_array(chunk,rp%data,ndata)

else if (in_header) then
      !
      ! There should not be any pcdata in header in this version...
      !write(std_out,*) "Header data:"
      !write(std_out,*) trim(chunk)
end if

end subroutine pcdata_chunk
!!***

#endif

end module m_xml_pseudo
!!***
