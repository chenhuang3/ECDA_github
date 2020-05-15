!{\src2tex{textfont=tt}}
!!****f* ABINIT/Lij
!! NAME
!! Lij
!!
!! FUNCTION
!! Routine which computes PAW onsite angular momentum expectation values
!!
!! COPYRIGHT
!! Copyright (C) 2005-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ntypat :: number of types of atoms
!!  type(pawrad_type) :: pawrad(ntypat) Data on PAW radial grid
!!  type(pawtab_type) :: pawtab(ntypat) Tabulated PAW data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  type(bfield_type) :: dtbfield. The onsite Lij values and Lij/r^3
!!                                 Lij/r^3 are stored in dtbfield
!!
!! NOTES
!!
!! PARENTS
!!      initorbmag
!!
!! CHILDREN
!!      pawrad_deducer0,simp_gen,slxyzs
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine Lij(dtbfield,ntypat,pawrad,pawtab)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_bfield

 use m_pawrad, only : pawrad_type, pawrad_deducer0, simp_gen
 use m_pawtab, only : pawtab_type
 use m_sphharm, only : slxyzs

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Lij'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: ntypat
 type(bfield_type),intent(inout) :: dtbfield

!arrays
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: idir, itypat, ilmn,il,im,iln,ilm, jlmn,jl,jm,jlm,jln,j0lmn
 integer :: klmn,kln, mesh_size
 real(dp) :: intg,intgr3
 complex(dpc) :: lms
!arrays
 integer, ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp),allocatable :: ff(:)

! *************************************************************************

!loop over types of atoms in cell
 do itypat = 1, ntypat
   indlmn => pawtab(itypat)%indlmn
   mesh_size=pawrad(itypat)%mesh_size
   ABI_ALLOCATE(ff,(mesh_size))

!  loop over basis state pairs for this type
   do jlmn=1,pawtab(itypat)%lmn_size
     jl=indlmn(1,jlmn)
     jm=indlmn(2,jlmn)
     jlm=indlmn(4,jlmn)
     jln=indlmn(5,jlmn)
     j0lmn=jlmn*(jlmn-1)/2
     do ilmn=1,jlmn
       il=indlmn(1,ilmn)
       im=indlmn(2,ilmn)
       iln=indlmn(5,ilmn)
       ilm=indlmn(4,ilmn)
       klmn=j0lmn+ilmn
       kln = pawtab(itypat)%indklmn(2,klmn)

!      Computation of <phi_i|phi_j>- <tphi_i|tphi_j> radial integral
!      this is NOT the same as the sij non-local overlap, because that also
!      involves an angular integral over the S_i*S_j spherical harmonics
       ff(2:mesh_size)=pawtab(itypat)%phiphj(2:mesh_size,kln)-&
&       pawtab(itypat)%tphitphj(2:mesh_size,kln)
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))

       ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln)-&
&       pawtab(itypat)%tphitphj(2:mesh_size,kln))/&
&       pawrad(itypat)%rad(2:mesh_size)**3
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intgr3,ff,pawrad(itypat))

       do idir = 1, 3

         call slxyzs(il,im,idir,jl,jm,lms)
         dtbfield%Lij(1,klmn,itypat,idir)=intg*real(lms)
         dtbfield%Lij(2,klmn,itypat,idir)=intg*aimag(lms)
         dtbfield%Lijr3(1,klmn,itypat,idir)=intgr3*real(lms)
         dtbfield%Lijr3(2,klmn,itypat,idir)=intgr3*aimag(lms)

       end do

     end do ! end loop over ilmn
   end do ! end loop over jlmn


   ABI_DEALLOCATE(ff)

 end do ! end loop over atom types

 dtbfield%has_Lij = 2
 dtbfield%has_Lijr3 = 2

end subroutine Lij
!!***
