!{\src2tex{textfont=tt}}
!!****f* ABINIT/mklocl_wavelets
!!
!! NAME
!! mklocl_wavelets
!!
!! FUNCTION
!! Compute the ionic local potential when the pseudo-potentials are GTH, using
!! the special decomposition of these pseudo. The resulting potential is computed with
!! free boundary conditions. It gives the same result than mklocl_realspace for the
!! GTH pseudo only with a different way to compute the potential.
!!
!! Optionally compute :
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  xcart(3,natom)=cartesian atomic coordinates.
!!
!! OUTPUT
!!  (if option==1) vpsp(nfft)=the potential resulting from the ionic
!!                 density of charge.
!!  (if option==2) grtn(3,natom)=grads of Etot wrt tn. These gradients are in
!!                 reduced coordinates. Multiply them by rprimd to get
!!                 gradients in cartesian coordinates.
!!
!! PARENTS
!!      mklocl,wvl_wfsinp_scratch
!!
!! CHILDREN
!!      createionicpotential,local_forces,psolver_hartree,wrtout
!!      wvl_rhov_abi2big,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mklocl_wavelets(efield, icoulomb, grtn, mpi_enreg, natom, nfft, &
     & nscforder, nspden, option, rhor, rprimd, vpsp, wvl_den, wvl_descr, xcart)

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_abi2big, only : wvl_rhov_abi2big
 use m_profiling_abi
 use m_xmpi
 use m_errors
#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only : createIonicPotential, local_forces
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mklocl_wavelets'
 use interfaces_14_hidewrite
 use interfaces_62_poisson
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option, natom, nfft, nspden, nscforder, icoulomb
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(wvl_internal_type), intent(in) :: wvl_descr
!arrays
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),efield(3)
 real(dp),intent(inout) :: grtn(3,natom)
 real(dp), intent(inout) :: vpsp(nfft)
 real(dp),intent(inout) :: xcart(3,natom)

!Local variables-------------------------------
!scalars
 integer :: i,i1,i2,i3,ia,ierr,igeo,me,nproc,shift,comm
 real(dp) :: charge
 real(dp) :: energy
 real(dp) :: locstrten(6,4)
 character(len=500) :: message
!arrays
 real(dp) :: epot(3),elecfield(3)
 real(dp),allocatable :: gxyz(:,:),vhartr(:)
! type(wvl_coulomb_operator)::kernel
 character(len=20) :: atomnames(100)

! *********************************************************************

 elecfield=zero !not used here

#if defined HAVE_DFT_BIGDFT
 comm=mpi_enreg%comm_wvl
 nproc=xcomm_size(comm)
 me=xcomm_rank(comm)


 shift = 1 + wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2) * &
& wvl_den%denspot%dpbox%nscatterarr(me, 4)
 if (option == 1) then
!!  We get the kernel for the Poisson solver (used to go from the ionic
!!  charge to the potential).
!!  If the kernel is uncomputed, it does it now.
!   call psolver_kernel(wvl_den%denspot%dpbox%hgrids, 2, icoulomb, me, kernel, &
!&   comm, wvl_den%denspot%dpbox%ndims, nproc, nscforder)
!   if (.not. associated(kernel%co%kernel)) then
!     call psolver_kernel(wvl_den%denspot%dpbox%hgrids, 1, icoulomb, me, kernel, &
!&     comm, wvl_den%denspot%dpbox%ndims, nproc, nscforder)
!   end if

   write(message, '(a,a)' ) ch10,&
&   ' mklocl_wavelets: Create local potential from ions.'
   call wrtout(std_out,message,'COLL')

!  Call the BigDFT routine...
   call createIonicPotential(wvl_descr%atoms%astruct%geocode, me, nproc, (me == 0), wvl_descr%atoms, &
&   xcart, wvl_den%denspot%dpbox%hgrids(1), wvl_den%denspot%dpbox%hgrids(2), &
&   wvl_den%denspot%dpbox%hgrids(3), &
&   elecfield, wvl_descr%Glr%d%n1, wvl_descr%Glr%d%n2, wvl_descr%Glr%d%n3, &
&   wvl_den%denspot%dpbox%n3pi, wvl_den%denspot%dpbox%nscatterarr(me, 3) + 1, &
&   wvl_den%denspot%dpbox%ndims(1), wvl_den%denspot%dpbox%ndims(2), &
&   wvl_den%denspot%dpbox%ndims(3), wvl_den%denspot%pkernel, vpsp(shift), 0.d0,wvl_descr%rholoc)

!  copy vpsp into bigdft object:
   call wvl_rhov_abi2big(nfft,1,1,vpsp,wvl_den%denspot%v_ext)
   wvl_den%denspot%v_ext = reshape(vpsp, shape(wvl_den%denspot%v_ext))

   if (maxval(efield) > zero) then
     write(message, '(a,a)' ) ch10,&
&     ' mklocl_wavelets: Add the electric field.'
     call wrtout(std_out,message,'COLL')

!    We add here the electric field since in BigDFT, the field must be on x...
     epot(:) = real(0.5, dp) * efield(:) * wvl_den%denspot%dpbox%hgrids(:)
     do i3 = 1, wvl_den%denspot%dpbox%n3pi, 1
       ia = (i3 - 1) * wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2)
       do i2 = -14, 2 * wvl_descr%Glr%d%n2 + 16, 1
         i = ia + (i2 + 14) * wvl_den%denspot%dpbox%ndims(1)
         do i1 = -14, 2 * wvl_descr%Glr%d%n1 + 16, 1
           i = i + 1
           vpsp(shift + i) = vpsp(shift + i) + &
&           epot(1) * real(i1 - wvl_descr%Glr%d%n1, dp) + &
&           epot(2) * real(i2 - wvl_descr%Glr%d%n2, dp) + &
&           epot(3) * real(i3 - wvl_descr%Glr%d%n3, dp)
         end do
       end do
     end do
   end if

 else if (option == 2) then
!  Dummy arguments
   atomnames(:) = "'Unknown name'"

!  Compute forces
   write(message, '(a)' ) 'mklocl_wavelets: compute local forces.'
   call wrtout(std_out,message,'COLL')

!  Compute Hartree's potential from rhor.
   ABI_ALLOCATE(vhartr,(nfft))
   call psolver_hartree(energy, wvl_den%denspot%dpbox%hgrids, icoulomb, me, comm, nfft, &
&   wvl_den%denspot%dpbox%ndims, nproc, nscforder, nspden, rhor, vhartr, 1)

!  Allocate temporary array for forces.
   ABI_ALLOCATE(gxyz,(3, natom))

!  calculate local part of the forces grtn (BigDFT routine)
   call local_forces(me, wvl_descr%atoms, xcart, &
&   wvl_den%denspot%dpbox%hgrids(1), wvl_den%denspot%dpbox%hgrids(2), &
&   wvl_den%denspot%dpbox%hgrids(3), &
&   wvl_descr%Glr%d%n1, wvl_descr%Glr%d%n2, wvl_descr%Glr%d%n3, &
&   wvl_den%denspot%dpbox%nscatterarr(me, 2), wvl_den%denspot%dpbox%nscatterarr(me, 3) + 1, &
&   wvl_den%denspot%dpbox%ndims(1), wvl_den%denspot%dpbox%ndims(2), rhor(shift, 1), &
&   vhartr(shift),gxyz,locstrten,charge)
   ABI_DEALLOCATE(vhartr)

!  Pending: floc,locstrten and charge are not used here.
!  Pending: put mpi_enreg%nscatterarr... in object denspot, initialize object, etc.

   if (nproc > 1) then
     call xmpi_sum(gxyz, comm, ierr)
   end if

!  Forces should be in reduced coordinates.
   do ia = 1, natom, 1
     do igeo = 1, 3, 1
       grtn(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) - &
&       rprimd(2, igeo) * gxyz(2, ia) - &
&       rprimd(3, igeo) * gxyz(3, ia)
     end do
   end do

!  Deallocate local variables
   ABI_DEALLOCATE(gxyz)
 else ! option switch
   message = ' mklocl_wavelets : internal error, option should be 1 or 2.'
   MSG_ERROR(message)
 end if
 
#else
 BIGDFT_NOTENABLED_ERROR()
#endif

end subroutine mklocl_wavelets
!!***
