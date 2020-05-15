!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_kptrank
!! NAME
!! m_kptrank
!!
!! FUNCTION
!! This module deals with rank objects for hashing k-point vector lists
!!
!! COPYRIGHT
!! Copyright (C) 2010-2014 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_kptrank

 use defs_basis
 use m_profiling_abi
 use m_errors

 implicit none

 private
!!***

!!****t* m_kptrank/kptrank_type
!! NAME
!! kptrank_type
!! 
!! FUNCTION
!!  structure to contain a rank/inverse rank pair of arrays, with dimensions
!! 
!! SOURCE

 type,public :: kptrank_type
   integer :: max_linear_density
   integer :: max_rank
   integer :: npoints
   integer,allocatable :: invrank(:)
   integer,allocatable :: rank(:)
   integer,allocatable :: multipl(:)
 end type kptrank_type

 public :: mkkptrank       ! Sets up the kpt ranks for comparing kpts
 public :: get_rank_1kpt   ! Calculates the rank for one kpt
 public :: copy_kptrank    ! Copy the object
 public :: destroy_kptrank ! Free memory
 public :: dump_kptrank    ! Prints the arrays and dimensions of a kptrank_type structure
!!***

contains
!!***

!!****f* m_kptrank/mkkptrank
!!
!! NAME
!! mkkptrank
!!
!! FUNCTION
!! This routine sets up the kpt ranks for comparing kpts
!!
!! INPUTS
!!  npt = number of kpoints
!!  kpt = coordinates of kpoints
!!
!! OUTPUT
!!  krank = object containing ranking and inverse ranking
!!
!! PARENTS
!!      get_full_kgrid,m_double_grid,m_gamma,m_nesting,m_tetrahedron,mkfskgrid
!!      mkqptequiv,new_integrate_gamma,order_fs_kpts,outelph,printbxsf
!!      read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkkptrank (kpt,nkpt,krank,nsym,symrec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkkptrank'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 integer,intent(in), optional :: nsym
!arrays
 type(kptrank_type), intent(out) :: krank
 real(dp),intent(in) :: kpt(3,nkpt)
 integer,intent(in), optional :: symrec(3,3, *)

!Local variables -------------------------
!scalars
 integer :: ikpt, isym, symkptrank, irank
 real(dp) :: smallestlen
 character(len=500) :: msg
!arrays
 real(dp) :: symkpt(3)

! *********************************************************************

! find smallest linear length
 smallestlen = one
 do ikpt=1, nkpt
   if (abs(kpt(1,ikpt)) > tol10) &
&     smallestlen = min(smallestlen, abs(kpt(1,ikpt)))
   if (abs(kpt(2,ikpt)) > tol10) &
&     smallestlen = min(smallestlen, abs(kpt(2,ikpt)))
   if (abs(kpt(3,ikpt)) > tol10) &
&     smallestlen = min(smallestlen, abs(kpt(3,ikpt)))
 end do

 krank%max_linear_density = int(one/smallestlen)+1
 krank%max_rank = 2*krank%max_linear_density**3
 krank%npoints = nkpt

 ABI_ALLOCATE(krank%rank,(nkpt))
 ABI_CHECK_ALLOC("out of memory %rank")

 ABI_ALLOCATE(krank%invrank,(krank%max_rank))
 ABI_CHECK_ALLOC("out of memory %invrank")
 krank%invrank(:) = -1

!Ensure kpt(i)+one is positive, and the smallest
!difference between kpts should be larger than 1/100
!ie ngkpt < 100.
 do ikpt=1,nkpt
   call get_rank_1kpt (kpt(:,ikpt), krank%rank(ikpt), krank)

   if (krank%rank(ikpt) > krank%max_rank .or. krank%rank(ikpt) < 1) then
     write(msg,'(a,2i0)')" max rank exceeded or < 1, ikpt, rank ", ikpt, krank%rank(ikpt)
     MSG_ERROR(msg)
   end if
   krank%invrank(krank%rank(ikpt)) = ikpt
 end do
 
! if symrec is provided, fill invrank with appropriate irred kpt indices
! for symmetry completion: kptrank_t%invrank points to the irred k-point
! equivalent to the k-point whose rank is provided
 if (present(symrec)) then
   ABI_CHECK(present(nsym), "need both symrec and nsym arguments together")
   do ikpt=1,nkpt
     do isym = 1, nsym
       symkpt = matmul(symrec(:,:,isym), kpt(:, ikpt))
       
       call get_rank_1kpt (symkpt(:), symkptrank, krank)

       krank%invrank(symkptrank) = ikpt
     end do
   end do
 end if 

 ABI_ALLOCATE(krank%multipl,(nkpt))
 ABI_CHECK_ALLOC("out of memory %multipl")

 krank%multipl = 0
! find multiplicity of ikpt
 do irank = 1, krank%max_rank
   ikpt = krank%invrank(irank)
   if (ikpt > 0) krank%multipl(ikpt) = krank%multipl(ikpt) + 1
 end do

end subroutine mkkptrank
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/get_rank_1kpt
!!
!! NAME
!! get_rank_1kpt
!!
!! FUNCTION
!! This routine calculates the rank for one kpt
!!
!! INPUTS
!!  kpt = coordinates of kpoints
!!  krank = rank object for the k-grid we are using
!!
!! OUTPUT
!!  rank = rank of the kpoint
!!
!! PARENTS
!!      elphon,get_full_kgrid,integrate_gamma,integrate_gamma_alt,k_neighbors
!!      m_gamma,m_kptrank,m_nesting,m_tetrahedron,mkfskgrid,mkqptequiv
!!      printbxsf,read_el_veloc,read_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_rank_1kpt(kpt,rank,krank)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_rank_1kpt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: rank
 type(kptrank_type), intent(in) :: krank
!arrays
 real(dp),intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!arrays    
 real(dp) :: redkpt(3)

! *************************************************************************

! wrap to [0, 1[ -> replaced call to wrap2_zeroone inline, to encapsulate this module
 if (kpt(1)>zero) then
   redkpt(1)=mod((kpt(1)+tol12),one)-tol12
 else
   redkpt(1)=-mod(-(kpt(1)-one+tol12),one)+one-tol12
 end if
 if(abs(redkpt(1))<tol12)redkpt(1)=0.0_dp

 if (kpt(2)>zero) then
   redkpt(2)=mod((kpt(2)+tol12),one)-tol12
 else
   redkpt(2)=-mod(-(kpt(2)-one+tol12),one)+one-tol12
 end if
 if(abs(redkpt(2))<tol12)redkpt(2)=0.0_dp

 if (kpt(3)>zero) then
   redkpt(3)=mod((kpt(3)+tol12),one)-tol12
 else
   redkpt(3)=-mod(-(kpt(3)-one+tol12),one)+one-tol12
 end if
 if(abs(redkpt(3))<tol12)redkpt(3)=0.0_dp



! rank = int(real(krank%max_linear_density)*(redkpt(3)+half+tol8 +&
!&           real(krank%max_linear_density)*(redkpt(2)+half+tol8 +&
!&           real(krank%max_linear_density)*(redkpt(1)+half+tol8))))
 rank = int(real(krank%max_linear_density)*(redkpt(1)+half+tol8 +&
&           real(krank%max_linear_density)*(redkpt(2)+half+tol8 +&
&           real(krank%max_linear_density)*(redkpt(3)+half+tol8))))

 if (rank > krank%max_rank) then
   write(msg,'(a,i0)') ' rank should be inferior to ', krank%max_rank
   MSG_ERROR(msg)
 end if

end subroutine get_rank_1kpt
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/copy_kptrank
!!
!! NAME
!! copy_kptrank
!!
!! FUNCTION
!! Copy the object 
!!
!! INPUTS
!!
!! OUTPUT
!!  krank = object containing ranking and inverse ranking, to be deallocated
!!
!! PARENTS
!!      defs_elphon,elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_kptrank (krank_in, krank_out)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_kptrank'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kptrank_type), intent(in) :: krank_in
 type(kptrank_type), intent(out) :: krank_out

! *********************************************************************
 krank_out%max_linear_density = krank_in%max_linear_density
 krank_out%max_rank = krank_in%max_rank
 krank_out%npoints = krank_in%npoints
 
 ABI_ALLOCATE(krank_out%rank,(krank_out%npoints))
 krank_out%rank = krank_in%rank
 
 ABI_ALLOCATE(krank_out%invrank,(krank_out%max_rank))
 krank_out%invrank = krank_in%invrank

 ABI_ALLOCATE(krank_out%multipl,(krank_out%npoints))
 krank_out%multipl = krank_in%multipl
 
end subroutine copy_kptrank
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/destroy_kptrank
!!
!! NAME
!! destroy_kptrank
!!
!! FUNCTION
!! This routine deallocates the arrays in a kptrank_type structure
!!
!! INPUTS
!!  krank = object containing ranking and inverse ranking, to be deallocated
!!
!! PARENTS
!!      defs_elphon,get_full_kgrid,m_double_grid,m_gamma,m_nesting
!!      m_tetrahedron,mkfskgrid,mkqptequiv,new_integrate_gamma
!!      new_integrate_gamma_tr,new_integrate_gamma_tr_lova,order_fs_kpts
!!      outelph,printbxsf,read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_kptrank (krank)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_kptrank'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kptrank_type), intent(inout) :: krank

! *********************************************************************

 if (allocated(krank%rank))  then
   ABI_DEALLOCATE(krank%rank)
 end if
 if (allocated(krank%invrank))  then
   ABI_DEALLOCATE(krank%invrank)
 end if
 if (allocated(krank%multipl))  then
   ABI_DEALLOCATE(krank%multipl)
 end if

end subroutine destroy_kptrank
!!***

!----------------------------------------------------------------------

!!****f* m_kptrank/dump_kptrank
!!
!! NAME
!! dump_kptrank
!!
!! FUNCTION
!! This routine prints the arrays and dimensions of a kptrank_type structure
!!
!! INPUTS
!!  krank = object containing ranking and inverse ranking
!!  unout = unit for open file to print to
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dump_kptrank (krank, unout)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dump_kptrank'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: unout
!arrays
 type(kptrank_type), intent(in) :: krank

! *********************************************************************

  write(unout, *) 
  write(unout, '(a)') ' Dump of the contents of a kptrank_type structure with k-point rank information' 
  write(unout, '(a,I8)') ' max linear density of points in 3 directions: max_linear_density = ',  krank%max_linear_density
  write(unout, '(a,I8)') ' maximum rank for any point in grid: max_rank = ',  krank%max_rank
  write(unout, '(a,I8)') ' number of points in input grid: npoints = ',  krank%npoints
  write(unout, *) 
  write(unout, '(a)') ' invrank array = '
  write(unout, '(I4)') krank%invrank(:)
  write(unout, '(a)') ' rank array = '
  write(unout, '(I4)') krank%rank(:)
  write(unout, '(a)') ' multiplicity array = '
  write(unout, '(I4)') krank%multipl(:)
  write(unout, *) 

end subroutine dump_kptrank
!!***

!----------------------------------------------------------------------

end module m_kptrank
!!***
