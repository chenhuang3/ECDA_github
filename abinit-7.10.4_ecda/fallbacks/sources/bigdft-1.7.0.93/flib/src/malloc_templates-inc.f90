!> @file
!! Include fortran file for allocation templates
!! file included in module dynamic_memory.f90
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

subroutine i1_all(array,m)
  use metadata_interfaces, metadata_address => geti1
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer, dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i1_all

subroutine i1_all_free(array)
  use metadata_interfaces, metadata_address => geti1
  implicit none
  integer, dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine i1_all_free

subroutine i2_all(array,m)
  use metadata_interfaces, metadata_address => geti2
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer, dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i2_all

subroutine i2_all_free(array)
  use metadata_interfaces, metadata_address => geti2
  implicit none
  integer, dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine i2_all_free

subroutine i3_all(array,m)
  use metadata_interfaces, metadata_address => geti3
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i3_all

subroutine i3_all_free(array)
  use metadata_interfaces, metadata_address => geti3
  implicit none
  integer, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine i3_all_free

subroutine i4_all(array,m)
  use metadata_interfaces, metadata_address => geti4
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer, dimension(:,:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i4_all

subroutine i4_all_free(array)
  use metadata_interfaces, metadata_address => geti4
  implicit none
  integer, dimension(:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine i4_all_free

!character arrays
subroutine c1_all(array,m)
  use metadata_interfaces, metadata_address => getc1
  implicit none
  type(malloc_information_str_all), intent(in) :: m
  character(len=m%len), dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-c-inc.f90'
end subroutine c1_all

!subroutine c1_all_free(length,array)
subroutine f_free_str(length,array)
  use metadata_interfaces, metadata_address => getc1
  implicit none
  integer, intent(in) :: length
  character(len=length), dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-c-inc.f90' 
end subroutine f_free_str


subroutine l1_all(array,m)
  use metadata_interfaces, metadata_address => getl1
  implicit none
  type(malloc_information_all), intent(in) :: m
  logical, dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine l1_all

subroutine l1_all_free(array)
  use metadata_interfaces, metadata_address => getl1
  implicit none
  logical, dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine l1_all_free


subroutine l2_all(array,m)
  use metadata_interfaces, metadata_address => getl2
  implicit none
  type(malloc_information_all), intent(in) :: m
  logical, dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine l2_all

subroutine l2_all_free(array)
  use metadata_interfaces, metadata_address => getl2
  implicit none
  logical, dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine l2_all_free

subroutine l3_all(array,m)
  use metadata_interfaces, metadata_address => getl3
  implicit none
  type(malloc_information_all), intent(in) :: m
  logical, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine l3_all

subroutine l3_all_free(array)
  use metadata_interfaces, metadata_address => getl3
  implicit none
  logical, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine l3_all_free

subroutine r1_all(array,m)
  use metadata_interfaces, metadata_address => getr1
  implicit none
  type(malloc_information_all), intent(in) :: m
  real, dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine r1_all

subroutine r1_all_free(array)
  use metadata_interfaces, metadata_address => getr1
  implicit none
  real, dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine r1_all_free

subroutine r2_all(array,m)
  use metadata_interfaces, metadata_address => getr2
  implicit none
  type(malloc_information_all), intent(in) :: m
  real, dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine r2_all

subroutine r2_all_free(array)
  use metadata_interfaces, metadata_address => getr2
  implicit none
  real, dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine r2_all_free

subroutine r3_all(array,m)
  use metadata_interfaces, metadata_address => getr3
  implicit none
  type(malloc_information_all), intent(in) :: m
  real, dimension(:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90' 
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine r3_all

subroutine r3_all_free(array)
  use metadata_interfaces, metadata_address => getr3
  implicit none
  real, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine r3_all_free

subroutine d1_all(array,m)
  use metadata_interfaces, metadata_address => getdp1
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d1_all

subroutine d1_all_free(array)
  use metadata_interfaces, metadata_address => getdp1
  implicit none
  double precision, dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine d1_all_free


subroutine d2_all(array,m)
  use metadata_interfaces, metadata_address => getdp2
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d2_all

subroutine d2_all_free(array)
  use metadata_interfaces, metadata_address => getdp2
  implicit none
  double precision, dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine d2_all_free

subroutine d3_all(array,m)
  use metadata_interfaces, metadata_address => getdp3
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90' 
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d3_all

subroutine d3_all_free(array)
  use metadata_interfaces, metadata_address => getdp3
  implicit none
  double precision, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine d3_all_free

subroutine d4_all(array,m)
  use metadata_interfaces, metadata_address => getdp4
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90' 
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3),&
       m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d4_all

subroutine d4_all_free(array)
  use metadata_interfaces, metadata_address => getdp4
  implicit none
  double precision, dimension(:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine d4_all_free

subroutine d5_all(array,m)
  use metadata_interfaces, metadata_address => getdp5
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90' 
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3),&
       m%lbounds(4):m%ubounds(4),m%lbounds(5):m%ubounds(5)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d5_all

subroutine d5_all_free(array)
  use metadata_interfaces, metadata_address => getdp5
  implicit none
  double precision, dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine d5_all_free

subroutine d6_all(array,m)
  use metadata_interfaces, metadata_address => getdp6
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90' 
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3),&
       m%lbounds(4):m%ubounds(4),m%lbounds(5):m%ubounds(5),&
       m%lbounds(6):m%ubounds(6)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d6_all

subroutine d6_all_free(array)
  use metadata_interfaces, metadata_address => getdp6
  implicit none
  double precision, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine d6_all_free



!test to see if this is convenient
subroutine d1_all_free_multi(arrayA,arrayB,arrayC,arrayD,arrayE,arrayF,arrayG,arrayH)
  implicit none
  double precision, dimension(:), allocatable, intent(inout) :: arrayA
  double precision, dimension(:), allocatable, intent(inout) :: arrayB
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayC
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayD
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayE
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayF
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayG
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayH

  call d1_all_free(arrayA)
  call d1_all_free(arrayB)
  if (present(arrayC)) then
     call d1_all_free(arrayC)
  end if
  if (present(arrayD)) then
     call d1_all_free(arrayD)
  end if
  if (present(arrayE)) then
     call d1_all_free(arrayE)
  end if
  if (present(arrayF)) then
     call d1_all_free(arrayF)
  end if
  if (present(arrayG)) then
     call d1_all_free(arrayG)
  end if
  if (present(arrayH)) then
     call d1_all_free(arrayH)
  end if
end subroutine d1_all_free_multi

!pointers
subroutine d1_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp1ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:), pointer, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)

  include 'allocate-inc.f90'
end subroutine d1_ptr

subroutine d1_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp1ptr
  implicit none
  double precision, dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-inc.f90'
  nullify(array)
end subroutine d1_ptr_free

subroutine i1_ptr(array,m)
  use metadata_interfaces, metadata_address => geti1ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  integer, dimension(:), pointer, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)

  include 'allocate-inc.f90'
end subroutine i1_ptr

subroutine i1_ptr_free(array)
  use metadata_interfaces, metadata_address => geti1ptr
  implicit none
  integer, dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-inc.f90'
  nullify(array)
end subroutine i1_ptr_free


subroutine d2_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp2ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d2_ptr

subroutine d2_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp2ptr
  implicit none
  double precision, dimension(:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-inc.f90'
  nullify(array)
end subroutine d2_ptr_free

subroutine i2_ptr(array,m)
  use metadata_interfaces, metadata_address => geti2ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  integer, dimension(:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i2_ptr

subroutine i2_ptr_free(array)
  use metadata_interfaces, metadata_address => geti2ptr
  implicit none
  integer, dimension(:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-inc.f90'
  nullify(array)
end subroutine i2_ptr_free


subroutine d3_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp3ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:,:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d3_ptr

subroutine d3_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp3ptr
  implicit none
  double precision, dimension(:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-inc.f90'
  nullify(array)
end subroutine d3_ptr_free

subroutine d4_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp4ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:,:,:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d4_ptr

subroutine d4_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp4ptr
  implicit none
  double precision, dimension(:,:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-inc.f90'
  nullify(array)
end subroutine d4_ptr_free

subroutine d5_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp5ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:,:,:,:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4),&
       m%lbounds(5):m%ubounds(5)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d5_ptr

subroutine d5_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp5ptr
  implicit none
  double precision, dimension(:,:,:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-inc.f90'
  nullify(array)
end subroutine d5_ptr_free

!character arrays
subroutine c1_ptr(array,m)
  use metadata_interfaces, metadata_address => getc1ptr
  implicit none
  type(malloc_information_str_ptr), intent(in) :: m
  character(len=m%len), dimension(:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-c-inc.f90'
end subroutine c1_ptr

!subroutine c1_ptr_free(length,array)
subroutine f_free_str_ptr(length,array)
  use metadata_interfaces, metadata_address => getc1ptr
  implicit none
  integer, intent(in) :: length
  character(len=length), dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-c-inc.f90'
  nullify(array)
end subroutine f_free_str_ptr
