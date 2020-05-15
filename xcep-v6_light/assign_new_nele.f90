!========================================
! notify new number of electrons to 
! subsystems 
! the subroutine will make sure the 
! subsystems have received the orders
!
!========================================
subroutine assign_new_nele(nsys,sub_ele)

  implicit none

  integer,intent(in) :: nsys
  real(kind=8),intent(in) :: sub_ele(nsys)
 
  ! local vars ......
  !==================
  integer :: s
  character(len=500) :: message

  ! >>>>>>>> FUNCTIONS <<<<<<<<< !

  do s=1,nsys
   ! write(message,'(a,f16.10)')'new_nele ',sub_ele(s)
   ! call write_msg(s,message)
   ! write(6,'(a,f12.8,a,i2)') & 
   !   'assign_new_nele(): new_ele =',sub_ele(s),' <= subsys:',s
  enddo

  ! wait for the nodes reading the new electron number
  do s=1,nsys
    call wait_for_system(s)
  enddo

end subroutine assign_new_nele
