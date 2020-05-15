!
! read zion.dat file from global system 
!
! created by Chen Huang 1/10/2020
subroutine read_system_zion(natom,zion)

  implicit none 
  integer :: natom, iatom
  real(8) :: zion(natom)


  open(file='global_system/zion.dat',action='read',unit=111,form='unformatted')
  do iatom = 1,natom   
     read(111)zion(iatom)
  enddo 
  close(111)

end subroutine 
