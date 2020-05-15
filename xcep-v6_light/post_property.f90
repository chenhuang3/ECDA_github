
! 
! calculate all properties here if needed
! 
subroutine post_property(nsys,nspin,n1,n2,n3,dvol,rhor)

  implicit none

  ! exteranl vars 
  integer :: nsys, nspin, n1,n2,n3
  real(8) :: rhor(n1,n2,n3,nspin,nsys), & 
             dvol

  ! local vars 
  integer :: isub

  if (nspin==2) then 
    print *,' '
    write(6,'(a)')      & 
    '   -------- post_property.f90 ---------'
    do isub=1,nsys
    write(6,'(a,i4,a,f12.4)') & 
    '   magnetism of subsystem ',isub," :",sum(rhor(:,:,:,1,isub))*dvol-sum(rhor(:,:,:,2,isub))*dvol
    enddo
    write(6,'(a,f12.4)')& 
    '   total magnetism: ',sum(rhor(:,:,:,1,:))*dvol-sum(rhor(:,:,:,2,:))*dvol
    print *,''
  endif


end subroutine post_property
