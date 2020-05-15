

subroutine dipole_moment(logf,natom,nspin,n1,n2,n3,dvol,points,rho,xcart,zion,dipole)
  
  implicit none 

  integer :: n1,n2,n3,nspin,logf,natom 
 
  real(8) :: points(3,n1*n2*n3), & 
             rho(n1*n2*n3,nspin), & 
             dipole(3), dvol, & 
             xcart(3,natom), & 
             zion(natom), val_q


  integer :: counter, i, j, k
  
  ! compute dipole
  
  counter = 0
  dipole = 0.0d0

  do k=1,n3
     do j=1,n2
        do i=1,n1
          counter = counter + 1
          dipole(1) = dipole(1) + points(1,counter) * sum(rho(counter,:)) * dvol
          dipole(2) = dipole(2) + points(2,counter) * sum(rho(counter,:)) * dvol
          dipole(3) = dipole(3) + points(3,counter) * sum(rho(counter,:)) * dvol
        enddo
     enddo
  enddo

  do i=1,natom 
     val_q = zion(i) 
     dipole(1) = dipole(1) - val_q*xcart(1,i)
     dipole(2) = dipole(2) - val_q*xcart(2,i)
     dipole(3) = dipole(3) - val_q*xcart(3,i)
  enddo

  dipole = - dipole 

  ! unti conversion 
  dipole = dipole * 1.602E-19 * 0.5291772016E-10 ! to SI first 
  dipole = dipole / 3.33564E-30 ! to Debye  http://cccbdb.nist.gov/debye.asp

  write(logf,'(a,3f9.4,a,f9.4,a)')'dipole moment (x,y,z): ',dipole(1:3), & 
     ' (Debye)   total dipole: ',sqrt(sum(dipole**2)),' (Debye)'
end subroutine dipole_moment




