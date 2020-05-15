subroutine make_points(rprimd,nfft,n1,n2,n3,points)

  implicit none 

  integer :: nfft,n1,n2,n3
  real(8) :: rprimd(3,3),aa,bb,cc
  real(8) :: points(3,nfft)

  integer :: q, i,j,k, counter 

  ! a b c are stored in column wise way 
  ! rprimd = | a1 b1 c1|
  !          | a2 b2 c2|
  !          | a3 b3 c3|

  counter  = 0
  do k=1,n3
      do j=1,n2
          do i=1,n1
            counter  = counter + 1

            aa = (i-1)/dble(n1)
            bb = (j-1)/dble(n2)
            cc = (k-1)/dble(n3)
            
            do q=1,3
               points(q,counter) = aa*rprimd(q,1) + bb*rprimd(q,2) + cc*rprimd(q,3)
            enddo

          enddo
      enddo
  enddo

end subroutine 
