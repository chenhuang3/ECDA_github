!=========================================
! make the q-space vectors
!=========================================
subroutine make_q_vector(n1,n2,n3,gprim,qvec)
  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=8),intent(out) :: qvec(3,(n1/2+1)*n2*n3)
  real(kind=8),intent(in)  :: gprim(3,3)   ! this should be stored as [g1|g2|g3]
  
  ! Local vars
  integer      :: i,j,k,index,dim1
  real(kind=8) :: ivect(3)
  
  dim1 = n1/2+1   ! FFTW stores the fourier component 
                  ! onto half of the 1st dimension 
  do k=1,n3
     ivect(3)=dble(k)-1.0d0
     if (k-1 > n3/2) ivect(3)=ivect(3)-dble(n3)
     do j=1,n2
        ivect(2)=dble(j)-1.0d0
        if (j-1 > n2/2) ivect(2)=ivect(2)-dble(n2)
        do i=1,dim1
          ivect(1)=dble(i)-1.0d0
          ! follow OFDFT's routine
          ! gprim is [g1|g2|g3]
          index = i + (j-1)*dim1 + (k-1)*dim1*n2
          qvec(:,index) = MATMUL(gprim,ivect)
        enddo
     enddo
  enddo

  return
end subroutine make_q_vector
