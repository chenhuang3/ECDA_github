
!=======================================
! just to debug the fft() subroutine 
!
!  Chen Huang  Jan/2011
!=======================================
subroutine debug_fft_code
 implicit none 
 !! DEBUG FFT
 integer,parameter :: dp=8
 real(kind=8) :: arr(24),c_arr(2,18) ! dimension is 4x3x2
                                 ! the complex one is 3x3x2
 arr = (/2.0,3.0,2.0,6.0, & 
         7.0,3.0,4.0,9.0, & 
         5.0,3.0,2.0,7.0, &
         6.0,2.0,4.0,3.0, &
         9.0,3.0,1.0,9.0, &
         4.0,2.0,4.0,7.0/)

 call fft(24,4,3,2,arr,c_arr,1)
 print *, 'arr =========' 
 print *, arr
 print *, 'REAL PART c_arr =========='
 print *, c_arr(1,:)
 print *, 'IMAG PART c_arr =========='
 print *, c_arr(2,:)

 print *,'back transform, has set arr=0'
 arr = 0.0
 call fft(24,4,3,2,arr,c_arr,-1)
 print *, 'arr =========' 
 print *, arr
 print *, 'REAL PART c_arr =========='
 print *, c_arr(1,:)
 print *, 'IMAG PART c_arr =========='
 print *, c_arr(2,:)
 stop
 !!ENDDEBUG OF FFT
end subroutine  debug_fft_code



!====================================
! DEBUG grad_lap 
!
!====================================
subroutine debug_grad_lap(nfft,a,b,c,gprim,n1,n2,n3,dvol)
! DEBUG gradient and laplacian code
! 
  implicit none
  integer :: i,ii,nfft,n1,n2,n3
  real(kind=8)  :: a(3),b(3),c(3),gprim(3,3)
  real(kind=8)  :: u(nfft),tmp3d(3,nfft),dvol
  real(kind=8)  :: qvec(3,(n1/2+1)*n2*n3),lap(nfft),vw_energy,gnorm

!!  open(file='den1.dat',unit=111,action='read')
!!  open(file='./h_dimer/totden.dat',unit=111,action='read')
  open(file='./sum1d',unit=111,action='read')
  print *,' readin sum1d '
  do ii=1,nfft
    read(111,*) u(ii)
  enddo
  close(111)
  print *,'max/min den1.dat=', minval(u),maxval(u)
  print *,'sum (den) = ', sum(u) * dvol
  call make_gprim(a,b,c,gprim)
  call make_q_vector(nfft,n1,n2,n3,gprim,qvec)
  print *,'max/min qvec=', minval(qvec),maxval(qvec)
!  call gradient(nfft,n1,n2,n3,u,qvec,tmp3d)
!  call laplacian(nfft,n1,n2,n3,u,qvec,lap)
print *,'computing gradient ...'
  call gradient(nfft,n1,n2,n3,sqrt(u),qvec,tmp3d)
  print *,' grd =',minval(tmp3d),maxval(tmp3d)
print *,'computing laplacian ...'
  call laplacian(nfft,n1,n2,n3,sqrt(u),qvec,lap)
  print *,' lap =',minval(lap),maxval(lap)

 vw_energy = 0.d0
 do i=1,nfft
   gnorm = sqrt(dot_product(tmp3d(:,i),tmp3d(:,i)))
!   vw_energy = vw_energy + 0.5d0 * gnorm**2 
   vw_energy = vw_energy - 0.5d0 * sqrt(u(i))*lap(i) 
 enddo
 vw_energy = vw_energy * dvol

  print *,' vw = ', vw_energy ,' stop in debug_grad_lap()'
  stop  
! ENDDEBUG  
end subroutine debug_grad_lap

