!
! compute the dot production of two vectors in fourier space  
! we have the option to add a metric in the fourier space 
!
function dotprod_four(n1,n2,n3,vec1,vec2,ucvol,qvec,use_damp_metric) 

 implicit none  

 integer :: n1,n2,n3
 real(8) :: vec1(n1,n2,n3)
 real(8) :: vec2(n1,n2,n3),  &  
            dotprod_four, & 
            ucvol, & 
            qvec(3,(n1/2+1)*n2*n3) 
 logical :: use_damp_metric

 ! 
 real(8)    :: metric(n1/2+1,n2,n3), & 
               qnorm(n1/2+1,n2,n3), &
               v1(n1,n2,n3), & 
               v2(n1,n2,n3), & 
               w,q0, & 
               len_scale=3.0d0, & ! ! more weight for wave length over 3 bohr
               Amax = 5.d0, &    ! max weight for q=0
               aa, bb 

 complex(8) :: fft1((n1/2+1),n2,n3)
 complex(8) :: fft2((n1/2+1),n2,n3)

 ! function 
 call fft(n1,n2,n3,reshape(vec1,(/n1,n2,n3/)),fft1,1)
 call fft(n1,n2,n3,reshape(vec2,(/n1,n2,n3/)),fft2,1)

 if (use_damp_metric) then 
   qnorm = reshape(sum(qvec**2,1),(/n1/2+1,n2,n3/))
   qnorm(1,1,1) = 1e8
   q0 = minval(qnorm)
   qnorm(1,1,1) = 0.d0
   w = (2.d0*3.1415926d0/len_scale)**2 ! more weight for wave length over 3 bohr
   !
   ! The model for metric is 
   ! 
   ! metric = 1 + aa * w / (qnorm + bb)
   ! 
   ! metric -> 0           as qnorm -> infinit
   ! metric -> Amax/2      as qnorm -> w
   !
   bb = (Amax/2.d0*w-w)*2.d0/3.d0/Amax
   aa = (Amax-1.d0)*bb/w
   metric = 1.d0 + aa*w/(qnorm+bb) 
   fft1 = fft1 / sqrt(metric)
   fft2 = fft2 / sqrt(metric)
   call fft(n1,n2,n3,v1,fft1,-1)
   call fft(n1,n2,n3,v2,fft2,-1)
 else 
   metric = 1.d0 
 endif 

 dotprod_four = sum(v1*v2)

 return                                                                          

 !! DEBUG 
 !test_fft(1,1,1) = 1;
 !test_fft(2,1,1) = 2;
 !test_fft(3,1,1) = 8;
 !test_fft(1,2,1) = 3;
 !test_fft(2,2,1) = 4;
 !test_fft(3,2,1) = 4;
 !test_fft(1,1,2) = 5;
 !test_fft(2,1,2) = 6;
 !test_fft(3,1,2) = 6; 
 !test_fft(1,2,2) = 7;
 !test_fft(2,2,2) = 8;
 !test_fft(3,2,2) = 8; 
 !print *, 'dot in real space:', sum(test_fft*test_fft)
 !if (myrank==0) then 
 !  call fft(3,2,2,reshape(test_fft,(/3,2,2/)),fft2,1)
 !  print *,fft2
 !  print *,sum(fft2(1,1:2,1:2)**2) & 
 !        + sum(fft2(2:3/2+1,1:2,1:2)*conjg(fft2(2:3/2+1,1:2,1:2)))*2.d0
 !
 !  print *, 'dot in real space:', sum(global_vks(:,1)**2)
 !  call fft(n1,n2,n3,reshape(global_vks(:,1),(/n1,n2,n3/)),fft1,1)
 !  print *,sum(fft1(1,1:n2,1:n3)**2) & 
 !        + sum(fft1(2:n1/2+1,1:n2,1:n3)*conjg(fft1(2:n1/2+1,1:n2,1:n3)))*2.d0
 ! 
 !
 !  stop
 !endif 


! print *,fft2(1,1,1)**2+sum(fft1(2:n1/2+1,2:n2,2:n3)*conjg(fft1(2:n1/2+1,2:n2,2:n3)))*2.d0

endfunction 
