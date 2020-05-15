! regularize the Ax=b equation 
! three options: 1) regularize |x|
!                2) regularize grad(x)
!                3) regularize |x| and grad(x)
!
! see code for details 
!
! On exit, Ax is updated
!
subroutine regularization(reg_type,reg_lambda,nspin,nfft,cell_nfft,qvec,Ax,x)

 implicit none 
 integer :: reg_type, nfft, nspin, & 
            cell_nfft(3), & 
            isp 
 real(8) :: Ax(nfft,nspin), & 
            x(nfft,nspin),  & 
            reg_lambda, & 
            qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), & 
            lap(nfft)

 select case(reg_type) 
 case(1)
   !
   ! Regularizing x (Tikhonov regularization)
   ! regularize the norm <x,x>
   ! f(x) = 1/2<x,Ax> - <b,x> + 1/2*reg_lambda*<x,x>
   ! df/dx --> (A+reg_lambda*I)x = b
   ! which means that the nearly zero eigenvalues of A are added with reg_lambda
   !
   Ax = Ax + x*reg_lambda
 case(2)
   !
   ! regularize gradient 
   ! the objective function is 1/2*<x,Ax> - <b,x> + 1/2<grad x,grad x>*lambda
   ! the functional derivative of 1/2*<grad_x, grad_x> is - lambda*\laplacian x
   ! in the end we solve
   !       (A - \laplacian*lambda)|x> = |b>
   !
   do isp=1,nspin
      call laplacian(nfft,cell_nfft(1),cell_nfft(2),cell_nfft(3),x(:,isp),qvec,lap)
      Ax(:,isp) = Ax(:,isp) - reg_lambda*lap
   enddo
 case(3)
   !
   ! regulate both norm and gradient 
   ! norm => remove the null space 
   ! gradient => enforce smoothness
   !
   do isp=1,nspin
      call laplacian(nfft,cell_nfft(1),cell_nfft(2),cell_nfft(3),x(:,isp),qvec,lap)
      Ax(:,isp) = Ax(:,isp) - reg_lambda*lap + reg_lambda*x(:,isp)
   enddo
 endselect


endsubroutine regularization
