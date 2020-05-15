! compute L-BFGS direction  
! 
!  mhist:  number of historical iterations 
!  x(:,1:m+1):   store m+1 previous x
!  g(:,1:m+1):   store m+1 previous g
!  x0: current x
!  g0: current g
!
! output: 
!  q: L-BFGS direction 
!  x: histroy x are updated on exit 
!  g: history g are updated on exit 
!
! based on 
!   https://github.com/GuipengLi/optLBFGS/blob/master/optLBFGS.m
!
!   wikipedia (seems to have a bug, but helps me understand L-BFGS)
!
!   See Alg. 3 in "A Stochastic Quasi-Newton Method for Online Convex Optimization", N.N. Schraudolph et al. 
!     Proceedings of Machine Learning Research, 2:436-443, (2007)
!     http://proceedings.mlr.press/v2/schraudolph07a.html
!
! Created on 1/21/2019 by Chen Huang 
!
subroutine  lbfgs_dir(nspin,nfft,mhist,iter,g0,x0,g,x,q)

  implicit none 

  integer, intent(in) :: nfft, mhist, iter, nspin
  real(8)             :: g(nspin*nfft,mhist+1), & 
                         x(nspin*nfft,mhist+1), & 
                         g0(nspin*nfft), & 
                         x0(nspin*nfft)

  real(8), intent(out) :: q(nspin*nfft) 


  ! local vars 
  integer :: j, i, m
  real(8) :: s(nspin*nfft), & 
             y(nspin*nfft), &
             rho, alpha(mhist), beta


  ! >>>> FUNCTION BEGINS <<<<

  ! following wikipedia 
  ! https://en.wikipedia.org/wiki/Limited-memory_BFGS#Algorithm

  ! update historical x and g
  ! 1,2,...,m are historical x and g
  do i=1,mhist
    x(:,i) = x(:,i+1)
    g(:,i) = g(:,i+1)
  enddo
  ! m+1 is the current x and g
  x(:,mhist+1) = x0
  g(:,mhist+1) = g0


  if (mhist >= iter-1) then 
    m = iter-1
  else 
    m = mhist
  endif 


  if (m==0) then 
    q = -g0
    return 
  endif 

  q = -g0

  ! first loop
  do i = 1, m
    s = x(:,mhist-i+2)-x(:,mhist-i+1)  ! s_i
    y = g(:,mhist-i+2)-g(:,mhist-i+1)  ! y_i
    rho = 1.d0/sum(y*s)      ! define rho 
    alpha(i) = rho*sum(s*q)  ! define alpha
    q = q - alpha(i)*y          
  enddo 

  s = x0 - x(:,mhist)
  y = g0 - g(:,mhist)
  q = sum(s*y)/sum(y*y)*q

  ! second loop
  do i = m,1,-1
    s = x(:,mhist-i+2)-x(:,mhist-i+1)  ! s_i
    y = g(:,mhist-i+2)-g(:,mhist-i+1)  ! y_i
    rho = 1.d0/sum(y*s)      ! define rho 
    beta = rho*sum(y*q)     ! define beta 
    q = q + s*(alpha(i)-beta)
  enddo 


end subroutine lbfgs_dir
