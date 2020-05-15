 !---------------------------------------------------------------------------
 !
 ! The A matrix is A = I - |1-weight><dfermi/dvks|
 ! Note that A is not symmetric, therefore CG method cannot be used. 
 ! We solve A^T x = b using the CGLS method.
 ! Also, note the transform of matrix A above.
 ! the matrix can be non-symmetric. The details can be found here: 
 !
 ! http://web.stanford.edu/group/SOL/software/cgls/
 !
 ! The code following 
 ! http://web.stanford.edu/group/SOL/software/cgls/matlab/cgls.m
 !
 !---------------------------------------------------------------------------
 subroutine apply_Ainv(nfft,ucvol,trans_A,weight,dfermi_dvks,b,x)
 
  use mpi

  implicit none

  logical :: trans_A   ! true:  solve A'x=b
                       ! false: solve  Ax=b
  integer :: nfft, iter_cg
  real(8) :: ucvol, & 
             weight(nfft),  & 
             dfermi_dvks(nfft), & 
             b(nfft), &  ! the b vector in Ax=b
             x(nfft), & 
             cg_r(nfft), & 
             Av(nfft), & 
             cg_p(nfft), &
             cg_q(nfft), &
             cg_s(nfft), & 
             gamma, gamma1, & 
             reg_lambda = 0.0d0, & 
             dvol, delta, alpha, beta, dtmp

  ! MPI
  integer :: ierr, myrank


  ! >>> function begins <<< !

  call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, ierr )
  dvol = ucvol / dble(nfft)


  if (myrank==0) then 
    print *,'-------- apply_Ainv() -------'
    print *,'solve Ax=b with b: ',minval(b),maxval(b)
    print *,'dfermi/dvks: ',minval(dfermi_dvks),maxval(dfermi_dvks)
    print *,'weight:      ',minval(weight),maxval(weight)
  endif 

  x = 0.0d0
  cg_r = b ! r = b - Ax

  ! get s = A'*r
  if (trans_A) then 
    Av = cg_r - (1.d0-weight)*sum(dfermi_dvks*cg_r)*dvol
  else 
    Av = cg_r - dfermi_dvks*sum((1.d0-weight)*cg_r)*dvol
  endif 
  cg_s = Av - reg_lambda*x 

  ! set p=s
  cg_p  = cg_s
  gamma = sum(cg_s**2)*dvol

  iter_cg = 0

  !-----------------------------------------------------
  ! CGLS (conjugate gradient for least squared method)
  !-----------------------------------------------------
  do while (.true.) 

    ! converge test 
    dtmp = sum(cg_r**2)*dvol
    if ( dtmp<1e-16 ) then 
      if (myrank==0) then 
        write(*,'(a)')'apply_Ainv is done.'
        write(*,*)'final x: ',minval(x),maxval(x)
        ! compute residual directly Ax-b
        if (trans_A) then 
           Av = x - dfermi_dvks*sum((1.d0-weight)*x)*dvol
        else 
           Av = x - (1.d0-weight)*sum(dfermi_dvks*x)*dvol
        endif 
        dtmp = dvol*sum((Av-b)**2)
        write(*,*)'direct check resid: ',dtmp
      endif
      exit
    endif 

    iter_cg = iter_cg +  1

    ! compute q=A*p
    if (trans_A) then 
      cg_q = cg_p - dfermi_dvks*sum((1.d0-weight)*cg_p)*dvol
    else
      cg_q = cg_p - (1.d0-weight)*sum(dfermi_dvks*cg_p)*dvol
    endif 
    delta = sum(cg_q**2)*dvol + reg_lambda*sum(cg_p**2)*dvol;
    alpha = gamma / delta;

    x    = x    + alpha*cg_p
    cg_r = cg_r - alpha*cg_q 

    ! compute s = A'*r - shift*x
    if (trans_A) then 
      Av = cg_r - (1.d0-weight)*sum(dfermi_dvks*cg_r)*dvol
    else 
      Av = cg_r - dfermi_dvks*sum((1.d0-weight)*cg_r)*dvol
    endif 
    cg_s = Av - reg_lambda*x

    gamma1 = gamma
    gamma  = sum(cg_s**2)*dvol
    beta   = gamma / gamma1

    ! update p = s + beta*p
    cg_p   = cg_s + beta * cg_p

    dtmp = sum(cg_r**2)*dvol

    if (myrank==0) then
      write(*,'(a,i3,a,es8.2,a,es10.2,a,es10.2,a,2f12.4)') & 
       "iter: ",iter_cg,"  |resid|: ",dtmp,' alpha:',alpha,' beta: ',beta,' x:',minval(x),maxval(x)
    endif


    if (iter_cg>10000) then 
      print *,'after 10000 iterations in apply_Ainv.f90, CGLS cannot converge. stop, for rank: ',myrank
      stop
    endif 
  enddo ! end of CG of solving global_vxc

  !!stop
 end subroutine apply_Ainv



