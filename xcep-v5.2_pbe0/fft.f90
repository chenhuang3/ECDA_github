!=================================================
! forward/backward FFT, using FFTW3 libraray
!  d: the incoming/output real array
!  c: the outcoming/input complex array 
!  is: 1  forward
!     -1  backward
!
!  Note: (1) FFTW3 does not do normalization in FFT
!            I did it manually in the backward transformation part
!       
!        (2) The forward FFT here gives the same results
!            as the MATLAB's fftn(X) subroutine
! 
!  Chen Huang, Jan/2011
!
!
! Note on using FFTW:
!
!   (1) your real(kind=8) array is defined as L x N x M
!       the FORTRAN interface from FFTW will do 
!       the conversion from column-major to row-major
!       for you, so just use your FOTRAN-style array as input.
!
!   (2) The transformed array must be defined as 
!       (L/2+1) x N x M, the 1st dim is halved then add one, 
!       to save memory in FFTW.
!
!   (3) The FFT subroutine in this code, can only do 
!       Real => Complex and Complex => Real.
!
subroutine  fft(n1,n2,n3,din,cin,is)

  implicit none

  include 'fftw3.f'

  ! external vars ---------

  integer,intent(in)   :: is                              ! is =1, forward FFT, is=-1, backward
  integer,intent(in)   :: n1,n2,n3                        ! dimensions
  real(kind=8)   ,intent(inout) :: din(n1,n2,n3)         ! the data REAL to be transformed
  complex(kind=8),intent(inout) :: cin((n1/2+1),n2,n3)   ! real part is in the 1st dim,imaginary in 2nd

  ! local vars ---------

  integer*8        :: plan                ! integer*8 is fftw_plan type, see http://www.fftw.org/doc/Fortran-Examples.html
  integer          :: dim1                ! the halved 1st dimension
  double precision :: d(n1,n2,n3)
  double complex   :: c(n1/2+1,n2,n3)

  ! check consistence ----------

!  if ( mod(n1,2)/=0 ) then
!    print *,' n1 must be even number for using FFTW, stop!'
!    stop
!  endif

  d = din
  c = cin
  
  ! fourier-transform --------------

  dim1 = n1/2+1         ! to work with FFTW, which saves memory by 
                        ! halve the first dimension and add 1
  if (is==1) then 

    ! without the FFTW_ESTIMATE flag, the d and c will 
    ! be overwritten during the planning

    call dfftw_plan_dft_r2c_3d(plan,n1,n2,n3,d,c,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan,d,c)

  else if (is==-1) then

    ! backward fourier transform -------------

    ! without the FFTW_ESTIMATE flag, the d3d and c3d will 
    ! be overwritten during the planning

    call dfftw_plan_dft_c2r_3d(plan,n1,n2,n3,c,d,FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan,c,d)

    ! FFTW DOES NOT PERFORM NORMALIZATION
    ! DO THE NORMALIZATION 
    d = d/dble(n1*n2*n3)

  else
    print *,'undefined IS in FFT(), error in fft, stop '
    stop
  endif

  din = d
  cin = c

  call dfftw_destroy_plan(plan)
  return 


end subroutine fft




