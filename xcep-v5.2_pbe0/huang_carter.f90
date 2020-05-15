!============================================
! implement of the Huang-Carter 2010 KEDF
! fully interpolated version 
!
!   Chen Huang
!    Feb/2011
!
!============================================
subroutine  huang_carter(nfft,n1,n2,n3,rho,qvec,dvol,pot,ke)

  implicit none

  ! external vars
  !=================
  integer,intent(in)  :: nfft,n1,n2,n3
  real(kind=8)        :: rho(nfft),qvec(3,(n1/2+1)*n2*n3),dvol

  real(kind=8), intent(out) :: ke
  real(kind=8), intent(out) :: pot(nfft)

  ! local vars
  !=============
  integer,parameter :: numEta=50001
  real(kind=8), parameter :: pi  = 3.1415926d0, &
                       cTF = 2.87123400018819d0 , & 
                       const = 8.d0*3.d0*pi**2, & 
                       ratio = 1.2d0, & 
                       rhoS  = 1.e-2, & 
                       beta  = 0.51d0

  integer :: ii, &
             nref, &           ! size of refxi array
             upperpower, & 
             downpower, & 
             allocatestat

  real(kind=8) :: midxi, maxxi, minxi

  real(kind=8),dimension(n1,n2,n3) :: & 
    rho3d, xi

  real(kind=8),dimension(3,n1/2+1,n2,n3) :: & 
    qvec3d

  real(kind=8),allocatable :: &  
    refxi(:)
  

  ! function begins
  !================
  print *,''
  print *,' ---- enter huang_carter() -----'
  print *,' beta   =',beta
  print *,' ratio  =',ratio

!! debug
!  open(file='dump_den',action='read',unit=111)
!  read(111,*)rho
!  print *,'load dump_rho'
!  print *,'max/min of rho=',minval(rho),maxval(rho)
!  close(111)
!! end debug

  rho3d = reshape(rho,(/n1,n2,n3/))
  qvec3d  = reshape(qvec,(/3,n1/2+1,n2,n3/))
  xi = (3.d0*pi**2*rho3d)**0.33333333333d0

  ! make refxi array
  !==================
  maxxi = maxval(xi)
  minxi = minval(xi)
  midXi = (3.d0*pi**2.d0*rhoS)**(1.d0/3.d0) 
  upperPower = ceiling(LOG(maxXi/midXi) / LOG(ratio))
  downPower  = floor  (LOG(minXi/midXi) / LOG(ratio))
  nref = upperpower-downpower+1
  print *,'nref = ',nref
  allocate (refxi(nref), stat = allocatestat)
  if (allocatestat /= 0) then
    write(*,*) 'error in allocating rhobnd() in dividedensity (), stop'
    stop
  endif
  do ii = upperPower - downPower + 1, 1, -1
    refxi(ii) = midxi * ratio**(downPower + ii-1)
  enddo

  ! compute kinetic energy
  !========================
  call ke_HC01(n1,n2,n3,rho3d,xi,nref,refxi,qvec3d,beta,dvol,ke)
  
  ! compute kinetic energy potential
  !=================================
  call kepot_HC01(n1,n2,n3,rho3d,xi,nref,refxi,qvec3d,beta,pot)

  deallocate(refxi)

  return
end subroutine  huang_carter


!=======================================
! compute HC01 kinetic energy 
! 
!  Chen Huang
!   Feb/2011
!=======================================
subroutine ke_HC01(n1,n2,n3,rho,xi,nref,refxi,qvec3d,beta,dvol,ke)

  use mathfunc, only: & 
    stepfun, & 
    h00,h10,h01,h11

  implicit none
 
  ! external vars
  !================
  integer                 :: n1,n2,n3,nref
  real(kind=8),intent(in) :: rho(n1,n2,n3),xi(n1,n2,n3),beta,refxi(nref)
  real(kind=8),intent(in) :: qvec3d(n1/2+1,n2,n3)
  real(kind=8),intent(in) :: dvol
  real(kind=8),intent(out) :: ke
  
  ! local vars
  !============
  integer :: jj
  real(kind=8),dimension(n1,n2,n3) :: & 
    R_alp, & 
    tmp3d, & 
    mask,t

  complex(kind=8),dimension(n1/2+1,n2,n3) :: & 
    FFTx, &
    ffttmp, &
    kernel0, kernel1, & 
    kernel0a,kernel1a

  real(kind=8) :: & 
    alpha, & 
    cTF = 2.87123400018819d0 , & 
    c   = 8.d0 * 3.d0 * 3.14159265d0**2 ,& 
    h
   
  !=================
  ! function begins
  !=================

  print *,'enter ke_HC01...'
  
  alpha = 8.d0/3.d0 - beta    ! adsorb (2kF)**3 to Rho_alp
  R_alp = rho**alpha/(2.d0*xi)**3

  FFTx  = 0.d0
  DO jj = 1, nref-1
    
    call make_HC01_kernel(refxi(jj),  n1,n2,n3,qvec3d,kernel0,kernel1)
    call make_HC01_kernel(refxi(jj+1),n1,n2,n3,qvec3d,kernel0a,kernel1a)

    h = refXi(jj+1) - refXi(jj)
    t = (xi - refXi(jj)) / h
    mask = stepfun(n1,n2,n3,xi-refXi(jj)) * stepfun(n1,n2,n3,refxi(jj+1)-xi)

    call fft(n1,n2,n3,R_alp*h00(n1,n2,n3,t)*mask,ffttmp,1)
    FFTx = FFTx + kernel0 * ffttmp

    call fft(n1,n2,n3,R_alp*h10(n1,n2,n3,t)*mask,ffttmp,1)
    FFTx = FFTx + kernel1 * h * ffttmp

    call fft(n1,n2,n3,R_alp*h01(n1,n2,n3,t)*mask,ffttmp,1)
    FFTx = FFTx + kernel0a * ffttmp

    call fft(n1,n2,n3,R_alp*h11(n1,n2,n3,t)*mask,ffttmp,1)
    FFTx = FFTx + kernel1a * h * ffttmp

!    FFTx = FFTx &
!       + kernel0(:,:,:) * FFT( R_alp*h00(t)*mask )      &
!       + kernel1(:,:,:) * h * FFT( R_alp*h10(t)*mask )  &
!       + kernel0a(:,:,:)* FFT( R_alp*h01(t)*mask )      &
!       + kernel1a(:,:,:)* h * FFT( R_alp*h11(t)*mask ) 
  ENDDO

  call fft(n1,n2,n3,tmp3d,FFTx,-1)
   
  ke = cTF * C * SUM( rho**beta * tmp3d ) * dvol

  print *, 'HC01 nonlocal part =',ke
!  stop

  RETURN
end subroutine  ke_HC01



!=======================================
! compute HC01 kinetic energy potential
!
!  Chen Huang
!   Feb/2011
!=======================================
subroutine kepot_HC01(n1,n2,n3,rho,xi,nref,refxi,qvec3d,beta,kepot)

  use mathfunc, only : & 
    stepfun , & 
    h00, h01, & 
    h10, h11

  implicit none

  ! external vars
  !================
  integer       :: n1,n2,n3,nref
  real(kind=8)  :: & 
    beta, &
    refxi(nref)
  real(kind=8),dimension(n1,n2,n3)  :: & 
    rho,xi,kepot
  real(kind=8),dimension(3,n1/2+1,n2,n3) :: & 
    qvec3d
  
  ! local vars
  !==============
  integer :: jj 

  real(kind=8) :: & 
    cTF = 2.87123400018819d0 , & 
    c   = 8.d0 * 3.d0 * 3.14159265d0**2 

  real(kind=8) :: & 
    gradrho(n1,n2,n3,3), & 
    h, & 
    alpha

  real(kind=8),dimension(n1,n2,n3) :: & 
    Q_term,N_term , & 
    R_alp,  & 
    t,t2,& 
    mask, & 
    rho_alp_int, &
    tmppot,      & 
    tmp3d, &
    h00p,h10p,h01p,h11p
    
  complex(kind=8) :: & 
    imag
  
  complex(kind=8),dimension(n1/2+1,n2,n3) :: &
    FFTtmp, & 
    FFTx, FFTb,  & 
    kernel0,kernel0a, & 
    kernel1,kernel1a
  
  print *,'enter kepot_HC01 ...'
  print *,'nref = ',nref
!  print *,'refxi',refxi

!! debug
!  open(file='dump_den',action='read',unit=111)
!  read(111,*)rho
!  print *,'load dump_rho'
!  print *,'max/min of rho=',minval(rho),maxval(rho)
!  close(111)
!! end debug
 
  ! grad(rho)
  !============
  imag = cmplx(0.d0,1.d0)

  call fft(n1,n2,n3,rho,FFTtmp,1)
  call fft(n1,n2,n3,gradrho(:,:,:,1),imag*qvec3d(1,:,:,:)*FFTtmp,-1)
  call fft(n1,n2,n3,gradrho(:,:,:,2),imag*qvec3d(2,:,:,:)*FFTtmp,-1)
  call fft(n1,n2,n3,gradrho(:,:,:,3),imag*qvec3d(3,:,:,:)*FFTtmp,-1)

!  gradRho(:,:,:,1) = FFT(imag*qVectors(:,:,:,1)*FFTtmp)
!  gradRho(:,:,:,2) = FFT(imag*qVectors(:,:,:,2)*FFTtmp)
!  gradRho(:,:,:,3) = FFT(imag*qVectors(:,:,:,3)*FFTtmp)

  alpha = 8.d0/3.d0-beta
  R_alp = rho**alpha/(2.d0*xi)**3
  call fft(n1,n2,n3,rho**beta,FFTb,1)
!  FFTb  = FFT(rho**beta)
  FFTx  = 0.d0
  Q_term= 0.d0
  N_term= 0.d0

  !===== loop over bins ==========
  DO jj = 1, nref-1

!!    print *,'jj=',jj

    call make_HC01_kernel(refxi(jj),  n1,n2,n3,qvec3d,kernel0, kernel1)
    call make_HC01_kernel(refxi(jj+1),n1,n2,n3,qvec3d,kernel0a,kernel1a)
!    CALL MakeKernel(refXi(jj),kernel0,kernel1)
!    CALL MakeKernel(refXi(jj+1),kernel0a,kernel1a)

    h = refXi(jj+1) - refXi(jj)
    t = (xi - refXi(jj)) / h
    mask = stepfun(n1,n2,n3,xi-refXi(jj)) * stepfun(n1,n2,n3,refXi(jj+1)-xi)

    !============ PART I: term Q in notes (on wiki) ==========
    call fft(n1,n2,n3,tmp3d,kernel0*FFTb,-1)
    Q_term = Q_term + mask * h00(n1,n2,n3,t)*tmp3d

    call fft(n1,n2,n3,tmp3d,kernel1*FFTb,-1)
    Q_term = Q_term + mask * h * h10(n1,n2,n3,t)*tmp3d

    call fft(n1,n2,n3,tmp3d,kernel0a*FFTb,-1)
    Q_term = Q_term + mask * h01(n1,n2,n3,t)*tmp3d

    call fft(n1,n2,n3,tmp3d,kernel1a*FFTb,-1)
    Q_term = Q_term + mask * h * h11(n1,n2,n3,t)*tmp3d

!    Q_term = Q_term +  mask * & 
!     ( h00(n1,n2,n3,t) * FFT(kernel0(:,:,:)  * FFTb) + &
!       h10(n1,n2,n3,t) * FFT(kernel1(:,:,:)  * h * FFTb) + & 
!       h01(n1,n2,n3,t) * FFT(kernel0a(:,:,:) * FFTb)     + & 
!       h11(n1,n2,n3,t) * FFT(kernel1a(:,:,:) * h * FFTb )) 

    !============ PART II: term M in notes (on wiki) ============
    call fft(n1,n2,n3,R_alp*h00(n1,n2,n3,t)*mask,ffttmp,1)
    FFTx = FFTx + kernel0 * ffttmp

    call fft(n1,n2,n3,R_alp*h10(n1,n2,n3,t)*mask,ffttmp,1)
    FFTx = FFTx + kernel1 * h * ffttmp

    call fft(n1,n2,n3,R_alp*h01(n1,n2,n3,t)*mask,ffttmp,1)
    FFTx = FFTx + kernel0a * ffttmp

    call fft(n1,n2,n3,R_alp*h11(n1,n2,n3,t)*mask,ffttmp,1)
    FFTx = FFTx + kernel1a * h * ffttmp
!    FFTx = FFTx &
!      + kernel0(:,:,:) * FFT( R_alp * h00(t) * mask )                  &
!      + kernel1(:,:,:) * h * FFT( R_alp * h10(t) * mask )  &
!      + kernel0a(:,:,:)* FFT( R_alp * h01(t) * mask )                  &
!      + kernel1a(:,:,:)* h * FFT( R_alp * h11(t) * mask ) 

    !=========== PART III: term N in notes (on wiki) ==========
    ! d (h(t))/ dt
    t2 = t**2
    h00p = 6.d0*t2 - 6.d0*t
    h10p = 3.d0*t2 - 4.d0*t + 1.d0
    h01p = -6.d0*t2 + 6.d0*t
    h11p = 3.d0*t2 - 2.d0*t
 
    call fft(n1,n2,n3,tmp3d,kernel0*FFTb,-1)
    N_term = N_term + mask * 1.d0/h * h00p * tmp3d

    call fft(n1,n2,n3,tmp3d,kernel1*FFTb,-1)
    N_term = N_term + mask * h10p * tmp3d

    call fft(n1,n2,n3,tmp3d,kernel0a*FFTb,-1)
    N_term = N_term + mask * 1.d0/h * h01p * tmp3d

    call fft(n1,n2,n3,tmp3d,kernel1a*FFTb,-1)
    N_term = N_term + mask * h11p * tmp3d
    
!    N_term = N_term + mask * & 
!      ( 1.d0/h * h00p * FFT(kernel0(:,:,:)*FFTb) + & 
!        h10p * FFT(kernel1(:,:,:)*FFTb) + & 
!        1.d0/h * h01p * FFT(kernel0a(:,:,:)*FFTb) + & 
!        h11p * FFT(kernel1a(:,:,:)*FFTb) )

  ENDDO ! jj

  !============= pot1 in notes (on wiki)======
  kepot = c*cTF*alpha*rho**(alpha-1)/(2.d0*xi)**3*Q_term + & 
          c*cTF*(-3.d0/8.d0)*rho**alpha/xi**4*xi/(3.d0*rho)*Q_term 

  !============= pot2 in notes (on wiki)=======
  call fft(n1,n2,n3,rho_alp_int,FFTx,-1)
  
  ! remove instability where rho < 1e-5
  tmppot = rho**(beta-1.d0)
  where (rho<1e-5) 
    tmppot = 1e-5**(beta-1.d0)
  endwhere
  tmppot = c*cTF*beta*tmppot*rho_alp_int

  kepot  = kepot + tmppot

  !============= pot3 in notes (on wiki)=======
  tmppot = c*cTF*R_alp*xi/(3.d0*rho)*N_term 
  kepot  = kepot + tmppot

  write(6,'(a,2es12.4)')'leave kepot_HC01. kepot=>',minval(kepot),maxval(kepot)

!  print *,'kepot =>',minval(kepot),maxval(kepot)
!  stop
  return
end subroutine kepot_HC01



!============================================
! prepare kernenls with given fermi vector
!
!  Chen Huang
!   Feb/2011
!
!============================================
subroutine make_HC01_kernel(kF,n1,n2,n3,qvec3d,k1,k2)
  implicit none
  
  ! external vars
  !===============
  integer,intent(in) :: n1,n2,n3
  real(kind=8), intent(in) :: kF
  real(kind=8), intent(in) :: qvec3d(3,n1/2+1,n2,n3)
  complex(kind=8), intent(out) :: & 
       k1(n1/2+1,n2,n3), &      ! F[(2kF)**3*w(r)]
       k2(n1/2+1,n2,n3)         ! F[ d{(2kF)**3*w(r)}/d{deltaR}*deltaR ]

  ! local vars
  !=============
  integer,parameter :: numEta = 50001
  integer :: ii,jj,kk
  integer :: left_index,right_index
  real(kind=8) :: w(numEta,2),wp(numEta),eta_int
  real(kind=8) :: q0,eta,coeff
   
  ! function begins
  !================

  ! make kernel w(eta)
  ! we just read it from file: huang_carter_kernel.dat
  !===================================================
!  print *,'reading HC01.dat file ...'
  open(file='HC01.dat',unit=111,action='read')
  do ii=1,numEta
    read(111,*) w(ii,1),w(ii,2),wp(ii)  ! eta and w(eta)
    wp(ii) = wp(ii) * w(ii,1)           ! wp = w'*eta
  enddo
  close(111)
  eta_int = w(numEta,1)/dble((size(w,1))-1)
!  print *,'max eta  = ',w(numEta,1)
!  print *,'eta_int  = ',eta_int
!  print *,'# of eta = ',size(w,1)

  ! make kernel3d array in q space
  !=================================
!  print *,'making kernel3d ...'
  k1 = -100.d0
  k2 = -100.d0
  do kk=1,n3
    do jj=1,n2
      do ii=1,n1/2+1
        q0 = sqrt(sum(qvec3d(:,ii,jj,kk)**2,1))
        eta = q0 / (2.d0*kF)
        if (eta > w(size(w,1),1)) then 
          k1(ii,jj,kk) = w(size(w,1),2)
          k2(ii,jj,kk) = wp(size(w,1))
        else
          left_index  = floor(eta/eta_int) + 1
          right_index = left_index + 1
          coeff       = eta/eta_int - floor(eta/eta_int)
          ! linear interpolation , must do interpolation for q=0 point
          k1(ii,jj,kk) = cmplx( w(left_index,2)*(1.d0-coeff)+coeff*w(right_index,2),0.d0)
          k2(ii,jj,kk) = cmplx(-(wp(left_index)*(1.d0-coeff)+coeff*wp(right_index))/kF,0.d0)
        endif
!        print *,q0,eta,real(kernel3d(ii,jj,kk))
      enddo
    enddo
  enddo
!  print *,'k1 max/min=',maxval(real(k1)),minval(real(k1))
!  print *,'k2 max/min=',maxval(real(k2)),minval(real(k2))
!  print *,'done kernel3d .'

  return

end subroutine make_HC01_kernel

