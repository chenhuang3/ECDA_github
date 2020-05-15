program posdopspectra

!! Script to analyze the electron-positron momentum distributions
!! included in the _DOPPLER file. Calculated 1D distributions in
!! (001), (011) and (111) directions, convolutes the data with
!! a Gaussian and interpolates them on an uniform grid with
!! 0.1 mrad spacing. Calculates the S and W parameters 
!! of the Doppler broadening of annihilation radiation.

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
!scalars
 integer, parameter:: dp=kind(0.d0) 
 integer :: i1,i2,i3,ig1,ig2,ig3,ii,ikpt,iz
 integer :: nfft,nkpt,nd1,nd2,nd3,nd4,nd5,nd6,np1,np2,np3,nspline
 real(dp) :: der1,der2,der3,der4,der5,der6,fwhm,lambda,S_dop
 real(dp) :: zero,one,two,three,two_pi,InvFineStruct
 real(dp) :: vec,vecg,ucvol,W_dop
 logical :: file_exists
 character(len=264) :: filnam
!arrays
 real(dp),allocatable :: der001(:),der011(:),der111(:),gauss(:)
 real(dp),allocatable :: mesh_spline(:),mesh001(:),mesh011(:),mesh111(:)
 real(dp),allocatable :: pcart(:,:,:),rho_moment(:,:),rho_001(:),rho001(:),rho001conv(:)
 real(dp),allocatable :: rho_011(:),rho011(:),rho011conv(:),rho_111(:),rho111(:),rho111conv(:)
 character(len=500) :: msg

!******************************************************************
 zero=0._dp;one=1._dp;two=2._dp;three=3._dp
 two_pi=8*atan(1.d0)
 InvFineStruct=137.035999679_dp  ! Inverse of fine structure constant
! Get name of density file
 write(*,*)
 write(*,*) ' What is the name of the 3D electron-positron momentum distribution file?'
 read(*,'(a)')filnam
 write(*,*)
!  Checking the existence of data file
 INQUIRE(FILE=filnam, EXIST=file_exists)   ! file_exists will be TRUE if the file
 if (.not.file_exists) then
   write(*,*) 'Missing data file: '//TRIM(filnam)
   stop
 end if

!  if (open(unit=19,file=filnam,form='unformatted',status='old') /=0) then
!    MSG_ERROR(msg)
!  end if 
 write(*,'(a,a,a,i4)')'  posdopspectra : read file ',trim(filnam)
 write(*,*)
 write(*,'(a,a)') 'Opening file ', filnam
 open(unit=19,file=filnam,form='unformatted',status='old')
 read(19) nfft,nkpt,np1,np2,np3,vec,ucvol
 ALLOCATE(pcart(3,nfft,nkpt))
 ALLOCATE(rho_moment(nfft,nkpt))
 do ikpt=1,nkpt
   read(19) pcart(1:3,1:nfft,ikpt),rho_moment(1:nfft,ikpt)
 end do
 write(*,*) ' => Your momentum distribution file is : ',trim(filnam)
 write(*,*) ' Choose FWHM (in mrad) for the convolution.'
 read(*,*) fwhm
! Normalize the momentum distibution to one (sum(rho001)*vec=1)
! and transform rho_moment to 10^3(m_0c)^-1

!  lambda=sum(rho_moment(1:nfft,1:nkpt))*(nkpt**(one/three))*(two_pi)**three/(ucvol)
!  rho_moment(1:nfft,1:nkpt)=rho_moment(1:nfft,1:nkpt)*InvFineStruct/(lambda*1000_dp)

 nd1=-(np1+1)/2-(np2+1)/2-(np3+1)/2+3;nd2=np1/2+np2/2+np3/2
 nd3=          -(np2+1)/2-(np3+1)/2+2;nd4=      np2/2+np3/2
 nd5=                    -(np3+1)/2+1;nd6=            np3/2

 ALLOCATE(rho111(nd1:nd2))
 ALLOCATE(rho011(nd3:nd4))
 ALLOCATE(rho001(nd5:nd6))
 ALLOCATE(mesh111(nd1:nd2))
 ALLOCATE(mesh011(nd3:nd4))
 ALLOCATE(mesh001(nd5:nd6))

 do iz=nd1,nd2
   mesh111(iz)=iz*vec*two_pi*1000_dp/(InvFineStruct*sqrt(three))
 end do
 do iz=nd3,nd4
   mesh011(iz)=iz*vec*two_pi*1000_dp/(InvFineStruct*sqrt(two))
 end do
 do iz=nd5,nd6
   mesh001(iz)=iz*vec*two_pi*1000_dp/InvFineStruct
 end do
 rho111(:)=zero;rho011(:)=zero;rho001(:)=zero
! Compute rho along (111), (011) and (001) axis through integration

 do ikpt=1,nkpt
   do ii=1,nfft
     ig1=nint(pcart(1,ii,ikpt)/vec)
     ig2=nint(pcart(2,ii,ikpt)/vec)
     ig3=nint(pcart(3,ii,ikpt)/vec)
!      vecg=sqrt(pcart(1,ii,ikpt)**2+pcart(2,ii,ikpt)**2+pcart(3,ii,ikpt)**2)*two_pi*1000_dp/InvFineStruct
!      if (vecg<70) then
     rho111(ig1+ig2+ig3)=rho111(ig1+ig2+ig3)+rho_moment(ii,ikpt)
     rho011(    ig2+ig3)=rho011(    ig2+ig3)+rho_moment(ii,ikpt)
     rho001(        ig3)=rho001(        ig3)+rho_moment(ii,ikpt)
!      end if
   end do
 end do
 DEALLOCATE(pcart)
 DEALLOCATE(rho_moment)
 rho111(:)=rho111(:)*sqrt(three)*((two_pi)**two)/((ucvol/real(nkpt))**(two/three))
 rho011(:)=rho011(:)*sqrt(two  )*((two_pi)**two)/((ucvol/real(nkpt))**(two/three))
 rho001(:)=rho001(:)*(            (two_pi)**two)/((ucvol/real(nkpt))**(two/three))

! Convolve the 1d spectrum with a Gaussian (now FWHM=3,71 10^-3 m0c)

 ALLOCATE(rho111conv(nd1:nd2))
 ALLOCATE(rho011conv(nd3:nd4))
 ALLOCATE(rho001conv(nd5:nd6))
 rho111conv=zero;rho011conv=zero; rho001conv=zero
 ALLOCATE(gauss(nd1:nd2))
 do iz=nd1,nd2
   gauss(iz)=exp(-dble(mesh111(iz)*mesh111(iz))/(two*((fwhm/2.35482_dp)**2)))
 end do
 do iz=nd1,0
   rho111conv(iz) = zero
   i1=nd1
   i2=iz-nd1
   do while (i2 >= nd1)
     rho111conv(iz) = rho111conv(iz) + rho111(i1)*gauss(i2)
     i1 = i1+1
     i2 = i2-1
   end do
 end do
 do iz=1,nd2
   rho111conv(iz) = zero
   i1=iz-nd2
   i2=nd2
   do while (i1 <= nd2)
     rho111conv(iz) = rho111conv(iz) + rho111(i1)*gauss(i2)
     i1 = i1+1
     i2 = i2-1
   end do
 end do
! Normalize to one
 lambda=sum(rho111conv(0:nd2))*mesh111(1)
 rho111conv(:)=rho111conv(:)/lambda

 gauss=zero

 do iz=nd3,nd4
   gauss(iz)=exp(-dble(mesh011(iz)*mesh011(iz))/(two*((fwhm/2.35482_dp)**2)))
 end do
 do iz=nd3,0
   rho011conv(iz) = zero
   i1=nd3
   i2=iz-nd3
   do while (i2 >= nd3)
     rho011conv(iz) = rho011conv(iz) + rho011(i1)*gauss(i2)
     i1 = i1+1
     i2 = i2-1
   end do
 end do
 do iz=1,nd4
   rho011conv(iz) = zero
   i1=iz-nd4
   i2=nd4
   do while (i1 <= nd4)
     rho011conv(iz) = rho011conv(iz) + rho011(i1)*gauss(i2)
     i1 = i1+1
     i2 = i2-1
   end do
 end do
! Normalize to one
 lambda=sum(rho011conv(0:nd4))*mesh011(1)
 rho011conv(:)=rho011conv(:)/lambda

 gauss=zero

 do iz=nd5,nd6
   gauss(iz)=exp(-dble(mesh001(iz)*mesh001(iz))/(two*((fwhm/2.35482_dp)**2)))
 end do
 do iz=nd5,0
   rho001conv(iz) = zero
   i1=nd5
   i2=iz-nd5
   do while (i2 >= nd5)
     rho001conv(iz) = rho001conv(iz) + rho001(i1)*gauss(i2)
     i1 = i1+1
     i2 = i2-1
   end do
 end do
 do iz=1,nd6
   rho001conv(iz) = zero
   i1=iz-nd6
   i2=nd6
   do while (i1 <= nd6)
     rho001conv(iz) = rho001conv(iz) + rho001(i1)*gauss(i2)
     i1 = i1+1
     i2 = i2-1
   end do
 end do
! Normalize to one
 lambda=sum(rho001conv(0:nd6))*mesh001(1)
 rho001conv(:)=rho001conv(:)/lambda

 DEALLOCATE(gauss)

 open(unit=114, file='doppler_out', status='replace')

 write(114,*)'Raw momentum distribution in (001) direction:'
 write(114,*)
 write(114,*)'  p001 (mrad):             Probalility:'
 do iz=0,nd6
   write(114,*) mesh001(iz),rho001(iz)
 end do
 write(114,*)
 write(114,*)'Raw momentum distribution in (011) direction:'
 write(114,*)
 write(114,*)'  p011 (mrad):             Probalility:'
 do iz=0,nd4
   write(114,*) mesh011(iz),rho011(iz)
 end do
 write(114,*)
 write(114,*)'Raw momentum distribution in (111) direction:'
 write(114,*)
 write(114,*)'  p111 (mrad):             Probalility:'
 do iz=0,nd2
   write(114,*) mesh111(iz),rho111(iz)
 end do
 write(114,*)
 write(114,*)'Convoluted and normalized momentum distribution in (001) direction:'
 write(114,*)
 write(114,*)'  p001 (mrad):             Probalility:'
 do iz=0,nd6
   write(114,*) mesh001(iz),rho001conv(iz)
 end do
 write(114,*)
 write(114,*)'Convoluted and normalized momentum distribution in (011) direction:'
 write(114,*)
 write(114,*)'  p011 (mrad):             Probalility:'
 do iz=0,nd4
   write(114,*) mesh011(iz),rho011conv(iz)
 end do
 write(114,*)
 write(114,*)'Convoluted and normalized momentum distribution in (111) direction:'
 write(114,*)
 write(114,*)'  p111 (mrad):             Probalility:'
 do iz=0,nd2
   write(114,*) mesh111(iz),rho111conv(iz)
 end do
 write(114,*)
! Interpolate results and put them on the same grid (0,01 10^-3m_0c spacing)
! this is for testing only at the moment
 der1=rho111conv(nd1);der2=rho111conv(nd2)
 der3=rho011conv(nd3);der4=rho011conv(nd4)
 der5=rho001conv(nd5);der6=rho001conv(nd6)
 ALLOCATE(der111(nd1:nd2))
 ALLOCATE(der011(nd3:nd4))
 ALLOCATE(der001(nd5:nd6))

 call spline(mesh111,rho111conv,nd2-nd1,der1,der2,der111)
 call spline(mesh011,rho011conv,nd4-nd3,der3,der4,der011)
 call spline(mesh001,rho001conv,nd6-nd5,der5,der6,der001)

 nspline=1001
 ALLOCATE(rho_111(nspline))

 ALLOCATE(rho_011(nspline))

 ALLOCATE(rho_001(nspline))

 ALLOCATE(mesh_spline(nspline))

 do iz=1,nspline
   mesh_spline(iz)=0.1_dp*(iz-1)
 end do

 call splint(nd2-nd1,mesh111,rho111conv,der111,nspline,mesh_spline,rho_111)
 call splint(nd4-nd3,mesh011,rho011conv,der011,nspline,mesh_spline,rho_011)
 call splint(nd6-nd5,mesh001,rho001conv,der001,nspline,mesh_spline,rho_001)

 DEALLOCATE(rho111conv)
 DEALLOCATE(rho011conv)
 DEALLOCATE(rho001conv)
 DEALLOCATE(der111)
 DEALLOCATE(der011)
 DEALLOCATE(der001)
 DEALLOCATE(mesh111)
 DEALLOCATE(mesh011)
 DEALLOCATE(mesh001)

 open(unit=111, file='rho_111', status='replace')
 do iz=1,nspline
   write(111,'(f10.2,es24.15)') mesh_spline(iz),rho_111(iz)
 end do
 close(111)
  open(unit=112, file='rho_011', status='replace')
 do iz=1,nspline
   write(112,'(f10.2,es24.15)') mesh_spline(iz),rho_011(iz)
 end do
 close(112)
 open(unit=113, file='rho_001', status='replace')
 do iz=1,nspline
   write(113,'(f10.2,es24.15)') mesh_spline(iz),rho_001(iz)
 end do
 close(113)

 write(*,*)
 write(*,*) 'Convoluted and normalized spectra on a uniform grid&
& have been written to rho_001, rho_011 and rho_111 files.' 
 write(*,*)

 S_dop=zero;W_dop=zero
 do iz=2,30
   S_dop=S_dop+rho_111(iz)+rho_011(iz)+rho_001(iz)
 end do
 S_dop=(S_dop+(rho_111(1)+rho_011(1)+rho_001(1))/two)*0.1_dp/three
 do iz=107,274
   W_dop=W_dop+rho_111(iz)+rho_011(iz)+rho_001(iz)
 end do
 W_dop=W_dop*0.1_dp/three

 DEALLOCATE(mesh_spline)
 DEALLOCATE(rho_111)
 DEALLOCATE(rho_011)
 DEALLOCATE(rho_001)

!Write results
 write(114,*)'S parameter:',S_dop
 write(114,*)'W parameter:',W_dop
 close(114)

 write(*,*)
 write(*,*) 'S and W parameters and momentum distribution spectra have been written doppler_out file.' 
 write(*,*)
 DEALLOCATE(rho111)
 DEALLOCATE(rho011)
 DEALLOCATE(rho001)

end program posdopspectra

!   subroutines spline and splint copied from the ABINIT source
!   from the m_splines module
!   (/src/28_numeric_noabirule/m_splines.F90)

subroutine spline( t, y, n, ybcbeg, ybcend, ypp )
!
!  Author:
!
!    John Burkardt
!    (XGonze got it from http://www.psc.edu/~burkardt/src/spline/spline.html)
!
!  Parameters:
!
!    Input, integer N, the number of data points; N must be at least 2. 
!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the 
!    spline will actually be linear. 
!
!    Input, double precision T(N), the knot values, that is, the points where data
!    is specified.  The knot values should be distinct, and increasing.
!
!    Input, double precision Y(N), the data values to be interpolated.
!
!    Input, double precision YBCBEG, YBCEND, the values to be used in the boundary
!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!
!    Output, double precision YPP(N), the second derivatives of the cubic spline.


  implicit none

  integer, parameter:: dp=kind(0.d0) 
  integer, intent(in) :: n
  real(dp), intent(in) :: t(n)
  real(dp), intent(in) :: y(n)
  real(dp), intent(in) :: ybcbeg
  real(dp), intent(in) :: ybcend

  real(dp), intent(out) :: ypp(n)

  integer :: ibcbeg
  integer :: ibcend
  integer :: i,k
  real(dp) :: ratio,pinv
  real(dp), allocatable :: tmp(:)

  ALLOCATE(tmp(n))

!
!  XG041127
  ibcbeg=1 ; ibcend=1
  if(ybcbeg>1.0d+30)ibcbeg=0
  if(ybcend>1.0d+30)ibcend=0
!
!  Set the first and last equations.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.d0
    tmp(1) = 0.d0
  else if ( ibcbeg == 1 ) then
    ypp(1) = -0.5d0
    tmp(1) = (3.d0/(t(2)-t(1)))*((y(2)-y(1))/(t(2)-t(1))-ybcbeg)
  end if
  if ( ibcend == 0 ) then
    ypp(n) = 0.d0
    tmp(n) = 0.d0
  else if ( ibcend == 1 ) then
    ypp(n) = 0.5d0
    tmp(n) = (3.d0/(t(n)-t(n-1)))*(ybcend-(y(n)-y(n-1))/(t(n)-t(n-1)))
  end if

!
!  Set the intermediate equations.
!
  do i=2,n-1
   ratio=(t(i)-t(i-1))/(t(i+1)-t(i-1))
   pinv = 1.0d0/(ratio*ypp(i-1) + 2.0d0)
   ypp(i) = (ratio-1.0d0)*pinv
   tmp(i)=(6.0d0*((y(i+1)-y(i))/(t(i+1)-t(i))-(y(i)-y(i-1)) &
&    /(t(i)-t(i-1)))/(t(i+1)-t(i-1))-ratio*tmp(i-1))*pinv
   if (abs(tmp(i))<1.d5*tiny(0.d0)) tmp(i)=0.d0   !MT20050927
  enddo

! Solve the equations
  ypp(n) = (tmp(n)-ypp(n)*tmp(n-1))/(ypp(n)*ypp(n-1)+1.0d0)
  do k=n-1,1,-1
   ypp(k)=ypp(k)*ypp(k+1)+tmp(k)
  enddo

  DEALLOCATE(tmp)

  return
end subroutine spline
!!***
!!***

!----------------------------------------------------------------------

!!****f* m_splines/splint
!! NAME
!!  splint
!!
!! FUNCTION
!!  Compute spline interpolation. There is no hypothesis
!!  about the spacing of the input grid points.
!!
!! INPUTS
!!  nspline: number of grid points of input mesh
!!  xspline(nspline): input mesh
!!  yspline(nspline): function on input mesh
!!  ysplin2(nspline): second derivative of yspline on input mesh
!!  nfit: number of points of output mesh
!!  xfit(nfit): output mesh
!!
!! OUTPUT
!!  yfit(nfit): function on output mesh
!!  [ierr]=A non-zero value is used to signal that some points in xfit exceed xspline(nspline).
!!    The input value is incremented by the number of such points.

subroutine splint(nsplin,xspline,yspline,ysplin2,nfit,xfit,yfit)

 implicit none

 integer, parameter:: dp=kind(0.d0) 
 integer, intent(in) :: nfit, nsplin
 real(dp), intent(in) :: xspline(nsplin)
 real(dp), intent(in) :: yspline(nsplin)
 real(dp), intent(in) :: ysplin2(nsplin)
 real(dp), intent(in) :: xfit(nfit)

 real(dp), intent(out) :: yfit(nfit)

!local
 integer :: left,i,k,right,my_err
 real(dp) :: delarg,invdelarg,aa,bb

!source
 my_err=0

 left = 1
 do i=1, nfit
   yfit(i)=0.d0  ! Initialize for the unlikely event that rmax exceed r(mesh)
   !
   do k=left+1, nsplin

     if(xspline(k) >= xfit(i)) then
       if(xspline(k-1) <= xfit(i)) then
         right = k
         left = k-1
       end if
       delarg= xspline(right) - xspline(left)
       invdelarg= 1.0d0/delarg
       aa= (xspline(right)-xfit(i))*invdelarg
       bb= (xfit(i)-xspline(left))*invdelarg

       yfit(i) = aa*yspline(left) + bb*yspline(right)    &
&               +( (aa*aa*aa-aa)*ysplin2(left) +         &
&                  (bb*bb*bb-bb)*ysplin2(right) ) *delarg*delarg/6.0d0
       exit
     end if
   end do ! k
   !

   if (k==nsplin+1) my_err=my_err+1 ! xfit not found 
 end do ! i

end subroutine splint
