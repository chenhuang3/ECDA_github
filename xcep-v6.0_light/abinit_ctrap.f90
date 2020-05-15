!
!
! Do corrected trapezoidal integral on uniform grid of spacing hh.
! from ABINIT code  ctrap.F90 file
!
!
!INPUTS
!
!  imax=highest index of grid=grid point number of upper limit
!  ff(imax)=integrand values
!  hh=spacing between x points
!
!OUTPUT
!
!  ans=resulting integral by corrected trapezoid
!

subroutine ctrap(imax,ff,hh,ans)

 implicit none

 integer,parameter :: dp=8
 integer,parameter :: eight=8.d0, nine=9.d0, two=2.d0, four=4.d0, three=3.d0

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: imax
 real(dp),intent(in) :: hh
 real(dp),intent(out) :: ans
!arrays
 real(dp),intent(in) :: ff(imax)

!Local variables-------------------------------
!scalars
 integer :: ir,ir2
 real(dp) :: endpt,sum

! *************************************************************************

 if (imax>=10)then

!  endpt=end point correction terms (low and high ends)
   endpt  = (23.75d0*(ff(1)+ff(imax  )) &
&   + 95.10d0*(ff(2)+ff(imax-1)) &
&   + 55.20d0*(ff(3)+ff(imax-2)) &
&   + 79.30d0*(ff(4)+ff(imax-3)) &
&   + 70.65d0*(ff(5)+ff(imax-4)))/ 72.d0
   ir2 = imax - 5
   sum=0.00d0
   if (ir2 > 5) then
     do ir=6,ir2
       sum = sum + ff(ir)
     end do
   end if
   ans = (sum + endpt ) * hh

 else if (imax>=8)then
   endpt  = (17.0d0*(ff(1)+ff(imax  )) &
&   + 59.0d0*(ff(2)+ff(imax-1)) &
&   + 43.0d0*(ff(3)+ff(imax-2)) &
&   + 49.0d0*(ff(4)+ff(imax-3)) )/ 48.d0
   sum=0.0d0
   if(imax==9)sum=ff(5)
   ans = (sum + endpt ) * hh

 else if (imax==7)then
   ans = (17.0d0*(ff(1)+ff(imax  )) &
&   + 59.0d0*(ff(2)+ff(imax-1)) &
&   + 43.0d0*(ff(3)+ff(imax-2)) &
&   + 50.0d0* ff(4)                )/ 48.d0  *hh

 else if (imax==6)then
   ans = (17.0d0*(ff(1)+ff(imax  )) &
&   + 59.0d0*(ff(2)+ff(imax-1)) &
&   + 44.0d0*(ff(3)+ff(imax-2)) )/ 48.d0  *hh

 else if (imax==5)then
   ans = (     (ff(1)+ff(5)) &
&   + four*(ff(2)+ff(4)) &
&   + two * ff(3)         )/ three  *hh

 else if (imax==4)then
   ans = (three*(ff(1)+ff(4)) &
&   + nine *(ff(2)+ff(3))  )/ eight  *hh

 else if (imax==3)then
   ans = (     (ff(1)+ff(3)) &
&   + four* ff(2)         )/ three *hh

 else if (imax==2)then
   ans = (ff(1)+ff(2))/ two  *hh

 else if (imax==1)then
   ans = ff(1)*hh

 end if

end subroutine ctrap