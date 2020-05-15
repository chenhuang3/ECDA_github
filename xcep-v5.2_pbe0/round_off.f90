!=============================
! round off input r 
!  e.g. r = 4.3 then return 4
!       r = 1.7 then return 2
!=============================
function round_off(r)
 implicit none
 real(kind=8)  :: r
 integer ::  round_off
 round_off =  floor ((r*10.d0 + 5.d0)/10.d0 )
 return 
end function round_off

