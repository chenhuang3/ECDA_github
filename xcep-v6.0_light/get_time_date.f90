

subroutine get_time_date(date_string)

   character(len=100) :: date_string
   integer*4 :: today(3), now(3)
   
!!   call idate(today(1),today(2),today(3))   ! today(1)=day, (2)=month, (3)=year
   call idate(today)   ! today(1)=day, (2)=month, (3)=year
   call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
   write ( date_string, 1000 )  today(2), today(1), today(3), now
   1000 format ( 'date ', i2.2, '/', i2.2, '/', i4.4, ' time ', &
            i2.2, ':', i2.2, ':', i2.2 )


end subroutine 
