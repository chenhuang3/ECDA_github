subroutine init_uemb(nspin,nfft,start_u,u)

 implicit none 

 integer :: nspin, nfft, start_u, isp, s
 real(8) :: u(nfft,nspin)


 ! initialize embedding potential 
 if (start_u==0) then   
   u = 0.d0
 else if (start_u==1) then 
   print *,'random initilization of uemb has been disabled at this point! '
   stop
!   do s=1,nfft
!     do isp=1,nspin
!       u(s,isp) = (rand(0)-0.5d0)*0.1d0
!     enddo
!   enddo 
 else if ( start_u==2 ) then 
   ! restart by loading u.restart file 
   print *,''
   print *,'restart by loading u.restart file.'
   print *,''
   open(file='u.restart',unit=111,action='read',form='unformatted')
   read(111)u
   close(111)
 endif


end subroutine init_uemb
