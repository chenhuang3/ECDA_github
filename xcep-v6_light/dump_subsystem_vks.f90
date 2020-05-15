

subroutine dump_subsystem_vks(nsys,nfft,nspin,tmp_vks)

 implicit none 

 integer :: nsys,nfft,nspin
 real(8) :: tmp_vks(nfft,nspin,nsys)

 ! internal vars 
 character(len=500) :: ss,fname,fname2
 integer :: ii,jj

 !
 ! write new KS potential to subsystems 
 !
 do ii=1,nsys
   write(ss,*)ii
   !
   ! dump vks.dat
   ! 
   if (nspin==1) then 
     fname = './subsys'//trim(adjustl(ss))//'/bare_ks_pot.dat'
     open(file=fname,action='write',form='formatted',unit=111)
     do jj=1,nfft
       write(111,*)tmp_vks(jj,1,ii)
     enddo
     close(111)
   else
     fname =  './subsys'//trim(adjustl(ss))//'/bare_ks_pot_up.dat'
     fname2 = './subsys'//trim(adjustl(ss))//'/bare_ks_pot_down.dat'
     open(file=fname, action='write',form='formatted',unit=111)
     open(file=fname2,action='write',form='formatted',unit=112)
     do jj=1,nfft
       write(111,*)tmp_vks(jj,1,ii)
       write(112,*)tmp_vks(jj,2,ii)
     enddo
     close(111)
     close(112)
   endif
 enddo

end subroutine dump_subsystem_vks
