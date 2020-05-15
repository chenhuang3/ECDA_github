!
!
!
subroutine wait_for_system_end_xcpatch(j_atom,sys)

 implicit none
 logical           :: file_ext
 integer           :: sys,stat,dtmp,j_atom,scf_code=1
 character(len=500):: fname,msg,end_msg,stmp,stmp2,fname2

 !  write(6,'(a,i4,a)') 'waiting for j_atom ',j_atom,' to end ...'
   write(stmp,*) j_atom

   ! ------------ cluster --------
   if (sys==1) fname  = 'subsys'//trim(adjustl(stmp))//'/done_end.dat'
   if (sys==1) fname2 = 'subsys'//trim(adjustl(stmp))//'/done2_end.dat'

   ! ------------ env --------------
   if (sys==2) fname  = 'subsys'//trim(adjustl(stmp))//'_env/done_end.dat'
   if (sys==2) fname2 = 'subsys'//trim(adjustl(stmp))//'_env/done2_end.dat'


   ! check if scf is done.
   file_ext = .false.
   do while (.not. file_ext)
     inquire(file=trim(fname),exist=file_ext)
     call sleep(1)
   enddo  

   file_ext = .false.
   do while (.not. file_ext)
     inquire(file=trim(fname2),exist=file_ext)
     call sleep(1)
   enddo

   !!write(6,'(a,i3,a,i3,a)')'subsystem for atom ',j_atom,' ended.'

end subroutine  
