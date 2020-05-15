
subroutine  set_clu_env_q(nspin,j_atom,q_clu,q_env,natom_clu,natom_env)

 use comm_data 


 implicit none 

 integer :: j_atom, nspin, & 
            natom_clu, & 
            natom_env 
 real(8) :: q_clu(nspin), q_env(nspin)

 ! local vars 
 logical :: file_ext 
 integer :: sys,calc_mode, ib, ispin, stat
 character(len=200) fname, ss, str1, str2, stmp, fname2


 write(715,'(a)') NEW_LINE('a')//'enter set_clu_env_q().'
 write(715,'(a,i4)') 'natom_clu: ',natom_clu
 write(715,'(a,i4)') 'natom_env: ',natom_env
 call flush(715)


 ! launch abinit 

 ! ============= cluster ============
 write(str1,*)j_atom
 str2='./subsys'//trim(adjustl(str1))//'/'
 call chdir(trim(str2))

 open(file='new_q.dat',unit=111,form='formatted',action='write')
 if (nspin==1) write(111,'(f16.8)')  q_clu(1:nspin)
 if (nspin==2) write(111,'(2f16.8)') q_clu(1:nspin)
 close(111)

 ! remove done_end.dat and rho.dat 
 open(file='done_end.dat',unit=111,iostat=stat,status='old')
 if (stat==0) close(111,status='delete')
 open(file='rho.dat',unit=111,iostat=stat,status='old')
 if (stat==0) close(111,status='delete')

 open(file='control.dat',unit=111,action='write',form='formatted')
 write(111,*) 1  ! tell abinit that we are doing subsystem 
 calc_mode = 600   ! for cluster xc patching
 write(111,*) calc_mode, -1 
 write(111,*) -1   ! do_Exx
 close(111)

 !message_from_dfet and message_from_dfet2
 open(file='message_from_dfet2',unit=111,action='write',form='formatted')
 close(111,status='delete')
 inquire(file='message_from_dfet2',exist=file_ext)
 if (file_ext .eqv. .true.) then 
    print *,'[set_clu_env] failed to remove message_from_dfet2'
    call flush(6)
    stop
 endif 
 !
 open(file='message_from_dfet',unit=111,action='write',form='formatted')
 write(111,'(a)') 'RUN'
 close(111)
 open(file='message_from_dfet2',unit=111,action='write',form='formatted')
 write(111,'(a)') 'RUN'
 close(111)
 call chdir("../")
 

 !
 ! ============ env ==============
 !
 if ( (do_envOF .eqv. .false.) .and. (natom_env > 0)) then 

   str2='./subsys'//trim(adjustl(str1))//'_env/'
   call chdir(trim(str2))

   open(file='new_q.dat',unit=111,form='formatted',action='write')
   write(111,'(2f16.8)') q_env(1:nspin)
   close(111)

   ! remove done_end.dat and rho.dat 
   open(file='done_end.dat',unit=111,iostat=stat,status='old')
   if (stat==0) close(111,status='delete')
   open(file='rho.dat',unit=111,iostat=stat,status='old')
   if (stat==0) close(111,status='delete')

   open(file='control.dat',unit=111,action='write',form='formatted')
   write(111,*) 1  ! tell abinit that we are doing subsystem 
   calc_mode = 600   ! for environment (xc patch)
   write(111,*) calc_mode, -1
   write(111,*) -1   ! do_Exx
   close(111)

   !message_from_dfet and message_from_dfet2
   open(file='message_from_dfet2',unit=111,action='write',form='formatted')
   close(111,status='delete')
   inquire(file='message_from_dfet2',exist=file_ext)
   if (file_ext .eqv. .true.) then 
     print *,'[set_clu_env] failed to remove message_from_dfet2'
     call flush(6)
     stop
   endif 
   !
   open(file='message_from_dfet',unit=111,action='write',form='formatted')
   write(111,'(a)') 'RUN'
   close(111)
   open(file='message_from_dfet2',unit=111,action='write',form='formatted')
   write(111,'(a)') 'RUN'
   close(111)
   call chdir("../")
 endif ! natom_env > 0



 !==================================
 ! wait cluster and env to finish
 !==================================
 do sys=1,2 

   ! no atom in environment, just skip 
   if (sys==2 .and. natom_env==0) cycle 
   if (sys==2 .and. do_envOF) exit

   write(stmp,*) j_atom 
   ! ------------ cluster --------
   if (sys==1) fname  = 'subsys'//trim(adjustl(stmp))//'/done_end.dat'
   if (sys==1) fname2 = 'subsys'//trim(adjustl(stmp))//'/done2_end.dat'

   ! ------------ env --------------
   if (sys==2) fname  = 'subsys'//trim(adjustl(stmp))//'_env/done_end.dat'
   if (sys==2) fname2 = 'subsys'//trim(adjustl(stmp))//'_env/done2_end.dat'

   ! check if abinit is done.
   file_ext = .false.
   do while (.not. file_ext)
     inquire(file=trim(fname),exist=file_ext)
   enddo  
   file_ext = .false.
   do while (.not. file_ext)
     inquire(file=trim(fname2),exist=file_ext)
   enddo

   ! remove done_end.dat
   open(file='done_end.dat',unit=111,iostat=stat,status='old')
   if (stat==0) close(111,status='delete')
 enddo

  
 write(715,'(a)')'done set_clu_env_q()'
 call flush(715)

end subroutine 
