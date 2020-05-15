!
! Tell ABINIT to update the cooridates (xred) when doing 
! geometry relaxation 
!
! Created (1/13/2019) Chen Huang
!
subroutine update_coords_subsystem(logf,myrank,natom,natom_clu,natom_env, & 
                                   atom_class,xcart,xcart_clu,xcart_env)
 use comm_data 

 implicit none 

 ! Declare the interface for POSIX fsync function
 interface
   function fsync (fd) bind(c,name="fsync")
   use iso_c_binding, only: c_int
     integer(c_int), value :: fd
     integer(c_int) :: fsync
   end function fsync
 end interface



 character(len=200) :: fname_clu, fname_env, ss, str2
 logical :: in_cluster, file_exist
 integer :: logf, ret, natom, myrank, ia, kk, stat, & 
            calc_mode,  atom_class(natom), & 
            natom_clu, natom_env
 real(8) :: xcart(3,natom), & 
            xcart_clu(3,natom), & 
            xcart_env(3,natom) 

 
 ! >>>>>>>>>>> function begins <<<<<<<<<!

 write(logf,'(a)')'enter update_coords_subsystem().'
 write(715,'(a)')'enter update_coords_subsystem().'
 call flush(logf)
 call flush(715)



 ! send new geometry to clusters and environments to updated their 
 ! nonlocal projectors (for NCPP) case or PAW terms (for PAW cases)

 write(ss,*) myrank+1

 ! cluster 
 fname_clu = './subsys'//trim(adjustl(ss))//'/new_coords.dat'
 open(file=fname_clu,unit=111,action='write',form='unformatted')
 do ia=1,natom
    if (atom_class(ia)==1) then 
       write(111) xcart(:,ia)
       write(715,'(a,3f16.8,a)') 'cluster atom coord: ',xcart(:,ia)*0.52917721067d0,' ang'
       call flush(715)
    endif
 enddo 
 close(111)


 ! env 
 fname_env = './subsys'//trim(adjustl(ss))//'_env/new_coords.dat'
 open(file=fname_env,unit=111,action='write',form='unformatted')
 do ia=1,natom
   if (atom_class(ia)==2) then ! env
     write(111) xcart(:,ia)
     write(715,'(a,3f16.8,a)') 'env atom coord: ',xcart(:,ia)*0.52917721067d0,' ang'
     call flush(715)
   endif  
 enddo 
 close(111)



 !========================================
 ! inform cluster to update atom coords 
 !========================================
 
 ! ============ cluster =============
 str2='./subsys'//trim(adjustl(ss))//'/'
 call chdir(trim(str2),status=stat)

 ! remove done_end.dat
 open(file='done_end.dat',unit=111,iostat=stat,status='old')
 close(111,status='delete')
 open(file='control.dat',unit=111,action='write',form='formatted')
 write(111,'(i4)') 1  ! tell abinit that we are doing subsystem 
 calc_mode = 11       ! tell abinit update atom coord for cluster
 write(111,'(i4)') calc_mode, 0
 write(111,*) -1
 close(111)

 
 ! message_from_dfet and message_from_dfet2
 ! remove message_from_dfet2 first
 open(file='message_from_dfet2',unit=111,action='write',form='formatted')
 close(111,status='delete')
 inquire(file='message_from_dfet2',exist=file_exist)
 if (file_exist .eqv. .true.) then 
   print *,'failed to remove message_from_dfet2 in run_subsystem_xcpatch.f90, myrank: ',myrank
   call flush(6)
   stop
 endif 
 open(file='message_from_dfet',unit=111,action='write',form='formatted')
 write(111,'(a)') 'RUN'
 call flush(111)
 ret = fsync(fnum(111))
 close(111)
 open(file='message_from_dfet2',unit=111,action='write',form='formatted')
 write(111,'(a)') 'RUN'
 call flush(111)
 ret = fsync(fnum(111))
 close(111)

 call chdir("../",status=stat)
 if (stat/=0) stop 'ERROR, cannot change dir in subsystem'
 




 if (.not. do_envOF) then 
   ! ============ env ==============
   str2='./subsys'//trim(adjustl(ss))//'_env/'
   call chdir(trim(str2),status=stat)
   if (stat/=0) stop 'ERROR, cannot change dir in subsystem'

   ! remove done_end.dat
   open(file='done_end.dat',unit=111,iostat=stat,status='old')
   if (stat==0) close(111,status='delete')
   if (stat/=0) close(111)

   open(file='control.dat',unit=111,action='write',form='formatted')
   write(111,*) 1  ! tell abinit that we are doing subsystem 
   calc_mode = 11  ! for environment to update its atom coords 
   write(111,'(i4)') calc_mode, 0
   write(111,'(i4)') -1
   close(111)


   !message_from_dfet and message_from_dfet2
   !
   ! remove message_from_dfet2 first
   open(file='message_from_dfet2',unit=111,action='write',form='formatted')
   close(111,status='delete')
   inquire(file='message_from_dfet2',exist=file_exist)
   if (file_exist .eqv. .true.) then 
     print *,'failed to remove message_from_dfet2 in run_subsystem_xcpatch.f90, myrank:',myrank
     call flush(6)
     stop
   endif 
   open(file='message_from_dfet',unit=111,action='write',form='formatted')
   write(111,'(a)') 'RUN'
   call flush(111)
   ret = fsync(fnum(111))
   close(111)
   open(file='message_from_dfet2',unit=111,action='write',form='formatted')
   write(111,'(a)') 'RUN'
   call flush(111)
   ret = fsync(fnum(111))
   close(111)

   call chdir("../",status=stat)
   if (stat/=0) stop 'ERROR, cannot change dir in subsystem'
 endif 


 !=======================
 ! wait abinit programs 
 !=======================
 call wait_for_system_end_xcpatch(myrank+1,1) ! wait for cluster 
 if (.not. do_envOF)  & 
   call wait_for_system_end_xcpatch(myrank+1,2) ! wait for env



 write(logf,'(a)')'leave update_coords_subsystem().'
 write(715,'(a)')'leave update_coords_subsystem().'
 write(logf,*)''
 call flush(logf)
 call flush(715)

end subroutine       
