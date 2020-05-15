!
! inform ABINIT to compute the new density 
! with the new occupation numbers 
! 
subroutine new_occ_abinit_density(j_atom,mxband,nband,nfft,nspin,occ,rho,eigenvalues,fermi)

 implicit none 

 integer :: nspin, mxband, nband(2), j_atom, nfft
 real(8) :: occ(mxband,nspin,2), & 
            fermi(nspin), & 
            eigenvalues(mxband,nspin,2), & 
            rho(nfft,nspin,2)

 ! local vars 

 logical :: file_ext 
 integer :: sys,calc_mode, ib, ispin, stat
 character(len=200) fname, ss, str1, str2, stmp, fname2


 ! write new occ to cluster and env 
 ! --------------------------------

 ! loop over cluster and env 
 do sys=1,2
   write(ss,*)j_atom
   ! write eigenvalues of the subsystem 
   if (sys==1) fname = 'subsys'//trim(adjustl(ss))//'/new_occ.xcpp'
   if (sys==2) fname = 'subsys'//trim(adjustl(ss))//'_env/new_occ.xcpp'
   open(file=fname,unit=111,action='write',form='formatted')
   ! 
   ! fermi levels
   !
   if (nspin==2)  write(111,'(a,2es22.12)')'fermi: ',fermi(1:nspin)
   if (nspin==1)  write(111,'(a, es22.12)')'fermi: ',fermi(1:nspin)
   !
   ! write the occupation numbers
   !
   do ispin=1,nspin 
     do ib=1,nband(sys)
       write(111,'(es22.12,es22.12,a)') occ(ib,ispin,sys),eigenvalues(ib,ispin,sys),' occ, eigen '
     enddo 
   enddo
   close(111)
 enddo



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! launch abinit 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! ============= cluster ============

 write(str1,*)j_atom
 str2='./subsys'//trim(adjustl(str1))//'/'
 call chdir(trim(str2))

 ! remove done_end.dat and rho.dat 
 open(file='done_end.dat',unit=111,iostat=stat,status='old')
 if (stat==0) close(111,status='delete')
 open(file='rho.dat',unit=111,iostat=stat,status='old')
 if (stat==0) close(111,status='delete')

 open(file='control.dat',unit=111,action='write',form='formatted')
 write(111,*) 1  ! tell abinit that we are doing subsystem 
 calc_mode = 100            ! for cluster xc patching
 write(111,*) calc_mode, -1 ! calc_mode, conv_unocc_orb 
 write(111,*)'-1  do_exx'
 close(111)

 !message_from_dfet and message_from_dfet2
 open(file='message_from_dfet2',unit=111,action='write',form='formatted')
 close(111,status='delete')
 inquire(file='message_from_dfet2',exist=file_ext)
 if (file_ext .eqv. .true.) then 
    print *,'failed to remove message_from_dfet2 in new_occ_abinit_density.f90'
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
 
 ! ============ env ==============

 str2='./subsys'//trim(adjustl(str1))//'_env/'
 call chdir(trim(str2))

 ! remove done_end.dat and rho.dat 
 open(file='done_end.dat',unit=111,iostat=stat,status='old')
 if (stat==0) close(111,status='delete')
 open(file='rho.dat',unit=111,iostat=stat,status='old')
 if (stat==0) close(111,status='delete')

 open(file='control.dat',unit=111,action='write',form='formatted')
 write(111,*) 1  ! tell abinit that we are doing subsystem 
 calc_mode = 100            ! for environment (xc patch)
 write(111,*) calc_mode, -1 ! calc_mode, conv_unocc_orb
 write(111,*) -1
 close(111)

 !message_from_dfet and message_from_dfet2
 open(file='message_from_dfet2',unit=111,action='write',form='formatted')
 close(111,status='delete')
 inquire(file='message_from_dfet2',exist=file_ext)
 if (file_ext .eqv. .true.) then 
    print *,'failed to remove message_from_dfet2 in new_occ_abinit_density.f90'
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
 write(6,'(a,i3,a)') '[new_occ_abinit_density] wait for atom:',j_atom,' to finish'
 call flush(6)


 ! wait abinit to finish
 ! ---------------------

 do sys=1,2 
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


 ! collect density 
 ! ---------------

 do sys=1,2 
   write(ss,*) j_atom 
   if (sys==1) fname = 'subsys'//trim(adjustl(ss))//'/rho.dat'
   if (sys==2) fname = 'subsys'//trim(adjustl(ss))//'_env/rho.dat'
   open(file=fname,unit=111,action='read',form='unformatted')
   read(111) rho(:,1,sys)  
   if ( nspin==2) then 
     read(111) rho(:,2,sys)
   endif
   close(111)
 enddo

 print *,'finish new_occ_abinit_density.f90'
 call flush(6)

end subroutine 
