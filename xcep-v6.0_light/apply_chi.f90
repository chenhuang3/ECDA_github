!
! For a given KS potential, solve the KS equation, and compute the 
! perturbed density for the perturbing potential vperb
! On exit, rho_perb is the perturbed density 
!
subroutine apply_chi(info,system_type,j_atom,nfft,nspin,vks,vperb,rho_perb)

  implicit none 

  integer :: system_type, &   ! 0: total system 
                              ! 1: cluster
                              ! 2: env 
             j_atom, nfft, nspin, &
             type_of_system

  real(8) :: vks(nfft,nspin), & 
             vperb(nfft,nspin) , &   ! perturbing potnetial 
             rho_perb(nfft,nspin)


  integer            :: stat, info, ierror
  logical            :: file_ext
  character(len=255) :: cwd,stmp
  character(len=500) :: str1, str2,str3,fname,fname2


  ! >>>>>>>>>>>>>> function begins <<<<<<<<<<<<<<<<<<!

  if (system_type==0) then 
     ! total system 
     str2='./global_system/'
     call chdir(trim(str2))

     open(file='done_end.dat',unit=111,iostat=stat,status='old')
     if (stat==0) close(111,status='delete',iostat=ierror)
     if (stat/=0) close(111)

     type_of_system = 2
     print *,'launching ./global_system/ ...'

  elseif (system_type == 1) then 
     ! ========= cluster ==========
     write(str1,*)j_atom
     str2='./subsys'//trim(adjustl(str1))//'/'
     call chdir(trim(str2))

     ! remove done_end.dat
     open(file='done_end.dat',unit=111,iostat=stat,status='old')
     if (stat==0) close(111,status='delete')

     type_of_system = 1

  elseif (system_type == 2) then 
     ! ============ env ==============
     write(str1,*)j_atom
     str2='./subsys'//trim(adjustl(str1))//'_env/'
     call chdir(trim(str2))

     ! remove done_end.dat
     open(file='done_end.dat',unit=111,iostat=stat,status='old')
     if (stat==0) close(111,status='delete')

     type_of_system = 1
  endif


  !
  ! Write control file 
  !
  open(file='control.dat',unit=111,action='write',form='formatted')
  write(111,*) type_of_system  ! tell abinit that we are doing subsystem 
  if (info==501) then 
    write(111,*) 501,1  ! calc_mode (speical case, consider Fermi level equilibrated between both clustr and env)
  else 
    write(111,*) 500,1  ! calc_mode (do DFPT)
  endif 
  write(111,*) -1
  close(111)

  open(file='vks_q.dat',unit=111,action='write',form='unformatted')
  write(111) vks
  write(111) vperb
  close(111)
  



  !message_from_dfet and message_from_dfet2
  open(file='message_from_dfet2',unit=111,action='write',form='formatted')
  close(111,status='delete')
  inquire(file='message_from_dfet2',exist=file_ext)
  if (file_ext .eqv. .true.) then 
    print *,'failed to remove message_from_dfet2 in apply_chi.f90'
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



  !
  ! ======= wait jobs to finish =============
  !
  fname  = 'done_end.dat'
  fname2 = 'done2_end.dat'
 ! write(6,*)'files to monitor are: ',trim(fname)
 ! write(6,*)'files to monitor are: ',trim(fname2)
  call flush(6)
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

  ! ========= load in cluster density =============
  open(file='rho_perb.dat',action='read',unit=111,form='unformatted')
  read(111) rho_perb

  ! return to the working dir 
  call chdir("../")


end subroutine
