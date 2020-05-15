!
! For a given KS potential, solve the KS equation, and compute the 
! perturbed density for the perturbing potential vperb
! On exit, rho_perb is the perturbed density 
!
subroutine calc_subsystem_dfpt(natom,natom_clu,natom_env,j_atom,nfft,nspin,ucvol, & 
                               cluster_weight,env_weight,qvec,sub_rhor,cell_nfft,& 
                               dfermi_dvks,vks_cluster,vks_env,vperb,rho_perb)

  use comm_data

  implicit none 

  integer :: j_atom, nfft, nspin, & 
             natom, & 
             cell_nfft(3), & 
             natom_clu, & 
             natom_env 

  real(8) :: ucvol, & 
             sub_rhor(nfft,nspin,2), & 
             qvec(3,(cell_nfft(1)/2+1)*cell_nfft(2)*cell_nfft(3)), & 
             dfermi_dvks(nfft,nspin,2), & 
             vks_cluster(nfft,nspin), & 
             vks_env(nfft,nspin), & 
             vperb(nfft,nspin) , &   ! perturbing potnetial 
             rho_clu(nfft,nspin), &  ! perturbed density of cluster 
             rho_env(nfft,nspin), & 
             rho_perb(nfft,nspin), & 
             cluster_weight(nfft), & 
             env_weight(nfft)

  ! local vars 
  logical :: file_ext, trans_A
  character(len=255) :: cwd,stmp
  character(len=500) :: str1, str2,str3,fname,fname2
  real(8) :: vec(nfft), &
             OF_chempot, & 
             dmu
  integer  :: stat,isp


  !================================
  !  cluster 
  !================================
  write(str1,*)j_atom
  str2='./subsys'//trim(adjustl(str1))//'/'
  call chdir(trim(str2))

  ! remove done_end.dat
  open(file='done_end.dat',unit=111,iostat=stat,status='old')
  if (stat==0) close(111,status='delete')

  open(file='control.dat',unit=111,action='write',form='formatted')
  write(111,*) 1       ! tell abinit that we are doing subsystem 
  write(111,*) 500,1   ! calc_mode, conv_unocc_orb
  write(111,*) -1      ! do_EXX
  close(111)

  open(file='vks_q.dat',unit=111,action='write',form='unformatted')
  write(111) vks_cluster
  write(111) vperb
  close(111)  

  !message_from_dfet and message_from_dfet2
  open(file='message_from_dfet2',unit=111,action='write',form='formatted')
  close(111,status='delete')
  inquire(file='message_from_dfet2',exist=file_ext)
  if (file_ext .eqv. .true.) then 
    print *,'failed to remove message_from_dfet2 in calc_subsystem_dfpt.f90'
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
!  write(6,'(a,i3,a,2f12.4)')'(dfpt) subsystem for atom ',j_atom,' launched'
!  call flush(6)




  !===============================================
  !  environment 
  !===============================================
  if (do_envOF) then 
    !
    ! env is treated by OFDFT =============
    !
    if (natom_env>0) then 
      write(logOF,'(a)')'In calc_subsystem_dfpt() -> call OF_dfpt(()'
      do isp=1,nspin 
        call OF_dfpt(nfft,cell_nfft,ucvol,qvec, & 
          sub_rhor(:,isp,2),vperb(:,isp),rho_env(:,isp),dmu)
      enddo
    else 
      rho_env = 0.d0 
    endif 
  else 
    !
    ! env is treated by KS-DFT ============
    !
    if (natom_env>0) then 
      str2='./subsys'//trim(adjustl(str1))//'_env/'
      call chdir(trim(str2))

      ! remove done_end.dat
      open(file='done_end.dat',unit=111,iostat=stat,status='old')
      if (stat==0) close(111,status='delete')

      open(file='control.dat',unit=111,action='write',form='formatted')
      write(111,*) 1       ! tell abinit that we are doing subsystem 
      write(111,*) 500,1   ! calc_mode, conv_unocc_orb
      write(111,*) -1
      close(111)

      open(file='vks_q.dat',unit=111,action='write',form='unformatted')
      write(111) vks_env
      write(111) vperb
      close(111)

      !message_from_dfet and message_from_dfet2
      open(file='message_from_dfet2',unit=111,action='write',form='formatted')
      close(111,status='delete')
      inquire(file='message_from_dfet2',exist=file_ext)
      if (file_ext .eqv. .true.) then 
        print *,'failed to remove message_from_dfet2 in calc_subsystem_dfpt.f90 '
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
  endif 




  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !  wait for job to finish 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  write(stmp,*) j_atom

  ! ------------ cluster --------
  fname  = 'subsys'//trim(adjustl(stmp))//'/done_end.dat'
  fname2 = 'subsys'//trim(adjustl(stmp))//'/done2_end.dat'
  call flush(6)

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
  print *,'cluster ',j_atom,' is done'
  call flush(6)

  ! ========= load in cluster density =============
  str2='./subsys'//trim(adjustl(str1))//'/'
  call chdir(trim(str2))
  open(file='rho_perb.dat',action='read',unit=111,form='unformatted')
  read(111) rho_clu
  call chdir("../")




  ! ------------ env --------------
  if (.not. do_envOF) then 
    if (natom_env > 0) then 
      fname  = 'subsys'//trim(adjustl(stmp))//'_env/done_end.dat'
      fname2 = 'subsys'//trim(adjustl(stmp))//'_env/done2_end.dat' 
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
      print *,'env ',j_atom,' is done'
      call flush(6)
      ! ========= load in env  density =============
      str2='./subsys'//trim(adjustl(str1))//'_env/'
      call chdir(trim(str2))
      open(file='rho_perb.dat',action='read',unit=111,form='unformatted')
      read(111)rho_env
      call chdir("../")
    else 
      rho_env = 0.d0 
    endif 
  endif 


  ! for local fermi scheme, by considering the change of cluster and env's 
  ! Fermi levels, we compute y^T Aclu^{-1} (\chi_clu A_clu^{-1} + \chi_env A_env^{-1})
  rho_perb = rho_clu + rho_env

end subroutine
