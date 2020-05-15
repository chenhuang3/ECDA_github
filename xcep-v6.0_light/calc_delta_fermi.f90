!
! If global_fermi>0, we consider the charge transfer between cluster and env 
! in that case, cluster and env share the same Fermi level and we need to 
! modify the standard density functional perturbation theory.
!
! To compute the last term (due to the change of the Fermi level) in the Eq. 75 in 
! "Phonons and related crystal properties from density-functional 
! perturbation theory" by Stefano Baroni, Stefano de Gironcoli, and Andrea Dal Corso 
! Reivews of Modern Physics, 73, 515 (2001), this subroutine simply sums the change 
! of the Fermi levels from cluster and env together.
! 
! This subroutine computes the change of the common Fermi level with respect 
! to the change of the embedding potential. It simpily lets cluster and env
! computes their own fermi level change and then sum the two together. 
!
!

subroutine calc_delta_fermi(j_atom,nspin,nfft,vperb)

  implicit none 

  logical :: file_ext
  integer :: j_atom, & 
             nspin, & 
             nfft, isp, & 
             stat, & 
             system_type
  real(8) :: delta_Ef(nspin,2), & ! the last index for cluster and env
             vperb(nfft,nspin) 
  character(len=200) :: str1, str2, fname, fname2


  !
  ! loop over cluster and environment to get their fermi level change
  !
  do system_type =1,2
  
    ! cluster 
    if (system_type == 1) then 
      write(str1,*)j_atom
      str2='./subsys'//trim(adjustl(str1))//'/'
      call chdir(trim(str2))  
      ! remove done_end.dat
      open(file='done_end.dat',unit=111,iostat=stat,status='old')
      if (stat==0) close(111,status='delete')
    endif 

    ! env 
    if (system_type == 2) then 
      write(str1,*)j_atom
      str2='./subsys'//trim(adjustl(str1))//'_env/'
      call chdir(trim(str2)) 
      ! remove done_end.dat
      open(file='done_end.dat',unit=111,iostat=stat,status='old')
      if (stat==0) close(111,status='delete')
    endif


    !~~~~~~~~~~~~~~~~~~~~~~~~
    ! write control file 
    !~~~~~~~~~~~~~~~~~~~~~~~~
    open(file='control.dat',unit=111,action='write',form='formatted')
    write(111,*) 1      ! tell abinit that we are doing subsystem 
    write(111,*) 502,-1 ! calc_mode (just compute Delta_Ef)
    write(111,*) -1
    close(111)

    open(file='vperb.dat',unit=111,action='write',form='unformatted')
    write(111) vperb
    close(111)


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! run ABINIT jobs 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    open(file='message_from_dfet2',unit=111,action='write',form='formatted')
    close(111,status='delete')
    inquire(file='message_from_dfet2',exist=file_ext)
    if (file_ext .eqv. .true.) then 
      print *,'failed to remove message_from_dfet2 in apply_chi.f90'
      call flush(6)
      stop
    endif 
    open(file='message_from_dfet',unit=111,action='write',form='formatted')
    write(111,'(a)') 'RUN'
    close(111)
    open(file='message_from_dfet2',unit=111,action='write',form='formatted')
    write(111,'(a)') 'RUN'
    close(111)


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! wait jobs to finish
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fname  = 'done_end.dat'
    fname2 = 'done2_end.dat'
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

    ! read in delta_Ef 
    open(file='my_delta_Ef.dat',action='read',unit=111,form='formatted')
    do isp=1,nspin 
      read(111,*) delta_Ef(isp,system_type)
    enddo 
    close(111)

    ! return to the working dir 
    call chdir("../")
  enddo  ! loop over cluster and env 




  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! write Delta fermi to cluster and env 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do system_type =1,2

    if (system_type == 1) then 
      ! cluster 
      write(str1,*)j_atom
      str2='./subsys'//trim(adjustl(str1))//'/'
      call chdir(trim(str2))  
    endif 

    if (system_type == 2) then 
      ! env 
      write(str1,*)j_atom
      str2='./subsys'//trim(adjustl(str1))//'_env/'
      call chdir(trim(str2)) 
    endif

    open(file='delta_fermi.dat',unit=111,action='write',form='formatted')
    do isp=1,nspin 
      write(111,*) sum(delta_Ef(isp,:))
    enddo 
    close(111)

    call chdir("../")
  enddo


end subroutine calc_delta_fermi
