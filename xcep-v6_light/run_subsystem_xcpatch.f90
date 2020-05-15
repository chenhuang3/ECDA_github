!
!  run each subsystem 
!
subroutine run_subsystem_xcpatch(j_atom, do_exx, load_new_occ, do_xcep, last_iter)

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

  integer :: s,calc_mode,j_atom, do_exx,  & 
             load_new_occ, do_xcep, & 
             conv_unocc_orb

  logical :: last_iter
  logical :: there
  character(len=255) :: cwd
  character(len=100) :: str1, str2,str3

  integer :: stat, ret
  logical :: file_exist

  conv_unocc_orb = -1

  !
  ! ============= cluster ============
  !
  write(str1,*)j_atom
  str2='./subsys'//trim(adjustl(str1))//'/'
  call chdir(trim(str2),status=stat)
  if (stat/=0) then 
    print *,trim(str2)
    stop 'ERROR, cannot change dir in subsystem'
  endif 

  ! remove done_end.dat
  open(file='done_end.dat',unit=111,iostat=stat,status='old')
  if (stat==0) close(111,status='delete')
  if (stat/=0) close(111)

  open(file='control.dat',unit=111,action='write',form='formatted')
  write(111,'(i4)') 1  ! tell abinit that we are doing subsystem 

  ! final step in embedding potential solver 
  if (load_new_occ>0 .and. do_exx>0) then 
    calc_mode = 200 
  else 
    calc_mode = 4   ! for cluster xc patching
  endif

  ! compute EXX-related quantities for XCEP 
  if (do_xcep>0) then 
    calc_mode = 7
  endif 
  if (last_iter) conv_unocc_orb = 1
  write(111,'(i4)') calc_mode, conv_unocc_orb
  if (do_exx>0) then 
    write(111,*) 1 
  else
    write(111,*) -1
  endif 
  close(111)

  !message_from_dfet and message_from_dfet2
  !
  ! remove message_from_dfet2 first
  open(file='message_from_dfet2',unit=111,action='write',form='formatted')
  close(111,status='delete')
  inquire(file='message_from_dfet2',exist=file_exist)
  if (file_exist .eqv. .true.) then 
    print *,'failed to remove message_from_dfet2 in run_subsystem_xcpatch.f90 '
    call flush(6)
    stop
  endif 
  !
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
  !write(6,'(a,i3,a,2f12.4)')'subsystem for atom ',j_atom,' launched'
  


  !================================
  !  env
  !================================
  if (.not. do_envOF) then 
    str2='./subsys'//trim(adjustl(str1))//'_env/'
    call chdir(trim(str2),status=stat)
    if (stat/=0) stop 'ERROR, cannot change dir in subsystem'

    ! remove done_end.dat
    open(file='done_end.dat',unit=111,iostat=stat,status='old')
    if (stat==0) close(111,status='delete')
    if (stat/=0) close(111)

    open(file='control.dat',unit=111,action='write',form='formatted')
    write(111,*) 1  ! tell abinit that we are doing subsystem 
    calc_mode = 5   ! for environment (xc patch)
    if (last_iter) conv_unocc_orb = 1
    write(111,'(i4)') calc_mode, conv_unocc_orb
    write(111,'(i4)') -1
    close(111)

    !message_from_dfet and message_from_dfet2
    !
    ! remove message_from_dfet2 first
    open(file='message_from_dfet2',unit=111,action='write',form='formatted')
    close(111,status='delete')
    inquire(file='message_from_dfet2',exist=file_exist)
    if (file_exist .eqv. .true.) then 
      print *,'failed to remove message_from_dfet2 in run_subsystem_xcpatch.f90 '
      call flush(6)
      stop
    endif 
    !
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
    call flush(6)
  endif ! do_envOF


end subroutine  run_subsystem_xcpatch






