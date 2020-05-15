
!========================================
! write potential u(r) to the disk
!========================================

subroutine dump_uemb_xcpatch(j_atom,u,nfft,nspin)

  use comm_data 

  implicit none

  integer,intent(in)      :: nspin,nfft,j_atom
  real(kind=8),intent(in) :: u(nfft,nspin)
  character(len=100)      :: fil,str1,ss



  ! ============= cluster ==============

  write(ss,*)j_atom 
  ss = 'subsys'//trim(adjustl(ss))
  call chdir(ss)

  ! write u.dat file to subsystem
  fil="uemb.dat"
  open(file=trim(fil),unit=111,action='write',form='unformatted')
  write(111)u
  close(111)
  call chdir('../')

  
  ! ============ env =============
  if (do_envOF .eqv. .false.) then 
    write(ss,*)j_atom 
    ss = 'subsys'//trim(adjustl(ss))//'_env'
    call chdir(ss)
    ! write u.dat file to subsystem
    fil="uemb.dat"
    open(file=trim(fil),unit=111,action='write',form='unformatted')
    write(111)u
    close(111)
    call chdir('../')
  endif 


end subroutine dump_uemb_xcpatch
