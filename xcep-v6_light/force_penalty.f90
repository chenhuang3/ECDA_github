!
!  compute force due to \partial Penalty / \partial R 
!
subroutine force_penalty(nspin,nfft,n1,n2,n3,natom,qvec, & 
                        vhart,vks,vpsp,pen_coeff,force_pen)

  implicit none 

  integer :: nfft, nspin, natom, n1,n2,n3
  real(8) :: vhart(nfft), vks(nfft,nspin), vpsp(nfft), & 
             pen_coeff, force_pen(3,natom), & 
             qvec(3,(n1/2+1)*n2*n3)

  ! local vars 
  character(len=600) :: fname, fname2
  integer :: isp, stat
  logical :: file_ext
  real(8) :: vd(nfft,nspin)
  real(8) :: lap_vd(nfft,nspin)


  ! dump global ks potential
  do isp=1,nspin
    vd(:,isp) =  vks(:,isp) - vpsp - vhart
    call laplacian(nfft,n1,n2,n3,vd(:,isp),qvec,lap_vd(:,isp))
  enddo 
  open(file='./global_system/vd.dat', action='write',form='unformatted',unit=111)
  write(111)vd
  close(111)
  open(file='./global_system/lap_vd.dat', action='write',form='unformatted',unit=111)
  write(111)lap_vd
  close(111)
  open(file='./global_system/pen_coeff.dat', action='write',form='unformatted',unit=111)
  write(111)pen_coeff
  close(111)
  

  ! write control.dat file 
  open(file='./global_system/control.dat',unit=111,action='write')
  write(111,*) 2       ! type of system, global system KS
  write(111,*) 12,1    ! calc mode ask ABINIT for dExc/dvks; conv_unocc_orb
  write(111,*) -1      ! ask ABINIT for EXX energy density and potential
  close(111)
  

  ! launch global KS DFT (non-self-consistent)
  call chdir("./global_system")
  ! remove done_end.dat
  open(file='done_end.dat',unit=111,iostat=stat,status='old')
  if (stat==0) close(111,status='delete')
  if (stat/=0) close(111)


  !message_from_dfet and message_from_dfet2
  open(file='message_from_dfet2',unit=111,action='write',form='formatted')
  close(111,status='delete')
  !
  inquire(file='message_from_dfet2',exist=file_ext)
  if (file_ext .eqv. .true.) then 
    print *,'failed to remove message_from_dfet2 in run_global_ks.f90'
    call flush(6)
    stop
  endif 
  !
  open(file='message_from_dfet',unit=111,action='write',form='formatted')
  write(111,'(a)')'RUN'
  close(111)
  open(file='message_from_dfet2',unit=111,action='write',form='formatted')
  write(111,'(a)')'RUN'
  close(111)


  ! wait it to finish 
  fname  = 'done_end.dat'
  fname2 = 'done2_end.dat'

  ! check if one shot Kohn-Sham is done.
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

  ! load forces 
  open(file='force_pen.dat',unit=111,action='read',form='unformatted')
  read(111) force_pen
  close(111)

  call chdir("../")

endsubroutine  force_penalty
