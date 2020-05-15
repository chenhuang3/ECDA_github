

subroutine check_bfgs_flag(task, new_x, nfft, nspin, u, iter_umin, etot)

  implicit none 

  integer           :: iter_umin, nfft, nspin  ! BFGS iterations 
  logical           :: new_x
  character(len=60) :: task
  real(kind=8)      :: etot, u(nfft,nspin)

  print *,''
  print *,'    ------- check_bfgs_flag.f90 ------'
  print *,''

  new_x = .FALSE.

  if (task(1:5).eq.'NEW_X') then 
    
    if (nspin==1) then 
      write(6,'(a,i3,a,2es12.4,a,es12.4,a,f14.6)')& 
       'BFGS: TASK=NEW_X, iter-> ',iter_umin,' min/max(u): ',minval(u),maxval(u), & 
       ' amp(u): ',maxval(u)-minval(u),' W: ',etot
    else if (nspin==2) then 
      write(6,'(a,i3,a,2es12.4,a,es12.4)')& 
       'BFGS: TASK=NEW_X, iter-> ',iter_umin,' min/max(u_alpha): ',minval(u(:,1)),maxval(u(:,1)), & 
       ' amp(u_alpha): ',maxval(u(:,1))-minval(u(:,1))
      write(6,'(a,i3,a,2es12.4,a,es12.4,a,f14.6)')& 
       'BFGS: TASK=NEW_X, iter-> ',iter_umin,' min/max(u_beta) : ',minval(u(:,2)),maxval(u(:,2)), & 
       ' amp(u_beta ): ',maxval(u(:,2))-minval(u(:,2)),' W: ',etot
    endif
    write(6,*)''
    new_x = .TRUE.

  elseif (task(1:2) .EQ. 'FG') then

    if (nspin==1) then 
      write(6,'(a,2es12.4,a,f14.6)')'BFGS: TASK=FG, min/max(u): ',minval(u),maxval(u),' W: ',etot
    else if(nspin==2) then 
      write(6,'(a,2es12.4,a,f14.6)')'BFGS: TASK=FG, min/max(u_alpha): ',minval(u(:,1)),maxval(u(:,1))
      write(6,'(a,2es12.4,a,f14.6)')'BFGS: TASK=FG, min/max(u_beta ): ',minval(u(:,2)),maxval(u(:,2)),' W: ',etot
    endif 
    write(715,'(a,es16.8)')'TASK=FG  cgW: ',etot
    call flush(715)
    print *,''

  elseif( task(1:5) .eq. 'START') then

    if (nspin==1) then 
      write(6,'(a,i3,a,2es12.4,a,es12.4,a,f14.6)')& 
       'BFGS: TASK=START iter-> ',0,' min/max(u): ',minval(u),maxval(u), & 
       ' amp(u): ',maxval(u)-minval(u),' W: ',etot
    else if (nspin==2) then 
      write(6,'(a,i3,a,2es12.4,a,es12.4)')& 
       'BFGS: TASK=START iter-> ',0,' min/max(u_alpha): ',minval(u(:,1)),maxval(u(:,1)), & 
       ' amp(u_alpha): ',maxval(u(:,1))-minval(u(:,1))
      write(6,'(a,i3,a,2es12.4,a,es12.4,a,f14.6)')& 
       'BFGS: TASK=START iter-> ',0,' min/max(u_beta) : ',minval(u(:,2)),maxval(u(:,2)), & 
       ' amp(u_beta ): ',maxval(u(:,2))-minval(u(:,2)),' W: ',etot
    endif
    write(6,*)''
  else  

    write(6,'(a,a)')'unknown task. task from bfgs => ',task(1:10)
    new_x = .TRUE.

  endif

end subroutine check_bfgs_flag
