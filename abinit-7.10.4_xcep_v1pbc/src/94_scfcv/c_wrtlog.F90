subroutine c_wrtlog(msg)

  implicit none
  character(len=*), intent(in) :: msg;
  INTEGER :: me, ierr
  INTEGER,save :: master = 0

!  me=xcomm_rank(spaceComm)
!  if ( me ==0) then

      ! only master node write to invKS_log file
      ! open file  ...
      open(unit=2222,action='write',file='invKS_log',form='formatted',access='append');
      
      IF (msg(1:13) .eq. 'WELCOMEMSG_BF') THEN
        WRITE(2222,*)'-----------------------------------------------------------------'
        write(2222,*)'      Revised version based on ABINIT v7.10.4                    '
        write(2222,*)'                                                                 '
        write(2222,*)'-    This program solve wavefuntion for a given KS potential    -'
        write(2222,*)'-----------------------------------------------------------------'
      ELSE
        write (2222, '(A)') trim(msg);
      ENDIF
      

      CLOSE (2222)

!  endif

end subroutine c_wrtlog
