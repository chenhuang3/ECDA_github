subroutine run_global_ks(n1,n2,n3,nfft,nspin,natom,qvec,ucvol,do_global_exx,do_xcep, &
         global_vks,ref_rho,ehart,ke,nlpsp,TS,elocal,ewald,corepsp,total_fermi,dfermi_dvks, & 
         band_energy,exc_eps,exc_exx,dExc_dvks,force_abinit,vhart)

    implicit none 

    ! external vars 
    integer :: nfft,nspin,n1,n2,n3,do_global_exx,do_xcep,natom
    real(8) :: ref_rho(nfft,nspin),etotal,ucvol,& 
               global_vks(nfft,nspin),qvec(3,(n1/2+1)*n2*n3), & 
               dfermi_dvks(nfft,nspin), & 
               ehart,ke,nlpsp,TS,elocal,ewald,corepsp,& 
               exc_eps(nfft), exc_exx, total_fermi(nspin), & 
               dExc_dvks(nfft,nspin), & 
               force_abinit(3,natom)

    ! local vars
    integer :: jj, stat
    logical :: file_ext, there
    real(8) :: vhart(nfft), & 
               band_energy
    character(len=255) :: cwd
    character(len=100) :: str1, str2,str3,fname,fname2


    print *,'enter run_global_ks()'

    ! dump global ks potential
    open(file='./global_system/vks.dat', action='write',form='unformatted',unit=111)
    write(111)global_vks
    close(111)
    

    ! write control.dat file 
    open(file='./global_system/control.dat',unit=111,action='write')
    write(111,*) 2  ! type of system, global system KS
    if (do_global_exx>0) then 
      write(111,*) 7,1     ! calc mode ask ABINIT for dExc/dvks; conv_unocc_orb
      write(111,*) 1       ! ask ABINIT for EXX energy density and potential
    else 
      write(111,*) 1,1     ! calc_mode, we do normal KS DFT; conv_unocc_orb
      write(111,*) -1  ! ask ABINIT for EXX energy density and potential
    endif 
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

    print *,'(run_global_ks) Kohn-Sham DFT on the global system is activated. Waiting ...'
    call flush(6)
    
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

    write(6,'(a,i3,a,i3,a)')'non-self consistent Kohn-Sham on the global system is done.'
    call flush(6)

    ! load global electron density 
    ! ----------------------------
    open(file='rho.dat',unit=111,action='read',form='unformatted')
    read(111) ref_rho(:,1)
    if (nspin==2) read(111) ref_rho(:,2)
    close(111)

    ! For atom i, its y(r) is defined as 
    !
    ! y(r) = \int dr' w_atom(r') * d eps(r')/dv_KS(r)
    ! 
    ! For global system, in ABINIT, we just set w_atom = 1.0
    if (do_global_exx>0) then  
      open(file='yvec.dat',unit=111,action='read',form='unformatted')
      read(111) exc_exx
      read(111) exc_eps
      read(111) dExc_dvks
      close(111)
    endif 

    open(file='dfermi_dvks.dat',unit=111,action='read',form='unformatted')
    read(111) dfermi_dvks 
    close(111)

    open(file='forces.dat',unit=111,action='read',form='unformatted')
    read(111) force_abinit
    close(111)

    ! load results & energies 
    fname = 'result.dat'
    open(file=fname,unit=111,action='read',form='formatted')
    read(111,*) etotal   ! this is actually ke + nlpsp - TS + \int (bare_vks*rhor)
    read(111,*) ke
    read(111,*) nlpsp
    read(111,*) TS
    read(111,*) elocal
    read(111,*) ewald
    read(111,*) corepsp
    read(111,*) band_energy
    read(111,*) total_fermi(1:nspin)
    read(111,*) ehart
    close(111)


    open(file='vhartr.dat',unit=111,action='read',form='unformatted')
    read(111) vhart
    close(111)

    print *,''
    write(6,*)'------------- results from global KS --------------'
    print *,''
    write(6,'(a,f16.6)')'[run_global_ks] ke:      ',ke
    write(6,'(a,f16.6)')'[run_global_ks] nlps:    ',nlpsp
    write(6,'(a,f16.6)')'[run_global_ks] -TS:     ',TS
    write(6,'(a,f16.6)')'[run_global_ks] elocal:  ',elocal
    write(6,'(a,f16.6)')'[run_global_ks] ehartree:',ehart
    write(6,'(a,f16.6)')'[run_global_ks] ewald:   ',ewald
    write(6,'(a,f16.6)')'[run_global_ks] corepsp: ',corepsp
    print *,''
    print *,'exit run_global_ks()'
    print *,''

!    call chdir(trim(cwd))
    call chdir("../")

end subroutine run_global_ks
