  subroutine make_global_vks_debug(nfft,nspin,global_vks,int_exc)

    integer :: nfft,nspin
    integer :: ii,isp
    real(8) :: tot_rho(nfft,nspin), & 
               tmp_exc,tmp_ehart, & 
               vhart(nfft),vxc(nfft,nspin), & 
               global_vks(nfft,nspin), &
               tmp_vpsp(nfft), int_exc
    character(len=400) :: fname,ss

    print *,''
    print *,''
    print *,'enter make_global_vks()'

    ! total LDA vxc
    tot_rho = ref_rho
    call hartree(n1,n2,n3,nspin,ucvol,qvec,tot_rho,vhart,tmp_ehart)
    select case(global_ixc)
    case(1)
      call calculate_xc_pw92LSDA(n1,n2,n3,nspin,tot_rho,dvol,vxc,tmp_exc);
    case(2)
      call calculate_xc_pbe(n1,n2,n3,nspin,tot_rho,qvec,dvol,vxc,tmp_exc);
    end select
    print *,'[make_global_vks] exc_LDA(rho_total): ',tmp_exc,' hartr: ',tmp_ehart
    int_exc = tmp_exc

    do isp=1,nspin
      global_vks(:,isp) = vhart + vxc(:,isp)
    enddo

    ! get global vpsp 
    do ii=1,nsys
      write(ss,*)ii
      ! local local potential
      fname = './subsys'//trim(adjustl(ss))//'/vpsp.dat'
      open(file=fname,unit=111,action='read',form='unformatted')
      read(111)tmp_vpsp
      close(111)
      do isp=1,nspin
        global_vks(:,isp) = global_vks(:,isp) + tmp_vpsp
      enddo
    enddo

!    open(file='vks_up.dat',action='read',unit=111)
!    do ii=1,nfft
!      read(111,*)global_vks(ii,1)
!    enddo
!    close(111)
!    open(file='vks_down.dat',action='read',unit=111)
!    do ii=1,nfft
!      read(111,*)global_vks(ii,2)
!    enddo
!    close(111)
!
!    return

    ! subsytem exc in LDA or GGA
    do ii=1,nsys
      select case(global_ixc)
      case (1)
        call calculate_xc_pw92LSDA(n1,n2,n3,nspin,sub_rhor(:,:,ii),dvol,vxc,tmp_exc);
      case (2)
        call calculate_xc_pbe(n1,n2,n3,nspin,sub_rhor(:,:,ii),qvec,dvol,vxc,tmp_exc)
      end select
      int_exc = int_exc - tmp_exc
      global_vks = global_vks - vxc 
      write(6,'(a,f16.6,a,i4)')'[make_global_vks] exc_LDA(rho_subsystem): ',tmp_exc,' subsystem: ',ii
    enddo

    ! subsystem vxc (they are dumpped by ABINIT)
    do ii=1,nsys
      write(ss,*)ii
      fname = './subsys'//trim(adjustl(ss))//'/vxc.dat'
      open(file=fname,unit=111,action='read',form='unformatted')
      read(111)tmp_exc
      read(111)vxc
      close(111)
      write(6,'(a,f16.6,a,i4,a)')'[make_global_vks] exc(rho_subsystem): ', & 
      tmp_exc,' subsystem: ',ii,' (Read from abinit_s vxc.dat file)'
      int_exc = int_exc + tmp_exc
      global_vks = global_vks + vxc
    enddo

    print *,'[make_global_vks] int_exc: ',int_exc
    if (nspin==2) then 
      print *,' min/max(global_vks): ',minval(global_vks(:,1)),maxval(global_vks(:,1)),'  spin_up'
      print *,' min/max(global_vks): ',minval(global_vks(:,2)),maxval(global_vks(:,2)),'  spin_dn'
    else 
      print *,' min/max(global_vks): ',minval(global_vks(:,1)),maxval(global_vks(:,1)),' '
    endif
    print *,'exit make_global_vks()'
    print *,''
    print *,''
    print *,''
  end subroutine make_global_vks_debug

