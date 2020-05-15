!
!  collect data from subsystems,
!  energies, density in real space
!  and in fourier space
!   
!  Chen Huang, Jan/2011
!
subroutine get_subsystem_data_xcpatch(j_atom,sys,nfft,dvol, & 
             nspin,etotal,nlpsp,TS,ke,rhor,fermi, & 
             mxband,eigenvalues,occ,nband, & 
             do_exx,exc_exx,cluster_exc_eps,yvec,load_dEdVks,drho_dN,deps_dN)

 implicit none

 ! external vars ..........
 integer      :: sys, &   ! 1: cluster, 2: env
                 nfft,nspin,j_atom,do_exx,load_dEdVks

 integer      :: mxband, nband 
 real(kind=8) :: eigenvalues(mxband,nspin), &   ! eigenvalues from the ABINIT program 
                 occ(mxband,nspin)              ! occupation number of KS oribtals 

 real(kind=8) :: etotal, dvol, &      ! etotal: total energies for each subsystem 
                 TS, nlpsp,ke, exc, ehart, elocal, & !
                 cluster_exc_eps(nfft,nspin), & 
                 exc_exx, &
                 drho_dN(nfft,nspin), & 
                 yvec(nfft,nspin), & 
                 deps_dN(nfft,nspin), & 
                 fermi(nspin), &      ! fermi levels for cluster and env
                 rhor(nfft,nspin)     ! subsystem electron density 

 ! local vars ............
 character(len=500) :: fname,ss,sc, string
 integer :: isp, itmp, ispin, ib


 ! function .....................

 write(ss,*)j_atom
 if (sys==1) fname = 'subsys'//trim(adjustl(ss))//'/result.dat'
 if (sys==2) fname = 'subsys'//trim(adjustl(ss))//'_env/result.dat'

! write(6,'(a,a)')'get_data: reading file:  ', trim(fname)
 open(file=fname,unit=111,action='read',form='formatted')
 read(111,*) etotal   ! this is actually ke + nlpsp - TS + \int (bare_vks*rhor)
 read(111,*) ke
 read(111,*) nlpsp
 read(111,*) TS
 close(111)

 ! load density 
 if (sys==1) fname = 'subsys'//trim(adjustl(ss))//'/rho.dat'
 if (sys==2) fname = 'subsys'//trim(adjustl(ss))//'_env/rho.dat'
 open(file=fname,unit=111,action='read',form='unformatted')
 read(111) rhor(:,1)
 if ( nspin==2) then 
   read(111) rhor(:,2)
 endif
 close(111)

 ! load eigenvalues of the subsystem 
 if (sys==1) fname = 'subsys'//trim(adjustl(ss))//'/eigen.dat'
 if (sys==2) fname = 'subsys'//trim(adjustl(ss))//'_env/eigen.dat'
 open(file=fname,unit=111,action='read',form='formatted')
 read(111,*) string, nband, string, itmp, string
 !
 occ = 0.d0 
 eigenvalues= 0.d0
 ! read the eigen.dat file
 do ispin=1,nspin 
   do ib=1,nband 
     read(111,*) string, itmp, itmp, itmp, occ(ib,ispin), eigenvalues(ib,ispin)
   enddo 
 enddo
 close(111)


 ! read fermi level 
 if (sys==1) fname = 'subsys'//trim(adjustl(ss))//'/fermi.dat'
 if (sys==2) fname = 'subsys'//trim(adjustl(ss))//'_env/fermi.dat'
 open(file=fname,unit=111,action='read',form='unformatted')
 read(111) fermi
 close(111)


! if (sys==1) fname = 'subsys'//trim(adjustl(ss))//'/dfermi_dvks.dat'
! if (sys==2) fname = 'subsys'//trim(adjustl(ss))//'_env/dfermi_dvks.dat'
! open(file=fname,unit=111,action='read',form='unformatted')
! read(111) dfermi_dvks 
! close(111)

 if (sys==1) fname = 'subsys'//trim(adjustl(ss))//'/drho_dN.dat'
 if (sys==2) fname = 'subsys'//trim(adjustl(ss))//'_env/drho_dN.dat'
 open(file=fname,unit=111,action='read',form='unformatted')
 read(111) drho_dN
 close(111)


 ! only load deps_dN for cluster 
 if (load_dEdVks==1 .and. sys==1) then 
   write(ss,*)j_atom
   open(file='subsys'//trim(adjustl(ss))//'/deps_dN.dat',unit=111,action='read',form='unformatted')
   read(111) deps_dN
   close(111)
 endif


 ! ======= load dEdV for XCEP ==========
 if (load_dEdVks==1 .and. sys==1) then 
   write(ss,*)j_atom
   fname = 'subsys'//trim(adjustl(ss))//'/yvec.dat'
   open(file=fname,unit=111,action='read',form='unformatted')
   read(111) exc_exx
   read(111) cluster_exc_eps
   read(111) yvec
   close(111)
   write(9944,'(a,i3)')'get yvec.dat from cluster: ',j_atom
 endif 

 ! 
 !
 !
! print *,''
! write(6,'(a,i2,a,i3)')'subsystem ',sys
! print *,'total energy is defined as ke + nlpsp - TS + \int(rho*bare_vks)'
! write(6,'(a,es20.12,a,i3,a,i3)')  '  etotal   = ', etotal, ' for sub: ',sys
! write(6,'(a,es20.12,a)')'  kinetic  = ', ke
! write(6,'(a,es20.12)')  '  TS       = ', TS
! write(6,'(a,es20.12)')  '  nlpsp    = ', nlpsp
! if (nspin==1) then 
!   write(6,'(a,f10.6)')  '  Q        = ', dvol*sum(rhor(:,1))
!   write(6,*) ' min/max(rho): ',minval(rhor),maxval(rhor)
! else
!  print *,''
!   !
!   ! display information
!   !
!   write(6,'(a,f10.6,a,f12.6)')'Q[alpha] = ', dvol*sum(rhor(:,1))
!   write(6,'(a,f10.6,a,f12.6)')'Q[beta]  = ', dvol*sum(rhor(:,2))
!   write(6,'(a,2e16.6)')'min/max(rho_alpha): ',minval(rhor(:,1)),maxval(rhor(:,1))
!   write(6,'(a,2e16.6)')'min/max(rho_beta):  ',minval(rhor(:,2)),maxval(rhor(:,2))
! endif

! print *,'left get_subsystem_data_xcpatch()'
! print *,''
! call flush(6)

 return 

end subroutine get_subsystem_data_xcpatch
