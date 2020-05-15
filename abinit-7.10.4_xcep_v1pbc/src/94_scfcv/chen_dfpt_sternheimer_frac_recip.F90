!
!  compute orbital perb (NOT orbital shift) for v_perb, 
!  by solving a modified (DFPT based) Sternheimer equation 
!  This code works fractional occupations
!
subroutine chen_dfpt_sternheimer_frac_recip(isp,iband,fermi,eigen,occ,occ_tol, & 
       n_homo,e_window,cg_fixphase,wfr,vlocal,orb_perb,vperb,npwarr,kg,dtset,ucvol,psps,xred,rprimd, & 
       rmet,gmet,gprimd,ylm,ylmgr,dtfil,atindx,atindx1,mpi_enreg,nattyp)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters
 use m_cgtools
 use m_hamiltonian
 use m_pawcprj
 use m_pawtab 
 use m_paw_ij

 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_65_nonlocal
 use interfaces_66_wfs
 use interfaces_67_common
 use interfaces_56_xc, only: rhohxc
 use m_fock,           only : fock_type, fock_get_getghc_call

 implicit none 

 ! ----------- external varaibles ----------------

 type(dataset_type),intent(in)         :: dtset
 type(datafiles_type),intent(in)       :: dtfil
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),intent(inout)          :: mpi_enreg
 integer, intent(in)                   :: isp,npwarr(dtset%nkpt),atindx(dtset%natom),iband,n_homo, & 
                                          atindx1(dtset%natom),nattyp(psps%ntypat),kg(3,npwarr(1))
 real(8), intent(in)                   :: occ_tol, & 
                                          e_window ! used to define subspace for removing the singularity of (H-e)
 real(kind=8),intent(inout)            :: vlocal(dtset%nfft)
 real(kind=8),intent(in)               :: rprimd(3,3),rmet(3,3),ucvol,gmet(3,3),gprimd(3,3),xred(3,dtset%natom), & 
                                          cg_fixphase(2,dtset%mpw,dtset%nband(1),dtset%nkpt,dtset%nspden), & ! wave function in the plane-wave basis set  
                                                                                                             ! I explicitely reshape it and fixed its phase, 
                                                                                                             ! wave function is real in the space
                                          ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm), & 
                                          ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm), & 
                                          vperb(dtset%nfft), & 
                                          wfr(dtset%nfft,dtset%nband(1)), & 
                                          fermi, & 
                                          occ(dtset%nband(1)), & 
                                          eigen(dtset%nband(1))  
 real(kind=8), intent(out)             :: orb_perb(dtset%nfft)

 ! ------------ local vars -----------------------

 character(len=1500) :: message
 type(fock_type),pointer   :: fock
 integer                   :: ia, iatom, ik, matblk, ndegen, & 
                              ider, idir, dimffnl, nkpg, ib, & 
                              nband_occupied   ! number of the band occupied 

 type(paw_ij_type)         :: paw_ij(dtset%natom)
 type(pawtab_type)         :: pawtab(psps%ntypat*psps%usepaw)
 type(gs_hamiltonian_type) :: gs_hamk
 type(pawcprj_type)        :: cwaveprj(dtset%natom,0)  ! usepaw=0, the last dim is 0

 real(kind=8)              :: dvol,tsmear, & 
                              shift=0.2d0,  &  ! remove the sigularity of the sternheimer equation 
                              dotr, doti, dotr1 ,dotr2, doti1, doti2, & 
                              factor 

 real(kind=8)              :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom), &
                              gsqcut,boxcut,res2,vres2,vxcavg, & 
                              vlocal_getghc(dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),1), &
                              kinpw(npwarr(1)),kpoint(3),arg,eval,& 
                              ghc(2,npwarr(1)),gvnlc(2,npwarr(1)), & 
                              gsc_dummy(2,npwarr(1)*dtset%nspinor*1*((0+1)/2))

 real(kind=8),allocatable  :: ffnl(:,:,:,:), & 
                              kpg_k(:,:), & 
                              ph3d(:,:,:)

 ! For solving the AX=b, with CG 
 ! transformation from recip space to real space
 logical      :: init_cg_x_zero
 integer      :: cg_step, fft_dir
 real(kind=8) :: cg_x(2,npwarr(1)), &       ! guess vector for solving AX=b  with CG
                 ghc_real(dtset%nfft),  &   ! ghc is transformed from recip space to real space
                 cg_resid(2,npwarr(1)), &    ! residual vector in solving Ax=b with CG 
                 cg_resid_old(2,npwarr(1)), &! residual vector in solving Ax=b with CG          
                 cg_p(2,npwarr(1)) , &       ! the p vector in CG algrithm on wiki
                 cg_Ap(2,npwarr(1)), &       ! A*p in the CG algrithm
                 cg_b(2,npwarr(1)), &        ! the vector b in Ax=b
                 cg_b_real(dtset%nfft), &    ! the vector b in Ax=b
                 wfg(2,npwarr(1)), &         ! 
                 vperb_wf_g(2,npwarr(1)), &  ! 
                 wfg_tmp(2,npwarr(1)), &     ! 
                 cg_beta, &                 ! beta parameter in CG solver
                 cg_alpha, &                ! alpha parameter in CG solver
                 aa, dotv, &                ! alpha_k in Eq. 72 of Rev. Mod. Phys. 73, 515
                 xx, rr, beta, en, em

 nullify(fock)
 call init_ham()

 dvol = ucvol/dtset%nfft
 tsmear = dtset%tsmear

 if (dtset%nsppol==1) factor = 2.d0 
 if (dtset%nsppol==2) factor = 1.d0 

 write(*,*)''
 write(*,*)'------ get in chen_dfpt_sternheimer_frac_recip() -------'
 write(*,'(a,i4)')"doing band: ",iband
 write(6,'(a,f18.10)')'incoming eigevalue :', eigen(iband)
 write(*,'(a,f12.4)')"e_window: ",e_window
 write(6,'(a,f12.4)')'tsmear: ',tsmear
 write(6,'(a,f12.4)')'fermi: ',fermi
 write(*,'(a,2es12.4)') "vperb: ",minval(vperb),maxval(vperb)
 call flush(6)


 ! we first convert cg_p to recip space 
! call recip_to_real_NEW(-1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
!   dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),wfg,wfr(:,iband),.false.)  
 wfg = cg_fixphase(:,:,iband,1,isp)

 ! we first convert uxc to recip space 
 call recip_to_real_NEW(-1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
   dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),vperb_wf_g,vperb*wfr(:,iband),.false.)  

 ! convert vlocal to the vlocal_getghc for getghc() subroutine
 call fftpac(1,mpi_enreg,dtset%nspden,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3), &
     dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),dtset%ngfft,vlocal,vlocal_getghc,2)




 !===============================
 !          check codes 
 !===============================
 if ( abs(theta(fermi,eigen(iband),tsmear) - occ(iband)/factor)>1e-4) then 
    print *,'abs(theta(fermi,eigen(iband),tsmear) - occ(iband))>1e-4 in chen_dfpt_sternheimer_frac.F90'
    print *,'theta(fermi,eigen(iband),tsmear)',theta(fermi,eigen(iband),tsmear)
    print *,'occ(iband)/factor: ',occ(iband)/factor
    stop 
 endif 
 ! no need to normalize the dotr 
 call dotprod_g(dotr,doti,dtset%istwfk(1),npwarr(1),2,wfg,wfg,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
 if ( abs(dotr-1.d0)>1e-5 ) then 
    print *,'error <wf,wf>/= 1'
    stop 
 endif
 call getghc(-1,wfg,cwaveprj,dimffnl,ffnl,dtfil%filstat,ghc,gsc_dummy,gs_hamk,gvnlc,kg,&
&    kinpw,eval,mpi_enreg,dtset%natom,1,npwarr(1),dtset%nspinor,dtset%paral_kgb,ph3d, & 
     dtset%prtvol,0,0,0,vlocal_getghc,fock)
 call dotprod_g(dotr,doti,dtset%istwfk(1),npwarr(1),2,wfg,ghc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
 print *,'<wf_g|H|wf_g>: ',dotr
 if ( abs(dotr-eigen(iband))>1e-5 ) then 
    print *,'error <wf,wf>/= incoming eigenvalue'
    stop 
 endif 



 !======================
 !      make cg_b
 !======================

 cg_b = -occ(iband)/factor*vperb_wf_g
 !
 ! apply P operator Eq.(72) in Rev. Mod. Phys. 73, 515
 !
 do ib=1,dtset%nband(1)
 !do ib=1,n_homo
   en = eigen(iband)
   em = eigen(ib)

!   call recip_to_real_NEW(-1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
!     dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),wfg_tmp,wfr(:,ib),.false.)  
   wfg_tmp = cg_fixphase(:,:,ib,1,isp)

   call dotprod_g(dotr,doti,dtset%istwfk(1),npwarr(1),2,wfg_tmp, & 
                 vperb_wf_g,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   dotv = dotr 
   !!dotv = sum(wfr(:,ib)*vperb*wfr(:,iband))*dvol

   ! first part 
   cg_b =  cg_b + (occ(iband)/factor*theta(en,em,tsmear)+occ(ib)/factor*theta(em,en,tsmear))*wfg_tmp*dotv
   
   ! second part 
   ! We only consider the bands <  e + e_window
   !
   ! Note that A must positive definite to make CG work, so all the bands lower than 
   ! e + e_window will be considered 
   !
   if ( em < en + e_window) then 
     write(6,*)'added band:',ib
     if (ib/=iband) then 
       rr = (theta(fermi,en,tsmear)-theta(fermi,em,tsmear))/(en-em)
     else 
       ! use the derivative of the Fermi-Dirac distribution 
       xx = (fermi-em)/tsmear
       rr = -1.d0/tsmear*0.5d0*1.d0/(1.d0+cosh(xx))
     endif 
     ! by setting aa below, (H+Q-e) is always non-singular
     aa = en - em + shift
     cg_b = cg_b + aa*rr*theta(em,en,tsmear)*wfg_tmp*dotv
   endif 
 enddo
 cg_resid  = cg_b      ! since cg_x = 0.d0
 cg_p      = cg_resid

 ! convert cg_b to real space
 call recip_to_real_NEW(1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
   dtset%nspinor,dtset%mgfft,dtset%nfft,gs_hamk%ucvol,gs_hamk%ngfft,kg,npwarr(1),cg_b,cg_b_real,.false.)  
 print *, 'cg_b_real:', minval(cg_b_real),maxval(cg_b_real)




 ! =============================================================
 !       cg loop (for solving ax=b for orbital shift)
 ! =============================================================
 cg_x = 0.d0  
 print *,'cg_step   resid '
 cg_step = 0

 do while (.true.) 
   cg_step = cg_step + 1 

   ! compute residual 
   call dotprod_g(dotr1,doti1,dtset%istwfk(1),npwarr(1),2, & 
     cg_resid,cg_resid,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

   if ( mod(cg_step-1,1)==0 ) then 
     write(6,'(i5,1x,es12.4)')cg_step,dotr1
     call flush(6)
   endif
   if (dotr1 < 1e-12) then 
     write(6,*)'|cg_resid| < 1e-12.'
     call flush(6)
     exit
   endif
   if (cg_step>=10 .and. dotr1 < 1e-9) then 
     write(6,*)'|cg_resid| converged to 1e-9.'
     call flush(6)
     exit
   endif

   ! now compute <g|H|cg_p>
   call getghc(-1,cg_p,cwaveprj,dimffnl,ffnl,dtfil%filstat,cg_Ap,gsc_dummy,gs_hamk,gvnlc,kg,&
&    kinpw,eval,mpi_enreg,dtset%natom,1,npwarr(1),dtset%nspinor, & 
     dtset%paral_kgb,ph3d,dtset%prtvol,0,0,0,vlocal_getghc,fock)

   ! A|cg_p> = (H - ei)|cg_p>
   cg_Ap = cg_Ap - eigen(iband)*cg_p

   ! Apply the Q operator in Eq. 72 in Rev. Mod. Phys. 73, 515
   ! basically speakinng, we add in all the wave functions that are within 
   ! the small window centering at eigen(iband)
   !
   do ib=1,dtset%nband(1)
   !do ib = 1, n_homo
     ! Note that A must positive definite to make CG work, so all the bands lower than 
     ! e + e_window will be considered 
     !
     if (eigen(ib)<eigen(iband)+e_window) then 
       aa = eigen(iband)-eigen(ib) + shift

       ! convert wfr(:,ib) to recip space 
       !call recip_to_real_NEW(-1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
       !  dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),wfg_tmp,wfr(:,ib),.false.)  
       wfg_tmp = cg_fixphase(:,:,ib,1,isp)
       
       ! compute <wfr,cg_p>
       call dotprod_g(dotr,doti,dtset%istwfk(1),npwarr(1),2,cg_p,wfg_tmp,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
       cg_Ap = cg_Ap + aa*wfg_tmp*dotr
     endif 
   enddo

   call dotprod_g(dotr1,doti1,dtset%istwfk(1),npwarr(1),2,cg_resid,cg_resid,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   call dotprod_g(dotr2,doti2,dtset%istwfk(1),npwarr(1),2,cg_p,cg_Ap,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

   cg_alpha = dotr1 / dotr2

   if (cg_alpha<0) then 
     write(message,'(a)') 'WARNING!!!! chen_dfpt_sternheimer_frac_recip() => cg_alpha<0, exit dfpt()'
     call c_wrtlog(message)
     write(message,'(a,2e12.4)') 'chen_dfpt_sternheimer_frac_recip() => <cg_resid,cg_resid>: ',dotr1,doti1
     call c_wrtlog(message)
     write(message,'(a,2e12.4)') 'chen_dfpt_sternheimer_frac_recip() => <cg_p,cg_Ap>: ',dotr2,doti2
     call c_wrtlog(message)
     exit 
   endif 

   cg_x     = cg_x + cg_alpha * cg_p

   cg_resid_old = cg_resid
   cg_resid     = cg_resid - cg_alpha * cg_Ap
 
   call dotprod_g(dotr1,doti1,dtset%istwfk(1),npwarr(1),2,cg_resid,cg_resid,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   call dotprod_g(dotr2,doti2,dtset%istwfk(1),npwarr(1),2,cg_resid_old,cg_resid_old,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

   cg_beta = dotr1 / dotr2
   cg_p    = cg_resid + cg_beta * cg_p

 enddo
 ! ========== end of cg loop =======================


 call recip_to_real_NEW(1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
   dtset%nspinor,dtset%mgfft,dtset%nfft,gs_hamk%ucvol,gs_hamk%ngfft,kg,npwarr(1),cg_x,orb_perb,.false.)  

 print *,'orb_perb: ',minval(orb_perb),maxval(orb_perb)
 print *,'chen_dfpt_sternheimer_frac_recip() done! '
 print *,''
 call flush(6)

 deallocate(ffnl,kpg_k,ph3d)






contains 

  ! -------------------------------------------------------
  ! initialize all abinit routines for following operations
  ! -------------------------------------------------------
  subroutine init_ham()
    !
    kpoint = 0.0d0           ! we have only one kpoint at this moment

    ! initialize gs_hamk for getghc()
    call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nspden, & 
                          dtset%natom,dtset%ntypat,dtset%typat,xred,dtset%nfft,dtset%mgfft, & 
                          dtset%ngfft,rprimd,dtset%nloalg)

    gs_hamk%istwf_k = dtset%istwfk(1)
    gs_hamk%npw     = npwarr(1)
    gs_hamk%kpoint  = 0.d0              ! only consider one kpoint in the 
    call sphereboundary(gs_hamk%gbound,dtset%istwfk(1),kg,dtset%mgfft,npwarr(1))

    ! make ffnl, nonlocal form factors, for each type of atom up to ntypat
    ! and for each angular momentum.
    ider=0;idir=0;dimffnl=1
    nkpg=3*dtset%optforces*dtset%nloalg(5)  ! for norm-conserving nloalg(5)=0, and nkpg=0
    allocate(ffnl(npwarr(1),dimffnl,psps%lmnmax,dtset%ntypat))
    allocate(kpg_k(npwarr(1),nkpg))
    call mkkpg(kg,kpg_k,kpoint,nkpg,npwarr(1))  ! make kpg_k
    ! make ffnl
    call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
        gmet,gprimd,ider,idir,psps%indlmn,kg,kpg_k,kpoint,psps%lmnmax,&
        psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
        npwarr(1),dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
        !!psps%usepaw,psps%useylm,ylm_k,ylmgr)
        psps%usepaw,psps%useylm,ylm,ylmgr)  ! we only have one kpoint, replace ylm_k with ylm
 
    ! make kinpw
    call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg,kinpw,kpoint,npwarr(1))

    ! make the rest of gs_hamk 
    ! compute phkxred and eventually ph3d.
    do ia = 1,dtset%natom
      iatom = atindx(ia)
      kpoint = 0.d0
      arg = two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
      gs_hamk%phkxred(1,iatom)=cos(arg)
      gs_hamk%phkxred(2,iatom)=sin(arg)
    end do
 
    ! compute ph3d and ph1d
    matblk = dtset%natom
    allocate(ph3d(2,npwarr(1),matblk))
    call getph(atindx,dtset%natom,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),ph1d,xred)
    call ph1d3d(1,dtset%natom,kg,matblk,dtset%natom,npwarr(1), & 
            dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),gs_hamk%phkxred,ph1d,ph3d)

    ! print *,'ffnl=',minval(ffnl),maxval(ffnl),ffnl(1,1,1,1),sum(ffnl)
    ! print *,'ph1d=',minval(gs_hamk%ph1d),maxval(gs_hamk%ph1d)
    ! print *,'ph3d=',minval(ph3d),maxval(ph3d),sum(ph3d)
    ! print *,'kinpw=',minval(kinpw),maxval(kinpw),sum(kinpw)
    ! print *,'gs_hamk%gbound=',gs_hamk%gbound
    ! stop
 
    !
    ! You need to check very carefully in vtorho() subroutine and in init_hamiltonian() 
    ! to see how the gs_hamk is initialized in ABINIT.
    ! gs_hamk delivers all the information to the getghc() surboutine, very important!
    ! if nonlocal psp is used, gs_hamk%ekb will be initalized with init_hamiltonian()
    !
    gs_hamk%matblk = matblk
    gs_hamk%typat  = dtset%typat
    call getcut(boxcut,dtset%ecut,gmet,gsqcut,dtset%iboxcut,std_out,kpoint,dtset%ngfft)


  end subroutine init_ham


  !
  ! Fermi-Dirac occupation numbers 
  ! for details, see the text below Eq. 69 in Baroni et al. 
  ! "Phonons and related properties of extended systems from density-functional perturbation theory", Rev. Mod. Phys. 73, 515
  !
  function theta(e1,e2,tsmear)
     implicit none 
     real(8) :: theta
     real(8) :: e1, e2, tsmear ! two energy levels 
     ! fermi-dirac smearing 
     theta = 1.d0/(1.d0+exp(-(e1-e2)/tsmear))
  end function theta 



end subroutine chen_dfpt_sternheimer_frac_recip
