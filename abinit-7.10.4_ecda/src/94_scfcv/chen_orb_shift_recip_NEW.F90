!!
!!  compute orbital shift for a band and a spin
!!
!! NOTE: it is STRONLGY recommend that you solve the orbital shift equaiton (that is sterheimer equation)
!! in the reciprocal space. If you do it in real space, the back and forth FFT will cause 
!! small errors in the gradient and finally will give you the wrong solution to the Ax=b equation
!! I spent nearly three days to trace down this problem. So do it in reciprocal space! 
!!  
!!
subroutine chen_orb_shift_recip_NEW(isp,iband,eigen,fermi,cg_fixphase,wfr,vlocal,orb_shift,uxc,vxc,npwarr,& 
                          kg,dtset,ucvol,psps,xred,rprimd,rmet,gmet,gprimd, & 
                          ylm,ylmgr,dtfil,atindx,atindx1,mpi_enreg,nattyp)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters
 use m_hamiltonian
 use m_pawcprj
 use m_pawtab 
 use m_paw_ij
 use m_cgtools   ! for dotprod_g

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
 integer, intent(in)                   :: isp, &    ! spin channel to work on 
                                          iband, &  ! band to work on 
                                          npwarr(dtset%nkpt),atindx(dtset%natom), & 
                                          atindx1(dtset%natom),nattyp(psps%ntypat),kg(3,npwarr(1))
 real(kind=8), intent(in)              :: cg_fixphase(2,dtset%mpw,dtset%nband(1),dtset%nkpt,dtset%nspden)  ! wave function in the plane-wave basis set  
                                                                                                           ! I explicitely reshape it and fixed its phase, 
                                                                                                           ! wave function is real in the space.
 real(kind=8),intent(inout)            :: rprimd(3,3),rmet(3,3),ucvol,gmet(3,3),gprimd(3,3),xred(3,dtset%natom), & 
                                          ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm), & 
                                          ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm), & 
                                          vlocal(dtset%nfft), & 
                                          uxc(dtset%nfft), & 
                                          vxc(dtset%nfft), & 
                                          wfr(dtset%nfft,dtset%nband(1)), & 
                                          eigen(dtset%nband(1))                                          
 real(kind=8), intent(out)             :: orb_shift(dtset%nfft)
 real(8)                               :: fermi

 ! ------------ local vars -----------------------

 type(fock_type),pointer   :: fock
 integer                   :: ia, iatom, ik, ndegen, & 
                              ider, idir, dimffnl, nkpg, ib, rss
 type(paw_ij_type)         :: paw_ij(dtset%natom)
 type(pawtab_type)         :: pawtab(psps%ntypat*psps%usepaw)
 type(gs_hamiltonian_type) :: gs_hamk
 type(pawcprj_type)        :: cwaveprj(dtset%natom,0)  ! usepaw=0, the last dim is 0

 real(8)                   :: dvol 
 real(kind=8)              :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom), &
                              e_window = 0.1d0, & 
                              e_thr = 0.1d0, & 
                              dtmp, & 
                              bar_uxc_ij, & 
                              gsqcut,boxcut,res2,vres2,vxcavg, & 
                              vlocal_getghc(dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),1), &
                              kinpw(npwarr(1)),kpoint(3),arg,eval,bar_uxc,bar_vxc, & 
                              ghc(2,npwarr(1)),gvnlc(2,npwarr(1)), & 
                              gsc_dummy(2,npwarr(1)*dtset%nspinor*1*((0+1)/2))

 real(kind=8),allocatable  :: ffnl(:,:,:,:), & 
                              kpg_k(:,:), & 
                              ph3d(:,:,:)

 ! For solving the AX=b, with CG 
 ! transformation from recip space to real space
 logical      :: init_cg_x_zero
 integer      :: cg_step, fft_dir, i
 real(kind=8) :: cg_x(2,npwarr(1)), &        ! guess vector for solving AX=b  with CG
                 dotr1, dotr2, doti1, doti2, dotr, doti,& 
                 wfg(2,npwarr(1)), &
                 wfg_tmp(2,npwarr(1)), &
                 wfg_tmp_conj(2,npwarr(1)), &
                 uxc_g(2,npwarr(1)), &
                 cg_Ax(2,npwarr(1)),  &   ! ghc is transformed from recip space to real space
                 cg_resid(2,npwarr(1)), &    ! residual vector in solving Ax=b with CG 
                 cg_resid_old(2,npwarr(1)), &! residual vector in solving Ax=b with CG          
                 cg_p(2,npwarr(1)) , &       ! the p vector in CG algrithm on wiki
                 cg_alpha, &                ! alpha parameter in CG solver
                 cg_Ap(2,npwarr(1)), &       ! A*p in the CG algrithm
                 cg_beta, &                 ! beta parameter in CG solver
                 cg_b(2,npwarr(1)), &        ! the vector b in Ax=b
                 cg_b_real(dtset%nfft), &        ! the vector b in Ax=b
                 fftg_out(2,npwarr(1))      ! store the output ( recip space) from recip_to_real()

 nullify(fock)

 write(*,*)'>>>> enter chen_orb_shift_recip_NEW() <<<<'
 print *,'fermi:  ',fermi
 print *,'doing isp: ',isp
 write(*,'(a,i4)')"doing band: ",iband
 write(6,'(a,f18.10)')'incoming eigevalue :', eigen(iband)
 print *,'uxc:    ',minval(uxc),maxval(uxc)
 print *,'vlocal: ',minval(vlocal),maxval(vlocal)
 call flush(6)

 call init_abinit_arrays()
 dvol = ucvol/dtset%nfft

! call recip_to_real_NEW(-1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
!   dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),wfg,wfr(:,iband),.false.) 
!
! NOTE: 
! you cannot directly use the cg from ABINIT, but has to use this phase-fixed wave funtion cg_fixphase
! why? becauase uxc is computed using the orbitals whose phase are fixed. if you use the ABINIT's cg
! whose phases are not fixed. the results will be wrong.
!
 wfg = cg_fixphase(:,:,iband,1,isp)  ! consider only one kpoint at this point of time 

 !
 ! we first convert uxc to recip space, 
 ! NOTE: uxc already contains wf
 !
 call recip_to_real_NEW(-1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
   dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),uxc_g,uxc,.false.)  

 ! convert vlocal to the vlocal_getghc for getghc() subroutine
 call fftpac(1,mpi_enreg,dtset%nspden,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3), &
     dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),dtset%ngfft,vlocal,vlocal_getghc,2)


 !================
 ! check the code
 !================
 ! no need to normalize the dotr 
 call dotprod_g(dotr,doti,dtset%istwfk(1),npwarr(1),2,wfg,wfg,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
 if ( abs(dotr-1.d0)>1e-5 ) then 
    print *,'error <wf,wf>/= 1'
    print *,'<wf,wf>: ',dotr,doti
    stop 
 endif
 if (npwarr(1)/=dtset%mpw) then
   print *,'npwarr(1)/=dtset%mpw, they are'
   print *,npwarr(1),dtset%mpw
   stop
 endif 

 call getghc(-1,wfg,cwaveprj,dimffnl,ffnl,dtfil%filstat,ghc,gsc_dummy,gs_hamk,gvnlc,kg,&
&    kinpw,eval,mpi_enreg,dtset%natom,1,npwarr(1),dtset%nspinor,dtset%paral_kgb,ph3d, & 
     dtset%prtvol,0,0,0,vlocal_getghc,fock)
 
 call dotprod_g(dotr,doti,dtset%istwfk(1),npwarr(1),2,wfg,ghc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
 print *,'<wf_g|H|wf_g>: ',dotr
 if ( abs(dotr-eigen(iband))>1e-5 ) then 
    print *,'error <wf,H wf>/= incoming eigenvalue'
    print *,'<wf,H wf>: ',dotr,doti
    stop 
 endif 

 !=============
 !  make cg_b
 !=============
 cg_x = 0.d0  
 cg_b = uxc_g
 call dotprod_g(dotr,doti,dtset%istwfk(1),npwarr(1),2,wfg,uxc_g,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
 bar_uxc = dotr
 cg_b    = cg_b - bar_uxc*wfg

 !
 ! We do (psi is the orbital shift)
 !   [H_KS - e_i + G] \psi
 ! G is constructed to remove the signularity of [H_KS - e_i]
 !
 ! Define the G operator as 
 !   G = \sum_alpha |phi_alpha> <phi_alpha| (e_i - e_alpha + e_thr)
 !   with   e_i - e_window < e_alpha < e_i + e_window 
 !
 ! NOTE: for conjugate gradient to work, A must be positive, 
 ! therefore, the final window is e_alpha < e_i + e_winidow
 !
 ! Addition term on the RHS of orbital shift equation due to the G|psi>
 !
 print *,'cg_b ',minval(cg_b),maxval(cg_b),' before G operator'
 print *,'shift up all the bands below currrent band'
 do ib=1,dtset%nband(1)
   if (ib==iband) cycle 
 
   ! NOTE: for conjugate gradient to work, A must be positive, 
   ! therefore, the final window is e_alpha < e_i + e_winidow
   !
   if ( eigen(ib)<eigen(iband)+e_window) then 
     write(6,'(a,i4,a)')'+ Add |psi><psi|, ib:',ib,''
     dtmp = eigen(iband)-eigen(ib)

     !call recip_to_real_NEW(-1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
     !  dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),wfg_tmp,wfr(:,ib),.false.)  
     wfg_tmp = cg_fixphase(:,:,ib,1,isp)

     call dotprod_g(dotr,doti,dtset%istwfk(1),npwarr(1),2,wfg_tmp,uxc_g,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
     bar_uxc_ij = dotr

     cg_b = cg_b - wfg_tmp*bar_uxc_ij*(dtmp+e_thr)/dtmp
   endif 
 enddo
 print *,'G operator is applied to remove possible degeneracy'

 cg_resid  = cg_b      ! since cg_x = 0.d0
 cg_p      = cg_resid

 ! convert <g|H|cg_p> to real space
 call recip_to_real_NEW(1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
   dtset%nspinor,dtset%mgfft,dtset%nfft,gs_hamk%ucvol,gs_hamk%ngfft,kg,npwarr(1),cg_b,cg_b_real,.false.)  

 print *,'cg_b ',minval(cg_b_real),maxval(cg_b_real)
 write(6,'(a,2f12.6)')'bar_vxc, bar_uxc:', bar_vxc,bar_uxc
 print *,'cg_step   resid   cg_alpha  cg_beta'
 call flush(6)


 !==================================================
 ! CG loop (for solving Ax=b for orbital shift)
 !===================================================

 do cg_step = 1,10000

   if ( mod(cg_step-1,1)==0 ) then 
     call dotprod_g(dotr1,doti1,dtset%istwfk(1),npwarr(1),2,cg_resid,cg_resid,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
     write(6,'(i5,1x,es12.4,1x,2es12.4)') cg_step,dotr1,cg_alpha,cg_beta
     call flush(6)
   endif

   call dotprod_g(dotr2,doti2,dtset%istwfk(1),npwarr(1),2,cg_resid,cg_resid,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   if ( dotr2<1e-12 ) then 
     write(6,*)'|cg_resid|<1e-12, exit'
     exit 
   endif 
   if ( cg_step>=10 .and. dotr2 < 1e-6) then 
     write(6,*)'(chen_cgwf_orb_shift) cg_step > 10 and CG converged to 1e-6.'
     exit
   endif

   cg_resid_old = cg_resid

   ! compute A*p
   ! now compute <g|H|cg_p>
   !
   call getghc(-1,cg_p,cwaveprj,dimffnl,ffnl,dtfil%filstat,cg_Ap,gsc_dummy,gs_hamk,gvnlc,kg,&
&    kinpw,eval,mpi_enreg,dtset%natom,1,npwarr(1),dtset%nspinor,dtset%paral_kgb,ph3d, & 
     dtset%prtvol,0,0,0,vlocal_getghc,fock)

   ! obtain the complete A|cg_p> = (H - ei + |ei><ei|)|cg_p>
   cg_Ap = cg_Ap - eigen(iband)*cg_p

   ! Apply the G operator 
   do ib=1,dtset%nband(1)
     
     ! NOTE: for conjugate gradient to work, A must be positive, 
     ! therefore, the final window is e_alpha < e_i + e_winidow
     !
     if ( eigen(ib)<eigen(iband)+e_window ) then 
       ! convert wfr(:,ib) to recip space 
       !call recip_to_real_NEW(-1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
       !  dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),wfg_tmp,wfr(:,ib),.false.)  
       wfg_tmp = cg_fixphase(:,:,ib,1,isp)
       
       ! compute <wfr,cg_p>
       call dotprod_g(dotr,doti,dtset%istwfk(1),npwarr(1),2,wfg_tmp,cg_p,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

       wfg_tmp_conj(:,1) =  wfg_tmp(:,1)
       wfg_tmp_conj(:,2) = -wfg_tmp(:,2)

       cg_Ap = cg_Ap + wfg_tmp*dotr*(eigen(iband)-eigen(ib)+e_thr)
     endif 
   enddo 

   call dotprod_g(dotr1,doti1,dtset%istwfk(1),npwarr(1),2,cg_resid,cg_resid,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   call dotprod_g(dotr2,doti2,dtset%istwfk(1),npwarr(1),2,cg_p,cg_Ap,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

   cg_alpha = dotr1 / dotr2

   if (cg_alpha<0) then 
     print *,'<cg_resid,cg_resid>: ',dotr1,doti1
     print *,'<cg_p,cg_Ap>: ',dotr2,doti2
     print *,'cg_alpha<0, BUG!!'
     call flush(6)
     stop
   endif 

   cg_x     = cg_x + cg_alpha * cg_p
   cg_resid = cg_resid - cg_alpha * cg_Ap
 
   call dotprod_g(dotr1,doti1,dtset%istwfk(1),npwarr(1),2,cg_resid,cg_resid,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   call dotprod_g(dotr2,doti2,dtset%istwfk(1),npwarr(1),2,cg_resid_old,cg_resid_old,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

   cg_beta = dotr1 / dotr2
   cg_p    = cg_resid + cg_beta * cg_p

 enddo ! cg loop

 ! convert <g|H|cg_p> to real space
 call recip_to_real_NEW(1,mpi_enreg,dtset%paral_kgb,1,dtset%mpw,npwarr(1),0,dtset%istwfk(1), & 
   dtset%nspinor,dtset%mgfft,dtset%nfft,gs_hamk%ucvol,gs_hamk%ngfft,kg,npwarr(1),cg_x,orb_shift,.false.)  

 ! project out the current KS wave function
 orb_shift = orb_shift - sum(orb_shift*wfr(:,iband))*dvol*wfr(:,iband)

 write(6,'(a,2es18.8)')'orbital shift: ',minval(orb_shift),maxval(orb_shift)

 deallocate(ffnl,kpg_k,ph3d)

 write(6,'(a)')'done chen_orb_shift_recip_NEW()'
 write(6,'(a)')''

contains 



  ! -------------------------------------------------------
  ! initialize all abinit routines for following operations
  ! -------------------------------------------------------
  subroutine init_abinit_arrays()
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
    allocate(ph3d(2,npwarr(1),gs_hamk%matblk))

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
    call getph(atindx,dtset%natom,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),ph1d,xred)
    call ph1d3d(1,dtset%natom,kg,gs_hamk%matblk,dtset%natom,npwarr(1), & 
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
    gs_hamk%typat  = dtset%typat
    call getcut(boxcut,dtset%ecut,gmet,gsqcut,dtset%iboxcut,std_out,kpoint,dtset%ngfft)

  end subroutine init_abinit_arrays

end subroutine 
