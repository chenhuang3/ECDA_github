
!
!  compute orbital shift for a band and a spin
!
subroutine chen_orb_shift_gmres(iband,eigen,fermi,wfr,vlocal,orb_shift,uxc,vxc,npwarr,& 
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
 integer, intent(in)                   :: npwarr(dtset%nkpt),atindx(dtset%natom),iband, & 
                                          atindx1(dtset%natom),nattyp(psps%ntypat),kg(3,npwarr(1))
 real(kind=8),intent(inout)            :: rprimd(3,3),rmet(3,3),ucvol,gmet(3,3),gprimd(3,3),xred(3,dtset%natom), & 
                                          ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm), & 
                                          ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm), & 
                                          vlocal(dtset%nfft),uxc(dtset%nfft),vxc(dtset%nfft), & 
                                          wfr(dtset%nfft,dtset%nband(1)),eigen(dtset%nband(1))                                          
 real(kind=8), intent(out)             :: orb_shift(dtset%nfft)
 real(8)                               :: fermi

 ! ------------ local vars -----------------------

 type(fock_type),pointer   :: fock
 integer                   :: ia, iatom, ik, matblk, ndegen, & 
                              ider, idir, dimffnl, nkpg, ib
 type(paw_ij_type)         :: paw_ij(dtset%natom)
 type(pawtab_type)         :: pawtab(psps%ntypat*psps%usepaw)
 type(gs_hamiltonian_type) :: gs_hamk
 type(pawcprj_type)        :: cwaveprj(dtset%natom,0)  ! usepaw=0, the last dim is 0

 real(8)                   :: dvol 
 real(kind=8)              :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom), &
                              e_window = 0.5d0, & 
                              e_thr = 0.5d0, & 
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
 integer      :: cg_step, fft_dir
 real(kind=8) :: cg_x(dtset%nfft), &        ! guess vector for solving AX=b  with CG
                 cg_Ax(dtset%nfft),  &   ! ghc is transformed from recip space to real space
                 ghc_real(dtset%nfft),  &   ! ghc is transformed from recip space to real space
                 cg_resid(dtset%nfft), &    ! residual vector in solving Ax=b with CG 
                 cg_resid_old(dtset%nfft), &! residual vector in solving Ax=b with CG          
                 cg_p(dtset%nfft) , &       ! the p vector in CG algrithm on wiki
                 cg_alpha, &                ! alpha parameter in CG solver
                 cg_Ap(dtset%nfft), &       ! A*p in the CG algrithm
                 cg_beta, &                 ! beta parameter in CG solver
                 cg_b(dtset%nfft), &        ! the vector b in Ax=b
                 fftg_out(2,npwarr(1))      ! store the output ( recip space) from recip_to_real()

 ! GMRES
 integer,parameter :: gmres_mgmres = 0, & 
                      gmres_j = 5
 integer :: gmres_maxits,gmres_n,gmres_iflag
 logical :: gmres_oktest
 character(len=3) :: gmres_stc
 real(8) :: gmres_work(dtset%nfft,0:2*gmres_j+gmres_mgmres+2-1), & 
            gmres_eps, gmres_resid

 nullify(fock)

 call init_abinit_arrays()
 dvol = ucvol/dtset%nfft

 write(*,*)''
 print *,'fermi:  ',fermi
 write(*,'(a,i4)')    "oep (chen_cgwf_orb_shift) doing band: ",iband
 write(6,'(a,f18.10)')'oep (chen_cgwf_orb_shift) incoming eigevalue :', eigen(iband)
 print *,'oep => uxc: ',minval(uxc),maxval(uxc)
 print *,'oep => vxc: ',minval(vxc),maxval(vxc)
 print *,'vlocal: ',minval(vlocal),maxval(vlocal)

 cg_x = 0.d0  
 cg_b = -(vxc-uxc)*wfr(:,iband)
 bar_vxc = dvol*sum(wfr(:,iband)*wfr(:,iband)*vxc)
 bar_uxc = dvol*sum(wfr(:,iband)*wfr(:,iband)*uxc)
 cg_b    = cg_b + (bar_vxc-bar_uxc)*wfr(:,iband)

 !
 ! We do (psi is the orbital shift)
 !   [H_KS - e_i + G] \psi
 ! G is constructed to remove the signularity of [H_KS - e_i]
 ! Define the G operator as 
 !   G = \sum_alpha |phi_alpha> <phi_alpha| (e_i - e_alpha + e_thr)
 !   with   e_i - e_window < e_alpha < e_i + e_window 
 !
 !
 ! Addition term on the RHS of orbital shift equation due to the G|psi>
 !
 print *,'cg_b ',minval(cg_b),maxval(cg_b),' before G operator'
 do ib=1,dtset%nband(1)
   if (ib==iband) cycle 
   if (eigen(ib)>eigen(iband)-e_window  .and.  &  
       eigen(ib)<eigen(iband)+e_window) then 
     print *,'add |psi><psi|, ib:',ib
     dtmp = eigen(iband)-eigen(ib)
     bar_uxc_ij = dvol*sum(wfr(:,iband)*wfr(:,ib)*uxc)
     cg_b = cg_b - wfr(:,ib)*bar_uxc_ij*(dtmp+e_thr)/dtmp
     !!cg_b = cg_b - wfr(:,ib)*bar_uxc_ij*(fermi-eigen(ib)+e_thr)/(eigen(iband)-eigen(ib))
   endif 
 enddo
 print *,'G operator is applied to remove possible degeneracy'
 cg_resid  = cg_b      ! since cg_x = 0.d0
 cg_p      = cg_resid

 print *,'cg_b ',minval(cg_b),maxval(cg_b)
 write(6,'(a,2f12.6)')'bar_vxc, bar_uxc:', bar_vxc,bar_uxc
 print *,'cg_step   resid     min(cg_x)    max(cg_x)'
 call flush(6)


 gmres_stc = 'rel'
 gmres_maxits = 100
 gmres_n = dtset%nfft
 gmres_oktest = .true.
 gmres_eps = 1e-4

 call gmresr(gmres_oktest,gmres_n,gmres_j,gmres_mgmres,cg_b,cg_x, & 
    gmres_work,gmres_eps,gmres_stc,gmres_maxits,gmres_resid,chen_matvec,gmres_iflag)


 orb_shift = cg_x 

 !
 ! project out the current KS wave function
 ! this is a must, since current KS orbital is in the null space of (H-eigen)
 !
 orb_shift = orb_shift - sum(orb_shift*wfr(:,iband))*dvol*wfr(:,iband)

 print *,'oep <cg_x,wfr> ',sum(cg_x*wfr(:,iband))*dvol,' should be close to zero'
 print *,'oep orbital shift: ',minval(orb_shift),maxval(orb_shift)

 deallocate(ffnl,kpg_k,ph3d)
 write(6,'(a)')'(chen_cgwf_orb_shift) finished CG.'

contains 


  subroutine chen_matvec(x,Ax,nfft)

   integer :: nfft
   real(8) :: x(nfft), Ax(nfft)

   ! compute Ap
   ! we first convert cg_p to recip space 
   call recip_to_real(-1,mpi_enreg,dtset%paral_kgb,1,npwarr(1),0,dtset%istwfk(1), & 
     dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),fftg_out,x)  

   ! convert vlocal to the vlocal_getghc for getghc() subroutine
   call fftpac(1,mpi_enreg,dtset%nspden,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3), &
     dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),dtset%ngfft,vlocal,vlocal_getghc,2)
 
   ! now compute <g|H|cg_p>
   call getghc(-1,fftg_out,cwaveprj,dimffnl,ffnl,dtfil%filstat,ghc,gsc_dummy,gs_hamk,gvnlc,kg,&
&    kinpw,eval,mpi_enreg,dtset%natom,1,npwarr(1),dtset%nspinor, & 
     dtset%paral_kgb,ph3d,dtset%prtvol,0,0,0,vlocal_getghc,fock)

   ! convert <g|H|cg_p> to real space
   call recip_to_real(1,mpi_enreg,dtset%paral_kgb,1,npwarr(1),0,dtset%istwfk(1), & 
     dtset%nspinor,dtset%mgfft,dtset%nfft,ucvol,dtset%ngfft,kg,npwarr(1),ghc,Ax)  

   ! obtain the complete A|cg_p> = (H - ei + |ei><ei|)|cg_p>
   Ax = Ax - eigen(iband) * x

   ! Apply the G operator 
   do ib=1,dtset%nband(1)
     ! energy window 
     if (eigen(ib)>eigen(iband)-e_window .and. & 
         eigen(ib)<eigen(iband)+e_window) then 

       Ax = Ax + wfr(:,ib)*sum(dvol*wfr(:,ib)*x)*(eigen(iband)-eigen(ib)+e_thr)
     !!  cg_Ap = cg_Ap + wfr(:,ib)*sum(dvol*wfr(:,ib)*cg_p)*(fermi-eigen(ib)+e_thr)
     endif 
   enddo 

  end subroutine chen_matvec


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


  end subroutine init_abinit_arrays


end subroutine chen_orb_shift_gmres
