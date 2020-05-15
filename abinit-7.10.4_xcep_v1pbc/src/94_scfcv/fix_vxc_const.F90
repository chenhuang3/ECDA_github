!
! We fix the constant in the Vxc by making a variation of the electron number.
! the rate of dE_xc / dN is then equated to the \int V_xc * d rho /d N 
! see my notes and paper for more detail. This method is extermely simple 
! and is designed for ensemble (F-S smearing) system only. 
!
! Chen Huang, November 2016
!
subroutine fix_vxc_const(dexx_docc,chen_occ,nfft,nspin,n_homo,norb_extra,& 
                         nband,eigenvalues,fermi,phir,tsmear,dvol,vxc_kli)

  implicit none 

  ! >>>>>>>>>> external vars <<<<<<<<<<!
  integer :: nfft, nspin, nband, n_homo(nspin),norb_extra
  real(8) :: dexx_docc(nspin,nband), tsmear, & 
             vxc_kli(nfft,nspin), dvol, &
             chen_occ(nspin,nband), & 
             fermi(nspin), & 
             eigenvalues(nspin,nband),&
             phir(nfft,nspin,nband) 

  ! >>>>>>>> internal vars <<<<<<<<<<<<!
  character(len=200) :: message
  integer :: isp, ib1, ib2, & 
             method_docc_dN ! 1: docc_dN will be computed by performing analytic derivatives
                            ! 2: docc_dN will be computed using finite difference 

  real(8) :: vxc_const, sum_docc_dmu, & 
             fd_N , q, & 
             dExc_dN, &
             docc_dmu(nband), & 
             tol, denom, dtmp, & 
             docc_dN_fd(nband), &  ! compute using finite difference 
             docc_dN, vxc_bar , & 
             occ_tmp1(nband), occ_tmp2(nband),& 
             fermi_tmp1, fermi_tmp2, & 
             occi, occj, ei_i, ei_j, & 
             contrib

  logical :: second_order_ft = .false.

  call c_wrtlog("(fix_vxc_const) determine the constant based on the change of EXX energy w.r.t. the electron number. ")
  write(message,'(a,2f12.4)')'(fix_vxc_const) fermi: ',fermi(1:2)
  call c_wrtlog(message)
  if (.not. second_order_ft) then 
    write(message,'(a,2f12.4)')'dExc_dN is take from the electron-sufficient side => Ionization energy match.'
    call c_wrtlog(message)  
  endif

  !print *,'eigen: ',eigenvalues(1,:)
  !print *,'eigen: ',eigenvalues(2,:)

  method_docc_dN = 2

  do isp=1,nspin
      
      !
      ! Employ the finite difference way to compute docc_dN 
      ! one potential benifit of the computing docc_dN using finite difference
      ! method is that it is simple, 
      ! when the gap is large and temperature is low, the fermi energy can be 
      ! hard to determined, This finite difference method does not rely on fermi energy. 
      !
      if (method_docc_dN==2) then 

        q = sum(chen_occ(isp,:))

        if (second_order_ft) then 

          fd_N = 1e-3
          call chen_fermi_occ(tsmear, nband, eigenvalues(isp,:), q+fd_N, occ_tmp1, fermi_tmp1)
          call chen_fermi_occ(tsmear, nband, eigenvalues(isp,:), q-fd_N, occ_tmp2, fermi_tmp2)

          docc_dN_fd = (occ_tmp1-occ_tmp2)/2.d0/fd_N

          write(message,'(a,es12.4,a,i3)')'(fix_vxc_const) 2nd order fintie difference to compute docc_dN, q:',q,' for spin ',isp
          call c_wrtlog(message)
        else
          fd_N = 1e-3
          call chen_fermi_occ(tsmear, nband, eigenvalues(isp,:), q,      occ_tmp1, fermi_tmp1)
          call chen_fermi_occ(tsmear, nband, eigenvalues(isp,:), q-fd_N, occ_tmp2, fermi_tmp2)

          docc_dN_fd = (occ_tmp1-occ_tmp2)/fd_N

          write(message,'(a,es12.4,a,i3)')'(fix_vxc_const) same ionization energy to determine docc_dN, q:',q,' for spin ',isp
          call c_wrtlog(message)
        endif


      endif 
      
      !
      ! computing vxc_constant 
      !
      vxc_const = 0.0d0 
      dExc_dN = 0.0d0

      do ib1=1,n_homo(isp) + norb_extra


        if (method_docc_dN==1) then 
          !
          ! new way of computing docc_dN 
          !
          !    df_i/dN = f_i(1-f_i) / sum(f_j(1-f_j))
          !
          ! we do not follow the above expression in fact, we do the following 
          !
          !    df_i/dN = 1 / (1+sum_{j\neq i} f_j(1-f_j)/(f_i*(1-f_i)))
          !
          ! we consider the cases that f_i is very close to zero or one, which will make the 
          ! above denominator (f_i*(1-f_i)) to zero. Ok, the four caese are
          !
          !   1) f_i -> 0 and f_j -> 0
          !   2) f_i -> 0 and f_j -> 1
          !   3) f_i -> 0 and f_j is not close to 0 or 1
          !   4) f_i -> 1 and f_j -> 0
          !   5) f_i -> 1 and f_j -> 1
          !   6) f_i -> 1 and f_j is not close to 0 or 1
          !   7) both f_i and f_j are not close to 0 or 1
          ! 
          ! apparently, we cannot directly use the occupation numbers (f_i), since they will be 
          ! just 1 or 0. We need to work on the eigenvalues and chemical potential.  
          !
          ! note that for Fermi-Dirac statistics f=1/(1+x) with x=exp[(E-\mu)*beta]
          ! we have 
          !         for f->0 => x->\infinity => f ~ 1/x
          !         for f->1 => x->0         => f ~ x
          !
          denom = 1.0d0;

          do ib2 = 1, n_homo(isp)+norb_extra

             if (ib2==ib1) cycle
             ! 
             ! we consider four cases that f_i and f_j are close to 
             ! zero or one.
             !
             tol = 1e-6  

             occi = chen_occ(isp,ib1)  ! f_i 
             occj = chen_occ(isp,ib2)  ! f_j 
             ei_i = eigenvalues(isp,ib1) 
             ei_j = eigenvalues(isp,ib2)
             !
             ! the expression below is based on the expression denom = f_j(1-f_j) / (f_i(1-f_i)) 
             !
             if (occi<tol) then 
               !
               ! for f_i ~ 0
               !
               if (occj<tol) then 
                   dtmp =  ( ei_i - ei_j )/tsmear
                   denom = denom + exp(dtmp)

               elseif (abs(1.d0-occj)<tol ) then 
                   dtmp =  ( ei_i + ei_j - 2.d0*fermi(isp))/tsmear
                   denom = denom + exp(dtmp)

               elseif (occj>tol .and. abs(1.d0-occj)>tol) then 
                   dtmp = (ei_i - fermi(isp))/tsmear
                   denom = denom + exp(dtmp)*occj*(1.0d0-occj)
               endif 

             elseif ( abs(1.d0-occi)<tol) then 
               ! 
               ! for f_i ~ 1
               ! 
               if (occj<tol) then
                  dtmp =  -( ei_i + ei_j - 2.d0*fermi(isp))/tsmear
                  denom = denom + exp(dtmp)

               elseif (abs(1.d0-occj)<tol) then  
                  dtmp =  ( ei_j - ei_i )/tsmear
                  denom = denom + exp(dtmp)

               elseif (occj>tol .and. abs(1.d0-occj)>tol) then 
                  dtmp  = (ei_i - fermi(isp))/tsmear
                  denom = denom + exp(-dtmp)*occj*(1.d0-occj)
               endif
             else
               !
               ! for the case that f_i is neither 1 or 0
               !
               denom = denom + occj*(1.d0-occj)/(occi*(1.d0-occi))
             endif

           enddo ! loop over all the bands

           ! compute d f_i / d N
           docc_dN = 1.d0 / denom;
        endif 

        !
        ! finite difference way for computing docc_dN 
        !
        if ( method_docc_dN == 2) then 
           docc_dN = docc_dN_fd(ib1)
        endif

        !print *,'docc_dN : ',docc_dN
        dExc_dN = dExc_dN + dexx_docc(isp,ib1)*docc_dN
        vxc_bar = dvol*sum(phir(:,isp,ib1)**2*vxc_kli(:,isp))
        contrib = dexx_docc(isp,ib1)*docc_dN - vxc_bar*docc_dN
        vxc_const = vxc_const + contrib

        write(message,'(a,es12.4,a,es12.4,a,i3,a,i2,a,f12.4)') & 
        '   docc_dN ',docc_dN,'   contrib: ',contrib,' ib: ',ib1,' isp: ',isp,' <phi|vxc|phi>: ',vxc_bar
        call c_wrtlog(message)
        
        !
        ! the vxc_shift from the highest band considered here should be small enough
        !
        if ( ib1==n_homo(isp)+norb_extra .and. abs(contrib) > 1e-3) then
           print *,'the last contrib is too large!, contrib:',contrib
           print *,'n_homo: ',n_homo(isp)
           print *,'n_homo+norb_extra: ',n_homo(isp)+norb_extra
           stop 
        endif

       ! write(message,'(a,es12.4,a,i3,a,i3)')'(fix_vxc_const) vxc_bar ',vxc_bar,' for spin: ',isp,' ib:',ib1
       ! call c_wrtlog(message)
       ! write(message,'(a,es12.4,a,i3,a,i3)')'(fix_vxc_const) vxc_const ',vxc_const,' for spin: ',isp,' ib:',ib1
       ! call c_wrtlog(message)


      enddo  ! loop over all the bands 
  
      write(message,'(a,es12.4,a,i3,a,f12.4)')'==> vxc_const: ',vxc_const,'   spin: ',isp,' dExc_dN: ',dExc_dN
      call c_wrtlog(message)
      vxc_kli(:,isp) = vxc_kli(:,isp) + vxc_const

  enddo ! loop over spin 

  call c_wrtlog("(fix_vxc_const) vxc_kli is shifted based on EXX energy change")

end subroutine fix_vxc_const


