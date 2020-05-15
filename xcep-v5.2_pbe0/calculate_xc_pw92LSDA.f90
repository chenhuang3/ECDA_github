
!
! wrapper for libxc pw92LSDA
!

subroutine  calculate_xc_pw92LSDA(n1,n2,n3,nspin,rho,dvol,xc_pot,xc_energy)

  use xc_f90_types_m 
  use xc_f90_lib_m

  implicit none 

  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  integer :: func_id

  integer                 :: n1,n2,n3,ixc,nspin
  real(kind=8),intent(in) :: rho(n1*n2*n3,nspin), & 
                             dvol

  real(kind=8),intent(out) :: xc_pot(n1*n2*n3,nspin), & 
                              xc_energy

  ! for libxc 
  real(kind=8) :: exc(n1*n2*n3) ! The energy per unit particle array always has dimensions exc[np]

  ! If the functional was initialized with nspin=XC_UNPOLARIZED, 
  ! the spin indices should be dropped from the previous expressions, 
  ! resulting in arrays of size np. 
  ! Otherwise, the parameters have dimensions 
  ! rho[2*np], vxc[2*np], sigma[3*np], vsigma[3*np], 
  ! v2rho2[3*np], v2rhosigma[6*np], v2sigma2[6*np]. 

  real(kind=8) :: rho_tmp(nspin), &
                  exc_tmp(1),     & 
                  vrho_tmp(nspin)

  integer :: nfft, j

  nfft = n1*n2*n3
  xc_energy = 0.0d0
  xc_pot = 0.d0

  ! DEFINITIONS IN LIBXC 
  !  The quantity returned by libxc is the energy per unit particle, 
  !  However that the derivatives returned by libxc are 
  !  with respect to the energy per unit volume.

  ! ============================
  ! === NON-SPIN POLARIZED   ===
  ! ============================

  if (nspin==1) then 

    ! exchange ===========
    call xc_f90_func_init(xc_func,xc_info,XC_LDA_X,XC_UNPOLARIZED)
    do j=1,nfft
      rho_tmp(1) = rho(j,1)
      call xc_f90_lda_exc(xc_func,1,rho_tmp(1),exc_tmp(1))
      call xc_f90_lda_vxc(xc_func,1,rho_tmp(1),vrho_tmp(1))
      xc_energy   = xc_energy + exc_tmp(1)*rho_tmp(1)*dvol
      xc_pot(j,1) = xc_pot(j,1) + vrho_tmp(1)
      !
      exc(j) = exc_tmp(1)*rho_tmp(1)
    enddo
    call xc_f90_func_end(xc_func)

    ! correlation =============
    call xc_f90_func_init(xc_func,xc_info,XC_LDA_C_PW,XC_UNPOLARIZED)
    do j=1,nfft
      rho_tmp(1) = rho(j,1)
      call xc_f90_lda_exc(xc_func,1,rho_tmp(1),exc_tmp(1))
      call xc_f90_lda_vxc(xc_func,1,rho_tmp(1),vrho_tmp(1))
      xc_energy   = xc_energy + exc_tmp(1)*rho_tmp(1)*dvol
      xc_pot(j,1) = xc_pot(j,1) + vrho_tmp(1)
      !
      exc(j) = exc(j) + exc_tmp(1)*rho_tmp(1)
    enddo
    call xc_f90_func_end(xc_func)

  endif


  ! ============================
  ! === SPIN POLARIZED       ===
  ! ============================
  if (nspin==2) then 

    ! exchange 
    ! --------
    call xc_f90_func_init(xc_func,xc_info,XC_LDA_X,XC_POLARIZED)
    do j=1,nfft
      rho_tmp = rho(j,:)
      call xc_f90_lda_exc_vxc(xc_func,1,rho_tmp(1),exc_tmp(1),vrho_tmp(1))
      xc_energy   = xc_energy + exc_tmp(1)*(rho(j,1)+rho(j,2))*dvol
      xc_pot(j,1) = vrho_tmp(1)
      xc_pot(j,2) = vrho_tmp(2)
      !
      exc(j) = exc_tmp(1)*(rho_tmp(1)+rho_tmp(2))
    enddo
    call xc_f90_func_end(xc_func)

    ! correlation 
    ! -----------
    call xc_f90_func_init(xc_func,xc_info,XC_LDA_C_PW,XC_POLARIZED)
    do j=1,nfft
      rho_tmp = rho(j,:)
      call xc_f90_lda_exc_vxc(xc_func,1,rho_tmp(1),exc_tmp(1),vrho_tmp(1))
      xc_energy   = xc_energy + exc_tmp(1)*(rho(j,1)+rho(j,2))*dvol
      xc_pot(j,1) = xc_pot(j,1) + vrho_tmp(1)
      xc_pot(j,2) = xc_pot(j,2) + vrho_tmp(2)
      !
      exc(j) = exc(j) + exc_tmp(1)*(rho_tmp(1)+rho_tmp(2))
    enddo
    call xc_f90_func_end(xc_func)

  endif

  ! temperarily for xc patching 
  open(file='tmp_pw92_density.dat',unit=111,action='write',form='unformatted')
  write(111) exc
  close(111)
  
  print *,''
  print *,'   --------- calculate_xc_pw92LSDA.f90 -------------'
  print *,''
  if (nspin==1) & 
    write(6,'(a,f12.6)')   ' Q        : ',sum(rho(:,1))*dvol
  if (nspin==2) then 
    write(6,'(a,f12.6)')   ' Q_alpha  : ',sum(rho(:,1))*dvol
    write(6,'(a,f12.6)')   ' Q_beta   : ',sum(rho(:,2))*dvol
  endif
  write(6,'(a,f14.8,a)')   ' xc_energy: ',xc_energy,' Ha'
  if (nspin==1) then 
    write(6,'(a,2f14.8,a)')' xc_pot   : ',minval(xc_pot(:,1)),maxval(xc_pot(:,1)),' Ha'
  endif
  if (nspin==2) then 
    write(6,'(a,2f14.8,a)')' xc_pot   : ',minval(xc_pot(:,1)),maxval(xc_pot(:,1)),' Ha [alpha]'
    write(6,'(a,2f14.8,a)')' xc_pot   : ',minval(xc_pot(:,2)),maxval(xc_pot(:,2)),' Ha [beta]'
  endif
  
!  call LSDAPW92(n1,n2,n3,rho,xc_pot,xc_energy,dvol)
!  write(6,'(a,f14.8,a)') 'calculate_xc_pw92LSDA: xc_energy: ',xc_energy,' Ha'
!  write(6,'(a,2f14.8,a)')'calculate_xc_pw92LSDA: xc_pot   : ',minval(xc_pot(:,1)),maxval(xc_pot(:,1)), ' Ha'
!  if (nspin==2) &
!  write(6,'(a,2f14.8,a)')'calculate_xc_pw92LSDA: xc_pot   : ',minval(xc_pot(:,2)),maxval(xc_pot(:,2)),' [beta]'

end subroutine calculate_xc_pw92LSDA
