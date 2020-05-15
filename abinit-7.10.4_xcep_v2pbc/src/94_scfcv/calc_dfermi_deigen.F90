!
! Compute d(fermi level)/d(eigen) for a given band
! The formula for this is very simple, 
!
!  d(fermi_level) / d(eigen_m) = 1 / [1 + sum_{m/=n} fm(1-fm)/fn(1-fn)] 
!
! why do we need to write such a lengthy code to do it? 
! because, the occupation numbers can get very close to 1 and zero. 
! therefore making the denominator go to zero, therefore you get NaN.
!
! INPUT/OUTPUT
! ============
! iband:  the band to be computed
! isp: the spin channel to be computed
! nband: number of bands
! nspin: number of spin
!
! 12/27/2017 (Chen Huang)
!
subroutine  calc_dfermi_deigen(iband,isp,nband,nspin,tsmear,fermi,occ,eigen,dmu_deigen)

  implicit none 
                         
  integer :: method = 2  ! 1: a simple code just compute fn(1-fn)/sum(fm(1-fm)) 
                         ! 2: a lengthy code
  integer :: iband,nband,nspin,isp
  real(8) :: fermi(nspin), tsmear, & 
             eigen(nspin,nband), & 
             occ(nspin,nband), & 
             dmu_deigen

  ! local vars 
  integer :: ib 
  real(8) :: zz, fn, fm, en, em, expo, dsum
  real(8) :: tiny_value=1e-6,  & 
             log_fn_part, log_fm_part, & 
             factor 

 
  if (nspin==2) factor = 1.d0 
  if (nspin==1) factor = 2.d0 


  !=============================================
  !  a simple version of the code 
  !=============================================
  if (method==1) then 
     dsum = 0.d0 
     do ib=1,nband
        fm = occ(isp,ib)/factor 
        dsum = dsum + fm*(1.d0-fm)
     enddo 
     fn = occ(isp,iband) / factor 
     dmu_deigen = fn*(1.d0-fn)/dsum
  endif 


  !===========================
  ! a more lengthy code
  !===========================
  if (method==2) then 
    fn = occ(isp,iband) / factor 
    en = eigen(isp,iband)

    ! instead of computing fm(1-fm)/fn(1-fm) 
    ! we compute log(fm(1-fm)) - log(fn(1-fn))
    !
    ! various situations  
    if (fn<tiny_value) then 
      log_fn_part = (fermi(isp)-en)/tsmear
    elseif (fn>tiny_value .and. fn<1-tiny_value) then
      log_fn_part = log(fn*(1.d0-fn))
    elseif (fn>1.d0-tiny_value) then 
      log_fn_part = (en-fermi(isp))/tsmear 
    endif 

    !
    ! d(fermi_level)/d(eigen_n) = 1 / [1 + sum_{m/=n} fm(1-fm)/fn(1-fn)]
    !
    zz = 0.0d0
    do ib=1,nband
      if (ib==iband) cycle

      fm = occ(isp,ib)/factor 
      em = eigen(isp,ib)

      ! different situatoins 
      if (fm<tiny_value) then 
        log_fm_part = (fermi(isp)-em)/tsmear
      elseif (fm>tiny_value .and. fm<1.d0-tiny_value) then 
        log_fm_part = log(fm*(1.d0-fm))
      elseif (fm>1.d0-tiny_value) then 
        log_fm_part = (em-fermi(isp))/tsmear
      endif

      expo = log_fm_part - log_fn_part 
      expo = limit_expo(expo)

      zz = zz + exp(expo)
    enddo 
    dmu_deigen = 1.d0/(1.d0+zz)

  endif 

  
  



  contains 

  function limit_expo(expo)
      implicit none 
      real(8) :: expo
      real(8) :: limit_expo

      if (expo>50.d0) then 
         limit_expo =  50.d0
      elseif (expo<-50.d0) then 
         limit_expo = -50.d0
      else 
         limit_expo = expo
      endif 

  endfunction 
  
end subroutine 
