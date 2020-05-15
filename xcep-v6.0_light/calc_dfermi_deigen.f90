!
! Compute d(fermi level)/d(eigen)
! 
! The formula for this is very simple, 
!
!  d(fermi_level) / d(eigen_m) = 1 / [1 + sum_{m/=n} fm(1-fm)/fn(1-fn)] 
!
! why do we need to write such a lengthy code to do it? 
! because, the occupation numbers can get very close to 1 and zero. 
! therefore making the denominator go to zero, therefore you get NaN.
!
! input/output
! ============
! iband:  the band to be computed
! isp: the spin channel to be computed
! nband: number of bands
! nspin: number of spin
!
! 12/27/2017 (Chen Huang)
!
subroutine  calc_dfermi_deigen(j_atom,mxband,nband,nfft,nspin,occ,eigenvalues,tsmear,fermi)

  implicit none 

  integer :: j_atom, mxband, & 
             nband(2),  &  ! 1: cluster, 2: env 
             nfft, nspin

  real(8) :: occ(mxband,nspin,2), & 
             tsmear, & 
             eigenvalues(mxband,nspin,2), & 
             fermi(nspin), & 
             dfermi_deigen(mxband,nspin,2)  ! last index for cluster/env

  ! ABINIT vars 
  integer :: iband,isp
  real(8) :: ab_eigen(nspin,sum(nband)), & 
             ab_occ(nspin,sum(nband)), & 
             dmu_deigen

  ! local vars 
  character(len=200) :: fname, ss
  integer :: ib 
  real(8) :: zz, fn, fm, en, em, expo
  real(8) :: tiny_value=1e-6,  & 
             log_fn_part, & 
             log_fm_part

  
  ! >>>>>>>>>>> function begins <<<<<<<<<<<<!

  ! assign XCEP variables to ABINIT
  ! this save my effort to debug the program 

  do isp=1,nspin 
    ab_eigen(isp,1:nband(1)) = eigenvalues(1:nband(1),isp,1) ! cluster 
    ab_occ  (isp,1:nband(1)) = occ(1:nband(1),isp,1)         ! cluster 

    ab_eigen(isp,nband(1)+1:nband(1)+nband(2)) = eigenvalues(1:nband(2),isp,2) ! env
    ab_occ  (isp,nband(1)+1:nband(1)+nband(2)) = occ(1:nband(2),isp,2)         ! env
  enddo

  !
  ! compute d(Fermi level) / d(eigenvalue) for all 
  ! the bands in cluster and environment 
  !
  dfermi_deigen = 0.d0 

  do isp=1,nspin 
    do iband = 1,sum(nband)

      fn = ab_occ(isp,iband)
      en = ab_eigen(isp,iband)

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
    !  print *,'log_fn_part:' ,log_fn_part
  
      !
      ! d(fermi_level)/d(eigen_{n}) = 1 / [1 + sum_{m/=n} fm(1-fm)/fn(1-fn)]
      !
      zz = 0.0d0
      do ib=1,sum(nband)
        if (ib==iband) cycle

        fm = ab_occ(isp,ib)
        em = ab_eigen(isp,ib)

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

     !    print *,'log_fm_part:' ,log_fm_part
     !    print *,'expo: ',expo

        zz = zz + exp(expo)
      enddo 

      dmu_deigen = 1.d0/(1.d0+zz)
  
      if (iband>=nband(1)+1) then 
        dfermi_deigen(iband-nband(1),isp,2) = dmu_deigen ! cluster 
      else 
        dfermi_deigen(iband,isp,1) = dmu_deigen ! env 
      endif 
  
    enddo ! end of band 
  enddo ! end of spin 


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! write dfermi_deigen to cluster and env 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  write(ss,'(i4)')j_atom
  fname = 'subsys'//trim(adjustl(ss))//'/dfermi_deigen.dat'
  open(file=fname,unit=111,action='write',form='formatted')
  write(111,'(2f12.6)') fermi(1:2)
  do iband = 1,nband(1)
     write(111,'(2es16.8)')dfermi_deigen(iband,1:2,1)
  enddo 
  close(111)

  fname = 'subsys'//trim(adjustl(ss))//'_env/dfermi_deigen.dat'
  open(file=fname,unit=111,action='write',form='formatted')
  write(111,'(2f12.6)')fermi(1:2)
  do iband = 1,nband(1)
     write(111,'(2es16.8)')dfermi_deigen(iband,1:2,2)
  enddo 
  close(111)
  
!  print *,'dmu_deigen: ',dmu_deigen
!  stop




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


end subroutine calc_dfermi_deigen
