!
! compute the entropy based on occ 
!
subroutine  entropy(mxband, nspin, occ, ss)

  implicit none 

  integer :: mxband, nspin, isp, ib
  real(8) :: occ(mxband,nspin), tsmear, ss, docc, factor 

  ss = 0.d0 
  if (nspin==2) factor = 1.d0 
  if (nspin==1) factor = 2.d0

  do isp = 1,nspin 
      do ib = 1,mxband 
         docc = occ(ib,isp)/factor 
         if ( docc>1e-6 .and. docc<1.d0-1e-6) then 
           ss = ss + docc*log(docc) + (1.d0-docc)*log(1.d0-docc)
         endif 
      enddo 
  enddo 

  ss = -ss
  if (nspin==1) ss = ss*2.d0 

end subroutine
