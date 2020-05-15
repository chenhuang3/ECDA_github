!
! for given eigenvalues and electron number (q), compute the 
! occupation numbers according to the F-D smearing 
!
subroutine chen_fermi_occ(tsmear, nband, ei, q, occ, fermi)
 
   implicit none 

   integer :: nband,iband,iter
   real(8) :: tsmear, q, ei(nband), occ(nband), totalQ
   real(8) :: th0, th1, th2, aaa, fermi, occ_old(nband), max_diff

   !
   ! find the fermi level by bisecting 
   ! 
   th0 = minval(ei)-10.0;
   th2 = maxval(ei)+10.0;

   print *,'enter chen_fermi_occ() ...'
   print *,"   lower bound: ",th0," upper bound: ",th2
   print *,'   tsmear:',tsmear,' q:',q,' nband: ',nband

   th1 = (th0+th2)/2.d0
   iter = 0

   ! bisecting .........
   do while (.true.)

     iter = iter + 1

     totalQ = 0.d0

     do iband=1,nband
         aaa = (ei(iband) - th1)/tsmear
         if (aaa>40.d0) then 
           occ(iband) = 0.d0 
         else if (aaa<-40.d0) then 
           occ(iband) = 1.d0 
         else
           occ(iband) = 1.d0 / (exp(aaa) + 1.d0) 
         endif
         !print *,'occ(iband): ',occ(iband),' aaa: ',aaa
         totalQ = occ(iband) + totalQ
     enddo

     if (iter>3) then 
       max_diff = maxval( abs(occ_old-occ) )
     endif
     occ_old = occ 

     if (totalQ < q) then 
       th0 = th1
     else
       th2 = th1
     endif
     th1 = (th0+th2)/2.d0

     if ( abs(th0-th2)<1e-10 .and. abs(totalQ-q)<1e-12 .and. max_diff<1e-6 ) exit

!     print *,'th0, th2: ',th0,th2, ' abs(totalQ-q):',abs(totalQ-q),' abs(th0-th2):',abs(th0-th2),' totalQ:',totalQ

!!!!!!!!!!!!! debug 
!     iter = iter +1 
!     if (iter>100) then 
!         stop
!     endif

   enddo 

   fermi = th1

   print *,'done chen_fermi_occ(). fermi level: ',th1

end subroutine 
