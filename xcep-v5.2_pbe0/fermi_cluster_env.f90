!
!
!  based on the eigenvalues from ALL subsystem, 
!  recompute the occupation numbers
!  new occupation number will be written for subsystems via the file new_occ.tmp
!
!
subroutine fermi_cluster_env(mxband,file_dfet_out,nspin,nelectr, & 
                             tsmear,nband,eigenvalue,occ_new,fermi_level)

 implicit none

 ! we have two subsystem here: cluster + env 
 integer,parameter :: nsys = 2 
 integer   :: mxband,nspin,file_dfet_out,nband( nsys )
 real(8)   :: eigenvalue(mxband,nspin,nsys), &  ! the last index indicates cluster or env
              occ_new(mxband,nspin,nsys),    & 
              nelectr(nspin), & 
              tsmear

 ! local vars ............

 character(len=500) :: fname,ss,sc,stemp
 integer :: iband,isp,is,system,counter,total_nband,jj,itmp,nkpt(nsys)
 real(8) :: th1,th2,th0,kbT,occ,totalQ, & 
            entropy_tmp(nspin),fermi_level(nspin),aaa,ent,occ_tmp,ftmp,ab_tsmear


 ! function .....................

 print *,'enter fermi_cluster_env().'

 kbT = tsmear 
 print *,'kbT:   ',kbT,' hartree'
 print *,'totalQ:',nelectr(1:nspin)
 print *,'bisecing to compute the fermi levels ...'
 call flush(6)

 !
 ! find the fermi level by bisecting 
 ! 
 do isp=1,nspin

   th0 = min(minval(eigenvalue(1:nband(1),isp,1)),minval(eigenvalue(1:nband(2),isp,2))) - 10.0
   th2 = max(maxval(eigenvalue(1:nband(1),isp,1)),maxval(eigenvalue(1:nband(2),isp,2))) + 10.0

   print *,"lower bound: ",th0," upper bound: ",th2

   th1 = (th0+th2)/2.d0

   ! bisecting .........
   do while (.true.)
     totalQ = 0.d0
     !
     ! loop over cluster and env 
     !
     do system=1,nsys
       ! 
       do iband=1,nband(system)
         aaa = (eigenvalue(iband,isp,system)-th1)/kBT
         if (aaa>40.d0) then 
           occ = 0.d0 
         else if (aaa<-40.d0) then 
           if (nspin==2) occ = 1.d0 
           if (nspin==1) occ = 2.d0 
         else
           if (nspin==2) occ = 1.d0 / (exp(aaa) + 1.d0) 
           if (nspin==1) occ = 2.d0 / (exp(aaa) + 1.d0) 
         endif
         occ_new(iband,isp,system) = occ
         totalQ = occ + totalQ
       enddo
     enddo
     if (totalQ < nelectr(isp)) then 
       th0 = th1
     else
       th2 = th1
     endif
     th1 = (th0+th2)/2.d0
   !!  print *,'totalQ,  target ',totalQ,nelectr(isp)
   !!  if ( abs(totalQ-nelectr(isp))<1e-12 ) exit
     if ( abs(th0-th2)<1e-6 .and. abs(totalQ-nelectr(isp))<1e-12 ) exit
   enddo 
!   write(file_dfet_out,'(a,f12.6,a,i2,a)') & 
!   '[fermi_level_across_subsystem] fine new fermi level: ',th1,' (spin: ',isp,')'
   print *,'fermi level:',th1
   fermi_level(isp) = th1
 enddo 

! print *,'cluster eigen:'
! print *,'spin 1:',eigenvalue(1:nband(1),1,1)
! print *,'spin 2:',eigenvalue(1:nband(1),2,1)
!
! print *,'env eigen:'
! print *,'spin 1:',eigenvalue(1:nband(1),1,2)
! print *,'spin 2:',eigenvalue(1:nband(1),2,2)
!
! print *,'cluster occ_new:'
! print *,'spin 1:',occ_new(1:nband(1),1,1)
! print *,'spin 2:',occ_new(1:nband(1),2,1)
! print *,'env occ_new:'
! print *,'spin 1:',occ_new(1:nband(1),1,2)
! print *,'spin 2:',occ_new(1:nband(1),2,2)
!
! print *,'done bisecing fermi levels.'
! print *,'left fermi_cluster_env().'

end subroutine
