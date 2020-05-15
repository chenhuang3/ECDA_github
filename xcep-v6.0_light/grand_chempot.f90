!
! chemical potential for grand canonical potential
!
subroutine grand_chempot(nsys,nspin,q_alpha,q_beta,eigen,nband, & 
 nconfig_max,nconfig,config,grand_temp,sub_etotal,weight,fermi)

 implicit none 

 integer :: nsys,nconfig(nsys), nspin,&
            nconfig_max, & 
            nband(nsys)
 
 real(8) :: weight(nconfig_max,nsys), fermi(nspin), & 
            eigen(nconfig_max,2,100,nsys), & 
            sub_etotal(nconfig_max,nsys), & 
            dsum, grand_temp, min_energy, & 
            config(2,nconfig_max,nsys)  ! number of spin up and spin down electrons

 double precision :: & 
   low, high, middle, dd, & 
   q_alpha, q_beta, workQ,& 
   min_eigen, max_eigen

 integer :: iconf, s, isp, ib



 print *,'enter chempot'

 ! get the lowest/highest eigenvalues 
 do isp=1,nspin
   do s = 1,nsys
     do ib=1,nband(s)
       if ( (isp==1) .and. (s==1) .and. (ib==1)) then 
         min_eigen = eigen(1,1,1,1)
         max_eigen = eigen(1,1,1,1)
       else if ( min_eigen > eigen(1,isp,ib,s) ) then 
         min_eigen = eigen(1,isp,ib,s)
       else if ( max_eigen < eigen(1,isp,ib,s) ) then 
         max_eigen = eigen(1,isp,ib,s)
       endif
     enddo
   enddo
 enddo

 print *,'chempot: nband: ',nband
 print *,'chempot: min_eigen: ',min_eigen
 print *,'chempot: max_eigen: ',max_eigen

 grand_temp = 0.0036749 ! hartree (=0.1 eV)
! grand_temp = 0.0036749*2.0d0 ! hartree (=0.2 eV)

 do isp=1,nspin
   low  =  min_eigen
   high =  max_eigen
   if (isp==1) workQ=q_alpha
   if (isp==2) workQ=q_beta
   !
   do while (.true.) 
     middle = (low+high)/2.d0
     dd=0.d0
     do s=1,nsys
       do ib=1,nband(s)
         dd = dd + 1.d0/(1.d0 + exp((eigen(1,isp,ib,s)-middle)/grand_temp))
       enddo
     enddo
     !print *,"dd=",dd,' middle: ',middle
     if (dd > workQ) then 
       high = middle 
     else
       low = middle 
     endif
     if ( (isp==1) .and. (abs(q_alpha-dd)<1e-8) ) then 
       fermi(1) = middle 
       write(6,'(a,f12.4,a)')'chempot: spin-up fermi level:   ',middle,' Hartree'
       exit 
     endif
     if ( (isp==2) .and. (abs(q_beta-dd)<1e-8) ) then 
       fermi(2) = middle 
       write(6,'(a,f12.4,a)')'chempot: spin-down fermi level: ',middle,' Hartree'
       exit 
     endif
   enddo
 enddo

 !
 ! get charges on each subsystem 
 !
 do s=1,nsys
   write(6,'(a,i3,a)',advance='no')'chempot: subsys[',s,'] q_alpha, q_beta, total_q: '
   do isp=1,nspin
     middle = (low+high)/2.d0
     dd=0.d0
     do ib=1,nband(s)
       dd = dd + 1.d0/(1.d0+exp((eigen(1,isp,ib,s)-fermi(isp))/grand_temp))
     enddo
     config(isp,1,s) = dd
     write(6,'(f16.6)',advance='no') dd
   enddo
   write(6,'(f16.6)',advance='no') sum(config(:,1,s))
   print *,''
 enddo


end subroutine grand_chempot
