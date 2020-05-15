!
! calculate kinetic energy and potential for give density rhor, as 
! well as its nonlocal psp energy and entropy energy
!
subroutine calculate_kedf(isys,n1,n2,n3,nspin,ike,rhor,qvec,dvol,ts_pot,ts_energy)

 implicit none

 integer      :: n1,n2,n3,nspin,ike,isys
 real(kind=8) :: rhor(n1*n2*n3,nspin), & 
                 qvec(3,(n1/2+1)*n2*n3), & 
                 dvol, & 
                 ts_pot(n1*n2*n3,nspin), & 
                 ts_energy

 ! local vars 
 
 character(len=50) :: str1,str2
 character(len=100) :: path
 integer      :: nfft, isp
 real(kind=8) :: tmp_energy, & 
                 tmp_rho(n1*n2*n3), & 
                 tmp_pot(n1*n2*n3), & 
                 pot_up(n1*n2*n3), & 
                 pot_dn(n1*n2*n3)

 ! function begins 

 ts_energy   = 0.0d0
 ts_pot      = 0.0d0
 nfft        = n1*n2*n3

 ! make density positive everywhere
 where (rhor<0.D0) 
   rhor = 1e-12
 endwhere

 !
 !  compute kinetic energy 
 !  for spin polarized case 
 !  T = 0.5*T[2*n_alpha] + 0.5*T[2*n_beta]
 !
 select case(ike)
 case(1) 
   do isp=1,nspin
     call TF(nfft,rhor(:,isp)*dble(nspin),dvol,tmp_pot,tmp_energy)
     ts_energy     = ts_energy + 1.d0/dble(nspin)*tmp_energy
     ts_pot(:,isp) = tmp_pot
   enddo

 case(2) 
   do isp=1,nspin
     call VW(nfft,n1,n2,n3,rhor(:,isp)*dble(nspin),qvec,dvol,tmp_pot,tmp_energy)
     tmp_energy = 1.d0/dble(nspin)*tmp_energy
     print *,'vw: ',tmp_energy
     ts_energy     = ts_energy + tmp_energy
     ts_pot(:,isp) = tmp_pot
   enddo

 case(3)
   do isp=1,nspin
     ! TF + 1/9 VW 
     call TF(nfft,rhor(:,isp)*dble(nspin),dvol,tmp_pot,tmp_energy)
     ts_energy     = ts_energy     + 1.0d0/dble(nspin)*tmp_energy
     ts_pot(:,isp) = ts_pot(:,isp) + tmp_pot
     call VW(nfft,n1,n2,n3,rhor(:,isp)*dble(nspin),qvec,dvol,tmp_pot,tmp_energy)
     ts_energy     = ts_energy     + 1.0d0/dble(nspin)*tmp_energy*1.0D0/9.0D0
     ts_pot(:,isp) = ts_pot(:,isp) + tmp_pot*1.0D0/9.0D0
   enddo

 case(4)
   do isp=1,nspin
     ! TF + VW + Huang-Carter nonlocal KEDF term
     call TF(nfft,rhor(:,isp)*dble(nspin),dvol,tmp_pot,tmp_energy)
     ts_energy     = ts_energy     + 1.d0/dble(nspin)*tmp_energy
     ts_pot(:,isp) = ts_pot(:,isp) + tmp_pot
     call VW(nfft,n1,n2,n3,rhor(:,isp)*dble(nspin),qvec,dvol,tmp_pot,tmp_energy)
     ts_energy     = ts_energy     + 1.d0/dble(nspin)*tmp_energy
     ts_pot(:,isp) = ts_pot(:,isp) + tmp_pot
     call huang_carter(nfft,n1,n2,n3,rhor(:,isp)*dble(nspin),qvec,dvol,tmp_pot,tmp_energy)
     ts_energy     = ts_energy     + 1.d0/dble(nspin)*tmp_energy
     ts_pot(:,isp) = ts_pot(:,isp) + tmp_pot
   enddo

 case(5)
   if (nspin==1) then 
     call LLP_KEDF_spin(nfft,n1,n2,n3,sum(rhor,2)/2.d0,sum(rhor,2)/2.d0,qvec,dvol,pot_up,pot_dn,ts_energy)
     ts_pot(:,1) = (pot_up + pot_dn)/2.D0
   endif
   if (nspin==2) then 
     call LLP_KEDF_spin(nfft,n1,n2,n3,rhor(:,1),rhor(:,2),qvec,dvol,pot_up,pot_dn,ts_energy)
     ts_pot(:,1) = pot_up
     ts_pot(:,2) = pot_dn
   endif

 case default
   print *,'ike=',ike,' is not implemented yet! error stop, compute_int_pot.f90'
   stop

 end select

 print *,''
 print *,'  ---------- calculate_kedf.f90 ---------'
 write(6,'(a,f12.6)')      '   total Q     : ',sum(rhor)*dvol
 write(6,'(a,f12.6)')      '   total Q_alpha: ',sum(rhor(:,1))*dvol
 write(6,'(a,f12.6)')      '   total Q_beta : ',sum(rhor(:,2))*dvol
 write(6,'(a,f14.8,a)')    '   ke          : ',ts_energy,' hartree'
 if (nspin==1) &
   write(6,'(a,2es12.4  )')'   ke potential: ',minval(ts_pot(:,1)),maxval(ts_pot(:,1))
 if (nspin==2) then 
   write(6,'(a,2es12.4,a)')'   ke potential: ',minval(ts_pot(:,1)),maxval(ts_pot(:,1)),' alpha'
   write(6,'(a,2es12.4,a)')'               : ',minval(ts_pot(:,2)),maxval(ts_pot(:,2)),' beta'
 endif
 print *,''

end subroutine
