!
! this is spin polarized KEDF from Lee-Lee-Parr (PRA 1991)
! if the system is non-spin-polarized, then 
! you need to make rho_up   = rho/2
!                  rho_down = rho/2
!                    
! then the dT/drho = 1/2(pot_up + pot_dn)
!
!   Chen Huang
!   Feb/28/2011
!
subroutine  LLP_KEDF_spin(nfft,n1,n2,n3,rho_up,rho_dn,qvec,dvol,pot_up,pot_dn,ke)

 implicit none 
 
 ! external vars
 integer,intent(in) :: nfft,n1,n2,n3
 real(kind=8),intent(in ) :: rho_up(nfft),rho_dn(nfft),dvol
 real(kind=8),intent(in)  :: qvec(3,(n1/2+1)*n2*n3)
 real(kind=8),intent(out) :: pot_up(nfft),pot_dn(nfft),ke

 ! local vars
 integer :: ii,isp

 real(kind=8) , parameter :: & 
    gamma=0.0253d0, &   ! parameter in LLP KEDF
!    alpha=0.d0  , &  ! parameter in LLP KEDF
    alpha=4.4188e-3, &  ! parameter in LLP KEDF
    frt = 4.d0/3.d0, & 
    fvt = 5.d0/3.d0, & 
    tt = 2.d0/3.d0 , & 
    ot = 1.d0/3.d0, &
    cF = 2.8712d0

 real(kind=8) :: g(3,nfft), &                 ! gradient
           Gx(nfft),   &                ! the G(x) in paper
           dGx_dx(nfft), & 
           x(nfft), & 
           pot_tmp(nfft), & 
           rho_tmp(nfft), &
           norm_grho(nfft), &
           temp(nfft), &
           q3d(3,n1/2+1,n2,n3)

 complex(kind=8) :: & 
   ffttmp(n1/2+1,n2,n3)


 !===================
 ! functionl begins
 !===================
 print *,''
 print *,'enter  LLP()'
 print *,'rho_up=',minval(rho_up),maxval(rho_up)
 print *,'rho_dn=',minval(rho_dn),maxval(rho_dn)

 q3d = reshape(qvec,(/3,n1/2+1,n2,n3/))

 !compute potential
 !===================
 ! for spin up
 !
 ke = 0.d0
 do isp=1, 2
   if (isp==1)  rho_tmp = rho_up
   if (isp==2)  rho_tmp = rho_dn
   call gradient (nfft,n1,n2,n3,rho_tmp,qvec,g)

   ! compute x
   do ii=1,nfft
     norm_grho(ii) = sqrt(g(1,ii)**2+g(2,ii)**2+g(3,ii)**2)
     if (rho_tmp(ii)>1e-4) then 
       x(ii) = norm_grho(ii) / rho_tmp(ii)**frt
     else 
       x(ii) = norm_grho(ii) / 1e-4**frt
     endif
   enddo

   print *,'norm_g_rho =',minval(norm_grho),maxval(norm_grho),'<- spin:',isp

   ! compute Gx
!   Gx = x**2 / (1.d0 + gamma * asinh(x))

   ! compute dGx / dx
!   dGx_dx = 2.d0*x / (1.d0+gamma*asinh(x)) - & 
!            x**2* gamma/sqrt(1.d0+x**2)/(1.d0+gamma*asinh(x))**2
   
   ! pot A
   pot_tmp = 2.d0**tt * CF * fvt * rho_tmp**tt * (1.d0+alpha*Gx) ! pot A
   print *,'potA =',minval(pot_tmp),maxval(pot_tmp),'<- spin:',isp

   ! pot B
   pot_tmp = pot_tmp + 2.d0**tt * CF * alpha * dGx_dx * (-4.0d0/3.d0) * x * rho_tmp**tt
   print *,'+potB=',minval(pot_tmp),maxval(pot_tmp),'<- spin:',isp
   
   ! pot C
   do ii=1,3
     where (rho_tmp>1e-4) 
!       temp = alpha*(  2.d0/(1.d0+gamma*asinh(x)) & 
!                   - x*gamma/(1.d0+gamma*asinh(x))**2/sqrt(1.d0+x**2)) / rho_tmp
     elsewhere
!       temp = alpha*(  2.d0/(1.d0+gamma*asinh(x)) & 
!                   - x*gamma/(1.d0+gamma*asinh(x))**2/sqrt(1.d0+x**2)) / 1e-4
     endwhere
     temp = temp * g(ii,:)
     call fft(n1,n2,n3,temp,ffttmp,1)
     ffttmp = dcmplx(0.d0,-q3d(ii,:,:,:))*ffttmp
     call fft(n1,n2,n3,temp,ffttmp,-1)
     pot_tmp = pot_tmp - 2.d0**tt*CF*temp
     print *,'+potC=',minval(pot_tmp),maxval(pot_tmp),'<- spin:',isp
   enddo 

   if (isp==1) pot_up = pot_tmp
   if (isp==2) pot_dn = pot_tmp
   
   ! kinetic energy
   ke = ke + 2.d0**tt * CF * dvol * sum(rho_tmp**fvt*(1.d0+alpha*Gx)) 

 enddo

 print *,'LLP: ke=', ke
 print *,'LLP: pot_up =',minval(pot_up),maxval(pot_up)
 print *,'LLP: pot_dn =',minval(pot_dn),maxval(pot_dn)

 print *,''
 print *,'leave LLP()'
 print *,''

! stop

end subroutine LLP_KEDF_spin
