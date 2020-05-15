!
! code from PROFESS, compute the xc of Perdew and Wang 1992
!
SUBROUTINE LSDAPW92(n1,n2,n3,rho,LDAPotential,LDAEnergy,dVol)

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  integer,parameter :: dp=8
  integer,intent(in)        :: n1,n2,n3 
  real(kind=dp),intent(in)  :: rho(n1,n2,n3,2),dVol
  real(kind=dp),intent(out) :: LDAPotential(n1,n2,n3,2)
  real(kind=dp),intent(out) :: LDAEnergy

                     !>> INTERNAL VARIABLES <<!

  real(kind=dp), parameter :: &
    pi = 3.1415926d0, &
    p75vpi = 0.75d0/pi, &
    ax = -0.7385588d0, &
    fzz = 1.709921d0, &        
    gamma = 0.5198421d0, &
    one =1.d0, &
    two = 2.d0, &
    three = 3.d0, &
    four = 4.d0

  REAL(kind=DP) :: &
    rs, &                              ! (3/(4.pi.rho))^1/3
    zet, &                             ! Spin-polarization
    ex(2), vx(2),exc, &                           ! exchange energy
    ec, vc(2), &                              ! correlation energy
    eu,f,z4,ecrs,eurs,eczet,comm,fz,d, &  ! work variables
    ep,eprs,alfrsm,alfm,ac2,third         ! work variables

  INTEGER :: &
    ix, iy, iz                           ! Dummy counters

  !! we ignore spin now
  !! 
  write(6,'(a)')' enter LSDAPW92()'

  third = one/three
  ac2 = four/three
  exc= 0.d0
  DO iz=1, n3
     DO iy=1, n2
        DO ix=1, n1
           ! exchange
           d=two*rho(ix,iy,iz,1)
           ex(1)= ax*d**third
           vx(1)=ex(1)*ac2
           d=two*rho(ix,iy,iz,2)
           ex(2)= ax*d**third
           vx(2)=ex(2)*ac2
           ! local correlation
           d=sum(rho(ix,iy,iz,:))
           rs=(p75vpi/d)**third
           call spn_gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0, &
                1.6382d0,0.49294d0,1.00d0,rs,eu,eurs)
           call spn_gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0, &
                3.3662d0,0.62517d0,1.00D0,rs,ep,eprs)
           call spn_gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0, &
                0.88026d0,0.49671d0,1.00d0,rs,alfm,alfrsm)
           zet=(rho(ix,iy,iz,1)-rho(ix,iy,iz,2))/d
           f = ((one+zet)**ac2+(one-zet)**ac2-two)/gamma
           z4 = zet**4
           ec = eu*(one-f*z4)+ep*f*z4-alfm*f*(one-z4)/fzz
           ecrs = eurs*(one-f*z4)+eprs*f*z4-alfrsm*f*(one-z4)/fzz
           fz = ac2*((one+zet)**third-(one-zet)**third)/gamma
           eczet = four*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu &
                -(one-z4)*alfm/fzz)
           comm = ec -rs*ecrs/three-zet*eczet
           vc(1) = comm + eczet
           vc(2) = comm - eczet

!           vxc(i,1)=two*(vx(1)+vc(1))
!           vxc(i,2)=two*(vx(2)+vc(2))
           ec=(rho(ix,iy,iz,1)*ex(1)+rho(ix,iy,iz,2)*ex(2)+ec*d)
           exc = exc + ec
           LDAPotential(ix,iy,iz,1)=(vx(1)+vc(1)) !- 0.5d0 * emf
           LDAPotential(ix,iy,iz,2)=(vx(2)+vc(2)) !+ 0.5d0 * emf
        ENDDO
     ENDDO
  ENDDO
  LDAEnergy = dVol*exc

  write(6,'(a,f16.10,a)')' leave LSDAPW92(), e_lda=',LDAEnergy,' hartree'
END SUBROUTINE LSDAPW92

!---------------------------------------------------------------
!
! This subroutine computes the local correlation energy and
! potential for the Perdew-Wang exchange-correlation scheme.
!
!---------------------------------------------------------------
SUBROUTINE spn_gcor(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs)

  implicit none
  !
  ! Input/Output variables:
  !
  integer,parameter :: dp=8
  real(dp), intent(in) :: a,a1,b1,b2,b3,b4,p,rs
  real(dp), intent(out) :: gg,ggrs
  !
  ! Work variables:
  !
  real(dp) :: p1,q0,rsp,q1,q2,q3,rs12,rs32,two, one,three
  !---------------------------------------------------------------
  one = 1.d0
  two = 2.d0
  three = 3.d0
  p1 = p + 1.d0
  q0 = -two*a*(one+a1*rs)
  rs12 = sqrt(rs)
  rs32 = rs12**3
  rsp = rs**p
  q1 = two*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
  q2 = log(one+one/q1)
  gg = q0*q2
  q3 = a*(b1/rs12+two*b2+three*b3*rs12+two*b4*p1*rsp)
  ggrs = -two*a*a1*q2-q0*q3/(q1**2+q1)

END SUBROUTINE spn_gcor

