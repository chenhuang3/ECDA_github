
! ================================================
!  this subroutine is adapted from PROFESS code 
! ================================================

!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is a general version of PBE using libxc and work for both spin
!   unpolarized and polarized cases. The algo is the same as ABINIT libxc PBE
!   subroutine.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!      libxc and algorithm in ABINIT rhohxc,drivexc,libxc_functional, xcpot, xcmult
!------------------------------------------------------------------------------
! REVISION LOG:
!  02/28/2014 : Created By Steven Xia
!  03/04/2014 : Ported to new FFT, general optimizations (JMD)
!------------------------------------------------------------------------------

subroutine calculate_x_pbe_eps(n1,n2,n3,nspin,rho,qvec,dvol,xc_pot,xc_eps)

     use xc_f90_types_m
     use xc_f90_lib_m
     use libxc_funcs_m

     implicit none

     integer      :: n1,n2,n3,nspin
     real(kind=8) :: dvol
     real(kind=8),intent(in) :: rho(n1,n2,n3,nspin)     ! you can just pass rho(n1*n2*n3,nspin)
     real(kind=8) :: qvec(3,(n1/2+1)*n2*n3)
     real(kind=8) :: xc_pot(n1,n2,n3,nspin)  ! you can just pass xc_pot(n1*n2*n3,nspin) 
     real(kind=8) :: xc_eps(n1,n2,n3)

     !working variables
     !gradup2 = gradup * gradup
     !sigma(:,:,:,1) = gradup*gradup
     !sigma(2)=gradup*graddn
     !simga(3)=graddn*graddn
     !vsigma = dE/dsigma
     !depsxc(:,:,:,1) = dE/dup
     !depsxc(:,:,:,2) = dE/ddn
     !depsxc(:,:,:,3) = 1/|gradup|*dE/d|gradup|
     !depsxc(:,:,:,4) = 1/|graddn|*dE/d|graddn|
     !depsxc(:,:,:,5) = 1/|gradtot|*dE/d|gradtot|
     !rhonow(:,:,:,1,1) = up
     ! if numspin = 1
     !rhonow(:,:,:,1,2) = gradtot/|gradtot|*dE/d|gradtot|
     ! if numspin = 2
     !rhonow(:,:,:,1,2) = (gradup/|gradup|*dE/d|gradup| + gradtot/|gradtot| * dE/d|gradtot|)xdir
     !rhonow(:,:,:,1,3) = (gradup/|gradup|*dE/d|gradup| + gradtot/|gradtot| * dE/d|gradtot|)ydir
     !rhonow(:,:,:,1,4) = (gradup/|gradup|*dE/d|gradup| + gradtot/|gradtot| * dE/d|gradtot|)zdir
     !rhonow(:,:,:,2,:) is for spin_dn and the fifth dimension is the same as described above

     !rhotmp,vctmp,vxtmp,sigmatmp,vxsigmatmp,vcsigmatmp,extmp,ectmp are temp variables for libxc calling

     ! local variables (-> please at some point move them to heap [JMD])
     COMPLEX(kind=8), DIMENSION(n1/2+1,n2,n3) :: rhoRecip,trafo1,trafo2,trafo3
     real(kind=8) :: gtmp(3,n1*n2*n3)
     real(kind=8) :: q3d(3,n1/2+1,n2,n3)
     real(kind=8),dimension(n1,n2,n3):: rhotot,up,dn,excpbe
     real(kind=8),dimension(n1,n2,n3,3):: gradup,graddn,gradtot,sigma
     real(kind=8),dimension(n1,n2,n3,5):: depsxc
     real(kind=8),dimension(n1,n2,n3,nspin,4):: rhonow
     real(kind=8):: rhotmp(nspin),vctmp(nspin),vxtmp(nspin)
     real(kind=8):: sigmatmp(3),vxsigmatmp(3),vcsigmatmp(3)
     real(kind=8):: extmp,ectmp

     integer:: i,j,k,ispin

     TYPE(xc_f90_pointer_t) :: x_func, c_func
     TYPE(xc_f90_pointer_t) :: x_info, c_info

     real(8) :: correlation(n1,n2,n3), exchange(n1,n2,n3)

!======================================function body=====================================

  !   write(*,*) "USING NEW SPIN LIBXC!!!!"

     !first get parameter needed for libxc
     If(nspin==1) then
       rhotot = rho(:,:,:,1)
     Else
       up = rho(:,:,:,1)
       dn = rho(:,:,:,2)
       rhotot = up + dn
     Endif

     q3d(1,:,:,:) = reshape(qvec(1,:),(/n1/2+1,n2,n3/))
     q3d(2,:,:,:) = reshape(qvec(2,:),(/n1/2+1,n2,n3/))
     q3d(3,:,:,:) = reshape(qvec(3,:),(/n1/2+1,n2,n3/))

     !compute gradtot
     call gradient(n1*n2*n3,n1,n2,n3,rhotot,qvec,gtmp)
     gradtot(:,:,:,1) = reshape(gtmp(1,:),(/n1,n2,n3/))
     gradtot(:,:,:,2) = reshape(gtmp(2,:),(/n1,n2,n3/))
     gradtot(:,:,:,3) = reshape(gtmp(3,:),(/n1,n2,n3/))

     If(nspin==1) then
         sigma(:,:,:,1) = gradtot(:,:,:,1)*gradtot(:,:,:,1) &
                         +gradtot(:,:,:,2)*gradtot(:,:,:,2) &
                         +gradtot(:,:,:,3)*gradtot(:,:,:,3)
         sigma(:,:,:,2) = 0.d0
         sigma(:,:,:,3) = 0.d0
     Else
         ! compute gradup, gradup2
         call gradient(n1*n2*n3,n1,n2,n3,up,qvec,gtmp)
         gradup(:,:,:,1) = reshape(gtmp(1,:),(/n1,n2,n3/))
         gradup(:,:,:,2) = reshape(gtmp(2,:),(/n1,n2,n3/))
         gradup(:,:,:,3) = reshape(gtmp(3,:),(/n1,n2,n3/))

         sigma(:,:,:,1) = gradup(:,:,:,1)*gradup(:,:,:,1) &
                         +gradup(:,:,:,2)*gradup(:,:,:,2) &
                         +gradup(:,:,:,3)*gradup(:,:,:,3)

         ! compute graddn, graddn2
         call gradient(n1*n2*n3,n1,n2,n3,dn,qvec,gtmp)
         graddn(:,:,:,1) = reshape(gtmp(1,:),(/n1,n2,n3/))
         graddn(:,:,:,2) = reshape(gtmp(2,:),(/n1,n2,n3/))
         graddn(:,:,:,3) = reshape(gtmp(3,:),(/n1,n2,n3/))

         sigma(:,:,:,3) = graddn(:,:,:,1)*graddn(:,:,:,1) &
                         +graddn(:,:,:,2)*graddn(:,:,:,2) &
                         +graddn(:,:,:,3)*graddn(:,:,:,3)

         !compute gradup*graddn
         sigma(:,:,:,2) = gradup(:,:,:,1)*graddn(:,:,:,1) &
                         +gradup(:,:,:,2)*graddn(:,:,:,2) &
                         +gradup(:,:,:,3)*graddn(:,:,:,3)
     Endif

     !use libxc to get dE/drho_up/dn and dE/dsigma
     call xc_f90_func_init(x_func, x_info, XC_GGA_X_PBE, nspin)
     call xc_f90_func_init(c_func, c_info, XC_GGA_C_PBE, nspin)

     exchange    = 0.d0 
     correlation = 0.0d0 

     ! JMD: it may be possible to get more speed out of this by reordering the loops.
     do k = 1,n3
         do j = 1,n2
             do i = 1,n1
                 rhotmp = rho(i,j,k,:)
                 sigmatmp = sigma(i,j,k,:)
                 call xc_f90_gga_exc_vxc(x_func,1,rhotmp(1),sigmatmp(1),extmp,vxtmp(1),vxsigmatmp(1))
                 call xc_f90_gga_exc_vxc(c_func,1,rhotmp(1),sigmatmp(1),ectmp,vctmp(1),vcsigmatmp(1))

                 !excpbe(i,j,k)   =  extmp + ectmp
                 excpbe(i,j,k)   =  extmp
                 !depsxc(i,j,k,1) = vxtmp(1) + vctmp(1)
                 depsxc(i,j,k,1) =  vxtmp(1)

                 If(nspin==1) then
                     !depsxc(i,j,k,5) =  2.d0 * (vcsigmatmp(1) + vxsigmatmp(1))
                     depsxc(i,j,k,5) =  2.d0 * vxsigmatmp(1)
                 Else
                     !depsxc(i,j,k,2) =  vxtmp(2) + vctmp(2)
                     !depsxc(i,j,k,3) =  2.d0*(vxsigmatmp(1)+vcsigmatmp(1)) - vxsigmatmp(2) - vcsigmatmp(2)
                     !depsxc(i,j,k,4) =  2.d0*(vxsigmatmp(3)+vcsigmatmp(3)) - vxsigmatmp(2) - vcsigmatmp(2)
                     !depsxc(i,j,k,5) =  vxsigmatmp(2) + vcsigmatmp(2)

                     depsxc(i,j,k,2) =  vxtmp(2)
                     depsxc(i,j,k,3) =  2.d0*(vxsigmatmp(1)) - vxsigmatmp(2)
                     depsxc(i,j,k,4) =  2.d0*(vxsigmatmp(3)) - vxsigmatmp(2)
                     depsxc(i,j,k,5) =  vxsigmatmp(2)
                 Endif
             enddo
         enddo
     enddo

     call xc_f90_func_end(x_func)
     call xc_f90_func_end(c_func)

     ! now already get depsxc, calculate rhonow
     If(nspin==1) then
         rhonow(:,:,:,1,1) = rhotot
         rhonow(:,:,:,1,2) = depsxc(:,:,:,5) * gradtot(:,:,:,1)
         rhonow(:,:,:,1,3) = depsxc(:,:,:,5) * gradtot(:,:,:,2)
         rhonow(:,:,:,1,4) = depsxc(:,:,:,5) * gradtot(:,:,:,3)
     Else
         rhonow(:,:,:,1,1) = up(:,:,:)
         rhonow(:,:,:,2,1) = dn(:,:,:)

         rhonow(:,:,:,1,2) = gradup(:,:,:,1)*depsxc(:,:,:,3)+gradtot(:,:,:,1)*depsxc(:,:,:,5)
         rhonow(:,:,:,1,3) = gradup(:,:,:,2)*depsxc(:,:,:,3)+gradtot(:,:,:,2)*depsxc(:,:,:,5)
         rhonow(:,:,:,1,4) = gradup(:,:,:,3)*depsxc(:,:,:,3)+gradtot(:,:,:,3)*depsxc(:,:,:,5)

         rhonow(:,:,:,2,2) = graddn(:,:,:,1)*depsxc(:,:,:,4)+gradtot(:,:,:,1)*depsxc(:,:,:,5)
         rhonow(:,:,:,2,3) = graddn(:,:,:,2)*depsxc(:,:,:,4)+gradtot(:,:,:,2)*depsxc(:,:,:,5)
         rhonow(:,:,:,2,4) = graddn(:,:,:,3)*depsxc(:,:,:,4)+gradtot(:,:,:,3)*depsxc(:,:,:,5)
     Endif

     ! now already get rhonow, do divergence calculation
     do ispin =1,nspin

         call fft(n1,n2,n3,rhonow(:,:,:,ispin,2),trafo1,1)
         call fft(n1,n2,n3,rhonow(:,:,:,ispin,3),trafo2,1)
         call fft(n1,n2,n3,rhonow(:,:,:,ispin,4),trafo3,1)

         trafo1 = DCMPLX(0.d0,q3d(1,:,:,:))*trafo1 + &
                  DCMPLX(0.d0,q3d(2,:,:,:))*trafo2 + & 
                  DCMPLX(0.d0,q3d(3,:,:,:))*trafo3 
 
         call fft(n1,n2,n3,xc_pot(:,:,:,ispin),trafo1,-1)
         xc_pot(:,:,:,ispin) = - xc_pot(:,:,:,ispin) + depsxc(:,:,:,ispin) 

     enddo

     !now energy
     xc_eps = rhotot*excpbe

!     print *,'pbe: exchange energy: ',sum(rhotot*exchange)*dvol 
!     print *,'pbe: correlation energy: ',sum(rhotot*correlation)*dvol
!
!     print *,''
!     print *,               '   -------- calculate_xc_pbe.f90 --------'
!     write(6,'(a,f14.6)')   '       xc_energy: ',xc_energy
!     if (nspin==1) then 
!       write(6,'(a,2f16.6)')'       xc_pot: ', & 
!       minval(xc_pot(:,:,:,1)),maxval(xc_pot(:,:,:,1))
!     endif
!     if (nspin==2) then 
!       write(6,'(a,2f16.6)')'       xc_pot[alpha]: ', & 
!       minval(xc_pot(:,:,:,1)),maxval(xc_pot(:,:,:,1))
!       write(6,'(a,2f16.6)')'       xc_pot[beta ]: ', & 
!       minval(xc_pot(:,:,:,2)),maxval(xc_pot(:,:,:,2))
!     endif

end subroutine 
