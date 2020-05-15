!=========================================
! compute VW kinetic energy and potential
! with incoming electron density
!
!  Chen Huang  Jan/2011 
!  
!=========================================
subroutine VW(nfft,n1,n2,n3,rhor,qvec,dvol,vw_reg,vw_pot,vw_energy)
 implicit none

 integer,intent(in) :: nfft,n1,n2,n3
 real(kind=8),intent(in)  :: rhor(nfft), & 
                             vw_reg, & 
                             qvec(3,(n1/2+1)*n2*n3),dvol
 real(kind=8),intent(out) :: vw_pot(nfft),vw_energy

 ! local vars
 integer :: i
 real(kind=8) :: gnorm,g(3,nfft),lap(nfft),sq_rhor(nfft)


 !=========================================================================
 ! NOTE:
 !=========================================================================
 ! When electron density is low, VW potential can be large for small density
 ! VW functional is regularized by defining 
 !
 !  phi = sqrt(rho+vw_reg)
 !
 !  VW = phi*(-1/2*\nabla^2)*phi
 !
 ! I set vw_reg = 1e-5 in general. If we set vw_reg = 0, 
 ! we restore the original VW energy functional 
 !
 ! The potential is then 
 !
 ! VW_pot = \delta VW/delta rho = - \nabla^2 phi / (2*phi)
 !
 !=========================================================================

 sq_rhor = sqrt(rhor+vw_reg)   ! regularized 

!! call gradient (nfft,n1,n2,n3,sq_rhor,qvec,g)
 call laplacian(nfft,n1,n2,n3,sq_rhor,qvec,lap)

 vw_energy = sum(sq_rhor*(-0.5d0)*lap)*dvol 
 vw_pot = - 0.5d0 * lap / sq_rhor
 
! write(6,'(a,2es12.4)')'min/max rho    = ',minval(rhor),maxval(rhor)
! write(6,'(a,2es12.4)')'min/max qvec   = ',minval(qvec),maxval(qvec)
! write(6,'(a,2es12.4)')'min/max sqrho  = ',minval(sq_rhor),maxval(sq_rhor)
! write(6,'(a,2es12.4)')'min/max g      = ',minval(g),maxval(g)
! write(6,'(a,2es12.4)')'min/max lap    = ',minval(lap),maxval(lap)
! write(6,'(a,es20.12)')'vw_energy      = ',vw_energy
! write(6,'(a,2es12.4)')'min/max vw_pot = ',minval(vw_pot),maxval(vw_pot)
! print *,'leave vW() '
 return 
end subroutine VW
