!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function obtained from ABINIT version 8.10.2
! The bare_vqg() function for ABINIT version 7 does not compute error function and Spencer-Alavi method 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!****f* ABINIT/bare_vqg
!! NAME
!! bare_vqg
!!
!! FUNCTION
!! Compute bare coulomb term in G-space on the FFT mesh i.e. 4pi/(G+q)**2
!!
!! INPUTS
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  gsqcut=cutoff value on G**2 for sphere inside fft box. (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  divgq0= value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq. Used if q = Gamma
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  izero=if 1, unbalanced components of V(q,g) are set to zero
!!  hyb_mixing=hybrid mixing coefficient for the Fock contribution
!!  hyb_mixing_sr=hybrid mixing coefficient for the short-range Fock contribution
!!  hyb_range_fock=hybrid range for separation
!!  nfft=Total number of FFT grid points.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  vqg(nfft)=4pi/(G+q)**2, G=0 component is set to divgq0/pi if q = Gamma.
!!
!! NOTES
!!  This routine operates on the full FFT mesh. DO NOT PASS MPI_TYPE
!!  One can easily implemente MPI-FFT by just calling this routine and then
!!  extracting the G-vectors treated by the node.
!!
!! PARENTS
!!      m_fock_getghc
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine bare_vqg_abinit_v8(qphon,gsqcut,gmet,izero,hyb_mixing,hyb_mixing_sr,hyb_range_fock,nfft,nkpt_bz,ngfft,ucvol,vqg)


!!This section has been created automatically by the script Abilint (TD).
!!Do not modify the following lines by hand.
!#undef ABI_FUNC
!#define ABI_FUNC 'bare_vqg'
!!End of the abilint section


 use interfaces_53_ffts, only: zerosym

 implicit none

 !!!!!!!!!!!!!!!!!!! chen !!!!!!!!!!!!!!!
 integer, parameter :: dp = 8
 real(8), parameter :: pi = 3.14159265359d0
 real(8):: two=2.d0, & 
           one=1.d0, & 
           four_pi = 4.d0*pi, & 
           half=0.5d0,  & 
           three=3.d0,  & 
           piinv = 1.d0/pi, & 
           two_pi = 2.d0*pi, & 
           zero = 0d0, & 
           tol12 = 1e-12, & 
           tol8 = 1e-8
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,nfft,nkpt_bz
 real(dp),intent(in) :: gsqcut,hyb_mixing,hyb_mixing_sr,hyb_range_fock,ucvol
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),qphon(3)
 real(dp),intent(out) ::  vqg(nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex1=1
 integer :: i1,i2,i23,i3,id1,id2,id3
 integer :: ig,ig1min,ig1,ig1max,ig2,ig2min,ig2max,ig3,ig3min,ig3max
 integer :: ii,ii1,ing,n1,n2,n3,qeq0,qeq05
 real(dp),parameter :: tolfix=1.000000001e0_dp ! Same value as the one used in hartre
 real(dp) :: cutoff,den,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3,rcut,divgq0
 character(len=100) :: msg
!arrays
 integer :: id(3)
 real(dp),allocatable :: gq(:,:)

! *************************************************************************

 if (abs(hyb_mixing_sr)>tol8.and.abs(hyb_range_fock)<tol8) then
   msg='SR mixing<>0 while range separation=0!'
 !  MSG_BUG(msg)
 end if

!Treatment of the divergence at q+g=zero
!For the time being, only Spencer-Alavi scheme...
 rcut= (three*nkpt_bz*ucvol/four_pi)**(one/three)
 divgq0= two_pi*rcut**two
!divgq0=zero
!Initialize a few quantities
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 cutoff=gsqcut*tolfix
 vqg=zero

!Some peculiar values of q
 qeq0=0; if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1
 qeq05=0
 if (qeq0==0) then
   if (abs(abs(qphon(1))-half)<tol12.or.abs(abs(qphon(2))-half)<tol12.or. &
&   abs(abs(qphon(3))-half)<tol12) qeq05=1
 end if

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 allocate(gq(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qphon(ii)
   end do
 end do
 ig1max=-1;ig2max=-1;ig3max=-1
 ig1min=n1;ig2min=n2;ig3min=n3

 id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

 ! Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
     gqgm12=gq(2,i2)*gmet(1,2)*2
     gqg2p3=gqgm13+gqgm12

     i23=n1*(i2-1 +(n2)*(i3-1))
     ! Do the test that eliminates the Gamma point outside of the inner loop
     ii1=1
     if (i23==0 .and. qeq0==1  .and. ig2==0 .and. ig3==0) then
       ii1=2
       ! value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq
       vqg(1+i23)=hyb_mixing*divgq0

!      Note the combination of Spencer-Alavi and Erfc screening
       if (abs(hyb_range_fock)>tol8)then
         vqg(1+i23)=vqg(1+i23)+hyb_mixing_sr*(pi/hyb_range_fock**2)
!        This would give a combination of Spencer-Alavi and Erfc screening,
!        unfortunately, it modifies also the tests for pure HSE06, so was not retained.
!        vqg(1+i23)=vqg(1+i23)+hyb_mixing_sr*min(divgq0,pi/(hyb_range_fock**2))
       endif

     end if

     ! Final inner loop on the first dimension (note the lower limit)
     do i1=ii1,n1
       gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
       ii=i1+i23

       if(gs<=cutoff)then
         ! Identify min/max indexes (to cancel unbalanced contributions later)
         ! Count (q+g)-vectors with similar norm
         if ((qeq05==1).and.(izero==1)) then
           ig1=i1-(i1/id1)*n1-1
           ig1max=max(ig1max,ig1); ig1min=min(ig1min,ig1)
           ig2max=max(ig2max,ig2); ig2min=min(ig2min,ig2)
           ig3max=max(ig3max,ig3); ig3min=min(ig3min,ig3)
         end if

         den=piinv/gs

!        Spencer-Alavi screening
         if (abs(hyb_mixing)>tol8)then
           vqg(ii)=vqg(ii)+hyb_mixing*den*(one-cos(rcut*sqrt(four_pi/den)))
!&         vqg(ii)=vqg(ii)+hyb_mixing*den
         endif
!        Erfc screening
         if (abs(hyb_mixing_sr)>tol8) then
           vqg(ii)=vqg(ii)+hyb_mixing_sr*den*(one-exp(-pi/(den*hyb_range_fock**2)))
!          This other possibility combines Erfc and Spencer-Alavi screening in case rcut is too small or hyb_range_fock too large
!          if(divgq0<pi/(hyb_range_fock**2))then
!            vqg(ii)=vqg(ii)+hyb_mixing_sr*den*&
!&             (one-exp(-pi/(den*hyb_range_fock**2)))*(one-cos(rcut*sqrt(four_pi/den)))
!          endif
         endif

       end if ! Cut-off
     end do ! End loop on i1
   end do ! End loop on i2
 end do ! End loop on i3

 if (izero==1) then
   ! Set contribution of unbalanced components to zero
   if (qeq0==1) then !q=0
     call zerosym(vqg,cplex1,n1,n2,n3)
   else if (qeq05==1) then
     !q=1/2; this doesn't work in parallel
     ig1=-1;if (mod(n1,2)==0) ig1=1+n1/2
     ig2=-1;if (mod(n2,2)==0) ig2=1+n2/2
     ig3=-1;if (mod(n3,2)==0) ig3=1+n3/2
     if (abs(abs(qphon(1))-half)<tol12) then
       if (abs(ig1min)<abs(ig1max)) ig1=abs(ig1max)
       if (abs(ig1min)>abs(ig1max)) ig1=n1-abs(ig1min)
     end if
     if (abs(abs(qphon(2))-half)<tol12) then
       if (abs(ig2min)<abs(ig2max)) ig2=abs(ig2max)
       if (abs(ig2min)>abs(ig2max)) ig2=n2-abs(ig2min)
     end if
     if (abs(abs(qphon(3))-half)<tol12) then
       if (abs(ig3min)<abs(ig3max)) ig3=abs(ig3max)
       if (abs(ig3min)>abs(ig3max)) ig3=n3-abs(ig3min)
     end if
     call zerosym(vqg,cplex1,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3)
   end if
 end if

 DEALLOCATE(gq)

end subroutine bare_vqg_abinit_v8
!!***
