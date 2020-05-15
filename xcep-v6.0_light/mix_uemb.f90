!
! Mix embedding potential to accelerate the convergence of 
! embedding potential, pulay mixing is used 
!
! iopt should start from 1
!
! (1) if iopt >= 2 start doing pulay mixing
! (2) if iopt == 1 we simply mixing 
!

subroutine mix_uemb(npt,iopt,uemb,uemb_new,work,npulay)

 implicit none 

 integer :: npt
 integer :: iopt, npulay, ii, jj, neff, lapack_info
 real(kind=8),intent(in)    :: uemb_new(npt)
 real(kind=8),intent(inout) :: uemb(npt), & 
                               work(npt,npulay,2)     ! the work(:,1,:) is residual 
                                                      ! the work(:,2,:) is the vector
 ! >>>>>>> local vars <<<<<<<<< !

 real(kind=8),allocatable :: pulay_matrix(:,:), & 
                             lapack_work(:)
 integer, allocatable     :: lapack_ipiv(:)                             
 real(kind=8)             :: pulay_beta 

 ! Pulay mixing ===============================

 print *,''
 print *,'    ------------------------- mix_uemb() --------------------------'
 print *,''

 ! push into history array for pulay mixing 
 do ii=1,npulay-1
   work(:,ii,1) = work(:,ii+1,1)
   work(:,ii,2) = work(:,ii+1,2)
 enddo
 work(:,npulay,1) = uemb_new - uemb ! store new residual 
 work(:,npulay,2) = uemb            ! store vector

 write(6,'(a,2g16.8,a,g16.8)')'min/max(uemb_old): ', minval(uemb),maxval(uemb),& 
  ' amp(uemb_old): ', maxval(uemb)-minval(uemb)
 write(6,'(a,2g16.8,a,g16.8)')'min/max(uemb_new): ', minval(uemb_new),maxval(uemb_new), & 
  ' amp(uemb_new): ', maxval(uemb_new)-minval(uemb_new)

! debug 
   ! simple mixing -------------
   print *,'simple mixing of uemb.'
   uemb = (uemb_new + uemb)/2.d0
   print *,'after mixing, min/max(uemb): ', minval(uemb),maxval(uemb),' amp(uemb): ', maxval(uemb)-minval(uemb)
   return 
! end of debug

 if ( iopt==1 ) then 

   ! simple mixing -------------
   print *,'simple mixing of uemb.'
   uemb = (uemb_new + uemb)/2.d0

 else

   ! pulay mixing ---------------
   print *,'pulay mixing of uemb. npulay:',npulay
   if ( iopt>=npulay ) then 
     neff = npulay
   else
     neff = iopt
   endif

   allocate(pulay_matrix(neff+1,neff+1))
   allocate(lapack_ipiv(neff+1))
   allocate(lapack_work(neff+1))

   ! make pulay matrix
   pulay_matrix(neff+1,:)      = -1.0d0
   pulay_matrix(:,neff+1)      = -1.0d0
   pulay_matrix(neff+1,neff+1) =  0.0d0

   do ii=npulay-neff+1,npulay
     do jj=npulay-neff+1,npulay
       pulay_matrix(ii-npulay+neff,jj-npulay+neff) = sum(work(:,ii,1)*work(:,jj,1))
     enddo
   enddo

   ! invert pulay_matrix 
   call dgetrf(neff+1,neff+1,pulay_matrix,neff+1,lapack_ipiv,lapack_info)
   call dgetri(neff+1,pulay_matrix,neff+1,lapack_ipiv,lapack_work,neff+1,lapack_info)
   write(*,*),'pulay coeffs: ',pulay_matrix(1:neff,neff+1)

   ! pulay mixing
   uemb = 0.d0
   do ii=1,neff
     ! if pulay_beta = 1.0 then we are mixing F[n] 
     ! if pulay_beta < 1.0 then we are mixing beta*F[n] + (1-beta)*n
     pulay_beta = 0.2d0
     uemb = uemb - & 
       (work(:,npulay-neff+ii,2) + pulay_beta*work(:,npulay-neff+ii,1))*pulay_matrix(ii,neff+1)
   enddo

   deallocate(pulay_matrix)
   deallocate(lapack_ipiv)
   deallocate(lapack_work)

 endif

 write(6,'(a,2g16.8,a,g16.8)') & 
  'after mixing, min/max(uemb): ', minval(uemb),maxval(uemb), & 
  ' amp(uemb): ', maxval(uemb)-minval(uemb)

end subroutine 
