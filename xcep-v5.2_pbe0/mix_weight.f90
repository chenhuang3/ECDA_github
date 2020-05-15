!
! Mix embedding potential to accelerate the convergence of 
! embedding potential, pulay mixing is used 
!
! iopt should start from 1
!
! we iterate   x=F[x]
!  x is current trial 
!  x_new is the F[x]
! 
! on output 
!   x_new will be updated
!
! (1) if iopt >= 2 start doing pulay mixing
! (2) if iopt == 1 we simply mixing 
!

subroutine mix_weight(npt,iopt,uemb,uemb_new,work,npulay)

 implicit none 

 integer :: npt
 integer :: iopt, npulay, ii, jj, neff, lapack_info
 real(kind=8),intent(inout) :: uemb_new(npt)
 real(kind=8),intent(in)    :: uemb(npt)
 real(8) ::                    work(npt,npulay,2)     ! the work(:,1,:) is residual 
                                                      ! the work(:,2,:) is the vector
 ! >>>>>>> local vars <<<<<<<<< !

 real(kind=8),allocatable :: pulay_matrix(:,:), & 
                             lapack_work(:), matb(:)
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

 print *,'pulay_work array: <F[u], u> '
 do ii=1,npulay
   do jj=1,npt
     write(6,'(f12.6,a,f12.6)',advance='no')work(jj,ii,1),' ',work(jj,ii,2)
   enddo
   write(6,'(a)',advance="YES")""
 enddo

 ! display ==================
 write(6,'(a)',advance='no')'mix: old(sub_etotal): '
 do ii=1,npt
   write(6,'(f14.6)',advance='no')uemb(ii)
 enddo
 write(6,*)''
 write(6,'(a)',advance='no')'mix: new(sub_etotal): '
 do ii=1,npt
   write(6,'(f14.6)',advance='no')uemb_new(ii)
 enddo
 write(6,*)""
 write(6,*)'mix: norm(diff_etotal): ',sqrt(sum((uemb-uemb_new)**2))

! debug 
!   ! simple mixing -------------
!   print *,'simple mixing of uemb.'
!   uemb = (uemb_new + uemb)/2.d0
!   print *,'after mixing, min/max(uemb): ', minval(uemb),maxval(uemb),' amp(uemb): ', maxval(uemb)-minval(uemb)
!   return 
! end of debug

 if ( iopt==1 ) then 

   ! simple mixing -------------
   print *,'simple mixing of uemb.'
   uemb_new = (uemb_new + uemb)/2.d0

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
   pulay_matrix(neff+1,:)      = 1.0d0
   pulay_matrix(:,neff+1)      = 1.0d0
   pulay_matrix(neff+1,neff+1) = 0.0d0

   print *,'pulay_matrix: '
   do ii=npulay-neff+1,npulay
     do jj=npulay-neff+1,npulay
       pulay_matrix(ii-npulay+neff,jj-npulay+neff) = sum(work(:,ii,1)*work(:,jj,1))
       write(6,'(f16.6)',advance='no') pulay_matrix(ii-npulay+neff,jj-npulay+neff)
     enddo
     write(6,'(a)',advance='yes')""
   enddo

   !
   ! solve for AX=b using lapack
   !
   ! see equation (5) in
   ! http://scitation.aip.org/content/aip/journal/jcp/137/5/10.1063/1.4740249
   !
   allocate(matb(neff+1))
   matb(1:neff) = 0.0d0
   matb(neff+1) = 1.0d0
   !call DPOSV('U',neff+1,1,pulay_matrix,neff+1,matb,neff+1,lapack_info)
   call DGESV(neff+1,1,pulay_matrix,neff+1,lapack_ipiv,matb,neff+1,lapack_info)
   pulay_matrix(:,neff+1) = matb(:)
   deallocate(matb)
   write(6,'(a)',advance='no')'mix: pulay coeffs: '
   do ii=1,neff
     write(6,'(f12.4)',advance='no')pulay_matrix(ii,neff+1)
   enddo
   write(6,'(a,f16.8)',advance='yes')'  sum(c):',sum(pulay_matrix(1:neff,neff+1))

   ! pulay mixing
   uemb_new = 0.d0
   do ii=1,neff
     ! if pulay_beta = 1.0 then we are mixing F[n] 
     ! if pulay_beta < 1.0 then we are mixing beta*F[n] + (1-beta)*n
     pulay_beta = 0.1d0
     uemb_new = uemb_new + & 
       (work(:,npulay-neff+ii,2) + pulay_beta*work(:,npulay-neff+ii,1))*pulay_matrix(ii,neff+1)
   enddo

   deallocate(pulay_matrix)
   deallocate(lapack_ipiv)
   deallocate(lapack_work)

 endif

 write(6,'(a)',advance='no')'mix: after mixing: '
 do ii=1,npt
   write(6,'(f14.6)',advance='no')uemb_new(ii)
 enddo
 write(6,'(a)',advance="YES")""

end subroutine 
