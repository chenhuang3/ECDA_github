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

subroutine pulay_mix_den_metric(npt,iopt,rhor,uemb,uemb_new,work,npulay,pulay_beta)

 use mpi 

 implicit none 

 integer :: npt, rank
 integer :: iopt, npulay, ii, jj, neff, lapack_info
 real(kind=8),intent(inout) :: uemb_new(npt)
 real(kind=8),intent(in)    :: uemb(npt), rhor(npt)
 real(8) ::                    work(npt,npulay,2)     ! the work(:,1,:) is residual 
                                                      ! the work(:,2,:) is the vector
 ! >>>>>>> local vars <<<<<<<<< !

 real(kind=8),allocatable :: pulay_matrix(:,:), & 
                             lapack_work(:), matb(:)
 integer, allocatable     :: lapack_ipiv(:)
 integer                  :: ierr, i
 real(kind=8)             :: pulay_beta,  metric(npt)

 call MPI_Comm_rank ( MPI_COMM_WORLD, rank, ierr )


 ! Pulay mixing ===============================

 metric = (rhor + 1e-4/10) / (rhor + 1e-4)

 ! push into history array for pulay mixing 
 do ii=1,npulay-1
   work(:,ii,1) = work(:,ii+1,1)
   work(:,ii,2) = work(:,ii+1,2)
 enddo
 work(:,npulay,1) = uemb_new - uemb ! store new residual 
 work(:,npulay,2) = uemb            ! store vector


 if ( iopt==1 ) then 
   ! simple mixing -------------
   uemb_new = uemb_new*pulay_beta + uemb*(1.d0-pulay_beta)
 else
   ! pulay mixing ---------------
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

   do ii=npulay-neff+1,npulay
     do jj=npulay-neff+1,npulay
       pulay_matrix(ii-npulay+neff,jj-npulay+neff) = sum(work(:,ii,1)*work(:,jj,1)*metric)
     enddo
   enddo

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

   ! pulay mixing
   uemb_new = 0.d0
   do ii=1,neff
     ! if pulay_beta = 1.0 then we are mixing F[n] 
     ! if pulay_beta < 1.0 then we are mixing beta*F[n] + (1-beta)*n
     uemb_new = uemb_new + & 
       (work(:,npulay-neff+ii,2) + & 
        pulay_beta*work(:,npulay-neff+ii,1))*pulay_matrix(ii,neff+1)
   enddo

   deallocate(pulay_matrix)
   deallocate(lapack_ipiv)
   deallocate(lapack_work)

 endif

! write(6,'(a)',advance='no')'mix: after mixing: '
! do ii=1,npt
!   write(6,'(f14.6)',advance='no')uemb_new(ii)
! enddo
! write(6,'(a)',advance="yes")""

end subroutine 
