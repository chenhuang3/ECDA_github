!
! Mix embedding potential to accelerate the convergence of 
! embedding potential, pulay mixing is used 
!
! Based on paper:
!   "An efficient and robust technique for achieving self consistency in electronic structure calculations"
!   Bowler and Gillan, Chemical Physics Letters 325 2000 473–476
!
! we iterate   x=F[x]
!  x is current trial 
!  x_new is the F[x]
!
! on first entry
!  task = 'START'
!   
! on output 
!   uemb_new is updated
!

subroutine pulay_mix_zp_gr(task,npt,iwork,uemb,uemb_new,work,npulay,pulay_beta)

 use mpi 

 implicit none 

 character(len=7)           :: task 
 integer                    :: npt, rank, & 
                               npulay, ii, jj, neff, lapack_info
 real(kind=8),intent(inout) :: uemb_new(npt)
 real(kind=8),intent(in)    :: uemb(npt)
 real(8)                    :: delta_u(npt)
 integer ::                    iwork(5) 
 real(8) ::                    work(npt,npulay,2)     ! the work(:,1,:) is residual 
                                                      ! the work(:,2,:) is the vector
 ! >>>>>>> local vars <<<<<<<<< !

 real(kind=8),allocatable :: pulay_matrix(:,:), & 
                             lapack_work(:), matb(:)
 integer, allocatable     :: lapack_ipiv(:)
 integer                  :: ierr 
 real(kind=8)             :: pulay_beta, & 
                             mix_coeff, amp, amp_new

 call MPI_Comm_rank ( MPI_COMM_WORLD, rank, ierr )

 print *,'--------- pulay_mix_gr() -----------'

 !===================
 ! start GR-pulay 
 !===================

 if (task(1:5) == 'START') then 
   work(:,npulay,1) = uemb_new - uemb ! store new residual 
   work(:,npulay,2) = uemb            ! store vector

   iwork(1)=1 ! number of history vectors and residules

   print *,'simple mixing ........'
   delta_u = uemb_new - uemb
   ! make sure that the norm of uemb does not change over 1%
   mix_coeff = 0.01d0 / maxval(delta_u)
   print *,'mix_coeff: ',mix_coeff
   uemb_new = uemb + mix_coeff*delta_u

   task = 'NEW_RES'
   return 
 endif 


 !
 ! update work with optimal solution and its residule 
 !
 if (task(1:5) == 'OPT_X') then 
   work(:,npulay,1) = uemb_new - uemb ! store new residual 
   work(:,npulay,2) = uemb            ! store vector

   delta_u = uemb_new - uemb
   mix_coeff = min(0.01d0, 0.005d0/maxval(delta_u))
   uemb_new = uemb + delta_u*mix_coeff

   task = 'NEW_RES'
   return 
 endif 


 !
 ! pulay mixing 
 !
 if (task(1:7)=='NEW_RES') then 

   ! push into history array for pulay mixing 
   do ii=1,npulay-1
     work(:,ii,1) = work(:,ii+1,1)
     work(:,ii,2) = work(:,ii+1,2)
   enddo
   work(:,npulay,1) = uemb_new - uemb ! store new residual 
   work(:,npulay,2) = uemb            ! store vector

   iwork(1) = iwork(1)+1

   task = 'OPT_X' ! this is the optimal X 

   if (rank==0) then 
     print *,'mix: norm(new_vec-old_vec): ',sqrt(sum((uemb-uemb_new)**2))
     print *,'pulay mixing of uemb. npulay:',npulay
   endif 

   ! pulay mixing ---------------
   if ( iwork(1)>=npulay ) then 
     neff = npulay
   else
     neff = iwork(1)
   endif
  
   allocate(pulay_matrix(neff+1,neff+1))
   allocate(lapack_ipiv(neff+1))
   allocate(lapack_work(neff+1))
  
   ! make pulay matrix
   pulay_matrix(neff+1,:)      = 1.0d0
   pulay_matrix(:,neff+1)      = 1.0d0
   pulay_matrix(neff+1,neff+1) = 0.0d0
  
   if (rank==0) print *,'pulay_matrix: '
   do ii=npulay-neff+1,npulay
     do jj=npulay-neff+1,npulay
       pulay_matrix(ii-npulay+neff,jj-npulay+neff) = sum(work(:,ii,1)*work(:,jj,1))
       if(rank==0) write(6,'(f16.6)',advance='no') pulay_matrix(ii-npulay+neff,jj-npulay+neff)
     enddo
     if(rank==0)write(6,'(a)',advance='yes')""
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


   if(rank==0) write(6,'(a)',advance='no')'pulay mixing coeffs: '
   do ii=1,neff
     if(rank==0) write(6,'(f12.4)',advance='no')pulay_matrix(ii,neff+1)
   enddo
   if(rank==0) write(6,'(a,f16.8)',advance='yes')'  sum(c):',sum(pulay_matrix(1:neff,neff+1))
 

   ! pulay mixing
   uemb_new = 0.d0
   do ii=1,neff
     ! if pulay_beta = 1.0 then we are mixing F[n] 
     ! if pulay_beta < 1.0 then we are mixing beta*F[n] + (1-beta)*n
     pulay_beta = 0.0d0
     uemb_new = uemb_new + & 
     ( work(:,npulay-neff+ii,2) + & 
       pulay_beta*work(:,npulay-neff+ii,1))*pulay_matrix(ii,neff+1)
   enddo

   deallocate(pulay_matrix)
   deallocate(lapack_ipiv)
   deallocate(lapack_work)

 endif 



end subroutine pulay_mix_zp_gr
