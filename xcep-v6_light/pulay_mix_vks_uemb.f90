subroutine pulay_mix_vks_uemb(npt,iopt,vks,vks_new,uemb,uemb_new, & 
                              work,work_emb,npulay,pulay_beta)

 use mpi 

 implicit none 

 integer :: npt, rank
 integer :: iopt, npulay

 real(kind=8),intent(inout) :: uemb_new(npt),vks_new(npt)
 real(kind=8),intent(in)    :: uemb(npt),vks(npt)

 real(8) :: work(npt,npulay,2)     ! the work(:,1,:) is residual for vks
                                   ! the work(:,2,:) is the vector for vks 
 real(8) :: work_emb(npt,npulay,2) ! the work_emb(:,1,:) is residual for uemb
                                   ! the work_emb(:,2,:) is the vector for uemb


 ! >>>>>>> local vars <<<<<<<<< !
 real(kind=8),allocatable :: pulay_matrix(:,:), & 
                             lapack_work(:), matb(:)
 integer, allocatable     :: lapack_ipiv(:)
 integer                  :: ierr, lapack_info, neff, ii, jj
 real(kind=8)             :: pulay_beta, res_vks, res_emb

 call MPI_Comm_rank ( MPI_COMM_WORLD, rank, ierr )


 ! push into history array for pulay mixing 
 do ii=1,npulay-1
   work(:,ii,1) = work(:,ii+1,1)
   work(:,ii,2) = work(:,ii+1,2)
   !
   work_emb(:,ii,1) = work_emb(:,ii+1,1)
   work_emb(:,ii,2) = work_emb(:,ii+1,2)
 enddo
 !
 work(:,npulay,1) = vks_new - vks ! store new residual 
 work(:,npulay,2) = vks           ! store vector
 !
 work_emb(:,npulay,1) = uemb_new - uemb ! store new residual 
 work_emb(:,npulay,2) = uemb            ! store vector


 if ( iopt==1 ) then 
   ! simple mixing -------------
   uemb_new = uemb_new*pulay_beta + uemb*(1.d0-pulay_beta)
   vks_new  = vks_new *pulay_beta + vks *(1.d0-pulay_beta)
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
   !
   ! note that: the vector is [vks, uemb_1, uemb_2, uemb_3, ... uemb_n]
   !
   do ii=npulay-neff+1,npulay
     do jj=npulay-neff+1,npulay

       res_vks = sum(work(:,ii,1)*work(:,jj,1))         ! residule for vks
       res_emb = sum(work_emb(:,ii,1)*work_emb(:,jj,1)) ! residule for embedding potential 
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,res_emb,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

       pulay_matrix(ii-npulay+neff,jj-npulay+neff) = res_vks + res_emb
     enddo
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
!   if(rank==0) write(6,'(a)',advance='no')'pulay mixing coeffs: '
   do ii=1,neff
!     if(rank==0) write(6,'(f12.4)',advance='no')pulay_matrix(ii,neff+1)
   enddo
!   if(rank==0) write(6,'(a,f16.8)',advance='yes')'  sum(c):',sum(pulay_matrix(1:neff,neff+1))
   !
   ! pulay mixing
   !
   vks_new  = 0.d0 
   uemb_new = 0.d0
   do ii=1,neff
     ! if pulay_beta = 1.0 then we are mixing F[n] 
     ! if pulay_beta < 1.0 then we are mixing beta*F[n] + (1-beta)*n

     uemb_new = uemb_new + & 
       (work_emb(:,npulay-neff+ii,2) + & 
        pulay_beta*work_emb(:,npulay-neff+ii,1))*pulay_matrix(ii,neff+1)

     vks_new  = vks_new + & 
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