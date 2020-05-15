!================================================================
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
!================================================================

subroutine pulay_mix_kerker_dfet(n1,n2,n3,npt,iopt,uemb,uemb_new,work,npulay,pulay_beta,qvec,ucvol)

 use mpi 

 implicit none 

 integer :: npt, rank, n1,n2,n3
 integer :: iopt, npulay, ii, jj, neff, lapack_info
 real(kind=8),intent(inout) :: uemb_new(npt)
 real(kind=8),intent(in)    :: uemb(npt), ucvol, & 
                               pulay_beta, & 
                               qvec(3,(n1/2+1)*n2*n3)
 real(8) ::                    work(npt,npulay,2)     ! the work(:,1,:) is residual 
                                                      ! the work(:,2,:) is the vector
 ! >>>>>>> local vars <<<<<<<<< !
 real(kind=8),allocatable :: pulay_matrix(:,:), & 
                             lapack_work(:), matb(:)
 integer, allocatable     :: lapack_ipiv(:)
 integer                  :: ierr, ix, iy, iz
 real(8)                  :: qq((n1/2+1)*n2*n3), & 
                             qq_min, qq_max, qq1, factor, qq_opt 
 real(8)                  :: coeff, vec(npt), resid(npt), px(npt)



 call MPI_Comm_rank ( MPI_COMM_WORLD, rank, ierr )


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
       
       !pulay_matrix(ii-npulay+neff,jj-npulay+neff) = sum(work(:,ii,1)*work(:,jj,1))
       
       !=============================================================
       ! Wave-vector depending metric --- long wavelength is focused 
       !
       qq = sum(qvec**2,1)
       qq_max = maxval(qq)

       ! get the smallest nonzero qq
       qq_min = 1e10
       do ix=1,((n1/2)+1)*n2*n3
         qq1 = qq(ix)
         if ( qq1>1e-8 .and. qq_min>qq1 ) then 
           qq_min = qq1
         endif 
       enddo

       ! compute optimal qq to make long-wave is factor=20 times larger than short-wave 
       factor = 20.d0
       qq_opt = (factor-1.d0)*qq_max*qq_min/(qq_max-factor*qq_min)
       !!print *,'qq_opt: ',qq_opt
       call precond_vks(n1,n2,n3,n1*n2*n3,qvec,qq_opt,qq_min,work(:,jj,1),px)

       pulay_matrix(ii-npulay+neff,jj-npulay+neff) = sum(work(:,ii,1)*px)
     enddo
   enddo



   !=======================================
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

   !if(rank==0) write(6,'(a)',advance='no')'pulay mixing coeffs: '
   !do ii=1,neff
   !  if(rank==0) write(6,'(f12.4)',advance='no')pulay_matrix(ii,neff+1)
   !enddo
   !if(rank==0) write(6,'(a,f16.8)',advance='yes')'  sum(c):',sum(pulay_matrix(1:neff,neff+1))


   
   !===================================
   ! mixing 
   uemb_new = 0.d0
   do ii=1,neff
     ! if pulay_beta = 1.0 then we are mixing F[n] 
     ! if pulay_beta < 1.0 then we are mixing beta*F[n] + (1-beta)*n

     vec   = work(:,npulay-neff+ii,2)  ! vecor 
     resid = work(:,npulay-neff+ii,1)  ! residule 
     coeff = pulay_matrix(ii,neff+1)   ! pulay coefficient 

     ! precondion resid
     !call precond_kerker(n1,n2,n3,n1*n2*n3,qvec,1.d0,resid,px)
     !uemb_new = uemb_new + (vec + pulay_beta*px)*coeff

     uemb_new = uemb_new + (vec + pulay_beta*resid)*coeff
   enddo

   deallocate(pulay_matrix)
   deallocate(lapack_ipiv)
   deallocate(lapack_work)

 endif



contains 


 !===========================================================
 ! precondition x
 !
 subroutine precond_vks(n1,n2,n3,nfft,qvec,q0,qq_min,x,Px)
  implicit none 

  integer :: n1,n2,n3,nfft, dim1, ix, iy, iz
  real(8) :: x(nfft), Px(nfft), & 
             q3d(3,n1/2+1,n2,n3), & 
             qvec(3,(n1/2+1)*n2*n3), & 
             factor_k, q0, qq_min, & 
             rtmp(n1,n2,n3), metric, & 
             qq((n1/2+1),n2,n3),qq1

  complex(kind=8) :: fft1(n1/2+1,n2,n3), & 
                     fft2(n1/2+1,n2,n3)

  dim1=(n1/2+1)
  call FFT(n1,n2,n3,reshape(x,(/n1,n2,n3/)),fft1,1)

  ! convert qvec to 3-dimension 
  q3d(1,:,:,:) = reshape(qvec(1,:),(/dim1,n2,n3/))
  q3d(2,:,:,:) = reshape(qvec(2,:),(/dim1,n2,n3/))
  q3d(3,:,:,:) = reshape(qvec(3,:),(/dim1,n2,n3/))

  qq = sum(q3d**2,1)

  do ix=1,dim1
    do iy=1,n2
      do iz=1,n3
        qq1 = qq(ix,iy,iz)

        if ( qq1>1e-8 ) then 
          metric = (qq1+q0)/qq1
        else 
          metric = (qq_min+q0)/qq_min
        endif 
        fft2(ix,iy,iz) = fft1(ix,iy,iz) / metric
       !! fft2(ix,iy,iz) = fft1(ix,iy,iz) 
       ! fft2(ix,iy,iz) = fft1(ix,iy,iz) * qq1/(qq1+q0)
      enddo
    enddo
  enddo

  call FFT(n1,n2,n3,rtmp,fft2,-1)
  Px = reshape(rtmp,(/n1*n2*n3/))
 end subroutine 

end subroutine 
