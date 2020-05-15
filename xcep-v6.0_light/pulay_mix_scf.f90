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
!
! The code follows the idea in paper:
!
! "Preconditioning of self-consistent-field cycles in density-functional theory: The extrapolar method"
! Anglade and Gonze (PRB 78, 045126 2008)
! 
! See equations (14) and (15) in above paper.
! 
! The residule V_out - V_in is preconditioned by the inverse of the dielectric function before
! being used by Pulay mixing. 
!
! The model dielectric function follows ABINIT 
! https://docs.abinit.org/variables/gstate/#dielng
! 
!                 1 + L^2 * q^2
!  diel(q^2)= ------------------------
!              (1/diemac + L^2 * q^2)
!
!  Note that: at q=0, diel(q) is imposed to 1. 
!
!
subroutine pulay_mix_scf(n1,n2,n3,npt,iopt,uemb,uemb_new,work,npulay,pulay_beta,qvec,ucvol)

 use mpi 

 implicit none 

 integer :: npt, rank, n1,n2,n3
 integer :: iopt, npulay, ii, jj, neff, lapack_info

 real(8)                    :: diemac = 10.d0 
 real(8)                    :: dielen = 5.d0  ! bohr 
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

 ! precondition residual 
 call precond_scf_residual(n1,n2,n3,npt,qvec,diemac,dielen,uemb_new-uemb,px)
 work(:,npulay,1) = px     ! store new preconditioned residual 
 work(:,npulay,2) = uemb   ! store new vector


 if ( iopt==1 ) then 
   ! simple mixing -------------
   uemb_new = uemb_new*0.1d0 + uemb*0.9d0
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
       pulay_matrix(ii-npulay+neff,jj-npulay+neff) = sum(work(:,ii,1)*work(:,jj,1))
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

   
   !===================================
   ! mixing 
   uemb_new = 0.d0
   do ii=1,neff
     ! if pulay_beta = 1.0 then we are mixing F[n] 
     ! if pulay_beta < 1.0 then we are mixing beta*F[n] + (1-beta)*n
     vec   = work(:,npulay-neff+ii,2)  ! vecor 
     resid = work(:,npulay-neff+ii,1)  ! residule 
     coeff = pulay_matrix(ii,neff+1)   ! pulay coefficient 

     uemb_new = uemb_new + (vec + pulay_beta*resid)*coeff
   enddo

   deallocate(pulay_matrix)
   deallocate(lapack_ipiv)
   deallocate(lapack_work)

 endif


end subroutine 
