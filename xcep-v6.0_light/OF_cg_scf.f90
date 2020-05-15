!
! Use L-BFGS or conjugate gradient method to direct minimize total energy w.r.t. electron density 
! see ECDA/OFenv paper for details 
!
subroutine OF_cg_scf(logOF,n1,n2,n3,nfft,nspin,vw_lam,vw_reg, & 
                 qvec,ucvol,vext,rho,user_init,chempot)
 
 use mpi

 implicit none 

 logical :: user_init

 integer :: n1,n2,n3, & 
            nfft, nspin, & 
            logOF, & 
            isp 

 real(8) :: vext(nfft,nspin), & 
            cTF = 3.d0/10.d0*(3.d0*3.1415926d0**2)**(2.d0/3.d0), & 
            rho(nfft,nspin), & 
            new_sq_rho(nfft), & 
            vOF(nfft),vw_lam,vw_reg,  &
            qvec(3,(n1/2+1)*n2*n3), & 
            resid, ucvol, dvol, chempot, & 
            q(nspin), q_tmp, & 
            vw_pot(nfft), vw_energy, & 
            tf_pot(nfft), tf_energy, & 
            phi(nfft), phi_tmp(nfft), phi_reg(nfft), &  
            f(nfft), dE_dphi(nfft), & 
            cg_beta, cg_alpha

 logical :: print_info = .false., & 
            do_precond

 integer :: opt_method, & 
            iter,k, & 
            ipt, ierr

 real(8) :: dE_dtheta, etotal, etotal2, & 
            wf(nfft), etotal_hist(2), etol, & 
            etotal_tmp, & 
            etotal_tmp2, & 
            dtmp, theta0,  & 
            arr_tmp(nfft),& 
            AA, BB, orth 


 integer,parameter :: neig=1
 integer :: ndiag, nline, niter, m, ii, ib, ib2, j, lapack_info 
 real(8) :: ke, lap(nfft), & 
            v(nfft,neig), coeff, lam, dconv, theta, & 
            v0(nfft), Ax(nfft), pg_old(nfft), & 
            pg(nfft), g(nfft), g_old(nfft), d_proj(nfft), d(nfft), & 
            Ad(nfft), R(neig,neig), cg_tol, & 
            Ax_hist(nfft,neig), a(neig), eigenvalues(neig), & 
            lapack_w(neig), lapack_WORK(3*neig-1)

 

 ! >>>>>>>>>>>>>> function <<<<<<<<<<<<!

 write(logOF,'(a)')NEW_LINE('a')//'enter OF_solver_scf()'
 dvol = ucvol / dble(nfft)


 write(logOF,*)"OF_cg() is only partially coded, does not work, stop"
 write(6,*)"OF_cg() is only partially coded, does not work, stop"
 call flush(logOF)
 call flush(6)
 stop
 
 do isp=1,nspin 
   q(isp) = dvol*sum(rho(:,isp))
 enddo
 write(logOF,'(a,f12.5,a,f8.4,a,es12.4)') & 
   'q: ',q(1:nspin),' vw_lambda:',vw_lam,' vw_reg: ',vw_reg
 write(logOF,'(a,2es12.4)')'min/max(vext): ',minval(vext),maxval(vext)
 call flush(logOF)


 niter = 10000 ! number of density mixing 
 ndiag = 200   ! number of diag steps 
 nline = 100   ! number of line searches 

 ! code based on DOI: 10.1504/IJCSE.2006.012774
 ! "Conjugate-gradient eigenvalue solvers in computing electronic properties of 
 !  nanostructure architectures", Lin-Wang Wang and coworkers 
 ! International Journal of Computational Science and Engineering, 2006 Vol 2, pp.205
 !
 do isp=1,nspin 

   ! initialize v (randomize)
   do m = 1,neig 
     if (m==1 .and. user_init) then 
       v(:,m) = sqrt(rho(:,isp)+vw_reg)
     else 
       do ii=1,nfft
         v(ii,m) = rand()
       enddo
     endif
     v(:,m) = v(:,m)/sqrt(sum(v(:,m)**2)*dvol)
     ! orthorgalize v(m) to v(1:m-1)
     do ib=1,m-1
       coeff = sum(v(:,m)*v(:,ib))*dvol
       v(:,m) = v(:,m) - coeff*v(:,ib)
     enddo
     v(:,m) = v(:,m)/sqrt(sum(v(:,m)**2)*dvol)
     !!write(logOF,*)'|v|',dvol*sum(v(:,m)**2)
   enddo


   !==================================
   ! mixing density 
   do k=1,niter 

     ! check norm of v
     do m=1,neig
       dtmp = dvol*sum(v(:,m)**2)
       if ( abs(dtmp-1.0d0)>1e-6 ) then
         write(logOF,*)'for m:',m,' norm of v is not 1, it is',dtmp
         stop
       endif 
     enddo

     !=====================================
     ! solve the eigenvectors  
     ! Outloop: diagonalize the subspace 
     do iter=1,ndiag 
         
       write(logof,'(a,i3,a)') '  ==> optimizing all eigenvectors. iter: ',iter,' <=='
       cg_tol = 1e-5
            
       !=======================================================
       ! conjugate gradient for optimizing each eigenvector
       do m = 1,neig

          write(logOF,'(a,i3,a)')NEW_LINE('a')//'optimizing eigenvector #',m; 
  
          ! compute kinetic energy 
          phi_tmp = sqrt(rho(:,isp)+vw_reg)
          phi_tmp = phi_tmp/sqrt(sum(dvol*phi_tmp**2))
          call laplacian(nfft,n1,n2,n3,phi_tmp,qvec,lap)
          ke = -sum(dvol*phi_tmp*lap)/2.d0
          write(logOF,'(a,f12.4)')'ke for this band: ',ke
          write(logOF,*)'qq_max: ',maxval(sum(qvec**2,1))

          do_precond = .true.

          !===============================
          ! cg performs this eigenvector 
          theta = 0.0d0
          do j=1,nline   

             call compute_Mx(nfft,rho(:,isp),vext(:,isp),v(:,m),Ax)
             lam = sum(v(:,m)*Ax)*dvol 
             g = Ax - lam*v(:,m)
             dconv = sqrt(sum((Ax-lam*v(:,m))**2)*dvol)

             ! make gradient orthogonal to all v
             do ib=1,neig
               coeff = sum(g*v(:,ib))*dvol
               g = g - coeff*v(:,ib)
             enddo

             ! precondition g => pg
             !================================
             if (do_precond) then 
               call precond_teter(n1,n2,n3,nfft,qvec,ke,g,pg)

               ! make pg orthorgnal to all bands 
               do ib=1,neig
                 coeff = sum(pg*v(:,ib))*dvol
                 pg = pg - coeff*v(:,ib)
               enddo
             endif 

             ! make conjugate gradient dir
             !================================
             if (j==1) then
               if (do_precond) d = -pg
               if (.not. do_precond) d = -g
             else
               if (do_precond) then 
                 cg_beta = sum(pg*g)/sum(pg_old*g_old)
                 d = -pg + cg_beta*d
               else 
                 cg_beta = sum(g*g)/sum(g_old*g_old)
                 d = -g + cg_beta*d
               endif 
             endif 

             ! make cg direction orth to all v
             !================================
             d_proj = d
             do ib=1,neig
               coeff = sum(d_proj*v(:,ib))*dvol
               d_proj = d_proj - coeff*v(:,ib)
             enddo
             ! normalize d
             d_proj = d_proj/sqrt(sum(d_proj*d_proj*dvol))

             ! exit? 
             if ( dconv<cg_tol .and. j>=3 ) then 
               write(logOF,'(a,es12.4,a)')'|Ax-lambda*x|<',cg_tol,', done cg'
               exit
             endif
             write(logOF,'(a,i4,a,f13.7,a,es9.2,a,es10.2,a,L2)') & 
              '(cg) iter: ',j,' <v,Av>:',lam,' |Ax-lam*x|:',dconv, &
              ' theta: ',theta,'   precond: ',do_precond
             call flush(logOF)

             g_old  = g
             pg_old = pg

!! DEBUG 
!             do ii=1,1000
!               theta = 0.001*(ii-1)
!               v0 = v(:,m)*cos(theta) + sin(theta)*d
!               call compute_Mx(nfft,rho(:,isp),vext(:,isp),v0,Ad)
!               write(logOF,*)'  -- theta: ',theta,sum(v0*Ad)*dvol
!             enddo
!! END OF DEBUG

             call compute_Mx(nfft,rho(:,isp),vext(:,isp),d_proj,Ad)
             theta = 0.5d0*atan(2.d0*sum(d_proj*Ax*dvol)/(lam-sum(d_proj*Ad*dvol)))

             ! theta in the window [0,pi/2] is the solution 
             ! if theta is < 0, wrap it to the [0,pi/2]
             if (theta<0.0d0) theta = theta + 3.1415926d0/2.d0

             v(:,m) = v(:,m)*cos(theta) + sin(theta)*d_proj
             v(:,m) = v(:,m)/sqrt(sum(v(:,m)**2*dvol)) ! enforce normalization 

             Ax = cos(theta)*Ax + sin(theta)*Ad
           enddo ! nline

           Ax_hist(:,m) = Ax

        enddo ! loop over all eigenvalues


        ! ensure that vectors are orthorgal to each other
        ! (even this step is not needed in principle) 
        do ib=1,neig
          do ib2=1,ib
            coeff = sum(v(:,ib)*v(:,ib2))*dvol
            write(logOF,'(es14.5)',advance='no')coeff
        !    v(:,ib) = v(:,ib) - coeff*v(:,ib2)
          enddo
          write(logOF,'(a)',advance='yes')''
        !  v(:,ib) = v(:,ib)/sqrt(sum(v(:,ib)**2)*dvol)
        enddo


        !
        ! Rayleigh-Ritz on span{X}
        ! https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Ritz_method
        !
        write(logOF,'(a)') NEW_LINE('a')//'Rayleigh-Ritz, making R matrix ...'
        do m=1,neig
          do j=1,neig
            R(m,j) = sum(v(:,j)*Ax_hist(:,m))*dvol
            write(logOF,'(f14.5)',advance='no')R(m,j)
          enddo
          write(logOF,'(a)',advance='yes')''
        enddo
        ! symmetrize R matrix (R in fact is symmetric in principle)
        R = (R + transpose(R))/2.d0
     
        !call check_R_symmetry(neig,R)
     
        ! diagonlize matrix R
        call DSYEV('V','U',neig,R,neig,lapack_w,lapack_WORK,3*neig-1,lapack_INFO )
     
        ! get new eigenvectors 
        v = MATMUL(v,R)
     
        write(logOF,'(a)')'        eigenvalues   '
        do m=1,neig
          a(m) = lapack_w(m)
          write(logOF,'(i4,f14.6)') m,a(m) 
        enddo
        write(logOF,*)'lowest eigenvalue: ',a(1)
     
        ! ensure that vectors are orthorgal to each other
        ! (even this step is not needed in principle) 
        do ib=1,neig
          do ib2=1,ib-1
            coeff = sum(v(:,ib)*v(:,ib2))*dvol
            v(:,ib) = v(:,ib) - coeff*v(:,ib2)
          enddo
          v(:,ib) = v(:,ib)/sqrt(sum(v(:,ib)**2)*dvol)
        enddo

      enddo ! iter 

      eigenvalues = a
      stop


      !============================
      ! mix density 

   enddo  ! density mixing 

 enddo ! spin


 write(logOF,'(a)')'leave OF_solver_scf() '//NEW_LINE('a')
 call flush(logOF)





 contains 

 !
 ! wf_reg is the solution of the equation
 ! -1/2*nabla^2 |wf_reg> + vscf|wf_reg> = mu|wf_reg>
 ! mu = <wf_ref,H wf_reg>/<wf_reg,wf_reg>
 !
 subroutine compute_Mx(n,rho,vext,x,Mx)
   implicit none 
   integer :: n
   real(8) :: rho(n), vext(n), & 
              x(n),ke,Mx(n) 

   ! local vars 
   real(8) :: lap(n), vscf(n), vtf(n)

   vtf = cTF*5.d0/3.d0*rho**(2.d0/3.d0)
   vscf = (vext + vtf)/vw_lam   
   call laplacian(nfft,n1,n2,n3,x,qvec,lap)
   Mx = -0.5d0*lap + vscf*x
   !mu = vw_lam*dvol*sum(wf_reg*arr_tmp)/(dvol*sum(wf_reg**2))
   !mu = mu/vw_lam

 end subroutine compute_Mx



end subroutine  OF_cg_scf

