
subroutine make_weights_pbc(file_dfet_out,natom,nspin,nfft,nsys,melement,atom_info, & 
                        cell_nfft,xcart,cell_acell,znucl,rprimd,cov_radius,watom)

 use mpi
! use interface_funcs

 implicit none 

 integer :: shape_func=3  ! 1: Gaussian function 
                          ! 2: Bell function defined as 1/( 1 + ((x-R)/a)^2b  )
                          ! 3: Becke weights

 integer :: nfft, nspin, nsys, natom,  ib, & 
            melement, file_dfet_out, cell_nfft(3)

 real(8) :: cov_radius(melement), & 
            watom_tmp(nfft,nsys), & 
            watom(nfft,nsys), & 
            xcart(3,natom), r_scale,  rc, & 
            cell_acell(3,3), & 
            znucl(natom) , & 
            rprimd(3,3), & 
            atom_info(5,natom)

 ! local vars 
 logical :: wgauss_too_small
 real(8) :: wgauss
 character(len=300) :: filename, fname, msg, str2
 integer :: nimag1,nimag2,nimag3, & 
            iatom,ipt, im1, im2, im3, & 
            ix,iy,iz,i,j, & 
            imag1,imag2,imag3, & 
            n1,n2,n3
 real(8) :: distance, atom_x, atom_y, atom_z, gau_width, expo, xtmp(3)
 integer :: myrank, ierr, span

 real(8) :: aa,bb,f(3), mu, & 
            atom_x2, atom_y2, atom_z2, dist2, ff,ss, dist3


 call mpi_comm_rank ( mpi_comm_world, myrank, ierr )

 n1 =  cell_nfft(1)
 n2 =  cell_nfft(2)
 n3 =  cell_nfft(3)

 ! define the guassian width as the covalent radius
 if (myrank==0) call system("rm  wsubsys_*xsf")
 watom_tmp = 0.d0


 if (myrank==0) then 
   write(*,'(a,a)')''
   write(*,'(a,a)')'making subsystem weight (PBC) in parallel ...'
 endif
 !
 ! nimag should be determined automatically by the code based on acell 
 ! and cov_radius
 !
 nimag1 = 1
 nimag2 = 1
 nimag3 = 1

 span = nfft/10.d0;

 ! ------- loop over points in the cell ----
 ipt = 0
 do iz =1,n3
   do iy=1,n2
     do ix=1,n1

       ipt = ipt + 1

       ! show progress .................
       if ( mod(ipt,span) == 0  .and.  myrank==0) then
         write(6,'(a,f6.1,a)')'processor 0: ', dble(ipt/span)/10.d0*100.d0,' %'; call flush(6)
       endif 

       ! make the cooridate of the points 
       ! fractional coord.
       f(1) = dble(ix-1)/dble(n1)
       f(2) = dble(iy-1)/dble(n2)
       f(3) = dble(iz-1)/dble(n3)

       ! lattice constants are stored in col-wise way [a,b,c]
       ! cartesion coord in bohr 
       xtmp = matmul(rprimd,f)
       
       ! ======= make weight for iatom (consider its images) =======
       ! distribute the job to each processor 
       !!!do iatom=1,natom

         iatom = myrank + 1 ! set iatom for this processor 

         do im1 = -1,1
           do im2 = -1,1
             do im3 = -1,1

               wgauss = 1.0d0

               atom_x = xcart(1,iatom) + im1*cell_acell(1,1) + im2*cell_acell(2,1) + im3*cell_acell(3,1)
               atom_y = xcart(2,iatom) + im1*cell_acell(1,2) + im2*cell_acell(2,2) + im3*cell_acell(3,2)
               atom_z = xcart(3,iatom) + im1*cell_acell(1,3) + im2*cell_acell(2,3) + im3*cell_acell(3,3)

               ! distance is in bohr
               distance = sqrt( (xtmp(1)-atom_x)**2 + &
                                (xtmp(2)-atom_y)**2 + &
                                (xtmp(3)-atom_z)**2  )

               ! becke's method   
               ! loop over other atomm in the system (we just loop over cells +/- 1 in x, y, z direction
               ! constructing S_AB function at the r point 
               wgauss_too_small = .false. 

               do ib=1,natom 
                 do imag1 = -1,1
                   do imag2 = -1,1
                     do imag3 = -1,1

                       if (ib==iatom .and. im1==imag1 .and. im2==imag2 .and. im3==imag3) cycle 

                       atom_x2 = xcart(1,ib) + imag1*cell_acell(1,1) + imag2*cell_acell(2,1) + imag3*cell_acell(3,1)
                       atom_y2 = xcart(2,ib) + imag1*cell_acell(1,2) + imag2*cell_acell(2,2) + imag3*cell_acell(3,2)
                       atom_z2 = xcart(3,ib) + imag1*cell_acell(1,3) + imag2*cell_acell(2,3) + imag3*cell_acell(3,3)

                       dist2 = sqrt( (atom_x2 - xtmp(1))**2 + & 
                                     (atom_y2 - xtmp(2))**2 + & 
                                     (atom_z2 - xtmp(3))**2 ) 
  
                       dist3 = sqrt( (atom_x2 - atom_x)**2 + & 
                                     (atom_y2 - atom_y)**2 + & 
                                     (atom_z2 - atom_z)**2 ) 

                       mu = (distance - dist2)/dist3
                       if ( abs(mu)>1.d0 + 1e-6) then 
                         print *,'mu cannot be larger than 1, error stop'
                         print *,'mu: ',mu
                         stop
                       endif 
                 
                       ff = becke_func(mu)
                       !! the less ff() the softer
                       !!  ff = becke_func(ff) 
                       !!  ff = becke_func(ff) 
                       ss = 0.5d0 * (1.d0 - ff)
                       ! ss is from [0,1]
                       wgauss = ss*wgauss
                     enddo 
                   enddo
                 enddo
               enddo
               !========== end of loop of second atom in system ==========
             
               watom_tmp(ipt,iatom) = watom_tmp(ipt,iatom) + wgauss

             enddo ! image 1
           enddo ! image 2
         enddo ! image 3

         ! exit for this processor 
       !!!enddo ! iatom
       !============== end of loop of first atom ============


     enddo ! loop over points 
   enddo ! loop over points 
 enddo ! loop over atoms


 call mpi_allreduce(MPI_IN_PLACE,watom_tmp,natom*nfft,mpi_double_precision,mpi_sum, mpi_comm_world, ierr)


 if (myrank==0) write(*,'(a,a)')'done making weights for atoms.'


 !
 ! normalize atom weights
 !
 do j=1,natom 
   watom(:,j) = watom_tmp(:,j)/sum(watom_tmp,2)
 enddo 


 !
 ! dump atomic weight functions for plot
 !
 if (myrank==0) then 
   print *,'writing atom_weight files ... '
   call flush(6)

   do j=1,natom
     write(filename,'(i4)') j
     str2='subsys'//trim(adjustl(filename))//'/'
     filename = trim(str2)//"atom_weight.dat"
     write(6,'(a)')'writing file: '//trim(filename)
     open(file=filename,unit=111,action='write',form='unformatted')
     write(111) watom(:,j) 
     close(111)

     ! write the XSF files 
     write(filename,'(i4)') j
     filename = "atom_weight_"//trim(adjustl(filename))//".xsf"
     print *,'writing file: ',trim(filename),' ...'
     call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,watom(:,j))
     write(file_dfet_out,'(a,a)')trim(filename),' dumped.'
   enddo

   write(*,*)''
 endif 
 call flush(6)
 call mpi_barrier(mpi_comm_world,ierr)

contains 

 function  becke_func(x)
   real(8) :: becke_func
   real(8) :: x 
   becke_func = (x/2.d0)*(3.0d0-x*x)
 end function 


end subroutine 






