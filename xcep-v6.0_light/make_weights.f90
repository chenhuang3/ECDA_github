
subroutine make_weights(file_dfet_out,natom,nspin,nfft,nsys,melement,atom_info, & 
                        cell_nfft,xcart,cell_acell,znucl,rprimd,cov_radius,watom)

 use mpi
! use interface_funcs

 implicit none 

 integer :: shape_func=3  ! 1: Gaussian function 
                          ! 2: Bell function defined as 1/( 1 + ((x-R)/a)^2b  )
                          ! 3: Becke weights

 integer :: nfft, nspin, nsys, natom,  ib, subsys, & 
            melement, file_dfet_out, cell_nfft(3)

 real(8) :: cov_radius(melement), & 
            wsubsys(nfft,nsys), & 
            watom(nfft,nsys), & 
            xcart(3,natom), r_scale,  rc, & 
            cell_acell(3,3), & 
            znucl(natom) , & 
            rprimd(3,3), & 
            atom_info(5,natom)

 ! local vars 
 real(8),allocatable :: wgauss(:)
 character(len=300) :: filename, fname, msg, str2
 integer :: nimag1,nimag2,nimag3, & 
            iatom,counter, & 
            ix,iy,iz,i,j, & 
            imag1,imag2,imag3, & 
            n1,n2,n3
 real(8) :: distance, atom_x, atom_y, atom_z, gau_width, expo, xtmp(3)
 integer :: myrank, ierr

 real(8) :: aa,bb,f(3), mu, & 
            atom_x2, atom_y2, atom_z2, dist2, ff,ss, dist3


 call mpi_comm_rank ( mpi_comm_world, myrank, ierr )

 n1 =  cell_nfft(1)
 n2 =  cell_nfft(2)
 n3 =  cell_nfft(3)
 !
 ! define the guassian width as the covalent radius
 !
 allocate(wgauss(nfft))
 if (myrank==0) call system("rm  wsubsys_*xsf")
 wsubsys = 0.d0
 !
 if (myrank==0) then 
   write(*,'(a,a)')''
   write(*,'(a,a)')'making subsystem weight...'
 endif
 !
 ! nimag should be determined automatically by the code based on acell 
 ! and cov_radius
 !
 nimag1 = 0
 nimag2 = 0
 nimag3 = 0

! if (myrank==0) then 
!   write(file_dfet_out,*)'WARNING: NIMAG1: 0, might not be large enough'
!   write(file_dfet_out,*)'WARNING: NIMAG2: 0, might not be large enough'
!   write(file_dfet_out,*)'WARNING: NIMAG3: 0, might not be large enough'
! endif

 ! 
 ! first, we make atom-centered gaussian weights 
 !
 do iatom=1,natom
   !!!print *,'atom: ',iatom
   wgauss = 1.0d0
   if (myrank==0) write(6,'(a,i3,a,i3)')'doing weight for atom:',iatom,' / ',natom
   call flush(6)
   !
   ! set the wgauss for each point in the cell
   do imag1 = -nimag1,  nimag1
     do imag2 = -nimag2,  nimag2
       do imag3 = -nimag3,  nimag3
         !
         atom_x = xcart(1,iatom) + & 
                  imag1*cell_acell(1,1) + imag2*cell_acell(2,1) + imag3*cell_acell(3,1)
         atom_y = xcart(2,iatom) + & 
                  imag1*cell_acell(1,2) + imag2*cell_acell(2,2) + imag3*cell_acell(3,2)
         atom_z = xcart(3,iatom) + & 
                  imag1*cell_acell(1,3) + imag2*cell_acell(2,3) + imag3*cell_acell(3,3)
         !
         ! ------- loop over points ----
         !
         counter = 0
         do iz =1,n3
           do iy=1,n2
             do ix=1,n1

               counter = counter + 1

               ! make the cooridate of the points 
               ! fractional coord.
               f(1) = dble(ix-1)/dble(n1)
               f(2) = dble(iy-1)/dble(n2)
               f(3) = dble(iz-1)/dble(n3)

               ! lattice constants are stored in col-wise way [a,b,c]
               ! cartesion coord in bohr 
               xtmp = matmul(rprimd,f)

               ! distance is in bohr
               distance = sqrt( (xtmp(1)-atom_x)**2 + &
                                (xtmp(2)-atom_y)**2 + &
                                (xtmp(3)-atom_z)**2  )

               if (shape_func==1) then 
                 !
                 gau_width = cov_radius(int(atom_info(2,iatom)))/0.5291772106d0
                 gau_width = gau_width/r_scale
                 !
                 expo = -0.5d0 * ( distance/gau_width )**2
                 if (exp(expo) /= exp(expo) ) then 
                   print *,'NaN for exp(expo) in make_weight.f90 file stop'
                   stop
                 endif 
                 if (expo<-60.d0) then
                   wgauss(counter) = wgauss(counter) + 0.0d0
                 else
                   wgauss(counter) = wgauss(counter) + exp( expo )
                 endif
               else if (shape_func==2) then 
                 !
                 ! bell function 
                 ! see
                 ! http://www.bindichen.co.uk/post/Fundamentals/bell-shaped-function.html
                 !
                 bb = 4.0
                 rc = cov_radius(int(atom_info(2,iatom)))
                 rc = rc/0.5291772106d0 ! convert to bohr 
                 aa = rc 
                 wgauss(counter) = wgauss(counter) + 1.d0/(1.d0+(distance/aa)**(2*bb))

               else if (shape_func==3) then 
                 ! becke's method   

                 ! loop over other atomm 
                 ! constructing S_AB function at the r point 
                 do ib=1,natom 
                   if (ib==iatom) cycle 
                   !
                   atom_x2 = xcart(1,ib)
                   atom_y2 = xcart(2,ib)
                   atom_z2 = xcart(3,ib) 

                   dist2 = sqrt( (atom_x2-xtmp(1))**2 + & 
                                 (atom_y2-xtmp(2))**2 + & 
                                 (atom_z2-xtmp(3))**2 ) 

                   dist3 = sqrt(  (atom_x2 - atom_x)**2 + & 
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
                   wgauss(counter) = ss*wgauss(counter)
                 enddo
               endif ! shape_function 

             enddo
           enddo
         enddo
         !------ loop over points -----
        ! print *,'min/max wgauss: ',minval(wgauss),maxval(wgauss)
       enddo
     enddo
   enddo

   !
   ! ====== make atom weights =====
   !
   subsys = atom_info(1,iatom)
   wsubsys(:,subsys) = wgauss
   !
 enddo ! loop over atoms
 deallocate(wgauss)

 if (myrank==0) write(*,'(a,a)')'done making weights for atoms.'



 ! normalize atom weights  
 do j=1,natom 
   watom(:,j) = wsubsys(:,j)/sum(wsubsys,2)
 enddo 


 !
 ! dump subsystem weights for plot
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

!     ! write the XSF files 
!     write(filename,'(i4)') j
!     filename = "atom_weight_"//trim(adjustl(filename))//".xsf"
!     print *,'writing file: ',trim(filename),' ...'
!     call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,watom(:,j))
!     write(file_dfet_out,'(a,a)')trim(filename),' dumped.'
   enddo

   write(*,*)''
 endif 
! call flush(6)
!
 call mpi_barrier(mpi_comm_world,ierr)

contains 

 function  becke_func(x)
   real(8) :: becke_func
   real(8) :: x 
   becke_func = (x/2.d0)*(3.0d0-x*x)
 end function 


end subroutine 






