
subroutine make_weights_recip(file_dfet_out,natom,nspin,nfft,nsys,melement,atom_info, & 
                              cell_nfft,xcart,cell_acell,znucl,rprimd,cov_radius,wsubsys)

 use mpi

 implicit none 

 integer :: nfft, nspin, nsys, natom, & 
            melement, file_dfet_out, cell_nfft(3)

 real(8) :: cov_radius(melement), & 
            wsubsys(nfft,nsys), & 
            xcart(3,natom), & 
            cell_acell(3,3), & 
            znucl(natom) , & 
            rprimd(3,3), & 
            atom_info(5,natom)

 ! local vars 
 real(8),allocatable :: wgauss(:)
 character(len=300) :: filename 
 integer :: nimag1,nimag2,nimag3, & 
            iatom,counter, & 
            ix,iy,iz,i,j, & 
            imag1,imag2,imag3, & 
            n1,n2,n3, myrank, ierr
 real(8) :: distance, atom_x, atom_y, atom_z, gau_width, expo, xtmp(3)


 complex(kind=8) :: fft1(cell_nfft(1)/2+1,cell_nfft(2),cell_nfft(3))
 complex(kind=8) :: fft2(cell_nfft(1)/2+1,cell_nfft(2),cell_nfft(3))

 call mpi_comm_rank ( mpi_comm_world, myrank, ierr )

 n1 =  cell_nfft(1)
 n2 =  cell_nfft(2)
 n3 =  cell_nfft(3)
 !
 ! define the guassian width as the covalent radius
 !
 allocate(wgauss(nfft))
 if (myrank==0) then 
    call system("rm  wsubsys_*xsf")
 endif
 wsubsys = 0.d0
 !
 write(file_dfet_out,'(a,a)')''
 write(file_dfet_out,'(a,a)')'making subsystem weight...'
 !
 ! nimag should be determined automatically by the code based on acell 
 ! and cov_radius
 !
 nimag1 = 1
 nimag2 = 1
 nimag3 = 1
 write(file_dfet_out,*)'WARNING: NIMAG1: 1, might not be large enough'
 write(file_dfet_out,*)'WARNING: NIMAG2: 1, might not be large enough'
 write(file_dfet_out,*)'WARNING: NIMAG3: 1, might not be large enough'
 ! 
 ! first, we make atom-centered gaussian weights 
 !
 do iatom=1,natom
   !!!print *,'atom: ',iatom
   wgauss =0.d0
   write(6,'(a,i3,a,i3)')'doing weight for atom:',iatom,' / ',natom
   call flush(6)
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
               xtmp(1) = rprimd(1,1)/dble(n1)*dble(ix-1)
               xtmp(2) = rprimd(2,2)/dble(n2)*dble(iy-1)
               xtmp(3) = rprimd(3,3)/dble(n3)*dble(iz-1)
               ! distance is in bohr
               distance = sqrt( (xtmp(1)-atom_x)**2 + &
                                (xtmp(2)-atom_y)**2 + &
                                (xtmp(3)-atom_z)**2  )
               !
               gau_width = cov_radius(int(atom_info(2,iatom)))/0.529177d0
               gau_width = gau_width/2.d0
              ! gau_width = gau_width/3.d0
               !
               expo = -0.5d0 * ( distance/gau_width )**2
               if (expo<-30.d0) then
                 wgauss(counter) = wgauss(counter) + exp(-30.d0)
               else
                 wgauss(counter) = wgauss(counter) + exp( expo )
               endif
             enddo
           enddo
         enddo
         !------ loop over points -----
       enddo
     enddo
   enddo
   !
   ! ====== make subsystem weights =====
   !
   do j=1,nsys
     ! check if the atom is in the subsystem 
     if (atom_info(1,iatom)==j) then 
       !print *,'add atom',iatom,' to subsystem weight ',j
       wsubsys(:,j) = wsubsys(:,j) + wgauss
     endif
   enddo
   !
 enddo ! loop over atoms
 deallocate(wgauss)
 write(file_dfet_out,'(a,a)')'done making subsystem weight.'

 !
 ! dump subsystem weights for plot
 !
 do j=1,nsys
   write(filename,'(i4)') j
   filename = "wsubsys_"//trim(adjustl(filename))//".xsf"
   call wrt_xsf(nfft,cell_nfft,natom,xcart,znucl,rprimd,filename,wsubsys(:,j))
   write(file_dfet_out,'(a,a)')trim(filename),' dumped.'
 enddo
 write(file_dfet_out,'(a,a)')''
 
end subroutine 
