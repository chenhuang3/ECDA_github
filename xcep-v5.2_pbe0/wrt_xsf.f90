!
! write .xsf files for VESTA and xcrysden
! input: dat(nfft) is the data to be dump
! dat_type=1 : density 
! dat_type=2 : potential 
! fil: file name to be dumped 
! natom_dump: number of atoms to be dumpped
! xcart_dump: cartesion coordinates of atoms to be dumpped
! znucl_dump: zions of the atoms to be dumped
!

subroutine wrt_xsf(nfft, cell_nfft,& 
                   natom_dump, & 
                   xcart_dump, & 
                   znucl_dump, & 
                   rprimd,fil,dat)
  
  implicit none 
  
  character(len=*) :: fil
  integer :: n1,n2,n3, i_eff, j_eff, k_eff, & 
             ipt, ix,iy,iz, & 
             i,j,k,jj, kk, it, mu, iatom, nfft, & 
             natom_dump,cell_nfft(3)  ! 1: density, 2: potential

  real(8) :: xcart_dump(3,natom_dump), & 
             rprimd(3,3), & 
             dat(cell_nfft(1),cell_nfft(2),cell_nfft(3)),& 
             dat_aug(cell_nfft(1)+1,cell_nfft(2)+1,cell_nfft(3)+1),&
             znucl_dump(natom_dump)

  !do iatom=1,dtset%natom
  !  do mu=1,3
  !    xcart(mu,iatom)=rprimd(mu,1)*xred(1,iatom) & 
  !        +rprimd(mu,2)*xred(2,iatom) & 
  !        +rprimd(mu,3)*xred(3,iatom)
  !  end do
  !end do

  n1 = cell_nfft(1)
  n2 = cell_nfft(2)
  n3 = cell_nfft(3)

  
  ! map dat -> dat_aug
  do i=1,n1+1
    do j=1,n2+1
      do k=1,n3+1
        i_eff = i
        j_eff = j 
        k_eff = k
        if (i==n1+1) i_eff = 1
        if (j==n2+1) j_eff = 1
        if (k==n3+1) k_eff = 1
        dat_aug(i,j,k) = dat(i_eff,j_eff,k_eff)
      enddo
    enddo
  enddo


  open(file=fil,unit=222,action='write',form='formatted')
  write(222,*)'DIM-GROUP'
  write(222,*)'3    1'
  write(222,*)'PRIMVEC'
  write(222,'(3f16.4)') rprimd(1:3,1)*0.529177
  write(222,'(3f16.4)') rprimd(1:3,2)*0.529177
  write(222,'(3f16.4)') rprimd(1:3,3)*0.529177
  write(222,*)'PRIMCOORD'
  write(222,'(a,i5,a)')'  ',natom_dump,'  1'
  do jj=1,natom_dump
    !it = dtset%typat(jj)
    write(222,'(a,i5,a,3f18.4)')'  ',int(znucl_dump(jj)),'  ',xcart_dump(1:3,jj)*0.529177
  enddo
  write(222,*)'ATOMS'
  do jj=1,natom_dump
    !it = dtset%typat(jj)
    write(222,'(a,i5,a,3f18.4)')'  ',int(znucl_dump(jj)),'  ',xcart_dump(1:3,jj)*0.529177
  enddo
  write(222,*)'BEGIN_BLOCK_DATAGRID3D'
  write(222,*)'datagrids'
  write(222,*)'DATAGRID_3D_DENSITY'
  
  ! Note: XSF uses general lattice which includes the boundaries in all direction 
  write(222,'(3i6)')cell_nfft(1)+1,cell_nfft(2)+1,cell_nfft(3)+1
  write(222,*)'0.0   0.0   0.0'
  write(222,'(3f16.4)') rprimd(1:3,1)*0.529177
  write(222,'(3f16.4)') rprimd(1:3,2)*0.529177
  write(222,'(3f16.4)') rprimd(1:3,3)*0.529177
  !==================
  ! data
  !==================
  ipt = 0
  do iz=1,n3+1
     do iy=1,n2+1
        do ix=1,n1+1 
          ipt = ipt + 1
          write(222,'(f20.8)',ADVANCE = "NO")dat_aug(ix,iy,iz)
          if (mod(ipt,6)==0 .and. ipt/=nfft) write(222,*)
        enddo 
     enddo 
  enddo
  write(222,*)
  write(222,*)'END_DATAGRID_3D'
  write(222,*)'END_BLOCK_DATAGRID3D'
  close(222)

end subroutine wrt_xsf


