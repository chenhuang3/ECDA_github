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

subroutine wrt_xsf_geo(nfft, cell_nfft,& 
                   natom_dump, & 
                   xcart_dump, & 
                   znucl_dump, & 
                   rprimd,fil)
  
  implicit none 
  
  character(len=*) :: fil
  integer :: jj, kk, it, mu, iatom, nfft, & 
             natom_dump,cell_nfft(3)  ! 1: density, 2: potential

  real(8) :: xcart_dump(3,natom_dump), & 
             rprimd(3,3), & 
             dat(nfft), & 
             znucl_dump(natom_dump)

  !do iatom=1,dtset%natom
  !  do mu=1,3
  !    xcart(mu,iatom)=rprimd(mu,1)*xred(1,iatom) & 
  !        +rprimd(mu,2)*xred(2,iatom) & 
  !        +rprimd(mu,3)*xred(3,iatom)
  !  end do
  !end do

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
  close(222)

end subroutine wrt_xsf_geo


