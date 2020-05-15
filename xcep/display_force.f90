subroutine display_force(logf,natom,force_xc,force_xc1,force_xc2,force_xcN,force_pen,force_abinit,force_total)

  implicit none 

  integer :: natom,ia,logf
  real(8) :: force_xc(3,natom)
  real(8) :: force_xc1(3,natom)
  real(8) :: force_xc2(3,natom)
  real(8) :: force_xcN(3,natom)
  real(8) :: force_abinit(3,natom)
  real(8) :: force_pen(3,natom)
  real(8) :: force_total(3,natom)

  ! display force 
  write(logf,*)''
  write(logf,'(a)')'total forces: '
  do ia=1,natom 
    write(logf,'(a,i4,3f14.6)')'  atom: ',ia,force_total(1:3,ia)
  enddo
  write(logf,*)''
  write(logf,'(a,es12.4)')'max(abs(total force)): ',maxval(abs(force_total))
  write(logf,*)''
  write(logf,'(a)')'ion-ion + e-ion forces from ABINIT:'
  do ia=1,natom 
    write(logf,'(a,i4,3f14.6)')'  atom: ',ia,force_abinit(1:3,ia)
  enddo 
  write(logf,*)''
  write(logf,'(a)')'non-Hellmann-Feynman XC force:'
  do ia=1,natom 
    write(logf,'(a,i4,3f14.6)')'  atom: ',ia,force_xc(1:3,ia)
  enddo
  write(logf,'(a)')'non-Hellmann-Feynman penalty force:'
  do ia=1,natom 
    write(logf,'(a,i4,3f14.6)')'  atom: ',ia,force_pen(1:3,ia)
  enddo


!=================================
! details of the forces 
!=================================

if (.false.) then 
  write(logf,*)''
  write(logf,'(a)')'--- decomposition of the non-Hellmann-Feynman XC force) --- '
  write(logf,*)'[force_xc1 (due to change of atomic weights)]:'
  do ia=1,natom 
    write(logf,'(a,i4,3f14.6)')'atom: ',ia,force_xc1(1:3,ia)
  enddo 
  write(logf,'(a)')'[force_xc2 (due to change of cluster KS potential)]:'
  do ia=1,natom 
    write(logf,'(a,i4,3f14.6)')'atom: ',ia,force_xc2(1:3,ia)
  enddo 
  write(logf,'(a)')'[force_xcN (due to change of cluster electron number)]:'
  do ia=1,natom 
    write(logf,'(a,i4,3f14.6)')'atom: ',ia,force_xcN(1:3,ia)
  enddo 
endif   



  write(logf,*)''
  call flush(logf)

end subroutine  
