!
! compute the derivative of atom j's weight w.r.t. the coordinate of atom i (R_i)
! using finite difference method 
!
! output: 
!   dwatom_dR: derivative of atom weights with respect to atom iR 's coordinate
!   dwcluster_dR: derivative of clusters' weights w.r.t. atom iR's cooridates
!
subroutine calc_dWatom_dR(logf,natom,nspin,nfft,nsys,melement,atom_info, & 
                          cell_nfft,xcart,cell_acell,znucl,rprimd,cov_radius, & 
                          mlist,nlist,ilist,iR,dwatom_dR,dwcluster_dR)
               
implicit none 

integer :: iR, &   ! index of the atom to move 
           nfft, & 
           nspin, nsys, natom, logf, & 
           melement, cell_nfft(3), & 
           mlist, nlist(natom), & 
           ilist(mlist,natom)

real(8) :: dwatom_dR(3,nfft,natom), & 
           dwcluster_dR(3,nfft,natom), &
           cov_radius(melement), & 
           watom1(nfft,nsys), & 
           watom2(nfft,nsys), & 
           xcart(3,natom), & 
           cell_acell(3,3), & 
           znucl(natom) , & 
           rprimd(3,3), & 
           atom_info(5,natom), & 
           xcart_new(3,natom), & 
           dstep 

! local vars 
integer  :: i, method, ia, kk, k


!>> function begin <<!

method = 1   ! 1: finite difference 
             ! 2: analytical derivative  


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! d(atom_weight)/dR for all the atoms
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (method==1) then 
  dstep = 0.001d0

  ! X, Y, Z directions  
  do i = 1,3
    ! (+) 
    xcart_new = xcart  
    xcart_new(i,iR) = xcart(i,iR) + dstep 
    call make_weights(logf,natom,nspin,nfft,nsys,melement,atom_info,& 
                      cell_nfft,xcart_new,cell_acell,znucl,rprimd,cov_radius,watom1)
    ! (-)
    xcart_new = xcart 
    xcart_new(i,iR) = xcart(i,iR) - dstep 
    call make_weights(logf,natom,nspin,nfft,nsys,melement,atom_info,& 
                      cell_nfft,xcart_new,cell_acell,znucl,rprimd,cov_radius,watom2)

    dwatom_dR(i,:,:) = (watom1-watom2)/2.d0/dstep
  enddo 
endif 



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! get d(cluster_weight) / d(R) for all clusters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
do ia = 1,natom  ! loop over clusters
  dwcluster_dR(:,:,ia) = dwatom_dR(:,:,ia)
  ! all buffer atoms 
  do kk=1,nlist(ia)
    dwcluster_dR(:,:,ia) = dwcluster_dR(:,:,ia) + dwatom_dR(:,:,ilist(kk,ia))
  enddo       
enddo


end subroutine calc_dwatom_dR
