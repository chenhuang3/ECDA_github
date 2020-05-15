!
! convert data cell_data -> atom_data 
! they are on different meshes 
!
subroutine cell2atom(s,cell_data,cell_nfft,cell_acell, & 
                     atbox_acell,atbox_npt,xatom,atom_data)
  
  use bspline 

  implicit none 

  integer :: s, &            ! atom index
             cell_nfft(3), &
             atbox_npt       ! do not double count boundary points

  real(8) :: cell_data(cell_nfft(1),cell_nfft(2),cell_nfft(3)), & 
             cell_acell(3), & 
             atbox_acell, & 
             atom_data(atbox_npt,atbox_npt,atbox_npt), & 
             xatom(3) ! coordinate of atom in the cell

  ! internal vars 

  integer :: q,j1,j2,j3,i1,i2,i3,xpo,ypo,zpo
  real(8) :: pt(3),space(3),xd,yd,zd, atSpace, & 
             c00,c01,c10,c11,c0,c1

  ! for interpolation
  integer,parameter :: kxord=4,kyord=4,kzord=4  ! third order piecewise polynomial
  real(8) :: xx,yy,zz,ret 
  integer :: nx,ny,nz,ldf, mdf
  real(8) :: xknot(kxord+cell_nfft(1)+1),yknot(kzord+cell_nfft(2)+1),zknot(kyord+cell_nfft(3)+1)
  real(8) :: BCOEF((cell_nfft(1)+1)*(cell_nfft(2)+1)*(cell_nfft(3)+1))              
  real(8) :: xvec(cell_nfft(1)+1),yvec(cell_nfft(2)+1),zvec(cell_nfft(3)+1)
  real(8) :: xyzdata( cell_nfft(1)+1,cell_nfft(2)+1,cell_nfft(3)+1 )

  ! >>>>>>>>>>> function begins <<<<<<<<<<!

  print *,'enter cell2atom().'
  write(6,'(a,3f12.4)')'xatom: ',xatom

  space = cell_acell/cell_nfft

  ! subsystem performs aperiodic calculations
  ! the total number of points are atbox_npt,which include the 
  ! all the points in boundaries.

  atSpace = atbox_acell/dble(atbox_npt-1)

  !---------- 3d interpolation -------------
  !
  ! example can be found at
  ! http://docs.roguewave.com/imsl/fortran/6.0/math/default.htm?turl=bs3in.htm
  ! a good ref: https://www.cs.unc.edu/~dm/UNC/COMP258/Papers/bsplbasic.pdf
  !
  nx = cell_nfft(1)+1
  ny = cell_nfft(2)+1
  nz = cell_nfft(3)+1
  do j1 = 1,nx
    xvec(j1) = j1 
  enddo
  do j1 = 1,ny
    yvec(j1) = j1
  enddo
  do j1=1,nz
    zvec(j1) = j1 
  enddo 
  ! make xyzdata, consider periodic boundary condition
  do j1=1,nx
    do j2=1,ny
      do j3=1,nz
        xpo=j1
        ypo=j2
        zpo=j3
        if (j1==cell_nfft(1)+1) xpo=1
        if (j2==cell_nfft(2)+1) ypo=1
        if (j3==cell_nfft(3)+1) zpo=1
        xyzdata(j1,j2,j3)=cell_data(xpo,ypo,zpo)
      enddo
    enddo
  enddo
  !  generate knots
  call dbsnak (nx, xvec, kxord, xknot)
  call dbsnak (ny, yvec, kyord, yknot)
  call dbsnak (nz, zvec, kzord, zknot)
  ldf = nx 
  mdf = ny 
  call dbs3in(nx,xvec,ny,yvec,nz,zvec,xyzdata,ldf,mdf,kxord,kyord,kzord,xknot,yknot,zknot, & 
              bcoef)  ! output
              

  ! loop over atom box 
  do j1=1,atbox_npt
    do j2=1,atbox_npt
      do j3=1,atbox_npt

        ! In Octopus, the atom is centered in the box, 
        ! atbox_npt double counts the boundary points
        ! Find the coordinate of the mesh point j1,j2,j3 

        pt(1) = - atbox_acell/2.d0 + atSpace*(j1-1)
        pt(2) = - atbox_acell/2.d0 + atSpace*(j2-1)
        pt(3) = - atbox_acell/2.d0 + atSpace*(j3-1)

        ! adjust pt using atom coordinate
        pt = pt + xatom

        ! wrap pt into the box
        do q=1,3
          do while (pt(q)<0.d0)
            pt(q) = pt(q) + cell_acell(q)
          enddo
          do while (pt(q)>cell_acell(q))
            pt(q) = pt(q) - cell_acell(q)
          enddo
        enddo

        ! locate the small mesh of the cell that enclose this point
        xx = pt(1)/space(1)+1
        yy = pt(2)/space(2)+1
        zz = pt(3)/space(3)+1
        ! cubic bspline 
        atom_data(j1,j2,j3) = dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nx,ny,nz,bcoef)
      enddo
    enddo
  enddo

  print *,'left cell2atom().'

end subroutine cell2atom
