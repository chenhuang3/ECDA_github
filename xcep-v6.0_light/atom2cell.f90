
!
! convert data cell_data -> atom_data 
! they are on different meshes 
!

subroutine atom2cell(s,cell_data,cell_nfft,cell_acell,ucvol,& 
                     atbox_acell,atbox_npt,xatom,atom_data)

  use bspline

  implicit none 

  integer :: s, &            ! atom index
             cell_nfft(3), &
             atbox_npt       ! do not double count boundary points

  real(8) :: cell_data(cell_nfft(1),cell_nfft(2),cell_nfft(3)), & 
             cell_acell(3), ucvol,&
             atbox_acell, & 
             atom_data(atbox_npt,atbox_npt,atbox_npt), & 
             xatom(3) ! coordinate of atom in the cell

  ! internal vars 

  logical :: out_range

  integer :: q,j1,j2,j3,i1,i2,i3,xpo,ypo,zpo,n1,n2,n3, &
             nimag1, nimag2, nimag3, ix, iy, iz, counter

  real(8) :: pt(3),space(3),xd,yd,zd, atSpace, & 
             c00,c01,c10,c11,c0,c1,offset(3),xatom_tmp(3)


  ! for interpolation
  real(8) :: xx,yy,zz,ret 
  integer,parameter :: kxord=4,kyord=4,kzord=4
  integer :: nx,ny,nz,ldf, mdf
  real(8) :: xknot(kxord+atbox_npt),yknot(kzord+atbox_npt),zknot(kyord+atbox_npt)
  real(8) :: BCOEF(atbox_npt*atbox_npt*atbox_npt)              
  real(8) :: xvec(atbox_npt),yvec(atbox_npt),zvec(atbox_npt)
  real(8) :: xyzdata(atbox_npt,atbox_npt,atbox_npt)


  ! >>>>>>>>>>> function begins <<<<<<<<<<!

  n1 = cell_nfft(1)
  n2 = cell_nfft(2)
  n3 = cell_nfft(3)

  print *,'enter atom2cell().'

  write(6,'(a,3f12.4)')'xatom: ',xatom

  space = cell_acell/dble(cell_nfft)

  atSpace = atbox_acell/dble(atbox_npt-1)

  cell_data = 0.d0

  ! how many images do we need? 
  nimag1 = int(atbox_acell/2.d0 / cell_acell(1))+2
  nimag2 = int(atbox_acell/2.d0 / cell_acell(2))+2
  nimag3 = int(atbox_acell/2.d0 / cell_acell(3))+2

  write(6,'(a,3i4)')'nimag: ',nimag1,nimag2,nimag3

  counter = 0
  
  !---------- 3d interpolation -------------
  !
  ! example can be found at
  ! http://docs.roguewave.com/imsl/fortran/6.0/math/default.htm?turl=bs3in.htm
  ! a good ref: https://www.cs.unc.edu/~dm/UNC/COMP258/Papers/bsplbasic.pdf
  !
  nx = atbox_npt 
  ny = atbox_npt
  nz = atbox_npt 
  do j1 = 1,nx 
    xvec(j1) = j1 
    yvec(j1) = j1 
    zvec(j1) = j1 
  enddo 
  xyzdata = atom_data
  !  generate knots
  call dbsnak (nx, xvec, kxord, xknot)
  call dbsnak (ny, yvec, kyord, yknot)
  call dbsnak (nz, zvec, kzord, zknot)
  ldf = nx 
  mdf = ny 
  call dbs3in(nx,xvec,ny,yvec,nz,zvec,xyzdata,ldf,mdf,kxord,kyord,kzord,xknot,yknot,zknot, & 
              bcoef)  ! output 


  ! loop over all points in the cell 
  do j1=1,n1
    do j2=1,n2
      do j3=1,n3

        ! the cooridnate 
        pt(1) = space(1)*dble(j1-1)
        pt(2) = space(2)*dble(j2-1)
        pt(3) = space(3)*dble(j3-1)
        
        ! loop over images 
        do ix = -nimag1,nimag1
          do iy= -nimag2,nimag2
            do iz= -nimag3,nimag3

              xatom_tmp(1) = xatom(1) + dble(ix)*cell_acell(1)
              xatom_tmp(2) = xatom(2) + dble(iy)*cell_acell(2)
              xatom_tmp(3) = xatom(3) + dble(iz)*cell_acell(3)
              
              ! test if pt is outside of atom box
              out_range = .false.
              if (  abs(pt(1)-xatom_tmp(1))>=atbox_acell/2.d0 .or. & 
                    abs(pt(2)-xatom_tmp(2))>=atbox_acell/2.d0 .or. & 
                    abs(pt(3)-xatom_tmp(3))>=atbox_acell/2.d0 ) then
                out_range = .true.
              endif

              ! in the atom box 
              if ( .not. out_range ) then 
                ! find the small element encloses pt 
                offset = pt + 1e-10 - (xatom_tmp - atbox_acell/2.d0)
                xx = offset(1)/atSpace+1
                yy = offset(2)/atSpace+1
                zz = offset(3)/atSpace+1
                ! cubic bspline
                ret = dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nx,ny,nz,bcoef)
                cell_data(j1,j2,j3) = ret
                counter = counter+ 1
              endif

            enddo
          enddo
        enddo
        ! end of loop over images 

      enddo  ! loop over cell points 
    enddo
  enddo

  write(6,'(a,g12.5,a)')'atom2cell: total_Q (on cell mesh): ', & 
  sum(cell_data)*ucvol/dble(n1*n2*n3),' [check it]'
  write(6,'(a,g12.5,a)')'atom2cell: total_Q (on atom mesh): ', & 
  sum(atom_data)*atbox_acell**3/dble((atbox_npt-1)**3),' [check it]'
  print *,'left atom2cell().'


end subroutine atom2cell
