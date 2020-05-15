!> @file
!!  Routines to manipulate the grid
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculates the overall size of the simulation cell 
!! and shifts the atoms such that their position is the most symmetric possible.
!! Assign these values to the global localisation region descriptor.
subroutine system_size(iproc,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Glr,shift)
   use module_base
   use module_types
   use yaml_output
   implicit none
   type(atoms_data), intent(inout) :: atoms
   integer, intent(in) :: iproc
   real(gp), intent(in) :: crmult,frmult
   real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
   real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
   real(gp), intent(inout) :: hx,hy,hz
   type(locreg_descriptors), intent(out) :: Glr
   real(gp), dimension(3), intent(out) :: shift
   !Local variables
   !character(len=*), parameter :: subname='system_size'
   integer, parameter :: lupfil=14
   real(gp), parameter ::eps_mach=1.e-12_gp
   integer :: iat,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
   real(gp) :: rad,cxmin,cxmax,cymin,cymax,czmin,czmax,alatrue1,alatrue2,alatrue3

   !check the geometry code with the grid spacings
   if (atoms%astruct%geocode == 'F' .and. (hx/=hy .or. hx/=hz .or. hy/=hz)) then
      write(*,'(1x,a)')'ERROR: The values of the grid spacings must be equal in the Free BC case'
      stop
   end if

   !calculate the extremes of the boxes taking into account the spheres around the atoms
   cxmax=-1.e10_gp 
   cxmin=1.e10_gp

   cymax=-1.e10_gp 
   cymin=1.e10_gp

   czmax=-1.e10_gp 
   czmin=1.e10_gp

   do iat=1,atoms%astruct%nat

      rad=radii_cf(atoms%astruct%iatype(iat),1)*crmult

      cxmax=max(cxmax,rxyz(1,iat)+rad) 
      cxmin=min(cxmin,rxyz(1,iat)-rad)

      cymax=max(cymax,rxyz(2,iat)+rad) 
      cymin=min(cymin,rxyz(2,iat)-rad)

      czmax=max(czmax,rxyz(3,iat)+rad) 
      czmin=min(czmin,rxyz(3,iat)-rad)
   enddo

   !eliminate epsilon form the grid size calculation
   !!  cxmax=cxmax+eps_mach 
   !!  cymax=cymax+eps_mach  
   !!  czmax=czmax+eps_mach  
   !!
   !!  cxmin=cxmin-eps_mach
   !!  cymin=cymin-eps_mach
   !!  czmin=czmin-eps_mach

   !define the box sizes for free BC, and calculate dimensions for the fine grid with ISF
   if (atoms%astruct%geocode == 'F') then
      atoms%astruct%cell_dim(1)=(cxmax-cxmin)
      atoms%astruct%cell_dim(2)=(cymax-cymin)
      atoms%astruct%cell_dim(3)=(czmax-czmin)

      ! grid sizes n1,n2,n3
      n1=int(atoms%astruct%cell_dim(1)/hx)
      !if (mod(n1,2)==1) n1=n1+1
      n2=int(atoms%astruct%cell_dim(2)/hy)
      !if (mod(n2,2)==1) n2=n2+1
      n3=int(atoms%astruct%cell_dim(3)/hz)
      !if (mod(n3,2)==1) n3=n3+1
      alatrue1=real(n1,gp)*hx
      alatrue2=real(n2,gp)*hy
      alatrue3=real(n3,gp)*hz

      n1i=2*n1+31
      n2i=2*n2+31
      n3i=2*n3+31

   else if (atoms%astruct%geocode == 'P') then 
      !define the grid spacings, controlling the FFT compatibility
      call correct_grid(atoms%astruct%cell_dim(1),hx,n1)
      call correct_grid(atoms%astruct%cell_dim(2),hy,n2)
      call correct_grid(atoms%astruct%cell_dim(3),hz,n3)
      alatrue1=(cxmax-cxmin)
      alatrue2=(cymax-cymin)
      alatrue3=(czmax-czmin)

      n1i=2*n1+2
      n2i=2*n2+2
      n3i=2*n3+2

   else if (atoms%astruct%geocode == 'S') then
      call correct_grid(atoms%astruct%cell_dim(1),hx,n1)
      atoms%astruct%cell_dim(2)=(cymax-cymin)
      call correct_grid(atoms%astruct%cell_dim(3),hz,n3)

      alatrue1=(cxmax-cxmin)
      n2=int(atoms%astruct%cell_dim(2)/hy)
      alatrue2=real(n2,gp)*hy
      alatrue3=(czmax-czmin)

      n1i=2*n1+2
      n2i=2*n2+31
      n3i=2*n3+2

   end if

   !balanced shift taking into account the missing space
   cxmin=cxmin+0.5_gp*(atoms%astruct%cell_dim(1)-alatrue1)
   cymin=cymin+0.5_gp*(atoms%astruct%cell_dim(2)-alatrue2)
   czmin=czmin+0.5_gp*(atoms%astruct%cell_dim(3)-alatrue3)

   !correct the box sizes for the isolated case
   if (atoms%astruct%geocode == 'F') then
      atoms%astruct%cell_dim(1)=alatrue1
      atoms%astruct%cell_dim(2)=alatrue2
      atoms%astruct%cell_dim(3)=alatrue3
   else if (atoms%astruct%geocode == 'S') then
      cxmin=0.0_gp
      atoms%astruct%cell_dim(2)=alatrue2
      czmin=0.0_gp
   else if (atoms%astruct%geocode == 'P') then
      !for the moment we do not put the shift, at the end it will be tested
      !here we should put the center of mass
      cxmin=0.0_gp
      cymin=0.0_gp
      czmin=0.0_gp
   end if

   !assign the shift to the atomic positions
   shift(1)=cxmin
   shift(2)=cymin
   shift(3)=czmin

   !here we can put a modulo operation for periodic directions
   do iat=1,atoms%astruct%nat
      rxyz(1,iat)=rxyz(1,iat)-shift(1)
      rxyz(2,iat)=rxyz(2,iat)-shift(2)
      rxyz(3,iat)=rxyz(3,iat)-shift(3)
   enddo

   ! fine grid size (needed for creation of input wavefunction, preconditioning)
   if (atoms%astruct%nat == 0) then
      !For homogeneous gaz, we fill the box with the fine grid
      nfl1=0 
      nfl2=0 
      nfl3=0

      nfu1=n1
      nfu2=n2
      nfu3=n3
   else
      !we start with nfl max to find th emin and nfu min to find the max
      nfl1=n1 
      nfl2=n2 
      nfl3=n3

      nfu1=0 
      nfu2=0 
      nfu3=0
   end if

   do iat=1,atoms%astruct%nat
      rad=radii_cf(atoms%astruct%iatype(iat),2)*frmult
      if (rad > 0.0_gp) then
         nfl1=min(nfl1,ceiling((rxyz(1,iat)-rad)/hx - eps_mach))
         nfu1=max(nfu1,floor((rxyz(1,iat)+rad)/hx + eps_mach))

         nfl2=min(nfl2,ceiling((rxyz(2,iat)-rad)/hy - eps_mach))
         nfu2=max(nfu2,floor((rxyz(2,iat)+rad)/hy + eps_mach))

         nfl3=min(nfl3,ceiling((rxyz(3,iat)-rad)/hz - eps_mach)) 
         nfu3=max(nfu3,floor((rxyz(3,iat)+rad)/hz + eps_mach))
      end if
   enddo

   !correct the values of the delimiter if they go outside the box
   if (nfl1 < 0 .or. nfu1 > n1) then
      nfl1=0
      nfu1=n1
   end if
   if (nfl2 < 0 .or. nfu2 > n2) then
      nfl2=0
      nfu2=n2
   end if
   if (nfl3 < 0 .or. nfu3 > n3) then
      nfl3=0
      nfu3=n3
   end if

   !correct the values of the delimiter if there are no wavelets
   if (nfl1 == n1 .and. nfu1 == 0) then
      nfl1=n1/2
      nfu1=n1/2
   end if
   if (nfl2 == n2 .and. nfu2 == 0) then
      nfl2=n2/2
      nfu2=n2/2
   end if
   if (nfl3 == n3 .and. nfu3 == 0) then
      nfl3=n3/2
      nfu3=n3/2
   end if

   if (iproc == 0) then
      if (atoms%astruct%ntypes > 0) then
         call yaml_comment('Atom Positions',hfill='-')
         call yaml_open_sequence('Atomic positions within the cell (Atomic and Grid Units)')
         do iat=1,atoms%astruct%nat
            call yaml_sequence(advance='no')
            call yaml_open_map(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),flow=.true.)
              call yaml_map('AU',rxyz(1:3,iat),fmt='(1pg12.5)')
              call yaml_map('GU',(/rxyz(1,iat)/hx,rxyz(2,iat)/hy,rxyz(3,iat)/hz/),fmt='(1pg12.5)')
            call yaml_close_map(advance='no')
            call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
         enddo
         call yaml_close_sequence()
         call yaml_map('Rigid Shift Applied (AU)',(/-cxmin,-cymin,-czmin/),fmt='(1pg12.5)')
      end if
      call yaml_comment('Grid properties',hfill='-')
      call yaml_map('Box Grid spacings',(/hx,hy,hz/),fmt='(f7.4)')
      call yaml_open_map('Sizes of the simulation domain')
        call yaml_map('AU',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),atoms%astruct%cell_dim(3)/),fmt='(1pg12.5)')
        call yaml_map('Angstroem',(/atoms%astruct%cell_dim(1)*Bohr_Ang,&
             atoms%astruct%cell_dim(2)*Bohr_Ang,atoms%astruct%cell_dim(3)*Bohr_Ang/),fmt='(1pg12.5)')
        call yaml_map('Grid Spacing Units',(/n1,n2,n3/),fmt='(i4)')
        call yaml_open_map('High resolution region boundaries (GU)',flow=.false.)
          call yaml_map('From',(/nfl1,nfl2,nfl3/),fmt='(i4)')
          call yaml_map('To',(/nfu1,nfu2,nfu3/),fmt='(i4)')
        call yaml_close_map()
      call yaml_close_map()
!!$      write(*,'(1x,a,19x,a)') 'Shifted atomic positions, Atomic Units:','grid spacing units:'
!!$      do iat=1,atoms%astruct%nat
!!$         write(*,'(1x,i5,1x,a6,3(1x,1pe12.5),3x,3(1x,0pf9.3))') &
!!$            &   iat,trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),&
!!$            &   (rxyz(j,iat),j=1,3),rxyz(1,iat)/hx,rxyz(2,iat)/hy,rxyz(3,iat)/hz
!!$      enddo
!!$      write(*,'(1x,a,3(1x,1pe12.5),a,3(1x,0pf7.4))') &
!!$         &   '   Shift of=',-cxmin,-cymin,-czmin,' H grids=',hx,hy,hz
!!$      write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
!!$         &   '  Box Sizes=',atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),atoms%astruct%cell_dim(3),n1,n2,n3
!!$      write(*,'(1x,a,3x,3(3x,i4,a1,i0))')&
!!$         &   '      Extremes for the high resolution grid points:',&
!!$         &   nfl1,'<',nfu1,nfl2,'<',nfu2,nfl3,'<',nfu3
   endif

   !assign the values
   Glr%geocode=atoms%astruct%geocode
   Glr%d%n1  =n1  
   Glr%d%n2  =n2  
   Glr%d%n3  =n3  
   Glr%d%n1i =n1i 
   Glr%d%n2i =n2i 
   Glr%d%n3i =n3i 
   Glr%d%nfl1=nfl1
   Glr%d%nfl2=nfl2
   Glr%d%nfl3=nfl3
   Glr%d%nfu1=nfu1
   Glr%d%nfu2=nfu2
   Glr%d%nfu3=nfu3

   Glr%ns1=0
   Glr%ns2=0
   Glr%ns3=0
   Glr%nsi1=0
   Glr%nsi2=0
   Glr%nsi3=0

   !while using k-points this condition should be disabled
   !evaluate if the condition for the hybrid evaluation if periodic BC hold
   Glr%hybrid_on=                   (nfu1-nfl1+lupfil < n1+1)
   Glr%hybrid_on=(Glr%hybrid_on.and.(nfu2-nfl2+lupfil < n2+1))
   Glr%hybrid_on=(Glr%hybrid_on.and.(nfu3-nfl3+lupfil < n3+1))

  !allocate projflg
!   allocate(Glr%projflg(atoms%astruct%nat),stat=i_stat)
!   call memocc(i_stat,Glr%projflg,'Glr%projflg',subname)
!   Glr%projflg = 1 
   
   !OCL convolutions not compatible with hybrid boundary conditions
   if (OCLConv) Glr%hybrid_on = .false.

!!$   if (Glr%hybrid_on) then
!!$      if (iproc == 0) write(*,*)'wavelet localization is ON'
!!$   else
!!$      if (iproc == 0) write(*,*)'wavelet localization is OFF'
!!$   endif
   if (iproc==0) call yaml_map('High Res. box is treated separately',Glr%hybrid_on)
END SUBROUTINE system_size


!>   Here the dimensions should be corrected in order to 
!!   allow the fft for the preconditioner and for Poisson Solver
subroutine correct_grid(a,h,n)
   use module_base
   use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
   implicit none
   real(gp), intent(in) :: a
   integer, intent(inout) :: n
   real(gp), intent(inout) :: h
   !local variables
   integer :: m,m2,nt

   n=ceiling(a/h)-1
   nt=n+1
   do
      !correct the direct dimension
      call fourier_dim(nt,m)

      !control if the double of this dimension is compatible with the FFT
      call fourier_dim(2*m,m2)
      !if this check is passed both the preconditioner and the PSolver works
      if (m2==2*m .and. mod(m,2) ==0) exit !only even dimensions are considered so far

      nt=m+1
   end do
   n=m-1

   !!!  !here the dimensions should be corrected in order to 
   !!!  !allow the fft for the preconditioner
   !!!  m=2*n+2
   !!!  do 
   !!!     call fourier_dim(m,m)
   !!!     if ((m/2)*2==m) then
   !!!        n=(m-2)/2
   !!!        exit
   !!!     else
   !!!        m=m+1
   !!!     end if
   !!!  end do

   h=a/real(n+1,gp)

END SUBROUTINE correct_grid


!> Calculates the length of the keys describing a wavefunction data structure
subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
   implicit none
   integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3
   logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid 
   integer, intent(out) :: mseg,mvctr
   !local variables
   logical :: plogrid
   integer :: i1,i2,i3,nsrt,nend,nsrti,nendi,mvctri
   mvctr=0
   nsrt=0
   nend=0
   !$omp parallel default(private) shared(nl3,nu3,nl2,nu2,nl1,nu1,logrid,mvctr,nsrt,nend)
   mvctri=0
   nsrti=0
   nendi=0
   !$omp do  
   do i3=nl3,nu3 
      do i2=nl2,nu2
         plogrid=.false.
         do i1=nl1,nu1
            if (logrid(i1,i2,i3)) then
               mvctri=mvctri+1
               if (.not. plogrid) then
                  nsrti=nsrti+1
               endif
            else
               if (plogrid) then
                  nendi=nendi+1
               endif
            endif
            plogrid=logrid(i1,i2,i3)
         enddo
         if (plogrid) then
            nendi=nendi+1
         endif
      enddo
   enddo
   !$omp enddo
   !$omp critical
   mvctr=mvctr+mvctri
   nsrt=nsrt+nsrti
   nend=nend+nendi
   !$omp end critical
   !$omp end parallel
   if (nend /= nsrt) then 
      write(*,*)' ERROR: nend <> nsrt',nend,nsrt
      stop 
   endif
   mseg=nend

END SUBROUTINE num_segkeys


!>   Calculates the keys describing a wavefunction data structure
subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
   implicit none
   integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,mseg
   logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid  
   integer, dimension(mseg), intent(out) :: keyv
   integer, dimension(2,mseg), intent(out) :: keyg
   !local variables
   logical :: plogrid
   integer :: mvctr,nsrt,nend,i1,i2,i3,ngridp

   mvctr=0
   nsrt=0
   nend=0
   do i3=nl3,nu3 
      do i2=nl2,nu2
         plogrid=.false.
         do i1=nl1,nu1
            ngridp=i3*((n1+1)*(n2+1)) + i2*(n1+1) + i1+1
            if (logrid(i1,i2,i3)) then
               mvctr=mvctr+1
               if (.not. plogrid) then
                  nsrt=nsrt+1
                  keyg(1,nsrt)=ngridp
                  keyv(nsrt)=mvctr
               endif
            else
               if (plogrid) then
                  nend=nend+1
                  keyg(2,nend)=ngridp-1
               endif
            endif
            plogrid=logrid(i1,i2,i3)
         enddo
         if (plogrid) then
            nend=nend+1
            keyg(2,nend)=ngridp
         endif
      enddo
   enddo
   if (nend /= nsrt) then 
      write(*,*) 'nend , nsrt',nend,nsrt
      stop 'nend <> nsrt'
   endif
   !mseg=nend
END SUBROUTINE segkeys


!> Set up an array logrid(i1,i2,i3) that specifies whether the grid point
!! i1,i2,i3 is the center of a scaling function/wavelet
subroutine fill_logrid(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
      &   ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
   use module_base
   implicit none
   !Arguments
   character(len=*), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
   integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
   real(gp), intent(in) :: rmult,hx,hy,hz
   integer, dimension(nat), intent(in) :: iatype
   real(gp), dimension(ntypes), intent(in) :: radii
   real(gp), dimension(3,nat), intent(in) :: rxyz
   logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid
   !local variables
   real(kind=8), parameter :: eps_mach=1.d-12
   integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3
   real(gp) :: dx,dy2,dz2,rad

   !some checks
   if (geocode(1:1) /= 'F') then
      !the nbuf value makes sense only in the case of free BC
      if (nbuf /=0) then
         write(*,'(1x,a)')'ERROR: a nonzero value of nbuf is allowed only for Free BC (tails)'
         stop
      end if
   else
      !The grid spacings must be the same
      if (hx /= hy .or. hy /= hz .or. hx /= hz) then
         !        write(*,'(1x,a)')'ERROR: For Free BC the grid spacings must be the same'
      end if
   end if

   if (geocode(1:1) == 'F') then
      do i3=nl3,nu3 
         do i2=nl2,nu2 
            do i1=nl1,nu1
               logrid(i1,i2,i3)=.false.
            enddo
         enddo
      enddo
   else !
      !Special case if no atoms (homogeneous electron gas): all points are used (TD)
      if (nat == 0) then
         do i3=0,n3 
            do i2=0,n2 
               do i1=0,n1
                  logrid(i1,i2,i3)=.true.
               enddo
            enddo
         enddo
      else
         do i3=0,n3 
            do i2=0,n2 
               do i1=0,n1
                  logrid(i1,i2,i3)=.false.
               enddo
            enddo
         enddo
      end if
   end if

   do iat=1,nat
      rad=radii(iatype(iat))*rmult+real(nbuf,gp)*hx
      if (rad /= 0.0_gp) then
         ml1=ceiling((rxyz(1,iat)-rad)/hx - eps_mach)  
         ml2=ceiling((rxyz(2,iat)-rad)/hy - eps_mach)   
         ml3=ceiling((rxyz(3,iat)-rad)/hz - eps_mach)   
         mu1=floor((rxyz(1,iat)+rad)/hx + eps_mach)
         mu2=floor((rxyz(2,iat)+rad)/hy + eps_mach)
         mu3=floor((rxyz(3,iat)+rad)/hz + eps_mach)

         !for Free BC, there must be no incoherences with the previously calculated delimiters
         if (geocode(1:1) == 'F') then
           if (ml1 < nl1) then
               write(*,'(a,i0,3x,i0)')  'ERROR: ml1 < nl1  ', ml1, nl1
               stop
           end if
           if (ml2 < nl2) then
               write(*,'(a,i0,3x,i0)')  'ERROR: ml2 < nl2  ', ml2, nl2
               stop
           end if
           if (ml3 < nl3) then
               write(*,'(a,i0,3x,i0)')  'ERROR: ml3 < nl3  ', ml3, nl3
               stop
           end if

           if (mu1 > nu1) then
               write(*,'(a,i0,3x,i0)')  'ERROR: mu1 > nu1  ', mu1, nu1
               stop
           end if
           if (mu2 > nu2) then
               write(*,'(a,i0,3x,i0)')  'ERROR: mu2 > nu2  ', mu2, nu2
               stop
           end if
           if (mu3 > nu3) then
               write(*,'(a,i0,3x,i0)')  'ERROR: mu3 > nu3  ', mu3, nu3
               stop
           end if
         end if
         !what follows works always provided the check before
         !$omp parallel default(shared) private(i3,dz2,j3,i2,dy2,j2,i1,j1,dx)
         !$omp do
         do i3=max(ml3,-n3/2-1),min(mu3,n3+n3/2+1)
            dz2=(real(i3,gp)*hz-rxyz(3,iat))**2
            j3=modulo(i3,n3+1)
            do i2=max(ml2,-n2/2-1),min(mu2,n2+n2/2+1)
               dy2=(real(i2,gp)*hy-rxyz(2,iat))**2
               j2=modulo(i2,n2+1)
               do i1=max(ml1,-n1/2-1),min(mu1,n1+n1/2+1)
                  j1=modulo(i1,n1+1)
                  dx=real(i1,gp)*hx-rxyz(1,iat)
                  if (dx**2+(dy2+dz2)-eps_mach <= rad**2) then 
                     logrid(j1,j2,j3)=.true.
                  endif
               enddo
            enddo
         enddo
         !$omp enddo
         !$omp end parallel
      end if
   enddo

END SUBROUTINE fill_logrid


subroutine make_bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
   implicit none
   integer, intent(in) :: n1,n2,n3
   logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid
   integer, dimension(2,0:n2,0:n3), intent(out) :: ibyz
   integer, dimension(2,0:n1,0:n3), intent(out) :: ibxz
   integer, dimension(2,0:n1,0:n2), intent(out) :: ibxy
   !local variables
   integer :: i1,i2,i3

   do i3=0,n3 
      do i2=0,n2 
         ibyz(1,i2,i3)= 1000
         ibyz(2,i2,i3)=-1000

         loop_i1s: do i1=0,n1
            if (logrid(i1,i2,i3)) then 
               ibyz(1,i2,i3)=i1
               exit loop_i1s
            endif
         enddo loop_i1s

         loop_i1e: do i1=n1,0,-1
            if (logrid(i1,i2,i3)) then 
               ibyz(2,i2,i3)=i1
               exit loop_i1e
            endif
         enddo loop_i1e
      end do
   end do


   do i3=0,n3 
      do i1=0,n1
         ibxz(1,i1,i3)= 1000
         ibxz(2,i1,i3)=-1000

         loop_i2s: do i2=0,n2 
            if (logrid(i1,i2,i3)) then 
               ibxz(1,i1,i3)=i2
               exit loop_i2s
            endif
         enddo loop_i2s

         loop_i2e: do i2=n2,0,-1
            if (logrid(i1,i2,i3)) then 
               ibxz(2,i1,i3)=i2
               exit loop_i2e
            endif
         enddo loop_i2e

      end do
   end do

   do i2=0,n2 
      do i1=0,n1 
         ibxy(1,i1,i2)= 1000
         ibxy(2,i1,i2)=-1000

         loop_i3s: do i3=0,n3
            if (logrid(i1,i2,i3)) then 
               ibxy(1,i1,i2)=i3
               exit loop_i3s
            endif
         enddo loop_i3s

         loop_i3e: do i3=n3,0,-1
            if (logrid(i1,i2,i3)) then 
               ibxy(2,i1,i2)=i3
               exit loop_i3e
            endif
         enddo loop_i3e
      end do
   end do

END SUBROUTINE make_bounds
