!
! Code are based on ABINIT psp5lo.F90 and mklocl_recipspace.F90 
! pspcore of the total system is 
!   N_{tot} * sum_i(epsatm_i) / (cell volume), where i loops over all atoms
!              

subroutine make_subloc(nsys, psps, sub_vloc, gmet, nfft, cell_nfft, cell_acell, xatom, epsatm, ucvol)

 use mpi
 use bspline

 implicit none 

 character(len=100) :: psps(nsys)
 integer :: nsys,cell_nfft(3)
 real(8) :: sub_vloc(nfft,nsys), & 
            gmet(3,3), epsatm(nsys),  & 
            ucvol, xatom(3,nsys), cell_acell(3)

 ! internal vars 

 integer,parameter  :: mqgrid = 5000
 character(len=100) :: scmd 
 integer :: pspcode,j,jj,ii,ir,mmax,mmax2,lloc,nstart,nend,dd,iq, & 
            nfft, n1,n2,n3,i1,i2,i3,ia,ia2,ia1,id1,id2,id3,ig1,ig2,ig3,im=2,re=1

 integer :: ider,irmu,irn,mesh_mult,mmax_new
 integer,parameter :: NPT_IN_2PI=200
 real(8) :: xx,amesh,amesh_new,r0,fp1,fpn,qmesh

 real(8),allocatable :: vloc(:),rad(:),rad_new(:),work(:), & 
                        rvlpz(:),rvlpz_new(:),sprvlpz(:,:)

 real(8) :: rad_expo, result, arg, pi=3.1415926d0, test, & 
            rmtoin, zero=0.0d0, zion, ztor1, ratio, tolfix=1.0000001d0, & 
            qgrid(mqgrid), vlspl(mqgrid,2), q2vq(mqgrid), & 
            work_space(mqgrid), work_spl(mqgrid), work1(2,nfft),& 
            gmag, gsqcut, vion1, qmax, & 
            yp1, ypn, diff, aa, bb, cc, cutoff,  & 
            dq, dq2div6, dqdiv6, dqm1, gsquar

 real(8) :: gsq,gsq_max,dtmp,sfi,sfr

 integer :: ierr, myrank
 complex(8) :: g_fftw(cell_nfft(1)/2+1,cell_nfft(2),cell_nfft(3))


 ! ========================================
 ! make gcut^2 and strucutre factor 
 ! ========================================


 ! based on file mklocl_recipspace.F90
 ! Define G^2 based on G space metric gmet.
 gsq(i1,i2,i3)=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
             & dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
             & dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)

 print *,'enter make_subvext()'

 n1 = cell_nfft(1)
 n2 = cell_nfft(2)
 n3 = cell_nfft(3)
 nfft = n1*n2*n3

 ! get the maximum gsq 
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 gsq_max = 0.0d0

 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     do i1=1,n1
       ig1=i1-(i1/id1)*n1-1
       gsquar=gsq(ig1,ig2,ig3)
       if ( gsq_max < gsquar ) then 
         gsq_max = gsquar
       endif
     end do
   end do
 end do
 ! our definition of gmet is (2pi)^2 larger than ABINIT's definition
 gsq_max = gsq_max / (2.d0*pi)**2
 !print *,'gsq_max: ',gsq_max


 !
 ! setup qmax, mqgrid, q2vq, vlspl arrays
 !
 !Set up uniform q grids, make qmax 20% larger than largest expected:
 qmax=1.2d0 * sqrt(gsq_max)
 dq=qmax/dble(mqgrid-1)
 do ii=1,mqgrid
   qgrid(ii)=(ii-1)*dq
 end do
 !print *,'qmax: ',qmax
 !print *,'dq: ',dq

 !
 ! load the psps and make the vlocal for this atom 
 ! in the cell (G=0 is set to zero)
 !
 do j=1,nsys

   ! ==================================
   ! load the local part of *.fhi 
   ! ==================================

   write(6,'(a,a)')' psps file:  ',trim(psps(j))

   ! pspfile code 
   if (myrank==0) then 
     write(scmd,'(a,a,a)')"head -n 3 ",trim(psps(j))," | tail -n 1 | awk '{print $1}' > .tmp_oo"
     call system( scmd )
   endif
   call mpi_barrier(mpi_comm_world, ierr)
   open(file='.tmp_oo',action='read',unit=111,form='formatted')
   read(111,*) pspcode
   close(111)
   write(6,'(a,i4)')' pspscode: ',pspcode

   ! ------------------------------------
   !  pspcode = 6, *.cpi file of ABIINT 
   ! ------------------------------------

   if (pspcode==6) then 
     ! array length of *.fhi file
     if (myrank==0) then
       write(scmd,'(a,a,a)')"head -n 3 ",trim(psps(j))," | tail -n 1 | awk '{print $5}' > .tmp_oo"
       call system( scmd )
     endif 
     call mpi_barrier(mpi_comm_world, ierr)
     open(file='.tmp_oo',action='read',unit=111,form='formatted')
     read(111,*)mmax
     close(111)
     write(6,'(a,i4)')' mmax: ',mmax

     ! zion from fhi file 
     if (myrank==0) then 
       write(scmd,'(a,a,a)')"head -n 2 ",trim(psps(j))," | tail -n 1 | awk '{print $2}' > .tmp_oo"
       call system( scmd )
     endif
     call mpi_barrier(mpi_comm_world, ierr)
     open(file='.tmp_oo',action='read',unit=111,form='formatted')
     read(111,*)zion
     close(111)
     write(6,'(a,f8.4)')' zion: ',zion

     ! which l is local 
     if (myrank==0) then 
       write(scmd,'(a,a,a)')"head -n 3 ",trim(psps(j))," | tail -n 1 | awk '{print $4}' > .tmp_oo"
       call system( scmd )
     endif
     call mpi_barrier(mpi_comm_world, ierr)
     open(file='.tmp_oo',action='read',unit=111,form='formatted')
     read(111,*)lloc
     close(111)
     write(6,'(a,i4)')' local ang. mom.: ',lloc

     ! load local potential
     allocate(vloc(mmax),rad(mmax),work(mmax))
     nstart = 18 + (lloc+1)*(mmax+1)
     nend = mmax
     if (myrank==0) then 
       write(scmd,'(a,i4,a,a,a,i4,a)')"head -n ",nstart, & 
         " ",trim(psps(j))," | tail -n ",nend," > .tmp_oo"
       call system( scmd )
     endif
     call mpi_barrier(mpi_comm_world, ierr)
     open(file='.tmp_oo',action='read',unit=111,form='formatted')
     do jj=1,mmax
       read(111,*)dd,rad(jj),dtmp,vloc(jj)
     enddo 
     close(111)
     ratio = rad(mmax)/rad(1)
     rad_expo = log(ratio)/dble(mmax-1)
     print *,'rad_expo: ',rad_expo

     ! ==================================================
     !   compute q^2 v(q) for each type of atom. 
     !==================================================
 
     ! code below is based on psp5lo.F90 in ABINIT
 
     ! Do q=0 separately (compute epsatm)
     ! Do integral from 0 to r1
     ztor1=(zion/2.0d0+rad(1)*vloc(1)/3.d0)*rad(1)**2
 
     !Set up integrand for q=0: $ \int[r^2 (V(r)+\frac{Zv}{r}) dr]$
     !with extra factor of r to convert to uniform grid in exponent
     do ir=1,mmax
       !  First handle tail region
       test=vloc(ir)+zion/rad(ir)
       !  DEBUG
       !  write(std_out,*)ir,rad(ir),test
       !  ENDDEBUG
       !  Ignore small contributions, or impose a cut-off in the case
       !  the pseudopotential data are in single precision.
       !  (it is indeed expected that vloc is very close to zero beyond 20,
       !  so a value larger than 2.0d-8 is considered anomalous)
       if (abs(test)<1.0d-20 .or. (rad(ir)>20.0d0 .and. abs(test)>2.0d-8) ) then
         work(ir)=zero
       else
         work(ir)=(rad(ir)*rad(ir))*(rad(ir)*vloc(ir)+zion)
       end if
     end do
  
     !Do integral from r(1) to r(max)
     call ctrap(mmax,work,rad_expo,result)
     !Do integral from r(mmax) to infinity
     !compute decay length lambda at r(mmax)
     !$\lambda=-\log((rad(im1)*vloc(im1)+zion)$/ &
     !$(rad(imat)*vloc(imat)+zion))/(rad(im1)-rad(imat))$
     !rmtoin=$(rad(mmax)*vloc(mmax)+zion)*(rad(mmax)+1.d0/\lambda)/\lambda$
     !Due to inability to fit exponential decay to r*V(r)+Zv
     !in tail, NO TAIL CORRECTION IS APPLIED
     !(numerical trouble might be removed if atomic code is
     !cleaned up in tail region)
     rmtoin=0.0d0
     epsatm(j)=4.d0*pi*(result+ztor1+rmtoin)
     print *,'epsatm: ',epsatm(j)
     q2vq(1)=-zion/pi
  
     !Loop over the rest of q values
     do iq=2,mqgrid
       arg=2.d0*pi*qgrid(iq)
  !  ztor1=$ -Zv/\pi+2q \int_0^{r1}[\sin(2\pi q r)(rV(r)+Zv) dr]$
       ztor1=(vloc(1)*sin(arg*rad(1))/arg-(rad(1)*vloc(1)+zion)* &
       &   cos(arg*rad(1)) )/pi
  
       !  set up integrand
       do  ir=1,mmax
         test=vloc(ir)+zion/rad(ir)
  !    Ignore contributions within decade of machine precision
  !       if ((scale+abs(test)).eq.scale) then
         if ( abs(test) < 1e-8 ) then
           work(ir)=zero
         else
           work(ir)=rad(ir)*sin(arg*rad(ir))*(rad(ir)*vloc(ir)+zion)
         end if
       end do
  !  do integral from r(1) to r(mmax)
       call ctrap(mmax,work,rad_expo,result)
  
  !  do integral from r(mmax) to infinity
  !  rmtoin=(r(mmax)*vr(mmax)+zion)*(lambda*sin(arg*r(mmax))+
  !  arg*cos(arg*r(mmax)))/(arg**2+lambda**2)
  !  See comment above; no tail correction
       rmtoin=0.0d0
  
  !  store q^2 v(q)
       q2vq(iq)=ztor1+2.d0*qgrid(iq)*(result+rmtoin)
  
     end do

     ! ==============================================================
     !   compute derivative of q^2 v(q) at two ends for later spline
     ! ==============================================================
  
     !Compute derivatives of q^2 v(q) at ends of interval
     yp1=0.0d0
     !ypn=$ 2\int_0^\infty[(\sin(2\pi qmax r)+(2\pi qmax r)*\cos(2\pi qmax r)(r V(r)+Z) dr]$
     !integral from 0 to r1
     arg=2.0d0*pi*qgrid(mqgrid)
     ztor1=zion*rad(1)*sin(arg*rad(1))
     ztor1=ztor1+ 3.d0*rad(1)*vloc(1)*cos(arg*rad(1))/arg + &
     & (rad(1)**2-1.0d0/arg**2)*vloc(1)*sin(arg*rad(1))
     !integral from r(mmax) to infinity is overkill; ignore
     !set up integrand
     do ir=1,mmax
       test=vloc(ir)+zion/rad(ir)    
     !  Ignore contributions within decade of machine precision
     !  if ((scale+abs(test)).eq.scale) then
       if ( abs(test) < 1e-8 ) then
         work(ir)=0.0d0
       else
         work(ir)=rad(ir)*(sin(arg*rad(ir))+arg*rad(ir)*cos(arg*rad(ir))) * &
  &      (rad(ir)*vloc(ir)+zion)
       end if
     end do
     call ctrap(mmax,work,rad_expo,result)
     ypn=2.0d0 * (ztor1 + result)

   endif 

   ! ------------------------------
   !  pspcode = 8 local psp 
   ! ------------------------------

   if (pspcode == 8) then 

     print *,''
     print *,' --------- local psp (pspcode=8) -----------'
     print *,''
     ! array length of *.fhi file
     if (myrank==0) then 
       write(scmd,'(a,a,a)')"head -n 3 ",trim(psps(j))," | tail -n 1 | awk '{print $5}' > .tmp_oo"
       call system( scmd )
     endif
     call mpi_barrier(mpi_comm_world, ierr)
     open(file='.tmp_oo',action='read',unit=111,form='formatted')
     read(111,*)mmax
     close(111)
     write(6,'(a,i4)')' mmax: ',mmax

     ! zion from fhi file 
     if (myrank==0) then 
       write(scmd,'(a,a,a)')"head -n 2 ",trim(psps(j))," | tail -n 1 | awk '{print $2}' > .tmp_oo"
       call system( scmd )
     endif
     call mpi_barrier(mpi_comm_world, ierr)
     open(file='.tmp_oo',action='read',unit=111,form='formatted')
     read(111,*)zion
     close(111)
     write(6,'(a,f8.4)')' zion: ',zion

     ! load local potential
     allocate(vloc(mmax),rad(mmax),work(mmax),rvlpz(mmax))
     nstart = 8+mmax
     nend = mmax
     if (myrank==0) then 
       write(scmd,'(a,i4,a,a,a,i4,a)')"head -n ",nstart, & 
         " ",trim(psps(j))," | tail -n ",nend," > .tmp_oo"
       call system( scmd )
     endif
     call mpi_barrier(mpi_comm_world, ierr)
     open(file='.tmp_oo',action='read',unit=111,form='formatted')
     do jj=1,mmax
       read(111,*)dd,rad(jj),vloc(jj)
     enddo 
     close(111)
     
     amesh = rad(2)-rad(1)
     !Do q=0 separately (compute epsatm)
     ztor1=(zion/2.0d0+rad(1)*vloc(1)/3.d0)*rad(1)**2
     !Set up integrand for q=0: $ \int[r^2 (V(r)+\frac{Zv}{r}) dr]$
     do ir=1,mmax
       rvlpz(ir)=rad(ir)*vloc(ir)+zion
       work(ir)=rad(ir)*rvlpz(ir)
     end do

     !Do integral from zero to r(max)
     call ctrap(mmax,work,amesh,result)

     epsatm(j)=4.d0*pi*result
     q2vq(1)=-zion/pi

     print *,'make_subloc.f90: epsatm-> ',epsatm(j)

     !Find r mesh spacing necessary for accurate integration at qmax
     amesh_new=2.d0*pi/(NPT_IN_2PI*qgrid(mqgrid))

     !Choose submultiple of input mesh
     mesh_mult=int(amesh/amesh_new) + 1
     mmax_new=mesh_mult*(mmax-1)+1
     amesh_new=amesh/dble(mesh_mult)

     allocate(rad_new(mmax_new))
     allocate(rvlpz_new(mmax_new))

     if(mesh_mult==1) then
       rad_new(:)=rad(:)
       rvlpz_new(:)=rvlpz(:)
     else
       !  Set up spline and interpolate to finer mesh.
       !  First, compute derivatives at end points
       fp1=(-50.d0*rvlpz(1)+96.d0*rvlpz(2)-72.d0*rvlpz(3)+32.d0*rvlpz(4)&
  &      -6.d0*rvlpz(5))/(24.d0*amesh)
       fpn=(6.d0*rvlpz(mmax-4)-32.d0*rvlpz(mmax-3)+72.d0*rvlpz(mmax-2)&
  &      -96.d0*rvlpz(mmax-1)+50.d0*rvlpz(mmax))/(24.d0*amesh)
       allocate(sprvlpz(mmax,2))
       work(:)=zero
  
       !  Spline fit
       call abinit_spline(rad, rvlpz,mmax,fp1,fpn,sprvlpz(:,2))
       sprvlpz(:,1)=rvlpz(:)

       !  Set up new radial mesh
       irn=1
       do ir=1,mmax-1
         do irmu=0,mesh_mult-1
           rad_new(irn)=rad(ir)+dble(irmu)*amesh_new
           irn=irn+1
         end do
       end do
       rad_new(mmax_new)=rad(mmax)

       ider=0
       call abinit_splfit(rad,work,sprvlpz,ider,rad_new,rvlpz_new,mmax,mmax_new)

       deallocate(sprvlpz,work)
       allocate(work(mmax_new))
     end if

     !Loop over q values
     do iq=2,mqgrid
       arg=2.d0*pi*qgrid(iq)

       !  Set up integrand
       do  ir=1,mmax_new
         work(ir)=sin(arg*rad_new(ir))*rvlpz_new(ir)
       end do

       !  Do integral from zero to rad(mmax)
       call ctrap(mmax_new,work,amesh_new,result)

       !  Store q^2 v(q)
       q2vq(iq)=q2vq(1)+2.d0*qgrid(iq)*result
     end do

     !Compute derivatives of q^2 v(q) at ends of interval
     qmesh=qgrid(2)-qgrid(1)
     yp1=(-50.d0*q2vq(1)+96.d0*q2vq(2)-72.d0*q2vq(3)+32.d0*q2vq(4)&
  &       -6.d0*q2vq(5))/(24.d0*qmesh)
     ypn=(6.d0*q2vq(mqgrid-4)-32.d0*q2vq(mqgrid-3)+72.d0*q2vq(mqgrid-2)&
  &       -96.d0*q2vq(mqgrid-1)+50.d0*q2vq(mqgrid))/(24.d0*qmesh)
  
      deallocate(rad_new,rvlpz_new)
      deallocate(rvlpz) !! by MM
   endif

  
   ! =====================
   !  spline q2vq 
   ! =====================

   !Fit spline to q^2 V(q) (Numerical Recipes subroutine)
   vlspl(:,1) = q2vq(:)
   call abinit_spline (qgrid,vlspl(:,1),mqgrid,yp1,ypn,work_spl)
   vlspl(:,2) = work_spl(:)

   ! ==================================
   !  construct vlocal for each atom 
   ! ==================================

   ! based on file mklocl_recipspace.F90

   !Zero out array to permit accumulation over atom types below:
   work1(:,:)=zero

   dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
   dqm1=1.0d0/dq
   dqdiv6=dq/6.0d0
   dq2div6=dq**2/6.0d0
   cutoff=gsqcut*tolfix
   id1=n1/2+2
   id2=n2/2+2
   id3=n3/2+2
   ia1=1

   ii=0
   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     do i2=1,n2
        ig2=i2-(i2/id2)*n2-1
        do i1=1,n1
           ig1=i1-(i1/id1)*n1-1

           ii=ii+1
!          ***     GET RID OF THIS THESE IF STATEMENTS (if they slow code)
!          Skip G=0:
           if (ii==1) cycle
            
           gsquar=gsq(ig1,ig2,ig3)
!
! chen: our gmet = gmet_ABIINT * (2*pi)**2
! in fact, our gmet is defined based on the k in k*r, not the q in 2 pi q*r
! FFTW performs (2 pi q r), that is why we should divide gmet by 4pi^2
!
           gsquar=gsquar / (2.d0*pi)**2
!          Skip G**2 outside cutoff:
!           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)

!            Compute vion(G) for given type of atom
             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Evaluate spline fit from q^2 V(q) to get V(q):
!            (p. 86 Numerical Recipes, Press et al;
!            NOTE error in book for sign
!            of "aa" term in derivative; also see splfit routine).

             bb = diff*dqm1
             aa = 1.0d0-bb
             cc = aa*(aa**2-1.0d0)*dq2div6
             dd = bb*(bb**2-1.0d0)*dq2div6

             vion1 = (aa*vlspl(jj,1)+bb*vlspl(jj+1,1) +&
        &             cc*vlspl(jj,2)+dd*vlspl(jj+1,2) ) / gsquar

!
!  structure factor to 1.d0 for a orthogonal cube here
!
             dtmp = 2.d0*pi/cell_acell(1)*ig1*xatom(1,j) +  & 
                    2.d0*pi/cell_acell(2)*ig2*xatom(2,j) + & 
                    2.d0*pi/cell_acell(3)*ig3*xatom(3,j)

! the minus sign is due to the e(-i 2 pi q*r) definition in FFTW

             sfr = cos(-dtmp)
             sfi = sin(-dtmp)

!            Multiply structure factor times vion:
             work1(re,ii)=sfr*vion1
             work1(im,ii)=sfi*vion1
!            End skip G**2 outside cutoff:
!           end if

!          End loop on n1, n2, n3. There is a "cycle" inside the loop
        end do
     end do
   end do


!  convert work1 (abinit) array to fftw3 array
   call array_abinit2fftw3(work1,n1,n2,n3,g_fftw)

   call fft(n1,n2,n3,sub_vloc(:,j),g_fftw,-1)

   deallocate(vloc,rad,work)

!  convert sub_vloc to correct unit. This is tricky, generally speaking, in ABINIT 
!  we divide sub_vloc by the cell volume, however in our fft()
!  function, we also divided by the (n1*n2*n3), therefore, 
!  we first times (n1*n2*n3) back and then divide it by ucvol
!  now everything is consistent with ABINIT's mklocl_recipspace.F90 file 

   sub_vloc(:,j)=sub_vloc(:,j)*dble(n1*n2*n3)/ucvol

   print *,'sub_vloc: ',minval(sub_vloc),maxval(sub_vloc)
   print *,''

 enddo

 print *,'leave make_subloc()'

end subroutine make_subloc
