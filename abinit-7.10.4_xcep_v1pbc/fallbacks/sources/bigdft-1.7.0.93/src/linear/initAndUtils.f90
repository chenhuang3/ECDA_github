!> @file 
!!   Initializations
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 

subroutine allocateBasicArraysInputLin(lin, ntypes)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(linearInputParameters),intent(inout) :: lin
  integer, intent(in) :: ntypes
  
  ! Local variables
  integer :: istat
  character(len=*),parameter :: subname='allocateBasicArrays'
  
  allocate(lin%norbsPerType(ntypes), stat=istat)
  call memocc(istat, lin%norbsPerType, 'lin%norbsPerType', subname)

  allocate(lin%potentialPrefac_ao(ntypes), stat=istat)
  call memocc(istat, lin%potentialPrefac_ao, 'lin%potentialPrefac_ao', subname)

  allocate(lin%potentialPrefac_lowaccuracy(ntypes), stat=istat)
  call memocc(istat, lin%potentialPrefac_lowaccuracy, 'lin%potentialPrefac_lowaccuracy', subname)

  allocate(lin%potentialPrefac_highaccuracy(ntypes), stat=istat)
  call memocc(istat, lin%potentialPrefac_highaccuracy, 'lin%potentialPrefac_highaccuracy', subname)

  allocate(lin%locrad_type(ntypes),stat=istat)
  call memocc(istat,lin%locrad_type,'lin%locrad_type',subname)

  allocate(lin%kernel_cutoff(ntypes), stat=istat)
  call memocc(istat, lin%kernel_cutoff, 'lin%kernel_cutoff', subname)

end subroutine allocateBasicArraysInputLin

subroutine deallocateBasicArraysInput(lin)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(linearinputParameters),intent(inout) :: lin
  
  ! Local variables
  integer :: i_stat,i_all
  character(len=*),parameter :: subname='deallocateBasicArrays'
 
  if(associated(lin%potentialPrefac_ao)) then
    i_all = -product(shape(lin%potentialPrefac_ao))*kind(lin%potentialPrefac_ao)
    deallocate(lin%potentialPrefac_ao,stat=i_stat)
    call memocc(i_stat,i_all,'lin%potentialPrefac_ao',subname)
    nullify(lin%potentialPrefac_ao)
  end if 
  if(associated(lin%potentialPrefac_lowaccuracy)) then
    i_all = -product(shape(lin%potentialPrefac_lowaccuracy))*kind(lin%potentialPrefac_lowaccuracy)
    deallocate(lin%potentialPrefac_lowaccuracy,stat=i_stat)
    call memocc(i_stat,i_all,'lin%potentialPrefac_lowaccuracy',subname)
    nullify(lin%potentialPrefac_lowaccuracy)
  end if 
  if(associated(lin%potentialPrefac_highaccuracy)) then
    i_all = -product(shape(lin%potentialPrefac_highaccuracy))*kind(lin%potentialPrefac_highaccuracy)
    deallocate(lin%potentialPrefac_highaccuracy,stat=i_stat)
    call memocc(i_stat,i_all,'lin%potentialPrefac_highaccuracy',subname)
    nullify(lin%potentialPrefac_highaccuracy)
  end if 

  if(associated(lin%norbsPerType)) then
    i_all = -product(shape(lin%norbsPerType))*kind(lin%norbsPerType)
    deallocate(lin%norbsPerType,stat=i_stat)
    call memocc(i_stat,i_all,'lin%norbsPerType',subname)
    nullify(lin%norbsPerType)
  end if 

  if(associated(lin%locrad)) then
    i_all = -product(shape(lin%locrad))*kind(lin%locrad)
    deallocate(lin%locrad,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad',subname)
    nullify(lin%locrad)
  end if 

  if(associated(lin%locrad_lowaccuracy)) then
    i_all = -product(shape(lin%locrad_lowaccuracy))*kind(lin%locrad_lowaccuracy)
    deallocate(lin%locrad_lowaccuracy,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad_lowaccuracy',subname)
    nullify(lin%locrad_lowaccuracy)
  end if 

  if(associated(lin%locrad_highaccuracy)) then
    i_all = -product(shape(lin%locrad_highaccuracy))*kind(lin%locrad_highaccuracy)
    deallocate(lin%locrad_highaccuracy,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad_highaccuracy',subname)
    nullify(lin%locrad_highaccuracy)
  end if 

  if(associated(lin%locrad_type)) then
    i_all = -product(shape(lin%locrad_type))*kind(lin%locrad_type)
    deallocate(lin%locrad_type,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad_type',subname)
    nullify(lin%locrad_type)
  end if 

  if(associated(lin%kernel_cutoff)) then
    i_all = -product(shape(lin%kernel_cutoff))*kind(lin%kernel_cutoff)
    deallocate(lin%kernel_cutoff,stat=i_stat)
    call memocc(i_stat,i_all,'lin%kernel_cutoff',subname)
    nullify(lin%kernel_cutoff)
  end if 

end subroutine deallocateBasicArraysInput



! lzd%llr already allocated, locregcenter and locrad already filled - could tidy this!
subroutine initLocregs(iproc, nproc, lzd, hx, hy, hz, astruct, orbs, Glr, locregShape, lborbs)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initLocregs
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(inout) :: lzd
  real(kind=8),intent(in) :: hx, hy, hz
  type(atomic_structure),intent(in) :: astruct
  type(orbitals_data),intent(in) :: orbs
  type(locreg_descriptors),intent(in) :: Glr
  character(len=1),intent(in) :: locregShape
  type(orbitals_data),optional,intent(in) :: lborbs
  
  ! Local variables
  integer :: istat, jorb, jjorb, jlr, iall
  character(len=*),parameter :: subname='initLocregs'
  logical,dimension(:),allocatable :: calculateBounds

  
  allocate(calculateBounds(lzd%nlr), stat=istat)
  call memocc(istat, calculateBounds, 'calculateBounds', subname)
  calculateBounds=.false.
  
  do jorb=1,orbs%norbp
     jjorb=orbs%isorb+jorb
     jlr=orbs%inWhichLocreg(jjorb)
     calculateBounds(jlr)=.true.
  end do

  if(present(lborbs)) then
     do jorb=1,lborbs%norbp
        jjorb=lborbs%isorb+jorb
        jlr=lborbs%inWhichLocreg(jjorb)
        calculateBounds(jlr)=.true.
     end do
  end if
  
  if(locregShape=='c') then
      stop 'locregShape c is deprecated'
  else if(locregShape=='s') then
      call determine_locregSphere_parallel(iproc, nproc, lzd%nlr, hx, hy, hz, &
           astruct, orbs, Glr, lzd%Llr, calculateBounds)
  end if
  
  iall=-product(shape(calculateBounds))*kind(calculateBounds)
  deallocate(calculateBounds, stat=istat)
  call memocc(istat, iall, 'calculateBounds', subname)
  
  !DEBUG
  !do ilr=1,lin%nlr
  !    if(iproc==0) write(*,'(1x,a,i0)') '>>>>>>> zone ', ilr
  !    if(iproc==0) write(*,'(3x,a,4i10)') 'nseg_c, nseg_f, nvctr_c, nvctr_f', lin%Llr(ilr)%wfd%nseg_c, lin%Llr(ilr)%wfd%nseg_f, lin%Llr(ilr)%wfd%nvctr_c, lin%Llr(ilr)%wfd%nvctr_f
  !    if(iproc==0) write(*,'(3x,a,3i8)') 'lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i', lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i
  !    if(iproc==0) write(*,'(a,6i8)') 'lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3',&
  !    lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3
  !end do
  !END DEBUG
  
  lzd%linear=.true.

end subroutine initLocregs

function megabytes(bytes)
  implicit none
  
  integer,intent(in) :: bytes
  integer :: megabytes
  
  megabytes=nint(dble(bytes)/1048576.d0)
  
end function megabytes


subroutine init_foe(iproc, nproc, lzd, astruct, input, orbs_KS, orbs, foe_obj, reset)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(atomic_structure),intent(in) :: astruct
  type(input_variables),intent(in) :: input
  type(orbitals_data),intent(in) :: orbs_KS, orbs
  type(foe_data),intent(out) :: foe_obj
  logical, intent(in) :: reset
  
  ! Local variables
  integer :: iorb, iiorb, jjorb, istat, iseg, ilr, jlr
  integer :: iwa, jwa, itype, jtype, ierr, iall
  logical :: seg_started
  real(kind=8) :: tt, cut
  logical,dimension(:,:),allocatable :: kernel_locreg
  character(len=*),parameter :: subname='initMatrixCompression'
!  integer :: ii, iseg
  
  call timing(iproc,'init_matrCompr','ON')

  if (reset) then
     foe_obj%ef=0.d0
     foe_obj%evlow=input%lin%evlow
     foe_obj%evhigh=input%lin%evhigh
     foe_obj%bisection_shift=1.d-1
     foe_obj%fscale=input%lin%fscale
     foe_obj%ef_interpol_det=input%lin%ef_interpol_det
     foe_obj%ef_interpol_chargediff=input%lin%ef_interpol_chargediff
     foe_obj%charge=0.d0
     do iorb=1,orbs_KS%norb
          foe_obj%charge=foe_obj%charge+orbs_KS%occup(iorb)
     end do
  end if

  call nullify_foe(foe_obj)

  ! Initialize kernel_locreg
  if (input%lin%scf_mode==LINEAR_FOE) then ! otherwise don't need to allocate just nullify as above
     allocate(kernel_locreg(orbs%norbp,orbs%norb), stat=istat)
     call memocc(istat, kernel_locreg, 'kernel_locreg', subname)
     allocate(foe_obj%kernel_nseg(orbs%norb), stat=istat)
     call memocc(istat, foe_obj%kernel_nseg, 'foe_obj%kernel_nseg', subname)
     call to_zero(orbs%norb, foe_obj%kernel_nseg(1))
     do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        ilr=orbs%inwhichlocreg(iiorb)
        iwa=orbs%onwhichatom(iiorb)
        itype=astruct%iatype(iwa)
        foe_obj%kernel_nseg(iiorb)=0
        seg_started=.false.
        do jjorb=1,orbs%norb
           jlr=orbs%inwhichlocreg(jjorb)
           jwa=orbs%onwhichatom(jjorb)
           jtype=astruct%iatype(jwa)
           tt = (lzd%llr(ilr)%locregcenter(1)-lzd%llr(jlr)%locregcenter(1))**2 + &
                (lzd%llr(ilr)%locregcenter(2)-lzd%llr(jlr)%locregcenter(2))**2 + &
                (lzd%llr(ilr)%locregcenter(3)-lzd%llr(jlr)%locregcenter(3))**2
           cut = input%lin%kernel_cutoff(itype)+input%lin%kernel_cutoff(jtype)
           tt=sqrt(tt)
           if (tt<=cut) then
              kernel_locreg(iorb,jjorb)=.true.
              if (.not.seg_started) then
                 foe_obj%kernel_nseg(iiorb)=foe_obj%kernel_nseg(iiorb)+1
              end if
              seg_started=.true.
           else
              kernel_locreg(iorb,jjorb)=.false.
              seg_started=.false.
           end if
        end do
     end do
     call mpiallred(foe_obj%kernel_nseg(1), orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)

     allocate(foe_obj%kernel_segkeyg(2,maxval(foe_obj%kernel_nseg),orbs%norb), stat=istat)
     call memocc(istat, foe_obj%kernel_segkeyg, 'foe_obj%kernel_segkeyg', subname)
     call to_zero(2*maxval(foe_obj%kernel_nseg)*orbs%norb, foe_obj%kernel_segkeyg(1,1,1))
     do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        iseg=0
        seg_started=.false.
        do jjorb=1,orbs%norb
           if(kernel_locreg(iorb,jjorb)) then
              if (.not.seg_started) then
                 iseg=iseg+1
                 foe_obj%kernel_segkeyg(1,iseg,iiorb)=jjorb
              end if
              seg_started=.true.
           else
              if (seg_started) then
                 foe_obj%kernel_segkeyg(2,iseg,iiorb)=jjorb-1
              end if
              seg_started=.false.
           end if
        end do
        if (seg_started) then
           foe_obj%kernel_segkeyg(2,iseg,iiorb)=orbs%norb
        end if
     end do
     call mpiallred(foe_obj%kernel_segkeyg(1,1,1), 2*maxval(foe_obj%kernel_nseg)*orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)

     iall = -product(shape(kernel_locreg))*kind(kernel_locreg) 
     deallocate(kernel_locreg,stat=istat)
     call memocc(istat,iall,'kernel_locreg',subname)
  end if

  call timing(iproc,'init_matrCompr','OF')


end subroutine init_foe


subroutine check_linear_and_create_Lzd(iproc,nproc,linType,Lzd,atoms,orbs,nspin,rxyz)
  use module_base
  use module_types
  use module_xc
  implicit none

  integer, intent(in) :: iproc,nproc,nspin
  type(local_zone_descriptors), intent(inout) :: Lzd
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data),intent(inout) :: orbs
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  integer, intent(in) :: linType
!  real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
  !Local variables
  character(len=*), parameter :: subname='check_linear_and_create_Lzd'
  logical :: linear
  integer :: iat,ityp,nspin_ig,i_all,i_stat,ilr
  real(gp), dimension(:), allocatable :: locrad
  logical,dimension(:),allocatable :: calculateBounds

  !default variables
  Lzd%nlr = 1

  if (nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=nspin
  end if

  linear  = .true.
  if (linType == INPUT_IG_FULL) then
     Lzd%nlr=atoms%astruct%nat
     allocate(locrad(Lzd%nlr+ndebug),stat=i_stat)
     call memocc(i_stat,locrad,'locrad',subname)
     ! locrad read from last line of  psppar
     do iat=1,atoms%astruct%nat
        ityp = atoms%astruct%iatype(iat)
        locrad(iat) = atoms%rloc(ityp,1)
     end do  
     call timing(iproc,'check_IG      ','ON')
     call check_linear_inputguess(iproc,Lzd%nlr,rxyz,locrad,&
          Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
          Lzd%Glr,linear) 
     call timing(iproc,'check_IG      ','OF')
     if(nspin >= 4) linear = .false. 
  end if

  ! If we are using cubic code : by choice or because locregs are too big
  Lzd%linear = .true.
  if (linType == INPUT_IG_LIG .or. linType == INPUT_IG_OFF .or. .not. linear) then
     Lzd%linear = .false.
     Lzd%nlr = 1
  end if


  if(linType /= INPUT_IG_TMO) then
     allocate(Lzd%Llr(Lzd%nlr+ndebug))
     do ilr=1,Lzd%nlr
        Lzd%Llr(ilr)=locreg_null()
     end do
     !for now, always true because we want to calculate the hamiltonians for all locregs
     if(.not. Lzd%linear) then
        Lzd%lintyp = 0
        !copy Glr to Llr(1)
        call nullify_locreg_descriptors(Lzd%Llr(1))
        call copy_locreg_descriptors(Lzd%Glr,Lzd%Llr(1),subname)
     else 
        Lzd%lintyp = 1
        ! Assign orbitals to locreg (for LCAO IG each orbitals corresponds to an atomic function. WILL NEED TO CHANGE THIS)
        call assignToLocreg(iproc,nproc,orbs%nspinor,nspin_ig,atoms,orbs,Lzd)

        ! determine the localization regions
        ! calculateBounds indicate whether the arrays with the bounds (for convolutions...) shall also
        ! be allocated and calculated. In principle this is only necessary if the current process has orbitals
        ! in this localization region.
        allocate(calculateBounds(lzd%nlr),stat=i_stat)
        call memocc(i_stat,calculateBounds,'calculateBounds',subname)
        calculateBounds=.true.
!        call determine_locreg_periodic(iproc,Lzd%nlr,rxyz,locrad,hx,hy,hz,Lzd%Glr,Lzd%Llr,calculateBounds)
        call determine_locreg_parallel(iproc,nproc,Lzd%nlr,rxyz,locrad,&
             Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),Lzd%Glr,Lzd%Llr,&
             orbs,calculateBounds)  
        i_all = -product(shape(calculateBounds))*kind(calculateBounds) 
        deallocate(calculateBounds,stat=i_stat)
        call memocc(i_stat,i_all,'calculateBounds',subname)
        i_all = -product(shape(locrad))*kind(locrad)
        deallocate(locrad,stat=i_stat)
        call memocc(i_stat,i_all,'locrad',subname)

        ! determine the wavefunction dimension
        call wavefunction_dimension(Lzd,orbs)
     end if
  else
     Lzd%lintyp = 2
  end if
  
!DEBUG
!!if(iproc==0)then
!!print *,'###################################################'
!!print *,'##        General information:                   ##'
!!print *,'###################################################'
!!print *,'Lzd%nlr,linear, ndimpotisf :',Lzd%nlr,Lzd%linear,Lzd%ndimpotisf
!!print *,'###################################################'
!!print *,'##        Global box information:                ##'
!!print *,'###################################################'
!!write(*,'(a24,3i4)')'Global region n1,n2,n3:',Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3
!!write(*,*)'Global fine grid: nfl',Lzd%Glr%d%nfl1,Lzd%Glr%d%nfl2,Lzd%Glr%d%nfl3
!!write(*,*)'Global fine grid: nfu',Lzd%Glr%d%nfu1,Lzd%Glr%d%nfu2,Lzd%Glr%d%nfu3
!!write(*,*)'Global inter. grid: ni',Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i
!!write(*,'(a27,f6.2,f6.2,f6.2)')'Global dimension (1x,y,z):',Lzd%Glr%d%n1*hx,Lzd%Glr%d%n2*hy,Lzd%Glr%d%n3*hz
!!write(*,'(a17,f12.2)')'Global volume: ',Lzd%Glr%d%n1*hx*Lzd%Glr%d%n2*hy*Lzd%Glr%d%n3*hz
!!print *,'Global wfd statistics:',Lzd%Glr%wfd%nseg_c,Lzd%Glr%wfd%nseg_f,Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f
!!print *,'###################################################'
!!print *,'##        Local boxes information:               ##'
!!print *,'###################################################'
!!do i_stat =1, Lzd%nlr
!!   write(*,*)'=====> Region:',i_stat
!!   write(*,'(a24,3i4)')'Local region n1,n2,n3:',Lzd%Llr(i_stat)%d%n1,Lzd%Llr(i_stat)%d%n2,Lzd%Llr(i_stat)%d%n3
!!   write(*,*)'Local fine grid: nfl',Lzd%Llr(i_stat)%d%nfl1,Lzd%Llr(i_stat)%d%nfl2,Lzd%Llr(i_stat)%d%nfl3
!!   write(*,*)'Local fine grid: nfu',Lzd%Llr(i_stat)%d%nfu1,Lzd%Llr(i_stat)%d%nfu2,Lzd%Llr(i_stat)%d%nfu3
!!   write(*,*)'Local inter. grid: ni',Lzd%Llr(i_stat)%d%n1i,Lzd%Llr(i_stat)%d%n2i,Lzd%Llr(i_stat)%d%n3i
!!   write(*,'(a27,f6.2,f6.2,f6.2)')'Local dimension (1x,y,z):',Lzd%Llr(i_stat)%d%n1*hx,Lzd%Llr(i_stat)%d%n2*hy,&
!!            Lzd%Llr(i_stat)%d%n3*hz
!!   write(*,'(a17,f12.2)')'Local volume: ',Lzd%Llr(i_stat)%d%n1*hx*Lzd%Llr(i_stat)%d%n2*hy*Lzd%Llr(i_stat)%d%n3*hz
!!   print *,'Local wfd statistics:',Lzd%Llr(i_stat)%wfd%nseg_c,Lzd%Llr(i_stat)%wfd%nseg_f,Lzd%Llr(i_stat)%wfd%nvctr_c,&
!!            Lzd%Llr(i_stat)%wfd%nvctr_f
!!end do
!!end if
!!call mpi_finalize(i_stat)
!!stop
!END DEBUG

end subroutine check_linear_and_create_Lzd

subroutine create_LzdLIG(iproc,nproc,nspin,linearmode,hx,hy,hz,Glr,atoms,orbs,rxyz,Lzd)
  use module_base
  use module_types
  use module_xc
  implicit none

  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: Glr
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data),intent(inout) :: orbs
  integer, intent(in) :: linearmode
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  type(local_zone_descriptors), intent(inout) :: Lzd
!  real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
  !Local variables
  character(len=*), parameter :: subname='check_linear_and_create_Lzd'
  logical :: linear
  integer :: iat,ityp,nspin_ig,i_all,i_stat,ilr
  real(gp), dimension(:), allocatable :: locrad
  logical,dimension(:),allocatable :: calculateBounds

  !default variables
  Lzd%nlr = 1

  Lzd%hgrids(1)=hx
  Lzd%hgrids(2)=hy
  Lzd%hgrids(3)=hz

  if (nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=nspin
  end if

  linear  = .true.
  if (linearmode == INPUT_IG_LIG .or. linearmode == INPUT_IG_FULL) then
     Lzd%nlr=atoms%astruct%nat
     allocate(locrad(Lzd%nlr+ndebug),stat=i_stat)
     call memocc(i_stat,locrad,'locrad',subname)
     ! locrad read from last line of  psppar
     do iat=1,atoms%astruct%nat
        ityp = atoms%astruct%iatype(iat)
        locrad(iat) = atoms%rloc(ityp,1)
     end do  
     call timing(iproc,'check_IG      ','ON')
     call check_linear_inputguess(iproc,Lzd%nlr,rxyz,locrad,hx,hy,hz,&
          Glr,linear) 
     call timing(iproc,'check_IG      ','OF')
     if(nspin >= 4) linear = .false. 
  end if

  ! If we are using cubic code : by choice or because locregs are too big
  if (linearmode == INPUT_IG_OFF .or. .not. linear) then
     linear = .false.
     Lzd%nlr = 1
  end if

  Lzd%linear = .true.
  if (.not. linear)  Lzd%linear = .false.

!  print *,'before Glr => Lzd%Glr'
  call nullify_locreg_descriptors(Lzd%Glr)
  call copy_locreg_descriptors(Glr,Lzd%Glr,subname)

  if(linearmode /= INPUT_IG_TMO) then
     allocate(Lzd%Llr(Lzd%nlr+ndebug),stat=i_stat)
     do ilr=1,Lzd%nlr
        Lzd%Llr(ilr)=locreg_null()
     end do
     !for now, always true because we want to calculate the hamiltonians for all locregs

     if(.not. Lzd%linear) then
        Lzd%lintyp = 0
        !copy Glr Lzd%Llr(1)
        call nullify_locreg_descriptors(Lzd%Llr(1))
!        print *,'before Glr => Lzd%Llr(1)'
        call copy_locreg_descriptors(Glr,Lzd%Llr(1),subname)
     else 
        Lzd%lintyp = 1
        ! Assign orbitals to locreg (for LCAO IG each orbitals corresponds to an atomic function. WILL NEED TO CHANGE THIS)
        call assignToLocreg(iproc,nproc,orbs%nspinor,nspin_ig,atoms,orbs,Lzd)

        ! determine the localization regions
        ! calculateBounds indicate whether the arrays with the bounds (for convolutions...) shall also
        ! be allocated and calculated. In principle this is only necessary if the current process has orbitals
        ! in this localization region.
        allocate(calculateBounds(lzd%nlr),stat=i_stat)
        call memocc(i_stat,calculateBounds,'calculateBounds',subname)
        calculateBounds=.true.
!        call determine_locreg_periodic(iproc,Lzd%nlr,rxyz,locrad,hx,hy,hz,Glr,Lzd%Llr,calculateBounds)
        call determine_locreg_parallel(iproc,nproc,Lzd%nlr,rxyz,locrad,&
             hx,hy,hz,Glr,Lzd%Llr,&
             orbs,calculateBounds)  
        i_all = -product(shape(calculateBounds))*kind(calculateBounds) 
        deallocate(calculateBounds,stat=i_stat)
        call memocc(i_stat,i_all,'calculateBounds',subname)
        i_all = -product(shape(locrad))*kind(locrad)
        deallocate(locrad,stat=i_stat)
        call memocc(i_stat,i_all,'locrad',subname)

        ! determine the wavefunction dimension
        call wavefunction_dimension(Lzd,orbs)
     end if
  else
     Lzd%lintyp = 2
  end if

!DEBUG
!!if(iproc==0)then
!!print *,'###################################################'
!!print *,'##        General information:                   ##'
!!print *,'###################################################'
!!print *,'Lzd%nlr,linear, Lpsidimtot, ndimpotisf, Lnprojel:',Lzd%nlr,Lzd%linear,Lzd%ndimpotisf
!!print *,'###################################################'
!!print *,'##        Global box information:                ##'
!!print *,'###################################################'
!!write(*,'(a24,3i4)')'Global region n1,n2,n3:',Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3
!!write(*,*)'Global fine grid: nfl',Lzd%Glr%d%nfl1,Lzd%Glr%d%nfl2,Lzd%Glr%d%nfl3
!!write(*,*)'Global fine grid: nfu',Lzd%Glr%d%nfu1,Lzd%Glr%d%nfu2,Lzd%Glr%d%nfu3
!!write(*,*)'Global inter. grid: ni',Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i
!!write(*,'(a27,f6.2,f6.2,f6.2)')'Global dimension (1x,y,z):',Lzd%Glr%d%n1*hx,Lzd%Glr%d%n2*hy,Lzd%Glr%d%n3*hz
!!write(*,'(a17,f12.2)')'Global volume: ',Lzd%Glr%d%n1*hx*Lzd%Glr%d%n2*hy*Lzd%Glr%d%n3*hz
!!print *,'Global wfd statistics:',Lzd%Glr%wfd%nseg_c,Lzd%Glr%wfd%nseg_f,Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f
!!write(*,'(a17,f12.2)')'Global volume: ',Lzd%Glr%d%n1*input%hx*Lzd%Glr%d%n2*input%hy*Lzd%Glr%d%n3*input%hz
!!print *,'Global wfd statistics:',Lzd%Glr%wfd%nseg_c,Lzd%Glr%wfd%nseg_f,Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f
!!print *,'###################################################'
!!print *,'##        Local boxes information:               ##'
!!print *,'###################################################'
!!do i_stat =1, Lzd%nlr
!!   write(*,*)'=====> Region:',i_stat
!!   write(*,'(a24,3i4)')'Local region n1,n2,n3:',Lzd%Llr(i_stat)%d%n1,Lzd%Llr(i_stat)%d%n2,Lzd%Llr(i_stat)%d%n3
!!   write(*,*)'Local fine grid: nfl',Lzd%Llr(i_stat)%d%nfl1,Lzd%Llr(i_stat)%d%nfl2,Lzd%Llr(i_stat)%d%nfl3
!!   write(*,*)'Local fine grid: nfu',Lzd%Llr(i_stat)%d%nfu1,Lzd%Llr(i_stat)%d%nfu2,Lzd%Llr(i_stat)%d%nfu3
!!   write(*,*)'Local inter. grid: ni',Lzd%Llr(i_stat)%d%n1i,Lzd%Llr(i_stat)%d%n2i,Lzd%Llr(i_stat)%d%n3i
!!   write(*,'(a27,f6.2,f6.2,f6.2)')'Local dimension (1x,y,z):',Lzd%Llr(i_stat)%d%n1*hx,Lzd%Llr(i_stat)%d%n2*hy,&
!!            Lzd%Llr(i_stat)%d%n3*hz
!!   write(*,'(a17,f12.2)')'Local volume: ',Lzd%Llr(i_stat)%d%n1*hx*Lzd%Llr(i_stat)%d%n2*hy*Lzd%Llr(i_stat)%d%n3*hz
!!   print *,'Local wfd statistics:',Lzd%Llr(i_stat)%wfd%nseg_c,Lzd%Llr(i_stat)%wfd%nseg_f,Lzd%Llr(i_stat)%wfd%nvctr_c,&
!!            Lzd%Llr(i_stat)%wfd%nvctr_f
!!end do
!!end if
!call mpi_finalize(i_stat)
!stop
!END DEBUG

end subroutine create_LzdLIG




integer function optimalLength(totalLength, value)
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: totalLength, value
  
  optimalLength=totalLength-ceiling(log10(dble(value+1)+1.d-10))

end function optimalLength

subroutine init_orbitals_data_for_linear(iproc, nproc, nspinor, input, astruct, rxyz, lorbs)
  use module_base
  use module_types
  use module_interfaces, except_this_one => init_orbitals_data_for_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nspinor
  type(input_variables),intent(in) :: input
  type(atomic_structure),intent(in) :: astruct
  real(kind=8),dimension(3,astruct%nat),intent(in) :: rxyz
  type(orbitals_data),intent(out) :: lorbs
  
  ! Local variables
  integer :: norb, norbu, norbd, ityp, iat, ilr, istat, iall, iorb, nlr
  integer,dimension(:),allocatable :: norbsPerLocreg, norbsPerAtom
  real(kind=8),dimension(:,:),allocatable :: locregCenter
  character(len=*),parameter :: subname='init_orbitals_data_for_linear'

  call timing(iproc,'init_orbs_lin ','ON')
  
  call nullify_orbitals_data(lorbs)
 
  ! Count the number of basis functions.
  allocate(norbsPerAtom(astruct%nat), stat=istat)
  call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)
  norb=0
  nlr=0
  do iat=1,astruct%nat
      ityp=astruct%iatype(iat)
      norbsPerAtom(iat)=input%lin%norbsPerType(ityp)
      norb=norb+input%lin%norbsPerType(ityp)
      nlr=nlr+input%lin%norbsPerType(ityp)
  end do

  ! Distribute the basis functions among the processors.
  norbu=norb
  norbd=0
!!$  call orbitals_descriptors_forLinear(iproc, nproc, norb, norbu, norbd, input%nspin, nspinor,&
!!$       input%nkpt, input%kpt, input%wkpt, lorbs)
!!$  call repartitionOrbitals(iproc, nproc, lorbs%norb, lorbs%norb_par,&
!!$       lorbs%norbp, lorbs%isorb_par, lorbs%isorb, lorbs%onWhichMPI)
 
  call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, nspinor,&
       input%gen_nkpt, input%gen_kpt, input%gen_wkpt, lorbs,.true.) !simple repartition

  allocate(locregCenter(3,nlr), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)
  
  ilr=0
  do iat=1,astruct%nat
      ityp=astruct%iatype(iat)
      do iorb=1,input%lin%norbsPerType(ityp)
          ilr=ilr+1
          locregCenter(:,ilr)=rxyz(:,iat)
          ! DEBUGLR write(10,*) iorb,locregCenter(:,ilr)
      end do
  end do
 
  allocate(norbsPerLocreg(nlr), stat=istat)
  call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname)
  norbsPerLocreg=1 !should be norbsPerLocreg
    
  iall=-product(shape(lorbs%inWhichLocreg))*kind(lorbs%inWhichLocreg)
  deallocate(lorbs%inWhichLocreg, stat=istat)
  call memocc(istat, iall, 'lorbs%inWhichLocreg', subname)
  call assignToLocreg2(iproc, nproc, lorbs%norb, lorbs%norb_par, astruct%nat, nlr, &
       input%nspin, norbsPerLocreg, locregCenter, lorbs%inwhichlocreg)

  iall=-product(shape(lorbs%onwhichatom))*kind(lorbs%onwhichatom)
  deallocate(lorbs%onwhichatom, stat=istat)
  call memocc(istat, iall, 'lorbs%onwhichatom', subname)
  call assignToLocreg2(iproc, nproc, lorbs%norb, lorbs%norb_par, astruct%nat, astruct%nat, &
       input%nspin, norbsPerAtom, rxyz, lorbs%onwhichatom)
  
  allocate(lorbs%eval(lorbs%norb), stat=istat)
  call memocc(istat, lorbs%eval, 'lorbs%eval', subname)
  lorbs%eval=-.5d0
  
  iall=-product(shape(norbsPerLocreg))*kind(norbsPerLocreg)
  deallocate(norbsPerLocreg, stat=istat)
  call memocc(istat, iall, 'norbsPerLocreg', subname)
  
  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter, stat=istat)
  call memocc(istat, iall, 'locregCenter', subname)

  iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
  deallocate(norbsPerAtom, stat=istat)
  call memocc(istat, iall, 'norbsPerAtom', subname)


  call timing(iproc,'init_orbs_lin ','OF')

end subroutine init_orbitals_data_for_linear


! initializes locrad and locregcenter
subroutine lzd_init_llr(iproc, nproc, input, astruct, rxyz, orbs, lzd)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(input_variables),intent(in) :: input
  type(atomic_structure),intent(in) :: astruct
  real(kind=8),dimension(3,astruct%nat),intent(in) :: rxyz
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(inout) :: lzd
  
  ! Local variables
  integer :: iat, ityp, ilr, istat, iorb, iall
  real(kind=8),dimension(:,:),allocatable :: locregCenter
  character(len=*),parameter :: subname='lzd_init_llr'
  real(8):: t1, t2

  call timing(iproc,'init_locregs  ','ON')
  t1=mpi_wtime()
  
  nullify(lzd%llr)

  lzd%nlr=orbs%norb

  allocate(locregCenter(3,lzd%nlr), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)
  
  ilr=0
  do iat=1,astruct%nat
      ityp=astruct%iatype(iat)
      do iorb=1,input%lin%norbsPerType(ityp)
          ilr=ilr+1
          locregCenter(:,ilr)=rxyz(:,iat)
      end do
  end do

  ! Allocate the array of localisation regions
  allocate(lzd%Llr(lzd%nlr),stat=istat)
  do ilr=1,lzd%nlr
     lzd%Llr(ilr)=locreg_null()
  end do
  do ilr=1,lzd%nlr
      lzd%llr(ilr)%locrad=input%lin%locrad(ilr)
      lzd%llr(ilr)%locregCenter=locregCenter(:,ilr)
  end do

  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter, stat=istat)
  call memocc(istat, iall, 'locregCenter', subname)
  
  t2=mpi_wtime()
  !if(iproc==0) write(*,*) 'in lzd_init_llr: time',t2-t1
  call timing(iproc,'init_locregs  ','OF')

end subroutine lzd_init_llr


subroutine update_locreg(iproc, nproc, nlr, locrad, locregCenter, glr_tmp, &
           useDerivativeBasisFunctions, nscatterarr, hx, hy, hz, astruct, input, &
           orbs_KS, orbs, lzd, npsidim_orbs, npsidim_comp, lbcomgp, lbcollcom, lfoe, lbcollcom_sr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => update_locreg
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nlr
  integer,intent(out) :: npsidim_orbs, npsidim_comp
  logical,intent(in) :: useDerivativeBasisFunctions
  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(kind=8),intent(in) :: hx, hy, hz
  type(atomic_structure),intent(in) :: astruct
  type(input_variables),intent(in) :: input
  real(kind=8),dimension(nlr),intent(in) :: locrad
  type(orbitals_data),intent(in) :: orbs_KS, orbs
  real(kind=8),dimension(3,nlr),intent(in) :: locregCenter
  type(locreg_descriptors),intent(in) :: glr_tmp
  type(local_zone_descriptors),intent(inout) :: lzd
  type(p2pComms),intent(inout) :: lbcomgp
  type(foe_data),intent(inout),optional :: lfoe
  type(collective_comms),intent(inout) :: lbcollcom
  type(collective_comms),intent(inout),optional :: lbcollcom_sr

  
  ! Local variables
  integer :: iorb, ilr, npsidim, istat
  character(len=*),parameter :: subname='update_locreg'

  call timing(iproc,'updatelocreg1','ON') 
  if (present(lfoe)) call nullify_foe(lfoe)
  call nullify_collective_comms(lbcollcom)
  if (present(lbcollcom_sr)) then
      call nullify_collective_comms(lbcollcom_sr)
  end if
  call nullify_p2pComms(lbcomgp)
  call nullify_local_zone_descriptors(lzd)
  !!tag=1

  ! Allocate the array of localisation regions
  lzd%nlr=nlr
  allocate(lzd%Llr(lzd%nlr),stat=istat)
  do ilr=1,lzd%nlr
     lzd%Llr(ilr)=locreg_null()
  end do
  do ilr=1,lzd%nlr
      lzd%llr(ilr)%locrad=locrad(ilr)
      lzd%llr(ilr)%locregCenter=locregCenter(:,ilr)
  end do
  call timing(iproc,'updatelocreg1','OF') 
  call initLocregs(iproc, nproc, lzd, hx, hy, hz, astruct, orbs, glr_tmp, 's')!, llborbs)
  call timing(iproc,'updatelocreg1','ON') 
  call nullify_locreg_descriptors(lzd%glr)
  call copy_locreg_descriptors(glr_tmp, lzd%glr, subname)
  lzd%hgrids(1)=hx
  lzd%hgrids(2)=hy
  lzd%hgrids(3)=hz

  npsidim = 0
  do iorb=1,orbs%norbp
   ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
   npsidim = npsidim + lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
  end do

  npsidim_orbs=max(npsidim,1)
  ! set npsidim_comp here too?!
  npsidim_comp=1

!  ! don't really want to keep this unless we do so right from the start, but for now keep it to avoid updating refs
!  orbs%eval=-.5d0

  call timing(iproc,'updatelocreg1','OF') 

  if (present(lfoe)) call init_foe(iproc, nproc, lzd, astruct, input, orbs_KS, orbs, lfoe, .false.)

  call init_collective_comms(iproc, nproc, npsidim_orbs, orbs, lzd, lbcollcom)
  if (present(lbcollcom_sr)) then
      call init_collective_comms_sumro(iproc, nproc, lzd, orbs, nscatterarr, lbcollcom_sr)
  end if

  call initialize_communication_potential(iproc, nproc, nscatterarr, orbs, lzd, lbcomgp)
  call allocateCommunicationsBuffersPotential(lbcomgp, subname)

end subroutine update_locreg


subroutine update_ldiis_arrays(tmb, subname, ldiis)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(DFT_wavefunction),intent(in) :: tmb
  character(len=*),intent(in) :: subname
  type(localizedDIISParameters),intent(inout) :: ldiis

  ! Local variables
  integer :: iall, istat, ii, iorb, ilr

  iall=-product(shape(ldiis%phiHist))*kind(ldiis%phiHist)
  deallocate(ldiis%phiHist, stat=istat)
  call memocc(istat, iall, 'ldiis%phiHist', subname)
  iall=-product(shape(ldiis%hphiHist))*kind(ldiis%hphiHist)
  deallocate(ldiis%hphiHist, stat=istat)
  call memocc(istat, iall, 'ldiis%hphiHist', subname)

  ii=0
  do iorb=1,tmb%orbs%norbp
      ilr=tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb)
      ii=ii+ldiis%isx*(tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f)
  end do

  allocate(ldiis%phiHist(ii), stat=istat)
  call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
  allocate(ldiis%hphiHist(ii), stat=istat)
  call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)

end subroutine update_ldiis_arrays


subroutine allocate_auxiliary_basis_function(npsidim, subname, lphi, lhphi)
  use module_base
  implicit none

  ! Calling arguments
  integer,intent(in) :: npsidim
  real(kind=8),dimension(:),pointer,intent(out) :: lphi, lhphi
  character(len=*),intent(in) :: subname

  ! Local variables
  integer :: istat

  allocate(lphi(npsidim), stat=istat)
  call memocc(istat, lphi, 'lphi', subname)
  allocate(lhphi(npsidim), stat=istat)
  call memocc(istat, lhphi, 'lhphi', subname)

  call to_zero(npsidim, lphi(1))
  call to_zero(npsidim, lhphi(1))

end subroutine allocate_auxiliary_basis_function


subroutine deallocate_auxiliary_basis_function(subname, lphi, lhphi)
  use module_base
  implicit none

  ! Calling arguments
  real(kind=8),dimension(:),pointer :: lphi, lhphi
  character(len=*),intent(in) :: subname

  ! Local variables
  integer :: istat, iall

  iall=-product(shape(lphi))*kind(lphi)
  deallocate(lphi, stat=istat)
  call memocc(istat, iall, 'lphi', subname)
  iall=-product(shape(lhphi))*kind(lhphi)
  deallocate(lhphi, stat=istat)
  call memocc(istat, iall, 'lhphi', subname)

end subroutine deallocate_auxiliary_basis_function



subroutine destroy_new_locregs(iproc, nproc, tmb)
  use module_base
  use module_types
  use module_interfaces, except_this_one => destroy_new_locregs
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(inout) :: tmb

  ! Local variables
  character(len=*),parameter :: subname='destroy_new_locregs'

  !!call wait_p2p_communication(iproc, nproc, tmb%comgp)
  call synchronize_onesided_communication(iproc, nproc, tmb%comgp)
 ! call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
  call deallocate_p2pComms(tmb%comgp, subname)

  call deallocate_local_zone_descriptors(tmb%lzd, subname)
  call deallocate_orbitals_data(tmb%orbs, subname)
  call deallocate_foe(tmb%foe_obj, subname)
  !call deallocate_sparseMatrix(tmb%linmat%denskern, subname)
  !call deallocate_sparseMatrix(tmb%linmat%ham, subname)
  !call deallocate_sparseMatrix(tmb%linmat%ovrlp, subname)

  call deallocate_collective_comms(tmb%collcom, subname)
  call deallocate_collective_comms(tmb%collcom_sr, subname)

end subroutine destroy_new_locregs


subroutine destroy_DFT_wavefunction(wfn)
  use module_base
  use module_types
  use module_interfaces, except_this_one => destroy_DFT_wavefunction
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(DFT_wavefunction),intent(inout) :: wfn

  ! Local variables
  integer :: istat, iall
  character(len=*),parameter :: subname='destroy_DFT_wavefunction'

  iall=-product(shape(wfn%psi))*kind(wfn%psi)
  deallocate(wfn%psi, stat=istat)
  call memocc(istat, iall, 'wfn%psi', subname)

  call deallocate_p2pComms(wfn%comgp, subname)
  call deallocate_foe(wfn%foe_obj, subname)
  call deallocate_sparseMatrix(wfn%linmat%denskern, subname)
  call deallocate_sparseMatrix(wfn%linmat%inv_ovrlp, subname)
  call deallocate_sparseMatrix(wfn%linmat%ovrlp, subname)
  call deallocate_sparseMatrix(wfn%linmat%ham, subname)

  call deallocate_orbitals_data(wfn%orbs, subname)
  !call deallocate_communications_arrays(wfn%comms, subname)
  call deallocate_collective_comms(wfn%collcom, subname)
  call deallocate_collective_comms(wfn%collcom_sr, subname)
  call deallocate_local_zone_descriptors(wfn%lzd, subname)

  if (associated(wfn%coeff)) then
      iall=-product(shape(wfn%coeff))*kind(wfn%coeff)
      deallocate(wfn%coeff, stat=istat)
      call memocc(istat, iall, 'wfn%coeff', subname)
  end if

  !call deallocate_p2pComms(wfn%ham_descr%comgp, subname)
  !call deallocate_local_zone_descriptors(wfn%ham_descr%lzd, subname)
  !call deallocate_matrixDescriptors_foe(wfn%ham_descr%mad, subname)
  !call deallocate_collective_comms(wfn%ham_descr%collcom, subname)

end subroutine destroy_DFT_wavefunction


subroutine update_wavefunctions_size(lzd,npsidim_orbs,npsidim_comp,orbs,iproc,nproc)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer, intent(in) :: iproc, nproc
  integer, intent(out) :: npsidim_orbs, npsidim_comp

  ! Local variables
  integer :: npsidim, ilr, iorb
  integer :: nvctr_tot,jproc,istat,ierr,iall
  integer, allocatable, dimension(:) :: ncntt 
  integer, allocatable, dimension(:,:) :: nvctr_par
  character(len = *), parameter :: subname = "update_wavefunctions_size"

  npsidim = 0
  do iorb=1,orbs%norbp
   ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
   npsidim = npsidim + lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
  end do
  npsidim_orbs=max(npsidim,1)


  nvctr_tot = 1
  do iorb=1,orbs%norbp
     ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
     nvctr_tot = max(nvctr_tot,lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f)
  end do
  if (nproc>1) call mpiallred(nvctr_tot, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)

  allocate(nvctr_par(0:nproc-1,1),stat=istat)
  call memocc(istat,nvctr_par,'nvctr_par',subname)

  call kpts_to_procs_via_obj(nproc,1,nvctr_tot,nvctr_par)

  allocate(ncntt(0:nproc-1+ndebug),stat=istat)
  call memocc(istat,ncntt,'ncntt',subname)

  ncntt(:) = 0
  do jproc=0,nproc-1
     ncntt(jproc)=ncntt(jproc)+&
          nvctr_par(jproc,1)*orbs%norbp*orbs%nspinor
  end do

  npsidim_comp=sum(ncntt(0:nproc-1))

  iall=-product(shape(nvctr_par))*kind(nvctr_par)
  deallocate(nvctr_par,stat=istat)
  call memocc(istat,iall,'nvctr_par',subname) 

  iall=-product(shape(ncntt))*kind(ncntt)
  deallocate(ncntt,stat=istat)
  call memocc(istat,iall,'ncntt',subname)  

end subroutine update_wavefunctions_size


subroutine create_large_tmbs(iproc, nproc, KSwfn, tmb, denspot, input, at, rxyz, lowaccur_converged)
  use module_base
  use module_types
  use module_interfaces, except_this_one => create_large_tmbs
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(DFT_Wavefunction),intent(inout):: KSwfn, tmb
  type(DFT_local_fields),intent(in):: denspot
  type(input_variables),intent(in):: input
  type(atoms_data),intent(in):: at
  real(8),dimension(3,at%astruct%nat),intent(in):: rxyz
  logical,intent(in):: lowaccur_converged

  ! Local variables
  integer:: iorb, ilr, istat, iall
  real(8),dimension(:),allocatable:: locrad_tmp
  real(8),dimension(:,:),allocatable:: locregCenter
  character(len=*),parameter:: subname='create_large_tmbs'

  allocate(locregCenter(3,tmb%lzd%nlr), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)
  allocate(locrad_tmp(tmb%lzd%nlr), stat=istat)
  call memocc(istat, locrad_tmp, 'locrad_tmp', subname)

  do iorb=1,tmb%orbs%norb
      ilr=tmb%orbs%inwhichlocreg(iorb)
      locregCenter(:,ilr)=tmb%lzd%llr(ilr)%locregCenter
  end do
  do ilr=1,tmb%lzd%nlr
      locrad_tmp(ilr)=tmb%lzd%llr(ilr)%locrad+8.d0*tmb%lzd%hgrids(1)
  end do

  !temporary,  moved from update_locreg
  tmb%orbs%eval=-0.5_gp
  call update_locreg(iproc, nproc, tmb%lzd%nlr, locrad_tmp, locregCenter, tmb%lzd%glr, &
       .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
       at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%ham_descr%lzd, tmb%ham_descr%npsidim_orbs, tmb%ham_descr%npsidim_comp, &
       tmb%ham_descr%comgp, tmb%ham_descr%collcom)

  call allocate_auxiliary_basis_function(max(tmb%ham_descr%npsidim_comp,tmb%ham_descr%npsidim_orbs), subname, &
       tmb%ham_descr%psi, tmb%hpsi)

  tmb%ham_descr%can_use_transposed=.false.
  nullify(tmb%ham_descr%psit_c)
  nullify(tmb%ham_descr%psit_f)
  allocate(tmb%confdatarr(tmb%orbs%norbp), stat=istat)

  if(.not.lowaccur_converged) then
      call define_confinement_data(tmb%confdatarr,tmb%orbs,rxyz,at,&
           tmb%ham_descr%lzd%hgrids(1),tmb%ham_descr%lzd%hgrids(2),tmb%ham_descr%lzd%hgrids(3),&
           4,input%lin%potentialPrefac_lowaccuracy,tmb%ham_descr%lzd,tmb%orbs%onwhichatom)
  else
      call define_confinement_data(tmb%confdatarr,tmb%orbs,rxyz,at,&
           tmb%ham_descr%lzd%hgrids(1),tmb%ham_descr%lzd%hgrids(2),tmb%ham_descr%lzd%hgrids(3),&
           4,input%lin%potentialPrefac_highaccuracy,tmb%ham_descr%lzd,tmb%orbs%onwhichatom)
  end if

  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter, stat=istat)
  call memocc(istat, iall, 'locregCenter', subname)
  iall=-product(shape(locrad_tmp))*kind(locrad_tmp)
  deallocate(locrad_tmp, stat=istat)
  call memocc(istat, iall, 'locrad_tmp', subname)

end subroutine create_large_tmbs

