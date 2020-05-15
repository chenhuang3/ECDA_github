!{\src2tex{textfont=tt}}
!!****f* ABINIT/qmc_prep_ctqmc
!! NAME
!! qmc_prep_ctqmc
!!
!! FUNCTION
!! Prepare and call the qmc subroutines
!!
!! COPYRIGHT
!! Copyright (C) 1999-2011 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  hu <type(hu_type)>= U interaction
!!  paw_dmft <type(paw_dmft_type)>= DMFT data structure
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawprtvol = drive the amount of writed data.
!!  weiss <type(green_type)>= weiss function
!!
!! OUTPUT
!!  green <type(green_type)>= green function
!!
!! NOTES
!!
!! PARENTS
!!      impurity_solve
!!
!! CHILDREN
!!      add_matlu,compute_levels,copy_green,copy_matlu,ctqmc_printgreen
!!      ctqmcinterface_finalize,ctqmcinterface_init,ctqmcinterface_run
!!      ctqmcinterface_setopts,destroy_green,destroy_matlu,destroy_oper
!!      diag_matlu,diff_matlu,fac_matlu,flush_unit
!!      hybridization_asymptotic_coefficient,identity_matlu,init_green
!!      init_matlu,init_oper,int_fct,inverse_oper,nullify_matlu,occup_green_tau
!!      print_green,print_matlu,printocc_green,printplot_matlu,prod_matlu
!!      rotate_matlu,shift_matlu,slm2ylm_matlu,sym_matlu,testcode_ctqmc,wrtout
!!      xginv,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine qmc_prep_ctqmc(cryst_struc,green,self,hu,paw_dmft,pawang,pawprtvol,weiss)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_xmpi

 use m_pawang, only : pawang_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type,occup_green_tau,print_green,printocc_green,spline_fct,copy_green,init_green,destroy_green,&
& int_fct,greenldacompute_green
 use m_paw_dmft, only : paw_dmft_type
 use m_abilasi,         only : xginv
 use m_oper, only : oper_type,destroy_oper,init_oper,inverse_oper
 use m_self, only : self_type
 use m_matlu, only : matlu_type,sym_matlu, print_matlu, &
& diag_matlu,init_matlu,destroy_matlu,nullify_matlu,rotate_matlu,checkdiag_matlu, &
& copy_matlu, diff_matlu, slm2ylm_matlu, shift_matlu, prod_matlu,fac_matlu,add_matlu,printplot_matlu,identity_matlu
 use m_hu, only : hu_type,rotatevee_hu
 use m_Ctqmc
 use m_CtqmcInterface
 use m_GreenHyb
 use m_data4entropyDMFT
 !use m_self, only : self_type,initialize_self,destroy_self,print_self,rw_self
 use m_io_tools, only : flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'qmc_prep_ctqmc'
 use interfaces_14_hidewrite
 use interfaces_68_dmft, except_this_one => qmc_prep_ctqmc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(inout) :: green  ! MGNAG: This fix the problem with v7[27:29] on nag@petrus
 type(hu_type), intent(in) :: hu(cryst_struc%ntypat)
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type), intent(in) :: pawang
 integer, intent(in) :: pawprtvol
 type(green_type), intent(inout) :: weiss
 type(self_type), intent(in) :: self

!Local variables ------------------------------
 character(len=500) :: message
 character(len=2) :: gtau_iter,iatomnb
 integer :: iatom,ierr,if1,if2,iflavor,iflavor1,iflavor2,ifreq,im,im1,ispinor,ispinor1,isppol,itau,itypat,im2,ispinor2
 integer :: lpawu,master,mbandc,natom,nflavor,nkpt,nspinor,nsppol,nsppol_imp,tndim,ispa,ispb,ima,imb
 integer :: nproc,myproc,spacecomm,opt_diag,opt_nondiag,testcode,testrot,dmft_nwlo,opt_fk,useylm,nomega,opt_rot
 complex(dpc) :: omega_current,integral(2,2),xsum
 real(dp) :: Doccsum,Noise,omega
 real(dp) :: facd,facnd
 real(dp) :: tsec(2)
! arrays
 real(dp), allocatable :: docc(:,:)
 real(dp), allocatable :: gtmp(:,:), levels_ctqmc(:)
 complex(dpc), allocatable :: levels_ctqmc_nd(:,:)
 complex(dpc), allocatable :: hybri_limit(:,:)
 real(dp), allocatable :: gtmp_nd(:,:,:)
 real(dp) :: umod(2,2)
 character(len=4) :: tag_proc
 character(len=30) :: tmpfil
 integer :: unitnb,ier
 complex(dpc), allocatable :: fw1(:,:),gw_tmp(:,:)
 complex(dpc), allocatable :: gw_tmp_nd(:,:,:)
 complex(dpc), allocatable :: fw1_nd(:,:,:)
 complex(dpc), allocatable :: gw1_nd(:,:,:)
 complex(dpc), allocatable :: shift(:)
 integer,parameter :: optdb=0
 !type(coeff2_type), allocatable :: udens_atoms(:)
! Type    -----------------------------------------
 type(coeff2c_type), allocatable :: eigvectmatlu(:,:)
 type(green_type)  :: weiss_for_rot
 type(matlu_type), allocatable :: dmat_diag(:)
 type(matlu_type), allocatable :: matlu1(:)
 type(matlu_type), allocatable :: matlu2(:)
 type(matlu_type), allocatable :: matlu3(:)
 type(matlu_type), allocatable :: matlu4(:)
 type(matlu_type), allocatable :: identity(:)
 type(matlu_type), allocatable :: level_diag(:)
 type(oper_type)  :: energy_level
 !type(self_type) :: self
! type(green_type) :: gw_loc
 type(CtqmcInterface) :: hybrid   !!! WARNING THIS IS A BACKUP PLAN
 type(green_type) :: greenlda
 type(matlu_type), allocatable  :: hybri_coeff(:)
! ************************************************************************
 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 nsppol=paw_dmft%nsppol
 natom=paw_dmft%natom
 nspinor=paw_dmft%nspinor
 greenlda%whichgreen="LDA"

 call init_green(weiss_for_rot,paw_dmft,opt_oper_ksloc=2)
! weiss_for_rot=>weiss
! call init_green(gw_loc,paw_dmft)
 call copy_green(weiss,weiss_for_rot,opt_tw=2)

!=======================================================================
!== Use one QMC solver   ===============================================
!=======================================================================
 write(message,'(2a)') ch10,'  ===  CT-QMC solver === '
 call wrtout(std_out,message,'COLL')

! Initialise for compiler
 omega_current=czero

! Initialise nproc
 nproc=paw_dmft%nproc

! ======================================
! Allocations: diagonalization and eigenvectors
! ======================================
 !ABI_DATATYPE_ALLOCATE(udens_atoms,(natom))
 ABI_DATATYPE_ALLOCATE(eigvectmatlu,(natom,nsppol))
 ABI_DATATYPE_ALLOCATE(dmat_diag,(natom))
 ABI_DATATYPE_ALLOCATE(identity,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,dmat_diag)
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,identity)
 call identity_matlu(identity,natom)
 do iatom=1,cryst_struc%natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     tndim=nspinor*(2*lpawu+1)
     do isppol=1,nsppol
       ABI_ALLOCATE(eigvectmatlu(iatom,isppol)%value,(tndim,tndim))
     end do
     !ABI_ALLOCATE(udens_atoms(iatom)%value,(2*(2*lpawu+1),2*(2*lpawu+1)))
     dmat_diag(iatom)%mat=czero
   end if
 end do

! ___________________________________________________________________________________
!
!  FIRST PART: DIAGONALISATION AND ROTATIONS.
! ___________________________________________________________________________________
! 

! =================================================================
! Choose to diagonalize and how to do it
! =================================================================

! =================================================================
! Impose diago of density matrix
! =================================================================

! =================================================================
! Impose diago of levels and Ylm basis if opt_nondiag=1
! =================================================================
! opt_diag=2 ! 2: diago density matrix (can be used for historical reasons)
! opt_diag=3 ! 3: weiss function is diagonalized: for test, bad and should NOT be used
! opt_diag=1 ! 1: diago the levels (The best choice).

 if(nspinor==2) then
   opt_diag=1    ! because levels are the quantity to diagonalise for CTQMC and NOT dmat
   useylm=1      ! to avoid complex G(tau)
   opt_nondiag=1 ! mandatory for soc.
 end if
 useylm=0
 opt_diag=2
 opt_nondiag=0 ! impose it now for diagonal CTQMC code

 if(useylm==0) then
   write(std_out,*) " Slm basis is used"
 else
   write(std_out,*) " Ylm basis is used for CTQMC"
 end if
 if(opt_diag==1) then
   write(std_out,*) " The atomic levels are diagonalized"
 else if(opt_diag==2) then
   write(std_out,*) " The correlated occupation matrix is diagonalized"
 end if
! if(useylm==1.and.opt_diag/=1) MSG_ERROR("useylm==1 and opt_diag/=0 is not possible")
 if(hu(1)%jpawu_zero.and.nsppol==2) nsppol_imp=2 ! J=0 and nsppol=2
 if(.not.hu(1)%jpawu_zero.or.nsppol/=2) nsppol_imp=1  ! J/=0 ou nsppol=1
! =================================================================
! Compute LDA Green's function to compare to weiss_for_rot (check)
! =================================================================
! call init_green(greenlda,paw_dmft,opt_oper_ksloc=3)
! call greenldacompute_green(cryst_struc,greenlda,pawang,paw_dmft)
!! call copy_green(greenlda,weiss_for_rot,2)

! =================================================================
! Compute atomic levels
! =================================================================
 call init_oper(paw_dmft,energy_level,opt_ksloc=3)

 ! Compute atomic levels in Slm basis
 ! ----------------------------------
 call compute_levels(cryst_struc,energy_level,self%hdc,pawang,paw_dmft)


! =================================================================
! First rotate to Ylm basis the atomic levels
! =================================================================

 if(useylm==1) then 

   ! Rotate from Slm to Ylm the atomic levels
   ! ----------------------------------------
   call slm2ylm_matlu(energy_level%matlu,natom,1,2)

   ! Print atomic energy levels in Ylm basis
   ! --------------------------------
   if(pawprtvol>=3) then
     write(message,'(a,a)') ch10, " == Print Energy levels in Ylm basis"
     call wrtout(std_out,message,'COLL')
     call print_matlu(energy_level%matlu,natom,1)
   end if

 end if ! useylm

! ===========================================================================================
! Start for diagonalization of levels/density matrix according to opt_diag
! ===========================================================================================
 !opt_rot=2 ! do it one time before CTQMC
 opt_rot=1 ! do all the rotations successively on all different quantities.
 if(opt_diag==1.or.opt_diag==0) then


   if(opt_diag==1) then
! =================================================================
! Diagonalize atomic levels
! =================================================================
     ABI_DATATYPE_ALLOCATE(level_diag,(natom))
     call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,level_diag)

     ! Diagonalise atomic levels (opt_real is necessary, because
     ! rotation must be real in order for the occupations and Green's
     ! function to be real
     ! ---------------------------------------------------------------
     call diag_matlu(energy_level%matlu,level_diag,natom,&
&     prtopt=3,eigvectmatlu=eigvectmatlu,nsppol_imp=nsppol_imp,optreal=1)

     if(opt_rot==1) call copy_matlu(level_diag,energy_level%matlu,natom)

! do iatom=1,cryst_struc%natom
!   lpawu=paw_dmft%lpawu(iatom)
!   if(lpawu/=-1) then
!     tndim=nspinor*(2*lpawu+1)
!     do isppol=1,nsppol
!        do im=1,tndim
!        do im1=1,tndim
!        xsum=czero
!        do im2=1,tndim
!          xsum=xsum+ conjg(eigvectmatlu(iatom,isppol)%value(im,im2))*eigvectmatlu(iatom,isppol)%value(im1,im2)
!        enddo
!        write(6,*) im,im1,xsum
!        enddo
!        enddo
!     end do
!   end if
! end do
! write(message,'(a,2x,a,f13.5)') ch10," == Print I before rot"
! call wrtout(std_out,message,'COLL')
! call print_matlu(identity,natom,1,compl=1,opt_exp=2)
! call rotate_matlu(identity,eigvectmatlu,natom,3,1)
! write(message,'(a,2x,a,f13.5)') ch10," == Print I after rot"
! call wrtout(std_out,message,'COLL')
! call print_matlu(identity,natom,1,compl=1,opt_exp=2)

     call destroy_matlu(level_diag,natom)
     call nullify_matlu(level_diag,natom)
     ABI_DATATYPE_DEALLOCATE(level_diag)

     ! Print diagonalized levels
     ! --------------------------
     if(pawprtvol>=3) then
       write(message,'(a,2x,a,f13.5)') ch10,&
&       " == Print Diagonalized Energy levels for Fermi Level=",paw_dmft%fermie
       call wrtout(std_out,message,'COLL')
       call print_matlu(energy_level%matlu,natom,1,compl=1,opt_exp=1)
     else
       write(message,'(a,2x,a,f13.5)') ch10,&
&       " == Energy levels Diagonalized for Fermi Level=",paw_dmft%fermie
       call wrtout(std_out,message,'COLL')
     end if
   end if

 else if(opt_diag==2) then
! =================================================================
! Diagonalizes density matrix and keep eigenvectors in eigvectmatlu
! =================================================================

   ! Print density matrix before diagonalization
   ! -------------------------------------------
   if(pawprtvol>=3) then
     write(message,'(a,2x,a)') ch10,        " == Density Matrix ="
     call wrtout(std_out,message,'COLL')
     !MGNAG: This call is wrong if green has intent(out), now we use intent(inout)
     call print_matlu(green%occup%matlu,natom,1)
   end if

!!  checkstop: we can have two different diagonalisation basis for the up and dn 
!!  but one use the same basis, unless the error is really to large(>0.1)

   ! Diagonalize density matrix 
   ! ---------------------------
   call diag_matlu(green%occup%matlu,dmat_diag,natom,&
&   prtopt=4,eigvectmatlu=eigvectmatlu,nsppol_imp=nsppol_imp,checkstop=.false.)

   ! Print diagonalized density matrix
   ! ----------------------------------
   if(pawprtvol>=3) then
     write(message,'(a,2x,a)') ch10,&
&     " == Diagonalized Density Matrix in the basis used for QMC ="
     call wrtout(std_out,message,'COLL')
     call print_matlu(dmat_diag,natom,1)

     !write(message,'(2a,i3,13x,a)') ch10,'    ==  Rotation of interaction matrix =='
     !call wrtout(std_out,message,'COLL')
   end if
   
   !if (.not.hu(1)%jpawu_zero) & 
   !MSG_WARNING("In qmc_prep_ctqmc J/=0 and rotation matrix not rotated")
!  Rotate interaction.
!   call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms)

 else if(opt_diag==3) then
! =================================================================
! Diagonalizes weiss function
! =================================================================

   ! Diagonalize Weiss function (not recommanded)
   ! --------------------------------------------
   call diag_matlu(weiss_for_rot%oper(1)%matlu,dmat_diag,natom,&
&   prtopt=5,eigvectmatlu=eigvectmatlu,nsppol_imp=nsppol_imp)

   ! Print Weiss
   ! -----------
   if(pawprtvol>=3) then
     write(message,'(a,2x,a,f13.5)') ch10,&
&     " == Print Diagonalized Weiss function for first frequency and Fermi level=",paw_dmft%fermie
     call wrtout(std_out,message,'COLL')
     call print_matlu(dmat_diag,natom,1)
   end if

 end if
! ===========================================================================================
! END Of diagonalization 
! ===========================================================================================

 call flush_unit(std_out)

! ===========================================================================================
! Broadcast matrix of rotation from processor 0 to the other
! In case of degenerate levels, severals rotations are possible. Here we
! choose the rotation of proc 0. It is arbitrary.
! ===========================================================================================
 do iatom=1,cryst_struc%natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     tndim=nspinor*(2*lpawu+1)
     do isppol=1,nsppol
       call xmpi_bcast(eigvectmatlu(iatom,isppol)%value,0,paw_dmft%spacecomm,ier)
     end do
   end if
 end do


     !unitnb=300000+paw_dmft%myproc
     !call int2char4(paw_dmft%myproc,tag_proc)
     !tmpfil = 'eigvectmatluaftermpi'//tag_proc
     !open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
     !do iflavor1=1,14
     !  do iflavor2=1,14
     !    write(unitnb,*) iflavor1,iflavor2,eigvectmatlu(1,1)%value(iflavor1,iflavor2)
     !  enddo
     !enddo

! ===========================================================================================
! Now rotate various quantities in the new basis
! ===========================================================================================

!=======================================================
! Allocate, Compute, and Rotate atomic levels for CTQMC
!=======================================================

   ! If levels not rotated, rotate them
   ! -----------------------------------
 if((opt_diag==2.or.opt_diag==3).and.opt_rot==1) call rotate_matlu(energy_level%matlu,eigvectmatlu,natom,3,1)

   ! Print atomic levels
   ! -------------------
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels after rotation"
   call wrtout(std_out,message,'COLL')
   call print_matlu(energy_level%matlu,natom,1)
 else 
   write(message,'(a,2x,a,f13.5)') ch10," == CT-QMC Energy levels rotated"
   call wrtout(std_out,message,'COLL')
 end if

!====================================================================
! If levels were diagonalized before, then rotate density matrix for
! information.
!====================================================================
 if(opt_diag==1) then
   ABI_DATATYPE_ALLOCATE(matlu1,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu1)
   call copy_matlu(green%occup%matlu,matlu1,natom)

   ! 1) rotate density matrix to Ylm basis
   ! --------------------------------------
   if(useylm==1) then 
     call slm2ylm_matlu(matlu1,natom,1,2)
     if(pawprtvol>=3) then
       write(message,'(a,a)') ch10, " == Print occupations in Ylm basis"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu1,natom,1)
     end if
   end if

   ! 2) rotate density matrix to rotated basis
   ! -------------------------------------------
   if(opt_rot==1.or.opt_rot==2) call rotate_matlu(matlu1,eigvectmatlu,natom,3,1)
   write(message,'(a,2x,a,f13.5)') ch10," == Rotated occupations (for information)"
   call wrtout(std_out,message,'COLL')
   call print_matlu(matlu1,natom,1,compl=1)
   call destroy_matlu(matlu1,natom)
   call nullify_matlu(matlu1,natom)
   ABI_DATATYPE_DEALLOCATE(matlu1)

 end if

 call flush_unit(std_out)


! =================================================================
! Rotate weiss function according to eigenvectors.
! =================================================================
!!!stop
  ! Rotate Weiss function first in Ylm basis 
  ! -------------------------------------------------------------------
 if(useylm==1) then
   write(message,'(a,2x,a)') ch10, " == Rotation of weiss and greenlda in the Ylm Basis="
   call wrtout(std_out,message,'COLL')
   do ifreq=1,paw_dmft%dmft_nwlo
     call slm2ylm_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,1,0)
     ! call slm2ylm_matlu(greenlda%oper(ifreq)%matlu,natom,1,0)
   end do
 end if

 if(pawprtvol>=3) then
   !   write(message,'(a,2x,a,f13.5)') ch10,& ! debug
   !   " == Print weiss for small freq 1 before rot" ! debug
   !   call wrtout(std_out,message,'COLL') ! debug
   !   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1) !  debug

    ! Print Weiss function
    ! --------------------
   write(message,'(a,2x,a,f13.5)') ch10,& ! debug
&  " == Print weiss for 1st freq before rot" ! debug
   call wrtout(std_out,message,'COLL') ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1,compl=1) !  debug
   write(message,'(a,2x,a,f13.5)') ch10,& ! debug
&  " == Print weiss for last freq before rot" ! debug
   call wrtout(std_out,message,'COLL') ! debug
   call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1) !  debug
!    write(message,'(a,2x,a,f13.5)') ch10,& ! debug
!&   " == Print LDA G for 1st freq before rot" ! debug
!    call wrtout(std_out,message,'COLL') ! debug
!    call print_matlu(greenlda%oper(1)%matlu,natom,1,compl=1,opt_exp=2) !  debug
!    write(message,'(a,2x,a,f13.5)') ch10,& ! debug
!&   " == Print LDA G for last freq before rot" ! debug
!    call wrtout(std_out,message,'COLL') ! debug
!    call print_matlu(greenlda%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1,opt_exp=2) !  debug
 end if

 if(opt_diag/=0) then
   ! Rotate Weiss function first in Ylm basis then in the rotated basis.
   ! -------------------------------------------------------------------
   write(message,'(a,2x,a)') ch10, " == Rotation of weiss ="
   call wrtout(std_out,message,'COLL')
   do ifreq=1,paw_dmft%dmft_nwlo
     if(opt_rot==1) call rotate_matlu(weiss_for_rot%oper(ifreq)%matlu,eigvectmatlu,natom,3,1)
!    call checkdiag_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,tol6)
   end do

   call flush_unit(std_out)
   if(pawprtvol>=3) then
     write(message,'(a,2x,a,f13.5)') ch10,& ! debug
&    " == Print weiss for small freq 1 after rot" ! debug
     call wrtout(std_out,message,'COLL') ! debug
     call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1,compl=1) !  debug
     write(message,'(a,2x,a,f13.5)') ch10,&   ! debug
&    " == Print weiss for last freq after rot"   ! debug
     call wrtout(std_out,message,'COLL')   ! debug
     call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1) ! debug
   end if

!   ! Rotate LDA Green's function first in Ylm basis then in the rotated basis and compare to weiss_for_rot
!   ! -----------------------------------------------------------------------------------------------------
!   write(message,'(a,2x,a)') ch10, " == Rotation of greenlda ="
!   call wrtout(std_out,message,'COLL')
!   do ifreq=1,paw_dmft%dmft_nwlo
!     if(opt_rot==1) call rotate_matlu(greenlda%oper(ifreq)%matlu,eigvectmatlu,natom,3,1)
!     call diff_matlu("Weiss_for_rot","greenlda",weiss_for_rot%oper(ifreq)%matlu,greenlda%oper(ifreq)%matlu,natom,1,tol14)
!!    call checkdiag_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,tol6)
!   end do
!   if(pawprtvol>=3) then
!     write(message,'(a,2x,a,f13.5)') ch10,& ! debug
!&    " == Print greenlda for small freq 1 after rot" ! debug
!     call wrtout(std_out,message,'COLL') ! debug
!     call print_matlu(greenlda%oper(1)%matlu,natom,1,compl=1,opt_exp=2) !  debug
!     write(message,'(a,2x,a,f13.5)') ch10,&   ! debug
!&    " == Print greenlda for last freq after rot"   ! debug
!     call wrtout(std_out,message,'COLL')   ! debug
!     call print_matlu(greenlda%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1,opt_exp=2) ! debug
!   end if
!   call flush_unit(std_out)
 end if

! =================================================================
! Compute analytic limit of hybridization and rotate it
! =================================================================
 ABI_DATATYPE_ALLOCATE(hybri_coeff,(paw_dmft%natom))
 call init_matlu(paw_dmft%natom,paw_dmft%nspinor,paw_dmft%nsppol,paw_dmft%lpawu,hybri_coeff)
 !write(6,*)"hybri1",hybri_coeff(1)%mat(1,1,1,1,1),paw_dmft%natom,cryst_struc%natom

 ! Compute analytical C_ij such that F_ij -> C_ij/iw_n
 ! ---------------------------------------
 call hybridization_asymptotic_coefficient(cryst_struc,paw_dmft,pawang,hybri_coeff)
 write(message,'(a,2x,a)') ch10," == Coeff analytical C_ij such that F -> C_ij/iw_n for large frequency"
 call wrtout(std_out,message,'COLL')

 ! Print analytical C_ij (not rotated)
 ! ---------------------------------------
 call print_matlu(hybri_coeff,natom,1)

 ! Rotate analytical C_ij in Ylm basis
 ! ---------------------------------------
 if(useylm==1) call slm2ylm_matlu(hybri_coeff,natom,1,2)
 if(opt_diag/=0)  then

 ! Rotate analytical C_ij in rotated basis
 ! ---------------------------------------
   if(opt_rot==1.or.opt_rot==2) call rotate_matlu(hybri_coeff,eigvectmatlu,natom,3,1)

 ! Print analytical C_ij (rotated)
 ! ---------------------------------------
   write(message,'(a,2x,a)') ch10," == Coeff analytical C_ij such that F -> C_ij/iw_n after rotation"
   call wrtout(std_out,message,'COLL')
   call print_matlu(hybri_coeff,natom,1,compl=1,opt_exp=1)
 end if

! =================================================================
! Check if rotation is properly done.
! =================================================================
 if(3==4) then
   write(message,'(a,2x,a)') ch10,&
&   " == Print  dmat before rot"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup%matlu,natom,1)
   if(useylm==1) call slm2ylm_matlu(green%occup%matlu,natom,1,2)
   if(opt_rot==1) call rotate_matlu(green%occup%matlu,eigvectmatlu,natom,3,1)
   write(message,'(a,2x,a)') ch10,&
&   " == Print  dmat after rot"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup%matlu,natom,1)

   write(message,'(2a)') ch10,' QMC STOP: DEBUG'
   call wrtout(std_out,message,'COLL')
   MSG_ERROR(message)
 end if
! =================================================================
! Check
! =================================================================

! write(message,'(a,2x,a,f13.5)') ch10,&
!&   " == Print weiss for small tau"
! call wrtout(std_out,message,'COLL')
! call print_matlu(weiss%oper(1)%matlu,natom,1)
! write(message,'(a,2x,a,f13.5)') ch10,&
!&   " == Print weiss for large tau"
! call wrtout(std_out,message,'COLL')
! call print_matlu(weiss%oper(paw_dmft%dmft_nwlo)%matlu,natom,1)
! call flush_unit(std_out)
! write(message,'(2a)') ch10,' Check weiss_for_rot(last freq)'
! call wrtout(std_out,message,'COLL')
! call checkdiag_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,tol6,opt=nspinor)
! call flush_unit(std_out)
! write(message,'(2a)') ch10,' Check weiss_for_rot(ifreq=1)'
! call wrtout(std_out,message,'COLL')
! call checkdiag_matlu(weiss_for_rot%oper(1)%matlu,natom,tol6,opt=nspinor)
! call flush_unit(std_out)

 master=0

! =================================================================
! Print out
! =================================================================

! Print Weiss
! -------------
 if(paw_dmft%dmft_prgn==1) then
   call print_green('Weiss_diag',weiss_for_rot,1,paw_dmft,pawprtvol=1,opt_wt=1,opt_decim=1)
 end if

 write(message,'(a,2x,a,f13.5)') ch10,&
& " == Preparing data for CTQMC"
 call wrtout(std_out,message,'COLL')

! Print Rotate Weiss for 1st and last frequencies
! ------------------------------------------------
 if (pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print rotated weiss function for small freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1,compl=1)  ! debug
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print rotated weiss function for largest freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1)  ! debug
 end if

! =================================================================
!  VARIABLES FOR CTQMC TESTS
 testcode = 0 
 testrot  = 0
 opt_fk=0 ! for developpers to check Fourier transform and computes G0(tau)
 opt_fk=1 ! usual case: for real calculations
! =================================================================

! ___________________________________________________________________________________
!
!  SECOND PART : BUILT HYBRIDIZATION FROM G0
! ___________________________________________________________________________________
! 
! ===========================================================================================
! Compute inverse of weiss  and compute hybridization
! ===========================================================================================

! Compute inverse of weiss  for each Frequency
! ----------------------------------------------
 do ifreq=1,paw_dmft%dmft_nwlo
   ABI_DATATYPE_ALLOCATE(matlu1,(natom))
   ABI_DATATYPE_ALLOCATE(matlu2,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu1)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu2)

   call copy_matlu(weiss_for_rot%oper(ifreq)%matlu,matlu1,natom)

   ! Print G_0(iw_n)
   ! ----------------
   if(optdb==1) call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,paw_dmft%omega_lo(ifreq),"go",60000,imre=1)

   ! Compute G_0^-1
   ! -------------------------------------------
   ! if opt_fk=1 or testcode/=0  Do the inversion
   ! if opt_fk=0                 Do not inverse. 
   ! If testcode=2 and opt_fk=0  Do the inversion 
   ! If testcode=1 and opt_fk=0  Do the inversion but no effect, because it will nevertheless be erased
   ! If opt_fk=1                 Do the inversion
   ! -------------------------------------------
   if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"weiss",12000,imre=1)
   if(opt_fk==1.or.testcode/=0) call inverse_oper(weiss_for_rot%oper(ifreq),option=1,prtopt=1)

   ! Print G_0^-1(iw_n)
   ! ----------------
   if(optdb==1) call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,paw_dmft%omega_lo(ifreq),"goinv",70000,imre=1)

   if(pawprtvol>=4.or.ifreq==paw_dmft%dmft_nwlo) then
     if(opt_fk==1.or.testcode/=0) then
      ! Check inversion : do the product
      ! ----------------------------------------------
       call prod_matlu(weiss_for_rot%oper(ifreq)%matlu,matlu1,matlu2,natom)
       write(message,'(a,2x,a,i7)') ch10,&  ! debug
&      " == Print product of  weiss times invers for freq",ifreq
       call wrtout(std_out,message,'COLL')  ! debug
       call print_matlu(matlu2,natom,1)  ! debug
     end if
   end if

   call destroy_matlu(matlu1,natom)
   call nullify_matlu(matlu1,natom)
   call destroy_matlu(matlu2,natom)
   call nullify_matlu(matlu2,natom)
   ABI_DATATYPE_DEALLOCATE(matlu1)
   ABI_DATATYPE_DEALLOCATE(matlu2)
 end do

 ! Copy weiss_for_rot into weiss
 ! -------------------------------
 !call copy_matlu(weiss_for_rot%oper(ifreq)%matlu,weiss%oper(ifreq)%matlu,natom)


 ! Print G_0^-1 for 1st and last frequencies.
 ! -----------------------------------------
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print G_0^-1 for small freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1)  ! debug
   write(message,'(a,2x,a,e18.10,a)') ch10,&   ! debug
&  " == Print G_0^-1 for last freq in the rotated basis (last freq=", paw_dmft%omega_lo(paw_dmft%dmft_nwlo),")"  ! debug
   call wrtout(std_out,message,'COLL')   ! debug
   call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1) ! debug
 end if

! Substract frequency from diagonal part
! ======================================

 ABI_ALLOCATE(shift,(natom))
 do ifreq=1,paw_dmft%dmft_nwlo
   shift(:)=cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)

!  write(5555,'(400e17.4)') paw_dmft%omega_lo(ifreq),((((((weiss_for_rot%oper(ifreq)%matlu(1)%mat&
!  & (im,im1,isppol,ispinor,ispinor1)-cmplx(0.d0,paw_dmft%omega_lo(ifreq),kind=dp)),im=1,2*3+1),&
!&      im1=1,2*3+1),isppol=1,nsppol),ispinor=1,nspinor),ispinor1=1,nspinor)

  ! Compute G_0^-1-iw_n
  ! --------------------
   if(opt_fk==1) call shift_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,shift)

  ! Compute -G_0^-1+iw_n
  ! --------------------
   if(opt_fk==1) call fac_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,-cone)

  ! Print -G_0^-1+iw_n
  ! --------------------
   if(optdb==1) then
     call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,paw_dmft%omega_lo(ifreq),"G0inv_minus_omega",20000,imre=1)
   end if
 end do

 ! Print -G_0^+1-iw_n=(F-levels) for last freq in the rotated basis" 
 ! ------------------------------------------------------------------
 ABI_DEALLOCATE(shift)
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print G_0^-1-iw_n=-(F-levels) for last freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')   ! debug
   call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1) ! debug
 end if
 
! Check numerical limit of F(i_wn)*iw_n (can be used also to compute F )
! ======================================

 if(pawprtvol>=6) then
   ABI_DATATYPE_ALLOCATE(matlu1,(natom))
   ABI_DATATYPE_ALLOCATE(matlu2,(natom))
   ABI_DATATYPE_ALLOCATE(matlu3,(natom))
   ABI_DATATYPE_ALLOCATE(matlu4,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu1)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu2)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu3)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu4)

   write(message,'(a,2x,a,f13.5)') ch10,  " == energy_levels"
   call wrtout(std_out,message,'COLL')   
   call print_matlu(energy_level%matlu,natom,1,opt_exp=2,compl=1) 

   do ifreq=paw_dmft%dmft_nwlo,1,-1 ! necessary to have matlu4 computed for the max frequency and available for all frequency.
   !do ifreq=paw_dmft%dmftqmc_l,1,-1 ! necessary to have matlu4 computed for the max frequency and available for all frequency.
      ! Compute F (substract levels) for max frequency
      ! -----------------------------------------------
     call add_matlu(weiss_for_rot%oper(ifreq)%matlu,energy_level%matlu,matlu1,natom,-1)

      ! Print F(iw_n)=-(G_0^-1-iw_n+levels)  for last frequency.
      ! --------------------------------------------------------
     if(ifreq==paw_dmft%dmft_nwlo.or.ifreq==paw_dmft%dmftqmc_l) then
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
&       " == Print F(iw_n)=-(G_0^-1-iw_n+levels) for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')   
       call print_matlu(matlu1,natom,1,opt_exp=1,compl=1) 
     end if
     if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"Hybridization",10000,imre=1)

      ! Put F in weiss_for_rot -> CTQMC
      ! -------------------------------
     if(opt_rot==2) call rotate_matlu(weiss_for_rot%oper(ifreq)%matlu,eigvectmatlu,natom,3,1)
!   The following line will produce directly the weiss function for the CTQMC code
     !if(opt_fk==1) call copy_matlu(matlu1,weiss_for_rot%oper(ifreq)%matlu,natom) 


      ! Multiply F by frequency
      ! ------------------------
     call copy_matlu(matlu1,matlu2,natom) 
     call fac_matlu(matlu1,natom,cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp))
     if(ifreq==paw_dmft%dmft_nwlo.or.ifreq==paw_dmft%dmftqmc_l) then
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
&       " == Print numerical C_ij = F(iw_n)*iw_n for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')   
       call print_matlu(matlu1,natom,1,opt_exp=1,compl=1) 
     end if
     if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"cij",72800,imre=1)
     !call rotate_matlu(matlu1,eigvectmatlu,natom,3,1)

     if(ifreq==paw_dmft%dmft_nwlo.or.ifreq==paw_dmft%dmftqmc_l) then
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
&       " == Print numerical after back rotation C_ij = F(iw_n)*iw_n for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')   
       call print_matlu(matlu1,natom,1,opt_exp=1,compl=1) 
     end if
     if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"cij_rotated",72900,imre=1)

      ! Built C_ij/iw_n
      ! ------------------------
     call copy_matlu(hybri_coeff,matlu1,natom)
     call fac_matlu(matlu1,natom,1.d0/cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp))
     if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"cij_over_omega",72000)
    ! if(ifreq==paw_dmft%dmft_nwlo) then
    !   write(message,'(a,2x,a,f13.5)') ch10,  " == Print numerical C_ij/iw_n for frequency",paw_dmft%omega_lo(ifreq)
    !   call wrtout(std_out,message,'COLL')   
    !   call print_matlu(matlu1,natom,1,opt_exp=1,compl=1) 
    ! endif

      ! For test: put C_ij/i_wn into weiss_for_rot
      ! --------------------------------------------
     !call copy_matlu(matlu1,weiss_for_rot%oper(ifreq)%matlu,natom,opt_non_diag=1) 

      ! Compute Hybri - C_ij/iw_n
      ! ------------------------
     call add_matlu(matlu2,matlu1,matlu3,natom,-1)

      ! Print Hybri - C_ij/iw_n
      ! ------------------------
     if(optdb==1) call printplot_matlu(matlu3,natom,paw_dmft%omega_lo(ifreq),"hybri_minus_asymp",74000,imre=1)

      ! Multiply (F-C_ij/i_wn) by (iw_n)**2 to find D_ij such that (F-C_ij/i_wn) -> D_ij/(iw_n)^2 only for last frequency.
      ! ------------------------------------------------------------------------------------------------------------------
     call copy_matlu(matlu3,matlu2,natom)
     call fac_matlu(matlu2,natom,cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)**2)
     if(optdb==1) call printplot_matlu(matlu2,natom,paw_dmft%omega_lo(ifreq),"fminuscijtimesw2",75000,imre=1)
     if(ifreq==paw_dmft%dmft_nwlo.or.ifreq==paw_dmft%dmftqmc_l) then
       call copy_matlu(matlu2,matlu4,natom)
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
&       " == Print numerical (F(iw_n)-C_ij/iw_n)%iw_n^2 for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')   
       call print_matlu(matlu4,natom,1) 
     end if

      ! Built C_ij/iw_n+D_ij/(iw_n)^2
      ! ------------------------
     call copy_matlu(matlu4,matlu3,natom,opt_re=1)
     call fac_matlu(matlu3,natom,1.d0/cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)**2)
     call add_matlu(matlu1,matlu3,matlu2,natom,1)
     if(optdb==1) call printplot_matlu(matlu2,natom,paw_dmft%omega_lo(ifreq),"cij_w_plus_dij_w2",72700,imre=1)
      ! For test: put C_ij/i_wn +D_ij/(iw_n)^2 into weiss_for_rot
      ! --------------------------------------------
     !call copy_matlu(matlu2,weiss_for_rot%oper(ifreq)%matlu,natom,opt_non_diag=1) 


   end do

   call destroy_matlu(matlu1,natom)
   call destroy_matlu(matlu2,natom)
   call destroy_matlu(matlu3,natom)
   call destroy_matlu(matlu4,natom)
   call nullify_matlu(matlu1,natom)
   call nullify_matlu(matlu2,natom)
   call nullify_matlu(matlu3,natom)
   call nullify_matlu(matlu4,natom)
   ABI_DATATYPE_DEALLOCATE(matlu1)
   ABI_DATATYPE_DEALLOCATE(matlu2)
   ABI_DATATYPE_DEALLOCATE(matlu3)
   ABI_DATATYPE_DEALLOCATE(matlu4)
 end if

! =========================================================================================
! Start big loop to compute hybridization
! =========================================================================================
 do iatom=1,cryst_struc%natom
   green%ecorr_qmc(iatom)=zero
   itypat=cryst_struc%typat(iatom)
   lpawu=paw_dmft%lpawu(iatom)
   tndim=2*lpawu+1
   if(lpawu/=-1) then

     nflavor=2*(tndim)
     if(testcode>=1) then
       nflavor=2
       if(testcode==2) then
         ispa=1
         ispb=2
         if(nspinor==1) ispb=1
         ima=1
         imb=1
         if(tndim>4) then
           ima=5 ! row
           imb=4 ! column
         end if
       end if
     end if

     !allocate(correl_loc(nflavor,nflavor))
     !ABI_ALLOCATE(f_with_k,(MIN(paw_dmft%dmft_nwli,paw_dmft%dmftqmc_l),nflavor))
     ABI_ALLOCATE(fw1,(paw_dmft%dmft_nwlo,nflavor))
     ABI_ALLOCATE(fw1_nd,(paw_dmft%dmft_nwlo,nflavor,nflavor))
     ABI_ALLOCATE(levels_ctqmc,(nflavor))
     ABI_ALLOCATE(levels_ctqmc_nd,(nflavor,nflavor))
     levels_ctqmc_nd=czero
     ABI_ALLOCATE(hybri_limit,(nflavor,nflavor))
     hybri_limit=czero
     fw1_nd=czero
     fw1=czero
     !allocate(fw2(paw_dmft%dmft_nwli))
! =================================================================
! Compute Hybridization
! =================================================================

     if (testcode==0) then
       iflavor1=0
       iflavor2=0
       do isppol=1,nsppol 
         do ispinor1=1,nspinor
           do ispinor2=1,nspinor
             do im1=1,tndim
               do im2=1,tndim
                 ! first diagonal terms whatever opt_nondiag
                 iflavor1=im1+tndim*(ispinor1-1)+tndim*(isppol-1)
                 iflavor2=im2+tndim*(ispinor2-1)+tndim*(isppol-1)
                 !write(6,*) isppol,ispinor1,ispinor2,im1,im2
                 !write(6,*) iflavor1,iflavor2

                 if ( iflavor1==iflavor2 ) then

!              Do spline of weiss function for im and isppol
!              Construction of fw1
!             =================================================================

                   do ifreq=1,paw_dmft%dmft_nwlo
                     if(opt_fk==1) then
                       fw1(ifreq,iflavor1)= weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1,im1,isppol,ispinor1,ispinor1) 
                     else if (opt_fk==0) then
                       fw1(ifreq,iflavor1)= weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1,im1,isppol,ispinor1,ispinor1) 
                     end if
                   end do
                   fw1_nd(:,iflavor1,iflavor1)=fw1(:,iflavor1)

                   levels_ctqmc(iflavor1)=real(energy_level%matlu(iatom)%mat(im1,im1,isppol,ispinor1,ispinor1),kind=dp)
                   hybri_limit(iflavor1,iflavor1)=hybri_coeff(iatom)%mat(im1,im1,isppol,ispinor1,ispinor1)


                   if(nsppol==1.and.nspinor==1) then
                     !f_with_k(:,iflavor+tndim)=f_with_k(:,iflavor)
                     fw1(:,iflavor1+tndim)=fw1(:,iflavor1)
                     fw1_nd(:,iflavor1+tndim,iflavor1+tndim)=fw1(:,iflavor1)
                     levels_ctqmc(iflavor1+tndim)=levels_ctqmc(iflavor1)
                     hybri_limit(iflavor1+tndim,iflavor1+tndim)=hybri_limit(iflavor1,iflavor1)
                   end if

                 else 

                   do ifreq=1,paw_dmft%dmft_nwlo
                     if(opt_fk==1) then
                       fw1_nd(ifreq,iflavor1,iflavor2)= &
&                       weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)
                     else if (opt_fk==0) then
                       fw1_nd(ifreq,iflavor1,iflavor2)= &
&                       weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)
                     end if

        ! omega=pi*paw_dmft%temp*(two*float(ifreq)-1)
         !write(3333,*) omega,weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)
         !write(4444,*) omega,fw1_nd(ifreq,iflavor1,iflavor2),"#",iflavor1,iflavor2
        ! if(iflavor1/=iflavor2)write(5555,*) omega,imag(fw1_nd(ifreq,iflavor1,iflavor2)),"#",iflavor1,iflavor2

                   end do
                   hybri_limit(iflavor1,iflavor2)=hybri_coeff(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)

        ! write(3333,*)
        ! write(4444,*)
        ! write(5555,*)

                   if(nsppol==1.and.nspinor==1) then
                     fw1_nd(:,iflavor1+tndim,iflavor2+tndim) = fw1_nd(:,iflavor1,iflavor2)
                     hybri_limit(iflavor1+tndim,iflavor2+tndim)=hybri_limit(iflavor1,iflavor2)
                   end if

                 end if

! <  / HACK >
               end do !im2
             end do !im1
           end do  !ispinor2
         end do  !ispinor1
       end do  !isppol
! < HACK >
       ! JB. On 1000 cpus this can not work since all CPU try to open/write the files
       ! Action : Don't print it or check only one cpu does it.
       if(pawprtvol>=10000000) then
         write(message,'(a,2x,a)') ch10,  " == Hybri for all flavors for CTQMC "
         call wrtout(std_out,message,'COLL')
         do iflavor1=1,nflavor
           write(message,'(4x,14(2e14.5,2x))') (hybri_limit(iflavor1,iflavor2),iflavor2=1,nflavor)
           call wrtout(std_out,message,'COLL')
         end do

         open (unit=111,file='Hybri_cijoveromega',status='unknown',form='formatted')
         open (unit=112,file='Hybri',status='unknown',form='formatted')
         do ifreq=1,paw_dmft%dmft_nwlo
!              weiss_for_rot is G_0^-1-iw_n=-(F-levels)
           if(optdb==1) call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,paw_dmft%omega_lo(ifreq),"weissbefore112",30000)
         end do
         do iflavor1=1,nflavor
           do iflavor2=1,nflavor
             do ifreq=1,paw_dmft%dmft_nwlo
               omega=pi*paw_dmft%temp*(two*float(ifreq)-1)
!              fw1_nd is -G_0^+1-iw_n=(F-levels)
               write(111,'(300e16.5)') paw_dmft%omega_lo(ifreq)&
&               ,fw1_nd(ifreq,iflavor1,iflavor2)-hybri_limit(iflavor1,iflavor2)/cmplx(0.d0,paw_dmft%omega_lo(ifreq),kind=dp)
               write(112,'(300e16.5)') paw_dmft%omega_lo(ifreq),fw1_nd(ifreq,iflavor1,iflavor2)
               !write(1111,*) omega,real(fw1_nd(ifreq,iflavor1,iflavor2))
               !write(1112,*) omega,imag(fw1_nd(ifreq,iflavor1,iflavor2))
             end do
             write(111,*)
             write(112,*)
            ! write(1111,*) 
            ! write(1112,*) 
           end do
         end do
         close(111)
         close(112)
       end if
     end if ! testcode
   ! </ HACK >

! ====================================================================================
!  TEST
!  For testing purpose, built ultra simple hybridization (constant in
!  imaginary time or very simple) or extract some part of the calculated hybridization
! ====================================================================================

     if(testcode>=1) then
       dmft_nwlo=paw_dmft%dmft_nwlo
       paw_dmft%dmft_nwlo=paw_dmft%dmftqmc_l
       ABI_ALLOCATE(gw1_nd,(paw_dmft%dmft_nwlo,nflavor,nflavor))
       gw1_nd=czero

       if (testcode==1) then
         call testcode_ctqmc(paw_dmft%dmftqmc_l,fw1_nd,fw1,gtmp_nd,gw_tmp_nd,&
&         levels_ctqmc,hybri_limit,nflavor,1,paw_dmft%temp,testrot,testcode,umod)
       else if (testcode==2) then
         facnd=0.8d0
         facd=1.0d0
         !write(6,*) "fac",facnd,facd
         levels_ctqmc_nd(2,2)   = energy_level%matlu(iatom)%mat(imb,imb,1,ispb,ispb)
         levels_ctqmc_nd(1,1)   = energy_level%matlu(iatom)%mat(ima,ima,1,ispa,ispa)
         levels_ctqmc(2)   = real(energy_level%matlu(iatom)%mat(imb,imb,1,ispb,ispb),kind=dp)
         levels_ctqmc(1)   = real(energy_level%matlu(iatom)%mat(ima,ima,1,ispa,ispa),kind=dp)
         if(opt_diag/=1) then
           levels_ctqmc_nd(1,2)   = energy_level%matlu(iatom)%mat(ima,imb,1,ispa,ispb)
           levels_ctqmc_nd(2,1)   = energy_level%matlu(iatom)%mat(imb,ima,1,ispb,ispa)
         end if 
         hybri_limit(1,1)  = facd*hybri_coeff(iatom)%mat(ima,ima,1,ispa,ispa)
         hybri_limit(2,2)  = facd*hybri_coeff(iatom)%mat(imb,imb,1,ispb,ispb)
         hybri_limit(1,2)  = facnd*hybri_coeff(iatom)%mat(ima,imb,1,ispa,ispb)
         hybri_limit(2,1)  = facnd*hybri_coeff(iatom)%mat(imb,ima,1,ispb,ispa)
         !write(6,*) "hybri_limit",hybri_limit
         !write(6,*) "levels_ctqmc",levels_ctqmc
         umod=zero

         tmpfil = 'fw1_nd_re'
         open (unit=777,file=trim(tmpfil),status='unknown',form='formatted')
         tmpfil = 'fw1_nd_im'
         open (unit=888,file=trim(tmpfil),status='unknown',form='formatted')
         write(std_out,*) "testcode==2",ispa,ispb,ima,imb
         write(std_out,*) "opt_fk==",opt_fk
         do ifreq=1,paw_dmft%dmftqmc_l
           if (opt_fk==1) then
             fw1_nd(ifreq,1,1) = facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,ima,1,ispa,ispa)
             fw1_nd(ifreq,2,2) = facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,imb,1,ispb,ispb) 
             !fw1_nd(ifreq,1,2) =  weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,imb,1,ispa,ispb)
             !fw1_nd(ifreq,2,1) =  weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,ima,1,ispb,ispa) 
             fw1_nd(ifreq,1,2) = facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,imb,1,ispa,ispb)
             fw1_nd(ifreq,2,1) = facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,ima,1,ispb,ispa) 
             omega=pi*paw_dmft%temp*(two*float(ifreq)-1)
           else if (opt_fk==0) then
             fw1_nd(ifreq,1,1) =  facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,ima,1,ispa,ispa) 
             fw1_nd(ifreq,2,2) =  facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,imb,1,ispb,ispb) 
             fw1_nd(ifreq,1,2) =  facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,imb,1,ispa,ispb)
             fw1_nd(ifreq,2,1) =  facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,ima,1,ispb,ispa) 
             call xginv(fw1_nd(ifreq,:,:),2)
           end if
         end do
         close(777)
         close(888)
       end if

       do if1=1,2
         do if2=1,2
           do ifreq=1,paw_dmft%dmftqmc_l
             omega=pi*paw_dmft%temp*(two*float(ifreq)-1)
             if(if1==if2) then
               gw1_nd(ifreq,if1,if2) =  (cmplx(0.d0,omega,kind=dp)-fw1_nd(ifreq,if1,if2))
             else 
               gw1_nd(ifreq,if1,if2) =  (-fw1_nd(ifreq,if1,if2))
             end if
           end do
         end do
       end do
       do ifreq=1,paw_dmft%dmftqmc_l
         call xginv(gw1_nd(ifreq,:,:),2)
       end do
       write(std_out,*) " testctqmc high frequency limit of hybridization",fw1_nd(paw_dmft%dmftqmc_l,:,:)

       do if1=1,2
         do if2=1,2
           do ifreq=1,paw_dmft%dmftqmc_l
             omega=pi*paw_dmft%temp*(two*float(ifreq)-1)
             write(999,*) omega,gw1_nd(ifreq,if1,if2)
           end do
           write(999,*) 
           call int_fct(gw1_nd(:,if1,if2),(if1==if2),2,paw_dmft,integral(if1,if2))  ! test_1
          ! write(std_out,*) "testctqmc occupations of input Green's function",(integral(if1,if2))
         end do
       end do
       do if1=1,2
         do if2=1,2
          ! write(std_out,*) "testctqmc occupations of input Green's function",(integral(if1,if2)+conjg(integral(if2,if1)))/two
         end do
       end do
       write(std_out,*) "Occupation of model in matrix form"
       do if1=1,2
         write(std_out,'(2(2f13.5,3x))') ((integral(if1,if2)+conjg(integral(if2,if1)))/two,if2=1,2)
       end do
       write(std_out,*) "Limit of hybridization "
       do if1=1,2
         write(std_out,'(2(2f13.5,3x))') (hybri_limit(if1,if2),if2=1,2)
       end do

       ABI_DEALLOCATE(gw1_nd)
       paw_dmft%dmft_nwlo=dmft_nwlo
     end if

     call flush_unit(std_out)
! =================================================================

! ___________________________________________________________________________________
! 
!  THIRD PART : CALL CTQMC
! ___________________________________________________________________________________

! =================================================================
!    Main calls to CTQMC code
! =================================================================
     write(message,'(a,2x,a)') ch10,&
&     " == Initializing CTQMC"
     call wrtout(std_out,message,'COLL')

!    Initialisation
! =================================================================
     nomega=paw_dmft%dmftqmc_l
     call CtqmcInterface_init(hybrid,paw_dmft%dmftqmc_seed,paw_dmft%dmftqmc_n, &
&     paw_dmft%dmftqmc_therm, paw_dmft%dmftctqmc_meas,nflavor,paw_dmft%dmftqmc_l,one/paw_dmft%temp,zero,&
&     std_out,paw_dmft%spacecomm)

!  for non diagonal code
!     call CtqmcInterface_init(hybrid,paw_dmft%dmftqmc_seed,paw_dmft%dmftqmc_n, &
!&         paw_dmft%dmftqmc_therm, paw_dmft%dmftctqmc_meas,nflavor,&
!&         paw_dmft%dmftqmc_l,one/paw_dmft%temp,zero,&
!&         std_out,paw_dmft%spacecomm,opt_nondiag) 

!    options
! =================================================================
     call CtqmcInterface_setOpts(hybrid,&
     opt_Fk      =opt_fk,& 
&     opt_order   =paw_dmft%dmftctqmc_order ,&
&     opt_movie   =paw_dmft%dmftctqmc_mov   ,&
&     opt_analysis=paw_dmft%dmftctqmc_correl,&
&     opt_check   =paw_dmft%dmftctqmc_check ,&
&     opt_noise   =paw_dmft%dmftctqmc_grnns ,&
&     opt_spectra =paw_dmft%dmftctqmc_mrka  ,&
&     opt_gmove   =paw_dmft%dmftctqmc_gmove )  
     write(message,'(a,2x,2a)') ch10,&
&     " == Initialization CTQMC done", ch10
     call wrtout(std_out,message,'COLL')

     ABI_ALLOCATE(gw_tmp,(paw_dmft%dmft_nwlo,nflavor+1))
     ABI_ALLOCATE(gw_tmp_nd,(paw_dmft%dmft_nwlo,nflavor,nflavor+1))

!     use  gw_tmp to put freq
     do ifreq=1,paw_dmft%dmft_nwlo
       gw_tmp(ifreq,nflavor+1)=cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
       gw_tmp_nd(ifreq,nflavor,nflavor+1)=cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
     end do
     ABI_ALLOCATE(gtmp,(paw_dmft%dmftqmc_l,nflavor))
     ! THIS IS A BACKUP PLAN. USING paw_dmft%hybrid makes a segfault on TIKAL
     ! PSC with MPI only (and max2_open64). paw_dmf%hybrid is corrupted
     ! somewhere but I could not find the place in all DMFT routines
     ABI_ALLOCATE(gtmp_nd,(paw_dmft%dmftqmc_l,nflavor,nflavor))
     call flush_unit(std_out)

     if(testcode==0) then
      !unitnb=100000+paw_dmft%myproc
      !call int2char4(paw_dmft%myproc,tag_proc)
      !tmpfil = 'hybrilimit'//tag_proc
      !open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
      !do iflavor1=1,nflavor
      !  do iflavor2=1,nflavor
      !    write(unitnb,*) iflavor1,iflavor2,hybri_limit(iflavor1,iflavor2)
      !  enddo
      !enddo

! =================================================================
!    CTQMC run
! =================================================================
     ABI_ALLOCATE(docc,(1:nflavor,1:nflavor))
     docc(:,:) = zero
       call CtqmcInterface_run(hybrid,fw1(1:paw_dmft%dmftqmc_l,:),Gtau=gtmp,&
         &       Gw=gw_tmp,D=docc(:,:),E=green%ecorr_qmc(iatom),&
&       matU=hu(itypat)%udens,opt_levels=levels_ctqmc)
     call data4entropyDMFT_setDocc(paw_dmft%forentropyDMFT,iatom,docc)
     ABI_DEALLOCATE(docc)

!  for non diagonal code:
!       call CtqmcInterface_run(hybrid,fw1_nd(1:paw_dmft%dmftqmc_l,:,:),Gtau=gtmp_nd,&
!&       Gw=gw_tmp_nd,D=Doccsum,E=green%ecorr_qmc(iatom),&
!&       Noise=Noise,matU=hu(itypat)%udens,opt_levels=levels_ctqmc,hybri_limit=hybri_limit)

     else if (testcode>=1) then

       write(std_out,*) "testcode=",testcode
! =================================================================
!    CTQMC run for tests
! =================================================================
       write(std_out,*) "nomega,dmftqmc_l",nomega,paw_dmft%dmftqmc_l
       call CtqmcInterface_run(hybrid,fw1(1:nomega,:),Gtau=gtmp,&
&       Gw=gw_tmp,E=green%ecorr_qmc(iatom),&
&       matU=umod,opt_levels=levels_ctqmc)

! for non diagonal code
!       call CtqmcInterface_run(hybrid,fw1_nd(1:nomega,:,:),Gtau=gtmp_nd,&
!&       Gw=gw_tmp_nd,D=Doccsum,E=green%ecorr_qmc(iatom),&
!&       Noise=Noise,matU=umod,opt_levels=levels_ctqmc,hybri_limit=hybri_limit)

     end if

! =================================================================
!  TEST
!  If test of the code is activated, and testrot =1 rotate back green's function
!  and stop the code.
! =================================================================
     if(testcode==1) then

       call testcode_ctqmc(paw_dmft%dmftqmc_l,fw1_nd,fw1,gtmp_nd,gw_tmp_nd,&
&       levels_ctqmc,hybri_limit,nflavor,2,paw_dmft%temp,testrot,testcode,umod)

     end if
     if(testcode==2) then
       write(message,'(2a)') ch10,' testcode 2 end of test calculation'
       MSG_ERROR(message)
     end if
! =================================================================
     ABI_DEALLOCATE(fw1)
     ABI_DEALLOCATE(fw1_nd)

     !----------------------------------------
     ! <DEBUG>
     !----------------------------------------
     ! Construct UNIT
     if(paw_dmft%idmftloop < 10) then
       write(gtau_iter,'("0",i1)') paw_dmft%idmftloop
     elseif(paw_dmft%idmftloop >= 10 .and. paw_dmft%idmftloop < 100) then
       write(gtau_iter,'(i2)') paw_dmft%idmftloop
     else
       gtau_iter="xx"
     end if
     if(iatom < 10) then
       write(iatomnb,'("0",i1)') iatom
     elseif(iatom >= 10 .and. iatom < 100) then
       write(iatomnb,'(i2)') iatom
     else
       iatomnb='xx'
     end if

     if(paw_dmft%myproc .eq. mod(nproc+1,nproc)) then
! < HACK >
       open(unit=4242, file=trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_Gtau_"//gtau_iter//".dat")
       call Ctqmc_printGreen(paw_dmft%hybrid(iatom)%hybrid,4242)
       close(4242)
       !open(unit=4243, file=trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_F_"//gtau_iter//".dat")
       !call BathOperator_printF(paw_dmft%hybrid(iatom)%hybrid%bath,4243)
       !close(4243)
       open(unit=4242, file=trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_Gw_"//gtau_iter//".dat")
       do ifreq=1,paw_dmft%dmft_nwlo
         write(4242,'(6f21.14)') (/ (gw_tmp_nd(ifreq,iflavor,iflavor), iflavor=1, nflavor) /) 
       end do
       close(4242)
! </ HACK >
     end if
     !----------------------------------------
     ! </DEBUG>
     !----------------------------------------
     write(message,'(a,2x,a)') ch10,&
&     " == Destroy CTQMC"
     call wrtout(std_out,message,'COLL')
     call CtqmcInterface_finalize(hybrid)
     write(message,'(a,2x,a)') ch10,&
&     " == Destroy CTQMC done"
     call wrtout(std_out,message,'COLL')
     !ABI_DEALLOCATE(f_with_k)
     ABI_DEALLOCATE(hybri_limit)
     ABI_DEALLOCATE(levels_ctqmc_nd)
     ABI_DEALLOCATE(levels_ctqmc)

! ___________________________________________________________________________________
!
!  FOURTH PART : USE OUTPUT OF CTQMC AND DO BACK ROTATION
! ___________________________________________________________________________________
! 

     do itau=1,paw_dmft%dmftqmc_l
       green%oper_tau(itau)%matlu(iatom)%mat(:,:,:,:,:)=czero
     end do
     green%occup_tau%matlu(iatom)%mat(:,:,:,:,:)=czero

     do ifreq=1,paw_dmft%dmft_nwlo
       green%oper(ifreq)%matlu(iatom)%mat(:,:,:,:,:)=czero
     end do
     green%occup%matlu(iatom)%mat(:,:,:,:,:)=czero

!   built time and frequency green's function from output of CTQMC
! =================================================================
     if(opt_nondiag==1) then
       do isppol=1,nsppol
         do ispinor1=1,nspinor
           do im1=1,tndim
             iflavor1=im1+tndim*(ispinor1-1)+tndim*(isppol-1)
             do ispinor2=1,nspinor
               do im2=1,tndim
                 iflavor2=im2+tndim*(ispinor2-1)+tndim*(isppol-1)
                  ! NNtodo here: nsppol=1 above, only, with symetrisation
                  ! automatic
!               iflavor1=im+(isppol-1)*tndim
                 do itau=1,paw_dmft%dmftqmc_l
                   green%oper_tau(itau)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
&                   gtmp_nd(itau,iflavor1,iflavor2)
                   if(nsppol==1.and.nspinor==1) then
                     green%oper_tau(itau)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
&                     (gtmp_nd(itau,iflavor1,iflavor2)+gtmp_nd(itau,iflavor1+tndim,iflavor2+tndim))/two
                   end if
                  ! NNtodo here: isppol above should be one and symetrized
                  ! gtmp
                 end do  !itau
!               ifreq2=0
                 do ifreq=1,paw_dmft%dmft_nwlo
!                 if(paw_dmft%select_log(ifreq)==1) then
!                   ifreq2=ifreq2+1
                   green%oper(ifreq)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
&                   gw_tmp_nd(ifreq,iflavor1,iflavor2)
                   if(nsppol==1.and.nspinor==1) then
                     green%oper(ifreq)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
&                     (gw_tmp_nd(ifreq,iflavor1,iflavor2)+&
&                     gw_tmp_nd(ifreq,iflavor1+tndim,iflavor2+tndim))/two
                   end if
                  ! NNtodo here: isppol above should be one and symetrized
                  ! gw_tmp
!                 endif
                 end do ! ifreq
               end do  ! im2
             end do  ! ispinor2
           end do  ! im1
         end do  ! ispinor
       end do ! isppol
     else
       iflavor=0
       do isppol=1,nsppol
         do ispinor=1,nspinor
           do im=1,tndim
              ! NNtodo here: nsppol=1 above, only, with symetrisation
              ! automatic
!           iflavor=im+(isppol-1)*tndim
             iflavor=iflavor+1
             do itau=1,paw_dmft%dmftqmc_l
               green%oper_tau(itau)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=gtmp(itau,iflavor)
               if(nsppol==1.and.nspinor==1) then
                 green%oper_tau(itau)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=&
&                 (gtmp(itau,iflavor)+gtmp(itau,iflavor+tndim))/two
               end if
              ! NNtodo here: isppol above should be one and symetrized
              ! gtmp
             end do
!           ifreq2=0
             do ifreq=1,paw_dmft%dmft_nwlo
!             if(paw_dmft%select_log(ifreq)==1) then
!               ifreq2=ifreq2+1
               green%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=gw_tmp(ifreq,iflavor)
               if(nsppol==1.and.nspinor==1) then
                 green%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=&
&                 (gw_tmp(ifreq,iflavor)+gw_tmp(ifreq,iflavor+tndim))/two
               end if
              ! NNtodo here: isppol above should be one and symetrized
              ! gw_tmp
!             endif
             end do
           end do
         end do
       end do
     end if
     ABI_DEALLOCATE(gw_tmp)
     ABI_DEALLOCATE(gw_tmp_nd)
     ABI_DEALLOCATE(gtmp)
     ABI_DEALLOCATE(gtmp_nd)
     if(nsppol==1.and.nspinor==1) then
       write(message,'(a,2x,a,f13.5)') ch10,& 
&       " == nsppol==1 and nspden==1: Green functions from CTQMC have been symetrized over spin"
       call wrtout(std_out,message,'COLL')  
     end if
     !write(message,'(i3,4x,2e21.14)') 5,weiss_for_rot%oper(1)%matlu(1)%mat(1,1,1,1,1)
     !call wrtout(std_out,message,'COLL')  ! debug
!     do im=1,tndim
!       do itau=1,paw_dmft%dmftqmc_l
!         gtt=(green%oper_tau(itau)%matlu(iatom)%mat(im,im,1,1,1)+&
!&         green%oper_tau(itau)%matlu(iatom)%mat(im,im,2,1,1))/two
!         green%oper_tau(itau)%matlu(iatom)%mat(im,im,1,1,1)=gtt
!         green%oper_tau(itau)%matlu(iatom)%mat(im,im,2,1,1)=gtt
!       enddo
!     enddo
!     write(6,*)" SYMETRISATION OVER ISPPOL"
     !deallocate(correl_loc,f_with_k,gtmp,fw1,fw2) 
   end if
 end do ! iatom


 if(paw_dmft%dmft_prgn==1) then
   call print_green('QMC_diag_notsym',green,1,paw_dmft,pawprtvol=1,opt_wt=2)
   call print_green('QMC_diag_notsym',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
 end if
 !write(message,'(i3,4x,2e21.14)') 6,weiss_for_rot%oper(1)%matlu(1)%mat(1,1,1,1,1)
 !call wrtout(std_out,message,'COLL')  ! debug
! =================================================================
! Inverse Weiss, then
! Copy Weiss_for_rot into weiss and rotate back weiss to the original basis
! =================================================================

! ABI_ALLOCATE(shift,(natom))
! do ifreq=1,paw_dmft%dmft_nwlo
!  ! First weiss_for_rot contains -G_0^-1+iw_n
!  ! -------------------------------------------
!  ! Compute G_0^-1-iw_n
!  ! --------------------
!       write(6,*) "1"
!  if(opt_fk==1) call fac_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,-cone)
!
!
!       write(6,*) "2"
!  ! Compute G_0^-1
!  ! --------------------
!  shift(:)=cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
!  if(opt_fk==1) call shift_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,shift,signe=1)
!
!       write(6,*) "3"
!  ! Compute G_0
!  ! --------------------
!   call inverse_oper(weiss_for_rot%oper(ifreq),option=1,prtopt=1)
!   ! No need to copy if weiss_for_rot is a pointer to weiss ...
!!   if(useylm==1) call slm2ylm_matlu(weiss%oper(ifreq)%matlu,natom,2,0)
!!   if(opt_diag/=0) call rotate_matlu(weiss%oper(ifreq)%matlu,eigvectmatlu,natom,3,0)
!
!  ! Compute G_0 in the original basis
!  ! --------------------
!   call rotate_matlu(weiss_for_rot%oper(ifreq)%matlu,eigvectmatlu,natom,3,0)
! end do
! ABI_DEALLOCATE(shift)

! =================================================================
! Here compute Self energy from Dyson and print it
! Warning : Weiss_for_rot is inversed inside dyson
! =================================================================
! call initialize_self(self,paw_dmft)
! call dyson(green,paw_dmft,self,weiss_for_rot,opt_weissself=2)
! call rw_self(self,mpi_enreg,paw_dmft,prtopt=2,opt_rw=2,opt_char="diag")
! call destroy_self(self)
 !write(message,'(i3,4x,2e21.14)') 7,weiss%oper(1)%matlu(1)%mat(1,1,1,1,1)
 !call wrtout(std_out,message,'COLL')  ! debug

! =================================================================
! Rotate back green function to original basis (non-diagonal)
!  (and Weiss for further use: might be useful if an back Fourier
!     transformation is done).  
! =================================================================
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print green function for small tau after CTQMC"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper_tau(1)%matlu,natom,1)  ! debug
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print green function for small freq after CTQMC"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper(1)%matlu,natom,1)  ! debug
 end if

 if(pawprtvol>=3) then
!  === Compute non rotated Occupations in green%occup_tau
   call occup_green_tau(green)
   write(message,'(a,2x,a,f13.5)') ch10," == Occupation from G(tau) in the original basis" 
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup_tau%matlu,natom,1)
 end if 

 write(message,'(a,2x,a,f13.5)') ch10,&  
& " == Rotate Green function to original basis "
 call wrtout(std_out,message,'COLL') 
 !write(message,'(i3,4x,2e21.14)') 8,weiss%oper(1)%matlu(1)%mat(1,1,1,1,1)
 !call wrtout(std_out,message,'COLL')  ! debug
 do itau=1,paw_dmft%dmftqmc_l
   if(useylm==1) call slm2ylm_matlu(green%oper_tau(itau)%matlu,natom,2,0)
   if(opt_diag/=0) call rotate_matlu(green%oper_tau(itau)%matlu,eigvectmatlu,natom,3,0)
 end do
 !write(message,'(i3,4x,2e21.14)') 9,weiss%oper(1)%matlu(1)%mat(1,1,1,1,1)
 !call wrtout(std_out,message,'COLL')  ! debug
 do ifreq=1,paw_dmft%dmft_nwlo
   if(useylm==1) call slm2ylm_matlu(green%oper(ifreq)%matlu,natom,2,0)
   if(opt_diag/=0) call rotate_matlu(green%oper(ifreq)%matlu,eigvectmatlu,natom,3,0)
 end do
 !write(message,'(i3,4x,2e21.14)') 10,weiss%oper(1)%matlu(1)%mat(1,1,1,1,1)
 !call wrtout(std_out,message,'COLL')  ! debug

 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&                  ! debug
&  " == Print green function for small time after rotation (in the original basis)" ! debug
   call wrtout(std_out,message,'COLL')                  ! debug
   call print_matlu(green%oper_tau(1)%matlu,natom,1)  ! debug
   write(message,'(a,2x,a,f13.5)') ch10,&                  ! debug
&  " == Print green function for small freq after rotation (in the original basis)" ! debug
   call wrtout(std_out,message,'COLL')                  ! debug
   call print_matlu(green%oper(1)%matlu,natom,1)  ! debug
   !< HACK >
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print diagonalized weiss_for_rot function after rotation for small freq in the original basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1)  ! debug
   !</ HACK >
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print weiss function for small freq in the original basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss%oper(1)%matlu,natom,1)  ! debug
 end if

!  === Compute rotated Occupations in green%occup_tau
 call occup_green_tau(green)

 
 ABI_DATATYPE_ALLOCATE(matlu1,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu1)
 call copy_matlu(green%occup_tau%matlu,matlu1,natom)
 call sym_matlu(cryst_struc,matlu1,pawang)

 write(message,'(a,2x,a,f13.5)') ch10," == Occupation from G(tau) in the original basis" 
 call wrtout(std_out,message,'COLL')
 call print_matlu(green%occup_tau%matlu,natom,1)

 write(message,'(a,2x,a,f13.5)') ch10," == Symetrized occupations"
 call wrtout(std_out,message,'COLL')
 call print_matlu(matlu1,natom,1)

 call diff_matlu("CTQMC Occup","CTQMC Occup symetrized",green%occup_tau%matlu,matlu1,natom,0,tol4,ierr)
 call destroy_matlu(matlu1,natom)
 call nullify_matlu(matlu1,natom)
 ABI_DATATYPE_DEALLOCATE(matlu1)

! =================================================================
! Symetrise green function G(tau) and G(ifreq) to recover symetry 
! artificially broken by QMC
! =================================================================
 write(message,'(a,2x,a,f13.5)') ch10,&  
& " == Symetrise green function after QMC "
 call wrtout(std_out,message,'COLL')  
 do itau=1,paw_dmft%dmftqmc_l
   call sym_matlu(cryst_struc,green%oper_tau(itau)%matlu,pawang)
 end do
 do ifreq=1,paw_dmft%dmft_nwlo
   call sym_matlu(cryst_struc,green%oper(ifreq)%matlu,pawang)
 end do
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print green function for small time after symetrisation"  !  debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper_tau(1)%matlu,natom,1)  ! debug
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print green function for small freq after symetrisation"  !  debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper(1)%matlu,natom,1)  ! debug
 end if
 if(paw_dmft%dmft_prgn==1) then
   call print_green('QMC_sym',green,1,paw_dmft,pawprtvol=1,opt_wt=2)
   call print_green('QMC_sym',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
 end if

!  === Compute Occupations  (Symetrized from oper_tau)
 call occup_green_tau(green)


!  === Print occupations
 call printocc_green(green,6,paw_dmft,3)

 call destroy_oper(energy_level)
 call destroy_matlu(dmat_diag,natom)
 call nullify_matlu(dmat_diag,natom)
 ABI_DATATYPE_DEALLOCATE(dmat_diag)
 call destroy_matlu(identity,natom)
 call nullify_matlu(identity,natom)
 ABI_DATATYPE_DEALLOCATE(identity)
 do iatom=1,cryst_struc%natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     do isppol=1,nsppol
       ABI_DEALLOCATE(eigvectmatlu(iatom,isppol)%value)
       !ABI_DEALLOCATE(udens_atoms(iatom))
     end do
   end if
 end do
 !ABI_DATATYPE_DEALLOCATE(udens_atoms)
 ABI_DATATYPE_DEALLOCATE(eigvectmatlu)
 call destroy_green(weiss_for_rot)
! call destroy_green(gw_loc)
! call destroy_green(greenlda)

!  destroy limit of hybridization
 call destroy_matlu(hybri_coeff,paw_dmft%natom)
 call nullify_matlu(hybri_coeff,natom)
 ABI_DATATYPE_DEALLOCATE(hybri_coeff)

end subroutine qmc_prep_ctqmc
!!***
