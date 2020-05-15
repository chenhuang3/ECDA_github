!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_psps
!! NAME
!!  m_psps
!!
!! FUNCTION
!!  This module provides method to allocate/free/initialize the 
!!  pseudopotential_type object.
!!
!! COPYRIGHT
!!  Copyright (C) 2014 ABINIT group (XG,DC,MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_psps
    
 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_io_tools,      only : open_file
 use defs_datatypes,  only : pspheader_type, pseudopotential_type, pseudopotential_gth_type 
 use defs_abitypes,   only : dataset_type

 implicit none

 private

 public :: psps_init_global        ! Allocate and init all part of psps structure that are independent of a given dataset.
 public :: psps_init_from_dtset    ! Allocate and init all part of psps structure that are dependent of a given dataset.
 public :: psps_free               ! Deallocate all memory of psps structure.
 public :: psps_print
 !public :: psps_plot
!!***

contains 

!!****f* m_psps/psps_init_global
!! NAME
!! psps_init_global
!!
!! FUNCTION
!! Allocate and initialise all part of psps structure that are independent of a given dataset.
!!
!! INPUTS
!! npsp=the number of read pseudo files.
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! psps=<type pseudopotential_type>the pseudopotentials description
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine psps_init_global(mtypalch, npsp, psps, pspheads)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psps_init_global'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mtypalch,npsp
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 integer :: ii, mpsang, n1xccc

! *************************************************************************

!Allocation of some arrays independent of the dataset
 ABI_ALLOCATE(psps%filpsp,(npsp))
 ABI_ALLOCATE(psps%pspcod,(npsp))
 ABI_ALLOCATE(psps%pspdat,(npsp))
 ABI_ALLOCATE(psps%pspso,(npsp))
 ABI_ALLOCATE(psps%pspxc,(npsp))
 ABI_ALLOCATE(psps%title,(npsp))
 ABI_ALLOCATE(psps%zionpsp,(npsp))
 ABI_ALLOCATE(psps%znuclpsp,(npsp))
 call psp2params_init(psps%gth_params, npsp)

 psps%filpsp(1:npsp)=pspheads(1:npsp)%filpsp
 psps%pspcod(1:npsp)=pspheads(1:npsp)%pspcod
 psps%pspdat(1:npsp)=pspheads(1:npsp)%pspdat
 psps%pspso(1:npsp)=pspheads(1:npsp)%pspso
 psps%pspxc(1:npsp)=pspheads(1:npsp)%pspxc
 psps%title(1:npsp)=pspheads(1:npsp)%title
 psps%zionpsp(1:npsp)=pspheads(1:npsp)%zionpsp
 psps%znuclpsp(1:npsp)=pspheads(1:npsp)%znuclpsp

!Set values independant from dtset
 psps%npsp   = npsp
!Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
!mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! might not work with HP compiler
!n1xccc=maxval(pspheads(1:npsp)%xccc)
 mpsang=1
 n1xccc=pspheads(1)%xccc
 do ii=1,psps%npsp
   mpsang=max(pspheads(ii)%lmax+1,mpsang)
   n1xccc=max(pspheads(ii)%xccc,n1xccc)
 end do
 psps%mpsang = mpsang
 psps%n1xccc = n1xccc
!Determine here whether the calculation is PAW
 psps%usepaw  =0
 if (pspheads(1)%pspcod==7.or.pspheads(1)%pspcod==17) psps%usepaw=1  ! If paw, all pspcod necessarily are 7 or 17 (see iofn2)
 psps%mtypalch = mtypalch

end subroutine psps_init_global
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psps_init_from_dtset
!! NAME
!! psps_init_from_dtset
!!
!! FUNCTION
!! Allocate and initialise all part of psps structure that are dependent of a given dataset.
!!
!! INPUTS
!! dtset=<type dataset_type>a given dataset
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! psps=<type pseudopotential_type>the pseudopotentials description
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine psps_init_from_dtset(dtset, idtset, psps, pspheads)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psps_init_from_dtset'
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_57_iovars
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idtset
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 type(pspheader_type),intent(in) :: pspheads(psps%npsp)

!Local variables-------------------------------
!scalars
 integer,save :: dimekb_old=-1,lmnmax_old=-1,lnmax_old=-1,mqgridff_old=0
 integer,save :: mqgridvl_old=0,ntypat_old=-1,usepaw_old=-1
 integer :: ipsp,lmnmax,lmnmaxso,lnmax,lnmaxso,newmqgrid,newmqgriddg,nptsgvec
 real(dp) :: gprimd_orig(3,3)

! *************************************************************************

 psps%optnlxccc   = dtset%optnlxccc
!Determine the number of points needed in reciprocal space to represent the
!pseudopotentials (either set by hand from input variable or set automatically by abinit)
 nptsgvec         = 200 !This has to be chosen one and for all or else ??
 newmqgrid        = dtset%mqgrid
 newmqgriddg      = dtset%mqgriddg
 !JB:Which image to use ? I guess 1 always works
 call matr3inv(dtset%rprimd_orig(:,:,1),gprimd_orig)
 if ( dtset%usewvl == 0) then
   call setmqgrid(newmqgrid,newmqgriddg,dtset%ecut*dtset%dilatmx**2,&
&       dtset%pawecutdg*dtset%dilatmx**2,gprimd_orig,nptsgvec,psps%usepaw)
 else
   call setmqgrid(newmqgrid,newmqgriddg,one,one,gprimd_orig,nptsgvec,psps%usepaw)
 end if
 psps%mqgrid_ff   = newmqgrid
 if (psps%usepaw == 1) then
   psps%mqgrid_vl = newmqgriddg
 else
   psps%mqgrid_vl = newmqgrid
 end if

!Determine the maximum number of projectors, for the set of pseudo atom
 call getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,dtset%mixalch_orig,dtset%nimage,psps%npsp,dtset%npspalch,&
& dtset%ntypat,dtset%ntypalch,pspheads)

 psps%npspalch = dtset%npspalch
 psps%ntypat   = dtset%ntypat
 psps%ntypalch = dtset%ntypalch
 psps%ntyppure = dtset%ntyppure

!Set the flag for reciprocal space or real space calculations
 psps%vlspl_recipSpace = (dtset%icoulomb /= 1)
!changed by RShaltaf
 psps%positron = dtset%positron
 psps%useylm   = dtset%useylm

!Added by T. Rangel for WVL+PAW
 psps%usewvl   = dtset%usewvl

 if (idtset > 1) then
   if (allocated(psps%algalch))  then
     ABI_DEALLOCATE(psps%algalch)
   end if
   if (allocated(psps%mixalch))  then
     ABI_DEALLOCATE(psps%mixalch)
   end if
 end if
 ABI_ALLOCATE(psps%algalch,(psps%ntypalch))
 ABI_ALLOCATE(psps%mixalch,(psps%npspalch,psps%ntypalch))
 psps%algalch(1:psps%ntypalch)=dtset%algalch(1:psps%ntypalch)
!This value will be overwritten elsewhere in case there are different images ...
 psps%mixalch(1:psps%npspalch,1:psps%ntypalch)=dtset%mixalch_orig(1:psps%npspalch,1:psps%ntypalch,1)

!Set mpspso and psps%pspso
!Warning : mpspso might be different for each dataset.
 psps%mpspso=1
 do ipsp=1,dtset%npsp
   if(dtset%nspinor==1)then
     psps%pspso(ipsp)=0
   else
     if(dtset%so_psp(ipsp)/=1)then
       psps%pspso(ipsp)=dtset%so_psp(ipsp)
     else
       psps%pspso(ipsp)=pspheads(ipsp)%pspso
     end if
     if(psps%pspso(ipsp)/=0)psps%mpspso=2
   end if
!  Ideally the following line should not exist, but at present, the space has to be booked
   if(pspheads(ipsp)%pspso/=0)psps%mpspso=2
 end do

!Set mpssoang, lmnmax, lnmax
 if(psps%mpspso==1)then
   psps%mpssoang=psps%mpsang
   psps%lmnmax  =lmnmax
   psps%lnmax   =lnmax
 else
   psps%mpssoang=2*psps%mpsang-1
   psps%lmnmax=lmnmaxso
   psps%lnmax=lnmaxso
 end if
!T. Rangel: for wvl + paw do not change psps%lmnmax
 if (psps%useylm==0 .and. psps%usepaw/=1 ) then
   psps%lmnmax=psps%lnmax
 end if

!Set dimekb
 if (psps%usepaw==0) then
   psps%dimekb=psps%lnmax
 else
   psps%dimekb=psps%lmnmax*(psps%lmnmax+1)/2
 end if

!The following arrays are often not deallocated before the end of the dtset loop
!and might keep their content from one dataset to the other,
!if the conditions are fulfilled
 if(dimekb_old/=psps%dimekb .or. ntypat_old/=dtset%ntypat .or. usepaw_old/=psps%usepaw) then
   if(idtset/=1) then
     if (allocated(psps%ekb))  then
       ABI_DEALLOCATE(psps%ekb)
     end if
   end if
   ABI_ALLOCATE(psps%ekb,(psps%dimekb,dtset%ntypat*(1-psps%usepaw)))
   dimekb_old=psps%dimekb
 end if
 if(lmnmax_old/=psps%lmnmax .or. ntypat_old/=dtset%ntypat)then
   if(idtset/=1) then
     if (allocated(psps%indlmn))  then
       ABI_DEALLOCATE(psps%indlmn)
     end if
   end if
   ABI_ALLOCATE(psps%indlmn,(6,psps%lmnmax,dtset%ntypat))
   lmnmax_old=psps%lmnmax
 end if
 if(mqgridff_old/=psps%mqgrid_ff .or. lnmax_old/=psps%lnmax .or. ntypat_old/=dtset%ntypat)then
   if(idtset/=1) then
     if (allocated(psps%ffspl))  then
       ABI_DEALLOCATE(psps%ffspl)
     end if
     if (allocated(psps%qgrid_ff))  then
       ABI_DEALLOCATE(psps%qgrid_ff)
     end if
   end if
   ABI_ALLOCATE(psps%ffspl,(psps%mqgrid_ff,2,psps%lnmax,dtset%ntypat))
   ABI_ALLOCATE(psps%qgrid_ff,(psps%mqgrid_ff))
   mqgridff_old=psps%mqgrid_ff
   lnmax_old=psps%lnmax
 end if
 if(mqgridvl_old/=psps%mqgrid_vl .or. ntypat_old/=dtset%ntypat)then
   if(idtset/=1) then
     if (allocated(psps%qgrid_vl))  then
       ABI_DEALLOCATE(psps%qgrid_vl)
     end if
     if (allocated(psps%vlspl))  then
       ABI_DEALLOCATE(psps%vlspl)
     end if
   end if
   if (idtset/=1 .and. .not.psps%vlspl_recipSpace) then
     if (allocated(psps%dvlspl))  then
       ABI_DEALLOCATE(psps%dvlspl)
     end if
   end if
   ABI_ALLOCATE(psps%vlspl,(psps%mqgrid_vl,2,dtset%ntypat))
   ABI_ALLOCATE(psps%qgrid_vl,(psps%mqgrid_vl))
   if (.not.psps%vlspl_recipSpace) then
     ABI_ALLOCATE(psps%dvlspl,(psps%mqgrid_vl,2,dtset%ntypat))
   end if
   mqgridvl_old=psps%mqgrid_vl
 end if
 if(ntypat_old/=dtset%ntypat.or. usepaw_old/=psps%usepaw)then
   if(idtset/=1) then
     if (allocated(psps%xccc1d))  then
       ABI_DEALLOCATE(psps%xccc1d)
     end if
   end if
   ABI_ALLOCATE(psps%xccc1d,(psps%n1xccc*(1-psps%usepaw),6,dtset%ntypat))
   usepaw_old=psps%usepaw
 end if
 if(ntypat_old/=dtset%ntypat)then
   if(idtset/=1) then
     if (allocated(psps%xcccrc))  then
       ABI_DEALLOCATE(psps%xcccrc)
     end if
     if (allocated(psps%ziontypat))  then
       ABI_DEALLOCATE(psps%ziontypat)
     end if
     if (allocated(psps%znucltypat))  then
       ABI_DEALLOCATE(psps%znucltypat)
     end if
   end if
   ABI_ALLOCATE(psps%xcccrc,(dtset%ntypat))
   ABI_ALLOCATE(psps%znucltypat,(dtset%ntypat))
   ABI_ALLOCATE(psps%ziontypat,(dtset%ntypat))
   ntypat_old=dtset%ntypat
 end if
 psps%ziontypat(:)=dtset%ziontypat(:)

end subroutine psps_init_from_dtset
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psps_free
!! NAME
!! psps_free
!!
!! FUNCTION
!! Deallocate all memory of psps structure.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! psps=<type pseudopotential_type>the pseudopotentials description
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine psps_free(psps)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psps_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(pseudopotential_type),intent(inout) :: psps

! *************************************************************************

!Allocation of some arrays independent of the dataset
 if (allocated(psps%filpsp))  then
   ABI_DEALLOCATE(psps%filpsp)
 end if
 if (allocated(psps%pspcod))  then
   ABI_DEALLOCATE(psps%pspcod)
 end if
 if (allocated(psps%pspdat))  then
   ABI_DEALLOCATE(psps%pspdat)
 end if
 if (allocated(psps%pspso))  then
   ABI_DEALLOCATE(psps%pspso)
 end if
 if (allocated(psps%pspxc))  then
   ABI_DEALLOCATE(psps%pspxc)
 end if
 if (allocated(psps%title))  then
   ABI_DEALLOCATE(psps%title)
 end if
 if (allocated(psps%zionpsp))  then
   ABI_DEALLOCATE(psps%zionpsp)
 end if
 if (allocated(psps%znuclpsp))  then
   ABI_DEALLOCATE(psps%znuclpsp)
 end if

 call psp2params_free(psps%gth_params)

 if (allocated(psps%algalch))  then
   ABI_DEALLOCATE(psps%algalch)
 end if
 if (allocated(psps%mixalch))  then
   ABI_DEALLOCATE(psps%mixalch)
 end if
 if (allocated(psps%ekb))  then
   ABI_DEALLOCATE(psps%ekb)
 end if
 if (allocated(psps%indlmn))  then
   ABI_DEALLOCATE(psps%indlmn)
 end if
 if (allocated(psps%ffspl))  then
   ABI_DEALLOCATE(psps%ffspl)
 end if
 if (allocated(psps%qgrid_ff))  then
   ABI_DEALLOCATE(psps%qgrid_ff)
 end if
 if (allocated(psps%qgrid_vl))  then
   ABI_DEALLOCATE(psps%qgrid_vl)
 end if
 if (allocated(psps%vlspl))  then
   ABI_DEALLOCATE(psps%vlspl)
 end if

 if (.not.psps%vlspl_recipSpace) then
   if (allocated(psps%dvlspl))  then
     ABI_DEALLOCATE(psps%dvlspl)
   end if
 end if
 if (allocated(psps%xccc1d))  then
   ABI_DEALLOCATE(psps%xccc1d)
 end if
 if (allocated(psps%xcccrc))  then
   ABI_DEALLOCATE(psps%xcccrc)
 end if
 if (allocated(psps%ziontypat))  then
   ABI_DEALLOCATE(psps%ziontypat)
 end if
 if (allocated(psps%znucltypat))  then
   ABI_DEALLOCATE(psps%znucltypat)
 end if

end subroutine psps_free
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psps_print
!! NAME
!! psps_print
!!
!! FUNCTION
!!  Method to print the content of a pseudopotential_type derived type
!!
!! INPUTS
!!  psps=<type pseudopotential_type>
!!  unit(optional)=unit number for output
!!  prtvol(optional)=verbosity level
!!  mode_paral(optional): either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only writing 
!!
!! NOTES
!!  Should add information coming from pspheads
!!
!! PARENTS
!!      bethe_salpeter,pspini,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine psps_print(psps,unit,prtvol,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psps_print'
 use interfaces_14_hidewrite
 use interfaces_57_iovars
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(pseudopotential_type),intent(in) :: psps

!Local variables-------------------------------
!scalars
 integer :: ierr,ips,ipsp_alch,ityp_alch,itypat,unt,verb
 character(len=4) :: mode
 character(len=500) :: msg
!arrays
 integer :: cond_values(3)
 character(len=9) :: cond_string(3)

! *************************************************************************

 ! Some initialisations
 verb=0     ; if (PRESENT(prtvol)) verb=prtvol
 unt=std_out; if (PRESENT(unit)) unt=unit
 mode='COLL'; if (PRESENT(mode_paral)) mode=mode_paral
 ierr=0 ; cond_string(1:3)=' ' ; cond_values(1:3)=(/0,0,0/)
 !
 ! /** General info including spin-orbit **/
 write(msg,'(2a)')ch10,' ==== Info on pseudopotentials ==== '
 call wrtout(unt,msg,mode)

 SELECT CASE (psps%usepaw) 
 CASE (0)
   write(msg,'(a)')'  Norm-conserving pseudopotentials '
   call wrtout(unt,msg,mode)
   write(msg,'(a,i4)')'  Max number of Kleinman-Bylander energies ',psps%dimekb
   call wrtout(unt,msg,mode)
   !do itypat=1,psps%ntypat 
   ! write(msg,'(a,i4,a,f9.4)')' Type ',itypat,' K-B energies ',(psps%ekb(ikbe,itypat),ikbe=1,psps%dimekb)
   !end do
 CASE (1)
   write(msg,'(a)')'  PAW calculation'
   call wrtout(unt,msg,mode)
   write(std_out,*)'  Max number of D_ij coefficients ',psps%dimekb
 CASE DEFAULT 
   call chkint_eq(0,0,cond_string,cond_values,ierr,'usepaw',psps%usepaw,2,(/0,1/),unt)
 END SELECT

 !integer :: dimekb
 ! Dimension of Ekb
 ! ->Norm conserving : Max. number of Kleinman-Bylander energies
 !                     for each atom type
 !                     dimekb=lnmax (lnmax: see this file)
 ! ->PAW : Max. number of Dij coefficients connecting projectors
 !                     for each atom type
 !                     dimekb=lmnmax*(lmnmax+1)/2 (lmnmax: see this file)

 !real(dp), pointer :: ekb(:,:)
  ! ekb(dimekb,ntypat*(1-usepaw))
  !  ->NORM-CONSERVING PSPS ONLY:
  !    (Real) Kleinman-Bylander energies (hartree)
  !           for number of basis functions (l,n) (lnmax)
  !           and number of atom types (ntypat)
  ! NOTE (MT) : ekb (norm-conserving) is now diagonal (one dimension
  !             lnmax); it would be easy to give it a second
  !             (symmetric) dimension by putting
  !             dimekb=lnmax*(lnmax+1)/2
  !             in the place of dimekb=lmnmax.
 SELECT CASE (psps%positron)
 CASE (0) 
   !write(std_out,*)' Standard Electron Calculation '
 CASE (1,2)
   write(msg,'(a,i3)')'  Positron Calculation with positron .. ',psps%positron 
   call wrtout(unt,msg,mode)
 CASE DEFAULT
   call chkint_eq(0,0,cond_string,cond_values,ierr,'positron',psps%positron,3,(/0,1,2/),unt)
 END SELECT

 write(msg,'(a,i4,2a,i4)')&
&  '  Number of pseudopotentials .. ',psps%npsp,ch10,&
&  '  Number of types of atoms   .. ',psps%ntypat 
 call wrtout(unt,msg,mode)

 SELECT CASE (psps%mpspso) 
 CASE (1) 
   write(msg,'(a)')'  Calculation without spin-orbit '
   call wrtout(unt,msg,mode)
 CASE (2)
   write(msg,'(3a,i3)')&
&    '  Calculation with spin-orbit coupling ',ch10,&
&    '  Max number of channels (spin-orbit included) ',psps%mpssoang
   call wrtout(unt,msg,mode)
   do itypat=1,psps%ntypat 
     if (psps%pspso(itypat)==2) then 
       write(msg,'(a,i4,a)')'  - Atom type ',itypat,' has spin-orbit characteristics'
       call wrtout(unt,msg,mode)
     end if 
   end do
 CASE DEFAULT
   call chkint_eq(0,0,cond_string,cond_values,ierr,'mpspso',psps%mpspso,2,(/1,2/),unt)
 END SELECT
 !
 ! /** Info on nonlocal part **/
 !
 SELECT CASE (psps%useylm)
 CASE (0)
   write(msg,'(a)')'  Nonlocal part applied using Legendre polynomials '
 CASE (1)
   write(msg,'(a)')'  Nonlocal part applied using real spherical harmonics '
 CASE DEFAULT
   call chkint_eq(0,0,cond_string,cond_values,ierr,'psps%useylm',psps%useylm,2,(/0,1/),unt)
 END SELECT
 call wrtout(unt,msg,mode)

 !FIXME this does not work, it seems it is always 0 , except for HGH
 !write(msg,'(a,i3)')' Max number of non-local projectors over l and type ',psps%mproj 
 !if (psps%mproj==0) then 
 ! write(msg,'(a)')TRIM(msg)//' (All local) '
 !end if
 !call wrtout(unt,msg,mode)
 write(msg,'(a,i3,2a,i3,2a,i3)')&
& '  Highest angular momentum +1 ....... ',psps%mpsang,ch10,&
& '  Max number of (l,n)   components .. ',psps%lnmax, ch10,&
& '  Max number of (l,m,n) components .. ',psps%lmnmax
 call wrtout(unt,msg,mode)
  !integer :: lnmax
  !  Max. number of (l,n) components over all type of psps
  !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
  !  so, it is equal to the max of lnprojso, see pspheader_type
 !integer :: lmnmax
  !  If useylm=0, max number of (l,m,n) comp. over all type of psps (lnproj)
  !  If useylm=1, max number of (l,n)   comp. over all type of psps (lmnproj)
  !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
  !  so, it is equal to the max of lmnprojso or lnprojso, see pspheader_type

!$integer, pointer :: indlmn(:,:,:)
! indlmn(6,lmnmax,ntypat)
! For each type of psp,
! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
!                                or i=lmn (if useylm=1)

 !FIXME for paw n1xccc==1
 !
 ! /** Non-linear Core correction **/
 if (psps%n1xccc/=0) then 
  write(msg,'(3a,2(a,i4,a),2a)')ch10,&
&   ' *** Pseudo-Core Charge Info *** ',ch10,&
&   '  Number of radial points for pseudo-core charge .. ',psps%n1xccc,ch10,&
&   '  XC core-correction treatment (optnlxccc) ........ ',psps%optnlxccc,ch10,&
&   '  Radius for pseudo-core charge for each type ..... ',ch10
  call wrtout(unt,msg,mode)
  do itypat=1,psps%ntypat 
    write(msg,'(a,i4,a,f7.4)')'  - Atom type ',itypat,' has pseudo-core radius .. ',psps%xcccrc(itypat)
    call wrtout(unt,msg,mode)
  end do
 end if
 !
 ! /** Alchemical mixing **/
 if (psps%mtypalch/=0) then 
  write(msg,'(3a,3(a,i4,a))')ch10,&
&   ' *** Calculation with alchemical mixing *** ',ch10,&
&   '  Number of pure pseudoatoms .... ',psps%ntyppure,ch10,&
&   '  Number of pseudos for mixing .. ',psps%npspalch,ch10,&
&   '  Alchemical pseudoatoms ........ ',psps%ntypalch,ch10
  call wrtout(unt,msg,mode)
  do ipsp_alch=1,psps%npspalch 
    do ityp_alch=1,psps%ntypalch 
      write(std_out,*)' mixalch ',psps%mixalch(ipsp_alch,ityp_alch)
    end do
  end do
  do ityp_alch=1,psps%ntypalch 
    write(msg,'(a,i4,a,i4)')' For alchemical atom no. ',ityp_alch,' algalch is .. ',psps%algalch(ityp_alch)
    call wrtout(unt,msg,mode)
  end do
 end if
 !integer :: mtypalch
  ! Maximum number of alchemical pseudo atoms. If non-zero,
  ! the mechanism to generate mixing of pseudopotentials is activated
 !integer :: ntypat
  ! Number of types of atoms (might be alchemy wrt pseudopotentials)
 !integer :: ntyppure
  ! Number of types of pure pseudoatoms
 !integer :: ntypalch
  ! Number of types of alchemical pseudoatoms
 !integer :: npspalch
  ! Number of types of pseudopotentials use for alchemical purposes
 !integer, pointer :: algalch(:)   ! algalch(ntypalch)
  ! For each type of pseudo atom, the algorithm to mix the pseudopotentials
 !real(dp), pointer :: mixalch(:,:)
  ! mixalch(npspalch,ntypalch)
  ! Mixing coefficients to generate alchemical pseudo atoms

 !
 ! /** Info in Q-grid for spline of form factors **/
 !
 write(msg,'(3a,a,i6,a,a,i6)')ch10,&
&  ' *** Info on the Q-grid used for form factors in spline form *** ',ch10,&
&  '  Number of q-points for radial functions ffspl .. ',psps%mqgrid_ff,ch10,&
&  '  Number of q-points for vlspl ................... ',psps%mqgrid_vl 
 call wrtout(unt,msg,mode)
 if (psps%vlspl_recipSpace) then 
   call wrtout(unt,'  vlspl is computed in Reciprocal Space ',mode)
 else 
   call wrtout(unt,'  vlsp is computed in Real Space ',mode)
 end if
 !TODO additional stuff tbat might be printed

 !real(dp), pointer :: ffspl(:,:,:,:)
  ! ffspl(mqgrid_ff,2,lnmax,ntypat)
  ! Gives, on the radial grid, the different non-local projectors,
  ! in both the norm-conserving case, and the PAW case
 !real(dp), pointer :: qgrid_ff(:)
  ! qgrid_ff(mqgrid_ff)
  ! The coordinates of all the points of the radial grid for the nl form factors
 !real(dp), pointer :: qgrid_vl(:)
   ! qgrid_vl(mqgrid_vl)
   ! The coordinates of all the points of the radial grid for the local part of psp
 !real(dp), pointer :: vlspl(:,:,:)
  ! vlspl(mqgrid_vl,2,ntypat)
  ! Gives, on the radial grid, the local part of each type of psp.
 ! real(dp), pointer :: dvlspl(:,:,:)
  ! dvlspl(mqgrid_vl,2,ntypat)
  ! Gives, on the radial grid, the first derivative of the local
  ! part of each type of psp (computed when the flag 'vlspl_recipSpace'
  ! is true).
 !real(dp), pointer :: xccc1d(:,:,:)
  ! xccc1d(n1xccc*(1-usepaw),6,ntypat)
  ! Norm-conserving psps only
  ! The component xccc1d(n1xccc,1,ntypat) is the pseudo-core charge
  ! for each type of atom, on the radial grid. The components
  ! xccc1d(n1xccc,ideriv,ntypat) give the ideriv-th derivative of the
  ! pseudo-core charge with respect to the radial distance.

  !write(msg,'(2a)')ch10,' Z_ion pseudo Z_at ' 
  !call wrtout(unt,msg,mode)
  !do ips=1,psps%npsp
  ! write(std_out,*)psps%zionpsp(ips),psps%znuclpsp(ips)
  !end do

   !real(dp), pointer :: zionpsp(:)
   ! zionpsp(npsp)
   ! For each pseudopotential, the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)
   !real(dp), pointer :: ziontypat(:)
   ! ziontypat(ntypat)
   !  For each type of atom (might be alchemy wrt psps), the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)
   !real(dp), pointer :: znuclpsp(:)
   ! znuclpsp(npsp)
   ! The atomic number of each pseudopotential
   !real(dp), pointer :: znucltypat(:)
   ! znucltypat(ntypat)
   ! The atomic number of each type of atom (might be alchemy wrt psps)

  do itypat=1,psps%ntypat 
    write(msg,'(a,i3,a,i3)')' XC functional for type ',itypat,' is ',psps%pspxc(itypat)
    call wrtout(unt,msg,mode)
    !write(std_out,*)psps%ziontypat(itypat),psps%znucltypat(itypat)
  end do
  !integer, pointer :: pspxc(:)
   ! pspxc(ntypat)
   ! For each type of psp, the XC functional that was used to generate it,
   ! as given by the psp file

  if (verb>=3) then 
    do ips=1,psps%npsp
      write(std_out,*)' Pseudo number   ',ips,' read from ',TRIM(psps%filpsp(ips))
      write(std_out,*)' Format or Code  ',psps%pspcod(ips)
      write(std_out,*)' Generation Date ',psps%pspdat(ips)
      write(std_out,*)' Content of first line ',TRIM(psps%title(ips))
    end do
  end if

  !character(len=fnlen), pointer :: filpsp(:)
   ! filpsp(ntypat)
   ! The filename of the pseudopotential
  !character(len=fnlen), pointer :: title(:)
   ! title(ntypat)
   ! The content of first line read from the psp file
!  integer, pointer :: pspdat(:)
   ! pspdat(ntypat)
   ! For each type of psp, the date of psp generation, as given by the psp file
  !integer, pointer :: pspcod(:)
   ! pspcod(npsp)
   ! For each type of psp, the format -or code- of psp generation,
   !  as given by the psp file


! Types for pseudo-potentials that are based on parameters. Currently, only
! GTH are supported (see pseudopotential_gth_type). To add one, one should
! create an initialisation method and a destruction method in 02psp (see psp2params.F90). 
! These methods are called in driver().
!TODO this is still missing
!  type(pseudopotential_gth_type) :: gth_params

 ! If there was a problem, then stop.
 if (ierr/=0) then 
    MSG_ERROR("Fatal error")
 end if

end subroutine psps_print
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psps_plot
!! NAME
!! psps_plot
!!
!! FUNCTION
!!  Writes on external files some of the arrays defined in the pseudopotential_type.
!!
!! INPUTS
!!   root_filename=File prefix.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine psps_plot(psps,root_filename)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psps_plot'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in),optional :: root_filename
 type(pseudopotential_type),intent(in) :: psps

!Local variables-------------------------------
!scalars
 integer :: ider,iln,iq,ir,ityp,unt
 character(len=100) :: fmt
 character(len=500) :: msg
 character(len=fnlen) :: fname,root

! *************************************************************************

 root='PSPS'; if (present(root_filename)) root=root_filename

 !TODO most of the pointer are not nullified, 
 !this part will print a lot of quantities that are not used

 if (allocated(psps%vlspl)) then 
   fname=trim(root)//'_VLSPL'
   if (open_file(fname,msg,newunit=unt,status='new',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write(unt,*)' Local part of each type of atom '  
   write(unt,*)' q-mesh and, for each type, v_loc(q)? and second derivative '
   write(fmt,*)'(',1+2*psps%ntypat,'(es17.9,1x))'

   do iq=1,psps%mqgrid_vl
     write(unt,fmt)psps%qgrid_vl(iq),((psps%vlspl(iq,ider,ityp),ider=1,2),ityp=1,psps%ntypat)
   end do

   close(unt)
 end if

 if (allocated(psps%ffspl)) then 
   !TODO write error handler for open
   !here I need the pseudo_header to avoid writing columns made of zero 
   fname=trim(root)//'_FFSPL'
   if (open_file(fname,msg,newunit=unt,status='new',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write(unt,*)' Form factors for each type of atom '  
   write(unt,*)' q-mesh and, for each type and each (l,n) channel, ffnl(q) and second derivative '
   write(fmt,*)'(',1+2*psps%lnmax*psps%ntypat,'(es17.9,1x))'

   do iq=1,psps%mqgrid_ff
     write(unt,fmt)psps%qgrid_ff(iq),&
&      (((psps%ffspl(iq,ider,iln,ityp),ider=1,2),iln=1,psps%lnmax),ityp=1,psps%ntypat)
   end do

   close(unt)
 end if

 if (allocated(psps%xccc1d)) then 
   !TODO write error handler for open

   fname=trim(root)//'_PSCC'
   if (open_file(fname,msg,newunit=unt,status='new',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write(unt,*)' Pseudo-core charge for each type of atom, on the radial grid. '
   write(unt,*)' radial-mesh and, for each type rho_pscore(r) '
   write(fmt,*)'(',psps%ntypat,'(es17.9,1x))'

   ! TODO Grid is missing
   do ir=1,psps%n1xccc
     !write(unt,fmt)((psps%xccc1d(ir,1,itypat)),itypat=1,psps%ntypat)
   end do

   close(unt)
 end if

end subroutine psps_plot 
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psp2params_init
!! NAME
!! psp2params_init
!!
!! FUNCTION
!! Allocate and initialise the data structure holding parameters for the GTH
!! pseudo-potentials.
!!
!!  MJV note: this should be renamed: psp2 suggests it relates to pspcod 2,
!!     whereas it is actually 3 
!!    the parameters would also be better off separated into C and h arrays
!!
!! INPUTS
!!  npsp=number of true pseudo used (not alchemy).
!!
!! OUTPUT
!!  gth_params <type (pseudopotential_gth_type)>=the values to allocate and initialise.
!!
!! PARENTS
!!      m_psps
!!
!! CHILDREN
!!
!! SOURCE

subroutine psp2params_init(gth_params, npsp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp2params_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npsp
 type(pseudopotential_gth_type),intent(out) :: gth_params

! *********************************************************************

!Check array, no params are currently set.
 ABI_ALLOCATE(gth_params%set,(npsp))
 gth_params%set(:) = .false.

!Check array, have geometric informations been filled?
 ABI_ALLOCATE(gth_params%hasGeometry,(npsp))
 gth_params%hasGeometry(:) = .false.

!Coefficients for local part and projectors
 ABI_ALLOCATE(gth_params%psppar,(0:4, 0:6, npsp))
 gth_params%psppar = zero

!Coefficients for spin orbit part
 ABI_ALLOCATE(gth_params%psp_k_par,(1:4, 1:3, npsp))
 gth_params%psp_k_par = zero

!Different radii
 ABI_ALLOCATE(gth_params%radii_cov,(npsp))
 ABI_ALLOCATE(gth_params%radii_cf,(npsp, 3))

!Number of semicore electrons
 ABI_ALLOCATE(gth_params%semicore,(npsp))

end subroutine psp2params_init
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psp2params_free
!! NAME
!! psp2params_free
!!
!! FUNCTION
!! Deallocate a previously allocated data structure for storage of GTH parameters.
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  gth_params <type (pseudopotential_gth_type)>=the values to deallocate.
!!
!! PARENTS
!!      m_psps
!!
!! CHILDREN
!!
!! SOURCE

subroutine psp2params_free(gth_params)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp2params_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(pseudopotential_gth_type),intent(inout) :: gth_params

! *********************************************************************

!Check arrays.
 if (allocated(gth_params%set))  then
   ABI_DEALLOCATE(gth_params%set)
 end if
 if (allocated(gth_params%hasGeometry))  then
   ABI_DEALLOCATE(gth_params%hasGeometry)
 end if

!Coefficients for local part and projectors
 if (allocated(gth_params%psppar))  then
   ABI_DEALLOCATE(gth_params%psppar)
 end if

!Coefficients for spin orbit part
 if (allocated(gth_params%psp_k_par))  then
   ABI_DEALLOCATE(gth_params%psp_k_par)
 end if

!Different radii
 if (allocated(gth_params%radii_cov))  then
   ABI_DEALLOCATE(gth_params%radii_cov)
 end if
 if (allocated(gth_params%radii_cf))  then
   ABI_DEALLOCATE(gth_params%radii_cf)
 end if

!Number of semicore electrons
 if (allocated(gth_params%semicore))  then
   ABI_DEALLOCATE(gth_params%semicore)
 end if

end subroutine psp2params_free
!!***

end module m_psps
!!***
