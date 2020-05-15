!{\src2tex{textfont=tt}}
!!****m* ABINIT/libxc_functionals
!! NAME
!!  libxc_functionals
!!
!! FUNCTION
!!  Module containing interfaces to the LibXC library, for exchange
!!  correlation potentials and energies. The interfacing between
!!  the ABINIT and LibXC formats and datastructures happens here.
!!  Also contains basic container datatype for LibXC interfacing.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (MOliveira,LHH,FL,GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  libxc_functionals.F90 uses a global structure (funcs) to store the XC parameters.
!!  This structure is initialized in driver with the value of ixc specified by the user in the input file. 
!!  In order to change the value of ixc at run-time, we have to reinitialize the global structure 
!!  with the new value of ixc before computing XC quantities.
!!  Moreover one has to reinstate the old functional before returning so that the other routines
!!  will continue to used the previous ixc. This task can be accomplished with the following pseudocode
!!  
!!  #ifdef HAVE_DFT_LIBXC
!!   ! Reinitialize the libxc module with the overriden values
!!   if (old_ixc<0) call libxc_functionals_end()
!!   if (new_ixc<0) call libxc_functionals_init(new_ixc,nspden)
!!  #endif
!!
!!  ! Compute XC stuff here.
!!  
!!  #ifdef HAVE_DFT_LIBXC
!!   ! Revert libxc module to the original settings
!!   if (new_ixc<0) call libxc_functionals_end()
!!   if (old_ixc<0) call libxc_functionals_init(old_ixc,nspden)
!!  #endif
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

module libxc_functionals

 use defs_basis
 use m_profiling_abi
 use m_errors

#if defined HAVE_DFT_LIBXC
 use xc_f90_types_m
 use libxc_funcs_m
 use xc_f90_lib_m
#endif

 implicit none

#if defined HAVE_DFT_LIBXC
 type libxc_functional
   private
   integer                :: family  ! LDA, GGA, etc.
   integer                :: id      ! identifier
   integer                :: nspin   ! # of spin components
   logical                :: has_fxc ! TRUE is fxc is available for the functional
   type(xc_f90_pointer_t) :: conf    ! the pointer used to call the library
   type(xc_f90_pointer_t) :: info    ! information about the functional
 end type libxc_functional

 type(libxc_functional) :: funcs(2)

 private

 integer,save :: ABI_IXC = HUGE(0)

 public :: libxc_functionals_init        ! Initialize the desired XC functional, from LibXC.
 public :: libxc_functionals_fullname    ! Return full name of the XC functional
 public :: libxc_functionals_getid       ! Return identifer of a XC functional from its name
 public :: libxc_functionals_global_ixc  ! The value of ixc used to initialize the global structure funcs
 public :: libxc_functionals_getvxc      ! Return XC potential and energy, from input density (event gradient etc...)
 public :: libxc_functionals_isgga
 public :: libxc_functionals_ismgga
 public :: libxc_functionals_has_kxc
 public :: libxc_functionals_nspin       ! The number of spin components for the XC functionals
 public :: libxc_functionals_end         ! End usage of LibXC functional.

contains
!!***

!!****f* libxc_functionals/libxc_functionals_init
!! NAME
!!  libxc_functionals_init
!!
!! FUNCTION
!!  Initialize the desired XC functional, from LibXC.
!!  * Call the LibXC initializer
!!  * Fill preliminary fields in module structures.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_vhxc_me,driver,m_kxc
!!
!! CHILDREN
!!      wrtout,xc_f90_gga_exc_vxc,xc_f90_gga_fxc,xc_f90_gga_vxc
!!      xc_f90_lda_exc_vxc,xc_f90_lda_fxc,xc_f90_lda_kxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc,xc_f90_mgga_x_tb09_set_par
!!
!! SOURCE

 subroutine libxc_functionals_init(ixc,nspden)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_init'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nspden
 integer, intent(in) :: ixc

!Local variables-------------------------------
!scalars
 integer :: ii,jj,nspden_eff
 character(len=500) :: message
 type(xc_f90_pointer_t) :: str

! *************************************************************************
 ! Save abinit value for reference 
 ABI_IXC = ixc

 funcs(1)%id = -ixc/1000
 funcs(2)%id = -ixc - funcs(1)%id*1000

 nspden_eff=min(nspden,2)
 funcs(1)%nspin = nspden_eff
 funcs(2)%nspin = nspden_eff

 do ii = 1, 2
   if (funcs(ii)%id == 0) then
     funcs(ii)%family = 0
     funcs(ii)%has_fxc= .false.
     cycle
   end if

   ! Get XC functional family
   funcs(ii)%family = xc_f90_family_from_id(funcs(ii)%id)
   select case (funcs(ii)%family)
   case (XC_FAMILY_LDA, XC_FAMILY_GGA,XC_FAMILY_HYB_GGA,XC_FAMILY_MGGA)
     call xc_f90_func_init(funcs(ii)%conf,funcs(ii)%info,funcs(ii)%id,nspden_eff)
     funcs(ii)%has_fxc=(iand(xc_f90_info_flags(funcs(ii)%info),XC_FLAGS_HAVE_FXC)>0)
   case default
     write(message, '(a,i8,2a,i8,6a)' )&
&      'Invalid IXC = ',ixc,ch10,&
&      'The LibXC functional family ',funcs(ii)%family,&
&      'is currently unsupported by ABINIT',ch10,&
&      '(-1 means the family is unknown to the LibXC itself)',ch10,&
&      'Please consult the LibXC documentation',ch10
     MSG_ERROR(message)
   end select

   if (funcs(ii)%id == XC_LDA_C_XALPHA) then
     call xc_f90_lda_c_xalpha_set_par(funcs(ii)%conf,zero)
   end if

   ! Dump functional information
   call xc_f90_info_name(funcs(ii)%info,message)
   call wrtout(std_out,message,'COLL')
   jj = 0
   call xc_f90_info_refs(funcs(ii)%info,jj,str,message)
   do while (jj >= 0)
     call wrtout(std_out,message,'COLL')
     call xc_f90_info_refs(funcs(ii)%info,jj,str,message)
   end do
 end do

end subroutine libxc_functionals_init
!!***

!----------------------------------------------------------------------

!!****f* libxc_functionals/libxc_functionals_end
!! NAME
!!  libxc_functionals_end
!!
!! FUNCTION
!!  End usage of LibXC functional. Call LibXC end function,
!!  and deallocate module contents.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_vhxc_me,driver,m_kxc
!!
!! CHILDREN
!!      wrtout,xc_f90_gga_exc_vxc,xc_f90_gga_fxc,xc_f90_gga_vxc
!!      xc_f90_lda_exc_vxc,xc_f90_lda_fxc,xc_f90_lda_kxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc,xc_f90_mgga_x_tb09_set_par
!!
!! SOURCE
 subroutine libxc_functionals_end()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_end'
!End of the abilint section

 implicit none

!Local variables-------------------------------
 integer :: ii

! *************************************************************************
 ABI_IXC = HUGE(0)

 do ii = 1, 2
   if (funcs(ii)%id == 0) cycle
   call xc_f90_func_end(funcs(ii)%conf)
 end do

 end subroutine libxc_functionals_end
!!***

!----------------------------------------------------------------------

!!****f* libxc_functionals/libxc_functionals_fullname
!! NAME
!!  libxc_functionals_fullname
!!
!! FUNCTION
!!  Return full name of the XC functional
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
 function libxc_functionals_fullname()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_fullname'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=100) :: libxc_functionals_fullname

!Local variables-------------------------------
 character(len=100) :: xcname

! *************************************************************************

 if (funcs(1)%id == 0) then
   if (funcs(2)%id /= 0) then
     call xc_f90_info_name(funcs(2)%info,libxc_functionals_fullname)
   else
     libxc_functionals_fullname='No XC functional'
   end if
 else if (funcs(2)%id == 0) then
   if (funcs(1)%id /= 0) then
     call xc_f90_info_name(funcs(1)%info,libxc_functionals_fullname)
   else
     libxc_functionals_fullname='No XC functional'
   end if
 else
   call xc_f90_info_name(funcs(1)%info,libxc_functionals_fullname)
   call xc_f90_info_name(funcs(2)%info,xcname)
   libxc_functionals_fullname=trim(libxc_functionals_fullname)//'+'//trim(xcname)
 end if

end function libxc_functionals_fullname
!!***

!----------------------------------------------------------------------

!!****f* libxc_functionals/libxc_functionals_getid
!! NAME
!!  libxc_functionals_getid
!!
!! FUNCTION
!!  Return identifer of a XC functional from its name
!!  Return -1 if undefined 
!!
!! INPUTS
!!  xcname= string containing the name of a XC functional
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
 function libxc_functionals_getid(xcname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_getid'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: libxc_functionals_getid
 character(len=*),intent(in) :: xcname

 !Local variables-------------------------------
 integer,external :: xc_f90_functional_get_number
 character(len=256) :: str

! *************************************************************************

 str=xcname
 if (xcname(1:3)=="XC_".or.xcname(1:3)=="xc_") str=xcname(4:)

 libxc_functionals_getid=xc_f90_functional_get_number(trim(str))

end function libxc_functionals_getid
!!***

!----------------------------------------------------------------------

!!****f* libxc_functionals/libxc_functionals_global_ixc
!! NAME
!!  libxc_functionals_global_ixc
!!
!! FUNCTION
!!  Return the value of ixc used to initialize the global structure funcs
!!  Return HUGE(0) if undefined 
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function libxc_functionals_global_ixc()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_global_ixc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: libxc_functionals_global_ixc

! *************************************************************************

 libxc_functionals_global_ixc = ABI_IXC

end function libxc_functionals_global_ixc
!!***

!----------------------------------------------------------------------

!!****f* libxc_functionals/libxc_functionals_isgga
!! NAME
!!  libxc_functionals_isgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a GGA or not
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function libxc_functionals_isgga()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_isgga'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 logical :: libxc_functionals_isgga

! *************************************************************************

 if (any(funcs%family == XC_FAMILY_GGA) .or. any(funcs%family == XC_FAMILY_HYB_GGA)) then
   libxc_functionals_isgga = .true.
 else
   libxc_functionals_isgga = .false.
 end if

end function libxc_functionals_isgga
!!***

!----------------------------------------------------------------------

!!****f* libxc_functionals/libxc_functionals_ismgga
!! NAME
!!  libxc_functionals_ismgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a Meta-GGA or not
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function libxc_functionals_ismgga()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_ismgga'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 logical :: libxc_functionals_ismgga

! *************************************************************************

 if (any(funcs%family == XC_FAMILY_MGGA)) then
   libxc_functionals_ismgga = .true.
 else
   libxc_functionals_ismgga = .false.
 end if

end function libxc_functionals_ismgga
!!***

!----------------------------------------------------------------------

!!****f* libxc_functionals/libxc_functionals_has_kxc
!! NAME
!!  libxc_functionals_has_kxc
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  provides Kxc or not fxc in the libXC convention)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function libxc_functionals_has_kxc()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_has_kxc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 logical :: libxc_functionals_has_kxc

!Local variables-------------------------------
 integer :: ii

! *************************************************************************

 libxc_functionals_has_kxc = .true.

 do ii=1,2
   if (funcs(ii)%id/=0) then
     if (.not.funcs(ii)%has_fxc) libxc_functionals_has_kxc = .false.
   end if
 end do

end function libxc_functionals_has_kxc
!!***

!----------------------------------------------------------------------

!!****f* libxc_functionals/libxc_functionals_nspin
!! NAME
!!  libxc_functionals_nspin
!!
!! FUNCTION
!!  Returns the number of spin components for the XC functionals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function libxc_functionals_nspin()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_nspin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: libxc_functionals_nspin

! *************************************************************************

 if (any(funcs%nspin == XC_POLARIZED)) then
   libxc_functionals_nspin = 2
 else
   libxc_functionals_nspin = 1
 end if

end function libxc_functionals_nspin
!!***

!----------------------------------------------------------------------

!!****f* libxc_functionals/libxc_functionals_getvxc
!! NAME
!!  libxc_functionals_getvxc
!!
!! FUNCTION
!!  Return XC potential and energy, from input density (event gradient etc...)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!      wrtout,xc_f90_gga_exc_vxc,xc_f90_gga_fxc,xc_f90_gga_vxc
!!      xc_f90_lda_exc_vxc,xc_f90_lda_fxc,xc_f90_lda_kxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc,xc_f90_mgga_x_tb09_set_par
!!
!! SOURCE

 subroutine libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxc,&
&                 grho2,vxcgr,lrho,vxclrho,tau,vxctau,dvxc,d2vxc,xc_tb09_c) ! Optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_getvxc'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: ndvxc,nd2vxc,npts,nspden,order
 real(dp),intent(in)  :: rho(npts,nspden)
 real(dp),intent(out) :: vxc(npts,nspden), exc(npts)
 real(dp),intent(in),optional :: grho2(npts,2*min(nspden,2)-1)
 real(dp),intent(out),optional :: vxcgr(npts,3)
 real(dp),intent(in),optional :: lrho(npts,nspden)
 real(dp),intent(out),optional :: vxclrho(npts,nspden)
 real(dp),intent(in),optional :: tau(npts,nspden)
 real(dp),intent(out),optional :: vxctau(npts,nspden)
 real(dp),intent(out),optional :: dvxc(npts,ndvxc)
 real(dp),intent(out),optional :: d2vxc(npts,nd2vxc)
 real(dp),intent(in),optional :: xc_tb09_c

!Local variables -------------------------------
 integer  :: i, ipts
 real(dp) :: c, rhotmp(nspden), exctmp, sigma(3), vsigma(3), vxctmp(nspden)
 real(dp) :: v2rho2(3),v2rhosigma(6),v2sigma2(6),v3rho3(4)
 real(dp) :: lrhotmp(nspden), tautmp(nspden), vxclrhotmp(nspden), vxctautmp(nspden)
 real(dp), allocatable :: gnon(:)
 character(len=500) :: message

! *************************************************************************

 ! Inititalize all relevant arrays to zero
 vxc=zero
 exc=zero
 vxctmp=zero
 exctmp=zero


!LHH,FL,GMR
 v2rho2=zero
 v2rhosigma=zero
 v2sigma2=zero
 v3rho3=zero
 if (order**2 >1) dvxc=zero
 if (order**2 >4) d2vxc=zero
!LHH,FL,GMR

 if (any(funcs%family == XC_FAMILY_GGA) .or. any(funcs%family == XC_FAMILY_HYB_GGA)) vxcgr=zero
 if (any(funcs%family == XC_FAMILY_MGGA)) then
   vxcgr=zero
   vxclrho=zero
   vxctau=zero
 end if

 !The TB09 MGGA functional requires an extra quantity
 if (any(funcs%id == XC_MGGA_X_TB09)) then
   if (PRESENT(xc_tb09_c)) then
     c = xc_tb09_c
   else
     ABI_ALLOCATE(gnon,(npts))
     do ipts = 1, npts
       if (sum(rho(ipts, :)) <= 1e-7_dp) then
         gnon(ipts) = zero
       else
         if (nspden == 1) then
           gnon(ipts) = sqrt(grho2(ipts,1))/rho(ipts, 1)
         else
           gnon(ipts) = sqrt(grho2(ipts,3))/sum(rho(ipts, :))
         end if
       end if
     end do
     c = -0.012_dp + 1.023_dp*sqrt(sum(gnon)/npts)
     ABI_DEALLOCATE(gnon)
   end if
   do i = 1, 2
     if (funcs(i)%id == XC_MGGA_X_TB09) then
       call xc_f90_mgga_x_tb09_set_par(funcs(i)%conf, c)
       if (PRESENT(xc_tb09_c)) then
         write(message, '(2a,f9.6)' ) ch10,&
              &     ' In the functional TB09 c is fixed by the user (xc_tb09_c set in input file) and is equal to ', c
       else
         write(message, '(2a,f9.6)' ) ch10,&
              &     ' In the functional TB09 c = ', c
       end if
       call wrtout(std_out,message,'COLL')
     end if
   end do
 end if

 !Loop over points
 do ipts = 1, npts

   ! Convert the quantities provided by ABINIT to the ones needed by libxc
   if (nspden == 1) then
     ! ABINIT passes rho_up in the spin-unpolarized case, while the libxc
     ! expects the total density
     rhotmp(1:nspden) = two*rho(ipts,1:nspden)
   else
     rhotmp(1:nspden) = rho(ipts,1:nspden)
   end if
   if (any(funcs%family == XC_FAMILY_GGA) .or. any(funcs%family == XC_FAMILY_HYB_GGA)&
&      .or. any(funcs%family == XC_FAMILY_MGGA)) then
     sigma=zero
     if (nspden==1) then
       ! ABINIT passes |grho_up|^2 while Libxc needs |grho_tot|^2
       sigma(1) = four*grho2(ipts,1)
     else
       ! ABINIT passes |grho_up|^2, |grho_dn|^2, and |grho_tot|^2
       ! while Libxc needs |grho_up|^2, grho_up.grho_dn, and |grho_dn|^2
       sigma(1) = grho2(ipts,1)
       sigma(2) = (grho2(ipts,3) - grho2(ipts,1) - grho2(ipts,2))/two
       sigma(3) = grho2(ipts,2)
     end if
   end if
   if (any(funcs%family == XC_FAMILY_MGGA)) then
     if (nspden==1) then
       lrhotmp(1:nspden) = two*lrho(ipts,1:nspden)
       tautmp(1:nspden) = two*tau(ipts,1:nspden)
     else
       lrhotmp(1:nspden) = lrho(ipts,1:nspden)
       tautmp(1:nspden) = tau(ipts,1:nspden)
     end if
   end if

   !Loop over functionals
   do i = 1,2
     if (funcs(i)%id == 0) cycle

     !Get the potential (and possibly the energy)
     if (iand(xc_f90_info_flags(funcs(i)%info), XC_FLAGS_HAVE_EXC) .ne. 0) then
       select case (funcs(i)%family)
       case (XC_FAMILY_LDA)
         call xc_f90_lda_exc_vxc(funcs(i)%conf,1,rhotmp(1),exctmp,vxctmp(1))
         if (order**2 > 1) then
           call xc_f90_lda_fxc(funcs(i)%conf,1,rhotmp(1),v2rho2(1))
         endif
         if (order**2 > 4) then
           call xc_f90_lda_kxc(funcs(i)%conf,1,rhotmp(1),v3rho3(1))
         endif
       case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
         call xc_f90_gga_exc_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),exctmp,vxctmp(1),vsigma(1))
         if (order**2 > 1) then
           call xc_f90_gga_fxc(funcs(i)%conf,1,rhotmp(1),sigma(1),v2rho2(1),v2rhosigma(1),v2sigma2(1))
         endif

       case (XC_FAMILY_MGGA)
         call xc_f90_mgga_exc_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),lrhotmp(1),&
                     tautmp(1),exctmp,vxctmp(1),vsigma(1),vxclrhotmp(1),vxctautmp(1))
       end select

     else
       exctmp=zero
       select case (funcs(i)%family)
       case (XC_FAMILY_LDA)
         call xc_f90_lda_vxc(funcs(i)%conf,1,rhotmp(1),vxctmp(1))
       case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
         call xc_f90_gga_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),vxctmp(1),vsigma(1))


       case (XC_FAMILY_MGGA)
         call xc_f90_mgga_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),lrhotmp(1),&
                     tautmp(1),vxctmp(1),vsigma(1),vxclrhotmp(1),vxctautmp(1))
       end select
     end if

     exc(ipts) = exc(ipts) + exctmp
     vxc(ipts,1:nspden) = vxc(ipts,1:nspden) + vxctmp(1:nspden)

!LHH,FL,GMR: deal with fxc and kxc
     if (order**2>1) then
       select case (funcs(i)%family)
       case (XC_FAMILY_LDA)
         if (nspden==1) then
           if(order>=2) then
             dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
             if(order==3) then
               d2vxc(ipts,1)=d2vxc(ipts,1)+v3rho3(1)
             endif
           else
             dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
             dvxc(ipts,2)=dvxc(ipts,2)+v2rho2(1)
           endif
         else
           dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
           dvxc(ipts,2)=dvxc(ipts,2)+v2rho2(2)
           dvxc(ipts,3)=dvxc(ipts,3)+v2rho2(3)
           if(order==3) then
             d2vxc(ipts,1)=d2vxc(ipts,1)+v3rho3(1)
             d2vxc(ipts,2)=d2vxc(ipts,2)+v3rho3(2)
             d2vxc(ipts,3)=d2vxc(ipts,3)+v3rho3(3)
             d2vxc(ipts,4)=d2vxc(ipts,4)+v3rho3(4)
           endif
         endif
       case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
         select case(xc_f90_info_kind(funcs(i)%info))
         case(XC_EXCHANGE)
           if (nspden==1) then
             dvxc(ipts,1)=v2rho2(1)*two
             dvxc(ipts,2)=dvxc(ipts,1)
             dvxc(ipts,3)=two*two*vsigma(1)
             dvxc(ipts,4)=dvxc(ipts,3)
             dvxc(ipts,5)=four*two*v2rhosigma(1)
             dvxc(ipts,6)=dvxc(ipts,5)
             dvxc(ipts,7)=two*four*four*v2sigma2(1)
             dvxc(ipts,8)=dvxc(ipts,7)
           else
             dvxc(ipts,1)=v2rho2(1)
             dvxc(ipts,2)=v2rho2(3)
             dvxc(ipts,3)=two*vsigma(1)
             dvxc(ipts,4)=two*vsigma(3)
             dvxc(ipts,5)=two*v2rhosigma(1)
             dvxc(ipts,6)=two*v2rhosigma(6)
             dvxc(ipts,7)=four*v2sigma2(1)
             dvxc(ipts,8)=four*v2sigma2(6)
           end if
         case(XC_CORRELATION)
           if (nspden==1) then
             dvxc(ipts,9)=v2rho2(1)
             dvxc(ipts,10)=dvxc(ipts,9)
             dvxc(ipts,11)=dvxc(ipts,9)
             dvxc(ipts,12)=two*vsigma(1)
             dvxc(ipts,13)=two*v2rhosigma(1)
             dvxc(ipts,14)=dvxc(ipts,13)
             dvxc(ipts,15)=four*v2sigma2(1)
           else
             dvxc(ipts,9)=v2rho2(1)
             dvxc(ipts,10)=v2rho2(2)
             dvxc(ipts,11)=v2rho2(3)
             dvxc(ipts,12)=two*vsigma(1)
             dvxc(ipts,13)=two*v2rhosigma(1)
             dvxc(ipts,14)=two*v2rhosigma(6)
             dvxc(ipts,15)=four*v2sigma2(1)
           end if
         case(XC_EXCHANGE_CORRELATION)
           message=' KXC is not available for GGA XC_EXCHANGE_CORRELATION functionals from LibCXC'
           MSG_ERROR(message)
         end select
       end select
     end if

     if (any(funcs%family == XC_FAMILY_GGA) .or. any(funcs%family == XC_FAMILY_HYB_GGA)&
&        .or. any(funcs%family == XC_FAMILY_MGGA)) then
       !Convert the quantities returned by Libxc to the ones needed by ABINIT
       if (nspden == 1) then
         vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(1)*two
       else
         vxcgr(ipts,1) = vxcgr(ipts,1) + two*vsigma(1) - vsigma(2)
         vxcgr(ipts,2) = vxcgr(ipts,2) + two*vsigma(3) - vsigma(2)
         vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(2)
       end if
     end if
     if (any(funcs%family == XC_FAMILY_MGGA)) then
       vxclrho(ipts,1:nspden) = vxclrho(ipts,1:nspden) + vxclrhotmp(1:nspden)
       vxctau(ipts,1:nspden) = vxctau(ipts,1:nspden) + vxctautmp(1:nspden)
     end if

   end do

 end do

end subroutine libxc_functionals_getvxc
!----------------------------------------------------------------------

#endif

end module
!!***
