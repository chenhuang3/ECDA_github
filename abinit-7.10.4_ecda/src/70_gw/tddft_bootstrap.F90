!{\src2tex{textfont=tt}}
!!****f* ABINIT/tddft_bootstrap
!! NAME
!! tddft_bootstrap
!!
!! FUNCTION
!!  Compute the TDDFT Bootstrap kernel. Calculate RPA $\tilde\epsilon^{-1}$
!!  Use a special treatment of non-Analytic behavior of heads and wings in reciprocal space
!!  calculating these quantities for different small q-directions specified by the user
!!  (Not yet operative)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2014 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  iqibz=index of the q-point in the array Vcp%qibz where epsilon^-1 has to be calculated
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  npwe=Number of G-vectors in chi0.
!!  nI,nJ=Number of rows/columns in chi0_ij (1,1 in collinear case)
!!  nomega=Number of frequencies.
!!  omega(nomega)
!!  dim_wing=Dimension of the wings (0 or 3 if q-->0)
!!  chi0_head(dim_wing,dim_wing)=Head of of chi0 (only for q-->0)
!!  chi0_lwing(npwe*nI,dim_wing)=Lower wings of chi0 (only for q-->0)
!!  chi0_uwing(npwe*nJ,dim_wing)=Upper wings of chi0 (only for q-->0)
!!  comm=MPI communicator.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chi0(npwe*nI,npwe*nJ): in input the irreducible polarizability, in output 
!!   the symmetrized inverse dielectric matrix.
!! PARENTS
!!
!! CHILDREN
!!      atddft_symepsm1,rpa_symepsm1,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine tddft_bootstrap(iqibz,Vcp,npwe,nI,nJ,nomega,omega,chi0,nsteps,tolerance,&
&  my_nqlwl,dim_wing,chi0_head,chi0_lwing,chi0_uwing,epsm_lf,epsm_nlf,eelf,conv_err,comm)

 use defs_basis
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_vcoul
 use m_screening

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tddft_bootstrap'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nI,nJ,npwe,nomega,dim_wing,nsteps,my_nqlwl,comm
 real(dp),intent(in) :: tolerance
 real(dp),intent(out) :: conv_err
 type(vcoul_t),intent(in) :: Vcp
!arrays
 complex(dpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ,nomega)
 complex(dpc),intent(inout) :: chi0_lwing(npwe*nI,nomega,dim_wing)
 complex(dpc),intent(inout) :: chi0_uwing(npwe*nJ,nomega,dim_wing)
 complex(dpc),intent(inout) :: chi0_head(dim_wing,dim_wing,nomega)
 real(dp),intent(out) :: eelf(nomega,my_nqlwl)
 complex(dpc),intent(out) :: epsm_lf(nomega,my_nqlwl),epsm_nlf(nomega,my_nqlwl)

!Local variables-------------------------------
!scalars
 integer :: step,iw,g1,g2,comm_w,nprocs,ierr,option_test
 complex(gwpc) :: chi0w0_head
 logical :: isconverged
 !character(len=500) :: msg
!arrays
 complex(gwpc),allocatable :: fxc_boot(:,:),chi0_save(:,:,:)

! *************************************************************************

 nprocs = xcomm_size(comm)
 !
 ! MPI parallelization over frequencies.
 comm_w = comm

 chi0w0_head = chi0(1,1,1)

 ABI_MALLOC(chi0_save,(npwe*nI,npwe*nJ,nomega))
 ABI_CHECK_ALLOC("out of memory in fxc_boot")
 chi0_save = chi0

 ABI_MALLOC(fxc_boot,(npwe*nI,npwe*nJ))
 ABI_CHECK_ALLOC("out of memory in fxc_boot")
 !
 ! Compute RPA e^{-1} from chi0.
 do iw=1,nomega
   call rpa_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,iw),my_nqlwl,dim_wing,&
&    chi0_head(:,:,iw),chi0_lwing(:,iw,:),chi0_uwing(:,iw,:),epsm_lf(iw,:),epsm_nlf(iw,:),eelf(iw,:),comm_w)
 end do

 do step=1,nsteps
   !
   ! Build bootstrap Kernel.
   fxc_boot = chi0(:,:,1) / chi0w0_head
   !
   ! Compute new e^{-1} within TDDFT.
   chi0 = chi0_save

   epsm_lf  = zero
   epsm_nlf = zero
   eelf     = zero

   do iw=1,nomega
     option_test=0 ! TESTPARTICLE
     call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,iw),fxc_boot,option_test,my_nqlwl,dim_wing,omega(iw),&
&      chi0_head(:,:,iw),chi0_lwing(:,iw,:),chi0_uwing(:,iw,:),epsm_lf(iw,:),epsm_nlf(iw,:),eelf(iw,:),comm_w)
   end do

   if (nprocs > 1) then
     call xmpi_sum(epsm_lf,comm,ierr)
     call xmpi_sum(epsm_nlf,comm,ierr)
     call xmpi_sum(eelf,comm,ierr)
   end if
   !
   ! Check for convergence.
   conv_err = smallest_real
   do iw=1,nomega
     do g2=1,npwe*nJ
       do g1=1,npwe*nI
         conv_err= MAX(conv_err, ABS(chi0(g1,g2,iw) - chi0_save(g1,g2,iw)) )
       end do 
     end do
   end do
   isconverged = (conv_err <= tolerance)
   !
   ! Compute optical properties.
   if (isconverged .or. step == nsteps) then
     ! Gather final epsm-1 on each node.
     do iw=1,nomega
       call xmpi_sum(chi0(:,:,iw),comm,ierr)
     end do
     exit
   end if
 end do

 if (.not. isconverged) then

 else 
   conv_err = -conv_err
 end if

 ABI_FREE(chi0_save)
 ABI_FREE(fxc_boot)

end subroutine tddft_bootstrap
!!***
