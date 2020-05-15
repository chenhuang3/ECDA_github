!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abi2big
!! NAME
!!  m_abi2big
!!
!! FUNCTION
!!  Module to copy objects from ABINIT to BigDFT and viceversa.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel,D. Caliste)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_abi2big
    
 use defs_basis
 use m_errors
 use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'm_abi2big'
!End of the abilint section

 implicit none
 
 private
 
 public :: wvl_vtrial_abi2big 
 !to copy vtrial to wvl_den%rhov and viceversa.

 public :: wvl_rho_abi2big
 !to copy a density from ABINIT to BigDFT or viceversa

 public :: wvl_rhov_abi2big
 !generic routine to copy a density or potential from/to 
 !ABINIT to/from BigDFT.

 public :: wvl_vxc_abi2big
 ! to copy Vxc from ABINIT to BigDFT and viceversa

 public :: wvl_vhartr_abi2big
 ! to copy V_hartree from ABINIT to BigDFT and viceversa

 public :: wvl_occ_abi2big
 ! to copy occupations from/to ABINIT to/from BigDFT
 
 public :: wvl_eigen_abi2big
 ! to copy eigenvalues from/to ABINIT to/from BigDFT

 public :: wvl_occopt_abi2big
 ! maps occupation method in ABINIT and BigDFT


 logical,parameter::hmem=.false. !high memory
!!  Set hmem=.false. if memory is limited. It will copy element by element.
!!  If  hmem=.true. all elements are copied at once.

contains
!!***

!!****f* m_abi2big/wvl_vtrial_abi2big
!! NAME
!!  wvl_vtrial_abi2big
!!
!! FUNCTION
!!  Copies vtrial in ABINIT to BigDFT objects and viceversa.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel, D. Caliste)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nfft= number of points (real space points) in vtrial.
!!  nspden= number of spin polarization densities.
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  vtrial(nfft,nspden)= trial potential (ABINIT array)
!!  wvl_den= density-potential BigDFT object
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!! vtrial is copied to wvl_den, or viceversa, depending on "opt" (see above).
!! It verifies that (or sets) wvl_den%rhov_is = KS_POTENTIAL.
!!
!! NOTES
!! It uses the generic routine wvl_rhov_abi2big.
!!
!! PARENTS
!!      afterscfloop,newvtr,rhotov,setvtr,wvl_psitohpsi
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_vtrial_abi2big(nfft,nspden,opt,vtrial,wvl_den)
    
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : KS_POTENTIAL
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_vtrial_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: nfft,nspden,opt
 real(dp) , intent(inout)  :: vtrial(nfft,nspden)
 type(wvl_denspot_type), intent(inout) :: wvl_den

!Local variables-------------------------------
! integer::ii ,ifft,ispden
 character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
 DBG_ENTER("COLL")

! write(message, '(a,a,a,a)' ) ch10, ' wvl_vtrial_abi2big : but why are you copying me :..o('
! call wrtout(std_out,message,'COLL')

#if defined HAVE_DFT_BIGDFT
 if(opt==1) then !ABINIT -> BIGDFT

   call wvl_rhov_abi2big(nfft,nspden,opt,vtrial,wvl_den%denspot%rhov)
   wvl_den%denspot%rhov_is = KS_POTENTIAL


 elseif(opt==2) then !BigDFT -> ABINIT

   if(wvl_den%denspot%rhov_is .ne. KS_POTENTIAL) then
     write(message, '(a,a)' ) ch10,&
&     ' wvl_vtrial_abi2big : rhov should contain the KS_POTENTIAL.'
     MSG_ERROR(message)
   end if

   call wvl_rhov_abi2big(nfft,nspden,opt,vtrial,wvl_den%denspot%rhov)

 else
   write(message, '(a,a)' ) ch10,&
&   ' wvl_vtrial_abi2big : wrong option.'
   MSG_ERROR(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
#endif

 DBG_EXIT("COLL")
end subroutine wvl_vtrial_abi2big
!!***

!!****f* m_abi2big/wvl_rho_abi2big
!! NAME
!!  wvl_rho_abi2big
!!
!! FUNCTION
!!  Copies the density from ABINIT to BigDFT, or viceversa.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel, D. Caliste)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nfft= number of points (real space points) in vtrial.
!!  nspden= number of spin polarization densities.
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  rhor(nfft,nspden)= trial potential (ABINIT array)
!!  wvl_den= density-potential BigDFT object
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! density copied from ABINIT to BigDFT or viceversa.
!! It verifies that (or sets) wvl_den%rhov_is= ELECTRONIC_DENSITY.
!!
!! NOTES
!! It uses the generic routine wvl_rhov_abi2big.
!!
!! PARENTS
!!      newrho,vtorho,wvl_mkrho
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_rho_abi2big(nfft,nspden,opt,rhor,wvl_den)
    
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : ELECTRONIC_DENSITY
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rho_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: nfft,nspden,opt
 real(dp) , intent(inout)  :: rhor(nfft,nspden)
 type(wvl_denspot_type), intent(inout) :: wvl_den

!Local variables-------------------------------
 !integer::ifft
 !real(dp)::tmp_up,tmp_dn,tmp_tot
 character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 DBG_ENTER("COLL")

! write(message, '(a,a,a,a)' ) ch10, ' wvl_rho_abi2big : but why are you copying me :..o('
! call wrtout(std_out,message,'COLL')

 if(nspden==4) then
   write(message, '(a)' ) ch10,&
&   ' wvl_rho_abi2big : nspden=4 not coded yet'
   MSG_ERROR(message)
 end if

#if defined HAVE_DFT_BIGDFT
 if(opt==1) then !ABINIT -> BIGDFT

   call wvl_rhov_abi2big(nfft,nspden,opt,rhor,wvl_den%denspot%rhov)
   wvl_den%denspot%rhov_is = ELECTRONIC_DENSITY


 elseif(opt==2) then !BigDFT -> ABINIT

   if(wvl_den%denspot%rhov_is .ne. ELECTRONIC_DENSITY) then
     write(message, '(a,a)' ) ch10,&
&     ' wvl_rho_abi2big : rhov should contain the ELECTRONIC_DENSITY.'
     MSG_ERROR(message)
   end if

   call wvl_rhov_abi2big(nfft,nspden,opt,rhor,wvl_den%denspot%rhov)
 else
   write(message, '(a,a)' ) ch10,&
&   ' wvl_rho_abi2big : wrong option.'
   MSG_ERROR(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
#endif


 DBG_EXIT("COLL")
end subroutine wvl_rho_abi2big
!!***

!!****f* m_abi2big/wvl_rhov_abi2big
!! NAME
!!  wvl_rhov_abi2big
!!
!! FUNCTION
!!  Copies a density/potential from ABINIT to BigDFT or viceversa.
!!  This is a generic routine to copy objects.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nfft= number of points (real space points) in vtrial.
!!  nspden= number of spin polarization densities.
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  rhov_abi(nfft,nspden) = density/potential array in ABINIT
!!  rhov_big(nfft,nspden) = density/potential array in BigDFT
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! At output rhov_abi is copied to rhov_abinit, or viceversa
!!
!! NOTES
!! This routine is duplicated:
!!                  This option is faster but it requires more memory.
!! Notice that we cannot point the variables since the spin convention is not
!! the same in BigDFT and ABINIT.
!! In ABINIT: index 1 is for the total spin (spin up + spin down) and index 2 is for spin up.
!! In BigDFT: indices 1 and 2 are for spin up and down, respectively.
!!
!! PARENTS
!!      m_abi2big,mklocl_wavelets,psolver_rhohxc
!!
!! CHILDREN
!!
!! SOURCE


subroutine wvl_rhov_abi2big(nfft,nspden,opt,rhov_abi,rhov_big)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rhov_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: nfft,nspden,opt
 real(dp):: rhov_big(nfft,nspden)
 real(dp):: rhov_abi(nfft,nspden)

!Local variables-------------------------------
 real(dp) :: tmpUp,tmpDown,tmpTot
 integer::ifft
 character(len=500) :: message                   ! to be uncommented, if needed
 real(dp),allocatable::rhoup(:),rhodn(:),rhotot(:)
 
! *************************************************************************
 DBG_ENTER("COLL")

! write(message, '(a,a,a,a)' ) ch10, ' wvl_rhov_abi2big : but why are you copying me :..o('
! call wrtout(std_out,message,'COLL')

 if(size(rhov_big)==1 .and. size(rhov_abi)==0) then
  return
  !no objects to copy
  ! In BigDFT by default they have size of 1!
 end if

 if(hmem) then
   if(nspden==2) then
     if(opt==1) then
       ABI_ALLOCATE(rhoup,(nfft))
       ABI_ALLOCATE(rhodn,(nfft))
     elseif(opt==2)  then
       ABI_ALLOCATE(rhotot,(nfft))
     end if
   end if
 end if

#if defined HAVE_DFT_BIGDFT
 if(opt==1) then !ABINIT -> BIGDFT
  if(nspden==2) then
    if(hmem) then
      rhoup(:)=rhov_abi(:,2)
      rhodn(:)=rhov_abi(:,1)-rhoup(:)
      rhov_big(:,1)=rhoup(:)
      rhov_big(:,2)=rhodn(:)
    else
      do ifft = 1, nfft
!       We change convention for BigDFT
        tmpDown=rhov_abi(ifft,1)-rhov_abi(ifft,2)
        tmpUp  =rhov_abi(ifft,2)
        rhov_big(ifft,1)=tmpUp
        rhov_big(ifft,2)=tmpDown
      end do
    end if !hmem
  else !nspden==1
    if(hmem) then
      rhov_big=rhov_abi
    else
      do ifft = 1, nfft
       rhov_big(ifft,1)=rhov_abi(ifft,1)
      end do
    end if!hmem
  end if !nspden

 elseif(opt==2) then !BigDFT -> ABINIT
  if(nspden==2) then
    if(hmem) then
      rhotot(:)=rhov_big(:,1)+rhov_big(:,2)
      rhov_abi(:,1)=rhotot(:)
      rhov_abi(:,2)=rhov_big(:,1)
    else 
      do ifft = 1, nfft
!       We change convention for BigDFT
        tmpTot=rhov_big(ifft,1)+rhov_big(ifft,2)
        rhov_abi(ifft,1)=tmpTot
        rhov_abi(ifft,2)=rhov_big(ifft,1) !Spin Up
      end do
    end if !hmem
  else if(nspden==1) then
    if(hmem) then
      rhov_abi=rhov_big
    else
      do ifft=1,nfft
        rhov_abi(ifft,1)=rhov_big(ifft,1)
      end do
    end if !hmem
  end if !nspden

 else
   write(message, '(a,a)' ) ch10,&
&   ' wvl_rhov_abi2big : wrong option.'
   MSG_ERROR(message)
 end if


 if(hmem) then
   if(nspden==2) then
     if(opt==1) then
       ABI_DEALLOCATE(rhoup)
       ABI_DEALLOCATE(rhodn)
     else if(opt==2) then
       ABI_DEALLOCATE(rhotot)
     end if !opt
   end if !nspden
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
#endif


 DBG_EXIT("COLL")
end subroutine wvl_rhov_abi2big
!!***

!!****f* m_abi2big/wvl_vxc_abi2big
!! NAME
!!  wvl_vxc_abi2big
!!
!! FUNCTION
!!  It copies the Vxc potential from ABINIT  to BigDFT or viceversa.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel,D. Caliste)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nfft= number of points (real space points) in vtrial.
!!  nspden= number of spin polarization densities.
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  vxc(nfft,nspden)= trial potential (ABINIT array)
!!  wvl_den= density-potential BigDFT object
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! vxc is copied to wvl_den, or viceversa, depending on "opt" (see above).
!!
!! NOTES
!! Vxc object in BigDFT has 4 indices.
!! The spin convention is not the same in both codes (see comments in wvl_rhov_abi2big).
!!
!! PARENTS
!!      psolver_rhohxc,wvl_psitohpsi
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_vxc_abi2big(nfft,nspden,opt,vxc,wvl_den)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_vxc_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: nfft,nspden,opt
 real(dp) , intent(inout)  :: vxc(nfft,nspden)
 type(wvl_denspot_type), intent(inout) :: wvl_den

!Local variables-------------------------------
 integer:: dim1,dim2,dim3,dim4,ii,i1,i2,i3
 integer:: size_abi, size_big
 real(dp):: vxc_up,vxc_dn,vxc_tot
 character(len=500) :: message
 
! *************************************************************************
 DBG_ENTER("COLL")

#if defined HAVE_DFT_BIGDFT

 dim1=size(wvl_den%denspot%v_xc,1)
 dim2=size(wvl_den%denspot%v_xc,2)
 dim3=size(wvl_den%denspot%v_xc,3)
 dim4=size(wvl_den%denspot%v_xc,4)

 size_abi=nspden*nfft
 size_big=dim1*dim2*dim3*dim4
 
 if(size_abi .ne. size_big) then
   if(size_abi == 0 .and. size_big == 1) then
     return
     !no objects to copy
     ! In BigDFT by default they have size of 1!
   end if
   write(message, '(a,a)') ch10,&
&    'wvl_vxc_abi2big: ABINIT and BigDFT objects do not have the same size'
   MSG_ERROR(message)
 end if
 if(nspden==4) then
   write(message, '(a,a)') ch10,&
&    'wvl_vxc_abi2big: nspden=4 not yet supported'
   MSG_ERROR(message)
 end if

 if(opt==1) then !ABINIT -> BIGDFT

  !wvl_den%denspot%V_XC=reshape(vxc,shape(wvl_den%denspot%V_XC))
  if(nspden==2) then
    ii=0
    do i3=1,dim3
      do i2=1,dim2
        do i1=1,dim1
          ii=ii+1
          vxc_tot=vxc(ii,1)
          vxc_up=vxc(ii,2)
          vxc_dn=vxc_tot-vxc_up
          wvl_den%denspot%v_xc(i1,i2,i3,1)=vxc_up
          wvl_den%denspot%v_xc(i1,i2,i3,2)=vxc_dn
        end do
      end do
    end do
  elseif(nspden==1) then
   ii=0
   do i3=1,dim3
     do i2=1,dim2
       do i1=1,dim1
         ii=ii+1 
         wvl_den%denspot%v_xc(i1,i2,i3,1)=vxc(ii,1)
       end do
     end do
   end do
  end if

 elseif(opt==2) then !BigDFT -> ABINIT


  if(nspden==2) then
    ii=0
    do i3=1,dim3
      do i2=1,dim2
        do i1=1,dim1
          ii=ii+1 
          vxc_up=wvl_den%denspot%v_xc(i1,i2,i3,1)
          vxc_dn=wvl_den%denspot%v_xc(i1,i2,i3,2)
          vxc_tot=vxc_dn+vxc_up
          vxc(ii,1)=vxc_tot
          vxc(ii,2)=vxc_up
        end do
      end do
    end do

  elseif(nspden==1) then

   !vxc=reshape(wvl_den%denspot%V_XC,shape(vxc))
   ii=0
   do i3=1,dim3
     do i2=1,dim2
       do i1=1,dim1
         ii=ii+1 !i1+dim1*(i2-1+dim2*(i3-1))
         vxc(ii,1)=wvl_den%denspot%v_xc(i1,i2,i3,1)
       end do
     end do
   end do

  else

    write(message, '(a,a,1x,i5)') ch10,&
&     'wvl_vxc_abi2big: invalid nspden=',nspden
    MSG_ERROR(message)

  end if !nspden

 else
   write(message, '(a,a)' ) ch10,&
&   ' wvl_vxc_abi2big : wrong option.'
   MSG_ERROR(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
#endif


 DBG_EXIT("COLL")
end subroutine wvl_vxc_abi2big
!!***

!!****f* m_abi2big/wvl_vhartr_abi2big
!! NAME
!!  wvl_vhartr_abi2big
!!
!! FUNCTION
!!  Copies vhartree in ABINIT to BigDFT objects and viceversa.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nfft= number of points (real space points) in vtrial.
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  vhartr(nfft)= Hartree potential (ABINIT array)
!!  wvl_den= density-potential BigDFT object
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!! vhartr is copied to wvl_den, or viceversa, depending on "opt" (see above).
!! It verifies that (or sets) wvl_den%rhov_is = HARTREE_POTENTIAL
!!
!! NOTES
!! It uses the generic routine wvl_rhov_abi2big.
!!
!! PARENTS
!!      psolver_rhohxc
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_vhartr_abi2big(nfft,opt,vhartr,wvl_den)
    
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : HARTREE_POTENTIAL
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_vhartr_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: nfft,opt
 real(dp) , intent(inout)  :: vhartr(nfft)
 type(wvl_denspot_type), intent(inout) :: wvl_den

!Local variables-------------------------------
! integer::ii !,ifft,ispden
 character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
 DBG_ENTER("COLL")

! write(message, '(a,a,a,a)' ) ch10, ' wvl_vhartr_abi2big : but why are you copying me :..o('
! call wrtout(std_out,message,'COLL')

#if defined HAVE_DFT_BIGDFT
 if(opt==1) then !ABINIT -> BIGDFT

   call wvl_rhov_abi2big(nfft,1,opt,vhartr,wvl_den%denspot%rhov(1:nfft))
   wvl_den%denspot%rhov_is = HARTREE_POTENTIAL


 elseif(opt==2) then !BigDFT -> ABINIT

   if(wvl_den%denspot%rhov_is .ne. HARTREE_POTENTIAL) then
     write(message, '(a,a)' ) ch10,&
&     ' wvl_vhartr_abi2big : rhov should contain the HARTREE_POTENTIAL.'
     MSG_ERROR(message)
   end if

   call wvl_rhov_abi2big(nfft,1,opt,vhartr,wvl_den%denspot%rhov(1:nfft))

 else
   write(message, '(a,a)' ) ch10,&
&   ' wvl_vhartr_abi2big : wrong option.'
   MSG_ERROR(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
#endif

 DBG_EXIT("COLL")
end subroutine wvl_vhartr_abi2big
!!***

!!****f* ABINIT/wvl_occ_abi2big
!! NAME
!!  wvl_occ_abi2big
!!
!! FUNCTION
!!  Copies occupations in ABINIT to BigDFT objects and viceversa.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  nsppol= number of spin polarization
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!! occ is copied to wfs%ks%orbs%occup, or viceversa, depending on "opt" (see above).
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop,gstate,vtorho,wvl_wfsinp_disk,wvl_wfsinp_scratch
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_occ_abi2big(mband,nkpt,nsppol,occ,opt,wvl_wfs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_occ_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: mband,nkpt,nsppol,opt
 real(dp) , intent(inout)  :: occ(mband*nkpt*nsppol)
 type(wvl_wf_type), intent(inout) :: wvl_wfs

!Local variables-------------------------------
 integer :: norb,norbd,norbu,ii
 character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
 DBG_ENTER("COLL")

!PENDING: I am not sure this will work for nsppol==2
!check also the parallel case.

#if defined HAVE_DFT_BIGDFT
 norbu=wvl_wfs%ks%orbs%norbu
 norbd=wvl_wfs%ks%orbs%norbd
 norb =wvl_wfs%ks%orbs%norb
 if(opt==1) then !ABINIT -> BIGDFT
   if (nsppol == 1) then
    do ii=1,norb
      wvl_wfs%ks%orbs%occup(ii)=occ(ii)
    end do
   else
     wvl_wfs%ks%orbs%occup(1:norbu)=occ(1:norbu)
     wvl_wfs%ks%orbs%occup(norbu + 1:norb)= &
&       occ(mband + 1:mband + norbd)
   end if
 elseif(opt==2) then !BigDFT -> ABINIT
   if (nsppol == 1) then
     do ii=1,norb
       occ=wvl_wfs%ks%orbs%occup
     end do
   else
     occ(1:norbu) = wvl_wfs%ks%orbs%occup(1:norbu)
     occ(mband + 1:mband + norbd) = &
&     wvl_wfs%ks%orbs%occup(norbu + 1:norb)
   end if
 else
   write(message, '(a,a)' ) ch10,&
&   ' wvl_occ_abi2big : wrong option.'
   MSG_ERROR(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
#endif

 DBG_EXIT("COLL")
end subroutine wvl_occ_abi2big
!!***

!!****f* ABINIT/wvl_eigen_abi2big
!! NAME
!!  wvl_eigen_abi2big
!!
!! FUNCTION
!!  Copies eigenvalues in ABINIT to BigDFT objects and viceversa.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  nsppol= number of spin polarization
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!! occ is copied to wfs%ks%orbs%occup, or viceversa, depending on "opt" (see above).
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop,vtorho
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_eigen_abi2big(mband,nkpt,nsppol,eigen,opt,wvl_wfs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_eigen_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: mband,nkpt,nsppol,opt
 real(dp) , intent(inout)  :: eigen(mband*nkpt*nsppol)
 type(wvl_wf_type), intent(inout) :: wvl_wfs

!Local variables-------------------------------
 integer :: ii,norb,norbd,norbu
 character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
 DBG_ENTER("COLL")

!PENDING: I am not sure this will work for nsppol==2
!check also the parallel case.

#if defined HAVE_DFT_BIGDFT
 norbu=wvl_wfs%ks%orbs%norbu
 norbd=wvl_wfs%ks%orbs%norbd
 norb =wvl_wfs%ks%orbs%norb
 if(opt==1) then !ABINIT -> BIGDFT
   if (nsppol == 1) then
    wvl_wfs%ks%orbs%eval=eigen
   else
     wvl_wfs%ks%orbs%eval(1:norbu)=eigen(1:norbu)
     wvl_wfs%ks%orbs%eval(norbu + 1:norb)= &
&       eigen(mband + 1:mband + norbd)
   end if
 elseif(opt==2) then !BigDFT -> ABINIT
   if (nsppol == 1) then
     do ii=1,norb
       eigen(ii)=wvl_wfs%ks%orbs%eval(ii)
     end do
   else
     eigen(1:norbu) = wvl_wfs%ks%orbs%eval(1:norbu)
     eigen(mband + 1:mband + norbd) = &
&     wvl_wfs%ks%orbs%eval(norbu + 1:norb)
   end if
 else
   write(message, '(a,a)' ) ch10,&
&   ' wvl_eigen_abi2big : wrong option.'
   MSG_ERROR(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
#endif

 DBG_EXIT("COLL")
end subroutine wvl_eigen_abi2big
!!***

!!****f* ABINIT/wvl_occopt_abi2big
!! NAME
!!  wvl_occopt_abi2big
!!
!! FUNCTION
!!  Copies occopt in ABINIT to BigDFT objects and viceversa.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Several smearing schemes do not exists in both codes such 
!! as the SMEARING_DIST_ERF in BigDFT.
!!
!! PARENTS
!!      vtorho,wvl_wfsinp_scratch
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

subroutine wvl_occopt_abi2big(occopt_abi,occopt_big,opt)
    
#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only : &
&  SMEARING_DIST_FERMI, SMEARING_DIST_COLD1, SMEARING_DIST_COLD2,&
&  SMEARING_DIST_METPX, SMEARING_DIST_ERF
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_occopt_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(inout)  :: occopt_abi,occopt_big
 integer , intent(in)     :: opt

!Local variables-------------------------------
 character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 
 DBG_ENTER("COLL")


#if defined HAVE_DFT_BIGDFT

 if(opt==1) then !ABINIT -> BIGDFT
   if(occopt_abi==3) then
     occopt_big=SMEARING_DIST_FERMI
   elseif(occopt_abi==4) then
     occopt_big=SMEARING_DIST_COLD1 
   elseif(occopt_abi==5) then
     occopt_big=SMEARING_DIST_COLD2
   elseif(occopt_abi==6) then
     occopt_big=SMEARING_DIST_METPX
   else
     write(message, '(4a)' ) ch10,&
&     ' wvl_occopt_abi2big : occopt does not have a corresponding option in BigDFT.',ch10,&
&     ' Action: change the value of occopt to a number between 3 and 6'
     MSG_ERROR(message)
   end if
 elseif(opt==2) then !BigDFT -> ABINIT
   if(occopt_big==SMEARING_DIST_FERMI) then
     occopt_abi=3
   elseif(occopt_big==SMEARING_DIST_COLD1) then
     occopt_abi=4
   elseif(occopt_big==SMEARING_DIST_COLD2) then
     occopt_abi=5
   elseif(occopt_big==SMEARING_DIST_METPX) then
     occopt_abi=6
   else
!    One should never get here.
     write(message, '(4a)' ) ch10,&
&     ' wvl_occopt_abi2big : occopt in BigDFT does not have a corresponding option in ABINIT.',ch10,&
&     ' Action: contact the ABINIT group'
     MSG_ERROR(message)
   end if
 else
   write(message, '(a,a)' ) ch10,&
&   ' wvl_occopt_abi2big : wrong option. Contact the ABINIT group'
   MSG_ERROR(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
#endif

 DBG_EXIT("COLL")
end subroutine wvl_occopt_abi2big

end module m_abi2big
!!***

