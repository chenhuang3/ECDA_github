!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_screening
!! NAME
!!  m_screening
!!
!! FUNCTION
!!  This module contains the definition of the object used to deal
!!  with the inverse dielectric matrix as well as related methods.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_screening

 use defs_basis
 use defs_abitypes
 use m_blas
 use m_linalg_interfaces
 use m_xmpi
 use m_errors
 use m_copy
 use m_splines
 use m_profiling_abi
 use m_lebedev

 use m_gwdefs,          only : GW_TOLQ0, czero_gw
 use m_fstrings,        only : toupper
 use m_io_tools,        only : open_file
 use m_numeric_tools,   only : print_arr, hermitianize
 use m_geometry,        only : normv
 use m_abilasi,         only : xginv
 use m_crystal,         only : crystal_t
 use m_bz_mesh,         only : kmesh_t, get_BZ_item, box_len
 use m_fft_mesh,        only : g2ifft
 use m_gsphere,         only : gsphere_t
 use m_vcoul,           only : vcoul_t
 use m_io_screening,    only : free_scrhdr, scr_hdr_io, read_screening, write_screening, &
&                              copy_scrhdr, HSCR_LATEST_HEADFORM, scrhdr_type, read_pole_screening
 use m_model_screening, only : re_and_im_screening_with_phase
 use m_spectra,         only : spectra_type, init_spectra, destroy_spectra, repr_dielconst
 use m_sphharm,         only : ylmc
 use m_mpinfo,          only : destroy_mpi_enreg 

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_screening/epsilonm1_results
!! NAME
!! epsilonm1_results
!!
!! FUNCTION
!! For the GW part of ABINIT, the epsilonm1_results structured datatype
!! gather the results of screening : the inverse dielectric matrix,
!! and the omega matrices .
!!
!! SOURCE

 type,public :: Epsilonm1_results

  integer :: ID                          ! Matrix identifier: O if not yet defined, 1 for chi0,
                                         ! 2 for chi, 3 for epsilon, 4 for espilon^{-1}, 5 for W.
  integer :: ikxc                        ! Kxc kernel used, 0 for None (RPA), >0 for static TDDFT (=ixc), <0 for TDDFT
  integer :: fform                       ! File format: 1002 for SCR|SUSC files. 2002 for pole fit
  integer :: mqmem                       ! =0 for out-of-core solution, =nqibz if entire matrix is stored in memory.
  integer :: nI,nJ                       ! Number of components (rows,columns) in chi|eps^-1. (1,1) if collinear.
  integer :: nqibz                       ! Number of q-points in the IBZ used.
  integer :: nqlwl                       ! Number of point used for the treatment of the long wave-length limit.
  integer :: nomega                      ! Number of frequencies used.
  integer :: nomega_i                    ! Number of purely imaginary frequencies used.
  integer :: nomega_r                    ! Number of real frequencies used.
  integer :: nomega_c                    ! Number of frequencies in the complex plane.
  integer :: npoles                      ! Number of poles for pole-fit screening
  integer :: ncoeff                      ! Number of coefficients = npoles*3+1 for phase
  integer :: npwe                        ! Number of G vectors used.
  integer :: test_type                   ! 0 for None, 1 for TEST-PARTICLE, 2 for TEST-ELECTRON (only for TDDFT)
  integer :: Tordering                   ! 0 if not defined, 1 for Time-Ordered, 2 for Advanced, 3 for Retarded.

  character(len=fnlen) :: fname          ! Name of the file from which epsm1 is read.

!arrays
  integer,allocatable :: gvec(:,:)
  ! gvec(3,npwe)
  ! G-vectors used to describe the two-point function (r.l.u.).

  real(dp),allocatable :: qibz(:,:)
  ! qibz(3,nqibz)
  ! q-points in reduced coordinates

  real(dp),allocatable :: qlwl(:,:)
  ! qlwl(3,nqlwl)
  ! q-points used for the long wave-length limit treatment.

  real(gwp),allocatable :: epsm1_pole(:,:,:,:)  
  ! epsm1(npwe,npwe,nomega,nqibz)
  ! Contains the two-point function $\epsilon_{G,Gp}(q,omega)$ in frequency and reciprocal space.

  complex(gwpc),allocatable :: epsm1(:,:,:,:) 
  ! epsm1(npwe,npwe,nomega,nqibz)
  ! Contains the two-point function $\epsilon_{G,Gp}(q,omega)$ in frequency and reciprocal space.

  complex(dpc),allocatable :: lwing(:,:,:) 
  ! lwing(npwe,nomega,nqlwl)
  ! Lower wings for the different q"s -->0

  complex(dpc),allocatable :: omega(:)
  ! omega(nomega)
  ! Frequencies used both along the real and the imaginary axis.

  complex(dpc),allocatable :: uwing(:,:,:)
  ! uwing(npwe,nomega,nqlwl)
  ! Upper wings for the different q"s -->0

  type(ScrHdr_type) :: Hscr
  ! The header reported in the _SCR of _SUSC file.
  ! This object contains information on the susceptibility or the inverse dielectric matrix
  ! as stored in the external file. These quantities do *NOT* correspond to the quantities
  ! used during the GW calculation since some parameters might differ, actually they might be smaller.
  ! For example, the number of G-vectors used can be smaller than the number of G"s stored on file.

 end type Epsilonm1_results

 public :: destroy_epsilonm1_results      ! Free all associated pointers
 public :: print_epsilonm1_results        ! Print basic info
 public :: Epsm1_symmetrizer              ! Symmetrize two-point function at a q-point in the BZ.
 public :: Epsm1_symmetrizer_inplace      ! In-place version of the above
 public :: Epsm1_pole_symmetrizer         ! Symmetrize pole fit screening
 public :: Epsm1_pole_symmetrizer_inplace ! In-place version of the above
 public :: init_Er_from_file              ! Initialize the object from file
 public :: mkdump_Er                      ! Dump the object to a file.
 public :: get_epsm1
 public :: get_pole_epsm1
 public :: decompose_epsm1
 public :: make_epsm1_driver
 public :: make_W                         ! Calculate W from the data stored in Er (in place, content of Er% is changed).
 public :: mkem1_q0
 public :: screen_mdielf                  ! Calculates W_{G,G'}(q,w) for a given q-point in the BZ using a model dielectric function.
 public :: recalculate_epsm1_freq_grid
 public :: rpa_symepsm1
!!***

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_screening/destroy_epsilonm1_results
!! NAME
!! destroy_epsilonm1_results
!!
!! FUNCTION
!! Deallocate all the pointers in Er that result to be associated.
!! Perform also a cleaning of the Header.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening,mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_epsilonm1_results(Er)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_epsilonm1_results'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_results),intent(inout) :: Er
! *************************************************************************

 !@Epsilonm1_results
 !integer
 if (allocated(Er%gvec)) then
   ABI_FREE(Er%gvec)
 end if

 !real
 if (allocated(Er%qibz)) then
   ABI_FREE(Er%qibz)
 end if
 if (allocated(Er%qlwl)) then
   ABI_FREE(Er%qlwl)
 end if
 if (allocated(Er%epsm1_pole)) then
   ABI_FREE(Er%epsm1_pole)
 end if

 !complex
 if (allocated(Er%epsm1)) then
   ABI_FREE(Er%epsm1)
 end if
 if (allocated(Er%lwing)) then
   ABI_FREE(Er%lwing)
 end if
 if (allocated(Er%omega)) then
   ABI_FREE(Er%omega)
 end if
 if (allocated(Er%uwing)) then
   ABI_FREE(Er%uwing)
 end if

 !datatypes
 call free_scrhdr(Er%Hscr)

end subroutine destroy_epsilonm1_results
!!***

!----------------------------------------------------------------------

!!****f* m_screening/print_epsilonm1_results
!! NAME
!!  print_epsilonm1_results
!!
!! FUNCTION
!! Print the basic dimensions and the most important
!! quantities reported in the Epsilonm1_results data type.
!!
!! INPUTS
!!  Er<Epsilonm1_results>=The data type.
!!  unit[optional]=the unit number for output.
!!  prtvol[optional]=verbosity level.
!!  mode_paral[optional]=either COLL or PERS.
!!
!! OUTPUT
!!  Only printing.
!!
!! PARENTS
!!      m_screening,mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_epsilonm1_results(Er,unit,prtvol,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_epsilonm1_results'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 type(Epsilonm1_results),intent(in) :: Er

!Local variables-------------------------------
 integer :: iomega,iqibz,iqlwl,unt,verbose,rdwr
 character(len=50) :: rfname,rforder,rfapprox,rftest,kxcname
 character(len=500) :: msg
 character(len=4) :: mode
! *************************************************************************

 unt    =std_out; if (PRESENT(unit      )) unt    =unit
 verbose=0      ; if (PRESENT(prtvol    )) verbose=prtvol
 mode   ='COLL' ; if (PRESENT(mode_paral)) mode   =mode_paral

 ! === chi0 or \epsilon^{-1} ? ===
 SELECT CASE (Er%ID)
 CASE (0)
   rfname='Undefined'
 CASE (1)
   rfname='Irreducible Polarizability'
 CASE (2)
   rfname='Polarizability'
 CASE (3)
   rfname='Symmetrical Dielectric Matrix'
 CASE (4)
   rfname='Symmetrical Inverse Dielectric Matrix'
 CASE DEFAULT
   write(msg,'(a,i3)') 'Wrong value of Er%ID= ',Er%ID
   MSG_BUG(msg)
 END SELECT

 ! === For chi, \espilon or \epsilon^{-1}, define the approximation ===
 rfapprox='None'
 if (Er%ID>=2.or.Er%ID<=4) then
   if (Er%ikxc==0) then
     rfapprox='RPA'
   else if (Er%ikxc>0) then
     rfapprox='Static TDDFT'
   else
     rfapprox='TDDFT'
   end if
 end if

 ! === If TDDFT and \epsilon^{-1}, define the type ===
 rftest='None'
! if (Er%ID==0) then
!  if (Er%test_type==0) then
!   rftest='TEST-PARTICLE'
!  else if (Er%test_type==1) then
!   rftest='TEST-ELECTRON'
!  else
!   write(msg,'(4a,i3)')ch10,&
!&   ' print_epsilonm1_results : BUG - ',ch10,&
!&   ' Wrong value of Er%test_type = ',Er%test_type
!   MSG_ERROR(msg)
!  end if
! end if

 ! === Define time-ordering ===
 rforder='Undefined'
 if (Er%Tordering==1) then
   rforder='Time-Ordered'
 else if (Er%Tordering==2) then
   rforder='Advanced'
 else if (Er%Tordering==3) then
   rforder='Retarded'
 else
   write(msg,'(a,i3)')'Wrong value of Er%Tordering= ',Er%Tordering
   MSG_BUG(msg)
 end if

 kxcname='None'
 if (Er%ikxc/=0) then
   !TODO Add function to retrieve kxc name
   MSG_ERROR('Add function to retrieve kxc name')
   kxcname='XXXXX'
 end if

 write(msg,'(6a,5(3a))')ch10,&
&  ' ==== Info on the Response Function ==== ',ch10,&
&  '  Associated File ................  ',TRIM(Er%fname),ch10,&
&  '  Response Function Type .......... ',TRIM(rfname),ch10,&
&  '  Type of Approximation ........... ',TRIM(rfapprox),ch10,&
&  '  XC kernel used .................. ',TRIM(kxcname),ch10,&
&  '  Type of probing particle ........ ',TRIM(rftest),ch10,&
&  '  Time-Ordering ................... ',TRIM(rforder),ch10
 call wrtout(unt,msg,mode)
 write(msg,'(a,2i4,a,3(a,i4,a),a,3i4,2a,i4,a)')&
&  '  Number of components ............ ',Er%nI,Er%nJ,ch10,&
&  '  Number of q-points in the IBZ ... ',Er%nqibz,ch10,&
&  '  Number of q-points for q-->0 .... ',Er%nqlwl,ch10,&
&  '  Number of G-vectors ............. ',Er%npwe,ch10,&
&  '  Number of frequencies ........... ',Er%nomega,Er%nomega_r,Er%nomega_i,ch10,&
&  '  Value of mqmem .................. ',Er%mqmem,ch10
 call wrtout(unt,msg,mode)
 if (Er%npoles>0) then ! We have a pole-fit
   write(msg,'(3a,i4,2a,i4,3a)')&
&    '  - THIS IS A POLE-FIT SCREENING -  ',ch10,&
&    '  Number of poles ................. ',Er%npoles,ch10,&
&    '  Number of coefficients .......... ',Er%ncoeff,ch10,&
&    '  (npoles*3 + 1 value for phase)... ',ch10
   call wrtout(unt,msg,mode)
 end if

 if (Er%nqlwl/=0) then
   write(msg,'(a,i3)')' q-points for long wavelength limit: ',Er%nqlwl
   call wrtout(unt,msg,mode)
   do iqlwl=1,Er%nqlwl
     write(msg,'(1x,i5,a,3es16.8)')iqlwl,') ',Er%qlwl(:,iqlwl)
     call wrtout(unt,msg,mode)
   end do
 end if

 if (verbose>0) then ! Print out head and wings in the long-wavelength limit.
   ! TODO add additional stuff.
   if (Er%nqlwl>0) then
     write(msg,'(1x,2a)')' Heads and wings of chi0(G,G'')',ch10
     call wrtout(unt,msg,mode)
     do iqlwl=1,Er%nqlwl
       write(msg,'(1x,a,i2,a)')' chi0(qlwl =',iqlwl,')'
       call wrtout(unt,msg,mode)
       do iomega=1,Er%nomega
         write(msg,'(2x,a,i4,a,2f9.4,a)')&
&          ' Upper and lower wings at the ',iomega,' th omega',Er%omega(iomega)*Ha_eV,' [eV]'
         call wrtout(unt,msg,mode)
         call print_arr(Er%uwing(:,iomega,iqlwl),max_r=9,unit=unt)
         call print_arr(Er%lwing(:,iomega,iqlwl),max_r=9,unit=unt)
       end do
     end do
   end if

   write(msg,'(a,i4)')' Calculated Frequencies: ',Er%nomega
   call wrtout(unt,msg,mode)
   do iomega=1,Er%nomega
     write(msg,'(i4,es14.6)')iomega,Er%omega(iomega)*Ha_eV
     call wrtout(unt,msg,mode)
   end do

   write(msg,'(a,i4)')' Calculated q-points: ',Er%nqibz
   call wrtout(unt,msg,mode)
   do iqibz=1,Er%nqibz
     write(msg,'(1x,i4,a,3es16.8)')iqibz,') ',Er%qibz(:,iqibz)
     call wrtout(unt,msg,mode)
   end do

   rdwr=4
   !%call hdr_io_int(Er%fform,Er%Hscr%Hdr,rdwr,unt)
 end if ! verbose>0

end subroutine print_epsilonm1_results
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_symmetrizer
!! NAME
!!  Epsm1_symmetrizer
!!
!! FUNCTION
!!  Symmetrize the inverse dielectric matrix, namely calculate epsilon^{-1} at a generic
!!  q-point in the BZ starting from the knowledge of the matrix at a q-point in the IBZ.
!!  The procedure is quite generic and can be used for every two-point function which has
!!  the same symmetry as the crystal.
!!
!! INPUTS
!!  nomega=Number of frequencies required. All frequencies from 1 up to nomega are symmetrized.
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  remove_exchange=If .TRUE., return e^{-1}-1 namely remove the exchange part.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!    %grottb
!!    %phmSGt
!!  Qmesh<kmesh_t>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi)
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required.
!!
!! OUTPUT
!!  epsm1_qbz(npwc,npwc,nomega)=The inverse dielectric matrix at the q-point defined by iq_bz.
!!   Exchange part can be subtracted out.
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!!
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine Epsm1_symmetrizer(iq_bz,nomega,npwc,Er,Gsph,Qmesh,remove_exchange,epsm1_qbz)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Epsm1_symmetrizer'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npwc
 logical,intent(in) :: remove_exchange
 type(Epsilonm1_results),intent(in) :: Er
 type(gsphere_t),target,intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
!arrays
 complex(gwpc),intent(out) :: epsm1_qbz(npwc,npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: iomega,ii,jj,iq_ibz,itim_q,isym_q,iq_loc,sg1,sg2
 complex(gwpc) :: phmsg1t,phmsg2t_star
!arrays
 real(dp) :: qbz(3)

! *********************************************************************

 ABI_CHECK(Er%nomega>=nomega,'Too many frequencies required')
 ABI_CHECK(Er%npwe  >=npwc , 'Too many G-vectors required')

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled.
 iq_loc=iq_ibz; if (Er%mqmem==0) iq_loc=1

 ! MG: grottb is a 1-1 mapping, hence we can collapse the loops (false sharing is not an issue here).
 !grottb => Gsph%rottb (1:npwc,itim_q,isym_q)
 !phmSgt => Gsph%phmSGt(1:npwc,isym_q)

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(sg2,sg1,phmsg1t,phmsg2t_star)
 do iomega=1,nomega
   do jj=1,npwc
     sg2 = Gsph%rottb(jj,itim_q,isym_q)  
     phmsg2t_star = CONJG(Gsph%phmSGt(jj,isym_q))
     do ii=1,npwc
       sg1 = Gsph%rottb(ii,itim_q,isym_q)  
       phmsg1t = Gsph%phmSGt(ii,isym_q)
       epsm1_qbz(sg1,sg2,iomega) = Er%epsm1(ii,jj,iomega,iq_loc) * phmsg1t * phmsg2t_star
       !epsm1_qbz(sg1,sg2,iomega) = Er%epsm1(ii,jj,iomega,iq_loc) * phmSgt(ii) * CONJG(phmSgt(jj))
     end do
   end do
 end do
 !
 ! === Account for time-reversal ===
 !epsm1_qbz(:,:,iomega)=TRANSPOSE(epsm1_qbz(:,:,iomega))
 if (itim_q==2) then
!$OMP PARALLEL DO IF (nomega > 1)
   do iomega=1,nomega
     call sqmat_itranspose(npwc,epsm1_qbz(:,:,iomega))
   end do
 end if

 if (remove_exchange) then
   ! === Subtract the exchange contribution ===
   ! If it's a pole screening, the exchange contribution is already removed
!$OMP PARALLEL DO IF (nomega > 1)
   do iomega=1,nomega
     do ii=1,npwc
       epsm1_qbz(ii,ii,iomega)=epsm1_qbz(ii,ii,iomega)-CMPLX(1.0_gwp,0.0_gwp)
     end do
   end do
 endif

end subroutine Epsm1_symmetrizer
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_symmetrizer_inplace
!! NAME
!!  Epsm1_symmetrizer_inplace
!!
!! FUNCTION
!!  Same function as Epsm1_symmetrizer, ecept now the array Ep%epsm1 is modified inplace
!!  thorugh an auxiliary work array of dimension (npwc,npwc)
!!
!! INPUTS
!!  nomega=Number of frequencies required. All frequencies from 1 up to nomega are symmetrized.
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  remove_exchange=If .TRUE., return e^{-1}-1 namely remove the exchange part.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!    %grottb
!!    %phmSGt
!!  Qmesh<kmesh_t>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi)
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required.
!!
!! OUTPUT
!!  Er%epsm1(npwc,npwc,nomega,iq_loc) symmetrised
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!!
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine Epsm1_symmetrizer_inplace(iq_bz,nomega,npwc,Er,Gsph,Qmesh,remove_exchange)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Epsm1_symmetrizer_inplace'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npwc
 logical,intent(in) :: remove_exchange
 type(Epsilonm1_results) :: Er
 type(gsphere_t),target,intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh

!Local variables-------------------------------
!scalars
 integer :: iomega,ii,jj,iq_ibz,itim_q,isym_q,iq_loc,sg1,sg2
!arrays
 real(dp) :: qbz(3)
 complex(gwpc) :: phmsg1t,phmsg2t_star
 complex(gwpc),allocatable :: work(:,:)

! *********************************************************************

 ABI_CHECK(Er%nomega>=nomega,'Too many frequencies required')
 ABI_CHECK(Er%npwe  >=npwc , 'Too many G-vectors required')

 ABI_MALLOC(work,(npwc,npwc))

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled.
 iq_loc=iq_ibz; if (Er%mqmem==0) iq_loc=1

 !grottb => Gsph%rottb (1:npwc,itim_q,isym_q)
 !phmSgt => Gsph%phmSGt(1:npwc,isym_q)

!$OMP PARALLEL DO PRIVATE(sg2,sg1,phmsg1t,phmsg2t_star) IF (nomega > 1)
 do iomega=1,nomega
   do jj=1,npwc
     sg2 = Gsph%rottb(jj,itim_q,isym_q)  
     phmsg2t_star = CONJG(Gsph%phmSGt(jj,isym_q))
     do ii=1,npwc
       sg1 = Gsph%rottb(ii,itim_q,isym_q)  
       phmsg1t = Gsph%phmSGt(ii,isym_q)
       work(sg1,sg2) = Er%epsm1(ii,jj,iomega,iq_loc) * phmsg1t * phmsg2t_star
       !work(grottb(ii),grottb(jj))=Er%epsm1(ii,jj,iomega,iq_loc)*phmSgt(ii)*CONJG(phmSgt(jj))
     end do
   end do
   Er%epsm1(:,:,iomega,iq_loc) = work(:,:)
 end do
 !
 ! === Account for time-reversal ===
 !Er%epsm1(:,:,iomega,iq_loc)=TRANSPOSE(Er%epsm1(:,:,iomega,iq_loc))
 if (itim_q==2) then
!$OMP PARALLEL DO IF (nomega > 1)
   do iomega=1,nomega
     call sqmat_itranspose(npwc,Er%epsm1(:,:,iomega,iq_loc))
   end do
 end if

 ! === Subtract the exchange contribution ===
 if (remove_exchange) then
!$OMP PARALLEL DO IF (nomega > 1)
   do iomega=1,nomega
     do ii=1,npwc
       Er%epsm1(ii,ii,iomega,iq_loc)=Er%epsm1(ii,ii,iomega,iq_loc)-1.0_gwp
     end do
   end do
 endif

 ABI_FREE(work)

end subroutine Epsm1_symmetrizer_inplace
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_pole_symmetrizer
!! NAME
!!  Epsm1_symmetrizer
!!
!! FUNCTION
!!  Symmetrize the inverse dielectric matrix, namely calculate epsilon^{-1} at a generic
!!  q-point in the BZ starting from the knowledge of the matrix at a q-point in the IBZ.
!!  The procedure is quite generic and can be used for every two-point function which has
!!  the same symmetry as the crystal. This routine is meant for pole-fit screening
!!
!! INPUTS
!!  ncoeff=Number of coefficients in pole-fit + phase value
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!    %grottb
!!    %phmSGt
!!  Qmesh<kmesh_t>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi)
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required.
!!
!! OUTPUT
!!  epsm1_pole_qbz(npwc,npwc,nomega)=The inverse dielectric matrix at the q-point defined by iq_bz.
!!   Exchange part can be subtracted out.
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!!
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! SOURCE

subroutine Epsm1_pole_symmetrizer(iq_bz,ncoeff,npwc,Er,Gsph,Qmesh,epsm1_pole_qbz)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Epsm1_pole_symmetrizer'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,ncoeff,npwc
 type(Epsilonm1_results),intent(in) :: Er
 type(gsphere_t),target,intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
!arrays
 real(gwp),intent(out) :: epsm1_pole_qbz(npwc,npwc,ncoeff)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,iq_ibz,itim_q,isym_q,iq_loc,icoeff
 real(gwp) :: phase
 complex   :: cphase
!arrays
 integer,pointer :: grottb(:)
 real(dp) :: qbz(3)
 complex(gwpc),pointer :: phmSgt(:)

! *********************************************************************

 ABI_CHECK(Er%npwe>=npwc , 'Too many G-vectors required')

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 grottb => Gsph%rottb (1:npwc,itim_q,isym_q)
 phmSgt => Gsph%phmSGt(1:npwc,isym_q)

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled.
 iq_loc=iq_ibz ; if (Er%mqmem==0) iq_loc=1

 do jj=1,npwc
   do ii=1,npwc
     do icoeff=1,ncoeff-1,3
       epsm1_pole_qbz(grottb(ii),grottb(jj),icoeff)=Er%epsm1_pole(ii,jj,icoeff,iq_loc)
       epsm1_pole_qbz(grottb(ii),grottb(jj),icoeff+1)=Er%epsm1_pole(ii,jj,icoeff+1,iq_loc)
       epsm1_pole_qbz(grottb(ii),grottb(jj),icoeff+2)=Er%epsm1_pole(ii,jj,icoeff+2,iq_loc)
     end do
     ! Add phase
     cphase = phmSgt(ii)*CONJG(phmSgt(jj))
     phase  = ATAN(AIMAG(cphase)/REAL(cphase))
     epsm1_pole_qbz(grottb(ii),grottb(jj),ncoeff)=Er%epsm1_pole(ii,jj,ncoeff,iq_loc)-phase
   end do
 end do
 !
 ! === Account for time-reversal ===
 if (itim_q==2) then
   do icoeff=1,ncoeff
     !epsm1_pole_qbz(:,:,icoeff)=TRANSPOSE(epsm1_pole_qbz(:,:,icoeff))
     call sqmat_itranspose(npwc,epsm1_pole_qbz(:,:,icoeff))
   end do
 end if

end subroutine Epsm1_pole_symmetrizer
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_pole_symmetrizer_inplace
!! NAME
!!  Epsm1_symmetrizer
!!
!! FUNCTION
!!  Symmetrize the inverse dielectric matrix, namely calculate epsilon^{-1} at a generic
!!  q-point in the BZ starting from the knowledge of the matrix at a q-point in the IBZ.
!!  The procedure is quite generic and can be used for every two-point function which has
!!  the same symmetry as the crystal. This routine is meant for pole-fit screening
!!
!! INPUTS
!!  ncoeff=Number of coefficients in pole-fit + phase value
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!    %grottb
!!    %phmSGt
!!  Qmesh<kmesh_t>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi)
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required.
!!
!! OUTPUT
!!  epsm1_pole_qbz(npwc,npwc,nomega)=The inverse dielectric matrix at the q-point defined by iq_bz.
!!   Exchange part can be subtracted out.
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!!
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! SOURCE

subroutine Epsm1_pole_symmetrizer_inplace(iq_bz,ncoeff,npwc,Er,Gsph,Qmesh)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Epsm1_pole_symmetrizer_inplace'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,ncoeff,npwc
 type(Epsilonm1_results),intent(inout) :: Er
 type(gsphere_t),target,intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
!arrays

!Local variables-------------------------------
!scalars
 integer :: ii,jj,iq_ibz,itim_q,isym_q,iq_loc,icoeff
 real(gwp) :: phase
 complex(gwpc) :: cphase
!arrays
 integer,pointer :: grottb(:)
 real(dp) :: qbz(3)
 complex(gwpc),pointer :: phmSgt(:)
 real(gwp), allocatable :: work(:,:)

! *********************************************************************

 ABI_CHECK(Er%npwe>=npwc, 'Too many G-vectors required')

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 ABI_MALLOC(work,(npwc,npwc))

 grottb => Gsph%rottb (1:npwc,itim_q,isym_q)
 phmSgt => Gsph%phmSGt(1:npwc,isym_q)

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled.
 iq_loc=iq_ibz ; if (Er%mqmem==0) iq_loc=1

 do icoeff=1,ncoeff-1
   do jj=1,npwc
     do ii=1,npwc
         work(grottb(ii),grottb(jj))=Er%epsm1_pole(ii,jj,icoeff,iq_loc)
     end do
   end do
   Er%epsm1_pole(:,:,icoeff,iq_loc)=work(:,:)
 end do
 ! Add phase
 do jj=1,npwc
   do ii=1,npwc
     cphase = phmSgt(ii)*CONJG(phmSgt(jj))
     phase  = ATAN(AIMAG(cphase)/REAL(cphase))
     work(grottb(ii),grottb(jj))=Er%epsm1_pole(ii,jj,ncoeff,iq_loc)+phase
   end do
 end do
 Er%epsm1_pole(:,:,ncoeff,iq_loc) = work(:,:)
 !
 ! === Account for time-reversal ===
 if (itim_q==2) then
   do icoeff=1,ncoeff
     !Er%epsm1_pole(:,:,icoeff,iq_loc)=TRANSPOSE(Er%epsm1_pole(:,:,icoeff,iq_loc))
     call sqmat_itranspose(npwc,Er%epsm1_pole(:,:,icoeff,iq_loc))
   end do
 end if

 ABI_FREE(work)

end subroutine Epsm1_pole_symmetrizer_inplace
!!***

!----------------------------------------------------------------------

!!****f* m_screening/init_Er_from_file
!! NAME
!!  init_Er_from_file
!!
!! FUNCTION
!!  Initialize basic dimensions and the important (small) arrays in an Epsilonm1_results data type
!!  starting from a file containing either epsilon^{-1} (_SCR) or chi0 (_SUSC).
!!
!! INPUTS
!!  accesswff=Option defining the file format of the external file.
!!  fname=The name of the external file used to read the matrix.
!!  mqmem=0 for out-of-core solution, /=0 if entire matrix has to be stored in memory.
!!  npwe_asked=Number of G-vector to be used in the calculation, if <=0 use Max allowed number.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Er<Epsilonm1_results>=The structure initialized with basic dimensions and arrays.
!!
!! PARENTS
!!      m_screening,mrgscr,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_Er_from_file(Er,fname,mqmem,npwe_asked,accesswff,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_Er_from_file'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: mqmem,accesswff,npwe_asked,comm
 character(len=*),intent(in) :: fname
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: iomega,fform,rdwr,my_rank,master,unt,unclassified
 real(dp) :: re, im, tol
 character(len=500) :: msg

! *********************************************************************

 DBG_ENTER("COLL")

 !@Epsilonm1_results
 my_rank = xcomm_rank(comm)
 master=0

 ! === Open file ===
 !if (my_rank==master.or.localrdwf==1) then
 write(msg,'(3a)')' init_Er_from_file : testing file ',TRIM(fname),ch10
 call wrtout(std_out,msg,'COLL')
 if (open_file(fname,msg,newunit=unt,form="unformatted",status="old",action="read") /= 0) then
   MSG_ERROR(msg)
 end if
 !end if

 rdwr=5
 call scr_hdr_io(fform,rdwr,unt,comm,master,accesswff,Er%Hscr)
 !if (my_rank==master.or.localrdwf==1) close(unt)
 close(unt)

 ! === Master echoes the header ===
 if (my_rank==master) then
   rdwr=4
   call scr_hdr_io(fform,rdwr,std_out,comm,master,accesswff,Er%Hscr)
 end if

 ! === Generic Info ===
 Er%ID         =0       ! Not yet initialized as epsm1 is calculated in mkdump_Er.F90
 Er%fname      =fname
 Er%fform      =fform
 if (fform==2002) then
   Er%npoles     =Er%Hscr%npoles
   Er%ncoeff     =Er%Hscr%ncoeff
 else
   Er%npoles     =0
   Er%ncoeff     =0
 end if

 Er%Tordering=Er%Hscr%Tordering

!TODO these quantitities should be checked and initiliazed in mkdump_Er
!BEGIN HARCODED
 Er%nI       = 1
 Er%nJ       = 1
 Er%ikxc     = 0
 Er%test_type=-1

 Er%Hscr%headform=HSCR_LATEST_HEADFORM   ! XG20090912
!END HARDCODED

 Er%nqibz=Er%Hscr%nqibz
 Er%mqmem=mqmem ; if (mqmem/=0) Er%mqmem=Er%nqibz
 ABI_MALLOC(Er%qibz,(3,Er%nqibz))
 Er%qibz(:,:)=Er%Hscr%qibz(:,:)

 Er%nqlwl=Er%Hscr%nqlwl
 ABI_MALLOC(Er%qlwl,(3,Er%nqlwl))
 Er%qlwl(:,:)=Er%Hscr%qlwl(:,:)

 Er%nomega=Er%Hscr%nomega
 ABI_MALLOC(Er%omega,(Er%nomega))
 Er%omega(:)=Er%Hscr%omega(:)

 ! Count number of real, imaginary, and complex frequencies.
 Er%nomega_r = 1
 Er%nomega_i = 0
 Er%nomega_c = 0
 if (Er%nomega == 2) then
   Er%nomega_i = 1
 else
   unclassified = 0
   tol = tol6*Ha_eV
   do iomega = 2, Er%nomega
     re =  REAL(Er%omega(iomega))
     im = AIMAG(Er%omega(iomega))
     if ((re > tol) .AND. (im < tol)) then
       Er%nomega_r = iomega ! Real freqs are packed in the first locations.
     else if ((re < tol) .AND. (im > tol)) then
       Er%nomega_i = Er%nomega_i + 1
     else if ((re > tol) .AND. (im > tol)) then
       Er%nomega_c = Er%nomega_c + 1
     else
       unclassified = unclassified + 1
     end if
   end do
   if (unclassified > 0) then
     write(msg,'(3a,i6)')&
&      'Some complex frequencies are too small to qualify as real or imaginary.',ch10,&
&      'Number of unidentified frequencies = ', unclassified
     MSG_WARNING(msg)
   end if
 end if

 ! === Get G-vectors ===
 Er%npwe=Er%Hscr%npwe
 if (npwe_asked>0) then
   if (npwe_asked>Er%Hscr%npwe) then
     write(msg,'(a,i8,2a,i8)')&
&     'Number of G-vectors saved on file is less than the value required = ',npwe_asked,ch10,&
&     'Calculation will proceed with Max available npwe = ',Er%Hscr%npwe
     MSG_WARNING(msg)
   else  ! Redefine the no. of G"s for W.
     Er%npwe=npwe_asked
   end if
 end if

 ! pointer to Er%Hscr%gvec ?
 ABI_MALLOC(Er%gvec,(3,Er%npwe))
 Er%gvec=Er%Hscr%gvec(:,1:Er%npwe)

 DBG_EXIT("COLL")

end subroutine init_Er_from_file
!!***

!----------------------------------------------------------------------

!!****f* m_screening/mkdump_Er
!! NAME
!!  mkdump_Er
!!
!! FUNCTION
!!  Dump the content of an Epsilonm1_results data type on file.
!!
!! INPUTS
!!  id_required=Identifier of the matrix to be calculated
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  ngfft(18)=Info on the FFT mesh.
!!  nfftot=Total number of point on the FFT mesh.
!!  gvec(3,npwe)=Reduced coordinates of plane waves for the response functions
!!  npwe=Number of plane waves.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!      mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkdump_Er(Er,Vcp,npwe,gvec,nkxc,kxcg,id_required,approx_type,&
&                    ikxc_required,option_test,fname_dump,accesswff,&
&                    nfftot,ngfft,comm,fxc_ADA,reconstruct_scr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkdump_Er'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: id_required,approx_type,option_test,ikxc_required,nkxc
 integer,intent(in) :: accesswff,nfftot,npwe,comm
 integer,intent(in),optional :: reconstruct_scr
 type(Epsilonm1_results),intent(inout) :: Er
 type(vcoul_t),intent(in) :: Vcp
 character(len=*),intent(in) :: fname_dump
!arrays
 integer,intent(in) :: ngfft(18),gvec(3,npwe)
 complex(gwpc),intent(in) :: kxcg(nfftot,nkxc)
 complex(gwpc),intent(in), optional :: fxc_ADA(Er%npwe*Er%nI,Er%npwe*Er%nJ,Er%nqibz)

!Local variables-------------------------------
!scalars
 integer :: dim_wing,iqibz,is_qeq0,mqmem_,npwe_asked
 integer :: unt_dump,fform,rdwr,ig1,ig2
 integer :: master,my_rank,comm_self
 real(dp) :: ucvol
 character(len=500) :: msg
 type(ScrHdr_type) :: Hscr_cp
 type(Spectra_type) :: Spectra
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 complex(gwpc),allocatable :: epsm1(:,:,:)
 complex(dpc),allocatable :: dummy_lwing(:,:,:),dummy_uwing(:,:,:),dummy_head(:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(id_required==4,'Value of id_required not coded')
 ABI_CHECK(npwe==Er%npwe,"mismatch in npwe")

 my_rank = xcomm_rank(comm)
 master=0
 comm_self = xmpi_self

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 !! if (Er%ID/=0) then
 !!   call reset_Epsilonm1(Er)
 !! end if
 Er%ID=id_required

 write(std_out,*) 'Er%ID:',Er%ID
 write(std_out,*) 'Er%Hscr%ID:',Er%Hscr%ID

 if (Er%ID==Er%Hscr%ID) then
   ! === The two-point function we are asking for is already stored on file ===
   ! * According to mqmem either read and store the entire matrix in memory or do nothing.
   !
   if (Er%mqmem>0) then ! In-core solution.
     ABI_MALLOC(Er%lwing,(Er%npwe,Er%nomega,Er%nqlwl))
     ABI_MALLOC(Er%uwing,(Er%npwe,Er%nomega,Er%nqlwl))
     if (Er%nqlwl>0) then
       Er%lwing(:,:,:)=Er%Hscr%lwing(1:Er%npwe,1:Er%nomega,1:Er%nqlwl)
       Er%uwing(:,:,:)=Er%Hscr%uwing(1:Er%npwe,1:Er%nomega,1:Er%nqlwl)
     end if

     if (Er%fform/=2002) then ! Check for pole-fit scr
       write(msg,'(a,f12.1,a)')' Memory needed for Er%epsm1 = ',2*gwpc*Er%npwe**2*Er%nomega*Er%nqibz*b2Mb,' [Mb]'
       call wrtout(std_out,msg,'PERS')
       ABI_MALLOC(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,Er%nqibz))
       ABI_CHECK_ALLOC("Out-of-memory in Er%epsm1 (in-core)")

#ifdef HAVE_MPI_IO
       call wrtout(std_out,ABI_FUNC//"read_screening with MPI_IO","COLL")
       !call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,IO_MODE_MPI,comm)
       call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,accesswff,comm)
#else
       call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,accesswff,comm)
#endif

     else
       !write(std_out,*) Er%fform
       !write(std_out,*) Er%npwe,Er%npwe,Er%ncoeff,Er%nqibz
       ABI_MALLOC(Er%epsm1_pole,(Er%npwe,Er%npwe,Er%ncoeff,Er%nqibz))
       ABI_CHECK_ALLOC('Out-of-memory in Er%epsm1_pole (in-core)')

       call read_pole_screening(Er%fname,Er%npwe,Er%nqibz,Er%ncoeff,Er%epsm1_pole,accesswff,comm)

       if (present(reconstruct_scr)) then
         if (reconstruct_scr == 1) then
           ! We are supposed to reconstruct the screening file from the
           ! pole-fit for each omega as if it was all read from file.
           ABI_MALLOC(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,Er%nqibz))
           ABI_CHECK_ALLOC('Out-of-memory in Er%epsm1 (in-core, reconstruct screening)')
           do iqibz=1,Er%nqibz
             do ig2=1,Er%npwe
               do ig1=1,Er%npwe
                 call re_and_im_screening_with_phase(Er%omega,Er%epsm1(ig1,ig2,:,iqibz),Er%nomega, &
&                 Er%epsm1_pole(ig1,ig2,:,iqibz),Er%ncoeff)
                 if (ig1==ig2) Er%epsm1(ig1,ig2,:,iqibz) = &
&                 Er%epsm1(ig1,ig2,:,iqibz) + CMPLX(1.0_gwp,0.0_gwp)
               end do
             end do
           end do
           if (allocated(Er%epsm1_pole))  then
             ABI_FREE(Er%epsm1_pole)
           end if
           Er%fform = 1002
           Er%npoles=0
           Er%ncoeff=0
           msg = "--- NOTE: Reconstructed full screening from pole-fit."
           MSG_COMMENT(msg)
         end if
       end if
     end if
     !
   else  ! Out-of-core solution ===
     msg = " mqmem==0 => allocating a single q-slice of (W|chi0) (slower but less memory)."
     MSG_COMMENT(msg)
     continue
   end if

   RETURN

 else
   ! === The matrix stored on file do not correspond to the quantity required ===
   ! * Presently only the transformation chi0 => e^-1 is coded
   ! * According to Er%mqmem either calculate e^-1 dumping the result to a file
   !   for a subsequent use or calculate e^-1 keeping everything in memory.
   if (Er%mqmem==0) then
     ! === Open file and write the header for the SCR file ===
     ! * For the moment only master works.
     if (my_rank==master) then
       !
       write(msg,'(3a)')ch10,&
&       'mkdump_Er : calculating and writing epsilon^-1 matrix on file ',TRIM(fname_dump)
       call wrtout(std_out,msg,'COLL')

       if (open_file(fname_dump,msg,newunit=unt_dump,form="unformatted",status="unknown",action="write") /= 0) then
         MSG_ERROR(msg)
       end if

       ! === Update the entries in the header that have been modified ===
       ! TODO, write function to return title, just for info
       call copy_scrhdr(Er%Hscr,Hscr_cp)
       Hscr_cp%ID        = id_required
       Hscr_cp%ikxc      = ikxc_required
       Hscr_cp%test_type = option_test
       Hscr_cp%title(1)  = 'SCR file: epsilon^-1'
       Hscr_cp%title(2)  = 'TESTPARTICLE'

       ! Treat the case in which a smaller matrix is used.
       Hscr_cp%npwe      = Er%npwe
       if (Hscr_cp%nqlwl>0) then ! Reallocate wing with correct dimension keeping previous values.
         ABI_MALLOC(dummy_lwing,(Er%npwe*Er%nI,Er%nomega,Hscr_cp%nqlwl))
         dummy_lwing(:,:,:) = Hscr_cp%lwing(1:Er%npwe*Er%nI,:,:)
         ABI_FREE(Hscr_cp%lwing)
         ABI_MALLOC(Hscr_cp%lwing,(Hscr_cp%npwe,Hscr_cp%nomega,Hscr_cp%nqlwl))
         Hscr_cp%lwing = dummy_lwing
         dummy_lwing(:,:,:) = Hscr_cp%uwing(1:Er%npwe*Er%nI,:,:)
         ABI_FREE(Hscr_cp%uwing)
         ABI_MALLOC(Hscr_cp%uwing,(Hscr_cp%npwe,Hscr_cp%nomega,Hscr_cp%nqlwl))
         Hscr_cp%uwing = dummy_lwing
         ABI_FREE(dummy_lwing)
       end if

       rdwr=2; fform=Hscr_cp%fform
       call scr_hdr_io(fform,rdwr,unt_dump,comm_self,master,accesswff,Hscr_cp)
       call free_scrhdr(Hscr_cp)

       ABI_MALLOC(epsm1,(Er%npwe,Er%npwe,Er%nomega))
       ABI_CHECK_ALLOC("out of memory in epsm1")

       do iqibz=1,Er%nqibz
         is_qeq0=0
         if (normv(Er%qibz(:,iqibz),gmet,'G')<GW_TOLQ0) is_qeq0=1

#ifdef HAVE_MPI_IO
         ! FIXME there's a problem with SUSC files and MPI-IO
         !call read_screening(Er%fname,Er%npwe,1,Er%nomega,epsm1,IO_MODE_MPI,comm_self,iqiA=iqibz)
         call read_screening(Er%fname,Er%npwe,1,Er%nomega,epsm1,accesswff,comm_self,iqiA=iqibz)
#else
         call read_screening(Er%fname,Er%npwe,1,Er%nomega,epsm1,accesswff,comm_self,iqiA=iqibz)
#endif

         dim_wing=0; if (is_qeq0==1) dim_wing=3
         ABI_MALLOC(dummy_lwing,(Er%npwe*Er%nI,Er%nomega,dim_wing))
         ABI_MALLOC(dummy_uwing,(Er%npwe*Er%nJ,Er%nomega,dim_wing))
         ABI_MALLOC(dummy_head,(dim_wing,dim_wing,Er%nomega))

         if (approx_type<2) then
           MSG_WARNING('Entering out-of core RPA or Kxc branch')
           call make_epsm1_driver(iqibz,dim_wing,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                    approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                    dummy_lwing,dummy_uwing,epsm1,Spectra,comm_self)
         else
           MSG_WARNING('Entering out-of core fxc_ADA branch')
           call make_epsm1_driver(iqibz,dim_wing,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                    approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                    dummy_lwing,dummy_uwing,epsm1,Spectra,comm_self,fxc_ADA(:,:,iqibz))
         end if

         ABI_FREE(dummy_head)
         ABI_FREE(dummy_uwing)
         ABI_FREE(dummy_lwing)

         if (is_qeq0==1) then
          call repr_dielconst(Spectra,msg)
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out,msg,'COLL')
         end if

         call destroy_spectra(Spectra)

         call write_screening(unt_dump,accesswff,Er%npwe,Er%nomega,epsm1)
       end do

       close(unt_dump)
       ABI_FREE(epsm1)
     end if !master
     !
     ! A synchronization is required here, else the other procs start to read the
     ! file _SCR before it is written by the master
     call xmpi_barrier(comm)

     ! Now Er% "belongs" to the file "fname_dump", thus
     ! each proc has to destroy and re-initialize the object.
     call destroy_Epsilonm1_results(Er)

     mqmem_=Er%mqmem; npwe_asked=Er%npwe
     call init_Er_from_file(Er,fname_dump,mqmem_,npwe_asked,accesswff,comm)

     !Now Er% has been reinitialized and ready-to-use.
     Er%ID=id_required
     call print_epsilonm1_results(Er)
   else
     ! ========================
     ! === In-core solution ===
     ! ========================
     ABI_MALLOC(Er%lwing,(Er%npwe,Er%nomega,Er%nqlwl))
     ABI_MALLOC(Er%uwing,(Er%npwe,Er%nomega,Er%nqlwl))
     ABI_MALLOC(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,Er%nqibz))
     ABI_CHECK_ALLOC('Out-of-memory in Er%epsm1 (in-core)')

#ifdef HAVE_MPI_IO
     ! FIXME there's a problem with SUSC files and MPI-IO
     !call wrtout(std_out,ABI_FUNC//"read_screening with MPI_IO","COLL")
     !call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,IO_MODE_MPI,comm)
     call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,accesswff,comm)
#else
     call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,accesswff,comm)
#endif

     do iqibz=1,Er%nqibz
       is_qeq0=0
       if (normv(Er%qibz(:,iqibz),gmet,'G')<GW_TOLQ0) is_qeq0=1

       dim_wing=0; if (is_qeq0==1) dim_wing=3 ! FIXME
       ABI_MALLOC(dummy_lwing,(Er%npwe*Er%nI,Er%nomega,dim_wing))
       ABI_MALLOC(dummy_uwing,(Er%npwe*Er%nJ,Er%nomega,dim_wing))
       ABI_MALLOC(dummy_head,(dim_wing,dim_wing,Er%nomega))

       if (approx_type<2) then
         MSG_WARNING('Entering in-core RPA and Kxc branch')
         call make_epsm1_driver(iqibz,dim_wing,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                  approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                  dummy_lwing,dummy_uwing,Er%epsm1(:,:,:,iqibz),Spectra,comm)
       else
         MSG_WARNING('Entering in-core fxc_ADA branch')
         call make_epsm1_driver(iqibz,dim_wing,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                  approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                  dummy_lwing,dummy_uwing,Er%epsm1(:,:,:,iqibz),Spectra,comm,fxc_ADA=fxc_ADA(:,:,iqibz))
       end if

       ABI_FREE(dummy_lwing)
       ABI_FREE(dummy_uwing)
       ABI_FREE(dummy_head)

       if (is_qeq0==1) then
         call repr_dielconst(Spectra,msg)
         call wrtout(std_out,msg,'COLL')
         call wrtout(ab_out,msg,'COLL')
       end if

       call destroy_spectra(Spectra)
     end do

     Er%ID=id_required
     call print_epsilonm1_results(Er)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine mkdump_Er
!!***

!----------------------------------------------------------------------

!!****f* m_screening/get_epsm1
!! NAME
!!  get_epsm1
!!
!! FUNCTION
!!  Work in progress but the main is idea is as follows:
!!
!!  Return the symmetrized inverse dielectric matrix.
!!  This method implements both in-core and the out-of-core solution
!!  In the later, epsilon^-1 or chi0 are read from file.
!!  It is possible to specify options to retrieve (RPA |TDDDT, [TESTCHARGE|TESTPARTICLE]).
!!  All dimensions are already initialized in the Er% object, this method
!!  should act as a wrapper around rdscr and make_epsm1_driver. A better
!!  implementation will be done in the following once the coding of file handlers is completed.
!!
!! INPUTS
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  iqibzA[optional]=Index of the q-point to be read from file (only for out-of-memory solutions)
!!  accesswff=option definig the file format.
!!  option_test
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Er%epsm1
!!
!! TODO
!!  Remove this routine. Now everything should be done with mkdump_Er
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_epsm1(Er,Vcp,approx_type,option_test,accesswff,comm,iqibzA)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_epsm1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,option_test,approx_type,comm
 integer,optional,intent(in) :: iqibzA
 type(vcoul_t),intent(in) :: Vcp
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: my_approx_type,my_option_test,ng
 !character(len=500) :: msg

! *********************************************************************

 DBG_ENTER("COLL")

 my_approx_type = approx_type
 my_option_test = option_test

 ! Vcp not yet used.
 ng = Vcp%ng

 select case (Er%mqmem)

 case (0) !  Out-of-core solution
   if (allocated(Er%lwing))  then
     ABI_FREE(Er%lwing)
   end if
   if (allocated(Er%uwing))  then
     ABI_FREE(Er%uwing)
   end if
   if (allocated(Er%epsm1))  then
     ABI_FREE(Er%epsm1)
   end if
   ABI_MALLOC(Er%lwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_MALLOC(Er%uwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_MALLOC(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,1))
   ABI_CHECK_ALLOC('Out-of-memory in Er%epsm1 (out-of-core)')

#ifdef HAVE_MPI_IO
   ! FIXME there's a problem with SUSC files and MPI-IO
   !call wrtout(std_out,ABI_FUNC//": calling read_screening with MPI_IO","COLL")
   !call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,IO_MODE_MPI,comm,iqiA=iqibzA)
   call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,accesswff,comm,iqiA=iqibzA)
#else
   call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,accesswff,comm,iqiA=iqibzA)
#endif

   if (Er%ID==4) then
     ! === If q-slice of epsilon^-1 has been read then return ===
     !call print_epsilonm1_results(Er)
     RETURN
   else
     MSG_ERROR('Wrong Er%ID')
   end if

 case default
   ! ========================
   ! === In-core solution ===
   ! ========================
   MSG_ERROR("you should not be here")
 end select

 DBG_EXIT("COLL")

end subroutine get_epsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/get_pole_epsm1
!! NAME
!!  get_epsm1
!!
!! FUNCTION
!!  Work in progress but the main is idea is as follows:
!!
!!  Return the symmetrized inverse dielectric matrix.
!!  This method implements both in-core and the out-of-core solution
!!  In the later, epsilon^-1 or chi0 are read from file.
!!  It is possible to specify options to retrieve (RPA |TDDDT, [TESTCHARGE|TESTPARTICLE]).
!!  All dimensions are already initialized in the Er% object, this method
!!  should act as a wrapper around rdscr and make_epsm1_driver. A better
!!  implementation will be done in the following once the coding of file handlers is completed.
!!
!! INPUTS
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  iqibzA[optional]=Index of the q-point to be read from file (only for out-of-memory solutions)
!!  accesswff=option definig the file format.
!!  option_test
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Er%epsm1
!!
!! TODO
!!  Remove this routine. Now everything should be done with mkdump_Er
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! SOURCE

subroutine get_pole_epsm1(Er,Vcp,approx_type,option_test,accesswff,comm,iqibzA,reconstruct_scr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_pole_epsm1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,option_test,approx_type,comm
 integer,optional,intent(in) :: iqibzA,reconstruct_scr
 type(vcoul_t),intent(in) :: Vcp
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: my_approx_type,my_option_test,ng,ig1,ig2
 character(len=500) :: msg

! *********************************************************************

 DBG_ENTER("COLL")

 my_approx_type = approx_type
 my_option_test = option_test

 ! Vcp not yet used.
 ng = Vcp%ng

 select case (Er%mqmem)

 case (0) !  Out-of-core solution
   if (allocated(Er%lwing))  then
     ABI_FREE(Er%lwing)
   end if
   if (allocated(Er%uwing))  then
     ABI_FREE(Er%uwing)
   end if
   if (allocated(Er%epsm1))  then
     ABI_FREE(Er%epsm1)
   end if
   if (allocated(Er%epsm1_pole))  then
     ABI_FREE(Er%epsm1_pole)
   end if
   ABI_MALLOC(Er%lwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_MALLOC(Er%uwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_MALLOC(Er%epsm1_pole,(Er%npwe,Er%npwe,Er%ncoeff,1))
   ABI_CHECK_ALLOC('Out-of-memory in Er%epsm1_pole (out-of-core)')

   call read_pole_screening(Er%fname,Er%npwe,Er%nqibz,Er%ncoeff,Er%epsm1_pole,accesswff,comm,iqiA=iqibzA)

   if (present(reconstruct_scr)) then
     if (reconstruct_scr==1) then
       ABI_MALLOC(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,1))
       ABI_CHECK_ALLOC('Out-of-memory in Er%epsm1 (out-of-core, reconstruct screening)')
       do ig2=1,Er%npwe
         do ig1=1,Er%npwe
           call re_and_im_screening_with_phase(Er%omega,Er%epsm1(ig1,ig2,:,1),Er%nomega, &
&           Er%epsm1_pole(ig1,ig2,:,1),Er%ncoeff)
           if (ig1==ig2) Er%epsm1(ig1,ig2,:,1) = &
&            Er%epsm1(ig1,ig2,:,1) + CMPLX(1.0_gwp,0.0_gwp)
         end do
       end do
       if (allocated(Er%epsm1_pole))  then
         ABI_FREE(Er%epsm1_pole)
       end if
       msg = "--- NOTE: Reconstructed full screening from pole-fit. (out-of-core)"
       MSG_COMMENT(msg)
     end if
   end if

   if (Er%ID==4) then
     ! === If q-slice of epsilon^-1 has been read then return ===
     !call print_epsilonm1_results(Er)
     RETURN
   else
     MSG_ERROR('Wrong Er%ID')
   end if

 case default
   ! ========================
   ! === In-core solution ===
   ! ========================
   MSG_ERROR("you should not be here")
 end select

 DBG_EXIT("COLL")

end subroutine get_pole_epsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/decompose_epsm1
!! NAME
!! decompose_epsm1
!!
!! FUNCTION
!! Decompose the complex symmetrized dielectric
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine decompose_epsm1(Er,iqibz,eigenvalues)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'decompose_epsm1'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz
 type(Epsilonm1_results),intent(in) :: Er
!arrays
 complex(dpc),intent(out) :: eigenvalues(Er%npwe,Er%nomega)

!Local variables-------------------------------
!scalars
 integer :: info,lwork,iomega,negw,ig1,ig2,idx,sdim
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: ww(:),rwork(:)
 complex(dpc),allocatable :: work(:),Adpp(:),eigvec(:,:),Afull(:,:),vs(:,:),wwc(:)
 logical,allocatable :: bwork(:)
 logical :: sortcplx !BUG in abilint

! *********************************************************************

 ABI_CHECK(Er%mqmem/=0,'mqmem==0 not implemented')

 do iomega=1,Er%nomega

   if (ABS(REAL(Er%omega(iomega)))>0.00001) then ! Eigenvalues for a generic complex matrix ===
     !if (.TRUE.) then
     lwork=4*2*Er%npwe
     ABI_MALLOC(wwc,(Er%npwe))
     ABI_MALLOC(work,(lwork))
     ABI_MALLOC(rwork,(Er%npwe))
     ABI_MALLOC(bwork,(Er%npwe))
     ABI_MALLOC(vs,(Er%npwe,Er%npwe))
     ABI_MALLOC(Afull,(Er%npwe,Er%npwe))

     Afull=Er%epsm1(:,:,iomega,iqibz)

     !for the moment no sort, maybe here I should sort using the real part?
     call ZGEES('V','N',sortcplx,Er%npwe,Afull,Er%npwe,sdim,wwc,vs,Er%npwe,work,lwork,rwork,bwork,info)
     if (info/=0) then
       write(msg,'(2a,i10)')' decompose_epsm1 : Error in ZGEES, diagonalizing complex matrix, info = ',info
       call wrtout(std_out,msg,'COLL')
     end if

     eigenvalues(:,iomega)=wwc(:)

     ABI_FREE(wwc)
     ABI_FREE(work)
     ABI_FREE(rwork)
     ABI_FREE(bwork)
     ABI_FREE(vs)
     ABI_FREE(Afull)

   else ! Hermitian version.

     lwork=2*Er%npwe-1
     ABI_MALLOC(ww,(Er%npwe))
     ABI_MALLOC(work,(lwork))
     ABI_MALLOC(rwork,(3*Er%npwe-2))
     ABI_MALLOC(eigvec,(Er%npwe,Er%npwe))
     ABI_MALLOC(Adpp,(Er%npwe*(Er%npwe+1)/2))
     ABI_CHECK_ALLOC('out of memory in Adpp')

     idx=0 ! Pack the matrix
     do ig2=1,Er%npwe
       do ig1=1,ig2
         idx=idx+1
         Adpp(idx)=Er%epsm1(ig1,ig2,iomega,iqibz)
       end do
     end do

     ! For the moment we require also the eigenvectors.
     call ZHPEV('V','U',Er%npwe,Adpp,ww,eigvec,Er%npwe,work,rwork,info)
     if (info/=0) then
       write(msg,'(a,i5)')' decompose_epsm1 : Error in ZHPEV, diagonalizing matrix, info = ',info
       MSG_ERROR(msg)
     end if

     negw=(COUNT((REAL(ww)<tol6)))
     if (negw/=0) then
       write(msg,'(a,i5,a,i3,a,f8.4)')&
&        'Found negative eigenvalues. No. ',negw,' at iqibz= ',iqibz,' minval= ',MINVAL(REAL(ww))
       MSG_WARNING(msg)
     end if

     eigenvalues(:,iomega)=ww(:)

     ABI_FREE(ww)
     ABI_FREE(work)
     ABI_FREE(rwork)
     ABI_FREE(eigvec)
     ABI_FREE(Adpp)
   end if
 end do !iomega

! contains
! function sortcplx(carg) result(res)
!  implicit none
!  complex(dpc),intent(in) :: carg
!  logical :: res
!  res=.TRUE.
! end function sortcplx

end subroutine decompose_epsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/make_epsm1_driver
!! NAME
!! make_epsm1_driver
!!
!! FUNCTION
!!  Driver routine to calculate the inverse symmetrical dielectric matrix starting
!!  from the irreducible polarizability. The routine considers a single q-point, and
!!  performs the following tasks:
!!
!!  1) Calculate $\tilde\epsilon^{-1}$ using different approximations:
!!      * RPA
!!      * ALDA within TDDFT
!!
!!  2) Use a special treatment of non-Analytic behavior of heads and wings in reciprocal space
!!     calculating these quantities for different small q-directions specified by the user
!!     (Not yet operative)
!!
!!  3) Output the electron energy loss function and the macroscopic dielectric function with and
!!     without local field effects (only if non-zero real frequencies are available)
!!
!! INPUTS
!!  iqibz=index of the q-point in the array Vcp%qibz where epsilon^-1 has to be calculated
!!  npwe=Number of G-vectors in chi0.
!!  nI,nJ=Number of rows/columns in chi0_ij (1,1 in collinear case)
!!  nomega=Number of frequencies.
!!  dim_wing=Dimension of the wings (0 or 3 if q-->0)
!!  approx_type=Integer flag defining the type of approximation
!!   == 0 for RPA   ==
!!   == 1 for TDDFT ==
!!  option_test=Only for TDDFT:
!!   == 0 for TESTPARTICLE ==
!!   == 1 for TESTELECTRON ==
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!   %nqibz=Number of q-points.
!!   %qibz(3,nqibz)=q-points in the IBZ.
!!  nkxc=Integer defining the dimension of the kernel in reciprocal space
!!  kxcg(nfftot,nkxc)=TDDFT kernel in reciprocal space on the FFT mesh. Needed only if approx_type==1
!!  comm=MPI communicator.
!!  ngfft(18)=Info on the FFT mesh.
!!  nfftot=Total number of points in the FFT mesh.
!!  chi0_lwing(npwe*nI,nomega,dim_wing)=Lower wings of chi0 (only for q-->0)
!!  chi0_uwing(npwe*nJ,nomega,dim_wing)=Upper wings of chi0 (only for q-->0)
!!  chi0_head(dim_wing,dim_wing,nomega)=Head of of chi0 (only for q-->0)
!!
!! OUTPUT
!!  Different files are written according to the type of calculation
!!  See also side effects
!!
!! SIDE EFFECTS
!!  chi0(npwe*nI,npwe*nJ,nomega): in input the irreducible polarizability, in output
!!   the symmetrized inverse dielectric matrix.
!!
!! PARENTS
!!      m_screen,m_screening,screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_epsm1_driver(iqibz,dim_wing,npwe,nI,nJ,nomega,omega,&
& approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,chi0_head,&
& chi0_lwing,chi0_uwing,chi0,Spectra,comm,fxc_ADA)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_epsm1_driver'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nI,nJ,npwe,nomega,dim_wing,approx_type,option_test,nkxc,nfftot,comm
 type(vcoul_t),intent(in) :: Vcp
 type(Spectra_type),intent(out) :: Spectra
!arrays
 integer,intent(in) :: ngfft(18),gvec(3,npwe)
 complex(gwpc),intent(in) :: kxcg(nfftot,nkxc)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(inout) :: chi0_lwing(:,:,:)   !(npwe*nI,nomega,dim_wing)
 complex(dpc),intent(inout) :: chi0_uwing(:,:,:)   !(npwe*nJ,nomega,dim_wing)
 complex(dpc),intent(inout) :: chi0_head(:,:,:)   !(dim_wing,dim_wing,nomega)
 complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ,nomega)
 complex(gwpc),intent(in), optional :: fxc_ADA(npwe*nI,npwe*nJ)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ig1,ig2,io,ierr,irank,master,my_nqlwl !iqlwl
 integer :: nor,my_rank,nprocs,comm_self,g1mg2_idx
 real(dp) :: ucvol
 logical :: is_qeq0,use_MPI
 character(len=500) :: msg
!arrays
 integer :: omega_distrb(nomega)
 integer,allocatable :: istart(:),istop(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: eelf(:,:),tmp_eelf(:)
 !complex(gwpc),allocatable :: chitmp(:,:)
 !complex(gwpc),pointer :: vc_sqrt(:)
 complex(dpc),allocatable :: epsm_lf(:,:),epsm_nlf(:,:),tmp_lf(:),tmp_nlf(:) 
 complex(dpc),allocatable :: buffer_lwing(:,:),buffer_uwing(:,:)
 complex(gwpc),allocatable :: kxcg_mat(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 ABI_UNUSED(chi0_head(1,1,1))

 if (nI/=1.or.nJ/=1) then
   MSG_ERROR("nI or nJ=/1 not yet implemented")
 end if

 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master  = 0

 ! MG TODO We use comm_self for the inversion as the single precision version is not yet available
 comm_self = xmpi_self

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 is_qeq0 = (normv(Vcp%qibz(:,iqibz),gmet,'G')<GW_TOLQ0)

 omega_distrb = my_rank
 use_MPI = .FALSE.
 use_MPI = (nprocs>=nomega)  ! Parallelism is not used

 if (use_MPI) then
   ! * Initialize distribution table for frequencies.
   ABI_MALLOC(istart,(nprocs))
   ABI_MALLOC(istop,(nprocs))
   call xmpi_split_work2_i4b(nomega,nprocs,istart,istop,msg,ierr)
   omega_distrb(:)=xmpi_undefined_rank
   do irank=0,nprocs-1
     i1 = istart(irank+1)
     i2 = istop (irank+1)
     if (i1<=i2) omega_distrb(i1:i2) = irank
   end do
   ABI_FREE(istart)
   ABI_FREE(istop)
 end if
 !
 ! * Initialize container for spectral results
 do nor=1,nomega
   if (ABS(AIMAG(omega(nor)))>1.e-3) EXIT
 end do
 nor=nor-1; if (nor==0) nor = 1 ! only imag !?

 if (dim_wing==3) then
   call wrtout(std_out,' Analyzing long wavelength limit for several q','COLL')
   call init_spectra(Spectra,nor,REAL(omega(1:nor)),Vcp%nqlwl,Vcp%qlwl)
   my_nqlwl = 1
   !my_nqlwl = dim_wing ! TODO
   !ABI_CHECK(dim_wing==SIZE(Vcp%vcqlwl_sqrt,DIM=2),"WRONG DIMS")
 else
   call init_spectra(Spectra,nor,REAL(omega(1:nor)),1,Vcp%qibz(:,iqibz))
   my_nqlwl = 1
 end if
 !
 ! NOTE: all processors have to perform this operation in order to have the
 !       epsm1 matrix when performing a sigma calculation starting with the file _SUS
 !
 ! Temporary arrays to store spectra.
 ABI_MALLOC(epsm_lf,(nomega,my_nqlwl))
 ABI_MALLOC(epsm_nlf,(nomega,my_nqlwl))
 ABI_MALLOC(eelf,(nomega,my_nqlwl))
 epsm_lf=czero; epsm_nlf=czero; eelf=zero

 ! Temporary arrays used to store output results.
 ABI_MALLOC(tmp_lf, (my_nqlwl))
 ABI_MALLOC(tmp_nlf, (my_nqlwl))
 ABI_MALLOC(tmp_eelf, (my_nqlwl))

 SELECT CASE (approx_type)

 CASE (0)
   ! * RPA: \tepsilon=1 - Vc^{1/2} chi0 Vc^{1/2}
   ! * vc_sqrt contains vc^{1/2}(q,G), complex-valued to allow for a possible cutoff.
   do io=1,nomega
     if (omega_distrb(io) == my_rank) then
       !write(std_out,*)"dim_wing",dim_wing
       call rpa_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),my_nqlwl,dim_wing,&
&        chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),&
&        tmp_lf,tmp_nlf,tmp_eelf,comm_self)

         ! Store results.
         epsm_lf(io,:) = tmp_lf
         epsm_nlf(io,:) = tmp_nlf
         eelf(io,:) = tmp_eelf
     end if
   end do ! nomega

 CASE (1) ! Vertex correction from Adiabatic TDDFT. chi_{G1,G2} = [\delta -\chi0 (vc+kxc)]^{-1}_{G1,G3} \chi0_{G3,G2}

   ABI_CHECK(Vcp%nqlwl==1,"nqlwl/=1 not coded")
   ABI_CHECK(nkxc==1,"nkxc/=1 not coded")

   ! Make kxcg_mat(G1,G2) = kxcg(G1-G2) from kxcg defined on the FFT mesh.
   ABI_MALLOC(kxcg_mat,(npwe,npwe))
   ABI_CHECK_ALLOC("out-of-memory kxcg_mat")

   ierr=0
   do ig2=1,npwe
     do ig1=1,npwe
       g1mg2_idx = g2ifft(gvec(:,ig1)-gvec(:,ig2),ngfft)
       if (g1mg2_idx>0) then
         kxcg_mat(ig1,ig2) = kxcg(g1mg2_idx,1)
       else
         ierr=ierr+1
         kxcg_mat(ig1,ig2) = czero
       end if
     end do
   end do

   if (ierr/=0) then
     write(msg,'(a,i4,3a)')&
&     'Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
&     'Enlarge the FFT mesh to get rid of this problem. '
     MSG_WARNING(msg)
   end if

   !FIXME "recheck TDDFT code and parallel"
   ABI_CHECK(nkxc==1,"nkxc/=1 not coded")
   do io=1,nomega
     if (omega_distrb(io) == my_rank) then
       call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),kxcg_mat,option_test,my_nqlwl,dim_wing,omega(io),&
&        chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm)

       ! Store results.
       epsm_lf(io,:) = tmp_lf
       epsm_nlf(io,:) = tmp_nlf
       eelf(io,:) = tmp_eelf
     end if
   end do

   ABI_FREE(kxcg_mat)
   do io=1,nomega
     write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

 CASE(2) ! ADA nonlocal vertex correction contained in fxc_ADA
   MSG_WARNING('Entered fxc_ADA branch: EXPERIMENTAL!')
   ! Test that argument was passed
   if (.NOT.present(fxc_ADA)) then
     MSG_ERROR('make_epsm1_driver was not called with optional argument fxc_ADA')
   end if

   ABI_CHECK(Vcp%nqlwl==1,"nqlwl/=1 not coded")

   do io=1,nomega
     if (omega_distrb(io) == my_rank) then
       call atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0(:,:,io),fxc_ADA,option_test,my_nqlwl,dim_wing,omega(io),&
&        chi0_head(:,:,io),chi0_lwing(:,io,:),chi0_uwing(:,io,:),tmp_lf,tmp_nlf,tmp_eelf,comm)

       ! Store results.
       epsm_lf(io,:) = tmp_lf
       epsm_nlf(io,:) = tmp_nlf
       eelf(io,:) = tmp_eelf
     end if
   end do

   do io=1,nomega
     write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

 CASE DEFAULT
   write(msg,'(a,i0)')'Wrong value for approx_type= ',approx_type
   MSG_BUG(msg)
 END SELECT

 if (use_MPI) then 
   ! Collect results on each node.
   ABI_MALLOC(buffer_lwing, (size(chi0_lwing,dim=1), size(chi0_lwing, dim=3)))
   ABI_MALLOC(buffer_uwing, (size(chi0_uwing,dim=1), size(chi0_uwing, dim=3)))

   do io=1,nomega
     if (omega_distrb(io)/=my_rank) then
       ! Zero arrays.
       chi0(:,:,io) = czero_gw
       if(dim_wing>0) then
          chi0_lwing(:,io,:) = zero
          chi0_uwing(:,io,:) = zero
          chi0_head(:,:,io)  = czero
       endif
       epsm_lf(io,:) = czero
       epsm_nlf(io,:) = czero
       eelf(io,:) = zero
     end if

     call xmpi_sum(chi0(:,:,io), comm,ierr)

     if(dim_wing>0) then
        ! Build contiguous arrays
        buffer_lwing = chi0_lwing(:,io,:)
        buffer_uwing = chi0_uwing(:,io,:)
        call xmpi_sum(buffer_lwing,comm,ierr)
        call xmpi_sum(buffer_uwing,comm,ierr)
        chi0_lwing(:,io,:) = buffer_lwing
        chi0_uwing(:,io,:) = buffer_uwing
        if (size(chi0_head(:,:,io))/= zero) then
          call xmpi_sum(chi0_head(:,:,io),comm,ierr)
        end if
     end if

   end do ! iomega

   call xmpi_sum(epsm_lf, comm,ierr )
   call xmpi_sum(epsm_nlf,comm,ierr)
   call xmpi_sum(eelf,    comm,ierr)
   ABI_FREE(buffer_lwing)
   ABI_FREE(buffer_uwing) 
 end if
 !
 ! * Save results in Spectra%, mind the slicing.
 Spectra%emacro_nlf(:,:) = epsm_nlf(1:nor,:)
 Spectra%emacro_lf (:,:) = epsm_lf (1:nor,:)
 Spectra%eelf      (:,:) = eelf    (1:nor,:)

 ABI_FREE(epsm_lf)
 ABI_FREE(epsm_nlf)
 ABI_FREE(eelf)

 ABI_FREE(tmp_lf)
 ABI_FREE(tmp_nlf)
 ABI_FREE(tmp_eelf)

 DBG_EXIT("COLL")

end subroutine make_epsm1_driver
!!***

!----------------------------------------------------------------------

!!****f* m_screening/rpa_symepsm1
!! NAME
!! rpa_symepsm1
!!
!! FUNCTION
!!  Calculate RPA $\tilde\epsilon^{-1}$
!!
!!  Use a special treatment of non-Analytic behavior of heads and wings in reciprocal space
!!  calculating these quantities for different small q-directions specified by the user
!!  (Not yet operative)
!!
!! INPUTS
!!  iqibz=index of the q-point in the array Vcp%qibz where epsilon^-1 has to be calculated
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  npwe=Number of G-vectors in chi0.
!!  nI,nJ=Number of rows/columns in chi0_ij (1,1 in collinear case)
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
!!
!! PARENTS
!!      m_screening,tddft_bootstrap
!!
!! CHILDREN
!!
!! SOURCE

subroutine rpa_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0,my_nqlwl,dim_wing,chi0_head,chi0_lwing,chi0_uwing,epsm_lf,epsm_nlf,eelf,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rpa_symepsm1'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nI,nJ,npwe,dim_wing,my_nqlwl,comm
 type(vcoul_t),target,intent(in) :: Vcp
!arrays
 complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ)
 complex(dpc),intent(inout) :: chi0_lwing(:,:) !(npwe*nI,dim_wing)
 complex(dpc),intent(inout) :: chi0_uwing(:,:) !(npwe*nJ,dim_wing)
 complex(dpc),intent(inout) :: chi0_head(:,:) !(dim_wing,dim_wing)
 real(dp),intent(out) :: eelf(my_nqlwl)
 complex(dpc),intent(out) :: epsm_lf(my_nqlwl),epsm_nlf(my_nqlwl)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,master,iqlwl,my_rank,nprocs
 real(dp) :: ucvol
 logical :: is_qeq0
 !character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 complex(gwpc),pointer :: vc_sqrt(:)
 complex(gwpc),allocatable :: chi0_save(:,:)

! *************************************************************************

 ABI_UNUSED(chi0_head(1,1))

 if (nI/=1.or.nJ/=1) then
   MSG_ERROR("nI or nJ=/1 not yet implemented")
 end if

 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master  = 0

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 is_qeq0 = (normv(Vcp%qibz(:,iqibz),gmet,'G')<GW_TOLQ0)
 if (is_qeq0) then
   ABI_CHECK(iqibz==1,"q is 0 but iq_ibz /= 1")
 end if
 !
 if (my_nqlwl>1) then
   ABI_MALLOC(chi0_save,(npwe*nI,npwe*nJ))
   chi0_save = chi0
 end if
 !
 ! Symmetrized RPA epsilon = 1 - Vc^{1/2} chi0 Vc^{1/2}
 !   * vc_sqrt contains vc^{1/2}(q,G), complex-valued to allow for a possible cutoff.
 !
 ! * Loop over small q"s (if any) to treat the nonanalytical behavior.
 do iqlwl=my_nqlwl,1,-1
   !
   if (my_nqlwl>1) then
     chi0(:,:) = chi0_save           ! restore pristine polarizability
     chi0(:,1) = chi0_lwing(:,iqlwl) ! change the wings
     chi0(1,:) = chi0_uwing(:,iqlwl)
   end if

   if (iqibz==1) then
     vc_sqrt => Vcp%vcqlwl_sqrt(:,iqlwl)  ! Use Coulomb term for q-->0
   else
     vc_sqrt => Vcp%vc_sqrt(:,iqibz)
   end if

   do ig2=1,npwe*nJ
     do ig1=1,npwe*nI
       chi0(ig1,ig2)=-vc_sqrt(ig1)*chi0(ig1,ig2)*vc_sqrt(ig2)
     end do
     chi0(ig2,ig2)=one+chi0(ig2,ig2)
   end do

   epsm_nlf(iqlwl)=chi0(1,1) ! * chi0, now contains \tepsilon.

   if (.FALSE.) then
     call wrtout(std_out,' Symmetrical epsilon(G,G'') ','COLL')
     call print_arr(chi0)
   end if
   !
   ! === Invert tepsilon and calculate macroscopic dielectric constant ===
   ! * epsm_lf(w)=1/epsm1(G=0,Gp=0,w).
   ! * Since G=Gp=0 there is no difference btw symmetrical and not symmetrical.
   !
   call xginv(chi0,npwe,comm=comm)

   epsm_lf(iqlwl) = one/chi0(1,1)
   eelf(iqlwl) = -AIMAG(chi0(1,1))

   if (.FALSE.) then
     call wrtout(std_out," Symmetrical epsilon^-1(G,G'')",'COLL')
     call print_arr(chi0)
   end if
   !
   ! Save wings of e^-1 overwriting input values.
   if (dim_wing>0.and..FALSE.) then
     chi0_lwing(:,iqlwl) = chi0(:,1)
     chi0_uwing(:,iqlwl) = chi0(1,:)
   end if

 end do !iqlwl

 if (allocated(chi0_save))  then
   ABI_FREE(chi0_save)
 end if

end subroutine rpa_symepsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/atddft_symepsm1
!! NAME
!! atddft_symepsm1
!!
!! FUNCTION
!!  Calculate $\tilde\epsilon^{-1}$ using ALDA within TDDFT
!!
!!  2) Use a special treatment of non-Analytic behavior of heads and wings in reciprocal space
!!     calculating these quantities for different small q-directions specified by the user
!!     (Not yet operative)
!!
!!  Output the electron energy loss function and the macroscopic dielectric function with and
!!  without local field effects (only if non-zero real frequencies are available)
!!
!! INPUTS
!!  iqibz=index of the q-point in the array Vcp%qibz where epsilon^-1 has to be calculated
!!  npwe=Number of G-vectors in chi0.
!!  nI,nJ=Number of rows/columns in chi0_ij (1,1 in collinear case)
!!  dim_wing=Dimension of the wings (0 or 3 if q-->0)
!!  option_test=Only for TDDFT:
!!   == 0 for TESTPARTICLE ==
!!   == 1 for TESTELECTRON ==
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!   %nqibz=Number of q-points.
!!   %qibz(3,nqibz)=q-points in the IBZ.
!!  comm=MPI communicator.
!!  chi0_lwing(npwe*nI,dim_wing)=Lower wings of chi0 (only for q-->0)
!!  chi0_uwing(npwe*nJ,dim_wing)=Upper wings of chi0 (only for q-->0)
!!  chi0_head(dim_wing,dim_wing)=Head of of chi0 (only for q-->0)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chi0(npwe*nI,npwe*nJ): in input the irreducible polarizability, in output
!!   the symmetrized inverse dielectric matrix.
!!
!! PARENTS
!!      m_screening,tddft_bootstrap
!!
!! CHILDREN
!!
!! SOURCE

subroutine atddft_symepsm1(iqibz,Vcp,npwe,nI,nJ,chi0,kxcg_mat,option_test,my_nqlwl,dim_wing,omega,&
& chi0_head,chi0_lwing,chi0_uwing,epsm_lf,epsm_nlf,eelf,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'atddft_symepsm1'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nI,nJ,npwe,dim_wing,my_nqlwl
 integer,intent(in) :: option_test,comm
 type(vcoul_t),target,intent(in) :: Vcp
!arrays
 complex(gwpc),intent(in) :: kxcg_mat(npwe,npwe)
 complex(dpc),intent(in) :: omega
 complex(dpc),intent(inout) :: chi0_lwing(npwe*nI,dim_wing)
 complex(dpc),intent(inout) :: chi0_uwing(npwe*nJ,dim_wing)
 complex(dpc),intent(inout) :: chi0_head(dim_wing,dim_wing)
 complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ)
 real(dp),intent(out) :: eelf(my_nqlwl)
 complex(dpc),intent(out) :: epsm_lf(my_nqlwl),epsm_nlf(my_nqlwl)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,master,my_rank,nprocs
 real(dp) :: ucvol
 logical :: is_qeq0
 character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 complex(gwpc),allocatable :: chitmp(:,:)
 complex(gwpc),pointer :: vc_sqrt(:)

! *************************************************************************

 ABI_UNUSED(chi0_head(1,1))
 ABI_UNUSED(chi0_lwing(1,1))
 ABI_UNUSED(chi0_uwing(1,1))

 if (nI/=1.or.nJ/=1) then
   MSG_ERROR("nI or nJ=/1 not yet implemented")
 end if

 ABI_CHECK(Vcp%nqlwl==1,"nqlwl/=1 not coded")
 ABI_CHECK(my_nqlwl==1,"my_nqlwl/=1 not coded")

 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master  = 0

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 is_qeq0 = (normv(Vcp%qibz(:,iqibz),gmet,'G')<GW_TOLQ0)

 if (iqibz==1) then
   !%vc_sqrt => Vcp%vcqlwl_sqrt(:,iqlwl)  ! Use Coulomb term for q-->0
   vc_sqrt => Vcp%vcqlwl_sqrt(:,1)  ! TODO add treatment of non-Analytic behavior
 else
   vc_sqrt => Vcp%vc_sqrt(:,iqibz)
 end if

 write(msg,'(a,f8.2,a)')" chitmp requires: ",npwe**2*gwpc*b2Mb," Mb"
 ABI_MALLOC(chitmp,(npwe,npwe))
 ABI_CHECK_ALLOC(msg)
 !
 ! * Calculate chi0*fxc.
 chitmp = MATMUL(chi0,kxcg_mat)
 ! * First calculate the NLF contribution
 do ig1=1,npwe
   do ig2=1,npwe
     chitmp(ig1,ig2)=-chitmp(ig1,ig2)
   end do
   chitmp(ig1,ig1)=chitmp(ig1,ig1)+one
 end do

 call xginv(chitmp,npwe,comm=comm)

 chitmp = MATMUL(chitmp,chi0)
 !if (.not. ABS(REAL(omega))> tol3) then
 !  call hermitianize(chitmp,"All")
 !end if
 chitmp(1,1)=-vc_sqrt(1)*chitmp(1,1)*vc_sqrt(1)
 chitmp(1,1)=chitmp(1,1)+one

 epsm_nlf(1)=chitmp(1,1)

 chitmp = MATMUL(chi0,kxcg_mat)
 ! * Calculate (1-chi0*Vc-chi0*Kxc) and put it in chitmp.
 do ig1=1,npwe
   do ig2=1,npwe
     chitmp(ig1,ig2)=-chitmp(ig1,ig2)-chi0(ig1,ig2)*vc_sqrt(ig2)**2
   end do
   chitmp(ig1,ig1)=chitmp(ig1,ig1)+one
 end do

 ! * Invert (1-chi0*Vc-chi0*Kxc) and Multiply by chi0.
 call xginv(chitmp,npwe,comm=comm)
 chitmp=MATMUL(chitmp,chi0)

 ! * Save result, now chi0 contains chi.
 chi0=chitmp

 SELECT CASE (option_test)

 CASE (0) ! Symmetrized TESTPARTICLE epsilon^-1
   call wrtout(std_out,' Calculating TESTPARTICLE epsilon^-1(G,G") = 1 + Vc*chi','COLL')
   do ig1=1,npwe
     chi0(ig1,:)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:)
     chi0(ig1,ig1)=one+chi0(ig1,ig1)
   end do

 CASE (1) ! Symmetrized TESTELECTRON epsilon^-1
   call wrtout(std_out,' Calculating TESTELECTRON epsilon^-1(G,G") = 1 + (Vc + fxc)*chi',"COLL")
   chitmp=MATMUL(kxcg_mat,chi0)
   ! Perform hermitianization, only valid along the imaginary axis.
   if (.not. ABS(REAL(omega))> tol3) then
     call hermitianize(chitmp,"All")
   end if
   do ig1=1,npwe
     chi0(ig1,:)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:)+chitmp(ig1,:)
     chi0(ig1,ig1)=one+chi0(ig1,ig1)
   end do

 CASE DEFAULT
   write(msg,'(a,i3)')'Wrong value for option_test= ',option_test
   MSG_BUG(msg)
 END SELECT

 ABI_FREE(chitmp)
 !
 ! === chi0 now contains symmetrical epsm1 ===
 ! * Calculate macroscopic dielectric constant epsm_lf(w)=1/epsm1(G=0,Gp=0,w).
 epsm_lf(1) =  one/chi0(1,1)
 eelf   (1) = -AIMAG(chi0(1,1))

 !write(msg,'(a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at omega',omega*Ha_eV,' [eV]'
 !call wrtout(std_out,msg,'COLL')
 !call print_arr(chi0(:,:,io),mode_paral='PERS')

end subroutine atddft_symepsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/make_W
!! NAME
!! make_W
!!
!! FUNCTION
!!  calculate the symmetrical inverse dielectric matrix starting
!!  from the irreducible polarizability. The routine considers a single q-point, and
!!  performs the following tasks:
!!
!! INPUTS
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!   %nqibz=Number of q-points.
!!   %qibz(3,nqibz)=q-points in the IBZ.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix. See also SIDE EFFECTS.
!!
!! SIDE EFFECTS
!!  Er%epsm1(npwe*nI,npwe*nJ,nomega): in input the symmetryzed inverse dieletric matrix,
!!   in output the screened interaction (including a possible cutoff in real space)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_W(Er,Vcp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_W'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(vcoul_t),target,intent(in) :: Vcp
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: iq_ibz,ig1,ig2,io
 real(dp) :: ucvol
 logical :: is_qeq0
 character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 complex(gwpc),ABI_CONTIGUOUS pointer :: vc_sqrt(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (Er%nI/=1.or.Er%nJ/=1) then
   MSG_ERROR("nI or nJ=/1 not yet implemented")
 end if

 if (Er%ID/=4) then
   write(msg,'(a,i0,a)')" found Er%ID = ",Er%ID," while it should be 4"
   MSG_ERROR(msg)
 end if

 !TODO mqmem==0, change info ER%, and update Hscr
 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 do iq_ibz=1,Er%nqibz
   is_qeq0 = (normv(Er%qibz(:,iq_ibz),gmet,'G')<GW_TOLQ0) ! Check if q==0

   if (is_qeq0) then
     vc_sqrt => Vcp%vcqlwl_sqrt(:,1)  ! Use Coulomb term for q-->0, only first Q is used, shall we average if nqlwl>1?
   else
     vc_sqrt => Vcp%vc_sqrt(:,iq_ibz)
   end if

   do io=1,Er%nomega
     do ig1=1,Er%npwe
       do ig2=1,Er%npwe
         Er%epsm1(ig1,ig2,io,iq_ibz) = Er%epsm1(ig1,ig2,io,iq_ibz) * vc_sqrt(ig1) * vc_sqrt(ig2)
       end do
     end do
   end do
   !
 end do ! nqibz

end subroutine make_W
!!***

!----------------------------------------------------------------------

!!****f* m_screening/mkem1_q0
!! NAME
!! mkem1_q0
!!
!! FUNCTION
!!   This routine construct the microscopic dieletric matrix for q-->0 starting from the heads, wings and the body
!!   of the irreducible polarizability. Afterwards it calculates the symmetrized inverse dieletric matrix
!!   via a block wise inversion thus obtaining the heads and the wings of e^{-1} that can be
!!   used to describe the non-analytic behavior for q-->0.
!!
!! INPUTS
!! npwe=Number of Gs used to describe chi0
!! nomega=Number of frequencies in chi0.
!! n1,n2=Factors used to define the same of the chi0 matrix (1,1 if collinear, the typical case)
!! Cryst<crystal_t>=Structure describing the crystal structure.
!! Vcp<vcoul_t>=datatypes gathering info on the Coulomb term
!! gvec(3,npwe)=G-vector for chi0 in reduced coordinates.
!!
!! OUTPUT
!! eps_head(3,3,nomega)=The macroscopic dieletric tensor in reduced coordinates.
!!   The dieletric matrix along versor \hat q can be obtained with
!!     e(\hat q) = \hat q.eps_head \hat q if all quantities are given in Cartesian coordinates.
!!
!! SIDE EFFECTS
!! chi0(npwe*n1,npwe*n2,nomega)= Input: polarizability. output: inverse dieletric matrix (only the body is used)
!! chi0_lwing(npwe*n1,nomega,3)
!! chi0_uwing(npwe*n2,nomega,3)  Input:  the lower and upper wings of the polarizability
!!                               Output: the "lower" and "upper" wings of the inverse dieletric matrix. See notes below.
!! chi0_head(3,3,nomega)= Input: the polarizability tensor in Cartesian coordinates.
!!                        Ouput: The "head" of the inverse dieletric matrix. See notes below.
!!
!! NOTES
!!  Matrix inversion in block form.
!!
!!         1  n-1
!!  M =  | c  u^t| 1     ==>   M^{-1} =  |  1/k          -u^t A^{-1}/k                    |
!!       | v  A  | n-1                   | -A^{-1} v/k    A^{-1} + (A^{-1}v u^t A^{-1})/k |
!!
!!                             where k = c - u^t A^{-1} v
!!
!!  Let q be a versor in reciprocal space, the symmetrized dielectric matrix with bare coulomb interaction
!!  can be written as
!!
!!  \tilde\epsilon = | q.Eq      q.Wl(G2) |  where   E_ij = \delta_ij -4\pi chi0_head_ij
!!                   | q.Wl(G1)  B(G1,G2  |          Wl(G1) = -4\pi chi0_lwing(G1)
!!                                                   Wu(G2) = -4\pi chi0_uwing(G1)
!!  therefore, in Cartesian coordinates, we have:
!!
!!  1) e^{-1}_{00}(q) = [ q_i q_j (E_{ij} - \sum_{GG'} Wu_i(G)a B_{GG'}^{-1} Wl_j(G')) ]^{-1} = 1/(q.Lq)
!!
!!  2) e^{-1}_{0G'}(q) = -e^{-1}_{00}(q) [ \sum_{iG} q_i Wu_i(G)a B_{GG'}^{-1} ] = (q.Su) /(q.Lq)
!!
!!  3) e^{-1}_{G0}(q)  = -e^{-1}_{00}(q) [ \sum_{iG'} q_i B_{GG'}^{-1} Wl_i(G') ] = (q.Sl) /(q.Lq)
!!
!!  4) e^{-1}_{GG'}(q) = B_{GG'}^{-1} +
!!     [ \sum_{ij} q_i q_j ( \sum_T B^{-1}_{GT}^{-1} Wl_i(T)) (\sum_T' Wu_j(T') B^{-1}_{T'G'}^{-1} ] / (q.Lq)
!!
!!  where Su(G,3) and Sl(G,3) are the "upper" and "lower" wings of the inverse dielectric matrix and
!!  L is the inverse dielectric tensor. Similar equations hold even if vectors and tensors are given in terms
!!  of the reciprocal lattice vectors provided that the metric tensor is taken into account.
!!  The main difference is in the expression for the tensor as only one metric tensor can be
!!  absorbed in the scalar product, the second metric multiplies one of the wings.
!!
!!  *) The present implementation assumes that no cutoff technique is used in the Coulomb term.
!!
!!  *) Once the tensor in know it is possible to average the quadratic form on the sphere exactly.
!!  In Cartesian coordinates one obtains.
!!
!!    \dfrac{1}{4\pi} \int v.Tv d\Omega = Trace(T)/3
!!
!!  For the inverse dielectric matrix we have to resort to a numerical integration
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkem1_q0(npwe,n1,n2,nomega,Cryst,Vcp,gvec,chi0_head,chi0_lwing,chi0_uwing,chi0,eps_head)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkem1_q0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe,nomega,n1,n2
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
!arrays
 integer,intent(in) :: gvec(3,npwe)
 complex(gwpc),intent(inout) :: chi0(npwe*n1,npwe*n2,nomega)
 complex(dpc),intent(inout) :: chi0_lwing(npwe*n1,nomega,3)
 complex(dpc),intent(inout) :: chi0_uwing(npwe*n2,nomega,3)
 complex(dpc),intent(inout) :: chi0_head(3,3,nomega)
 complex(dpc),intent(out) :: eps_head(3,3,nomega)

!Local variables ------------------------------
!scalars
 integer :: iomega,ig,ig1,ig2,idir,jdir
!arrays
 real(dp),allocatable :: modg_inv(:)
 complex(dpc),allocatable :: eps_lwing(:,:),eps_uwing(:,:),eps_body(:,:),cvec(:)

!************************************************************************

 ABI_CHECK(npwe/=1,"npwe must be >1")

 ! Precompute 1/|G|.
 ABI_MALLOC(modg_inv,(npwe-1))
 do ig=1,npwe-1
   modg_inv(ig) = one/normv(gvec(:,ig+1),Cryst%gmet,'G')
 end do

 ABI_MALLOC(eps_uwing,((npwe-1)*n1,3))
 ABI_MALLOC(eps_lwing,((npwe-1)*n2,3))
 ABI_MALLOC(eps_body,(npwe-1,npwe-1))
 ABI_MALLOC(cvec,(npwe-1))

 do iomega=1,nomega
   !
   ! Head and wings of the symmetrized epsilon.
   eps_head(:,:,iomega) = -four_pi*chi0_head(:,:,iomega)
   do idir=1,3
     eps_head(idir,idir,iomega) = one + eps_head(idir,idir,iomega)
     eps_lwing(:,idir) = -four_pi * modg_inv * chi0_lwing(2:,iomega,idir)
     eps_uwing(:,idir) = -four_pi * modg_inv * chi0_uwing(2:,iomega,idir)
     !eps_lwing(:,idir) = -chi0_lwing(2:,iomega,idir) * SQRT(four_pi) * Vcp%vcqlwl_sqrt(2:npwe,1)
     !eps_uwing(:,idir) = -chi0_uwing(2:,iomega,idir) * SQRT(four_pi) * Vcp%vcqlwl_sqrt(2:npwe,1)
   end do

   write(std_out,*)" espilon head"
   call print_arr(eps_head(:,:,iomega))
   !
   ! Construct the body of the symmetrized epsilon then invert it.
   do ig2=1,npwe-1
     do ig1=1,npwe-1
       eps_body(ig1,ig2) = -four_pi * modg_inv(ig1)*chi0(ig1+1,ig2+1,iomega )*modg_inv(ig2)
       !eps_body(ig1,ig2) = -Vcp%vcqlwl_sqrt(ig1+1,1)*chi0(ig1+1,ig2+1,iomega)* Vcp%vcqlwl_sqrt(ig2+1,1)
     end do
     eps_body(ig2,ig2) = one + eps_body(ig2,ig2)
   end do

   call xginv(eps_body,npwe-1)
   !
   ! Overwrite chi0_head and chi0_wings with the head and the wings of the inverse dielectric matrix.
   do jdir=1,3
     !
     ! Head.
     cvec=czero
     do idir=1,3
       cvec = cvec + two_pi*Cryst%gmet(jdir,idir)*MATMUL(eps_body,eps_lwing(:,idir)) ! as we work in reciprocal coords.
     end do
     !cvec = MATMUL(eps_body,eps_lwing(:,jdir))
     do idir=1,3
       chi0_head(idir,jdir,iomega) = eps_head(idir,jdir,iomega) - xdotu(npwe-1,eps_uwing(:,idir),1,cvec,1)
     end do
     !
     ! Now the wings.
     chi0_uwing(2:,iomega,jdir) = -MATMUL(eps_uwing(:,jdir),eps_body)
     chi0_lwing(2:,iomega,jdir) = -MATMUL(eps_body,eps_lwing(:,jdir))
     !
   end do !jdir

   write(std_out,*)"espilon^1 head after block inversion"
   call print_arr(chi0_head(:,:,iomega))
   !
   ! Change the body but do not add the corrections due to the head and the wings.
   ! since they can be obtained on the fly from eps_body and the wings of eps^{-1}.
   !%chi0(2:,2:,iomega) = eps_body
 end do !iomega

 ABI_FREE(modg_inv)
 ABI_FREE(cvec)
 ABI_FREE(eps_lwing)
 ABI_FREE(eps_uwing)
 ABI_FREE(eps_body)

 RETURN
 ABI_UNUSED(Vcp%ng)

end subroutine mkem1_q0
!!***

!----------------------------------------------------------------------

!!****f* m_screening/lebedev_laikov_int
!! NAME
!!  lebedev_laikov_int
!!
!! FUNCTION
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

subroutine lebedev_laikov_int()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lebedev_laikov_int'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
!scalars
 integer :: on,npts,ii,ll,mm,lmax,leb_idx !ierr,
 real(dp) :: accuracy
 complex(dpc) :: ang_int
!arrays
 real(dp) :: cart_vpt(3) !,real_pars(0)
 real(dp),allocatable :: vx(:),vy(:),vz(:),ww(:)
 complex(dpc) :: tensor(3,3),cplx_pars(9)
 complex(dpc),allocatable :: ref_func(:),expd_func(:) !tmp_momenta(:)

! *************************************************************************

 MSG_ERROR("lebedev_laikov_int is still under development")

 !tensor=RESHAPE((/4.0,2.0,4.0,0.5,2.1,0.0,5.4,2.1,5.0/),(/3,3/))
 tensor=RESHAPE((/4.0,0.0,0.0,0.0,4.0,0.0,0.0,0.0,5.0/),(/3,3/))
 !tensor=RESHAPE((/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/),(/3,3/))

 npts=26
 ABI_MALLOC(vx,(npts))
 ABI_MALLOC(vy,(npts))
 ABI_MALLOC(vz,(npts))
 ABI_MALLOC(ww,(npts))

 call LD0026(vx,vy,vz,ww,on)

 ang_int=czero
 do ii=1,npts
   cart_vpt = (/vx(ii),vy(ii),vz(ii)/)
   ang_int = ang_int + ww(ii)*DOT_PRODUCT(cart_vpt,MATMUL(tensor,cart_vpt))
 end do

 !write(std_out,*)"quadratic form associated to tensor=",tensor
 write(std_out,*)"on ang_int",on,ang_int

 ABI_FREE(vx)
 ABI_FREE(vy)
 ABI_FREE(vz)
 ABI_FREE(ww)

 call init_lebedev_gridset()
 cplx_pars = RESHAPE(tensor,(/9/)); accuracy=tol10

 ! This is the function to be expanded evaluated on the lebedev_laikov grid of index leb_idx
 leb_idx=3; npts=lebedev_npts(leb_idx)
 ABI_MALLOC(ref_func,(npts))
 do ii=1,npts
   cart_vpt = Lgridset(leb_idx)%versor(:,ii)
   ref_func(ii) = one/DOT_PRODUCT(cart_vpt,MATMUL(tensor,cart_vpt))
 end do

 ! Calculate the expansion in angular momenta of 1/{q.Tq}.
 ! Only even l-components contribute thanks to the parity of the integrand.
 ! tol6 seems to be an acceptable error, convergence wrt lmax is very slow even for simple tensors.
 ABI_MALLOC(expd_func,(npts))
 expd_func=czero
 lmax=10
 do ll=0,lmax,2
   !allocate(tmp_momenta(-ll:ll))
   do mm=-ll,ll
     ! MG: Commented becase it causes problems with the new version of abilint
     !call lebedev_quadrature(ylmstar_over_qTq,(/ll,mm/),real_pars,cplx_pars,ang_int,ierr,accuracy)
     write(std_out,*)ll,mm,ang_int
     !tmp_momenta(mm) = ang_int
     do ii=1,npts
       cart_vpt = Lgridset(leb_idx)%versor(:,ii)
       expd_func(ii) = expd_func(ii) + four_pi*ang_int*ylmc(ll,mm,cart_vpt)
     end do
   end do
   !deallocate(tmp_momenta)
   write(std_out,*)"Error in angular expansion at l=",ll," is ",MAXVAL(ABS(expd_func-ref_func))
 end do

!BEGINDEBUG
 do ii=1,npts
   write(777,*)ref_func(ii)
   write(778,*)expd_func(ii)
 end do
!ENDDEBUG

 ABI_FREE(expd_func)
 ABI_FREE(ref_func)
 call destroy_lebedev_gridset()

 MSG_ERROR("Exiting from lebedev_laikov_int")

end subroutine lebedev_laikov_int
!!***

!----------------------------------------------------------------------

!!****f* m_screening/ylmstar_over_qTq
!! NAME
!!  ylmstar_over_qTq
!!
!! FUNCTION
!!  Return Ylm(q)^*/(q.Tq) where q is a versor in Cartesian coordinates.
!!  and Ylm is a complex spherical Harmonics whose index (l,m) are
!!  passed via int_pars(1:2). T is a tensore in Cartesian coordinates
!!  passed via cplx_pars(1:9).
!!
!! INPUTS
!!  cart_vers(3)=Cartesian components of the versor
!!  int_pars(1:2)=(l,m) indeces in Ylm. l>=0 and  m \in [-l,l]
!!  cplx_pars(1:9)=Tensor T in Cartesian coordinates.
!!  real_pars=Not used.
!!
!! OUTPUT
!!  Value of Ylm(q)^*/(q.Tq)
!!
!! PARENTS
!!
!! SOURCE

function ylmstar_over_qTq(cart_vers,int_pars,real_pars,cplx_pars)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ylmstar_over_qTq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: cart_vers(3)
 integer,intent(in) :: int_pars(:)
 real(dp),intent(in) :: real_pars(:)
 complex(dpc),intent(in) :: cplx_pars(:)
 complex(dpc) :: ylmstar_over_qTq
!arrays

!Local variables-------------------------------
!scalars
 integer :: ll,mm
!arrays
 complex(dpc) :: tensor(3,3)

! *************************************************************************

 tensor = RESHAPE(cplx_pars(1:9),(/3,3/))
 ll = int_pars(1) ! ll starts from zero.
 mm = int_pars(2) ! m \in [-l,l]

 ylmstar_over_qTq = CONJG(ylmc(ll,mm,cart_vers))/DOT_PRODUCT(cart_vers,MATMUL(tensor,cart_vers))

 RETURN
 ABI_UNUSED(real_pars(1))

end function ylmstar_over_qTq
!!***

!----------------------------------------------------------------------

!!****f* m_screening/ylmstar_wtq_over_qTq
!! NAME
!!  ylmstar_wtq_over_qTq
!!
!! FUNCTION
!!  Return Ylm(q)^* weight(q)/(q.Tq) where q is a versor in Cartesian coordinates.
!!  Ylm is a complex spherical Harmonics whose index (l,m) are
!!  passed via int_pars(1:2). T is a tensor in Cartesian coordinates
!!  passed via cplx_pars(1:9). weight(q) is the weighting function giving
!!  the length of the vector parallel to versor q that connects the origin
!!  of the lattice to one of the boundaries of the small cell centered at Gamma
!!
!! INPUTS
!!  cart_vers(3)=Cartesian components of the versor
!!  int_pars(1:2)=(l,m) indeces in Ylm. l>=0 and  m \in [-l,l]
!!  cplx_pars(1:9)=Tensor T in Cartesian coordinates.
!!  real_pars(1:9)=The Cartesian vectors defining the small box centered around gamma point
!!    when referred to this vectors the points in the box are given by {(x,y,z) | x,y,z \in [-1,1]}.
!!
!! OUTPUT
!!  Value of Ylm(q)^* weigh(q)/(q.Tq)
!!
!! PARENTS
!!
!! SOURCE

function ylmstar_wtq_over_qTq(cart_vers,int_pars,real_pars,cplx_pars)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ylmstar_wtq_over_qTq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: cart_vers(3)
 integer,intent(in) :: int_pars(:)
 real(dp),intent(in) :: real_pars(:)
 complex(dpc),intent(in) :: cplx_pars(:)
 complex(dpc) :: ylmstar_wtq_over_qTq
!arrays

!Local variables-------------------------------
!scalars
 integer :: ll,mm
 real(dp) :: wtq
!arrays
 real(dp) :: gprimd(3,3),rprimd(3,3),red_vers(3)
 complex(dpc) :: tensor(3,3)

! *************************************************************************

 MSG_ERROR("Work in progress")
 ! box_len has to be tested

 gprimd = RESHAPE(real_pars(1:9),(/3,3/))
 red_vers = MATMUL(rprimd,cart_vers)
 wtq = box_len(red_vers,gprimd)

 tensor = RESHAPE(cplx_pars(1:9),(/3,3/))
 ll = int_pars(1) ! true ll i.e. not shifted
 mm = int_pars(2)

 ylmstar_wtq_over_qTq = CONJG(ylmc(ll,mm,cart_vers))*wtq/DOT_PRODUCT(cart_vers,MATMUL(tensor,cart_vers))

end function ylmstar_wtq_over_qTq
!!***

!----------------------------------------------------------------------

!!****f* m_screening/k_fermi
!! NAME
!!  k_fermi
!!
!! FUNCTION
!!  Returns the Fermi wave vector corresponding to the local value of the real space density rhor.
!!
!! INPUTS
!!  rhor=Local density in real space.
!!
!! PARENTS
!!
!! SOURCE

elemental function k_fermi(rhor)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'k_fermi'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rhor
 real(dp) :: k_fermi
!arrays

!Local variables-------------------------------
!scalars
 real(dp),parameter :: pisq=pi**2

! *************************************************************************

 k_fermi = (three*pisq*rhor)**third

end function k_fermi
!!***

!----------------------------------------------------------------------

!!****f* m_screening/k_thfermi
!! NAME
!!  k_thfermi
!!
!! FUNCTION
!!  Returns the Thomas-Fermi wave vector corresponding to the local value of the real space density rhor.
!!
!! INPUTS
!!  rhor=Local density in real space.
!!
!! PARENTS
!!
!! SOURCE

elemental function k_thfermi(rhor)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'k_thfermi'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rhor
 real(dp) :: k_thfermi

!Local variables-------------------------------
!scalars
 real(dp),parameter :: pisq=pi**2

! *************************************************************************

 k_thfermi = SQRT(four*k_fermi(rhor)*piinv)

end function k_thfermi
!!***

!----------------------------------------------------------------------

!!****f* m_screening/mdielf_bechstedt
!! NAME
!!  mdielf_bechstedt
!!
!! FUNCTION
!!  Calculates the model dielectric function for the homogeneous system
!!  as proposed by F. Bechstedt, in Solid State Commun. 84, 765 1992.
!!
!! INPUTS
!!  eps_inf=Dielectric constant of the material
!!  qnrm=The modulus of the q-point.
!!  rhor=The local value of the density
!!
!! PARENTS
!!
!! SOURCE

elemental function mdielf_bechstedt(eps_inf,qnrm,rhor) result(mdielf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mdielf_bechstedt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: eps_inf,qnrm,rhor
 real(dp) :: mdielf

! *************************************************************************

 mdielf = one + &
&  one / ( one/(eps_inf-one) + (qnrm/k_thfermi(rhor))**2 + (three*qnrm**4)/(four*k_fermi(rhor)**2 * k_thfermi(rhor)**2) )

end function mdielf_bechstedt
!!***

!----------------------------------------------------------------------

!!****f* m_screening/screen_mdielf
!! NAME
!!  screen_mdielf
!!
!! FUNCTION
!!  Calculates W_{G,G'}(q,w) for a given q-point in the BZ using a model dielectric function.
!!
!! INPUTS
!!  iq_bz=The index of the q-point in the BZ where W(q) is calculated.
!!  npw=Number of plane waves for W
!!  nomega=Number of frequency points.
!!  model_type=Flag defining the model.
!!  eps_inf=Dielectric constant of the material.
!!  Cryst<crystal_t>=Info on the unit cell
!!  Qmesh<kmesh_t>=Info on the set of q-points.
!!  Vcp<vcoul_t datatype>= containing information on the cutoff technique
!!  Gsph<Gsphere>=The G-sphere for W.
!!  nspden=Number of spin density components of the density.
!!  nfft=Number of FFT points on the dense FFT mesh
!!  ngfft(18)=contain all needed information about 3D FFT.
!!  rhor(nfft,nspden)=Electron density in real space (The PAW AE term is included)
!!  which= Set it to "EM1" if the symmetrized inverse dielectric matrix is wanted.
!!   By default the routines returns W.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  w_qbz(npw,npw,nomega)
!!
!! NOTES
!!   W_{G1,G2} =  1/2 {
!!     v(q+G1) \int em1(|q+G1|,r) e^{-i(G1-G2).r} dr  +
!!     v(q+G2) \int em1(|q+G2|,r) e^{-i(G1-G2).r} dr } / \Omega
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!
!! SOURCE

subroutine screen_mdielf(iq_bz,npw,nomega,model_type,eps_inf,Cryst,Qmesh,Vcp,Gsph,nspden,nfft,ngfft,rhor,which,w_qbz,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_mdielf'
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nomega,nfft,nspden,iq_bz,comm,model_type
 real(dp),intent(in) :: eps_inf
 character(len=*),intent(in) :: which
 type(kmesh_t),intent(in) :: Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),target,intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,nspden)
 complex(gwpc),intent(out) :: w_qbz(npw,npw,nomega)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0,paral_kgb0=0,cplex1=1
 integer :: my_gstart,my_gstop,iq_ibz,ig,itim_q,isym_q
 integer :: ig1,ig2,g1mg2_fft,iw,ii,ierr,nprocs,isg,ifft !,row ,col
 real(dp) :: qpg2_nrm
 complex(dpc) :: ph_mqbzt
 logical :: is_qeq0,isirred
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: umklp(3)
 integer,allocatable :: igfft(:),g1mg2(:,:)
 real(dp) :: qpg2(3),qpt_bz(3)
 real(dp),allocatable :: em1_qpg2r(:),fofg(:,:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: vc_sqrt_ibz(:)
 complex(gwpc),allocatable :: vc_qbz(:),ctmp(:,:)
 logical,allocatable :: mask(:)

! *************************************************************************

 ABI_CHECK(nomega==1,"screen_mdielf does not support nomega>1")

!Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft(2),ngfft(3),'all')

 nprocs = xcomm_size(comm)
 call xmpi_split_work(npw,comm,my_gstart,my_gstop,msg,ierr)

 call get_bz_item(Qmesh,iq_bz,qpt_bz,iq_ibz,isym_q,itim_q,ph_mqbzt,umklp,isirred)

 !if (itim_q/=1.or.isym_q/=1.or.ANY(umklp/=0) ) then
 !  MSG_ERROR("Bug in mdielf_bechstedt")
 !end if
 !
 ! Symmetrize Vc in the full BZ.
 is_qeq0 = (normv(qpt_bz,Cryst%gmet,'G')<GW_TOLQ0) ! Check if q==0
 if (is_qeq0) then
   vc_sqrt_ibz => Vcp%vcqlwl_sqrt(:,1)  ! Use Coulomb term for q-->0, only first Q is used, shall we average if nqlwl>1?
 else
   vc_sqrt_ibz => Vcp%vc_sqrt(:,iq_ibz)
 end if

 ABI_MALLOC(vc_qbz,(npw))
 do ig=1,npw
   isg = Gsph%rottb(ig,itim_q,isym_q)
   vc_qbz(isg) = vc_sqrt_ibz(ig)**2
 end do

 ABI_MALLOC(igfft,(npw))
 ABI_MALLOC(g1mg2,(3,npw))
 ABI_MALLOC(fofg,(2,nfft))
 ABI_MALLOC(em1_qpg2r,(nfft))
 ABI_MALLOC(mask,(npw))

 w_qbz=czero
 do ig2=my_gstart,my_gstop
   !
   ! Compute the index of G-G2 wave in the FFT grid.
   do ii=1,npw
     g1mg2(:,ii) = Gsph%gvec(:,ii) - Gsph%gvec(:,ig2)
   end do
   call kgindex(igfft,g1mg2,mask,MPI_enreg_seq,ngfft,npw)

   ! TODO can use zero-padding FFT to speed up the transform.
   !call sphereboundary(gbound,istwfk1,g1mg2,mgfft,npw)

   ! Evaluate em1_qpg2r = \int em1(|q+G2|,r) e^{-i(G1-G2).r} dr }.
   qpg2 = qpt_bz + Gsph%gvec(:,ig2)
   qpg2_nrm = normv(qpg2,Cryst%gmet,"G")

   do iw=1,nomega
     !
     select case (model_type)
     case (1)
       do ifft=1,nfft
         em1_qpg2r(ifft) = one / mdielf_bechstedt(eps_inf,qpg2_nrm,rhor(ifft,1))
       end do
     case default
       write(msg,'(a,i0)')" Unknown value for model_type ",model_type
       MSG_ERROR(msg)
     end select

     call fourdp(cplex1,fofg,em1_qpg2r,-1,MPI_enreg_seq,nfft,ngfft,paral_kgb0,tim_fourdp0)
     !
     ! Here, unlike the other parts of the code, the unsymmetrized e^{-1} is used.
     do ig1=1,npw
       g1mg2_fft = igfft(ig1)
       w_qbz(ig1,ig2,iw) = DCMPLX(fofg(1,g1mg2_fft), fofg(2,g1mg2_fft)) * vc_qbz(ig2) !/ Cryst%ucvol
     end do
   end do ! iw
   !
 end do ! ig2

 ABI_FREE(em1_qpg2r)
 ABI_FREE(fofg)
 ABI_FREE(igfft)
 ABI_FREE(g1mg2)
 ABI_FREE(mask)
 !
 ! W = 1/2 * (A + A^H)
 ! The MPI sum is done inside the loop to avoid problems with the size of the packet.
 ABI_MALLOC(ctmp,(npw,npw))
 ABI_CHECK_ALLOC("out of memory in ctmp")

 do iw=1,nomega
   !ctmp = TRANSPOSE(CONJG(w_qbz(:,:,iw)))
   ctmp = GWPC_CONJG(w_qbz(:,:,iw))
   call sqmat_itranspose(npw,ctmp)
   w_qbz(:,:,iw) = half * (ctmp + w_qbz(:,:,iw))
   call xmpi_sum(w_qbz(:,:,iw),comm,ierr)
 end do
 !
 ! Calculate the symmetrized Em1. W = vc(G1)^{1/2} \tilde Em1 vc(G2)^{1/2} -------------------------
 if (toupper(which)=="EM1") then
   do ig=1,npw
     isg = Gsph%rottb(ig,itim_q,isym_q)
     vc_qbz(isg) = vc_sqrt_ibz(ig)  ! Workspace storing vc*{1/2}(q_BZ,G).
   end do

   do ig2=1,npw
     do ig1=1,npw
       ctmp(ig1,ig2) =  one / (vc_qbz(ig1) * vc_qbz(ig2))
     end do
   end do

   do iw=1,nomega
     w_qbz(:,:,iw) = w_qbz(:,:,iw) * ctmp(:,:)
   end do
 end if

 call destroy_mpi_enreg(MPI_enreg_seq)

 ABI_FREE(vc_qbz)
 if (allocated(ctmp)) then
   ABI_FREE(ctmp)
 end if

end subroutine screen_mdielf
!!***

!!****f* m_screening/recalculate_epsm1_freq_grid
!! NAME
!! recalculate_epsm1_freq_grid
!!
!! FUNCTION
!!  Recalculate the frequency gridpoints in the Epsilonm1_results structure.
!!  This is useful when a pole-fit screening is used and the grid can be
!!  arbitrarily chaged for the sigma calculation.
!!
!! INPUTS
!!  Er        - The Epsilonm1_results structure
!!  nfreqre   - Number of real frequencies
!!  nfreqim   - Number of imaginary frequencies
!!  freqremax - Maximum real frequency
!!  freqremin - Minimum real frequency
!!
!! OUTPUT
!!   The Epsilonm1_results structure with changed grid parameters
!!
!! NOTES
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine recalculate_epsm1_freq_grid(Er,nfreqre,nfreqim,freqremin,freqremax,ppmfrq,freqim_alpha)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'recalculate_epsm1_freq_grid'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfreqre,nfreqim
 real(dp),intent(in) :: freqremin,freqremax,ppmfrq,freqim_alpha
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: iomega
 real(dp) :: domegareal

! *************************************************************************

 Er%nomega   = nfreqre+nfreqim
 Er%nomega_i = nfreqim
 Er%nomega_r = nfreqre
 Er%nomega_c = 0
 if (allocated(Er%omega))  then
   ABI_FREE(Er%omega)
 end if
 ABI_MALLOC(Er%omega,(Er%nomega))

! Real frequencies
 Er%omega(1)=CMPLX(freqremin,zero,kind=dpc)

 if (Er%nomega_r>1) then ! Avoid division by zero.
   domegareal=(freqremax-freqremin)/(Er%nomega_r-1)
   do iomega=2,Er%nomega_r
     Er%omega(iomega)=CMPLX(freqremin+(iomega-1)*domegareal,zero,kind=dpc)
   end do
 end if

! Imaginary frequencies
 do iomega=1,Er%nomega_i
   Er%omega(Er%nomega_r+iomega)=CMPLX(zero,ppmfrq/(freqim_alpha-two)*(EXP(two/(Er%nomega_i+1)&
&   *LOG(freqim_alpha-one)*iomega)-one),kind=dpc)
 end do

end subroutine recalculate_epsm1_freq_grid
!!***

END MODULE m_screening
!!***
