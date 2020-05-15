!{\src2tex{textfont=tt}}
!!****f* ABINIT/pspatm_pspio
!! NAME
!!  pspatm_pspio
!!
!! FUNCTION
!!  Open atomic pseudopotential data file for a given atom,
!!  read the three first lines, make some checks, then
!!  call appropriate subroutine for the reading of
!!  V(r) and wf R(r) data for each angular momentum, and subsequent
!!  Fourier and Bessel function transforms for local and nonlocal potentials.
!!  Close psp file at end.
!!
!!  Handles pseudopotential files produced by (pspcod=1 or 4) Teter code,
!!  or from the Goedecker-Teter-Hutter paper (pspcod=2),
!!  or from the Hartwigsen-Goedecker-Hutter paper (pspcod=3 or 10)
!!  or "Phoney pseudopotentials" (Hamman grid in real space) (pspcod=5)
!!  or "Troullier-Martins pseudopotentials" from the FHI (pspcod=6)
!!  or "XML format" (pspcod=9)
!!  or "UPF PWSCF format" (pspcod=11)
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (YP)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dq= spacing of the q-grid
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | ixc=exchange-correlation choice from main routine data file
!!   | pawxcdev=choice of XC development in PAW formalism
!!   | usexcnhat=choice for use of nhat in Vxc in PAW formalism
!!   | xclevel= XC functional level
!!  ipsp=id in the array of the currently read pseudo.
!!
!! OUTPUT
!!  ekb(dimekb)=
!!    ->NORM-CONSERVING PSPS ONLY (pspcod/=7):
!!      (Real) Kleinman-Bylander energies (hartree)
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for number of basis functions (l,n) (dimekb=lnmax)
!!             If any, spin-orbit components begin at l=mpsang+1
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!  xccc1d(n1xccc*(1-usepaw),6)=1D core charge function and five derivatives,
!!                              from psp file (used in NC only)
!!
!! SIDE EFFECTS
!!  Input/Output :
!!   psps <type(pseudopotential_type)>=at output, values depending on the read
!!                                     pseudo are set.
!!    | dimekb(IN)=dimension of ekb (see module defs_psp.f)
!!    | filpsp(IN)=name of formatted external file containing atomic psp data.
!!    | lmnmax(IN)=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!    |           =if useylm=0, max number of (l,n)   comp. over all type of psps
!!    | lnmax(IN)=max. number of (l,n) components over all type of psps
!!    |           angular momentum of nonlocal pseudopotential
!!    | mpsang(IN)= 1+maximum angular momentum for nonlocal pseudopotentials
!!    | mpssoang(IN)= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!    | mqgrid_ff(IN)=dimension of q (or G) grid for nl form factors (array ffspl)
!!    | mqgrid_vl(IN)=dimension of q (or G) grid for Vloc (array vlspl)
!!    | n1xccc(IN)=dimension of xccc1d ; 0 if no XC core correction is used
!!    | optnlxccc(IN)=option for nl XC core correction
!!    | positron(IN)=0 if electron GS calculation
!!    |              1 if positron GS calculation
!!    |              2 if electron GS calculation in presence of the positron
!!    | pspso(INOUT)=spin-orbit characteristics, govern the content of ffspl and ekb
!!    |          if =0 : this input requires NO spin-orbit characteristics of the psp
!!    |          if =2 : this input requires HGH characteristics of the psp
!!    |          if =3 : this input requires HFN characteristics of the psp
!!    |          if =1 : this input will be changed at output to 1, 2, 3, according
!!    |                  to the intrinsic characteristics of the psp file
!!    | qgrid_ff(mqgrid_ff)(IN)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!    | qgrid_vl(mqgrid_vl)(IN)=values of q on grid from 0 to qmax (bohr^-1) for Vloc
!!    | usepaw(IN)= 0 for non paw calculation; =1 for paw calculation
!!    | useylm(IN)=governs the way the nonlocal operator is to be applied:
!!    |            1=using Ylm, 0=using Legendre polynomials
!!    | vlspl_recipSpace(IN)=.true. if pseudo are expressed in reciprocal space.
!!    | znuclpsp(IN)=atomic number of atom as specified in input file to main routine
!!
!! NOTES
!!  Format expected for the three first lines of pseudopotentials
!!  (1) title (character) line
!!  (2) znucl,zion,pspdat
!!  (3) pspcod,pspxc,lmax,lloc,mmax,r2well
!!
!!  Dimensions of form factors and Vloc q grids must be the same in Norm-Conserving case
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#if defined HAVE_TRIO_LIBPSPIO
#include "pspio_common.h"
#endif

#include "abi_common.h"

subroutine pspatm_pspio(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&
&  psps,vlspl,dvlspl,xcccrc,xccc1d,comm_mpi)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
#if defined HAVE_TRIO_LIBPSPIO
 use pspio_f90_types_m
 use pspio_f90_lib_m
#endif

 use m_pawrad, only : pawrad_type
 use m_pawtab, only : pawtab_type
 use m_pawpsp, only : pawpsp_bcast

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspatm_pspio'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_65_psp, except_this_one => pspatm_pspio
!End of the abilint section

  implicit none

!Arguments ---------------------------------------------
!scalars
!no_abirules
  integer,intent(in) :: ipsp
  integer, optional,intent(in) :: comm_mpi
  real(dp),intent(in) :: dq
  real(dp),intent(out) :: epsatm,xcccrc
  type(dataset_type),intent(in) :: dtset
  type(datafiles_type),intent(in) :: dtfil
  type(pawrad_type),intent(out) :: pawrad
  type(pawtab_type),intent(out) :: pawtab
  type(pseudopotential_type),intent(inout) :: psps
!arrays
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
  real(dp),intent(inout) :: ekb(psps%dimekb*(1-psps%usepaw)) !vz_i
  real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax) !vz_i
  real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
  real(dp),intent(inout) :: xccc1d(psps%n1xccc*(1-psps%usepaw),6) !vz_i

!Local variables ---------------------------------------
!scalars
!no_abirules
  integer :: ii,il,ilmn,iln,iln0,lloc,lmax,me,mmax,paral_mode
  integer :: pio_fmt,pspcod,pspdat,psp_x,psp_c,pspxc
  real(dp) :: accuracy_p,qchrg,r2well,zion,znucl
  character(len=500) :: msg
  character(len=70) :: testxml
  character(len=fnlen) :: title
  character(len=fnlen) :: filnam
!arrays
  integer,allocatable :: nproj(:)
  real(dp) :: tsec(2)
  real(dp),allocatable :: e990(:),e999(:),ekb1(:),ekb2(:),epspsp(:),rcpsp(:)
  real(dp),allocatable :: rms(:)
#if defined HAVE_TRIO_LIBPSPIO
  type(pspheader_type)       :: pspheads
  type(pspio_f90_pspdata_t) :: pspdata
  type(pspio_f90_xc_t) :: pspio_xc
#endif

! ******************************************************************************

#if defined HAVE_TRIO_LIBPSPIO

 if (psps%usepaw == 1) then
   MSG_ERROR('PAW is not implemented in Libpspio.')
 end if

!paral_mode defines how we access to the psp file
!  paral_mode=0: all processes access to the file (sequentially)
!  paral_mode=1: only proc. 0 access to the file and then broadcast
 paral_mode=0
 if (present(comm_mpi)) then
   if (psps%usepaw==1.and.xcomm_size(comm_mpi)>1) paral_mode=1
 end if
 me=0;if (paral_mode==1) me=xcomm_rank(comm_mpi)

 if ( me == 0 ) then
!  Dimensions of form factors and Vloc q grids must be the same in
!  Norm-Conserving case
   if ( (psps%usepaw == 0) .and. (psps%mqgrid_ff /= psps%mqgrid_vl) ) then
     write(msg, '(a,a,a,a,a)' )&
&     'Dimension of q-grid for nl form factors (mqgrid_ff)',ch10,&
&     'is different from dimension of q-grid for Vloc (mqgrid_vl) !',ch10,&
&     'This is not allowed for norm-consevring psp.'
     MSG_ERROR(msg)
   end if

!  Allocate nproj here: can be read in now for UPF
   ABI_ALLOCATE(nproj, (psps%mpssoang))
   nproj(:) = 0

!  ----------------------------------------------------------------------------
!  Read psp info using Libpspio
   write(msg, '(a,t38,a)' )'- pspatm_pspio: opening atomic psp file',trim(psps%filpsp(ipsp))
   call wrtout(ab_out,  msg,'COLL')
   call wrtout(std_out,  msg,'COLL')

!  Init pspio data structure and parse file
   pio_fmt = PSPIO_FMT_UNKNOWN
   call pspio_check_error(pspio_f90_pspdata_init(pspdata))
   call pspio_check_error(pspio_f90_pspdata_read(pspdata, pio_fmt, &
&   trim(psps%filpsp(ipsp))))

!  Transfer relevant information to Abinit
   pspcod = 0
   r2well = 0
   call pspio_f90_pspdata_get_z(pspdata, znucl)
   call pspio_f90_pspdata_get_zvalence(pspdata, zion)
   call pspio_f90_pspdata_get_l_max(pspdata, lmax)
   call pspio_f90_pspdata_get_l_local(pspdata, lloc)
   call pspio_f90_pspdata_get_xc(pspdata, pspio_xc)

!  FIXME: ensure consistency of XC indices
   call pspio_f90_xc_get_exchange(pspio_xc, psp_x)
   call pspio_f90_xc_get_correlation(pspio_xc, psp_c)
   pspxc = - psp_x - 1000 * psp_c

!  FIXME: find better strategy
   select case ( pio_fmt )
     case (PSPIO_FMT_ABINIT_1)
       pspcod = 1
     case (PSPIO_FMT_ABINIT_2)
       pspcod = 2
     case (PSPIO_FMT_ABINIT_3)
       pspcod = 3
     case (PSPIO_FMT_ABINIT_4)
       pspcod = 4
     case (PSPIO_FMT_ABINIT_5)
       pspcod = 5
     case (PSPIO_FMT_ABINIT_6)
       pspcod = 6
     case (PSPIO_FMT_ABINIT_7)
       pspcod = 7
     case (PSPIO_FMT_ABINIT_8)
       pspcod = 8
     case (PSPIO_FMT_ABINIT_9)
       pspcod = 9
     case (PSPIO_FMT_ABINIT_10)
       pspcod = 10
     case (PSPIO_FMT_ABINIT_11)
       pspcod = 11
     case (PSPIO_FMT_ABINIT_17)
       pspcod = 17
   end select

!  ------------------------------------------------------------------------------
!  Check data for consistency against main routine input

!  Does required spin-orbit characteristics agree with format
   if ( psps%pspso(ipsp) /= 0 ) then
     write(msg, '(a,a,a,i3,a,a,a)' )&
&     'Libpspio cannot extract spin-orbit characteristics,',ch10,&
&     'while pspso(itypat)=',psps%pspso(ipsp),'.',ch10,&
&     'Action : contact ABINIT Group.'
     MSG_ERROR(msg)
   end if

!  Does nuclear charge znuclpsp agree with psp input znucl
   if ( abs(psps%znuclpsp(ipsp)-znucl) > tol8 ) then
     write(msg, '(a,f10.5,2a,f10.5,5a)' )&
&     'Pseudopotential file znucl=',znucl,ch10,&
&     'does not equal input znuclpsp=',psps%znuclpsp(ipsp),' better than 1e-08 .',ch10,&
&     'znucl is read from the psp file in pspatm_pspio, while',ch10,&
&     'znuclpsp is read in iofn2.'
     MSG_BUG(msg)
   end if

!  Is the highest angular momentum within limits?
!  Recall mpsang is 1+highest l for nonlocal correction.
!  Nonlocal corrections for s, p, d, and f are supported.
   if ( lmax >= psps%mpsang ) then
     write(msg, '(a,i10,a,i10,a,a)' )&
&     'input lmax+1=',lmax+1,' exceeds mpsang=',psps%mpsang,ch10,&
&     'indicates input lmax too large for dimensions.'
     MSG_BUG(msg)
   end if

!  Check several choices for ixc against pspxc
!  ixc is from ABINIT code; pspxc is from atomic psp file
   if ( dtset%ixc == 0 ) then
     msg = 'Note that input ixc=0 => no xc is being used.'
     MSG_WARNING(msg)
   else if ( dtset%ixc /= pspxc ) then
     write(msg,'(a,i0,3a,i0,a,a,a,a,a,a,a,a,a,a)' )&
&     'Pseudopotential file pspxc=',pspxc,',',ch10,&
&     'not equal to input ixc=',dtset%ixc,'.',ch10,&
&     'These parameters must agree to get the same xc ',ch10,&
&     'in ABINIT code as in psp construction.',ch10,&
&     'Action : check psp design or input file.',ch10,&
&     'Assume experienced user. Execution will continue.',ch10
     MSG_WARNING(msg)
   end if

   if (lloc>lmax .and. pspcod/=4 .and. pspcod/=8 .and. pspcod/=10) then
     write(msg, '(a,2i12,a,a,a,a)' )&
&     'lloc,lmax=',lloc,lmax,ch10,&
&     'chosen l of local psp exceeds range from input data.',ch10,&
&     'Action : check pseudopotential input file.'
     MSG_ERROR(msg)
   end if

!  Does the pspcod agree with type of calculation (paw or not)?
   if ( ((pspcod /= 7) .and. (pspcod /= 17) .and. (psps%usepaw == 1)) .or. &
&   ((pspcod == 7) .or. (pspcod == 17)) .and. (psps%usepaw == 0) ) then
     write(msg, '(a,i2,a,a,i1,a)' )&
&     'In reading atomic psp file, finds pspcod=',pspcod,ch10,&
&     'This is not an allowed value with usepaw=',psps%usepaw,'.'
     MSG_BUG(msg)
   end if

   if ( (.not. psps%vlspl_recipSpace) .and. (pspcod /= 2) .and. &
&   (pspcod /= 3) .and. (pspcod /= 10) .and. (pspcod /= 7) ) then
     write(msg, '(a,i2,a,a)' )&
&     'In reading atomic psp file, finds pspcod=',pspcod,ch10,&
&     'This is not an allowed value with real space computation.'
     MSG_BUG(msg)
   end if


!  -----------------------------------------------------------------------
!  Set various terms to 0 in case not defined below
   ABI_ALLOCATE(e990,(psps%mpssoang))
   ABI_ALLOCATE(e999,(psps%mpssoang))
   ABI_ALLOCATE(rcpsp,(psps%mpssoang))
   ABI_ALLOCATE(rms,(psps%mpssoang))
   ABI_ALLOCATE(epspsp,(psps%mpssoang))
   ABI_ALLOCATE(ekb1,(psps%mpssoang))
   ABI_ALLOCATE(ekb2,(psps%mpssoang))
   e990(:)   = zero
   e999(:)   = zero
   rcpsp(:)  = zero
   rms(:)    = zero
   ekb1(:)   = zero
   ekb2(:)   = zero
   epspsp(:) = zero
   qchrg     = 0

!  ----------------------------------------------------------------------
!  FIXME: extract psp data

   if ( pspcod == 6 ) then

     if (positron==1.and.abs(fchrg)<=tol14) then
       write(message,'(6a)')&
&       'You can only perform positronic ground-state calculations', &
&       ' (positron=1)',ch10,&
&       'using fhi pseudopotentials with a core density (fchrg>0)',ch10,&
&       'Action: change your psp file (add fchrg>0).'
       MSG_ERROR(message)
     end if

     !Will now proceed at the reading of pots and wfs

     !rad(:)=radial grid r(i)
     !vpspll(:,1),...,vpspll(:,4)=nonlocal pseudopotentials
     !wfll(:,1),...,wfll(:,4)=reference config. wavefunctions
     ABI_ALLOCATE(rad,(mmax))
     ABI_ALLOCATE(vpspll,(mmax,mpsang))
     ABI_ALLOCATE(wfll,(mmax,mpsang))

     !Read atomic pseudopotential for each l, filling up arrays vpspll
     !and wfll. Also set up rad array (actually read more than once)
     !Note: put each l into vpspll(:,l+1)
     do ipsang=1,lmax+1
       nproj(ipsang)=1
       read(tmp_unit,*)mmax2,amesh
       if(ipsang==1)then
         write(message, '(f10.6,t20,a)' ) amesh,' amesh (Hamman grid)'
         al_announced=log(amesh)
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')
       end if
       do irad=1,mmax
         read(tmp_unit,*)jj,rad(irad),wfll(irad,ipsang),vpspll(irad,ipsang)
         !    DEBUG
         !    Maybe the normalization is different
         !    wfll(irad,ipsang)=wfll(irad,ipsang)/rad(irad)
         !    ENDDEBUG
       end do
     end do


     !Generate core charge function and derivatives, if needed
     if(fchrg>tol14)then

       if (positron==1) then
         call psp6cc(mmax,n1xccc,rchrg,xccc1d,znucl,vh_tnzc=vpspll(:,lloc+1))
       else if(optnlxccc==1)then
         call psp6cc(mmax,n1xccc,rchrg,xccc1d,znucl)
       else if(optnlxccc==2)then
         call psp6cc_drh(mmax,n1xccc,rchrg,xccc1d)
       end if

       !The core charge function for pspcod=6
       !becomes zero beyond rchrg. Thus xcccrc must be set
       !equal to rchrg .
       xcccrc=rchrg
     else
       xccc1d(:,:)=zero
       xcccrc=zero
     end if

     !Compute in real(dp) al : the announced amesh is inaccurate.
     ratio=rad(mmax)/rad(1)
     al=log(ratio)/dble(mmax-1)

     !vloc(:)=Vlocal(r), lloc=0, 1, or 2 or -1 for avg.
     ABI_ALLOCATE(vloc,(mmax))
     !Copy appropriate nonlocal psp for use as local one
     vloc( 1:mmax ) = vpspll( 1:mmax , lloc+1 )

     !--------------------------------------------------------------------
     !Carry out calculations for local (lloc) pseudopotential.
     !Obtain Fourier transform (1-d sine transform)
     !to get q^2 V(q).

     call psp5lo(al,epsatm,mmax,mqgrid,qgrid,&
&     vlspl(:,1),rad,vloc,yp1,ypn,zion)

     !Fit spline to q^2 V(q) (Numerical Recipes subroutine)
     ABI_ALLOCATE(work_space,(mqgrid))
     ABI_ALLOCATE(work_spl,(mqgrid))
     call spline (qgrid,vlspl(:,1),mqgrid,yp1,ypn,work_spl)
     vlspl(:,2)=work_spl(:)
     ABI_DEALLOCATE(work_space)
     ABI_DEALLOCATE(work_spl)

     !--------------------------------------------------------------------
     !Take care of non-local part

     ABI_ALLOCATE(ekb_tmp,(mpsang))
     ABI_ALLOCATE(ffspl_tmp,(mqgrid,2,mpsang))

     !Zero out all Kleinman-Bylander energies to initialize
     ekb_tmp(:)=zero
     ekb(:)=zero

     !Allow for option of no nonlocal corrections (lloc=lmax=0)
     if (lloc==0.and.lmax==0) then
       write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     else

     !----------------------------------------------------------------------
     !Compute KB form factors and fit splines

       call psp5nl(al,ekb_tmp,ffspl_tmp,lmax,mmax,mpsang,mqgrid,qgrid,rad,&
&       vloc,vpspll,wfll)

     end if

     jj=0;index=0;indlmn(:,:)=0
     do ipsang=1,lmax+1
       !nproj had been set at 1, by default
       if(abs(ekb_tmp(ipsang))<tol10)then
         nproj(ipsang)=0
       end if
       !Possible values for nproj in this routine : 0 or 1.
       if(nproj(ipsang)==1)then
         if (useylm==1) then
           jj=jj+1
           do mm=1,2*ipsang-1
             index=index+1
             indlmn(1,index)=ipsang-1
             indlmn(2,index)=mm-ipsang
             indlmn(3,index)=1
             indlmn(4,index)=mm+(ipsang-1)*(ipsang-1)
             indlmn(5,index)=jj
             indlmn(6,index)=1
           end do
         else
           jj=jj+1
           index=index+1
           indlmn(1,index)=ipsang-1
           indlmn(2,index)=0
           indlmn(3,index)=1
           indlmn(4,index)=ipsang+(ipsang-1)*(ipsang-1)
           indlmn(5,index)=jj
           indlmn(6,index)=1
         end if
       end if
     end do
     !Transfer ekb and ffspl to their definitive location
     jpsang=1
     do ipsang=1,lmax+1
       if(nproj(ipsang)/=0)then
         ekb(jpsang)=ekb_tmp(ipsang)
         ffspl(:,:,jpsang)=ffspl_tmp(:,:,ipsang)
         jpsang=jpsang+1
         if(jpsang>lnmax)then
           write(message,'(3a,2i6)')&
&           'Problem with the dimension of the ekb and ffspl arrays.',ch10,&
&           'ipsang,lnmax=',ipsang,lnmax
         end if
       end if
     end do

     ABI_DEALLOCATE(ekb_tmp)
     ABI_DEALLOCATE(ffspl_tmp)
     ABI_DEALLOCATE(vpspll)
     ABI_DEALLOCATE(rad)
     ABI_DEALLOCATE(vloc)
     ABI_DEALLOCATE(wfll)

   end if

!----------------------------------------------------------------------
   if (pspcod==2 .or. pspcod==3 .or. pspcod==10)then
     write(msg, '(a,a,a,a,a,a,a,a,a,a)' )ch10,&
&     ' pspatm_pspio : COMMENT -',ch10,&
&     '  the projectors are not normalized,',ch10,&
&     '  so that the KB energies are not consistent with ',ch10,&
&     '  definition in PRB44, 8503 (1991). ',ch10,&
&     '  However, this does not influence the results obtained hereafter.'
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,  msg,'COLL')
   end if

   if (pspcod/=7.and.pspcod/=17) then
     write(msg, '(a,f14.8,a,a)' ) ' pspatm_pspio: epsatm=',epsatm,ch10,&
&     '         --- l  ekb(1:nproj) -->'

     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,  msg,'COLL')
     iln0=0
     do ilmn=1,psps%lmnmax
       iln=indlmn(5,ilmn)
       if (iln>iln0) then
         il=indlmn(1,ilmn)
         if (indlmn(6,ilmn)==1) then
           iln0=iln0+nproj(il+1)
           write(msg, '(13x,i1,4f12.6)' ) il,&
&           (ekb(iln+ii),ii=0,nproj(il+1)-1)
         else
           iln0=iln0+nproj(il+psps%mpsang)
           write(msg, '(2x,a,i1,4f12.6)' ) 'spin-orbit ',il,&
           (ekb(iln+ii),ii=0,nproj(il+psps%mpsang)-1)
         end if
         call wrtout(ab_out,msg,'COLL')
         call wrtout(std_out,  msg,'COLL')
       end if
     end do
   end if

   write(msg,'(3a)') ' pspatm_pspio: atomic psp has been read ',' and splines computed',ch10
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')

   ABI_DEALLOCATE(e990)
   ABI_DEALLOCATE(e999)
   ABI_DEALLOCATE(rcpsp)
   ABI_DEALLOCATE(rms)
   ABI_DEALLOCATE(ekb1)
   ABI_DEALLOCATE(ekb2)
   ABI_DEALLOCATE(epspsp)
   ABI_DEALLOCATE(nproj)

   if (dtset%prtvol > 9 .and. psps%usepaw==0 .and. psps%lmnmax>3) then
     write (filnam, '(a,i0,a)') trim(dtfil%fnameabo_pspdata), ipsp, ".dat"
     open (unit=dtfil%unt_pspdat, file=filnam)
     write (dtfil%unt_pspdat,*) '# Pseudopotential data in reciprocal space as used by ABINIT'
     write (dtfil%unt_pspdat,'(a)', ADVANCE='NO') '# index       vlocal   '
     if (psps%lnmax > 0) &
&     write (dtfil%unt_pspdat,'(a,I3)', ADVANCE='NO')   '           1st proj(l=', indlmn(1,1)
     if (psps%lnmax > 1) &
&     write (dtfil%unt_pspdat,'(a,I3)', ADVANCE='NO')   ')            2nd(l=', indlmn(1,2)
     if (psps%lnmax > 2) &
&     write (dtfil%unt_pspdat,'(a,I3,a)', ADVANCE='NO') ')            3rd(l=', indlmn(1,3), ')'
     write (dtfil%unt_pspdat,*)

     do ii = 1, psps%mqgrid_vl
       write(dtfil%unt_pspdat, '(I5,E24.16)', ADVANCE='NO') ii, vlspl(ii,1)
       if (psps%lnmax > 0) write(dtfil%unt_pspdat, '(E24.16)', ADVANCE='NO') ffspl(ii,1,1)
       if (psps%lnmax > 1) write(dtfil%unt_pspdat, '(E24.16)', ADVANCE='NO') ffspl(ii,1,2)
       if (psps%lnmax > 2) write(dtfil%unt_pspdat, '(E24.16)', ADVANCE='NO') ffspl(ii,1,3)
       write(dtfil%unt_pspdat, *)
     end do
     close(dtfil%unt_pspdat)

     write (filnam, '(a,i0,a)') trim(dtfil%fnameabo_nlcc_derivs), ipsp, ".dat"
     open (unit=dtfil%unt_nlcc_derivs, file=filnam)
     write (dtfil%unt_nlcc_derivs,*) '# Non-linear core corrections'
     write (dtfil%unt_nlcc_derivs,*) '#  r, pseudocharge, 1st, 2nd, 3rd, 4th, 5th derivatives'
     do ii = 1, psps%n1xccc
       write (dtfil%unt_nlcc_derivs,*) xcccrc*(ii-1)/(psps%n1xccc-1), xccc1d(ii,1), xccc1d(ii,2), &
&       xccc1d(ii,3), xccc1d(ii,4), &
&       xccc1d(ii,5), xccc1d(ii,6)
     end do
     close(dtfil%unt_nlcc_derivs)
   end if

 end if !me=0

 if (paral_mode==1) then
   call timab(48,1,tsec)
   call pawpsp_bcast(comm_mpi,epsatm,ffspl,pawrad,pawtab,vlspl,xcccrc)
   call timab(48,2,tsec)
 end if

#endif

 contains
!!***

!!****f* pspatm_pspio/pspio_check_error
!! NAME
!!  pspio_check_error
!!
!! FUNCTION
!!  Basic error handler for Libpspio.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2014 ABINIT group (YP)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ierr= error code
!!
!! OUTPUT
!!  Only printing.
!!
!! NOTES
!!  Stops Abinit upon error.
!! 
!! PARENTS
!!      pspatm_pspio
!!
!! CHILDREN
!!
!! SOURCE

  subroutine pspio_check_error(ierr)

    use defs_basis
#if defined HAVE_TRIO_LIBPSPIO
    use pspio_f90_lib_m
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspio_check_error'
!End of the abilint section

    implicit none

!Arguments ---------------------------------------------
!scalars
    integer, intent(in) :: ierr

!Local variables ---------------------------------------
!scalars
    character(len=512) :: msg
    integer :: itmp

! ******************************************************************************

! ******************************************************************************

#if defined HAVE_TRIO_LIBPSPIO

   if ( ierr /= PSPIO_SUCCESS ) then
     itmp = pspio_f90_error_flush()
     if ( itmp /= PSPIO_SUCCESS ) then
       write(msg, '(a,a,a,a)' )&
&       '  could not fetch error message from Libpspio.', ch10, ch10, &
&       '  Please report to Yann Pouillon <yann.pouillon@ehu.es>.'
       MSG_ERROR(msg)
     end if
     write(msg, '(a,a,a,a)' )&
&     '  possibly unsupported format.', ch10, ch10, &
&     '  Please check your pseudopotential file.'
     MSG_ERROR(msg)
   end if

#endif
 end subroutine pspio_check_error

end subroutine pspatm_pspio
!!***
