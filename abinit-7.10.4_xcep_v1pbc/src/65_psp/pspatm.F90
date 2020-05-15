!{\src2tex{textfont=tt}}
!!****f* ABINIT/pspatm
!! NAME
!!  pspatm
!!
!! FUNCTION
!!  Wrapper routine to read pseudopotentials using Libpspio or Abinit.
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
!!      pspini
!!
!! CHILDREN
!!      pspatm_abinit,pspatm_pspio
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pspatm(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&
&  psps,vlspl,dvlspl,xcccrc,xccc1d,comm_mpi)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_pawrad, only : pawrad_type
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspatm'
 use interfaces_65_psp, except_this_one => pspatm
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipsp
 integer, optional,intent(in) :: comm_mpi
 real(dp),intent(in) :: dq
 real(dp),intent(out) :: epsatm,xcccrc
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(pawrad_type),intent(inout) :: pawrad !vz_i
 type(pawtab_type),intent(inout) :: pawtab !vz_i
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(inout) :: indlmn(6,psps%lmnmax) !vz_i
 real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
 real(dp),intent(inout) :: ekb(psps%dimekb*(1-psps%usepaw))   !vz_i
 real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax) !vz_i
 real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
 real(dp),intent(inout) :: xccc1d(psps%n1xccc*(1-psps%usepaw),6) !vz_i

!Local variables ---------------------------------------

! ******************************************************************************

#if defined HAVE_TRIO_LIBPSPIO
 if ( present(comm_mpi) ) then
   call pspatm_pspio(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&
&   psps,vlspl,dvlspl,xcccrc,xccc1d,comm_mpi)
 else
   call pspatm_pspio(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&
&   psps,vlspl,dvlspl,xcccrc,xccc1d)
 end if
#else
 if ( present(comm_mpi) ) then
   call pspatm_abinit(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&
&   psps,vlspl,dvlspl,xcccrc,xccc1d,comm_mpi)
 else
   call pspatm_abinit(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&
&   psps,vlspl,dvlspl,xcccrc,xccc1d)
 end if
#endif

end subroutine pspatm
!!***
