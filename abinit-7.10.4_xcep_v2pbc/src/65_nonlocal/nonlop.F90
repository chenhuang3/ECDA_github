!{\src2tex{textfont=tt}}
!!****f* ABINIT/nonlop
!! NAME
!! nonlop
!!
!! FUNCTION
!! This routine is a driver to compute:
!! * Application of a nonlocal operator Vnl in order to get:
!!    - contracted elements (energy, forces, stresses, ...), if signs=1
!!    - a function in reciprocal space (|out> = Vnl|in>), if signs=2
!! * Optionally, in case of PAW calculation:
!!   - Application of the overlap matrix in reciprocal space
!!     (<in|S|in> or (I+S)|in>).
!!   - Application of (Vnl-lambda.S) in reciprocal space
!!     (<in|Vnl-lambda.S|in> and derivatives or (Vnl-lambda.S)|in>).
!! According to user s choice, the routine calls a subroutine, computing all quantities:
!!   - using Legendre Polynomials Pl (Norm-conserving psps only)
!!   - using Spherical Harmonics Ylm (N-conserving or PAW ; compulsory for PAW)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  choice: chooses possible output:
!!    choice=0 => do nothing (only compute WF projected with NL projectors)
!!          =1 => non-local energy contribution
!!          =2 => 1st derivative(s) with respect to atomic position(s)
!!          =3 => 1st derivative(s) with respect to strain(s)
!!          =23=> 1st derivative(s) with respect to atm. pos. and strain(s)
!!          =4 => 2nd derivative(s) with respect to 2 atomic pos.
!!          =24=> 1st derivative(s) with respect to atm. pos. and
!!                2nd derivative(s) with respect to 2 atomic pos.
!!          =5 => 1st derivative(s) with respect to k wavevector, typically
!!                $\sum_{ij}|p_i\rangle D_{ij}\langle dp_j/dk| + |dp_i/dk\rangle D_{ij}\langle p_j|$
!!          =51 =>right 1st derivative(s) with respect to k wavevector, typically
!!                $\sum_{ij}|p_i\rangle D_{ij}\langle dp_j/dk|$
!!          =52 =>left 1st derivative(s) with respect to k wavevector, typically
!!                $\sum_{ij}|dp_i/dk\rangle D_{ij}\langle p_j|$
!!          =53 =>twist 1st derivative(s) with respect to k, typically
!!                $\sum_{ij}|dp_i/dk_(idir+1)\rangle D_{ij}\langle dp_j/dk_(idir-1)| -
!!                 |dp_i/dk_(idir-1)\rangle D_{ij}\langle dp_j/dk_(idir+1)|$
!!          =54=> mixed 2nd derivative(s) with respect to atomic pos. and left k wavevector
!!          =7 => apply operator $\sum_{i}|p_i\rangle \langle p_i|$,
!!                same as overlap operator with s_ij=identity (paw_opt==3 only)
!!          =8 => 2nd derivatives with respect to 2 k wavevectors
!!    Only choices 1,2,3,23,4,5,6 are compatible with useylm=0.
!!    Only choices 1,2,3,5,7,51,52 and 53 are compatible with signs=2
!!    Choices 51,52 are not compatible with signs=1
!!  cpopt=flag defining the status of cprjin%cp(:)=<Proj_i|Cnk> scalars (see below, side effects)
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  dimffnlin=second dimension of ffnlin (1+number of derivatives)
!!  dimffnlout=second dimension of ffnlout (1+number of derivatives)
!!  enl(dimenl1,dimenl2,nspinortot**2)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!  ->PAW :             ==== when paw_opt=1 or 4 ====
!!                      (Real, symmetric) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!                      These are complex numbers if cplex_enl=2
!!                        enl(:,:,1) contains Dij^up-up
!!                        enl(:,:,2) contains Dij^dn-dn
!!                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!  ffnlin(npwin,dimffnlin,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |in> vector
!!  ffnlout(npwout,dimffnlout,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |out> vector
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,51,or 52 and signs=2)
!!                          for choice 53 signs=2, cross derivatives are in idir-1 and idir+1 directions
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                                  or i=lmn (if useylm=1)
!!  istwf_k=option parameter that describes the storage of wfs
!!  kgin(3,npwin)=integer coords of planewaves in basis sphere, for the |in> vector
!!  kgout(3,npwout)=integer coords of planewaves in basis sphere, for the |out> vector
!!  kpgin(npw,npkgin)= (k+G) components and related data, for the |in> vector  (only if useylm=1)
!!  kpgout(npw,nkpgout)=(k+G) components and related data, for the |out> vector (only if useylm=1)
!!  kptin(3)=k point in terms of recip. translations, for the |in> vector
!!  kptout(3)=k point in terms of recip. translations, for the |out> vector
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!         Typically lambda is the eigenvalue (or its guess)
!!  lmnmax=max. number of (l,m,n) components over all types of atoms
!!  matblk=dimension of the arrays ph3din and ph3dout
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+max(spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  ndat=number of wavefunctions on which to apply nonlop
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpgin,nkpgout=second sizes of arrays kpgin/kpgout
!!  nloalg(5)=governs the choice of the algorithm for nonlocal operator
!!  nnlout=dimension of enlout (when signs=1 and choice>0):
!!         ==== if paw_opt=0, 1 or 2 ====
!!         choice   nnlout     |  choice   nnlout
!!              1   1          |       5   3
!!              2   3*natom    |      53   3
!!              3   6          |      54   9*natom
!!              4   6*natom    |       6   36+18*natom
!!             23   6+3*natom  |       8   6
!!             24   9*natom    |
!!         ==== if paw_opt=3 ====
!!         choice   nnlout
!!              1   1
!!              2   3*natom
!!             54   9*natom
!!              7   1
!!         ==== if paw_opt=4 ====
!!         not available
!!  npwin=number of planewaves for given k point, for the |in> vector
!!  npwout=number of planewaves for given k point, for the |out> vector
!!  nspinor=total number of spinorial components of the wavefunctions
!!  nspinortot=number of spinorial components of the wavefunctions on current proc
!!  ntypat=number of types of atoms in cell
!!  only_SO=flag to calculate only the SO part in nonlop
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  phkxredin(2,natom)=phase factors exp(2 pi kptin.xred)
!!  phkxredout(2,natom)=phase factors exp(2 pi kptout.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
!!  ph3din(2,npwin,matblk)=3D structure factors, for each atom and plane wave (in)
!!  ph3dout(2,npwout,matblk)=3-dim structure factors, for each atom and plane wave (out)
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  sij(dimenl1,ntypat)=overlap matrix components (only if paw_opt=2, 3 or 4)
!!  tim_nonlop=timing code of the calling routine (can be set to 0 if not attributed)
!!  ucvol=unit cell volume (bohr^3)
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  vectin(2,npwin*nspinor*ndat)=input cmplx wavefunction coefficients <G|Cnk>
!!
!! OUTPUT
!! ==== if (signs==1) ====
!! --If (paw_opt==0, 1 or 2)
!!    enlout(nnlout)= contribution to the non-local part of the following properties:
!!      if choice=1 : enlout(1)             -> the energy
!!      if choice=2 : enlout(3*natom)       -> 1st deriv. of energy wrt atm. pos (forces)
!!      if choice=3 : enlout(6)             -> 1st deriv. of energy wrt strain (stresses)
!!      if choice=4 : enlout(6*natom)       -> 2nd deriv. of energy wrt 2 atm. pos (dyn. mat.)
!!      if choice=23: enlout(6+3*natom)     -> 1st deriv. of energy wrt atm. pos (forces) and
!!                                             1st deriv. of energy wrt strain (stresses)
!!      if choice=24: enlout(9*natom)       -> 1st deriv. of energy wrt atm. pos (forces) and
!!                                             2nd deriv. of energy wrt 2 atm. pos (dyn. mat.)
!!      if choice=5 : enlout(3)             -> 1st deriv. of energy wrt k
!!      if choice=53: enlout(3)             -> 1st deriv. (twist) of energy wrt k
!!      if choice=54: enlout(9*natom)       -> 2nd deriv. of energy wrt atm. pos and left  k (Born eff. charge)
!!      if choice=6 : enlout(36+18*natom)   -> 2nd deriv. of energy wrt 2 strains (elast. tensor) and
!!                                             2nd deriv. of energy wrt to atm. pos and strain (internal strain)
!!      if choice=8 : enlout(6)             -> 2nd deriv. of energy wrt 2 k
!! --If (paw_opt==3)
!!      if choice=1 : enlout(1)             -> contribution to <c|S|c> (note: not including <c|c>)
!!      if choice=2 : enlout(3*natom)       -> contribution to <c|dS/d_atm.pos|c>
!!      if choice=54: enlout(9*natom)       -> 2nd deriv. of energy wrt atm. pos and left  k (Born eff. charge)
!!      if choice=7 : enlout(1)             -> contribution to <c|sum_i[p_i><p_i]|c>
!! --If (paw_opt==4)
!!      not available
!! ==== if (signs==2) ====
!! --if (paw_opt=0, 1 or 4)
!!    vectout(2,npwout*nspinor*ndat)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.:
!!      if (choice=1) <G|V_nonlocal|vect_in>
!!      if (choice=2) <G|dV_nonlocal/d(atm. pos)|vect_in>
!!      if (choice=3) <G|dV_nonlocal/d(strain)|vect_in>
!!      if (choice=5) <G|dV_nonlocal/d(k)|vect_in>
!!      if (choice=51) <G|d(right)V_nonlocal/d(k)|vect_in>
!!      if (choice=52) <G|d(left)V_nonlocal/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)V_nonlocal/d(k)|vect_in>
!!  if (paw_opt=2)
!!    vectout(2,npwout*nspinor*ndat)=final vector in reciprocal space:
!!      if (choice=1) <G|V_nonlocal-lamdba.(I+S)|vect_in> (note: not including <G|I|c>)
!!      if (choice=2) <G|d[V_nonlocal-lamdba.(I+S)]/d(atm. pos)|vect_in>
!!      if (choice=3) <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_in>
!!      if (choice=5) <G|d[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=51) <G|d(right)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=52) <G|d(left)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!! --if (paw_opt=3 or 4)
!!    svectout(2,npwout*nspinor*ndat)=result of the aplication of Sij (overlap matrix)
!!                  or one of its derivatives to the input vect.:
!!      if (choice=1) <G|I+S|vect_in> (note: not including <G|I|c>)
!!      if (choice=2) <G|dS/d(atm. pos)|vect_in>
!!      if (choice=3) <G|dS/d(strain)|vect_in>
!!      if (choice=5) <G|dS/d(k)|vect_in>
!!      if (choice=51) <G|d(right)S/d(k)|vect_in>
!!      if (choice=52) <G|d(left)S/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)S/d(k)|vect_in>
!!      if (choice=3) <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_in>
!!      if (choice=7) <G|sum_i[p_i><p_i]|vect_in>
!!
!! SIDE EFFECTS
!!  ==== ONLY IF useylm=1
!!  cprjin(natom,nspinor*ndat) <type(pawcprj_type)>=projected input wave function |in> on non-local projectors
!!                                  =<p_lmn|in> and derivatives
!!                    Treatment depends on cpopt parameter:
!!                     if cpopt=-1, <p_lmn|in> (and derivatives)
!!                                  are computed here (and not saved)
!!                     if cpopt= 0, <p_lmn|in> are computed here and saved
!!                                  derivatives are eventually computed but not saved
!!                     if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 2  <p_lmn|in> are already in memory;
!!                                  first (and 2nd) derivatives are computed here and not saved
!!                     if cpopt= 3  <p_lmn|in> are already in memory;
!!                                  first derivatives are computed here and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!!                                  other derivatives are not computed
!!                                  This option is not compatible with choice=4,24 or 6
!!                     (if useylm=0, should have cpopt=-1)
!!
!! NOTES
!! * See nonlop_pl and nonlop_ylm to have more comments...
!! * In the case signs=1, the array vectout is not used, nor modified
!!   so that the same array as vectin can be used as a dummy argument;
!!   the same is true for the pairs npwin-npwout, ffnlin-ffnlout,
!!   kgin-kgout, ph3din-ph3dout, phkredin-phkxredout).
!!
!! PARENTS
!!      d2frnl,eltfrnl3,energy,forstrnps,getgh1c,getghc,getgsc,ladielmt,lavnl
!!      lobpcgIIwf,lobpcgccIIIwf,lobpcgccIIwf,m_invovl,m_lobpcgIIIwf
!!      make_grad_berry,nstwf4,prctfvw1,prctfvw2,prep_nonlop,resp3dte,vtowfk
!!      wfd_vnlpsi
!!
!! CHILDREN
!!      gemm_nonlop,nonlop_gpu,nonlop_pl,nonlop_ylm,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine nonlop(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,&
&                 use_gpu_cuda) !optional argument

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_gemm_nonlop

 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nonlop'
 use interfaces_18_timing
 use interfaces_65_nonlocal, except_this_one => nonlop
!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimffnlin,dimffnlout,idir
 integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
 integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
 integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm
 integer,optional,intent(in) :: use_gpu_cuda
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 !arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kgin(3,npwin)
 integer,intent(in) :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(5)
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
 real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat),gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpgin(npwin,nkpgin*useylm)
 real(dp),intent(in) :: kpgout(npwout,nkpgout*useylm),kptin(3),kptout(3)
 real(dp),intent(in) :: lambda(ndat)
 !real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),phkxredin(2,natom) !vz_d
 real(dp),intent(in) :: phkxredin(2,natom) !vz_d
 real(dp),intent(in) :: ph1d(2,*) !vz_d
 real(dp),intent(in) :: phkxredout(2,natom),sij(dimenl1,ntypat*((paw_opt+1)/3))
 real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
 real(dp),intent(inout),target :: vectin(2,npwin*nspinor*ndat)
 real(dp),intent(out),target :: enlout(nnlout*ndat)
 real(dp),intent(out),target :: svectout(2,npwout*nspinor*(paw_opt/3)*ndat)
 real(dp),intent(inout),target :: vectout(2,npwout*nspinor*ndat) !vz_i
 type(pawcprj_type),intent(inout),target :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

 !Local variables-------------------------------
 !scalars
 integer :: use_gpu_cuda_
 logical :: use_gemm_nonlop
 character(len=500) :: message
 !arrays
 real(dp) :: tsec(2)
 !idat pointers
 integer :: idat
 type(pawcprj_type),pointer :: cprjin_idat(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: vectin_idat(:,:), vectout_idat(:,:), svectout_idat(:,:), enlout_idat(:)

 DBG_ENTER("COLL")

! **********************************************************************

!Keep track of time spent in this routine. Note the selection of
!different slots for different choices.
 call timab(220+tim_nonlop,1,tsec)

 use_gpu_cuda_=0;if (present(use_gpu_cuda)) use_gpu_cuda_=use_gpu_cuda

 !Error on bad input
 if (useylm==0) then
   if(paw_opt>0) then
     message = ' when paw_opt>0 you must use ylm version of nonlop! Set useylm 1'
     MSG_BUG(message)
   end if
   if(cpopt/=-1) then
     message = 'if useylm=0, ie no PAW, then cpopt/=-1 is not allowed !'
     MSG_BUG(message)
   end if
   if(use_gpu_cuda_/=0) then
     message = 'when use_gpu_cuda/=0 you must use ylm version of nonlop! Set useylm 1'
     MSG_BUG(message)
   end if
 end if

!A specific version of nonlop based on BLAS3 can be used
!But there are several restrictions
 use_gemm_nonlop= ( &
& gemm_nonlop_use_gemm &
& .and. signs == 2 .and. paw_opt /= 2 .and. &
& (choice < 2 .or. choice == 7) .and. nspinor == 1 &
& .and. cpopt < 3 .and. useylm /= 0 )

 if(use_gemm_nonlop) then
   call gemm_nonlop(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimffnlin,dimffnlout,&
&   enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&   kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&   mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&   nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&   phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&   tim_nonlop,ucvol,useylm,vectin,vectout,use_gpu_cuda)
 else
   do idat=1, ndat
     vectin_idat => vectin(:,1+npwin*nspinor*(idat-1):npwin*nspinor*idat)     
     if(choice /= 0) then
       vectout_idat => vectout(:,1+npwout*nspinor*(idat-1):npwout*nspinor*idat)
     else
       vectout_idat => vectout
     end if    
     if(paw_opt >= 3) then
       svectout_idat => svectout(:,1+npwout*nspinor*(idat-1):npwout*nspinor*idat)
     else
       svectout_idat => svectout
     end if     
     if(cpopt >= 0) then
       cprjin_idat => cprjin(:, nspinor*(idat-1)+1:nspinor*(idat))
     else
       cprjin_idat => cprjin
     end if    
     if (nnlout>0) then
       enlout_idat => enlout((idat-1)*nnlout+1:(idat*nnlout))
     else
       enlout_idat => enlout
     end if
     if (useylm==0) then
       call nonlop_pl (choice,dimenl1,dimenl2,dimffnlin,dimffnlout,&
&       enl,enlout_idat,ffnlin,ffnlout,gmet,gprimd,idir,&
&       indlmn,istwf_k,kgin,kgout,kpgin,kpgout,kptin,kptout,lmnmax,matblk,&
&       mgfft,mpi_enreg,mpsang,mpssoang,natom,nattyp,ngfft,nkpgin,nkpgout,&
&       nloalg,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,phkxredin,&
&       phkxredout,ph1d,ph3din,ph3dout,signs,ucvol,vectin_idat,vectout_idat)
     else if (use_gpu_cuda_==0) then
       call nonlop_ylm(atindx1,choice,cpopt,cprjin_idat,dimenl1,dimenl2,dimffnlin,dimffnlout,&
&       enl,enlout_idat,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&
&       kgin,kgout,kpgin,kpgout,kptin,kptout,lambda(idat),lmnmax,matblk,mgfft,&
&       mpi_enreg,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&
&       npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&
&       ph3din,ph3dout,signs,sij,svectout_idat,ucvol,vectin_idat,vectout_idat)
     else
       call nonlop_gpu(atindx1,choice,cpopt,cprjin_idat,dimenl1,dimenl2,dimffnlin,dimffnlout,&
&       enl,enlout_idat,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&
&       kgin,kgout,kpgin,kpgout,kptin,kptout,lambda(idat),lmnmax,matblk,mgfft,&
&       mpi_enreg,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&
&       npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&
&       ph3din,ph3dout,signs,sij,svectout_idat,ucvol,vectin_idat,vectout_idat)
     end if
   end do
 end if

 call timab(220+tim_nonlop,2,tsec)

 DBG_EXIT("COLL")

end subroutine nonlop
!!***
