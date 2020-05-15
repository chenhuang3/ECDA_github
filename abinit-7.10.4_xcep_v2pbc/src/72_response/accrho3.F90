!{\src2tex{textfont=tt}}
!!****f* ABINIT/accrho3
!!
!! NAME
!! accrho3
!!
!! FUNCTION
!! Response function calculation only:
!!  Accumulate contribution to first-order density due do current (k,band)
!!  Also accumulate zero-order potential part of the 2nd-order total energy (if needed)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  counter=counter for status file
!!  cplex=1 if 1st-order density is real, 2 if 1st-order density is complex
!!  cwave0(2,npw*nspinor)=GS wavefunction at k, in reciprocal space
!!  cwave1(2,npw1*nspinor)=1st-order wavefunction at k,q, in reciprocal space
!!  cwavef(2,npw1*nspinor)=1st-order wavefunction at k,q, in reciprocal space, without correction due to occupation change
!!  cwaveprj0(natom,nspinor*usecprj)= GS wave function at k projected with nl projectors
!!  cwaveprj1(natom,nspinor*usecprj)= 1st-order wave function at k,q projected with nl projectors
!!  dimcprj(natom*usepaw)=array of dimensions of arrays cprj, cprjq (ordered by atom-type)
!!  dimffnlk=second dimension of ffnlk (1+number of derivatives)
!!  dimphkxred=second dimension of phkxred
!!  ffnlk(npw_k,dimffnlk,lmnmax,1+usepaw(ntypat-1))=nonloc form factors at k, for the displaced atom.
!!  filstat=name of the status file
!!  gbound(2*mgfft+8,2)=G sphere boundary
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  iband=index of current band
!!  idir=direction of the current perturbation
!!  ipert=type of the perturbation
!!  isppol=1 index of current spin component
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kpg_k(npw_k,nkpg)= (k+G) components at k (only if useylm=1)
!!  kpt(3)=reduced coordinates of k points.
!!  kptopt=option for the generation of k points
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband_k=number of bands at this k point for that spin polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
!!  npw_k=number of planewaves in basis sphere at k
!!  npw1_k=number of planewaves in basis sphere at k+q
!!  nspinor=number of spinorial components of the wavefunctions
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  option= 1: accumulate 1st-order density,
!!          2: accumulate 0-order potential part of the 2nd-order total energy
!!          3: accumulate both
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  phkxred(2,dimphkxred)=phase factors exp(2 pi kpoint.xred) at k
!!  prtvol=control print volume and debugging output
!!  tim_fourwf= timing code for fourwf (5 from vtowfk3, 18 from nstwf3)
!!  usecprj= 1 if cwaveprj0 array is stored in memory
!!  vlocal(n4,n5,n6)= GS local potential in real space, on the augmented fft grid
!!  wf_corrected=flag put to 1 if cwave1 is different from cwavef (if there is a contribution from occ. change)
!!  wtk_k=weight assigned to the k point.
!!
!! OUTPUT
!!  ====== if option=2 or option=3 =====
!!  eloc0_k=zero-order local contribution to 2nd-order total energy for current band and k
!!
!! SIDE EFFECTS
!!  ====== if option=1 or option=3 =====
!!    rhoaug1(cplex*n4,n5,n6)= density in electrons/bohr**3,
!!    ==== if gs_hamkq%usepaw=1 =====
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!                                            (cumulative, so input as well as output)
!!
!! NOTES
!!  In this part of the treatment of one band, one has to
!!  perform Fourier transforms, and to treat separately the
!!  two spinorial components of the wavefunction.
!!  Was part of vtowfk3 before.
!!
!! PARENTS
!!      nstpaw3,vtowfk3,wfkfermi3
!!
!! CHILDREN
!!      fourwf,get_my_atmtab,getcprj,pawaccrhoij,pawcprj_alloc,pawcprj_destroy
!!      status
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine accrho3(counter,cplex,cwave0,cwave1,cwavef,cwaveprj0,cwaveprj1,&
&                  dimcprj,dimffnlk,ffnlk,dimphkxred,eloc0_k,filstat,gbound,&
&                  gs_hamkq,iband,idir,ipert,isppol,kg_k,kg1_k,kpg_k,kpt,kptopt,&
&                  mpi_enreg,natom,nband_k,ncpgr,nkpg,npw_k,npw1_k,nspinor,occ_k,&
&                  option,paral_kgb,pawrhoij1,ph3d,phkxred,prtvol,rhoaug1,tim_fourwf,&
&                  usecprj,vlocal,wf_corrected,wtk_k,mpi_comm_atom,mpi_atmtab)


 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_pawrhoij,      only : pawrhoij_type
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_destroy
 use m_paral_atom,    only : get_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'accrho3'
 use interfaces_32_util
 use interfaces_53_ffts
 use interfaces_65_nonlocal
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: counter,cplex,dimffnlk,dimphkxred,iband,idir,ipert,isppol,kptopt,natom,nband_k
 integer,intent(in) :: ncpgr,nkpg,npw_k,npw1_k,nspinor,option,paral_kgb,prtvol,tim_fourwf
 integer,intent(in) :: usecprj,wf_corrected
 integer,optional,intent(in) :: mpi_comm_atom
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: wtk_k
 real(dp),intent(out) :: eloc0_k
 character(len=fnlen),intent(in) :: filstat
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: dimcprj(natom*gs_hamkq%usepaw),gbound(2*gs_hamkq%mgfft+8,2)
 integer,intent(in) :: kg1_k(3,npw1_k),kg_k(3,npw_k)
 real(dp),intent(in),target :: cwave0(2,npw_k*nspinor),cwave1(2,npw1_k*nspinor),cwavef(2,npw1_k*nspinor)
 real(dp),intent(in) :: ffnlk(npw_k,dimffnlk,gs_hamkq%lmnmax,1+gs_hamkq%usepaw*(gs_hamkq%ntypat-1))
 real(dp),intent(in) :: kpg_k(npw_k,nkpg),kpt(3),occ_k(nband_k),phkxred(2,dimphkxred)
 real(dp),intent(in) :: vlocal(gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6*(option/2))
 real(dp),intent(inout) :: ph3d(2,npw1_k,gs_hamkq%matblk)
 real(dp),intent(inout) :: rhoaug1(cplex*gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6)
 type(pawcprj_type),intent(in) :: cwaveprj0(natom,nspinor*usecprj),cwaveprj1(natom,nspinor*gs_hamkq%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=14
 integer :: choice,cplex_cprj,i1,i2,i3,iexit,ispinor,my_comm_atom,my_natom,n1,n2,n3,option_rhoij
 logical :: my_atmtab_allocated,paral_atom
 logical :: usetimerev
 real(dp) :: im0,im1,re0,re1,valuer,weight
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: dummy(2,1)
 real(dp),allocatable :: rhoaug(:,:,:),wfraug(:,:,:,:),wfraug1(:,:,:,:)
 real(dp),pointer :: cwavef_sp(:,:)
 type(pawcprj_type),allocatable :: cwaveprj_tmp(:,:)


! *********************************************************************

 DBG_ENTER("COLL")

 if (option/=1.and.option/=2.and.option/=3) return

!Initializations
 ABI_ALLOCATE(rhoaug,(gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
 ABI_ALLOCATE(wfraug1,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
 n1=gs_hamkq%ngfft(1);n2=gs_hamkq%ngfft(2);n3=gs_hamkq%ngfft(3)
 if (option==2.or.option==3) eloc0_k=zero

!Loop on spinorial components
 do ispinor=1,nspinor

   if (prtvol>=10) then
     call status(counter,filstat,iexit,level,'density update')
   end if

!  Part devoted to the accumulation of the 0-order potential part of the 2nd-order total energy
!  --------------------------------------------------------------------------------------------

!  Fourier transform of cwavef. Here, rhoaug1 is a dummy variable.
!  NOTE : should take into account nspinor
   if (wf_corrected==0.or.option==2.or.option==3) then
     if (ispinor==1) then
       cwavef_sp => cwavef(:,1:npw1_k)
     else
       cwavef_sp => cwavef(:,1+npw1_k:2*npw1_k)
     end if
     call fourwf(cplex,rhoaug,cwavef_sp,dummy,wfraug1,gs_hamkq%gbound,gs_hamkq%gbound,&
&     gs_hamkq%istwf_k,kg1_k,kg1_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,npw1_k,&
&     1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,paral_kgb,tim_fourwf,weight,weight,&
&     use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     nullify(cwavef_sp)
   end if

!  Compute contribution of this band to zero-order potential part of the 2nd-order total energy
   if (option==2.or.option==3) then
     valuer=zero
     do i3=1,n3
       do i2=1,n2
         do i1=1,n1
           valuer=valuer+vlocal(i1,i2,i3)*(wfraug1(1,i1,i2,i3)**2+wfraug1(2,i1,i2,i3)**2)
         end do
       end do
     end do

!    Local potential energy of this band
     eloc0_k=eloc0_k+two*valuer/dble(gs_hamkq%nfft)
   end if ! option

!  Part devoted to the accumulation of the 1st-order density
!  ---------------------------------------------------------
   if (option==1.or.option==3) then

!    Compute 1st-order WF in real space
!    One needs the Fourier transform of cwave1. However, only the one of
!    cwavef is available. If cwavef and cwave1 differs, this Fourier
!    transform must be computed. In both case the result is in wfraug1.
     if (wf_corrected==1) then
       if (ispinor==1) then
         cwavef_sp => cwave1(:,1:npw1_k)
       else
         cwavef_sp => cwave1(:,1+npw1_k:2*npw1_k)
       end if
       call fourwf(cplex,rhoaug,cwavef_sp,dummy,wfraug1,gs_hamkq%gbound,gs_hamkq%gbound,&
&       gs_hamkq%istwf_k,kg1_k,kg1_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,npw1_k,&
&       1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,paral_kgb,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
       nullify(cwavef_sp)
     end if

!    Compute 0-order WF in real space
     ABI_ALLOCATE(wfraug,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
     if (ispinor==1) then
       cwavef_sp => cwave0(:,1:npw_k)
     else
       cwavef_sp => cwave0(:,1+npw_k:2*npw_k)
     end if
     call fourwf(1,rhoaug,cwavef_sp,dummy,wfraug,gbound,gbound,gs_hamkq%istwf_k,kg_k,kg_k,&
&     gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,npw_k,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,&
&     0,paral_kgb,tim_fourwf,weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     nullify(cwavef_sp)

!    The factor 2 is not the spin factor (see Eq.44 of PRB55,10337 (1997))
     weight=two*occ_k(iband)*wtk_k/gs_hamkq%ucvol

!    Accumulate 1st-order density
     if (cplex==2) then
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             re0=wfraug(1,i1,i2,i3)  ; im0=wfraug(2,i1,i2,i3)
             re1=wfraug1(1,i1,i2,i3) ; im1=wfraug1(2,i1,i2,i3)
             rhoaug1(2*i1-1,i2,i3)=rhoaug1(2*i1-1,i2,i3)+weight*(re0*re1+im0*im1)
             rhoaug1(2*i1  ,i2,i3)=rhoaug1(2*i1  ,i2,i3)+weight*(re0*im1-im0*re1)
           end do
         end do
       end do
     else
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             rhoaug1(i1,i2,i3)=rhoaug1(i1,i2,i3) &
&             +weight*(wfraug(1,i1,i2,i3)*wfraug1(1,i1,i2,i3) &
&             +wfraug(2,i1,i2,i3)*wfraug1(2,i1,i2,i3))
           end do
         end do
       end do
     end if
     ABI_DEALLOCATE(wfraug)
   end if ! option

 end do ! Loop on spinorial components

 ABI_DEALLOCATE(rhoaug)
 ABI_DEALLOCATE(wfraug1)

!Part devoted to the accumulation of the 1st-order occupation matrix in PAW case
!-------------------------------------------------------------------------------


 if ((option==1.or.option==3).and.gs_hamkq%usepaw==1) then

!  Set up parallelism over atoms
   my_natom=natom; if(gs_hamkq%usepaw==1) my_natom=size(pawrhoij1)
   paral_atom=(present(mpi_comm_atom).and.(my_natom/=natom))
   my_comm_atom=xmpi_self;if (present(mpi_comm_atom)) my_comm_atom=mpi_comm_atom
   nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
   call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

   cplex_cprj=2;if (gs_hamkq%istwf_k>1) cplex_cprj=1
   option_rhoij=2;usetimerev=(kptopt>0.and.kptopt<3)

   if (usecprj==1) then
     call pawaccrhoij(gs_hamkq%atindx,cplex_cprj,cwaveprj0,cwaveprj1,ipert,isppol,&
&     my_natom,natom,nspinor,occ_k(iband),option_rhoij,pawrhoij1,usetimerev,wtk_k,&
&     mpi_comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
   else
     ABI_DATATYPE_ALLOCATE(cwaveprj_tmp,(natom,nspinor))
     call pawcprj_alloc(cwaveprj_tmp,ncpgr,dimcprj)
     choice=2
     call getcprj(choice,0,cwave0,cwaveprj_tmp,gs_hamkq%dimekb1,gs_hamkq%dimekb2,&
&     dimffnlk,gs_hamkq%ekb(:,:,1),ffnlk,idir,gs_hamkq%indlmn,gs_hamkq%istwf_k,&
&     kg_k,kpg_k,kpt,gs_hamkq%lmnmax,gs_hamkq%matblk,gs_hamkq%mgfft,mpi_enreg,natom,&
&     gs_hamkq%nattyp,gs_hamkq%ngfft,nkpg,gs_hamkq%nloalg,npw_k,nspinor,gs_hamkq%ntypat,phkxred,&
&     gs_hamkq%ph1d,ph3d,gs_hamkq%ucvol,gs_hamkq%usepaw,gs_hamkq%useylm)
     call pawaccrhoij(gs_hamkq%atindx,cplex_cprj,cwaveprj_tmp,cwaveprj1,ipert,isppol,&
&     my_natom,natom,nspinor,occ_k(iband),option_rhoij,pawrhoij1,usetimerev,wtk_k, &
&     mpi_comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
     call pawcprj_destroy(cwaveprj_tmp)
     ABI_DATATYPE_DEALLOCATE(cwaveprj_tmp)
   end if

 end if

 DBG_EXIT("COLL")

end subroutine accrho3
!!***
