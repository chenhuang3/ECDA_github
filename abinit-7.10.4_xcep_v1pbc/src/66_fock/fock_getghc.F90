!{\src2tex{textfont=tt}}
!!****f* ABINIT/fock_getghc
!! NAME
!!  fock_getghc
!!
!! FUNCTION
!!  Compute the matrix elements <G|Vx|psi> of the Fock operator.
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (CMartins)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cwavef(2,npw*nspinor*ndat)= planewave coefficients of wavefunctions on which Fock operator is applied.
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  gbound(2*mgfft+8,2)= G sphere boundary for zero-padded FFT.
!!  gmet(3,3)= reciprocal space metric tensor in Bohr**-2
!!  istwf_k= option parameter that describes the storage of the input wavefunctions (k-dependent)
!!  kg_k(3,npw)= reduced coordinates of the G-vectors.
!!  mgfft= maximum size for 1D FFTs 
!!  mpi_enreg= information about MPI parallelization
!!  n4,n5,n6= same as ngfft(4:6)
!!  nfft= number of FFT grid points.
!!  ngfft(18)=contain all needed information about 3D FFT,
!!  npw= number of planewaves in basis for given k point.
!!  paral_kgb
!!  use_gpu_cuda= governs wheter we do the hamiltonian calculation on gpu (1) or not
!!
!! SIDE EFFECTS
!!  ghc(2,npw*ndat)= matrix elements <G|H|C> or <G|H-lambda.S|C> (if sij_opt>=0 or =-1 in getghc)
!!                   contains the fock exchange term for cwavef at the end.
!!
!! NOTES
!!  The current version assumes that :
!!   * nspinor = 1
!!   * no "my_nspinor"
!!   * no restriction to the value of istwfk_bz (but must be tested in all case)
!!   * all the data for the occupied states (cgocc_bz) are the same as those for the current states (cg)
!!
!! PARENTS
!!      getghc
!!
!! CHILDREN
!!      fourdp,fourwf,hartre,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fock_getghc(cwavef,fock,gbound,ghc,gmet,istwf_k,kpoint_i,kg_k,&
&  mgfft,mpi_enreg,n4,n5,n6,nfft,ngfft,npw,paral_kgb,use_gpu_cuda)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_fock
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_getghc'
 use interfaces_18_timing
 use interfaces_53_ffts
 use interfaces_56_xc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! Scalars
 integer, intent(in) :: istwf_k,mgfft,n4,n5,n6,nfft,npw,paral_kgb
 integer, intent(in),optional :: use_gpu_cuda
 type(fock_type),pointer,intent(inout) :: fock
 type(MPI_type),intent(inout) :: mpi_enreg
! Arrays
 integer,intent(in) :: kg_k(3,npw),gbound(2*mgfft+8,2),ngfft(18)
 real(dp),intent(inout) :: cwavef(2,npw),ghc(2,npw)
 real(dp),intent(in) :: gmet(3,3),kpoint_i(3)

!Local variables-------------------------------
! Scalars
 integer,parameter :: tim_fourwf0=0,tim_fourdp0=0,ndat1=1
 integer :: bdtot_jindex,i1,i2,i3,ier,ind,ipw,ifft
 integer :: jband,jbg,jkpt,my_jsppol,nband_k
 real(dp),parameter :: weight1=one
 real(dp) :: eigen,imcwf,imcwocc,imvloc,recwf,recwocc,revloc,occ,wtk
! Arrays
 real(dp) :: qvec_j(3),dummytab2(2,1),tsec(2)
 real(dp), allocatable :: dummytab3(:,:,:),cwavef_r(:,:,:,:),ghc1(:,:)
 real(dp), allocatable :: rhog_munu(:,:),vlocpsi_r(:,:,:,:),work_tmp3(:),vqg(:)
 !character(len=500) :: msg  
 real(dp), ABI_CONTIGUOUS  pointer :: cwaveocc_r(:,:,:,:)

! real(dp), pointer :: cwavef_r(:,:,:,:),dummytab2(:,:),dummytab3(:,:,:)
! real(dp), pointer :: rhog_munu(:,:),vlocpsi_r(:,:,:,:),work_tmp3(:) 

! *************************************************************************
 
 call timab(1504,1,tsec)
 call timab(1505,1,tsec)

 ABI_CHECK(associated(fock),"fock must be associated")

! ===========================
! === Initialize pointers ===
! ===========================

!* Initialization of local pointers
!   cwavef_r   => fock%cwavef_r   !* cwavef_r = current wavefunction in r-space  
!   dummytab2  => fock%dummytab2  !* dummytab2 = variables for fourwf
!   dummytab3  => fock%dummytab3  !* dummytab3 = variables for fourwf
!   rhog_munu  => fock%rhog_munu  !* rhogmunu = overlap matrix between cwavef and (jkpt,mu) in G-space
!   vlocpsi_r  => fock%vlocpsi_r  !* vlocpsi_r = partial local Fock operator (jkpt,mu) in r-space
!   work_tmp3  => fock%work_tmp3  !* work_tmp3 = temporary variables (fftpac)

!* Initializtion of the array cwavef_r
!* cwavef_r = current wavefunction in r-space  
   ABI_ALLOCATE(cwavef_r,(2,n4,n5,n6))

   ABI_ALLOCATE(rhog_munu,(2,nfft))
!* rhogmunu = overlap matrix between cwavef and (jkpt,mu) in G-space
   ABI_ALLOCATE(dummytab3,(n4,n5,n6))
!* dummytab3 = variables for fourwf
!   ABI_ALLOCATE(dummytab2,(2,ndat)) ! max(ndat,ndat_occ)
!* dummytab2 = variables for fourwf
   ABI_ALLOCATE(work_tmp3,(2*nfft))
!* work_tmp3 = temporary variables (fftpac)

   ABI_ALLOCATE(vqg, (nfft))

!* Initialization of the array ghc1
!* ghc1 will contain the exact exchange contribution to the Hamiltonian
   ABI_ALLOCATE(ghc1,(2,npw))

!* Initialization of the array vlocpsi_r
!* vlocpsi_r = partial local Fock operator applied to cwavef in r-space and summed over all occupied (jkpt,mu) 
   ABI_ALLOCATE(vlocpsi_r,(2,n4,n5,n6))
   vlocpsi_r=zero

! ==========================================
! === Get cwavef in real space using FFT ===
! ==========================================
!   fock%cwavef_r=zero
   call fourwf(0,dummytab3,cwavef,dummytab2,cwavef_r,gbound,gbound,istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat1,&
&       ngfft,npw,1,n4,n5,n6,0,paral_kgb,tim_fourwf0,weight1,weight1,use_gpu_cuda=use_gpu_cuda)

! ! =====================================================
! ! === Select the states in cgocc_bz with the same spin ===
! ! =====================================================
! !* Initialization of the indices/shifts, according to the value of isppol
   bdtot_jindex=0
! !* bdtot_jindex = shift to be applied on the location of data in the array occ_bz ?
   jbg=0
! !* jbg = shift to be applied on the location of data in the array cprj/occ

   my_jsppol=fock%isppol
   if((fock%isppol==2).and.(mpi_enreg%nproc_kpt/=1)) my_jsppol=1

   call timab(1505,2,tsec)
   call timab(1506,1,tsec)

! ===================================
! === Loop on the k-points in IBZ ===
! ===================================
   do jkpt=1,fock%mkpt
     nband_k=fock%nbandocc_bz(jkpt,my_jsppol)
!* nband_k = number of bands at point k_j
     wtk=fock%wtk_bz(jkpt)
!* wtk = weight in BZ of this k point

! ======================================
! === Calculate the vector q=k_i-k_j ===
! ======================================
!* Evaluation of kpoint_j, the considered k-point in reduced coordinates 
!     kpoint_j(:)=fock%kptns_bz(:,jkpt)
!* the vector qvec is expressed in reduced coordinates.
!     qvec(:)=kpoint_i(:)-kpoint_j(:)
     qvec_j(:)=kpoint_i(:)-fock%kptns_bz(:,jkpt)

     call bare_vqg(qvec_j,fock%gsqcut,fock%divgq0,gmet,fock%usepaw,nfft,ngfft,vqg) 
     !call fock_vqg_fftbox(qvec_j, nfft, vqg)

! =================================================
! === Loop on the band indices jband of cgocc_k ===
! =================================================
     do jband=1,nband_k

! ==============================================
! === Get cwaveocc_r in real space using FFT ===
! ==============================================
       cwaveocc_r => fock%cwaveocc_bz(:,:,:,:,jband+jbg,my_jsppol)
       occ=fock%occ_bz(jband+bdtot_jindex,my_jsppol)
!* occ = occupancy of jband at this k point

! ================================================
! === Get the overlap density matrix rhog_munu ===
! ================================================
!* Calculate the overlap density matrix in real space = conj(cwaveocc_r)*cwavef_r
!* work_tmp3 will contain the overlap density matrix.
       call timab(1508,1,tsec)
       do i3=1,ngfft(3)
         do i2=1,ngfft(2)
           do i1=1,ngfft(1)
             ind=i1+ngfft(1)*(i2-1+ngfft(2)*(i3-1))

             recwf=cwavef_r(1,i1,i2,i3) ; imcwf=cwavef_r(2,i1,i2,i3)
             recwocc=cwaveocc_r(1,i1,i2,i3) ; imcwocc=cwaveocc_r(2,i1,i2,i3)
             work_tmp3(2*ind-1)= recwocc*recwf+imcwocc*imcwf
             work_tmp3(2*ind)= recwocc*imcwf-imcwocc*recwf
           end do ! i1
         end do ! i2
       end do ! i3
       call timab(1508,2,tsec)

! ===========================================
! === Compute compensation charge density ===
! ===========================================
!       ider=0;izero=0;nhat12_grdim=0
!       nfftf=fock%pawfgr%nfft
!       ngfftf=fock%pawfgr%ngfft
!       natom=fock%natom

!       ABI_ALLOCATE(grnhat12,(2,nfftf,nspinor**2,3*nhat12_grdim))
!       ABI_ALLOCATE(nhat12,(2,nfftf,nspinor**2))

!       call pawmknhat_psipsi(cprjocc(:,???),cwaveprj,ider,izero,natom,natom,nfftf,ngfftf,&
!&        nhat12_grdim,nspinor,fock%ntypat,fock%pawang,fock%pawfgrtab,grnhat12,nhat12,fock%pawtab)

!!Transfer pseudo density from coarse grid to fine grid
!       call transgrid(cplex,mpi_enreg,nspden,+1,1,0,paral_kgb,pawfgr,rhopsg,rhodum,rhopsr,rhor)


! if (present(pawnhat)) then
!   pawnhat_ptr => pawnhat
! else
!   ABI_ALLOCATE(pawnhat_ptr,(pawfgr%nfft,nspden))
! end if
! if (present(pawrhoij0)) then
!   pawrhoij0_ptr => pawrhoij0
! else
!   pawrhoij0_ptr => pawrhoij_ptr
! end if
! call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,natom,&
!& pawfgr%nfft,pawfgr%ngfft,ider,nspden,ntypat,pawang,pawfgrtab,&
!& rhodum,pawnhat_ptr,pawrhoij_ptr,pawrhoij0_ptr,pawtab,qphon,rprimd,ucvol,usewvl,xred,&
!& mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
!& mpi_comm_fft=mpi_enreg%comm_fft,paral_kgb=paral_kgb,me_g0=mpi_enreg%me_g0,&
!& distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)

! call pawmknhat_psipsi(cprjocc,cwaveprj,ider,izero,my_natom,natom,pawfgr%nfft,pawfgr%ngfft,nhat12_grdim,&
!&          nspinor,ntypat,pawang,pawfgrtab,grnhat12,nhat12,pawtab)

!!Transfer pseudo density from coarse grid to fine grid
! if(usewvl==0) then
!   call transgrid(cplex,mpi_enreg,nspden,+1,1,0,paral_kgb,pawfgr,rhopsg,rhodum,rhopsr,rhor)
! end if

!!Add pseudo density and compensation charge density (on fine grid)
! rhor(:,:)=rhor(:,:)+pawnhat_ptr(:,:)

!!Free temporary memory spaces
! if (.not.present(pawnhat))  then
!   ABI_DEALLOCATE(pawnhat_ptr)
! end if

!* Perform an FFT using fourwf to get rhog_munu = FFT^-1(work_tmp3)
       call timab(1509,1,tsec)
       call fourdp(2,rhog_munu,work_tmp3,-1,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp0)
       call timab(1509,2,tsec)

! ===================================================
! === Calculate the local potential vfockloc_munu ===
! ===================================================
!* Apply the Poisson solver to "rhog_munu" while taking into account the effect of the vector "qvec"
!* This is precisely what is done in the subroutine hartre, with option cplex=2.
!* work_tmp3 will contain the local Fock potential, the result of hartre routine.
!* work_tmp3 = FFT( rhog_munu/|g+qvec|^2 )


       call timab(1510,1,tsec)
#if 0 
!         rhog_munu=rhog_munu*fock%wtk_bz(jkpt)
       call hartre(2,gmet,fock%gsqcut,fock%usepaw,mpi_enreg,nfft,ngfft,&
&        paral_kgb,qvec_j,rhog_munu,work_tmp3,divgq0=fock%divgq0)
#else
       do ifft=1,nfft
         rhog_munu(1,ifft) = rhog_munu(1,ifft) * vqg(ifft)
         rhog_munu(2,ifft) = rhog_munu(2,ifft) * vqg(ifft)
       end do
       call fourdp(2,rhog_munu,work_tmp3,+1,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp0)
#endif
       call timab(1510,2,tsec)

!* CMartins : the variable izero is set to gs_ham%usepaw as in all other routines.

!Transfer pseudo density from coarse grid to fine grid
!  if PAW  call transgrid(cplex,mpi_enreg,nspden,-1,0,0,paral_kgb,pawfgr,rhopsg,rhodum,rhopsr,rhor)

! =============================================================
! === Apply the local potential vfockloc_munu to cwaveocc_r ===
! =============================================================
       call timab(1507,1,tsec)
       do i3=1,ngfft(3) 
         do i2=1,ngfft(2) 
           do i1=1,ngfft(1) 
             ind=i1+ngfft(1)*(i2-1+ngfft(2)*(i3-1))
             revloc=work_tmp3(2*ind-1) ; imvloc=work_tmp3(2*ind)

             recwocc=cwaveocc_r(1,i1,i2,i3) ; imcwocc=cwaveocc_r(2,i1,i2,i3)
             vlocpsi_r(1,i1,i2,i3)=vlocpsi_r(1,i1,i2,i3)-(revloc*recwocc-imvloc*imcwocc)*occ*wtk
             vlocpsi_r(2,i1,i2,i3)=vlocpsi_r(2,i1,i2,i3)-(revloc*imcwocc+imvloc*recwocc)*occ*wtk
           end do
         end do
       end do
       call timab(1507,2,tsec)
   
! ==========================================================================
! === Evaluate ghc_munu, the contribution to ghc of the occupied band mu ===
! ==========================================================================
!* initialize ghc_munu to zero
!       ghc_munu=zero
!* ghc_munu will contain the contribution to the Hamiltonian of each band nu
!       weight=one
!* Perform the FFT using fourwf, ghc_munu = FFT^-1(vlocpsi_r)
!       call fourwf(0,dummytab3,dummytab2,ghc_munu,vlocpsi_r,&
!&        gbound,gbound,istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat1,ngfft,1,npw,n4,n5,n6,3,&
!&        paral_kgb,tim_fourwf0,weight,weight,use_gpu_cuda=use_gpu_cuda)

! =========================================
! === Sum the contribution of each band ===
! =========================================
!* The value of ghc_munu is added to ghc ; 
!* it is weighted by the occupancy of the band occ_bz(jband) and the weight of kpoint_j wtk_bz(jkpt)
!       ghc=ghc-ghc_munu*fock%occ_bz(jband+bdtot_jindex,my_jsppol)*fock%wtk_bz(jkpt)

     end do ! jband

! ================================================
! === End : update of shifts and deallocations ===
! ================================================
!* Update of the shifts to be applied (reminder : mkmem is not 0, nspinor=1)
     jbg=jbg+nband_k
     bdtot_jindex=bdtot_jindex+nband_k
   end do ! jkpt

   call timab(1506,2,tsec)
   call timab(1511,1,tsec)

!* Perform an FFT using fourwf to get ghc1 = FFT^-1(vlocpsi_r)
   call fourwf(0,dummytab3,dummytab2,ghc1,vlocpsi_r,gbound,gbound,istwf_k,kg_k,kg_k,mgfft,mpi_enreg,ndat1,&
&       ngfft,1,npw,n4,n5,n6,3,paral_kgb,tim_fourwf0,weight1,weight1,use_gpu_cuda=use_gpu_cuda)

!* If the calculation is parallelized, perform an MPI_allreduce to sum all the contributions in the array ghc
   ghc(:,:)=ghc(:,:)/mpi_enreg%nproc_hf + fock%alpha*ghc1(:,:)
   call xmpi_sum(ghc,mpi_enreg%comm_hf,ier)
   
   call timab(1511,2,tsec)

! ============================================
! === Calculate the contribution to energy ===
! ============================================
!* Only the contribution when cwavef=cgocc_bz are calculated, in order to cancel exactly the self-interaction 
!* at each convergence step. (consistent definition with the defintion of hartree energy)
   if (fock%ieigen/=0) then
     eigen=0.d0
!* Dot product of cwavef and ghc 
!* inspired from the routine 53_spacepar/meanvalue_g but without the reference to parallelism and filtering
     if(istwf_k==2) then 
       eigen=half*cwavef(1,1)*ghc1(1,1)
     else
       eigen=cwavef(1,1)*ghc1(1,1)+cwavef(2,1)*ghc1(2,1)
     end if
     do ipw=2,npw
       eigen=eigen+cwavef(1,ipw)*ghc1(1,ipw)+cwavef(2,ipw)*ghc1(2,ipw)
     end do
     if(istwf_k>=2) eigen=two*eigen
     call xmpi_sum(eigen,mpi_enreg%comm_hf,ier)
     fock%eigen_ikpt(fock%ieigen)= eigen
     fock%ieigen = 0 
   end if

! ===============================
! === Deallocate local arrays ===
! ===============================
   ABI_DEALLOCATE(cwavef_r)
   ABI_DEALLOCATE(ghc1)
   ABI_DEALLOCATE(rhog_munu)
   ABI_DEALLOCATE(vlocpsi_r)
   ABI_DEALLOCATE(dummytab3)
   ABI_DEALLOCATE(work_tmp3)
   ABI_DEALLOCATE(vqg)

 call timab(1504,2,tsec)

end subroutine fock_getghc
!!***
