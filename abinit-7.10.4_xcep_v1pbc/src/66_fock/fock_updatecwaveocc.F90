!{\src2tex{textfont=tt}}
!!****f* ABINIT/fock_updatecwaveocc
!! NAME
!!  fock_updatecwaveocc
!!
!! FUNCTION
!!  Update in the fock datastructure the fields relative to the occupied states.
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (CMartins)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mcg)= Input wavefunctions
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!  istep=index of the number of steps in the routine scfcv
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  occ(mband*nkpt*nsppol)= occupation number for each band (often 2) at each k point
!!
!!  fock_energy= ???
!!  nsppol=number of independent spin polarizations
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_exactX = Fock contribution to the total energy (Hartree)
!!
!! OUTPUT
!!  none
!!
!! SIDE EFFECTS
!!   The field fock%cgocc_bz contains the table cg at the end.
!!   The fields kg_bz, occ_bz are simultaneously updated. 
!!
!! NOTES
!!
!!  ############################
!!  ### Not fully tested yet ###
!!  ############################
!!
!! May be improved by selecting only the occupied states with the same spin isppol.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      fourwf,timab,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fock_updatecwaveocc(cg,dtset,fock,fock_energy,istep,mcg,mpi_enreg,npwarr,occ)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_fock
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_updatecwaveocc'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!scalars
 integer, intent(in) :: istep,mcg
 real(dp) :: fock_energy
 type(dataset_type),intent(in) :: dtset
 type(fock_type),intent(inout),pointer :: fock
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer, intent(in) :: npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg),occ(dtset%mband*dtset%nkpt*dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf0=0
 integer :: iband,ibg,icg,ier,ikpt,isppol,jbg,jcg,jkg,jkpt,jpw,jstwfk
 integer :: mband,mkpt,mpw,my_jsppol,my_jband,my_jkpt
 integer :: n4,n5,n6,nkpt_bz,npwj,nsppol
 real(dp),parameter :: weight1=one
 real(dp) :: cgre,cgim 
 character(len=500) :: message               
! arrays
 integer, ABI_CONTIGUOUS pointer :: gbound_k(:,:),kg_k(:,:)
 real(dp) :: tsec(2),tsec2(2)
 real(dp),allocatable :: cgocc_tmp(:),cgocc(:,:),dummytab2(:,:),dummytab3(:,:,:),phase_jkpt(:,:)

! *************************************************************************
 
! DEBUG
 write (std_out,*) ' fock_updatecwaveocc : enter'
! ENDDEBUG

 call timab(1502,1,tsec)

 ABI_CHECK(associated(fock),"fock must be associated")

 if (associated(fock)) then

   if (mod(istep-1,fock%nnsclo_hf)==0) then 

! Local variables = useful dimensions
     mband=fock%mband
     mkpt=fock%mkpt
     mpw=dtset%mpw
     nkpt_bz=fock%nkpt_bz
     nsppol=dtset%nsppol

! Local variables : useful arrays
     ABI_ALLOCATE(cgocc,(2,mpw))
     cgocc=zero
     ABI_ALLOCATE(cgocc_tmp,(2*mpw+1))
     cgocc_tmp=zero

! Local variables to perform FFT
     n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
     ABI_ALLOCATE(dummytab3,(n4,n5,n6))
     ABI_ALLOCATE(dummytab2,(2,1))

     if(ANY(fock%calc_phase(:)/=0)) then
       ABI_ALLOCATE(phase_jkpt,(2,mpw))
       phase_jkpt=zero
     end if

! =======================================================
! === Update the data relative to the occupied states ===
! =======================================================
!* The arrays cgocc_bz, kg_bz, occ_bz and npwarr_bz are already allocated with the maximal size. 
!     if ((dtset%kptopt>=1).and.(dtset%kptopt<=4)) then
!       if (dtset%kptopt/=3) then
     do isppol=1,nsppol
       jbg=0 ; jcg=0 ; jkg=0
          
       my_jsppol=isppol
       if ((isppol==2).and.(mpi_enreg%nproc_kpt/=1)) my_jsppol=1
!* Both spins are treated on the same proc., only in the case where nproc_kpt=1; 
!* otherwise each proc. treats only one spin.

       ! MG: This loop is not effient!
       ! Ok the number of k-points in the BZ is usually `small` when hybrids are used
       ! but what happens if we use a 12x12x12.
       ! One should loop over the IBZ, broadcast and reconstruct the star of the k-point.
       my_jkpt=0
       do jkpt=1,nkpt_bz
           
       if (proc_distrb_cycle(mpi_enreg%distrb_hf,jkpt,1,mband,1,mpi_enreg%me_hf)) cycle
!* In this case, the processor does not calculate the exchange with any occupied state on jkpt.

!               if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,jkpt,1,dtset%nbandhf,jsppol,mpi_enreg%me_kpt))) then
!* The state (jkpt,jband,jsppol) is stored in the array cg of this processor and copied in cgocc_tmp.
!                 icg=icg+dtset%nband(jkpt)*npwj
!               end if
!               ibg=ibg+dtset%nband(jkpt)
!               cycle
!             end if
       my_jkpt=my_jkpt+1
      
       ikpt=fock%tab_ikpt(my_jkpt)
!* ikpt = the point of IBZ that jkpt is an image of in BZ
       npwj=npwarr(ikpt)
!* npwj= number of plane wave in basis for the wavefunction 
       jstwfk=fock%istwfk_bz(my_jkpt)     
!* jstwfk= how is stored the wavefunction 
       ibg=fock%tab_ibg(my_jkpt,my_jsppol)
!* ibg = shift to be applied on the location of data in the array cprj/occ
       icg=fock%tab_icg(my_jkpt,my_jsppol)
!* icg = shift to be applied on the location of data in the array cg
       gbound_k => fock%gbound_bz(:,:,my_jkpt)
!* boundary of the basis sphere of G vectors
       kg_k => fock%kg_bz(:,1+jkg:npwj+jkg)
!* reduced plean wave coordinates
       if (fock%calc_phase(my_jkpt)==1) then
         phase_jkpt(:,1:npwj)=fock%phase(:,1+jkg:npwj+jkg)
       end if
!* phase factor at k-point j

!* Initialize the band counter
       my_jband=0
       do iband=1,dtset%nband(ikpt)
         cgocc_tmp=zero
             
         if(ABS(occ(iband+ibg))>tol8) then
!* If the band is occupied

!* To avoid segmentation fault, my_jband should not be greater than nbandhf
           if ((my_jband+1)>mband) then
             write(message,*) 'The number of occupied band at k-point',ikpt,' is greater than the value of nbandhf'
             MSG_ERROR(message)
           end if

!* If the processor does not calculate the exchange with the occupied state (jkpt,my_jband), cycle
           if (mpi_enreg%distrb_hf(jkpt,(my_jband+1),1)/=mpi_enreg%me_hf) cycle
!                   if (mpi_enreg%proc_distrb(jkpt,jband,jsppol)==mpi_enreg%me_kpt) then
!* The state (jkpt,jband,jsppol) is stored in the array cg of this processor ; shift are incremented.
!                     icg=icg+npwj
!                   end if
!                   ibg=ibg+1
!* Skip the end of the loop
!                   cycle
!                 end if

!* increment the number of occupied bands treated on this processor
           my_jband = my_jband+1

!* In this case, the processor calculates the exchange with the occupied state (jkpt,my_jband). 
           if (mpi_enreg%proc_distrb(ikpt,iband,isppol)==mpi_enreg%me_kpt) then
!* The state (ikpt,iband,isppol) is stored in the array cg of this processor and copied in cgocc_tmp.
             if(icg==-1) then
               write(100,*) 'icg=-1',mpi_enreg%me,isppol,my_jsppol,jkpt,my_jkpt,ikpt,iband
             end if
             ! MG: Why packing re and im part?
             cgocc_tmp(1)=occ(iband+ibg)
             cgocc_tmp(2:npwj+1)=cg(1,1+(iband-1)*npwj+icg:iband*npwj+icg)
             cgocc_tmp(npwj+2:2*npwj+1)=cg(2,1+(iband-1)*npwj+icg:iband*npwj+icg)
           end if
!* Broadcast the state (ikpt,iband,isppol) to all the processors of comm_kpt.
           call timab(1503,1,tsec2)
           call xmpi_bcast(cgocc_tmp,mpi_enreg%proc_distrb(ikpt,iband,isppol),mpi_enreg%comm_kpt,ier)
           call timab(1503,2,tsec2)

!* Keep the processors in %comm_kpt which needs the values in cgocc_tmp to build their own %cwaveocc and %occ_bz. 
           if ((mpi_enreg%nproc_kpt/=1).and.(nsppol==2)) then
             if (fock%timerev(my_jkpt)==mpi_enreg%my_isppoltab(isppol)) cycle
!* In the case of a parallel spin-polarized calculation 
!* when time reversal symmetry is applied at this k-point (timrev==1), only the processors with the opposite spin (my_isppoltab==0) are kept.
!* when time reversal symmetry is not applied at this k-point (timrev==0), only the processors with the same spin (my_isppoltab==1) are kept.

!             if (fock%timerev(my_jkpt)==1)) then
!               if (mpi_enreg%my_isppoltab(isppol)==1) cycle
!* In the case of a parallel spin-polarized calculation and when time reversal symmetry is applied at this k-point,
!* only the processors with the opposite spin are kept.
!             else
!               if (mpi_enreg%my_isppoltab(isppol)==0) cycle
!* only the processors with isppol are kept.
!             end if
           end if

!* Copy the values of cgocc_tmp in the arrays cgocc and %occ_bz
           fock%occ_bz(my_jband+jbg,my_jsppol) = cgocc_tmp(1)
           cgocc(1,1:npwj) = cgocc_tmp(2:npwj+1)
           cgocc(2,1:npwj) = cgocc_tmp(npwj+2:2*npwj+1)

!* calculate cg and store it in cgocc_bz
           if (fock%calc_phase(my_jkpt)==1) then
             do jpw=1,npwj
               cgre=cgocc(1,jpw) ; cgim=cgocc(2,jpw)
               cgocc(1,jpw) = phase_jkpt(1,jpw)*cgre - phase_jkpt(2,jpw)*cgim
               cgocc(2,jpw) = phase_jkpt(1,jpw)*cgim + phase_jkpt(2,jpw)*cgre
             end do
           end if ! phase

!* apply time reversal symmetry if necessary
           if (fock%timerev(my_jkpt)==1) then 
             cgocc(2,:) = - cgocc(2,:)
             if((mpi_enreg%nproc_kpt==1).and.(nsppol==2)) my_jsppol=mod(my_jsppol,2)+1
!* exchange spin (1 ->2 ; 2-> 1) in the sequential case.
           end if
                
!* apply FFT to get cwaveocc in real space
           call fourwf(0,dummytab3,cgocc,dummytab2,fock%cwaveocc_bz(:,:,:,:,my_jband+jbg,my_jsppol), &
&            gbound_k,gbound_k,jstwfk,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&            npwj,1,n4,n5,n6,tim_fourwf0,dtset%paral_kgb,0,weight1,weight1,use_gpu_cuda=dtset%use_gpu_cuda)

!               else
!* The band is empty ; the array cgocc_tmp remains equal to 0.d0.
!                 if (mpi_enreg%proc_distrb(jkpt,jband,jsppol)==mpi_enreg%me_kpt) then
!* The state (jkpt,jband,jsppol) is stored in the array cg of this processor ; shift are incremented.
!                   icg=icg+npwj
!                 end if

         end if ! band occupied

!* update the shift to apply to occ in all case because this array is not distributed among the proc.
!               ibg=ibg+1

       end do ! iband

!* Save the true number of occupied bands in the array %nbandocc_bz
       fock%nbandocc_bz(my_jkpt,my_jsppol) = my_jband

!* update the shifts to apply
       jbg=jbg+my_jband
       jcg=jcg+npwj*my_jband
       jkg=jkg+npwj
       end do ! ikpt
     end do ! isppol
     ABI_DEALLOCATE(cgocc_tmp)
     ABI_DEALLOCATE(cgocc)
     if(allocated(phase_jkpt)) then 
       ABI_DEALLOCATE(phase_jkpt)
     end if
     ABI_DEALLOCATE(dummytab3)
     ABI_DEALLOCATE(dummytab2)

! Restricted or unrestricted HF
     if (nsppol==1) then 
!* Update the array %occ_bz => May be limited to the occupied states only 
      fock%occ_bz(:,:)=half*fock%occ_bz(:,:)
! If nsppol=1, this is a restricted Hartree-Fock calculation.
! If nsppol=2, this is an unrestricted Hartree-Fock calculation.
     end if

   end if ! istep

!* Set the Fock contribution to total energy to zero
   fock_energy=zero

 end if

 call timab(1502,2,tsec)
! DEBUG
write (std_out,*) ' fock_updatecwaveocc : exit'
! stop
! ENDDEBUG

end subroutine fock_updatecwaveocc
!!***
