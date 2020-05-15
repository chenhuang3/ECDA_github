!> @file
!! H|psi> and orthonormalization (linear scaling version)
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
           ldiis, fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, alpha_mean, alpha_max, &
           energy_increased, tmb, lhphiold, overlap_calculated, &
           energs, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, &
           energy_only, hpsi_small, experimental_mode, ksorbs, hpsi_noprecond)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, it
  type(DFT_wavefunction), target, intent(inout):: tmb
  type(localizedDIISParameters), intent(inout) :: ldiis
  real(kind=8), dimension(tmb%orbs%norb), intent(inout) :: fnrmOldArr
  real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha
  real(kind=8), intent(out):: trH, fnrm, fnrmMax, alpha_mean, alpha_max
  real(kind=8), intent(inout):: trHold
  logical,intent(out) :: energy_increased
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout):: lhphiold
  logical,intent(inout):: overlap_calculated
  type(energy_terms), intent(in) :: energs
  real(kind=8), dimension(:), pointer:: hpsit_c, hpsit_f
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint
  logical, intent(in) :: energy_only, experimental_mode
  real(kind=8), dimension(tmb%npsidim_orbs), intent(out) :: hpsi_small
  type(orbitals_data),intent(in) :: ksorbs
  real(kind=8), dimension(tmb%npsidim_orbs), optional,intent(out) :: hpsi_noprecond

  ! Local variables
  integer :: iorb, iiorb, ilr, ncount, ierr, ist, ncnt, istat, iall, ii, jjorb, i, jorb
  integer :: matrixindex_in_compressed, lwork, info
  real(kind=8) :: ddot, tt, gnrmArr, fnrmOvrlp_tot, fnrm_tot, fnrmold_tot
  !real(kind=8) :: eval_zero
  character(len=*), parameter :: subname='calculate_energy_and_gradient_linear'
  real(kind=8), dimension(:), pointer :: hpsittmp_c, hpsittmp_f
  real(kind=8), dimension(:), allocatable :: fnrmOvrlpArr, fnrmArr, work
  real(kind=8), dimension(:), allocatable :: hpsi_conf, hpsi_tmp
  real(kind=8), dimension(:), pointer :: kernel_compr_tmp
  !type(sparseMatrix) :: lagmat
  real(kind=8), dimension(:), allocatable :: prefac
  real(wp), dimension(2) :: garray
  real(dp) :: gnrm,gnrm_zero,gnrmMax,gnrm_old ! for preconditional2, replace with fnrm eventually, but keep separate for now
  real(8),dimension(:),allocatable :: prefacarr, dphi, dpsit_c, dpsit_f
  real(kind=8),dimension(:,:),allocatable :: SK, KS, HK, KHK, KSKHK, KHKSK , Q
  integer,dimension(:),allocatable :: ipiv
  real(kind=8) :: fnrm_low, fnrm_high, fnrm_in, fnrm_out, rx, ry, rz, rr, hh, fnrm_tot2
  integer :: iseg, isegf, j0, jj, j1, i1, i2, i3, i0, istart, iold, inew, ind_ham, ind_denskern, iorbp
  real(kind=8),dimension(3) :: noise

  if (target_function==TARGET_FUNCTION_IS_HYBRID) then
      allocate(hpsi_conf(tmb%npsidim_orbs), stat=istat)
      call memocc(istat, hpsi_conf, 'hpsi_conf', subname)
      call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%hpsi, hpsi_conf)
      call timing(iproc,'eglincomms','ON')
      ist=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          tt=ddot(ncount, hpsi_conf(ist), 1, tmb%psi(ist), 1)
          call daxpy(ncount, -tt, tmb%psi(ist), 1, hpsi_conf(ist), 1)
          ist=ist+ncount
      end do
      call timing(iproc,'eglincomms','OF')
  end if

  ! by default no quick exit
  energy_increased=.false.


  allocate(hpsittmp_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, hpsittmp_c, 'hpsittmp_c', subname)
  allocate(hpsittmp_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, hpsittmp_f, 'hpsittmp_f', subname)

  if(target_function==TARGET_FUNCTION_IS_ENERGY .or. &
     target_function==TARGET_FUNCTION_IS_HYBRID) then

      if(sum(tmb%ham_descr%collcom%nrecvcounts_c)>0) &
          call dcopy(sum(tmb%ham_descr%collcom%nrecvcounts_c), hpsit_c(1), 1, hpsittmp_c(1), 1)
      if(sum(tmb%ham_descr%collcom%nrecvcounts_f)>0) &
          call dcopy(7*sum(tmb%ham_descr%collcom%nrecvcounts_f), hpsit_f(1), 1, hpsittmp_f(1), 1)

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call timing(iproc,'eglincomms','ON')
          allocate(kernel_compr_tmp(tmb%linmat%denskern%nvctr), stat=istat)
          call memocc(istat, kernel_compr_tmp, 'kernel_compr_tmp', subname)
          call vcopy(tmb%linmat%denskern%nvctr, tmb%linmat%denskern%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
          !ii=0
          !do iseg=1,tmb%linmat%denskern%nseg
          !    do jorb=tmb%linmat%denskern%keyg(1,iseg), tmb%linmat%denskern%keyg(2,iseg)
          !        ii=ii+1
          !        iiorb = (jorb-1)/tmb%orbs%norb + 1
          !        jjorb = jorb - (iiorb-1)*tmb%orbs%norb
              do ii=1,tmb%linmat%denskern%nvctr
                      iiorb = tmb%linmat%denskern%orb_from_index(1,ii)
                      jjorb = tmb%linmat%denskern%orb_from_index(2,ii)
                  if(iiorb==jjorb) then
                  !if(iiorb==jjorb .or. mod(iiorb-1,9)+1>4 .or.  mod(jjorb-1,9)+1>4) then
                  !if(iiorb==jjorb .or. mod(iiorb-1,9)+1>4) then
                      tmb%linmat%denskern%matrix_compr(ii)=0.d0
                  else
                      tmb%linmat%denskern%matrix_compr(ii)=kernel_compr_tmp(ii)
                  end if
              end do
          !end do

          ist=1
          do iorb=tmb%orbs%isorb+1,tmb%orbs%isorb+tmb%orbs%norbp
              ilr=tmb%orbs%inwhichlocreg(iorb)
              !ii=0
              !do iseg=1,tmb%linmat%denskern%nseg
              !    do jorb=tmb%linmat%denskern%keyg(1,iseg), tmb%linmat%denskern%keyg(2,iseg)
              !        ii=ii+1
              !        iiorb = (jorb-1)/tmb%orbs%norb + 1
              !        jjorb = jorb - (iiorb-1)*tmb%orbs%norb
              do ii=1,tmb%linmat%denskern%nvctr
                      iiorb = tmb%linmat%denskern%orb_from_index(1,ii)
                      jjorb = tmb%linmat%denskern%orb_from_index(2,ii)
                      if(iiorb==jjorb .and. iiorb==iorb) then
                          ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
                          call dscal(ncount, kernel_compr_tmp(ii), tmb%hpsi(ist), 1)
                          ist=ist+ncount
                      end if
                  !end do
              end do
          end do
          call timing(iproc,'eglincomms','OF')
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
          call build_linear_combination_transposed(tmb%ham_descr%collcom, &
               tmb%linmat%denskern, hpsittmp_c, hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
          ! copy correct kernel back
          call vcopy(tmb%linmat%denskern%nvctr, kernel_compr_tmp(1), 1, tmb%linmat%denskern%matrix_compr(1), 1)
          iall=-product(shape(kernel_compr_tmp))*kind(kernel_compr_tmp)
          deallocate(kernel_compr_tmp, stat=istat)
          call memocc(istat, iall, 'kernel_compr_tmp', subname)
      else
          call build_linear_combination_transposed(tmb%ham_descr%collcom, &
               tmb%linmat%denskern, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
      end if
  end if

  !!! EXPERIMENTAL: correction for co- / contravariant ===============================================================
  !!! Calculate the overlap matrix, can be optimized ############################
  !!! Use ham since it has the correct SHAMOP pattern
  !!!if(.not.tmb%ham_descr%can_use_transposed) then
  !!    if(.not.associated(tmb%ham_descr%psit_c)) then
  !!        allocate(tmb%ham_descr%psit_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
  !!        call memocc(istat, tmb%ham_descr%psit_c, 'tmb%ham_descr%psit_c', subname)
  !!    end if
  !!    if(.not.associated(tmb%ham_descr%psit_f)) then
  !!        allocate(tmb%ham_descr%psit_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
  !!        call memocc(istat, tmb%ham_descr%psit_f, 'tmb%ham_descr%psit_f', subname)
  !!    end if
  !!!end if
  !!call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
  !!     tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd)
  !!tmb%ham_descr%can_use_transposed=.true.

  !!call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, tmb%ham_descr%psit_c, &
  !!     tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%psit_f, tmb%linmat%ham)

  !!!!%%! invert overlap
  !!!!%%allocate(tmb%linmat%ham%matrix(tmb%orbs%norb,tmb%orbs%norb))
  !!!!%%call uncompressMatrix(iproc, tmb%linmat%ham)
  !!!!%%allocate(ipiv(tmb%orbs%norb))
  !!!!%%lwork=10*tmb%orbs%norb
  !!!!%%allocate(work(lwork))
  !!!!%%call dgetrf(tmb%orbs%norb, tmb%orbs%norb, tmb%linmat%ham%matrix, tmb%orbs%norb, ipiv, info)
  !!!!%%call dgetri(tmb%orbs%norb, tmb%linmat%ham%matrix, tmb%orbs%norb, ipiv, work, lwork, info)
  !!!!%%call compress_matrix_for_allreduce(iproc,tmb%linmat%ham)
  !!!!%%deallocate(ipiv)
  !!!!%%deallocate(work)
  !!!!%%deallocate(tmb%linmat%ham%matrix)


  !!do ii=1,tmb%linmat%ham%nvctr
  !!   !iorb = tmb%linmat%ham%orb_from_index(1,ii)
  !!   !jorb = tmb%linmat%ham%orb_from_index(2,ii)
  !!   !if (iproc==0) write(333,'(a,2i8,es16.6)') 'iorb, jorb, matrix', iorb, jorb, tmb%linmat%ham%matrix_compr(ii)
  !!   !!if (iorb==jorb) then
  !!   !!    tmb%linmat%ham%matrix_compr(ii)=1.d0
  !!   !!else
  !!   !!    tmb%linmat%ham%matrix_compr(ii)=0.d0
  !!   !!end if
  !!end do
  !!!if (iproc==0) write(333,'(a)') '========================================='

  !!!! ###########################################################################

  !!hpsittmp_c=hpsit_c
  !!hpsittmp_f=hpsit_f
  !!call build_linear_combination_transposed(tmb%ham_descr%collcom, tmb%linmat%ham, &
  !!     hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
  !!! END EXPERIMENTAL ===============================================================================================

  iall=-product(shape(hpsittmp_c))*kind(hpsittmp_c)
  deallocate(hpsittmp_c, stat=istat)
  call memocc(istat, iall, 'hpsittmp_c', subname)
  iall=-product(shape(hpsittmp_f))*kind(hpsittmp_f)
  deallocate(hpsittmp_f, stat=istat)
  call memocc(istat, iall, 'hpsittmp_f', subname)




  !!! EXPERIMENTAL: add the term stemming from the derivative of the kernel with respect to the support funtions #################

  !!if (iproc==0) write(*,*) 'kernel term...'

  !!! Calculate the matrix Q = 3KHK - 2KSKHK - 2KHKSK
  !!if(.not. tmb%ham_descr%can_use_transposed) then
  !!    allocate(tmb%ham_descr%psit_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
  !!    !call memocc(istat, tmb%ham_descr%psit_c, 'tmb%ham_descr%psit_c', subname)

  !!    allocate(tmb%ham_descr%psit_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
  !!    !call memocc(istat, tmb%ham_descr%psit_f, 'tmb%ham_descr%psit_f', subname)

  !!    call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, tmb%ham_descr%psi, &
  !!         tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd)
  !!    !can_use_transposed=.true.
  !!end if

  !!! It is assumed that this routine is called with the transposed gradient ready if it is associated...
  !!if(.not.associated(hpsit_c)) then
  !!    allocate(hpsit_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
  !!    !call memocc(istat, hpsit_c, 'hpsit_c', subname)
 
  !!    allocate(hpsit_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
  !!    !call memocc(istat, hpsit_f, 'hpsit_f', subname)
 
  !!   call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, tmb%hpsi, &
  !!        hpsit_c, hpsit_f, tmb%ham_descr%lzd)
  !!end if

  !!if(.not.tmb%can_use_transposed) then
  !!    if(.not.associated(tmb%psit_c)) then
  !!        allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
  !!        call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
  !!    end if
  !!    if(.not.associated(tmb%psit_f)) then
  !!        allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
  !!        call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
  !!    end if
  !!    call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
  !!         tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
  !!    tmb%can_use_transposed=.true.
  !!end if

  !!call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
  !!     tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%ovrlp)
  !!call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, tmb%ham_descr%psit_c, hpsit_c, &
  !!     tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%ham)
  !!!!call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, tmb%ham_descr%psit_c, &
  !!!!     tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%psit_f, tmb%linmat%ovrlp)

  !!allocate(tmb%linmat%ham%matrix(tmb%orbs%norb,tmb%orbs%norb))
  !!call uncompressMatrix(iproc,tmb%linmat%ham)
  !!allocate(tmb%linmat%denskern%matrix(tmb%orbs%norb,tmb%orbs%norb))
  !!call uncompressMatrix(iproc,tmb%linmat%denskern)
  !!allocate(tmb%linmat%ovrlp%matrix(tmb%orbs%norb,tmb%orbs%norb))
  !!call uncompressMatrix(iproc,tmb%linmat%ovrlp)

  !!allocate(SK(tmb%orbs%norb,tmb%orbs%norb))
  !!allocate(KS(tmb%orbs%norb,tmb%orbs%norb))
  !!allocate(HK(tmb%orbs%norb,tmb%orbs%norb))
  !!allocate(KHK(tmb%orbs%norb,tmb%orbs%norb))
  !!allocate(KSKHK(tmb%orbs%norb,tmb%orbs%norb))
  !!allocate(KHKSK(tmb%orbs%norb,tmb%orbs%norb))
  !!allocate(Q(tmb%orbs%norb,tmb%orbs%norb))

  !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, tmb%linmat%ovrlp%matrix, tmb%orbs%norb, &
  !!     tmb%linmat%denskern%matrix, tmb%orbs%norb, 0.d0, SK, tmb%orbs%norb)
  !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, tmb%linmat%denskern%matrix, tmb%orbs%norb, &
  !!     tmb%linmat%ovrlp%matrix, tmb%orbs%norb, 0.d0, KS, tmb%orbs%norb)
  !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, tmb%linmat%ham%matrix, tmb%orbs%norb, &
  !!     tmb%linmat%denskern%matrix, tmb%orbs%norb, 0.d0, HK, tmb%orbs%norb)
  !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, tmb%linmat%denskern%matrix, tmb%orbs%norb, &
  !!     HK, tmb%orbs%norb, 0.d0, KHK, tmb%orbs%norb)
  !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, KS, tmb%orbs%norb, &
  !!     KHK, tmb%orbs%norb, 0.d0, KSKHK, tmb%orbs%norb)
  !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, KHK, tmb%orbs%norb, &
  !!     SK, tmb%orbs%norb, 0.d0, KHKSK, tmb%orbs%norb)

  !!Q = 3*KHK - 2*KSKHK -2*KHKSK

  !!! Store the matrix Q temporaily in tmb%linmat%ham
  !!tmb%linmat%ham%matrix=Q
  !!call compress_matrix_for_allreduce(iproc,tmb%linmat%ham)

  !!call build_linear_combination_transposed(tmb%ham_descr%collcom, tmb%linmat%ham, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, &
  !!     .false., hpsit_c, hpsit_f, iproc)

  !!call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
  !!     hpsit_c, hpsit_f, tmb%hpsi, tmb%ham_descr%lzd)

  !!deallocate(tmb%linmat%ham%matrix)
  !!deallocate(tmb%linmat%denskern%matrix)
  !!deallocate(tmb%linmat%ovrlp%matrix)

  !!deallocate(SK)
  !!deallocate(KS)
  !!deallocate(HK)
  !!deallocate(KHK)
  !!deallocate(KSKHK)
  !!deallocate(KHKSK)
  !!deallocate(Q)

  !!deallocate(hpsit_c)
  !!deallocate(hpsit_f)

  !!! not sure about this
  !!if(.not. tmb%ham_descr%can_use_transposed) then
  !!    deallocate(tmb%ham_descr%psit_c)
  !!    deallocate(tmb%ham_descr%psit_f)
  !!end if

  !!! END EXPERIMENTAL ###########################################################################################################


  ! make lagmat a structure with same sparsity as h
  !call nullify_sparsematrix(lagmat)
  !call sparse_copy_pattern(tmb%linmat%ham, lagmat, iproc, subname)
  !allocate(lagmat%matrix_compr(lagmat%nvctr), stat=istat)
  !call memocc(istat, lagmat%matrix_compr, 'lagmat%matrix_compr', subname)

  ! Calculate the overlap matrix, can be optimized ############################
  if (correction_orthoconstraint==0) then
      !if(.not.tmb%can_use_transposed) then
          if(.not.associated(tmb%psit_c)) then
              allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
              call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
          end if
          if(.not.associated(tmb%psit_f)) then
              allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
              call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
          end if
      !end if
      call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
           tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
      tmb%can_use_transposed=.true.

      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
           tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%ovrlp)
  end if
  ! ###########################################################################



  !!! Calculate the derivative basis functions
  !!allocate(dphi(3*tmb%ham_descr%npsidim_orbs))
  !!call get_derivative_supportfunctions(tmb%ham_descr%npsidim_orbs, tmb%ham_descr%lzd%hgrids(1), tmb%ham_descr%lzd, tmb%orbs, tmb%ham_descr%psi, dphi)
  !!if(.not. associated(tmb%hpsi)) then
  !!    allocate(hpsit_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
  !!    allocate(hpsit_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
  !!    call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
  !!end if
  !!allocate(dpsit_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
  !!allocate(dpsit_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
  !!ist=1
  !!do i=1,3
  !!    call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, dphi(ist), dpsit_c, dpsit_f, tmb%ham_descr%lzd)
  !!    call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, dpsit_c, hpsit_c, dpsit_f, hpsit_f, tmb%linmat%ham)
  !!    noise(i)=0.d0
  !!    do iorb=1,tmb%orbs%norb
  !!       ii=matrixindex_in_compressed(tmb%linmat%ham,iorb,iorb)
  !!       noise(i) = noise(i) + tmb%linmat%ham%matrix_compr(ii)
  !!    end do
  !!    ist=ist+tmb%ham_descr%npsidim_orbs
  !!end do
  !!!!if (nproc>1) then
  !!!!   call mpiallred(noise(1), 3, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !!!!end if
  !!if (iproc==0) write(*,'(a,3es12.3)') 'noise', noise


  !if (iproc==0) write(*,*) 'correction_orthoconstraint',correction_orthoconstraint
  call orthoconstraintNonorthogonal(iproc, nproc, tmb%ham_descr%lzd, tmb%ham_descr%npsidim_orbs, tmb%ham_descr%npsidim_comp, &
       tmb%orbs, tmb%ham_descr%collcom, tmb%orthpar, correction_orthoconstraint, tmb%linmat, tmb%ham_descr%psi, tmb%hpsi, &
       tmb%linmat%ham, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, hpsit_c, hpsit_f, tmb%ham_descr%can_use_transposed, &
       overlap_calculated, experimental_mode)

  !!EXPERIMENTAL
  !!call calculate_residue_ks(iproc, nproc, num_extra, ksorbs, tmb, hpsit_c, hpsit_f)
  !!call calculate_residue_ks(iproc, nproc, 0, ksorbs, tmb, hpsit_c, hpsit_f)
  !!END EXPERIMENTAL

  call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
       tmb%orbs, tmb%hpsi, hpsi_small)




  !!! Gradient in the outer shell
  !!hh=(tmb%lzd%hgrids(1)+tmb%lzd%hgrids(2)+tmb%lzd%hgrids(3))/3.d0
  !!fnrm_in=0.d0
  !!fnrm_out=0.d0

  !!istart=0
  !!iold=0
  !!inew=0

  !!do iorb=1,tmb%orbs%norbp

  !!    iiorb=tmb%orbs%isorb+iorb
  !!    ilr=tmb%orbs%inwhichlocreg(iiorb)
  !!    ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

  !!    do iseg=1,tmb%lzd%llr(ilr)%wfd%nseg_c
  !!        jj=tmb%lzd%llr(ilr)%wfd%keyvloc(iseg)
  !!        j0=tmb%lzd%llr(ilr)%wfd%keygloc(1,iseg)
  !!        j1=tmb%lzd%llr(ilr)%wfd%keygloc(2,iseg)
  !!        ii=j0-1
  !!        i3=ii/((tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1))
  !!        ii=ii-i3*(tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1)
  !!        i2=ii/(tmb%lzd%llr(ilr)%d%n1+1)
  !!        i0=ii-i2*(tmb%lzd%llr(ilr)%d%n1+1)
  !!        i1=i0+j1-j0
  !!        do i=i0,i1
  !!            rx=(tmb%lzd%llr(ilr)%ns1+i)*tmb%lzd%hgrids(1)
  !!            ry=(tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)
  !!            rz=(tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)
  !!            rr = sqrt((rx-tmb%lzd%llr(ilr)%locregcenter(1))**2 + &
  !!                      (ry-tmb%lzd%llr(ilr)%locregcenter(2))**2 + &
  !!                      (rz-tmb%lzd%llr(ilr)%locregcenter(3))**2)
  !!            if (rr<tmb%lzd%llr(ilr)%locrad-8*hh) then
  !!                fnrm_in=fnrm_in+hpsi_small(istart+i-i0+jj)**2
  !!            else
  !!                fnrm_out=fnrm_out+hpsi_small(istart+i-i0+jj)**2
  !!            end if
  !!        end do
  !!    end do


  !!    isegf=tmb%lzd%llr(ilr)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr)%wfd%nseg_f)
  !!    do iseg=isegf,isegf+tmb%lzd%llr(ilr)%wfd%nseg_f-1
  !!        jj=tmb%lzd%llr(ilr)%wfd%keyvloc(iseg)
  !!        j0=tmb%lzd%llr(ilr)%wfd%keygloc(1,iseg)
  !!        j1=tmb%lzd%llr(ilr)%wfd%keygloc(2,iseg)
  !!        ii=j0-1
  !!        i3=ii/((tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1))
  !!        ii=ii-i3*(tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1)
  !!        i2=ii/(tmb%lzd%llr(ilr)%d%n1+1)
  !!        i0=ii-i2*(tmb%lzd%llr(ilr)%d%n1+1)
  !!        i1=i0+j1-j0
  !!        do i=i0,i1
  !!            rx=(tmb%lzd%llr(ilr)%ns1+i)*tmb%lzd%hgrids(1)
  !!            ry=(tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)
  !!            rz=(tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)
  !!            rr = sqrt((rx-tmb%lzd%llr(ilr)%locregcenter(1))**2 + &
  !!                      (ry-tmb%lzd%llr(ilr)%locregcenter(2))**2 + &
  !!                      (rz-tmb%lzd%llr(ilr)%locregcenter(3))**2)
  !!            if (rr<tmb%lzd%llr(ilr)%locrad-8*hh) then
  !!                fnrm_in=fnrm_in+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+1)**2
  !!                fnrm_in=fnrm_in+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+2)**2
  !!                fnrm_in=fnrm_in+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+3)**2
  !!                fnrm_in=fnrm_in+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+4)**2
  !!                fnrm_in=fnrm_in+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+5)**2
  !!                fnrm_in=fnrm_in+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+6)**2
  !!                fnrm_in=fnrm_in+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+7)**2
  !!            else
  !!                fnrm_out=fnrm_out+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+1)**2
  !!                fnrm_out=fnrm_out+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+2)**2
  !!                fnrm_out=fnrm_out+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+3)**2
  !!                fnrm_out=fnrm_out+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+4)**2
  !!                fnrm_out=fnrm_out+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+5)**2
  !!                fnrm_out=fnrm_out+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+6)**2
  !!                fnrm_out=fnrm_out+hpsi_small(istart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i-i0+jj-1)+7)**2
  !!            end if
  !!        end do
  !!    end do

  !!    istart=istart+ncount

  !!end do

  !!call mpiallred(fnrm_in, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !!call mpiallred(fnrm_out, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !!fnrm_tot2=fnrm_in+fnrm_out
  !!fnrm_in=sqrt(fnrm_in/dble(tmb%orbs%norb))
  !!fnrm_out=sqrt(fnrm_out/dble(tmb%orbs%norb))
  !!fnrm_tot2=sqrt(fnrm_tot2/dble(tmb%orbs%norb))
  !!if (iproc==0) write(*,'(a,3es16.4)') 'fnrm_in, fnrm_out, fnrm_tot2', fnrm_in, fnrm_out, fnrm_tot2





  if (present(hpsi_noprecond)) call dcopy(tmb%npsidim_orbs, hpsi_small, 1, hpsi_noprecond, 1)

  ! Calculate trace (or band structure energy, resp.)
  trH=0.d0
  call timing(iproc,'eglincomms','ON')
  do iorb=1,tmb%orbs%norb
     ii=matrixindex_in_compressed(tmb%linmat%ham,iorb,iorb)
     trH = trH + tmb%linmat%ham%matrix_compr(ii)
     !!if (iproc==0) write(*,*) 'iorb, value', iorb, tmb%linmat%ham%matrix_compr(ii)
  end do
  call timing(iproc,'eglincomms','OF')
  !call deallocate_sparseMatrix(lagmat,subname)

  ! trH is now the total energy (name is misleading, correct this)
  ! Multiply by 2 because when minimizing trace we don't have kernel
  !if(iproc==0)print *,'trH,energs',trH,energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
  if(tmb%orbs%nspin==1 .and. target_function/= TARGET_FUNCTION_IS_ENERGY) trH=2.d0*trH
  trH=trH-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp


  ! Cycle if the trace increased (steepest descent only)
  if(.not. ldiis%switchSD .and. ldiis%isx==0) then
      if(trH > ldiis%trmin+1.d-12*abs(ldiis%trmin)) then !1.d-12 is here to tolerate some noise...
          !!if(iproc==0) write(*,'(1x,a,es18.10,a,es18.10)') &
          !!    'WARNING: the target function is larger than its minimal value reached so far:',trH,' > ', ldiis%trmin
          if (iproc==0) then
              call yaml_newline()
              call yaml_warning('target function larger than its minimal value reached so far, &
                  &D='//trim(yaml_toa(trH-ldiis%trmin,fmt='(1es10.3)')))!//'. &
                  !&Decrease step size and restart with previous TMBs')
          end if
          !if(iproc==0) write(*,'(1x,a)') 'Decrease step size and restart with previous TMBs'
          energy_increased=.true.
      end if
  end if

  allocate(fnrmOvrlpArr(tmb%orbs%norb), stat=istat)
  call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)
  allocate(fnrmArr(tmb%orbs%norb), stat=istat)
  call memocc(istat, fnrmArr, 'fnrmArr', subname)

  ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
  ! of the previous iteration (fnrmOvrlpArr).
  call timing(iproc,'eglincomms','ON')
  ist=1
  fnrm_low=0.d0
  fnrm_high=0.d0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      if(it>1) then
         fnrmOvrlpArr(iorb)=ddot(ncount, hpsi_small(ist), 1, lhphiold(ist), 1)
         !fnrmOldArr(iorb)=ddot(ncount, lhphiold(ist), 1, lhphiold(ist), 1)
      end if
      fnrmArr(iorb)=ddot(ncount, hpsi_small(ist), 1, hpsi_small(ist), 1)
      ist=ist+ncount
      !!!if (mod(iiorb-1,9)+1<=4) then
      !!!    fnrm_low=fnrm_low+fnrmArr(iorb)
      !!!else
      !!!    fnrm_high=fnrm_high+fnrmArr(iorb)
      !!!end if
  end do

  !!!call mpiallred(fnrm_low, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !!!call mpiallred(fnrm_high, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !!!fnrm_low=sqrt(fnrm_low/(4.d0/9.d0*dble(tmb%orbs%norb)))
  !!!fnrm_high=sqrt(fnrm_high/(5.d0/9.d0*dble(tmb%orbs%norb)))
  !!!if (iproc==0) write(*,'(a,2es16.6)') 'fnrm_low, fnrm_high', fnrm_low, fnrm_high


  ! Determine the gradient norm and its maximal component. In addition, adapt the
  ! step size for the steepest descent minimization (depending on the angle 
  ! between the current gradient and the one from the previous iteration).
  ! This is of course only necessary if we are using steepest descent and not DIIS.
  ! if newgradient is true, the angle criterion cannot be used and the choice whether to
  ! decrease or increase the step size is only based on the fact whether the trace decreased or increased.

  ! TEMPORARY #################################
  ! This is just for tests, can be improved
  fnrmOvrlp_tot=0.d0
  fnrm_tot=0.d0
  fnrmOld_tot=0.d0
  do iorb=1,tmb%orbs%norbp
      if (it>1) fnrmOvrlp_tot=fnrmOvrlp_tot+fnrmOvrlpArr(iorb)
      fnrm_tot=fnrm_tot+fnrmArr(iorb)
      if (it>1) fnrmOld_tot=fnrmOld_tot+fnrmOldArr(iorb)
  end do
  if (it>1) call mpiallred(fnrmOvrlp_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  call mpiallred(fnrm_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  call mpiallred(fnrmOld_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  ! ###########################################

  fnrm=0.d0
  fnrmMax=0.d0
  do iorb=1,tmb%orbs%norbp
      fnrm=fnrm+fnrmArr(iorb)
      if(fnrmArr(iorb)>fnrmMax) fnrmMax=fnrmArr(iorb)
      if(it>1 .and. ldiis%isx==0 .and. .not.ldiis%switchSD) then
      ! Adapt step size for the steepest descent minimization.
          if (experimental_mode) then
              if (iproc==0 .and. iorb==1) then
                  !write(*,*) 'WARNING: USING SAME STEP SIZE'
                  call yaml_warning('Using same step size for all TMBs')
              end if
              tt=fnrmOvrlp_tot/sqrt(fnrm_tot*fnrmOld_tot)
          else
              tt=fnrmOvrlpArr(iorb)/sqrt(fnrmArr(iorb)*fnrmOldArr(iorb))
          end if
          ! apply thresholds so that alpha never goes below around 1.d-2 and above around 2
          if(tt>.6d0 .and. trH<trHold .and. alpha(iorb)<1.8d0) then
              alpha(iorb)=alpha(iorb)*1.1d0
          else if (alpha(iorb)>1.7d-3) then
              alpha(iorb)=alpha(iorb)*.6d0
          end if
          !!alpha(iorb)=min(alpha(iorb),1.5d0)
      end if
  end do
  call mpiallred(fnrm, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  call mpiallred(fnrmMax, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
  fnrm=sqrt(fnrm/dble(tmb%orbs%norb))
  fnrmMax=sqrt(fnrmMax)

  call dcopy(tmb%orbs%norb, fnrmArr(1), 1, fnrmOldArr(1), 1)

  iall=-product(shape(fnrmOvrlpArr))*kind(fnrmOvrlpArr)
  deallocate(fnrmOvrlpArr, stat=istat)
  call memocc(istat, iall, 'fnrmOvrlpArr', subname)

  iall=-product(shape(fnrmArr))*kind(fnrmArr)
  deallocate(fnrmArr, stat=istat)
  call memocc(istat, iall, 'fnrmArr', subname)

  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  alpha_mean=tt/dble(tmb%orbs%norb)
  alpha_max=maxval(alpha)
  call mpiallred(alpha_max, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)

  ! Copy the gradient (will be used in the next iteration to adapt the step size).
  call dcopy(tmb%npsidim_orbs, hpsi_small, 1, lhphiold, 1)
  call timing(iproc,'eglincomms','OF')

  ! if energy has increased or we only wanted to calculate the energy, not gradient, we can return here
  ! rather than calculating the preconditioning for nothing
  if ((energy_increased .or. energy_only) .and. target_function/=TARGET_FUNCTION_IS_HYBRID) return

  !!! Precondition the gradient.
  !!if(iproc==0) then
  !!    write(*,'(a)',advance='no') 'Preconditioning... '
  !!end if
 

  !if (target_function==TARGET_FUNCTION_IS_HYBRID) then
  !    allocate(hpsi_tmp(tmb%npsidim_orbs), stat=istat)
  !    call memocc(istat, hpsi_conf, 'hpsi_conf', subname)
  !end if
  !ist=1
  !do iorb=1,tmb%orbs%norbp
  !    iiorb=tmb%orbs%isorb+iorb
  !    ilr = tmb%orbs%inWhichLocreg(iiorb)
  !    ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  !    if(target_function==TARGET_FUNCTION_IS_HYBRID) then
  !        tt=ddot(ncnt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
  !        tt=tt/ddot(ncnt, hpsi_conf(ist), 1, hpsi_conf(ist), 1)
  !        do i=ist,ist+ncnt-1
  !            hpsi_tmp(i)=tt*hpsi_conf(i)
  !        end do
  !        call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
  !             tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
  !             nit_precond, hpsi_tmp(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
  !             tmb%confdatarr(iorb)%prefac, iorb, eval_zero)
  !        call daxpy(ncnt, -tt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
  !        call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
  !             tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
  !             nit_precond, hpsi_small(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
  !             0.d0, iorb, eval_zero)
  !        call daxpy(ncnt, 1.d0, hpsi_tmp(ist), 1, hpsi_small(ist), 1)
  !    else
  !    !    call choosePreconditioner2(iproc, nproc, tmb%orbs, tmb%lzd%llr(ilr), &
  !    !         tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
  !    !         nit_precond, hpsi_small(ist:ist+ncnt-1), tmb%confdatarr(iorb)%potorder, &
  !    !         tmb%confdatarr(iorb)%prefac, iorb, eval_zero)
  !    end if
  !    ist=ist+ncnt
  !end do

  !!if (iproc==0) write(*,*) 'HACK S.M.: precond'
  if(target_function==TARGET_FUNCTION_IS_HYBRID) then
     allocate(hpsi_tmp(tmb%npsidim_orbs), stat=istat)
     call memocc(istat, hpsi_tmp, 'hpsi_tmp', subname)
     ist=1
     do iorb=1,tmb%orbs%norbp
        iiorb=tmb%orbs%isorb+iorb
        ilr = tmb%orbs%inWhichLocreg(iiorb)
        ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

        tt=ddot(ncnt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
        tt=tt/ddot(ncnt, hpsi_conf(ist), 1, hpsi_conf(ist), 1)
        do i=ist,ist+ncnt-1
           hpsi_tmp(i)=tt*hpsi_conf(i)
        end do
        call daxpy(ncnt, -tt, hpsi_conf(ist), 1, hpsi_small(ist), 1)

        ist=ist+ncnt
     end do

     !!if (ldiis%isx>0) then
     !!    if (iproc==0) write(*,*) 'HACK precond: 10*conf'
     !!    tmb%confdatarr(:)%prefac=10.d0*tmb%confdatarr(:)%prefac
     !!end if
     !if (iproc==0) write(*,*) 'HACK precond: max(prefac,1.d-3)'
     !allocate(prefacarr(tmb%orbs%norbp))
     !prefacarr(:)=tmb%confdatarr(:)%prefac
     !tmb%confdatarr(:)%prefac=max(tmb%confdatarr(:)%prefac,1.d-3)
     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_tmp,tmb%confdatarr,gnrm,gnrm_zero)
     !tmb%confdatarr(:)%prefac=prefacarr(:)
     !deallocate(prefacarr)
     !!if (ldiis%isx>0) then
     !!    tmb%confdatarr(:)%prefac=0.1d0*tmb%confdatarr(:)%prefac
     !!end if

     ! temporarily turn confining potential off...
     allocate(prefac(tmb%orbs%norbp),stat=istat)
     call memocc(istat, prefac, 'prefac', subname)
     prefac(:)=tmb%confdatarr(:)%prefac
     tmb%confdatarr(:)%prefac=0.0d0
     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero) ! prefac should be zero
     call daxpy(tmb%npsidim_orbs, 1.d0, hpsi_tmp(1), 1, hpsi_small(1), 1)
     ! ...revert back to correct value
     tmb%confdatarr(:)%prefac=prefac

     iall=-product(shape(prefac))*kind(prefac)
     deallocate(prefac, stat=istat)
     call memocc(istat, iall, 'prefac', subname)
     iall=-product(shape(hpsi_conf))*kind(hpsi_conf)
     deallocate(hpsi_conf, stat=istat)
     call memocc(istat, iall, 'hpsi_conf', subname)
     iall=-product(shape(hpsi_tmp))*kind(hpsi_tmp)
     deallocate(hpsi_tmp, stat=istat)
     call memocc(istat, iall, 'hpsi_tmp', subname)
  else
     !!if (ldiis%isx>0) then
         !if (iproc==0) write(*,*) 'HACK precond: max(prefac,1.d-4)'
         !tmb%confdatarr(:)%prefac=10.d0*tmb%confdatarr(:)%prefac
         !allocate(prefacarr(tmb%orbs%norbp))
         !prefacarr(:)=tmb%confdatarr(:)%prefac
         !tmb%confdatarr(:)%prefac=max(tmb%confdatarr(:)%prefac,1.d-4)
     !!end if
     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero)
     !!if (ldiis%isx>0) then
         !tmb%confdatarr(:)%prefac=0.1d0*tmb%confdatarr(:)%prefac
         !tmb%confdatarr(:)%prefac=prefacarr(:)
         !deallocate(prefacarr)
     !!end if
  end if

  if (iproc==0) then
      call yaml_map('Preconditioning',.true.)
  end if


  !sum over all the partial residues
  if (nproc > 1) then
      garray(1)=gnrm
      garray(2)=gnrm_zero
      call mpiallred(garray(1),2,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
      gnrm     =garray(1)
      gnrm_zero=garray(2)
  end if

  !if (target_function==TARGET_FUNCTION_IS_HYBRID) then
  !    iall=-product(shape(hpsi_conf))*kind(hpsi_conf)
  !    deallocate(hpsi_conf, stat=istat)
  !    call memocc(istat, iall, 'hpsi_conf', subname)
  !    iall=-product(shape(hpsi_tmp))*kind(hpsi_tmp)
  !    deallocate(hpsi_tmp, stat=istat)
  !    call memocc(istat, iall, 'hpsi_tmp', subname)
  !end if

  !!if(iproc==0) then
  !!    write(*,'(a)') 'done.'
  !!end if

  call timing(iproc,'eglincomms','ON')
  ist=1
  gnrm_old=gnrm
  gnrm=0.d0
  gnrmMax=0.d0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      gnrmArr=ddot(ncount, hpsi_small(ist), 1, hpsi_small(ist), 1)
      gnrm=gnrm+gnrmArr
      if(gnrmArr>gnrmMax) gnrmMax=gnrmArr
      ist=ist+ncount
  end do

  call mpiallred(gnrm, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  call mpiallred(gnrmMax, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
  gnrm=sqrt(gnrm/dble(tmb%orbs%norb))
  gnrmMax=sqrt(gnrmMax)
  !if (iproc==0) write(*,'(a,3es16.6)') 'AFTER: gnrm, gnrmmax, gnrm/gnrm_old',gnrm,gnrmmax,gnrm/gnrm_old
  call timing(iproc,'eglincomms','OF')


end subroutine calculate_energy_and_gradient_linear

subroutine calculate_residue_ks(iproc, nproc, num_extra, ksorbs, tmb, hpsit_c, hpsit_f)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, num_extra
  type(dft_wavefunction), intent(in) :: tmb
  type(orbitals_data), intent(in) :: ksorbs
  real(kind=8),dimension(:),pointer :: hpsit_c, hpsit_f

  integer :: iorb, jorb, istat, iall, ierr
  real(kind=8) :: ksres_sum
  real(kind=8), dimension(:), allocatable :: ksres
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, grad_coeff
  type(sparseMatrix) :: grad_ovrlp
  character(len=256) :: subname='calculate_residue_ks'


  ! want to calculate the residue of the KS states here, not just the tmbs
  ! for now just occupied, eventually would want occupied+num_extra
  ! probably better to calculate c_i^a|g_a> first but copying and pasting for now (INEFFICIENT but just testing)
  !!if(associated(hpsit_c)) then
  !!    iall=-product(shape(hpsit_c))*kind(hpsit_c)
  !!    deallocate(hpsit_c, stat=istat)
  !!    call memocc(istat, iall, 'hpsit_c', subname)
  !!end if
  !!if(associated(hpsit_f)) then
  !!    iall=-product(shape(hpsit_f))*kind(hpsit_f)
  !!    deallocate(hpsit_f, stat=istat)
  !!    call memocc(istat, iall, 'hpsit_f', subname)
  !!end if
  !!allocate(hpsit_c(sum(collcom%nrecvcounts_c)), stat=istat)
  !!call memocc(istat, hpsit_c, 'hpsit_c', subname)
  !!allocate(hpsit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
  !!call memocc(istat, hpsit_f, 'hpsit_f', subname)

  ! should already be done in orthoconstraintnonorthogonal so can use directly
  !call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
  !     tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
  !!can_use_transposed=.true.

  call nullify_sparsematrix(grad_ovrlp)
  call sparse_copy_pattern(tmb%linmat%ham, grad_ovrlp, iproc, subname)
  allocate(grad_ovrlp%matrix_compr(grad_ovrlp%nvctr), stat=istat)
  call memocc(istat, grad_ovrlp%matrix_compr, 'grad_ovrlp%matrix_compr', subname)
  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, hpsit_c, hpsit_c, &
       hpsit_f, hpsit_f, grad_ovrlp)

  allocate(grad_coeff(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, grad_coeff, 'grad_coeff', subname)

  allocate(coeff_tmp(tmb%orbs%norbp,max(tmb%orbs%norb,1)), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)

  allocate(grad_ovrlp%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, grad_ovrlp%matrix, 'grad_ovrlp%matrix', subname)
  call uncompressMatrix(iproc,grad_ovrlp)

  ! can change this so only go up to ksorbs%norb...
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, 1.d0, grad_ovrlp%matrix(tmb%orbs%isorb+1,1), &
          tmb%orbs%norb, tmb%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp, tmb%orbs%norbp)
     call dgemm('t', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norbp, 1.d0, tmb%coeff(tmb%orbs%isorb+1,1), &
          tmb%orbs%norb, coeff_tmp, tmb%orbs%norbp, 0.d0, grad_coeff, tmb%orbs%norb)
  else
     call to_zero(tmb%orbs%norb**2, grad_coeff(1,1))
  end if

  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=istat)
  call memocc(istat,iall,'coeff_tmp',subname)

  if (nproc>1) then
      call mpiallred(grad_coeff(1,1), tmb%orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  ! now calculate sqrt(<g_i|g_i>) and mean value
  allocate(ksres(ksorbs%norb+num_extra), stat=istat)
  call memocc(istat, ksres, 'ksres', subname)
  
  ksres_sum=0.0d0
  do iorb=1,ksorbs%norb+num_extra
    ksres(iorb)=dsqrt(grad_coeff(iorb,iorb))
    !ksres_sum=ksres_sum+ksres(iorb)
    ksres_sum=ksres_sum+grad_coeff(iorb,iorb)
    if (iproc==0) write(*,*) 'KS residue',iorb,ksres(iorb)!,tmb%orbs%occup(iorb)
  end do
  if (iproc==0) write(*,*) 'Average KS residue',sqrt(ksres_sum/real(ksorbs%norb+num_extra,gp))


  !call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
  !     tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ham)
  ! calculate Tr[Kg]  (not recalculating kernel as don't have the correct occs here)
  !call calculate_kernel_and_energy(iproc,nproc,denskern,grad_coeff,ksres_sum,tmb%coeff,orbs,tmb%orbs,.true.)
  call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,grad_ovrlp,ksres_sum,tmb%coeff,tmb%orbs,tmb%orbs,.false.)
  if (iproc==0) write(*,*) 'KS residue from trace',dsqrt(ksres_sum)/real(tmb%orbs%norb,gp) ! should update normalization as would only be occ here not extra?

  call deallocate_sparseMatrix(grad_ovrlp, subname)

  iall=-product(shape(grad_coeff))*kind(grad_coeff)
  deallocate(grad_coeff,stat=istat)
  call memocc(istat,iall,'grad_coeff',subname)

  iall=-product(shape(ksres))*kind(ksres)
  deallocate(ksres,stat=istat)
  call memocc(istat,iall,'ksres',subname)


end subroutine calculate_residue_ks



subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb,  &
           lphiold, alpha, trH, alpha_mean, alpha_max, alphaDIIS, hpsi_small, ortho, psidiff, &
           experimental_mode)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => hpsitopsi_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  type(localizedDIISParameters), intent(inout) :: ldiis
  type(DFT_wavefunction), target,intent(inout) :: tmb
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: lphiold
  real(kind=8), intent(in) :: trH, alpha_mean, alpha_max
  real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha, alphaDIIS
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: hpsi_small
  real(kind=8), dimension(tmb%npsidim_orbs), optional,intent(out) :: psidiff
  logical, intent(in) :: ortho, experimental_mode
  
  ! Local variables
  integer :: istat, iall, i
  character(len=*), parameter :: subname='hpsitopsi_linear'
  real(kind=8),dimension(:),allocatable :: psittmp_c, psittmp_f
  real(kind=8), dimension(:),allocatable :: norm
  real(kind=8) :: ddot

  call DIISorSD(iproc, it, trH, tmb, ldiis, alpha, alphaDIIS, lphiold)
  if(iproc==0) then
      call yaml_newline()
      call yaml_open_map('Optimization',flow=.true.)
      if(ldiis%isx>0) then
          !!write(*,'(1x,3(a,i0))') 'DIIS informations: history length=',ldiis%isx, ', consecutive failures=', &
          !!    ldiis%icountDIISFailureCons, ', total failures=', ldiis%icountDIISFailureTot
          call yaml_map('algorithm','DIIS')
          call yaml_map('history length',ldiis%isx)
          call yaml_map('consecutive failures',ldiis%icountDIISFailureCons)
          call yaml_map('total failures',ldiis%icountDIISFailureTot)
      else
          !!write(*,'(1x,2(a,es9.3),a,i0)') 'SD informations: mean alpha=', alpha_mean, ', max alpha=', alpha_max,&
          !!', consecutive successes=', ldiis%icountSDSatur
          call yaml_map('algorithm','SD')
          call yaml_map('mean alpha',alpha_mean,fmt='(es9.3)')
          call yaml_map('max alpha',alpha_max,fmt='(es9.3)')
          call yaml_map('consecutive successes',ldiis%icountSDSatur)
      end if
      call yaml_close_map()
      call yaml_newline()
  end if

  ! Improve the orbitals, depending on the choice made above.
  if (present(psidiff)) call dcopy(tmb%npsidim_orbs, tmb%psi, 1, psidiff, 1)
  if(.not.ldiis%switchSD) then
      call improveOrbitals(iproc, tmb, ldiis, alpha, hpsi_small, experimental_mode)
  else
      !if(iproc==0) write(*,'(1x,a)') 'no improvement of the orbitals, recalculate gradient'
      if (iproc==0) then
          call yaml_warning('no improvement of the orbitals, recalculate gradient')
          call yaml_newline()
      end if
  end if
  if (present(psidiff)) then
      do i=1,tmb%npsidim_orbs
          psidiff(i)=tmb%psi(i)-psidiff(i)
      end do 
  end if

  ! The transposed quantities can now not be used any more...
  if(tmb%can_use_transposed) then
      iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
      deallocate(tmb%psit_c, stat=istat)
      call memocc(istat, iall, 'tmb%psit_c', subname)
      iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
      deallocate(tmb%psit_f, stat=istat)
      call memocc(istat, iall, 'tmb%psit_f', subname)
      tmb%can_use_transposed=.false.
  end if


  !!!!  ! EXPERIMENTAL -- DON'T USE IT ########################################################
!!!!
!!!!  write(*,*) 'before: iproc, ddot', iproc, ddot(tmb%npsidim_orbs, tmb%psi, 1, tmb%psi, 1)
!!!!
!!!!  call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
!!!!       tmb%orbs, tmb%psi, tmb%ham_descr%psi)
!!!!
!!!!  !if (.not.tmb%ham_descr%can_use_transposed) then
!!!!  if (.not.associated(tmb%ham_descr%psit_c)) then
!!!!      allocate(tmb%ham_descr%psit_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
!!!!      call memocc(istat, tmb%ham_descr%psit_c, 'tmb%ham_descr%psit_c', subname)
!!!!  end if
!!!!  if (.not.associated(tmb%ham_descr%psit_f)) then
!!!!      allocate(tmb%ham_descr%psit_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
!!!!      call memocc(istat, tmb%ham_descr%psit_f, 'tmb%ham_descr%psit_f', subname)
!!!!  end if
!!!!  if (.not.associated(tmb%ham_descr%psi)) then
!!!!      write(*,*) 'ALLOCATE'
!!!!      allocate(tmb%ham_descr%psi(tmb%ham_descr%npsidim_orbs), stat=istat)
!!!!      call memocc(istat, tmb%ham_descr%psi, 'tmb%ham_descr%psi', subname)
!!!!  end if
!!!!  allocate(psittmp_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
!!!!  call memocc(istat, psittmp_c, 'psittmp_c', subname)
!!!!  allocate(psittmp_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
!!!!  call memocc(istat, psittmp_f, 'psittmp_f', subname)
!!!!
!!!!  call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
!!!!       tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd)
!!!!  call dcopy(sum(tmb%ham_descr%collcom%nrecvcounts_c), tmb%ham_descr%psit_c, 1, psittmp_c, 1)
!!!!  call dcopy(7*sum(tmb%ham_descr%collcom%nrecvcounts_f), tmb%ham_descr%psit_f, 1, psittmp_f, 1)
!!!!  !tmb%linmat%denskern%matrix_compr=1.d0
!!!!  write(*,*) 'iproc, ddot trans c', iproc, &
!!!!      ddot(sum(tmb%ham_descr%collcom%nrecvcounts_c), tmb%ham_descr%psit_c, 1, psittmp_c, 1)
!!!!  write(*,*) 'iproc, ddot trans f', iproc, &
!!!!      ddot(7*sum(tmb%ham_descr%collcom%nrecvcounts_f), tmb%ham_descr%psit_f, 1, psittmp_f, 1)
!!!!  call build_linear_combination_transposed(tmb%ham_descr%collcom, &
!!!!       tmb%linmat%denskern, psittmp_c, psittmp_f, .true., tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, iproc)
!!!!  write(*,*) 'iproc, ddot trans c after', iproc, &
!!!!      ddot(sum(tmb%ham_descr%collcom%nrecvcounts_c), tmb%ham_descr%psit_c, 1, psittmp_c, 1)
!!!!  write(*,*) 'iproc, ddot trans f after', iproc, &
!!!!      ddot(7*sum(tmb%ham_descr%collcom%nrecvcounts_f), tmb%ham_descr%psit_f, 1, psittmp_f, 1)
!!!! !!call build_linear_combination_transposed(tmb%ham_descr%collcom, &
!!!! !!     tmb%linmat%denskern, hpsittmp_c, hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
!!!!  call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
!!!!       tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%psi, tmb%ham_descr%lzd)
!!!!  write(*,*) 'after untranspose: iproc, ddot', iproc, ddot(tmb%ham_descr%npsidim_orbs, tmb%ham_descr%psi, 1, tmb%ham_descr%psi, 1)
!!!!
!!!!          !!call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
!!!!          !!     tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
!!!!          !!call build_linear_combination_transposed(tmb%ham_descr%collcom, &
!!!!          !!     tmb%linmat%denskern, hpsittmp_c, hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
!!!!
!!!!  call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
!!!!       tmb%orbs, tmb%ham_descr%psi, tmb%psi)
!!!!  write(*,*) 'after: iproc, ddot', iproc, ddot(tmb%npsidim_orbs, tmb%psi, 1, tmb%psi, 1)
!!!!
!!!!  iall=-product(shape(tmb%ham_descr%psit_c))*kind(tmb%ham_descr%psit_c)
!!!!  deallocate(tmb%ham_descr%psit_c, stat=istat)
!!!!  call memocc(istat, iall, 'tmb%ham_descr%psit_c', subname)
!!!!  iall=-product(shape(tmb%ham_descr%psit_f))*kind(tmb%ham_descr%psit_f)
!!!!  deallocate(tmb%ham_descr%psit_f, stat=istat)
!!!!  call memocc(istat, iall, 'tmb%ham_descr%psit_f', subname)
!!!!  iall=-product(shape(psittmp_c))*kind(psittmp_c)
!!!!  deallocate(psittmp_c, stat=istat)
!!!!  call memocc(istat, iall, 'psittmp_c', subname)
!!!!  iall=-product(shape(psittmp_f))*kind(psittmp_f)
!!!!  deallocate(psittmp_f, stat=istat)
!!!!  call memocc(istat, iall, 'psittmp_f', subname)
!!!!  tmb%ham_descr%can_use_transposed=.false.
!!!!  ! END EXPERIMENTAL ###################################################################

  if (.not.ortho .and. iproc==0) then
      call yaml_map('Orthogonalization',.false.)
  end if

  if(.not.ldiis%switchSD.and.ortho) then
      !!if(iproc==0) then
      !!     write(*,'(1x,a)',advance='no') 'Orthonormalization... '
      !!end if

      ! CHEATING here and passing tmb%linmat%denskern instead of tmb%linmat%inv_ovrlp
      call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
           tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, &
           tmb%can_use_transposed)
      if (iproc == 0) then
          call yaml_map('Orthogonalization',.true.)
      end if
  else if (experimental_mode) then
      ! Wasteful to do it transposed...
      !if (iproc==0) write(*,*) 'normalize...'
      if (iproc==0) call yaml_map('normalization',.true.)
      if(associated(tmb%psit_c)) then
          iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
          deallocate(tmb%psit_c, stat=istat)
          call memocc(istat, iall, 'tmb%psit_c', subname)
      end if
      if(associated(tmb%psit_f)) then
          iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
          deallocate(tmb%psit_f, stat=istat)
          call memocc(istat, iall, 'tmb%psit_f', subname)
      end if
      allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
      allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
      tmb%can_use_transposed=.true.

      call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
      allocate(norm(tmb%orbs%norb))
      call normalize_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, tmb%psit_f, norm)
      deallocate(norm)
      call untranspose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, tmb%psit_c, tmb%psit_f, tmb%psi, tmb%lzd)
      if (iproc == 0) then
          call yaml_map('Normalization',.true.)
      end if
  end if

  ! Emit that new wavefunctions are ready.
  if (tmb%c_obj /= 0) then
     call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
  end if

end subroutine hpsitopsi_linear

