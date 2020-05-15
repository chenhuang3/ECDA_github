!> @file
!! Optimize the coefficients
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!
!! @note
!!  Coefficients are defined for Ntmb KS orbitals so as to maximize the number
!!  of orthonormality constraints. This should speedup the convergence by
!!  reducing the effective number of degrees of freedom.
subroutine optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, fnrm_crit, itmax, energy, sd_fit_curve, &
    factor, itout, it_scc, it_cdft, reorder, num_extra)
  use module_base
  use module_types
  use module_interfaces, fake_name => optimize_coeffs
  use diis_sd_optimization
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, itmax, itout, it_scc, it_cdft
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  type(DIIS_obj), intent(inout) :: ldiis_coeff
  real(kind=gp),intent(in):: fnrm_crit
  real(kind=gp),intent(out):: fnrm
  real(kind=gp), intent(inout) :: energy
  logical, intent(in) :: sd_fit_curve
  real(kind=gp), intent(in) :: factor
  integer, optional, intent(in) :: num_extra
  logical, optional, intent(in) :: reorder

  ! Local variables
  integer:: iorb, jorb, iiorb, ierr, it, itlast
  real(kind=gp),dimension(:,:),allocatable:: grad, grad_cov_or_coeffp !coeffp, grad_cov
  real(kind=gp) :: tt, ddot, energy0, pred_e

  call f_routine(id='optimize_coeffs')

  if (ldiis_coeff%idsx == 0 .and. sd_fit_curve) then
     ! calculate initial energy for SD line fitting and printing (maybe don't need to (re)calculate kernel here?)
     call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,tmb%linmat%ham,energy0,&
          tmb%coeff,orbs,tmb%orbs,.true.)
  else
     energy0=energy
  end if

  if (present(num_extra)) then
     grad=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='grad')
     grad_cov_or_coeffp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='grad_cov_or_coeffp')
  else
     grad=f_malloc((/tmb%orbs%norb,orbs%norbp/), id='grad')
     grad_cov_or_coeffp=f_malloc((/tmb%orbs%norb,orbs%norbp/), id='grad_cov_or_coeffp')
  end if

  if (iproc==0) then
      call yamL_newline()
      call yaml_open_sequence('expansion coefficients optimization',label=&
           'it_coeff'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'_'//&
           trim(adjustl(yaml_toa(it_cdft,fmt='(i3.3)')))//&
           '_'//trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)'))))
  end if


  do it=1,itmax

      if (iproc==0) then
          call yaml_newline()
          call yaml_sequence(advance='no')
          call yaml_open_map(flow=.true.)
          call yaml_comment('it coeff:'//yaml_toa(it,fmt='(i6)'),hfill='-')
      end if

     if (present(num_extra)) then
        call calculate_coeff_gradient_extra(iproc,nproc,num_extra,tmb,orbs,grad_cov_or_coeffp(1,1),grad(1,1))
     else
        call calculate_coeff_gradient(iproc,nproc,tmb,orbs,grad_cov_or_coeffp,grad)
     end if

     ! Precondition the gradient (only making things worse...)
     !call precondition_gradient_coeff(tmb%orbs%norb, orbs%norbp, tmb%linmat%ham%matrix, tmb%linmat%ovrlp%matrix, grad)

     call timing(iproc,'dirmin_sddiis','ON')

     !For fnrm, we only sum on the occupied KS orbitals
     tt=0.d0
     if (present(num_extra)) then
        do iorb=1,tmb%orbs%norbp
            tt=tt+ddot(tmb%orbs%norb, grad_cov_or_coeffp(1,iorb), 1, grad(1,iorb), 1)
        end do
     else
        do iorb=1,orbs%norbp
            tt=tt+ddot(tmb%orbs%norb, grad_cov_or_coeffp(1,iorb), 1, grad(1,iorb), 1)
        end do
     end if
     call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
     fnrm=2.0_gp*tt

     !scale the gradient (not sure if we always want this or just fragments/constrained!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
     if (present(num_extra)) then
        call dscal(tmb%orbs%norb*tmb%orbs%norbp,factor,grad,1)
     else
        call dscal(tmb%orbs%norb*orbs%norbp,factor,grad,1)
     end if

     if (ldiis_coeff%idsx > 0) then !do DIIS
        !TO DO: make sure DIIS works
        ldiis_coeff%mids=mod(ldiis_coeff%ids,ldiis_coeff%idsx)+1
        ldiis_coeff%ids=ldiis_coeff%ids+1

        if (present(num_extra)) then
           call dcopy(tmb%orbs%norb*tmb%orbs%norbp,tmb%coeff(1,tmb%orbs%isorb+1),1,grad_cov_or_coeffp,1)

           call diis_opt(iproc,nproc,1,0,1,(/iproc/),(/tmb%orbs%norb*tmb%orbs%norbp/),tmb%orbs%norb*tmb%orbs%norbp,&
                grad_cov_or_coeffp,grad,ldiis_coeff) 
        else
           call dcopy(tmb%orbs%norb*orbs%norbp,tmb%coeff(1,orbs%isorb+1),1,grad_cov_or_coeffp,1)

           call diis_opt(iproc,nproc,1,0,1,(/iproc/),(/tmb%orbs%norb*orbs%norbp/),tmb%orbs%norb*orbs%norbp,&
                grad_cov_or_coeffp,grad,ldiis_coeff) 
        end if
     else  !steepest descent with curve fitting for line minimization
        call timing(iproc,'dirmin_sddiis','OF')
        if (present(num_extra)) then   
           if (sd_fit_curve) call find_alpha_sd(iproc,nproc,ldiis_coeff%alpha_coeff,tmb,tmb%orbs,&
                grad_cov_or_coeffp,grad,energy0,fnrm,pred_e)
           call timing(iproc,'dirmin_sddiis','ON')
           do iorb=1,tmb%orbs%norbp
              iiorb = tmb%orbs%isorb + iorb
              do jorb=1,tmb%orbs%norb
                 grad_cov_or_coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-ldiis_coeff%alpha_coeff*grad(jorb,iorb)
              end do
           end do
        else
           if (sd_fit_curve) call find_alpha_sd(iproc,nproc,ldiis_coeff%alpha_coeff,tmb,orbs,&
                grad_cov_or_coeffp,grad,energy0,fnrm,pred_e)
           call timing(iproc,'dirmin_sddiis','ON')
           do iorb=1,orbs%norbp
              iiorb = orbs%isorb + iorb
              do jorb=1,tmb%orbs%norb
                 grad_cov_or_coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-ldiis_coeff%alpha_coeff*grad(jorb,iorb)
              end do
           end do
        end if
     end if

     call timing(iproc,'dirmin_sddiis','OF')

     call timing(iproc,'dirmin_allgat','ON')
     if (present(num_extra)) then  
        if(nproc > 1) then 
           call mpi_allgatherv(grad_cov_or_coeffp, tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, tmb%coeff, &
              tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, mpi_double_precision, &
              bigdft_mpi%mpi_comm, ierr)
        else
           call dcopy(tmb%orbs%norb*tmb%orbs%norb,grad_cov_or_coeffp(1,1),1,tmb%coeff(1,1),1)
        end if

        call timing(iproc,'dirmin_allgat','OF')

        fnrm=sqrt(fnrm/dble(orbs%norb+num_extra))
     else 
        if(nproc > 1) then 
           call mpi_allgatherv(grad_cov_or_coeffp, tmb%orbs%norb*orbs%norbp, mpi_double_precision, tmb%coeff, &
              tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, &
              bigdft_mpi%mpi_comm, ierr)
        else
           call dcopy(tmb%orbs%norb*orbs%norb,grad_cov_or_coeffp(1,1),1,tmb%coeff(1,1),1)
        end if

        call timing(iproc,'dirmin_allgat','OF')

        fnrm=sqrt(fnrm/dble(orbs%norb))
     end if

     !! experimenting with calculating cHc and diagonalizing
     !if (present(num_extra).and.present(reorder)) then
     !   call reorthonormalize_coeff(iproc, nproc, orbs%norb+num_extra, -8, -8, 0, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
     !   !call reorthonormalize_coeff(iproc, nproc, orbs%norb+num_extra, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
     !   call reordering_coeffs(iproc, nproc, num_extra, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
     !   if (reorder) call reordering_coeffs(iproc, nproc, num_extra, orbs, tmb%orbs, &
     !        tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
     !else if (present(reorder)) then
     !   call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 0, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
     !   !call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
     !   call reordering_coeffs(iproc, nproc, 0, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
     !   if (reorder) call reordering_coeffs(iproc, nproc, 0, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
     !else
     !   call reordering_coeffs(iproc, nproc, 0, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, .false.)
     !end if

     ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
     ! instead of twice could add some criterion to check accuracy?
     if (present(num_extra)) then
        call reorthonormalize_coeff(iproc, nproc, orbs%norb+num_extra, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
        call reorthonormalize_coeff(iproc, nproc, orbs%norb+num_extra, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
     else
        call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
        call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
     end if
     !!!!!!!!!!!!!!!!!!!!!!!!
     !can't put coeffs directly in ksorbs%eval as intent in, so change after - problem with orthonormality of coeffs so adding extra
     !call find_eval_from_coeffs(iproc, nproc, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, &
     !     tmb%coeff, tmb%orbs%eval, .true., .true.)
     !call order_coeffs_by_energy(orbs%norb,tmb%orbs%norb,tmb%coeff,tmb%orbs%eval)
     !!!!!!!!!!!!!!!!!!!!!!!!

     call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,tmb%linmat%ham,energy,&
          tmb%coeff,orbs,tmb%orbs,.true.)
     !write(127,*) ldiis_coeff%alpha_coeff,energy
     !close(127)

     ! can't check ebs only, need to check Etot in linearScaling, but if we do it>1 without curve fitting will still need to update here
     !if (ldiis_coeff%idsx == 0 .and. (.not. sd_fit_curve) .and. energy0/=0.0_gp) then ! only update alpha after first iteration
     !   ! apply a cap so that alpha_coeff never goes below around 1.d-2 or above 2
     !   if ((energy-energy0)<0.d0 .and. ldiis_coeff%alpha_coeff < 1.8d0) then
     !      ldiis_coeff%alpha_coeff=1.1d0*ldiis_coeff%alpha_coeff
     !   else if (ldiis_coeff%alpha_coeff > 1.7d-3) then
     !      ldiis_coeff%alpha_coeff=0.5d0*ldiis_coeff%alpha_coeff
     !   end if
     !   if (iproc==0) print*,'EBSdiff,alpha',energy-energy0,ldiis_coeff%alpha_coeff,energy,energy0
     !end if

     if (iproc==0) write(*,*) ''
     if (sd_fit_curve .and. ldiis_coeff%idsx == 0) then
        !!if (iproc==0) write(*,'(a,I4,2x,6(ES16.6e3,2x))')'DminSD: it, fnrm, ebs, ebsdiff, alpha, pred E, diff',&
        !!     it,fnrm,energy0,energy-energy0,ldiis_coeff%alpha_coeff,pred_e,pred_e-energy
        if (iproc==0) then
            call yaml_map('method','DminSD')
            call yaml_newline()
            call yaml_map('iter',it)
            call yaml_map('fnrm',fnrm,fmt='(es9.2)')
            call yaml_map('eBS',energy0,fmt='(es24.17)')
            call yaml_map('D',energy-energy0,fmt='(es10.3)')
            call yaml_map('alpha',ldiis_coeff%alpha_coeff,fmt='(es10.3)')
            call yaml_map('predicted energy',pred_e,fmt='(es24.17)')
            call yaml_map('D',pred_e-energy,fmt='(es10.3)')
        end if
     else if (ldiis_coeff%idsx == 0) then
        !!if (iproc==0) write(*,'(a,I4,2x,4(ES16.6e3,2x))')'DminSD: it, fnrm, ebs, ebsdiff, alpha',&
        !!     it,fnrm,energy0,energy-energy0,ldiis_coeff%alpha_coeff
        if (iproc==0) then
            call yaml_map('method','DminSD')
            call yaml_newline()
            call yaml_map('iter',it)
            call yaml_map('fnrm',fnrm,fmt='(es9.2)')
            call yaml_map('eBS',energy0,fmt='(es24.17)')
            call yaml_map('D',energy-energy0,fmt='(es10.3)')
            call yaml_map('alpha',ldiis_coeff%alpha_coeff,fmt='(es10.3)')
        end if
     else
        !!if (iproc==0) write(*,'(a,I4,2x,3(ES16.6e3,2x))')'DminDIIS: it, fnrm, ebs, ebsdiff',&
        !!     it,fnrm,energy0,energy-energy0
        if (iproc==0) then
            call yaml_map('method','DminDIIS')
            call yaml_newline()
            call yaml_map('iter',it)
            call yaml_map('fnrm',fnrm,fmt='(es9.2)')
            call yaml_map('eBS',energy0,fmt='(es24.17)')
            call yaml_map('D',energy-energy0,fmt='(es10.3)')
        end if
     end if

     energy0=energy

     if (iproc==0) then
         call yaml_close_map()
         call bigdft_utils_flush(unit=6)
     end if


     itlast=it
     if (fnrm<fnrm_crit) exit

  end do

  !!if (iproc==0) then
  !!    call yaml_newline()
  !!    call yaml_sequence(label='final_coeff'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'_'//&
  !!         trim(adjustl(yaml_toa(it_cdft,fmt='(i3.3)')))//'_'//&
  !!         trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)'))),advance='no')
  !!    call yaml_open_map(flow=.true.)
  !!    call yaml_comment('iter:'//yaml_toa(itlast,fmt='(i6)'),hfill='-')
  !!    call yaml_map('iter',itlast,fmt='(i6)')
  !!    call yaml_map('fnrm',fnrm,fmt='(es9.2)')
  !!    call yaml_map('eBS',energy0,fmt='(es24.17)')
  !!    call yaml_map('D',energy-energy0,fmt='(es10.3)')
  !!    call yaml_close_map()
  !!    call bigdft_utils_flush(unit=6)
  !!end if


  if (iproc==0) then
      call yaml_close_sequence()
      call yaml_newline()
  end if


  call f_free(grad_cov_or_coeffp)
  call f_free(grad)

  call f_release_routine()

end subroutine optimize_coeffs

! subset of reordering coeffs - need to arrange this routines better but taking the lazy route for now
! (also assuming we have no extra - or rather number of extra bands come from input.mix not input.lin)
subroutine coeff_weight_analysis(iproc, nproc, input, ksorbs, tmb, ref_frags)
  use module_base
  use module_types
  use module_interfaces
  use module_fragments
  use constrained_dft
  use yaml_output
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc
  type(orbitals_data), intent(in) :: ksorbs
  type(dft_wavefunction), intent(inout) :: tmb
  type(input_variables),intent(in) :: input
  type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags

  integer :: iorb, jorb, istat, iall, ierr, itmb, jtmb, ifrag
  integer, dimension(2) :: ifrag_charged
  real(kind=8), dimension(:,:,:), allocatable :: weight_coeff
  type(sparseMatrix) :: weight_matrix
  character(len=256) :: subname='coeff_weight_analysis'

  call nullify_sparsematrix(weight_matrix)
  call sparse_copy_pattern(tmb%linmat%ham, weight_matrix, iproc, subname)
  allocate(weight_matrix%matrix_compr(weight_matrix%nvctr), stat=istat)
  call memocc(istat, weight_matrix%matrix_compr, 'weight_matrix%matrix_compr', subname)

  allocate(weight_coeff(ksorbs%norb,ksorbs%norb,input%frag%nfrag), stat=istat)
  call memocc(istat, weight_coeff, 'weight_coeff', subname)

  do ifrag=1,input%frag%nfrag
     ifrag_charged(1)=ifrag
     call calculate_weight_matrix_lowdin(weight_matrix,1,ifrag_charged,tmb,input,ref_frags,.false.)
     allocate(weight_matrix%matrix(weight_matrix%full_dim1,weight_matrix%full_dim1), stat=istat)
     call memocc(istat, weight_matrix%matrix, 'weight_matrix%matrix', subname)
     call uncompressmatrix(iproc,weight_matrix)
     call calculate_coeffMatcoeff(weight_matrix%matrix,tmb%orbs,ksorbs,tmb%coeff,weight_coeff(1,1,ifrag))
     iall=-product(shape(weight_matrix%matrix))*kind(weight_matrix%matrix)
     deallocate(weight_matrix%matrix,stat=istat)
     call memocc(istat,iall,'weight_matrix%matrix',subname)
  end do

  !if (iproc==0) write(*,*) 'Weight analysis:'
  if (iproc==0) call yaml_open_sequence('Weight analysis',flow=.true.)
  if (iproc==0) call yaml_newline()
  !if (iproc==0) write(*,*) 'coeff, occ, eval, frac for each frag'
  if (iproc==0) call yaml_comment ('coeff, occ, eval, frac for each frag')
  ! only care about diagonal elements
  do iorb=1,ksorbs%norb
     !if (iproc==0) write(*,'(i4,2x,f6.4,1x,f10.6,2x)',ADVANCE='no') iorb,KSorbs%occup(iorb),tmb%orbs%eval(iorb)
     if (iproc==0) then
         call yaml_open_map(flow=.true.)
         call yaml_map('iorb',iorb,fmt='(i4)')
         call yaml_map('occ',KSorbs%occup(iorb),fmt='(f6.4)')
         call yaml_map('eval',tmb%orbs%eval(iorb),fmt='(f10.6)')
     end if
     do ifrag=1,input%frag%nfrag
        !if (iproc==0) write(*,'(f6.4,2x)',ADVANCE='no') weight_coeff(iorb,iorb,ifrag)
        if (iproc==0) call yaml_map('frac',weight_coeff(iorb,iorb,ifrag),fmt='(f6.4)')
     end do
     !if (iproc==0) write(*,*) ''
     if (iproc==0) call yaml_close_map()
     if (iproc==0) call yaml_newline()
  end do
  if (iproc==0) call yaml_close_sequence()

  call deallocate_sparseMatrix(weight_matrix, subname)

  iall=-product(shape(weight_coeff))*kind(weight_coeff)
  deallocate(weight_coeff,stat=istat)
  call memocc(istat,iall,'weight_coeff',subname)

end subroutine coeff_weight_analysis


! subset of reordering coeffs - need to arrange this routines better but taking the lazy route for now
! (also assuming we have no extra - or rather number of extra bands come from input.mix not input.lin)
subroutine find_eval_from_coeffs(iproc, nproc, ksorbs, basis_orbs, ham, ovrlp, coeff, eval, calc_overlap, diag)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc
  type(orbitals_data), intent(in) :: basis_orbs, ksorbs
  type(sparseMatrix),intent(in) :: ham, ovrlp
  real(kind=8),dimension(basis_orbs%norb,ksorbs%norb),intent(inout) :: coeff
  real(kind=8),dimension(ksorbs%norb),intent(inout) :: eval
  logical, intent(in) :: diag, calc_overlap

  integer :: iorb, jorb, istat, iall, ierr, itmb, jtmb
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, ham_coeff,  ovrlp_coeff
  real(kind=8) :: offdiagsum, offdiagsum2, coeff_orthog_threshold
  character(len=256) :: subname='reordering_coeffs'

  coeff_orthog_threshold=1.0d-3

  allocate(ham_coeff(ksorbs%norb,ksorbs%norb), stat=istat)
  call memocc(istat, ham_coeff, 'ham_coeff', subname)

  call calculate_coeffMatcoeff(ham%matrix,basis_orbs,ksorbs,coeff,ham_coeff)

  if (calc_overlap) then
     allocate(ovrlp_coeff(ksorbs%norb,ksorbs%norb), stat=istat)
     call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)
     call calculate_coeffMatcoeff(ovrlp%matrix,basis_orbs,ksorbs,coeff,ovrlp_coeff)
  end if

  ! above is overkill, actually just want diagonal elements but print off as a test out of curiosity
  offdiagsum=0.0d0
  offdiagsum2=0.0d0
  do iorb=1,ksorbs%norb
     do jorb=1,ksorbs%norb
        if (iorb==jorb) then
           eval(iorb)=ham_coeff(iorb,iorb)
        else
           offdiagsum=offdiagsum+abs(ham_coeff(iorb,jorb))
           if (calc_overlap) offdiagsum2=offdiagsum2+abs(ovrlp_coeff(iorb,jorb))
        end if
     end do
  end do
  offdiagsum=offdiagsum/(ksorbs%norb**2-ksorbs%norb)
  if (calc_overlap) offdiagsum2=offdiagsum2/(ksorbs%norb**2-ksorbs%norb)
  if (iproc==0) print*,''
  if (calc_overlap) then
     if (iproc==0) print*,'offdiagsum (ham,ovrlp):',offdiagsum,offdiagsum2
  else
     if (iproc==0) print*,'offdiagsum (ham):',offdiagsum
  end if

  ! if coeffs are too far from orthogonality
  if (calc_overlap .and. offdiagsum2>coeff_orthog_threshold) then
     call reorthonormalize_coeff(iproc, nproc, ksorbs%norb, -8, -8, 0, basis_orbs, ovrlp, coeff, ksorbs)
  end if

  if (diag.or.offdiagsum>1.0d-2) then
     ! diagonalize within the space of occ+extra
     if (.not.calc_overlap) then
        ! assume ovrlp_coeff is orthgonal for now
        allocate(ovrlp_coeff(ksorbs%norb,ksorbs%norb), stat=istat)
        call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)
        call calculate_coeffMatcoeff(ovrlp%matrix,basis_orbs,ksorbs,coeff,ovrlp_coeff)
        do iorb=1,ksorbs%norb
           do jorb=1,ksorbs%norb
              if (iorb==jorb) then
                 ovrlp_coeff(iorb,jorb)=1.0d0
              else
                 ovrlp_coeff(iorb,jorb)=0.0d0
              end if
           end do
        end do
     end if
     ! diagonalize within the space of occ+extra
     call diagonalizeHamiltonian2(iproc, ksorbs%norb, ham_coeff, ovrlp_coeff, eval)
     coeff_tmp=f_malloc((/basis_orbs%norb,ksorbs%norb/),id='coeff_tmp')

     ! multiply new eigenvectors by coeffs
     call dgemm('n', 'n', basis_orbs%norb, ksorbs%norb, ksorbs%norb, 1.d0, coeff(1,1), &
          basis_orbs%norb, ham_coeff, ksorbs%norb, 0.d0, coeff_tmp, basis_orbs%norb)
 
     call dcopy(basis_orbs%norb*(ksorbs%norb),coeff_tmp(1,1),1,coeff(1,1),1)
     call f_free(coeff_tmp)
  end if

  if (calc_overlap.or.diag) then
     iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
     deallocate(ovrlp_coeff,stat=istat)
     call memocc(istat,iall,'ovrlp_coeff',subname)
  end if

  iall=-product(shape(ham_coeff))*kind(ham_coeff)
  deallocate(ham_coeff,stat=istat)
  call memocc(istat,iall,'ham_coeff',subname)

end subroutine find_eval_from_coeffs


subroutine calculate_coeffMatcoeff(matrix,basis_orbs,ksorbs,coeff,mat_coeff)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data), intent(in) :: basis_orbs, ksorbs
  real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(in) :: matrix
  real(kind=8),dimension(basis_orbs%norb,ksorbs%norb),intent(inout) :: coeff
  real(kind=8),dimension(ksorbs%norb,ksorbs%norb),intent(inout) :: mat_coeff

  integer :: iall, istat, ierr
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp
  character(len=256) :: subname='calculate_coeffMatcoeff'

  allocate(coeff_tmp(basis_orbs%norbp,ksorbs%norb), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)

  if (basis_orbs%norbp>0) then
     call dgemm('n', 'n', basis_orbs%norbp, ksorbs%norb, basis_orbs%norb, 1.d0, matrix(basis_orbs%isorb+1,1), &
          basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
     call dgemm('t', 'n', ksorbs%norb, ksorbs%norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
          basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, mat_coeff, ksorbs%norb)
  else
     call to_zero(ksorbs%norb**2, mat_coeff(1,1))
  end if

  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=istat)
  call memocc(istat,iall,'coeff_tmp',subname)

  if (bigdft_mpi%nproc>1) then
      call mpiallred(mat_coeff(1,1), ksorbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

end subroutine calculate_coeffMatcoeff


 
  ! not really fragment related so prob should be moved - reorders coeffs by eval
  subroutine order_coeffs_by_energy(nstate,ntmb,coeff,eval)
    use module_base
    use module_types
    implicit none
    integer, intent(in) :: nstate, ntmb
    real(kind=gp), dimension(ntmb,nstate), intent(inout) :: coeff
    real(kind=gp), dimension(nstate), intent(inout) :: eval

    integer :: itmb, jorb
    integer, allocatable, dimension(:) :: ipiv
    real(gp), dimension(:), allocatable :: tmp_array
    real(gp), dimension(:,:), allocatable :: tmp_array2

    ipiv=f_malloc(nstate,id='coeff_final')
    tmp_array=f_malloc(nstate,id='tmp_array')

    do itmb=1,nstate
       tmp_array(itmb)=-eval(itmb)
    end do

    call sort_positions(nstate,tmp_array,ipiv)

    do itmb=1,nstate
       eval(itmb)=-tmp_array(ipiv(itmb))
    end do

    call f_free(tmp_array)

    tmp_array2=f_malloc((/ntmb,nstate/),id='tmp_array2')

    do jorb=1,nstate
       do itmb=1,ntmb
          tmp_array2(itmb,jorb)=coeff(itmb,jorb)
       end do
    end do

    do jorb=1,nstate
       do itmb=1,ntmb
          coeff(itmb,jorb)=tmp_array2(itmb,ipiv(jorb))
       end do
    end do

    call f_free(tmp_array2)
    call f_free(ipiv)

  end subroutine order_coeffs_by_energy

! experimental - for num_extra see if extra states are actually lower in energy and reorder (either by sorting or diagonalization)
! could improve efficiency by only calculating cSc or cHc up to norb (i.e. ksorbs%norb+extra)
subroutine reordering_coeffs(iproc, nproc, num_extra, ksorbs, basis_orbs, ham, ovrlp, coeff, reorder)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, num_extra
  type(orbitals_data), intent(in) :: basis_orbs, ksorbs
  type(sparseMatrix),intent(in) :: ham, ovrlp
  real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(inout) :: coeff
  logical, intent(in) :: reorder

  integer :: iorb, jorb, istat, iall, ierr, itmb, jtmb
  integer, allocatable, dimension(:) :: ipiv
  real(gp), dimension(:), allocatable :: tmp_array, eval
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, ham_coeff, tmp_array2, ham_coeff_small, ovrlp_coeff_small, ovrlp_coeff
  real(kind=8) :: offdiagsum, offdiagsum2
  character(len=256) :: subname='reordering_coeffs'

  allocate(ham_coeff(basis_orbs%norb,basis_orbs%norb), stat=istat)
  call memocc(istat, ham_coeff, 'ham_coeff', subname)

  allocate(coeff_tmp(basis_orbs%norbp,max(basis_orbs%norb,1)), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)

  if (basis_orbs%norbp>0) then
     call dgemm('n', 'n', basis_orbs%norbp, basis_orbs%norb, basis_orbs%norb, 1.d0, ham%matrix(basis_orbs%isorb+1,1), &
          basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
     call dgemm('t', 'n', basis_orbs%norb, basis_orbs%norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
          basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, ham_coeff, basis_orbs%norb)
  else
     call to_zero(basis_orbs%norb**2, ham_coeff(1,1))
  end if

  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=istat)
  call memocc(istat,iall,'coeff_tmp',subname)

  if (nproc>1) then
      call mpiallred(ham_coeff(1,1), basis_orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  allocate(ovrlp_coeff(basis_orbs%norb,basis_orbs%norb), stat=istat)
  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)

  allocate(coeff_tmp(basis_orbs%norbp,max(basis_orbs%norb,1)), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)

  if (basis_orbs%norbp>0) then
     call dgemm('n', 'n', basis_orbs%norbp, basis_orbs%norb, basis_orbs%norb, 1.d0, ovrlp%matrix(basis_orbs%isorb+1,1), &
          basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
     call dgemm('t', 'n', basis_orbs%norb, basis_orbs%norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
          basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, ovrlp_coeff, basis_orbs%norb)
  else
     call to_zero(basis_orbs%norb**2, ovrlp_coeff(1,1))
  end if

  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=istat)
  call memocc(istat,iall,'coeff_tmp',subname)

  if (nproc>1) then
      call mpiallred(ovrlp_coeff(1,1), basis_orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  ! above is overkill, actually just want diagonal elements but print off as a test out of curiosity
  offdiagsum=0.0d0
  offdiagsum2=0.0d0
  do iorb=1,ksorbs%norb+num_extra!basis_orbs%norb
     do jorb=1,ksorbs%norb+num_extra!basis_orbs%norb
        if (iorb==jorb) cycle
        offdiagsum=offdiagsum+abs(ham_coeff(iorb,jorb))
        offdiagsum2=offdiagsum2+abs(ovrlp_coeff(iorb,jorb))
     end do
  end do
  offdiagsum=offdiagsum/(basis_orbs%norb**2-basis_orbs%norb)
  offdiagsum2=offdiagsum2/(basis_orbs%norb**2-basis_orbs%norb)
  if (iproc==0) print*,'offdiagsum (ham,ovrlp):',offdiagsum,offdiagsum2

  ! sort the states - really need just ks+extra not all, otherwise sloshing!
  ipiv=f_malloc(ksorbs%norb+num_extra,id='coeff_final')
  tmp_array=f_malloc(ksorbs%norb+num_extra,id='tmp_array')
  do itmb=1,ksorbs%norb+num_extra
     tmp_array(itmb)=-ham_coeff(itmb,itmb)
  end do

  do itmb=1,ksorbs%norb+num_extra
     ipiv(itmb)=itmb
  end do
  !call sort_positions(ksorbs%norb+num_extra,tmp_array,ipiv)
  do jtmb=1,ksorbs%norb+num_extra
     ham_coeff(jtmb,jtmb)=-tmp_array(ipiv(jtmb))
  end do

  tmp_array2=f_malloc((/basis_orbs%norb,ksorbs%norb+num_extra/),id='tmp_array2')
  do jtmb=1,ksorbs%norb+num_extra
     do itmb=1,basis_orbs%norb
        tmp_array2(itmb,jtmb)=coeff(itmb,jtmb)
     end do
  end do

  do jtmb=1,ksorbs%norb+num_extra
     do itmb=1,basis_orbs%norb
        coeff(itmb,jtmb)=tmp_array2(itmb,ipiv(jtmb))
     end do
  end do
  call f_free(tmp_array2)

  if (reorder) then
     ! diagonalize within the space of occ+extra
     eval=f_malloc((/ksorbs%norb+num_extra/),id='eval')
     ham_coeff_small=f_malloc((/ksorbs%norb+num_extra,ksorbs%norb+num_extra/),id='ham_coeff_small')
     ovrlp_coeff_small=f_malloc((/ksorbs%norb+num_extra,ksorbs%norb+num_extra/),id='ovrlp_coeff_small')
     ! assume ovrlp_coeff is orthgonal for now
     do iorb=1,ksorbs%norb+num_extra
        do jorb=1,ksorbs%norb+num_extra
           ham_coeff_small(iorb,jorb)=ham_coeff(iorb,jorb)
           if (iorb==jorb) then
              ovrlp_coeff_small(iorb,jorb)=1.0d0
           else
              ovrlp_coeff_small(iorb,jorb)=0.0d0
           end if
        end do
     end do
     ! diagonalize within the space of occ+extra
     call diagonalizeHamiltonian2(iproc, ksorbs%norb+num_extra, ham_coeff_small, ovrlp_coeff_small, eval)
     call f_free(ovrlp_coeff_small)
     coeff_tmp=f_malloc((/basis_orbs%norb,ksorbs%norb+num_extra/),id='coeff_tmp')

     ! multiply new eigenvectors by coeffs
     call dgemm('n', 'n', basis_orbs%norb, ksorbs%norb+num_extra, ksorbs%norb+num_extra, 1.d0, coeff(1,1), &
          basis_orbs%norb, ham_coeff_small, ksorbs%norb+num_extra, 0.d0, coeff_tmp, basis_orbs%norb)
 
     call f_free(ham_coeff_small)
     call dcopy(basis_orbs%norb*(ksorbs%norb+num_extra),coeff_tmp(1,1),1,coeff(1,1),1)
     call f_free(coeff_tmp)

     do iorb=1,basis_orbs%norb
        if (iorb<=ksorbs%norb) then
           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
                ksorbs%occup(iorb),ham_coeff(iorb,iorb),eval(iorb)
        else if (iorb<=ksorbs%norb+num_extra) then
           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
                0.0d0,ham_coeff(iorb,iorb),eval(iorb)
        !else
        !   if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,ham_coeff(iorb,iorb),basis_orbs%occup(iorb),&
        !        0.0d0,ham_coeff(iorb,iorb)
        end if
     end do
     call f_free(eval)
   else
      do iorb=1,basis_orbs%norb
        if (iorb<=ksorbs%norb) then
           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
                ksorbs%occup(iorb),ham_coeff(iorb,iorb)
        else if (iorb<=ksorbs%norb+num_extra) then
           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
                0.0d0,ham_coeff(iorb,iorb)
        !else
        !   if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,ham_coeff(iorb,iorb),basis_orbs%occup(iorb),&
        !        0.0d0,ham_coeff(iorb,iorb)
        end if
     end do
  end if

  call f_free(ipiv)
  call f_free(tmp_array)
  iall=-product(shape(ham_coeff))*kind(ham_coeff)
  deallocate(ham_coeff,stat=istat)
  call memocc(istat,iall,'ham_coeff',subname)

end subroutine reordering_coeffs



subroutine find_alpha_sd(iproc,nproc,alpha,tmb,orbs,coeffp,grad,energy0,fnrm,pred_e)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc, nproc
  real(kind=gp), intent(inout) :: alpha
  type(DFT_wavefunction) :: tmb
  type(orbitals_data), intent(in) :: orbs
  real(kind=gp), dimension(tmb%orbs%norb,orbs%norbp), intent(inout) :: coeffp
  real(kind=gp), dimension(tmb%orbs%norb,orbs%norbp), intent(in) :: grad
  real(kind=gp), intent(in) :: energy0, fnrm
  real(kind=gp), intent(out) :: pred_e
  integer :: iorb, iiorb, jorb, ierr
  real(kind=gp) :: energy1, a, b, c, alpha_old
  real(kind=gp),dimension(:,:),allocatable :: coeff_tmp

  call timing(iproc,'dirmin_sdfit','ON')

  ! take an initial step to get 2nd point
  coeff_tmp=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_tmp')
  do iorb=1,orbs%norbp
     iiorb = orbs%isorb + iorb
     do jorb=1,tmb%orbs%norb
        coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-alpha*grad(jorb,iorb)
     end do
  end do

  if(nproc > 1) then 
     call mpi_allgatherv(coeffp, tmb%orbs%norb*orbs%norbp, mpi_double_precision, coeff_tmp, &
        tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call dcopy(tmb%orbs%norb*orbs%norb,coeffp(1,1),1,coeff_tmp(1,1),1)
  end if

  ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
  ! instead of twice could add some criterion to check accuracy?
  call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, coeff_tmp, orbs)
  call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, coeff_tmp, orbs)
  call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,tmb%linmat%ham,energy1,&
       coeff_tmp,orbs,tmb%orbs,.true.)
  call f_free(coeff_tmp)

  ! find ideal alpha using both points
  alpha_old=alpha
  a=fnrm/alpha_old+(energy1-energy0)/alpha_old**2
  b=-fnrm
  c=energy0
  alpha=-0.5_gp*b/a
  ! don't update if we have found a maximum, or negative alpha is predicted
  ! do something better here - don't just want to do the same thing twice, so at least check if energy has decreased
  if (alpha<0.0_gp .or. a<0.0_gp) alpha=alpha_old
  pred_e=a*alpha**2+b*alpha+c

  !open(127)
  !write(127,*) '#',a,b,c,(energy1-energy0)/alpha_old,b-(energy1-energy0)/alpha_old
  !write(127,*) 0.0_gp,energy0
  !write(127,*) alpha_old,energy1

  call timing(iproc,'dirmin_sdfit','OF')

end subroutine find_alpha_sd


subroutine calculate_kernel_and_energy(iproc,nproc,denskern,ham,energy,coeff,orbs,tmb_orbs,calculate_kernel)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc, nproc
  type(sparseMatrix), intent(in) :: ham
  type(sparseMatrix), intent(inout) :: denskern
  logical, intent(in) :: calculate_kernel
  real(kind=gp), intent(out) :: energy
  type(orbitals_data), intent(in) :: orbs, tmb_orbs
  real(kind=gp), dimension(tmb_orbs%norb,tmb_orbs%norb), intent(in) :: coeff

  integer :: iorb, jorb, ind_ham, ind_denskern, ierr, iorbp
  integer :: matrixindex_in_compressed

  if (calculate_kernel) then 
     call calculate_density_kernel(iproc, nproc, .true., orbs, tmb_orbs, coeff, denskern)
  end if

  energy=0.0_gp
  do iorbp=1,tmb_orbs%norbp
     iorb=iorbp+tmb_orbs%isorb
     do jorb=1,tmb_orbs%norb
        ind_ham = matrixindex_in_compressed(ham,iorb,jorb)
        ind_denskern = matrixindex_in_compressed(denskern,jorb,iorb)
        if (ind_ham==0.or.ind_denskern==0) cycle
        energy = energy + denskern%matrix_compr(ind_denskern)*ham%matrix_compr(ind_ham)
     end do
  end do
  if (nproc>1) then
     call mpiallred(energy, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if


end subroutine calculate_kernel_and_energy


! calculate grad_cov_i^a = f_i (I_ab - S_ag K^gb) H_bg c_i^d
! then grad=S^-1grad_cov
subroutine calculate_coeff_gradient(iproc,nproc,tmb,KSorbs,grad_cov,grad)
  use module_base
  use module_types
  implicit none

  integer, intent(in) :: iproc, nproc
  type(DFT_wavefunction), intent(inout) :: tmb
  type(orbitals_data), intent(in) :: KSorbs
  real(gp), dimension(tmb%orbs%norb,KSorbs%norbp), intent(out) :: grad_cov, grad  ! could make grad_cov KSorbs%norbp

  integer :: iorb, iiorb, info, ierr
  real(gp),dimension(:,:),allocatable :: sk, skh, skhp, inv_ovrlp
  integer,dimension(:),allocatable:: ipiv
  real(kind=gp), dimension(:,:), allocatable:: grad_full
  character(len=*),parameter:: subname='calculate_coeff_gradient'

  call f_routine(id='calculate_coeff_gradient')
  call timing(iproc,'dirmin_lagmat1','ON')

  ! we have the kernel already, but need it to not contain occupations so recalculate here
  ! don't want to lose information in the compress/uncompress process - ideally need to change sparsity pattern of kernel
  !call calculate_density_kernel(iproc, nproc, .false., KSorbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern)
  tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='denskern')
  !call uncompressMatrix(iproc,tmb%linmat%denskern)
  call calculate_density_kernel_uncompressed(iproc, nproc, .false., KSorbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern%matrix)

  sk=f_malloc0((/tmb%orbs%norbp,tmb%orbs%norb/), id='sk')

  ! calculate I-S*K - first set sk to identity
  do iorb=1,tmb%orbs%norbp
     iiorb=tmb%orbs%isorb+iorb
     sk(iorb,iiorb) = 1.d0
  end do 

  if (tmb%orbs%norbp>0) then
     call dgemm('t', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, -1.d0, &
          tmb%linmat%ovrlp%matrix(1,tmb%orbs%isorb+1), tmb%orbs%norb, &
          tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 1.d0, sk, tmb%orbs%norbp)
  end if

  ! coeffs and therefore kernel will change, so no need to keep it
  call f_free_ptr(tmb%linmat%denskern%matrix)

  skhp=f_malloc((/tmb%orbs%norb,max(tmb%orbs%norbp,1)/), id='skhp')

  ! multiply by H to get (I_ab - S_ag K^gb) H_bd, or in this case the transpose of the above
  if (tmb%orbs%norbp>0) then
     call dgemm('t', 't', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, 1.d0, tmb%linmat%ham%matrix(1,1), &
          tmb%orbs%norb, sk(1,1), tmb%orbs%norbp, 0.d0, skhp(1,1), tmb%orbs%norb)
  end if

  call f_free(sk)

  skh=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/), id='skh')

  call timing(iproc,'dirmin_lagmat1','OF')
  call timing(iproc,'dirmin_lagmat2','ON')

  ! gather together
  if(nproc > 1) then
     call mpi_allgatherv(skhp(1,1), tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, skh(1,1), &
        tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
        mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call dcopy(tmb%orbs%norbp*tmb%orbs%norb,skhp(1,1),1,skh(1,1),1)
  end if

  call timing(iproc,'dirmin_lagmat2','OF')
  call timing(iproc,'dirmin_lagmat1','ON')

  call f_free(skhp)

  ! calc for i on this proc: (I_ab - S_ag K^gb) H_bg c_i^d
  if (KSorbs%norbp>0) then
     call dgemm('t', 'n', tmb%orbs%norb, KSorbs%norbp, tmb%orbs%norb, 1.d0, skh(1,1), &
          tmb%orbs%norb, tmb%coeff(1,KSorbs%isorb+1), tmb%orbs%norb, 0.d0, grad_cov(1,1), tmb%orbs%norb)
  end if

  call f_free(skh)

  ! multiply by f_i to get grad_i^a
  do iorb=1,KSorbs%norbp
     iiorb=KSorbs%isorb+iorb
     grad_cov(:,iorb)=grad_cov(:,iorb)*KSorbs%occup(iiorb)
  end do

  call timing(iproc,'dirmin_lagmat1','OF')
  call timing(iproc,'dirmin_dgesv','ON') !lr408t

  info = 0 ! needed for when some processors have orbs%norbp=0
  ! Solve the linear system ovrlp*grad=grad_cov
  if(tmb%orthpar%blocksize_pdsyev<0) then
     !! keep the covariant gradient to calculate fnrm correctly
     !call dcopy(tmb%orbs%norb*KSorbs%norbp,grad_cov,1,grad,1)
     !if (KSorbs%norbp>0) then
     !   ipiv=f_malloc(tmb%orbs%norb,id='ipiv')
     !   call dgesv(tmb%orbs%norb, KSorbs%norbp, tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, ipiv(1), &
     !        grad(1,1), tmb%orbs%norb, info)
     !   call f_free(ipiv)
     !end if
     inv_ovrlp=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
     call overlapPowerMinusOne(iproc, nproc, 1, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, inv_ovrlp)
     if (KSorbs%norbp>0) then
        call dgemm('n', 'n', tmb%orbs%norb, KSorbs%norbp, tmb%orbs%norb, 1.d0, inv_ovrlp(1,1), &
             tmb%orbs%norb, grad_cov(1,1), tmb%orbs%norb, 0.d0, grad(1,1), tmb%orbs%norb)
     else
        call dcopy(tmb%orbs%norb*KSorbs%norbp,grad_cov,1,grad,1)
     end if
     call f_free(inv_ovrlp)
  else
      grad_full=f_malloc((/tmb%orbs%norb,KSorbs%norb/),id='grad_full')
      ! do allgather instead of allred so we can keep grad as per proc
      if(nproc > 1) then 
         call mpi_allgatherv(grad_cov, tmb%orbs%norb*KSorbs%norbp, mpi_double_precision, grad_full, &
            tmb%orbs%norb*KSorbs%norb_par(:,0), tmb%orbs%norb*KSorbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
         call dcopy(tmb%orbs%norb*KSorbs%norb,grad_cov(1,1),1,grad_full(1,1),1)
      end if
      !call mpiallred(grad(1,1), tmb%orbs%norb*KSorbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)

      call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
           tmb%orbs%norb, KSorbs%norb, tmb%linmat%ovrlp%matrix, tmb%orbs%norb, grad_full, tmb%orbs%norb, info)

      call dcopy(tmb%orbs%norb*KSorbs%norbp,grad_full(1,KSorbs%isorb+1),1,grad(1,1),1)

      call f_free(grad_full)
  end if

  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
      stop
  end if


  call timing(iproc,'dirmin_dgesv','OF') !lr408t
  call f_release_routine()

end subroutine calculate_coeff_gradient

! calculate grad_cov_i^a = f_i (I_ab - S_ag K^gb) H_bg c_i^d
! then grad=S^-1grad_cov
subroutine calculate_coeff_gradient_extra(iproc,nproc,num_extra,tmb,KSorbs,grad_cov,grad)
  use module_base
  use module_types
  implicit none

  integer, intent(in) :: iproc, nproc, num_extra
  type(DFT_wavefunction), intent(inout) :: tmb
  type(orbitals_data), intent(in) :: KSorbs
  real(gp), dimension(tmb%orbs%norb,tmb%orbs%norbp), intent(out) :: grad_cov, grad  ! could make grad_cov KSorbs%norbp

  integer :: iorb, iiorb, info, ierr
  real(gp),dimension(:,:),allocatable :: sk, skh, skhp, inv_ovrlp
  integer :: matrixindex_in_compressed
  integer,dimension(:),allocatable:: ipiv
  real(kind=gp), dimension(:), allocatable:: occup_tmp
  real(kind=gp), dimension(:,:), allocatable:: grad_full
  character(len=*),parameter:: subname='calculate_coeff_gradient'

  call f_routine(id='calculate_coeff_gradient')
  call timing(iproc,'dirmin_lagmat1','ON')

  occup_tmp=f_malloc(tmb%orbs%norb,id='occup_tmp')
  call dcopy(tmb%orbs%norb,tmb%orbs%occup(1),1,occup_tmp(1),1)

  call razero(tmb%orbs%norb,tmb%orbs%occup(1))
  do iorb=1,KSorbs%norb+num_extra
     tmb%orbs%occup(iorb)=1.0d0
  end do

  ! we have the kernel already, but need it to not contain occupations so recalculate here
  !call calculate_density_kernel(iproc, nproc, .true., tmb%orbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern)
  tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='denskern')
  !call uncompressMatrix(iproc,tmb%linmat%denskern)
  call calculate_density_kernel_uncompressed (iproc, nproc, .true., tmb%orbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern%matrix)

  call dcopy(tmb%orbs%norb,occup_tmp(1),1,tmb%orbs%occup(1),1)
  call f_free(occup_tmp)

  sk=f_malloc0((/tmb%orbs%norbp,tmb%orbs%norb/), id='sk')

  ! calculate I-S*K - first set sk to identity
  do iorb=1,tmb%orbs%norbp
     iiorb=tmb%orbs%isorb+iorb
     sk(iorb,iiorb) = 1.d0
  end do 

  if (tmb%orbs%norbp>0) then
     call dgemm('t', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, -1.d0, &
          tmb%linmat%ovrlp%matrix(1,tmb%orbs%isorb+1), tmb%orbs%norb, &
          tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 1.d0, sk, tmb%orbs%norbp)
  end if

  ! coeffs and therefore kernel will change, so no need to keep it
  call f_free_ptr(tmb%linmat%denskern%matrix)

  skhp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='skhp')

  ! multiply by H to get (I_ab - S_ag K^gb) H_bd, or in this case the transpose of the above
  if (tmb%orbs%norbp>0) then
     call dgemm('t', 't', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, 1.d0, tmb%linmat%ham%matrix(1,1), &
          tmb%orbs%norb, sk(1,1), tmb%orbs%norbp, 0.d0, skhp(1,1), tmb%orbs%norb)
  end if

  call f_free(sk)

  skh=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/), id='skh')

  call timing(iproc,'dirmin_lagmat1','OF')
  call timing(iproc,'dirmin_lagmat2','ON')

  ! gather together
  if(nproc > 1) then
     call mpi_allgatherv(skhp(1,1), tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, skh(1,1), &
        tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
        mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call dcopy(tmb%orbs%norbp*tmb%orbs%norb,skhp(1,1),1,skh(1,1),1)
  end if

  call timing(iproc,'dirmin_lagmat2','OF')
  call timing(iproc,'dirmin_lagmat1','ON')

  call f_free(skhp)

  ! calc for i on this proc: (I_ab - S_ag K^gb) H_bg c_i^d
  if (tmb%orbs%norbp>0) then
     call dgemm('t', 'n', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, 1.d0, skh(1,1), &
          tmb%orbs%norb, tmb%coeff(1,tmb%orbs%isorb+1), tmb%orbs%norb, 0.d0, grad_cov(1,1), tmb%orbs%norb)
  end if

  call f_free(skh)

  ! multiply by f_i to get grad_i^a
  do iorb=1,tmb%orbs%norbp
     iiorb=tmb%orbs%isorb+iorb
     grad_cov(:,iorb)=grad_cov(:,iorb)*tmb%orbs%occup(iiorb)
  end do

  call timing(iproc,'dirmin_lagmat1','OF')
  call timing(iproc,'dirmin_dgesv','ON') !lr408t

  info = 0 ! needed for when some processors have orbs%norbp=0
  ! Solve the linear system ovrlp*grad=grad_cov
  if(tmb%orthpar%blocksize_pdsyev<0) then
     !! keep the covariant gradient to calculate fnrm correctly
     !call dcopy(tmb%orbs%norb*tmb%orbs%norbp,grad_cov,1,grad,1)
     !if (tmb%orbs%norbp>0) then
     !   ipiv=f_malloc(tmb%orbs%norb,id='ipiv')
     !   call dgesv(tmb%orbs%norb, tmb%orbs%norbp, tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, ipiv(1), &
     !        grad(1,1), tmb%orbs%norb, info)
     !   call f_free(ipiv)
     !end if
     inv_ovrlp=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
     call overlapPowerMinusOne(iproc, nproc, 1, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, inv_ovrlp)
     if (tmb%orbs%norbp>0) then
        call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, 1.d0, inv_ovrlp(1,1), &
             tmb%orbs%norb, grad_cov(1,1), tmb%orbs%norb, 0.d0, grad(1,1), tmb%orbs%norb)
     else
        call dcopy(tmb%orbs%norb*tmb%orbs%norbp,grad_cov,1,grad,1)
     end if
     call f_free(inv_ovrlp)
  else
      grad_full=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='grad_full')
      ! do allgather instead of allred so we can keep grad as per proc
      if(nproc > 1) then 
         call mpi_allgatherv(grad_cov, tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, grad_full, &
            tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
         call dcopy(tmb%orbs%norb*tmb%orbs%norb,grad_cov(1,1),1,grad_full(1,1),1)
      end if
      !call mpiallred(grad(1,1), tmb%orbs%norb*tmb%orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)

      call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
           tmb%orbs%norb, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, tmb%orbs%norb, grad_full, tmb%orbs%norb, info)

      call dcopy(tmb%orbs%norb*tmb%orbs%norbp,grad_full(1,tmb%orbs%isorb+1),1,grad(1,1),1)

      call f_free(grad_full)
  end if

  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
      stop
  end if


  call timing(iproc,'dirmin_dgesv','OF') !lr408t
  call f_release_routine()

end subroutine calculate_coeff_gradient_extra

subroutine precondition_gradient_coeff(ntmb, norb, ham, ovrlp, grad)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: ntmb, norb
  real(8),dimension(ntmb,ntmb),intent(in):: ham, ovrlp
  real(8),dimension(ntmb,norb),intent(inout):: grad
  
  ! Local variables
  integer:: iorb, itmb, jtmb, info, istat, iall
  complex(8),dimension(:,:),allocatable:: mat
  complex(8),dimension(:,:),allocatable:: rhs
  integer,dimension(:),allocatable:: ipiv
  character(len=*),parameter:: subname='precondition_gradient_coeff'
  
  allocate(mat(ntmb,ntmb), stat=istat)
  !call memocc(istat, mat, 'mat', subname)
  allocate(rhs(ntmb,norb), stat=istat)
  !call memocc(istat, mat, 'mat', subname)
  
  ! Build the matrix to be inverted
  do itmb=1,ntmb
      do jtmb=1,ntmb
          mat(jtmb,itmb) = cmplx(ham(jtmb,itmb)+.5d0*ovrlp(jtmb,itmb),0.d0,kind=8)
      end do
      mat(itmb,itmb)=mat(itmb,itmb)+cmplx(0.d0,-1.d-1,kind=8)
      !mat(itmb,itmb)=mat(itmb,itmb)-cprec
  end do
  do iorb=1,norb
      do itmb=1,ntmb
          rhs(itmb,iorb)=cmplx(grad(itmb,iorb),0.d0,kind=8)
      end do
  end do
  
  allocate(ipiv(ntmb), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)
  
  call zgesv(ntmb, norb, mat(1,1), ntmb, ipiv, rhs(1,1), ntmb, info)
  if(info/=0) then
      stop 'ERROR in dgesv'
  end if
  !call dcopy(nel, rhs(1), 1, grad(1), 1)
  do iorb=1,norb
      do itmb=1,ntmb
          grad(itmb,iorb)=real(rhs(itmb,iorb))
      end do
  end do
  
  iall=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv, stat=istat)
  call memocc(istat, iall, 'ipiv', subname)
  
  iall=-product(shape(mat))*kind(mat)
  deallocate(mat, stat=istat)
  !call memocc(istat, iall, 'mat', subname)
  
  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  !call memocc(istat, iall, 'rhs', subname)

end subroutine precondition_gradient_coeff



subroutine DIIS_coeff(iproc, orbs, tmb, grad, coeff, ldiis)
  use module_base
  use module_types
  use module_interfaces, except_this_one => DIIS_coeff
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(in):: tmb
  real(8),dimension(tmb%orbs%norb*tmb%orbs%norbp),intent(in):: grad
  real(8),dimension(tmb%orbs%norb*tmb%orbs%norb),intent(inout):: coeff
  type(localizedDIISParameters),intent(inout):: ldiis
  
  ! Local variables
  integer:: iorb, jorb, ist, ncount, jst, i, j, mi, ist1, ist2, istat, lwork, info
  integer:: mj, jj, k, jjst, isthist, iall
  real(8):: ddot
  real(8),dimension(:,:),allocatable:: mat
  real(8),dimension(:),allocatable:: rhs, work
  integer,dimension(:),allocatable:: ipiv
  character(len=*),parameter:: subname='DIIS_coeff'
  
  !!call timing(iproc,'optimize_DIIS ','ON')
  
  ! Allocate the local arrays.
  allocate(mat(ldiis%isx+1,ldiis%isx+1), stat=istat)
  call memocc(istat, mat, 'mat', subname)
  allocate(rhs(ldiis%isx+1), stat=istat)
  call memocc(istat, rhs, 'rhs', subname)
  allocate(ipiv(ldiis%isx+1), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)
  
  mat=0.d0
  rhs=0.d0
  call to_zero((ldiis%isx+1)**2, mat(1,1))
  call to_zero(ldiis%isx+1, rhs(1))
  
  ncount=tmb%orbs%norb

  ! Copy coeff and grad to history.
  ist=1
  do iorb=1,tmb%orbs%norbp
      jst=1
      do jorb=1,iorb-1
          jst=jst+ncount*ldiis%isx
      end do
      jst=jst+(ldiis%mis-1)*ncount
      call dcopy(ncount, coeff(ist+tmb%orbs%isorb*tmb%orbs%norb), 1, ldiis%phiHist(jst), 1)
      call dcopy(ncount, grad(ist), 1, ldiis%hphiHist(jst), 1)
      ist=ist+ncount
  end do
  
  do iorb=1,tmb%orbs%norbp
      ! Shift the DIIS matrix left up if we reached the maximal history length.
      if(ldiis%is>ldiis%isx) then
         do i=1,ldiis%isx-1
            do j=1,i
               ldiis%mat(j,i,iorb)=ldiis%mat(j+1,i+1,iorb)
            end do
         end do
      end if
  end do
  
  do iorb=1,tmb%orbs%norbp
      ! Calculate a new line for the matrix.
      i=max(1,ldiis%is-ldiis%isx+1)
      jst=1
      ist1=1
      do jorb=1,iorb-1
          jst=jst+ncount*ldiis%isx
          ist1=ist1+ncount
      end do
      do j=i,ldiis%is
         mi=mod(j-1,ldiis%isx)+1
         ist2=jst+(mi-1)*ncount
         if(ist2>size(ldiis%hphiHist)) then
             write(*,'(a,7i8)') 'ERROR ist2: iproc, iorb, ldiis%is, mi, ncount, ist2, size(ldiis%hphiHist)', iproc, iorb, ldiis%is,&
                                 mi, ncount, ist2, size(ldiis%hphiHist)
         end if
         ldiis%mat(j-i+1,min(ldiis%isx,ldiis%is),iorb)=ddot(ncount, grad(ist1), 1, ldiis%hphiHist(ist2), 1)
         ist2=ist2+ncount
      end do
  end do
  
  ist=1+tmb%orbs%isorb*tmb%orbs%norb
  do iorb=1,tmb%orbs%norbp
      ! Copy the matrix to an auxiliary array and fill with the zeros and ones.
      do i=1,min(ldiis%isx,ldiis%is)
          mat(i,min(ldiis%isx,ldiis%is)+1)=1.d0
          rhs(i)=0.d0
          do j=i,min(ldiis%isx,ldiis%is)
              mat(i,j)=ldiis%mat(i,j,iorb)
          end do
      end do
      mat(min(ldiis%isx,ldiis%is)+1,min(ldiis%isx,ldiis%is)+1)=0.d0
      rhs(min(ldiis%isx,ldiis%is)+1)=1.d0
   
      ! Solve the linear system
      !!do istat=1,ldiis%isx+1
          !!do iall=1,ldiis%isx+1
              !!if(iproc==0) write(500,*) istat, iall, mat(iall,istat)
          !!end do
      !!end do

      if(ldiis%is>1) then
         lwork=-1   !100*ldiis%isx
         allocate(work(1000), stat=istat)
         call memocc(istat, work, 'work', subname)
         call dsysv('u', min(ldiis%isx,ldiis%is)+1, 1, mat, ldiis%isx+1,  & 
              ipiv, rhs(1), ldiis%isx+1, work, lwork, info)
         lwork=nint(work(1))
         iall=-product(shape(work))*kind(work)
         deallocate(work,stat=istat)
         call memocc(istat,iall,'work',subname)
         allocate(work(lwork), stat=istat)
         call memocc(istat, work, 'work', subname)
         call dsysv('u', min(ldiis%isx,ldiis%is)+1, 1, mat, ldiis%isx+1,  & 
              ipiv, rhs(1), ldiis%isx+1, work, lwork, info)
         iall=-product(shape(work))*kind(work)
         deallocate(work, stat=istat)
         call memocc(istat, iall, 'work', subname)
         
         if (info /= 0) then
            write(*,'(a,i0)') 'ERROR in dsysv (DIIS_coeff), info=', info
            stop
         end if
      else
         rhs(1)=1.d0
      endif
    
      ! Make a new guess for the orbital.
      call razero(ncount, coeff(ist))
      isthist=max(1,ldiis%is-ldiis%isx+1)
      jj=0
      jst=0
      do jorb=1,iorb-1
          jst=jst+ncount*ldiis%isx
      end do
      do j=isthist,ldiis%is
          jj=jj+1
          mj=mod(j-1,ldiis%isx)+1
          jjst=jst+(mj-1)*ncount
          do k=1,ncount
              coeff(ist+k-1) = coeff(ist+k-1) + rhs(jj)*(ldiis%phiHist(jjst+k)-ldiis%hphiHist(jjst+k))
          end do
      end do
      ist=ist+ncount
  end do
    
  iall=-product(shape(mat))*kind(mat)
  deallocate(mat, stat=istat)
  call memocc(istat, iall, 'mat', subname)
  
  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  call memocc(istat, iall, 'rhs', subname)

  iall=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv, stat=istat)
  call memocc(istat, iall, 'ipiv', subname)
  
  !!call timing(iproc,'optimize_DIIS ','OF')

end subroutine DIIS_coeff


subroutine initialize_DIIS_coeff(isx, ldiis)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: isx
  type(localizedDIISParameters),intent(inout):: ldiis
  
  ! Local variables
  character(len=*),parameter:: subname='initialize_DIIS_coeff'
    
  ldiis%isx=isx
  ldiis%is=0
  ldiis%switchSD=.false.
  ldiis%trmin=1.d100
  ldiis%trold=1.d100
  ldiis%alpha_coeff=0.1d0

end subroutine initialize_DIIS_coeff


subroutine allocate_DIIS_coeff(tmb, ldiis)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(DFT_wavefunction),intent(in):: tmb
  type(localizedDIISParameters),intent(inout):: ldiis
  
  ! Local variables
  integer:: ii, istat
  character(len=*),parameter:: subname='allocate_DIIS_coeff'

  allocate(ldiis%mat(ldiis%isx,ldiis%isx,tmb%orbs%norbp),stat=istat)
  call memocc(istat, ldiis%mat, 'ldiis%mat', subname)

  ii=ldiis%isx*tmb%orbs%norb*tmb%orbs%norbp
  allocate(ldiis%phiHist(ii), stat=istat)
  call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
  allocate(ldiis%hphiHist(ii), stat=istat)
  call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)

end subroutine allocate_DIIS_coeff

