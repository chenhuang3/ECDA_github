!> @file
!! @brief External routines of the BigDFT library.
!! @details
!! To be documented in detail once stabilized
!! All the call to BigDFT code should be performed from these routines
!! No interface should be required to manipulate these routines
!! Non-intrinsic objects should be mapped to addresses which have to be manipulated
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Routine which initalizes the BigDFT environment
subroutine bigdft_init(mpi_info,nconfig,run_id,ierr)
  use BigDFT_API
  implicit none
  integer, dimension(4), intent(out) :: mpi_info !< first entry: id of MPI task in the groups,
                                                 !! 2nd: number of MPI tasks, third: id of group, fourth: total number of taskgroups
  integer, intent(out) :: nconfig                !< if negative, run_is is a list_posinp, otherwise comes from the taskgroups
  character(len=*), intent(out) :: run_id        !< radical of the taskgroups or list_posinp name
  integer, intent(out) :: ierr                   !< error code
  !local variables
  logical :: exist_list
  integer :: nconfig_file,mpi_groupsize
  character(len=60) :: posinp_file,radical

  !Initalize the global mpi environment
  call bigdft_mpi_init(ierr)
  if (ierr /= MPI_SUCCESS) return

  call bigdft_init_errors()

  call command_line_information(mpi_groupsize,posinp_file,radical,ierr)

  call bigdft_init_mpi_env(mpi_info, mpi_groupsize, ierr)

  !minimum number of different configurations dictated by ngroups
  nconfig=bigdft_mpi%ngroup
  write(run_id,'(a)')trim(radical)

  if (len_trim(posinp_file) > 0) then
     inquire(file=trim(posinp_file),exist=exist_list)
     if (exist_list) then
        write(run_id,'(a)')posinp_file
        open(54,file=trim(run_id))
        read(54,*) nconfig_file
        close(54)
        nconfig=-nconfig_file
     else
        write(*,'(a)')'ERROR (bigdft_init): runs-file absent'
        stop
     end if
  end if
end subroutine bigdft_init

subroutine bigdft_init_mpi_env(mpi_info,mpi_groupsize, ierr)
  use BigDFT_API
  implicit none
  
  integer, dimension(4), intent(out) :: mpi_info
  integer, intent(in) :: mpi_groupsize
  integer, intent(out) :: ierr
  !local variables
  integer :: iproc,nproc,ngroup_size

  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  if (ierr /= MPI_SUCCESS) return

  !set the memory limit for the allocation library
  call f_malloc_set_status(memory_limit=memorylimit,iproc=iproc)
  !call memocc_set_memory_limit(memorylimit)

!!$  print *,'list_posinp',trim(posinp_file),'iproc',iproc
!!$  print *,'run_id',trim(radical),'iproc',iproc
!!$  print *,'mpi_groupsize',mpi_groupsize,'iproc',iproc

  !if the taskgroup size is not a divisor of nproc do not create taskgroups
  if (nproc >1 .and. mpi_groupsize > 0 .and. mpi_groupsize < nproc .and.&
       mod(nproc,mpi_groupsize)==0) then
     ngroup_size=mpi_groupsize
  else
     ngroup_size=nproc
  end if
  call mpi_environment_set(bigdft_mpi,iproc,nproc,MPI_COMM_WORLD,ngroup_size)

  !final values
  mpi_info(1)=bigdft_mpi%iproc
  mpi_info(2)=bigdft_mpi%nproc
  mpi_info(3)=bigdft_mpi%igroup
  mpi_info(4)=bigdft_mpi%ngroup
end subroutine bigdft_init_mpi_env

subroutine bigdft_init_mpi_force(igroup, ngroup)
  use BigDFT_API
  implicit none
  integer, intent(in) :: igroup, ngroup

  if (igroup >= 0) bigdft_mpi%igroup = igroup
  if (ngroup >= 0) bigdft_mpi%ngroup = ngroup
END SUBROUTINE bigdft_init_mpi_force

subroutine bigdft_finalize(ierr)
  use BigDFT_API
  implicit none
  integer, intent(out) :: ierr
  
  !here a routine to free the environment should be called
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call mpi_environment_free(bigdft_mpi)
   !wait all processes before finalisation
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)
end subroutine bigdft_finalize


subroutine bigdft_get_run_ids(nconfig,run_id,arr_radical,arr_posinp,ierr)
  use BigDFT_API
  use yaml_output
  implicit none
  integer, intent(in) :: nconfig !< number of configurations, if negative; starts from the radical if positive
  character(len=*), intent(in) :: run_id !< posinp file if nconfig < 0 else radical info
  character(len=*), dimension(abs(nconfig)), intent(out) :: arr_radical,arr_posinp
  integer, intent(out) :: ierr !< error code
  !local variables
  logical :: exist_list
  integer :: iconfig,nconfig_file
  integer, external :: bigdft_error_ret

  ierr=BIGDFT_SUCCESS
  !case without list_posinp file
  if (nconfig > 0) then
     if (nconfig > 1) then
        do iconfig=1,nconfig
           write(arr_radical(iconfig),'(a)')trim(run_id)//&
                trim(adjustl(yaml_toa(iconfig,fmt='(i3)')))
        end do
     else
        write(arr_radical(1),'(a)')trim(run_id)
     end if
     !input positions
     if (trim(run_id)=="input") then
        if (nconfig > 1) then
           do iconfig=1,nconfig
              write(arr_posinp(iconfig),'(a)')'posinp'//&
                   trim(adjustl(yaml_toa(iconfig,fmt='(i3)')))
           end do
        else
           write(arr_posinp(1),'(a)')'posinp'
        end if
     else
        do iconfig=1,nconfig
           write(arr_posinp(iconfig),'(a)')arr_radical(iconfig)
        end do
     end if
  !case for posinp_list file
  else
     inquire(file=trim(run_id),exist=exist_list)
     !if should be present, otherwise raise an error
     if (exist_list) then
        open(54,file=trim(run_id))
        read(54,*) nconfig_file
        if (nconfig_file /= -nconfig) then
           ierr=bigdft_error_ret(BIGDFT_INCONSISTENCY,&
                'bigdft_get_run_ids: inconsistency in list_posinp file')
           close(54)
           return
        end if
        do iconfig=1,nconfig_file
           read(54,*) arr_posinp(iconfig)
           write(arr_radical(iconfig),'(a)')arr_posinp(iconfig)
        enddo
        close(54)
     else
        ierr=bigdft_error_ret(BIGDFT_INCONSISTENCY,&
                'bigdft_get_run_ids: list_posinp file absent')
        return
     endif
  end if

end subroutine bigdft_get_run_ids


function bigdft_error_ret(err_signal,err_message) result (ierr)
  implicit none
  character(len=*), intent(in) :: err_message
  integer, intent(in) :: err_signal
  integer :: ierr

  ierr=err_signal
  
end function bigdft_error_ret

!accessors for external programs
!> Get the number of orbitals of the run in rst
function bigdft_get_number_of_atoms(atoms) result(nat)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms !> BigDFT restart variables. call_bigdft already called
  integer :: nat !> Number of atoms

  nat=atoms%astruct%nat

  if (f_err_raise(nat < 0 ,'Number of atoms unitialized')) return

end function bigdft_get_number_of_atoms

!> Get the number of orbitals of the run in rst
function bigdft_get_number_of_orbitals(rst,istat) result(norb)
  use module_base
  use module_types
  implicit none
  type(restart_objects), intent(in) :: rst !> BigDFT restart variables. call_bigdft already called
  integer :: norb !> Number of orbitals of run in rst
  integer, intent(out) :: istat

  istat=BIGDFT_SUCCESS

  norb=rst%KSwfn%orbs%norb
  if (norb==0) istat = BIGDFT_UNINITIALIZED

end function bigdft_get_number_of_orbitals

!> Fill the array eval with the number of orbitals of the last run
subroutine bigdft_get_eigenvalues(rst,eval,istat)
  use module_base
  use module_types
  implicit none
  type(restart_objects), intent(in) :: rst !> BigDFT restart variables. call_bigdft already called
  real(gp), dimension(*), intent(out) :: eval !> Buffer for eigenvectors. Should have at least dimension equal to bigdft_get_number_of_orbitals(rst,istat)
  integer, intent(out) :: istat !> Error code
  !local variables
  integer :: norb,bigdft_get_number_of_orbitals

  norb=bigdft_get_number_of_orbitals(rst,istat)

  if (istat /= BIGDFT_SUCCESS) return

  if (.not. associated(rst%KSwfn%orbs%eval)) then
     istat = BIGDFT_UNINITIALIZED
     return
  end if

  if (product(shape(rst%KSwfn%orbs%eval)) < norb) then
     istat = BIGDFT_INCONSISTENCY
     return
  end if

  call vcopy(norb,rst%KSwfn%orbs%eval(1),1,eval(1),1)

end subroutine bigdft_get_eigenvalues

subroutine bigdft_severe_abort()
  use module_base
  implicit none
  integer :: ierr

  !the MPI_ABORT works only in MPI_COMM_WORLD
  call f_malloc_dump_status()
  !call f_dump_last_error()
  call f_dump_all_errors()
  call MPI_ABORT(MPI_COMM_WORLD,816437,ierr)
  if (ierr/=0) stop 'Problem in MPI_ABORT'

end subroutine bigdft_severe_abort
