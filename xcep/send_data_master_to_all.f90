!
! send result on master to all the other processors using mpi_allreduce()
!
subroutine send_data_master_to_all(nsize,myrank,arr)

 use mpi
 implicit none 

 integer :: ierr,nsize, myrank
 real(8) :: arr(nsize)
 if (myrank/=0) arr = 0.d0
 call mpi_allreduce(MPI_IN_PLACE,arr,nsize,mpi_double_precision,mpi_sum, mpi_comm_world, ierr)
 if (ierr /=0) then 
   print *,'ierr/=0,  mpi_allreduce() error in main().'
   print *,'ierr: ',ierr
   stop
 endif 

end subroutine send_data_master_to_all
