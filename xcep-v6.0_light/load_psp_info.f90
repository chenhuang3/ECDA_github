subroutine load_psp_info(npsp, psps_file) 

 use interface_funcs
 use mpi

 implicit none 

 character(len=pspfile_len),dimension(max_psp_file) :: psps_file
 character(len=100) :: string
 integer            :: npsp, file_unit=111, ii, iost


 integer :: myrank, ierr
 call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, ierr )


 open(file='param.in',unit=file_unit,action='read',iostat=iost)

 do while (.true.) 
    read(file_unit,*,iostat=iost) string 
    if (iost<0) exit                 ! end of the file 
    if (string(1:9)=='psps_file' ) then 
      do ii=1,npsp
        read(file_unit,*) psps_file(ii)
        if (myrank==0) print *,'psp filename is ',adjustl(trim(psps_file(ii)))
      enddo
    endif
 enddo 


 close(file_unit)


end subroutine load_psp_info
