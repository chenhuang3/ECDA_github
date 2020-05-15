integer, parameter :: MPI_COMM_NULL=2,MPI_SUCCESS=0,MPI_SUM=1, MPI_COMM_WORLD=1
integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1
integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1
integer, parameter :: MPI_MIN=1, MPI_MAX=1, MPI_CHARACTER=1, MPI_REAL8=1
integer, parameter :: MPI_MAX_PROCESSOR_NAME=10, MPI_STATUS_SIZE=1,MPI_LAND=1
integer, parameter :: MPI_REQUEST_NULL=1,MPI_STATUS_IGNORE=1
integer, parameter :: mpi_tag_ub=1,mpi_address_kind=8,mpi_info_null=0
integer, parameter :: mpi_mode_noprecede=0
real(kind=8), external :: mpi_wtime
