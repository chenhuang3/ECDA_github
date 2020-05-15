!
! make the prime vector in the reciprocal space.
! 
! refer to metric.F90 of ABINIT 
!
! Note:
!   My gprim = 2*pi*grimd_ABINIT
!   My gmet = (2*pi)^2*gmet_ABINIT
!
! Compute first dimensional primitive translation vectors in reciprocal space
! gprim from rprimd, and eventually writes out.
! Then, computes metrics for real and recip space rmet and gmet using length
! dimensional primitive translation vectors in columns of rprimd(3,3) and gprim(3,3).
! gprim is the inverse transpose of rprimd.
!       $ gmet_{i,j}= \sum_k ( gprim_{k,i}*gprim_{k,j} )  $
!
!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!  gprim(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!

subroutine make_gprim(cell_acell,gprim,gmet,abinit_gmet)

 use mpi 


 implicit none

 integer :: myrank,  ierr
 real(8),parameter :: pi = 3.14159265359d0
 real(8) :: cell_acell(3,3)
 real(kind=8) :: acell(3,3),res(3,3)
 real(kind=8) :: gprim(3,3),gmet(3,3),abinit_gmet(3,3)



 call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)


 ! this is the OFDFT's cellreal definition
 ! acell  = [a|b|c]
 if (myrank==0) then 
   print *,''
   print *,'only consider orthogonal box now'
   print *,''
 endif

 !
 ! cell_acell is in row-wise for a b c
 !
 acell = transpose(cell_acell)

 
 ! then |g1|
 !      |g2| * [a|b|c] = 2PI
 !      |g3|
 ! finally , we want [g1|g2|g3] in column wise
 ! this is the OFDFT's definition
 
 call inverse3(acell,res)
 gprim = TRANSPOSE(res) * 2.d0 * pi

 if (myrank==0) then 
   print *,''
   print *,'gprim in q-space :'
   write(6,'(a,3es18.8,a)') '  gprim_a = ',gprim(:,1),' 1/bohr'
   write(6,'(a,3es18.8,a)') '  gprim_b = ',gprim(:,2),' 1/bohr'
   write(6,'(a,3es18.8,a)') '  gprim_c = ',gprim(:,3),' 1/bohr'
   print *,''
 endif


 !Compute reciprocal space metric.
 gmet = MATMUL(TRANSPOSE(gprim),gprim)

 abinit_gmet = gmet / (2.d0 * pi)**2

 if (myrank==0) then 
   print *,''
   print *,'gmet:'
   write(6,'(a,3es18.8,a)') ' ',gmet(:,1),' '
   write(6,'(a,3es18.8,a)') ' ',gmet(:,2),' '
   write(6,'(a,3es18.8,a)') ' ',gmet(:,3),' '
   print *,''
 endif 
 
 return
end subroutine make_gprim
