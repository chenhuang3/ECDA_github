!
! Based on the fftw3_c2r_op() function in ABINIT
! 
! AIM: 
!
!   In ABINIT, array in Fourier space is defined as work(2,nfft), which
!   is not the input for FFTW3, which uses work_fft(2,(n1/2+1)*n2*n3))
!   This is due to the fact, work array is corresponding to the Fourier transformation 
!   of the real data, therefore half of work() array is complex conjugate of 
!   the other half of work() array.
!
! INPUT: 
!
!   abinit_arr:  Fourire transformtion in q-space from ABINIT
!   nx,ny,nz  :  grid mesh in real space
!
! OUTPUT:
!   g_arr_fftw:  g_arr_fftw(2,(n1/2+1)*n2*n3)) to be provided to FFTW3 for transformation to real space.
!

subroutine array_abinit2fftw3(abinit_arr,nx,ny,nz,g_fftw)
 
  implicit none 

  integer :: nx,ny,nz,m1,m2,m3
  real(8) :: abinit_arr(2,nx*ny*nz)         ! ABINIT array, The complex array to be transformed.
  complex(8) :: g_fftw((nx/2+1)*ny*nz)      ! The backwards real array from FFTW3

  ! local vars 
  !scalars
  integer   :: ldx,ldy,ldz, j, & 
               nhp,my_flags,padx,i2,i3,igp,igf,idat,padatf,padatp,idist,odist,stride
  !arrays
  real(8) :: im,re

  ! Fill the Hermitian part: Hermitian redundancy: out[i] is the conjugate of out[n-i]

  padx = (nx/2+1)

  do i3=1,nz
    do i2=1,ny

      igf = (i3-1)*nx*ny   + (i2-1)*nx  
      igp = (i3-1)*padx*ny + (i2-1)*padx 

      do j = 1,padx
        im = abinit_arr(1,igf+j)
        re = abinit_arr(2,igf+j)
        g_fftw(igp+j) = cmplx(im,re)
      enddo

    end do
  end do

end subroutine array_abinit2fftw3
