function check_nan(n,var)

  logical :: check_nan
  integer :: n, i
  real(8) :: var(n)

  check_nan = .false.
  do i=1,n
    if (isnan(var(i))) then 
      check_nan = .true. 
      return 
    endif 
  enddo

  return  

end function 
