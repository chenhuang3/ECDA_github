!!@LICENSE
!
!     MODULE m_io
!
! Copyright Alberto Garcia, 1996, 1997, 1998
!
! This module implements an interface to the FORTRAN logical unit
! system. Based on code by Richard Maine.
!
! Alberto Garcia, December 30, 1996
! Rewritten as a single subroutine 
! with multiple entry points, March 7, 1998
! Converted to a module by J.M.Soler. Aug. 2009
!---------------------------------------------------------------
!
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


      MODULE m_io
!
!-----------------------------------------------------------------
!
!     Used module procedures
!
      use m_errors,   only : die

      implicit none
!

!-----------------------------------------------------------------
!
!     Public procedures provided by this module
!
      public :: &
     &  io_seterr,  & ! Set standard error unit
     &  io_setout,  & ! Set standard output unit
     &  io_geterr,  & ! Get standard error unit
     &  io_getout,  & ! Get standard output unit
     &  io_assign,  & ! Get some available IO unit and reserve it
     &  io_reserve, & ! Reserve a specific IO unit
     &  io_close,   & ! Close and free a given IO unit
     &  io_status     ! Print all used IO units

! Nothing is declared public below this point
     private
!
!----------------------------------------------------------------
!
!     Module variables
!
!     Logical unit management. Units 0 to min_lun-1 are "reserved",
!     since most of the "typical" files (output, etc) use them.
!
!     Logical units min_lun to min_max are managed by this module.
!
      integer, parameter:: min_lun = 10
      integer, parameter:: max_lun = 99
      integer, parameter:: nunits = max_lun-min_lun+1
      integer, save:: stdout = 6
      integer, save:: stderr = 0
      logical, save:: lun_is_free(min_lun:max_lun) = .true.
!
!-----------------------------------------------------------------
!
!     Internal and dummy variables
!
      integer  :: i, iostat
      logical  :: used, named, opened
      character:: filename*50, form*11
!
      CONTAINS
!
!-----------------------------------------------------------------
!
!     Simple interfaces to modify standard units
!
      subroutine io_seterr(unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io_seterr'
!End of the abilint section

      integer,intent(in):: unit
      stderr = unit
      end subroutine io_seterr
!
!-----------------------------------------------------------------
!
      subroutine io_setout(unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io_setout'
!End of the abilint section

      integer,intent(in):: unit
      stdout = unit
      end subroutine io_setout
!
!-----------------------------------------------------------------
!
      subroutine io_geterr(unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io_geterr'
!End of the abilint section

      integer,intent(out):: unit
      unit = stderr
      end subroutine io_geterr
!
!-----------------------------------------------------------------
!
      subroutine io_getout(unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io_getout'
!End of the abilint section

      integer,intent(out):: unit
      unit = stdout
      end subroutine io_getout
!
!------------------------------------------------------------------     
!
!     Logical unit management
!
      subroutine io_assign(lun)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io_assign'
!End of the abilint section

      integer,intent(out):: lun
!
!     Looks for a free unit and assigns it to lun
!
      do lun= min_lun, max_lun
         if (lun_is_free(lun)) then
            inquire(unit=lun, opened=used, iostat=iostat)
            if (iostat .ne. 0) used = .true.
            lun_is_free(lun) = .false.
            if (.not. used) return
         endif
      enddo
      call die("No luns available in io_assign")

      end subroutine io_assign
!
!------------------------------------------------------------------     
!
      subroutine io_reserve(lun)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io_reserve'
!End of the abilint section

      integer,intent(in):: lun
!
!     Useful to specify that one needs to use a particular unit number
!
!     For example, assume some legacy code expects to work with unit 15:
!
!     call io_reserve(15)   ! this call at the beginning of the program
!     ...
!     open(15,....)
!
      inquire(unit=lun, opened=used, iostat=iostat)
      if (iostat .ne. 0) used = .true.
      if (used) then
        call die("Cannot reserve unit. Already connected")
      end if
      if (lun .ge. min_lun .and. lun .le. max_lun)     &
     &                      lun_is_free(lun) = .false.

      end subroutine io_reserve
!
!------------------------------------------------------------------     
!
      subroutine io_close(lun)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io_close'
!End of the abilint section

      integer,intent(in):: lun
!
!     Use this routine instead of a simple close!!
!
      close(lun)
      if (lun .ge. min_lun .and. lun .le. max_lun)   &
     &                     lun_is_free(lun) = .true.

      end subroutine io_close
!
!------------------------------------------------------------------     
!
      subroutine io_status
!
!     Prints a list of the connected logical units and the names of
!     the associated files
!

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'io_status'
!End of the abilint section

      write(stdout,'(a)') '******** io_status ********'
      do i = 0, max_lun
         inquire(i,opened=opened,named=named,name=filename,   &
     &           form=form,iostat=iostat)
         if (iostat .eq. 0) then
            if (opened) then
               if (named) then
                  write(stdout,9000) i, form, filename
               else
                  write(stdout,9000) i, form, 'No name available'
               endif
            endif
         else
            write(stdout,9000) i, 'Iostat error'
         endif
      enddo
      write(stdout,'(a)') '********           ********'

 9000 format(i4,5x,a,5x,a)
      end subroutine io_status

      END MODULE m_io

