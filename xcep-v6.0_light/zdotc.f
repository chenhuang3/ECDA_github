       COMPLEX*16 FUNCTION zdotc(N,ZX,INCX,ZY,INCY)
       INTEGER incx,incy,n
       COMPLEX*16 zx(*),zy(*)
       COMPLEX*16 ztemp
       INTEGER i,ix,iy
       INTRINSIC dconjg
       ztemp = (0.0d0,0.0d0)
       zdotc = (0.0d0,0.0d0)
       IF (n.LE.0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
! *
! *        code for both increments equal to 1
! *
          DO i = 1,n
             ztemp = ztemp + dconjg(zx(i))*zy(i)
          END DO
       ELSE
! *
! *        code for unequal increments or equal increments
! *          not equal to 1
! *
          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
             ztemp = ztemp + dconjg(zx(ix))*zy(iy)
             ix = ix + incx
             iy = iy + incy
          END DO
       END IF
       zdotc = ztemp
       RETURN
       END
