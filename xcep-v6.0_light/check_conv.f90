subroutine  check_conv(iopt,ucvol,npt,error,error2,rho_tot_old,rho_tot_new,thres,tolv,bExit)

 implicit none 

 logical :: bExit 
 integer :: tolv,&  ! tolv>0, convergence is based on potential 
                    ! tolv<0, convergence is based on density
            iopt,&  ! iteratoin number 
            npt

 real(kind=8) :: & 
    ucvol,  & 
    error, & 
    error2, & 
    thres, & 
    rho_tot_new(npt), & 
    rho_tot_old(npt)

 print *,''    
 print *,' ----------------- check_conv.f90 ---------------------'    
 print *,''    

 if (tolv>0) write(6,'(a,es12.4)')'Convergence is based on potential, with threshold: ', thres
 if (tolv<=0)write(6,'(a,es12.4)')'Convergence is based on density, with threshold: ', thres

 error2 = error
 error = sqrt(sum((rho_tot_old-rho_tot_new)**2)*ucvol/dble(npt))

 if (iopt>=2) then 
   write(6,*)'min/max(rho_tot_old-rho_tot_new):',minval(rho_tot_old-rho_tot_new),maxval(rho_tot_old-rho_tot_new)
   if (iopt==2) write(6,'(a, es12.4)')'last iteration converges to    : ', error
   if (iopt>=3) write(6,'(a,2es12.4)')'last two iterations converge to: ', error, error2
 endif 

 bExit = .false.
 if (iopt>=2) then 
   if (error < thres .and. error2 < thres) bExit = .true.
 endif

 rho_tot_old = rho_tot_new
end subroutine 
