subroutine  wrap_rhorsp(ndim,npad,aa,bb,cc,alphaQ,betaQ)

  implicit none 

  integer :: isp, ndim(3), npad(3), itmp, & 
             shift1, shift2, shift3, i, j, k, p1, p2,p3, gms_nfft
  character(len=300) :: stemp
  real(kind=8) :: rhor(ndim(1)*npad(1),ndim(2)*npad(2),ndim(3)*npad(3)), & 
                  wrap_rhor(ndim(1),ndim(2),ndim(3)),dtmp, &
                  aa(3),bb(3),cc(3), det, dv, alphaQ, betaQ


  gms_nfft = ndim(1)*npad(1)*ndim(2)*npad(2)*ndim(3)*npad(3)

  det = aa(1)*bb(2)*cc(3)+aa(2)*bb(3)*cc(1)+aa(3)*bb(1)*cc(2)-& 
        aa(3)*bb(2)*cc(1)-aa(1)*bb(3)*cc(2)-aa(3)*bb(1)*cc(3)

  det = det * npad(1) * npad(2) * npad(3)

  dv = abs(det)/dble( gms_nfft )

  print *,'# of grids in GMS : ',gms_nfft
  print *,'volume (bohr)     : ',det
  print *,'dv (bohr)         : ',dv

  do isp=1,2
     
     ! ---------- load rhor files from GMS -------

     if (isp==1) then 
       write(stemp,'(a)')"/home/chenh/work/gamess_spinden_NewEMB/scr/mcscf_ci/RHOR_ALPH.dat"
       open(file=stemp,unit=111,action='read',form='formatted',status='old')
       print *,'loading RHOR_ALPH.dat ...'
     endif
     if (isp==2) then
       write(stemp,'(a)')"/home/chenh/work/gamess_spinden_NewEMB/scr/mcscf_ci/RHOR_BETA.dat"
       open(file=stemp,unit=111,action='read',form='formatted',status='old')
       print *,'loading RHOR_BETA.dat ...'
     endif
     do k=1,ndim(3)*npad(3)
       do j=1,ndim(2)*npad(2)
         do i=1,ndim(1)*npad(1)
           read(111,*)itmp,dtmp,dtmp,dtmp,rhor(i,j,k)
         enddo
       enddo
     enddo
     close(111)

     !---------- wrap to unit cell --------
     
     wrap_rhor = 0.d0

     if (isp==1) print *,'wrapping rhor_alpha to unit cell ...'
     if (isp==2) print *,'wrapping rhor_beta  to unit cell ...'

     do p3=1,npad(3)
       do p2=1,npad(2)
         do p1=1,npad(1)

            ! ---------- loop over wrap rho array -------
            do k=1,ndim(3)
              do j=1,ndim(2)
                do i=1,ndim(1)
                  shift1 = (p1-1)*ndim(1)
                  shift2 = (p2-1)*ndim(2)
                  shift3 = (p3-1)*ndim(3)
                  wrap_rhor(i,j,k)=wrap_rhor(i,j,k)+rhor(i+shift1,j+shift2,k+shift3)
                enddo
              enddo
            enddo

         enddo
       enddo
     enddo

     !--------- store to disk -------------
     print *,'min/max(wrap_rhor): ', minval(wrap_rhor),maxval(wrap_rhor)
     if (isp==1) then 
        print *,'Q_alpha = ',sum(wrap_rhor)*dv, ' is normalized to: ',alphaQ
        wrap_rhor = wrap_rhor / (sum(wrap_rhor)*dv) * alphaQ
        open(file='wrap_rho_alph.dat',unit=111,action='write',form='unformatted')
        print *,'rhor is wrapped to unit cell and save to wrap_rho_alph.dat'
     endif
     if (isp==2) then 
        print *,'Q_beta  = ',sum(wrap_rhor)*dv, ' is normalized to: ',betaQ
        open(file='wrap_rho_beta.dat',unit=111,action='write',form='unformatted')
        print *,'rhor is wrapped to unit cell and save to wrap_rho_beta.dat'
        wrap_rhor = wrap_rhor / (sum(wrap_rhor)*dv) * betaQ
     endif
     write(111) wrap_rhor
     close(111)
  enddo



end subroutine 
