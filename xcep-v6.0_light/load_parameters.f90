!
! load parameters from param.in file
!
!  Chen Huang , Jan/2011
!
subroutine load_parameters(restart,file_dfet_out,natom,nspin,npsps,cell_nfft,nsys,global_fermi,mix_gvks, & 
                           start_u,dft_flavor,nopt_uemb,npulay_gvks,pulay_gvks_beta,cgW_tol,&
                           cell_acell,ucvol,rprimd,rmet,tsmear,relax_geo,only_lda, & 
                           reg_type,zp_alpha,reg_vxc,zp_coeff,zp_conv,vac,plot_sub_rho,relax_nscf, & 
                           exx_method,hse_range)

 use mpi 
 use comm_data

 implicit none 
 
 ! external vars
 logical              :: plot_sub_rho,only_lda
 integer,parameter    :: dp=8
 integer              :: restart,nopt_uemb,global_fermi,mix_gvks,exx_method 
 integer              :: file_dfet_out,dft_flavor,relax_geo, & 
                         natom,nspin,cell_nfft(3),nsys,start_u,npsps,npulay_gvks, & 
                         relax_nscf,  reg_type 
 real(dp),intent(out) :: ucvol,cell_acell(3,3),rprimd(3,3),rmet(3,3), cgW_tol, & 
                         pulay_gvks_beta,tsmear, hse_range, & 
                         zp_alpha, reg_vxc, zp_coeff, zp_conv, vac

 ! local vars
 logical :: set_qalpha=.false., & 
            set_qbeta=.false., &
            set_startu = .false., & 
            set_acell = .false.

 real(dp) :: dtmp
 integer :: iost, file_unit, myrank, ierr,itmp
 character(len=100)  :: string
 integer,allocatable :: ntmp(:)


 !=============== default values ===========
 precond_vks = .false.
 do_envOF = .false.
 only_lda = .false. 
 vac = 0.5d0 ! default value 
 plot_sub_rho = .false.
 relax_nscf = -1
 exx_method = -99
 reg_type = -1 
 zp_coeff = -1
 zp_conv  = -1
 tsmear = -1
 file_unit = 111
 cgW_tol = -1.0d0
 npulay_gvks=-1
 mix_gvks = -1;
 pulay_gvks_beta=-1.d0
 ucvol = -100.0d0
 npsps = -100
 nsys = -100
 natom = -100
 nspin = -100
 reg_vxc = -100
 global_fermi = -999
 nopt_uemb = -1 
 restart = 999

 call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, ierr )

 open(file='param.in',unit=file_unit,action='read',iostat=iost)
 if (iost/=0) then
   print *,'error on opening param.in file, iostat=',iost
   stop
   close(file_unit)
 endif


 do while (.true.) 
   read(file_unit,*,iostat=iost) string 
   if (iost<0) exit                 ! end of the file 
   if (string(1:1) .eq. '#') cycle  ! this line is comment 
   if (string(1:6) .eq. 'vw_lam' ) then 
     backspace file_unit
     read(file_unit,*) string, vw_lam
     write(file_dfet_out,'(a,f12.4)')'[param] vw_lam:        ',vw_lam
   endif
   if (string(1:7) .eq. 'nsubsys' ) then 
     backspace file_unit
     read(file_unit,*) string, nsys
   endif
   if (string(1:11) .eq. 'precond_vks' ) then 
     backspace file_unit
     read(file_unit,*) string,dtmp
     if (dtmp>0) precond_vks = .true.
     write(file_dfet_out,'(a)')'[param] precond_vks is TRUE'
   endif
   if (string(1:9) .eq. 'env_OFDFT' ) then 
     backspace file_unit
     read(file_unit,*) string,dtmp
     if (dtmp>0) do_envOF = .true. 
   endif
   if (string(1:10) .eq. 'npsps_file' ) then 
     backspace file_unit
     read(file_unit,*) string,npsps
   endif
   if (string(1:8) .eq. 'zp_alpha' ) then 
     backspace file_unit
     read(file_unit,*) string,zp_alpha
   endif
   if (string(1:8) .eq. 'reg_type' ) then 
     backspace file_unit
     read(file_unit,*) string,reg_type
     write(file_dfet_out,'(a,i4)')'[param] reg_type:      ',reg_type
   endif
   if (string(1:10) .eq. 'exx_method' ) then 
     backspace file_unit
     read(file_unit,*) string,exx_method,hse_range
     write(file_dfet_out,'(a,i4)')     '[param] exx_method:',exx_method
     write(file_dfet_out,'(a,f12.4,a)')'[param] rc:        ',hse_range,' bohr'
     ! the input is the inverse of the cutoff length 
     hse_range = 1.d0/hse_range 
   endif
   if (string(1:7) .eq. 'reg_vxc' ) then 
     backspace file_unit
     read(file_unit,*) string,reg_vxc
     write(file_dfet_out,'(a,es12.4)')'[param] reg_vxc:    ',reg_vxc
   endif
   if (string(1:10) .eq. 'restart' ) then 
     backspace file_unit
     read(file_unit,*) string,restart
   endif
   if (string(1:10) .eq. 'dft_flavor' ) then 
     backspace file_unit
     read(file_unit,*) string,dft_flavor
     write(file_dfet_out,'(a,i4,a)')&
     '[param] dft_flavor:   ', dft_flavor ,' '
   endif
   if ((string(1:7) .eq. 'cgw_tol') .or. (string(1:7) .eq. 'cgW_tol')  ) then 
     backspace file_unit
     read(file_unit,*) string, cgW_tol
     write(file_dfet_out,'(a,es12.4)') &
       '[param] cgW_tol:    ',cgW_tol
   endif
   if (string(1:12) .eq. 'global_fermi' ) then
     backspace file_unit
     read(file_unit,*) string, global_fermi
     write(file_dfet_out,'(a,i3)')     &
       '[param] global_fermi: ',global_fermi
   endif
   if (string(1:8) .eq. 'mix_gvks' ) then
     backspace file_unit
     read(file_unit,*) string,mix_gvks
   endif
   if (string(1:15) .eq. 'pulay_gvks_beta' ) then
     backspace file_unit
     read(file_unit,*) string, pulay_gvks_beta
     write(file_dfet_out,'(a,f12.4)') & 
       '[param] pulay_gvks_beta: ',pulay_gvks_beta
   endif
   if (string(1:11) .eq. 'npulay_gvks' ) then
     backspace file_unit
     read(file_unit,*) string,npulay_gvks
     write(file_dfet_out,'(a,i3)')'[param] npulay_gvks:',npulay_gvks
   endif
   if (string(1:11) .eq. 'relax_nscf' ) then
     backspace file_unit
     read(file_unit,*) string, relax_nscf
     write(file_dfet_out,'(a,i3)')'[param] relax_nscf:', relax_nscf
   endif
   if (string(1:8) .eq. 'zp_coeff' ) then
     backspace file_unit
     read(file_unit,*) string,zp_coeff
     write(file_dfet_out,'(a,f12.4)')'[param] zp_coeff:',zp_coeff
   endif
   if (string(1:8) .eq. 'zp_conv' ) then
     backspace file_unit
     read(file_unit,*) string,zp_conv
     write(file_dfet_out,'(a,es12.4)')'[param] zp_conv:',zp_conv
   endif
   if (string(1:6) .eq. 'tsmear' ) then
     backspace file_unit
     read(file_unit,*) string, tsmear
     write(file_dfet_out,'(a,f12.4,a)')'[param] tsmear: ',tsmear,' eV'
     tsmear = tsmear / 27.2114
   endif
   if (string(1:5) .eq. 'natom' ) then 
     backspace file_unit
     read(file_unit,*) string,natom
   endif
   if (string(1:7) .eq. 'start_u' ) then 
     backspace file_unit
     set_startu=.true.
     read(file_unit,*) string,start_u
   endif
   if (string(1:9) .eq. 'relax_geo' ) then 
     backspace file_unit
     read(file_unit,*) string,relax_geo
   endif
   if (string(1:15) .eq. 'cell_acell_bohr') then 
     set_acell = .true.
     backspace file_unit
     read(file_unit,*) string
     read(file_unit,*) cell_acell(1,:)
     read(file_unit,*) cell_acell(2,:)
     read(file_unit,*) cell_acell(3,:)
   endif 
   if (string(1:16) .eq. 'cell_acell_angst') then 
     set_acell = .true.
     backspace file_unit
     read(file_unit,*) string
     read(file_unit,*) cell_acell(1,:)
     read(file_unit,*) cell_acell(2,:)
     read(file_unit,*) cell_acell(3,:)
     cell_acell = cell_acell / 0.5291772106D0
   endif 
   if (string(1:9) .eq. 'nopt_uemb') then 
     backspace file_unit
     read(file_unit,*) string,nopt_uemb
     if (myrank==0)  print *,'# of iterations for optimizing embedding potential',nopt_uemb
   endif 
   if (string(1:12) .eq. 'plot_sub_rho') then 
     backspace file_unit
     read(file_unit,*) string,dtmp
     plot_sub_rho = .false.
     if (dtmp>0) plot_sub_rho = .true.
   endif 
   if (string(1:8) .eq. 'only_lda') then 
     backspace file_unit
     read(file_unit,*) string,dtmp
     if (dtmp>0) only_lda = .true.
   endif 
   if (string(1:3) .eq. 'vac') then 
     backspace file_unit
     read(file_unit,*) string,vac
   endif 
   if (string(1:9) .eq. 'cell_nfft') then 
     backspace file_unit
     read(file_unit,*) string,cell_nfft(1:3)
   endif 
   if (string(1:5) .eq. 'nspin') then 
     backspace file_unit
     read(file_unit,*) string, nspin
   endif
   if (string(1:10) .eq. 'share_atom') then 
     backspace file_unit
     read(file_unit,*) string,dtmp
     if (dtmp>0) share_atom = .true.
   endif
 enddo
 write(file_dfet_out,'(a,L)')'[param] plot_sub_rho:   ',plot_sub_rho
 write(file_dfet_out,'(a,f12.4)') '[param] vac:   ',vac
 write(file_dfet_out,'(a,es12.4)')'[param] zp_alpha:   ',zp_alpha

 close(file_unit)

 ! check 
 ! =====
 if (cgw_tol < 0.0d0 ) then 
    print *,'give cgw_tol in param.ini file! stop'
    stop
 endif
 if (npulay_gvks<0 .or. mix_gvks<0) then 
   print *,'npulay_gvks or mix_gvks are not given. stop '
   stop
 endif
 if ( mix_gvks /=1 .and. mix_gvks/=2) then 
   print *,'undefined mix_gvks value: ',mix_gvks
   stop
 endif
 if (nspin<0) then 
   print *,'no nspin found '
   stop
 endif
 if (tsmear<0) then 
   print *,'tsmear is not given'
   stop
 endif
 if (natom<0) then 
   print *,'no natom key found in param.in'
   stop
 endif
 if (nsys<0) then 
   print *,'no nsys paramter found, stop'
   stop
 endif
 if (npsps<0) then 
   print *,'no npsps found, stop!'
   stop
 endif
 if (zp_coeff<0) then 
   print *,'zp_coeff is not set, stop'
   stop
 endif 
 if (zp_conv<0) then 
   print *,'zp_conv is not set, stop'
   stop
 endif 
 if (relax_nscf<0) then 
   print *,' relax_nscf shoud be given, stop'
   stop
 endif 
 if (pulay_gvks_beta<0.d0) then 
   print *,'no pulay_gvks_beta is given, stop'
   stop
 endif
 if ( global_fermi == -999 ) then 
   print *,'no global_fermi is given, stop'
   stop
 endif
 if (restart==999) then 
   print *,'restart keyword is not specifiied! '
   stop
 endif 
 if (reg_vxc<0) then 
   print *,'reg_vxc is not given, stop'
   stop
 endif 
 if (reg_type<0) then 
   print *,'reg_type is not given, stop'
   stop
 endif 
 if (exx_method==-99) then 
   print *,'exx_method is not given'
   stop
 endif 


 rprimd = transpose(cell_acell)

 ucvol=rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3))+&
& rprimd(2,1)*(rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3))+&
& rprimd(3,1)*(rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3))

! rprimd is in column wise from ABINIT's webste 

 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)  ! ABINIT's metric.F90 file

 if (myrank==0 ) then 
 print *,''
 print *,'Parameters'
 print *,'----------'
 print *,' '
 write(6,'(a,i4)')       '  nspin         : ',nspin
 write(6,'(a,i4)')       '  nsubsystem    : ',nsys
 write(6,'(a,i4)')       '  start_u       : ',start_u
 write(6,'(a,3i4)')      '  nfft_mesh     : ',cell_nfft
 write(6,'(a,3f12.4,a)') '  latt. box  a= : ',cell_acell(1,:),' bohr'
 write(6,'(a,3f12.4,a)') '             b= : ',cell_acell(2,:),' bohr'
 write(6,'(a,3f12.4,a)') '             c= : ',cell_acell(3,:),' bohr'
 write(6,'(a,f12.4,a)')  '  ucvol         : ',ucvol,' bohr^3'
 write(6,'(a,f12.4,a)')  '  ucvol         : ',ucvol*0.52917721067**3,' angstr^3'
 print *,""
 print *,""
 endif

 return 
end subroutine load_parameters

