!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_haydock
!! NAME
!! m_haydock
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2014 ABINIT group (M.Giantomassi, Y. Gillet, L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_haydock

 use defs_basis
 use m_profiling_abi
 use m_bs_defs
 use m_xmpi
 use m_errors
 use m_ncfile
 use m_haydock_io
 use m_linalg_interfaces
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif

 use m_fstrings,          only : indent, strcat, sjoin, itoa
 use m_io_tools,          only : file_exists, open_file
 use defs_abitypes,       only : Hdr_type
 use defs_datatypes,      only : ebands_t, pseudopotential_type
 use m_geometry,          only : normv
 use m_blas,              only : xdotc, xgemv
 use m_numeric_tools,     only : print_arr, symmetrize, hermitianize, continued_fract, wrap2_pmhalf
 use m_fft_mesh,          only : calc_ceigr
 use m_crystal,           only : crystal_t 
 use m_crystal_io,        only : crystal_ncwrite
 use m_bz_mesh,           only : kmesh_t, findqg0, get_bz_item
 use m_double_grid,       only : double_grid_t, get_kpt_from_indices_coarse, compute_corresp
 use m_interp,            only : interpolator_t, interpolator_init, interpolator_free, interpolator_normalize, INTERP_YG
 use m_paw_commutator,    only : HUr_commutator
 use m_wfs,               only : wfd_t, wfd_sym_ur, wfd_get_ur, wfd_change_ngfft
 use m_bse_io,            only : exc_read_rcblock, exc_write_optme
 use m_pawtab,            only : pawtab_type
 use m_vcoul,             only : vcoul_t

 implicit none

 private 
!!***

 public :: exc_haydock_driver     ! Driver for the Haydock method (main entry point for client code). 

!----------------------------------------------------------------------

!!!!!****t* m_haydock/hexc_t
!!!!! NAME
!!!!! hexc_t
!!!!! 
!!!!! FUNCTION
!!!!! 
!!!!! SOURCE
!!!
!!! type,private :: hexc_t
!!!
!!!   integer :: nsppol 
!!!   ! Number of spins
!!!
!!!   !integer :: hsize_coarse(nsppol)
!!!   !integer :: hsize_dense(nsppol)
!!!
!!!   integer :: nvert=8
!!!   ! Number of vertices for interpolation.
!!!
!!!   integer :: comm
!!!   ! MPI communicator
!!!
!!!   integer,allocatable :: corresp(:,:)
!!!   ! corresp(hsize,8)
!!!   ! mapping between coarse points and neighbours
!!!
!!!   integer,allocatable :: indices(:,:)
!!!   ! indices(BSp%nreh(1),grid%ndiv)
!!!   ! Table (t_coarse, index_in_fine_box) --> index of the k-point in Trans_interp
!!!
!!!   real(dp),allocatable :: interp_factors(:,:,:)
!!!   ! interp_factors(BSP%nreh(spin),8,grid%ndiv)
!!!   ! Interpolation factors.
!!!
!!!   complex(dpc),allocatable :: all_hmat(:,:,:)
!!!   ! all_hmat,(hsize,hsize,8))
!!!   ! Coarse excitonic matrix in a format suitable for interpolation in k-space
!!!
!!!   complex(dpc),allocatable :: all_acoeffs(:,:,:)
!!!   ! all_acoeffs(hsize,hsize,8))
!!!   ! a coefficients in a format suitable for interpolation in k-space
!!!
!!!   complex(dpc),allocatable :: all_bcoeffs(:,:,:)
!!!   ! all_bcoeffs(hsize,hsize,8))
!!!   ! b coefficients in a format suitable for interpolation in k-space
!!!
!!!   complex(dpc),allocatable :: all_ccoeffs(:,:,:)
!!!   ! all_ccoeffs(hsize,hsize,8))
!!!   ! c coefficients in a format suitable for interpolation in k-space
!!!
!!!   complex(gwpc),allocatable :: overlaps(:,:,:) 
!!!   ! overlaps,(BSp%nreh_interp(spin1),BSp%maxnbndv*BSp%maxnbndc,8))
!!!   ! Overlap matrix
!!!
!!!   ! Pointers to datatypes that are already in memory.
!!!   type(excparam),pointer :: bsp => null()
!!!   type(crystal_t),pointer :: crystal => null()
!!!   type(kmesh_t),pointer :: kmesh_dense => null()
!!!   type(kmesh_t),pointer :: kmesh_coarse => null()
!!!   type(double_grid_t),pointer :: double_grid => null()
!!!   type(vcoul_t),pointer :: vcp_dense => null()
!!! end type hexc_t
!!!
!!! !public :: hexc_init       ! Construct the object.
!!! !public :: hexc_free       ! Release memory.
!!! !public :: hexc_build      ! Interpolate the Hamiltonian and store it in memory
!!! !public :: hexc_gemv       ! Matrix-vector multiplication.
!!!!!***

CONTAINS  !=======================================================================
!!***

!!****f* m_haydock/exc_haydock_driver
!! NAME
!! exc_haydock_driver
!!
!! FUNCTION
!!  Calculate the imaginary part of the macroscopic dielectric function with the Haydock recursive method.
!!
!! INPUTS
!! BSp<type(excparam)=The parameter for the Bethe-Salpeter run.
!! BS_files<excparam>=Files associated to the bethe_salpeter code.
!! Cryst<crystal_t>=Info on the crystalline structure.
!! Kmesh<type(kmesh_t)>=The list of k-points in the BZ, IBZ and symmetry tables.
!! Cryst<type(crystal_t)>=Info on the crystalline structure.
!! Hdr_bse
!! KS_BSt=The KS energies.
!! QP_BSt=The QP energies.
!! Wfd<wfd_t>=Wavefunction descriptor (input k-mesh)
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials.
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data.
!! Hur(Cryst%natom*usepaw)<type(HUr_commutator)>=Only for PAW and LDA+U, quantities used to evaluate the commutator [H_u,r].
!!
!! OUTPUT
!!  The imaginary part of the macroscopic dielectric function is written on the external file _EXC_MDF
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine exc_haydock_driver(BSp,BS_files,Cryst,Kmesh,Hdr_bse,KS_BSt,QP_Bst,Wfd,Psps,Pawtab,Hur,&
& Kmesh_dense, KS_BSt_dense, QP_BSt_dense, Wfd_dense, Vcp_dense, grid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_haydock_driver'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_69_wfdesc
 use interfaces_71_bse
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(Hdr_type),intent(in) :: Hdr_bse
 type(wfd_t),intent(inout) :: Wfd
 type(pseudopotential_type),intent(in) :: Psps
 type(ebands_t),intent(in) :: KS_BSt,QP_Bst
!Interp@BSE
 type(double_grid_t),intent(in),optional :: grid
 type(kmesh_t),intent(in),optional :: Kmesh_dense
 type(wfd_t),intent(inout),optional :: Wfd_dense
 type(ebands_t),intent(in),optional :: KS_BSt_dense, QP_Bst_dense
 type(vcoul_t),intent(in),optional :: Vcp_dense
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(HUr_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: io,my_rank,iq,itt,ierr
 integer :: hsize,comm,my_t1,my_t2,nsppol,nkets,nproc
 integer :: spin,spad,ik_bz,iv,ic,trans_idx,lomo_min,max_band
 integer :: max_r,max_c
 !Interp@BSE
 integer :: hsize_dense
 real(dp) :: omegaev,rand_phi !,norm
 complex(dpc) :: ks_avg,gw_avg,exc_avg
 logical :: is_resonant,use_mpio,diago_is_real,prtdos
 character(len=500) :: msg
 character(len=fnlen) :: hreso_fname,hcoup_fname !,ome_fname
 type(ncfile_t) :: ncf
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: dos(:),dos_gw(:),dos_ks(:)
 complex(dpc),allocatable :: green(:,:),hreso(:,:),hcoup(:,:),test(:,:)
 complex(dpc),allocatable :: opt_cvk(:,:,:,:,:),kets(:,:)
 complex(dpc),allocatable :: eps_rpanlf(:,:),eps_gwnlf(:,:)
 complex(dpc),allocatable :: tensor_cart(:,:),tensor_cart_rpanlf(:,:),tensor_cart_gwnlf(:,:)
 complex(dpc),allocatable :: tensor_red(:,:),tensor_red_rpanlf(:,:),tensor_red_gwnlf(:,:)

!************************************************************************

 call timab(690,1,tsec) ! exc_haydock_driver
 call timab(691,1,tsec) ! exc_haydock_driver(read)

 if (BSp%have_complex_ene) then
   MSG_ERROR("Complex energies are not supported yet")
 end if

 my_rank = Wfd%my_rank
 comm    = Wfd%comm
 nsppol  = Wfd%nsppol
 nproc   = Wfd%nproc

 use_mpio=.FALSE.
#ifdef HAVE_MPI_IO
 use_mpio = (nproc > 1)
 !use_mpio = .TRUE. 
#endif
 use_mpio=.FALSE.
 !use_mpio = .TRUE. 

 ! Hsize refers to the size of the individual blocks (resonant and coupling). 
 ! Thanks to the symmetry property of the starting vector, the Haydock method 
 ! can be reformulated in terms of matrix-vector multiplication involving the 
 ! blocks thus avoiding to allocation of the full matrix ( R   C )
 !                                                        -C* -R*)
 hsize=SUM(BSp%nreh)

 ! Divide the columns of the Hamiltonian among the nodes.
 call xmpi_split_work(hsize,comm,my_t1,my_t2,msg,ierr)
 if (ierr/=0) then
   MSG_WARNING(msg)
 end if

 ABI_CHECK(my_t2-my_t1+1>0,"found processor with 0 rows")
                                                             
 ABI_MALLOC(hreso,(hsize,my_t1:my_t2))
 ABI_CHECK_ALLOC("out of memory in hreso")

 ! Read the resonant block from file.
 if (BS_files%in_hreso /= BSE_NOFILE) then
   hreso_fname = BS_files%in_hreso
 else 
   hreso_fname = BS_files%out_hreso
 end if

 is_resonant=.TRUE.; diago_is_real=(.not.BSp%have_complex_ene)
 call exc_read_rcblock(hreso_fname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,my_t1,my_t2,hreso,use_mpio,comm)

 if (BSp%use_interp) then
   ABI_CHECK(nsppol==1,"Interpolation not yet coded with spin")
   ! No parallelization at present ...
   hsize_dense = SUM(BSp%nreh_interp)
 end if

 !call hermitianize(hreso,"All")

!BEGIN DEBUG
 if (use_mpio) then
   MSG_WARNING("Testing MPI-IO routines")
   ABI_MALLOC(test,(hsize,my_t1:my_t2))
   ABI_CHECK_ALLOC("out of memory in hreso")
   diago_is_real=(.not.BSp%have_complex_ene)
   call exc_read_rcblock(hreso_fname,Bsp,is_resonant,diago_is_real,nsppol,Bsp%nreh,hsize,my_t1,my_t2,test,.FALSE.,comm)
   test = test-hreso
   write(std_out,*)"DEBUG: Diff MPI-IO - Fortran ",MAXVAL(ABS(test))
   max_r=20; max_c=10
   write(std_out,*)" **** Testing resonant block **** "
   call print_arr(test,max_r=max_r,max_c=max_c,unit=std_out)
   if (nsppol==2) then
     write(std_out,*)" **** D down down ****"
     call print_arr(test(hsize/2+1:,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
     write(std_out,*)" **** V up down ****"
     call print_arr(test(1:hsize/2,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
     write(std_out,*)" **** V down up ****"
     call print_arr(test(hsize/2+1:,1:hsize/2),max_r=max_r,max_c=max_c,unit=std_out)
   end if
   ABI_FREE(test)
 end if
!END DEBUG
 !
 ! Read coupling block.
 if (BSp%use_coupling>0) then 
   ABI_CHECK(.not. Bsp%use_interp,"interpolation with coupling not coded!")
   if (BS_files%in_hcoup /= BSE_NOFILE) then
     hcoup_fname = BS_files%in_hcoup
   else 
     hcoup_fname = BS_files%out_hcoup
   end if

   ABI_MALLOC(hcoup,(hsize,my_t1:my_t2))
   ABI_CHECK_ALLOC("out of memory in hcoup")
   is_resonant=.FALSE.; diago_is_real=.FALSE.
   call exc_read_rcblock(hcoup_fname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,my_t1,my_t2,hcoup,use_mpio,comm)
   !call symmetrize(hcoup,"ALL")

   if (use_mpio) then
     MSG_WARNING("Testing MPI-IO routines")
     ABI_MALLOC(test,(hsize,my_t1:my_t2))
     ABI_CHECK_ALLOC("out of memory in text")
     diago_is_real=.FALSE.
     call exc_read_rcblock(hcoup_fname,Bsp,is_resonant,diago_is_real,nsppol,Bsp%nreh,hsize,my_t1,my_t2,test,.FALSE.,comm)
     test = test-hcoup
     write(std_out,*)"DEBUG: Diff MPI-IO - Fortran ",MAXVAL(ABS(test))
     max_r=20; max_c=10
     write(std_out,*)" **** Testing coupling block **** "
     call print_arr(test,max_r=max_r,max_c=max_c,unit=std_out)
     if (nsppol==2) then
       write(std_out,*)" **** D down down ****"
       call print_arr(test(hsize/2+1:,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
       write(std_out,*)" **** V up down ****"
       call print_arr(test(1:hsize/2,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
       write(std_out,*)" **** V down up ****"
       call print_arr(test(hsize/2+1:,1:hsize/2),max_r=max_r,max_c=max_c,unit=std_out)
     end if
     ABI_FREE(test)
   end if
 end if

 call timab(691,2,tsec) ! exc_haydock_driver(read)
 call timab(692,1,tsec) ! exc_haydock_driver(prep)
 !
 ! Prepare the starting vectors for the Lanczos chain.
 nkets=Bsp%nq

 prtdos=.FALSE. !prtdos=.TRUE.
 if (prtdos) then
   nkets=nkets+1
   if (Bsp%use_coupling>0) then 
     MSG_ERROR("DOS with coupling not coded")
     nkets=nkets+1
   end if
 end if

 if (BSp%use_interp) then
   ABI_MALLOC(kets,(hsize_dense,nkets))
 else
   ABI_MALLOC(kets,(hsize,nkets))
 end if
 ABI_CHECK_ALLOC("out of memory in kets")
 kets=czero
 !
 ! Prepare the kets for the macroscopic dielectric function.
 lomo_min=Bsp%lomo_min; max_band=Bsp%nbnds

 if (BSp%use_interp) then
   ABI_MALLOC(opt_cvk,(lomo_min:max_band,lomo_min:max_band,BSp%nkbz_interp,Wfd%nsppol,BSp%nq))
 else
   ABI_MALLOC(opt_cvk,(lomo_min:max_band,lomo_min:max_band,BSp%nkbz,Wfd%nsppol,BSp%nq))
 end if
 ABI_CHECK_ALLOC("out of memory in opt_cvk")

 do iq=1,Bsp%nq

   if (BSp%use_interp) then
     ! Use dense mesh for the oscillator matrix elements.
     ! Note KS_BSt is used here to calculate the commutator.
     call calc_optical_mels(Wfd_dense,Kmesh_dense,KS_BSt_dense,Cryst,Psps,Pawtab,Hur, &
&       BSp%inclvkb,BSp%lomo_spin,lomo_min,max_band,BSp%nkbz_interp,BSp%q(:,iq),opt_cvk(:,:,:,:,iq))

     ! Fill ket0 using the same ordering for the indeces as the one used for the excitonic Hamiltonian.
     ! Note that only the resonant part is used here.
     do spin=1,nsppol
       spad=(spin-1)*BSp%nreh_interp(1)
       do ik_bz=1,BSp%nkbz_interp
         do iv=BSp%lomo_spin(spin),BSp%homo_spin(spin)
           do ic=BSp%lumo_spin(spin),BSp%nbnds
             trans_idx = BSp%vcks2t_interp(iv,ic,ik_bz,spin)
             if (trans_idx>0) kets(trans_idx+spad,iq)=opt_cvk(ic,iv,ik_bz,spin,iq)
           end do
         end do
       end do
     end do
   else
     !
     ! KS_BSt is used here to calculate the commutator.
     call calc_optical_mels(Wfd,Kmesh,KS_BSt,Cryst,Psps,Pawtab,Hur,BSp%inclvkb,Bsp%lomo_spin,lomo_min,max_band,&
&                           BSp%nkbz,BSp%q(:,iq),opt_cvk(:,:,:,:,iq))
     !
     ! Fill ket0 using the same ordering for the indeces as the one used for the excitonic Hamiltonian.
     ! Note that only the resonant part is used here.
     do spin=1,nsppol
       spad=(spin-1)*BSp%nreh(1)
       do ik_bz=1,BSp%nkbz
         do iv=BSp%lomo_spin(spin),BSp%homo_spin(spin)
           do ic=BSp%lumo_spin(spin),BSp%nbnds
             trans_idx = BSp%vcks2t(iv,ic,ik_bz,spin)
             if (trans_idx>0) kets(trans_idx+spad,iq)=opt_cvk(ic,iv,ik_bz,spin,iq)
           end do
         end do
       end do
     end do
    end if

 end do

 call timab(692,2,tsec) ! exc_haydock_driver(prep)
 call timab(693,1,tsec) ! exc_haydock_driver(wo lf) - that is, without local field
 !
 ! ========================================================
 ! === Write the Optical Matrix Elements to NetCDF file ===
 ! ========================================================

 !if (.false.) then
 !  ome_fname='test_OME.nc'
 !  call exc_write_optme(ome_fname,minb,maxb,BSp%nkbz,Wfd%nsppol,BSp%nq,opt_cvk,ierr)
 !end if

 ! =======================================================
 ! === Make EPS RPA and GW without local-field effects ===
 ! =======================================================
 ABI_MALLOC(eps_rpanlf,(BSp%nomega,BSp%nq))
 ABI_MALLOC(dos_ks,(BSp%nomega))
 ABI_MALLOC(eps_gwnlf ,(BSp%nomega,BSp%nq))
 ABI_MALLOC(dos_gw,(BSp%nomega))
 
 if (BSp%use_interp) then
   call wrtout(std_out," Calculating Interpolated RPA NLF and QP NLF epsilon","COLL")

   call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,BSp%lomo_min,BSp%homo_spin,Kmesh_dense,KS_BSt_dense,BSp%nq,nsppol,&
&    opt_cvk,Cryst%ucvol,BSp%broad,BSp%nomega,BSp%omega,eps_rpanlf,dos_ks)

   call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,BSp%lomo_min,BSp%homo_spin,Kmesh_dense,QP_BSt_dense,BSp%nq,nsppol,&
&    opt_cvk,Cryst%ucvol,Bsp%broad,BSp%nomega,BSp%omega,eps_gwnlf,dos_gw)

 else
   call wrtout(std_out," Calculating RPA NLF and QP NLF epsilon","COLL")

   call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,BSp%lomo_min,BSp%homo_spin,Kmesh,KS_BSt,BSp%nq,nsppol,opt_cvk,&
&    Cryst%ucvol,BSp%broad,BSp%nomega,BSp%omega,eps_rpanlf,dos_ks)

   call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,BSp%lomo_min,BSp%homo_spin,Kmesh,QP_BSt,BSp%nq,nsppol,opt_cvk,&
&    Cryst%ucvol,Bsp%broad,BSp%nomega,BSp%omega,eps_gwnlf,dos_gw)
 end if

 if (my_rank==master) then ! Only master works.
   !
   ! Master node writes final results on file.
   call exc_write_data(BSp,BS_files,"RPA_NLF_MDF",eps_rpanlf,dos=dos_ks)

   call exc_write_data(BSp,BS_files,"GW_NLF_MDF",eps_gwnlf,dos=dos_gw)

   ! Computing and writing tensor in files

   ! RPA_NLF
   ABI_MALLOC(tensor_cart_rpanlf,(BSp%nomega,6))
   ABI_MALLOC(tensor_red_rpanlf,(BSp%nomega,6))

   call wrtout(std_out," Calculating RPA NLF dielectric tensor","COLL")
   call haydock_mdf_to_tensor(BSp,Cryst,eps_rpanlf,tensor_cart_rpanlf, tensor_red_rpanlf, ierr)

   if(ierr == 0) then
      ! Writing tensor
      call exc_write_tensor(BSp,BS_files,"RPA_NLF_TSR_CART",tensor_cart_rpanlf)
      call exc_write_tensor(BSp,BS_files,"RPA_NLF_TSR_RED",tensor_red_rpanlf)
   else 
      write(msg,'(3a)')&
&       'The RPA_NLF dielectric complex tensor cannot be computed',ch10,&
&       'There must be 6 different q-points in long wavelength limit (see gw_nqlwl)'
      MSG_COMMENT(msg)
   end if

   ABI_FREE(tensor_cart_rpanlf)
   ABI_FREE(tensor_red_rpanlf)

   ! GW_NLF
   ABI_MALLOC(tensor_cart_gwnlf,(BSp%nomega,6))
   ABI_MALLOC(tensor_red_gwnlf,(BSp%nomega,6))

   call wrtout(std_out," Calculating GW NLF dielectric tensor","COLL")

   call haydock_mdf_to_tensor(BSp,Cryst,eps_gwnlf,tensor_cart_gwnlf, tensor_red_gwnlf, ierr)

   if(ierr == 0) then
      ! Writing tensor
      call exc_write_tensor(BSp,BS_files,"GW_NLF_TSR_CART",tensor_cart_gwnlf)
      call exc_write_tensor(BSp,BS_files,"GW_NLF_TSR_RED",tensor_red_gwnlf)
   else
      write(msg,'(3a)')&
&       'The GW_NLF dielectric complex tensor cannot be computed',ch10,&
&       'There must be 6 different q-points in long wavelength limit (see gw_nqlwl)'
      MSG_COMMENT(msg)
   end if

   ABI_FREE(tensor_cart_gwnlf)
   ABI_FREE(tensor_red_gwnlf)
 
   !call wrtout(std_out," Checking Kramers Kronig on Excitonic Macroscopic Epsilon","COLL")
   !call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_exc(:,1))

   !call wrtout(std_out," Checking Kramers Kronig on RPA NLF Macroscopic Epsilon","COLL")
   !call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_rpanlf(:,1))

   !call wrtout(std_out," Checking Kramers Kronig on GW NLF Macroscopic Epsilon","COLL")
   !call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_gwnlf(:,1))

   !call wrtout(std_out," Checking f-sum rule on Excitonic Macroscopic Epsilon","COLL")

   !if (BSp%exchange_term>0) then 
   !  MSG_COMMENT(' f-sum rule should be checked without LF')
   !end if
   !call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_exc(:,1)),drude_plsmf)

   !call wrtout(std_out," Checking f-sum rule on RPA NLF Macroscopic Epsilon","COLL")
   !call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_rpanlf(:,1)),drude_plsmf)

   !call wrtout(std_out," Checking f-sum rule on GW NLF Macroscopic Epsilon","COLL")
   !call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_gwnlf(:,1)),drude_plsmf)
 end if ! my_rank==master

 ABI_FREE(opt_cvk)
 !call xmpi_barrier(comm)
 !
 ! The ket for the approximated DOS.
 if (prtdos) then 
   MSG_WARNING("Calculating DOS with Haydock method")
   ABI_CHECK(BSp%use_coupling==0,"DOS with coupling not coded")
   iq = BSp%nq + 1
   if (my_rank==master) then
     !call random_seed()
     do itt=1,SUM(Bsp%nreh)
       call RANDOM_NUMBER(rand_phi)
       rand_phi = two_pi*rand_phi
       kets(itt,iq) = CMPLX( COS(rand_phi), SIN(rand_phi) )
     end do
     ! Normalize the vector.
     !norm = SQRT( DOT_PRODUCT(kets(:,iq), kets(:,iq)) ) 
     !kets(:,iq) = kets(:,iq)/norm
   end if
   call xmpi_bcast(kets(:,iq),master,comm,ierr)
 end if

 call timab(693,2,tsec) ! exc_haydock_driver(wo lf    - that is, without local field
 call timab(694,1,tsec) ! exc_haydock_driver(apply

 ABI_MALLOC(green,(BSp%nomega,nkets))

 if (BSp%use_coupling==0) then 
   if (BSp%use_interp) then
      call haydock_herm_interp(BSp,BS_files,Cryst,Psps,Pawtab,Hdr_bse,hsize,hsize_dense,my_t1,my_t2,hreso,&
&       nkets,kets,grid,Wfd,Wfd_dense,Kmesh,Kmesh_dense,Vcp_dense,green,comm)
   else
      call haydock_herm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hreso,nkets,kets,green,comm)
   end if
 else
   if (BSp%use_interp) then
     MSG_ERROR("BSE Interpolation with coupling is not supported")
   else
     call haydock_psherm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hreso,hcoup,nkets,kets,green,comm)
   end if
 end if
 !
 ! Add 1 to have the real part right.
 green = one + green

 ABI_FREE(kets)

 call timab(694,2,tsec) ! exc_haydock_driver(apply
 call timab(695,1,tsec) ! exc_haydock_driver(end)

 if (my_rank==master) then ! Master writes the final results.
   !
   if (prtdos) then
     ABI_MALLOC(dos,(BSp%nomega))
     dos = -AIMAG(green(:,BSp%nq+1))
     call exc_write_data(BSp,BS_files,"EXC_MDF",green,dos=dos)
     ABI_FREE(dos)
   else 
     call exc_write_data(BSp,BS_files,"EXC_MDF",green)
   end if
   !
   ! =========================
   ! === Write out Epsilon ===
   ! =========================

   ABI_MALLOC(tensor_cart,(BSp%nomega,6))
   ABI_MALLOC(tensor_red,(BSp%nomega,6))

   call wrtout(std_out," Calculating EXC dielectric tensor","COLL")
   call haydock_mdf_to_tensor(BSp,Cryst,green,tensor_cart,tensor_red,ierr)

   if (ierr == 0) then
       ! Writing tensor
       call exc_write_tensor(BSp,BS_files,"EXC_TSR_CART",tensor_cart)
       call exc_write_tensor(BSp,BS_files,"EXC_TSR_RED",tensor_red)
   else
       write(msg,'(3a)')&
&        'The EXC dielectric complex tensor cannot be computed',ch10,&
&        'There must be 6 different q-points in long wavelength limit (see gw_nqlwl)'
       MSG_COMMENT(msg)
   end if

   ABI_FREE(tensor_cart)
   ABI_FREE(tensor_red)
   !
   ! This part will be removed when fldiff will be able to compare two mdf files.
   write(ab_out,*)" "
   write(ab_out,*)"Macroscopic dielectric function:"
   write(ab_out,*)"omega [eV] <KS_RPA_nlf>  <GW_RPA_nlf>  <BSE> "
   do io=1,MIN(BSp%nomega,10)
     omegaev = REAL(BSp%omega(io))*Ha_eV
     ks_avg  = SUM( eps_rpanlf(io,:)) / Bsp%nq
     gw_avg  = SUM( eps_gwnlf (io,:)) / Bsp%nq
     exc_avg = SUM( green     (io,:)) / BSp%nq
     write(ab_out,'(7f9.4)')omegaev,ks_avg,gw_avg,exc_avg
   end do
   write(ab_out,*)" "

   ! Write MDF file with the final results.
   ! FIXME: It won't work if prtdos == True
#ifdef HAVE_TRIO_ETSF_IO
     NCF_CHECK(ncfile_create(ncf,TRIM(BS_files%out_basename)//"_MDF.nc", NF90_CLOBBER), "Creating MDF file")
     call crystal_ncwrite(Cryst,ncf%ncid)
     call mdfs_ncwrite(ncf%ncid, Bsp, green, eps_rpanlf, eps_gwnlf)
     NCF_CHECK(ncfile_close(ncf),"Closing MDF file")
#else
     ABI_UNUSED(ncf%ncid)
#endif
 end if 

 ABI_FREE(green)
 ABI_FREE(eps_rpanlf)
 ABI_FREE(eps_gwnlf)
 ABI_FREE(dos_ks)
 ABI_FREE(dos_gw)

 ABI_FREE(hreso)
 if (allocated(hcoup)) then
   ABI_FREE(hcoup)
 end if

 call timab(695,2,tsec) ! exc_haydock_driver(end)
 call timab(690,2,tsec) ! exc_haydock_driver

end subroutine exc_haydock_driver
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_herm_interp
!! NAME
!! haydock_herm_interp
!!
!! FUNCTION
!!  Reads the excitonic Hamiltonian from file and construct the Lanczos set of vectors 
!!  by iterative matrix-vector multiplications.
!!
!! INPUTS
!! BSp<excparam>=Parameters for the Bethe-Salpeter calculation.
!! BS_files<excparam>=Files associated to the bethe_salpeter code.
!! Cryst<crystal_t>=Info on the crystalline structure.
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials.
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data.
!! hize=Size of the excitonic matrix.
!! my_t1,my_t2=First and last columns treated by this node.
!! hmat(hsize,my_t1:my_t2)=Excitonic matrix.
!! nkets=Number of starting vectors for Haydock method.
!! kets(hsize,nkets)=The kets in the eh representation.
!! comm=MPI communicator.
!!
!! OUTPUT
!!  green(BSp%nomega,nkets)=
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine haydock_herm_interp(BSp,BS_files,Cryst,Psps,Pawtab,Hdr_bse,hsize,hsize_dense,my_t1,my_t2,hmat,&
& nkets,kets,grid,Wfd,Wfd_dense,Kmesh_coarse,Kmesh_dense,Vcp_dense,green,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_herm_interp'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize,hsize_dense,my_t1,my_t2,nkets,comm
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh_coarse,Kmesh_dense
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(pseudopotential_type),intent(in) :: Psps
 type(Hdr_type),intent(in) :: Hdr_bse
 type(double_grid_t),intent(in) :: grid
 type(wfd_t),intent(inout) :: Wfd,Wfd_dense
!arrays
 complex(dp),intent(out) :: green(BSp%nomega,nkets)
 complex(dpc),intent(in) :: hmat(hsize,my_t1:my_t2),kets(hsize_dense,nkets)
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
!DBYG
 type(vcoul_t),intent(in) :: Vcp_dense

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0,spin1=1
 integer :: inn,itt,nproc,my_rank,ierr,iovlp,ix,iy,iz,ii ! itp,itp_
 integer :: niter_file,niter_max,niter_done,nsppol,iq,my_nt,term_type,n_all_omegas
 real(dp) :: norm,nfact
 logical :: can_restart,is_converged,is_resonant,diago_is_real,use_mpio
 complex(dpc) :: factor
 character(len=500) :: msg
 character(len=fnlen),parameter :: tag_file="_HAYDR_SAVE"
 character(len=fnlen) :: restart_file,out_file,tmpfname,hreso_fname
 type(haydock_type) :: haydock_file
 type(interpolator_t) :: interpolator
!arrays
 integer,allocatable :: div2kdense(:,:), kdense2div(:) 
 real(dp),allocatable :: bb_file(:)
 real(dp),allocatable :: bb(:)
 complex(dpc),allocatable :: aa(:),phi_nm1(:),phi_n(:),hphi_n(:),hphi_nm1(:)
 complex(dpc),allocatable :: work_coeffs(:,:)
 complex(dpc),allocatable :: aa_file(:),phi_n_file(:),phi_nm1_file(:)
 complex(dpc),allocatable :: ket0(:),all_omegas(:),green_temp(:,:)
 complex(dpc),allocatable :: all_hmat(:,:,:),all_acoeffs(:,:,:),all_bcoeffs(:,:,:),all_ccoeffs(:,:,:)
 complex(dpc),allocatable :: diag_coarse(:,:),diag_dense(:),hinterp(:,:)
 logical :: check(2)

!************************************************************************

 ABI_CHECK(Bsp%nsppol==1,"nsppol > 1 not implemented yet")

 nproc  = xcomm_size(comm); my_rank= xcomm_rank(comm)
 nsppol = Hdr_bse%nsppol

 if (BSp%use_interp) then
   MSG_COMMENT("No parallelization in Interpolation")
   my_nt = hsize_dense
 else
   my_nt = my_t2-my_t1+1
 end if

 ABI_CHECK(nproc == 1,"Parallelization not available in interpolation")
 ABI_CHECK(my_nt>0,"One of the processors has zero columns")


 call interpolator_init(interpolator, grid, Wfd_dense, Wfd, Kmesh_dense, Kmesh_coarse, BSp, Cryst, Psps, Pawtab, INTERP_YG)

 if(BSp%sum_overlaps) then
   call interpolator_normalize(interpolator)
 end if

 ! Table of correspondance between kdense and idiv needed for interpolation
 ! TODO : allocate in memory only if needed (method 1 & method 3)
 ABI_MALLOC(kdense2div,(grid%nbz_dense))
 ABI_MALLOC(div2kdense,(grid%nbz_coarse,grid%ndiv))
 call compute_corresp(grid,div2kdense,kdense2div)

 if (any(BSp%interp_mode == [2,3])) then
   ! Read a, b, c coefficient matrices from file.
   ! For the time being, we read the full matrix in a temporary array, and 
   ! then we store the data in a form suitable for the interpolation.
   is_resonant=.TRUE.; diago_is_real=(.not.BSp%have_complex_ene); use_mpio=.FALSE.

   if (BS_files%in_hreso /= BSE_NOFILE) then
     hreso_fname = BS_files%in_hreso
   else 
     hreso_fname = BS_files%out_hreso
   end if

   tmpfname = hreso_fname; ii = LEN_TRIM(hreso_fname)

   ! Allocate workspace array
   ABI_MALLOC(work_coeffs,(hsize,my_t1:my_t2))
   ABI_CHECK_ALLOC("out of memory in work_coeffs")

   ! TODO: Write new IO routines to read MPI-distributed data in a format suitable for the interpolation
   tmpfname(ii-2:ii+1) = 'ABSR'
   call exc_read_rcblock(tmpfname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,my_t1,my_t2,work_coeffs,use_mpio,comm)

   ABI_MALLOC(all_acoeffs,(hsize,hsize,interpolator%nvert))
   ABI_CHECK_ALLOC("out of memory in all_acoeffs")
   do iovlp=1,interpolator%nvert
     all_acoeffs(:,:,iovlp) = work_coeffs(:,interpolator%corresp(:,iovlp,spin1))
   end do

   tmpfname(ii-2:ii+1) = 'BBSR'
   call exc_read_rcblock(tmpfname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,my_t1,my_t2,work_coeffs,use_mpio,comm)

   ABI_MALLOC(all_bcoeffs,(hsize,hsize,interpolator%nvert))
   ABI_CHECK_ALLOC("out of memory in all_bcoeffs")
   do iovlp=1,interpolator%nvert
     all_bcoeffs(:,:,iovlp) = work_coeffs(:,interpolator%corresp(:,iovlp,spin1))
   end do

   tmpfname(ii-2:ii+1) = 'CBSR'
   call exc_read_rcblock(tmpfname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,my_t1,my_t2,work_coeffs,use_mpio,comm)

   ABI_MALLOC(all_ccoeffs,(hsize,hsize,interpolator%nvert))
   ABI_CHECK_ALLOC("out of memory in all_ccoeffs")
   do iovlp=1,interpolator%nvert
     all_ccoeffs(:,:,iovlp) = work_coeffs(:,interpolator%corresp(:,iovlp,spin1))
   end do

   ABI_FREE(work_coeffs)
 end if 

 ! TODO: Use vector
 !call wrtout(std_out," Computing diagonal term on the coarse mesh","COLL")
 ABI_MALLOC(diag_coarse,(hsize,hsize))
 ABI_CHECK_ALLOC("out of memory in diag_coarse")
 diag_coarse = czero
 do itt=1,BSp%nreh(spin1) ! 1 is for spin 1
   diag_coarse(itt,itt) = Bsp%Trans(itt,spin1)%en
 end do
 !
 ! Compute overlaps & compute all hmat
 ABI_MALLOC(all_hmat,(hsize,hsize,interpolator%nvert))
 ABI_CHECK_ALLOC("out of memory in all_hmat")

 do iovlp = 1,interpolator%nvert
   ix = (iovlp-1)/4
   iy = (iovlp-ix*4-1)/2
   iz = (iovlp-ix*4-iy*2-1)

   all_hmat(:,:,iovlp) = hmat(:,interpolator%corresp(:,iovlp,spin1)) - diag_coarse(:,interpolator%corresp(:,iovlp,spin1))

   !do itp=1,hsize
   !  itp_interp = corresp(itp,iovlp)
   !  do itt=1,hsize
   !    !all_hmat(itt,itp,iovlp) = hmat(itt,itp_interp) - diag_coarse(itt,itp_interp)  
   !    all_hmat(itt,itp,iovlp) = hmat(itt,itp_interp) 
   !    if (itt == itp_interp) all_hmat(itt,itp,iovlp) = all_hmat(itt,itp,iovlp) - Bsp%Trans(itt,spin1)%en
   !  end do
   !end do
 end do

 ABI_FREE(diag_coarse)

 if (any(BSp%interp_mode == [2,3])) then 
   ! If mode = divergence abc, not yet interpolated product ...
   write(msg,"(a,f8.1,a)")"Memory needed for hinterp = ",one*(hsize_dense**2)*2*dpc*b2Mb," Mb"
   call wrtout(std_out,msg,"COLL")

   ABI_MALLOC(hinterp,(hsize_dense,hsize_dense))
   ABI_CHECK_ALLOC('Out of memory in hinterp')

   call compute_hinterp(BSp,hsize,hsize_dense,all_hmat,grid,&
&    BSp%maxnbndv*BSp%maxnbndc,interpolator,kdense2div,&
&    all_acoeffs,all_bcoeffs,all_ccoeffs,Kmesh_dense,Vcp_dense,cryst%gmet,hinterp)

   ABI_FREE(all_acoeffs)
   ABI_FREE(all_bcoeffs)
   ABI_FREE(all_ccoeffs)
 end if

 write(msg,'(a,i0)')' Haydock algorithm with MAX number of iterations: ',BSp%niter
 call wrtout(std_out,msg,"COLL")
 !
 ! Select the terminator for the continued fraction.
 term_type=0; if (Bsp%hayd_term>0) term_type=1
 call wrtout(std_out,sjoin("Using terminator type: ",itoa(term_type)),"COLL")
 !
 ! Check for presence of the restart file.
 can_restart=.FALSE.
 if ( BS_files%in_haydock_basename /= BSE_NOFILE) then
   restart_file = TRIM(BS_files%in_haydock_basename)//TRIM(tag_file)
   if (file_exists(restart_file) ) then
     can_restart=.TRUE.
     msg = " Restarting Haydock calculation from file: "//TRIM(restart_file)
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")
   else
     can_restart=.FALSE.
     call wrtout(ab_out," WARNING: cannot find restart file: "//TRIM(restart_file),"COLL")
   end if
 end if
 ABI_CHECK(.not.can_restart,"restart not yet implemented")
 
 ! Open the file and write basic dimensions and info.
 if (my_rank==master) then
   out_file = TRIM(BS_files%out_basename)//TRIM(tag_file)
   call open_haydock(out_file,haydock_file) 
   haydock_file%hsize = hsize_dense
   haydock_file%use_coupling = Bsp%use_coupling
   haydock_file%op = BSE_HAYD_IMEPS
   haydock_file%nq = nkets
   haydock_file%broad = Bsp%broad
   call write_dim_haydock(haydock_file)
 end if

 ! Compute diagonal part of the dense Ham
 ABI_MALLOC(diag_dense,(hsize_dense))
 do itt=1,BSp%nreh_interp(spin1) ! 1 is for spin 1
   diag_dense(itt) = Bsp%Trans_interp(itt,spin1)%en
 end do
 !
 ! Calculate green(w) for the different starting points.
 green=czero
 do iq=1,nkets
   ABI_MALLOC(ket0,(hsize_dense))
   ket0=kets(:,iq)
   !
   niter_file=0
   if (can_restart) then
     call haydock_restart(BSp,restart_file,BSE_HAYD_IMEPS,iq,hsize_dense,&
&      niter_file,aa_file,bb_file,phi_nm1_file,phi_n_file,comm)
   end if 
   !
   ! For n>1, we have:
   !  1) a_n = <n|H|n>
   !  2) b_n = || H|n> - a_n|n> -b_{n-1}|n-1> ||
   !  3) |n+1> = [H|n> -a_n|n> -b_{n-1}|n-1>]/b_n
   !
   ! The sequences starts with |1> normalized to 1 and b_0 =0, therefore:
   !  a_1 = <1|H|1>
   !  b_1 = || H|1> - a_1|1> ||
   !  |2> = [H|1> - a_1|1>]/b_1
   !
   ABI_MALLOC(hphi_n,(hsize_dense))
   ABI_MALLOC(hphi_nm1,(hsize_dense))
   ABI_MALLOC(phi_nm1,(hsize_dense))
   ABI_MALLOC(phi_n,(hsize_dense))

   niter_max = niter_file + Bsp%niter
   ABI_MALLOC(aa,(niter_max))
   ABI_MALLOC(bb,(niter_max))
   aa=czero; bb=zero

   if (niter_file==0) then       ! Calculation from scratch.
     phi_nm1=ket0(:)   ! Select the slice treated by this node.
     norm = DZNRM2(hsize_dense,ket0,1) ! Normalization  
     phi_nm1=phi_nm1/norm

     ! hphi_n = MATMUL(hmat,phi_nm1)
     if (any(BSp%interp_mode == [2,3])) then
       call haydock_interp_matmul(BSp,hsize,hsize_dense,all_hmat,diag_dense,phi_nm1,hphi_n,grid,&
&        BSp%maxnbndv*BSp%maxnbndc,interpolator,div2kdense,kdense2div,hinterp)
     else
       call haydock_interp_matmul(BSp,hsize,hsize_dense,all_hmat,diag_dense,phi_nm1,hphi_n,grid,&
&        BSp%maxnbndv*BSp%maxnbndc,interpolator,div2kdense,kdense2div)
     end if
     !temp_phi = hphi_n

     !call xgemv('N',hsize_dense,my_nt,cone,diag_dense,hsize_dense,phi_nm1,1,czero,hphi_n,1)
     !call xmpi_sum(hphi_n,comm,ierr)

     aa(1)=xdotc(my_nt,phi_nm1,1,hphi_n(:),1)
     call xmpi_sum(aa(1:1),comm,ierr)

     phi_n = hphi_n(:) - aa(1)*phi_nm1

     bb(1) = xdotc(my_nt,phi_n,1,phi_n,1)
     call xmpi_sum(bb(1:1),comm,ierr)
     bb(1) = SQRT(bb(1))

     phi_n = phi_n/bb(1)
     niter_done=1

   else ! Use the previous a and b.
     niter_done=niter_file
     aa(1:niter_done) = aa_file
     bb(1:niter_done) = bb_file
     phi_nm1=phi_nm1_file(:)   ! Select the slice treated by this node.
     phi_n  =phi_n_file  (:)
   end if

   if (can_restart) then
     ABI_FREE(aa_file)
     ABI_FREE(bb_file)
     ABI_FREE(phi_nm1_file)
     ABI_FREE(phi_n_file)
   end if

   ! Multiplicative factor (k-point sampling and unit cell volume)  
   ! TODO be careful with the spin here
   ! TODO four_pi comes from the coulomb term 1/|q| is already included in the 
   ! oscillators hence the present approach wont work if a cutoff interaction is used.
   nfact = -four_pi/(Cryst%ucvol*BSp%nkbz_interp)
   if (nsppol==1) nfact=two*nfact

   factor = nfact*(DZNRM2(hsize_dense,ket0,1)**2)

   ! Which quantity should be checked for convergence?
   check = (/.TRUE.,.TRUE./)
   if (ABS(Bsp%haydock_tol(2)-one)<tol6) check = (/.TRUE. ,.FALSE./)
   if (ABS(Bsp%haydock_tol(2)-two)<tol6) check = (/.FALSE.,.TRUE./)

   ! Create new frequencies "mirror" in negative range to add 
   ! their contributions. Can be improved by computing only once
   ! zero frequency, but loosing clearness
   n_all_omegas = 2*BSp%nomega

   ABI_MALLOC(all_omegas,(n_all_omegas))
   ! Put all omegas with frequency > 0 in table
   all_omegas(BSp%nomega+1:n_all_omegas) = BSp%omega
   ! Put all omegas with frequency < 0
   ! Warning, the broadening must be kept positive
   all_omegas(1:BSp%nomega) = -DBLE(BSp%omega(BSp%nomega:1:-1)) + j_dpc*AIMAG(BSp%omega(BSp%nomega:1:-1))

   ABI_MALLOC(green_temp,(n_all_omegas,nkets))

   if (any(BSp%interp_mode == [2,3])) then
     call haydock_herm_algo_interp(BSp,niter_done,niter_max,n_all_omegas,all_omegas,BSp%haydock_tol(1),check,hsize,hsize_dense,&
&      my_t1,my_t2,all_hmat,diag_dense,grid,factor,term_type,aa,bb,phi_nm1,phi_n,&
&      green_temp(:,iq),inn,is_converged,BSp%maxnbndv*BSp%maxnbndc,interpolator,div2kdense,kdense2div,&
&      comm,hinterp)
   else
     call haydock_herm_algo_interp(BSp,niter_done,niter_max,n_all_omegas,all_omegas,BSp%haydock_tol(1),check,hsize,hsize_dense,&
&      my_t1,my_t2,all_hmat,diag_dense,grid,factor,term_type,aa,bb,phi_nm1,phi_n,&
&      green_temp(:,iq),inn,is_converged,BSp%maxnbndv*BSp%maxnbndc,interpolator,div2kdense,kdense2div,&
&      comm)
   end if

   ! Computing result from two ranges of frequencies
   ! The real part is added, the imaginary part is substracted
   green(:,iq) = green_temp(BSp%nomega+1:n_all_omegas,iq)+CONJG(green_temp(BSp%nomega:1:-1,iq))

   ABI_FREE(all_omegas)
   ABI_FREE(green_temp)
   !
   ! Save the a"s and the b"s for possible restarting.
   ! 1) Info on the Q.
   ! 2) Number of iterations performed.
   ! 3) do iter=1,niter_performed 
   !      aa(iter),bb(iter)
   !    end do
   ! 4) |n-1>
   !    |n>
   !
   hphi_nm1 = czero
   hphi_nm1(1:hsize_dense) = phi_nm1
   call xmpi_sum_master(hphi_nm1,master,comm,ierr)

   hphi_n = czero
   hphi_n(1:hsize_dense) = phi_n
   call xmpi_sum_master(hphi_n,master,comm,ierr)

   if (my_rank==master) then 
     ! Write data for restarting
     call write_haydock(haydock_file, hsize_dense, Bsp%q(:,iq), aa, bb, hphi_n, hphi_nm1, MIN(inn,niter_max), factor)
   end if

   ABI_FREE(hphi_n)
   ABI_FREE(hphi_nm1)
   ABI_FREE(phi_nm1)
   ABI_FREE(phi_n)
   ABI_FREE(aa)
   ABI_FREE(bb)
   ABI_FREE(ket0)
 end do ! iq

 ABI_FREE(diag_dense)
 ABI_FREE(all_hmat)
 ABI_FREE(div2kdense)
 ABI_FREE(kdense2div)

 call interpolator_free(interpolator)

 if (allocated(hinterp)) then
   ABI_FREE(hinterp)
 end if

 if (my_rank==master) call close_haydock(haydock_file)

 call xmpi_barrier(comm)

end subroutine haydock_herm_interp
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_herm_algo_interp
!! NAME
!! haydock_herm_algo_interp
!!
!! FUNCTION
!!
!! INPUTS
!!  niter_done=Number of iterations already performed (0 if the run starts from scratch).
!!  niter_max=Max number of iterations. Always > niter_done
!!  nomega=Number of Frequency points for the evaluation of the matrix element.
!!  omega(nomega)=Frequency set (imaginary part is already included).
!!  tol_iter=Tolerance used to stop the algorithm.
!!  check(2)=Logical flags to specify where both the real and the imaginary part of the 
!!    matrix elements of the Green functions have to be checked for convergence. 
!!  hsize=Size of the blocks.
!!  my_t1,my_t2=Indices of the first and last column stored treated by this done.
!!  term_type=0 if no terminator is used, 1 otherwise.
!!  hmat(hsize,my_t1:my_t2)=The columns of the block.
!!  factor
!!  ntrans = Number of transitions
!!  nbnd_coarse = Product of number of conduction and number of valences
!!  corresp = mapping between coarse points and neighbours
!!  overlaps = overlaps of wavefunctions between dense k-point coarse neighbours and bands
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  green(nomega)=Output matrix elements.
!!  inn=Last iteration performed.
!!  is_converged=.TRUE. of the algorithm converged.
!!
!! SIDE EFFECTS
!!  phi_nm1(my_t2-my_t1+1), phi_n(my_t2-my_t1+1)
!!    input: vectors used to initialize the iteration
!!    output: the vectors obtained in the last iteration 
!!  aa(niter_max) and bb(niter_max)
!!    if niter_done>0: aa(1:niter_done), bb(1:niter_done) store the coefficients of the previous run.
!!    when the routine returns aa(1:inn) and bb(1:inn) contain the matrix elements of the tridiagonal form.
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine haydock_herm_algo_interp(BSp,niter_done,niter_max,nomega,omega,tol_iter,check,hsize,hsize_dense,&
& my_t1,my_t2,hmat,diag_dense,grid,factor,term_type,aa,bb,phi_nm1,phi_n,&
& green,inn,is_converged,nbnd_coarse,interpolator,div2kdense,kdense2div,comm,hinterp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_herm_algo_interp'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter_max,niter_done,nomega
 integer,intent(in) :: hsize,my_t1,my_t2,term_type
 integer,intent(in) :: hsize_dense,nbnd_coarse,comm
 integer,intent(out) :: inn
 logical,intent(out) :: is_converged
 real(dp),intent(in) :: tol_iter
 complex(dpc),intent(in) :: factor
 type(excparam),intent(in) :: BSp
 type(double_grid_t),intent(in) :: grid
 type(interpolator_t),intent(in) :: interpolator
!arrays
 integer,intent(in) :: div2kdense(grid%nbz_coarse,grid%ndiv), kdense2div(grid%nbz_dense) 
 real(dp),intent(inout) :: bb(niter_max)
 complex(dpc),intent(out) :: green(nomega)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(inout) :: aa(niter_max)
 complex(dpc),intent(in) :: hmat(hsize,hsize,8)
 complex(dpc),intent(in) :: diag_dense(hsize_dense)
 complex(dpc),intent(inout) :: phi_nm1(hsize_dense)
 complex(dpc),intent(inout) :: phi_n  (hsize_dense)
 complex(dpc),optional,intent(in) :: hinterp(hsize_dense,hsize_dense)
 logical,intent(in) :: check(2)

!Local variables ------------------------------
!scalars
 integer :: ierr,my_nt,niter_min,nconv
 character(len=500) :: msg
 logical,parameter :: force_real=.TRUE.
!arrays
 real(dp) :: abs_err(nomega,2) !,rel_err(nomega,2)
 complex(dpc),allocatable :: oldg(:),newg(:)
 complex(dpc),allocatable :: phi_np1(:),hphi_n(:),cfact(:)
 logical :: test(2)

!************************************************************************

 ! The sequences starts with |1> normalized to 1 and b_0 =0, therefore:
 !  a_1 = <1|H|1>
 !  b_1 = || H|1> - a_1|1> ||
 !  |2> = [H|1> - a_1|1>]/b_1
 !
 ! For n>1 we have
 !  1) a_n = <n|H|n>
 !  2) b_n = || H|n> - a_n|n> -b_{n-1}|n-1> ||
 !  3) |n+1> = [H|n> -a_n|n> -b_{n-1}|n-1>]/b_n
 !
 my_nt = my_t2-my_t1+1
 my_nt = hsize_dense

 ABI_MALLOC(hphi_n,(hsize_dense))
 ABI_CHECK_ALLOC("out-of-memory hphi_n")

 ABI_MALLOC(phi_np1,(my_nt))

 ABI_MALLOC(oldg,(nomega))
 oldg=czero
 ABI_MALLOC(newg,(nomega))
 newg=czero
 ABI_MALLOC(cfact,(nomega))
 cfact=czero

 nconv=0
 do inn=niter_done+1,niter_max
   !
   ! hphi_n = MATMUL(hmat,phi_n)
   if (any(Bsp%interp_mode == [2,3])) then
     call haydock_interp_matmul(BSp,hsize,hsize_dense,hmat,diag_dense,phi_n,hphi_n,grid,&
&      nbnd_coarse,interpolator,div2kdense,kdense2div,hinterp)
   else
     call haydock_interp_matmul(BSp,hsize,hsize_dense,hmat,diag_dense,phi_n,hphi_n,grid,&
&      nbnd_coarse,interpolator,div2kdense,kdense2div)
   end if

   !call xgemv('N',hsize_dense,my_nt,cone,diag_dense,hsize_dense,phi_n,1,czero,hphi_n,1)
   !call xmpi_sum(hphi_n,comm,ierr)

   aa(inn) = xdotc(my_nt,phi_n,1,hphi_n(:),1)
   call xmpi_sum(aa(inn:inn),comm,ierr)
   if (force_real) aa(inn) = DBLE(aa(inn)) ! Matrix is Hermitian.

   ! |n+1> = H|n> - A(n)|n> - B(n-1)|n-1>
   phi_np1 = hphi_n(:) - aa(inn)*phi_n - bb(inn-1)*phi_nm1

   bb(inn) = xdotc(my_nt,phi_np1,1,phi_np1,1)
   call xmpi_sum(bb(inn),comm,ierr)
   bb(inn) = SQRT(bb(inn))

   phi_np1 = phi_np1/bb(inn)

   phi_nm1 = phi_n
   phi_n   = phi_np1

   write(msg,'(a,i0,a,3es12.4)')' Iteration number ',inn,', b_i RE(a_i) IM(a_i) ',bb(inn),REAL(aa(inn)),AIMAG(aa(inn))
   call wrtout(std_out,msg,"COLL")

   call continued_fract(inn,term_type,aa,bb,nomega,omega,cfact)

   newg= factor*cfact
   !
   ! Avoid spurious convergence.
   niter_min=4; if (niter_done>1) niter_min=niter_done+1
   if (inn>niter_min) then
     test=.TRUE.
     abs_err(:,1) = ABS(DBLE (newg-oldg))
     abs_err(:,2) = ABS(AIMAG(newg-oldg))
     !
     if (tol_iter>zero) then
       ! Test on the L1 norm.
       if (check(1)) test(1) = SUM(abs_err(:,1)) < tol_iter*SUM(ABS(DBLE (newg)))
       if (check(2)) test(2) = SUM(abs_err(:,2)) < tol_iter*SUM(ABS(AIMAG(newg)))

     else
       ! Stringent test for each point.
       if (check(1)) test(1) = ALL( abs_err(:,1) < -tol_iter*ABS(DBLE (newg)))
       if (check(2)) test(2) = ALL( abs_err(:,2) < -tol_iter*ABS(AIMAG(newg)))
     end if
     !
     if (ALL(test)) then
       nconv = nconv+1
     else
       nconv = 0
     end if
     if (nconv==2) then
       write(msg,'(a,es10.2,a,i0,a)')&
&        " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after ",inn," iterations."
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       EXIT
     end if
   end if

   oldg = newg
 end do ! inn

 green = newg
 if (nconv/=2) then
   write(msg,'(a,es10.2,a,i0,a)')&
&    " WARNING: Haydock algorithm did not converge within ",tol_iter," after ",niter_max," iterations."
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 is_converged = (nconv==2)

 ABI_FREE(oldg)
 ABI_FREE(newg)
 ABI_FREE(cfact)
 ABI_FREE(hphi_n)
 ABI_FREE(phi_np1)

end subroutine haydock_herm_algo_interp
!!***

!-------------------------------------------------------------------------------------

!!****f* m_haydock/compute_subhinterp
!! NAME
!! compute_subhinterp
!!
!! FUNCTION
!! TODO
!!
!! INPUTS
!! TODO
!!
!! OUTPUT
!! TODO
!!
!! PARENTS
!!      haydock
!!
!! SOURCE

subroutine compute_subhinterp(BSp,grid,nbnd_coarse,&
&  interpolator,kdense2div,&
&  work_coeffs,Cmat,ikp_dense,&
&  overlaps)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compute_subhinterp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnd_coarse
 integer,intent(in) :: ikp_dense
 type(excparam),intent(in) :: BSp
 type(double_grid_t),intent(in) :: grid
 type(interpolator_t),intent(in) :: interpolator
!arrays
 integer,intent(in) :: kdense2div(grid%nbz_dense)
 complex(gwpc),intent(in) :: overlaps(interpolator%mband_coarse,interpolator%mband_dense,interpolator%nvert)
 complex(dpc),intent(in) :: work_coeffs(nbnd_coarse,interpolator%nvert)
 complex(dpc),intent(out) :: Cmat(nbnd_coarse)

!Local variables ------------------------------
!scalars
 integer,parameter :: spin1=1,spin2=1
 integer :: iv1,ic1
 integer :: icp,ivp,idivp,ibndp_coarse,ibndp_coarse1,ineighbourp
 integer :: indwithnb
 integer :: lumo2,lomo2,humo2,homo2
 complex(dpc) :: tmp_val, tmp2, tmp4
!arrays
 complex(dpc),allocatable :: btemp(:),ctemp(:)

!*********************************************************************

 ABI_MALLOC(btemp,(interpolator%nvert*nbnd_coarse))
 btemp = czero

 ABI_MALLOC(ctemp,(interpolator%nvert*nbnd_coarse))
 ctemp = czero

 lumo2 = BSp%lumo_spin(spin2)
 lomo2 = BSp%lomo_spin(spin2)
 humo2 = BSp%humo_spin(spin2)
 homo2 = BSp%homo_spin(spin2)

 Cmat = czero
     
 idivp = kdense2div(ikp_dense)

 do ineighbourp = 1,interpolator%nvert

   btemp(((ineighbourp-1)*nbnd_coarse+1):(ineighbourp*nbnd_coarse)) = &
&       interpolator%interp_factors(ineighbourp,idivp)*work_coeffs(:,ineighbourp)

 end do !ineighbourp

 ! Loop over the (c', v') part of the right transition
 do ivp = lomo2,homo2
   do icp = lumo2,humo2

     ibndp_coarse = (ivp-lomo2)*BSp%maxnbndc+(icp-lumo2+1)
     ! Now we now it_dense, and itp_dense

     do ineighbourp = 1,interpolator%nvert

       do iv1 = lomo2, homo2
         ! BSp%lumo_spin(spin2),BSp%humo_spin(spin2)

         tmp4 = overlaps(iv1,ivp,ineighbourp)

         do ic1 = lumo2,humo2
           tmp2 = GWPC_CONJG(overlaps(ic1,icp,ineighbourp))
           ! BSp%lomo_spin(spin2),BSp%homo_spin(spin2)

           ibndp_coarse1 = (iv1-lomo2)*BSp%maxnbndc+(ic1-lumo2+1)
           indwithnb = (ineighbourp-1)*nbnd_coarse+ibndp_coarse1

           ctemp(indwithnb) = &
&             tmp4 &
&            *tmp2

         end do ! iv1
       end do !ic1

     end do !ineighbourp

     tmp_val = xdotc(interpolator%nvert*nbnd_coarse,ctemp,1,btemp,1)
     !tmp_val = DOT_PRODUCT(ctemp,btemp)

     Cmat(ibndp_coarse) = tmp_val
   end do !ivp
 end do !icp

 ABI_FREE(btemp)
 ABI_FREE(ctemp)

end subroutine compute_subhinterp
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/compute_hinterp
!! NAME
!! compute_hinterp
!!
!! FUNCTION
!! Compute interpolated matrix elements for methods 2 and 3
!!
!! INPUTS
!! BSp<type(excparam)=The parameter for the Bethe-Salpeter run.
!!  hsize_coarse=Size of the coarse Hamiltonian
!!  hsize_dense=Size of the dense Hamiltonian
!!  hmat(hsize_coarse,hsize_coarse,8)=Excitonic matrix 
!!  grid <double_grid_t> = Correspondence between coarse and dense k-mesh.
!!  ntrans
!!  nbnd_coarse
!!  corresp(hsize_coarse,8)=mapping between coarse points and neighbours
!!  overlaps(ntrans,nbnd_coarse,8)=overlaps of wavefunctions between dense k-point coarse neighbours and bands
!!  interp_factors(BSp%nreh(1),8,grid%ndiv) = Interpolation factors
!!  indices(BSp%nreh(1),grid%ndiv)
!!  acoeffs(hsize_coarse,hsize_coarse,8)
!!  bcoeffs(hsize_coarse,hsize_coarse,8)
!!  ccoeffs(hsize_coarse,hsize_coarse,8)
!! Kmesh_dense<type(kmesh_t)>=The list of k-points in the BZ, IBZ and symmetry tables.
!! Vcp<vcoul_t>=Coulomb interation in G-space on the dense Q-mesh
!! gmet(3,3)=Metric tensor in G-space
!!
!! OUTPUT
!!   hinterp(hsize_dense,hsize_dense)
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine compute_hinterp(BSp,hsize_coarse,hsize_dense,hmat,grid,nbnd_coarse,&
&  interpolator,kdense2div,&
&  acoeffs,bcoeffs,ccoeffs,Kmesh_dense,Vcp_dense,gmet,hinterp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compute_hinterp'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize_coarse,hsize_dense,nbnd_coarse !,ntrans
 type(excparam),intent(in) :: BSp
 type(double_grid_t),intent(in) :: grid
 type(vcoul_t),intent(in) :: Vcp_dense
 type(kmesh_t),intent(in) :: Kmesh_dense
 type(interpolator_t),intent(in) :: interpolator
!arrays
 integer,intent(in) :: kdense2div(grid%nbz_dense)
 real(dp),intent(in) :: gmet(3,3)
 complex(dpc),intent(in) :: hmat(hsize_coarse,hsize_coarse,interpolator%nvert)
 complex(dpc),intent(in) :: acoeffs(hsize_coarse,hsize_coarse,interpolator%nvert)
 complex(dpc),intent(in) :: bcoeffs(hsize_coarse,hsize_coarse,interpolator%nvert)
 complex(dpc),intent(in) :: ccoeffs(hsize_coarse,hsize_coarse,interpolator%nvert)
 complex(dpc),intent(out) :: hinterp(hsize_dense,hsize_dense)

!Local variables ------------------------------
!scalars
 integer,parameter :: spin1=1, spin2=1
 integer :: ic,iv,iv1,ic1,ik_dense,ik_coarse,it_coarse,it_dense,idiv,ibnd_coarse,ibnd_coarse1,ineighbour
 integer :: icp,ivp,ikp_dense,ikp_coarse,itp_coarse,itp_dense,idivp,ibndp_coarse,ibndp_coarse1,ineighbourp,itp_coarse1
 integer :: dump_unt,itc,it_dense1,indwithnb
 real(dp) :: factor,vc_sqrt_qbz,qnorm
 complex(dpc) :: term,http
 logical :: writeh
 character(len=500) :: msg
 logical :: newway
!arrays
 real(dp) :: kmkp(3),q2(3),shift(3),qinred(3),tsec(2)
 complex(dpc),allocatable :: btemp(:),ctemp(:),Cmat(:,:,:) !Temp matrices for optimized version
 complex(dpc),allocatable :: tmp_Cmat(:)
 complex(dpc),allocatable :: work_coeffs(:,:)
 integer,allocatable :: band2it(:)

!************************************************************************

 !TODO: This is most CPU-expensive part!!!!!

 newway = .True.

 if (any(BSp%interp_mode == [2,3])) then
   if(Vcp_dense%mode /= 'CRYSTAL' .and. Vcp_dense%mode /= 'AUXILIARY_FUNCTION') then
     MSG_BUG('Vcp_dense%mode not implemented yet !')
   end if
 end if
                                                                                                 
 if(BSp%nsppol > 1) then
   MSG_BUG("nsppol > 1 not yet implemented")
 end if

 call timab(696,1,tsec)
 factor = one/grid%ndiv

 hinterp = czero; term = czero

 ABI_MALLOC(btemp,(interpolator%nvert*nbnd_coarse))
 btemp = czero

 ABI_MALLOC(ctemp,(interpolator%nvert*nbnd_coarse))
 ctemp = czero

 ABI_MALLOC(Cmat,(nbnd_coarse,nbnd_coarse,interpolator%nvert))
 ABI_CHECK_ALLOC("out of memory in Cmat")
 Cmat = czero

 ABI_MALLOC(band2it,(nbnd_coarse))

 if(newway) then
   ABI_MALLOC(tmp_Cmat,(nbnd_coarse))
   ABI_CHECK_ALLOC("out of memory in tmp_Cmat")
   ABI_MALLOC(work_coeffs,(nbnd_coarse,interpolator%nvert))
   ABI_CHECK_ALLOC("out of memory in work_coeffs")
 end if

 do ik_dense = 1,grid%nbz_dense
   write(std_out,*) "Kdense = ",ik_dense,"/",grid%nbz_dense
   ik_coarse = grid%dense_to_coarse(ik_dense)
   do ikp_dense = 1,grid%nbz_dense
     ikp_coarse = grid%dense_to_coarse(ikp_dense)

     do iv1 = BSp%lomo_spin(spin2),BSp%homo_spin(spin2)
       do ic1 = BSp%lumo_spin(spin2),BSp%humo_spin(spin2)
         itp_coarse1 = BSp%vcks2t(iv1,ic1,ikp_coarse,spin2)
         ibndp_coarse1 = (iv1-BSp%lomo_spin(spin2))*BSp%maxnbndc+(ic1-BSp%lumo_spin(spin2)+1)

         band2it(ibndp_coarse1) = itp_coarse1 
       end do
     end do


     if (any(BSp%interp_mode == [2,3])) then
       ! Check if we are along the diagonal
       kmkp = Kmesh_dense%bz(:,ik_dense) - Kmesh_dense%bz(:,ikp_dense)

       call wrap2_pmhalf(kmkp(:),q2(:),shift(:))
       qinred = MATMUL(grid%kptrlatt_coarse,q2)

       ! We are outside the diagonal
       if (BSp%interp_mode==3 .and. ANY((ABS(qinred)-tol7) > one)) cycle

       qnorm = two_pi*SQRT(DOT_PRODUCT(q2,MATMUL(gmet,q2)))

       if(ALL(ABS(q2(:)) < 1.e-3)) then
         vc_sqrt_qbz = SQRT(Vcp_dense%i_sz)
       else
         vc_sqrt_qbz = SQRT(four_pi/qnorm**2)
       end if

       !!DEBUG CHK !
       !!COMPUTE Qpoint
       !call findqg0(iq_bz,g0,kmkp,Qmesh_dense%nbz,Qmesh_dense%bz,BSp%mG0)

       !! * Get iq_ibz, and symmetries from iq_bz
       !call get_BZ_item(Qmesh_dense,iq_bz,qbz,iq_ibz,isym_q,itim_q)

       !if(iq_ibz > 1 .and. ABS(vc_sqrt_qbz - Vcp_dense%vc_sqrt(1,iq_ibz)) > 1.e-3) then
       !   write(*,*) "vc_sqrt_qbz = ",vc_sqrt_qbz
       !   write(*,*) "Vcp_dense%vc_sqrt(1,iq_ibz) = ",Vcp_dense%vc_sqrt(1,iq_ibz)
       !   MSG_ERROR("vcp are not the same !")
       !else if(iq_ibz == 1 .and. ABS(vc_sqrt_qbz - SQRT(Vcp_dense%i_sz)) > 1.e-3) then
       !   write(*,*) "vc_sqrt_qbz = ",vc_sqrt_qbz
       !   write(*,*) "SQRT(Vcp_dense%i_sz) = ",SQRT(Vcp_dense%i_sz)
       !   MSG_ERROR("vcp are not the same !")
       !end if
       !!END DEBUG CHK !
     end if

     if(newway) then

       Cmat = czero

       work_coeffs = czero

       do ineighbour = 1,interpolator%nvert

         ! Loop over the (c, v) part of the left transition
         do iv = BSp%lomo_spin(spin1),BSp%homo_spin(spin1)
           do ic = BSp%lumo_spin(spin1),BSp%humo_spin(spin1)
        
             it_dense = BSp%vcks2t_interp(iv,ic,ik_dense,spin1)
             it_coarse = BSp%vcks2t(iv,ic,ik_coarse,spin1)
             ibnd_coarse = (iv-BSp%lomo_spin(spin1))*BSp%maxnbndc+(ic-BSp%lumo_spin(spin1)+1)

             itc = interpolator%corresp(it_coarse,ineighbour,spin1) 

             if (any(BSp%interp_mode == [1,3])) then
               work_coeffs(:,:) = hmat(itc,band2it(:),:)

               call compute_subhinterp(BSp,grid,nbnd_coarse,&
&        interpolator,kdense2div,&
&        work_coeffs,tmp_Cmat,ikp_dense,&
&        interpolator%overlaps(:,:,:,ikp_dense,spin2))

               if(BSp%interp_mode == 1) then
                 Cmat(ibnd_coarse,:,ineighbour) = tmp_Cmat
               else if (BSp%interp_mode == 3) then
                 Cmat(ibnd_coarse,:,ineighbour) = -tmp_Cmat
               end if
             end if


             if (any(BSp%interp_mode == [2,3])) then
               work_coeffs(:,:) = acoeffs(itc,band2it(:),:)

               call compute_subhinterp(BSp,grid,nbnd_coarse,&
&        interpolator,kdense2div,&
&        work_coeffs,tmp_Cmat,ikp_dense,&
&        interpolator%overlaps(:,:,:,ikp_dense,spin2))

               tmp_Cmat = tmp_Cmat * (vc_sqrt_qbz**2)
               Cmat(ibnd_coarse,:,ineighbour) = Cmat(ibnd_coarse,:,ineighbour) + tmp_Cmat
             end if


             if (any(BSp%interp_mode == [2,3])) then
               work_coeffs(:,:) = bcoeffs(itc,band2it(:),:)

               call compute_subhinterp(BSp,grid,nbnd_coarse,&
&        interpolator,kdense2div,&
&        work_coeffs,tmp_Cmat,ikp_dense,&
&        interpolator%overlaps(:,:,:,ikp_dense,spin2))

               tmp_Cmat = tmp_Cmat * (vc_sqrt_qbz)
               Cmat(ibnd_coarse,:,ineighbour) = Cmat(ibnd_coarse,:,ineighbour) + tmp_Cmat
             end if
               
               
             if (any(BSp%interp_mode == [2,3])) then
               work_coeffs(:,:) = ccoeffs(itc,band2it(:),:)

               call compute_subhinterp(BSp,grid,nbnd_coarse,&
&        interpolator,kdense2div,&
&        work_coeffs,tmp_Cmat,ikp_dense,&
&        interpolator%overlaps(:,:,:,ikp_dense,spin2))

               Cmat(ibnd_coarse,:,ineighbour) = Cmat(ibnd_coarse,:,ineighbour) + tmp_Cmat
             end if
           end do ! ic
         end do ! iv
       end do ! ineighbour

     else
       ! Loop over the (c, v) part of the left transition
       do iv = BSp%lomo_spin(spin1),BSp%homo_spin(spin1)
         do ic = BSp%lumo_spin(spin1),BSp%humo_spin(spin1)

           it_dense = BSp%vcks2t_interp(iv,ic,ik_dense,spin1)
           it_coarse = BSp%vcks2t(iv,ic,ik_coarse,spin1)
           ibnd_coarse = (iv-BSp%lomo_spin(spin1))*BSp%maxnbndc+(ic-BSp%lumo_spin(spin1)+1)

           ! Loop over the (c', v') part of the right transition
           do ivp = BSp%lomo_spin(spin2),BSp%homo_spin(spin2)
             do icp = BSp%lumo_spin(spin2),BSp%humo_spin(spin2)

               itp_dense = BSp%vcks2t_interp(ivp,icp,ikp_dense,spin2)
               itp_coarse = BSp%vcks2t(ivp,icp,ikp_coarse,spin2)
               ibndp_coarse = (ivp-Bsp%lomo_spin(spin2))*BSp%maxnbndc+(icp-BSp%lumo_spin(spin2)+1)
               ! Now we now it_dense, and itp_dense

               idivp = kdense2div(ikp_dense)

               btemp = czero; ctemp = czero

               ! MG TODO: This way of looping is not optimal
               do ineighbour = 1,interpolator%nvert
                 itc = interpolator%corresp(it_coarse,ineighbour,spin1)            
 
                 do ineighbourp = 1,interpolator%nvert

                   do iv1 = BSp%lomo_spin(spin2),BSp%homo_spin(spin2)
                     do ic1 = BSp%lumo_spin(spin2),BSp%humo_spin(spin2)

                       ibndp_coarse1 = (iv1-BSp%lomo_spin(spin2))*BSp%maxnbndc+(ic1-BSp%lumo_spin(spin2)+1)
                       indwithnb = (ineighbourp-1)*nbnd_coarse+ibndp_coarse1
                       itp_coarse1 = BSp%vcks2t(iv1,ic1,ikp_coarse,spin2)

                       select case (BSp%interp_mode)
                       case (1)
                         btemp(indwithnb) = hmat(itc,itp_coarse1,ineighbourp)
                       case (2)
                         btemp(indwithnb) = acoeffs(itc,itp_coarse1,ineighbourp)*(vc_sqrt_qbz**2) &
&                                         + bcoeffs(itc,itp_coarse1,ineighbourp)*(vc_sqrt_qbz) &
&                                         + ccoeffs(itc,itp_coarse1,ineighbourp)
                       case (3)
                         ! Diff between divergence and hmat
                         btemp(indwithnb) = acoeffs(itc,itp_coarse1,ineighbourp)*(vc_sqrt_qbz**2) &
&                                         + bcoeffs(itc,itp_coarse1,ineighbourp)*(vc_sqrt_qbz) &
&                                         + ccoeffs(itc,itp_coarse1,ineighbourp) &
&                                         - hmat(itc,itp_coarse1,ineighbourp)
                       case default
                         MSG_ERROR("Wrong Bsp%interp_mode")
                       end select

                       ctemp(indwithnb) = &
&                        interpolator%overlaps(iv1,ivp,ineighbourp,ikp_dense,spin2) &
&                        * GWPC_CONJG(interpolator%overlaps(ic1,icp,ineighbourp,ikp_dense,spin2)) &
&                        *interpolator%interp_factors(ineighbourp,idivp)
                     end do ! ic1
                   end do !iv1

                 end do !ineighbourp
                 Cmat(ibnd_coarse,ibndp_coarse,ineighbour) = xdotc(interpolator%nvert*nbnd_coarse,ctemp,1,btemp,1)
               end do !ineighbour

             end do !icp
           end do !ivp

         end do !ic
       end do !iv

     end if

     do iv = BSp%lomo_spin(spin1),BSp%homo_spin(spin1)
       do ic = BSp%lumo_spin(spin1),BSp%humo_spin(spin1)
         it_dense = BSp%vcks2t_interp(iv,ic,ik_dense,spin1)
         it_coarse = BSp%vcks2t(iv,ic,ik_coarse,spin1)
         ibnd_coarse = (iv-BSp%lomo_spin(spin1))*BSp%maxnbndc+(ic-BSp%lumo_spin(spin1)+1)

         idiv = kdense2div(ik_dense)

         do ivp = BSp%lomo_spin(spin2),BSp%homo_spin(spin2)
           do icp = BSp%lumo_spin(spin2),BSp%humo_spin(spin2)
             itp_dense = BSp%vcks2t_interp(ivp,icp,ikp_dense,spin2)
             itp_coarse = BSp%vcks2t(ivp,icp,ikp_coarse,spin2)
             ibndp_coarse = (ivp-Bsp%lomo_spin(spin2))*BSp%maxnbndc+(icp-BSp%lumo_spin(spin2)+1)

             btemp = czero; ctemp = czero

             do ineighbour = 1,interpolator%nvert
               do iv1 = BSp%lomo_spin(spin1),BSp%homo_spin(spin1)
                 do ic1 = BSp%lumo_spin(spin1),BSp%humo_spin(spin1)
                   ibnd_coarse1 = (iv1-BSp%lomo_spin(spin1))*BSp%maxnbndc+(ic1-BSp%lumo_spin(spin1)+1)
                   it_dense1 = BSp%vcks2t_interp(iv1,ic1,ik_dense,spin1)
                   indwithnb = (ineighbour-1)*nbnd_coarse+ibnd_coarse1

                   btemp(indwithnb) = Cmat(ibnd_coarse1,ibndp_coarse,ineighbour)

                   ctemp(indwithnb) = GWPC_CONJG(interpolator%overlaps(iv1,iv,ineighbour,ik_dense,spin1)) &
&                                    *interpolator%overlaps(ic1,ic,ineighbour,ik_dense,spin1) &
&                                    *interpolator%interp_factors(ineighbour,idiv)
                 end do !ic1
               end do !iv1
             end do !ineighbour

             ! Save interpolated value.
             hinterp(it_dense,itp_dense) = xdotc(interpolator%nvert*nbnd_coarse,ctemp,1,btemp,1)
             !DOT_PRODUCT(ctemp,btemp)

           end do !icp
         end do !ivp

       end do !ic
     end do !iv

   end do !ikp
 end do !ik

 ABI_FREE(btemp)
 ABI_FREE(ctemp)
 ABI_FREE(Cmat)
 ABI_FREE(band2it)

 if(newway) then
   ABI_FREE(tmp_Cmat)
   ABI_FREE(work_coeffs)
 end if

 hinterp = hinterp*factor
 
 call timab(696,2,tsec)

 !DBYG
 writeh = .False.
 if (writeh) then
   dump_unt = 991+BSp%interp_mode
   call wrtout(dump_unt,'Interpolated Reasonant Hamiltonian matrix elements: ',"PERS")
   call wrtout(dump_unt,'    k v  c  s     k" v" c" s"       H',"PERS")
   do itp_dense=1,BSp%nreh_interp(spin1)      
     ikp_dense = Bsp%Trans_interp(itp_dense,spin1)%k
     ivp       = Bsp%Trans_interp(itp_dense,spin1)%v
     icp       = Bsp%Trans_interp(itp_dense,spin1)%c
     do it_dense=1,BSp%nreh_interp(spin2)
       ik_dense = Bsp%Trans_interp(it_dense,spin2)%k
       iv       = Bsp%Trans_interp(it_dense,spin2)%v
       ic       = Bsp%Trans_interp(it_dense,spin2)%c
       http = hinterp(it_dense,itp_dense)
       !if (ABS(http) > tol3) then
       write(msg,'(2(i0,1x),2(i5,3i3,3x),2f24.20)')it_dense,itp_dense,ik_dense,iv,ic,spin1,ikp_dense,ivp,icp,spin2,http
       call wrtout(dump_unt,msg,"PERS")
       !end if
     end do
   end do
 end if
 !ENDDBYG

end subroutine compute_hinterp
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_interp_matmul
!! NAME
!! haydock_interp_matmul
!!
!! FUNCTION
!! Compute matrix-vector product Hmat * phi by interpolating coarse Hmat
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the BS calculation
!!  hsize_coarse = Size of the coarse hamiltonian
!!  hsize_dense = Size of the dense hamiltonian
!!  hmat(hsize_coarse,hsize_coarse,8) = coarse hamiltonian
!!  diag_dense(hsize_dense) = Diagonal with the interpolated transition energies
!!  phi(hsize_dense) = ket on which apply the matrix
!!  grid <double_grid_t> = Correspondence between coarse and dense k-mesh.
!!  nbnd_coarse
!!  corresp(hsize_coarse,8)
!!  overlaps(ntrans,nbnd_coarse,8) = Overlap coefficients.
!!  interp_factors(BSp%nreh(1),8,grid%ndiv) = Interpolation factors
!!  indices(BSp%nreh(1),grid%ndiv)
!!  hinterp(hsize_dense,hsize_dense) = Interpolated Hamiltonian
!!
!! OUTPUT
!!  hphi(hsize_dense) = Interp(hmat)*phi
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine haydock_interp_matmul(BSp,hsize_coarse,hsize_dense,hmat,diag_dense,phi,hphi,grid,&
&   nbnd_coarse,interpolator,div2kdense,kdense2div,hinterp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_interp_matmul'
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize_coarse,hsize_dense,nbnd_coarse !,ntrans
 type(excparam),intent(in) :: BSp
 type(double_grid_t),intent(in) :: grid
 type(interpolator_t),intent(in) :: interpolator
!arrays
! integer,intent(in) :: corresp(hsize_coarse,8)
! integer,intent(in) :: indices(BSp%nreh(1),grid%ndiv)
! real(dp),intent(in) :: interp_factors(BSp%nreh(1),8,grid%ndiv)
 integer,intent(in) :: div2kdense(grid%nbz_coarse,grid%ndiv), kdense2div(grid%nbz_dense) 
 complex(dpc),intent(in) :: phi(hsize_dense)
 complex(dpc),intent(in) :: hmat(hsize_coarse,hsize_coarse,interpolator%nvert)
 complex(dpc),intent(in) :: diag_dense(hsize_dense)
! complex(gwpc),intent(in) :: overlaps(ntrans,nbnd_coarse,8)
 complex(dpc),intent(out) :: hphi(hsize_dense)
 complex(dpc),optional,intent(in) :: hinterp(hsize_dense,hsize_dense)

!Local variables ------------------------------
!scalars
 integer :: itt,ik_dense,ik_coarse,it_coarse
 integer :: ic,iv,iv1,ic1, ibnd_coarse
 integer :: ibnd_coarse1
 integer :: ineighbour,idense,ikpt
 integer :: my_k1,my_k2,ind_with_nb,is, is1
 real(dp) :: factor
 complex(dpc) :: tmp
 logical,parameter :: use_blas=.True.
!arrays
 integer :: allindices(nbnd_coarse)
 real(dp) :: tsec(2)
 complex(dpc) :: allp(hsize_coarse,interpolator%nvert), test(hsize_coarse)
 complex(dpc) :: ophi(grid%nbz_dense,interpolator%nvert,nbnd_coarse)
 complex(dpc),allocatable :: b(:), c(:),A(:,:)
 complex(dpc),allocatable :: tmp_array(:), tmp_array2(:,:)

!************************************************************************

 call timab(697,1,tsec)

 factor = one/grid%ndiv

 ! Apply the diagonal part of the matrix.
 hphi = diag_dense*phi

 if (any(BSp%interp_mode == [2,3])) then
   ! hphi = hphi + MATMUL(hinterp,phi)
   call xgemv('N',hsize_dense,hsize_dense,cone,hinterp,hsize_dense,phi,1,cone,hphi,1)
   if (BSp%interp_mode == 2) then
     call timab(697,2,tsec); return ! We are done
   end if
 end if

 ! Outer index : k point in the dense zone
 ! Sum over vc
 ! Index of result : k point in the dense zone, v2,c2,neighbour

 ! Parallelization on nbz in the coarse mesh !
 my_k1 = 1
 my_k2 = grid%nbz_coarse

 ABI_MALLOC(A,(interpolator%nvert*nbnd_coarse,nbnd_coarse))
 ABI_MALLOC(b,(nbnd_coarse))
 ABI_MALLOC(c,(interpolator%nvert*nbnd_coarse))

 c = czero; ophi = czero

 !$OMP PARALLEL DO DEFAULT(none) &
 !$OMP PRIVATE(A,b,c,is1,iv1,ic1,ibnd_coarse,ibnd_coarse1,itt,allindices,is,iv,ic,idense,ineighbour,ind_with_nb) &
 !$OMP SHARED(BSp,phi,nbnd_coarse,ophi,grid,interpolator)
 do ik_dense = 1,grid%nbz_dense
   ! if( ik_dense is not in my set of k-points)
   !   ! continue
   !
   do is1 = 1, BSp%nsppol
     do iv1 = BSp%lomo_spin(is1),Bsp%homo_spin(is1)
       do ic1 = BSp%lumo_spin(is1),Bsp%humo_spin(is1)
         ibnd_coarse = (iv1-BSp%lomo_spin(is1))*BSp%maxnbndc+(ic1-BSp%lumo_spin(is1)+1)
         itt = BSp%vcks2t_interp(iv1,ic1,ik_dense,is1)
         allindices(ibnd_coarse) = itt
       end do !ic1
     end do !iv1
   end do !is1

   b(:) = phi(allindices(:))
  
   do is = 1, BSp%nsppol
     do iv = BSp%lomo_spin(is),Bsp%homo_spin(is)
       do ic = BSp%lumo_spin(is),Bsp%humo_spin(is)
         ibnd_coarse = (iv-BSp%lomo_spin(is))*BSp%maxnbndc+(ic-BSp%lumo_spin(is)+1)
         idense = Bsp%vcks2t_interp(iv,ic,ik_dense,is)

         do ineighbour = 1,interpolator%nvert
           ind_with_nb = (ineighbour-1)*(nbnd_coarse)+ibnd_coarse

           !A(ind_with_nb,:) = overlaps(allindices(:),ibnd_coarse,ineighbour)

           ! Should be optimized !!!
           do iv1 = BSp%lomo_spin(is),Bsp%homo_spin(is)
             do ic1 = BSp%lumo_spin(is),Bsp%humo_spin(is)
               ibnd_coarse1 = (iv1-BSp%lomo_spin(is))*BSp%maxnbndc+(ic1-BSp%lumo_spin(is)+1)
               A(ind_with_nb,ibnd_coarse1) = GWPC_CONJG(interpolator%overlaps(iv,iv1,ineighbour,ik_dense,is)) &
&                                          *interpolator%overlaps(ic,ic1,ineighbour,ik_dense,is)
             end do !ic1
           end do !iv1
         end do !ineighbour
       end do !ic
     end do !iv
   end do !is

   if(use_blas) then
     call xgemv('N',interpolator%nvert*nbnd_coarse,nbnd_coarse,cone,A,interpolator%nvert*nbnd_coarse,b,1,czero,c,1)
   else
     c = MATMUL(A,b)
   end if

   do is = 1, BSp%nsppol
     do iv = BSp%lomo_spin(is),BSp%homo_spin(is)
       do ic = BSp%lumo_spin(is),BSp%humo_spin(is)
         ibnd_coarse = (iv-BSp%lomo_spin(is))*BSp%maxnbndc+(ic-BSp%lumo_spin(is)+1)
         do ineighbour = 1,interpolator%nvert
           ind_with_nb = (ineighbour-1)*(nbnd_coarse)+ibnd_coarse
           ophi(ik_dense,ineighbour,ibnd_coarse) = c(ind_with_nb)
         end do !ineighbour
       end do !ic
     end do !iv
   end do !is

 end do !ik_dense
 !$OMP END PARALLEL DO

 ABI_FREE(A)
 ABI_FREE(b)
 ABI_FREE(c)

 !call xmpi_sum_(ophi,comm,ierr)

 ! Outer index : k,v,c in the coarse zone, ineighbour
 ! Sum over all k-dense relative to one coarse point
 ! Index of result : k,v,c in the coarse zone, ineighbour

 ABI_MALLOC(b,(grid%ndiv))
 ABI_MALLOC(c,(grid%ndiv))

 allp = czero

 !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
 !$OMP PRIVATE(is, it_coarse, ibnd_coarse, ineighbour, b, c, tmp,ik_coarse) &
 !$OMP SHARED(allp,BSp,ophi,interpolator,div2kdense)
 do is = 1, BSp%nsppol
   do ineighbour = 1,interpolator%nvert

     do it_coarse = 1, BSp%nreh(is)
       ibnd_coarse = (Bsp%trans(it_coarse,is)%v-BSp%lomo_spin(is))*BSp%maxnbndc+&
&            (BSp%Trans(it_coarse,is)%c-BSp%lumo_spin(is)+1)
       ik_coarse = BSp%trans(it_coarse,is)%k
       !b(:) = interp_factors(it_coarse,ineighbour,:) 
       b(:) = interpolator%interp_factors(ineighbour,:) 
       !c(:) = ophi(indices(it_coarse,:),ineighbour,ibnd_coarse)
       c(:) = ophi(div2kdense(ik_coarse,:),ineighbour,ibnd_coarse)
       tmp = DOT_PRODUCT(b,c)
       allp(it_coarse,ineighbour) = tmp     
     end do

   end do
 end do
 !$OMP END PARALLEL DO

 !call xmpi_sum_(allp,comm,ierr)

 ABI_FREE(b)
 ABI_FREE(c)

 ABI_MALLOC(tmp_array,(hsize_coarse))
 ABI_MALLOC(tmp_array2,(hsize_coarse,hsize_coarse))
 tmp_array(:) = czero
 tmp_array2(:,:) = czero

 test = czero

 ! Second step : Multiplication by hmat
 ! Note: OMP is deactivated since this would require large copies in stack of
 ! each thread !
 !!!$OMP PARALLEL DO DEFAULT(none) &
 !!!$OMP PRIVATE(ineighbour,tmp_array,tmp_array2) 
 !!!$OMP SHARED(factor,hmat,allp,hsize_coarse,std_out) reduction(+:test)
 do ineighbour = 1,interpolator%nvert
   if(use_blas) then
     !call xgemv('N',hsize_coarse,hsize_coarse,cone,factor*(hmat(:,:,ineighbour)),hsize_coarse,allp(:,ineighbour),1,czero,tmp_array,1)
     tmp_array2 = hmat(:,:,ineighbour)
     tmp_array2 = factor*tmp_array2
     call xgemv('N',hsize_coarse,hsize_coarse,cone,tmp_array2,hsize_coarse,allp(:,ineighbour),1,czero,tmp_array,1)
     test = test + tmp_array 
   else 
     test = test+MATMUL(factor*(hmat(:,:,ineighbour)),allp(:,ineighbour))
   end if
 end do
 !!!$OMP END PARALLEL DO
 
 ABI_FREE(tmp_array)
 ABI_FREE(tmp_array2)

 ! Outer index : ineighbour
 ! Sum over all v c
 ! Index of result : ineighbour, k_dense, v,c
 ABI_MALLOC(A,(nbnd_coarse,nbnd_coarse))
 ABI_MALLOC(b,(nbnd_coarse))
 ABI_MALLOC(c,(nbnd_coarse))
 c = czero

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
!$OMP PRIVATE(ineighbour,ik_dense,ibnd_coarse1,is1,iv1,ic1,is,iv,ic,A,b,c,ibnd_coarse,ik_coarse,itt,idense) &
!$OMP SHARED(test,ophi,BSp,grid,nbnd_coarse,interpolator)
 do ineighbour = 1,interpolator%nvert
   do ik_dense = 1,grid%nbz_dense

     do is1 = 1, Bsp%nsppol
       do iv1 = Bsp%lomo_spin(is1),Bsp%homo_spin(is1)
         do ic1 = BSp%lumo_spin(is1), Bsp%humo_spin(is1)
           ibnd_coarse = (iv1-BSp%lomo_spin(is1))*BSp%maxnbndc+(ic1-BSp%lumo_spin(is1)+1)

           ik_coarse = grid%dense_to_coarse(ik_dense)
           itt = BSp%vcks2t(iv1,ic1,ik_coarse,is1)
           b(ibnd_coarse) = test(interpolator%corresp(itt,ineighbour,is1))
         end do ! ic1
       end do ! iv1
     end do ! is1

     do is = 1, BSp%nsppol
       do iv = BSp%lomo_spin(is),Bsp%homo_spin(is)
         do ic = BSp%lumo_spin(is),BSp%humo_spin(is)
           ibnd_coarse = (iv-BSp%lomo_spin(is))*Bsp%maxnbndc+(ic-BSp%lumo_spin(is)+1)
           idense = BSp%vcks2t_interp(iv,ic,ik_dense,is)

           !A(ibnd_coarse,:) = CONJG(overlaps(idense,:,ineighbour))

           ! Should be optimized !!!
           do iv1 = BSp%lomo_spin(is),Bsp%homo_spin(is)
             do ic1 = BSp%lumo_spin(is),Bsp%humo_spin(is)
               ibnd_coarse1 = (iv1-BSp%lomo_spin(is))*BSp%maxnbndc+(ic1-BSp%lumo_spin(is)+1)
               A(ibnd_coarse,ibnd_coarse1) = (interpolator%overlaps(iv1,iv,ineighbour,ik_dense,is)) &
&                                          *GWPC_CONJG(interpolator%overlaps(ic1,ic,ineighbour,ik_dense,is))
             end do !ic1
           end do !iv1 
         end do ! ic
       end do !iv
     end do !is

     if(use_blas) then
       call xgemv('N',nbnd_coarse,nbnd_coarse,cone,A,nbnd_coarse,b,1,czero,c,1)
     else
       c = MATMUL(A,b)
     end if

     do is = 1, BSp%nsppol
       do iv = BSp%lomo_spin(is),Bsp%homo_spin(is)
         do ic = BSp%lumo_spin(is),BSp%humo_spin(is)
           ibnd_coarse = (iv-BSp%lomo_spin(is))*BSp%maxnbndc+(ic-BSp%lumo_spin(is)+1)
           idense = Bsp%vcks2t_interp(iv,ic,ik_dense,is)
           !ophi(ik_dense,ineighbour,ibnd_coarse) = c(idense)
           ophi(ik_dense,ineighbour,ibnd_coarse) = c(ibnd_coarse)
         end do
       end do
     end do

   end do ! ik_dense
 end do ! ineighbour
!$OMP END PARALLEL DO

 !call xmpi_sum_(ophi,comm,ierr)
 
 ABI_FREE(A)
 ABI_FREE(b)
 ABI_FREE(c)

 ! Outer indices : it_dense 
 ! Sum over neighbours
 ! Index of result : it_dense (ik,ic,iv dense)

 ABI_MALLOC(b,(interpolator%nvert))
 ABI_MALLOC(c,(interpolator%nvert))

 do is = 1, BSp%nsppol
!!!Disable OpenMP since it leads to segmentation faults !
!!!OMP PARALLEL DO DEFAULT(none) &
!!!$OMP PRIVATE(is,itt,ik_dense,ic,iv,ik_coarse,it_coarse,ibnd_coarse,ix,ikpt,b,c) &
!!!$OMP SHARED(BSp,grid,kdense2div,ophi,hphi,std_out) 
   do itt = 1,BSp%nreh_interp(is)
    ! From itt -> ik_ibz,ic,iv
    ik_dense = BSp%Trans_interp(itt,is)%k
    ic = BSp%Trans_interp(itt,is)%c
    iv = BSp%Trans_interp(itt,is)%v

    ! From ik_ibz in the dense mesh -> indices_dense
    ik_coarse = grid%dense_to_coarse(ik_dense)
    it_coarse = BSp%vcks2t(iv,ic,ik_coarse,is)

    ibnd_coarse = (iv-BSp%lomo_spin(is))*BSp%maxnbndc+(ic-BSp%lumo_spin(is)+1)

    ikpt = kdense2div(ik_dense)
    !ikpt = -1
    !do ix = 1,grid%ndiv
    !  if (indices(it_coarse,ix) == ik_dense) then
    !    ikpt = ix
    !    exit 
    !  end if
    !end do
    !ABI_CHECK(ikpt/=-1,"Cannot find ik_dense")

    !b = interp_factors(it_coarse,:,ikpt)
    b = interpolator%interp_factors(:,ikpt)
    c =  ophi(ik_dense,:,ibnd_coarse)

    hphi(itt) = hphi(itt) + xdotc(interpolator%nvert, b, 1, c, 1)
   end do 
!!$OMP END PARALLEL DO
 end do

 ABI_FREE(b)
 ABI_FREE(c)

 call timab(697,2,tsec)

end subroutine haydock_interp_matmul
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_herm
!! NAME
!! haydock_herm
!!
!! FUNCTION
!!  Reads the excitonic Hamiltonian from file and construct the Lanczos set of vectors 
!!  by iterative matrix-vector multiplications.
!!
!! INPUTS
!! BSp<excparam>=Parameters for the Bethe-Salpeter calculation.
!! BS_files<excparam>=Files associated to the bethe_salpeter code.
!! Cryst<crystal_t>=Info on the crystalline structure.
!! hize=Size of the excitonic matrix.
!! my_t1,my_t2=First and last columns treated by this node.
!! hmat(hsize,my_t1:my_t2)=Excitonic matrix.
!! nkets=Number of starting vectors for Haydock method.
!! kets(hsize,nkets)=The kets in the eh representation.
!! comm=MPI communicator.
!!
!! OUTPUT
!!  green(BSp%nomega,nkets)=
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine haydock_herm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hmat,nkets,kets,green,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_herm'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize,my_t1,my_t2,nkets,comm
 type(crystal_t),intent(in) :: Cryst
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
!arrays
 complex(dp),intent(out) :: green(BSp%nomega,nkets)
 complex(dpc),intent(in) :: hmat(hsize,my_t1:my_t2),kets(hsize,nkets)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: inn,nproc,my_rank,ierr
 integer :: niter_file,niter_max,niter_done,nsppol,iq,my_nt,term_type
 integer :: n_all_omegas
 real(dp) :: norm,nfact
 logical :: can_restart,is_converged
 complex(dpc) :: factor
 character(len=500) :: msg
 character(len=fnlen),parameter :: tag_file="_HAYDR_SAVE"
 character(len=fnlen) :: restart_file,out_file 
 type(haydock_type) :: haydock_file
!arrays
 real(dp),allocatable :: bb_file(:)
 real(dp),allocatable :: bb(:)
 complex(dpc),allocatable :: aa(:),phi_nm1(:),phi_n(:),hphi_n(:),hphi_nm1(:)
 complex(dpc),allocatable :: aa_file(:),phi_n_file(:),phi_nm1_file(:)
 complex(dpc),allocatable :: ket0(:)
 complex(dpc),allocatable :: all_omegas(:)
 complex(dpc),allocatable :: green_temp(:,:)
 logical :: check(2)
 
!************************************************************************

 nproc  = xcomm_size(comm)
 my_rank= xcomm_rank(comm)
 nsppol = Hdr_bse%nsppol

 my_nt = my_t2-my_t1+1
 ABI_CHECK(my_nt>0,"One of the processors has zero columns")

 write(msg,'(a,i0)')' Haydock algorithm with MAX number of iterations: ',BSp%niter
 call wrtout(std_out,msg,"COLL")
 !
 ! Select the terminator for the continued fraction.
 term_type=0; if (Bsp%hayd_term>0) term_type=1 
 call wrtout(std_out,sjoin("Using terminator type: ",itoa(term_type)),"COLL")
 !
 ! Check for presence of the restart file.
 can_restart=.FALSE.
 if ( BS_files%in_haydock_basename /= BSE_NOFILE) then
   restart_file = TRIM(BS_files%in_haydock_basename)
   if (file_exists(restart_file) ) then
     can_restart=.TRUE.
     msg = " Restarting Haydock calculation from file: "//TRIM(restart_file)
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")
   else 
     can_restart=.FALSE.
     MSG_WARNING("Cannot find restart file: "//TRIM(restart_file))
   end if
 end if
 !
 ! Open the file and write basic dimensions and info.
 if (my_rank==master) then
   out_file = TRIM(BS_files%out_basename)//TRIM(tag_file)
   call open_haydock(out_file,haydock_file) 
   haydock_file%hsize = hsize
   haydock_file%use_coupling = Bsp%use_coupling
   haydock_file%op = BSE_HAYD_IMEPS
   haydock_file%nq = nkets
   haydock_file%broad = Bsp%broad
   call write_dim_haydock(haydock_file)
 end if
 !
 ! Calculate green(w) for the different starting points.
 green=czero
 do iq=1,nkets
   ABI_MALLOC(ket0,(hsize))
   ket0=kets(:,iq)
   !
   niter_file=0
   if (can_restart) then
     call haydock_restart(BSp,restart_file,BSE_HAYD_IMEPS,iq,hsize,&
&      niter_file,aa_file,bb_file,phi_nm1_file,phi_n_file,comm)
   end if 
   !
   ! For n>1, we have:
   !  1) a_n = <n|H|n>
   !  2) b_n = || H|n> - a_n|n> -b_{n-1}|n-1> ||
   !  3) |n+1> = [H|n> -a_n|n> -b_{n-1}|n-1>]/b_n
   !
   ! The sequences starts with |1> normalized to 1 and b_0 =0, therefore:
   !  a_1 = <1|H|1>
   !  b_1 = || H|1> - a_1|1> ||
   !  |2> = [H|1> - a_1|1>]/b_1
   !
   ABI_MALLOC(hphi_n,(hsize))
   ABI_MALLOC(hphi_nm1,(hsize))
   ABI_MALLOC(phi_nm1,(my_nt))
   ABI_MALLOC(phi_n,(my_nt))

   niter_max = niter_file + Bsp%niter
   ABI_MALLOC(aa,(niter_max))
   ABI_MALLOC(bb,(niter_max))
   aa=czero; bb=zero

   if (niter_file==0) then       ! Calculation from scratch.
     phi_nm1=ket0(my_t1:my_t2)   ! Select the slice treated by this node.
     norm = DZNRM2(hsize,ket0,1) ! Normalization  
     phi_nm1=phi_nm1/norm      
                                                                                
     ! hphi_n = MATMUL(hmat,phi_nm1)
     call xgemv('N',hsize,my_nt,cone,hmat,hsize,phi_nm1,1,czero,hphi_n,1)
     call xmpi_sum(hphi_n,comm,ierr)

     aa(1)=xdotc(my_nt,phi_nm1,1,hphi_n(my_t1:),1)
     call xmpi_sum(aa(1:1),comm,ierr)

     phi_n = hphi_n(my_t1:my_t2) - aa(1)*phi_nm1

     bb(1) = xdotc(my_nt,phi_n,1,phi_n,1)
     call xmpi_sum(bb(1:1),comm,ierr)
     bb(1) = SQRT(bb(1))

     phi_n = phi_n/bb(1)
     niter_done=1

   else ! Use the previous a and b.
     niter_done=niter_file
     aa(1:niter_done) = aa_file
     bb(1:niter_done) = bb_file
     phi_nm1=phi_nm1_file(my_t1:my_t2)   ! Select the slice treated by this node.
     phi_n  =phi_n_file  (my_t1:my_t2)   
   end if

   if (can_restart) then
     ABI_FREE(aa_file)
     ABI_FREE(bb_file)
     ABI_FREE(phi_nm1_file)
     ABI_FREE(phi_n_file)
   end if

   ! Multiplicative factor (k-point sampling and unit cell volume)  
   ! TODO be careful with the spin here
   ! TODO four_pi comes from the coulomb term 1/|q| is already included in the 
   ! oscillators hence the present approach wont work if a cutoff interaction is used.
   nfact = -four_pi/(Cryst%ucvol*BSp%nkbz)
   if (nsppol==1) nfact=two*nfact 

   factor = nfact*(DZNRM2(hsize,ket0,1)**2)

   ! Which quantity should be checked for convergence?
   check = (/.TRUE.,.TRUE./) 
   if (ABS(Bsp%haydock_tol(2)-one)<tol6) check = (/.TRUE. ,.FALSE./) 
   if (ABS(Bsp%haydock_tol(2)-two)<tol6) check = (/.FALSE.,.TRUE./) 

   ! Create new frequencies "mirror" in negative range to add 
   ! their contributions. Can be improved by computing only once
   ! zero frequency, but loosing clearness
   n_all_omegas = 2*BSp%nomega

   ABI_MALLOC(all_omegas,(n_all_omegas))
   ! Put all omegas with frequency > 0 in table
   all_omegas(BSp%nomega+1:n_all_omegas) = BSp%omega
   ! Put all omegas with frequency < 0
   ! Warning, the broadening must be kept positive
   all_omegas(1:BSp%nomega) = -DBLE(BSp%omega(BSp%nomega:1:-1)) + j_dpc*AIMAG(BSp%omega(BSp%nomega:1:-1))   

   ABI_MALLOC(green_temp,(n_all_omegas,nkets))

   ! Calling haydock_herm_algo with green_temp with full range of frequencies
   call haydock_herm_algo(niter_done,niter_max,n_all_omegas,all_omegas,BSp%haydock_tol(1),check,hsize,&
&    my_t1,my_t2,hmat,factor,term_type,aa,bb,phi_nm1,phi_n,green_temp(:,iq),inn,is_converged,comm)

   ! Computing result from two ranges of frequencies
   ! The real part is added, the imaginary part is substracted
   green(:,iq) = green_temp(BSp%nomega+1:n_all_omegas,iq)+CONJG(green_temp(BSp%nomega:1:-1,iq))

   ABI_FREE(all_omegas)
   ABI_FREE(green_temp)
   !
   ! Save the a"s and the b"s for possible restarting.
   ! 1) Info on the Q.
   ! 2) Number of iterations performed.
   ! 3) do iter=1,niter_performed 
   !      aa(iter),bb(iter)
   !    end do
   ! 4) |n-1>
   !    |n>
   !
   hphi_nm1 = czero
   hphi_nm1(my_t1:my_t2) = phi_nm1
   call xmpi_sum_master(hphi_nm1,master,comm,ierr)

   hphi_n = czero
   hphi_n(my_t1:my_t2) = phi_n
   call xmpi_sum_master(hphi_n,master,comm,ierr)

   if (my_rank==master) then 
     ! Write data for restarting
     call write_haydock(haydock_file, hsize, Bsp%q(:,iq), aa, bb, hphi_n, hphi_nm1, MIN(inn,niter_max), factor)
   end if
  
   ABI_FREE(hphi_n)
   ABI_FREE(hphi_nm1)
   ABI_FREE(phi_nm1)
   ABI_FREE(phi_n)
   ABI_FREE(aa)
   ABI_FREE(bb)
   ABI_FREE(ket0)
 end do ! iq

 if (my_rank==master) call close_haydock(haydock_file)

 call xmpi_barrier(comm)

end subroutine haydock_herm
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_herm_algo
!! NAME
!! haydock_herm_algo
!!
!! FUNCTION
!!  Haydock algorithm for Hermitian matrix
!!
!! INPUTS
!!  niter_done=Number of iterations already performed (0 if the run starts from scratch).
!!  niter_max=Max number of iterations. Always > niter_done
!!  nomega=Number of Frequency points for the evaluation of the matrix element.
!!  omega(nomega)=Frequency set (imaginary part is already included).
!!  tol_iter=Tolerance used to stop the algorithm.
!!  check(2)=Logical flags to specify where both the real and the imaginary part of the 
!!    matrix elements of the Green functions have to be checked for convergence. 
!!  hsize=Size of the blocks.
!!  my_t1,my_t2=Indices of the first and last column stored treated by this done.
!!  term_type=0 if no terminator is used, 1 otherwise.
!!  hmat(hsize,my_t1:my_t2)=The columns of the block.
!!  factor
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  green(nomega)=Output matrix elements.
!!  inn=Last iteration performed.
!!  is_converged=.TRUE. of the algorithm converged.
!!
!! SIDE EFFECTS
!!  phi_nm1(my_t2-my_t1+1), phi_n(my_t2-my_t1+1)
!!    input: vectors used to initialize the iteration
!!    output: the vectors obtained in the last iteration 
!!  aa(niter_max) and bb(niter_max)
!!    if niter_done>0: aa(1:niter_done), bb(1:niter_done) store the coefficients of the previous run.
!!    when the routine returns aa(1:inn) and bb(1:inn) contain the matrix elements of the tridiagonal form.
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine haydock_herm_algo(niter_done,niter_max,nomega,omega,tol_iter,check,hsize,my_t1,my_t2,hmat,&
&  factor,term_type,aa,bb,phi_nm1,phi_n,green,inn,is_converged,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_herm_algo'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter_max,niter_done,nomega
 integer,intent(in) :: comm,hsize,my_t1,my_t2,term_type
 integer,intent(out) :: inn
 logical,intent(out) :: is_converged
 real(dp),intent(in) :: tol_iter
 complex(dpc),intent(in) :: factor
!arrays
 real(dp),intent(inout) :: bb(niter_max)
 complex(dpc),intent(out) :: green(nomega)
 complex(dpc),intent(in) :: omega(nomega) 
 complex(dpc),intent(inout) :: aa(niter_max)
 complex(dpc),intent(in) :: hmat(hsize,my_t1:my_t2)
 complex(dpc),intent(inout) :: phi_nm1(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_n  (my_t2-my_t1+1)
 logical,intent(in) :: check(2)

!Local variables ------------------------------
!scalars
 integer :: ierr,my_nt,niter_min,nconv
 character(len=500) :: msg
 logical,parameter :: force_real=.TRUE.
!arrays
 real(dp) :: abs_err(nomega,2) !,rel_err(nomega,2)
 complex(dpc),allocatable :: oldg(:),newg(:)
 complex(dpc),allocatable :: phi_np1(:),hphi_n(:),cfact(:)
 logical :: test(2)

!************************************************************************

 ! The sequences starts with |1> normalized to 1 and b_0 =0, therefore:
 !  a_1 = <1|H|1>
 !  b_1 = || H|1> - a_1|1> ||
 !  |2> = [H|1> - a_1|1>]/b_1
 !
 ! For n>1 we have
 !  1) a_n = <n|H|n>
 !  2) b_n = || H|n> - a_n|n> -b_{n-1}|n-1> ||
 !  3) |n+1> = [H|n> -a_n|n> -b_{n-1}|n-1>]/b_n
 !
 my_nt = my_t2-my_t1+1

 ABI_MALLOC(hphi_n,(hsize))
 ABI_MALLOC(phi_np1,(my_nt))

 ABI_MALLOC(oldg,(nomega))
 ABI_MALLOC(newg,(nomega))
 ABI_MALLOC(cfact,(nomega))
 oldg=czero; newg=czero; cfact=czero 

 nconv=0
 do inn=niter_done+1,niter_max
   !
   ! hphi_n = MATMUL(hmat,phi_n)
   call xgemv('N',hsize,my_nt,cone,hmat,hsize,phi_n,1,czero,hphi_n,1)
   call xmpi_sum(hphi_n,comm,ierr)

   aa(inn) = xdotc(my_nt,phi_n,1,hphi_n(my_t1:),1)
   call xmpi_sum(aa(inn:inn),comm,ierr)
   if (force_real) aa(inn) = DBLE(aa(inn)) ! Matrix is Hermitian.

   ! |n+1> = H|n> - A(n)|n> - B(n-1)|n-1>
   phi_np1 = hphi_n(my_t1:my_t2) - aa(inn)*phi_n - bb(inn-1)*phi_nm1

   bb(inn) = xdotc(my_nt,phi_np1,1,phi_np1,1)
   call xmpi_sum(bb(inn),comm,ierr)
   bb(inn) = SQRT(bb(inn))

   phi_np1 = phi_np1/bb(inn)
   
   phi_nm1 = phi_n
   phi_n   = phi_np1

   write(msg,'(a,i0,a,3es12.4)')' Iteration number ',inn,', b_i RE(a_i) IM(a_i) ',bb(inn),REAL(aa(inn)),AIMAG(aa(inn)) 
   call wrtout(std_out,msg,"COLL")

   call continued_fract(inn,term_type,aa,bb,nomega,omega,cfact)

   newg= factor*cfact
   !
   ! Avoid spurious convergence.
   niter_min=4; if (niter_done>1) niter_min=niter_done+1
   if (inn>niter_min) then
     test=.TRUE.
     abs_err(:,1) = ABS(DBLE (newg-oldg))
     abs_err(:,2) = ABS(AIMAG(newg-oldg))
     !
     if (tol_iter>zero) then 
       ! Test on the L1 norm.
       if (check(1)) test(1) = SUM(abs_err(:,1)) < tol_iter*SUM(ABS(DBLE (newg)))
       if (check(2)) test(2) = SUM(abs_err(:,2)) < tol_iter*SUM(ABS(AIMAG(newg)))

     else 
       ! Stringent test for each point.
       if (check(1)) test(1) = ALL( abs_err(:,1) < -tol_iter*ABS(DBLE (newg))) 
       if (check(2)) test(2) = ALL( abs_err(:,2) < -tol_iter*ABS(AIMAG(newg)))
     end if
     !
     if (ALL(test)) then 
       nconv = nconv+1
     else 
       nconv = 0
     end if
     if (nconv==2) then 
       write(msg,'(a,es10.2,a,i0,a)')&
&        " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after ",inn," iterations." 
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       EXIT
     end if
   end if

   oldg = newg
 end do ! inn

 green = newg
 if (nconv/=2) then
   write(msg,'(a,es10.2,a,i0,a)')&
&    " WARNING: Haydock algorithm did not converge within ",tol_iter," after ",niter_max," iterations."
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write(yaml_out,"(8a)")ch10,&
&    "--- !HaydockConvergenceWarning",ch10,&
&    "message: | ",ch10,TRIM(indent(msg)),ch10,&
&    "..."
 end if

 is_converged = (nconv==2)

 ABI_FREE(oldg)
 ABI_FREE(newg)
 ABI_FREE(cfact)
 ABI_FREE(hphi_n)
 ABI_FREE(phi_np1)

end subroutine haydock_herm_algo
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_restart
!! NAME
!! haydock_restart
!!
!! FUNCTION
!! Restart the Haydock method from file reading the data produced in a previous run.
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the Bethe-Salpeter calculation.
!!    omega(BSp%nomega)=Frequency mesh for the macroscopic dielectric function (broadening is already included).
!!  iq_search=The index of the q-point to be searched.
!!  hsize
!!  comm=MPI communicator.
!!  nsppol
!!  restart_file
!!
!! OUTPUT
!!  niter_file=Number of iterations already performed. 0 to signal that an error occurred during the reading
!!  bb_file(:)
!!  aa_file(:)
!!  phi_n_file(:)
!!  phi_nm1_file(:)
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine haydock_restart(BSp,restart_file,ftype,iq_search,hsize,niter_file,aa_file,bb_file,phi_nm1_file,phi_n_file,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_restart'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,hsize,iq_search,ftype
 integer,intent(out) :: niter_file
 character(len=*),intent(in) :: restart_file
 type(excparam),intent(in) :: BSp
!arrays
 real(dp),allocatable,intent(out) :: bb_file(:)
 complex(dpc),allocatable,intent(out) :: aa_file(:),phi_n_file(:),phi_nm1_file(:)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: nproc,my_rank,ierr,op_file
 integer :: hsize_file,use_coupling_file
 complex(dpc) :: factor_file
 character(len=500) :: msg
 type(haydock_type) :: haydock_file

!************************************************************************

 nproc = xcomm_size(comm); my_rank= xcomm_rank(comm)

 if (my_rank==master) then
   call open_haydock(restart_file, haydock_file)

   call read_dim_haydock(haydock_file)

   if (haydock_file%op/=ftype) then
     write(msg,"(2(a,i0))")" Expecting restart file with filetype: ",ftype," but found ",op_file
     MSG_ERROR(msg)
   end if

   if (haydock_file%hsize/=hsize) then
     write(msg,"(2(a,i0))")&
&      " Rank of H_exc read from file: ",hsize_file," differs from the one used in this run: ",hsize
     MSG_ERROR(msg)
   end if

   if (haydock_file%use_coupling /= BSp%use_coupling) then
     write(msg,'(2(a,i0))')&
&      " use_coupling_file: ",use_coupling_file," differs from input file value: ",BSp%use_coupling
     MSG_ERROR(msg)
   end if

   call read_haydock(haydock_file, Bsp%q(:,iq_search), aa_file, bb_file, &
&                   phi_n_file, phi_nm1_file, niter_file, factor_file)

   if (niter_file == 0) then
     write(msg,"(a,3f8.4,3a)")&
&      " Could not find q-point: ",BSp%q(:,iq_search)," in file ",TRIM(restart_file),&
&      " Cannot restart Haydock iterations for this q-point"
     MSG_COMMENT(msg)
   else
     write(msg,'(a,i0)')" Number of iterations already performed: ",niter_file
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")

     if ( ABS(haydock_file%broad - BSp%broad) > tol6) then
       write(msg,'(2a,2(a,f8.4),a)')&
&        " Restart file has been produced with a different Lorentzian broadening: ",ch10,&
&        " broad_file: ",haydock_file%broad," input broadening: ",BSp%broad," Continuing anyway. "
       MSG_WARNING(msg)
     end if

     call close_haydock(haydock_file)
   end if
 end if
 !
 ! Master broadcasts the data.
 call xmpi_bcast(niter_file,master,comm,ierr)

 if (my_rank/=master) then 
   ABI_MALLOC(aa_file,(niter_file))
   ABI_MALLOC(bb_file,(niter_file))
   ABI_MALLOC(phi_nm1_file,(hsize))
   ABI_MALLOC(phi_n_file,(hsize))
 end if

 call xmpi_bcast(aa_file,master,comm,ierr)
 call xmpi_bcast(bb_file,master,comm,ierr)
 call xmpi_bcast(phi_nm1_file,master,comm,ierr)
 call xmpi_bcast(phi_n_file,master,comm,ierr)

end subroutine haydock_restart
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_mdf_to_tensor
!! NAME
!! haydock_mdf_to_tensor
!!
!! FUNCTION
!! Transform macroscopic dielectric function from green function to each components of the tensor in red and cart coord.
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the Bethe-Salpeter calculation.
!!    omega(BSp%nomega)=Frequency mesh for the macroscopic dielectric function (broadening is already included).
!!  Cryst=Parameters of the crystal
!!  eps(BSp%nomega,BSp%nq) = Macroscopic dielectric function to be written.
!!
!! OUTPUT
!!  tensor_cart(BSp%nomega,6) = dielectric tensor for each frequency, order (11,22,33,12,13,23) in cart. coord.
!!  tensor_red(BSp%nomega, 6) = idem in reduced coordinated
!!  ierr = 0 if the tensors have been successfully computed
!!      \= 0 if the system is ill-posed in terms of q-points (not enough or not independent q-points)
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine haydock_mdf_to_tensor(BSp,Cryst,eps,tensor_cart,tensor_red,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_mdf_to_tensor'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ierr
 type(excparam),intent(in) :: BSp
 type(crystal_t),intent(in) :: Cryst
!arrays
 complex(dpc),intent(in) :: eps(BSp%nomega,BSp%nq)
 complex(dpc),intent(out) :: tensor_cart(BSp%nomega,6), tensor_red(BSp%nomega,6)

!Local variables ------------------------------
!scalars
 integer :: iq,info
 real(dp) :: normqcart, normqred
!arrays
 integer,allocatable :: ipiv(:)
 real(dp) :: qcart(3), qtmet(3)
 real(dp) :: qred2cart(3,3),qcart2red(3,3)
 complex(dpc) :: qqcart(BSp%nq,6), qqred(BSp%nq,6)
 complex(dpc) :: b(6,BSP%nomega)

!************************************************************************

 ! Error flag
 ierr = 0

 if(BSp%nq /= 6) then
    ierr = -1
    return
 end if

 ! Transformation matrices from reduced coordinates to cartesian coordinates
 qred2cart = two_pi*Cryst%gprimd
 qcart2red = qred2cart
 call matrginv(qcart2red,3,3)
 do iq = 1, 6

   ! Computing cartesian q-vector
   qcart = MATMUL(qred2cart, BSp%q(:,iq))

   ! Computing product 'metric - qred' to form quadratic form
   qtmet = (two_pi**2)*MATMUL(Cryst%gmet, BSp%q(:,iq))
 
   ! squared norms
   normqcart = qcart(1)**2+qcart(2)**2+qcart(3)**2
   normqred = (normv(BSp%q(:,iq),Cryst%gmet,"G"))**2

   ! Compute line 'iq' for matrix in cartesian coord
   qqcart(iq,1) = (qcart(1))**2
   qqcart(iq,2) = (qcart(2))**2
   qqcart(iq,3) = (qcart(3))**2
   qqcart(iq,4) = 2*(qcart(1)*qcart(2))
   qqcart(iq,5) = 2*(qcart(1)*qcart(3))
   qqcart(iq,6) = 2*(qcart(2)*qcart(3))

   ! Compute line 'iq' for matrix in reduced coord
   qqred(iq,1) = (qtmet(1))**2
   qqred(iq,2) = (qtmet(2))**2
   qqred(iq,3) = (qtmet(3))**2
   qqred(iq,4) = 2*(qtmet(1)*qtmet(2))
   qqred(iq,5) = 2*(qtmet(1)*qtmet(3))
   qqred(iq,6) = 2*(qtmet(2)*qtmet(3))

   ! Renormalize line
   qqcart(iq,:) = qqcart(iq,:)/normqcart
   qqred(iq,:) = qqred(iq,:)/normqred
 end do

 ABI_MALLOC(ipiv,(6))

 ! Solving linear system
 b = TRANSPOSE(eps)
 call ZGESV(6,BSp%nomega,qqcart,6,ipiv,b,6,info)
 tensor_cart = TRANSPOSE(b)

 if(info /= 0) then
   ! Skipping the rest of the routine
   ierr = info
   ABI_FREE(ipiv)
   return
 end if

 b = TRANSPOSE(eps)
 call ZGESV(6,BSp%nomega,qqred,6,ipiv,b,6,info)
 tensor_red = TRANSPOSE(b)

 if(info /= 0) ierr = info

 ABI_FREE(ipiv)

end subroutine haydock_mdf_to_tensor
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_psherm
!! NAME
!! haydock_psherm
!!
!! FUNCTION
!!  Reads the excitonic Hamiltonian from file and construct the Lanczos set of vectors 
!!  by iterative matrix-vector multiplications.
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the Bethe-Salpeter calculation.
!!    omega(BSp%nomega)=Frequency mesh for the macroscopic dielectric function (broadening is already included).
!! hize
!! my_t1,my_t2
!! hreso(hsize,my_t1:my_t2)
!! hcoup(hsize,my_t1:my_t2)
!! nkets
!! kets(hsize,nkets)
!! comm=MPI communicator.
!!
!! OUTPUT
!!  green(BSp%nomega)=The imaginary part of the macroscopic dielectric function.
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine haydock_psherm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hreso,hcoup,nkets,kets,green,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_psherm'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize,my_t1,my_t2,nkets,comm
 type(crystal_t),intent(in) :: Cryst
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
!arrays
 complex(dp),intent(out) :: green(BSp%nomega,BSp%nq)
 complex(dpc),intent(in) :: hreso(hsize,my_t1:my_t2) 
 complex(dpc),intent(in) :: kets(hsize,nkets)
 complex(dpc),intent(in) :: hcoup(hsize,my_t1:my_t2)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: inn,itt,out_unt,nproc,my_rank,ierr
 integer :: niter_file,niter_max,niter_done,nsppol,iq,my_nt,term_type
 real(dp) :: ket0_hbar_norm,nfact
 logical :: can_restart,is_converged
 complex(dpc) :: factor
 character(len=fnlen),parameter :: tag_file="_HAYDC_SAVE"
 character(len=500) :: msg
 character(len=fnlen) :: restart_file,out_file
!arrays
 real(dp),allocatable :: bb_file(:)
 real(dp),allocatable :: bb(:)
 complex(dpc),allocatable :: aa(:),cc(:),phi_np1(:),phi_n(:),phi_nm1(:),cbuff(:)
 complex(dpc),allocatable :: aa_file(:),phi_n_file(:),phi_np1_file(:),cc_file(:)
 complex(dpc),allocatable :: ket0(:)
 logical :: check(2)

!************************************************************************

 MSG_WARNING("Haydock + coupling is still under development")

 nproc  = xcomm_size(comm)
 my_rank= xcomm_rank(comm)
 nsppol = Hdr_bse%nsppol

 my_nt = my_t2-my_t1+1
 ABI_CHECK(my_nt>0,"One of the processors has zero columns")

 ! Multiplicative factor (k-point sampling and unit cell volume)  
 ! TODO be careful with the spin here
 ! TODO four_pi comes from the coulomb term 1/|q| is already included in the 
 ! oscillators hence the present approach wont work if a cutoff interaction is used.
 nfact = four_pi/(Cryst%ucvol*BSp%nkbz)
 if (nsppol==1) nfact=two*nfact 

 write(msg,'(a,i0)')' Haydock algorithm with MAX number of iterations: ',BSp%niter
 call wrtout(std_out,msg,"COLL")
 !
 ! Check for presence of the restart file.
 can_restart=.FALSE.

 if ( BS_files%in_haydock_basename /= BSE_NOFILE) then
   restart_file = strcat(BS_files%in_haydock_basename,tag_file)
   if (file_exists(restart_file) ) then
     can_restart=.TRUE.
     msg = strcat(" Restarting Haydock calculation from file: ",restart_file)
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")
     MSG_ERROR("Restart is not tested")
   else 
     can_restart=.FALSE.
     MSG_WARNING(strcat("Cannot find restart file: ",restart_file))
   end if
 end if
 !
 ! Open the file and writes basic dimensions and info.
 if (my_rank==master) then 
   out_file = TRIM(BS_files%out_basename)//TRIM(tag_file)
   if (open_file(out_file,msg,newunit=out_unt,form="unformatted") /= 0) then
     MSG_ERROR(msg)
   end if
   ! write header TODO: standardize this part.
   write(out_unt)hsize,Bsp%use_coupling,BSE_HAYD_IMEPS,nkets,Bsp%broad
 end if
 !
 ! Select the terminator for the continued fraction.
 term_type=0 !; if (Bsp%hayd_term>0) term_type=2
 call wrtout(std_out,sjoin("Using terminator type: ",itoa(term_type)),"COLL")
 !
 ! Calculate green(w) for the different starting kets.
 green=czero
 do iq=1,nkets
   ABI_MALLOC(ket0,(my_nt))
   ket0 = kets(my_t1:my_t2,iq)
   !
   niter_file=0

   if (can_restart) then
     call haydock_restart(BSp,restart_file,BSE_HAYD_IMEPS,iq,hsize,&
&      niter_file,aa_file,bb_file,phi_np1_file,phi_n_file,comm)
   end if 
   !
   ABI_MALLOC(phi_nm1,(my_nt))
   ABI_MALLOC(phi_n,(my_nt))
   ABI_MALLOC(phi_np1,(my_nt))
   !
   ! TODO: Note the different convention used for the coefficients
   ! Should use the same convention in the Hermitian case.
   niter_max = niter_file + Bsp%niter
   ABI_MALLOC(aa,(niter_max))
   ABI_MALLOC(bb,(niter_max+1))
   ABI_MALLOC(cc,(niter_max+1))
   aa=czero; bb=czero; cc=czero

   if (niter_file==0) then ! Calculation from scratch.
     phi_n   = ket0
     phi_np1 = MATMUL(hreso,ket0) - MATMUL(hcoup,CONJG(ket0)) 
     ket0_hbar_norm = SQRT(two*DBLE(DOT_PRODUCT(phi_n,phi_np1)))  
     phi_n   = phi_n  /ket0_hbar_norm
     phi_np1 = phi_np1/ket0_hbar_norm
     !ket0    = ket0/ket0_hbar_norm
     cc(1)=zero ! <P|F|P>
     !cc(1) =  DOT_PRODUCT(ket0,phi_np1)
     write(std_out,*)" cc(1), ket0_hbar_norm =",cc(1),ket0_hbar_norm  

     phi_nm1 = czero
     niter_done=0  ! TODO Be careful here

   else ! Use the previously calculates a and b.
     niter_done=niter_file
     MSG_ERROR("Restart not coded")
     !aa(1:niter_done) = aa_file
     !bb(1:niter_done) = bb_file
     !phi_np1=phi_np1_file(my_t1:my_t2)   ! Select the slice treated by this node.
     !phi_n  =phi_n_file  (my_t1:my_t2)   
   end if

   if (can_restart) then
     ABI_FREE(aa_file)
     ABI_FREE(bb_file)
     ABI_FREE(cc_file)
     ABI_FREE(phi_np1_file)
     ABI_FREE(phi_n_file)
   end if

   ! This factor gives the correct results
   factor = -nfact*ket0_hbar_norm / SQRT(two)

   ! Which quantity should be checked for convergence?
   check = (/.TRUE.,.TRUE./) 
   if (ABS(Bsp%haydock_tol(2)-one)<tol6) check = (/.TRUE. ,.FALSE./) 
   if (ABS(Bsp%haydock_tol(2)-two)<tol6) check = (/.FALSE.,.TRUE./) 

   call haydock_psherm_optalgo(niter_done,niter_max,BSp%nomega,BSp%omega,BSp%haydock_tol(1),check,hsize,&
&    my_t1,my_t2,hreso,hcoup,factor,term_type,aa,bb,cc,ket0,ket0_hbar_norm,phi_nm1,phi_n,phi_np1,green(:,iq),inn,is_converged,comm)
   !
   ! Save the a"s and the b"s for possible restarting.
   ! 1) Info on the Q.
   ! 2) Number of iterations performed.
   ! 3) do iter=1,niter_performed 
   !      aa(iter),bb(iter)
   !    end do
   ! 4) |n-1>
   !    |n>
   !    |n+1>
   !
   if (my_rank==master) then ! Open the file and writes basic dimensions and info.
     write(out_unt)Bsp%q(:,iq)
     write(out_unt)MIN(inn,niter_max)  ! NB: if the previous loop completed inn=niter_max+1
     do itt=1,MIN(inn,niter_max)        !     if we exited then inn is not incremented by one.
       write(out_unt)itt,aa(itt),bb(itt)
     end do
   end if
   !
   ! cbuff is used as workspace to gather |n-1>, |n> and |n+1>.
   ABI_MALLOC(cbuff,(hsize))
   cbuff=czero; cbuff(my_t1:my_t2) = phi_nm1
   call xmpi_sum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n-1>

   cbuff=czero; cbuff(my_t1:my_t2) = phi_n
   call xmpi_sum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n>

   cbuff=czero; cbuff(my_t1:my_t2) = phi_np1
   call xmpi_sum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n+1>

   ABI_FREE(phi_nm1)
   ABI_FREE(phi_n)
   ABI_FREE(phi_np1)
   ABI_FREE(cbuff)
   ABI_FREE(aa)
   ABI_FREE(bb)
   ABI_FREE(cc)
   ABI_FREE(ket0)
 end do ! iq

 if (my_rank==master) close(out_unt)

 call xmpi_barrier(comm)

end subroutine haydock_psherm
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_psherm_optalgo
!! NAME
!! haydock_psherm_optalgo
!!
!! FUNCTION
!!  Haydock algorithm for pseudo-hermitian matrix
!!
!! INPUTS
!!  niter_done=Number of iterations already performed (0 if the run starts from scratch).
!!  niter_tot=Max number of iterations. Always > niter_done
!!  nomega=Number of Frequency points for the evaluation of the matrix element.
!!  omega(nomega)=Frequency set (imaginary part is already included).
!!  tol_iter=Tollerance used to stop the the algorithm.
!!  check(2)=Logical flags to specify where both the real and the imaginary part of the 
!!    matrix elements of the Green functions have to be checked for convergence. 
!!  hsize=Size of the blocks.
!!  my_t1,my_t2=Indeces of the first and last column stored treated by this done.
!!  term_type=0 if no terminator is used, 1 otherwise.
!!  hreso(hsize,my_t1:my_t2)=The columns of the resonant block.
!!  hcoup(hsize,my_t1:my_t2)=The columns of the coupling block.
!!  factor
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  green(nomega)=Output matrix elements.
!!  inn=Last iteration performed.
!!  is_converged=.TRUE. of the algorithm converged.
!!
!! SIDE EFFECTS
!!  phi_nm1(my_t2-my_t1+1), phi_n(my_t2-my_t1+1)
!!    input: vectors used to initialize the iteration
!!    output: the vectors obtained in the last iteration 
!!  aa(niter_tot) and bb(niter_tot+1)
!!    if niter_done>0: aa(1:niter_done), bb(1:niter_done) store the coefficients of the previous run.
!!    when the routine returns aa(1:inn) and bb(1:inn) contain the matrix elements of the tridiagonal form.
!!  cc(niter_tot+1)
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

subroutine haydock_psherm_optalgo(niter_done,niter_tot,nomega,omega,tol_iter,check,hsize,my_t1,my_t2,hreso,hcoup,&
&  factor,term_type,aa,bb,cc,ket0,ket0_hbar_norm,phi_nm1,phi_n,phi_np1,green,inn,is_converged,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_psherm_optalgo'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter_tot,niter_done,nomega,comm,hsize,my_t1,my_t2,term_type
 integer,intent(out) :: inn
 logical,intent(out) :: is_converged
 real(dp),intent(in) :: tol_iter,ket0_hbar_norm
 complex(dpc),intent(in) :: factor
!arrays
 real(dp),intent(inout) :: bb(niter_tot+1)
 complex(dpc),intent(out) :: green(nomega)
 complex(dpc),intent(in) :: omega(nomega) 
 complex(dpc),intent(inout) :: aa(niter_tot),cc(niter_tot+1)
 complex(dpc),intent(in) :: hreso(hsize,my_t1:my_t2)
 complex(dpc),intent(in) :: hcoup(hsize,my_t1:my_t2)
 complex(dpc),intent(in) :: ket0(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_nm1(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_n  (my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_np1(my_t2-my_t1+1)
 logical,intent(in) :: check(2)

!Local variables ------------------------------
!scalars
 integer :: my_nt,niter_min,nconv,parity,ii,jj,tdim !,ierr
 integer :: row_max,col_max,nlev
 character(len=500) :: msg
 real(dp) :: max_err,mean_err,mean_err2,std_dev,err
 logical :: keep_vectors=.TRUE.
!arrays
 real(dp) :: abs_err(nomega,2) !,ww_err(nomega,2)
 complex(dpc) :: gn0(nomega,niter_tot)
 complex(dpc),allocatable :: oldg(:),newg(:) 
 complex(dpc),allocatable :: hphi_n(:),save_phi(:,:)
 complex(dpc),allocatable ::  alpha(:,:),beta(:,:),ovlp(:,:)
 complex(dpc),allocatable :: phi_test(:),phi_test2(:),g00(:)
 logical :: test(2)

!************************************************************************

 ABI_UNUSED(ket0_hbar_norm)

 my_nt = my_t2-my_t1+1

 ABI_MALLOC(oldg,(nomega))
 ABI_MALLOC(newg,(nomega))
 ABI_MALLOC(g00,(nomega))
 oldg=czero; newg=czero; g00=czero 
 nconv=0

 keep_vectors = (keep_vectors.and.xcomm_size(comm)==1)
 if (keep_vectors) then 
   ABI_MALLOC(save_phi,(my_t2-my_t1+1,niter_tot))
   ABI_CHECK_ALLOC("out of memory in save_phi")        
   save_phi=czero
 end if

 ABI_MALLOC(hphi_n,(hsize))

 do inn=niter_done+1,niter_tot
   !
   ! a(n) = <Vn+1|F|Vn+1> = <Vn|HFH|Vn>) = 0 by symmetry.
   aa(inn)=zero

   ! |n+1> = |n+1> - a(n)|Vn> - a(n)|n-1>
   phi_np1 = phi_np1 - bb(inn)*phi_nm1
   !
   ! |n-1> = |n> 
   ! |n>   = |n+1> 
   phi_nm1 = phi_n
   phi_n   = phi_np1
   !
   !|n+1> = H |n> using resonant eh components.
   parity = (-1)**(inn+1)
   phi_np1 = MATMUL(hreso,phi_n) + parity * MATMUL(hcoup,CONJG(phi_n))
   !call xmpi_sum(hphi_np1,comm,ierr)
   !
   ! B(n+1)= <n|F|n+1>^(1/2) = <n|FH|n>^(1/2))= (2*Re(<n|V+1>))^(1/2) 
   ! by symmetry, where the dot_product is done in the resonant eh sub-space. 
   !
   bb(inn+1)=SQRT(two*DBLE(DOT_PRODUCT(phi_n,phi_np1)))
   !bb(inn+1)=two*DBLE(DOT_PRODUCT(phi_n,phi_np1))
   !call xmpi_sum(bb(inn+1),comm,ierr)
   !bb(inn+1)=SQRT(bb(inn+1)
   !
   !|n+1> =|n+1>/B(n+1)
   phi_n   = phi_n  /bb(inn+1)
   phi_np1 = phi_np1/bb(inn+1)

   if (keep_vectors) save_phi(:,inn) = phi_n

   parity = (-1)**(inn+1) 
   !if (parity==-1) then 
   !  cc(inn+1)=czero
   !else 
     cc(inn+1)=DOT_PRODUCT(ket0,phi_n) + parity * DOT_PRODUCT(phi_n,ket0)
   !end if
   !call xmpi_sum(cc(inn+1),comm,ierr)

   write(msg,'(a,i0,a,3es12.4)')' Iteration number ',inn,', b_i RE(c_i+1) IM(c_i+1) ',bb(inn),REAL(cc(inn+1)),AIMAG(cc(inn+1)) 
   call wrtout(std_out,msg,"COLL")

   call continued_fract(inn,term_type,aa,bb(2:),nomega,omega,g00)
   gn0(:,1) = g00

   if (.FALSE.) then
     gn0(:,2) = (one - omega(:)*g00(:))/bb(2)
     do ii=3,inn
       gn0(:,ii) = -(-bb(ii)*gn0(:,ii-2) -omega(:)*gn0(:,ii-1))/bb(ii+1)
     end do
   else 
     do ii=2,inn
       nlev = inn-ii
       call continued_fract(nlev,term_type,aa,bb(ii+1:),nomega,omega,g00)
       gn0(:,ii) = +bb(ii+1) * g00 * gn0(:,ii-1)
     end do
   end if

   newg=czero
   do ii=1,inn
     newg(:) = newg + cc(ii)* gn0(:,ii)
   end do
   newg = factor*newg
   !
   ! Avoid spurious convergence.
   niter_min=4; if (niter_done>1) niter_min=niter_done+1
   if (inn>niter_min) then
     test=.TRUE.
     abs_err(:,1) = ABS(DBLE (newg-oldg))
     abs_err(:,2) = ABS(AIMAG(newg-oldg))
     !
     if (tol_iter>zero) then 
       ! Test on the L1 norm.
       if (check(1)) test(1) = SUM(abs_err(:,1)) < tol_iter*SUM(ABS(DBLE (newg)))
       if (check(2)) test(2) = SUM(abs_err(:,2)) < tol_iter*SUM(ABS(AIMAG(newg)))
     else
       ! Stringent test for each point.
       if (check(1)) test(1) = ALL( abs_err(:,1) < -tol_iter*ABS(DBLE (newg))) 
       if (check(2)) test(2) = ALL( abs_err(:,2) < -tol_iter*ABS(AIMAG(newg)))
     end if
     !
     if (ALL(test)) then 
       nconv = nconv+1
     else 
       nconv = 0
     end if
     if (nconv==2) then 
       write(msg,'(a,es10.2,a,i0,a)')&
&        " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after ",inn," iterations." 
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       EXIT
     end if
   end if
   !
   oldg = newg
 end do ! inn

 green = newg
 if (nconv/=2) then
   write(msg,'(a,es10.2,a,i0,a)')&
&    " WARNING: Haydock algorithm did not converge within ",tol_iter," after ",niter_tot," iterations."
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 is_converged = (nconv==2)

 ABI_FREE(oldg)
 ABI_FREE(newg)
 ABI_FREE(g00)
 ABI_FREE(hphi_n)

 if (keep_vectors) then
   tdim = MIN(inn,niter_tot)
   ABI_MALLOC(ovlp,(tdim,tdim))

   ABI_MALLOC(phi_test,(hsize))
   ABI_MALLOC(phi_test2,(hsize))

   max_err=smallest_real; mean_err=zero; mean_err2=zero; row_max=-1
   do ii=1,tdim
     parity = (-1)**(ii+1) 
     phi_test  = save_phi(:,ii)
     phi_test2 = MATMUL(hreso,phi_test) + parity * MATMUL(hcoup,CONJG(phi_test))
     ovlp(ii,ii) = DOT_PRODUCT(phi_test,phi_test2) + DOT_PRODUCT(phi_test2,phi_test) 
     err = ABS(ovlp(ii,ii)-cone)
     mean_err  = mean_err + err
     mean_err2 = mean_err2 + err**2
     if (err > max_err) then
       max_err = err 
       row_max = ii
     end if
   end do
   mean_err = mean_err/tdim
   std_dev = mean_err2/tdim -mean_err**2
   write(std_out,'(a,i0,1x,3es14.6)')&
&   " Error in normalization (ii, max_err,mean,std_dev): ",row_max,max_err,mean_err,std_dev

   ABI_FREE(phi_test)
   ABI_FREE(phi_test2)
                                             
   ABI_MALLOC(alpha,(hsize,tdim))
   alpha = MATMUL(hreso,save_phi(:,1:tdim))

   do ii=1,tdim
     parity = (-1)**(ii+1) 
     alpha(:,ii) =  alpha(:,ii) + parity*MATMUL(hcoup,CONJG(save_phi(:,ii))) 
   end do

   ovlp = MATMUL(TRANSPOSE(CONJG(save_phi(:,1:tdim))),alpha)

   ABI_MALLOC(beta,(hsize,tdim))
   do ii=1,tdim
     parity = (-1)**(ii+1) 
     beta(:,ii)  =  parity*save_phi(:,ii)
     alpha(:,ii) = -parity*alpha(:,ii)
   end do

   ovlp = ovlp - MATMUL(TRANSPOSE(CONJG(beta)),alpha)

   max_err=smallest_real; row_max=-1; col_max=-1
   mean_err=zero; mean_err2=zero
   do jj=1,tdim
     do ii=1,jj
       err = ABS(ovlp(ii,jj))
       if (ii==jj) err = ABS(err - one)
       mean_err  = mean_err + err
       mean_err2 = mean_err2 + err**2
       if (err > max_err) then
         max_err = err
         row_max=ii
         col_max=jj
       end if
     end do
   end do

   mean_err = mean_err/(tdim*(tdim+1)/2)
   std_dev = mean_err2/(tdim*(tdim+1)/2) - mean_err**2
   write(std_out,'(a,2(i0,1x),3es14.6)')&
&     " Error in Hbar-ortho (i,j), max_err, mean, std_dev ",row_max,col_max,max_err,mean_err,std_dev
   !call print_arr(ovlp,max_r=185,max_c=10,unit=std_out)

   ABI_FREE(alpha)
   ABI_FREE(beta)
   ABI_FREE(ovlp)
   ABI_FREE(save_phi)
 end if

end subroutine haydock_psherm_optalgo
!!***

!----------------------------------------------------------------------

end module m_haydock
!!***
