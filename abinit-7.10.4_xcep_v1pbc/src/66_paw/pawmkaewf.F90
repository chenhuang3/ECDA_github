!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmkaewf
!! NAME
!! pawmkaewf
!!
!! FUNCTION
!! Construct complete AE wave functions on the fine FFT grid adding onsite PAW corrections.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! dimcprj(natom)=array of dimensions of array cprj (not ordered)
!! mband=maximum number of bands
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mkmem=number of k points treated by this node.
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! mpi_comm_atom=--optional-- MPI communicator over atoms
!! mpw=maximum dimensioned size of npw.
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in cell
!! ntypat=number of types of atoms in the cell
!! nkpt=Total number of k-points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! unks=unit number for G vectors.
!! nband(nkpt*nsppol)=Number of bands for each k-point and spin.
!! istwfk(nkpt)=Storage mode at each k-point.
!! paral_kgb=Option for kgb parallelism
!! Pawfgrtab(natom) <type(pawfgrtab_type)> : data about the fine grid around each atom
!! Pawrad(ntypat) <type(pawrad_type)> : radial mesh data for each type of atom
!! Pawtab(ntypat) <type(pawtab_type)> : PAW functions around each type of atom
!! Psps <type(pseudopotential_type)> : basic pseudopotential data
!! Dtfil <type(datafiles_type)>=variables related to files
!! cg(2,mcg)=planewave coefficients of wavefunctions.
!! Cprj(natom,nspinor*mband*mkmem*nsppol)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!   and each |p_lmn> non-local projector
!! npwarr(nkpt)=Number of plane waves at each k-point
!! ngfftf(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  Note that ngfftf refers to the fine mesh.
!! kg(3,mpw*mkmem)=reduced planewave coordinates
!! eigen(mband*nkpt*nsppol)=eigenvalues (hartree) for all bands at each k point
!! occ(mband*nkpt*nsppol)=occupations for all bands at each k point
!! Hdr<hdr_type>=the header of wf, den and pot files
!! kpt(3,nkpt)=reduced coordinates of k points.
!!
!! OUTPUT
!!  ierr=Status error
!!  Main output is written on file (NETCDF file format).
!!
!! NOTES
!! In PAW calculations, the pseudized wavefunction us represented
!! on a relatively small plane wave basis set and is not normalized
!! as it does not include the on-site PAW contributions which is described
!! in terms of real spherical harmonics and radial functions.
!! For post-processing and proper visualization, it is necessary
!! to use the full electronic wave function, which is what this subroutine constructs.
!! Specifically, it computes the pseudo part by doing an FFT from G- to r-space
!! using the dense mesh defined by pawecutdg. The on-site PAW terms are also
!! computed in real space inside each sphere and added to the pseudo part.
!! Notice that this formula is expressed on the fine grid, and requires
!! interpolating the PAW radial functions onto this grid, as well as calling
!! initylmr in order to get the angular functions on the grid points.
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      abi_etsf_dims_init,abi_etsf_electrons_put,abi_etsf_geo_put
!!      etsf_io_data_init,etsf_io_electrons_put,etsf_io_file_merge
!!      etsf_io_low_close,etsf_io_low_open_modify,etsf_io_low_set_write_mode
!!      etsf_io_main_def,etsf_io_main_put,flush_unit,fourwf,free_my_atmtab
!!      get_my_atmtab,hdr_io_etsf,ini_wf_etsf,int2char4,nhatgrid
!!      paw_pwaves_lmn_free,paw_pwaves_lmn_init,pawcprj_alloc,pawcprj_destroy
!!      pawcprj_diskinit_r,pawfgrtab_destroy,pawfgrtab_init,pawfgrtab_print
!!      sphereboundary,wrap2_zero_one,wrtout,xcart2xred,xmpi_barrier,xmpi_max
!!      xmpi_sum,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawmkaewf(Dtset,my_natom,natom,mpw,mband,mcg,mcprj,nkpt,mkmem,nsppol,ntypat,nband,istwfk,npwarr,kpt,&
& paral_kgb,ngfftf,kg,dimcprj,Pawfgrtab,Pawrad,Pawtab,gmet,rprimd,ucvol,&
& Psps,Hdr,Dtfil,eigen,occ,cg,Cprj,MPI_enreg,ierr,pseudo_norms,set_k,set_band , &
& mpi_atmtab,mpi_comm_atom) ! Optional arguments


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_xmpi
 use m_blas
 use m_splines
 use m_wffile
 use m_errors

 use m_io_tools,       only : get_unit, flush_unit
 use m_fstrings,       only : int2char4
 use m_numeric_tools,  only : imax_loc, wrap2_zero_one
 use m_pptools,        only : printxsf
 use m_header,         only : hdr_io_etsf, hdr_skip, hdr_io
 use m_crystal,        only : crystal_free, crystal_t
 use m_crystal_io,     only : crystal_from_hdr
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgrtab,      only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_destroy, pawfgrtab_print
 use m_pawcprj,        only : pawcprj_type, pawcprj_diskinit_r, pawcprj_alloc, pawcprj_get, pawcprj_destroy
 use m_paw_pwaves_lmn, only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_paral_atom,     only : get_my_atmtab, free_my_atmtab

#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
 use etsf_io_file,     only : etsf_io_file_merge
 use m_abi_etsf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmkaewf'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_66_paw, except_this_one => pawmkaewf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,ntypat,mband,mcg,mcprj,mkmem,mpw,nsppol,paral_kgb,nkpt
 integer,intent(in),optional :: mpi_comm_atom,set_k,set_band
 integer,intent(out) :: ierr
 type(Datafiles_type),intent(in) :: Dtfil
 type(pseudopotential_type),intent(in) :: Psps
 type(MPI_type),intent(inout) :: MPI_enreg
 type(hdr_type),intent(inout) :: Hdr
 type(dataset_type),intent(in) :: Dtset
!arrays
 integer,intent(in) :: nband(nkpt*nsppol),istwfk(nkpt),npwarr(nkpt),dimcprj(natom)
 integer,intent(in) :: ngfftf(18),kg(3,mpw*mkmem)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),target,intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),target,intent(in) :: occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: kpt(3,nkpt)
 real(dp),intent(in) :: gmet(3,3),rprimd(3,3),ucvol
 real(dp),optional,intent(out) :: pseudo_norms(nsppol,nkpt,mband)
 type(pawfgrtab_type),intent(in) :: Pawfgrtab(my_natom)
 type(pawrad_type),intent(in) :: Pawrad(ntypat)
 type(pawtab_type),intent(in) :: Pawtab(ntypat)
 type(pawcprj_type),intent(in) :: Cprj(natom,mcprj)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf0=0,tim_rwwf0=0
 integer :: bdtot_index,iband,icg,mgfftf
 integer :: iatom,iatom_tot,ifgd,ifftsph,ifft,itypat,ispinor,ipw,ndat,ii,i1,i2,i3
 integer :: jl,jm,jlmn
 integer :: max_nfgd,nfgd,ln_size,lmn_size,my_comm_atom,option
 integer :: iorder_cprj,spaceComm,rank,ibsp,ibg,isppol,ikpt,nband_k,cplex
 integer :: n1,n2,n3,n4,n5,n6,ikg,npwout,istwf_k,npw_k
 integer :: indx,nfftot,my_spin,nprocs,tmp_unt
 integer :: optcut,optgr0,optgr1,optgr2,optrad,start_band,start_kpt,stop_kpt,stop_band
 logical :: my_atmtab_allocated,paral_atom
 real(dp),parameter :: weight1=one
 real(dp) :: phj,tphj,re_p,im_p,norm,norm_rerr,max_rerr,imur,reur,arg
 character(len=500) :: msg
 !character(len=fnlen) :: my_basename , xsf_fname
 character(len=fnlen) :: my_basenames(4)
!arrays
 integer,allocatable :: l_size_atm(:)
 integer,pointer :: my_atmtab(:)
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 integer,allocatable :: atindx(:),atindx1(:),nattyp(:)
 integer,allocatable,target :: my_kpoints(:),my_spins(:)
 integer,allocatable :: my_kstable(:,:),my_nkpt(:)
 real(dp),allocatable :: r0shift(:,:,:),phk_atm(:,:,:)
 real(dp) :: red(3),shift(3)
 real(dp) :: rfft(3)
 real(dp) :: kpoint(3),cp_fact(2)
 real(dp),allocatable :: buf_tmp(:,:,:),fofgin(:,:),fofgin_down(:,:),fofgout(:,:)
 real(dp),allocatable :: denpot(:,:,:)
 real(dp),allocatable :: fofr(:,:,:,:),fofr_down(:,:,:,:),phkr(:,:)
 real(dp),allocatable,target :: ur(:,:), ur_pw(:,:)
 real(dp),allocatable,target :: ur_ae_onsite(:,:),ur_ps_onsite(:,:)
 real(dp),allocatable :: ur_mask(:),xcart(:,:),dummy_1d(:)
 real(dp),allocatable :: rsph_red(:,:)
 real(dp),pointer :: rsph_cart(:,:)
 type(pawcprj_type),allocatable :: Cprj_k(:,:)
 type(pawfgrtab_type) :: local_pawfgrtab(my_natom)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

#if defined HAVE_TRIO_ETSF_IO
 integer :: fform,rdwr,var_main,irank
 integer :: ncids(4)
 logical :: kdep,lstat
 character(len=10) :: tag
 !character(len=fnlen) :: out_file !,fname_part
 character(len=fnlen) :: my_fnames(4), out_files(4)
 character(len=256),allocatable :: merge_files(:,:)
 type(etsf_dims) :: Dims
 type(etsf_groups_flags) :: Groups_flags
 !type(etsf_groups) :: Groups
 type(etsf_split) :: Split
 type(etsf_main) :: MainFolder
 type(etsf_io_low_error) :: Error_data
 type(etsf_electrons) :: Electrons_folder
 type(Wvl_wf_type) :: Dummy_wfs
#endif

! ************************************************************************

 DBG_ENTER("COLL")

!Init parallelism
 spaceComm=MPI_enreg%comm_cell
 nprocs=xcomm_size(spaceComm)
 rank=MPI_enreg%me_kpt

!Compatibility tests
 ABI_CHECK(paral_kgb==0,"paral_kgb/=0 not coded")
 ABI_CHECK(MPI_enreg%paral_kgb/=1,"mode_para=b not coded")
 ABI_CHECK(SIZE(dimcprj)>0,"dimcprj should be allocated")
 ABI_CHECK(mpi_enreg%paral_spinor==0,"parallelisation over spinors not implemented")
 ABI_CHECK(nprocs==1,"k spin parallelism not yet active")

 ABI_ALLOCATE(atindx,(natom))
 ABI_ALLOCATE(atindx1,(natom))
 ABI_ALLOCATE(nattyp,(ntypat))

!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_self;if (present(mpi_comm_atom)) my_comm_atom=mpi_comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)
 ABI_ALLOCATE(l_size_atm,(my_natom))

 indx=1
 do itypat=1,ntypat
   nattyp(itypat)=0
   do iatom=1,natom
     if (Dtset%typat(iatom)==itypat) then
       atindx (iatom )=indx
       atindx1(indx  )=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

 ABI_ALLOCATE(xcart,(3,Dtset%natom))
 call xred2xcart(Dtset%natom,Hdr%rprimd,xcart,Hdr%xred)

!If collection of pseudo norms is enabled, make sure the array is initialised
 if (present(pseudo_norms)) pseudo_norms = zero

!use a local copy of pawfgrtab to make sure we use the correction in the paw spheres
!the usual pawfgrtab uses r_shape which may not be the same as r_paw
 if (my_natom>0) then
   do iatom = 1, my_natom
     itypat=pawfgrtab(iatom)%itypat
     l_size_atm(iatom) = pawtab(itypat)%lcut_size
   end do
   if (paral_atom) then
     call pawfgrtab_init(local_pawfgrtab,Pawfgrtab(1)%cplex,l_size_atm,Dtset%nspden,Dtset%typat,&
&     mpi_atmtab=my_atmtab,mpi_comm_atom=my_comm_atom)
   else
     call pawfgrtab_init(local_pawfgrtab,Pawfgrtab(1)%cplex,l_size_atm,Dtset%nspden,Dtset%typat)
   end if
 end if
 optcut = 1 ! use rpaw to construct local_pawfgrtab
 optgr0 = 0; optgr1 = 0; optgr2 = 0 ! dont need gY terms locally
 optrad = 1 ! do store r-R

 if (paral_atom) then
   call nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfftf,ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,local_pawfgrtab,pawtab,rprimd,Dtset%typat,ucvol,Hdr%xred,&
&   mpi_comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfftf,ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,local_pawfgrtab,pawtab,rprimd,Dtset%typat,ucvol,Hdr%xred)
 end if
!now local_pawfgrtab is ready to use

 max_nfgd=MAXVAL(local_pawfgrtab(:)%nfgd) ! MAX no. of points in the fine grid for this PAW sphere
 ABI_ALLOCATE(r0shift,(3,max_nfgd,my_natom))
 ABI_ALLOCATE(phk_atm,(2,max_nfgd,my_natom))

 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   nfgd=local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere
   ABI_ALLOCATE(rsph_red,(3,nfgd))
   ABI_ALLOCATE(rsph_cart,(3,nfgd))
!  rsph_cart => pawfgrtab(iatom)%rfgd
   do ifgd=1,nfgd
     rsph_cart(:,ifgd) = local_pawfgrtab(iatom)%rfgd(:,ifgd) + xcart(:,iatom_tot)
   end do
   call xcart2xred(nfgd,rprimd,rsph_cart,rsph_red) ! we work in reduced coordinates.
   do ifgd=1,nfgd
     call wrap2_zero_one(rsph_red(1,ifgd),red(1),shift(1)) ! num = red + shift
     call wrap2_zero_one(rsph_red(2,ifgd),red(2),shift(2))
     call wrap2_zero_one(rsph_red(3,ifgd),red(3),shift(3))
     r0shift(:,ifgd,iatom) = shift
     if (ANY( ABS(shift) > tol12)) then
!      MSG_WARNING("rmR_red is outside the first unit cell.")
!      write(std_out,*)rsph_red(:,ifgd),shift
     end if
   end do
   ABI_DEALLOCATE(rsph_red)
   ABI_DEALLOCATE(rsph_cart)
 end do

 if ((.not.paral_atom).and.my_natom>0) then
   call pawfgrtab_print(local_pawfgrtab,natom=natom,unit=std_out,&
&   prtvol=Dtset%prtvol,mode_paral="COLL")
 end if

 ierr=0
#ifndef HAVE_TRIO_ETSF_IO
 ierr=-1
 write(msg,'(3a)')&
& " ETSF-IO support must be enabled in order to output AE PAW wavefunction. ",ch10,&
& " No output will be produced, use --enable-etsf-io at configure-time. "
 MSG_WARNING(msg)
 RETURN
!These statements are necessary to avoid the compiler complain about unused variables:
 ii=Dtset%usepaw;ii=Dtfil%unpaw;ii=Hdr%usepaw
#endif


 if (xmpi_paral==1)then
   call wrtout(std_out,'pawmkaewf: loop on k-points and spins done in parallel','COLL')
   call xmpi_barrier(spaceComm)
 end if

 mgfftf=MAXVAL(ngfftf(1:3))


!=== Calculate my list of k-points and spin ===
!* my_kstable gives the sequential index for each k-point treated by rank.
!* cannot check for MPI_enreg%proc_distrb if nprocs ==1
 ABI_ALLOCATE(my_kstable,(nkpt,nsppol))
 my_kstable=0

 if (nprocs==1) then

   ii=0
   do ikpt=1,nkpt
     ii=ii+1; my_kstable(ikpt,:) = ii
   end do

   ABI_ALLOCATE(my_spins,(nsppol))
   do isppol=1,nsppol
     my_spins(isppol)=isppol
   end do

   ABI_ALLOCATE(my_kpoints,(nkpt))
   do ikpt=1,nkpt
     my_kpoints(ikpt) = ikpt
   end do

 else ! parallelism over k and spin.

   do isppol=1,nsppol
     ii=0
     do ikpt=1,nkpt
       nband_k = nband(ikpt+(isppol-1)*nkpt)
       if (ALL(MPI_enreg%proc_distrb(ikpt,1:nband_k,isppol)==rank)) then
         ii=ii+1
         my_kstable(ikpt,isppol)=ii
       end if
     end do
   end do

   ABI_ALLOCATE(my_nkpt,(nsppol))
   do isppol=1,nsppol
     my_nkpt(isppol) = COUNT(my_kstable(:,isppol)>0)
   end do

!  Each node has to deal with a single spin.
   if (nsppol>1 .and. ALL(my_nkpt>0)) then
     msg =' Non optimal distribution, some wave functions won''t be correctly initialized.'
     MSG_ERROR(msg)
   end if
   my_spin = imax_loc(my_nkpt)

   ABI_ALLOCATE(my_spins,(1))
   my_spins(1)=my_spin

   ABI_ALLOCATE(my_kpoints,(my_nkpt(my_spin)))
   ii=0
   do ikpt=1,nkpt
     if (my_kstable(ikpt,my_spin)/=0) then
       ii=ii+1
       my_kpoints(ii) = ikpt
     end if
   end do

   ABI_DEALLOCATE(my_nkpt)

 end if ! nprocs==1


#if defined HAVE_TRIO_ETSF_IO

!=== Initialize NETCDF files ===
 my_basenames(1)=dtfil%fnameabo_ae_wfk
 my_basenames(2)=dtfil%fnameabo_ps_wfk
 my_basenames(3)=dtfil%fnameabo_ae_onsite_wfk
 my_basenames(4)=dtfil%fnameabo_ps_onsite_wfk

!* For parallel case: the index of the processor must be appended.
!XG 100108 : One would better have it done in dtfil_init1 !!!!
 if (nprocs>1) then
   ABI_ALLOCATE(merge_files,(4, nprocs))
   out_files(1)=TRIM(dtfil%fnameabo_ae_wfk)//"-etsf.nc"
   out_files(2)=TRIM(dtfil%fnameabo_ps_wfk)//"-etsf.nc"
   out_files(3)=TRIM(dtfil%fnameabo_ae_onsite_wfk)//"-etsf.nc"
   out_files(4)=TRIM(dtfil%fnameabo_ps_onsite_wfk)//"-etsf.nc"
   do ii = 1, 4
     do irank=1,nprocs
       call int2char4(irank,tag)
       ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
       !merge_files(ii,irank)=TRIM(Dtfil%filnam_ds(4))//'_P-'//trim(tag)//"_AE_WFK-etsf.nc"
       merge_files(ii,irank)=TRIM(my_basenames(ii))//'_P-'//trim(tag)//"-etsf.nc"
       if (irank==rank+1) then
         my_basenames(ii)=TRIM(my_basenames(ii))//'_P-'//trim(tag)
         !my_basenames(1)=TRIM(Dtfil%filnam_ds(4))//'_P-'//trim(tag)//"_AE_WFK"
       end if
     end do
   end do
!  my_basenames(1)=TRIM(Dtfil%filnam_ds(4))//'_P-'//trim(tag)//"_AE_WFK"
!  my_basenames(1)=merge_files(rank+1)
 end  if

 do ii = 1, 4
   my_fnames(ii)=TRIM(my_basenames(ii))//"-etsf.nc"
 end do

 write(msg,'(2a)')' Opening file for AE PAW wave functions: ',TRIM(my_fnames(1))
 call wrtout(std_out,msg,'PERS')
 call wrtout(ab_out,msg,'PERS')
 call xmpi_barrier(spaceComm)


!Initialize Dims, remeber the hacking with with Dims
!call abi_etsf_dims_init(Dtset,2,Hdr%lmn_size,Psps,Dummy_wfs)
 call abi_etsf_dims_init(Dims,Dtset,2,Psps,Dummy_wfs)

!Change some values since we work in real space on the dense FFT mesh.
 Dims%max_number_of_coefficients      = etsf_no_dimension
 Dims%max_number_of_basis_grid_points = etsf_no_dimension
 Dims%number_of_localization_regions  = etsf_no_dimension
 Dims%real_or_complex_coefficients    = etsf_no_dimension

 Dims%real_or_complex_wavefunctions   = 2

 Dims%number_of_grid_points_vector1  = ngfftf(1)
 Dims%number_of_grid_points_vector2  = ngfftf(2)
 Dims%number_of_grid_points_vector3  = ngfftf(3)

!Dimensions for variables that can be splitted.
 Dims%my_number_of_kpoints = SIZE(my_kpoints) !etsf_no_dimension
 Dims%my_number_of_spins   = SIZE(my_spins)   !etsf_no_dimension

!Split data using k-points and spins
 if (nprocs>1) then
   !BEGIN DEBUG
   write(msg,'(a,i4)') 'number of kpoints:', SIZE(my_kpoints)
   call wrtout(std_out, msg,'PERS')
   write(msg, *) my_kpoints(:)
   call wrtout(std_out, msg,'PERS')
   !END DEBUG
   nullify(Split%my_kpoints)
   Split%my_kpoints => my_kpoints(:)
   nullify(Split%my_spins)
   if (nsppol>1) Split%my_spins => my_spins(:)
 else
   nullify(Split%my_kpoints)
   nullify(Split%my_spins)
 end if

!=== Set-up the variables ===
!* These mandatory values are always written by the hdr_io_etsf() routine.
 Groups_flags%geometry  = etsf_geometry_all
 Groups_flags%electrons = etsf_electrons_all - etsf_electrons_x_functional - etsf_electrons_c_functional
 Groups_flags%kpoints   = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights

!These variables may be written depending on prt<something> input variables.
 Groups_flags%basisdata = etsf_basisdata_basis_set
 if (Dtset%usewvl==0) then
   Groups_flags%basisdata= Groups_flags%basisdata + &
&   etsf_basisdata_kin_cutoff + etsf_basisdata_n_coeff + etsf_basisdata_red_coord_pw
 else
   Groups_flags%basisdata= Groups_flags%basisdata + etsf_basisdata_coord_grid + etsf_basisdata_n_coeff_grid
 end if

!Groups_flags%basisdata = etsf_basisdata_none

!Wavefunctions in real space.
 Groups_flags%main = etsf_main_wfs_rsp

!=== Create the file ===
!* If the group contains main, we remove it for a while to be sure to
!add it at the end, after ABINIT private variables.
 var_main = Groups_flags%main
 Groups_flags%main = etsf_main_none

 write(std_out,*)"Before  etsf_io_data_init"
 kdep=.TRUE.

 call etsf_io_data_init(my_fnames(1),Groups_flags,Dims,&
& 'PAW AE wavefunction given in real space',    &
& 'PAW AE Wavefunction File generated by ABINIT with ETSF_IO',&
& lstat,Error_data,k_dependent=kdep,overwrite=.TRUE.,Split_definition=Split)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_data_init(my_fnames(2),Groups_flags,Dims,&
& 'PAW PS wavefunction given in real space',    &
& 'PAW PS Wavefunction File generated by ABINIT with ETSF_IO',&
& lstat,Error_data,k_dependent=kdep,overwrite=.TRUE.,Split_definition=Split)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_data_init(my_fnames(3),Groups_flags,Dims,&
& 'PAW AE Onsite wavefunction given in real space',    &
& 'PAW AE Onsite Wavefunction File generated by ABINIT with ETSF_IO',&
& lstat,Error_data,k_dependent=kdep,overwrite=.TRUE.,Split_definition=Split)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_data_init(my_fnames(4),Groups_flags,Dims,&
& 'PAW PS Onsite wavefunction given in real space',    &
& 'PAW PS Onsite Wavefunction File generated by ABINIT with ETSF_IO',&
& lstat,Error_data,k_dependent=kdep,overwrite=.TRUE.,Split_definition=Split)
 ETSF_CHECK_ERROR(lstat,Error_data)

 write(std_out,*)" my_number_of_kpoints ",Dims%my_number_of_kpoints
 write(std_out,*)" my_number_of_spins ",Dims%my_number_of_spins


 do ii = 1, 4

  !* Add the private ABINIT variables.
   call etsf_io_low_open_modify(ncids(ii),my_fnames(ii),lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
   
   call ini_wf_etsf(ncids(ii),Dtset%usewvl,Hdr%lmn_size,Psps%npsp,Psps%ntypat)
   
  !Add the main part as last variables in the ETSF file.
   write(std_out,*)"Before  etsf_io_main_def"
   call etsf_io_main_def(ncids(ii),lstat,Error_data,k_dependent=kdep,flags=var_main, Split=Split)
   ETSF_CHECK_ERROR(lstat,Error_data)
   
  !Close the file.
   call etsf_io_low_close(ncids(ii),lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
   
  !Complete the geometry information with missing values from hdr_io().
   call abi_etsf_geo_put(Dtset,my_basenames(ii),Psps)
   
  !To use the following statements, do not forget to declare:
  !timrev(integer), Cryst(crystal_t)
  !timrev=2
  !call crystal_from_hdr(Cryst,Hdr,timrev)
  !call abi_crystal_put(Cryst,my_fnames(ii))
  !call destroy_crystal(Cryst)
   
   call abi_etsf_electrons_put(Dtset,my_basenames(ii))
   
  !We open again for further additions
   call etsf_io_low_open_modify(ncids(ii),my_fnames(ii),lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
   
  !* Write the header.
  !FIXME problem in the case of splitting over k-points due to SIZE mismatch
  !in hdr%npwarr(number_of_kpoints) and number_of_coefficients(my_mkpt)
   
   fform=502; rdwr=2
   call hdr_io_etsf(fform,Hdr,rdwr,ncids(ii))
   
   call etsf_io_low_close(ncids(ii),lstat,Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
   
  !=== Prepare the writing of the results ===
  !
  !1) Open file for writing
   call etsf_io_low_open_modify(ncids(ii),my_fnames(ii),lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
   
  !2) Switch to write mode.
   call etsf_io_low_set_write_mode(ncids(ii),lstat,Error_data=Error_data)
   ETSF_CHECK_ERROR(lstat,Error_data)
 end do
#endif

!Init structure storing phi_{nlm} and tphi_(nlm} on the dense FFT points located in the PAW spheres.
 ABI_DATATYPE_ALLOCATE(Paw_onsite,(natom))
!DEBUG
 write(std_out,*)' pawmkaewf : before call init_paw_pwaves_lmn, rprimd=',rprimd
!ENDDEBUG
 if (paral_atom) then
   call paw_pwaves_lmn_init(Paw_onsite,my_natom,natom,ntypat,rprimd,xcart,Pawtab,Pawrad,local_pawfgrtab, &
&   mpi_comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call paw_pwaves_lmn_init(Paw_onsite,my_natom,natom,ntypat,rprimd,xcart,Pawtab,Pawrad,local_pawfgrtab)
 end if

!FIXME check ordering in cprj and Eventually in external file
!why is iorder_cprj not stored in the file for crosschecking purpose?
!Here Im assuming cprj are not ordered!

 iorder_cprj=0
 call pawcprj_diskinit_r(atindx1,natom,iorder_cprj,mkmem,natom,0,dimcprj,dtset%nspinor,Dtfil%unpaw)

!n4,n5,n6 are FFT dimensions, modified to avoid cache trashing
 n1=ngfftf(1); n2=ngfftf(2); n3=ngfftf(3)
 n4=ngfftf(4); n5=ngfftf(5); n6=ngfftf(6)
 nfftot=PRODUCT(ngfftf(1:3))

 !my_basename = trim(Dtfil%fnameabo_pawwf)

!=== Loop over spin and k points ===
 bdtot_index=0; icg=0; ibg=0; norm_rerr=smallest_real
 do isppol=1,nsppol

   ikg=0
   start_kpt=1
   stop_kpt=nkpt
!  Check if k-point was specified (only serial)
   if (present(set_k).AND.nprocs==1) then
     if (set_k/=0) then
       start_kpt = set_k
       stop_kpt = set_k
     end if
   end if

   do ikpt=start_kpt,stop_kpt

     kpoint  = kpt(:,ikpt)
     nband_k = nband(ikpt+(isppol-1)*nkpt)
     npw_k   = npwarr(ikpt)
     istwf_k = istwfk(ikpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,rank)) then
       bdtot_index=bdtot_index+nband_k; CYCLE
     end if

     ABI_ALLOCATE(phkr,(2,nfftot))
     do i3=0,n3-1
       rfft(3)=DBLE(i3)/n3
       do i2=0,n2-1
         rfft(2)=DBLE(i2)/n2
         do i1=0,n1-1
           rfft(1)=DBLE(i1)/n1
           ifft = 1 +i1 +i2*n1 +i3*n1*n2
           phkr(1,ifft) = COS(two_pi*dot_product(kpoint,rfft))
           phkr(2,ifft) = SIN(two_pi*dot_product(kpoint,rfft))
         end do
       end do
     end do
!    phkr(1,:)=one
!    phkr(2,:)=zero

!    Calculate the phase for the onsite PAW contributions.
     do iatom=1,my_natom
       nfgd=local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere
       do ifgd=1,nfgd
         arg = -two_pi* dot_product(r0shift(:,ifgd,iatom),kpoint)
         phk_atm(1,ifgd,iatom) = COS(arg)
         phk_atm(2,ifgd,iatom) = SIN(arg)
       end do
     end do

     ABI_DATATYPE_ALLOCATE(Cprj_k,(natom,dtset%nspinor*nband_k))
     call pawcprj_alloc(Cprj_k,0,dimcprj)

!    Extract cprj for this k-point.
     ibsp=0
     do iband=1,nband_k
       do ispinor=1,dtset%nspinor
         ibsp=ibsp+1
         do iatom=1,natom
           Cprj_k(iatom,ibsp)%cp(:,:)=Cprj(iatom,ibsp+ibg)%cp(:,:)
         end do
       end do
     end do

     ABI_ALLOCATE(gbound,(2*mgfftf+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))

!    Do i/o as needed
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     call sphereboundary(gbound,istwf_k,kg_k,mgfftf,npw_k)

#if defined HAVE_TRIO_ETSF_IO
!    === Write eigenvalues and occupations ===
!    write(std_out,*)eigen(bdtot_index+1:bdtot_index+nband_k)
     Electrons_folder%eigenvalues__number_of_states = Dtset%mband
     Electrons_folder%eigenvalues%data1D         => eigen(bdtot_index+1:bdtot_index+nband_k)
     Electrons_folder%eigenvalues__kpoint_access = my_kstable(ikpt,isppol) !ikpt
     Electrons_folder%eigenvalues__spin_access   = isppol
     if (nprocs>1) Electrons_folder%eigenvalues__spin_access = 1

     Electrons_folder%occupations__number_of_states = Dtset%mband
     Electrons_folder%occupations%data1D         => occ(bdtot_index+1:bdtot_index+nband_k)
     Electrons_folder%occupations__kpoint_access = my_kstable(ikpt,isppol) !ikpt
     Electrons_folder%occupations__spin_access   = isppol
     if (nprocs>1) Electrons_folder%occupations__spin_access = 1

     write(std_out,*)"rank ",rank," about to write",my_kstable(ikpt,isppol)

     call etsf_io_electrons_put(ncids(1),Electrons_folder,lstat,Error_data)
     ETSF_CHECK_ERROR(lstat,Error_data)
#endif

     start_band = 1
     stop_band = nband_k
!    If a single band is requested, neuter the loop (only serial)
     if (present(set_band).AND.nprocs==1) then
       if (set_band/=0) then
         start_band = set_band
         stop_band = set_band
       end if
     end if
     do iband=start_band,stop_band ! Loop over bands.

!      * Fourier transform on the real fft box of the smooth part.
       ndat=Dtset%nspinor
       ABI_ALLOCATE(fofgin,(2,npw_k*ndat))
       ABI_ALLOCATE(fofr,(2,n4,n5,n6*ndat))

       do ipw=1,npw_k*dtset%nspinor
         fofgin(:,ipw)=cg(:,ipw+(iband-1)*npw_k*dtset%nspinor+icg)
       end do

!      Complex can be set to 0 with this option(0) of fourwf
       option=0; cplex=0; npwout=1
       ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
       ABI_ALLOCATE(fofgout,(2,npwout*ndat))

       call fourwf(cplex,denpot,fofgin(:,1:npw_k),fofgout,fofr(:,:,:,1:n6),gbound,gbound,istwf_k,kg_k,kg_k,&
&       mgfftf,MPI_enreg,1,ngfftf,npw_k,npwout,n4,n5,n6,option,paral_kgb,tim_fourwf0,weight1,weight1,&
&       use_gpu_cuda=Dtset%use_gpu_cuda)

!      Here I do not know if fourwf works in the case of spinors,
!      It seems that not all fftalg option support ndata! should check!
!      Do not forget to declare real(dp)::fofgin_down(:,:) to use the following statements
       if (Dtset%nspinor==2) then
         ABI_ALLOCATE(fofgin_down,(2,npw_k))
         ABI_ALLOCATE(fofr_down,(2,n4,n5,n6))
         fofgin_down(:,:)=fofgin(:,1+npw_k:2*npw_k)
!        Complex can be set to 0 with this option(0) of fourwf
!        cplex=1; option=1; npwout=1; ndat=1
!        NOTE: fofr_down can NOT be replaced by fofr(:,:,:,n6+1:2*n6), or else
!        the data in fofr(:,:,:,1:n6) will be the same with fofr(:,:,:,n6+1:2*n6)
         call fourwf(cplex,denpot,fofgin_down,fofgout,fofr_down,gbound,gbound,istwf_k,kg_k,kg_k,&
&         mgfftf,MPI_enreg,1,ngfftf,npw_k,npwout,n4,n5,n6,option,paral_kgb,tim_fourwf0,weight1,weight1)
         ABI_DEALLOCATE(fofgin_down)
       end if

       ABI_ALLOCATE(ur,(2,n1*n2*n3*ndat))
       ABI_ALLOCATE(ur_ae_onsite,(2,n1*n2*n3))
       ABI_ALLOCATE(ur_ps_onsite,(2,n1*n2*n3))
       ABI_ALLOCATE(ur_pw,(2,n1*n2*n3*ndat))
       ABI_ALLOCATE(ur_mask,(n1*n2*n3))

       ur=zero;ur_ae_onsite=zero;ur_ps_onsite=zero;ur_pw=zero;ur_mask=zero

!      * Add phase e^{ikr} since it is contained in cprj.
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             ii = i1 + n1*(i2-1)+ n1*n2*(i3-1)
             ur_pw(:,ii)=fofr(:,i1,i2,i3) ! Save pw part separately without the phase.
             ur(1,ii)= fofr(1,i1,i2,i3) * phkr(1,ii) - fofr(2,i1,i2,i3) * phkr(2,ii)
             ur(2,ii)= fofr(1,i1,i2,i3) * phkr(2,ii) + fofr(2,i1,i2,i3) * phkr(1,ii)
             if(Dtset%nspinor==2) then
               ur_pw(:,ii+n1*n2*n3)=fofr_down(:,i1,i2,i3) ! Save pw part separately without the phase.
               ur(1,ii+n1*n2*n3)= fofr_down(1,i1,i2,i3) * phkr(1,ii) - fofr_down(2,i1,i2,i3) * phkr(2,ii)
               ur(2,ii+n1*n2*n3)= fofr_down(1,i1,i2,i3) * phkr(2,ii) + fofr_down(2,i1,i2,i3) * phkr(1,ii)
             end if
           end do
         end do
       end do
       ABI_DEALLOCATE(fofr)
       if(Dtset%nspinor==2) then
         ABI_DEALLOCATE(fofr_down)
       end if

!      === Add onsite term on the augmented FFT mesh ===
       do iatom=1,my_natom
         itypat  =local_pawfgrtab(iatom)%itypat
         lmn_size=Pawtab(itypat)%lmn_size
         ln_size =Pawtab(itypat)%basis_size ! no. of nl elements in PAW basis
         nfgd    =local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere

         ibsp=(iband-1)*dtset%nspinor
         do ispinor=1,dtset%nspinor
           ibsp=ibsp+1
           do jlmn=1,lmn_size
             jl=Pawtab(itypat)%indlmn(1,jlmn)
             jm=Pawtab(itypat)%indlmn(2,jlmn)
             cp_fact(1) = Cprj_k(iatom,ibsp)%cp(1,jlmn) *sqrt(ucvol) ! Magic factor
             cp_fact(2) = Cprj_k(iatom,ibsp)%cp(2,jlmn) *sqrt(ucvol)

             do ifgd=1,nfgd ! loop over fine grid points in current PAW sphere.
               ifftsph = local_pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid
               phj  = Paw_onsite(iatom)% phi(ifgd,jlmn)
               tphj = Paw_onsite(iatom)%tphi(ifgd,jlmn)
!              old code
!              re_p = cp_fact(1)
!              im_p = cp_fact(2)
!              apply the phase
               re_p = cp_fact(1) * phk_atm(1,ifgd,iatom) - cp_fact(2) * phk_atm(2,ifgd,iatom)
               im_p = cp_fact(1) * phk_atm(2,ifgd,iatom) + cp_fact(2) * phk_atm(1,ifgd,iatom)

               ur(1,ifftsph+(ispinor-1)*nfftot) = ur(1,ifftsph+(ispinor-1)*nfftot) + re_p * (phj-tphj)
               ur(2,ifftsph+(ispinor-1)*nfftot) = ur(2,ifftsph+(ispinor-1)*nfftot) + im_p * (phj-tphj)
               ur_ae_onsite(1,ifftsph) = ur_ae_onsite(1,ifftsph) + re_p * phj
               ur_ae_onsite(2,ifftsph) = ur_ae_onsite(2,ifftsph) + im_p * phj
               ur_ps_onsite(1,ifftsph) = ur_ps_onsite(1,ifftsph) + re_p * tphj
               ur_ps_onsite(2,ifftsph) = ur_ps_onsite(2,ifftsph) + im_p * tphj
               ur_mask(ifftsph) = one
             end do

           end do !jlmn
         end do !ispinor
       end do !iatom

       if (paral_atom) then
         ABI_ALLOCATE(buf_tmp,(2,n1*n2*n3,3))
         buf_tmp(:,:,1)=ur(:,:);buf_tmp(:,:,2)=ur_ae_onsite(:,:);buf_tmp(:,:,3)=ur_ps_onsite(:,:)
         call xmpi_sum(buf_tmp,my_comm_atom,ierr)
         ur(:,:)=buf_tmp(:,:,1);ur_ae_onsite(:,:)=buf_tmp(:,:,2);ur_ps_onsite(:,:)=buf_tmp(:,:,3)
         ABI_DEALLOCATE(buf_tmp)
       end if

!      * Remove the phase e^{ikr}, we store u(r).
#if 1
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             ii = i1 + n1*(i2-1)+ n1*n2*(i3-1)
             reur=ur(1,ii)
             imur=ur(2,ii)
             ur(1,ii)=  reur * phkr(1,ii) + imur * phkr(2,ii)
             ur(2,ii)= -reur * phkr(2,ii) + imur * phkr(1,ii)
             if(Dtset%nspinor==2) then
               reur=ur(1,ii+nfftot)    ! Important!
               imur=ur(2,ii+nfftot)
               ur(1,ii+nfftot)=  reur * phkr(1,ii) + imur * phkr(2,ii)
               ur(2,ii+nfftot)= -reur * phkr(2,ii) + imur * phkr(1,ii)
             end if
             reur=ur_ae_onsite(1,ii)
             imur=ur_ae_onsite(2,ii)
             ur_ae_onsite(1,ii)=  reur * phkr(1,ii) + imur * phkr(2,ii)
             ur_ae_onsite(2,ii)= -reur * phkr(2,ii) + imur * phkr(1,ii)
             reur=ur_ps_onsite(1,ii)
             imur=ur_ps_onsite(2,ii)
             ur_ps_onsite(1,ii)=  reur * phkr(1,ii) + imur * phkr(2,ii)
             ur_ps_onsite(2,ii)= -reur * phkr(2,ii) + imur * phkr(1,ii)
           end do
         end do
       end do
#endif

       norm=zero
       do ii=1,npw_k*Dtset%nspinor
         norm=norm+fofgin(1,ii)**2+fofgin(2,ii)**2
       end do
       write(std_out,'(a,2i5,f22.16)',advance='no') 'ikpt,iband, norm (G,PSWF)=',ikpt,iband,norm
       norm=zero
       do ifft=1,nfftot*Dtset%nspinor
         norm = norm + ur(1,ifft)**2+ur(2,ifft)**2
       end do
       norm=norm/nfftot
       norm_rerr = MAX((ABS(norm-one))*100,norm_rerr)
       write(std_out,*)"norm (R,AEWF)= ",norm
       call flush_unit(std_out)

!      MS: Various testing and debugging options
       if (.TRUE..and.nprocs==1) then

!        Dump results to .xsf files if running in serial
         tmp_unt=get_unit()

         if (present(pseudo_norms)) then
!          Check the supposedly zero overlap |\tilde{Psi_n}-\tilde{Psi_n^1}|^2
           ABI_ALLOCATE(dummy_1d,(n1*n2*n3))
           dummy_1d = zero
           norm = zero
           do ifft = 1, nfftot
             dummy_1d(ifft) = ((ur_pw(1,ifft)-ur_ps_onsite(1,ifft))**2 &
&             +  (ur_pw(2,ifft)-ur_ps_onsite(2,ifft))**2) * ur_mask(ifft)
             norm = norm + dummy_1d(ifft)
           end do
           norm = norm / nfftot
           pseudo_norms(isppol,ikpt,iband) = norm
!           if (Dtset%prtvol>9) then
!             write(std_out,'(a,3(a,I0),a,F14.9,a)') ch10,' State sp',isppol,' kpt',ikpt,' bd',iband,&
!&             ' |\tilde{Psi_n}-\tilde{Psi_n^1}|^2 norm:',norm,ch10
!             write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_PSImTPSI2','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!             open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!             call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&             Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!             close(tmp_unt)
!           end if
           ABI_DEALLOCATE(dummy_1d)
         end if


!         ABI_ALLOCATE(dummy_1d,(n1*n2*n3))
!         dummy_1d=zero
!
!         if (Dtset%prtvol>9) then
!!          Onsite AE part
!           dummy_1d = ur_ae_onsite(1,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Re_AE_onsite','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!           dummy_1d = ur_ae_onsite(2,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Im_AE_onsite','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!!          Calculate norm
!           ur_ae_onsite(1,:) = ur_ae_onsite(1,:)**2 + ur_ae_onsite(2,:)**2
!           dummy_1d = ur_ae_onsite(1,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Abs_AE_onsite','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!
!
!!          Onsite PS part-
!           dummy_1d = ur_ps_onsite(1,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Re_PS_onsite','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!           dummy_1d = ur_ps_onsite(2,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Im_PS_onsite','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!!          Calculate norm
!           ur_ps_onsite(1,:) = ur_ps_onsite(1,:)**2 + ur_ps_onsite(2,:)**2
!           dummy_1d = ur_ps_onsite(1,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Abs_PS_onsite','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!
!
!!          PW part
!           dummy_1d = ur_pw(1,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Re_PW','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!           dummy_1d = ur_pw(2,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Im_PW','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!!          Calculate norm
!           ur_pw(1,:) = ur_pw(1,:)**2 + ur_pw(2,:)**2
!           dummy_1d = ur_pw(1,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Abs_PW','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!
!!          Use ur_pw as dummy array
!           ur_pw = ur
!
!!          Full all-electron
!           dummy_1d = ur_pw(1,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Re_AE','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!           dummy_1d = ur_pw(2,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Im_AE','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!!          Calculate norm
!           ur_pw(1,:) = ur_pw(1,:)**2 + ur_pw(2,:)**2
!           dummy_1d = ur_pw(1,1:nfftot)
!           write(xsf_fname,'(2a,3(a,I0),a)') trim(my_basename),'_Abs_AE','_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!           close(tmp_unt)
!
!           if (iband==1.AND.ikpt==1) then ! This is printed only once
!!            masking function - 1 in spheres, zero outside
!             write(xsf_fname,'(a)') 'PAW_mask.xsf'
!             open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!             call printxsf(n1,n2,n3,ur_mask,Hdr%rprimd,(/zero,zero,zero/),&
!&             Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!             close(tmp_unt)
!           end if
!
!         end if ! prtvol>9
!
!         ABI_DEALLOCATE(dummy_1d)

       else
         write(msg,'(5a)')&
&         " The option to print PAW all-electron wavefunctions is on, but execution ",ch10,&
&         " is in parallel on two or more processors. XcrysDen files with individual con-",ch10,&
&         " tributions will not be written. In order to enable this you must run in serial."
!        MG
         MSG_WARNING(msg)
       end if ! Check if serial run

#if defined HAVE_TRIO_ETSF_IO
       MainFolder%wfs_rsp__spin_access   =  isppol !this is wrong if para!
       if (nprocs>1) MainFolder%wfs_rsp__spin_access = 1
       MainFolder%wfs_rsp__kpoint_access = my_kstable(ikpt,isppol) !ikpt
       MainFolder%wfs_rsp__state_access  = iband
!      main_folder%wfs_coeff__number_of_coefficients = npw * dtset%nspinor

!      We use the group level write routine.
       MainFolder%real_space_wavefunctions%data2D  => ur
       call etsf_io_main_put(ncids(1),MainFolder,lstat,Error_data=Error_data)
       ETSF_CHECK_ERROR(lstat,Error_data)

       MainFolder%real_space_wavefunctions%data2D  => ur_pw
       call etsf_io_main_put(ncids(2),MainFolder,lstat,Error_data=Error_data)
       ETSF_CHECK_ERROR(lstat,Error_data)

       MainFolder%real_space_wavefunctions%data2D  => ur_ae_onsite
       call etsf_io_main_put(ncids(3),MainFolder,lstat,Error_data=Error_data)
       ETSF_CHECK_ERROR(lstat,Error_data)

       MainFolder%real_space_wavefunctions%data2D  => ur_ps_onsite
       call etsf_io_main_put(ncids(4),MainFolder,lstat,Error_data=Error_data)
       ETSF_CHECK_ERROR(lstat,Error_data)


#endif

       ABI_DEALLOCATE(ur)
       ABI_DEALLOCATE(ur_ae_onsite)
       ABI_DEALLOCATE(ur_ps_onsite)
       ABI_DEALLOCATE(ur_pw)
       ABI_DEALLOCATE(ur_mask)

       ABI_DEALLOCATE(fofgin)
       ABI_DEALLOCATE(fofgout)
       ABI_DEALLOCATE(denpot)
     end do !nband_k

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
       ibg=ibg+dtset%nspinor*nband_k
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(phkr)

     call pawcprj_destroy(Cprj_k)
     ABI_DATATYPE_DEALLOCATE(Cprj_k)

   end do !ikpt
 end do !nsppol

!* Free augmentation waves.
 call paw_pwaves_lmn_free(Paw_onsite)
 ABI_DATATYPE_DEALLOCATE(Paw_onsite)

 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(l_size_atm)

 ABI_DEALLOCATE(my_kpoints)
 ABI_DEALLOCATE(my_spins)
 ABI_DEALLOCATE(my_kstable)

!* Maximum relative error over CPUs.
 call xmpi_max(norm_rerr,max_rerr,spaceComm,ierr)
 write(std_out,*)"max_rerr=",max_rerr
 if (max_rerr>ten) then
   write(msg,'(7a)')&
&   " Inaccuracy on the normalization of the wave funtions exceeds 10%. ",ch10,&
&   " Likely due to the use of a too coarse FFT mesh or unconverged wavefunctions. ",ch10,&
&   " Numerical values inside the augmentation regions might be inaccurate. ",ch10,&
&   " Action: increase pawecutdg in your input file. "
   MSG_COMMENT(msg)
 end if

#if defined HAVE_TRIO_ETSF_IO
!=== Merge partial files ===
 if (nprocs>1) then
   call xmpi_barrier(spaceComm)
   if (rank==0) then
     do ii = 1, 4
       write(msg,'(2a)')'Master node is merging NETCDF partial files into: ',TRIM(out_files(ii))
       call wrtout(std_out, msg,'COLL')
       call etsf_io_file_merge(out_files(ii),merge_files(ii,:),lstat,Error_data)
       ETSF_CHECK_ERROR(lstat,Error_data)
     end do
   end if
   call xmpi_barrier(spaceComm)
 end if

 if (allocated(merge_files))  then
   ABI_DEALLOCATE(merge_files)
 end if
#endif

 ABI_DEALLOCATE(r0shift)
 ABI_DEALLOCATE(phk_atm)
 ABI_DEALLOCATE(xcart)
 call pawfgrtab_destroy(local_pawfgrtab)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawmkaewf
!!***
