!{\src2tex{textfont=tt}}
!!****f* ABINIT/build_spectra
!! NAME
!!  build_spectra
!!
!! FUNCTION
!!  Driver routine for the computation of optical spectra.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  usepaw=1 for PAW calculations, 0 otherwise.
!!  drude_plsmf=Drude plasma frequency.
!!  Bsp<excparam>=Data type gathering the paramenters used for the Bethe-Salpeter calculation.
!!    inclvkb=If different from 0, [Vnl,r] is included in the calculation of the matrix elements of the velocity operator.
!!  BS_files<excfiles>=filenames used in the Bethe-Salpeter part.
!!  Kmesh<kmesh_t>=the k-point sampling for the wave functions.
!!  Cryst<crystal_t>=Structure defining the crystalline structure.
!!  KS_BSt=The KS energies.
!!  QP_BSt=The QP energies.
!!  Psps <pseudopotential_type>=variables related to pseudopotentials.
!!  Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data
!!  Hur(Cryst%natom*usepaw)<HUr_commutator>=Only for PAW and LDA+U, quantities used to evaluate the commutator [H_u,r].
!!  Wfd<wfd_t>=Handler for the wavefunctions.
!!    nsppol=Number of independent spin polarizations.
!!    nspinor=Number of spinorial components.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  No output. The routine calls specialized routines where the computation and the output of the spectra is done.
!!
!! PARENTS
!!      m_exc_diago
!!
!! CHILDREN
!!      c_f_pointer,etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine build_spectra(BSp,BS_files,Cryst,Kmesh,KS_BSt,QP_BSt,Psps,Pawtab,Wfd,Hur,drude_plsmf,comm)

 use defs_basis 
 use m_bs_defs
 use defs_datatypes   
 use defs_abitypes   
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_ncfile
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif

 use m_crystal,         only : crystal_t 
 use m_crystal_io,      only : crystal_ncwrite
 use m_bz_mesh,         only : kmesh_t
 use m_ebands,          only : ebands_ncwrite
 use m_commutator_vkbr, only : kb_potential
 use m_pawtab,          only : pawtab_type
 use m_paw_commutator,  only : HUr_commutator
 use m_wfs,             only : wfd_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'build_spectra'
 use interfaces_14_hidewrite
 use interfaces_69_wfdesc
 use interfaces_71_bse, except_this_one => build_spectra
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 real(dp),intent(in) :: drude_plsmf
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(pseudopotential_type),intent(in) :: Psps
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(ebands_t),intent(in) :: KS_BSt,QP_BSt
 type(wfd_t),intent(inout) :: Wfd
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(HUr_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)
 
!Local variables ------------------------------
!scalars
 integer :: my_rank,master,iq,io,nsppol,lomo_min,max_band
 real(dp) :: omegaev
 complex(dpc) :: ks_avg,gw_avg,exc_avg
 !character(len=500) :: msg
 type(ncfile_t) :: ncf
!arrays
 real(dp),allocatable :: dos_exc(:),dos_gw(:),dos_ks(:)
 complex(dpc),allocatable :: eps_rpanlf(:,:),eps_gwnlf(:,:)
 complex(dpc),allocatable :: eps_exc(:,:),opt_cvk(:,:,:,:,:)
      
!************************************************************************

 my_rank = Wfd%my_rank
 master  = Wfd%master
 nsppol  = Wfd%nsppol
 !
 ! =====================================================
 ! === Calculate fcv(k)=<c k s|e^{-iqr}|v k s> in BZ ===
 ! =====================================================
 lomo_min=Bsp%lomo_min; max_band=Bsp%nbnds
 ABI_MALLOC(opt_cvk,(lomo_min:max_band,lomo_min:max_band,BSp%nkbz,nsppol,BSp%nq))

 do iq=1,BSp%nq 
   call calc_optical_mels(Wfd,Kmesh,KS_BSt,Cryst,Psps,Pawtab,Hur,BSp%inclvkb,Bsp%lomo_spin,lomo_min,max_band,&
&                         BSp%nkbz,BSp%q(:,iq),opt_cvk(:,:,:,:,iq))
 end do
 !
 ! ============================
 ! ==== Make EPS EXCITONIC ====
 ! ============================
 if (my_rank==master) then ! Only master works.

   ABI_MALLOC(eps_exc,(BSp%nomega,BSp%nq))
   ABI_MALLOC(dos_exc,(BSp%nomega))

   if (BSp%use_coupling==0) then
     call exc_eps_resonant(BSp,BS_files,lomo_min,max_band,BSp%nkbz,nsppol,opt_cvk,Cryst%ucvol,BSp%nomega,BSp%omega,eps_exc,dos_exc)
   else
     call exc_eps_coupling(Bsp,BS_files,lomo_min,max_band,BSp%nkbz,nsppol,opt_cvk,Cryst%ucvol,BSp%nomega,BSp%omega,eps_exc,dos_exc)
   end if
   !
   ! =======================================================
   ! === Make EPS RPA and GW without local-field effects ===
   ! =======================================================
   call wrtout(std_out," Calculating RPA NLF and QP NLF epsilon","COLL")

   ABI_MALLOC(eps_rpanlf,(BSp%nomega,BSp%nq))
   ABI_MALLOC(dos_ks,(BSp%nomega))

   call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,Bsp%lomo_min,BSp%homo_spin,Kmesh,KS_BSt,BSp%nq,nsppol,opt_cvk,&
&    Cryst%ucvol,BSp%broad,BSp%nomega,BSp%omega,eps_rpanlf,dos_ks)

   ABI_MALLOC(eps_gwnlf ,(BSp%nomega,BSp%nq))
   ABI_MALLOC(dos_gw,(BSp%nomega))

   call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,Bsp%lomo_min,BSp%homo_spin,Kmesh,QP_BSt,BSp%nq,nsppol,opt_cvk,&
&    Cryst%ucvol,Bsp%broad,BSp%nomega,BSp%omega,eps_gwnlf,dos_gw)
   !
   ! =========================
   ! === Write out Epsilon ===
   ! =========================
   !this is just for the automatic tests, It will be removed when fldiff 
   !will be able to compare two optical spectral
   write(ab_out,*)" "
   write(ab_out,*)"Macroscopic dielectric function:"
   write(ab_out,*)"omega [eV] <KS_RPA_nlf>  <GW_RPA_nlf>  <BSE> "
   do io=1,MIN(10,BSp%nomega)
     omegaev = REAL(BSp%omega(io))*Ha_eV
     ks_avg  = SUM( eps_rpanlf(io,:)) / Bsp%nq
     gw_avg  = SUM( eps_gwnlf (io,:)) / Bsp%nq
     exc_avg = SUM( eps_exc   (io,:)) / Bsp%nq
     write(ab_out,'(7f9.4)')omegaev,ks_avg,gw_avg,exc_avg
   end do
   write(ab_out,*)" "
   !
   ! Master node writes final results on file.
   call exc_write_data(BSp,BS_files,"RPA_NLF_MDF",eps_rpanlf,dos=dos_ks)

   call exc_write_data(BSp,BS_files,"GW_NLF_MDF",eps_gwnlf,dos=dos_gw)

   call exc_write_data(BSp,BS_files,"EXC_MDF",eps_exc,dos=dos_exc)
 
   call wrtout(std_out," Checking Kramers Kronig on Excitonic Macroscopic Epsilon","COLL")
   call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_exc(:,1))

   call wrtout(std_out," Checking Kramers Kronig on RPA NLF Macroscopic Epsilon","COLL")
   call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_rpanlf(:,1))

   call wrtout(std_out," Checking Kramers Kronig on GW NLF Macroscopic Epsilon","COLL")
   call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_gwnlf(:,1))

   call wrtout(std_out," Checking f-sum rule on Excitonic Macroscopic Epsilon","COLL")

   if (BSp%exchange_term>0) then 
     MSG_COMMENT(' f-sum rule should be checked without LF')
   end if
   call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_exc(:,1)),drude_plsmf)

   call wrtout(std_out," Checking f-sum rule on RPA NLF Macroscopic Epsilon","COLL")
   call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_rpanlf(:,1)),drude_plsmf)

   call wrtout(std_out," Checking f-sum rule on GW NLF Macroscopic Epsilon","COLL")
   call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_gwnlf(:,1)),drude_plsmf)

#ifdef HAVE_TRIO_ETSF_IO
     NCF_CHECK(ncfile_create(ncf,TRIM(BS_files%out_basename)//"_MDF.nc", NF90_CLOBBER), "Creating MDF file")
     call crystal_ncwrite(Cryst,ncf%ncid)
     ! TODO: Write bands but try to avoid having to pass dtset!
     !call ebands_ncwrite(KS_BSt, dtset%nshiftk_orig,dtset%shiftk_orig,dtset%nshiftk,dtset%shiftk,dtset%ngkpt,dtset%kptrlatt,ncf%ncid)
     call mdfs_ncwrite(ncf%ncid, Bsp, eps_exc,eps_rpanlf,eps_gwnlf)
     NCF_CHECK(ncfile_close(ncf),"Closing MDF file")
#else
     ABI_UNUSED(ncf%ncid)
#endif

   ABI_FREE(eps_rpanlf)
   ABI_FREE(eps_gwnlf)
   ABI_FREE(eps_exc)
   ABI_FREE(dos_exc)
   ABI_FREE(dos_ks)
   ABI_FREE(dos_gw)
 end if ! my_rank==master

 ABI_FREE(opt_cvk)

 call xmpi_barrier(comm)

end subroutine build_spectra
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/exc_write_data
!! NAME
!!  exc_write_data
!!
!! FUNCTION
!!  This routine drives the writing of the files produced by the Bethe-Salpeter code. 
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! BSp<excparam>=Bethe-Salpeter Parameters.
!! what= "EXC_MDF"
!!       "RPA_NLF_MDF"
!!       "GW_NLF_MDF"
!! [dos(nomega)]
!!
!! OUTPUT
!!  Only writing.
!!
!! SIDE EFFECTS
!!  eps(BSp%nomega,BSp%nq) = Macroscopic dielectric function to be written.
!!
!! PARENTS
!!      exc_spectra,m_haydock
!!
!! CHILDREN
!!      c_f_pointer,etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine exc_write_data(BSp,BS_files,what,eps,dos)

 use defs_basis
 use m_bs_defs
 use m_errors
 use m_profiling_abi

 use m_io_tools,        only : open_file
 use m_fstrings,        only : toupper, strcat
 use m_numeric_tools,   only : simpson_int

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_write_data'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: what
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
!arrays
 real(dp),optional,intent(in) :: dos(BSp%nomega)
 complex(dpc),intent(in) :: eps(BSp%nomega,BSp%nq)

!Local variables ------------------------------
!scalars
 integer :: io,iq,funt 
 real(dp) :: omegaev,step
 !real(dp),parameter :: SMALL=5.0d-99
!arrays
 real(dp) :: int_dos(BSp%nomega)
 real(dp) :: tmp_eps(2,BSp%nq)
 character(len=500) :: lf_type,block_type,wgg_type,frm,str_type,msg 
 character(len=fnlen) :: fname

!************************************************************************

 fname = strcat(BS_files%out_basename,'_',toupper(what))

 if (open_file(fname,msg,newunit=funt,form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 select case (toupper(what))
 case ("EXC_MDF")
   call wrtout(ab_out," Writing EXC Macroscopic dielectric function to file: "//trim(fname),"COLL")

   write(funt,'("# Macroscopic dielectric function obtained with the BS equation.")')

   lf_type = 'WITHOUT LOCAL FIELD EFFECTS'
   if (BSp%exchange_term>0) lf_type='LOCAL FIELD EFFECTS INCLUDED'
   call bsp_calctype2str(Bsp,str_type)
   write(funt,'("# ",a,"     " ,a)') TRIM(str_type), TRIM(lf_type)

   block_type = 'RESONANT-ONLY calculation'
   if (BSp%use_coupling>0) block_type = 'RESONANT+COUPLING calculation'
   write(funt,'("# ",a)') TRIM(block_type)

   if (BSp%use_coulomb_term) then
     wgg_type = "Coulomb term constructed with full W(G1,G2)"
     if ( BSp%use_diagonal_Wgg ) wgg_type = "Coulomb term constructed with diagonal approximation W(G1,G1)"
     write(funt,'("# ",a)') TRIM(wgg_type)
   end if

   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%soenergy*Ha_eV,' [eV]'

 case ("RPA_NLF_MDF")
   call wrtout(ab_out," Writing KS-RPA macroscopic dielectric function without local fields to file: "//trim(fname),"COLL")
   write(funt,'("# RPA macroscopic dielectric function without local fields")')

 case ("GW_NLF_MDF")
   call wrtout(ab_out," Writing GW-RPA macroscopic dielectric function without local fields to file: "//trim(fname),"COLL")

   write(funt,'("# GW Macroscopic dielectric function without local field effects ")')
   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%soenergy*Ha_eV,' [eV]'

 case default
   MSG_ERROR("Unknown value for what: "//trim(what))
 end select
 !
 ! Paramaters common to the different calculations.
 if (BSp%algorithm /= BSE_ALGO_HAYDOCK) then
   write(funt,'(a,i0)')"# nstates included in the diagonalization = ",BSp%nstates    
 end if
 
 if (BSp%algorithm == BSE_ALGO_HAYDOCK) then
   write(funt,'(a,2f7.4)')'# Tolerance = ',BSp%haydock_tol
 end if

 write(funt,'(a,i0)')"# npweps  = ",BSp%npweps
 write(funt,'(a,i0)')"# npwwfn  = ",BSp%npwwfn
 write(funt,'(a,i0)')"# nbands  = ",BSp%nbnds
 write(funt,'(a,i0)')"# loband  = ",BSp%lomo_spin(1)
 if (Bsp%nsppol==2) write(funt,'(a,i0)')"# loband(spin=2) = ",BSp%lomo_spin(2)
 write(funt,'(a,i0)')"# nkibz   = ",BSp%nkibz
 write(funt,'(a,i0)')"# nkbz    = ",BSp%nkbz
 write(funt,'(a,f7.4,a)')'# Lorentzian broadening = ',BSp%broad*Ha_eV,' [eV]'
 !
 ! Write the list of q-points.
 write(funt,'(a)')"# List of q-points for the optical limit:"
 do iq=1,BSp%nq
   write(funt,'(a,3(f9.6,","),a)')'# q = ',BSp%q(:,iq),' [Reduced coords] '
 end do
 !
 ! Write spectra.
 if (.not.PRESENT(dos)) then
   write(funt,'(a)')"# omega [eV]    RE(eps(q=1)) IM(eps(q=1) RE(eps(q=2) ) ... "
   !write(frm,*)'(f7.3,',2*BSp%nq,'es12.4)'
   write(frm,*)'(f7.3,',2*BSp%nq,'(1x,f9.4))'
   do io=1,BSp%nomega
     omegaev = DBLE(BSp%omega(io))*Ha_eV
     tmp_eps(1,:) = REAL (eps(io,:))
     tmp_eps(2,:) = AIMAG(eps(io,:))
     !where (ABS(tmp_eps) < SMALL) ! this to improve the portability of the automatic tests.
     !  tmp_eps = zero
     !end where
     write(funt,frm) omegaev,(tmp_eps(:,iq), iq=1,BSp%nq)
   end do

 else 
   write(funt,'(a)')"# omega [eV]    RE(eps(q=1)) IM(eps(q=1) RE(eps(q=2) ) ... DOS   IDOS"
   step = DBLE(BSp%omega(2) - BSp%omega(1))
   if ( ABS( step - DBLE((BSp%omega(BSp%nomega) - BSp%omega(BSp%nomega-1)))) > tol6 ) then
     MSG_WARNING("Frequency mesh must be linear for using simpson_int")
   end if
   call simpson_int(Bsp%nomega,step,dos,int_dos)
   !write(frm,*)'(f7.3,',2*BSp%nq,'es12.4,2es12.4)'
   write(frm,*)'(f7.3,',2*BSp%nq,'(1x,f9.4,1x,f9.4,1x,f9.4))'
   do io=1,BSp%nomega
     omegaev = DBLE(BSp%omega(io))*Ha_eV
     tmp_eps(1,:) = REAL (eps(io,:))
     tmp_eps(2,:) = AIMAG(eps(io,:))
     !where (ABS(tmp_eps) < SMALL) ! this to improve the portability of the automatic tests.
     !  tmp_eps = zero
     !end where
     !write(funt,frm) omegaev,(eps(io,iq), iq=1,BSp%nq), dos(io), int_dos(io)
     write(funt,frm) omegaev,(tmp_eps(:,iq), iq=1,BSp%nq), dos(io), int_dos(io)
   end do
 end if

 close(funt)

end subroutine exc_write_data
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/exc_eps_rpa
!! NAME
!!  exc_eps_rpa
!!
!! FUNCTION
!!  Make epsilon within RPA and GW. 
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nkbz=Number of points in the BZ 
!! nbnds=Number of bands
!! lomo_spin(nsppol)
!! lomo_min=Lowest occupied state
!! homo=Number of occupied states.
!! homo_spin(nsppol)
!! nsppol=Number of independent spin polarizations.
!! nomega=Number of frequencies
!! omega(nomega)=Frequency mesh.
!! ucvol=Unit cell volume.
!! broad=Broadening used for the DOS.
!! opt_cvk(nbnds,nbnds,nkbz)=Matrix elements <b k|e^{-iqr}|b" k> for a given q in the full BZ.
!!
!! OUTPUT
!!  eps_rpa(nomega)=RPA spectrum without local-field effects.
!!  dos(nomega)=The DOS.
!!
!! PARENTS
!!      exc_spectra,m_haydock
!!
!! CHILDREN
!!      c_f_pointer,etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine exc_eps_rpa(nbnds,lomo_spin,lomo_min,homo_spin,Kmesh,Bst,nq,nsppol,opt_cvk,ucvol,broad,nomega,omega,eps_rpa,dos)

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_profiling_abi

 use m_bz_mesh,         only : kmesh_t
 use m_special_funcs,   only : gaussian

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_eps_rpa'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,lomo_min,nsppol,nomega,nq
 real(dp),intent(in) :: ucvol,broad
 type(kmesh_t),intent(in) :: Kmesh
 type(ebands_t),intent(in) :: BSt
!arrays
 integer,intent(in) :: lomo_spin(nsppol),homo_spin(nsppol)
 real(dp),intent(out) :: dos(nomega)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(in) :: opt_cvk(lomo_min:nbnds,lomo_min:nbnds,Kmesh%nbz,nsppol,nq)
 complex(dpc),intent(out) :: eps_rpa(nomega,nq)

!Local variables ------------------------------
!scalars
 integer :: iw,ib_v,ib_c,ik_bz,ik_ibz,spin,iq
 real(dp) :: fact,arg,ediff
 complex(dpc) :: ctemp

!************************************************************************

 ! TODO: four_pi comes from the bare Coulomb term hence the 
 ! present implementation is not compatible with the cutoff technique.
 fact=four_pi/(ucvol*Kmesh%nbz)
 if (nsppol==1) fact=two*fact ! two accounts for the occupation factors.

 eps_rpa=czero; dos=zero

 !write(std_out,*)nsppol,Kmesh%nbz,lomo_min,homo,nbnds
 !
 ! Sum over all QP transitions.
 do spin=1,nsppol
   do ik_bz=1,Kmesh%nbz
     ik_ibz = Kmesh%tab(ik_bz)
     do ib_v=lomo_spin(spin),homo_spin(spin)
       do ib_c=homo_spin(spin)+1,nbnds
         !
         ! TODO here energies are always assumed to be real.
         ediff = BSt%eig(ib_c,ik_ibz,spin) - BSt%eig(ib_v,ik_ibz,spin)
         !
         do iq=1,nq
           ctemp = opt_cvk(ib_c,ib_v,ik_bz,spin,iq)
           do iw=1,nomega
             eps_rpa(iw,iq) = eps_rpa(iw,iq)  + ctemp * CONJG(ctemp) * (one/(ediff-omega(iw)) + one/(ediff+omega(iw)))
           end do
         end do
         !
         ! The JDOS at q=0
         !if (ediff*Ha_eV < 0.3) then
         !  write(std_out,*)"Small transition ",ik_ibz,ib_v,ib_c
         !end if

         do iw=1,nomega
           arg = DBLE(omega(iw)) - ediff
           dos(iw) = dos(iw) + gaussian(arg,broad)        
         end do
         !
       end do !ib_c
     end do !ib_v
   end do !ik_bz
 end do !spin

 dos = dos/Kmesh%nbz
 eps_rpa = cone + fact*eps_rpa

end subroutine exc_eps_rpa
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/exc_eps_resonant
!! NAME
!!  exc_eps_resonant
!!
!! FUNCTION
!!  This routine calculates the macroscopic dielectric function with excitonic effects.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Bsp
!! lomo_min,max_band
!! nkbz=Number of points in the BZ 
!! nsppol=Number of independent polarizations.
!! nomega=Number of frequencies
!! omega(nomega)=frequency mesh (complex shift is already included)
!! ucvol=Volume of the unit cell.
!! opt_cvk(lomo_min:max_band,mib:max_band,nkbz,nsppol,Bsp%nq)=Matrix elements <b k|e^{-iqr}|b" k> for a given q in the full BZ.
!!
!! OUTPUT
!!  eps_exc(nomega,Bsp%nq)=Macroscopic dielectric function with excitonic effects.
!!  dos_exc(nomega)=The DOS of the excitonic Hamiltonian
!!
!! PARENTS
!!      exc_spectra
!!
!! CHILDREN
!!      c_f_pointer,etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE


subroutine exc_eps_resonant(Bsp,BS_files,lomo_min,max_band,nkbz,nsppol,opt_cvk,ucvol,nomega,omega,eps_exc,dos_exc)

 use defs_basis
 use m_bs_defs
 use m_errors
 use m_profiling_abi

 use defs_abitypes,    only : hdr_type
 use m_io_tools,       only : open_file
 use m_fstrings,       only : strcat
 use m_special_funcs,  only : gaussian
 use m_header,         only : hdr_free
 use m_bse_io,         only : exc_amplitude

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_eps_resonant'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lomo_min,max_band,nkbz,nomega,nsppol
 real(dp),intent(in) :: ucvol
 type(excfiles),intent(in) :: BS_files
 type(excparam),intent(in) :: BSp
!arrays
 real(dp),intent(out) :: dos_exc(nomega)
 complex(dpc),intent(in) :: opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,nsppol,BSp%nq),omega(nomega)
 complex(dpc),intent(out) :: eps_exc(nomega,BSp%nq)

!Local variables ------------------------------
!scalars
 integer :: ll,it,iw,ib_v,ib_c,ik_bz,neig_read,eig_unt,exc_size,iq !fform,
 integer :: spin,spad,hsize_read,nstates,ost_unt
 real(dp) :: fact,arg
 character(len=500) :: msg,frm
 character(len=fnlen) :: filbseig,ost_fname
 !type(Hdr_type) :: tmp_Hdr
!arrays
 real(dp),allocatable :: exc_ene(:),ostrength(:,:)
 complex(dpc) :: ctemp(BSp%nq)
 complex(dpc),allocatable :: exc_ene_cplx(:),exc_state(:)

!************************************************************************

 call wrtout(std_out," Calculating excitonic epsilon with antiresonant","COLL")

 if (nsppol==2) then
   MSG_WARNING("nsppol==2 still under development")
 end if

 exc_size = SUM(BSp%nreh) 
 nstates  = BSp%nstates

 if (ANY(Bsp%nreh/=Bsp%nreh(1))) then
   write(msg,'(a,2(i0,1x))')"BSE does not support different number of transitions for the two spin channels. nreh: ",Bsp%nreh
   MSG_WARNING(msg)
 end if
 !
 ! TODO: 
 ! four_pi comes from the bare Coulomb term hence the 
 ! present implementation is not compatible with the cutoff technique.
 fact=four_pi/(ucvol*nkbz); if (nsppol==1) fact=two*fact ! two to account for the occupation numbers.

 if (BS_files%in_eig /= BSE_NOFILE) then
   filbseig = BS_files%in_eig
 else 
   filbseig = BS_files%out_eig
 end if

 call wrtout(std_out," Reading excitonic eigenstates from file: "//TRIM(filbseig),"COLL")
 if (open_file(filbseig,msg,newunit=eig_unt,form="unformatted",status="old",action="read") /= 0) then
   MSG_ERROR(msg)
 end if

 !$call hdr_io_int(fform,tmp_Hdr,1,eig_unt)
 !$call hdr_free(tmp_Hdr)

 read(eig_unt)hsize_read,neig_read

 if (hsize_read /= exc_size) then
   write(msg,'(2(a,i0))')" Wrong size of the Hamiltonian: read: ",hsize_read," expected= ",exc_size
   MSG_ERROR(msg)
 end if

 if (neig_read /= nstates) then
   write(msg,'(2(a,i0))')" Wrong number of eigenstates: read: ",neig_read," expected= ",nstates
   MSG_ERROR(msg)
 end if
 !
 ! Read eigenvalues, ignore possibly small imaginary part.
 ABI_MALLOC(exc_ene_cplx,(neig_read))
 read(eig_unt) exc_ene_cplx 

 ABI_MALLOC(exc_ene,(neig_read))
 exc_ene = exc_ene_cplx
 ABI_FREE(exc_ene_cplx)
 !
 ! Calculate oscillator strength.
 ABI_MALLOC(exc_state,(exc_size))
 ABI_MALLOC(ostrength,(neig_read,BSp%nq))

 do ll=1,neig_read ! Loop over excitonic eigenstates reported on file.
   read(eig_unt) exc_state(:)

   ctemp(:) = czero 
   do spin=1,nsppol
     spad=(spin-1)*BSp%nreh(1) ! Loop over spin channels.
     do it=1,BSp%nreh(spin)    ! Loop over resonant transition t = (k,v,c,s)
       ik_bz = Bsp%Trans(it,spin)%k
       ib_v  = Bsp%Trans(it,spin)%v
       ib_c  = Bsp%Trans(it,spin)%c
       do iq=1,BSp%nq
         ctemp(iq) = ctemp(iq) + CONJG(opt_cvk(ib_c,ib_v,ik_bz,spin,iq)) * exc_state(it+spad)
       end do
     end do ! it
   end do
   ostrength(ll,:) = ctemp(:)*CONJG(ctemp(:))
 end do ! ll

 close(eig_unt)
 ABI_FREE(exc_state)

 eps_exc = one
 do ll=1,neig_read ! Sum over all excitonic eigenstates read from file.
   do iq=1,BSp%nq
      do iw=1,nomega
        eps_exc(iw,iq) = eps_exc(iw,iq) +  &
&         fact * ostrength(ll,iq) * (one/(exc_ene(ll) - omega(iw)) - one/(-exc_ene(ll) - omega(iw)))
      end do
   end do !ll
 end do !iw
 !
 ! The excitonic DOS.
 dos_exc=zero
 do ll=1,neig_read ! Sum over the calculate excitonic eigenstates.
   do iw=1,nomega
     arg = ( DBLE(omega(iw)) - exc_ene(ll))
     dos_exc(iw) = dos_exc(iw) + gaussian(arg,Bsp%broad)
   end do
 end do
 !
 ! Write the oscillator strengths of file.
 ost_fname = strcat(BS_files%out_basename,"_EXC_OST")

 if (open_file(ost_fname,msg,newunit=ost_unt,form="formatted",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(ost_unt,'("# Oscillator strengths of the excitonic states for the different q-polarizations.")')
 !
 ! Write the list of q-points.
 write(ost_unt,*)"# List of q-points for the optical limit"
 do iq=1,BSp%nq
   write(ost_unt,'(a,3(f9.6,","),a)')'# q = ',BSp%q(:,iq),' [Reduced coords] '
 end do

 write(ost_unt,*)"# E_lambda [eV]     ostrength(q=1) ostrength(q=2) .... "
 write(frm,*)'(f8.4,',BSp%nq,'es12.4)'
 do ll=1,neig_read 
   write(ost_unt,frm)exc_ene(ll)*Ha_eV,(ostrength(ll,iq), iq=1,BSp%nq)
 end do

 close(ost_unt)

 ABI_FREE(ostrength)
 ABI_FREE(exc_ene)

 !call exc_amplitude(Bsp,filbseig,1,(/(ll,ll=1,10)/),"TEST_AMPLITUDE")
 !call exc_amplitude(Bsp,filbseig,1,(/30/),"TEST_AMPLITUDE")

end subroutine exc_eps_resonant
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/exc_eps_coupling
!! NAME
!!  exc_eps_coupling
!!
!! FUNCTION
!!  Make epsilon EXCITONIC with full COUPLING.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Bsp
!! nkbz=Number of points in the BZ 
!! lomo_min,max_band
!! nomega=Number of frequencies
!! omega(nomega)=frequency mesh.
!! nsppol=Number of independent spin polarizations.
!! ucvol=Unit cell volume.
!! BS_files<excfiles>File names used in the Bethe-Salpeter code.
!! opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,nsppol)=Matrix elements <b k|e^{-iqr}|b" k> for a given q in the full BZ.
!!
!! OUTPUT
!!  eps_exc(nomega)=Macroscopic dielectric function with excitonic effects calculated including the COUPLING.
!!  dos_exc(nomega)=The DOS of the excitonic Hamiltonian
!!
!! PARENTS
!!      exc_spectra
!!
!! CHILDREN
!!      c_f_pointer,etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE


subroutine exc_eps_coupling(Bsp,BS_files,lomo_min,max_band,nkbz,nsppol,opt_cvk,ucvol,nomega,omega,eps_exc,dos_exc)

 use defs_basis
 use m_bs_defs
 use m_errors
 use m_profiling_abi

 use m_io_tools,       only : open_file
 use m_fstrings,       only : strcat
 use m_special_funcs,  only : gaussian
 use m_blas,           only : xdotu
 use defs_abitypes,    only : hdr_type
 use m_header,         only : hdr_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_eps_coupling'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lomo_min,max_band,nkbz,nomega,nsppol
 real(dp),intent(in) :: ucvol
 type(excfiles),intent(in) :: BS_files
 type(excparam),intent(in) :: BSp
!arrays
 real(dp),intent(out) :: dos_exc(nomega)
 complex(dpc),intent(in) :: opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,nsppol,BSp%nq),omega(nomega)
 complex(dpc),intent(out) :: eps_exc(nomega,BSp%nq)

!Local variables ------------------------------
!scalars
 integer :: mi,it,ii,ib_v,ib_c,ik_bz,exc_size_read,nstates_read,eig_unt !,fform
 integer :: exc_size,iq,spin,tr_idx,tar_idx,nstates,iw,ll
 real(dp) :: fact,arg
 complex(dpc) :: eps,fam,famp
 character(len=500) :: msg
 character(len=fnlen) :: filbseig
 !type(Hdr_type) :: tmp_Hdr
!arrays
 complex(dpc),allocatable :: Ami(:),exc_ene(:),Sm1mi(:)
 complex(dpc),allocatable :: msfap(:,:),fa(:,:),fap(:,:)

!************************************************************************

 call wrtout(std_out," Calculating absorption strength with full coupling","COLL")

 if (nsppol==2) then 
   MSG_WARNING("nsppol==2 is still under development")
 end if

 ! Rank of the entire excitonic Hamiltonian including the coupling block.
 exc_size = 2*SUM(BSp%nreh); if (nsppol==2) exc_size = 2*(SUM(BSp%nreh) + BSp%nreh(2)) 
 nstates  = BSp%nstates

 ! TODO: four_pi comes from the bare Coulomb term hence the 
 ! present implementation is not compatible with the cutoff technique.
 ! factor two is due to the occupation factors.
 fact=four_pi/(ucvol*nkbz); if (nsppol==1) fact=two*fact

 if (BS_files%in_eig /= BSE_NOFILE) then
   filbseig = BS_files%in_eig
 else 
   filbseig = BS_files%out_eig
 end if

 call wrtout(std_out," Reading excitonic eigenstates from file: "//trim(filbseig),"COLL")
 if (open_file(filbseig,msg,newunit=eig_unt,form="unformatted", status="old", action="read") /= 0) then
   MSG_ERROR(msg)
 end if
 !$call hdr_io_int(fform,tmp_Hdr,1,eig_unt)
 !$call hdr_free(tmp_Hdr)

 read(eig_unt) exc_size_read, nstates_read
 ABI_CHECK(exc_size_read==exc_size,"wrong file")
 ABI_CHECK(nstates_read==nstates,"Partial diago not supported yet")
 !
 ! Read eigenvalues
 ABI_MALLOC(exc_ene,(nstates))
 read(eig_unt) exc_ene(:)

 ABI_MALLOC(Ami,(exc_size))
 ABI_MALLOC(fa,(nstates,BSp%nq))
 ABI_MALLOC(fap,(nstates,BSp%nq))
 ABI_CHECK_ALLOC(" out-of-memory Ami")

 do mi=1,nstates ! Loop on excitonic eigenvalues mi
   read(eig_unt) Ami(:)
   !
   do iq=1,BSp%nq
     !
     fam  = czero
     famp = czero
     do spin=1,nsppol
       do it=1,BSp%nreh(spin) ! Loop over transition t = (k,v,c)
         ik_bz = Bsp%Trans(it,spin)%k
         ib_v  = Bsp%Trans(it,spin)%v
         ib_c  = Bsp%Trans(it,spin)%c
         tr_idx  = it + (spin-1)*Bsp%nreh(1)
         if (nsppol==1) then
           tar_idx = it + Bsp%nreh(1)
         else
           if (spin==1) tar_idx = it + SUM(Bsp%nreh)
           if (spin==2) tar_idx = it + 2*Bsp%nreh(1)+Bsp%nreh(2)
         end if

         fam = fam + CONJG(opt_cvk(ib_c,ib_v,ik_bz,spin,iq)) * Ami(tr_idx) &
&                  + CONJG(opt_cvk(ib_v,ib_c,ik_bz,spin,iq)) * Ami(tar_idx)

         famp = famp - opt_cvk(ib_c,ib_v,ik_bz,spin,iq) * CONJG(Ami(tr_idx)) &
&                    + opt_cvk(ib_v,ib_c,ik_bz,spin,iq) * CONJG(Ami(tar_idx))
       end do
     end do
     ! Save results.
     fa (mi,iq) = fam
     fap(mi,iq) = famp
   end do
 end do ! mi

 ABI_FREE(Ami)

 !$call hdr_io_int(fform,tmp_Hdr,1,eig_unt)
 !$call hdr_free(tmp_Hdr)

 ! Read O{-1} and sum over the eigenstates.
 ABI_MALLOC(Sm1mi,(nstates))
 ABI_MALLOC(msfap,(nstates,BSp%nq))
 ABI_CHECK_ALLOC("out-of-memory Sm1mi")

 do mi=1,nstates
   read(eig_unt) Sm1mi
   Sm1mi = DCONJG(Sm1mi) ! This gives the row since O^{-1} is Hermitian.
   do iq=1,BSp%nq
     msfap(mi,iq) = xdotu(exc_size,Sm1mi,1,fap(:,iq),1)
   end do
 end do

 ABI_FREE(Sm1mi)

 close(eig_unt)
 !
 ! === Calculate excitonic epsilon with coupling ===
 do iq=1,BSp%nq
   !
   do ii=1,nomega
     eps = czero
     do mi=1,nstates ! sum over all exciton eigenstates
       eps = eps - fa(mi,iq) * msfap(mi,iq) / (exc_ene(mi) - omega(ii))
     end do
     eps_exc(ii,iq) = one + fact * eps
   end do
   !
 end do

 ABI_FREE(fa)
 ABI_FREE(msfap)
 ABI_FREE(fap)
 !
 ! The excitonic DOS.
 dos_exc=zero
 do ll=1,nstates ! Sum over the calculate excitonic eigenstates.
   do iw=1,nomega
     arg = DBLE(omega(iw) - exc_ene(ll))
     dos_exc(iw) = dos_exc(iw) + gaussian(arg,Bsp%broad)
   end do
 end do

 ABI_FREE(exc_ene)

end subroutine exc_eps_coupling
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/exc_write_tensor
!! NAME
!!  exc_write_tensor
!!
!! FUNCTION
!!  This routine drives the writing of complex dielectric tensor
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (Y. Gillet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! BSp<excparam>=Bethe-Salpeter Parameters.
!! what= "EXC_TSR_CART" or "EXC_TSR_RED"
!!       "RPA_NLF_TSR_CART" or "RPA_NLF_TSR_RED"
!!       "GW_NLF_TSR_CART" or "GW_NLF_TSR_RED"
!!
!! OUTPUT
!!  Only writing.
!!
!! SIDE EFFECTS
!!  tensor(BSp%nomega,6) = Complex dielectric tensor to be written
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      c_f_pointer,etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine exc_write_tensor(BSp,BS_files,what,tensor)

 use defs_basis
 use m_bs_defs
 use m_errors
 use m_profiling_abi

 use m_io_tools,     only : open_file
 use m_fstrings,     only : toupper, strcat

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_write_tensor'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: what
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
!arrays
 complex(dpc),intent(in) :: tensor(BSp%nomega,6)

!Local variables ------------------------------
!scalars
 integer :: io,iq,funt
 real(dp) :: omegaev
!arrays
 character(len=500) :: lf_type,block_type,wgg_type,frm,str_type, msg
 character(len=fnlen) :: fname

!************************************************************************

 fname = strcat(BS_files%out_basename,'_',toupper(what))
 if (open_file(fname,msg,newunit=funt,form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 select case (toupper(what))
 case ("EXC_TSR_CART")

   write(funt,'("# Complex dielectric tensor (cart. coord.) obtained with the BS equation.")')

   lf_type = 'WITHOUT LOCAL FIELD EFFECTS'
   if (BSp%exchange_term>0) lf_type='LOCAL FIELD EFFECTS INCLUDED'
   call bsp_calctype2str(Bsp,str_type)
   write(funt,'("# ",a,"     " ,a)') TRIM(str_type), TRIM(lf_type)

   block_type = 'RESONANT-ONLY calculation'
   if (BSp%use_coupling>0) block_type = 'RESONANT+COUPLING calculation'
   write(funt,'("# ",a)') TRIM(block_type)

   if (BSp%use_coulomb_term) then
     wgg_type = "Coulomb term constructed with full W(G1,G2)"
     if ( BSp%use_diagonal_Wgg ) wgg_type = "Coulomb term constructed with diagonal approximation W(G1,G1)"
     write(funt,'("# ",a)') TRIM(wgg_type)
   end if

   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%soenergy*Ha_eV,' [eV]'
 
 case ("EXC_TSR_RED")

   write(funt,'("# Complex dielectric tensor (red. coord.) obtained with the BS equation.")')

   lf_type = 'WITHOUT LOCAL FIELD EFFECTS'
   if (BSp%exchange_term>0) lf_type='LOCAL FIELD EFFECTS INCLUDED'
   call bsp_calctype2str(Bsp,str_type)
   write(funt,'("# ",a,"     " ,a)') TRIM(str_type), TRIM(lf_type)

   block_type = 'RESONANT-ONLY calculation'
   if (BSp%use_coupling>0) block_type = 'RESONANT+COUPLING calculation'
   write(funt,'("# ",a)') TRIM(block_type)

   if (BSp%use_coulomb_term) then
     wgg_type = "Coulomb term constructed with full W(G1,G2)"
     if ( BSp%use_diagonal_Wgg ) wgg_type = "Coulomb term constructed with diagonal approximation W(G1,G1)"
     write(funt,'("# ",a)') TRIM(wgg_type)
   end if

   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%soenergy*Ha_eV,' [eV]'

 case ("RPA_NLF_TSR_CART")

   write(funt,'("# RPA complex dielectric tensor (cart. coord.) without local fields")')

 case ("RPA_NLF_TSR_RED")
   
   write(funt,'("# RPA complex dielectric tensor (red. coord.) without local fields")')

 case ("GW_NLF_TSR_CART")

   write(funt,'("# GW complex dielectric tensor (cart. coord.) without local field effects ")')
   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%soenergy*Ha_eV,' [eV]'

 case ("GW_NLF_TSR_RED")

   write(funt,'("# GW complex dielectric tensor (red. coord.) without local field effects ")')
   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%soenergy*Ha_eV,' [eV]'

 case default
   MSG_ERROR("Unknown value for what: "//TRIM(what))
 end select
 !
 ! Paramaters common to the different calculations.
 if (BSp%algorithm /= BSE_ALGO_HAYDOCK) then
   write(funt,'(a,i0)')"# nstates included in the diagonalization = ",BSp%nstates    
 end if
 
 if (BSp%algorithm == BSE_ALGO_HAYDOCK) then
   write(funt,'(a,2f7.4)')'# Tolerance = ',BSp%haydock_tol
 end if

 write(funt,'(a,i0)')"# npweps  = ",BSp%npweps
 write(funt,'(a,i0)')"# npwwfn  = ",BSp%npwwfn
 write(funt,'(a,i0)')"# nbands  = ",BSp%nbnds
 write(funt,'(a,i0)')"# loband  = ",BSp%lomo_spin(1)
 if (Bsp%nsppol==2) write(funt,'(a,i0)')"# loband(spin=2) = ",BSp%lomo_spin(2)
 write(funt,'(a,i0)')"# nkibz   = ",BSp%nkibz
 write(funt,'(a,i0)')"# nkbz    = ",BSp%nkbz
 write(funt,'(a,f7.4,a)')'# Lorentzian broadening = ',BSp%broad*Ha_eV,' [eV]'

 !
 ! Write tensor.
 write(funt,'(3a)') "# omega [eV] RE(eps_11) IM(eps_11) RE(eps_22)", &
&  "IM(eps_22) RE(eps_33) IM(eps_33) RE(eps_12) IM(eps_12)", &
&  "RE(eps_13) IM(eps_13) RE(eps_23) IM(eps_23))"
 write(frm,*) '(f7.3,12es14.6)'
 do io=1,BSp%nomega
   omegaev = DBLE(BSp%omega(io))*Ha_eV
   write(funt,frm) omegaev,(tensor(io,iq), iq=1,6)
 end do

 close(funt)

end subroutine exc_write_tensor
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/mdfs_ncwrite
!! NAME
!! mdfs_ncwrite
!!
!! FUNCTION
!!  Writes the MDF.nc file with the final results.
!!
!! INPUTS
!!  ncid =NC file handle
!!  Bsp<excparam>=Data type gathering the paramenters used for the Bethe-Salpeter calculation.
!!  eps_exc = Excitonic MDF
!!  eps_rpanlf = KS-RPA MDF without local-field effects.
!!  eps_gwnlf = GW-RPA MDF without local-field effects.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      exc_spectra,m_haydock
!!
!! CHILDREN
!!      c_f_pointer,etsf_io_low_def_var,etsf_io_low_set_define_mode
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!
!! SOURCE

subroutine mdfs_ncwrite(ncid,Bsp,eps_exc,eps_rpanlf,eps_gwnlf)

 use defs_basis
 use m_bs_defs
 use defs_datatypes   
 use defs_abitypes   
 use m_profiling_abi
 use m_errors
 use m_ncfile
 use iso_c_binding
#ifdef HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mdfs_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid 
 type(excparam),intent(in) :: BSp
!arrays
 complex(dpc),target,intent(in) :: eps_exc(BSp%nomega,BSp%nq)
 complex(dpc),target,intent(in) :: eps_rpanlf(BSp%nomega,BSp%nq)
 complex(dpc),target,intent(in) :: eps_gwnlf(BSp%nomega,BSp%nq)

!Local variables-------------------------------
!scalars
#ifdef HAVE_TRIO_ETSF_IO
 logical :: lstat
 type(ETSF_io_low_error) :: Error_data
 real(dp),ABI_CONTIGUOUS pointer :: real_ptr(:,:,:)

! *************************************************************************
 ! =========================
 ! === Write the dimensions 
 ! =========================
 call etsf_io_low_set_define_mode(ncid, lstat, Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'two',2,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'three',3,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_qpoints',Bsp%nq,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'number_of_frequencies',Bsp%nomega,lstat,Error_data=Error_data)  
 ETSF_CHECK_ERROR(lstat,Error_data)

! Define variables.
! scalars
! call etsf_io_low_def_var(ncid,'deltae',etsf_io_low_double,lstat,Error_data=Error_data)
! ETSF_CHECK_ERROR(lstat,Error_data)

!arrays
 call etsf_io_low_def_var(ncid,'qpoints',etsf_io_low_double,&
&  (/pad('three'),pad('number_of_qpoints')/),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'wmesh',etsf_io_low_double,(/'number_of_frequencies'/),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'exc_mdf',etsf_io_low_double,&
&  (/pad('two'), pad('number_of_frequencies'), pad('number_of_qpoints')/),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'rpanlf_mdf',etsf_io_low_double,&
&  (/pad('two'), pad('number_of_frequencies'), pad('number_of_qpoints')/),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'gwnlf_mdf',etsf_io_low_double,&
&  (/pad('two'), pad('number_of_frequencies'), pad('number_of_qpoints')/),lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

! Write data.
 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'qpoints',Bsp%q,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! Write frequency in mesh in eV.
 call etsf_io_low_write_var(ncid,'wmesh',REAL(Bsp%omega)*Ha_eV,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 ! Associate complex pointers with real inputs via the C pointers
 call C_F_pointer(C_loc(eps_exc), real_ptr, shape=[2,Bsp%nomega,BSp%nq])
 call etsf_io_low_write_var(ncid,'exc_mdf',real_ptr,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call C_F_pointer(C_loc(eps_rpanlf), real_ptr, shape=[2,Bsp%nomega,BSp%nq])
 call etsf_io_low_write_var(ncid,'rpanlf_mdf',real_ptr,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

 call C_F_pointer(C_loc(eps_gwnlf), real_ptr, shape=[2,Bsp%nomega,BSp%nq])
 call etsf_io_low_write_var(ncid,'gwnlf_mdf',real_ptr,lstat,Error_data=Error_data)
 ETSF_CHECK_ERROR(lstat,Error_data)

#else 
 MSG_ERROR("ETSF-IO support is not activated.")
#endif

end subroutine mdfs_ncwrite
!!***

