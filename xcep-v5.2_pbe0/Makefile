FC=mpif90

#===================================================
#  compile with libxc 2.1.0   
#  for >2.2.0, need to modify the Makefile 
#===================================================

libxc_instpath=/gpfs/home/chuang3/XCEP-orb-free/libxc-2.1.0/install
fftw=/gpfs/home/chuang3/codes/fftw-3.3.4

libxc_lib=$(libxc_instpath)/lib/libxc.a

LIBXCMOD=$(libxc_instpath)/include  # where libxc modules are 

FFLAGS=-c -O2 -g -I$(LIBXCMOD)  -fbounds-check  # gfortran 
#FFLAGS=-c -O2 -g -I$(LIBXCMOD) -check bounds

#-I/panfs/panasas1/users/chuang3/work_fsu/codes/pfet/bspllib_22/

LDFLAGS=-L$(fftw)/install/lib  $(libxc_lib)  -lfftw3 -lm   -llapack -lblas

OBJECTS=interface_funcs.o \
    zdotc.o \
    dPackgmres.o \
    comm_data.o \
	read_system_zion.o \
	send_data_master_to_all.o \
    write_cluster_env_structure_XSF.o\
	drho2dphi.o\
	lbfgs_dir.o\
	force_penalty.o\
	apply_hess_vwreg1.o \
	compute_OF_dEdrho_vwreg1.o\
	compute_OF_dEdrho_vwreg3.o\
	compute_OF_dEdrho_vwreg2.o\
	compute_OF_dEdrho_complex_vwreg1.o\
	compute_OF_dEdrho_complex_vwreg3.o\
	compute_grad_vwreg1.o\
	compute_etotal_vwreg1.o\
	compute_etotal_vwreg3.o\
	compute_grad_complex_vwreg2.o\
	compute_grad_complex_vwreg3.o\
	compute_etotal_vwreg2.o\
	compute_OF_E.o \
	compute_OF_g.o \
	OF_fixpoint_solver.o\
	OF_drho_dN_vfix.o \
	OF_cg_scf.o \
	OF_cgmin.o \
	OF_cgmin_tot.o \
	OF_dfpt.o \
	OF_dfpt_inv.o \
	precond.o\
	precond_scf_residual.o \
	precond_teter.o\
	precond_kerker.o\
	precond_teter_inv.o\
	dotprod_four.o \
	load_xcart_cluster_env.o\
	update_coords_subsystem.o \
    print_newline.o \
	apply_Ainv.o \
    bspline90_22.o \
    entropy.o\
	make_points.o \
    dipole_moment.o \
	remove_mean.o\
	display_force.o\
	check_NaN.o\
	calc_delta_fermi.o \
	calc_dfermi_deigen.o \
	dfet_zp.o \
	get_reg_matrix.o \
	regularization.o \
	calc_dwatom_dR.o \
	system_mem_usage.o \
	sigma_vec.o\
	set_clu_env_q.o\
	new_occ_abinit_density.o \
	make_weights.o \
	make_weights_pbc.o \
	init_uemb.o \
	print_optimize_w_header.o \
	init_conv_radius.o \
	load_init_guess_vks.o \
	print_header_optimize_global_system.o \
	print_bfgs_conv_header.o \
	print_monitor_electron_numbers.o \
	calc_subsystem_dfpt.o \
	apply_chi.o \
	load_init_guess_refrho.o \
	check_file_length.o \
	run_subsystem_xcpatch.o \
	wait_for_system_end_xcpatch.o \
	get_subsystem_data_xcpatch.o \
	proj_dump_cluster_env_vks.o\
	fermi_cluster_env.o \
	run_global_ks.o \
	solve_z.o\
	solve_z_gmres.o\
	solve_z_lbfgs.o\
	solve_z_lbfgs_exlsr.o\
	solve_z_trans.o\
    solve_total_oep.o\
    solve_total_oep_regvxc.o\
    solve_total_oep_regvxc_trans.o\
	solve_total_oep_lbfgs.o\
	solve_total_oep_prec.o\
    check_subsystem_charges.o \
	file_len.o \
	dump_subsystem_vks.o \
	wrt_xsf.o \
	wrt_xsf_geo.o \
	pulay_mix.o\
	pulay_mix_scf.o\
	pulay_mix_vks_uemb.o \
	pulay_mix_kerker.o \
	pulay_mix_kerker_dfet.o \
	pulay_mix_gr.o\
	pulay_mix_zp_gr.o\
	pulay_mix_dfet_metric.o\
	pulay_mix_den_metric.o \
	load_psp_info.o\
	abinit_splfit.o \
	abinit_abi_derfc.o \
	abinit_ewald.o\
        abinit_spline.o  \
	abinit_ctrap.o \
	grand_chempot.o \
	load_atom_info.o \
        array_abinit2fftw3.o \
        get_time_date.o \
	make_subloc.o \
        cell2atom.o \
	atom2cell.o \
        welcome_message.o \
        hartree.o \
        hartree_screened.o \
        mathfunc.o \
		patch_xc.o \
        main.o \
        load_parameters.o \
        post_property.o \
        make_gprim.o \
        make_q_vector.o\
        Lbfgsb_chen.o \
        fft.o \
        display_info.o \
        check_bfgs_flag.o \
        det3.o \
        inverse3.o \
        round_off.o \
        debug_code.o \
        pw92lsda.o \
        tf.o \
        vw.o \
    	dump_uemb_xcpatch.o \
        gradient.o laplacian.o \
        calculate_xc_pbe.o\
		calculate_x_pbe_eps.o \
		calculate_c_pbe_eps.o \
        calculate_c_pbe.o\
        calculate_x_pbe.o\
        calculate_xc_pbe_xcpatch.o\
        calculate_xc_pw92LSDA_eps.o \
        calculate_xc_pw92LSDA.o \
        calculate_kedf.o \
        huang_carter.o \
        llp.o \
        check_conv.o \
        info_display.o \
        wrap_rhorsp.o


all: $(OBJECTS)
	$(FC) -o xcep $(OBJECTS) $(LDFLAGS) 

%.o : %.f90
	$(FC) $(FFLAGS) $< 


%.o : %.f
	$(FC) $(FFLAGS) $< 

clean: 
	rm -rf *.o  pfet  *.mod
