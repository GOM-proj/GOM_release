#=============================================================================#
# This makefile is originally written by Jungwoo Lee                           
# This makefile is for the release version.
#=============================================================================#
FC = gfortran
FCFLAGS = -O2
FCFLAGS += -fopenmp
#FCFLAGS_LINK = -fopenmp

# Activate debugging options -------------------------------------------------#
#FCFLAGS = -g

VPATH=../ 

objects = \
	mod_global_variables.o \
	mod_file_definition.o \
	mod_function_library.o \
	allocate_variables.o \
	ana_depth_adjustment.o \
	ana_windp.o \
	bilinear_interpolation.o \
	calculate_analytical_density.o \
	calculate_bottom_friction.o \
	calculate_coriolis_factor.o \
	calculate_density_full.o \
	calculate_density_linear.o \
	calculate_Q.o \
	calculate_reference_diff_days.o \
	calculate_velocity_at_node.o \
	calculate_voronoi_center.o \
	calculate_wind_stress.o \
	check_intersection.o \
	closing_message.o \
	compute_boundary_eta.o \
	compute_boundary_velocity.o \
	coordinate_conversion.o \
	create_regular_grid.o \
	ELM_backtrace_v1.o \
	find_element_and_face.o \
	find_hurricane_start_end_jday.o \
	find_neighbour_nodes.o \
	get_windp_at_face.o \
	holland_storm_surge.o \
	IDW_5.o \
	IDW_windp_ser.o \
	initialize_transport_model.o \
	julian.o \
	julian_to_localtime.o \
	model_termination_check.o \
	open_2D_out_files.o \
	open_3D_out_files.o \
	open_dump2D.o \
	open_dump3D.o \
	open_tser_out_files.o \
	open_tser_outfiles_tec.o \
	open_tser_outfiles_vtk.o \
	open_tser_outfiles_sub.o \
	pre_conj_grad.o \
	prepare_gom.o \
	read_bottom_roughness_inp.o \
	read_cell_inp.o \
	read_eta_ser.o \
	read_harmonic_ser.o \
	read_hurricane_ser_0.o \
	read_hurricane_ser_1.o \
	read_hurricane_ser_2.o \
	read_main_inp.o \
	read_node_inp.o \
	read_q_ser.o \
	read_salt_init.o \
	read_salt_ser.o \
	read_temp_init.o \
	read_temp_ser.o \
	read_windp_ser.o \
	scan_eta_ser.o \
	scan_gom.o \
	scan_harmonic_ser.o \
	scan_q_ser.o \
	scan_salt_ser.o \
	scan_temp_ser.o \
	scan_windp_ser.o \
	search_straight_line.o \
	seek_card.o \
	set_geometry_1.o \
	set_geometry_2.o \
	setup_cold_start.o \
	setup_hot_start.o \
	setup_spinup.o \
	skip_header_lines.o \
	solve_free_surface_equation.o \
	solve_momentum_equation.o \
	solve_nonlinear_advection.o \
	solve_transport_equation_v19.o \
	solve_velocities.o \
	string_to_lower.o \
	tridiagonal_solver.o \
	update_flood_map.o \
	update_variables.o \
	update_vertical_layers_0.o \
	update_vertical_layers_1.o \
	welcome_message.o \
	write_dump2D.o \
	write_dump3D.o \
	write_flood_map_tec.o \
	write_flood_map_vtk.o \
	write_grid_checking_files.o \
	write_output_files.o \
	write_restart.o \
	write_tec2D_body.o \
	write_tec2D_head.o \
	write_tec3D_full_body_zto_sigma.o \
	write_tec3D_full_head_zto_sigma.o \
	write_tec3D_surf_body.o \
	write_tec3D_surf_head.o \
	write_tser_out_files.o \
	write_vtk2D.o \
	write_vtk3D.o \
	zto_sigma_node_values.o \
	zto_sigma_cell_values.o \
	solve_transport_equation_ELM_v8.o \
	find_openboundary_salt_temp_v3.o \
	write_tec3D_full_head_z.o \
	write_tec3D_full_body_z.o \
	heat_term_by_term.o \
	heat_equilibrium_temperature.o \
	scan_air_ser.o \
	read_air_ser.o \
	calculate_vertically_averaged_values.o \
	calculate_Q_WR.o \
	calculate_Q_SS.o \
	main.o

all: $(objects)
	$(FC) $(FCFLAGS) ${FCFLAGS_LINK} $(objects) $(MY_LIBRARY_PATH)$(LIBS) -o run_release.exe 

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<	

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod *.MOD *.exe

# End of makefile ============================================================#
#=============================================================================#