!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! Collection of Global Variables
!!
module mod_global_variables
	implicit none
	!! define precision control parameters ===================================!!
	integer, parameter :: sp = kind(0.0) 				! single precision
	integer, parameter :: dp = kind(0.d0) 				! double precision: this is for f90 version
	! integer, parameter :: dp = selected_real_kind(16) ! double precision: this is for f95 version
	
	!! === cpu time variables ================================================!!
	real(dp),save :: start, finish						! for cpu time
	real(dp),save :: elapsed_time = 0.0_dp 			! [sec], initialize to 0.0 second
	integer, save :: gom_start_year, gom_start_month, gom_start_day, gom_start_hour, gom_start_minute, gom_start_second			! for wall clock
	integer, save :: gom_finish_year, gom_finish_month, gom_finish_day, gom_finish_hour, gom_finish_minute, gom_finish_second	! for wall clock
	
	!! === Parameters ========================================================!!
   real(dp),parameter :: small_06	=	1.e-6
	
	! Calculation of pi()	
	! (1) this is the original method; use this approach.
   real(dp),parameter :: pi = 3.141592653589793_dp
   real(dp),parameter :: deg2rad = pi/180.0_dp		! degree to radian
   real(dp),parameter :: rad2deg = 180.0_dp/pi		! radian to degree
   
   ! (2) this is the second method; don't use this approach.
	! real(dp),parameter :: pi = 4*atan(1.0)
	! real(dp),parameter :: deg2rad = pi/180_dp

	! real(dp), parameter :: earth_omega = 7.292d-5	! Earth's angular velocity [rad/sec]: (2pi/[23h 56min 4.09sec])
	
   !! === main.inp variables ================================================!!
	!!========================================================================!!
	!! C0 ~ C19: General Model Configuration Information							  !!
	!!========================================================================!!

	! C00  TITLE FOR RUN ------------------------------------------------------!
	character(len=200), save :: project_name
	
	!! C01	Define general parameters ---------------------------------------!!
 	real(dp),save :: gravity									! gravitational acceleration, g [m/s2]
   real(dp),save :: rho_o										! reference density of seawater [kg/m3]
   real(dp),save :: rho_a										! reference density of air [kg/m3] @ 0 degree
   integer, save :: max_no_neighbor_node					! max allowed number of neighbor nodes and elements

	!! C1	Read grid information -----------------------------------------------!
	integer, save :: maxnod, maxele							! node number & cell number in a model grid
	real(dp),save :: h_node_adjust							! Additive adjustment to the node depth provided in "node.inp".
	integer, save :: ana_depth									! analytical depth adjustment on/off option
	integer, save :: coordinate_system						! 1: Cartician, 2: Lon/Lat
	real(dp),save :: lon, lat									! Longitude/Latitude of the model domain [degree]
	integer, save :: voronoi									! switch to whether Voronoi center will be use or not
	integer, save :: node_mirr, cell_mirr
	
	! additional variables:
	integer, save :: utm_projection_zone					! UTM project zone for the grid.
	real(dp),save :: lon_mid, lat_mid						! model grid's mid lon/lat, when using coordinate_system == 2 (i.e., Longitude/Latitude system)
	real(dp),save :: standard_meridian						! standard meridan of the model domain, this will be used if heat transport is activated
	
	!! C2 Read vertical layer information
	integer, save :: maxlayer									! Maximum vertical layers
	real(dp),save :: MSL											! Mean sea level from the arbitary bottom
	! real(dp),save :: z_level(0)								! this is defined at C2_1
	
	!! C2_1
   real(dp),save,allocatable,dimension(:) ::	&
   &			delta_z,									&			! (maxlayer); delta_z of each vertical layer; allocate_variables.f90, 
   &			z_level												! (0:maxlayer); level of vertical grid (each vertical layer's top face location), z_level(0) is the bottom elevation of the 1st z-layer; allocate_variables.f90
	
	! additional variables:
   real(dp),save :: delta_z_min								! minimum delta_z; it will be used in solve_nonlinear_advectioin.f90 (for ELM)
	
	!! C3 Read time related variables -----------------------------------------!
	real(dp),save :: dt											! simulation time step (delta t)
	integer, save :: ndt											! total number of simulation time steps
	integer, save :: data_start_year							! All data starts at data_start_year/January/1st 00:00:00
	integer, save :: data_time_shift							! All data will be sifted: data_start_year/January/1st 00:00:00 + data_time_shift*86400.0
   integer, save :: start_year, start_month, start_day, start_hour, start_minute, start_second	! simulation start year/month/day/hour/minute/second

	! additional variables:
   integer, save :: it = 0										! iteration number, set to zero at the start	
   real(dp),save :: current_jday_1900, air_p_julian_day_1900_2, wind_julian_day_1900_2   
   
   ! jday is the julian day from the beginning of the model simulation year (i.e., it is fixed), but
   ! julian_day is the time marching julian day (i.e., julian_day = jday + it*dt)
   ! i.e., jday is the actual julian day (which starts from "1") minus 1; and this is just one time value for the model simulation start day
   !       julian_day is the time marching from jday
   ! so, whenever need interpolation, I have to use "julian_day" not "jday"
   real(dp),save :: julian_day , jday, jday_1900		
   real(dp),save :: elapsed_time_hr
   integer, save :: year, month, day, hour, minute
   integer, save :: reference_diff_days
   
	!! C4	Restart & screen show control ---------------------------------------!
	integer, save :: restart_in, restart_out, restart_freq, ishow, ishow_frequency

   !! C5_01 Momentum equation variables I - Propagation/Nonlinear advection ---!
   real(dp),save :: theta
   integer, save :: advection_flag
   real(dp),save :: adv_onoff_depth
   integer, save :: ELM_backtrace_flag,	ELM_sub_iter, ELM_min_iter, ELM_max_iter

	!! C5_02	Momentum equation variables II - Bottom friction
   integer, save :: bf_flag
   integer, save :: bf_varying
   real(dp),save,allocatable,dimension(:) :: manning		! (maxface), Manning's n at face, allocate_variables.f90
   real(dp),save,allocatable,dimension(:) :: bf_height	! (maxface), bottom roughness height in [m] at face, allocate_variables.f90
   real(dp),save :: von_Karman									! von Karman constant for log law
	
	!! C5_03	Momentum equation variables III - Baroclinic gradient/Coriolis
   integer, save :: baroclinic_flag
   integer, save :: ana_density
   real(dp),save :: ref_salt, ref_temp					! reference salinity & temperature
   real(dp),save :: dry_depth
   integer, save :: Coriolis_option
   real(dp),save :: Coriolis								! constant coriolis
	! additional variables:
   real(dp),save,allocatable,dimension(:) :: 	&
   &			coriolis_factor								! (maxface), variable coriolis (it will be calculated) [1/s], allocated in: allocate_variables.f90
	
	!! C5_04 Momentum equation variables IV - Vertical & horizontal diffusion --!
   real(dp),save :: Ah_0									! Background horizontal eddy viscosity [m2/s]; constant add value
   real(dp),save :: Kh_0									! Background horizontal eddy diffusivity [m2/s]; constant add value
   real(dp),save :: Smagorinsky_parameter
	real(dp),save,allocatable,dimension(:,:) :: &	! horizontal diffusivity is defined at face
	&			Kh													! (maxlayer,maxface); horizontal diffusivity [m2/s]; allocate_variables.f90
	
	real(dp),save :: Av_0									! Background vertical eddy viscosity [m2/s]; constat add value
	real(dp),save :: Kv_0									! Background vertical eddy diffusivity [m2/s]; constant add value
   real(dp),save,allocatable,dimension(:,:) ::	&
   ! Note 1: Even though GOM uses constant vertical eddy viscosity,
   ! it is set as an arrayed variable to clarify the code.
   ! Note 2: Even though vertical eddy viscosity terms are defined [0:maxlayer,j] in vertical,
   ! the values at the top and the bottom (Av[maxlayer,j] and Av[0,j]) will not be used
   ! since they will be substituded with surface and bottom friction terms.
   ! However, I define these values as [0:maxlayer,j] to match with schematic diagram.
   ! Thus, these values' actual dimension, used in the code, is [1:maxlayer-1,maxface] 
   ! Note, Av (for momentuum equation) is defined at face center
   ! 		  Kv (for transport equation) is defined at element center
   &			Av													! (0:maxlayer,maxface); vertical eddy viscosity [m2/s]; allocate_variables.f90
   
	real(dp),save,allocatable,dimension(:,:) :: 	&
	&			Kv													! (0:maxlayer,maxele); vertical eddy diffusivity [m2/s]; allocate_variables.f90
! 	&			Kv_f,											&	! (0:maxlayer,maxface); vertical eddy diffusivity at face; allocate_variables.f90; I am not using this
! 	&			Kv_n												! (0:maxlayer,maxnod); vertical eddy diffusivity at node; allocate_variables.f90; I am not using this

	!! C6		Matrix solver option: OpneMP version of the Jacobi Preconditioned Conjugate Gradient method:
	integer, save :: max_iteration_pcg
	real(dp),save :: error_tolerance
	integer, save :: pcg_result_show



	!! C7		off/on tranport equations & selection of the solver -----------------!	
	integer,save :: transport_flag, transport_solver

	!! C7_1	TVD solver -------------------------------------------------------!
	integer,save :: trans_sub_iter, h_flux_limiter, v_flux_limiter
	
	!! C7_2	Turn off/on transport variables ----------------------------------!
	! in the future, I will increase these array as adding other trasportive materials
	! here, 2 = maxtran
	integer,save :: is_tran(2)=0
	integer,save :: num_tran_ser(2)=0
	integer,save :: tran_ser_shape(2)=0
	
	! additional variables
   integer,save :: maxtran = 2							! maximum transport variables allowed in GOM; currently I have two materials: salt & temp
	integer,save :: maxtran2								! number of active transport material; this will be used only for TVD algorithm
	integer,save,allocatable,dimension(:) :: &
	&			tran_id											! (maxtran2), scan_gom.f90
	
	! additional variables in salt_ser.inp
	integer,save :: is_salt, num_salt_ser
   integer,save :: salt_ser_shape						! salt_ser1.inp or salt_ser2.inp
	integer,save :: max_salt_ser, max_salt_data_num
	integer,save,allocatable,dimension(:) :: 	&
	&			salt_ser_data_num								! (max_salt_ser)
	
	! following variables are set to local variables:	
	! real(dp) :: salt_time_conv,	salt_time_adjust,	salt_unit_conv, salt_adjust
	
	real(dp),save,allocatable,dimension(:,:) :: 	&
	&			salt_ser_time,								&	! (max_salt_data_num,max_salt_ser)
	&			salt_ser_salt									! (max_salt_data_num,max_salt_ser)
	
	! additional variables for salt transport
	real(dp),save,allocatable,dimension(:,:) :: salt_at_obck	! (maxlayer,num_ob_cell)
	real(dp),save,allocatable,dimension(:) :: salt_at_Qbc 	! (num_Qb_cell)
	
	!! C7_3	Heat transport control variables 1 -------------------------------!
	integer,save :: heat_option, sol_swr
	real(dp),save:: fWz_a, fWz_b, fWz_c, wind_height
	
	!! C7_4 	Heat transport control variables 2 -------------------------------!
	real(dp),save:: light_extinction, sol_absorb, sed_water_exchange, T_sed, sed_temp_coeff
	
	! additional variable for heat transport
	real(dp),save,allocatable,dimension(:) :: &
	&			phi_n,  									& 	! (maxele), net surface heat flux
	&			phi_sw										! (maxele), sediment heat exchange at the bottom layer
	real(dp),save,allocatable,dimension(:,:):: &
	&			phi_sz										! (maxlayer,maxele), absorved short-wave solar radiation at each vertical layer
		
	! additional variables for air_ser.inp
	integer,save :: max_air_data_num
	real(dp),save,allocatable,dimension(:) ::	&
	&			air_ser_time,							&	! (max_air_data_num)
	&			T_air,									&	! (max_air_data_num)
	&			T_dew,									&	! (max_air_data_num)
	&			air_wind_speed,						&	! (max_air_data_num), wind_speed in air_ser.inp; this is to distinguish the general wind_speed from wind stress calculation
	&			solar_radiation,						&	! (max_air_data_num)
	&			cloud,									&	! (max_air_data_num)
	&			rain,										&	! (max_air_data_num)
	&			evaporation									! (max_air_data_num)
	
	! additional variables for temp transport
	real(dp),save,allocatable,dimension(:,:) :: temp_at_obck	! (maxlayer,num_ob_cell)
	real(dp),save,allocatable,dimension(:) :: temp_at_Qbc 	! (num_Qb_cell)

	! additional variables in temp_ser.inp
	integer,save :: is_temp, num_temp_ser
   integer,save :: temp_ser_shape						! temp_ser1.inp or temp_ser2.inp
	integer,save :: max_temp_ser, max_temp_data_num
	integer,save,allocatable,dimension(:) :: 	&
	&			temp_ser_data_num								! (max_temp_ser)
	
	real(dp),save,allocatable,dimension(:,:) :: 	&
	&			temp_ser_time,								&	! (max_temp_data_num,max_temp_ser)
	&			temp_ser_temp									! (max_temp_data_num,max_temp_ser)

   real(dp),save,allocatable,dimension(:) ::		&
   &			temp_ob											! (max_ob_element); interpolated temperature from temp_ser.inp; allocate_variables.f90

	
	!! C8		Model termination control ----------------------------------------!
	integer, save :: terminate_check_freq
	real(dp),save :: eta_min_terminate,	eta_max_terminate, uv_terminate, salt_terminate, temp_terminate
	
		      
	
   
	
	
	!!========================================================================!!
	!! C20 ~ C29, Tidal & River Boundary Conditions									  !!
	!!========================================================================!!
   ! C20 Surface elevation (tidal) open boundary information -----------------!
   integer, save :: num_ob_cell							! total number of open boundary cells = (num_tide_bc + num_eta_bc)
   integer, save :: ob_info_option						! open boundary information providing option
   integer, save :: num_tide_bc							! number of (harmonic) tidal boundary cells
   integer, save :: num_harmonic_ser					! number of harmonic tidal series
   integer, save :: check_tide							! 
   integer, save :: num_eta_bc							! number of surface elevation boundary cells
   integer, save :: num_eta_ser							! number of surface elevation time series
   integer, save :: check_eta
   integer, save :: eta_ser_shape						! eta_ser1.inp or eta_ser2.inp
   
   ! additional variables: ---------------------------------------------------!
   ! eta_ser.inp related variables:
   integer, save :: eta_ser_frequency
   integer, save :: max_eta_data_num, max_eta_ser
   real(dp),save,allocatable :: 						&
   &			eta_ser_time(:,:),						&	! (max_eta_data_num, max_eta_ser)
   &			eta_ser_eta(:,:)								! (max_eta_data_num, max_eta_ser
	
	integer,save,allocatable ::						&
	&			eta_ser_data_num(:)							! (max_eta_ser)
	
   real(dp),save,allocatable,dimension(:) ::		&
   &			eta_ob											! (max_ob_element); interpolated water surface elevation from eta_ser.inp; allocate_variables.f90
   
   real(dp),save,allocatable,dimension(:) ::		&	! (max_ob_element); interpolated water surface elevation from eta_ser.inp; allocate_variables.f90
   &			eta_at_ob_old,								&	! (max_ob_element); interpolated water surface elevation from eta_ser.inp; allocate_variables.f90
   &			eta_at_ob_new									
   
   ! C20_1 open boundary information -----------------------------------------!
	integer, save, allocatable, dimension(:,:) ::&
	&			ob_nodes 										! (max_ob_element,2)
   ! additional variables:
   integer, save, allocatable, dimension(:) :: 	&
   &			ob_cell_id,									&	! (max_ob_element)
   &			ob_face_id										! (max_ob_element)
	
	integer, save, allocatable, dimension(:) :: 	&
	&			ob_eta_type, 								&	! (max_ob_element); allocate_variables.f90
   &			harmonic_ser_id,							&	! (max_ob_element); allocate_variables.f90
   &			eta_ser_id,									&	! (max_ob_element); allocate_variables.f90 
   &			salt_ser_id,								&	! (max_ob_element); allocate_variables.f90
   &			temp_ser_id										! (max_ob_element); allocate_variables.f90
   
   
   real(dp),save :: eta_ser_dt

   	
   integer, save :: max_ob_node							! number of nodes which belongs to the open boundary
   integer, save :: max_ob_element						! number of elements which belongs to the open boundary
   integer, save :: max_ob_face							! number of facess which belongs to the open boundary (= max_ob_element)

   integer, save,allocatable :: 						&
   &			num_ob_element(:,:), 					&	! (max_ob_segment,max_ob_element); element ID numbers at each OB segment; allocate_variables.f90
   &			no_ob_elements(:)								! (max_ob_segment); element numbers at each OB segment; allocate_variables.f90
	
	! harmonic_ser.inp variables
   integer, save, allocatable, dimension(:) ::  &
   &			no_tidal_constituent							! (num_harmonic_ser); allocate_variables.f90
   integer, save :: max_no_tidal_constituent			! MAX(no_tidal_constituent)
   real(dp),save, allocatable, dimension(:) :: 	&
   &			tidal_phase_shift								! (num_harmonic_ser); allocate_variables.f90

   real(dp),save, allocatable, dimension(:,:)::	&
   &			tidal_amplitude,   						&	! (max_no_tidal_constituent,num_harmonic_ser); allocate_variables.f90
   &  		tidal_phase, 								&	! (max_no_tidal_constituent,num_harmonic_ser); allocate_variables.f90
   &			tidal_period,								&	! (max_no_tidal_constituent,num_harmonic_ser); allocate_variables.f90
   &			tidal_nodal_factor,						&	! (max_no_tidal_constituent,num_harmonic_ser); allocate_variables.f90
   &			equilibrium_argument							! (max_no_tidal_constituent,num_harmonic_ser); allocate_variables.f90


   ! C21 Setup spinup function for tide --------------------------------------!
   integer, save :: tide_spinup, baroclinic_spinup	
   real(dp),save :: tide_spinup_period, baroclinic_spinup_period

   
   ! C22 Include River Boundary ----------------------------------------------!
   integer, save :: num_Qb_cell, Qbc_info_option, num_Q_ser, check_Q, Q_ser_shape
   real(dp),save, allocatable, dimension(:) :: 	&
   &			Qu_boundary, 								&	! (max_Q_bc); allocate_variables.f90, face normal velocity [m/s] generated by river inflow
   &			Qv_boundary										! (max_Q_bc); allocate_variables.f90, face tangential velocity [m/s] generated by river inflow
   
   ! additional variables:
   integer, save :: max_Q_bc, max_Q_ser
   
   ! C22_1 -------------------------------------------------------------------!
	integer, save, allocatable :: 					&
	&			Q_nodes(:,:), 								&	! (max_Q_bc,2); allocate_variables.f90
	&			Q_ser_id(:),								&	! (max_Q_bc); allocate_variables.f90
	&			Q_salt_ser_id(:),							&	! (max_Q_bc); allocate_variables.f90
	&			Q_temp_ser_id(:)								! (max_Q_bc); allocate_variables.f90
	real(dp),save, allocatable :: 					&
	&			Q_portion(:)									! (max_Q_bc); allocate_variables.f90
   
   ! additional variables
   integer, save, allocatable :: 					&
   &			Q_boundary(:,:)								! (max_Q_bc,2); Q_boundary(i,1) = element_id, Q_boundary(i,2) = face_id; allocate_variables.f90
   real(dp),save, allocatable :: 					&
   &			Q_add(:)											! (max_Q_bc); allocate_variables.f90
   real(dp),save, allocatable :: 					&
   &			q_ser_time(:,:), 							&	! (max_q_data_num,max_Q_ser); allocate_variables.f90
   &			q_ser_Q(:,:)									! (max_q_data_num,max_Q_ser); allocate_variables.f90
   integer, save, allocatable :: 					&
   &			q_data_num(:)									! (max_Q_ser); allocate_variables.f90
   
   integer, save :: q_interp_method
   integer, save :: face_id, element_id
	integer, save :: max_q_data_num 						! maximum number of q time series data num
	
	! C23	Include Withdrawal/Return (WR) point boundary -----------------------!
	! integer, save :: num_WR_cell,	num_WR_ser, check_WR, WR_ser_shape, WR_layer
	integer, save :: num_WR_cell

	! additional variables
	integer, save :: max_WR_bc
	integer, save, allocatable ::	&
	&			WR_boundary(:,:),		&	! (max_WR_bc,2)
	&			isflowside4(:),		&	! (maxface)
	&			WR_element_flag(:)		! (maxele)
   real(dp),save, allocatable, dimension(:) :: 	&
   &			WRu_boundary, 								&	! (max_WR_bc); allocate_variables.f90, face normal velocity [m/s] generated by Withdrawl/Return
   &			WRv_boundary									! (max_WR_bc); allocate_variables.f90, face tangential velocity [m/s] generated by Withdrawl/Return
   real(dp),save, allocatable :: 					&
   &			Q_add_WR(:)										! (max_WR_bc); allocate_variables.f90

	
	! C23_1	Withdraw/Return (WR) point boundary information ----------------!
	integer, save, allocatable :: 	&
	&			WR_nodes(:,:),				&	! (max_WR_bc,2); allocate_variables.f90
	&			WR_layer(:),				&	! (max_WR_bc)
	&			WR_Q_ser_id(:),			&	! (max_WR_bc)
	&			WR_salt_ser_id(:),		&	! (max_WR_bc)
	&			WR_temp_ser_id(:)				! (max_WR_bc)
	real(dp),save, allocatable ::		&
	&			WR_portion(:)					! (max_WR_bc)
	
		
	! C24	Include Source/Sink (SS) point boundary -----------------------------!
	integer, save :: num_SS_cell

	! additional variables
	integer, save :: max_SS_bc
	integer, save, allocatable :: SS_element_flag(:)	! (maxele)
!   real(dp),save, allocatable, dimension(:) :: 	&
!   &			SSw_boundary 									! (max_SS_bc); allocate_variables.f90, w velocity [m/s] generated by Source/Sink
   real(dp),save, allocatable :: 						&
   &			Q_add_SS(:)											! (max_SS_bc); allocate_variables.f90

	! C24_1	Source/Sink (SS) point boundary information ----------------------!
	integer, save, allocatable ::		&
	&			SS_cell(:),					&	! (max_SS_bc)
	&			SS_layer(:),				&	! (max_SS_bc)
	&			SS_Q_ser_id(:),			&	! (max_SS_bc)
	&			SS_salt_ser_id(:),		&	! (max_SS_bc)
	&			SS_temp_ser_id(:)				! (max_SS_bc)
	real(dp),save, allocatable ::		&
	&			SS_portion(:)					! (max_SS_bc)
	
	
   ! C25 Include earth equilibrium tidal potential ---------------------------!
   integer, save :: no_Etide_species 
   real(dp),save :: Etide_cutoff_depth
   
   ! additional variables:
   integer, save :: max_Etide_species
   
   ! C25_1 -------------------------------------------------------------------!
   integer, save, allocatable :: 					&
   &			Etide_species(:) 								! (max_Etide_species); allocate_variables.f90
   real(dp),save, allocatable, dimension(:) :: 	&
   &			Etide_amplitude, 							&	! (max_Etide_species); allocate_variables.f90
   &			Etide_frequency, 							&	! (max_Etide_species); allocate_variables.f90
   &			Etide_nodal_factor, 						&	! (max_Etide_species); allocate_variables.f90
   &			Etide_astro_arg_degree						! (max_Etide_species); allocate_variables.f90

	! additional variables:
	real(dp),save,allocatable, dimension(:,:)::	&
	&			Etide_species_coef_at_node,			&	! (0:2,maxnod), allocate_variables.f90
	&			Etide_species_coef_at_element				! (0:2,maxele), allocate_variables.f90
	   
!=============================================================================!
! C23_1 is the end of the second set.														!
! C30 ~ C39, Air/Sea Boundary Conditions													!
!=============================================================================!
   ! C30 Wind & pressure data type #1 ----------------------------------------!
   integer, save :: wind_flag, airp_flag, num_windp_ser, wind_formular, wind_spinup
   real(dp),save :: wind_spinup_period
	
	! additional variables   
	! windp_ser.inp variables:
	integer, save :: max_windp_data_num					! maximum number of data in windp_ser.inp, scan_gom.f90
	integer, save :: max_windp_station					! maximum number of windp station (or number of windp series)
	! real(dp),save :: windp_time_conv					! Time unit conversion factor to [sec], this will set as a local variable in read_windp_ser.f90
	integer, save,allocatable,dimension(:) ::		&
	&			windp_ser_data_num,						&	! (max_windp_data_num); Total number of data points in each station; allocate_variables.f90
	&			windp_station_node							! (max_windp_station); closest node number for each wind station; allocate_variables.f90
	real(dp),save,allocatable,dimension(:,:):: 	&
	&			windp_ser_time,							&	! (max_windp_data_num,max_windp_station), windp time in [Julian day]
	&			windp_u,										&	! (max_windp_data_num,max_windp_station), Wind speed in x-direction [m/s]
	&			windp_v,										&	! (max_windp_data_num,max_windp_station), Wind speed in y-direction [m/s]
	&			windp_p											! (max_windp_data_num,max_windp_station), Pressure in [Millibar]
	! end of windp_ser.inp
	
	! windp_model.inp variables:
! 	real(dp),save :: windp_model_start_time			! wind model start time in [julian day]
! 	real(dp),save :: windp_model_end_time				! wind model end time in [Julian day]
! 	integer, save :: windp_model_delta_t				! wind model time increment in [sec]
! 	integer, save :: max_windp_model						! maximum wind model time series
	
! 	integer, save,allocatable,dimension(:)  :: 	&
! 	&	windp_model_time										! (max_windp_model)
	
! 	real(dp),save,allocatable,dimension(:,:):: 	&
! 	&	windp_model_u,										&	! (maxnod,max_windp_model); Wind speed in x-direction [m/s]
! 	&	windp_model_v,										&	! (maxnod,max_windp_model); Wind speed in y-direction [m/s]
! 	&	windp_model_p											! (maxnod,max_windp_model); Pressure in [Millibar]
	! end of windp_model.inp
	
	
   real(dp),save,allocatable,dimension(:)	::		&
   &			wind_u_at_node,							&	! (maxnod); wind speed in x-direction at each node [m/s]; allocate_variables.f90
   &			wind_v_at_node,							&	! (maxnod); wind speed in y-direction at each node [m/s]; allocate_variables.f90
   &			wind_u_at_face, 							&	! (maxface); wind speed in x-direction at each face [m/s]; allocate_variables.f90
   &			wind_v_at_face,							&	! (maxface); wind speed in y-direction at each face [m/s]; allocate_variables.f90
   &			wind_stress_normal,						&	! (maxface); normal wind stress at face; allocate_variables.f90
   &			wind_stress_tangnt							! (maxface); tangential wind stress at face; allocate_variables.f90
	
	real(dp),save,allocatable,dimension(:) ::		&
	&			airp_at_node,								&	! (maxnod); air pressure at node [Millibars]
	&			airp_at_face									! (maxface); air pressure at face [Millibars]
	   

   ! C31 Wind & pressure data type #2 ----------------------------------------!
   integer, save :: hurricane_flag															! Hurricane data set on/off option
   integer, save :: hurricane_data_type, hurricane_interp_method					! 
   integer, save :: hurricane_dt																! hurricane data time step in [sec]
	
	! C31_1 hurricane_ser.inp related variables -------------------------------!
	integer, save :: hurricane_start_year, hurricane_start_month, hurricane_start_day
	integer, save :: hurricane_end_year, hurricane_end_month, hurricane_end_day
	
	! additional variables:
	real(dp),save :: hurricane_start_jday, hurricane_end_jday
   
	! hurricane_ser.inp variables:
   real(dp),save,allocatable,dimension(:) :: 	&
   &			wind_u0, 									&	! (maxnod), first wind set from hurricane_ser.inp, allocate_variables.f90
   &			wind_u1, 									&	! (maxnod), interpolated value, allocate_variables.f90
   &			wind_u2,										&	! (maxnod), second wind set from hurricane_ser.inp, allocate_variables.f90
   &			wind_v0,										&	! (maxnod), first wind set from hurricane_ser.inp, allocate_variables.f90
   &			wind_v1, 									&	! (maxnod), interpolated value, allocate_variables.f90
   &			wind_v2,										&	! (maxnod), second wind set from hurricane_ser.inp, allocate_variables.f90
   &			air_p0, 										&	! (maxnod), first airp set from hurricane_ser.inp, allocate_variables.f90
   &			air_p1, 										&	! (maxnod), interpolated value, allocate_variables.f90
   &			air_p2											! (maxnod), second airp set from hurricane_ser.inp, allocate_variables.f90
	! additional variables for space interpolation: 
	! these are originally in read_hurricane_ser.f90, but moved here to avoid allocation/deallocation every time
   real(dp),save,allocatable,dimension(:) :: 	&	
   &			wind_u1_new, 								&	! (maxnod), space interpolated value for the 1st snapshot: wind_u, read_main_inp.f90
   &			wind_v1_new, 								&	! (maxnod), space interpolated value for the 1st snapshot: wind_v, read_main_inp.f90
   &			wind_u2_new, 								&	! (maxnod), space interpolated value for the 1st snapshot: wind_u, read_main_inp.f90
   &			wind_v2_new, 								&	! (maxnod), space interpolated value for the 1st snapshot: wind_v, read_main_inp.f90
   &			air_p1_new, 								&	! (maxnod), space interpolated value for the 1st snapshot: airp, read_main_inp.f90
   &			air_p2_new,									&	! (maxnod), space interpolated value for the 1st snapshot: airp, read_main_inp.f90
   &			shiftx, 										&	! (maxnod), space shifted node_x information from the 1st snapshot, read_main_inp.f90
   &			shifty											! (maxnod), space shifted node_y information from the 1st snapshot, read_main_inp.f90
	
	! hurricane_center.inp variables:
   integer, save :: &
   &			hurric_year_1, hurric_month_1, hurric_day_1, hurric_hour_1, hurric_minute_1, &
   &			hurric_year_2, hurric_month_2, hurric_day_2, hurric_hour_2, hurric_minute_2
   
   real(dp),save :: &
   &			hurric_x, hurric_y , hurric_latitude, 			&
   &  	   hurric_delta_pressure, hurric_mwr,				&
   &			hurric_x_1, hurric_y_1 , hurric_latitude_1,  &
   &  	   hurric_delta_pressure_1, hurric_mwr_1,			&
   &			hurric_x_2, hurric_y_2 , hurric_latitude_2,	&
   &  	   hurric_delta_pressure_2, hurric_mwr_2
   real(dp),save :: &
   &			hurric_julian_day_1     , hurric_julian_day_2,		&
   &  	   hurric_julian_day_1900_1, hurric_julian_day_1900_2	  
   real(dp),save :: hurric_coriolis, hurric_beta
   real(dp),save :: hurricane_time_1, hurricane_time_2
   real(dp),save :: hurricane_center_x, hurricane_center_y
	
	! additional variables:	
   integer, save :: hurricane_read_interval
	real(dp),save, allocatable :: regular_grid_xc(:), regular_grid_yc(:)
	integer, save, allocatable :: regular_grid_node_count(:,:), regular_grid(:,:,:)
	integer, save :: regular_grid_xi, regular_grid_yi ! total grid number on each axis
	real(dp),save:: regular_grid_half_dx, regular_grid_half_dy
	
	integer,save, allocatable :: hurricane_search_nodes(:,:), hurricane_search_count(:)

   ! C32 Wind & pressure data type #3 ----------------------------------------!
   integer, save :: holland_flag
   real(dp),save :: holland_start_jday, holland_end_jday

	! additional variables:
   real(dp),save,allocatable,dimension(:) :: 	&
   &			eta_from_Holland_at_ob						! (max_ob_element), allocate_variables.f90
   
!=============================================================================!
! C32 is the end of the third set.															!
! C40 ~ C49: Output Control Options															!
!=============================================================================!

   ! C40 Check grid ----------------------------------------------------------!
   integer, save :: check_grid_2D, check_grid_2DO, check_grid_3D, check_grid_format, check_grid_info
   real(dp),save :: check_grid_unit_conv

   ! additional variables:
	integer, save :: cell_connectivity_size_2D	! total number of data in 'CELLS' for vtk file format, for 2D
	integer, save :: cell_connectivity_size_3D	! total number of data in 'CELLS' for vtk file format, for 3D

   
   ! C41 ---------------------------------------------------------------------!
   integer, save :: tser_station_num, tser_info_option, tser_hloc, tser_time, tser_frequency, tser_format
   
   ! additional variables:
   real(dp),save :: tser_time_conv ! additiona, time conversion factor
   integer, save :: max_tser_station
	real(dp),save,allocatable,dimension(:) :: 	&
	&			eta_from_bottom								! (tser_station_num), allocated in: allocate_variables.f90
   
	! C41_1 -------------------------------------------------------------------!
   integer, save, allocatable, dimension(:) :: 	&
   &			tser_station_cell,						&	! (tser_station_num), allocate_variables.f90
   &			tser_station_node								! (tser_station_num), allocate_variables.f90
   character(25), save,allocatable, dimension(:) :: &
   &			tser_station_name								! (tser_station_num), allocate_variables.f90

   ! C41_2 -------------------------------------------------------------------!
   integer, save :: tser_eta, tser_H, tser_u, tser_v, tser_salt, tser_temp, tser_airp
   
	! C42 ---------------------------------------------------------------------!
	integer, save :: IS2D_switch, IS2D_time, IS2D_frequency, IS2D_start, IS2D_end, IS2D_format, IS2D_binary, IS2D_File_freq
	real(dp),save :: IS2D_unit_conv
	real(dp),save :: IS2D_time_conv 						! additional, time conversion factor 
	integer, save :: IS2D_File_num, IS2D_vtk_num
	character(len = 100), save :: IS2D_File_name
	integer, save :: IS2D_flood_map

	integer, save :: IS2D_buff = 0, IS2D_buff2 = 0 ! additional variables
	

	! flood map variables
	real(dp),save,allocatable,dimension(:) :: max_eta_node, max_flood_time	! (maxnod), allocate_variables.f90
	integer, save,allocatable,dimension(:) :: flood_id								! (maxnod), allocate_variables.f90

	! C42_1		2D output file (2D_tec.out) variable selection
	integer, save :: IS2D_variable(6) = 0
	integer, save :: zone_num_2D							! tecplot zone number
	
	
	! C43	3D output file (3D_tec.out) control; only for the surface elevation (eta)
	integer, save :: IS3D_surf_switch, IS3D_full_switch, IS3D_time, IS3D_frequency, IS3D_start, IS3D_end, &
	&					  IS3D_format, IS3D_grid_format, IS3D_binary, IS3D_File_freq
	real(dp),save :: IS3D_unit_conv
	real(dp),save :: IS3D_time_conv ! additional, time conversion factor 
	integer, save :: zone_num_3D_surf, zone_num_3D_full
	integer, save :: IS3D_File_num_surf, IS3D_File_num_full, IS3D_vtk_num
	character(len = 100), save :: IS3D_File_name_surf, IS3D_File_name_full

 	integer, save :: IS3D_buff = 0, IS3D_buff2 = 0 ! additional variables

	! C43_1		3D output file (3D_tec.out) variable selection
	integer, save :: IS3D_variable(6) = 0
	
	! C44	2D Dump option: data dump at horizontal 2D grid points
	integer, save :: IS2D_dump_switch, IS2D_dump_binary, IS2D_dump_time, IS2D_dump_File_freq, &
	&					  IS2D_dump_frequency, IS2D_dump_start, IS2D_dump_end
	integer, save :: IS2D_dump_File_num
	character(len = 100), save :: IS2D_dump_File_name
	real(dp),save :: IS2D_dump_time_conv 						! additional, time conversion factor 

 	integer, save :: IS2D_dump_buff = 0, IS2D_dump_buff2 = 0 ! additional variables


	! C45	3D Dump option
	integer, save :: IS3D_dump_switch, IS3D_dump_binary, IS3D_dump_time, IS3D_dump_File_freq, &
	&					  IS3D_dump_frequency, IS3D_dump_start, IS3D_dump_end
	integer, save :: IS3D_dump_File_num
	character(len = 100), save :: IS3D_dump_File_name
	real(dp),save :: IS3D_dump_time_conv 						! additional, time conversion factor 

 	integer, save :: IS3D_dump_buff = 0, IS3D_dump_buff2 = 0 ! additional variables
	
	! C46 NetCDF dump file
	integer, save :: netcdf_switch, netcdf_File_freq, netcdf_frequency, netcdf_start, netcdf_end
	integer, save :: netcdf_File_num
	character(len = 100), save :: netcdf_File_name

	integer, save :: netcdf_buff = 0, netcdf_buff2 = 0 ! additional variables
	integer, save :: netcdf_buff3 = 0 ! 
	integer, save :: netcdf_maxtime ! total data number (timeseries) in each netcdf output file
	
	! C46_1 NetCDF output variable selection at cell
	integer, save :: netcdf_variable_cell(6) = 0
	character(len=10), dimension(6) :: netcdf_variable_name_cell	! variables names are defined in allocate_variables.f90

	! C46_2	NetCDF output variable selection at node
	integer, save :: netcdf_variable_node(15) = 0
	character(len=10), dimension(15) :: netcdf_variable_name_node	! variables names are defined in allocate_variables.f90

	! C46_3	NetCDF output variable selection at face
	integer, save :: netcdf_variable_face(2) = 0
	character(len=10), dimension(2) :: netcdf_variable_name_face	! variables names are defined in allocate_variables.f90
	
	
	! C47 diagnostic files
	integer, save :: dia_advection, dia_momentum, dia_freesurface, dia_eta_at_ob
	integer, save :: dia_bottom_friction, dia_face_velocity, dia_node_velocity
   !=== End of main.inp variables ============================================!
   
   
	!==========================================================================!
   ! Additional variables ====================================================!
   real(dp),save :: 	&
   &			wtime1, wtime2, wtratio, windx1_s, windx2_s, windy1_s, windy2_s,	&
   &      	spinup_function_wind
      

   ! allocatable variables ===================================================!
   ! Velocity related variables ==============================================!
   real(dp),save,allocatable,dimension(:,:) ::	&
	! at face (at vertical layer): face normal velocities; these are essential velocities for momentum equation:
	! Note: I use un_, vn_ at face since they are face normal and tangent velocities not true u & v in xy coordinate
   &			un_ELM, 										&	! (maxlayer,maxface), face normal velocity at face for ELM; allocate_variables.f90
   &			vn_ELM,										&	! (maxlayer,maxface), face tangent velocity at face for ELM; allocate_variables.f90, 
   &			un_face,										&	! (maxlayer,maxface), normal velocity at face; allocate_variables.f90
   &			vn_face,										&	! (maxlayer,maxface), tangent velocity at face; allocate_variables.f90
   &			un_face_new,								&	! (maxlayer,maxface), normal velocity at face; allocate_variables.f90
   &			vn_face_new,								&	! (maxlayer,maxface), tangent velocity at face; allocate_variables.f90
   ! at cell (at vertical level): face normal velocities at prism's top and bottom faces; these are essential velocities for momentum equation:
   ! Note: I also use wn_ since they are face normal (in vertical coordinate) velocities and they are located at prism's top and bottom faces
   &			wn_cell,										&	! (0:maxlayer,maxele), velocity w at prism face, allocate_variables.f90
   &			wn_cell_new,								&	! (0:maxlayer,maxele), velocity w at prism face, allocate_variables.f90
   ! at node (at vertical level): these variables will be used when calculating ELM: ELM_backtrace.f90
	&			u_node,										&	! (0:maxlayer,maxnod), velocity u (true east/west) at node (at each vertical level); allocate_variables.f90
	&			v_node,										&	! (0:maxlayer,maxnod), velocity v (true north/south) at node (at each vertical level); allocate_variables.f90
	&			w_node											! (0:maxlayer,maxnod), velocity w at node (at each vertical level); allocate_variables.f90

	! at face: true face velocities at each vertical level: true u & v in xy coordinate
	! This velocity is to calculate nodal velocity and it is only used in calculate_velocity_at_node, i.e., it is a temporary velocities.
	! However, I set it as a global variable since it is a large array and thus allocation takes a long time when I put this in the subroutine
   real(dp),save,allocatable,dimension(:,:) :: 	&
   &			u_face_level,								&	! (0:maxlayer,maxface)
   &			v_face_level									! (0:maxlayer,maxface)
   
   real(dp),save,allocatable,dimension(:) :: 	&
   &			u_boundary, 								&	! (num_ob_cell); allocate_variables.f90
   &			v_boundary										! (num_ob_cell); allocate_variables.f90
   
   ! jw, I have to check this again.
   real(dp),save,allocatable,dimension(:,:) :: 	&
   &			velo_u_transport, 						&	! (0:maxlayer+1,maxnod), allocate_variables.f90
   &			velo_v_transport, 						&	! (0:maxlayer+1,maxnod), allocate_variables.f90
   &			velo_w_transport								! (0:maxlayer+1,maxnod), allocate_variables.f90
   
   ! 2D vertically averaged variables (only at node, not at cell)
   real(dp),save,allocatable,dimension(:) :: 	&
   &			ubar_node,									&	! (maxnod), allocate_variables.f90, vertically averaged velocity u (true east/wets) at node
   &			vbar_node,									&	! (maxnod), allocate_variables.f90, vertically averaged velocity v (true north/south) at node
   &			sbar_node,									&	! (maxnod), allocate_variables.f90, vertically averaged salinity at node
   &			tbar_node,									&	! (maxnod), allocate_variables.f90, vertically averaged temperature at node
   &			rbar_node										! (maxnod), allocate_variables.f90, vertically averaged water density (rho) at node
   ! End of velocity related variables =======================================!
   
   ! Geometry related variables ==============================================!
   integer, save :: maxface								! total number of faces, it will be calculated
   real(dp),save :: xn_min, yn_min 						! origin of model domain, minimums of x_node and y_node, for tecplot outputs
   integer, save :: start_end_node(4,4,3) 			! defines start and end node number at each face of an element, see more details in set_geometry_1.f90
   
   ! x, y, z coordinates -----------------------------------------------------!
   real(dp),save,allocatable,dimension(:) :: 	&
   ! at node
   &			x_node, 										&	! (maxnod), x coordinate of each node [m], scan_gom.f90 -> read_node_inp.f90,
   &			y_node, 										&	! (maxnod), y coordinate of each node [m], scan_gom.f90 -> read_node_inp.f90
   &			lon_node,									&	! (maxnod), longitude of each node [degree], scan_gom.f90 -> read_node_inp.f90
   &			lat_node,									&	! (maxnod), latitude of each node [degree], scan_gom.f90 -> read_node_inp.f90
   &			h_node, 										&	! (maxnod), depth at node [m], scan_gom.f90 -> read_node_inp.f90
   &			eta_node,									&	! (maxnod), water surface elevation at node [m], allocate_variables.f90
   ! at cell center
   &			x_cell,       								&	! (maxele), allocated in: allocate_variables.f90, x coordinate of each element [m]
   &			y_cell,       								&	! (maxele), allocated in: allocate_variables.f90, y coordinate of each element [m]
   &			lon_cell, 									&	! (maxele), allocated in: allocate_variables.f90, longitude at cell center [degree]
   &			lat_cell,									&	! (maxele), allocated in: allocate_variables.f90, latitude at cell center [degree]
	&			h_cell,										&	! (maxele), allocated in: allocate_variables.f90, depth at cell center [m]
	&			bed_elev,									&	! (maxele), allocated in: allocate_variables.f90, bed elevation - it is not yet used but I will keep this for the future.
   &			eta_cell,									&	! (maxele), allocated in: allocate_variables.f90, water surface elevation at cell center
	&			eta_cell_new,								&	! (maxele), allocated in: allocate_variables.f90, water surface elevation at cell center
	! at face center
   &			face_length,                 			&	! (maxface), allocated in: allocate_variables.f90
   &  		x_face, 										&	! (maxface), allocated in: allocate_variables.f90, x coordinate at face center [m]
   &			y_face,              					&	! (maxface), allocated in: allocate_variables.f90, y coordinate at face center [m]
   &			h_face											! (maxface), allocated in: allocate_variables.f90, depth at face center [m]
	
	! delta z -----------------------------------------------------------------!
	! dz_ is defined between lower and upper levels
	! dzhalf_ is defined between current and upper layer centers
	! Thus, dz_ will have    [1:maxlayer,  maxnod/maxele/maxface]
	!       dzhalf will have [1:maxlayer-1,maxnod/maxele/maxface]
	! However, I will define array sizes following layer and level notation for convenience:
	!       dz_					 [maxlayer,   maxnod/maxele/maxface]
	!       dzhalf_			 [0:maxlayer, maxnod/maxele/maxface]
	! Actually dz_node and dzhalf_node are not required for calculation, but just keep them.
	! dz_node is only required when calculating water pressure at nodal point in "calculate_density.f90", and it can be easily substituted with other calculation.
	! But, I just keep it. 
	real(dp),save,allocatable,dimension(:,:) ::	&
	! at node
	&			dz_node,										&	! (maxlayer,  maxnod), allocate_variables.f90, dz at each nodal point (at layer); it is required when calculating water pressure at nodal point in "calculate_density.f90"
	&			dzhalf_node,								&	! (0:maxlayer,maxnod), allocate_variables.f90, dz at each nodal point (between layers); actually, it is not required, but keep it.
	&			dz_node_new,								&	! (maxlayer,  maxnod)
	&			dzhalf_node_new,							&	! (0:maxlayer,maxnod)
   ! at cell
   &			dz_cell,										&	! (maxlayer,  maxele), allocate_variables.f90, dz at each cell (at layer)
   &			dzhalf_cell,								&	! (0:maxlayer,maxele), allocate_variables.f90, dz at each cell (between layers)
   &			dz_cell_new,								&	! (maxlayer,  maxele), this is required for solve_transport_equation.f90
   &			dzhalf_cell_new,							&	! (0:maxlayer,maxele), this is required for solve_transport_equation.f90
   ! at face
   &			dz_face,										&	! (maxlayer,  maxface), allocated in: allocate_variables.f90, dz_face is the distance between each level
   &			dzhalf_face,								&	! (0:maxlayer,maxface), allocated in: allocate_variables.f90, 
   &			dz_face_new,								&	! (maxlayer,  maxface), this is required for solve_transport_equation.f90
   &			dzhalf_face_new								! (0:maxlayer,maxface), this is required for solve_transport_equation.f90
   																! dzhalf_face is not dz_face/2 but dz_face(k+1,j)/2 + dz_face(k,j)/2, i.e., dz_face at each level
   																! Thus, it requires (maxlayer + 1) arrays
	! Top and Bottom layers ---------------------------------------------------!
   integer, save, allocatable, dimension(:) ::	&
   ! at node
   &			top_layer_at_node,						&	! (maxnod), allocated in: allocate_variables.f90, top layer number at each node, this will be updated at each time step
   &			bottom_layer_at_node,					&	! (maxnod), allocated in: allocate_variables.f90, bottom layer number at each node, this is fixed at the beginning
	! at cell
	&			top_layer_at_element, 					&	! (0:maxele), allocated in: allocate_variables.f90, include dummy room (0), top layer number at each element
   &			bottom_layer_at_element,				&	! (0:maxele), allocated in: allocate_variables.f90, include dummy room (0), bottom layer number at each element
	! at face
   &			top_layer_at_face, 						&	! (maxface), allocated in: allocate_variables.f90
   &			bottom_layer_at_face							! (maxface), allocated in: allocate_variables.f90

	! Adjacent (neighbour) nodes & cells & faces ------------------------------!
	! at node
   integer, save, allocatable, dimension(:) ::	&
   &			adj_cells_at_node,						&	! (maxnod), allocated in: scan_gom.f90 -> set_geometry_1.f90, total adjacent element numbers which sharing this node
   &			adj_nodes_at_node,						&	! (maxnod), allocated in: scan_gom.f90 -> set_geometry_1.f90, total adjacent node numbers at each node
   &			initial_wetdry_node,						&	! (maxnod), allocated in: scan_gom.f90 -> read_node_inp.f90, wet (=0) or dry (=1) index at each node, for initial condition
   &			wetdry_node 									! (maxnod), allocated in: allocate_variables.f90, wet (=0) or dry (= 1) index at each node	
	integer, save, allocatable, dimension(:,:)::	&  
	&			adj_cellnum_at_node,	   				&	! (max_no_neighbor_node,maxnod), allocated in: scan_gom.f90 -> set_geometry_1.f90, adjacent cell ID number at each node
	&			adj_nodenum_at_node,   					&	! (max_no_neighbor_node,maxnod), allocated in: scan_gom.f90 -> set_geometry_1.f90, adjacent node ID number at each node
	&			node_count_each_element						! (max_no_neighbor_node,maxnod), allocated in: scan_gom.f90 -> set_geometry_1.f90, the location of this node appears in elements which include this node (see more details in set_geometry_1.f90)
	! at cell
   integer, save, allocatable, dimension(:,:):: &
   &			nodenum_at_cell,							& 	! (4,maxele), allocated in: scan_gom.f90 -> read_cell_inp.f90, node number at each element (i.e., element constructing node IDs): 3 nodes in triangular cell and 4 nodes in quadrilateral cell
  	&			nodenum_at_cell_tec,						&	! (4,maxele), allocated in: scan_gom.f90 -> read_cell_inp.f90, node number at each element for tecplot
   &			adj_cellnum_at_cell,						&	! (4,maxele), allocated in: scan_gom.f90 -> set_geometry_1.f90, adjacent element ID number for face number
	&			facenum_at_cell								! (4,maxele), allocated in: allocate_variables.f90, face IDs at cell
	! at face
   integer, save,allocatable,dimension(:,:) :: 	&
   &			adj_cellnum_at_face, 					&	! (2,maxface), allocated in: set_geometry_1.f90, each face can have two adjacent cells (IDs) except at the boundary face
   &			nodenum_at_face								! (2,maxface), allocated in: set_geometry_1.f90, two node numbers (IDs) which consist each face 

	! Additional geometry -----------------------------------------------------!
	! at cell
   real(dp),save,allocatable,dimension(:) :: 	&
   &			area   	                      			! (maxele), allocated in: scan_gom.f90 -> read_cell_inp.f90, area [m^2] of each element
   integer, save, allocatable, dimension(:) ::	&
   &			tri_or_quad										! (maxele), allocated in: scan_gom.f90 -> read_cell_inp.f90, element geometry identifier (triangle = 3, rectangle = 4)
	integer,save,allocatable,dimension(:,:) ::	&
   &			sign_in_outflow								! (4,maxele), allocated in: set_geometry_1.f90, flow direction indicator at each element's face
	! at face
   real(dp),save,allocatable,dimension(:) :: 	&
   &			delta_j,   									&	! (maxface), allocated in: set_geometry_1.f90, distance between the centers of two adjacent elements at each face
   &			cos_theta,									&	! (maxface), allocated in: set_geometry_1.f90, cosine theta between each face and true xy coordinate
   &			sin_theta,									&	! (maxface), allocated in: set_geometry_1.f90, sine theta between each face and true xy coordinate
   &			cos_theta2,									&	! (maxface), allocated in: set_geometry_1,f90, angle difference from true face-normal to the fake face-normal
   &			sin_theta2										! (maxface), allocated in: set_geometry_1.f90, angle difference from true face-normal to the fake face-normal
   ! End of geometry related variables =======================================!

	! Open Boundary related variables =========================================!
	! at node
   integer, save, allocatable :: 					&
   &			num_serial_ob_node(:) 						! (max_ob_segment), allocated in: allocate_variables.f90, total number of open boundary nodes at each OB segment -1
! 	&  		ob_facenum(:,:)								! (max_ob_segment,max_ob_face), allocated in: allocate_variables.f90, face number (ID) at each open boundary segment
   ! at element
   integer,save,allocatable,dimension(:) ::		&
   &			ob_element_flag,							&	! (maxele), allocated in: allocate_variables.f90, ith ob element
   &			Qb_element_flag								! (maxele), allocated in: allocate_variables.f90, ith ob river
   ! at face
   integer,save,allocatable,dimension(:) :: 		&
   &			boundary_type_of_face,					&	! (maxface), allocated in: allocate_variables.f90
   &			isflowside,		 							&	! (maxface), allocated in: allocate_variables.f90, ith tidal flow boundary face
   &			isflowside2,								&	! (maxface), allocated in: allocate_variables.f90, ith radiation boundary face
   &			isflowside3										! (maxface), allocated in: allocate_variables.f90, ith river boundary (Qbc) face


   integer, save :: i_sponge_layer_flag
   real(dp),save :: spinup_function_tide, spinup_function_baroclinic
      
   integer, save, allocatable, dimension(:) ::  &
   &			i_tidal_boundary_temperature_type, 	&	! (max_ob_segment), allocated in: allocate_variables.f90
   &			i_tidal_boundary_salinity_type			! (max_ob_segment), allocated in: allocate_variables.f90

   character(len = 10), save, allocatable, dimension(:) :: 	&
   &			tidal_constituent_name_at_ob,			&	! (max_no_tidal_constituent), allocated in: allocate_variables.f90
   &  		tide_name										! (max_no_tidal_constituent), allocated in: allocate_variables.f90

   real(dp),save, allocatable, dimension(:) :: 	&
   &			tide_Q, 										& 	! (max_ob_segment), allocated in: allocate_variables.f90
   &			tth, 											& 	! (max_ob_segment), allocated in: allocate_variables.f90
   &			sth, 											& 	! (max_ob_segment), allocated in: allocate_variables.f90
   &			ath		  										! (max_ob_segment), allocated in: allocate_variables.f90   
	! End of open boundary related variables ==================================!


	! Momentum equation variables =============================================!   
   real(dp),save,allocatable,dimension(:,:) ::	&	! matrix indicies for momentum equations
   &			AinvDeltaZ1, 								&	! (maxlayer,maxface), allocated in: allocate_variables.f90, inv(A)dz_face for normal velocity
   &			AinvDeltaZ2,								&	! (maxlayer,maxface), allocated in: allocate_variables.f90, inv(A)dz_face for tangential velocity
   &			AinvG1, 										&	! (maxlayer,maxface), allocated in: allocate_variables.f90, G vector for normal velocity
   &			AinvG2											! (maxlayer,maxface), allocated in: allocate_variables.f90, G vector for tangent velocity

   

	
	! ELM variables ===========================================================!
	integer, save, allocatable, dimension(:,:)::	&  
	&			num_sub_elm_iteration						! (maxlayer,maxnod), allocated in: allocate_variables.f90
   
   real(dp),save :: denominator_min_for_matrix		! jw, check this number again

	! bottom_friction_variables ===============================================!
   real(dp),save,allocatable,dimension(:) :: 	&
   &			bottom_roughness, 						&	! (maxface), allocated in: allocate_variables.f90
   &			bottom_drag_coefficient,				&	! (maxface), allocated in: allocate_variables.f90
   &			Gamma_B,										&	! (maxface), allocated in: allocate_variables.f90
   &			Gamma_T,										&	! (maxface), allocated in: allocate_variables.f90 ! this is surface friction term (not bottom)
   &			Cdb												! (maxface), allocated in: allocate_variables.f90, bottom drag coefficient
   
   ! transport model variables ===============================================!
   ! Note that these variables at node and face are at vertically center position (not in each level position)
   ! That is why following variables only have maxlayer (not maxlayeer + 1) in vertical 
   real(dp),save,allocatable,dimension(:,:) :: 	&
   &			salt_cell,									&	! (maxlayer,maxele), salinity at cell center, allocated in: allocate_variables.f90
   &			salt_cell_new,								&	! (maxlayer,maxele), new salinity at cell center
   &			salt_node,									&	! (maxlayer,maxnod), salinity at horizontal nodal point at each veritcal layer, allocated in: allocate_variables.f90
   &			salt_face,									&	! (maxlayer,maxface), salinity at face at each vertical layer, allocated in: allocate_variables.f90
   &			temp_cell,									&	! (maxlayer,maxele), temperature at cell center, allocated in: allocate_variables.f90
   &			temp_cell_new,								&	! (maxlayer,maxele), new salinity at cell center
   &			temp_node,									&	! (maxlayer,maxnod), temperature at node, allocated in: allocate_variables.f90
   &			temp_face										! (maxlayer,maxface), temperature at face, allocated in: allocate_variables.f90
   
!   real(dp),save,allocatable,dimension(:,:,:) ::&	! currently, we have two transport variables: 1 = salinity, 2 = temperature
!   &			con_cell,									&	! (maxlayer,maxele,maxtran2), concentration at cell center, allocated in: allocate_variables.f90
!   &			con_cell_new,								&	! (maxlayer,maxele,maxtran2), new concentration at cell center
!   &			con_node_new,								&	! (maxlayer,maxnod,maxtran2), concentration at horizontal nodal point at each veritcal layer, allocated in: allocate_variables.f90
!   &			con_face_new,								&	! (maxlayer,maxface,maxtran2), concentration at face at each vertical layer, allocated in: allocate_variables.f90
!   &			r_jk,											&	! (maxlayer,maxface,maxtran2), horizontal consecutive gradient, allocated in: allocate_variables.f90
!   &			r_ik												! (0:maxlayer,maxele,maxtran2), vertical consecutive gradient, allocated in: allocate_variables.f90

	
   real(dp),save,allocatable,dimension(:) :: 	&
   &			air_temperature1, 						&	! (maxnod), allocated in: allocate_variables.f90
   &			air_temperature2,							&	! (maxnod), allocated in: allocate_variables.f90
   &			srad, 										&	! (maxnod), allocated in: allocate_variables.f90
   &			hradu, 										&	! (maxnod), allocated in: allocate_variables.f90
   &			hradd, 										&	! (maxnod), allocated in: allocate_variables.f90
   &			shum1, 										&	! (maxnod), allocated in: allocate_variables.f90
   &			shum2, 										&	! (maxnod), allocated in: allocate_variables.f90
   &			sflux, 										&	! (maxnod), allocated in: allocate_variables.f90
   &			fluxsu, 										&	! (maxnod), allocated in: allocate_variables.f90
   &			fluxlu 											! (maxnod), allocated in: allocate_variables.f90
   
      
   ! Densities are calculated at each horizontal points (node, face, and cell) at veritcal prism center; as other trasportive variables such as salinity and temperature.
   ! Thus, they are requied only [maxnod/maxele/maxface, maxlayer]
   real(dp),save,allocatable,dimension(:,:) :: 	&	! densities at face, cell, and node
   &			rho_face,									&	! (maxlayer,maxface), rho at face [kg/m3], allocate_variables.f90
   &			rho_cell,									&	! (maxlayer,maxele), rho at element center [kg/m3], allocate_variables.f90
   &			rho_node											! (maxlayer,maxnod), rho at node [kg/m3], allocate_variables.f90
   
   real(dp),save :: salinity_min, salinity_max, temperature_min, temperature_max
   integer, save :: i_density_flag


   character(len=2), save :: i_turbulence_model_name, stability_function 
   integer, save :: i_turbulence_flag
   real(dp),save :: qd, qd2, vd, td
      
   real(dp),save :: schk, schpsi
   real(dp),save :: cmiu0, cpsi1, cpsi2, rpub, rmub, rnub, psimin, eps_min, bgdiff, h1_pp, h2_pp,   &
   &       	tdmin_pp, vdmax_pp2, vdmin_pp2, vdmax_pp1, vdmin_pp1 , q2min 

   real(dp),save,allocatable,dimension(:) :: 	&
   &			sponge_relax, 								& 	! (maxnod), allocated in: allocate_variables.f90
   &			etaic												! (maxele), allocated in: allocate_variables.f90
   real(dp),save :: ttt, qq 
	
	! general_flags -----------------------------------------------------------!
   integer, save :: 										&
   &        heat_model_flag, 							&
   &        initial_temperature_field_flag, initial_salinity_field_flag,      &
   &        i_transport_model_flag
   integer, save :: ifile
	! End of global variables =================================================!
	!==========================================================================!
	


	! Transport variables for ELM ---------------------------------------------!
! 	real(dp),save,allocatable,dimension(:,:) ::	&
! 	&			tsd,											&	! (maxlayer,maxface)
! 	&			ssd,											&	! (maxlayer,maxface)
! 	&			tnd,											&	! (maxlayer,maxnod)
! 	&			snd,											&	! (maxlayer,maxnod)
	real(dp),save,allocatable,dimension(:,:) ::	&
	&			tem0,											&	! (maxlayer,maxnod)
	&			sal0,											&	! (maxlayer,maxnod)
	&			q2,											&	! (0:maxlayer,maxface)
	&			xl,											&	! (0:maxlayer,maxface)
! 	&			velo_u_transport,							&	! (0:maxlayer,maxnod)
! 	&			velo_v_transport,							&	! (0:maxlayer,maxnod)
! 	&			velo_w_transport,							&	! (0:maxlayer,maxnod)
	&			srho,											&	! (maxlayer,maxface)
	&			erho,											&	! (maxlayer,maxele)
	&			prho												! (maxlayer,maxnod)
	
end module mod_global_variables

