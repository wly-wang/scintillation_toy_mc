// this file contains parameters used in simulation

namespace parameters {
	
	// output file
	const char *output_file_name = "testing.root";

	// events to generate
	// number
	const int number_events =  5000;
	
	// type

	// energy
	const double energy = 25;	// MeV
	
	// position
	// random position ranges
	const double x_position_range[2] {10,350};	// cm
	const double y_position_range[2] {-500,500};	// cm
	const double z_position_range[2] {400,1000};	// cm


	// semi-analytic simulation options
	// visible light
	const bool include_reflected = true;

	// timings
	const bool include_timings = true;
	const double timing_discretisation_step_size = 1.0;	// cm


	// photon detection system properties
	const double quantum_efficiency = 0.035;				// arapuca QE
	const double wireplane_factor = 0.7;	
	const double vuv_transmission = 0.5;	// efficiency of any detector to VUV photons (128 or 174 nm)
	const double vis_transmission = 0.7; 	// efficiency of any detector to Visible photons
	
	const double opdet_fraction_vuv_only = 0.0;
	const double opdet_fraction_visible_only = 0.0; 
	const double opdet_fraction_both = 1.0 - opdet_fraction_vuv_only - opdet_fraction_visible_only;
	
	const double cathode_tpb_frac = 0.8;		


	// scintillation properties
	const int scintillation_yield = 24000; 		// 24000 photons/MeV at 500 V/m
	const double t_singlet = 0.000000006; 		// 6ns 
	const double t_triplet = 0.0000015; 		// 1.5 us
	const double scint_time_window = 0.00001; 	// 10 us
	const double particle_type = 0;			// ionising particle: 0 - electron, 1 - alpha
		
}
