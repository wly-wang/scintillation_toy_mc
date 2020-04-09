// this file contains parameters used in simulation MERGED FILE
#include "TRandom3.h"
namespace parameters {
	
	// output file
	const char *output_file_name = "event_file_test_merge_Po210.root";

        ///-------------------------------------
        //--------WHAT events to generate?------
        ///-------------------------------------
        bool fixed_energy = false; //If gamma activity this will be set to true in the analyze_light.cpp
	double fixedE = -1.0; //Set to -1.0 to check if fixed value has been correctly assigned
        bool supernova = false;
        bool solar = false;
        bool gen_hep = false;
        bool gen_argon = false;
        bool gen_Co60B = false;
        bool gen_Co60G1 = false; //fixed energy = 1.173 Mev
        bool gen_Co60G2 = false; //fixed energy = 1.332 Mev
        bool gen_Ar42 = false;
        bool gen_K42 = false;
        bool gen_40KB = false;
        bool gen_40KG = false; //fixed energy = 1.46 MeV
        bool gen_Kr85B1 = false;
        bool gen_Kr85B2 = false;
        bool gen_Kr85G1 = false; //fixed energy = 0.151 Mev
        bool gen_Kr85G2 = false; //fixed energy = 0.305 Mev
        bool gen_radon = false;
	bool gen_Pb214 = false;
	//////-------------------------------///////
	//////--Radon alpha chain daughters--///////
	//////-------------------------------///////
	bool gen_Rn222 = false;
	bool gen_Po218 = false;
	bool gen_Po214 = false;
	bool gen_Po210 = true;
        bool gen_Bi214 = false;

        ///-------------------------------------
        ////----LAr & Ar39 properties--------------------
        /////-------------------------------------
        const double MassE = 0.510998910;       // mass electron - MeV/c^2
        const double Q_Ar = 0.565;                      //Q value of decay - Ar39
        const double Q_Co60B = 0.318;           //End point of beta decay of Co60[MeV]
        const double Q_Ar42 = 0.599;            //End point of beta decay of Ar42[MeV]
        const double Q_K42 = 3.525;
        const double Q_40KB = 1.35;             //End point of beta decay of 40K[MeV]
        const double Q_Kr85B1 = 0.84;           //End point of beta decay of Kr85[MeV]
        const double Q_Kr85B2 = 0.687;
        const double Q_Kr85G1 = 0.151;
        const double Q_Kr85G2 = 0.305;
        const double Q_Pb214 = 1.03;
        const double Q_Bi214 = 3.2;

        const double activity_Ar = 1.01;

        ///-------------------------------------
        ////----SN properties---------------------
        /////-------------------------------------
        const double Eav = 20.;                         // Average energy for SN spectrum
        const double expected_sn = 2.8;         // For poisson weighting
        ///-------------------------------------

        ///-------------------------------------
        ////----Radon properties---------------------
        /////-------------------------------------
        const int scint_yield_alpha = 16800;    // SY of alpha particles at 500 V/cm - from larsoft
        const double activity_Rn = 0.000021;    // Bq/kg
        const double massAlpha = 3727.3794;     // alpha particle mass - MeV/c^2
        const double Q_Rn222 = 5.490;                              // deposited energy from a radon decay - Rn-222 --> Po-218
        const double Q_Po218 = 6.0;
        const double Q_Po214 = 7.7;
        const double Q_Po210 = 5.3;
        ///-------------------------------------

	// position
	// random position ranges
	// Patrick Dimensions "Middle third" of a module
	const double entire_x_position_range[2] {5,350};	// cm
	const double entire_y_position_range[2] {-600,600};	// cm
	//const double entire_z_position_range[2] {400,1000};	// cm
        const double entire_z_position_range[2] {0,1400};   // cm

        //chrisflynn dimensions
        //whole volume
        //Ar39, Ar40, Bi214, K42, Kr85B1, Kr85B2, Kr85G1, Kr85G2, Pb214, Rn222
        /*const double entire_x_position_range[2] {10,363};      // cm
        const double entire_y_position_range[2] {-658,600};    // cm
        const double entire_z_position_range[2] {0,1400};    // cm
        */
        //40KB,40KG 
        const double K_x_position_range[2] {3.495e2,3.505e2};      // cm
        const double K_y_position_range[2] {-600,600};    // cm
        const double K_z_position_range[2] {0,1395};    // cm

        //Co60B
        const double Co_x_position_range[2] {-5e-1,5e-1};      // cm
        const double Co_y_position_range[2] {-600,600};    // cm
        const double Co_z_position_range[2] {0,1395};    // cm

        //Po210
        const double Po_x_position_range[2] {4.77e-1,1.477};      // cm
        const double Po_y_position_range[2] {-600,600};    // cm
        const double Po_z_position_range[2] {0,1395};    // cm

        

        //chrisflynn mass
        const double entire_mass = 0.001396*(entire_x_position_range[1]-entire_x_position_range[0])*(entire_y_position_range[1]-entire_y_position_range[0])*(entire_z_position_range[1]-entire_z_position_range[0]); //density = 1.396 g/cm^3

        //chrisflynn random variables
        const double frame_time = 0.0025; // 2.5 miliseconds; equivalent to the drift time of electrons in DUNE
        const double time_window = 1; //readout window (running time of the simulation in seconds)
        const double time_frames = int(time_window/frame_time); // number of frames registered
        


        ///-------------------------------------
        ////----Number of Events------------------
        /////-------------------------------------
        // Fixed energy (electron like) events:
        const int max_events_FE = 1;
	//const int number_events = max_events_FE;
	//const double energy = 25;     //MeV

        // Ar-39 events:
        //const int max_events_Ar = 10;
        const int max_events_Ar = activity_Ar * entire_mass * time_window;//FULL volume for 1 TPC
        const int Ar_decays_per_sec = activity_Ar* entire_mass; // decay rate in one TPC

        // Radon events:
        //const int max_events_Rn = 1000000;//To check the Fprompt, you need more Rn events, hope fully it would give you 6000000 photons
        const int max_events_Rn = 5.584e-5 * (entire_x_position_range[1]-entire_x_position_range[0])*(entire_y_position_range[1]-entire_y_position_range[0])*(entire_z_position_range[1]-entire_z_position_range[0]) * time_window * 0.01; //the 0.01 is included because of design changes to reduce alpha see https://indico.fnal.gov/event/20794/contribution/0/material/slides/0.pdf
        const int max_events_Rn222 = max_events_Rn * 0.25;
	const int max_events_Po218 = max_events_Rn * 0.25; //0.25 comes from 4 alpha decays in chain from Rn222, each contributes 10Bq/kg * 0.01 with new Rn222 condition
        const int max_events_Po214 = max_events_Rn * 0.25;
	//Will want a new variable for Po210 here if including events from Rn222 decay chain 
        const int max_events_Po210 = 5e-6 * (Po_x_position_range[1]-Po_x_position_range[0])*(Po_y_position_range[1]-Po_y_position_range[0])*(Po_z_position_range[1]-Po_z_position_range[0]) * time_window;//Activity from Po210 in the APAs
        const int max_events_Pb214 = max_events_Rn * 0.25;  
        const int max_events_Bi214 = max_events_Rn * 0.25;  
        //const int max_events_Rn = activity_Rn * mass * time_window;//Half volume for 1 (NOTE: for a small time window, this will probably return 0)
        const double Rn_decays_per_sec = activity_Rn* entire_mass; // decay rate in one TPC

        // Supernova events:
         const int max_events_SN = time_frames;
        //int max_events_SN = utility::poisson(expected_sn,gRandom->Uniform(1.),1.);

        // Solar neutrino events:
        const int max_events_SO = 10000;
        //const int max_events_SO = 0.00962286; //chrisflynn solarneutrino calculated for 10kt DUNE Module
        //int max_events_SO = utility::poisson(expected_sn,gRandom->Uniform(1.),1.);

        const int max_events_Co60B = 8.2e-5 * (Co_x_position_range[1]-Co_x_position_range[0])*(Co_y_position_range[1]-Co_y_position_range[0])*(Co_z_position_range[1]-Co_z_position_range[0]) * time_window;//For beta decay, near 100% would have beta decay and 200% gamma decay
        const int max_events_Co60G1 = max_events_Co60B;
        const int max_events_Co60G2 = max_events_Co60B;
        const int max_events_Ar42 = 1.283768e-7 * (entire_x_position_range[1]-entire_x_position_range[0])*(entire_y_position_range[1]-entire_y_position_range[0])*(entire_z_position_range[1]-entire_z_position_range[0]) * time_window;//1.41e-3 Bq*cm^-3 from mcc11 simulation
        const int max_events_K42 = max_events_Ar42 * 0.819;
        const int max_events_40KB = 2.7195e-3 * (K_x_position_range[1]-K_x_position_range[0])*(K_y_position_range[1]-K_y_position_range[0])*(K_z_position_range[1]-K_z_position_range[0]) * time_window * 0.8928;//89.28% ratio for beta decay
        const int max_events_40KG = 2.7195e-3 * (K_x_position_range[1]-K_x_position_range[0])*(K_y_position_range[1]-K_y_position_range[0])*(K_z_position_range[1]-K_z_position_range[0]) * time_window * 0.1072;//10.72% ratio for gamma decay
        const int max_events_Kr85B1 = 1.6e-4 * (entire_x_position_range[1]-entire_x_position_range[0])*(entire_y_position_range[1]-entire_y_position_range[0])*(entire_z_position_range[1]-entire_z_position_range[0]) * time_window * 0.785;//near 100% of Kr82 would have beta decay
        const int max_events_Kr85B2 = 1.6e-4 * (entire_x_position_range[1]-entire_x_position_range[0])*(entire_y_position_range[1]-entire_y_position_range[0])*(entire_z_position_range[1]-entire_z_position_range[0]) * time_window * 0.14;
        const int max_events_Kr85G1 = 1.6e-4 * (entire_x_position_range[1]-entire_x_position_range[0])*(entire_y_position_range[1]-entire_y_position_range[0])*(entire_z_position_range[1]-entire_z_position_range[0]) * time_window * 0.785 * 0.752;
        const int max_events_Kr85G2 = max_events_Kr85B2;
        const int max_events_hep = 100;





	// semi_analytic hits
	// set to false
	//const bool use_crowns_model = false;	// crowns model for visible hits, preliminary

	// timings
	const bool include_timings = true;
	const double timing_discretisation_step_size = 1.0;	// cm

	// visible light
	const bool include_reflected = true;


	// photon detection system properties
	const double quantum_efficiency = 0.035;				// arapuca QE
	const double wireplane_factor = 0.7;	
	const double vuv_transmission = 0.5;
	const double vis_transmission = 0.7;
 
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
