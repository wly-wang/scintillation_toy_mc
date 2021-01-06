// main driver for scintillation toy mc simulation code

#include<string>
#include<iostream>
#include<fstream>
#include<chrono>
#include <sstream>//necessary?
#include <vector>//necessary?
#include <algorithm>//necessary?
#include "TH1.h"

#include "TRandom.h"
#include "TVector3.h"

#include "data_output.h"
#include "semi_analytic_hits.h"
#include "time_parameterisation.h"
#include "utility_functions.h"


// include parameter file
#include "simulation_parameters.h"
using namespace std;

int main() {

	gRandom->SetSeed(0);

	// ------- Initialise output class ----------
	data_output output_file(parameters::output_file_name, parameters::include_timings, parameters::include_reflected);

	// -------- Initialise semi-analytic hits class ---------
	semi_analytic_hits hits_model;

	// -------- Initialise timing parametrisation class ---------
	time_parameterisation times_model(parameters::timing_discretisation_step_size);

	// -------- Initialise utility/energy spectrum class ---------
	utility_functions utility;
	utility.initalise_scintillation_function(parameters::t_singlet, parameters::t_triplet, parameters::scint_time_window, parameters::particle_type);
        utility.initalise_scintillation_function_alpha(parameters::t_singlet, parameters::t_triplet, parameters::scint_time_window, 1); //added to allow for alpha and gammas in same .root

	// ------- Read photon detector positions and types --------
	std::vector<std::vector<int>> opdet_type;
	std::vector<std::vector<double>> opdet_position;

	/////////////-------------INCLUDING BACKGROUND BOOLS---------------////////////////////

        double Q_beta_endpoint = 0;     
        if(parameters::gen_argon == true) {Q_beta_endpoint = parameters::Q_Ar;}
        if(parameters::gen_Co60B == true) {Q_beta_endpoint = parameters::Q_Co60B;}
        if(parameters::gen_Ar42 == true) {Q_beta_endpoint = parameters::Q_Ar42;}
        if(parameters::gen_K42 == true) {Q_beta_endpoint = parameters::Q_K42;}
        if(parameters::gen_40KB == true) {Q_beta_endpoint = parameters::Q_40KB;}
        if(parameters::gen_Kr85B1 == true) {Q_beta_endpoint = parameters::Q_Kr85B1;}
        if(parameters::gen_Kr85B2 == true) {Q_beta_endpoint = parameters::Q_Kr85B2;}
        if(parameters::gen_Pb214 == true) {Q_beta_endpoint = parameters::Q_Pb214;}
        if(parameters::gen_Bi214 == true) {Q_beta_endpoint = parameters::Q_Bi214;}

        if(parameters::gen_Co60G1 == true) {parameters::fixed_energy = true; parameters::fixedE = 1.173;}
        else if(parameters::gen_Co60G2 == true) {parameters::fixed_energy = true; parameters::fixedE = 1.332;}
        else if(parameters::gen_40KG == true) {parameters::fixed_energy = true; parameters::fixedE = 1.46;}
        else if(parameters::gen_Kr85G1 == true) {parameters::fixed_energy = true; parameters::fixedE = 0.151;}
        else if(parameters::gen_Kr85G2 == true) {parameters::fixed_energy = true; parameters::fixedE = 0.305;}


        TF1 *fSpectrum = new TF1("fSpectrum",utility_functions::SpectrumFunction,0,Q_beta_endpoint,1);
        TF1 *flandau_sn = new TF1("flandau_sn",utility_functions::fsn, 0, 50, 1);
        TF1 *flandau_so = new TF1("flandau_so",utility_functions::fso, 0, 16.56, 1); 
        TF1 *flandau_hep = new TF1("flandau_hep",utility_functions::fhep, 0, 18.78, 1); 

        flandau_sn->SetParameter(0, parameters::Eav);
        flandau_so->SetParameter(0, parameters::Eav);
        flandau_hep->SetParameter(0, parameters::Eav);

        TRandom3 *fGauss = new TRandom3();

        double a_gamma_distance;
	double gamma_length = 12.0;
	double added_x;
	double added_y;
	double added_z;

        TFile *f_alpha = new TFile("../spectra/real_alpha_positions_363.root");
	TH1D *alpha_gamma = (TH1D*)f_alpha->Get("histo");

	//////-----NOT SURE IF THESE VARIABLES ARE RELEVANT-----//////
        int max_events;
        int scint_yield;
        string particle;

        // GAMMA backgrounds

        if(parameters::fixed_energy == true) {
           max_events = parameters::max_events_FE;
        if(parameters::gen_Co60G1 == true) {max_events = parameters::max_events_Co60G1;}
        if(parameters::gen_Co60G2 == true) {max_events = parameters::max_events_Co60G2;}
        if(parameters::gen_40KG == true) {max_events = parameters::max_events_40KG;}
        if(parameters::gen_Kr85G1 == true) {max_events = parameters::max_events_Kr85G1;}
        if(parameters::gen_Kr85G2 == true) {max_events = parameters::max_events_Kr85G2;}

           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << " events, of fixed energy: " << parameters::fixedE << " MeV." << endl;
        }
        if(parameters::gen_argon == true) {
           fSpectrum->SetParameter(0, parameters::Q_Ar);
           //utility.SpectrumFunction->SetParameter(0, parameters::Q_Ar);//Setting the endpoint to Ar39 beta decay endpoint

           max_events = parameters::max_events_Ar;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", Ar 39 decays in time window: " << parameters::time_window << " seconds." << endl;
           cout << "This is equal to " << parameters::time_frames << " PMT readout frames." << endl;
           cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
           cout << "***NOTE. In ONE TPC, We expect to see: " << parameters::Ar_decays_per_sec << " Ar 39 decays each second.***" << endl;
           cout << "//////////////////////////////////////////////////////////////////////////////" << endl << endl;
        }
        if(parameters::gen_Co60B == true) {
           fSpectrum->SetParameter(0, parameters::Q_Co60B);
           //utility.SpectrumFunction(double *x, double *par)->SetParameter(0, parameters::Q_Co60B);//Setting the endpoint to Ar39 beta decay endpoint
           max_events = parameters::max_events_Co60B;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", Co60 beta decays in time window: " << parameters::time_window << " seconds." << endl;
           cout << "This is equal to " << parameters::time_frames << " PMT readout frames." << endl;
           cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
           cout << "//////////////////////////////////////////////////////////////////////////////" << endl << endl;
        }
        if(parameters::gen_Ar42 == true) {
           fSpectrum->SetParameter(0, parameters::Q_Ar42);
           //utility.SpectrumFunction(double *x, double *par)->SetParameter(0, parameters::Q_Ar42);//Setting the endpoint to Ar39 beta decay endpoint
           max_events = parameters::max_events_Ar42;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", Ar42 beta decays in time window: " << parameters::time_window << " seconds." << endl;
           cout << "This is equal to " << parameters::time_frames << " PMT readout frames." << endl;
           cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
           cout << "//////////////////////////////////////////////////////////////////////////////" << endl << endl;
        }
        if(parameters::gen_40KB == true) {
           fSpectrum->SetParameter(0, parameters::Q_40KB);
           max_events = parameters::max_events_40KB;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", 40K beta decays in time window: " << parameters::time_window << " seconds." << endl;
           cout << "This is equal to " << parameters::time_frames << " PMT readout frames." << endl;
           cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
           cout << "//////////////////////////////////////////////////////////////////////////////" << endl << endl;
        }
        if(parameters::gen_Kr85B1 == true) {
           fSpectrum->SetParameter(0, parameters::Q_Kr85B1);
           max_events = parameters::max_events_Kr85B1;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", Kr85B1 beta decays in time window: " << parameters::time_window << " seconds." << endl;
           cout << "This is equal to " << parameters::time_frames << " PMT readout frames." << endl;
           cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
           cout << "//////////////////////////////////////////////////////////////////////////////" << endl << endl;
        }
        if(parameters::gen_Kr85B2 == true) {
           fSpectrum->SetParameter(0, parameters::Q_Kr85B2);
           max_events = parameters::max_events_Kr85B2;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", Kr85B2 beta decays in time window: " << parameters::time_window << " seconds." << endl;
           cout << "This is equal to " << parameters::time_frames << " PMT readout frames." << endl;
           cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
           cout << "//////////////////////////////////////////////////////////////////////////////" << endl << endl;
        }
        if(parameters::supernova == true){
           max_events = parameters::max_events_SN;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", supernova events.\n";
        }
        if(parameters::solar == true){
           max_events = parameters::max_events_SO;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", solar events.\n";
        }
        if(parameters::gen_Pb214 == true){
           fSpectrum->SetParameter(0, parameters::Q_Pb214);
           max_events = parameters::max_events_Pb214;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", Pb214 events.\n";
        }
        if(parameters::gen_Bi214 == true){
           fSpectrum->SetParameter(0, parameters::Q_Bi214);
           max_events = parameters::max_events_Bi214;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", Bi214 events.\n";
        }
        if(parameters::gen_K42 == true){
           fSpectrum->SetParameter(0, parameters::Q_K42);
           max_events = parameters::max_events_K42;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", K42 events.\n";
        }
        if(parameters::gen_hep == true){
           max_events = parameters::max_events_hep;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", hep events.\n";
        }
        if(parameters::gen_radon == true){ //gen_radon == true
           max_events = parameters::max_events_Rn;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Radon 222 decays in time window: " << parameters::time_window << " seconds." << endl;
           cout << "This is equal to " << parameters::time_frames << " PMT readout frames." << endl;
           cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
           cout << "***NOTE. In ONE TPC, we expect to see: " << parameters::Rn_decays_per_sec << " radon-222 decays each second.***" << endl;
           cout << "/////////////////////////////////////////////////////////////////////////////" << endl << endl;
        }
        if(parameters::gen_Po210 == true){ //gen_radon == true
           max_events = parameters::max_events_Po210;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Po210 decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Rn222 == true){ //gen_radon == true
           max_events = parameters::max_events_Rn222;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Rn222 decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Po218 == true){ //gen_radon == true
           max_events = parameters::max_events_Po218;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Po218 decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Po214 == true){ //gen_radon == true
           max_events = parameters::max_events_Po214;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Po214 decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_alpha_gamma == true){ //gen_radon == true
           max_events = parameters::max_events_alpha_gamma;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Alpha-gamma decays in time window: " << parameters::time_window << " seconds." << endl;
        }



	///////----------INCLUDED THIS max_events variable TO KEEP CONSISTENT BETWEEN VERSIONS, can delete one-------////////
	//parameters::number_events = max_events;

	std::cout << "Loading Photon Detector positions..." << std::endl;
        std::ifstream detector_positions_file;
        detector_positions_file.open("optical_detectors_dune1x2x6.txt");
        if(detector_positions_file.is_open()) std::cout << "File opened successfully" << std::endl;
        else {std::cout << "File not found." << std::endl; exit(1);}
        while(!detector_positions_file.eof()) {
        int num_opdet, type_opdet; double x_opdet, y_opdet, z_opdet;
        if(detector_positions_file >> num_opdet >> x_opdet >> y_opdet >> z_opdet >> type_opdet) {
            std::vector<int> type({num_opdet, type_opdet});
            std::vector<double> position({x_opdet, y_opdet, z_opdet});         
            opdet_type.push_back(type);
            opdet_position.push_back(position);
        }
        else{ break; }
    }
    detector_positions_file.close();
    int number_opdets = opdet_type.size();
    std::cout << "Positions Loaded: " << number_opdets << " optical detectors." << std::endl << std::endl;

	// ----------- Create Events ------------
	// generate event positions and energies, storing information in output file	

	std::vector<double> energy_list; 
	//energy_list.reserve(parameters::number_events);
	energy_list.reserve(max_events);
	//std::vector<std::vector<double>> position_list(parameters::number_events, std::vector<double>(3,0.0)); 
	std::vector<std::vector<double>> position_list(max_events, std::vector<double>(3,0.0));

	std::cout << "Generating events..." << std::endl;
	for (int event = 0; event < max_events; event++){

		// output completion %
        if ( (event != 0) && (max_events >= 10) &&  (event % (max_events/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", event / (max_events/10));
        }

        // determine energy of event
        double energy;
        if(parameters::fixed_energy == true) {energy = parameters::fixedE;}
        if(parameters::gen_argon == true) {energy = fSpectrum->GetRandom();}    
        if(parameters::gen_Ar42 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_K42 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Bi214 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Pb214 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Co60B == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_40KB == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Kr85B1 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Kr85B2 == true) {energy = fSpectrum->GetRandom();}   
        if(parameters::supernova == true) {energy = flandau_sn->GetRandom();}   
        if(parameters::solar == true) {energy = flandau_so->GetRandom();}   
        if(parameters::gen_hep == true) {energy = flandau_hep->GetRandom();}
        if(parameters::gen_radon == true) {energy = fGauss->Gaus(parameters::Q_Rn222, 0.05);}  
        if(parameters::gen_Po210 == true) {energy = fGauss->Gaus(parameters::Q_Po210, 0.05);}
        if(parameters::gen_Rn222 == true) {energy = fGauss->Gaus(parameters::Q_Rn222, 0.05);}
        if(parameters::gen_Po218 == true) {energy = fGauss->Gaus(parameters::Q_Po218, 0.05);}
        if(parameters::gen_Po214 == true) {energy = fGauss->Gaus(parameters::Q_Po214, 0.05);}
	//alpha-gamma energy
	if(parameters::gen_alpha_gamma == true) {energy = fGauss->Gaus(parameters::Q_Rn222, 0.05);}
	    if(event % 2 ==0) { //Alpha event
		//scint_yield = parameters::scint_yield_alpha;
                energy = fGauss->Gaus(parameters::Q_Rn222, 0.05);
	    }
            else { //Gamma event
		//energy = fGauss->Gaus(15, 0.1);
		//scint_yield = parameters::scintillation_yield;
		energy = fGauss->Gaus(15.0, 2.9);
	    }
        energy_list.push_back(energy);

        // determine position of event
        // choose random position in 3D range
        if(parameters::fixed_energy == true || parameters::gen_argon == true || parameters::gen_Ar42 == true || parameters::gen_K42 == true || parameters::gen_Bi214 == true || parameters::gen_Pb214 == true || parameters::gen_Kr85B1 == true || parameters::gen_Kr85B2 == true || parameters::supernova == true || parameters::solar == true || parameters::gen_hep == true || parameters::gen_radon == true || parameters::gen_Kr85G1 == true || parameters::gen_Kr85G2 == true){
           position_list[event][0] = gRandom->Uniform(parameters::entire_x_position_range[0],parameters::entire_x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(parameters::entire_y_position_range[0],parameters::entire_y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(parameters::entire_z_position_range[0],parameters::entire_z_position_range[1]);
        }
        else if(parameters::gen_40KB == true || parameters::gen_40KG == true){
           position_list[event][0] = gRandom->Uniform(parameters::K_x_position_range[0],parameters::K_x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(parameters::K_y_position_range[0],parameters::K_y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(parameters::K_z_position_range[0],parameters::K_z_position_range[1]);
        }
        else if(parameters::gen_Co60B == true){
           position_list[event][0] = gRandom->Uniform(parameters::Co_x_position_range[0],parameters::Co_x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(parameters::Co_y_position_range[0],parameters::Co_y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(parameters::Co_z_position_range[0],parameters::Co_z_position_range[1]);
        }
        else if(parameters::gen_Po210 == true){
           position_list[event][0] = gRandom->Uniform(parameters::Po_x_position_range[0],parameters::Po_x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(parameters::Po_y_position_range[0],parameters::Po_y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(parameters::Po_z_position_range[0],parameters::Po_z_position_range[1]);
        }
        else {
           position_list[event][0] = gRandom->Uniform(parameters::entire_x_position_range[0],parameters::entire_x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(parameters::entire_y_position_range[0],parameters::entire_y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(parameters::entire_z_position_range[0],parameters::entire_z_position_range[1]);
        }

        //Event generation for alpha-gammas from Rn222 chain
        if(parameters::gen_alpha_gamma == true) {
	    if(event % 2 ==0) { //Alpha event
                position_list[event][0] = alpha_gamma->GetRandom();
                position_list[event][1] = gRandom->Uniform(parameters::entire_y_position_range[0],parameters::entire_y_position_range[1]);
                position_list[event][2] = gRandom->Uniform(parameters::entire_z_position_range[0],parameters::entire_z_position_range[1]);
	    }
	    else { //Gamma event
		a_gamma_distance = gRandom->Exp(gamma_length);
		TRandom obj;
		obj.SetSeed(0);
		obj.Sphere(added_x, added_y, added_z, a_gamma_distance);
		//std::cout << "r is " << a_gamma_distance << endl;
		//std::cout << "x is " << added_x << "  y is " << added_y << "  z is " << added_z << endl;
		position_list[event][0] = position_list[event-1][0] + added_x;
                position_list[event][1] = position_list[event-1][1] + added_y;
                position_list[event][2] = position_list[event-1][2] + added_z;
            }
	}
                //These statements are to stop gammas exiting through the edges of the TPC
                if (position_list[event][0] > 363.0) {
                    position_list[event][0] = 363.0;
                }
                if (position_list[event][0] < 1.0) {
                    position_list[event][0] = 1.0;
                }
                if (position_list[event][1] > 600.0) {
                    position_list[event][1] = 600.0;
                }
                if (position_list[event][1] < -600.0) {
                    position_list[event][1] = -600.0;
                }
                if (position_list[event][2] > 1400.0) {
                    position_list[event][2] = 1400.0;
                }
                if (position_list[event][2] < 0.0) {
                    position_list[event][2] = 0.0;
                }


	//---------SET THESE FOR FIXED POSITIONS----------//
        //position_list[event][0] = gRandom->Uniform(310.0,parameters::entire_x_position_range[1]);
        //position_list[event][0] = 5.0; //Fixed x position
        //position_list[event][1] = 0.0; //Fixed y position
        //position_list[event][2] = 200.0; //Fixed z position

        // add event properties to output file
        output_file.add_event(event, energy_list[event], position_list[event]);

	}//End of event loop

	std::cout << "Event generation complete." << std::endl << std::endl;
	
	// --------- Calculate hits and times ----------
	std::cout << "Determining number of photon hits..." << std::endl;

	//Efficiency testing
	/*
	std::chrono::steady_clock::time_point t_all_i; std::chrono::steady_clock::time_point t_all_f;
        std::chrono::steady_clock::time_point t_VUVTime_i; std::chrono::steady_clock::time_point t_VUVTime_f;
        std::chrono::steady_clock::time_point t_ReflTime_i; std::chrono::steady_clock::time_point t_ReflTime_f;
        std::chrono::steady_clock::time_point t_VUVHits_i; std::chrono::steady_clock::time_point t_VUVHits_f;
        std::chrono::steady_clock::time_point t_ReflHits_i; std::chrono::steady_clock::time_point t_ReflHits_f;
        std::chrono::duration<double> timespan_all (0.);
        std::chrono::duration<double> timespan_VUV_times (0.);
        std::chrono::duration<double> timespan_Refl_times (0.);
        std::chrono::duration<double> timespan_VUV_hits (0.);
        std::chrono::duration<double> timespan_Refl_hits (0.);

        t_all_i = std::chrono::steady_clock::now();
*/
	// loop each event in events list
	//for(int event = 0; event < parameters::number_events; event++) {
	for(int event = 0; event < max_events; event++) {

		// output completion %
        /*if ( (event != 0) && (parameters::number_events >= 10) &&  (event % (parameters::number_events/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", event / (parameters::number_events/10));
        }*/
        if ( (event != 0) && (max_events >= 10) &&  (event % (max_events/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", event / (max_events/10));
        }

        if(event % 2 ==0) { //Alpha event

        // calculate total scintillation yield from the event
        int number_photons = utility.poisson(static_cast<double>(parameters::scint_yield_alpha) * energy_list.at(event), gRandom->Uniform(1.), energy_list.at(event));

        // loop over each optical channel
        for(int op_channel = 0; op_channel < number_opdets; op_channel++) {

                // get optical detector type - rectangular or disk aperture
                int op_channel_type = opdet_type[op_channel][1];

                // get scintillation point and detection channel coordinates (in cm)
            TVector3 ScintPoint(position_list[event][0],position_list[event][1],position_list[event][2]);
            TVector3 OpDetPoint(opdet_position[op_channel][0],opdet_position[op_channel][1],opdet_position[op_channel][2]);

            // determine number of hits on optical channel via semi-analytic model:
            // VUV
            int num_VUV = 0;
            // incident photons
            int num_VUV_geo = hits_model.VUVHits(number_photons, ScintPoint, OpDetPoint, op_channel_type);       // calculate hits       
            // apply additional factors QE etc.            
            for(int i = 0; i < num_VUV_geo; i++) {
                if (gRandom->Uniform(1.) <= parameters::quantum_efficiency * parameters::wireplane_factor * parameters::vuv_transmission * (parameters::opdet_fraction_both + parameters::opdet_fraction_vuv_only)) num_VUV++;
            }
            // Visible
            int num_VIS = 0;
            if (parameters::include_reflected) {
                // incident photons
                int num_VIS_geo = 0;
                num_VIS_geo = hits_model.VisHits(number_photons, ScintPoint, OpDetPoint, op_channel_type);  // calculate hits with hotspot model        
                // apply additional factors QE etc.
                for(int j = 0; j < num_VIS_geo; j++) {
                        if (gRandom->Uniform(1.) <= parameters::quantum_efficiency * parameters::wireplane_factor * parameters::cathode_tpb_frac * parameters::vis_transmission * (parameters::opdet_fraction_both +  parameters::opdet_fraction_visible_only)) num_VIS++;
                }
            }

            // if no photons from this event for this optical channel, go to the next channel.
            if(num_VUV+num_VIS == 0) { continue; } // forces the next iteration

            // calculate timings
            std::vector<double> total_time_vuv; total_time_vuv.reserve(num_VUV);
            std::vector<double> total_time_vis; total_time_vis.reserve(num_VIS);
            if (parameters::include_timings){
                // VUV                  
                if(num_VUV > 0) {
                        // transport times
                        double distance_to_pmt = (OpDetPoint-ScintPoint).Mag();
                        double cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance_to_pmt;
                                        double theta = acos(cosine)*180./3.14159;
                                        int angle_bin = theta/45;       // 45 deg bins    

                        std::vector<double> transport_time_vuv = times_model.getVUVTime(distance_to_pmt, angle_bin, num_VUV);

                        // total times
                        for(auto& x: transport_time_vuv) {
                                double total_time = (x*0.001 + utility.get_scintillation_time_alpha()*1000000. + 2.5*0.001); // in microseconds
                                total_time_vuv.push_back(total_time);
                        }
                }
                // VIS
                if (num_VIS > 0 && parameters::include_reflected) {
                        // transport times
                        std::vector<double> transport_time_vis = times_model.getVisTime(ScintPoint, OpDetPoint, num_VIS);
                        // total times
                        for(auto& y: transport_time_vis) {
                                double total_time = (y*0.001 + utility.get_scintillation_time_alpha()*1000000. + 2.5*0.001); // in microseconds
                                total_time_vis.push_back(total_time);
                        }
                }
            }

            // fill data trees for each photon detected
            if (parameters::include_timings) output_file.add_data(event, op_channel, num_VUV, num_VIS, ScintPoint, total_time_vuv, total_time_vis);
            else output_file.add_data(event, op_channel, num_VUV, num_VIS, ScintPoint);

        } // end of optical channel loop

	} // end of alpha event

	else { //Gamma event
		
        // calculate total scintillation yield from the event
        int number_photons = utility.poisson(static_cast<double>(parameters::scintillation_yield) * energy_list.at(event), gRandom->Uniform(1.), energy_list.at(event));

        // loop over each optical channel
        for(int op_channel = 0; op_channel < number_opdets; op_channel++) { 

        	// get optical detector type - rectangular or disk aperture
        	int op_channel_type = opdet_type[op_channel][1];

        	// get scintillation point and detection channel coordinates (in cm)
            TVector3 ScintPoint(position_list[event][0],position_list[event][1],position_list[event][2]);
            TVector3 OpDetPoint(opdet_position[op_channel][0],opdet_position[op_channel][1],opdet_position[op_channel][2]);

            // determine number of hits on optical channel via semi-analytic model:
            // VUV
            int num_VUV = 0;
            // incident photons
            int num_VUV_geo = hits_model.VUVHits(number_photons, ScintPoint, OpDetPoint, op_channel_type);       // calculate hits       
            // apply additional factors QE etc.            
            for(int i = 0; i < num_VUV_geo; i++) {
                if (gRandom->Uniform(1.) <= parameters::quantum_efficiency * parameters::vuv_transmission * (parameters::opdet_fraction_both + parameters::opdet_fraction_vuv_only)) num_VUV++;   
            }
            // Visible
            int num_VIS = 0;
            if (parameters::include_reflected) {
            	// incident photons
            	int num_VIS_geo = 0;                
                num_VIS_geo = hits_model.VisHits(number_photons, ScintPoint, OpDetPoint, op_channel_type);  // calculate hits with hotspot model        
                // apply additional factors QE etc.
            	for(int j = 0; j < num_VIS_geo; j++) {
                	if (gRandom->Uniform(1.) <= parameters::quantum_efficiency * parameters::cathode_tpb_frac * parameters::vis_transmission * (parameters::opdet_fraction_both +  parameters::opdet_fraction_visible_only)) num_VIS++;
                }
            }

            // if no photons from this event for this optical channel, go to the next channel.
            if(num_VUV+num_VIS == 0) { continue; } // forces the next iteration

            // calculate timings
            std::vector<double> total_time_vuv; total_time_vuv.reserve(num_VUV);
            std::vector<double> total_time_vis; total_time_vis.reserve(num_VIS);
            if (parameters::include_timings){
            	// VUV            	
            	if(num_VUV > 0) {
            		// transport times
            		double distance_to_pmt = (OpDetPoint-ScintPoint).Mag();
            		double cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance_to_pmt;
  					double theta = acos(cosine)*180./3.14159;
  					int angle_bin = theta/45;	// 45 deg bins	  
            		
            		std::vector<double> transport_time_vuv = times_model.getVUVTime(distance_to_pmt, angle_bin, num_VUV);
            		
            		// total times
            		for(auto& x: transport_time_vuv) {
            			double total_time = (x*0.001 + utility.get_scintillation_time()*1000000. + 2.5*0.001); // in microseconds
            			total_time_vuv.push_back(total_time);
            		}
            	}
            	// VIS
            	if (num_VIS > 0 && parameters::include_reflected) {
            		// transport times
            		std::vector<double> transport_time_vis = times_model.getVisTime(ScintPoint, OpDetPoint, num_VIS);
            		// total times
            		for(auto& y: transport_time_vis) {
            			double total_time = (y*0.001 + utility.get_scintillation_time()*1000000. + 2.5*0.001); // in microseconds
            			total_time_vis.push_back(total_time);
            		}
            	}
            }

            // fill data trees for each photon detected
            if (parameters::include_timings) output_file.add_data(event, op_channel, num_VUV, num_VIS, ScintPoint, total_time_vuv, total_time_vis);
            else output_file.add_data(event, op_channel, num_VUV, num_VIS, ScintPoint);

        } // end of optical channel loop

        } //end of gamma event

	} // end of event loop

	f_alpha->Close();
	// write output root file
	output_file.write_output_file();

	std::cout << "Program finished." << std::endl;
}
