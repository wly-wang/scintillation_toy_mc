// main driver for scintillation toy mc simulation code

#include<string>
#include<iostream>
#include<fstream>
#include<chrono>

#include "TRandom.h"
#include "TVector3.h"

#include "data_output.h"
#include "semi_analytic_hits.h"
#include "time_parameterisation.h"
#include "utility_functions.h"


// include parameter file
#include "simulation_parameters.h"

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

	// ------- Read photon detector positions and types --------
	std::vector<std::vector<int>> opdet_type;
	std::vector<std::vector<double>> opdet_position;

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
	energy_list.reserve(parameters::number_events);
	std::vector<std::vector<double>> position_list(parameters::number_events, std::vector<double>(3,0.0)); 

	std::cout << "Generating events..." << std::endl;
	for (int event = 0; event < parameters::number_events; event++){

		// output completion %
        if ( (event != 0) && (parameters::number_events >= 10) &&  (event % (parameters::number_events/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", event / (parameters::number_events/10));
        }

        // determine energy of event
        energy_list.push_back(parameters::energy);

        // determine position of event
        // choose random position in 3D range
        position_list[event][0] = gRandom->Uniform(parameters::x_position_range[0],parameters::x_position_range[1]);
        position_list[event][1] = gRandom->Uniform(parameters::y_position_range[0],parameters::y_position_range[1]);
        position_list[event][2] = gRandom->Uniform(parameters::z_position_range[0],parameters::z_position_range[1]);

        // add event properties to output file
        output_file.add_event(event, energy_list[event], position_list[event]);
	}
	std::cout << "Event generation complete." << std::endl << std::endl;
	
	// --------- Calculate hits and times ----------
	std::cout << "Determining number of photon hits..." << std::endl;

	// loop each event in events list
	for(int event = 0; event < parameters::number_events; event++) {

		// output completion %
        if ( (event != 0) && (parameters::number_events >= 10) &&  (event % (parameters::number_events/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", event / (parameters::number_events/10));
        }

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
                if (gRandom->Uniform(1.) <= parameters::quantum_efficiency * parameters::mesh_factor * parameters::vuv_tpb_transmission * parameters::opdet_tpb_frac) num_VUV++;   
            }
            // Visible
            int num_VIS = 0;
            if (parameters::include_reflected) {
            	// incident photons
            	int num_VIS_geo = 0;                
                num_VIS_geo = hits_model.VisHits(number_photons, ScintPoint, OpDetPoint, op_channel_type);  // calculate hits with hotspot model        
                // apply additional factors QE etc.
            	for(int j = 0; j < num_VIS_geo; j++) {
                	if (gRandom->Uniform(1.) <= parameters::quantum_efficiency * parameters::mesh_factor * parameters::cathode_tpb_frac * (parameters::vis_tpb_transmission*parameters::opdet_tpb_frac + (1-parameters::opdet_tpb_frac))) num_VIS++;
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
            		std::vector<double> transport_time_vuv = times_model.getVUVTime(distance_to_pmt, num_VUV);
            		
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

	} // end of event loop

	// write output root file
	output_file.write_output_file();

	std::cout << "Program finished." << std::endl;
}