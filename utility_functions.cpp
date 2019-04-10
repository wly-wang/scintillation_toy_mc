// implementation of utility functions class

#include "utility_functions.h"

#include <cmath>

// destructor
utility_functions::~utility_functions(){
	//if fScintillation_function != nullptr {
		delete fScintillation_function;
	//}
}

// poisson distribution function
int utility_functions::poisson(const double mean, const double draw, const double eng) const {
	int number = 0;
	const int border = 16;
	double limit = 2e9;

	if(mean <= border) {
		double position = draw;
		double poissonValue = std::exp(-mean);
		double poissonSum = poissonValue;

		while(poissonSum <= position) {
			number++;
			poissonValue *= mean/number;
			poissonSum += poissonValue;
		}
		return number;
	} // the case of mean <= 16

	double value, t, y;
	t = std::sqrt(-2*std::log(draw));
	y = 2*3.141592654*eng;
	t *= std::cos(y);
	value = mean + t*std::sqrt(mean) + 0.5;
	if(value <= 0) {return 0; }
	if(value >= limit) { return limit; }
	return value;
}

// scintillation time spectrum function
double utility_functions::scintillation_function(const double *t, const double *par) {
	double time = *t;
	double t_singlet = par[0];
	double t_triplet = par[1];
	double type = par[2]; 		// type will be defined at 0 or 1, 0 is an electron, 1 is an alpha particle
	double singlet_part;
	double triplet_part;

	if(type == 0){ // particle is an electron
	  singlet_part = 0.30;
	  triplet_part = 0.70;
	}

	if(type == 1){ // particle is an alpha
	  singlet_part = 0.70;
	  triplet_part = 0.30;
	}

	double Scintillation = exp(-(time/t_singlet))*singlet_part/t_singlet;  + exp(-(time/t_triplet))*triplet_part/t_triplet;

	return Scintillation;
}

// function to create scintillation function TF1 with required parameters
void utility_functions::initalise_scintillation_function(const double t_singlet, const double t_triplet, const double scint_time_window, const double particle_type) {

	// create scintillation spectrum
	fScintillation_function = new TF1("Scintillation Timing", scintillation_function, 0, scint_time_window, 3);
	fScintillation_function->SetParameter(0, t_singlet); 
    fScintillation_function->SetParameter(1, t_triplet);  
    if(particle_type == 0) fScintillation_function->FixParameter(2, 0); 	// electron  
    if(particle_type == 1) fScintillation_function->FixParameter(2, 1);		// alpha
}