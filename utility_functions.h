#ifndef UTILITY_FUNCTION_H
#define UTILITY_FUNCTION_H

// class containing energy spectrums and other utility functions required by main code

#include "TF1.h"

class utility_functions {

private:
	TF1 *fScintillation_function = nullptr;

public:

	// constructor
	utility_functions(){};

	// destructor 
	~utility_functions();

	// poisson distribution function
	int poisson(const double mean, const double draw, const double eng) const;

	// scintillation function
	static double scintillation_function( const double *t, const double *par);
	void initalise_scintillation_function(const double t_singlet, const double t_triplet, const double scint_time_window, const double particle_type);
	double get_scintillation_time() { return fScintillation_function->GetRandom(); }

        // Spectrum Function (Beta Decay) /////chrisflynn
        static double SpectrumFunction(double *x, double *par);

        //Supernova Spectrum
        static double fsn(double *x, double *par);


        // HEP spectrum ////chrisflynn
        static double fhep(double *x, double *par);

        //Solar spectrum ////chrisflynn
        static double fso(double *x, double *par);

        //Rn Spectrum ////chrisflynn
        static double Rn_function(double *x, double *par);


};

#endif
