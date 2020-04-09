// implementation of utility functions class

#include "utility_functions.h"
#include <cmath>
#include "TMath.h" //chrisflynn
#include "TVector3.h" //chriflynn
#include "TF1.h" //chriflynn
#include "TRandom.h" //chriflynn
#include <chrono> //chriflynn
#include <iostream> //chriflynn
#include "TSpline.h" //chriflynn
#include "TFile.h" //chriflynn

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

	double Scintillation = exp(-(time/t_singlet))*singlet_part/t_singlet + exp(-(time/t_triplet))*triplet_part/t_triplet;

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

//Beta Decay Function ///chrisflynn
double utility_functions::SpectrumFunction(double *x, double *par)
{
       double KE = *x;
       double Q = par[0];
       double MassE = 0.510998910; //mass electron - Mev/c^2

       double N = std::sqrt(pow(KE,2) + 2*KE*MassE) * std::pow((Q-KE),2) * (KE+MassE);

       return N;
}

//Supernova neutrino energy spectrum ///chrisflynn
//Spectrum energy comes froma 1kpc galactic supernova
double utility_functions::fsn(double *x, double *par)
{
       double E = *x;
       double Eav = par[0];

       double f_nu = std::pow(E,3)*std::exp(-4*E/Eav);

       return f_nu;
}

//Solar neutrino energy spectrum ////chrisflynn
double utility_functions::fso(double *x, double *par)
{
       double E = *x;
       double Eav = par[0];

       TFile *f = new TFile("solar_neu_energy_spline.root");
       TSpline3 *spline = (TSpline3*)f->Get("Spline3");
       double f_s_neu = spline->Eval(E);
       f->Close();

       return f_s_neu;
}


//hep spectrum ///chrisflynn
double utility_functions::fhep(double *x, double *par)
{
      double E = *x;
      double Eav = par[0];

      TFile *f1 = new TFile("hep_neu_energy_spline.root");
      TSpline3 *spline = (TSpline3*)f1->Get("Spline3");
      double f_hep_neu = spline->Eval(E);
      f1->Close();

      return f_hep_neu;
}

//Radon-22 decay energy spectrum (Gaussian about the alpha decay energy) //chrisflynn
double utility_functions::Rn_function(double *x, double *par)
{
      double E = *x;
      double Q_Rn = par[0];

      double sigma = 0.01;
      //double sigma = (E-Q_Rn) / std::sgrt(1.3863); //1.3863 = ln(4)
      double sigma_sq = sigma * sigma;
      double gauss = 1/(sigma*std::sqrt(2*3.1416)) * std::exp((-1*std::pow((E-Q_Rn),2))/(2*sigma_sq));
      
      return gauss;
      }
      
