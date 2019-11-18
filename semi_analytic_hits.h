#ifndef SOLID_ANGLE_H
#define SOLID_ANGLE_H

// class containing functions for calculating number of hits on each optical channel via the solid angle method
// calculates number of VUV hits, using gaisser-hillas corrections for Rayleigh scattering

// calculates number of visible hits via calculating number of hits on the cathode corrected with gaisser-hillas curves then the number
// of hits from the cathode for each optical channel using correction analogous to gaisser-hillas curves

#include <vector>
#include <string>

#include "TF1.h"
#include "TVector3.h"

class semi_analytic_hits {

private:

	// useful constants
	const double pi = 3.141592653589793;

	// detector type flag
	//const std::string flagDetector;

	// *************************************************************************************************
	//                    	OPTICAL DETECTOR SHAPE/SIZE
	// *************************************************************************************************	
	// Supercells: type = 1;
	double y_dimension_detector = 9.3;	// cm
	double z_dimension_detector = 46.8;	// cm
	// PMTs: type = 0; 
	double radius = 8*2.54/2.;	//	8" PMT diameter  // cm

	// structure definition for solid angle of rectangle function
	struct acc{
		// ax,ay,az = centre of rectangle; w = y dimension; h = z dimension
  		double ax, ay, az, w, h; 
	};

	bool _mathmore_loaded_ = false;


	// *************************************************************************************************
	//                    NUMBER OF VUV HITS PARAMETRIZATION
	// *************************************************************************************************
	// VUV hits Gaisser-Hillas Rayleigh scattering correction
	double fGHVUVPars[4][9];
	double fVUVBorderCorr[2];
	
	//  DUNE-SP Gaisser-Hillas							
	const double GH_RS60cm_SP[4][9] = { {1.37688, 1.35055, 1.29736, 1.21255, 1.0884, 0.985292, 0.828479, 0.686076, 0.531654}, 
										{105.721, 108.404, 114.289, 123.753, 142.937, 154.483, 179.057, 190.186, 188.029}, 
										{80, 80, 80, 80, 80, 80, 80, 80, 80}, 
										{-300, -300, -300, -300, -300, -300, -300, -300, -300} };

	// DUNE-SP border effects correction
	const double BORDER_RS60cm_SP[2] = {9.58514e-05, 0.0858229};

	// bin size of the offset angle used in the parametrization in degrees
	const double delta_angulo = 10.;
	// LAr absorption length in cm
	const double L_abs = 2000.;

	// active volume corner
	const double fYactive_corner = 1200/2; 	// cm
	const double fZactive_corner = 1400/2;	// cm
	double fReference_to_corner;


	// *************************************************************************************************
	//                    NUMBER OF VIS HITS PARAMETRIZATION
	// *************************************************************************************************
	// properties of LAr for reflected light path calculation
	// refractive indices in LAr
	const double n_LAr_VUV = 2.632;     // effective index due to group vel.
	const double n_LAr_vis = 1.23;

	// Detector properties: 
	
	// plane depth
	// Dune
	const double plane_depth = 363.38405;	// cm 
	
	// Cathode plane covered by foils
	// Dune
	// size
	const double y_dimension_foils = 1204.7255 + 5.466;	// cm		// 2 panels y height 602.36275 with 5.466cm gaps between
	const double z_dimension_foils = 1359.144 + 35.196;	// cm		// 6 panels of z width 226.524cm with 5.866cm gaps between them
	
	// centre coordinates
	// Dune
	const double x_foils = 363.38405; const double y_foils = 0; const double z_foils = 696.294;	// cm

	// Visible hits correction
	TF1* VIS_pol[9];
	const double delta_angle = 10.;

	// DUNE SP geometric correction [Arapuca/supercell, front window only]
	const double VIS_RS60cm_SP[6][9] = { {1.53686,1.46322,1.34083,1.14279,0.915688,0.648431,0.41744,0.41744,0.41744},
		{-0.0160284,-0.0145082,-0.0124793,-0.00975818,-0.0072178,-0.00464946,-0.00279529,-0.00279529,-0.00279529},
		{0.000112998,9.87981e-05,8.42566e-05,6.65798e-05,5.35068e-05,3.65524e-05,2.144e-05,2.144e-05,2.144e-05},
		{-4.75727e-07,-4.05306e-07,-3.49779e-07,-2.83337e-07,-2.52561e-07,-1.83123e-07,-1.04058e-07,-1.04058e-07,-1.04058e-07},
		{1.08826e-09,9.16466e-10,8.10421e-10,6.77421e-10,6.61672e-10,5.08153e-10,2.89974e-10,2.89974e-10,2.89974e-10},
		{-1.03393e-12,-8.73396e-13,-7.94945e-13,-6.86547e-13,-7.14112e-13,-5.73443e-13,-3.34853e-13,-3.34853e-13,-3.34853e-13} };
		
	/*
	// ** 7m drift **
	// RS60cm
	const double VIS_RS60cm_SP[6][9] = { {1.40583,1.15359,1.25194,1.54649,1.70694,1.99575,1.99575,1.99575,1.99575},
	{-0.0073011,-0.00344252,-0.00411537,-0.00622707,-0.00588044,-0.00674468,-0.00674468,-0.00674468,-0.00674468},
	{3.21944e-05,1.0809e-05,1.46532e-05,2.42112e-05,2.14304e-05,2.87278e-05,2.87278e-05,2.87278e-05,2.87278e-05},
	{-9.17456e-08,-3.63564e-08,-4.78928e-08,-7.18351e-08,-6.77508e-08,-9.91708e-08,-9.91708e-08,-9.91708e-08,-9.91708e-08},
	{1.3461e-10,6.7504e-11,8.3892e-11,1.15271e-10,1.1713e-10,1.72574e-10,1.72574e-10,1.72574e-10,1.72574e-10},
	{-7.49764e-14,-4.42768e-14,-5.28968e-14,-6.92522e-14,-7.38904e-14,-1.074e-13,-1.074e-13,-1.074e-13,-1.074e-13} };	
	*/

	// DUNE SP border correction [Arapuca/supercell, front window only]
	const std::vector<double> vDistances_x = {10, 40, 85, 130, 175, 220, 265, 310, 340, 355};			// cm	[10]
	const std::vector<double> vDistances_r = {122, 451, 647, 790};										// cm	[4]

	const std::vector<std::vector<std::vector<double>>> VIS_RS60_SP_Borders = { 
		{ {1.36687,1.22296,1.16675,1.14379,1.12301,1.10163,1.10099,1.09356,1.10187,1.12753},
		  {1.83254,1.61419,1.49286,1.43927,1.3912,1.35437,1.31543,1.27008,1.24884,1.25714},
		  {2.05931,1.85453,1.78544,1.79056,1.79707,1.8195,1.80777,1.7556,1.70645,1.70289},
		  {2.23968,2.15344,2.18786,2.36734,2.56407,2.75233,2.89681,2.90468,2.82028,2.79517}
		},
		{ {1.34339,1.20667,1.13506,1.11603,1.08991,1.07892,1.07119,1.06193,1.07188,1.09697},
		  {1.81072,1.59609,1.46559,1.4195,1.37173,1.32814,1.28887,1.2501,1.23295,1.24069},
		  {2.03827,1.82543,1.75276,1.74577,1.77848,1.78771,1.78426,1.73086,1.69444,1.68646},
		  {1.95378,1.84539,1.89026,2.02207,2.17497,2.33397,2.43329,2.42094,2.33427,2.31921}
		},
		{ {1.30212,1.1679,1.09268,1.06124,1.04809,1.02856,1.01896,1.00673,1.01583,1.03687},
		  {1.77623,1.53716,1.40967,1.36582,1.31276,1.28379,1.24822,1.21333,1.19201,1.19782},
		  {1.97939,1.75382,1.67763,1.70338,1.70174,1.71431,1.71404,1.67669,1.64612,1.64257},
		  {1.8617,1.72862,1.73262,1.86505,1.96976,2.09492,2.18476,2.16275,2.08759,2.06053}
		},
		{ {1.2217,1.09217,1.02034,1.0015,0.982541,0.969169,0.964904,0.960506,0.967135,0.985449},
		  {1.76632,1.55709,1.42782,1.38663,1.35201,1.32472,1.30128,1.2761,1.26203,1.26602},
		  {1.85359,1.6475,1.55838,1.56052,1.57037,1.57436,1.56395,1.53414,1.49935,1.49745},
		  {1.70973,1.60428,1.57908,1.68053,1.78143,1.87097,1.94755,1.94401,1.89018,1.874}
		},
		{ {1.19709,1.0421,0.975426,0.959783,0.943776,0.928057,0.925317,0.927927,0.939079,0.954253},
		  {1.80677,1.56675,1.45455,1.43069,1.40725,1.39217,1.38509,1.38628,1.39264,1.39857},
		  {1.78879,1.58522,1.50806,1.5289,1.54999,1.56731,1.58403,1.59104,1.58466,1.5856},
		  {1.51917,1.39067,1.3775,1.43522,1.52172,1.59911,1.65968,1.67423,1.63938,1.62868}
		},
		{ {1.31939,1.14616,1.06629,1.04966,1.03393,1.02082,1.0237,1.03905,1.05748,1.07774},
		  {1.44279,1.2267,1.11587,1.08317,1.0596,1.04128,1.03517,1.0428,1.05148,1.06037},
		  {1.49286,1.2971,1.20787,1.20624,1.20363,1.20469,1.20909,1.22173,1.2211,1.22693},
		  {1.35806,1.22853,1.19188,1.24884,1.30861,1.37408,1.43273,1.46065,1.45178,1.44975}
		},
		{ {1.67305,1.45578,1.35305,1.32855,1.32184,1.32226,1.33484,1.36954,1.40883,1.44865},
		  {1.539,1.30279,1.15584,1.11571,1.09565,1.08999,1.09669,1.11903,1.14344,1.16256},
		  {1.41221,1.20621,1.09857,1.08457,1.09274,1.10523,1.12587,1.15326,1.17331,1.18866},
		  {1.12174,0.991628,0.935007,0.961442,1.01713,1.07786,1.13943,1.18148,1.19394,1.2022}
		},
		{ {1.67305,1.45578,1.35305,1.32855,1.32184,1.32226,1.33484,1.36954,1.40883,1.44865},
		  {1.37135,1.13834,0.992489,0.951638,0.935147,0.933131,0.942681,0.970196,0.993699,1.01489},
		  {1.15366,0.96338,0.857661,0.841519,0.846182,0.86541,0.893007,0.92589,0.949376,0.968278},
		  {0.884385,0.759848,0.707666,0.724935,0.767402,0.824071,0.883741,0.935287,0.958996,0.97325}
		},
		{ {1.67305,1.45578,1.35305,1.32855,1.32184,1.32226,1.33484,1.36954,1.40883,1.44865},
		  {1.37135,1.13834,0.992489,0.951638,0.935147,0.933131,0.942681,0.970196,0.993699,1.01489},
		  {1.15366,0.96338,0.857661,0.841519,0.846182,0.86541,0.893007,0.92589,0.949376,0.968278},
		  {0.884385,0.759848,0.707666,0.724935,0.767402,0.824071,0.883741,0.935287,0.958996,0.97325}
		}
	};


public:	
	// constructor 
	semi_analytic_hits();

	// destructor
	~semi_analytic_hits(){};

	// hits calculating functions
	int VUVHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type);
	int VisHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type);

	// gaisser-hillas function
	static Double_t GaisserHillas(double x, double *par);

	// solid angle of rectangular aperture calculation functions
	double omega(const double &a, const double &b, const double &d) const;
	double solid(const acc& out, const TVector3 &v) const;

	// solid angle of circular aperture calculation functions
	double Disk_SolidAngle(double *x, double *p);
	double Disk_SolidAngle(double d, double h, double b);

	// linear interpolation function
	double interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate );

};

#endif