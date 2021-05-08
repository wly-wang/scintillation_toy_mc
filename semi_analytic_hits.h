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
	// Note: wires are now incorporated into the parameterisations, no longer need to add scaling factor
	
	// DUNE-SP Gaisser-Hillas
	// Argon, RS = 99.9cm, flat PDs (Arapucas/Supercells)
	const double fGHVUVPars_flat[4][9] = { {1.23286, 1.20491, 1.1597, 1.08148, 0.986607, 0.868075, 0.725637, 0.633551, 0.469717},
  								{150.325, 150.109, 156.294, 169.271, 179.642, 212.9, 232.173, 226.557, 215.166},
  								{19.0558, 18.8478, 18.8757, 18.7149, 50.8662, 49.2247, 57.6717, 157.92, 172.494},
  								{-3000, -3000, -3000, -3000, -1000, -1000, -1000, -100, -100} };
	std::vector<double> angulo_flat = {5, 15, 25, 35, 45, 55, 65, 75, 85};
	std::vector<double> slopes1_flat = {-9.63855e-05, -6.82604e-05, -9.63478e-05, -0.000121181, -0.000126611, -0.000115481, -8.61492e-05, -0.000112594, -7.80935e-05};
	std::vector<double> slopes2_flat = {-0.0662469, -0.0504497, -0.0596321, -0.0418021, -0.0342462, -0.0531668, -0.0522639, -0.0578887, -0.0591081};
	std::vector<double> slopes3_flat = {-0.00593207, -0.00672713, -0.0020843, -0.00216374, 0.00901291, 0.00385402, 0.0066081, 0.0341547, 0.0446519};
    
    // Xenon, flat PDs (Arapucas/Supercells)
	/*							  
    const double fGHVUVPars_flat[4][9] = { {0.959653, 0.92975, 0.874812, 0.790038, 0.685581, 0.566253, 0.433698, 0.31073, 0.197959},
  								{94.2397, 94.5496, 104.245, 139.412, 189.763, 192.723, 340.598, 501.523, 503.954},
  								{385.308, 413.922, 405.968, 371.839, 330.268, 484.33, 513.097, 512.926, 253.092},
  								{-400, -400, -400, -400, -400, -400, -400, -400, -400} };
	std::vector<double> angulo_flat = {5, 15, 25, 35, 45, 55, 65, 75, 85};
	std::vector<double> slopes1_flat = {0.000212425, 0.000200659, 0.000179874, 0.000165088, 0.00015168, 0.000122913, 8.94813e-05, 3.98645e-05, 1.24325e-13};
	std::vector<double> slopes2_flat = {0.151348, 0.168286, 0.153929, 0.12419, 0.0834517, 0.0955464, -0.0175567, -0.112032, -1.98752e-10};
	std::vector<double> slopes3_flat = {-0.114771, -0.19272, -0.182644, -0.130611, -0.115369, -0.279817, -0.276752, -0.294378, 0.0783223};

	*/
	
	// bin size of the offset angle used in the parametrization in degrees
	const double delta_angle = 10.;
	
	// LAr absorption length in cm
	const double L_abs = 2000.;

	// *************************************************************************************************
	//                    NUMBER OF VIS HITS PARAMETRIZATION
	// *************************************************************************************************
	// Note: wires are now incorporated into the parameterisations, no longer need to add scaling factor
	
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

	// DUNE-SP corrections
	// Argon, flat PDs (Arapucas/Supercells)
	const std::vector<double> vDistances_x_flat = {3, 13, 38, 63, 88, 113, 138, 163, 188, 213, 238, 263, 288, 313, 338, 353};			// cm	[16]
	const std::vector<double> vDistances_r_flat = {54.1115, 154.253, 269.706, 367.262, 496.739, 605.581, 693.216, 782.765, 862.208};		// cm	[9]

	const std::vector<std::vector<std::vector<double>>> fVISPars_flat = { 
		{ {1.67995, 1.5023, 1.20657, 1.02106, 0.882746, 0.776639, 0.695785, 0.629938, 0.582714, 0.545081, 0.515464, 0.493708, 0.472181, 0.45963, 0.439441, 0.419622},
		  {1.68963, 1.51982, 1.23126, 1.03816, 0.900409, 0.796292, 0.713189, 0.649865, 0.602912, 0.564819, 0.539057, 0.517199, 0.505427, 0.492089, 0.481653, 0.449229},
		  {1.73304, 1.55607, 1.26007, 1.07054, 0.934781, 0.832277, 0.759501, 0.704314, 0.658422, 0.627825, 0.606923, 0.591168, 0.584167, 0.580226, 0.563891, 0.529764},
		  {1.78193, 1.60802, 1.31866, 1.1289, 0.9946, 0.893599, 0.81836, 0.762132, 0.722593, 0.690561, 0.67056, 0.655973, 0.654003, 0.647708, 0.640756, 0.605073},
		  {1.96176, 1.77259, 1.47397, 1.2754, 1.13531, 1.02732, 0.942267, 0.880779, 0.828975, 0.794916, 0.771843, 0.751323, 0.74329, 0.744061, 0.733977, 0.70005},
		  {2.21549, 2.04521, 1.72647, 1.48965, 1.31053, 1.17082, 1.06559, 0.98202, 0.921138, 0.873365, 0.841385, 0.81795, 0.807915, 0.803652, 0.795182, 0.753378},
		  {2.3808, 2.20671, 1.89092, 1.64262, 1.4475, 1.29706, 1.17708, 1.08489, 1.01635, 0.962341, 0.926461, 0.904984, 0.894703, 0.894324, 0.890895, 0.844522},
		  {2.7039, 2.5078, 2.18003, 1.93111, 1.71334, 1.54146, 1.39972, 1.29637, 1.19956, 1.13696, 1.08587, 1.04728, 1.03061, 1.03128, 1.04224, 1.06682},
		  {3.28429, 3.14942, 2.81891, 2.46231, 2.13464, 1.86717, 1.63899, 1.46508, 1.32431, 1.22158, 1.13801, 1.07975, 1.04774, 1.0268, 1.00872, 0.934299}	  
		},
		{ {1.61717, 1.45536, 1.18141, 1.00703, 0.877535, 0.778688, 0.700656, 0.640282, 0.593421, 0.554809, 0.526049, 0.504949, 0.486008, 0.471742, 0.454623, 0.438439},
		  {1.63238, 1.47453, 1.21037, 1.02925, 0.896804, 0.79514, 0.717278, 0.658347, 0.611546, 0.575992, 0.548571, 0.52959, 0.516161, 0.505923, 0.495831, 0.470395},
		  {1.67016, 1.50472, 1.23162, 1.05456, 0.927059, 0.833322, 0.761788, 0.709941, 0.667458, 0.638414, 0.61695, 0.600871, 0.595867, 0.589934, 0.576962, 0.547778},
		  {1.72845, 1.56667, 1.29773, 1.1202, 0.992758, 0.898139, 0.826881, 0.772625, 0.73407, 0.70565, 0.685211, 0.673677, 0.669364, 0.665694, 0.661327, 0.62945},
		  {1.93974, 1.75853, 1.4716, 1.28379, 1.14714, 1.04343, 0.963162, 0.902172, 0.851864, 0.820813, 0.796031, 0.78107, 0.77271, 0.768238, 0.767237, 0.738134},
		  {2.11276, 1.94999, 1.65909, 1.44723, 1.28647, 1.16348, 1.0694, 0.995028, 0.937634, 0.896204, 0.86791, 0.850113, 0.843451, 0.841015, 0.839559, 0.805772},
		  {2.27089, 2.10995, 1.82097, 1.60253, 1.42969, 1.2978, 1.18599, 1.10159, 1.03826, 0.991185, 0.957985, 0.940432, 0.930111, 0.936145, 0.929783, 0.895972},
		  {2.675, 2.49434, 2.18241, 1.9411, 1.7415, 1.57744, 1.43963, 1.33051, 1.24657, 1.17447, 1.13066, 1.09834, 1.08741, 1.0859, 1.10706, 1.12759},
		  {2.96947, 2.85863, 2.57822, 2.2702, 1.98558, 1.75909, 1.56001, 1.40263, 1.27867, 1.18196, 1.11067, 1.05927, 1.03217, 1.0172, 0.9976, 0.934535}	  
		},
		{ {1.50825, 1.36889, 1.13392, 0.98191, 0.866922, 0.77688, 0.707548, 0.651111, 0.606337, 0.572121, 0.544768, 0.525048, 0.509095, 0.494043, 0.47836, 0.45858},
		  {1.52415, 1.38959, 1.16104, 1.00319, 0.886744, 0.794913, 0.724086, 0.670173, 0.625644, 0.593862, 0.570216, 0.551808, 0.53998, 0.532302, 0.522777, 0.495825},
		  {1.56948, 1.42471, 1.18729, 1.03115, 0.919158, 0.834466, 0.770977, 0.72419, 0.686219, 0.660066, 0.641763, 0.630002, 0.623834, 0.621193, 0.608531, 0.576298},
		  {1.64288, 1.50027, 1.26387, 1.10682, 0.992256, 0.906822, 0.84274, 0.793871, 0.758439, 0.733316, 0.71634, 0.706797, 0.70469, 0.703601, 0.697613, 0.662949},
		  {1.86287, 1.70058, 1.44341, 1.27538, 1.15264, 1.05872, 0.986552, 0.931587, 0.88889, 0.860201, 0.839109, 0.827028, 0.82099, 0.820279, 0.818245, 0.783072},
		  {1.96894, 1.83178, 1.5845, 1.40403, 1.26573, 1.15779, 1.07672, 1.01293, 0.961107, 0.927483, 0.903816, 0.891158, 0.88809, 0.89041, 0.887703, 0.850857},
		  {2.11486, 1.98106, 1.74011, 1.55249, 1.40211, 1.28869, 1.18925, 1.11334, 1.0575, 1.01646, 0.989677, 0.978543, 0.972021, 0.97899, 0.974276, 0.938631},
		  {2.4422, 2.2967, 2.05021, 1.85131, 1.67819, 1.53385, 1.41887, 1.3196, 1.24542, 1.18909, 1.14934, 1.12965, 1.12177, 1.13114, 1.14473, 1.17223},
		  {2.5988, 2.52471, 2.32597, 2.09048, 1.86465, 1.67534, 1.50487, 1.3751, 1.26983, 1.18503, 1.1275, 1.08634, 1.06218, 1.05617, 1.042, 0.977655}	  
		},
		{ {1.36203, 1.25143, 1.06416, 0.940405, 0.846393, 0.770026, 0.710644, 0.663008, 0.625575, 0.597437, 0.572537, 0.555241, 0.5426, 0.528512, 0.515126, 0.496204},
		  {1.37792, 1.27147, 1.08773, 0.961108, 0.864395, 0.789027, 0.729774, 0.681964, 0.647106, 0.619295, 0.599427, 0.585466, 0.575829, 0.57062, 0.562177, 0.533077},
		  {1.43685, 1.32193, 1.12809, 1.00146, 0.908067, 0.837362, 0.783724, 0.745501, 0.715154, 0.694964, 0.681247, 0.673038, 0.671752, 0.672051, 0.6613, 0.626777},
		  {1.56602, 1.4446, 1.24389, 1.10849, 1.00957, 0.935897, 0.880299, 0.839168, 0.809274, 0.788491, 0.775715, 0.770102, 0.771824, 0.774011, 0.770752, 0.734038},
		  {1.65725, 1.52999, 1.32736, 1.19666, 1.10039, 1.0278, 0.972371, 0.929854, 0.898858, 0.878205, 0.86608, 0.860904, 0.861646, 0.866037, 0.86557, 0.831242},
		  {1.78054, 1.67683, 1.48465, 1.34081, 1.22965, 1.14298, 1.07772, 1.0273, 0.983856, 0.959071, 0.943603, 0.93805, 0.9416, 0.947582, 0.949896, 0.911353},
		  {1.86241, 1.76882, 1.59574, 1.45292, 1.33642, 1.24684, 1.16578, 1.10885, 1.06611, 1.03589, 1.01871, 1.01663, 1.01979, 1.03478, 1.04109, 0.997536},
		  {2.02194, 1.92338, 1.76362, 1.63216, 1.522, 1.41935, 1.33986, 1.27154, 1.22612, 1.18624, 1.16734, 1.16245, 1.16453, 1.18303, 1.20928, 1.24244},
		  {2.1917, 2.15394, 2.03866, 1.88331, 1.71926, 1.57749, 1.44916, 1.34803, 1.26762, 1.20227, 1.16045, 1.13188, 1.12252, 1.12368, 1.11628, 1.04992}	  
		},
		{ {1.21192, 1.12501, 0.981162, 0.886385, 0.81395, 0.755101, 0.708174, 0.671753, 0.641806, 0.61953, 0.602772, 0.590021, 0.580067, 0.570257, 0.558218, 0.538113},
		  {1.23759, 1.15185, 1.01135, 0.912632, 0.836856, 0.778859, 0.73174, 0.695823, 0.668569, 0.649014, 0.634208, 0.625161, 0.621106, 0.617743, 0.612021, 0.582028},
		  {1.3372, 1.24048, 1.08385, 0.981411, 0.906339, 0.850626, 0.809228, 0.780127, 0.757369, 0.744707, 0.737064, 0.735092, 0.739165, 0.743367, 0.733825, 0.697567},
		  {1.39625, 1.3, 1.14429, 1.04039, 0.966035, 0.911477, 0.870675, 0.841821, 0.823572, 0.812433, 0.809316, 0.811214, 0.81902, 0.827248, 0.827171, 0.788705},
		  {1.44179, 1.34299, 1.19297, 1.09814, 1.0303, 0.980126, 0.943691, 0.917849, 0.901886, 0.893587, 0.891415, 0.897052, 0.906638, 0.916486, 0.92216, 0.88791},
		  {1.49443, 1.42303, 1.29061, 1.1906, 1.11335, 1.054, 1.01118, 0.977035, 0.952696, 0.941875, 0.937759, 0.942588, 0.955226, 0.97001, 0.976097, 0.937184},
		  {1.59488, 1.52823, 1.40851, 1.31078, 1.22903, 1.1642, 1.1103, 1.07635, 1.05194, 1.03894, 1.03624, 1.0454, 1.0604, 1.08511, 1.09897, 1.0563},
		  {1.6671, 1.60158, 1.50068, 1.41958, 1.35017, 1.28541, 1.23678, 1.2002, 1.1771, 1.16522, 1.15902, 1.17343, 1.19414, 1.22399, 1.26039, 1.28191},
		  {1.80691, 1.78974, 1.73667, 1.63844, 1.53272, 1.43474, 1.34692, 1.27759, 1.22442, 1.1885, 1.16391, 1.15584, 1.16103, 1.1786, 1.18017, 1.11401}	  
		},
		{ {1.15524, 1.07988, 0.958385, 0.879462, 0.818869, 0.77138, 0.734666, 0.705303, 0.6833, 0.667096, 0.655951, 0.649252, 0.643192, 0.63591, 0.623073, 0.601985},
		  {1.1885, 1.11803, 0.996096, 0.912411, 0.849796, 0.800743, 0.762943, 0.734867, 0.71462, 0.699669, 0.692586, 0.688702, 0.689522, 0.690546, 0.68734, 0.654688},
		  {1.21943, 1.14146, 1.01304, 0.930509, 0.872149, 0.830513, 0.800824, 0.780524, 0.769998, 0.764595, 0.76575, 0.771854, 0.779683, 0.788188, 0.783455, 0.74605},
		  {1.20788, 1.1374, 1.016, 0.939093, 0.884635, 0.847011, 0.820966, 0.805591, 0.798819, 0.798217, 0.804501, 0.816042, 0.829053, 0.845455, 0.849412, 0.807885},
		  {1.22609, 1.15651, 1.04232, 0.974567, 0.928931, 0.897108, 0.877528, 0.866212, 0.863744, 0.866935, 0.876967, 0.892472, 0.910572, 0.927446, 0.932572, 0.902381},
		  {1.2411, 1.19224, 1.09978, 1.03027, 0.978537, 0.939776, 0.916201, 0.895056, 0.892352, 0.893211, 0.902722, 0.91799, 0.938479, 0.962924, 0.974229, 0.933753},
		  {1.26582, 1.22115, 1.14423, 1.08429, 1.0345, 0.996373, 0.971323, 0.956652, 0.952058, 0.956109, 0.968782, 0.993272, 1.02211, 1.05529, 1.07656, 1.03333},
		  {1.30958, 1.26597, 1.20323, 1.15558, 1.11544, 1.08319, 1.06375, 1.04858, 1.04532, 1.05349, 1.06848, 1.0941, 1.12702, 1.17656, 1.21581, 1.2369},
		  {1.4037, 1.39983, 1.37779, 1.32183, 1.25879, 1.19901, 1.1496, 1.11269, 1.0877, 1.07653, 1.07314, 1.08548, 1.10827, 1.1378, 1.15578, 1.08895}	  
		},
		{ {1.27426, 1.19175, 1.05934, 0.973038, 0.909823, 0.861483, 0.821608, 0.792506, 0.771765, 0.757368, 0.747217, 0.742414, 0.737547, 0.731485, 0.71678, 0.692145},
		  {1.20255, 1.13089, 1.00916, 0.927784, 0.866592, 0.821308, 0.786168, 0.762307, 0.744728, 0.735681, 0.730353, 0.731412, 0.734882, 0.737942, 0.735619, 0.700012},
		  {1.11775, 1.04616, 0.93068, 0.857827, 0.808068, 0.773642, 0.75082, 0.738697, 0.73243, 0.732565, 0.738754, 0.748695, 0.762565, 0.776194, 0.769628, 0.732771},
		  {1.07179, 1.00731, 0.90292, 0.837423, 0.793408, 0.763976, 0.745886, 0.737359, 0.736678, 0.742324, 0.753317, 0.769762, 0.790495, 0.807173, 0.813678, 0.775096},
		  {1.03601, 0.974732, 0.88131, 0.827123, 0.790807, 0.768682, 0.757436, 0.752386, 0.755774, 0.765491, 0.780436, 0.799704, 0.82227, 0.842884, 0.853788, 0.823337},
		  {1.03307, 0.991846, 0.91548, 0.859529, 0.820516, 0.793129, 0.777652, 0.768362, 0.770054, 0.779399, 0.79407, 0.814441, 0.839751, 0.866729, 0.880961, 0.845997},
		  {1.03008, 0.992874, 0.929713, 0.880949, 0.843799, 0.817422, 0.802006, 0.796185, 0.79869, 0.808495, 0.827176, 0.856709, 0.884164, 0.920383, 0.943264, 0.906534},
		  {1.02415, 0.989935, 0.938482, 0.901237, 0.870795, 0.848773, 0.834976, 0.830254, 0.834168, 0.845807, 0.865375, 0.892356, 0.929334, 0.969198, 1.00822, 1.02599},
		  {1.03957, 1.03493, 1.01538, 0.973857, 0.928619, 0.889948, 0.858456, 0.837881, 0.827436, 0.826027, 0.834931, 0.852386, 0.879559, 0.910994, 0.936875, 0.913504}	  
		},
		{ {1.27426, 1.19175, 1.05934, 0.973038, 0.909823, 0.861483, 0.821608, 0.792506, 0.771765, 0.757368, 0.747217, 0.742414, 0.737547, 0.731485, 0.71678, 0.692145},
		  {1.20255, 1.13089, 1.00916, 0.927784, 0.866592, 0.821308, 0.786168, 0.762307, 0.744728, 0.735681, 0.730353, 0.731412, 0.734882, 0.737942, 0.735619, 0.700012},
		  {1.27204, 1.18739, 1.05069, 0.967695, 0.90854, 0.865532, 0.834186, 0.818909, 0.812163, 0.810419, 0.819609, 0.833564, 0.843202, 0.852653, 0.848203, 0.808108},
		  {1.16399, 1.0912, 0.971812, 0.897956, 0.847812, 0.814007, 0.793112, 0.782349, 0.77974, 0.785083, 0.797298, 0.813708, 0.836949, 0.850748, 0.858756, 0.823327},
		  {1.0813, 1.01036, 0.908115, 0.84808, 0.807807, 0.782758, 0.769239, 0.763294, 0.761293, 0.773114, 0.787471, 0.805788, 0.826063, 0.847037, 0.862569, 0.827133},
		  {1.02045, 0.976022, 0.894052, 0.833797, 0.791417, 0.761919, 0.742856, 0.734295, 0.73183, 0.737728, 0.749916, 0.76965, 0.79265, 0.817232, 0.830053, 0.807644},
		  {0.975507, 0.937157, 0.869956, 0.817599, 0.777362, 0.749731, 0.72941, 0.719385, 0.718641, 0.725071, 0.738938, 0.762324, 0.788287, 0.822578, 0.847537, 0.840179},
		  {0.966222, 0.927697, 0.869503, 0.827392, 0.792519, 0.765055, 0.749054, 0.737671, 0.736818, 0.741364, 0.756918, 0.782049, 0.816928, 0.862733, 0.911826, 0.947701},
		  {0.941417, 0.931757, 0.901292, 0.852774, 0.802177, 0.759511, 0.726692, 0.705822, 0.699634, 0.699843, 0.718528, 0.749023, 0.787584, 0.845292, 0.91141, 0.966078}	  
		},
		{ {1.27426, 1.19175, 1.05934, 0.973038, 0.909823, 0.861483, 0.821608, 0.792506, 0.771765, 0.757368, 0.747217, 0.742414, 0.737547, 0.731485, 0.71678, 0.692145},
		  {1.20255, 1.13089, 1.00916, 0.927784, 0.866592, 0.821308, 0.786168, 0.762307, 0.744728, 0.735681, 0.730353, 0.731412, 0.734882, 0.737942, 0.735619, 0.700012},
		  {1.27204, 1.18739, 1.05069, 0.967695, 0.90854, 0.865532, 0.834186, 0.818909, 0.812163, 0.810419, 0.819609, 0.833564, 0.843202, 0.852653, 0.848203, 0.808108},
		  {1.16399, 1.0912, 0.971812, 0.897956, 0.847812, 0.814007, 0.793112, 0.782349, 0.77974, 0.785083, 0.797298, 0.813708, 0.836949, 0.850748, 0.858756, 0.823327},
		  {1.0813, 1.01036, 0.908115, 0.84808, 0.807807, 0.782758, 0.769239, 0.763294, 0.761293, 0.773114, 0.787471, 0.805788, 0.826063, 0.847037, 0.862569, 0.827133},
		  {1.02045, 0.976022, 0.894052, 0.833797, 0.791417, 0.761919, 0.742856, 0.734295, 0.73183, 0.737728, 0.749916, 0.76965, 0.79265, 0.817232, 0.830053, 0.807644},
		  {0.975507, 0.937157, 0.869956, 0.817599, 0.777362, 0.749731, 0.72941, 0.719385, 0.718641, 0.725071, 0.738938, 0.762324, 0.788287, 0.822578, 0.847537, 0.840179},
		  {0.966222, 0.927697, 0.869503, 0.827392, 0.792519, 0.765055, 0.749054, 0.737671, 0.736818, 0.741364, 0.756918, 0.782049, 0.816928, 0.862733, 0.911826, 0.947701},
		  {0.941417, 0.931757, 0.901292, 0.852774, 0.802177, 0.759511, 0.726692, 0.705822, 0.699634, 0.699843, 0.718528, 0.749023, 0.787584, 0.845292, 0.91141, 0.966078}	  
		}
	};

	// Xenon, flat PDs (Arapucas/Supercells)
	/*
	const std::vector<double> vDistances_x_flat = {3, 13, 38, 63, 88, 113, 138, 163, 188, 213, 238, 263, 288, 313, 338, 353};			// cm	[16]
	const std::vector<double> vDistances_r_flat = {54.1115, 154.253, 269.706, 367.262, 496.739, 605.581, 693.216, 782.765, 862.208};		// cm	[9]

	const std::vector<std::vector<std::vector<double>>> fVISPars_flat = { 
		{ {1.42292, 1.31942, 1.14314, 1.03354, 0.954893, 0.895668, 0.850742, 0.813135, 0.784575, 0.769604, 0.758446, 0.747137, 0.748827, 0.763408, 0.776546, 0.769917},
		  {1.41836, 1.32156, 1.15245, 1.03549, 0.955398, 0.893329, 0.844977, 0.803317, 0.780003, 0.762924, 0.754932, 0.749316, 0.753423, 0.764773, 0.775107, 0.73847},
		  {1.43064, 1.32918, 1.15538, 1.04562, 0.957868, 0.902, 0.858058, 0.82298, 0.798712, 0.783286, 0.772047, 0.764398, 0.775728, 0.790624, 0.797297, 0.75866},
		  {1.46201, 1.36552, 1.20406, 1.09302, 1.01223, 0.951456, 0.905415, 0.868787, 0.844916, 0.826713, 0.817526, 0.824628, 0.829448, 0.83788, 0.832992, 0.824314},
		  {1.61987, 1.52623, 1.38658, 1.22113, 1.1464, 1.07976, 1.02709, 0.986773, 0.953123, 0.935845, 0.913943, 0.91304, 0.915525, 0.925434, 0.93213, 0.955269},
		  {1.89154, 1.78597, 1.64064, 1.48965, 1.38724, 1.29546, 1.21908, 1.16388, 1.11387, 1.08341, 1.06077, 1.04836, 1.05, 1.05101, 1.04239, 0.983244},
		  {2.05023, 1.96241, 1.79936, 1.66012, 1.5446, 1.45309, 1.36128, 1.30031, 1.24959, 1.21221, 1.18766, 1.16757, 1.16475, 1.17336, 1.15028, 1.11875},
		  {2.40085, 2.26858, 2.10522, 1.95488, 1.82863, 1.72549, 1.63033, 1.55232, 1.49184, 1.43598, 1.41357, 1.38856, 1.39596, 1.37445, 1.40755, 1.47731},
		  {3.00131, 2.96722, 2.87576, 2.69863, 2.50406, 2.32312, 2.15486, 2.01486, 1.8907, 1.79851, 1.72108, 1.65956, 1.61751, 1.59252, 1.5179, 1.49579}
		},
		{ {1.37407, 1.2769, 1.11892, 1.02321, 0.954931, 0.900042, 0.86009, 0.828382, 0.803265, 0.789315, 0.775509, 0.768375, 0.775618, 0.786733, 0.805713, 0.810876},
		  {1.36715, 1.28218, 1.12903, 1.02615, 0.954457, 0.896786, 0.854357, 0.817711, 0.797061, 0.784298, 0.773302, 0.773805, 0.780595, 0.796419, 0.810985, 0.788883},
		  {1.37248, 1.28618, 1.12822, 1.03165, 0.953403, 0.905577, 0.865954, 0.833164, 0.812087, 0.797174, 0.790532, 0.78481, 0.794702, 0.811046, 0.826019, 0.800091},
		  {1.41499, 1.33017, 1.18144, 1.08634, 1.01416, 0.959179, 0.918611, 0.886187, 0.862766, 0.847704, 0.843878, 0.850882, 0.85918, 0.868345, 0.872051, 0.872198},
		  {1.60827, 1.5209, 1.40127, 1.23055, 1.16318, 1.09904, 1.05255, 1.01438, 0.986219, 0.967248, 0.94471, 0.946333, 0.954531, 0.96617, 0.985107, 1.00816},
		  {1.79109, 1.69186, 1.56698, 1.43749, 1.34995, 1.27419, 1.20955, 1.16172, 1.12054, 1.09672, 1.08068, 1.07198, 1.07992, 1.07807, 1.0831, 1.04207},
		  {1.94105, 1.86262, 1.72586, 1.61135, 1.51152, 1.43401, 1.35988, 1.30749, 1.26334, 1.23276, 1.2112, 1.19643, 1.19967, 1.20709, 1.18601, 1.1724},
		  {2.38598, 2.2534, 2.10385, 1.96647, 1.85609, 1.75645, 1.67094, 1.59531, 1.53349, 1.4914, 1.46947, 1.45301, 1.4647, 1.44839, 1.4737, 1.55239},
		  {2.72183, 2.68984, 2.63433, 2.50323, 2.34891, 2.2024, 2.06186, 1.93903, 1.8393, 1.75606, 1.69391, 1.6483, 1.60708, 1.58162, 1.50839, 1.48248}
		},
		{ {1.2794, 1.20006, 1.07694, 1.00277, 0.950365, 0.908544, 0.87576, 0.851775, 0.83527, 0.825862, 0.818268, 0.815444, 0.823338, 0.835547, 0.858962, 0.859531},
		  {1.27703, 1.20877, 1.08646, 1.00402, 0.948941, 0.906946, 0.870997, 0.842256, 0.829672, 0.820007, 0.815418, 0.818627, 0.828569, 0.844322, 0.862867, 0.833061},
		  {1.2874, 1.21515, 1.09059, 1.01232, 0.950562, 0.915758, 0.885219, 0.862458, 0.845484, 0.836697, 0.833105, 0.831037, 0.846697, 0.862853, 0.884124, 0.848749},
		  {1.34102, 1.27211, 1.15161, 1.07866, 1.02081, 0.97755, 0.94591, 0.920072, 0.903687, 0.893782, 0.892264, 0.905507, 0.916919, 0.927907, 0.932861, 0.928055},
		  {1.53413, 1.45896, 1.35954, 1.22871, 1.16754, 1.12509, 1.08659, 1.05639, 1.03529, 1.02295, 1.00668, 1.0114, 1.02304, 1.03636, 1.05397, 1.08081},
		  {1.65867, 1.59076, 1.49496, 1.39836, 1.3354, 1.27703, 1.22807, 1.19247, 1.16281, 1.14552, 1.13656, 1.13351, 1.14515, 1.15131, 1.15634, 1.11231},
		  {1.80453, 1.74518, 1.65345, 1.57666, 1.5019, 1.44094, 1.38424, 1.34313, 1.30982, 1.28811, 1.27395, 1.2624, 1.2718, 1.28459, 1.26347, 1.24585},
		  {2.10153, 2.05328, 1.97521, 1.90129, 1.82371, 1.75283, 1.68732, 1.63021, 1.58068, 1.54824, 1.52926, 1.5136, 1.50177, 1.52156, 1.55705, 1.62856},
		  {2.36231, 2.37399, 2.39017, 2.3408, 2.23382, 2.14055, 2.03151, 1.94209, 1.86046, 1.79827, 1.74587, 1.70999, 1.67968, 1.66205, 1.59837, 1.57015}
		},
		{ {1.1587, 1.10258, 1.01763, 0.972936, 0.942476, 0.918824, 0.90139, 0.890344, 0.883709, 0.88271, 0.886015, 0.890102, 0.902704, 0.920961, 0.950356, 0.959914},
		  {1.15581, 1.11106, 1.0278, 0.974191, 0.941421, 0.916465, 0.897497, 0.88092, 0.878494, 0.878512, 0.88354, 0.891103, 0.908239, 0.932995, 0.959165, 0.932838},
		  {1.18233, 1.13225, 1.04602, 0.996879, 0.954006, 0.937118, 0.921234, 0.909001, 0.903469, 0.905103, 0.90932, 0.91506, 0.936128, 0.962579, 0.989865, 0.955045},
		  {1.27987, 1.22599, 1.13806, 1.09147, 1.05338, 1.02533, 1.00669, 0.991285, 0.984225, 0.983176, 0.989949, 1.01319, 1.02911, 1.0466, 1.05565, 1.05465},
		  {1.37098, 1.31428, 1.25887, 1.16969, 1.13438, 1.11293, 1.09334, 1.082, 1.07498, 1.07414, 1.06876, 1.08114, 1.10154, 1.12102, 1.14582, 1.17267},
		  {1.50236, 1.46554, 1.41301, 1.36147, 1.32546, 1.29033, 1.26487, 1.24705, 1.23106, 1.2259, 1.22713, 1.23504, 1.25331, 1.26746, 1.28144, 1.23209},
		  {1.60315, 1.5783, 1.53919, 1.50024, 1.47074, 1.44497, 1.40436, 1.389, 1.372, 1.3651, 1.36882, 1.37046, 1.38737, 1.41173, 1.41223, 1.37083},
		  {1.75648, 1.73054, 1.71033, 1.69857, 1.67488, 1.65189, 1.62651, 1.60748, 1.59238, 1.58493, 1.57967, 1.59267, 1.61092, 1.62596, 1.65529, 1.74053},
		  {2.00571, 2.04056, 2.13422, 2.15003, 2.11537, 2.07883, 2.01903, 1.96674, 1.92093, 1.88048, 1.85433, 1.83814, 1.83025, 1.81992, 1.76921, 1.72067}
		},
		{ {1.04658, 1.00778, 0.95966, 0.942102, 0.936297, 0.934151, 0.935516, 0.940309, 0.949998, 0.962494, 0.977172, 0.990667, 1.01269, 1.03264, 1.07541, 1.09656},
		  {1.05348, 1.02322, 0.974812, 0.949217, 0.941289, 0.936595, 0.935076, 0.933571, 0.948089, 0.959977, 0.977878, 0.996721, 1.02375, 1.05695, 1.09578, 1.07113},
		  {1.11199, 1.07678, 1.02098, 0.99726, 0.977216, 0.981083, 0.982685, 0.987138, 0.995204, 1.00864, 1.02642, 1.0376, 1.07245, 1.10808, 1.1523, 1.11488},
		  {1.1551, 1.11778, 1.07282, 1.04944, 1.03799, 1.02988, 1.02786, 1.03494, 1.04407, 1.05898, 1.07942, 1.09781, 1.13716, 1.16616, 1.19515, 1.18285},
		  {1.20648, 1.16921, 1.15214, 1.10204, 1.08995, 1.09863, 1.10272, 1.11084, 1.12357, 1.14067, 1.15096, 1.17381, 1.20943, 1.23644, 1.2713, 1.30183},
		  {1.28456, 1.26706, 1.26798, 1.25408, 1.25225, 1.25116, 1.25023, 1.25701, 1.26068, 1.27704, 1.2955, 1.31538, 1.34945, 1.37064, 1.39389, 1.34066},
		  {1.38944, 1.38761, 1.39433, 1.40126, 1.40151, 1.40766, 1.40087, 1.40847, 1.41943, 1.43345, 1.45508, 1.48215, 1.51808, 1.55559, 1.57311, 1.5031},
		  {1.46756, 1.46635, 1.4887, 1.5255, 1.54733, 1.56291, 1.5749, 1.59327, 1.60871, 1.62682, 1.65092, 1.68253, 1.72136, 1.75329, 1.78396, 1.8622},
		  {1.66989, 1.72621, 1.86738, 1.94471, 1.96866, 1.98785, 1.97897, 1.97116, 1.96262, 1.96178, 1.96451, 1.97127, 1.98368, 2.00068, 1.9676, 1.88032}
		},
		{ {1.03118, 1.00281, 0.976428, 0.980169, 0.991869, 1.00883, 1.02906, 1.05058, 1.0795, 1.11169, 1.14202, 1.17375, 1.20765, 1.23382, 1.28736, 1.3177},
		  {1.04496, 1.02017, 0.997629, 0.995814, 1.00113, 1.01321, 1.03076, 1.05255, 1.0761, 1.1072, 1.13805, 1.17251, 1.21393, 1.25785, 1.31348, 1.29418},
		  {1.04761, 1.02269, 0.996304, 0.995023, 1.0007, 1.01443, 1.03109, 1.05143, 1.07711, 1.10826, 1.1405, 1.17231, 1.21302, 1.25633, 1.3085, 1.27693},
		  {1.03766, 1.01826, 0.994187, 0.991873, 0.998316, 1.0141, 1.03499, 1.05694, 1.0847, 1.11766, 1.15486, 1.20303, 1.24369, 1.28379, 1.31056, 1.30615},
		  {1.05731, 1.03858, 1.01847, 1.02548, 1.05105, 1.06778, 1.09361, 1.12376, 1.15584, 1.19305, 1.22445, 1.2704, 1.30898, 1.34238, 1.38562, 1.409},
		  {1.10877, 1.10754, 1.14076, 1.15307, 1.18504, 1.20609, 1.22933, 1.25754, 1.28835, 1.32514, 1.36309, 1.39769, 1.45247, 1.49321, 1.52081, 1.46434},
		  {1.1436, 1.15103, 1.18693, 1.22099, 1.25743, 1.29632, 1.31765, 1.35473, 1.39353, 1.4329, 1.48052, 1.53162, 1.58169, 1.64218, 1.67934, 1.59459},
		  {1.19421, 1.20256, 1.25152, 1.317, 1.36939, 1.41774, 1.46407, 1.51032, 1.55923, 1.60923, 1.66674, 1.71655, 1.77754, 1.83076, 1.8659, 1.93086},
		  {1.34172, 1.40313, 1.56572, 1.68537, 1.75519, 1.81807, 1.85753, 1.89438, 1.93073, 1.96484, 2.00243, 2.04032, 2.08798, 2.12176, 2.12229, 2.00869}
		},
		{ {1.16989, 1.14166, 1.11952, 1.12862, 1.15072, 1.1794, 1.21092, 1.24722, 1.28841, 1.33695, 1.38336, 1.42854, 1.47133, 1.49804, 1.54741, 1.60228},
		  {1.10141, 1.07998, 1.06145, 1.06426, 1.08325, 1.10846, 1.13749, 1.16747, 1.20553, 1.25047, 1.29551, 1.34354, 1.40275, 1.45218, 1.51533, 1.50297},
		  {1.01096, 0.991787, 0.971516, 0.976379, 0.991447, 1.0196, 1.05113, 1.08405, 1.12178, 1.1661, 1.21459, 1.25266, 1.3159, 1.36625, 1.42783, 1.38305},
		  {0.969217, 0.952436, 0.942125, 0.952499, 0.972033, 0.996329, 1.02925, 1.06522, 1.10686, 1.1524, 1.20293, 1.26395, 1.32181, 1.3664, 1.3952, 1.38077},
		  {0.947213, 0.934058, 0.939794, 0.941679, 0.973627, 1.00591, 1.04379, 1.08653, 1.13184, 1.1822, 1.22609, 1.28252, 1.33969, 1.38291, 1.42353, 1.43162},
		  {0.974022, 0.978218, 1.02104, 1.04427, 1.08599, 1.11956, 1.15665, 1.20023, 1.24537, 1.2971, 1.34998, 1.39711, 1.46284, 1.52256, 1.55725, 1.50591},
		  {0.983942, 0.997837, 1.03974, 1.0837, 1.12872, 1.17995, 1.21721, 1.26786, 1.3201, 1.37515, 1.43636, 1.50599, 1.57261, 1.639, 1.69817, 1.60824},
		  {0.986637, 1.00021, 1.05809, 1.12605, 1.18612, 1.24586, 1.30626, 1.36646, 1.42805, 1.49162, 1.55908, 1.62912, 1.70098, 1.77127, 1.82632, 1.83638},
		  {1.05328, 1.11283, 1.26775, 1.3846, 1.46585, 1.54778, 1.6025, 1.66216, 1.7191, 1.77506, 1.83249, 1.89588, 1.95564, 2.01094, 2.03729, 1.92352}
		},
		{ {1.16989, 1.14166, 1.11952, 1.12862, 1.15072, 1.1794, 1.21092, 1.24722, 1.28841, 1.33695, 1.38336, 1.42854, 1.47133, 1.49804, 1.54741, 1.60228},
		  {1.10141, 1.07998, 1.06145, 1.06426, 1.08325, 1.10846, 1.13749, 1.16747, 1.20553, 1.25047, 1.29551, 1.34354, 1.40275, 1.45218, 1.51533, 1.50297},
		  {1.18138, 1.14739, 1.12995, 1.13489, 1.15539, 1.18368, 1.21695, 1.25575, 1.29879, 1.35524, 1.41206, 1.46715, 1.5353, 1.58196, 1.64623, 1.61357},
		  {1.07715, 1.05718, 1.04135, 1.0538, 1.07493, 1.10168, 1.13983, 1.17895, 1.22699, 1.2772, 1.33161, 1.41493, 1.47309, 1.51852, 1.54152, 1.53889},
		  {1.02121, 1.00515, 1.03757, 1.0002, 1.03372, 1.07104, 1.11477, 1.15716, 1.20884, 1.26404, 1.30206, 1.36415, 1.43529, 1.4845, 1.53161, 1.54125},
		  {0.991434, 0.993547, 1.04104, 1.06174, 1.09845, 1.13584, 1.17396, 1.22239, 1.26504, 1.32346, 1.38022, 1.4284, 1.49953, 1.55411, 1.60144, 1.5439},
		  {0.958746, 0.973306, 1.01366, 1.05779, 1.09917, 1.14531, 1.18026, 1.23175, 1.28569, 1.33957, 1.40024, 1.45703, 1.52317, 1.60065, 1.64575, 1.56922},
		  {0.958867, 0.972587, 1.02301, 1.08662, 1.14244, 1.19739, 1.25631, 1.31383, 1.37222, 1.43563, 1.50117, 1.56715, 1.63567, 1.70922, 1.76777, 1.7753},
		  {0.986332, 1.04235, 1.18522, 1.29253, 1.36492, 1.43449, 1.4859, 1.53936, 1.59273, 1.64912, 1.70134, 1.7577, 1.81235, 1.86855, 1.89174, 1.80892}
		},
		{ {1.16989, 1.14166, 1.11952, 1.12862, 1.15072, 1.1794, 1.21092, 1.24722, 1.28841, 1.33695, 1.38336, 1.42854, 1.47133, 1.49804, 1.54741, 1.60228},
		  {1.10141, 1.07998, 1.06145, 1.06426, 1.08325, 1.10846, 1.13749, 1.16747, 1.20553, 1.25047, 1.29551, 1.34354, 1.40275, 1.45218, 1.51533, 1.50297},
		  {1.18138, 1.14739, 1.12995, 1.13489, 1.15539, 1.18368, 1.21695, 1.25575, 1.29879, 1.35524, 1.41206, 1.46715, 1.5353, 1.58196, 1.64623, 1.61357},
		  {1.07715, 1.05718, 1.04135, 1.0538, 1.07493, 1.10168, 1.13983, 1.17895, 1.22699, 1.2772, 1.33161, 1.41493, 1.47309, 1.51852, 1.54152, 1.53889},
		  {1.02121, 1.00515, 1.03757, 1.0002, 1.03372, 1.07104, 1.11477, 1.15716, 1.20884, 1.26404, 1.30206, 1.36415, 1.43529, 1.4845, 1.53161, 1.54125},
		  {0.991434, 0.993547, 1.04104, 1.06174, 1.09845, 1.13584, 1.17396, 1.22239, 1.26504, 1.32346, 1.38022, 1.4284, 1.49953, 1.55411, 1.60144, 1.5439},
		  {0.958746, 0.973306, 1.01366, 1.05779, 1.09917, 1.14531, 1.18026, 1.23175, 1.28569, 1.33957, 1.40024, 1.45703, 1.52317, 1.60065, 1.64575, 1.56922},
		  {0.958867, 0.972587, 1.02301, 1.08662, 1.14244, 1.19739, 1.25631, 1.31383, 1.37222, 1.43563, 1.50117, 1.56715, 1.63567, 1.70922, 1.76777, 1.7753},
		  {0.986332, 1.04235, 1.18522, 1.29253, 1.36492, 1.43449, 1.4859, 1.53936, 1.59273, 1.64912, 1.70134, 1.7577, 1.81235, 1.86855, 1.89174, 1.80892}
		}
	};
*/		
	

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

	// solid angle of dome (PMTs)
	double Omega_Dome_Model(const double distance, const double theta) const;

	// linear interpolation function
	double interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate );

};

#endif
