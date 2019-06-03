#include "semi_analytic_hits.h"

// implementation of semi-analytic model for number of incident photons

#include <iostream>
#include <cmath>

#include "TRandom.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFormula.h"
#include "Math/SpecFuncMathMore.h"

using namespace std;

// constructor
semi_analytic_hits::semi_analytic_hits() {

  // load mathmore library
  gSystem->Load("libMathMore.so");
  if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
  _mathmore_loaded_ = true;

  std::cout<<"Light simulation for DUNE Single Phase detector"<<std::endl;

  // initialise gaisser hillas functions for VUV Rayleigh scattering correction
  std::cout<<"Loading Gaisser Hillas corrections for VUV semi-analytic model."<< std::endl;
  double pars_ini[4] = {0,0,0,0};
  for(int bin = 0; bin < 9; bin++) {
    GH[bin] =  new TF1("GH",GaisserHillas,0.,2000,4);      
    for(int j=0; j < 4; j++) {
      pars_ini[j] = GH_RS60cm_SP[j][bin];		
    }
    GH[bin]->SetParameters(pars_ini);
  }
  

  // initialise pol5 functions for visible hits correction
  std::cout<<"Loading corrections for visible semi-analytic model."<< std::endl;
  double pars_ini_vis[6] = {0,0,0,0,0,0};
  for (int bin = 0; bin < 9; bin++) {
    VIS_pol[bin] = new TF1 ("pol", "pol5", 0, 2000);
    for (int j = 0; j < 6; j++){
      pars_ini_vis[j] = VIS_RS60cm_SP[j][bin];
    }
    VIS_pol[bin]->SetParameters(pars_ini_vis);
  }

  // initialise pol7 functions for visible hits correction, crowns model
  double pars_ini_vis_crowns[8] = {0,0,0,0,0,0,0,0};
  for (int bin = 0; bin < 9; bin++) {
    VIS_pol_crowns[bin] = new TF1 ("pol", "pol7", 0, 2000);
    for (int j = 0; j < 8; j++){
      pars_ini_vis_crowns[j] = VIS_RS60cm_SP_Crowns[j][bin];
    }
    VIS_pol_crowns[bin]->SetParameters(pars_ini_vis_crowns);
  }

  std::cout << std::endl;

}

// VUV hits calculation
int semi_analytic_hits::VUVHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type) {
  gRandom->SetSeed(0);

  // distance and angle between ScintPoint and OpDetPoint
  double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
  double cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance;
  double theta = acos(cosine)*180./pi;

  // calculate solid angle:
  double solid_angle = 0;
  // rectangular aperture
  if (optical_detector_type == 1) {
    // set Arapuca geometry struct for solid angle function
    acc detPoint; 
    detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
    detPoint.w = y_dimension_detector; detPoint.h = z_dimension_detector; // width and height in cm of arapuca active window

    // get scintillation point coordinates relative to arapuca window centre
    TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;  

    // calculate solid angle
    solid_angle = solid(detPoint, ScintPoint_rel);
  }
  // disk aperture
  else if (optical_detector_type == 0) {
    // offset in z-y plane
    double d = sqrt(pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
    // drift distance (in x)
    double h =  sqrt(pow(ScintPoint[0] - OpDetPoint[0],2));
    // Solid angle of a disk
    solid_angle = Disk_SolidAngle(d, h, radius);
  }
  else {
    std::cout << "Error: Invalid optical detector type." << endl;
    exit(1);
  }  

  // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
  double hits_geo = exp(-1.*distance/L_abs) * (solid_angle / (4*pi)) * Nphotons_created;

  // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
  // offset angle bin
  int j = (theta/delta_angulo);
  double hits_rec = gRandom->Poisson( GH[j]->Eval(distance)*hits_geo/cosine );

  // round to integer value, cannot have non-integer number of hits
  int hits_vuv = std::round(hits_rec);

  return hits_vuv;
}

// Visible hits calculation
int semi_analytic_hits::VisHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type) {
  gRandom->SetSeed(0);

  // 1). calculate total number of hits of VUV photons on reflective foils via solid angle + Gaisser-Hillas corrections:

  // set cathode plane struct for solid angle function
  acc cathode_plane; 
  cathode_plane.ax = x_foils; cathode_plane.ay = y_foils; cathode_plane.az = z_foils;   // centre coordinates of cathode plane
  cathode_plane.w = y_dimension_foils; cathode_plane.h = z_dimension_foils;                        // width and height in cm

  // get scintpoint coords relative to centre of cathode plane
  TVector3 cathodeCentrePoint(x_foils,y_foils,z_foils);
  TVector3 ScintPoint_relative = ScintPoint - cathodeCentrePoint; 

  // calculate solid angle of cathode from the scintillation point
  double solid_angle_cathode = solid(cathode_plane, ScintPoint_relative);

  // calculate distance and angle between ScintPoint and hotspot
  // vast majority of hits in hotspot region directly infront of scintpoint,therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
  double distance_cathode = std::abs(x_foils - ScintPoint[0]); 

  double cosine_cathode = 1;
  double theta_cathode = 0;

  // calculate hits on cathode plane via geometric acceptance
  double cathode_hits_geo = exp(-1.*distance_cathode/L_abs) * (solid_angle_cathode / (4.*pi)) * Nphotons_created;

  // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
  // offset angle bin
  int j = (theta_cathode/delta_angulo);
  double cathode_hits_rec = GH[j]->Eval(distance_cathode)*cathode_hits_geo/cosine_cathode;
  
  // 2). calculate number of these hits which reach the optical channel from the hotspot via solid angle 
  
  // calculate hotspot location  
  TVector3 v_to_wall(x_foils - ScintPoint[0],0,0);        
  TVector3 hotspot = ScintPoint + v_to_wall;

  // solid angle :
  double solid_angle_detector = 0;
  // rectangular aperture
  if (optical_detector_type == 1) {
    // set Arapuca geometry struct for solid angle function
    acc detPoint; 
    detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
    detPoint.w = y_dimension_detector; detPoint.h = z_dimension_detector; // width and height in cm of arapuca active window

    
    // get hotspot coordinates relative to detpoint
    TVector3 emission_relative = hotspot - OpDetPoint;

    // calculate solid angle of optical channel
    solid_angle_detector = solid(detPoint, emission_relative);
  }  
  // disk aperture
  else if (optical_detector_type == 0) {
    // offset in z-y plane
    double d = sqrt(pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
    // drift distance (in x)
    double h =  sqrt(pow(hotspot[0] - OpDetPoint[0],2));
    // Solid angle of a disk
    solid_angle_detector = Disk_SolidAngle(d, h, radius);
  }
  else {
    std::cout << "Error: Invalid optical detector type." << endl;
    exit(1);
  }

  // calculate number of hits via geometeric acceptance  
  double hits_geo = (solid_angle_detector / (2*pi)) * cathode_hits_rec;

  // distance to hotspot
  double distance_vuv = sqrt(pow(ScintPoint[0] - hotspot[0],2) + pow(ScintPoint[1] - hotspot[1],2) + pow(ScintPoint[2] - hotspot[2],2));
  // distance from hotspot to arapuca
  double distance_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2) + pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
  // angle between hotspot and arapuca
  double cosine_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2)) / distance_vis;
  double theta_vis = acos(cosine_vis)*180./pi;

  // apply correction curves, 5th order polynomial 
  int k = (theta_vis/delta_angle);
  double hits_rec = VIS_pol[k]->Eval(distance_vuv)*hits_geo/cosine_vis;

  // apply border correction
  // interpolate in x for each y bin
  std::vector<double> interp_vals = {0,0,0,0,0,0};
  for (int i = 0; i < 6; i++){
    interp_vals[i] = interpolate(vDistances_x, VIS_RS60_SP_Borders[k][i], std::abs(ScintPoint[0]), false);
  }
  // interpolate in y
  double border_correction = interpolate(vDistances_y, interp_vals, std::abs(ScintPoint[1]), false);
  // apply correction
  double hits_rec_borders = border_correction * hits_rec / cosine_vis;

  // round final result
  int hits_vis = std::round(hits_rec_borders);

  return hits_vis;
}

// Visible hits calculation using updated crowns method
int semi_analytic_hits::VisHits_crowns(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type) {

  gRandom->SetSeed(0);

  // utilising concentric square rings corresponding to approximately 0-10, 10-20, etc. degrees on cathode plane to model number of hits geometrically
  // GH corrections applied to account for effects of Rayleigh scattering for each concentric ring
  // hits on optical detector calculated via solid angle from four points corresponding to middle top, bottom, left, right of ring
  // corrections applied to account for reflections/absorption from borders 

  // 1). calculate total number of hits of VUV photons on reflective foils via solid angle + Gaisser-Hillas corrections:
  
  // determine boundaries of cathode plane
  double top_boundary = y_foils + y_dimension_foils/2;
  double bottom_boundary = y_foils - y_dimension_foils/2;
  double left_boundary = z_foils + z_dimension_foils/2;
  double right_boundary = z_foils - z_dimension_foils/2;


  // determine square dimensions corresponding to each angular bin, calculating solid angle of each square
  // containers to store dimensions and positions of each square
  // full size of square
  double size_cathode_squares[9];
  // centre and coordinates one boundary walls taken into account
  TVector3 square_centres[9];
  double top_cathode_squares[9];
  double bottom_cathode_squares[9];
  double left_cathode_squares[9];
  double right_cathode_squares[9];  

  // solid angle of each square  
  double solid_angle_cathode_squares[9];
  
  // distance to cathode from scintillation point
  double distance_cathode = std::abs(plane_depth - ScintPoint[0]);
  
  // loop over angle bins determining dimensions, coordinates and solid angle of each square
  for (int bin = 0; bin < 9; bin++) {

    // angle
    double theta_boundary;
    if (bin == 8) theta_boundary = (bin+1)*10 - 0.01;  // special case for 80-90 degree bin, take maximum as 89.99 (90 impossible)
    else theta_boundary = (bin+1) * 10;
    
    // size of square
    size_cathode_squares[bin] = distance_cathode * tan(theta_boundary * (pi/180));
    
    // dimensions of squares checking against boundaries:
    // top
    if (ScintPoint[1] + size_cathode_squares[bin] < top_boundary) top_cathode_squares[bin] = ScintPoint[1] + size_cathode_squares[bin];
    else top_cathode_squares[bin] = top_boundary;
    // bottom
    if (ScintPoint[1] - size_cathode_squares[bin] > bottom_boundary) bottom_cathode_squares[bin] = ScintPoint[1] - size_cathode_squares[bin];
    else bottom_cathode_squares[bin] = bottom_boundary;
    // left
    if (ScintPoint[2] + size_cathode_squares[bin] < left_boundary) left_cathode_squares[bin] = ScintPoint[2] + size_cathode_squares[bin];
    else left_cathode_squares[bin] = left_boundary;
    // right
    if (ScintPoint[2] - size_cathode_squares[bin] > right_boundary) right_cathode_squares[bin] = ScintPoint[2] - size_cathode_squares[bin];
    else right_cathode_squares[bin] = right_boundary;

    // set square geometry struct for solid angle function
    acc square; 
    square.ax = plane_depth; square.ay = (top_cathode_squares[bin] + bottom_cathode_squares[bin])/2; square.az = (right_cathode_squares[bin] + left_cathode_squares[bin])/2;  // centre coordinates of square
    square.h = std::abs(left_cathode_squares[bin] - right_cathode_squares[bin]); square.w = std::abs(top_cathode_squares[bin] - bottom_cathode_squares[bin]); // width and height of square

    // get scintillation point coordinates relative to centre of square
    TVector3 square_centre(square.ax, square.ay, square.az); square_centres[bin] = square_centre;
    TVector3 ScintPoint_rel = square_centre - ScintPoint;  

    // calculate solid angle of square
    solid_angle_cathode_squares[bin] = solid(square, ScintPoint_rel);

    // debugging
    //cout << "Square solid angle [" << bin*10 << ", " << (bin+1)*10 << "]: " << solid_angle_cathode_squares[bin] << endl;
  }


  // use solid angle squares to determine solid angle of angular bin ring
  // calculate number of hits on each ring applying appropriate GH correction

  // hits on each ring
  double cathode_hits_rec[9];

  // find maximum distance to boundary wall from scintpoint in y,z directions 
  double distances_to_boundary[] = {abs(top_boundary - ScintPoint[1]), abs(bottom_boundary - ScintPoint[1]), abs(left_boundary - ScintPoint[2]), abs(right_boundary - ScintPoint[2])};
  double max_distance_to_boundary = *(std::max_element(distances_to_boundary, distances_to_boundary+4));
  
  // loop over angular bins
  for (int bin = 0; bin < 9; bin ++) {
    
    // solid angle of ring for this bin + distance for GH correction
    double bin_solid_angle;
    double bin_distance;
    // centre special case, use full square solid angle & distance to hotspot directly infront of scintillation point
    if (bin == 0){
      bin_solid_angle = solid_angle_cathode_squares[bin];
      bin_distance = distance_cathode;
    }
    else {
      // solid angle is difference between this square and previous square
      bin_solid_angle = solid_angle_cathode_squares[bin] - solid_angle_cathode_squares[bin - 1];
      // distance to the middle of this bin, checking against boundary of cathode     
      double distance_square;
      if (size_cathode_squares[bin] > max_distance_to_boundary) distance_square = size_cathode_squares[bin - 1] + ((max_distance_to_boundary - size_cathode_squares[bin - 1]) / 2);
      //if (size_cathode_squares[bin] > max_distance_to_boundary) distance_square = size_cathode_squares[bin - 1];
      else distance_square = size_cathode_squares[bin - 1] + ((size_cathode_squares[bin] - size_cathode_squares[bin - 1]) / 2);
      bin_distance = sqrt(pow(distance_cathode, 2) + pow(distance_square, 2));
    }

    // debugging
    //cout << "Bin: " << bin << "   bin solid angle = " << bin_solid_angle << "   bin distance = " << bin_distance << "   size square: " << size_cathode_squares[bin]; 

    // calculate hits on this ring via geometric acceptance
    double cathode_hits_geo = exp(-1.*bin_distance/L_abs) * (bin_solid_angle / (4.*pi)) * Nphotons_created;
    
    // apply Gaisser-Hillas correction for Rayleigh scattering
    // offset angle bin
    double theta_boundary = 10 + bin*10;
    int j = (theta_boundary-5) / delta_angulo;  // minus 5 degrees - middle of bin
    double cosine_cathode = cos((theta_boundary - 5) * (pi/180));
    // corrected hits for this bin
    cathode_hits_rec[bin] = (GH[j]->Eval(bin_distance))*cathode_hits_geo/cosine_cathode;  

    // debugging
    //cout << "   cathode hits = " << cathode_hits_rec[bin] << endl;
  }


  // 2). calculate number of hit from each ring which reach the optical detector via solid angle 

  // set Arapuca geometry struct for solid angle function
  acc detPoint; 
  detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
  detPoint.w = y_dimension_detector; detPoint.h = z_dimension_detector;                   // width and height in cm of arapuca active window

  // calculate hits from each ring via geometric acceptance 
  double hits_geo[9];
  // loop over rings
  for (int bin = 0; bin < 9; bin++) {
        
    // first bin, hits from centre of first square
    if (bin == 0) {
      double solid_angle_detector = 0;
      // rectangle
      if (optical_detector_type == 1) {
        TVector3 emission_relative = square_centres[0] - OpDetPoint;
        solid_angle_detector = solid(detPoint, emission_relative);
      }
      // disk
      else if (optical_detector_type == 0) {
        // offset in z-y plane
        double d = sqrt(pow(square_centres[0][1] - OpDetPoint[1],2) + pow(square_centres[0][2] - OpDetPoint[2],2));
        // drift distance (in x)
        double h =  sqrt(pow(square_centres[0][0] - OpDetPoint[0],2));
        // Solid angle of a disk
        solid_angle_detector = Disk_SolidAngle(d, h, radius);
      }
      else {
        std::cout << "Erorr: Invalid optical detector type." << std::endl;
        exit(1);
      }
      double hits_geo_ring = (solid_angle_detector / (2*pi)) * cathode_hits_rec[bin];
      // store
      hits_geo[bin] = hits_geo_ring;
    }
    // for remaining bins, take four points around ring weighted by width of ring in each direction to account for rings cut-off by cathode boundaries
    else {      
      
      // widths of ring in each direction acting as weighting factor for solid angle hits from each spot
      double widths[4];
      widths[0] = abs(top_cathode_squares[bin] - top_cathode_squares[bin - 1]);       // top width
      widths[1] = abs(bottom_cathode_squares[bin] - bottom_cathode_squares[bin - 1]); // bottom width
      widths[2] = abs(left_cathode_squares[bin] - left_cathode_squares[bin - 1]);     // left width
      widths[3] = abs(right_cathode_squares[bin] - right_cathode_squares[bin - 1]);   // right width
      // sum for normalisation of weighting factor
      double width_sum = widths[0] + widths[1] + widths[2] + widths[3];

      // coordinates of points on ring
      TVector3 spot[4];
      spot[0] = square_centres[bin]; spot[0][1] = top_cathode_squares[bin] - (widths[0] / 2);
      spot[1] = square_centres[bin]; spot[1][1] = bottom_cathode_squares[bin] + (widths[1] / 2);
      spot[2] = square_centres[bin]; spot[2][2] = left_cathode_squares[bin] - (widths[2] / 2);
      spot[3] = square_centres[bin]; spot[3][2] = right_cathode_squares[bin] + (widths[3] / 2);

      // calculate hits from each spot via geometric acceptance adding to total from ring
      hits_geo[bin] = 0;
      // loop over the fours spots
      for (int i = 0; i < 4; i++) {
        // solid angle from spot
        double solid_angle_detector = 0;
        // recantular
        if (optical_detector_type == 1) {
          TVector3 emission_relative = spot[i] - OpDetPoint;
          double solid_angle_detector = solid(detPoint, emission_relative);
        }
        // disk
        else if (optical_detector_type == 0) {
          // offset in z-y plane
          double d = sqrt(pow(spot[i][1] - OpDetPoint[1],2) + pow(spot[i][2] - OpDetPoint[2],2));
          // drift distance (in x)
          double h =  sqrt(pow(spot[i][0] - OpDetPoint[0],2));
          // Solid angle of a disk
          solid_angle_detector = Disk_SolidAngle(d, h, radius);
        }
        else {
           std::cout << "Erorr: Invalid optical detector type." << std::endl;
           exit(1);
        }
      
        // store hits from spot applying weighting factor
        if (width_sum == 0) hits_geo[bin] += 0;    // special case to prevent issues in case of zero width ring (i.e. previous ring already up to boundary)     
        else hits_geo[bin] += ((solid_angle_detector / (2*pi)) * cathode_hits_rec[bin]) * widths[i]/width_sum;  // hits multiplied by weighting factor

      }
    }
  }
  
  // sum total hits
  double hits_geo_total = 0; 
  for (int bin = 0; bin < 9; bin++) {
    hits_geo_total += hits_geo[bin];
  }
  
  // apply corrections to analytic model:
  TVector3 hotspot(plane_depth, ScintPoint[1], ScintPoint[2]);
  double distance_vuv = sqrt(pow(ScintPoint[0] - hotspot[0],2) + pow(ScintPoint[1] - hotspot[1],2) + pow(ScintPoint[2] - hotspot[2],2));
  double distance_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2) + pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
  double distance_fullpath = distance_vuv + distance_vis;
  double cosine_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2)) / distance_vis;
  double theta_vis = acos(cosine_vis)*180./pi;
  int k = (theta_vis/delta_angle);

  // geometric correction
  double hits_rec = VIS_pol_crowns[k]->Eval(distance_vuv)*hits_geo_total/cosine_vis;

  // border effect correction
  // interpolate in x for each y bin
  std::vector<double> interp_vals = {0,0,0,0,0,0};
  for (int i = 0; i < 6; i++){
    interp_vals[i] = interpolate(vDistances_x, VIS_RS60_SP_Crowns_Borders[k][i], std::abs(ScintPoint[0]), false);
  }
  // interpolate in y
  double border_correction = interpolate(vDistances_y, interp_vals, std::abs(ScintPoint[1]), false);
 
  // apply correction
  double hits_rec_borders = border_correction * hits_rec / cosine_vis; 
  
  // round final result
  int hits_vis = std::round(hits_rec_borders);
  
  return hits_vis;
}


// gaisser-hillas function definition
Double_t semi_analytic_hits::GaisserHillas(double *x,double *par) {
  //This is the Gaisser-Hillas function
  Double_t X_mu_0=par[3];
  Double_t Normalization=par[0];
  Double_t Diff=par[1]-X_mu_0;
  Double_t Term=pow((*x-X_mu_0)/Diff,Diff/par[2]);
  Double_t Exponential=TMath::Exp((par[1]-*x)/par[2]);
  
  return ( Normalization*Term*Exponential);
}


// solid angle of rectanglular aperture calculation functions

double semi_analytic_hits::omega(const double &a, const double &b, const double &d) const{

  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));

}

double semi_analytic_hits::solid(const acc& out, const TVector3 &v) const{

  //v is the position of the track segment with respect to 
  //the center position of the arapuca window 
 
  // arapuca plane fixed in x direction	

  if( v.Y()==0.0 && v.Z()==0.0){
    return omega(out.w,out.h,v.X());
  }
  
  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  
  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(b-B),d)+omega(2*A,2*(b-B),d)+omega(2*(a-A),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(b-B),d)-omega(2*A,2*(b-B),d)+omega(2*(A+a),2*B,d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(B+b),d)-omega(2*(a-A),2*B,d)+omega(2*A,2*(B+b),d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}


// solid angle of circular aperture
double semi_analytic_hits::Disk_SolidAngle(double *x, double *p) {
  const double d = x[0];
  const double h = x[1];
  const double b = p[0];
  if(b <= 0. || d < 0. || h <= 0.) return 0.; 
  const double aa = TMath::Sqrt(h*h/(h*h+(b+d)*(b+d)));
  if(d == 0) {
    return 2.*TMath::Pi()*(1.-aa);
  }
  const double bb = TMath::Sqrt(4*b*d/(h*h+(b+d)*(b+d)));
  const double cc = 4*b*d/((b+d)*(b+d));

  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  if(TMath::Abs(ROOT::Math::comp_ellint_1(bb) - bb) < 1e-10 && TMath::Abs(ROOT::Math::comp_ellint_3(cc,bb) - cc) <1e-10) {
    throw(std::runtime_error("please do gSystem->Load(\"libMathMore.so\") before running Disk_SolidAngle for the first time!"));
  }
  if(d < b) {
    return 2.*TMath::Pi() - 2.*aa*(ROOT::Math::comp_ellint_1(bb) + TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb));
  }
  if(d == b) {
    return TMath::Pi() - 2.*aa*ROOT::Math::comp_ellint_1(bb);
  }
  if(d > b) {
    return 2.*aa*(TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb) - ROOT::Math::comp_ellint_1(bb));
  }

  return 0.;
}

double semi_analytic_hits::Disk_SolidAngle(double d, double h, double b) {
  double x[2] = { d, h };
  double p[1] = { b };
  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  return Disk_SolidAngle(x,p);
}

double semi_analytic_hits::interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate ) {
  int size = xData.size();
  int i = 0;                                          // find left end of interval for interpolation
  if ( x >= xData[size - 2] )                         // special case: beyond right end
    {
      i = size - 2;
    }
  else
    {
      while ( x > xData[i+1] ) i++;
    }
  double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1]; // points on either side (unless beyond ends)
  if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
    {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
    }
  double dydx = ( yR - yL ) / ( xR - xL );            // gradient
  return yL + dydx * ( x - xL );                      // linear interpolation
}