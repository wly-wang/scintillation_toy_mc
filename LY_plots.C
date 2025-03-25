#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <sstream>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

void LY(const char *infilename) {

	const double plane_depth = 363.38405;

	// read output file from 7mdrift toy MC
	TFile *filein = new TFile(infilename);

	// set branch addresses 
	TTree *data_tree = nullptr;
	filein->GetObject("data_tree", data_tree);
	int data_event;
	data_tree->SetBranchAddress("data_event", &data_event);

	// TTree *data_tree_vis = nullptr;
	// filein->GetObject("data_tree_vis", data_tree_vis);
	// int data_event_vis;
	// data_tree_vis->SetBranchAddress("data_event_vis", &data_event_vis);

	TTree *data_tree_vuv = nullptr;
	filein->GetObject("data_tree_vuv", data_tree_vuv);
	int data_event_vuv, data_pmt_vuv;
	data_tree_vuv->SetBranchAddress("data_event_vuv", &data_event_vuv);
    data_tree_vuv->SetBranchAddress("data_pmt_vuv", &data_pmt_vuv);

	TTree *event_tree = nullptr;
	filein->GetObject("event_tree", event_tree);
	int event_no; double event_x_pos, event_y_pos, event_z_pos, event_E;
	event_tree->SetBranchAddress("event_no", &event_no);
	event_tree->SetBranchAddress("event_x_pos", &event_x_pos);
	event_tree->SetBranchAddress("event_y_pos", &event_y_pos);
	event_tree->SetBranchAddress("event_z_pos", &event_z_pos);
	event_tree->SetBranchAddress("event_E", &event_E); 

    
	// loop through data trees counting entries to get total hits for each event
	int number_events = event_tree->GetEntries();
	vector<int> hits_both(number_events,0);
	vector<int> hits_vuv(number_events,0);
	vector<int> hits_vuv_anode(number_events, 0);
	vector<int> hits_vuv_cathode(number_events, 0);
	// vector<int> hits_vis(number_events,0);

	int size_data_tree = data_tree->GetEntries();
	for (int i = 0; i < size_data_tree; i++) {
		data_tree->GetEntry(i);
		hits_both[data_event]++;
	}

	cout << "data tree complete" << endl;

	int size_data_tree_vuv = data_tree_vuv->GetEntries();
    int number_of_opDet_anode = 480;  // 0-479 anode 480-912 cathode PDs
	for (int i = 0; i < size_data_tree_vuv; i++) {
		data_tree_vuv->GetEntry(i);
		hits_vuv[data_event_vuv]++;
        if (data_pmt_vuv < number_of_opDet_anode){
			hits_vuv_anode[data_event_vuv]++;
		} else {
			hits_vuv_cathode[data_event_vuv]++;
		}
	}
	cout << "vuv only data tree complete" << endl;

	// int size_data_tree_vis = data_tree_vis->GetEntries();
	// for (int i = 0; i < size_data_tree_vis; i++) {
	// 	data_tree_vis->GetEntry(i);
	// 	hits_vis[data_event_vis]++;
	// }

	// cout << "vis only data tree complete" << endl;


	//vector<double> LY_both;
	//vector<double> LY_vis;
	vector<double> LY_vuv;
    vector<double> LY_vuv_anode;
    vector<double> LY_vuv_cathode;
	vector<double> x_bin;
	vector<double> err_x;
	//vector<double> err_y_both;
	//vector<double> err_y_vis;
	vector<double> err_y_vuv;
    vector<double> err_y_vuv_anode;
    vector<double> err_y_vuv_cathode;


	double bin_step = 30.;
	double drift_size = 718.;
	
	// loop over bins in x, calculating light yield
	for (int bin = 0; bin < drift_size/bin_step; bin++){
		//vector<double> LY_vals_both; LY_vals_both.reserve(1e6);
		vector<double> LY_vals_vuv; LY_vals_vuv.reserve(1e6);
        vector<double> LY_vals_vuv_anode; LY_vals_vuv_anode.reserve(1e6);
        vector<double> LY_vals_vuv_cathode; LY_vals_vuv_cathode.reserve(1e6);
		//vector<double> LY_vals_vis; LY_vals_vis.reserve(1e6);

		double xbin_min = bin*bin_step;
		double xbin_max = (bin+1)*bin_step;
		x_bin.push_back((xbin_max + xbin_min)/2  - 359.);
		err_x.push_back(bin_step/2);
		// loop over events selecting those within current bin
		for (int i = 0; i < number_events; i++){
			event_tree->GetEntry(i);
			// LY of events occuring within bin
			if (event_x_pos > xbin_min - 359. && event_x_pos <= xbin_max - 359.) {
				//double LY_val_both = hits_both[event_no]/event_E;
				//double LY_val_vis = hits_vis[event_no]/event_E;
				double LY_val_vuv = hits_vuv[event_no]/event_E;
				double LY_val_vuv_anode = hits_vuv_anode[event_no]/event_E;
				double LY_val_vuv_cathode = hits_vuv_cathode[event_no]/event_E;

				//LY_vals_both.push_back(LY_val_both);
				LY_vals_vuv.push_back(LY_val_vuv);
				//LY_vals_vis.push_back(LY_val_vis);
				LY_vals_vuv_anode.push_back(LY_val_vuv_anode);
				LY_vals_vuv_cathode.push_back(LY_val_vuv_cathode);
                
			}
		}
		// calculate mean and stdev of LY for this bin for plotting 
		double m_both, m_vuv, m_vuv_anode, m_vuv_cathode; // ,m_vis; 
        double stdev_both, stdev_vuv, stdev_vuv_anode, stdev_vuv_cathode; // ,stdev_vis;		
		
		// // both
		// // mean
		// double sum_both = std::accumulate(std::begin(LY_vals_both), std::end(LY_vals_both), 0.0);
		// m_both =  sum_both / LY_vals_both.size();
		// // stdev
		// double accum_both = 0.0;
		// std::for_each (std::begin(LY_vals_both), std::end(LY_vals_both), [&](const double d) { accum_both += (d - m_both) * (d - m_both); } );
		// stdev_both = sqrt(accum_both / (LY_vals_both.size()-1));

		// vuv
		// mean
		double sum_vuv = std::accumulate(std::begin(LY_vals_vuv), std::end(LY_vals_vuv), 0.0);
		m_vuv =  sum_vuv / LY_vals_vuv.size();
		// stdev
		double accum_vuv = 0.0;
		std::for_each (std::begin(LY_vals_vuv), std::end(LY_vals_vuv), [&](const double d) { accum_vuv += (d - m_vuv) * (d - m_vuv); } );
		stdev_vuv = sqrt(accum_vuv / (LY_vals_vuv.size()-1));

        // vuv anode
        // mean
        double sum_vuv_anode = std::accumulate(std::begin(LY_vals_vuv_anode), std::end(LY_vals_vuv_anode), 0.0);
        m_vuv_anode =  sum_vuv_anode / LY_vals_vuv_anode.size();
        // stdev
        double accum_vuv_anode = 0.0;
        std::for_each (std::begin(LY_vals_vuv_anode), std::end(LY_vals_vuv_anode), [&](const double d) { accum_vuv_anode += (d - m_vuv_anode) * (d - m_vuv_anode); } );
        stdev_vuv_anode = sqrt(accum_vuv_anode / (LY_vals_vuv_anode.size()-1));

        // vuv cathode
        // mean
        double sum_vuv_cathode = std::accumulate(std::begin(LY_vals_vuv_cathode), std::end(LY_vals_vuv_cathode), 0.0);
        m_vuv_cathode =  sum_vuv_cathode / LY_vals_vuv_cathode.size();
        // stdev
        double accum_vuv_cathode = 0.0;
        std::for_each (std::begin(LY_vals_vuv_cathode), std::end(LY_vals_vuv_cathode), [&](const double d) { accum_vuv_cathode += (d - m_vuv_cathode) * (d - m_vuv_cathode); } );
        stdev_vuv_cathode = sqrt(accum_vuv_cathode / (LY_vals_vuv_cathode.size()-1));

		// // vis
		// // mean
		// double sum_vis = std::accumulate(std::begin(LY_vals_vis), std::end(LY_vals_vis), 0.0);
		// m_vis =  sum_vis / LY_vals_vis.size();
		// // stdev
		// double accum_vis = 0.0;
		// std::for_each (std::begin(LY_vals_vis), std::end(LY_vals_vis), [&](const double d) { accum_vis += (d - m_vis) * (d - m_vis); } );
		// stdev_vis = sqrt(accum_vis / (LY_vals_vis.size()-1));

		// // add mean and stdev to vectors for plotting
		// LY_both.push_back(m_both);
		// LY_vis.push_back(m_vis);
		LY_vuv.push_back(m_vuv);
		LY_vuv_anode.push_back(m_vuv_anode);
        LY_vuv_cathode.push_back(m_vuv_cathode);
        
        //err_y_both.push_back(stdev_both);
        // err_y_vis.push_back(stdev_vis);
        err_y_vuv.push_back(stdev_vuv);
        err_y_vuv_anode.push_back(stdev_vuv_anode);
        err_y_vuv_cathode.push_back(stdev_vuv_cathode);

	}

    // create output .root file
    unique_ptr<TFile> myFile(TFile::Open("LY_plots.root", "RECREATE"));
    cout << "An empty .root file created for storing output plots." << endl;
	// create plot
    TCanvas *c1 = new TCanvas("c1","",200,10,1080,1080);
    auto mg = new TMultiGraph();
    auto legend = new TLegend(0.65,0.7,0.9,0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

	// TGraph *gr1 = new TGraphErrors(x_bin.size(), &(x_bin[0]),&(LY_both[0]),&(err_x[0]),&(err_y_both[0]));
	TGraph *gr2 = new TGraphErrors(x_bin.size(), &(x_bin[0]),&(LY_vuv[0]),&(err_x[0]),&(err_y_vuv[0]));
    TGraph *gr3 = new TGraphErrors(x_bin.size(), &(x_bin[0]),&(LY_vuv_anode[0]),&(err_x[0]),&(err_y_vuv_anode[0]));
    TGraph *gr4 = new TGraphErrors(x_bin.size(), &(x_bin[0]),&(LY_vuv_cathode[0]),&(err_x[0]),&(err_y_vuv_cathode[0]));
	// TGraph *gr3 = new TGraphErrors(x_bin.size(), &(x_bin[0]),&(LY_vis[0]),&(err_x[0]),&(err_y_vis[0]));

	//gr1->SetMarkerColor(1);
	gr2->SetMarkerColor(1);
    gr3->SetMarkerColor(2);
    gr4->SetMarkerColor(4);
	//gr3->SetMarkerColor(2);

	//gr1->SetLineColor(1);
	gr2->SetLineColor(1);
    gr3->SetLineColor(2);
    gr4->SetLineColor(4);
	//gr3->SetLineColor(2);

	//gr1->SetMarkerStyle(20);
	gr2->SetMarkerStyle(22);
    gr3->SetMarkerStyle(23);
    gr4->SetMarkerStyle(24);
	//gr3->SetMarkerStyle(23);

	///gr1->SetMarkerSize(2);
	gr2->SetMarkerSize(2.5);
    gr3->SetMarkerSize(1);
    gr4->SetMarkerSize(1);
	//gr3->SetMarkerSize(2);

	//legend->AddEntry(gr1, "Total" ,"p");
	legend->AddEntry(gr2, "Direct" ,"p");
    legend->AddEntry(gr3, "Anode" ,"p");
    legend->AddEntry(gr4, "Cathode" ,"p");
	//legend->AddEntry(gr3, "Reflected" ,"p");

	//mg->Add(gr1);
	mg->Add(gr2);
    mg->Add(gr3);
    mg->Add(gr4);
	//mg->Add(gr3);

	// draw graph
	mg->GetXaxis()->SetTitle("x [cm]");
	mg->GetYaxis()->SetTitle("Light Yield [PE/MeV]");
	mg->GetHistogram()->SetTitle("");

	mg->GetXaxis()->SetRangeUser(-359.,359.);
	mg->GetYaxis()->SetRangeUser(0.,50.);
	mg->GetYaxis()->SetTitleSize(0.045);
    mg->GetXaxis()->SetTitleSize(0.045);

	mg->Draw("ap");
	legend->Draw("same");
	c1->Update();

	TLatex ltx;
	// ltx.SetTextSize(0.036);
	// ltx.DrawLatex(100,42,"#splitline{Wavelength-shifting}{reflector foils}");
    c1 -> Print("plots.pdf");
  	
}