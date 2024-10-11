#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <algorithm>
#include <tuple>

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"

int reco_pion_purity() {  	
	// List of files
	const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/1March24/";
	// tuple: type, run, file path, run weight
	const std::vector<std::string> files = {
		"/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_73_run4a_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_73_run4b_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4c_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4d_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5.root"
	};

	// create a 2D histogram
	TH2F *h1 = new TH2F("h1", "True Pion Energy vs Completeness", 50, 0, 0.5, 10, 0, 1);
	TH2F *h2 = new TH2F("h2", "True Pion Energy vs Purity", 50, 0, 0.5, 10, 0, 1);
	

	for (const auto& filePath : files) {
		std::cout << "Processing file: " << filePath << std::endl;

		TFile *f1 = new TFile(filePath.c_str());
		if(!f1 || !f1->IsOpen()) {
			std::cout << "Could not open input file!" << std::endl; 
			exit(1);
		}

		TTree *t1 = (TTree*)f1->Get("nuselection/NeutrinoSelectionFilter");
		if(!t1) {
			std::cout << "Could not find TTree 'nuselection/NeutrinoSelectionFilter' in file!" << std::endl; 
			exit(1);
		}

		// set addresses for variables
		std::vector<int> *mc_pdg = nullptr;
		std::vector<float> *backtracked_purity = nullptr;
		std::vector<int> *backtracked_pdg = nullptr;
		std::vector<float> *mc_E = nullptr;
		std::vector<float> *backtracked_completeness = nullptr;
		std::vector<float> *mc_px = nullptr; // added mc_px
		std::vector<float> *mc_py = nullptr; // added mc_py
		std::vector<float> *mc_pz = nullptr; // added mc_pz
		std::vector<float> *mc_completeness = nullptr; // added mc_completeness
		std::vector<float> *mc_purity = nullptr; // added mc_purity
		t1->SetBranchAddress("mc_pdg", &mc_pdg);
		t1->SetBranchAddress("backtracked_purity", &backtracked_purity);
		t1->SetBranchAddress("backtracked_pdg", &backtracked_pdg);
		t1->SetBranchAddress("mc_E", &mc_E);
		t1->SetBranchAddress("backtracked_completeness", &backtracked_completeness);
		t1->SetBranchAddress("mc_px", &mc_px); // added mc_px
		t1->SetBranchAddress("mc_py", &mc_py); // added mc_py
		t1->SetBranchAddress("mc_pz", &mc_pz); // added mc_pz
		t1->SetBranchAddress("mc_completeness", &mc_completeness); // added mc_completeness
		t1->SetBranchAddress("mc_purity", &mc_purity); // added mc_purity
		// loop over events
		// std::cout << "WARNING! Only using 5\% of events for testing purposes" << std::endl;
		int n_entries = t1->GetEntries();
		std::cout << "Number events: " << n_entries << std::endl;

		for (int e = 0; e < n_entries; e++) {
			// get current entry
			t1->GetEntry(e);

			if ( (e != 0) && (n_entries >= 10) &&  (e % (n_entries/10) == 0) ) {
				std::cout << Form("%i0%% Completed...\n", e / (n_entries/10));
			}

			// count the number of entries that have an absolute value of 211, 13, 2212, and other
			int count_211 = 0;
			int count_13 = 0;
			int count_2212 = 0;
			int count_other = 0;

			for (int pdg : *mc_pdg) {
				int abs_pdg = std::abs(pdg);
				if (abs_pdg == 211) {
					count_211++;
				} else if (abs_pdg == 13) {
					count_13++;
				} else if (abs_pdg == 2212) {
					count_2212++;
				} else if (abs_pdg == 11 || abs_pdg == 22 || abs_pdg == 111 || abs_pdg == 321) {
					count_other++;
				}
			}

			// check if there is zero other and exactly one 13 (muon) and one pion (211)
			if (count_other == 0 && count_13 == 1 && count_211 == 1) {
				int reco_pion_index = -1;
				double max_completeness = -1;
				std::cout << "Found muon and pion" << std::endl;
				for (size_t i = 0; i < backtracked_pdg->size(); i++) {
					if (std::abs((*mc_pdg)[i]) == 211) {
						// std::cout << "Found pion with completeness: " << (*mc_completeness)[i] << " and purity: " << (*mc_purity)[i] << std::endl;
						if ((*backtracked_completeness)[i] > max_completeness) {
							max_completeness = (*backtracked_completeness)[i];
							reco_pion_index = i;
						}
						if(reco_pion_index != -1) throw std::runtime_error("More than one pion found in mc_pdg");
						reco_pion_index = i;
					}
				}

				// if (reco_pion_index == -1) {
				// 	continue;
				// }

				float pionMomentum = -1;
				float purity = -1;
				float completeness = -1;
				if (reco_pion_index != -1) {
					float px = (*mc_px)[reco_pion_index];
					float py = (*mc_py)[reco_pion_index];
					float pz = (*mc_pz)[reco_pion_index];
					pionMomentum = std::sqrt(px*px + py*py + pz*pz);
					purity = (*mc_purity)[reco_pion_index];
					completeness = (*mc_completeness)[reco_pion_index];
				}
				else
				{
					throw std::runtime_error("No pion found in mc_pdg");
				}

				// fill a plot with the event
				h1->Fill(pionMomentum, completeness);
				h2->Fill(pionMomentum, purity);
			}
		}

		f1->Close();
		delete f1;
	}

	// // save the histogram
	// TFile *fout = new TFile("plots/pionPurityVsEnergy.root", "RECREATE");
	// h2->Write();

	// create a canvas
	TCanvas *c1 = new TCanvas("c1", "True Pion Energy vs Purity", 800, 600);

	// remove the stats box
	gStyle->SetOptStat(0);

	// label the x and z axis
	h2->GetXaxis()->SetTitle("True Pion Energy [GeV]");
	h2->GetZaxis()->SetTitle("Purity");

	// draw the histogram on the canvas
	h2->Draw("COLZ");

	// save the canvas as a PNG file
	c1->SaveAs("plots/pionPurityVsEnergy.png");

	// fout->Close();
	delete c1;
	// delete fout;


	TCanvas *c2 = new TCanvas("c2", "True Pion Energy vs Completeness", 800, 600);

	// remove the stats box
	gStyle->SetOptStat(0);

	// label the x and z axis
	h1->GetXaxis()->SetTitle("True Pion Energy [GeV]");
	h1->GetZaxis()->SetTitle("Completeness");
	
	// draw the histogram on the canvas#
	h1->Draw("COLZ");

	// save the canvas as a PNG file
	c2->SaveAs("plots/pionCompletenessVsEnergy.png");

	delete c2;

	return 0;
}