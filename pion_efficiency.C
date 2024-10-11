#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <algorithm>
#include <tuple>

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"


float rangeFitFunc(float r) {
	// Pion from range fit
	const float a = 0.25798;
	const float b = 0.0024088;
	const float c = 0.18828;
	const float d = 0.11687;
	// Work out the maximum straight-line range a particle can have, so we can set sensible limits on the fit
	const auto minRange = 0.f;
	const float lowX  = 0.f;
	const float highX = 256.35f;
	const float lowY  = -116.5f;
	const float highY = 116.5f;
	const float lowZ  = 0.f;
	const float highZ = 1036.8f;
	const auto maxRange = std::pow( std::pow(highX - lowX, 2) + std::pow(highY - lowY, 2) + std::pow(highZ - lowZ, 2) , 0.5f);

	// std::shared_ptr<TF1> pFunc(new TF1(("fitFunc_0", "[0] + [1]*x - [2]*pow(x, -[3])", minRange, maxRange));

	if(r<minRange || r>maxRange) return -1.f;
	return a + b * r - c * std::pow(r, -d);
};

void MakePlot2D(TH2F* hist, const std::string& name, const bool drawDiagonal = false, const bool axisTicks = true, const bool drawText = false)
{
    // Set the color and line thickness of the histograms
    hist->SetLineColor(kGreen);
    hist->SetLineWidth(3);

    // Remove the stats box
    hist->SetStats(0);

    // Draw the histogram
    TCanvas* c1 = new TCanvas(name.c_str(), name.c_str(), 800, 600);

    if(!axisTicks)
    {
        hist->GetXaxis()->SetNdivisions(hist->GetNbinsX(), kFALSE);
        hist->GetYaxis()->SetNdivisions(hist->GetNbinsY(), kFALSE);
    }

    if (drawText)
    {
        // Set the text format to 2 digits after the decimal point
        gStyle->SetPaintTextFormat("4.2f");
        // Increase the text size
        gStyle->SetTextSize(0.2);
        // Draw the histogram with text
        hist->Draw("COLZ TEXT");
    }
    else
    {
        hist->Draw("COLZ");
    }

    if (drawDiagonal)
    {
        // Add a diagonal line
        double x1 = hist->GetXaxis()->GetXmin();
        double y1 = hist->GetYaxis()->GetXmin();
        double x2 = hist->GetXaxis()->GetXmax();
        double y2 = hist->GetYaxis()->GetXmax();
        TLine *line = new TLine(x1, y1, x2, y2);
        line->SetLineColor(kBlack);  // Set line color to black
        line->SetLineStyle(2);       // Set line style to dashed
        line->SetLineWidth(2);       // Set line width
        line->Draw("same");          // Draw line on the same canvas
    }

    // Add more space on the right side of the plot for the z-axis labels
    gStyle->SetPadRightMargin(0.15);

    c1->SaveAs(("plots/pion_efficiency_" + name + "_testingOnly_lowPiMomThreshold.pdf").c_str());
	c1->SaveAs(("plots/pion_efficiency_" + name + "_testingOnly_lowPiMomThreshold.png").c_str());
	c1->SaveAs(("plots/pion_efficiency_" + name + "_testingOnly_lowPiMomThreshold.C").c_str());

    // Delete the TCanvas
    delete c1;
}

int pion_efficiency() {  	
	// List of files
	const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/1March24/";
	// tuple: type, run, file path, run weight
	const std::vector<std::string> files = {
		"/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu.root",
		"/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_73_run4a_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_73_run4b_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4c_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4d_nu.root",
		// "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5.root"
	};

	// create a 2D histogram
	TH2F *h2 = new TH2F("h2", "True vs Reco Pion Momentum", 30, 0, 0.6, 30, 0, 0.6);
	double binEdges[] = {-1.1, 0, 0.1, 0.6};
	TH2F *h2_min_bins = new TH2F("h2_min_bins", "True vs Reco Pion Momentum", 3, binEdges, 3, binEdges);

	TH2F *h3 = new TH2F("h3", "Reco Pion Momentum vs Purity; Reco Pion Momentum; Purity", 50, 0, 0.5, 10, 0, 1);
	TH2F *h4 = new TH2F("h4", "Reco Pion Momentum vs Completeness; Reco Pion Momentum; Completeness", 50, 0, 0.5, 10, 0, 1);
	TH2F *h5 = new TH2F("h5", "True Pion Momentum vs Purity; True Pion Momentum; Purity", 50, 0, 0.5, 10, 0, 1);
	TH2F *h5_min_bins = new TH2F("h5_min_bins", "True Pion Momentum vs Purity; True Pion Momentum; Purity", 50, 0, 0.5, 4, 0, 1);
	TH2F *h6 = new TH2F("h6", "True Pion Momentum vs Completeness; True Pion Momentum; Completeness", 50, 0, 0.5, 10, 0, 1);
	TH2F *h7 = new TH2F("h7", "True Pion Momentum vs Reco Pion Momentum; True Pion Momentum; Reco Pion Momentum", 50, 0, 0.5, 50, 0, 0.5);


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

		// Disable all branches
        t1->SetBranchStatus("*", 0);
		// Activate only the branches we need
		t1->SetBranchStatus("mc_pdg", 1);
		t1->SetBranchStatus("backtracked_*", 1);
		t1->SetBranchStatus("mc_E", 1);
		t1->SetBranchStatus("backtracked_completeness", 1);
		t1->SetBranchStatus("mc_px", 1); // added mc_px
		t1->SetBranchStatus("mc_py", 1); // added mc_py
		t1->SetBranchStatus("mc_pz", 1); // added mc_pz
		t1->SetBranchStatus("mc_completeness", 1); // added mc_completeness
		t1->SetBranchStatus("mc_purity", 1); // added mc_purity
		t1->SetBranchStatus("trk_len_v", 1); // added mc_purity

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
		std::vector<float> *trk_len_v = nullptr; // added mc_purity
		std::vector<float> *backtracked_px = nullptr;
		std::vector<float> *backtracked_py = nullptr;
		std::vector<float> *backtracked_pz = nullptr;
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
		t1->SetBranchAddress("trk_len_v", &trk_len_v); // added mc_purity
		t1->SetBranchAddress("backtracked_px", &backtracked_px);
		t1->SetBranchAddress("backtracked_py", &backtracked_py);
		t1->SetBranchAddress("backtracked_pz", &backtracked_pz);
		// loop over events
		// std::cout << "WARNING! Only using 20\% of events for testing purposes" << std::endl;
		int n_entries = t1->GetEntries();
		std::cout << "Number events: " << n_entries << std::endl;

		std::vector<int> michelCandidates;
		int cc1piCounter = 0;

		for (int e = 0; e < n_entries; e++) {
			// get current entry
			t1->GetEntry(e);

			if ( n_entries >= 10 &&  (e % (n_entries/10) == 0) ) {
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
					break; // to accelerate the process
				}
			}

			// check if there is zero other and exactly one 13 (muon) and one pion (211)
			if (count_other == 0 && count_13 == 1 && count_211 == 1) {

				// int true_pion_index = -1;
				// for (size_t i = 0; i < mc_pdg->size(); i++) {
				// 	if (std::abs((*mc_pdg)[i]) == 211) {
				// 		true_pion_index = i;
				// 		break;
				// 	}
				// }

				// if(true_pion_index == -1) {
				// 	throw std::runtime_error("Error: Missing pion in backtracked_pdg or mc_pdg");
				// }

				// std::cout << "DEBUG Point 8" << std::endl;

				// float px = (*mc_px)[true_pion_index];
				// float py = (*mc_py)[true_pion_index];
				// float pz = (*mc_pz)[true_pion_index];
				// const auto truePionMomentum = std::sqrt(px*px + py*py + pz*pz);

				int backtracked_pion_index = -1;
				float highest_completenss = 0;
				// int match_count = 0;
				for (size_t i = 0; i < backtracked_pdg->size(); i++) {
					if (std::abs((*backtracked_pdg)[i]) == 211) 
					{
						if(backtracked_completeness->at(i) > highest_completenss) {
							highest_completenss = backtracked_completeness->at(i);
							backtracked_pion_index = i;
						}
					}
				}

				if(backtracked_pion_index == -1)
					continue;

				// std::cout << "DEBUG Point 9" << std::endl;
				float recoPionMomentumFromRange = -1;

				if(backtracked_pion_index != -1)
				{
					const float backtrackedPionRange = (*trk_len_v)[backtracked_pion_index];
					recoPionMomentumFromRange = rangeFitFunc(backtrackedPionRange);
				}

				const auto px = (*backtracked_px)[backtracked_pion_index];
				const auto py = (*backtracked_py)[backtracked_pion_index];
				const auto pz = (*backtracked_pz)[backtracked_pion_index];
				const auto truePionMomentum = std::sqrt(px*px + py*py + pz*pz);

				// std::cout << "DEBUG Point 10" << std::endl;

				// // fill a plot with the event
				// h2->Fill(recoPionMomentumFromRange, truePionMomentum);
				// h2_min_bins->Fill(recoPionMomentumFromRange, truePionMomentum);

				h3->Fill(recoPionMomentumFromRange, backtracked_purity->at(backtracked_pion_index));
				h4->Fill(recoPionMomentumFromRange, backtracked_completeness->at(backtracked_pion_index));

				h5->Fill(truePionMomentum, backtracked_purity->at(backtracked_pion_index));
				h5_min_bins->Fill(truePionMomentum, backtracked_purity->at(backtracked_pion_index));
				h6->Fill(truePionMomentum, backtracked_completeness->at(backtracked_pion_index));
				h7->Fill(truePionMomentum, recoPionMomentumFromRange);

				// std::cout << "DEBUG Point 11" << std::endl;
			}
		}

		f1->Close();
		delete f1;
	}

	// // save the histogram
	// TFile *fout = new TFile("plots/pionPurityVsEnergy.root", "RECREATE");
	// h2->Write();

	// // create a canvas
	// TCanvas *c1 = new TCanvas("c1", "True Pion Momentum vs Backtracked Reco Pion Momentum", 800, 600);

	// // remove the stats box
	// gStyle->SetOptStat(0);

	// // label the x and z axis
	// h2->GetXaxis()->SetTitle("Reco Pion Momentum [GeV/c]");
	// h2->GetYaxis()->SetTitle("True Pion Momentum [GeV/c]");

	// // draw the histogram on the canvas
	// h2->Draw("COLZ");

	// // save the canvas as a PNG file
	// c1->SaveAs("plots/truePionMomVsBacktrackedRecoPionMom.pdf");
	// c1->SaveAs("plots/truePionMomVsBacktrackedRecoPionMom.png");
	// c1->SaveAs("plots/truePionMomVsBacktrackedRecoPionMom.C");

	// // fout->Close();
	// delete c1;
	// // delete fout;

	// // The bin sizes in the TH2F are not uniform h2_min_bins
	// // Create a new TH2F version of h2_min_bins with all bins sized as 1
	// TH2F *h2_min_bins_uniform = new TH2F("h2_min_bins_uniform", "True vs Reco Pion Momentum", 3, -1, 2, 3, -1, 2);
	// for (int i = 1; i <= h2_min_bins->GetNbinsX(); i++) {
	// 	for (int j = 1; j <= h2_min_bins->GetNbinsY(); j++) {
	// 		h2_min_bins_uniform->Fill(i, j, h2_min_bins->GetBinContent(i, j));
	// 	}
	// }


	// MakePlot2D(h2_min_bins, "truePionMomVsBacktrackedRecoPionMom_min_bins", false, false, true);
	// MakePlot2D(h2_min_bins_uniform, "truePionMomVsBacktrackedRecoPionMom_min_bins_uniform", false, false, true);

	MakePlot2D(h3, "purityVsRecoPionMom", false, true, false);
	MakePlot2D(h4, "completenessVsRecoPionMom", false, true, false);
	MakePlot2D(h5, "purityVsTruePionMom", false, true, false);
	MakePlot2D(h5_min_bins, "purityVsTruePionMom_min_bins", false, true, false);
	MakePlot2D(h6, "completenessVsTruePionMom", false, true, false);
	MakePlot2D(h7, "truePionMomVsRecoPionMom", true, true, false);

	return 0;
}