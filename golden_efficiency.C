#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <algorithm>
#include <tuple>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

int golden_efficiency() {  	
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        std::make_tuple("nu mc", "1",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 0.13011 * 2.0), // *2 because sample is only half size
        std::make_tuple("nu mc", "2",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root", 0.25750 * 2.0),
        std::make_tuple("nu mc", "3",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root", 0.20113 * 2.0),
        std::make_tuple("nu mc", "4bcd", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root", 0.13074 * 2.0),
        std::make_tuple("nu mc", "5",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196 * 2.0),
    };

    TH1D *generic_selected_signal = new TH1D("h1", "True Pion Momentum", 30, 0.1, 1);
    TH1D *golden_selected_signal = new TH1D("h2", "True Pion Momentum", 30, 0.1, 1);
    TH1D *all_signal = new TH1D("h3", "True Pion Momentum", 30, 0.1, 1);

    TH1D *true_golden_signal = new TH1D("h4", "True Pion Momentum", 30, 0.1, 1);
    TH1D *true_non_golden_signal = new TH1D("h5", "True Pion Momentum", 30, 0.1, 1);

	for (auto file : files) {
        // std::string type = std::get<0>(file);
        // std::string run = std::get<1>(file);
        const std::string filePath = std::get<2>(file);
        const float runWeight = std::get<3>(file);
		std::cout << "Processing file: " << filePath << std::endl;

		TFile *f1 = new TFile(filePath.c_str());
		if(!f1 || !f1->IsOpen()) {
			std::cout << "Could not open input file!" << std::endl; 
			exit(1);
		}

		TTree *t1 = (TTree*)f1->Get("stv_tree");
		if(!t1) {
			std::cout << "Could not find TTree 'nuselection/NeutrinoSelectionFilter' in file!" << std::endl; 
			exit(1);
		}

		// Disable all branches
        t1->SetBranchStatus("*", 0);
		// Activate only the branches we need
		t1->SetBranchStatus("cc1pi_signal", 1);
        t1->SetBranchStatus("cc1pi_selected_generic", 1);
        t1->SetBranchStatus("cc1pi_selected_golden", 1);
		t1->SetBranchStatus("cc1pi_truth_pionMomentum", 1);
        t1->SetBranchStatus("spline_weight", 1);
        t1->SetBranchStatus("tuned_cv_weight", 1);
        t1->SetBranchStatus("category", 1);

        // set addresses for variables
        Bool_t cc1pi_signal = 0;
        Bool_t cc1pi_selected_generic = 0;
        Bool_t cc1pi_selected_golden = 0;
        Float_t cc1pi_truth_pionMomentum = 0;
        Float_t spline_weight = 0;
        Float_t tuned_cv_weight = 0;
        Int_t category = 0;
        t1->SetBranchAddress("cc1pi_signal", &cc1pi_signal);
        t1->SetBranchAddress("cc1pi_selected_generic", &cc1pi_selected_generic);
        t1->SetBranchAddress("cc1pi_selected_golden", &cc1pi_selected_golden);
        t1->SetBranchAddress("cc1pi_truth_pionMomentum", &cc1pi_truth_pionMomentum);
        t1->SetBranchAddress("spline_weight", &spline_weight);
        t1->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);
        t1->SetBranchAddress("category", &category);

		// loop over events
		// std::cout << "WARNING! Only using 20\% of events for testing purposes" << std::endl;
		int n_entries = t1->GetEntries();
		std::cout << "Number events: " << n_entries << std::endl;

		for (int e = 0; e < n_entries; e++) {
			// get current entry
			t1->GetEntry(e);
            auto eventWeight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;
            eventWeight *= runWeight; // Apply the run weight

			if ( n_entries >= 10 &&  (e % (n_entries/10) == 0) ) {
				std::cout << Form("%i0%% Completed...\n", e / (n_entries/10));
			}

            // Fill the histograms with the true pion momentum is they are cc1pi_signal and their true pion momentum is > 0.1 
            if (cc1pi_signal && cc1pi_truth_pionMomentum > 0.1) {
                all_signal->Fill(cc1pi_truth_pionMomentum, eventWeight);
                if (category == 1)
                {
                    true_non_golden_signal->Fill(cc1pi_truth_pionMomentum, eventWeight);
                } else if (category == 2)
                {
                    true_golden_signal->Fill(cc1pi_truth_pionMomentum, eventWeight);
                } else throw std::runtime_error("Category is not 1 or 2");

                if (cc1pi_selected_generic) {
                    generic_selected_signal->Fill(cc1pi_truth_pionMomentum, eventWeight);
                }
                if (cc1pi_selected_golden) {
                    golden_selected_signal->Fill(cc1pi_truth_pionMomentum, eventWeight);
                }
            }

		}

		f1->Close();
		delete f1;
	}

    // Create a final histogram by taking the ratio of generic_selected_signal and golden_selected_signal with all_signal and plotting these two ratios on one plot
    TH1D *generic_selected_signal_ratio = (TH1D*)generic_selected_signal->Clone("generic_selected_signal_ratio");
    generic_selected_signal_ratio->Divide(all_signal);
    TH1D *golden_selected_signal_ratio = (TH1D*)golden_selected_signal->Clone("golden_selected_signal_ratio");
    golden_selected_signal_ratio->Divide(all_signal);

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    gStyle->SetOptStat(0); // Remove the stats box
    generic_selected_signal_ratio->SetLineColor(kBlack);
    // generic_selected_signal_ratio->SetLineStyle(2); // Dashed line
    generic_selected_signal_ratio->SetLineWidth(2); // Make line thicker
    golden_selected_signal_ratio->SetLineColor(kBlue);
    generic_selected_signal_ratio->GetYaxis()->SetRangeUser(0, generic_selected_signal_ratio->GetMaximum()*1.3); // Start y-axis at 0 and scale to 10% above max
    golden_selected_signal_ratio->SetLineWidth(2); // Make line thicker
    generic_selected_signal_ratio->SetTitle("");
    generic_selected_signal_ratio->Draw("HIST");
    golden_selected_signal_ratio->Draw("HIST same");

    generic_selected_signal_ratio->GetXaxis()->SetTitle("p_{#pi^{#pm}} (GeV/c)"); // Set x axis label
    generic_selected_signal_ratio->GetXaxis()->CenterTitle(true); // Center the title horizontally
    generic_selected_signal_ratio->GetXaxis()->SetTitleOffset(1.1); // Adjust the position to center it horizontally
    generic_selected_signal_ratio->GetYaxis()->SetTitle("Signal Selection Efficiency"); // Set y axis label

    // // Add vertical lines to indicate the selection cuts
    // double linePositions[] = {0.1, 0.16, 0.19, 0.22, 0.6};
    // for (double position : linePositions) {
    //     TLine *line = new TLine(position, 0, position, generic_selected_signal_ratio->GetMaximum()*1.3);
    //     line->SetLineColor(kBlack);
    //     line->SetLineStyle(7); // Dotted line
    //     line->Draw();
    // }

    // Add a legend
    TLegend* legend = new TLegend(0.6, 0.7, 0.95, 0.9);
    legend->SetTextSize(0.04); // Increase font size of the legend title and entries
    legend->SetHeader("Selection", "C"); // Centered title
    legend->AddEntry(generic_selected_signal_ratio, "General", "l");
    legend->AddEntry(golden_selected_signal_ratio, "Unscattered Enhanced", "l");
    legend->Draw();

    // Save the plot as pdf and .C
    c1->SaveAs("plots/generic_vs_golden_efficiency.pdf");
    c1->SaveAs("plots/generic_vs_golden_efficiency.C");    


    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);

    gStyle->SetOptStat(0); // Remove the stats box
    generic_selected_signal->SetLineColor(kBlack);
    generic_selected_signal->SetLineWidth(2); // Make line thicker
    golden_selected_signal->SetLineColor(kBlue);
    golden_selected_signal->SetLineWidth(2); // Make line thicker
    // all_signal->SetLineColor(kRed); // Set color for all_signal
    // all_signal->SetLineWidth(2); // Make line thicker for all_signal
    generic_selected_signal->GetYaxis()->SetRangeUser(0, generic_selected_signal->GetMaximum()*1.1); // Start y-axis at 0 and scale to 30% above max
    generic_selected_signal->SetTitle("");
    generic_selected_signal->Draw("HIST");
    golden_selected_signal->Draw("HIST same");
    // all_signal->Draw("HIST same"); // Draw all_signal on the same canvas
    generic_selected_signal->GetXaxis()->SetTitle("True p_{#pi^{#pm}} (GeV/c)"); // Set x axis label

    generic_selected_signal->GetXaxis()->CenterTitle(true); // Center the title horizontally
    generic_selected_signal->GetXaxis()->SetTitleOffset(1.1); // Adjust the position to center it horizontally
    generic_selected_signal->GetYaxis()->SetTitle("Predicted Signal Events"); // Set y axis label

    // Add a legend
    TLegend* legend2 = new TLegend(0.5, 0.7, 0.9, 0.9);
    legend2->SetTextSize(0.04); // Increase font size of the legend title and entries
    legend2->SetHeader("Selection", "C"); // Centered title
    // legend2->AddEntry(all_signal, "All", "l"); // Add all_signal to the legend
    legend2->AddEntry(generic_selected_signal, "General", "l");
    legend2->AddEntry(golden_selected_signal, "Unscattered Enhanced", "l");
    legend2->Draw();

    // Save the plot as pdf and .C
    c2->SaveAs("plots/generic_vs_golden_event_rate.pdf");
    c2->SaveAs("plots/generic_vs_golden_event_rate.C");


    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    
    gStyle->SetOptStat(0); // Remove the stats box
    
    // Set line colors and increase line thickness
    true_golden_signal->SetLineColor(kGreen);
    true_golden_signal->SetLineWidth(4); // Increase line thickness
    true_non_golden_signal->SetLineColor(kGreen+3);
    true_non_golden_signal->SetLineWidth(5); // Increase line thickness
    true_non_golden_signal->SetLineStyle(2); // Dashed line
    
    // Scale y-axis
    true_golden_signal->GetYaxis()->SetRangeUser(0, std::max(true_golden_signal->GetMaximum(), true_non_golden_signal->GetMaximum()) * 1.1);
    
    // Remove the title
    true_golden_signal->SetTitle("");
    
    // Increase the axis label and title font size
    true_golden_signal->GetXaxis()->SetLabelSize(0.05);
    true_golden_signal->GetYaxis()->SetLabelSize(0.05);
    true_golden_signal->GetXaxis()->SetTitleSize(0.05);
    true_golden_signal->GetYaxis()->SetTitleSize(0.05);
    
    // Set axis labels and center the x-axis title
    true_golden_signal->GetXaxis()->SetTitle("True p_{#pi^{#pm}} (GeV/c)");
    true_golden_signal->GetXaxis()->CenterTitle(true);
    true_golden_signal->GetXaxis()->SetTitleOffset(1.1);
    true_golden_signal->GetYaxis()->SetTitle("Predicted Signal Events");
    
    // Adjust margins to ensure everything fits on screen
    c3->SetLeftMargin(0.15);
    c3->SetRightMargin(0.05);
    c3->SetTopMargin(0.1);
    c3->SetBottomMargin(0.15);
    
    // Draw histograms
    true_golden_signal->Draw("HIST");
    true_non_golden_signal->Draw("HIST same"); // Draw on the same canvas
    
    // Add a legend
    TLegend* legend3 = new TLegend(0.5, 0.7, 0.9, 0.9);
    legend3->SetTextSize(0.04); // Increase font size of the legend title and entries
    legend3->SetHeader("Signal Type", "C"); // Centered title
    legend3->AddEntry(true_golden_signal, "Unscattered pion", "l");
    legend3->AddEntry(true_non_golden_signal, "Reinteracting pion", "l");
    legend3->Draw();
    
    // Save the plot as pdf and .C
    c3->SaveAs("plots/generic_vs_golden_truth_total.pdf");
    c3->SaveAs("plots/generic_vs_golden_truth_total.C");



	return 0;
}