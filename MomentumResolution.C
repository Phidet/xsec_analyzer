#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>


void FillHistogram(TTree* tree, TH1D* hist, const std::string& variable, const std::string& condition, const float runWeight)
{

    TH1D* hist_tmp = (TH1D*)hist->Clone("hist");

    hist_tmp->Reset();

    std::string hist_tmp_name = hist_tmp->GetName();
    std::string weight = "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)*" + std::to_string(runWeight);

    // Fill the histograms with weights
    tree->Project(hist_tmp_name.c_str(), variable.c_str(), ("(" + condition + ")*" + weight).c_str());

    hist_tmp->Sumw2();

    // Add the tmp histograms to the main histograms
    hist->Add(hist_tmp);

    delete hist_tmp;
}

void MakePlot(TH1D* hist, const std::string& name)
{

    // Set the color and line thickness of the histograms
    hist->SetLineColor(kGreen);
    hist->SetLineWidth(3);

    // hist->SetStats(0);

    // Draw the histogram
    TCanvas* c1 = new TCanvas(("c1_" + name).c_str(), ("c1_" + name).c_str(), 800, 500);
    hist->Draw("E hist");

    // Add a legend
    // TLegend* legend = new TLegend(0.2, 0.1);
    // legend->AddEntry(hist, "CC0pi", "l");
    // legend->AddEntry(cc1pi_hist, "CC1pi", "l");
    // legend->Draw();

    c1->SaveAs(("plots/MomentumResolution_" + name + "_testingOnly_lowPiMomThreshold.pdf").c_str());

    // Delete the TCanvas
    delete c1;
}


// Create plots with fitted Gaussian drawn on top
auto MakePlotWithFit = [](TH1D* hist, const std::string& name) {
    // Set the color and line thickness of the histograms
    hist->SetLineColor(kGreen);
    hist->SetLineWidth(3);

    // Draw the histogram
    TCanvas* c = new TCanvas(("c_fit_" + name).c_str(), ("c_fit_" + name).c_str(), 800, 500);
    hist->Draw("E hist");

    // Fit a Gaussian to the histogram
    TF1* fitFunc = new TF1(("fitFunc_" + name).c_str(), "gaus", -1, 1);
    hist->Fit(fitFunc, "Q");
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("same");
    std::cout << "Fit parameters for " << name << " Mean: " << fitFunc->GetParameter(1) << " Sigma: " << fitFunc->GetParameter(2) << std::endl;

    // Save the canvas
    c->SaveAs(("plots/MomentumResolution_" + name + "_withFit.pdf").c_str());

    // Delete the TCanvas and fit function
    delete c;
    delete fitFunc;
};

void MomentumResolution() 
{
    // // List of files
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        std::make_tuple("nu mc", "1",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
        std::make_tuple("nu mc", "2",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root", 0.25750*2.0),
        std::make_tuple("nu mc", "3",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root", 0.20113*2.0),
        std::make_tuple("nu mc", "4bcd", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root", 0.13074*2.0),
        std::make_tuple("nu mc", "5",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196*2.0),
    };

    // Define a 1D histogram for contained muon momentum resolution and three for pion momentum resolution (all, scattered, unscattered)
    TH1D* h_muMom_contained_res = new TH1D("h_muMom_contained_res", "Contained #mu Momentum Resolution;#frac{p_{#mu}^{reco} - p_{#mu}^{true}}{p_{#mu}^{true}};Events", 99, -1, 1);
    TH1D* h_muMom_uncontained_res = new TH1D("h_muMom_uncontained_res", "Uncontained #mu Momentum Resolution;#frac{p_{#mu}^{reco} - p_{#mu}^{true}}{p_{#mu}^{true}};Events", 99, -1, 1);
    TH1D* h_piMom_res = new TH1D("h_piMom_res", "All #pi Momentum Resolution;#frac{p_{#pi}^{reco} - p_{#pi}^{true}}{p_{#pi}^{true}};Events", 99, -1, 1);
    TH1D* h_piMom_unscattered_res = new TH1D("h_piMom_unscattered_res", "Unscattered #pi Momentum Resolution;#frac{p_{#pi}^{reco} - p_{#pi}^{true}}{p_{#pi}^{true}};Events", 99, -1, 1);
    TH1D* h_piCosTheta_res = new TH1D("h_piCosTheta_res", "#cos(#theta_{#pi}) Resolution;#frac{#cos(#theta_{#pi}^{reco}) - #cos(#theta_{#pi}^{true})}{#cos(#theta_{#pi}^{true})};Events", 99, -1, 1);


    // Use this to calculate the pion and muon momentum resolution values (reco - true)/true
    // cc1pi_truth_pionMomentum
    // cc1pi_reco_pionMomentum
    // cc1pi_truth_muonMomentum
    // cc1pi_reco_muonMomentum

    // Use this condition to select which events to include (this is for pion momentum all; I'll add the remaining three later):
    // cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1

    // Plot the distribution of the (reco - true)/true values for all events for the given selection filling a histogram
    // Loop over the files
    for (int f = 0; f < files.size(); f++)
    {
        const auto [sampleType, run, filePath, runWeight] = files.at(f);
        std::cout<<"DEBUG - filePath: "<<filePath<<std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(filePath.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");

        // Fill the histograms
        // FillHistogram(tree, h_muMom_contained_res, "(cc1pi_reco_muonMomentum - cc1pi_truth_muonMomentum)/cc1pi_truth_muonMomentum", "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truthMuon_IsContained", runWeight);
        // FillHistogram(tree, h_piMom_res, "(cc1pi_reco_pionMomentum - cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum", "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1", runWeight);
        // FillHistogram(tree, h_piMom_scattered_res, "(cc1pi_reco_pionMomentum - cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum", "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && !true_golden_cc1pi", runWeight);
        // FillHistogram(tree, h_piMom_unscattered_res, "(cc1pi_reco_pionMomentum - cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum", "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && true_golden_cc1pi", runWeight);

        FillHistogram(tree, h_muMom_contained_res, "(cc1pi_reco_muonMomentum - cc1pi_truth_muonMomentum)/cc1pi_truth_muonMomentum", "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained && (cc1pi_backtracked_muonPDG == 13 || cc1pi_backtracked_muonPDG == -13)", runWeight);
        FillHistogram(tree, h_muMom_uncontained_res, "(cc1pi_reco_muonMomentum - cc1pi_truth_muonMomentum)/cc1pi_truth_muonMomentum", "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained && (cc1pi_backtracked_muonPDG == 13 || cc1pi_backtracked_muonPDG == -13)", runWeight);
        FillHistogram(tree, h_piMom_res, "(cc1pi_reco_pionMomentum - cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum", "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && (cc1pi_backtracked_protonPDG == 211 || cc1pi_backtracked_protonPDG == -211)", runWeight); // This variable name is wrong and protn containes the information about the pion
        FillHistogram(tree, h_piMom_unscattered_res, "(cc1pi_reco_pionMomentum - cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum", "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_selected_golden && cc1pi_reco_pionMomentum >= 0.1 && (cc1pi_backtracked_protonPDG == 211 || cc1pi_backtracked_protonPDG == -211)", runWeight);
        FillHistogram(tree, h_piCosTheta_res, "(cc1pi_reco_pionCosTheta - cc1pi_truth_pionCosTheta)/cc1pi_truth_pionCosTheta", "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && (cc1pi_backtracked_protonPDG == 211 || cc1pi_backtracked_protonPDG == -211)", runWeight);

        tFile->Close();
    }


    // Make the plots with fitted Gaussian
    MakePlotWithFit(h_muMom_contained_res, "muMom_contained_res");
    MakePlotWithFit(h_muMom_uncontained_res, "muMom_uncontained_res");
    MakePlotWithFit(h_piMom_res, "piMom_res");
    MakePlotWithFit(h_piMom_unscattered_res, "piMom_unscattered_res");
    MakePlotWithFit(h_piCosTheta_res, "piCosTheta_res");

    // Make the plots
    MakePlot(h_muMom_contained_res, "muMom_contained_res");
    MakePlot(h_muMom_uncontained_res, "muMom_uncontained_res");
    MakePlot(h_piMom_res, "piMom_res");
    MakePlot(h_piMom_unscattered_res, "piMom_unscattered_res");
    MakePlot(h_piCosTheta_res, "piCosTheta_res");

    // Create a custom plot showing both h_piMom_res and h_piMom_unscattered_res as line plots without fill
    TCanvas* c2 = new TCanvas("c2_piMom_res_comparison", "c2_piMom_res_comparison", 800, 500);

    // Adjust margins to ensure everything fits on screen
    c2->SetLeftMargin(0.17);
    c2->SetRightMargin(0.05);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.15);

    // Subtract h_piMom_unscattered_res from h_piMom_res
    TH1D* h_piMom_res_subtracted = (TH1D*)h_piMom_res->Clone("h_piMom_res_subtracted");
    h_piMom_res_subtracted->Add(h_piMom_unscattered_res, -1);

    // Set properties for the subtracted histogram
    h_piMom_res_subtracted->SetLineColor(kBlue);
    h_piMom_res_subtracted->SetLineWidth(4);
    h_piMom_res_subtracted->SetLineStyle(2);
    h_piMom_res_subtracted->SetFillStyle(0); // No fill

    // Set properties for the unscattered histogram
    h_piMom_unscattered_res->SetLineColor(kGreen);
    h_piMom_unscattered_res->SetLineWidth(4);
    h_piMom_unscattered_res->SetLineStyle(1);
    h_piMom_unscattered_res->SetFillStyle(0); // No fill

    // Rescale y axis to fit both plots
    double max_y = std::max(h_piMom_res_subtracted->GetMaximum(), h_piMom_unscattered_res->GetMaximum());
    h_piMom_res_subtracted->SetMaximum(max_y * 1.1);
    h_piMom_unscattered_res->SetMaximum(max_y * 1.1);

    // Draw the histograms
    h_piMom_res_subtracted->Draw("hist");
    h_piMom_unscattered_res->Draw("hist same");

    // Set axis titles and properties
    h_piMom_res_subtracted->SetTitle("");
    h_piMom_res_subtracted->GetXaxis()->SetTitle("(p_{#pi}^{reco} - p_{#pi}^{truth})/p_{#pi}^{truth}");
    h_piMom_res_subtracted->GetYaxis()->SetTitle("#splitline{Selected Simulated Signal Events}{with Correct Pion Candidate}");
    h_piMom_res_subtracted->GetXaxis()->SetLabelSize(0.05);
    h_piMom_res_subtracted->GetYaxis()->SetLabelSize(0.05);
    h_piMom_res_subtracted->GetXaxis()->SetTitleSize(0.05);
    h_piMom_res_subtracted->GetYaxis()->SetTitleSize(0.05);
    h_piMom_res_subtracted->GetXaxis()->SetLabelOffset(0.012);
    h_piMom_res_subtracted->GetYaxis()->SetLabelOffset(0.012);
    h_piMom_res_subtracted->GetXaxis()->CenterTitle();
    h_piMom_res_subtracted->GetXaxis()->SetTitleOffset(1.2);
    h_piMom_res_subtracted->GetYaxis()->SetTitleOffset(1.45);

    // Add a legend
    TLegend* legend = new TLegend(0.62, 0.75, 0.95, 0.95);
    legend->AddEntry(h_piMom_unscattered_res, "#splitline{Events in unscattered-}{enhanced subset}", "l");
    legend->AddEntry(h_piMom_res_subtracted, "#splitline{Additional events only}{in general selection}", "l");
    legend->Draw();

    // Remove stats box
    h_piMom_res_subtracted->SetStats(0);

    // Add LaTeX text to the plot
    TLatex latex1;
    latex1.SetTextSize(0.04);
    latex1.DrawLatexNDC(0.68, 0.96, "MicroBooNE Simulation");

    // Save the canvas
    c2->SaveAs("plots/MomentumResolution_piMom_comparison.pdf");

    // Delete the TCanvas
    delete c2;

    // Create a custom plot showing both h_muMom_contained_res and h_muMom_uncontained_res as line plots without fill
    TCanvas* c3 = new TCanvas("c3_muMom_res_comparison", "c3_muMom_res_comparison", 800, 500);

    // Adjust margins to ensure everything fits on screen
    c3->SetLeftMargin(0.17);
    c3->SetRightMargin(0.05);
    c3->SetTopMargin(0.05);
    c3->SetBottomMargin(0.15);


    // Set properties for the uncontained histogram
    h_muMom_uncontained_res->SetLineColor(kRed);
    h_muMom_uncontained_res->SetLineWidth(4);
    h_muMom_uncontained_res->SetLineStyle(2);
    h_muMom_uncontained_res->SetFillStyle(0); // No fill

    // Set properties for the contained histogram
    h_muMom_contained_res->SetLineColor(kCyan);
    h_muMom_contained_res->SetLineWidth(4);
    h_muMom_contained_res->SetLineStyle(1);
    h_muMom_contained_res->SetFillStyle(0); // No fill

    // Rescale y axis to fit both plots
    max_y = std::max(h_muMom_contained_res->GetMaximum(), h_muMom_uncontained_res->GetMaximum());
    h_muMom_contained_res->SetMaximum(max_y * 1.1);
    h_muMom_uncontained_res->SetMaximum(max_y * 1.1);

    // Draw the histograms
    h_muMom_uncontained_res->Draw("hist");
    h_muMom_contained_res->Draw("hist  same");

    // Set axis titles and properties
    h_muMom_uncontained_res->SetTitle("");
    h_muMom_uncontained_res->GetXaxis()->SetTitle("(p_{#mu}^{reco} - p_{#mu}^{truth})/p_{#mu}^{truth}");
    h_muMom_uncontained_res->GetYaxis()->SetTitle("#splitline{Selected Simulated Signal Events}{with Correct Muon Candidate}");
    h_muMom_uncontained_res->GetXaxis()->SetLabelSize(0.05);
    h_muMom_uncontained_res->GetYaxis()->SetLabelSize(0.05);
    h_muMom_uncontained_res->GetXaxis()->SetTitleSize(0.05);
    h_muMom_uncontained_res->GetYaxis()->SetTitleSize(0.05);
    h_muMom_uncontained_res->GetXaxis()->SetLabelOffset(0.012);
    h_muMom_uncontained_res->GetYaxis()->SetLabelOffset(0.012);
    h_muMom_uncontained_res->GetXaxis()->CenterTitle();
    h_muMom_uncontained_res->GetXaxis()->SetTitleOffset(1.2);
    h_muMom_uncontained_res->GetYaxis()->SetTitleOffset(1.45);

    // Add a legend
    TLegend* legend_mu = new TLegend(0.62, 0.75, 0.95, 0.95);
    legend_mu->AddEntry(h_muMom_contained_res, "#splitline{Events in contained}{muon subset}", "l");
    legend_mu->AddEntry(h_muMom_uncontained_res, "#splitline{Additional events only}{in general selection}", "l");
    legend_mu->Draw();

    // Remove stats box
    h_muMom_uncontained_res->SetStats(0);

    // Add LaTeX text to the plot
    TLatex latex2;
    latex2.SetTextSize(0.04);
    latex2.DrawLatexNDC(0.68, 0.96, "MicroBooNE Simulation");

    // Save the canvas
    c3->SaveAs("plots/MomentumResolution_muMom_comparison.pdf");

    // Delete the TCanvas
    delete c3;
}