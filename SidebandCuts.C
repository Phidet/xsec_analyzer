#include <TChain.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>


void FillHistograms(TTree* tree, TH1F* cc0pi_hist, TH1F* cc1pi_hist, const std::string& x_axis_variable1, const std::string& x_axis_variable2, const float runWeight, bool additionalCuts = true)
{
    // // Get the names of the histograms
    // std::string cc0pi_hist_name = cc0pi_hist->GetName();
    // std::string cc1pi_hist_name = cc1pi_hist->GetName();

    TH1F* cc0pi_hist_tmp = (TH1F*)cc0pi_hist->Clone((cc0pi_hist->GetName()+x_axis_variable1+x_axis_variable2).c_str());
    TH1F* cc1pi_hist_tmp = (TH1F*)cc1pi_hist->Clone((cc1pi_hist->GetName()+x_axis_variable1+x_axis_variable2).c_str());

    cc0pi_hist_tmp->Reset();
    cc1pi_hist_tmp->Reset();

    std::string cc0pi_hist_tmp_name = cc0pi_hist_tmp->GetName();
    std::string cc1pi_hist_tmp_name = cc1pi_hist_tmp->GetName();
    
    // Define the weight condition
    // std::string weight = "(weight_splines_general_Spline*weight_TunedCentralValue_UBGenie >= 0 && weight_splines_general_Spline*weight_TunedCentralValue_UBGenie <= 30) ? weight_splines_general_Spline*weight_TunedCentralValue_UBGenie : 1";
    // std::string weight = "(weight_splines_general_Spline*weight_TunedCentralValue_UBGenie >= 0 && weight_splines_general_Spline*weight_TunedCentralValue_UBGenie <= 30 && std::isfinite(weight_splines_general_Spline*weight_TunedCentralValue_UBGenie)) ? weight_splines_general_Spline*weight_TunedCentralValue_UBGenie : 1";
    std::string weight = "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)*" + std::to_string(runWeight);
    // std::string weight = "(std::isfinite(spline_weight) ? spline_weight : 1)*(std::isfinite(tuned_cv_weight) ? tuned_cv_weight : 1) *" + std::to_string(runWeight);
    // std::string weight = "1";

    std::string additionalCutString = additionalCuts ? " && cc0pi_leadingProton_protonBDTScore < 0.1 && cc0pi_leadingProton_muonBDTScore > -0.4" : "";

    // Fill the histograms with weights
    // tree->Project(cc0pi_hist_tmp_name.c_str(), x_axis_variable1.c_str(), ("(true_cc0pi && cc0pi_selected_generic && cc0pi_leadingProton_protonBDTScore < 0.4 && cc0pi_leadingProton_muonBDTScore > -0.5)*" + weight).c_str());
    tree->Project(cc0pi_hist_tmp_name.c_str(), x_axis_variable1.c_str(), ("(true_cc0pi && cc0pi_selected_generic && cc0pi_reco_protonMomentum > 0.1" + additionalCutString + ")*" + weight).c_str());
    // tree->Project(cc0pi_hist_tmp_name.c_str(), x_axis_variable1.c_str(), ("(true_cc0pi && cc0pi_selected_generic)*" + weight).c_str());
    tree->Project(cc1pi_hist_tmp_name.c_str(), x_axis_variable2.c_str(), ("(true_cc0pi && cc1pi_selected_generic && cc1pi_reco_pionMomentum > 0.1)*" + weight).c_str());

    // chain->Project(cc0pi_hist_name.c_str(), x_axis_variable1.c_str(), ("(true_cc0pi && cc0pi_selected_generic && cc0pi_leadingProton_protonBDTScore < 0.10 && cc0pi_leadingProton_muonBDTScore > -0.4 && cc0pi_backtracked_protonMomentum > 0)*" + weight).c_str());
    // chain->Project(cc1pi_hist_name.c_str(), x_axis_variable2.c_str(), ("(true_cc0pi && cc1pi_selected_generic && cc1pi_backtracked_protonMomentum > 0)*" + weight).c_str());

    // Add the tmp histograms to the main histograms
    cc0pi_hist->Add(cc0pi_hist_tmp);
    cc1pi_hist->Add(cc1pi_hist_tmp);

    delete cc0pi_hist_tmp, cc1pi_hist_tmp;
}

void MakePlot(TH1F* cc0pi_hist, TH1F* cc1pi_hist, const std::string& name, const bool normalize = true, bool additionalCuts = true)
{

    // Normalize the histograms
    if(normalize)
    {
        const double epsilon = 1e-8;
        double cc0piIntegral = cc0pi_hist->Integral();
        double cc1piIntegral = cc1pi_hist->Integral();

        if(cc0piIntegral < epsilon || cc1piIntegral < epsilon) {
            throw std::runtime_error("Cannot normalize histograms as one of the integrals is zero. [cc0piIntegral: " + std::to_string(cc0piIntegral) + ", cc1piIntegral: " + std::to_string(cc1piIntegral) + "]");
        }

        cc0pi_hist->Scale(1.0 / cc0piIntegral);
        cc1pi_hist->Scale(1.0 / cc1piIntegral);
    }

    // Set the color and line thickness of the histograms
    if(additionalCuts)
    {
        cc0pi_hist->SetLineColor(kGreen + 2);
    }
    else
    {
        cc0pi_hist->SetLineColor(kBlue);
    }
    cc0pi_hist->SetLineWidth(3);
    cc1pi_hist->SetLineColor(kBlack);
    cc1pi_hist->SetFillColor(kGray);
    cc1pi_hist->SetFillStyle(3001); // light grey, semi-transparent
    cc1pi_hist->SetLineWidth(3);

    cc1pi_hist->GetXaxis()->SetTitleOffset(cc1pi_hist->GetXaxis()->GetTitleOffset() + 0.1);

    // Remove the stats box
    cc1pi_hist->SetStats(0);
    cc0pi_hist->SetStats(0);

    std::cout << "Calculating maximum y value with error..." << std::endl;
    double maxYWithError = 0.0;
    for (int i = 1; i <= cc1pi_hist->GetNbinsX(); ++i) {
        double val = cc1pi_hist->GetBinContent(i);
        double err = cc1pi_hist->GetBinError(i);
        if (val + err > maxYWithError)
            maxYWithError = val + err;
    }
    for (int i = 1; i <= cc0pi_hist->GetNbinsX(); ++i) {
        double val = cc0pi_hist->GetBinContent(i);
        double err = cc0pi_hist->GetBinError(i);
        if (val + err > maxYWithError)
            maxYWithError = val + err;
    }
    cc1pi_hist->GetYaxis()->SetRangeUser(0, maxYWithError * 1.3);

    // Draw the histogram
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    cc1pi_hist->Draw("E hist");
    cc0pi_hist->Draw("E hist same");

    // Add a legend
    TLegend* legend = new TLegend(0.15, 0.75, 0.70, 0.89);
    legend->AddEntry(cc1pi_hist, "True CC0#pi Xp, X#geq 1 in CC1#pi^{#pm} Selection", "l");
    legend->AddEntry(cc0pi_hist, "True CC0#pi Xp, X#geq 1 in Sideband Selection", "l");
    legend->Draw();

    c1->SaveAs(("plots/sidebandCuts_true_" + name + ".pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void SidebandCuts() 
{
    // List of files
    const std::vector<std::string> files = {
        // "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root",
        // "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root",
        // "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root",
        // "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/overlay_peleeTuple_uboone_v08_00_00_73_run4a_nu_ubcc1pi.root"
        // "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root"

        "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root",
        "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root",
        "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root",
        "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root",
        "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root",
    };

    // const std::vector<float> overlayPOT {1.28278e+21, 1.01592e+21, 1.31355e+21, 2.48772e+20, 7.64399e+20, 4.64842e+20, 8.66958e+20, 9.70631e+20}; // run 1, 2, 3, 4a, 4b, 4c, 4d, 5
    // const std::vector<float> dataBNBTor860WCut {1.669e+20, 2.616e+20, 2.562e+20, 3.607e+19, 1.39e+20, 8.586e+19, 4.919e+19, 1.23e+20};
    // std::vector<float> ratioPOT;
    // for (size_t i = 0; i < overlayPOT.size(); i++)
    // {
    //     ratioPOT.push_back(dataBNBTor860WCut[i] / overlayPOT[i]);
    // }

    std::vector<float> ratioPOT{
        0.13011*2.0, // run 1
        0.25750*2.0, // run 2
        0.20113*2.0, // run 3
        0.13074*2.0, // run 4bcd
        0.15196*2.0  // run 5
    };

    const std::vector<float> weights = ratioPOT;

    const std::map<std::string, std::vector<double>> binEdges{
        {"muonMomentum", {10, 0, 1.5}},
        {"muonPhi", {15, -3.141592654, 3.141592654}},
        {"muonCosTheta", {15, -1, 1}},
        {"nProtons", {4, 0, 4}},
        {"protonMomentum", {15, 0.2, 2.0}},
        {"protonPhi", {15, -3.141592654, 3.141592654}},
        {"protonCosTheta", {15, -1, 1}},
        {"muonProtonAngle", {15, 0, 3.141592654}},
        // {"muonPionAngle", {0, 0.49, 0.93, 1.26, 1.57, 1.88, 2.21, 2.65}},
    };

    // std::vector<std::string> variables1 = {"cc0pi_truth_muonMomentum", "cc0pi_truth_muonPhi", "cc0pi_truth_muonCosTheta", "cc0pi_truth_nProtons", "cc0pi_backtracked_protonMomentum", "cc0pi_backtracked_protonPhi", "cc0pi_backtracked_protonCosTheta"};
    // std::vector<std::string> variables2 = {"cc0pi_truth_muonMomentum", "cc0pi_truth_muonPhi", "cc0pi_truth_muonCosTheta", "cc0pi_truth_nProtons", "cc1pi_backtracked_protonMomentum", "cc1pi_backtracked_protonPhi", "cc1pi_backtracked_protonCosTheta"}; // Backtracked variables are different
    
    std::vector<std::string> variables1 = {"cc0pi_truth_muonMomentum", "cc0pi_truth_muonPhi", "cc0pi_truth_muonCosTheta", "cc0pi_truth_nProtons", "cc0pi_truth_protonMomentum", "cc0pi_truth_protonPhi", "cc0pi_truth_protonCosTheta", "cc0pi_truth_muonProtonAngle"};
    std::vector<std::string> variables2 = variables1;

    std::vector<std::string> names = {"muonMomentum", "muonPhi", "muonCosTheta", "nProtons", "protonMomentum", "protonPhi", "protonCosTheta", "muonProtonAngle"};
    std::vector<std::string> latexNames = {"True p_{#mu} / GeV", "True #phi_{#mu} / rad", "True cos(#theta_{#mu})", "True Number of protons", "True p_{p} / GeV", "True #phi_{p} / rad", "True cos(#theta_{p})", "True #theta_{#mu p} / rad"};

    if(variables1.size() != variables2.size()) {
        throw std::runtime_error("variables1 and variables2 must have the same size");
    }

    std::map<std::string, TH1F*> cc0pi_histograms;
    std::map<std::string, TH1F*> cc1pi_histograms;
    std::map<std::string, TH1F*> cc0pi_histograms_noAdditionalCuts;
    std::map<std::string, TH1F*> cc1pi_histograms_noAdditionalCuts;

    for (size_t i = 0; i < names.size(); ++i) {
        cc0pi_histograms[names.at(i)] = new TH1F(("cc0pi_" + names.at(i)).c_str(), (";"+ latexNames.at(i) + ";" + "Area normalised event rate").c_str(), binEdges.at(names.at(i)).at(0), binEdges.at(names.at(i)).at(1), binEdges.at(names.at(i)).at(2));
        cc1pi_histograms[names.at(i)] = new TH1F(("cc1pi_" + names.at(i)).c_str(), (";"+ latexNames.at(i) + ";" + "Area normalised event rate").c_str(), binEdges.at(names.at(i)).at(0), binEdges.at(names.at(i)).at(1), binEdges.at(names.at(i)).at(2));
        cc0pi_histograms.at(names.at(i))->Sumw2();
        cc1pi_histograms.at(names.at(i))->Sumw2();

        // Create histograms without additional cuts
        cc0pi_histograms_noAdditionalCuts[names.at(i)] = new TH1F(("cc0pi_noAdditionalCuts_" + names.at(i)).c_str(), (";"+ latexNames.at(i) + ";" + "Area normalised event rate").c_str(), binEdges.at(names.at(i)).at(0), binEdges.at(names.at(i)).at(1), binEdges.at(names.at(i)).at(2));
        cc1pi_histograms_noAdditionalCuts[names.at(i)] = new TH1F(("cc1pi_noAdditionalCuts_" + names.at(i)).c_str(), (";"+ latexNames.at(i) + ";" + "Area normalised event rate").c_str(), binEdges.at(names.at(i)).at(0), binEdges.at(names.at(i)).at(1), binEdges.at(names.at(i)).at(2));
        cc0pi_histograms_noAdditionalCuts.at(names.at(i))->Sumw2();
        cc1pi_histograms_noAdditionalCuts.at(names.at(i))->Sumw2();
    }

    // Loop over the files
    for (unsigned int f=0; f<files.size(); f++)  
    {
        const auto fileName = files.at(f); 
        const auto weight = weights.at(f);
        std::cout << "Processing file " << fileName << " with weight " << weight << std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(fileName.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");

        // Disable all branches
        tree->SetBranchStatus("*", 0);

        // Enable only the branches you need
        tree->SetBranchStatus("cc0pi_*", 1);
        tree->SetBranchStatus("cc1pi_*", 1);
        tree->SetBranchStatus("true_cc0pi", 1);
        tree->SetBranchStatus("spline_weight", 1);
        tree->SetBranchStatus("tuned_cv_weight", 1);

        // Fill and the histograms
        for (size_t i = 0; i < variables1.size(); ++i) {
            FillHistograms(tree, cc0pi_histograms.at(names.at(i)), cc1pi_histograms.at(names.at(i)), variables1.at(i), variables2.at(i), weight);
            FillHistograms(tree, cc0pi_histograms_noAdditionalCuts.at(names.at(i)), cc1pi_histograms_noAdditionalCuts.at(names.at(i)), variables1.at(i), variables2.at(i), weight, false);
            std::cout << "Done filling histogram for " << names.at(i) << std::endl;
        }
        tFile->Close();
    }

    // Make the plots
    std::cout << "Making plots..." << std::endl;
    const auto normalize = true;
    for (size_t i = 0; i < names.size(); ++i) {
        std::cout << "Making plot for " << names.at(i) << std::endl;
        MakePlot(cc0pi_histograms.at(names.at(i)), cc1pi_histograms.at(names.at(i)), names.at(i), normalize);
        std::cout << "Making plot for " << names.at(i) << "_noAdditionalCuts" << std::endl;
        MakePlot(cc0pi_histograms_noAdditionalCuts.at(names.at(i)), cc1pi_histograms_noAdditionalCuts.at(names.at(i)), names.at(i) + "_noAdditionalCuts", normalize, false);
    }
}