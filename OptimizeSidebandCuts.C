#include <TChain.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>


void FillHistograms(TTree* tree, TH1F* cc0pi_hist, TH1F* cc1pi_hist, const std::string& x_axis_variable1, const std::string& x_axis_variable2, const float runWeight, const float upperProtonBDTLimit = 0.1, const float lowerMuonBDTLimit = -0.4)
{
    std::cout<<"DEBUG FillHistograms: upperProtonBDTLimit = "<<upperProtonBDTLimit<<", lowerMuonBDTLimit = "<<lowerMuonBDTLimit<<" runWeight = "<<runWeight<<" x_axis_variable1 = "<<x_axis_variable1<<" x_axis_variable2 = "<<x_axis_variable2<<std::endl;
    // // Get the names of the histograms
    // std::string cc0pi_hist_name = cc0pi_hist->GetName();
    // std::string cc1pi_hist_name = cc1pi_hist->GetName();

    TH1F* cc0pi_hist_tmp = (TH1F*)cc0pi_hist->Clone(("cc0pi_hist_tmp_name"+x_axis_variable1+x_axis_variable2).c_str());
    TH1F* cc1pi_hist_tmp = (TH1F*)cc1pi_hist->Clone(("cc1pi_hist_tmp_name"+x_axis_variable1+x_axis_variable2).c_str());

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

    // Fill the histograms with weights
    tree->Project(cc0pi_hist_tmp_name.c_str(), x_axis_variable1.c_str(), ("(true_cc0pi && cc0pi_selected_generic && cc0pi_leadingProton_protonBDTScore < " + std::to_string(upperProtonBDTLimit) + " && cc0pi_leadingProton_muonBDTScore > " + std::to_string(lowerMuonBDTLimit) + ")*" + weight).c_str());
    tree->Project(cc1pi_hist_tmp_name.c_str(), x_axis_variable2.c_str(), ("(true_cc0pi && cc1pi_selected_generic)*" + weight).c_str());

    // Add the tmp histograms to the main histograms
    cc0pi_hist->Add(cc0pi_hist_tmp);
    cc1pi_hist->Add(cc1pi_hist_tmp);


    delete cc0pi_hist_tmp, cc1pi_hist_tmp;
}

void MakePlot(TH1F* cc0pi_hist, TH1F* cc1pi_hist, const std::string& name, const bool normalize = true)
{

    // Normalize the histograms
    if(normalize)
    {
        const double epsilon = 1e-8;
        double cc0piIntegral = cc0pi_hist->Integral();
        double cc1piIntegral = cc1pi_hist->Integral();

        if(cc0piIntegral < epsilon || cc1piIntegral < epsilon) {
            throw std::runtime_error("Cannot normalize histograms as one of the integrals is zero");
        }

        cc0pi_hist->Scale(1.0 / cc0piIntegral);
        cc1pi_hist->Scale(1.0 / cc1piIntegral);
    }

    // Set the color and line thickness of the histograms
    cc0pi_hist->SetLineColor(kGreen);
    cc0pi_hist->SetLineWidth(3);
    cc1pi_hist->SetLineColor(kBlack);
    cc1pi_hist->SetLineWidth(3);

    // Get the maximum y value from both histograms
    double maxY = std::max(cc0pi_hist->GetMaximum(), cc1pi_hist->GetMaximum());

    // Set the y-axis range to be slightly more than the maximum y value
    cc1pi_hist->GetYaxis()->SetRangeUser(0, maxY * 1.3);

    // Remove the stats box
    cc1pi_hist->SetStats(0);
    cc0pi_hist->SetStats(0);

    // Draw the histogram
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    cc1pi_hist->Draw("E hist");
    cc0pi_hist->Draw("E hist same");

    // Add a legend
    // TLegend* legend = new TLegend(0.2, 0.1);
    // legend->AddEntry(cc0pi_hist, "CC0pi", "l");
    // legend->AddEntry(cc1pi_hist, "CC1pi", "l");
    // legend->Draw();

    c1->SaveAs(("plots/optimizeSidebandCuts_true_" + name + ".pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void OptimizeSidebandCuts() 
{
    // List of files
    const std::vector<std::string> files = {
        "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/100Percent_10/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root",
        // "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root",
        // "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root",
        // "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/overlay_peleeTuple_uboone_v08_00_00_73_run4a_nu_ubcc1pi.root"
    };

    const std::vector<float> overlayPOT {1.28278e+21, 1.01592e+21, 1.31355e+21, 2.48772e+20, 7.64399e+20, 4.64842e+20, 8.66958e+20, 9.70631e+20}; // run 1, 2, 3, 4a, 4b, 4c, 4d, 5
    const std::vector<float> dataBNBTor860WCut {1.669e+20, 2.616e+20, 2.562e+20, 3.607e+19, 1.39e+20, 8.586e+19, 4.919e+19, 1.23e+20};
    std::vector<float> ratioPOT;
    for (size_t i = 0; i < overlayPOT.size(); i++)
    {
        ratioPOT.push_back(dataBNBTor860WCut[i] / overlayPOT[i]);
    }
    const std::vector<float> weights = {ratioPOT.at(0)};//, ratioPOT.at(1), ratioPOT.at(2)};

    const std::map<std::string, std::vector<double>> binEdges{
        {"muonMomentum", {10, 0, 1.5}},
        {"muonPhi", {15, -3.141592654, 3.141592654}},
        {"muonCosTheta", {15, -1, 1}},
        {"nProtons", {4, 0, 4}},
        {"pionMomentum", {4, 0.3, 1.5}},
        {"pionPhi", {15, -3.141592654, 3.141592654}},
        {"pionCosTheta", {15, -1, 1}},
        // {"muonPionAngle", {0, 0.49, 0.93, 1.26, 1.57, 1.88, 2.21, 2.65}},
    };

    std::vector<std::string> variables1 = {"cc0pi_truth_muonMomentum", "cc0pi_truth_muonPhi", "cc0pi_truth_muonCosTheta", "cc0pi_truth_nProtons", "cc0pi_backtracked_protonMomentum", "cc0pi_backtracked_protonPhi", "cc0pi_backtracked_protonCosTheta"};
    std::vector<std::string> variables2 = {"cc0pi_truth_muonMomentum", "cc0pi_truth_muonPhi", "cc0pi_truth_muonCosTheta", "cc0pi_truth_nProtons", "cc1pi_backtracked_protonMomentum", "cc1pi_backtracked_protonPhi", "cc1pi_backtracked_protonCosTheta"}; // Backtracked variables are different
    
    // std::vector<std::string> variables1 = {"cc0pi_truth_muonMomentum", "cc0pi_truth_muonPhi", "cc0pi_truth_muonCosTheta", "cc0pi_truth_nProtons", "cc0pi_truth_protonMomentum", "cc0pi_truth_protonPhi", "cc0pi_truth_protonCosTheta"};
    // std::vector<std::string> variables2 = variables1;
    std::vector<std::string> names = {"muonMomentum", "muonPhi", "muonCosTheta", "nProtons", "pionMomentum", "pionPhi", "pionCosTheta"};
    if(variables1.size() != variables2.size()) {
        throw std::runtime_error("variables1 and variables2 must have the same size");
    }

    std::map<std::string, TH1F*> cc0pi_histograms;
    std::map<std::string, TH1F*> cc1pi_histograms;

    for (size_t i = 0; i < names.size(); ++i) {
        cc0pi_histograms[names.at(i)] = new TH1F(("cc0pi_" + names.at(i)).c_str(), (names.at(i) + " of true CC0pi events selected by CC0pi selection").c_str(), binEdges.at(names.at(i)).at(0), binEdges.at(names.at(i)).at(1), binEdges.at(names.at(i)).at(2));
        cc1pi_histograms[names.at(i)] = new TH1F(("cc1pi_" + names.at(i)).c_str(), (names.at(i) + " of true CC0pi events selected by CC1pi selection").c_str(), binEdges.at(names.at(i)).at(0), binEdges.at(names.at(i)).at(1), binEdges.at(names.at(i)).at(2));
        cc0pi_histograms.at(names.at(i))->Sumw2();
        cc1pi_histograms.at(names.at(i))->Sumw2();
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

        // auto muonCutValue = [](const int l) { return  -0.65f+l*0.05f; };
        // auto protonCutValue = [](const int k) { return 0.025f+k*0.0125f; };

        auto muonCutValue = [](const int l) { return  -0.65f; };
        auto protonCutValue = [](const int k) { return 0.025f; };

        // Fill and the histograms
        // for (size_t i = 0; i < variables1.size(); ++i) {
        const size_t v = 4;
        std::cout<<"----- Optimizing for "<<names.at(v)<<" -----"<<std::endl;

        float lowestScore = std::numeric_limits<float>::max();
        int lowestScore_l = 0;
        int lowestScore_k = 0;
        for(int k = 0; k < 2; k++)
        {
            for(int l = 0; l < 2; l++)
            {
                std::cout<<"protonCutValue = "<<protonCutValue(k)<<", muonCutValue = "<<muonCutValue(l)<<": "<<std::flush;
                const auto cc0pi_hist = new TH1F(("cc0pi_" + std::to_string(k) + "_" + std::to_string(l) + "_" + names.at(v)).c_str(), (names.at(v) + " of true CC0pi events selected by CC0pi selection").c_str(), binEdges.at(names.at(v)).at(0), binEdges.at(names.at(v)).at(1), binEdges.at(names.at(v)).at(2));
                const auto cc1pi_hist = new TH1F(("cc1pi_" + std::to_string(k) + "_" + std::to_string(l) + "_" + names.at(v)).c_str(), (names.at(v) + " of true CC0pi events selected by CC1pi selection").c_str(), binEdges.at(names.at(v)).at(0), binEdges.at(names.at(v)).at(1), binEdges.at(names.at(v)).at(2));
                FillHistograms(tree, cc0pi_hist, cc1pi_hist, variables1.at(v), variables2.at(v), weight, protonCutValue(k), muonCutValue(l));

                // Print out the integral values
                std::cout<<"----> cc0pi: "<<cc0pi_hist->Integral()<<" cc1pi: "<<cc1pi_hist->Integral()<<std::endl;

                // Normalise the histograms first by dividing through their integrals
                cc0pi_hist->Scale(1.0 / cc0pi_hist->Integral());
                cc1pi_hist->Scale(1.0 / cc1pi_hist->Integral());

                // Calculate a difference metric
                // const auto diff = cc0pi_hist - cc1pi_hist;
                // const auto squaredDiff = diff * diff;
                // float score = squaredDiff.Integral();
                const auto score = cc0pi_hist->Chi2Test(cc1pi_hist, "WW CHI2/NDF");
                if(score<lowestScore)
                {
                    lowestScore = score;
                    lowestScore_l = l;
                    lowestScore_k = k;
                }

                std::cout<<score<<" "<<std::flush;


                // TH1F* cc0pi_hist = cc0pi_hist;
                // TH1F* cc1pi_hist = cc1pi_hist;

                int nBins = cc0pi_hist->GetNbinsX();
                for (int bin = 1; bin <= nBins; ++bin) {
                    float cc0pi_binContent = cc0pi_hist->GetBinContent(bin);
                    float cc1pi_binContent = cc1pi_hist->GetBinContent(bin);

                    std::cout << "Bin " << bin << ": "
                            << "cc0pi = " << cc0pi_binContent << ", "
                            << "cc1pi = " << cc1pi_binContent << std::endl;
                }


                cc0pi_hist->Reset();
                cc1pi_hist->Reset();
                delete cc0pi_hist, cc1pi_hist;
            }
            std::cout<<std::endl;
        }

        std::cout<<"----- Result -----"<<std::endl;
        std::cout<<"Lowest score: "<<lowestScore<<" with protonCutValue = "<<protonCutValue(lowestScore_k)<<" and muonCutValue = "<<muonCutValue(lowestScore_l)<<std::endl;
        
        for (size_t i = 0; i < variables1.size(); ++i) {
            std::cout<<"Filling histograms for "<<names.at(i)<<std::endl;
            FillHistograms(tree, cc0pi_histograms.at(names.at(v)), cc1pi_histograms.at(names.at(v)), variables1.at(v), variables2.at(v), weight, protonCutValue(lowestScore_k), muonCutValue(lowestScore_l));
        }

        // }
        tFile->Close();
    }

    // Make the plots
    const auto normalize = true;
    for (size_t i = 0; i < names.size(); ++i) {
        MakePlot(cc0pi_histograms.at(names.at(i)), cc1pi_histograms.at(names.at(i)), names.at(i), normalize);
    }
}