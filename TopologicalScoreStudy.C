#include <map>
#include <string>

void TopologicalScoreStudy() 
{
    const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/1March24/";
    const std::string rootPath2 = "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/";

    // List of files; tuple: file type, run, path, fileWeight
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        // std::make_tuple("beam on", "1",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi.root", 1.f),
        std::make_tuple("beam on", "2",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "3",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "4bcd", rootPath + "bnb_beam_on_peleeTuple_uboone_run4bcd_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "5",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root", 1.f),

        // std::make_tuple("beam off", "1",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi.root", 0.56421),
        std::make_tuple("beam off", "2",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 0.40202),
        // std::make_tuple("beam off", "3",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 0.30657),
        // std::make_tuple("beam off", "4bcd", rootPath + "bnb_beam_off_peleeTuple_uboone_run4bcd_ubcc1pi.root", 0.30175),
        // std::make_tuple("beam off", "5",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root", 0.32807),

        // std::make_tuple("nu mc", "1",  rootPath2 + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
        std::make_tuple("nu mc", "2",  rootPath2 + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root", 0.25750*2.0),
        // std::make_tuple("nu mc", "3",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root", 0.20113*2.0),
        // std::make_tuple("nu mc", "4bcd", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root", 0.13074*2.0),
        // std::make_tuple("nu mc", "5",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196*2.0),

        // std::make_tuple("dirt mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root", 0.52806),
        std::make_tuple("dirt mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root", 0.27521),
        // std::make_tuple("dirt mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root", 0.80892),
        // std::make_tuple("dirt mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi.root", 0.39701),
        // std::make_tuple("dirt mc", "5",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_ubcc1pi.root", 0.41280),
    };

    // 2D TH hist
    // TH2D* mc2DHist = new TH2D("mc2DHist", "mc2DHist", 15, -3.1415, 3.1415, 10, 0, 1);
    // TH2D* data2DHist = new TH2D("data2DHist", "data2DHist", 15, -3.1415, 3.1415, 10, 0, 1);

    // Define the bin edges for the Y-axis
    // double yBinEdges[] = {0.2, 0.67, 1.0};
    double yBinEdges[] = {0.0, 0.2, 1.0};
    
    // 2D TH hist with variable bin widths for the Y-axis
    TH2D* mc2DHist = new TH2D("mc2DHist", "MC 2D Histogram;Muon Phi;Topological Score", 15, -3.1415, 3.1415, 2, yBinEdges);
    TH2D* data2DHist = new TH2D("data2DHist", "Data 2D Histogram;Muon Phi;Topological Score", 15, -3.1415, 3.1415, 2, yBinEdges);

    // Loop over the files
    for (const auto &[fileType, run, path, fileWeight] : files)
    {
        std::cout << "Processing file " << path <<std::endl;
        // const auto weight = weights.at(f);
        // std::cout << "Processing file " << path << " with weight " << weight << std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(path.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");

        // const std::string weightString = "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)*" + std::to_string(fileWeight);

        // Define variables to hold branch values
        Float_t spline_weight, tuned_cv_weight, muonPhi, topologicalScore;
        Bool_t isTrueCC1pi, isTrueGoldenPion, passedMuonNotInGap;
        Int_t nUncontained;

        // Set branch addresses
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);
        tree->SetBranchAddress("passed_muonNotInGap", &passedMuonNotInGap);
        // tree->SetBranchAddress("passed_2NonProtons", &passedMuonNotInGap);

        tree->SetBranchAddress("event_cutValue_muonPhi", &muonPhi);
        tree->SetBranchAddress("event_cutValue_topologicalScore", &topologicalScore);

        tree->SetBranchAddress("event_cutValue_nUncontained", &nUncontained);


        // std::cout << "\033[1;31mWARNING - Only processing 10% of events!!!\033[0m" << std::endl;
        const auto nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++)
        {
            tree->GetEntry(i);
            if(!passedMuonNotInGap) continue;
            if(nUncontained != 0) continue;

            if(i%(nEntries/100)==0) std::cout<<"\r"<<(100*i)/nEntries<<"%"<<std::flush;

            double weight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;
            weight *= fileWeight;

            if(fileType == "beam on") data2DHist->Fill(muonPhi, topologicalScore, weight);
            else mc2DHist->Fill(muonPhi, topologicalScore, weight);

        } // End of loop over events
    } // End of loop over files

    // Define the number of colors and the stops for the gradient
    const Int_t NRGBs = 3;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = {0.00, 0.50, 1.00}; // Positions of the colors along the gradient
    Double_t red[NRGBs]   = {1.00, 1.00, 0.00}; // Red at 0.00, then transitioning to white (1.00) and finally to 0.00 at the end
    Double_t green[NRGBs] = {0.00, 1.00, 0.00}; // Starts at 0.00, white in the middle (1.00), back to 0.00
    Double_t blue[NRGBs]  = {0.00, 1.00, 1.00}; // Starts at 0.00, transitions to white (1.00), and ends at blue (1.00)

    // Create the custom color palette
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // Draw and save the MC 2D histogram
    TCanvas *c1 = new TCanvas("c1", "MC 2D Histogram", 800, 600);
    mc2DHist->SetStats(false); // Disable the stats box
    mc2DHist->Draw("COLZ"); // Draw with color palette
    mc2DHist->GetXaxis()->SetTitle("Muon Phi");
    mc2DHist->GetYaxis()->SetTitle("Topological Score");
    mc2DHist->GetZaxis()->SetTitle("Counts");
    mc2DHist->GetZaxis()->SetRangeUser(0.5, 1.5);
    c1->SaveAs("plots/mc2DHist_muonPhiVsToplogicalScore.pdf");
    
    // Draw and save the Data 2D histogram
    TCanvas *c2 = new TCanvas("c2", "Data 2D Histogram", 800, 600);
    data2DHist->SetStats(false); // Disable the stats box
    data2DHist->Draw("COLZ"); // Draw with color palette
    data2DHist->GetXaxis()->SetTitle("Muon Phi");
    data2DHist->GetYaxis()->SetTitle("Topological Score");
    data2DHist->GetZaxis()->SetTitle("Counts");
    data2DHist->GetZaxis()->SetRangeUser(0.5, 1.5);
    c2->SaveAs("plots/data2DHist_muonPhiVsToplogicalScore.pdf");

    // Clone the MC histogram to use as a template for the ratio histogram
    TH2D* ratioHist = (TH2D*)mc2DHist->Clone("ratioHist");
    ratioHist->SetTitle("Ratio Histogram;Muon Phi;Topological Score;Data/MC Ratio");

    // Calculate the ratio in each bin
    for (int xBin = 1; xBin <= mc2DHist->GetNbinsX(); ++xBin) {
        for (int yBin = 1; yBin <= mc2DHist->GetNbinsY(); ++yBin) {
            double mcValue = mc2DHist->GetBinContent(xBin, yBin);
            double dataValue = data2DHist->GetBinContent(xBin, yBin);
            double ratio = 0;
            if (mcValue != 0) { // Avoid division by zero
                ratio = dataValue / mcValue;
            }
            ratioHist->SetBinContent(xBin, yBin, ratio);
        }
    }

    // Draw and save the ratio histogram
    TCanvas *c3 = new TCanvas("c3", "Ratio Histogram", 800, 600);
    ratioHist->SetStats(false); // Disable the stats box
    ratioHist->Draw("COLZ"); // Draw with color palette
    ratioHist->GetXaxis()->SetTitle("Muon Phi");
    ratioHist->GetYaxis()->SetTitle("Topological Score");
    ratioHist->GetZaxis()->SetTitle("Data/MC Ratio");
    ratioHist->GetZaxis()->SetRangeUser(0.5, 1.5);
    c3->SaveAs("plots/ratioHist_muonPhiVsToplogicalScore.pdf");



    // Clean up
    delete c1;
    delete c2;
    delete c3;
}