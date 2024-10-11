void testRecoVsTruth() {
    // Open the file
    TFile* file = TFile::Open("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/100Percent_10/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root");

    // Get the TTree
    TTree* tree = (TTree*)file->Get("stv_tree");

    // Create histograms
    TH1F* h1 = new TH1F("h1", "Muon CosTheta", 10, -1, 1);
    TH1F* h2 = new TH1F("h2", "Muon CosTheta", 10, -1, 1);

    // Define the weight
    std::string weight = "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)";

    // Fill the histograms with weights
    tree->Project("h1", "cc1pi_reco_muonCosTheta", ("(cc1pi_signal && cc1pi_selected_generic)*" + weight).c_str());
    tree->Project("h2", "cc1pi_truth_muonCosTheta", ("(cc1pi_signal && cc1pi_selected_generic)*" + weight).c_str());

    // Draw the histograms
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    h1->Draw();
    c1->SaveAs("plots/testRecoVsTruth_muonCosTheta_reco.png");

    TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
    h2->Draw();
    c2->SaveAs("plots/testRecoVsTruth_muonCosTheta_truth.png");

    // Clean up
    delete c1;
    delete h1;
    delete c2;
    delete h2;
    delete tree;
    delete file;
}