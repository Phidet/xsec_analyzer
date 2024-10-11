void test() {
  TFile *file = new TFile("/uboone/data/users/jdetje/ubcc1piOutput/3Percent_2/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root");
  TTree *tree = (TTree*)file->Get("stv_tree");

  TH1D *histogram1 = new TH1D("histogram1", "Muon Momentum Histogram", 10, 0, 1); // Adjust binning as needed
  TH1D *histogram2 = new TH1D("histogram2", "Muon Momentum Histogram", 10, 0, 1); // Adjust binning as needed

//   tree->Draw("cc1pi_reco_muonMomentum>>histogram", "cc1pi_selected_generic");
//   tree->Draw("cc1pi_reco_pionMomentum>>histogram", "cc1pi_selected_generic");
    // tree->Draw("cc1pi_reco_muonCosTheta>>histogram", "cc1pi_selected_generic");

    // tree->Draw("cc1pi_reco_muonCosTheta>>histogram", "");
    // tree->Draw("cc1pi_truth_muonCosTheta>>histogram", "cc1pi_signal && cc1pi_selected_generic");

    // tree->Draw("cc1pi_truth_muonCosTheta>>histogram1", "cc1pi_signal", "same");
    // tree->Draw("cc1pi_truth_muonCosTheta>>histogram2", "cc1pi_selected_generic && cc1pi_signal");

    // tree->Draw("cc1pi_signal>>histogram1", "cc1pi_signal");
    // tree->Draw("cc1pi_signal>>histogram2", "cc1pi_selected_generic && cc1pi_signal");

    tree->Draw("event_cutValue_goldenPionBDT>>histogram1", "");
    // tree->Draw("cc1pi_reco_muonPhi>>histogram2", "cc1pi_selected_generic");

    //   tree->Draw("cc1pi_truth_muonMomentum>>histogram", "cc1pi_signal");

    histogram1->SetTitle("Muon Momentum");
    histogram1->GetXaxis()->SetTitle("Muon Momentum");
    histogram1->GetYaxis()->SetTitle("Number of Events");
    histogram1->Draw();

    // histogram2->SetTitle("Muon Momentum");
    // histogram2->GetXaxis()->SetTitle("Muon Momentum");
    // histogram2->GetYaxis()->SetTitle("Number of Events");
    // histogram2->SetLineStyle(2); // Set dotted line style
    // histogram2->Draw("same");
}
