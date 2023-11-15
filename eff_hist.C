// Efficiency histogram plotter

void plot_eff_hist( TTree& stv_tree, const std::string& branch,
  const std::string& signal_cuts, const std::string& selection_cuts,
  const std::string& variable_name, const std::string& hist_name,
  const std::string& unit_name, int num_bins, double x_min, double x_max, const std::string additionalTitle = "")
{
  // For computing efficiencies and purities, we need to only use MC events.
  // Unconditionally add this requirement to the cuts defined above.
  std::string signal = signal_cuts + " && is_mc ";
  std::string selection = selection_cuts + " && is_mc ";

  std::string unit_part;
  if ( !unit_name.empty() ) unit_part = " [" + unit_name + ']';
  std::string eff_title = "Efficiency in true " + variable_name + " bins" + additionalTitle + "; "
    + variable_name + unit_part + "; Efficiency";
  TH1D* eff_hist = new TH1D( hist_name.c_str(), eff_title.c_str(), num_bins,
    x_min, x_max );

  for ( int b = 1; b <= num_bins; ++b ) {
    double xlow = eff_hist->GetBinLowEdge( b );
    double xhigh = eff_hist->GetBinLowEdge( b + 1 );

    // Define a new set of cuts requiring the event to fall inside the current
    // bin
    std::stringstream ss;
    ss << branch << " >= " << xlow << " && "
      << branch << " < " << xhigh;
    std::string bin_cut = ss.str();

    // Combine the bin cuts with the signal and selection cuts
    std::string bin_signal = signal + " && " + bin_cut;
    std::string bin_selection = selection + " && " + bin_cut;

    // These are actually integer counts, but since we will use them below to
    // compute ratios, intrinsically cast them to double-precision values for
    // convenience.
    double num_signal = stv_tree.Draw( "", bin_signal.c_str(), "goff" );
    double num_selected_signal = stv_tree.Draw( "",
      (bin_signal + " && " + selection).c_str(), "goff" );

    // Compute the efficiency for the current bin
    double eff = 0.;
    double eff_stat_err = 0.;
    if ( num_signal > 0. && num_selected_signal > 0. ) {
      eff = num_selected_signal / num_signal;
      eff_stat_err = eff * std::sqrt( (1. / num_selected_signal)
        + (1. / num_signal) );
    }

    eff_hist->SetBinContent( b, eff );
    eff_hist->SetBinError( b, eff_stat_err );
  }

  TCanvas* c1 = new TCanvas;
  eff_hist->SetStats( false );
  eff_hist->GetYaxis()->SetRangeUser( 0., 1. );

  c1->SetBottomMargin(0.15);

  eff_hist->SetLineWidth(4);
  eff_hist->SetLineColor(kBlack);
  eff_hist->GetYaxis()->SetRangeUser(0., 1.);
  eff_hist->GetYaxis()->CenterTitle(true);
  eff_hist->GetYaxis()->SetTitleOffset(0.95);
  eff_hist->GetYaxis()->SetTitleSize(0.05);

  eff_hist->GetXaxis()->SetTitleOffset(1.2);
  eff_hist->GetXaxis()->CenterTitle(true);
  eff_hist->GetXaxis()->SetTitleSize(0.05);

  eff_hist->SetStats( false );
  eff_hist->Draw( );

  c1->SaveAs( ( "plots/eff_hist_" + hist_name + ".jpg").c_str() );
  c1->Close();

}

void plot_pur_hist( TTree& stv_tree, const std::string& branch,
  const std::string& signal_cuts, const std::string& selection_cuts,
  const std::string& variable_name, const std::string& hist_name,
  const std::string& unit_name, int num_bins, double x_min, double x_max, const std::string additionalTitle = "")
{
  // For computing efficiencies and purities, we need to only use MC events.
  // Unconditionally add this requirement to the cuts defined above.
  std::string signal = signal_cuts + " && is_mc ";
  std::string selection = selection_cuts + " && is_mc ";

  std::string unit_part;
  if ( !unit_name.empty() ) unit_part = " [" + unit_name + ']';
  std::string eff_title = "Purity in true " + variable_name + " bins" + additionalTitle + "; "
    + variable_name + unit_part + "; Purity";
  TH1D* eff_hist = new TH1D( hist_name.c_str(), eff_title.c_str(), num_bins,
    x_min, x_max );

  for ( int b = 1; b <= num_bins; ++b ) {
    double xlow = eff_hist->GetBinLowEdge( b );
    double xhigh = eff_hist->GetBinLowEdge( b + 1 );

    // Define a new set of cuts requiring the event to fall inside the current
    // bin
    std::stringstream ss;
    ss << branch << " >= " << xlow << " && "
      << branch << " < " << xhigh;
    std::string bin_cut = ss.str();

    // Combine the bin cuts with the signal and selection cuts
    std::string bin_signal = signal + " && " + bin_cut;
    std::string bin_selection = selection + " && " + bin_cut;

    // These are actually integer counts, but since we will use them below to
    // compute ratios, intrinsically cast them to double-precision values for
    // convenience.
    double num_signal = stv_tree.Draw( "", bin_signal.c_str(), "goff" );
    double num_selected = stv_tree.Draw( "", bin_selection.c_str(), "goff" );
    double num_selected_signal = stv_tree.Draw( "",
      (bin_signal + " && " + selection).c_str(), "goff" );

    // Compute the efficiency for the current bin
    double eff = 0.;
    double eff_stat_err = 0.;
    // std::cout<<"bin "<<b<<" - num_selected: "<<num_selected<<" num_selected_signal: "<<num_selected_signal<<std::endl;
    // std::cout<<"DEBUG: bin_signal: "<<bin_signal<<" - bin_selection: "<<bin_selection<<" - selection: "<<selection<<std::endl;
    if ( num_selected > 0. && num_selected_signal > 0. ) {
      eff = num_selected_signal / num_selected;
      eff_stat_err = eff * std::sqrt( (1. / num_selected_signal)
        + (1. / num_selected) );
    }

    eff_hist->SetBinContent( b, eff );
    eff_hist->SetBinError( b, eff_stat_err );
  }

  TCanvas* c1 = new TCanvas;
  eff_hist->SetStats( false );
  eff_hist->GetYaxis()->SetRangeUser( 0., 1. );

  c1->SetBottomMargin(0.15);

  eff_hist->SetLineWidth(4);
  eff_hist->SetLineColor(kBlack);
  eff_hist->GetYaxis()->SetRangeUser(0., 1.);
  eff_hist->GetYaxis()->CenterTitle(true);
  eff_hist->GetYaxis()->SetTitleOffset(0.95);
  eff_hist->GetYaxis()->SetTitleSize(0.05);

  eff_hist->GetXaxis()->SetTitleOffset(1.2);
  eff_hist->GetXaxis()->CenterTitle(true);
  eff_hist->GetXaxis()->SetTitleSize(0.05);

  eff_hist->SetStats( false );
  eff_hist->Draw( );

  c1->SaveAs( ( "plots/pur_hist_" + hist_name + ".jpg").c_str() );
  c1->Close();
}


void eff_hist() {

  TChain stv_ch( "stv_tree" );
  stv_ch.Add( "/uboone/data/users/jdetje/ubcc1piOutput/3Percent_1/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root" );
  stv_ch.Add( "/uboone/data/users/jdetje/ubcc1piOutput/3Percent_1/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root" );
  stv_ch.Add( "/uboone/data/users/jdetje/ubcc1piOutput/3Percent_1/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root" );

  //// Fake data samples (NuWro)
  //stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/stv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run1_NuWro_reco2_reco2.root" );
  //stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/stv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run2_NuWro_reco2_reco2.root" );
  //stv_ch.Add( "/uboone/data/users/gardiner/ntuples-stv/stv-high_stat_prodgenie_bnb_nu_overlay_DetVar_Run3_NuWro_reco2_reco2.root" );

//  plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "cc1pi_signal", "passed_particleTrackScore", "p_{#mu}^{true}", "pmutrue", "GeV", 20, 0.15, 1.5);
//  plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "cc1pi_signal", "passed_startNearVertex", "p_{#mu}^{true}", "pmutrue", "GeV", 20, 0.15, 1.5);


// // // cc1pi_signal vs true_cc1pi
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_particleTrackScore", "p_{#mu}^{true}", "pmutrue_passed_particleTrackScore", "GeV", 20, 0.15, 1.5, " - Passed particleTrackScore");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_particleVertexDistance", "p_{#mu}^{true}", "pmutrue_passed_particleVertexDistance", "GeV", 20, 0.15, 1.5, " - particleVertexDistance");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_particleGeneration", "p_{#mu}^{true}", "pmutrue_passed_particleGeneration", "GeV", 20, 0.15, 1.5, " - particleGeneration");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_particleTrackLength", "p_{#mu}^{true}", "pmutrue_passed_particleTrackLength", "GeV", 20, 0.15, 1.5, " - particleTrackLength");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_particleProtonChi2", "p_{#mu}^{true}", "pmutrue_passed_particleProtonChi2", "GeV", 20, 0.15, 1.5, " - particleProtonChi2");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_particleMuonChi2", "p_{#mu}^{true}", "pmutrue_passed_particleMuonChi2", "GeV", 20, 0.15, 1.5, " - particleMuonChi2");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_particleProtonChi2OverMuonChi2", "p_{#mu}^{true}", "pmutrue_passed_particleProtonChi2OverMuonChi2", "GeV", 20, 0.15, 1.5, " - particleProtonChi2OverMuonChi2");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_pandoraNuPDGIsNumu", "p_{#mu}^{true}", "pmutrue_passed_pandoraNuPDGIsNumu", "GeV", 20, 0.15, 1.5, " - pandoraNuPDGIsNumu");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_daughterVerticesContained", "p_{#mu}^{true}", "pmutrue_passed_daughterVerticesContained", "GeV", 20, 0.15, 1.5, " - daughterVerticesContained");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_nuVertexFiducial", "p_{#mu}^{true}", "pmutrue_passed_nuVertexFiducial", "GeV", 20, 0.15, 1.5, " - nuVertexFiducial");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_topologicalOrFlashMatch", "p_{#mu}^{true}", "pmutrue_passed_topologicalOrFlashMatch", "GeV", 20, 0.15, 1.5, " - topologicalOrFlashMatch");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_min2Tracks", "p_{#mu}^{true}", "pmutrue_passed_min2Tracks", "GeV", 20, 0.15, 1.5, " - min2Tracks");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_max1Uncontained", "p_{#mu}^{true}", "pmutrue_passed_max1Uncontained", "GeV", 20, 0.15, 1.5, " - max1Uncontained");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_2NonProtons", "p_{#mu}^{true}", "pmutrue_passed_2NonProtons", "GeV", 20, 0.15, 1.5, " - 2NonProtons");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_pionHasValiddEdx", "p_{#mu}^{true}", "pmutrue_passed_pionHasValiddEdx", "GeV", 20, 0.15, 1.5, " - pionHasValiddEdx");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_pionNotInGap", "p_{#mu}^{true}", "pmutrue_passed_pionNotInGap", "GeV", 20, 0.15, 1.5, " - pionNotInGap");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_muonNotInGap", "p_{#mu}^{true}", "pmutrue_passed_muonNotInGap", "GeV", 20, 0.15, 1.5, " - muonNotInGap");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_openingAngle", "p_{#mu}^{true}", "pmutrue_passed_openingAngle", "GeV", 20, 0.15, 1.5, " - openingAngle");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_topologicalScore", "p_{#mu}^{true}", "pmutrue_passed_topologicalScore", "GeV", 20, 0.15, 1.5, " - topologicalScore");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_startNearVertex", "p_{#mu}^{true}", "pmutrue_passed_startNearVertex", "GeV", 20, 0.15, 1.5, " - startNearVertex");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonMomentum", "true_cc1pi", "passed_likelyGoldenPion", "p_{#mu}^{true}", "pmutrue_passed_likelyGoldenPion", "GeV", 20, 0.15, 1.5, " - likelyGoldenPion");


// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_particleTrackScore", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_particleTrackScore", "rad", 20, 0, 3.1415, " - particleTrackScore");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_particleVertexDistance", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_particleVertexDistance", "rad", 20, 0, 3.1415, " - particleVertexDistance");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_particleGeneration", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_particleGeneration", "rad", 20, 0, 3.1415, " - particleGeneration");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_particleTrackLength", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_particleTrackLength", "rad", 20, 0, 3.1415, " - particleTrackLength");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_particleProtonChi2", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_particleProtonChi2", "rad", 20, 0, 3.1415, " - particleProtonChi2");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_particleMuonChi2", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_particleMuonChi2", "rad", 20, 0, 3.1415, " - particleMuonChi2");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_particleProtonChi2OverMuonChi2", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_particleProtonChi2OverMuonChi2", "rad", 20, 0, 3.1415, " - particleProtonChi2OverMuonChi2");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_pandoraNuPDGIsNumu", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_pandoraNuPDGIsNumu", "rad", 20, 0, 3.1415, " - pandoraNuPDGIsNumu");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_daughterVerticesContained", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_daughterVerticesContained", "rad", 20, 0, 3.1415, " - daughterVerticesContained");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_nuVertexFiducial", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_nuVertexFiducial", "rad", 20, 0, 3.1415, " - nuVertexFiducial");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_topologicalOrFlashMatch", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_topologicalOrFlashMatch", "rad", 20, 0, 3.1415, " - topologicalOrFlashMatch");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_min2Tracks", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_min2Tracks", "rad", 20, 0, 3.1415, " - min2Tracks");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_max1Uncontained", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_max1Uncontained", "rad", 20, 0, 3.1415, " - max1Uncontained");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_2NonProtons", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_2NonProtons", "rad", 20, 0, 3.1415, " - 2NonProtons");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_pionHasValiddEdx", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_pionHasValiddEdx", "rad", 20, 0, 3.1415, " - pionHasValiddEdx");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_pionNotInGap", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_pionNotInGap", "rad", 20, 0, 3.1415, " - pionNotInGap");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_muonNotInGap", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_muonNotInGap", "rad", 20, 0, 3.1415, " - muonNotInGap");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_openingAngle", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_openingAngle", "rad", 20, 0, 3.1415, " - openingAngle");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_topologicalScore", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_topologicalScore", "rad", 20, 0, 3.1415, " - topologicalScore");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_startNearVertex", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_startNearVertex", "rad", 20, 0, 3.1415, " - startNearVertex");
// plot_eff_hist( stv_ch, "cc1pi_truth_muonPionAngle", "true_cc1pi", "passed_likelyGoldenPion", "#theta_{#mu#pi}^{true}", "openingangletrue_passed_likelyGoldenPion", "rad", 20, 0, 3.1415, " - likelyGoldenPion");


// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_particleTrackScore", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_particleTrackScore", "rad", 20, 0, 3.1415, " - particleTrackScore");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_particleVertexDistance", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_particleVertexDistance", "rad", 20, 0, 3.1415, " - particleVertexDistance");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_particleGeneration", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_particleGeneration", "rad", 20, 0, 3.1415, " - particleGeneration");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_particleTrackLength", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_particleTrackLength", "rad", 20, 0, 3.1415, " - particleTrackLength");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_particleProtonChi2", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_particleProtonChi2", "rad", 20, 0, 3.1415, " - particleProtonChi2");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_particleMuonChi2", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_particleMuonChi2", "rad", 20, 0, 3.1415, " - particleMuonChi2");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_particleProtonChi2OverMuonChi2", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_particleProtonChi2OverMuonChi2", "rad", 20, 0, 3.1415, " - particleProtonChi2OverMuonChi2");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_pandoraNuPDGIsNumu", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_pandoraNuPDGIsNumu", "rad", 20, 0, 3.1415, " - pandoraNuPDGIsNumu");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_daughterVerticesContained", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_daughterVerticesContained", "rad", 20, 0, 3.1415, " - daughterVerticesContained");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_nuVertexFiducial", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_nuVertexFiducial", "rad", 20, 0, 3.1415, " - nuVertexFiducial");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_topologicalOrFlashMatch", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_topologicalOrFlashMatch", "rad", 20, 0, 3.1415, " - topologicalOrFlashMatch");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_min2Tracks", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_min2Tracks", "rad", 20, 0, 3.1415, " - min2Tracks");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_max1Uncontained", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_max1Uncontained", "rad", 20, 0, 3.1415, " - max1Uncontained");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_2NonProtons", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_2NonProtons", "rad", 20, 0, 3.1415, " - 2NonProtons");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_pionHasValiddEdx", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_pionHasValiddEdx", "rad", 20, 0, 3.1415, " - pionHasValiddEdx");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_pionNotInGap", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_pionNotInGap", "rad", 20, 0, 3.1415, " - pionNotInGap");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_muonNotInGap", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_muonNotInGap", "rad", 20, 0, 3.1415, " - muonNotInGap");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_openingAngle", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_openingAngle", "rad", 20, 0, 3.1415, " - openingAngle");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_topologicalScore", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_topologicalScore", "rad", 20, 0, 3.1415, " - topologicalScore");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_startNearVertex", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_startNearVertex", "rad", 20, 0, 3.1415, " - startNearVertex");
// plot_pur_hist( stv_ch, "cc1pi_reco_muonPionAngle", "true_cc1pi", "passed_likelyGoldenPion", "#theta_{#mu#pi}^{reco}", "openingangletrue_passed_likelyGoldenPion", "rad", 20, 0, 3.1415, " - likelyGoldenPion");



plot_eff_hist( stv_ch, "event_cutValue_openingAngle", "true_cc1pi", "passed_muonNotInGap", "#theta_{#mu#pi}^{reco}", "pmutrue_passed_muonNotInGap", "rad", 48, 0.0, 3.15, " - muonNotInGap cut");
plot_eff_hist( stv_ch, "event_cutValue_openingAngle", "true_cc1pi", "passed_openingAngle", "#theta_{#mu#pi}^{reco}", "pmutrue_passed_openingAngle", "rad", 48, 0.0, 3.15, " - openingAngle - cut");

plot_pur_hist( stv_ch, "event_cutValue_openingAngle", "true_cc1pi", "passed_muonNotInGap", "#theta_{#mu#pi}^{reco}", "pmutrue_passed_muonNotInGap", "rad", 48, 0.0, 3.15, " - muonNotInGap - cut");
plot_pur_hist( stv_ch, "event_cutValue_openingAngle", "true_cc1pi", "passed_openingAngle", "#theta_{#mu#pi}^{reco}", "pmutrue_passed_openingAngle", "rad", 48, 0.0, 3.15, " - openingAngle - cut");
}
