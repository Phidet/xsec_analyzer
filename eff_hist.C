// Efficiency histogram plotter

const std::string signal_cuts = "mc_is_signal";
const std::string selection_cuts = "sel_CCNp0pi";

void plot_eff_hist( TTree& stv_tree, const std::string& branch,
  const std::string& variable_name, const std::string& hist_name,
  int num_bins, double x_min, double x_max)
{
  // For computing efficiencies and purities, we need to only use MC events.
  // Unconditionally add this requirement to the cuts defined above.
  std::string signal = signal_cuts + " && is_mc && genie_ok ";
  std::string selection = selection_cuts + " && is_mc && genie_ok ";

  std::string eff_title = hist_name + "; " + variable_name + "; efficiency";
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
    if ( num_signal > 0. ) eff = num_selected_signal / num_signal;

    eff_hist->SetBinContent( b, eff );
  }

  TCanvas* c1 = new TCanvas;
  eff_hist->SetStats( false );
  //eff_hist->GetYaxis()->SetRangeUser( 0., 1. );
  eff_hist->Draw();

}

void eff_hist() {

  TChain stv_ch( "stv_tree" );
  stv_ch.Add( "*stv.root" );

  //plot_eff_hist( stv_ch, "mc_p3_mu.Mag()", "p_{#mu}^{true}", "pmutrue", 10, 0., 2.);
  //plot_eff_hist( stv_ch, "mc_p3_mu.CosTheta()", "cos#theta_{#mu}^{true}", "cthmutrue", 10, -1., 1.);
  //plot_eff_hist( stv_ch, "mc_p3_mu.Phi()", "#phi_{#mu}^{true}", "phimutrue", 10, 0., M_PI);

  //plot_eff_hist( stv_ch, "mc_p3_lead_p.Mag()", "p_{lead p}^{true}", "pptrue", 10, 0., 2.);
  //plot_eff_hist( stv_ch, "mc_p3_lead_p.CosTheta()", "cos#theta_{lead p}^{true}", "cthptrue", 10, -1., 1.);
  //plot_eff_hist( stv_ch, "mc_p3_lead_p.Phi()", "#phi_{lead p}^{true}", "phiptrue", 10, 0., M_PI);

  plot_eff_hist( stv_ch, "mc_p3_lead_p.X() / mc_p3_lead_p.Mag()", "px_{lead p}^{true} / p_{lead p}^{true}", "pxp", 10., 0., 1.);
  plot_eff_hist( stv_ch, "mc_p3_lead_p.Y() / mc_p3_lead_p.Mag()", "py_{lead p}^{true} / p_{lead p}^{true}", "pyp", 10., 0., 1.);
  plot_eff_hist( stv_ch, "mc_p3_lead_p.Z() / mc_p3_lead_p.Mag()", "pz_{lead p}^{true} / p_{lead p}^{true}", "pzp", 10., 0., 1.);

  //plot_eff_hist( stv_ch, "mc_delta_pT", "#deltap_{T}", "deltapT", 10, 0., 2.);
  //plot_eff_hist( stv_ch, "mc_delta_phiT", "#delta#phi_{T}", "deltaphiT", 10, 0., M_PI);
  //plot_eff_hist( stv_ch, "mc_delta_alphaT", "#delta#alpha_{T}", "deltaalphaT", 10, 0., M_PI);
  //plot_eff_hist( stv_ch, "mc_delta_pL", "#deltap_{L}", "deltapL", 10, 0., 2.);
  //plot_eff_hist( stv_ch, "mc_pn", "p_{n}", "pn", 10, 0., 2.);

}