// Standard library includes
#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "PlotUtils.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"

using NFT = NtupleFileType;

// #define USE_FAKE_DATA "yes"

void get_slice_uncertainties() {

  std::cout<<"DEBUG tutorial_slice_plots Point 0"<<std::endl;
  #ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro
    // MC ntuples as if they were data
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "nuwro_file_properties.txt" );
    auto* syst_ptr = new MCC9SystematicsCalculator(
      "/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_10/univmake_output_nuwro_with_sideband_overflow_all_13Jan23.root",
      "systcalc_unfold_fd_min.conf" );
    std::string nameExtension = "_fd";
  #else
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "file_properties_testingOnly.txt" );
    auto* syst_ptr = new MCC9SystematicsCalculator(
    // "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_run1234bcd5_3Mar24.root",
    "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_run1234bcd5_20Mar24_testingOnly.root",
    "systcalc.conf" );
    std::string nameExtension = "_bnb_testingOnly";
    // std::string nameExtension = "_bnb_xsec_testingOnly";
    // std::string nameExtension = "_bnb_detVar_testingOnly";
  #endif

  std::cout<<"DEBUG tutorial_slice_plots Point 1"<<std::endl;
  auto& syst = *syst_ptr;

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  std::cout<<"DEBUG tutorial_slice_plots Point 2"<<std::endl;

  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Create new histograms for signal and background
  TH1D* mc_signal_hist = new TH1D("signal", "Signal", category_hist->GetNbinsX(), category_hist->GetXaxis()->GetXmin(), category_hist->GetXaxis()->GetXmax());
  TH1D* mc_background_hist = new TH1D("background", "Background", category_hist->GetNbinsX(), category_hist->GetXaxis()->GetXmin(), category_hist->GetXaxis()->GetXmax());

  // Combine categories for signal
  for (int bin = 1; bin <= category_hist->GetNbinsX(); ++bin) {
    double content = category_hist->GetBinContent(bin, 1) + category_hist->GetBinContent(bin, 2);
    mc_signal_hist->SetBinContent(bin, content);
  }

  // Combine categories for background
  for (int bin = 1; bin <= category_hist->GetNbinsX(); ++bin) {
    double content = 0;
    for (int category = 3; category <= 12; ++category) {
      content += category_hist->GetBinContent(bin, category);
    }
    mc_background_hist->SetBinContent(bin, content);
  }

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );


  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* sb_ptr = new SliceBinning( "ubcc1pi_neutral_slice_config.txt" );
  auto& sb = *sb_ptr;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
    const auto& slice = sb.slices_.at( sl_idx );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Show fractional uncertainties computed using these covariance matrices
    // in the ROOT plot. All configured fractional uncertainties will be
    // included in the output pgfplots file regardless of whether they appear
    // in this vector.
    
    std::vector< std::string > cov_mat_keys = { "total", "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets", "MCstats", "EXTstats"};
    const auto total_name = "total";
    // std::vector< std::string > cov_mat_keys = { "xsec_total", "xsec_multi", "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH", "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad", "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie"};
    // const auto total_name = "xsec_total";
    // std::vector< std::string > cov_mat_keys = {"detVar_total", "detVarLYatten", "detVarLYdown", "detVarLYrayl", "detVarRecomb2", "detVarSCE", "detVarWMAngleXZ", "detVarWMAngleYZ", "detVarWMX", "detVarWMYZ"};
    // const auto total_name = "detVar_total";
    // std::vector< std::string > cov_mat_keys = { "total", "detVar_total", "flux", "reint", "xsec_multi", "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH", "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad", "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie", "POT", "numTargets", "MCstats", "EXTstats"};

    #ifdef USE_FAKE_DATA
    // cov_mat_keys = { "total", "xsec_total", "MCstats", "EXTstats"};
    cov_mat_keys = { "total", "MCstats", "EXTstats", "xsec_multi", "xsec_unisim", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie"};
    // cov_mat_keys = { "total", "xsec_multi", "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH", "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad", "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie", "MCstats", "EXTstats"};
    #endif

    // Loop over the various systematic uncertainties
    int color = 0;
    for ( const auto& pair : matrix_map ) {

      const auto& key = pair.first;
      const auto& cov_matrix = pair.second;

      SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, &cov_matrix );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
        double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
        double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        double frac = 0.;
        if ( y > 0. ) frac = err / y;
        slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
        slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
      }

      // Check whether the current covariance matrix name is present in
      // the vector defined above this loop. If it isn't, don't bother to
      // plot it, and just move on to the next one.
      auto cbegin = cov_mat_keys.cbegin();
      auto cend = cov_mat_keys.cend();
      auto iter = std::find( cbegin, cend, key );
      if ( iter == cend )
      {
        std::cout << "DEBUG skipping " << key << std::endl;
        continue;
      }
      else
      {
        std::cout<<"DEBUG not skipping "<<key<<std::endl;
      }

      frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

      // if ( color <= 9 ) ++color;
      // if ( color == 5 ) ++color;
      if ( color >= 9 ) color += 7;
      else color++;
      //   slice_for_syst->hist_->SetLineStyle( 3 );
      // if(color > 20)
      //   slice_for_syst->hist_->SetLineStyle( 4 );

      slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );
    }

    std::cout<<"DEBUG tutorial_slice_plots Point 7"<<std::endl;

    TCanvas* c2 = new TCanvas;
    // TLegend* lg2 = new TLegend( 0.7, 0.7, 0.9, 0.9 );
    TLegend* lg2 = new TLegend( 0.2, 0.3);

    std::cout<<"DEBUG tutorial_slice_plots Point 7.1"<<std::endl;

    auto* total_frac_err_hist = frac_uncertainty_hists.at( total_name );
    std::cout<<"DEBUG tutorial_slice_plots Point 7.2"<<std::endl;
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
      total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineStyle( 2 );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->Draw( "hist" );
    total_frac_err_hist->SetTitle("Fractional Uncertainty of Selected #nu_{#mu}CC1#pi^{#pm}Xp, X #geq 0 Events");
    total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");

    std::cout<<"DEBUG tutorial_slice_plots Point 7.3"<<std::endl;

    // const auto frac_ymax = 0.35;
    // total_frac_err_hist->GetYaxis()->SetRangeUser( 0., frac_ymax);

    lg2->AddEntry( total_frac_err_hist, total_name, "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      std::cout<<"DEBUG tutorial_slice_plots Point 7.4"<<std::endl;
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == total_name ) continue;
      if (name.size() >= 5 && name.substr(name.size() - 5) == "stats") 
      {
        hist->SetLineStyle( 2 );
      }
      std::cout<<"DEBUG tutorial_slice_plots Point 7.5"<<std::endl;

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );


      std::cout << name << " frac err in bin #1 = "
        << hist->GetBinContent( 1 )*100. << "%\n";
    }

    lg2->Draw( "same" );

    std::string frac_out_name = "plots/get_slice_uncertainties_slice";
    if ( sl_idx < 10 ) frac_out_name += "0";
    frac_out_name += std::to_string( sl_idx ) + nameExtension;
    c2->SaveAs( (frac_out_name + ".pdf").c_str() );
    c2->SaveAs( (frac_out_name + ".png").c_str() );
    c2->SaveAs( (frac_out_name + ".C").c_str() );

    std::cout << "Total frac error in bin #1 = "
      << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";

  } // End of loop over slices

  std::cout<<"--------- All Done ---------"<<std::endl;

} // End of function get_slice_uncertainties
