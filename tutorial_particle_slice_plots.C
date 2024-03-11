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

void scale_by_bin_width(SliceHistogram* pSlice)
{
    int num_slice_bins = pSlice->hist_->GetNbinsX();
    TMatrixD trans_mat( num_slice_bins, num_slice_bins );
    for ( int b = 0; b < num_slice_bins; ++b ) {
      const auto width = pSlice->hist_->GetBinWidth( b + 1 );
      // width *= other_var_width;
      trans_mat( b, b ) = 1 / width;
    }
    pSlice->transform(trans_mat);
}

void slice_plots(const bool normaliseByBinWidth) {

  std::cout<<"DEBUG tutorial_slice_plots Point 0"<<std::endl;
  #ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro
    // MC ntuples as if they were data
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "nuwro_file_properties.txt" );
    auto* syst_ptr = new MCC9SystematicsCalculator(
      "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_run1234bcd5_6Mar24.root", // <-- Yes the name is wrong and should say nuwro
      "systcalc_unfold_fd_min.conf" );
    std::string nameExtension = "_fd";
  #else
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "file_properties.txt" );
    auto* syst_ptr = new MCC9SystematicsCalculator(
    "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/bdt_input_variables_truncMeandEdx_run1234bcd5.root", // golden pion cut plot with full uncertainties
    "systcalc.conf" );
    // std::string nameExtension = "_bnb";
    std::string nameExtension = "_bnb_truncMeandEdx_run1234bcd5";
  #endif

  std::cout<<"DEBUG tutorial_slice_plots Point 1"<<std::endl;
  auto& syst = *syst_ptr;

  auto* sb_ptr = new SliceBinning( "bdt_input_slice_config_truncMeandEdx.txt" );

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  // #ifdef USE_FAKE_DATA
  //   // Add the EXT to the "data" when working with fake data
  //   reco_bnb_hist->Add( reco_ext_hist );
  // #endif

  std::cout<<"DEBUG tutorial_slice_plots Point 2"<<std::endl;

  // TH2D* category_hist = syst.cv_universe().hist_categ_.get();

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

  auto& sb = *sb_ptr;

  for (const auto& pair : matrix_map) {
    const std::string& key = pair.first;
    std::cout << "Key: " << key << std::endl;
  }

  std::cout<<"DEBUG tutorial_slice_plots Point 3"<<std::endl;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
    std::cout<<"DEBUG tutorial_slice_plots Point 3.1 sl_idx: "<<sl_idx<<std::endl;
    // if(sl_idx!=2) continue;

    const auto& slice = sb.slices_.at( sl_idx );
    std::cout<<"DEBUG tutorial_slice_plots Point 3.2"<<std::endl;

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    std::cout<<"DEBUG tutorial_slice_plots Point 3.21"<<std::endl;

    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );

    std::cout<<"DEBUG tutorial_slice_plots Point 3.22"<<std::endl;

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    std::cout<<"DEBUG tutorial_slice_plots Point 3.23"<<std::endl;

    // auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    // std::cout << "Slice " << sl_idx << ": \u03C7\u00b2 = "
    //  << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
    //  << " p-value = " << chi2_result.p_value_ << '\n';

    // Prepare the plot legend
    // TLegend* lg = new TLegend( 0.3, 0.4 );
    TLegend* lg = new TLegend( 0.75, 0.7, 0.99, 0.9);
    std::cout<<"DEBUG tutorial_slice_plots Point 3.3"<<std::endl;

    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    const auto& eci = EventCategoryInterpreter::Instance();
    eci.set_ext_histogram_style( slice_ext->hist_.get() );

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );
    if(normaliseByBinWidth) scale_by_bin_width(slice_ext);
    slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

    // const auto& cat_map = eci.label_map();
    // // Go in reverse so that signal ends up on top. Note that this index is
    // // one-based to match the ROOT histograms
    // auto cat_bin_index = cat_map.size();
    const auto total_events = slice_mc_plus_ext->hist_->Integral();

    std::cout<<"DEBUG tutorial_slice_plots Point 4"<<std::endl;

    // for ( auto iter = cat_map.begin(); iter != cat_map.end(); ++iter )
    // {
    //   EventCategory cat = iter->first;
    //   std::cout<<"DEBUG std::to_string( cat ): "<<std::to_string( cat )<<" vs cat_bin_index: "<<cat_bin_index<<std::endl;
    //   TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
    //     cat+1, cat+1 );
    //   temp_mc_hist->SetDirectory( nullptr );

    //   SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
    //     *temp_mc_hist, slice  );

    //   eci.set_mc_histogram_style( cat, temp_slice_mc->hist_.get() );

    //   const auto label = iter->second;

    //   if(normaliseByBinWidth) scale_by_bin_width(temp_slice_mc);

    //   if(cat == kExternal)
    //   {
    //     // const auto events_in_category = slice_ext->hist_->Integral();
    //     // const auto category_percentage = events_in_category * 100. / total_events;
    //     // const auto cat_pct_label = Form( "%.2f%#%", category_percentage );
    //     lg->AddEntry( slice_ext->hist_.get(), (label).c_str(), "f" ); // + ", " + cat_pct_label
    //   }
    //   else if(cat != kUnknown)
    //   {
    //     // const auto events_in_category = temp_slice_mc->hist_->Integral();
    //     // const auto category_percentage = events_in_category * 100. / total_events;
    //     // const auto cat_pct_label = Form( "%.2f%#%", category_percentage );
    //     lg->AddEntry( temp_slice_mc->hist_.get(), (label).c_str(), "f" ); // + ", " + cat_pct_label
    //   }
    //   slice_pred_stack->Add( temp_slice_mc->hist_.get() );

    //   --cat_bin_index;
    // }

    std::cout<<"DEBUG tutorial_slice_plots Point 5"<<std::endl;

    TCanvas* c1 = new TCanvas;
    c1->SetRightMargin(0.252); // Allow space for the legend
    c1->SetLeftMargin(0.12); // Allow a bit more space for the y axis label on the left
    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 2 );
    // slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    // slice_bnb->hist_->SetMarkerSize( 0.8 );
    slice_bnb->hist_->SetStats( false );
    if(normaliseByBinWidth) scale_by_bin_width(slice_bnb);
    if(normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext);
    double ymax = 0;
    for (int i = 1; i <= slice_bnb->hist_->GetNbinsX(); ++i) {
        double binContent = slice_bnb->hist_->GetBinContent(i);
        double binError = slice_bnb->hist_->GetBinError(i);
        if (binContent + binError > ymax) {
            ymax = binContent + binError;
        }
    }
    for (int i = 1; i <= slice_mc_plus_ext->hist_->GetNbinsX(); ++i) {
        double binContent = slice_mc_plus_ext->hist_->GetBinContent(i);
        double binError = slice_mc_plus_ext->hist_->GetBinError(i);
        if (binContent + binError > ymax) {
            ymax = binContent + binError;
        }
    }
    ymax *= 1.07;
    slice_bnb->hist_->GetYaxis()->SetRangeUser(0., ymax);

    slice_bnb->hist_->Draw( "e" );
    // slice_pred_stack->Draw( "hist same" );

    slice_mc_plus_ext->hist_->SetLineWidth( 3 );
    slice_mc_plus_ext->hist_->SetFillColor(kGray + 1);
    slice_mc_plus_ext->hist_->SetFillStyle(3244);
    slice_mc_plus_ext->hist_->Draw( "same E2" );

    // Clone slice_mc_plus_ext and draw it as a line
    TH1F *slice_mc_plus_ext_line = (TH1F*)slice_mc_plus_ext->hist_->Clone();
    slice_mc_plus_ext_line->SetFillColor(0);
    slice_mc_plus_ext_line->SetFillStyle(0);
    slice_mc_plus_ext_line->SetLineColor(kGray + 1);
    slice_mc_plus_ext_line->SetLineWidth(2);
    slice_mc_plus_ext_line->Draw("same hist");

    lg->AddEntry( slice_mc_plus_ext_line, "MC & EXT", "l" );

    // // Create a dummy histogram with the same style as slice_mc_plus_ext
    // TH1F *dummy = new TH1F(*(TH1F*)slice_mc_plus_ext->hist_.get());
    // dummy->SetFillColor(kGray + 2);
    // dummy->SetFillStyle(3003);
    // lg->AddEntry(dummy, "MC & EXT Uncertainty", "f");

    slice_bnb->hist_->Draw( "same e" );
    lg->AddEntry( slice_bnb->hist_.get(), "Data (beam on)", "lp" );
    slice_bnb->hist_->SetTitle("");
    const std::string y_title = normaliseByBinWidth ? "# Reconstructed Particles / Bin width" : "# Reconstructed Particles";
    slice_bnb->hist_->GetYaxis()->SetTitle(y_title.c_str());

    // std::ostringstream oss;
    // auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    // oss << "#splitline{#chi^{2} = " << std::setprecision( 3 ) << chi2_result.chi2_ << " / "
    // << chi2_result.num_bins_ << " bin";
    // if ( chi2_result.num_bins_ > 1 ) oss << "s";
    // oss<<"}{";
    // if(chi2_result.num_bins_ > 1) oss<<"p-value = " << chi2_result.p_value_<<"}";
    // else oss<<"}";
    // const auto title =  oss.str();

    // // std::string legend_title = get_legend_title( pot_on );
    // lg->SetHeader( title.c_str(), "C" );

    lg->SetBorderSize( 0 );
    // Increase the font size for the legend header
    // (see https://root-forum.cern.ch/t/tlegend-headers-font-size/14434)
    TLegendEntry* lg_header = dynamic_cast< TLegendEntry* >(
      lg->GetListOfPrimitives()->First() );
    lg_header->SetTextSize( 0.03 );
    lg->Draw( "same" );

    std::string out_pdf_name = "plots/plot_slice_";
    if ( sl_idx < 10 ) out_pdf_name += "0";
    out_pdf_name += std::to_string( sl_idx ) + nameExtension;
    out_pdf_name += normaliseByBinWidth ? "_norm.pdf" : ".pdf";
    c1->SaveAs( out_pdf_name.c_str() );

    std::cout<<"DEBUG tutorial_slice_plots Point 6"<<std::endl;

    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    TH1* slice_hist = dynamic_cast< TH1* >(
      slice.hist_->Clone("slice_hist") );

    slice_hist->SetDirectory( nullptr );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Show fractional uncertainties computed using these covariance matrices
    // in the ROOT plot. All configured fractional uncertainties will be
    // included in the output pgfplots file regardless of whether they appear
    // in this vector.
    
    std::vector< std::string > cov_mat_keys = { "total", "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets", "MCstats", "EXTstats"};
    // std::vector< std::string > cov_mat_keys = { "total", "detVar_total", "flux", "reint", "xsec_multi", "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH", "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad", "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie", "POT", "numTargets", "MCstats", "EXTstats", "BNBstats"};

    #ifdef USE_FAKE_DATA
    // cov_mat_keys = { "total", "xsec_total", "MCstats", "EXTstats", "BNBstats"};
    cov_mat_keys = { "total", "MCstats", "EXTstats", "BNBstats", "xsec_multi", "xsec_unisim", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie"};
    // cov_mat_keys = { "total", "xsec_multi", "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH", "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad", "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie", "MCstats", "EXTstats", "BNBstats"};
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

      if ( color <= 9 ) ++color;
      if ( color == 5 ) ++color;
      if ( color >= 10 ) color += 10;

      slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );
    }

    std::cout<<"DEBUG tutorial_slice_plots Point 7"<<std::endl;

    TCanvas* c2 = new TCanvas;
    // TLegend* lg2 = new TLegend( 0.7, 0.7, 0.9, 0.9 );
    TLegend* lg2 = new TLegend( 0.2, 0.3);

    std::cout<<"DEBUG tutorial_slice_plots Point 7.1"<<std::endl;

    auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
    std::cout<<"DEBUG tutorial_slice_plots Point 7.2"<<std::endl;
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
      total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineStyle( 9 );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->Draw( "hist" );
    total_frac_err_hist->SetTitle("Fractional Uncertainty");
    total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");

    std::cout<<"DEBUG tutorial_slice_plots Point 7.3"<<std::endl;

    // const auto frac_ymax = 0.35;
    // total_frac_err_hist->GetYaxis()->SetRangeUser( 0., frac_ymax);

    lg2->AddEntry( total_frac_err_hist, "total", "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      std::cout<<"DEBUG tutorial_slice_plots Point 7.4"<<std::endl;
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == "total" ) continue;
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

    std::string frac_out_pdf_name = "plots/plot_frac_slice_";
    if ( sl_idx < 10 ) frac_out_pdf_name += "0";
    frac_out_pdf_name += std::to_string( sl_idx ) + nameExtension +".pdf";
    c2->SaveAs( frac_out_pdf_name.c_str() );

    std::cout << "Total frac error in bin #1 = "
      << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";

  } // slices

  std::cout<<"DEBUG tutorial_slice_plots Point 8"<<std::endl;

}

int tutorial_particle_slice_plots() {
  slice_plots(false);
  return 0;
}
