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

// #define USE_FAKE_DATA "no"

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

  // ################################################
  // Select one of the slice configurations here  
  // ################################################
  constexpr unsigned int configIndex = 0;
  // ################################################
  // ################################################

  const std::vector<std::string> rootFiles = {
    "bdt_input_variables_logBragg_pToMIP_run1.root",
    "",
    "",
    "",
    "",
    "",
  };

  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 0"<<std::endl;
  #ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro
    // MC ntuples as if they were data
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "file_properties_run1.txt" );
    auto* syst_ptr = new MCC9SystematicsCalculator(
      "/uboone/data/users/jdetje/ubcc1pi_univmake/" + rootFiles.at(configIndex),
      "systcalc_unfold_fd_min.conf" );
    std::string nameExtension = "_fd";
  #else
    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 0.1"<<std::endl;
    auto& fpm = FilePropertiesManager::Instance();
    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 0.2"<<std::endl;
    fpm.load_file_properties( "file_properties_run1.txt" );
    auto* syst_ptr = new MCC9SystematicsCalculator(
      "/uboone/data/users/jdetje/ubcc1pi_univmake/" + rootFiles.at(configIndex),
      "systcalc_noDetVar.conf" );
    std::string nameExtension = "_bnb";
  #endif

  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 0.3"<<std::endl;

  const std::vector<std::pair<std::string, std::string>> sliceConfigs = {
    {"bdt_input_slice_config_logBragg_pToMIP.txt", "logBragg_pToMIP"},
    {"bdt_input_slice_config_logBragg_piToMIP.txt", "logBragg_piToMIP"},
    {"bdt_input_slice_config_nDescendents.txt", "nDescendents"},
    {"bdt_input_slice_config_trackScore.txt", "trackScore"},
    {"bdt_input_slice_config_truncMeandEdx.txt", "truncMeandEdx"},
    {"bdt_input_slice_config_wiggliness.txt", "wiggliness"}
  };

  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 0.4"<<std::endl;

  nameExtension += "_" + sliceConfigs.at(configIndex).second;
  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 0.5"<<std::endl;
  const auto path = "./" + sliceConfigs.at(configIndex).first;
  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 0.6"<<std::endl;
  auto* sb_ptr = new SliceBinning( path );
  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 0.7"<<std::endl;
  auto& sb = *sb_ptr;

  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 1"<<std::endl;
  auto& syst = *syst_ptr;

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  // #ifdef USE_FAKE_DATA
  //   // Add the EXT to the "data" when working with fake data
  //   reco_bnb_hist->Add( reco_ext_hist );
  // #endif

  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 2"<<std::endl;

  // TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 2.1"<<std::endl;

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 2.2"<<std::endl;

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 2.2.1"<<std::endl;
  auto& matrix_map = *matrix_map_ptr;

  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 2.3"<<std::endl;

  for (const auto& pair : matrix_map) {
    const std::string& key = pair.first;
    std::cout << "Key: " << key << std::endl;
  }

  std::cout<<"DEBUG bdt_input_variable_slice_plots Point 3"<<std::endl;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 3.1 sl_idx: "<<sl_idx<<std::endl;
    // if(sl_idx!=2) continue;

    const auto& slice = sb.slices_.at( sl_idx );
    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 3.2"<<std::endl;

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 3.21"<<std::endl;

    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );

    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 3.22"<<std::endl;

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 3.23"<<std::endl;

    // auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    // std::cout << "Slice " << sl_idx << ": \u03C7\u00b2 = "
    //  << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
    //  << " p-value = " << chi2_result.p_value_ << '\n';

    // Prepare the plot legend
    // TLegend* lg = new TLegend( 0.3, 0.4 );
    TLegend* lg = new TLegend( 0.75, 0.1, 0.99, 0.9);
    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 3.3"<<std::endl;

    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    const auto& eci = EventCategoryInterpreter::Instance();
    eci.set_ext_histogram_style( slice_ext->hist_.get() );

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );
    if(normaliseByBinWidth) scale_by_bin_width(slice_ext);
    slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

    const auto& cat_map = eci.label_map();
    // Go in reverse so that signal ends up on top. Note that this index is
    // one-based to match the ROOT histograms
    auto cat_bin_index = cat_map.size();
    const auto total_events = slice_mc_plus_ext->hist_->Integral();

    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 4"<<std::endl;

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

    std::cout<<"DEBUG bdt_input_variable_slice_plots Point 5"<<std::endl;

    TCanvas* c1 = new TCanvas;
    c1->SetRightMargin(0.252); // Allow space for the legend
    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 3 );
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 0.8 );
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
    slice_pred_stack->Draw( "hist same" );

    slice_mc_plus_ext->hist_->SetLineWidth( 3 );
    slice_mc_plus_ext->hist_->SetFillColor(kGray + 2);
    slice_mc_plus_ext->hist_->SetFillStyle(3003);
    slice_mc_plus_ext->hist_->Draw( "same E2" );
    // Create a dummy histogram with the same style as slice_mc_plus_ext
    TH1F *dummy = new TH1F(*(TH1F*)slice_mc_plus_ext->hist_.get());
    dummy->SetFillColor(kGray + 2);
    dummy->SetFillStyle(3003);
    lg->AddEntry(dummy, "MC & EXT Uncertainties", "f");

    slice_bnb->hist_->Draw( "same e" );
    lg->AddEntry( slice_bnb->hist_.get(), "Data (beam on)", "lp" );
    slice_bnb->hist_->SetTitle("Selected #nu_{#mu}CC1#pi^{#pm}Xp, N #geq 0 Events");
    const std::string y_title = normaliseByBinWidth ? "Events / Bin width" : "Events";
    slice_bnb->hist_->GetYaxis()->SetTitle(y_title.c_str());

    std::ostringstream oss;
    auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
    oss << "#splitline{#chi^{2} = " << std::setprecision( 3 ) << chi2_result.chi2_ << " / "
    << chi2_result.num_bins_ << " bin";
    if ( chi2_result.num_bins_ > 1 ) oss << "s";
    oss<<"}{";
    if(chi2_result.num_bins_ > 1) oss<<"p-value = " << chi2_result.p_value_<<"}";
    else oss<<"}";
    const auto title =  oss.str();

    // std::string legend_title = get_legend_title( pot_on );
    lg->SetHeader( title.c_str(), "C" );

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
  } // End of loop over slices
} // End of function slice_plots

int bdt_input_variable_slice_plots() {
  slice_plots(false);
  return 0;
}
