#include <algorithm>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

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

void scale_by_bin_width(SliceHistogram* pSlice)
{
    int num_slice_bins = pSlice->hist_->GetNbinsX();
    TMatrixD trans_mat( num_slice_bins, num_slice_bins );
    for ( int b = 0; b < num_slice_bins; ++b ) {
        const auto width = pSlice->hist_->GetBinWidth( b + 1 );
        trans_mat( b, b ) = 1 / width;
    }
    pSlice->transform(trans_mat);
}

struct inputFiles
{
  std::string rootFile;
  std::string fileList;
  std::string config;
  std::string sliceConfig;
  std::string nameExtension;
};

void make_slice_plots(const bool normaliseByBinWidth) {
  // inputFiles input{
  //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_run1234bcd5_5Feb25_bdts_full_uncertainty.root",
  //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
  //     "systcalc.conf",
  //     "ubcc1pi_slice_config_bdt.txt",
  //     "_bdts_full_uncertainty"
  // };

  // // 22 and 25 bin binning to align with the cut lines; random sampling of particles
  // inputFiles input{
  //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_16July25_bdts_all_quality_cuts_22And25bins_onlyAppliedParticles_randomIndex.root",
  //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore_randomIndex.txt",
  //     "systcalc.conf",
  //     "ubcc1pi_slice_config_bdt_all_quality_cuts_22And25bins_highestMuonBDTScore.txt",
  //     "_bdts_full_uncertainty_highestMuonBDTScore"
  // };

  // // 22 and 25 bin binning to align with the cut lines; random sampling of particles; includes all later quality cuts.
  // inputFiles input{
  //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_17July25_bdts_all_quality_cuts_22And25bins_onlyAppliedParticles_randomIndex.root",
  //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore_randomIndex.txt",
  //     "systcalc.conf",
  //     "ubcc1pi_slice_config_bdt_all_quality_cuts_22And25bins_highestMuonBDTScore.txt",
  //     "_bdts_full_uncertainty_highestMuonBDTScore"
  // };


  // 22 and 25 bin binning to align with the cut lines; random sampling of particles; includes all later quality cuts; fixed vertex distance variable
  inputFiles input{
      "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_21July25_bdts_all_quality_cuts_22And25bins_onlyAppliedParticles_randomIndex.root",
      "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore_randomIndex_fixedTrackDistance.txt",
      "systcalc.conf",
      "ubcc1pi_slice_config_bdt_all_quality_cuts_22And25bins_highestMuonBDTScore.txt",
      "_bdts_full_uncertainty_highestMuonBDTScore_fixedTrackDistance"
  };

  auto &fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties(input.fileList);
  auto* syst_ptr = new MCC9SystematicsCalculator(input.rootFile, input.config);
  syst_ptr->set_syst_mode(syst_ptr->SystMode::VaryBackgroundAndSignalDirectly);
  std::string nameExtension = input.nameExtension;
  auto& syst = *syst_ptr;
  auto* sb_ptr = new SliceBinning( input.sliceConfig );
  auto& sb = *sb_ptr;
  
  // Get the required histograms
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
  
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >( reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );
  
  // Retrieve the covariance matrices (only needed for stacking)
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  // Loop over slices and produce stacked histogram plots
  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
    const auto& slice = sb.slices_.at( sl_idx );
    
    // Create slice histograms for data, EXT, and MC+EXT
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(*reco_bnb_hist, slice, &matrix_map.at("BNBstats"));
    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(*reco_ext_hist, slice, &matrix_map.at("EXTstats"));
    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(*reco_mc_plus_ext_hist, slice, &matrix_map.at("total"));
    
    // Build a stack: first add EXT and then MC groups
    THStack* slice_pred_stack = new THStack("mc+ext", "");
    if (normaliseByBinWidth)
      scale_by_bin_width(slice_ext);
    slice_pred_stack->Add(slice_ext->hist_.get());
    
    // Group definitions for MC components
    std::map<std::string, std::vector<EventCategory>> group_map = {
      { "#nu_{#mu}/#bar{#nu}_{#mu} CC1#pi^{#pm} (unscattered #pi^{#pm})", {kNumuCC1PiChargedGolden} },
      { "#nu_{#mu}/#bar{#nu}_{#mu} CC1#pi^{#pm} (scattered #pi^{#pm})", {kNumuCC1PiChargedNonGolden} },
      { "#nu_{#mu}/#bar{#nu}_{#mu} CC1#pi^{#pm} Non-Signal", {kNumuCC1PiNonSignal} },
      { "#nu_{#mu}/#bar{#nu}_{#mu} CC0#pi", {kNumuCC0PiSignal, kNumuCC0Pi} },
      { "#nu_{#mu}/#bar{#nu}_{#mu} CC1#pi^{0}", {kNumuCC1PiZero} },
      { "Other #nu_{#mu}/#bar{#nu}_{#mu} CC", {kNumuCCOther} },
      { "NC", {kNC} },
      { "Other (Non-Fiducial + Dirt + #nu_{e}/#bar{#nu}_{e})", {kNonFiducial, kDirt, kNue} }
    };
    std::vector<std::string> group_order = {
      "Other (Non-Fiducial + Dirt + #nu_{e}/#bar{#nu}_{e})",
      "NC",
      "Other #nu_{#mu}/#bar{#nu}_{#mu} CC",
      "#nu_{#mu}/#bar{#nu}_{#mu} CC1#pi^{0}",
      "#nu_{#mu}/#bar{#nu}_{#mu} CC0#pi",
      "#nu_{#mu}/#bar{#nu}_{#mu} CC1#pi^{#pm} Non-Signal",
      "#nu_{#mu}/#bar{#nu}_{#mu} CC1#pi^{#pm} (scattered #pi^{#pm})",
      "#nu_{#mu}/#bar{#nu}_{#mu} CC1#pi^{#pm} (unscattered #pi^{#pm})"
    };
    
    // Legend entries (added only for first slice)
    std::vector<std::pair<TH1D*, std::string>> legend_entries;
    
    const auto& eci = EventCategoryInterpreter::Instance();
    eci.set_ext_histogram_style(slice_ext->hist_.get());
    
    for (const auto& group_name : group_order) {
      const auto& group = group_map[group_name];
      TH1D* group_hist = nullptr;
      for (const auto& cat : group) {
          TH1D* temp_mc_hist = syst.cv_universe().hist_categ_->ProjectionY("temp_mc_hist", cat + 1, cat + 1);
          if (!group_hist)
              group_hist = (TH1D*)temp_mc_hist->Clone();
          else
              group_hist->Add(temp_mc_hist);
      }
      if (!group_hist) continue;
      
      SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(*group_hist, slice);
      eci.set_mc_histogram_style(group.front(), temp_slice_mc->hist_.get());
      if (normaliseByBinWidth)
          scale_by_bin_width(temp_slice_mc);
      slice_pred_stack->Add(temp_slice_mc->hist_.get());
      if (sl_idx == 0 && group_name != "External")
          legend_entries.push_back(std::make_pair((TH1D*)temp_slice_mc->hist_.get(), group_name));
      delete group_hist;
    }
    
    // Create and draw the stack plot
    TCanvas* c1 = new TCanvas("c1", "Stacked Histogram", 600, 400);
    c1->SetRightMargin(0.05);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.12);
    c1->SetTopMargin(0.07);
    
    if (normaliseByBinWidth) {
      scale_by_bin_width(slice_bnb);
      scale_by_bin_width(slice_mc_plus_ext);
    }
    
    slice_bnb->hist_->SetLineColor(kBlack);
    slice_bnb->hist_->SetLineWidth(3);
    slice_bnb->hist_->SetMarkerStyle(kFullCircle);
    slice_bnb->hist_->SetMarkerSize(0.8);
    slice_bnb->hist_->SetStats(false);
    
    // Draw in order: data, stack, and MC+EXT uncertainty band
    slice_bnb->hist_->Draw("e");
    slice_pred_stack->Draw("hist same");
    slice_mc_plus_ext->hist_->SetLineWidth(3);
    slice_mc_plus_ext->hist_->SetFillColor(kGray + 1);
    slice_mc_plus_ext->hist_->SetFillStyle(3244);
    slice_mc_plus_ext->hist_->Draw("same E2");
    slice_bnb->hist_->Draw("same e");
    
    // // Create legend only for the first slice
    // if(sl_idx == 0) {
    //   TLegend* lg = new TLegend(0.09, 0.9, 0.99, 0.99);
    //   for (auto it = legend_entries.rbegin(); it != legend_entries.rend(); ++it)
    //     lg->AddEntry(it->first, it->second.c_str(), "f");
    //   lg->AddEntry(slice_ext->hist_.get(), "External", "f");
    //   lg->AddEntry(slice_bnb->hist_.get(), "Data", "lp");
    //   lg->Draw();
    // }
    
    // Set plot aesthetics
    slice_bnb->hist_->SetTitle("");
    const std::string y_title = normaliseByBinWidth ? "# Events / Bin width" : "# Events";
    slice_bnb->hist_->GetYaxis()->SetTitle(y_title.c_str());
    
    c1->Update();
    
    // Save the stacked histogram plot
    std::string out_pdf_name = "plots/plot_slice_";
    out_pdf_name += (sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + nameExtension;
    out_pdf_name += normaliseByBinWidth ? "_norm.pdf" : ".pdf";

    // Set y-axis minimum to 0 and maximum to (max value + uncertainty) * 1.1
    double max_val = 0.0;
    for (int bin = 1; bin <= slice_mc_plus_ext->hist_->GetNbinsX(); ++bin) {
      double val = slice_mc_plus_ext->hist_->GetBinContent(bin);
      double err = slice_mc_plus_ext->hist_->GetBinError(bin);
      if (val + err > max_val)
        max_val = val + err;
    }
    slice_bnb->hist_->GetYaxis()->SetRangeUser(0, max_val * 1.1);

    c1->SaveAs(out_pdf_name.c_str());
    
  } // slices loop
  
  std::cout << "-------- All done --------" << std::endl;
}

int slice_plots_bdt() {
  make_slice_plots(false);
  return 0;
}
