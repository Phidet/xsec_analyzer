// Standard library includes
#include <iostream>

// ROOT includes
#include "TFile.h"
#include "TH1D.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"

using NFT = NtupleFileType;

void print_uncertainties() {
    // Define input files and configurations
    std::string rootFile = "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_run1234bcd5_31Jan25_all_three_total_xsecs.root";
    std::string fileList = "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt";
    // std::string config = "systcalc.conf";
    std::string config = "systcalc_unfold.conf";
    std::string sliceConfig = "ubcc1pi_neutral_slice_config_all_three_total_xsecs.txt";

    // Initialize the FilePropertiesManager and load file properties
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties(fileList);

    // Initialize the SystematicsCalculator
    auto* syst_ptr = new MCC9SystematicsCalculator(rootFile, config);
    syst_ptr->set_syst_mode(syst_ptr->SystMode::VaryBackgroundAndSignalDirectly);
    auto& syst = *syst_ptr;

    // Initialize the SliceBinning
    auto* sb_ptr = new SliceBinning(sliceConfig);
    auto& sb = *sb_ptr;

    // Get the relevant histograms
    TH1D* reco_bnb_hist = syst.data_hists_.at(NFT::kOnBNB).get();
    TH1D* reco_ext_hist = syst.data_hists_.at(NFT::kExtBNB).get();
    TH1D* reco_mc_hist = syst.cv_universe().hist_reco_.get();
    TH2D* category_hist = syst.cv_universe().hist_categ_.get();

    // Histogram for MC prediction + ext
    TH1D* reco_mc_with_ext_hist = dynamic_cast<TH1D*>(reco_mc_hist->Clone("reco_mc_with_ext_hist"));
    reco_mc_with_ext_hist->Add(reco_ext_hist);
    reco_mc_with_ext_hist->SetDirectory(nullptr);

    // Create a histogram for the background (non-signal MC + ext)
    TH1D* reco_mc_background_hist = dynamic_cast<TH1D*>(reco_ext_hist->Clone("reco_mc_background_hist"));
    reco_mc_background_hist->SetDirectory(nullptr);
    // Create an empty histogram for the signal
    TH1D* reco_mc_signal_hist = dynamic_cast<TH1D*>(reco_mc_hist->Clone("reco_mc_signal_hist"));
    reco_mc_signal_hist->Reset();
    reco_mc_signal_hist->SetDirectory(nullptr);


    // Add the background categories to the background histogram
    const auto& eci = EventCategoryInterpreter::Instance();
    const auto& cat_map = eci.label_map();
    std::cout << "Before:" <<reco_mc_background_hist->GetBinContent(60) << std::endl;
    for (const auto& pair : cat_map) {
        EventCategory cat = pair.first;
        if(cat != kExternal)
        {
            if (cat != kNumuCC1PiChargedGolden && cat != kNumuCC1PiChargedNonGolden) {
                TH1D* temp_mc_hist = category_hist->ProjectionY("temp_mc_hist", cat + 1, cat + 1);
                reco_mc_background_hist->Add(temp_mc_hist);
                std::cout << " cat:" << cat << " " <<reco_mc_background_hist->GetBinContent(60) << std::endl;
                delete temp_mc_hist;
            }
            else {
                TH1D* temp_mc_hist = category_hist->ProjectionY("temp_mc_hist", cat + 1, cat + 1);
                reco_mc_signal_hist->Add(temp_mc_hist);
                delete temp_mc_hist;
            }
        }
    }

    // Get the covariance matrices
    auto* matrix_map_ptr = syst.get_covariances().release();
    auto& matrix_map = *matrix_map_ptr;

    // Loop over slices
    for (size_t sl_idx = 0; sl_idx < sb.slices_.size(); ++sl_idx) {
        const auto& slice = sb.slices_.at(sl_idx);

        // Create SliceHistograms for data and MC+EXT
        SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(*reco_bnb_hist, slice, &matrix_map.at("DataStats"));
        // SliceHistogram* slice_mc_background = SliceHistogram::make_slice_histogram(*reco_mc_background_hist, slice, &matrix_map.at("total"));
        // SliceHistogram* slice_mc_signal = SliceHistogram::make_slice_histogram(*reco_mc_signal_hist, slice, &matrix_map.at("total")); // matrix_map is just placeholder


        SliceHistogram* slice_mc_with_ext = SliceHistogram::make_slice_histogram(*reco_mc_with_ext_hist, slice, &matrix_map.at("total"));

        // Print absolute uncertainties for each bin
        std::cout << "Slice " << sl_idx << ":\n";
        for (int bin = 1; bin <= slice_bnb->hist_->GetNbinsX(); ++bin) {
            double data_value = slice_bnb->hist_->GetBinContent(bin);
            double data_stat_uncertainty = slice_bnb->hist_->GetBinError(bin);
            // double mc_background = slice_mc_background->hist_->GetBinContent(bin);
            // double mc_background_uncertainty = slice_mc_background->hist_->GetBinError(bin);
            // double mc_signal = slice_mc_signal->hist_->GetBinContent(bin);
            double mc_with_ext = slice_mc_with_ext->hist_->GetBinContent(bin);
            double mc_with_ext_uncertainty = slice_mc_with_ext->hist_->GetBinError(bin);
            std::cout << "  Bin " << bin << ": Data = " << data_value << " +- " << data_stat_uncertainty << 
                // ", Predicted Signal = " << mc_signal << ", Predicted Background = " << mc_background << " +- " << mc_background_uncertainty << " (syst + stat)" <<
                ", Predicted Total = " << mc_with_ext << " +- " << mc_with_ext_uncertainty << " (syst + stat)\n";
        }
    }

    std::cout << "-------- All done --------" << std::endl;
}

int main() {
    print_uncertainties();
    return 0;
}
