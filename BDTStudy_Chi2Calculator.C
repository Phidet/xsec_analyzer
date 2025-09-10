#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>

#include "MCC9SystematicsCalculator.hh"
#include "FilePropertiesManager.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"

#include "TFile.h"
#include "TH1D.h"

void BDTStudy_Chi2Calculator() 
{
    // Define input files structure and set inputs
    struct inputFiles {
        std::string rootFile;
        std::string fileList;
        std::string config;
        std::string sliceConfig;
        std::string nameExtension;
    };
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_run1234bcd5_5Feb25_bdts_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_removedZeroBins.txt",
    //     // "ubcc1pi_slice_config_bdt_removedZeroBins_experimental.txt",
    //     "_bdts_full_uncertainty"
    // };

    // // Coarser binning and all quality cuts applied
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_20May25_bdts_all_quality_cuts_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts.txt",
    //     "_bdts_full_uncertainty"
    // };

    // // Even coarser binning and all quality cuts applied
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_21May25_bdts_all_quality_cuts_17bins_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_17bins.txt",
    //     "_bdts_full_uncertainty"
    // };

    // // Even coarser binning and all quality cuts applied + sight binning edge change for proton BDT
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_21May25_bdts_all_quality_cuts_17bins_changedProtonBDTBinning_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_17bins.txt",
    //     "_bdts_full_uncertainty"
    // };


    // Coarser binning (14) and all quality cuts applied
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_27May25_bdts_all_quality_cuts_14bins_changedProtonBDTBinning_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_14bins.txt",
    //     "_bdts_full_uncertainty"
    // };


    // // Coarser binning (17) and all quality cuts applied + only including events and particles that the BDTs are applied to (fixed version)
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_28May25_bdts_all_quality_cuts_17bins_onlyAppliedParticles_fixed.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_17bins_highestMuonBDTScore.txt",
    //     "_bdts_full_uncertainty_highestMuonBDTScore"
    // };

    // // Coarser binning (15) only including events and particles that the BDTs are applied to (fixed version)
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_29May25_bdts_all_quality_cuts_15bins_onlyAppliedParticles_fixed.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_15bins_highestMuonBDTScore.txt",
    //     "_bdts_full_uncertainty_highestMuonBDTScore"
    // };

    // // Coarser binning (15) only including events and particles that the BDTs are applied to (fixed version) + only muon BDT plot using only the first particle in each event
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_30May25_bdts_all_quality_cuts_15bins_onlyAppliedParticles_randomSampling.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_15bins_randomSampling.txt",
    //     "_bdts_full_uncertainty_highestMuonBDTScore"
    // };

    // 22 and 25 bin binning to align with the cut lines 
    inputFiles input{
        "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_30May25_bdts_all_quality_cuts_22And25bins_onlyAppliedParticles_fixed.root",
        "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore.txt",
        "systcalc.conf",
        "ubcc1pi_slice_config_bdt_all_quality_cuts_22And25bins_highestMuonBDTScore.txt",
        "_bdts_full_uncertainty_highestMuonBDTScore"
    };

    // Load file properties
    auto &fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties(input.fileList);

    // Create the systematics calculator and slice binning object
    auto* syst_ptr = new MCC9SystematicsCalculator(input.rootFile, input.config);
    syst_ptr->set_syst_mode(syst_ptr->SystMode::VaryBackgroundAndSignalDirectly);
    auto& syst = *syst_ptr;
    auto* sb_ptr = new SliceBinning(input.sliceConfig);
    auto& sb = *sb_ptr;

    // Get the covariance matrices (release pointer for later deletion)
    auto* matrix_map_ptr = syst.get_covariances().release();
    auto& matrix_map = *matrix_map_ptr;

    // Get the data histogram and the MC+EXT prediction histogram
    TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
    TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
    TH1D* reco_mc_plus_ext_hist = dynamic_cast<TH1D*>(reco_ext_hist->Clone("reco_mc_plus_ext_hist"));
    reco_mc_plus_ext_hist->SetDirectory(nullptr);
    reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

    // Loop over the three BDTs and compute the chi2 values.
    std::vector<std::string> bdts = {"muonBDTScore", "protonBDTScore", "goldenPionBDTScore"};
    for (size_t i = 0; i < bdts.size(); i++)
    {
        const std::string &bdt = bdts[i];
        // Use corresponding slice: index 0 for muon, 1 for proton, 2 for golden pion.
        auto &slice = sb.slices_[i];

        // Create slice histograms using the covariance matrices.
        SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(*reco_bnb_hist, slice, &matrix_map.at("BNBstats"));
        SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(*reco_mc_plus_ext_hist, slice, &matrix_map.at("total"));

        auto chi2_result = slice_mc_plus_ext->get_chi2(*slice_bnb);

        // Print the chi2 result.
        std::cout << "BDT: " << bdt
                  << " -> chi²: " << std::fixed << std::setprecision(2) << chi2_result.chi2_
                  << " / " << chi2_result.num_bins_ << " bins, p-value: " << chi2_result.p_value_
                  << std::endl;

        // Print the modified chi2 result.
        if(i==1)
        {
            const std::vector<int> excluded_bins = {20,21,22};
            auto chi2_result_excluded_bins = slice_mc_plus_ext->get_chi2_excluded_bins(*slice_bnb, excluded_bins);

            std::cout << "BDT: " << bdt
                      << " -> chi² (excluded bins): " << std::fixed << std::setprecision(2) << chi2_result_excluded_bins.chi2_
                      << " / " << chi2_result_excluded_bins.num_bins_ << " bins, p-value: " << chi2_result_excluded_bins.p_value_
                      << std::endl;
        }


        delete slice_bnb;
        delete slice_mc_plus_ext;
    }

    // Clean up
    delete matrix_map_ptr;
    delete sb_ptr;
    delete syst_ptr;
}