// Standard library includes
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"

// STV analysis includes
#include "../FilePropertiesManager.hh"
#include "../DAgostiniUnfolder.hh"
#include "../FiducialVolume.hh"
#include "../MatrixUtils.hh"
#include "../MCC9SystematicsCalculator.hh"
#include "../NormShapeCovMatrix.hh"
#include "../PGFPlotsDumpUtils.hh"
#include "../SliceBinning.hh"
#include "../SliceHistogram.hh"
#include "../WienerSVDUnfolder.hh"
#include "../utils.hh"
#include "../UnfolderHelper.hh"

constexpr double BIG_DOUBLE = 1e300;

void multiply_1d_hist_by_matrix(TMatrixD *mat, TH1 *hist)
{
    // Copy the histogram contents into a column vector
    int num_bins = mat->GetNcols();
    TMatrixD hist_mat(num_bins, 1);
    for (int r = 0; r < num_bins; ++r)
    {
        hist_mat(r, 0) = hist->GetBinContent(r + 1);
    }

    // Multiply the column vector by the input matrix
    // TODO: add error handling here related to matrix dimensions
    TMatrixD hist_mat_transformed(*mat, TMatrixD::EMatrixCreatorsOp2::kMult,
                                  hist_mat);

    // Update the input histogram contents with the new values
    for (int r = 0; r < num_bins; ++r)
    {
        double val = hist_mat_transformed(r, 0);
        hist->SetBinContent(r + 1, val);
    }
}


std::string toLatexScientific(double value) {
    std::stringstream stream;
    stream << std::scientific << std::setprecision(2) << value;
    std::string str = stream.str();
    size_t pos = str.find('e');
    if (pos != std::string::npos) {
        str.replace(pos, 1, " \\times 10^{");
        str += "}";
    }
    return str;
}
// const std::string SAMPLE_NAME = "MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC";

struct TruthFileInfo
{
    TruthFileInfo() {}
    TruthFileInfo(const std::string &file_name, int color, int style)
        : file_name_(file_name), color_(color), style_(style) {}

    std::string file_name_;
    int color_;
    int style_;
};

// Keys are generator names and versions, values are TruthFileInfo objects
// describing nuiscomp output files containing the differential cross-section
// predictions in each true bin
std::map<std::string, TruthFileInfo> truth_file_map = {};

std::map<std::string, std::string> samples_to_hist_names{
    {"Unfolded Selection", "UnfData"},
    {"MicroBooNE Tune", "uBTune"},
    {"NuWro Truth", "FakeData"}};

struct SampleInfo
{
    SampleInfo() {}

    SampleInfo(const std::string &resp, const std::string &width,
               const std::string &sb) : respmat_file_(resp), widths_file_(width),
                                        sb_file_(sb) {}

    // File containing the output of respmat.C with the same true binning
    // as was used in NUISANCE to make theoretical predictions
    std::string respmat_file_;

    // NUISANCE data file in which the bin widths (along each relevant axis) are
    // stored. These files are used by the preliminary NUISANCE implementation of
    // this cross-section measurement.
    std::string widths_file_;

    // SliceBinning configuration file (needs to be consistent with the bin
    // definitions used in the other two files)
    std::string sb_file_;
};

// Keys are NUISANCE sample names, values are SampleInfo objects that give
// the corresponding input file paths needed for this script
std::map<std::string, SampleInfo> sample_info_map = {};

void dump_slice_errors(const std::string &hist_col_prefix,
                       const Slice &slice, const std::map<std::string, std::unique_ptr<SliceHistogram>> &slice_hist_cov_matrix_map,
                       std::map<std::string, std::vector<double>> &pgf_plots_hist_table)
{
    for (const auto &pair : slice_hist_cov_matrix_map)
    {
        std::string err_name = pair.first;
        std::string err_col_name = hist_col_prefix + '_' + err_name + "_error";
        pgf_plots_hist_table[err_col_name] = std::vector<double>();
    }

    for (const auto &bin_pair : slice.bin_map_)
    {
        // TODO: revisit for multi-dimensional slices
        int global_bin_idx = bin_pair.first;

        for (const auto &err_pair : slice_hist_cov_matrix_map)
        {

            std::string err_name = err_pair.first;
            std::string err_col_name = hist_col_prefix + '_' + err_name + "_error";

            const auto *hist = err_pair.second->hist_.get();
            double err = hist->GetBinError(global_bin_idx);

            pgf_plots_hist_table.at(err_col_name).push_back(err);
        }

    } // slice bins

    // Add a (presumably empty) overflow bin to get certain PGFPlots styles to
    // look right.
    for (const auto &err_pair : slice_hist_cov_matrix_map)
    {
        std::string err_name = err_pair.first;
        std::string err_col_name = hist_col_prefix + '_' + err_name + "_error";

        pgf_plots_hist_table.at(err_col_name).push_back(0.);
    }
}

// Helper function that dumps a lot of the results to simple text files.
// The events_to_xsec_factor is a constant that converts expected true event
// counts to a total cross section (10^{-39} cm^2 / Ar) via multiplication.
void dump_overall_results(const UnfoldedMeasurement &result,
                          const std::map<std::string, std::shared_ptr<TMatrixD>> &unf_cov_matrix_map,
                          double events_to_xsec_factor, const TMatrixD &genie_cv_true_events,
                          const TMatrixD &fake_data_true_events,
                          const std::map<std::string, TMatrixD *> &generator_truth_map,
                          bool using_fake_data)
{
    // Dump the unfolded flux-averaged total cross sections (by converting
    // the units on the unfolded signal event counts)
    TMatrixD unf_signal = *result.unfolded_signal_;
    unf_signal *= events_to_xsec_factor;
    dump_text_column_vector("dump/vec_table_unfolded_signal.txt", unf_signal);

    // Dump similar tables for each of the theoretical predictions (and the fake
    // data truth if applicable). Note that this function expects that the
    // additional smearing matrix A_C has not been applied to these predictions.
    TMatrixD temp_genie_cv = genie_cv_true_events;
    temp_genie_cv *= events_to_xsec_factor;
    dump_text_column_vector("dump/vec_table_uBTune.txt", temp_genie_cv);

    if (using_fake_data)
    {
        TMatrixD temp_fake_truth = fake_data_true_events;
        temp_fake_truth *= events_to_xsec_factor;
        dump_text_column_vector("dump/vec_table_FakeData.txt", temp_fake_truth);
    }

    for (const auto &gen_pair : generator_truth_map)
    {
        std::string gen_short_name = samples_to_hist_names.at(gen_pair.first);
        TMatrixD temp_gen = *gen_pair.second;
        temp_gen *= events_to_xsec_factor;
        dump_text_column_vector("dump/vec_table_" + gen_short_name + ".txt",
                                temp_gen);
    }

    // No unit conversions are necessary for the unfolding, error propagation,
    // and additional smearing matrices since they are dimensionless
    dump_text_matrix("dump/mat_table_unfolding.txt", *result.unfolding_matrix_);
    dump_text_matrix("dump/mat_table_err_prop.txt", *result.err_prop_matrix_);
    dump_text_matrix("dump/mat_table_add_smear.txt", *result.add_smear_matrix_);

    // Convert units on the covariance matrices one-by-one and dump them
    for (const auto &cov_pair : unf_cov_matrix_map)
    {
        const auto &name = cov_pair.first;
        TMatrixD temp_cov_matrix = *cov_pair.second;
        // Note that we need to square the unit conversion factor for the
        // covariance matrix elements
        temp_cov_matrix *= std::pow(events_to_xsec_factor, 2);
        dump_text_matrix("dump/mat_table_cov_" + name + ".txt", temp_cov_matrix);
    }

    // Finally, dump a summary table of the flux-averaged total cross section
    // measurements and their statistical and total uncertainties
    TMatrixD temp_stat_cov = *unf_cov_matrix_map.at("DataStats");
    TMatrixD temp_total_cov = *unf_cov_matrix_map.at("total");
    temp_stat_cov *= std::pow(events_to_xsec_factor, 2);
    temp_total_cov *= std::pow(events_to_xsec_factor, 2);

    // Open the output file and set up the output stream so that full numerical
    // precision is preserved in the ascii text representation
    std::ofstream out_summary_file("dump/xsec_summary_table.txt");
    out_summary_file << std::scientific
                     << std::setprecision(std::numeric_limits<double>::max_digits10);

    int num_bins = unf_signal.GetNrows();
    out_summary_file << "numXbins " << num_bins;

    for (int bin = 0; bin < num_bins; ++bin)
    {
        double xsec = unf_signal(bin, 0);
        double stat_err = std::sqrt(std::max(0., temp_stat_cov(bin, bin)));
        double total_err = std::sqrt(std::max(0., temp_total_cov(bin, bin)));
        out_summary_file << '\n'
                         << bin << "  " << xsec << "  " << stat_err
                         << "  " << total_err;
    }
}

TH1D* get_generator_hist(const TString& filePath, const unsigned int sl_idx, const float scaling = 1.f )
{
    // Open the file
    TFile* file = new TFile(filePath, "readonly");

    // Check if the file was successfully opened
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return nullptr;
    }

    // Define the plot names
    std::vector<TString> plotNames = {
        "TrueMuonCosThetaPlot",
        "TrueMuonPhiPlot",
        "TrueMuonMomentumPlot",
        "TruePionCosThetaPlot",
        "TruePionPhiPlot",
        "TruePionMomentumPlot",
        "TrueMuonPionOpeningAnglePlot",
        "TrueTotalPlot"
    };

    // Check if the index is valid
    if (sl_idx >= plotNames.size()) {
        std::cerr << "Invalid slice index: " << sl_idx << std::endl;
        return nullptr;
    }

    // Get the histogram from the file
    TH1D* hist = (TH1D*)file->Get(plotNames[sl_idx]);

    // Check if the histogram was successfully retrieved
    if (!hist) {
        std::cerr << "Failed to retrieve histogram: " << plotNames[sl_idx] << std::endl;
        return nullptr;
    }

    // Scale the histogram
    hist->Scale(scaling);

    // Apply the scaling and selection
    // hist->SetLineWidth(4);
    // hist->GetXaxis()->SetNdivisions(8);
    // hist->GetYaxis()->SetNdivisions(6);
    // hist->GetYaxis()->SetTitle("Cross Section [10^{-38} cm^{2}/Ar]");

    return hist;
}

struct inputFiles
{
  std::string rootFile;
  std::string fileList;
  std::string config;
  std::string sliceConfig;
  std::string nameExtension;
};

void unfold_nuwro()
{
    // inputFiles input{ // All nuwro cross-sections with only contained muons for the muon momentum cross-section
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_nuwro_run1234bcd5_09Apr24_testingOnly_lowPiMomThreshold_containedMuon.root",
    //     "../nuwro_file_properties_testingOnly_lowPiMomThreshold.txt",
    //     "../systcalc_unfold_fd_min.conf",
    //     "../ubcc1pi_neutral_slice_config.txt", // <-- This version includes the _lowPiMomThreshold underflow bin removal
    //     "_fd_testingOnly_lowPiMomThreshold_containedMuon_withGenerators"
    // };

    // inputFiles input{ // All nuwro cross-sections with only contained muons for the muon momentum cross-section and modified bin size (muon momentum ending at 1.2GeV rather than 1.5)
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_nuwro_run1234bcd5_15Apr24_testingOnly_lowPiMomThreshold_containedMuon_muon1200MeV.root",
    //     "../nuwro_file_properties_testingOnly_lowPiMomThreshold.txt",
    //     "../systcalc_unfold_fd_min.conf",
    //     "../ubcc1pi_neutral_slice_config_muon1200MeV.txt", // <-- This verion includes the _lowPiMomThreshold underflow bin removal & muon momentum endind at 1.2GeV rather than 1.5
    //     "_fd_testingOnly_lowPiMomThreshold_containedMuon_muon1200MeV"
    // };

    // inputFiles input{ // All nuwro cross-sections
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_nuwro_run1234bcd5_28Mar24_testingOnly_lowPiMomThreshold.root",
    //     "../nuwro_file_properties_testingOnly_lowPiMomThreshold.txt",
    //     "../systcalc_unfold_fd_min.conf",
    //     "../ubcc1pi_neutral_slice_config.txt", // <-- This version includes the _lowPiMomThreshold underflow bin removal
    //     "_fd_testingOnly_lowPiMomThreshold_withGenerators"
    // };

    // inputFiles input{ // All nuwro cross-sections; removed the nProton cross-section
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_nuwro_run1234bcd5_07May24_testingOnly_lowPiMomThreshold.root",
    //     "../nuwro_file_properties_testingOnly_lowPiMomThreshold.txt",
    //     "../systcalc_unfold_fd_min.conf",
    //     "../ubcc1pi_neutral_slice_config.txt", // <-- This version includes the _lowPiMomThreshold underflow bin removal and removes the nProton cross-section
    //     "_fd_testingOnly_lowPiMomThreshold_withGenerators_08May24"
    // };

    // inputFiles input{ // All nuwro cross-sections with full uncertainties like real data
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_nuwro_run1234bcd5_29May24_testingOnly_lowPiMomThreshold_allUncertainties.root",
    //     "../nuwro_file_properties_testingOnly_lowPiMomThreshold_allUncertainties.txt",
    //     "../systcalc_unfold.conf",
    //     "../ubcc1pi_neutral_slice_config.txt", // <-- This version includes the _lowPiMomThreshold underflow bin removal and remove nProton cross-section
    //     "_fd_testingOnly_lowPiMomThreshold_withGenerators_30May24_allUncertainties"
    // };

    // ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    // #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
    // ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    // ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    // # # # # # # # # # # # # # # # #

    // ##########################################################
    // Remember to consider how to handle the overflow in truth 
    // ##########################################################
    // std::cout<<"\n\n\n\n\nWarning: The overflow bins are excluded both truth and reco here\n\n\n\n\n"<<std::endl;
    // inputFiles input{ // The bin definition used for the reco file ensures that CC1pi events with p_pi < 100MeV are not double counted + overflow bins have been removed from the bin and slice definitions
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_nuwro_run1234bcd5_28Jun24_testingOnly_lowPiMomThreshold_allUncertainties_fixedBackground_noOverflow.root",
    //     "../nuwro_file_properties_testingOnly_lowPiMomThreshold_allUncertainties.txt",
    //     "../systcalc_unfold.conf",
    //     "../ubcc1pi_neutral_slice_config_noOverflow.txt",
    //     "_nuwro_fixedBackground_noOverflow"
    // };

    // inputFiles input{ // The bin definition used for the reco file ensures that CC1pi events with p_pi < 100MeV are not double counted
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_nuwro_run1234bcd5_26Jun24_testingOnly_lowPiMomThreshold_allUncertainties_fixedBackground.root",
    //     "../nuwro_file_properties_testingOnly_lowPiMomThreshold_allUncertainties.txt",
    //     "../systcalc_unfold_fd_min.conf",
    //     "../ubcc1pi_neutral_slice_config.txt",
    //     "_nuwro_fixedBackground"
    // };

    inputFiles input{ // The bin definition used for the reco file ensures that CC1pi events with p_pi < 100MeV are not double counted + overflow bins have been removed from the bin and slice definitions
        "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_nuwro_run1234bcd5_02Jul24_testingOnly_lowPiMomThreshold_allUncertainties_fixedBackground_mergedOverflow.root",
        "../nuwro_file_properties_testingOnly_lowPiMomThreshold_allUncertainties.txt",
        "../systcalc_unfold_fd_min.conf",
        "../ubcc1pi_neutral_slice_config_mergedOverflow.txt",
        "_nuwro_fixedBackground_mergedOverflow"
    };

    // Initialize the FilePropertiesManager and tell it to treat the NuWro MC ntuples as if they were data
    auto &fpm = FilePropertiesManager::Instance();

    // fpm.load_file_properties("../nuwro_file_properties.txt");
    // fpm.load_file_properties("../nuwro_file_properties_testingOnly.txt");
    // fpm.load_file_properties("../nuwro_file_properties_testingOnly_lowPiMomThreshold.txt");
    fpm.load_file_properties( input.fileList );

    // const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_10/univmake_output_nuwro_with_sideband_noOverflow_13Jan23.root");
    // const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_10/univmake_output_nuwro_with_sideband_noOverflow_noGolden_13Jan23.root");
    // const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_10/univmake_output_nuwro_with_sideband_overflow_all_13Jan23.root");
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_run1234bcd5_6Mar24.root"); // <-- Yes the name is wrong and should say nuwro
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_run1234bcd5_14Mar24_testingOnly.root"); // <-- Yes the name is wrong and should say nuwro
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_nuwro_run1234bcd5_28Mar24_testingOnly_lowPiMomThreshold.root");
    const std::string respmat_file_name( input.rootFile );

    // const std::string postfix = "_run1_reduced";
    // const std::string postfix = "_with_Overflow_onlyGolden";
    // const std::string postfix = "_fd_testingOnly";
    // const std::string postfix = "_fd_testingOnly_lowPiMomThreshold";
    const std::string postfix = input.nameExtension;

    // Plot slices of the unfolded result
    // auto *sb_ptr = new SliceBinning("../ubcc1pi_neutral_slice_config_withPionUnderflow.txt");
    // auto *sb_ptr = new SliceBinning("../ubcc1pi_neutral_slice_config.txt");
    auto *sb_ptr = new SliceBinning( input.sliceConfig );
    auto &sb = *sb_ptr;

    // Do the systematics calculations in preparation for unfolding
    // const auto *syst_ptr = new MCC9SystematicsCalculator(respmat_file_name, "../systcalc_unfold_fd_min.conf");
    const auto *syst_ptr = new MCC9SystematicsCalculator(respmat_file_name, input.config);
    const auto &syst = *syst_ptr;

    const std::string genPath = "/exp/uboone/app/users/jdetje/BuildEventGenerators/FlatTreeAnalyzer/OutputFiles/";
    // Format: path, lineColor, lineStyle, lineWidth, name, scaling
    constexpr double fluxIntNuMu = 7.3762291e-10; // This is technically flux/POT
	constexpr double fluxIntNuMuBar = 4.5589190e-11;
    std::vector<GeneratorInfo> generators = {
        // {genPath + "FlatTreeAnalyzerOutput_GENIE.root", kBlue+2, 0, 1, "GENIE NUMU 2.12.10", 1.f}, // Old numu only sample from wiki
        // {genPath + "FlatTreeAnalyzerOutput_NEUT.root", kRed+1, 0, 1, "NEUT NUMU 5.4.0.1", 1.f}, // Old numu only sample from wiki
        // {genPath + "FlatTreeAnalyzerOutput_GiBUU.root", kOrange+7, 0, 1, "GiBUU NUMU ?3.0.6?", 100.f}, // Old numu only sample from wiki
        // {genPath + "FlatTreeAnalyzerOutput_NuWro.root", kGreen+1, 0, 1, "NuWro NUMU 19.02.1", 1.f}, // Old numu only sample from wiki

        // {genPath + "FlatTreeAnalyzerOutput_GENIE_NuMu_NuMuBar.root", kBlue+1, 1, 1, "GENIE 3.04.02", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_NuWro_NuMu_NuMuBar.root", kGreen+3, 1, 1, "NuWro 21.09.2", 1.f}, // <----------- Use these
        // {genPath + "FlatTreeAnalyzerOutput_NuWro_NuMu_NuMuBar_19_02_1.root", kGreen, 1, 1, "NuWro 19.02.1", 1.f},

        // {genPath + "FlatTreeAnalyzerOutput_NuWro_numu_numubar_19_02_1_events_2000000_each.root", kRed, 1, 1, "NuWro 19.02.1 equal numu/bar", 1.f}, // Doesn't change anything

        // {genPath + "FlatTreeAnalyzerOutput_NuWro_NuMu_19_02_1.root", kRed+1, 1, 1, "NuWro NUMU with Scaling 19.02.1", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_GENIE_NuMu.root", kRed+2, 1, 1, "GENIE NUMU with Scaling 3.04.02", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_GENIE_NuMu_noScaling.root", kBlue+2, 2, 1, "GENIE NUMU NoScaling", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_NuWro_NuMu_noScaling.root", kGreen+2, 2, 1, "NuWro NUMU NoScaling", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_NuWro_NuMu_19_02_1_noScaling.root", kGreen+7, 2, 1, "NuWro NUMU 19.02.1 NoScaling", 1.f},

        // {genPath + "FlatTreeAnalyzerOutput_Genie_numu_numubar_CC_v3_4_2_AR23_20i_00_000.root", kBlue+4, 1, 1, "GENIE 3.04.02 DUNE/SBN Tune", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_Genie_numu_numubar_CC_v3_4_2_G18_10a_02_11a.root", kBlue+2, 1, 1, "GENIE 3.04.02 G18_10a_02_11a", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_Genie_numu_CC_v3_4_2_G18_10a_02_11a.root", kBlue, 1, 1, "GENIE 3.04.02 G18_10a_02_11a NuMu", (fluxIntNuMu + fluxIntNuMuBar)/fluxIntNuMu},
        // {genPath + "FlatTreeAnalyzerOutput_Genie_numubar_CC_v3_4_2_G18_10a_02_11a.root", kBlue+1, 1, 1, "GENIE 3.04.02 G18_10a_02_11a NuMuBar", (fluxIntNuMu + fluxIntNuMuBar)/fluxIntNuMuBar},
        // {genPath + "FlatTreeAnalyzerOutput_Genie_numu_CC_v3_4_2_G18_10a_02_11a.root", kBlue+3, 1, 1, "GENIE 3.04.02 G18_10a_02_11a NuMu Scaled", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_Genie_numubar_CC_v3_4_2_G18_10a_02_11a.root", kBlue+4, 1, 1, "GENIE 3.04.02 G18_10a_02_11a NuMuBar Scaled", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_NuWro_numu_numubar_CC_v21_09_2.root", kGreen+3, 1, 1, "NuWro 21.09.2", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_NuWro_numu_numubar_CC_v19_02_1.root", kRed, 1, 1, "NuWro 19.02.1", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_NuWro_numu_CC_v19_02_1.root", kRed + 1, 1, 1, "NuWro 19.02.1 NuMu", (fluxIntNuMu + fluxIntNuMuBar)/fluxIntNuMu},
        // {genPath + "FlatTreeAnalyzerOutput_NuWro_numubar_CC_v19_02_1.root", kRed + 3, 1, 1, "NuWro 19.02.1 NuMuBar", (fluxIntNuMu + fluxIntNuMuBar)/fluxIntNuMuBar},


        {genPath + "FlatTreeAnalyzerOutput_Genie_numu_numubar_CC_v3_4_2_G18_10a_02_11a.root", kAzure+1, 1, 1, "#splitline{GENIE 3.04.02}{G18_10a_02_11a}", 1.f},
        {genPath + "FlatTreeAnalyzerOutput_Genie_numu_numubar_CC_v3_4_2_ModAR23_20i_00_000.root", kAzure+3, 1, 1, "#splitline{GENIE 3.04.02}{#splitline{AR23_20i_00_000}{Nuclear deexcitation off}}", 1.f},
        {genPath + "FlatTreeAnalyzerOutput_NuWro_numu_numubar_CC_v21_09_2.root", kTeal-7, 1, 1, "NuWro 21.09.2", 1.f},
        {genPath + "FlatTreeAnalyzerOutput_GiBUU2023_patch3_numu_150_numubar_15.root", kOrange, 1, 1, "GiBUU 2023 Patch 3", 1.f},
        // {genPath + "FlatTreeAnalyzerOutput_NuWro_numu_numubar_CC_v19_02_1.root", kTeal+6, 1, 1, "NuWro 19.02.1", 1.f},
        {genPath + "FlatTreeAnalyzerOutput_NuWro_numu_numubar_CC_v19_02_1_cthorpe.root", kGreen+4, 1, 1, "NuWro (v19 + Hyperon Models*)", 1.f},
    };

    // Get the tuned GENIE CV prediction in each true bin (including the
    // background true bins)
    TH1D *genie_cv_truth = syst.cv_universe().hist_true_.get();
    int num_true_bins = genie_cv_truth->GetNbinsX();

    // While we're at it, clone the histogram and zero it out. We'll fill this
    // one with our unfolded result for easy comparison
    TH1D *unfolded_events = dynamic_cast<TH1D *>(genie_cv_truth->Clone("unfolded_events"));
    unfolded_events->Reset();

    auto getBinCount = [](auto &bins, auto type, auto condition)
    {
        return std::count_if(bins.begin(), bins.end(), [&](const auto &bin)
                             { return bin.type_ == type == condition; });
    };

    // If present, then get the fake data event counts in each true bin
    // (including the background true bins). We hope to approximately reproduce
    // these event counts in the signal true bins via unfolding the fake data.
    const auto &fake_data_univ = syst.fake_data_universe();
    TH1D *fake_data_truth_hist = fake_data_univ ? fake_data_univ->hist_true_.get() : nullptr;
    const auto using_fake_data = (fake_data_univ != nullptr);
    std::cout<<"DEBUG using_fake_data: "<<using_fake_data<<std::endl;

    auto num_ordinary_reco_bins = getBinCount(syst.reco_bins_, kSidebandRecoBin, false);
    auto num_sideband_reco_bins = syst.reco_bins_.size() - num_ordinary_reco_bins;
    auto num_true_signal_bins = getBinCount(syst.true_bins_, kSignalTrueBin, true);

    std::cout << "NUM ORDINARY RECO BINS = " << num_ordinary_reco_bins << '\n';
    std::cout << "NUM TRUE SIGNAL BINS = " << num_true_signal_bins << '\n';

    auto *matrix_map_ptr = syst.get_covariances().release();
    auto &matrix_map = *matrix_map_ptr;
    auto *cov_mat = matrix_map.at("total").cov_matrix_.get();


    // // Define your custom labels and intervals
    // const Int_t n = 8;
    // const std::vector<int> bins{0, 11, 26, 32, 39, 49, 54, 61, 62};
    // const Int_t overUnderFlowN = 2;
    // const std::vector<int> overUnderFlowBins{31, 53};
    // const Char_t *labels[n] = {"cos(#theta_{#mu})", "#phi_{#mu}", "p_{#mu}", "cos(#theta_{#pi})", "#phi_{#pi}", "p_{#pi}^{**}", "#theta_{#pi #mu}", "Total"};

    TMatrixD* corr = new TMatrixD(util::CovarianceMatrixToCorrelationMatrix(util::TH2DToTMatrixD(*cov_mat)));
    TH2D* corr_hist = new TH2D(util::TMatrixDToTH2D(*corr, "corr", "Correlation Matrix", 0, corr->GetNrows(), 0, corr->GetNcols()));
    delete corr; // Delete the dynamically allocated memory
    UnfolderHelper::plot_entire_matrix(corr_hist, "Correlation Matrix", "Reco Bins", "Reco Bins", "corr"+postfix);

    // std::unique_ptr<Unfolder> unfolder(new DAgostiniUnfolder(DAgostiniUnfolder::ConvergenceCriterion::FigureOfMerit, 0.025));
    std::unique_ptr<Unfolder> unfolder(new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kFirstDeriv));
    // std::unique_ptr<Unfolder> unfolder(new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv));
    // constexpr int NUM_DAGOSTINI_ITERATIONS = 50;
    // std::unique_ptr<Unfolder> unfolder(new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS ));
    constexpr bool USE_ADD_SMEAR = true;

    // // Extract the inputs needed for the unfolding procedure from the
    // // supplied SystematicsCalculator object
    // const auto &smearcept = syst.get_cv_smearceptance_matrix();
    // const auto &true_signal = syst.get_cv_true_signal();

    // // const auto &reco_selected = syst.get_cv_reco_selected();
    // // const auto data_signal = syst.get_cv_reco_signal();
    // const auto &meas = syst.get_measured_events();
    // const auto &data_signal = meas.reco_signal_;
    // const auto &data_covmat = meas.cov_matrix_;

    const auto result = unfolder->unfold(syst);

    TH2D* error_prop_hist = new TH2D(util::TMatrixDToTH2D(*result.err_prop_matrix_, "err_prop_matrix", "Error Propagation Matrix", 0, result.err_prop_matrix_->GetNrows(), 0, result.err_prop_matrix_->GetNcols()));
    UnfolderHelper::plot_entire_matrix(error_prop_hist, "Error Propagation Matrix", "Truth Bins", "Reco Bins", "err_prop"+postfix, false, -5, 10);//, true, -10, 10);
    delete error_prop_hist;

    // Propagate all defined covariance matrices through the unfolding procedure
    TMatrixD err_prop_tr(TMatrixD::kTransposed, *result.err_prop_matrix_);
    std::map<std::string, std::shared_ptr<TMatrixD>> unfolded_cov_matrix_map;
    std::map<std::string, std::shared_ptr<TMatrixD>> cov_TMatrix_map;
    for (const auto &matrix_pair : matrix_map)
    {
        unfolded_cov_matrix_map[matrix_pair.first] = std::make_unique<TMatrixD>(
            *result.err_prop_matrix_, TMatrixD::EMatrixCreatorsOp2::kMult,
            TMatrixD(*matrix_pair.second.get_matrix(), TMatrixD::EMatrixCreatorsOp2::kMult, err_prop_tr));

        cov_TMatrix_map[matrix_pair.first] = std::make_unique<TMatrixD>(*matrix_pair.second.get_matrix());
    }

    // Decompose the block-diagonal pieces of the total covariance matrix
    // into normalization, shape, and mixed components (for later plotting
    // purposes)
    NormShapeCovMatrix bd_ns_covmat = make_block_diagonal_norm_shape_covmat(
        *result.unfolded_signal_, *result.cov_matrix_, syst.true_bins_);

    // // Assuming result_cov_matrix is a std::unique_ptr<TMatrixD>
    // int nRows = result.cov_matrix_->GetNrows();
    // int nCols = result.cov_matrix_->GetNcols();

    // // Create a new TH2D object
    // TH2D *h2 = new TH2D("h2", "Covariance Matrix", nCols, 0, nCols, nRows, 0, nRows);

    // // Fill the TH2D object with the values from the TMatrixD
    // for (int i = 1; i <= nRows; ++i) {
    //     for (int j = 1; j <= nCols; ++j) {
    //         h2->SetBinContent(j, i, (*result.cov_matrix_)(i-1, j-1));
    //     }
    // }

    // Draw the TH2D object on the canvas
    TCanvas *cm1 = new TCanvas("cm1", "Canvas", 800, 600);

    TMatrixD* corr_unfolded = new TMatrixD(util::CovarianceMatrixToCorrelationMatrix(*result.cov_matrix_));
    TH2D* corr_unfolded_hist = new TH2D(util::TMatrixDToTH2D(*corr_unfolded, "corr_unfolded", "Correlation Matrix", 0, corr_unfolded->GetNrows(), 0, corr_unfolded->GetNcols()));
    delete corr_unfolded; // Delete the dynamically allocated memory
    UnfolderHelper::plot_entire_matrix(corr_unfolded_hist, "Correlation matrix after blockwise unfolding", "Truth Bins", "Truth Bins", "unfolded_corr"+postfix);
        TH2D* cov_unfolded_hist = new TH2D(util::TMatrixDToTH2D(*result.cov_matrix_, "cov_unfolded", "Covariance Matrix", 0, result.cov_matrix_->GetNrows(), 0, result.cov_matrix_->GetNcols()));

    // Add the blockwise decomposed matrices into the map
    unfolded_cov_matrix_map["total_blockwise_norm"] = std::make_unique<TMatrixD>(bd_ns_covmat.norm_);
    unfolded_cov_matrix_map["total_blockwise_shape"] = std::make_unique<TMatrixD>(bd_ns_covmat.shape_);
    unfolded_cov_matrix_map["total_blockwise_mixed"] = std::make_unique<TMatrixD>(bd_ns_covmat.mixed_);

    // Set the event counts in each bin of the histogram that displays the
    // unfolded result. Note that we don't care about the background true bins
    // (which are assumed to follow all of the signal true bins) since we've
    // subtracted out an estimate of the background before unfolding.
    for (int t = 0; t < num_true_bins; ++t)
    {
        const auto evts = (t < num_true_signal_bins) ? result.unfolded_signal_->operator()(t, 0) : 0.;
        const auto error = (t < num_true_signal_bins) ? std::sqrt(std::max(0., result.cov_matrix_->operator()(t, t))) : 0.;

        // We need to use one-based indices while working with TH1D bins
        unfolded_events->SetBinContent(t + 1, evts);
        unfolded_events->SetBinError(t + 1, error);
    }

    unfolded_events->SetStats(false);
    unfolded_events->SetLineColor(kBlack);
    unfolded_events->SetLineWidth(5);
    unfolded_events->GetXaxis()->SetRangeUser(0, num_true_signal_bins);

    // Save the fake data truth (before A_C multiplication) using a column vector
    // of event counts
    TMatrixD fake_data_truth(num_true_signal_bins, 1);
    TMatrixD fake_data_truth_cov(num_true_signal_bins, num_true_signal_bins);
    if (using_fake_data)
    {
        for (int b = 0; b < num_true_signal_bins; ++b)
        {
            const auto true_evts = fake_data_truth_hist->GetBinContent(b + 1);
            const auto true_evts_err = fake_data_truth_hist->GetBinError(b + 1);
            fake_data_truth(b, 0) = true_evts;
            fake_data_truth_cov(b, b) = true_evts_err * true_evts_err;
        }
    }

    // Save the GENIE CV model (before A_C multiplication) using a column vector
    // of event counts
    TMatrixD genie_cv_truth_vec(num_true_signal_bins, 1);
    TMatrixD genie_cv_truth_err_vec(num_true_signal_bins, num_true_signal_bins);
    for (int b = 0; b < num_true_signal_bins; ++b)
    {
        const auto true_evts = genie_cv_truth->GetBinContent(b + 1);
        const auto true_evts_err = genie_cv_truth->GetBinError(b + 1);
        genie_cv_truth_vec(b, 0) = true_evts;
        genie_cv_truth_err_vec(b, b) = true_evts_err * true_evts_err;
    }

    // Multiply the truth-level GENIE prediction histogram by the additional
    // smearing matrix
    TMatrixD *A_C = result.add_smear_matrix_.get();

    // Define the colors and stops for the gradient
    constexpr int nColors = 4;
    Double_t stops[nColors] = {0.0, 1/4.0, 3/4.0, 1.0};
    Double_t red[nColors] = {1.0, 1.0, 0.0, 0.0};
    Double_t green[nColors] = {0.0, 1.0, 1.0, 0.0};
    Double_t blue[nColors] = {0.0, 1.0, 1.0, 1.0};

    // Create the gradient color table
    // const int paletteNumber = 
    TColor::CreateGradientColorTable(nColors, stops, red, green, blue, 20);
    // gStyle->SetPalette(paletteNumber);

    // Convert TMatrixD to TH2D using the TMatrixDToTH2D function
    const TH2D* h_A_C = new TH2D(util::TMatrixDToTH2D(*A_C, "h_A_C", "Additional smearing matrix", 0, A_C->GetNcols(), 0, A_C->GetNrows()));
    UnfolderHelper::plot_entire_matrix(h_A_C, "Additional smearing matrix", "Truth Bins", "Regularized Truth Bins", "ac"+postfix, true, -0.5, 1.5);

    if(USE_ADD_SMEAR) multiply_1d_hist_by_matrix(A_C, genie_cv_truth); // Can be directly multiplied for all slices because A_C is a block diagonal matrix 

    genie_cv_truth->SetStats(false);
    genie_cv_truth->SetLineColor(kRed);
    genie_cv_truth->SetLineWidth(3);
    genie_cv_truth->SetLineStyle(9);

    unfolded_events->Draw("E1");
    genie_cv_truth->Draw("hist same");

    if (using_fake_data)
    {

        // Multiply the fake data truth histogram by the additional smearing matrix
        if(USE_ADD_SMEAR) multiply_1d_hist_by_matrix(A_C, fake_data_truth_hist);

        fake_data_truth_hist->SetStats(false);
        fake_data_truth_hist->SetLineColor(kBlue);
        fake_data_truth_hist->SetLineWidth(1);
        fake_data_truth_hist->SetLineStyle(2);
        fake_data_truth_hist->Draw("hist same");
    }

    // TLegend* lg = new TLegend( 0.15, 0.7, 0.3, 0.85 ); // Top right
    // TLegend* lg = new TLegend( 0.7, 0.7, 0.85, 0.85 ); // Top left
    // TLegend* lg = new TLegend( 0.7, 0.15, 0.85, 0.3 ); // Bottom left
    // TLegend* lg = new TLegend( 0.15, 0.15, 0.3, 0.3 ); // Bottom right
    TLegend *lg = new TLegend(0.6, 0.2);
    lg->AddEntry(unfolded_events, "unfolded", "l");
    lg->AddEntry(genie_cv_truth, "uB tune", "l");
    if (using_fake_data)
    {
        lg->AddEntry(fake_data_truth_hist, "Fake Data Truth (NuWro 19.02.1)", "l");
    }

    lg->Draw("same");


    // Get the factors needed to convert to cross-section units
    double total_pot = syst.total_bnb_data_pot_;
    double integ_flux = integrated_numu_flux_in_FV(total_pot);
    double num_Ar = num_Ar_targets_in_FV();

    std::cout << "TOTAL POT = " << total_pot << '\n';
    std::cout << "INTEGRATED numu FLUX = " << integ_flux << '\n';
    std::cout << "NUM Ar atoms in fiducial volume = " << num_Ar << '\n';

    // Retrieve the true-space expected event counts from NUISANCE output files
    // for each available generator model
    double conv_factor = (num_Ar * integ_flux) / 1e39;
    std::map<std::string, TMatrixD *> generator_truth_map = {}; // get_true_events_nuisance( sample_info, conv_factor );

    // Dump overall results to text files. Total cross section units (10^{-39}
    // cm^2 / Ar) will be used throughout. Do this before adjusting the
    // truth-level prediction TMatrixD objects via multiplication by A_C
    dump_overall_results(result, unfolded_cov_matrix_map, 1.0 / conv_factor, genie_cv_truth_vec, fake_data_truth, generator_truth_map, using_fake_data);

    if (USE_ADD_SMEAR)
    {
        // Get access to the additional smearing matrix
        const TMatrixD &A_C = *result.add_smear_matrix_;
        const TMatrixD A_C_T = TMatrixD(TMatrixD::kTransposed, A_C);

        // Start with the fake data truth if present
        if (using_fake_data)
        {
            fake_data_truth = TMatrixD(A_C, TMatrixD::kMult, fake_data_truth);
            fake_data_truth_cov = TMatrixD(TMatrixD(A_C, TMatrixD::kMult, fake_data_truth_cov), TMatrixD::kMult, A_C_T); // A_c * cov * A_c^T
        }

        // Also transform the GENIE CV model
        genie_cv_truth_vec = TMatrixD(A_C, TMatrixD::kMult, genie_cv_truth_vec);
        genie_cv_truth_err_vec = TMatrixD(A_C, TMatrixD::kMult, genie_cv_truth_err_vec);

        // Now do the other generator predictions
        for (const auto &pair : generator_truth_map)
            *pair.second = TMatrixD(A_C, TMatrixD::kMult, *pair.second);
    }

    const auto cv_hist_2d = syst.get_cv_hist_2d();

    // ###########################################################
    // ###      ####  ########  ####      ####      ####      ######
    // ###  ########  ########  ####  ########  ########  ############
    // ###      ####  ########  ####  ########      ####      ##########
    // #######  ####  ########  ####  ########  ############  ########
    // ###      ####      ####  ####      ####      ####      ######
    // ###########################################################
    for (size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx)
    {
        const auto &slice = sb.slices_.at(sl_idx);

        // Make a histogram showing the unfolded true event counts in the current slice
        SliceHistogram *slice_unf = SliceHistogram::make_slice_histogram(
            *result.unfolded_signal_, slice, result.cov_matrix_.get());

        // Make a histogram showing the unfolded true event counts in the current slice using only the DataStats covariance matrix
        SliceHistogram *slice_unf_stats_only = SliceHistogram::make_slice_histogram(
                    *result.unfolded_signal_, slice, unfolded_cov_matrix_map.at("DataStats").get());


        auto meas = syst.get_measured_events();
        const auto& data_signal = meas.reco_signal_;
        const auto cov_matrix_map = syst.get_covariances();


        // Make unfolded fractional uncertainty plots
        std::vector< std::string > cov_mat_keys_total = { "total", "xsec_total", "MCstats", "BNBstats"};
        std::vector< std::string > cov_mat_keys_xsec = { "xsec_total", "xsec_multi", "xsec_unisim", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie"};

        // TMatrixD* bkgnd_subtracted_reco_data_hist = new TMatrixD(meas.reco_signal_);
        UnfolderHelper::fractional_uncertainty_plot(slice, cov_TMatrix_map, data_signal.get(), sl_idx, cov_mat_keys_total, "total", postfix + "_reco_total");
        UnfolderHelper::fractional_uncertainty_plot(slice, cov_TMatrix_map, data_signal.get(), sl_idx, cov_mat_keys_xsec, "xsec_total", postfix + "_reco_xsec");
        // delete bkgnd_subtracted_reco_data_hist;

        // This assumes a continuous range of individual bins
        // int t_start = slice.bin_map_.begin()->first;
        // int t_stop = slice.bin_map_.rbegin()->first;
        size_t start = std::numeric_limits<size_t>::max();
        size_t stop = 0;
        for (const auto& entry : slice.bin_map_) {
            const auto& set = entry.second;
            if (set.size() != 1) {
                throw std::runtime_error("Error: set in bin_map_ has more or less than 1 entry");
            }
            start = std::min(start, *set.begin());
            stop = std::max(stop, *set.rbegin());
        }

        TMatrixD cv_hist_2d_slice(stop - start + 1, stop - start + 1);
        TMatrixD fake_hist_2d_slice(stop - start + 1, stop - start + 1);
        TMatrixD corr_unfolded_hist_slice(stop - start + 1, stop - start + 1);
        TMatrixD cov_unfolded_hist_slice(stop - start + 1, stop - start + 1);
        TMatrixD corr_hist_slice(stop - start + 1, stop - start + 1);
        TMatrixD cov_hist_slice(stop - start + 1, stop - start + 1);
        for (int i = start; i <= stop; i++) {
            for (int j = start; j <= stop; j++) {
                cv_hist_2d_slice(i - start, j - start) = cv_hist_2d->GetBinContent(i+1, j+1);
                if(using_fake_data) fake_hist_2d_slice(i - start, j - start) = fake_data_univ->hist_2d_->GetBinContent(i+1, j+1);
                corr_unfolded_hist_slice(i - start, j - start) = corr_unfolded_hist->GetBinContent(i+1, j+1);
                cov_unfolded_hist_slice(i - start, j - start) = cov_unfolded_hist->GetBinContent(i+1, j+1);
                corr_hist_slice(i - start, j - start) = corr_hist->GetBinContent(i+1, j+1);
                cov_hist_slice(i - start, j - start) = cov_mat->GetBinContent(i+1, j+1);
            }
        }

        // Add over and underflow bins to additional smearing matrix plots
        if(sl_idx==2 || sl_idx==5) stop++; // Overflow bin
        // if(sl_idx==5) start--; // Underflow bin // Removed and replaced withh 100 MeV phase space restriction
        TMatrixD ac_hist_slice(stop - start + 1, stop - start + 1);
        for (int i = start; i <= stop; i++) {
            for (int j = start; j <= stop; j++) {
                ac_hist_slice(i - start, j - start) = A_C->operator()(i, j);
            }
        }

        const auto cv_confusion_mat = util::CountsToConfusionMatrix(cv_hist_2d_slice, "row");
        const auto fake_confusion_mat = util::CountsToConfusionMatrix(fake_hist_2d_slice, "row");
        std::string title = slice_unf->hist_->GetXaxis()->GetTitle();
        // std::regex pattern("true (.*)");
        // title = std::regex_replace(title, pattern, "$1");
        const auto conf_title = std::string("Confusion Matrix") + title +";Truth Bins;Reco Bins";
        const auto ac_title = std::string("Additional Smearing Matrix") + title +";Truth Bins;Regularized Truth Bins";
        const auto corr_unf_title = std::string("Correlation Matrix for Unfolded Bins") + title +";Truth Bins;Truth Bins";
        const auto cov_unf_title = std::string("Covariance Matrix for Unfolded Bins") + title +";Truth Bins;Truth Bins";
        const auto corr_title = std::string("Correlation Matrix for Reco Bins") + title +";Reco Bins;Reco Bins";
        const auto cov_title = std::string("Covariance Matrix for Reco Bins") + title +";Reco Bins;Reco Bins";

        const int num_bins = slice.hist_->GetNbinsX();
        const int num_bins_ac = stop - start + 1;
        // Double_t edges[num_bins+1];
        // for (int i = 0; i < num_bins; i++) {
        //     edges[i] = slice.hist_->GetBinLowEdge(i+1);
        // }
        // edges[num_bins] = slice.hist_->GetBinLowEdge(num_bins) + slice.hist_->GetBinWidth(num_bins);

        // Create a TH2D histogram
        TH2D *cv_confusion_hist = new TH2D(("cv_confusion slice "+std::to_string(sl_idx)).c_str(),("Genie " + conf_title).c_str(),
                                num_bins, 0, num_bins,
                                num_bins, 0, num_bins);

        // Create a TH2D histogram
        TH2D *fake_confusion_hist = new TH2D(("fake_confusion slice "+std::to_string(sl_idx)).c_str(),("NuWro " + conf_title).c_str(),
                                num_bins, 0, num_bins,
                                num_bins, 0, num_bins);

        // Create a TH2D histogram
        TH2D *ac_hist = new TH2D(("ac_hist slice "+std::to_string(sl_idx)).c_str(),(ac_title).c_str(),
                                num_bins_ac, 0, num_bins_ac,
                                num_bins_ac, 0, num_bins_ac);
                                
        TH2D *trimmed_corr_unf_hist = new TH2D(("trimmed_corr_unf_hist slice "+std::to_string(sl_idx)).c_str(),(corr_unf_title).c_str(),
                                num_bins, 0, num_bins,
                                num_bins, 0, num_bins);

        TH2D *trimmed_cov_unf_hist = new TH2D(("trimmed_cov_unf_hist slice "+std::to_string(sl_idx)).c_str(),(cov_unf_title).c_str(),
                                num_bins, 0, num_bins,
                                num_bins, 0, num_bins);

        TH2D *trimmed_corr_hist = new TH2D(("trimmed_corr_hist slice "+std::to_string(sl_idx)).c_str(), corr_title.c_str(),
                                num_bins, 0, num_bins,
                                num_bins, 0, num_bins);

        TH2D *trimmed_cov_hist = new TH2D(("trimmed_cov_hist slice "+std::to_string(sl_idx)).c_str(), cov_title.c_str(),
                                num_bins, 0, num_bins,
                                num_bins, 0, num_bins);

        // Fill histogram with values from TMatrixD
        for (int i = 0; i < cv_confusion_mat.GetNrows(); i++) {
            for (int j = 0; j < cv_confusion_mat.GetNcols(); j++) {
                cv_confusion_hist->SetBinContent(i+1, j+1, cv_confusion_mat(i, j));
            }
        }

        if(using_fake_data)
        {
            // Fill fake histogram with values from TMatrixD
            for (int i = 0; i < fake_confusion_mat.GetNrows(); i++) {
                for (int j = 0; j < fake_confusion_mat.GetNcols(); j++) {
                    fake_confusion_hist->SetBinContent(i+1, j+1, fake_confusion_mat(i, j));
                }
            }
        }

        for (int i = 0; i < ac_hist_slice.GetNrows(); i++) {
            for (int j = 0; j < ac_hist_slice.GetNcols(); j++) {
                ac_hist->SetBinContent(i+1, j+1, ac_hist_slice.operator()(i, j));
            }
        }

        for(int i = 0; i < corr_unfolded_hist_slice.GetNrows(); i++) {
            for(int j = 0; j < corr_unfolded_hist_slice.GetNcols(); j++) {
                trimmed_corr_unf_hist->SetBinContent(i+1, j+1, corr_unfolded_hist_slice.operator()(i, j));
            }
        }
        
        std::cout<<"DEBUG trimmed_cov_unf_hist:"<<std::endl;
        for(int i = 0; i < cov_unfolded_hist_slice.GetNrows(); i++) {
            for(int j = 0; j < cov_unfolded_hist_slice.GetNcols(); j++) {
                trimmed_cov_unf_hist->SetBinContent(i+1, j+1, cov_unfolded_hist_slice.operator()(i, j));
                // std::cout<<"["<<i<<", "<<j<<"] = "<<cov_unfolded_hist_slice.operator()(i, j)<<std::endl;
            }
        }

        for(int i = 0; i < corr_hist_slice.GetNrows(); i++) {
            for(int j = 0; j < corr_hist_slice.GetNcols(); j++) {
                trimmed_corr_hist->SetBinContent(i+1, j+1, corr_hist_slice.operator()(i, j));
            }
        }

        std::cout<<"DEBUG trimmed_cov_hist:"<<std::endl;
        for(int i = 0; i < cov_hist_slice.GetNrows(); i++) {
            for(int j = 0; j < cov_hist_slice.GetNcols(); j++) {
                trimmed_cov_hist->SetBinContent(i+1, j+1, cov_hist_slice.operator()(i, j));
                // std::cout<<"["<<i<<", "<<j<<"] = "<<cov_hist_slice.operator()(i, j)<<std::endl;
            }
        }


        TCanvas* c_conf_1 = new TCanvas(("c_conf_1 slice "+std::to_string(sl_idx)).c_str(), "c confusion matrix", 800, 600);
        cv_confusion_hist->SetStats(false);
        cv_confusion_hist->SetMinimum(0);
        cv_confusion_hist->SetMaximum(1.0);
        // gStyle->SetPalette(util::CreateWhiteToBlueColorPalette(20));
        gStyle->SetPalette();

        // Set x and y axis labels to bin numbers
        for (int i = 1; i <= cv_confusion_hist->GetNbinsX(); i++) {
            // const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 1 || i == 6));
            // const auto binLabel = overFlow ? Form("%d*", i) : Form("%d", i);
            cv_confusion_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
        }
        for (int i = 1; i <= cv_confusion_hist->GetNbinsY(); i++) {
            // const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 1 || i == 6));
            // const auto binLabel = overFlow ? Form("%d*", i) : Form("%d", i);
            cv_confusion_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
        }
        // Increase the font size of the axis labels
        cv_confusion_hist->GetXaxis()->SetLabelSize(0.05);
        cv_confusion_hist->GetYaxis()->SetLabelSize(0.05);

        cv_confusion_hist->Draw("colz");
        // cv_hist_2d_slice.Draw("same TEXT");

        // Fill histogram with slice data
        for (int i = 0; i < num_bins; i++) {
            for (int j = 0; j < num_bins; j++) {
                double bin_content = cv_hist_2d_slice(i, j);
                TLatex* latex = new TLatex(cv_confusion_hist->GetXaxis()->GetBinCenter(i+1), cv_confusion_hist->GetYaxis()->GetBinCenter(j+1), Form("#splitline{%.1f%%}{(%.1f)}", 100*cv_confusion_hist->GetBinContent(i+1, j+1), bin_content));
                latex->SetTextFont(42);
                latex->SetTextSize(0.02);
                latex->SetTextAlign(22);
                latex->Draw();
            }
        }

        c_conf_1->SaveAs(("plots/cv_confusion_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro" + postfix + ".pdf").c_str());


        if(using_fake_data)
        {
            TCanvas* c_conf_fake = new TCanvas(("c_conf_fake slice "+std::to_string(sl_idx)).c_str(), "c confusion matrix", 800, 600);
            fake_confusion_hist->SetStats(false);
            fake_confusion_hist->SetMinimum(0);
            fake_confusion_hist->SetMaximum(1.0);
            fake_confusion_hist->Draw("colz");

            // Set x and y axis labels to bin numbers
            for (int i = 1; i <= fake_confusion_hist->GetNbinsX(); i++) {
                fake_confusion_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
            }
            for (int i = 1; i <= fake_confusion_hist->GetNbinsY(); i++) {
                fake_confusion_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
            }
            // Increase the font size of the axis labels
            fake_confusion_hist->GetXaxis()->SetLabelSize(0.05);
            fake_confusion_hist->GetYaxis()->SetLabelSize(0.05);

            // Fill histogram with slice data
            for (int i = 0; i < num_bins; i++) {
                for (int j = 0; j < num_bins; j++) {
                    double bin_content = fake_hist_2d_slice(i, j);
                    TLatex* latex = new TLatex(fake_confusion_hist->GetXaxis()->GetBinCenter(i+1), fake_confusion_hist->GetYaxis()->GetBinCenter(j+1), Form("#splitline{%.1f%%}{(%.1f)}", 100*fake_confusion_hist->GetBinContent(i+1, j+1), bin_content));
                    latex->SetTextFont(42);
                    latex->SetTextSize(0.02);
                    latex->SetTextAlign(22);
                    latex->Draw();
                }
            }

            c_conf_fake->SaveAs(("plots/fake_confusion_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro" + postfix + ".pdf").c_str());
        }


        TCanvas* c_ac_slice = new TCanvas(("c_ac_slice slice "+std::to_string(sl_idx)).c_str(), "c ac matrix", 800, 600);
        ac_hist->SetTitleSize(0.05, "t");
        ac_hist->SetStats(false);
        ac_hist->GetZaxis()->SetRangeUser(-0.5, 1.5); // Set the z range
        ac_hist->Draw("colz");

        // Print out number of bins
        // std::cout<<"DEBUG ac_hist.GetNbinsX() = "<<ac_hist->GetNbinsX()<<std::endl;

        // Set x and y axis labels to bin numbers
        for (int i = 1; i <= ac_hist->GetNbinsX(); i++) {
            // const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 1 || i == 6));
            const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 5));
            const auto binLabel = overFlow ? Form("%d*", i) : Form("%d", i);
            ac_hist->GetXaxis()->SetBinLabel(i, binLabel);
        }
        for (int i = 1; i <= ac_hist->GetNbinsY(); i++) {
            // const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 1 || i == 6));
            const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 5));
            const auto binLabel = overFlow ? Form("%d*", i) : Form("%d", i);
            ac_hist->GetYaxis()->SetBinLabel(i, binLabel);
        }
        // Increase the font size of the axis labels
        ac_hist->GetXaxis()->SetLabelSize(0.05);
        ac_hist->GetYaxis()->SetLabelSize(0.05);

        TColor::CreateGradientColorTable(nColors, stops, red, green, blue, 20);
        c_ac_slice->SaveAs(("plots/additional_smearing_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro" + postfix + ".pdf").c_str());

        // Plot for trimmed_corr_unf_hist
        TCanvas* c_corr_unf_slice = new TCanvas("c_corr_unf_slice", "c correlation matrix", 800, 600);
        util::CreateRedToBlueColorPalette(20);
        trimmed_corr_unf_hist->SetTitleSize(0.05, "t");
        trimmed_corr_unf_hist->GetZaxis()->SetRangeUser(-1, 1);
        trimmed_corr_unf_hist->SetStats(false);
        trimmed_corr_unf_hist->GetZaxis()->SetRangeUser(-1.0, 1.0); // Set the z range
        trimmed_corr_unf_hist->Draw("colz");

        // Set x and y axis labels to bin numbers
        for (int i = 1; i <= trimmed_corr_unf_hist->GetNbinsX(); i++) {
            trimmed_corr_unf_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
        }
        for (int i = 1; i <= trimmed_corr_unf_hist->GetNbinsY(); i++) {
            trimmed_corr_unf_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
        }
        // Increase the font size of the axis labels
        trimmed_corr_unf_hist->GetXaxis()->SetLabelSize(0.05);
        trimmed_corr_unf_hist->GetYaxis()->SetLabelSize(0.05);

        c_corr_unf_slice->SaveAs(("plots/unfolded_correlation_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro" + postfix + ".pdf").c_str());

        // Plot for trimmed_cov_unf_hist
        TCanvas* c_cov_unf_slice = new TCanvas("c_cov_unf_slice", "c covariance matrix", 800, 600);
        util::CreateRedToBlueColorPalette(20);
        trimmed_cov_unf_hist->SetTitleSize(0.05, "t");
        trimmed_cov_unf_hist->SetStats(false);
        trimmed_cov_unf_hist->Draw("colz");

        // Set x and y axis labels to bin numbers
        for (int i = 1; i <= trimmed_cov_unf_hist->GetNbinsX(); i++) {
            trimmed_cov_unf_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
        }
        for (int i = 1; i <= trimmed_cov_unf_hist->GetNbinsY(); i++) {
            trimmed_cov_unf_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
        }

        // Increase the font size of the axis labels
        trimmed_cov_unf_hist->GetXaxis()->SetLabelSize(0.05);
        trimmed_cov_unf_hist->GetYaxis()->SetLabelSize(0.05);

        c_cov_unf_slice->SaveAs(("plots/unfolded_covariance_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + postfix + ".pdf").c_str());

        // Plot for trimmed_corr_hist
        TCanvas* c_corr_slice = new TCanvas(("c_corr_slice slice "+std::to_string(sl_idx)).c_str(), "c correlation matrix", 800, 600);
        util::CreateRedToBlueColorPalette(20);
        trimmed_corr_hist->SetTitleSize(0.05, "t");
        trimmed_corr_hist->GetZaxis()->SetRangeUser(-1, 1);
        trimmed_corr_hist->SetStats(false);
        trimmed_corr_hist->Draw("colz");

        // Set x and y axis labels to bin numbers
        for (int i = 1; i <= trimmed_corr_hist->GetNbinsX(); i++) {
            trimmed_corr_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
        }
        for (int i = 1; i <= trimmed_corr_hist->GetNbinsY(); i++) {
            trimmed_corr_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
        }
        // Increase the font size of the axis labels
        trimmed_corr_hist->GetXaxis()->SetLabelSize(0.05);
        trimmed_corr_hist->GetYaxis()->SetLabelSize(0.05);

        c_corr_slice->SaveAs(("plots/correlation_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro" + postfix + ".pdf").c_str());

        // Plot for trimmed_cov_hist
        TCanvas* c_cov_slice = new TCanvas(("c_cov_slice slice "+std::to_string(sl_idx)).c_str(), "c covariance matrix", 800, 600);
        util::CreateRedToBlueColorPalette(20);
        trimmed_cov_hist->SetTitleSize(0.05, "t");
        // trimmed_cov_hist->GetZaxis()->SetRangeUser(-1, 1);
        trimmed_cov_hist->SetStats(false);
        trimmed_cov_hist->Draw("colz");

        // Set x and y axis labels to bin numbers
        for (int i = 1; i <= trimmed_cov_hist->GetNbinsX(); i++) {
            trimmed_cov_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
        }
        for (int i = 1; i <= trimmed_cov_hist->GetNbinsY(); i++) {
            trimmed_cov_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
        }
        // Increase the font size of the axis labels
        trimmed_cov_hist->GetXaxis()->SetLabelSize(0.05);
        trimmed_cov_hist->GetYaxis()->SetLabelSize(0.05);

        c_cov_slice->SaveAs(("plots/covariance_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + postfix + ".pdf").c_str());

        // Also use the GENIE CV model to do the same
        SliceHistogram *slice_cv = SliceHistogram::make_slice_histogram(
            genie_cv_truth_vec, slice, nullptr);

        // If present, also use the truth information from the fake data to do the same
        SliceHistogram *slice_truth = using_fake_data ? SliceHistogram::make_slice_histogram(fake_data_truth, slice, &fake_data_truth_cov) : nullptr;

        // Keys are legend labels, values are SliceHistogram objects containing true-space predictions from the corresponding generator models
        auto *slice_gen_map_ptr = new std::map<std::string, SliceHistogram *>();
        auto &slice_gen_map = *slice_gen_map_ptr;

        slice_gen_map["Unfolded Selection"] = slice_unf;
        slice_gen_map["MicroBooNE Tune"] = slice_cv;

        if (using_fake_data)
            slice_gen_map["NuWro Truth"] = slice_truth;

        for (const auto &pair : generator_truth_map)
            slice_gen_map[pair.first] = SliceHistogram::make_slice_histogram(*pair.second, slice, nullptr);

        auto *slice_gen_map_alt_ptr = new std::map<std::string, SliceHistogram *>(); // Used for stats only histograms
        auto &slice_gen_map_alt = *slice_gen_map_alt_ptr;
        slice_gen_map_alt["Unfolded Selection Stats Only"] = slice_unf_stats_only;

        int var_count = 0;
        double other_var_width = 1.;
        std::string diff_xsec_denom, diff_xsec_units_denom, diff_xsec_denom_latex, diff_xsec_units_denom_latex;

        for (const auto &ov_spec : slice.other_vars_)
        {
            double high = ov_spec.high_bin_edge_;
            double low = ov_spec.low_bin_edge_;
            const auto &var_spec = sb.slice_vars_.at(ov_spec.var_index_);

            if (high != low && std::abs(high - low) < BIG_DOUBLE)
            {
                ++var_count;
                other_var_width *= (high - low);
                diff_xsec_denom += 'd' + var_spec.name_;
                diff_xsec_denom_latex += " d" + var_spec.latex_name_;

                if (!var_spec.units_.empty())
                {
                    diff_xsec_units_denom += " / " + var_spec.units_;
                    diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;
                }
            }
        }

        for (size_t av_idx : slice.active_var_indices_)
        {
            const auto &var_spec = sb.slice_vars_.at(av_idx);
            if (var_spec.name_ != "true bin number")
            {
                var_count += slice.active_var_indices_.size();
                diff_xsec_denom += 'd' + var_spec.name_;
                diff_xsec_denom_latex += " d" + var_spec.latex_name_;

                if (!var_spec.units_.empty())
                {
                    diff_xsec_units_denom += " / " + var_spec.units_;
                    diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;
                }
            }
        }

        if(sl_idx == 7) var_count = 0; // TODO: Find a better way to do this


        if(other_var_width != 1.) throw std::runtime_error("Error: other_var_width is not 1.0 and has value " + std::to_string(other_var_width)); // Should not happen for this cross-section

        int num_slice_bins = slice_unf->hist_->GetNbinsX();
        TMatrixD trans_mat(num_slice_bins, num_slice_bins);
        for (int b = 0; b < num_slice_bins; ++b)
        {
            double width = slice_unf->hist_->GetBinWidth(b + 1) * other_var_width;
            trans_mat(b, b) = 1e39 / (width * integ_flux * num_Ar);
        }

        std::string slice_y_title = var_count > 0 ? "d" : "#sigma";
        std::string slice_y_latex_title = var_count > 0 ? "{$d" : "\\sigma";

        if (var_count > 1)
        {
            std::string var_count_str = "^{" + std::to_string(var_count) + "}";
            slice_y_title += var_count_str;
            slice_y_latex_title += var_count_str;
        }

        if (var_count > 0)
        {
            slice_y_title += "#sigma/" + diff_xsec_denom;
            slice_y_latex_title += "\\sigma / " + diff_xsec_denom_latex;
        }

        slice_y_title += " (10^{-39} cm^{2}" + diff_xsec_units_denom + " / Ar)";
        slice_y_latex_title += "\\text{ }(10^{-39}\\text{ cm}^{2}" + diff_xsec_units_denom_latex + " / \\mathrm{Ar})$}";

        for (auto &pair : slice_gen_map)
        {
            auto *slice_h = pair.second;
            slice_h->transform(trans_mat);
            slice_h->hist_->GetYaxis()->SetTitle(slice_y_title.c_str());
        }

        for (auto &pair : slice_gen_map_alt)
        {
            auto *slice_h = pair.second;
            slice_h->transform(trans_mat);
            slice_h->hist_->GetYaxis()->SetTitle(slice_y_title.c_str());
        }

        // for (auto &sh_cov_pair : sh_cov_map)
        // {
        //     auto &slice_h = sh_cov_pair.second;
        //     slice_h->transform(trans_mat);
        // }

        std::map<std::string, SliceHistogram::Chi2Result> chi2_map;
        for (const auto &pair : slice_gen_map)
        {
            const auto &name = pair.first;
            const auto *slice_h = pair.second;
            chi2_map[name] = (name == "Unfolded Selection") ? SliceHistogram::Chi2Result() : slice_h->get_chi2(*slice_gen_map.at("Unfolded Selection"));
        }

        TCanvas *c1 = new TCanvas("c1", "c1");//, 800, 600); // 800x600 pixels
        gStyle->SetLegendBorderSize(0);
        const auto rightMargin = 0.25;

        const auto noRatioPlot = (sl_idx == 7);

        TPad* pad2 = new TPad(("pad2 slice "+std::to_string(sl_idx)).c_str(), "", 0, 0.0, 1.0, 0.3);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.35);
        pad2->SetRightMargin(rightMargin);
        if(!noRatioPlot) pad2->Draw();
        pad2->cd();

        // Create the ratio histogram
        TH1D* h_ratio = dynamic_cast<TH1D*>(slice_unf->hist_.get()->Clone("h_ratio"));
        h_ratio->SetStats(false);

        // Create a temporary clone of slice_cv->hist_ with zero errors
        TH1D* h_cv_no_errors = dynamic_cast<TH1D*>(slice_cv->hist_.get()->Clone("h_cv_no_errors"));
        for (int i = 1; i <= h_cv_no_errors->GetNbinsX(); ++i) {
            h_cv_no_errors->SetBinError(i, 0.0);
        }


        // Draw a horizontal dashed line at ratio == 1
        TLine* line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
        line->SetLineColor(kAzure - 7);
        line->SetLineWidth(2);
        line->SetLineStyle(5);


        h_ratio->Divide(h_cv_no_errors);
        h_ratio->SetLineWidth(2);
        h_ratio->SetLineColor(kBlack);
        h_ratio->SetMarkerStyle(kFullCircle);
        h_ratio->SetMarkerSize(0.7);
        h_ratio->SetTitle("");


        // x-axis
        h_ratio->GetXaxis()->SetTitle(slice_unf->hist_->GetXaxis()->GetTitle());
        h_ratio->GetXaxis()->SetTitleOffset(1.3); // Hide x-axis        
        h_ratio->GetXaxis()->CenterTitle(true);
        h_ratio->GetXaxis()->SetLabelSize(0.1);
        h_ratio->GetXaxis()->SetTitleSize(0.08);
        h_ratio->GetXaxis()->SetTickLength(0.05);
        // h_ratio->GetXaxis()->SetTitleOffset(0.9);

        if(sl_idx == 7)
        {
            h_ratio->GetXaxis()->SetLabelOffset(999); // Hide x-axis labels
            h_ratio->GetXaxis()->SetTickLength(0); // Hide x-axis ticks
        }

        // y-axis
        h_ratio->GetYaxis()->SetTitle("Ratio");
        h_ratio->GetYaxis()->CenterTitle(true);
        h_ratio->GetYaxis()->SetLabelSize(0.08);
        h_ratio->GetYaxis()->SetTitleSize(0.08);
        h_ratio->GetYaxis()->SetTitleOffset(0.4);

        TH1D* slice_unf_stats_only_ratio = dynamic_cast<TH1D*>(slice_unf_stats_only->hist_.get()->Clone("slice_unf_stats_only_ratio"));
        // Print out the bin sizes for slice_unf_stats_only_ratio and h_cv_no_errors separately to check wethether they are the same
        std::cout << "Bin information for slice_unf_stats_only_ratio:" << std::endl;
        for (int i = 1; i <= slice_unf_stats_only_ratio->GetNbinsX(); ++i) {
            std::cout << "Bin " << i << ": Lower edge = " << slice_unf_stats_only_ratio->GetBinLowEdge(i)
                    << ", Upper edge = " << slice_unf_stats_only_ratio->GetBinLowEdge(i) + slice_unf_stats_only_ratio->GetBinWidth(i)
                    << ", Width = " << slice_unf_stats_only_ratio->GetBinWidth(i) << std::endl;
        }

        std::cout << "\nBin information for h_cv_no_errors:" << std::endl;
        for (int i = 1; i <= h_cv_no_errors->GetNbinsX(); ++i) {
            std::cout << "Bin " << i << ": Lower edge = " << h_cv_no_errors->GetBinLowEdge(i)
                    << ", Upper edge = " << h_cv_no_errors->GetBinLowEdge(i) + h_cv_no_errors->GetBinWidth(i)
                    << ", Width = " << h_cv_no_errors->GetBinWidth(i) << std::endl;
        }
        slice_unf_stats_only_ratio->Divide(h_cv_no_errors);
        slice_unf_stats_only_ratio->SetLineWidth(4);
        slice_unf_stats_only_ratio->SetLineColor(kBlack);
        slice_unf_stats_only_ratio->SetMarkerStyle(kFullCircle);
        slice_unf_stats_only_ratio->SetMarkerSize(0.7);

        // Set y-axis range
        double min = std::numeric_limits<double>::max();
        double max = - std::numeric_limits<double>::max();

        for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
            double binContent = h_ratio->GetBinContent(i);
            double binError = h_ratio->GetBinError(i);
            if (binContent - binError < min) {
                min = binContent - binError;
            }
            if (binContent + binError > max) {
                max = binContent + binError;
            }
        }

        h_ratio->GetYaxis()->SetRangeUser(0.95*min, 1.05*max);
        h_ratio->Draw("e");

        for(const auto& generator : generators)
        {
            const auto gen_hist = get_generator_hist(generator.path, sl_idx, 1.0/generator.scaling);
            if (gen_hist) {

                if(USE_ADD_SMEAR) 
                {
                    multiply_1d_hist_by_matrix(&ac_hist_slice, gen_hist);
                }

                // Normalise by bin width
                for (int i = 1; i <= gen_hist->GetNbinsX(); ++i) {
                    double bin_content = gen_hist->GetBinContent(i);
                    double bin_error = gen_hist->GetBinError(i);
                    // double bin_width = gen_hist->GetXaxis()->GetBinWidth(i);
                    double bin_width = slice_unf->hist_->GetXaxis()->GetBinWidth(i); // Get the bin width from the slice_unf histogram
                    gen_hist->SetBinContent(i, bin_content / bin_width);
                    gen_hist->SetBinError(i, bin_error / bin_width);
                }

                // Clone gen_hist and create a 'ratio' version and divide that by h_cv_no_errors
                TH1D* gen_hist_ratio = dynamic_cast<TH1D*>(gen_hist->Clone("gen_hist_ratio"));
                gen_hist_ratio->Divide(h_cv_no_errors); // Warnings about different bin limits come from floating point differences in bin sizes

                // Update min and max based on gen_hist_ratio
                for (int i = 1; i <= gen_hist_ratio->GetNbinsX(); ++i) {
                    double binContent = gen_hist_ratio->GetBinContent(i);
                    double binError = gen_hist_ratio->GetBinError(i);
                    if (binContent - binError < min) {
                        min = binContent - binError;
                    }
                    if (binContent + binError > max) {
                        max = binContent + binError;
                    }
                }
                drawHistogramWithBand(gen_hist_ratio, generator, 0.3, true);
            }
        }

        h_ratio->GetYaxis()->SetRangeUser(0.95 * std::min(min, 1.0), 1.05 * std::max(max, 1.0));

        if (using_fake_data)
        {
            TH1D* h_ratio_truth = dynamic_cast<TH1D*>(slice_truth->hist_.get()->Clone("h_ratio_truth"));

            h_ratio_truth->SetStats(false);
            h_ratio_truth->Divide(h_cv_no_errors);
            h_ratio_truth->SetLineColor(kGreen);
            h_ratio_truth->SetLineWidth(2);
            h_ratio_truth->SetLineStyle(5);
            min = TMath::Min(min, 0.9 * h_ratio_truth->GetMinimum());
            max = TMath::Max(max, 1.1 * h_ratio_truth->GetMaximum());
            h_ratio_truth->GetYaxis()->SetRangeUser(min, max);
            h_ratio_truth->Draw("hist same");
        }

        line->Draw();
        h_ratio->Draw("e same"); // draw again to be over everything

        slice_unf_stats_only_ratio->Draw("EX0 same");

        // Go back to the main canvas before creating the second pad
        c1->cd();

        TPad* pad1 = new TPad(("pad1 slice "+std::to_string(sl_idx)).c_str(), "", 0.0, noRatioPlot ? 0.1 : 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.03);
        pad1->SetTopMargin(0.1);
        pad1->SetRightMargin(rightMargin);
        pad1->Draw();
        pad1->cd();


        slice_unf->hist_->SetLineColor(kBlack);
        // slice_unf->hist_->SetLineWidth(5);
        slice_unf->hist_->SetMarkerStyle(kFullCircle);
        slice_unf->hist_->SetMarkerSize(0.7);
        slice_unf->hist_->SetStats(false);
        slice_unf->hist_->SetLineWidth(2);


        double ymax = -DBL_MAX;
        slice_unf->hist_->Draw("e");

        if(sl_idx == 7)
        {
            slice_unf->hist_->GetXaxis()->SetLabelOffset(999); // Hide x-axis labels
            slice_unf->hist_->GetXaxis()->SetTickLength(0); // Hide x-axis ticks
        }

        for (const auto &pair : slice_gen_map)
        {
            const auto &name = pair.first;
            const auto *slice_h = pair.second;

            double max = slice_h->hist_->GetMaximum();
            // if (max > ymax)
            //     ymax = max;

            for (int i = 1; i <= slice_h->hist_->GetNbinsX(); ++i) {
                double binContent = slice_h->hist_->GetBinContent(i);
                double binError = (name == "Unfolded Selection") ? slice_h->hist_->GetBinError(i) : 0;
                if (binContent + binError > ymax) {
                    ymax = binContent + binError;
                }
            }

            if (name == "Unfolded Selection" || name == "NuWro Truth" || name == "MicroBooNE Tune")
                continue;

            const auto &file_info = truth_file_map.at(name);
            slice_h->hist_->SetLineColor(file_info.color_);
            slice_h->hist_->SetLineStyle(file_info.style_);
            slice_h->hist_->SetLineWidth(3);
            slice_h->hist_->Draw("hist same");
        }

        drawHistogramWithBand(slice_cv->hist_.get(), kAzure - 7, 2, 5, 0.3, true);

        if (using_fake_data)
        {
            const auto& slice_truth_cov_mat_ptr = slice_truth->cmat_.get_matrix(); // Get the pointer first
            if (!slice_truth_cov_mat_ptr) 
            {
                throw std::runtime_error("Error: slice_truth_cov_mat_ptr is nullptr");
            }
            const auto& slice_truth_cov_mat = *slice_truth_cov_mat_ptr; // Now dereference

            //Set the slice_truth histogram errors to the sqrt of the diagonal slice_truth->cmat_ values
            for (int i = 1; i <= slice_truth->hist_->GetNbinsX(); ++i) {
                double error = sqrt(slice_truth_cov_mat(i-1, i-1));
                slice_truth->hist_->SetBinError(i, error);
            }

            drawHistogramWithBand(slice_truth->hist_.get(), kGreen, 2, 5, 0.3, true);
            // drawHistogramWithBand(h_slice_truth_with_error, kGreen, 2, 5, 0.2, false);
        }

        std::map<std::string, SliceHistogram::Chi2Result> gen_metrics;
        for(const auto& generator : generators)
        {
            const auto gen_hist = get_generator_hist(generator.path, sl_idx, 1.0/generator.scaling);
            if (gen_hist) {
                
                if(USE_ADD_SMEAR) 
                {
                    multiply_1d_hist_by_matrix(&ac_hist_slice, gen_hist);
                }

                // Normalise by bin width
                for (int i = 1; i <= gen_hist->GetNbinsX(); ++i) {
                    double bin_content = gen_hist->GetBinContent(i);
                    double bin_error = gen_hist->GetBinError(i);
                    // double bin_width = gen_hist->GetXaxis()->GetBinWidth(i);
                    double bin_width = slice_unf->hist_->GetXaxis()->GetBinWidth(i); // Get the bin width from the slice_unf histogram
                    gen_hist->SetBinContent(i, bin_content / bin_width);
                    gen_hist->SetBinError(i, bin_error / bin_width);
                }

                // Create the generator slice with nullptr convariance matrix
                SliceHistogram *gen_slice_h = SliceHistogram::slice_histogram_from_histogram(*gen_hist);

                if (gen_slice_h->hist_->GetNbinsX() != slice_gen_map.at("Unfolded Selection")->hist_->GetNbinsX()) {
                    std::ostringstream oss;
                    oss << "Error: gen_slice_h->hist_->GetNbinsX() (" << gen_slice_h->hist_->GetNbinsX() << ") != slice_gen_map.at(\"Unfolded NuWro Reco\")->hist_->GetNbinsX() (" << slice_gen_map.at("Unfolded Selection")->hist_->GetNbinsX() << ")";
                    throw std::runtime_error(oss.str());
                }

                for (int i = 1; i <= gen_slice_h->hist_->GetNbinsX(); ++i) {
                    double binContent = gen_slice_h->hist_->GetBinContent(i);
                    double binError = 0;
                    if (binContent + binError > ymax) {
                        ymax = binContent + binError;
                    }
                }

                // Compute chi2 to unfolded NuWro Reco
                const auto metrics = gen_slice_h->get_chi2(*slice_gen_map.at("Unfolded Selection"));
                gen_metrics[generator.name] = metrics;
                delete gen_slice_h;

                drawHistogramWithBand(gen_hist, generator, 0.3, true);
            }
        }

        slice_unf->hist_->GetYaxis()->SetRangeUser(0., ymax * 1.05);
        slice_unf->hist_->Draw("E same");
        // slice_unf->hist_->SetTitle("Unfolded NuWro CC1#pi^{#pm}Xp");
        slice_unf->hist_->SetTitle(""); // No title
        slice_unf->hist_->GetXaxis()->SetLabelOffset(999); // Hide x-axis
        slice_unf->hist_->GetXaxis()->SetTitleOffset(999); // Hide x-axis title

        gStyle->SetEndErrorSize(3);
        slice_unf_stats_only->hist_->SetLineColor(kBlack);
        slice_unf_stats_only->hist_->SetMarkerStyle(kFullCircle);
        slice_unf_stats_only->hist_->SetMarkerSize(0.7);
        slice_unf_stats_only->hist_->SetStats(false);
        slice_unf_stats_only->hist_->SetLineWidth(4);
        slice_unf_stats_only->hist_->Draw("EX0 same");

        c1->cd();

        // TLegend *lg = new TLegend(0.6, 0.2);
        TLegend *lg = new TLegend(1 - rightMargin + 0.02, 0.09, 1 - 0.02, 0.93);

        // lg->AddEntry((TObject*)0, "Generator Predictions", "");
        for(const auto& generator : generators)
        {
            const auto gen_hist = get_generator_hist(generator.path, sl_idx, 1.0/generator.scaling);
            if (gen_hist) {
                gen_hist->SetLineColor(generator.lineColor); // Set the line color
                gen_hist->SetLineWidth(generator.lineWidth); // Set the line width
                gen_hist->SetLineStyle(generator.lineStyle); // Set the line style

                const auto metrics = gen_metrics[generator.name];
                // Format the chi2 values and append to the generator name
                std::ostringstream oss;
                oss << "#splitline{" << generator.name << "}{" 
                    << "#splitline{#chi^{2} = " << (metrics.chi2_>= 0.01 && metrics.chi2_ < 100 ? std::fixed : std::scientific) << std::setprecision(2) << metrics.chi2_ << " / " << metrics.num_bins_ << " bin" << (metrics.num_bins_ > 1 ? "s" : "") << "}{"; 
                
                if (metrics.num_bins_ > 1)
                    oss << "p = " << metrics.p_value_;
                
                oss << "}}";
                const std::string label = oss.str();

                lg->AddEntry(gen_hist, label.c_str(), "l");
            }
        }

        for (const auto &pair : slice_gen_map)
        {
            const auto &name = pair.first;
            const auto *slice_h = pair.second;

            // std::string label = name;
            const auto &chi2_result = chi2_map.at(name);
            // std::ostringstream oss;
            // oss << std::setprecision(3) << chi2_result.chi2_ << " / " << chi2_result.num_bins_ << " bin";
            // if (chi2_result.num_bins_ > 1)
            //     oss << "s; p-value = " << chi2_result.p_value_;

            // label += (name != "Unfolded Selection") ? ": #chi^{2} = " + oss.str() : "";

            std::ostringstream oss;
            if((name != "Unfolded Selection"))
            {
                oss << "#splitline{";
                if (name == "NuWro Truth") oss << "#splitline{Fake Data Truth}{(NuWro 19.02.1)}";
                else if (name == "MicroBooNE Tune") oss << "#splitline{GENIE 3.0.6}{G18_10a_02_11a}";
                else oss << name;

                oss << "}{" << "#splitline{#chi^{2} = " << (chi2_result.chi2_>= 0.01 && chi2_result.chi2_ < 100 ? std::fixed : std::scientific) << std::setprecision(2) << chi2_result.chi2_ << " / " << chi2_result.num_bins_ << " bin" << (chi2_result.num_bins_ > 1 ? "s" : "") << "}{"; 
                
                if (chi2_result.num_bins_ > 1)
                    oss << "p = " << chi2_result.p_value_;
                
                oss << "}}";
            }
            else
            {
                oss << name;
            }
            const std::string label = oss.str();

            lg->AddEntry(slice_h->hist_.get(), label.c_str(), (name == "Unfolded Selection") ? "lp" : "l");
        }

        std::ostringstream headerss;
        headerss << "#splitline{NuWro Fake Data}"
                         << "{" << toLatexScientific(total_pot) << " POT}";
                 
        lg->SetHeader(headerss.str().c_str());
        // Increase the font size for the legend header
        // (see https://root-forum.cern.ch/t/tlegend-headers-font-size/14434)
        TLegendEntry* lg_header = dynamic_cast< TLegendEntry* >(
            lg->GetListOfPrimitives()->First() );
        lg_header->SetTextSize( 0.03 );
        lg->Draw("same");
        c1->Update();

        std::string out_pdf_name = "plots/plot_unfolded_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro" + postfix + ".pdf";
        c1->SaveAs(out_pdf_name.c_str());


        // A separate plot for the event rate
        SliceHistogram *slice_unf_event_rate_stats_only = SliceHistogram::make_slice_histogram(
                    *result.unfolded_signal_, slice, unfolded_cov_matrix_map.at("DataStats").get());

        // Selected reco before unfolding
        SliceHistogram *slice_reco_selected_event_rate = SliceHistogram::make_slice_histogram(
                    *syst.data_hists_.at( NFT::kOnBNB ), slice, nullptr);


        SliceHistogram *slice_reco_bkgd_subtracted_selected_event_rate = SliceHistogram::make_slice_histogram(
                       *data_signal, slice, cov_matrix_map->at("DataStats").get_matrix().get());

        UnfolderHelper::event_rate_plot(slice_reco_selected_event_rate->hist_.get(), sl_idx, "_selected_reco"+postfix);
        UnfolderHelper::event_rate_plot(slice_reco_bkgd_subtracted_selected_event_rate->hist_.get(), sl_idx, "_selected_reco_bkgd_subtracted"+postfix);
        UnfolderHelper::event_rate_plot(slice_unf_event_rate_stats_only->hist_.get(), sl_idx, "_unfolded"+postfix);


        // Dump the unfolded results to text files compatible with PGFPlots
        std::map<std::string, std::vector<double>> slice_hist_table;
        std::map<std::string, std::string> slice_params_table;

        dump_slice_variables(sb, sl_idx, slice_params_table);

        for (const auto &pair : slice_gen_map)
        {
            const auto hist_name = samples_to_hist_names.at(pair.first);
            const auto *slice_hist = pair.second;
            bool include_coords_and_error = (hist_name == "UnfData");

            dump_slice_histogram(hist_name, *slice_hist, slice, slice_hist_table,
                                 include_coords_and_error, include_coords_and_error);
        }


        dump_slice_plot_limits(*slice_unf, *slice_cv, slice, slice_params_table);

        // dump_slice_errors("UnfData", slice, sh_cov_map, slice_hist_table);

        // Dump the chi^2 test results
        for (const auto &chi2_pair : chi2_map)
        {
            const auto hist_name = samples_to_hist_names.at(chi2_pair.first);
            const auto &chi2_result = chi2_pair.second;

            // Comparing the data histogram to itself is trivial, so skip it
            if (hist_name != "UnfData")
            {
                slice_params_table[hist_name + "_chi2"] = std::to_string(chi2_result.chi2_);
                slice_params_table[hist_name + "_pvalue"] = std::to_string(chi2_result.p_value_);
            }
        }

        // Dump the total data POT and number of bins in the slice
        slice_params_table["bnb_data_pot"] = std::to_string(total_pot);
        slice_params_table["num_bins"] = std::to_string(num_slice_bins);

        // Dump a LaTeX title for the y-axis
        slice_params_table["y_axis_title"] = slice_y_latex_title;

        // Prepare the output file prefix
        std::string output_file_prefix = "dump/pgfplots_slice_" + std::string(3 - std::to_string(sl_idx).length(), '0') + std::to_string(sl_idx);

        write_pgfplots_files(output_file_prefix, slice_hist_table, slice_params_table);

        // Make unfolded fractional uncertainty plots
        UnfolderHelper::fractional_uncertainty_plot(slice, unfolded_cov_matrix_map, result.unfolded_signal_.get(), sl_idx, cov_mat_keys_total, "total", postfix + "_unfolded_total");
        UnfolderHelper::fractional_uncertainty_plot(slice, unfolded_cov_matrix_map, result.unfolded_signal_.get(), sl_idx, cov_mat_keys_xsec, "xsec_total", postfix + "_unfolded_xsec");
    } // slices
    delete corr_unfolded_hist, corr_hist, cov_unfolded_hist;
    std::cout << "---- All done ----" << std::endl;
    return 0;
}

int main()
{
    unfold_nuwro();
    return 0;
}
