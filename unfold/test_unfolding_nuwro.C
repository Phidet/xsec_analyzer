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
// #include "../DAgostiniUnfolder.hh"
#include "../FiducialVolume.hh"
#include "../MatrixUtils.hh"
#include "../MCC9SystematicsCalculator.hh"
#include "../NormShapeCovMatrix.hh"
#include "../PGFPlotsDumpUtils.hh"
#include "../SliceBinning.hh"
#include "../SliceHistogram.hh"
#include "../WienerSVDUnfolder.hh"
#include "../utils.hh"

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
    {"Unfolded NuWro Reco", "UnfData"},
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
                          const std::map<std::string, std::unique_ptr<TMatrixD>> &unf_cov_matrix_map,
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

void test_unfolding_nuwro()
{
    // Initialize the FilePropertiesManager and tell it to treat the NuWro MC ntuples as if they were data
    auto &fpm = FilePropertiesManager::Instance();
  
    // fpm.load_file_properties("../nuwro_file_properties_closure.txt");
    // const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/univmake_output_nuwro_v5_noext_nodirt_run1.root");
    // // Do the systematics calculations in preparation for unfolding
    // std::cout << "DEBUG test_unfolding_nuwro - Point 0" << std::endl;
    // const auto *syst_ptr = new MCC9SystematicsCalculator(respmat_file_name, "../systcalc_unfold_fd_closure.conf");

    fpm.load_file_properties("../nuwro_file_properties.txt");
    const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_9/univmake_output_nuwro_sidband_reduced_28Dec23.root");
    // const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_9/univmake_output_nuwro_with_sidband_22Dec23.root");
    // Do the systematics calculations in preparation for unfolding
    const auto *syst_ptr = new MCC9SystematicsCalculator(respmat_file_name, "../systcalc_unfold_fd_min.conf");

    const auto &syst = *syst_ptr;

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
    std::cout<<"DEBUG U-1 V3"<<std::endl;
    auto &matrix_map = *matrix_map_ptr;
    std::cout<<"DEBUG U-0.9"<<std::endl;
    auto *cov_mat = matrix_map.at("total").cov_matrix_.get();
    std::cout<<"DEBUG U0"<<std::endl;




    // Define your custom labels and intervals
    const Int_t n = 9;
    Double_t bins[n+1] = {0, 11, 26, 32, 39, 49, 55, 62, 65, 66};
    const Int_t overUnderFlowN = 3;
    Double_t overUnderFlowBins[overUnderFlowN] = {31, 49, 54};
    const Char_t *labels[n] = {"cos(#theta_{#mu})", "#phi_{#mu}", "p_{#mu}", "cos(#theta_{#pi})", "#phi_{#pi}", "p_{#pi}^{**}", "#theta_{#pi #mu}", "N_{p}", "Total"};

    // Set the color palette
    util::CreateRedToBlueColorPalette(20);
    // gStyle->SetPalette(util::CreateRedToBlueColorPalette(20));
    gStyle->SetTitleFontSize(0.05); // Set the title font size to 0.05

    TMatrixD* corr = new TMatrixD(util::CovarianceMatrixToCorrelationMatrix(util::TH2DToTMatrixD(*cov_mat)));
    TH2D* corr_hist = new TH2D(util::TMatrixDToTH2D(*corr, "corr", "Correlation Matrix", 0, corr->GetNrows(), 0, corr->GetNcols()));
    delete corr; // Delete the dynamically allocated memory

    TCanvas *cm2 = new TCanvas("cm2", "Canvas", 800, 600);
    corr_hist->SetTitleSize(0.05, "t"); // Set the title font size to 0.05
    corr_hist->SetTitle("Correlation matrix");
    corr_hist->GetZaxis()->SetRangeUser(-1, 1);
    corr_hist->GetXaxis()->SetRangeUser(0, 66); // Set the range of the x-axis
    corr_hist->GetXaxis()->SetTitle("Reco Bins"); // Set x-axis label
    corr_hist->GetYaxis()->SetTitle("Reco Bins"); // Set y-axis label
    corr_hist->SetStats(kFALSE);
    corr_hist->Draw("COLZ");

    // Draw vertical and horizontal lines at the bin edges
    for (Int_t i = 1; i < n; i++) {
        TLine *vline = new TLine(bins[i], 0, bins[i], corr_hist->GetNbinsY());
        vline->SetLineColor(kBlack);
        vline->Draw();

        TLine *hline = new TLine(0, bins[i], corr_hist->GetNbinsX(), bins[i]);
        hline->SetLineColor(kBlack);
        hline->Draw();
    }

    for (Int_t i = 1; i <= n; i++) {
        // Draw white dotted lines from bins[i-1] to bins[i]
        TLine *vline_dotted1 = new TLine(bins[i], bins[i-1], bins[i], bins[i]);
        vline_dotted1->SetLineColor(kWhite);
        vline_dotted1->SetLineStyle(2); // Set line style to dotted
        vline_dotted1->Draw();

        TLine *hline_dotted1 = new TLine(bins[i-1], bins[i], bins[i], bins[i]);
        hline_dotted1->SetLineColor(kWhite);
        hline_dotted1->SetLineStyle(2); // Set line style to dotted
        hline_dotted1->Draw();

        if(i<n)
        {
            TLine *vline_dotted2 = new TLine(bins[i], bins[i+1], bins[i], bins[i]);
            vline_dotted2->SetLineColor(kWhite);
            vline_dotted2->SetLineStyle(2); // Set line style to dotted
            vline_dotted2->Draw();

            TLine *hline_dotted2 = new TLine(bins[i+1], bins[i], bins[i], bins[i]);
            hline_dotted2->SetLineColor(kWhite);
            hline_dotted2->SetLineStyle(2); // Set line style to dotted
            hline_dotted2->Draw();
        }
    }

    // Add labels in the middle of the intervals
    for (Int_t i = 0; i < n; i++) {
        Double_t midPoint = i==n-1 ? bins[i+1]+1 : (bins[i] + bins[i+1]) / 2.0;
        TLatex *text = new TLatex(midPoint, 1.03*corr_hist->GetNbinsY(), labels[i]);
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->Draw();
        // delete text;
    }

    // Add asterisk for over-/underrflow bins
    for (Int_t i = 0; i < overUnderFlowN; i++) {
        Double_t midPoint = overUnderFlowBins[i] + 0.5;
        TLatex *text = new TLatex(midPoint, 1.0*corr_hist->GetNbinsY(), "*");
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->SetTextColor(kGray+3); // Set text color to grey
        text->Draw();
        // delete text;
    }

    // Add footnote
    TLatex *footnote2 = new TLatex(0, -0.1*corr_hist->GetNbinsY(), "* Under-/overflow bin; ** Selection subset");
    footnote2->SetTextSize(0.02); // Set text size to something smaller
    footnote2->SetTextColor(kGray+3);
    // footnote->SetTextAlign(22); // Center alignment
    footnote2->Draw();

    cm2->SaveAs("plots/plot_slice_entire_corr_nuwro.pdf");

    std::cout<<"DEBUG U0.1"<<std::endl;

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

    auto result = unfolder->unfold(syst);
    std::cout<<"DEBUG U2"<<std::endl;

    // Propagate all defined covariance matrices through the unfolding procedure
    TMatrixD err_prop_tr(TMatrixD::kTransposed, *result.err_prop_matrix_);

    // std::map<std::string, std::unique_ptr<TMatrixD>> unfolded_cov_matrix_map;

    // for (const auto &matrix_pair : matrix_map)
    // {
    //     std::cout << "DEBUG matrix_pair.first: " << matrix_pair.first << std::endl;
    //     unfolded_cov_matrix_map[matrix_pair.first] = std::make_unique<TMatrixD>(
    //         *result.err_prop_matrix_, TMatrixD::EMatrixCreatorsOp2::kMult,
    //         TMatrixD(*matrix_pair.second.get_matrix(), TMatrixD::EMatrixCreatorsOp2::kMult, err_prop_tr));
    // }

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

    corr_unfolded_hist->SetTitleSize(0.05, "t"); // Set the title font size to 0.05
    corr_unfolded_hist->SetTitle("Correlation matrix after blockwise unfolding");
    corr_unfolded_hist->GetZaxis()->SetRangeUser(-1, 1);
    corr_unfolded_hist->GetXaxis()->SetTitle("Truth Bins"); // Set x-axis label
    corr_unfolded_hist->GetYaxis()->SetTitle("Truth Bins"); // Set y-axis label
    corr_unfolded_hist->SetStats(kFALSE);
    corr_unfolded_hist->Draw("COLZ");

    // Draw vertical and horizontal lines at the bin edges
    for (Int_t i = 1; i < n; i++) {
        TLine *vline = new TLine(bins[i], 0, bins[i], corr_unfolded_hist->GetNbinsY());
        vline->SetLineColor(kBlack);
        vline->Draw();

        TLine *hline = new TLine(0, bins[i], corr_unfolded_hist->GetNbinsX(), bins[i]);
        hline->SetLineColor(kBlack);
        hline->Draw();
    }

    for (Int_t i = 1; i <= n; i++) {
        // Draw white dotted lines from bins[i-1] to bins[i]
        TLine *vline_dotted1 = new TLine(bins[i], bins[i-1], bins[i], bins[i]);
        vline_dotted1->SetLineColor(kWhite);
        vline_dotted1->SetLineStyle(2); // Set line style to dotted
        vline_dotted1->Draw();

        TLine *hline_dotted1 = new TLine(bins[i-1], bins[i], bins[i], bins[i]);
        hline_dotted1->SetLineColor(kWhite);
        hline_dotted1->SetLineStyle(2); // Set line style to dotted
        hline_dotted1->Draw();

        if(i<n)
        {
            TLine *vline_dotted2 = new TLine(bins[i], bins[i+1], bins[i], bins[i]);
            vline_dotted2->SetLineColor(kWhite);
            vline_dotted2->SetLineStyle(2); // Set line style to dotted
            vline_dotted2->Draw();

            TLine *hline_dotted2 = new TLine(bins[i+1], bins[i], bins[i], bins[i]);
            hline_dotted2->SetLineColor(kWhite);
            hline_dotted2->SetLineStyle(2); // Set line style to dotted
            hline_dotted2->Draw();
        }
    }

    // Add labels in the middle of the intervals
    for (Int_t i = 0; i < n; i++) {
        Double_t midPoint = i==n-1 ? bins[i+1] + 1 : (bins[i] + bins[i+1]) / 2.0;
        TLatex *text = new TLatex(midPoint, 1.03*corr_unfolded_hist->GetNbinsY(), labels[i]);
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->Draw();
        // delete text;
    }

    // Add asterisk for over-/underrflow bins
    for (Int_t i = 0; i < overUnderFlowN; i++) {
        Double_t midPoint = overUnderFlowBins[i] + 0.5;
        TLatex *text = new TLatex(midPoint, 1.0*corr_unfolded_hist->GetNbinsY(), "*");
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->SetTextColor(kGray+3); // Set text color to grey
        text->Draw();
        // delete text;
    }

    // Add footnote
    TLatex *footnote1 = new TLatex(0, -0.1*corr_unfolded_hist->GetNbinsY(), "* Under-/overflow bin; ** Selection subset");
    footnote1->SetTextSize(0.02); // Set text size to something smaller
    footnote1->SetTextColor(kGray+3);
    // footnote->SetTextAlign(22); // Center alignment
    footnote1->Draw();

    cm1->SaveAs("plots/plot_unfolded_slice_entire_corr_nuwro.pdf");

    // Add the blockwise decomposed matrices into the map
    // unfolded_cov_matrix_map["total_blockwise_norm"] = std::make_unique<TMatrixD>(bd_ns_covmat.norm_);
    // unfolded_cov_matrix_map["total_blockwise_shape"] = std::make_unique<TMatrixD>(bd_ns_covmat.shape_);
    // unfolded_cov_matrix_map["total_blockwise_mixed"] = std::make_unique<TMatrixD>(bd_ns_covmat.mixed_);

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
    if (using_fake_data)
    {
        for (int b = 0; b < num_true_signal_bins; ++b)
        {
            double true_evts = fake_data_truth_hist->GetBinContent(b + 1);
            fake_data_truth(b, 0) = true_evts;
        }
    }

    // Save the GENIE CV model (before A_C multiplication) using a column vector
    // of event counts
    TMatrixD genie_cv_truth_vec(num_true_signal_bins, 1);
    for (int b = 0; b < num_true_signal_bins; ++b)
    {
        double true_evts = genie_cv_truth->GetBinContent(b + 1);
        genie_cv_truth_vec(b, 0) = true_evts;
    }

    std::cout<<"DEBUG U2"<<std::endl;
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
    TH2D h_A_C = util::TMatrixDToTH2D(*A_C, "h_A_C", "Additional smearing matrix", 0, A_C->GetNcols(), 0, A_C->GetNrows());

    TCanvas *c_ac = new TCanvas("c_ac","A_C Matrix",200,10,700,500);
    c_ac->SetRightMargin(0.15);
    h_A_C.SetStats(0); // Disable the statistics box
    h_A_C.GetZaxis()->SetRangeUser(-0.5, 1.5); // Set the z range
    h_A_C.Draw("colz");
    // h_A_C.GetXaxis()->SetTitle("Truth Bins");
    // h_A_C.GetYaxis()->SetTitle("Truth Bins");
    c_ac->Update();

    // Draw vertical and horizontal lines at the bin edges
    for (Int_t i = 1; i < n; i++) {
        TLine *vline = new TLine(bins[i], 0, bins[i], h_A_C.GetNbinsY());
        vline->SetLineColor(kBlack);
        vline->Draw();

        TLine *hline = new TLine(0, bins[i], h_A_C.GetNbinsX(), bins[i]);
        hline->SetLineColor(kBlack);
        hline->Draw();
    }

    for (Int_t i = 1; i <= n; i++) {
        // Draw white dotted lines from bins[i-1] to bins[i]
        TLine *vline_dotted1 = new TLine(bins[i], bins[i-1], bins[i], bins[i]);
        vline_dotted1->SetLineColor(kWhite);
        vline_dotted1->SetLineStyle(2); // Set line style to dotted
        vline_dotted1->Draw();

        TLine *hline_dotted1 = new TLine(bins[i-1], bins[i], bins[i], bins[i]);
        hline_dotted1->SetLineColor(kWhite);
        hline_dotted1->SetLineStyle(2); // Set line style to dotted
        hline_dotted1->Draw();

        if(i<n)
        {
            TLine *vline_dotted2 = new TLine(bins[i], bins[i+1], bins[i], bins[i]);
            vline_dotted2->SetLineColor(kWhite);
            vline_dotted2->SetLineStyle(2); // Set line style to dotted
            vline_dotted2->Draw();

            TLine *hline_dotted2 = new TLine(bins[i+1], bins[i], bins[i], bins[i]);
            hline_dotted2->SetLineColor(kWhite);
            hline_dotted2->SetLineStyle(2); // Set line style to dotted
            hline_dotted2->Draw();
        }
    }

    // Add labels in the middle of the intervals
    for (Int_t i = 0; i < n; i++) {
        Double_t midPoint = i==n-1 ? bins[i+1] + 1 : (bins[i] + bins[i+1]) / 2.0;
        TLatex *text = new TLatex(midPoint, 1.03*h_A_C.GetNbinsY(), labels[i]);
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->Draw();
        // delete text;
    }

    // Add asterisk for over-/underrflow bins
    for (Int_t i = 0; i < overUnderFlowN; i++) {
        Double_t midPoint = overUnderFlowBins[i] + 0.5;
        TLatex *text = new TLatex(midPoint, 1.0*h_A_C.GetNbinsY(), "*");
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->SetTextColor(kGray+3); // Set text color to grey
        text->Draw();
        // delete text;
    }

    // Add footnote
    TLatex *footnote3 = new TLatex(0, -0.1*h_A_C.GetNbinsY(), "* Under-/overflow bin; ** Selection subset");
    footnote3->SetTextSize(0.02); // Set text size to something smaller
    footnote3->SetTextColor(kGray+3);
    // footnote->SetTextAlign(22); // Center alignment
    footnote3->Draw();


    c_ac->SaveAs("plots/plot_entire_additional_smearing_matrix_nuwro.pdf");

    std::cout<<"DEBUG U3"<<std::endl;

    // std::cout << "A_C elements:\n";
    // for (int i = 0; i < A_C->GetNrows(); ++i) {
    //     for (int j = 0; j < A_C->GetNcols(); ++j) {
    //         std::cout << (*A_C)(i, j) << " ";
    //     }
    //     std::cout << std::endl;
    // }

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
        lg->AddEntry(fake_data_truth_hist, "NuWro Truth", "l");
    }

    lg->Draw("same");

    // Plot slices of the unfolded result
    auto *sb_ptr = new SliceBinning("../ubcc1pi_neutral_slice_config_reduced.txt");
    auto &sb = *sb_ptr;

    // Get the factors needed to convert to cross-section units
    double total_pot = syst.total_bnb_data_pot_;
    double integ_flux = integrated_numu_flux_in_FV(total_pot);
    double num_Ar = num_Ar_targets_in_FV();

    std::cout << "INTEGRATED numu FLUX = " << integ_flux << '\n';
    std::cout << "NUM Ar atoms in fiducial volume = " << num_Ar << '\n';

    // Retrieve the true-space expected event counts from NUISANCE output files
    // for each available generator model
    double conv_factor = (num_Ar * integ_flux) / 1e39;
    std::map<std::string, TMatrixD *> generator_truth_map = {}; // get_true_events_nuisance( sample_info, conv_factor );
    std::cout << "DEBUG test_unfolding_nuwro - Point 2" << std::endl;

    // Dump overall results to text files. Total cross section units (10^{-39}
    // cm^2 / Ar) will be used throughout. Do this before adjusting the
    // truth-level prediction TMatrixD objects via multiplication by A_C
    // dump_overall_results(result, unfolded_cov_matrix_map, 1.0 / conv_factor, genie_cv_truth_vec, fake_data_truth, generator_truth_map, using_fake_data);

    std::cout << "DEBUG test_unfolding_nuwro - Point 3" << std::endl;

    if (USE_ADD_SMEAR)
    {
        // Get access to the additional smearing matrix
        const TMatrixD &A_C = *result.add_smear_matrix_;

        // Start with the fake data truth if present
        if (using_fake_data)
            fake_data_truth = TMatrixD(A_C, TMatrixD::kMult, fake_data_truth);

        // Also transform the GENIE CV model
        genie_cv_truth_vec = TMatrixD(A_C, TMatrixD::kMult, genie_cv_truth_vec);

        // Now do the other generator predictions
        for (const auto &pair : generator_truth_map)
            *pair.second = TMatrixD(A_C, TMatrixD::kMult, *pair.second);
    }

    std::cout << "DEBUG test_unfolding_nuwro - Point 4" << std::endl;

    const auto cv_hist_2d = syst.get_cv_hist_2d();

    // Assuming cv_hist_2d is already defined and filled
    TCanvas* c_cv_hist_2d = new TCanvas("c_cv_hist_2d", "CV Histogram 2D", 800, 600);
    cv_hist_2d->Draw("colz");
    gPad->SetLogz();
    util::CreateWhiteToBlueColorPalette(20);
    c_cv_hist_2d->SaveAs("plots/entire_cv_hist_2d_nuwro.pdf");

    for (size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx)
    {
        // if (sl_idx != 4)
        //     continue; // TODO: remove !!!!!!!!!!!!!!!!!!!!!
        std::cout << "DEBUG test_unfolding_nuwro - Point 5" << std::endl;
        const auto &slice = sb.slices_.at(sl_idx);
        std::cout << "DEBUG test_unfolding_nuwro - Point 6" << std::endl;

        // Make a histogram showing the unfolded true event counts in the current slice
        SliceHistogram *slice_unf = SliceHistogram::make_slice_histogram(
            *result.unfolded_signal_, slice, result.cov_matrix_.get());

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

        std::cout<<"DEBUG truth from "<<start<<" to "<<stop<<std::endl;
        TMatrixD cv_hist_2d_slice(stop - start + 1, stop - start + 1);
        TMatrixD fake_hist_2d_slice(stop - start + 1, stop - start + 1);
        TMatrixD corr_unfolded_hist_slice(stop - start + 1, stop - start + 1);
        TMatrixD corr_hist_slice(stop - start + 1, stop - start + 1);
        for (int i = start; i <= stop; i++) {
            for (int j = start; j <= stop; j++) {
                cv_hist_2d_slice(i - start, j - start) = cv_hist_2d->GetBinContent(i+1, j+1);
                if(using_fake_data) fake_hist_2d_slice(i - start, j - start) = fake_data_univ->hist_2d_->GetBinContent(i+1, j+1);
                corr_unfolded_hist_slice(i - start, j - start) = corr_unfolded_hist->GetBinContent(i+1, j+1);
                corr_hist_slice(i - start, j - start) = corr_hist->GetBinContent(i+1, j+1);
            }
        }

        // Add over and underflow bins to additional smearing matrix plots
        if(sl_idx==2 || sl_idx==5) stop++; // Overflow bin
        if(sl_idx==5) start--; // Underflow bin
        TMatrixD ac_hist_slice(stop - start + 1, stop - start + 1);
        for (int i = start; i <= stop; i++) {
            for (int j = start; j <= stop; j++) {
                ac_hist_slice(i - start, j - start) = A_C->operator()(i, j);
            }
        }
        // We print out the dimensions of ac_hist_slice
        std::cout<<"DEBUG ac_hist_slice.GetNrows() = "<<ac_hist_slice.GetNrows()<<std::endl;

        const auto cv_confusion_mat = util::CountsToConfusionMatrix(cv_hist_2d_slice, "row");
        const auto fake_confusion_mat = util::CountsToConfusionMatrix(fake_hist_2d_slice, "row");
        std::string title = slice_unf->hist_->GetXaxis()->GetTitle();
        // std::regex pattern("true (.*)");
        // title = std::regex_replace(title, pattern, "$1");
        const auto conf_title = std::string("Confusion Matrix") + title +";Truth Bins;Reco Bins";
        const auto ac_title = std::string("Additional Smearing Matrix") + title +";Truth Bins;Reco Bins";
        const auto corr_unf_title = std::string("Correlation Matrix for Unfolded Bins") + title +";Truth Bins;Truth Bins";
        const auto corr_title = std::string("Correlation Matrix for Reco Bins") + title +";Reco Bins;Reco Bins";


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

        TH2D *trimmed_corr_hist = new TH2D(("trimmed_corr_hist slice "+std::to_string(sl_idx)).c_str(), corr_title.c_str(),
                                num_bins, 0, num_bins,
                                num_bins, 0, num_bins);

        // // Create a TH2D histogram
        // TH2D *cv_confusion_hist = new TH2D("cv_confusion",("Genie " + conf_title).c_str(),
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges,
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges);

        // // Create a TH2D histogram
        // TH2D *fake_confusion_hist = new TH2D("fake_confusion",("NuWro " + conf_title).c_str(),
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges,
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges);

        // // Create a TH2D histogram
        // TH2D *ac_hist = new TH2D("ac_hist",(ac_title).c_str(),
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges,
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges);
                                
        // TH2D *trimmed_corr_unf_hist = new TH2D("trimmed_corr_unf_hist",(corr_unf_title).c_str(),
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges,
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges);

        // TH2D *trimmed_corr_hist = new TH2D("trimmed_corr_hist", corr_title.c_str(),
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges,
        //                         sizeof(edges) / sizeof(Double_t) - 1, edges);

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

        for(int i = 0; i < corr_hist_slice.GetNrows(); i++) {
            for(int j = 0; j < corr_hist_slice.GetNcols(); j++) {
                trimmed_corr_hist->SetBinContent(i+1, j+1, corr_hist_slice.operator()(i, j));
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

        c_conf_1->SaveAs(("plots/cv_confusion_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro.pdf").c_str());


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

            c_conf_fake->SaveAs(("plots/fake_confusion_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro.pdf").c_str());
        }


        TCanvas* c_ac_slice = new TCanvas(("c_ac_slice slice "+std::to_string(sl_idx)).c_str(), "c ac matrix", 800, 600);
        ac_hist->SetTitleSize(0.05, "t");
        ac_hist->SetStats(false);
        ac_hist->GetZaxis()->SetRangeUser(-0.5, 1.5); // Set the z range
        ac_hist->Draw("colz");

        // Print out number of bins
        std::cout<<"DEBUG ac_hist.GetNbinsX() = "<<ac_hist->GetNbinsX()<<std::endl;

        // Set x and y axis labels to bin numbers
        for (int i = 1; i <= ac_hist->GetNbinsX(); i++) {
            const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 1 || i == 6));
            const auto binLabel = overFlow ? Form("%d*", i) : Form("%d", i);
            ac_hist->GetXaxis()->SetBinLabel(i, binLabel);
        }
        for (int i = 1; i <= ac_hist->GetNbinsY(); i++) {
            const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 1 || i == 6));
            const auto binLabel = overFlow ? Form("%d*", i) : Form("%d", i);
            ac_hist->GetYaxis()->SetBinLabel(i, binLabel);
        }
        // Increase the font size of the axis labels
        ac_hist->GetXaxis()->SetLabelSize(0.05);
        ac_hist->GetYaxis()->SetLabelSize(0.05);

        TColor::CreateGradientColorTable(nColors, stops, red, green, blue, 20);
        c_ac_slice->SaveAs(("plots/additional_smearing_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro.pdf").c_str());

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

        c_corr_unf_slice->SaveAs(("plots/unfolded_correlation_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro.pdf").c_str());

        // Plot for trimmed_corr_hist
        TCanvas* c_corr_slice = new TCanvas(("c_corr_slice slice "+std::to_string(sl_idx)).c_str(), "c correlation matrix", 800, 600);
        util::CreateRedToBlueColorPalette(20);
        trimmed_corr_hist->SetTitleSize(0.05, "t");
        trimmed_corr_hist->GetZaxis()->SetRangeUser(-1, 1);
        trimmed_corr_hist->SetStats(false);
        trimmed_corr_hist->GetZaxis()->SetRangeUser(-1.0, 1.0); // Set the z range
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

        c_corr_slice->SaveAs(("plots/correlation_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro.pdf").c_str());


        // TMatrixD* matrix = result.cov_matrix_.get();
        // std::cout<<"DEBUG result.cov_matrix_: "<<std::endl;
        // for (int i = 0; i < matrix->GetNrows(); i++) {
        //     for (int j = 0; j < matrix->GetNcols(); j++) {
        //         std::cout << (*matrix)(i, j) << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // TH2D* hist2D = slice_unf->cmat_.cov_matrix_.get();
        // std::cout<<"DEBUG slice_unf->cmat_.cov_matrix_: "<<std::endl;
        // for (int i = 1; i <= hist2D->GetNbinsX(); i++) {
        //     for (int j = 1; j <= hist2D->GetNbinsY(); j++) {
        //         std::cout << hist2D->GetBinContent(i, j) << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // const auto matrix2 = slice_unf->cmat_.get_matrix();
        // std::cout<<"DEBUG slice_unf->cmat_.get_matrix: "<<std::endl;
        // for (int i = 0; i < matrix2->GetNrows(); i++) {
        //     for (int j = 0; j < matrix2->GetNcols(); j++) {
        //         std::cout << (*matrix2)(i, j) << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // Temporary copies of the unfolded true event count slices with different covariance matrices
        // std::map<std::string, std::unique_ptr<SliceHistogram>> sh_cov_map;
        // for (const auto &uc_pair : unfolded_cov_matrix_map)
        // {
        //     sh_cov_map[uc_pair.first].reset(
        //         SliceHistogram::make_slice_histogram(*result.unfolded_signal_, slice,
        //                                              uc_pair.second.get()));
        // }

        // Also use the GENIE CV model to do the same
        SliceHistogram *slice_cv = SliceHistogram::make_slice_histogram(
            genie_cv_truth_vec, slice, nullptr);

        // If present, also use the truth information from the fake data to do the same
        SliceHistogram *slice_truth = using_fake_data ? SliceHistogram::make_slice_histogram(fake_data_truth, slice, nullptr) : nullptr;

        // Keys are legend labels, values are SliceHistogram objects containing true-space predictions from the corresponding generator models
        auto *slice_gen_map_ptr = new std::map<std::string, SliceHistogram *>();
        auto &slice_gen_map = *slice_gen_map_ptr;

        slice_gen_map["Unfolded NuWro Reco"] = slice_unf;
        slice_gen_map["MicroBooNE Tune"] = slice_cv;

        if (using_fake_data)
            slice_gen_map["NuWro Truth"] = slice_truth;

        for (const auto &pair : generator_truth_map)
            slice_gen_map[pair.first] = SliceHistogram::make_slice_histogram(*pair.second, slice, nullptr);

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

        if(sl_idx == 8) var_count = 0; // TODO: Find a better way to do this

        std::cout << "DEBUG test_unfolding_nuwro - Point 7" << std::endl;

        int num_slice_bins = slice_unf->hist_->GetNbinsX();
        TMatrixD trans_mat(num_slice_bins, num_slice_bins);
        for (int b = 0; b < num_slice_bins; ++b)
        {
            double width = slice_unf->hist_->GetBinWidth(b + 1) * other_var_width;
            trans_mat(b, b) = 1e39 / (width * integ_flux * num_Ar);
            std::cout << "DEBUG at " << b << " trans_mat: " << trans_mat(b, b) << " width: " << width << " integ_flux: " << integ_flux << " num_Ar: " << num_Ar << std::endl;
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

        // for (auto &sh_cov_pair : sh_cov_map)
        // {
        //     auto &slice_h = sh_cov_pair.second;
        //     slice_h->transform(trans_mat);
        // }

        std::map<std::string, SliceHistogram::Chi2Result> chi2_map;
        std::cout << '\n';

        for (const auto &pair : slice_gen_map)
        {
            const auto &name = pair.first;
            const auto *slice_h = pair.second;

            if (name == "Unfolded NuWro Reco")
            {
                chi2_map[name] = SliceHistogram::Chi2Result();
                continue;
            }

            const auto &chi2_result = chi2_map[name] = slice_h->get_chi2(*slice_gen_map.at("Unfolded NuWro Reco"));

            double chi2_alt = slice_h->hist_->Chi2Test(slice_gen_map.at("Unfolded NuWro Reco")->hist_.get(), "UW CHI2");
            std::cout<<"DEBUG (slice "<<sl_idx<<") Chi-square alt: "<<chi2_alt<<std::endl;


            std::cout << "Slice " << sl_idx << ", " << name << ": \u03C7\u00b2 = "
                      << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bin"
                      << (chi2_result.num_bins_ > 1 ? "s" : "") << ", p-value = " << chi2_result.p_value_ << '\n';
        }

        std::cout << "DEBUG Unfolding Point 4" << std::endl;











        TCanvas *c1 = new TCanvas;

        TPad* pad2 = new TPad(("pad2 slice "+std::to_string(sl_idx)).c_str(), "", 0, 0.0, 1.0, 0.3);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();
        pad2->cd();

        // Create the ratio histogram
        TH1D* h_ratio = dynamic_cast<TH1D*>(slice_unf->hist_.get()->Clone("h_ratio"));
        h_ratio->SetStats(false);

        // Create a temporary clone of slice_cv->hist_ with zero errors
        TH1D* h_cv_no_errors = dynamic_cast<TH1D*>(slice_cv->hist_.get()->Clone("h_cv_no_errors"));
        for (int i = 1; i <= h_cv_no_errors->GetNbinsX(); ++i) {
            h_cv_no_errors->SetBinError(i, 0.0);
        }

        double chi2_1 = slice_unf->hist_->Chi2Test(slice_cv->hist_.get(), "WU CHI2");
        double chi2_2 = slice_unf->hist_->Chi2Test(slice_cv->hist_.get(), "WW CHI2");
        std::cout << "DEBUG (slice "<<sl_idx<<") Chi-square unf 1: " << chi2_1 << " Chi-square unf 2: "<<chi2_2<<std::endl;
        std::cout<<"DEBUG unfolded vs genie truth: ";
        for (int i = 1; i <= slice_unf->hist_->GetNbinsX(); ++i) {
            double value = slice_unf->hist_->GetBinContent(i);
            double error = slice_unf->hist_->GetBinError(i);
            double genieTruthValue = slice_cv->hist_->GetBinContent(i);
            std::cout << value << " +- " << error << "("<<genieTruthValue<<")" << ", ";
        }
        std::cout<<std::endl;

        // Draw a horizontal dashed line at ratio == 1
        TLine* line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
        line->SetLineColor(kAzure - 7);
        line->SetLineWidth(2);
        line->SetLineStyle(5);

        h_ratio->Divide(h_cv_no_errors);
        h_ratio->SetLineWidth(2);
        h_ratio->SetLineColor(kBlack);
        h_ratio->SetMarkerStyle(kFullCircle);
        h_ratio->SetMarkerSize(0.8);
        h_ratio->SetTitle("");

        // x-axis
        h_ratio->GetXaxis()->SetTitle(slice_unf->hist_->GetXaxis()->GetTitle());
        // h_ratio->GetXaxis()->CenterTitle(true);
        h_ratio->GetXaxis()->SetLabelSize(0.1);
        h_ratio->GetXaxis()->SetTitleSize(0.08);
        h_ratio->GetXaxis()->SetTickLength(0.05);
        // h_ratio->GetXaxis()->SetTitleOffset(0.9);

        // y-axis
        h_ratio->GetYaxis()->SetTitle("ratio");
        h_ratio->GetYaxis()->CenterTitle(true);
        h_ratio->GetYaxis()->SetLabelSize(0.08);
        h_ratio->GetYaxis()->SetTitleSize(0.15);
        // h_ratio->GetYaxis()->SetTitleOffset(0.35);

        // Set y-axis range
        double min = 0.9 * (h_ratio->GetBinContent(h_ratio->GetMinimumBin()) - h_ratio->GetBinError(h_ratio->GetMinimumBin()));
        double max = 1.1 * (h_ratio->GetBinContent(h_ratio->GetMaximumBin()) + h_ratio->GetBinError(h_ratio->GetMaximumBin()));

        h_ratio->GetYaxis()->SetRangeUser(min, max);
        h_ratio->Draw("e");

        if (using_fake_data)
        {
            TH1D* h_ratio_truth = dynamic_cast<TH1D*>(slice_truth->hist_.get()->Clone("h_ratio_truth"));
            double chi2_1 = slice_unf->hist_->Chi2Test(slice_truth->hist_.get(), "WU CHI2");
            double chi2_2 = slice_unf->hist_->Chi2Test(slice_truth->hist_.get(), "WW CHI2");
            std::cout << "DEBUG (slice "<<sl_idx<<") Chi-square fake truth 1: " << chi2_1 << " Chi-square fake truth 2: "<<chi2_2<<std::endl;
            std::cout<<"DEBUG unfolded vs nuwro truth: ";
            for (int i = 1; i <= slice_unf->hist_->GetNbinsX(); ++i) {
                double value = slice_unf->hist_->GetBinContent(i);
                double error = slice_unf->hist_->GetBinError(i);
                double nuwroTruthValue = slice_truth->hist_->GetBinContent(i);
                std::cout << value << " +- " << error << "("<<nuwroTruthValue<<")" << ", ";
            }
            std::cout<<std::endl;
            h_ratio_truth->SetStats(false);
            h_ratio_truth->Divide(h_cv_no_errors);
            h_ratio_truth->SetLineColor(kGreen);
            h_ratio_truth->SetLineWidth(2);
            min = TMath::Min(min, 0.9 * h_ratio_truth->GetMinimum());
            max = TMath::Max(max, 1.1 * h_ratio_truth->GetMaximum());
            h_ratio_truth->GetYaxis()->SetRangeUser(min, max);
            h_ratio_truth->Draw("hist same");
        }

        // TH1D* statUncertHist = dynamic_cast<TH1D*>(sh_cov_map.at("DataStats")->hist_.get());
        // TH1D* statUncertHist_ratio = dynamic_cast<TH1D*>(statUncertHist->Clone("statUncertHist_ratio"));
        // statUncertHist_ratio->Divide(h_cv_no_errors);
        // statUncertHist_ratio->SetLineWidth(2);
        // statUncertHist_ratio->SetLineColor(kBlack);
        // statUncertHist_ratio->Draw("E1 same");

        line->Draw();


        // Go back to the main canvas before creating the second pad
        c1->cd();

        TPad* pad1 = new TPad(("pad1 slice "+std::to_string(sl_idx)).c_str(), "", 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.03);
        pad1->SetTopMargin(0.1);
        pad1->Draw();
        pad1->cd();


        slice_unf->hist_->SetLineColor(kBlack);
        // slice_unf->hist_->SetLineWidth(5);
        slice_unf->hist_->SetMarkerStyle(kFullCircle);
        slice_unf->hist_->SetMarkerSize(0.8);
        slice_unf->hist_->SetStats(false);
        slice_unf->hist_->SetLineWidth(2);

        std::cout << "DEBUG Unfolding Point 5" << std::endl;

        double ymax = -DBL_MAX;
        slice_unf->hist_->Draw("e");

        for (const auto &pair : slice_gen_map)
        {
            const auto &name = pair.first;
            const auto *slice_h = pair.second;

            double max = slice_h->hist_->GetMaximum();
            if (max > ymax)
                ymax = max;

            if (name == "Unfolded NuWro Reco" || name == "NuWro Truth" || name == "MicroBooNE Tune")
                continue;

            const auto &file_info = truth_file_map.at(name);
            slice_h->hist_->SetLineColor(file_info.color_);
            slice_h->hist_->SetLineStyle(file_info.style_);
            slice_h->hist_->SetLineWidth(3);
            slice_h->hist_->Draw("hist same");
        }

        slice_cv->hist_->SetStats(false);
        slice_cv->hist_->SetLineColor(kAzure - 7);
        slice_cv->hist_->SetLineWidth(2);
        slice_cv->hist_->SetLineStyle(5);
        slice_cv->hist_->Draw("hist same");

        if (using_fake_data)
        {
            slice_truth->hist_->SetStats(false);
            slice_truth->hist_->SetLineColor(kGreen);
            slice_truth->hist_->SetLineWidth(2);
            slice_truth->hist_->Draw("hist same");
        }

        slice_unf->hist_->GetYaxis()->SetRangeUser(0., ymax * 1.07);
        slice_unf->hist_->Draw("3 same");
        slice_unf->hist_->SetTitle("Unfolded NuWro CC1#pi^{#pm}Np");
        slice_unf->hist_->GetXaxis()->SetLabelOffset(999); // Hide x-axis
        slice_unf->hist_->GetXaxis()->SetTitleOffset(999); // Hide x-axis title

        // auto graph1_sys = new TGraphErrors(5, x, py1, zero, ey_sys1);
        // graph1_sys->Draw("[]");

        // for (const auto& key : {std::string("total_blockwise_norm"), std::string("total_blockwise_shape"), std::string("total_blockwise_mixed")}) {
        //     TH1D* hist = dynamic_cast<TH1D*>(sh_cov_map.at(key)->hist_.get());
        //     // Draw the histogram
        //     if (key == std::string("total_blockwise_norm")) {
        //         hist->SetFillColor(kBlue);
        //         hist->SetLineWidth(2); // Set line width for this histogram
        //     } else if (key == std::string("total_blockwise_shape")) {
        //         hist->SetFillColor(kRed);
        //         hist->SetLineWidth(3); // Set a different line width for this histogram
        //     } else if (key == std::string("total_blockwise_mixed")) {
        //         hist->SetFillColor(kGreen);
        //         hist->SetLineWidth(4); // Set yet another line width for this histogram
        //     }
        //     hist->Draw("e same");
        // }
        
        // statUncertHist->SetLineWidth(2);
        // statUncertHist->SetLineColor(kBlack);
        // statUncertHist->Draw("E1 same");


        std::cout << "DEBUG Unfolding Point 6" << std::endl;

        TLegend *lg = new TLegend(0.6, 0.2);
        for (const auto &pair : slice_gen_map)
        {
            const auto &name = pair.first;
            const auto *slice_h = pair.second;

            std::string label = name;
            const auto &chi2_result = chi2_map.at(name);
            std::ostringstream oss;
            oss << std::setprecision(3) << chi2_result.chi2_ << " / " << chi2_result.num_bins_ << " bin";
            if (chi2_result.num_bins_ > 1)
                oss << "s; p-value = " << chi2_result.p_value_;

            label += (name != "Unfolded NuWro Reco") ? ": #chi^{2} = " + oss.str() : "";

            lg->AddEntry(slice_h->hist_.get(), label.c_str(), "l");
        }

        lg->Draw("same");

        c1->Update();
        std::cout << "DEBUG test_unfolding_nuwro - Point 8" << std::endl;

        std::string out_pdf_name = "plots/plot_unfolded_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro.pdf";
        c1->SaveAs(out_pdf_name.c_str());



        // // Assuming slice.bin_map_ is a std::map<int, std::set<size_t>>
        // size_t indexXMin = std::numeric_limits<size_t>::max();
        // size_t indexXMax = 0;

        // for (const auto &pair : slice.bin_map_) { // ATTN: This only works correctly for a continuous set of individual bins
        //     if (!pair.second.empty()) {
        //         size_t min_val = *pair.second.begin();
        //         size_t max_val = *pair.second.rbegin();

        //         if (min_val < indexXMin) {
        //             indexXMin = min_val;
        //         }
        //         if (max_val > indexXMax) {
        //             indexXMax = max_val;
        //         }
        //     }
        // }

        // Now indexXMin and indexXMax hold the overall smallest and largest values respectively
        // You can use these variables in your code

        // indexXMin++;
        // indexXMax++;
        // const auto indexYMin = indexXMin;
        // const auto indexYMax = indexXMax;

        // // Assuming cov_mat is a std::unique_ptr<TH2D>
        // // Create a new histogram with the dimensions of the subset
        // TH2D *cov_matrix = new TH2D("cov_matrix", "Subset of original histogram", 
        //                         indexXMax-indexXMin+1, indexXMin-1, indexXMax, 
        //                         indexYMax-indexYMin+1, indexYMin-1, indexYMax);

        // // Fill the new histogram with the values from the original histogram
        // for (int i = indexXMin; i <= indexXMax; ++i) {
        //     for (int j = indexYMin; j <= indexYMax; ++j) {
        //         double bin_content = cov_mat->GetBinContent(i, j);
        //         cov_matrix->SetBinContent(i-indexXMin+1, j-indexYMin+1, bin_content);
        //     }
        // }

        // // const auto cov_matrix = slice_hist->cmat_.cov_matrix_.get();
        // TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
        // cov_matrix->Draw("COLZ");
        // cov_matrix->SetTitle("Covariance matrix");
        // std::string out_pdf_name_cov = "plots/plot_slice_cov_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + ".pdf";
        // c2->SaveAs(out_pdf_name_cov.c_str());

        // const auto unfolded_cov_matrix = slice_gen_map.at("Unfolded NuWro Reco")->cmat_.cov_matrix_.get();
        // TCanvas *c3 = new TCanvas("c3", "Canvas", 800, 600);
        // unfolded_cov_matrix->Draw("COLZ");
        // unfolded_cov_matrix->SetTitle("Unfolded covariance matrix");
        // std::string out_pdf_name_unf_cov = "plots/plot_unfolded_slice_cov_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + ".pdf";
        // c3->SaveAs(out_pdf_name_unf_cov.c_str());
        // delete cov_matrix;


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

        std::cout << "DEBUG Unfolding Point 7" << std::endl;

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

        std::cout << "DEBUG test_unfolding_nuwro - Point 9" << std::endl;
        write_pgfplots_files(output_file_prefix, slice_hist_table, slice_params_table);
        std::cout << "DEBUG test_unfolding_nuwro - Point 10" << std::endl;
    } // slices
    std::cout << "DEBUG test_unfolding_nuwro - Point 11" << std::endl;
    return 0;
}

int main()
{
    std::cout << "DEBUG Unfolding Point main 0" << std::endl;
    test_unfolding_nuwro();
    std::cout << "DEBUG Unfolding Point main 1" << std::endl;
    return 0;
}
