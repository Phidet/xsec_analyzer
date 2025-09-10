// Standard library includes
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

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
  TMatrixD hist_mat_transformed(*mat, TMatrixD::EMatrixCreatorsOp2::kMult, hist_mat);

  // Update the input histogram contents with the new values
  for (int r = 0; r < num_bins; ++r)
  {
    double val = hist_mat_transformed(r, 0);
    hist->SetBinContent(r + 1, val);
  }
}

// const std::string SAMPLE_NAME = "MicroBooNE_CC1MuXp_XSec_2D_PpCosp_nu_MC";

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
    {"unfolded data", "UnfData"},
    {"MicroBooNE Tune", "uBTune"},
    {"truth", "FakeData"}};

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
// counts to a total cross section (10^{-38} cm^2 / Ar) via multiplication.
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

void test_unfolding_simpler()
{

  //// Initialize the FilePropertiesManager and tell it to treat the NuWro
  //// MC ntuples as if they were data
  auto &fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties("../nuwro_file_properties.txt");

  //   const auto& sample_info = sample_info_map.at( SAMPLE_NAME );
  // const auto& respmat_file_name = sample_info.respmat_file_;

  const std::string respmat_file_name(
      "/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_4/univmake_output_NuWro.root");

  // Do the systematics calculations in preparation for unfolding
  auto *syst_ptr = new MCC9SystematicsCalculator(respmat_file_name, "../systcalc_unfold_fd_modified.conf");
  // auto* syst_ptr = new MCC9SystematicsCalculator( respmat_file_name, "../systcalc.conf" );
  auto &syst = *syst_ptr;

  // Get the tuned GENIE CV prediction in each true bin (including the
  // background true bins)
  TH1D *genie_cv_truth = syst.cv_universe().hist_true_.get();
  int num_true_bins = genie_cv_truth->GetNbinsX();

  // While we're at it, clone the histogram and zero it out. We'll fill this
  // one with our unfolded result for easy comparison
  TH1D *unfolded_events = dynamic_cast<TH1D *>(
      genie_cv_truth->Clone("unfolded_events"));
  unfolded_events->Reset();

  // If present, then get the fake data event counts in each true bin
  // (including the background true bins). We hope to approximately reproduce
  // these event counts in the signal true bins via unfolding the fake data.
  const auto &fake_data_univ = syst.fake_data_universe();
  TH1D *fake_data_truth_hist = nullptr;

  bool using_fake_data = false;
  if (fake_data_univ)
  {
    using_fake_data = true;
    fake_data_truth_hist = fake_data_univ->hist_true_.get();
  }

  std::cout << "\n\n############\nDEBUG: using_fake_data is " << (using_fake_data ? "true" : "false") << "\n############\n\n"
            << std::endl;

  int num_ordinary_reco_bins = 0;
  int num_sideband_reco_bins = 0;
  for (int b = 0; b < syst.reco_bins_.size(); ++b)
  {
    const auto &rbin = syst.reco_bins_.at(b);
    if (rbin.type_ == kSidebandRecoBin)
      ++num_sideband_reco_bins;
    else
      ++num_ordinary_reco_bins;
  }

  int num_true_signal_bins = 0;
  for (int t = 0; t < syst.true_bins_.size(); ++t)
  {
    const auto &tbin = syst.true_bins_.at(t);
    if (tbin.type_ == kSignalTrueBin)
      ++num_true_signal_bins;
  }

  std::cout << "NUM ORDINARY RECO BINS = " << num_ordinary_reco_bins << '\n';
  std::cout << "NUM TRUE SIGNAL BINS = " << num_true_signal_bins << '\n';

  auto *matrix_map_ptr = syst.get_covariances().release();
  auto &matrix_map = *matrix_map_ptr;

  auto *cov_mat = matrix_map.at("total").cov_matrix_.get();

  constexpr int NUM_DAGOSTINI_ITERATIONS = 2;
  constexpr bool USE_ADD_SMEAR = true;

  std::unique_ptr<Unfolder> unfolder(
      // new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS )
      new DAgostiniUnfolder(DAgostiniUnfolder::ConvergenceCriterion ::FigureOfMerit, 0.025)
      // new WienerSVDUnfolder( true,
      //  WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv )
  );

  std::cout << "DEBUG Pre unfolding" << std::endl;
  UnfoldedMeasurement result = unfolder->unfold(syst);
  std::cout << "DEBUG Post unfolding" << std::endl;

  // For real data only, add some new covariance matrices in which only the
  // signal response or the background is varied. We could calculate these
  // for the fake data, but it seems unnecessary at this point.
  if (!using_fake_data)
  {
    syst.set_syst_mode(MCC9SystematicsCalculator ::SystMode::VaryOnlyBackground);
    auto *bkgd_matrix_map_ptr = syst.get_covariances().release();
    auto &bkgd_matrix_map = *bkgd_matrix_map_ptr;

    syst.set_syst_mode(MCC9SystematicsCalculator ::SystMode::VaryOnlySignalResponse);
    auto *sigresp_matrix_map_ptr = syst.get_covariances().release();
    auto &sigresp_matrix_map = *sigresp_matrix_map_ptr;

    for (const auto &m_pair : bkgd_matrix_map)
    {
      auto &my_temp_cov_mat = matrix_map["bkgd_only_" + m_pair.first];
      my_temp_cov_mat += m_pair.second;
    }

    for (const auto &m_pair : sigresp_matrix_map)
    {
      auto &my_temp_cov_mat = matrix_map["sigresp_only_" + m_pair.first];
      my_temp_cov_mat += m_pair.second;
    }
  }

  // Propagate all defined covariance matrices through the unfolding procedure
  const TMatrixD &err_prop = *result.err_prop_matrix_;
  TMatrixD err_prop_tr(TMatrixD::kTransposed, err_prop);

  std::map<std::string, std::unique_ptr<TMatrixD>> unfolded_cov_matrix_map;

  for (const auto &matrix_pair : matrix_map)
  {
    const std::string &matrix_key = matrix_pair.first;
    auto temp_cov_mat = matrix_pair.second.get_matrix();

    TMatrixD temp_mat(*temp_cov_mat, TMatrixD::EMatrixCreatorsOp2::kMult,
                      err_prop_tr);

    unfolded_cov_matrix_map[matrix_key] = std::make_unique<TMatrixD>(
        err_prop, TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat);
  }

  // Decompose the block-diagonal pieces of the total covariance matrix
  // into normalization, shape, and mixed components (for later plotting
  // purposes)
  NormShapeCovMatrix bd_ns_covmat = make_block_diagonal_norm_shape_covmat(
      *result.unfolded_signal_, *result.cov_matrix_, syst.true_bins_);

  // Add the blockwise decomposed matrices into the map
  unfolded_cov_matrix_map["total_blockwise_norm"] = std::make_unique<TMatrixD>(bd_ns_covmat.norm_);

  unfolded_cov_matrix_map["total_blockwise_shape"] = std::make_unique<TMatrixD>(bd_ns_covmat.shape_);

  unfolded_cov_matrix_map["total_blockwise_mixed"] = std::make_unique<TMatrixD>(bd_ns_covmat.mixed_);

  //// Test against RooUnfold implementation of the D'Agostini method
  // auto test_response = get_test_response( mcc9 );
  // RooUnfoldBayes roounfold_unfolder( test_response.get(),
  //   test_response->Hmeasured(), NUM_DAGOSTINI_ITERATIONS );
  // auto measured = mcc9.get_measured_events();
  // roounfold_unfolder.SetMeasuredCov( *measured.cov_matrix_ );
  // roounfold_unfolder.IncludeSystematics( 1 );

  // auto* roo_hist = roounfold_unfolder.Hreco(
  //   RooUnfold::ErrorTreatment::kCovariance );

  // auto roo_cov = roounfold_unfolder.Ereco(
  //   RooUnfold::ErrorTreatment::kCovariance );

  // for ( int t = 0; t < num_true_signal_bins; ++t ) {
  //   double mine = result.unfolded_signal_->operator()( t, 0 );
  //   double my_err = std::sqrt( std::max(0.,
  //     result.cov_matrix_->operator()( t, t )) );
  //   double theirs = roo_hist->GetBinContent( t + 1 );
  //   double their_err = std::sqrt( std::max(0., roo_cov(t, t)) );

  //  std::cout << "t = " << t << ' ' << mine << " ± " << my_err
  //    << "  " << theirs << " ± " << their_err << "  "
  //    << (mine - theirs) / theirs << " ± "
  //    << (my_err - their_err) / their_err << '\n';
  //}

  // auto smearcept = mcc9.get_cv_smearceptance_matrix();
  // auto true_signal = syst.get_cv_true_signal();

  //// Sanity check (these were equal to each other as expected)
  //// int num_reco_bins = smearcept->GetNrows();
  ////TMatrixD reco_signal( *smearcept, TMatrixD::kMult, *true_signal );
  ////auto data_sig = mcc9.get_cv_ordinary_reco_signal();
  ////for ( int r = 0; r < num_reco_bins; ++r ) {
  ////  double r1 = reco_signal( r, 0 );
  ////  double r2 = data_sig->operator()( r, 0 );
  ////  std::cout << "r = " << r << ", r1 = " << r1 << ", r2 = " << r2
  ////    << ", diff = " << r1 - r2 << '\n';
  ////}

  // auto meas = mcc9.get_measured_events();
  //// ASIMOV TEST: USE CV RECO SIGNAL PREDICTION INSTEAD OF
  //// BACKGROUND-SUBTRACTED DATA
  ////const auto& data_signal = meas.reco_signal_;
  // auto data_signal = mcc9.get_cv_ordinary_reco_signal();
  // const auto& data_covmat = meas.cov_matrix_;

  // auto result = unfolder->unfold( *data_signal, *data_covmat,
  //   *smearcept, *true_signal );

  ////// Test with original implementation
  ////TVectorD true_sig_vec( num_true_bins );
  ////for ( int t = 0; t < num_true_bins; ++t ) {
  ////  true_sig_vec( t ) = true_signal->operator()( t, 0 );
  ////}

  ////TVectorD data_sig_vec( num_ordinary_reco_bins );
  ////for ( int r = 0; r < num_ordinary_reco_bins; ++r ) {
  ////  data_sig_vec( r ) = data_signal->operator()( r, 0 );
  ////}

  ////TMatrixD A_C_svd( num_true_bins, num_true_bins );
  ////TVectorD WF_svd( num_true_bins );
  ////TMatrixD UnfoldCov_svd( num_true_bins, num_true_bins );

  ////TVectorD wsvd_unfolded_signal = WienerSVD( *smearcept, true_sig_vec,
  ////  data_sig_vec, *data_covmat, 0, 0, A_C_svd, WF_svd, UnfoldCov_svd, 0. );

  // for ( int t = 0; t < num_true_bins; ++t ) {
  //   double t_evts = true_signal->operator()( t, 0 );
  //   double evts = result.unfolded_signal_->operator()( t, 0 );
  //   double error = std::sqrt( std::max(0.,
  //     result.cov_matrix_->operator()( t, t )) );
  //   //double evts_svd = wsvd_unfolded_signal( t );
  //   //double error_svd = std::sqrt( std::max(0.,
  //   //  UnfoldCov_svd( t, t )) );
  //   std::cout << "t = " << t << ", t_evts = " << t_evts
  //     << ", evts = " << evts << " ± " << error
  //     << ", diff = " << t_evts - evts << '\n';
  // }

  // Set the event counts in each bin of the histogram that displays the
  // unfolded result. Note that we don't care about the background true bins
  // (which are assumed to follow all of the signal true bins) since we've
  // subtracted out an estimate of the background before unfolding.
  for (int t = 0; t < num_true_bins; ++t)
  {
    double evts = 0.;
    double error = 0.;
    if (t < num_true_signal_bins)
    {
      evts = result.unfolded_signal_->operator()(t, 0);
      error = std::sqrt(std::max(0., result.cov_matrix_->operator()(t, t)));
    }

    // We need to use one-based indices while working with TH1D bins
    unfolded_events->SetBinContent(t + 1, evts);
    unfolded_events->SetBinError(t + 1, error);
  }

  unfolded_events->SetStats(false);
  unfolded_events->SetLineColor(kBlack);
  unfolded_events->SetLineWidth(3);
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

  // Multiply the truth-level GENIE prediction histogram by the additional
  // smearing matrix
  TMatrixD *A_C = result.add_smear_matrix_.get();
  multiply_1d_hist_by_matrix(A_C, genie_cv_truth);

  genie_cv_truth->SetStats(false);
  genie_cv_truth->SetLineColor(kRed);
  genie_cv_truth->SetLineWidth(3);
  genie_cv_truth->SetLineStyle(9);

  unfolded_events->Draw("e");
  genie_cv_truth->Draw("hist same");

  if (using_fake_data)
  {

    // Multiply the fake data truth histogram by the additional smearing matrix
    multiply_1d_hist_by_matrix(A_C, fake_data_truth_hist);

    fake_data_truth_hist->SetStats(false);
    fake_data_truth_hist->SetLineColor(kBlue);
    fake_data_truth_hist->SetLineWidth(3);
    fake_data_truth_hist->SetLineStyle(2);
    fake_data_truth_hist->Draw("hist same");
  }

  // TLegend* lg = new TLegend( 0.15, 0.7, 0.3, 0.85 ); // Top right
  // TLegend* lg = new TLegend( 0.7, 0.7, 0.85, 0.85 ); // Top left
  // TLegend* lg = new TLegend( 0.7, 0.15, 0.85, 0.3 ); // Bottom left
  // TLegend* lg = new TLegend( 0.15, 0.15, 0.3, 0.3 ); // Bottom right
  TLegend *lg = new TLegend(0.5, 0.2);
  lg->AddEntry(unfolded_events, "unfolded", "l");
  lg->AddEntry(genie_cv_truth, "uB tune", "l");
  if (using_fake_data)
  {
    lg->AddEntry(fake_data_truth_hist, "truth", "l");
  }

  lg->Draw("same");

  // Plot slices of the unfolded result
  auto *sb_ptr = new SliceBinning("../ubcc1pi_slice_config.txt");
  auto &sb = *sb_ptr;

  // Get the factors needed to convert to cross-section units
  double total_pot = syst.total_bnb_data_pot_;
  double integ_flux = integrated_numu_flux_in_FV(total_pot);
  double num_Ar = num_Ar_targets_in_FV();

  std::cout << "INTEGRATED numu FLUX = " << integ_flux << '\n';
  std::cout << "NUM Ar atoms in fiducial volume = " << num_Ar << '\n';

  // Retrieve the true-space expected event counts from NUISANCE output files
  // for each available generator model
  double conv_factor = (num_Ar * integ_flux) / 1e38;
  std::map<std::string, TMatrixD *> generator_truth_map = {}; // get_true_events_nuisance( sample_info, conv_factor );

  // Dump overall results to text files. Total cross section units (10^{-38}
  // cm^2 / Ar) will be used throughout. Do this before adjusting the
  // truth-level prediction TMatrixD objects via multiplication by A_C
  dump_overall_results(result, unfolded_cov_matrix_map, 1.0 / conv_factor,
                       genie_cv_truth_vec, fake_data_truth, generator_truth_map,
                       using_fake_data);

  if (USE_ADD_SMEAR)
  {

    // Get access to the additional smearing matrix
    const TMatrixD &A_C = *result.add_smear_matrix_;

    // Start with the fake data truth if present
    if (using_fake_data)
    {
      TMatrixD ac_truth(A_C, TMatrixD::kMult, fake_data_truth);
      fake_data_truth = ac_truth;
    }

    // Also transform the GENIE CV model
    TMatrixD genie_cv_temp(A_C, TMatrixD::kMult, genie_cv_truth_vec);
    genie_cv_truth_vec = genie_cv_temp;

    // Now do the other generator predictions
    for (const auto &pair : generator_truth_map)
    {
      const auto &model_name = pair.first;
      TMatrixD *truth_mat = pair.second;

      TMatrixD ac_temp(A_C, TMatrixD::kMult, *truth_mat);
      *truth_mat = ac_temp;
    }
  }

  for (size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx)
  {

    // if( sl_idx != 0 ) continue;

    const auto &slice = sb.slices_.at(sl_idx);

    // Make a histogram showing the unfolded true event counts in the current
    // slice
    SliceHistogram *slice_unf = SliceHistogram::make_slice_histogram(
        *result.unfolded_signal_, slice, result.cov_matrix_.get());
    std::cout << "DEBUG at " << sl_idx << " result.cov_matrix_.get() is " << (result.cov_matrix_.get() == nullptr ? "" : "not") << " nullptr" << std::endl;

    // Temporary copies of the unfolded true event count slices with
    // different covariance matrices
    std::map<std::string, std::unique_ptr<SliceHistogram>> sh_cov_map;
    for (const auto &uc_pair : unfolded_cov_matrix_map)
    {
      const auto &uc_name = uc_pair.first;
      const auto &uc_matrix = uc_pair.second;

      auto &uc_ptr = sh_cov_map[uc_name];
      uc_ptr.reset(
          SliceHistogram::make_slice_histogram(*result.unfolded_signal_, slice,
                                               uc_matrix.get()));
    }

    // Also use the GENIE CV model to do the same
    SliceHistogram *slice_cv = SliceHistogram::make_slice_histogram(
        genie_cv_truth_vec, slice, nullptr);

    // If present, also use the truth information from the fake data to do the
    // same
    SliceHistogram *slice_truth = nullptr;
    if (using_fake_data)
    {
      slice_truth = SliceHistogram::make_slice_histogram(fake_data_truth,
                                                         slice, nullptr);
    }

    // Keys are legend labels, values are SliceHistogram objects containing
    // true-space predictions from the corresponding generator models
    auto *slice_gen_map_ptr = new std::map<std::string, SliceHistogram *>();
    auto &slice_gen_map = *slice_gen_map_ptr;

    slice_gen_map["unfolded data"] = slice_unf;
    if (using_fake_data)
    {
      slice_gen_map["truth"] = slice_truth;
    }
    slice_gen_map["MicroBooNE Tune"] = slice_cv;

    for (const auto &pair : generator_truth_map)
    {
      const auto &model_name = pair.first;
      TMatrixD *truth_mat = pair.second;

      SliceHistogram *temp_slice = SliceHistogram::make_slice_histogram(
          *truth_mat, slice, nullptr);

      slice_gen_map[model_name] = temp_slice;
    }

    int var_count = 0;
    std::string diff_xsec_denom;
    std::string diff_xsec_units_denom;
    std::string diff_xsec_denom_latex;
    std::string diff_xsec_units_denom_latex;
    double other_var_width = 1.;
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
        const std::string &temp_units = var_spec.units_;
        if (!temp_units.empty())
        {
          diff_xsec_units_denom += " / " + temp_units;
          diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;
        }
      }
    }

    for (size_t av_idx : slice.active_var_indices_)
    {
      const auto &var_spec = sb.slice_vars_.at(av_idx);
      const std::string &temp_name = var_spec.name_;
      if (temp_name != "true bin number")
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

    // NOTE: This currently assumes that each slice is a 1D histogram
    // TODO: revisit as needed
    int num_slice_bins = slice_unf->hist_->GetNbinsX();
    TMatrixD trans_mat(num_slice_bins, num_slice_bins);
    for (int b = 0; b < num_slice_bins; ++b)
    {
      double width = slice_unf->hist_->GetBinWidth(b + 1);
      width *= other_var_width;
      trans_mat(b, b) = 1e38 / (width * integ_flux * num_Ar);
      std::cout << "DEBUG at " << b << " trans_mat: " << trans_mat(b, b) << " width: " << width << " integ_flux: " << integ_flux << " num_Ar: " << num_Ar << std::endl;
    }

    std::string slice_y_title;
    std::string slice_y_latex_title;
    if (var_count > 0)
    {
      slice_y_title += "d";
      slice_y_latex_title += "{$d";
      if (var_count > 1)
      {
        slice_y_title += "^{" + std::to_string(var_count) + "}";
        slice_y_latex_title += "^{" + std::to_string(var_count) + "}";
      }
      slice_y_title += "#sigma/" + diff_xsec_denom;
      slice_y_latex_title += "\\sigma / " + diff_xsec_denom_latex;
    }
    else
    {
      slice_y_title += "#sigma";
      slice_y_latex_title += "\\sigma";
    }
    slice_y_title += " (10^{-38} cm^{2}" + diff_xsec_units_denom + " / Ar)";
    slice_y_latex_title += "\\text{ }(10^{-38}\\text{ cm}^{2}" + diff_xsec_units_denom_latex + " / \\mathrm{Ar})$}";

    // Convert all slice histograms from true event counts to differential
    // cross-section units
    for (auto &pair : slice_gen_map)
    {
      auto *slice_h = pair.second;
      slice_h->transform(trans_mat);
      slice_h->hist_->GetYaxis()->SetTitle(slice_y_title.c_str());
    }

    // Also transform all of the unfolded data slice histograms which have
    // specific covariance matrices
    for (auto &sh_cov_pair : sh_cov_map)
    {
      auto &slice_h = sh_cov_pair.second;
      slice_h->transform(trans_mat);
    }

    // Keys are generator legend labels, values are the results of a chi^2
    // test compared to the unfolded data (or, in the case of the unfolded
    // data, to the fake data truth)
    std::map<std::string, SliceHistogram::Chi2Result> chi2_map;
    std::cout << '\n';
    for (const auto &pair : slice_gen_map)
    {
      const auto &name = pair.first;
      const auto *slice_h = pair.second;

      // Decide what other slice histogram should be compared to this one,
      // then calculate chi^2
      SliceHistogram *other = nullptr;
      // We don't need to compare the unfolded data to itself, so just skip to
      // the next SliceHistogram and leave a dummy Chi2Result object in the map
      if (name == "unfolded data")
      {
        chi2_map[name] = SliceHistogram::Chi2Result();
        continue;
      }
      // Compare all other distributions to the unfolded data
      else
      {
        other = slice_gen_map.at("unfolded data");
        std::cout << "DEBUG other is histogram at 'unfolded data'" << std::endl;
      }

      // Store the chi^2 results in the map
      const auto &chi2_result = chi2_map[name] = slice_h->get_chi2(*other);

      std::cout << "Slice " << sl_idx << ", " << name << ": \u03C7\u00b2 = "
                << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bin";
      if (chi2_result.num_bins_ > 1)
        std::cout << 's';
      std::cout << ", p-value = " << chi2_result.p_value_ << '\n';
    }

    std::cout << "DEBUG Unfolding Point 4" << std::endl;

    TCanvas *c1 = new TCanvas;
    slice_unf->hist_->SetLineColor(kBlack);
    slice_unf->hist_->SetLineWidth(3);
    slice_unf->hist_->SetMarkerStyle(kFullCircle);
    slice_unf->hist_->SetMarkerSize(0.8);
    slice_unf->hist_->SetStats(false);

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

      if (name == "unfolded data" || name == "truth" || name == "MicroBooNE Tune")
        continue;

      const auto &file_info = truth_file_map.at(name);
      slice_h->hist_->SetLineColor(file_info.color_);
      slice_h->hist_->SetLineStyle(file_info.style_);
      slice_h->hist_->SetLineWidth(3);

      slice_h->hist_->Draw("hist same");
    }

    slice_cv->hist_->SetStats(false);
    slice_cv->hist_->SetLineColor(kAzure - 7);
    slice_cv->hist_->SetLineWidth(3);
    slice_cv->hist_->SetLineStyle(5);
    slice_cv->hist_->Draw("hist same");

    if (using_fake_data)
    {
      slice_truth->hist_->SetStats(false);
      slice_truth->hist_->SetLineColor(kOrange);
      slice_truth->hist_->SetLineWidth(3);
      slice_truth->hist_->Draw("hist same");
    }

    slice_unf->hist_->GetYaxis()->SetRangeUser(0., ymax * 1.07);
    slice_unf->hist_->Draw("e same");
    slice_unf->hist_->SetTitle("Unfolded NuWro CC1#pi^{#pm}Xp");

    std::cout << "DEBUG Unfolding Point 6" << std::endl;

    // TLegend* lg = new TLegend( 0.45, 0.7, 0.9, 0.88 ); // Top right
    // TLegend* lg = new TLegend( 0.1, 0.7, 0.55, 0.88 ); // Top Left
    // TLegend* lg = new TLegend( 0.1, 0.2, 0.55, 0.38 ); // Bottom Left
    // TLegend* lg = new TLegend( 0.45, 0.2, 0.9, 0.38 ); // Bottom Right
    TLegend *lg = new TLegend(0.5, 0.2);
    for (const auto &pair : slice_gen_map)
    {
      const auto &name = pair.first;
      const auto *slice_h = pair.second;

      std::string label = name;
      std::ostringstream oss;
      const auto &chi2_result = chi2_map.at(name);
      oss << std::setprecision(3) << chi2_result.chi2_ << " / "
          << chi2_result.num_bins_ << " bin";
      if (chi2_result.num_bins_ > 1)
        oss << "s; p-value = " << chi2_result.p_value_;

      if (name != "unfolded data")
      {
        label += ": #chi^{2} = " + oss.str();
      }

      lg->AddEntry(slice_h->hist_.get(), label.c_str(), "l");
    }

    lg->Draw("same");

    std::string out_pdf_name = "plots/plot_unfolded_slice_";
    if (sl_idx < 10)
      out_pdf_name += "0";
    out_pdf_name += std::to_string(sl_idx) + ".pdf";
    c1->SaveAs(out_pdf_name.c_str());

    // Dump the unfolded results to text files compatible with PGFPlots
    std::map<std::string, std::vector<double>> slice_hist_table;
    std::map<std::string, std::string> slice_params_table;

    dump_slice_variables(sb, sl_idx, slice_params_table);

    for (const auto &pair : slice_gen_map)
    {
      const auto hist_name = samples_to_hist_names.at(pair.first);
      const auto *slice_hist = pair.second;
      bool include_x_coords = (hist_name == "UnfData");
      bool include_y_error = include_x_coords;
      dump_slice_histogram(hist_name, *slice_hist, slice, slice_hist_table,
                           include_y_error, include_x_coords);
    }

    std::cout << "DEBUG Unfolding Point 7" << std::endl;

    dump_slice_plot_limits(*slice_unf, *slice_cv, slice, slice_params_table);

    dump_slice_errors("UnfData", slice, sh_cov_map, slice_hist_table);

    // Dump the chi^2 test results
    for (const auto &chi2_pair : chi2_map)
    {
      const auto hist_name = samples_to_hist_names.at(chi2_pair.first);
      const auto &chi2_result = chi2_pair.second;

      // Comparing the data histogram to itself is trivial, so skip it
      if (hist_name == "UnfData")
        continue;
      else
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

    // Before moving on to the next slice, dump information about the
    // current one to new pgfplots files that can be used for offline plotting
    std::string output_file_prefix = "dump/pgfplots_slice_";
    // Use at least three digits for numbering the slice output files
    if (sl_idx < 10)
      output_file_prefix += '0';
    if (sl_idx < 100)
      output_file_prefix += '0';
    output_file_prefix += std::to_string(sl_idx);

    write_pgfplots_files(output_file_prefix, slice_hist_table,
                         slice_params_table);

    std::cout << "DEBUG Unfolding Point 8" << std::endl;
    // break; // TODO: remove !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  } // slices
  std::cout << "DEBUG Unfolding Point 9" << std::endl;
  return;

  std::cout << "DEBUG Unfolding Point 10" << std::endl;
  // ******* Also look at reco-space results
  TH1D *reco_data_hist = dynamic_cast<TH1D *>(
      syst.data_hists_.at(NFT::kOnBNB)->Clone("reco_data_hist"));
  TH1D *reco_ext_hist = syst.data_hists_.at(NFT::kExtBNB).get();
  const auto &cv_univ = syst.cv_universe();
  int num_reco_bins = reco_data_hist->GetNbinsX();

  // Clone the reco data hist twice. We will fill the clones with the CV
  // MC+EXT prediction and the constrained one
  TH1D *reco_mc_and_ext_hist = dynamic_cast<TH1D *>(
      reco_data_hist->Clone("reco_mc_and_ext_hist"));
  reco_mc_and_ext_hist->Reset();
  reco_mc_and_ext_hist->Add(reco_ext_hist);
  reco_mc_and_ext_hist->Add(cv_univ.hist_reco_.get());

  TH1D *reco_constrained_hist = dynamic_cast<TH1D *>(
      reco_data_hist->Clone("reco_constrained_hist"));
  reco_constrained_hist->Reset();

  // Get the post-constraint event counts and covariance matrix in the
  // signal region
  auto meas = syst.get_measured_events();

  for (int rb = 0; rb < num_reco_bins; ++rb)
  {

    double mcc9_err = std::sqrt(
        std::max(0., cov_mat->GetBinContent(rb + 1, rb + 1)));
    reco_mc_and_ext_hist->SetBinError(rb + 1, mcc9_err);

    if (rb >= num_ordinary_reco_bins)
    {
      double data_evts = reco_data_hist->GetBinContent(rb + 1);
      reco_constrained_hist->SetBinContent(rb + 1, data_evts);
      reco_constrained_hist->SetBinError(rb + 1, 0.);
    }
    else
    {
      double constr_pred = meas.reco_mc_plus_ext_->operator()(rb, 0);
      double constr_err = std::sqrt(
          std::max(0., meas.cov_matrix_->operator()(rb, rb)));

      reco_constrained_hist->SetBinContent(rb + 1, constr_pred);
      reco_constrained_hist->SetBinError(rb + 1, constr_err);
    }
  }

  TCanvas *c2 = new TCanvas;

  reco_data_hist->SetLineColor(kBlack);
  reco_data_hist->SetLineWidth(3);

  reco_mc_and_ext_hist->SetLineColor(kRed);
  reco_mc_and_ext_hist->SetLineStyle(2);
  reco_mc_and_ext_hist->SetLineWidth(3);

  reco_constrained_hist->SetLineColor(kBlue);
  reco_constrained_hist->SetLineStyle(9);
  reco_constrained_hist->SetLineWidth(3);

  reco_data_hist->Draw("e");
  reco_mc_and_ext_hist->Draw("same hist e");
  reco_constrained_hist->Draw("same hist e");

  reco_data_hist->Draw("same e");

  TLegend *lg2 = new TLegend(0.15, 0.7, 0.3, 0.85);
  lg2->AddEntry(reco_data_hist, using_fake_data ? "fake data" : "data",
                "l");
  lg2->AddEntry(reco_mc_and_ext_hist, "uB tune + EXT", "l");
  lg2->AddEntry(reco_constrained_hist, "post-constraint", "l");

  lg2->Draw("same");
}

int main()
{
  std::cout << "DEBUG Unfolding Point main 0" << std::endl;
  test_unfolding_simpler();
  std::cout << "DEBUG Unfolding Point main 1" << std::endl;
  return 0;
}
