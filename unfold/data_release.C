// Standard library includes
#include <iostream>
#include <memory>
#include <vector>
#include <stdexcept>
#include <algorithm> // added
#include <cmath>     // added

// ROOT includes
#include "TFile.h"
#include "TMatrixD.h"
#include "TH1D.h"

// UBCC1PI includes
#include "../FilePropertiesManager.hh"
#include "../FiducialVolume.hh" // added: integrated_numu_flux_in_FV, num_Ar_targets_in_FV
#include "../MCC9SystematicsCalculator.hh"
#include "../WienerSVDUnfolder.hh"
#include "../utils.hh"
#include "../SliceBinning.hh" // ADDED: to get per-bin widths from the slice config

// Use the same input format as in unfold_bnb.C
struct inputFiles
{
  std::string rootFile;
  std::string fileList;
  std::string config;
  std::string sliceConfig;   // not used here but kept for consistency
  std::string nameExtension; // not used here
};

// Zero-based index ranges [first, second) to remove rows/cols (phi and total bins)
static TMatrixD *removeMatrixBins(const TMatrixD &mat,
                                  const std::vector<std::pair<int, int>> &rowsToRemove,
                                  const std::vector<std::pair<int, int>> &colsToRemove)
{
    const int nrows = mat.GetNrows();
    const int ncols = mat.GetNcols();

    auto count_removed = [](const std::vector<std::pair<int, int>> &ranges) {
        int n = 0;
        for (const auto &r : ranges) n += (r.second - r.first);
        return n;
    };

    const int keptRows = nrows - count_removed(rowsToRemove);
    const int keptCols = ncols - count_removed(colsToRemove);

    auto in_any_range = [](int idx, const std::vector<std::pair<int, int>> &ranges) {
        for (const auto &r : ranges)
            if (idx >= r.first && idx < r.second) return true;
        return false;
    };

    auto *reduced = new TMatrixD(keptRows, keptCols);

    int rOut = 0;
    for (int r = 0; r < nrows; ++r)
    {
        if (in_any_range(r, rowsToRemove)) continue;

        int cOut = 0;
        for (int c = 0; c < ncols; ++c)
        {
            if (in_any_range(c, colsToRemove)) continue;
            (*reduced)(rOut, cOut) = mat(r, c);
            ++cOut;
        }
        ++rOut;
    }

    return reduced;
}

int data_release()
{
    // Inputs matching the unfold_bnb.C configuration
    inputFiles input{
        "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_run1234bcd5_17Oct24_testingOnly_lowPiMomThreshold_fullDetVars_fixedBackground_mergedOverflow_containedMuXSec.root",
        "../file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
        "../systcalc_unfold.conf",
        "../ubcc1pi_neutral_slice_config_mergedOverflow_noSuperscripts.txt",
        "_bnb_fixedBackground_mergedOverflow_containedMuXSec"
    };

    // Load file properties as in unfold_bnb.C
    auto &fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties(input.fileList);

    // Build systematics calculator and unfolder (same as unfold_bnb defaults)
    MCC9SystematicsCalculator syst(input.rootFile, input.config);
    std::unique_ptr<Unfolder> unfolder(
        new WienerSVDUnfolder(true, WienerSVDUnfolder::RegularizationMatrixType::kFirstDeriv));

    auto result = unfolder->unfold(syst);
    const TMatrixD &A_C = *result.add_smear_matrix_;
    const TMatrixD &unfolded_signal = *result.unfolded_signal_;
    // Also need covariance for uncertainties
    const TMatrixD &cov = *result.cov_matrix_;

    // Sanity check: expect 60x60 before cropping in this analysis setup
    if (A_C.GetNrows() != A_C.GetNcols())
        std::cerr << "Warning: A_C is not square: " << A_C.GetNrows() << "x" << A_C.GetNcols() << std::endl;

    // Remove phi bins and the total xsec bin (0-based indices; end exclusive)
    // Matches TH2D cropping in unfold_bnb.C: remove [12..26], [39..48], and total bin 60.
    // Here we translate to 0-based: [11,26), [38,48), [59,60)
    const std::vector<std::pair<int, int>> removeRanges = {
        {11, 26}, // muon phi
        {38, 48}, // pion phi
        {59, 60}  // total xsec bin
    };

    std::unique_ptr<TMatrixD> A_C_cropped(removeMatrixBins(A_C, removeRanges, removeRanges));

    // Compute events->xsec scaling like in unfold_bnb.C
    const double total_pot = syst.total_bnb_data_pot_;
    const double integ_flux = integrated_numu_flux_in_FV(total_pot);
    const double num_Ar = num_Ar_targets_in_FV();
    const double events_to_xsec = 1.0 / ((num_Ar * integ_flux) / 1e38);

    // Crop the unfolded covariance matrix to match the differential bins and scale to xsec units
    std::unique_ptr<TMatrixD> cov_unfolded_cropped(removeMatrixBins(cov, removeRanges, removeRanges));
    (*cov_unfolded_cropped) *= (events_to_xsec * events_to_xsec);

    // removed TH2D conversion for A_C and covariance; keep TMatrixD outputs

    // Helper for range checks (added)
    auto in_any_range = [](int idx, const std::vector<std::pair<int, int>> &ranges) {
        for (const auto &r : ranges)
            if (idx >= r.first && idx < r.second) return true;
        return false;
    };

    // Build a TMatrixD (nKept x 1) for the differential cross section; no uncertainties
    const int nOrig = unfolded_signal.GetNrows();
    int nKept = 0;
    for (int i = 0; i < nOrig; ++i)
        if (!in_any_range(i, removeRanges)) ++nKept;

    TMatrixD unfolded_data_diff(nKept, 1);
    {
        int out = 0;
        for (int i = 0; i < nOrig; ++i)
        {
            if (in_any_range(i, removeRanges)) continue;
            const double val = unfolded_signal(i, 0) * events_to_xsec;
            unfolded_data_diff(out, 0) = val;
            ++out;
        }
    }

    // Build a TMatrixD (nKept x 1) for the cv GENIE ("uBTune") prediction with the same kept bins (no A_C)
    TH1D *genie_cv_truth = syst.cv_universe().hist_true_.get();
    TMatrixD uBTune_diff(nKept, 1);
    {
        int out = 0;
        for (int i = 0; i < nOrig; ++i)
        {
            if (in_any_range(i, removeRanges)) continue;
            const double val = genie_cv_truth->GetBinContent(i + 1) * events_to_xsec;
            uBTune_diff(out, 0) = val;
            ++out;
        }
    }

    // Build a single-bin histogram for the total cross section with A_C divided out, and scale to xsec units
    const int n = A_C.GetNrows();
    double totalRegFactor = 1.0;
    if (n == A_C.GetNcols() && n > 0)
    {
        totalRegFactor = A_C(n - 1, n - 1);
        if (totalRegFactor == 0.0)
        {
            std::cerr << "Warning: total regularization factor is zero; skipping division." << std::endl;
            totalRegFactor = 1.0;
        }
    }
    else
    {
        std::cerr << "Warning: A_C not square or empty; assuming regularization factor = 1." << std::endl;
    }

    const double total_val = (unfolded_signal(n - 1, 0) / totalRegFactor) * events_to_xsec;
    double total_err = 0.0;
    if (cov.GetNrows() == n && cov.GetNcols() == n)
        total_err = std::sqrt(std::max(0.0, cov(n - 1, n - 1))) / totalRegFactor * events_to_xsec;
    else
        std::cerr << "Warning: covariance matrix size mismatch; setting total error to 0." << std::endl;

    TH1D *unfolded_signal_total = new TH1D("unfolded_data_total", "Total cross section;bin;10^{-38} cm^{2}/Ar", 1, 0.0, 1.0);
    unfolded_signal_total->SetBinContent(1, total_val);
    unfolded_signal_total->SetBinError(1, total_err);

    // uBTune_total as a 1x1 TMatrixD (no uncertainties)
    // const double ubtune_total_val = genie_cv_truth->GetBinContent(n) * events_to_xsec; // TH1 is 1-based
    // TMatrixD uBTune_total(1, 1);
    // uBTune_total(0, 0) = ubtune_total_val;

    // Build a TMatrixD (nKept x 1) with the physical bin widths for unfolded_data_diff/uBTune_diff
    std::unique_ptr<SliceBinning> sb_ptr(new SliceBinning(input.sliceConfig));
    const auto &sb = *sb_ptr;

    std::vector<double> widths_all;
    widths_all.reserve(nOrig);
    for (const auto &slice : sb.slices_)
    {
        const int nb = slice.hist_->GetNbinsX();
        for (int b = 1; b <= nb; ++b)
            widths_all.push_back(slice.hist_->GetBinWidth(b));
    }

    if (static_cast<int>(widths_all.size()) != nOrig)
    {
        std::cerr << "Error: concatenated bin width count (" << widths_all.size()
                  << ") != number of global truth bins (" << nOrig << ")" << std::endl;
        return 1;
    }

    TMatrixD bin_widths_diff(nKept, 1);
    {
        int out = 0;
        for (int i = 0; i < nOrig; ++i)
        {
            if (in_any_range(i, removeRanges)) continue;
            bin_widths_diff(out, 0) = widths_all[i];
            ++out;
        }
    }

    // Write to file
    static TFile *outputFile = TFile::Open("dataRelease/data_release.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie())
    {
        std::cerr << "Error: failed to create dataRelease/data_release.root" << std::endl;
        return 1;
    }

    outputFile->cd();
    // Save TMatrixD objects and the single-bin total histogram
    A_C_cropped->Write("A_C_diff");
    unfolded_data_diff.Write("unfolded_data_diff");
    cov_unfolded_cropped->Write("cov_unfolded_diff");
    uBTune_diff.Write("uBTune_diff");
    bin_widths_diff.Write("bin_widths_diff");
    // uBTune_total.Write("uBTune_total");
    unfolded_signal_total->Write("unfolded_data_total");

    outputFile->Close();

    std::cout << "Saved A_C_diff (TMatrixD), unfolded_data_diff (TMatrixD), cov_unfolded_diff (TMatrixD), uBTune_diff (TMatrixD), bin_widths_diff (TMatrixD), unfolded_data_total (TH1D) to dataRelease/data_release.root" << std::endl;
    return 0;
}