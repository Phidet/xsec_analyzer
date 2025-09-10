// Standard library includes
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <set>

// ROOT includes
#include "TCanvas.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPad.h"
#include "TLine.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"
#include "utils.hh"
#include "UnfolderHelper.hh"
#include "MatrixUtils.hh"

struct inputFiles {
  std::string rootFile;
  std::string fileList;
  std::string config;
  std::string sliceConfig;
  std::string nameExtension;
};

// Function to plot a matrix with slice boundary lines
void plot_matrix_with_boundaries(TH2D* hist, const std::string& title, const std::string& xTitle, 
                                const std::string& yTitle, const std::string& outputName, 
                                const std::vector<int>& sliceBoundaries, 
                                bool setZRange = false, double zLow = 0, double zHigh = 0.1,
                                bool useScientificNotation = false, bool differentColorScheme = false) {
    TCanvas* canvas = new TCanvas("canvas", title.c_str(), 900, 800);
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.12);
    canvas->SetBottomMargin(0.12);
    canvas->SetTopMargin(0.1);
    
    // Set color scheme
    if (differentColorScheme) {
        gStyle->SetPalette(kBird);
    } else {
        util::CreateRedToBlueColorPalette(20);
    }
    
    // Configure histogram
    hist->SetTitle(title.c_str());
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetLabelSize(0.035);
    
    if (setZRange) {
        hist->GetZaxis()->SetRangeUser(zLow, zHigh);
    }
    
    if (useScientificNotation) {
        gStyle->SetPaintTextFormat("1.2e");
        hist->GetZaxis()->SetMoreLogLabels();
    }
    
    hist->SetStats(false);
    hist->Draw("colz");
    
    // Draw slice boundary lines
    std::vector<TLine*> lines;
    for (auto boundary : sliceBoundaries) {
        // Draw horizontal lines
        TLine* hLine = new TLine(0, boundary, hist->GetNbinsX(), boundary);
        hLine->SetLineStyle(2);
        hLine->SetLineColor(kBlack);
        hLine->SetLineWidth(1);
        hLine->Draw();
        lines.push_back(hLine);
        
        // Draw vertical lines
        TLine* vLine = new TLine(boundary, 0, boundary, hist->GetNbinsY());
        vLine->SetLineStyle(2);
        vLine->SetLineColor(kBlack);
        vLine->SetLineWidth(1);
        vLine->Draw();
        lines.push_back(vLine);
    }
    
    canvas->SaveAs(("plots/" + outputName + ".pdf").c_str());
    canvas->SaveAs(("plots/" + outputName + ".png").c_str());
    
    delete canvas;
    for (auto line : lines) {
        delete line;
    }
}

void plot_reco_cov_matrices() {
    // Input files configuration
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_20May25_bdts_all_quality_cuts_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc_unfold.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts.txt",
    //     "_bdts_full_uncertainty"
    // };

    // Even coarser binning and all quality cuts applied
    inputFiles input{
        "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_21May25_bdts_all_quality_cuts_17bins_full_uncertainty.root",
        "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
        "systcalc.conf",
        "ubcc1pi_slice_config_bdt_all_quality_cuts_17bins.txt",
        "_bdts_full_uncertainty"
    };

    const std::string postfix = input.nameExtension;

    // Initialize FilePropertiesManager
    auto &fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties(input.fileList);
    const std::string respmat_file_name(input.rootFile);

    // Initialize SliceBinning
    auto *sb_ptr = new SliceBinning(input.sliceConfig);
    auto &sb = *sb_ptr;

    // Initialize SystematicsCalculator
    const auto *syst_ptr = new MCC9SystematicsCalculator(respmat_file_name, input.config);
    const auto &syst = *syst_ptr;

    // Get covariance matrices
    auto *matrix_map_ptr = syst.get_covariances().release();
    auto &matrix_map = *matrix_map_ptr;
    auto *cov_mat = matrix_map.at("total").cov_matrix_.get();

    // Create output directory if it doesn't exist
    gSystem->Exec("mkdir -p plots");

    // Create correlation matrix from covariance matrix
    TMatrixD* corr = new TMatrixD(util::CovarianceMatrixToCorrelationMatrix(util::TH2DToTMatrixD(*cov_mat)));
    TH2D* corr_hist = new TH2D(util::TMatrixDToTH2D(*corr, "corr", "Correlation Matrix", 0, corr->GetNrows(), 0, corr->GetNcols()));
    
    // Determine slice boundaries
    std::vector<int> sliceBoundaries;
    
    for (size_t sl_idx = 0; sl_idx < sb.slices_.size(); ++sl_idx) {
        const auto &slice = sb.slices_.at(sl_idx);
        
        // Determine the bin range for this slice
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
        
        // Only add the end boundary of each slice (except the last one)
        if (sl_idx < sb.slices_.size() - 1) {
            sliceBoundaries.push_back(stop + 1); // +1 because histogram bins start at 1
        }
    }
    
    // Sort boundaries to ensure they're in order
    std::sort(sliceBoundaries.begin(), sliceBoundaries.end());
    
    // Plot the full correlation matrix with slice boundaries
    TH2D* cov_hist = (TH2D*)cov_mat->Clone("cov");
    cov_hist->SetTitle("Covariance Matrix");
    
    plot_matrix_with_boundaries(corr_hist, "Correlation Matrix", "Reconstructed Bins", 
                              "Reconstructed Bins", "correlation_matrix"+postfix, 
                              sliceBoundaries, true, -1, 1);
                              
    plot_matrix_with_boundaries(cov_hist, "Covariance Matrix", "Reconstructed Bins", 
                              "Reconstructed Bins", "covariance_matrix"+postfix, 
                              sliceBoundaries, false, 0, 0.1, true, true);
                              
    // Loop through each slice and create slice-specific covariance and correlation matrices
    for (size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx) {
        const auto &slice = sb.slices_.at(sl_idx);
        
        std::cout << "Processing slice " << sl_idx << std::endl;
        
        // Determine the bin range for this slice
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
        
        // Extract slice-specific matrices
        int slice_size = stop - start + 1;
        TMatrixD corr_slice(slice_size, slice_size);
        TMatrixD cov_slice(slice_size, slice_size);
        
        for (int i = start; i <= stop; i++) {
            for (int j = start; j <= stop; j++) {
                corr_slice(i - start, j - start) = corr_hist->GetBinContent(i+1, j+1);
                cov_slice(i - start, j - start) = cov_hist->GetBinContent(i+1, j+1);
            }
        }
        
        // Create histogram versions of the slices
        TH2D* trimmed_corr_hist = new TH2D(("trimmed_corr_hist slice "+std::to_string(sl_idx)).c_str(), 
                                          ("Correlation Matrix for Slice "+std::to_string(sl_idx)).c_str(),
                                          slice_size, 0, slice_size,
                                          slice_size, 0, slice_size);
                                          
        TH2D* trimmed_cov_hist = new TH2D(("trimmed_cov_hist slice "+std::to_string(sl_idx)).c_str(), 
                                         ("Covariance Matrix for Slice "+std::to_string(sl_idx)).c_str(),
                                         slice_size, 0, slice_size,
                                         slice_size, 0, slice_size);
        
        // Fill the histograms
        for (int i = 0; i < corr_slice.GetNrows(); i++) {
            for (int j = 0; j < corr_slice.GetNcols(); j++) {
                trimmed_corr_hist->SetBinContent(i+1, j+1, corr_slice(i, j));
                trimmed_cov_hist->SetBinContent(i+1, j+1, cov_slice(i, j));
            }
        }
        
        // Set up the histogram labels
        for (int i = 1; i <= trimmed_corr_hist->GetNbinsX(); i++) {
            trimmed_corr_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
            trimmed_corr_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
            
            trimmed_cov_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
            trimmed_cov_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
        }
        
        // Plot correlation matrix for this slice
        TCanvas* c_corr_slice = new TCanvas(("c_corr_slice "+std::to_string(sl_idx)).c_str(), "Correlation Matrix", 800, 600);
        util::CreateRedToBlueColorPalette(20);
        trimmed_corr_hist->SetTitleSize(0.05, "t");
        trimmed_corr_hist->GetZaxis()->SetRangeUser(-1, 1);
        trimmed_corr_hist->SetStats(false);
        trimmed_corr_hist->GetXaxis()->SetLabelSize(0.05);
        trimmed_corr_hist->GetYaxis()->SetLabelSize(0.05);
        trimmed_corr_hist->Draw("colz");
        c_corr_slice->SaveAs(("plots/correlation_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + postfix + ".pdf").c_str());
        
        // Plot covariance matrix for this slice
        TCanvas* c_cov_slice = new TCanvas(("c_cov_slice "+std::to_string(sl_idx)).c_str(), "Covariance Matrix", 800, 700);
        c_cov_slice->SetRightMargin(0.15);
        gStyle->SetPalette(kBird);
        trimmed_cov_hist->SetTitleSize(0.05, "t");
        trimmed_cov_hist->SetStats(false);
        trimmed_cov_hist->GetXaxis()->SetLabelSize(0.05);
        trimmed_cov_hist->GetYaxis()->SetLabelSize(0.05);
        trimmed_cov_hist->Draw("colz");
        
        // Add text to show the actual values
        for (int i = 0; i < slice_size; i++) {
            for (int j = 0; j < slice_size; j++) {
                TLatex* latex = new TLatex(trimmed_cov_hist->GetXaxis()->GetBinCenter(i+1), 
                                          trimmed_cov_hist->GetYaxis()->GetBinCenter(j+1), 
                                          Form("%.3g", trimmed_cov_hist->GetBinContent(i+1, j+1)));
                latex->SetTextFont(42);
                latex->SetTextSize(0.02);
                latex->SetTextAlign(22);
                latex->Draw();
            }
        }
        c_cov_slice->SaveAs(("plots/covariance_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + postfix + ".pdf").c_str());
        
        // Clean up
        delete trimmed_corr_hist;
        delete trimmed_cov_hist;
        delete c_corr_slice;
        delete c_cov_slice;
    }
    
    // Clean up
    delete corr;
    delete sb_ptr;
    delete syst_ptr;
    delete matrix_map_ptr;
    
    std::cout << "Finished creating all covariance and correlation matrix plots." << std::endl;
}

int main() {
    plot_reco_cov_matrices();
    return 0;
}
