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
#include "ConstrainedCalculator.hh"
#include "PlotUtils.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"

using NFT = NtupleFileType;

#define USE_FAKE_DATA "yes"

void scale_by_bin_width(SliceHistogram *pSlice)
{
    int num_slice_bins = pSlice->hist_->GetNbinsX();
    TMatrixD trans_mat(num_slice_bins, num_slice_bins);
    for (int b = 0; b < num_slice_bins; ++b)
    {
        const auto width = pSlice->hist_->GetBinWidth(b + 1);
        // width *= other_var_width;
        trans_mat(b, b) = 1 / width;
    }
    pSlice->transform(trans_mat);
}

void slice_plots(const bool normaliseByBinWidth)
{
    std::cout << "DEBUG constrained_slice_plots Point 0" << std::endl;
#ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro
    // MC ntuples as if they were data
    auto &fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties("nuwro_file_properties.txt");
    const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_9/univmake_output_nuwro_with_sidband_5Jan23.root");
    // const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_9/univmake_output_nuwro_sidband_reduced_28Dec23.root");
    // const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_9/univmake_output_nuwro_sidband_reduced_04Jan23.root");

    auto *syst_ptr = new MCC9SystematicsCalculator(
        respmat_file_name,
        "systcalc_unfold_fd_min.conf");
    // std::string nameExtension = "_fd_reduced";

    auto* syst_ptr_constr = new ConstrainedCalculator(
      respmat_file_name,
      "systcalc_unfold_fd_min.conf" );
    // std::string nameExtension = "_fd_reduced_constrained";
    const bool using_fake_data = true;
#else
    auto *syst_ptr = new MCC9SystematicsCalculator(
        "...",
        "systcalc.conf");
    std::string nameExtension = "_bnb";
    const bool using_fake_data = false;
#endif

    auto &syst = *syst_ptr;
    auto &syst_constr = *syst_ptr_constr;

    // Get the measured events from the systematics calculator
    const auto meas = syst.get_measured_events();
    const auto &data_signal = meas.reco_signal_;
    const auto &data_bkgd = meas.reco_bkgd_;
    const auto &data_mc_plus_ext = meas.reco_mc_plus_ext_;
    const auto &data_covmat = meas.cov_matrix_;
    const auto data_signal_plus_bkgd = std::make_unique<TMatrixD>(*data_signal + *data_bkgd);

    // Get the measured events from the constrained calculator
    const auto meas_constr = syst_constr.get_measured_events();
    const auto &data_signal_constr = meas_constr.reco_signal_;
    const auto &data_bkgd_constr = meas_constr.reco_bkgd_;
    const auto &data_mc_plus_ext_constr = meas_constr.reco_mc_plus_ext_;
    const auto &data_covmat_constr = meas_constr.cov_matrix_;
    const auto data_signal_plus_bkgd_constr = std::make_unique<TMatrixD>(*data_signal_constr + *data_bkgd_constr);
    
    TMatrixD data_covmat_placeholder = *data_covmat; // Copy the matrix
    data_covmat_placeholder = 999; // Set all values to 999

    auto* sb_ptr = new SliceBinning( "ubcc1pi_neutral_slice_config.txt" );
    auto& sb = *sb_ptr;

    for (size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx)
    {
        // continue;
        // if(sl_idx != 0) continue; // todo remove this line !!!!!!!!!!!!
        const auto &slice = sb.slices_.at(sl_idx);

        // Make histograms in slice space.
        SliceHistogram *slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
            *data_mc_plus_ext, slice, &data_covmat_placeholder);

        SliceHistogram *slice_reco_signal_plus_bkgd = SliceHistogram::make_slice_histogram(
            *data_signal_plus_bkgd, slice, data_covmat.get());

        // Make constrained versions of the histograms
        SliceHistogram *slice_mc_plus_ext_constr= SliceHistogram::make_slice_histogram(
            *data_mc_plus_ext_constr, slice, &data_covmat_placeholder);

        SliceHistogram *slice_reco_signal_plus_bkgd_constr = SliceHistogram::make_slice_histogram(
            *data_signal_plus_bkgd_constr, slice, data_covmat_constr.get());

        if(normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext);
        if(normaliseByBinWidth) scale_by_bin_width(slice_reco_signal_plus_bkgd);
        if(normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext_constr);
        if(normaliseByBinWidth) scale_by_bin_width(slice_reco_signal_plus_bkgd_constr);

        TCanvas* c = new TCanvas;
        // Set the color and line thickness of the histograms
        slice_mc_plus_ext->hist_->SetLineColor(kBlue);
        slice_mc_plus_ext->hist_->SetLineWidth(3);
        slice_reco_signal_plus_bkgd->hist_->SetLineColor(kBlack);
        slice_reco_signal_plus_bkgd->hist_->SetLineWidth(3);

        slice_mc_plus_ext_constr->hist_->SetLineColor(kGreen);
        slice_mc_plus_ext_constr->hist_->SetLineStyle(2);
        slice_reco_signal_plus_bkgd_constr->hist_->SetLineColor(kRed);
        slice_reco_signal_plus_bkgd_constr->hist_->SetLineStyle(2);

        // Draw the histograms
        slice_mc_plus_ext->hist_->Draw("hist");
        slice_reco_signal_plus_bkgd->hist_->Draw("hist same");
        slice_mc_plus_ext_constr->hist_->Draw("hist same");
        slice_reco_signal_plus_bkgd_constr->hist_->Draw("hist same");

        std::string out_pdf_name = "plots/constrained_plots_slice_";
        if ( sl_idx < 10 ) out_pdf_name += "0";
        out_pdf_name += std::to_string( sl_idx )+".pdf";
        c->SaveAs(out_pdf_name.c_str());
        delete c;
    }






    // ---------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------





    auto* matrix_map_ptr = syst_constr.get_covariances().release();
    auto& matrix_map = *matrix_map_ptr;

    auto* mcc9_map_ptr = syst.get_covariances().release();
    auto& mcc9_map = *mcc9_map_ptr;

    auto* cov_mat = matrix_map.at( "total" ).cov_matrix_.get();
    auto* mcc9_cov_mat = mcc9_map.at( "total" ).cov_matrix_.get();


    std::cout << "CC size = " << syst_constr.get_covariance_matrix_size()
    << ", MCC9 size = " << syst.get_covariance_matrix_size() << '\n';

    int num_ordinary_reco_bins = 0;
    int num_sideband_reco_bins = 0;
    for ( int b = 0; b < syst_constr.reco_bins_.size(); ++b ) {
    const auto& rbin = syst_constr.reco_bins_.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) ++num_sideband_reco_bins;
    else ++num_ordinary_reco_bins;
    }

    int num_true_signal_bins = 0;
    for ( int t = 0; t < syst_constr.true_bins_.size(); ++t ) {
    const auto& tbin = syst_constr.true_bins_.at( t );
    if ( tbin.type_ == kSignalTrueBin ) ++num_true_signal_bins;
    }

    std::cout << "NUM ORDINARY RECO BINS = " << num_ordinary_reco_bins << '\n';
    std::cout << "NUM TRUE SIGNAL BINS = " << num_true_signal_bins << '\n';

    TH1D* reco_data_hist = dynamic_cast< TH1D* >(
    syst_constr.data_hists_.at( NFT::kOnBNB )->Clone( "reco_data_hist" )
    );
    TH1D* reco_ext_hist = syst_constr.data_hists_.at( NFT::kExtBNB ).get();
    const auto& cv_univ = syst_constr.cv_universe();
    int num_reco_bins = reco_data_hist->GetNbinsX();

    // Clone the reco data hist twice. We will fill the clones with the CV
    // MC+EXT prediction and the constrained one
    TH1D* reco_mc_and_ext_hist = dynamic_cast< TH1D* >(
    reco_data_hist->Clone( "reco_mc_and_ext_hist" )
    );
    reco_mc_and_ext_hist->Reset();
    reco_mc_and_ext_hist->Add( reco_ext_hist );
    reco_mc_and_ext_hist->Add( cv_univ.hist_reco_.get() );

    TH1D* reco_constrained_hist = dynamic_cast< TH1D* >(
    reco_data_hist->Clone( "reco_constrained_hist" )
    );
    reco_constrained_hist->Reset();

    // Get the post-constraint event counts and covariance matrix in the
    // signal region
    // auto meas = syst_constr.get_measured_events();

    for ( int rb = 0; rb < num_reco_bins; ++rb ) {

    double mcc9_err = std::sqrt(
        std::max( 0., mcc9_cov_mat->GetBinContent(rb + 1, rb + 1) )
    );
    reco_mc_and_ext_hist->SetBinError( rb + 1, mcc9_err );

    if ( rb >= num_ordinary_reco_bins ) {
        double data_evts = reco_data_hist->GetBinContent( rb + 1 );
        reco_constrained_hist->SetBinContent( rb + 1, data_evts );
        reco_constrained_hist->SetBinError( rb + 1, 0. );
    }
    else {
        double constr_pred = meas.reco_mc_plus_ext_->operator()( rb, 0 );
        double constr_err = std::sqrt(
        std::max( 0., meas.cov_matrix_->operator()(rb, rb) )
        );

        reco_constrained_hist->SetBinContent( rb + 1, constr_pred );
        reco_constrained_hist->SetBinError( rb + 1, constr_err );
    }

    }

    TCanvas* c2 = new TCanvas;

    reco_data_hist->SetLineColor( kBlack );
    reco_data_hist->SetLineWidth( 1 );

    reco_mc_and_ext_hist->SetLineColor( kRed );
    reco_mc_and_ext_hist->SetLineStyle( 2 );
    reco_mc_and_ext_hist->SetLineWidth( 3 );

    reco_constrained_hist->SetLineColor( kBlue );
    reco_constrained_hist->SetLineStyle( 9 );
    reco_constrained_hist->SetLineWidth( 2 );

    // Set the x-axis range to the first num_ordinary_reco_bins bins
    reco_data_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    reco_mc_and_ext_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    reco_constrained_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);

    // Set the y-axis range to start at 0
    // const auto yMaxSideband = std::max({reco_data_hist->GetMaximum(), reco_mc_and_ext_hist->GetMaximum(), reco_constrained_hist->GetMaximum()});
    // reco_data_hist->GetYaxis()->SetRangeUser(0, 1.1*yMaxSideband);
    // reco_mc_and_ext_hist->GetYaxis()->SetRangeUser(0, 1.1*yMaxSideband);
    // reco_constrained_hist->GetYaxis()->SetRangeUser(0, 1.1*yMaxSideband);

    // Set the y-axis to a logarithmic scale
    c2->SetLogy();

    // Draw the histograms
    reco_data_hist->Draw( "e" );
    reco_mc_and_ext_hist->Draw( "same hist e" );
    reco_constrained_hist->Draw( "same hist e" );
    reco_data_hist->Draw( "same e" );

    TLegend* lg2 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    lg2->AddEntry( reco_data_hist, using_fake_data ? "fake data" : "data", "l" );
    lg2->AddEntry( reco_mc_and_ext_hist, "uB tune + EXT", "l" );
    lg2->AddEntry( reco_constrained_hist, "post-constraint", "l" );
    lg2->Draw( "same" );
    c2->SaveAs("plots/constrained_plots_reco_sideband.pdf");


    TCanvas* c3 = new TCanvas;

    reco_data_hist->SetLineColor( kBlack );
    reco_data_hist->SetLineWidth( 1 );

    reco_mc_and_ext_hist->SetLineColor( kRed );
    reco_mc_and_ext_hist->SetLineStyle( 2 );
    reco_mc_and_ext_hist->SetLineWidth( 3 );

    reco_constrained_hist->SetLineColor( kBlue );
    reco_constrained_hist->SetLineStyle( 9 );
    reco_constrained_hist->SetLineWidth( 2 );

    // Set the x-axis range to the first num_ordinary_reco_bins bins
    reco_data_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_mc_and_ext_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_constrained_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);

    // Set the y-axis range to start at 0
    // const auto yMax = std::max({reco_data_hist->GetMaximum(), reco_mc_and_ext_hist->GetMaximum(), reco_constrained_hist->GetMaximum()});
    // reco_data_hist->GetYaxis()->SetRangeUser(0, 1.1*yMax);
    // reco_mc_and_ext_hist->GetYaxis()->SetRangeUser(0, 1.1*yMax);
    // reco_constrained_hist->GetYaxis()->SetRangeUser(0, 1.1*yMax);

    // Set the y-axis to a logarithmic scale
    c3->SetLogy();

    // Draw the histograms
    reco_data_hist->Draw( "e" );
    reco_mc_and_ext_hist->Draw( "same hist e" );
    reco_constrained_hist->Draw( "same hist e" );
    reco_data_hist->Draw( "same e" );

    TLegend* lg3 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    lg3->AddEntry( reco_data_hist, using_fake_data ? "fake data" : "data", "l" );
    lg3->AddEntry( reco_mc_and_ext_hist, "uB tune + EXT", "l" );
    lg3->AddEntry( reco_constrained_hist, "post-constraint", "l" );
    lg3->Draw( "same" );
    c3->SaveAs("plots/constrained_plots_reco.pdf");

    // for (size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx)
    // {
    //     const auto &slice = sb.slices_.at(sl_idx);

    //     // Check if bin_map_ is not empty
    //     // Initialize the lowest and highest values to the first value in the map
    //     size_t lowestValue = *slice.bin_map_.begin()->second.begin();
    //     size_t highestValue = lowestValue;

    //     // Iterate over the map
    //     for (const auto &pair : slice.bin_map_) {
    //         for (const auto &value : pair.second) {
    //             // Update the lowest and highest values
    //             lowestValue = std::min(lowestValue, value);
    //             highestValue = std::max(highestValue, value);
    //         }
    //     }

    //     TCanvas* c3 = new TCanvas;

    //     reco_data_hist->SetLineColor( kBlack );
    //     reco_data_hist->SetLineWidth( 2 );

    //     reco_mc_and_ext_hist->SetLineColor( kRed );
    //     reco_mc_and_ext_hist->SetLineStyle( 2 );
    //     reco_mc_and_ext_hist->SetLineWidth( 4 );

    //     reco_constrained_hist->SetLineColor( kBlue );
    //     reco_constrained_hist->SetLineStyle( 9 );
    //     reco_constrained_hist->SetLineWidth( 3 );

    //     // Set the x-axis range to the first num_ordinary_reco_bins bins
    //     reco_data_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1+lowestValue, 2*num_ordinary_reco_bins+highestValue);
    //     reco_mc_and_ext_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1+lowestValue, 2*num_ordinary_reco_bins+highestValue);
    //     reco_constrained_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1+lowestValue, 2*num_ordinary_reco_bins+highestValue);

    //     // Set the y-axis range to start at 0
    //     // const auto yMaxSideband = std::max({reco_data_hist->GetMaximum(), reco_mc_and_ext_hist->GetMaximum(), reco_constrained_hist->GetMaximum()});
    //     // reco_data_hist->GetYaxis()->SetRangeUser(0, 1.1*yMaxSideband);
    //     // reco_mc_and_ext_hist->GetYaxis()->SetRangeUser(0, 1.1*yMaxSideband);
    //     // reco_constrained_hist->GetYaxis()->SetRangeUser(0, 1.1*yMaxSideband);

    //     // Draw the histograms
    //     reco_data_hist->Draw( "e" );
    //     reco_mc_and_ext_hist->Draw( "same hist e" );
    //     reco_constrained_hist->Draw( "same hist e" );
    //     reco_data_hist->Draw( "same e" );

    //     TLegend* lg3 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    //     lg3->AddEntry( reco_data_hist, using_fake_data ? "fake data" : "data", "l" );
    //     lg3->AddEntry( reco_mc_and_ext_hist, "uB tune + EXT", "l" );
    //     lg3->AddEntry( reco_constrained_hist, "post-constraint", "l" );
    //     lg3->Draw( "same" );
    //     c3->SaveAs(("plots/constrained_plots_reco_sideband_slice_0" + std::to_string(sl_idx) + ".pdf").c_str());


    //     TCanvas* c4 = new TCanvas;

    //     reco_data_hist->SetLineColor( kBlack );
    //     reco_data_hist->SetLineWidth( 2 );

    //     reco_mc_and_ext_hist->SetLineColor( kRed );
    //     reco_mc_and_ext_hist->SetLineStyle( 2 );
    //     reco_mc_and_ext_hist->SetLineWidth( 4 );

    //     reco_constrained_hist->SetLineColor( kBlue );
    //     reco_constrained_hist->SetLineStyle( 9 );
    //     reco_constrained_hist->SetLineWidth( 3 );

    //     // Set the x-axis range to the first num_ordinary_reco_bins bins
    //     reco_data_hist->GetXaxis()->SetRange(1+lowestValue, num_ordinary_reco_bins+highestValue);
    //     reco_mc_and_ext_hist->GetXaxis()->SetRange(1+lowestValue, num_ordinary_reco_bins+highestValue);
    //     reco_constrained_hist->GetXaxis()->SetRange(1+lowestValue, num_ordinary_reco_bins+highestValue);

    //     // Set the y-axis range to start at 0
    //     // const auto yMax = std::max({reco_data_hist->GetMaximum(), reco_mc_and_ext_hist->GetMaximum(), reco_constrained_hist->GetMaximum()});
    //     // reco_data_hist->GetYaxis()->SetRangeUser(0, 1.1*yMax);
    //     // reco_mc_and_ext_hist->GetYaxis()->SetRangeUser(0, 1.1*yMax);
    //     // reco_constrained_hist->GetYaxis()->SetRangeUser(0, 1.1*yMax);

    //     // Draw the histograms
    //     reco_data_hist->Draw( "e" );
    //     reco_mc_and_ext_hist->Draw( "same hist e" );
    //     reco_constrained_hist->Draw( "same hist e" );
    //     reco_data_hist->Draw( "same e" );

    //     TLegend* lg4 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    //     lg4->AddEntry( reco_data_hist, using_fake_data ? "fake data" : "data", "l" );
    //     lg4->AddEntry( reco_mc_and_ext_hist, "uB tune + EXT", "l" );
    //     lg4->AddEntry( reco_constrained_hist, "post-constraint", "l" );
    //     lg4->Draw( "same" );
    //     c4->SaveAs(("plots/constrained_plots_reco_slice_0" + std::to_string(sl_idx) + ".pdf").c_str());
    // }
}

int constrained_slice_plots()
{
    slice_plots(false);
    std::cout << "------------All Done------------" << std::endl;
    return 0;
}
