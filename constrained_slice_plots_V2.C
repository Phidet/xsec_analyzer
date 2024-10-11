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

#include "utils.hh"

using NFT = NtupleFileType;

// #define USE_FAKE_DATA "yes"

// ---------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------

void slice_plots(const bool normaliseByBinWidth)
{
    // ******************************************************************************************************************
    // Configure the input files ****************************************************************************************
    // ******************************************************************************************************************
#ifdef USE_FAKE_DATA
    // Initialize the FilePropertiesManager and tell it to treat the NuWro MC ntuples as if they were data
    auto &fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties("nuwro_file_properties_run1.txt");
    // const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_10/univmake_output_nuwro_sidband_run1_reduced_23Feb23_muonCosTheta.root");
    // const std::string respmat_file_name("/uboone/data/users/jdetje/ubcc1pi_univmake/100Percent_10/univmake_output_nuwro_sidband_run1_reduced_22Feb23_muonCosAndMom2.root");

    std::string nameExtension = "_fd_reduced_constrained";
    const bool using_fake_data = true;
#else
    auto &fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties("file_properties_testingOnly.txt");
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_reduced_muonCosTheta_run1_3Mar24.root");
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_reduced_muonCosTheta_run1234bcd5_5Mar24.root");
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_noExtraBDTCuts_reduced_muonCosTheta_run1234bcd5_5Mar24.root");
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_run1234bcd5_3Mar24.root");
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_run1234bcd5_3Mar24_gardiner.root");
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_noExtraBDTCuts_noGolden_run1234bcd5_12Mar24.root");
    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_noGolden_noProton_run1234bcd5_12Mar24.root");

    // const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_run1234bcd5_18Mar24_noExtraBDTCuts_noGolden_noProton.root");
    const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_run1234bcd5_18Mar24_noGolden_noProton.root");
    // std::string nameExtension = "_bnb_reduced_constrained";
    // std::string nameExtension = "_bnb_noExtraBDTCuts_noGolden_noProton";
    std::string nameExtension = "_bnb_noGolden_noProton";
    const bool using_fake_data = false;
#endif

    // auto* sb_ptr = new SliceBinning( "ubcc1pi_neutral_slice_config.txt" ); // all
    // auto* sb_ptr = new SliceBinning( "ubcc1pi_neutral_slice_config_reduced_muonCosTheta.txt" );
    // auto* sb_ptr = new SliceBinning( "ubcc1pi_neutral_slice_config_reduced_muonCosTheta.txt" ); // muonCosTheta and muonMomentum
    auto* sb_ptr = new SliceBinning( "ubcc1pi_neutral_slice_config_noGolden_noProton.txt" ); // no pion momentum and proton multipllcity
    auto& sb = *sb_ptr;

    // ******************************************************************************************************************
    // Get the unconstrained objects ************************************************************************************
    // ******************************************************************************************************************
    // Get the object constaining all the selection and uncertainty information
#ifdef USE_FAKE_DATA
    auto *syst_ptr = new MCC9SystematicsCalculator(
        respmat_file_name,
        "systcalc_unfold_fd_min.conf");
    auto &syst = *syst_ptr;
#else
    auto *syst_ptr = new MCC9SystematicsCalculator(
        respmat_file_name,
        "systcalc.conf");
    auto &syst = *syst_ptr;
#endif

    // Get the unconstrained reco ext & mc + ext histograms
    TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
    TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
    reco_mc_plus_ext_hist->SetDirectory( nullptr );
    reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

    TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();

    // Get the map of covariance matrices
    const auto* matrix_map_ptr = syst.get_covariances().release();
    const auto& matrix_map = *matrix_map_ptr;


    // ******************************************************************************************************************
    // Get the constrained objects **************************************************************************************
    // ******************************************************************************************************************
    // Get the object constaining all the selection and uncertainty information for the constrained calculator
#ifdef USE_FAKE_DATA
    auto* syst_ptr_constr = new ConstrainedCalculator(
    respmat_file_name,
    "systcalc_unfold_fd_min.conf" );
    auto &syst_constr = *syst_ptr_constr;
#else
    auto* syst_ptr_constr = new ConstrainedCalculator(
    respmat_file_name,
    "systcalc.conf" );
    auto &syst_constr = *syst_ptr_constr;
#endif

    // Get the measured events from the constrained calculator
    const auto meas_constr = syst_constr.get_measured_events();
    const auto &data_signal_constr = meas_constr.reco_signal_;
    const auto &data_bkgd_constr = meas_constr.reco_bkgd_;
    const auto &mc_plus_ext_constr = meas_constr.reco_mc_plus_ext_;
    const auto &mc_plus_ext_covmat_constr = meas_constr.cov_matrix_;
    const auto *data_signal_plus_bkgd_constr = new TMatrixD(*data_signal_constr, TMatrixD::EMatrixCreatorsOp2::kPlus, *data_bkgd_constr);

    // ******************************************************************************************************************
    // Get the total covariance matrices for the unconstrained **********************************************************
    // ******************************************************************************************************************
    const auto total_covmat = matrix_map.at("total").cov_matrix_.get();
    // Plot the matrix as a colz plot
    TCanvas* cCovTotal = new TCanvas;
    total_covmat->Draw("colz");
    // Remove x and y axis titles
    total_covmat->GetXaxis()->SetTitle("");
    total_covmat->GetYaxis()->SetTitle("");
    gStyle->SetOptStat(0); // Add this line to remove the stats box

    // Draw dotted lines to separate the sideband and the signal region
    // use the size of the TMatrixD total_covmat_tmatrixd and draw a vertical and a horizontal line between half the bins
    const auto nSignalBins = total_covmat->GetXaxis()->GetNbins() / 2;
    TLine *l_0_1 = new TLine(0, nSignalBins, 2*nSignalBins, nSignalBins);
    l_0_1->SetLineStyle(2);
    l_0_1->Draw();
    TLine *l_0_2 = new TLine(nSignalBins, 0, nSignalBins, 2*nSignalBins);
    l_0_2->SetLineStyle(2);
    l_0_2->Draw();

    // Add labels to the x-axis
    TLatex *ltx_0_X1 = new TLatex();
    ltx_0_X1->SetTextSize(0.03);
    ltx_0_X1->SetTextAlign(22);
    ltx_0_X1->DrawLatex(0.5*nSignalBins, -0.12*nSignalBins, "Signal reco bins");
    TLatex *ltx_0_X2 = new TLatex();
    ltx_0_X2->SetTextSize(0.03);
    ltx_0_X2->SetTextAlign(22);
    ltx_0_X2->DrawLatex(1.5*nSignalBins, -0.12*nSignalBins, "Sideband reco bins");

    // Add labels to the y-axis
    TLatex *ltx_0_Y1 = new TLatex();
    ltx_0_Y1->SetTextSize(0.03);
    ltx_0_Y1->SetTextAlign(22);
    ltx_0_Y1->SetTextAngle(90); // Rotate the text by 90 degrees
    ltx_0_Y1->DrawLatex(-0.12*nSignalBins, 0.5*nSignalBins, "Signal reco bins");
    TLatex *ltx_0_Y2 = new TLatex();
    ltx_0_Y2->SetTextSize(0.03);
    ltx_0_Y2->SetTextAlign(22);
    ltx_0_Y2->SetTextAngle(90); // Rotate the text by 90 degrees
    ltx_0_Y2->DrawLatex(-0.12*nSignalBins, 1.5*nSignalBins, "Sideband reco bins");
    

    // Add a title to the plot
    TPaveText *ptTotal = new TPaveText(0.1, 0.94, 0.9, 0.98, "brNDC");
    ptTotal->AddText("Unconstrained Total Covariance Matrix");
    ptTotal->Draw();

    std::string out_pdf_name_cov_total = "plots/cov_matrix_total_unconstrained" + nameExtension;
    cCovTotal->SaveAs((out_pdf_name_cov_total + ".pdf").c_str());
    cCovTotal->SaveAs((out_pdf_name_cov_total + ".png").c_str());
    cCovTotal->SaveAs((out_pdf_name_cov_total + ".C").c_str());
    

    // ******************************************************************************************************************
    // Get the total covariance matrices for the unconstrained case as a correlation matrix *****************************
    // ******************************************************************************************************************
    // Plot the matrix as a colz plot
    TCanvas* cCorrTotal = new TCanvas;
    const auto total_covmat_tmatrixd = util::TH2DToTMatrixD(*total_covmat);
    auto total_corrmat = util::CovarianceMatrixToCorrelationMatrix( total_covmat_tmatrixd );
    // total_corrmat.Draw("colz");
    gStyle->SetOptStat(0); // Add this line to remove the stats box

    // Convert the TMatrixT<double> to a TH2D
    int nRows = total_corrmat.GetNrows();
    int nCols = total_corrmat.GetNcols();
    TH2D h("h", "", nCols, 0, nCols, nRows, 0, nRows);
    for (int i = 1; i <= nRows; i++) {
        for (int j = 1; j <= nCols; j++) {
            h.SetBinContent(j, i, total_corrmat(i-1, j-1));
        }
    }

    // Set the range of the z-axis
    h.GetZaxis()->SetRangeUser(-1, 1);

    // Draw the histogram
    h.Draw("colz");

    // Draw dotted lines to separate the sideband and the signal region
    // use the size of the TMatrixD total_covmat_tmatrixd and draw a vertical and a horizontal line between half the bins
    TLine *l1 = new TLine(0, nSignalBins, 2*nSignalBins, nSignalBins);
    l1->SetLineStyle(2);
    l1->Draw();
    TLine *l2 = new TLine(nSignalBins, 0, nSignalBins, 2*nSignalBins);
    l2->SetLineStyle(2);
    l2->Draw();

    // Add labels to the x-axis
    TLatex *ltxX1 = new TLatex();
    ltxX1->SetTextSize(0.03);
    ltxX1->SetTextAlign(22);
    ltxX1->DrawLatex(0.5*nSignalBins, -0.12*nSignalBins, "Signal bins");
    TLatex *ltxX2 = new TLatex();
    ltxX2->SetTextSize(0.03);
    ltxX2->SetTextAlign(22);
    ltxX2->DrawLatex(1.5*nSignalBins, -0.12*nSignalBins, "Sideband bins");

    // Add labels to the y-axis
    TLatex *ltxY1 = new TLatex();
    ltxY1->SetTextSize(0.03);
    ltxY1->SetTextAlign(22);
    ltxY1->SetTextAngle(90); // Rotate the text by 90 degrees
    ltxY1->DrawLatex(-0.12*nSignalBins, 0.5*nSignalBins, "Signal bins");
    TLatex *ltxY2 = new TLatex();
    ltxY2->SetTextSize(0.03);
    ltxY2->SetTextAlign(22);
    ltxY2->SetTextAngle(90); // Rotate the text by 90 degrees
    ltxY2->DrawLatex(-0.12*nSignalBins, 1.5*nSignalBins, "Sideband bins");
    

    // Add a title to the plot
    TPaveText *ptTotalCorr = new TPaveText(0.1, 0.94, 0.9, 0.98, "brNDC");
    ptTotalCorr->AddText("Unconstrained Total Correlation Matrix");
    ptTotalCorr->Draw();

    std::string out_pdf_name_corr_total = "plots/corr_matrix_total_unconstrained" + nameExtension;
    cCorrTotal->SaveAs((out_pdf_name_corr_total + ".pdf").c_str());
    cCorrTotal->SaveAs((out_pdf_name_corr_total + ".png").c_str());
    cCorrTotal->SaveAs((out_pdf_name_corr_total + ".C").c_str());


    // ******************************************************************************************************************
    // Plot the histograms for the entire reco space ********************************************************************
    // ******************************************************************************************************************
    // Get the total covariance matrices for the unconstrained and constrained cases
    auto* cov_mat = matrix_map.at( "total" ).cov_matrix_.get();

    // Get the number of ordinary and sideband reco bins
    int num_ordinary_reco_bins = 0;
    int num_sideband_reco_bins = 0;
    for ( int b = 0; b < syst_constr.reco_bins_.size(); ++b ) 
    {
        const auto& rbin = syst_constr.reco_bins_.at( b );
        if ( rbin.type_ == kSidebandRecoBin ) ++num_sideband_reco_bins;
        else ++num_ordinary_reco_bins;
    }

    // Get the number of true signal bins
    int num_true_signal_bins = 0;
    for ( int t = 0; t < syst_constr.true_bins_.size(); ++t )
    {
        const auto& tbin = syst_constr.true_bins_.at( t );
        if ( tbin.type_ == kSignalTrueBin ) ++num_true_signal_bins;
    }


    int num_reco_bins = reco_bnb_hist->GetNbinsX();

    TH1D* reco_mc_and_ext_hist_constr = dynamic_cast< TH1D* >(
    reco_bnb_hist->Clone( "reco_mc_and_ext_hist_constr" )
    );
    reco_mc_and_ext_hist_constr->Reset();

    TH1D* reco_mc_and_ext_hist_no_constr = dynamic_cast< TH1D* >(
        reco_mc_plus_ext_hist->Clone( "reco_mc_and_ext_hist_no_constr" )
    );

    // Get the post-constraint event counts and covariance matrix in the signal region
    for ( int rb = 0; rb < num_reco_bins; ++rb ) {
        if(cov_mat->GetBinContent(rb + 1, rb + 1) < 0)
            throw std::runtime_error("Negative diagonal covariance matrix element");

        const double err = std::sqrt(
            std::max( 0., cov_mat->GetBinContent(rb + 1, rb + 1) )
        );
        reco_mc_and_ext_hist_no_constr->SetBinError( rb + 1, err );

        if ( rb >= num_ordinary_reco_bins )
        {
            double data_evts = reco_bnb_hist->GetBinContent( rb + 1 );
            reco_mc_and_ext_hist_constr->SetBinContent( rb + 1, data_evts );
            reco_mc_and_ext_hist_constr->SetBinError( rb + 1, 0. );
            // reco_mc_and_ext_hist_no_constr->SetBinContent( rb + 1, data_evts );
            // reco_mc_and_ext_hist_no_constr->SetBinError( rb + 1, 0 );
        }
        else 
        {
            double constr_pred = mc_plus_ext_constr->operator()( rb, 0 );
            double constr_err = std::sqrt(
            std::max( 0., mc_plus_ext_covmat_constr->operator()(rb, rb) )
            );

            reco_mc_and_ext_hist_constr->SetBinContent( rb + 1, constr_pred );
            reco_mc_and_ext_hist_constr->SetBinError( rb + 1, constr_err );
        }
    }




    // ******************************************************************************************************************
    // Plot the histograms for sideband region **************************************************************************
    // ******************************************************************************************************************
    TCanvas* c2 = new TCanvas;

    reco_bnb_hist->SetLineColor( kBlack );
    reco_bnb_hist->SetLineWidth( 1 );

    // reco_mc_and_ext_hist->SetLineColor( kRed );
    // reco_mc_and_ext_hist->SetLineStyle( 2 );
    // reco_mc_and_ext_hist->SetLineWidth( 3 );

    reco_mc_and_ext_hist_no_constr->SetLineColor( kOrange );
    // reco_mc_and_ext_hist_no_constr->SetLineStyle( 2 );
    reco_mc_and_ext_hist_no_constr->SetLineWidth( 4 );

    reco_mc_and_ext_hist_constr->SetLineColor( kGreen );
    // reco_mc_and_ext_hist_constr->SetLineStyle( 9 );
    reco_mc_and_ext_hist_constr->SetLineWidth( 2 );


    // Create copies of the histograms
    TH1* reco_mc_and_ext_hist_no_constr_copy = (TH1*)reco_mc_and_ext_hist_no_constr->Clone();
    TH1* reco_mc_and_ext_hist_constr_copy = (TH1*)reco_mc_and_ext_hist_constr->Clone();

    // Set the x-axis range to the first num_ordinary_reco_bins bins
    reco_bnb_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    // reco_mc_and_ext_hist_constr->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    // reco_mc_and_ext_hist_no_constr->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    reco_mc_and_ext_hist_no_constr_copy->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    reco_mc_and_ext_hist_constr_copy->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);

    // Apply fill styles to the copies
    reco_mc_and_ext_hist_no_constr_copy->SetFillStyle(3004);
    reco_mc_and_ext_hist_no_constr_copy->SetFillColor(kOrange);
    reco_mc_and_ext_hist_constr_copy->SetFillStyle(3005);
    reco_mc_and_ext_hist_constr_copy->SetFillColor(kGreen);

    // Draw the histograms
    reco_mc_and_ext_hist_no_constr_copy->Draw("e2");
    reco_mc_and_ext_hist_no_constr->Draw("same hist");
    reco_mc_and_ext_hist_constr->Draw("same hist");
    reco_bnb_hist->Draw("same E1");

    TLegend* lg2 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    lg2->AddEntry( reco_bnb_hist, using_fake_data ? "fake data" : "data", "l" );
    lg2->AddEntry( reco_mc_and_ext_hist_no_constr, "uB tune + EXT", "l" );
    lg2->AddEntry( reco_mc_and_ext_hist_constr, "uB tune + EXT (constrained)", "l" );
    lg2->Draw( "same" );
    c2->SaveAs(("plots/constrained_plots_reco_sideband" + nameExtension + ".pdf").c_str());
    c2->SaveAs(("plots/constrained_plots_reco_sideband" + nameExtension + ".png").c_str());
    c2->SaveAs(("plots/constrained_plots_reco_sideband" + nameExtension + ".C").c_str());
    std::cout << "########################## Saved " << "plots/constrained_plots_reco_sideband" + nameExtension + ".pdf" << " ##########################" << std::endl;


    // ******************************************************************************************************************
    // Plot the histograms for signal region ****************************************************************************
    // ******************************************************************************************************************
    TCanvas* c3 = new TCanvas;

    reco_mc_and_ext_hist_constr->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_mc_and_ext_hist_constr_copy->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_bnb_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_mc_and_ext_hist_no_constr->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_mc_and_ext_hist_no_constr_copy->GetXaxis()->SetRange(1, num_ordinary_reco_bins);

    // Draw the histograms
    reco_mc_and_ext_hist_no_constr_copy->Draw("e2");
    reco_mc_and_ext_hist_no_constr->Draw("same hist");
    reco_mc_and_ext_hist_constr_copy->Draw("same e2");
    reco_mc_and_ext_hist_constr->Draw("same hist");
    reco_bnb_hist->Draw("same E1");

    TLegend* lg3 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    lg3->AddEntry( reco_bnb_hist, using_fake_data ? "fake data" : "data", "l" );
    lg3->AddEntry( reco_mc_and_ext_hist_no_constr, "uB tune + EXT", "l" );
    // lg3->AddEntry( reco_mc_and_ext_hist_no_constr, "uB tune + EXT (safety check)", "l" );
    lg3->AddEntry( reco_mc_and_ext_hist_constr, "uB tune + EXT (constrained)", "l" );
    lg3->Draw( "same" );


    c3->SaveAs(("plots/constrained_plots_reco" + nameExtension + ".pdf").c_str());
    c3->SaveAs(("plots/constrained_plots_reco" + nameExtension + ".png").c_str());
    c3->SaveAs(("plots/constrained_plots_reco" + nameExtension + ".C").c_str());
    std::cout << "########################## Saved " << "plots/constrained_plots_reco" + nameExtension + ".pdf" << " ##########################" << std::endl;

    // Deleting dynamically allocated objects to avoid memory leaks
    delete data_signal_plus_bkgd_constr;
    delete syst_ptr;
    delete syst_ptr_constr;
    delete sb_ptr;
}

int constrained_slice_plots_V2()
{
    slice_plots(false);
    std::cout << "------------All Done------------" << std::endl;
    return 0;
}
