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


void plot_frac_slice( const Slice& slice, const std::map< std::string, CovMatrix >& matrix_map, TH1D* reco_mc_plus_ext_hist, const int sl_idx, const std::string nameExtension ){
    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    std::vector< std::string > cov_mat_keys = { "total", "detVar_total", "flux", "reint", "xsec_multi", "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH", "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad", "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie", "POT", "numTargets", "MCstats", "EXTstats", "BNBstats"};

    #ifdef USE_FAKE_DATA
    cov_mat_keys = { "total", "MCstats", "EXTstats", "BNBstats", "xsec_multi", "xsec_unisim", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie"};
    #endif

    // Loop over the various systematic uncertainties
    int color = 0;
    for ( const auto& pair : matrix_map ) {
      const auto& key = pair.first;
      const auto& cov_matrix = pair.second;
      SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, &cov_matrix );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
        double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
        double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        double frac = 0.;
        if ( y > 0. ) frac = err / y;
        slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
        slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
      }

      // Check whether the current covariance matrix name is present in
      // the vector defined above this loop. If it isn't, don't bother to
      // plot it, and just move on to the next one.
      auto cbegin = cov_mat_keys.cbegin();
      auto cend = cov_mat_keys.cend();
      auto iter = std::find( cbegin, cend, key );
      if ( iter == cend )
        continue;

      frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

      if ( color <= 9 ) ++color;
      if ( color == 5 ) ++color;
      if ( color >= 10 ) color += 10;

      slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );
    }

    TCanvas* c2 = new TCanvas;
    TLegend* lg2 = new TLegend( 0.2, 0.3);

    auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
      total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineStyle( 9 );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->Draw( "hist" );
    total_frac_err_hist->SetTitle("Fractional Uncertainty of Selected #nu_{#mu}CC1#pi^{#pm}Np, N #geq 0 Events");
    total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");

    lg2->AddEntry( total_frac_err_hist, "total", "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == "total" ) continue;
      if (name.size() >= 5 && name.substr(name.size() - 5) == "stats") 
      {
        hist->SetLineStyle( 2 );
      }

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );
    }
    lg2->Draw( "same" );

    std::string frac_out_pdf_name = "plots/plot_frac_slice_";
    if ( sl_idx < 10 ) frac_out_pdf_name += "0";
    frac_out_pdf_name += std::to_string( sl_idx ) + nameExtension;
    c2->SaveAs( (frac_out_pdf_name + ".pdf").c_str() );
    c2->SaveAs( (frac_out_pdf_name + ".png").c_str() );
    c2->SaveAs( (frac_out_pdf_name + ".C").c_str() );
    std::cout << "########################## Saved " << frac_out_pdf_name << " ##########################" << std::endl;
}

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
    fpm.load_file_properties("file_properties_run1.txt");   
    const std::string respmat_file_name("/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_sideband_reduced_muonCosTheta_run1_3Mar24.root");
    std::string nameExtension = "_bnb_reduced_constrained";
    const bool using_fake_data = false;
#endif

    // auto* sb_ptr = new SliceBinning( "ubcc1pi_neutral_slice_config_reduced_muonCosTheta.txt" );
    auto* sb_ptr = new SliceBinning( "ubcc1pi_neutral_slice_config_reduced_muonCosTheta.txt" ); // muonCosTheta and muonMomentum
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

    // Get the measured events from the systematics calculator
    const auto meas = syst.get_measured_events();
    const auto &data_signal = meas.reco_signal_; // Background-subtracted data event counts in the ordinary reco bins
    const auto &data_bkgd = meas.reco_bkgd_; // Background that was subtracted from each reco bin to form the signal measurement
    const auto &data_mc_plus_ext = meas.reco_mc_plus_ext_; // Total MC+EXT prediction in each reco bin
    const auto &data_covmat = meas.cov_matrix_; // Covariance matrix for the background-subtracted data
    const auto *data_signal_plus_bkgd = new TMatrixD(*data_signal, TMatrixD::EMatrixCreatorsOp2::kPlus, *data_bkgd);

    // Get the unconstrained reco ext & mc + ext histograms
    TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
    TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
    reco_mc_plus_ext_hist->SetDirectory( nullptr );
    reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

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

    // Get the constrained reco ext & mc + ext histograms
    TH1D* reco_ext_hist_constr = syst_constr.data_hists_.at( NFT::kExtBNB ).get();
    TH1D* reco_mc_plus_ext_hist_constr = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist_constr") );
    reco_mc_plus_ext_hist_constr->SetDirectory( nullptr );
    reco_mc_plus_ext_hist_constr->Add( syst_constr.cv_universe().hist_reco_.get() );

    // Get the measured events from the constrained calculator
    const auto meas_constr = syst_constr.get_measured_events();
    const auto &data_signal_constr = meas_constr.reco_signal_;
    const auto &data_bkgd_constr = meas_constr.reco_bkgd_;
    const auto &data_mc_plus_ext_constr = meas_constr.reco_mc_plus_ext_;
    const auto &data_covmat_constr = meas_constr.cov_matrix_;
    const auto *data_signal_plus_bkgd_constr = new TMatrixD(*data_signal_constr, TMatrixD::EMatrixCreatorsOp2::kPlus, *data_bkgd_constr);

    // Get the map of covariance matrices (identical to unconstrained map)
    const auto* matrix_map_constr_ptr = syst_constr.get_covariances().release();
    const auto& matrix_map_constr = *matrix_map_constr_ptr;


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
    total_corrmat.Draw("colz");
    gStyle->SetOptStat(0); // Add this line to remove the stats box

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
    // Loop over slices *************************************************************************************************
    // ******************************************************************************************************************
    for (size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx)
    {
        std::cout << "DEBUG constrained_slice_plots Point 3" << std::endl;
        const auto &slice = sb.slices_.at(sl_idx);

        // Make histograms in slice space.
        SliceHistogram *slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
            *data_mc_plus_ext, slice, data_covmat.get());

        SliceHistogram *slice_reco_signal_plus_bkgd = SliceHistogram::make_slice_histogram(
            *data_signal_plus_bkgd, slice, matrix_map.at("BNBstats").get_matrix().get());

        // Make constrained versions of the histograms
        SliceHistogram *slice_mc_plus_ext_constr= SliceHistogram::make_slice_histogram(
            *data_mc_plus_ext_constr, slice, data_covmat_constr.get());

        SliceHistogram *slice_reco_signal_plus_bkgd_constr = SliceHistogram::make_slice_histogram(
            *data_signal_plus_bkgd_constr, slice, matrix_map_constr.at("BNBstats").get_matrix().get() );


        // Create plot for the unconstrained covariance matrix
        TCanvas* cCov = new TCanvas;
        auto data_corrmat = util::CovarianceMatrixToCorrelationMatrix( *data_covmat );
        // Plot the TMatrixD data_covmat as a colz plot
        data_corrmat.Draw("colz");
        gStyle->SetOptStat(0); // Add this line to remove the stats box

        // Add a title to the plot
        TPaveText *pt = new TPaveText(0.1, 0.94, 0.9, 0.98, "brNDC");
        pt->AddText("Unconstrained Correlation Matrix");
        pt->Draw();

        std::string out_pdf_name_cov = "plots/corr_matrix_";
        if ( sl_idx < 10 ) out_pdf_name_cov += "0";
        out_pdf_name_cov += std::to_string( sl_idx )+"_unconstrained" + nameExtension;
        cCov->SaveAs((out_pdf_name_cov + ".pdf").c_str());
        cCov->SaveAs((out_pdf_name_cov + ".png").c_str());
        cCov->SaveAs((out_pdf_name_cov + ".C").c_str());

        //  Create plot for the constrained covariance matrix
        TCanvas* cCovConstr = new TCanvas;
        auto data_corrmat_constr = util::CovarianceMatrixToCorrelationMatrix( *data_covmat_constr );
        // Plot the TMatrixD data_covmat as a colz plot
        data_corrmat_constr.Draw("colz");
        gStyle->SetOptStat(0); // Add this line to remove the stats box

        // Add a title to the plot
        TPaveText *ptConstr = new TPaveText(0.1, 0.94, 0.9, 0.98, "brNDC");
        ptConstr->AddText("Constrained Correlation Matrix");
        ptConstr->Draw();

        std::string out_pdf_name_cov_constr = "plots/corr_matrix_";
        if ( sl_idx < 10 ) out_pdf_name_cov_constr += "0";
        out_pdf_name_cov_constr += std::to_string( sl_idx )+"_constrained" + nameExtension;
        cCovConstr->SaveAs((out_pdf_name_cov_constr + ".pdf").c_str());
        cCovConstr->SaveAs((out_pdf_name_cov_constr + ".png").c_str());
        cCovConstr->SaveAs((out_pdf_name_cov_constr + ".C").c_str());


        // Get chi2
        const auto chi2 = slice_reco_signal_plus_bkgd->get_chi2(*slice_mc_plus_ext);
        const auto chi2Const = slice_reco_signal_plus_bkgd_constr->get_chi2(*slice_mc_plus_ext_constr);

        std::cout << "DEBUG constrained_slice_plots Point 4" << std::endl;
        if(normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext);
        if(normaliseByBinWidth) scale_by_bin_width(slice_reco_signal_plus_bkgd);
        if(normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext_constr);
        if(normaliseByBinWidth) scale_by_bin_width(slice_reco_signal_plus_bkgd_constr);

        TCanvas* c = new TCanvas;
        // Set the color and line thickness of the histograms
        slice_mc_plus_ext->hist_->SetLineColor(kOrange);
        slice_mc_plus_ext->hist_->SetLineWidth(4);
        slice_reco_signal_plus_bkgd->hist_->SetLineColor(kBlack);
        slice_reco_signal_plus_bkgd->hist_->SetLineWidth(1);

        slice_mc_plus_ext_constr->hist_->SetLineColor(kGreen);
        slice_mc_plus_ext_constr->hist_->SetLineWidth(2);
        // slice_mc_plus_ext_constr->hist_->SetLineStyle(2);
        // slice_reco_signal_plus_bkgd_constr->hist_->SetLineColor(kRed);
        // slice_reco_signal_plus_bkgd_constr->hist_->SetLineStyle(2);

        // Set the minimum value of the y-axis to 0
        slice_mc_plus_ext->hist_->SetMinimum(0.0);
        slice_reco_signal_plus_bkgd->hist_->SetMinimum(0.0);
        slice_mc_plus_ext_constr->hist_->SetMinimum(0.0);
        // slice_reco_signal_plus_bkgd_constr->hist_->SetMinimum(0.0);

        // Create copies of the histograms
        TH1* slice_mc_plus_ext_hist_copy = (TH1*)slice_mc_plus_ext->hist_->Clone();
        TH1* slice_mc_plus_ext_constr_hist_copy = (TH1*)slice_mc_plus_ext_constr->hist_->Clone();

        // Apply fill styles to the copies
        slice_mc_plus_ext_hist_copy->SetFillStyle(3004);
        slice_mc_plus_ext_hist_copy->SetFillColor(kOrange);
        slice_mc_plus_ext_constr_hist_copy->SetFillStyle(3005);
        slice_mc_plus_ext_constr_hist_copy->SetFillColor(kGreen);

        // Draw the histograms
        slice_mc_plus_ext_hist_copy->Draw("E2");
        slice_mc_plus_ext->hist_->Draw("hist same");
        slice_mc_plus_ext_constr_hist_copy->Draw("E2 same");
        slice_mc_plus_ext_constr->hist_->Draw("hist same");
        slice_reco_signal_plus_bkgd->hist_->Draw("E1 same");

        // slice_reco_signal_plus_bkgd_constr->hist_->Draw("E hist same");

        // Convert chi2 values to strings
        std::string chi2Str = " - Chi2: " + std::to_string(chi2.chi2_);
        std::string chi2ConstStr = " - Chi2 (constrained): " + std::to_string(chi2Const.chi2_);

        // Add a legend
        TLegend* legendSlice = new TLegend(0.1,0.7,0.48,0.9);
        legendSlice->AddEntry(slice_mc_plus_ext->hist_.get(), ("MC+EXT" + chi2Str).c_str(), "l");
        legendSlice->AddEntry(slice_mc_plus_ext_constr->hist_.get(), ("MC+EXT (constr)" + chi2ConstStr).c_str(), "l");
        legendSlice->AddEntry(slice_reco_signal_plus_bkgd->hist_.get(), "Data", "l");
        // legendSlice->AddEntry(slice_reco_signal_plus_bkgd_constr->hist_.get(), "Data (constr)", "l");
        legendSlice->Draw();

        std::string out_pdf_name = "plots/constrained_plots_slice_";
        if ( sl_idx < 10 ) out_pdf_name += "0";
        out_pdf_name += std::to_string( sl_idx ) + nameExtension;
        c->SaveAs((out_pdf_name + ".pdf").c_str());
        c->SaveAs((out_pdf_name + ".png").c_str());
        c->SaveAs((out_pdf_name + ".C").c_str());
        std::cout << "########################## Saved " << out_pdf_name << " ##########################" << std::endl;
        delete c;

        // Plot the fractional uncertainties without and with constraint applied
        std::cout << "DEBUG constrained_slice_plots Point 5" << std::endl;
        plot_frac_slice(slice, matrix_map, reco_mc_plus_ext_hist, sl_idx, nameExtension+"_without_constr");
        std::cout << "DEBUG constrained_slice_plots Point 6" << std::endl;
        plot_frac_slice(slice, matrix_map_constr, reco_mc_plus_ext_hist_constr, sl_idx, nameExtension+"_with_constr");
        std::cout << "DEBUG constrained_slice_plots Point 7" << std::endl;
    } // End of slice loop


    // ******************************************************************************************************************
    // Plot the histograms for the entire reco space ********************************************************************
    // ******************************************************************************************************************
    // Get the total covariance matrices for the unconstrained and constrained cases
    auto* cov_mat = matrix_map.at( "total" ).cov_matrix_.get();
    auto* cov_mat_constr = matrix_map_constr.at( "total" ).cov_matrix_.get();

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

    TH1D* reco_data_hist = dynamic_cast< TH1D* >(
        syst_constr.data_hists_.at( NFT::kOnBNB )->Clone( "reco_data_hist" )
    );

    // TH1D* reco_ext_hist = syst_constr.data_hists_.at( NFT::kExtBNB ).get();
    const auto& cv_univ = syst_constr.cv_universe();
    int num_reco_bins = reco_data_hist->GetNbinsX();

    std::cout << "DEBUG constrained_slice_plots Point 8" << std::endl;
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

    // reco_mc_and_ext_hist_no_constr serves as a safety check that
    // reco_mc_and_ext_hist uncertainties do not change when set the same way as reco_constrained_hist
    TH1D* reco_mc_and_ext_hist_no_constr = dynamic_cast< TH1D* >(
        reco_mc_and_ext_hist->Clone( "reco_mc_and_ext_hist_no_constr" )
    );

    // Get the post-constraint event counts and covariance matrix in the signal region
    for ( int rb = 0; rb < num_reco_bins; ++rb ) {
        if(cov_mat_constr->GetBinContent(rb + 1, rb + 1) < 0 || cov_mat->GetBinContent(rb + 1, rb + 1) < 0)
            throw std::runtime_error("Negative diagonal covariance matrix element");

        // const double err_constr = std::sqrt(
        //     std::max( 0., cov_mat_constr->GetBinContent(rb + 1, rb + 1) )
        // );
        // reco_mc_and_ext_hist->SetBinError( rb + 1, err_constr );

        const double err = std::sqrt(
            std::max( 0., cov_mat->GetBinContent(rb + 1, rb + 1) )
        );
        reco_mc_and_ext_hist_no_constr->SetBinError( rb + 1, err );

        if ( rb >= num_ordinary_reco_bins )
        {
            double data_evts = reco_data_hist->GetBinContent( rb + 1 );
            reco_constrained_hist->SetBinContent( rb + 1, data_evts );
            reco_constrained_hist->SetBinError( rb + 1, 0. );
            // reco_mc_and_ext_hist_no_constr->SetBinContent( rb + 1, data_evts );
            // reco_mc_and_ext_hist_no_constr->SetBinError( rb + 1, 0 );
        }
        else 
        {
            double constr_pred = meas_constr.reco_mc_plus_ext_->operator()( rb, 0 );
            double constr_err = std::sqrt(
            std::max( 0., meas_constr.cov_matrix_->operator()(rb, rb) )
            );

            reco_constrained_hist->SetBinContent( rb + 1, constr_pred );
            reco_constrained_hist->SetBinError( rb + 1, constr_err );

            double pred = meas.reco_mc_plus_ext_->operator()( rb, 0 );
            double err = std::sqrt(
            std::max( 0., meas.cov_matrix_->operator()(rb, rb) )
            );
            reco_mc_and_ext_hist_no_constr->SetBinContent( rb + 1, pred );
            reco_mc_and_ext_hist_no_constr->SetBinError( rb + 1, err );
        }
    }


    // ******************************************************************************************************************
    // Plot the histograms for sideband region **************************************************************************
    // ******************************************************************************************************************
    TCanvas* c2 = new TCanvas;

    reco_data_hist->SetLineColor( kBlack );
    reco_data_hist->SetLineWidth( 1 );

    // reco_mc_and_ext_hist->SetLineColor( kRed );
    // reco_mc_and_ext_hist->SetLineStyle( 2 );
    // reco_mc_and_ext_hist->SetLineWidth( 3 );

    reco_mc_and_ext_hist_no_constr->SetLineColor( kOrange );
    // reco_mc_and_ext_hist_no_constr->SetLineStyle( 2 );
    reco_mc_and_ext_hist_no_constr->SetLineWidth( 4 );

    reco_constrained_hist->SetLineColor( kGreen );
    // reco_constrained_hist->SetLineStyle( 9 );
    reco_constrained_hist->SetLineWidth( 2 );


    // Create copies of the histograms
    TH1* reco_mc_and_ext_hist_no_constr_copy = (TH1*)reco_mc_and_ext_hist_no_constr->Clone();
    TH1* reco_constrained_hist_copy = (TH1*)reco_constrained_hist->Clone();

        // Set the x-axis range to the first num_ordinary_reco_bins bins
    reco_data_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    // reco_mc_and_ext_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    reco_constrained_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    reco_mc_and_ext_hist_no_constr->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    reco_mc_and_ext_hist_no_constr_copy->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);
    reco_constrained_hist_copy->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 2*num_ordinary_reco_bins);

    // Apply fill styles to the copies
    reco_mc_and_ext_hist_no_constr_copy->SetFillStyle(3004);
    reco_mc_and_ext_hist_no_constr_copy->SetFillColor(kOrange);
    reco_constrained_hist_copy->SetFillStyle(3005);
    reco_constrained_hist_copy->SetFillColor(kGreen);

    // Draw the histograms
    reco_mc_and_ext_hist_no_constr_copy->Draw("e2");
    reco_mc_and_ext_hist_no_constr->Draw("same hist");
    // reco_constrained_hist_copy->Draw("same e2");
    reco_constrained_hist->Draw("same hist");
    reco_data_hist->Draw("same E1");

    TLegend* lg2 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    lg2->AddEntry( reco_data_hist, using_fake_data ? "fake data" : "data", "l" );
    lg2->AddEntry( reco_mc_and_ext_hist_no_constr, "uB tune + EXT", "l" );
    lg2->AddEntry( reco_constrained_hist, "uB tune + EXT (constrained)", "l" );
    lg2->Draw( "same" );
    c2->SaveAs(("plots/constrained_plots_reco_sideband" + nameExtension + ".pdf").c_str());
    c2->SaveAs(("plots/constrained_plots_reco_sideband" + nameExtension + ".png").c_str());
    c2->SaveAs(("plots/constrained_plots_reco_sideband" + nameExtension + ".C").c_str());
    std::cout << "########################## Saved " << "plots/constrained_plots_reco_sideband" + nameExtension + ".pdf" << " ##########################" << std::endl;


    // ******************************************************************************************************************
    // Plot the histograms for signal region ****************************************************************************
    // ******************************************************************************************************************
    TCanvas* c3 = new TCanvas;

    // reco_data_hist->SetLineColor( kBlack );
    // reco_data_hist->SetLineWidth( 1 );

    // // reco_mc_and_ext_hist->SetLineColor( kGreen );
    // // reco_mc_and_ext_hist->SetLineStyle( 2 );
    // // reco_mc_and_ext_hist->SetLineWidth( 3 );

    // reco_mc_and_ext_hist_no_constr->SetLineColor( kOrange );
    // // reco_mc_and_ext_hist_no_constr->SetLineStyle( 2 );
    // reco_mc_and_ext_hist_no_constr->SetLineWidth( 4 );

    // reco_constrained_hist->SetLineColor( kGreen );
    // // reco_constrained_hist->SetLineStyle( 9 );
    // reco_constrained_hist->SetLineWidth( 2 );

    // Set the x-axis range to the first num_ordinary_reco_bins bins
    // reco_mc_and_ext_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    // reco_mc_and_ext_hist_no_constr->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_constrained_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_constrained_hist_copy->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_data_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_mc_and_ext_hist_no_constr->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_mc_and_ext_hist_no_constr_copy->GetXaxis()->SetRange(1, num_ordinary_reco_bins);

    // // Create copies of the histograms
    // TH1* reco_mc_and_ext_hist_no_constr_copy = (TH1*)reco_mc_and_ext_hist_no_constr->Clone();
    // TH1* reco_constrained_hist_copy = (TH1*)reco_constrained_hist->Clone();

    // // Apply fill styles to the copies
    // reco_mc_and_ext_hist_no_constr_copy->SetFillStyle(3004);
    // reco_mc_and_ext_hist_no_constr_copy->SetFillColor(kOrange);
    // reco_constrained_hist_copy->SetFillStyle(3005);
    // reco_constrained_hist_copy->SetFillColor(kGreen);

    // Draw the histograms
    reco_mc_and_ext_hist_no_constr_copy->Draw("e2");
    reco_mc_and_ext_hist_no_constr->Draw("same hist");
    reco_constrained_hist_copy->Draw("same e2");
    reco_constrained_hist->Draw("same hist");
    reco_data_hist->Draw("same E1");

    TLegend* lg3 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    lg3->AddEntry( reco_data_hist, using_fake_data ? "fake data" : "data", "l" );
    lg3->AddEntry( reco_mc_and_ext_hist_no_constr, "uB tune + EXT", "l" );
    // lg3->AddEntry( reco_mc_and_ext_hist_no_constr, "uB tune + EXT (safety check)", "l" );
    lg3->AddEntry( reco_constrained_hist, "uB tune + EXT (constrained)", "l" );
    lg3->Draw( "same" );


    c3->SaveAs(("plots/constrained_plots_reco" + nameExtension + ".pdf").c_str());
    c3->SaveAs(("plots/constrained_plots_reco" + nameExtension + ".png").c_str());
    c3->SaveAs(("plots/constrained_plots_reco" + nameExtension + ".C").c_str());
    std::cout << "########################## Saved " << "plots/constrained_plots_reco" + nameExtension + ".pdf" << " ##########################" << std::endl;

    delete data_signal_plus_bkgd, data_signal_plus_bkgd_constr, syst_ptr, syst_ptr_constr, sb_ptr;
}

int constrained_slice_plots()
{
    slice_plots(false);
    std::cout << "------------All Done------------" << std::endl;
    return 0;
}
