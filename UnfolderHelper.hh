#pragma once

struct GeneratorInfo {
    std::string path;
    int lineColor;
    int lineStyle;
    int lineWidth;
    std::string name;
    float scaling;
};

class UnfolderHelper {
public:
    // Define your custom labels and intervals
    static const Int_t n = 8;
    static const std::vector<int> bins;
    static const Int_t overUnderFlowN;
    static const std::vector<int> overUnderFlowBins;
    static const Char_t *labels[n];

    // Existing methods...

    static void fractional_uncertainty_plot(
        const Slice& slice, 
        const std::map<std::string, std::shared_ptr<TMatrixD>> unfolded_cov_matrix_map,
        const TMatrixD* unfolded_signal,
        int sl_idx,
        std::vector<std::string> cov_mat_keys,
        const std::string& totalName,
        const std::string& nameExtension
    );

    static void event_rate_plot(const TH1* hist_event_rate_original, const size_t sl_idx, const std::string& postfix);
    static void plot_entire_matrix(const TH2D* hist, const std::string title, const std::string xTitle, const std::string yTitle, const std::string name, const bool set_y_range = true, const float y_min = -1 , const float y_max = 1);
};

// Static member initialization

// With overflow bins
const std::vector<int> UnfolderHelper::bins = {0, 11, 26, 32, 39, 49, 54, 61, 62};
const std::vector<int> UnfolderHelper::overUnderFlowBins = {31, 53};
const Int_t UnfolderHelper::overUnderFlowN = 2;

// For thee bin + slice config without overflow bins
// const std::vector<int> UnfolderHelper::bins = {0, 11, 26, 31, 38, 48, 52, 59, 60};
// const std::vector<int> UnfolderHelper::overUnderFlowBins = {}; 
// const Int_t UnfolderHelper::overUnderFlowN = 0;

const Char_t* UnfolderHelper::labels[UnfolderHelper::n] = {"cos(#theta_{#mu})", "#phi_{#mu}", "p_{#mu}", "cos(#theta_{#pi})", "#phi_{#pi}", "p_{#pi}^{**}", "#theta_{#pi #mu}", "Total"};

void drawHistogramWithBand(TH1* hist, const int lineColor, const int lineWidth, const int lineStyle, const float alpha, const bool noBand = false)
{
    // if (!hist_base) {
    //     throw std::runtime_error("Error in drawHistogramWithBand: hist_base is nullptr");
    // }

    // Remove const qualifier since Clone() produces a modifiable copy
    // TH1D* hist = dynamic_cast<TH1D*>(hist_base->Clone());
    if (!hist) {
        throw std::runtime_error("Error in drawHistogramWithBand: dynamic_cast failed, possibly due to incorrect histogram type");
    }


    if(!noBand)
    {
        // Clone the histogram to draw the error band
        TH1D* hist_band = dynamic_cast<TH1D*>(hist->Clone(Form("%s_band", hist->GetName())));

        // Set the properties for the error band
        hist_band->SetFillColorAlpha(lineColor, alpha); // Set the fill color to be semi-transparent

        // Draw the error band
        hist_band->Draw("E2 same");
    }

    // Print out the bin values and errors
    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        double binValue = hist->GetBinContent(i);
        double binError = hist->GetBinError(i);
        std::cout << "Bin " << i << ": Value = " << binValue << ", Error = " << binError << std::endl;
    }

    // Set the properties for the line
    hist->SetStats(false);
    hist->SetLineColor(lineColor); // Set the line color
    hist->SetLineWidth(lineWidth); // Set the line width
    hist->SetLineStyle(lineStyle); // Set the line style

    // Draw the line over the error band
    hist->Draw("HIST same");
}


void drawHistogramWithBand(TH1* hist, GeneratorInfo genInfo, float alpha, const bool noBand = false)
{
    drawHistogramWithBand(hist, genInfo.lineColor, genInfo.lineWidth, genInfo.lineStyle, alpha, noBand);
}

void UnfolderHelper::fractional_uncertainty_plot(
    const Slice& slice, 
    const std::map<std::string, std::shared_ptr<TMatrixD>> unfolded_cov_matrix_map,
    const TMatrixD* unfolded_signal,
    int sl_idx,
    std::vector< std::string > cov_mat_keys,
    const std::string& totalName,
    const std::string& nameExtension
)
{
    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    TH1* slice_hist = dynamic_cast< TH1* >(
      slice.hist_->Clone("slice_hist") );

    slice_hist->SetDirectory( nullptr );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Loop over the various systematic uncertainties
    int color = 0;
    for ( const auto& pair : unfolded_cov_matrix_map ) {
        const auto& key = pair.first;
        const auto& cov_matrix_ptr = pair.second.get();

        SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
            *unfolded_signal, slice, cov_matrix_ptr );

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
        {
            // std::cout << "DEBUG skipping " << key << std::endl;
            continue;
        }
        // else
        // {
        //     std::cout<<"DEBUG not skipping "<<key<<std::endl;
        // }

        frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

        if ( color <= 9 ) ++color;
        if ( color == 5 ) ++color;
        if ( color >= 10 ) color += 10;

        slice_for_syst->hist_->SetLineColor( color );
        slice_for_syst->hist_->SetLineWidth( 3 );
    }


    TCanvas* c2 = new TCanvas;
    // Set right padding to allow for the legend
    c2->SetRightMargin( 0.21 );

    // TLegend* lg2 = new TLegend( 0.2, 0.3);
    TLegend* lg2 = new TLegend( 0.8, 0.1, 0.95, 0.9); // x1, y1, x2, y2


    auto* total_frac_err_hist = frac_uncertainty_hists.at( totalName );
    // total_frac_err_hist->SetTitle("Fractional Uncertainty of Selected #nu_{#mu}CC1#pi^{#pm}Xp, X #geq 0 Events");
    total_frac_err_hist->SetTitle("");
    total_frac_err_hist->SetTitle("");
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
    total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineStyle( 9 );
    total_frac_err_hist->SetLineWidth( 3 );
    total_frac_err_hist->Draw( "hist" );
    total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");

    // const auto frac_ymax = 0.35;
    // total_frac_err_hist->GetYaxis()->SetRangeUser( 0., frac_ymax);

    lg2->AddEntry( total_frac_err_hist, totalName.c_str(), "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the totalName one above
      if ( name == totalName ) continue;
      if (name.size() >= 5 && name.substr(name.size() - 5) == "stats") 
      {
        hist->SetLineStyle( 2 );
      }

      lg2->AddEntry( hist, name.c_str(), "l" );
      hist->Draw( "same hist" );


    //   std::cout << name << " frac err in bin #1 = "
    //     << hist->GetBinContent( 1 )*100. << "%\n";
    }

    lg2->Draw( "same" );

    std::string frac_out_pdf_name = "plots/plot_frac_slice_";
    if ( sl_idx < 10 ) frac_out_pdf_name += "0";
    frac_out_pdf_name += std::to_string( sl_idx ) + nameExtension +".pdf";
    c2->SaveAs( frac_out_pdf_name.c_str() );
}

// void UnfolderHelper::fractional_uncertainty_plot(
//     const Slice& slice, 
//     const std::map<std::string, std::shared_ptr<TMatrixD>> unfolded_cov_matrix_map,
//     const TMatrixD* unfolded_signal,
//     int sl_idx,
//     const std::string& nameExtension
// )
// {
//     // Get the binning and axis labels for the current slice by cloning the
//     // (empty) histogram owned by the Slice object
//     TH1* slice_hist = dynamic_cast< TH1* >(
//       slice.hist_->Clone("slice_hist") );

//     slice_hist->SetDirectory( nullptr );

//     // Keys are labels, values are fractional uncertainty histograms
//     auto* fr_unc_hists = new std::map< std::string, TH1* >();
//     auto& frac_uncertainty_hists = *fr_unc_hists;

//     // Show fractional uncertainties computed using these covariance matrices
//     // in the ROOT plot. All configured fractional uncertainties will be
//     // included in the output pgfplots file regardless of whether they appear
//     // in this vector.
    
//     // std::vector< std::string > cov_mat_keys = { "total", "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets", "MCstats", "EXTstats"};

//     std::vector< std::string > cov_mat_keys = { "total", "xsec_total", "MCstats", "BNBstats"};//, "EXTstats", "BNBstats"};
//     // cov_mat_keys = { "total", "MCstats", "xsec_multi", "xsec_unisim", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie"}; // removed ext + beam on stats
//     // cov_mat_keys = { "total", "MCstats", "EXTstats", "BNBstats", "xsec_multi", "xsec_unisim", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie"};
//     // cov_mat_keys = { "total", "xsec_multi", "xsec_AxFFCCQEshape", "xsec_DecayAngMEC", "xsec_NormCCCOH", "xsec_NormNCCOH", "xsec_RPA_CCQE", "xsec_ThetaDelta2NRad", "xsec_Theta_Delta2Npi", "xsec_VecFFCCQEshape", "xsec_XSecShape_CCMEC", "xsec_xsr_scc_Fa3_SCC", "xsec_xsr_scc_Fv3_SCC", "NuWroGenie", "MCstats", "EXTstats", "BNBstats"};

//     // Loop over the various systematic uncertainties
//     int color = 0;
//     for ( const auto& pair : unfolded_cov_matrix_map ) {
//         const auto& key = pair.first;
//         const auto& cov_matrix_ptr = pair.second.get();

//         SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
//             *unfolded_signal, slice, cov_matrix_ptr );

//         // The SliceHistogram object already set the bin errors appropriately
//         // based on the slice covariance matrix. Just change the bin contents
//         // for the current histogram to be fractional uncertainties. Also set
//         // the "uncertainties on the uncertainties" to zero.
//         // TODO: revisit this last bit, possibly assign bin errors here
//         for ( const auto& bin_pair : slice.bin_map_ ) {
//             int global_bin_idx = bin_pair.first;
//             double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
//             double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
//             double frac = 0.;
//             if ( y > 0. ) frac = err / y;
//             slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
//             slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
//         }


//         // Check whether the current covariance matrix name is present in
//         // the vector defined above this loop. If it isn't, don't bother to
//         // plot it, and just move on to the next one.
//         auto cbegin = cov_mat_keys.cbegin();
//         auto cend = cov_mat_keys.cend();
//         auto iter = std::find( cbegin, cend, key );
//         if ( iter == cend )
//         {
//             std::cout << "DEBUG skipping " << key << std::endl;
//             continue;
//         }
//         else
//         {
//             std::cout<<"DEBUG not skipping "<<key<<std::endl;
//         }

//         frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

//         if ( color <= 9 ) ++color;
//         if ( color == 5 ) ++color;
//         if ( color >= 10 ) color += 10;

//         slice_for_syst->hist_->SetLineColor( color );
//         slice_for_syst->hist_->SetLineWidth( 3 );
//     }


//     TCanvas* c2 = new TCanvas;
//     // c2->SetLogy(); // Use this for golden Pion Cut variable plots
//     // TLegend* lg2 = new TLegend( 0.7, 0.7, 0.9, 0.9 );
//     TLegend* lg2 = new TLegend( 0.2, 0.3);


//     auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
//     total_frac_err_hist->SetStats( false );
//     total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
//     total_frac_err_hist->GetMaximum() * 1.05 );
//     total_frac_err_hist->SetLineColor( kBlack );
//     total_frac_err_hist->SetLineStyle( 9 );
//     total_frac_err_hist->SetLineWidth( 3 );
//     total_frac_err_hist->Draw( "hist" );
//     total_frac_err_hist->SetTitle("Fractional Uncertainty of Selected #nu_{#mu}CC1#pi^{#pm}Xp, X #geq 0 Events");
//     total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");


//     // const auto frac_ymax = 0.35;
//     // total_frac_err_hist->GetYaxis()->SetRangeUser( 0., frac_ymax);

//     lg2->AddEntry( total_frac_err_hist, "total", "l" );

//     for ( auto& pair : frac_uncertainty_hists ) {
//       const auto& name = pair.first;
//       TH1* hist = pair.second;
//       // We already plotted the "total" one above
//       if ( name == "total" ) continue;
//       if (name.size() >= 5 && name.substr(name.size() - 5) == "stats") 
//       {
//         hist->SetLineStyle( 2 );
//       }

//       lg2->AddEntry( hist, name.c_str(), "l" );
//       hist->Draw( "same hist" );


//       std::cout << name << " frac err in bin #1 = "
//         << hist->GetBinContent( 1 )*100. << "%\n";
//     }

//     lg2->Draw( "same" );

//     std::string frac_out_pdf_name = "plots/plot_unfolded_frac_slice_";
//     if ( sl_idx < 10 ) frac_out_pdf_name += "0";
//     frac_out_pdf_name += std::to_string( sl_idx ) + nameExtension +".pdf";
//     c2->SaveAs( frac_out_pdf_name.c_str() );
// }

void UnfolderHelper::event_rate_plot(const TH1 *hist_event_rate_original, const size_t sl_idx, const std::string &postfix)
{
    // Clone hist_event_rate_original to avoid modifying the original histogram
    TH1* hist_event_rate = dynamic_cast<TH1*>(hist_event_rate_original->Clone("hist_event_rate"));

    TCanvas* c1_event_rate = new TCanvas((postfix + std::to_string(sl_idx)).c_str(), (postfix + std::to_string(sl_idx)).c_str(), 800, 600);

    hist_event_rate->SetLineColor(kBlack);
    // hist_event_rate->SetLineWidth(5);
    hist_event_rate->SetMarkerStyle(kFullCircle);
    hist_event_rate->SetMarkerSize(0.7);
    hist_event_rate->SetStats(false);
    hist_event_rate->SetLineWidth(2);


    double ymax = -DBL_MAX;
    hist_event_rate->Draw("e");


    if(sl_idx == 7)
    {
        hist_event_rate->GetXaxis()->SetLabelOffset(999); // Hide x-axis labels
        hist_event_rate->GetXaxis()->SetTickLength(0); // Hide x-axis ticks
    }

    // for (const auto &pair : slice_gen_map)
    // {
    //     const auto &name = pair.first;
    //     const auto *slice_h = pair.second;

    //     double max = slice_h->hist_->GetMaximum();
    //     // if (max > ymax)
    //     //     ymax = max;

    //     for (int i = 1; i <= slice_h->hist_->GetNbinsX(); ++i) {
    //         double binContent = slice_h->hist_->GetBinContent(i);
    //         double binError = (name == "Unfolded Selection") ? slice_h->hist_->GetBinError(i) : 0;
    //         if (binContent + binError > ymax) {
    //             ymax = binContent + binError;
    //         }
    //     }

    //     if (name == "Unfolded Selection" || name == "NuWro Truth" || name == "MicroBooNE Tune")
    //         continue;

    //     const auto &file_info = truth_file_map.at(name);
    //     slice_h->hist_->SetLineColor(file_info.color_);
    //     slice_h->hist_->SetLineStyle(file_info.style_);
    //     slice_h->hist_->SetLineWidth(3);
    //     slice_h->hist_->Draw("hist same");
    // }

    // drawHistogramWithBand(slice_cv->hist_.get(), kAzure - 7, 2, 5, 0.3, true);

    // if (using_fake_data)
    // {
    //     const auto& slice_truth_cov_mat_ptr = slice_truth->cmat_.get_matrix(); // Get the pointer first
    //     if (!slice_truth_cov_mat_ptr) 
    //     {
    //         throw std::runtime_error("Error: slice_truth_cov_mat_ptr is nullptr");
    //     }
    //     const auto& slice_truth_cov_mat = *slice_truth_cov_mat_ptr; // Now dereference

    //     //Set the slice_truth histogram errors to the sqrt of the diagonal slice_truth->cmat_ values
    //     for (int i = 1; i <= slice_truth->hist_->GetNbinsX(); ++i) {
    //         double error = sqrt(slice_truth_cov_mat(i-1, i-1));
    //         std::cout<<"DEBUG slice_truth->hist_->GetBinContent(i): "<<slice_truth->hist_->GetBinContent(i)<<" error^2: "<< slice_truth_cov_mat(i-1, i-1) <<" error: "<<error<<std::endl;
    //         slice_truth->hist_->SetBinError(i, error);
    //     }

    //     drawHistogramWithBand(slice_truth->hist_.get(), kGreen, 2, 5, 0.3, false);
    // }


    // hist_event_rate->GetYaxis()->SetRangeUser(0., ymax * 1.05);
    hist_event_rate->GetYaxis()->SetRangeUser(0., hist_event_rate->GetMaximum() * 1.1);
    hist_event_rate->Draw("E same");
    // hist_event_rate->SetTitle("Unfolded NuWro CC1#pi^{#pm}Xp");
    hist_event_rate->SetTitle(""); // No title
    hist_event_rate->GetXaxis()->SetLabelOffset(999); // Hide x-axis
    hist_event_rate->GetXaxis()->SetTitleOffset(999); // Hide x-axis title

    std::string out_pdf_name = "plots/plot_unfolded_event_rate_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + postfix + ".pdf";
    c1_event_rate->SaveAs(out_pdf_name.c_str());
}



void UnfolderHelper::plot_entire_matrix(const TH2D* hist_init, const std::string title, const std::string xTitle, const std::string yTitle, const std::string name, const bool set_y_range, const float y_min, const float y_max) {
    // Clone the histogram to avoid modifying the original histogram
    TH2D* hist = dynamic_cast<TH2D*>(hist_init->Clone("hist"));

    // Set the color palette
    util::CreateRedToBlueColorPalette(20);
    gStyle->SetTitleFontSize(0.05); // Set the title font size to 0.05

    TCanvas *cm2 = new TCanvas((title+name).c_str(), (title+name).c_str(), 800, 600);
    hist->SetTitleSize(0.05, "t"); // Set the title font size to 0.05
    hist->SetTitle(title.c_str());
    if(set_y_range) hist->GetZaxis()->SetRangeUser(y_min, y_max);
    hist->GetXaxis()->SetRangeUser(0, 66); // Set the range of the x-axis
    hist->GetXaxis()->SetTitle(xTitle.c_str()); // Set x-axis label
    hist->GetYaxis()->SetTitle(yTitle.c_str()); // Set y-axis label
    hist->SetStats(kFALSE);
    hist->Draw("COLZ");

    // Draw vertical and horizontal lines at the bin edges
    for (Int_t i = 1; i < UnfolderHelper::n; i++) {
        TLine *vline = new TLine(UnfolderHelper::bins[i], 0, UnfolderHelper::bins[i], hist->GetNbinsY());
        vline->SetLineColor(kBlack);
        vline->Draw();

        TLine *hline = new TLine(0, UnfolderHelper::bins[i], hist->GetNbinsX(), UnfolderHelper::bins[i]);
        hline->SetLineColor(kBlack);
        hline->Draw();
    }

    // Draw white dotted lines for visual separation
    // Note: The original logic for dotted lines seems incorrect. Adjusting to draw between bins.
    for (Int_t i = 1; i < UnfolderHelper::n; i++) {
        TLine *vline_dotted = new TLine(UnfolderHelper::bins[i], 0, UnfolderHelper::bins[i], hist->GetNbinsY());
        vline_dotted->SetLineColor(kWhite);
        vline_dotted->SetLineStyle(2); // Set line style to dotted
        vline_dotted->Draw();

        TLine *hline_dotted = new TLine(0, UnfolderHelper::bins[i], hist->GetNbinsX(), UnfolderHelper::bins[i]);
        hline_dotted->SetLineColor(kWhite);
        hline_dotted->SetLineStyle(2); // Set line style to dotted
        hline_dotted->Draw();
    }

    // Add labels in the middle of the intervals
    for (Int_t i = 0; i < UnfolderHelper::n - 1; i++) {
        Double_t midPoint = (UnfolderHelper::bins[i] + UnfolderHelper::bins[i+1]) / 2.0;
        TLatex *text = new TLatex(midPoint, 1.03*hist->GetNbinsY(), UnfolderHelper::labels[i]);
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->Draw();
    }

    // Add asterisk for over-/underrflow bins
    for (Int_t i = 0; i < UnfolderHelper::overUnderFlowN; i++) {
        Double_t midPoint = UnfolderHelper::overUnderFlowBins[i] + 0.5;
        TLatex *text = new TLatex(midPoint, 1.0*hist->GetNbinsY(), "*");
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->SetTextColor(kGray+3); // Set text color to grey
        text->Draw();
    }

    // Add footnote
    TLatex *footnote = new TLatex(0, -0.1*hist->GetNbinsY(), "* Under-/overflow bin; ** Selection subset");
    footnote->SetTextSize(0.02); // Set text size to something smaller
    footnote->SetTextColor(kGray+3);
    footnote->Draw();

    cm2->SaveAs(("plots/plot_entire_" + name + ".pdf").c_str());
}