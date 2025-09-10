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
    static const std::vector<std::string> labels;

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
    static void plot_entire_matrix(const TH2D* hist, const std::string title, const std::string xTitle, const std::string yTitle, const std::string name, 
        const bool setYRange = true, const float y_min = -1 , const float y_max = 1, const std::vector<int> binsX = UnfolderHelper::bins, const std::vector<int> binsY = UnfolderHelper::bins, 
        const std::vector<std::string> labelsX = UnfolderHelper::labels, const std::vector<std::string> labelsY = UnfolderHelper::labels, 
        const bool addText = false, const int xOffset = 0, const int yOffset = 0, const bool footnote = true, const bool differentStyle = false, const bool scientificNotation = false);
};

// Static member initialization

// With overflow bins
// const std::vector<int> UnfolderHelper::bins = {0, 11, 26, 32, 39, 49, 54, 61, 62};
// const std::vector<int> UnfolderHelper::overUnderFlowBins = {31, 53};
// const Int_t UnfolderHelper::overUnderFlowN = 2;

// For thee bin + slice config without overflow bins
const std::vector<int> UnfolderHelper::bins = {0, 11, 26, 31, 38, 48, 52, 59, 60};
const std::vector<int> UnfolderHelper::overUnderFlowBins = {}; 
const Int_t UnfolderHelper::overUnderFlowN = 0;

// const Char_t* UnfolderHelper::labels[UnfolderHelper::n] = {"cos(#theta_{#mu})", "#phi_{#mu}", "p_{#mu}", "cos(#theta_{#pi})", "#phi_{#pi}", "p_{#pi}^{**}", "#theta_{#pi #mu}", "Total"};
const std::vector<std::string> UnfolderHelper::labels = {"cos(#theta_{#mu})", "#phi_{#mu}", "p_{#mu}", "cos(#theta_{#pi})", "#phi_{#pi}", "p_{#pi}^{**}", "#theta_{#pi #mu}", "Total"};

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
    std::vector<std::string> cov_mat_keys,
    const std::string& totalName,
    const std::string& nameExtension
) {
    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    TH1* slice_hist = dynamic_cast<TH1*>(slice.hist_->Clone("slice_hist"));
    slice_hist->SetDirectory(nullptr);

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map<std::string, TH1*>();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Create the legend outside the loop
    TLegend* lg2 = new TLegend(0.1, 0.1, 0.9, 0.9); // Adjust legend size to fill the canvas
    lg2->SetHeader("Uncertainties", "C"); // Add a header to the legend

    // Loop over the various systematic uncertainties
    int color = 0;
    for (const auto& pair : unfolded_cov_matrix_map) {
        const auto& key = pair.first;
        const auto& cov_matrix_ptr = pair.second.get();

        SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
            *unfolded_signal, slice, cov_matrix_ptr);

        // The SliceHistogram object already set the bin errors appropriately
        // based on the slice covariance matrix. Just change the bin contents
        // for the current histogram to be fractional uncertainties. Also set
        // the "uncertainties on the uncertainties" to zero.
        for (const auto& bin_pair : slice.bin_map_) {
            int global_bin_idx = bin_pair.first;
            double y = slice_for_syst->hist_->GetBinContent(global_bin_idx);
            double err = slice_for_syst->hist_->GetBinError(global_bin_idx);
            double frac = 0.;
            if (y > 0.) frac = err / y;
            slice_for_syst->hist_->SetBinContent(global_bin_idx, frac);
            slice_for_syst->hist_->SetBinError(global_bin_idx, 0.);
        }

        // Check whether the current covariance matrix name is present in
        // the vector defined above this loop. If it isn't, don't bother to
        // plot it, and just move on to the next one.
        if (std::find(cov_mat_keys.begin(), cov_mat_keys.end(), key) == cov_mat_keys.end()) {
            continue;
        }

        frac_uncertainty_hists[key] = slice_for_syst->hist_.get();

        if (color <= 9) ++color;
        if (color == 5) ++color;
        if (color >= 10) color += 9;

        slice_for_syst->hist_->SetLineColor(color);
        slice_for_syst->hist_->SetLineWidth(3);

        // Add to the legend only for the first slice
        if (sl_idx == 0) {
            std::string label;
            if (key == "EXTstats") label = "Beam-off Statistical";
            else if (key == "total") label = "Total";
            else if (key == "MCstats") label = "MC Statistical";
            else if (key == "BNBstats") label = "Beam-on Statistical";
            else if (key == "POT") label = "POT";
            else if (key == "detVar_total") label = "Detector Variations";
            else if (key == "flux") label = "Flux";
            else if (key == "numTargets") label = "Target";
            else if (key == "reint") label = "Reinteraction";
            else if (key == "xsec_total") label = "Interaction Model";
            else label = key; // Default to the original name if not matched

            lg2->AddEntry(slice_for_syst->hist_.get(), label.c_str(), "l");
        }
    }

    // Create the canvas for the fractional uncertainty plot
    TCanvas* c2 = new TCanvas;
    c2->SetRightMargin(0.05); // Reduce the right margin to minimize whitespace

    auto* total_frac_err_hist = frac_uncertainty_hists.at(totalName);
    total_frac_err_hist->SetTitle("");
    total_frac_err_hist->SetStats(false);
    total_frac_err_hist->GetYaxis()->SetRangeUser(0., total_frac_err_hist->GetMaximum() * 1.05);
    total_frac_err_hist->SetLineColor(kBlack);
    total_frac_err_hist->SetLineStyle(9);
    total_frac_err_hist->SetLineWidth(3);
    total_frac_err_hist->Draw("hist");
    total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");

    for (auto& pair : frac_uncertainty_hists) {
        const auto& name = pair.first;
        TH1* hist = pair.second;
        if (name == totalName) continue;
        if (name.size() >= 5 && name.substr(name.size() - 5) == "stats") {
            hist->SetLineStyle(2);
        }
        hist->Draw("same hist");
    }

    std::string frac_out_pdf_name = "plots/plot_frac_slice_";
    if (sl_idx < 10) frac_out_pdf_name += "0";
    frac_out_pdf_name += std::to_string(sl_idx) + nameExtension + ".pdf";
    c2->SaveAs(frac_out_pdf_name.c_str());

    // Save the legend as a separate figure
    if (sl_idx == 0) {
        TCanvas* legend_canvas = new TCanvas("legend_canvas", "Legend", 400, 400); // Adjust canvas size to match legend
        legend_canvas->cd();
        lg2->Draw();
        legend_canvas->SaveAs((std::string("plots/plot_frac_uncertainties_legend_") + nameExtension + ".pdf").c_str());
        legend_canvas->SaveAs((std::string("plots/plot_frac_uncertainties_legend_") + nameExtension + ".C").c_str());
    }
}

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



void UnfolderHelper::plot_entire_matrix(const TH2D* hist_init, const std::string title, const std::string xTitle, const std::string yTitle, const std::string name, const bool setYRange, const float y_min, const float y_max, 
    const std::vector<int> binsX, const std::vector<int> binsY, const std::vector<std::string> labelsX, const std::vector<std::string> labelsY, const bool addText, const int xOffset, const int yOffset, const bool footnote, const bool differentStyle, const bool scientificNotation){
    std::cout<<"DEBUG: Removed bins point 1.0"<<std::endl;
    // Clone the histogram to avoid modifying the original histogram
    TH2D* hist = dynamic_cast<TH2D*>(hist_init->Clone("hist"));

    std::cout<<"DEBUG: Removed bins point 2.0"<<std::endl;
    TCanvas *cm2 = new TCanvas((title+name).c_str(), (title+name).c_str(), 800, 600);

    // Set the color palette
    if(!differentStyle)
    {
        util::CreateRedToBlueColorPalette(20);
        gStyle->SetTitleFontSize(0.04); // Set the title font size to 0.05
    }
    else
    {
        cm2->SetLogz();
        gStyle->SetPalette(57);
    }

    cm2->SetRightMargin(0.16); // Adjust the right margin to add whitespace
    cm2->SetTopMargin(0.12); // Adjust the top margin to add whitespace
    hist->SetTitleSize(0.05, "t"); // Set the title font size to 0.05
    hist->SetTitle(title.c_str());
    if(setYRange) hist->GetZaxis()->SetRangeUser(y_min, y_max);
    // hist->GetXaxis()->SetRangeUser(0, 66); // Set the range of the x-axis
    hist->GetXaxis()->SetTitle(xTitle.c_str()); // Set x-axis label
    hist->GetYaxis()->SetTitle(yTitle.c_str()); // Set y-axis label
    hist->SetStats(kFALSE);
    hist->Draw("COLZ");

    // Adjust the z-axis label offset to create more space
    // hist->GetZaxis()->SetLabelOffset(0.01);
    // hist->GetZaxis()->SetTitleOffset(1.3);
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
    if (palette) {
        palette->SetX1NDC(0.89); // Adjust this value as needed
        palette->SetX2NDC(0.94); // Adjust this value as needed
    }

    // Set custom bin labels for x-axis
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        hist->GetXaxis()->SetBinLabel(i, std::to_string(i + xOffset).c_str());
    }
    hist->GetXaxis()->SetLabelSize(hist->GetNbinsX() > 15 ? 0.04 : 0.06); // Adjust the label size as needed
    
    // Set custom bin labels for y-axis
    for (int i = 1; i <= hist->GetNbinsY(); ++i) {
        hist->GetYaxis()->SetBinLabel(i, std::to_string(i + yOffset).c_str());
    }
    hist->GetYaxis()->SetLabelSize(hist->GetNbinsY() > 15 ? 0.04 : 0.06); // Adjust the label size as needed

    std::cout<<"DEBUG: Removed bins point 3.0"<<std::endl;

    std::cout << "Size of hist: " << hist->GetNbinsX() << " x " << hist->GetNbinsY() << std::endl;

    const int nX = binsX.size();
    const int nY = binsY.size();

    // Draw white dotted lines for visual separation
    // Note: The original logic for dotted lines seems incorrect. Adjusting to draw between bins.
    for (Int_t i = 1; i < nX-1; i++) {
        TLine *vline = new TLine(binsX[i], 0, binsX[i], hist->GetNbinsY());
        vline->SetLineColor(kBlack);
        vline->Draw();

        TLine *vline_dotted = new TLine(binsX[i], 0, binsX[i], hist->GetNbinsY());
        vline_dotted->SetLineColor(kWhite);
        vline_dotted->SetLineStyle(2); // Set line style to dotted
        vline_dotted->Draw();
    }

    for (Int_t i = 0; i < nX-1; i++) {    
        Double_t midPoint = (binsX[i] + binsX[i+1]) / 2.0;
        TLatex *text = new TLatex(midPoint, 1.035*hist->GetNbinsY(), (labelsX[i]).c_str());
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->Draw();
    }

    for (Int_t i = 1; i < nY-1; i++) {
        TLine *hline = new TLine(0, binsY[i], hist->GetNbinsX(), binsY[i]);
        hline->SetLineColor(kBlack);
        hline->Draw();

        TLine *hline_dotted = new TLine(0, binsY[i], hist->GetNbinsX(), binsY[i]);
        hline_dotted->SetLineColor(kWhite);
        hline_dotted->SetLineStyle(2); // Set line style to dotted
        hline_dotted->Draw();
    }

    for (Int_t i = 0; i < nY-1; i++) {
        Double_t midPoint = (binsY[i] + binsY[i+1]) / 2.0;
        TLatex *text = new TLatex(1.035*hist->GetNbinsX(), midPoint, (labelsY[i]).c_str());
        text->SetTextSize(0.03); // Set text size to something smaller
        text->SetTextAlign(22); // Center alignment
        text->SetTextAngle(90); // Rotate text by 90 degrees
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

    std::cout<<"DEBUG: Removed bins point 5.0"<<std::endl;

    if(footnote)
    {
        // TLatex *footnote = new TLatex(0, -0.1*hist->GetNbinsY(), "* Under-/overflow bin; ** Selection subset");
        TLatex *footnote = new TLatex(0, -0.1*hist->GetNbinsY(), "* Selection subset");
        footnote->SetTextSize(0.03);
        footnote->SetTextColor(kGray+3);
        footnote->Draw();
    }


    if(addText)
    {
        const int num_bins_X = hist->GetNbinsX();
        const int num_bins_Y = hist->GetNbinsY();
        for (int i = 0; i < num_bins_X; i++) {
            for (int j = 0; j < num_bins_Y; j++) 
            {
                const auto text = scientificNotation ? Form("%.2e", hist->GetBinContent(i+1, j+1)) : Form("%.3f", hist->GetBinContent(i+1, j+1));
                TLatex* latex = new TLatex(hist->GetXaxis()->GetBinCenter(i+1), hist->GetYaxis()->GetBinCenter(j+1), text);
                latex->SetTextFont(42);
                latex->SetTextSize(0.02);
                latex->SetTextAlign(22);
                latex->Draw();
            }
        }
    }
    std::cout<<"DEBUG: Removed bins point 6.0"<<std::endl;

    cm2->SaveAs(("plots/plot_entire_" + name + ".pdf").c_str());
    std::cout<<"DEBUG: Removed bins point 7.0"<<std::endl;
}