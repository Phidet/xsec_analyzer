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

    #ifdef USE_FAKE_DATA
        std::cout << "########################## Using fake data ##########################" << std::endl;
        // Initialize the FilePropertiesManager and tell it to treat the NuWro
        // MC ntuples as if they were data
        auto &fpm = FilePropertiesManager::Instance();
        fpm.load_file_properties("closure_file_properties_run1.txt");
    #endif

    auto *syst_ptr = new MCC9SystematicsCalculator(
        "/uboone/data/users/jdetje/ubcc1pi_univmake/univmake_output_closure_run1_noDirt_noEXT_V2.root",
        "systcalc_unfold_fd_closure.conf"); // "systcalc_unfold_fd.conf");
    auto &syst = *syst_ptr;

    // Get access to the relevant histograms owned by the SystematicsCalculator
    // object. These contain the reco bin counts that we need to populate the
    // slices below.
    TH1D *reco_bnb_hist = syst.data_hists_.at(NFT::kOnBNB).get();
    TH1D *reco_ext_hist = syst.data_hists_.at(NFT::kExtBNB).get();

    // std::cout<<"WARNING Remove scaling !!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    // reco_ext_hist->Scale(0.); // Todo remove

    #ifdef USE_FAKE_DATA
        // Add the EXT to the "data" when working with fake data
        reco_bnb_hist->Add(reco_ext_hist);
    #endif

    TH2D *category_hist = syst.cv_universe().hist_categ_.get();

    // Total MC+EXT prediction in reco bin space. Start by getting EXT.
    TH1D *reco_mc_plus_ext_hist = dynamic_cast<TH1D *>(
        reco_ext_hist->Clone("reco_mc_plus_ext_hist"));
    reco_mc_plus_ext_hist->SetDirectory(nullptr);

    // Add in the CV MC prediction
    reco_mc_plus_ext_hist->Add(syst.cv_universe().hist_reco_.get());

    const int num_bins = reco_mc_plus_ext_hist->GetNbinsX();
    // Ensure both histograms have the same number of bins
    assert(num_bins == reco_bnb_hist->GetNbinsX());
    for (int i = 1; i <= num_bins; ++i)
    {
        double reco_mc_plus_ext_value = reco_mc_plus_ext_hist->GetBinContent(i);
        double reco_bnb_value = reco_bnb_hist->GetBinContent(i);

        std::cout << "Bin " << i << ": "
                  << "reco_mc_plus_ext = " << reco_mc_plus_ext_value << ", "
                  << "reco_bnb = " << reco_bnb_value << std::endl;
    }

    // Keys are covariance matrix types, values are CovMatrix objects that
    // represent the corresponding matrices
    auto *matrix_map_ptr = syst.get_covariances().release();
    auto &matrix_map = *matrix_map_ptr;

    auto *sb_ptr = new SliceBinning("ubcc1pi_slice_config.txt");
    auto &sb = *sb_ptr;

    for (size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx)
    {
        std::cout << "DEBUG Point 0" << std::endl;
        if (sl_idx != 0u)
            continue; // todo remove !!!!!!!!!!!!!!!!!!!

        std::cout << "DEBUG Point 1" << std::endl;

        const auto &slice = sb.slices_.at(sl_idx);

        std::cout << "DEBUG Point 1.1" << std::endl;
        // We now have all of the reco bin space histograms that we need as input.
        // Use them to make new histograms in slice space.
        SliceHistogram *slice_bnb = SliceHistogram::make_slice_histogram(
            *reco_bnb_hist, slice, &matrix_map.at("BNBstats"));

        std::cout << "DEBUG Point 1.2" << std::endl;
        SliceHistogram *slice_ext = SliceHistogram::make_slice_histogram(
            *reco_ext_hist, slice, &matrix_map.at("EXTstats"));

        std::cout << "DEBUG Point 1.3" << std::endl;
        SliceHistogram *slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
            *reco_mc_plus_ext_hist, slice, &matrix_map.at("total"));

        std::cout << "DEBUG Point 2" << std::endl;

        // auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
        // std::cout << "Slice " << sl_idx << ": \u03C7\u00b2 = "
        //  << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
        //  << " p-value = " << chi2_result.p_value_ << '\n';

        // Build a stack of categorized central-value MC predictions plus the
        // extBNB contribution in slice space
        const auto &eci = EventCategoryInterpreter::Instance();
        eci.set_ext_histogram_style(slice_ext->hist_.get());

        THStack *slice_pred_stack = new THStack("mc+ext", "");
        if (normaliseByBinWidth)
            scale_by_bin_width(slice_ext);

        slice_pred_stack->Add(slice_ext->hist_.get()); // extBNB

        const auto &cat_map = eci.label_map();

        // Go in reverse so that signal ends up on top. Note that this index is
        // one-based to match the ROOT histograms
        int cat_bin_index = cat_map.size();
        for (auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter)
        {
            EventCategory cat = iter->first;
            TH1D *temp_mc_hist = category_hist->ProjectionY("temp_mc_hist",
                                                            cat_bin_index, cat_bin_index);
            temp_mc_hist->SetDirectory(nullptr);

            SliceHistogram *temp_slice_mc = SliceHistogram::make_slice_histogram(
                *temp_mc_hist, slice);

            eci.set_mc_histogram_style(cat, temp_slice_mc->hist_.get());

            if (normaliseByBinWidth)
                scale_by_bin_width(temp_slice_mc);
            slice_pred_stack->Add(temp_slice_mc->hist_.get());

            std::string cat_col_prefix = "MC" + std::to_string(cat);

            --cat_bin_index;
        }

        TCanvas *c1 = new TCanvas;
        slice_bnb->hist_->SetLineColor(kBlack);
        slice_bnb->hist_->SetLineWidth(3);
        slice_bnb->hist_->SetMarkerStyle(kFullCircle);
        slice_bnb->hist_->SetMarkerSize(0.8);
        slice_bnb->hist_->SetStats(false);
        if (normaliseByBinWidth)
            scale_by_bin_width(slice_bnb);
        if (normaliseByBinWidth)
            scale_by_bin_width(slice_mc_plus_ext);
        double ymax = std::max(slice_bnb->hist_->GetMaximum(),
                               slice_mc_plus_ext->hist_->GetMaximum()) * 1.07;
        slice_bnb->hist_->GetYaxis()->SetRangeUser(0., ymax);

        slice_bnb->hist_->Draw("e");
        slice_pred_stack->Draw("hist same");

        slice_mc_plus_ext->hist_->SetLineWidth(3);
        slice_mc_plus_ext->hist_->SetFillColor(kGray + 1);
        slice_mc_plus_ext->hist_->SetFillStyle(3353);
        slice_mc_plus_ext->hist_->Draw("same E2");

        slice_bnb->hist_->Draw("same e");
        slice_bnb->hist_->SetTitle("Selected CC1#pi^{#pm}Xp Events");

        std::cout << "DEBUG Point 3" << std::endl;

        // slice_bnb->hist_->GetYaxis()->SetTitle("Counts");

        // TLegend* lg = new TLegend(0.6,0.7,0.9,0.9);
        // lg->AddEntry(slice_bnb->hist_.get(), "BNB", "l");
        // lg->AddEntry(slice_mc_plus_ext->hist_.get(), "MC + EXT", "f");
        // lg->AddEntry(slice_pred_stack->GetStack()->Last(), "Prediction", "f");
        // lg->Draw();

        std::ostringstream oss;
        auto chi2_result = slice_bnb->get_chi2(*slice_mc_plus_ext);
        oss << std::setprecision(3) << chi2_result.chi2_ << " / "
            << chi2_result.num_bins_ << " bin";
        if (chi2_result.num_bins_ > 1)
            oss << "s; p-value = " << chi2_result.p_value_;
        std::string label = ": #chi^{2} = " + oss.str();

        // create a TLatex object to add the label to the plot
        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextFont(42);
        latex->SetTextSize(0.04);

        // add the label to the plot
        latex->DrawLatex(0.12, 0.86, label.c_str());

        // // Prepare the plot legend
        // TLegend* lg = new TLegend( 0.64, 0.32, 0.94, 0.85 );

        // double pot_on = pot_map.at( NFT::kOnBNB );
        // std::string legend_title = get_legend_title( pot_on );
        // lg->SetHeader( "Legend", "C" );

        // lg->AddEntry( slice_bnb->hist_, "NuWro as data", "lp" );
        // lg->AddEntry( slice_mc_plus_ext->hist_, "MC", "f" );

        std::string out_pdf_name = "plots/plot_slice_";
        if (sl_idx < 10)
            out_pdf_name += "0";
        out_pdf_name += std::to_string(sl_idx);
        out_pdf_name += normaliseByBinWidth ? "_norm.pdf" : ".pdf";
        c1->SaveAs(out_pdf_name.c_str());

        // Get the binning and axis labels for the current slice by cloning the
        // (empty) histogram owned by the Slice object
        TH1 *slice_hist = dynamic_cast<TH1 *>(
            slice.hist_->Clone("slice_hist"));

        slice_hist->SetDirectory(nullptr);

        // Keys are labels, values are fractional uncertainty histograms
        auto *fr_unc_hists = new std::map<std::string, TH1 *>();
        auto &frac_uncertainty_hists = *fr_unc_hists;

        std::cout << "DEBUG Point 4" << std::endl;

        // Show fractional uncertainties computed using these covariance matrices
        // in the ROOT plot. All configured fractional uncertainties will be
        // included in the output pgfplots file regardless of whether they appear
        // in this vector.

        std::vector<std::string> cov_mat_keys = {
            "total", "detVar_total", "flux", "reint", "xsec_total", "POT", "numTargets",
            "MCstats", "EXTstats", "BNBstats"};

        #ifdef USE_FAKE_DATA
                cov_mat_keys = {"total", "xsec_total", "MCstats", "EXTstats"};
        #endif

        // Loop over the various systematic uncertainties
        int color = 0;
        std::cout << "DEBUG Point 5" << std::endl;
        for (const auto &pair : matrix_map)
        {

            const auto &key = pair.first;
            const auto &cov_matrix = pair.second;

            SliceHistogram *slice_for_syst = SliceHistogram::make_slice_histogram(
                *reco_mc_plus_ext_hist, slice, &cov_matrix);

            // The SliceHistogram object already set the bin errors appropriately
            // based on the slice covariance matrix. Just change the bin contents
            // for the current histogram to be fractional uncertainties. Also set
            // the "uncertainties on the uncertainties" to zero.
            // TODO: revisit this last bit, possibly assign bin errors here
            for (const auto &bin_pair : slice.bin_map_)
            {
                int global_bin_idx = bin_pair.first;
                double y = slice_for_syst->hist_->GetBinContent(global_bin_idx);
                double err = slice_for_syst->hist_->GetBinError(global_bin_idx);
                double frac = 0.;
                if (y > 0.)
                    frac = err / y;
                slice_for_syst->hist_->SetBinContent(global_bin_idx, frac);
                slice_for_syst->hist_->SetBinError(global_bin_idx, 0.);
            }

            // Check whether the current covariance matrix name is present in
            // the vector defined above this loop. If it isn't, don't bother to
            // plot it, and just move on to the next one.
            auto cbegin = cov_mat_keys.cbegin();
            auto cend = cov_mat_keys.cend();
            auto iter = std::find(cbegin, cend, key);
            if (iter == cend)
                continue;

            frac_uncertainty_hists[key] = slice_for_syst->hist_.get();

            if (color <= 9)
                ++color;
            if (color == 5)
                ++color;
            if (color >= 10)
                color += 10;
            if (color == 6)
                ++color; // Since we aare skipping total detvar

            slice_for_syst->hist_->SetLineColor(color);
            slice_for_syst->hist_->SetLineWidth(3);
        }
        std::cout << "DEBUG Point 6" << std::endl;

        TCanvas *c2 = new TCanvas;
        TLegend *lg2 = new TLegend(0.7, 0.7, 0.9, 0.9);

        auto *total_frac_err_hist = frac_uncertainty_hists.at("total");
        total_frac_err_hist->SetStats(false);
        total_frac_err_hist->GetYaxis()->SetRangeUser(0.,
                                                      total_frac_err_hist->GetMaximum() * 1.05);
        total_frac_err_hist->SetLineColor(kBlack);
        total_frac_err_hist->SetLineStyle(9);
        total_frac_err_hist->SetLineWidth(3);
        total_frac_err_hist->Draw("hist");
        total_frac_err_hist->SetTitle("Fractional Uncertainty Selected CC1#pi^{#pm}Xp Events");

        // const auto frac_ymax = 0.35;
        // total_frac_err_hist->GetYaxis()->SetRangeUser( 0., frac_ymax);

        lg2->AddEntry(total_frac_err_hist, "total", "l");
        std::cout << "DEBUG Point 7" << std::endl;

        for (auto &pair : frac_uncertainty_hists)
        {
            const auto &name = pair.first;
            TH1 *hist = pair.second;
            // We already plotted the "total" one above
            if (name == "total")
                continue;
            if (name.size() >= 5 && name.substr(name.size() - 5) == "stats")
            {
                hist->SetLineStyle(2);
            }

            lg2->AddEntry(hist, name.c_str(), "l");
            hist->Draw("same hist");

            std::cout << name << " frac err in bin #1 = "
                      << hist->GetBinContent(1) * 100. << "%\n";
        }

        lg2->Draw("same");

        std::string frac_out_pdf_name = "plots/plot_frac_slice_";
        if (sl_idx < 10)
            frac_out_pdf_name += "0";
        frac_out_pdf_name += std::to_string(sl_idx) + ".pdf";
        c2->SaveAs(frac_out_pdf_name.c_str());

        std::cout << "Total frac error in bin #1 = "
                  << total_frac_err_hist->GetBinContent(1) * 100. << "%\n";
        std::cout << "DEBUG Point 8" << std::endl;

    } // slices
    std::cout << "DEBUG Point 9" << std::endl;
}

int tutorial_slice_plots_closure()
{
    slice_plots(true);
    return 0;
}
