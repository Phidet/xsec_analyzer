// Standard library includes
#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <iostream>
#include <array>
#include <set>

// ROOT includes
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TPad.h"

// STV analysis includes
#include "EventCategory.hh"
#include "FiducialVolume.hh"
#include "FilePropertiesManager.hh"
#include "HistUtils.hh"
#include "PlotUtils.hh"

// Abbreviation to make using the enum class easier
using NFT = NtupleFileType;

#define USE_FAKE_DATA "yes"
#define FILE_PROPERTIES "nuwro_file_properties.txt"

void make_plots(const std::string& branchexpr, const std::string& selection,
                const std::string& signal, const std::set<int>& runs,
                std::vector<double> bin_low_edges,
                const std::string& x_axis_label = "",
                const std::string& y_axis_label = "",
                const std::string& title = "",
                const bool effPurToggle = true, // true = efficiency, false = purity
                const std::string& mc_event_weight = DEFAULT_MC_EVENT_WEIGHT)
{
    // Get the number of bins to use in histograms
    int Nbins = bin_low_edges.size() - 1;

    FilePropertiesManager& fpm = FilePropertiesManager::Instance();

    #ifdef USE_FAKE_DATA
      fpm.load_file_properties( FILE_PROPERTIES );
    #endif

    // Define the file type (NFT::kNumuMC)
    const NFT file_type = NFT::kNumuMC;

    // Prepare TChains needed to loop over the event ntuples to be analyzed. Also
    // prepare maps to keep track of the corresponding POT normalizations and
    // total number of triggers (the latter of these is actually used only for
    // data samples).
    std::map< NFT, std::unique_ptr<TChain> > tchain_map;
    tchain_map.emplace( std::make_pair(file_type, new TChain("stv_tree")) );

    // Prepare strings used by multiple histograms below
    std::string hist_name_prefix = "PurityEfficiencyPlots_";
    std::string plot_title = title + "; " + x_axis_label + (effPurToggle ? "; Efficiency" : "; Purity");

    const auto& ntuple_map = fpm.ntuple_file_map();

    // Loop over the runs
    for (const int& run : runs) {
        const auto& run_map = ntuple_map.at(run);
        const auto& ntuple_files = run_map.at(file_type);
        auto* tchain = tchain_map.at(file_type).get();  // Access the TChain for the current file type

        for (const auto& file_name : ntuple_files) {
            // Add the current file to the TChain
            tchain->Add(file_name.c_str());
        } // file names
    }

    // Initialize histograms for selected and signal events
    TH1D* selected_signal_hist = new TH1D("selected_signal_hist", plot_title.c_str(), Nbins, bin_low_edges.data());
    TH1D* denominator_hist = new TH1D("denominator_hist", plot_title.c_str(), Nbins, bin_low_edges.data());

    // Loop over the files and fill the histograms based on the selection and signal branches
    for (const int& run : runs) {
        const auto& ntuple_files = ntuple_map.at(run).at(file_type);
        for (const auto& file_name : ntuple_files) {
            TChain* mc_ch = tchain_map.at(file_type).get();
            mc_ch->Draw((branchexpr + " >> selected_signal_hist").c_str(), (mc_event_weight + "*(" + selection + " && " + signal + ")").c_str(), "goff");
            if(effPurToggle) mc_ch->Draw((branchexpr + " >> denominator_hist").c_str(), (mc_event_weight + "*(" + signal + ")").c_str(), "goff");
            else mc_ch->Draw((branchexpr + " >> denominator_hist").c_str(), (mc_event_weight + "*(" + selection + ")").c_str(), "goff");
        }
    }

    // Compute the efficiency as the ratio of selected events to signal events
    TH1D* eff_hist = new TH1D("eff_hist", plot_title.c_str(), Nbins, bin_low_edges.data());
    eff_hist->Divide(selected_signal_hist, denominator_hist, 1.0, 1.0, "B");

    // Set the style for the efficiency plot
    eff_hist->SetLineColor(kBlack);
    // eff_hist->SetMarkerStyle(20);
    // eff_hist->SetMarkerSize(0.7);
    // eff_hist->SetMarkerColor(kBlack);
    eff_hist->SetFillStyle(3004);
    eff_hist->SetFillColor(kGray);
    eff_hist->GetXaxis()->SetTitle(x_axis_label.c_str());
    eff_hist->GetYaxis()->SetTitle(y_axis_label.c_str());
    eff_hist->SetStats( false );

    // Draw the efficiency plot
    auto* c1 = new TCanvas;
    eff_hist->Draw("E2");

    TH1D* eff_hist_copy = (TH1D*)eff_hist->Clone(); // Create a copy of the histogram
    eff_hist_copy->SetFillStyle(0); // Set fill style to 0 for the copy
    eff_hist_copy->Draw("same hist");

    // Save the plot
    std::string plot_name = ("plots/" + hist_name_prefix + plot_title + ".pdf");
    // Remove ; from plot name
    plot_name.erase(std::remove(plot_name.begin(), plot_name.end(), ';'), plot_name.end());
    // Replace spaces with _
    std::replace(plot_name.begin(), plot_name.end(), ' ', '_');

    c1->SaveAs(plot_name.c_str());

    delete c1;
    delete eff_hist;
    delete eff_hist_copy;
    delete selected_signal_hist;
    delete denominator_hist;
}

// Overloaded version with constant-width binning
void make_plots( const std::string& branchexpr, const std::string& selection, const std::string& signal,
    const std::set<int>& runs, double xmin, double xmax, int Nbins,
    const std::string& x_axis_label = "", const std::string& y_axis_label = "",
    const std::string& title = "",
    const bool effPurToggle = true, // true = efficiency, false = purity
    const std::string& mc_event_weight = DEFAULT_MC_EVENT_WEIGHT )
{
    // Generates a vector of bin low edges equivalent to the approach used by the
    // TH1 constructor that takes xmin and xmax in addition to the number of bins
    auto bin_low_edges = get_bin_low_edges( xmin, xmax, Nbins );

    make_plots( branchexpr, selection, signal, runs, bin_low_edges, x_axis_label,
      y_axis_label, title, effPurToggle, mc_event_weight );
}


void PurityEfficiencyPlots() 
{
    const auto runSet = std::set<int>{1};
    const std::map<std::string, std::vector<double>> binEdges{
        {"muonCosTheta", {-1, -0.27, 0.29, 0.46, 0.58, 0.67, 0.77, 0.82, 0.88, 0.93, 0.97, 1}},
        {"muonPhi", {-3.141592654, -2.722713633, -2.303834613, -1.884955592, -1.466076572, -1.047197551, -0.6283185307, -0.2094395102, 0.2094395102, 0.6283185307, 1.047197551, 1.466076572, 1.884955592, 2.303834613, 2.722713633 , 3.141592654}},
        {"muonMomentum", {0.15, 0.23, 0.32, 0.45, 0.66, 1.5}},
        {"pionCosTheta", {-1, -0.47, 0, 0.39, 0.65, 0.84, 0.93, 1}},
        {"pionPhi", {-3.141592654, -2.513274123, -1.884955592, -1.256637061, -0.6283185307, 0, 0.6283185307, 1.256637061, 1.884955592, 2.513274123, 3.141592654}},
        {"pionMomentum", {0, 0.1, 0.16, 0.19, 0.22, 0.6}},
        {"muonPionAngle", {0, 0.49, 0.93, 1.26, 1.57, 1.88, 2.21, 2.65}},
        {"nProtons", {0, 1, 2, 3}}, // 3 should be labeled std::numeric_limits<double>::max()
        {"total", {-1, 1}},
    };

    // make_plots( /* branchexpr = */        "cc1pi_truth_muonCosTheta",
    //             /* selection = */         "cc1pi_selected_generic",
    //             /* signal = */            "cc1pi_signal",
    //             /* runs = */              runSet,
    //             /* xmin = */              binEdges.at("muonPhi").front(),
    //             /* xmax = */              binEdges.at("muonPhi").back(),
    //             /* Nbins = */             1,
    //             /* x_axis_label = */      "True muonCosTheta",
    //             /* y_axis_label = */      "# Events", 
    //             /* title = */             "CC1pi Selection Efficiency",
    //             /* effPurToggle = */      true); // true = efficiency, false = purity

    // make_plots( /* branchexpr = */        "cc1pi_truth_muonCosTheta",
    //             /* selection = */         "cc1pi_selected_generic",
    //             /* signal = */            "cc1pi_signal",
    //             /* runs = */              runSet,
    //             /* xmin = */              binEdges.at("muonPhi").front(),
    //             /* xmax = */              binEdges.at("muonPhi").back(),
    //             /* Nbins = */             1,
    //             /* x_axis_label = */      "True muonCosTheta",
    //             /* y_axis_label = */      "# Events", 
    //             /* title = */             "CC1pi Selection Purity",
    //             /* effPurToggle = */      false); // true = efficiency, false = purity

    make_plots( /* branchexpr = */        "true_cc1pi",
                /* selection = */         "cc1pi_selected_generic",
                /* signal = */            "true_cc1pi",
                /* runs = */              runSet,
                /* xmin = */              1,
                /* xmax = */              2,
                /* Nbins = */             1,
                /* x_axis_label = */      "True muonCosTheta",
                /* y_axis_label = */      "# Events", 
                /* title = */             "CC1pi Selection Efficiency",
                /* effPurToggle = */      true); // true = efficiency, false = purity

    make_plots( /* branchexpr = */        "cc1pi_selected_generic",
                /* selection = */         "cc1pi_selected_generic",
                /* signal = */            "true_cc1pi",
                /* runs = */              runSet,
                /* xmin = */              1,
                /* xmax = */              2,
                /* Nbins = */             1,
                /* x_axis_label = */      "True muonCosTheta",
                /* y_axis_label = */      "# Events", 
                /* title = */             "CC1pi Selection Purity",
                /* effPurToggle = */      false); // true = efficiency, false = purity

}