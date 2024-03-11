// Standard library includes
#include <iomanip>
#include <map>
#include <memory>
#include <string>

// ROOT includes
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
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
#define FILE_PROPERTIES "file_properties.txt"


void make_plots( const std::string& out_name,
    const std::string& branchexpr, const std::string& selection,
    const std::set<int>& runs, std::vector<double> bin_low_edges,
    const std::string& x_axis_label = "",
    const std::string& y_axis_label = "", const std::string& title = "",
    const std::set<float> &vertical_lines = std::set<float>{},
    const bool xAxisLog = false,
    const bool yAxisLog = false,
    const std::string& mc_event_weight = DEFAULT_MC_EVENT_WEIGHT )
{
    // Get the number of bins to use in histograms
    int Nbins = bin_low_edges.size() - 1;

    // Make a counter that iterates each time this function is called. We'll use
    // it to avoid duplicate histogram names (which can confuse ROOT).
    static long plot_counter = -1;
    ++plot_counter;

    // Get access to the singleton utility classes that we'll need
    const EventCategoryInterpreter& eci = EventCategoryInterpreter::Instance();
    FilePropertiesManager& fpm = FilePropertiesManager::Instance();

    #ifdef USE_FAKE_DATA
      fpm.load_file_properties( FILE_PROPERTIES );
    #endif

    // Consider samples for data taken with the beam on, data taken with the beam
    // off, and CV MC samples for numus, intrinsic nues, and dirt events
    constexpr std::array< NFT, 4 > file_types = { NFT::kOnBNB, NFT::kExtBNB,
      NFT::kNumuMC, NFT::kDirtMC }; // NFT::kIntrinsicNueMC

    // Similar array that includes only the CV MC samples
    constexpr std::array< NFT, 2 > mc_file_types = { NFT::kNumuMC,
      NFT::kDirtMC }; // NFT::kIntrinsicNueMC

    // Prepare TChains needed to loop over the event ntuples to be analyzed. Also
    // prepare maps to keep track of the corresponding POT normalizations and
    // total number of triggers (the latter of these is actually used only for
    // data samples).
    std::map< NFT, std::unique_ptr<TChain> > tchain_map;
    std::map< NFT, double > pot_map;
    std::map< NFT, long > trigger_map;
    for ( const auto& type : file_types ) {
      tchain_map.emplace( std::make_pair(type, new TChain("stv_tree")) );
      pot_map[ type ] = 0.;
      trigger_map[ type ] = 0;
    }

    // Add files for each of the selected runs to the appropriate TChain. Also
    // update the corresponding POT normalizations. Use the FilePropertiesManager
    // to find the right ntuple files for each run.
    const auto& ntuple_map = fpm.ntuple_file_map();
    const auto& data_norm_map = fpm.data_norm_map();

    for ( const int& run : runs ) {
      // Get the map storing the ntuple file names for the current run
      const auto& run_map = ntuple_map.at( run );

      for ( const auto& type : file_types ) {

        // Get the set of ntuple files for the current run and sample type
        const auto& ntuple_files = run_map.at( type );

        // Get access to the corresponding TChain, total POT value, and total
        // number of triggers that we want to use to handle these files
        auto* tchain = tchain_map.at( type ).get();
        double& total_pot = pot_map.at( type );
        long& total_triggers = trigger_map.at( type );

        for ( const auto& file_name : ntuple_files ) {
          // Add the current file to the appropriate TChain
          tchain->Add( file_name.c_str() );

          // For data samples, get normalization information from the
          // FilePropertiesManager and add it to the total (it's not stored in
          // the files themselves)
          if ( type == NFT::kOnBNB || type == NFT::kExtBNB ) {
            const auto& norm_info = data_norm_map.at( file_name );
            total_triggers += norm_info.trigger_count_;
            // This will just be zero for beam-off data. We will calculate an
            // effective value using the trigger counts below.
            total_pot += norm_info.pot_;
          }
          // For MC samples, extract the POT normalization from the TParameter
          // stored in the file
          else if ( type == NFT::kNumuMC || type == NFT::kIntrinsicNueMC
            || type == NFT::kDirtMC )
          {
            TFile temp_file( file_name.c_str(), "read" );
            TParameter<float>* temp_pot = nullptr;
            temp_file.GetObject( "summed_pot", temp_pot );
            double pot = temp_pot->GetVal();
            total_pot += pot;
          }
        } // file names
      } // ntuple types
    } // runs

    // Prepare strings used by multiple histograms below
    std::string hist_name_prefix = "hist_plot" + std::to_string( plot_counter );
    std::string plot_title = title + "; " + x_axis_label + "; " + y_axis_label;

    // Fill the beam-off data histogram using the matching TChain
    std::string off_data_hist_name = hist_name_prefix + "_ext";
    TH1D* off_data_hist = new TH1D( off_data_hist_name.c_str(),
      plot_title.c_str(), Nbins, bin_low_edges.data() );

    TChain* off_chain = tchain_map.at( NFT::kExtBNB ).get();
    off_chain->Draw( (branchexpr + " >> " + off_data_hist_name).c_str(),
      selection.c_str(), "goff" );
    //off_data_hist->SetDirectory( nullptr );

    // We need to scale the beam-off data to an effective POT based on the ratio
    // of the total trigger counts for beam-off and beam-on data. Do that here.
    double pot_on = pot_map.at( NFT::kOnBNB );
    double trigs_on = trigger_map.at( NFT::kOnBNB );
    double trigs_off = trigger_map.at( NFT::kExtBNB );

    // Compute the effective POT and store it in the map
    double ext_effective_pot = trigs_off * pot_on / trigs_on;
    pot_map[ NFT::kExtBNB ] = ext_effective_pot;

    // Scale the beam-off data based on the effective POT
    off_data_hist->Scale( pot_on / ext_effective_pot );

    eci.set_ext_histogram_style( off_data_hist );

    // Fill the beam-on data histogram using the matching TChain
    std::string on_data_hist_name = hist_name_prefix + "_on";
    TH1D* on_data_hist = new TH1D( on_data_hist_name.c_str(),
      plot_title.c_str(), Nbins, bin_low_edges.data() );

    TChain* on_chain = tchain_map.at( NFT::kOnBNB ).get();
    on_chain->Draw( (branchexpr + " >> " + on_data_hist_name).c_str(),
      selection.c_str(), "goff" );

    //on_data_hist->SetDirectory( nullptr );
    on_data_hist->Scale( 1. );

    eci.set_bnb_data_histogram_style( on_data_hist );

    // Initialize empty stacked histograms organized by MC event category
    std::map< EventCategory, TH1D* > mc_hists;
    // Loop over all MC event categories
    for ( const auto& pair : eci.label_map() ) {
      EventCategory cat = pair.first;
      std::string cat_label = pair.second;

      std::string temp_mc_hist_name = hist_name_prefix + "_mc"
        + std::to_string( cat );

      TH1D* temp_mc_hist = new TH1D( temp_mc_hist_name.c_str(),
        plot_title.c_str(), Nbins, bin_low_edges.data() );

      mc_hists[ cat ] = temp_mc_hist;

      //temp_mc_hist->SetDirectory( nullptr );
      eci.set_mc_histogram_style( cat, temp_mc_hist );
    }

    // Loop over the different MC samples and collect their contributions. We
    // have to handle them separately in order to get the POT normalization
    // correct.

    // Counter that avoids duplicate temporary MC histogram names. This is used
    // to avoid annoying ROOT warnings.
    static int dummy_counter = 0;
    for ( const auto& type : mc_file_types ) {

      // std::cout<<"DEBUG plot.C Point 16"<<std::endl;
      TChain* mc_ch = tchain_map.at( type ).get();
      double on_pot = pot_map.at( NFT::kOnBNB );
      double mc_pot = pot_map.at( type );
      // std::cout<<"DEBUG plot.C Point 17"<<std::endl;

      // Add this sample's contribution to the stacked histograms by MC event
      // category
      for ( const auto& pair : eci.label_map() ) {

        EventCategory ec = pair.first;
        // std::cout<<"DEBUG plot.C Point 17.1"<<std::endl;

        std::string temp_mc_hist_name = hist_name_prefix + "_temp_mc"
          + std::to_string( ec ) + "_number" + std::to_string( dummy_counter );

        ++dummy_counter;

        // std::cout<<"DEBUG plot.C Point 17.2"<<std::endl;

        TH1D* temp_mc_hist = new TH1D( temp_mc_hist_name.c_str(),
          plot_title.c_str(), Nbins, bin_low_edges.data() );

        mc_ch->Draw( (branchexpr + " >> " + temp_mc_hist_name).c_str(),
          (mc_event_weight + "*(" + selection + " && category == "
          + std::to_string(ec) + ')').c_str(), "goff" );

        // std::cout<<"DEBUG plot.C Point 17.3"<<std::endl;

        // Scale to the same exposure as the beam-on data
        temp_mc_hist->Scale( on_pot / mc_pot );

        // Add this histogram's contribution (now properly scaled) to the total
        // std::cout<<"DEBUG plot.C Point 18"<<std::endl;
        mc_hists.at( ec )->Add( temp_mc_hist );
        // std::cout<<"DEBUG plot.C Point 19"<<std::endl;

        // We don't need the temporary histogram anymore, so just get rid of it
        delete temp_mc_hist;

      } // event categories

    } // MC samples

    // All the input histograms are now ready. Prepare the plot.
    auto* c1 = new TCanvas("c1", "c1", 1200, 800); // Width is 1200 pixels, height is 800 pixels
    //c1->SetLeftMargin( 0.12 );
    //c1->SetBottomMargin( 1.49 );

    TPad* pad1 = new TPad( "pad1", "", 0.0, 0.23, 0.75, 1.0 );
    pad1->SetBottomMargin( 0.01 );
    pad1->SetRightMargin( 0.03 );
    pad1->SetLeftMargin( 0.13 );
    pad1->SetGridx();
    // Add bottom padding to the pad
    pad1->SetBottomMargin( 0.05 );
    pad1->Draw();
    pad1->cd();
    
    pad1->SetLogx(xAxisLog);
    pad1->SetLogy(yAxisLog);

    if (yAxisLog) on_data_hist->SetMinimum(1);

    on_data_hist->Draw( "E1" );

    // Stack of categorized MC predictions plus extBNB contribution
    THStack* stacked_hist = new THStack( "mc", "" );

    // Sum all contributions into this TH1D so that we can get the overall
    // statistical uncertainty easily
    TH1D* stat_err_hist = new TH1D(
      ("stat_err_hist_" + hist_name_prefix).c_str(), "",
      Nbins, bin_low_edges.data()
    );

    stacked_hist->Add( off_data_hist );
    stat_err_hist->Add( off_data_hist );

    for ( auto citer = mc_hists.crbegin(); citer != mc_hists.crend();
      ++citer )
    {
      TH1D* hist = citer->second;
      stacked_hist->Add( hist );
      stat_err_hist->Add( hist );
    }

    if (yAxisLog) stacked_hist->SetMinimum(1);
    stacked_hist->Draw( "hist same" );

    on_data_hist->Draw( "E1 same" );

    eci.set_stat_err_histogram_style( stat_err_hist );
    if (yAxisLog) stat_err_hist->SetMinimum(1);
    stat_err_hist->Draw( "E2 same" );

    // Adjust y-axis range for stacked plot. Check both the data and the
    // stacked histograms (via the combined stat_err_hist)
    double ymax = stat_err_hist->GetBinContent( stat_err_hist->GetMaximumBin() );
    double ymax2 = on_data_hist->GetBinContent( on_data_hist->GetMaximumBin() );
    if ( ymax < ymax2 ) ymax = ymax2;

    // Redraw the histograms with the updated y-axis range
    on_data_hist->GetYaxis()->SetRangeUser( 1, 1.2*ymax );
    on_data_hist->Draw( "E1 same" );

    for (const auto& x : vertical_lines) {
      TLine* line = new TLine( x, 0, x, 1.2*ymax );
      line->SetLineColor( kBlack );
      line->SetLineStyle( 9 ); // dashed
      line->SetLineWidth( 2 );
      line->Draw();
    }

    // Prepare the plot legend
    c1->cd(); // Go back from pad1 to main canvas c1
    TLegend* lg = new TLegend(0.75, 0.1, 1, 0.95); // Position the legend to the right of the pads

    std::string legend_title = get_legend_title( pot_on );
    lg->SetHeader( legend_title.c_str(), "C" );

    lg->AddEntry( on_data_hist, "Data (beam on)", "lp" );
    lg->AddEntry( stat_err_hist, "Statistical uncertainty", "f" );

    double total_events = stat_err_hist->Integral();
    for ( const auto& pair : eci.label_map() ) {
      EventCategory ec = pair.first;
      std::string label = pair.second;

      // Skip categories with no events
      if ( mc_hists.at( ec )->Integral() == 0 ) continue;

      // std::cout<<"DEBUG plot.C Point 20"<<std::endl;
      TH1* category_hist = mc_hists.at( ec );
      // std::cout<<"DEBUG plot.C Point 21"<<std::endl;

      // Use TH1::Integral() to account for CV reweighting correctly
      double events_in_category = category_hist->Integral();
      double category_percentage = events_in_category * 100. / total_events;

      std::string cat_pct_label = Form( "%.2f%#%", category_percentage );

      lg->AddEntry( category_hist, (label + ", " + cat_pct_label).c_str(), "f" );
    }

    double beam_off_events = off_data_hist->Integral();
    double beam_off_percentage = beam_off_events * 100. / total_events;

    std::string off_pct_label = Form( "%.2f%#%", beam_off_percentage );

    lg->AddEntry( off_data_hist, ("Data (beam off), "
      + off_pct_label).c_str(), "f" );

    lg->SetBorderSize( 0 );

    // Increase the font size for the legend header
    // (see https://root-forum.cern.ch/t/tlegend-headers-font-size/14434)
    TLegendEntry* lg_header = dynamic_cast< TLegendEntry* >(
      lg->GetListOfPrimitives()->First() );
    lg_header->SetTextSize( 0.03 );

    lg->Draw( "same" );

    // Ratio plot
    TPad* pad2 = new TPad( "pad2", "", 0, 0.01, 0.75, 0.23 );
    pad2->SetFrameFillStyle( 4000 );
    pad2->SetBottomMargin( 0.38 );
    pad2->SetTopMargin( 0.01 );
    pad2->SetRightMargin( 0.03 );
    pad2->SetLeftMargin( 0.13 );

    pad2->SetGridx();
    pad2->Draw();
    pad2->cd(); // change current pad to pad2
    pad2->SetLogx(xAxisLog);

    // Ratio plot
    TH1D* h_ratio = dynamic_cast<TH1D*>( on_data_hist->Clone("h_ratio") );
    h_ratio->SetStats( false );
    h_ratio->Divide( stat_err_hist );
    h_ratio->SetLineWidth( 2 );
    h_ratio->SetLineColor( kBlack );
    h_ratio->SetMarkerStyle( kFullCircle );
    h_ratio->SetMarkerSize( 0.8 );
    h_ratio->SetTitle( "" );

    // x-axis
    h_ratio->GetXaxis()->SetTitle( on_data_hist->GetXaxis()->GetTitle() );
    h_ratio->GetXaxis()->CenterTitle( true );
    h_ratio->GetXaxis()->SetLabelSize( 0.12 );
    h_ratio->GetXaxis()->SetTitleSize( 0.18 );
    h_ratio->GetXaxis()->SetTickLength( 0.05 );
    h_ratio->GetXaxis()->SetTitleOffset( 0.9 );

    // y-axis
    h_ratio->GetYaxis()->SetTitle( "Ratio" ); //"#frac{Beam ON}{Beam OFF + MC}" );
    h_ratio->GetYaxis()->CenterTitle( true );
    h_ratio->GetYaxis()->SetLabelSize( 0.08);
    h_ratio->GetYaxis()->SetTitleSize( 0.15 );
    h_ratio->GetYaxis()->SetTitleOffset( 0.35 );

    h_ratio->Draw( "E1" );

    gStyle->SetGridColor( 17 );

    // Adjust y-axis
    double ratio_max = -std::numeric_limits<double>::max();
    double ratio_min = std::numeric_limits<double>::max();

    int nBins = h_ratio->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
      double binContent = h_ratio->GetBinContent(i);
      if (binContent > 0) {
        if (binContent > ratio_max) {
          ratio_max = binContent;
        }
        if (binContent < ratio_min) {
          ratio_min = binContent;
        }
      }
    }

    h_ratio->SetMaximum( ratio_max + ratio_max*0.15 );
    h_ratio->SetMinimum( ratio_min - ratio_min*0.2 );

    gPad->Update();

    // Draw a horizontal dashed line at ratio == 1
    TLine* line = new TLine( h_ratio->GetXaxis()->GetXmin(), 1.0,
      h_ratio->GetXaxis()->GetXmax(), 1.0 );
    line->SetLineColor( kBlack );
    line->SetLineStyle( 9 ); // dashed
    line->Draw();

    c1->Update();

    // Save the plot to a file  named after selection and the runs included
    std::string runNames = "";
    for ( const auto& run : runs ) {
      runNames += std::to_string( run );
    }
    std::string plot_name = "plots/" + out_name + "_runs" + runNames;
    c1->SaveAs( (plot_name + ".pdf").c_str() );
    c1->SaveAs( (plot_name + ".png").c_str() );
    c1->SaveAs( (plot_name + ".C").c_str() );
}

// // Overloaded version with constant-width binning
// void make_plots( const std::string& branchexpr, const std::string& selection,
//     const std::set<int>& runs, double xmin, double xmax, int Nbins,
//     const std::string& x_axis_label = "", const std::string& y_axis_label = "",
//     const std::string& title = "",
//     const std::string& mc_event_weight = DEFAULT_MC_EVENT_WEIGHT )
// {
//     // Generates a vector of bin low edges equivalent to the approach used by the
//     // TH1 constructor that takes xmin and xmax in addition to the number of bins
//     auto bin_low_edges = get_bin_low_edges( xmin, xmax, Nbins );

//     make_plots( branchexpr, selection, runs, bin_low_edges, x_axis_label,
//       y_axis_label, title, mc_event_weight );
// }


std::vector<double> get_log_bin_edges( double xMin, double xMax, int nBins )
{
  std::vector<double> binEdges(nBins + 1);
  double logMin = TMath::Log10(xMin);
  double logMax = TMath::Log10(xMax);
  double binWidth = (logMax - logMin) / nBins;
  for (int i = 0; i <= nBins; i++) {
    binEdges[i] = TMath::Power(10, logMin + i * binWidth);
  }
  return binEdges;
}


void CutPlots() {
  const std::vector<std::set<int>> runSets = {
    // std::set<int>{1},
    // std::set<int>{2},
    // std::set<int>{3},
    std::set<int>{4},
    std::set<int>{5},
    std::set<int>{1,2,3,4,5},
  };

  const std::map<std::string, std::vector<double>> binEdges{
      {"nTracks", {1, 2, 3, 4, 5, 6, 7}},
      {"nUncontained", {0, 1, 2, 3, 4}},
      {"nNonProtons", {1, 2, 3, 4, 5, 6}},
      {"pionTruncatedMeandEdx", get_bin_low_edges( 0, 5, 30 )},
      {"pionPhi", get_bin_low_edges( -3.15, 3.15, 30 )},
      {"muonPhi", get_bin_low_edges( -3.15, 3.15, 30 )},
      // {"pionPhi", {0.f, 1.f, 2.f}},
      // {"muonPhi", {0.f, 1.f, 2.f}},
      {"topologicalScore", get_bin_low_edges( 0, 1, 30 )},      
      {"maxVertexDist", get_log_bin_edges( 0.1, 1000, 30 )},
      {"openingAngle", get_bin_low_edges( 0, 3.15, 30 )},
      {"goldenPionBDTScore", get_bin_low_edges(-0.8, 0.6, 30 )},
  };


  for (const auto& runSet : runSets)
  {
    make_plots( /* out_name = */          "cut1",
                /* branchexpr = */        "event_cutValue_nTracks",
                /* selection = */         "passed_topologicalScoreCC",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("nTracks"),
                /* x_axis_label = */      "Number of fitted tracks",
                /* y_axis_label = */      "# Events", 
                /* title = */             "Cut #1: At least two tracks",
                /* vertical_lines = */ std::set<float>{2.f} );

    make_plots( /* out_name = */         "cut2",
                /* branchexpr = */        "event_cutValue_nUncontained",
                /* selection = */         "passed_min2Tracks",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("nUncontained"),
                /* x_axis_label = */      "Number of uncontained particles",
                /* y_axis_label = */      "# Events", 
                /* title = */             "Cut #2: At most one uncontained",
                /* vertical_lines = */ std::set<float>{2.f} );

    make_plots( /* out_name = */          "cut3",
                /* branchexpr = */        "event_cutValue_nNonProtons",
                /* selection = */         "passed_max1Uncontained",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("nNonProtons"),
                /* x_axis_label = */      "Number of non-protons",
                /* y_axis_label = */      "# Events", 
                /* title = */             "Cut #3: Two non-protons",
                /* vertical_lines = */ std::set<float>{2.f, 3.f} );

    make_plots( /* out_name = */          "cut4",
                /* branchexpr = */        "event_cutValue_pionTruncatedMeandEdx",
                /* selection = */         "passed_2NonProtons",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("pionTruncatedMeandEdx"),
                /* x_axis_label = */      "Reconstructed dE/dx / MeV / cm",
                /* y_axis_label = */      "# Events", 
                /* title = */             "Cut #4: Pion candidate dE/dx is valid",
                /* vertical_lines = */ std::set<float>{1.0f} );

    make_plots( /* out_name = */          "cut5a",
                /* branchexpr = */        "event_cutValue_pionPhi",
                /* selection = */         "passed_pionHasValiddEdx",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("pionPhi"),
                /* x_axis_label = */      "Reconstructed #phi_{#pi} / rad",
                /* y_axis_label = */      "# Events", 
                /* title = */             "Before cut #5: With pion-candidates in gaps");

    make_plots( /* out_name = */          "cut5b",
                /* branchexpr = */        "event_cutValue_pionPhi",
                /* selection = */         "passed_pionNotInGap",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("pionPhi"),
                /* x_axis_label = */      "Reconstructed #phi_{#pi} / rad",
                /* y_axis_label = */      "# Events", 
                /* title = */             "After cut #5: Without pion-candidates in gaps");

    make_plots( /* out_name = */          "cut6a",
                /* branchexpr = */        "event_cutValue_muonPhi",
                /* selection = */         "passed_pionNotInGap",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("muonPhi"),
                /* x_axis_label = */      "Reconstructed #phi_{#mu} / rad",
                /* y_axis_label = */      "# Events",
                /* title = */             "Before cut #6: With muon-candidates in gaps");

    make_plots( /* out_name = */          "cut6b",
                /* branchexpr = */        "event_cutValue_muonPhi",
                /* selection = */         "passed_muonNotInGap",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("muonPhi"),
                /* x_axis_label = */      "Reconstructed #phi_{#mu} / rad",
                /* y_axis_label = */      "# Events", 
                /* title = */             "After cut #6: Without muon-candidates in gaps");

    make_plots( /* out_name = */    "cut7",
                /* branchexpr = */      "event_cutValue_topologicalScore",
                /* selection = */         "passed_muonNotInGap",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("topologicalScore"),
                /* x_axis_label = */      "Neutrino-vs-cosmic topological score",
                /* y_axis_label = */      "# Events", 
                /* title = */             "Cut #7: Topological score",
                /* vertical_lines = */ std::set<float>{0.67f},
                /* xAxisLog = */ false,
                /* yAxisLog = */ true );

    make_plots( /* out_name = */        "cut8",
                /* branchexpr = */      "TMath::Sqrt(event_cutValue_maxVertexDist)",
                /* selection = */         "passed_topologicalScore",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("maxVertexDist"),
                /* x_axis_label = */      "Furthest track-start to vertex distance / cm",
                /* y_axis_label = */      "# Events", 
                /* title = */             "Cut #8: All tracks start close to neutrino vertex",
                /* vertical_lines = */ std::set<float>{9.5f},
                /* xAxisLog = */ true );

    make_plots( /* out_name = */        "cut9",
                /* branchexpr = */      "event_cutValue_openingAngle",
                /* selection = */         "passed_startNearVertex",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("openingAngle"),
                /* x_axis_label = */      "Muon-pion opening angle, #theta_{#mu#pi}",
                /* y_axis_label = */      "# Events", 
                /* title = */             "Cut #9: Cut on opening angle",
                /* vertical_lines = */ std::set<float>{2.65f} );

    make_plots( /* out_name = */          "cut10",
                /* branchexpr = */        "event_cutValue_goldenPionBDT",
                /* selection = */         "passed_openingAngle",
                /* runs = */              runSet,
                /* bin_low_edges = */     binEdges.at("goldenPionBDTScore"),
                /* x_axis_label = */      "Golden pion BDT response",
                /* y_axis_label = */      "# Events", 
                /* title = */             "Cut #10: Golden-pion BDT",
                /* vertical_lines = */ std::set<float>{-0.03f} );
  }
}