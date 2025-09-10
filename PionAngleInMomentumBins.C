#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>


void FillHistogram(TTree* tree, TH1D* hist, const std::string& variable, const std::string& condition, const float runWeight)
{

    TH1D* hist_tmp = (TH1D*)hist->Clone("hist");

    hist_tmp->Reset();

    std::string hist_tmp_name = hist_tmp->GetName();
    std::string weight = "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)*" + std::to_string(runWeight);

    // Fill the histograms with weights
    tree->Project(hist_tmp_name.c_str(), variable.c_str(), ("(" + condition + ")*" + weight).c_str());

    hist_tmp->Sumw2();

    // Add the tmp histograms to the main histograms
    hist->Add(hist_tmp);

    delete hist_tmp;
}

void FillHistogram2D(TTree* tree, TH2D* hist, const std::string& variable1, const std::string& variable2, const std::string& condition, const float runWeight)
{

    TH2D* hist_tmp = (TH2D*)hist->Clone("hist2D");
    hist_tmp->Reset();

    std::string hist_tmp_name = hist_tmp->GetName();
    std::string weight = "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)*" + std::to_string(runWeight);

    // Fill the histograms with weights
    tree->Project(hist_tmp_name.c_str(), (variable2+":"+variable1).c_str(), ("(" + condition + ")*" + weight).c_str());

    // hist_tmp->Sumw2();

    // Add the tmp histograms to the main histograms
    hist->Add(hist_tmp);

    delete hist_tmp;
}

template <typename T>
void FillParticleHistogram(TTree* tree, TH1D* hist, const std::string& variable, const std::string& condition, const float runWeight)
{
    std::vector<T> *values = nullptr;
    // tree->SetBranchAddress(variable.c_str(), &values);
    tree->SetBranchAddress("trueCC1pi_recoParticle_muonBDTScore_vector", &values);

    // Check if values is null
    if (!values)
    {
        std::cerr << "Error: No branch named " << variable << " found in the tree." << std::endl;
        return;
    }

    // Loop over the entries in the tree
    for (Long64_t i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);

        std::cout << "DEBUG FillParticleHistogram Point values->size(): " << values->size() << std::endl;
        // Apply the condition
        if (condition.empty() || tree->GetLeaf(condition.c_str())->GetValue() != 0)
        {
            // Loop over the values in the vector and fill the histogram
            for (const auto& value : *values)
            {
                hist->Fill(value, runWeight);
            }
        }
    }
}

void MakePlot(TH1D* hist, const std::string& name)
{

    // Set the color and line thickness of the histograms
    hist->SetLineColor(kGreen);
    hist->SetLineWidth(3);

    // Remove the stats box
    hist->SetStats(0);

    // Draw the histogram
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    hist->Draw("E hist");

    // Add a legend
    // TLegend* legend = new TLegend(0.2, 0.1);
    // legend->AddEntry(hist, "CC0pi", "l");
    // legend->AddEntry(cc1pi_hist, "CC1pi", "l");
    // legend->Draw();

    c1->SaveAs(("plots/PionAngleInMomentumBins_" + name + "_testingOnly_lowPiMomThreshold.pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void MakePlot2D(TH2D* hist, const std::string& name, const bool drawDiagonal = false, const bool axisTicks = true, const bool drawText = false, const std::vector<float>& binEdges = {})
{
    // // Set the color and line thickness of the histograms
    // hist->SetLineColor(kGreen);
    // hist->SetLineWidth(3);

    // Remove the stats box
    hist->SetStats(0);

    // Draw the histogram
    TCanvas* c1 = new TCanvas(name.c_str(), name.c_str(), 800, 600);

    if(!axisTicks)
    {
        hist->GetXaxis()->SetNdivisions(hist->GetNbinsX(), kFALSE);
        hist->GetYaxis()->SetNdivisions(hist->GetNbinsY(), kFALSE);
    }

    if (drawText)
    {
        // Set the text format to 2 digits after the decimal point
        gStyle->SetPaintTextFormat("4.2f");
        // Increase the text size
        gStyle->SetTextSize(0.2);
        // Draw the histogram with text
        hist->Draw("COLZ TEXT");
    }
    else
    {
        hist->Draw("COLZ");
    }

    if (drawDiagonal)
    {
        // Add a diagonal line
        double x1 = hist->GetXaxis()->GetXmin();
        double y1 = hist->GetYaxis()->GetXmin();
        double x2 = hist->GetXaxis()->GetXmax();
        double y2 = hist->GetYaxis()->GetXmax();
        TLine *line = new TLine(x1, y1, x2, y2);
        line->SetLineColor(kBlack);  // Set line color to black
        line->SetLineStyle(2);       // Set line style to dashed
        line->SetLineWidth(2);       // Set line width
        line->Draw("same");          // Draw line on the same canvas
    }

    // Draw vertical and horizontal lines at bin edges when binEdges is not empty
    if (!binEdges.empty()) {
        auto previousEdge = binEdges.at(0);
        for (int i = 0; i < binEdges.size(); i++) {
            const auto edge = binEdges.at(i);
            // const auto nextEdge = binEdges.at(i+1);
            // Draw vertical line
            // TLine* vLine = new TLine(edge, previousEdge, edge, nextEdge);
            TLine* vLine = new TLine(edge, hist->GetYaxis()->GetXmin(), edge, hist->GetYaxis()->GetXmax()); 
            vLine->SetLineColor(kRed);
            vLine->SetLineStyle(2);
            vLine->SetLineWidth(2);
            vLine->Draw("same");

            // Draw horizontal line
            // TLine* hLine = new TLine(previousEdge, edge, nextEdge, edge);
            TLine* hLine = new TLine(hist->GetXaxis()->GetXmin(), edge, hist->GetXaxis()->GetXmax(), edge);
            hLine->SetLineColor(kRed);
            hLine->SetLineStyle(2);
            hLine->SetLineWidth(2);
            hLine->Draw("same");

            previousEdge = edge;
        }
    }

    gStyle->SetPadRightMargin(0.15);

    // Print out hist title
    std::cout << "DEBUG - hist title 0: " << hist->GetTitle() << std::endl;

    const auto title = hist->GetTitle();
    const auto xTitle = hist->GetXaxis()->GetTitle();
    const auto yTitle = hist->GetYaxis()->GetTitle();


    // Add plot title
    hist->SetTitle(title);
    // Add axis labels title
    hist->GetXaxis()->SetTitle(xTitle);
    hist->GetYaxis()->SetTitle(yTitle);

    // Set x axis label offset
    hist->GetXaxis()->SetTitleOffset(1.2);

    std::cout << "DEBUG - hist title 1: " << hist->GetTitle() << std::endl;

    c1->SaveAs(("plots/PionAngleInMomentumBins_" + name + "_testingOnly_lowPiMomThreshold.pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void PionAngleInMomentumBins() 
{
    // // List of files
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        std::make_tuple("nu mc", "1",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
        std::make_tuple("nu mc", "2",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root", 0.25750*2.0),
        std::make_tuple("nu mc", "3",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root", 0.20113*2.0),
        std::make_tuple("nu mc", "4bcd", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root", 0.13074*2.0),
        std::make_tuple("nu mc", "5",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196*2.0),
    };

    // 2d histograms
    std::vector<double> phi_bin_edges = {-3.141592654, -2.513274123, -1.884955592, -1.256637061, -0.6283185307, 0, 0.6283185307, 1.256637061, 1.884955592, 2.513274123, 3.141592654};
    TH2D* h_piPhi_0 = new TH2D("pi_phi_h0", "Truth vs Reco #phi_{#pi} for Selected Signal Events with true p_{#pi} < 0 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_1 = new TH2D("pi_phi_h1", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 100 MeV <= true p_{#pi} < 160 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_2 = new TH2D("pi_phi_h2", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 160 MeV <= true p_{#pi} < 190 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_3 = new TH2D("pi_phi_h3", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 190 MeV <= true p_{#pi} < 220 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_4 = new TH2D("pi_phi_h4", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 220 MeV <= true p_{#pi} < 600 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_5 = new TH2D("pi_phi_h5", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 600 MeV <= true p_{#pi}; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);

    TH2D* h_piPhi_reco_0 = new TH2D("pi_phi_reco_h0", "Truth vs Reco #phi_{#pi} for Selected Signal Events with reco p_{#pi} < 0 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_reco_1 = new TH2D("pi_phi_reco_h1", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 100 MeV <= reco p_{#pi} < 160 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_reco_2 = new TH2D("pi_phi_reco_h2", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 160 MeV <= reco p_{#pi} < 190 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_reco_3 = new TH2D("pi_phi_reco_h3", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 190 MeV <= reco p_{#pi} < 220 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_reco_4 = new TH2D("pi_phi_reco_h4", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 220 MeV <= reco p_{#pi} < 600 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_reco_5 = new TH2D("pi_phi_reco_h5", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 600 MeV <= reco p_{#pi}; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);

    TH2D* h_piPhi_fineBinning_true_00 = new TH2D("h_piPhi_fineBinning_true_00", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 0 MeV <= true p_{#pi} < 80 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_fineBinning_true_80 = new TH2D("h_piPhi_fineBinning_true_80", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 80 MeV <= true p_{#pi} < 100 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_fineBinning_true_100 = new TH2D("h_piPhi_fineBinning_true_100", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 100 MeV <= true p_{#pi} < 120 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_fineBinning_true_120 = new TH2D("h_piPhi_fineBinning_true_120", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 120 MeV <= true p_{#pi} < 140 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);

    
    TH2D* h_piPhi_fineBinning_reco_00 = new TH2D("h_piPhi_fineBinning_reco_00", "Truth vs Reco #phi_{#pi} for Selected Signal Events with reco p_{#pi} < 100 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_fineBinning_reco_100 = new TH2D("h_piPhi_fineBinning_reco_100", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 100 MeV <= reco p_{#pi} < 120 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_fineBinning_reco_120 = new TH2D("h_piPhi_fineBinning_reco_120", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 120 MeV <= reco p_{#pi} < 140 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_fineBinning_reco_140 = new TH2D("h_piPhi_fineBinning_reco_140", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 140 MeV <= reco p_{#pi} < 160 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);
    TH2D* h_piPhi_fineBinning_reco_160 = new TH2D("h_piPhi_fineBinning_reco_160", "Truth vs Reco #phi_{#pi} for Selected Signal Events with 160 MeV <= reco p_{#pi} < 180 MeV; True #phi_{#pi}; Reconstructed #phi_{#pi}", phi_bin_edges.size() - 1, &phi_bin_edges[0], phi_bin_edges.size() - 1, &phi_bin_edges[0]);

    // TH2D* h_piMom_recoGolden = new TH2D("h_piMom_recoGolden", "Truth vs reco p_{#pi} for signal events in unscattered-enhanced selection; True p_{#pi}; Reconstructed p_{#pi}", 100, 0, 0.6, 100, 0, 0.6);
    // TH2D* h_piMom_reco = new TH2D("h_piMom_reco", "Truth vs reco p_{#pi} for signal events in generic selection; True p_{#pi}; Reconstructed p_{#pi}", 100, 0, 0.6, 100, 0, 0.6);

    TH2D* h_piMom_recoGolden = new TH2D("h_piMom_recoGolden", "; True p_{#pi} (GeV); Reconstructed p_{#pi} (GeV)", 100, 0, 0.6, 100, 0, 0.6);
    TH2D* h_piMom_reco = new TH2D("h_piMom_reco", "; True p_{#pi} (GeV); Reconstructed p_{#pi} (GeV)", 100, 0, 0.6, 100, 0, 0.6);

    TH2D* h_piMom_trueGolden = new TH2D("h_piMom_trueGolden", "Truth vs Reco p_{#pi} for Golden Selection Golden Signal Events; True p_{#pi} (GeV); Reconstructed p_{#pi} (GeV)", 100, 0, 0.6, 100, 0, 0.6);

    // std::vector<double> cosTheta_bin_edges = {-1, -0.47, 0, 0.39, 0.65, 0.84, 0.93, 1};

    // Loop over the files
    for (int f = 0; f < files.size(); f++)
    {
        const auto [sampleType, run, filePath, runWeight] = files.at(f);
        std::cout<<"DEBUG - filePath: "<<filePath<<std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(filePath.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");

        FillHistogram2D(tree, h_piPhi_0,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.0 && cc1pi_truth_pionMomentum < 0.1", runWeight);

        FillHistogram2D(tree, h_piPhi_1, 
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi", 
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionMomentum < 0.16", runWeight);

        FillHistogram2D(tree, h_piPhi_2,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.16 && cc1pi_truth_pionMomentum < 0.19", runWeight);

        FillHistogram2D(tree, h_piPhi_3,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.19 && cc1pi_truth_pionMomentum < 0.22", runWeight);

        FillHistogram2D(tree, h_piPhi_4,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.22 && cc1pi_truth_pionMomentum < 0.6", runWeight);

        FillHistogram2D(tree, h_piPhi_5,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.6", runWeight);

        // **********************************
        // Reco pion momentum
        // **********************************

        FillHistogram2D(tree, h_piPhi_reco_0,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.0 && cc1pi_reco_pionMomentum < 0.1", runWeight);

        FillHistogram2D(tree, h_piPhi_reco_1,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_reco_pionMomentum < 0.16", runWeight);

        FillHistogram2D(tree, h_piPhi_reco_2,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.16 && cc1pi_reco_pionMomentum < 0.19", runWeight);

        FillHistogram2D(tree, h_piPhi_reco_3,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.19 && cc1pi_reco_pionMomentum < 0.22", runWeight);

        FillHistogram2D(tree, h_piPhi_reco_4,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.22 && cc1pi_reco_pionMomentum < 0.6", runWeight);

        FillHistogram2D(tree, h_piPhi_reco_5,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.6", runWeight);

        // **********************************
        // True momentum fine binnin around 100 MeV
        // **********************************

        FillHistogram2D(tree, h_piPhi_fineBinning_true_00,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.00 && cc1pi_truth_pionMomentum < 0.08", runWeight);

        FillHistogram2D(tree, h_piPhi_fineBinning_true_80,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.08 && cc1pi_truth_pionMomentum < 0.1", runWeight);

        FillHistogram2D(tree, h_piPhi_fineBinning_true_100,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionMomentum < 0.12", runWeight);

        FillHistogram2D(tree, h_piPhi_fineBinning_true_120,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.12 && cc1pi_truth_pionMomentum < 0.14", runWeight);

        // **********************************
        // Reco momentum fine binnin around 100 MeV
        // **********************************

        FillHistogram2D(tree, h_piPhi_fineBinning_reco_00,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum < 0.1", runWeight);

        FillHistogram2D(tree, h_piPhi_fineBinning_reco_100,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_reco_pionMomentum < 0.12", runWeight);

        FillHistogram2D(tree, h_piPhi_fineBinning_reco_120,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.12 && cc1pi_reco_pionMomentum < 0.14", runWeight);

        FillHistogram2D(tree, h_piPhi_fineBinning_reco_140,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.14 && cc1pi_reco_pionMomentum < 0.16", runWeight);

        FillHistogram2D(tree, h_piPhi_fineBinning_reco_160,
        /*variable1*/ "cc1pi_truth_pionPhi",
        /*variable2*/ "cc1pi_reco_pionPhi",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.16 && cc1pi_reco_pionMomentum < 0.18", runWeight);


        // **********************************
        // Golden Selection
        // **********************************

        FillHistogram2D(tree, h_piMom_recoGolden,
        /*variable1*/ "cc1pi_truth_pionMomentum",
        /*variable2*/ "cc1pi_reco_pionMomentum",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_golden", runWeight);

        FillHistogram2D(tree, h_piMom_reco,
        /*variable1*/ "cc1pi_truth_pionMomentum",
        /*variable2*/ "cc1pi_reco_pionMomentum",
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic", runWeight);

        FillHistogram2D(tree, h_piMom_trueGolden,
        /*variable1*/ "cc1pi_truth_pionMomentum",
        /*variable2*/ "cc1pi_reco_pionMomentum",
        /*Condition*/ "cc1pi_signal && true_golden_cc1pi && cc1pi_selected_golden", runWeight);


        // **********************************
        // CosTheta
        // **********************************

        // FillHistogram2D(tree, h_thetaPhi_0,
        // /*variable1*/ "cc1pi_reco_pionCosTheta",
        // /*variable2*/ "cc1pi_truth_pionCosTheta",
        // /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.0 && cc1pi_truth_pionMomentum < 0.1", runWeight);

        // FillHistogram2D(tree, h_thetaPhi_1,
        // /*variable1*/ "cc1pi_reco_pionCosTheta",
        // /*variable2*/ "cc1pi_truth_pionCosTheta",
        // /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionMomentum < 0.16", runWeight);

        // FillHistogram2D(tree, h_thetaPhi_2,
        // /*variable1*/ "cc1pi_reco_pionCosTheta",
        // /*variable2*/ "cc1pi_truth_pionCosTheta",
        // /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.16 && cc1pi_truth_pionMomentum < 0.19", runWeight);

        // FillHistogram2D(tree, h_thetaPhi_3,
        // /*variable1*/ "cc1pi_reco_pionCosTheta",
        // /*variable2*/ "cc1pi_truth_pionCosTheta",
        // /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.19 && cc1pi_truth_pionMomentum < 0.22", runWeight);

        // FillHistogram2D(tree, h_thetaPhi_4,
        // /*variable1*/ "cc1pi_reco_pionCosTheta",
        // /*variable2*/ "cc1pi_truth_pionCosTheta",
        // /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.22 && cc1pi_truth_pionMomentum < 0.6", runWeight);

        // FillHistogram2D(tree, h_thetaPhi_5,
        // /*variable1*/ "cc1pi_reco_pionCosTheta",
        // /*variable2*/ "cc1pi_truth_pionCosTheta",
        // /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.6", runWeight);


        tFile->Close();
    }
    std::cout<<"DEBUG - Done with file loop"<<std::endl;

    MakePlot2D(h_piPhi_0, "recoVsTruePiPhi_bin0", false, true, false);
    MakePlot2D(h_piPhi_1, "recoVsTruePiPhi_bin1", false, true, false);
    MakePlot2D(h_piPhi_2, "recoVsTruePiPhi_bin2", false, true, false);
    MakePlot2D(h_piPhi_3, "recoVsTruePiPhi_bin3", false, true, false);
    MakePlot2D(h_piPhi_4, "recoVsTruePiPhi_bin4", false, true, false);
    MakePlot2D(h_piPhi_5, "recoVsTruePiPhi_bin5", false, true, false);

    MakePlot2D(h_piPhi_reco_0, "recoVsTruePiPhi_recoBin0", false, true, false);
    MakePlot2D(h_piPhi_reco_1, "recoVsTruePiPhi_recoBin1", false, true, false);
    MakePlot2D(h_piPhi_reco_2, "recoVsTruePiPhi_recoBin2", false, true, false);
    MakePlot2D(h_piPhi_reco_3, "recoVsTruePiPhi_recoBin3", false, true, false);
    MakePlot2D(h_piPhi_reco_4, "recoVsTruePiPhi_recoBin4", false, true, false);
    MakePlot2D(h_piPhi_reco_5, "recoVsTruePiPhi_recoBin5", false, true, false);

    MakePlot2D(h_piPhi_fineBinning_true_00, "recoVsTruePiPhi_fineBinning_true_00", false, true, false);
    MakePlot2D(h_piPhi_fineBinning_true_80, "recoVsTruePiPhi_fineBinning_true_80", false, true, false);
    MakePlot2D(h_piPhi_fineBinning_true_100, "recoVsTruePiPhi_fineBinning_true_100", false, true, false);
    MakePlot2D(h_piPhi_fineBinning_true_120, "recoVsTruePiPhi_fineBinning_true_120", false, true, false);

    MakePlot2D(h_piPhi_fineBinning_reco_00, "recoVsTruePiPhi_fineBinning_reco_80", false, true, false);
    MakePlot2D(h_piPhi_fineBinning_reco_100, "recoVsTruePiPhi_fineBinning_reco_100", false, true, false);
    MakePlot2D(h_piPhi_fineBinning_reco_120, "recoVsTruePiPhi_fineBinning_reco_120", false, true, false);
    MakePlot2D(h_piPhi_fineBinning_reco_140, "recoVsTruePiPhi_fineBinning_reco_140", false, true, false);
    MakePlot2D(h_piPhi_fineBinning_reco_160, "recoVsTruePiPhi_fineBinning_reco_160", false, true, false);

    // const std::vector<float> binEdges = {0.1, 0.16, 0.19, 0.22}; // No 0.6 for this plot because of merged overflow bin
    const std::vector<float> binEdges = {0.1};
    MakePlot2D(h_piMom_recoGolden, "recoVsTruePiMom_goldenSelection", false, true, false, binEdges);
    MakePlot2D(h_piMom_reco, "recoVsTruePiMom_genericSelection", false, true, false, binEdges);
    MakePlot2D(h_piMom_trueGolden, "recoVsTruePiMom_goldenSelection_trueGolden", false, true, false, binEdges);

    // MakePlot2D(h_thetaPhi_0, "recoVsTruePiCosTheta_bin0", false, true, false);
    // MakePlot2D(h_thetaPhi_1, "recoVsTruePiCosTheta_bin1", false, true, false);
    // MakePlot2D(h_thetaPhi_2, "recoVsTruePiCosTheta_bin2", false, true, false);
    // MakePlot2D(h_thetaPhi_3, "recoVsTruePiCosTheta_bin3", false, true, false);
    // MakePlot2D(h_thetaPhi_4, "recoVsTruePiCosTheta_bin4", false, true, false);
    // MakePlot2D(h_thetaPhi_5, "recoVsTruePiCosTheta_bin5", false, true, false);
}