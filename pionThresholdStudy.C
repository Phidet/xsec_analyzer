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

    c1->SaveAs(("plots/pionThresholdStudy_" + name + "_testingOnly_lowPiMomThreshold.pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void MakePlot2D(TH2D* hist, const std::string& name, const bool drawDiagonal = false, const bool axisTicks = true, const bool drawText = false)
{
    // Set the color and line thickness of the histograms
    hist->SetLineColor(kGreen);
    hist->SetLineWidth(3);

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

    // Add more space on the right side of the plot for the z-axis labels
    gStyle->SetPadRightMargin(0.15);

    c1->SaveAs(("plots/pionThresholdStudy_" + name + "_testingOnly_lowPiMomThreshold.pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void pionThresholdStudy() 
{

    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        std::make_tuple("nu mc", "1",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
        // std::make_tuple("nu mc", "2",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root", 0.25750*2.0),
        // std::make_tuple("nu mc", "3",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root", 0.20113*2.0),
        // std::make_tuple("nu mc", "4bcd", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root", 0.13074*2.0),
        // std::make_tuple("nu mc", "5",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196*2.0),
    };

    // 1d histograms
    TH1D* reco_pion_momentum_selectedCC1pi = new TH1D("reco_pion_momentum_selectedCC1pi", "Selected CC1pi Events; (Reco-True)/True Pion Momentum; Area normalised event count", 30, -1, 1);
    TH1D* reco_pion_momentum_selectedCC1pi_0_100MeV = new TH1D("reco_pion_momentum_selectedCC1pi_0_100MeV", "Selected CC1pi Events - 0 to 100 MeV True Pion Momentum; (Reco-True)/True Pion Momentum; Area normalised event count", 30, -1, 1);
    TH1D* reco_pion_momentum_selectedCC1pi_100_200MeV = new TH1D("reco_pion_momentum_selectedCC1pi_100_200MeV 100", "Selected CC1pi Events - 100 to 200 MeV True Pion Momentum; (Reco-True)/True Pion Momentum; Area normalised event count", 30, -1, 1);
    TH1D* reco_pion_momentum_selectedCC1pi_over_200MeV = new TH1D("reco_pion_momentum_selectedCC1pi_over_200MeV", "Selected CC1pi Events - Over 200 MeV True Pion Momentum; (Reco-True)/True Pion Momentum; Area normalised event count", 30, -1, 1);

    // 1d histograms
    TH1D* reco_pion_momentum_selectedCC1pi_trueGolden = new TH1D("reco_pion_momentum_selectedCC1pi_trueGolden", "Selected CC1pi Events; (Reco-True)/True Pion Momentum; Area normalised event count", 30, -1, 1);
    TH1D* reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden = new TH1D("reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden", "Selected CC1pi Events - 0 to 100 MeV True Pion Momentum; (Reco-True)/True Pion Momentum; Area normalised event count", 30, -1, 1);
    TH1D* reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden = new TH1D("reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden 100", "Selected CC1pi Events - 100 to 200 MeV True Pion Momentum; (Reco-True)/True Pion Momentum; Area normalised event count", 30, -1, 1);
    TH1D* reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden = new TH1D("reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden", "Selected CC1pi Events - Over 200 MeV True Pion Momentum; (Reco-True)/True Pion Momentum; Area normalised event count", 30, -1, 1);


    // 2d histograms
    
    // Show the true vs reco pion momentum
    TH2D* recoPionMomentumVStruePionMomentum_selectedTrueCC1pi = new TH2D("recoPionMomentumVStruePionMomentum_selectedTrueCC1pi", "Selected True CC1pi Events; Reco Pion Momentum (GeV/c); True Pion Momentum (GeV/c)", 30, 0, 0.6, 30, 0, 0.6);
    // The same as above but reduced to 2 bins in each dimension separating in and out of phase space regions 
    const double binEdges[] = {0, 0.1, 1.5};
    TH2D* recoPionMomentumVStruePionMomentum_selectedTrueCC1pi_minBins = new TH2D("recoPionMomentumVStruePionMomentum_selectedTrueCC1pi_minBins", "Selected True CC1pi Events; Reco Pion Momentum (GeV/c); True Pion Momentum (GeV/c)", 2, binEdges, 2, binEdges);

    // Loop over the files
    for (int f = 0; f < files.size(); f++)
    {
        const auto [sampleType, run, filePath, runWeight] = files.at(f);
        std::cout<<"DEBUG - filePath: "<<filePath<<std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(filePath.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");

        std::cout<<"DEBUG - Point X0"<<std::endl;

        // Disable all branches
        tree->SetBranchStatus("*", 0);

        // Enable only the branches you need
        // tree->SetBranchStatus("cc0pi_*", 1);
        tree->SetBranchStatus("true*", 1);
        tree->SetBranchStatus("cc1pi_*", 1);
        // tree->SetBranchStatus("true_cc0pi", 1);
        tree->SetBranchStatus("spline_weight", 1);
        tree->SetBranchStatus("tuned_cv_weight", 1);

        // std::cout<<"DEBUG - Point X0.1"<<std::endl;

        FillHistogram(tree, reco_pion_momentum_selectedCC1pi,
        /*variable*/ "(cc1pi_reco_pionMomentum-cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum",
        /*Condition*/ "true_cc1pi && cc1pi_selected_generic", runWeight);

        FillHistogram(tree, reco_pion_momentum_selectedCC1pi_0_100MeV,
        /*variable*/ "(cc1pi_reco_pionMomentum-cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum",
        /*Condition*/ "true_cc1pi && cc1pi_selected_generic && cc1pi_truth_pionMomentum < 0.1", runWeight);

        FillHistogram(tree, reco_pion_momentum_selectedCC1pi_100_200MeV,
        /*variable*/ "(cc1pi_reco_pionMomentum-cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum",
        /*Condition*/ "true_cc1pi && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionMomentum < 0.2", runWeight);

        FillHistogram(tree, reco_pion_momentum_selectedCC1pi_over_200MeV,
        /*variable*/ "(cc1pi_reco_pionMomentum-cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum",
        /*Condition*/ "true_cc1pi && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.2", runWeight);

        FillHistogram(tree, reco_pion_momentum_selectedCC1pi_trueGolden,
        /*variable*/ "(cc1pi_reco_pionMomentum-cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum",
        /*Condition*/ "true_golden_cc1pi && cc1pi_selected_generic", runWeight);

        FillHistogram(tree, reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden,
        /*variable*/ "(cc1pi_reco_pionMomentum-cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum",
        /*Condition*/ "true_golden_cc1pi && cc1pi_selected_generic && cc1pi_truth_pionMomentum < 0.1", runWeight);

        FillHistogram(tree, reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden,
        /*variable*/ "(cc1pi_reco_pionMomentum-cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum",
        /*Condition*/ "true_golden_cc1pi && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionMomentum < 0.2", runWeight);

        FillHistogram(tree, reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden,
        /*variable*/ "(cc1pi_reco_pionMomentum-cc1pi_truth_pionMomentum)/cc1pi_truth_pionMomentum",
        /*Condition*/ "true_golden_cc1pi && cc1pi_selected_generic && cc1pi_truth_pionMomentum >= 0.2", runWeight);



        // FillHistogram2D(tree, recoPionMomentumVStruePionMomentum_selectedTrueCC1pi,
        // /*variable1*/ "cc1pi_reco_pionMomentum",
        // /*variable2*/ "cc1pi_truth_pionMomentum",
        // /*Condition*/ "true_cc1pi && cc1pi_selected_generic", runWeight);

        // FillHistogram2D(tree, recoPionMomentumVStruePionMomentum_selectedTrueCC1pi_minBins,
        // /*variable1*/ "cc1pi_reco_pionMomentum",
        // /*variable2*/ "cc1pi_truth_pionMomentum",
        // /*Condition*/ "true_cc1pi && cc1pi_selected_generic", runWeight);

        tFile->Close();
    }

    // MakePlot2D(muonBDTVSrange_trueMuon_trueCC1pi, "muonBDTVSrange_trueMuon_trueCC1pi");
    // MakePlot2D(muonBDTVSTrueRange_trueMuon_trueCC1pi, "muonBDTVSTrueRange_trueMuon_trueCC1pi");

    // MakePlot2D(muonBDTVSrange_truePion_trueCC1pi, "muonBDTVSrange_truePion_trueCC1pi");
    // MakePlot2D(muonBDTVSTrueRange_truePion_trueCC1pi, "muonBDTVSTrueRange_truePion_trueCC1pi");

    // MakePlot(histograms.at(names.at(i)), cc1pi_histograms.at(names.at(i)), names.at(i));

    // MakePlot(reco_pion_momentum_selectedCC1pi, "reco_pion_momentum_selectedCC1pi");

    // Set the color and line thickness of the histograms
    reco_pion_momentum_selectedCC1pi->SetLineColor(kBlack);
    reco_pion_momentum_selectedCC1pi->SetLineWidth(3);
    reco_pion_momentum_selectedCC1pi_0_100MeV->SetLineColor(kGreen);
    reco_pion_momentum_selectedCC1pi_0_100MeV->SetLineWidth(3);
    reco_pion_momentum_selectedCC1pi_100_200MeV->SetLineColor(kBlue);
    reco_pion_momentum_selectedCC1pi_100_200MeV->SetLineWidth(3);
    reco_pion_momentum_selectedCC1pi_over_200MeV->SetLineColor(kRed);
    reco_pion_momentum_selectedCC1pi_over_200MeV->SetLineWidth(3);

    // Area normalize all three histograms
    reco_pion_momentum_selectedCC1pi->Scale(1.0 / reco_pion_momentum_selectedCC1pi->Integral());
    reco_pion_momentum_selectedCC1pi_0_100MeV->Scale(1.0 / reco_pion_momentum_selectedCC1pi_0_100MeV->Integral());
    reco_pion_momentum_selectedCC1pi_100_200MeV->Scale(1.0 / reco_pion_momentum_selectedCC1pi_100_200MeV->Integral());
    reco_pion_momentum_selectedCC1pi_over_200MeV->Scale(1.0 / reco_pion_momentum_selectedCC1pi_over_200MeV->Integral());

    // Remove the stats box
    reco_pion_momentum_selectedCC1pi_0_100MeV->SetStats(0);

    // Draw the histogram
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    // Draw the first histogram
    reco_pion_momentum_selectedCC1pi_0_100MeV->Draw("E hist");

    // Find the maximum value across all histograms
    double max_val = reco_pion_momentum_selectedCC1pi_0_100MeV->GetBinContent(reco_pion_momentum_selectedCC1pi_0_100MeV->GetMaximumBin()) + reco_pion_momentum_selectedCC1pi_0_100MeV->GetBinError(reco_pion_momentum_selectedCC1pi_0_100MeV->GetMaximumBin());
    double temp_max = reco_pion_momentum_selectedCC1pi_100_200MeV->GetBinContent(reco_pion_momentum_selectedCC1pi_100_200MeV->GetMaximumBin()) + reco_pion_momentum_selectedCC1pi_100_200MeV->GetBinError(reco_pion_momentum_selectedCC1pi_100_200MeV->GetMaximumBin());
    if (temp_max > max_val) max_val = temp_max;
    temp_max = reco_pion_momentum_selectedCC1pi_over_200MeV->GetBinContent(reco_pion_momentum_selectedCC1pi_over_200MeV->GetMaximumBin()) + reco_pion_momentum_selectedCC1pi_over_200MeV->GetBinError(reco_pion_momentum_selectedCC1pi_over_200MeV->GetMaximumBin());
    if (temp_max > max_val) max_val = temp_max;
    temp_max = reco_pion_momentum_selectedCC1pi->GetBinContent(reco_pion_momentum_selectedCC1pi->GetMaximumBin()) + reco_pion_momentum_selectedCC1pi->GetBinError(reco_pion_momentum_selectedCC1pi->GetMaximumBin());
    if (temp_max > max_val) max_val = temp_max;

    // Set the y-axis range
    reco_pion_momentum_selectedCC1pi_0_100MeV->GetYaxis()->SetRangeUser(0, max_val * 1.1); // Add 10% for a little extra space

    // Draw the other histograms
    reco_pion_momentum_selectedCC1pi_100_200MeV->Draw("E hist same");
    reco_pion_momentum_selectedCC1pi_over_200MeV->Draw("E hist same");
    reco_pion_momentum_selectedCC1pi->Draw("E hist same");

    // Add a legend
    TLegend* legend = new TLegend(0.2, 0.2);
    legend->SetHeader("Selected True CC1pi Events", "C"); // "C" centers the header
    legend->AddEntry(reco_pion_momentum_selectedCC1pi_0_100MeV, "0 to 100 MeV", "l");
    legend->AddEntry(reco_pion_momentum_selectedCC1pi_100_200MeV, "100 to 200 MeV", "l");
    legend->AddEntry(reco_pion_momentum_selectedCC1pi_over_200MeV, "Over 200 MeV", "l");
    legend->AddEntry(reco_pion_momentum_selectedCC1pi, "All", "l");
    legend->Draw();

    c1->SaveAs("plots/pionThresholdStudy_reco_pion_momentum_selectedCC1pi_testingOnly_lowPiMomThreshold.pdf");

    // Delete the TCanvas
    delete c1;


    // Set the color and line thickness of the histograms
    reco_pion_momentum_selectedCC1pi_trueGolden->SetLineColor(kBlack);
    reco_pion_momentum_selectedCC1pi_trueGolden->SetLineWidth(3);
    reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->SetLineColor(kGreen);
    reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->SetLineWidth(3);
    reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden->SetLineColor(kBlue);
    reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden->SetLineWidth(3);
    reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden->SetLineColor(kRed);
    reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden->SetLineWidth(3);

    // Area normalize all four histograms
    reco_pion_momentum_selectedCC1pi_trueGolden->Scale(1.0 / reco_pion_momentum_selectedCC1pi_trueGolden->Integral());
    reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->Scale(1.0 / reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->Integral());
    reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden->Scale(1.0 / reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden->Integral());
    reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden->Scale(1.0 / reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden->Integral());

    // Remove the stats box
    reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->SetStats(0);

    // Draw the histogram
    TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
    // Draw the first histogram
    reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->Draw("E hist");

    // Find the maximum value across all histograms
    double max_val_trueGolden = reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->GetBinContent(reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->GetMaximumBin()) + reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->GetBinError(reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->GetMaximumBin());
    double temp_max_trueGolden = reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden->GetBinContent(reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden->GetMaximumBin()) + reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden->GetBinError(reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden->GetMaximumBin());
    if (temp_max_trueGolden > max_val_trueGolden) max_val_trueGolden = temp_max_trueGolden;
    temp_max_trueGolden = reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden->GetBinContent(reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden->GetMaximumBin()) + reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden->GetBinError(reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden->GetMaximumBin());
    if (temp_max_trueGolden > max_val_trueGolden) max_val_trueGolden = temp_max_trueGolden;
    temp_max_trueGolden = reco_pion_momentum_selectedCC1pi_trueGolden->GetBinContent(reco_pion_momentum_selectedCC1pi_trueGolden->GetMaximumBin()) + reco_pion_momentum_selectedCC1pi_trueGolden->GetBinError(reco_pion_momentum_selectedCC1pi_trueGolden->GetMaximumBin());
    if (temp_max_trueGolden > max_val_trueGolden) max_val_trueGolden = temp_max_trueGolden;

    // Set the y-axis range
    reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden->GetYaxis()->SetRangeUser(0, max_val_trueGolden * 1.1); // Add 10% for a little extra space

    // Draw the other histograms
    reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden->Draw("E hist same");
    reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden->Draw("E hist same");
    reco_pion_momentum_selectedCC1pi_trueGolden->Draw("E hist same");

    // Add a legend
    TLegend* legend_trueGolden = new TLegend(0.2, 0.2);
    legend_trueGolden->SetHeader("Selected True CC1pi Events", "C"); // "C" centers the header
    legend_trueGolden->AddEntry(reco_pion_momentum_selectedCC1pi_0_100MeV_trueGolden, "0 to 100 MeV", "l");
    legend_trueGolden->AddEntry(reco_pion_momentum_selectedCC1pi_100_200MeV_trueGolden, "100 to 200 MeV", "l");
    legend_trueGolden->AddEntry(reco_pion_momentum_selectedCC1pi_over_200MeV_trueGolden, "Over 200 MeV", "l");
    legend_trueGolden->AddEntry(reco_pion_momentum_selectedCC1pi_trueGolden, "All", "l");
    legend_trueGolden->Draw();

    c2->SaveAs("plots/pionThresholdStudy_reco_pion_momentum_selectedCC1pi_trueGolden_testingOnly_lowPiMomThreshold.pdf");

    // Delete the TCanvas
    delete c2;


    // MakePlot2D(recoPionMomentumVStruePionMomentum_selectedTrueCC1pi, "recoPionMomentumVStruePionMomentum_selectedTrueCC1pi", false, true, false);
    // MakePlot2D(recoPionMomentumVStruePionMomentum_selectedTrueCC1pi_minBins, "recoPionMomentumVStruePionMomentum_selectedTrueCC1pi_minBins", false, true, true);
}