#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>


void FillHistogram(TTree* tree, TH1D* hist, const std::string& variable, const std::string& condition, const float runWeight)
{
    TH1D* hist_tmp = (TH1D*)hist->Clone("hist_tmp");
    hist_tmp->Reset();

    std::string hist_tmp_name = hist_tmp->GetName();
    std::string weight = "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)*" + std::to_string(runWeight);

    // Fill the histograms with weights
    tree->Project(hist_tmp_name.c_str(), variable.c_str(), ("(" + condition + ")*" + weight).c_str());

    // Check if Sumw2 has already been called
    if (!hist_tmp->GetSumw2N()) {
        hist_tmp->Sumw2();
    }

    // Add the tmp histograms to the main histograms
    hist->Add(hist_tmp);

    delete hist_tmp;
}

void FillHistogram2D(TTree* tree, TH2D* hist, const std::string& variable1, const std::string& variable2, const std::string& condition, const float runWeight)
{
    TH2D* hist_tmp = (TH2D*)hist->Clone("hist2D_tmp");
    hist_tmp->Reset();

    std::string hist_tmp_name = hist_tmp->GetName();
    std::string weight = "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)*" + std::to_string(runWeight);

    // Fill the histograms with weights
    tree->Project(hist_tmp_name.c_str(), (variable2 + ":" + variable1).c_str(), ("(" + condition + ")*" + weight).c_str());

    // Check if Sumw2 has already been called
    if (!hist_tmp->GetSumw2N()) {
        hist_tmp->Sumw2();
    }

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

    c1->SaveAs(("plots/MuonMomentumVSAngleContainmentStudy_" + name + ".pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void MakeCombinedPlot(TH1D* histData, TH1D* histMC, const std::string& name)
{
    // Set the color and line thickness of the histograms
    histData->SetLineColor(kBlack);
    histData->SetLineWidth(2);
    histMC->SetLineColor(kRed);
    histMC->SetLineWidth(2);

    // Remove the stats box
    histData->SetStats(0);
    histMC->SetStats(0);

    // Draw the histograms
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    histData->Draw("E");
    histMC->Draw("E HIST SAME");

    // Add a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(histData, "Data", "l");
    legend->AddEntry(histMC, "MC", "l");
    legend->Draw();

    // Save the canvas
    c1->SaveAs(("plots/MuonMomentumVSAngleContainmentStudy_" + name + ".pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void MakePlot2D(TH2D* hist, const std::string& name, const bool drawDiagonal = false, const bool axisTicks = true, const bool drawText = false, const bool drawExponential = false, const bool altZContour = false, const std::vector<float>& binEdges = {})
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

    // Draw exponential function line if requested
    if (drawExponential) {
        TF1* expFunc = new TF1("expFunc", "0.0405714841414*exp(3.3*x) + 0.3", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
        expFunc->SetLineColor(kRed);   // Set line color to blue
        expFunc->SetLineStyle(1);       // Set line style to solid
        expFunc->SetLineWidth(2);       // Set line width
        expFunc->Draw("same");          // Draw the exponential function on the same canvas
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

    std::cout << "DEBUG - hist title 1: " << hist->GetTitle() << std::endl;

    if(altZContour)
    {
        const auto zMax = hist->GetMaximum();
        // Set the color palette for the z-axis
        const Int_t NRGBs = 5;
        const Int_t NCont = 255;
        Double_t stops[NRGBs] = { 0.00, 0.5/zMax, 1.0/zMax, 1.5/zMax, 1.00 };
        Double_t red[NRGBs]   = { 0.10, 0.10, 0.00, 1.00, 1.00 };
        Double_t green[NRGBs] = { 0.10, 0.10, 1.00, 0.20, 0.20 };
        Double_t blue[NRGBs]  = { 0.90, 0.90, 0.00, 0.00, 0.00 };;
        TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        // gStyle->SetNumberContours(NCont);
        hist->SetContour(NCont);
    }

    c1->SaveAs(("plots/MuonMomentumVSAngleContainmentStudy_" + name + ".pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void MuonMomentumVSAngleContainmentStudy()
{
    const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/1March24/";
    // // List of files
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {

        std::make_tuple("beam on", "1",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi.root", 1.f),
        std::make_tuple("beam on", "2",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 1.f),
        std::make_tuple("beam on", "3",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 1.f),
        std::make_tuple("beam on", "4bcd", rootPath + "bnb_beam_on_peleeTuple_uboone_run4bcd_ubcc1pi.root", 1.f),
        std::make_tuple("beam on", "5",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root", 1.f),

        std::make_tuple("beam off", "1",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi.root", 0.56421),
        std::make_tuple("beam off", "2",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 0.40202),
        std::make_tuple("beam off", "3",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 0.30657),
        std::make_tuple("beam off", "4bcd", rootPath + "bnb_beam_off_peleeTuple_uboone_run4bcd_ubcc1pi.root", 0.30175),
        std::make_tuple("beam off", "5",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root", 0.32807),

        std::make_tuple("dirt mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root", 0.52806),
        std::make_tuple("dirt mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root", 0.27521),
        std::make_tuple("dirt mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root", 0.80892),
        std::make_tuple("dirt mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi.root", 0.39701),
        std::make_tuple("dirt mc", "5",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root", 0.41280),

        std::make_tuple("nu mc", "1",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
        std::make_tuple("nu mc", "2",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root", 0.25750*2.0),
        std::make_tuple("nu mc", "3",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root", 0.20113*2.0),
        std::make_tuple("nu mc", "4bcd", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root", 0.13074*2.0),
        std::make_tuple("nu mc", "5",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196*2.0),
    };

    // 2d histograms
    const int nMomentumBins = 15;
    const int nAngleBins = 15;
    
    TH2D* h_muMomvsCosTheta_Contained_Data = new TH2D("muMomvsCosTheta_Contained_Data", "p_{#mu} vs cos(#theta_{#mu}) for Contained Selected Events (Data); Reco Cos(#theta_{#mu}; Reco p_{#mu}", nMomentumBins, -1, 1, nAngleBins, 0, 3);
    TH2D* h_muMomvsCosTheta_Uncontained_Data = new TH2D("muMomvsCosTheta_Uncontained_Data", "p_{#mu} vs cos(#theta_{#mu}) for Uncontained Selected Events (Data); Reco Cos(#theta_{#mu}; Reco p_{#mu}", nMomentumBins, -1, 1, nAngleBins, 0, 3);
    TH2D* h_muMomvsCosTheta_Contained_Data_expPhaseSpaceCut = new TH2D("muMomvsCosTheta_Contained_Data_expPhaseSpaceCut", "p_{#mu} vs cos(#theta_{##mu}) for Contained Selected Events (Data) with phase space cut; Reco Cos(#theta_{#mu}; Reco p_{#mu}", nMomentumBins, -1, 1, nAngleBins, 0, 3);
    TH2D* h_muMomvsCosTheta_Uncontained_Data_expPhaseSpaceCut = new TH2D("muMomvsCosTheta_Uncontained_Data_expPhaseSpaceCut", "p_{#mu} vs cos(#theta_{#mu}) for Uncontained Selected Events (Data) with phase space cut; Reco Cos(#theta_{#mu}; Reco p_{#mu}", nMomentumBins, -1, 1, nAngleBins, 0, 3);
    TH2D* h_muMomvsCosTheta_Contained_MC = new TH2D("muMomvsCosTheta_Contained_MC", "p_{#mu} vs cos(#theta_{#mu}) for Contained Selected Events (MC); Reco Cos(#theta_{#mu}; Reco p_{#mu}", nMomentumBins, -1, 1, nAngleBins, 0, 3);
    TH2D* h_muMomvsCosTheta_Uncontained_MC = new TH2D("muMomvsCosTheta_Uncontained_MC", "p_{#mu} vs cos(#theta_{#mu}) for Uncontained Selected Events (MC); Reco Cos(#theta_{#mu}; Reco p_{#mu}", nMomentumBins, -1, 1, nAngleBins, 0, 3);
    TH2D* h_muMomvsCosTheta_Contained_MC_expPhaseSpaceCut = new TH2D("muMomvsCosTheta_Contained_MC_expPhaseSpaceCut", "p_{#mu} vs cos(#theta_{#mu}) for Contained Selected Events (MC) with phase space cut; Reco Cos(#theta_{#mu}; Reco p_{#mu}", nMomentumBins, -1, 1, nAngleBins, 0, 3);
    TH2D* h_muMomvsCosTheta_Uncontained_MC_expPhaseSpaceCut = new TH2D("muMomvsCosTheta_Uncontained_MC_expPhaseSpaceCut", "p_{#mu} vs cos(#theta_{#mu}) for Uncontained Selected Events (MC) with phase space cut; Reco Cos(#theta_{#mu}; Reco p_{#mu}", nMomentumBins, -1, 1, nAngleBins, 0, 3);

    TH1D* h_muonPhi_Uncontained_Data = new TH1D("muonPhi_Uncontained_Data", "Muon Phi for Uncontained Events; Reco #phi_{#mu}", 10, -3.14159, 3.14159);
    TH1D* h_muonMom_Uncontained_Data = new TH1D("muonMom_Uncontained_Data", "Muon Momentum for Uncontained Events; Reco p_{#mu}", 20, 0, 3);
    TH1D* h_muonPhi_Uncontained_MC = new TH1D("muonPhi_Uncontained_MC", "Muon Phi for Uncontained Events; Reco #phi_{#mu}", 10, -3.14159, 3.14159);
    TH1D* h_muonMom_Uncontained_MC = new TH1D("muonMom_Uncontained_MC", "Muon Momentum for Uncontained Events; Reco p_{#mu}", 20, 0, 3);
    
    TH1D* h_muonPhi_Uncontained_Data_expPhaseSpaceCut = new TH1D("muonPhi_Uncontained_Data_expPhaseSpaceCut", "Muon Phi for Uncontained Events with phase space cut; Reco #phi_{#mu}", 10, -3.14159, 3.14159);
    TH1D* h_muonMom_Uncontained_Data_expPhaseSpaceCut = new TH1D("muonMom_Uncontained_Data_expPhaseSpaceCut", "Muon Momentum for Uncontained Events with phase space cut; Reco p_{#mu}", 20, 0, 3);
    TH1D* h_muonPhi_Uncontained_MC_expPhaseSpaceCut = new TH1D("muonPhi_Uncontained_MC_expPhaseSpaceCut", "Muon Phi for Uncontained Events with phase space cut; Reco #phi_{#mu}", 10, -3.14159, 3.14159);
    TH1D* h_muonMom_Uncontained_MC_expPhaseSpaceCut = new TH1D("muonMom_Uncontained_MC_expPhaseSpaceCut", "Muon Momentum for Uncontained Events with phase space cut; Reco p_{#mu}", 20, 0, 3);

    TH1D* h_muonPhi_Contained_Data = new TH1D("muonPhi_Contained_Data", "Muon Phi for Contained Events; Reco #phi_{#mu}", 10, -3.14159, 3.14159);
    TH1D* h_muonMom_Contained_Data = new TH1D("muonMom_Contained_Data", "Muon Momentum for Contained Events; Reco p_{#mu}", 20, 0, 3);
    TH1D* h_muonPhi_Concontained_MC = new TH1D("muonPhi_Concontained_MC", "Muon Phi for Contained Events; Reco #phi_{#mu}", 10, -3.14159, 3.14159);
    TH1D* h_muonMom_Contained_MC = new TH1D("muonMom_Contained_MC", "Muon Momentum for Contained Events; Reco p_{#mu}", 20, 0, 3);

    TH1D* h_muonPhi_Contained_Data_expPhaseSpaceCut = new TH1D("muonPhi_Contained_Data_expPhaseSpaceCut", "Muon Phi for Contained Events with phase space cut; Reco #phi_{#mu}", 10, -3.14159, 3.14159);
    TH1D* h_muonMom_Contained_Data_expPhaseSpaceCut = new TH1D("muonMom_Contained_Data_expPhaseSpaceCut", "Muon Momentum for Contained Events with phase space cut; Reco p_{#mu}", 20, 0, 3);
    TH1D* h_muonPhi_Contained_MC_expPhaseSpaceCut = new TH1D("muonPhi_Contained_MC_expPhaseSpaceCut", "Muon Phi for Contained Events with phase space cut; Reco #phi_{#mu}", 10, -3.14159, 3.14159);
    TH1D* h_muonMom_Contained_MC_expPhaseSpaceCut = new TH1D("muonMom_Contained_MC_expPhaseSpaceCut", "Muon Momentum for Contained Events with phase space cut; Reco p_{#mu}", 20, 0, 3);

    // Loop over the files
    for (int f = 0; f < files.size(); f++)
    {
        const auto [sampleType, run, filePath, runWeight] = files.at(f);

        // if(run != "1") continue;

        std::cout<<"DEBUG - filePath: "<<filePath<<std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(filePath.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");

        if(sampleType == "beam on")
        {
            FillHistogram2D(tree, h_muMomvsCosTheta_Contained_Data,
            /*variable1*/ "cc1pi_reco_muonCosTheta",
            /*variable2*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram2D(tree, h_muMomvsCosTheta_Uncontained_Data,
            /*variable1*/ "cc1pi_reco_muonCosTheta",
            /*variable2*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram2D(tree, h_muMomvsCosTheta_Contained_Data_expPhaseSpaceCut,
            /*variable1*/ "cc1pi_reco_muonCosTheta",
            /*variable2*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);

            FillHistogram2D(tree, h_muMomvsCosTheta_Uncontained_Data_expPhaseSpaceCut,
            /*variable1*/ "cc1pi_reco_muonCosTheta",
            /*variable2*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);


            FillHistogram(tree, h_muonPhi_Uncontained_Data,
            /*variable*/ "cc1pi_reco_muonPhi",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram(tree, h_muonMom_Uncontained_Data,
            /*variable*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram(tree, h_muonPhi_Uncontained_Data_expPhaseSpaceCut,
            /*variable*/ "cc1pi_reco_muonPhi",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);

            FillHistogram(tree, h_muonMom_Uncontained_Data_expPhaseSpaceCut,
            /*variable*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);
            

            FillHistogram(tree, h_muonPhi_Contained_Data,
            /*variable*/ "cc1pi_reco_muonPhi",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram(tree, h_muonMom_Contained_Data,
            /*variable*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram(tree, h_muonPhi_Contained_Data_expPhaseSpaceCut,
            /*variable*/ "cc1pi_reco_muonPhi",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);

            FillHistogram(tree, h_muonMom_Contained_Data_expPhaseSpaceCut,
            /*variable*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);
        }
        else
        {
            FillHistogram2D(tree, h_muMomvsCosTheta_Contained_MC,
            /*variable1*/ "cc1pi_reco_muonCosTheta",
            /*variable2*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram2D(tree, h_muMomvsCosTheta_Uncontained_MC,
            /*variable1*/ "cc1pi_reco_muonCosTheta",
            /*variable2*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram2D(tree, h_muMomvsCosTheta_Contained_MC_expPhaseSpaceCut,
            /*variable1*/ "cc1pi_reco_muonCosTheta",
            /*variable2*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);

            FillHistogram2D(tree, h_muMomvsCosTheta_Uncontained_MC_expPhaseSpaceCut,
            /*variable1*/ "cc1pi_reco_muonCosTheta",
            /*variable2*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);


            FillHistogram(tree, h_muonPhi_Uncontained_MC,
            /*variable*/ "cc1pi_reco_muonPhi",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram(tree, h_muonMom_Uncontained_MC,
            /*variable*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram(tree, h_muonPhi_Uncontained_MC_expPhaseSpaceCut,
            /*variable*/ "cc1pi_reco_muonPhi",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);

            FillHistogram(tree, h_muonMom_Uncontained_MC_expPhaseSpaceCut,
            /*variable*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && !cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);


            FillHistogram(tree, h_muonPhi_Concontained_MC,
            /*variable*/ "cc1pi_reco_muonPhi",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram(tree, h_muonMom_Contained_MC,
            /*variable*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained", runWeight);

            FillHistogram(tree, h_muonPhi_Contained_MC_expPhaseSpaceCut,
            /*variable*/ "cc1pi_reco_muonPhi",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);

            FillHistogram(tree, h_muonMom_Contained_MC_expPhaseSpaceCut,
            /*variable*/ "cc1pi_reco_muonMomentum",
            /*Condition*/ "cc1pi_selected_generic && cc1pi_reco_pionMomentum >= 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum < (0.0405714841414*exp(3.3 * cc1pi_reco_muonCosTheta) + 0.3)", runWeight);
        }

        tFile->Close();
    }
    std::cout<<"DEBUG - Done with file loop"<<std::endl;

    MakePlot2D(h_muMomvsCosTheta_Contained_Data, "muMomvsCosTheta_Contained_Data", false, true, true, true);
    MakePlot2D(h_muMomvsCosTheta_Uncontained_Data, "muMomvsCosTheta_Uncontained_Data", false, true, true, true);
    MakePlot2D(h_muMomvsCosTheta_Contained_Data_expPhaseSpaceCut, "muMomvsCosTheta_Contained_Data_expPhaseSpaceCut", false, true, true, true);
    MakePlot2D(h_muMomvsCosTheta_Uncontained_Data_expPhaseSpaceCut, "muMomvsCosTheta_Uncontained_Data_expPhaseSpaceCut", false, true, true, true);
    MakePlot2D(h_muMomvsCosTheta_Contained_MC, "muMomvsCosTheta_Contained_MC", false, true, true, true);
    MakePlot2D(h_muMomvsCosTheta_Uncontained_MC, "muMomvsCosTheta_Uncontained_MC", false, true, true, true);
    MakePlot2D(h_muMomvsCosTheta_Contained_MC_expPhaseSpaceCut, "muMomvsCosTheta_Contained_MC_expPhaseSpaceCut", false, true, true, true);
    MakePlot2D(h_muMomvsCosTheta_Uncontained_MC_expPhaseSpaceCut, "muMomvsCosTheta_Uncontained_MC_expPhaseSpaceCut", false, true, true, true);

    // Create a new histogram for the ratio
    TH2D* h_muMomvsCosTheta_Ratio_Data_Data = (TH2D*)h_muMomvsCosTheta_Contained_Data->Clone("muMomvsCosTheta_Ratio");
    h_muMomvsCosTheta_Ratio_Data_Data->SetTitle("Contained / Uncontained (Data); Reco Cos(#theta_{#mu}); Reco p_{#mu}");
    // Compute the element-wise ratio
    h_muMomvsCosTheta_Ratio_Data_Data->Divide(h_muMomvsCosTheta_Uncontained_Data);    
    // Optionally, you can plot the ratio histogram
    MakePlot2D(h_muMomvsCosTheta_Ratio_Data_Data, "muMomvsCosTheta_Ratio_Data_Data", false, true, true, true, true);

    // Create a new histogram for the ratio (MC with phase space cut)
    TH2D* h_muMomvsCosTheta_Ratio_MC_MC = (TH2D*)h_muMomvsCosTheta_Contained_MC->Clone("muMomvsCosTheta_Ratio_MC_MC");
    h_muMomvsCosTheta_Ratio_MC_MC->SetTitle("Contained / Uncontained (MC); Reco Cos(#theta_{#mu}); Reco p_{#mu}");
    // Compute the element-wise ratio
    h_muMomvsCosTheta_Ratio_MC_MC->Divide(h_muMomvsCosTheta_Uncontained_MC);
    // Optionally, you can plot the ratio histogram
    MakePlot2D(h_muMomvsCosTheta_Ratio_MC_MC, "muMomvsCosTheta_Ratio_MC_MC", false, true, true, true, true);
    
    // Create a new histogram for the ratio (Data/MC for contained)
    TH2D* h_muMomvsCosTheta_Ratio_Contained_Data_MC = (TH2D*)h_muMomvsCosTheta_Contained_Data->Clone("muMomvsCosTheta_Ratio_Contained_Data_MC");
    h_muMomvsCosTheta_Ratio_Contained_Data_MC->SetTitle("Contained (Data/MC); Reco Cos(#theta_{#mu}); Reco p_{#mu}");
    // Compute the element-wise ratio
    h_muMomvsCosTheta_Ratio_Contained_Data_MC->Divide(h_muMomvsCosTheta_Contained_MC);
    // Optionally, you can plot the ratio histogram
    MakePlot2D(h_muMomvsCosTheta_Ratio_Contained_Data_MC, "muMomvsCosTheta_Ratio_Contained_Data_MC", false, true, true, true, true);

    // Create a new histogram for the ratio (Data/MC for uncontained)
    TH2D* h_muMomvsCosTheta_Ratio_Uncontained_Data_MC = (TH2D*)h_muMomvsCosTheta_Uncontained_Data->Clone("muMomvsCosTheta_Ratio_Uncontained_Data_MC");
    h_muMomvsCosTheta_Ratio_Uncontained_Data_MC->SetTitle("Uncontained (Data/MC); Reco Cos(#theta_{#mu}); Reco p_{#mu}");
    // Compute the element-wise ratio
    h_muMomvsCosTheta_Ratio_Uncontained_Data_MC->Divide(h_muMomvsCosTheta_Uncontained_MC);
    // Optionally, you can plot the ratio histogram
    MakePlot2D(h_muMomvsCosTheta_Ratio_Uncontained_Data_MC, "muMomvsCosTheta_Ratio_Uncontained_Data_MC", false, true, true, true, true);

    // Create 1D plots for uncontained data and MC
    MakeCombinedPlot(h_muonPhi_Uncontained_Data, h_muonPhi_Uncontained_MC, "muonPhi_Uncontained");
    MakeCombinedPlot(h_muonMom_Uncontained_Data, h_muonMom_Uncontained_MC, "muonMom_Uncontained");

    // Create 1D plots for uncontained data and MC with exponential phase space cut
    MakeCombinedPlot(h_muonPhi_Uncontained_Data_expPhaseSpaceCut, h_muonPhi_Uncontained_MC_expPhaseSpaceCut, "muonPhi_Uncontained_expPhaseSpaceCut");
    MakeCombinedPlot(h_muonMom_Uncontained_Data_expPhaseSpaceCut, h_muonMom_Uncontained_MC_expPhaseSpaceCut, "muonMom_Uncontained_expPhaseSpaceCut");

    // Create 1D plots for contained data and MC
    MakeCombinedPlot(h_muonPhi_Contained_Data, h_muonPhi_Concontained_MC, "muonPhi_Contained");
    MakeCombinedPlot(h_muonMom_Contained_Data, h_muonMom_Contained_MC, "muonMom_Contained");

    // Create 1D plots for contained data and MC with exponential phase space cut
    MakeCombinedPlot(h_muonPhi_Contained_Data_expPhaseSpaceCut, h_muonPhi_Contained_MC_expPhaseSpaceCut, "muonPhi_Contained_expPhaseSpaceCut");
    MakeCombinedPlot(h_muonMom_Contained_Data_expPhaseSpaceCut, h_muonMom_Contained_MC_expPhaseSpaceCut, "muonMom_Contained_expPhaseSpaceCut");

    // Delete the 2D histograms
    delete h_muMomvsCosTheta_Contained_Data;
    delete h_muMomvsCosTheta_Uncontained_Data;
    delete h_muMomvsCosTheta_Contained_MC;
    delete h_muMomvsCosTheta_Uncontained_MC;
    delete h_muMomvsCosTheta_Contained_MC_expPhaseSpaceCut;
    delete h_muMomvsCosTheta_Uncontained_MC_expPhaseSpaceCut;
    
    // Delete the 1D histograms for uncontained data and MC
    delete h_muonPhi_Uncontained_Data;
    delete h_muonMom_Uncontained_Data;
    delete h_muonPhi_Uncontained_MC;
    delete h_muonMom_Uncontained_MC;
    
    // Delete the 1D histograms for uncontained data and MC with exponential phase space cut
    delete h_muonPhi_Uncontained_Data_expPhaseSpaceCut;
    delete h_muonMom_Uncontained_Data_expPhaseSpaceCut;
    delete h_muonPhi_Uncontained_MC_expPhaseSpaceCut;
    delete h_muonMom_Uncontained_MC_expPhaseSpaceCut;
    
    // Delete the 1D histograms for contained data and MC
    delete h_muonPhi_Contained_Data;
    delete h_muonMom_Contained_Data;
    delete h_muonPhi_Concontained_MC;
    delete h_muonMom_Contained_MC;
    
    // Delete the 1D histograms for contained data and MC with exponential phase space cut
    delete h_muonPhi_Contained_Data_expPhaseSpaceCut;
    delete h_muonMom_Contained_Data_expPhaseSpaceCut;
    delete h_muonPhi_Contained_MC_expPhaseSpaceCut;
    delete h_muonMom_Contained_MC_expPhaseSpaceCut;

    delete h_muMomvsCosTheta_Ratio_Data_Data;
    delete h_muMomvsCosTheta_Ratio_MC_MC;
    delete h_muMomvsCosTheta_Ratio_Contained_Data_MC;
    delete h_muMomvsCosTheta_Ratio_Uncontained_Data_MC;
}