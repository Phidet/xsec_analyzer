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

    hist_tmp->Sumw2();

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

    c1->SaveAs(("plots/muonMomentumStudy_" + name + "_testingOnly_lowPiMomThreshold.pdf").c_str());

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

    c1->SaveAs(("plots/muonMomentumStudy_" + name + "_testingOnly_lowPiMomThreshold.pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void MuonMomentumStudy() 
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
    // TH2D* trueVsRecoMuonMom_selectedTrueCC1pi = new TH2D("muonMomVSpionMom_selectedTrueCC1pi", "Selected True CC1pi Events (and muon candidate has backtracked truth particle); True Muon Momentum (GeV/c); True Pion Momentum (GeV/c)", 30, 0, 1.5, 30, 0, 1.5);
    // TH2D* trueVsRecoMuonMom_TrueCC1pi = new TH2D("trueVsRecoMuonMom_TrueCC1pi", "True CC1pi; True Muon Momentum (GeV/c); Reco Muon Momentum (GeV/c)", 50, 0, 2, 50, 0, 2);
    TH2D* trueVsRecoMuonMom_selectedTrueCC1pi = new TH2D("trueVsRecoMuonMom_selectedTrueCC1pi", "Selected True CC1pi; True Muon Momentum (GeV/c); Reco Muon Momentum (GeV/c)", 50, 0, 2, 50, 0, 2);


    // Loop over the files
    for (int f = 0; f < files.size(); f++)
    {
        const auto [sampleType, run, filePath, runWeight] = files.at(f);
        std::cout<<"DEBUG - filePath: "<<filePath<<std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(filePath.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");


        // FillHistogram2D(tree, trueVsRecoMuonMom_TrueCC1pi, 
        // /*variable1*/ "cc1pi_truth_muonMomentum",
        // /*variable2*/ "cc1pi_reco_muonMomentum", 
        // /*Condition*/ "cc1pi_signal  &&  cc1pi_truth_pionMomentum > 0.1", runWeight);

        // std::cout<<"DEBUG - Point X0.2"<<std::endl;

        FillHistogram2D(tree, trueVsRecoMuonMom_selectedTrueCC1pi, 
        /*variable1*/ "cc1pi_truth_muonMomentum",
        /*variable2*/ "cc1pi_reco_muonMomentum", 
        /*Condition*/ "cc1pi_signal  &&  cc1pi_truth_pionMomentum > 0.1 && cc1pi_selected_generic && cc1pi_reco_pionMomentum > 0.1", runWeight);


        tFile->Close();
    }
    std::cout<<"DEBUG - Done with file loop"<<std::endl;

    // MakePlot2D(trueVsRecoMuonMom_TrueCC1pi, "trueVsRecoMuonMom_TrueCC1pi");
    MakePlot2D(trueVsRecoMuonMom_selectedTrueCC1pi, "trueVsRecoMuonMom_selectedTrueCC1pi");

}