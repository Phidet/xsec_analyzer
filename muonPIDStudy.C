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

    c1->SaveAs(("plots/muonPIDStudy_" + name + ".pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void MakePlot2D(TH2D* hist, const std::string& name, const bool drawDiagonal = false, const bool axisTicks = true)
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

    hist->Draw("COLZ");

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

    c1->SaveAs(("plots/muonPIDStudy_" + name + ".pdf").c_str());

    // Delete the TCanvas
    delete c1;
}

void muonPIDStudy() 
{
    // List of files
    const std::string rootPath = "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/22Feb24/";
    // tuple: type, run, file path, run weight
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        std::make_tuple("nu mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root", 0.13011),
        std::make_tuple("nu mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root", 0.25750),
        std::make_tuple("nu mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root", 0.20113),
        std::make_tuple("nu mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi.root", 0.13074),
        std::make_tuple("nu mc", "5",  rootPath + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root", 0.15196),
    };

    // const std::vector<float> overlayPOT {1.28278e+21, 1.01592e+21, 1.31355e+21, 2.48772e+20, 7.64399e+20, 4.64842e+20, 8.66958e+20, 9.70631e+20}; // run 1, 2, 3, 4a, 4b, 4c, 4d, 5
    // const std::vector<float> dataBNBTor860WCut {1.669e+20, 2.616e+20, 2.562e+20, 3.607e+19, 1.39e+20, 8.586e+19, 4.919e+19, 1.23e+20};
    // std::vector<float> ratioPOT;
    // for (size_t i = 0; i < overlayPOT.size(); i++)
    // {
    //     ratioPOT.push_back(dataBNBTor860WCut[i] / overlayPOT[i]);
    // }
    // const std::vector<float> weights = {ratioPOT.at(0), ratioPOT.at(1), ratioPOT.at(2)};

    // 1d histograms
    // TH1D* pdg_backtrackedMuon = new TH1D("cc1pi_backtracked_muonPDG", "cc1pi_backtracked_muonPDG", 56, -14, 14);
    TH1D* trueCC1pi_muonBDT_trueMuon = new TH1D("trueCC1pi_muonBDT_trueMuon", "Particles in True CC1pi Events; Muon BDT Score; # Events (Area normalised)", 20, -1.1, 1.1);
    TH1D* trueCC1pi_muonBDT_trueProton = new TH1D("trueCC1pi_muonBDT_trueProton", "trueCC1pi_muonBDT_trueProton", 20, -1.1, 1.1);
    TH1D* trueCC1pi_muonBDT_truePion = new TH1D("trueCC1pi_muonBDT_truePion", "trueCC1pi_muonBDT_truePion", 20, -1.1, 1.1);

    TH1D* trueCC1pi_protonBDT_trueMuon = new TH1D("trueCC1pi_protonBDT_trueMuon", "Particles in True CC1pi Events; Proton BDT Score; # Events (Area normalised)", 20, -1.1, 1.1);
    TH1D* trueCC1pi_protonBDT_trueProton = new TH1D("trueCC1pi_protonBDT_trueProton", "trueCC1pi_protonBDT_trueProton", 20, -1.1, 1.1);
    TH1D* trueCC1pi_protonBDT_truePion = new TH1D("trueCC1pi_protonBDT_truePion", "trueCC1pi_protonBDT_truePion", 20, -1.1, 1.1);

    TH1D* trueCC1pi_muonBDT_trueMuon_contained = new TH1D("trueCC1pi_muonBDT_trueMuon_contained", "Contained Particles in True CC1pi Events; Muon BDT Score; # Events (Area normalised)", 20, -1.1, 1.1);
    TH1D* trueCC1pi_muonBDT_trueProton_contained = new TH1D("trueCC1pi_muonBDT_trueProton_contained", "trueCC1pi_muonBDT_trueProton_contained", 20, -1.1, 1.1);
    TH1D* trueCC1pi_muonBDT_truePion_contained = new TH1D("trueCC1pi_muonBDT_truePion_contained", "trueCC1pi_muonBDT_truePion_contained", 20, -1.1, 1.1);

    TH1D* trueCC1pi_muonBDT_trueMuon_uncontained = new TH1D("trueCC1pi_muonBDT_trueMuon_uncontained", "Escaping Particles in True CC1pi Events; Muon BDT Score; # Events (Area normalised)", 20, -1.1, 1.1);
    TH1D* trueCC1pi_muonBDT_trueProton_uncontained = new TH1D("trueCC1pi_muonBDT_trueProton_uncontained", "trueCC1pi_muonBDT_trueProton_uncontained", 20, -1.1, 1.1);
    TH1D* trueCC1pi_muonBDT_truePion_uncontained = new TH1D("trueCC1pi_muonBDT_truePion_uncontained", "trueCC1pi_muonBDT_truePion_uncontained", 20, -1.1, 1.1);


    TH1D* trueCC1pi_protonBDT_trueMuon_contained = new TH1D("trueCC1pi_protonBDT_trueMuon_contained", "Contained Particles in True CC1pi Events; Proton BDT Score; # Events (Area normalised)", 20, -1.1, 1.1);
    TH1D* trueCC1pi_protonBDT_trueProton_contained = new TH1D("trueCC1pi_protonBDT_trueProton_contained", "trueCC1pi_protonBDT_trueProton_contained", 20, -1.1, 1.1);
    TH1D* trueCC1pi_protonBDT_truePion_contained = new TH1D("trueCC1pi_protonBDT_truePion_contained", "trueCC1pi_protonBDT_truePion_contained", 20, -1.1, 1.1);

    // 2d histograms
    TH2D* muonMomVSpionMom_selectedTrueCC1pi = new TH2D("muonMomVSpionMom_selectedTrueCC1pi", "Selected True CC1pi Events (and muon candidate has backtracked truth particle); True Muon Momentum (GeV/c); True Pion Momentum (GeV/c)", 20, 0, 1.5, 20, 0, 1.5);
    TH2D* muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon = new TH2D("muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon", "Fraction of selected true CC1pi events with correctly identified muon; True Muon Momentum (GeV/c); True Pion Momentum (GeV/c)", 20, 0, 1.5, 20, 0, 1.5);

    TH2D* muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi = new TH2D("muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi", "Selected True CC1pi Events (and muon candidate has backtracked truth particle); Muon Track Length > Pion Track Length; Muon Is Contained", 2, 0, 2, 2, 0, 2);
    TH2D* muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon = new TH2D("muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon", "Fraction of selected true CC1pi events with correctly identified muon; Muon Track Length > Pion Track Length; Muon Is Contained", 2, 0, 2, 2, 0, 2);

    TH2D* muonBDTVSrange_trueMuon_trueCC1pi = new TH2D("muonBDTVSrange_trueMuon_trueCC1pi", "Backtracked Muon Particles in True CC1pi Events; Reconstructed track length / cm; Muon BDT Score", 15, 0, 500, 15, -1, 0.8);
    TH2D* muonBDTVSTrueRange_trueMuon_trueCC1pi = new TH2D("muonBDTVSTruerange_trueMuon_trueCC1pi", "Backtracked Muon Particles in True CC1pi Events; True track length / cm; Muon BDT Score", 15, 0, 500, 15, -1, 0.8);

    TH2D* muonBDTVSrange_truePion_trueCC1pi = new TH2D("muonBDTVSrange_truePion_trueCC1pi", "Backtracked Pion Particles in True CC1pi Events; Reconstructed track length / cm; Muon BDT Score", 15, 0, 150, 15, -1, 0.8);
    TH2D* muonBDTVSTrueRange_truePion_trueCC1pi = new TH2D("muonBDTVSTruerange_truePion_trueCC1pi", "Backtracked Pion Particles in True CC1pi Events; True track length / cm; Muon BDT Score", 15, 0, 150, 15, -1, 0.8);

    // Loop over the files
    for (int f = 0; f < files.size(); f++)
    {
        const auto [sampleType, run, filePath, runWeight] = files.at(f);
        // Open the file and get the tree
        TFile* tFile = TFile::Open(filePath.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");

        // Disable all branches
        tree->SetBranchStatus("*", 0);

        // Enable only the branches you need
        // tree->SetBranchStatus("cc0pi_*", 1);
        tree->SetBranchStatus("trueCC1pi_*", 1);
        tree->SetBranchStatus("cc1pi_*", 1);
        // tree->SetBranchStatus("true_cc0pi", 1);
        tree->SetBranchStatus("spline_weight", 1);
        tree->SetBranchStatus("tuned_cv_weight", 1);

        // FillHistogram(tree, pdg_backtrackedMuon,
        // /*variable*/ "cc1pi_backtracked_muonPDG",
        // /*Condition*/ "cc1pi_signal && cc1pi_selected_generic", runWeight);

        FillHistogram2D(tree, muonMomVSpionMom_selectedTrueCC1pi, 
        /*variable1*/ "cc1pi_truth_muonMomentum",
        /*variable2*/ "cc1pi_truth_pionMomentum", 
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_backtracked_muonPDG > -2147483648", runWeight);
        FillHistogram2D(tree, muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon, 
        /*variable1*/ "cc1pi_truth_muonMomentum",
        /*variable2*/ "cc1pi_truth_pionMomentum", 
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && (cc1pi_backtracked_muonPDG == 13 ||  cc1pi_backtracked_muonPDG == -13)", runWeight);

        FillHistogram2D(tree, muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi, 
        /*variable1*/ "cc1pi_truthMuon_TrackLength>cc1pi_truthPion_TrackLength",
        /*variable2*/ "cc1pi_truthMuon_IsContained", 
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && cc1pi_backtracked_muonPDG > -2147483648", runWeight);

        FillHistogram2D(tree, muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon,
        /*variable1*/ "cc1pi_truthMuon_TrackLength>cc1pi_truthPion_TrackLength",
        /*variable2*/ "cc1pi_truthMuon_IsContained", 
        /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && (cc1pi_backtracked_muonPDG == 13 ||  cc1pi_backtracked_muonPDG == -13)", runWeight);

        // FillHistogram2D(tree, muonBDTVSrange_trueMuon_trueCC1pi, 
        // /*variable1*/ "cc1pi_recoMuon_TrackLength",
        // /*variable2*/ "cc1pi_recoMuon_muonBDTScore", 
        // /*Condition*/ "cc1pi_signal && cc1pi_selected_generic && (cc1pi_backtracked_muonPDG == 13 ||  cc1pi_backtracked_muonPDG == -13)", runWeight);

        // Particle plots
        // trueCC1pi_muonBDT
        // FillParticleHistogram<double>(tree, trueCC1pi_muonBDT_trueMuon, "trueCC1pi_recoParticle_muonBDTScore_vector", "cc1pi_signal", runWeight);
        // FillParticleHistogram<double>(tree, trueCC1pi_muonBDT_trueProton, "trueCC1pi_recoParticle_muonBDTScore_vector", "cc1pi_signal", runWeight);
        // FillParticleHistogram<double>(tree, trueCC1pi_muonBDT_truePion, "trueCC1pi_recoParticle_muonBDTScore_vector", "cc1pi_signal", runWeight);



        std::vector<int> *pBacktracked_absTruePDG = nullptr;
        std::vector<bool> *pIsContained = nullptr;
        std::vector<bool> *pBacktracked_isContained = nullptr;
        std::vector<double> *pTrackLength = nullptr;
        std::vector<double> *pBacktracked_TrackLength = nullptr;
        std::vector<double> *pProtonBDTScore = nullptr;
        std::vector<double> *pMuonBDTScore = nullptr;

        tree->SetBranchAddress("trueCC1pi_recoParticle_backtracked_absTruePDG_vector", &pBacktracked_absTruePDG);
        tree->SetBranchAddress("trueCC1pi_recoParticle_isContained_vector", &pIsContained);
        tree->SetBranchAddress("trueCC1pi_recoParticle_backtracked_isContained_vector", &pBacktracked_isContained);
        tree->SetBranchAddress("trueCC1pi_recoParticle_TrackLength_vector", &pTrackLength);
        tree->SetBranchAddress("trueCC1pi_recoParticle_backtracked_TrackLength_vector", &pBacktracked_TrackLength);
        tree->SetBranchAddress("trueCC1pi_recoParticle_protonBDTScore_vector", &pProtonBDTScore);
        tree->SetBranchAddress("trueCC1pi_recoParticle_muonBDTScore_vector", &pMuonBDTScore);

        Float_t spline_weight, tuned_cv_weight;
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);

        // Check if values is null
        if (!pMuonBDTScore || !pProtonBDTScore) throw std::runtime_error("Error: Branchnot found in the tree.");

        std::cout << "Warning - Only using 3\% of the events for testing purposes." << std::endl;
        const auto nEvents = tree->GetEntries()/33;
        // Loop over the entries in the tree
        for (Long64_t i = 0; i < nEvents; i++)
        {
            tree->GetEntry(i);

            // Update the progress bar at every percent
            if (i % (nEntries / 100) == 0)
            {
                const auto progress = static_cast<int>(100.0 * i / nEntries);
                std::cout << "\r[" << std::string(progress, '|') << std::string(100 - progress, ' ') << "] " << progress << "%" << std::flush;
            }

            auto eventWeight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;
            eventWeight *= runWeight; // Apply the run weight

            // std::cout << "DEBUG FillParticleHistogram Point pMuonBDTScore->size(): " << pMuonBDTScore->size() << std::endl;
            // Apply the condition
            if (tree->GetLeaf("cc1pi_signal")->GetValue() == true)
            {
                // Loop over the values in the vector and fill the histogram
                for (Long64_t v = 0; v< pProtonBDTScore->size(); v++)
                {
                    if (pBacktracked_absTruePDG->at(v) == 2212) trueCC1pi_protonBDT_trueProton->Fill(pProtonBDTScore->at(v), eventWeight);
                    if (pBacktracked_absTruePDG->at(v) == 211) trueCC1pi_protonBDT_truePion->Fill(pProtonBDTScore->at(v), eventWeight);
                    if (pBacktracked_absTruePDG->at(v) == 13) trueCC1pi_protonBDT_trueMuon->Fill(pProtonBDTScore->at(v), eventWeight);

                    if(pBacktracked_isContained->at(v))
                    {
                        if (pBacktracked_absTruePDG->at(v) == 2212) trueCC1pi_protonBDT_trueProton_contained->Fill(pProtonBDTScore->at(v), eventWeight);
                        if (pBacktracked_absTruePDG->at(v) == 211) trueCC1pi_protonBDT_truePion_contained->Fill(pProtonBDTScore->at(v), eventWeight);
                        if (pBacktracked_absTruePDG->at(v) == 13) trueCC1pi_protonBDT_trueMuon_contained->Fill(pProtonBDTScore->at(v), eventWeight);
                    }
                }

                for (Long64_t v = 0; v< pMuonBDTScore->size(); v++)
                {
                    if (pBacktracked_absTruePDG->at(v) == 2212) trueCC1pi_muonBDT_trueProton->Fill(pMuonBDTScore->at(v), eventWeight);
                    if (pBacktracked_absTruePDG->at(v) == 211) trueCC1pi_muonBDT_truePion->Fill(pMuonBDTScore->at(v), eventWeight);
                    if (pBacktracked_absTruePDG->at(v) == 13) trueCC1pi_muonBDT_trueMuon->Fill(pMuonBDTScore->at(v), eventWeight);

                    if(pBacktracked_isContained->at(v))
                    {
                        if (pBacktracked_absTruePDG->at(v) == 2212) trueCC1pi_muonBDT_trueProton_contained->Fill(pMuonBDTScore->at(v), eventWeight);
                        if (pBacktracked_absTruePDG->at(v) == 211) trueCC1pi_muonBDT_truePion_contained->Fill(pMuonBDTScore->at(v), eventWeight);
                        if (pBacktracked_absTruePDG->at(v) == 13) trueCC1pi_muonBDT_trueMuon_contained->Fill(pMuonBDTScore->at(v), eventWeight);
                    }
                    else
                    {
                        if (pBacktracked_absTruePDG->at(v) == 2212) trueCC1pi_muonBDT_trueProton_uncontained->Fill(pMuonBDTScore->at(v), eventWeight);
                        if (pBacktracked_absTruePDG->at(v) == 211) trueCC1pi_muonBDT_truePion_uncontained->Fill(pMuonBDTScore->at(v), eventWeight);
                        if (pBacktracked_absTruePDG->at(v) == 13) trueCC1pi_muonBDT_trueMuon_uncontained->Fill(pMuonBDTScore->at(v), eventWeight);
                    }

                    if (pBacktracked_absTruePDG->at(v) == 13)
                    {
                        muonBDTVSrange_trueMuon_trueCC1pi->Fill(pTrackLength->at(v), pMuonBDTScore->at(v), eventWeight);
                        muonBDTVSTrueRange_trueMuon_trueCC1pi->Fill(pBacktracked_TrackLength->at(v), pMuonBDTScore->at(v), eventWeight);
                    }
                    else if (pBacktracked_absTruePDG->at(v) == 211)
                    {
                        muonBDTVSrange_truePion_trueCC1pi->Fill(pTrackLength->at(v), pMuonBDTScore->at(v), eventWeight);
                        muonBDTVSTrueRange_truePion_trueCC1pi->Fill(pBacktracked_TrackLength->at(v), pMuonBDTScore->at(v), eventWeight);
                    }
                }
            }
        }


        tFile->Close();
    }

    // Area normalise
    trueCC1pi_muonBDT_trueMuon->Scale(1.0 / trueCC1pi_muonBDT_trueMuon->Integral());
    trueCC1pi_muonBDT_trueProton->Scale(1.0 / trueCC1pi_muonBDT_trueProton->Integral());
    trueCC1pi_muonBDT_truePion->Scale(1.0 / trueCC1pi_muonBDT_truePion->Integral());

    // Plot the histograms together
    
    TCanvas* c01 = new TCanvas("c01", "c01", 800, 600);    

    // Set y axis range
    double maxBinContent_muonBDT = std::max({trueCC1pi_muonBDT_trueMuon->GetMaximum(), trueCC1pi_muonBDT_trueProton->GetMaximum(), trueCC1pi_muonBDT_truePion->GetMaximum()});
    trueCC1pi_muonBDT_trueMuon->GetYaxis()->SetRangeUser(0, maxBinContent_muonBDT * 1.1);
    trueCC1pi_muonBDT_trueProton->GetYaxis()->SetRangeUser(0, maxBinContent_muonBDT * 1.1);
    trueCC1pi_muonBDT_truePion->GetYaxis()->SetRangeUser(0, maxBinContent_muonBDT * 1.1);

    trueCC1pi_muonBDT_trueMuon->SetLineColor(kBlue);
    trueCC1pi_muonBDT_trueMuon->SetLineWidth(3);
    trueCC1pi_muonBDT_trueMuon->SetStats(0);
    trueCC1pi_muonBDT_trueMuon->Draw("E hist");
    trueCC1pi_muonBDT_trueProton->SetLineColor(kRed);
    trueCC1pi_muonBDT_trueProton->SetLineWidth(3);
    trueCC1pi_muonBDT_trueProton->Draw("E hist same");
    trueCC1pi_muonBDT_truePion->SetLineColor(kGreen);
    trueCC1pi_muonBDT_truePion->SetLineWidth(3);
    trueCC1pi_muonBDT_truePion->Draw("E hist same");
    
    TLegend* legend01 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend01->SetHeader("Backtracked Particle Type:");
    legend01->AddEntry(trueCC1pi_muonBDT_trueMuon, "Muon", "l");
    legend01->AddEntry(trueCC1pi_muonBDT_trueProton, "Proton", "l");
    legend01->AddEntry(trueCC1pi_muonBDT_truePion, "Pion", "l");
    legend01->Draw();
    c01->SaveAs("plots/muonPIDStudy_trueCC1pi_muonBDT.pdf");


    // Area normalise
    trueCC1pi_muonBDT_trueMuon_contained->Scale(1.0 / trueCC1pi_muonBDT_trueMuon_contained->Integral());
    trueCC1pi_muonBDT_trueProton_contained->Scale(1.0 / trueCC1pi_muonBDT_trueProton_contained->Integral());
    trueCC1pi_muonBDT_truePion_contained->Scale(1.0 / trueCC1pi_muonBDT_truePion_contained->Integral());

    // Plot the histograms together
    TCanvas* c011 = new TCanvas("c011", "c011", 800, 600);
    double maxBinContent_muonBDT_contained = std::max({trueCC1pi_muonBDT_trueMuon_contained->GetMaximum(), trueCC1pi_muonBDT_trueProton_contained->GetMaximum(), trueCC1pi_muonBDT_truePion_contained->GetMaximum()});
    trueCC1pi_muonBDT_trueMuon_contained->GetYaxis()->SetRangeUser(0, maxBinContent_muonBDT_contained * 1.1);
    trueCC1pi_muonBDT_trueProton_contained->GetYaxis()->SetRangeUser(0, maxBinContent_muonBDT_contained * 1.1);
    trueCC1pi_muonBDT_truePion_contained->GetYaxis()->SetRangeUser(0, maxBinContent_muonBDT_contained * 1.1);

    trueCC1pi_muonBDT_trueMuon_contained->SetLineColor(kBlue);
    trueCC1pi_muonBDT_trueMuon_contained->SetLineWidth(3);
    trueCC1pi_muonBDT_trueMuon_contained->SetStats(0);
    trueCC1pi_muonBDT_trueMuon_contained->Draw("E hist");
    trueCC1pi_muonBDT_trueProton_contained->SetLineColor(kRed);
    trueCC1pi_muonBDT_trueProton_contained->SetLineWidth(3);
    trueCC1pi_muonBDT_trueProton_contained->Draw("E hist same");
    trueCC1pi_muonBDT_truePion_contained->SetLineColor(kGreen);
    trueCC1pi_muonBDT_truePion_contained->SetLineWidth(3);
    trueCC1pi_muonBDT_truePion_contained->Draw("E hist same");
    
    TLegend* legend011 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend011->SetHeader("Backtracked Particle Type:");
    legend011->AddEntry(trueCC1pi_muonBDT_trueMuon_contained, "Muon", "l");
    legend011->AddEntry(trueCC1pi_muonBDT_trueProton_contained, "Proton", "l");
    legend011->AddEntry(trueCC1pi_muonBDT_truePion_contained, "Pion", "l");
    legend011->Draw();
    c011->SaveAs("plots/muonPIDStudy_trueCC1pi_muonBDT_contained.pdf");

    // Area normalise
    trueCC1pi_muonBDT_trueMuon_uncontained->Scale(1.0 / trueCC1pi_muonBDT_trueMuon_uncontained->Integral());
    trueCC1pi_muonBDT_trueProton_uncontained->Scale(1.0 / trueCC1pi_muonBDT_trueProton_uncontained->Integral());
    trueCC1pi_muonBDT_truePion_uncontained->Scale(1.0 / trueCC1pi_muonBDT_truePion_uncontained->Integral());

    // Plot the histograms together
    TCanvas* c012 = new TCanvas("c012", "c012", 800, 600);
    double maxBinContent_muonBDT_uncontained = std::max({trueCC1pi_muonBDT_trueMuon_uncontained->GetMaximum(), trueCC1pi_muonBDT_trueProton_uncontained->GetMaximum(), trueCC1pi_muonBDT_truePion_uncontained->GetMaximum()});
    trueCC1pi_muonBDT_trueMuon_uncontained->GetYaxis()->SetRangeUser(0, maxBinContent_muonBDT_uncontained * 1.1);
    trueCC1pi_muonBDT_trueProton_uncontained->GetYaxis()->SetRangeUser(0, maxBinContent_muonBDT_uncontained * 1.1);
    trueCC1pi_muonBDT_truePion_uncontained->GetYaxis()->SetRangeUser(0, maxBinContent_muonBDT_uncontained * 1.1);

    trueCC1pi_muonBDT_trueMuon_uncontained->SetLineColor(kBlue);
    trueCC1pi_muonBDT_trueMuon_uncontained->SetLineWidth(3);
    trueCC1pi_muonBDT_trueMuon_uncontained->SetStats(0);
    trueCC1pi_muonBDT_trueMuon_uncontained->Draw("E hist");
    trueCC1pi_muonBDT_trueProton_uncontained->SetLineColor(kRed);
    trueCC1pi_muonBDT_trueProton_uncontained->SetLineWidth(3);
    trueCC1pi_muonBDT_trueProton_uncontained->Draw("E hist same");
    trueCC1pi_muonBDT_truePion_uncontained->SetLineColor(kGreen);
    trueCC1pi_muonBDT_truePion_uncontained->SetLineWidth(3);
    trueCC1pi_muonBDT_truePion_uncontained->Draw("E hist same");
    TLegend* legend012 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend012->SetHeader("Backtracked Particle Type:");
    legend012->AddEntry(trueCC1pi_muonBDT_trueMuon_uncontained, "Muon", "l");
    legend012->AddEntry(trueCC1pi_muonBDT_trueProton_uncontained, "Proton", "l");
    legend012->AddEntry(trueCC1pi_muonBDT_truePion_uncontained, "Pion", "l");
    legend012->Draw();
    c012->SaveAs("plots/muonPIDStudy_trueCC1pi_muonBDT_uncontained.pdf");


    // Area normalise
    trueCC1pi_protonBDT_trueMuon->Scale(1.0 / trueCC1pi_protonBDT_trueMuon->Integral());
    trueCC1pi_protonBDT_trueProton->Scale(1.0 / trueCC1pi_protonBDT_trueProton->Integral());
    trueCC1pi_protonBDT_truePion->Scale(1.0 / trueCC1pi_protonBDT_truePion->Integral());

    // Plot the histograms together
    TCanvas* c021 = new TCanvas("c021", "c021", 800, 600);
    double maxBinContent_protonBDT = std::max({trueCC1pi_protonBDT_trueMuon->GetMaximum(), trueCC1pi_protonBDT_trueProton->GetMaximum(), trueCC1pi_protonBDT_truePion->GetMaximum()});
    trueCC1pi_protonBDT_trueMuon->GetYaxis()->SetRangeUser(0, maxBinContent_protonBDT * 1.1);
    trueCC1pi_protonBDT_trueProton->GetYaxis()->SetRangeUser(0, maxBinContent_protonBDT * 1.1);
    trueCC1pi_protonBDT_truePion->GetYaxis()->SetRangeUser(0, maxBinContent_protonBDT * 1.1);

    trueCC1pi_protonBDT_trueMuon->SetLineColor(kBlue);
    trueCC1pi_protonBDT_trueMuon->SetLineWidth(3);
    trueCC1pi_protonBDT_trueMuon->SetStats(0);
    trueCC1pi_protonBDT_trueMuon->Draw("E hist");
    trueCC1pi_protonBDT_trueProton->SetLineColor(kRed);
    trueCC1pi_protonBDT_trueProton->SetLineWidth(3);
    trueCC1pi_protonBDT_trueProton->Draw("E hist same");
    trueCC1pi_protonBDT_truePion->SetLineColor(kGreen);
    trueCC1pi_protonBDT_truePion->SetLineWidth(3);
    trueCC1pi_protonBDT_truePion->Draw("E hist same");
    TLegend* legend021 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend021->SetHeader("Backtracked Particle Type:");
    legend021->AddEntry(trueCC1pi_protonBDT_trueMuon, "Muon", "l");
    legend021->AddEntry(trueCC1pi_protonBDT_trueProton, "Proton", "l");
    legend021->AddEntry(trueCC1pi_protonBDT_truePion, "Pion", "l");
    legend021->Draw();
    c021->SaveAs("plots/muonPIDStudy_trueCC1pi_protonBDT.pdf");    


    // Area normalise
    trueCC1pi_protonBDT_trueMuon_contained->Scale(1.0 / trueCC1pi_protonBDT_trueMuon_contained->Integral());
    trueCC1pi_protonBDT_trueProton_contained->Scale(1.0 / trueCC1pi_protonBDT_trueProton_contained->Integral());
    trueCC1pi_protonBDT_truePion_contained->Scale(1.0 / trueCC1pi_protonBDT_truePion_contained->Integral());

    // Plot the histograms together
    TCanvas* c02 = new TCanvas("c02", "c02", 800, 600);
    double maxBinContent_protonBDT_contained = std::max({trueCC1pi_protonBDT_trueMuon_contained->GetMaximum(), trueCC1pi_protonBDT_trueProton_contained->GetMaximum(), trueCC1pi_protonBDT_truePion_contained->GetMaximum()});
    trueCC1pi_protonBDT_trueMuon_contained->GetYaxis()->SetRangeUser(0, maxBinContent_protonBDT_contained * 1.1);
    trueCC1pi_protonBDT_trueProton_contained->GetYaxis()->SetRangeUser(0, maxBinContent_protonBDT_contained * 1.1);
    trueCC1pi_protonBDT_truePion_contained->GetYaxis()->SetRangeUser(0, maxBinContent_protonBDT_contained * 1.1);

    trueCC1pi_protonBDT_trueMuon_contained->SetLineColor(kBlue);
    trueCC1pi_protonBDT_trueMuon_contained->SetLineWidth(3);
    trueCC1pi_protonBDT_trueMuon_contained->SetStats(0);
    trueCC1pi_protonBDT_trueMuon_contained->Draw("E hist");
    trueCC1pi_protonBDT_trueProton_contained->SetLineColor(kRed);
    trueCC1pi_protonBDT_trueProton_contained->SetLineWidth(3);
    trueCC1pi_protonBDT_trueProton_contained->Draw("E hist same");
    trueCC1pi_protonBDT_truePion_contained->SetLineColor(kGreen);
    trueCC1pi_protonBDT_truePion_contained->SetLineWidth(3);
    trueCC1pi_protonBDT_truePion_contained->Draw("E hist same");
    TLegend* legend02 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend02->SetHeader("Backtracked Particle Type:");
    legend02->AddEntry(trueCC1pi_protonBDT_trueMuon_contained, "Muon", "l");
    legend02->AddEntry(trueCC1pi_protonBDT_trueProton_contained, "Proton", "l");
    legend02->AddEntry(trueCC1pi_protonBDT_truePion_contained, "Pion", "l");
    legend02->Draw();
    c02->SaveAs("plots/muonPIDStudy_trueCC1pi_protonBDT_contained.pdf");



    // MakePlot(pdg_backtrackedMuon, "pdg_backtrackedMuon");
    // MakePlot(trueCC1pi_muonBDT_trueMuon, "trueCC1pi_muonBDT_trueMuon");

    MakePlot2D(muonMomVSpionMom_selectedTrueCC1pi, "muonMomVSpionMom_selectedTrueCC1pi");
    // MakePlot2D(muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon, "muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon");
    // Calculate ratio
    TH2D* muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon_ratio = (TH2D*)muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon->Clone("muonMomVSpionMom_trueCC1pi_correctRecoMuon_ratio");
    muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon_ratio->Divide(muonMomVSpionMom_selectedTrueCC1pi);
    MakePlot2D(muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon_ratio, "muonMomVSpionMom_selectedTrueCC1pi_correctRecoMuon_ratio", true);

    // Calculate ratio
    TH2D* muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio = (TH2D*)muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon->Clone("muonContainedVSMuonTrackLargerThanPionTrack_trueCC1pi_correctRecoMuon_ratio");
    muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio->Divide(muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi);
    MakePlot2D(muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio, "muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio", false, false);

    MakePlot2D(muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi, "muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi", false, false);

    // Print out the values of muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio
    std::cout<<"muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio:"<<std::endl;
    for (Int_t i = 1; i <= muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio->GetNbinsX(); ++i) {
        for (Int_t j = 1; j <= muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio->GetNbinsY(); ++j) {
            Double_t binContent = muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio->GetBinContent(i, j);
            std::cout << "Bin (" << i << ", " << j << ") content: " << binContent << std::endl;
        }
    }

    std::cout<<"muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi:"<<std::endl;
    // Print out the values of muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi_correctRecoMuon_ratio
    for (Int_t i = 1; i <= muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi->GetNbinsX(); ++i) {
        for (Int_t j = 1; j <= muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi->GetNbinsY(); ++j) {
            Double_t binContent = muonContainedVSMuonTrackLargerThanPionTrack_selectedTrueCC1pi->GetBinContent(i, j);
            std::cout << "Bin (" << i << ", " << j << ") content: " << binContent << std::endl;
        }
    }
    


    MakePlot2D(muonBDTVSrange_trueMuon_trueCC1pi, "muonBDTVSrange_trueMuon_trueCC1pi");
    MakePlot2D(muonBDTVSTrueRange_trueMuon_trueCC1pi, "muonBDTVSTrueRange_trueMuon_trueCC1pi");

    MakePlot2D(muonBDTVSrange_truePion_trueCC1pi, "muonBDTVSrange_truePion_trueCC1pi");
    MakePlot2D(muonBDTVSTrueRange_truePion_trueCC1pi, "muonBDTVSTrueRange_truePion_trueCC1pi");

    // MakePlot(histograms.at(names.at(i)), cc1pi_histograms.at(names.at(i)), names.at(i));
}