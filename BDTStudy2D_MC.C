#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

void BDTStudyBacktrackedProperties(const std::vector<std::tuple<std::string, std::string, std::string, float>>& files, 
                                   const std::map<std::string, std::tuple<int, double, double, bool>>& binningInfo, 
                                   const std::vector<std::string>& bdts)
{
    // Define 2D histograms for backtracked properties
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_Phi;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_CosTheta;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_Momentum;

    for (const auto& bdt : bdts)
    {
        const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
        
        histograms2D_Phi[bdt] = std::make_unique<TH2D>(
            ("h2D_Phi_" + bdt).c_str(),
            ("BDT Score vs. Backtracked Phi for " + bdt + "; " + bdt + " BDT Score; Backtracked #phi [rad]").c_str(),
            nBinsBDT, xMin, xMax, 50, -3.15, 3.15);

        histograms2D_CosTheta[bdt] = std::make_unique<TH2D>(
            ("h2D_CosTheta_" + bdt).c_str(),
            ("BDT Score vs. Backtracked cos(#theta) for " + bdt + "; " + bdt + " BDT Score; Backtracked cos(#theta)").c_str(),
            nBinsBDT, xMin, xMax, 50, -1.0, 1.0);

        histograms2D_Momentum[bdt] = std::make_unique<TH2D>(
            ("h2D_Momentum_" + bdt).c_str(),
            ("BDT Score vs. Backtracked Momentum for " + bdt + "; " + bdt + " BDT Score; Backtracked Momentum [GeV/c]").c_str(),
            nBinsBDT, xMin, xMax, 50, 0, 2.0);
    }

    // Loop over the files
    for (const auto& [sampleType, run, filePath, runWeight] : files)
    {
        // Only process MC files
        if (sampleType != "nu mc" && sampleType != "dirt mc") {
            continue;
        }

        TFile* tFile = TFile::Open(filePath.c_str());
        if (!tFile || tFile->IsZombie()) {
            std::cerr << "Error: Unable to open file: " << filePath << std::endl;
            continue;
        }

        TTree* tree = (TTree*)tFile->Get("stv_tree");
        if (!tree) {
            std::cerr << "Error: Unable to get tree from file: " << filePath << std::endl;
            tFile->Close();
            continue;
        }

        std::vector<float>* pReco_particle_ccinc_muonBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_protonBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_goldenPionBDTScore = nullptr;
        std::vector<int>* pReco_particle_ccinc_generation = nullptr;
        std::vector<bool>* pReco_particle_ccinc_isContained = nullptr;
        std::vector<float>* pReco_particle_ccinc_backtracked_phi = nullptr;
        std::vector<float>* pReco_particle_ccinc_backtracked_cosTheta = nullptr;
        std::vector<float>* pReco_particle_ccinc_backtracked_momentum = nullptr;

        tree->SetBranchAddress("reco_particle_ccinc_muonBDTScore", &pReco_particle_ccinc_muonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_protonBDTScore", &pReco_particle_ccinc_protonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_goldenPionBDTScore", &pReco_particle_ccinc_goldenPionBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_generation", &pReco_particle_ccinc_generation);
        tree->SetBranchAddress("reco_particle_ccinc_contained", &pReco_particle_ccinc_isContained);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_phi", &pReco_particle_ccinc_backtracked_phi);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_cosTheta", &pReco_particle_ccinc_backtracked_cosTheta);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_momentum", &pReco_particle_ccinc_backtracked_momentum);

        Float_t spline_weight, tuned_cv_weight;
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);

        const auto nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++)
        {
            if (i % (nEntries / 100) == 0)
            {
                const auto progress = static_cast<int>(100.0 * i / nEntries);
                std::cout << "\r"<< sampleType << " run " << run << "[" << std::string(progress, '|') << std::string(100 - progress, ' ') << "] " << progress << "%" << std::flush;
            }

            tree->GetEntry(i);
            
            // Apply event selection criteria
            if (!(tree->GetLeaf("passed_max1Uncontained")->GetValue() && 
                  tree->GetLeaf("event_cutValue_topologicalScore")->GetValue() >= 0.67 && 
                  tree->GetLeaf("event_cutValue_maxVertexDist")->GetValue() <= 9.5*9.5 && 
                  tree->GetLeaf("event_cutValue_maxVertexDist")->GetValue() >= 0))
            {
                continue;
            }

            for (size_t v = 0; v < pReco_particle_ccinc_generation->size(); v++)
            {
                if(pReco_particle_ccinc_generation->at(v) != 2 || !pReco_particle_ccinc_isContained->at(v)) 
                    continue; // Only consider generation 2 particles that are contained

                const auto phi = pReco_particle_ccinc_backtracked_phi->at(v);
                const auto cosTheta = pReco_particle_ccinc_backtracked_cosTheta->at(v);
                const auto momentum = pReco_particle_ccinc_backtracked_momentum->at(v);

                // Check if backtracked values are valid
                if (!std::isfinite(phi) || !std::isfinite(cosTheta) || !std::isfinite(momentum))
                    continue;

                for (const auto& bdt : bdts)
                {
                    const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
                    const auto bdtScore = (bdt == "muonBDTScore") ? pReco_particle_ccinc_muonBDTScore->at(v) :
                                          (bdt == "protonBDTScore") ? pReco_particle_ccinc_protonBDTScore->at(v) :
                                          pReco_particle_ccinc_goldenPionBDTScore->at(v);

                    // Ensure BDT scores are within histogram dimensions
                    if (bdtScore < xMin || bdtScore > xMax)
                        continue;

                    auto eventWeight = std::isfinite(spline_weight * tuned_cv_weight) && 
                                       spline_weight * tuned_cv_weight >= 0 && 
                                       spline_weight * tuned_cv_weight <= 30 ? 
                                       spline_weight * tuned_cv_weight : 1;
                    eventWeight *= runWeight;

                    histograms2D_Phi[bdt]->Fill(bdtScore, phi, eventWeight);
                    histograms2D_CosTheta[bdt]->Fill(bdtScore, cosTheta, eventWeight);
                    histograms2D_Momentum[bdt]->Fill(bdtScore, momentum, eventWeight);
                }
            }
        }

        tFile->Close();
    }

    // Find maximum values for color scaling
    double maxPhiValue = 0.0;
    double maxCosThetaValue = 0.0;
    double maxMomentumValue = 0.0;

    for (const auto& bdt : bdts)
    {
        maxPhiValue = std::max(maxPhiValue, histograms2D_Phi[bdt]->GetMaximum());
        maxCosThetaValue = std::max(maxCosThetaValue, histograms2D_CosTheta[bdt]->GetMaximum());
        maxMomentumValue = std::max(maxMomentumValue, histograms2D_Momentum[bdt]->GetMaximum());
    }

    // Save the 2D histograms
    for (const auto& bdt : bdts)
    {
        // Save Phi histograms
        auto cPhi = new TCanvas(("c2D_Phi_" + bdt).c_str(), "", 800, 600);
        histograms2D_Phi[bdt]->SetStats(0);
        histograms2D_Phi[bdt]->SetMinimum(0);
        histograms2D_Phi[bdt]->SetMaximum(maxPhiValue);
        histograms2D_Phi[bdt]->Draw("COLZ");
        cPhi->SaveAs(("plots/BDTStudy2D_MC_Phi_" + bdt + ".pdf").c_str());
        cPhi->SaveAs(("plots/BDTStudy2D_MC_Phi_" + bdt + ".C").c_str());
        delete cPhi;

        // Save CosTheta histograms
        auto cCosTheta = new TCanvas(("c2D_CosTheta_" + bdt).c_str(), "", 800, 600);
        histograms2D_CosTheta[bdt]->SetStats(0);
        histograms2D_CosTheta[bdt]->SetMinimum(0);
        histograms2D_CosTheta[bdt]->SetMaximum(maxCosThetaValue);
        histograms2D_CosTheta[bdt]->Draw("COLZ");
        cCosTheta->SaveAs(("plots/BDTStudy2D_MC_CosTheta_" + bdt + ".pdf").c_str());
        cCosTheta->SaveAs(("plots/BDTStudy2D_MC_CosTheta_" + bdt + ".C").c_str());
        delete cCosTheta;

        // Save Momentum histograms
        auto cMomentum = new TCanvas(("c2D_Momentum_" + bdt).c_str(), "", 800, 600);
        histograms2D_Momentum[bdt]->SetStats(0);
        histograms2D_Momentum[bdt]->SetMinimum(0);
        histograms2D_Momentum[bdt]->SetMaximum(maxMomentumValue);
        histograms2D_Momentum[bdt]->Draw("COLZ");
        cMomentum->SaveAs(("plots/BDTStudy2D_MC_Momentum_" + bdt + ".pdf").c_str());
        cMomentum->SaveAs(("plots/BDTStudy2D_MC_Momentum_" + bdt + ".C").c_str());
        delete cMomentum;
    }
}

void BDTStudy2D_MC()
{
    const std::string rootPathWithLengths = "/exp/uboone/data/users/jdetje/ubcc1piPelee/14May25_withRecoLengthAndUniversalVertDist/";

    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        std::make_tuple("nu mc", "5", rootPathWithLengths + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196 * 2.0),
        std::make_tuple("dirt mc", "5", rootPathWithLengths + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root", 0.41280),
    };

    const std::map<std::string, std::tuple<int, double, double, bool>> binningInfo = {
        {"muonBDTScore", std::make_tuple(32, -0.95, 0.65, false)},
        {"protonBDTScore", std::make_tuple(29, -0.85, 0.60, false)},
        {"goldenPionBDTScore", std::make_tuple(28, -0.90, 0.50, false)},
    };

    const std::vector<std::string> bdts = {"muonBDTScore", "protonBDTScore", "goldenPionBDTScore"};

    BDTStudyBacktrackedProperties(files, binningInfo, bdts);
}
