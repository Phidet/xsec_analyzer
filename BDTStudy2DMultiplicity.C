#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

void BDTStudy2DMultiplicity(const std::vector<std::tuple<std::string, std::string, std::string, float>>& files, const std::map<std::string, std::tuple<int, double, double, bool>>& binningInfo, const std::vector<std::string>& bdts, bool useLongestTracks = false, unsigned int nLongestTracks = 3)
{
    // Define 2D histograms for run 0
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_MC;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_BeamOn;

    // Define additional 2D histograms for reco_particle_ccinc_length
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_Length_MC;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_Length_BeamOn;

    // Define additional 2D histograms for reco_particle_ccinc_truncMeandEdx
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_dEdx_MC;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_dEdx_BeamOn;

    // Define additional 2D histograms for new variables
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_logBragg_pToMIP_MC;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_logBragg_pToMIP_BeamOn;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_logBragg_piToMIP_MC;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_logBragg_piToMIP_BeamOn;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_wiggliness_MC;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_wiggliness_BeamOn;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_trackScore_MC;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_trackScore_BeamOn;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_nDescendents_MC;
    std::map<std::string, std::unique_ptr<TH2D>> histograms2D_nDescendents_BeamOn;

    for (const auto& bdt : bdts)
    {
        const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
        histograms2D_MC[bdt] = std::make_unique<TH2D>(
            ("h2D_MC_" + bdt).c_str(),
            ("BDT Score vs. Cont. Prim. Particles (MC+Beam Off+Dirt) for " + bdt + "; " + bdt + "; Cont. Prim. Particles").c_str(),
            nBinsBDT, xMin, xMax, 10, 0, 10);

        histograms2D_BeamOn[bdt] = std::make_unique<TH2D>(
            ("h2D_BeamOn_" + bdt).c_str(),
            ("BDT Score vs. Cont. Prim. Particles (Beam On) for " + bdt + "; " + bdt + "; Cont. Prim. Particles").c_str(),
            nBinsBDT, xMin, xMax, 10, 0, 10);

        histograms2D_Length_MC[bdt] = std::make_unique<TH2D>(
            ("h2D_Length_MC_" + bdt).c_str(),
            ("BDT Score vs. Reco Length (MC+Beam Off+Dirt) for " + bdt + "; " + bdt + "; Reco Length [cm]").c_str(),
            nBinsBDT, xMin, xMax, 50, 0, 50);

        histograms2D_Length_BeamOn[bdt] = std::make_unique<TH2D>(
            ("h2D_Length_BeamOn_" + bdt).c_str(),
            ("BDT Score vs. Reco Length (Beam On) for " + bdt + "; " + bdt + "; Reco Length [cm]").c_str(),
            nBinsBDT, xMin, xMax, 50, 0, 50);

        histograms2D_dEdx_MC[bdt] = std::make_unique<TH2D>(
            ("h2D_dEdx_MC_" + bdt).c_str(),
            ("BDT Score vs. Trunc. Mean dE/dx (MC+Beam Off+Dirt) for " + bdt + "; " + bdt + "; Trunc. Mean dE/dx [MeV/cm]").c_str(),
            // nBinsBDT, xMin, xMax, 50, 0, 15);
            nBinsBDT, xMin, xMax, 5, -std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

        histograms2D_dEdx_BeamOn[bdt] = std::make_unique<TH2D>(
            ("h2D_dEdx_BeamOn_" + bdt).c_str(),
            ("BDT Score vs. Trunc. Mean dE/dx (Beam On) for " + bdt + "; " + bdt + "; Trunc. Mean dE/dx [MeV/cm]").c_str(),
            // nBinsBDT, xMin, xMax, 50, 0, 15);
            nBinsBDT, xMin, xMax, 5, -std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

        histograms2D_logBragg_pToMIP_MC[bdt] = std::make_unique<TH2D>(
            ("h2D_logBragg_pToMIP_MC_" + bdt).c_str(),
            ("BDT Score vs. log(Bragg p/MIP) (MC+Beam Off+Dirt) for " + bdt + "; " + bdt + "; log(Bragg p/MIP)").c_str(),
            // nBinsBDT, xMin, xMax, 50, -8, 8);
            nBinsBDT, xMin, xMax, 5, -std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

        histograms2D_logBragg_pToMIP_BeamOn[bdt] = std::make_unique<TH2D>(
            ("h2D_logBragg_pToMIP_BeamOn_" + bdt).c_str(),
            ("BDT Score vs. log(Bragg p/MIP) (Beam On) for " + bdt + "; " + bdt + "; log(Bragg p/MIP)").c_str(),
            // nBinsBDT, xMin, xMax, 50, -8, 8);
            nBinsBDT, xMin, xMax, 5, -std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

        histograms2D_logBragg_piToMIP_MC[bdt] = std::make_unique<TH2D>(
            ("h2D_logBragg_piToMIP_MC_" + bdt).c_str(),
            ("BDT Score vs. log(Bragg pi/MIP) (MC+Beam Off+Dirt) for " + bdt + "; " + bdt + "; log(Bragg pi/MIP)").c_str(),
            nBinsBDT, xMin, xMax, 50, -4, 7);

        histograms2D_logBragg_piToMIP_BeamOn[bdt] = std::make_unique<TH2D>(
            ("h2D_logBragg_piToMIP_BeamOn_" + bdt).c_str(),
            ("BDT Score vs. log(Bragg pi/MIP) (Beam On) for " + bdt + "; " + bdt + "; log(Bragg pi/MIP)").c_str(),
            nBinsBDT, xMin, xMax, 50, -4, 7);

        histograms2D_wiggliness_MC[bdt] = std::make_unique<TH2D>(
            ("h2D_wiggliness_MC_" + bdt).c_str(),
            ("BDT Score vs. Wiggliness (MC+Beam Off+Dirt) for " + bdt + "; " + bdt + "; Wiggliness").c_str(),
            nBinsBDT, xMin, xMax, 50, 0.0, 0.03);

        histograms2D_wiggliness_BeamOn[bdt] = std::make_unique<TH2D>(
            ("h2D_wiggliness_BeamOn_" + bdt).c_str(),
            ("BDT Score vs. Wiggliness (Beam On) for " + bdt + "; " + bdt + "; Wiggliness").c_str(),
            nBinsBDT, xMin, xMax, 50, 0.0, 0.03);

        histograms2D_trackScore_MC[bdt] = std::make_unique<TH2D>(
            ("h2D_trackScore_MC_" + bdt).c_str(),
            ("BDT Score vs. Track Score (MC+Beam Off+Dirt) for " + bdt + "; " + bdt + "; Track Score").c_str(),
            nBinsBDT, xMin, xMax, 50, 0, 1.0);

        histograms2D_trackScore_BeamOn[bdt] = std::make_unique<TH2D>(
            ("h2D_trackScore_BeamOn_" + bdt).c_str(),
            ("BDT Score vs. Track Score (Beam On) for " + bdt + "; " + bdt + "; Track Score").c_str(),
            nBinsBDT, xMin, xMax, 50, 0, 1.0);

        histograms2D_nDescendents_MC[bdt] = std::make_unique<TH2D>(
            ("h2D_nDescendents_MC_" + bdt).c_str(),
            ("BDT Score vs. nDescendents (MC+Beam Off+Dirt) for " + bdt + "; " + bdt + "; nDescendents").c_str(),
            nBinsBDT, xMin, xMax, 5, 0, 5);

        histograms2D_nDescendents_BeamOn[bdt] = std::make_unique<TH2D>(
            ("h2D_nDescendents_BeamOn_" + bdt).c_str(),
            ("BDT Score vs. nDescendents (Beam On) for " + bdt + "; " + bdt + "; nDescendents").c_str(),
            nBinsBDT, xMin, xMax, 5, 0, 5);
    }

    // Loop over the files
    for (const auto& [sampleType, run, filePath, runWeight] : files)
    {
        const bool isBeamOn = (sampleType == "beam on");
        const bool isBeamOff = (sampleType == "beam off");
        const bool isMC = (sampleType == "nu mc" || sampleType == "dirt mc");

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
        // std::vector<float>* pReco_particle_ccinc_length = nullptr;
        std::vector<float>* pReco_particle_ccinc_truncMeandEdx = nullptr; // Added for dE/dx
        std::vector<float>* pReco_particle_ccinc_logBragg_pToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_logBragg_piToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_wiggliness = nullptr;
        std::vector<float>* pReco_particle_ccinc_trackScore = nullptr;
        std::vector<int>* pReco_particle_ccinc_nDescendents = nullptr;
        std::vector<float>* pReco_particle_ccinc_backtracked_phi = nullptr; // Only for MC
        std::vector<float>* pReco_particle_ccinc_backtracked_cosTheta = nullptr; // Only for MC
        std::vector<float>* pReco_particle_ccinc_backtracked_momentum = nullptr; // Only for MC

        tree->SetBranchAddress("reco_particle_ccinc_muonBDTScore", &pReco_particle_ccinc_muonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_protonBDTScore", &pReco_particle_ccinc_protonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_goldenPionBDTScore", &pReco_particle_ccinc_goldenPionBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_generation", &pReco_particle_ccinc_generation);
        tree->SetBranchAddress("reco_particle_ccinc_contained", &pReco_particle_ccinc_isContained);
        // tree->SetBranchAddress("reco_particle_ccinc_length", &pReco_particle_ccinc_length);
        tree->SetBranchAddress("reco_particle_ccinc_truncMeandEdx", &pReco_particle_ccinc_truncMeandEdx);
        tree->SetBranchAddress("reco_particle_ccinc_logBragg_pToMIP", &pReco_particle_ccinc_logBragg_pToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_logBragg_piToMIP", &pReco_particle_ccinc_logBragg_piToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_wiggliness", &pReco_particle_ccinc_wiggliness);
        tree->SetBranchAddress("reco_particle_ccinc_trackScore", &pReco_particle_ccinc_trackScore);
        tree->SetBranchAddress("reco_particle_ccinc_nDescendents", &pReco_particle_ccinc_nDescendents);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_phi", &pReco_particle_ccinc_backtracked_phi);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_cosTheta", &pReco_particle_ccinc_backtracked_cosTheta);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_momentum", &pReco_particle_ccinc_backtracked_momentum);

        // std::vector<float>* pEvent_cutValue_topologicalScore = nullptr;
        // tree->SetBranchAddress("event_cutValue_topologicalScore", &pEvent_cutValue_topologicalScore);

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
            
            // if (!tree->GetLeaf("passed_topologicalScoreCC")->GetValue())
            // if (!tree->GetLeaf("passed_openingAngle")->GetValue())
            // if (!(tree->GetLeaf("passed_max1Uncontained")->GetValue() && tree->GetLeaf("event_cutValue_topologicalScore")->GetValue() >= 0.67))
            // if (!(tree->GetLeaf("passed_max1Uncontained")->GetValue() && tree->GetLeaf("event_cutValue_topologicalScore")->GetValue() >= 0.67 && tree->GetLeaf("event_cutValue_nUncontained")->GetValue() == 0))
            // std::cout<<"DEBUG: event_cutValue_maxVertexDist"<< tree->GetLeaf("event_cutValue_maxVertexDist")->GetValue()<<std::endl;
            // if (!(tree->GetLeaf("passed_max1Uncontained")->GetValue() && tree->GetLeaf("event_cutValue_topologicalScore")->GetValue() >= 0.67 && tree->GetLeaf("event_cutValue_maxVertexDist")->GetValue() <= 9.5*9.5 && tree->GetLeaf("event_cutValue_maxVertexDist")->GetValue() >= 0))
            if (!(tree->GetLeaf("passed_max1Uncontained")->GetValue() && tree->GetLeaf("event_cutValue_nUncontained")->GetValue() == 0)) // Contained events for the muon BDT
            {
                continue;
            }

            unsigned int mainParticles = 0;
            std::vector<std::pair<float, size_t>> trackLengths; // Pair of (length, index)


            for (size_t v = 0; v < pReco_particle_ccinc_generation->size(); v++)
            {
                if (pReco_particle_ccinc_generation->at(v) == 2 && pReco_particle_ccinc_isContained->at(v) && pReco_particle_ccinc_truncMeandEdx->at(v) >= 0 && pReco_particle_ccinc_truncMeandEdx->at(v) < 10000 )
                {
                    mainParticles++;
                    // if (useLongestTracks)
                    // {
                    //     trackLengths.emplace_back(pReco_particle_ccinc_length->at(v), v);
                    // }
                }
            }

            // If using longest tracks, sort and select the top nLongestTracks
            std::vector<size_t> selectedIndices;
            if (useLongestTracks)
            {
                std::sort(trackLengths.rbegin(), trackLengths.rend()); // Sort in descending order
                for (size_t i = 0; i < std::min(nLongestTracks, static_cast<unsigned int>(trackLengths.size())); i++)
                {
                    selectedIndices.push_back(trackLengths[i].second);
                }
            }

            for (size_t v = 0; v < pReco_particle_ccinc_generation->size(); v++)
            {
                if(pReco_particle_ccinc_generation->at(v) != 2 || !pReco_particle_ccinc_isContained->at(v)) 
                    continue; // Only consider generation 2 particles that are contained

                if (useLongestTracks && std::find(selectedIndices.begin(), selectedIndices.end(), v) == selectedIndices.end())
                {
                    continue; // Skip if not one of the n longest tracks
                }

                // const auto recoLength = pReco_particle_ccinc_length->at(v);
                const auto recoLength = 0;
                const auto recoDEdx = pReco_particle_ccinc_truncMeandEdx->at(v);
                const auto recoLogBragg_pToMIP = pReco_particle_ccinc_logBragg_pToMIP->at(v);
                const auto recoLogBragg_piToMIP = pReco_particle_ccinc_logBragg_piToMIP->at(v);
                const auto recoWiggliness = pReco_particle_ccinc_wiggliness->at(v);
                const auto recoTrackScore = pReco_particle_ccinc_trackScore->at(v);
                const auto recoNDescendents = pReco_particle_ccinc_nDescendents->at(v);

                for (const auto& bdt : bdts)
                {
                    const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
                    const auto bdtScore = (bdt == "muonBDTScore") ? pReco_particle_ccinc_muonBDTScore->at(v) :
                                          (bdt == "protonBDTScore") ? pReco_particle_ccinc_protonBDTScore->at(v) :
                                          pReco_particle_ccinc_goldenPionBDTScore->at(v);

                    // Ensure values are within histogram dimensions
                    if (bdtScore < xMin || bdtScore > xMax || recoLength < 0 || recoLength > 50)
                        continue;

                    auto eventWeight = std::isfinite(spline_weight * tuned_cv_weight) && spline_weight * tuned_cv_weight >= 0 && spline_weight * tuned_cv_weight <= 30 ? spline_weight * tuned_cv_weight : 1;
                    eventWeight *= runWeight;

                    if (isBeamOn)
                    {
                        histograms2D_BeamOn[bdt]->Fill(bdtScore, mainParticles, eventWeight);
                        histograms2D_Length_BeamOn[bdt]->Fill(bdtScore, recoLength, eventWeight);
                        histograms2D_dEdx_BeamOn[bdt]->Fill(bdtScore, recoDEdx, eventWeight);
                        histograms2D_logBragg_pToMIP_BeamOn[bdt]->Fill(bdtScore, recoLogBragg_pToMIP, eventWeight);
                        histograms2D_logBragg_piToMIP_BeamOn[bdt]->Fill(bdtScore, recoLogBragg_piToMIP, eventWeight);
                        histograms2D_wiggliness_BeamOn[bdt]->Fill(bdtScore, recoWiggliness, eventWeight);
                        histograms2D_trackScore_BeamOn[bdt]->Fill(bdtScore, recoTrackScore, eventWeight);
                        histograms2D_nDescendents_BeamOn[bdt]->Fill(bdtScore, recoNDescendents, eventWeight);
                    }
                    else if (isMC || isBeamOff)
                    {
                        histograms2D_MC[bdt]->Fill(bdtScore, mainParticles, eventWeight);
                        histograms2D_Length_MC[bdt]->Fill(bdtScore, recoLength, eventWeight);
                        histograms2D_dEdx_MC[bdt]->Fill(bdtScore, recoDEdx, eventWeight);
                        histograms2D_logBragg_pToMIP_MC[bdt]->Fill(bdtScore, recoLogBragg_pToMIP, eventWeight);
                        histograms2D_logBragg_piToMIP_MC[bdt]->Fill(bdtScore, recoLogBragg_piToMIP, eventWeight);
                        histograms2D_wiggliness_MC[bdt]->Fill(bdtScore, recoWiggliness, eventWeight);
                        histograms2D_trackScore_MC[bdt]->Fill(bdtScore, recoTrackScore, eventWeight);
                        histograms2D_nDescendents_MC[bdt]->Fill(bdtScore, recoNDescendents, eventWeight);
                    }
                }
            }
        }

        tFile->Close();
    }

    // Determine the maximum histogram cell value across all histograms
    double maxCellValue = 0.0;
    for (const auto& bdt : bdts)
    {
        maxCellValue = std::max(maxCellValue, histograms2D_MC[bdt]->GetMaximum());
        maxCellValue = std::max(maxCellValue, histograms2D_BeamOn[bdt]->GetMaximum());
    }

    // Determine the maximum histogram cell value for reco length histograms
    double maxCellValueLength = 0.0;
    for (const auto& bdt : bdts)
    {
        maxCellValueLength = std::max(maxCellValueLength, histograms2D_Length_MC[bdt]->GetMaximum());
        maxCellValueLength = std::max(maxCellValueLength, histograms2D_Length_BeamOn[bdt]->GetMaximum());
    }

    // Determine the maximum histogram cell value for dE/dx histograms
    double maxCellValueDEdx = 0.0;
    // Determine maximum cell values for new variables
    double maxCellValueLogBragg_pToMIP = 0.0;
    double maxCellValueLogBragg_piToMIP = 0.0;
    double maxCellValueWiggliness = 0.0;
    double maxCellValueTrackScore = 0.0;
    double maxCellValueNDescendents = 0.0;
    
    for (const auto& bdt : bdts)
    {
        maxCellValueDEdx = std::max(maxCellValueDEdx, histograms2D_dEdx_MC[bdt]->GetMaximum());
        maxCellValueDEdx = std::max(maxCellValueDEdx, histograms2D_dEdx_BeamOn[bdt]->GetMaximum());

        maxCellValueLogBragg_pToMIP = std::max(maxCellValueLogBragg_pToMIP, histograms2D_logBragg_pToMIP_MC[bdt]->GetMaximum());
        maxCellValueLogBragg_pToMIP = std::max(maxCellValueLogBragg_pToMIP, histograms2D_logBragg_pToMIP_BeamOn[bdt]->GetMaximum());
        
        maxCellValueLogBragg_piToMIP = std::max(maxCellValueLogBragg_piToMIP, histograms2D_logBragg_piToMIP_MC[bdt]->GetMaximum());
        maxCellValueLogBragg_piToMIP = std::max(maxCellValueLogBragg_piToMIP, histograms2D_logBragg_piToMIP_BeamOn[bdt]->GetMaximum());
        
        maxCellValueWiggliness = std::max(maxCellValueWiggliness, histograms2D_wiggliness_MC[bdt]->GetMaximum());
        maxCellValueWiggliness = std::max(maxCellValueWiggliness, histograms2D_wiggliness_BeamOn[bdt]->GetMaximum());
        
        maxCellValueTrackScore = std::max(maxCellValueTrackScore, histograms2D_trackScore_MC[bdt]->GetMaximum());
        maxCellValueTrackScore = std::max(maxCellValueTrackScore, histograms2D_trackScore_BeamOn[bdt]->GetMaximum());
        
        maxCellValueNDescendents = std::max(maxCellValueNDescendents, histograms2D_nDescendents_MC[bdt]->GetMaximum());
        maxCellValueNDescendents = std::max(maxCellValueNDescendents, histograms2D_nDescendents_BeamOn[bdt]->GetMaximum());
    }

    gStyle->SetNumberContours(100); // Increase the number of color contours for finer color granularity

    // Save the 2D histograms with consistent z-axis (color) range
    for (const auto& bdt : bdts)
    {
        auto c2D_MC = new TCanvas(("c2D_MC_" + bdt).c_str(), "", 800, 600);
        histograms2D_MC[bdt]->SetStats(0);
        histograms2D_MC[bdt]->SetMinimum(0);           // Set minimum z value (color scale)
        histograms2D_MC[bdt]->SetMaximum(maxCellValue); // Set maximum z value (color scale)
        histograms2D_MC[bdt]->Draw("COLZ");
        c2D_MC->SaveAs(("plots/BDTStudy2D_MC_run0_" + bdt + ".pdf").c_str());
        c2D_MC->SaveAs(("plots/BDTStudy2D_MC_run0_" + bdt + ".C").c_str());

        auto c2D_BeamOn = new TCanvas(("c2D_BeamOn_" + bdt).c_str(), "", 800, 600);
        histograms2D_BeamOn[bdt]->SetStats(0);
        histograms2D_BeamOn[bdt]->SetMinimum(0);              // Set minimum z value
        histograms2D_BeamOn[bdt]->SetMaximum(maxCellValue);   // Set maximum z value
        histograms2D_BeamOn[bdt]->Draw("COLZ");
        c2D_BeamOn->SaveAs(("plots/BDTStudy2D_BeamOn_run0_" + bdt + ".pdf").c_str());
        c2D_BeamOn->SaveAs(("plots/BDTStudy2D_BeamOn_run0_" + bdt + ".C").c_str());

        auto c2D_Length_MC = new TCanvas(("c2D_Length_MC_" + bdt).c_str(), "", 800, 600);
        histograms2D_Length_MC[bdt]->SetStats(0);
        histograms2D_Length_MC[bdt]->SetMinimum(0);                  // Set minimum z value
        histograms2D_Length_MC[bdt]->SetMaximum(maxCellValueLength); // Set maximum z value
        histograms2D_Length_MC[bdt]->Draw("COLZ");
        c2D_Length_MC->SaveAs(("plots/BDTStudy2D_Length_MC_run0_" + bdt + ".pdf").c_str());
        c2D_Length_MC->SaveAs(("plots/BDTStudy2D_Length_MC_run0_" + bdt + ".C").c_str());

        auto c2D_Length_BeamOn = new TCanvas(("c2D_Length_BeamOn_" + bdt).c_str(), "", 800, 600);
        histograms2D_Length_BeamOn[bdt]->SetStats(0);
        histograms2D_Length_BeamOn[bdt]->SetMinimum(0);                   // Set minimum z value
        histograms2D_Length_BeamOn[bdt]->SetMaximum(maxCellValueLength);  // Set maximum z value
        histograms2D_Length_BeamOn[bdt]->Draw("COLZ");
        c2D_Length_BeamOn->SaveAs(("plots/BDTStudy2D_Length_BeamOn_run0_" + bdt + ".pdf").c_str());
        c2D_Length_BeamOn->SaveAs(("plots/BDTStudy2D_Length_BeamOn_run0_" + bdt + ".C").c_str());

        // Save dE/dx 2D histograms
        auto c2D_dEdx_MC = new TCanvas(("c2D_dEdx_MC_" + bdt).c_str(), "", 800, 600);
        histograms2D_dEdx_MC[bdt]->SetStats(0);
        histograms2D_dEdx_MC[bdt]->SetMinimum(0);                  // Set minimum z value
        histograms2D_dEdx_MC[bdt]->SetMaximum(maxCellValueDEdx);   // Set maximum z value
        histograms2D_dEdx_MC[bdt]->Draw("COLZ");
        c2D_dEdx_MC->SaveAs(("plots/BDTStudy2D_dEdx_MC_run0_" + bdt + ".pdf").c_str());
        c2D_dEdx_MC->SaveAs(("plots/BDTStudy2D_dEdx_MC_run0_" + bdt + ".C").c_str());

        auto c2D_dEdx_BeamOn = new TCanvas(("c2D_dEdx_BeamOn_" + bdt).c_str(), "", 800, 600);
        histograms2D_dEdx_BeamOn[bdt]->SetStats(0);
        histograms2D_dEdx_BeamOn[bdt]->SetMinimum(0);                   // Set minimum z value
        histograms2D_dEdx_BeamOn[bdt]->SetMaximum(maxCellValueDEdx);    // Set maximum z value
        histograms2D_dEdx_BeamOn[bdt]->Draw("COLZ");
        c2D_dEdx_BeamOn->SaveAs(("plots/BDTStudy2D_dEdx_BeamOn_run0_" + bdt + ".pdf").c_str());
        c2D_dEdx_BeamOn->SaveAs(("plots/BDTStudy2D_dEdx_BeamOn_run0_" + bdt + ".C").c_str());

        // Save logBragg_pToMIP 2D histograms
        auto c2D_logBragg_pToMIP_MC = new TCanvas(("c2D_logBragg_pToMIP_MC_" + bdt).c_str(), "", 800, 600);
        histograms2D_logBragg_pToMIP_MC[bdt]->SetStats(0);
        histograms2D_logBragg_pToMIP_MC[bdt]->SetMinimum(0);
        histograms2D_logBragg_pToMIP_MC[bdt]->SetMaximum(maxCellValueLogBragg_pToMIP);
        histograms2D_logBragg_pToMIP_MC[bdt]->Draw("COLZ");
        c2D_logBragg_pToMIP_MC->SaveAs(("plots/BDTStudy2D_logBragg_pToMIP_MC_run0_" + bdt + ".pdf").c_str());
        c2D_logBragg_pToMIP_MC->SaveAs(("plots/BDTStudy2D_logBragg_pToMIP_MC_run0_" + bdt + ".C").c_str());

        auto c2D_logBragg_pToMIP_BeamOn = new TCanvas(("c2D_logBragg_pToMIP_BeamOn_" + bdt).c_str(), "", 800, 600);
        histograms2D_logBragg_pToMIP_BeamOn[bdt]->SetStats(0);
        histograms2D_logBragg_pToMIP_BeamOn[bdt]->SetMinimum(0);
        histograms2D_logBragg_pToMIP_BeamOn[bdt]->SetMaximum(maxCellValueLogBragg_pToMIP);
        histograms2D_logBragg_pToMIP_BeamOn[bdt]->Draw("COLZ");
        c2D_logBragg_pToMIP_BeamOn->SaveAs(("plots/BDTStudy2D_logBragg_pToMIP_BeamOn_run0_" + bdt + ".pdf").c_str());
        c2D_logBragg_pToMIP_BeamOn->SaveAs(("plots/BDTStudy2D_logBragg_pToMIP_BeamOn_run0_" + bdt + ".C").c_str());

        // Save logBragg_piToMIP 2D histograms
        auto c2D_logBragg_piToMIP_MC = new TCanvas(("c2D_logBragg_piToMIP_MC_" + bdt).c_str(), "", 800, 600);
        histograms2D_logBragg_piToMIP_MC[bdt]->SetStats(0);
        histograms2D_logBragg_piToMIP_MC[bdt]->SetMinimum(0);
        histograms2D_logBragg_piToMIP_MC[bdt]->SetMaximum(maxCellValueLogBragg_piToMIP);
        histograms2D_logBragg_piToMIP_MC[bdt]->Draw("COLZ");
        c2D_logBragg_piToMIP_MC->SaveAs(("plots/BDTStudy2D_logBragg_piToMIP_MC_run0_" + bdt + ".pdf").c_str());
        c2D_logBragg_piToMIP_MC->SaveAs(("plots/BDTStudy2D_logBragg_piToMIP_MC_run0_" + bdt + ".C").c_str());

        auto c2D_logBragg_piToMIP_BeamOn = new TCanvas(("c2D_logBragg_piToMIP_BeamOn_" + bdt).c_str(), "", 800, 600);
        histograms2D_logBragg_piToMIP_BeamOn[bdt]->SetStats(0);
        histograms2D_logBragg_piToMIP_BeamOn[bdt]->SetMinimum(0);
        histograms2D_logBragg_piToMIP_BeamOn[bdt]->SetMaximum(maxCellValueLogBragg_piToMIP);
        histograms2D_logBragg_piToMIP_BeamOn[bdt]->Draw("COLZ");
        c2D_logBragg_piToMIP_BeamOn->SaveAs(("plots/BDTStudy2D_logBragg_piToMIP_BeamOn_run0_" + bdt + ".pdf").c_str());
        c2D_logBragg_piToMIP_BeamOn->SaveAs(("plots/BDTStudy2D_logBragg_piToMIP_BeamOn_run0_" + bdt + ".C").c_str());

        // Save wiggliness 2D histograms
        auto c2D_wiggliness_MC = new TCanvas(("c2D_wiggliness_MC_" + bdt).c_str(), "", 800, 600);
        histograms2D_wiggliness_MC[bdt]->SetStats(0);
        histograms2D_wiggliness_MC[bdt]->SetMinimum(0);
        histograms2D_wiggliness_MC[bdt]->SetMaximum(maxCellValueWiggliness);
        histograms2D_wiggliness_MC[bdt]->Draw("COLZ");
        c2D_wiggliness_MC->SaveAs(("plots/BDTStudy2D_wiggliness_MC_run0_" + bdt + ".pdf").c_str());
        c2D_wiggliness_MC->SaveAs(("plots/BDTStudy2D_wiggliness_MC_run0_" + bdt + ".C").c_str());

        auto c2D_wiggliness_BeamOn = new TCanvas(("c2D_wiggliness_BeamOn_" + bdt).c_str(), "", 800, 600);
        histograms2D_wiggliness_BeamOn[bdt]->SetStats(0);
        histograms2D_wiggliness_BeamOn[bdt]->SetMinimum(0);
        histograms2D_wiggliness_BeamOn[bdt]->SetMaximum(maxCellValueWiggliness);
        histograms2D_wiggliness_BeamOn[bdt]->Draw("COLZ");
        c2D_wiggliness_BeamOn->SaveAs(("plots/BDTStudy2D_wiggliness_BeamOn_run0_" + bdt + ".pdf").c_str());
        c2D_wiggliness_BeamOn->SaveAs(("plots/BDTStudy2D_wiggliness_BeamOn_run0_" + bdt + ".C").c_str());

        // Save trackScore 2D histograms
        auto c2D_trackScore_MC = new TCanvas(("c2D_trackScore_MC_" + bdt).c_str(), "", 800, 600);
        histograms2D_trackScore_MC[bdt]->SetStats(0);
        histograms2D_trackScore_MC[bdt]->SetMinimum(0);
        histograms2D_trackScore_MC[bdt]->SetMaximum(maxCellValueTrackScore);
        histograms2D_trackScore_MC[bdt]->Draw("COLZ");
        c2D_trackScore_MC->SaveAs(("plots/BDTStudy2D_trackScore_MC_run0_" + bdt + ".pdf").c_str());
        c2D_trackScore_MC->SaveAs(("plots/BDTStudy2D_trackScore_MC_run0_" + bdt + ".C").c_str());

        auto c2D_trackScore_BeamOn = new TCanvas(("c2D_trackScore_BeamOn_" + bdt).c_str(), "", 800, 600);
        histograms2D_trackScore_BeamOn[bdt]->SetStats(0);
        histograms2D_trackScore_BeamOn[bdt]->SetMinimum(0);
        histograms2D_trackScore_BeamOn[bdt]->SetMaximum(maxCellValueTrackScore);
        histograms2D_trackScore_BeamOn[bdt]->Draw("COLZ");
        c2D_trackScore_BeamOn->SaveAs(("plots/BDTStudy2D_trackScore_BeamOn_run0_" + bdt + ".pdf").c_str());
        c2D_trackScore_BeamOn->SaveAs(("plots/BDTStudy2D_trackScore_BeamOn_run0_" + bdt + ".C").c_str());

        // Save nDescendents 2D histograms
        auto c2D_nDescendents_MC = new TCanvas(("c2D_nDescendents_MC_" + bdt).c_str(), "", 800, 600);
        histograms2D_nDescendents_MC[bdt]->SetStats(0);
        histograms2D_nDescendents_MC[bdt]->SetMinimum(0);
        histograms2D_nDescendents_MC[bdt]->SetMaximum(maxCellValueNDescendents);
        histograms2D_nDescendents_MC[bdt]->Draw("COLZ");
        c2D_nDescendents_MC->SaveAs(("plots/BDTStudy2D_nDescendents_MC_run0_" + bdt + ".pdf").c_str());
        c2D_nDescendents_MC->SaveAs(("plots/BDTStudy2D_nDescendents_MC_run0_" + bdt + ".C").c_str());

        auto c2D_nDescendents_BeamOn = new TCanvas(("c2D_nDescendents_BeamOn_" + bdt).c_str(), "", 800, 600);
        histograms2D_nDescendents_BeamOn[bdt]->SetStats(0);
        histograms2D_nDescendents_BeamOn[bdt]->SetMinimum(0);
        histograms2D_nDescendents_BeamOn[bdt]->SetMaximum(maxCellValueNDescendents);
        histograms2D_nDescendents_BeamOn[bdt]->Draw("COLZ");
        c2D_nDescendents_BeamOn->SaveAs(("plots/BDTStudy2D_nDescendents_BeamOn_run0_" + bdt + ".pdf").c_str());
        c2D_nDescendents_BeamOn->SaveAs(("plots/BDTStudy2D_nDescendents_BeamOn_run0_" + bdt + ".C").c_str());
    }
}

void BDTStudy2DMultiplicity()
{
    // // const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/1March24/";
    // const std::string rootPathWithLengths = "/exp/uboone/data/users/jdetje/ubcc1piPelee/14May25_withRecoLengthAndUniversalVertDist/";

    // const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
    //     std::make_tuple("beam on", "5", rootPathWithLengths + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root", 1.f),
    //     std::make_tuple("beam off", "5", rootPathWithLengths + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root", 0.32807),
    //     std::make_tuple("nu mc", "5", rootPathWithLengths + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196 * 2.0),
    //     std::make_tuple("dirt mc", "5", rootPathWithLengths + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root", 0.41280),
    // };

    const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/1March24/";

    // List of files; tuple: file type, run, path, fileWeight
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        // std::make_tuple("beam on", "1",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "2",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "3",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 1.f),
        std::make_tuple("beam on", "4bcd", rootPath + "bnb_beam_on_peleeTuple_uboone_run4bcd_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "5",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root", 1.f),

        // std::make_tuple("beam off", "1",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi.root", 0.56421),
        // std::make_tuple("beam off", "2",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 0.40202),
        // std::make_tuple("beam off", "3",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 0.30657),
        std::make_tuple("beam off", "4bcd", rootPath + "bnb_beam_off_peleeTuple_uboone_run4bcd_ubcc1pi.root", 0.30175),
        // std::make_tuple("beam off", "5",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root", 0.32807),

        // std::make_tuple("nu mc", "1",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
        // std::make_tuple("nu mc", "2",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root", 0.25750*2.0),
        // std::make_tuple("nu mc", "3",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root", 0.20113*2.0),
        std::make_tuple("nu mc", "4bcd", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root", 0.13074*2.0),
        // std::make_tuple("nu mc", "5",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196*2.0),

        // std::make_tuple("dirt mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root", 0.52806),
        // std::make_tuple("dirt mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root", 0.27521),
        // std::make_tuple("dirt mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root", 0.80892),
        std::make_tuple("dirt mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi.root", 0.329301), // old value: 0.39701),
        // std::make_tuple("dirt mc", "5",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root", 0.41280),
    };


    const std::map<std::string, std::tuple<int, double, double, bool>> binningInfo = {
        {"muonBDTScore", std::make_tuple(25, -0.95, 0.65, false)},
        {"protonBDTScore", std::make_tuple(22, -0.830, 0.582, false)},
        {"goldenPionBDTScore", std::make_tuple(25, -0.705, 0.500, false)},
    };

    const std::map<std::string, std::tuple<int, double, double, bool>> binningInfoAlt = {
        {"muonBDTScore", std::make_tuple(20, -0.95, 0.65, false)},
        {"protonBDTScore", std::make_tuple(20, -0.85, 0.60, false)},
        {"goldenPionBDTScore", std::make_tuple(20, -0.90, 0.50, false)},
    };

    const std::vector<std::string> bdts = {"muonBDTScore", "protonBDTScore", "goldenPionBDTScore"};

    BDTStudy2DMultiplicity(files, binningInfo, bdts, /*useLongestTracks*/ false, /*nLongestTracks*/ 1);
}
