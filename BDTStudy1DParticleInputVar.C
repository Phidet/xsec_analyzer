#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <map>
#include <string>
#include <memory>

void BDTStudy1DParticleInputVar() 
{
    // Enabke multi-threading
    ROOT::EnableImplicitMT();

    const std::string rootPath = "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/22Feb24/";
    // tuple: type, run, file path, run weight
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        std::make_tuple("nu mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root", 0.13011),
        std::make_tuple("nu mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root", 0.25750),
        std::make_tuple("nu mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root", 0.20113),
        std::make_tuple("nu mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi.root", 0.13074),
        std::make_tuple("nu mc", "5",  rootPath + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root", 0.15196),
    };


    const std::vector<std::string> runs {"0", "1", "2", "3", "4bcd", "5"}; // Here 0 is the full set of all runs
    const std::vector<std::string> variables = {"logBragg_pToMIP", "logBragg_piToMIP", "truncMeandEdx", "wiggliness", "trackScore", "nDescendents"};
    const std::vector<std::string> particles = {"P", "Mu", "Pi", "Golden Pi", "EXT"};
    const std::vector<bool> trainingOptions = {true, false};

    // map: variable, tuple(nBins, xMin, xMax, xAxisLog, yAxisLog)
    const std::map<std::string, std::tuple<int, double, double, bool, bool>> binningInfo = {
        {"logBragg_pToMIP", std::make_tuple(60, -8, 8, false, false)},
        {"logBragg_piToMIP", std::make_tuple(60, -4, 7, false, false)},
        {"truncMeandEdx", std::make_tuple(60, 0, 10.0, false, false)},
        {"wiggliness", std::make_tuple(60, 0.00005, 0.25, true, false)},
        {"trackScore", std::make_tuple(60, 0, 1.0, false, true)},
        {"nDescendents", std::make_tuple(4, 0, 4, false, false)}
    };

    // map: run, variable, particle
    std::map<std::string, std::map<std::string, std::map<std::string, std::unique_ptr<TH1D>>>> histograms1D;
    for (const auto& run : runs) {
        for (const auto& variable : variables) {
            const auto [nBinsVar, xMin, xMax, xAxisLog, yAxisLog] = binningInfo.at(variable);
            for (const auto& particle : particles) {
                    const auto name = "h_" + run + "_" + variable + "_" + particle;
                    // histograms1D[run][variable][particle] = std::make_unique<TH1D>((name + "_" + variable + "_" + particle).c_str(), (name + " reco Particles Backtracked to True " + particle + " in Signal Events; " + variable + " BDT Score").c_str(), nBinsVar, xMin, xMax);

                    if(xAxisLog)
                    {
                        // Create a vector to hold the bin edges
                        std::vector<double> binEdges(nBinsVar + 1);
                        double logMin = TMath::Log10(xMin);
                        double logMax = TMath::Log10(xMax);
                        double binWidth = (logMax - logMin) / nBinsVar;
                        for (int i = 0; i <= nBinsVar; i++) {
                            binEdges[i] = TMath::Power(10, logMin + i * binWidth);
                        }

                        // Create the histogram with logarithmic binning
                        histograms1D[run][variable][particle] = std::make_unique<TH1D>((name + "_" + variable + "_" + particle).c_str(), (name + " reco Particles Backtracked to True " + particle + " in Signal Events; " + variable + " BDT Score").c_str(), nBinsVar, binEdges.data());
                    }
                    else
                    {
                        histograms1D[run][variable][particle] = std::make_unique<TH1D>((name + "_" + variable + "_" + particle).c_str(), (name + " reco Particles Backtracked to True " + particle + " in Signal Events; " + variable + " BDT Score").c_str(), nBinsVar, xMin, xMax);
                    }
                    
                    histograms1D[run][variable][particle]->Sumw2();
            }
        }
    }

    // Loop over the files
    for (int f = 0; f < files.size(); f++)
    {
        const auto [sampleType, run, filePath, runWeight] = files.at(f);

        const auto isBeamOn = (sampleType == "beam on");
        const auto isBeamOff = (sampleType == "beam off");
        const auto isNuMC = (sampleType == "nu mc");
        const auto isMC = (sampleType == "nu mc" || sampleType == "dirt mc");
        const auto isDirt = (sampleType == "dirt mc");

        // Open the file and get the tree
        TFile* tFile = TFile::Open(filePath.c_str());
        if (!tFile || tFile->IsZombie()) {
            std::cerr << "Error: Unable to open file: " << filePath << std::endl;
            return;
        }

        TTree* tree = (TTree*)tFile->Get("stv_tree");
        if (!tree) {
            std::cerr << "Error: Unable to get tree from file: " << filePath << std::endl;
            return;
        }

        // std::vector<float>* pReco_particle_ccinc_muonBDTScore = nullptr;
        // std::vector<float>* pReco_particle_ccinc_protonBDTScore = nullptr;
        // std::vector<float>* pReco_particle_ccinc_goldenPionBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_logBragg_pToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_logBragg_piToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_truncMeandEdx = nullptr;
        std::vector<float>* pReco_particle_ccinc_wiggliness = nullptr;
        std::vector<float>* pReco_particle_ccinc_trackScore = nullptr;
        std::vector<int>* pReco_particle_ccinc_nDescendents = nullptr;
        std::vector<int>* pReco_particle_ccinc_generation = nullptr;
        std::vector<int>* pReco_particle_ccinc_backtracked_pdg = nullptr;
        // std::vector<float>* pReco_particle_ccinc_backtracked_momentum = nullptr;
        // std::vector<float>* pReco_particle_ccinc_backtracked_cosTheta = nullptr;
        // std::vector<float>* pReco_particle_ccinc_backtracked_phi = nullptr;
        std::vector<bool>* pReco_particle_ccinc_backtracked_goldenPion = nullptr;
        std::vector<bool>* pTrueCC1pi_recoParticle_isContained = nullptr;


        // tree->SetBranchAddress("reco_particle_ccinc_muonBDTScore", &pReco_particle_ccinc_muonBDTScore);
        // tree->SetBranchAddress("reco_particle_ccinc_protonBDTScore", &pReco_particle_ccinc_protonBDTScore);
        // tree->SetBranchAddress("reco_particle_ccinc_goldenPionBDTScore", &pReco_particle_ccinc_goldenPionBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_logBragg_pToMIP", &pReco_particle_ccinc_logBragg_pToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_logBragg_piToMIP", &pReco_particle_ccinc_logBragg_piToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_truncMeandEdx", &pReco_particle_ccinc_truncMeandEdx);
        tree->SetBranchAddress("reco_particle_ccinc_wiggliness", &pReco_particle_ccinc_wiggliness);
        tree->SetBranchAddress("reco_particle_ccinc_trackScore", &pReco_particle_ccinc_trackScore);
        tree->SetBranchAddress("reco_particle_ccinc_nDescendents", &pReco_particle_ccinc_nDescendents);
        tree->SetBranchAddress("reco_particle_ccinc_generation", &pReco_particle_ccinc_generation);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_pdg", &pReco_particle_ccinc_backtracked_pdg);
        // tree->SetBranchAddress("reco_particle_ccinc_backtracked_momentum", &pReco_particle_ccinc_backtracked_momentum);
        // tree->SetBranchAddress("reco_particle_ccinc_backtracked_cosTheta", &pReco_particle_ccinc_backtracked_cosTheta);
        // tree->SetBranchAddress("reco_particle_ccinc_backtracked_phi", &pReco_particle_ccinc_backtracked_phi);

        tree->SetBranchAddress("reco_particle_ccinc_backtracked_goldenPion", &pReco_particle_ccinc_backtracked_goldenPion);
        tree->SetBranchAddress("trueCC1pi_recoParticle_isContained_vector", &pTrueCC1pi_recoParticle_isContained);


        Float_t spline_weight, tuned_cv_weight;
        Bool_t isTrueCC1pi, isTrueGoldenPion;
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);
        tree->SetBranchAddress("true_cc1pi", &isTrueCC1pi);
        // tree->SetBranchAddress("true_golden_cc1pi", &isTrueGoldenPion);

        // Get the total number of entries
        std::cout<<"WARNING - Only using 1% of the entries!!!!!!!!!"<<std::endl;
        const auto nEntries = tree->GetEntries()/100;

        // Loop over the entries in the tree
        for (Long64_t i = 0; i < nEntries; i++)
        {
            tree->GetEntry(i);

            // Update the progress bar at every percent
            if (i % (nEntries / 100) == 0)
            {
                const auto progress = static_cast<int>(100.0 * i / nEntries);
                std::cout << "\r[" << std::string(progress, '|') << std::string(100 - progress, ' ') << "] " << progress << "%" << std::flush;
            }

            // Apply the condition
            if (tree->GetLeaf("passed_topologicalScoreCC")->GetValue())
            {
                if(!isTrueCC1pi) // Only true CC1pi events
                    continue;

                auto eventWeight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;
                eventWeight *= runWeight; // Apply the run weight

                if(pTrueCC1pi_recoParticle_isContained->size()!=0 && (pReco_particle_ccinc_generation->size() != pTrueCC1pi_recoParticle_isContained->size()))
                    throw std::runtime_error("pReco_particle_ccinc_generation" + std::to_string(pReco_particle_ccinc_generation->size()) + " and pTrueCC1pi_recoParticle_isContained" + std::to_string(pTrueCC1pi_recoParticle_isContained->size()) + " vectors are not the same size");

                // Loop over the values in the vector and fill the histogram
                for (Long64_t v = 0; v< pReco_particle_ccinc_generation->size(); v++)
                {
                    const auto generation = pReco_particle_ccinc_generation->at(v);
                    const auto isContained = pTrueCC1pi_recoParticle_isContained->at(v);
                    if(generation == 2 && isContained)
                    {
                        const auto pdg = std::abs(pReco_particle_ccinc_backtracked_pdg->at(v));
                        std::string particle;

                        switch (pdg)
                        {
                            case 13:
                                particle = "Mu";
                                break;
                            case 211:
                            {
                                const auto isGolden = pReco_particle_ccinc_backtracked_goldenPion->at(v);
                                particle = isGolden ? "Golden Pi" : "Pi";
                                break;
                            }
                            case 2212:
                                particle = "P";
                                break;
                            default:
                            {
                                if(pdg !=2147483647) std::cout << "DEBUG Default particle PDG: " << pdg << std::endl;
                                particle = "EXT";
                                break;
                            }
                        }

                        histograms1D[run]["logBragg_pToMIP"][particle]->Fill(pReco_particle_ccinc_logBragg_pToMIP->at(v), eventWeight);
                        histograms1D[run]["logBragg_piToMIP"][particle]->Fill(pReco_particle_ccinc_logBragg_piToMIP->at(v), eventWeight);
                        histograms1D[run]["truncMeandEdx"][particle]->Fill(pReco_particle_ccinc_truncMeandEdx->at(v), eventWeight);
                        histograms1D[run]["wiggliness"][particle]->Fill(pReco_particle_ccinc_wiggliness->at(v), eventWeight);
                        histograms1D[run]["trackScore"][particle]->Fill(pReco_particle_ccinc_trackScore->at(v), eventWeight);
                        histograms1D[run]["nDescendents"][particle]->Fill(pReco_particle_ccinc_nDescendents->at(v), eventWeight);

                        histograms1D["0"]["logBragg_pToMIP"][particle]->Fill(pReco_particle_ccinc_logBragg_pToMIP->at(v), eventWeight);
                        histograms1D["0"]["logBragg_piToMIP"][particle]->Fill(pReco_particle_ccinc_logBragg_piToMIP->at(v), eventWeight);
                        histograms1D["0"]["truncMeandEdx"][particle]->Fill(pReco_particle_ccinc_truncMeandEdx->at(v), eventWeight);
                        histograms1D["0"]["wiggliness"][particle]->Fill(pReco_particle_ccinc_wiggliness->at(v), eventWeight);
                        histograms1D["0"]["trackScore"][particle]->Fill(pReco_particle_ccinc_trackScore->at(v), eventWeight);
                        histograms1D["0"]["nDescendents"][particle]->Fill(pReco_particle_ccinc_nDescendents->at(v), eventWeight);
                    }
                } // End of loop over vector entries (aka particles)
            } // End of if signal
        } // End of loop over tree entries
    } // End of loop over files

    // Create a particle to color map
    std::map<std::string, int> particleColors = {
        {"P", kOrange+1},
        {"Mu", kBlue},
        {"Pi", kMagenta},
        {"Golden Pi", kGreen},
        {"EXT", kBlack},
    };

    // Map with max y values for each plot; map: run, variable
    std::map<std::string, std::map<std::string, double>> maxYValues;
    
    // Area normalise all plots
    for (const auto& run : runs)
    {
        for (const auto& variable : variables)
        {
            for (const auto& particle : particles)
            {
                auto h = histograms1D[run][variable][particle].get();
                const auto integral = h->Integral();
                h->Scale(1.0/integral);
                const auto maxY = histograms1D[run][variable][particle]->GetMaximum();
                if (maxY > maxYValues[run][variable])
                    maxYValues[run][variable] = maxY;
            }
        }
    }

    // Loop over the runs and variables in histograms1D and for each draw all particle types on one canvas
    for (const auto& run : runs)
    {
        for (const auto& variable : variables)
        {
            // Create canvas here and then plot the particle histograms on it
            auto firstHist = true;
            const auto outName = "BDTStudy1DParticleInputVariables_" + run + "_" + variable;
            auto c = new TCanvas(outName.c_str(), "", 800, 600);
            auto legend = variable == "trackScore" ? new TLegend(0.3, 0.7, 0.5, 0.9) : new TLegend(0.7, 0.7, 0.9, 0.9);
            const auto [nBinsVar, xMin, xMax, xAxisLog, yAxisLog] = binningInfo.at(variable);
            for (const auto& particle : particles)
            {
                auto h = histograms1D[run][variable][particle].get();
                h->SetStats(0); // Remove stats box
                h->SetLineColor(particleColors[particle]);
                h->SetLineWidth(2);
                std::string style = "HIST";
    

                if(!firstHist)
                {
                    style += " SAME";
                }
                else
                {
                    firstHist = false;
                    // Set the title
                    const std::string variableName = variable == "logBragg_pToMIP" ? "log(R_{p}/MIP)" : (variable == "logBragg_piToMIP" ? "log(R_{#pi}/MIP)" : (variable == "truncMeandEdx" ? "Truncated Mean dE/dx" : (variable == "wiggliness" ? "Wiggliness" : (variable == "trackScore" ? "Track Score" : (variable == "nDescendents" ? "# Descendents" : "Unknown Variable")))));
                    const std::string runName = run == "0" ? "All Runs" : ("Run " + run);
                    std::string title = variableName + " for " + runName;
                    h->SetTitle(("True CC1#pi Events Passing the #nu_#mu Preselection for " + runName).c_str());
                    // Set the y axis label
                    std::string yAxisLabel = "# Events (area normalised)";
                    h->SetYTitle(yAxisLabel.c_str());
                    // Set the x axis label
                    std::string xAxisLabel = variableName;
                    h->SetXTitle(xAxisLabel.c_str());
                    // Set y axis range
                    h->GetYaxis()->SetRangeUser(0, 1.2*maxYValues[run][variable]);
                }

                if(yAxisLog) h->SetMinimum(0.0001);
                h->Draw(style.c_str());

                // Make a copy of the histogram and set the fill color and style
                auto hCopy = (TH1D*)h->Clone();
                const auto transparentColor = TColor::GetColorTransparent(particleColors[particle], 0.3); // 30% transparency
                hCopy->SetFillColor(transparentColor); // Set the fill color
                std::string errorStyle = "E2 SAME";
                if(yAxisLog) hCopy->SetMinimum(0.0001);
                hCopy->Draw(errorStyle.c_str());

                // Add the entry to the legend
                const std::string particleName = particle == "Golden Pi" ? "Golden Pion" : (particle == "Pi" ? "Other Pion" : (particle == "Mu" ? "Muon" :(particle == "P" ? "Proton" : "EXT")));
                legend->AddEntry(h, particleName.c_str(), "l");
            }

            legend->Draw();
            c->SetLogx(xAxisLog);
            c->SetLogy(yAxisLog);
            c->SaveAs(("plots/" + outName + ".pdf").c_str());
            c->SaveAs(("plots/" + outName + ".png").c_str());
            c->SaveAs(("plots/" + outName + ".C").c_str());
        }
    }
}