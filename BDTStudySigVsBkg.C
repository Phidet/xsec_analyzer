#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <map>
#include <string>
#include <memory>

void BDTStudySigVsBkg() 
{

    const std::vector versions = {""};
    
    const auto version = versions.at(0);


    // List of files; map: run, path
    std::vector<std::tuple<std::string, std::string>> files = {
        std::make_tuple("1", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root"),
        std::make_tuple("2", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root"),
        std::make_tuple("3", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root"),
        std::make_tuple("4a", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "overlay_peleeTuple_uboone_v08_00_00_73_run4a_nu_ubcc1pi.root"),
        std::make_tuple("4b", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "overlay_peleeTuple_uboone_v08_00_00_73_run4b_nu_ubcc1pi.root"),
        std::make_tuple("4c", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "overlay_peleeTuple_uboone_v08_00_00_70_run4c_nu_ubcc1pi.root"),
        std::make_tuple("4d", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "overlay_peleeTuple_uboone_v08_00_00_70_run4d_nu_ubcc1pi.root"),
        std::make_tuple("5", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root"),
    };

    const std::vector runs {"1", "2", "3", "4a", "4b", "4c", "4d", "5"};

    const std::vector<std::string> bdts = {"muonBDTScore", "protonBDTScore", "goldenPionBDTScore"};
    const std::vector<std::string> variables = {"logBragg_pToMIP", "logBragg_piToMIP", "truncMeandEdx", "wiggliness", "trackScore", "nDescendents"};
    const std::vector<std::string> particles = {"Signal", "Background"};

    // map: variable, tuple(nBins, xMin, xMax, axisLog)
    const std::map<std::string, std::tuple<int, double, double, bool>> binningInfo = {
        {"muonBDTScore", std::make_tuple(20, -1, 1, false)},
        {"protonBDTScore", std::make_tuple(20, -1, 1, false)},
        {"goldenPionBDTScore", std::make_tuple(20, -1, 1, false)},
        {"logBragg_pToMIP", std::make_tuple(20, -10, 10, false)},
        {"logBragg_piToMIP", std::make_tuple(20, -10, 10, false)},
        {"truncMeandEdx", std::make_tuple(20, 0, 20.0, false)},
        {"wiggliness", std::make_tuple(20, 0.00001, 0.25, true)},
        {"trackScore", std::make_tuple(20, 0, 1.0, false)},
        {"nDescendents", std::make_tuple(4, 0, 4, false)}
    };

    // Loop over the files
    for (const auto& [run, filePath] : files) {
        // const auto weight = weights.at(f);
        // std::cout << "Processing file " << filePath << " with weight " << weight << std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(filePath.c_str());
        if (!tFile || tFile->IsZombie()) {
            std::cerr << "Error: Unable to open file: " << filePath << std::endl;
            return;
        }

        // std::cout<< "DEBUG FillParticleHistogram Point 2.1" << std::endl;

        TTree* tree = (TTree*)tFile->Get("stv_tree");
        if (!tree) {
            std::cerr << "Error: Unable to get tree from file: " << filePath << std::endl;
            return;
        }

        // std::cout<< "DEBUG FillParticleHistogram Point 2.2" << std::endl;

        std::vector<float>* pReco_particle_ccinc_muonBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_protonBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_goldenPionBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_logBragg_pToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_logBragg_piToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_truncMeandEdx = nullptr;
        std::vector<float>* pReco_particle_ccinc_wiggliness = nullptr;
        std::vector<float>* pReco_particle_ccinc_trackScore = nullptr;
        std::vector<int>* pReco_particle_ccinc_nDescendents = nullptr;
        std::vector<int>* pReco_particle_ccinc_generation = nullptr;
        std::vector<int>* pReco_particle_ccinc_backtracked_pdg = nullptr;
        std::vector<float>* pReco_particle_ccinc_backtracked_momentum = nullptr;
        std::vector<float>* pReco_particle_ccinc_backtracked_cosTheta = nullptr;
        std::vector<float>* pReco_particle_ccinc_backtracked_phi = nullptr;
        std::vector<bool>* pReco_particle_ccinc_backtracked_goldenPion = nullptr;

        // std::cout<< "DEBUG FillParticleHistogram Point 2.3" << std::endl;

        // std::cout<< "DEBUG FillParticleHistogram Point 2.4" << std::endl;

        tree->SetBranchAddress("reco_particle_ccinc_muonBDTScore", &pReco_particle_ccinc_muonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_protonBDTScore", &pReco_particle_ccinc_protonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_goldenPionBDTScore", &pReco_particle_ccinc_goldenPionBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_logBragg_pToMIP", &pReco_particle_ccinc_logBragg_pToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_logBragg_piToMIP", &pReco_particle_ccinc_logBragg_piToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_truncMeandEdx", &pReco_particle_ccinc_truncMeandEdx);
        tree->SetBranchAddress("reco_particle_ccinc_wiggliness", &pReco_particle_ccinc_wiggliness);
        tree->SetBranchAddress("reco_particle_ccinc_trackScore", &pReco_particle_ccinc_trackScore);
        tree->SetBranchAddress("reco_particle_ccinc_nDescendents", &pReco_particle_ccinc_nDescendents);
        tree->SetBranchAddress("reco_particle_ccinc_generation", &pReco_particle_ccinc_generation);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_pdg", &pReco_particle_ccinc_backtracked_pdg);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_momentum", &pReco_particle_ccinc_backtracked_momentum);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_cosTheta", &pReco_particle_ccinc_backtracked_cosTheta);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_phi", &pReco_particle_ccinc_backtracked_phi);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_goldenPion", &pReco_particle_ccinc_backtracked_goldenPion);

        Float_t spline_weight, tuned_cv_weight;
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);

        // map: bdt, variable, particle 
        std::map<std::string, std::map<std::string, std::map<std::string, std::unique_ptr<TH2D>>>> histograms2D;
        // map: bdt, particle
        std::map<std::string, std::map<std::string, std::unique_ptr<TH1D>>> histograms1D;
        for (const auto& bdt : bdts) {
            const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
            for (const auto& particle : particles) {
                for (const auto& variable : variables) {
                    std::string histName = run + "_" + bdt + "_" + variable + "_" + particle + "_" + version;
                    std::string title = run + " reco Particles Backtracked to True " + particle + " in Signal Events; " + bdt + " BDT Score; ";
                    const auto [nBinsVar, yMin, yMax, yAxisLog] = binningInfo.at(variable);
                    if(yAxisLog)
                    {
                        // Create a vector to hold the bin edges
                        std::vector<double> binEdges(nBinsVar + 1);
                        double logMin = TMath::Log10(yMin);
                        double logMax = TMath::Log10(yMax);
                        double binWidth = (logMax - logMin) / nBinsVar;
                        for (int i = 0; i <= nBinsVar; i++) {
                            binEdges[i] = TMath::Power(10, logMin + i * binWidth);
                        }

                        // Create the histogram with logarithmic binning
                        histograms2D[bdt][variable][particle] = std::make_unique<TH2D>(histName.c_str(), title.c_str(), nBinsBDT, xMin, xMax, nBinsVar, binEdges.data());
                    }
                    else
                    {
                        histograms2D[bdt][variable][particle] = std::make_unique<TH2D>(histName.c_str(), title.c_str(), nBinsBDT, xMin, xMax, nBinsVar, yMin, yMax);
                    }
                }
            }
        }

        // Get the total number of entries
        // std::cout<<"WARNING - Only using 5% of the entries!!!!!!!!!"<<std::endl;
        const auto nEntries = tree->GetEntries();
        // std::cout<< "DEBUG FillParticleHistogram Point 3" << std::endl;
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
            if (tree->GetLeaf("passed_topologicalScoreCC")->GetValue() && !tree->GetLeaf("isTrainingEvent")->GetValue())
            {
                const auto eventWeight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;
                // Loop over the values in the vector and fill the histogram
                for (Long64_t v = 0; v< pReco_particle_ccinc_protonBDTScore->size(); v++)
                {
                    const auto generation = pReco_particle_ccinc_generation->at(v);
                    if(generation == 2)
                    {
                        const auto pdg = std::abs(pReco_particle_ccinc_backtracked_pdg->at(v));
                        std::string particle;
                        
                        for (const auto& bdt : bdts) {
                            float bdtScore;
                            if (bdt == "muonBDTScore") {
                                bdtScore = pReco_particle_ccinc_muonBDTScore->at(v);
                                particle = std::abs(pdg) == 13 ? "Signal" :"Background";
                            }
                            else if (bdt == "protonBDTScore") {
                                bdtScore = pReco_particle_ccinc_protonBDTScore->at(v);
                                particle = std::abs(pdg) == 2212 ? "Signal" :"Background";
                            }
                            else if (bdt == "goldenPionBDTScore") {
                                const auto isGoldenPion = pReco_particle_ccinc_backtracked_goldenPion->at(v);
                                const auto isPion = (std::abs(pdg) == 211);
                                if(isGoldenPion && (!isPion))
                                    throw std::runtime_error("Golden pion is not a pion");
                                bdtScore = pReco_particle_ccinc_goldenPionBDTScore->at(v);
                                particle = isGoldenPion ? "Signal" :"Background";
                            }
                            else throw std::runtime_error("Unknown BDT");
    
                            for (const auto& variable : variables) {
                                float variableValue;
                                if(variable == "logBragg_pToMIP") {
                                    variableValue = pReco_particle_ccinc_logBragg_pToMIP->at(v);
                                }
                                else if(variable == "logBragg_piToMIP") {
                                    variableValue = pReco_particle_ccinc_logBragg_piToMIP->at(v);
                                }
                                else if(variable == "truncMeandEdx") {
                                    variableValue = pReco_particle_ccinc_truncMeandEdx->at(v);
                                }
                                else if(variable == "wiggliness") {
                                    variableValue = pReco_particle_ccinc_wiggliness->at(v);
                                }
                                else if(variable == "trackScore") {
                                    variableValue = pReco_particle_ccinc_trackScore->at(v);
                                }
                                else if(variable == "nDescendents") {
                                    variableValue = pReco_particle_ccinc_nDescendents->at(v);
                                }
                                else throw std::runtime_error("Unknown variable");

                                histograms2D[bdt][variable][particle]->Fill(bdtScore, variableValue, eventWeight);
                            } // End of loop over variables#

                        } // End of loop over BDTs
                    }
                } // End of loop over vector entries (aka particles)
            } // End of if signal
        } // End of loop over tree entries

        // Draw the histograms2D and save them
        gStyle->SetOptStat(0);
        for (const auto& bdt : bdts) {
            for (const auto& variable : variables) {
                for (const auto& particle : particles) {
                    auto& histogram = histograms2D[bdt][variable][particle];
                    TCanvas canvas;
                    const auto [nBinsVar, yMin, yMax, yAxisLog] = binningInfo.at(variable);
                    if(yAxisLog)
                        canvas.SetLogy();
                    // Disable stats box
                    histogram->Draw("colz");
                    // Set x, y and title labels
                    if(bdt == "muonBDTScore") histogram->GetXaxis()->SetTitle("Muon BDT Score");
                    else if(bdt == "protonBDTScore") histogram->GetXaxis()->SetTitle("Proton BDT Score");
                    else if(bdt == "goldenPionBDTScore") histogram->GetXaxis()->SetTitle("Golden Pion BDT Score");

                    if (variable == "logBragg_pToMIP") {
                        histogram->GetYaxis()->SetTitle("logBragg_pToMIP");
                    } else if (variable == "logBragg_piToMIP") {
                        histogram->GetYaxis()->SetTitle("logBragg_piToMIP");
                    } else if (variable == "truncMeandEdx") {
                        histogram->GetYaxis()->SetTitle("Truncated Mean dE/dx");
                    } else if (variable == "wiggliness") {
                        histogram->GetYaxis()->SetTitle("Wiggliness");
                    } else if (variable == "trackScore") {
                        histogram->GetYaxis()->SetTitle("Track Score");
                    } else if (variable == "nDescendents") {
                        histogram->GetYaxis()->SetTitle("Number of Descendents");
                    }

                    histogram->SetTitle("");

                    canvas.SaveAs(("plots/BDTStudy_SigVsBkg_nu_overlay_run" + run + "_" + bdt + "_" + variable + "_" + particle + ".pdf").c_str());
                    canvas.SaveAs(("plots/BDTStudy_SigVsBkg_nu_overlay_run" + run + "_" + bdt + "_" + variable + "_" + particle + ".png").c_str());
                    canvas.SaveAs(("plots/BDTStudy_SigVsBkg_nu_overlay_run" + run + "_" + bdt + "_" + variable + "_" + particle + ".C").c_str());
                } // End of loop over particle types
            } // End of loop over variables
        } // End of loop over BDTs
    } // End of loop over files
}