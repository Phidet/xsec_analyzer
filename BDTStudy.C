#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <map>
#include <string>
#include <memory>

void BDTStudy() 
{

    const std::vector versions = {"detVar_DecPelee2", // Full run 1-5 trained BDTs
                                    "detVar_DecPelee3", // Run 1 only trained BDTs
                                    "detVar_DecPelee3_withDetVar", // Run 1 only trained BDTs with additional detector variations (SCE & Recomb2)
                                    "detVar_DecPelee4", // Run 1 only trained BDTs
                                    "detVar_DecPelee4_withDetVar"}; // Run 1 only trained BDTs with additional detector variations (Recomb2 with weights scaled to overlay)
    const auto version = versions.at(3);

    // List of files
    std::vector<std::tuple<std::string, std::string>> files = {
        std::make_tuple("detVar_SCE", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "/detvar_run3_peleeTuple_uboone_v08_00_00_73_SCE_ubcc1pi.root"),
        std::make_tuple("detVar_CVextra", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "/detvar_run3_peleeTuple_uboone_v08_00_00_73_CVextra_ubcc1pi.root"),
        std::make_tuple("detVar_Recomb2", std::string("/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/") + version + "/detvar_run3_peleeTuple_uboone_v08_00_00_73_Recomb2_ubcc1pi.root"),
    };

    // const std::vector<float> overlayPOT {1.28278e+21, 1.01592e+21, 1.31355e+21, 2.48772e+20, 7.64399e+20, 4.64842e+20, 8.66958e+20, 9.70631e+20}; // run 1, 2, 3, 4a, 4b, 4c, 4d, 5
    // const std::vector<float> dataBNBTor860WCut {1.669e+20, 2.616e+20, 2.562e+20, 3.607e+19, 1.39e+20, 8.586e+19, 4.919e+19, 1.23e+20};
    // std::vector<float> ratioPOT;
    // for (size_t i = 0; i < overlayPOT.size(); i++)
    // {
    //     ratioPOT.push_back(dataBNBTor860WCut[i] / overlayPOT[i]);
    // }
    // const std::vector<float> weights = {ratioPOT.at(0), ratioPOT.at(1), ratioPOT.at(2)};
    
    // Map for the max y axis values to use for the 1D histograms
    // map: bdt
    std::map<std::string, double> yAxisRange = {{"muonBDTScore", 0.0}, {"protonBDTScore", 0.0}, {"goldenPionBDTScore", 0.0}};

    // Loop over the files
    for (const auto& [name, filePath] : files) {
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


        std::vector<std::string> bdts = {"muonBDTScore", "protonBDTScore", "goldenPionBDTScore"};
        std::vector<std::string> variables = {"logBragg_pToMIP", "logBragg_piToMIP", "truncMeandEdx", "wiggliness", "trackScore", "nDescendents"};
        std::vector<std::string> particles = {"P", "Pi", "Mu"};

        std::map<std::string, std::tuple<int, double, double>> binningInfo = {
            {"muonBDTScore", std::make_tuple(20, -1, 1)},
            {"protonBDTScore", std::make_tuple(20, -1, 1)},
            {"goldenPionBDTScore", std::make_tuple(20, -1, 1)},
            {"logBragg_pToMIP", std::make_tuple(20, -10, 10)},
            {"logBragg_piToMIP", std::make_tuple(20, -10, 10)},
            {"truncMeandEdx", std::make_tuple(20, 0, 7.0)},
            {"wiggliness", std::make_tuple(20, 0, 0.25)},
            {"trackScore", std::make_tuple(20, 0, 1.0)},
            {"nDescendents", std::make_tuple(5, 0, 4)}
        };

        // map: bdt, variable, particle 
        std::map<std::string, std::map<std::string, std::map<std::string, std::unique_ptr<TH2D>>>> histograms2D;
        // map: bdt, particle
        std::map<std::string, std::map<std::string, std::unique_ptr<TH1D>>> histograms1D;
        for (const auto& bdt : bdts) {
            auto [nBinsBDT, xMinBDT, xMaxBDT] = binningInfo[bdt];
            for (const auto& particle : particles) {
                for (const auto& variable : variables) {
                    std::string histName = name + "_" + bdt + "_" + variable + "_" + particle + "_" + version;
                    std::string title = name + " reco Particles Backtracked to True " + particle + " in Signal Events; " + bdt + " BDT Score; ";
                    auto [nBinsVar, xMinVar, xMaxVar] = binningInfo[variable];
                    histograms2D[bdt][variable][particle] = std::make_unique<TH2D>(histName.c_str(), title.c_str(), nBinsBDT, xMinBDT, xMaxBDT, nBinsVar, xMinVar, xMaxVar);
                }
                histograms1D[bdt][particle] = std::make_unique<TH1D>((name + "_" + bdt + "_" + particle).c_str(), (name + " reco Particles Backtracked to True " + particle + " in Signal Events; " + bdt + " BDT Score").c_str(), nBinsBDT, xMinBDT, xMaxBDT);
            }
        }

        // Get the total number of entries
        const auto nEntries = tree->GetEntries();
        // std::cout<< "DEBUG FillParticleHistogram Point 3" << std::endl;
        // Loop over the entries in the tree
        for (Long64_t i = 0; i < nEntries; i++)
        {
            tree->GetEntry(i);
            // std::cout<< "DEBUG FillParticleHistogram Point 3.1" << std::endl;

            // Update the progress bar at every percent
            if (i % (nEntries / 100) == 0)
            {
                const auto progress = static_cast<int>(100.0 * i / nEntries);
                std::cout << "\r[" << std::string(progress, '|') << std::string(100 - progress, ' ') << "] " << progress << "%" << std::flush;
            }
            // std::cout<< "DEBUG FillParticleHistogram Point 3.2" << std::endl;

            // std::cout<< "DEBUG FillParticleHistogram Point 4" << std::endl;
            // std::cout << "DEBUG FillParticleHistogram Point pMuonBDTScore->size(): " << pMuonBDTScore->size() << std::endl;
            // Apply the condition
            if (tree->GetLeaf("cc1pi_signal")->GetValue() && !tree->GetLeaf("isTrainingEvent")->GetValue())
            {
                // std::cout<< "DEBUG FillParticleHistogram Point pReco_particle_ccinc_protonBDTScore->size(): " << pReco_particle_ccinc_protonBDTScore->size() << std::endl;
                // Loop over the values in the vector and fill the histogram
                for (Long64_t v = 0; v< pReco_particle_ccinc_protonBDTScore->size(); v++)
                {
                    const auto generation = pReco_particle_ccinc_generation->at(v);
                    if(generation == 2)
                    {
                        const auto pdg = std::abs(pReco_particle_ccinc_backtracked_pdg->at(v));
                        std::string particle = "";
                        switch (pdg)
                        {
                            case 2212:
                                particle = "P";
                                break;
                            case 13:
                                particle = "Mu";
                                break;
                            case 211:
                                particle = "Pi";
                                break;
                        }
                        
                        for (const auto& bdt : bdts) {
                            float bdtScore;
                            if (bdt == "muonBDTScore") {
                                bdtScore = pReco_particle_ccinc_muonBDTScore->at(v);
                            }
                            else if (bdt == "protonBDTScore") {
                                bdtScore = pReco_particle_ccinc_protonBDTScore->at(v);
                            }
                            else if (bdt == "goldenPionBDTScore") {
                                bdtScore = pReco_particle_ccinc_goldenPionBDTScore->at(v);
                            }
                            else throw std::runtime_error("Unknown BDT");

                            if(particle != "") // For the 2D plots require a each particle to be identified
                            {      
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

                                    histograms2D[bdt][variable][particle]->Fill(bdtScore, variableValue);
                                } // End of loop over variables
                            } // End of if particle != ""
                            // 1D BDT Histograms
                            if(particle != "") histograms1D[bdt][particle]->Fill(bdtScore);
                        } // End of loop over BDTs
                    }
                } // End of loop over vector entries (aka particles)
            } // End of if signal
            // std::cout<< "DEBUG FillParticleHistogram Point 5" << std::endl;
        } // End of loop over tree entries

        // Draw the histograms2D and save them
        gStyle->SetOptStat(0);
        for (const auto& bdt : bdts) {
            for (const auto& variable : variables) {
                for (const auto& particle : particles) {
                    auto& histogram = histograms2D[bdt][variable][particle];
                    TCanvas canvas;
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

                    canvas.SaveAs(("plots/BDTStudy_" + name + "_" + bdt + "_" + variable + "_" + particle + "_" + version + ".pdf").c_str());
                    canvas.SaveAs(("plots/BDTStudy_" + name + "_" + bdt + "_" + variable + "_" + particle + "_" + version + ".png").c_str());
                    canvas.SaveAs(("plots/BDTStudy_" + name + "_" + bdt + "_" + variable + "_" + particle + "_" + version + ".C").c_str());
                }
            }
        }

        // Draw the 1D histograms and save them
        gStyle->SetOptStat(0);
        std::map<std::string, int> particleColors = {{"P", kRed}, {"Pi", kBlue}, {"Mu", kGreen}};
        for (const auto& bdt : bdts) {
            TCanvas canvas1D;
            bool isFirstHistogram = true;

            // Set x, y and title labels
            histograms1D[bdt]["P"]->GetXaxis()->SetTitle((bdt + " BDT Score").c_str());
            histograms1D[bdt]["P"]->GetYaxis()->SetTitle("Fraction of Particles");
            histograms1D[bdt]["P"]->SetTitle((name + " reco Particles Backtracked to True Particle in Signal Events; " + bdt + " BDT Score").c_str());

            for (const auto& particle : particles) {
                auto& histogram = histograms1D[bdt][particle];
                histogram->Scale(1.0 / histogram->Integral());
                histogram->SetLineColor(particleColors[particle]);
                histogram->SetMarkerColor(particleColors[particle]);
            }

            if(yAxisRange.at(bdt) == 0.0) // Only set for the first file. This assumes that all the files have roughly the same max value
            {
                double maxVal = 0;
                for (const auto& particle : particles) {
                    maxVal = std::max(maxVal, histograms1D[bdt][particle]->GetMaximum());
                }
                yAxisRange.at(bdt) = 1.1*maxVal;
            }

            for (const auto& particle : particles) {
                auto& histogram = histograms1D[bdt][particle];
                if (isFirstHistogram) {
                    histogram->GetYaxis()->SetRangeUser(0, yAxisRange.at(bdt));
                    histogram->Draw("Hist E");
                    isFirstHistogram = false;
                } else {
                    histogram->Draw("Hist E SAME");
                }
            }
            TLegend legend(0.7, 0.7, 0.9, 0.9);
            legend.SetHeader("Particle");
            legend.AddEntry(histograms1D[bdt]["P"].get(), "Proton", "l");
            legend.AddEntry(histograms1D[bdt]["Pi"].get(), "Pion", "l");
            legend.AddEntry(histograms1D[bdt]["Mu"].get(), "Muon", "l");
            legend.Draw();

            canvas1D.SaveAs(("plots/BDTStudy_" + name + "_" + bdt + "_" + version + ".pdf").c_str());
            canvas1D.SaveAs(("plots/BDTStudy_" + name + "_" + bdt + "_" + version + ".png").c_str());
            canvas1D.SaveAs(("plots/BDTStudy_" + name + "_" + bdt + "_" + version + ".C").c_str());
        }

        // // Reset all the histograms
        // for (const auto& bdt : bdts) {
        //     for (const auto& particle : particles) {
        //         for (const auto& variable : variables) {
        //             histograms2D[bdt][variable][particle].reset();
        //         }
        //         histograms1D[bdt][particle].reset();
        //     }
        // }

        
    } // End of loop over files
}