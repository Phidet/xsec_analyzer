#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <map>
#include <string>
#include <memory>

#include "EventCategory.hh" 

void WigglinessStudy1D() 
{
    const EventCategoryInterpreter& eventCategoryInterpreter = EventCategoryInterpreter::Instance();

    const std::string rootPath = "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/22Feb24/";
    // tuple: type, run, file path, run weight
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        // std::make_tuple("beam on", "1",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "2",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "3",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "4bcd", rootPath + "bnb_beam_on_peleeTuple_uboone_run4bcd_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "5",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root", 1.f),

        // std::make_tuple("beam off", "1",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi.root", 0.56421),
        // std::make_tuple("beam off", "2",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 0.40202),
        // std::make_tuple("beam off", "3",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 0.30657),
        // std::make_tuple("beam off", "4bcd", rootPath + "bnb_beam_off_peleeTuple_uboone_run4bcd_ubcc1pi.root", 0.30175),
        // std::make_tuple("beam off", "5",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root", 0.32807),

        std::make_tuple("nu mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root", 0.13011),
        std::make_tuple("nu mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root", 0.25750),
        std::make_tuple("nu mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root", 0.20113),
        std::make_tuple("nu mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi.root", 0.13074),
        std::make_tuple("nu mc", "5",  rootPath + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root", 0.15196),

        // std::make_tuple("dirt mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root", 0.52806),
        // std::make_tuple("dirt mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root", 0.27521),
        // std::make_tuple("dirt mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root", 0.80892),
        // std::make_tuple("dirt mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi.root", 0.39701),
        // std::make_tuple("dirt mc", "5",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_ubcc1pi.root", 0.41280),
    };

    const std::vector<int> nDescendentOptions {0, 1};
    const std::vector<std::string> bdts = {"wiggliness"};
    // const std::vector<std::string> types = {"EXT", "Dirt", "E", "Photon", "K", "P", "Mu", "Pi", "Golden Pi", "Beam On"};
    const std::vector<std::string> types = {"Pi", "Golden Pi"};
    const std::vector<bool> dataOptions = {true, false};

    // map: variable, tuple(nBins, xMin, xMax, axisLog)
    const std::map<std::string, std::tuple<int, double, double, bool>> binningInfo = {
        // {"muonBDTScore", std::make_tuple(60, -1, 1, false)},
        // {"protonBDTScore", std::make_tuple(60, -1, 1, false)},
        // {"goldenPionBDTScore", std::make_tuple(60, -1, 1, false)},
        // {"logBragg_pToMIP", std::make_tuple(60, -8, 8, false)},
        // {"logBragg_piToMIP", std::make_tuple(60, -4, 7, false)},
        // {"truncMeandEdx", std::make_tuple(60, 0, 10.0, false)},
        {"wiggliness", std::make_tuple(60, 0.00005, 0.25, true)},
        // {"trackScore", std::make_tuple(60, 0, 1.0, false)},
        // {"nDescendents", std::make_tuple(4, 0, 4, false)}
    };

    // map: nDescendents, bdt, particle
    std::map<int, std::map<std::string, std::map<std::string, std::unique_ptr<TH1D>>>> histograms1D;
    for (const auto& nDescendents : nDescendentOptions) {
        for (const auto& bdt : bdts) {
            const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
            for (const auto& type : types) {
                    const auto name = "h_" + std::to_string(nDescendents) + "_" + bdt + "_" + type;
                    if(xAxisLog)
                    {
                        std::cout << "DEBUG: Creating histogram with logarithmic binning for " << name << std::endl;
                        // Create a vector to hold the bin edges
                        std::vector<double> binEdges(nBinsBDT + 1);
                        double logMin = TMath::Log10(xMin);
                        double logMax = TMath::Log10(xMax);
                        double binWidth = (logMax - logMin) / nBinsBDT;
                        for (int i = 0; i <= nBinsBDT; i++) {
                            binEdges[i] = TMath::Power(10, logMin + i * binWidth);
                        }

                        // Create the histogram with logarithmic binning
                        histograms1D[nDescendents][bdt][type] = std::make_unique<TH1D>((name + "_" + bdt + "_" + type).c_str(), (name + " reco Particles Backtracked to True " + type + " in Signal Events; " + bdt + " BDT Score").c_str(), nBinsBDT, binEdges.data());
                    }
                    else
                    {
                        histograms1D[nDescendents][bdt][type] = std::make_unique<TH1D>((name + "_" + bdt + "_" + type).c_str(), (name + " reco Particles Backtracked to True " + type + " in Signal Events; " + bdt + " BDT Score").c_str(), nBinsBDT, xMin, xMax);
                    }
                    histograms1D[nDescendents][bdt][type]->Sumw2();
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

        std::vector<int>* pReco_particle_ccinc_generation = nullptr;
        std::vector<int>* pReco_particle_ccinc_backtracked_pdg = nullptr;
        std::vector<bool>* pReco_particle_ccinc_backtracked_goldenPion = nullptr;
        std::vector<bool>* pReco_particle_ccinc_isContained = nullptr;

        std::vector<int>* pReco_particle_ccinc_nDescendents = nullptr;
        std::vector<float>* pReco_particle_ccinc_wiggliness = nullptr;

        tree->SetBranchAddress("reco_particle_ccinc_generation", &pReco_particle_ccinc_generation);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_pdg", &pReco_particle_ccinc_backtracked_pdg);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_goldenPion", &pReco_particle_ccinc_backtracked_goldenPion);
        tree->SetBranchAddress("reco_particle_ccinc_contained", &pReco_particle_ccinc_isContained);

        tree->SetBranchAddress("reco_particle_ccinc_nDescendents", &pReco_particle_ccinc_nDescendents);
        tree->SetBranchAddress("reco_particle_ccinc_wiggliness", &pReco_particle_ccinc_wiggliness);

        Float_t spline_weight, tuned_cv_weight;
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);


        // Get the total number of entries
        std::cout<<"WARNING - Only using 10\% of the entries!!!!!!!!!"<<std::endl;
        const auto nEntries = tree->GetEntries()/10;

        // Loop over the entries in the tree
        for (Long64_t i = 0; i < nEntries; i++)
        {
            tree->GetEntry(i);

            // Update the progress bar at every percent
            if (i % (nEntries / 100) == 0)
            {
                const auto progress = static_cast<int>(100.0 * i / nEntries);
                std::cout << "\r"<< " " << sampleType << " run " << run <<"[" << std::string(progress, '|') << std::string(100 - progress, ' ') << "] " << progress << "%" << std::flush;
            }

            // Apply the condition
            if (tree->GetLeaf("passed_topologicalScoreCC")->GetValue())
            {
                // const auto isTrainingEvent = tree->GetLeaf("isTrainingEvent")->GetValue();
                // if(isTrainingEvent)
                //     continue;

                auto eventWeight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;
                eventWeight *= runWeight; // Apply the run weight

                // Loop over the values in the vector and fill the histogram
                for (Long64_t v = 0; v< pReco_particle_ccinc_generation->size(); v++)
                {
                    // Set the type to the type determined by the file first
                    std::string type;
                    

                    const auto generation = pReco_particle_ccinc_generation->at(v);
                    const auto isContained = pReco_particle_ccinc_isContained->at(v);
                    if(pReco_particle_ccinc_generation->size() != pReco_particle_ccinc_isContained->size())
                        throw std::runtime_error("ERROR: pReco_particle_ccinc_generation and pReco_particle_ccinc_isContained have different sizes");
                    if(generation == 2 && isContained)
                    {
                        if(type == "")
                        {
                            const auto pdg = std::abs(pReco_particle_ccinc_backtracked_pdg->at(v));

                            switch (pdg)
                            {
                                // case 13:
                                //     type = "Mu";
                                //     break;
                                case 211:
                                {
                                    const auto isGolden = pReco_particle_ccinc_backtracked_goldenPion->at(v);
                                    type = isGolden ? "Golden Pi" : "Pi";
                                    break;
                                }
                                // case 2212:
                                //     type = "P";
                                //     break;
                                // default:
                                //     type = "Other";
                            }
                        }
                        
                        if (type != "")
                        {
                            if(pReco_particle_ccinc_nDescendents->at(v) >= 0)
                            {
                                const auto nDescendents = pReco_particle_ccinc_nDescendents->at(v) == 0 ? 0 : 1;
                                histograms1D.at(nDescendents).at("wiggliness").at(type)->Fill(pReco_particle_ccinc_wiggliness->at(v), eventWeight);
                            }
                        }
                    } // End of if contained and generation 2
                } // End of loop over vector entries (aka types)
            } // End of if signal
        } // End of loop over tree entries
    } // End of loop over files


    std::cout << "DEBUG: Finished filling histograms" << std::endl;

    // Create a particle to color map
    std::map<std::string, int> particleColors = {
        // {"P", kOrange+1},
        // {"Mu", kBlue},
        {"Pi", kMagenta},
        // {"E", kCyan},
        // {"Photon", kYellow},
        // {"K", kMagenta},
        {"Golden Pi", kGreen},
        // {"EXT", kBlack},
        // {"Dirt", kOrange-6},
        // {"Beam On", kBlack},
        // {"Other", kGray}
    };


    std::cout << "DEBUG: Finished setting particle colors" << std::endl;
    for (const auto& nDescendents : nDescendentOptions)
    {
        for (const auto& bdt : bdts)
        {
            const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
            std::cout << "DEBUG Point Z0" << std::endl;
            const auto outName = "WigglinessStudy1D_nDescendents" + std::to_string(nDescendents) + "_" + bdt;
            auto c = new TCanvas(outName.c_str(), "", 800, 600);
            c->SetLogx(xAxisLog);
            auto legend = (bdt != "trackScore") ? new TLegend(0.75, 0.6, 0.9, 0.9) : new TLegend(0.45, 0.6, 0.6, 0.9);

            float maxYTotal = 0;
            std::cout << "DEBUG Point Z1" << std::endl;
            for (const auto& type : types)
            {
                std::cout << "DEBUG Point Z2" << std::endl;
                auto h = histograms1D.at(nDescendents).at(bdt).at(type).get();
                h->SetStats(0); // Remove stats box
                const auto maxY = h->GetMaximum();
                if(maxY > maxYTotal)
                    maxYTotal = maxY;

                // Set the line color
                h->SetLineColor(particleColors.at(type));

                // Draw the histogram as a line
                h->Draw("HIST SAME");

                // Add an entry to the legend for this histogram
                const std::string typeName = type == "Other" ? "Other (including EXT)" : (type == "Beam On" ? "Beam On" : (type == "P" ? "Proton" : (type == "Mu" ? "Muon" : (type == "Pi" ? "Other Pion" : "Golden Pion"))));
                legend->AddEntry(h, typeName.c_str(), "f");
            }

            std::cout << "DEBUG Point Z3" << std::endl;

            // // Set y axis title
            // c->GetYaxis()->SetTitle("# Reconstructed Particles");

            // // Set title
            // const std::string bdtTitle = bdt == "muonBDTScore" ? "Muon BDT Score" : (bdt == "protonBDTScore" ? "Proton BDT Score" : (bdt == "goldenPionBDTScore" ? "Golden Pion BDT Score" : (bdt == "logBragg_pToMIP" ? "log(R_{p}/MIP)" : (bdt == "logBragg_piToMIP" ? "log(R_{#pi}/MIP)" : (bdt == "truncMeandEdx" ? "Truncated Mean dE/dx" : (bdt == "wiggliness" ? "Wiggliness" : (bdt == "trackScore" ? "Track Score" : "Number of Descendents")))))));
            // const std::string nDescendentTitle = "nDescendents " + std::to_string(nDescendents);
            // c->SetTitle((bdtTitle + " for " + nDescendentTitle + "\nReconstructed & Contained Beam-on and MC+EXT Particles in Events that Pass the CC #nu_#mu Preselection").c_str());

            legend->Draw();
            c->SaveAs(("plots/" + outName + ".pdf").c_str());
        }
    }
}