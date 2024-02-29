#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <map>
#include <string>
#include <memory>

#include "EventCategory.hh" 

void BDTStudy1DVsData() 
{
    const EventCategoryInterpreter& eventCategoryInterpreter = EventCategoryInterpreter::Instance();

    const std::string rootPath = "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/22Feb24/";
    // tuple: type, run, file path, run weight
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

        std::make_tuple("nu mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root", 0.13011),
        std::make_tuple("nu mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root", 0.25750),
        std::make_tuple("nu mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root", 0.20113),
        std::make_tuple("nu mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi.root", 0.13074),
        std::make_tuple("nu mc", "5",  rootPath + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root", 0.15196),

        std::make_tuple("dirt mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root", 0.52806),
        std::make_tuple("dirt mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root", 0.27521),
        std::make_tuple("dirt mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root", 0.80892),
        std::make_tuple("dirt mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi.root", 0.39701),
        std::make_tuple("dirt mc", "5",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_ubcc1pi.root", 0.41280),
    };

    const std::vector<std::string> runs {"0", "1", "2", "3", "4bcd", "5"}; // Here 0 is the full set of all runs
    const std::vector<std::string> bdts = {"muonBDTScore", "protonBDTScore", "goldenPionBDTScore", "logBragg_pToMIP", "logBragg_piToMIP", "truncMeandEdx", "wiggliness", "trackScore", "nDescendents"};
    // const std::vector<std::string> types = {"EXT", "Dirt", "E", "Photon", "K", "P", "Mu", "Pi", "Golden Pi", "Beam On"};
    const std::vector<std::string> types = {"Other", "P", "Mu", "Pi", "Golden Pi", "Beam On"};
    const std::vector<bool> dataOptions = {true, false};

    // map: variable, tuple(nBins, xMin, xMax, axisLog)
    const std::map<std::string, std::tuple<int, double, double, bool>> binningInfo = {
        {"muonBDTScore", std::make_tuple(60, -1, 1, false)},
        {"protonBDTScore", std::make_tuple(60, -1, 1, false)},
        {"goldenPionBDTScore", std::make_tuple(60, -1, 1, false)},
        {"logBragg_pToMIP", std::make_tuple(60, -8, 8, false)},
        {"logBragg_piToMIP", std::make_tuple(60, -4, 7, false)},
        {"truncMeandEdx", std::make_tuple(60, 0, 10.0, false)},
        {"wiggliness", std::make_tuple(60, 0.00005, 0.25, true)},
        {"trackScore", std::make_tuple(60, 0, 1.0, false)},
        {"nDescendents", std::make_tuple(4, 0, 4, false)}
    };

    // map: run, bdt, particle
    std::map<std::string, std::map<std::string, std::map<std::string, std::unique_ptr<TH1D>>>> histograms1D;
    for (const auto& run : runs) {
        for (const auto& bdt : bdts) {
            const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
            for (const auto& type : types) {
                    const auto name = "h_" + run + "_" + bdt + "_" + type;
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
                        histograms1D[run][bdt][type] = std::make_unique<TH1D>((name + "_" + bdt + "_" + type).c_str(), (name + " reco Particles Backtracked to True " + type + " in Signal Events; " + bdt + " BDT Score").c_str(), nBinsBDT, binEdges.data());
                    }
                    else
                    {
                        histograms1D[run][bdt][type] = std::make_unique<TH1D>((name + "_" + bdt + "_" + type).c_str(), (name + " reco Particles Backtracked to True " + type + " in Signal Events; " + bdt + " BDT Score").c_str(), nBinsBDT, xMin, xMax);
                    }
                    histograms1D[run][bdt][type]->Sumw2();
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

        std::vector<float>* pReco_particle_ccinc_muonBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_protonBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_goldenPionBDTScore = nullptr;
        std::vector<int>* pReco_particle_ccinc_generation = nullptr;
        std::vector<int>* pReco_particle_ccinc_backtracked_pdg = nullptr;
        std::vector<bool>* pReco_particle_ccinc_backtracked_goldenPion = nullptr;
        std::vector<bool>* pReco_particle_ccinc_isContained = nullptr;

        std::vector<float>* pReco_particle_ccinc_logBragg_pToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_logBragg_piToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_truncMeandEdx = nullptr;
        std::vector<float>* pReco_particle_ccinc_wiggliness = nullptr;
        std::vector<float>* pReco_particle_ccinc_trackScore = nullptr;
        std::vector<int>* pReco_particle_ccinc_nDescendents = nullptr;

        tree->SetBranchAddress("reco_particle_ccinc_muonBDTScore", &pReco_particle_ccinc_muonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_protonBDTScore", &pReco_particle_ccinc_protonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_goldenPionBDTScore", &pReco_particle_ccinc_goldenPionBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_generation", &pReco_particle_ccinc_generation);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_pdg", &pReco_particle_ccinc_backtracked_pdg);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_goldenPion", &pReco_particle_ccinc_backtracked_goldenPion);
        tree->SetBranchAddress("reco_particle_ccinc_contained", &pReco_particle_ccinc_isContained);

        tree->SetBranchAddress("reco_particle_ccinc_logBragg_pToMIP", &pReco_particle_ccinc_logBragg_pToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_logBragg_piToMIP", &pReco_particle_ccinc_logBragg_piToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_truncMeandEdx", &pReco_particle_ccinc_truncMeandEdx);
        tree->SetBranchAddress("reco_particle_ccinc_wiggliness", &pReco_particle_ccinc_wiggliness);
        tree->SetBranchAddress("reco_particle_ccinc_trackScore", &pReco_particle_ccinc_trackScore);
        tree->SetBranchAddress("reco_particle_ccinc_nDescendents", &pReco_particle_ccinc_nDescendents);

        Float_t spline_weight, tuned_cv_weight;
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);

        // std::string metaType = isBeamOn ? "Beam On" : (isDirt ? "Dirt" : (isBeamOff ? "EXT" : ""));
        std::string metaType = isBeamOn ? "Beam On" : ((isDirt || isBeamOff) ? "Other" : "");

        // Get the total number of entries
        // std::cout<<"WARNING - Only using 3% of the entries!!!!!!!!!"<<std::endl;
        const auto nEntries = tree->GetEntries();

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

                // if(isNuMC)
                //     eventWeight *= 2; // Compensate for skipping the BDT training events and only using 50% of the sample

                // Loop over the values in the vector and fill the histogram
                for (Long64_t v = 0; v< pReco_particle_ccinc_generation->size(); v++)
                {
                    // Set the type to the type determined by the file first
                    std::string type = metaType; // Only when that is not set, use the particle type
                    

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
                                case 13:
                                    type = "Mu";
                                    break;
                                case 211:
                                {
                                    const auto isGolden = pReco_particle_ccinc_backtracked_goldenPion->at(v);
                                    type = isGolden ? "Golden Pi" : "Pi";
                                    break;
                                }
                                case 2212:
                                    type = "P";
                                    break;
                                default:
                                    type = "Other";


                                // case 11:
                                //     type = "E";
                                //     break;
                                // case 111: // pi0 are treated as photons due to quick decay
                                // case 22:
                                // {
                                //     std::cout << "DEBUG: photon PDG: " << pdg << " with weight: " << eventWeight << std::endl;
                                //     type = "Photon";
                                //     break;
                                // }
                                // case 321:
                                //     type = "K";
                                //     break;
                                // case 2147483647: // Could not find a backtracked pdg
                                // case 0:
                                //     type = "EXT";
                                //     break;
                                // default:
                                // {
                                //     std::cout << "DEBUG: switch default PDG: " << pdg << std::endl;
                                // }
                            }
                        }
                        
                        if (type != "")
                        {
                            // for (const auto& bdt : bdts) 
                            // {
                            //     float bdtScore;
                            //     bdtScore = bdt == "muonBDTScore" ? pReco_particle_ccinc_muonBDTScore->at(v) : (bdt == "protonBDTScore" ? pReco_particle_ccinc_protonBDTScore->at(v) : pReco_particle_ccinc_goldenPionBDTScore->at(v));
                            //     histograms1D.at(run).at(bdt).at(type)->Fill(bdtScore, eventWeight); // Fill the hist for this run
                            //     histograms1D.at("0").at(bdt).at(type)->Fill(bdtScore, eventWeight); // Fill the hist for the combined run = "0"
                            // } // End of loop over BDTs

                            // // Skip if muonBDT is not in the range // replacement for hasFeatures
                            // if(pReco_particle_ccinc_muonBDTScore->at(v) < -1 || pReco_particle_ccinc_muonBDTScore->at(v) > 1)
                            //     continue;
                            // if(pReco_particle_ccinc_protonBDTScore->at(v) < -1 || pReco_particle_ccinc_protonBDTScore->at(v) > 1)
                            //     continue;
                            // if(pReco_particle_ccinc_goldenPionBDTScore->at(v) < -1 || pReco_particle_ccinc_goldenPionBDTScore->at(v) > 1)
                            //     continue;

                            histograms1D.at(run).at("muonBDTScore").at(type)->Fill(pReco_particle_ccinc_muonBDTScore->at(v), eventWeight);
                            histograms1D.at(run).at("protonBDTScore").at(type)->Fill(pReco_particle_ccinc_protonBDTScore->at(v), eventWeight);
                            histograms1D.at(run).at("goldenPionBDTScore").at(type)->Fill(pReco_particle_ccinc_goldenPionBDTScore->at(v), eventWeight);

                            histograms1D.at("0").at("muonBDTScore").at(type)->Fill(pReco_particle_ccinc_muonBDTScore->at(v), eventWeight);
                            histograms1D.at("0").at("protonBDTScore").at(type)->Fill(pReco_particle_ccinc_protonBDTScore->at(v), eventWeight);
                            histograms1D.at("0").at("goldenPionBDTScore").at(type)->Fill(pReco_particle_ccinc_goldenPionBDTScore->at(v), eventWeight);

                            histograms1D.at(run).at("logBragg_pToMIP").at(type)->Fill(pReco_particle_ccinc_logBragg_pToMIP->at(v), eventWeight);
                            histograms1D.at(run).at("logBragg_piToMIP").at(type)->Fill(pReco_particle_ccinc_logBragg_piToMIP->at(v), eventWeight);
                            histograms1D.at(run).at("truncMeandEdx").at(type)->Fill(pReco_particle_ccinc_truncMeandEdx->at(v), eventWeight);
                            histograms1D.at(run).at("wiggliness").at(type)->Fill(pReco_particle_ccinc_wiggliness->at(v), eventWeight);
                            histograms1D.at(run).at("trackScore").at(type)->Fill(pReco_particle_ccinc_trackScore->at(v), eventWeight);
                            histograms1D.at(run).at("nDescendents").at(type)->Fill(pReco_particle_ccinc_nDescendents->at(v), eventWeight);

                            histograms1D.at("0").at("logBragg_pToMIP").at(type)->Fill(pReco_particle_ccinc_logBragg_pToMIP->at(v), eventWeight);
                            histograms1D.at("0").at("logBragg_piToMIP").at(type)->Fill(pReco_particle_ccinc_logBragg_piToMIP->at(v), eventWeight);
                            histograms1D.at("0").at("truncMeandEdx").at(type)->Fill(pReco_particle_ccinc_truncMeandEdx->at(v), eventWeight);
                            histograms1D.at("0").at("wiggliness").at(type)->Fill(pReco_particle_ccinc_wiggliness->at(v), eventWeight);
                            histograms1D.at("0").at("trackScore").at(type)->Fill(pReco_particle_ccinc_trackScore->at(v), eventWeight);
                            histograms1D.at("0").at("nDescendents").at(type)->Fill(pReco_particle_ccinc_nDescendents->at(v), eventWeight);
                        }
                        else
                        {
                            throw std::runtime_error("ERROR: type is empty");
                        }
                    } // End of if contained and generation 2
                } // End of loop over vector entries (aka types)
            } // End of if signal
        } // End of loop over tree entries
    } // End of loop over files


    std::cout << "DEBUG: Finished filling histograms" << std::endl;

    // Create a particle to color map
    std::map<std::string, int> particleColors = {
        {"P", kOrange+1},
        {"Mu", kBlue},
        {"Pi", kMagenta},
        // {"E", kCyan},
        // {"Photon", kYellow},
        // {"K", kMagenta},
        {"Golden Pi", kGreen},
        // {"EXT", kBlack},
        // {"Dirt", kOrange-6},
        {"Beam On", kBlack},
        {"Other", kGray}
    };


    std::cout << "DEBUG: Finished setting particle colors" << std::endl;
    for (const auto& run : runs)
    {
        for (const auto& bdt : bdts)
        {
            const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
            std::cout << "DEBUG Point Z0" << std::endl;
            const auto outName = "BDTStudy1DVsData_run" + run + "_" + bdt;
            auto c = new TCanvas(outName.c_str(), "", 800, 600);
            c->SetLogx(xAxisLog);
            auto legend = (bdt != "trackScore") ? new TLegend(0.75, 0.6, 0.9, 0.9) : new TLegend(0.45, 0.6, 0.6, 0.9);

            // Create two pads
            TPad *pad1 = new TPad((outName+"pad1").c_str(), "pad1", 0, 0.3, 1, 0.95);
            pad1->SetBottomMargin(0.05); // Upper pad, no space on the bottom
            pad1->Draw();             // Draw the upper pad
            pad1->cd();               // pad1 becomes the current pad
            gPad->SetLogx(xAxisLog);

            // Create a THStack
            auto hs = new THStack(outName.c_str(), "");

            // Create an initially empty histogram
            TH1D *sum = (TH1D*)histograms1D.at(run).at(bdt).at(types[0]).get()->Clone();
            sum->Reset(); // Reset all bins to 0

            float maxYTotal = 0;
            std::cout << "DEBUG Point Z1" << std::endl;
            for (const auto& type : types)
            {
                std::cout << "DEBUG Point Z2" << std::endl;
                auto h = histograms1D.at(run).at(bdt).at(type).get();
                h->SetStats(0); // Remove stats box
                const auto maxY = h->GetMaximum();
                if(maxY > maxYTotal)
                    maxYTotal = maxY;
                if(type == "Beam On")
                {
                    legend->AddEntry(h, "Beam On", "l");
                    continue;
                }

                if(type == "Other")//type == "EXT")
                    eventCategoryInterpreter.set_ext_histogram_style(h);
                else
                    h->SetFillColor(particleColors[type]);

                h->SetLineColor(kWhite); // Add this line to remove the outline
                h->SetLineWidth(0); // Add this line to remove the outline
                // Add the histogram to the stack
                hs->Add(h);

                // Add the histogram to the sum
                sum->Add(h);

                // Add an entry to the legend for this histogram
                const std::string typeName = type == "Other" ? "Other (including EXT)" : (type == "Beam On" ? "Beam On" : (type == "P" ? "Proton" : (type == "Mu" ? "Muon" : (type == "Pi" ? "Other Pion" : "Golden Pion"))));
                legend->AddEntry(h, typeName.c_str(), "f");
            }

            std::cout << "DEBUG Point Z3" << std::endl;

            hs->SetMaximum(maxYTotal*1.2);
            
            hs->Draw("HIST");
            // auto hsCopy = (THStack*)hs->Clone();
            sum->SetFillColor(TColor::GetColorTransparent(kBlack, 0.5));
            sum->Draw("E2 SAME");

            // Get the beam on histogram
            auto h = histograms1D.at(run).at(bdt).at("Beam On").get();
            eventCategoryInterpreter.set_bnb_data_histogram_style(h);
            h->SetStats(0); // Remove stats box
            h->Draw("E same");

            legend->Draw();

            // Set y axis title
            hs->GetYaxis()->SetTitle("# Reconstructed Particles");

            // Set title
            const std::string bdtTitle = bdt == "muonBDTScore" ? "Muon BDT Score" : (bdt == "protonBDTScore" ? "Proton BDT Score" : (bdt == "goldenPionBDTScore" ? "Golden Pion BDT Score" : (bdt == "logBragg_pToMIP" ? "log(R_{p}/MIP)" : (bdt == "logBragg_piToMIP" ? "log(R_{#pi}/MIP)" : (bdt == "truncMeandEdx" ? "Truncated Mean dE/dx" : (bdt == "wiggliness" ? "Wiggliness" : (bdt == "trackScore" ? "Track Score" : "Number of Descendents")))))));
            const std::string runTitle = run == "0" ? "All Runs" : ("Run " + run);
            hs->SetTitle((bdtTitle + " for " + runTitle + "\nReconstructed & Contained Beam-on and MC+EXT Particles in Events that Pass the CC #nu_#mu Preselection").c_str());
            
            // Draw the ratio plot in the lower pad
            c->cd(); // Go back to the main canvas before creating a new pad
            TPad *pad2 = new TPad((outName+"pad2").c_str(), "pad2", 0, 0.05, 1, 0.25);
            pad2->SetTopMargin(0);
            pad2->SetBottomMargin(0.3);
            pad2->Draw();
            pad2->cd(); // pad2 becomes the current pad
            gPad->SetLogx(xAxisLog);

            // Create the ratio plot
            TH1F *ratio = (TH1F*)h->Clone("ratio");

            for (int i = 1; i <= ratio->GetNbinsX(); i++)
            {
                std::cout << "Bin " << i << " orig ratio: " << ratio->GetBinContent(i) <<std::endl;
            }

            for (int i = 1; i <= sum->GetNbinsX(); i++)
            {
                std::cout << "Bin " << i << " sum: " << sum->GetBinContent(i) <<std::endl;
            }
            std::cout << "DEBUG Point Z4" << std::endl;

            ratio->SetLineColor(kBlack);
            ratio->Sumw2();
            ratio->Divide(sum);
            // Get the max and min y values and set the range
            double maxY = 1.2*ratio->GetMaximum();
            double minY = 0.8*ratio->GetMinimum();
            ratio->SetMinimum(minY);
            ratio->SetMaximum(maxY);

            // Adjust the font size of the x-axis labels and title
            ratio->GetXaxis()->SetTitleSize(0.12);
            ratio->GetXaxis()->SetLabelSize(0.12);

            // Adjust the font size of the y-axis labels and title
            ratio->GetYaxis()->SetTitleSize(0.12);
            ratio->GetYaxis()->SetLabelSize(0.12);

            ratio->Draw("E"); // Draw the ratio plot

            // Add a horizontal line at y=1 to the ratio plot (Bottom pad)
            TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 1, ratio->GetXaxis()->GetXmax(), 1);
            line->SetLineColor(kBlack);
            line->SetLineWidth(2);
            line->Draw();

            // No title for the ratio plot
            ratio->SetTitle("");

            // Set a nicely formatted y-axis title for the ratio plot
            ratio->GetYaxis()->SetTitle("Beam On / (MC + EXT)");
            // Set a nicely formatted x-axis title for the ratio plot
            ratio->GetXaxis()->SetTitle(bdtTitle.c_str());

            // Print out the bin values for the ratio plot
            for (int i = 1; i <= ratio->GetNbinsX(); i++)
            {
                std::cout << "Bin " << i << " ratio: " << ratio->GetBinContent(i) <<std::endl;
            }
            std::cout << "DEBUG Point Z5" << std::endl;

            c->SaveAs(("plots/" + outName + ".pdf").c_str());
            // c->SaveAs(("plots2/" + outName + ".png").c_str());
            // c->SaveAs(("plots2/" + outName + ".C").c_str());
        }
    }
}