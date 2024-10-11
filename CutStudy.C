#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <map>
#include <string>
#include <memory>

void CutStudy() 
{
    ROOT::EnableImplicitMT();

    const std::vector versions = {"detVar_DecPelee2", // Full run 1-5 trained BDTs
                                    "detVar_DecPelee3", // Run 1 only trained BDTs
                                    "detVar_DecPelee3_withDetVar", // Run 1 only trained BDTs with additional detector variations (SCE & Recomb2)
                                    "detVar_DecPelee4", // Run 1 only trained BDTs
                                    "detVar_DecPelee4_withDetVar"}; // Run 1 only trained BDTs with additional detector variations (Recomb2 with weights scaled to overlay)
    const auto version = versions.at(4);

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
    // std::map<std::string, double> yAxisRange = {{"muonBDTScore", 0.0}, {"protonBDTScore", 0.0}, {"goldenPionBDTScore", 0.0}};

    // // map: name, pelee name
    // std::map<std::string, std::string> variables = {{"muonMomentum", "cc1pi_truth_muonMomentum"},
    //                                                 {"muonPionAngle", "cc1pi_truth_muonPionAngle"},
    //                                                 {"total", "passed_particleTrackScore"}}; // Use this as a placeholder as the values fall in the -1 to 1 range 

    std::map<std::string, std::string> cuts = {   
                                        {"passed_topologicalScoreCC", "event_cutValue_nTracks"},
                                        {"passed_min2Tracks", "event_cutValue_nUncontained"},
                                        {"passed_max1Uncontained", "event_cutValue_nNonProtons"},
                                        {"passed_2NonProtons", "event_cutValue_pionTruncatedMeandEdx"},
                                        {"passed_pionHasValiddEdx", "passed_pionNotInGap"},
                                        {"passed_pionNotInGap", "passed_muonNotInGap"},
                                        {"passed_muonNotInGap", "event_cutValue_topologicalScore"},
                                        {"passed_topologicalScore", "event_cutValue_maxVertexDist"},
                                        {"passed_startNearVertex", "event_cutValue_openingAngle"},
                                        {"passed_openingAngle", "event_cutValue_goldenPionBDT"},
                                        /*{"passed_likelyGoldenPion", "cutValueVariable"},*/};


    // map: file name, cut, variable, isSignal
    std::map<std::string, std::map<std::string, std::map<std::string, std::map<bool, std::unique_ptr<TH1D>>>>> histograms1DTruth;

    auto generateBinEdges = [](float start, float stop, int numBins) {
        std::vector<float> binEdges;
        float step = (stop - start) / numBins;
        for (int i = 0; i <= numBins; ++i) {
            binEdges.push_back(start + i * step);
        }
        return binEdges;
    };

    std::map<std::string, std::vector<float>> binningInfo = {
        {"event_cutValue_nTracks", generateBinEdges(0, 8, 8)},
        {"event_cutValue_nUncontained", generateBinEdges(0, 4, 4)},
        {"event_cutValue_nNonProtons", generateBinEdges(0, 8, 8)},
        {"event_cutValue_pionTruncatedMeandEdx", generateBinEdges(0.0f, 10.0f, 11)},
        {"passed_pionNotInGap", generateBinEdges(0, 1, 2)},
        {"passed_muonNotInGap", generateBinEdges(0, 1, 2)},
        {"event_cutValue_topologicalScore", generateBinEdges(0.0f, 2.0f, 11)}, 
        {"event_cutValue_maxVertexDist", generateBinEdges(0, 20, 11)},
        {"event_cutValue_openingAngle", generateBinEdges(0, 3.14f, 11)}, 
        {"event_cutValue_goldenPionBDT", generateBinEdges(-1.0f, 1.0f, 11)}, 
    };


    for (const auto& [name, filePath] : files) {
        for (const auto& [cut, cutVariable] : cuts) {
            auto& binEdges = binningInfo.at(cutVariable);
            int nBinsVar = binEdges.size() - 1;

            for (const auto& isSignal : {true, false})
            {
                std::string histName = cutVariable + "_" + name + "_" + cut + (isSignal ? "_signal" : "_background");
                std::string title = cutVariable + "_" + name + "_" + cut + (isSignal ? "_signal" : "_background");
                histograms1DTruth[name][cut][cutVariable].insert({isSignal, std::make_unique<TH1D>(histName.c_str(), title.c_str(), nBinsVar, binEdges.data())});
            }
        }
    }

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


        std::cout<< "DEBUG FillParticleHistogram Point 2.3" << std::endl;
        // Get the total number of entries
        const auto nEntries = tree->GetEntries();
        std::cout<< "DEBUG FillParticleHistogram Point 3" << std::endl;
        // Loop over the entries in the tree
        std::cout << "Processing file " << filePath << " with " << nEntries << " entries" << std::endl;
        for (Long64_t i = 0; i < nEntries; i++)
        {
            tree->GetEntry(i);

            // Update the progress bar at every percent
            if (i % (nEntries / 100) == 0)
            {
                const auto progress = static_cast<int>(100.0 * i / nEntries);
                std::cout << "\r[" << std::string(progress, '~') << std::string(100 - progress, ' ') << "] " << progress << "%" << std::flush;
            }

            // Check events pass at least the CC inclusive selection
            if(tree->GetLeaf("passed_topologicalScoreCC")->GetValue() && !tree->GetLeaf("isTrainingEvent")->GetValue())
            {
                const auto isSignal = tree->GetLeaf("cc1pi_signal")->GetValue();
                for (const auto& [cut, cutVariable] : cuts) {
                    if(tree->GetLeaf(cut.c_str())->GetValue()) // Add events that pass the cut to the histogram
                    {
                        const auto fillValue = tree->GetLeaf(cutVariable.c_str())->GetValue();
                        histograms1DTruth.at(name).at(cut).at(cutVariable).at(true)->Fill(fillValue); // todo temporaryily save everything in here
                    }
                } // End of loop over cuts
            } // End of CC inclusive selection
        } // End of loop over tree entries
    } // End of loop over files

    for (const auto& [name, filePath] : files) {
        if(name != "detVar_CVextra")
        {
            for (const auto& [cut, cutVariable] : cuts) {
                // ###### Signal events with true variable plots ######
                TCanvas canvas1D("canvas1D", "", 800, 500);
                TLegend legend(0.7, 0.7, 0.9, 0.9);

                auto& hist = histograms1DTruth.at(name).at(cut).at(cutVariable).at(true);
                hist->SetLineColor(kGreen);
                auto& histCV = histograms1DTruth.at("detVar_CVextra").at(cut).at(cutVariable).at(true);
                histCV->SetLineColor(kBlue);
                histCV->SetLineWidth(2);
                histCV->SetLineStyle(9);
                legend.AddEntry(hist.get(), name.c_str(), "l");
                legend.AddEntry(histCV.get(), "CV", "l");
                hist->Draw("Hist E");
                histCV->Draw("Hist E SAME");
                legend.Draw();
                hist->SetStats(0);
                const auto histMax = std::max(hist->GetMaximum(), histCV->GetMaximum());
                const auto histMin = std::min(hist->GetMinimum(), histCV->GetMinimum());
                hist->GetYaxis()->SetRangeUser(histMin * 0.8, histMax * 1.2);
                hist->GetXaxis()->SetTitle(cutVariable.c_str());
                hist->GetYaxis()->SetTitle("# Events");
                hist->SetTitle(("Events passing cut: " + cut).c_str());

                canvas1D.SaveAs(("plots/CutStudy_signal_" + cutVariable + "_" + cut + "_" + name + "_" + version + ".png").c_str());
                canvas1D.SaveAs(("plots/CutStudy_signal_" + cutVariable + "_" + cut + "_" + name + "_" + version + ".pdf").c_str());
                canvas1D.SaveAs(("plots/CutStudy_signal_" + cutVariable + "_" + cut + "_" + name + "_" + version + ".C").c_str());
            } // End of loop over cuts
        } // End of if statement for CVextra
    } // End of loop over files

    // for (const auto& cut : cuts) {
    //     TCanvas canvas("canvas", "", 500, 500);

    //     canvas.SetTitle(("Events passing cut: " + cut).c_str());
    //     TPaveText title(0.1, 0.95, 0.9, 0.99, "brNDC");
    //     title.AddText(("Events passing cut: " + cut).c_str());
    //     title.SetFillColor(0); // make the background of the title box transparent
    //     title.SetBorderSize(0); // remove the border
    //     title.SetTextSize(0.07); // change the text size
    //     title.Draw();

    //     THStack stack("stack", "");

    //     // Create two empty TH1D histograms with three bins for signal and background
    //     std::unique_ptr<TH1D> signalHist = std::make_unique<TH1D>("signalHist", "signalHist", 3, 0, 3);
    //     std::unique_ptr<TH1D> backgroundHist = std::make_unique<TH1D>("backgroundHist", "backgroundHist", 3, 0, 3);

    //     // Set fill colors
    //     signalHist->SetFillColor(kGreen);
    //     backgroundHist->SetFillColor(kGray);

    //     // Fill the three bins of each histograms with the three histograms for each detector variation
    //     for (auto b = 1; b <= 3; b++)
    //     {   
    //         const auto& [name, filePath] = files.at(b-1);
    //         signalHist->SetBinContent(b, histograms1DTruth.at(name).at(cut).at("total").at(true)->Integral());
    //         backgroundHist->SetBinContent(b, histograms1DTruth.at(name).at(cut).at("total").at(false)->Integral());
    //     }

    //     // Add histograms to the stack
    //     stack.Add(backgroundHist.get());
    //     stack.Add(signalHist.get());

    //     // Draw the stack
    //     stack.Draw("B");
    //     stack.Draw("E SAME");


    //     // Set labels for the bins
    //     stack.GetXaxis()->SetBinLabel(1, "SCE");
    //     stack.GetXaxis()->SetBinLabel(2, "CVExtra");
    //     stack.GetXaxis()->SetBinLabel(3, "Recomb2");

    //     // Set y-axis label
    //     stack.GetYaxis()->SetTitle("# Events");

    //     // Save the canvas
    //     canvas.SaveAs(("plots/CutStudy_signal_and_background_total_" + cut + "_" + version + ".png").c_str());
    //     canvas.SaveAs(("plots/CutStudy_signal_and_background_total_" + cut + "_" + version + ".pdf").c_str());
    //     canvas.SaveAs(("plots/CutStudy_signal_and_background_total_" + cut + "_" + version + ".C").c_str());
    // }

}