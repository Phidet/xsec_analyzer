#include <map>
#include <string>

void SelectionTable_GenieFSIBug() 
{
    const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/1March24/";

    // List of files; tuple: file type, run, path, fileWeight
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        std::make_tuple("nu mc", "1", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 1.0), // Ignoring proper scaling
        // FSI fix sample (treating it as run 2 here)
        std::make_tuple("nu mc", "2", "/exp/uboone/data/users/jdetje/ubcc1piPelee/pionFSIFix/GENIE_PionFSIFix_50_Percent_run1_ubcc1pi.root", 1.0),
    };

    // Vector of pairs; pair: stop when cut is failed, cut formula
    std::vector<std::pair<bool, std::string>> cuts {
        {true, "all"}, // <-- All
        {true, "passed_particleTrackScore"},
        {true, "passed_particleVertexDistance"},
        {true, "passed_particleGeneration"},
        {true, "passed_particleTrackLength"},
        {true, "passed_particleProtonChi2"},
        {true, "passed_particleMuonChi2"},
        {true, "passed_particleProtonChi2OverMuonChi2"},
        //{true, "passed_pandoraNuPDGIsNumu"},
        {true, "passed_daughterVerticesContained"},
        {true, "passed_nuVertexFiducial"},
        {true, "passed_topologicalOrFlashMatch"},
        //{true, "passed_topologicalScoreCC"},
        {true, "passed_min2Tracks"},
        {true, "passed_max1Uncontained"},
        {true, "passed_2NonProtons"},
        {true, "passed_pionHasValiddEdx"},
        {true, "passed_pionNotInGap"},
        {true, "passed_muonNotInGap"},
        {true, "passed_topologicalScore"},
        {true, "passed_startNearVertex"},
        {true, "passed_openingAngle"},
        {true, "passed_openingAngle && cc1pi_reco_muonMomentum > 0.15"},
        {true, "cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1"}, // isnt this the same as just cc1pi_selected_generic && cc1pi_reco_pionMomentum > 0.1
        // Contained Muon Selection
        {false, "cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained"},
        // Golden pion selection
        {true, "cc1pi_selected_golden && passed_likelyGoldenPion && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1"}, // Opening angle cut is already applied to passed_likelyGoldenPion
    };

    // Vector of pairs; pair: show in latex, cut name
    std::vector<std::pair<bool, std::string>> cutsLatex {
        {true, "all"},
        {true, "trackScore"},
        {true, "vertexDistance"},
        {true, "generation"},
        {true, "trackLength"},
        {true, "protonChi2"},
        {true, "muonChi2"},
        {true, "protonChi2OverMuonChi2"},
        // {false, "\\sout{pandoraNuPDGIsNumu}"},
        {true, "daughterVerticesContained"},
        {true, "nuVertexFiducial"},
        {true, "topologicalScore"},
        // {true, "topological\\sout{OrFlashMatch}"},
        // {false, "\\sout{topologicalScore}"},
        {true, "min2Tracks"},
        {true, "max1Uncontained"},
        {true, "2NonProtons"},
        {true, "pionHasValiddEdx"},
        {true, "pionNotInGap"},
        {true, "muonNotInGap"},
        {true, "topologicalScore"},
        {true, "startNearVertex"},
        {true, "openingAngle"},
        {true, "muonMomentum"},
        {true, "pionMomentum"},
        // Contained Muon Selection
        {true, "containedMuon"},
        // Golden pion selection
        {true, "likelyGoldenPion"},
    };

    // Map for all events; map: event type, cut, run
    std::map<std::string, std::map<std::string, std::map<std::string, double>>> eventCountMap;

    // Map for signal events; map: cut, run
    std::map<std::string, std::map<std::string, double>> signalCountMap, goldenCountMap;
    std::vector<std::string> runs = {"0", "1", "2", "3", "4bcd", "5"}; // 0 is the total of all runs

    for (const auto& [stopOnFailedCut, cut] : cuts) {
        for (const auto& run : runs) {
            eventCountMap["data"][cut][run] = 0.0f;
            eventCountMap["prediction"][cut][run] = 0.0f;
            signalCountMap[cut][run] = 0.0f;
            goldenCountMap[cut][run] = 0.0f;
        }
    }

    double totalMCSum0 = 0.f;
    // Loop over the files
    for (const auto &[fileType, run, path, fileWeight] : files)
    {
        std::cout << "Processing file " << path <<std::endl;
        // const auto weight = weights.at(f);
        // std::cout << "Processing file " << path << " with weight " << weight << std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(path.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");

        const std::string eventType = (fileType == "beam on") ? "data" : "prediction";
        const bool isMC = (fileType == "nu mc" || fileType == "dirt mc");

        // std::string weightString = (fileType == "prediction") ? "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)*" : "";
        // weightString += std::to_string(fileWeight);

        // Define variables to hold branch values
        Float_t spline_weight, tuned_cv_weight;
        Bool_t isCC1PiSignal, isTrueGoldenPion;
        Float_t cc1piTruthPionMomentum;

        // Set branch addresses
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);
        tree->SetBranchAddress("cc1pi_signal", &isCC1PiSignal);
        tree->SetBranchAddress("true_golden_cc1pi", &isTrueGoldenPion);
        tree->SetBranchAddress("cc1pi_truth_pionMomentum", &cc1piTruthPionMomentum);

        // std::cout << "\033[1;31mWARNING - Only processing 10% of events!!!\033[0m" << std::endl;
        const auto nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++)
        {
            if(i%(nEntries/100)==0) std::cout<<"\r"<<(100*i)/nEntries<<"%"<<std::flush;
            tree->GetEntry(i);

            for ( const auto& [stopOnFailedCut, cut] : cuts)
            {
                bool passed = false;
                if (cut == "all")
                {
                    passed = true;
                }
                else
                {
                    TTreeFormula formula("formula", cut.c_str(), tree);
                    passed = formula.EvalInstance();
                }
                
                if (passed)
                // if (cut == "all" || tree->GetLeaf(cut.c_str())->GetValue())
                {
                    // Calculate weight
                    double weight = (isMC && std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30) ? spline_weight*tuned_cv_weight : 1;
                    weight *= fileWeight;
                    eventCountMap.at(eventType).at(cut).at(run) += weight;
                    eventCountMap.at(eventType).at(cut).at("0") += weight;
                    if(isCC1PiSignal && cc1piTruthPionMomentum >= 0.1)
                    {
                        signalCountMap.at(cut).at(run) += weight;
                        signalCountMap.at(cut).at("0") += weight;
                        if(isTrueGoldenPion)
                        {
                            goldenCountMap.at(cut).at(run) += weight;
                            goldenCountMap.at(cut).at("0") += weight;
                        }
                    }
                }
                else
                {
                    if(stopOnFailedCut)
                        break;
                }
            } // End of loop over cuts
        } // End of loop over events
    } // End of loop over files

    std::cout<<std::endl;

    // *********************************
    // Save the table again as latex
    // *********************************
    for (const auto& run : runs) {
        // Open the output file
        std::ofstream outputFileLatex("eventCountTable_run" + run + "_onlyTesting_lowPiMomThreshold_fixed_muonContainedOption_GENIEFSIBugFix.tex");

        // Write the header row
        std::string header = 
R"(
\begin{tabular}{lrrrrrrrr}
\large{\textbf{Generic Selection}} \vspace{1mm}\\
\textbf{Cuts for Run%s} & \rot{\textbf{Signal}} & \rot{\textbf{Background}} & \rot{\textbf{Efficiency} ($E$)} & \rot{\textbf{Purity} ($P$)} & \rot{$E \times P$} & \rot{\textbf{Golden Fraction}} & \rot{\textbf{Data}} & \rot{\textbf{Ratio}} \\
\midrule
\textbf{Particle-level Preselection} \\
\midrule
)";
        char buffer[512];
        std::string runStr = run == "0" ? "s 1-5" : " "+run;
        sprintf(buffer, header.c_str(), runStr.c_str());
        outputFileLatex << buffer;

        // Write the data rows
        for (int i = 0; i < cuts.size(); i++)
        {
            const auto cut = cuts.at(i).second;
            const auto cutLatexLabel = cutsLatex.at(i).second;
            const auto cutLatexShow = cutsLatex.at(i).first;
            const auto dataMCRatio = eventCountMap.at("prediction").at(cut).at(run) != 0 ? eventCountMap.at("data").at(cut).at(run)/eventCountMap.at("prediction").at(cut).at(run) : 0;
            const auto efficiency = signalCountMap.at(cut).at(run) / signalCountMap.at("all").at(run);
            const auto purity = signalCountMap.at(cut).at(run) / eventCountMap.at("prediction").at(cut).at(run);
            const auto effpur =  efficiency*purity;
            const auto goldenFraction = goldenCountMap.at(cut).at(run) / signalCountMap.at(cut).at(run);
            outputFileLatex  << cutLatexLabel;

            if(cutLatexShow)
            {
                outputFileLatex << " & " << std::round(10*signalCountMap.at(cut).at(run))/10.0 // Round to 1 decimal place
                << " & " << std::round(10*(eventCountMap.at("prediction").at(cut).at(run) - signalCountMap.at(cut).at(run)))/10.0 // Round to 1 decimal place
                << " & " << std::round(1000*efficiency)/10.0 << "\\%" // Convert to percentage and round to 1 decimal place
                << " & " << std::round(1000*purity)/10.0 << "\\%" // Convert to percentage and round to 1 decimal place
                << " & " << std::round(1000*effpur)/1000.0 // Round to 3 decimal places
                << " & " << std::round(1000*goldenFraction)/10.0 << "\\%" // Convert to percentage and round to 1 decimal place
                << " & " << eventCountMap.at("data").at(cut).at(run)
                << " & " << std::round(100*dataMCRatio)/100.0; // Round to 2 decimal places
                outputFileLatex << " \\\\\n";
                // outputFileLatex << "\\\\\n\\midrule\n";
            }
            else
            {
                outputFileLatex << " & & & & & & & & \\\\\n";
            }

            if(cut == "passed_particleProtonChi2OverMuonChi2")
            {
                outputFileLatex << "\\midrule\n";
                outputFileLatex << "\\textbf{Event-level Preselection}";
                outputFileLatex << " \\\\\n\\midrule\n";
            } else if(cut == "passed_topologicalOrFlashMatch")
            // else if(cut == "passed_topologicalScoreCC")
            {
                outputFileLatex << "\\midrule\n";
                outputFileLatex << "\\textbf{Charged Pion Selection}";
                outputFileLatex << " \\\\\n\\midrule\n";
            } else if(cut == "passed_startNearVertex")
            {
                outputFileLatex << "\\midrule\n";
                outputFileLatex << "\\textbf{Phase Space Cuts}";
                outputFileLatex << " \\\\\n\\midrule\n";
            } else if(cut == "cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1")
            {
                outputFileLatex << "\\\\ [3mm] \\midrule\n"; // Add 5mm of vertical space and a midrule
                outputFileLatex << "\\textbf{Subset 1: Muon Momentum}\\\\\n";
                outputFileLatex << "Cut applied to generic selection";
                outputFileLatex << " \\\\\\midrule\n";
            } else if(cut == "cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained")
            {
                outputFileLatex << "\\\\ [1mm] \\midrule\n"; // Add 5mm of vertical space and a midrule
                outputFileLatex << "\\textbf{Subset 2: Pion Momentum}\\\\\n";
                outputFileLatex << "Cut applied to generic selection";
                outputFileLatex << " \\\\\\midrule\n";
            }


            // else if(cut == "passed_openingAngle")
            // {
            //     outputFileLatex << "\\midrule\n";
            //     outputFileLatex << "\\textbf{Golden Pion Selection}";
            //     outputFileLatex << " \\\\\n\\midrule\n";
            // }
        }

        // End the tabular environment
        outputFileLatex << "\\end{tabular}\n";

        // Close the output file
        outputFileLatex.close();
    }

    // Print out eventCountMap.at("prediction").at(cut).at(run) and signalCountMap.at(cut).at(run) for all runs at cut "all"
    double totalMCSum = 0;
    double totalSignalSum = 0;
    for (const auto& run : runs) {
        std::cout << "Event count for run " << run << " at cut \"all\": " << eventCountMap.at("prediction").at("all").at(run) << std::endl;
        std::cout << "Signal count for run " << run << " at cut \"all\": " << signalCountMap.at("all").at(run) << std::endl;
        if (run != "0") {
            totalMCSum += eventCountMap.at("prediction").at("all").at(run);
            totalSignalSum += signalCountMap.at("all").at(run);
        }
    }
    std::cout << "Total MC sum runs 1-5: " << totalMCSum << std::endl;
    std::cout << "Total Signal sum runs 1-5: " << totalSignalSum << std::endl;

    std::cout<<"Done! Output written to eventCountTable_run*_GENIEFSIBugFix.tex"<<std::endl;
}