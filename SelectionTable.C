#include <map>
#include <string>

void SelectionTable() 
{
    const std::string rootPath = "/pnfs/uboone/persistent/users/jdetje/ubcc1piPelee/22Feb24/";

    // List of files; tuple: file type, run, path, fileWeight
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

        std::make_tuple("nu mc", "1",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
        std::make_tuple("nu mc", "2",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root", 0.25750*2.0),
        std::make_tuple("nu mc", "3",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root", 0.20113*2.0),
        std::make_tuple("nu mc", "4bcd", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root", 0.13074*2.0),
        std::make_tuple("nu mc", "5",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196*2.0),

        std::make_tuple("dirt mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root", 0.52806),
        std::make_tuple("dirt mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root", 0.27521),
        std::make_tuple("dirt mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root", 0.80892),
        std::make_tuple("dirt mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi.root", 0.39701),
        std::make_tuple("dirt mc", "5",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_ubcc1pi.root", 0.41280),
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
        {true, "passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1"},
        // Contained Muon Selection
        {false, "passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained"},
        // Golden pion selection
        {true, "passed_likelyGoldenPion && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1"}, // Opening angle cut is already applied to passed_likelyGoldenPion
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
            eventCountMap["bnb"][cut][run] = 0.0f;
            eventCountMap["mc"][cut][run] = 0.0f;
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

        const std::string eventType = (fileType == "beam on") ? "bnb" : "mc";

        const std::string weightString = "(std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1)*" + std::to_string(fileWeight);

        // Define variables to hold branch values
        Float_t spline_weight, tuned_cv_weight;
        Bool_t isTrueCC1pi, isTrueGoldenPion;

        // Set branch addresses
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);
        tree->SetBranchAddress("true_cc1pi", &isTrueCC1pi);
        tree->SetBranchAddress("true_golden_cc1pi", &isTrueGoldenPion);

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
                    double weight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;
                    weight *= fileWeight;
                    eventCountMap.at(eventType).at(cut).at(run) += weight;
                    eventCountMap.at(eventType).at(cut).at("0") += weight;
                    if(isTrueCC1pi)
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
        // std::cout<<"\n 0---->"<<eventCountMap.at("mc").at("all").at("0")<<std::endl;
        // std::cout<<" run---->"<<eventCountMap.at("mc").at("all").at(run)<<std::endl;
        // std::cout<<" diff---->"<<eventCountMap.at("mc").at("all").at("0")-eventCountMap.at("mc").at("all").at(run)<<std::endl;
        // if(eventType == "mc") totalMCSum0 += eventCountMap.at("mc").at("all").at(run);
        // std::cout<<" tot---->"<<totalMCSum0<<"\n"<<std::endl;
    } // End of loop over files

    std::cout<<std::endl;

    // *********************************
    // Save the table again as latex
    // *********************************
    for (const auto& run : runs) {
        // Open the output file
        std::ofstream outputFileLatex("eventCountTable_run" + run + "_onlyTesting_lowPiMomThreshold_fixed_muonContainedOption.tex");

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
            const auto dataMCRatio = eventCountMap.at("mc").at(cut).at(run) != 0 ? eventCountMap.at("bnb").at(cut).at(run)/eventCountMap.at("mc").at(cut).at(run) : 0;
            const auto efficiency = signalCountMap.at(cut).at(run) / signalCountMap.at("all").at(run);
            const auto purity = signalCountMap.at(cut).at(run) / eventCountMap.at("mc").at(cut).at(run);
            const auto effpur =  efficiency*purity;
            const auto goldenFraction = goldenCountMap.at(cut).at(run) / signalCountMap.at(cut).at(run);
            outputFileLatex  << cutLatexLabel;

            if(cutLatexShow)
            {
                outputFileLatex << " & " << std::round(10*signalCountMap.at(cut).at(run))/10.0 // Round to 1 decimal place
                << " & " << std::round(10*(eventCountMap.at("mc").at(cut).at(run) - signalCountMap.at(cut).at(run)))/10.0 // Round to 1 decimal place
                << " & " << std::round(1000*efficiency)/10.0 << "\\%" // Convert to percentage and round to 1 decimal place
                << " & " << std::round(1000*purity)/10.0 << "\\%" // Convert to percentage and round to 1 decimal place
                << " & " << std::round(1000*effpur)/1000.0 // Round to 3 decimal places
                << " & " << std::round(1000*goldenFraction)/10.0 << "\\%" // Convert to percentage and round to 1 decimal place
                << " & " << eventCountMap.at("bnb").at(cut).at(run)
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
            } else if(cut == "passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1")
            {
                outputFileLatex << "\\\\ [3mm] \\midrule\n"; // Add 5mm of vertical space and a midrule
                outputFileLatex << "\\textbf{Subset 1: Muon Momentum}\\\\\n";
                outputFileLatex << "Cut applied to generic selection";
                outputFileLatex << " \\\\\\midrule\n";
            } else if(cut == "passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained")
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

    // Print out eventCountMap.at("mc").at(cut).at(run) and signalCountMap.at(cut).at(run) for all runs at cut "all"
    double totalMCSum = 0;
    double totalSignalSum = 0;
    for (const auto& run : runs) {
        std::cout << "Event count for run " << run << " at cut \"all\": " << eventCountMap.at("mc").at("all").at(run) << std::endl;
        std::cout << "Signal count for run " << run << " at cut \"all\": " << signalCountMap.at("all").at(run) << std::endl;
        if (run != "0") {
            totalMCSum += eventCountMap.at("mc").at("all").at(run);
            totalSignalSum += signalCountMap.at("all").at(run);
        }
    }
    std::cout << "Total MC sum runs 1-5: " << totalMCSum << std::endl;
    std::cout << "Total Signal sum runs 1-5: " << totalSignalSum << std::endl;

    std::cout<<"Done! Output written to eventCountTable_run*.tex"<<std::endl;
}