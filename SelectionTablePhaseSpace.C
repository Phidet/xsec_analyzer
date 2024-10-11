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

        // std::make_tuple("nu mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root", 0.13011),
        // std::make_tuple("nu mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root", 0.25750),
        // std::make_tuple("nu mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root", 0.20113),
        // std::make_tuple("nu mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi.root", 0.13074),
        // std::make_tuple("nu mc", "5",  rootPath + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root", 0.15196),

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

    std::vector<std::string> cuts {
        "cc1pi_selected_generic && cc1pi_reco_pionMomentum > 0.1",
        "cc1pi_selected_golden && cc1pi_reco_pionMomentum > 0.1",
    };

    std::vector<std::pair<bool, std::string>> cutsLatex {
        {true, "Generic selection"},
        {true, "Golden selection"}
    };

    // Map for all events; map: event type, cut, run
    std::map<std::string, std::map<std::string, std::map<std::string, float>>> eventCountMap;

    // Map for signal events; map: cut, run
    std::map<std::string, std::map<std::string, float>> signalCountMap, goldenCountMap;
    std::vector<std::string> runs = {"0", "1", "2", "3", "4bcd", "5"}; // 0 is the total of all runs

    for (const auto& cut : cuts) {
        for (const auto& run : runs) {
            eventCountMap["bnb"][cut][run] = 0.0f;
            eventCountMap["mc"][cut][run] = 0.0f;
            signalCountMap[cut][run] = 0.0f;
            goldenCountMap[cut][run] = 0.0f;
        }
    }

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

        std::cout << "\033[1;31mWARNING - Only processing 1% of events!!!\033[0m" << std::endl;
        const auto nEntries = tree->GetEntries()/100;
        for (Long64_t i = 0; i < nEntries; i++)
        {
            std::cout<<"\r"<<(100*i)/nEntries<<"%"<<std::flush;
            // Print progress bar
            // if (i % (nEntries / 100) == 0) { // Update progress bar for every 1% progress
            //     std::cout << "\r" // Move cursor to the beginning of the line
            //               << "[" // Start of progress bar
            //               << std::string(i / (nEntries / 100), '=') // Progress
            //               << std::string(100 - i / (nEntries / 100), ' ') // Remaining
            //               << "]" // End of progress bar
            //               << " " << i / (nEntries / 100) << "%" // Percentage
            //               << std::flush; // Flush output
            // }

            tree->GetEntry(i);

            // Calculate weight
            double weight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;
            // std::cout << "DEBUG entry for run: " << run << " of fileType: " << fileType << " with event weight: " << weight << " and file weight: " << fileWeight << " and total weight: " << weight*fileWeight << " and spline_weight: " << spline_weight << " and tuned_cv_weight: " << tuned_cv_weight << std::endl;
            weight *= fileWeight;

            for ( const auto& cut : cuts)
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
        std::ofstream outputFileLatex("eventCountTable_with_phase_space_run" + run + "_onlyTesting_lowPiMomThreshold.tex");

        // Write the header row
        std::string header = 
R"(
\begin{tabular}{lcccccccc}
\textbf{Cuts for Run%s} & \rotatebox{90}{\textbf{Signal}} & \rotatebox{90}{\textbf{Background}} & \rotatebox{90}{\textbf{Efficiency} ($E$)} & \rotatebox{90}{\textbf{Purity} ($P$)} & \rotatebox{90}{$E \times P$} & \rotatebox{90}{\textbf{Golden Fraction}} & \rotatebox{90}{\textbf{Data}} & \rotatebox{90}{\textbf{Data/MC Ratio}} \\
\midrule
\textbf{CC Inclusive Preselection} \\
\textbf{Particle-level Cuts} & & & & & & & &\\
\midrule
)";
        char buffer[512];
        std::string runStr = run == "0" ? "s 1-5" : " "+run;  
        sprintf(buffer, header.c_str(), runStr.c_str());
        outputFileLatex << buffer;

        // Write the data rows
        for (int i = 0; i < cuts.size(); i++)
        {
            const auto cut = cuts.at(i);
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

        //     if(cut == "passed_particleProtonChi2OverMuonChi2")
        //     {
        //         outputFileLatex << "\\midrule\n";
        //         outputFileLatex << "\\textbf{CC Inclusive Preselection}";
        //         outputFileLatex << " \\\\\n\\textbf{Event-level Cuts} & & & & & & & &";
        //         outputFileLatex << " \\\\\n\\midrule\n";
        //     } else if(cut == "passed_topologicalScoreCC")
        //     {
        //         outputFileLatex << "\\midrule\n";
        //         outputFileLatex << "\\textbf{Charged Pion Selection}";
        //         outputFileLatex << " \\\\\n\\midrule\n";
        //     } else if(cut == "passed_openingAngle")
        //     {
        //         outputFileLatex << "\\midrule\n";
        //         outputFileLatex << "\\textbf{Golden Pion Selection}";
        //         outputFileLatex << " \\\\\n\\midrule\n";
        //     }
        // }

        // End the tabular environment
        outputFileLatex << "\\end{tabular}\n";

        // Close the output file
        outputFileLatex.close();
    }

    std::cout<<"Done! Output written to eventCountTable_run*.tex"<<std::endl;
}