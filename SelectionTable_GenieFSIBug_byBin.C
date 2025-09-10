#include <map>
#include <string>

void SelectionTable_GenieFSIBug_byBin() 
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
        {false, "all"}, // <-- Alll
        // General Selection
        // Total
        {true, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1"}, 
        
        // Muon cos theta
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= -1 && cc1pi_truth_muonCosTheta < -0.27    &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= -1 && cc1pi_reco_muonCosTheta < -0.27"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= -0.27 && cc1pi_truth_muonCosTheta < 0.29  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= -0.27 && cc1pi_reco_muonCosTheta < 0.29"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.29 && cc1pi_truth_muonCosTheta < 0.46   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= 0.29 && cc1pi_reco_muonCosTheta < 0.46"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.46 && cc1pi_truth_muonCosTheta < 0.58   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= 0.46 && cc1pi_reco_muonCosTheta < 0.58"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.58 && cc1pi_truth_muonCosTheta < 0.67   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= 0.58 && cc1pi_reco_muonCosTheta < 0.67"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.67 && cc1pi_truth_muonCosTheta < 0.77   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= 0.67 && cc1pi_reco_muonCosTheta < 0.77"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.77 && cc1pi_truth_muonCosTheta < 0.82   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= 0.77 && cc1pi_reco_muonCosTheta < 0.82"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.82 && cc1pi_truth_muonCosTheta < 0.88   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= 0.82 && cc1pi_reco_muonCosTheta < 0.88"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.88 && cc1pi_truth_muonCosTheta < 0.93   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= 0.88 && cc1pi_reco_muonCosTheta < 0.93"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.93 && cc1pi_truth_muonCosTheta < 0.97   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= 0.93 && cc1pi_reco_muonCosTheta < 0.97"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.97 && cc1pi_truth_muonCosTheta <= 1     &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonCosTheta >= 0.97 && cc1pi_reco_muonCosTheta <= 1"},

        // Pion cos theta
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= -1 && cc1pi_truth_pionCosTheta < -0.47   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionCosTheta >= -1    && cc1pi_reco_pionCosTheta < -0.47"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= -0.47 && cc1pi_truth_pionCosTheta < 0    &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionCosTheta >= -0.47 && cc1pi_reco_pionCosTheta < 0"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0 && cc1pi_truth_pionCosTheta < 0.39     &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionCosTheta >= 0     && cc1pi_reco_pionCosTheta < 0.39"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0.39 && cc1pi_truth_pionCosTheta < 0.65  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionCosTheta >= 0.39  && cc1pi_reco_pionCosTheta < 0.65"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0.65 && cc1pi_truth_pionCosTheta < 0.84  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionCosTheta >= 0.65  && cc1pi_reco_pionCosTheta < 0.84"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0.84 && cc1pi_truth_pionCosTheta < 0.93  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionCosTheta >= 0.84  && cc1pi_reco_pionCosTheta < 0.93"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0.93 && cc1pi_truth_pionCosTheta <= 1    &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionCosTheta >= 0.93  && cc1pi_reco_pionCosTheta <= 1"},

        // Muon pion angle
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 0 && cc1pi_truth_muonPionAngle < 0.49      &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonPionAngle >= 0    && cc1pi_reco_muonPionAngle < 0.49"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 0.49 && cc1pi_truth_muonPionAngle < 0.93   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonPionAngle >= 0.49 && cc1pi_reco_muonPionAngle < 0.93"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 0.93 && cc1pi_truth_muonPionAngle < 1.26   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonPionAngle >= 0.93 && cc1pi_reco_muonPionAngle < 1.26"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 1.26 && cc1pi_truth_muonPionAngle < 1.57   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonPionAngle >= 1.26 && cc1pi_reco_muonPionAngle < 1.57"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 1.57 && cc1pi_truth_muonPionAngle < 1.88   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonPionAngle >= 1.57 && cc1pi_reco_muonPionAngle < 1.88"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 1.88 && cc1pi_truth_muonPionAngle < 2.21   &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonPionAngle >= 1.88 && cc1pi_reco_muonPionAngle < 2.21"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 2.21 && cc1pi_truth_muonPionAngle <= 2.65  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_muonPionAngle >= 2.21 && cc1pi_reco_muonPionAngle <= 2.65"},
        
        // Contained Muon Selection
        // Muon momentum
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.15 && cc1pi_truth_muonMomentum < 0.23  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum  >= 0.15 && cc1pi_reco_muonMomentum < 0.23"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.23 && cc1pi_truth_muonMomentum < 0.32  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum  >= 0.23 && cc1pi_reco_muonMomentum < 0.32"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.32 && cc1pi_truth_muonMomentum < 0.45  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum  >= 0.32 && cc1pi_reco_muonMomentum < 0.45"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.45 && cc1pi_truth_muonMomentum < 0.66  &&  cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum  >= 0.45 && cc1pi_reco_muonMomentum < 0.66"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.66 &&                                      cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained && cc1pi_reco_muonMomentum  >= 0.66"},

        // Golden pion selection
        // Pion momentum
        {false, "cc1pi_signal && true_golden_cc1pi && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionMomentum < 0.16   &&  cc1pi_selected_golden && passed_likelyGoldenPion && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionMomentum  >= 0.1  &&  cc1pi_reco_pionMomentum  < 0.16"},
        {false, "cc1pi_signal && true_golden_cc1pi && cc1pi_truth_pionMomentum >= 0.16 && cc1pi_truth_pionMomentum < 0.19  &&  cc1pi_selected_golden && passed_likelyGoldenPion && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionMomentum  >= 0.16 &&  cc1pi_reco_pionMomentum  < 0.19"},
        {false, "cc1pi_signal && true_golden_cc1pi && cc1pi_truth_pionMomentum >= 0.19 && cc1pi_truth_pionMomentum < 0.22  &&  cc1pi_selected_golden && passed_likelyGoldenPion && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionMomentum  >= 0.19 &&  cc1pi_reco_pionMomentum  < 0.22"},
        {false, "cc1pi_signal && true_golden_cc1pi && cc1pi_truth_pionMomentum >= 0.22 &&                                      cc1pi_selected_golden && passed_likelyGoldenPion && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_reco_pionMomentum  >= 0.22"},
    };

    // Hard-coded truth cuts corresponding to each reco cut in the cuts vector
    std::vector<std::pair<bool, std::string>> truthCuts {
        // "all" truth cut
        {true, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1"},
        
        // Generic selection truth cut
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1"},
        
        // Muon cos theta bins truth cuts
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= -1 && cc1pi_truth_muonCosTheta < -0.27"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= -0.27 && cc1pi_truth_muonCosTheta < 0.29"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.29 && cc1pi_truth_muonCosTheta < 0.46"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.46 && cc1pi_truth_muonCosTheta < 0.58"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.58 && cc1pi_truth_muonCosTheta < 0.67"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.67 && cc1pi_truth_muonCosTheta < 0.77"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.77 && cc1pi_truth_muonCosTheta < 0.82"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.82 && cc1pi_truth_muonCosTheta < 0.88"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.88 && cc1pi_truth_muonCosTheta < 0.93"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.93 && cc1pi_truth_muonCosTheta < 0.97"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonCosTheta >= 0.97 && cc1pi_truth_muonCosTheta <= 1"},
        
        // Pion cos theta bins truth cuts
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= -1 && cc1pi_truth_pionCosTheta < -0.47"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= -0.47 && cc1pi_truth_pionCosTheta < 0"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0 && cc1pi_truth_pionCosTheta < 0.39"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0.39 && cc1pi_truth_pionCosTheta < 0.65"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0.65 && cc1pi_truth_pionCosTheta < 0.84"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0.84 && cc1pi_truth_pionCosTheta < 0.93"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionCosTheta >= 0.93 && cc1pi_truth_pionCosTheta <= 1"},
        
        // Muon pion angle bins truth cuts
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 0 && cc1pi_truth_muonPionAngle < 0.49"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 0.49 && cc1pi_truth_muonPionAngle < 0.93"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 0.93 && cc1pi_truth_muonPionAngle < 1.26"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 1.26 && cc1pi_truth_muonPionAngle < 1.57"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 1.57 && cc1pi_truth_muonPionAngle < 1.88"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 1.88 && cc1pi_truth_muonPionAngle < 2.21"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonPionAngle >= 2.21 && cc1pi_truth_muonPionAngle <= 2.65"},
        
        // Contained muon momentum bins truth cuts
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.15 && cc1pi_truth_muonMomentum < 0.23"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.23 && cc1pi_truth_muonMomentum < 0.32"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.32 && cc1pi_truth_muonMomentum < 0.45"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.45 && cc1pi_truth_muonMomentum < 0.66"},
        {false, "cc1pi_signal && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_muonMomentum >= 0.66"},
        
        // Golden pion momentum bins truth cuts
        {false, "cc1pi_signal && true_golden_cc1pi && cc1pi_truth_pionMomentum >= 0.1 && cc1pi_truth_pionMomentum < 0.16"},
        {false, "cc1pi_signal && true_golden_cc1pi && cc1pi_truth_pionMomentum >= 0.16 && cc1pi_truth_pionMomentum < 0.19"},
        {false, "cc1pi_signal && true_golden_cc1pi && cc1pi_truth_pionMomentum >= 0.19 && cc1pi_truth_pionMomentum < 0.22"},
        {false, "cc1pi_signal && true_golden_cc1pi && cc1pi_truth_pionMomentum >= 0.22"}
    };

    // Vector of pairs; pair: show in latex, cut name
    std::vector<std::pair<bool, std::string>> cutsLatex {
        {true, "all"},
        {true, "genericSelection"},
        {true, "muonCosTheta bin 1: -1 to -0.27"},
        {true, "muonCosTheta bin 2: -0.27 to 0.29"},
        {true, "muonCosTheta bin 3: 0.29 to 0.46"},
        {true, "muonCosTheta bin 4: 0.46 to 0.58"},
        {true, "muonCosTheta bin 5: 0.58 to 0.67"},
        {true, "muonCosTheta bin 6: 0.67 to 0.77"},
        {true, "muonCosTheta bin 7: 0.77 to 0.82"},
        {true, "muonCosTheta bin 8: 0.82 to 0.88"},
        {true, "muonCosTheta bin 9: 0.88 to 0.93"},
        {true, "muonCosTheta bin 10: 0.93 to 0.97"},
        {true, "muonCosTheta bin 11: 0.97 to 1"},
        
        // Pion cos theta bins
        {true, "pionCosTheta bin 1: -1 to -0.47"},
        {true, "pionCosTheta bin 2: -0.47 to 0"},
        {true, "pionCosTheta bin 3: 0 to 0.39"},
        {true, "pionCosTheta bin 4: 0.39 to 0.65"},
        {true, "pionCosTheta bin 5: 0.65 to 0.84"},
        {true, "pionCosTheta bin 6: 0.84 to 0.93"},
        {true, "pionCosTheta bin 7: 0.93 to 1"},
        
        // Muon pion angle bins
        {true, "muonPionAngle bin 1: 0 to 0.49"},
        {true, "muonPionAngle bin 2: 0.49 to 0.93"},
        {true, "muonPionAngle bin 3: 0.93 to 1.26"},
        {true, "muonPionAngle bin 4: 1.26 to 1.57"},
        {true, "muonPionAngle bin 5: 1.57 to 1.88"},
        {true, "muonPionAngle bin 6: 1.88 to 2.21"},
        {true, "muonPionAngle bin 7: 2.21 to 2.65"},
        
        // Contained Muon Selection
        {true, "containedMuonMomentum bin 1: 0.15 to 0.23"},
        {true, "containedMuonMomentum bin 2: 0.23 to 0.32"},
        {true, "containedMuonMomentum bin 3: 0.32 to 0.45"},
        {true, "containedMuonMomentum bin 4: 0.45 to 0.66"},
        {true, "containedMuonMomentum bin 5: 0.66+"},
        
        // Golden pion selection
        {true, "goldenPionMomentum bin 1: 0.1 to 0.16"},
        {true, "goldenPionMomentum bin 2: 0.16 to 0.19"},
        {true, "goldenPionMomentum bin 3: 0.19 to 0.22"},
        {true, "goldenPionMomentum bin 4: 0.22+"}
    };

    // Map for all events; map: event type, cut, run
    std::map<std::string, std::map<std::string, std::map<std::string, double>>> eventCountMap;

    // Map for signal events; map: cut, run
    std::map<std::string, std::map<std::string, double>> signalCountMap, goldenCountMap;
    
    // Map for truth bin events; map: cut index, run
    std::map<int, std::map<std::string, double>> truthBinCountMap;
    
    std::vector<std::string> runs = {"0", "1", "2", "3", "4bcd", "5"}; // 0 is the total of all runs

    for (int i = 0; i < cuts.size(); i++) {
        const auto& [stopOnFailedCut, cut] = cuts[i];
        for (const auto& run : runs) {
            eventCountMap["data"][cut][run] = 0.0f;
            eventCountMap["prediction"][cut][run] = 0.0f;
            signalCountMap[cut][run] = 0.0f;
            goldenCountMap[cut][run] = 0.0f;
            truthBinCountMap[i][run] = 0.0f;
        }
    }

    double totalMCSum0 = 0.f;
    // Loop over the files
    for (const auto &[fileType, run, path, fileWeight] : files)
    {
        std::cout << "Processing file " << path <<std::endl;
        // Open the file and get the tree
        TFile* tFile = TFile::Open(path.c_str());
        TTree* tree = (TTree*)tFile->Get("stv_tree");

        const std::string eventType = (fileType == "beam on") ? "data" : "prediction";
        const bool isMC = (fileType == "nu mc" || fileType == "dirt mc");

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

        // Create TTreeFormula objects for truth cuts
        std::vector<TTreeFormula*> truthFormulas;
        for (const auto& [stopOnFailedCut, truthCut] : truthCuts) {
            truthFormulas.push_back(new TTreeFormula(("truth_" + truthCut).c_str(), truthCut.c_str(), tree));
        }

        const auto nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++)
        {
            if(i%(nEntries/100)==0) std::cout<<"\r"<<(100*i)/nEntries<<"%"<<std::flush;
            tree->GetEntry(i);

            // Calculate weight
            double weight = (isMC && std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30) ? spline_weight*tuned_cv_weight : 1;
            weight *= fileWeight;
            
            // Count events in truth bins for denominator in efficiency calculation
            for (int j = 0; j < truthFormulas.size(); j++) {
                if (truthFormulas[j]->EvalInstance()) {
                    truthBinCountMap[j][run] += weight;
                    truthBinCountMap[j]["0"] += weight;
                }
                else {
                    if(truthCuts[j].first) 
                        break;
                }
            }

            for (int j = 0; j < cuts.size(); j++)
            {
                const auto& [stopOnFailedCut, cut] = cuts[j];
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
                {
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
        
        // Clean up TTreeFormula objects
        for (auto formula : truthFormulas) {
            delete formula;
        }
        truthFormulas.clear();
        
        tFile->Close();
    } // End of loop over files

    std::cout<<std::endl;

    // *********************************
    // Save the table again as latex
    // *********************************
    for (const auto& run : runs) {
        // Open the output file
        std::ofstream outputFileLatex("eventCountTable_run" + run + "_onlyTesting_lowPiMomThreshold_fixed_muonContainedOption_GENIEFSIBugFixByBin_withTotalSignal.tex");

        // Write the header row
        std::string header = 
R"(
\begin{tabular}{lrrrrrrrrr}
\large{\textbf{Generic Selection}} \vspace{1mm}\\
\textbf{Cuts for Run%s} & \rot{\textbf{Signal}} & \rot{\textbf{Background}} & \rot{\textbf{Total Signal}} & \rot{\textbf{Efficiency} ($E$)} & \rot{\textbf{Purity} ($P$)} & \rot{$E \times P$} & \rot{\textbf{Golden Fraction}} & \rot{\textbf{Data}} & \rot{\textbf{Ratio}} \\
\midrule
\textbf{Selections} \\
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
            
            // Modified efficiency calculation using truth bin count as denominator
            double efficiency = 0.0;
            if (truthBinCountMap[i][run] > 0) {
                efficiency = signalCountMap.at(cut).at(run) / truthBinCountMap[i][run];
            }
            
            const auto purity = signalCountMap.at(cut).at(run) / eventCountMap.at("prediction").at(cut).at(run);
            const auto effpur = efficiency*purity;
            const auto goldenFraction = goldenCountMap.at(cut).at(run) / signalCountMap.at(cut).at(run);
            outputFileLatex << cutLatexLabel;

            if(cutLatexShow)
            {
                outputFileLatex << " & " << std::round(10*signalCountMap.at(cut).at(run))/10.0 // Round to 1 decimal place
                << " & " << std::round(10*(eventCountMap.at("prediction").at(cut).at(run) - signalCountMap.at(cut).at(run)))/10.0 // Round to 1 decimal place
                << " & " << std::round(10*truthBinCountMap[i][run])/10.0 // Add Total Signal column (rounded to 1 decimal place)
                << " & " << std::round(1000*efficiency)/10.0 << "\\%" // Convert to percentage and round to 1 decimal place
                << " & " << std::round(1000*purity)/10.0 << "\\%" // Convert to percentage and round to 1 decimal place
                << " & " << std::round(1000*effpur)/1000.0 // Round to 3 decimal places
                << " & " << std::round(1000*goldenFraction)/10.0 << "\\%" // Convert to percentage and round to 1 decimal place
                << " & " << eventCountMap.at("data").at(cut).at(run)
                << " & " << std::round(100*dataMCRatio)/100.0; // Round to 2 decimal places
                outputFileLatex << " \\\\\n";
            }
            else
            {
                outputFileLatex << " & & & & & & & & & \\\\\n";
            }

            if(cut == "cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1")
            {
                outputFileLatex << "\\midrule\n";
                outputFileLatex << "\\textbf{Contained Muon Selection}\\\\\n";
                outputFileLatex << "Cut applied to generic selection";
                outputFileLatex << " \\\\\\midrule\n";
            } else if(cut == "cc1pi_selected_generic && passed_openingAngle && cc1pi_reco_muonMomentum > 0.15 && cc1pi_reco_pionMomentum > 0.1 && cc1pi_recoMuon_IsContained")
            {
                outputFileLatex << "\\midrule\n";
                outputFileLatex << "\\textbf{Golden Pion Selection}\\\\\n";
                outputFileLatex << "Cut applied to generic selection";
                outputFileLatex << " \\\\\\midrule\n";
            }
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
        std::cout << "Truth bin count for run " << run << " at cut \"all\": " << truthBinCountMap[0][run] << std::endl;
        if (run != "0") {
            totalMCSum += eventCountMap.at("prediction").at("all").at(run);
            totalSignalSum += signalCountMap.at("all").at(run);
        }
    }
    std::cout << "Total MC sum runs 1-5: " << totalMCSum << std::endl;
    std::cout << "Total Signal sum runs 1-5: " << totalSignalSum << std::endl;

    std::cout<<"Done! Output written to eventCountTable_run*_GENIEFSIBugFixByBin_withTotalSignal.tex"<<std::endl;
}