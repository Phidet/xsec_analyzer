#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLeaf.h>
#include <TParameter.h> // Added for TParameter
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <tuple>

/**
 * Adds a new branch called reco_event_ccinc_highestMuonBDTScore to ROOT files.
 * This branch contains the highest value of reco_particle_ccinc_muonBDTScore 
 * (between -1 and 1) among contained generation 2 particles.
 * 
 * @param inputFile Path to the input ROOT file
 * @param outputFile Path to the output ROOT file. If empty, will append "_updated" to the input filename.
 */
void ProcessFile(const std::string& inputFile, const std::string& outputFile = "") {
    // Determine the output filename
    std::string outFileName = outputFile;
    if (outFileName.empty()) {
        size_t dotPos = inputFile.find_last_of('.');
        if (dotPos != std::string::npos) {
            outFileName = inputFile.substr(0, dotPos) + "_updated" + inputFile.substr(dotPos);
        } else {
            outFileName = inputFile + "_updated";
        }
    }
    
    std::cout << "Processing file: " << inputFile << std::endl;
    std::cout << "Output will be written to: " << outFileName << std::endl;
    
    // Open the input file
    TFile* inFile = TFile::Open(inputFile.c_str(), "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open input file: " << inputFile << std::endl;
        return;
    }
    
    // Check for and retrieve "summed_pot" if it exists
    TObject* summedPotObj = nullptr;
    if (inFile->Get("summed_pot")) {
        summedPotObj = inFile->Get("summed_pot");
        std::cout << "Found 'summed_pot' in input file." << std::endl;
    } else {
        std::cout << "'summed_pot' not found in input file." << std::endl;
    }

    // Get the tree
    TTree* inTree = (TTree*)inFile->Get("stv_tree");
    if (!inTree) {
        std::cerr << "Error: Could not find stv_tree in the input file" << std::endl;
        inFile->Close();
        return;
    }
    
    // Create output file and clone the tree
    TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
    TTree* outTree = inTree->CloneTree(0);
    
    // Add the new branch
    Float_t highestMuonBDTScore;
    TBranch* newBranch = outTree->Branch("reco_event_ccinc_highestMuonBDTScore", &highestMuonBDTScore, "reco_event_ccinc_highestMuonBDTScore/F");
    
    // Set up branch addresses for the input tree
    std::vector<float>* pReco_particle_ccinc_muonBDTScore = nullptr;
    std::vector<int>* pReco_particle_ccinc_generation = nullptr;
    // std::vector<bool>* pReco_particle_ccinc_isContained = nullptr;
    
    // Only set the needed branches to active
    // inTree->SetBranchStatus("*", 0); // Disable all branches
    // inTree->SetBranchStatus("passed_topologicalScoreCC", 1); // Enable only the needed branches
    // inTree->SetBranchStatus("reco_particle_ccinc_muonBDTScore", 1);
    // inTree->SetBranchStatus("reco_particle_ccinc_generation", 1);
    // inTree->SetBranchStatus("reco_particle_ccinc_contained", 1);

    inTree->SetBranchAddress("reco_particle_ccinc_muonBDTScore", &pReco_particle_ccinc_muonBDTScore);
    inTree->SetBranchAddress("reco_particle_ccinc_generation", &pReco_particle_ccinc_generation);
    // inTree->SetBranchAddress("reco_particle_ccinc_contained", &pReco_particle_ccinc_isContained);
    
    // Get total number of entries
    Long64_t nEntries = inTree->GetEntries();
    
    // Loop over all entries
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        // std::cout << "DEBUG Point -2" << std::endl;
        // Display progress
        if (entry % 10000 == 0 || entry == nEntries - 1) {
            std::cout << "\rProcessing entry " << entry << " of " << nEntries 
                      << " (" << (100.0 * entry / nEntries) << "%)" << std::flush;
        }

        // std::cout << "DEBUG Point -1" << std::endl;
        inTree->GetEntry(entry);
        // std::cout << "DEBUG Point 0" << std::endl;
        
        // Find the highest muon BDT score for contained generation 2 particles
        highestMuonBDTScore = -std::numeric_limits<float>::max();
        
        if(inTree->GetLeaf("passed_topologicalScoreCC")->GetValue())
        {
            for (size_t i = 0; i < pReco_particle_ccinc_generation->size(); i++) {
                const auto generation = pReco_particle_ccinc_generation->at(i);
                // const auto isContained = pReco_particle_ccinc_isContained->at(i);
                if (generation == 2){// && isContained) {
                    float score = pReco_particle_ccinc_muonBDTScore->at(i);
                    if (score > highestMuonBDTScore && score >= -1 && score <= 1) {
                        highestMuonBDTScore = score;
                    }
                }
            }
        }

        // std::cout << "DEBUG Point 1" << std::endl;

        // If no valid muon BDT score was found, we keep
        
        outTree->Fill();
        // std::cout << "DEBUG Point 2" << std::endl;
    }
    
    std::cout << std::endl << "Writing output tree..." << std::endl;
    outTree->Write();
    
    // Write "summed_pot" to the output file if it was found
    if (summedPotObj) {
        outFile->cd(); // Ensure we are in the context of the output file
        summedPotObj->Write("summed_pot");
        std::cout << "Copied 'summed_pot' to output file." << std::endl;
    }
    
    outFile->Close();
    inFile->Close();
    
    std::cout << "Processing complete: " << inputFile << " → " << outFileName << std::endl;
}

/**
 * Process a batch of files from the standard locations used in the analysis.
 */
void ProcessAllFiles() {
    const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/17July2025_fixedTrackDistance/complete/";
    const std::string rootPathMC = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/completeFixedCategoryForMC/";
    const std::string outputDir = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/fixedCategoryAndHighestMuonBDTScore/";   

    // List of files to process; tuple: file type, run, path
    const std::vector<std::tuple<std::string, std::string, std::string>> files = {
        // // Beam on files
        std::make_tuple("beam on", "1",    rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi.root"),
        std::make_tuple("beam on", "2",    rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root"),
        std::make_tuple("beam on", "3",    rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root"),
        std::make_tuple("beam on", "4bcd", rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_73_run4bcd_ubcc1pi.root"),
        std::make_tuple("beam on", "5",    rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root"),

        // // Beam off files
        std::make_tuple("beam off", "1",    rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi.root"),
        std::make_tuple("beam off", "2",    rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root"),
        std::make_tuple("beam off", "3",    rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root"),
        std::make_tuple("beam off", "4bcd", rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_73_run4bcd_ubcc1pi.root"),
        std::make_tuple("beam off", "5",    rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root"),

        // // Nu MC files
        std::make_tuple("nu mc", "1",    rootPathMC + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root"),
        std::make_tuple("nu mc", "2",    rootPathMC + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root"),
        std::make_tuple("nu mc", "3",    rootPathMC + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root"),
        std::make_tuple("nu mc", "4bcd", rootPathMC + "overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_nu_ubcc1pi.root"),
        std::make_tuple("nu mc", "5",    rootPathMC + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root"),

        // // Dirt MC files
        std::make_tuple("dirt mc", "1",    rootPathMC + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root"),
        std::make_tuple("dirt mc", "2",    rootPathMC + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root"),
        std::make_tuple("dirt mc", "3",    rootPathMC + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root"),
        std::make_tuple("dirt mc", "4bcd", rootPathMC + "overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_dirt_ubcc1pi.root"),
        std::make_tuple("dirt mc", "5",    rootPathMC + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root"),

        // DetVars files
        std::make_tuple("DetVars", "CVextra",        rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_CVextra_ubcc1pi.root"),
        std::make_tuple("DetVars", "CV",             rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_CV_virtualParentsFix_ubcc1pi.root"),
        std::make_tuple("DetVars", "LYAtteNuation",  rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_LYAttenuation_ubcc1pi.root"),
        std::make_tuple("DetVars", "LYDown",         rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_LYDown_ubcc1pi.root"),
        std::make_tuple("DetVars", "LYRayleigh",     rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_LYRayleigh_virtualParentsFix_ubcc1pi.root"),
        std::make_tuple("DetVars", "Recomb2",        rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_Recomb2_ubcc1pi.root"),
        std::make_tuple("DetVars", "SCE",            rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_SCE_ubcc1pi.root"),
        std::make_tuple("DetVars", "WireModThetaXZ", rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModThetaXZ_ubcc1pi.root"),
        std::make_tuple("DetVars", "WireModThetaYZ", rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModThetaYZ_ubcc1pi.root"),
        std::make_tuple("DetVars", "WireModX",       rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModX_virtualParentsFix_ubcc1pi.root"),
        std::make_tuple("DetVars", "WireModYZ",      rootPathMC + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModYZ_ubcc1pi.root"),

        // NuWro files
        std::make_tuple("NuWro", "1",    rootPathMC + "nuwro_fakedata_peleeTuple_uboone_run1_ubcc1pi.root"),
        std::make_tuple("NuWro", "2",    rootPathMC + "nuwro_fakedata_peleeTuple_uboone_run2_ubcc1pi.root"),
        std::make_tuple("NuWro", "3",    rootPathMC + "nuwro_fakedata_peleeTuple_uboone_run3_ubcc1pi.root"),
        std::make_tuple("NuWro", "4bcd", rootPathMC + "high_stat_prodgenie_bnb_nu_nuwro_overlay_run4_pelee_ubcc1pi.root"),
        std::make_tuple("NuWro", "5",    rootPathMC + "high_stat_prodgenie_bnb_nu_nuwro_overlay_run5_pelee_ubcc1pi.root")

    };
    
    // Process each file
    for (const auto& [type, run, path] : files) {
        // Extract the original filename from the full path
        std::string originalFileName = path.substr(path.find_last_of('/') + 1);
        
        // Construct the output path using the fixed directory and original filename
        std::string outputPath = outputDir + originalFileName;
        
        // Append "_highestMuonBDT" before the file extension
        size_t dotPos = outputPath.find_last_of('.');
        if (dotPos != std::string::npos) {
            outputPath = outputPath.substr(0, dotPos) + "_highestMuonBDT" + outputPath.substr(dotPos);
        } else {
            outputPath += "_highestMuonBDT";
        }
        
        ProcessFile(path, outputPath);
    }
}

/**
 * Main entry point that processes a single file or all files.
 * @param inputFile Path to a single input file to process, or empty to process all files
 * @param outputFile Optional path for the output file
 */
void AddHighestMuonBDTScoreBranch(const char* inputFile = "", const char* outputFile = "") {
    if (strlen(inputFile) > 0) {
        // Process a single file
        ProcessFile(inputFile, outputFile);
    } else {
        // Process all files from the standard locations
        ProcessAllFiles();
    }
}
