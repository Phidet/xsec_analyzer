// Necessary ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"

void get_testing_events() {
    // Input and output files
    // std::string input_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/highestMuonBDTScoreAndRandomIndex/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_highestMuonBDT_randomIndex.root";
    // std::string input_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/highestMuonBDTScoreAndRandomIndex/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_highestMuonBDT_randomIndex.root";
    // std::string input_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/highestMuonBDTScoreAndRandomIndex/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_highestMuonBDT_randomIndex.root";
    // std::string input_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/highestMuonBDTScoreAndRandomIndex/overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_nu_ubcc1pi_highestMuonBDT_randomIndex.root";
    std::string input_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/highestMuonBDTScoreAndRandomIndex/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_highestMuonBDT_randomIndex.root";

    // std::string output_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/onlyTestEvents/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root";
    // std::string output_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/onlyTestEvents/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root";
    // std::string output_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/onlyTestEvents/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root";
    // std::string output_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/onlyTestEvents/overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_nu_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root";
    std::string output_file_path = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/onlyTestEvents/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root";


    // Check if the output file already exists
    if (gSystem->AccessPathName(output_file_path.c_str()) == 0) {
        std::cerr << "Output file already exists!" << std::endl;
        return;
    }

    // Open the root file
    TFile* file = new TFile(input_file_path.c_str(), "READ");

    // Create a new root file
    TFile* new_file = new TFile(output_file_path.c_str(), "NEW");

    // Access the TParameter<float> object
    TParameter<float>* summed_pot_param = (TParameter<float>*)file->Get("summed_pot");

    // Check if the parameter exists
    if (summed_pot_param) {
        // View the current value
        float current_value = summed_pot_param->GetVal();
        std::cout << "Current summed_pot value: " << current_value << std::endl;

        // Change the value (for example, set it to a new value)
        float new_value = current_value/2.0; // New value you want to set
        summed_pot_param->SetVal(new_value);

        // Make new_file the current directory
        new_file->cd();

        // Write the changes to the new file
        summed_pot_param->Write();

        std::cout << "summed_pot value has been changed to: " << new_value << std::endl;
    } else {
        std::cerr << "summed_pot parameter not found in the file!" << std::endl;
    }

    // Access the stv_tree TTree object
    TTree* tree = (TTree*)file->Get("stv_tree");

    // Check if the TTree exists
    if (tree) {
        // Create a new TTree that is a subset of the original TTree where isTrainingEvent is false
        TTree* new_tree = tree->CopyTree("isTrainingEvent == 0");

        // Print the number of copied events
        int num_copied_events = new_tree->GetEntries();
        std::cout << "Copied " << num_copied_events << " events out of " << tree->GetEntries() << " events." << std::endl;

        // Make new_file the current directory
        new_file->cd();

        // Write the new TTree to the new file
        new_tree->Write();
    } else {
        std::cerr << "stv_tree TTree not found in the file!" << std::endl;
    }

    // Close the files
    file->Close();
    new_file->Close();
    delete file;
    delete new_file;
}
