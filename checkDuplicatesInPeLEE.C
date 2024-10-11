#include <iostream>
#include <TFile.h>
#include <TTree.h>

void printBranchValues(const char* filename) {
    // Open the ROOT file
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Failed to open file " << filename << std::endl;
        return;
    }

    // Access the TTree
    TTree* tree = dynamic_cast<TTree*>(file->Get("nuselection/NeutrinoSelectionFilter")); // Replace "tree_name" with the actual name of your TTree
    if (!tree) {
        std::cerr << "Error: Failed to retrieve TTree from file" << std::endl;
        file->Close();
        return;
    }

    // Variables to store branch values
    int run, sub, evt;

    // Set branch addresses
    tree->SetBranchAddress("run", &run);
    tree->SetBranchAddress("sub", &sub);
    tree->SetBranchAddress("evt", &evt);

    // Loop over events
    Long64_t numEvents = tree->GetEntries(); // Only check the first 10% of the events
    for (Long64_t i = 0; i < numEvents; ++i) {
        if (i % 1000 == 0) std::cerr << "Processing event " << i << " of " << numEvents << std::endl;
        tree->GetEntry(i);
        std::cout <<"run = " << run << ", sub = " << sub << ", evt = " << evt << std::endl;
    }

    // Close the file
    file->Close();
}

void checkDuplicatesInPeLEE() {
    const char* filename = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5.root"; // Replace "your_file.root" with the actual path to your .root file
    printBranchValues(filename);
}