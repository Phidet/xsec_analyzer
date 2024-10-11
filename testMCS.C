#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <iostream>

void printVector(const std::vector<float>& vec) {
    for (const auto& value : vec) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}

int testMCS() {
    // Open the first ROOT file
    TFile* file1 = new TFile("/uboone/data/users/jdetje/ubcc1piVSpelee/pelee/neutrinoselection_filt_0_4k.root", "READ");
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error: Could not open the first file." << std::endl;
        return 1;
    }

    // Open the second ROOT file
    TFile* file2 = new TFile("/uboone/data/users/jdetje/ubcc1piVSpelee/ubcc1pi/ubcc1piAnalysis_0_4k.root", "READ");
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error: Could not open the second file." << std::endl;
        return 1;
    }

    // Access the TTree in the first file
    TTree* tree1 = dynamic_cast<TTree*>(file1->Get("nuselection/NeutrinoSelectionFilter"));
    if (!tree1) {
        std::cerr << "Error: Could not find TTree 'nuselection' in the first file." << std::endl;
        return 1;
    }

    // Access the TTree in the second file
    TTree* tree2 = dynamic_cast<TTree*>(file2->Get("events"));
    if (!tree2) {
        std::cerr << "Error: Could not find TTree 'events' in the second file." << std::endl;
        return 1;
    }

    // Set branch addresses for the first file
    std::vector<float>* trk_mcs_muon_mom_v = nullptr;
    TBranch* branch1 = tree1->GetBranch("trk_mcs_muon_mom_v");
    if (branch1) {
        branch1->SetAddress(&trk_mcs_muon_mom_v);
    } else {
        std::cerr << "Error: Could not find TBranch 'trk_mcs_muon_mom_v' in the first file." << std::endl;
        return 1;
    }

    // Set branch addresses for the second file
    std::vector<float>* mcsMomentumForwardMuon = nullptr;
    std::vector<float>* mcsMomentumBackwardMuon = nullptr;
    TBranch* branch2_forward = tree2->GetBranch("reco_particle_mcsMomentumForwardMuon_vect");
    TBranch* branch2_backward = tree2->GetBranch("reco_particle_mcsMomentumBackwardMuon_vect");
    if (branch2_forward && branch2_backward) {
        branch2_forward->SetAddress(&mcsMomentumForwardMuon);
        branch2_backward->SetAddress(&mcsMomentumBackwardMuon);
    } else {
        std::cerr << "Error: Could not find TBranches in the second file." << std::endl;
        return 1;
    }

    // Check if both trees have the same number of entries
    if (tree1->GetEntries() != tree2->GetEntries()) {
        throw std::runtime_error("The trees have different numbers of entries.");
    }

    // Determine the number of entries to print
    Long64_t n = 20;
    Long64_t nEntries = std::min(n, tree1->GetEntries());

    // Loop over the first n entries in both files
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree1->GetEntry(i);
        tree2->GetEntry(i);

        std::cout << "Event " << i << ": ";

        if (trk_mcs_muon_mom_v) {
            std::cout << "Tree1 Vector = ";
            printVector(*trk_mcs_muon_mom_v);
        }

        if (mcsMomentumForwardMuon) {
            std::cout << "Tree2 Forward MCS Momentum = ";
            printVector(*mcsMomentumForwardMuon);
        }

        if (mcsMomentumBackwardMuon) {
            std::cout << "Tree2 Backward MCS Momentum = ";
            printVector(*mcsMomentumBackwardMuon);
        }

        std::cout << std::endl;
    }

    // Clean up
    file1->Close();
    file2->Close();

    return 0;
}
