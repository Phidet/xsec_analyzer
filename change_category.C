#include <string>
#include <vector>
#include <iostream>
#include <map>

#include "TFile.h"
#include "TTree.h"

void change_category() {
    const std::string iPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/17July2025_fixedTrackDistance/complete/";
    const std::string oPath = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/completeFixedCategoryForMC/";
    // List of input and output file path pairs
    std::vector<std::pair<std::string, std::string>> files = {
    // {iPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi.root",
    //     oPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi.root"},
    // {iPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root",
    //     oPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root"},
    // {iPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root",
    //     oPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root"},
    // {iPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root",
    //     oPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root"},
    // {iPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_73_run4bcd_ubcc1pi.root",
    //     oPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_73_run4bcd_ubcc1pi.root"},
    // {iPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi.root",
    //     oPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi.root"},
    // {iPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root",
    //     oPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root"},
    // {iPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root",
    //     oPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root"},
    // {iPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_73_run4bcd_ubcc1pi.root",
    //     oPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_73_run4bcd_ubcc1pi.root"},
    // {iPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root",
    //     oPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_CV_virtualParentsFix_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_CV_virtualParentsFix_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_CVextra_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_CVextra_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_LYAttenuation_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_LYAttenuation_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_LYDown_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_LYDown_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_LYRayleigh_virtualParentsFix_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_LYRayleigh_virtualParentsFix_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_Recomb2_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_Recomb2_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_SCE_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_SCE_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModThetaXZ_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModThetaXZ_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModThetaYZ_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModThetaYZ_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModX_virtualParentsFix_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModX_virtualParentsFix_ubcc1pi.root"},
    {iPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModYZ_ubcc1pi.root",
        oPath + "detvar_run345_peleeTuple_uboone_v08_00_00_73_WireModYZ_ubcc1pi.root"},
    {iPath + "high_stat_prodgenie_bnb_nu_nuwro_overlay_run4_pelee_ubcc1pi.root",
        oPath + "high_stat_prodgenie_bnb_nu_nuwro_overlay_run4_pelee_ubcc1pi.root"},
    {iPath + "high_stat_prodgenie_bnb_nu_nuwro_overlay_run5_pelee_ubcc1pi.root",
        oPath + "high_stat_prodgenie_bnb_nu_nuwro_overlay_run5_pelee_ubcc1pi.root"},
    {iPath + "nuwro_fakedata_peleeTuple_uboone_run1_ubcc1pi.root",
        oPath + "nuwro_fakedata_peleeTuple_uboone_run1_ubcc1pi.root"},
    {iPath + "nuwro_fakedata_peleeTuple_uboone_run2_ubcc1pi.root",
        oPath + "nuwro_fakedata_peleeTuple_uboone_run2_ubcc1pi.root"},
    {iPath + "nuwro_fakedata_peleeTuple_uboone_run3_ubcc1pi.root",
        oPath + "nuwro_fakedata_peleeTuple_uboone_run3_ubcc1pi.root"},
    {iPath + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root",
        oPath + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root"},
    {iPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root",
        oPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root"},
    {iPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root",
        oPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root"},
    {iPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root",
        oPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root"},
    {iPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root",
        oPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root"},
    {iPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root",
        oPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root"},
    {iPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root",
        oPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root"},
    {iPath + "overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_dirt_ubcc1pi.root",
        oPath + "overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_dirt_ubcc1pi.root"},
    {iPath + "overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_nu_ubcc1pi.root",
        oPath + "overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_nu_ubcc1pi.root"},
    {iPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root",
        oPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root"},
    };

    for (const auto& filePair : files) {
        const char* inputFile1 = filePair.first.c_str();
        const char* outputFile = filePair.second.c_str();

        TFile *f1 = new TFile(inputFile1, "read");
        if(!f1->IsOpen()) {std::cout << "Could not open input file 1!" << std::endl; exit(1);}

        // load tree
        TTree *t1 = (TTree*)f1->Get("stv_tree");

        // load the summed_pot parameter
        TParameter<float> *summed_pot = (TParameter<float>*)f1->Get("summed_pot");

        // create an output file at outputFile that copies the entire input file
        TFile *f2 = new TFile(outputFile, "recreate");

        // write the summed_pot parameter to the output file
        summed_pot->Write();

        TTree *t2 = t1->CloneTree(0);  // clone the structure of the tree, but not the entries

        // then modify the content of the category variable in the output file stv_tree TTree
        float cc1pi_truth_pionMomentum;
        bool true_cc1pi;
        int category;

        t1->SetBranchAddress("cc1pi_truth_pionMomentum", &cc1pi_truth_pionMomentum);
        t1->SetBranchAddress("true_cc1pi", &true_cc1pi);
        t1->SetBranchAddress("category", &category);
        t2->Branch("category", &category, "category/I");

        Long64_t nentries = t1->GetEntries();
        int progress = 0;
        for (Long64_t i=0;i<nentries;i++) {
            t1->GetEntry(i);
            if (true_cc1pi && cc1pi_truth_pionMomentum < 0.1) {
                category = 11;
            }
            t2->Fill();

            // Update the progress bar
            int newProgress = (i * 100) / nentries;
            if (newProgress > progress) {
                progress = newProgress;
                std::cout << "\rProgress: " << progress << "%" << std::flush;
            }
        }

        // Finally save the output file with the modified category variable
        t2->Write();
        delete f2;
        delete f1;
    }
    return;
}