// Code by Patrick Green
// script to fake eventweight structure for NuMI dirt files for use with STV code
// reads real structure from input file with proper weights (e.g. nu overlay or nue overlay)
// sets all these real weights to 1 i.e. unchanged from CV then writes to output dirt file

#include <string>
#include <vector>
#include <iostream>
#include <map>

#include "TFile.h"
#include "TTree.h"

int fake_weights(const char* inputFile1, const char* inputFile2, const char* outputFile) {  	
  	std::cout << "Adding fake weights to file: " << inputFile2 << ", using structure from: " << inputFile1 << std::endl;

  	// load first file containing real weights
	TFile *f1 = new TFile(inputFile1);
	if(!f1->IsOpen()) {std::cout << "Could not open input file 1!" << std::endl; exit(1);}

	// load tree
	TTree *t1 = (TTree*)f1->Get("nuselection/NeutrinoSelectionFilter");

	// set addresses to weights map
	std::map<std::string, std::vector<double>> *weights = nullptr;
	t1->SetBranchAddress("weights", &weights);

	// get first entry to populate weights
	t1->GetEntry(0);
	
	// number of weights
	std::cout <<  "Number of weight pairs in map: " << weights->size() << std::endl;

	// create new set of fake weights to populate
	std::map<std::string, std::vector<double>> fake_weights;

	// loop through weights
	std::cout << std::endl;
	for (auto pair : *weights) {
		std::cout << "Name: " << pair.first << ", Size: " << (pair.second).size() << std::endl;
		// reset all weights to 1
		std::fill((pair.second).begin(), (pair.second).end(), 1);
		// add to new vector
		fake_weights.insert({pair.first, pair.second});
	}
	std::cout << std::endl;

	f1->Close();

	// open second file weights are to be added to
	TFile *f2 = new TFile(inputFile2);
	if(!f2->IsOpen()) {std::cout << "Could not open input file 2!" << std::endl; exit(1);}
	TTree *t2_nu = (TTree*)f2->Get("nuselection/NeutrinoSelectionFilter");
    TTree *t2_pot = (TTree*)f2->Get("nuselection/SubRun");

	// open third file to write out, cloning tree from second file
	TFile *f3 = new TFile(outputFile,"recreate");
	f3->cd();
    f3->mkdir("nuselection");
    f3->cd("/nuselection/");

    // create new tree, disabling partially filled weights branch
    t2_nu->SetBranchStatus("weights",0);
    TTree *t3_nu = t2_nu->CloneTree(0);
    TTree *t3_pot = t2_pot->CloneTree();

    // create new branch weights
    std::map<std::string, std::vector<double>> new_weights;
    t3_nu->Branch( "weights", "std::map<std::string, std::vector<double>>", &new_weights );

    // loop over events adding fake weights to new branch
    int n_entries = t2_nu->GetEntries();
	std::cout << "Number events: " << n_entries << std::endl;

	for (int e = 0; e < n_entries; e++) {

		// get current entry
	  	t2_nu->GetEntry(e);    	

	    if ( (e != 0) && (n_entries >= 10) &&  (e % (n_entries/10) == 0) ) {
	      std::cout << Form("%i0%% Completed...\n", e / (n_entries/10));
	    }

	    new_weights = fake_weights;
	    t3_nu->Fill();
	}

	// write output
    t3_nu->Write("NeutrinoSelectionFilter");
    t3_pot->Write("SubRun");

    f2->Close();
    f3->Close();

    delete f1;
    delete f2;
    delete f3;

	return 0;
}