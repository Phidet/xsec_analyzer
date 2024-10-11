void testParticles() {
    // Define the particle types you want to separate
    std::set<int> particleTypes = {211, 2212, 11, 13, 22};  // Add more if needed

    // Define the fill colors for each particle type
    std::map<int, int> fillColors = {{211, kRed}, {2212, kBlue}, {11, kGreen}, {13, kMagenta}, {22, kOrange}}; // Add more if needed

    // Create a map to hold the histograms for each particle type
    std::map<int, TH1D*> histograms;
    std::cout<<"DEBUG Point 0"<<std::endl;

    for (int pdg : particleTypes) {
        TString histName = Form("histogram_%d", pdg);
        TH1D *histogram = new TH1D(histName, Form("Particle Type %d", pdg), 10, 0, 1); // Adjust binning as needed
        histogram->SetLineColor(pdg); // Set a unique line color for each particle type
        histogram->SetFillColor(fillColors[pdg]); // Set a custom fill color for each particle type
        histograms[pdg] = histogram;
    }

    std::cout<<"DEBUG Point 1"<<std::endl;

    std::vector<std::string> inputFiles{
        "/uboone/data/users/jdetje/ubcc1piOutput/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root"
    };

    std::cout<<"DEBUG Point 2"<<std::endl;

    std::vector<int>* particleTypesVec;
    std::vector<float>* trackScores;
    std::cout<<"DEBUG Point 6"<<std::endl;

    for(const auto inputFile:inputFiles){

        std::cout<<"DEBUG Point 3"<<std::endl;
        TFile *file = new TFile(inputFile.c_str());
        std::cout<<"DEBUG Point 4"<<std::endl;
        TTree *tree = (TTree*)file->Get("stv_tree");
        std::cout<<"DEBUG Point 5"<<std::endl;
        // Loop over the tree entries and fill the histograms based on particle type
        int nEntries = tree->GetEntries();

        tree->SetBranchAddress("track_bestMatchedTruthPDG", &particleTypesVec);
        tree->SetBranchAddress("particle_cutValue_trackScore", &trackScores);
        std::cout<<"DEBUG Point 7"<<std::endl;

        for (int entry = 0; entry < nEntries; entry++) {
            std::cout<<"DEBUG Point 7.1"<<std::endl;
            tree->GetEntry(entry);
            std::cout<<"DEBUG Point 7.2"<<std::endl;

            for (unsigned int i = 0; i < particleTypesVec->size(); i++) {
                std::cout<<"DEBUG Point 7.3"<<std::endl;
                int pdg = std::abs(particleTypesVec->at(i));
                std::cout<<"DEBUG Point 7.4"<<std::endl;
                std::cout << "DEBUG pdg: " << pdg << std::endl;
                if (histograms.find(pdg) != histograms.end()) {
                    histograms[pdg]->Fill(trackScores->at(i));
            }
            }
        }
        std::cout<<"DEBUG Point 8"<<std::endl;

        delete file;
    } // end loop over files
    std::cout<<"DEBUG Point 9"<<std::endl;


    // Create a stacked histogram
    THStack *stackedHist = new THStack("stackedHist", "Stacked Histogram");
    for (int pdg : particleTypes) {
        if (histograms.find(pdg) != histograms.end()) {
        stackedHist->Add(histograms[pdg]);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Muon Momentum Histogram", 800, 600);

    // Draw the stacked histogram
    stackedHist->Draw(); // "nostack"
    stackedHist->SetTitle("Muon Momentum");
    stackedHist->GetXaxis()->SetTitle("Muon Momentum");
    stackedHist->GetYaxis()->SetTitle("Number of Events");

    // Add a legend to distinguish the particle types
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (int pdg : particleTypes) {
    if (histograms.find(pdg) != histograms.end()) {
        legend->AddEntry(histograms[pdg], Form("Type %d", pdg), "l");
    }
    }
    legend->Draw();

    c1->BuildLegend();
}