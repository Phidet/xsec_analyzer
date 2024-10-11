void processFile(const std::string& filePath) {
    // Open the file
    TFile *file = TFile::Open(filePath.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return;
    }

    // Define the x-axis variables and categories
    // std::vector<std::string> xAxisVariables = {"TrueMuonCosThetaPlot", "TrueMuonPhiPlot", "TrueMuonMomentumPlot", "TruePionCosThetaPlot", "TruePionPhiPlot", "TruePionMomentumPlot", "TrueMuonPionOpeningAnglePlot", "TrueTotalPlot"};
    std::vector<std::string> xAxisVariables = {"TrueNuEPlot", "TrueNuEPlot_noPhaseSpaceRestriction", "TrueNuEPlot_muonMomentumRestricted", "TrueNuEPlot_muonAndPionMomentumRestricted"};

    std::vector<std::string> categories = {"QE", "MEC", "RES", "DIS", "COH"}; // "" is the total
    
    // Define colors for each category
    int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan};
    
    // Loop over each x-axis variable
    for (const auto& xAxisVar : xAxisVariables) {
        // Create a canvas for the current x-axis variable
        TCanvas *c = new TCanvas(xAxisVar.c_str(), xAxisVar.c_str(), 800, 600);
    
        // Create a stack for the current x-axis variable
        THStack *hs = new THStack((xAxisVar + "_stack").c_str(), (xAxisVar + ";X;Events").c_str());
    
        // Loop over each category
        for (size_t i = 0; i < categories.size(); ++i) {
            const auto& category = categories[i];
            std::string histName = category + xAxisVar; // Construct the histogram name
    
            // Retrieve the histogram
            TH1D *hist = dynamic_cast<TH1D*>(file->Get(histName.c_str()));
            if (!hist) {
                std::cerr << "Error: Could not find histogram " << histName << std::endl;
                continue;
            }
    
            // Prepare the histogram
            hist->SetFillColor(colors[i % (sizeof(colors) / sizeof(colors[0]))]); // Use colors array for fill color
            hist->SetFillStyle(1001); // Solid fill
    
            // Add the histogram to the stack
            hs->Add(hist);
        }
    
        // Set the y-axis maximum
        hs->SetMaximum(2.4);
    
        // Draw the stack
        hs->Draw("hist");
    
        // Add a legend
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        for (size_t i = 0; i < categories.size(); ++i) {
            const auto& category = categories[i];
            std::string histName = category + xAxisVar; // Construct the histogram name
    
            TH1D *hist = dynamic_cast<TH1D*>(file->Get(histName.c_str()));
            if (hist) {
                legend->AddEntry(hist, category.c_str(), "f");
            }
        }
        legend->Draw();
    
        // Save the canvas
        std::string outputFileName = filePath.substr(filePath.find_last_of("/\\") + 1);
        outputFileName = "plots/" + outputFileName.substr(0, outputFileName.find_last_of(".")) + "_" + xAxisVar + ".pdf";
        c->SaveAs(outputFileName.c_str());
    }
    
    // Clean up
    file->Close();
}

void signal_model_breakdown_plots() {
    std::vector<std::string> filePaths = {
        "/exp/uboone/app/users/jdetje/BuildEventGenerators/FlatTreeAnalyzer/OutputFiles/FlatTreeAnalyzerOutput_Genie_numu_numubar_CC_v3_4_2_AR23_20i_00_000.root",
        "/exp/uboone/app/users/jdetje/BuildEventGenerators/FlatTreeAnalyzer/OutputFiles/FlatTreeAnalyzerOutput_Genie_numu_numubar_CC_v3_4_2_G18_10a_02_11a.root"
    };

    for (const auto& filePath : filePaths) {
        processFile(filePath);
    }
}