void model_breakdown_plots() {
    const std::string new_path = "/exp/uboone/data/users/jdetje/generatorOutput/good_samples_july24/";

    std::vector<std::pair<std::string, std::string>> samples = {
        {new_path + "genie_2000000_14_1000180400_CC_v3_4_2_AR23_20i_00_000.flat.root", "genie_14_v3_4_2_AR23_20i_00_000"},
        {new_path + "genie_200000_-14_1000180400_CC_v3_4_2_AR23_20i_00_000.flat.root", "genie_-14_v3_4_2_AR23_20i_00_000"},
        {new_path + "genie_2000000_14_1000180400_CC_v3_4_2_G18_10a_02_11a.flat.root", "genie_14_v3_4_2_G18_10a_02_11a"},
        {new_path + "genie_200000_-14_1000180400_CC_v3_4_2_G18_10a_02_11a.flat.root", "genie_-14_v3_4_2_G18_10a_02_11a"}
    };

    for (const auto &[filePath, postfix] : samples)
    {
        // Open the file
        TFile *file = TFile::Open(filePath.c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open file " << filePath << std::endl;
            return; // Or handle the error as appropriate
        }
        
        // Get the tree
        TTree *tree = dynamic_cast<TTree*>(file->Get("FlatTree_VARS"));
        if (!tree) {
            std::cerr << "Error: Could not find TTree 'FlatTree_VARS' in file " << filePath << std::endl;
            file->Close(); // Clean up before exiting
            return; // Or handle the error as appropriate
        }
        
        // Define the grouped mode categories
        const int nCategories = 5;
        const char* categories[nCategories] = {"QE", "MEC", "RES", "DIS", "COH"};
        int modeValues[nCategories][7] = {
            {1},                // QE
            {2},                // MEC
            {10, 11, 12, 13, 17, 22, 23},  // RES
            {21, 26},           // DIS
            {16}                // COH
        };

        // Create an array to hold the histograms
        TH1F* hists[nCategories];

        // Create a canvas
        TCanvas *c1 = new TCanvas("c1", "Grouped Stacked Histogram", 800, 600);

        // Create a stack to hold the histograms
        THStack *hs = new THStack("hs", "Grouped Stacked Enu_true Histogram");

        // Define the histogram parameters
        int nBins = 50; // Number of bins in the histogram
        double xMin = 0; // Minimum x value
        double xMax = 3; // Maximum x value (adjust this according to your data)

        // Define an array of colors for the histograms
        int colors[nCategories] = {kRed, kBlue, kGreen, kMagenta, kCyan};

        // Loop over the categories and create histograms
        for (int i = 0; i < nCategories; ++i) {
            // Create a name and title for the histogram
            TString histName = Form("hist_%s", categories[i]);
            TString histTitle = Form("Enu_true for %s", categories[i]);

            // Create the histogram
            hists[i] = new TH1F(histName, histTitle, nBins, xMin, xMax);

            // Set the fill color for the histogram
            hists[i]->SetFillColor(colors[i]);
            hists[i]->SetLineStyle(0);

            // Construct the selection string
            TString selection = "(";
            for (int j = 0; j < sizeof(modeValues[i])/sizeof(modeValues[i][0]); ++j) {
                if (modeValues[i][j] == 0) break;  // Skip if the value is zero (not used)
                if (j != 0) selection += " || ";
                selection += Form("TMath::Abs(Mode) == %d", modeValues[i][j]);
            }
            selection += ")";

            // Draw the histogram with the selection
            tree->Draw(Form("Enu_true>>%s", histName.Data()), selection, "goff");

            // Add the histogram to the stack
            hs->Add(hists[i]);
        }


        // Draw the stack
        hs->Draw("hist");

        // Add axis labels and a legend
        hs->GetXaxis()->SetTitle("Enu_true");
        hs->GetYaxis()->SetTitle("Events");

        // Replace c1->BuildLegend(); with the following:
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust these coordinates as needed
        for (int i = 0; i < nCategories; ++i) {
            legend->AddEntry(hists[i], categories[i], "f");
        }
        legend->Draw();

        // Save the canvas as an image
        c1->SaveAs(("plots/model_breakdown_plot_" + postfix +".pdf").c_str());

        // Clean up
        file->Close();
    }
}