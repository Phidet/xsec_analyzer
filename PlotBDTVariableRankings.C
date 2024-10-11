#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <vector>
#include <string>
#include <memory>

void PlotBDTVariableRankings() 
{
    std::vector<std::pair<std::string, float>> goldenPionSeparation{ 
        {"ln(R_{p}/MIP)", 2.240e-01},
        {"Track score", 1.665e-01},
        {"Truncated mean dE/dx", 1.538e-01},
        {"Wiggliness", 6.130e-02},
        {"# Descendents", 3.464e-02},
        {"ln(R_{#pi}/MIP)", 2.393e-02}
    };

    std::vector<std::pair<std::string, float>> goldenPionImportance{
        {"ln(R_{p}/MIP)", 2.373e-01},
        {"Track score", 2.349e-01},
        {"Truncated mean dE/dx", 1.651e-01},
        {"ln(R_{#pi}/MIP)", 1.606e-01},
        {"Wiggliness", 1.115e-01},
        {"# Descendents", 9.069e-02}
    };

    std::vector<std::pair<std::string, float>> protonSeparation{
        {"ln(R_{p}/MIP)", 6.058e-01},
        {"Truncated mean dE/dx", 5.225e-01},
        {"Track score", 1.701e-01},
        {"ln(R_{#pi}/MIP)", 1.648e-01},
        {"Wiggliness", 6.342e-02}
    };

    std::vector<std::pair<std::string, float>> protonImportance{
        {"ln(R_{p}/MIP)", 2.699e-01},
        {"Truncated mean dE/dx", 2.359e-01},
        {"Track score", 2.157e-01},
        {"Wiggliness", 1.483e-01},
        {"ln(R_{#pi}/MIP)", 1.302e-01}
    };

    std::vector<std::pair<std::string, float>> muonSeparation{
        {"Track score", 3.739e-01},
        {"ln(R_{p}/MIP)", 3.194e-01},
        {"Truncated mean dE/dx", 3.000e-01},
        {"ln(R_{#pi}/MIP)", 2.122e-01},
        {"Wiggliness", 3.634e-02},
        {"# Descendents", 3.753e-03}
    };

    std::vector<std::pair<std::string, float>> muonImportance{
        {"ln(R_{p}/MIP)", 3.232e-01},
        {"Truncated mean dE/dx", 1.859e-01},
        {"Track score", 1.800e-01},
        {"ln(R_{#pi}/MIP)", 1.619e-01},
        {"Wiggliness", 1.026e-01},
        {"# Descendents", 4.649e-02}
    };

    std::map<std::string, std::map<std::string, std::vector<std::pair<std::string, float>>>> allRankings{
        {"goldenPion", {
            {"Separation", goldenPionSeparation},
            {"Importance", goldenPionImportance}
        }},
        {"proton", {
            {"Separation", protonSeparation},
            {"Importance", protonImportance}
        }},
        {"muon", {
            {"Separation", muonSeparation},
            {"Importance", muonImportance}
        }}
    };

    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------

    std::vector<std::string> bdtTypes = {"goldenPion", "muon", "proton"};
    std::vector<std::string> metrics = {"Separation", "Importance"};

    for (const auto& bdtType : bdtTypes) {
        for (const auto& metric : metrics) {
            // Get the data for this particle type and metric
            const auto& data = allRankings.at(bdtType).at(metric);

            std::cout<<"DEBUG data.size() = "<<data.size()<<std::endl;

            // Create a histogram
            TH1D h(("h "+bdtType+metric).c_str(), ("h "+bdtType+metric).c_str(), data.size(), 0, 1);

            // Set the title
            std::string bdtName = bdtType == "goldenPion" ? "Golden Pion" : (bdtType == "proton" ? "Proton" : "Muon");
            std::string title = bdtName + " BDT Input Variable Ranking";
            h.SetTitle(title.c_str());

            // Set the y axis label
            std::string yAxisLabel = metric == "Separation" ? "Separation score" : "Importance score";
            h.SetYTitle(yAxisLabel.c_str());

            // Fill the histogram in reverse order
            int i = 1;
            for (auto it = data.rbegin(); it != data.rend(); ++it) {
                const auto& [key, value] = *it;
                h.SetBinContent(i, value);
                h.GetXaxis()->SetBinLabel(i, key.c_str());
                ++i;
            }

            // // Center the axis labels on the bars
            // h.GetXaxis()->SetLabelOffset(-0.05);

            // Disable the stats box
            h.SetStats(0);

            // // Remove the y-axis ticks
            // h.GetYaxis()->SetNdivisions(0);

            // Access the y-axis and disable ticks
            h.GetXaxis()->SetTickLength(0); // Hides the ticks
            // h.GetXaxis()->SetLabelOffset(999); // Hides the labels

            // Create a canvas
            TCanvas c("c", "", 800, 600);

            // Disable the y-axis grid lines
            c.SetGridx(1);

            // Increase the left margin
            c.SetLeftMargin(0.2);

            // Draw the histogram as a horizontal bar chart
            h.SetBarWidth(0.5);
            h.SetBarOffset(0.1);

            if (metric == "Importance") {
                h.SetFillColor(kGreen);
            } else {
                h.SetFillColor(kBlue);
            }

            // Set the line color to black
            h.SetLineColor(kBlack);

            h.Draw("hbar");

            // Update the canvas
            c.Update();

            // Save the plot
            std::string filename = "plots/BDTVariableRankings_" + bdtType + "_" + metric;
            c.SaveAs((filename + ".pdf").c_str());
            c.SaveAs((filename + ".png").c_str());
            c.SaveAs((filename + ".C").c_str());                     
        }
    }
}