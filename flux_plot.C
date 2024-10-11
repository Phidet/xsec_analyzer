#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

void plotHistograms() {
    // Open the ROOT file
    TFile* file = new TFile("/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463_hist/MCC9_FluxHist_volTPCActive.root");

    // Get the histograms
    TH1* hEnumu_cv = (TH1*)file->Get("hEnumu_cv");
    TH1* hEnumubar_cv = (TH1*)file->Get("hEnumubar_cv");

    hEnumu_cv->Scale(1.11425e21 / (4997. * 5e8) / (256.35 * 233.));
    hEnumubar_cv->Scale(1.11425e21 / (4997. * 5e8) / (256.35 * 233.));

    // Set line colors and thickness
    hEnumu_cv->SetLineColor(kRed); // Set line color to red
    hEnumubar_cv->SetLineColor(kBlue); // Set line color to blue
    hEnumu_cv->SetLineWidth(2); // Set line thickness to 2
    hEnumubar_cv->SetLineWidth(2); // Set line thickness to 2

    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Histograms", 800, 600);

    hEnumu_cv->SetTitle("Predicted Active Volume Flux");
    hEnumu_cv->GetXaxis()->SetTitle("E_{#nu} / GeV");
    hEnumu_cv->GetYaxis()->SetTitle("#phi(E_{#nu}) / cm^{2}");
    hEnumu_cv->GetXaxis()->SetRangeUser(0, 5.5);

    // Draw the histograms on one plot
    hEnumu_cv->Draw("hist");
    hEnumubar_cv->Draw("hist same");

    // Remove the stats box
    gStyle->SetOptStat(0);

    // Add a legend
    TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    legend->SetHeader("Neutrino Type", "C");
    legend->SetTextSize(0.05);
    legend->AddEntry(hEnumu_cv, "#nu_{#mu}", "l");
    legend->AddEntry(hEnumubar_cv, "#bar{#nu_{#mu}}", "l");
    legend->Draw();

    // Save the plot as a PDF file
    canvas->SaveAs("plots/flux_histogram.pdf");

    // Clean up
    delete legend;
    delete canvas;
    file->Close();
    delete file;
}

void flux_plot() {
    plotHistograms();
}