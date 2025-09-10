#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

// File list with scaling factor (use 1.0 for all files)
struct FileEntry {
    int id;
    std::string path;
    float scale;
};

void ccInclusiveChi2PIDPlots() {
    // Disable the stats box for all histograms
    gStyle->SetOptStat(0);
    
    // Define the list of ROOT files
    std::vector<FileEntry> files = {
        {1, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu.root", 1.0},
        {2, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu.root", 1.0},
        {3, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu.root", 1.0},
        {5, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_73_run4b_nu.root", 1.0},
        {6, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4c_nu.root", 1.0},
        {7, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4d_nu.root", 1.0},
        {8, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5.root", 1.0}
    };

    // Define 1D histograms for chi2 distributions per particle category
    // Categories: "p", "#mu^{#pm}", "#pi^{#pm}", "Other"
    const int bins = 50;
    const float xminP = 0.0;
    const float xmaxP = 350.0;

    const float xminMu = 0.0;
    const float xmaxMu = 70.0;

    // Histograms for chi2ForwardProtonW
    std::map<std::string, TH1F*> hChi2Proton;
    hChi2Proton["p"] = new TH1F("hChi2Proton_proton", "chi2ForwardProtonW for Proton;#chi_{p}^{2};Particle Count", bins, xminP, xmaxP);
    hChi2Proton["#mu^{#pm}"]   = new TH1F("hChi2Proton_muon",   "chi2ForwardProtonW for Muon;#chi_{p}^{2};Particle Count", bins, xminP, xmaxP);
    hChi2Proton["#pi^{#pm}"]   = new TH1F("hChi2Proton_pion+",   "chi2ForwardProtonW for Charged Pion;#chi_{p}^{2};Particle Count", bins, xminP, xmaxP);
    // hChi2Proton["#pi^0"]     = new TH1F("hChi2Proton_pion0",  "chi2ForwardProtonW for Pion0;#chi_{p}^{2};Particle Count", bins, xminP, xmaxP);
    // hChi2Proton["e^#pm"]     = new TH1F("hChi2Proton_electron",  "chi2ForwardProtonW for Electron;#chi_{p}^{2};Particle Count", bins, xminP, xmaxP);
    // hChi2Proton["#gamma"]    = new TH1F("hChi2Proton_photon",  "chi2ForwardProtonW for Photon;#chi_{p}^{2};Particle Count", bins, xminP, xmaxP);
    hChi2Proton["Other"]     = new TH1F("hChi2Proton_other",  "chi2ForwardProtonW for Other;#chi_{p}^{2};Particle Count", bins, xminP, xmaxP);

    // Histograms for chi2ForwardMuonW
    std::map<std::string, TH1F*> hChi2Muon;
    hChi2Muon["p"] = new TH1F("hChi2Muon_proton", "chi2ForwardMuonW for Proton;#chi_{#mu}^{2};Particle Count", bins, xminMu, xmaxMu);
    hChi2Muon["#mu^{#pm}"]   = new TH1F("hChi2Muon_muon",   "chi2ForwardMuonW for Muon;#chi_{#mu}^{2};Particle Count", bins, xminMu, xmaxMu);
    hChi2Muon["#pi^{#pm}"]   = new TH1F("hChi2Muon_pion+",   "chi2ForwardMuonW for Charged Pion;#chi_{#mu}^{2};Particle Count", bins, xminMu, xmaxMu);
    // hChi2Muon["#pi^0"]     = new TH1F("hChi2Muon_pion0",  "chi2ForwardMuonW for Pion0;#chi_{#mu}^{2};Particle Count", bins, xminMu, xmaxMu);
    // hChi2Muon["e^#pm"]     = new TH1F("hChi2Muon_electron",  "chi2ForwardMuonW for Electron;#chi_{#mu}^{2};Particle Count", bins, xminMu, xmaxMu);
    // hChi2Muon["#gamma"]    = new TH1F("hChi2Muon_photon",  "chi2ForwardMuonW for Photon;#chi_{#mu}^{2};Particle Count", bins, xminMu, xmaxMu);
    hChi2Muon["Other"]     = new TH1F("hChi2Muon_other",  "chi2ForwardMuonW for Other;#chi_{#mu}^{2};Particle Count", bins, xminMu, xmaxMu);

    // Loop over files and fill histograms
    for (const auto& file : files) {
        std::cout << "Processing file: " << file.path << std::endl;
        TFile *f = TFile::Open(file.path.c_str());
        if (!f || !f->IsOpen()) {
            std::cerr << "Error opening file: " << file.path << std::endl;
            continue;
        }
        // Assumed TTree name
        TTree* t = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
        if (!t) {
            std::cerr << "TTree not found in file: " << file.path << std::endl;
            f->Close();
            continue;
        }
        
        // Branch variables
        std::vector<int>* pdgBacktracked = nullptr;
        std::vector<float>* chi2ForwardProtonW = nullptr;
        std::vector<float>* chi2ForwardMuonW = nullptr;
        std::vector<float>* pfp_generation_v = nullptr;
        // New branch variables for event weights and pfp_generation_v
        const bool isMC = true;
        // Change type from double to float to match branch types
        float spline_weight = 0;
        float tuned_cv_weight = 0;

        t->SetBranchAddress("backtracked_pdg", &pdgBacktracked);
        t->SetBranchAddress("trk_pid_chipr_v", &chi2ForwardProtonW);
        t->SetBranchAddress("trk_pid_chimu_v", &chi2ForwardMuonW);
        t->SetBranchAddress("pfp_generation_v", &pfp_generation_v);
        // Set branch addresses for weight variables
        t->SetBranchAddress("weightSpline", &spline_weight);
        t->SetBranchAddress("weightTune", &tuned_cv_weight);

        const int nEntries = t->GetEntries();
        std::cout << "Number of entries: " << nEntries << std::endl;
        for (int i = 0; i < nEntries; i++) {
            std::cout << "Processing entry: " << i << std::endl;
            t->GetEntry(i);
            std::cout << "Obtained entry: " << i << std::endl;

            // Calculate event weight (convert product to double for precision if needed)
            double weight = ( isMC 
                              && std::isfinite(spline_weight*tuned_cv_weight) 
                              && spline_weight*tuned_cv_weight >= 0 
                              && spline_weight*tuned_cv_weight <= 30 ) ?
                              spline_weight*tuned_cv_weight : 1;
            if (!pdgBacktracked || !chi2ForwardProtonW || !chi2ForwardMuonW) continue;
            size_t nParticles = pdgBacktracked->size();
            for (size_t j = 0; j < nParticles; j++) {
                const auto generation = (*pfp_generation_v)[j];
                if (generation != 2) continue;

                int pdg = (*pdgBacktracked)[j];
                float valProton = (*chi2ForwardProtonW)[j];
                float valMuon   = (*chi2ForwardMuonW)[j];
                std::string category;
                int absPdg = std::abs(pdg);
                if (absPdg == 2212) {
                    category = "p";
                } else if (absPdg == 13) {
                    category = "#mu^{#pm}";
                } else if (absPdg == 211) {
                    category = "#pi^{#pm}";
                // } else if (absPdg == 111) {
                //     category = "#pi^0";
                // } else if (absPdg == 11) {
                //     category = "e^#pm";
                // } else if (absPdg == 22) {
                //     category = "#gamma";
                } else {
                    category = "Other";
                }
                hChi2Proton[category]->Fill(valProton, weight);
                hChi2Muon[category]->Fill(valMuon, weight);
            }
        }
        f->Close();
        delete f;
    }

    // Set line colors for visual clarity
    hChi2Proton["p"]->SetLineColor(kRed);
    hChi2Proton["#mu^{#pm}"]->SetLineColor(kBlue);
    hChi2Proton["#pi^{#pm}"]->SetLineColor(kGreen+2);
    // hChi2Proton["#pi^0"]->SetLineColor(kOrange);
    // hChi2Proton["e^#pm"]->SetLineColor(kCyan);
    // hChi2Muon["gamma"]->SetLineColor(kYellow);
    hChi2Proton["Other"]->SetLineColor(kMagenta);

    hChi2Muon["p"]->SetLineColor(kRed);
    hChi2Muon["#mu^{#pm}"]->SetLineColor(kBlue);
    hChi2Muon["#pi^{#pm}"]->SetLineColor(kGreen+2);
    // hChi2Muon["#pi^0"]->SetLineColor(kOrange);
    // hChi2Muon["e^#pm"]->SetLineColor(kCyan);
    // hChi2Muon["gamma"]->SetLineColor(kYellow);
    hChi2Muon["Other"]->SetLineColor(kMagenta);

    // Create canvas and legend for chi2ForwardProtonW plot
    TCanvas* c1 = new TCanvas("c1", "chi2ForwardProtonW Distributions", 800, 600);
    bool first = true;
    TLegend* leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    double ymax = 0;
    for (auto &kp : hChi2Proton) {
        kp.second->SetLineWidth(2);
        // kp.second->Scale(1.0 / kp.second->Integral());
        ymax = std::max(ymax, kp.second->GetMaximum());
    }
    for (auto &kp : hChi2Proton) {
        kp.second->GetYaxis()->SetRangeUser(0, ymax * 1.1);
        if (first) 
        { 
            kp.second->SetTitle("");
            kp.second->Draw("HIST"); first = false; 
        }
        else
        {
            kp.second->Draw("HIST SAME"); 
        }
        leg1->AddEntry(kp.second, kp.first.c_str(), "l");
    } 

    leg1->Draw();
    c1->Update();  // Force the canvas to update so that legend appears

    c1->SaveAs("plots/chi2ForwardProtonW_Distributions.pdf");
    c1->SaveAs("plots/chi2ForwardProtonW_Distributions.png");

    // Create canvas and legend for chi2ForwardMuonW plot
    TCanvas* c2 = new TCanvas("c2", "chi2ForwardMuonW Distributions", 800, 600);
    first = true;
    TLegend* leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    ymax = 0;
    for (auto &km : hChi2Muon) {
        km.second->SetLineWidth(2);
        // km.second->Scale(1.0 / km.second->Integral());
        ymax = std::max(ymax, km.second->GetMaximum());
    }
    for (auto &km : hChi2Muon) {
        km.second->GetYaxis()->SetRangeUser(0, ymax * 1.1);
        // Force scientific notation with exponent (e.g. "40" and x10^3)
        // km.second->GetYaxis()->SetMaxDigits(3);
        if (first) {
            km.second->SetTitle("");
            km.second->Draw("HIST");
            first = false; 
        }
        else {
            km.second->Draw("HIST SAME");
        }
        leg2->AddEntry(km.second, km.first.c_str(), "l");
    }
    leg2->Draw();
    c2->Update();  // Force the canvas to update so that legend appears

    c2->SaveAs("plots/chi2ForwardMuonW_Distributions.pdf");
    c2->SaveAs("plots/chi2ForwardMuonW_Distributions.png");

    // Cleanup
    delete c1;
    delete c2;
}

#ifndef __CINT__
int main() {
    ccInclusiveChi2PIDPlots();
    return 0;
}
#endif