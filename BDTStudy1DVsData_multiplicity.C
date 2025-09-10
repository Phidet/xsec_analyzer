#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <memory>

void BDTStudy1DVsData_multiplicity()
{
    const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/1March24/";

    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        // std::make_tuple("beam on", "1", rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "2", rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "3", rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 1.f),
        // std::make_tuple("beam on", "4bcd", rootPath + "bnb_beam_on_peleeTuple_uboone_run4bcd_ubcc1pi.root", 1.f),
        std::make_tuple("beam on", "5", rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root", 1.f),

        // std::make_tuple("beam off", "1", rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi.root", 0.56421),
        // std::make_tuple("beam off", "2", rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 0.40202),
        // std::make_tuple("beam off", "3", rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 0.30657),
        // std::make_tuple("beam off", "4bcd", rootPath + "bnb_beam_off_peleeTuple_uboone_run4bcd_ubcc1pi.root", 0.30175),
        std::make_tuple("beam off", "5", rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root", 0.32807),

        // std::make_tuple("nu mc", "1", rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root", 0.13011),
        // std::make_tuple("nu mc", "2", rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root", 0.25750),
        // std::make_tuple("nu mc", "3", rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root", 0.20113),
        // std::make_tuple("nu mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi.root", 0.13074),
        std::make_tuple("nu mc", "5", rootPath + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi.root", 0.15196),

        // std::make_tuple("dirt mc", "1", rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root", 0.52806),
        // std::make_tuple("dirt mc", "2", rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root", 0.27521),
        // std::make_tuple("dirt mc", "3", rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root", 0.80892),
        // std::make_tuple("dirt mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi.root", 0.329301),
        std::make_tuple("dirt mc", "5", rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root", 0.41280)
    };

    const std::vector<std::string> types = {"nu mc", "beam off", "dirt mc", "beam on"};

    // Create histograms for each type
    std::map<std::string, std::unique_ptr<TH1D>> histograms;
    for (const auto& type : types)
    {
        histograms[type] = std::make_unique<TH1D>((type + "_hist").c_str(), ("Number of Particles (" + type + ");Number of Particles;Events").c_str(), 10, 0, 10);
        histograms[type]->Sumw2();
    }

    // Loop over the files
    for (const auto& [sampleType, run, filePath, runWeight] : files)
    {
        std::cout << "Processing file: " << filePath << std::endl;
        TFile* tFile = TFile::Open(filePath.c_str());
        if (!tFile || tFile->IsZombie())
        {
            std::cerr << "Error: Unable to open file: " << filePath << std::endl;
            continue;
        }

        TTree* tree = (TTree*)tFile->Get("stv_tree");
        if (!tree)
        {
            std::cerr << "Error: Unable to get tree from file: " << filePath << std::endl;
            tFile->Close();
            continue;
        }

        std::vector<int>* pReco_particle_ccinc_generation = nullptr;
        std::vector<bool>* pReco_particle_ccinc_isContained = nullptr;
        bool passed_topologicalScoreCC = false;
        bool passed_max1Uncontained = false;
        tree->SetBranchAddress("reco_particle_ccinc_generation", &pReco_particle_ccinc_generation);
        tree->SetBranchAddress("reco_particle_ccinc_contained", &pReco_particle_ccinc_isContained);
        tree->SetBranchAddress("passed_topologicalScoreCC", &passed_topologicalScoreCC);
        tree->SetBranchAddress("passed_max1Uncontained", &passed_max1Uncontained);

        Float_t spline_weight, tuned_cv_weight;
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);
        
        std::vector<float>* pReco_particle_ccinc_muonBDTScore = nullptr;
        tree->SetBranchAddress("reco_particle_ccinc_muonBDTScore", &pReco_particle_ccinc_muonBDTScore);

        const auto nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++)
        {
            tree->GetEntry(i);
            if (!passed_topologicalScoreCC) // !passed_max1Uncontained
            {
                continue;
            }

            auto eventWeight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;
            eventWeight *= runWeight; // Apply the run weight

            int nCountedParticles = 0;
            bool hasInvalidBDTScore = false;
            for (size_t v = 0; v < pReco_particle_ccinc_generation->size(); v++)
            {                
                if (pReco_particle_ccinc_generation->at(v) == 2 && pReco_particle_ccinc_isContained->at(v))
                {
                    float bdtScore = pReco_particle_ccinc_muonBDTScore->at(v);
                    if (bdtScore >= -1.f && bdtScore <= 1.f && std::isfinite(bdtScore))
                        hasInvalidBDTScore = true;

                    nCountedParticles++;
                }
            }

            if(!hasInvalidBDTScore)
                continue;

            histograms[sampleType]->Fill(nCountedParticles, eventWeight);
        }

        tFile->Close();
    }

    // Create a canvas and legend
    auto c = new TCanvas("multiplicity_canvas", "Particle Multiplicity", 800, 600);
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Create a THStack for the stacked histograms
    auto hs = new THStack("hs", "Particle Multiplicity;Number of Particles;Events");

    // Add histograms to the stack
    for (const auto& type : {"nu mc", "beam off", "dirt mc"})
    {
        std::string typeStr(type); // Convert to std::string
        histograms[typeStr]->SetFillColor(typeStr == "nu mc" ? kBlue : (typeStr == "beam off" ? kRed : kGreen));
        histograms[typeStr]->SetLineColor(kBlack);
        hs->Add(histograms[typeStr].get());
        legend->AddEntry(histograms[typeStr].get(), typeStr.c_str(), "f");
    }

    // Determine the maximum bin content across all histograms
    double maxBinContent = 0.0;
    for (const auto& [type, hist] : histograms)
    {
        double localMax = hist->GetMaximum();
        if (localMax > maxBinContent)
        {
            maxBinContent = localMax;
        }
    }

    // Scale the y-axis to ensure all data points are visible
    hs->SetMaximum(1.2 * maxBinContent); // Add a 20% margin above the maximum value

    // Draw the stack
    hs->Draw("HIST");
    hs->GetXaxis()->SetTitle("Number of Particles");
    hs->GetYaxis()->SetTitle("Events");

    // Overlay the beam on data
    histograms["beam on"]->SetMarkerStyle(20);
    histograms["beam on"]->SetMarkerColor(kBlack);
    histograms["beam on"]->Draw("E SAME");
    legend->AddEntry(histograms["beam on"].get(), "Beam On Data", "lep");

    // Draw the legend
    legend->Draw();

    // Save the canvas
    c->SaveAs("multiplicity_plot.pdf");
}