void makeAllPlots() {
    // Creates the cut-by-cut plots tables for all runs separately and together
    gROOT->ProcessLine(".x SelectionTable.C");
    
    // Creates 2d plots for the bdt input variables vs the bdt score for each run
    // Includes all events that pass the CC inclusive step
    gROOT->ProcessLine(".x BDTStudySigVsBkg.C");

    // Creates particle-level BDT score plots (including trained vs untrained nu overlay)
    // for all runs separately and together
    gROOT->ProcessLine(".x BDTStudy1DParticle.C");

    // Creates particle-level BDT input variable plots for all runs separately and together
    gROOT->ProcessLine(".x BDTStudy1DParticleInputVar.C");

    // Creates stacked particle-level BDT score plots for all runs separately and together
    // showing real data vs MC
    gROOT->ProcessLine(".x BDTStudy1DVsData.C");

    // Creates 2d pion vs muon momentum plots
    gROOT->ProcessLine(".x MuonPIDStudy.C");

    // Creates stacked histograms of the selection cuts
    gROOT->ProcessLine(".x CutPlots.C");

    // Creates multiple 2D reco vs true pion phi plots in different pion momentum regions
    gROOT->ProcessLine(".x PionAngleInMomentumBins.C");

    // MuonMomentumVSAngleContainementStudy.C

    // // Creates some plots that show compare the efficiency of the genric, contained muon and golden selection
    // gROOT->ProcessLine(".x mcEfficiencyPlots.C"); // Not used
}