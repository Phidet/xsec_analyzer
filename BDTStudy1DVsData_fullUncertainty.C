#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "MCC9SystematicsCalculator.hh"
#include "TMatrixD.h"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"

#include <map>
#include <string>
#include <memory>

#include "EventCategory.hh"

std::string toLatexScientific(double value) {
    std::stringstream stream;
    stream << std::scientific << std::setprecision(2) << value;
    std::string str = stream.str();
    size_t pos = str.find('e');
    if (pos != std::string::npos) {
        str.replace(pos, 1, " \\times 10^{");
        str += "}";
    }
    if (str.find('+') != std::string::npos) {
        str.replace(str.find('+'), 1, "");
    }
    return str;
}

void BDTStudy1DVsData_fullUncertainty() 
{
    const EventCategoryInterpreter& eventCategoryInterpreter = EventCategoryInterpreter::Instance();

    const std::string rootPath = "/exp/uboone/data/users/jdetje/ubcc1piPelee/1March24/";
    const std::string rootPathWithLengths = "/exp/uboone/data/users/jdetje/ubcc1piPelee/14May25_withRecoLengthAndUniversalVertDist/";

    const bool addUncertaintyBands = true; // Flag to toggle uncertainty bands
    const bool showRatioPlots = false;     // Flag to toggle showing ratio plots
  
    const bool useLongestTracks = false;    // Flag to toggle using only the longest tracks
    const unsigned int nLongestTracks = 2;  // Number of longest tracks to consider when useLongestTracks is true

    const bool limitToNParticleEvents = false; // If true, skip all particles with more or less gen=2 reco particles
    const int nParticleEvent = 5;

    const bool fractionalWeighting = false; // Flag to weight tracks as a fraction of the total number of tracks in event 

    const bool limitToAppliedParticles = true; // Only fill the BDT plots with the particles and events that the BDTs are applied to

    const bool useFirstParticle = false; // Only use the first particle in each event
    
    const bool showCutLines = true; // Flag to toggle showing cut lines on BDT plots

    const bool onlyChiNoP = false; // Only show the calculated chi^2 value and not the p-value

    const bool useRandomIndex = false; // Uses the randomIndex branch

    if(useFirstParticle && useRandomIndex)
    {
        // Throw an error that only one can be activated
        throw std::runtime_error("ERROR: useFirstParticle and useRandomIndex cannot be used together");
    }

    // std::string uniqName = "onlySelected_2July";
    std::string uniqName = "";
    // std::string uniqName = "noUncertainty";
    // std::string uniqName = "randomIndex";
    // std::string uniqName = "randomIndex_allQualityCuts_fixedTrackDistance";
    
    // // List of files; tuple: file type, run, path, fileWeight
    // const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
    //     // std::make_tuple("beam on", "1",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi.root", 1.f),
    //     // std::make_tuple("beam on", "2",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 1.f),
    //     // std::make_tuple("beam on", "3",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 1.f),
    //     // std::make_tuple("beam on", "4bcd", rootPath + "bnb_beam_on_peleeTuple_uboone_run4bcd_ubcc1pi.root", 1.f),
    //     // std::make_tuple("beam on", "5",  rootPathWithLengths + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root", 1.f),
    //     std::make_tuple("beam on", "5",  rootPath + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi.root", 1.f),

    //     // std::make_tuple("beam off", "1",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi.root", 0.56421),
    //     // std::make_tuple("beam off", "2",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi.root", 0.40202),
    //     // std::make_tuple("beam off", "3",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi.root", 0.30657),
    //     // std::make_tuple("beam off", "4bcd", rootPath + "bnb_beam_off_peleeTuple_uboone_run4bcd_ubcc1pi.root", 0.30175),
    //     // std::make_tuple("beam off", "5",  rootPathWithLengths + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root", 0.32807),
    //     std::make_tuple("beam off", "5",  rootPath + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi.root", 0.32807),

    //     // std::make_tuple("nu mc", "1",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
    //     // std::make_tuple("nu mc", "2",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing.root", 0.25750*2.0),
    //     // std::make_tuple("nu mc", "3",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing.root", 0.20113*2.0),
    //     // std::make_tuple("nu mc", "4bcd", "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing.root", 0.13074*2.0),
    //     // std::make_tuple("nu mc", "5",  rootPathWithLengths + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196*2.0), // recoVertex distance is defined differently (not dependent on passed_2NonProtons cut)
    //     std::make_tuple("nu mc", "5",  "/exp/uboone/data/users/jdetje/ubcc1piPelee/27March24_pionMomThreshold/overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing.root", 0.15196*2.0),

    //     // std::make_tuple("dirt mc", "1",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi.root", 0.52806),
    //     // std::make_tuple("dirt mc", "2",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi.root", 0.27521),
    //     // std::make_tuple("dirt mc", "3",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi.root", 0.80892),
    //     // std::make_tuple("dirt mc", "4bcd", rootPath + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi.root", 0.329301), // old value: 0.39701),
    //     // std::make_tuple("dirt mc", "5",  rootPathWithLengths + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root", 0.41280),
    //     std::make_tuple("dirt mc", "5",  rootPath + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi.root", 0.41280),        
    // };

    // const std::string rootPathWithHighestMuonBDTScore = "/pnfs/uboone/persistent/users/jdetje/ntupplesWithHighestMuonBDTScoreBranch/";
    // // List of files; tuple: file type, run, path, fileWeight
    // const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
    //     std::make_tuple("beam on", "1",  rootPathWithHighestMuonBDTScore + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi_highestMuonBDT.root", 1.f),
    //     std::make_tuple("beam on", "2",  rootPathWithHighestMuonBDTScore + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi_highestMuonBDT.root", 1.f),
    //     std::make_tuple("beam on", "3",  rootPathWithHighestMuonBDTScore + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi_highestMuonBDT.root", 1.f),
    //     std::make_tuple("beam on", "4bcd", rootPathWithHighestMuonBDTScore + "bnb_beam_on_peleeTuple_uboone_run4bcd_ubcc1pi_highestMuonBDT.root", 1.f),
    //     std::make_tuple("beam on", "5",  rootPathWithHighestMuonBDTScore + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi_highestMuonBDT.root", 1.f),

    //     std::make_tuple("beam off", "1",  rootPathWithHighestMuonBDTScore + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi_highestMuonBDT.root", 0.56421),
    //     std::make_tuple("beam off", "2",  rootPathWithHighestMuonBDTScore + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi_highestMuonBDT.root", 0.40202),
    //     std::make_tuple("beam off", "3",  rootPathWithHighestMuonBDTScore + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi_highestMuonBDT.root", 0.30657),
    //     std::make_tuple("beam off", "4bcd", rootPathWithHighestMuonBDTScore + "bnb_beam_off_peleeTuple_uboone_run4bcd_ubcc1pi_highestMuonBDT.root", 0.30175),
    //     std::make_tuple("beam off", "5",  rootPathWithHighestMuonBDTScore + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi_highestMuonBDT.root", 0.32807),

    //     std::make_tuple("nu mc", "1",  rootPathWithHighestMuonBDTScore + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing_highestMuonBDT.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
    //     std::make_tuple("nu mc", "2",  rootPathWithHighestMuonBDTScore + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing_highestMuonBDT.root", 0.25750*2.0),
    //     std::make_tuple("nu mc", "3",  rootPathWithHighestMuonBDTScore + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing_highestMuonBDT.root", 0.20113*2.0),
    //     std::make_tuple("nu mc", "4bcd", rootPathWithHighestMuonBDTScore + "overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing_highestMuonBDT.root", 0.13074*2.0),
    //     std::make_tuple("nu mc", "5",  rootPathWithHighestMuonBDTScore + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing_highestMuonBDT.root", 0.15196*2.0),

    //     std::make_tuple("dirt mc", "1",  rootPathWithHighestMuonBDTScore + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi_highestMuonBDT.root", 0.52806),
    //     std::make_tuple("dirt mc", "2",  rootPathWithHighestMuonBDTScore + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi_highestMuonBDT.root", 0.27521),
    //     std::make_tuple("dirt mc", "3",  rootPathWithHighestMuonBDTScore + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi_highestMuonBDT.root", 0.80892),
    //     std::make_tuple("dirt mc", "4bcd", rootPathWithHighestMuonBDTScore + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi_highestMuonBDT.root", 0.329301), // old value: 0.39701),
    //     std::make_tuple("dirt mc", "5",  rootPathWithHighestMuonBDTScore + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi_highestMuonBDT.root", 0.41280),
    // };

    const std::string rootPathWithRandomIndex = "/pnfs/uboone/persistent/users/jdetje/ntupplesWithRandomIndexBranch/";
    // List of files; tuple: file type, run, path, fileWeight
    const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
        std::make_tuple("beam on", "1",  rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),
        std::make_tuple("beam on", "2",  rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),
        std::make_tuple("beam on", "3",  rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),
        std::make_tuple("beam on", "4bcd", rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_run4bcd_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),
        std::make_tuple("beam on", "5",  rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),

        std::make_tuple("beam off", "1",  rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi_highestMuonBDT_randomIndex.root", 0.56421),
        std::make_tuple("beam off", "2",  rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi_highestMuonBDT_randomIndex.root", 0.40202),
        std::make_tuple("beam off", "3",  rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi_highestMuonBDT_randomIndex.root", 0.30657),
        std::make_tuple("beam off", "4bcd", rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_run4bcd_ubcc1pi_highestMuonBDT_randomIndex.root", 0.30175),
        std::make_tuple("beam off", "5",  rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi_highestMuonBDT_randomIndex.root", 0.32807),

        std::make_tuple("nu mc", "1",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_only_testing_highestMuonBDT_randomIndex.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
        std::make_tuple("nu mc", "2",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_only_testing_highestMuonBDT_randomIndex.root", 0.25750*2.0),
        std::make_tuple("nu mc", "3",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_only_testing_highestMuonBDT_randomIndex.root", 0.20113*2.0),
        std::make_tuple("nu mc", "4bcd", rootPathWithRandomIndex + "overlay_peleeTuple_uboone_run4bcd_nu_ubcc1pi_only_testing_highestMuonBDT_randomIndex.root", 0.13074*2.0),
        std::make_tuple("nu mc", "5",  rootPathWithRandomIndex + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_only_testing_highestMuonBDT_randomIndex.root", 0.15196*2.0),

        std::make_tuple("dirt mc", "1",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi_highestMuonBDT_randomIndex.root", 0.52806),
        std::make_tuple("dirt mc", "2",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi_highestMuonBDT_randomIndex.root", 0.27521),
        std::make_tuple("dirt mc", "3",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi_highestMuonBDT_randomIndex.root", 0.80892),
        std::make_tuple("dirt mc", "4bcd", rootPathWithRandomIndex + "overlay_peleeTuple_uboone_run4bcd_dirt_ubcc1pi_highestMuonBDT_randomIndex.root", 0.329301), // old value: 0.39701),
        std::make_tuple("dirt mc", "5",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi_highestMuonBDT_randomIndex.root", 0.41280),
    };

    // const std::string rootPathWithRandomIndex = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/highestMuonBDTScoreAndRandomIndex/";
    // const std::string rootPathWithRandomIndexTestEventsOnly = "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/17July2025_fixedTrackDistance/onlyTestEvents/";
    // // List of files; tuple: file type, run, path, fileWeight
    // const std::vector<std::tuple<std::string, std::string, std::string, float>> files = {
    //     std::make_tuple("beam on", "1",  rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),
    //     std::make_tuple("beam on", "2",  rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),
    //     std::make_tuple("beam on", "3",  rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),
    //     std::make_tuple("beam on", "4bcd", rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_v08_00_00_73_run4bcd_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),
    //     std::make_tuple("beam on", "5",  rootPathWithRandomIndex + "bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5_ubcc1pi_highestMuonBDT_randomIndex.root", 1.f),

    //     std::make_tuple("beam off", "1",  rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1_ubcc1pi_highestMuonBDT_randomIndex.root", 0.56421),
    //     std::make_tuple("beam off", "2",  rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2_ubcc1pi_highestMuonBDT_randomIndex.root", 0.40202),
    //     std::make_tuple("beam off", "3",  rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3_ubcc1pi_highestMuonBDT_randomIndex.root", 0.30657),
    //     std::make_tuple("beam off", "4bcd", rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_v08_00_00_73_run4bcd_ubcc1pi_highestMuonBDT_randomIndex.root", 0.30175),
    //     std::make_tuple("beam off", "5",  rootPathWithRandomIndex + "bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5_ubcc1pi_highestMuonBDT_randomIndex.root", 0.32807),

    //     std::make_tuple("nu mc", "1",  rootPathWithRandomIndexTestEventsOnly + "overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root", 0.13011*2.0), // Times two because the scaling is for the full MC and this is only half
    //     std::make_tuple("nu mc", "2",  rootPathWithRandomIndexTestEventsOnly + "overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root", 0.25750*2.0),
    //     std::make_tuple("nu mc", "3",  rootPathWithRandomIndexTestEventsOnly + "overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root", 0.20113*2.0),
    //     std::make_tuple("nu mc", "4bcd", rootPathWithRandomIndexTestEventsOnly + "overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_nu_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root", 0.13074*2.0),
    //     std::make_tuple("nu mc", "5",  rootPathWithRandomIndexTestEventsOnly + "overlay_nu_peleeTuple_uboone_v08_00_00_73_weightFix_run5_ubcc1pi_highestMuonBDT_randomIndex_onlyTestEvents.root", 0.15196*2.0),

    //     std::make_tuple("dirt mc", "1",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt_ubcc1pi_highestMuonBDT_randomIndex.root", 0.52806),
    //     std::make_tuple("dirt mc", "2",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt_ubcc1pi_highestMuonBDT_randomIndex.root", 0.27521),
    //     std::make_tuple("dirt mc", "3",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt_ubcc1pi_highestMuonBDT_randomIndex.root", 0.80892),
    //     std::make_tuple("dirt mc", "4bcd", rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run4bcd_dirt_ubcc1pi_highestMuonBDT_randomIndex.root", 0.329301), // old value: 0.39701),
    //     std::make_tuple("dirt mc", "5",  rootPathWithRandomIndex + "overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt_with_fake_weights_ubcc1pi_highestMuonBDT_randomIndex.root", 0.41280),
    // };

    const std::vector<std::string> runs = {"0", "1", "2", "3", "4bcd", "5"}; // Here 0 is the full set of all runs
    const std::vector<std::string> bdts = {"muonBDTScore", "protonBDTScore", "goldenPionBDTScore", "logBragg_pToMIP", "logBragg_piToMIP", "truncMeandEdx", "wiggliness", "trackScore", "nDescendents", "length"};
    // const std::vector<std::string> types = {"EXT", "Dirt", "E", "Photon", "K", "P", "Mu", "Pi+", "Unscattered Pi+", "Beam On"};
    const std::vector<std::string> types = {"Other", "P", "Mu", "Pi+", "Unscattered Pi+", "Beam On"};
    const std::vector<bool> dataOptions = {true, false};

    // map: variable, tuple(nBins, xMin, xMax, axisLog)
    const std::map<std::string, std::tuple<int, double, double, bool>> binningInfo = {
        // {"muonBDTScore", std::make_tuple(17, -0.95, 0.65, false)},
        // {"protonBDTScore", std::make_tuple(17, -0.85, 0.60, false)},
        // {"goldenPionBDTScore", std::make_tuple(17, -0.90, 0.50, false)},

        // limitToAppliedParticles = true
        {"muonBDTScore", std::make_tuple(25, -0.95, 0.65, false)},
        {"protonBDTScore", std::make_tuple(22, -0.83, 0.582, false)},
        {"goldenPionBDTScore", std::make_tuple(25, -0.705, 0.5, false)},

        // {"muonBDTScore", std::make_tuple(80, -1, 1, false)},
        // {"protonBDTScore", std::make_tuple(80, -1, 1, false)},
        // {"goldenPionBDTScore", std::make_tuple(80, -1, 1, false)},

        {"logBragg_pToMIP", std::make_tuple(60, -8, 8, false)},
        {"logBragg_piToMIP", std::make_tuple(60, -4, 7, false)},
        {"truncMeandEdx", std::make_tuple(60, 0, 10.0, false)},
        {"wiggliness", std::make_tuple(60, 0.00005, 0.25, true)},
        {"trackScore", std::make_tuple(60, 0, 1.0, false)},
        {"nDescendents", std::make_tuple(4, 0, 4, false)},
        {"length", std::make_tuple(50, 0, 300, false)}

        // {"muonBDTScore", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
        // {"protonBDTScore", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
        // {"goldenPionBDTScore", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
        // {"logBragg_pToMIP", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
        // {"logBragg_piToMIP", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
        // {"truncMeandEdx", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
        // {"wiggliness", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
        // {"trackScore", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
        // {"nDescendents", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
        // {"length", std::make_tuple(1, -std::numeric_limits<float>::max()/2, std::numeric_limits<float>::max()/2, false)},
    };

    float outOfBoundsMuonBDTScoresForData = 0;
    float outOfBoundsMuonBDTScoresForMC = 0;
    float totalMuonBDTScoresForData = 0;
    float totalMuonBDTScoresForMC = 0;

    // map: run, bdt, particle
    std::map<std::string, std::map<std::string, std::map<std::string, std::unique_ptr<TH1D>>>> histograms1D;
    for (const auto& run : runs) {
        for (const auto& bdt : bdts) {
            const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
            for (const auto& type : types) {
                    const auto name = "h_" + run + "_" + bdt + "_" + type;
                    if(xAxisLog)
                    {
                        std::cout << "DEBUG: Creating histogram with logarithmic binning for " << name << std::endl;
                        // Create a vector to hold the bin edges
                        std::vector<double> binEdges(nBinsBDT + 1);
                        double logMin = TMath::Log10(xMin);
                        double logMax = TMath::Log10(xMax);
                        double binWidth = (logMax - logMin) / nBinsBDT;
                        for (int i = 0; i <= nBinsBDT; i++) {
                            binEdges[i] = TMath::Power(10, logMin + i * binWidth);
                        }

                        // Create the histogram with logarithmic binning
                        histograms1D[run][bdt][type] = std::make_unique<TH1D>((name + "_" + bdt + "_" + type).c_str(), (name + " reco Particles Backtracked to True " + type + " in Signal Events; " + bdt + " BDT Score").c_str(), nBinsBDT, binEdges.data());
                    }
                    else
                    {
                        histograms1D[run][bdt][type] = std::make_unique<TH1D>((name + "_" + bdt + "_" + type).c_str(), (name + " reco Particles Backtracked to True " + type + " in Signal Events; " + bdt + " BDT Score").c_str(), nBinsBDT, xMin, xMax);
                    }
                    histograms1D[run][bdt][type]->Sumw2();
            }
        }
    }

    // Loop over the files
    for (int f = 0; f < files.size(); f++)
    {
        const auto [sampleType, run, filePath, runWeight] = files.at(f);

        const auto isBeamOn = (sampleType == "beam on");
        const auto isBeamOff = (sampleType == "beam off");
        const auto isNuMC = (sampleType == "nu mc");
        const auto isMC = (sampleType == "nu mc" || sampleType == "dirt mc");
        const auto isDirt = (sampleType == "dirt mc");

        // Open the file and get the tree
        TFile* tFile = TFile::Open(filePath.c_str());
        if (!tFile || tFile->IsZombie()) {
            std::cerr << "Error: Unable to open file: " << filePath << std::endl;
            return;
        }

        TTree* tree = (TTree*)tFile->Get("stv_tree");
        if (!tree) {
            std::cerr << "Error: Unable to get tree from file: " << filePath << std::endl;
            return;
        }

        std::vector<float>* pReco_particle_ccinc_muonBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_protonBDTScore = nullptr;
        std::vector<float>* pReco_particle_ccinc_goldenPionBDTScore = nullptr;
        std::vector<int>* pReco_particle_ccinc_generation = nullptr;
        std::vector<int>* pReco_particle_ccinc_backtracked_pdg = nullptr;
        std::vector<bool>* pReco_particle_ccinc_backtracked_goldenPion = nullptr;
        std::vector<bool>* pReco_particle_ccinc_isContained = nullptr;
        std::vector<float>* pReco_particle_ccinc_length = nullptr;

        std::vector<float>* pReco_particle_ccinc_logBragg_pToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_logBragg_piToMIP = nullptr;
        std::vector<float>* pReco_particle_ccinc_truncMeandEdx = nullptr;
        std::vector<float>* pReco_particle_ccinc_wiggliness = nullptr;
        std::vector<float>* pReco_particle_ccinc_trackScore = nullptr;
        std::vector<int>* pReco_particle_ccinc_nDescendents = nullptr;

        tree->SetBranchAddress("reco_particle_ccinc_muonBDTScore", &pReco_particle_ccinc_muonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_protonBDTScore", &pReco_particle_ccinc_protonBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_goldenPionBDTScore", &pReco_particle_ccinc_goldenPionBDTScore);
        tree->SetBranchAddress("reco_particle_ccinc_generation", &pReco_particle_ccinc_generation);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_pdg", &pReco_particle_ccinc_backtracked_pdg);
        tree->SetBranchAddress("reco_particle_ccinc_backtracked_goldenPion", &pReco_particle_ccinc_backtracked_goldenPion);
        tree->SetBranchAddress("reco_particle_ccinc_contained", &pReco_particle_ccinc_isContained);
        // tree->SetBranchAddress("reco_particle_ccinc_length", &pReco_particle_ccinc_length);
        
        // Set length to dummy value since only run 5 samples have it
        pReco_particle_ccinc_length = new std::vector<float>(pReco_particle_ccinc_generation->size(), 0.0f); // Initialize with zeros

        tree->SetBranchAddress("reco_particle_ccinc_logBragg_pToMIP", &pReco_particle_ccinc_logBragg_pToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_logBragg_piToMIP", &pReco_particle_ccinc_logBragg_piToMIP);
        tree->SetBranchAddress("reco_particle_ccinc_truncMeandEdx", &pReco_particle_ccinc_truncMeandEdx);
        tree->SetBranchAddress("reco_particle_ccinc_wiggliness", &pReco_particle_ccinc_wiggliness);
        tree->SetBranchAddress("reco_particle_ccinc_trackScore", &pReco_particle_ccinc_trackScore);
        tree->SetBranchAddress("reco_particle_ccinc_nDescendents", &pReco_particle_ccinc_nDescendents);

        Float_t spline_weight, tuned_cv_weight;
        tree->SetBranchAddress("spline_weight", &spline_weight);
        tree->SetBranchAddress("tuned_cv_weight", &tuned_cv_weight);

        // std::string metaType = isBeamOn ? "Beam On" : (isDirt ? "Dirt" : (isBeamOff ? "EXT" : ""));
        std::string metaType = isBeamOn ? "Beam On" : ((isDirt || isBeamOff) ? "Other" : "");

        // Get the total number of entries
        // std::cout<<"WARNING - Only using 3% of the entries!!!!!!!!!"<<std::endl;
        const auto nEntries = tree->GetEntries();

        // Loop over the entries in the tree
        for (Long64_t i = 0; i < nEntries; i++)
        {
            tree->GetEntry(i);

            // Update the progress bar at every percent
            if (i % (nEntries / 100) == 0)
            {
                const auto progress = static_cast<int>(100.0 * i / nEntries);
                std::cout << "\r"<< " " << sampleType << " run " << run <<"[" << std::string(progress, '|') << std::string(100 - progress, ' ') << "] " << progress << "%" << std::flush;
            }

            // Apply the condition
            // if (tree->GetLeaf("passed_topologicalScoreCC")->GetValue())
            // if (tree->GetLeaf("passed_openingAngle")->GetValue())
            // if (tree->GetLeaf("passed_max1Uncontained")->GetValue())
            // if (tree->GetLeaf("passed_max1Uncontained")->GetValue() && tree->GetLeaf("event_cutValue_topologicalScore")->GetValue() >= 0.67)
            if (tree->GetLeaf("passed_max1Uncontained")->GetValue() && tree->GetLeaf("event_cutValue_topologicalScore")->GetValue() >= 0.67 && tree->GetLeaf("event_cutValue_maxVertexDist")->GetValue() <= 9.5*9.5 && tree->GetLeaf("event_cutValue_maxVertexDist")->GetValue() >= 0)
            // if (tree->GetLeaf("cc1pi_selected_generic")->GetValue() && tree->GetLeaf("cc1pi_reco_pionMomentum")->GetValue() > 0.1)
            {
                // const auto isTrainingEvent = tree->GetLeaf("isTrainingEvent")->GetValue();
                // if(isTrainingEvent)
                //     continue;

                // Count the number of generation 2 particles
                int nCountedParticles = 0;
                for (Long64_t v = 0; v < pReco_particle_ccinc_generation->size(); v++)
                {
                    if (pReco_particle_ccinc_generation->at(v) == 2) // && pReco_particle_ccinc_isContained->at(v))
                    {
                        nCountedParticles++;
                    }
                }
 
                float fractionalWeight = 1.0;
                if(fractionalWeighting)
                {                    
                    // Calculate the fractional weight
                    if (nCountedParticles > 0) {
                        fractionalWeight = 1.0 / nCountedParticles;
                    }
                }

                if( limitToNParticleEvents && nCountedParticles != nParticleEvent) continue;

                auto eventWeight = std::isfinite(spline_weight*tuned_cv_weight) && spline_weight*tuned_cv_weight >= 0 && spline_weight*tuned_cv_weight <= 30 ? spline_weight*tuned_cv_weight : 1;                
                eventWeight *= runWeight; // Apply the run weight
                eventWeight *= fractionalWeight;

                // if(isNuMC)
                //     eventWeight *= 2; // Compensate for skipping the BDT training events and only using 50% of the sample

                // If we're using longest tracks, identify the N longest contained generation 2 tracks
                std::vector<std::pair<float, Long64_t>> trackLengths; // Pair of (length, index)
                std::vector<Long64_t> selectedIndices;
            

                if (useLongestTracks) {
                    for (Long64_t v = 0; v < pReco_particle_ccinc_generation->size(); v++) {
                        if (pReco_particle_ccinc_generation->at(v) == 2){ // } && pReco_particle_ccinc_isContained->at(v)) {
                            trackLengths.emplace_back(pReco_particle_ccinc_length->at(v), v);
                        }
                    }
                    
                    // Sort tracks by length in descending order
                    std::sort(trackLengths.rbegin(), trackLengths.rend());
                    
                    // Select the top N longest tracks
                    for (size_t i = 0; i < std::min(nLongestTracks, static_cast<unsigned int>(trackLengths.size())); i++) {
                        selectedIndices.push_back(trackLengths[i].second);
                    }
                }

                // // If we're limiting to applied particles, we need to get the index of the particle with the lowest muon BDT score for fully contained events
                // unsigned int highestMuonBDTIndex = std::numeric_limits<unsigned int>::max();
                // if(limitToAppliedParticles && tree->GetLeaf("event_cutValue_nUncontained")->GetValue() == 0)
                // {
                //     float highestMuonBDTScore = -std::numeric_limits<float>::max();
                //     for (unsigned int v = 0; v < pReco_particle_ccinc_generation->size(); ++v)
                //     {
                //         // std::cout << "\n DEBUG: Event:"<< i <<": pReco_particle_ccinc_generation->at("<<v<<"): " << pReco_particle_ccinc_generation->at(v) << std::endl;
                //         if (pReco_particle_ccinc_generation->at(v) == 2)
                //         {
                //             float score = pReco_particle_ccinc_muonBDTScore->at(v);
                //             // std::cout << "\n DEBUG: Event:"<< i <<": pReco_particle_ccinc_muonBDTScore->at("<<v<<"): " << score << std::endl;
                //             if (score > highestMuonBDTScore && score > -1.f && score < 1.f)
                //             {
                //                 highestMuonBDTScore = score;
                //                 highestMuonBDTIndex = v;
                //             }
                //         }
                //     }
                // }

                unsigned int nPionValues = 0;
                unsigned int excludedMuonValues = 0;
                // If we're using the random index, get a random index
                int randomIndex = -2; // Events with index -1 have no particles
                if (useRandomIndex) {
                    randomIndex = static_cast<int>(tree->GetLeaf("random_index")->GetValue());
                }


                for (Long64_t v = 0; v< pReco_particle_ccinc_generation->size(); v++)
                {
                    if((useFirstParticle && v>0) || (useRandomIndex && randomIndex==-1)) break;
                    if (useRandomIndex && v != randomIndex) continue;
                    
                    // Set the type to the type determined by the file first
                    std::string type = metaType; // Only when that is not set, use the particle type
                    
                    const auto generation = pReco_particle_ccinc_generation->at(v);
                    const auto isContained = pReco_particle_ccinc_isContained->at(v);
                    if(pReco_particle_ccinc_generation->size() != pReco_particle_ccinc_isContained->size())
                        throw std::runtime_error("ERROR: pReco_particle_ccinc_generation and pReco_particle_ccinc_isContained have different sizes");

                    if(generation == 2 && isContained)
                    {
                        // Skip if using longest tracks and this isn't one of the N longest
                        if (useLongestTracks && std::find(selectedIndices.begin(), selectedIndices.end(), v) == selectedIndices.end()) {
                            continue;
                        }

                        if(type == "")
                        {
                            const auto pdg = std::abs(pReco_particle_ccinc_backtracked_pdg->at(v));

                            switch (pdg)
                            {
                                case 13:
                                    type = "Mu";
                                    break;
                                case 211:
                                {
                                    const auto isUnscattered = pReco_particle_ccinc_backtracked_goldenPion->at(v);
                                    type = isUnscattered ? "Unscattered Pi+" : "Pi+";
                                    break;
                                }
                                case 2212:
                                    type = "P";
                                    break;
                                default:
                                    type = "Other";
                            }
                        }

                        
                        if (type != "")
                        {

                            if(!limitToAppliedParticles || tree->GetLeaf("event_cutValue_nUncontained")->GetValue() == 0)
                            {
                                histograms1D.at(run).at("muonBDTScore").at(type)->Fill(pReco_particle_ccinc_muonBDTScore->at(v), eventWeight);
                                histograms1D.at("0").at("muonBDTScore").at(type)->Fill(pReco_particle_ccinc_muonBDTScore->at(v), eventWeight);
                            }

                            // If limited to applied particles, only fill the proton BDT score with the scores of the non-muons
                            // we are already limited to contained particles, which excludes the exiting muons. 
                            // For fully contained events, we also exclude the particle with the lowest muon BDT score                            
                            if(!limitToAppliedParticles|| (tree->GetLeaf("event_cutValue_nUncontained")->GetValue() == 1 || (tree->GetLeaf("event_cutValue_nUncontained")->GetValue() == 0 
                            && std::abs(tree->GetLeaf("reco_event_ccinc_highestMuonBDTScore")->GetValue()) <= 1
                            && std::abs(tree->GetLeaf("reco_event_ccinc_highestMuonBDTScore")->GetValue() - pReco_particle_ccinc_muonBDTScore->at(v)) > 0.00000001 )))
                            {
                                histograms1D.at(run).at("protonBDTScore").at(type)->Fill(pReco_particle_ccinc_protonBDTScore->at(v), eventWeight);
                                histograms1D.at("0").at("protonBDTScore").at(type)->Fill(pReco_particle_ccinc_protonBDTScore->at(v), eventWeight);
                            }
                            // else
                            // {
                            //     std::cout << "DEBUG: Not filling protonBDTScore for v=" << v << " in event " << i << std::endl;
                            //     std::cout << "  event_cutValue_nUncontained = " << tree->GetLeaf("event_cutValue_nUncontained")->GetValue() << std::endl;
                            //     std::cout << "  reco_event_ccinc_highestMuonBDTScore = " << tree->GetLeaf("reco_event_ccinc_highestMuonBDTScore")->GetValue() << std::endl;
                            //     std::cout << "  pReco_particle_ccinc_muonBDTScore->at(v) = " << pReco_particle_ccinc_muonBDTScore->at(v) << std::endl;
                            //     std::cout << "  fabs(reco_event_ccinc_highestMuonBDTScore) <= 1: " << (std::abs(tree->GetLeaf("reco_event_ccinc_highestMuonBDTScore")->GetValue()) <= 1) << std::endl;
                            //     std::cout << "  fabs(reco_event_ccinc_highestMuonBDTScore - pReco_particle_ccinc_muonBDTScore->at(v)) > 1e-8: " 
                            //               << (std::abs(tree->GetLeaf("reco_event_ccinc_highestMuonBDTScore")->GetValue() - pReco_particle_ccinc_muonBDTScore->at(v)) > 0.00000001) << std::endl;
                            // }


                            if(tree->GetLeaf("event_cutValue_nUncontained")->GetValue() == 0 
                            && std::abs(tree->GetLeaf("reco_event_ccinc_highestMuonBDTScore")->GetValue()) <= 1
                            && std::abs(tree->GetLeaf("reco_event_ccinc_highestMuonBDTScore")->GetValue() - pReco_particle_ccinc_muonBDTScore->at(v)) <= 0.00000001 )
                            {
                                excludedMuonValues++;
                                if(excludedMuonValues > 1)
                                {
                                    std::cout << "Event " << i << " has multiple muon BDT values with limitToAppliedParticles = true:" << std::endl;
                                    std::cout << "All muonBDTScore values for contained generation 2 particles in this event:" << std::endl;
                                    for (unsigned int vv = 0; vv < pReco_particle_ccinc_generation->size(); ++vv)
                                    {
                                        if (pReco_particle_ccinc_generation->at(vv) == 2 && pReco_particle_ccinc_isContained->at(vv))
                                        {
                                            std::cout << "  Index " << vv << ": muonBDTScore = " << pReco_particle_ccinc_muonBDTScore->at(vv) << std::endl;
                                        }
                                    }
                                    std::cout << "reco_event_ccinc_highestMuonBDTScore = " << tree->GetLeaf("reco_event_ccinc_highestMuonBDTScore")->GetValue() << std::endl;
                                    throw std::runtime_error("ERROR: More than one muon value added to muonBDTScore despite limitToAppliedParticles = true");
                                }
                            }

                            // if(tree->GetLeaf("event_cutValue_nUncontained")->GetValue() == 0 && highestMuonBDTIndex == std::numeric_limits<unsigned int>::max())
                            // {
                            //     throw std::runtime_error("ERROR: highestMuonBDTIndex is not set");
                            // }

                            // If limited to applied particles, only fill the golden BDT score with the scores of the selected pion candidates
                            if(!limitToAppliedParticles 
                                || (tree->GetLeaf("passed_openingAngle")->GetValue() 
                                    && std::abs(tree->GetLeaf("event_cutValue_goldenPionBDT")->GetValue()) <= 1.f 
                                    && tree->GetLeaf("cc1pi_reco_pionMomentum")->GetValue() >= 0.1
                                    && std::abs(pReco_particle_ccinc_goldenPionBDTScore->at(v) - tree->GetLeaf("event_cutValue_goldenPionBDT")->GetValue()) < 0.00000001)
                            )
                            {
                                nPionValues++;
                                // std::cout << "\n DEBUG: Event:"<< i <<": pReco_particle_ccinc_goldenPionBDTScore->at("<<v<<"): " << pReco_particle_ccinc_goldenPionBDTScore->at(v) << " event_cutValue_goldenPionBDT: " << tree->GetLeaf("event_cutValue_goldenPionBDT")->GetValue() << std::endl;
                                if(limitToAppliedParticles && nPionValues > 1) throw std::runtime_error("ERROR: More than one golden value added to goldenPionBDTScore despite limitToAppliedParticles = true");

                                histograms1D.at(run).at("goldenPionBDTScore").at(type)->Fill(pReco_particle_ccinc_goldenPionBDTScore->at(v), eventWeight);
                                histograms1D.at("0").at("goldenPionBDTScore").at(type)->Fill(pReco_particle_ccinc_goldenPionBDTScore->at(v), eventWeight);
                            }

                            histograms1D.at(run).at("logBragg_pToMIP").at(type)->Fill(pReco_particle_ccinc_logBragg_pToMIP->at(v), eventWeight);
                            histograms1D.at(run).at("logBragg_piToMIP").at(type)->Fill(pReco_particle_ccinc_logBragg_piToMIP->at(v), eventWeight);
                            histograms1D.at(run).at("truncMeandEdx").at(type)->Fill(pReco_particle_ccinc_truncMeandEdx->at(v), eventWeight);
                            histograms1D.at(run).at("wiggliness").at(type)->Fill(pReco_particle_ccinc_wiggliness->at(v), eventWeight);
                            histograms1D.at(run).at("trackScore").at(type)->Fill(pReco_particle_ccinc_trackScore->at(v), eventWeight);
                            histograms1D.at(run).at("nDescendents").at(type)->Fill(pReco_particle_ccinc_nDescendents->at(v), eventWeight);
                            // Temporarily disabled since only in run 5 variables
                            // histograms1D.at(run).at("length").at(type)->Fill(pReco_particle_ccinc_length->at(v), eventWeight);

                            histograms1D.at("0").at("logBragg_pToMIP").at(type)->Fill(pReco_particle_ccinc_logBragg_pToMIP->at(v), eventWeight);
                            histograms1D.at("0").at("logBragg_piToMIP").at(type)->Fill(pReco_particle_ccinc_logBragg_piToMIP->at(v), eventWeight);
                            histograms1D.at("0").at("truncMeandEdx").at(type)->Fill(pReco_particle_ccinc_truncMeandEdx->at(v), eventWeight);
                            histograms1D.at("0").at("wiggliness").at(type)->Fill(pReco_particle_ccinc_wiggliness->at(v), eventWeight);
                            histograms1D.at("0").at("trackScore").at(type)->Fill(pReco_particle_ccinc_trackScore->at(v), eventWeight);
                            histograms1D.at("0").at("nDescendents").at(type)->Fill(pReco_particle_ccinc_nDescendents->at(v), eventWeight);
                            // Temporarily disabled since only in run 5 variables
                            // histograms1D.at("0").at("length").at(type)->Fill(pReco_particle_ccinc_length->at(v), eventWeight);


                            if(type == "Beam On")
                            {
                                totalMuonBDTScoresForData += eventWeight;
                            }
                            else
                            {
                                totalMuonBDTScoresForMC += eventWeight;
                            }
                            if(!(pReco_particle_ccinc_muonBDTScore->at(v) >= -1 && pReco_particle_ccinc_muonBDTScore->at(v) <= 1))
                            {

                                if(type == "Beam On")
                                {
                                    outOfBoundsMuonBDTScoresForData += eventWeight;
                                }
                                else
                                {
                                    outOfBoundsMuonBDTScoresForMC += eventWeight;
                                }
                            }
                        }
                        else
                        {
                            throw std::runtime_error("ERROR: type is empty");
                        }
                    } // End of if contained and generation 2
                } // End of loop over vector entries (aka types)
            } // End of if signal
        } // End of loop over tree entries
    } // End of loop over files

    std::cout << "\nDEBUG: Finished filling histograms" << std::endl;

    // Create a particle to color map
    std::map<std::string, int> particleColors = {
        {"P", kOrange+1},
        {"Mu", kBlue},
        {"Pi+", kGreen+3},
        // {"E", kCyan},
        // {"Photon", kYellow},
        // {"K", kMagenta},
        {"Unscattered Pi+", kGreen},
        // {"EXT", kBlack},
        // {"Dirt", kOrange-6},
        {"Beam On", kBlack},
        {"Other", kGray}
    };

    struct inputFiles
    {
      std::string rootFile;
      std::string fileList;
      std::string config;
      std::string sliceConfig;
      std::string nameExtension;
    };

    // std::cout << "DEBUG: Point T7" << std::endl;
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_output_bnb_run1234bcd5_5Feb25_bdts_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_removedZeroBins.txt",
    //     "_bdts_full_uncertainty"
    // };

    // // Coarser binning and all quality cuts applied
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_20May25_bdts_all_quality_cuts_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts.txt",
    //     "_bdts_full_uncertainty"
    // };

    // // Even coarser binning and all quality cuts applied
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_21May25_bdts_all_quality_cuts_17bins_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_17bins.txt",
    //     "_bdts_full_uncertainty"
    // };

    // Coarser binning (14) and all quality cuts applied
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_27May25_bdts_all_quality_cuts_14bins_changedProtonBDTBinning_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_14bins.txt",
    //     "_bdts_full_uncertainty"
    // };

    // // Coarser binning (17) and all quality cuts applied + only including events and particles that the BDTs are applied to 
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_28May25_bdts_all_quality_cuts_17bins_onlyAppliedParticles_full_uncertainty.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_17bins_highestMuonBDTScore.txt",
    //     "_bdts_full_uncertainty_highestMuonBDTScore"
    // };

    // Coarser binning (17) and all quality cuts applied + only including events and particles that the BDTs are applied to (fixed version)
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_28May25_bdts_all_quality_cuts_17bins_onlyAppliedParticles_fixed.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_17bins_highestMuonBDTScore.txt",
    //     "_bdts_full_uncertainty_highestMuonBDTScore"
    // };

    // // Coarser binning (15) and all quality cuts applied + only including events and particles that the BDTs are applied to (fixed version)
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_29May25_bdts_all_quality_cuts_15bins_onlyAppliedParticles_fixed.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_15bins_highestMuonBDTScore.txt",
    //     "_bdts_full_uncertainty_highestMuonBDTScore"
    // };

    // // Coarser binning (15) only including events and particles that the BDTs are applied to (fixed version) + only muon BDT plot using only the first particle in each event
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_30May25_bdts_all_quality_cuts_15bins_onlyAppliedParticles_randomSampling.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_15bins_randomSampling.txt",
    // };


    // 22 and 25 bin binning to align with the cut lines <--------------- Current main input
    inputFiles input{
        "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_30May25_bdts_all_quality_cuts_22And25bins_onlyAppliedParticles_fixed.root",
        "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore.txt",
        "systcalc.conf",
        "ubcc1pi_slice_config_bdt_all_quality_cuts_22And25bins_highestMuonBDTScore.txt",
        "_bdts_full_uncertainty_highestMuonBDTScore"
    };

    // // 22 and 25 bin binning to align with the cut lines; random sampling of particles
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_16July25_bdts_all_quality_cuts_22And25bins_onlyAppliedParticles_randomIndex.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore_randomIndex.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_22And25bins_highestMuonBDTScore.txt",
    //     "_bdts_full_uncertainty_highestMuonBDTScore"
    // };


    // // 22 and 25 bin binning to align with the cut lines; random sampling of particles; includes all later quality cuts.
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_17July25_bdts_all_quality_cuts_22And25bins_onlyAppliedParticles_randomIndex.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore_randomIndex.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_22And25bins_highestMuonBDTScore.txt",
    //     "_bdts_full_uncertainty_highestMuonBDTScore"
    // };

    // // 22 and 25 bin binning to align with the cut lines; random sampling of particles; includes all later quality cuts; fixed vertex distance variable
    // inputFiles input{
    //     "/exp/uboone/data/users/jdetje/ubcc1pi_univmake/22Feb24/univmake_bdt_21July25_bdts_all_quality_cuts_22And25bins_onlyAppliedParticles_randomIndex.root",
    //     "file_properties_testingOnly_lowPiMomThreshold_fullDetvars_highestMuonBDTScore_randomIndex_fixedTrackDistance.txt",
    //     "systcalc.conf",
    //     "ubcc1pi_slice_config_bdt_all_quality_cuts_22And25bins_highestMuonBDTScore.txt",
    //     "_bdts_full_uncertainty_highestMuonBDTScore_fixedTrackDistance"
    // };

    auto &fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties(input.fileList);
    auto* syst_ptr = new MCC9SystematicsCalculator(input.rootFile, input.config);
    syst_ptr->set_syst_mode(syst_ptr->SystMode::VaryBackgroundAndSignalDirectly);
    std::string nameExtension = input.nameExtension;
    auto& syst = *syst_ptr;
    auto* sb_ptr = new SliceBinning( input.sliceConfig );
    auto& sb = *sb_ptr;
    auto* matrix_map_ptr = syst.get_covariances().release();
    auto& matrix_map = *matrix_map_ptr;


    // Get the factors needed to convert to cross-section units
    const double total_pot = syst.total_bnb_data_pot_;

    // std::cout << "DEBUG: Point T8" << std::endl;

    TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
    TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
    // Total MC+EXT prediction in reco bin space. Start by getting EXT.
    TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
        reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
    reco_mc_plus_ext_hist->SetDirectory( nullptr );

    // Add in the CV MC prediction
    reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

    std::cout << "DEBUG: Finished setting particle colors" << std::endl;
    float maxYScale = 1.2;
    for (const auto& run : runs)
    {
        // std::cout << "DEBUG Point Y - run " << run << std::endl;
        for (const auto& bdt : bdts)
        {
            // std::cout << "DEBUG Point X - bdt " << bdt << std::endl;
            const auto [nBinsBDT, xMin, xMax, xAxisLog] = binningInfo.at(bdt);
            // std::cout << "DEBUG Point Z0" << std::endl;
            const auto outName = "BDTStudy1DVsData_fullUncertainty_run" + run + "_" + bdt + "_" + uniqName;
            auto c = new TCanvas(outName.c_str(), "", 800, 600);
            c->SetLogx(xAxisLog);
            // Set up legend position and style using a switch statement for BDT variable
            TLegend *legend = nullptr;
            if (bdt == "trackScore") {
                legend = new TLegend(0.45, 0.6, 0.6, 0.9);
            } else if (bdt == "protonBDTScore") {
                legend = new TLegend(0.15, 0.6, 0.40, 0.9);
                // maxYScale = 1.35;
                maxYScale = 1.2;
            } else if (bdt == "muonBDTScore") {
                // legend = new TLegend(0.55, 0.6, 0.8, 0.9);
                legend = new TLegend(0.50, 0.6, 0.75, 0.9);
                maxYScale = 1.05;
            } else if (bdt == "goldenPionBDTScore") {
                legend = new TLegend(0.15, 0.6, 0.40, 0.9);
                // maxYScale = 1.35;
                maxYScale = 1.4;
            } else {
                legend = new TLegend(0.75, 0.6, 0.9, 0.9);
            }

            // Create pad layout based on showRatioPlots flag
            TPad *pad1;
            if (showRatioPlots) {
                pad1 = new TPad((outName+"pad1").c_str(), "pad1", 0, 0.3, 1, 0.95);
                pad1->SetBottomMargin(0.05); // Upper pad, no space on the bottom
            } else {
                pad1 = new TPad((outName+"pad1").c_str(), "pad1", 0, 0.05, 1, 0.95);
                pad1->SetBottomMargin(0.15); // Bottom margin for x-axis labels when no ratio plot
            }
            pad1->SetLeftMargin(0.12); // Increase left margin to prevent y-axis label cutoff
            pad1->Draw();             // Draw the upper pad
            pad1->cd();               // pad1 becomes the current pad
            gPad->SetLogx(xAxisLog);

            // Create a THStack
            auto hs = new THStack(outName.c_str(), "");

            // Create an initially empty histogram
            TH1D *sum = (TH1D*)histograms1D.at(run).at(bdt).at(types[0]).get()->Clone();
            sum->Reset(); // Reset all bins to 0

            float maxYTotal = 0;
            // std::cout << "DEBUG Point Z1" << std::endl;
            for (const auto& type : types)
            {
                // std::cout << "DEBUG Point Z2" << std::endl;
                auto h = histograms1D.at(run).at(bdt).at(type).get();
                h->SetStats(0); // Remove stats box
                const auto maxY = h->GetMaximum();
                if(maxY > maxYTotal)
                    maxYTotal = maxY;
                if(type == "Beam On")
                {
                    legend->AddEntry(h, "Beam On", "l");
                    continue;
                }

                if(type == "Other")//type == "EXT")
                    eventCategoryInterpreter.set_ext_histogram_style(h);
                else
                    h->SetFillColor(particleColors[type]);

                h->SetLineColor(kWhite); // Add this line to remove the outline
                h->SetLineWidth(0); // Add this line to remove the outline
                // Add the histogram to the stack
                hs->Add(h);

                // Add the histogram to the sum
                sum->Add(h);

                // Add an entry to the legend for this histogram
                const std::string typeName = type == "Other" ? "Other" : (type == "Beam On" ? "Beam On" : (type == "P" ? "Proton" : (type == "Mu" ? "Muon" : (type == "Pi+" ? "Scattered Charged Pion" : "Unscattered Charged Pion"))));
                legend->AddEntry(h, typeName.c_str(), "f");
            }

            std::ostringstream legendStream; // For chi2 and p-value

            // Replace the uncertainty block with the following, applied only for run "0"
            if(run=="0" && (bdt=="muonBDTScore" || bdt=="protonBDTScore" || bdt=="goldenPionBDTScore") && addUncertaintyBands)
            {
                const auto& slice = (bdt=="muonBDTScore") ? sb.slices_[0] : 
                                     (bdt=="protonBDTScore") ? sb.slices_[1] : 
                                     sb.slices_[2];
                SliceHistogram* slice_sum = SliceHistogram::make_slice_histogram(*sum, slice, &((*syst.get_covariances())["total"]));
                for (int bin = 1; bin <= sum->GetNbinsX(); bin++)
                {
                    double fullErr = slice_sum->hist_->GetBinError(bin);
                    sum->SetBinError(bin, fullErr);
                }

                // After applying uncertainties, check if any bin content + error exceeds current maxYTotal
                for (int bin = 1; bin <= sum->GetNbinsX(); bin++) {
                    double binContentPlusError = sum->GetBinContent(bin) + sum->GetBinError(bin);
                    if (binContentPlusError > maxYTotal) {
                        maxYTotal = binContentPlusError;
                    }
                }

                SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
                    *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );
            
                SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
                    *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

                auto chi2_result = slice_mc_plus_ext->get_chi2(*slice_bnb);

                legendStream << "MicroBooNE in the BNB, " << toLatexScientific(total_pot) << " POT, "
                            << "#chi^{2} = " << std::fixed << std::setprecision(2) << chi2_result.chi2_
                            << " / " << chi2_result.num_bins_ << " bin";
                if (chi2_result.num_bins_ > 1) legendStream << "s";
                if(!onlyChiNoP) legendStream << ", p = " << std::fixed << std::setprecision(2) << chi2_result.p_value_;

                delete slice_sum;
            }

            // std::cout << "DEBUG Point Z3" << std::endl;

            hs->SetMaximum(maxYTotal*maxYScale);
            
            hs->Draw("HIST");
            // auto hsCopy = (THStack*)hs->Clone();
            // sum->SetFillColor(TColor::GetColorTransparent(kBlack, 0.5));
            sum->SetLineWidth(3);
            sum->SetFillColor(kGray + 1);
            sum->SetFillStyle(3244);
            sum->Draw("E2 SAME");

            // Get the beam on histogram
            auto h = histograms1D.at(run).at(bdt).at("Beam On").get();
            eventCategoryInterpreter.set_bnb_data_histogram_style(h);
            h->SetStats(0); // Remove stats box
            h->Draw("E same");

            legend->Draw();

            // Set y axis title
            hs->GetYaxis()->SetTitle("# Reconstructed Particles");

            // Set title
            const std::string bdtTitle = bdt == "muonBDTScore" ? "Muon BDT Score" : (bdt == "protonBDTScore" ? "Proton BDT Score" : (bdt == "goldenPionBDTScore" ? "Unscattered Pion BDT Score" : (bdt == "logBragg_pToMIP" ? "log(R_{p}/MIP)" : (bdt == "logBragg_piToMIP" ? "log(R_{#pi}/MIP)" : (bdt == "truncMeandEdx" ? "Truncated Mean dE/dx" : (bdt == "wiggliness" ? "Wiggliness" : (bdt == "trackScore" ? "Track Score" : (bdt == "nDescendents" ? "Number of Descendents" : "Track Length (cm)"))))))));
            const std::string runTitle = run == "0" ? "All Runs" : ("Run " + run);
            hs->SetTitle((bdtTitle + " for " + runTitle + "\nReconstructed & Contained Beam-on and MC+EXT Particles in Events that Pass the CC #nu_#mu Preselection").c_str());
            
            // If we're not showing ratio plots, set the x-axis title here
            if (!showRatioPlots) {
                hs->GetXaxis()->SetTitle(bdtTitle.c_str());
                hs->GetXaxis()->SetTitleSize(0.045);
                hs->GetXaxis()->SetLabelSize(0.045);
                hs->GetXaxis()->SetTitleOffset(1.1);
            }
            
            // Add cut lines and labels if enabled
            if (showCutLines) {
                if (bdt == "protonBDTScore") {
                    pad1->cd(); // Make sure we're in the top pad
                    
                    // Draw vertical dashed line at cut value -0.06
                    const float protonCutValue = -0.06;
                    TLine *protonCutLine = new TLine(protonCutValue, 0, protonCutValue, maxYTotal*maxYScale);
                    protonCutLine->SetLineColor(kBlack);
                    protonCutLine->SetLineStyle(2); // Dashed line
                    protonCutLine->SetLineWidth(2);
                    protonCutLine->Draw();
                    
                    // Draw text label with right-pointing triangle
                    TLatex *protonLabel = new TLatex();
                    protonLabel->SetTextSize(0.04);
                    protonLabel->SetTextAlign(12); // Left aligned
                    protonLabel->DrawLatex(protonCutValue + 0.01, maxYTotal*maxYScale*0.95, "Proton candidates");
                }
                else if (bdt == "goldenPionBDTScore") {
                    pad1->cd(); // Make sure we're in the top pad
                    
                    // Draw vertical dashed line at cut value -0.03
                    const float pionCutValue = -0.03;
                    TLine *pionCutLine = new TLine(pionCutValue, 0, pionCutValue, maxYTotal*maxYScale);
                    pionCutLine->SetLineColor(kBlack);
                    pionCutLine->SetLineStyle(2); // Dashed line
                    pionCutLine->SetLineWidth(2);
                    pionCutLine->Draw();
                    
                    // Draw text label with right-pointing triangle
                    TLatex *pionLabel = new TLatex();
                    pionLabel->SetTextSize(0.04);
                    pionLabel->SetTextAlign(12); // Left aligned
                    pionLabel->DrawLatex(pionCutValue + 0.01, maxYTotal*maxYScale*0.95, "Unscattered pion subset");
                }
            }
            
            // Only create ratio plot if enabled
            if (showRatioPlots) {
                // Draw the ratio plot in the lower pad
                c->cd(); // Go back to the main canvas before creating a new pad
                TPad *pad2 = new TPad((outName+"pad2").c_str(), "pad2", 0, 0.05, 1, 0.25);
                pad2->SetTopMargin(0);
                pad2->SetBottomMargin(0.3);
                pad2->Draw();
                pad2->cd(); // pad2 becomes the current pad
                gPad->SetLogx(xAxisLog);

                // Create the ratio plot
                TH1F *ratio = (TH1F*)h->Clone("ratio");

                for (int i = 1; i <= ratio->GetNbinsX(); i++)
                {
                    std::cout << "Bin " << i << " orig ratio: " << ratio->GetBinContent(i) <<std::endl;
                }

                for (int i = 1; i <= sum->GetNbinsX(); i++)
                {
                    std::cout << "Bin " << i << " sum: " << sum->GetBinContent(i) <<std::endl;
                }

                ratio->SetLineColor(kBlack);
                ratio->Sumw2();
                ratio->Divide(sum);
                // Get the max and min y values and set the range
                double maxY = 1.2*ratio->GetMaximum();
                double minY = 0.8*ratio->GetMinimum();
                ratio->SetMinimum(minY);
                ratio->SetMaximum(maxY);

                // Adjust the font size of the x-axis labels and title
                ratio->GetXaxis()->SetTitleSize(0.12);
                ratio->GetXaxis()->SetLabelSize(0.12);

                // Adjust the font size of the y-axis labels and title
                ratio->GetYaxis()->SetTitleSize(0.12);
                ratio->GetYaxis()->SetLabelSize(0.12);

                ratio->Draw("E"); // Draw the ratio plot

                // Add a horizontal line at y=1 to the ratio plot (Bottom pad)
                TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 1, ratio->GetXaxis()->GetXmax(), 1);
                line->SetLineColor(kBlack);
                line->SetLineWidth(2);
                line->Draw();

                // No title for the ratio plot
                ratio->SetTitle("");

                // Set a nicely formatted y-axis title for the ratio plot
                ratio->GetYaxis()->SetTitle("Beam On / (MC + EXT)");
                // Set a nicely formatted x-axis title for the ratio plot
                ratio->GetXaxis()->SetTitle(bdtTitle.c_str());
            }

            if (run == "0" && (bdt == "muonBDTScore" || bdt == "protonBDTScore" || bdt == "goldenPionBDTScore") && addUncertaintyBands)
            {
                // Add LaTeX text to the plot
                TLatex latex;
                latex.SetNDC(); // Use normalized device coordinates
                latex.SetTextSize(0.037); // Adjust the font size
                latex.SetTextAlign(22); // Center the text both horizontally and vertically
                latex.DrawLatexNDC(0.5, 0.94, legendStream.str().c_str()); // Adjust the position (NDC coordinates)
                c->Update();
            }

            // Append the ratio plot status to the filename for clarity
            std::string ratioStatus = showRatioPlots ? "_withRatio" : "_noRatio";
            c->SaveAs(("plots/" + outName + "_onlyTesting" + ratioStatus + ".pdf").c_str());
            c->SaveAs(("plots/" + outName + "_onlyTesting" + ratioStatus + ".C").c_str());
        }
    }
    
    std::cout<<"Out of bounds muon BDT events - Data: "<< outOfBoundsMuonBDTScoresForData << " out of " << totalMuonBDTScoresForData << " - as fraction: "<<outOfBoundsMuonBDTScoresForData/totalMuonBDTScoresForData<<std::endl;
    std::cout<<"Out of bounds muon BDT events - MC: "<< outOfBoundsMuonBDTScoresForMC << " out of " << totalMuonBDTScoresForMC << " - as fraction: "<<outOfBoundsMuonBDTScoresForMC/totalMuonBDTScoresForMC<<std::endl;
}