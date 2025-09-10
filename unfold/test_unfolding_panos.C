/*Author: P. Englezos <p.englezos@physics.rutgers.edu>
 *
 * Usage: ./chi_square_cc0pi files.txt
 *
 */

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <set>
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"
#include "EventCategory.hh"
#include "TreeUtils.hh"
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TApplication.h"
#include "TPaveLabel.h"
#include <cassert>
#include <set>
#include <vector>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TPad.h>
#include "TColor.h"
#include "TInterpreter.h"
#include <algorithm>
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "PlotUtils.hh"
#include "TMatrixT.h"
#include "WienerSVDUnfolder.hh"
#include "FiducialVolume.hh"
#include <iomanip>
#include "TMatrixD.h"
#include "MatrixUtils.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"
#include "TLatex.h"
#include "HistUtils.hh"
#include "RooStats/RooStatsUtils.h"

void multiply_1d_hist_by_matrix(TMatrixD *mat, TH1 *hist)
{
   int num_bins = mat->GetNcols();
    TMatrixD hist_mat(num_bins, 1);
    for (int r = 0; r < num_bins; ++r)
    {
        hist_mat(r, 0) = hist->GetBinContent(r + 1);
    }

    TMatrixD hist_mat_transformed(*mat, TMatrixD::EMatrixCreatorsOp2::kMult, hist_mat);

    for (int r = 0; r < num_bins; ++r)
    {
        double val = hist_mat_transformed(r, 0);
        hist->SetBinContent(r + 1, val);
    }
}


void chi_square_all_gens(std::string infiles) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  //std::vector< double > nbins = {0.1,0.275,0.375,0.525,0.7,0.9,1.1,2.};
  std::vector< double > nbins = {-1.,-0.525,-0.375,-0.225,-0.1,0.025,0.15,0.275,0.45,0.6,0.775,0.925,1.};


  //-------Most up-to-date root file for Genie Closure test----------//
  //auto* mcc9 = new MCC9SystematicsCalculator("/uboone/data/users/englezos/jointxsec/ru_github/stv-analysis-new/stv-analysis-new/univmake_genie_closure_mc_total-stats_run1_v2.root", "systcalc.conf" );
    auto* mcc9 = new MCC9SystematicsCalculator("./univmake_genie_angle_run1_only.root", "systcalc.conf" );

  //-------Most up-to-date root file for NuWro fake data test (without AltCV at xsec uncertainties)----------//
  //auto* mcc9 = new MCC9SystematicsCalculator("/uboone/data/users/englezos/jointxsec/ru_github/stv-analysis-new/stv-analysis-new/univmake_nuwro_xsec_mc_total-stats_run1_v2.root", "systcalc.conf" );

   //-------Most up-to-date root file for NuWro fake data test (including AltCV at xsec uncertainties)----------//
   //auto* mcc9 = new MCC9SystematicsCalculator("/uboone/data/users/englezos/jointxsec/ru_github/stv-analysis-new/stv-analysis-new/univmake_nuwro_fake_runs1-3.root","systcalc.conf" );
   //auto* mcc9 = new MCC9SystematicsCalculator("/uboone/data/users/englezos/jointxsec/ru_github/stv-analysis-new/stv-analysis-new/univmake_nuwro_fake_runs1-3_muon_angle.root","systcalc.conf" );

   //------------Genie with all uncertainties------------//
   //auto* mcc9 = new MCC9SystematicsCalculator("/uboone/data/users/englezos/jointxsec/ru_github/stv-analysis-new/stv-analysis-new/univmake_genie_with_uncertainties_runs1-3_muon_angle.root","systcalc.conf" );

  //std::vector< double > nbins = {-1.,-0.525,-0.375,-0.225,-0.1,0.025,0.15,0.275,0.45,0.6,0.775,0.925,1.};
  //auto* mcc9 = new MCC9SystematicsCalculator("./output_files/univmake_corrected-muon-angle_xsec_stat_nuwro_uncertainties.root", "systcalc.conf" );

  TH1D* genie_cv_truth_vals = new TH1D("genie_cv_truth_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());
  //TH1D* unfolded_events_vals = dynamic_cast< TH1D* >(genie_cv_truth_vals->Clone("unfolded_events_vals") );
  //TH1D* fake_data_truth_vals = dynamic_cast< TH1D* >(genie_cv_truth_vals->Clone("fake_data_truth_vals") );
  TH1D* fake_data_truth_vals = new TH1D("fake_data_truth_vals", ";p_{muon}; Scaled Events", nbins.size() - 1, nbins.data());
  TH1D* unfolded_events_vals = new TH1D("unfolded_events_vals", ";cos#theta_{muon}; Scaled Events", nbins.size() - 1, nbins.data());

  std::cout<<"Now examining real data."<<std::endl;

  double total_pot = mcc9->total_bnb_data_pot_;
  double integ_flux = integrated_numu_flux_in_FV( total_pot );
  double num_Ar = num_Ar_targets_in_FV();
  double conv_factor = (num_Ar * integ_flux)/1e38;

  TH1D* genie_cv_truth = mcc9->cv_universe().hist_true_.get();
  int num_true_bins = genie_cv_truth->GetNbinsX();
  TH1D* unfolded_events = dynamic_cast< TH1D* >(genie_cv_truth->Clone("unfolded_events") );
  unfolded_events->Reset();

  const auto& fake_data_univ = mcc9->fake_data_universe();
  TH1D* fake_data_truth = fake_data_univ->hist_true_.get();

  std::unique_ptr< Unfolder > unfolder (new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv ) );

  auto result = unfolder->unfold( *mcc9 );

for ( int t = 0; t < num_true_bins; ++t ) {
      double evts = 0.;
      double error = 0.;
      if ( t < nbins.size() - 1) {
         evts = result.unfolded_signal_->operator()( t, 0 );
         error = std::sqrt( std::max(0., result.cov_matrix_->operator()( t, t )) );

         unfolded_events->SetBinContent( t + 1, evts );
         unfolded_events->SetBinError( t + 1, error );
         unfolded_events_vals->SetBinContent( t + 1, evts/conv_factor/(nbins[t+1] - nbins[t]) );
         unfolded_events_vals->SetBinError( t + 1, error/conv_factor/(nbins[t+1] - nbins[t]) );
      }
 }

  TMatrixD *A_C = result.add_smear_matrix_.get();
  multiply_1d_hist_by_matrix(A_C, genie_cv_truth);
  multiply_1d_hist_by_matrix(A_C, fake_data_truth);

  for ( int t = 0; t < nbins.size() - 1; ++t ) {
      fake_data_truth_vals->SetBinContent( t + 1, fake_data_truth->GetBinContent(t + 1)/conv_factor/(nbins[t+1] - nbins[t]) );
      genie_cv_truth_vals->SetBinContent( t + 1, genie_cv_truth->GetBinContent(t + 1)/conv_factor/(nbins[t+1] - nbins[t]) );
  }

  for ( int t = 0; t < result.cov_matrix_->GetNcols() ; ++t ) {
      for ( int u = 0; u < result.cov_matrix_->GetNrows() ; ++u ) {
          result.cov_matrix_->operator()(u,t) = result.cov_matrix_->operator()(u,t)/conv_factor/conv_factor/(nbins[u+1] -  nbins[u])/(nbins[t+1] -  nbins[t]);
      }
  }

  auto inv_cov_mat = invert_matrix(*result.cov_matrix_, 1e-4 );

  std::vector<std::string> filenames;
  std::ifstream myfile(infiles);
  std::copy(std::istream_iterator<std::string>(myfile),
            std::istream_iterator<std::string>(),
            std::back_inserter(filenames));
  std::cout << "File Count: " << filenames.size() << std::endl;

  TCanvas* c2 = new TCanvas("c2");
  //TLegend *leg=new TLegend(0.65,0.65,0.9,0.85);
  TLegend *leg=new TLegend(0.25,0.55,0.5,0.85);

  std::cout<<"Now examining generator true data"<<std::endl;

  for (int i = 0; i < filenames.size(); i++) {

  std::cout << "File: " << filenames.at(i) << std::endl;
  TFile* file = TFile::Open((filenames.at(i)).c_str());
  TTree *tree = (TTree*)file->Get("FlatTree_VARS");

  TH1D* h_muons_gen = new TH1D("h_muons_gen", "; cos#theta_{#mu}; Number of Events", nbins.size() - 1, nbins.data());

  int nfsp, ninitp, mode, pdg[100], pdg_init[100], pdg_nu;
  char cc;
  float  px[100], py[100], pz[100];
  double fSF;

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("nfsp", 1);
  tree->SetBranchStatus("cc", 1);
  tree->SetBranchStatus("pdg", 1);
  tree->SetBranchStatus("px", 1);
  tree->SetBranchStatus("py", 1);
  tree->SetBranchStatus("pz", 1);
  tree->SetBranchStatus("fScaleFactor", 1);
  tree->SetBranchStatus("PDGnu", 1);
  tree->SetBranchStatus("Mode", 1);
  tree->SetBranchStatus("pdg_init", 1);
  tree->SetBranchStatus("ninitp", 1);

  tree->SetBranchAddress("nfsp", &nfsp);
  tree->SetBranchAddress("cc", &cc);
  tree->SetBranchAddress("pdg", &pdg);
  tree->SetBranchAddress("px", &px);
  tree->SetBranchAddress("py", &py);
  tree->SetBranchAddress("pz", &pz);
  tree->SetBranchAddress("fScaleFactor", &fSF);
  tree->SetBranchAddress("PDGnu", &pdg_nu);
  tree->SetBranchAddress("Mode", &mode);
  tree->SetBranchAddress("pdg_init", &pdg_init);
  tree->SetBranchAddress("ninitp", &ninitp);

 for (int i_tree=0; i_tree<tree->GetEntries(); i_tree++) {//Loop over the entries.
       tree->GetEntry(i_tree);

       if ( !(int)cc || pdg_nu != 14 ) continue;//Cut which keeps only events with numu CC interaction

       bool flag_pi0 = true;
       bool flag_picharged = true;

       for (int j = 0; j < nfsp; j++) {//Loop over the final state particles.
            if (pdg[j] == 111){//If a neutral pion is found, the event is rejected.
               flag_pi0 = false;
               break;
            }
            if (abs(pdg[j]) == 211){//Search for charged pion.
                float picharged_ptot = sqrt(px[j]*px[j] + py[j]*py[j] + pz[j]*pz[j]);
                if (picharged_ptot > 0.07){//If charged pion's total momentum is over 70 MeV/c, the event is rejected.
                    flag_picharged = false;
                    break;
                 }
             }
        }

       if (!flag_pi0 || !flag_picharged) continue;//The event is rejected because one of the above cuts is not satisfied.

       bool flag_additional_mu = false;
       double muon_fill = -1.;
       int muon_count = 0;
       for (int n = 0; n < nfsp; n++) {
           if (pdg[n] == 13){//Search for muon

               if (muon_count > 0){
                  std::cout<<"An additional muon has been found. Event is rejected?"<<std::endl;
                  flag_additional_mu = true;
                  break;
               }

               double imuon_ptot = sqrt(px[n]*px[n] + py[n]*py[n] + pz[n]*pz[n]);
               if (imuon_ptot >= 0.1 && imuon_ptot <= 2.){
                   //muon_fill = imuon_ptot;
                   muon_fill = pz[n]/imuon_ptot;
                   muon_count += 1;
               }
            }
       }
       if (!flag_additional_mu && muon_fill > -1.){
          h_muons_gen->Fill(muon_fill);
       }
    }

    multiply_1d_hist_by_matrix(A_C, h_muons_gen);
    for ( int t = 0; t < nbins.size() - 1; ++t ) {
        h_muons_gen->SetBinContent( t + 1, h_muons_gen->GetBinContent(t + 1)*40.*fSF*1e38/(nbins[t+1] - nbins[t]));
    }

   h_muons_gen->SetStats(0);
   h_muons_gen->SetLineWidth(3);

   if (i==0){
         h_muons_gen->GetYaxis()->SetTitle("d#sigma/dp_{#mu} [10^{-38} cm^{2}/(GeV/c)/Ar]");
        //h_muons_gen->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}} [ 10^{-38} #frac{cm^{2}}{GeV/c Ar} ]");

       //h_muons_gen->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#mu} [10^{-38} cm^{2}/(GeV/c)/Ar]");
       //h_muons_gen->SetAxisRange(1.,5.,"Y");
       h_muons_gen->SetLineColor( kOrange + 2);
       h_muons_gen->Draw();
   }
   else{
       if (i==1) h_muons_gen->SetLineColor(kCyan - 3);
       if (i==2) h_muons_gen->SetLineColor(kRed - 4);
       if (i==3) h_muons_gen->SetLineColor(kGreen + 3);
       h_muons_gen->Draw("same");
    }
       if (i==0) leg->AddEntry(h_muons_gen,"NEUT 5.0.4","l");
       if (i==1) leg->AddEntry(h_muons_gen,"Genie 3.0.6","l");
       if (i==2) leg->AddEntry(h_muons_gen,"Genie 2.12.10","l");
       if (i==3) leg->AddEntry(h_muons_gen,"NuWro 19.02","l");

  double chi_square = 0;
  for (int k=0; k<inv_cov_mat->GetNrows(); k++) {
       for (int j=0; j<inv_cov_mat->GetNcols(); j++) {
           chi_square += (unfolded_events_vals->GetBinContent(k+1)-h_muons_gen->GetBinContent(k+1))*(inv_cov_mat->operator()(k,j))*(unfolded_events_vals->GetBinContent(j+1)-h_muons_gen->GetBinContent(j+1));
       }
   }

  double dof = inv_cov_mat->GetNrows();
  double p_value = TMath::Prob(chi_square, dof);
  //double sigma = TMath::Sqrt( TMath::ChisquareQuantile( 1-p_value, dof ) );
  double sigma = RooStats::PValueToSignificance(p_value);
  std::cout<<"chi_square is: "<<chi_square<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  leg->AddEntry((TObject*)0, TString::Format("#chi^{2} / d.o.f = %g / %g", chi_square, dof), "");
  //leg->AddEntry((TObject*)0, TString::Format("p = %g || #sigma = %g ", p_value, sigma), "");
  leg->AddEntry((TObject*)0, TString::Format("p = %g", p_value), "");
}

  double chi_square_cv = 0;
  for (int k=0; k<inv_cov_mat->GetNrows(); k++) {
       for (int j=0; j<inv_cov_mat->GetNcols(); j++) {
           chi_square_cv += (unfolded_events_vals->GetBinContent(k+1) - genie_cv_truth_vals->GetBinContent(k+1))*(inv_cov_mat->operator()(k,j))*(unfolded_events_vals->GetBinContent(j+1) - genie_cv_truth_vals->GetBinContent(j+1));
       }
   }

  leg->AddEntry(genie_cv_truth_vals,"MicroBooNE Tune","l");
  double dof = inv_cov_mat->GetNrows();
  double p_value = TMath::Prob(chi_square_cv, dof);
  std::cout<<"chi_square is: "<<chi_square_cv<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
  double sigma = RooStats::PValueToSignificance(p_value);
  leg->AddEntry(genie_cv_truth_vals, TString::Format("#chi^{2} / d.o.f = %g / %g", chi_square_cv, dof, p_value), "");
  leg->AddEntry(genie_cv_truth_vals, TString::Format("p = %g",  p_value), "");
  //leg->AddEntry(genie_cv_truth_vals, TString::Format("p = %g || #sigma = %g",  p_value, sigma), "");

  double chi_square_fake = 0;
  for (int k=0; k<inv_cov_mat->GetNrows(); k++) {
       for (int j=0; j<inv_cov_mat->GetNcols(); j++) {
           chi_square_fake += (unfolded_events_vals->GetBinContent(k+1) - fake_data_truth_vals->GetBinContent(k+1))*(inv_cov_mat->operator()(k,j))*(unfolded_events_vals->GetBinContent(j+1) - fake_data_truth_vals->GetBinContent(j+1));
        }
   }

  //leg->AddEntry(fake_data_truth_vals,"Truth","l");
  p_value = TMath::Prob(chi_square_fake, dof);
  sigma = RooStats::PValueToSignificance(p_value);
  std::cout<<"chi_square is: "<<chi_square_fake<<", d.o.f. is: "<<dof<<" and p-value is: "<<p_value<<std::endl;
 // leg->AddEntry(fake_data_truth_vals, TString::Format("#chi^{2} / d.o.f = %g / %g ", chi_square_fake, dof), "");
  //leg->AddEntry(fake_data_truth_vals, TString::Format("p = %g || #sigma = %g", p_value, sigma), "");
  //leg->AddEntry(fake_data_truth_vals, TString::Format("p = %g", p_value), "");

  unfolded_events_vals->SetStats(0);
  unfolded_events_vals->SetLineWidth(3);
  unfolded_events_vals->SetLineColor(kBlack);
  unfolded_events_vals->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}} [ 10^{-38} #frac{cm^{2}}{GeV/c Ar} ]");
  unfolded_events_vals->Draw("e same");
  //unfolded_events_vals->Draw("e");
  fake_data_truth_vals->SetLineColor( kBlue );
  fake_data_truth_vals->SetLineWidth( 3 );
  fake_data_truth_vals->SetLineStyle( 2 );
  //fake_data_truth_vals->Draw( "hist same" );
  genie_cv_truth_vals->SetLineColor( kMagenta - 3 );
  genie_cv_truth_vals->SetLineWidth( 3 );
  genie_cv_truth_vals->SetLineStyle( 2 );
  genie_cv_truth_vals->Draw( "hist same" );
  leg->AddEntry(unfolded_events_vals,"Fake beam-on data","l");
  leg->Draw();
  c2->SaveAs("overlay_chi_sqr_CC0pi.pdf");
}

int main(int argc, char* argv[]) {
   chi_square_all_gens(argv[1]);
   return 0;
}