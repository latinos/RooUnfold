#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include <iomanip>
#include <iostream>
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TRandom.h"


//------------------------------------------------------------------------------
// prepare aMC@NLO for plotting
//------------------------------------------------------------------------------


void prepareAMCatNLO(TString inputFile = "", TString outputFile = "", Int_t   jetChannel = 0, Bool_t jetGenVeto = 0 ) {


  TH1::SetDefaultSumw2();

 
  //----------------------------------------------------------------------------
  // Input files
  //----------------------------------------------------------------------------

  TChain* tree = new TChain("hwwgenanalysis/hwwGenAnalysis");

  tree->Add(inputFile);

  //----------------------------------------------------------------------------
  // Define functions
  //----------------------------------------------------------------------------
  Float_t smear (Float_t xt); 

   //----------------------------------------------------------------------------
  // Output files
  //----------------------------------------------------------------------------
  
  TFile* output = new TFile(outputFile, "recreate");


  // Defining binning
  //----------------------------------------------------------------------------

  //Double_t pt1bins[6] = {25,50,100,150,200,400};
  
  Double_t GENpt1bins[7] = {10,20,50,100,150,200,210};

  const Int_t pt1Nbin = 9;
  const Int_t ptllNbin = 8; 
  const Int_t mllNbin = 9;
  const Int_t dphiNbin = 13;
  const Int_t jetEtNbin = 10;

  Double_t pt1bins[pt1Nbin] = {20,40,60,80,100,125,150,175,200}; 
  Double_t ptllbins[ptllNbin] = {30,40,50,60,70,85,120,150};
  Double_t mllbins[mllNbin] = {20,40,60,80,100,125,150,175,200};
  Double_t dphibins[dphiNbin] = {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};//{0,0.5,1,1.5,2,2.5,3};
  Double_t jetEtbins[jetEtNbin] = {30,40,50,60,70,80,90,100,110,120}; 


  // Pt, Dilepton, DeltaPhi, Mll

  // GEN level ( phase space)  differential histograms 
  //----------------------------------------------------------------------------

  TH1F* hPtLepton1_GEN  = new TH1F("hPtLepton1_GEN",       "", pt1Nbin-1, pt1bins);

  TH1F* hDilepton_GEN  = new TH1F("hDilepton_GEN",       "", ptllNbin-1, ptllbins);

  TH1F* hmll_GEN  = new TH1F("hmll_GEN",       "", mllNbin-1, mllbins);

  TH1F* hdphi_GEN  = new TH1F("hdphi_GEN",       "", dphiNbin-1, dphibins);

  TH1F* hjetEt_GEN  = new TH1F("hjetEt_GEN",       "", jetEtNbin-1, jetEtbins);

  TH1F* hInclusive_GEN  = new TH1F("hInclusive_GEN",       "", 3,0,3);



  // Declaration of leaf types
  //----------------------------------------------------------------------------  

  // Apply NNLL resummation 
  Float_t nllW = 1; //tree->SetBranchAddress("nllW", &nllW);

// GEN info... 


//Define Status1 leptons 

  Float_t lepGenpid1, lepGenpid2;
  tree->SetBranchAddress("id1", &lepGenpid1);
  tree->SetBranchAddress("id2", &lepGenpid2);

  Float_t lepGenpt1, lepGenpt2;
  tree->SetBranchAddress("pt1", &lepGenpt1);
  tree->SetBranchAddress("pt2", &lepGenpt2);

  Float_t lepGeneta1, lepGeneta2;
  tree->SetBranchAddress("eta1", &lepGeneta1);
  tree->SetBranchAddress("eta2", &lepGeneta2);

  Float_t jetGen1_pt, jetGen2_pt;
  tree->SetBranchAddress("jetpt1", &jetGen1_pt);
  tree->SetBranchAddress("jetpt2", &jetGen2_pt);

  Float_t jetGen1_eta, jetGen2_eta;
  tree->SetBranchAddress("jeteta1", &jetGen1_eta);
  tree->SetBranchAddress("jeteta2", &jetGen2_eta);

  Float_t dileptonGenPt;
  Float_t mllGen;
  Float_t dphiGen;
  tree->SetBranchAddress("ptll", &dileptonGenPt);
  tree->SetBranchAddress("mll", &mllGen);
  tree->SetBranchAddress("dphill", &dphiGen);
 

 // Set the channel
  //----------------------------------------------------------------------------
  Float_t SelectedChannel = -999;

  /*  if      (flavorChannel == "MuMu") SelectedChannel =  0;
  else if (flavorChannel == "EE"  ) SelectedChannel =  1;
  else if (flavorChannel == "EMu" ) SelectedChannel =  2;
  else if (flavorChannel == "MuE" ) SelectedChannel =  3;
  else if (flavorChannel == "All" ) SelectedChannel = -1;
  */

  int kk = 0;

 //----------------------------------------------------------------------------
  // Loop
  //----------------------------------------------------------------------------
  

  for (int ievent=0; ievent<tree->GetEntriesFast(); ievent++) {
  //for (int ievent=0; ievent<Nentries; ievent++) {
  //for (int ievent=Nentries; ievent<tree->GetEntriesFast(); ievent++) {
   
    tree->GetEntry(ievent);

    Double_t mybaseW = 5984.0/100400; // aMC@NLO

    Float_t luminosity = 19.365;

    Double_t totalWGen = mybaseW * luminosity * nllW ; //  * puW
  
    // The GEN selection begins here
    //--------------------------------------------------------------------------
    
    /// ---> 1) Need status 1 leptons to define the same fiducial region
    /// ---> 2) Count how many GEN leptons we have in each bin, applying the fidual region cuts
    /// ---> 3) Apply also, OF, jetbin and opposite-charged cuts.
    

    bool genEvent = false; 

    if (lepGenpt1 <= 20) continue;  
    if (lepGenpt2 <= 20) continue;

    if ( fabs(lepGenpid1) == fabs(lepGenpid2) ) continue;
    
    if ( (fabs(lepGenpid1) == 13 && fabs(lepGeneta1) >= 2.4) || 
	 (fabs(lepGenpid1) == 11 && fabs(lepGeneta1) >= 2.5)) continue;

    if ( (fabs(lepGenpid2) == 13 && fabs(lepGeneta2) >= 2.4) || 
	 (fabs(lepGenpid2) == 11 && fabs(lepGeneta2) >= 2.5)) continue;
    


    // If jet veto at GEN level
    //--------------------------------------------------------------------------

    Int_t nGenJets = 0, nGenJet1 = 0, nGenJet2 = 0;
  
    if ( jetGen1_pt>=30 ) nGenJet1++;
    if ( jetGen2_pt>=30 ) nGenJet2++;
   
    nGenJets = nGenJet1 + nGenJet2;
    
    if ( jetGenVeto && nGenJets > 0 )  continue;

    if ( jetChannel && nGenJets != 1 ) continue;
  
    Float_t Genpt1S = smear(lepGenpt1);    
    
   
    hPtLepton1_GEN->Fill(lepGenpt1, totalWGen);//*baseW*luminosity*0.00300652); // leading pt ---> which pt should I store here? 
    
    hDilepton_GEN->Fill(dileptonGenPt,totalWGen); // ptll 
    
    hmll_GEN->Fill(mllGen,totalWGen); // mll
      
    hdphi_GEN->Fill(dphiGen,totalWGen); // deltaPhi

    hjetEt_GEN->Fill(jetGen1_pt, totalWGen); 

    hInclusive_GEN->Fill(1, totalWGen);
 
  

    
 }


  
  // Save the histograms
  //----------------------------------------------------------------------------
  output->cd();
  output->Write("", TObject::kOverwrite);
  output->Close();

  


  
}



//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Float_t smear (Float_t xt)
{

  Float_t cutdummy= -99999.0;

  Float_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Float_t x= gRandom->Rndm();
  //  if (x>xeff) return cutdummy;
  Float_t xsmear= gRandom->Gaus(25,10);     // bias and smear
  return xt+xsmear;
}


//==============================================================================
// Linear Weight
//==============================================================================

Double_t linearW (Double_t xt, Double_t slope)
{

  Double_t weight = 1;

  Double_t inter = 1+ (xt-50)*slope; 

  if ( inter*inter > 0.1) { 

    weight = inter*inter ;
  } else {
    weight = 0.1;
  }
  
     
  return xt*weight;
}






