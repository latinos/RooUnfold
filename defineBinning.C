#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include <iomanip>
#include <iostream>
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TRandom.h"

gSystem->Load("RooUnfold/libRooUnfold.so");

#include "RooUnfold/src/RooUnfoldResponse.h"



//------------------------------------------------------------------------------
// defineBinning
//------------------------------------------------------------------------------


void defineBinning(Int_t   jetChannel = 0 , TString theSample = "WW") {

  TH1::SetDefaultSumw2();

  //TString path = Form("rootfiles/%djet/%s/", jetChannel, flavorChannel.Data());
 
  //gSystem->mkdir(path, kTRUE);
 
 
  
  TString path = Form("_GEN_%djet_full.root",  jetChannel); 
  
  TFile* output = new TFile( theSample+path, "recreate");


  // Defining binning
  //----------------------------------------------------------------------------

  //Double_t pt1bins[6] = {25,50,100,150,200,400};
  
  Double_t RECOpt1bins[11] = {10,20,40,60,80,100,125,150,175,200,210};
  Double_t GENpt1bins[7] = {10,20,50,100,150,200,210};



  Double_t pt1bins[9] = {20,40,60,80,100,125,150,175,200}; 
  Double_t ptllbins[7] = {30,40,60,80,100,150,200};
  Double_t mllbins[9] = {20,40,60,80,100,125,150,175,200};
  Double_t dphibins[13] = {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};



  // Pt, Dilepton, DeltaPhi, Mll

  // GEN level ( phase space)  differential histograms 
  //----------------------------------------------------------------------------

  TH1F* hPtLepton1_GEN  = new TH1F("hPtLepton1_GEN",       "", 8, pt1bins);
  TH1F* hPtLepton1_RECO  = new TH1F("hPtLepton1_RECO",       "", 8, pt1bins);
  TH2F* hPtLepton1_RECO_GEN =  new TH2F("hPtLepton1_RECO_GEN", "", 8, pt1bins, 8, pt1bins);

  TH1F* hDilepton_GEN  = new TH1F("hDilepton_GEN",       "", 6, ptllbins);
  TH1F* hDilepton_RECO  = new TH1F("hDilepton_RECO",       "", 6, ptllbins);
  TH2F* hDilepton_RECO_GEN =  new TH2F("hDilepton_RECO_GEN", "", 6, ptllbins, 6, ptllbins);

  TH1F* hmll_GEN  = new TH1F("hmll_GEN",       "", 8, mllbins);
  TH1F* hmll_RECO  = new TH1F("hmll_RECO",       "",8,  mllbins);
  TH2F* hmll_RECO_GEN =  new TH2F("hmll_RECO_GEN", "", 8, mllbins, 8, mllbins);


  TH1F* hdphi_GEN  = new TH1F("hdphi_GEN",       "", 12, dphibins);
  TH1F* hdphi_RECO  = new TH1F("hdphi_RECO",       "",12,  dphibins);
  TH2F* hdphi_RECO_GEN =  new TH2F("hdphi_RECO_GEN", "", 12, dphibins, 12, dphibins);



  // WW level differential histograms 
  //----------------------------------------------------------------------------

  TH1F* hPtLepton1WWLevel_RECO  = new TH1F("hPtLepton1WWLevel_RECO",       "", 8, pt1bins);
  TH1F* hPtLepton1WWLevel_GEN  = new TH1F("hPtLepton1WWLevel_GEN",       "", 8, pt1bins);
  TH2F* hPtLepton1WWLevel_RECO_GEN =  new TH2F("hPtLepton1WWLevel_RECO_GEN", "", 8, pt1bins, 8, pt1bins);

 
  TH1F* hDileptonWWLevel_RECO  = new TH1F("hDileptonWWLevel_RECO",       "", 6, ptllbins);
  TH1F* hDileptonWWLevel_GEN  = new TH1F("hDileptonWWLevel_GEN",       "", 6, ptllbins);
  TH2F* hDileptonWWLevel_RECO_GEN =  new TH2F("hDileptonWWLevel_RECO_GEN", "", 6, ptllbins, 6, ptllbins);

  TH1F* hmllWWLevel_RECO  = new TH1F("hmllWWLevel_RECO",       "", 8, mllbins);
  TH1F* hmllWWLevel_GEN  = new TH1F("hmllWWLevel_GEN",       "", 8, mllbins);
  TH2F* hmllWWLevel_RECO_GEN =  new TH2F("hmllWWLevel_RECO_GEN", "", 8, mllbins,8, mllbins);

  TH1F* hdphiWWLevel_RECO  = new TH1F("hdphiWWLevel_RECO",       "", 12, dphibins);
  TH1F* hdphiWWLevel_GEN  = new TH1F("hdphiWWLevel_GEN",       "", 12, dphibins);
  TH2F* hdphiWWLevel_RECO_GEN =  new TH2F("hdphiWWLevel_RECO_GEN", "", 12, dphibins,12, dphibins);


  TH1D* hPtLepton1WWLevel_nonselected = new TH1D("hPtLepton1WWLevel_nonselected", "hPtLepton1WWLevel_nonselected",8, pt1bins);

  TH1F* hdeltaR = new TH1F("hdeltaR", "hdeltaR", 50,0,0.1);


  // WW level  define response matrix
  
  RooUnfoldResponse responsePtLepton1(hPtLepton1WWLevel_RECO, hPtLepton1WWLevel_GEN);

  RooUnfoldResponse responseDilepton(hDileptonWWLevel_RECO, hDileptonWWLevel_GEN);

  RooUnfoldResponse responseMll(hmllWWLevel_RECO, hmllWWLevel_GEN );

  RooUnfoldResponse responseDphi(hdphiWWLevel_RECO , hdphiWWLevel_GEN);


  //----------------------------------------------------------------------------
  // Input files
  //----------------------------------------------------------------------------

  TString filesPath;

  filesPath = "/gpfs/csic_projects/tier3data/LatinosSkims/ReducedTrees/DiferentialXSection/";

  TChain* tree = new TChain("latino", "latino");

  //tree->Add(filesPath + "latino_000_WWJets2LMad_OF.root");
  //tree->Add(filesPath + "latinostep3_latinosYieldSkim_MC_WWmg.root");
  tree->Add(filesPath + "latino_000_WWJets_madgraph.root");

  // Declaration of leaf types
  //----------------------------------------------------------------------------

  Float_t baseW;        tree->SetBranchAddress("baseW"       , &baseW);
  Float_t channel;      tree->SetBranchAddress("channel"     , &channel);
  Float_t chmet;        tree->SetBranchAddress("chmet"       , &chmet);
  Float_t dataset;      tree->SetBranchAddress("dataset"     , &dataset);
  Float_t dphill;       tree->SetBranchAddress("dphill"      , &dphill);
  Float_t dphilljet;    tree->SetBranchAddress("dphilljet"   , &dphilljet);
  Float_t dphilljetjet; tree->SetBranchAddress("dphilljetjet", &dphilljetjet);
  Float_t drll;         tree->SetBranchAddress("drll"        , &drll);
  Float_t effW;         tree->SetBranchAddress("effW"        , &effW);
  Float_t jeteta1;      tree->SetBranchAddress("jeteta1"     , &jeteta1);
  Float_t jeteta2;      tree->SetBranchAddress("jeteta2"     , &jeteta2);
  Float_t jeteta3;      tree->SetBranchAddress("jeteta3"     , &jeteta3);
  Float_t jetpt1;       tree->SetBranchAddress("jetpt1"      , &jetpt1);
  Float_t jetpt2;       tree->SetBranchAddress("jetpt2"      , &jetpt2);
  Float_t jetpt3;       tree->SetBranchAddress("jetpt3"      , &jetpt3);
  Float_t jettche1;     tree->SetBranchAddress("jettche1"    , &jettche1);
  Float_t jettche2;     tree->SetBranchAddress("jettche2"    , &jettche2);
  Float_t mctruth;      tree->SetBranchAddress("mctruth"     , &mctruth);
  Float_t mll;          tree->SetBranchAddress("mll"         , &mll);
  Float_t mpmet;        tree->SetBranchAddress("mpmet"       , &mpmet); 
  Float_t mth;          tree->SetBranchAddress("mth"         , &mth);
  Float_t nbjet;        tree->SetBranchAddress("nbjet"       , &nbjet);
  Float_t nbjettche;    tree->SetBranchAddress("nbjettche"   , &nbjettche);
  Float_t nextra;       tree->SetBranchAddress("nextra"      , &nextra);
  Float_t njet;         tree->SetBranchAddress("njet"        , &njet);
  Float_t nvtx;         tree->SetBranchAddress("nvtx"        , &nvtx);
  Float_t pchmet;       tree->SetBranchAddress("pchmet"      , &pchmet);
  Float_t pfmet;        tree->SetBranchAddress("pfmet"       , &pfmet);
  Float_t ppfmet;       tree->SetBranchAddress("ppfmet"      , &ppfmet);
  Float_t isomva1;      tree->SetBranchAddress("isomva1"     , &isomva1);
  Float_t isomva2;      tree->SetBranchAddress("isomva2"     , &isomva2);
  Float_t pt1;          tree->SetBranchAddress("pt1"         , &pt1);
  Float_t pt2;          tree->SetBranchAddress("pt2"         , &pt2);
  Float_t pt3;          tree->SetBranchAddress("pt3"         , &pt3);
  Float_t phi1;         tree->SetBranchAddress("phi1"         , &phi1);
  Float_t dymva1;       tree->SetBranchAddress("dymva1"         , &dymva1);
  Float_t phi2;         tree->SetBranchAddress("phi2"         , &phi2);
  Float_t eta1;         tree->SetBranchAddress("eta1"         , &eta1);
  Float_t eta2;         tree->SetBranchAddress("eta2"         , &eta2);
  Float_t ch1;          tree->SetBranchAddress("ch1"         , &ch1);
  Float_t ch2;          tree->SetBranchAddress("ch2"         , &ch2);
  Float_t ptll;         tree->SetBranchAddress("ptll"        , &ptll);
  Float_t softtche;     tree->SetBranchAddress("softtche"    , &softtche);
  Float_t trigger;      tree->SetBranchAddress("trigger"     , &trigger);
  Float_t triggW;       tree->SetBranchAddress("triggW"      , &triggW);
  Int_t   bveto;        tree->SetBranchAddress("bveto"       , &bveto);
  Int_t   bveto_ip;     tree->SetBranchAddress("bveto_ip"    , &bveto_ip);
  Int_t   bveto_mu;     tree->SetBranchAddress("bveto_mu"    , &bveto_mu);
  Int_t   bveto_nj30;   tree->SetBranchAddress("bveto_nj30"  , &bveto_nj30);
  Int_t   dphiveto;     tree->SetBranchAddress("dphiveto"    , &dphiveto);
  Int_t   sameflav;     tree->SetBranchAddress("sameflav"    , &sameflav);
  Int_t   zveto;        tree->SetBranchAddress("zveto"       , &zveto);
  UInt_t  event;        tree->SetBranchAddress("event"       , &event);
  UInt_t  lumi;         tree->SetBranchAddress("lumi"        , &lumi);
  UInt_t  run;          tree->SetBranchAddress("run"         , &run);
  Float_t puW;          tree->SetBranchAddress("puW"         , &puW);
  



// GEN info... 


//Define Status1 leptons 

  Float_t lepGenpt1, lepGenpt2, lepGenpt3;
  tree->SetBranchAddress("genVV_S1lepton1_pt", &lepGenpt1);
  tree->SetBranchAddress("genVV_S1lepton2_pt", &lepGenpt2);
  tree->SetBranchAddress("genVV_S1lepton3_pt", &lepGenpt3);

  Float_t lepGeneta1, lepGeneta2, lepGeneta3;
  tree->SetBranchAddress("genVV_S1lepton1_eta", &lepGeneta1);
  tree->SetBranchAddress("genVV_S1lepton2_eta", &lepGeneta2);
  tree->SetBranchAddress("genVV_S1lepton3_eta", &lepGeneta3);

  Float_t lepGenphi1, lepGenphi2, lepGenphi3;
  tree->SetBranchAddress("genVV_S1lepton1_phi", &lepGenphi1);
  tree->SetBranchAddress("genVV_S1lepton2_phi", &lepGenphi2);
  tree->SetBranchAddress("genVV_S1lepton3_phi", &lepGenphi3);

  Float_t lepGenM1, lepGenM2, lepGenM3;
  tree->SetBranchAddress("genVV_S1lepton1_oVpid", &lepGenM1); 
  tree->SetBranchAddress("genVV_S1lepton2_oVpid", &lepGenM2); 
  tree->SetBranchAddress("genVV_S1lepton3_oVpid", &lepGenM3); 

  Float_t lepGenimTau1, lepGenimTau2, lepGenimTau3;
  tree->SetBranchAddress("genVV_S1lepton1_imTau", &lepGenimTau1); 
  tree->SetBranchAddress("genVV_S1lepton2_imTau", &lepGenimTau2); 
  tree->SetBranchAddress("genVV_S1lepton3_imTau", &lepGenimTau3); 

  Float_t lepGenpid1, lepGenpid2, lepGenpid3;
  tree->SetBranchAddress("genVV_S1lepton1_pid", &lepGenpid1);
  tree->SetBranchAddress("genVV_S1lepton2_pid", &lepGenpid2);
  tree->SetBranchAddress("genVV_S1lepton3_pid", &lepGenpid3);

  Float_t lepGenS3pid1, lepGenS3pid2, lepGenS3pid3;
  tree->SetBranchAddress("genVV_lepton1_pid", &lepGenS3pid1);
  tree->SetBranchAddress("genVV_lepton2_pid", &lepGenS3pid2);
  tree->SetBranchAddress("genVV_lepton3_pid", &lepGenS3pid3);

  Float_t lepGenS3M1, lepGenS3M2, lepGenS3M3;
  tree->SetBranchAddress("genVV_lepton1_oVpid", &lepGenS3M1); 
  tree->SetBranchAddress("genVV_lepton2_oVpid", &lepGenS3M2); 
  tree->SetBranchAddress("genVV_lepton3_oVpid", &lepGenS3M3); 

  // 

 

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
  

  // Float_t Nentries = 179545; 
  //Float_t Nentries = 187525;
  Float_t Nentries = 5000;


  for (int ievent=0; ievent<tree->GetEntriesFast(); ievent++) {
  //for (int ievent=0; ievent<Nentries; ievent++) {
  //for (int ievent=Nentries; ievent<tree->GetEntriesFast(); ievent++) {
   
    tree->GetEntry(ievent);
  


    Float_t luminosity = 19.468;

    Double_t efficiencyW =  effW * triggW ;
    // Double_t totalW = effW * triggW * baseW * efficiencyW * luminosity;

    Double_t totalW = 1;//efficiencyW;

    Double_t totalWGen = 1;//puW;
    
    Double_t totalWReco = 1;//puW * effW * triggW;



    // The GEN selection begins here
    //--------------------------------------------------------------------------
    
    /// ---> 1) Need status 1 leptons to define the same fiducial region
    /// ---> 2) Count how many GEN leptons we have in each bin, applying the fidual region cuts
    /// ---> 3) Apply also, OF, jetbin and opposite-charged cuts.
    

    bool genEvent = false; 

    
    if ( fabs(lepGenpid1) > 20 ) continue;
    if ( fabs(lepGenpid2) > 20 ) continue;
    
    //Select the pair of leptons coming from the two Ws   

    if (fabs(lepGenM1)!= 24) continue; 
    if (fabs(lepGenM2)!= 24) continue; 


    if (lepGenpt1 <= 20) continue;  
    if (lepGenpt2 <= 20) continue;

    if ( fabs(lepGenpid1) == fabs(lepGenpid2) ) continue;
    
    if ( (fabs(lepGenpid1) == 13 && fabs(lepGeneta1) >= 2.4) || 
	 (fabs(lepGenpid1) == 11 && fabs(lepGeneta1) >= 2.5)) continue;

    if ( (fabs(lepGenpid2) == 13 && fabs(lepGeneta2) >= 2.4) || 
    (fabs(lepGenpid2) == 11 && fabs(lepGeneta2) >= 2.5)) continue;


    Float_t dileptonGenPt;
    Float_t mllGen;
    Float_t dphiGen;

    TLorentzVector leptonGen1p4;
    TLorentzVector leptonGen2p4;
    leptonGen1p4.SetPtEtaPhiM(lepGenpt1, lepGeneta1, lepGenphi1, 0.0);
    leptonGen2p4.SetPtEtaPhiM(lepGenpt2, lepGeneta2, lepGenphi2, 0.0);   

    
    dileptonGenPt = (leptonGen1p4+leptonGen2p4).Pt();
    mllGen = (leptonGen1p4+leptonGen2p4).M();
    dphiGen = fabs(leptonGen1p4.DeltaPhi(leptonGen2p4));
   
    hPtLepton1_GEN->Fill(lepGenpt1, totalWGen);//*baseW*luminosity); // leading pt ---> which pt should I store here? 
    
    hDilepton_GEN->Fill(dileptonGenPt,totalWGen); // ptll 
    
    hmll_GEN->Fill(mllGen,totalWGen); // mll
      
    hdphi_GEN->Fill(dphiGen,totalWGen); // deltaPhi
    


    /// ---> 3) Going to apply the selection analysis cuts on RECO objects
 
  
    // The RECO selection begins here. The RECO leptons are supposed to pass the ID+ISO selection already? 
    //----------------------------------------------------------------------------------------------------

    Int_t dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
    
    Float_t metvar = (njet <= 1) ? mpmet : pfmet;
    
    Float_t jetbin = njet;
    
    Float_t dyMVA = ( !sameflav || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) );

    
  
    //--- define migration matrix without any selection efficiency propagated, but fiducial region as reco one.

   
    TLorentzVector lepton1p4;
    TLorentzVector lepton2p4;
    lepton1p4.SetPtEtaPhiM(pt1, eta1, phi1, 0.0);
    lepton2p4.SetPtEtaPhiM(pt2, eta2, phi2, 0.0);
    
    Float_t dileptonPt = (lepton1p4+lepton2p4).Pt();
    Float_t mll = (lepton1p4+lepton2p4).M();
    Float_t deltaphill = fabs(lepton1p4.DeltaPhi(lepton2p4));

 
    hPtLepton1_RECO_GEN->Fill(pt1, lepGenpt1, totalW);
    hPtLepton1_RECO->Fill(pt1, totalWReco);

    hDilepton_RECO->Fill(dileptonPt, totalW);
    hDilepton_RECO_GEN->Fill(dileptonPt, dileptonGenPt, totalW);

    hmll_RECO->Fill(mll, totalW);
    hmll_RECO_GEN->Fill(mll, mllGen , totalW);

    hdphi_RECO->Fill(dphill, totalW);
    hdphi_RECO_GEN->Fill(dphill,dphiGen , totalW);
 

      Bool_t isMatched = true; 
      Float_t deltaR = 9999.9;
    
    if(  pt1 > 20 && pt2 > 20 && !sameflav && ch1*ch2 < 0 &&
	 trigger == 1                        &&
	 nextra == 0                         && 
	 pfmet > 20                          &&
	 mll > 12                            &&
	 (zveto==1 || !sameflav)             &&
	 (mpmet > 20  && dyMVA)              &&
	 (dphiv || !sameflav)                &&
	 bveto_mu                            &&
	 ptll>30 && (!sameflav || ptll>45)   &&
	 jetbin == jetChannel                && 
	 (bveto_ip==1 &&  nbjettche==0)      
	   ){
       

      // DEFINE MATCHING RECO - GEN
      //--------------------------------------------------------------------------
    
      if ( lepton1p4.DeltaR(leptonGen1p4) >= 0.15 && lepton1p4.DeltaR(leptonGen2p4) >= 0.15 )  isMatched = false;
      // (**) Option 1
      //if (isMatched)  cout << 	ievent << "  " << phi1 << "  " << lepGenphi1 << "  " << lepGenpt3  << "  " << isMatched << endl;
      //if ( lepton1p4.DeltaR(leptonGen1p4) >= 0.15 )  isMatched = false;    // (**) Option 2

      //--------------------------------------------------------------------------
   


      Float_t pt1S = smear(pt1);

      //Float_t pt1S = linearW(pt1, 0.006);
    

      //cout << totalWReco  << endl;

      hPtLepton1WWLevel_RECO->Fill(pt1,    totalWReco); 
      hPtLepton1WWLevel_GEN->Fill(lepGenpt1,totalW);
      hPtLepton1WWLevel_RECO_GEN->Fill(pt1,lepGenpt1, totalW);
 
      hDileptonWWLevel_RECO->Fill(dileptonPt,        totalWReco); 
      hDileptonWWLevel_GEN->Fill(dileptonGenPt,totalW);
      hDileptonWWLevel_RECO_GEN->Fill(dileptonPt,dileptonGenPt , totalW);
     
      hmllWWLevel_RECO->Fill(mll,        totalWReco); 
      hmllWWLevel_GEN->Fill(mllGen,totalW);
      hmllWWLevel_RECO_GEN->Fill(mll , mllGen, totalW);
 
      hdphiWWLevel_RECO->Fill(dphill,        totalWReco); 
      hdphiWWLevel_GEN->Fill(dphiGen,totalW);
      hdphiWWLevel_RECO_GEN->Fill(dphill,dphiGen, totalW);


      //---- Fill response matrix 

      if (isMatched) {
	responsePtLepton1.Fill(pt1, lepGenpt1, totalWReco);
	responseDilepton.Fill(dileptonPt, dileptonGenPt, totalWReco);
	responseMll.Fill(mll, mllGen, totalWReco);
	responseDphi.Fill(dphill, dphiGen, totalWReco);
      }	else {

	responsePtLepton1.Fake(pt1, totalWReco);
	responseDilepton.Fake(dileptonPt, totalWReco);
	responseMll.Fake(mll,  totalWReco);
	responseDphi.Fkae(dphill, totalWReco);

	cout << "Not found matched GEN !!! " << endl;

	//	responsePtLepton1.Miss(lepGenpt1, totalWGen);


      }
	
    }

    else {
	
      hPtLepton1WWLevel_nonselected->Fill(lepGenpt1, totalW);

      //---- Fill response matrix with efficiency 

      responsePtLepton1.Miss(lepGenpt1, totalWGen);
      responseDilepton.Miss(dileptonGenPt, totalWGen);
      responseMll.Miss(mllGen, totalWGen);
      responseDphi.Miss(dphiGen, totalWGen);

      }
    
    
  }


  // Save the histograms
  //----------------------------------------------------------------------------
  output->cd();
  responsePtLepton1.Write();
  responseDilepton.Write();
  responseMll.Write();
  responseDphi.Write();
  output->Write("", TObject::kOverwrite);
  output->Close();



  // Define binning 

  


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






