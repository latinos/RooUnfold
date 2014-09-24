#include <iostream>
using std::cout;
using std::endl;

#include "TSVDUnfold.h"
#include "TRandom.h"
#include "TH1D.h"





gSystem->Load("libRooUnfold");

#include "src/RooUnfold.h"
#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"
#include "src/RooUnfoldTUnfold.h"
#include "src/RooUnfoldErrors.h"



void doUnfoldSignal(Int_t   jetChannel=0, Int_t plot=2, Int_t kreg =5)
{



  gInterpreter->LoadMacro("../draw/draw.C+");
  gInterpreter->LoadMacro("../draw/ChargeRatioStyle.C");

  gInterpreter->LoadMacro("RunToys.C+");


  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);


//----------------------------------------------------------------------------
// Selected options 
//----------------------------------------------------------------------------
 
 bool findRegularization = 0;
 bool doScale = false;
 bool doToyMCPulls = 0;
 bool printErrors = 0;
 bool doResidualStudy = 0;

 bool useOvf = 1;

 //----------------------------------------------------------------------------
 // RooUnfold::SetVerbose(level): 0=warnings, 1=verbose (default, as before), 2=debug, 3=detailed
 //----------------------------------------------------------------------------
  
 Int_t verbose = 0;
 
 //----------------------------------------------------------------------------
 // Input files
 //----------------------------------------------------------------------------
 
  
  TFile* inputWW_GEN     = new TFile("files/WW_GEN_0jet_pow_full_NNLL_JetGenVeto.root");
 
  TFile* inputWW_RECO     = new TFile("files/WW_GEN_0jet_pow_full_NNLL_JetGenVeto.root");

  TFile* inputWW_DATA     = new TFile("WW_0jet_Data.root");

  TString response = "";
  TString mcTruth = ""; 
  TString mcReco = ""; 
  TString dataReco = "";

  Int_t nbin = 0;

  TString title ="";

  if (plot == 0)  {
    
    response = "hPtLepton1WWLevel_RECO_hPtLepton1WWLevel_GEN";
    mcTruth = "hPtLepton1_GEN";
    mcReco = "hPtLepton1WWLevel_RECO";
    title = "Leading Lepton Pt (GeV)";
    dataReco = "hPtLepton1WWLevel_Diff";
  }


 if (plot == 1)  {
    
   response = "hDileptonWWLevel_RECO_hDileptonWWLevel_GEN";
    mcTruth = "hDilepton_GEN";
    mcReco = "hDileptonWWLevel_RECO";
    title = "Dilepton Pt (GeV)";
    dataReco = "hDileptonWWLevel_Diff";
  }


 if (plot == 2)  {
    
   response = "hmllWWLevel_RECO_hmllWWLevel_GEN";
    mcTruth = "hmll_GEN";
    mcReco = "hmllWWLevel_RECO";
    title = "mll (GeV)";
    dataReco = "hmllWWLevel_Diff";
  }


 if (plot == 3)  {
    
   response = "hdphiWWLevel_RECO_hdphiWWLevel_GEN";
    mcTruth = "hdphi_GEN";
    mcReco = "hdphiWWLevel_RECO";
    title = "dphi (rad)";
    dataReco = "hdphiWWLevel_Diff";
  }

 if (plot == 4)  {

    response = "hjetEt_RECO_hjetEt_GEN";
    mcTruth = "hjetEt_GEN";
    mcReco = "hjetEt_RECO";
    title = "Leading Jet Et (GeV)";
    dataReco = "hPtLepton1WWLevel_Diff";
  }

if (plot == 5)  {

    response = "hNjet_RECO_hNjet_GEN";
    mcTruth = "hNjet_GEN";
    mcReco = "hNjet_RECO";
    title = "# jets (Pt>30 GeV)";
    dataReco = "hPtLepton1WWLevel_Diff";
  }


  TH1F *h_mcTruth;
  h_mcTruth = (TH1F*) inputWW_GEN->Get(mcTruth);
  h_mcTruth->SetName("h_mcTruth");

  TH1F *h_mcReco;
  h_mcReco = (TH1F*) inputWW_RECO->Get(mcReco);
  h_mcReco->SetName("h_mcReco");

  if ( h_mcTruth->GetNbinsX() == h_mcReco->GetNbinsX() )  nbin = h_mcTruth->GetNbinsX();

  cout << nbin << endl;

  // ----> Substract backgrounds from data histrogram.
  //==============================================================================
  //  hData_bkgSub= (TH1D*) hData->Clone();
  //hData_bkgSub->Add(hBkg,-1.0);

 
  if ( doScale) {

    cout << "==================================== TEST ====================================" << endl;
    
    // ---->  Scale up the measured distribution, 30% up

    for (int ib =1; ib < nbin; ib++ ){
      
      Float_t entry = h_mcReco->GetBinContent(ib);
      Float_t entryE = h_mcReco->GetBinError(ib);
      
      h_mcReco->SetBinContent(ib, entry);
      h_mcReco->SetBinError(ib, entry*0.3);

	}
  }
    


 
  cout << "==================================== TRAIN ====================================" << endl;

  const RooUnfoldResponse *responseFinal =  inputWW_GEN->Get(response);

 
  /// To be included also the underflow bin
  if (useOvf ) { 
    responseFinal->UseOverflow();    
  }

  

  //==============================================================================
  // Unfold
  //==============================================================================

  if (verbose>=0) cout << "Create RooUnfold object for method " << endl;


  // ---> initialize SVD unfolding
 
  RooUnfoldSvd   unfold(responseFinal, h_mcReco,1);   
  unfold.SetRegParm(kreg);
 

  //RooUnfoldBayes unfold(responseFinal, h_mcReco,5); 
  
  //RooUnfoldTUnfold  unfold(responseFinal, h_mcReco); 
  // RooUnfoldBinByBin unfold(responseFinal, h_mcReco);

  //RooUnfoldInvert unfold(responseFinal, h_mcReco);

  unfold.SetVerbose (verbose);

 
  
 
  // ----> Choose error treatment::  
  //       0: Errors are the square root of the bin content
  //       1: Errors from the diagonals of the covariance matrix given by the unfolding 
  //          (variance values) 
  //       2: Errors from the covariance matrix given by the unfolding
  //       3: Errors from the covariance matrix from the variation of the results in toy MC tests
  //==============================================================================
 

  //unfold.IncludeSystematics(2); // Default 1, propagates both statistical+systematics, =2 no errors included.
 
  
  unfold.SetNToys(1000);

  RooUnfold::ErrorTreatment *doerror = 1;//RooUnfold::kCovariance;

  if (verbose>=0) cout << "Error treatment" << doerror << endl;


  if ( printErrors) {
    TH1* t= 0;
    RooUnfoldErrors *unfoldErr = new RooUnfoldErrors(1000,unfold,t);
    TH1 *hToyErr=unfoldErr->RMSResiduals();
    TH1 *hUnfErr=unfoldErr->UnfoldingError();
    hUnfErr->Draw("hist");
    hToyErr->Draw("P SAME");
  }

  if (verbose>=0)
    unfold.PrintTable (cout, h_mcTruth, (RooUnfold::ErrorTreatment)doerror);


  // ---->  choose regularization parameter as in top paper 
  //==============================================================================
 
  TSVDUnfold *myTSVD = (TSVDUnfold*) unfold->Impl();

  TH2D *covBM = (TH2D*) myTSVD->GetBCov();
 //TH2D *covM = (TH2D*) myTSVD->GetUnfoldCovMatrix(covBM, 1000);
  TH2D *covM = (TH2D*) myTSVD->GetXinv();

  //covM->DrawCopy("COLZTEXT");


  if (findRegularization) {

    
    RooUnfoldParms *parms= new RooUnfoldParms(&unfold,(RooUnfold::ErrorTreatment)doerror, h_mcTruth);
  
    parms->SetMinParm(1);
    parms->SetMaxParm(nbin);
    parms->SetStepSizeParm(1);

    TProfile *hParmChi2= parms->GetChi2();
    TProfile *hParmErr= parms->GetRMSError();  // Mean values of errors in each bin. 
    TProfile *hParmRes= parms->GetMeanResiduals();
    TH1      *hParmRms= parms->GetRMSResiduals();


    // -----> Account for mean correlation
  
    TH1D *hParmCorr = new TH1D("","", nbin, 1, nbin+1);

    Int_t min = 1;
    Int_t max = nbin;
    Int_t step = 1;
  
    for ( int i =  min; i <= max; i=i+step) {

      RooUnfoldSvd   unfoldTest(responseFinal, h_mcReco,i);   

      TMatrixD m_covMatTest (nbin,nbin);       
      m_covMatTest = unfoldTest.Ereco((RooUnfold::ErrorTreatment)doerror); 
 
  
      TH2D *hCorrTest = CorrelationHist ( m_covMatTest,    
					  "corr", "Unfolded correlation matrix",
					  responseFinal->Hresponse()->GetYaxis()->GetXmin(),
					  responseFinal->Hresponse()->GetYaxis()->GetXmax(), false);    
  

    
      Double_t sum = 0; 
      
      for ( int nbi=1; nbi <= nbin; nbi++) { 
   
      Double_t getBinMax = -999; 
    
      for ( int nbj=1; nbj <= nbin; nbj++) { 

	if (nbi==nbj)  continue; 
 
	Double_t getBin = fabs(hCorrTest->GetBinContent(nbi,nbj));
	
	if ( getBin > getBinMax )  getBinMax = getBin; 

       }
      
      sum += getBinMax*getBinMax;     
     
    }
   
      hParmCorr->SetBinContent(i,sqrt(sum)/(nbin));  
      if(i==7) cout << sqrt(sum)/(nbin) << endl;
    }
      




  hParmCorr->SetMarkerStyle(4);
  hParmCorr->SetTitle("Global correlation");
  
  TCanvas * param = new TCanvas("param", "param", 550, 550);
 
  /*  
  param->Divide(3,2);
  param->cd(1);
  hParmChi2->Draw("P");    
  param->cd(2);
  hParmErr->Draw("P"); 
  param->cd(3);
  hParmRes->Draw("P"); 
  param->cd(4);
  hParmRms->Draw("P"); 
  param->cd(5);
  //  hParmCorr->Draw("P"); 
  param->cd(6);
 */ 
  
  // ---->   Returns d vector (for choosing appropriate regularisation)
  //==============================================================================
  
  if (useOvf) nbin = nbin+2;

  TSVDUnfold *myTSVD = (TSVDUnfold*) unfold->Impl();

  TH1D *svVector = myTSVD->GetSV();
  
  TH1D *dVector = myTSVD->GetD();
  TH1D *diVector = (TH1D*) dVector->Clone();
  TH1D *dzVector = (TH1D*) dVector->Clone();

  double tau = svVector->GetBinContent(kreg+1);

  for (int i=1; i < nbin; i++) {

    double Si = svVector->GetBinContent(i);  
    double scale = (Si*Si)/((Si*Si)+(tau*tau));
    double di = dVector->GetBinContent(i);
  
    dzVector->SetBinContent(i, di*scale);
    diVector->SetBinContent(i, di);
  }

  
  dVector->GetYaxis()->SetTitle("d_{i}");
  dVector->GetXaxis()->SetTitle("i");
  
  if (plot == 0) dVector->SetTitle("SVD unfolding d_{i} for p_{T}^{max}");

  dVector->SetMarkerStyle(4);
  dVector->DrawCopy("hist");  
  
  dzVector->SetLineStyle(2);
  dzVector->SetLineColor(kRed);
  dzVector->Draw("histSAME");  
  
 
   TLine().DrawLine(dVector->GetBinLowEdge(1), 1.0, dVector->GetBinLowEdge(nbin+1), 1.0);  // draw a line at y=0;
  
  
  
  }
  

  // ---->  get unfolded histogram 
  //==============================================================================
  TH1F* h_mcReco_unfolded = (TH1F*) unfold.Hreco((RooUnfold::ErrorTreatment)doerror); 
  h_mcReco_unfolded->SetName("h_mcReco_unfolded");


  // ---->  create matrix to store covariance matrix
  //==============================================================================

  if (useOvf) nbin = nbin+2;

  TMatrixD m_covMat (nbin,nbin);
  m_covMat = unfold.Ereco((RooUnfold::ErrorTreatment)doerror); 
  //m_covMat->Draw("COLTEXT");
 
  
  TMatrixD m_errMat(nbin,nbin);

  for (Int_t i=0; i<nbin; i++) { 
    for (Int_t j=0; j<nbin; j++) { 
      
      m_errMat(i,j)= m_covMat(i,j)>=0 ? sqrt(m_covMat(i,j)) : -sqrt(-m_covMat(i,j));
    }
  }

  if (verbose>=0) PrintMatrix(m_errMat,"","covariance matrix",10);


  TH2D *hCorr = CorrelationHist ( m_covMat,    
				  "corr", "Unfolded correlation matrix",
				  responseFinal->Hresponse()->GetYaxis()->GetXmin(),
				  responseFinal->Hresponse()->GetYaxis()->GetXmax(), useOvf);




  // ---->  Calculate pulls and residuals
  //==============================================================================
  TH1F *hRes = (TH1F*) h_mcReco_unfolded->Clone("res");
  hRes  ->Reset();
  hRes  ->SetTitle ("Residuals");

  TH1F *hPulls= (TH1F*) h_mcReco_unfolded->Clone("pulls");
  hPulls->Reset();
  hPulls->SetTitle ("Pulls");


  for (Int_t i= 1; i<=nbin; i++) {
    if ((h_mcReco_unfolded->GetBinContent(i)!=0.0 || (doerror && h_mcReco_unfolded->GetBinError(i)>0.0)) &&
        (h_mcTruth->GetBinContent(i)!=0.0 || (doerror && h_mcTruth->GetBinError(i)>0.0))) {
      Double_t res= h_mcReco_unfolded->GetBinContent(i) - h_mcTruth->GetBinContent(i);
      Double_t err= h_mcReco_unfolded->GetBinError  (i);
      hRes->SetBinContent (i, res);
      hRes->SetBinError   (i, err);
      if (err>0.0) {
        hPulls->SetBinContent (i, res/err);
        hPulls->SetBinError   (i, 1.0);
      }
    }
  }



  // ---->  Get response matrix
  //==============================================================================
  TMatrixD outResponse (nbin,nbin);
  TMatrixD outResponse = responseFinal->Mresponse();
  //outResponse->Draw("COLZTEXT");




  TVectorD outMeasured = responseFinal->Emeasured();
  //outMeasured->Draw("TEXT");

  TVectorD outThruth = responseFinal->Etruth();
  
  TH1* t= 0;
  //RooUnfoldErrors errors= new RooUnfoldErrors(100,unfold,t);
  

 
  // ----> Plot results
  //==============================================================================
 
  
  TCanvas * unfolding = new TCanvas("unfolding", "unfolding", 1400, 600);
  unfolding->Divide(3,1);
  
  unfolding->cd(1);

 
  //h_mcTruth->SetTitle ("Unfold");
  h_mcTruth->GetXaxis()->SetTitle(title);

  h_mcTruth->Draw("histE1");
  h_mcTruth->SetLineColor(kRed);
  h_mcTruth->SetLineWidth(2);
  h_mcTruth->SetMarkerStyle(20);
  h_mcTruth->SetMarkerSize(1.0);
  h_mcTruth->SetMarkerColor(kRed);

  h_mcReco_unfolded->Draw("sameE1");
  h_mcReco_unfolded->SetMarkerStyle(20);
  h_mcReco_unfolded->SetMarkerColor(kBlue);
  h_mcReco_unfolded->SetMarkerSize(1.0);
  h_mcReco_unfolded->SetLineWidth(2);


  h_mcReco->Draw("SAMEhistE1");
  h_mcReco->SetLineColor(8);
  h_mcReco->SetLineWidth(2);
  h_mcReco->SetMarkerStyle(20);
  h_mcReco->SetMarkerSize(1.0);
  h_mcReco->SetMarkerColor(8);

  unfolding->cd(2);

  hRes->SetMarkerStyle(kFullDotLarge);
  hRes->Draw();
  TLine().DrawLine(hRes->GetBinLowEdge(1), 0.0, hRes->GetBinLowEdge(nbin+1), 0.0);  // draw a line at y=0;
  
  cout << hRes->GetMean() << endl;

  unfolding->cd(3);

  hPulls->SetMarkerStyle(kFullDotLarge);
  hPulls->Draw("P");
  TLine().DrawLine(hPulls->GetBinLowEdge(1), 0.0, hPulls->GetBinLowEdge(nbin+1), 0.0);  // draw a line at pull=0;


cout << hPulls->GetMean() << endl;

  TCanvas * unfoldingAN = new TCanvas("unfoldingAN", "unfoldingAN", 550, 1.2*600);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  pad1->SetTopMargin   (0.05);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.31); 
  pad2->SetTopMargin   (0.08);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();

  pad1->cd();

  //h_mcTruth->SetTitle ("Unfold");
  h_mcTruth->GetXaxis()->SetTitle(title);

  h_mcTruth->Draw("histE1");
  h_mcReco_unfolded->Draw("sameE1");
  h_mcReco->Draw("SAMEhistE1");


  if (plot!=2 && plot!=3) {
    DrawLegend(0.53, 0.74 , h_mcTruth, " mcTruth",    "lp", 0.035, 0.2, 0.05);
    DrawLegend(0.53, 0.64 , h_mcReco_unfolded, " mcUnfolded",    "lp", 0.035, 0.2, 0.05);
      DrawLegend(0.53, 0.54 , h_mcReco, " mcReco",    "lp", 0.035, 0.2, 0.05);
  }
  
  if (plot==2) {
    DrawLegend(0.63, 0.84 , h_mcTruth, " mcTruth",    "lp", 0.035, 0.2, 0.05);
    DrawLegend(0.63, 0.74 , h_mcReco_unfolded, " mcUnfolded",    "lp", 0.035, 0.2, 0.05);
    DrawLegend(0.63, 0.64 , h_mcReco, " mcReco",    "lp", 0.035, 0.2, 0.05);
  }
  if (plot==3) {
    DrawLegend(0.23, 0.74 , h_mcTruth, " mcTruth",    "lp", 0.035, 0.2, 0.05);
    DrawLegend(0.23, 0.64 , h_mcReco_unfolded, " mcUnfolded",    "lp", 0.035, 0.2, 0.05);
    DrawLegend(0.23, 0.54 , h_mcReco, " mcReco",    "lp", 0.035, 0.2, 0.05);
  }

  pad2->cd();
  
  TH1F* ratio       = h_mcReco_unfolded->Clone("ratio");
  TH1F* uncertainty = h_mcTruth->Clone("uncertainty");

  for (UInt_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) {

    Double_t mcValue = h_mcTruth->GetBinContent(ibin);
    Double_t mcError = h_mcTruth->GetBinError  (ibin);
    
    Double_t dtValue = ratio->GetBinContent(ibin);
    Double_t dtError = ratio->GetBinError  (ibin);
    
    Double_t ratioValue       = (mcValue > 0) ? dtValue/mcValue : 0.0;
    Double_t ratioError       = (mcValue > 0) ? dtError/mcValue : 0.0;
    Double_t uncertaintyError = (mcValue > 0) ? mcError/mcValue : 0.0;
    
    ratio->SetBinContent(ibin, ratioValue);
    ratio->SetBinError  (ibin, ratioError);
    
    uncertainty->SetBinContent(ibin, 1.0);
    uncertainty->SetBinError  (ibin, uncertaintyError);
  }

    
    uncertainty->Draw("e2");
    ratio      ->Draw("ep,same");
    ratio      ->SetLineColor(kBlue);
    ratio      ->SetMarkerSize(1.0);
    ratio      ->SetLineWidth(2);
    ratio      ->SetMarkerStyle(20);
    uncertainty->SetLineColor(kBlue);
    uncertainty->SetFillColor  (kGray+2);
    uncertainty->SetFillStyle  (   3345);
    uncertainty->SetLineColor  (kGray+2);
    uncertainty->SetMarkerColor(kGray+2);
    uncertainty->SetMarkerSize (      0);

    uncertainty->GetYaxis()->SetRangeUser(0.4, 1.6);

    Pad2TAxis(uncertainty, h_mcTruth->GetXaxis()->GetTitle(), "unfolded / generated");




  TCanvas * corrMatrix = new TCanvas("corrMatrix", "corrMatrix", 600, 600);
  corrMatrix->cd();

  //  hCorr->GetXaxis()->SetTitle("Reco");
  //  hCorr->GetYaxis()->SetTitle("Truth");
  hCorr->GetYaxis()->SetTitleOffset(1.2);
  hCorr->Draw("COLZTEXT");
  

  //TCanvas * responseMatrix = new TCanvas("responseMatrix", "responseMatrix", 600, 600);
  //responseMatrix->cd();
  // outResponse->Draw("COLZ");


  if (doResidualStudy ) { 

    TCanvas * residual = new TCanvas("residual", "residual", 550, 550);
    residual->cd();

    // ---->  Calculate pulls and residuals
    //==============================================================================
  
    TH1F *hRes = (TH1F*) h_mcReco_unfolded->Clone("res");
    hRes  ->Reset();
    hRes  ->SetTitle ("Residuals");

    TH1F *hRes_up1s   = (TH1F*) h_mcTruth->Clone("res");
    TH1F *hRes_up2s   = (TH1F*) h_mcTruth->Clone("res");
    TH1F *hRes_down1s = (TH1F*) h_mcTruth->Clone("res");
    TH1F *hRes_down2s = (TH1F*) h_mcTruth->Clone("res");

    TH1F *hPulls= (TH1F*) h_mcReco_unfolded->Clone("pulls");
    hPulls->Reset();
    hPulls->SetTitle ("Pulls");

  

    for (Int_t i= 1; i<=nbin; i++) {
      if ((h_mcReco_unfolded->GetBinContent(i)!=0.0 || (doerror && h_mcReco_unfolded->GetBinError(i)>0.0)) &&
	  ( h_mcTruth->GetBinContent(i)!=0.0 || (doerror &&  h_mcTruth->GetBinError(i)>0.0))) {
      
	Double_t res= h_mcReco_unfolded->GetBinContent(i) -  h_mcTruth->GetBinContent(i);
	Double_t err= h_mcReco_unfolded->GetBinError  (i);
      
	Double_t Err = h_mcTruth->GetBinError  (i);
	Double_t Bias = Err; 


	hRes->SetBinContent (i, res);
	hRes->SetBinError   (i, err);

	hRes_up1s->SetBinContent(i, Bias);
      	hRes_up2s->SetBinContent(i, 2*Bias);
	hRes_down1s->SetBinContent(i, -Bias);
      	hRes_down2s->SetBinContent(i, -2*Bias);


	if (err>0.0) {
	  hPulls->SetBinContent (i, res/err);
	  hPulls->SetBinError   (i, 1.0);
	}
      }
    }
    
     hRes->Draw();
    
    hRes_up2s->Draw("histsame");
    hRes_down2s->Draw("histsame");
    hRes_up2s->SetFillColor(kYellow);
    hRes_up2s->SetFillStyle(3004);
    hRes_up2s->SetLineColor(kYellow);
    hRes_down2s->SetFillColor(kYellow);
    hRes_down2s->SetFillStyle(3004);
    hRes_down2s->SetLineColor(kYellow);


    hRes_up1s->Draw("histsame");
    hRes_down1s->Draw("histsame");
    hRes_up1s->SetFillColor(kGreen);
    hRes_up1s->SetFillStyle(3002);
    hRes_up1s->SetLineColor(kGreen);
    hRes_down1s->SetFillColor(kGreen);
    hRes_down1s->SetFillStyle(3002);
    hRes_down1s->SetFillStyle(3002);
    hRes_down1s->SetLineColor(kGreen);
  
    hRes->SetMarkerStyle(kFullDotLarge);
    hRes->Draw("same");

    TLine().DrawLine(hRes->GetBinLowEdge(1), 0.0, hRes->GetBinLowEdge(nbin+1), 0.0);  // draw a line at y=0;
    
 }


  // ----> MC closure test:: Generate pulls in each bin for N toy MC 
  //==============================================================================

  if (doToyMCPulls) { 

    if (useOvf) nbin = nbin-2;

    Int_t NToys = 50; 
 
    const Int_t nB = nbin; 

    TH1F *hpull[nB]; 

    for (int n = 0; n < nbin; n++)  {
      delete gDirectory->FindObject(Form("Bin%.1i ",  n+1));
      hpull[n] = new TH1F(Form("Bin%.1i ",  n+1), Form("Bin%.1i ",  n+1), 50, -0.06, 0.06);
      hpull[n]->SetStats(1);
      hpull[n]->GetXaxis()->SetTitle("X_{unfold}-X_{truth} / #sigma");
      hpull[n]->SetMarkerStyle(kFullCircle);
    }

    

    for (Int_t k=1; k<NToys; k++){
  
      delete gDirectory->FindObject("newUnfold");
      delete gDirectory->FindObject("newUnfoldRes");

      RooUnfoldBayes* newUnfold = new RooUnfoldBayes("newUnfold", "");
      const RooUnfoldResponse* newUnfoldRes = new RooUnfoldResponse("newUnfoldRes", "");
      
      newUnfold->SetMeasured(h_mcReco); 
      
      newUnfoldRes->Setup( RunToys() );
      
      newUnfold->SetResponse(newUnfoldRes);
      
      newUnfold->SetVerbose (0);
     
      
      // ---->  get unfolded histogram 
      //==============================================================================
      delete gDirectory->FindObject("");
     			    
      TH1F* h_unfolded = (TH1F*) newUnfold.Hreco((RooUnfold::ErrorTreatment)doerror); 
      
      //cout << h_unfolded->GetBinContent(2) << "  " << h_mcTruth->GetBinContent(2) << endl;
     
      for (Int_t i= 1; i<nbin+1; i++) {
	if ((h_unfolded->GetBinContent(i)!=0.0 || (doerror && h_unfolded->GetBinError(i)>0.0)) &&
        (h_mcTruth->GetBinContent(i)!=0.0 || (doerror && h_mcTruth->GetBinError(i)>0.0))) {
	  Double_t res= h_unfolded->GetBinContent(i) - h_mcTruth->GetBinContent(i);
	  Double_t err= h_unfolded->GetBinError  (i);
	    
	     if (err>0.0) {
	    hpull[i-1]->Fill(res/err);
	   
	  }
	}
      }
      
    }

    gStyle->SetOptStat(1111);

    Int_t nPads = nB/2; 

    if( (nPads*2 - nbin)  <  0 ) nPads = nPads+1;
    //hpull[0]->Draw("P");
    
    TCanvas * pullHistos = new TCanvas("pullHistos", "pullHistos", 1200, 600);
    pullHistos->Divide(nPads, 2);
  
    for (int p=0; p < nbin; p++ ) {
      pullHistos->cd(p+1);
      hpull[p]->Draw("P");
    }
    
   
  }
}

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}



/*
//==============================================================================
// Add bin correlations to measured errors
//==============================================================================

void SetMeasuredCov ()
{
  if (bincorr==0.0) return;
  TMatrixD cov= unfold->GetMeasuredCov();  // initially diagonal
  Double_t corr= bincorr;
  for (Int_t k=1; k<nmbins; k++) {
    for (Int_t i=k; i<nmbins; i++)
      cov(i,i-k)= cov(i-k,i)= corr*cov(i,i);
    corr *= bincorr;
    if (corr==0.0) break;
  }
  unfold->SetMeasuredCov (cov);
  hMeasCorr= CorrelationHist (cov, "measCor", "Measured correlation matrix",
                              response->Hresponse()->GetXaxis()->GetXmin(),
                              response->Hresponse()->GetXaxis()->GetXmax());
}
*/



//==============================================================================
//Tool for printing covariance matrix 
//==============================================================================

void PrintMatrix(const TMatrixD& m, const char* format,
                                       const char* name, Int_t cols_per_sheet)
{

  // Print the matrix as a table of elements.
   // Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
   // By default the format "%11.4g" is used to print one element.
   // One can specify an alternative format with eg
   //  format ="%6.2f  "

   if (!m.IsValid()) {
     m.Error("PrintMatrix","%s is invalid",name);
     return;
   }

   const Int_t ncols  = m.GetNcols();
   const Int_t nrows  = m.GetNrows();
   const Int_t collwb = m.GetColLwb();
   const Int_t rowlwb = m.GetRowLwb();

   if (!(format && format[0])) format= "%11.4g ";
   char topbar[1000];
   snprintf(topbar,1000,format,123.456789);
   Int_t nch = strlen(topbar)+1;
   if (nch > 18) nch = 18;
   char ftopbar[20];
   for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
   Int_t nk = 1 + Int_t(log10(ncols));
   snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
   Int_t nch2 = strlen(ftopbar);
   for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
   ftopbar[nch] = '|';
   ftopbar[nch+1] = 0;

   printf("\n%dx%d %s is as follows",nrows,ncols,name);

   if (cols_per_sheet <= 0) {
     cols_per_sheet = 5;
     if (nch <= 8) cols_per_sheet =10;
   }
   nk = 5+nch*(cols_per_sheet<ncols ? cols_per_sheet : ncols);
   for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
   topbar[nk] = 0;
   for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
      printf("\n\n     |");
      for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
         printf(ftopbar,j+collwb-1);
      printf("\n%s\n",topbar);
      if (m.GetNoElements() <= 0) continue;
      for (Int_t i = 1; i <= nrows; i++) {
         printf("%4d |",i+rowlwb-1);
         for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
            printf(format,m(i+rowlwb-1,j+collwb-1));
         printf("\n");
      }
   }
   printf("\n");

}


TH2D* CorrelationHist (const TMatrixD& cov,
		       const char* name, const char* title,
		       Double_t lo, Double_t hi, Bool_t _useOvf)
{

  // correlation coef == cov(xy)/(sigmay*sigmay);

  Int_t nb= cov.GetNrows();

  TH2D* h; 

  h = new TH2D (name, title, nb, lo, hi, nb, lo, hi);

  if (_useOvf ) h = new TH2D (name, title, nb-2, lo, hi, nb-2, lo, hi);
  
  h->SetAxisRange (-1.0, 1.0, "Z");
  
  if (_useOvf) nb = nb-2;

  for(int i=0; i < nb; i++) {
  
    for(int j=0; j < nb; j++) {
      
      if (_useOvf ) {
      
	Double_t Viijj= cov(i+1,i+1)*cov(j+1,j+1);
	h->SetBinContent (i+1, j+1, cov(i+1,j+1)/sqrt(Viijj));
      } else {
	Double_t Viijj= cov(i,i)*cov(j,j);
	h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
      }
    }
  }

  return h;
}



//------------------------------------------------------------------------------
// Pad2TAxis
//------------------------------------------------------------------------------
void Pad2TAxis(TH1* hist, TString xtitle, TString ytitle)
{
  TAxis* xaxis = (TAxis*)hist->GetXaxis();
  TAxis* yaxis = (TAxis*)hist->GetYaxis();

  xaxis->SetLabelFont  (    42);
  xaxis->SetLabelOffset( 0.025);
  xaxis->SetLabelSize  (   0.1);
  xaxis->SetNdivisions (   505);
  xaxis->SetTitle      (xtitle);
  xaxis->SetTitleFont  (    42);
  xaxis->SetTitleOffset(  1.35);
  xaxis->SetTitleSize  (  0.11);

  yaxis->CenterTitle   (      );
  yaxis->SetLabelFont  (    42);
  yaxis->SetLabelOffset(  0.02);
  yaxis->SetLabelSize  (   0.1);
  yaxis->SetNdivisions (   505);
  yaxis->SetTitle      (ytitle);
  yaxis->SetTitleFont  (    42);
  yaxis->SetTitleOffset(  0.75);
  yaxis->SetTitleSize  (  0.11);
}



