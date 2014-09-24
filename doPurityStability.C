void doPurityStability(Int_t   jetChannel = 0, Int_t plot = 1, Int_t compare = 1) {


  gInterpreter->LoadMacro   ("../draw/draw.C+");

  gStyle      ->SetOptStat  (0);


  TString pathPow = Form("WW_GEN_0jet_pow_full_NNLL.root",  jetChannel); 
  TString pathMad = Form("WW_mad_test.root",  jetChannel); 
  TString pathMCNLO = Form("WW_pow_test.root",  jetChannel); 


  TFile *entryFilePow = new TFile(pathPow);
  TFile *entryFileMad = new TFile(pathMad);
  TFile *entryFileMCNLO = new TFile(pathMCNLO);

  entryFilePow->cd();

  bool doPt = 0;  
  bool doPtll = 0;
  bool doMinv =0;
  bool doPhi = 0;
  bool dojetEt = 0;
  bool doNjet = 0;

  if (plot == 0) doPt = 1; 
  if (plot == 1) doPtll = 1; 
  if (plot == 2) doMinv = 1; 
  if (plot == 3) doPhi = 1; 
  if (plot == 4) dojetEt = 1; 
  if (plot == 5) doNjet = 1;


  TLegend* legend = new TLegend(0.73,
			        0.70,
                                0.73+0.15,
				0.70+0.15);
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (   0.035);


  TLegend* legendMC = new TLegend(0.73,
			        0.70,
                                0.73+0.15,
				0.70+0.15);
  legendMC->SetBorderSize(    0);
  legendMC->SetFillColor (    0);
  legendMC->SetTextAlign (   12);
  legendMC->SetTextFont  (   42);
  legendMC->SetTextSize  (   0.035);


 


  //===================================================================
  // Do Pt 
  //===================================================================

  if (doPt) { 


    TH2F *fAij  = (TH2F *) hPtLepton1_RECO_GEN->Clone();

    TH1F *fgen_WW = (TH1F *) hPtLepton1WWLevel_GEN->Clone();
    TH2F *fAij_WW = (TH2F *) hPtLepton1WWLevel_RECO_GEN->Clone();
  

    TH2F *fAij_NormGEN = (TH2F *) hPtLepton1WWLevel_RECO_GEN->Clone();
    TH1F *f_NormGEN = (TH1F *) hPtLepton1WWLevel_GEN->Clone();

    TH2F *fAij_NormData = (TH2F *) hPtLepton1WWLevel_RECO_GEN->Clone();
    TH1F *f_NormData = (TH1F *) hPtLepton1WWLevel_RECO->Clone();


    TH1F *f_RECO = (TH1F *) hPtLepton1_RECO->Clone();
    TH1F *f_RECO_WW = (TH1F *) hPtLepton1WWLevel_RECO->Clone();
    
    TH1F *f_Selection_Eff = (TH1F *) hPtLepton1WWLevel_RECO->Clone();
    f_Selection_Eff->Divide(f_RECO_WW,f_RECO, 1.,1.,"");
    
    TH2F *fAij_Efficiency  = MultiplyEfficiency(fAij, f_Selection_Eff);


  //---------------------------------------------------------------------
  // -------- Do stability 
  //---------------------------------------------------------------------

   f_NormGEN = getStability(fAij);


   TCanvas* cStability = new TCanvas("stability" ,
                                "STABILITY + PURITY" ,750,750);

   cStability->cd();
   
   f_NormGEN->SetLineColor(kRed);
   f_NormGEN->SetMarkerColor(kRed);
   f_NormGEN->SetMarkerStyle(22);
   f_NormGEN->SetMarkerSize(2);
   f_NormGEN->Draw("E1");
   
   legend->AddEntry(f_NormGEN, "Stability");



 
  //---------------------------------------------------------------------
  // -------- Do purity
  //---------------------------------------------------------------------

  f_NormData = getPurity(fAij);


  cStability->cd();

  f_NormData->SetLineColor(kBlue);
  f_NormData->SetMarkerColor(kBlue);
  f_NormData->SetMarkerStyle(23);
  f_NormData->SetMarkerSize(2);
  f_NormData->Draw("sameE1");

  legend->AddEntry(f_NormData, "Purity");


  f_Selection_Eff->SetLineColor(kGreen+2);
  f_Selection_Eff->SetMarkerColor(kGreen+2);
  f_Selection_Eff->SetMarkerStyle(21);
  f_Selection_Eff->SetMarkerSize(2);
  f_Selection_Eff->Draw("sameE1");

  legend->AddEntry(f_Selection_Eff, "Efficiency");

  legend->Draw();

  } //----  end doPt 

  //===================================================================
  //===================================================================




  //===================================================================
  // Do Ptll
  //===================================================================

  if (doPtll) { 


    TH2F *fAij  = (TH2F *) hDilepton_RECO_GEN->Clone();

    TH1F *fgen_WW = (TH1F *) hDileptonWWLevel_GEN->Clone();
    TH2F *fAij_WW = (TH2F *) hDileptonWWLevel_RECO_GEN->Clone();
  

    TH2F *fAij_NormGEN = (TH2F *) hDileptonWWLevel_RECO_GEN->Clone();
    TH1F *f_NormGEN = (TH1F *) hDileptonWWLevel_GEN->Clone();

    TH2F *fAij_NormData = (TH2F *) hDileptonWWLevel_RECO_GEN->Clone();
    TH1F *f_NormData = (TH1F *) hDileptonWWLevel_RECO->Clone();


    TH1F *f_RECO = (TH1F *) hDilepton_RECO->Clone();
    TH1F *f_RECO_WW = (TH1F *) hDileptonWWLevel_RECO->Clone();
    
    TH1F *f_Selection_Eff = (TH1F *) hDileptonWWLevel_RECO->Clone();
    f_Selection_Eff->Divide(f_RECO_WW,f_RECO, 1.,1.,"");
    
    TH2F *fAij_Efficiency  = MultiplyEfficiency(fAij, f_Selection_Eff);


  //---------------------------------------------------------------------
  // -------- Do stability 
  //---------------------------------------------------------------------

   f_NormGEN = getStability(fAij_WW);


   TCanvas* cStability = new TCanvas("stability" ,
                                "STABILITY + PURITY" ,750,750);

   cStability->cd();
   
   f_NormGEN->SetLineColor(kRed);
   f_NormGEN->SetMarkerColor(kRed);
   f_NormGEN->SetMarkerStyle(22);
   f_NormGEN->SetMarkerSize(2);
   f_NormGEN->Draw("E1");
   
   legend->AddEntry(f_NormGEN, "Stability");


 
  //---------------------------------------------------------------------
  // -------- Do purity
  //---------------------------------------------------------------------

  f_NormData = getPurity(fAij_WW);


  cStability->cd();

  f_NormData->SetLineColor(kBlue);
  f_NormData->SetMarkerColor(kBlue);
  f_NormData->SetMarkerStyle(23);
  f_NormData->SetMarkerSize(2);
  f_NormData->Draw("sameE1");

  legend->AddEntry(f_NormData, "Purity");


  f_Selection_Eff->SetLineColor(kGreen+2);
  f_Selection_Eff->SetMarkerColor(kGreen+2);
  f_Selection_Eff->SetMarkerStyle(21);
  f_Selection_Eff->SetMarkerSize(2);
  f_Selection_Eff->Draw("sameE1");

  legend->AddEntry(f_Selection_Eff, "Efficiency");

  legend->Draw();

  } //----  end doPt 

  //===================================================================
  //===================================================================





  //===================================================================
  // Do Mll
  //===================================================================

  if (doMinv) { 


    TH2F *fAij  = (TH2F *) hmll_RECO_GEN->Clone();

    TH1F *fgen_WW = (TH1F *) hmllWWLevel_GEN->Clone();
    TH2F *fAij_WW = (TH2F *) hmllWWLevel_RECO_GEN->Clone();
  

    TH2F *fAij_NormGEN = (TH2F *) hmllWWLevel_RECO_GEN->Clone();
    TH1F *f_NormGEN = (TH1F *) hmllWWLevel_GEN->Clone();

    TH2F *fAij_NormData = (TH2F *) hmllWWLevel_RECO_GEN->Clone();
    TH1F *f_NormData = (TH1F *) hmllWWLevel_RECO->Clone();


    TH1F *f_RECO = (TH1F *) hmll_RECO->Clone();
    TH1F *f_RECO_WW = (TH1F *) hmllWWLevel_RECO->Clone();
    
    TH1F *f_Selection_Eff = (TH1F *) hmllWWLevel_RECO->Clone();
    f_Selection_Eff->Divide(f_RECO_WW,f_RECO, 1.,1.,"");
    
    TH2F *fAij_Efficiency  = MultiplyEfficiency(fAij, f_Selection_Eff);


  //---------------------------------------------------------------------
  // -------- Do stability 
  //---------------------------------------------------------------------

   f_NormGEN = getStability(fAij_WW);


   TCanvas* cStability = new TCanvas("stability" ,
                                "STABILITY + PURITY" ,750,750);

   cStability->cd();
   
   f_NormGEN->SetLineColor(kRed);
   f_NormGEN->SetMarkerColor(kRed);
   f_NormGEN->SetMarkerStyle(22);
   f_NormGEN->SetMarkerSize(2);
   f_NormGEN->Draw("E1");
   
   legend->AddEntry(f_NormGEN, "Stability");


 
  //---------------------------------------------------------------------
  // -------- Do purity
  //---------------------------------------------------------------------

  f_NormData = getPurity(fAij_WW);


  cStability->cd();

  f_NormData->SetLineColor(kBlue);
  f_NormData->SetMarkerColor(kBlue);
  f_NormData->SetMarkerStyle(23);
  f_NormData->SetMarkerSize(2);
  f_NormData->Draw("sameE1");

  legend->AddEntry(f_NormData, "Purity");


  f_Selection_Eff->SetLineColor(kGreen+2);
  f_Selection_Eff->SetMarkerColor(kGreen+2);
  f_Selection_Eff->SetMarkerStyle(21);
  f_Selection_Eff->SetMarkerSize(2);
  f_Selection_Eff->Draw("sameE1");

  legend->AddEntry(f_Selection_Eff, "Efficiency");

  legend->Draw();

  } //----  end doPt 

  //===================================================================
  //===================================================================






  //===================================================================
  // Do DeltaPhi
  //===================================================================

  if (doPhi) { 


    TH2F *fAij  = (TH2F *) hdphi_RECO_GEN->Clone();

    TH1F *fgen_WW = (TH1F *) hdphiWWLevel_GEN->Clone();
    TH2F *fAij_WW = (TH2F *) hdphiWWLevel_RECO_GEN->Clone();
  

    TH2F *fAij_NormGEN = (TH2F *) hdphiWWLevel_RECO_GEN->Clone();
    TH1F *f_NormGEN = (TH1F *) hdphiWWLevel_GEN->Clone();

    TH2F *fAij_NormData = (TH2F *) hdphiWWLevel_RECO_GEN->Clone();
    TH1F *f_NormData = (TH1F *) hdphiWWLevel_RECO->Clone();


    TH1F *f_RECO = (TH1F *) hdphi_RECO->Clone();
    TH1F *f_RECO_WW = (TH1F *) hdphiWWLevel_RECO->Clone();
    
    TH1F *f_Selection_Eff = (TH1F *) hdphiWWLevel_RECO->Clone();
    f_Selection_Eff->Divide(f_RECO_WW,f_RECO, 1.,1.,"");
    
    TH2F *fAij_Efficiency  = MultiplyEfficiency(fAij, f_Selection_Eff);


  //---------------------------------------------------------------------
  // -------- Do stability 
  //---------------------------------------------------------------------

   f_NormGEN = getStability(fAij_WW);


   TCanvas* cStability = new TCanvas("stability" ,
                                "STABILITY + PURITY" ,750,750);

   cStability->cd();
   
   f_NormGEN->SetLineColor(kRed);
   f_NormGEN->SetMarkerColor(kRed);
   f_NormGEN->SetMarkerStyle(22);
   f_NormGEN->SetMarkerSize(2);
   f_NormGEN->Draw("E1");
   
   legend->AddEntry(f_NormGEN, "Stability");


 
  //---------------------------------------------------------------------
  // -------- Do purity
  //---------------------------------------------------------------------

  f_NormData = getPurity(fAij_WW);


  cStability->cd();

  f_NormData->SetLineColor(kBlue);
  f_NormData->SetMarkerColor(kBlue);
  f_NormData->SetMarkerStyle(23);
  f_NormData->SetMarkerSize(2);
  f_NormData->Draw("sameE1");

  legend->AddEntry(f_NormData, "Purity");


  f_Selection_Eff->SetLineColor(kGreen+2);
  f_Selection_Eff->SetMarkerColor(kGreen+2);
  f_Selection_Eff->SetMarkerStyle(21);
  f_Selection_Eff->SetMarkerSize(2);
  f_Selection_Eff->Draw("sameE1");

  legend->AddEntry(f_Selection_Eff, "Efficiency");

  legend->Draw();

  } //----  end doPt 

  //===================================================================
  //===================================================================




  //===================================================================
  // Do JetEt
  //===================================================================

  if (dojetEt) { 


    TH2F *fAij  = (TH2F *) hjetEt_RECO_GEN->Clone();

    TH1F *fgen_WW = (TH1F *) hjetEtWWLevel_GEN->Clone();
    TH2F *fAij_WW = (TH2F *) hjetEtWWLevel_RECO_GEN->Clone();
  

    TH2F *fAij_NormGEN = (TH2F *) hjetEtWWLevel_RECO_GEN->Clone();
    TH1F *f_NormGEN = (TH1F *) hjetEtWWLevel_GEN->Clone();

    TH2F *fAij_NormData = (TH2F *) hjetEtWWLevel_RECO_GEN->Clone();
    TH1F *f_NormData = (TH1F *) hjetEtWWLevel_RECO->Clone();


    TH1F *f_RECO = (TH1F *) hjetEt_RECO->Clone();
    TH1F *f_RECO_WW = (TH1F *) hjetEtWWLevel_RECO->Clone();
    
    TH1F *f_Selection_Eff = (TH1F *) hjetEtWWLevel_RECO->Clone();
    f_Selection_Eff->Divide(f_RECO_WW,f_RECO, 1.,1.,"");
    
    TH2F *fAij_Efficiency  = MultiplyEfficiency(fAij, f_Selection_Eff);


  //---------------------------------------------------------------------
  // -------- Do stability 
  //---------------------------------------------------------------------

   f_NormGEN = getStability(fAij);


   TCanvas* cStability = new TCanvas("stability" ,
                                "STABILITY + PURITY" ,750,750);

   cStability->cd();
   
   f_NormGEN->SetLineColor(kRed);
   f_NormGEN->SetMarkerColor(kRed);
   f_NormGEN->SetMarkerStyle(22);
   f_NormGEN->SetMarkerSize(2);
   f_NormGEN->Draw("E1");
   
   legend->AddEntry(f_NormGEN, "Stability");



 
  //---------------------------------------------------------------------
  // -------- Do purity
  //---------------------------------------------------------------------

  f_NormData = getPurity(fAij);


  cStability->cd();

  f_NormData->SetLineColor(kBlue);
  f_NormData->SetMarkerColor(kBlue);
  f_NormData->SetMarkerStyle(23);
  f_NormData->SetMarkerSize(2);
  f_NormData->Draw("sameE1");

  legend->AddEntry(f_NormData, "Purity");


  f_Selection_Eff->SetLineColor(kGreen+2);
  f_Selection_Eff->SetMarkerColor(kGreen+2);
  f_Selection_Eff->SetMarkerStyle(21);
  f_Selection_Eff->SetMarkerSize(2);
  // f_Selection_Eff->Draw("sameE1");

  //legend->AddEntry(f_Selection_Eff, "Efficiency");

  legend->Draw();

  } //----  end dojetEt 

  //===================================================================
  //===================================================================





  //===================================================================
  // Do Njets
  //===================================================================

  if (doNjet) { 


    TH2F *fAij  = (TH2F *) hNjet_RECO_GEN->Clone();

    TH1F *f_RECO = (TH1F *) hNjet_RECO->Clone();
    TH1F *f_RECO_WW = (TH1F *) hNjet_RECO->Clone();
    
    TH1F *f_Selection_Eff = (TH1F *) hNjet_RECO->Clone();
    //f_Selection_Eff->Divide(f_RECO_WW,9.568, 1.,1.,""); 

    f_Selection_Eff->Scale(1./9568);

  //---------------------------------------------------------------------
  // -------- Do stability 
  //---------------------------------------------------------------------

   f_NormGEN = getStability(fAij);


   TCanvas* cStability = new TCanvas("stability" ,
                                "STABILITY + PURITY" ,750,750);

   cStability->cd();
   
   f_NormGEN->SetLineColor(kRed);
   f_NormGEN->SetMarkerColor(kRed);
   f_NormGEN->SetMarkerStyle(22);
   f_NormGEN->SetMarkerSize(2);
   f_NormGEN->Draw("E1");
   
   legend->AddEntry(f_NormGEN, "Stability");



 
  //---------------------------------------------------------------------
  // -------- Do purity
  //---------------------------------------------------------------------

  f_NormData = getPurity(fAij);


  cStability->cd();

  f_NormData->SetLineColor(kBlue);
  f_NormData->SetMarkerColor(kBlue);
  f_NormData->SetMarkerStyle(23);
  f_NormData->SetMarkerSize(2);
  f_NormData->Draw("sameE1");

  legend->AddEntry(f_NormData, "Purity");

  f_Selection_Eff->SetLineColor(kGreen+2);
  f_Selection_Eff->SetMarkerColor(kGreen+2);
  f_Selection_Eff->SetMarkerStyle(21);
  f_Selection_Eff->SetMarkerSize(2);
  //f_Selection_Eff->Draw("sameE1");

  //  legend->AddEntry(f_Selection_Eff, "Efficiency");

  legend->Draw();



  legend->Draw();






  } //----  end doNjet





  //---------------------------------------------------------------------
  // -------- Compare Purity/Stability with diffferent MCs 
  //---------------------------------------------------------------------
  
  if (compare) {

    TString distribution = ""; 

    if (doPt) distribution = "hPtLepton1_RECO_GEN";
    if (doMinv) distribution = "hmll_RECO_GEN";
    if (doPhi) distribution = "hdphi_RECO_GEN";
    if (doPtll) distribution = "hDilepton_RECO_GEN";
    if (doNjet) distribution = "hNjet_RECO_GEN";
  

    TH2F *fAij_Pow  = (TH2F *) entryFilePow->Get(distribution);
    TH2F *fAij_Mad  = (TH2F *) entryFileMad->Get(distribution);
    TH2F *fAij_MCNLO  = (TH2F *) entryFileMCNLO->Get(distribution);

    TH1F *h_PurityPow = getPurity(fAij_Pow);
    TH1F *h_PurityMad = getPurity(fAij_Mad);
    TH1F *h_PurityMCNLO = getPurity(fAij_MCNLO);

    TH1F *h_StabilityPow = getStability(fAij_Pow);
    TH1F *h_StabilityMad = getStability(fAij_Mad);
    TH1F *h_StabilityMCNLO = getStability(fAij_MCNLO);


   
    TH1F *h_RECO_WW_pow= (TH1F *) entryFilePow->Get("hdphiWWLevel_RECO");
    TH1F *h_RECO_pow= (TH1F *) entryFilePow->Get("hdphi_GEN");
    TH1F *plotEff_pow = h_RECO_pow;
    //plotEff_pow->Scale(1./8995);
    plotEff_pow->Divide(h_RECO_WW_pow,h_RECO_pow, 1.,1.,"");

    TH1F *h_RECO_WW_mad= (TH1F *) entryFileMad->Get("hdphiWWLevel_RECO");
    TH1F *h_RECO_mad= (TH1F *) entryFileMad->Get("hdphi_GEN");
    TH1F *plotEff_mad = h_RECO_mad;
    //plotEff_mad->Scale(1./9033);
    plotEff_mad->Divide(h_RECO_WW_mad,h_RECO_mad, 1.,1.,"");
    TH1F *h_RECO_WW_mcnlo= (TH1F *) entryFileMCNLO->Get("hdphiWWLevel_RECO");
    TH1F *h_RECO_mcnlo= (TH1F *) entryFileMCNLO->Get("hdphi_GEN");
    TH1F *plotEff_mcnlo = h_RECO_mcnlo; 
    //plotEff_mcnlo->Scale(1./8727);
    plotEff_mcnlo->Divide(h_RECO_WW_mcnlo,h_RECO_mcnlo, 1.,1.,"");
  

    TCanvas* cPurityMCs = new TCanvas("PurityMC" ,
                                "PURITY" ,750,750);

    h_PurityPow->SetLineColor(kBlue+2);
    h_PurityPow->SetMarkerColor(kBlue+2);
    h_PurityPow->SetMarkerStyle(20);
    h_PurityPow->SetMarkerSize(2);
    h_PurityPow->Draw("E1");

    h_PurityMad->SetLineColor(kRed+2);
    h_PurityMad->SetMarkerColor(kRed+2);
    h_PurityMad->SetMarkerStyle(21);
    h_PurityMad->SetMarkerSize(2);
    h_PurityMad->Draw("E1same");

    h_PurityMCNLO->SetLineColor(kGreen+2);
    h_PurityMCNLO->SetMarkerColor(kGreen+2);
    h_PurityMCNLO->SetMarkerStyle(22);
    h_PurityMCNLO->SetMarkerSize(2);
    h_PurityMCNLO->Draw("E1same");


    legendMC->AddEntry( h_PurityPow, "Powheg");
    legendMC->AddEntry( h_PurityMad, "Madgraph");
    legendMC->AddEntry( h_PurityMCNLO, "MCNLO");
 
    legendMC->Draw();


    TCanvas* cStabilityMCs = new TCanvas("StabilityMC" ,
                                "STABILITY" ,750,750);
    h_StabilityPow->SetLineColor(kBlue+2);
    h_StabilityPow->SetMarkerColor(kBlue+2);
    h_StabilityPow->SetMarkerStyle(20);
    h_StabilityPow->SetMarkerSize(2);
     h_StabilityPow->Draw("E1");

    h_StabilityMad->SetLineColor(kRed+2);
    h_StabilityMad->SetMarkerColor(kRed+2);
    h_StabilityMad->SetMarkerStyle(21);
    h_StabilityMad->SetMarkerSize(2);
    h_StabilityMad->Draw("E1same");

    h_StabilityMCNLO->SetLineColor(kGreen+2);
    h_StabilityMCNLO->SetMarkerColor(kGreen+2);
    h_StabilityMCNLO->SetMarkerStyle(22);
    h_StabilityMCNLO->SetMarkerSize(2);
    h_StabilityMCNLO->Draw("E1same");

    legendMC->Draw();


    TCanvas* cEfficiency = new TCanvas("cEfficiency" ,
				       "EFFICIENCY" ,750,750);
    
    plotEff_pow->SetLineColor(kBlue+2);
    plotEff_pow->SetMarkerColor(kBlue+2);
    plotEff_pow->SetMarkerStyle(20);
    plotEff_pow->SetMarkerSize(2);
    plotEff_pow->Draw("E1");

    plotEff_mad->SetLineColor(kRed+2);
    plotEff_mad->SetMarkerColor(kRed+2);
    plotEff_mad->SetMarkerStyle(21);
    plotEff_mad->SetMarkerSize(2);
    plotEff_mad->Draw("E1same");

    plotEff_mcnlo->SetLineColor(kGreen+2);
    plotEff_mcnlo->SetMarkerColor(kGreen+2);
    plotEff_mcnlo->SetMarkerStyle(22);
    plotEff_mcnlo->SetMarkerSize(2);
    plotEff_mcnlo->Draw("E1same");

    legendMC->Draw();

  }




}


TH1F* getStability(TH2F* hAij) {

  cout << "  " << endl;
  cout << "--------- STABILITY ---------" << endl;
  cout << "  " << endl;

  Int_t pt1_Aij_YBins =  hAij->GetNbinsY()+1;
  Int_t pt1_Aij_XBins =  hAij->GetNbinsX()+1;

  TH2F* h2DStability = (TH2F*) hAij->Clone();
  TH1F* h1DStability = (TH1F*) hAij->ProjectionX()->Clone();
   


 //------ DEFINE STABILITY -----//

  for ( int j = 1; j < pt1_Aij_YBins; j++) { 

    Float_t sumUpGEN = hAij->ProjectionY()->GetBinContent(j);
    Float_t sumUpGENerror = hAij->ProjectionY()->GetBinError(j);
    
    if ( sumUpGEN == 0) continue; 

    for ( int i = 1; i < pt1_Aij_XBins; i++) { 


      Float_t binData = hAij->GetBinContent(i,j); 
      Float_t binDataError = hAij->GetBinError(i,j);


     

      h2DStability->SetBinContent(i,j,binData/sumUpGEN);

      Float_t errorA = binDataError/sumUpGEN;
      Float_t errorB = (binData/(sumUpGEN*sumUpGEN))*sumUpGENerror;

      Float_t computeError = sqrt ( errorA*errorA + errorB*errorB);

      h2DStability->SetBinError(i,j,computeError);


      if ( i==j)  { 
	h1DStability->SetBinContent(j, binData/sumUpGEN);
	h1DStability->SetBinError(j, computeError);

	//	cout << sumUpGEN << "+-"  << sumUpGENerror << "  " << binData << "+-" << binDataError << endl;

	//ut << computeError << endl;
      }

    }

  }

  /*
  TCanvas* cStability2D = new TCanvas("stability2D" ,
                                "STABILITY 2D" ,750,750);

  cStability2D->cd();

  h2DStability->Draw("COLZ TEXT");
   */
  
  return h1DStability;

}


TH1F* getPurity(TH2F* hAij) {


  cout << "  " << endl;
  cout << "--------- PURITY ---------" << endl;
  cout << "  " << endl;



    Int_t pt1_Aij_YBins =  hAij->GetNbinsY()+1;
    Int_t pt1_Aij_XBins =  hAij->GetNbinsX()+1;

    TH2F* h2DPurity = (TH2F*) hAij->Clone();
    TH1F* h1DPurity = (TH1F*) hAij->ProjectionX()->Clone();
   

    //------ DEFINE PURITY -----//

    for ( int i = 1; i < pt1_Aij_XBins; i++) { 

      Float_t sumUpRECO = hAij->ProjectionX()->GetBinContent(i);
      Float_t sumUpRECOerror = hAij->ProjectionX()->GetBinError(i);
    
      if ( sumUpRECO == 0) continue; 

      for ( int j = 1; j < pt1_Aij_YBins; j++) { 


    	Float_t binData = hAij->GetBinContent(i,j); 
	Float_t binDataError = hAij->GetBinError(i,j);


	//cout << sumUpRECO << "  " << binData << endl;

	h2DPurity->SetBinContent(i,j,binData/sumUpRECO);

	Float_t errorA = binDataError/sumUpRECO;
	Float_t errorB = (binData/(sumUpRECO*sumUpRECO))*sumUpRECOerror;

	Float_t computeError = sqrt ( errorA*errorA + errorB*errorB);
      
	h2DPurity->SetBinError(i,j,computeError);

	if ( i==j)  { 
	  h1DPurity->SetBinContent(i, binData/sumUpRECO);
	  h1DPurity->SetBinError(i, computeError);
	}

      }

    }

  /*  
  TCanvas* cPurity = new TCanvas("purity" ,
                                "PURITY" ,750,750);

  cPurity->cd();
  h2DPurity->Draw("COLZ TEXT");
  */

    return h1DPurity;
}
