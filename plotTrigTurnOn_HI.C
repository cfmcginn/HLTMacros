#include <TFile.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
#include <TStyle.h>
#include "TH1F.h"
#include "TMath.h"
#include "TDatime.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"

#include <string>
#include <vector>
#include <fstream>

//ONLY EDITING SHOULD OCCUR HERE
const Int_t maxPlotTrig = 10;
const Int_t trigColors[10] = {1, kBlue, kRed, kYellow+1, kMagenta, kGreen+3, kCyan+2, kGray+1, kMagenta+2, kCyan-9};

void claverCanvasSaving(TCanvas* c, TString s,TString format="gif"){
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
  return;
}

int plotTrigTurnOn_HI(const std::string inHistFile, const std::string inTrigFileName, const std::string inOptionalEffFile = "")
{
  TH1::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  std::string buffer;
  std::vector<std::string> listOfTrig;
  std::vector<Int_t> listOfThresh;
  std::vector<std::string> listOfIsOpt;
  Int_t nTrigTypeTemp = 0;
  int nLines = 0;
  ifstream* inTrigFile = new ifstream(inTrigFileName.data());

  std::cout << inTrigFileName << std::endl;

  if(!inTrigFile->is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return 1;
  }
  else{
    while(true){
      *inTrigFile >> buffer;
      if(inTrigFile->eof()) break;
      if(std::string::npos== buffer.find("HLT") && std::string::npos== buffer.find("L1")) nTrigTypeTemp++;
      listOfTrig.push_back(buffer);
      nLines++;
    }
  }
  delete inTrigFile;
  std::cout << "Trigger List Loaded" << std::endl;

  const Int_t nTrigType = nTrigTypeTemp;
  std::string trigType[nTrigType];
  Int_t nTrigPlotTypeTemp = 0;
  std::string trigPlotType[nTrigType];
  std::string trigPlotType2[nTrigType];
  Int_t trigPlotThresh[nTrigType][2];
  Int_t trigPlotPos[nTrigType];
  Int_t trigPlotCol[nTrigType];

  for(Int_t iter = 0; iter < nTrigType; iter++){
    trigPlotPos[iter] = -1;
  }
  Int_t trigTypeCount[nTrigType];

  nLines = 0;

  for(Int_t iter = 0; iter < (Int_t)listOfTrig.size(); iter++){
    if(std::string::npos == listOfTrig[iter].find("HLT") && std::string::npos == listOfTrig[iter].find("L1")){

      std::size_t strIndex = 0;
      while(true){
        strIndex = listOfTrig[iter].find(",");
        if(strIndex == std::string::npos) break;
        trigPlotType[nLines] = listOfTrig[iter].substr(strIndex+1);
	Bool_t matchBool = false;
	for(Int_t iter2 = 0; iter2 < nLines; iter2++){
	  if(!strcmp(trigPlotType[nLines].c_str(), trigPlotType[iter2].c_str())){
	    matchBool = true;
	    break;
	  }
	}
	if(!matchBool) nTrigPlotTypeTemp++;
        listOfTrig[iter].replace(strIndex, std::string::npos, "");
      }

      trigType[nLines] = listOfTrig[iter];
      trigTypeCount[nLines] = 0;
      nLines++;
      listOfThresh.push_back(0);
      listOfIsOpt.push_back("");
    }
    else{
      trigTypeCount[nLines-1]++;

      std::size_t strIndex = 0;
      std::size_t strIndex2 = 0;
      while(true){
        strIndex = listOfTrig[iter].find(",");
        if(strIndex == std::string::npos) break;

	strIndex2 = listOfTrig[iter].substr(strIndex+1).find(",");
        listOfThresh.push_back(std::stoi(listOfTrig[iter].substr(strIndex+1, strIndex2)));
	listOfIsOpt.push_back(listOfTrig[iter].substr(strIndex+strIndex2+2));

        listOfTrig[iter].replace(strIndex, std::string::npos, "");
      }
    }
  }

  Int_t maxNTrig = -1;
  Int_t maxNum = 0;
  for(Int_t iter = 0; iter < nTrigType; iter++){
    if(trigTypeCount[iter] > maxNTrig) maxNTrig = trigTypeCount[iter];
  }

  const Int_t maxNTrig2 = maxNTrig;

  std::string trigName[nTrigType][maxNTrig2];
  Int_t trigThresh[nTrigType][maxNTrig2];
  std::string trigSpectIsOpt[nTrigType][maxNTrig2];

  nLines = 0;
  Int_t tempPosIter = 0;

  for(Int_t iter = 1; iter < (Int_t)(listOfTrig.size()); iter++){
    if(std::string::npos== listOfTrig[iter].find("HLT") && std::string::npos== listOfTrig[iter].find("L1")){
      nLines++;
      tempPosIter = 0;
    }
    else{
      trigName[nLines][tempPosIter] = listOfTrig[iter];
      trigThresh[nLines][tempPosIter] = listOfThresh[iter];
      trigSpectIsOpt[nLines][tempPosIter] = listOfIsOpt[iter];

      tempPosIter++;
    }
  }

  for(Int_t iter = 0; iter < nTrigType; iter++){
    if(trigPlotPos[iter] == -1){
      trigPlotPos[iter] = maxNum;
      trigPlotType2[maxNum] = trigPlotType[iter];
      trigPlotThresh[maxNum][0] = trigThresh[iter][0];
      trigPlotThresh[maxNum][1] = trigThresh[iter][trigTypeCount[iter] - 1];
      trigPlotCol[iter] = 1;
      maxNum++;
    }
    Int_t tempColPos = 1;
    for(Int_t iter2 = iter+1; iter2 < nTrigType; iter2++){
      if(!strcmp(trigPlotType[iter].c_str(), trigPlotType[iter2].c_str())){
	if(trigPlotPos[iter2] == -1) trigPlotCol[iter2] = trigColors[tempColPos];
	trigPlotPos[iter2] = trigPlotPos[iter];

	if(trigPlotThresh[maxNum][0] > trigThresh[iter][0]) trigPlotThresh[maxNum][0] = trigThresh[iter][0];
	if(trigPlotThresh[maxNum][1] < trigThresh[iter][trigTypeCount[iter] - 1]) trigPlotThresh[maxNum][1] = trigThresh[iter][trigTypeCount[iter] - 1];

	tempColPos++;
      }
    }
  }

  Bool_t isOpt = false;
  TFile* inOptFile_p;
  TGraphAsymmErrors* aPt_Opt_p[nTrigType][maxNTrig2];
  TH1F* maxMissedPt99_Opt_p[nTrigType];
  TH1F* maxMissedPt100_Opt_p[nTrigType];
  TH1F* leadingPtBoundSpect_Opt_p[nTrigType][maxNTrig2];
  TH1F* trkPtBoundSpect_Opt_p[nTrigType][maxNTrig2];
  if(strcmp(inOptionalEffFile.c_str(), "") != 0){
    isOpt = true;
    inOptFile_p = new TFile(inOptionalEffFile.c_str(), "READ");

    for(Int_t iter = 0; iter < nTrigType; iter++){
      std::cout << trigType[iter] << std::endl;
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(iter2 == maxPlotTrig) break;
	aPt_Opt_p[iter][iter2] = (TGraphAsymmErrors*)inOptFile_p->Get(Form("%sDir/%s_%s_pt_asymm", trigType[iter].c_str(), trigName[iter][iter2].c_str(), trigType[iter].c_str()));
	leadingPtBoundSpect_Opt_p[iter][iter2] = (TH1F*)inOptFile_p->Get(Form("%sDir/leadingPt_%s_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]));

	if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())) trkPtBoundSpect_Opt_p[iter][iter2] = (TH1F*)inOptFile_p->Get(Form("%sDir/trkPt_%s_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]));
      }

      maxMissedPt99_Opt_p[iter] = (TH1F*)inOptFile_p->Get(Form("%sDir/maxMissedPt99_%s_%d_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter] - 1]));
      maxMissedPt100_Opt_p[iter] = (TH1F*)inOptFile_p->Get(Form("%sDir/maxMissedPt100_%s_%d_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter] - 1]));
    }
  }

  TFile* inFile_p = new TFile(inHistFile.c_str(), "READ");
  inFile_p->cd();

  std::string canvPtName[nTrigType];
  TCanvas* trigCanvPt_p[nTrigType];
  TGraphAsymmErrors* aPt_p[nTrigType][maxNTrig2];

  const Int_t nTrigPlotType = nTrigPlotTypeTemp;
  std::string canvRateName[nTrigPlotType];
  TCanvas* trigCanvRate_p[nTrigPlotType];
  std::string canvSpectName[nTrigType];
  TCanvas* ptSpectBoundCanv_p[nTrigType];
  std::string canvSpectNameEff[nTrigType];
  TCanvas* ptSpectBoundCanv_Eff_p[nTrigType];
  std::string canvSpectNameEffLeadCut99[nTrigType];
  TCanvas* ptSpectBoundCanv_EffLeadCut99_p[nTrigType];
  std::string canvSpectNameEffLeadCut100[nTrigType];
  TCanvas* ptSpectBoundCanv_EffLeadCut100_p[nTrigType];

  std::string canvTrkSpectName[nTrigType];
  TCanvas* ptTrkSpectBoundCanv_p[nTrigType];
  std::string canvTrkSpectPrsclName[nTrigType];
  TCanvas* ptTrkSpectPrsclBoundCanv_p[nTrigType];
  std::string canvTrkSpectRelErrName[nTrigType];
  TCanvas* ptTrkSpectRelErrBoundCanv_p[nTrigType];

  std::string canvMaxMissedName[nTrigType];
  TCanvas* maxMissedPtCanv_p[nTrigType];
  std::string canvR9Name[nTrigType];
  TCanvas* R9Canv_p[nTrigType];

  TH1F* hEmptyRate[nTrigPlotType];
  Float_t maxXRate[nTrigPlotType];
  Float_t minXRate[nTrigPlotType];
  Float_t maxYRate[nTrigPlotType];
  Float_t minYRate[nTrigPlotType];
  Float_t binWidth[nTrigPlotType];
  TH1F* ratesUnmatched_p[nTrigType];
  TH1F* leadingPtBoundSpect_p[nTrigType][maxNTrig2];
  TH1F* leadingPtBoundSpect_Eff_p[nTrigType][maxNTrig2];
  TH1F* leadingPtBoundSpect_EffLeadCut99_p[nTrigType][maxNTrig2];
  TH1F* leadingPtBoundSpect_EffLeadCut100_p[nTrigType][maxNTrig2];
  TH1F* trkPtBoundSpect_p[nTrigType][maxNTrig2];
  TH1F* trkPtBoundSpectPrscl_p[nTrigType][maxNTrig2];
  TH1F* trkPtBoundSpectRelErr_p[nTrigType][maxNTrig2];
  TH1F* maxMissedPt99_p[nTrigType];
  TH1F* maxMissedPt100_p[nTrigType];
  TH1F* prescale_p[nTrigType];
  TH1F* e1R9_p[nTrigType];
  TH1F* e2R9_p[nTrigType];

  Float_t rateFract[nTrigType][maxNTrig];

  for(Int_t iter = 0; iter < nTrigPlotType; iter++){
    maxXRate[iter] = -1;
    minXRate[iter] = 100000;
    maxYRate[iter] = -1;
    minYRate[iter] = 10000000;
    binWidth[iter] = 10000000;
  }

  TLegend* trigLeg_p[nTrigType];
  TLegend* rateLeg_p[nTrigPlotType];
  TLegend* spectLeg_p[nTrigType];
  TLegend* spectLeg2_p[nTrigType];

  for(Int_t iter = 0; iter < nTrigType; iter++){      
    canvPtName[iter] = Form("%s_%d_%d_pt_c", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter] - 1]);
    trigCanvPt_p[iter] = new TCanvas(canvPtName[iter].c_str(), canvPtName[iter].c_str(), 700, 700);
    trigCanvPt_p[iter]->cd();
    trigCanvPt_p[iter]->SetTopMargin(0.01);
    trigCanvPt_p[iter]->SetRightMargin(0.01);
    trigLeg_p[iter] = new TLegend(0.65, 0.20, 0.98, 0.45);
    trigLeg_p[iter]->SetFillColor(0);
    trigLeg_p[iter]->SetTextFont(43);
    trigLeg_p[iter]->SetTextSize(16);

    std::cout << trigType[iter] << std::endl;
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      if(iter2 == maxPlotTrig) break;
      aPt_p[iter][iter2] = (TGraphAsymmErrors*)inFile_p->Get(Form("%sDir/%s_%s_pt_asymm", trigType[iter].c_str(), trigName[iter][iter2].c_str(), trigType[iter].c_str()));

      aPt_p[iter][iter2]->SetLineColor(trigColors[iter2]);
      aPt_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);

      rateFract[iter][iter2] = 0;

      if(iter2 == 0){
	TH1F* hEmpty = new TH1F("hEmpty", ";p_{T}^{reco};Efficiency", aPt_p[iter][iter2]->GetN(), aPt_p[iter][iter2]->GetXaxis()->GetXmin(), aPt_p[iter][iter2]->GetXaxis()->GetXmax());
	hEmpty->GetXaxis()->CenterTitle();
	hEmpty->GetYaxis()->CenterTitle();
	hEmpty->GetYaxis()->SetTitleOffset(1.4);
	hEmpty->SetMaximum(1.1);
	hEmpty->SetMinimum(0.0);
	hEmpty->DrawCopy();
	trigLeg_p[iter]->Draw("SAME");
	delete hEmpty;
      }

      aPt_p[iter][iter2]->Draw("P E");
      trigLeg_p[iter]->AddEntry(aPt_p[iter][iter2], Form("%s, %d", trigType[iter].c_str(), trigThresh[iter][iter2]), "P L");
    }

    TLine* oneLine_p = new TLine(aPt_p[iter][0]->GetXaxis()->GetXmin(), 1, aPt_p[iter][0]->GetXaxis()->GetXmax(), 1);
    oneLine_p->SetLineStyle(2);
    oneLine_p->DrawClone();

    delete oneLine_p;
  }

  for(Int_t iter = 0; iter < nTrigPlotType; iter++){
    canvRateName[iter] =  Form("%s_%d_%d_Rate_c", trigPlotType2[iter].c_str(), trigPlotThresh[iter][0], trigPlotThresh[iter][1]);

    trigCanvRate_p[iter] = new TCanvas(canvRateName[iter].c_str(), canvRateName[iter].c_str(), 700, 700);
    trigCanvRate_p[iter]->SetTopMargin(0.01);
    trigCanvRate_p[iter]->SetRightMargin(0.01);

    rateLeg_p[iter] = new TLegend(0.17, 0.15, 0.58, 0.40);
    rateLeg_p[iter]->SetFillColor(0);
    rateLeg_p[iter]->SetTextFont(43);
    rateLeg_p[iter]->SetTextSize(16);
    rateLeg_p[iter]->SetBorderSize(0);
  }

  for(Int_t iter = 0; iter < nTrigType; iter++){
    std::cout << "type1: " << trigType[iter] << std::endl;
    
    ratesUnmatched_p[iter] = (TH1F*)inFile_p->Get(Form("%sDir/ratesUnmatched_%s_%d_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]));

    if(maxXRate[trigPlotPos[iter]] < ratesUnmatched_p[iter]->GetXaxis()->GetXmax()) maxXRate[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetXaxis()->GetXmax();
    if(minXRate[trigPlotPos[iter]] > ratesUnmatched_p[iter]->GetXaxis()->GetXmin()) minXRate[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetXaxis()->GetXmin();

    if(maxYRate[trigPlotPos[iter]] < ratesUnmatched_p[iter]->GetMaximum()) maxYRate[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetMaximum();
    if(minYRate[trigPlotPos[iter]] > ratesUnmatched_p[iter]->GetMinimum()) minYRate[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetMinimum();

    if(binWidth[trigPlotPos[iter]] > ratesUnmatched_p[iter]->GetXaxis()->GetBinWidth(1)/2.0) binWidth[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetXaxis()->GetBinWidth(1)/2.0;

    prescale_p[iter] = (TH1F*)inFile_p->Get(Form("%sDir/prescale_%s_%d_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter] - 1]));

    e1R9_p[iter] = (TH1F*)inFile_p->Get(Form("%sDir/e1R9_%s_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][0]));
    e2R9_p[iter] = (TH1F*)inFile_p->Get(Form("%sDir/e2R9_%s_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][0]));

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      inFile_p->cd();
      if(isOpt && !strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Opt")){
	leadingPtBoundSpect_p[iter][iter2] = (TH1F*)leadingPtBoundSpect_Opt_p[iter][iter2]->Clone(Form("leadingPt_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));

	leadingPtBoundSpect_Eff_p[iter][iter2] = (TH1F*)leadingPtBoundSpect_Opt_p[iter][iter2]->Clone(Form("leadingPt_Eff_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));
	
	leadingPtBoundSpect_EffLeadCut99_p[iter][iter2] = (TH1F*)leadingPtBoundSpect_Opt_p[iter][iter2]->Clone(Form("leadingPt_EffLeadCut99_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));

	leadingPtBoundSpect_EffLeadCut100_p[iter][iter2] = (TH1F*)leadingPtBoundSpect_Opt_p[iter][iter2]->Clone(Form("leadingPt_EffLeadCut100_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));
	if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())){
	  trkPtBoundSpect_p[iter][iter2] = (TH1F*)trkPtBoundSpect_Opt_p[iter][iter2]->Clone(Form("trkPt_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));
	  trkPtBoundSpectPrscl_p[iter][iter2] = (TH1F*)trkPtBoundSpect_p[iter][iter2]->Clone(Form("trkPt_Prscl_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));
	  trkPtBoundSpectRelErr_p[iter][iter2] = (TH1F*)trkPtBoundSpect_p[iter][iter2]->Clone(Form("trkPt_RelErr_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));
	}
      }
      else{
	leadingPtBoundSpect_p[iter][iter2] = (TH1F*)inFile_p->Get(Form("%sDir/leadingPt_%s_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]));
	leadingPtBoundSpect_Eff_p[iter][iter2] = (TH1F*)leadingPtBoundSpect_p[iter][iter2]->Clone(Form("leadingPt_Eff_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));
	leadingPtBoundSpect_EffLeadCut99_p[iter][iter2] = (TH1F*)leadingPtBoundSpect_p[iter][iter2]->Clone(Form("leadingPt_EffLeadCut99_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));
	leadingPtBoundSpect_EffLeadCut100_p[iter][iter2] = (TH1F*)leadingPtBoundSpect_p[iter][iter2]->Clone(Form("leadingPt_EffLeadCut100_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));

	if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())){
	  trkPtBoundSpect_p[iter][iter2] = (TH1F*)inFile_p->Get(Form("%sDir/trkPt_%s_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]));
	  trkPtBoundSpectPrscl_p[iter][iter2] = (TH1F*)trkPtBoundSpect_p[iter][iter2]->Clone(Form("trkPt_Prscl_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));
	  trkPtBoundSpectRelErr_p[iter][iter2] = (TH1F*)trkPtBoundSpect_p[iter][iter2]->Clone(Form("trkPt_RelErr_%s_%d", trigType[iter].c_str(), trigThresh[iter][iter2]));
	}
      }

      if(leadingPtBoundSpect_p[iter][iter2]->Integral() != 0){
	if(isOpt && strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Not") != 0){
	  for(Int_t binIter = 0; binIter < aPt_Opt_p[iter][iter2]->GetN(); binIter++){
	    Double_t tempX = 0.0001;
	    Double_t tempEff = 0.0001;
	    aPt_Opt_p[iter][iter2]->GetPoint(binIter, tempX, tempEff);
	    
	    if(tempEff < 0.0001) tempEff = 0.0001;
	    
	    if(tempX < leadingPtBoundSpect_p[iter][iter2]->GetXaxis()->GetXmin() || tempX > leadingPtBoundSpect_p[iter][iter2]->GetXaxis()->GetXmax()) continue;
	    
	    Int_t bin = leadingPtBoundSpect_p[iter][iter2]->FindBin(tempX);
	    
	    leadingPtBoundSpect_Eff_p[iter][iter2]->SetBinContent(bin, leadingPtBoundSpect_Eff_p[iter][iter2]->GetBinContent(bin)*tempEff);
	    leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetBinContent(bin, leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->GetBinContent(bin)*tempEff);
	    leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetBinContent(bin, leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->GetBinContent(bin)*tempEff);
	  }
	}		
	inFile_p->cd();
	
	for(Int_t binIter = 0; binIter < leadingPtBoundSpect_p[iter][iter2]->GetNbinsX(); binIter++){
          if(leadingPtBoundSpect_p[iter][iter2]->GetBinLowEdge(binIter+2) > trigThresh[iter][iter2]) break;

          leadingPtBoundSpect_p[iter][iter2]->SetBinContent(binIter+1, 0);
        }

	Float_t rate = ratesUnmatched_p[iter]->GetBinContent(ratesUnmatched_p[iter]->FindBin(trigThresh[iter][iter2]));
	
	Int_t prscl = prescale_p[iter]->GetBinContent(prescale_p[iter]->FindBin(trigThresh[iter][iter2]));

	if(prscl >= 30) rate *= 800000;
	else rate *= 800000/3;

	leadingPtBoundSpect_p[iter][iter2]->Scale(1/leadingPtBoundSpect_p[iter][iter2]->Integral());
	leadingPtBoundSpect_p[iter][iter2]->Scale(rate);
	
	std::cout << "CHECK THIS: " << rate << ", " << leadingPtBoundSpect_p[iter][iter2]->Integral() << std::endl;

	leadingPtBoundSpect_Eff_p[iter][iter2]->Scale(1/leadingPtBoundSpect_Eff_p[iter][iter2]->Integral());
	leadingPtBoundSpect_Eff_p[iter][iter2]->Scale(rate);
	
	leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->Scale(1/leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->Integral());
	leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->Scale(rate);
	
	leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->Scale(1/leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->Integral());
	leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->Scale(rate);
    

	if(isOpt){
	  Float_t maxCut99 = maxMissedPt99_Opt_p[iter]->GetBinContent(maxMissedPt99_Opt_p[iter]->FindBin(trigThresh[iter][iter2]));
	  Int_t maxBin99 = leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->FindBin(maxCut99);

	  if(leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->Integral() > .1) rateFract[iter][iter2] = leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->Integral(maxBin99, leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->GetNbinsX()+1)/leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->Integral();
	  else rateFract[iter][iter2] = 0;

	  if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())) std::cout << rateFract[iter][iter2] << std::endl;

	  for(Int_t binIter = 0; binIter < maxBin99; binIter++){
	    leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetBinContent(binIter+1, 0);
	  }
	  
	  Float_t maxCut100 = maxMissedPt100_Opt_p[iter]->GetBinContent(maxMissedPt100_Opt_p[iter]->FindBin(trigThresh[iter][iter2]));
	  Int_t maxBin = leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->FindBin(maxCut100);
	  for(Int_t binIter = 0; binIter < maxBin; binIter++){
	    leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetBinContent(binIter+1, 0);
	  }
	}

	if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())){
	  trkPtBoundSpect_p[iter][iter2]->Scale(rate*rateFract[iter][iter2]);
	  trkPtBoundSpectPrscl_p[iter][iter2]->Scale(rate*rateFract[iter][iter2]);

	  Float_t prscl = prescale_p[iter]->GetBinContent(prescale_p[iter]->FindBin(trigThresh[iter][iter2]));
	  for(Int_t binIter = 0; binIter < trkPtBoundSpectPrscl_p[iter][iter2]->GetNbinsX(); binIter++){
	    trkPtBoundSpectPrscl_p[iter][iter2]->SetBinContent(binIter+1, trkPtBoundSpectPrscl_p[iter][iter2]->GetBinContent(binIter+1)*prscl);
	  }

	  std::cout << "rateFract: " <<  rateFract[iter][iter2] << std::endl;
	}
	
      }
    }

    if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())){
      const Int_t nPrsclBins = trkPtBoundSpectPrscl_p[iter][0]->GetNbinsX();
      for(Int_t binIter = 0; binIter < nPrsclBins; binIter++){
	Float_t tempMax = 0;
	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Not") != 0){
	    if(trkPtBoundSpectPrscl_p[iter][iter2]->GetBinContent(binIter+1) > tempMax){
	      trkPtBoundSpectRelErr_p[iter][iter2]->SetBinContent(binIter+1, 1./TMath::Sqrt(trkPtBoundSpect_p[iter][iter2]->GetBinContent(binIter+1)));
	      tempMax = trkPtBoundSpectPrscl_p[iter][iter2]->GetBinContent(binIter+1);
	      for(Int_t iter3 = 0; iter3 < trigTypeCount[iter]; iter3++){
		if(iter3 == iter2) continue;
		trkPtBoundSpectRelErr_p[iter][iter3]->SetBinContent(binIter+1, 0);
	      }
	    }
	    trkPtBoundSpectRelErr_p[iter][iter2]->SetBinError(binIter+1, 0);
	  }
	}

      }
    }

    maxMissedPt99_p[iter] = (TH1F*)inFile_p->Get(Form("%sDir/maxMissedPt99_%s_%d_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter] - 1]));
    maxMissedPt100_p[iter] = (TH1F*)inFile_p->Get(Form("%sDir/maxMissedPt100_%s_%d_%d", trigType[iter].c_str(), trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter] - 1]));
  }

  for(Int_t iter = 0; iter < nTrigPlotType; iter++){
    trigCanvRate_p[iter]->cd();

    maxYRate[iter] += 10*TMath::Sqrt(maxYRate[iter]);
    if(minYRate[iter] < .00001) minYRate[iter] = 0.1;
    else if(minYRate[iter] > 1) minYRate[iter] -= TMath::Sqrt(minYRate[iter]);
    else minYRate[iter] -= minYRate[iter]*minYRate[iter];

    if(minYRate[iter] > 1) minYRate[iter] = 1;


    hEmptyRate[iter] = new TH1F(Form("hEmptyRate%d", iter), Form(";p_{T,%s}^{trig};Rate (Hz)", trigPlotType2[iter].c_str()), 10, minXRate[iter], maxXRate[iter]);
    hEmptyRate[iter]->GetXaxis()->CenterTitle();
    hEmptyRate[iter]->GetYaxis()->CenterTitle();
    hEmptyRate[iter]->GetYaxis()->SetTitleOffset(1.4);
    hEmptyRate[iter]->SetMaximum(maxYRate[iter]);
    hEmptyRate[iter]->SetMinimum(minYRate[iter]);

    hEmptyRate[iter]->DrawCopy();
    rateLeg_p[iter]->Draw("SAME");
    delete hEmptyRate[iter];
    gPad->SetLogy();
  }

  for(Int_t iter = 0; iter < nTrigType; iter++){
    trigCanvRate_p[trigPlotPos[iter]]->cd();

    ratesUnmatched_p[iter]->SetMarkerStyle(20);
    ratesUnmatched_p[iter]->SetMarkerSize(1);
    ratesUnmatched_p[iter]->SetLineColor(trigPlotCol[iter]);
    ratesUnmatched_p[iter]->SetMarkerColor(trigPlotCol[iter]);

    gStyle->SetErrorX(.5*binWidth[trigPlotPos[iter]]/(ratesUnmatched_p[iter]->GetXaxis()->GetBinWidth(1)/2.0));
    ratesUnmatched_p[iter]->DrawCopy("E1 SAME");

    rateLeg_p[trigPlotPos[iter]]->AddEntry(ratesUnmatched_p[iter], Form("%s", trigType[iter].c_str()), "P L");

    TLine* tenLine_p = new TLine(minXRate[trigPlotPos[iter]], 10, maxXRate[trigPlotPos[iter]], 10);
    tenLine_p->SetLineStyle(2);
    if(10 > minYRate[trigPlotPos[iter]] && 10 < maxYRate[trigPlotPos[iter]]) tenLine_p->DrawClone();

    if(1 > minYRate[trigPlotPos[iter]] && 1 < maxYRate[trigPlotPos[iter]]) tenLine_p->DrawLine(minXRate[trigPlotPos[iter]], 1, maxXRate[trigPlotPos[iter]], 1);

    delete tenLine_p;
  }

  for(Int_t iter = 0; iter < nTrigType; iter++){
    canvSpectName[iter] =  Form("leadingPtBoundSpect_%s_c", trigType[iter].c_str());
    canvSpectNameEff[iter] =  Form("leadingPtBoundSpect_Eff_%s_c", trigType[iter].c_str());
    canvSpectNameEffLeadCut99[iter] =  Form("leadingPtBoundSpect_EffLeadCut99_%s_c", trigType[iter].c_str());
    canvSpectNameEffLeadCut100[iter] =  Form("leadingPtBoundSpect_EffLeadCut100_%s_c", trigType[iter].c_str());

    canvTrkSpectName[iter] =  Form("trkPtBoundSpect_%s_c", trigType[iter].c_str());
    canvTrkSpectPrsclName[iter] =  Form("trkPtBoundSpectPrscl_%s_c", trigType[iter].c_str());
    canvTrkSpectRelErrName[iter] =  Form("trkPtBoundSpectRelErr_%s_c", trigType[iter].c_str());

    ptSpectBoundCanv_p[iter] = new TCanvas(canvSpectName[iter].c_str(), canvSpectName[iter].c_str(), 700, 700);
    ptSpectBoundCanv_p[iter]->SetTopMargin(0.01);
    ptSpectBoundCanv_p[iter]->SetRightMargin(0.01);

    ptSpectBoundCanv_p[iter]->cd();
    gPad->SetLogy();

    Float_t tempMax = -1;
    Float_t tempLowest = 100000000000;
    Bool_t foundLowBin = false;
    Int_t tempLowBin = 1000000;

    spectLeg_p[iter] = new TLegend(0.25, 0.20, 0.45, 0.45);
    spectLeg_p[iter]->SetFillColor(0);
    spectLeg_p[iter]->SetTextFont(43);
    spectLeg_p[iter]->SetTextSize(16);
    spectLeg_p[iter]->SetFillStyle(0);
    spectLeg_p[iter]->SetBorderSize(0);

    spectLeg2_p[iter] = new TLegend(0.15, 0.50, 0.45, 0.85);
    spectLeg2_p[iter]->SetFillColor(0);
    spectLeg2_p[iter]->SetTextFont(43);
    spectLeg2_p[iter]->SetTextSize(16);
    spectLeg2_p[iter]->SetFillStyle(0);
    spectLeg2_p[iter]->SetBorderSize(0);

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      if(strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Not") != 0){

	if(iter2 != 0) leadingPtBoundSpect_p[iter][iter2]->Add(leadingPtBoundSpect_p[iter][iter2-1]);

	spectLeg_p[iter]->AddEntry(leadingPtBoundSpect_p[iter][iter2], Form("%s %d", trigType[iter].c_str(), trigThresh[iter][iter2]), "F");

	for(Int_t binIter = 0; binIter < leadingPtBoundSpect_p[iter][iter2]->GetNbinsX(); binIter++){
	  if(leadingPtBoundSpect_p[iter][iter2]->GetBinLowEdge(binIter+2) > trigThresh[iter][iter2]) break;

	  leadingPtBoundSpect_p[iter][iter2]->SetBinContent(binIter+1, 0);
	}
	
      	if(leadingPtBoundSpect_p[iter][iter2]->GetMaximum() > tempMax) tempMax = leadingPtBoundSpect_p[iter][iter2]->GetMaximum();


	leadingPtBoundSpect_p[iter][iter2]->SetMinimum(1);
	leadingPtBoundSpect_p[iter][iter2]->GetXaxis()->CenterTitle();
	leadingPtBoundSpect_p[iter][iter2]->GetYaxis()->CenterTitle();
	leadingPtBoundSpect_p[iter][iter2]->GetYaxis()->SetTitleOffset(leadingPtBoundSpect_p[iter][iter2]->GetYaxis()->GetTitleOffset()+.10);
      
	if(iter2 == 0){
	  leadingPtBoundSpect_p[iter][iter2]->SetMarkerColor(17);
	  leadingPtBoundSpect_p[iter][iter2]->SetLineColor(1);
	  leadingPtBoundSpect_p[iter][iter2]->SetFillColor(17);
	}
	else{
	  leadingPtBoundSpect_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);
	  leadingPtBoundSpect_p[iter][iter2]->SetLineColor(1);
	  leadingPtBoundSpect_p[iter][iter2]->SetFillColor(trigColors[iter2]);
	}
      }
    }

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){

	Int_t nextMaxBin = -1;
	Bool_t isFound = false;
	if(!foundLowBin){
	  tempLowBin = leadingPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetNbinsX();
	  foundLowBin = true;
	}



	for(Int_t binIter = 0; binIter < tempLowBin; binIter++){
	  if(leadingPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < tempLowest && leadingPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-iter2].c_str(), "Not") != 0){

            Bool_t contBool = false;

            for(Int_t iter3 = iter2+1; iter3 < trigTypeCount[iter]-1; iter3++){
              if(leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter3]->GetBinContent(binIter+1)){
                contBool = true;
                break;
              }
            }

            if(contBool) continue;
	    tempLowest = leadingPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1);
	  }

	  if(leadingPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && !isFound){
	    nextMaxBin = binIter+1;
	    isFound = true;
	  }
	}

	tempLowBin = nextMaxBin-1;
      }
    }

    Bool_t isDrawn = false;
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){
	leadingPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->SetMaximum(tempMax);

	if(!isDrawn) leadingPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST");
	else leadingPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST SAME");

	isDrawn = true;
      }
    }
    spectLeg_p[iter]->Draw("SAME");
    TLine* line_p = new TLine(leadingPtBoundSpect_p[iter][0]->GetXaxis()->GetXmin(), tempLowest, leadingPtBoundSpect_p[iter][0]->GetXaxis()->GetXmax(), tempLowest);
    line_p->SetLineStyle(2);
    line_p->DrawLine(leadingPtBoundSpect_p[iter][0]->GetXaxis()->GetXmin(), tempLowest, leadingPtBoundSpect_p[iter][0]->GetXaxis()->GetXmax(), tempLowest);
    tempLowest = 100000000;
    foundLowBin = false;
    tempLowBin = 1000000;


    if(isOpt){
      ptSpectBoundCanv_Eff_p[iter] = new TCanvas(canvSpectNameEff[iter].c_str(), canvSpectNameEff[iter].c_str(), 700, 700);
      ptSpectBoundCanv_Eff_p[iter]->SetTopMargin(0.01);
      ptSpectBoundCanv_Eff_p[iter]->SetRightMargin(0.01);
      
      ptSpectBoundCanv_Eff_p[iter]->cd();
      gPad->SetLogy();
      
      tempMax = -1;
      
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Not") != 0){
	  
	  if(iter2 != 0) leadingPtBoundSpect_Eff_p[iter][iter2]->Add(leadingPtBoundSpect_Eff_p[iter][iter2-1]);
	  
	  for(Int_t binIter = 0; binIter < leadingPtBoundSpect_Eff_p[iter][iter2]->GetNbinsX(); binIter++){
	    if(leadingPtBoundSpect_Eff_p[iter][iter2]->GetBinLowEdge(binIter+2) > trigThresh[iter][iter2]*.8) break;
	    
	    leadingPtBoundSpect_Eff_p[iter][iter2]->SetBinContent(binIter+1, 0);
	  }
	  
	  if(leadingPtBoundSpect_Eff_p[iter][iter2]->GetMaximum() > tempMax) tempMax = leadingPtBoundSpect_Eff_p[iter][iter2]->GetMaximum();
	  
	  
	  //Continue here
	  leadingPtBoundSpect_Eff_p[iter][iter2]->SetMinimum(1);
	  leadingPtBoundSpect_Eff_p[iter][iter2]->GetXaxis()->CenterTitle();
	  leadingPtBoundSpect_Eff_p[iter][iter2]->GetYaxis()->CenterTitle();
	  leadingPtBoundSpect_Eff_p[iter][iter2]->GetYaxis()->SetTitleOffset(leadingPtBoundSpect_Eff_p[iter][iter2]->GetYaxis()->GetTitleOffset()+.10);
	  
	  if(iter2 == 0){
	    leadingPtBoundSpect_Eff_p[iter][iter2]->SetMarkerColor(17);
	    leadingPtBoundSpect_Eff_p[iter][iter2]->SetLineColor(1);
	    leadingPtBoundSpect_Eff_p[iter][iter2]->SetFillColor(17);
	  }
	  else{
	    leadingPtBoundSpect_Eff_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);
	    leadingPtBoundSpect_Eff_p[iter][iter2]->SetLineColor(1);
	    leadingPtBoundSpect_Eff_p[iter][iter2]->SetFillColor(trigColors[iter2]);
	  }
	}
      }

      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){

	  Int_t nextMaxBin = -1;
	  Bool_t isFound = false;
	  if(!foundLowBin){
	    tempLowBin = leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->GetNbinsX();
	    foundLowBin = true;
	  }

	  for(Int_t binIter = 0; binIter < tempLowBin; binIter++){
	    if(leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < tempLowest && leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-iter2].c_str(), "Not") != 0){
	      Bool_t contBool = false;

	      for(Int_t iter3 = iter2+1; iter3 < trigTypeCount[iter]-1; iter3++){
		if(leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter3]->GetBinContent(binIter+1)){
		  contBool = true;
		  break;
		}
	      }

	      if(contBool) continue;
	      tempLowest = leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1);

	    }

	    if(leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && !isFound){
	      nextMaxBin = binIter+1;
	      isFound = true;
	    }
	  }

	  tempLowBin = nextMaxBin-1;
	}
      }
      
      isDrawn = false;
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){
	  leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->SetMaximum(tempMax);
	  
	if(!isDrawn) leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST");
	else leadingPtBoundSpect_Eff_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST SAME");
	isDrawn = true;
	}
      }
      spectLeg_p[iter]->Draw("SAME");      
      line_p->DrawLine(leadingPtBoundSpect_Eff_p[iter][0]->GetXaxis()->GetXmin(), tempLowest, leadingPtBoundSpect_Eff_p[iter][0]->GetXaxis()->GetXmax(), tempLowest);
      std::cout << trigType[iter] << ", " << tempLowest << std::endl;
      tempLowest = 10000000000;
      foundLowBin = false;
      tempLowBin = 1000000;



      ptSpectBoundCanv_EffLeadCut99_p[iter] = new TCanvas(canvSpectNameEffLeadCut99[iter].c_str(), canvSpectNameEffLeadCut99[iter].c_str(), 700, 700);
      ptSpectBoundCanv_EffLeadCut99_p[iter]->SetTopMargin(0.01);
      ptSpectBoundCanv_EffLeadCut99_p[iter]->SetRightMargin(0.01);
      
      ptSpectBoundCanv_EffLeadCut99_p[iter]->cd();
      gPad->SetLogy();
      
      tempMax = -1;
      
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Not") != 0){
	  
	  if(iter2 != 0) leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->Add(leadingPtBoundSpect_EffLeadCut99_p[iter][iter2-1]);
	  
	  Float_t maxCut99 = maxMissedPt99_Opt_p[iter]->GetBinContent(maxMissedPt99_Opt_p[iter]->FindBin(trigThresh[iter][iter2]));
	  Int_t maxBin99 = leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->FindBin(maxCut99);

	  for(Int_t binIter = 0; binIter < maxBin99; binIter++){
	    leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetBinContent(binIter+1, 0);
	  }
	  
	  if(leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->GetMaximum() > tempMax) tempMax = leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->GetMaximum();
	  
	  leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetMinimum(1);
	  leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->GetXaxis()->CenterTitle();
	  leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->GetYaxis()->CenterTitle();
	  leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->GetYaxis()->SetTitleOffset(leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->GetYaxis()->GetTitleOffset()+.10);
	  
	  if(iter2 == 0){
	    leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetMarkerColor(17);
	    leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetLineColor(1);
	    leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetFillColor(17);
	  }
	  else{
	    leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);
	    leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetLineColor(1);
	    leadingPtBoundSpect_EffLeadCut99_p[iter][iter2]->SetFillColor(trigColors[iter2]);
	  }
	}
      }

      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){

	  Int_t nextMaxBin = -1;
	  Bool_t isFound = false;
	  if(!foundLowBin){
	    tempLowBin = leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2]->GetNbinsX();
	    foundLowBin = true;
	  }

	  for(Int_t binIter = 0; binIter < tempLowBin; binIter++){
	    if(leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < tempLowest && leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-iter2].c_str(), "Not") != 0){
	      Bool_t contBool = false;

	      for(Int_t iter3 = iter2+1; iter3 < trigTypeCount[iter]-1; iter3++){
		if(leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter3]->GetBinContent(binIter+1)){
		  contBool = true;
		  break;
		}
	      }

	      if(contBool) continue;
	      tempLowest = leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1);

	    
	      if(iter2 != trigTypeCount[iter]-1) std::cout << "     " << binIter << ", " << tempLowest << ", " << leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2-1]->GetBinContent(binIter+1)  << std::endl;
	      else std::cout << "     " << binIter << ", " << tempLowest << std::endl;
	    }

	    if(leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && !isFound){
	      nextMaxBin = binIter+1;
	      isFound = true;
	    }
	  }

	  tempLowBin = nextMaxBin-1;
	}
      }

      
      isDrawn = false;
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){
	  leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2]->SetMaximum(tempMax);
	  
	  if(!isDrawn) leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST");
	  else leadingPtBoundSpect_EffLeadCut99_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST SAME");
	  isDrawn = true;
	}
      }
      spectLeg_p[iter]->Draw("SAME");
      line_p->DrawLine(leadingPtBoundSpect_EffLeadCut99_p[iter][0]->GetXaxis()->GetXmin(), tempLowest, leadingPtBoundSpect_EffLeadCut99_p[iter][0]->GetXaxis()->GetXmax(), tempLowest);
      std::cout << trigType[iter] << ", " << tempLowest << std::endl;
      tempLowest = 10000000000;
      foundLowBin = false;
      tempLowBin = 1000000;

      
      ptSpectBoundCanv_EffLeadCut100_p[iter] = new TCanvas(canvSpectNameEffLeadCut100[iter].c_str(), canvSpectNameEffLeadCut100[iter].c_str(), 700, 700);
      ptSpectBoundCanv_EffLeadCut100_p[iter]->SetTopMargin(0.01);
      ptSpectBoundCanv_EffLeadCut100_p[iter]->SetRightMargin(0.01);
      
      ptSpectBoundCanv_EffLeadCut100_p[iter]->cd();
      gPad->SetLogy();
      
      tempMax = -1;
      
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Not") != 0){
	  
	  if(iter2 != 0) leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->Add(leadingPtBoundSpect_EffLeadCut100_p[iter][iter2-1]);
	  
	  
	  Float_t maxCut100 = maxMissedPt100_Opt_p[iter]->GetBinContent(maxMissedPt100_Opt_p[iter]->FindBin(trigThresh[iter][iter2]));
	  Int_t maxBin = leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->FindBin(maxCut100);
	  for(Int_t binIter = 0; binIter < maxBin; binIter++){
	    leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetBinContent(binIter+1, 0);
	  }
	  
	  if(leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->GetMaximum() > tempMax) tempMax = leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->GetMaximum();
	  
	  leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetMinimum(1);
	  leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->GetXaxis()->CenterTitle();
	  leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->GetYaxis()->CenterTitle();
	  leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->GetYaxis()->SetTitleOffset(leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->GetYaxis()->GetTitleOffset()+.10);
	  
	  if(iter2 == 0){
	    leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetMarkerColor(17);
	    leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetLineColor(1);
	    leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetFillColor(17);
	  }
	  else{
	    leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);
	    leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetLineColor(1);
	    leadingPtBoundSpect_EffLeadCut100_p[iter][iter2]->SetFillColor(trigColors[iter2]);
	  }
	}
      }

      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){

	  Int_t nextMaxBin = -1;
	  Bool_t isFound = false;
	  if(!foundLowBin){
	    tempLowBin = leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2]->GetNbinsX();
	    foundLowBin = true;
	 } 

	  std::cout << "  " << trigType[iter] << ", " << tempLowBin << std::endl;

	  for(Int_t binIter = 0; binIter < tempLowBin; binIter++){
	    if(leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < tempLowest && leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-iter2].c_str(), "Not") != 0){
	      Bool_t contBool = false;

	      for(Int_t iter3 = iter2+1; iter3 < trigTypeCount[iter]-1; iter3++){
		if(leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter3]->GetBinContent(binIter+1)){
		  contBool = true;
		  break;
		}
	      }

	      if(contBool) continue;
	      tempLowest = leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1);

	    
	      if(iter2 != trigTypeCount[iter]-1) std::cout << "     " << binIter << ", " << tempLowest << ", " << leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2-1]->GetBinContent(binIter+1)  << std::endl;
	      else std::cout << "     " << binIter << ", " << tempLowest << std::endl;
	    }

	    if(leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && !isFound){
	      nextMaxBin = binIter+1;
	      isFound = true;
	    }
	  }

	  tempLowBin = nextMaxBin-1;
	}
      }

      
      isDrawn = false;
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){
	  leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2]->SetMaximum(tempMax);
	  
	  if(!isDrawn) leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST");
	  else leadingPtBoundSpect_EffLeadCut100_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST SAME");
	  isDrawn = true;
	}
      }
      spectLeg_p[iter]->Draw("SAME");
      line_p->DrawLine(leadingPtBoundSpect_EffLeadCut100_p[iter][0]->GetXaxis()->GetXmin(), tempLowest, leadingPtBoundSpect_EffLeadCut100_p[iter][0]->GetXaxis()->GetXmax(), tempLowest);
      std::cout << trigType[iter] << ", " << tempLowest << std::endl;
      tempLowest = 10000000000;
      foundLowBin = false;
      tempLowBin = 1000000;




      if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())){
	ptTrkSpectBoundCanv_p[iter] = new TCanvas(canvTrkSpectName[iter].c_str(), canvTrkSpectName[iter].c_str(), 700, 700);
	ptTrkSpectBoundCanv_p[iter]->SetTopMargin(0.01);
	ptTrkSpectBoundCanv_p[iter]->SetRightMargin(0.01);
	
	ptTrkSpectBoundCanv_p[iter]->cd();
	gPad->SetLogy();
	gPad->SetLogx();
	
	tempMax = -1;
	
	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Not") != 0){
	    
	    if(iter2 != 0) trkPtBoundSpect_p[iter][iter2]->Add(trkPtBoundSpect_p[iter][iter2-1]);
	    
	    if(trkPtBoundSpect_p[iter][iter2]->GetMaximum() > tempMax) tempMax = trkPtBoundSpect_p[iter][iter2]->GetMaximum();
	    
	    trkPtBoundSpect_p[iter][iter2]->SetMinimum(1);
	    trkPtBoundSpect_p[iter][iter2]->GetXaxis()->CenterTitle();
	    trkPtBoundSpect_p[iter][iter2]->GetYaxis()->CenterTitle();
	    trkPtBoundSpect_p[iter][iter2]->GetYaxis()->SetTitleOffset(trkPtBoundSpect_p[iter][iter2]->GetYaxis()->GetTitleOffset()+.20);
	    trkPtBoundSpect_p[iter][iter2]->GetXaxis()->SetTitleOffset(trkPtBoundSpect_p[iter][iter2]->GetXaxis()->GetTitleOffset()+.20);
	    
	    if(iter2 == 0){
	      trkPtBoundSpect_p[iter][iter2]->SetMarkerColor(17);
	      trkPtBoundSpect_p[iter][iter2]->SetLineColor(1);
	      trkPtBoundSpect_p[iter][iter2]->SetFillColor(17);
	    }
	    else{
	      trkPtBoundSpect_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);
	      trkPtBoundSpect_p[iter][iter2]->SetLineColor(1);
	      trkPtBoundSpect_p[iter][iter2]->SetFillColor(trigColors[iter2]);
	    }
	  }
	}
	
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){

	  Int_t nextMaxBin = -1;
	  Bool_t isFound = false;
	  if(!foundLowBin){
	    tempLowBin = trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetNbinsX();
	    foundLowBin = true;
	  }

	  std::cout << "  " << trigType[iter] << ", " << tempLowBin << std::endl;

	  for(Int_t binIter = 0; binIter < tempLowBin; binIter++){
	    if(trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < tempLowest && trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0){
	      Bool_t contBool = false;

	      for(Int_t iter3 = iter2+1; iter3 < trigTypeCount[iter]-1; iter3++){
		if(trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter3]->GetBinContent(binIter+1)){
		  contBool = true;
		  break;
		}
	      }

	      if(contBool) continue;
	      tempLowest = trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1);

	    
	      if(iter2 != trigTypeCount[iter]-1) std::cout << "     " << binIter << ", " << tempLowest << ", " << trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2-1]->GetBinContent(binIter+1)  << std::endl;
	      else std::cout << "     " << binIter << ", " << tempLowest << std::endl;
	    }

	    if(trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && !isFound){
	      nextMaxBin = binIter+1;
	      isFound = true;
	    }
	  }

	  tempLowBin = nextMaxBin-1;
	}
      }


	Bool_t isDrawn = false;
	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){
	    trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->SetMaximum(tempMax);
	    
	    if(!isDrawn) trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST");
	    else trkPtBoundSpect_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST SAME");
	    
	    isDrawn = true;
	  }
	}
	spectLeg_p[iter]->Draw("SAME");
	line_p->DrawLine(trkPtBoundSpect_p[iter][0]->GetXaxis()->GetXmin(), tempLowest, trkPtBoundSpect_p[iter][0]->GetXaxis()->GetXmax(), tempLowest);
	std::cout << trigType[iter] << ", " << tempLowest << std::endl;
	tempLowest = 10000000000;
	foundLowBin = false;
	tempLowBin = 1000000;

	

	ptTrkSpectPrsclBoundCanv_p[iter] = new TCanvas(canvTrkSpectPrsclName[iter].c_str(), canvTrkSpectPrsclName[iter].c_str(), 700, 700);
	ptTrkSpectPrsclBoundCanv_p[iter]->SetTopMargin(0.01);
	ptTrkSpectPrsclBoundCanv_p[iter]->SetRightMargin(0.01);
	
	ptTrkSpectPrsclBoundCanv_p[iter]->cd();
	gPad->SetLogy();
	gPad->SetLogx();
	
	tempMax = -1;
	
	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Not") != 0){
	    
	    if(iter2 != 0) trkPtBoundSpectPrscl_p[iter][iter2]->Add(trkPtBoundSpectPrscl_p[iter][iter2-1]);
	    
	    if(trkPtBoundSpectPrscl_p[iter][iter2]->GetMaximum() > tempMax) tempMax = trkPtBoundSpectPrscl_p[iter][iter2]->GetMaximum();
	    
	    trkPtBoundSpectPrscl_p[iter][iter2]->SetMinimum(1);
	    trkPtBoundSpectPrscl_p[iter][iter2]->GetXaxis()->CenterTitle();
	    trkPtBoundSpectPrscl_p[iter][iter2]->GetYaxis()->CenterTitle();
	    trkPtBoundSpectPrscl_p[iter][iter2]->GetYaxis()->SetTitleOffset(trkPtBoundSpectPrscl_p[iter][iter2]->GetYaxis()->GetTitleOffset()+.20);
	    trkPtBoundSpectPrscl_p[iter][iter2]->GetXaxis()->SetTitleOffset(trkPtBoundSpectPrscl_p[iter][iter2]->GetXaxis()->GetTitleOffset()+.20);
	    
	    if(iter2 == 0){
	      trkPtBoundSpectPrscl_p[iter][iter2]->SetMarkerColor(17);
	      trkPtBoundSpectPrscl_p[iter][iter2]->SetLineColor(1);
	      trkPtBoundSpectPrscl_p[iter][iter2]->SetFillColor(17);
	    }
	    else{
	      trkPtBoundSpectPrscl_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);
	      trkPtBoundSpectPrscl_p[iter][iter2]->SetLineColor(1);
	      trkPtBoundSpectPrscl_p[iter][iter2]->SetFillColor(trigColors[iter2]);
	    }
	  }
	}
	
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){

	  Int_t nextMaxBin = -1;
	  Bool_t isFound = false;
	  if(!foundLowBin){
	    tempLowBin = trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2]->GetNbinsX();
	    foundLowBin = true;
	  }

	  std::cout << "  " << trigType[iter] << ", " << tempLowBin << std::endl;

	  for(Int_t binIter = 0; binIter < tempLowBin; binIter++){
	    if(trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < tempLowest && trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0){
	      Bool_t contBool = false;

	      for(Int_t iter3 = iter2+1; iter3 < trigTypeCount[iter]-1; iter3++){
		if(trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) < trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter3]->GetBinContent(binIter+1)){
		  contBool = true;
		  break;
		}
	      }

	      if(contBool) continue;
	      tempLowest = trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1);

	    
	      if(iter2 != trigTypeCount[iter]-1) std::cout << "     " << binIter << ", " << tempLowest << ", " << trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2-1]->GetBinContent(binIter+1)  << std::endl;
	      else std::cout << "     " << binIter << ", " << tempLowest << std::endl;
	    }

	    if(trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2]->GetBinContent(binIter+1) != 0 && !isFound){
	      nextMaxBin = binIter+1;
	      isFound = true;
	    }
	  }

	  tempLowBin = nextMaxBin-1;
	}
      }


	isDrawn = false;
	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){
	    trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2]->SetMaximum(tempMax);
	    
	    if(!isDrawn) trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST");
	    else trkPtBoundSpectPrscl_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 HIST SAME");
	    
	    isDrawn = true;
	  }
	}
	spectLeg_p[iter]->Draw("SAME");
	line_p->DrawLine(trkPtBoundSpectPrscl_p[iter][0]->GetXaxis()->GetXmin(), tempLowest, trkPtBoundSpectPrscl_p[iter][0]->GetXaxis()->GetXmax(), tempLowest);

	ptTrkSpectRelErrBoundCanv_p[iter] = new TCanvas(canvTrkSpectRelErrName[iter].c_str(), canvTrkSpectRelErrName[iter].c_str(), 700, 700);
	ptTrkSpectRelErrBoundCanv_p[iter]->SetTopMargin(0.01);
	ptTrkSpectRelErrBoundCanv_p[iter]->SetRightMargin(0.01);
	
	ptTrkSpectRelErrBoundCanv_p[iter]->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	
	tempMax = -1;
	Double_t tempMin = 10000;
	
	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(strcmp(trigSpectIsOpt[iter][iter2].c_str(), "Not") != 0){
	    
	    if(trkPtBoundSpectRelErr_p[iter][iter2]->GetMaximum() > tempMax) tempMax = trkPtBoundSpectRelErr_p[iter][iter2]->GetMaximum();

	    for(Int_t binIter = 0; binIter < trkPtBoundSpectRelErr_p[iter][iter2]->GetNbinsX(); binIter++){
	      if(trkPtBoundSpectRelErr_p[iter][iter2]->GetBinContent(binIter+1) != 0 && trkPtBoundSpectRelErr_p[iter][iter2]->GetBinContent(binIter+1) < tempMin) tempMin = trkPtBoundSpectRelErr_p[iter][iter2]->GetBinContent(binIter+1);
	    }
	    
	    trkPtBoundSpectRelErr_p[iter][iter2]->GetXaxis()->CenterTitle();
	    trkPtBoundSpectRelErr_p[iter][iter2]->GetYaxis()->CenterTitle();
	    trkPtBoundSpectRelErr_p[iter][iter2]->GetYaxis()->SetTitleOffset(trkPtBoundSpectRelErr_p[iter][iter2]->GetYaxis()->GetTitleOffset()+.20);
	    trkPtBoundSpectRelErr_p[iter][iter2]->GetXaxis()->SetTitleOffset(trkPtBoundSpectRelErr_p[iter][iter2]->GetXaxis()->GetTitleOffset()+.20);
	    
	    trkPtBoundSpectRelErr_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);
	    trkPtBoundSpectRelErr_p[iter][iter2]->SetLineColor(1);
	    trkPtBoundSpectRelErr_p[iter][iter2]->SetFillColor(trigColors[iter2]);
	    trkPtBoundSpectRelErr_p[iter][iter2]->SetMarkerSize(1);
	    trkPtBoundSpectRelErr_p[iter][iter2]->SetMarkerStyle(20);

	    trkPtBoundSpectRelErr_p[iter][iter2]->GetYaxis()->SetTitle("Relative Error");

	    spectLeg2_p[iter]->AddEntry(trkPtBoundSpectRelErr_p[iter][iter2], Form("%s %d", trigType[iter].c_str(), trigThresh[iter][iter2]), "P");
	  }
	}

	std::cout << trigType[iter] << " tempMin " << tempMin << std::endl;

	isDrawn = false;
	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(strcmp(trigSpectIsOpt[iter][trigTypeCount[iter]-1-iter2].c_str(), "Not") != 0){
	    trkPtBoundSpectRelErr_p[iter][trigTypeCount[iter]-1-iter2]->SetMaximum(tempMax*2);
	    trkPtBoundSpectRelErr_p[iter][trigTypeCount[iter]-1-iter2]->SetMinimum(tempMin/2);


	    trkPtBoundSpectRelErr_p[iter][trigTypeCount[iter]-1-iter2]->Print("ALL") ;

	    if(!isDrawn) trkPtBoundSpectRelErr_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 P");
	    else trkPtBoundSpectRelErr_p[iter][trigTypeCount[iter]-1-iter2]->DrawCopy("E1 P SAME");
	    
	    isDrawn = true;
	  }
	}
	spectLeg2_p[iter]->Draw("SAME");
      }
    }
  }

  TLegend* maxLeg_p[nTrigType];

  for(Int_t iter = 0; iter < nTrigType; iter++){
    canvMaxMissedName[iter] =  Form("maxMissedPt_%s_c", trigType[iter].c_str());

    maxMissedPtCanv_p[iter] = new TCanvas(canvMaxMissedName[iter].c_str(), canvMaxMissedName[iter].c_str(), 700, 700);
    maxMissedPtCanv_p[iter]->SetTopMargin(0.01);
    maxMissedPtCanv_p[iter]->SetRightMargin(0.01);

    maxMissedPtCanv_p[iter]->cd();

    maxMissedPt100_p[iter]->GetXaxis()->CenterTitle();
    maxMissedPt100_p[iter]->GetYaxis()->CenterTitle();
    maxMissedPt100_p[iter]->GetYaxis()->SetTitleOffset(maxMissedPt100_p[iter]->GetYaxis()->GetTitleOffset()+.10);
    maxMissedPt100_p[iter]->SetMarkerColor(kRed);
    maxMissedPt100_p[iter]->SetLineColor(kRed);
    maxMissedPt100_p[iter]->SetTitle("");
    maxMissedPt100_p[iter]->GetXaxis()->SetTitle("HLT p_{T}");
    maxMissedPt100_p[iter]->GetYaxis()->SetTitle("Max Missed p_{T}");

    maxMissedPt100_p[iter]->DrawCopy("E1 HIST");

    maxMissedPt99_p[iter]->DrawCopy("E1 HIST SAME");


    maxLeg_p[iter] = new TLegend(0.17, 0.60, 0.58, 0.80);
    maxLeg_p[iter]->SetFillColor(0);
    maxLeg_p[iter]->SetFillStyle(0);
    maxLeg_p[iter]->SetTextFont(43);
    maxLeg_p[iter]->SetTextSize(16);
    maxLeg_p[iter]->SetBorderSize(0);

    maxLeg_p[iter]->AddEntry(maxMissedPt99_p[iter], "99%", "L P");
    maxLeg_p[iter]->AddEntry(maxMissedPt100_p[iter], "100%", "L P");
    maxLeg_p[iter]->Draw("SAME");
  }


  for(Int_t iter = 0; iter < nTrigType; iter++){
    canvR9Name[iter] =  Form("r9_%s_c", trigType[iter].c_str());

    R9Canv_p[iter] = new TCanvas(canvR9Name[iter].c_str(), canvR9Name[iter].c_str(), 700, 700);
    R9Canv_p[iter]->SetTopMargin(0.01);
    R9Canv_p[iter]->SetRightMargin(0.01);

    R9Canv_p[iter]->cd();

    gPad->SetLogy();

    Float_t maxR9 = e1R9_p[iter]->GetMaximum();
    if(e2R9_p[iter]->GetMaximum() > maxR9) maxR9 = e2R9_p[iter]->GetMaximum();

    e1R9_p[iter]->SetMaximum(maxR9+TMath::Sqrt(maxR9));
    e2R9_p[iter]->SetMaximum(maxR9+TMath::Sqrt(maxR9));

    e1R9_p[iter]->GetXaxis()->CenterTitle();
    e1R9_p[iter]->GetYaxis()->CenterTitle();
    e1R9_p[iter]->GetYaxis()->SetTitleOffset(e1R9_p[iter]->GetYaxis()->GetTitleOffset()+.10);
    e1R9_p[iter]->SetMarkerColor(kRed);
    e1R9_p[iter]->SetLineColor(kRed);
    e1R9_p[iter]->DrawCopy("E1 HIST");

    e2R9_p[iter]->SetMarkerColor(kBlue);
    e2R9_p[iter]->SetLineColor(kBlue);
    e2R9_p[iter]->DrawCopy("E1 HIST SAME");
  }

  std::string outName = inHistFile;
  const std::string inString = "HIST";
  const std::string outString = "PLOT";
  std::size_t strIndex = 0;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString);
  }


  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");

  for(Int_t iter = 0; iter < nTrigType; iter++){
    //    if(strcmp("GAMMA", trigType[iter].substr(0,5).c_str()) != 0) continue;

    trigCanvPt_p[iter]->RedrawAxis();
    trigCanvPt_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(trigCanvPt_p[iter], Form("pdfDir/%s", canvPtName[iter].c_str()), "pdf");

    //    ptSpectBoundCanv_p[iter]->RedrawAxis();
    //    ptSpectBoundCanv_p[iter]->Write("", TObject::kOverwrite);
    //    claverCanvasSaving(ptSpectBoundCanv_p[iter], Form("pdfDir/%s", canvSpectName[iter].c_str()), "pdf");

    if(isOpt){
      ptSpectBoundCanv_Eff_p[iter]->RedrawAxis();
      ptSpectBoundCanv_Eff_p[iter]->Write("", TObject::kOverwrite);
      claverCanvasSaving(ptSpectBoundCanv_Eff_p[iter], Form("pdfDir/%s", canvSpectNameEff[iter].c_str()), "pdf");
      
      ptSpectBoundCanv_EffLeadCut99_p[iter]->RedrawAxis();
      ptSpectBoundCanv_EffLeadCut99_p[iter]->Write("", TObject::kOverwrite);
      claverCanvasSaving(ptSpectBoundCanv_EffLeadCut99_p[iter], Form("pdfDir/%s", canvSpectNameEffLeadCut99[iter].c_str()), "pdf");

      ptSpectBoundCanv_EffLeadCut100_p[iter]->RedrawAxis();
      claverCanvasSaving(ptSpectBoundCanv_EffLeadCut100_p[iter], Form("pdfDir/%s", canvSpectNameEffLeadCut100[iter].c_str()), "pdf");

      if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())){
	ptTrkSpectBoundCanv_p[iter]->RedrawAxis();
	ptTrkSpectBoundCanv_p[iter]->Write("", TObject::kOverwrite);
	claverCanvasSaving(ptTrkSpectBoundCanv_p[iter], Form("pdfDir/%s", canvTrkSpectName[iter].c_str()), "pdf");
	ptTrkSpectPrsclBoundCanv_p[iter]->RedrawAxis();
	ptTrkSpectPrsclBoundCanv_p[iter]->Write("", TObject::kOverwrite);
	claverCanvasSaving(ptTrkSpectPrsclBoundCanv_p[iter], Form("pdfDir/%s", canvTrkSpectPrsclName[iter].c_str()), "pdf");

	ptTrkSpectRelErrBoundCanv_p[iter]->RedrawAxis();
	ptTrkSpectRelErrBoundCanv_p[iter]->Write("", TObject::kOverwrite);
	claverCanvasSaving(ptTrkSpectRelErrBoundCanv_p[iter], Form("pdfDir/%s", canvTrkSpectRelErrName[iter].c_str()), "pdf");
      }      
    }

    maxMissedPtCanv_p[iter]->RedrawAxis();
    maxMissedPtCanv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(maxMissedPtCanv_p[iter], Form("pdfDir/%s", canvMaxMissedName[iter].c_str()), "pdf");

    R9Canv_p[iter]->RedrawAxis();
    R9Canv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(R9Canv_p[iter], Form("pdfDir/%s", canvR9Name[iter].c_str()), "pdf");

  }

  for(Int_t iter = 0; iter < nTrigPlotType; iter++){
    trigCanvRate_p[iter]->RedrawAxis();
    trigCanvRate_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(trigCanvRate_p[iter], Form("pdfDir/%s", canvRateName[iter].c_str()), "pdf");
  }

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    delete trigCanvPt_p[iter];
    delete ptSpectBoundCanv_p[iter];
    if(isOpt){
      delete ptSpectBoundCanv_Eff_p[iter];
      delete ptSpectBoundCanv_EffLeadCut99_p[iter];
      delete ptSpectBoundCanv_EffLeadCut100_p[iter];
      if(!strcmp("TRACK", trigType[iter].substr(0,5).c_str())){
	delete ptTrkSpectBoundCanv_p[iter];
	delete ptTrkSpectPrsclBoundCanv_p[iter];
	delete ptTrkSpectRelErrBoundCanv_p[iter];
      }
    }
    delete maxMissedPtCanv_p[iter];
  }

  for(Int_t iter = 0; iter < nTrigPlotType; iter++){
    delete trigCanvRate_p[iter];
  }

  inFile_p->Close();
  delete inFile_p;

  if(isOpt){
    inOptFile_p->Close();
    delete inOptFile_p;
  }

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 3 && argc != 4){
    std::cout << "Usage: plotTrigTurnOn_HI <inHistFile> <inTrigFileName> <inOptionalEffFile>" << std::endl;
    std::cout << "argNum: " << argc << std::endl;
    for(Int_t iter = 0; iter < argc; iter++){
      std::cout << "arg " << iter << ": " << argv[iter] << std::endl;
    }

    return -1;
  }

  int rStatus = -1;

  if(argc == 3) rStatus = plotTrigTurnOn_HI(argv[1], argv[2]);
  else if(argc == 4) rStatus = plotTrigTurnOn_HI(argv[1], argv[2], argv[3]);

  return rStatus;
}
