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
const Int_t maxPlotTrig = 8;
const Int_t trigColors[8] = {1, kBlue, kRed, kYellow+1, kMagenta, kGreen+3, kCyan+2, kGray+1};

void claverCanvasSaving(TCanvas* c, TString s,TString format="gif"){
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
  return;
}

int plotTrigTurnOn_HI(const std::string inHistFile, const std::string inTrigFileName)
{
  gStyle->SetOptStat(0);

  std::string buffer;
  std::vector<std::string> listOfTrig;
  std::vector<Int_t> listOfThresh;
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
    }
    else{
      trigTypeCount[nLines-1]++;

      std::size_t strIndex = 0;
      while(true){
        strIndex = listOfTrig[iter].find(",");
        if(strIndex == std::string::npos) break;
        listOfThresh.push_back(std::stoi(listOfTrig[iter].substr(strIndex+1)));
        listOfTrig[iter].replace(strIndex, std::string::npos, "");
      }
    }
  }

  Int_t maxNTrig = -1;
  Int_t maxNum = 0;
  for(Int_t iter = 0; iter < nTrigType; iter++){
    if(trigTypeCount[iter] > maxNTrig) maxNTrig = trigTypeCount[iter];
  }

  std::cout << "A" << std::endl;
  
  const Int_t maxNTrig2 = maxNTrig;

  std::string trigName[nTrigType][maxNTrig2];
  Int_t trigThresh[nTrigType][maxNTrig2];

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

  std::cout << "B" << std::endl;

  TFile* inFile_p = new TFile(inHistFile.c_str(), "READ");

  std::string canvPtName[nTrigType];
  TCanvas* trigCanvPt_p[nTrigType];
  TGraphAsymmErrors* aPt_p[nTrigType][maxNTrig2];

  const Int_t nTrigPlotType = nTrigPlotTypeTemp;
  std::string canvRateName[nTrigPlotType];
  TCanvas* trigCanvRate_p[nTrigPlotType];
  TH1F* hEmptyRate[nTrigPlotType];
  Float_t maxXRate[nTrigPlotType];
  Float_t minXRate[nTrigPlotType];
  Float_t maxYRate[nTrigPlotType];
  Float_t minYRate[nTrigPlotType];
  Float_t binWidth[nTrigPlotType];
  TH1F* ratesUnmatched_p[nTrigType];

  for(Int_t iter = 0; iter < nTrigPlotType; iter++){
    maxXRate[iter] = -1;
    minXRate[iter] = 100000;
    maxYRate[iter] = -1;
    minYRate[iter] = 10000000;
    binWidth[iter] = 10000000;
  }

  TLegend* trigLeg_p[nTrigType];
  TLegend* rateLeg_p[nTrigPlotType];

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

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      if(iter2 == maxPlotTrig) break;
      aPt_p[iter][iter2] = (TGraphAsymmErrors*)inFile_p->Get(Form("%s_%s_pt_asymm", trigName[iter][iter2].c_str(), trigType[iter].c_str()));
      std::cout << Form("%s_%s_pt_asymm", trigName[iter][iter2].c_str(), trigType[iter].c_str()) << std::endl;
      aPt_p[iter][iter2]->SetLineColor(trigColors[iter2]);
      aPt_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);

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

    rateLeg_p[iter] = new TLegend(0.57, 0.73, 0.98, 0.98);
    rateLeg_p[iter]->SetFillColor(0);
    rateLeg_p[iter]->SetTextFont(43);
    rateLeg_p[iter]->SetTextSize(16);
  }

  for(Int_t iter = 0; iter < nTrigType; iter++){
    ratesUnmatched_p[iter] = (TH1F*)inFile_p->Get(Form("ratesUnmatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]));

    std::cout << Form("ratesUnmatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]) << std::endl;

    if(maxXRate[trigPlotPos[iter]] < ratesUnmatched_p[iter]->GetXaxis()->GetXmax()) maxXRate[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetXaxis()->GetXmax();
    if(minXRate[trigPlotPos[iter]] > ratesUnmatched_p[iter]->GetXaxis()->GetXmin()) minXRate[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetXaxis()->GetXmin();

    if(maxYRate[trigPlotPos[iter]] < ratesUnmatched_p[iter]->GetMaximum()) maxYRate[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetMaximum();
    if(minYRate[trigPlotPos[iter]] > ratesUnmatched_p[iter]->GetMinimum()) minYRate[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetMinimum();

    if(binWidth[trigPlotPos[iter]] > ratesUnmatched_p[iter]->GetXaxis()->GetBinWidth(1)/2.0) binWidth[trigPlotPos[iter]] = ratesUnmatched_p[iter]->GetXaxis()->GetBinWidth(1)/2.0;
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
    trigCanvPt_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(trigCanvPt_p[iter], Form("pdfDir/%s", canvPtName[iter].c_str()), "pdf");
  }

  for(Int_t iter = 0; iter < nTrigPlotType; iter++){
    trigCanvRate_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(trigCanvRate_p[iter], Form("pdfDir/%s", canvRateName[iter].c_str()), "pdf");
  }

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    delete trigCanvPt_p[iter];
  }

  for(Int_t iter = 0; iter < nTrigPlotType; iter++){
    delete trigCanvRate_p[iter];
  }

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 3){
    std::cout << "Usage: plotTrigTurnOn_HI <inHistFile> <inTrigFileName>" << std::endl;
    std::cout << "argNum: " << argc << std::endl;
    for(Int_t iter = 0; iter < argc; iter++){
      std::cout << "arg " << iter << ": " << argv[iter] << std::endl;
    }

    return -1;
  }

  int rStatus = -1;

  rStatus = plotTrigTurnOn_HI(argv[1], argv[2]);

  return rStatus;
}
