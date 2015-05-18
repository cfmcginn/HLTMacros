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
      if(std::string::npos== buffer.find("HLT")) nTrigTypeTemp++;
      listOfTrig.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "Trigger List Loaded" << std::endl;

  const Int_t nTrigType = nTrigTypeTemp;
  std::string trigType[nTrigType];
  Int_t trigTypeCount[nTrigType];

  delete inTrigFile;
  inTrigFile = new ifstream(inTrigFileName.data());
  nLines = 0;

  while(true){
    *inTrigFile >> buffer;
    if(inTrigFile->eof()) break;
    if(std::string::npos== buffer.find("HLT")){
      trigType[nLines] = buffer;
      nLines++;
    }
  }

  delete inTrigFile;
  nLines = 0;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    trigTypeCount[iter] = 0;
  }

  for(Int_t iter = 1; iter < (Int_t)(listOfTrig.size()); iter++){
    std::cout << listOfTrig[iter] << std::endl;

    if(std::string::npos== listOfTrig[iter].find("HLT")) nLines++;
    else trigTypeCount[nLines]++;
  }

  Int_t maxNTrig = -1;
  for(Int_t iter = 0; iter < nTrigType; iter++){
    if(trigTypeCount[iter] > maxNTrig) maxNTrig = trigTypeCount[iter];
  }
  const Int_t maxNTrig2 = maxNTrig;

  std::string trigName[nTrigType][maxNTrig2];

  nLines = 0;
  Int_t tempPosIter = 0;

  for(Int_t iter = 1; iter < (Int_t)(listOfTrig.size()); iter++){
    if(std::string::npos== listOfTrig[iter].find("HLT")){
      nLines++;
      tempPosIter = 0;
    }
    else{
      trigName[nLines][tempPosIter] = listOfTrig[iter];
      tempPosIter++;
    }
  }

  TFile* inFile_p = new TFile(inHistFile.c_str(), "READ");

  std::string canvPtName[nTrigType];
  TCanvas* trigCanvPt_p[nTrigType];
  TGraphAsymmErrors* aPt_p[nTrigType][maxNTrig2];
  //ONLY EDITING SHOULD OCCUR HERE
  const Int_t maxPlotTrig = 5;
  const Int_t trigColors[5] = {1, kBlue, kRed, kYellow+1, kMagenta};
  TLegend* trigLeg_p[nTrigType];

  for(Int_t iter = 0; iter < nTrigType; iter++){      
    canvPtName[iter] = Form("%s_%s_pt_c", trigType[iter].c_str(), trigName[iter][0].c_str());
    trigCanvPt_p[iter] = new TCanvas(canvPtName[iter].c_str(), canvPtName[iter].c_str(), 700, 700);
    trigCanvPt_p[iter]->cd();
    trigLeg_p[iter] = new TLegend(0.65, 0.20, 0.98, 0.45);
    trigLeg_p[iter]->SetFillColor(0);
    trigLeg_p[iter]->SetTextFont(43);
    trigLeg_p[iter]->SetTextSize(16);

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      if(iter2 == maxPlotTrig) break;
      aPt_p[iter][iter2] = (TGraphAsymmErrors*)inFile_p->Get(Form("%s_%s_pt_asymm", trigName[iter][iter2].c_str(), trigType[iter].c_str()));
      aPt_p[iter][iter2]->SetLineColor(trigColors[iter2]);
      aPt_p[iter][iter2]->SetMarkerColor(trigColors[iter2]);

      if(iter2 == 0){
	TH1F* hEmpty = new TH1F("hEmpty", ";p_{T}^{reco};Efficiency", aPt_p[iter][iter2]->GetN(), aPt_p[iter][iter2]->GetXaxis()->GetXmin(), aPt_p[iter][iter2]->GetXaxis()->GetXmax());
	hEmpty->SetMaximum(1.1);
	hEmpty->SetMinimum(0.0);
	hEmpty->DrawCopy();
	trigLeg_p[iter]->Draw("SAME");
	delete hEmpty;
      }

      aPt_p[iter][iter2]->Draw("P E");
      trigLeg_p[iter]->AddEntry(aPt_p[iter][iter2], trigName[iter][iter2].c_str(), "P L");
    }

    TLine* oneLine_p = new TLine(aPt_p[iter][0]->GetXaxis()->GetXmin(), 1, aPt_p[iter][0]->GetXaxis()->GetXmax(), 1);
    oneLine_p->SetLineStyle(2);
    oneLine_p->DrawClone();

    delete oneLine_p;
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

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    delete trigCanvPt_p[iter];
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
