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
      if(std::string::npos== buffer.find("HLT")) nTrigTypeTemp++;
      listOfTrig.push_back(buffer);
      nLines++;
    }
  }
  delete inTrigFile;
  std::cout << "Trigger List Loaded" << std::endl;

  const Int_t nTrigType = nTrigTypeTemp;
  std::string trigType[nTrigType];
  Int_t trigTypeCount[nTrigType];

  nLines = 0;

  for(Int_t iter = 0; iter < (Int_t)listOfTrig.size(); iter++){
    if(std::string::npos == listOfTrig[iter].find("HLT")){
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
  for(Int_t iter = 0; iter < nTrigType; iter++){
    if(trigTypeCount[iter] > maxNTrig) maxNTrig = trigTypeCount[iter];
  }
  const Int_t maxNTrig2 = maxNTrig;

  std::string trigName[nTrigType][maxNTrig2];
  Int_t trigThresh[nTrigType][maxNTrig2];

  nLines = 0;
  Int_t tempPosIter = 0;

  for(Int_t iter = 1; iter < (Int_t)(listOfTrig.size()); iter++){
    if(std::string::npos== listOfTrig[iter].find("HLT")){
      nLines++;
      tempPosIter = 0;
    }
    else{
      trigName[nLines][tempPosIter] = listOfTrig[iter];
      trigThresh[nLines][tempPosIter] = listOfThresh[iter];
      tempPosIter++;
    }
  }

  TFile* inFile_p = new TFile(inHistFile.c_str(), "READ");

  std::string canvPtName[nTrigType];
  TCanvas* trigCanvPt_p[nTrigType];
  TGraphAsymmErrors* aPt_p[nTrigType][maxNTrig2];

  std::string canvRateName[nTrigType];
  TCanvas* trigCanvRate_p[nTrigType];
  TH1F* ratesUnmatched_p[nTrigType];

  //ONLY EDITING SHOULD OCCUR HERE
  const Int_t maxPlotTrig = 5;
  const Int_t trigColors[5] = {1, kBlue, kRed, kYellow+1, kMagenta};
  TLegend* trigLeg_p[nTrigType];

  for(Int_t iter = 0; iter < nTrigType; iter++){      
    canvPtName[iter] = Form("%s_%d_%d_pt_c", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter] - 1]);
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
      trigLeg_p[iter]->AddEntry(aPt_p[iter][iter2], Form("%s, %d", trigType[iter].c_str(), trigThresh[iter][iter2]), "P L");
    }

    TLine* oneLine_p = new TLine(aPt_p[iter][0]->GetXaxis()->GetXmin(), 1, aPt_p[iter][0]->GetXaxis()->GetXmax(), 1);
    oneLine_p->SetLineStyle(2);
    oneLine_p->DrawClone();

    delete oneLine_p;
  }

  for(Int_t iter = 0; iter < nTrigType; iter++){
    canvRateName[iter] =  Form("%s_%d_%d_Rate_c", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter] - 1]);

    trigCanvRate_p[iter] = new TCanvas(canvRateName[iter].c_str(), canvRateName[iter].c_str(), 700, 700);
    trigCanvRate_p[iter]->cd();

    ratesUnmatched_p[iter] = (TH1F*)inFile_p->Get(Form("ratesUnmatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]));


    TH1F* hEmpty = new TH1F("hEmpty", Form(";p_{T,%s}^{reco};Rate (Hz)", trigType[iter].c_str()), ratesUnmatched_p[iter]->GetNbinsX(), ratesUnmatched_p[iter]->GetXaxis()->GetXmin(), ratesUnmatched_p[iter]->GetXaxis()->GetXmax());
    hEmpty->GetXaxis()->CenterTitle();
    hEmpty->GetYaxis()->CenterTitle();
    hEmpty->SetMaximum(ratesUnmatched_p[iter]->GetMaximum() + 10*TMath::Sqrt(ratesUnmatched_p[iter]->GetMaximum()));
    if(hEmpty->GetMaximum() < 100) hEmpty->SetMaximum(100);
    hEmpty->SetMinimum(0.5);
    hEmpty->DrawCopy();
    delete hEmpty;

    ratesUnmatched_p[iter]->DrawCopy("E1 SAME");
    gPad->SetLogy();

    TLine* tenLine_p = new TLine(ratesUnmatched_p[iter]->GetXaxis()->GetXmin(), 10, ratesUnmatched_p[iter]->GetXaxis()->GetXmax(), 10);
    tenLine_p->SetLineStyle(2);
    tenLine_p->DrawClone();

    tenLine_p->DrawLine(ratesUnmatched_p[iter]->GetXaxis()->GetXmin(), 1, ratesUnmatched_p[iter]->GetXaxis()->GetXmax(), 1);

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
    trigCanvRate_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(trigCanvPt_p[iter], Form("pdfDir/%s", canvPtName[iter].c_str()), "pdf");
    claverCanvasSaving(trigCanvRate_p[iter], Form("pdfDir/%s", canvRateName[iter].c_str()), "pdf");
  }

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    delete trigCanvPt_p[iter];
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
