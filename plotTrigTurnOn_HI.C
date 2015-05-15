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

int plotTrigTurnOn_HI(const std::string inHistFile, const std::string inTrigFileName, const std::string outFile)
{
  std::string buffer;
  std::vector<std::string> listOfTrig;
  int nLines = 0;
  ifstream inTrigFile(inTrigFileName.data());

  std::cout << inTrigFileName << std::endl;
  std::cout << inTrigFile.is_open() << std::endl;

  if(!inTrigFile.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return 1;
  }
  else{
    while(true){
      inTrigFile >> buffer;
      if(inTrigFile.eof()) break;
      listOfTrig.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "Trigger List Loaded" << std::endl;

  const Int_t nTrigType = 4;
  const std::string trigType[4] = {"TRK", "3JET", "4JET", "GAMMA"};
  Int_t trigTypeCount[nTrigType] = {0, 0, 0, 0};
  Bool_t trigTypeBool[nTrigType] = {false, false, false, false};

  for(Int_t iter = 0; iter < (Int_t)(listOfTrig.size()); iter++){
    std::cout << listOfTrig[iter] << std::endl;

    Bool_t typeSwitch = false;

    for(Int_t typeIter = 0; typeIter < nTrigType; typeIter++){
      if(!strcmp(listOfTrig[iter].c_str(), trigType[typeIter].c_str())){
	
	for(Int_t typeIter2 = 0; typeIter2 < nTrigType; typeIter2++){
	  trigTypeBool[typeIter2] = false;
	}

	trigTypeBool[typeIter] = true;
	typeSwitch = true;
	break;
      }
    }

    if(typeSwitch) continue;

    for(Int_t typeIter = 0; typeIter < nTrigType; typeIter++){
      if(trigTypeBool[typeIter]){
	trigTypeCount[typeIter]++;
	break;
      }
    }
  }

  Int_t maxNTrig = -1;
  for(Int_t iter = 0; iter < nTrigType; iter++){
    if(trigTypeCount[iter] > maxNTrig) maxNTrig = trigTypeCount[iter];
  }
  const Int_t maxNTrig2 = maxNTrig;

  std::string trigName[nTrigType][maxNTrig2];

  Int_t tempPosIter = 0;

  for(Int_t iter = 0; iter < (Int_t)(listOfTrig.size()); iter++){
    Bool_t typeSwitch = false;
    for(Int_t typeIter = 0; typeIter < nTrigType; typeIter++){
      if(!strcmp(listOfTrig[iter].c_str(), trigType[typeIter].c_str())){
        for(Int_t typeIter2 = 0; typeIter2 < nTrigType; typeIter2++){
          trigTypeBool[typeIter2] = false;
        }

        trigTypeBool[typeIter] = true;
	typeSwitch = true;
	tempPosIter = 0;
        break;
      }
    }

    if(typeSwitch) continue;

    for(Int_t typeIter = 0; typeIter < nTrigType; typeIter++){
      if(trigTypeBool[typeIter]){
        trigName[typeIter][tempPosIter] = listOfTrig[iter];
	tempPosIter++;
        break;
      }
    }
  }

  Int_t error = 0;
  std::cout << error << std::endl;
  error++;

  TFile* inFile_p = new TFile(inHistFile.c_str(), "READ");

  std::string canvPtName[nTrigType];
  TCanvas* trigCanvPt_p[nTrigType];
  TGraphAsymmErrors* aPt_p[nTrigType][maxNTrig2];
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
	TH1F* hEmpty = new TH1F("hEmpty", "p_{T}^{reco};Efficiency", aPt_p[iter][iter2]->GetN(), aPt_p[iter][iter2]->GetXaxis()->GetXmin(), aPt_p[iter][iter2]->GetXaxis()->GetXmax());
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


  std::cout << error << std::endl;
  error++;

  TFile* outFile_p = new TFile(outFile.c_str(), "UPDATE");

  for(Int_t iter = 0; iter < nTrigType; iter++){
    trigCanvPt_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(trigCanvPt_p[iter], Form("pdfDir/%s", canvPtName[iter].c_str()), "pdf");
  }

  std::cout << error << std::endl;
  error++;

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    delete trigCanvPt_p[iter];
  }

  std::cout << error << std::endl;
  error++;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 4){
    std::cout << "Usage: plotTrigTurnOn_HI <inHistFile> <inTrigFileName> <outFile>" << std::endl;
    std::cout << "argNum: " << argc << std::endl;
    for(Int_t iter = 0; iter < argc; iter++){
      std::cout << "arg " << iter << ": " << argv[iter] << std::endl;
    }

    return -1;
  }

  int rStatus = -1;

  rStatus = plotTrigTurnOn_HI(argv[1], argv[2], argv[3]);

  return rStatus;
}
