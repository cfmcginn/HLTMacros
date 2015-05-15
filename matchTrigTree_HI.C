#include "EventMatchingCMS.h"
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <iostream>
#include <TStyle.h>
#include "TH1F.h"
#include "TMath.h"

#include <string>
#include <vector>
#include <fstream>

const TString Ana3CaloTreename = "akPu3CaloJetAnalyzer/t";
const TString Ana4CaloTreename = "akPu4CaloJetAnalyzer/t";
const TString AnaPhotonTreename = "multiPhotonAnalyzer/photon";
const TString AnaHITreename = "hiEvtAnalyzer/HiTree";
const TString AnaSkimTreename = "skimanalysis/HltTree";
const TString AnaTrkTreename = "anaTrack/trackTree";

const TString HLTFilename = "openHLT_20150508_HIMinBias502_740F.root";

const int nBins = 200;
const double maxpt = 200;

int matchTrigTree_HI(const std::string inHLTFile, const std::string inForestFile, const std::string inTrigFileName, const std::string outFile)
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


  TFile *HLTFile_p =  new TFile(inHLTFile.c_str(), "READ");
  TTree *HLTTree = (TTree*)HLTFile_p->Get("hltbitanalysis/HltTree");

  ULong64_t hlt_event;
  Int_t hlt_run, hlt_lumi;

  Int_t maxNTrig = -1;
  for(Int_t iter = 0; iter < nTrigType; iter++){
    if(trigTypeCount[iter] > maxNTrig) maxNTrig = trigTypeCount[iter];
  }

  const Int_t maxNTrig2 = maxNTrig;

  std::string trigName[nTrigType][maxNTrig2];
  Int_t trigVal[nTrigType][maxNTrig2];
  Int_t nTrigFire[nTrigType][maxNTrig2];

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
	trigVal[typeIter][tempPosIter] = 0;
	nTrigFire[typeIter][tempPosIter] = 0;
	tempPosIter++;
        break;
      }
    }
  }

  HLTTree->SetBranchStatus("*", 0);

  HLTTree->SetBranchStatus("Event", 1);
  HLTTree->SetBranchStatus("Run", 1);
  HLTTree->SetBranchStatus("LumiBlock", 1);

  HLTTree->SetBranchAddress("Event", &hlt_event);
  HLTTree->SetBranchAddress("Run", &hlt_run);
  HLTTree->SetBranchAddress("LumiBlock", &hlt_lumi);


  for(Int_t iter = 0; iter < nTrigType; iter++){
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      HLTTree->SetBranchStatus(trigName[iter][iter2].c_str(), 1);
      HLTTree->SetBranchAddress(trigName[iter][iter2].c_str(), &(trigVal[iter][iter2]));
    }
  }

  TFile *AnaFile_p = new TFile(inForestFile.c_str(), "READ");
  TTree *Ana3CaloTree = (TTree*)AnaFile_p->Get(Ana3CaloTreename); 
  TTree *Ana4CaloTree = (TTree*)AnaFile_p->Get(Ana4CaloTreename); 
  TTree *AnaPhotonTree = (TTree*)AnaFile_p->Get(AnaPhotonTreename); 
  TTree *AnaHITree = (TTree*)AnaFile_p->Get(AnaHITreename); 
  TTree *AnaSkimTree = (TTree*)AnaFile_p->Get(AnaSkimTreename);
  TTree *AnaTrkTree = (TTree*)AnaFile_p->Get(AnaTrkTreename);

  Int_t ana_event, ana_lumi;//, ana_run, ana_lumi;
  Int_t hiBin;

  Int_t pcollisionEventSelection;

  const Int_t maxTrk = 50000;
  Int_t nTrk;
  Bool_t trkFake[maxTrk];
  Float_t trkPt[maxTrk], trkPhi[maxTrk], trkEta[maxTrk];

  const Int_t maxJt = 500;
  Int_t n3Caloref;
  Float_t jt3Calopt[maxJt], jt3Caloeta[maxJt], jt3Calophi[maxJt];
  Float_t ref3Calopt[maxJt], ref3Caloeta[maxJt], ref3Calophi[maxJt];
  Int_t n3Calogen;
  Float_t gen3Calopt[maxJt], gen3Caloeta[maxJt], gen3Calophi[maxJt];

  Int_t n4Caloref;
  Float_t jt4Calopt[maxJt], jt4Caloeta[maxJt], jt4Calophi[maxJt];
  Float_t ref4Calopt[maxJt], ref4Caloeta[maxJt], ref4Calophi[maxJt];
  Int_t n4Calogen;
  Float_t gen4Calopt[maxJt], gen4Caloeta[maxJt], gen4Calophi[maxJt];

  const Int_t maxGamma = 50;
  Int_t nPhotons;
  Float_t photonPt[maxGamma], photonEta[maxGamma], photonPhi[maxGamma];   //[nPhotons]

  AnaHITree->SetBranchStatus("*", 0);
  AnaHITree->SetBranchStatus("evt", 1);
  AnaHITree->SetBranchStatus("lumi", 1);
  AnaHITree->SetBranchStatus("hiBin", 1);
  AnaHITree->SetBranchAddress("evt", &ana_event);
  AnaHITree->SetBranchAddress("lumi", &ana_lumi);
  AnaHITree->SetBranchAddress("hiBin", &hiBin);

  AnaSkimTree->SetBranchStatus("*", 0);
  AnaSkimTree->SetBranchStatus("pcollisionEventSelection", 1);
  AnaSkimTree->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection);


  AnaTrkTree->SetBranchStatus("*", 0);
  AnaTrkTree->SetBranchStatus("nTrk", 1);
  AnaTrkTree->SetBranchStatus("trkFake", 1);
  AnaTrkTree->SetBranchStatus("trkPt", 1);
  AnaTrkTree->SetBranchStatus("trkPhi", 1);
  AnaTrkTree->SetBranchStatus("trkEta", 1);
  AnaTrkTree->SetBranchAddress("nTrk", &nTrk);
  AnaTrkTree->SetBranchAddress("trkFake", trkFake);
  AnaTrkTree->SetBranchAddress("trkPt", trkPt);
  AnaTrkTree->SetBranchAddress("trkPhi", trkPhi);
  AnaTrkTree->SetBranchAddress("trkEta", trkEta);


  Ana3CaloTree->SetBranchStatus("*", 0);
  Ana3CaloTree->SetBranchStatus("nref", 1);
  Ana3CaloTree->SetBranchStatus("jtpt", 1);
  Ana3CaloTree->SetBranchStatus("jteta", 1);
  Ana3CaloTree->SetBranchStatus("jtphi", 1);
  Ana3CaloTree->SetBranchStatus("refpt", 1);
  Ana3CaloTree->SetBranchStatus("refeta", 1);
  Ana3CaloTree->SetBranchStatus("refphi", 1);
  Ana3CaloTree->SetBranchStatus("ngen", 1);
  Ana3CaloTree->SetBranchStatus("genpt", 1);
  Ana3CaloTree->SetBranchStatus("genphi", 1);
  Ana3CaloTree->SetBranchStatus("geneta", 1);
  Ana3CaloTree->SetBranchAddress("nref", &n3Caloref);
  Ana3CaloTree->SetBranchAddress("jtpt", jt3Calopt);
  Ana3CaloTree->SetBranchAddress("jteta", jt3Caloeta);
  Ana3CaloTree->SetBranchAddress("jtphi", jt3Calophi);
  Ana3CaloTree->SetBranchAddress("refpt", ref3Calopt);
  Ana3CaloTree->SetBranchAddress("refeta", ref3Caloeta);
  Ana3CaloTree->SetBranchAddress("refphi", ref3Calophi);
  Ana3CaloTree->SetBranchAddress("ngen", &n3Calogen);
  Ana3CaloTree->SetBranchAddress("genpt", gen3Calopt);
  Ana3CaloTree->SetBranchAddress("genphi", gen3Calophi);
  Ana3CaloTree->SetBranchAddress("geneta", gen3Caloeta);

  Ana4CaloTree->SetBranchStatus("*", 0);
  Ana4CaloTree->SetBranchStatus("nref", 1);
  Ana4CaloTree->SetBranchStatus("jtpt", 1);
  Ana4CaloTree->SetBranchStatus("jteta", 1);
  Ana4CaloTree->SetBranchStatus("jtphi", 1);
  Ana4CaloTree->SetBranchStatus("refpt", 1);
  Ana4CaloTree->SetBranchStatus("refeta", 1);
  Ana4CaloTree->SetBranchStatus("refphi", 1);
  Ana4CaloTree->SetBranchStatus("ngen", 1);
  Ana4CaloTree->SetBranchStatus("genpt", 1);
  Ana4CaloTree->SetBranchStatus("genphi", 1);
  Ana4CaloTree->SetBranchStatus("geneta", 1);
  Ana4CaloTree->SetBranchAddress("nref", &n4Caloref);
  Ana4CaloTree->SetBranchAddress("jtpt", jt4Calopt);
  Ana4CaloTree->SetBranchAddress("jteta", jt4Caloeta);
  Ana4CaloTree->SetBranchAddress("jtphi", jt4Calophi);
  Ana4CaloTree->SetBranchAddress("refpt", ref4Calopt);
  Ana4CaloTree->SetBranchAddress("refeta", ref4Caloeta);
  Ana4CaloTree->SetBranchAddress("refphi", ref4Calophi);
  Ana4CaloTree->SetBranchAddress("ngen", &n4Calogen);
  Ana4CaloTree->SetBranchAddress("genpt", gen4Calopt);
  Ana4CaloTree->SetBranchAddress("genphi", gen4Calophi);
  Ana4CaloTree->SetBranchAddress("geneta", gen4Caloeta);

  AnaPhotonTree->SetBranchStatus("*", 0);
  AnaPhotonTree->SetBranchStatus("nPhotons", 1);
  AnaPhotonTree->SetBranchStatus("pt", 1);
  AnaPhotonTree->SetBranchStatus("eta", 1);
  AnaPhotonTree->SetBranchStatus("phi", 1);
  AnaPhotonTree->SetBranchAddress("nPhotons", &nPhotons);
  AnaPhotonTree->SetBranchAddress("pt", photonPt);
  AnaPhotonTree->SetBranchAddress("eta", photonEta);
  AnaPhotonTree->SetBranchAddress("phi", photonPhi);

  const Int_t nPtBins[nTrigType] = {100, 100, 100, 100};
  const Int_t maxPt[nTrigType] = {100, 150, 150, 100};
  const Int_t nEtaBins = 50;
  TH1F *histsPt_p[nTrigType][maxNTrig2+1], *histsEta_p[nTrigType][maxNTrig2+1];

  //edit here
  for(Int_t iter = 0; iter < nTrigType; iter++){
    histsPt_p[iter][0] = new TH1F(Form("leading%s_pt", trigType[iter].c_str()), Form("leading%s_pt", trigType[iter].c_str()), nPtBins[iter], 0.0, maxPt[iter]);

    histsEta_p[iter][0] = new TH1F(Form("leading%s_eta", trigType[iter].c_str()), Form("leading%s_eta", trigType[iter].c_str()), nEtaBins, -5.0, 5.0);

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      histsPt_p[iter][iter2+1] = (TH1F*)histsPt_p[iter][0]->Clone(Form("%s_%s_pt", trigName[iter][iter2].c_str(), trigType[iter].c_str()));

      histsEta_p[iter][iter2+1] = (TH1F*)histsEta_p[iter][0]->Clone(Form("%s_%s_eta", trigName[iter][iter2].c_str(), trigType[iter].c_str()));
    }
  }

  //book histos

  std::cout << "Events in HLT file: " << HLTTree->GetEntries() << std::endl;
  std::cout << "Events in Ana file: " << Ana3CaloTree->GetEntries() << std::endl;

  //make map
  EventMatchingCMS *matcher = new EventMatchingCMS();

  for(Long64_t entry = 0; entry < HLTTree->GetEntries(); ++entry){
    HLTTree->GetEntry(entry);
    
    matcher->addEvent(hlt_event, hlt_lumi, 0, entry);
  }

  // analysis loop
  std::cout << "COMMENCE LOOP" << std::endl;
  int matched = 0;
  for(Long64_t entry = 0; entry < Ana3CaloTree->GetEntries(); ++entry){
    if(entry%10000 == 0) std::cout << entry << std::endl;

    Ana3CaloTree->GetEntry(entry);
    Ana4CaloTree->GetEntry(entry);
    AnaPhotonTree->GetEntry(entry);
    AnaHITree->GetEntry(entry);
    AnaSkimTree->GetEntry(entry);
    AnaTrkTree->GetEntry(entry);

    long long hlt_entry = matcher->retrieveEvent(ana_event, ana_lumi, 0);
    if(hlt_entry == -1) continue;

    HLTTree->GetEntry(hlt_entry);
    matched++;

    for(Int_t iter = 0; iter < nTrigType; iter++){
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	if(trigVal[iter][iter2]) nTrigFire[iter][iter2]++;
      }
    }

    Double_t maxTrkPt = -1;
    Double_t maxTrkEta = -100;
  
    for(int i = 0; i < nTrk; i++){
      if(fabs(trkEta[i]) > 2.4) continue;
      if(trkFake[i]) continue;
      if(trkPt[i] > maxTrkPt){
	maxTrkPt = trkPt[i];
	maxTrkEta = trkEta[i];
      }
    }
  
    Double_t max3CaloAnaPt = -1;
    Double_t max3CaloAnaEta = -100;
    
    for(int i = 0; i < n3Caloref; ++i){
      if(fabs(jt3Caloeta[i]) > 2.0) continue;
      if(jt3Calopt[i] > max3CaloAnaPt){
	max3CaloAnaPt = jt3Calopt[i];
	max3CaloAnaEta = jt3Caloeta[i];
      }
    }

    Double_t max4CaloAnaPt = -1;
    Double_t max4CaloAnaEta = -100;

    for(int i = 0; i < n4Caloref; ++i){
      if(fabs(jt4Caloeta[i]) > 2.0) continue;
      if(jt4Calopt[i] > max4CaloAnaPt){
	max4CaloAnaPt = jt4Calopt[i];
	max4CaloAnaEta = jt4Caloeta[i];
      }
    }

    Double_t maxPhotonAnaPt = -1;
    Double_t maxPhotonAnaEta = -100;

    for(int i = 0; i < nPhotons; ++i){
      if(fabs(photonEta[i]) > 2.0) continue;
      if(photonPt[i] > maxPhotonAnaPt){
	maxPhotonAnaPt = photonPt[i];
	maxPhotonAnaEta = photonEta[i];
      }
    }

    Double_t trigOfflinePt[nTrigType] = {maxTrkPt, max3CaloAnaPt, max4CaloAnaPt, maxPhotonAnaPt};
    Double_t trigOfflineEta[nTrigType] = {maxTrkEta, max3CaloAnaEta, max4CaloAnaEta, maxPhotonAnaEta};

    for(Int_t iter = 0; iter < nTrigType; iter++){
      if(trigOfflinePt[iter] > 0){
	histsPt_p[iter][0]->Fill(trigOfflinePt[iter]);
	histsEta_p[iter][0]->Fill(trigOfflineEta[iter]);

	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(trigVal[iter][iter2]){
	    histsPt_p[iter][iter2+1]->Fill(trigOfflinePt[iter]);
	    histsEta_p[iter][iter2+1]->Fill(trigOfflineEta[iter]);
	  }
	}
      }
    }
  }

  std::cout << std::endl;
  std::cout << "Matched events: " << matched << std::endl;
  std::cout << "Trigger fires: " << std::endl;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    std::cout << "  Trigger Type: " << trigType[iter] << std::endl;
    
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      std::cout << "    " << trigName[iter][iter2] << ": " << nTrigFire[iter][iter2] << std::endl;
    }
  }

  TGraphAsymmErrors *aPt_p[nTrigType][maxNTrig2], *aEta_p[nTrigType][maxNTrig2];

  for(Int_t iter = 0; iter < nTrigType; iter++){
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      aPt_p[iter][iter2] = new TGraphAsymmErrors();
      aPt_p[iter][iter2]->BayesDivide(histsPt_p[iter][iter2+1],histsPt_p[iter][0]);
      aPt_p[iter][iter2]->SetName(Form("%s_%s_pt_asymm", trigName[iter][iter2].c_str(), trigType[iter].c_str()));
      aEta_p[iter][iter2] = new TGraphAsymmErrors();
      aEta_p[iter][iter2]->BayesDivide(histsEta_p[iter][iter2+1],histsEta_p[iter][0]);
      aEta_p[iter][iter2]->SetName(Form("%s_%s_eta_asymm", trigName[iter][iter2].c_str(), trigType[iter].c_str()));
    }
  }

  TFile* outFile_p = new TFile(outFile.c_str(), "UPDATE");

  for(Int_t iter = 0; iter < nTrigType; iter++){
    histsPt_p[iter][0]->Write("", TObject::kOverwrite);
    histsEta_p[iter][0]->Write("", TObject::kOverwrite);

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      histsPt_p[iter][iter2+1]->Write("", TObject::kOverwrite);
      histsEta_p[iter][iter2+1]->Write("", TObject::kOverwrite);

      aPt_p[iter][iter2]->Write("", TObject::kOverwrite);
      aEta_p[iter][iter2]->Write("", TObject::kOverwrite);
    }
  }

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    delete histsPt_p[iter][0];
    delete histsEta_p[iter][0];

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      delete aPt_p[iter][iter2];
      delete aEta_p[iter][iter2];
    }
  }

  AnaFile_p->Close();
  delete AnaFile_p;

  HLTFile_p->Close();
  delete HLTFile_p;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 5){
    std::cout << "Usage: matchTrigTree_HI <inHLTFile> <inForestFile> <inTrigFileName> <outFile>" << std::endl;
    std::cout << "argNum: " << argc << std::endl;
    for(Int_t iter = 0; iter < argc; iter++){
      std::cout << "arg " << iter << ": " << argv[iter] << std::endl;
    }

    return -1;
  }

  int rStatus = -1;

  rStatus = matchTrigTree_HI(argv[1], argv[2], argv[3], argv[4]);

  return rStatus;
}
