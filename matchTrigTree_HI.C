#include "EventMatchingCMS.h"
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <iostream>
#include <iomanip>
#include <TStyle.h>
#include "TH1F.h"
#include "TMath.h"

#include <string>
#include <vector>
#include <fstream>

const TString AnaPu3CaloTreename = "akPu3CaloJetAnalyzer/t";
const TString AnaPu4CaloTreename = "akPu4CaloJetAnalyzer/t";
const TString AnaPu3PFTreename = "akPu3PFJetAnalyzer/t";
const TString AnaPu4PFTreename = "akPu4PFJetAnalyzer/t";
const TString AnaVs3CaloTreename = "akVs3CaloJetAnalyzer/t";
const TString AnaVs4CaloTreename = "akVs4CaloJetAnalyzer/t";
const TString AnaPhotonTreename = "multiPhotonAnalyzer/photon";
const TString AnaHITreename = "hiEvtAnalyzer/HiTree";
const TString AnaSkimTreename = "skimanalysis/HltTree";
const TString AnaTrkTreename = "anaTrack/trackTree";
const TString AnaGenTreename = "HiGenParticleAna/hi";

const TString HLTFilename = "openHLT_20150508_HIMinBias502_740F.root";

int matchTrigTree_HI(const std::string inHLTFile, const std::string inForestFile, const std::string inTrigFileName)
{
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
    std::cout << "Gettng Trig List" << std::endl;
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
  Int_t trigTypeCount[nTrigType];

  nLines = 0;

  for(Int_t iter = 0; iter < (Int_t)listOfTrig.size(); iter++){
    if(std::string::npos == listOfTrig[iter].find("HLT") && std::string::npos == listOfTrig[iter].find("L1")){

      std::size_t strIndex = 0;
      while(true){
        strIndex = listOfTrig[iter].find(",");
        if(strIndex == std::string::npos) break;
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
  Int_t trigThresh[nTrigType][maxNTrig2];
  Int_t trigVal[nTrigType][maxNTrig2];
  Int_t nTrigFire[nTrigType][maxNTrig2];
  Float_t highPtMiss[nTrigType][maxNTrig2];

  Float_t nZTrigFireNum[nTrigType][maxNTrig2];
  Float_t nZTrigFireDenom[nTrigType][maxNTrig2];

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
      trigVal[nLines][tempPosIter] = 0;
      nTrigFire[nLines][tempPosIter] = 0;
      highPtMiss[nLines][tempPosIter] = 0;
      nZTrigFireNum[nLines][tempPosIter] = 0;
      nZTrigFireDenom[nLines][tempPosIter] = 0;
      tempPosIter++;
    }
  }

  HLTTree->SetBranchStatus("*", 0);

  HLTTree->SetBranchStatus("Event", 1);
  HLTTree->SetBranchStatus("Run", 1);
  HLTTree->SetBranchStatus("LumiBlock", 1);

  HLTTree->SetBranchAddress("Event", &hlt_event);
  HLTTree->SetBranchAddress("Run", &hlt_run);
  HLTTree->SetBranchAddress("LumiBlock", &hlt_lumi);

  //EDIT HERE FOR TRIGGER ANDING
  //prev -8
  for(Int_t iter = 0; iter < nTrigType; iter++){
    //    if(iter == 12 || iter == 13 || iter == 14) continue;
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      HLTTree->SetBranchStatus(trigName[iter][iter2].c_str(), 1);
      HLTTree->SetBranchAddress(trigName[iter][iter2].c_str(), &(trigVal[iter][iter2]));
      std::cout << trigName[iter][iter2].c_str() << std::endl;
    }
  }

  TFile *AnaFile_p = new TFile(inForestFile.c_str(), "READ");
  TTree *AnaPu3CaloTree = (TTree*)AnaFile_p->Get(AnaPu3CaloTreename); 
  TTree *AnaPu4CaloTree = (TTree*)AnaFile_p->Get(AnaPu4CaloTreename); 
  TTree *AnaPu3PFTree = (TTree*)AnaFile_p->Get(AnaPu3PFTreename); 
  TTree *AnaPu4PFTree = (TTree*)AnaFile_p->Get(AnaPu4PFTreename); 
  TTree *AnaVs3CaloTree = (TTree*)AnaFile_p->Get(AnaVs3CaloTreename); 
  TTree *AnaVs4CaloTree = (TTree*)AnaFile_p->Get(AnaVs4CaloTreename); 
  TTree *AnaPhotonTree = (TTree*)AnaFile_p->Get(AnaPhotonTreename); 
  TTree *AnaHITree = (TTree*)AnaFile_p->Get(AnaHITreename); 
  TTree *AnaSkimTree = (TTree*)AnaFile_p->Get(AnaSkimTreename);
  TTree *AnaTrkTree = (TTree*)AnaFile_p->Get(AnaTrkTreename);
  TTree *AnaGenTree = (TTree*)AnaFile_p->Get(AnaGenTreename);
  

  Int_t ana_event, ana_lumi;//, ana_run, ana_lumi;
  Int_t hiBin;

  Int_t pcollisionEventSelection;

  const Int_t maxTrk = 50000;
  Int_t nTrk;
  Bool_t trkFake[maxTrk];
  Float_t trkPt[maxTrk], trkPhi[maxTrk], trkEta[maxTrk];
  
  const Int_t maxGen = 100000;
  Int_t nGen;
  Float_t genPt[maxGen], genPhi[maxGen], genEta[maxGen];
  Int_t genPDG[maxGen], genStatus[maxGen];
  
  const Int_t maxJt = 500;
  Int_t nPu3Caloref;
  Float_t jtPu3Calopt[maxJt], jtPu3Caloeta[maxJt], jtPu3Calophi[maxJt];
  Float_t refPu3Calopt[maxJt], refPu3Caloeta[maxJt], refPu3Calophi[maxJt];
  Int_t nPu3Calogen;
  Float_t genPu3Calopt[maxJt], genPu3Caloeta[maxJt], genPu3Calophi[maxJt];

  Int_t nPu4Caloref;
  Float_t jtPu4Calopt[maxJt], jtPu4Caloeta[maxJt], jtPu4Calophi[maxJt];
  Float_t refPu4Calopt[maxJt], refPu4Caloeta[maxJt], refPu4Calophi[maxJt];
  Int_t nPu4Calogen;
  Float_t genPu4Calopt[maxJt], genPu4Caloeta[maxJt], genPu4Calophi[maxJt];

  Int_t nPu3PFref;
  Float_t jtPu3PFpt[maxJt], jtPu3PFeta[maxJt], jtPu3PFphi[maxJt];
  Float_t refPu3PFpt[maxJt], refPu3PFeta[maxJt], refPu3PFphi[maxJt];
  Int_t nPu3PFgen;
  Float_t genPu3PFpt[maxJt], genPu3PFeta[maxJt], genPu3PFphi[maxJt];

  Int_t nPu4PFref;
  Float_t jtPu4PFpt[maxJt], jtPu4PFeta[maxJt], jtPu4PFphi[maxJt];
  Float_t refPu4PFpt[maxJt], refPu4PFeta[maxJt], refPu4PFphi[maxJt];
  Int_t nPu4PFgen;
  Float_t genPu4PFpt[maxJt], genPu4PFeta[maxJt], genPu4PFphi[maxJt];

  Int_t nVs3Caloref;
  Float_t jtVs3Calopt[maxJt], jtVs3Caloeta[maxJt], jtVs3Calophi[maxJt];
  Float_t refVs3Calopt[maxJt], refVs3Caloeta[maxJt], refVs3Calophi[maxJt];
  Int_t nVs3Calogen;
  Float_t genVs3Calopt[maxJt], genVs3Caloeta[maxJt], genVs3Calophi[maxJt];

  Int_t nVs4Caloref;
  Float_t jtVs4Calopt[maxJt], jtVs4Caloeta[maxJt], jtVs4Calophi[maxJt];
  Float_t refVs4Calopt[maxJt], refVs4Caloeta[maxJt], refVs4Calophi[maxJt];
  Int_t nVs4Calogen;
  Float_t genVs4Calopt[maxJt], genVs4Caloeta[maxJt], genVs4Calophi[maxJt];

  const Int_t maxGamma = 50;
  Int_t nPhotons;
  Float_t photonPt[maxGamma], photonEta[maxGamma], photonPhi[maxGamma];   //[nPhotons]
  Float_t seedTime[maxGamma], swissCrx[maxGamma], sigmaIphiIphi[maxGamma], sigmaIetaIeta[maxGamma];
  Int_t isEle[maxGamma];

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
 
  AnaGenTree->SetBranchStatus("*", 0);
  
  AnaGenTree->SetBranchStatus("mult", 1);
  AnaGenTree->SetBranchStatus("pt", 1);
  AnaGenTree->SetBranchStatus("phi", 1);
  AnaGenTree->SetBranchStatus("eta", 1);
  AnaGenTree->SetBranchStatus("pdg", 1);
  AnaGenTree->SetBranchStatus("sta", 1);

  AnaGenTree->SetBranchAddress("mult", &nGen);
  AnaGenTree->SetBranchAddress("pt", genPt);
  AnaGenTree->SetBranchAddress("phi", genPhi);
  AnaGenTree->SetBranchAddress("eta", genEta);
  AnaGenTree->SetBranchAddress("pdg", genPDG);
  AnaGenTree->SetBranchAddress("sta", genStatus);
  

  AnaPu3CaloTree->SetBranchStatus("*", 0);
  AnaPu3CaloTree->SetBranchStatus("nref", 1);
  AnaPu3CaloTree->SetBranchStatus("jtpt", 1);
  AnaPu3CaloTree->SetBranchStatus("jteta", 1);
  AnaPu3CaloTree->SetBranchStatus("jtphi", 1);
  AnaPu3CaloTree->SetBranchStatus("refpt", 1);
  AnaPu3CaloTree->SetBranchStatus("refeta", 1);
  AnaPu3CaloTree->SetBranchStatus("refphi", 1);
  AnaPu3CaloTree->SetBranchStatus("ngen", 1);
  AnaPu3CaloTree->SetBranchStatus("genpt", 1);
  AnaPu3CaloTree->SetBranchStatus("genphi", 1);
  AnaPu3CaloTree->SetBranchStatus("geneta", 1);
  AnaPu3CaloTree->SetBranchAddress("nref", &nPu3Caloref);
  AnaPu3CaloTree->SetBranchAddress("jtpt", jtPu3Calopt);
  AnaPu3CaloTree->SetBranchAddress("jteta", jtPu3Caloeta);
  AnaPu3CaloTree->SetBranchAddress("jtphi", jtPu3Calophi);
  AnaPu3CaloTree->SetBranchAddress("refpt", refPu3Calopt);
  AnaPu3CaloTree->SetBranchAddress("refeta", refPu3Caloeta);
  AnaPu3CaloTree->SetBranchAddress("refphi", refPu3Calophi);
  AnaPu3CaloTree->SetBranchAddress("ngen", &nPu3Calogen);
  AnaPu3CaloTree->SetBranchAddress("genpt", genPu3Calopt);
  AnaPu3CaloTree->SetBranchAddress("genphi", genPu3Calophi);
  AnaPu3CaloTree->SetBranchAddress("geneta", genPu3Caloeta);

  AnaPu4CaloTree->SetBranchStatus("*", 0);
  AnaPu4CaloTree->SetBranchStatus("nref", 1);
  AnaPu4CaloTree->SetBranchStatus("jtpt", 1);
  AnaPu4CaloTree->SetBranchStatus("jteta", 1);
  AnaPu4CaloTree->SetBranchStatus("jtphi", 1);
  AnaPu4CaloTree->SetBranchStatus("refpt", 1);
  AnaPu4CaloTree->SetBranchStatus("refeta", 1);
  AnaPu4CaloTree->SetBranchStatus("refphi", 1);
  AnaPu4CaloTree->SetBranchStatus("ngen", 1);
  AnaPu4CaloTree->SetBranchStatus("genpt", 1);
  AnaPu4CaloTree->SetBranchStatus("genphi", 1);
  AnaPu4CaloTree->SetBranchStatus("geneta", 1);
  AnaPu4CaloTree->SetBranchAddress("nref", &nPu4Caloref);
  AnaPu4CaloTree->SetBranchAddress("jtpt", jtPu4Calopt);
  AnaPu4CaloTree->SetBranchAddress("jteta", jtPu4Caloeta);
  AnaPu4CaloTree->SetBranchAddress("jtphi", jtPu4Calophi);
  AnaPu4CaloTree->SetBranchAddress("refpt", refPu4Calopt);
  AnaPu4CaloTree->SetBranchAddress("refeta", refPu4Caloeta);
  AnaPu4CaloTree->SetBranchAddress("refphi", refPu4Calophi);
  AnaPu4CaloTree->SetBranchAddress("ngen", &nPu4Calogen);
  AnaPu4CaloTree->SetBranchAddress("genpt", genPu4Calopt);
  AnaPu4CaloTree->SetBranchAddress("genphi", genPu4Calophi);
  AnaPu4CaloTree->SetBranchAddress("geneta", genPu4Caloeta);

  AnaPu3PFTree->SetBranchStatus("*", 0);
  AnaPu3PFTree->SetBranchStatus("nref", 1);
  AnaPu3PFTree->SetBranchStatus("jtpt", 1);
  AnaPu3PFTree->SetBranchStatus("jteta", 1);
  AnaPu3PFTree->SetBranchStatus("jtphi", 1);
  AnaPu3PFTree->SetBranchStatus("refpt", 1);
  AnaPu3PFTree->SetBranchStatus("refeta", 1);
  AnaPu3PFTree->SetBranchStatus("refphi", 1);
  AnaPu3PFTree->SetBranchStatus("ngen", 1);
  AnaPu3PFTree->SetBranchStatus("genpt", 1);
  AnaPu3PFTree->SetBranchStatus("genphi", 1);
  AnaPu3PFTree->SetBranchStatus("geneta", 1);
  AnaPu3PFTree->SetBranchAddress("nref", &nPu3PFref);
  AnaPu3PFTree->SetBranchAddress("jtpt", jtPu3PFpt);
  AnaPu3PFTree->SetBranchAddress("jteta", jtPu3PFeta);
  AnaPu3PFTree->SetBranchAddress("jtphi", jtPu3PFphi);
  AnaPu3PFTree->SetBranchAddress("refpt", refPu3PFpt);
  AnaPu3PFTree->SetBranchAddress("refeta", refPu3PFeta);
  AnaPu3PFTree->SetBranchAddress("refphi", refPu3PFphi);
  AnaPu3PFTree->SetBranchAddress("ngen", &nPu3PFgen);
  AnaPu3PFTree->SetBranchAddress("genpt", genPu3PFpt);
  AnaPu3PFTree->SetBranchAddress("genphi", genPu3PFphi);
  AnaPu3PFTree->SetBranchAddress("geneta", genPu3PFeta);

  AnaPu4PFTree->SetBranchStatus("*", 0);
  AnaPu4PFTree->SetBranchStatus("nref", 1);
  AnaPu4PFTree->SetBranchStatus("jtpt", 1);
  AnaPu4PFTree->SetBranchStatus("jteta", 1);
  AnaPu4PFTree->SetBranchStatus("jtphi", 1);
  AnaPu4PFTree->SetBranchStatus("refpt", 1);
  AnaPu4PFTree->SetBranchStatus("refeta", 1);
  AnaPu4PFTree->SetBranchStatus("refphi", 1);
  AnaPu4PFTree->SetBranchStatus("ngen", 1);
  AnaPu4PFTree->SetBranchStatus("genpt", 1);
  AnaPu4PFTree->SetBranchStatus("genphi", 1);
  AnaPu4PFTree->SetBranchStatus("geneta", 1);
  AnaPu4PFTree->SetBranchAddress("nref", &nPu4PFref);
  AnaPu4PFTree->SetBranchAddress("jtpt", jtPu4PFpt);
  AnaPu4PFTree->SetBranchAddress("jteta", jtPu4PFeta);
  AnaPu4PFTree->SetBranchAddress("jtphi", jtPu4PFphi);
  AnaPu4PFTree->SetBranchAddress("refpt", refPu4PFpt);
  AnaPu4PFTree->SetBranchAddress("refeta", refPu4PFeta);
  AnaPu4PFTree->SetBranchAddress("refphi", refPu4PFphi);
  AnaPu4PFTree->SetBranchAddress("ngen", &nPu4PFgen);
  AnaPu4PFTree->SetBranchAddress("genpt", genPu4PFpt);
  AnaPu4PFTree->SetBranchAddress("genphi", genPu4PFphi);
  AnaPu4PFTree->SetBranchAddress("geneta", genPu4PFeta);

  AnaVs3CaloTree->SetBranchStatus("*", 0);
  AnaVs3CaloTree->SetBranchStatus("nref", 1);
  AnaVs3CaloTree->SetBranchStatus("jtpt", 1);
  AnaVs3CaloTree->SetBranchStatus("jteta", 1);
  AnaVs3CaloTree->SetBranchStatus("jtphi", 1);
  AnaVs3CaloTree->SetBranchStatus("refpt", 1);
  AnaVs3CaloTree->SetBranchStatus("refeta", 1);
  AnaVs3CaloTree->SetBranchStatus("refphi", 1);
  AnaVs3CaloTree->SetBranchStatus("ngen", 1);
  AnaVs3CaloTree->SetBranchStatus("genpt", 1);
  AnaVs3CaloTree->SetBranchStatus("genphi", 1);
  AnaVs3CaloTree->SetBranchStatus("geneta", 1);
  AnaVs3CaloTree->SetBranchAddress("nref", &nVs3Caloref);
  AnaVs3CaloTree->SetBranchAddress("jtpt", jtVs3Calopt);
  AnaVs3CaloTree->SetBranchAddress("jteta", jtVs3Caloeta);
  AnaVs3CaloTree->SetBranchAddress("jtphi", jtVs3Calophi);
  AnaVs3CaloTree->SetBranchAddress("refpt", refVs3Calopt);
  AnaVs3CaloTree->SetBranchAddress("refeta", refVs3Caloeta);
  AnaVs3CaloTree->SetBranchAddress("refphi", refVs3Calophi);
  AnaVs3CaloTree->SetBranchAddress("ngen", &nVs3Calogen);
  AnaVs3CaloTree->SetBranchAddress("genpt", genVs3Calopt);
  AnaVs3CaloTree->SetBranchAddress("genphi", genVs3Calophi);
  AnaVs3CaloTree->SetBranchAddress("geneta", genVs3Caloeta);

  AnaVs4CaloTree->SetBranchStatus("*", 0);
  AnaVs4CaloTree->SetBranchStatus("nref", 1);
  AnaVs4CaloTree->SetBranchStatus("jtpt", 1);
  AnaVs4CaloTree->SetBranchStatus("jteta", 1);
  AnaVs4CaloTree->SetBranchStatus("jtphi", 1);
  AnaVs4CaloTree->SetBranchStatus("refpt", 1);
  AnaVs4CaloTree->SetBranchStatus("refeta", 1);
  AnaVs4CaloTree->SetBranchStatus("refphi", 1);
  AnaVs4CaloTree->SetBranchStatus("ngen", 1);
  AnaVs4CaloTree->SetBranchStatus("genpt", 1);
  AnaVs4CaloTree->SetBranchStatus("genphi", 1);
  AnaVs4CaloTree->SetBranchStatus("geneta", 1);
  AnaVs4CaloTree->SetBranchAddress("nref", &nVs4Caloref);
  AnaVs4CaloTree->SetBranchAddress("jtpt", jtVs4Calopt);
  AnaVs4CaloTree->SetBranchAddress("jteta", jtVs4Caloeta);
  AnaVs4CaloTree->SetBranchAddress("jtphi", jtVs4Calophi);
  AnaVs4CaloTree->SetBranchAddress("refpt", refVs4Calopt);
  AnaVs4CaloTree->SetBranchAddress("refeta", refVs4Caloeta);
  AnaVs4CaloTree->SetBranchAddress("refphi", refVs4Calophi);
  AnaVs4CaloTree->SetBranchAddress("ngen", &nVs4Calogen);
  AnaVs4CaloTree->SetBranchAddress("genpt", genVs4Calopt);
  AnaVs4CaloTree->SetBranchAddress("genphi", genVs4Calophi);
  AnaVs4CaloTree->SetBranchAddress("geneta", genVs4Caloeta);

  AnaPhotonTree->SetBranchStatus("*", 0);
  AnaPhotonTree->SetBranchStatus("nPhotons", 1);
  AnaPhotonTree->SetBranchStatus("pt", 1);
  AnaPhotonTree->SetBranchStatus("eta", 1);
  AnaPhotonTree->SetBranchStatus("phi", 1);
  AnaPhotonTree->SetBranchStatus("seedTime", 1);
  AnaPhotonTree->SetBranchStatus("swissCrx", 1);
  AnaPhotonTree->SetBranchStatus("sigmaIphiIphi", 1);
  AnaPhotonTree->SetBranchStatus("sigmaIetaIeta", 1);
  AnaPhotonTree->SetBranchStatus("isEle", 1);
  
  AnaPhotonTree->SetBranchAddress("nPhotons", &nPhotons);
  AnaPhotonTree->SetBranchAddress("pt", photonPt);
  AnaPhotonTree->SetBranchAddress("eta", photonEta);
  AnaPhotonTree->SetBranchAddress("phi", photonPhi);
  AnaPhotonTree->SetBranchAddress("seedTime", seedTime);
  AnaPhotonTree->SetBranchAddress("swissCrx", swissCrx);
  AnaPhotonTree->SetBranchAddress("sigmaIphiIphi", sigmaIphiIphi);
  AnaPhotonTree->SetBranchAddress("sigmaIetaIeta", sigmaIetaIeta);
  AnaPhotonTree->SetBranchAddress("isEle", isEle);

  //FOR ADDITIONAL OFFLINE OBJECT MATCHING, EDIT HERE (1 of 2)


  //  const Int_t nPtBins[nTrigType] = {200, 200, 200, 200};
  //  const Int_t maxPt[nTrigType] = {200, 200, 200, 200};

  const Int_t nPtBins[nTrigType] = {20};
  const Int_t maxPt[nTrigType] = {20};

  //  const Int_t nPtBins[nTrigType] = {120, 60, 60, 60, 60, 30, 100, 60, 160, 160, 160, 100, 100, 100, 100, 160, 160, 160, 100, 100, 100, 100, 100, 100, 160, 30, 160, 160, 100};
  //  const Int_t maxPt[nTrigType] = {120, 60, 60, 60, 60, 30, 100, 60, 160, 160, 160, 100, 100, 100, 100, 160, 160, 160, 100, 100, 100, 100, 100, 100, 160, 30, 160, 160, 100};

  //  const Int_t nPtBins[nTrigType] = {100, 200, 200, 200, 200, 200, 200, 100, 100, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};
  //  const Int_t maxPt[nTrigType] = {100, 200, 200, 200, 200, 200, 200, 100, 100, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};
  //const Int_t nPtBins[nTrigType] = {120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120};
  //const Int_t maxPt[nTrigType] = {130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130};

  const Int_t nEtaBins = 50;
  TH1F *histsPt_p[nTrigType][maxNTrig2+1], *histsEta_p[nTrigType][maxNTrig2], *histsEtaTrig_p[nTrigType][maxNTrig2];
  TH1F* histsSpectPt_p[nTrigType];

  TH1F* ratesMatched_p[nTrigType];
  TH1F* ratesUnmatched_p[nTrigType];

  for(Int_t iter = 0; iter < nTrigType; iter++){
    histsPt_p[iter][0] = new TH1F(Form("leading%s_pt", trigType[iter].c_str()), Form("leading%s_pt", trigType[iter].c_str()), nPtBins[iter], 0.0, maxPt[iter]);

    histsSpectPt_p[iter] = new TH1F(Form("spect%s_pt", trigType[iter].c_str()), Form("spect%s_pt", trigType[iter].c_str()), nPtBins[iter], 0.0, maxPt[iter]);

    Int_t max = trigThresh[iter][trigTypeCount[iter]-1];
    Int_t min = trigThresh[iter][0];
    Float_t minDiff = trigThresh[iter][1] - trigThresh[iter][0];

    for(Int_t iter2 = 1; iter2 < trigTypeCount[iter]-1; iter2++){
      if(trigThresh[iter][iter2+1] - trigThresh[iter][iter2] < minDiff) minDiff = trigThresh[iter][iter2+1] - trigThresh[iter][iter2];
    }

    minDiff /= 2;
    if(min < minDiff) minDiff = min;

    Int_t nBins = (max+minDiff - (min-minDiff))/(minDiff*2);
    if(trigTypeCount[iter] == 1){
      nBins = 1;
      minDiff = 20;
      if(min < minDiff) minDiff = min;
      min = trigThresh[iter][0] - minDiff;
      max = trigThresh[iter][0] + minDiff;
    }
    else{
      min -= minDiff;
      max += minDiff;
    }

    ratesMatched_p[iter] = new TH1F(Form("ratesMatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), Form("ratesMatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), nBins, min, max);
    ratesUnmatched_p[iter] = new TH1F(Form("ratesUnmatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), Form("ratesUnmatched_%s_%d_%d", trigType[iter].c_str(), trigThresh[iter][0], trigThresh[iter][trigTypeCount[iter]-1]), nBins, min, max);

    std::cout << "Min, max: " << min << ", " << max << std::endl;

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      histsPt_p[iter][iter2+1] = (TH1F*)histsPt_p[iter][0]->Clone(Form("%s_%s_%d_pt", trigName[iter][iter2].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]));

      histsEta_p[iter][iter2] = new TH1F(Form("leading_%s_%d_eta", trigType[iter].c_str(), trigThresh[iter][iter2]), Form("leading_%s_%d_eta", trigType[iter].c_str(), trigThresh[iter][iter2]), nEtaBins, -5.0, 5.0);

      histsEtaTrig_p[iter][iter2] = new TH1F(Form("%s_%s_%d_eta", trigName[iter][iter2].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]), Form("%s_%s_%d_eta", trigName[iter][iter2].c_str(), trigType[iter].c_str(), trigThresh[iter][iter2]), nEtaBins, -5.0, 5.0);
    }
  }

  //book histos

  std::cout << "Events in HLT file: " << HLTTree->GetEntries() << std::endl;
  std::cout << "Events in Ana file: " << AnaPu3CaloTree->GetEntries() << std::endl;

  //make map
  EventMatchingCMS *matcher = new EventMatchingCMS();

  for(Long64_t entry = 0; entry < HLTTree->GetEntries(); ++entry){
    HLTTree->GetEntry(entry);

    for(Int_t iter = 0; iter < nTrigType; iter++){
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	//EDIT HERE FOR ANDING
	/*	
	if(iter >= nTrigType-8 && iter < nTrigType-6) trigVal[iter][iter2] = (trigVal[iter-4][0] && trigVal[3][iter2]);
	else if(iter >= nTrigType-6 && iter < nTrigType-4) trigVal[iter][iter2] = (trigVal[iter-2][0] && trigVal[3][iter2]);
	else if(iter >= nTrigType-4 && iter < nTrigType-2) trigVal[iter][iter2] = (trigVal[iter-6][0] && trigVal[3][iter2]);
	else if(iter >= nTrigType-2) trigVal[iter][iter2] = (trigVal[iter-2][0] && trigVal[3][iter2]);

	if(iter >= nTrigType-6 && iter < nTrigType-4 && trigVal[3][trigTypeCount[3] - 1]) trigVal[iter][iter2] = 0;
	else if(iter >= nTrigType-2 && trigVal[3][trigTypeCount[3] - 1]) trigVal[iter][iter2] = 0;

	if(iter == 12 || iter == 13 || iter == 14) trigVal[iter][iter2] =  trigVal[iter-3][iter2];
	if((iter == 12 || iter == 13 || iter == 14) && trigVal[3][trigTypeCount[3] - 1])  trigVal[iter][iter2] = 0;
	*/

	if(trigVal[iter][iter2]) nTrigFire[iter][iter2]++;
      }
    }
    
    matcher->addEvent(hlt_event, hlt_lumi, 0, entry);
  }

  // analysis loop
  std::cout << "COMMENCE LOOP" << std::endl;
  int matched = 0;
  for(Long64_t entry = 0; entry < AnaPu3CaloTree->GetEntries(); ++entry){
    if(entry%10000 == 0) std::cout << entry << std::endl;

    AnaPu3CaloTree->GetEntry(entry);
    AnaPu4CaloTree->GetEntry(entry);
    AnaPu3PFTree->GetEntry(entry);
    AnaPu4PFTree->GetEntry(entry);
    AnaVs3CaloTree->GetEntry(entry);
    AnaVs4CaloTree->GetEntry(entry);
    AnaPhotonTree->GetEntry(entry);
    AnaHITree->GetEntry(entry);
    AnaSkimTree->GetEntry(entry);
    AnaTrkTree->GetEntry(entry);
    AnaGenTree->GetEntry(entry);

    long long hlt_entry = matcher->retrieveEvent(ana_event, ana_lumi, 0);
    if(hlt_entry == -1) continue;

    HLTTree->GetEntry(hlt_entry);
    matched++;
    
    /*
    for(Int_t iter = 0; iter < nTrigType; iter++){
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	
        //EDIT HERE FOR ANDING                                                                       
        if(iter >= nTrigType-8 && iter < nTrigType-6) trigVal[iter][iter2] = (trigVal[iter-4][0] && trigVal[3][iter2]);
        else if(iter >= nTrigType-6 && iter < nTrigType-4) trigVal[iter][iter2] = (trigVal[iter-2][0] && trigVal[3][iter2]);
        else if(iter >= nTrigType-4 && iter < nTrigType-2) trigVal[iter][iter2] = (trigVal[iter-6][0]&& trigVal[3][iter2]);
        else if(iter >= nTrigType-2) trigVal[iter][iter2] = (trigVal[iter-2][0] && trigVal[3][iter2]);

        if(iter >= nTrigType-6 && iter < nTrigType-4 && trigVal[3][trigTypeCount[3] - 1]) trigVal[iter][iter2] = 0;
        else if(iter >= nTrigType-2 && trigVal[3][trigTypeCount[3] - 1]) trigVal[iter][iter2] = 0;

        if(iter == 12 || iter == 13 || iter == 14) trigVal[iter][iter2] =  trigVal[iter-3][iter2];
        if((iter == 12 || iter == 13 || iter == 14) && trigVal[3][trigTypeCount[3] - 1])  trigVal[iter][iter2] = 0;
	
      }
    }
	*/

    Double_t maxTrkPt = -1;
    //    Double_t maxTrkEta = -100;
  
    for(int i = 0; i < nTrk; i++){
      if(fabs(trkEta[i]) > 2.4) continue;
      if(trkFake[i]) continue;
      if(trkPt[i] > maxTrkPt){
	maxTrkPt = trkPt[i];
	//	maxTrkEta = trkEta[i];
      }
    }

    
    Double_t maxMuPt = -1;
    //    Double_t maxMuEta = -100;

    for(int i = 0; i < nGen; i++){
      if(TMath::Abs(genPDG[i]) != 13) continue;
      if(fabs(genEta[i]) > 2.4) continue;
      if(genPt[i] > maxMuPt){
	maxMuPt = genPt[i];
	//	maxMuEta = genEta[i];
      }
    }
    

    Double_t maxDPt = -1;
    Double_t maxDEta = -100;

    for(int i = 0; i < nGen; i++){
      if(TMath::Abs(genPDG[i]) != 421) continue;
      if(fabs(genEta[i]) > 2.0) continue;
      if(genPt[i] > maxDPt){
	maxDPt = genPt[i];
	maxDEta = genEta[i];
      }
    }  


    Double_t maxPu3CaloAnaPt = -1;
    //Double_t maxPu3CaloAnaEta = -100;
    
    for(int i = 0; i < nPu3Caloref; ++i){
      if(fabs(jtPu3Caloeta[i]) > 2.0) continue;
      if(jtPu3Calopt[i] > maxPu3CaloAnaPt){
	maxPu3CaloAnaPt = jtPu3Calopt[i];
	//	maxPu3CaloAnaEta = jtPu3Caloeta[i];
      }
    }

    Double_t maxPu4CaloAnaPt = -1;
    //    Double_t maxPu4CaloAnaEta = -100;

    Double_t maxMidPu4CaloAnaPt = -1;
    //    Double_t maxMidPu4CaloAnaEta = -100;

    Double_t maxDiPu4CaloEta1p1AnaPt = -1;
    //    Double_t maxDiPu4CaloEta1p1AnaEta = -100;

    Double_t twoDiPu4CaloEta1p1AnaPt = -1;
    //    Double_t twoDiPu4CaloAnaEta = -100;

    Double_t maxDiPu4CaloEta0p7AnaPt = -1;
    //    Double_t maxDiPu4CaloEta0p7AnaEta = -100;

    Double_t twoDiPu4CaloEta0p7AnaPt = -1;
    //    Double_t twoDiPu4CaloAnaEta = -100;

    Double_t maxTriPu4CaloAnaPt = -1;
    //    Double_t maxTriPu4CaloAnaEta = -100;

    Double_t twoTriPu4CaloAnaPt = -1;
    //    Double_t twoTriPu4CaloAnaEta = -100;

    Double_t threeTriPu4CaloAnaPt = -1;
    //    Double_t threeTriPu4CaloAnaEta = -100;

    for(int i = 0; i < nPu4Caloref; ++i){
      if(fabs(jtPu4Caloeta[i]) > 2.0) continue;
      if(jtPu4Calopt[i] > maxPu4CaloAnaPt){
	maxPu4CaloAnaPt = jtPu4Calopt[i];
	//	maxPu4CaloAnaEta = jtPu4Caloeta[i];
      }

      if(fabs(jtPu4Caloeta[i]) < 2.0){
	if(jtPu4Calopt[i] > maxMidPu4CaloAnaPt){
	  maxMidPu4CaloAnaPt = jtPu4Calopt[i];
	  //	  maxMidPu4CaloAnaEta = jtPu4Caloeta[i];
	}
      }

      if(fabs(jtPu4Caloeta[i]) < 0.6){
	if(jtPu4Calopt[i] > maxDiPu4CaloEta0p7AnaPt){
	  twoDiPu4CaloEta0p7AnaPt = maxDiPu4CaloEta0p7AnaPt;
	  //	  twoDiPu4CaloAnaEta = maxDiPu4CaloAnaEta;

	  maxDiPu4CaloEta0p7AnaPt = jtPu4Calopt[i];
	  //	  maxDiPu4CaloEta0p7AnaEta = jtPu4Caloeta[i];
	}
	else if(jtPu4Calopt[i] > twoDiPu4CaloEta0p7AnaPt){
          twoDiPu4CaloEta0p7AnaPt = jtPu4Calopt[i];
	  //	  twoDiPu4CaloAnaEta = jtPu4Caloeta[i];
	}
      }

      if(fabs(jtPu4Caloeta[i]) < 1.0){
	if(jtPu4Calopt[i] > maxDiPu4CaloEta1p1AnaPt){
	  twoDiPu4CaloEta1p1AnaPt = maxDiPu4CaloEta1p1AnaPt;
	  //	  twoDiPu4CaloAnaEta = maxDiPu4CaloAnaEta;

	  maxDiPu4CaloEta1p1AnaPt = jtPu4Calopt[i];
	  //	  maxDiPu4CaloEta1p1AnaEta = jtPu4Caloeta[i];
	}
	else if(jtPu4Calopt[i] > twoDiPu4CaloEta1p1AnaPt){
          twoDiPu4CaloEta1p1AnaPt = jtPu4Calopt[i];
	  //	  twoDiPu4CaloAnaEta = jtPu4Caloeta[i];
	}
      }

      if(jtPu4Calopt[i] > maxTriPu4CaloAnaPt){
	threeTriPu4CaloAnaPt = twoTriPu4CaloAnaPt;
	//	threeTriPu4CaloAnaEta = twoTriPu4CaloAnaEta;

	twoTriPu4CaloAnaPt = maxTriPu4CaloAnaPt;
	//	twoTriPu4CaloAnaEta = maxTriPu4CaloAnaEta;
	
	maxTriPu4CaloAnaPt = jtPu4Calopt[i];
	//       	maxTriPu4CaloAnaEta = jtPu4Caloeta[i];
      }
      else if(jtPu4Calopt[i] > twoTriPu4CaloAnaPt){
        threeTriPu4CaloAnaPt = twoTriPu4CaloAnaPt;
	//	threeTriPu4CaloAnaEta = twoTriPu4CaloAnaEta;

	twoTriPu4CaloAnaPt = jtPu4Calopt[i];
	//	twoTriPu4CaloAnaEta = jtPu4Caloeta[i];
      }
      else if(jtPu4Calopt[i] > threeTriPu4CaloAnaPt){
        threeTriPu4CaloAnaPt = jtPu4Calopt[i];
	//	threeTriPu4CaloAnaEta = jtPu4Caloeta[i];                                                  
      }
    }

    Double_t maxPu3PFAnaPt = -1;
    //    Double_t maxPu3PFAnaEta = -100;
    
    for(int i = 0; i < nPu3PFref; ++i){
      if(fabs(jtPu3PFeta[i]) > 2.0) continue;
      if(jtPu3PFpt[i] > maxPu3PFAnaPt){
	maxPu3PFAnaPt = jtPu3PFpt[i];
	//	maxPu3PFAnaEta = jtPu3PFeta[i];
      }
    }

    Double_t maxPu4PFAnaPt = -1;
    //    Double_t maxPu4PFAnaEta = -100;

    for(int i = 0; i < nPu4PFref; ++i){
      if(fabs(jtPu4PFeta[i]) > 2.0) continue;
      if(jtPu4PFpt[i] > maxPu4PFAnaPt){
	maxPu4PFAnaPt = jtPu4PFpt[i];
	//	maxPu4PFAnaEta = jtPu4PFeta[i];
      }
    }

    /*
    Double_t maxVs3CaloAnaPt = -1;
    Double_t maxVs3CaloAnaEta = -100;
    
    for(int i = 0; i < nVs3Caloref; ++i){
      if(fabs(jtVs3Caloeta[i]) > 2.0) continue;
      if(jtVs3Calopt[i] > maxVs3CaloAnaPt){
	maxVs3CaloAnaPt = jtVs3Calopt[i];
	maxVs3CaloAnaEta = jtVs3Caloeta[i];
      }
    }

    Double_t maxVs4CaloAnaPt = -1;
    Double_t maxVs4CaloAnaEta = -100;

    for(int i = 0; i < nVs4Caloref; ++i){
      if(fabs(jtVs4Caloeta[i]) > 2.0) continue;
      if(jtVs4Calopt[i] > maxVs4CaloAnaPt){
	maxVs4CaloAnaPt = jtVs4Calopt[i];
	maxVs4CaloAnaEta = jtVs4Caloeta[i];
      }
    }
    */

    Double_t maxPhotonAnaPt = -1;
    //    Double_t maxPhotonAnaEta = -100;

    Double_t maxPhotonMidAnaPt = -1;
    //    Double_t maxPhotonMidAnaEta = -100;

    Double_t maxElectronAnaPt = -1;
    Double_t maxElectronAnaEta = -100;

    Double_t twoElectronAnaPt = -1;
    Double_t twoElectronAnaEta = -100;

    for(int i = 0; i < nPhotons; ++i){
      if(fabs(photonEta[i]) > 2.1) continue;
      if(TMath::Abs(seedTime[i]) > 3) continue;
      if(swissCrx[i] > 0.9) continue;
      if(sigmaIetaIeta[i] < 0.002) continue;
      if(sigmaIphiIphi[i] < 0.002) continue;
      //      if(isEle[i]) continue;
      if(photonPt[i] > maxPhotonAnaPt){	

	maxPhotonAnaPt = photonPt[i];
	//	maxPhotonAnaEta = photonEta[i];

	if(TMath::Abs(photonEta[i]) < 2.0 && photonPt[i] > maxPhotonMidAnaPt){
	  maxPhotonMidAnaPt = photonPt[i];
	}
      }
    }

    for(int i = 0; i < nPhotons; ++i){
      if(fabs(photonEta[i]) > 1.44) continue;
      if(!isEle[i]) continue;
      if(photonPt[i] > maxElectronAnaPt){	
	twoElectronAnaPt = maxElectronAnaPt;
	//	twoElectronAnaEta = maxElectronAnaEta;

	maxElectronAnaPt = photonPt[i];
	//	maxElectronAnaEta = photonEta[i];
      }
      else if(photonPt[i] > twoElectronAnaPt){
        twoElectronAnaPt = photonPt[i];
	//	twoElectronAnaEta = photonEta[i];
      }
    }
    
    maxElectronAnaPt = -1;
    twoElectronAnaPt = -1;

    //    maxElectronAnaEta = -100;
    //    twoElectronAnaEta = -100;

    for(int i = 0; i < nGen; i++){
      if(genStatus[i] != 1) continue;
      if(TMath::Abs(genPDG[i]) != 11) continue;
      //      if(fabs(genEta[i]) > 3.0) continue;
      if(genPt[i] > maxElectronAnaPt){
	twoElectronAnaPt = maxElectronAnaPt;
	//     	twoElectronAnaEta = maxElectronAnaEta;

	maxElectronAnaPt = genPt[i];
	//     	maxElectronAnaEta = genEta[i];
      }
      else if(genPt[i] > twoElectronAnaPt){
        twoElectronAnaPt = genPt[i];
	//	twoElectronAnaEta = genEta[i];
      }
    }
    /*
    Double_t trueZAnaPt = -1;
    //    Double_t trueZAnaEta = -1;
    
    for(int i = 0; i < nGen; i++){
      if(TMath::Abs(genPDG[i]) == 23){
		trueZAnaPt = genPt[i];
	//	trueZAnaEta = genEta[i];
      }
    }
    */
    Double_t maxMidPhotonAnaPt = -1;
    //    Double_t maxMidPhotonAnaEta = -100;

    for(int i = 0; i < nPhotons; ++i){
      if(fabs(photonEta[i]) > 1.24) continue;
      if(TMath::Abs(seedTime[i]) > 3) continue;
      if(swissCrx[i] > 0.9) continue;
      if(sigmaIetaIeta[i] < 0.002) continue;
      if(sigmaIphiIphi[i] < 0.002) continue;
      if(photonPt[i] > maxMidPhotonAnaPt){
	maxMidPhotonAnaPt = photonPt[i];
	//	maxMidPhotonAnaEta = photonEta[i];
      }
    }

    //FOR ADDITIONAL OFFLINE OBJECT MATCHING, EDIT HERE (2 of 2)
    //    const Int_t nOfflineObj = 4;
    const Int_t nOfflineObj = 1;
    //current
    //    const Int_t nOfflineObj = 29;
    //    const Int_t nOfflineObj = 27;
    //    const Int_t nOfflineObj = 23;
    if(nOfflineObj != nTrigType){
      std::cout << "ERROR: OFFLINE OBJECT NUMBER MUST MATCH NUMBER OF TRIGGER 'TYPES' IN INPUT TEXT FILE; RETURN 1" << std::endl;
      std::cout << nTrigType << std::endl;
      return 1;
    }
    /*
    Double_t trigOfflinePt[nOfflineObj] = {maxPu4PFAnaPt, maxPu4PFAnaPt, maxPu4PFAnaPt, maxPu4PFAnaPt};
    Double_t trigOfflineEta[nOfflineObj] = {maxPu4PFAnaEta, maxPu4PFAnaEta, maxPu4PFAnaEta, maxPu4PFAnaEta};
    Bool_t trigCond[nOfflineObj] = {true, true, true, true};
    */
    Double_t trigOfflinePt[nOfflineObj] = {maxDPt};
    Double_t trigOfflineEta[nOfflineObj] = {maxDEta};
    Bool_t trigCond[nOfflineObj] = {true};

    //    Double_t trigOfflinePt[nOfflineObj] = {maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPhotonMidAnaPt, maxPhotonMidAnaPt, maxMuPt, maxPu4CaloAnaPt, maxPhotonMidAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPhotonMidAnaPt, maxPhotonMidAnaPt, maxPhotonAnaPt, maxPhotonAnaPt, maxDiPu4CaloEta1p1AnaPt, maxDiPu4CaloEta0p7AnaPt, maxTriPu4CaloAnaPt, trueZAnaPt, trueZAnaPt, trueZAnaPt, trueZAnaPt, trueZAnaPt, trueZAnaPt, maxPu4PFAnaPt, maxMuPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPhotonMidAnaPt};

    //    Double_t trigOfflineEta[nOfflineObj] = {maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPhotonMidAnaEta, maxPhotonMidAnaEta, maxMuEta, maxPu4CaloAnaEta, maxPhotonMidAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPhotonMidAnaEta, maxPhotonMidAnaEta, maxPhotonAnaEta, maxPhotonAnaEta, maxDiPu4CaloEta1p1AnaEta, maxDiPu4CaloEta0p7AnaEta, maxTriPu4CaloAnaEta, trueZAnaEta, trueZAnaEta, trueZAnaEta, trueZAnaEta, trueZAnaEta, trueZAnaEta, maxPu4PFAnaEta, maxMuEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPhotonMidAnaEta};

    //    Bool_t trigCond[nOfflineObj] = {true, hiBin > 100, hiBin < 100 && hiBin > 60, true, hiBin > 100, true, maxMuPt > 10, maxMuPt > 10, true, hiBin > 100, hiBin < 100 && hiBin > 60, true, hiBin > 100, true, hiBin > 100, twoDiPu4CaloEta1p1AnaPt > 50, twoDiPu4CaloEta0p7AnaPt > 50, twoTriPu4CaloAnaPt > 60 && threeTriPu4CaloAnaPt > 60, TMath::Abs(maxElectronAnaEta) < 1.44 && TMath::Abs(twoElectronAnaEta) < 1.44, TMath::Abs(maxElectronAnaEta) < 1.44 && TMath::Abs(twoElectronAnaEta) < 1.44, TMath::Abs(maxElectronAnaEta) < 1.44 && TMath::Abs(twoElectronAnaEta) < 1.44, TMath::Abs(maxElectronAnaEta) < 2.0 && TMath::Abs(twoElectronAnaEta) < 2.0, TMath::Abs(maxElectronAnaEta) < 2.0 && TMath::Abs(twoElectronAnaEta) < 2.0, TMath::Abs(maxElectronAnaEta) < 2.0 && TMath::Abs(twoElectronAnaEta) < 2.0, true, true, maxMuPt > 10, maxPhotonMidAnaPt > 20, maxMuPt > 10};
    /*
    Double_t trigOfflinePt[nOfflineObj] = {maxTrkPt, maxPu3CaloAnaPt, maxPu3CaloAnaPt,  maxPu4CaloAnaPt, maxMidPu4CaloAnaPt, maxPu4PFAnaPt, maxPu4CaloAnaPt, maxPhotonAnaPt, maxMidPhotonAnaPt, maxDiPu4CaloAnaPt, maxDiPu4CaloAnaPt, maxTriPu4CaloAnaPt, maxDiPu4CaloAnaPt, maxDiPu4CaloAnaPt, maxTriPu4CaloAnaPt, maxMuPt, maxMuPt, maxMuPt, maxMuPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt};
    Double_t trigOfflineEta[nOfflineObj] = {maxTrkEta, maxPu3CaloAnaEta, maxPu3CaloAnaEta, maxPu4CaloAnaEta, maxMidPu4CaloAnaEta, maxPu4PFAnaEta, maxPu4CaloAnaEta, maxPhotonAnaEta, maxMidPhotonAnaEta, maxDiPu4CaloAnaEta, maxDiPu4CaloAnaEta, maxTriPu4CaloAnaEta, maxDiPu4CaloAnaEta, maxDiPu4CaloAnaEta, maxTriPu4CaloAnaEta, maxMuEta, maxMuEta, maxMuEta, maxMuEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta};
    Bool_t trigCond[nOfflineObj] = {true, true, true, true, true, true, true, true, true, twoDiPu4CaloAnaPt > 50.0, twoDiPu4CaloAnaPt > 50.0, twoTriPu4CaloAnaPt > 65.0 && threeTriPu4CaloAnaPt > 65.0, twoDiPu4CaloAnaPt > 50.0, twoDiPu4CaloAnaPt > 50.0, twoTriPu4CaloAnaPt > 65.0 && threeTriPu4CaloAnaPt > 65.0, true, true, true, true, maxMuPt > 3, maxMuPt > 3, maxMuPt > 3, maxMuPt > 3, maxMuPt > 10, maxMuPt > 10, maxMuPt > 10, maxMuPt > 10};
    */
    //    Double_t trigOfflinePt[nOfflineObj] = {maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt};
    //    Double_t trigOfflineEta[nOfflineObj] = {maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta};
    //    Bool_t trigCond[nOfflineObj] = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};

    for(Int_t iter = 0; iter < nTrigType; iter++){
      if(trigOfflinePt[iter] > 0 && trigCond[iter]){
	histsPt_p[iter][0]->Fill(trigOfflinePt[iter]);

	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  if(trigOfflinePt[iter] > trigThresh[iter][iter2]) histsEta_p[iter][iter2]->Fill(trigOfflineEta[iter]);

	  if(trigVal[iter][iter2]){
	    histsPt_p[iter][iter2+1]->Fill(trigOfflinePt[iter]);
	    if(trigOfflinePt[iter] > trigThresh[iter][iter2]) histsEtaTrig_p[iter][iter2]->Fill(trigOfflineEta[iter]);
	  }
	  else{
	    if(highPtMiss[iter][iter2] < trigOfflinePt[iter]) highPtMiss[iter][iter2] = trigOfflinePt[iter];
	  }
	}
      }
      
      if(TMath::Abs(maxElectronAnaEta) < 1.44 && TMath::Abs(twoElectronAnaEta) < 1.44){
	for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
	  nZTrigFireDenom[iter][iter2]++;	  

	  if(trigVal[iter][iter2]) nZTrigFireNum[iter][iter2]++;
	}
      }
      
    }
  }

  std::cout << std::endl;
  std::cout << "Matched events: " << matched << std::endl;
  std::cout << "Trigger fires: " << std::endl;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    std::cout << "  Trigger Type, raw #, rate matched (Hz), rate unmatched (Hz), highestMiss: " << trigType[iter] << std::endl;
    
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      std::cout << "    " << trigName[iter][iter2] << ": " << nTrigFire[iter][iter2] << ", " << std::setprecision(5) << nTrigFire[iter][iter2]*30000./matched << ", " << std::setprecision(5) << nTrigFire[iter][iter2]*30000./HLTTree->GetEntries() << ", " << std::setprecision(5) << highPtMiss[iter][iter2] << std::endl;

      Int_t bin = ratesMatched_p[iter]->FindBin(trigThresh[iter][iter2]);
      Float_t num = nTrigFire[iter][iter2]*30000.;
      Float_t numErr = (nTrigFire[iter][iter2] + TMath::Sqrt(nTrigFire[iter][iter2]))*30000.;

      ratesMatched_p[iter]->SetBinContent(bin, num/matched);
      ratesMatched_p[iter]->SetBinError(bin, numErr/matched - num/matched);
      ratesUnmatched_p[iter]->SetBinContent(bin, num/HLTTree->GetEntries());
      ratesUnmatched_p[iter]->SetBinError(bin, numErr/HLTTree->GetEntries() - num/HLTTree->GetEntries());
    }
  }

  std::cout << std::endl;
  std::cout << "Z eff: " << std::endl;
  for(Int_t iter = 0; iter < nTrigType; iter++){
    std::cout << "  Trigger Type, raw #, denom #, eff" << trigType[iter] << std::endl;

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      std::cout << "    " << trigName[iter][iter2] << ": " << nZTrigFireNum[iter][iter2] << ", " << nZTrigFireDenom[iter][iter2] << ", " << std::setprecision(5) << nZTrigFireNum[iter][iter2]/nZTrigFireDenom[iter][iter2] << std::endl;
    }
  }

  TGraphAsymmErrors *aPt_p[nTrigType][maxNTrig2], *aEta_p[nTrigType][maxNTrig2];

  for(Int_t iter = 0; iter < nTrigType; iter++){
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      //      if(iter2 == 0) histsSpectPt_p[iter] = (TH1F*)(histsPt_p[iter][iter2+1]->Clone());
      histsSpectPt_p[iter]->Add(histsPt_p[iter][iter2+1]);

      aPt_p[iter][iter2] = new TGraphAsymmErrors();
      aPt_p[iter][iter2]->BayesDivide(histsPt_p[iter][iter2+1],histsPt_p[iter][0]);
      aPt_p[iter][iter2]->SetName(Form("%s_%s_pt_asymm", trigName[iter][iter2].c_str(), trigType[iter].c_str()));

      aEta_p[iter][iter2] = new TGraphAsymmErrors();
      aEta_p[iter][iter2]->BayesDivide(histsEtaTrig_p[iter][iter2],histsEta_p[iter][iter2]);
      aEta_p[iter][iter2]->SetName(Form("%s_%s_eta_asymm", trigName[iter][iter2].c_str(), trigType[iter].c_str()));
    }
  }

  std::string outName = inHLTFile;
  const std::string inString = ".root";
  const std::string outString = "_HIST.root";
  std::size_t strIndex = 0;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString); 
  }

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");

  for(Int_t iter = 0; iter < nTrigType; iter++){
    histsPt_p[iter][0]->Write("", TObject::kOverwrite);
    histsSpectPt_p[iter]->Write("", TObject::kOverwrite);

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      histsPt_p[iter][iter2+1]->Write("", TObject::kOverwrite);
      histsEta_p[iter][iter2]->Write("", TObject::kOverwrite);
      histsEtaTrig_p[iter][iter2]->Write("", TObject::kOverwrite);

      aPt_p[iter][iter2]->Write("", TObject::kOverwrite);
      aEta_p[iter][iter2]->Write("", TObject::kOverwrite);
    }
    ratesMatched_p[iter]->Write("", TObject::kOverwrite);
    ratesUnmatched_p[iter]->Write("", TObject::kOverwrite);
  }

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    delete histsPt_p[iter][0];
    delete histsSpectPt_p[iter];

    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      delete histsEta_p[iter][iter2];
      delete histsEtaTrig_p[iter][iter2];

      delete aPt_p[iter][iter2];
      delete aEta_p[iter][iter2];
    }
    delete ratesMatched_p[iter];
    delete ratesUnmatched_p[iter];
  }

  AnaFile_p->Close();
  delete AnaFile_p;

  HLTFile_p->Close();
  delete HLTFile_p;

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 4){
    std::cout << "Usage: matchTrigTree_HI <inHLTFile> <inForestFile> <inTrigFileName>" << std::endl;
    std::cout << "argNum: " << argc << std::endl;
    for(Int_t iter = 0; iter < argc; iter++){
      std::cout << "arg " << iter << ": " << argv[iter] << std::endl;
    }

    return -1;
  }

  int rStatus = -1;

  rStatus = matchTrigTree_HI(argv[1], argv[2], argv[3]);

  return rStatus;
}
