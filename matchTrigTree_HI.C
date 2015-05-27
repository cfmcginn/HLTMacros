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
      trigVal[nLines][tempPosIter] = 0;
      nTrigFire[nLines][tempPosIter] = 0;
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
  for(Int_t iter = 0; iter < nTrigType-8; iter++){
    if(iter == 11 || iter == 12 || iter == 13) continue;
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      HLTTree->SetBranchStatus(trigName[iter][iter2].c_str(), 1);
      HLTTree->SetBranchAddress(trigName[iter][iter2].c_str(), &(trigVal[iter][iter2]));
    }
  }

  TFile *AnaFile_p = new TFile(inForestFile.c_str(), "READ");
  TTree *AnaPu3CaloTree = (TTree*)AnaFile_p->Get(AnaPu3CaloTreename); 
  TTree *AnaPu4CaloTree = (TTree*)AnaFile_p->Get(AnaPu4CaloTreename); 
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
  Int_t genPDG[maxGen];

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

  const Int_t maxGamma = 50;
  Int_t nPhotons;
  Float_t photonPt[maxGamma], photonEta[maxGamma], photonPhi[maxGamma];   //[nPhotons]
  Float_t seedTime[maxGamma], swissCrx[maxGamma], sigmaIphiIphi[maxGamma], sigmaIetaIeta[maxGamma];

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
  AnaGenTree->SetBranchAddress("mult", &nGen);
  AnaGenTree->SetBranchAddress("pt", genPt);
  AnaGenTree->SetBranchAddress("phi", genPhi);
  AnaGenTree->SetBranchAddress("eta", genEta);
  AnaGenTree->SetBranchAddress("pdg", genPDG);

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

  AnaPhotonTree->SetBranchStatus("*", 0);
  AnaPhotonTree->SetBranchStatus("nPhotons", 1);
  AnaPhotonTree->SetBranchStatus("pt", 1);
  AnaPhotonTree->SetBranchStatus("eta", 1);
  AnaPhotonTree->SetBranchStatus("phi", 1);
  AnaPhotonTree->SetBranchStatus("seedTime", 1);
  AnaPhotonTree->SetBranchStatus("swissCrx", 1);
  AnaPhotonTree->SetBranchStatus("sigmaIphiIphi", 1);
  AnaPhotonTree->SetBranchStatus("sigmaIetaIeta", 1);
  
  AnaPhotonTree->SetBranchAddress("nPhotons", &nPhotons);
  AnaPhotonTree->SetBranchAddress("pt", photonPt);
  AnaPhotonTree->SetBranchAddress("eta", photonEta);
  AnaPhotonTree->SetBranchAddress("phi", photonPhi);
  AnaPhotonTree->SetBranchAddress("seedTime", seedTime);
  AnaPhotonTree->SetBranchAddress("swissCrx", swissCrx);
  AnaPhotonTree->SetBranchAddress("sigmaIphiIphi", sigmaIphiIphi);
  AnaPhotonTree->SetBranchAddress("sigmaIetaIeta", sigmaIetaIeta);

  //FOR ADDITIONAL OFFLINE OBJECT MATCHING, EDIT HERE (1 of 2)

  const Int_t nPtBins[nTrigType] = {100, 200, 200, 200, 200, 200, 100, 100, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};
  const Int_t maxPt[nTrigType] = {100, 200, 200, 200, 200, 200, 100, 100, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};
  //const Int_t nPtBins[nTrigType] = {120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120};
  //const Int_t maxPt[nTrigType] = {130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130};

  const Int_t nEtaBins = 50;
  TH1F *histsPt_p[nTrigType][maxNTrig2+1], *histsEta_p[nTrigType][maxNTrig2], *histsEtaTrig_p[nTrigType][maxNTrig2];

  TH1F* ratesMatched_p[nTrigType];
  TH1F* ratesUnmatched_p[nTrigType];

  for(Int_t iter = 0; iter < nTrigType; iter++){
    histsPt_p[iter][0] = new TH1F(Form("leading%s_pt", trigType[iter].c_str()), Form("leading%s_pt", trigType[iter].c_str()), nPtBins[iter], 0.0, maxPt[iter]);

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
	
	if(iter >= nTrigType-8 && iter < nTrigType-6) trigVal[iter][iter2] = (trigVal[iter-4][0] && trigVal[3][iter2]);
	else if(iter >= nTrigType-6 && iter < nTrigType-4) trigVal[iter][iter2] = (trigVal[iter-2][0] && trigVal[3][iter2]);
	else if(iter >= nTrigType-4 && iter < nTrigType-2) trigVal[iter][iter2] = (trigVal[iter-6][0] && trigVal[3][iter2]);
	else if(iter >= nTrigType-2) trigVal[iter][iter2] = (trigVal[iter-2][0] && trigVal[3][iter2]);

	if(iter >= nTrigType-6 && iter < nTrigType-4 && trigVal[3][trigTypeCount[3] - 1]) trigVal[iter][iter2] = 0;
	else if(iter >= nTrigType-2 && trigVal[3][trigTypeCount[3] - 1]) trigVal[iter][iter2] = 0;

	if(iter == 11 || iter == 12 || iter == 13) trigVal[iter][iter2] =  trigVal[iter-3][iter2];
	if((iter == 11 || iter == 12 || iter == 13) && trigVal[3][trigTypeCount[3] - 1])  trigVal[iter][iter2] = 0;
	

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
    AnaPhotonTree->GetEntry(entry);
    AnaHITree->GetEntry(entry);
    AnaSkimTree->GetEntry(entry);
    AnaTrkTree->GetEntry(entry);
    AnaGenTree->GetEntry(entry);

    long long hlt_entry = matcher->retrieveEvent(ana_event, ana_lumi, 0);
    if(hlt_entry == -1) continue;

    HLTTree->GetEntry(hlt_entry);
    matched++;
    
    for(Int_t iter = 0; iter < nTrigType; iter++){
      for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
        //EDIT HERE FOR ANDING                                                                       
        if(iter >= nTrigType-8 && iter < nTrigType-6) trigVal[iter][iter2] = (trigVal[iter-4][0] && trigVal[3][iter2]);
        else if(iter >= nTrigType-6 && iter < nTrigType-4) trigVal[iter][iter2] = (trigVal[iter-2][0] && trigVal[3][iter2]);
        else if(iter >= nTrigType-4 && iter < nTrigType-2) trigVal[iter][iter2] = (trigVal[iter-6][0]&& trigVal[3][iter2]);
        else if(iter >= nTrigType-2) trigVal[iter][iter2] = (trigVal[iter-2][0] && trigVal[3][iter2]);

        if(iter >= nTrigType-6 && iter < nTrigType-4 && trigVal[3][trigTypeCount[3] - 1]) trigVal[iter][iter2] = 0;
        else if(iter >= nTrigType-2 && trigVal[3][trigTypeCount[3] - 1]) trigVal[iter][iter2] = 0;

        if(iter == 11 || iter == 12 || iter == 13) trigVal[iter][iter2] =  trigVal[iter-3][iter2];
        if((iter == 11 || iter == 12 || iter == 13) && trigVal[3][trigTypeCount[3] - 1])  trigVal[iter][iter2] = 0;
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

    Double_t maxMuPt = -1;
    Double_t maxMuEta = -100;

    for(int i = 0; i < nGen; i++){
      if(TMath::Abs(genPDG[i]) != 13) continue;
      if(fabs(genEta[i]) > 2.4) continue;
      if(genPt[i] > maxMuPt){
	maxMuPt = genPt[i];
       	maxMuEta = genEta[i];
      }
    }

  
    Double_t maxPu3CaloAnaPt = -1;
    Double_t maxPu3CaloAnaEta = -100;
    
    for(int i = 0; i < nPu3Caloref; ++i){
      if(fabs(jtPu3Caloeta[i]) > 2.0) continue;
      if(jtPu3Calopt[i] > maxPu3CaloAnaPt){
	maxPu3CaloAnaPt = jtPu3Calopt[i];
	maxPu3CaloAnaEta = jtPu3Caloeta[i];
      }
    }

    Double_t maxPu4CaloAnaPt = -1;
    Double_t maxPu4CaloAnaEta = -100;

    Double_t maxMidPu4CaloAnaPt = -1;
    Double_t maxMidPu4CaloAnaEta = -100;

    Double_t maxDiPu4CaloAnaPt = -1;
    Double_t maxDiPu4CaloAnaEta = -100;

    Double_t twoDiPu4CaloAnaPt = -1;
    //    Double_t twoDiPu4CaloAnaEta = -100;

    Double_t maxTriPu4CaloAnaPt = -1;
    Double_t maxTriPu4CaloAnaEta = -100;

    Double_t twoTriPu4CaloAnaPt = -1;
    //    Double_t twoTriPu4CaloAnaEta = -100;

    Double_t threeTriPu4CaloAnaPt = -1;
    //    Double_t threeTriPu4CaloAnaEta = -100;

    for(int i = 0; i < nPu4Caloref; ++i){
      if(fabs(jtPu4Caloeta[i]) > 5.0) continue;
      if(jtPu4Calopt[i] > maxPu4CaloAnaPt){
	maxPu4CaloAnaPt = jtPu4Calopt[i];
	maxPu4CaloAnaEta = jtPu4Caloeta[i];
      }

      if(fabs(jtPu4Caloeta[i]) < 2.0){
	if(jtPu4Calopt[i] > maxMidPu4CaloAnaPt){
	  maxMidPu4CaloAnaPt = jtPu4Calopt[i];
	  maxMidPu4CaloAnaEta = jtPu4Caloeta[i];
	}
      }

      if(fabs(jtPu4Caloeta[i]) < 0.6){
	if(jtPu4Calopt[i] > maxDiPu4CaloAnaPt){
	  twoDiPu4CaloAnaPt = maxDiPu4CaloAnaPt;
	  //	  twoDiPu4CaloAnaEta = maxDiPu4CaloAnaEta;

	  maxDiPu4CaloAnaPt = jtPu4Calopt[i];
	  maxDiPu4CaloAnaEta = jtPu4Caloeta[i];
	}
	else if(jtPu4Calopt[i] > twoDiPu4CaloAnaPt){
          twoDiPu4CaloAnaPt = jtPu4Calopt[i];
	  //	  twoDiPu4CaloAnaEta = jtPu4Caloeta[i];
	}
      }

      if(jtPu4Calopt[i] > maxTriPu4CaloAnaPt){
	threeTriPu4CaloAnaPt = twoTriPu4CaloAnaPt;
	//	threeTriPu4CaloAnaEta = twoTriPu4CaloAnaEta;

	twoTriPu4CaloAnaPt = maxTriPu4CaloAnaPt;
	//	twoTriPu4CaloAnaEta = maxTriPu4CaloAnaEta;
	
	maxTriPu4CaloAnaPt = jtPu4Calopt[i];
	maxTriPu4CaloAnaEta = jtPu4Caloeta[i];
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

    Double_t maxPhotonAnaPt = -1;
    Double_t maxPhotonAnaEta = -100;

    for(int i = 0; i < nPhotons; ++i){
      if(fabs(photonEta[i]) > 1.8) continue;
      if(TMath::Abs(seedTime[i]) > 3) continue;
      if(swissCrx[i] > 0.9) continue;
      if(sigmaIetaIeta[i] < 0.002) continue;
      if(sigmaIphiIphi[i] < 0.002) continue;
      if(photonPt[i] > maxPhotonAnaPt){	

	maxPhotonAnaPt = photonPt[i];
	maxPhotonAnaEta = photonEta[i];
      }
    }

    Double_t maxMidPhotonAnaPt = -1;
    Double_t maxMidPhotonAnaEta = -100;

    for(int i = 0; i < nPhotons; ++i){
      if(fabs(photonEta[i]) > 1.24) continue;
      if(TMath::Abs(seedTime[i]) > 3) continue;
      if(swissCrx[i] > 0.9) continue;
      if(sigmaIetaIeta[i] < 0.002) continue;
      if(sigmaIphiIphi[i] < 0.002) continue;
      if(photonPt[i] > maxMidPhotonAnaPt){
	maxMidPhotonAnaPt = photonPt[i];
	maxMidPhotonAnaEta = photonEta[i];
      }
    }

    //FOR ADDITIONAL OFFLINE OBJECT MATCHING, EDIT HERE (2 of 2)
    const Int_t nOfflineObj = 26;
    //    const Int_t nOfflineObj = 23;
    if(nOfflineObj != nTrigType){
      std::cout << "ERROR: OFFLINE OBJECT NUMBER MUST MATCH NUMBER OF TRIGGER 'TYPES' IN INPUT TEXT FILE; RETURN 1" << std::endl;
      return 1;
    }

    Double_t trigOfflinePt[nOfflineObj] = {maxTrkPt, maxPu3CaloAnaPt, maxPu3CaloAnaPt, maxPu4CaloAnaPt, maxMidPu4CaloAnaPt, maxPu4CaloAnaPt, maxPhotonAnaPt, maxMidPhotonAnaPt, maxDiPu4CaloAnaPt, maxDiPu4CaloAnaPt, maxTriPu4CaloAnaPt, maxDiPu4CaloAnaPt, maxDiPu4CaloAnaPt, maxTriPu4CaloAnaPt, maxMuPt, maxMuPt, maxMuPt, maxMuPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt, maxPu4CaloAnaPt};
    Double_t trigOfflineEta[nOfflineObj] = {maxTrkEta, maxPu3CaloAnaEta, maxPu3CaloAnaEta, maxPu4CaloAnaEta, maxMidPu4CaloAnaEta, maxPu4CaloAnaEta, maxPhotonAnaEta, maxMidPhotonAnaEta, maxDiPu4CaloAnaEta, maxDiPu4CaloAnaEta, maxTriPu4CaloAnaEta, maxDiPu4CaloAnaEta, maxDiPu4CaloAnaEta, maxTriPu4CaloAnaEta, maxMuEta, maxMuEta, maxMuEta, maxMuEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta, maxPu4CaloAnaEta};
    Bool_t trigCond[nOfflineObj] = {true, true, true, true, true, true, true, true, twoDiPu4CaloAnaPt > 50.0, twoDiPu4CaloAnaPt > 50.0, twoTriPu4CaloAnaPt > 65.0 && threeTriPu4CaloAnaPt > 65.0, twoDiPu4CaloAnaPt > 50.0, twoDiPu4CaloAnaPt > 50.0, twoTriPu4CaloAnaPt > 65.0 && threeTriPu4CaloAnaPt > 65.0, true, true, true, true, maxMuPt > 3, maxMuPt > 3, maxMuPt > 3, maxMuPt > 3, maxMuPt > 10, maxMuPt > 10, maxMuPt > 10, maxMuPt > 10};

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
	}
      }
    }
  }

  std::cout << std::endl;
  std::cout << "Matched events: " << matched << std::endl;
  std::cout << "Trigger fires: " << std::endl;

  for(Int_t iter = 0; iter < nTrigType; iter++){
    std::cout << "  Trigger Type, raw #, rate matched (Hz), rate unmatched (Hz): " << trigType[iter] << std::endl;
    
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
      std::cout << "    " << trigName[iter][iter2] << ": " << nTrigFire[iter][iter2] << ", " << std::setprecision(5) << nTrigFire[iter][iter2]*30000./matched << ", " << std::setprecision(5) << nTrigFire[iter][iter2]*30000./HLTTree->GetEntries() << std::endl;

      Int_t bin = ratesMatched_p[iter]->FindBin(trigThresh[iter][iter2]);
      Float_t num = nTrigFire[iter][iter2]*30000.;
      Float_t numErr = (nTrigFire[iter][iter2] + TMath::Sqrt(nTrigFire[iter][iter2]))*30000.;

      ratesMatched_p[iter]->SetBinContent(bin, num/matched);
      ratesMatched_p[iter]->SetBinError(bin, numErr/matched - num/matched);
      ratesUnmatched_p[iter]->SetBinContent(bin, num/HLTTree->GetEntries());
      ratesUnmatched_p[iter]->SetBinError(bin, numErr/HLTTree->GetEntries() - num/HLTTree->GetEntries());
    }
  }

  TGraphAsymmErrors *aPt_p[nTrigType][maxNTrig2], *aEta_p[nTrigType][maxNTrig2];

  for(Int_t iter = 0; iter < nTrigType; iter++){
    for(Int_t iter2 = 0; iter2 < trigTypeCount[iter]; iter2++){
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
