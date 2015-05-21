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

const TString Ana3CaloTreename = "akPu3CaloJetAnalyzer/t";
const TString Ana4CaloTreename = "akPu4CaloJetAnalyzer/t";
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

  //FOR ADDITIONAL OFFLINE OBJECT MATCHING, EDIT HERE (1 of 2)

  const Int_t nPtBins[nTrigType] = {100, 200, 200, 200, 200, 200, 100, 100, 200, 200, 200, 200};
  const Int_t maxPt[nTrigType] = {100, 200, 200, 200, 200, 200, 100, 100, 200, 200, 200, 200};
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

    Int_t nBins = (max+minDiff - (min-minDiff))/(minDiff*2);
    if(trigTypeCount[iter] == 1){
      nBins = 1;
      min = trigThresh[iter][0] - 20;
      max = trigThresh[iter][0] + 20;
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
	//EDIT HERE FOR ANDING
	if(iter >= nTrigType-2) trigVal[iter][iter2] = (trigVal[iter][trigTypeCount[iter]-1] && trigVal[2][iter2]);

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

    Double_t maxDi4CaloAnaPt = -1;
    Double_t maxDi4CaloAnaEta = -100;

    Double_t twoDi4CaloAnaPt = -1;
    //    Double_t twoDi4CaloAnaEta = -100;

    Double_t maxTri4CaloAnaPt = -1;
    Double_t maxTri4CaloAnaEta = -100;

    Double_t twoTri4CaloAnaPt = -1;
    //    Double_t twoTri4CaloAnaEta = -100;

    Double_t threeTri4CaloAnaPt = -1;
    //    Double_t threeTri4CaloAnaEta = -100;

    for(int i = 0; i < n4Caloref; ++i){
      if(fabs(jt4Caloeta[i]) > 5.0) continue;
      if(jt4Calopt[i] > max4CaloAnaPt){
	max4CaloAnaPt = jt4Calopt[i];
	max4CaloAnaEta = jt4Caloeta[i];
      }

      if(fabs(jt4Caloeta[i]) < 0.5){
	if(jt4Calopt[i] > maxDi4CaloAnaPt){
	  twoDi4CaloAnaPt = maxDi4CaloAnaPt;
	  //	  twoDi4CaloAnaEta = maxDi4CaloAnaEta;

	  maxDi4CaloAnaPt = jt4Calopt[i];
	  maxDi4CaloAnaEta = jt4Caloeta[i];
	}
	else if(jt4Calopt[i] > twoDi4CaloAnaPt){
          twoDi4CaloAnaPt = jt4Calopt[i];
	  //          twoDi4CaloAnaEta = jt4Caloeta[i];
	}
      }

      if(jt4Calopt[i] > maxTri4CaloAnaPt){
	threeTri4CaloAnaPt = twoTri4CaloAnaPt;
	//	  threeTri4CaloAnaEta = twoTri4CaloAnaEta;

	twoTri4CaloAnaPt = maxTri4CaloAnaPt;
	//	  twoTri4CaloAnaEta = maxTri4CaloAnaEta;
	
	maxTri4CaloAnaPt = jt4Calopt[i];
	maxTri4CaloAnaEta = jt4Caloeta[i];
      }
      else if(jt4Calopt[i] > twoTri4CaloAnaPt){
        threeTri4CaloAnaPt = twoTri4CaloAnaPt;
	//          threeTri4CaloAnaEta = twoTri4CaloAnaEta;

	twoTri4CaloAnaPt = jt4Calopt[i];
	//          twoTri4CaloAnaEta = jt4Caloeta[i];
      }
      else if(jt4Calopt[i] > threeTri4CaloAnaPt){
        threeTri4CaloAnaPt = jt4Calopt[i];
        //          threeTri4CaloAnaEta = jt4Caloeta[i];                                                  
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

    //FOR ADDITIONAL OFFLINE OBJECT MATCHING, EDIT HERE (2 of 2)
    const Int_t nOfflineObj = 12;
    if(nOfflineObj != nTrigType){
      std::cout << "ERROR: OFFLINE OBJECT NUMBER MUST MATCH NUMBER OF TRIGGER 'TYPES' IN INPUT TEXT FILE; RETURN 1" << std::endl;
      return 1;
    }

    Double_t trigOfflinePt[nOfflineObj] = {maxTrkPt, max3CaloAnaPt, max3CaloAnaPt, max4CaloAnaPt, max4CaloAnaPt, max4CaloAnaPt, maxPhotonAnaPt, maxPhotonAnaPt, maxDi4CaloAnaPt, maxTri4CaloAnaPt, maxMuPt, maxMuPt};
    Double_t trigOfflineEta[nOfflineObj] = {maxTrkEta, max3CaloAnaEta, max3CaloAnaEta, max4CaloAnaEta, max4CaloAnaEta, max4CaloAnaEta, maxPhotonAnaEta, maxPhotonAnaEta, maxDi4CaloAnaEta, maxTri4CaloAnaEta, maxMuEta, maxMuEta};
    Bool_t trigCond[nOfflineObj] = {true, true, true, true, true, true, true, true, twoDi4CaloAnaPt > 55.0, twoTri4CaloAnaPt > 65.0 && threeTri4CaloAnaPt > 65.0, true, true};

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
