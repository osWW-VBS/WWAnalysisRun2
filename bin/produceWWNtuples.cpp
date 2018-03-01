#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include <TClonesArray.h>           

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"

// HEADER FILE FOR JES 
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// HEADER FOR B-TAG
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "../BtagUnc.hh"

#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"
#include "../interface/readJSONFile.h"
#include "../interface/Utils.hh"

using namespace std;
void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, 
		  TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22,
		  double& costheta1,  double& costheta2,  double& Phi,
		  double& costhetastar, double& Phi1);

double GetSFs_Lepton(double pt, double eta, TH1F* h1);
double GetMin(double x, double y);
double GetMax(double x, double y);
//float getPUPPIweight(float puppipt, float puppieta );
float getPUPPIweight(TF1* puppisd_corrGEN, TF1* puppisd_corrRECO_cen, TF1* puppisd_corrRECO_for, float puppipt, float puppieta );
double func( double pt, double eta, JetCorrectionUncertainty *fJetUnc); 
double GetPt_MET(double pfMET,  double phi, double pz);
double GetEta_MET(double pfMET, double phi, double pz);

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  
  int t0 = time(NULL);
  
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  int isMC = atoi(argv[3]);
  std::string cluster = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string TotalNumberOfEntries = argv[8];
  std::string TotalNumberOfNegativeEntries = argv[9];
  int applyTrigger = atoi(argv[10]);
  std::string jsonFileName = argv[11];
  int isLocal = atoi(argv[12]);
  int VBFSel  = atoi(argv[13]);
  
  std::string leptonName;

  if ( VBFSel==1)	cout<<"==> VBF selection method : Select two highest pT jets"<<endl;
  else if ( VBFSel==2)	cout<<"==> VBF selection method : Select pair with highest mjj..."<<endl;
  else if ( VBFSel==3)	cout<<"==> VBF selection method : Select pair with highest DeltaEta..."<<endl;
  else {	cout<<"\n\nERROR:	Enter valid vbf selection criteria....\n\n"<<endl;
  		exit(0);  
	}
  
  //applyTrigger=true;
  std::string iHLTFile="${CMSSW_BASE}/src/BaconAna/DataFormats/data/HLTFile_25ns";
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = iHLTFile.find(cmssw_base_env);
  if(start_pos != std::string::npos) {
  	iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }

  const baconhep::TTrigger triggerMenu(iHLTFile);  
  std::cout<<"Apply trigger: "<<applyTrigger<<std::endl;

  TLorentzVector W_type0,W_type0_jes_up, W_type0_jes_dn, W_type0_jer_up, W_type0_jer_dn, W_type2, W_run2,W_puppi_type2, W_puppi_type0, W_puppi_run2, LEP1, LEP2, SJ1_PuppiAK8, SJ2_PuppiAK8, SJ1, SJ2;
  TLorentzVector NU0,NU1,NU2,NU0_puppi,NU1_puppi,NU2_puppi;
  TLorentzVector NU0_jes_up, NU0_jes_dn;
  TLorentzVector NU0_puppi_jes_up, NU0_puppi_jes_dn;
  TLorentzVector NU0_jer_up, NU0_jer_dn;
  TLorentzVector JET, JET_PuppiAK8, AK4;
  TLorentzVector JET_jes_up, JET_jes_dn, JET_PuppiAK8_jes_up, JET_PuppiAK8_jes_dn;
  TLorentzVector AK4_JET1,AK4_JET2;
  TLorentzVector AK4_JET1_jes_up, AK4_JET1_jes_dn;
  TLorentzVector AK4_JET2_jes_up, AK4_JET2_jes_dn;
  TLorentzVector PuppiAK4_JET1,PuppiAK4_JET2;
  TLorentzVector PuppiAK4_JET1_jes_up, PuppiAK4_JET1_jes_dn;
  TLorentzVector PuppiAK4_JET2_jes_up, PuppiAK4_JET2_jes_dn;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector VBF1_jes_up, VBF1_jes_dn, VBF2_jes_up, VBF2_jes_dn;
  TLorentzVector ELE,MU;

  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;  

  int ok=0, total=0;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  	= new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen	= new baconhep::TGenEventInfo();
  TClonesArray *muonArr    	= new TClonesArray("baconhep::TMuon");
  TClonesArray *electronArr	= new TClonesArray("baconhep::TElectron");
  TClonesArray *vertexArr	= new TClonesArray("baconhep::TVertex");
  TClonesArray *jetArr		= new TClonesArray("baconhep::TJet");
  TClonesArray *vjetArrPuppi	= new TClonesArray("baconhep::TJet");
  TClonesArray *vjetAddArrPuppi	= new TClonesArray("baconhep::TAddJet");
  TClonesArray *lheWgtArr	= new TClonesArray("baconhep::TLHEWeight");
  

  char command1[3000];
  if ( cluster == "lxplus")
  	sprintf(command1, "eos find -f %s  | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://eoscms.cern.ch/\"$1}' > listTemp_%s.txt", (inputFolder).c_str(), outputFile.c_str());	// NF in awk command skips the blank line
  else 
	sprintf(command1,"xrdfs root://cmseos.fnal.gov ls %s | awk '{print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());
	//sprintf(command1,"eos root://cmseos.fnal.gov find -f %s | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());	// WORKS ONLY WITH INTERACTIVE NODE

  std::cout<<command1<<std::endl;
  system(command1);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);
  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);

  int fileCounter=0;

  vector<TString>  sampleName; 

  while (!rootList.eof())
  {
	char iRun_tW[700];
	rootList >> iRun_tW;
	if(!rootList.good())break;
	sampleName.push_back(iRun_tW);
	fileCounter++;
  }

  TFile *infile=0;
  TTree *eventTree=0;
  int cutEff[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  // Read and add pileup in histogram
  TFile* pileupFileMC = TFile::Open("puWeights_80x_37ifb.root");
  TH1D* puWeights = (TH1D*)pileupFileMC->Get("puWeights");
  TH1D* puWeightsUp = (TH1D*)pileupFileMC->Get("puWeightsUp");
  TH1D* puWeightsDown = (TH1D*)pileupFileMC->Get("puWeightsDown");
  puWeights->SetBins(75,0,75);
  puWeightsUp->SetBins(75,0,75);
  puWeightsDown->SetBins(75,0,75);

  TFile* file = TFile::Open( "puppiCorr.root","READ");
  TF1* puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
  TF1* puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
  TF1* puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");


  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  setOutputTree* WWTree = new setOutputTree(outTree);

  int nEvents=0;	
  Long64_t jentry2=0;
  int count_genEvents=0;

  int nInputFiles = sampleName.size();

  if (isLocal==1) nInputFiles = 3;
  cout<<"==> Total number of input files : "<<nInputFiles<<endl;

  TH1D *MCpu = new TH1D("MCpu","",75,0,75);
  TH1D *MCpu_up = new TH1D("MCpu_up","",75,0,75);
  TH1D *MCpu_down = new TH1D("MCpu_down","",75,0,75);
  
  Long64_t TotalNumberOfEvents = 0, nNegEvents = 0; 


  // Set up b-tag scale factor readers
  BTagCalibration calib("csvv2", "CSVv2_Moriond17_B_H.csv");
  BTagCalibrationReader bTagReader(BTagEntry::OP_LOOSE,  // working point: can be OP_LOOSE, OP_MEDIUM, OP_TIGHT 
                                   "central",             // label for the central value (see the scale factor file)
                                   {"up","down"});        // vector of labels for systematics
  bTagReader.load(calib, BTagEntry::FLAV_B, "comb");      // use the "comb" measurements for b-jets
  bTagReader.load(calib, BTagEntry::FLAV_C, "comb");      // use the "comb" measurements for c-jets
  bTagReader.load(calib, BTagEntry::FLAV_UDSG, "incl");   // use the "incl" measurements for light jets
  

  // Loop on input files
  for(int i=0;i<nInputFiles;i++)
  {
     infile = TFile::Open(sampleName[i]);
     eventTree = (TTree*)infile->Get("Events");
     
     TotalNumberOfEvents+=eventTree->GetEntries();
     if(isMC)
     { 
        eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  	TBranch *genBr=0;
     	eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
	for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
	{
	  //eventTree->GetEntry(jentry);
	    genBr->GetEntry(jentry);
    	    infoBr->GetEntry(jentry);	    
	    MCpu->Fill(info->nPUmean);
	    MCpu_up->Fill(info->nPUmeanp);
	    MCpu_down->Fill(info->nPUmeanm);
	    if (jentry2%50000 == 0) std::cout << "\t File no. " << i << "; Neg Event Count; read entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
	    if (gen->weight<0)	nNegEvents++;
	}
     }
     delete infile;
     infile=0, eventTree=0;
  }
  
  
  cout<<"==> Total number of events : "<<TotalNumberOfEvents<<endl;
  cout<<"==> Total number of negative events : "<<nNegEvents<<endl;

  //float weight = std::atof(xSecWeight.c_str())/TotalNumberOfEvents;
  float weight = std::atof(xSecWeight.c_str())/(std::atof(TotalNumberOfEntries.c_str()) - 2*std::atof(TotalNumberOfNegativeEntries.c_str()));
  cout<<"Weight of cross-sec/events = "<<weight<<endl;
  int totalEntries=0;



  //---------start loop on events------------
  std::cout << "---------start loop on events------------" << std::endl;
  jentry2=0;
  int countNans=0;
  for(int i=0;i<nInputFiles;i++)
  {
  cout<<"\n\n=====	Processing File Number : "<<i<<"/"<<nInputFiles<<"\n\t"<<sampleName[i]<<"\n-------"<<endl;

  infile = TFile::Open(sampleName[i]);
  eventTree = (TTree*)infile->Get("Events");
  
  totalEntries+=eventTree->GetEntries();

  nEvents=eventTree->GetEntries();

  cout<<"\t==> Entries = "<<nEvents<<endl;



  eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
  eventTree->SetBranchAddress("AK4CHS",   &jetArr); TBranch *jetBr = eventTree->GetBranch("AK4CHS");    
  eventTree->SetBranchAddress("AK8Puppi",   &vjetArrPuppi); TBranch *vjetBrPuppi = eventTree->GetBranch("AK8Puppi");  
  eventTree->SetBranchAddress("AddAK8Puppi",   &vjetAddArrPuppi); TBranch *vjetAddBrPuppi = eventTree->GetBranch("AddAK8Puppi");  
  TBranch *genBr=0,  *lhePartBr=0;
  if(isMC)
     { 
       eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
       if(eventTree->GetListOfBranches()->FindObject("LHEWeight"))
       {
       eventTree->SetBranchAddress("LHEWeight",&lheWgtArr); lhePartBr = eventTree->GetBranch("LHEWeight");	       }
     }

  for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
  {
    infoBr->GetEntry(jentry);	    

    int GenPassCut = 0;

    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();

    if (jentry2%10000 == 0) std::cout << "\tread entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
    
    //*********************************
    WWTree->initializeVariables(); //initialize all variables

    WWTree->run   = info->runNum;
    WWTree->event = info->evtNum;
    WWTree->lumi  = info->lumiSec;


    /////////////////MC Info
    if (isMC==1)
    {
      lheWgtArr->Clear();
      if(lhePartBr)
	{
	  lhePartBr->GetEntry(jentry);
	}
      genBr->GetEntry(jentry);

    	for (int i = 0; i<lheWgtArr->GetEntries();i++)     // Note that i is starting from 446.
	{
		const baconhep::TLHEWeight *lhe = (baconhep::TLHEWeight*)((*lheWgtArr)[i]);
		WWTree->LHEWeight[i] = lhe->weight;
	}
    }

    
    WWTree->issignal = 0;
    WWTree->wSampleWeight = weight; //xsec/TotalNumberOfEntries
    WWTree->pu_Weight = 1.; //temporary value
    WWTree->pu_Weight_up = 1.; //temporary value
    WWTree->pu_Weight_down = 1.; //temporary value
    WWTree->top1_NNLO_Weight = 1.;
    WWTree->top2_NNLO_Weight = 1.;
    WWTree->id_eff_Weight = 1.;
    WWTree->id_eff_Weight2 = 1.;
    WWTree->trig_eff_Weight = 1.;
    WWTree->trig_eff_Weight2 = 1.;

    if (gen->weight>0)
      WWTree->genWeight=1.;
    else if (gen->weight<0) {
      WWTree->genWeight=-1.;
    }
    cutEff[0]++;

    if (isMC==1)
    {
    if (GenPassCut == 1)   cutEff[1]++;
    }

    
    vertexArr->Clear();
    vertexBr->GetEntry(jentry);
    WWTree->nPV = vertexArr->GetEntries();
  
    //PILE-UP WEIGHT
    if (isMC==1) {
       if(int(info->nPUmean)<75){
           WWTree->pu_Weight = puWeights->GetBinContent(info->nPUmean); //our pu recipe
           WWTree->pu_Weight_up = puWeightsUp->GetBinContent(info->nPUmean); //our pu recipe
           WWTree->pu_Weight_down = puWeightsDown->GetBinContent(info->nPUmean); //our pu recipe
       }
       else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
           std::cout<<"Warning! n_pu too big"<<std::endl;
	   // throw logic_error("n_pu too big");
	   WWTree->pu_Weight = 0.;
	   WWTree->pu_Weight_up = 0.;
	   WWTree->pu_Weight_down = 0.;
       } 
    }


    if(applyTrigger==1)
      if(!(triggerMenu.pass("HLT_IsoMu24_v*",info->triggerBits) || triggerMenu.pass("HLT_IsoTkMu24_v*",info->triggerBits) ||  triggerMenu.pass("HLT_Ele27_WPTight_Gsf_v*",info->triggerBits))) continue;
  
    /////////////////THE SELECTED LEPTON
    int nTightEle=0, nLooseEle=0;
    int nTightMu=0, nLooseMu=0;
    double pt_cut = 20;
    double leadelept_cut = 30;
    double leadmupt_cut = 27;
    electronArr->Clear();
    electronBr->GetEntry(jentry);
    const baconhep::TElectron *leadele = NULL;
    const baconhep::TElectron *subele = NULL;
    //double iso = 1.5;
    for (int i=0; i<electronArr->GetEntries(); i++) {
      const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i]);
      if (ele->pt<=pt_cut) continue;
      if (fabs(ele->eta)>=2.5) continue;
      if(!passEleLooseSel(ele,info->rhoIso)) continue;
      nLooseEle++;
      if(!passEleTightSel(ele,info->rhoIso)) continue;
      ELE.SetPtEtaPhiE(ele->pt,ele->eta,ele->phi,ele->ecalEnergy);
      tightEle.push_back(ELE);
      nTightEle++;
      //iso = ele->chHadIso + TMath::Max( 0.0,(ele->gammaIso + ele->neuHadIso - info->rhoIso*eleEffArea(ele->eta)) );
      if(!leadele || ele->pt>leadele->pt)
	{
	  if(!(ele->pt>leadelept_cut)) continue;
	  subele = leadele;
	  leadele = ele;
	}
      else if (!subele || ele->pt > subele->pt)
	{
	  subele = ele;
	}
    }
    if(leadele)
      {
	WWTree->l_pt1  = leadele->pt;
      }
    if(subele)
      {
	WWTree->l_pt2  = subele->pt;
      }
    muonArr->Clear();
    muonBr->GetEntry(jentry);
    const baconhep::TMuon *leadmu = NULL;
    const baconhep::TMuon *submu = NULL;
    //double leadmue=-999, submue = -999;
    //iso = 1.5;
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);
      if (mu->pt<pt_cut) continue;
      if (fabs(mu->eta)>=2.4) continue;
      if(!passMuonLooseSel(mu)) continue;
      nLooseMu++;
      if(!passMuonTightSel(mu)) continue;
      nTightMu++;
      MU.SetPtEtaPhiM(mu->pt,mu->eta,mu->phi,0.1057);
      tightMuon.push_back(MU);
      //iso = mu->chHadIso + TMath::Max(mu->neuHadIso + mu->gammaIso - 0.5*(mu->puIso), double(0));
      if(!leadmu || mu->pt>leadmu->pt)
	{
	  if(!(mu->pt>leadmupt_cut)) continue;
	  submu = leadmu;
	  leadmu = mu;
	  //leadmue = MU.E();
	}
      else if (!submu || mu->pt > submu->pt)
	{
	  submu = mu;
	  //submue = MU.E();
	}
    }   
    if(leadmu)
      {
	WWTree->l_pt1  = leadmu->pt;
      }
    if(submu)
      {
	WWTree->l_pt2  = submu->pt;
      }
    if(!(WWTree->l_pt1>0)) continue;
    if ((nTightMu+nTightEle)==0) continue; //no leptons with required ID
    if((nLooseEle+nLooseMu)>2) continue;
    if(nTightMu>0 && nLooseEle>0) continue;
    if(nTightEle>0 && nLooseMu>0) continue;
    if(nTightMu==1 && nLooseMu>1) continue;
    if(nTightEle==1 && nLooseEle>1) continue;
    if(nTightMu>0){
      WWTree->type=0;
      leptonName = "mu";	// Added this part for neutrino pz calculation in case there is w-boson.
    }else{
      WWTree->type=1;
      leptonName = "el";
    }
    
    cutEff[2]++;


    ///////////THE FAT JET - AK8
    float tempMassW = 3000.0;
    
    ///////////THE FAT JET - PuppiAK8
    vjetArrPuppi->Clear();
    vjetBrPuppi->GetEntry(jentry);
    vjetAddArrPuppi->Clear();
    vjetAddBrPuppi->GetEntry(jentry);
    tempMassW = 3000.0;
    int nGoodPuppiAK8jets=0;
    
    for ( int i=0; i<vjetArrPuppi->GetEntries(); i++)
      {
	const baconhep::TJet *jet = (baconhep::TJet*)((*vjetArrPuppi)[i]);
	const baconhep::TAddJet *addjet = (baconhep::TAddJet*)((*vjetAddArrPuppi)[i]);
	TLorentzVector TempAK8;
	TempAK8.SetPtEtaPhiM(jet->pt,fabs(jet->eta),jet->phi,jet->mass);
	bool isCleanedJet = true;
	if (jet->pt<200 || fabs(jet->eta)>2.4)  continue; //be careful: this is not inside the synchntuple code
        //if (!passJetLooseSel(jet)) continue;

	//if (abs(jet->mass - 80.385) > tempMassW) continue; //save the jet closest to w-mass
	if (abs(addjet->mass_sd0 - 80.385) > tempMassW) continue; //save the jet closest to w-mass
      
	if (isCleanedJet==false) continue; //jet is overlapped with a lepton
	
	WWTree->ungroomed_PuppiAK8_jet_pt  = jet->pt;
	WWTree->ungroomed_PuppiAK8_jet_eta = jet->eta;
	WWTree->ungroomed_PuppiAK8_jet_phi = jet->phi;
	WWTree->ungroomed_PuppiAK8_jet_e   = TempAK8.E();
      
	WWTree->PuppiAK8_jet_mass  = jet->mass;
	WWTree->PuppiAK8_jet_mass_pr  = addjet->mass_prun;
	WWTree->PuppiAK8_jet_mass_so  = addjet->mass_sd0;
	WWTree->PuppiAK8_jet_mass_so_corr  = addjet->mass_sd0*getPUPPIweight(puppisd_corrGEN, puppisd_corrRECO_cen, puppisd_corrRECO_for, jet->pt, jet->eta);
	WWTree->PuppiAK8_jet_mass_tr  = addjet->mass_trim;
	WWTree->PuppiAK8_jet_tau2tau1 = addjet->tau2/addjet->tau1;
	WWTree->ungroomed_PuppiAK8_jet_charge = jet->q;
	WWTree->PuppiAK8_jet_sj1_pt   = addjet->sj1_pt;
	WWTree->PuppiAK8_jet_sj1_eta  = addjet->sj1_eta;
	WWTree->PuppiAK8_jet_sj1_phi  = addjet->sj1_phi;
	WWTree->PuppiAK8_jet_sj1_m    = addjet->sj1_m;
	WWTree->PuppiAK8_jet_sj1_q    = addjet->sj1_q;
	WWTree->PuppiAK8_jet_sj2_pt   = addjet->sj2_pt;
	WWTree->PuppiAK8_jet_sj2_eta  = addjet->sj2_eta;
	WWTree->PuppiAK8_jet_sj2_phi  = addjet->sj2_phi;
	WWTree->PuppiAK8_jet_sj2_m    = addjet->sj2_m;
	WWTree->PuppiAK8_jet_sj2_q    = addjet->sj2_q;
	WWTree->PuppiAK8_jetID_loose  = passJetLooseSel(jet);


  	WWTree->PuppiAK8jet_e3_b1  = addjet->e3_b1;
  	WWTree->PuppiAK8jet_e3_v1_b1 = addjet->e3_v1_b1;
  	WWTree->PuppiAK8jet_e3_v2_b1 = addjet->e3_v2_b1;
  	WWTree->PuppiAK8jet_e4_v1_b1 = addjet->e4_v1_b1;
  	WWTree->PuppiAK8jet_e4_v2_b1 = addjet->e4_v2_b1;
  	WWTree->PuppiAK8jet_e3_b2  = addjet->e3_b2;
  	WWTree->PuppiAK8jet_e3_v1_b2 = addjet->e3_v1_b2;
  	WWTree->PuppiAK8jet_e3_v2_b2 = addjet->e3_v2_b2;
  	WWTree->PuppiAK8jet_e4_v1_b2 = addjet->e4_v1_b2;
  	WWTree->PuppiAK8jet_e4_v2_b2 = addjet->e4_v2_b2;
  	
  	WWTree->PuppiAK8jet_e2_sdb1  = addjet->e2_sdb1 ;
  	WWTree->PuppiAK8jet_e3_sdb1  = addjet->e3_sdb1 ;
  	WWTree->PuppiAK8jet_e3_v1_sdb1 = addjet->e3_v1_sdb1  ;
  	WWTree->PuppiAK8jet_e3_v2_sdb1 = addjet->e3_v2_sdb1  ;
  	WWTree->PuppiAK8jet_e4_v1_sdb1 = addjet->e4_v1_sdb1  ;
  	WWTree->PuppiAK8jet_e4_v2_sdb1 = addjet->e4_v2_sdb1;    
  	
  	WWTree->PuppiAK8jet_e2_sdb2  = addjet->e2_sdb2 ;
  	WWTree->PuppiAK8jet_e3_sdb2  = addjet->e3_sdb2 ;
  	WWTree->PuppiAK8jet_e3_v1_sdb2 = addjet->e3_v1_sdb2  ;
  	WWTree->PuppiAK8jet_e3_v2_sdb2 = addjet->e3_v2_sdb2  ;
  	WWTree->PuppiAK8jet_e4_v1_sdb2 = addjet->e4_v1_sdb2  ;
  	WWTree->PuppiAK8jet_e4_v2_sdb2 = addjet->e4_v2_sdb2;    // Soft Dropped correlation function in puts beta=2
  	
  	
  	WWTree->PuppiAK8jet_qjet = addjet->qjet;
	
	tempMassW = abs(WWTree->PuppiAK8_jet_mass_so - 80.385);
	nGoodPuppiAK8jets++;
      }
    if (WWTree->ungroomed_PuppiAK8_jet_pt > 0.)
      {
	JET_PuppiAK8.SetPtEtaPhiM(WWTree->ungroomed_PuppiAK8_jet_pt,WWTree->ungroomed_PuppiAK8_jet_eta,WWTree->ungroomed_PuppiAK8_jet_phi,WWTree->PuppiAK8_jet_mass_so_corr);
	SJ1_PuppiAK8.SetPtEtaPhiM(WWTree->PuppiAK8_jet_sj1_pt, WWTree->PuppiAK8_jet_sj1_eta, WWTree->PuppiAK8_jet_sj1_phi, WWTree->PuppiAK8_jet_sj1_m);
	SJ2_PuppiAK8.SetPtEtaPhiM(WWTree->PuppiAK8_jet_sj2_pt, WWTree->PuppiAK8_jet_sj2_eta, WWTree->PuppiAK8_jet_sj2_phi, WWTree->PuppiAK8_jet_sj2_m);
      }
    
    
    // FAT JET SELECTION
    bool isGoodFatJet = true;
    if (nGoodPuppiAK8jets==0) isGoodFatJet = false; //not found a good hadronic W candidate
    if (WWTree->ungroomed_PuppiAK8_jet_pt<200) isGoodFatJet = false;
    if (!isGoodFatJet) continue;
    cutEff[5]++;
    
    
    cutEff[8]++;
    

    
    /////////VBF and b-tag section
    WWTree->njets=0;
    WWTree->nBTagJet_loose=0;
    WWTree->nBTagJet_medium=0;
    WWTree->nBTagJet_tight=0;
    
    WWTree->njets_unmerged=0;
    WWTree->nBTagJet_loose_unmerged=0;
    WWTree->nBTagJet_medium_unmerged=0;
    WWTree->nBTagJet_tight_unmerged=0;

    int OnlyTwoVBFTypeJets = 0;
    
    std::vector<int> indexGoodVBFJets;

    jetArr->Clear();
    jetBr->GetEntry(jentry);
    for ( int i=0; i<jetArr->GetEntries(); i++) //loop on AK4 jet
    {
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[i]);
      bool isCleanedFromFatJet = true;

      
      if (jet->pt<=20 ) continue;
      if (!passJetLooseSel(jet)) continue;

      //fill B-Tag info
      
      if (fabs(jet->eta) < 2.4 && jet->pt>20){
      		if (jet->csv>0.5426)  WWTree->nBTagJet_loose++;
      		if (jet->csv>0.8484)  WWTree->nBTagJet_medium++;
      		if (jet->csv>0.9535)  WWTree->nBTagJet_tight++;
      }
      
      //CLEANING FROM FAT JET
      if (WWTree->nGoodAK8jets > 0) {
        if (deltaR(WWTree->ungroomed_AK8jet_eta, WWTree->ungroomed_AK8jet_phi,
                   jet->eta,jet->phi) < 0.8 )
          isCleanedFromFatJet = false;
      } 
      
      
      if (isCleanedFromFatJet==false) continue;
      
      indexGoodVBFJets.push_back(i); //save index of the "good" vbf jets candidates
      
      WWTree->njets++;
      AK4.SetPtEtaPhiM(jet->pt,jet->eta,jet->phi,jet->mass);
      
      //------------------------------
      // !!! VBF non-Puppi missing !!!
      //------------------------------
    }

    int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
    if (indexGoodVBFJets.size()>=4) 
    {
      cutEff[9]++;
      float tempPtMax=0.;
      float DRvbf;
      nVBF1=-1; nVBF2=-1; //position of the two vbf jets
      
      for (std::size_t i=0; i<indexGoodVBFJets.size()-1; i++) {
        for ( std::size_t ii=i+1; ii<indexGoodVBFJets.size(); ii++) {
	  const baconhep::TJet *jet1 = (baconhep::TJet*)((*jetArr)[indexGoodVBFJets.at(i)]);
	  const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArr)[indexGoodVBFJets.at(ii)]);
	  if (jet1->pt < 30) continue;
	  if (jet2->pt < 30) continue;
          VBF1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
          VBF2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
          TOT = VBF1 + VBF2;

	  if (TOT.M()<500 ) continue;
	  if ( VBFSel==1)
	  {
		if (TOT.Pt() < tempPtMax) continue;
		tempPtMax = TOT.Pt(); //take the jet pair with largest Pt
	  }
	  else if ( VBFSel==2)
	  {
		if (TOT.M() < tempPtMax) continue;
		tempPtMax = TOT.M(); //take the jet pair with largest mjj
	  }
	  else if ( VBFSel==3)
	  {
		DRvbf = abs(VBF1.Eta()-VBF2.Eta());
		if (DRvbf < tempPtMax) continue;
		tempPtMax = DRvbf; //take the jet pair with largest dEtajj
	  }
	  else
	  {
	  	cout<<"\n\nERROR:	Enter valid vbf selection criteria....\n\n"<<endl;
		exit(0);
	  }
          nVBF1 = indexGoodVBFJets.at(i); //save position of the 1st vbf jet
          nVBF2 = indexGoodVBFJets.at(ii); //save position of the 2nd vbf jet
        }
      }
      if (nVBF1 == -1 ) continue;
      if (nVBF2 == -1 ) continue;
      if (nVBF1!=-1 && nVBF2 !=-1) OnlyTwoVBFTypeJets=1;
      
      if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
      {
        // nVBF1=0; nVBF2=1;
	//cout<<"nVBF1 = "<<nVBF1<<"\tnVBF2 = "<<nVBF2<<endl;
        const baconhep::TJet *jet1 = (baconhep::TJet*)((*jetArr)[nVBF1]);
	const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArr)[nVBF2]);
        VBF1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
        VBF2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
        TOT = VBF1 + VBF2;
	

		//cout<<"****** > "<<jet1->pt<<"\t"<<jet2->pt<<endl;
        	WWTree->vbf_maxpt_j1_pt = jet1->pt;
        	WWTree->vbf_maxpt_j1_eta = jet1->eta;
        	WWTree->vbf_maxpt_j1_phi = jet1->phi;
        	WWTree->vbf_maxpt_j1_e = VBF1.E();
        	WWTree->vbf_maxpt_j1_mass = VBF1.M();
        	WWTree->vbf_maxpt_j1_bDiscriminatorCSV = jet1->csv;
		WWTree->vbf_maxpt_j1_charge = jet1->q;
        	WWTree->vbf_maxpt_j2_pt = jet2->pt;
        	WWTree->vbf_maxpt_j2_eta = jet2->eta;
        	WWTree->vbf_maxpt_j2_phi = jet2->phi;
        	WWTree->vbf_maxpt_j2_e = VBF2.E();
        	WWTree->vbf_maxpt_j2_mass = VBF2.M();
        	WWTree->vbf_maxpt_j2_bDiscriminatorCSV = jet2->csv;
		WWTree->vbf_maxpt_j2_charge = jet2->q;


        WWTree->vbf_maxpt_jj_pt = TOT.Pt();
        WWTree->vbf_maxpt_jj_eta = TOT.Eta();
        WWTree->vbf_maxpt_jj_phi = TOT.Phi();
        WWTree->vbf_maxpt_jj_m = TOT.M();	
	WWTree->vbf_maxpt_jj_Deta = abs(VBF1.Eta() - VBF2.Eta());
	WWTree->AK4_DR_GENRECO_11 = abs(deltaR(WWTree->AK4_1_eta_gen, WWTree->AK4_1_phi_gen, WWTree->vbf_maxpt_j1_eta, WWTree->vbf_maxpt_j1_phi));
	WWTree->AK4_DR_GENRECO_12 = abs(deltaR(WWTree->AK4_1_eta_gen, WWTree->AK4_1_phi_gen, WWTree->vbf_maxpt_j2_eta, WWTree->vbf_maxpt_j2_phi));
	WWTree->AK4_DR_GENRECO_21 = abs(deltaR(WWTree->AK4_2_eta_gen, WWTree->AK4_2_phi_gen, WWTree->vbf_maxpt_j1_eta, WWTree->vbf_maxpt_j1_phi));
	WWTree->AK4_DR_GENRECO_22 = abs(deltaR(WWTree->AK4_2_eta_gen, WWTree->AK4_2_phi_gen, WWTree->vbf_maxpt_j2_eta, WWTree->vbf_maxpt_j2_phi));
      }
    //indexGoodVBFJets.clear();
    jetArr->Clear();
    jetBr->GetEntry(jentry);
    double ptBalanceForLepSearch = 0.0;
    //cout<<"DEBUG : 1 : nVBF1 = "<<nVBF1<<"\tnVBF2 = "<<nVBF2<<endl;
    double min = 99999999.0;
    for (std::size_t i=0; i<indexGoodVBFJets.size()-1; i++) {
      for ( std::size_t ii=i+1; ii<indexGoodVBFJets.size(); ii++) {
	if (indexGoodVBFJets.at(i) == nVBF1) continue;
	if (indexGoodVBFJets.at(ii) == nVBF1) continue;
	if (indexGoodVBFJets.at(i) == nVBF2) continue;
	if (indexGoodVBFJets.at(ii) == nVBF2) continue;
        const baconhep::TJet *jet1 = (baconhep::TJet*)((*jetArr)[indexGoodVBFJets.at(i)]);
        const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArr)[indexGoodVBFJets.at(ii)]);
        if (jet1->pt < 30) continue;
        if (jet2->pt < 30) continue;
        if (fabs(jet1->eta) > 2.5) continue;
        if (fabs(jet2->eta) > 2.5) continue;
        VBF1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
        VBF2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
        TOT = VBF1 + VBF2;
	ptBalanceForLepSearch = ( WWTree->ungroomed_PuppiAK8_jet_pt + WWTree->vbf_maxpt_j1_pt + WWTree->vbf_maxpt_j2_pt) - TOT.Pt(); 
	if (min > fabs(ptBalanceForLepSearch))
	{
		WWTree->l_pt1 = jet1->pt;
		WWTree->l_eta1 = jet1->eta; 
		WWTree->l_phi1 = jet1->phi;
		WWTree->l_e1 = jet1->mass;
		WWTree->l_charge1 = jet1->q;
		LEP1 = VBF1;

		WWTree->l_pt2 = jet2->pt;
		WWTree->l_eta2 = jet2->eta; 
		WWTree->l_phi2 = jet2->phi;
		WWTree->l_e2 = jet2->mass;
		WWTree->l_charge2 = jet2->q;
		LEP2 = VBF2;

		WWTree->dilep_pt  = (LEP1+LEP2).Pt();
		WWTree->dilep_eta = (LEP1+LEP2).Eta();
		WWTree->dilep_phi = (LEP1+LEP2).Phi();	
		WWTree->dilep_m = (LEP1+LEP2).M();	

		min = fabs(ptBalanceForLepSearch);
		//cout<<jentry<<"\tMin = "<<min<<endl;
	}
    	WWTree->dPtForLepMETBalance = min;

      }
    }
    }		// if (indexGoodVBFJets.size()>=4)

    
    //////////////////FOUR-BODY INVARIANT MASS
    WWTree->mass_lvj_type0 = (LEP1 + LEP2 + JET_PuppiAK8).M();
    WWTree->WWEta = (LEP1 + LEP2 + JET_PuppiAK8 ).Eta();
    WWTree->WWRapidity = (LEP1 + LEP2 + JET_PuppiAK8 ).Rapidity();
    WWTree->pt_lvj_type0_PuppiAK8 = (LEP1 + LEP2 + JET_PuppiAK8).Pt();
    WWTree->eta_lvj_type0_PuppiAK8 = (LEP1 + LEP2 + JET_PuppiAK8).Eta();
    WWTree->rapidity_lvj_type0_PuppiAK8 = (LEP1 + LEP2 + JET_PuppiAK8).Rapidity();
    WWTree->phi_lvj_type0_PuppiAK8 = (LEP1 + LEP2 + JET_PuppiAK8).Phi();
    WWTree->energy_lvj_type0_PuppiAK8 = (LEP1 + LEP2 + JET_PuppiAK8).E();
    WWTree->mt_lvj_type0_PuppiAK8 = (LEP1 + LEP2 + JET_PuppiAK8).Mt();
    WWTree->LepWEta = (LEP1 + LEP2 ).Eta();
    WWTree->LepWRapidity = (LEP1 + LEP2 ).Rapidity();
    WWTree->HadWEta = (JET_PuppiAK8 ).Eta();
    WWTree->HadWRapidity = (JET_PuppiAK8 ).Rapidity();
    

///////////////////////////////////////////////////////////////////////////////////////
    if (OnlyTwoVBFTypeJets == 1) WWTree->isVBF=1;
    if (OnlyTwoVBFTypeJets == 0 ) continue;
        cutEff[10]++;
    
    WWTree->totalEventWeight = WWTree->genWeight*WWTree->pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight*WWTree->id_eff_Weight;
    WWTree->totalEventWeight_2Lep = WWTree->genWeight*WWTree->pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight*WWTree->id_eff_Weight*WWTree->trig_eff_Weight2*WWTree->id_eff_Weight2;

    
    WWTree->nEvents = TotalNumberOfEvents;
    WWTree->nNegEvents = nNegEvents;
    WWTree->nTotEvents = std::atof(TotalNumberOfEntries.c_str());
    WWTree->nTotNegEvents = std::atof(TotalNumberOfNegativeEntries.c_str());


///////////////////////////////////////////////////
//
//	CHS ANGULAR VARIABLES
//
//////////////////////////////////////////////////
    if (WWTree->isVBF && nGoodPuppiAK8jets!=0){
    WWTree->PtBalance_type0 = ((JET_PuppiAK8+LEP1 + LEP2).Pt())/(JET_PuppiAK8.Pt()+(LEP1 + LEP2).Pt());

    WWTree->BosonCentrality_type0 = GetMin( GetMin(JET_PuppiAK8.Eta(),(LEP1+LEP2).Eta())-GetMin(WWTree->vbf_maxpt_j1_eta, WWTree->vbf_maxpt_j2_eta) , GetMax(WWTree->vbf_maxpt_j1_eta,WWTree->vbf_maxpt_j2_eta) - GetMax(JET_PuppiAK8.Eta(),(LEP1+LEP2).Eta())  );



    if (JET_PuppiAK8.Pt()>0){
    double a_costheta1, a_costheta2, a_costhetastar, a_Phi, a_Phi1;

    computeAngles( LEP1 + LEP2 + JET_PuppiAK8, LEP1 + LEP2, LEP1, LEP2, JET_PuppiAK8,  SJ1, SJ2, a_costheta1, a_costheta2, a_Phi, a_costhetastar, a_Phi1 );
    WWTree->costheta1_type0 = (float) a_costheta1;                
    if (LEP2.Beta() > 1.0) printf("------------- LEP2 beta = %.17g\n", LEP2.Beta());
    if (LEP2.Beta() > 1.0) std::cout << "#######    LEP2 beta = 1 + " << (LEP2.Beta() - 1.0) << std::endl;
    

    WWTree->costheta2_type0 = (float) a_costheta2;

    WWTree->costhetastar_type0 = (float) a_costhetastar;

    WWTree->phi_type0 = (float) a_Phi;

    WWTree->phi1_type0 = (float) a_Phi1;
    if ((isnan(WWTree->phi1_type0) == 1) || (isnan(WWTree->phi_type0) == 1) || (isnan(WWTree->costheta2_type0) == 1) || (isnan(WWTree->costheta1_type0) == 1) || (isnan(WWTree->costhetastar_type0) == 1))
    	{
	countNans++;
    	//cout<<jentry2<< "\tcostheta1_type0 is NaN" << endl;
	//if (LEP2.Beta()>1) cout<<"beta > 1"<<endl;
	//printf("------------- LEP2 beta = %.17g\n", LEP2.Beta());
	//std::cout << "#######    LEP2 beta = 1 + " << (LEP2.Beta() - 1.0) << std::endl;
	//printf(" Neutrino mass = %.17g\n", LEP2.M());
	}


    if (fabs(VBF1.Eta() - VBF2.Eta()) == 0.0)
    {  	 WWTree->VBSCentrality_type0 = -999.0;	}
    else
    {
    	WWTree->VBSCentrality_type0 = (fabs(VBF1.Eta()- (((LEP1 + LEP2).Eta()+JET_PuppiAK8.Eta()))- VBF2.Eta() ))/fabs(VBF1.Eta() - VBF2.Eta());
    }
     WWTree->RpT_type0 = (JET_PuppiAK8.Pt()*(LEP1 + LEP2).Pt())/(VBF1.Pt()*VBF2.Pt());
     WWTree->ZeppenfeldWH = JET_PuppiAK8.Eta() - (VBF1.Eta() + VBF2.Eta())/2.0;
     }
     WWTree->ZeppenfeldWL_type0 = (LEP1 + LEP2).Eta() - (VBF1.Eta() + VBF2.Eta())/2.0;
     WWTree->LeptonProjection_type0 = (LEP1.Pt()*cos(LEP1.Theta()-(LEP1 + LEP2).Theta()))/(LEP1 + LEP2).Pt();
    }

/////////////////////////////////////////////////// END	CHS ANGULAR VARIABLES


    outTree->Fill();

    }
    delete infile;
    infile=0, eventTree=0;
    /////////////////FILL THE TREE
  }
  //delete puWeight;	delete puWeight_up;	delete puWeight_down;
  delete MCpu;	delete MCpu_up;	delete MCpu_down;
  delete puWeightsDown;	delete puWeightsUp;	delete puWeights;
  //delete pileupHisto;
  //pileupFile->Close();
  pileupFileMC->Close();
  file->Close();
  std::cout << "---------end loop on events------------" << std::endl;
  std::cout << std::endl;
  std::cout << "GEN events = " << count_genEvents << std::endl;


  
  std::cout << "----------------------" << std::endl;
  std::cout << " SUMMARY" << std::endl;
  std::cout << "----------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"MC matching: "<<(float)ok/(float)total<<std::endl;
  std::cout<<"negative events: "<<nNegEvents<<std::endl;
  std::cout<<"Total number of events having NaNs = "<<countNans<<std::endl;
  std::cout << std::endl;
  std::cout<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<std::endl
	   <<"(2) tight lepton:      "<<cutEff[2]<<"\t:\t"<<((float)cutEff[2]*100.0)/(float)cutEff[0]<<std::endl
	   <<"(3) MET:               "<<cutEff[3]<<"\t:\t"<<((float)cutEff[3]*100.0)/(float)cutEff[2]<<std::endl
	   <<"(4) negative lep-MET:  "<<cutEff[4]<<"\t:\t"<<((float)cutEff[4]*100.0)/(float)cutEff[3]<<std::endl
	   <<"(5) 1 good AK8:        "<<cutEff[5]<<"\t:\t"<<((float)cutEff[5]*100.0)/(float)cutEff[4]<<std::endl
//	   <<"(6) 2 good AK4:        "<<cutEff[6]<<"\t:\t"<<((float)cutEff[6]*100.0)/(float)cutEff[5]<<std::endl
//	   <<"(7) 1 AK8 & 2 good AK4:"<<cutEff[7]<<"\t:\t"<<((float)cutEff[7]*100.0)/(float)cutEff[6]<<std::endl
	   <<"(8) m(WV) > 500:       "<<cutEff[8]<<"\t:\t"<<((float)cutEff[8]*100.0)/(float)cutEff[5]<<std::endl
	   <<"(9) >=2 good VBF jets: "<<cutEff[9]<<"\t:\t"<<((float)cutEff[9]*100.0)/(float)cutEff[8]<<std::endl
	   <<"(10) Found VBF jets:  "<<cutEff[10]<<"\t:\t"<<((float)cutEff[10]*100.)/(float)cutEff[9]<<std::endl;
  
 
  //--------close everything-------------
  delete info; delete gen;
  delete vertexArr; delete muonArr; delete electronArr; 
  delete jetArr; delete vjetArrPuppi; delete vjetAddArrPuppi;  delete lheWgtArr;
  outROOT->Write();
  outROOT->Close();
  int t1 = time(NULL);
  printf ("\n==> time to run this code = %0.3f min\n", (float)(t1 - t0)/60.0);
  return(0);
}

//////////////////////////////////
//Ref: https://github.com/ram1123/LHEAnalyzer/blob/LHEanalyzer/LHEanalyzer.cpp
//////////////////////////////////
void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1){
    
    ///////////////////////////////////////////////
    // check for z1/z2 convention, redefine all 4 vectors with convention
    ///////////////////////////////////////////////	
    TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
    p4H = thep4H;
    
    p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
    p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
    //// costhetastar
	TVector3 boostX = -(thep4H.BoostVector());
	TLorentzVector thep4Z1inXFrame( p4Z1 );
	TLorentzVector thep4Z2inXFrame( p4Z2 );	
	thep4Z1inXFrame.Boost( boostX );
	thep4Z2inXFrame.Boost( boostX );
	TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
	TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );    
    costhetastar = theZ1X_p3.CosTheta();
    
    //// --------------------------- costheta1
    TVector3 boostV1 = -(thep4Z1.BoostVector());
    TLorentzVector p4M11_BV1( p4M11 );
	TLorentzVector p4M12_BV1( p4M12 );	
    TLorentzVector p4M21_BV1( p4M21 );
	TLorentzVector p4M22_BV1( p4M22 );
    p4M11_BV1.Boost( boostV1 );
	p4M12_BV1.Boost( boostV1 );
	p4M21_BV1.Boost( boostV1 );
	p4M22_BV1.Boost( boostV1 );
    
    TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
    //// costheta1
    costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();
    
    //// --------------------------- costheta2
    TVector3 boostV2 = -(thep4Z2.BoostVector());
    TLorentzVector p4M11_BV2( p4M11 );
	TLorentzVector p4M12_BV2( p4M12 );	
    TLorentzVector p4M21_BV2( p4M21 );
	TLorentzVector p4M22_BV2( p4M22 );
    p4M11_BV2.Boost( boostV2 );
	p4M12_BV2.Boost( boostV2 );
	p4M21_BV2.Boost( boostV2 );
	p4M22_BV2.Boost( boostV2 );
    
    TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
    //// costheta2
    costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
    
    //// --------------------------- Phi and Phi1
    //    TVector3 boostX = -(thep4H.BoostVector());
    TLorentzVector p4M11_BX( p4M11 );
	TLorentzVector p4M12_BX( p4M12 );	
    TLorentzVector p4M21_BX( p4M21 );
	TLorentzVector p4M22_BX( p4M22 );	
    
	p4M11_BX.Boost( boostX );
	p4M12_BX.Boost( boostX );
	p4M21_BX.Boost( boostX );
	p4M22_BX.Boost( boostX );
    
    TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
    TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );    
    
    TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
    TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 
    
    //// Phi
    TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;    
    double tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
    double sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
    Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
    
    
    //////////////
    
    TVector3 beamAxis(0,0,1);
    TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();
    
    TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
    TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
    TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
    
    //// Phi1
    double tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
    double sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
    Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );    
    
    //    std::cout << "extractAngles: " << std::endl;
    //    std::cout << "costhetastar = " << costhetastar << ", costheta1 = " << costheta1 << ", costheta2 = " << costheta2 << std::endl;
    //    std::cout << "Phi = " << Phi << ", Phi1 = " << Phi1 << std::endl;    
    
}

  double GetSFs_Lepton(double pt, double eta, TH1F *h1){
	if (pt > h1->GetYaxis()->GetXmax())  // Check if pt is not ouside upper limit; if so then to assign same SF as the last bin to high pt resacle tempPt to less then maximum range of pt defined in histo.	
		pt = h1->GetYaxis()->GetXmax() - 1.0;	
	if (pt < h1->GetYaxis()->GetXmin()) // Check if pt is not ouside lower limit; if so then to assign same SF as the first bin to low pt resacle tempPt to greater then minimum range of pt defined in histo.
		pt = h1->GetYaxis()->GetXmin() + 1.0;

	return h1->GetBinContent(h1->GetXaxis()->FindFixBin(eta), h1->GetYaxis()->FindFixBin(pt));
  }	

  double GetMin(double x, double y){
  	if (x<y) return x;
	else	 return y;
  }
  double GetMax(double x, double y){
  	if (x>y) return x;
	else	 return y;
  }

float getPUPPIweight(TF1* puppisd_corrGEN, TF1* puppisd_corrRECO_cen, TF1* puppisd_corrRECO_for, float puppipt, float puppieta ){



  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
        
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ){
    recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  }
  else{
    recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  }
  
  totalWeight = genCorr * recoCorr;


  return totalWeight;
}
double func( double pt, double eta, JetCorrectionUncertainty *fJetUnc) { 
  fJetUnc->setJetPt ( pt  );
  fJetUnc->setJetEta( eta );
  return fJetUnc->getUncertainty(true);
}
double GetPt_MET(double pfMET, double phi, double pz){
	
	double px = pfMET*TMath::Cos(phi);
	double py = pfMET*TMath::Sin(phi);
	return TMath::Sqrt(px*px + py*py); 
}
double GetEta_MET(double pfMET, double phi, double pz){
	double px = pfMET*TMath::Cos(phi);
	double py = pfMET*TMath::Sin(phi);
	double p = TMath::Sqrt(px*px + py*py + pz*pz);
	return 0.5*log((p+pz)/(p-pz));
}
