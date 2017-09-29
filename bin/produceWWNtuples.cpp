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
  //std::string TotalNumberOfEntries = argv[8];
  //float LUMI = atof(argv[9]);
  int applyTrigger = atoi(argv[10]);
  std::string jsonFileName = argv[11];
  int isLocal = atoi(argv[12]);
  int VBFSel  = atoi(argv[13]);
  
  std::string leptonName;

  if ( VBFSel==1)
    {
      cout<<"==> VBF selection method : Select two highest pT jets"<<endl;
    }
  else if ( VBFSel==2)
    {
      cout<<"==> VBF selection method : Select pair with highest mjj..."<<endl;
    }
  else if ( VBFSel==3)
    {
      cout<<"==> VBF selection method : Select pair with highest DeltaEta..."<<endl;
    }
  else
    {
      cout<<"\n\nERROR:	Enter valid vbf selection criteria....\n\n"<<endl;
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
  std::cout<<"apply trigger: "<<applyTrigger<<std::endl;

  TLorentzVector W,W_puppi,LEP, LEP2, SJ1, SJ2;
  TLorentzVector NU0,NU1,NU2,NU0_puppi,NU1_puppi,NU2_puppi;
  TLorentzVector NU0_jes_up, NU0_jes_dn;
  TLorentzVector JET, JET_PuppiAK8, AK4;
  TLorentzVector JET_jes_up, JET_jes_dn, JET_PuppiAK8_jes_up, JET_PuppiAK8_jes_dn;
  TLorentzVector AK4_JET1,AK4_JET2;
  TLorentzVector AK4_JET1_jes_up, AK4_JET1_jes_dn;
  TLorentzVector AK4_JET2_jes_up, AK4_JET2_jes_dn;
  TLorentzVector PuppiAK4_JET1,PuppiAK4_JET2;
  TLorentzVector PuppiAK4_JET1_jes_up, PuppiAK4_JET1_jes_dn;
  TLorentzVector PuppiAK4_JET2_jes_up, PuppiAK4_JET2_jes_dn;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector ELE,MU;

  
  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;

  int ok=0, total=0;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen  = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr    = new TClonesArray("baconhep::TMuon");
  TClonesArray *electronArr    = new TClonesArray("baconhep::TElectron");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  TClonesArray *vjetArr  = new TClonesArray("baconhep::TJet");
  TClonesArray *vjetAddArr  = new TClonesArray("baconhep::TAddJet");
  TClonesArray *jetArr  = new TClonesArray("baconhep::TJet");
  TClonesArray *vjetArrPuppi  = new TClonesArray("baconhep::TJet");
  TClonesArray *vjetAddArrPuppi  = new TClonesArray("baconhep::TAddJet");
  TClonesArray *jetArrPuppi  = new TClonesArray("baconhep::TJet");
  TClonesArray *lheWgtArr = new TClonesArray("baconhep::TLHEWeight");
  

  char command1[3000];
  //exit(0);
  if ( cluster == "lxplus")
  	sprintf(command1, "eos find -f %s  | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://eoscms.cern.ch/\"$1}' > listTemp_%s.txt", (inputFolder).c_str(), outputFile.c_str());	// NF in awk command skips the blank line
  else 
  	sprintf(command1,"xrdfs root://cmseos.fnal.gov ls %s | awk '{print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());

  std::cout<<command1<<std::endl;
  system(command1);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);
  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);

  int fileCounter=0;
  Long64_t TotalNumberOfEvents=0;

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

  //--------pile up file -----------------
  //TFile* pileupFile = TFile::Open("PileUpData2016_23Sep2016ReReco_69200ub.root");  
  //TH1D* pileupHisto = (TH1D*)pileupFile->Get("pileup");

  TFile* pileupFileMC = TFile::Open("puWeights_80x_37ifb.root");
  TH1D* puWeights = (TH1D*)pileupFileMC->Get("puWeights");
  TH1D* puWeightsUp = (TH1D*)pileupFileMC->Get("puWeightsUp");
  TH1D* puWeightsDown = (TH1D*)pileupFileMC->Get("puWeightsDown");
  puWeights->SetBins(75,0,75);
  puWeightsUp->SetBins(75,0,75);
  puWeightsDown->SetBins(75,0,75);

  //---------------- Root Files for ID, ISO, Trigger, GSF correctiosn for ELE and MU both: Starts -------------
  TFile* IDIsoEle = TFile::Open("egammaEffi_EGM2D_TightCutBasedIDSF.root","READ");
  TH1F *hIDIsoEle = (TH1F*)IDIsoEle->Get("EGamma_SF2D");

  TFile* GSFCorrEle = TFile::Open("egammaEffi_SF2D_GSF_tracking.root","READ");
  TH1F *hGSFCorrEle = (TH1F*)GSFCorrEle->Get("EGamma_SF2D");

  TFile* TriggerEle = TFile::Open("ElectronTrigger_SF.root","READ");
  TH1F* hTriggerEle = (TH1F*)TriggerEle->Get("HLT_Ele27");

  TFile* IDMuA = TFile::Open("MuonID_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F *hIDMuA = (TH1F*)IDMuA->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile* IDMuB = TFile::Open("MuonID_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F *hIDMuB = (TH1F*)IDMuB->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile* IsoMuA = TFile::Open("MuonIso_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F *hIsoMuA = (TH1F*)IsoMuA->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

  TFile* IsoMuB = TFile::Open("MuonIso_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F *hIsoMuB = (TH1F*)IsoMuB->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

  TFile* TriggerMuA = TFile::Open("MuonTrigger_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F* hTriggerMuA = (TH1F*)TriggerMuA->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

  TFile* TriggerMuB = TFile::Open("MuonTrigger_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F* hTriggerMuB = (TH1F*)TriggerMuB->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

  //---------------- Root Files for ID, ISO, Trigger, GSF correctiosn for ELE and MU both: ENDS -------------


  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  setOutputTree* WWTree = new setOutputTree(outTree);

  int nEvents=0;	
  Long64_t jentry2=0;
  int count_genEvents=0;

  int nInputFiles = sampleName.size();

  if (isLocal==1) nInputFiles=1;
  cout<<"==> Total number of input files : "<<nInputFiles<<endl;

  TH1D *MCpu = new TH1D("MCpu","",75,0,75);
  TH1D *MCpu_up = new TH1D("MCpu_up","",75,0,75);
  TH1D *MCpu_down = new TH1D("MCpu_down","",75,0,75);
  
  int nNegEvents=0; 
  int RunEntry = 500000;
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
	//for (Long64_t jentry=0; jentry<RunEntry;jentry++,jentry2++)
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
  //puWeight->Divide(MCpu);

  float weight = std::atof(xSecWeight.c_str())/TotalNumberOfEvents;
  int totalEntries=0;

  //---------start loop on events------------
  std::cout << "---------start loop on events------------" << std::endl;
  jentry2=0;
  for(int i=0;i<nInputFiles;i++)
  {
  cout<<"\n\n=====	Processing File Number : "<<i<<"\n\t"<<sampleName[i]<<"\n-------"<<endl;

  infile = TFile::Open(sampleName[i]);
  eventTree = (TTree*)infile->Get("Events");
  
  totalEntries+=eventTree->GetEntries();

  nEvents=eventTree->GetEntries();

  cout<<"\t==> Entries = "<<nEvents<<endl;



  eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
  eventTree->SetBranchAddress("AK8CHS",   &vjetArr); TBranch *vjetBr = eventTree->GetBranch("AK8CHS");  
  eventTree->SetBranchAddress("AddAK8CHS",   &vjetAddArr); TBranch *vjetAddBr = eventTree->GetBranch("AddAK8CHS");  
  eventTree->SetBranchAddress("AK4CHS",   &jetArr); TBranch *jetBr = eventTree->GetBranch("AK4CHS");    
  eventTree->SetBranchAddress("AK8Puppi",   &vjetArrPuppi); TBranch *vjetBrPuppi = eventTree->GetBranch("AK8Puppi");  
  eventTree->SetBranchAddress("AddAK8Puppi",   &vjetAddArrPuppi); TBranch *vjetAddBrPuppi = eventTree->GetBranch("AddAK8Puppi");  
  eventTree->SetBranchAddress("AK4Puppi",   &jetArrPuppi); TBranch *jetBrPuppi = eventTree->GetBranch("AK4Puppi");  
  TBranch *genBr=0, *genPartBr=0, *lhePartBr=0;
  if(isMC)
     { 
       eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
       eventTree->SetBranchAddress("GenParticle",&genPartArr); genPartBr = eventTree->GetBranch("GenParticle");
       if(eventTree->GetListOfBranches()->FindObject("LHEWeight"))
       {
       eventTree->SetBranchAddress("LHEWeight",&lheWgtArr); lhePartBr = eventTree->GetBranch("LHEWeight");	       }
     }

  //for (Long64_t jentry=0; jentry<RunEntry;jentry++,jentry2++)
  for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
  {
    //eventTree->GetEntry(jentry);
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
      genPartArr->Clear();
      genBr->GetEntry(jentry);
      genPartBr->GetEntry(jentry);
      TLorentzVector hadW, lepW, VBFJ1, VBFJ2, VBFJ, temp;
      TLorentzVector genLep, genNeutrino, genWquarks, genVBFquarks;
      std::vector<TLorentzVector> v_GEN_hadW, v_GEN_lepW, v_GEN_VBFJ1, v_GEN_VBFJ2, v_GEN_VBFJ, v_GEN_temp;
      std::vector<TLorentzVector> v_genLep, v_genNeutrino, v_genWquarks, v_genVBFquarks;
      
      for (int i = 0; i<genPartArr->GetEntries();i++)
	{
	  const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
	  Int_t parentPdg=dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId;
	  if( (abs(genloop->pdgId) == 11 || abs(genloop->pdgId) == 13 ) && abs(parentPdg) == 24)
	    {
	      genLep.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genLep.push_back(genLep);
	    }
	  if( (abs(genloop->pdgId) == 12 || abs(genloop->pdgId) == 14 ) && abs(parentPdg) == 24)
	    {
	      genNeutrino.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genNeutrino.push_back(genNeutrino);
	    }
	  if( (abs(genloop->pdgId) == 1 || abs(genloop->pdgId) == 3 || abs(genloop->pdgId) == 5 || abs(genloop->pdgId) == 2 || abs(genloop->pdgId) == 4 || abs(genloop->pdgId) == 6) && abs(parentPdg) == 24)
	    {
	      genWquarks.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genWquarks.push_back(genWquarks);
	    }
	  if( (abs(genloop->pdgId) == 1 || abs(genloop->pdgId) == 3 || abs(genloop->pdgId) == 5 || abs(genloop->pdgId) == 2 || abs(genloop->pdgId) == 4 || abs(genloop->pdgId) == 6) && genloop->status == 23 && abs(parentPdg) != 24)
	    {
	      genVBFquarks.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genVBFquarks.push_back(genVBFquarks);
	    }
	}
      if (v_genLep.size()==1 && v_genNeutrino.size()==1 && v_genWquarks.size()==2 && v_genVBFquarks.size()==2)
	{
	  //cout<<"====found WWJJ event.... \n\n"<<endl;
	  WWTree->isGen           = 1;
	  WWTree->lep_pt_gen      = v_genLep[0].Pt();
	  WWTree->lep_eta_gen     = v_genLep[0].Eta();
	  WWTree->hadW_pt_gen     = (v_genWquarks[0]+v_genWquarks[1]).Pt();
	  WWTree->hadW_eta_gen    = (v_genWquarks[0]+v_genWquarks[1]).Eta();
	  WWTree->hadW_phi_gen    = (v_genWquarks[0]+v_genWquarks[1]).Phi();
	  WWTree->hadW_e_gen      = (v_genWquarks[0]+v_genWquarks[1]).E();
	  WWTree->hadW_m_gen      = (v_genWquarks[0]+v_genWquarks[1]).M();
	  
	  WWTree->lepW_pt_gen     = (v_genLep[0]+v_genNeutrino[0]).Pt();
	  WWTree->lepW_eta_gen    = (v_genLep[0]+v_genNeutrino[0]).Eta();
	  WWTree->lepW_phi_gen    = (v_genLep[0]+v_genNeutrino[0]).Phi();
	  WWTree->lepW_e_gen      = (v_genLep[0]+v_genNeutrino[0]).E();
	  WWTree->lepW_m_gen      = (v_genLep[0]+v_genNeutrino[0]).M();

	  WWTree->WW_mass_gen	= (v_genLep[0] + v_genNeutrino[0] + v_genWquarks[0] + v_genWquarks[1]).M();
	  WWTree->WW_mT_gen	= (v_genLep[0] + v_genNeutrino[0] + v_genWquarks[0] + v_genWquarks[1]).Mt();
	  WWTree->WW_pT_gen	= (v_genLep[0] + v_genNeutrino[0] + v_genWquarks[0] + v_genWquarks[1]).Pt();

	  WWTree->AK4_1_pt_gen	= v_genVBFquarks[0].Pt();
	  WWTree->AK4_1_eta_gen	= v_genVBFquarks[0].Eta();
	  WWTree->AK4_1_phi_gen	= v_genVBFquarks[0].Phi();
	  WWTree->AK4_1_e_gen	= v_genVBFquarks[0].E();
	  WWTree->AK4_1_mass_gen	= v_genVBFquarks[0].M();
	  
	  WWTree->AK4_2_pt_gen	= v_genVBFquarks[1].Pt();
	  WWTree->AK4_2_eta_gen	= v_genVBFquarks[1].Eta();
	  WWTree->AK4_2_phi_gen	= v_genVBFquarks[1].Phi();
	  WWTree->AK4_2_e_gen	= v_genVBFquarks[1].E();
	  WWTree->AK4_2_mass_gen	= v_genVBFquarks[1].M();
	  
	  WWTree->AK4_jj_DeltaEta_gen = abs(v_genVBFquarks[0].Eta() - v_genVBFquarks[1].Eta());
	  WWTree->AK4_jj_mass_gen = (v_genVBFquarks[0] + v_genVBFquarks[1]).M();
	  
	  count_genEvents++;
	}
	if ( WWTree->lep_pt_gen > 30 && abs(WWTree->lep_eta_gen) < 2.5 && WWTree->AK4_jj_DeltaEta_gen>2 && WWTree->hadW_pt_gen>200 && WWTree->AK4_1_pt_gen >30 && WWTree->AK4_2_pt_gen>30 && WWTree->AK4_jj_mass_gen > 500 ){
	        GenPassCut = 1;
	 }
    	for (int i = 0; i<lheWgtArr->GetEntries();i++)     // Note that i is starting from 446.
	{
		const baconhep::TLHEWeight *lhe = (baconhep::TLHEWeight*)((*lheWgtArr)[i]);
		WWTree->LHEid.push_back(lhe->id);
		WWTree->LHEWeight.push_back(lhe->weight);
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
    WWTree->trig_eff_Weight = 1.;

    if (gen->weight>0)
      WWTree->genWeight=1.;
    else if (gen->weight<0) {
      WWTree->genWeight=-1.;
      //nNegEvents++;
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
	WWTree->l_eta1 = leadele->eta;
	WWTree->l_phi1 = leadele->phi;	
	WWTree->l_e1 = leadele->ecalEnergy;	
	WWTree->l_charge1 = leadele->q;
      }
    if(subele)
      {
	WWTree->l_pt2  = subele->pt;
	WWTree->l_eta2 = subele->eta;
	WWTree->l_phi2 = subele->phi;	
	WWTree->l_e2 = subele->ecalEnergy;	
	WWTree->l_charge2 = subele->q;
      }
    muonArr->Clear();
    muonBr->GetEntry(jentry);
    const baconhep::TMuon *leadmu = NULL;
    const baconhep::TMuon *submu = NULL;
    double leadmue=-999, submue = -999;
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
      if(!leadmu || mu->pt>leadmu->pt)
	{
	  if(!(mu->pt>leadmupt_cut)) continue;
	  submu = leadmu;
	  leadmu = mu;
	  leadmue = MU.E();
	}
      else if (!submu || mu->pt > submu->pt)
	{
	  submu = mu;
	  submue = MU.E();
	}
    }   
    if(leadmu)
      {
	WWTree->l_pt1  = leadmu->pt;
	WWTree->l_eta1 = leadmu->eta;
	WWTree->l_phi1 = leadmu->phi;	
	WWTree->l_e1 = leadmue;	
	WWTree->l_charge1 = leadmu->q;
      }
    if(submu)
      {
	WWTree->l_pt2  = submu->pt;
	WWTree->l_eta2 = submu->eta;
	WWTree->l_phi2 = submu->phi;	
	WWTree->l_e2 = submue;	
	WWTree->l_charge2 = submu->q;
      }
    //std::cout << nTightMu << " " << nTightEle << " " << " " << nLooseEle << " " << nLooseMu << " " << std::endl; 
    //std::cout << leadele << " " << subele << " " << leadmu << " " << submu << std::endl;
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
    
    LEP.SetPtEtaPhiE(WWTree->l_pt1,WWTree->l_eta1,WWTree->l_phi1,WWTree->l_e1);
    LEP2.SetPtEtaPhiE(WWTree->l_pt2,WWTree->l_eta2,WWTree->l_phi2,WWTree->l_e2);
    if(WWTree->l_pt2>0)
      {
	WWTree->dilep_pt  = (LEP+LEP2).Pt();
	WWTree->dilep_eta = (LEP+LEP2).Eta();
	WWTree->dilep_phi = (LEP+LEP2).Phi();	
	WWTree->dilep_m = (LEP+LEP2).M();	
      }

    //ID & GSF efficiency SF for electrons (https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale)
    if (strcmp(leptonName.c_str(),"el")==0 && isMC==1) {
	//  apply ID, ISO SF's
	WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hIDIsoEle);	// Get Scale factor corresponding to the pt and eta.

	// apply GSF/RECO SF's for electrons
	WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hGSFCorrEle);
	
	WWTree->trig_eff_Weight = WWTree->trig_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hTriggerEle);
    }

    //ID&ISO efficiency SF for muons (https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults)
    if (strcmp(leptonName.c_str(),"mu")==0 && isMC==1) {
	//  apply ID SF's
	if (WWTree->run<278820)
		WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIDMuA);
	else
		WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIDMuB);

	//  apply ISO SF's
	if (WWTree->run<278820)
		WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIsoMuA);
	else
		WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIsoMuB);
	
	// apply Trigger SF's
	if (WWTree->run<278820)
		WWTree->trig_eff_Weight = WWTree->trig_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hTriggerMuA);
	else
		WWTree->trig_eff_Weight = WWTree->trig_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hTriggerMuB);
    }
	
    //////////////THE MET
    
    // //preselection on met
     //if (info->pfMET < 0) continue;
     cutEff[3]++;
    
    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*
    
    float Wmass = 80.385;
  
    TLorentzVector W_Met, W_Met_jes_up, W_Met_jes_dn;
  
    W_Met.SetPxPyPzE(info->pfMET * TMath::Cos(info->pfMETphi), info->pfMET * TMath::Sin(info->pfMETphi), 0., sqrt(info->pfMET*info->pfMET));
    //W_Met_jes_up.SetPxPyPzE(ReducedTree->METPtUp * TMath::Cos(ReducedTree->METPhiUp), ReducedTree->METPtUp * TMath::Sin(ReducedTree->METPhiUp), 0., sqrt(ReducedTree->METPtUp*ReducedTree->METPtUp));
    //W_Met_jes_dn.SetPxPyPzE(ReducedTree->METPtDown * TMath::Cos(ReducedTree->METPhiDown), ReducedTree->METPtDown * TMath::Sin(ReducedTree->METPhiDown), 0., sqrt(ReducedTree->METPtDown*ReducedTree->METPtDown));
  
    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    METzCalculator NeutrinoPz_type0_jes_up;
    METzCalculator NeutrinoPz_type0_jes_dn;
    METzCalculator_Run2 NeutrinoPz_run2;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(LEP);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_type0_jes_up.SetMET(W_Met_jes_up);
    NeutrinoPz_type0_jes_up.SetLepton(LEP);
    NeutrinoPz_type0_jes_up.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_type0_jes_dn.SetMET(W_Met_jes_dn);
    NeutrinoPz_type0_jes_dn.SetLepton(LEP);
    NeutrinoPz_type0_jes_dn.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(LEP);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());
  
    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    //double pz2_type0 = NeutrinoPz_type0.getOther();  // Default one
    
    double pz1_run2 = NeutrinoPz_run2.Calculate();
  
    //double pz1_type0_jes_up = NeutrinoPz_type0_jes_up.Calculate(); // Default one -> according to type0
    //double pz1_type0_jes_dn = NeutrinoPz_type0_jes_dn.Calculate(); // Default one -> according to type0
  
    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
    
    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0; 
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
    
    if(NeutrinoPz_type0.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(info->pfMETphi), nu_pt1 * TMath::Sin(info->pfMETphi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(info->pfMETphi), nu_pt2 * TMath::Sin(info->pfMETphi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }
    
    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(LEP);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());
    
    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    //double pz2_type2 = NeutrinoPz_type2.getOther();   // Default one
  
    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; 
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    
    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2; 
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    
    if (NeutrinoPz_type2.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(info->pfMETphi), nu_pt1 * TMath::Sin(info->pfMETphi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(info->pfMETphi), nu_pt2 * TMath::Sin(info->pfMETphi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }
    
    WWTree->pfMET = sqrt(info->pfMET*info->pfMET);
    //WWTree->pfMET_jes_up = sqrt(ReducedTree->METPtUp*ReducedTree->METPtUp);
    //WWTree->pfMET_jes_dn = sqrt(ReducedTree->METPtDown*ReducedTree->METPtDown);
    WWTree->pfMET_Phi = info->pfMETphi;
    WWTree->pfMET_Corr = info->pfMETC;
    WWTree->pfMET_Corr_phi = info->pfMETCphi;
    WWTree->pfMET_Corr_Cov00 = info->pfMETCCov00;
    WWTree->pfMET_Corr_Cov01 = info->pfMETCCov01;
    WWTree->pfMET_Corr_Cov11 = info->pfMETCCov11;
    WWTree->pfMET_Corr_jerup = info->pfMETCjerup;
    WWTree->pfMET_Corr_jerdn = info->pfMETCjerdn;
    WWTree->pfMET_Corr_jenup = info->pfMETCjenup;
    WWTree->pfMET_Corr_jendn = info->pfMETCjendn;
    WWTree->pfMET_Corr_uncup = info->pfMETCuncup;
    WWTree->pfMET_Corr_uncdn = info->pfMETCuncdn;
    WWTree->pfMET_Corr_jrsup = info->pfMETCjrsup;
    WWTree->pfMET_Corr_jrsdn = info->pfMETCjrsdn;
    WWTree->pfMET_Corr_phijerup = info->pfMETCphijerup;
    WWTree->pfMET_Corr_phijerdn = info->pfMETCphijerdn;
    WWTree->pfMET_Corr_phijenup = info->pfMETCphijenup;
    WWTree->pfMET_Corr_phijendn = info->pfMETCphijendn;
    WWTree->pfMET_Corr_phiuncup = info->pfMETCphiuncup;
    WWTree->pfMET_Corr_phiuncdn = info->pfMETCphiuncdn;
    WWTree->pfMET_Corr_phijrsup = info->pfMETCphijrsup;
    WWTree->pfMET_Corr_phijrsdn = info->pfMETCphijrsdn;

    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();
  
    
    /////////////////THE LEPTONIC W
  
    NU0.SetPxPyPzE(info->pfMET*TMath::Cos(info->pfMETphi),info->pfMET*TMath::Sin(info->pfMETphi),WWTree->nu_pz_type0,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type0*WWTree->nu_pz_type0));
    //NU0_jes_up.SetPxPyPzE(ReducedTree->METPtUp*TMath::Cos(ReducedTree->METPhiUp),ReducedTree->METPtUp*TMath::Sin(ReducedTree->METPhiUp),pz1_type0_jes_up,TMath::Sqrt(WWTree->pfMET_jes_up*WWTree->pfMET_jes_up+pz1_type0_jes_up*pz1_type0_jes_up));
    //NU0_jes_dn.SetPxPyPzE(ReducedTree->METPtDown*TMath::Cos(ReducedTree->METPhiDown),ReducedTree->METPtDown*TMath::Sin(ReducedTree->METPhiDown),pz1_type0_jes_dn,TMath::Sqrt(WWTree->pfMET_jes_dn*WWTree->pfMET_jes_dn+pz1_type0_jes_dn*pz1_type0_jes_dn));
    
    NU2.SetPxPyPzE(info->pfMET*TMath::Cos(info->pfMETphi),info->pfMET*TMath::Sin(info->pfMETphi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    NU1.SetPxPyPzE(info->pfMET*TMath::Cos(info->pfMETphi),info->pfMET*TMath::Sin(info->pfMETphi),WWTree->nu_pz_run2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_run2*WWTree->nu_pz_run2));

  
    W = LEP + NU2;
  
    WWTree->v_pt = W.Pt();
    WWTree->v_eta = W.Eta();
    WWTree->v_phi = W.Phi();
    WWTree->v_mt = TMath::Sqrt(2*LEP.Et()*NU2.Et()*(1-TMath::Cos(LEP.DeltaPhi(NU2))));
    WWTree->v_mass = W.M();
  
    //////////////THE PUPPI MET
  
    //preselection on met
    //if(info->pfMET<0 && info->puppET<0) continue;
    //cutEff[3]++;
    
    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*
    
    W_Met.SetPxPyPzE(info->puppET * TMath::Cos(info->puppETphi), info->puppET * TMath::Sin(info->puppETphi), 0., sqrt(info->puppET*info->puppET));
    //W_Met_jes_up.SetPxPyPzE(info->METpuppiPtUp * TMath::Cos(info->METpuppiPhiUp), info->METpuppiPtUp * TMath::Sin(info->METpuppiPhiUp), 0., sqrt(info->METpuppiPtUp*info->METpuppiPtUp));
    //W_Met_jes_dn.SetPxPyPzE(info->METpuppiPtDown * TMath::Cos(info->METpuppiPhiDown), info->METpuppiPtDown * TMath::Sin(info->METpuppiPhiDown), 0., sqrt(info->METpuppiPtDown*info->METpuppiPtDown));
    
    if(LEP.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }
    cutEff[4]++;
  
    // type0 calculation of neutrino pZ
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(LEP);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());
    
    //NeutrinoPz_type0_jes_up.SetMET(W_Met_jes_up);
    //NeutrinoPz_type0_jes_up.SetLepton(LEP);
    //NeutrinoPz_type0_jes_up.SetLeptonType(leptonName.c_str());
    //
    //NeutrinoPz_type0_jes_dn.SetMET(W_Met_jes_dn);
    //NeutrinoPz_type0_jes_dn.SetLepton(LEP);
    //NeutrinoPz_type0_jes_dn.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(LEP);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());
    
    pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    // pz2_type0 = NeutrinoPz_type0.getOther();  // Default one
    
    pz1_run2 = NeutrinoPz_run2.Calculate();
    
    //pz1_type0_jes_up = NeutrinoPz_type0_jes_up.Calculate(); // Default one -> according to type0
    //pz1_type0_jes_dn = NeutrinoPz_type0_jes_dn.Calculate(); // Default one -> according to type0
    
    // don't touch the neutrino pT
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
  
    // change the neutrino pT in case of complex solution in order to make it real
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
  
    if (NeutrinoPz_type0.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(info->puppETphi), nu_pt1 * TMath::Sin(info->puppETphi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(info->puppETphi), nu_pt2 * TMath::Sin(info->puppETphi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }
  
    // type2 calculation of neutrino pZ
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(LEP);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());
    
    pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    //double pz2_type2 = NeutrinoPz_type2.getOther();   // Default one
  
    // don't touch the neutrino pT
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
  
    // change the neutrino pT in case of complex solution in order to make it real
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    
    if (NeutrinoPz_type2.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(info->puppETphi), nu_pt1 * TMath::Sin(info->puppETphi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(info->puppETphi), nu_pt2 * TMath::Sin(info->puppETphi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }
    
    WWTree->pfMETpuppi = sqrt(info->puppET*info->puppET);
    //WWTree->pfMETpuppi_jes_up = sqrt(info->METpuppiPtUp*info->METpuppiPtUp);
    //WWTree->pfMETpuppi_jes_dn = sqrt(info->METpuppiPtDown*info->METpuppiPtDown);
    WWTree->pfMETpuppi_Phi = info->puppETphi;
    WWTree->pfMETpuppi_Cov00 = info->puppETCov00;
    WWTree->pfMETpuppi_Cov01 = info->puppETCov01;
    WWTree->pfMETpuppi_Cov11 = info->puppETCov11;
    WWTree->pfMETpuppi_Corr = info->puppETC;
    WWTree->pfMETpuppi_Corr_phi = info->puppETCphi;
    WWTree->pfMETpuppi_Corr_Cov00 = info->puppETCCov00;
    WWTree->pfMETpuppi_Corr_Cov01 = info->puppETCCov01;
    WWTree->pfMETpuppi_Corr_Cov11 = info->puppETCCov11;
    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();

    
    /////////////////THE LEPTONIC W PUPPI
    
    NU0_puppi.SetPxPyPzE(info->puppET*TMath::Cos(info->puppETphi),info->puppET*TMath::Sin(info->puppETphi),WWTree->nu_pz_type0,TMath::Sqrt(WWTree->pfMETpuppi*WWTree->pfMETpuppi+WWTree->nu_pz_type0*WWTree->nu_pz_type0));
    //NU0_jes_up.SetPxPyPzE(info->METpuppiPtUp*TMath::Cos(info->METpuppiPhiUp),info->METpuppiPtUp*TMath::Sin(info->METpuppiPhiUp),pz1_type0_jes_up,TMath::Sqrt(WWTree->pfMETpuppi_jes_up*WWTree->pfMETpuppi_jes_up+pz1_type0_jes_up*pz1_type0_jes_up));
    //NU0_jes_dn.SetPxPyPzE(info->METpuppiPtDown*TMath::Cos(info->METpuppiPhiDown),info->METpuppiPtDown*TMath::Sin(info->METpuppiPhiDown),pz1_type0_jes_dn,TMath::Sqrt(WWTree->pfMETpuppi_jes_dn*WWTree->pfMETpuppi_jes_dn+pz1_type0_jes_dn*pz1_type0_jes_dn));
    
    NU2_puppi.SetPxPyPzE(info->puppET*TMath::Cos(info->puppETphi),info->puppET*TMath::Sin(info->puppETphi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMETpuppi*WWTree->pfMETpuppi+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    NU1_puppi.SetPxPyPzE(info->puppET*TMath::Cos(info->puppETphi),info->puppET*TMath::Sin(info->puppETphi),WWTree->nu_pz_run2,TMath::Sqrt(WWTree->pfMETpuppi*WWTree->pfMETpuppi+WWTree->nu_pz_run2*WWTree->nu_pz_run2));
  
    W_puppi = LEP + NU2_puppi;
    WWTree->v_puppi_pt = W_puppi.Pt();
    WWTree->v_puppi_eta = W_puppi.Eta();
    WWTree->v_puppi_phi = W_puppi.Phi();
    WWTree->v_puppi_mt = TMath::Sqrt(2*LEP.Et()*NU2_puppi.Et()*(1-TMath::Cos(LEP.DeltaPhi(NU2_puppi))));
    WWTree->v_puppi_mass = W_puppi.M();
      
    ///////////THE FAT JET - AK8
    //float tempPt=0.;
    float tempTTbarMass=0.;
    float tempMassW = 3000.0;
    //int nGoodAK8jets=0;
    //int ttb_jet_position=-1; //position of AK8 jet in ttbar-topology
    vjetArr->Clear();
    vjetBr->GetEntry(jentry);
    vjetAddArr->Clear();
    vjetAddBr->GetEntry(jentry);
    for ( int i=0; i<vjetArr->GetEntries(); i++)
    {
	const baconhep::TJet *jet = (baconhep::TJet*)((*vjetArr)[i]);
	const baconhep::TAddJet *addjet = (baconhep::TAddJet*)((*vjetAddArr)[i]);
	TLorentzVector TempAK8;
	TempAK8.SetPtEtaPhiM(jet->pt,fabs(jet->eta),jet->phi,jet->mass);
	bool isCleanedJet = true;
	//if (jet->pt<200 )  continue;
	if (jet->pt<200 || fabs(jet->eta)>2.4)  continue;
	if (addjet->mass_prun>tempTTbarMass) {
	  if ( (jet->eta>0 && WWTree->l_eta1<0) || 
	      (jet->eta<0 && WWTree->l_eta1>0)) { //jet and lepton in opposite hemisphere for ttb
	    //ttb_jet_position=i; //save AK8 jet in ttbar topology
	    tempTTbarMass=addjet->mass_prun;
	  }
	}
	//if (jet->pt<=tempPt) continue; //save the jet with the largest pt
	//if (abs(addjet->mass_sd0 - 80.385) > tempMassW ) continue; //save the jet which is closest to W-mass
	if (abs(jet->mass - 80.385) > tempMassW ) continue; //save the jet which is closest to W-mass
	//if (ReducedTree->AK8Jets_AK8isLooseJetId[i]==false) continue; //fat jet must satisfy loose ID
	
	//CLEANING FROM LEPTONS
	for (std::size_t j=0; j<tightEle.size(); j++) {
	  if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
		     jet->eta, jet->phi) < 1.0)
	    isCleanedJet = false;
	}
	for ( std::size_t j=0; j<tightMuon.size(); j++) {
	  if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
		     jet->eta, jet->phi) < 1.0)
	    isCleanedJet = false;
	}
      
	if (isCleanedJet==false) continue; //jet is overlapped with a lepton
      
	WWTree->ungroomed_AK8jet_pt  = jet->pt;
	WWTree->ungroomed_AK8jet_eta = jet->eta;
	WWTree->ungroomed_AK8jet_phi = jet->phi;
	WWTree->ungroomed_AK8jet_e   = TempAK8.E();
	//WWTree->ungroomed_jet_pt_jes_up = (jet->pt[i]/ReducedTree->AK8Jets_AK8correction[i])*ReducedTree->AK8Jets_AK8correctionUp[i];
	//WWTree->ungroomed_jet_pt_jes_dn = (ReducedTree->AK8CHS_pt[i]/ReducedTree->AK8Jets_AK8correction[i])*ReducedTree->AK8Jets_AK8correctionDown[i];
      
	WWTree->AK8jet_mass     = jet->mass;
	WWTree->AK8jet_mass_pr  = addjet->mass_prun;
	WWTree->AK8jet_mass_so  = addjet->mass_sd0;
	WWTree->AK8jet_mass_tr  = addjet->mass_trim;
	WWTree->AK8jet_tau2tau1 = addjet->tau2/addjet->tau1;
	WWTree->AK8_jetID_loose = passJetLooseSel(jet);
	//WWTree->jet_mass_pr_jes_up = (ReducedTree->AddAK8CHS_mass_prun[i]/ReducedTree->AK8Jets_AK8massCorrection[i])*ReducedTree->AK8Jets_AK8massCorrectionUp[i];
	//WWTree->jet_mass_pr_jes_dn = (ReducedTree->AddAK8CHS_mass_prun[i]/ReducedTree->AK8Jets_AK8massCorrection[i])*ReducedTree->AK8Jets_AK8massCorrectionDown[i];
      
	//tempPt = WWTree->ungroomed_jet_pt;
	//tempMassW = abs(addjet->mass_sd0 - 80.385);
	tempMassW = abs(jet->mass - 80.385);
	WWTree->nGoodAK8jets++;
      }
    if(WWTree->ungroomed_AK8jet_pt > 0)
      {
	JET.SetPtEtaPhiE(WWTree->ungroomed_AK8jet_pt,WWTree->ungroomed_AK8jet_eta,WWTree->ungroomed_AK8jet_phi,WWTree->ungroomed_AK8jet_e);
	/* JET_jes_up.SetPtEtaPhiE(WWTree->ungroomed_jet_pt*(ReducedTree->AK8Jets_AK8correctionUp[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]),
	   WWTree->ungroomed_jet_eta,
	   WWTree->ungroomed_jet_phi,
	   WWTree->ungroomed_jet_e*(ReducedTree->AK8Jets_AK8correctionUp[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]));
	   JET_jes_dn.SetPtEtaPhiE(WWTree->ungroomed_jet_pt*(ReducedTree->AK8Jets_AK8correctionDown[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]),
	   WWTree->ungroomed_jet_eta,
	   WWTree->ungroomed_jet_phi,
	   WWTree->ungroomed_jet_e*(ReducedTree->AK8Jets_AK8correctionDown[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]));
	*/
      }
    
    
    ///////////THE FAT JET - PuppiAK8
    vjetArrPuppi->Clear();
    vjetBrPuppi->GetEntry(jentry);
    vjetAddArrPuppi->Clear();
    vjetAddBrPuppi->GetEntry(jentry);
    tempTTbarMass=0.;
    tempMassW = 3000.0;
    int nGoodPuppiAK8jets=0;
    //int hadWPuppiAK8pos = -1;
    int ttb_PuppiAK8_jet_position=-1; //position of AK8 jet in ttbar-topology
    
    for ( int i=0; i<vjetArrPuppi->GetEntries(); i++)
      {
	const baconhep::TJet *jet = (baconhep::TJet*)((*vjetArrPuppi)[i]);
	const baconhep::TAddJet *addjet = (baconhep::TAddJet*)((*vjetAddArrPuppi)[i]);
	TLorentzVector TempAK8;
	TempAK8.SetPtEtaPhiM(jet->pt,fabs(jet->eta),jet->phi,jet->mass);
	bool isCleanedJet = true;
	//if (jet->pt<200)  continue; //be careful: this is not inside the synchntuple code
	if (jet->pt<200 || fabs(jet->eta)>2.4)  continue; //be careful: this is not inside the synchntuple code
        //if (!passJetLooseSel(jet)) continue;
	if (addjet->mass_prun>tempTTbarMass) {
	  if ( (jet->eta>0 && WWTree->l_eta1<0) || 
	      (jet->eta<0 && WWTree->l_eta1>0)) { //jet and lepton in opposite hemisphere for ttb
	    ttb_PuppiAK8_jet_position=i; //save AK8 jet in ttbar topology
	    tempTTbarMass=addjet->mass_prun;
	  }
	}
	//if (jet->pt<=tempPt) continue; //save the jet with the largest pt
	if (abs(jet->mass - 80.385) > tempMassW) continue; //save the jet closest to w-mass
	//if (ReducedTree->PuppiAK8Jets_PuppiAK8isLooseJetId[i]==false) continue; //fat jet must satisfy loose ID
	
	//CLEANING FROM LEPTONS
	for ( std::size_t j=0; j<tightEle.size(); j++) {
	  if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
		     jet->eta, jet->phi) < 1.0)
	    isCleanedJet = false;
	}
	for ( std::size_t j=0; j<tightMuon.size(); j++) {
	  if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
		     jet->eta, jet->phi) < 1.0)
	    isCleanedJet = false;
	}
      
	if (isCleanedJet==false) continue; //jet is overlapped with a lepton
	
	WWTree->ungroomed_PuppiAK8_jet_pt  = jet->pt;
	WWTree->ungroomed_PuppiAK8_jet_eta = jet->eta;
	WWTree->ungroomed_PuppiAK8_jet_phi = jet->phi;
	WWTree->ungroomed_PuppiAK8_jet_e   = TempAK8.E();
	//WWTree->ungroomed_PuppiAK8_jet_pt_jes_up = (jet->pt[i]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[i])*ReducedTree->PuppiAK8Jets_PuppiAK8correctionUp[i];
	//WWTree->ungroomed_PuppiAK8_jet_pt_jes_dn = (jet->pt[i]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[i])*ReducedTree->PuppiAK8Jets_PuppiAK8correctionDown[i];
      
	WWTree->PuppiAK8_jet_mass  = jet->mass;
	WWTree->PuppiAK8_jet_mass_pr  = addjet->mass_prun;
	WWTree->PuppiAK8_jet_mass_so  = addjet->mass_sd0;
	WWTree->PuppiAK8_jet_mass_tr  = addjet->mass_trim;
	WWTree->PuppiAK8_jet_tau2tau1 = addjet->tau2/addjet->tau1;
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
	//     WWTree->PuppiAK8_jet_mass_pr_jes_up = (addjet->mass_prun[i]/ReducedTree->PuppiAK8Jets_PuppiAK8massCorrection[i])*ReducedTree->PuppiAK8Jets_PuppiAK8massCorrectionUp[i];
	//   WWTree->PuppiAK8_jet_mass_pr_jes_dn = (addjet->mass_prun[i]/ReducedTree->PuppiAK8Jets_PuppiAK8massCorrection[i])*ReducedTree->PuppiAK8Jets_PuppiAK8massCorrectionDown[i];
	
	//tempPt = WWTree->ungroomed_PuppiAK8_jet_pt;
	tempMassW = abs(jet->mass - 80.385);
	nGoodPuppiAK8jets++;
	//hadWPuppiAK8pos = i;
      }
    if (WWTree->ungroomed_PuppiAK8_jet_pt > 0.)
      {
	JET_PuppiAK8.SetPtEtaPhiE(WWTree->ungroomed_PuppiAK8_jet_pt,WWTree->ungroomed_PuppiAK8_jet_eta,WWTree->ungroomed_PuppiAK8_jet_phi,WWTree->ungroomed_PuppiAK8_jet_e);
	SJ1.SetPtEtaPhiM(WWTree->PuppiAK8_jet_sj1_pt, WWTree->PuppiAK8_jet_sj1_eta, WWTree->PuppiAK8_jet_sj1_phi, WWTree->PuppiAK8_jet_sj1_m);
	SJ2.SetPtEtaPhiM(WWTree->PuppiAK8_jet_sj2_pt, WWTree->PuppiAK8_jet_sj2_eta, WWTree->PuppiAK8_jet_sj2_phi, WWTree->PuppiAK8_jet_sj2_m);
	/* 
	   JET_PuppiAK8_jes_up.SetPtEtaPhiE(WWTree->ungroomed_PuppiAK8_jet_pt*(ReducedTree->PuppiAK8Jets_PuppiAK8correctionUp[hadWPuppiAK8pos]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[hadWPuppiAK8pos]),
	   WWTree->ungroomed_PuppiAK8_jet_eta,
	   WWTree->ungroomed_PuppiAK8_jet_phi,
	   WWTree->ungroomed_PuppiAK8_jet_e*(ReducedTree->PuppiAK8Jets_PuppiAK8correctionUp[hadWPuppiAK8pos]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[hadWPuppiAK8pos]));
	   JET_PuppiAK8_jes_dn.SetPtEtaPhiE(WWTree->ungroomed_PuppiAK8_jet_pt*(ReducedTree->PuppiAK8Jets_PuppiAK8correctionDown[hadWPuppiAK8pos]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[hadWPuppiAK8pos]),
	   WWTree->ungroomed_PuppiAK8_jet_eta,
	   WWTree->ungroomed_PuppiAK8_jet_phi,
	   WWTree->ungroomed_PuppiAK8_jet_e*(ReducedTree->PuppiAK8Jets_PuppiAK8correctionDown[hadWPuppiAK8pos]/ReducedTree->PuppiAK8Jets_PuppiAK8correction[hadWPuppiAK8pos]));
	*/
      }
    
    
    // FAT JET SELECTION
    bool isGoodFatJet = true;
    if (WWTree->nGoodAK8jets==0 && nGoodPuppiAK8jets==0) isGoodFatJet = false; //not found a good hadronic W candidate
    if (WWTree->ungroomed_AK8jet_pt<200 && WWTree->ungroomed_PuppiAK8_jet_pt<200) isGoodFatJet = false;
    if (!isGoodFatJet) continue;
    cutEff[5]++;
    
    
    #if 0
    //////////////////UNMERGED HADRONIC W
    jetArr->Clear();
    jetBr->GetEntry(jentry);
    float tempPt1 = 0.; int pos1 = -1;
    float tempPt2 = 0.; //int pos2 = -1;
    int nGoodAK4jets=0;
    for ( int i=0; i<jetArr->GetEntries(); i++) //loop on AK4 jet
    {
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[i]);
      TLorentzVector TempAK4;
      std::vector<TLorentzVector> TempAK4array;      
      TempAK4.SetPtEtaPhiM(jet->pt,fabs(jet->eta),jet->phi,jet->mass);
      TempAK4array.push_back(TempAK4);

      bool isCleanedJet = true;
      if (jet->pt<=30  || fabs(jet->eta)>=2.4)  continue;
      //if (ReducedTree->Jets_isLooseJetId[i]==false) continue;
      //if (jet->csv>0.890) continue;
      
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   jet->eta, jet->phi) < 0.3) {
          isCleanedJet = false;
        }
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   jet->eta, jet->phi) < 0.3) {
          isCleanedJet = false;
        }
      }
      
      if (isCleanedJet==false) continue;
      
      if (jet->pt<tempPt1 && jet->pt<tempPt2) continue;
      
      if (jet->pt>tempPt1)
      {
        WWTree->AK4_jet1_pt  = jet->pt;
        WWTree->AK4_jet1_eta = jet->eta;
        WWTree->AK4_jet1_phi = jet->phi;
        WWTree->AK4_jet1_e   = TempAK4array[i].E();
        //WWTree->AK4_jet1_pt_jes_up = (jet->pt[i]/ReducedTree->Jets_AK4correction[i])*ReducedTree->Jets_AK4correctionUp[i];
        //WWTree->AK4_jet1_pt_jes_dn = (jet->pt[i]/ReducedTree->Jets_AK4correction[i])*ReducedTree->Jets_AK4correctionDown[i];
        
        tempPt1 = WWTree->AK4_jet1_pt;
        
        if (pos1!=-1)
        {
	  const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArr)[pos1]);
          WWTree->AK4_jet2_pt  = jet2->pt;
          WWTree->AK4_jet2_eta = jet2->eta;
          WWTree->AK4_jet2_phi = jet2->phi;
          WWTree->AK4_jet2_e   = TempAK4array[pos1].E();
         // WWTree->AK4_jet2_pt_jes_up = (jet->pt[pos1]/ReducedTree->Jets_AK4correction[pos1])*ReducedTree->Jets_AK4correctionUp[pos1];
         // WWTree->AK4_jet2_pt_jes_dn = (jet->pt[pos1]/ReducedTree->Jets_AK4correction[pos1])*ReducedTree->Jets_AK4correctionDown[pos1];
          
          tempPt2 = WWTree->AK4_jet2_pt;
        }
        pos1 = i;
        nGoodAK4jets++;
      }
      else if (jet->pt>tempPt2)
      {
        WWTree->AK4_jet2_pt  = jet->pt;
        WWTree->AK4_jet2_eta = jet->eta;
        WWTree->AK4_jet2_phi = jet->phi;
        WWTree->AK4_jet2_e   = TempAK4array[i].E();
       // WWTree->AK4_jet2_pt_jes_up = (jet->pt[i]/ReducedTree->Jets_AK4correction[i])*ReducedTree->Jets_AK4correctionUp[i];
       // WWTree->AK4_jet2_pt_jes_dn = (jet->pt[i]/ReducedTree->Jets_AK4correction[i])*ReducedTree->Jets_AK4correctionDown[i];
        
        tempPt2 = WWTree->AK4_jet2_pt;
        //pos2 = i;
        nGoodAK4jets++;
      }
    }
    
    if (WWTree->AK4_jet1_pt > 0.)
    {
      AK4_JET1.SetPtEtaPhiE(WWTree->AK4_jet1_pt,WWTree->AK4_jet1_eta,WWTree->AK4_jet1_phi,WWTree->AK4_jet1_e);
      //AK4_JET1_jes_up.SetPtEtaPhiE(WWTree->AK4_jet1_pt*(ReducedTree->Jets_AK4correctionUp[pos1]/ReducedTree->Jets_AK4correction[pos1]),
      //                             WWTree->AK4_jet1_eta,
      //                             WWTree->AK4_jet1_phi,
      //                             WWTree->AK4_jet1_e*(ReducedTree->Jets_AK4correctionUp[pos1]/ReducedTree->Jets_AK4correction[pos1]));
      //AK4_JET1_jes_dn.SetPtEtaPhiE(WWTree->AK4_jet1_pt*(ReducedTree->Jets_AK4correctionDown[pos1]/ReducedTree->Jets_AK4correction[pos1]),
      //                             WWTree->AK4_jet1_eta,
      //                             WWTree->AK4_jet1_phi,
      //                             WWTree->AK4_jet1_e*(ReducedTree->Jets_AK4correctionDown[pos1]/ReducedTree->Jets_AK4correction[pos1]));
    }
    if (WWTree->AK4_jet2_pt > 0.)
    {
      AK4_JET2.SetPtEtaPhiE(WWTree->AK4_jet2_pt,WWTree->AK4_jet2_eta,WWTree->AK4_jet2_phi,WWTree->AK4_jet2_e);
      //AK4_JET2_jes_up.SetPtEtaPhiE(WWTree->AK4_jet2_pt*(ReducedTree->Jets_AK4correctionUp[pos2]/ReducedTree->Jets_AK4correction[pos2]),
      //                             WWTree->AK4_jet2_eta,
      //                             WWTree->AK4_jet2_phi,
      //                             WWTree->AK4_jet2_e*(ReducedTree->Jets_AK4correctionUp[pos2]/ReducedTree->Jets_AK4correction[pos2]));
      //AK4_JET2_jes_dn.SetPtEtaPhiE(WWTree->AK4_jet2_pt*(ReducedTree->Jets_AK4correctionDown[pos2]/ReducedTree->Jets_AK4correction[pos2]),
      //                             WWTree->AK4_jet2_eta,
      //                             WWTree->AK4_jet2_phi,
      //                             WWTree->AK4_jet2_e*(ReducedTree->Jets_AK4correctionDown[pos2]/ReducedTree->Jets_AK4correction[pos2]));
    }
    
    if (WWTree->AK4_jet2_pt>0) {
      WWTree->AK4_jetjet_pt = (AK4_JET1+AK4_JET2).Pt();
      WWTree->AK4_jetjet_mass = (AK4_JET1+AK4_JET2).M();
      WWTree->AK4_jetjet_deltaeta = deltaEta(AK4_JET1.Eta(),AK4_JET2.Eta());
      WWTree->AK4_jetjet_deltaphi = deltaPhi(AK4_JET1.Phi(),AK4_JET2.Phi());
      WWTree->AK4_jetjet_deltar = deltaR(AK4_JET1.Eta(),AK4_JET1.Phi(),AK4_JET2.Eta(),AK4_JET2.Phi());
    }

    jetArrPuppi->Clear();
    jetBrPuppi->GetEntry(jentry);
    tempPt1 = 0.; int pos1Puppi = -1;
    tempPt2 = 0.; //int pos2Puppi = -1;
    int nGoodPuppiAK4jets=0;
    for ( int i=0; i<jetArrPuppi->GetEntries(); i++) //loop on PuppiAK4 jet
    {
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArrPuppi)[i]);
      TLorentzVector TempAK4;
      std::vector<TLorentzVector> TempAK4array;      
      TempAK4.SetPtEtaPhiM(jet->pt,fabs(jet->eta),jet->phi,jet->mass);
      TempAK4array.push_back(TempAK4);
      bool isCleanedJet = true;
      if ( jet->pt<=30 || fabs(jet->eta)>=2.4)  continue;
      //if (ReducedTree->JetsPuppi_isLooseJetId[i]==false) continue;
      //if (jet->csv>0.890) continue;
      
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   jet->eta, jet->phi) < 0.3) {
          isCleanedJet = false;
        }
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   jet->eta, jet->phi) < 0.3) {
          isCleanedJet = false;
        }
      }
      
      if (isCleanedJet==false) continue;
      
      if (jet->pt<tempPt1 && jet->pt<tempPt2) continue;
      
      if (jet->pt>tempPt1)
      {
        WWTree->PuppiAK4_jet1_pt  = jet->pt;
        WWTree->PuppiAK4_jet1_eta = jet->eta;
        WWTree->PuppiAK4_jet1_phi = jet->phi;
        WWTree->PuppiAK4_jet1_e   = TempAK4array[i].E();
      //  WWTree->PuppiAK4_jet1_pt_jes_up = (jet->pt[i]/ReducedTree->JetsPuppi_AK4correction[i])*ReducedTree->JetsPuppi_AK4correctionUp[i];
      //  WWTree->PuppiAK4_jet1_pt_jes_dn = (jet->pt[i]/ReducedTree->JetsPuppi_AK4correction[i])*ReducedTree->JetsPuppi_AK4correctionDown[i];
        
        tempPt1 = WWTree->PuppiAK4_jet1_pt;
        
        if (pos1Puppi!=-1)
        {
	  const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArrPuppi)[pos1Puppi]);
          WWTree->PuppiAK4_jet2_pt  = jet2->pt;
          WWTree->PuppiAK4_jet2_eta = jet2->eta;
          WWTree->PuppiAK4_jet2_phi = jet2->phi;
          WWTree->PuppiAK4_jet2_e   = TempAK4array[pos1Puppi].E();
        //  WWTree->PuppiAK4_jet2_pt_jes_up = (jet->pt[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi])*ReducedTree->JetsPuppi_AK4correctionUp[pos1Puppi];
        //  WWTree->PuppiAK4_jet2_pt_jes_dn = (jet->pt[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi])*ReducedTree->JetsPuppi_AK4correctionDown[pos1Puppi];
          
          tempPt2 = WWTree->PuppiAK4_jet2_pt;
        }
        pos1Puppi = i;
        nGoodPuppiAK4jets++;
      }
      else if (jet->pt>tempPt2)
      {
        WWTree->PuppiAK4_jet2_pt  = jet->pt;
        WWTree->PuppiAK4_jet2_eta = jet->eta;
        WWTree->PuppiAK4_jet2_phi = jet->phi;
        WWTree->PuppiAK4_jet2_e   = TempAK4array[i].E();
        //WWTree->PuppiAK4_jet2_pt_jes_up = (jet->pt[i]/ReducedTree->JetsPuppi_AK4correction[i])*ReducedTree->JetsPuppi_AK4correctionUp[i];
        //WWTree->PuppiAK4_jet2_pt_jes_dn = (jet->pt[i]/ReducedTree->JetsPuppi_AK4correction[i])*ReducedTree->JetsPuppi_AK4correctionDown[i];
        
        tempPt2 = WWTree->PuppiAK4_jet2_pt;
        //pos2Puppi = i;
        nGoodPuppiAK4jets++;
      }
    }
    
    if (WWTree->PuppiAK4_jet1_pt > 0.)
    {
      PuppiAK4_JET1.SetPtEtaPhiE(WWTree->PuppiAK4_jet1_pt,WWTree->PuppiAK4_jet1_eta,WWTree->PuppiAK4_jet1_phi,WWTree->PuppiAK4_jet1_e);
      //PuppiAK4_JET1_jes_up.SetPtEtaPhiE(WWTree->PuppiAK4_jet1_pt*(ReducedTree->JetsPuppi_AK4correctionUp[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi]),
      //                                  WWTree->PuppiAK4_jet1_eta,
      //                                  WWTree->PuppiAK4_jet1_phi,
      //                                  WWTree->PuppiAK4_jet1_e*(ReducedTree->JetsPuppi_AK4correctionUp[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi]));
      //PuppiAK4_JET1_jes_dn.SetPtEtaPhiE(WWTree->PuppiAK4_jet1_pt*(ReducedTree->JetsPuppi_AK4correctionDown[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi]),
      //                                  WWTree->PuppiAK4_jet1_eta,
      //                                  WWTree->PuppiAK4_jet1_phi,
      //                                  WWTree->PuppiAK4_jet1_e*(ReducedTree->JetsPuppi_AK4correctionDown[pos1Puppi]/ReducedTree->JetsPuppi_AK4correction[pos1Puppi]));
    }
    if (WWTree->PuppiAK4_jet2_pt > 0.)
    {
      PuppiAK4_JET2.SetPtEtaPhiE(WWTree->PuppiAK4_jet2_pt,WWTree->PuppiAK4_jet2_eta,WWTree->PuppiAK4_jet2_phi,WWTree->PuppiAK4_jet2_e);
      //PuppiAK4_JET2_jes_up.SetPtEtaPhiE(WWTree->PuppiAK4_jet2_pt*(ReducedTree->Jets_AK4correctionUp[pos2Puppi]/ReducedTree->Jets_AK4correction[pos2Puppi]),
      //                                  WWTree->PuppiAK4_jet2_eta,
      //                                  WWTree->PuppiAK4_jet2_phi,
      //                                  WWTree->PuppiAK4_jet2_e*(ReducedTree->JetsPuppi_AK4correctionUp[pos2Puppi]/ReducedTree->JetsPuppi_AK4correction[pos2Puppi]));
      //PuppiAK4_JET2_jes_dn.SetPtEtaPhiE(WWTree->PuppiAK4_jet2_pt*(ReducedTree->JetsPuppi_AK4correctionDown[pos2Puppi]/ReducedTree->JetsPuppi_AK4correction[pos2Puppi]),
      //                                  WWTree->PuppiAK4_jet2_eta,
      //                                  WWTree->PuppiAK4_jet2_phi,
      //                                  WWTree->PuppiAK4_jet2_e*(ReducedTree->JetsPuppi_AK4correctionDown[pos2Puppi]/ReducedTree->JetsPuppi_AK4correction[pos2Puppi]));
    }
    
    if (WWTree->PuppiAK4_jet2_pt>0) {
      WWTree->PuppiAK4_jetjet_pt = (PuppiAK4_JET1+PuppiAK4_JET2).Pt();
      WWTree->PuppiAK4_jetjet_mass = (PuppiAK4_JET1+PuppiAK4_JET2).M();
      WWTree->PuppiAK4_jetjet_deltaeta = deltaEta(PuppiAK4_JET1.Eta(),PuppiAK4_JET2.Eta());
      WWTree->PuppiAK4_jetjet_deltaphi = deltaPhi(PuppiAK4_JET1.Phi(),PuppiAK4_JET2.Phi());
      WWTree->PuppiAK4_jetjet_deltar = deltaR(PuppiAK4_JET1.Eta(),PuppiAK4_JET1.Phi(),PuppiAK4_JET2.Eta(),PuppiAK4_JET2.Phi());
    }
    
    
    // UNMERGED JETS SELECTION
    bool isGoodUnmergedJets = true;
    if (nGoodAK4jets<2 && nGoodPuppiAK4jets<2) isGoodUnmergedJets = false; //not found a good hadronic W candidate
    if (WWTree->AK4_jetjet_pt<150 && WWTree->PuppiAK4_jetjet_pt<150) isGoodUnmergedJets = false;
    if (isGoodUnmergedJets) cutEff[6]++;
    
    if (isGoodFatJet && isGoodUnmergedJets)    cutEff[7]++;
    #endif
    
    
    //////////////////ANGULAR VARIABLES
    WWTree->deltaR_Wjet_GenReco = deltaR(WWTree->hadW_eta_gen,WWTree->hadW_phi_gen,JET.Eta(),JET.Phi());
    WWTree->deltaR_lak8jet = deltaR(JET.Eta(),JET.Phi(),LEP.Eta(),LEP.Phi());
    WWTree->deltaphi_METak8jet = deltaPhi(JET.Phi(),NU2.Phi());
    WWTree->deltaphi_Vak8jet = deltaPhi(JET.Phi(),W.Phi());
    WWTree->deltaR_lPuppiak8jet = deltaR(JET_PuppiAK8.Eta(),JET_PuppiAK8.Phi(),LEP.Eta(),LEP.Phi());
    WWTree->deltaphi_METPuppiak8jet = deltaPhi(JET_PuppiAK8.Phi(),NU2_puppi.Phi());
    WWTree->deltaphi_VPuppiak8jet = deltaPhi(JET_PuppiAK8.Phi(),W_puppi.Phi());
    if (WWTree->deltaR_lak8jet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak8jet)>2.0 && fabs(WWTree->deltaphi_Vak8jet)>2.0 && WWTree->nGoodAK8jets>0)
      WWTree->issignal=1;
    if (WWTree->deltaR_lPuppiak8jet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METPuppiak8jet)>2.0 && fabs(WWTree->deltaphi_VPuppiak8jet)>2.0 && nGoodPuppiAK8jets>0)
      WWTree->issignal_PuppiAK8=1;

    if (WWTree->AK4_jet2_pt>0) {
      WWTree->deltaR_lak4jetjet = deltaR((AK4_JET1+AK4_JET2).Eta(),(AK4_JET1+AK4_JET2).Phi(),LEP.Eta(),LEP.Phi());
      WWTree->deltaphi_METak4jetjet = deltaPhi((AK4_JET1+AK4_JET2).Phi(),NU2.Phi());
      WWTree->deltaphi_Vak4jetjet = deltaPhi((AK4_JET1+AK4_JET2).Phi(),W.Phi());
      if (WWTree->deltaR_lak4jetjet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak4jetjet)>2.0 && fabs(WWTree->deltaphi_Vak4jetjet)>2.0)
        WWTree->issignal_AK4jetjet=1;
    }
    if (WWTree->PuppiAK4_jet2_pt>0) {
      WWTree->deltaR_lPuppiak4jetjet = deltaR((PuppiAK4_JET1+PuppiAK4_JET2).Eta(),(PuppiAK4_JET1+PuppiAK4_JET2).Phi(),LEP.Eta(),LEP.Phi());
      WWTree->deltaphi_METPuppiak4jetjet = deltaPhi((PuppiAK4_JET1+PuppiAK4_JET2).Phi(),NU2.Phi());
      WWTree->deltaphi_VPuppiak4jetjet = deltaPhi((PuppiAK4_JET1+PuppiAK4_JET2).Phi(),W.Phi());
      if (WWTree->deltaR_lPuppiak4jetjet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METPuppiak4jetjet)>2.0 && fabs(WWTree->deltaphi_VPuppiak4jetjet)>2.0)
        WWTree->issignal_PuppiAK4jetjet=1;
    }

    
    
    
    //////////////////FOUR-BODY INVARIANT MASS
    WWTree->mass_lvj_type0 = (LEP + NU0 + JET).M();
    WWTree->mass_lvj_type2 = (LEP + NU2 + JET).M();
    WWTree->mass_lvj_run2  = (LEP + NU1 + JET).M();
    if (isnan(WWTree->mass_lvj_run2) == 1)
    	cout<<"==============> Run2 mass is NAN"<<"\t LEP mass = "<<LEP.M()<<"\tNUmass = "<<NU1.M()<<"\t"<<NU1.Px()<<"\t"<<NU1.Py()<<"\t"<<NU1.Pz()<<"\t"<<NU1.E()<<endl;
    WWTree->mass_lvj_type0_met_jes_up = (LEP + NU0_jes_up + JET_jes_up).M();
    WWTree->mass_lvj_type0_met_jes_dn = (LEP + NU0_jes_dn + JET_jes_dn).M();
    
    WWTree->LepWEta = (LEP + NU0_puppi ).Eta();
    WWTree->LepWRapidity = (LEP + NU0_puppi ).Rapidity();
    WWTree->HadWEta = (JET_PuppiAK8 ).Eta();
    WWTree->HadWRapidity = (JET_PuppiAK8 ).Rapidity();
    WWTree->WWEta = (LEP + NU0_puppi + JET_PuppiAK8 ).Eta();
    WWTree->WWRapidity = (LEP + NU0_puppi + JET_PuppiAK8 ).Rapidity();
    WWTree->pt_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Pt();
    WWTree->pt_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Pt();
    WWTree->pt_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Pt();
    WWTree->eta_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Eta();
    WWTree->eta_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Eta();
    WWTree->eta_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Eta();
    WWTree->rapidity_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Rapidity();
    WWTree->rapidity_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Rapidity();
    WWTree->rapidity_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Rapidity();
    WWTree->phi_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Phi();
    WWTree->phi_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Phi();
    WWTree->phi_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Phi();
    WWTree->energy_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).E();
    WWTree->energy_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).E();
    WWTree->energy_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).E();
    WWTree->mass_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).M();
    WWTree->mass_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).M();
    WWTree->mass_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).M();
    WWTree->mt_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).Mt();
    WWTree->mt_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).Mt();
    WWTree->mt_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).Mt();
    WWTree->mass_lvj_type0_met_PuppiAK8_jes_up = (LEP + NU0_jes_up + JET_PuppiAK8_jes_up).M();
    WWTree->mass_lvj_type0_met_PuppiAK8_jes_dn = (LEP + NU0_jes_dn + JET_PuppiAK8_jes_dn).M();
    
    WWTree->mass_lvjj_type0_AK4 = (LEP + NU0 + AK4_JET1 + AK4_JET2).M();
    WWTree->mass_lvjj_type2_AK4 = (LEP + NU2 + AK4_JET1 + AK4_JET2).M();
    WWTree->mass_lvjj_run2_AK4  = (LEP + NU1 + AK4_JET1 + AK4_JET2).M();
    WWTree->mass_lvjj_type0_met_jes_up_AK4 = (LEP + NU0_jes_up + AK4_JET1_jes_up + AK4_JET2_jes_up).M();
    WWTree->mass_lvjj_type0_met_jes_dn_AK4 = (LEP + NU0_jes_dn + AK4_JET1_jes_dn + AK4_JET2_jes_dn).M();
    
    WWTree->mass_lvjj_type0_PuppiAK4 = (LEP + NU0 + PuppiAK4_JET1 + PuppiAK4_JET2).M();
    WWTree->mass_lvjj_type2_PuppiAK4 = (LEP + NU2 + PuppiAK4_JET1 + PuppiAK4_JET2).M();
    WWTree->mass_lvjj_run2_PuppiAK4  = (LEP + NU1 + PuppiAK4_JET1 + PuppiAK4_JET2).M();
    WWTree->mass_lvjj_type0_met_jes_up_PuppiAK4 = (LEP + NU0_jes_up + PuppiAK4_JET1_jes_up + PuppiAK4_JET2_jes_up).M();
    WWTree->mass_lvjj_type0_met_jes_dn_PuppiAK4 = (LEP + NU0_jes_dn + PuppiAK4_JET1_jes_dn + PuppiAK4_JET2_jes_dn).M();
    WWTree->mass_llj_PuppiAK8 = (LEP + LEP2 + JET_PuppiAK8).M();

   // if (WWTree->mass_lvj_type0_PuppiAK8 < 500 || WWTree->mass_lvj_type2_PuppiAK8 < 500 || WWTree->mass_lvj_run2_PuppiAK8 < 500) continue;    // cut on WW invariant mass 
    cutEff[8]++;
    

    
    //--- ttbar topology ------
    if (ttb_PuppiAK8_jet_position>=0)
    {
      const baconhep::TJet *jet	= (baconhep::TJet*)((*vjetArrPuppi)[ttb_PuppiAK8_jet_position]);
      const baconhep::TAddJet *addjet = (baconhep::TAddJet*)((*vjetAddArrPuppi)[ttb_PuppiAK8_jet_position]);
      WWTree->ttb_ungroomed_jet_pt  = jet->pt;
      WWTree->ttb_ungroomed_jet_eta = jet->eta;
      WWTree->ttb_ungroomed_jet_phi = jet->phi;
      WWTree->ttb_jet_mass_pr       = addjet->mass_prun;
      WWTree->ttb_jet_mass_so       = addjet->mass_sd0;
      WWTree->ttb_jet_mass_tr       = addjet->mass_trim;
      WWTree->ttb_jet_tau2tau1      = addjet->tau2/addjet->tau1;
      WWTree->ttb_deltaeta_lak8jet = deltaEta(WWTree->ttb_ungroomed_jet_eta,WWTree->l_eta1);
    }
    
    
    /////////VBF and b-tag section
    WWTree->njets=0;
    WWTree->nBTagJet_loose=0;
    WWTree->nBTagJet_medium=0;
    WWTree->nBTagJet_tight=0;
    
    WWTree->njets_unmerged=0;
    WWTree->nBTagJet_loose_unmerged=0;
    WWTree->nBTagJet_medium_unmerged=0;
    WWTree->nBTagJet_tight_unmerged=0;

    float oldDeltaR = 1000.;
    float oldDeltaRLep = 1000.;
    int indexCloserJet = -1;
    int indexCloserJetLep = -1;
    float deltaRbtag_prev=100.;
    float deltaRbtag_prev_loose=100.;
    
    #if 0
    std::vector<int> indexGoodVBFJets;


    for ( int i=0; i<jetArr->GetEntries(); i++) //loop on AK4 jet
    {
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[i]);
      bool isCleaned = true;
      bool isCleanedFromFatJet = true;
      bool isCleanedFromUnmergedJets = true;
      
      if (jet->pt<=30 ) continue;
      //if (ReducedTree->Jets_isLooseJetId[i]==false) continue;
      
      //CLEANING FROM FAT JET
      if (WWTree->nGoodAK8jets > 0) {
        if (deltaR(WWTree->ungroomed_AK8jet_eta, WWTree->ungroomed_AK8jet_phi,
                   jet->eta,jet->phi) < 0.8 )
          isCleanedFromFatJet = false;
      } 
      //else if (nGoodAK4jets>0) {
      //  if (deltaR(WWTree->AK4_jet1_eta, WWTree->AK4_jet1_phi,
      //             jet->eta,jet->phi) < 0.4 )
      //    isCleanedFromUnmergedJets = false;
      //} else if (nGoodAK4jets>1) {
      //  if (deltaR(WWTree->AK4_jet2_eta, WWTree->AK4_jet2_phi,
      //             jet->eta,jet->phi) < 0.4 )
      //    isCleanedFromUnmergedJets = false;
      //}      
      
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   jet->eta,   jet->phi) < 0.3) {
          isCleaned = false;
        }
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   jet->eta,   jet->phi) < 0.3) {
          isCleaned = false;
        }
      }
      
      if (isCleaned==false) continue;
      
      
      if (isCleanedFromUnmergedJets==true && fabs(jet->eta)<2.4)
      {
        WWTree->njets_unmerged++;
        if (jet->csv>0.605) WWTree->nBTagJet_loose_unmerged++;
        if (jet->csv>0.890) WWTree->nBTagJet_medium_unmerged++;
        if (jet->csv>0.970) WWTree->nBTagJet_tight_unmerged++;
      }
      
      if (isCleanedFromFatJet==false) continue;
      
      indexGoodVBFJets.push_back(i); //save index of the "good" vbf jets candidates
      
      if (fabs(jet->eta)>=2.4) continue;
      
      WWTree->njets++;
      AK4.SetPtEtaPhiM(jet->pt,jet->eta,jet->phi,jet->mass);
      
      
      //fill B-Tag info
      if (jet->csv>0.605) {
        WWTree->nBTagJet_loose++;
      }
      
      if (jet->csv>0.890) {  
        WWTree->nBTagJet_medium++;
      }
      
      if (jet->csv>0.970) {
        WWTree->nBTagJet_tight++;
      }
      
      
      //------------------------------
      // !!! VBF non-Puppi missing !!!
      //------------------------------
    }
    #endif
    
    
    
    /////////VBF and b-tag section Puppi
    WWTree->njetsPuppi=0;
    WWTree->nBTagJetPuppi_loose=0;
    WWTree->nBTagJetPuppi_medium=0;
    WWTree->nBTagJetPuppi_tight=0;
    
    WWTree->njetsPuppi_unmerged=0;
    WWTree->nBTagJetPuppi_loose_unmerged=0;
    WWTree->nBTagJetPuppi_medium_unmerged=0;
    WWTree->nBTagJetPuppi_tight_unmerged=0;

    oldDeltaR = 1000.;
    oldDeltaRLep = 1000.;
    indexCloserJet = -1;
    indexCloserJetLep = -1;
    deltaRbtag_prev=100.;
    deltaRbtag_prev_loose=100.;
    
    std::vector<int> indexGoodVBFJetsPuppi;
    jetArrPuppi->Clear();
    jetBrPuppi->GetEntry(jentry);
    for ( int i=0; i<jetArrPuppi->GetEntries(); i++) //loop on PuppiAK4 jet
    {
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArrPuppi)[i]);
      bool isCleaned = true;
      bool isCleanedFromFatJet = true;
      bool isCleanedFromUnmergedJets = true;
      
      if (jet->pt<=30 || jet->pt<=20) continue;
      if (!passJetLooseSel(jet)) continue;
      
      //CLEANING FROM FAT JET
      if (nGoodPuppiAK8jets > 0) {
        if (deltaR(WWTree->ungroomed_PuppiAK8_jet_eta, WWTree->ungroomed_PuppiAK8_jet_phi,
                   jet->eta,jet->phi) < 0.8 )
          isCleanedFromFatJet = false;
      } 
      //else if (nGoodPuppiAK4jets>0) {
      //  if (deltaR(WWTree->PuppiAK4_jet1_eta, WWTree->PuppiAK4_jet1_phi,
      //             jet->eta,jet->phi) < 0.4 )
      //    isCleanedFromUnmergedJets = false;
      //} else if (nGoodPuppiAK4jets>1) {
      //  if (deltaR(WWTree->PuppiAK4_jet2_eta, WWTree->PuppiAK4_jet2_phi,
      //             jet->eta,jet->phi) < 0.4 )
      //    isCleanedFromUnmergedJets = false;
      //}      
      
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   jet->eta,   jet->phi) < 0.3) {
          isCleaned = false;
        }
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   jet->eta,   jet->phi) < 0.3) {
          isCleaned = false;
        }
      }
      
      if (isCleaned==false) continue;
      
      
      if (isCleanedFromUnmergedJets==true && fabs(jet->eta)<2.4)
      {
        WWTree->njetsPuppi_unmerged++;
        if (jet->csv>0.605) WWTree->nBTagJetPuppi_loose_unmerged++;
        if (jet->csv>0.890) WWTree->nBTagJetPuppi_medium_unmerged++;
        if (jet->csv>0.970) WWTree->nBTagJetPuppi_tight_unmerged++;
      }
      
      if (isCleanedFromFatJet==false) continue;
      
      indexGoodVBFJetsPuppi.push_back(i); //save index of the "good" vbf jets candidates
      WWTree->njetsPuppi++;
      
      if (fabs(jet->eta)>=2.4) continue;
      
      AK4.SetPtEtaPhiM(jet->pt,jet->eta,jet->phi,jet->mass);
      
      
      //fill B-Tag info
      if (jet->csv>0.605) { 
        WWTree->nBTagJetPuppi_loose++;
        float deltaRbtag = JET_PuppiAK8.DeltaR(AK4);
        if (deltaRbtag>0.8 && deltaRbtag<deltaRbtag_prev_loose) {
          WWTree->deltaR_AK8_closestBtagJet_loose = deltaRbtag;
          deltaRbtag_prev_loose = deltaRbtag;
        }	  
      }
      
      if (jet->csv>0.890) {  
        WWTree->nBTagJetPuppi_medium++;
        float deltaRbtag = JET_PuppiAK8.DeltaR(AK4);
        if (deltaRbtag>0.8 && deltaRbtag<deltaRbtag_prev) {
          WWTree->deltaR_AK8_closestBtagJet = deltaRbtag;
          deltaRbtag_prev = deltaRbtag;
        }	  
      }
      
      if (jet->csv>0.970) {
        WWTree->nBTagJetPuppi_tight++;
      }
      
      float deltaRlep = W.DeltaR(AK4);
      if (deltaRlep<oldDeltaRLep) indexCloserJetLep = i;
      
      float deltaR = JET_PuppiAK8.DeltaR(AK4);
      if (deltaR<0.8) continue; //the vbf jets must be outside the had W cone
      
      if (WWTree->njets!=0) {
        if (WWTree->jet2_pt!=0) {
          WWTree->jet3_pt=jet->pt;
          WWTree->jet3_eta=jet->eta;
          WWTree->jet3_phi=jet->phi;
          WWTree->jet3_e=AK4.E();
          WWTree->jet3_btag=jet->csv;
        }
        else {
          WWTree->jet2_pt=jet->pt;
          WWTree->jet2_eta=jet->eta;
          WWTree->jet2_phi=jet->phi;
          WWTree->jet2_e=AK4.E();
          WWTree->jet2_btag=jet->csv;
        }
      }	
      
      if (deltaR<oldDeltaR)  indexCloserJet = i; //index of the closest jet to the AK8
    }
    
    
    if (indexCloserJet>=0) { //fill hadronic top mass
      const baconhep::TJet *jet_temp = (baconhep::TJet*)((*jetArrPuppi)[indexCloserJet]);
      AK4.SetPtEtaPhiM(jet_temp->pt,jet_temp->eta,jet_temp->phi,jet_temp->mass);
      WWTree->mass_ungroomedjet_closerjet  = (JET_PuppiAK8 + AK4).M();
      WWTree->AK8_closerjet_pt = AK4.Pt();
      WWTree->AK8_closerjet_eta = AK4.Eta();
      WWTree->AK8_closerjet_phi = AK4.Phi();
      WWTree->AK8_closerjet_e = AK4.E();
    }
    if (indexCloserJetLep>=0) { //fill leptonic top mass
      const baconhep::TJet *jet_lep = (baconhep::TJet*)((*jetArrPuppi)[indexCloserJetLep]);
      AK4.SetPtEtaPhiM(jet_lep->pt,jet_lep->eta,jet_lep->phi,jet_lep->mass);	// Fix energy
      WWTree->mass_leptonic_closerjet  = (W + AK4).M();
    }
    
    
    int OnlyTwoVBFTypeJets = 0;
    if (indexGoodVBFJetsPuppi.size()>=2) 
    {
      cutEff[9]++;
      float tempPtMax=0.;
      float DRvbf;
      int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
      
      for (std::size_t i=0; i<indexGoodVBFJetsPuppi.size()-1; i++) {
        for ( std::size_t ii=i+1; ii<indexGoodVBFJetsPuppi.size(); ii++) {
	  const baconhep::TJet *jet1 = (baconhep::TJet*)((*jetArrPuppi)[indexGoodVBFJetsPuppi.at(i)]);
	  const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArrPuppi)[indexGoodVBFJetsPuppi.at(ii)]);
          VBF1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
          VBF2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
          TOT = VBF1 + VBF2;
	  if (TOT.M()<500) continue;
	  if ( VBFSel==1)
	  {
		if (TOT.Pt() < tempPtMax) continue;
		tempPtMax = TOT.Pt(); //take the jet pair with largest Pt
		//cout<<i<<"\t"<<ii<<"\t tempPtMax = "<<tempPtMax<<endl;
	  }
	  else if ( VBFSel==2)
	  {
		if (TOT.M() < tempPtMax) continue;
		tempPtMax = TOT.M(); //take the jet pair with largest Pt
		//cout<<"tempPtMax = "<<tempPtMax<<endl;
	  }
	  else if ( VBFSel==3)
	  {
	  	//DRvbf = abs(deltaR(VBF1.Eta(), VBF1.Phi(), VBF2.Eta(), VBF2.Phi()));
		DRvbf = abs(VBF1.Eta()-VBF2.Eta());
		if (DRvbf < tempPtMax) continue;
		tempPtMax = DRvbf; //take the jet pair with largest Pt
		//cout<<"tempPtMax = "<<tempPtMax<<endl;
	  }
	  else
	  {
	  	cout<<"\n\nERROR:	Enter valid vbf selection criteria....\n\n"<<endl;
		exit(0);
	  }
          nVBF1 = indexGoodVBFJetsPuppi.at(i); //save position of the 1st vbf jet
          nVBF2 = indexGoodVBFJetsPuppi.at(ii); //save position of the 2nd vbf jet
        }
      }
      if (nVBF1!=-1 && nVBF2 !=-1) OnlyTwoVBFTypeJets=1;
      
      if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
      {
        // nVBF1=0; nVBF2=1;
	//cout<<"nVBF1 = "<<nVBF1<<"\tnVBF2 = "<<nVBF2<<endl;
        const baconhep::TJet *jet1 = (baconhep::TJet*)((*jetArrPuppi)[nVBF1]);
	const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArrPuppi)[nVBF2]);
        VBF1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
        VBF2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
        TOT = VBF1 + VBF2;
	

//	if (WWTree->AK4_DR_GENRECO_11 < WWTree->AK4_DR_GENRECO_12)
//	{
        	WWTree->vbf_maxpt_j1_pt = jet1->pt;
        	WWTree->vbf_maxpt_j1_eta = jet1->eta;
        	WWTree->vbf_maxpt_j1_phi = jet1->phi;
        	WWTree->vbf_maxpt_j1_e = VBF1.E();
        	WWTree->vbf_maxpt_j1_bDiscriminatorCSV = jet1->csv;
        	WWTree->vbf_maxpt_j2_pt = jet2->pt;
        	WWTree->vbf_maxpt_j2_eta = jet2->eta;
        	WWTree->vbf_maxpt_j2_phi = jet2->phi;
        	WWTree->vbf_maxpt_j2_e = VBF2.E();
        	WWTree->vbf_maxpt_j2_bDiscriminatorCSV = jet2->csv;
/*	} else
	{
        	WWTree->vbf_maxpt_j2_pt = jet1->pt;
        	WWTree->vbf_maxpt_j2_eta = jet1->eta;
        	WWTree->vbf_maxpt_j2_phi = jet1->phi;
        	WWTree->vbf_maxpt_j2_e = VBF1.E();
        	WWTree->vbf_maxpt_j2_bDiscriminatorCSV = jet1->csv;
        	WWTree->vbf_maxpt_j1_pt = jet2->pt;
        	WWTree->vbf_maxpt_j1_eta = jet2->eta;
        	WWTree->vbf_maxpt_j1_phi = jet2->phi;
        	WWTree->vbf_maxpt_j1_e = VBF2.E();
        	WWTree->vbf_maxpt_j1_bDiscriminatorCSV = jet2->csv;
	}*/
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
    }

    if (OnlyTwoVBFTypeJets == 0) continue;
        cutEff[10]++;
    
    WWTree->totalEventWeight = WWTree->genWeight*WWTree->pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight*WWTree->id_eff_Weight;
   
    
    WWTree->nEvents = TotalNumberOfEvents;
    WWTree->nNegEvents = nNegEvents;

    double a_costheta1, a_costheta2, a_costhetastar, a_Phi, a_Phi1;

    computeAngles( LEP + NU0_puppi + JET_PuppiAK8, LEP + NU0_puppi, LEP, NU0_puppi, JET_PuppiAK8,  SJ1, SJ2, a_costheta1, a_costheta2, a_Phi, a_costhetastar, a_Phi1 );
    WWTree->costheta1_type0 = (float) a_costheta1;                
    WWTree->costheta2_type0 = (float) a_costheta2;
    WWTree->costhetastar_type0 = (float) a_costhetastar;
    WWTree->phi_type0 = (float) a_Phi;
    WWTree->phi1_type0 = (float) a_Phi1;

    computeAngles( LEP + NU2_puppi + JET_PuppiAK8, LEP + NU2_puppi, LEP, NU2_puppi, JET_PuppiAK8,  SJ1, SJ2, a_costheta1, a_costheta2, a_Phi, a_costhetastar, a_Phi1);
    WWTree->costheta1_type2 = (float) a_costheta1;                
    WWTree->costheta2_type2 = (float) a_costheta2;
    WWTree->costhetastar_type2 = (float) a_costhetastar;
    WWTree->phi_type2 = (float) a_Phi;
    WWTree->phi1_type2 = (float) a_Phi1;

    computeAngles( LEP + NU1_puppi + JET_PuppiAK8, LEP + NU1_puppi, LEP, NU1_puppi, JET_PuppiAK8,  SJ1, SJ2, a_costheta1, a_costheta2, a_Phi, a_costhetastar, a_Phi1);
    WWTree->costheta1_run2 = (float) a_costheta1;                
    WWTree->costheta2_run2 = (float) a_costheta2;                
    WWTree->costhetastar_run2 = (float) a_costhetastar;
    WWTree->phi_run2 = (float) a_Phi;
    WWTree->phi1_run2 = (float) a_Phi1;

    if (fabs(VBF1.Eta() - VBF2.Eta()) == 0.0)
    {
    	 WWTree->VBSCentrality_type0 = -999.0;
    	 WWTree->VBSCentrality_type2 = -999.0;
    	 WWTree->VBSCentrality_run2 = -999.0;
    }
    else
    {
    	WWTree->VBSCentrality_type0 = (fabs(VBF1.Eta()- (((LEP + NU0_puppi).Rapidity()+JET_PuppiAK8.Rapidity()))- VBF2.Eta() ))/fabs(VBF1.Eta() - VBF2.Eta());
    	WWTree->VBSCentrality_type2 = (fabs(VBF1.Eta()- (((LEP + NU2_puppi).Rapidity()+JET_PuppiAK8.Rapidity()))- VBF2.Eta() ))/fabs(VBF1.Eta() - VBF2.Eta());
    	WWTree->VBSCentrality_run2 = (fabs(VBF1.Eta()- (((LEP + NU1_puppi).Rapidity()+JET_PuppiAK8.Rapidity()))- VBF2.Eta() ))/fabs(VBF1.Eta() - VBF2.Eta());
    }
     WWTree->RpT_type0 = (JET_PuppiAK8.Pt()*(LEP + NU0_puppi).Pt())/(VBF1.Pt()*VBF2.Pt());
     WWTree->RpT_type2 = (JET_PuppiAK8.Pt()*(LEP + NU2_puppi).Pt())/(VBF1.Pt()*VBF2.Pt());
     WWTree->RpT_run2 = (JET_PuppiAK8.Pt()*(LEP + NU1_puppi).Pt())/(VBF1.Pt()*VBF2.Pt());
     WWTree->ZeppenfeldWH = JET_PuppiAK8.Rapidity() - (VBF1.Eta() + VBF2.Eta())/2.0;
     WWTree->ZeppenfeldWL_type0 = (LEP + NU0_puppi).Rapidity() - (VBF1.Eta() + VBF2.Eta())/2.0;
     WWTree->ZeppenfeldWL_type2 = (LEP + NU2_puppi).Rapidity() - (VBF1.Eta() + VBF2.Eta())/2.0;
     WWTree->ZeppenfeldWL_run2 = (LEP + NU1_puppi).Rapidity() - (VBF1.Eta() + VBF2.Eta())/2.0;
     WWTree->LeptonProjection_type0 = (LEP.Pt()*cos(LEP.Theta()-(LEP + NU0_puppi).Theta()))/(LEP + NU0_puppi).Pt();
     WWTree->LeptonProjection_type2 = (LEP.Pt()*cos(LEP.Theta()-(LEP + NU2_puppi).Theta()))/(LEP + NU2_puppi).Pt();
     WWTree->LeptonProjection_run2 = (LEP.Pt()*cos(LEP.Theta()-(LEP + NU1_puppi).Theta()))/(LEP + NU1_puppi).Pt();

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
  std::cout << "---------end loop on events------------" << std::endl;
  std::cout << std::endl;
  std::cout << "GEN events = " << count_genEvents << std::endl;


  
  std::cout << "----------------------" << std::endl;
  std::cout << " SUMMARY" << std::endl;
  std::cout << "----------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"MC matching: "<<(float)ok/(float)total<<std::endl;
  std::cout<<"negative events: "<<nNegEvents<<std::endl;
  std::cout << std::endl;
  std::cout<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<std::endl
  	   <<"(1) Gen Events:        "<<cutEff[1]<<"\t:\t"<<((float)cutEff[1]*100.0)/(float)cutEff[0]<<std::endl
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
  delete genPartArr; delete muonArr; delete electronArr; delete vertexArr;
  delete vjetArr; delete vjetAddArr; delete jetArr; delete vjetArrPuppi; delete vjetAddArrPuppi; delete jetArrPuppi; delete lheWgtArr;
  //TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  //outROOT->cd();
  //outTree->Write();
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
