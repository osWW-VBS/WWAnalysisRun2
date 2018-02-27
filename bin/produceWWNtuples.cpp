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
  TClonesArray *genPartArr 	= new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr    	= new TClonesArray("baconhep::TMuon");
  TClonesArray *electronArr	= new TClonesArray("baconhep::TElectron");
  TClonesArray *vertexArr	= new TClonesArray("baconhep::TVertex");
  TClonesArray *vjetArr		= new TClonesArray("baconhep::TJet");
  TClonesArray *vjetAddArr  	= new TClonesArray("baconhep::TAddJet");
  TClonesArray *jetArr		= new TClonesArray("baconhep::TJet");
  TClonesArray *vjetArrPuppi	= new TClonesArray("baconhep::TJet");
  TClonesArray *vjetAddArrPuppi	= new TClonesArray("baconhep::TAddJet");
  TClonesArray *jetArrPuppi	= new TClonesArray("baconhep::TJet");
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
  

  vector<const baconhep::TJet*> goodJetsv;
  
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


  JetCorrectorParameters paramAK4puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFPuppi.txt");
  JetCorrectorParameters paramAK4chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
  JetCorrectorParameters paramAK8chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");
  JetCorrectorParameters paramAK8puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFPuppi.txt");
  JetCorrectionUncertainty *fJetUnc_AK4chs = new JetCorrectionUncertainty(paramAK4chs);
  JetCorrectionUncertainty *fJetUnc_AK4puppi = new JetCorrectionUncertainty(paramAK4puppi);
  JetCorrectionUncertainty *fJetUnc_AK8chs = new JetCorrectionUncertainty(paramAK8chs);
  JetCorrectionUncertainty *fJetUnc_AK8puppi = new JetCorrectionUncertainty(paramAK8puppi);

  //---------start loop on events------------
  std::cout << "---------start loop on events------------" << std::endl;
  jentry2=0;
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
      genPartArr->Clear();
      genBr->GetEntry(jentry);
      genPartBr->GetEntry(jentry);
      TLorentzVector hadW, lepW, VBFJ1, VBFJ2, VBFJ, temp;
      TLorentzVector genLep, genNeutrino, genWquarks, genVBFquarks;
      std::vector<TLorentzVector> v_GEN_hadW, v_GEN_lepW, v_GEN_VBFJ1, v_GEN_VBFJ2, v_GEN_VBFJ, v_GEN_temp;
      std::vector<TLorentzVector> v_genLep, v_genNeutrino, v_genWquarks, v_genVBFquarks;

      v_GEN_hadW.clear();	v_GEN_lepW.clear();	v_GEN_VBFJ1.clear();	v_GEN_VBFJ2.clear();
      v_GEN_VBFJ.clear();	v_GEN_temp.clear();	v_genLep.clear();	v_genNeutrino.clear();
      v_genWquarks.clear();	v_genVBFquarks.clear();
      
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
	  WWTree->nu_pz_gen	  = v_genNeutrino[0].Pz();
	  WWTree->nu_pt_gen	  = v_genNeutrino[0].Pt();
	  WWTree->nu_eta_gen	  = v_genNeutrino[0].Eta();
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
		//WWTree->LHEid[i] = lhe->id;
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



    ///////////THE FAT JET - AK8
    float tempTTbarMass=0.;
    float tempMassW = 3000.0;
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
	if (jet->pt<200 || fabs(jet->eta)>2.4)  continue;
        //if (!passJetLooseSel(jet)) continue;
	if (addjet->mass_prun>tempTTbarMass) {
	  if ( (jet->eta>0 && WWTree->l_eta1<0) || 
	      (jet->eta<0 && WWTree->l_eta1>0)) { //jet and lepton in opposite hemisphere for ttb
	    //ttb_jet_position=i; //save AK8 jet in ttbar topology
	    tempTTbarMass=addjet->mass_prun;
	  }
	}
	if (abs(addjet->mass_sd0 - 80.385) > tempMassW ) continue; //save the jet which is closest to W-mass
	//if (abs(jet->mass - 80.385) > tempMassW ) continue; //save the jet which is closest to W-mass
	
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
	WWTree->ungroomed_AK8jet_charge = jet->q;
	//WWTree->ungroomed_jet_pt_jes_up = (jet->pt[i]/ReducedTree->AK8Jets_AK8correction[i])*ReducedTree->AK8Jets_AK8correctionUp[i];
	//WWTree->ungroomed_jet_pt_jes_dn = (ReducedTree->AK8CHS_pt[i]/ReducedTree->AK8Jets_AK8correction[i])*ReducedTree->AK8Jets_AK8correctionDown[i];
	
	WWTree->AK8jet_mass     = jet->mass;
	WWTree->AK8jet_mass_pr  = addjet->mass_prun;
	WWTree->AK8jet_mass_so  = addjet->mass_sd0;
	WWTree->AK8jet_mass_tr  = addjet->mass_trim;
	WWTree->AK8jet_tau2tau1 = addjet->tau2/addjet->tau1;
	WWTree->AK8jet_sj1_pt	= addjet->sj1_pt;
	WWTree->AK8jet_sj1_eta	= addjet->sj1_eta;
	WWTree->AK8jet_sj1_phi	= addjet->sj1_phi;
	WWTree->AK8jet_sj1_m	= addjet->sj1_m;
	WWTree->AK8jet_sj1_q	= addjet->sj1_q;
	WWTree->AK8jet_sj2_pt	= addjet->sj2_pt;
	WWTree->AK8jet_sj2_eta	= addjet->sj2_eta;
	WWTree->AK8jet_sj2_phi	= addjet->sj2_phi;
	WWTree->AK8jet_sj2_m	= addjet->sj2_m;
	WWTree->AK8jet_sj2_q	= addjet->sj2_q;
	WWTree->AK8_jetID_loose = passJetLooseSel(jet);

	WWTree->AK8jet_e3_b1	= addjet->e3_b1;
	WWTree->AK8jet_e3_v1_b1	= addjet->e3_v1_b1;
	WWTree->AK8jet_e3_v2_b1	= addjet->e3_v2_b1;
	WWTree->AK8jet_e4_v1_b1	= addjet->e4_v1_b1;
	WWTree->AK8jet_e4_v2_b1	= addjet->e4_v2_b1;
	WWTree->AK8jet_e3_b2	= addjet->e3_b2;
	WWTree->AK8jet_e3_v1_b2	= addjet->e3_v1_b2;
	WWTree->AK8jet_e3_v2_b2	= addjet->e3_v2_b2;
	WWTree->AK8jet_e4_v1_b2	= addjet->e4_v1_b2;
	WWTree->AK8jet_e4_v2_b2	= addjet->e4_v2_b2;
	
	WWTree->AK8jet_e2_sdb1	= addjet->e2_sdb1;
	WWTree->AK8jet_e3_sdb1	= addjet->e3_sdb1;
	WWTree->AK8jet_e3_v1_sdb1	= addjet->e3_v1_sdb1;
	WWTree->AK8jet_e3_v2_sdb1	= addjet->e3_v2_sdb1;
	WWTree->AK8jet_e4_v1_sdb1	= addjet->e4_v1_sdb1;
	WWTree->AK8jet_e4_v2_sdb1	= addjet->e4_v2_sdb1;    
	
	WWTree->AK8jet_e2_sdb2	= addjet->e2_sdb2;
	WWTree->AK8jet_e3_sdb2	= addjet->e3_sdb2;
	WWTree->AK8jet_e3_v1_sdb2	= addjet->e3_v1_sdb2;
	WWTree->AK8jet_e3_v2_sdb2	= addjet->e3_v2_sdb2;
	WWTree->AK8jet_e4_v1_sdb2	= addjet->e4_v1_sdb2;
	WWTree->AK8jet_e4_v2_sdb2	= addjet->e4_v2_sdb2;    // Soft Dropped correlation function in puts beta=2
	
	//WWTree->AK8jet_e2_sdb4	= addjet->e2_sdb4;
	//WWTree->AK8jet_e3_sdb4	= addjet->e3_sdb4;
	//WWTree->AK8jet_e3_v1_sdb4	= addjet->e3_v1_sdb4;
	//WWTree->AK8jet_e3_v2_sdb4	= addjet->e3_v2_sdb4;
	//WWTree->AK8jet_e4_v1_sdb4	= addjet->e4_v1_sdb4;
	//WWTree->AK8jet_e4_v2_sdb4	= addjet->e4_v2_sdb4;        // Soft Dropped correlation function in puts beta=4
	//
	//WWTree->AK8jet_e2_sdb05	= addjet->e2_sdb05;
	//WWTree->AK8jet_e3_sdb05	= addjet->e3_sdb05;
	//WWTree->AK8jet_e3_v1_sdb05	= addjet->e3_v1_sdb05;
	//WWTree->AK8jet_e3_v2_sdb05	= addjet->e3_v2_sdb05;
	//WWTree->AK8jet_e4_v1_sdb05	= addjet->e4_v1_sdb05;
	//WWTree->AK8jet_e4_v2_sdb05	= addjet->e4_v2_sdb05;  // Soft Dropped correlation function in puts beta=0.5

	WWTree->AK8jet_qjet	= addjet->qjet;

	//WWTree->jet_mass_pr_jes_up = (ReducedTree->AddAK8CHS_mass_prun[i]/ReducedTree->AK8Jets_AK8massCorrection[i])*ReducedTree->AK8Jets_AK8massCorrectionUp[i];
	//WWTree->jet_mass_pr_jes_dn = (ReducedTree->AddAK8CHS_mass_prun[i]/ReducedTree->AK8Jets_AK8massCorrection[i])*ReducedTree->AK8Jets_AK8massCorrectionDown[i];
      
	tempMassW = abs(addjet->mass_sd0 - 80.385);
	//tempMassW = abs(jet->mass - 80.385);
	WWTree->nGoodAK8jets++;
      }
    if(WWTree->ungroomed_AK8jet_pt > 0)
      {
	JET.SetPtEtaPhiE(WWTree->ungroomed_AK8jet_pt,WWTree->ungroomed_AK8jet_eta,WWTree->ungroomed_AK8jet_phi,WWTree->ungroomed_AK8jet_e);
	SJ1.SetPtEtaPhiM(WWTree->AK8jet_sj1_pt, WWTree->AK8jet_sj1_eta, WWTree->AK8jet_sj1_phi, WWTree->AK8jet_sj1_m);
	SJ2.SetPtEtaPhiM(WWTree->AK8jet_sj2_pt, WWTree->AK8jet_sj2_eta, WWTree->AK8jet_sj2_phi, WWTree->AK8jet_sj2_m);
	// calculate Up variation
        double unc = func( JET.Pt(), JET.Eta(), fJetUnc_AK8chs); 
	JET_jes_up.SetPtEtaPhiM((1.+unc)*JET.Pt(), JET.Eta(), JET.Phi(), (1.+unc)*JET.M());	
	// calculate Down variation
	JET_jes_dn.SetPtEtaPhiM((1.-unc)*JET.Pt(), JET.Eta(), JET.Phi(), (1.-unc)*JET.M());	
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
	if (jet->pt<200 || fabs(jet->eta)>2.4)  continue; //be careful: this is not inside the synchntuple code
        //if (!passJetLooseSel(jet)) continue;
	if (addjet->mass_prun>tempTTbarMass) {
	  if ( (jet->eta>0 && WWTree->l_eta1<0) || 
	      (jet->eta<0 && WWTree->l_eta1>0)) { //jet and lepton in opposite hemisphere for ttb
	    ttb_PuppiAK8_jet_position=i; //save AK8 jet in ttbar topology
	    tempTTbarMass=addjet->mass_prun;
	  }
	}
	//if (abs(jet->mass - 80.385) > tempMassW) continue; //save the jet closest to w-mass
	if (abs(addjet->mass_sd0 - 80.385) > tempMassW) continue; //save the jet closest to w-mass
	
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
  
  //WWTree->PuppiAK8jet_e2_sdb4  = addjet->e2_sdb4 ;
  //WWTree->PuppiAK8jet_e3_sdb4  = addjet->e3_sdb4 ;
  //WWTree->PuppiAK8jet_e3_v1_sdb4 = addjet->e3_v1_sdb4  ;
  //WWTree->PuppiAK8jet_e3_v2_sdb4 = addjet->e3_v2_sdb4  ;
  //WWTree->PuppiAK8jet_e4_v1_sdb4 = addjet->e4_v1_sdb4  ;
  //WWTree->PuppiAK8jet_e4_v2_sdb4 = addjet->e4_v2_sdb4;        // Soft Dropped correlation function in puts beta=4
  //
  //WWTree->PuppiAK8jet_e2_sdb05 = addjet->e2_sdb05  ;
  //WWTree->PuppiAK8jet_e3_sdb05 = addjet->e3_sdb05  ;
  //WWTree->PuppiAK8jet_e3_v1_sdb05  = addjet->e3_v1_sdb05 ;
  //WWTree->PuppiAK8jet_e3_v2_sdb05  = addjet->e3_v2_sdb05 ;
  //WWTree->PuppiAK8jet_e4_v1_sdb05  = addjet->e4_v1_sdb05 ;
  //WWTree->PuppiAK8jet_e4_v2_sdb05  = addjet->e4_v2_sdb05;  // Soft Dropped correlation function in puts beta=0.5

  WWTree->PuppiAK8jet_qjet = addjet->qjet;
	//     WWTree->PuppiAK8_jet_mass_pr_jes_up = (addjet->mass_prun[i]/ReducedTree->PuppiAK8Jets_PuppiAK8massCorrection[i])*ReducedTree->PuppiAK8Jets_PuppiAK8massCorrectionUp[i];
	//   WWTree->PuppiAK8_jet_mass_pr_jes_dn = (addjet->mass_prun[i]/ReducedTree->PuppiAK8Jets_PuppiAK8massCorrection[i])*ReducedTree->PuppiAK8Jets_PuppiAK8massCorrectionDown[i];
	
	//tempPt = WWTree->ungroomed_PuppiAK8_jet_pt;
	//tempMassW = abs(jet->mass - 80.385);
	tempMassW = abs(WWTree->PuppiAK8_jet_mass_so - 80.385);
	nGoodPuppiAK8jets++;
	//hadWPuppiAK8pos = i;
      }
    if (WWTree->ungroomed_PuppiAK8_jet_pt > 0.)
      {
	//JET_PuppiAK8.SetPtEtaPhiE(WWTree->ungroomed_PuppiAK8_jet_pt,WWTree->ungroomed_PuppiAK8_jet_eta,WWTree->ungroomed_PuppiAK8_jet_phi,WWTree->ungroomed_PuppiAK8_jet_e);
	JET_PuppiAK8.SetPtEtaPhiM(WWTree->ungroomed_PuppiAK8_jet_pt,WWTree->ungroomed_PuppiAK8_jet_eta,WWTree->ungroomed_PuppiAK8_jet_phi,WWTree->PuppiAK8_jet_mass_so_corr);
	SJ1_PuppiAK8.SetPtEtaPhiM(WWTree->PuppiAK8_jet_sj1_pt, WWTree->PuppiAK8_jet_sj1_eta, WWTree->PuppiAK8_jet_sj1_phi, WWTree->PuppiAK8_jet_sj1_m);
	SJ2_PuppiAK8.SetPtEtaPhiM(WWTree->PuppiAK8_jet_sj2_pt, WWTree->PuppiAK8_jet_sj2_eta, WWTree->PuppiAK8_jet_sj2_phi, WWTree->PuppiAK8_jet_sj2_m);
	// calculate Up variation
        double unc = func(JET_PuppiAK8.Pt(), JET_PuppiAK8.Eta(), fJetUnc_AK8puppi); 
	//JET_PuppiAK8_jes_up.SetPtEtaPhiM((1.+unc)*JET_PuppiAK8.Pt(), JET_PuppiAK8.Eta(), JET_PuppiAK8.Phi(), (1.+unc)*JET_PuppiAK8.M());	
	JET_PuppiAK8_jes_up.SetPtEtaPhiM((1.+unc)*JET_PuppiAK8.Pt(), JET_PuppiAK8.Eta(), JET_PuppiAK8.Phi(), (1.+unc)*WWTree->PuppiAK8_jet_mass_so_corr);	
	// calculate Down variation
	//JET_PuppiAK8_jes_dn.SetPtEtaPhiM((1.-unc)*JET_PuppiAK8.Pt(), JET_PuppiAK8.Eta(), JET_PuppiAK8.Phi(), (1.-unc)*JET_PuppiAK8.M());	
	JET_PuppiAK8_jes_dn.SetPtEtaPhiM((1.-unc)*JET_PuppiAK8.Pt(), JET_PuppiAK8.Eta(), JET_PuppiAK8.Phi(), (1.-unc)*WWTree->PuppiAK8_jet_mass_so_corr);	
      }
    
    
    // FAT JET SELECTION
    bool isGoodFatJet = true;
    if (WWTree->nGoodAK8jets==0 && nGoodPuppiAK8jets==0) isGoodFatJet = false; //not found a good hadronic W candidate
    if (WWTree->ungroomed_AK8jet_pt<200 && WWTree->ungroomed_PuppiAK8_jet_pt<200) isGoodFatJet = false;
    if (!isGoodFatJet) continue;
    cutEff[5]++;
    
    
    WWTree->ungroomed_AK8jet_pt_jes_up = JET_jes_up.Pt();
    WWTree->ungroomed_AK8jet_eta_jes_up = JET_jes_up.Eta();
    WWTree->ungroomed_AK8jet_phi_jes_up = JET_jes_up.Phi();
    WWTree->ungroomed_AK8jet_mass_jes_up = JET_jes_up.M();	
    WWTree->ungroomed_AK8jet_pt_jes_dn = JET_jes_dn.Pt();
    WWTree->ungroomed_AK8jet_eta_jes_dn = JET_jes_dn.Eta();
    WWTree->ungroomed_AK8jet_phi_jes_dn = JET_jes_dn.Phi();
    WWTree->ungroomed_AK8jet_mass_jes_dn = JET_jes_dn.M();	
    WWTree->ungroomed_PuppiAK8_jet_pt_jes_up = JET_PuppiAK8_jes_up.Pt();// calculated with softdrop corrected mass
    WWTree->ungroomed_PuppiAK8_jet_eta_jes_up = JET_PuppiAK8_jes_up.Eta();// calculated with softdrop corrected mass
    WWTree->ungroomed_PuppiAK8_jet_phi_jes_up = JET_PuppiAK8_jes_up.Phi();// calculated with softdrop corrected mass
    WWTree->ungroomed_PuppiAK8_jet_mass_jes_up = JET_PuppiAK8_jes_up.M();// calculated with softdrop corrected mass
    WWTree->ungroomed_PuppiAK8_jet_pt_jes_dn = JET_PuppiAK8_jes_dn.Pt();// calculated with softdrop corrected mass
    WWTree->ungroomed_PuppiAK8_jet_eta_jes_dn = JET_PuppiAK8_jes_dn.Eta();// calculated with softdrop corrected mass
    WWTree->ungroomed_PuppiAK8_jet_phi_jes_dn = JET_PuppiAK8_jes_dn.Phi();// calculated with softdrop corrected mass
    WWTree->ungroomed_PuppiAK8_jet_mass_jes_dn = JET_PuppiAK8_jes_dn.M();// calculated with softdrop corrected mass
    //////////////////ANGULAR VARIABLES
    WWTree->deltaR_Wjet_GenReco = deltaR(WWTree->hadW_eta_gen,WWTree->hadW_phi_gen,JET.Eta(),JET.Phi());
    

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
    int OnlyTwoVBFTypeJets = 0;
    
    std::vector<int> indexGoodVBFJets;


    jetArr->Clear();
    jetBr->GetEntry(jentry);
    for ( int i=0; i<jetArr->GetEntries(); i++) //loop on AK4 jet
    {
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[i]);
      bool isCleaned = true;
      bool isCleanedFromFatJet = true;

        double unc = func(jet->pt, jet->eta, fJetUnc_AK4chs); 
	// calculate Up variation
	VBF1_jes_up.SetPtEtaPhiM((1.+unc)*jet->pt, jet->eta, jet->phi, (1.+unc)*jet->mass);	
	// calculate Down variation
	VBF1_jes_dn.SetPtEtaPhiM((1.-unc)*jet->pt, jet->eta, jet->phi, (1.-unc)*jet->mass);	
      
      if (jet->pt<=20 && VBF1_jes_up.Pt()<=20 && VBF1_jes_dn.Pt()<=20 ) continue;
      if (!passJetLooseSel(jet)) continue;

      //fill B-Tag info
      
      if (fabs(jet->eta) < 2.4 && jet->pt>20){
      		if (jet->csv>0.5426)  WWTree->nBTagJet_loose++;
      		if (jet->csv>0.8484)  WWTree->nBTagJet_medium++;
      		if (jet->csv>0.9535)  WWTree->nBTagJet_tight++;

		goodJetsv.push_back(jet);
      }
      
      //CLEANING FROM FAT JET
      if (WWTree->nGoodAK8jets > 0) {
        if (deltaR(WWTree->ungroomed_AK8jet_eta, WWTree->ungroomed_AK8jet_phi,
                   jet->eta,jet->phi) < 0.8 )
          isCleanedFromFatJet = false;
      } 
      
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
      if (isCleanedFromFatJet==false) continue;
      
      indexGoodVBFJets.push_back(i); //save index of the "good" vbf jets candidates
      
      if (fabs(jet->eta)>=2.4) continue;
      
      WWTree->njets++;
      AK4.SetPtEtaPhiM(jet->pt,jet->eta,jet->phi,jet->mass);
      
      
      
      //------------------------------
      // !!! VBF non-Puppi missing !!!
      //------------------------------
    }

    int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
    if (indexGoodVBFJets.size()>=2) 
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

        double unc1 = func(jet1->pt, jet1->eta, fJetUnc_AK4chs); 
	// calculate Up variation
	VBF1_jes_up.SetPtEtaPhiM((1.+unc1)*jet1->pt, jet1->eta, jet1->phi, (1.+unc1)*jet1->mass);	
	// calculate Down variation
	VBF1_jes_dn.SetPtEtaPhiM((1.-unc1)*jet1->pt, jet1->eta, jet1->phi, (1.-unc1)*jet1->mass);	
	// calculate Up variation
        double unc2 = func(jet2->pt, jet2->eta, fJetUnc_AK4chs); 
	VBF2_jes_up.SetPtEtaPhiM((1.+unc2)*jet2->pt, jet2->eta, jet2->phi, (1.+unc2)*jet2->mass);	
	// calculate Down variation
	VBF2_jes_dn.SetPtEtaPhiM((1.-unc2)*jet2->pt, jet2->eta, jet2->phi, (1.-unc2)*jet2->mass);	
	  if (TOT.M()<500 && (VBF1_jes_up+VBF2_jes_up).M()<500 && (VBF1_jes_dn+VBF2_jes_dn).M()<500 ) continue;
	  if ( VBFSel==1)
	  {
		if (TOT.Pt() < tempPtMax) continue;
		tempPtMax = TOT.Pt(); //take the jet pair with largest Pt
	  }
	  else if ( VBFSel==2)
	  {
		if (TOT.M() < tempPtMax) continue;
		tempPtMax = TOT.M(); //take the jet pair with largest Pt
	  }
	  else if ( VBFSel==3)
	  {
		DRvbf = abs(VBF1.Eta()-VBF2.Eta());
		if (DRvbf < tempPtMax) continue;
		tempPtMax = DRvbf; //take the jet pair with largest Pt
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

	// calculate Up variation
        double unc1 = func(jet1->pt, jet1->eta, fJetUnc_AK4chs); 
	VBF1_jes_up.SetPtEtaPhiM((1.+unc1)*jet1->pt, jet1->eta, jet1->phi, (1.+unc1)*jet1->mass);	
	// calculate Down variation
	VBF1_jes_dn.SetPtEtaPhiM((1.-unc1)*jet1->pt, jet1->eta, jet1->phi, (1.-unc1)*jet1->mass);	
	// calculate Up variation
        double unc2 = func(jet2->pt, jet2->eta, fJetUnc_AK4chs); 
	VBF2_jes_up.SetPtEtaPhiM((1.+unc2)*jet2->pt, jet2->eta, jet2->phi, (1.+unc2)*jet2->mass);
	// calculate Down variation
	VBF2_jes_dn.SetPtEtaPhiM((1.-unc2)*jet2->pt, jet2->eta, jet2->phi, (1.-unc2)*jet2->mass);	

	WWTree->vbf_maxpt_j1_pt_jes_up 	= VBF1_jes_up.Pt();
	WWTree->vbf_maxpt_j1_eta_jes_up = VBF1_jes_up.Eta();
	WWTree->vbf_maxpt_j1_phi_jes_up = VBF1_jes_up.Phi();
	WWTree->vbf_maxpt_j1_mass_jes_up= VBF1_jes_up.M();
	WWTree->vbf_maxpt_j1_pt_jes_dn 	= VBF1_jes_dn.Pt();
	WWTree->vbf_maxpt_j1_eta_jes_dn = VBF1_jes_dn.Eta();
	WWTree->vbf_maxpt_j1_phi_jes_dn = VBF1_jes_dn.Phi();
	WWTree->vbf_maxpt_j1_mass_jes_dn= VBF1_jes_dn.M();
	WWTree->vbf_maxpt_j2_pt_jes_up 	= VBF2_jes_up.Pt();
	WWTree->vbf_maxpt_j2_eta_jes_up = VBF2_jes_up.Eta();
	WWTree->vbf_maxpt_j2_phi_jes_up = VBF2_jes_up.Phi();
	WWTree->vbf_maxpt_j2_mass_jes_up= VBF2_jes_up.M();
	WWTree->vbf_maxpt_j2_pt_jes_dn 	= VBF2_jes_dn.Pt();
	WWTree->vbf_maxpt_j2_eta_jes_dn = VBF2_jes_dn.Eta();
	WWTree->vbf_maxpt_j2_phi_jes_dn = VBF2_jes_dn.Phi();
	WWTree->vbf_maxpt_j2_mass_jes_dn= VBF2_jes_dn.M();

        WWTree->vbf_maxpt_jj_pt = TOT.Pt();
        WWTree->vbf_maxpt_jj_pt_jes_up = (VBF1_jes_up + VBF2_jes_up).Pt();
        WWTree->vbf_maxpt_jj_pt_jes_dn = (VBF1_jes_dn + VBF2_jes_dn).Pt();
        WWTree->vbf_maxpt_jj_eta = TOT.Eta();
        WWTree->vbf_maxpt_jj_phi = TOT.Phi();
        WWTree->vbf_maxpt_jj_m = TOT.M();	
        WWTree->vbf_maxpt_jj_m_jes_up = (VBF1_jes_up + VBF2_jes_up).M();
        WWTree->vbf_maxpt_jj_m_jes_dn = (VBF1_jes_dn + VBF2_jes_dn).M();
	WWTree->vbf_maxpt_jj_Deta = abs(VBF1.Eta() - VBF2.Eta());
	WWTree->vbf_maxpt_jj_Deta_jes_up = abs(VBF1_jes_up.Eta() - VBF2_jes_up.Eta());
	WWTree->vbf_maxpt_jj_Deta_jes_dn = abs(VBF1_jes_dn.Eta() - VBF2_jes_dn.Eta());
	WWTree->AK4_DR_GENRECO_11 = abs(deltaR(WWTree->AK4_1_eta_gen, WWTree->AK4_1_phi_gen, WWTree->vbf_maxpt_j1_eta, WWTree->vbf_maxpt_j1_phi));
	WWTree->AK4_DR_GENRECO_12 = abs(deltaR(WWTree->AK4_1_eta_gen, WWTree->AK4_1_phi_gen, WWTree->vbf_maxpt_j2_eta, WWTree->vbf_maxpt_j2_phi));
	WWTree->AK4_DR_GENRECO_21 = abs(deltaR(WWTree->AK4_2_eta_gen, WWTree->AK4_2_phi_gen, WWTree->vbf_maxpt_j1_eta, WWTree->vbf_maxpt_j1_phi));
	WWTree->AK4_DR_GENRECO_22 = abs(deltaR(WWTree->AK4_2_eta_gen, WWTree->AK4_2_phi_gen, WWTree->vbf_maxpt_j2_eta, WWTree->vbf_maxpt_j2_phi));
      }
    }
    indexGoodVBFJets.clear();
    jetArr->Clear();
    jetBr->GetEntry(jentry);
    double ptBalanceForLepSearch;
    //cout<<"DEBUG : 1 : nVBF1 = "<<nVBF1<<"\tnVBF2 = "<<nVBF2<<endl;
    for ( int i=0; i<jetArr->GetEntries(); i++) //loop on AK4 jet
    {
      if (i == nVBF1 || i == nVBF2) continue;
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[i]);
      if (!passJetLooseSel(jet)) continue;
      if (jet->pt<=30) continue;
      indexGoodVBFJets.push_back(i); //save index of the "good" vbf jets candidates
    }
    double min = 1000.0;
    for (std::size_t i=0; i<indexGoodVBFJets.size()-1; i++) {
      for ( std::size_t ii=i+1; ii<indexGoodVBFJets.size(); ii++) {
        const baconhep::TJet *jet1 = (baconhep::TJet*)((*jetArr)[indexGoodVBFJets.at(i)]);
        const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArr)[indexGoodVBFJets.at(ii)]);
        if (jet1->pt < 30) continue;
        if (jet2->pt < 50) continue;
        VBF1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
        VBF2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
        TOT = VBF1 + VBF2;
	ptBalanceForLepSearch = TOT.Pt() - ( WWTree->ungroomed_PuppiAK8_jet_pt + WWTree->vbf_maxpt_j1_pt + WWTree->vbf_maxpt_j2_pt); 
	if (min > ptBalanceForLepSearch)
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

		min = ptBalanceForLepSearch;
	}

      }
    }

    
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
    /////////////////////////////////////////////////////////////////////////////////////////
    //
    //		Calculate b-tag weight
    //		
    /////////////////////////////////////////////////////////////////////////////////////////

    	vector<float> btagSFv       = getJetSFs(goodJetsv, bTagReader, "central", "central");
    	vector<float> btagSFUpHFv   = getJetSFs(goodJetsv, bTagReader, "up",      "central");
    	vector<float> btagSFDownHFv = getJetSFs(goodJetsv, bTagReader, "down",    "central");
    	vector<float> btagSFUpLFv   = getJetSFs(goodJetsv, bTagReader, "central", "up");
    	vector<float> btagSFDownLFv = getJetSFs(goodJetsv, bTagReader, "central", "down"); 

    	WWTree->btag0Wgt       = getBtagEventReweightEtaBin(0, goodJetsv, btagSFv);
    	WWTree->btag1Wgt       = getBtagEventReweightEtaBin(1, goodJetsv, btagSFv);
    	WWTree->btag2Wgt       = getBtagEventReweightEtaBin(2, goodJetsv, btagSFv);

	//if (isnan(btagSFv) == 1) cout<<"********************* btag SF is nan"<<endl;
	if (isnan(WWTree->btag0Wgt) == 1 || WWTree->btag0Wgt == 0) { 
	cout<<"********************* btag weight is nan"<<endl;
	//for (unsigned int i=0; i<btagSFv.size(); i++)
	//cout<<"Vector = "<<btagSFv[i]<<"\t";
	//cout<<endl;
	}

    	WWTree->btag0WgtUpHF   = getBtagEventReweightEtaBin(0, goodJetsv, btagSFUpHFv);
    	WWTree->btag0WgtDownHF = getBtagEventReweightEtaBin(0, goodJetsv, btagSFDownHFv);
    	WWTree->btag0WgtUpLF   = getBtagEventReweightEtaBin(0, goodJetsv, btagSFUpLFv);
    	WWTree->btag0WgtDownLF = getBtagEventReweightEtaBin(0, goodJetsv, btagSFDownLFv);
    	WWTree->btag1WgtUpHF   = getBtagEventReweightEtaBin(1, goodJetsv, btagSFUpHFv);
    	WWTree->btag1WgtDownHF = getBtagEventReweightEtaBin(1, goodJetsv, btagSFDownHFv);
    	WWTree->btag1WgtUpLF   = getBtagEventReweightEtaBin(1, goodJetsv, btagSFUpLFv);
    	WWTree->btag1WgtDownLF = getBtagEventReweightEtaBin(1, goodJetsv, btagSFDownLFv);
    	
    	btagSFv.clear();	btagSFUpHFv.clear();	btagSFDownHFv.clear();
    	btagSFUpLFv.clear();	btagSFDownLFv.clear();
    /////////////////////////////////////////////////////////////////////////////////////////	END btag weight calculation
    
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
    
    if (isnan(WWTree->costheta1_type0) == 1)
    {
    	cout<<jentry2<< "\tcostheta1_type0 is NaN" << endl;
	if (LEP2.Beta()>1) cout<<"beta > 1"<<endl;
	printf("------------- LEP2 beta = %.17g\n", LEP2.Beta());
	std::cout << "#######    LEP2 beta = 1 + " << (LEP2.Beta() - 1.0) << std::endl;
	printf(" Neutrino mass = %.17g\n", LEP2.M());
    }
    WWTree->costheta2_type0 = (float) a_costheta2;
    if (isnan(WWTree->costheta2_type0) == 1)
    {	
    	cout<<jentry2<< "\tcostheta2_type0 is NaN" << endl;
	if (LEP2.Beta()>1) cout<<"beta > 1"<<endl;
	printf("------------- LEP2 beta = %.17g\n", LEP2.Beta());
	std::cout << "#######    LEP2 beta = 1 + " << (LEP2.Beta() - 1.0) << std::endl;
	printf(" Neutrino mass = %.17g\n", LEP2.M());
	}
    WWTree->costhetastar_type0 = (float) a_costhetastar;
    if (isnan(WWTree->costhetastar_type0) == 1)
    	{cout<< jentry2 << "\tWWTree->costhetastar_type0 is NaN" << endl;
	if (LEP2.Beta()>1) cout<<"beta > 1"<<endl;
	printf("------------- LEP2 beta = %.17g\n", LEP2.Beta());
	std::cout << "#######    LEP2 beta = 1 + " << (LEP2.Beta() - 1.0) << std::endl;
	printf(" Neutrino mass = %.17g\n", LEP2.M());
	}
    WWTree->phi_type0 = (float) a_Phi;
    if (isnan(WWTree->phi_type0) == 1)
    	{cout<<jentry2<< "\tphi_type0 is NaN" << endl;
	if (LEP2.Beta()>1) cout<<"beta > 1"<<endl;
	printf("------------- LEP2 beta = %.17g\n", LEP2.Beta());
	std::cout << "#######    LEP2 beta = 1 + " << (LEP2.Beta() - 1.0) << std::endl;
	printf(" Neutrino mass = %.17g\n", LEP2.M());
	}
    WWTree->phi1_type0 = (float) a_Phi1;
    if (isnan(WWTree->phi1_type0) == 1)
    	{cout<<jentry2<< "\tphi1_type0 is NaN" << endl;
	if (LEP2.Beta()>1) cout<<"beta > 1"<<endl;
	printf("------------- LEP2 beta = %.17g\n", LEP2.Beta());
	std::cout << "#######    LEP2 beta = 1 + " << (LEP2.Beta() - 1.0) << std::endl;
	printf(" Neutrino mass = %.17g\n", LEP2.M());
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

    goodJetsv.clear();
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
