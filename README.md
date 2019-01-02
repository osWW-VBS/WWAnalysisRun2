#####Table Of Content

* [Instructions](#instructions)
* [To Do List](#to-do-list)
* [POG Recipes for Moriond 2017](#pog-recipes-for-moriond-2017)
* [How To Make PU Distribution for data](#how-to-make-pu-distribution-for-data)
* [General command](#general-command)
* [Command to generate aQGC parametres summary from reweight cards](#command-to-generate-aQGC-parametres-summary-from-reweight-cards)

---
The package contains a code to produce ntuple for WW semileptonic final state.
It takes in input ntuples produced from miniAOD with the Bacon (https://github.com/ksung25/BaconProd )


## Instructions

	cmsrel CMSSW_8_0_26_patch1
	cd CMSSW_8_0_26_patch1/src
	cmsenv
	git clone git@github.com:ksung25/BaconAna.git
	cd BaconAna
	git checkout 33ffe39
        cd  ../
	mkdir WWAnalysis
	cd WWAnalysis
	git clone https://github.com/ram1123/WWAnalysisRun2.git;
	cd WWAnalysisRun2
	git checkout bacon_80x
	cd ../../
	scramv1 b -j 8
	python python/produceWWNtuples.py -i /store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/ -n WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8 -o WWTree_WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8 -w 0.9114 -no 1991227 -noNeg 0 -lumi 35900.0 --ismc 1 -trig 1 -c lpc -loc 1 

* To submit the batch job (**LXPLUS**):
	* Go to directory:

			CMSSW_8_0_26_patch1/src/WWAnalysis
			python WWAnalysisRun2/python/submit_on_lxbatch_MINIAODv2_MC.py
			python WWAnalysisRun2/python/submit_on_lxbatch_MINIAODv2_DataEle2.py
			python WWAnalysisRun2/python/submit_on_lxbatch_MINIAODv2_Data.py

* To submit the condor job (**LPC FNAL**):

		cd {...}/CMSSW_8_0_26_patch1/src/WWAnalysis/WWAnalysisRun2
		python python/submit_on_lpcCondor_MINIAODv2.py

This will give you two files named `runstep2condor.jdl` and `runstep2condor.sh`. To submit the condor job do

	voms-proxy-init --voms cms --valid 168:00  # if proxy was not set
	condor_submit runstep2condor.jdl

Monitor the status of jobs using:

		condor_q -submitter rasharma

* To see the various options availabe with **produceWWNtuples.py** do,

		python python/produceWWNtuples.py --help

## To Do List
- [ ] Clean the code

## POG Recipes for Moriond 2017

* Pile-up reweighting xsec = 69.2mb

Ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/POGRecipesICHEP2016

## How To Make PU Distribution for data

	pileupCalc.py -i MyAnalysisJSON.txt --inputLumiJSON pileup_latest.txt  --calcMode true --minBiasXsec 69200 --maxPileupBin 50 --numPileupBins 50  MyDataPileupHistogram.root
	pileupCalc.py -i Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --inputLumiJSON pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 70 --numPileupBins 70 MyDataPileupHistogram.root

where,
* MyAnalysisJSON.txt is the JSON file we are using.
* pileup_latest.txt : this is input json file. It can be found at link: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt

Reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II



## General command

	grep "time to run this code =" *.stdout | awk '{print $10,$7/60}'

	voms-proxy-init

	grep -r --exclude=\*.{root,o,exe,swp,bcup} genGravMass *

## Search data/mc on DAS

	dasgoclient --query="dataset=/*/RunIISpring16MiniAOD*/MINIAODSIM" --limit=0

There is a script to check many samples at once. Script name is `DasGoClientSummary.py`.

### How to use script DasGoClientSummary.py

1. Check if it runs fine:

		grep "Job Nubmer" <LogFileName>
   
   Lets say you are checking this for N samples. Then it should have all numbers from 1 to N.

2. Grab dataset name using this command:

		grep Dataset screenlog.0 | awk -F "/" '{print $2,"\t",$3}'

3. Grab number of events using command:

		grep nevents screenlog.0 | awk -F "," '{print $3}' | awk -F ":" '{print $2}'

4. To investigate we can paste output of step 2 and 3 in spreadsheet and look at it.
## Command to generate aQGC parametres summary from reweight cards

	grep launch aQGC_WMhadZlepJJ_EWK_LO_NPle1_mjj100pt10_reweight_card.dat | awk -F "=" '{print $2}' | awk -F "_" '{ gsub("p",".",$2); gsub("m","-",$2); print $1,$2}'

	grep launch aQGC_WMhadZlepJJ_EWK_LO_NPle1_mjj100pt10_reweight_card.dat | awk -F "=" '{print $2}' | awk -F "_" '{ gsub("p",".",$2); gsub("m","-",$2); print $1}'

	grep launch aQGC_WMhadZlepJJ_EWK_LO_NPle1_mjj100pt10_reweight_card.dat | awk -F "=" '{print $2}' | awk -F "_" '{ gsub("p",".",$2); gsub("m","-",$2); print $2}'
