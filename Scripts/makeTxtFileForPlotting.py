#!/usr/bin/python

import sys
import datetime
import os

from ConfigParser import RawConfigParser
config = RawConfigParser()   
config.optionxform = str       # Last two lines are done because ConfigParser will not preserve case
config.read("DataMCInfo.ini")

crossSections = dict([sample, float(xsec)] for sample, xsec in config.items('CrossSection'))
Name   = dict([sample, str(nPro)] for sample, nPro in config.items('name'))
Stack  = dict([sample, int(stit)] for sample, stit in config.items('StackIt'))
colorCode = dict([sample, int(color)] for sample, color in config.items('ColorCode'))

from pprint import pprint
#print "cross sections:" 
#pprint(crossSections)
#print "name of samples:"
#pprint(Name)

source = "/eos/uscms/store/user/rasharma/SecondStep/WWTree_IsoVarSaved_2018_02_26_13h34/"

os.system('xrdfs root://cmseos.fnal.gov/ mkdir '+source+'HaddedFiles/')

ifhaddOnly = 0

Arrayfilepath = []
for root, dirs, filenames in os.walk(source):
	for f in filenames:
		#filepath = root + os.sep + f
		filepath = f
		if filepath.endswith(".root"):
			Arrayfilepath.append(filepath)

#pprint(Arrayfilepath)
#pprint("====================================================")

import difflib
#print difflib.get_close_matches('WWTree_DYJetsToLL_M-50_amcatnlo_', Arrayfilepath, 100, 0.6)
#print difflib.get_close_matches('WWTree_SingleElectron', Arrayfilepath, 100,0.55)

import os
import ROOT as ROOT

ROOT.gROOT.SetBatch(True)

List = []
ListnEvents = []
ListnNegEvents = []
for sample in crossSections.keys():
	temp = []
	#temp.append("hadd "+sample+".root")
	temp.append(sample+".root")
	print "====================> ",sample
	nEvents=0
	nNegEvents=0
	for i,files in enumerate(Arrayfilepath):
		if files.find(sample) != -1:
			#temp.append(files)
			print files
			file = ROOT.TFile(source+"/"+files)
			if not file:
				print '\n==>Failed to open %s' % Arrayfilepath[i]
				print '\n'
				exit(0)
			tree = file.Get('otree')
			if tree:
				temp.append(files)
				if files.find('Single') != -1:
					nEvents=1
					nNegEvents=0
				else:
					if ifhaddOnly != 1:
						h1 = ROOT.TH1F("h1","", 999999999, 0 , 999999999)
						h2 = ROOT.TH1F("h2","", 999999999, 0 , 999999999)
						tree.Draw("nEvents>>h1","","",10)
						tree.Draw("nNegEvents>>h2","","",10)
						nEvents+=h1.GetMean()
						nNegEvents+=h2.GetMean()
						h1.Delete()
						h2.Delete()
					else:
						print "skip..."
				print files,nEvents,nNegEvents
			file.Delete()

	if len(temp)>1:
		List.append(temp)
		ListnEvents.append(nEvents)
		ListnNegEvents.append(nNegEvents)

pprint(List)

sampleInfo = []
for samples in List:
	temp=""
	temp="hadd -f "
	for i in range(0,len(samples)):
		if i == 0:
			temp+=source+'HaddedFiles/'+samples[i]+' '
		else:
			temp+=source+"/"+samples[i]+" "
	print temp
	os.system(temp)
	#for i in range(1,len(samples)):
	#	os.system('mv '+samples[i]+' '+source)
	sampleInfo.append(samples[0])
	print ""

if ifhaddOnly != 1:
	print "=============	MAKE SUMMARY	================"	
	#sampleInfo.sort()
	pprint(sampleInfo)
	
	OutPutFile = "DibosonBoostedElMuSamples13TeV_"+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M')+".txt";
	outScript = open(OutPutFile,"w");
	outScript.write("# name           file_location  xspb/lumipb  otherscale nMCevents       nMCNegEvents    colorcode       stackit\n")
	for i,files in enumerate(sampleInfo):
		for sample in crossSections.keys():
			#print i,files
			#print files,sample
			if files.find(sample) != -1:
				#print "===> ",files,sample
				print Name[sample],"\t",files,"\t",crossSections[sample],"\t1\t",int(ListnEvents[i]),"\t",int(ListnNegEvents[i]),"\t",colorCode[sample],"\t",Stack[sample]
				outScript.write(Name[sample]+"\t"+files+"\t"+str(crossSections[sample])+"\t1\t"+str(ListnEvents[i])+"\t"+str(ListnNegEvents[i])+"\t"+str(colorCode[sample])+"\t"+str(Stack[sample])+"\n")
