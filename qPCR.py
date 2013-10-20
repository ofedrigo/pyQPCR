#! /usr/bin/env python

import math,getopt,sys,sets,operator

LOG_FILENAME = 'qPCR.log'
logout=open(LOG_FILENAME,'w')

#/* FUNCTIONS */

def product(thisList): #prodcut
	result=1.0
	for item in thisList:
		result=result*float(item)
	return result

def geometric(thisList): #geometric mean
	return math.pow(product(thisList),1.0/float(len(thisList)))

def deltaCt(theseSamples,eff,theseControlGenes,whichRef): # calculate delta CT
	#calculate Quantities
	controlQ,sampleQ,minCT,allCTs={},{},{},[]
	# collect all Ct mean data and output the minium mean Ct per gene
	for geneName in theseSamples.keys(): # goes through all meantCt and stdvCt for the genes
		vector=[]
		for tissueName in theseSamples[geneName].keys(): # goes through tissuees
			for sampleID in theseSamples[geneName][tissueName].keys(): #goes through sample names
				vector.append(theseSamples[geneName][tissueName][sampleID]["mean"])
				allCTs.append(theseSamples[geneName][tissueName][sampleID]["mean"])
		minCT[geneName]=min(vector)
	for geneName in theseSamples.keys(): #goes though genes
		for tissueName in theseSamples[geneName].keys(): #goes through tissues
			for sampleID in theseSamples[geneName][tissueName].keys(): # goes through individuals
				#different ways to calculate the Ct mean reference
				if whichRef==0: # mean of all Ct means
					reference=mean(allCTs)
				elif whichRef==1: # minimum of all C tmeans
					reference=min(allCTs)
				deltaCT=reference-theseSamples[geneName][tissueName][sampleID]["mean"] # calculate delta Ct
				Q=math.pow(eff[geneName],deltaCT) # calcualte relative quantities, accounting for efficiency
				SDQ=Q*math.log(eff[geneName])*theseSamples[geneName][tissueName][sampleID]["stdv"] # propagate stdv
				# we need to seperate data if it is a gene of interest or a control gene
				if geneName in theseControlGenes: # control genes
					if (controlQ.has_key(tissueName))==False: controlQ[tissueName]={}
					if (controlQ[tissueName].has_key(sampleID))==False: controlQ[tissueName][sampleID]={}
					controlQ[tissueName][sampleID][geneName]={} #this assumes that in a given plate you should have only 1 gene/smapleId/tissue otherwise these are technical replicates and it has been already averaged
					controlQ[tissueName][sampleID][geneName]["mean"]=Q
					controlQ[tissueName][sampleID][geneName]["stdv"]=SDQ
				else: # gene of interest
					if (sampleQ.has_key(tissueName))==False: sampleQ[tissueName]={}
					if (sampleQ[tissueName].has_key(sampleID))==False: sampleQ[tissueName][sampleID]={}
					sampleQ[tissueName][sampleID][geneName]={} #this assumes that in a given plate you should have only 1 gene/smapleId/tissue otherwise these are technical replicates and it has been already averaged
					sampleQ[tissueName][sampleID][geneName]["mean"]=Q
					sampleQ[tissueName][sampleID][geneName]["stdv"]=SDQ
	return sampleQ,controlQ

def isFloat(f): #check if this is the string can be changed in float
	newf="".join(f.split("."))
	return newf.isdigit()
	
def checkSamples(samples,listOfGOI,listOfControls,interrunsFlag,interruns,interrunsGenes,interrunsControls,fullNames): #this is cleaning up the data when normalizer are missing, we need to eliminate the sample set associated with it
	newSample={}
	interrunData={}
	allTissues=sets.Set([])
	allSpecies=sets.Set([])
	for geneName in listOfGOI:
		for tissueName in samples[geneName].keys():
			for sampleID in samples[geneName][tissueName].keys():
				if tissueName!=interruns[0] and sampleID!=interruns[1]:
					allTissues.add(tissueName[1:])
					allSpecies.add(fullNames[tissueName[0]])
					cc=[]
					for controlGene in listOfControls:
						if samples.has_key(controlGene):
							if samples[controlGene].has_key(tissueName):
								if samples[controlGene][tissueName].has_key(sampleID): cc.append(controlGene)
					if len(cc)<len(listOfControls):
						line="Sample set ["+geneName+", "+tissueName+", "+sampleID+"] is missing control gene(s) "+",".join(filter(lambda x:x not in cc,listOfControls)+filter(lambda x:x not in listOfControls,cc))
						logout.write(line)
						#print line
					else:
						if (newSample.has_key(geneName))==False: newSample[geneName]={} # gene name
						if (newSample[geneName].has_key(tissueName))==False: newSample[geneName][tissueName]={} # tissue name
						if (newSample[geneName][tissueName].has_key(sampleID))==False: newSample[geneName][tissueName][sampleID]={} # individual name
						newSample[geneName][tissueName][sampleID]["mean"]=samples[geneName][tissueName][sampleID]["mean"]
						newSample[geneName][tissueName][sampleID]["stdv"]=samples[geneName][tissueName][sampleID]["stdv"]
						for controlGene in listOfControls:
							if (newSample.has_key(controlGene))==False: newSample[controlGene]={} # gene name
							if (newSample[controlGene].has_key(tissueName))==False: newSample[controlGene][tissueName]={} # tissue name
							if (newSample[controlGene][tissueName].has_key(sampleID))==False: newSample[controlGene][tissueName][sampleID]={} # individual name
							newSample[controlGene][tissueName][sampleID]["mean"]=samples[controlGene][tissueName][sampleID]["mean"]
							newSample[controlGene][tissueName][sampleID]["stdv"]=samples[controlGene][tissueName][sampleID]["stdv"]
	print "\t# genes of interest: "+",".join(listOfGOI)
	print "\t# control genes: "+",".join(listOfControls)
	print "\t# inter-plates genes of interest: "+",".join(interrunsGenes)
	print "\t# inter-plates control genes: "+",".join(interrunsControls)
	print "\t# tissues: "+",".join(list(allTissues))
	print "\t# species: "+",".join(list(allSpecies))
	for geneName in interrunsGenes:
		tissueName=interruns[0]
		sampleID=interruns[1]
		if (interrunData.has_key(geneName))==False: interrunData[geneName]={} # gene name
		if (interrunData[geneName].has_key(tissueName))==False: interrunData[geneName][tissueName]={} # tissue name
		if (interrunData[geneName][tissueName].has_key(sampleID))==False: interrunData[geneName][tissueName][sampleID]={} # individual name
		interrunData[geneName][tissueName][sampleID]["mean"]=samples[geneName][tissueName][sampleID]["mean"]
		interrunData[geneName][tissueName][sampleID]["stdv"]=samples[geneName][tissueName][sampleID]["stdv"]
		for controlGene in listOfControls:
			if (interrunData.has_key(controlGene))==False: interrunData[controlGene]={} # gene name
			if (interrunData[controlGene].has_key(tissueName))==False: interrunData[controlGene][tissueName]={} # tissue name
			if (interrunData[controlGene][tissueName].has_key(sampleID))==False: interrunData[controlGene][tissueName][sampleID]={} # individual name
			interrunData[controlGene][tissueName][sampleID]["mean"]=samples[controlGene][tissueName][sampleID]["mean"]
			interrunData[controlGene][tissueName][sampleID]["stdv"]=samples[controlGene][tissueName][sampleID]["stdv"]
	exitFlag=checkInterPlates(interrunData,interruns,interrunsGenes,interrunsControls)
	return newSample,interrunData,exitFlag

def checkInterPlates(samples,interruns,interrunsGenes,interrunsControl):
	exitFlag=False
	tissueName=interruns[0]
	sampleID=interruns[1]
	totalGenes=len(interrunsGenes)+len(interrunsControl)
	cc=0
	for geneName in (interrunsGenes+interrunsControl):
		if samples[geneName].has_key(tissueName):
			if samples[geneName][tissueName].has_key(sampleID):
				cc=cc+1
	if cc!=totalGenes: exitFlag=True
	return exitFlag

def getData(fileName): #get data and return data and effciency
	fileHandle=open(fileName,'r')
	myData=fileHandle.read().splitlines()
	fileHandle.close()
	thisData=[]
	for line in myData: thisData.append(line.strip(" "))
	amplEff_vector={} #get primer efficiencies
	index=thisData.index("gene	individual name	species/tissue name	CT")
	for datum in thisData[0:index]: # get efficiency and gene names
		items=datum.split("\t")
		geneName=items[0]
		efficiency=(float(items[1])/100.0)+1.0 #assume that efficiency will be between 0 and 100
		amplEff_vector[geneName]=efficiency
	technicals_data={}
	for datum in thisData[index+1:]:
		if datum[0]!="#": #ignore lines starting with "#"
			geneName,sampleID,tissueName,CTvalue=datum.split("\t")
			if isFloat(CTvalue): #deprecated: if CTvalue!="Undet.":
				if (technicals_data.has_key(tissueName))==False: technicals_data[tissueName]={} # tissues
				if (technicals_data[tissueName].has_key(geneName))==False: technicals_data[tissueName][geneName]={} # genename
				technicals_data[tissueName][geneName].setdefault(sampleID,[]).append(float(CTvalue)) # tehcnical replicates per individual
	return technicals_data,amplEff_vector

def same(string1,string2):
	if string1==string2 or string1.upper()==string2.upper():
		return True
	return False

def getParameters(parameterFile):
	fileNames,controls,GOIs,thresh,threshFlag,interrunsF,interrunsI,interrunsG,exprRef,outputFileName,fullNames,replicates,interrunsC=[],[],[],-1,False,False,[],[],0,"qPCR.out",{},False,[]
	fileHandle=open(parameterFile,'r')
	inpData=fileHandle.read().splitlines()
	fileHandle.close()
	cc=sets.Set([])
	for line in inpData:
		if line[0]!="#":
			subject,items=[x.strip(" ") for x in line.split("=")]
			if same(subject,"PLATES"):
				if items!="": fileNames=items.split(",")
				cc.add(subject)
			if same(subject,"CONTROLS"):
				if items!="": controls=items.split(",")
				cc.add(subject)
			if same(subject,"GOIs"):
				if items!="": GOIs=items.split(",")
				cc.add(subject)
			if same(subject,"STDV"):
				if items!="-1":
					thresh=float(items)
					threshFlag=True
				cc.add(subject)
			if same(subject,"INTERRUNS (TISSUE,INDIVIDUAL)"):
				if items!="": interrunsI=items.split(",")
				cc.add(subject)
			if same(subject,"INTERRUNS (GENES)"):
				if items!="": interrunsG=items.split(",")
				cc.add(subject)
			if same(subject,"INTERRUNS (CONTROLS)"):
				if items!="": interrunsC=items.split(",")
				cc.add(subject)
			if same(subject,"INTERRUNS"):
				interrunsF=same(items,"Yes")
				cc.add(subject)
			if same(subject,"EXPRESSION REFERENCE"):
				ExprRef=int(items)
				cc.add(subject)
			if same(subject,"OUTPUT FILE NAME"):
				outputFileName=items
				cc.add(subject)
			if same(subject,"FULL SPECIES NAMES"):
				if items!="": entries=items.split(",")
				for entry in entries:
					abreviation,fullName=entry.split(":")
					fullNames[abreviation]=fullName
				cc.add(subject)
			if same(subject,"REPLICATES"):
				replicates=same(items,"Yes")
				cc.add(subject)
	if len(list(cc))<12:
		print "Parameter file incomplete or contaning errors"
		sys.exit()
	return fileNames,controls,GOIs,thresh,threshFlag,interrunsF,interrunsI,interrunsG,ExprRef,outputFileName,fullNames,replicates,interrunsC

def ss(inlist): # sum of squares
	ss = 0
	for item in inlist:
		ss = ss + item*item
	return ss

def variance (inlist): # variance
	n = len(inlist)
	mn = mean(inlist)
	deviations = [0]*len(inlist)
	for i in range(len(inlist)):
		deviations[i] = inlist[i] - mn
	return ss(deviations)/float(n-1)

def stdev (inlist): # standard deviation
	if len(inlist)>1:
		return math.sqrt(variance(inlist))
	else:
		return 0

def mean(thisList): #calculate the mean
	return float(sum(thisList))/float(len(thisList))

def meanCt(thistechnicals,thresh,flag): #calculate mean CT
	samples_vector={}
	for tissueName in thistechnicals.keys(): # loop through tissuee names
		for geneName in thistechnicals[tissueName].keys():	# loop through gene name
			for sampleID in thistechnicals[tissueName][geneName].keys():	# loop through individal names
				meanCT=mean(thistechnicals[tissueName][geneName][sampleID]) # get the mean of 1 individual for all technical replicates
				stdvCT=stdev(thistechnicals[tissueName][geneName][sampleID])# get the stdv of 1 individual for all technical replicates
				if stdvCT>thresh and flag==True:
					line="Sample ["+tissueName+", "+geneName+", "+sampleID+"] eliminated because stdv > "+str(stdvCT)
					logout.write(line)
					#print line
				elif flag==False or stdvCT<=thresh: #to modifiy if we want to incllude everything
					# now we store mean and stdv per gene/tissue/individual
					if (samples_vector.has_key(geneName))==False: samples_vector[geneName]={} # gene name
					if (samples_vector[geneName].has_key(tissueName))==False: samples_vector[geneName][tissueName]={} # tissue name
					if (samples_vector[geneName][tissueName].has_key(sampleID))==False: samples_vector[geneName][tissueName][sampleID]={} # individual name
					samples_vector[geneName][tissueName][sampleID]["mean"]=meanCT
					samples_vector[geneName][tissueName][sampleID]["stdv"]=stdvCT
	return samples_vector

def getNF(Qvalues): #calculate normalizing factors from
	NFresults={}
	for sampleID in Qvalues.keys():
		vector=[] #create a vector of meanCt per individual/gene
		for geneName in Qvalues[sampleID].keys(): vector.append(Qvalues[sampleID][geneName]["mean"])
		thisNF=geometric(vector)
		SDNF=thisNF*math.pow(sum([math.pow(Qvalues[sampleID][geneName]["stdv"]/float(len(Qvalues[sampleID].keys())*Qvalues[sampleID][geneName]["mean"]),2) for geneName in Qvalues[sampleID].keys()]),0.5)
		NFresults[sampleID]={}
		NFresults[sampleID]["NF"]=thisNF
		NFresults[sampleID]["SDNF"]=SDNF
	return NFresults

def getAllNFs(controlQ):
	NFs={} # calculate normalizing factors per tissue/species
	for tissueName in controlQ.keys(): NFs[tissueName]=getNF(controlQ[tissueName])
	return NFs

def normalize(sampleQ,NFs): #normalize quantities
	#normalize quantities
	GOIs={}
	for tissueName in sampleQ.keys(): #goes through all tissues
		for sampleID in sampleQ[tissueName].keys(): # goes though all individuals
			for geneName in sampleQ[tissueName][sampleID].keys(): # goes through all genes
				normalized=sampleQ[tissueName][sampleID][geneName]["mean"]/NFs[tissueName][sampleID]["NF"]
				stdvNormalized=normalized*math.pow(math.pow(NFs[tissueName][sampleID]["SDNF"]/NFs[tissueName][sampleID]["NF"],2)+math.pow(sampleQ[tissueName][sampleID][geneName]["stdv"]/sampleQ[tissueName][sampleID][geneName]["mean"],2),0.5)
				if (GOIs.has_key(geneName))==False: GOIs[geneName]={} #we create a Gene Of Interest hash with gene/tissue/individual and store normalized mean and normalized stdv
				if (GOIs[geneName].has_key(tissueName))==False: GOIs[geneName][tissueName]={}
				GOIs[geneName][tissueName][sampleID]={}
				GOIs[geneName][tissueName][sampleID]["mean"]=normalized
				GOIs[geneName][tissueName][sampleID]["stdv"]=stdvNormalized
	return GOIs

def formatPrint(inData,title,legend):
	line="\n\n"+title+"\n"+legend
	for keyOne in inData.keys():
		for keyTwo in inData[keyOne].keys():
			for keyThree in inData[keyOne][keyTwo].keys():
				line=line+"\n"+keyOne+"\t"+keyTwo+"\t"+keyThree+"\t"+str(inData[keyOne][keyTwo][keyThree]["mean"])+"\t"+str(inData[keyOne][keyTwo][keyThree]["stdv"])
	return line

def interPlates(GOIs,interrunsFlag,interruns,interrunsGenes,interrunsControls,amplEff,exprRef,interrunData):
	# use inter plates normalizers
	interrun={}
	if interrunsFlag==True:
		print "\tInter-plates calibration is performed"
		interrunSampleQ,interrunControlQ=deltaCt(interrunData,amplEff,interrunsControls,exprRef)
		if len(interrunsControls)>0:
			print "\t   Used Hellemans et al., 2007 method to normalize (geometric mean of normalized relative quantities)"
			interrunNFs=getAllNFs(interrunControlQ) # calculate normalizing factors per tissue/species
			interrunGOIs=normalize(interrunSampleQ,interrunNFs) #normalize quantities
			allInterruns=[]
			for thisGenecontrol in interrunsGenes: allInterruns.append(interrunData[thisGenecontrol][interruns[0]][interruns[1]]["mean"])
			CF=geometric(allInterruns)
		else:
			print "\t   Normalize with the geometric mean of relative quantities" 
			allInterruns=[]
			for thisGenecontrol in interrunsGenes: allInterruns.append(interrunSampleQ[interruns[0]][interruns[1]][thisGenecontrol]["mean"])
			CF=geometric(allInterruns)
		for geneName in GOIs.keys():
			interrun[geneName]={}
			for tissueName in GOIs[geneName].keys():
					interrun[geneName][tissueName]={}
					for sampleID in GOIs[geneName][tissueName].keys():
						scaled=GOIs[geneName][tissueName][sampleID]["mean"]/CF
						interrun[geneName][tissueName][sampleID]={}
						interrun[geneName][tissueName][sampleID]["mean"]=scaled
						interrun[geneName][tissueName][sampleID]["stdv"]=GOIs[geneName][tissueName][sampleID]["stdv"]/CF
						finalVector.append(scaled)
	else:
		print "\nNo itner-plates calibration performed"
		for geneName in GOIs.keys():
			interrun[geneName]={}
			for tissueName in GOIs[geneName].keys():
					interrun[geneName][tissueName]={}
					for sampleID in GOIs[geneName][tissueName].keys():
						interrun[geneName][tissueName][sampleID]={}
						interrun[geneName][tissueName][sampleID]["mean"]=GOIs[geneName][tissueName][sampleID]["mean"]
						interrun[geneName][tissueName][sampleID]["stdv"]=GOIs[geneName][tissueName][sampleID]["stdv"]
						finalVector.append(GOIs[geneName][tissueName][sampleID]["mean"])
	return interrun,finalVector

fileNames_updated,interrun,finalVector=[],{},[]
parameterFile="aParameters.in"
plateNames,listOfControls,listOfGOI,thresh,threshFlag,interrunsFlag,interruns,interrunsGenes,exprRef,outputFileName,fullNames,replicates,interrunsControls=getParameters(parameterFile)
outputhandle=open(outputFileName,'w')
print""
for fileName in plateNames:
	print "\n"+fileName
	outputhandle.write("\n\n*** "+fileName+" ***")
	technicals,amplEff=getData(fileName) #get technical replicates[sampleName][[geneName][sampleId] and amplEfficiency for each gene
	samples=meanCt(technicals,thresh,threshFlag) # calculate meanCt across technical replicates and returns samples[geneName][sampleName][sampleId] mean and stdv
	outputhandle.write(formatPrint(samples,"Mean CT","Gene\tTissue\tIndividuals\tMean\tStdv"))
	#make new listOfControls, listOfGOI
	tempGOIs,tempControls=[],[]
	for geneName in samples.keys():
		if geneName in listOfGOI:
			tempGOIs.append(geneName)
		elif geneName in listOfControls:
			tempControls.append(geneName)
	listOfGOI=tempGOIs
	listOfControls=tempControls
	samples,interrunData,exitFlag=checkSamples(samples,listOfGOI,listOfControls,interrunsFlag,interruns,interrunsGenes,interrunsControls,fullNames) # check if the data is complete
	if exitFlag==True and interrunsFlag==True:
		line="Experiment "+fileName+" can't be analyzed with interplates normalizers"
		logout.write(line)
		print line
	else:
		outputhandle.write("\n"+fileName+" passed")
		fileNames_updated.append(fileName)
		sampleQ,controlQ=deltaCt(samples,amplEff,listOfControls,exprRef) #get quantities
		NFs=getAllNFs(controlQ) # calculate normalizing factors per tissue/species
		GOIs=normalize(sampleQ,NFs) #normalize quantities
		outputhandle.write(formatPrint(GOIs,"Normalized","Gene\tTissue\tIndividuals\tMean\tStdv"))
		interrun[fileName],finalVector=interPlates(GOIs,interrunsFlag,interruns,interrunsGenes,interrunsControls,amplEff,exprRef,interrunData) # uses inter plates normalizers
		outputhandle.write(formatPrint(interrun[fileName],"Interrun normalized","Gene\tTissue\tIndividuals\tMean\tStdv")) #output final results

# COMBINE ALL THE NORMALIZED DATA
if len(finalVector)==0:
	line="\nError: all the plates have been ignored because of missing inter-plates normalizer values\n"
	logout.write(line)
	print line
	sys.exit()
else:
	minExpress=min(finalVector) #normalized by the smallest scaled expression => lowest expressed==1
	globalTable={}
	for fileName in fileNames_updated: # consider only filenames that passed
		outputhandle.write("\n\n\t"+fileName+"\n")
		for geneName in interrun[fileName].keys():
			for tissueName in interrun[fileName][geneName].keys():
				for sampleID in interrun[fileName][geneName][tissueName].keys():
					if sampleID!=interruns[1]:
						speciesName=tissueName[0]
						experimentName=geneName+"_"+tissueName
						if (globalTable.has_key(experimentName))==False: globalTable[experimentName]={}
						if (globalTable[experimentName].has_key(speciesName))==False: globalTable[experimentName][speciesName]={}
						if (globalTable[experimentName][speciesName].has_key(sampleID))==False: globalTable[experimentName][speciesName][sampleID]={}
						globalTable[experimentName][speciesName][sampleID][fileName]={}
						globalTable[experimentName][speciesName][sampleID][fileName]["mean"]=interrun[fileName][geneName][tissueName][sampleID]["mean"]/minExpress
						globalTable[experimentName][speciesName][sampleID][fileName]["stdv"]=interrun[fileName][geneName][tissueName][sampleID]["stdv"]/minExpress
						outputhandle.write("\t\t"+fullNames[speciesName]+"\t"+geneName+"_"+tissueName[1:]+"\t"+sampleID+"\t"+str(interrun[fileName][geneName][tissueName][sampleID]["mean"]/minExpress)+"\t"+str(interrun[fileName][geneName][tissueName][sampleID]["stdv"]/minExpress)+"\n")
	
	newR={}
	if replicates==True: # if there are more replicates, average them
		print "\nAverage all replicates across plates"
		for experimentName in globalTable.keys():
			listTemp=sets.Set([])
			for speciesName in globalTable[experimentName].keys():
				for sampleID in globalTable[experimentName][speciesName].keys():
					for fileName in globalTable[experimentName][speciesName][sampleID].keys():
						listTemp.add(fileName)
			listOfPlates=list(listTemp)
			for speciesName in globalTable[experimentName].keys():
				for sampleID in globalTable[experimentName][speciesName].keys():
					line="\n"+sampleID
					val=[]
					val2=[]
					for fileName in listOfPlates:
						if globalTable[experimentName][speciesName][sampleID].has_key(fileName):
							line=line+"\t"+str(globalTable[experimentName][speciesName][sampleID][fileName])
							val.append(globalTable[experimentName][speciesName][sampleID][fileName]["mean"])
							val2.append(globalTable[experimentName][speciesName][sampleID][fileName]["stdv"])
					geneName,tissueName=experimentName.split("_")
					if newR.has_key(geneName)==False: newR[geneName]={}
					if newR[geneName].has_key(tissueName[1:])==False: newR[geneName][tissueName[1:]]={}
					if newR[geneName][tissueName[1:]].has_key(speciesName)==False: newR[geneName][tissueName[1:]][speciesName]={}
					newR[geneName][tissueName[1:]][speciesName][sampleID]=[str(mean(val))]
					"""
					newVal=[]
					for i in range(len(val2)):
						if val2[i]/val[i]<=0.3:
							newVal.append(val[i])
					if len(newVal)>0:
						newR[geneName][tissueName[1:]][speciesName][sampleID]=[str(mean(newVal))]
					else:
						newR[geneName][tissueName[1:]][speciesName][sampleID]=["None"]
					"""
	else:
		print "\nOutput all replicates across plates"
		for experimentName in globalTable.keys():
			listTemp=sets.Set([])	
			for speciesName in globalTable[experimentName].keys():
				for sampleID in globalTable[experimentName][speciesName].keys():
					for fileName in globalTable[experimentName][speciesName][sampleID].keys():
						listTemp.add(fileName)
			listOfPlates=list(listTemp)
			for speciesName in globalTable[experimentName].keys():
				for sampleID in globalTable[experimentName][speciesName].keys():
					line="\n"+sampleID
					geneName,tissueName=experimentName.split("_")
					if newR.has_key(geneName)==False: newR[geneName]={}
					if newR[geneName].has_key(tissueName[1:])==False: newR[geneName][tissueName[1:]]={}
					if newR[geneName][tissueName[1:]].has_key(speciesName)==False: newR[geneName][tissueName[1:]][speciesName]={}
					for fileName in listOfPlates:
						if globalTable[experimentName][speciesName][sampleID].has_key(fileName):
							line=line+"\t"+str(globalTable[experimentName][speciesName][sampleID][fileName])
							newR[geneName][tissueName[1:]][speciesName].setdefault(sampleID,[]).append(str(globalTable[experimentName][speciesName][sampleID][fileName]["mean"])+"("+str(globalTable[experimentName][speciesName][sampleID][fileName]["stdv"])+", "+str(globalTable[experimentName][speciesName][sampleID][fileName]["stdv"]/globalTable[experimentName][speciesName][sampleID][fileName]["mean"])+")")
	outputhandle.write("\n\nSpecies,Gene,Tissue,Individual,Expression")
	for geneName in newR.keys():
			outputhandle.write("\n\n*** "+geneName+" ***")
			for tissueName in newR[geneName].keys():
				outputhandle.write("\n\n"+tissueName)
				for speciesName in newR[geneName][tissueName].keys():
					for sampleID in newR[geneName][tissueName][speciesName].keys():
						outputhandle.write("\n"+"\t".join([fullNames[speciesName],geneName,tissueName,sampleID,"\t".join([str(x) for x in newR[geneName][tissueName][speciesName][sampleID]])]))
	outputhandle.close()
print "\nsee results in "+outputFileName
print "and possible errors in "+LOG_FILENAME+"\n"
logout.close()