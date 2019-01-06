#!/usr/bin/env python
from __future__ import absolute_import, division, print_function


####################################################################################-------BRIEF DESCRIPTION------###################################################################################################
#																																																					#
#Date: September 2018																																																#
#Purpose: This script is written to do some post processing on the SCRMshaw multiple offsets output and then Running peaks calling algorithm (MACs) in order to get more robust CRM predictions						#
#Input: It needs two input 																																															#
# (i) combined SCRMshaw_offset's output (including multiple training sets and multiple methods)																														#
# (ii) number of predictions to extract from each of the offset's prediction default = 5000																															#																																							
#Outputs of this script are individual BED formated PEAKs file for each of the training set and each method, which can be concatenated to one file and used as an input for the evaluation pipeline script			#
################################################-####################################################################################################################################################################
import os
import pybedtools
import statistics
import argparse
import sys
import shutil
from scipy import stats
import scipy.stats
import pprint
import csv
import subprocess
from collections import Counter
import numpy as np
import pandas as pd


##############################################------------FUNCTION---------------##########################################
# This functions will parse combined file of multiple outputs of scrmshaw to individual unique files for each of the training set and statistical method used and will return three methods dictionaries containing all the training sets as their keys
def parse_output(outfile,numparse): 
	# creating three separate dictionaries based on methods used  
	d_hexmcd={}
	d_imm={}
	d_pac={}
	x={''}
	global d_hexmcd_val
	global d_imm_val
	global d_pac_val
	d_hexmcd_val=['']
	d_imm_val=['']
	d_pac_val=['']

	with open (outfile) as file:
		rows=(line2.split('\t') for line2 in file)
		for row in rows:
		#based on the 14th column(names of different data sets) and 15th column (statistical method used) of scrmshaw_joined_output file giving values to each of the three method`s dictionaries
			if (row[15]=='hexmcd') and (int(row[16]) <= int(numparse)):
				#print(row[16])
				#print(numparse)
				if row[14] not in d_hexmcd:
					myRow = [] # create a new list to use
					myRow.append(row) # add my new row to our new list
					d_hexmcd[row[14]] = myRow  #create a new entry if it isn't in the dictionary already
				else:
					d_hexmcd.get(row[14]).append(row)
					#count_hexmcd=count_hexmcd+1
			elif (row[15]=='imm')and (int(row[16]) <= int(numparse)):
				if row[14] not in d_imm:
					myRow = []
					myRow.append(row)
					d_imm[row[14]] = myRow
				else:
					d_imm.get(row[14]).append(row)
			elif (row[15]=='pac') and (int(row[16]) <= int(numparse) ):
				if row[14] not in d_pac:
					myRow = []
					myRow.append(row)
					d_pac[row[14]] = myRow
				else:
					d_pac.get(row[14]).append(row)

			
		#calculating number of keys(datasets) each dictionary ends up having		
		for key in d_hexmcd.keys():
			d_hexmcd_val.append(key)
		for key in d_imm.keys():
			d_imm_val.append(key)
		for key in d_pac.keys():
			d_pac_val.append(key)
			
	#creating separate files for each method w.r.t datasets, using the above newly created three dictionaries and moving them to tmp (temporary folder).
	
	#These individual unique files(based on methods and their training sets) will go through downstream processes(calling peaks through MACs) one by one
	for key in d_hexmcd.keys():
		noOflines=len(d_hexmcd[key])
		
		with open(os.path.join(subdirectory,'hexmcd_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_hexmcd[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")
	
	for key in d_imm.keys():
		with open(os.path.join(subdirectory,'imm_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_imm[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")

	for key in d_pac.keys():
		with open(os.path.join(subdirectory,'pac_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_pac[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")
				
	return(d_imm_val,d_hexmcd_val,d_pac_val)


#------------------------------------------------------------------------------------------------------

# This function will extract user specified number of predictions from each of the offset for the given training set
def extract_topN_scrms(fullLengthFilePath,cutoff,method,TSET):
	#extract num of scrms from offset
	extractedFileName=str(cutoff)+'.'+method+'_'+TSET+'.bed'
	
	with open(fullLengthFilePath,'r') as infile, open(os.path.join(subdirectory,extractedFileName),'w') as outfile:
		for line in infile:

			col=line.split('\t')
			#print(col[10])
			rank=col[16].strip('\n')
			if int(rank) <= int(cutoff):
				outfile.write(line)
				
	#path=os.path.abspath(extractedFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==extractedFileName:
				path=os.path.abspath(os.path.join(root,name))
	
	return path
#---------------------------------------------------------------------------------------------------------------------------
#This function taked the tab delimited file and returns bed version of it to perform the functions of bed

def bed_conversion(tab_delimited_path):
	
	bed_tabdelimited=pybedtools.BedTool(tab_delimited_path)
		
	return bed_tabdelimited

#---------------------------------------------------------------------------------------------------------------------------
#This function will take in tab delimited file path and convert it into py bed version of that file(on which bedtools functions can be applied like sort/merge etc) and return its path
def bedtools_sorting(extractedScrmsPathBED,sortedFileName):
	
	extractedScrmsPathBED.sort().saveas(subdirectory+'/'+sortedFileName)	
	
	#path=os.path.abspath(sortedFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==sortedFileName:
				path=os.path.abspath(os.path.join(root,name))
	
	return path
	
	
	
#---------------------------------------------------------------------------------------------------------------------------
#This function will take in BED file and sum up the score for each of the 10 bp overlapping window and return every 10 bp window with its score(i.e summed)
def sum_of_score(sortedFilePath,sortedFileName,sumOfScoreFileName):
	#calculating sum of score and saving it to a dictionary
	diction={}

	with open(sortedFilePath,'r') as infile:

		for line in infile:
			col=line.split('\t')
			chrName=col[0]
			startCoord=col[1]
			endCoord=col[2]
			score=float(col[3])

			j=int(startCoord)
			tempI=int(startCoord)+10
			bp=range(tempI,int(endCoord)+10,10)

			for i in bp:
				keyName=chrName+':'+str(j)+'-'+str(i)
				#print(keyName)
				if keyName not in diction:
					#print('No')
					diction[keyName]=score

				else:
					#print('yes')
					diction[keyName]=diction[keyName]+score	
		
				j+=10	


	file2=sortedFileName.strip('bed')+'csv'

	#converting dictionary to csv file
	outfile=csv.writer(open('sumScored'+file2,'w'))
	for key, val in diction.items():
		outfile.writerow([key,val])			
	#pprint.pprint(diction)	


	#from csv to bed file
	file1='sumScored'+file2

	with open(file1,'r') as infile, open(os.path.join(subdirectory,sumOfScoreFileName),'w') as outfile: 

		for line in infile:
			col=line.split(',')
			coord=col[0]
			score=col[1]
			chrAndCoords=coord.split(':')
			chrName=chrAndCoords[0]
			Coords=chrAndCoords[1]

			bothCoords=Coords.split('-')
			startCoord= bothCoords[0]
			endCoord=bothCoords[1]

			outfile.write(chrName+'\t'+startCoord+'\t'+endCoord+'\t'+score)

	
	#path=os.path.abspath(sumOfScoreFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==sumOfScoreFileName:
				path=os.path.abspath(os.path.join(root,name))
	
	return path
#---------------------------------------------------------------------------------------------------------------------------
#This function will take in BED file and return whole genome coverage including any of the missing coordinates with the score value of 0.
def filling_missing_coords(sortedSumOfScoreFilePath,wholeGenomeFileName):
	with open(sortedSumOfScoreFilePath,'r') as infile,open(os.path.join(subdirectory,wholeGenomeFileName),'w') as outfile:
		prevEnd=0
		chrNames={}
		for line in infile:
			col=line.split("\t")

			chrName=col[0]
			start=col[1]
			end=col[2]
			score=col[3]
			if chrName not in chrNames:
				chrNames[chrName]=''
				prevEnd=0
			if prevEnd==start:
				outfile.write(line)

			elif int(start) > int(prevEnd) :
				outfile.write(chrName+'\t'+str(prevEnd)+'\t'+start+'\t'+str('0.0')+'\n')
				outfile.write(line)


			prevEnd=end
	
	#path=os.path.abspath(wholeGenomeFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==wholeGenomeFileName:
				path=os.path.abspath(os.path.join(root,name))
	
	
	
	return path
#---------------------------------------------------------------------------------------------------------------------------
#This function will run MACS program through system call.
def callingMACS(wholeGenomeFilePath,macsOutputName,cutoff):
	subprocess.call(["macs2","bdgpeakcall","-i",wholeGenomeFilePath,"-c",str(cutoff),"-o",str(macsOutputName)])
	
	path=os.path.abspath(macsOutputName)
	
	return path
#---------------------------------------------------------------------------------------------------------------------------
#This function is converting the output of MACs(narrowPeak) file to a BED version of it (more like a SCRMshaw output)
def peaksToScrms(macsOutputPath,peaksToScrmsName,TSET):

	with open(macsOutputPath,'r') as infile, open(os.path.join(subdirectory,peaksToScrmsName),'w') as outfile:
		for line in infile:
			if not line.startswith('track'):
				col=line.split('\t')
				if '_' not in TSET:
					setAndmeth=col[3].split('_')
					set=setAndmeth[0]
					meth=setAndmeth[1]
					amplitude=str(int(col[4])/10)
			
			#mapping1.adult_mesoderm_imm_narrowPeak5
				else:
					setAndmeth=col[3].split('_')
					set=setAndmeth[0]+'_'+setAndmeth[1]
					meth=setAndmeth[2]
					amplitude=str(int(col[4])/10)	
					
				outfile.write(col[0]+'\t'+col[1]+'\t'+col[2]+'\t'+amplitude+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+set+'\t'+meth+'\n')
				
	#path=os.path.abspath(peaksToScrmsName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==peaksToScrmsName:
				path=os.path.abspath(os.path.join(root,name))
	return path
	
#---------------------------------------------------------------------------------------------------------------------------
#This function will retrieve the SCRMshaw score of the peaks by intersecting it with SCRMshaw prediction file
def intersect_peaks_and_scrms(peaksToScrmsPathBED,extractedScrmsPathBED,intersectedFileName):
	peaksToScrmsPathBED.intersect(extractedScrmsPathBED,loj=True).merge(c=[4,20,6,7,8,9,10,11,12,13,14,15,16],o=['max','max','distinct','distinct','distinct','distinct','distinct','distinct','distinct','distinct','distinct','distinct','distinct']).saveas(subdirectory+'/'+intersectedFileName)
	
#	path=os.path.abspath(intersectedFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==intersectedFileName:
				path=os.path.abspath(os.path.join(root,name))

	return path
	
#---------------------------------------------------------------------------------------------------------------------------
#This function will sort the peaks output file according to the amplitude of their peaks and rank it based on that.
def sortAndRank_basedOnAmplitude(intersectedFilePath,finalPeaksFileName):
	data=pd.read_csv(intersectedFilePath,delimiter='\t',header=None)
	numOfPeaks=len(data)+1
	ranks=pd.Series(range(1,numOfPeaks))
	dataSorted=data.sort_values(by=3,ascending=False)
	dataSorted=dataSorted.reset_index(drop=True)
	dataSorted[17]=ranks
	#finalName='peaksFinal_'+TSET+'_'+method
	dataSorted.to_csv(finalPeaksFileName,sep='\t',index=False,header=False)
	
	path=os.path.abspath(finalPeaksFileName)
	return path,numOfPeaks



#############################################-------MAIN FUNCTION-----##########################################################

def main():
	#some of the variables are declared global because these are being used in main program as well as some functions
	global d  
	global d2
	global countd
	global subdirectory
	global patternRecovery
	global totalNumberOfCrmsKnownToCauseExpression
	
	totalNumberOfCrmsKnownToCauseExpression=0
	t=1 #count of sets done
	#Temporary directory in which all the intermediate files will be moved
	subdirectory='tmp' 

	#command line parsing files from users
	parser=argparse.ArgumentParser()
	parser.add_argument('-so','--scrmJoinedOutputFile',help='Scrmshaw Output file concatenated ',required=True)
	parser.add_argument('-num','--numOfScrms',help='Number of Scrms to start from, default is 5000',default=5000)
	args = parser.parse_args()
	scrmJoinedOutputFile=args.scrmJoinedOutputFile
	numOfScrms=args.numOfScrms
	numOfScrms=int(numOfScrms)
	my_path=os.getcwd()
	if not os.path.isdir(my_path+'/'+subdirectory):
		os.makedirs(my_path+'/'+subdirectory)
		
	#iterating through each keys (different training sets) of the three method's dictionary:
	methods=['imm','hexmcd','pac']
	
	#Creating list to iterate through three methods 
	methods_val=[None,None,None]
	num=0
	x=0
		
	#Parsing the output file into separate files for each training set and each method via creating three dictionaries for each method: keys being the names of training sets associated with that method 
	scrmJoinedOutputFile=os.path.abspath(scrmJoinedOutputFile)
	methods_val[0],methods_val[1],methods_val[2]=parse_output(scrmJoinedOutputFile,35000)

	#this loop is used to iterate through three methods dictionaries
	for num in range(len(methods)):
	
		#removing empty string from list
		while '' in methods_val[num]:
  	  		methods_val[num].remove('')
		
		print("Now method:"+methods[num])
		print("methods_Val of num")
		print(methods_val[num])
		#pprint.pprint(methods_val[num])
		
		#loop is used to iterate through all the training sets in the method		
		for x in methods_val[num]:
			TSET=x
			method= methods[num]
			print("Tset: "+TSET)
			print("Method: "+method)
			
			#individual training set's full length scrmshaw predictions file
			file= method+"_"+TSET+"_fullLength.bed"

			#getting the path of this file
			for root, dirs, files in os.walk(os.getcwd()):
				for name in files:
					if name==file:
						#print(name)
						my_path=os.path.abspath(os.path.join(root,name))
		
			#print("Path of your file: is"+ my_path)
			fileFullLengthPath=my_path		
		
			#extract num of scrms from offsets combined output
			extractedScrmsPath=extract_topN_scrms(fileFullLengthPath,numOfScrms,method,TSET)
						
			#calculating the min score to use as cutoff for macs from above file .
			process = subprocess.Popen(["sort","-k4",extractedScrmsPath], stdout=subprocess.PIPE)
			output = process.communicate()[0]
			tmp=output.split('\n')
			tmp1=tmp[0].split('\t')	 
			minScore= tmp1[3]
		
			#sorting the extracted individual offset combined file.
			#first converting it to bed
			extractedScrmsPathBED=bed_conversion(extractedScrmsPath)
			sortedFileName='sorted_'+str(numOfScrms)+'.'+method+'_'+TSET+'.bed'
			sortedFilePath=bedtools_sorting(extractedScrmsPathBED,sortedFileName)
			sumOfScoreFileName=str(numOfScrms)+'.'+'original_'+method+'_'+TSET+'.bdg'
			
			# 	calculating sum of score for each 10bp window 
			sumOfScoreFilePath=sum_of_score(sortedFilePath,sortedFileName,sumOfScoreFileName)

			#sorting the above created file
			sumOfScoreFilePathBED=bed_conversion(sumOfScoreFilePath)
			sumOfScoreFileName='sorted_'+sumOfScoreFileName
			sortedSumOfScoreFilePath=bedtools_sorting(sumOfScoreFilePathBED,sumOfScoreFileName)
			
			#filling out the missing windows across the genome with the score with 0.0 (required step for MACs)
			wholeGenomeFileName='whole_genome_coverage_'+sumOfScoreFileName
			wholeGenomeFilePath=filling_missing_coords(sortedSumOfScoreFilePath,wholeGenomeFileName)
			
			#calling MACs
			macsOutputName=TSET+'_'+method
			macsOutputPath=callingMACS(wholeGenomeFilePath,macsOutputName,minScore)
			

			#converting MACS output narrowPeak to BED file (like SCRMshaw prediction file)
			peaksToScrmsName='peaks_'+macsOutputName
			peaksToScrmsPath=peaksToScrms(macsOutputPath,peaksToScrmsName,TSET)
	

			#intersecting peaks output file with the SCRMshaw prediction file and get the SCRMshaw score for each of the peak.

			#converting both to pybed version to perform bedtools functions
			peaksToScrmsPathBED=bed_conversion(peaksToScrmsPath)
			extractedScrmsPathBED=bed_conversion(extractedScrmsPath)
			#intersecting
			intersectedFileName='Intersected_ScrmsAndPeaks_'+TSET+'_'+method+'.bed'
			intersectedFilePath=intersect_peaks_and_scrms(peaksToScrmsPathBED,extractedScrmsPathBED,intersectedFileName)
			
			#reading intersected file as pandas dataframe to sort based on amplitude	
			numOfPeaks=peaksToScrmsPathBED.count()
			finalPeaksFileName='scrmshawOutput_peaksCalled_'+TSET+'_'+method+'_'+str(numOfPeaks)+'_peaks.bed'
			finalPeaksFilePath,numOfpeaks2=sortAndRank_basedOnAmplitude(intersectedFilePath,finalPeaksFileName)
			print('Number of peaks for the set '+TSET+'_'+method+': '+str(numOfpeaks2))
		
			#moving the extra files to tmp directory
			for root, dirs, files in os.walk(os.getcwd(),topdown=False):
				for name in files:
					if name==TSET+'_'+method or name=='sumScoredsorted_'+str(numOfScrms)+'.'+method+'_'+TSET+'.csv':
						shutil.move(name, 'tmp/')
		
			
main()										
