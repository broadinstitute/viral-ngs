#! /usr/bin/env python

import pandas as pd
import numpy as np
import sys


def validate(row1, row2):
	'''
	This verifies that two libraries describe the same variant at a given position.
	If it is a SNP, we require the same two alleles, but in any order.
	If it is an indel, we require both to be indels, but do not check the alleles.
	'''
	
	alleles1 = [row1['ref'], row1['var']]
	alleles2 = [row2['ref'], row2['var']]
	alleles1 = [a[0].upper() for a in alleles1]
	alleles2 = [a[0].upper() for a in alleles2]
	return set(alleles1)==set(alleles2)
	

if __name__ == "__main__":

	'''
	This compares two replicate libraries and looks for positions that match.
	Used for validation of the batch2 Nextera library.
	'''	
	
	#get file names from the command line
	file1 = str(sys.argv[1])
	file2 = str(sys.argv[2])
	#get the files as dataframes
	rep1=pd.read_csv(file1,sep='\t')
	rep2=pd.read_csv(file2,sep='\t')
    
    #save list of positions
	names=set(rep1.pos)
	#save list of acceptable positions, to be filled in
	pos = set()
    
	for name in names:
    
    	#get all the rows with a particular position
		temp_rep1 = rep1[rep1.pos == name]
		temp_rep2 = rep2[rep2.pos == name]
		
		#create a list of patients
		patients = set(temp_rep1.patient.tolist())
		
		escape = False
		
		for p in patients:
		
			#get the rep1 row for a given patient
			rep1_row = temp_rep1[temp_rep1.patient == p]
			#ensure that we have exactly one row in the first library given this patient
			#assert rep1_row.shape[0]==1
			#save rep1_row as an object
			for index,row in rep1_row.iterrows():
				rep1_row = row
				
			#get the rep2 row for a given patient
			rep2_row = temp_rep2[temp_rep2.patient == p]
			
			#continue if we have a matching patient
			if rep2_row.shape[0]>0:
			
				#save rep2_row as an object
				for index,row in rep2_row.iterrows():
					rep2_row = row
			
				#look for matching alleles
				if validate(rep1_row, rep2_row):
					#if we find a position with matching alleles, add the position and exit
					pos.add(name)
					escape = True
					break
					
			if escape:
				break
		
	#now use pos to filter out positions that dont have duplicates
	rep1=rep1[(rep1.pos.isin(pos))]

	#save resulting dataframe to a file
	rep1.to_csv(file1.split('.')[0]+'_validated.txt',sep='\t',index=False)
    
    