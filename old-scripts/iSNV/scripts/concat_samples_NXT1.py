#! /usr/bin/env python

import pandas as pd
import numpy as np
import sys
import os.path

if __name__ == "__main__":
    
    #read in the patient shift file so we have a list of patients
    names=pd.read_csv('patients_all.txt',sep='\t')
    #save the names of the patients so we can loop through them
    names=names.patient
    
    #initialize the concatenated dataframe
    full=pd.DataFrame()
    
    #column names to use in the dataframe
    cols=['pos','var','ref','pval','type','freq','ct_1','ct_2','ct_3','ct_4','extra1','extra2']
    
    for name in names:
    
    	#clear data just in case they file doesn't exist in next round
		data=[]
    	
    	#first check if file exists then save to dataframe
		if os.path.exists(name+'.var.raw.txt.mapped.ref.vc.txt.filtered'):
			data=pd.read_csv(name+'.var.raw.txt.mapped.ref.vc.txt.filtered',sep='\t',header=None,names=cols)
		
			#add a column indicating the patient
			data.insert(0,'patient',name)
		
			#now concat the new patient to the whole dataframe	
			full=pd.concat([full,data],ignore_index=True)
    	
	#now save the final dataframe as a file
    full.to_csv('vphaser_iSNVs_combined_NXT1.txt',sep='\t',index=False)