#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 14:19:45 2025

@author: uqmzhou8
"""

#import all necessary packages
import pandas as pd
from cyvcf2 import VCF
import numpy as np
from missingpy import MissForest

#Function for encoding haplotypes
def vcf_read(namef):
    pos=[]
    genos=[]
    otherinfo= []
    for v in VCF(namef):
        otherinfo.append([v.CHROM, v.REF, v.ALT[0]])
        pos.append(v.POS)
        gts = v.genotype.array().astype(int)
    
        if gts.shape[1] != 3:
            raise IOError('VCF must be for diploids.')
        #genotypes are converted in the following format: 0|0=1 0|1=2 1|0=3 1|1=4
        #missing data will be -3
        gts[:, 2] = 0
        gts[gts[:, 0] == 0, 0] = 2
        gts[gts[:, 0] == 1, 0] = 4
        gts = gts.sum(axis=1)
        gts = gts -1
        genos.append(gts)
    genos=np.stack(genos)
    pos=np.array(pos, dtype= "int")
    otherinfo=np.stack(otherinfo)
    return genos, pos, otherinfo
 
def run_mfimpute(filen):
    
    geno, pos, otherinfo= vcf_read(filen)
    #run missForest
    clf= MissForest(missing_values = -3)
    new_geno= clf.fit_transform(geno, cat_vars = list(range(len(geno[0]))))
    tosto = pd.DataFrame(new_geno)
    
    #save the imputed results back into a VCF file
    header_1= VCF(filen).raw_header
    temp=header_1.split('\n')
    temp=temp[:-2]
    temp=temp+['##Imputed with missForest']
    df= pd.DataFrame([temp])
    df=df.T
    df.to_csv('imputed'+filen , sep="\t",header=None, index=None, quoting=3 ,escapechar="\n")
    columnheader=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    geinfo=np.zeros([len(pos),9])
    geinfo=pd.DataFrame(geinfo, columns= columnheader,dtype=int)
    geinfo=geinfo.replace(0,'.')
    geinfo['POS']= pos
    geinfo['FILTER']='PASS'
    geinfo['FORMAT']='GT'
    geinfo['#CHROM']=otherinfo[:,0]
    geinfo['REF']=otherinfo[:,1]
    geinfo['ALT']=otherinfo[:,2]
    geno2sto=pd.DataFrame(tosto)
    geno2sto=geno2sto.replace(1,'0|0')
    geno2sto=geno2sto.replace(2,'0|1')
    geno2sto=geno2sto.replace(3,'1|0')
    geno2sto=geno2sto.replace(4,'1|1')
    newdf= pd.concat([geinfo,geno2sto],axis=1)
    newdf.to_csv('imputed'+filen, mode='a', sep="\t", header=True, index=None, quoting=3 ,escapechar="\n")

