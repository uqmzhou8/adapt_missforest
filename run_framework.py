#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:54:17 2023

@author: uqmzhou8
"""

#Main script to perform preprocessing of VCF files to generate reference chimeric genomes

if __name__ == '__main__':
    from mfpipeline import *
    #enter file name for the samples to impute
    run_mfimpute('masked_test.vcf')
    
   