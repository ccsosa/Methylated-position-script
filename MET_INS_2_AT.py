# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 13:40:44 2020

@author: cami_
"""

import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import sys
import gc

#import copy

path = "E:/JAVERIANA/EPI"
os.chdir(path)
exons_charac = pd.read_excel('Epigenetic Aluminum tolerance regions_CCSA.xlsx',"CONC")
chromosomes = exons_charac['Chromosome'].unique().tolist()
#fn1 = path+"/"+"IR64-Control_mC"
#fn1 = path+"/"+"N22-Control_mC"
fn1 = path+"/"+"PK-Control_mC"

exons_charac["CG"]=0
exons_charac["CGH"]=0
exons_charac["CHH"]=0
#exons_charac2 = copy.deepcopy(exons_charac)
col=list(exons_charac.columns)
exons_charac = exons_charac.to_numpy()

#num_lines = sum(1 for line in open(fn1))


fp= open(fn1)
fp = fp.readlines()
gc.collect() 
for line in tqdm(fp):
    _data =  line.strip().split()
    _x_sub = np.array(np.where(exons_charac[:,3]==_data[0]))[0,]
    _x_sub = _x_sub[[(min(exons_charac[i,[5,6]])-1000) <= int(min(_data[1],_data[2])) <=max(exons_charac[i,[5,6]])+1000  for i in list(_x_sub)]]
    if len(_x_sub)>0:

        for i in range(len(_x_sub)):
            #print(str(_x_sub[i]))
            if _data[6]=="CG":
                exons_charac[_x_sub[i],8] +=1 
            elif _data[6]=="CHG":
                exons_charac[_x_sub[i],9] +=1 
            elif _data[6]=="CHH":
                exons_charac[_x_sub[i],10] +=1 

exons_charac2 = pd.DataFrame(exons_charac)
exons_charac2.columns = col
exons_charac2.to_csv('met_sum_PK_NEW.csv', index=False)
