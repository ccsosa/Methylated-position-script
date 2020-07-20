import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import sys
import copy
import gc

#import copy

path = "E:/JAVERIANA/EPI"
os.chdir(path)
exons_charac = pd.read_excel('Epigenetic Aluminum tolerance regions_CCSA.xlsx',"CONC")
chromosomes = exons_charac['Chromosome'].unique().tolist()
fn1 = path+"/"+"GSM1039489_sample03_BSseq.txt"

exons_charac["CG"]=0
exons_charac["CGH"]=0
exons_charac["CHH"]=0
exons_charac2 = copy.deepcopy(exons_charac)
col=list(exons_charac2.columns)
exons_charac = exons_charac.to_numpy()

#num_lines = sum(1 for line in open(fn1))
fp= open(fn1)

for line in tqdm(fp):

    _data =  fp.readline().strip().split()
    _x_sub = np.array(np.where(exons_charac[:,4]==int(_data[0])))[0,]
    _x_sub = _x_sub[[min(exons_charac[i,[5,6]])-1000 <= int(_data[1]) <=max(exons_charac[i,[5,6]])+1000  for i in list(_x_sub)]]
    if len(_x_sub)>0:
        for i in range(len(_x_sub)):

            if _data[3]=="CG":
                exons_charac[_x_sub[i],8] +=1 
            elif _data[3]=="CHG":
                exons_charac[_x_sub[i],9] +=1 
            elif _data[3]=="CHH":
                exons_charac[_x_sub[i],10] +=1 


exons_charac2 = pd.DataFrame(exons_charac)
exons_charac2.columns = col
exons_charac2.to_csv('met_sum_NIPPONBARE_NEW.csv', index=False)
