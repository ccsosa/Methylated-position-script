import numpy as np
import pandas as pd
import allel
from tqdm import tqdm
import gc
dir_i = "/users/ccsosaa/SNP/"



names = ["CX230.snp.vcf.gz","IRIS_313-8244.snp.vcf.gz","CX140.snp.vcf.gz"]
names2 = ["IR_64","PK","NIPPONBARE"]

#def VCF_COUNT(names,names2,exons_charac,dir_i):
for j in range(len(names)):
    exons_charac = pd.read_excel(dir_i+'Epigenetic Aluminum tolerance regions_CCSA.xlsx',"CONC")
    chromosomes = exons_charac['Chromosome'].unique().tolist()

    exons_charac["CHR"][exons_charac["CHR"]=="ChrUn"] = "chr01"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr2"] = "chr02"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr3"] = "chr03"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr4"] = "chr04"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr5"] = "chr05"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr6"] = "chr06"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr7"] = "chr07"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr8"] = "chr08"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr9"] = "chr09"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr10"] = "chr10"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr11"] = "chr11"
    exons_charac["CHR"][exons_charac["CHR"]=="Chr12"] = "chr12"
    chromosomes = exons_charac['CHR'].unique().tolist()
    exons_charac["AL_COUNT"]=0
    exons_charac["AL_SUM"]=0
    exons_charac["AL_HET_COUNT"]=0
    exons_charac["AL_HET_SUM"]=0
    

    """
    # iterating the columns 
    for col in exons_charac.columns: 
        print(col) 
    """

    PK_VCF = allel.read_vcf(dir_i+names[j])
    #sorted(PK_VCF.keys())
    gt = allel.GenotypeArray(PK_VCF['calldata/GT'])
    ac = gt.count_alleles()

    CHR = PK_VCF['variants/CHROM']
    POS = PK_VCF['variants/POS']
    COUNT = ac #np.sum(ac,axis=1)
    HET = gt.count_het(axis=1)

    del PK_VCF
    del gt
    del ac
    gc.collect()

    for i in tqdm(range(exons_charac.shape[0]),desc="gen"):
        _x = np.where(CHR == exons_charac.iloc[i,3])
        _ch = POS[_x]
        _ch3= np.where(np.logical_and(int(min(exons_charac.iloc[i,[5,6]]))<=_ch,_ch<=int(max(exons_charac.iloc[i,[5,6]]))))
        #POS[_ch3]
        exons_charac.iloc[i,8] = np.count_nonzero(COUNT[_ch3])  
        exons_charac.iloc[i,9] = np.sum(COUNT[_ch3])  
        exons_charac.iloc[i,10] = np.count_nonzero(HET[_ch3])  
        exons_charac.iloc[i,11] = np.sum(HET[_ch3]) 
    exons_charac.to_csv(dir_i+names2[j]+'.csv', index=False)
    print("DONE")

"""
for element in dir():
    if element[0:2] != "__":
        del globals()[element]

del element

"""