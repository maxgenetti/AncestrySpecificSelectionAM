#!/usr/bin/env python3
# Maximilian Genetti (mgenetti)

import pandas as pd
import numpy as np
import sys
from statsmodels.stats.multitest import multipletests

mthds = ['fdr_bh']
mthdDict = {'fdr_bh':"bh"}

fileIn = sys.argv[1]
if fileIn.startswith("summary") :
    quit()

p = pd.read_csv(fileIn+".tsv", sep = "\t")

for mthd in mthds :
    df = p[(p['GO-Term'] != "NA") & (p['Observed'] > 0)]

    left_tail = df['P_Enrichment'].to_list()
    bh_lt = multipletests(left_tail, method=mthd, is_sorted=False, returnsorted=False, alpha=0.05)

    temp = []
    for i in range(len(bh_lt[0])) :
        if bh_lt[0][i] == True :
            temp.append("\t".join([df.iloc[i]['GO-Term'],str(df.iloc[i]['Observed']),str(df.iloc[i]['P_Enrichment']),str(bh_lt[1][i])])+"\n")

    if temp != [] :
        fileO = open(f"{fileIn}.{mthdDict[mthd]}.full", "w")
        fileO.write("GO-Term\tObserved\tP-value\tFDR-Corrected\n")
        for line in temp:
            fileO.write(line)
        print(f"{fileIn}.{mthdDict[mthd]}.full", file = sys.stdout)
    else :
        print(f"{fileIn}.{mthdDict[mthd]}.full Skipped", file = sys.stderr)

    ##################

    df = p[(p['Observed'] > 1) & (p['Maximum'] > 10) &(p['GO-Term'] != "NA")]
    left_tail = df['P_Enrichment'].to_list()

    bh_lt = multipletests(left_tail, method=mthd, is_sorted=False, returnsorted=False, alpha=0.05)

    temp = []
    for i in range(len(bh_lt[0])) :
        if bh_lt[0][i] == True :
            temp.append("\t".join([df.iloc[i]['GO-Term'],str(df.iloc[i]['Observed']),str(df.iloc[i]['P_Enrichment']),str(bh_lt[1][i])])+"\n")

    if temp != [] :
        fileO = open(f"{fileIn}.{mthdDict[mthd]}.subset", "w")
        fileO.write("GO-Term\tObserved\tP-value\t\tFDR-Corrected\n")
        for line in temp:
            fileO.write(line)
        print(f"{fileIn}.{mthdDict[mthd]}.subset", file = sys.stdout)
    else :
        print(f"{fileIn}.{mthdDict[mthd]}.subset Skipped", file = sys.stderr)
