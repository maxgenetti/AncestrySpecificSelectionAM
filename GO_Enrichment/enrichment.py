#!/usr/bin/env python3
# Maximilian Genetti (mgenetti)
import pickle
from math import ceil
from math import floor
import os
import sys
import argparse
from collections import defaultdict
import time
import gzip
import numpy as np


class CountGOTerms :

    def __init__ (self) :
        '''constructor: saves attributes'''
        args = self.get_args() #reads command line
        self.chroms = {'1':[2714,29889910],'2':[2885,15493755],'3':[2562,13231344],'4':[3675,12714835],'5':[2555,14361741],'6':[3723,18471248],'7':[2567,13217654],'8':[2517,13497559],'9':[5606,10994606],'10':[2535,12913265],'11':[2665,14385295],'12':[4229,11897115],'13':[3172,10286922],'14':[2512,10244012],'15':[55735,9994945],'16':[3961,7151123]}
        self.pop = args.pop
        self.popName = f"../{args.pop}"
        self.sim = args.sim
        self.annotations = f"./annotations/hgd-ranges.{args.go}.txt"
        self.genes = defaultdict(int)
        self.goData = defaultdict(int)
        self.goSim = defaultdict(dict)
        self.v = args.verbose
        self.go = args.go
        self.direction = args.d
        self.window = args.window
        self.tests = 1000
        self.max = 0
        self.goMax = defaultdict(int)
        self.goPeaks = defaultdict(set)
        self.dup = args.dup
        self.nearest = args.nearest
        self.totalData = 0
        if self.nearest :
            self.out = f"./{self.go}/{args.pop}-{args.d}.{self.go}.{args.window}.tsv"
            self.infile = f'./{self.popName}/real-{self.direction}.{self.go}.1.txt'

        elif self.window == -1 :
            self.out = f"./{self.go}/{args.pop}-{args.d}{self.go}.range.tsv"
            self.infile = f'./{self.popName}/real-{self.direction}.{self.go}.2.txt'
        else :
            self.out = f"./{self.go}/{args.pop}-{args.d}.{self.go}.range{args.window}.tsv"
            self.infile = f'./{self.popName}/real-{self.direction}.{self.go}.0.txt'
        self.pkl = ''

    def get_args(self):
        parser = argparse.ArgumentParser(description='Read in files and paremeters')
        parser.add_argument('--pop', type=str, help='name of population directory')
        parser.add_argument('--sim', type=int, default=1000000, help='number in each sim')
        parser.add_argument('--dup', type=int, default=1, help='number of simes to read')
        parser.add_argument('--verbose', default=False, action='store_true', help='output GO terms found in data be default, or all go terms found in genome')
        parser.add_argument('--window', type=int, default=-2, help='window in KB used in newSim')
        parser.add_argument('--d', type=str, default='rev', help='direction')
        parser.add_argument('--go', action='store_const', const="mol", default='bio', help='option to thin on minimum selCoef')
        parser.add_argument('-n', '--nearest', action="store_true", help='Only report nearest gene in range')
        return parser.parse_args()

    def readGoTerms (self) :
        with open(self.annotations) as fileG :
            for line in fileG :
                chrom, p1,p2 = line.rstrip().split()[2:5]
                if int(p2)<self.chroms[chrom][0] or int(p1)>self.chroms[chrom][1] :
                    continue
                else :
                    terms = line.rstrip().split()[5:]
                    for term in terms :
                        self.goSim[term] = defaultdict(int) #makes goSim inclusive of desired terms
                        self.goMax[term] += 1
            fileG.close()

    def readGoSim (self) :
        j = 0
        for i in range(1,self.dup+1) :
            try : #read only first complete file
                if self.nearest :
                    self.pkl = f"./pickles/{self.pop}-{self.direction}.{self.go}.{self.window}.{i}.1.pkl"
                elif self.window == -1 :
                    self.pkl = f"./pickles/{self.pop}-{self.direction}.{self.go}.range.{i}.2.pkl"
                else :
                    self.pkl = f"./pickles/{self.pop}-{self.direction}.{self.go}.{self.window}.{i}.2.pkl"

                with open(self.pkl, 'rb') as fileP:
                    loaded_dict = pickle.load(fileP)
                    for goRaw, freqDict in loaded_dict.items() :
                        go = goRaw.rstrip()
                        if go in self.goSim.keys() : # included in desired goTerms
                            for freq in freqDict.keys() :
                                count = freqDict[freq]
                                self.goSim[go][freq] += count
                                self.max = max(freq,self.max)
                    fileP.close()
                j += 1
                if j == self.dup :
                    break
            except :
                continue
        if j < self.dup :
            print(f"Error reading simulations for {self.infile.split('/')[-1]} with window of {self.window}", file = sys.stderr )
            quit()
        for term in self.goSim.keys() :
            self.goSim[term].pop("0", None) #deleting zero instance just in case it's string formatted

    def readGoData (self) :
        with open(self.infile) as fileT:
            genes = []
            for line in fileT :
                l = line.rstrip().split("\t")
                if l[0] in genes :
                    continue
                elif self.window == -1 :
                    genes.append(l[0])
                    for term in l[7].split(','):
                        self.goData[term]+=1
                        self.goPeaks[term].add((l[2],l[3]))
                elif -self.window*1e3 <= int(l[4]) <= self.window*1e3 :
                    genes.append(l[0])
                    for term in l[7].split(','):
                        self.goData[term]+=1
                        self.goPeaks[term].add((l[2],l[3]))
                else :
                    print(f'Too Far: {line}') 
        fileT.close()

    def compare (self) :
        obs = ["SIM"+str(x) for x in list(range(0,self.max+1))]
        final_array = [['GO-Term','Observed','Peaks','P_Enrichment','Maximum']+obs]
        terms = list(self.goSim.keys())
        terms.sort()
        totalSimCount = 0
        for term in terms :
            simCount = sum(list(self.goSim[term].values())) #total sims term was observed
            totalSimCount += simCount
            value = self.goData[term]  #real data observations
            self.totalData += value
            self.goSim[term][0] = (self.sim*self.dup)-simCount  #add in zero observations to sims
            pairs = []
            maxObserved = max([int(x) for x in self.goSim[term].keys()])
            observed = list(range(0,max(value,maxObserved)+2))
            if simCount == 0 and value > 0 :
                print(f"check {term} in {self.pkl.split('/')[-1]}")
            if value in observed : #report range the frequnecy was observed in the newSim
                idx = observed.index(value)
                try:
                    first = sum([self.goSim[term][x] for x in observed[:idx+1]])/(self.sim*self.dup) #frequency of same or smaller
                    last  = sum([self.goSim[term][x] for x in observed[idx:]])/(self.sim*self.dup) #frequency of same or larger P_Enrichment
                except :
                    print(f'{term} : {self.goSim[term]}')
                    print(value)
                    print(observed)
                    print(maxObserved)
                    break
            else :
                print('fuck', file = sys.stderr)
            if self.goSim[term][self.goMax[term]+1] != 0 :
                print(f'{term} was seen {self.goMax[term]+1} times in {observed[self.goMax[term]+1]} sims, more than possible', file = sys.stderr)
                print(term)
                print("\t".join([term,str(value),str(len(self.goPeaks[term])),str(last),str(self.goMax[term])]+pairs))
                quit()
            pairs = [str(self.goSim[term][x]) for x in range(0, self.max+1)]
            final_array.append([term,str(value),str(len(self.goPeaks[term])),str(last),str(self.goMax[term])]+pairs)
        np.savetxt(self.out, np.array(final_array), delimiter='\t', fmt='%s')
        print(f'{self.totalData} {self.go} in {self.infile.split("/")[-1]} compared to {self.pkl.split("/")[-1]} with {totalSimCount/(self.sim*self.dup)} {self.go} per sim and printed to {self.out.split("/")[-1]}', file = sys.stderr)

def main() :
    test = CountGOTerms()
    test.readGoTerms()
    test.readGoData()
    #print(f'READ DATA after {time.time()-start_time}', file = sys.stderr)
    test.readGoSim()
    #print(f'READ SIMS after {time.time()-start_time}', file = sys.stderr)
    test.compare()
    #print(f'{test.infile.split("/")[-1]} compared to simulation {test.pkl.split("/")[-1]} and printed to {test.out.split("/")[-1]}', file = sys.stderr)
if __name__ == "__main__":
    main()
