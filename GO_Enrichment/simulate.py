#!/usr/bin/env python3
# Maximilian Genetti (mgenetti)
import gc
import pickle
import time
import os
import sys
import argparse
import math
import random
from collections import defaultdict
import intervaltree 

class PeakSimulation :
    '''
    Args:

    Returns:

    '''

    def __init__ (self) :
        '''constructor: saves attributes'''
        args = self.get_args() #reads command line
        self.popName = f"../{args.population}"
        self.group = args.group
        self.direction = args.direction
        self.annotations = f"./annotations/hgd-ranges.{args.go}.txt"
        self.sims = args.sims
        self.dist = args.dist*1000
        self.seperation = args.seperation
        self.nearest = args.nearest
        self.sites = {}
        self.totalSites = 0
        self.peaks = []
        self.chroms = []
        self.geneRanges = {}
        self.goTerms = {}
        self.summary = {}
        self.peak_widths = []
        self.siteIndex = {}
        self.sitesGenes = {} #only used when checking for nearest genes
        self.chromCount = defaultdict(int) #sampling check
        self.chroms = []
        self.sitesChroms = []

        if self.nearest :
            self.out = f"./pickles/{args.population}-{args.direction}.{args.go}.{args.dist}.{args.group}.1"
        else :
            if self.dist < 0 :
                self.out = f"./pickles/{args.population}-{args.direction}.{args.go}.range.{args.group}.2"
            else :
                self.out = f"./pickles/{args.population}-{args.direction}.{args.go}.{args.dist}.{args.group}.2"
        self.peaksName = f"peaks-{args.direction}.txt"

    def get_args(self):
        '''
        Notes:
            reads command line arguments
        '''
        parser = argparse.ArgumentParser(description='Read in files and paremeters')
        parser.add_argument('-p', '--population', type=str, help='name of population directory')
        parser.add_argument('-s', '--sims', type=int, default=1e6, help='number of simulations to run')
        parser.add_argument('-g', '--group', type=int, default=1, help='simulation group')
        parser.add_argument('-i', '--direction', type=str, default='rev', help='direction of introgression')
        parser.add_argument('-o', '--go', action='store_const', const="mol", default='bio', help='option to check for mol or bio')
        parser.add_argument('-d', '--dist', type=int, default=-1, help='Maximum distance from center of peak to check for genes in KB, must be greater than 1kb, ')
        parser.add_argument('-x', '--seperation', type=int, default=20, help='Minimum distance between peaks in AIMs')
        parser.add_argument('-n', '--nearest', action="store_true", help='Only report nearest gene in range, must specify --dist')
        return parser.parse_args()

    def readPanels (self) :
        '''
        Args:
        Sites used in ancestry_hmm input files
        Notes:

        '''
        for root, dirs, files in os.walk(self.popName) :
            for f in files :
                if f.endswith('.panel') :
                    k = 0
                    with open(self.popName+"/"+f) as fileP :
                        fileP.readline()
                        fileP.readline()
                        line = fileP.readline()
                        l = line.rstrip().split('\t')
                        chrom = l[0]
                        site = int(l[1])
                        sites = [site]
                        self.totalSites += 1
                        self.siteIndex[chrom] = {}
                        self.siteIndex[chrom][site] = 0
                        self.sitesChroms.append([chrom,site]) #list of all sites for randomization
                        for line in fileP :
                            l = line.rstrip().split('\t')
                            k += 1
                            site = int(l[1])
                            sites.append(site)
                            self.totalSites += 1
                            self.siteIndex[chrom][site] = k
                            self.sitesChroms.append([chrom,site]) #list of all sites for randomization
                        self.sites[chrom] = sites
                        self.chroms.append(chrom)
                        fileP.close()
        try:
            self.chroms.sort(key=int)
        except :
            self.chroms.sort()

    def readPeaks (self) :
        with open(f"{self.popName}/{self.peaksName}") as fileP :
            for line in fileP :
                l = line.rstrip().split('\t')
                peak = int(l[2])
                start= int(l[1])
                stop = int(l[3])
                self.peaks.append(stop-start+1) #inclusive width
                self.peak_widths.append([peak-start,stop-peak])
            fileP.close()
        if self.dist >=0 :
            self.peaks = len(self.peaks)*[self.dist*2+1] #inclusive width
            self.peak_widths = len(self.peaks)*[[self.dist,self.dist]]

    def readGenes (self) :
        for chrom in self.chroms :
            self.geneRanges[chrom] = intervaltree.IntervalTree()
        with open(self.annotations) as fileG :
            for line in fileG :
                l = line.rstrip().split('\t')
                geneID = l[0]
                chrom = l[2]
                pos1 = min(int(l[3]),int(l[4]))
                pos2 = max(int(l[3]),int(l[4]))
                terms = l[5:]
                self.geneRanges[chrom].addi(pos1, pos2+1, geneID) #inclusive store ending position
                self.goTerms[geneID] = terms
                for term in terms :
                    self.summary[term] = defaultdict(int)

    def simulate (self) : #store all positions each peak would cover
        percent = max(self.sims // 100, 1)
        start_time2 = time.time()
        max_attempts = 1000 #prevents getting stuck in second while loop if further peak selection becomes impossible.
        i = 0
        while i < self.sims : #should not get stuck, unless peaks are impossible.
            chromosomes_trees = {}
            peak_index = {}
            for chrom in self.chroms : #create empty range to prevent going over chromsomal distances, does not block selection of first or last site
                chromosomes_trees[chrom] = intervaltree.IntervalTree([intervaltree.Interval(0, self.sites[chrom][0]), intervaltree.Interval(self.sites[chrom][-1]+1, self.sites[chrom][-1]+2)]) 
                peak_index[chrom] = [] #dict to store indexes of peaks used and exclude self.seperation sites around them.
            genes = set()
            goCounts = defaultdict(int)
            retry = False
            for p1,p2 in self.peak_widths :
                attempts = 0
                while attempts < max_attempts :
                    chrom, position = random.choice(self.sitesChroms) #all sites are random
                    posIndex = self.siteIndex[chrom][position]
                    if posIndex in peak_index[chrom] : #check conflict in preselected aims for current simulation run
                        attempts += 1
                    else :
                        if not chromosomes_trees[chrom].overlaps(position-p1,position+p2+1): #check if not overlapping other peaks, inclusively
                            chromosomes_trees[chrom].addi(position-p1, position+p2+1)
                            self.chromCount[chrom]+=1
                            peak_index[chrom]+=(list(range(posIndex-self.seperation,posIndex+self.seperation+1))) #add distance in sites around peaks (from peak detection)
                            for gene in self.geneRanges[chrom][position-p1:position+p2+1] : #add overlapping genes, check inclusively
                                genes.add(gene.data)
                            break
                        else :
                            attempts += 1
                if attempts == max_attempts:
                    retry = True
                    break
            if not retry : #store simulation run as GO counts for run. Could check if gene never selected?
                i+=1
                for gene in genes :
                    for term in self.goTerms[gene] :
                        goCounts[term] += 1
                for key,value in goCounts.items() :
                    self.summary[key][value] += 1
                if i % percent == 0 :
                    print(f"{int(i/percent)}% : {time.time() - start_time2} seconds", file = sys.stderr)
                    random.shuffle(self.peaks)
        print(f"{int(i)} Simulations, Printing...", file = sys.stderr)
        with open(f'{self.out}.pkl', 'wb') as f:
            pickle.dump(self.summary, f)
            f.close()

    def sites2genes (self) :
        for chrom in self.chroms :
            self.sitesGenes[chrom] = {}
            for site in self.sites[chrom] :
                dist = self.dist
                bestGene = []
                for gene in self.geneRanges[chrom][site-self.dist:site+self.dist+1] :
                    if gene.begin<=site<=gene.end :
                        if dist == 0 : #multiple overlapping genes possible
                            bestGene.append(gene.data)
                        else :
                            bestGene = [gene.data]
                            dist = 0
                    elif dist > min(abs(site-gene.begin),abs(site-gene.end)) : #check if closer than current best
                        bestGene = [gene.data]
                        dist = min(abs(site-gene.begin),abs(site-gene.end))
                #save best gene(s), could be empty list
                self.sitesGenes[chrom][site] = bestGene #stores gene for each site

    def simulateNearest(self) :
        percent = max(self.sims // 100, 1)
        start_time2 = time.time()
        max_attempts = 10000 #prevents getting stuck in second while loop if further peak selection becomes impossible.
        i = 0
        while i < self.sims : #should not get stuck, unless peaks are impossible.
            peak_index = {}
            for chrom in self.chroms : #create empty intervaltrees
                peak_index[chrom] = [] #dict to store indexes of peaks used and exclude self.seperation sites around them.
            genes = set()
            goCounts = defaultdict(int)
            retry = False
            for p1,p2 in self.peak_widths :
                attempts = 0
                while attempts < max_attempts :
                    chrom, position = random.choice(self.sitesChroms) #all sites are random
                    posIndex = self.siteIndex[chrom][position]
                    if posIndex in peak_index[chrom] : #check conflict in preselected aims for current simulation run
                        attempts += 1
                    else :
                        self.chromCount[chrom]+=1
                        peak_index[chrom]+=(list(range(posIndex-self.seperation,posIndex+self.seperation+1)))
                        for gene in self.sitesGenes[chrom][position] :
                            genes.add(gene)
                        break
                if attempts == max_attempts:
                    retry = True
                    break
            if not retry : #store simulation run as GO counts for run. Could check if gene never selected?
                i+=1
                for gene in genes :
                    for term in self.goTerms[gene] :
                        goCounts[term] += 1
                for key,value in goCounts.items() :
                    self.summary[key][value] += 1
                if i % percent == 0 :
                    print(f"{int(i/percent)}% : {time.time() - start_time2} seconds", file = sys.stderr)
                    random.shuffle(self.peaks)
        print(f"{int(i)} Simulations, Printing...", file = sys.stderr)
        with open(f'{self.out}.pkl', 'wb') as f:
            pickle.dump(self.summary, f)
            f.close()

def main() :
    start_time = time.time()
    test = PeakSimulation()
    test.readPanels()
    print(f"{len(test.chroms)} Chromsomes containing {test.totalSites} AIMS", file = sys.stderr)
    test.readGenes()
    print(f"Annotation file contains {len(list(test.goTerms.keys()))} Genes with {len(set([j for sub in list(test.goTerms.values()) for j in sub]))} GO-terms", file = sys.stderr)
    test.readPeaks()
    print(f"{len(test.peaks)}  Peaks : {','.join([str(x) for x in test.peaks[::-1]])}", file = sys.stderr)
    if test.peaks == [] :
        print(f"No Peaks, runtime : {time.time() - start_time} seconds", file = sys.stderr)
        quit()
    elif test.nearest :
        print("Converting sites to closest gene(s) in range", file = sys.stderr)
        test.sites2genes()
        print(f"Starting {test.sims} simulations for {len(test.peaks)} peaks : {time.time() - start_time} seconds", file = sys.stderr)
        test.simulateNearest()
        print("Complete, runtime : %s seconds" %(time.time() - start_time), file = sys.stderr)
        for chrom in test.chroms :
            print(f'Sample Chromsome {chrom} {test.chromCount[chrom]} times from {len(test.sites[chrom])} AIMs', file = sys.stderr)
    else :
        print(f"Starting {test.sims} simulations for {len(test.peaks)} peaks : {time.time() - start_time} seconds", file = sys.stderr)
        test.simulate()
        print("Complete, runtime : %s seconds" %(time.time() - start_time), file = sys.stderr)
        for chrom in test.chroms :
            print(f'Sample Chromsome {chrom} {test.chromCount[chrom]} times from {len(test.sites[chrom])} AIMs', file = sys.stderr)

if __name__ == "__main__":
    main()
