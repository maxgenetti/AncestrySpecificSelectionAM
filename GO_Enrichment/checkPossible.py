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
from statsmodels.stats.multitest import multipletests
import statistics #delete after testing

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
        self.geneCounts = {} #delete after testing
        self.siteCounts = {} #delete after testing
        self.geneSites = {} #delete after testing

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
        parser.add_argument('-g', '--group', type=int, default=1, help='simulation block')
        parser.add_argument('-i', '--direction', type=str, default='rev', help='direction of introgression')
        parser.add_argument('-t', '--wide', action='store_const', const=".wide", default='', help='option to check unfiltered and untrimmed peaks')
        parser.add_argument('-d', '--dist', type=int, default=-1, help='Maximum distance from center of peak to check for genes in KB, must be greater than 1kb, ')
        parser.add_argument('-x', '--seperation', type=int, default=20, help='Minimum distance between peaks in AIMs')
        parser.add_argument('-n', '--nearest', action="store_true", help='Only report nearest gene in range, must specify --dist')
        return parser.parse_args()

    def readPanels (self) :
        '''
        Args:

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
                        chrom = line.split()[0]
                        site = int(line.split()[1])
                        sites = [site]
                        self.totalSites += 1
                        self.siteIndex[chrom] = {}
                        self.siteCounts[chrom] = {}
                        self.siteIndex[chrom][site] = 0
                        self.siteCounts[chrom][site] = 0 #delete after testing
                        self.sitesChroms.append([chrom,site]) #list of all sites for randomization
                        for line in fileP :
                            k += 1
                            site = int(line.split()[1])
                            sites.append(site)
                            self.totalSites += 1
                            self.siteIndex[chrom][site] = k
                            self.sitesChroms.append([chrom,site]) #list of all sites for randomization
                            self.siteCounts[chrom][site] = 0 #delete after testing
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
                peak = int(line.split()[2])
                start= int(line.split()[1])
                stop = int(line.split()[3])
                self.peaks.append(stop-start+1) #inclusive width
                self.peak_widths.append([peak-start,stop-peak])
            fileP.close()
        if self.dist >=0 :
            self.peaks = len(self.peaks)*[self.dist*2+1] #inclusive width
            self.peak_widths = len(self.peaks)*[[self.dist,self.dist]]

    def readGenes (self) :
        for chrom in self.chroms :
            self.geneRanges[chrom] = intervaltree.IntervalTree()
        with open("./annotations/hgd-ranges.txt") as fileG :
            for line in fileG :
                geneID = line.split("\t")[0]
                chrom = line.split("\t")[2]
                pos1 = min(int(line.split("\t")[3]),int(line.split("\t")[4]))
                pos2 = max(int(line.split("\t")[3]),int(line.split("\t")[4]))
                terms = line.split("\t")[5:]
                if self.dist == -1000 and (pos2 < self.sites[chrom][0] or pos1 > self.sites[chrom][-1]) :
                    continue
                elif self.dist != -1000 and (pos2 < self.sites[chrom][0] - self.dist or pos1 > self.sites[chrom][-1] + self.dist) :
                    continue
                else :
                    print(geneID)
                    self.geneRanges[chrom].addi(pos1, pos2+1, geneID) #inclusive store ending position
                self.goTerms[geneID] = terms
                self.geneCounts[geneID] = [chrom,0] #delete after testing
                self.geneSites[geneID] = [] #delete after testing
                for term in terms :
                    self.summary[term] = defaultdict(int)

    def sites2genes (self) :
        for chrom in self.chroms :
            print(chrom, file = sys.stderr)
            self.sitesGenes[chrom] = defaultdict(list)
            for site in self.sites[chrom] :
                for p1, p2 in self.peak_widths :
                    bestgene = []
                    dist = max(p1,p2)+2
                    if not self.nearest :
                        for gene in self.geneRanges[chrom][site-p1:site+p2+1] :
                            self.sitesGenes[chrom][site].append(gene.data) #stores gene for each site
                            self.geneSites[gene.data].append(site)
                    else :
                        for gene in self.geneRanges[chrom][site-p1:site+p2+1] :
                            if gene.begin<=site and gene.end>=site :
                                if dist == 0 :
                                    bestgene.append(gene.data)
                                else :
                                    bestgene = [gene.data]
                                    dist = 0
                            elif min(abs(site-gene.begin),abs(site-gene.end)) < dist:
                                bestgene = [gene.data]
                                dist = min(abs(site-gene.begin),abs(site-gene.end))
                            elif min(abs(site-gene.begin),abs(site-gene.end)) == dist :
                                bestgene.append(gene.data)
                            else :
                                continue
                        self.sitesGenes[chrom][site]+=bestgene #stores gene for each site
                        for gene in bestgene :
                            self.geneSites[gene].append(site)

def main() :
    start_time = time.time()
    test = PeakSimulation()
    test.readPanels()
    print(f"{len(test.chroms)} Chromsomes containing {test.totalSites} AIMS", file = sys.stderr)
    test.readGenes()
    print(f"Annotation file contains {len(list(test.goTerms.keys()))} Genes with {len(set([j for sub in list(test.goTerms.values()) for j in sub]))} GO-terms", file = sys.stderr)

if __name__ == "__main__":
    main()
