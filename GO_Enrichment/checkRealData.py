#!/usr/bin/env python3
# Maximilian Genetti (mgenetti)

import time
import os
import sys
import argparse
import random
from collections import defaultdict
import intervaltree 

'''
'''

class PeakSimulation :
    '''
    Args:

    Returns:

    '''

    def __init__ (self) :
        '''constructor: saves attributes'''
        args = self.get_args() #reads command line
        self.popName = f"../{args.pop}"
        self.direction = args.d
        self.out = f"real-{args.d}.{args.go}"
        self.peaksName = f"peaks-{args.d}.txt"
        self.annotations = f"./annotations/hgd-ranges.{args.go}.txt"
        self.goTerms = defaultdict(list)
        self.chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16']
        self.geneRanges = {}
        self.peaks = {}
        self.genes0 = {}
        self.genes1 = {}
        self.genes2 = {}
        for chrom in self.chroms :
            self.geneRanges[chrom] = intervaltree.IntervalTree()
            self.peaks[chrom] = defaultdict(list)
            self.genes0[chrom] = []
            self.genes1[chrom] = []
            self.genes2[chrom] = []

    def get_args(self):
        '''
        Notes:
            reads command line arguments
        '''
        parser = argparse.ArgumentParser(description='Read in files and paremeters')
        parser.add_argument('--pop', type=str, help='name of population directory')
        parser.add_argument('--d', type=str, default='fix later', help='name of peaks tsv')
        parser.add_argument('--go', action='store_const', const="mol", default='bio', help='option to use default peaks')
        return parser.parse_args()

    def readGenes (self) :
        with open(self.annotations) as fileG :
            for line in fileG :
                geneID = line.split("\t")[0]
                geneName = line.split("\t")[1]
                chrom = line.split("\t")[2]
                pos1 = min(int(line.split("\t")[3]),int(line.split("\t")[4]))
                pos2 = max(int(line.split("\t")[3]),int(line.split("\t")[4]))
                terms = line.rstrip().split("\t")[5:]
                self.geneRanges[chrom].addi(pos1, pos2+1, geneID) #inclusive store ending position
                self.goTerms[geneID] = [geneName, terms]
            fileG.close()

    def readPeaks (self) :
        with open(f"{self.popName}/{self.peaksName}") as fileP :
            for line in fileP :
                chrom = str(line.split()[0])
                start= int(line.split()[1])
                peak = int(line.split()[2])
                stop = int(line.split()[3])
                coef = float(line.split()[4])
                ll = float(line.split()[5])
                self.peaks[chrom][peak]=[peak, start, stop, coef, ll]
            fileP.close()

    def matchGoTerms (self) :
        for chrom in self.chroms : #iterate through chroms
            for pos1, peak in self.peaks[chrom].items() : #iterate through peaks
                pos, start, stop, coef, ll = peak
                for gene in self.geneRanges[chrom][pos-5e5:pos+5e5+1] :
                    if gene.begin <= pos <= gene.end :
                        dist = 0
                    else :
                        if abs(gene.begin-pos) > abs(gene.end-pos) :
                            dist = gene.end-pos
                        else :
                            dist = gene.begin-pos
                    self.genes0[chrom].append([gene.data,dist,pos])

                #closest to peak
                bestDist = 5e5
                dist = 5e5
                bestGene = ''
                zero = False
                for gene in self.geneRanges[chrom][pos-5e5:pos+5e5+1] :
                    if gene.begin <= pos <= gene.end : #within peak
                        bestDist = 0
                        self.genes1[chrom].append([gene.data,bestDist,pos])
                        zero = True
                    else:
                        if abs(gene.begin-pos) > abs(gene.end-pos) :
                            dist = gene.end-pos
                        else :
                            dist = gene.begin-pos
                        if abs(dist) < abs(bestDist) :
                            bestDist = dist
                            bestGene = gene.data
                        else :
                            continue
                if not zero :
                    self.genes1[chrom].append([bestGene,bestDist,pos])

                for gene in self.geneRanges[chrom][start:stop+1] :
                    if gene.begin <= pos <= gene.end :
                        dist = 0
                    else :
                        if abs(gene.begin-pos) > abs(gene.end-pos) :
                            dist = gene.end-pos
                        else :
                            dist = gene.begin-pos    
                    self.genes2[chrom].append([gene.data,dist,pos])

        f0 = open(f"{self.popName}/{self.out}.0.txt", "w") #all within range
        for chrom in self.chroms :
            for gene, dist, peak2 in self.genes0[chrom] :
                peak = self.peaks[chrom][peak2]    # [peak, start, stop, coef, ll])
                geneName, terms = self.goTerms[gene]
                newLine =f"{gene}\t{geneName}\t{chrom}\t{peak[0]}\t{dist}\t{peak[3]}\t{peak[4]}\t{','.join(terms)}\n"
                f0.write(newLine)
        f0.close()
        f1 = open(f"{self.popName}/{self.out}.1.txt", "w") #closest gene to peak only
        for chrom in self.chroms :
            for gene, dist, peak2 in self.genes1[chrom] :
                peak = self.peaks[chrom][peak2]    # [peak, start, stop, coef, ll])
                geneName, terms = self.goTerms[gene]
                newLine =f"{gene}\t{geneName}\t{chrom}\t{peak[0]}\t{dist}\t{peak[3]}\t{peak[4]}\t{','.join(terms)}\n"
                f1.write(newLine)
        f1.close()
        f2 = open(f"{self.popName}/{self.out}.2.txt", "w") #all within measured window 
        for chrom in self.chroms :
            for gene, dist, peak2 in self.genes2[chrom] :
                peak = self.peaks[chrom][peak2]    # [peak, start, stop, coef, ll])
                geneName, terms = self.goTerms[gene]
                newLine =f"{gene}\t{geneName}\t{chrom}\t{peak[0]}\t{dist}\t{peak[3]}\t{peak[4]}\t{','.join(terms)}\n"
                f2.write(newLine)
        f2.close()

def main() :

    test = PeakSimulation()
    test.readPeaks()
    test.readGenes()
    test.matchGoTerms()

if __name__ == "__main__":
    main()
