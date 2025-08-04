#!/usr/bin/env python3
# Maximilian Genetti (maxgenetti@gmail.com)

import numpy as np
import scipy
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from collections import defaultdict
import numpy as np 
import sys
import math

num_sites = 67329
populations = ["PR","MA"]
directions = ["rev","forw"]
dataDir = ''
if dataDir == '' :
    print("Must specify location of data",file = sys.stderr)
    quit()

chroms= {
        '1':29893408,
        '2':15549267,
        '3':13234341,
        '4':12718334,
        '5':14363272,
        '6':18472937,
        '7':13219345,
        '8':13546544,
        '9':11120453,
        '10':12965953,
        '11':14726556,
        '12':11902654,
        '13':10288499,
        '14':10253655,
        '15':10167229,
        '16':7207165
        }

lnl_cutoff = scipy.stats.chi2.ppf(1 - (0.05/num_sites),1)

for direction in directions :
    for population in populations :
        peaks = []
        offset = 0
        full_lnl = []
        full_sel = []
        full_pos = []
        
        # Generate sample data
        for chrom in chroms.keys() : #['1','14'] : #
            pos_values = []
            gen_pos = []
            sel_values = []
            lnl_values = []

            with open(f'/{dataDir}/{population}/ahmms-{direction}-chr{chrom}.out') as fileP :
                for line in fileP :
                    l = line.split()
                    pos_values.append(int(l[0]))
                    gen_pos.append(int(l[0])+offset)
                    sel_values.append(float(l[1]))
                    lnl_values.append(max(float(l[2]),0))
                fileP.close()

            ahmms_len=len(lnl_values)

            lnl_prom = 0.5*np.array(lnl_values)            

            p1, p2 = find_peaks(lnl_values, height=lnl_cutoff, prominence=lnl_prom, distance = 20, width = 21)

            for p in range(0,len(p1)) :
                diff_a = int((pos_values[int(p2["left_ips"][p])+1]- pos_values[int(p2["left_ips"][p])])* (p2["left_ips"][p]%1))
                diff_b = int((pos_values[int(p2["right_ips"][p])+1]-pos_values[int(p2["right_ips"][p])])*(p2["right_ips"][p]%1))
                a = int(p2["left_ips"][p])  #inclusive
                b = int(p2["right_ips"][p]) #inclusive

                if lnl_values[p1[p]] > lnl_cutoff :
                    rh1=1-(lnl_cutoff/lnl_values[p1[p]])
                    prom = p2['prominences'][p]
                    lnl = lnl_values[p1[p]]
                    rh=min(1-(lnl_cutoff-(lnl-prom))/prom, 5000)
                    pw = peak_widths(lnl_values, [p1[p]], rel_height=rh)
                    diff_c = int((pos_values[int(pw[2][0])+1]-pos_values[int(pw[2][0])]) * (pw[2][0]%1))
                    diff_d = int((pos_values[int(pw[3][0])+1]-pos_values[int(pw[3][0])]) * (pw[3][0]%1))
                    c = int(pw[2][0])
                    d = int(pw[3][0])

                    peaks.append([chrom,pos_values[p1[p]],sel_values[p1[p]],lnl_values[p1[p]],pos_values[a]+diff_a,pos_values[b]+diff_b,pos_values[c]+diff_c,pos_values[d]+diff_d])
        
            offset+=chroms[chrom]
            full_sel += sel_values
            full_lnl += lnl_values
            full_pos += gen_pos

        f = open(f'/{dataDir}/{population}/peaks-{direction}.txt', "w")
        print(f"{population}-{direction} : {len(peaks)}")
        for peak in peaks:
            f.write(f'{peak[0]}\t{peak[6]}\t{peak[1]}\t{peak[7]}\t{peak[2]:3f}\t{peak[3]}\n')
