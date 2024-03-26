#!/usr/bin/env python3
# Maximilian Genetti (maxgenetti@gmail.com)

pops = ['pr','ma']
sites = []
with open("sites.panels") as fileS :
        for line in fileS :
                chrom, pos = line.rstrip().split()
                sites.append([chrom,pos])
for pop in pops :
        values = []
        for i in range(0,len(sites)):
                pi = 0
                try:
                        with open(f'./{pop}/{i+1}.saf.sfs') as fileP :
                                line = fileP.readline()
                                l = [float(x) for x in line.rstrip().split()]
                                n = len(l)-1
                                pi = sum(i * (n - i) / (n * (n - 1)) * x for i, x in enumerate(l[1:-1]))/sum(l)
                                fileP.close()
                except :
                        print(f"broke on {pop}/{i+1}.saf.sfs")
                        quit()
                values.append(sites[i]+[str(pi)])
        with open(f"{pop}.pi", "w") as filePi :
                for value in values :
                        filePi.write('\t'.join(value)+"\n")
                filePi.close()
