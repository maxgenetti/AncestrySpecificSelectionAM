#!/usr/bin/env python3
# Maximilian Genetti (maxgenetti@gmail.com)

pops = ['pr','ma']

chromosomes = {
    "NC_007070.3": 29893408,
    "NC_007071.3": 15549267,
    "NC_007072.3": 13234341,
    "NC_007073.3": 12718334,
    "NC_007074.3": 14363272,
    "NC_007075.3": 18472937,
    "NC_007076.3": 13219345,
    "NC_007077.3": 13546544,
    "NC_007078.3": 11120453,
    "NC_007079.3": 12965953,
    "NC_007080.3": 14726556,
    "NC_007081.3": 11902654,
    "NC_007082.3": 10288499,
    "NC_007083.3": 10253655,
    "NC_007084.3": 10167229,
    "NC_007085.3": 7207165,
}
windows =[]
with open("sites.panels") as fileP:
    for line in fileP:
        chrom = (f'NC_00{int(line.split()[0])+7069}.3')
        pos = int(line.split()[1])
        windows.append((chrom,pos))

for pop in pops :
    j = 0
    for window in windows:
        j+=1
        start = max(window[1] - 100000, 0)
        end = min(window[1] + 100001, chromosomes[window[0]])
        output_prefix = f"{window[0]}_{window[1]}"
        command = f'/hb/home/mgenetti/bin/angsd/misc/realSFS {pop}.haploid.saf.idx -P 1 -fold 1 -r {window[0]}:{start}-{end} >>{pop}/{j}.saf.sfs'
        print(command)
