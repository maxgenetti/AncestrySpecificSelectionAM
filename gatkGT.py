import os
from collections import defaultdict
import argparse

def read_to_list(filename):
    with open(filename, 'r') as file:
        lines = [line.rstrip().split()[0] for line in file.readlines()]
    return lines

def writeCommands(chroms, first, ref, outdir) :
    x = first
    for chrom in chroms:
        gvcf = sample.split("/")[-1].split(".")[0]
        y = str(x)+'.sh'
        with open(y,"w") as f :
            f.write("gatk GenotypeGVCFs \\\n")
            f.write(f'-V gendb:/{work} \\\n')
            f.write(f"\t-R {ref} \\\n")
            f.write(f'\t-O {outdir}/{chrom} \n')
            x += 1

def main():
    parser = argparse.ArgumentParser(description='Read lines from a file into a list.')
    parser.add_argument('--chroms', type=str, help='Chromsomes to include in VCF')
    parser.add_argument('--format', type=str, help='Path to reference and working directory')
    parser.add_argument('--first', type=int, default = 1, help='starting number for job array')

    args = parser.parse_args()
    chroms = read_to_list(args.chroms)
    ref, work = read_to_list(args.format)
    writeCommands(samples, chroms, args.first, ref, work)

if __name__ == '__main__':
    main()