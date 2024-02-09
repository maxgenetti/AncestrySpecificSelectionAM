import os
from collections import defaultdict
import argparse

def read_to_list(filename):
    with open(filename, 'r') as file:
        lines = [line.rstrip().split()[0] for line in file.readlines()]
    return lines

def writeCommands(samples, chroms, first, tmp, work) :
    x = first
    for chrom in chroms:
        y = str(x)+'.sh'
        with open(y,"w") as f :
            f.write("gatk GenomicsDBImport \\\n")
            f.write(f"\t--genomicsdb-workspace-path {work}/{chrom.split(".")[0]} \\\n")
            for sample in samples:
                f.write(f'\t-V {sample} \\\n')
            f.write(f'\t--tmp-dir {tmp}/{chrom.split(".")[0]} \\\n')
            f.write('\t-L %s \n' %(chrom))
            x += 1

def main():
    parser = argparse.ArgumentParser(description='Read lines from a file into a list.')
    parser.add_argument('--samples', type=str, help='Samples with paths to include in VCF')
    parser.add_argument('--chroms', type=str, help='Chromsomes to include in VCF')
    parser.add_argument('--format', type=str, help='Path to tmp directory, working directory')
    parser.add_argument('--first', type=int, default = 1, help='starting number for job array')

    args = parser.parse_args()
    samples = read_to_list(args.samples)
    chroms = read_to_list(args.chroms)
    tmp, work = read_to_list(args.format)
    writeCommands(samples, chroms, args.first, tmp, work)

if __name__ == '__main__':
    main()