import os
from collections import defaultdict
import argparse

def read_to_list(filename):
    with open(filename, 'r') as file:
        lines = [line.rstrip().split()[0] for line in file.readlines()]
    return lines

def writeCommands(samples, chroms, first, ploidy, ref, outdir) :
    x = first
    for sample in samples:
        gvcf = sample.split("/")[-1].split(".")[0]
        y = str(x)+'.sh'
        with open(y,"w") as f :
            f.write("gatk HaplotypeCaller \\\n")
            f.write(f"\t-ERC GVCF \\\n")
            f.write(f"\t-ploidy {ploidy} \\\n")
            f.write(f"\t-R {ref} \\\n")
            f.write(f'\t-I {sample} \\\n')
            f.write(f'\t-O {outdir}/{gvcf}.g.vcf.gz \\\n')
            f.write(f'\t-L {chroms} \n')
            x += 1

def main():
    parser = argparse.ArgumentParser(description='Read lines from a file into a list.')
    parser.add_argument('--samples', type=str, help='Samples with paths to include in VCF')
    parser.add_argument('--chroms', type=str, help='Chromsomes to include in VCF')
    parser.add_argument('--format', type=str, help='Reference with path and output directory for gvcfs')
    parser.add_argument('--first', type=int, default = 1, help='starting number for job array')
    parser.add_argument('--ploidy', type=int, default = 2, help='ploidy')

    args = parser.parse_args()
    samples = read_to_list(args.samples)
    ref, outdir = read_to_list(args.format)
    writeCommands(samples, args.chroms, first, ploidy, ref, outdir)

if __name__ == '__main__':
    main()
