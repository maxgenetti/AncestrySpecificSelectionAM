import os
from collections import defaultdict
import argparse

def read_to_list(filename):
    with open(filename, 'r') as file:
        lines = [line.rstrip().split()[0] for line in file.readlines()]
    return lines

def read_to_dict(filename):
    with open(filename, 'r') as file:
        fileDict = defaultdict(list)
        for line in file :
            fileDict[line.rstrip().split()[0]].append(line.rstrip().split()[1],line.rstrip().split()[2])
    return fileDict

def writeCommands(samples, chroms, first, ref, outdir, tmp, minQ) :
    x = first
    for sample,fastqs in samples.items():
        r1 = " ".join([item[0] for item in fastqs])
        r2 = " ".join([item[1] for item in fastqs])
        with open(y,"w") as f :
            f.write(f"fastp -i <(zcat {r1}) -I <(zcat {r2}) -o {tmp}/{sample}.R1.fq.gz -O {tmp}/{sample}.R2.fq.gz --detect_adapter_for_pe --unpaired1 {tmp}/{sample}.unpaired.fq.gz --unpaired2 {tmp}/{sample}.unpaired.fq.gz \n")
            f.write(f"bwa mem -t {threads} -R '@RG\\tID:{sample}\\tSM:{sample}' {ref}  {tmp}/{sample}.unpaired.fq.gz | samtools sort -o {tmp}/{sample}.se.bam -@ {threads} \n")
            f.write(f"bwa mem -t {threads} -R '@RG\\tID:{sample}\\tSM:{sample}' {ref}  {tmp}/{sample}.R1.fq.gz {tmp}/{sample}.R2.fq.gz | samtools sort -o {tmp}/{sample}.pe.bam -@ {threads} \n")
            f.write(f"samtools cat {tmp}/{sample}.pe.bam {tmp}/{sample}.se.bam -@ {threads} | samtools sort -@ {threads} | samtools fixmate - -  -@ {threads} | sambamba markdup -r -t {threads} - {outdir}/{sample}.raw.bam \n")
            f.write(f"samtools index {outdir}/{sample}.raw.bam -Q {threads} \n")
            f.write(f"samtools view {outdir}/{sample}.raw.bam -L {chroms} -q {minQ} -o {outdir}/{sample}.filtered.bam \n")
            f.write(f"samtools index {outdir}/{sample}.filtered.bam -Q {threads} \n")
        x += 1

def main():
    parser = argparse.ArgumentParser(description='Read lines from a file into a list.')
    parser.add_argument('--fastq', type=str, help='file with path to fastq files, EXMAPLE:sample_name A_R1.fastq.gz A_R2.fastq.gz')
    parser.add_argument('--format', type=str, help='file with paths to reference and output directory')
    parser.add_argument('--first', type=int, default = 1, help='starting number for job array')
    parser.add_argument('--chroms', type=str, help='Chromsomes to include in filtered bam, bed format')
    parser.add_argument('--threads', type=int,default = 1, help='Chromsomes to include in VCF')  
    parser.add_argument('--minQ', type=int,default = 20, help='Chromsomes to include in VCF')  

    args = parser.parse_args()
    ref, outdir = read_to_list(args.format)
    samples = read_to_dict(args.samples)
    writeCommands(samples, chroms, args.first, args.threads, ref, outdir, args.minQ)

if __name__ == '__main__':
    main()


import os
from collections import defaultdict


    print('fastp -i <(zcat %s) -I <(zcat %s) -o out.R1.fq.gz -O out.R2.fq.gz --detect_adapter_for_pe --unpaired1 unpaired.fq.gz --unpaired2 unpaired.fq.gz' %(" ".join(locations1[sample]), " ".join(locations2[sample])))
    print("bwa mem -t 8 -R "+"'"+"@RG\\tID:%s\\tSM:%s' /scratch1/mgenetti/meowMix/reference/Oge-1_final.fasta.gz out.R1.fq.gz out.R2.fq.gz 2>./err/%s.err | samtools sort -o paired.bam" %(sample, sample,sample))
    print("bwa mem -t 8 -R "+"'"+"@RG\\tID:%s\\tSM:%s' /scratch1/mgenetti/meowMix/reference/Oge-1_final.fasta.gz unpaired.fq.gz 2>./err/%s.se.err | samtools sort -o unpaired.bam" %(sample, sample, sample))

    print("samtools cat paired.bam unpaired.bam | samtools sort --threads 8 | samtools fixmate - - | sambamba markdup -r -t 8 - ../bwa/%s.bam" %(sample))