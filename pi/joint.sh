#!/bin/bash
#SBATCH --job-name=joint
#SBATCH --output=joint.out
#SBATCH --error=joint.err
#SBATCH --ntasks=40
#SBATCH --mem 240GB
#SBATCH --oversubscribe
#SBATCH --account=256x20
#SBATCH --qos=256x20
#SBATCH --partition=256x20
#SBATCH -t 5-00:00:00

#/hb/home/mgenetti/bin/angsd/angsd -GL 1 -b pr.list -anc /hb/scratch/mgenetti/drones/reference/GCF_000002195.4_Amel_4.5_genomic.fna -P 20 -out pr.haploid -doSaf 1 -minMapQ 30 -minQ 20 -isHap 1 -sites inclusion_regions.bed &
#/hb/home/mgenetti/bin/angsd/angsd -GL 1 -b ma.list -anc /hb/scratch/mgenetti/drones/reference/GCF_000002195.4_Amel_4.5_genomic.fna -P 20 -out ma.haploid -doSaf 1 -minMapQ 30 -minQ 20 -isHap 1 -sites inclusion_regions.bed &
#wait
parallel -j 40 <joint.txt
