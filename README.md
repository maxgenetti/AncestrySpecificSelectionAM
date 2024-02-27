Directory PR contains the sample and panel files used to run Ancstry_HMM and AHMM-S on the Puerto Rican population of Africanized honey bee. The marginals are the population average of the Ancestry_HMM results. The ahmms-forw files are the AHMMS results for selection on African ancestry. The ahmms-rev files are the AHMMS results for selection on European ancestry.

Directory MA contains the sample and panel files used to run Ancstry_HMM and AHMM-S on the Mesoamerican population of Africanized honey bee. The marginals are the population average of the Ancestry_HMM results. The ahmms-forw files are the AHMMS results for selection on African ancestry. The ahmms-rev files are the AHMMS results for selection on European ancestry.


Example of how to generate PANEL files using vcf2aHMM.py :

python vcf2aHMM.py --vcf chr1.vcf --dist 2500 --pop ahmm.pop --distM 1e-4 --map apis.map -geno --minGTP 25 --minGTA 0.8 --minDif 0.10 1>chr1.panel 2>chr1.err


Example of how to run Ancestry_HMM :

ancestry_hmm -i chr1.panel -s PR.sample -e 1e-3 --ne 100000 -g -a 2 0.5 0.5 -p 0 10000 0.5 -p 1 -100 0.5 1>chr1.out 2>chr1.err


Example of how to run AHMM-S :

ahmm-s -i chr1.panel -g -s PR.sample -p 1 10000 0.5 -p 0 100 0.5 --window p 100 --gss 2 980000 1 0.001 0.5 --ne 100000 1>ahmms-forw-chr1.out 2>ahmms-forw-chr1.err


peakFinder.py can simply be run after specifying the populations, directions, and location of the data.
