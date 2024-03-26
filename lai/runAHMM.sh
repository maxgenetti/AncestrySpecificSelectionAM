#!/bin/sh
# Maximilian Genetti (maxgenetti@gmail.com)
# Example command: runAHMM.sh 1 2 3 format.tsv


IFS=$'\t'
while IFS=$'\t' read -r key value; do
    IFS=$'\t' read -r -a "$key" <<< "$value"
done < "$4"

workDir=$(pwd)
dataDir="~/scratch/mgenetti/vcfamel45"
sub="final"
threads=30

###Deletes command files and then starts
rm -f panels.${sub}.txt
rm -f panels-2.${sub}.txt
rm -f ahmm1.${sub}.txt
for dist in "${dists[@]}"; do
	for mist in "${mists[@]}"; do
		for diff in "${diffs[@]}"; do
			for mthd in "${mthds[@]}"; do
				### defines and removes old directory
				outDirPR="/scratch/mgenetti/results${sub}/PR_${dist}_${mist}_${diff}_${mthd}"
				outDirMA="/scratch/mgenetti/results${sub}/MA_${dist}_${mist}_${diff}_${mthd}"
				if [ "$1" = "1" ]; then
					#rm -r -f $outDirPR
					mkdir -p "${outDirPR}"
					mkdir -p "${outDirPR}/ahmm1"
					#rm -r -f $outDirMA
					mkdir -p "${outDirMA}"
					mkdir -p "${outDirMA}/ahmm1"
				fi
				if [ $mthd = "g" ]; then
					g="-g "
					geno="-geno "
				else
					g=" "
					geno=" "
				fi
				for chrom in "${chrm1[@]}"; do
					mkdir -p "${outDirPR}/ahmm1/chr${chrom}"
					mkdir -p "${outDirMA}/ahmm1/chr${chrom}"
					### makes panel files
					echo "vcf2aHMM.py --vcf ${dataDir}/chr${chrom}.vcf --dist ${dist} --pop ${dataDir}/ahmm.pop --distM ${mist} --map ${dataDir}/apis.map ${geno}--minGTP 25 --minGTA 0.8 --minDP 0.8 --minDif $diff 1>${outDirPR}/chr${chrom}.all.panel 2>$outDirPR/chr${chrom}.panel.err " >>panels.${sub}.txt
					echo "awk '{for (i=1; i<=67; i++) printf \$i \"\t\"; print \"\"}' ${outDirPR}/chr${chrom}.all.panel >${outDirPR}/chr${chrom}.panel" >>panels-2.${sub}.txt
					echo "awk '{for (i=1; i<=7; i++) printf \$i \"	\"; for (i=68; i<=NF; i++) printf \$i \"	\"; print \"\"}' ${outDirPR}/chr${chrom}.all.panel >${outDirMA}/chr${chrom}.panel" >>panels-2.${sub}.txt
					### makes panel files with swapped parentals
					echo "awk 'BEGIN{OFS=FS=\"	\"} { temp = \$3; \$3 = \$5; \$5 = temp; temp = \$4; \$4 = \$6; \$6 = temp; print }' ${outDirPR}/chr${chrom}.panel 1>${outDirPR}/chr${chrom}-2.panel" >>panels-3.${sub}.txt
					echo "awk 'BEGIN{OFS=FS=\"	\"} { temp = \$3; \$3 = \$5; \$5 = temp; temp = \$4; \$4 = \$6; \$6 = temp; print }' ${outDirMA}/chr${chrom}.panel 1>${outDirMA}/chr${chrom}-2.panel" >>panels-3.${sub}.txt
					### runs ancestry_hmm
					echo "cd ${outDirPR}/ahmm1/chr${chrom} && ancestry_hmm -i ${outDirPR}/chr${chrom}.panel -s ${dataDir}/PR.sample -e 1e-3 --ne 100000 ${g}-a 2 0.5 0.5 -p 0 10000 0.5 -p 1 -100 0.5 1>chr.out 2>chr.err" >>ahmm1.${sub}.txt
					echo "cd ${outDirMA}/ahmm1/chr${chrom} && ancestry_hmm -i ${outDirMA}/chr${chrom}.panel -s ${dataDir}/MA.sample -e 1e-3 --ne 100000 ${g}-a 2 0.5 0.5 -p 0 10000 0.5 -p 1 -100 0.5 1>chr.out 2>chr.err" >>ahmm1.${sub}.txt
				done
			done
		done
	done
done
wait

###Executes if arg added, defaults to DRY RUN
if [ "$1" = "1" ]; then
	parallel --joblog panels.${sub}.log -j $threads <panels.${sub}.txt 2>panels.${sub}.err
	wait
	cd $workDir
	parallel --joblog panels-2.${sub}.log -j $threads <panels-2.${sub}.txt 2>panels-2.${sub}.err
	wait
	cd $workDir
	parallel --joblog panels-3.${sub}.log -j $threads <panels-3.${sub}.txt 2>panels-3.${sub}.err
	wait
	echo 'Finished Making panels'
	cd $workDir
	parallel --joblog ahmm1.${sub}.log -j $threads <ahmm1.${sub}.txt 2>ahmm1.${sub}.err
	wait
	echo 'Finished ahmm1.${sub}.txt'
else
	echo "
	parallel --joblog panels.${sub}.log -j $threads <panels.${sub}.txt 2>panels.${sub}.err
	wait
	cd $workDir
	parallel --joblog panels-2.${sub}.log -j $threads <panels-2.${sub}.txt 2>panels-2.${sub}.err
	wait
	echo 'Finished Making panels'
	cd $workDir
	parallel --joblog ahmm1.${sub}.log -j $threads <ahmm1.${sub}.txt 2>ahmm1.${sub}.err
	wait
	echo 'Finished ahmm1.${sub}.txt'
	"
fi

#second run of ahmm
rm -f ahmm2.${sub}.txt
for dist in "${dists[@]}"; do
	for mist in "${mists[@]}"; do
		for diff in "${diffs[@]}"; do
			for mthd in "${mthds[@]}"; do
				outDirPR="/scratch/mgenetti/results${sub}/PR_${dist}_${mist}_${diff}_${mthd}"
				outDirMA="/scratch/mgenetti/results${sub}/MA_${dist}_${mist}_${diff}_${mthd}"
				mkdir -p "${outDirPR}/ahmm2"
				mkdir -p "${outDirMA}/ahmm2"
				if [ $mthd = "g" ]; then
					g="-g "
					geno="-geno "
				else
					g=" "
					geno=" "
				fi
				if [ "$2" == "2" ]; then
					pP=$(tail -qn +2 ${outDirPR}/ahmm1/chr*/*.posterior | awk '{sum += $3}END{printf "%.3f", sum/NR}')
					qP=$(awk "BEGIN {print 1 - ${pP}}")
					pM=$(tail -qn +2 ${outDirMA}/ahmm1/chr*/*.posterior | awk '{sum += $3}END{printf "%.3f", sum/NR}')
					qM=$(awk "BEGIN {print 1 - ${pM}}")
				fi
				for chrom in "${chrm1[@]}"; do
					mkdir -p "${outDirPR}/ahmm2/chr${chrom}"
					mkdir -p "${outDirMA}/ahmm2/chr${chrom}"
					echo "cd ${outDirPR}/ahmm2/chr${chrom} && ancestry_hmm -i ${outDirPR}/chr${chrom}.panel -s ${dataDir}/PR.sample -e 1e-3 --ne 100000 -a 2 ${pP} ${qP} -p 0 10000 ${pP} -p 1 -100 ${qP} 1>chr.out 2>chr.err" >>ahmm2.${sub}.txt
					echo "cd ${outDirMA}/ahmm2/chr${chrom} && ancestry_hmm -i ${outDirMA}/chr${chrom}.panel -s ${dataDir}/MA.sample -e 1e-3 --ne 100000 -a 2 ${pM} ${qM} -p 0 10000 ${pM} -p 1 -100 ${qM} 1>chr.out 2>chr.err" >>ahmm2.${sub}.txt
				done
			done
		done
	done
done

###Executes if arg added, defaults to DRY RUN
if [ "$2" = "2" ]; then
	cd $workDir
	parallel --joblog ahmm2.${sub}.log -j $threads <ahmm2.${sub}.txt 2>ahmm2.${sub}.err
	wait
	echo 'Finished ahmm2.${sub}.txt'
else
	echo "
	cd $workDir
	parallel --joblog ahmm2.${sub}.log -j $threads <ahmm2.${sub}.txt 2>ahmm2.${sub}.err
	wait
	echo 'Finished ahmm2.${sub}.txt'
	"
fi

#creates ahmms.${sub}.txt
rm -f ahmms.${sub}.txt
rm -f results.${sub}.txt
for demo in "${demos[@]}"; do
	for dist in "${dists[@]}"; do
		for mist in "${mists[@]}"; do
			for diff in "${diffs[@]}"; do
				for mthd in "${mthds[@]}"; do
					if [ $mthd = "g" ]; then
						g="-g "
					else
						g=" "
					fi
					outDir="/scratch/mgenetti/results${sub}/${demo}_${dist}_${mist}_${diff}_${mthd}"
					if [ "$3" = "3" ]; then
						parallel -j 12 --keep-order cat ${outDir}/ahmm2/chr{}/*posterior \| marginals.py 1 {} \>${outDir}/marginals-chr{}.txt ::: ${chrm2[@]} 2>${outDir}/summary-marginals.txt
						p=$(tail -qn +2 ${outDir}/ahmm2/chr*/*.posterior | awk '{sum += $3}END{printf "%.3f", sum/NR}')
						q=$(awk "BEGIN {print 1 - ${p}}")
						s=$(tail -qn +2 ${outDir}/ahmm2/chr*/*.posterior | cut -f 3 | awk '{delta = $1; sum += delta; sum_sq += delta^2} END {mean = sum/NR; variance = sum_sq/NR - mean^2; print sqrt(variance)}')
						g=$(tail -qn 1 ${outDir}/ahmm2/chr*/chr.out | cut -f 3 | awk '{sum += $1}END{print int(sum/NR)}')
						z=$(tail -qn 1 ${outDir}/ahmm2/chr*/chr.out | cut -f 3 | awk '{delta = $1; sum += delta; sum_sq += delta^2} END {mean = sum/NR; variance = sum_sq/NR - mean^2; print sqrt(variance)}')
					fi

					echo "chroms	${chrm2[@]}" >${outDir}/params.txt
					echo "generation	${g}" >>${outDir}/params.txt
					echo "pulse	${p}" >>${outDir}/params.txt
					echo "steps	1" >>${outDir}/params.txt
					#echo "pusle_SD	${s}" >>${outDir}/params.txt
					#echo "gen_SD	${z}" >>${outDir}/params.txt
					for morg in "${morgs[@]}"; do
						outDir2="/scratch/mgenetti/results${sub}/${demo}_${dist}_${mist}_${diff}_${mthd}/${morg//[ ]/_}"
						mkdir -p "${outDir2}"
						echo "${outDir2}" >>results.${sub}.txt

						for chrom in "${chrm2[@]}"; do
							if [ "$3" = "3" ]; then
								pos=$(tail -n 1 ${outDir}/chr${chrom}.panel | cut -f 1-2)
								echo "${pos}" >>${outDir}/contigs.txt
								#p=$(tail -qn +2 ${outDir}/ahmm2/chr${chrom}/*.posterior | awk '{sum += $3}END{print sum/NR}')
								#q=$(awk "BEGIN {print 1 - ${p}}")
							fi
							echo "ahmm-s -i ${outDir}/chr${chrom}.panel -s ${dataDir}/${demo}.sample -p 1 10000 ${p} -p 0 $g ${q} --window ${morg} --gss 2 980000 1 0.001 0.5 --ne 100000 1>${outDir2}/ahmms-forw-chr${chrom}.out 2>${outDir2}/ahmms-forw-chr${chrom}.err" >>${workDir}/ahmms.${sub}.txt
							echo "ahmm-s -i ${outDir}/chr${chrom}-2.panel -s ${dataDir}/${demo}.sample -p 1 10000 ${q} -p 0 $g ${p} --window ${morg} --gss 2 980000 1 0.001 0.5 --ne 100000 1>${outDir2}/ahmms-rev-chr${chrom}.out 2>${outDir2}/ahmms-rev-chr${chrom}.err" >>${workDir}/ahmms.${sub}.txt
						done
					done
				done
			done
		done
	done
done
wait

###Executes ahmms if arg added, defaults to DRY RUN
if [ "$3" = "4" ]; then
	cd $workDir
	conda activate ahmms
	parallel --joblog ahmms.${sub}.log -j $threads <ahmms.${sub}.txt 2>ahmms.${sub}.err
	wait
	echo 'DONE'
	conda deactivate
else
	echo "
	cd $workDir
	conda activate ahmms
	parallel --joblog ahmms.${sub}.log -j $threads <ahmms.${sub}.txt 2>ahmms.${sub}.err
	wait
	echo 'DONE'
	conda deactivate
	"
fi
