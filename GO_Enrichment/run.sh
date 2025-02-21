for G in {1..10}; do
    for P in PR MA; do
        for D in rev forw; do
            python3 simulate.py -g ${G} -p ${P} -s 100000 -i ${D} -d 50 2>./pickles/${P}-${D}.bio.${G}.window.err
	        python3 simulate.py -g ${G} -p ${P} -s 100000 -i ${D} -o -d 50 2>./pickles/${P}-${D}.mol.${G}.window.err
        done
    done
done

for P in MA PR; do
    for D in rev forw; do
        python3 checkRealData.py --pop ${P} --d ${D} --go &
        python3 checkRealData.py --pop ${P} --d ${D} &	
    done
done
wait

rm -f err.check
for P in PR MA; do
    for D in rev forw; do
        for W in 50; do
	    python3 enrichment.py --pop ${P} --d ${D} --sim 100000 --dup 10 --window ${W} --go 2>>err.check
	    python3 enrichment.py --pop ${P} --d ${D} --sim 100000 --dup 10 --window ${W} 2>>err.check
        done
    done
done

rm -f err.mol
rm -f err.bio
rm -f log.mol
rm -f log.bio
for filename in ./mol/*.tsv; do
	fname=${filename::-4}
	python3 fdr.py $fname 2>>err.mol 1>>log.mol&
done
wait
for filename in ./bio/*.tsv; do
	fname=${filename::-4}
	python3 fdr.py $fname 2>>err.bio 1>>log.bio&
done
wait
