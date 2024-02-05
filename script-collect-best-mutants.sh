#!/bin/bash

f="/project/olga-phylo/LPOR/2017.06.June/revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod"
f_in="results-collection"
f_out="results-REvolver-sim-1/best_mutant_per_SDR"

mkdir -p $f/$f_out
file="$f/TIGRFAM-best-mutants"

cat $file | while read -r line
do
	id="`echo "$line" | awk -F " " '{print $1}'`"
	mutant="`echo "$line" | awk -F " " '{print $3}'`"
	sp="`grep "root.$id-" $f/plots-joint/input/table-names-root-ids.txt | awk -F " " '{print $2}'`"
	folder="`ls -l $f/$f_in/* | grep "root.$id-" | awk -F " " '{print $9}'`"
	genus="`echo "$sp" | awk -F "_" '{print $1}'`"
	file_seq="$f/$f_in/*$genus*/$folder/sequences.$id.TIGRFAM.all"

	#echo "====================================================="
	#echo "genus: $genus"
	#echo "root.$id"
	#echo "species: $sp"
	#echo "mutant: $mutant"
	#echo "file: $file_seq"
	echo "`grep -A 1 "$mutant" $file_seq | sed 's/'$mutant'/'$sp-$mutant'/g' `" >> best_TIGRFAM_mutants_all.fa
done
