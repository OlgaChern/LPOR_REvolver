#!/bin/bash
f="/project/olga-phylo/LPOR/2017.06.June/revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod/results-collection"

for i in $f/*/root.*/
do
	id="`echo "$i" | awk -F "root." '{print $2}' | awk -F "-" '{print $1}'`"
	score="`head -n 1 $i/results.TIGRFAM.hmm.scores.ordered | awk -F " " '{print $9}'`"
	mutant="`head -n 1 $i/results.TIGRFAM.hmm.scores.ordered | awk -F " " '{print $1}'`"
	echo "$id $score $mutant"
done
