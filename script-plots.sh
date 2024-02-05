#!/bin/bash


f_main="/project/olga-phylo/LPOR/2017.03.March-final-dataset/revolver-analysis/input"
f="/project/olga-phylo/LPOR/2017.06.June/revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod"
f_min_max="$f/results/revolver.6"
f_in="$f/results-collection"

for folder in $f_in/*/root.*/
#for folder in /project/olga-phylo/LPOR/2017.06.June/revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod/results-collection/Agrobacterium-Rhizobium/root.11-WP_025418312.1_Rhizobium_leguminosarum 
#for folder in $f/results-collection/Caulobacter/root.35-WP_013079543.1_Caulobacter_segnis
#for folder in $f/results-collection/Agrobacterium-Rhizobium/root.13-WP_069610514.1_Rhizobium_sp._YK2
do
	id="`echo "$folder" | awk -F "root." '{print $2}' | awk -F "-" '{print $1}'`"
	echo "========================================"
	echo "Plotting figures for root $id.."
	echo "========================================"
	tigrfam_scores="$folder/results.TIGRFAM.hmm.scores"
	median="$folder/threshold.history.$id"
	min="$f_min_max/root.$id/generation_min.$id"
	max="$f_min_max/root.$id/generation_max.$id"
	gen_proteo="`awk -F " " '{print $2}' $folder/generation.cutoff.1e-186.proteo`"
        gen_cyano="`awk -F " " '{print $2}' $folder/generation.cutoff.1e-157.cyano | tail -n 1`"
	~/Dropbox/R.codes/LPOR-plots-median-min.r $tigrfam_scores $min $median $max $gen_cyano $gen_proteo


	#cp $f_main/root-seq-$id.fa $folder/initial-root-seq-$id.fa

done
