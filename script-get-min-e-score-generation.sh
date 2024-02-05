#!/bin/bash


#f_main="/project/olga-phylo/LPOR/2017.03.March-final-dataset/revolver-analysis"
f_out_end="revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod"
f_out="/project/olga-phylo/LPOR/2017.06.June/$f_out_end"
f_scratch=$f_out
#f_scratch="/scratch/olga/LPOR/2017.06.June/$f_out_end"
# for root 11
#f_scratch="$f_out/tar-all-results/revolver"

# project
#root_ids="10 19 2 25 27 33 34 36 38 51 53 56 59 63 67 68 71 73 8 80 83" #"1"
#root_ids="4"

# ritchie
#root_ids="48" #"14 15 16 21 24 28 29 30 37 39 40 41 44 48 49 50 58 6 62 65 66 72 74 75 76 79"

# fitch
#root_ids="78" #"31 47 54 78 82"

# zuse
#root_ids="12 18 32 35 43 5 52 61 69"
#root_ids="81"

# In separate locations
#root_ids="20" # franklin:scratch

# In tar-all
#root_ids="11"


# Filling gaps;) (from done-zuse.tar)
#root_ids="1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 24 25 26 27 28 29 30 31 33 34 35 36 37 38 39 4 40 41 42 43 44 45 46 47 48 49 5 50 51 53 54 56 59 6 61 62 63 64 65 66 67 68 7 71 72 73 74 75 78 79 8 80 81 82 83 9"

root_ids="13"


for id in $root_ids
do
	rev_test=6     # revolver model settings
	f_res_revolver="results/revolver.$rev_test/root.${id}"
	generation_min="$f_out/$f_res_revolver/generation_min.$id"
	mkdir -p $f_out/$f_res_revolver/
	threshold_history="$f_scratch/$f_res_revolver/threshold.history"
	# for root 11
	#threshold_history="$f_scratch/root.${id}/threshold.history"
	if [ -e $threshold_history ]
	then
		gen_final="`tail -n 1 $threshold_history | awk -F " " '{print $2}'`"
	else
		gen_final="`ls $f_res_revolver | grep "iteration" | wc -l | awk -F " " '{print $1}'`"
	fi
	echo "======================================="
	echo "Results for root $id..."
	echo "======================================="
	if [ -e $generation_min ]
	then
		echo "already done.."
	else
		echo "final generation $gen_final..."
		for s in $(seq 1 $gen_final)
		do
			gen_hmm_scores="$f_scratch/$f_res_revolver/iteration.$s/generation.hmm.scores"
			# for root 11
			#gen_hmm_scores="$f_scratch/root.${id}/iteration.$s/generation.hmm.scores" 
        		gen_min="$(~/Dropbox/R.codes/LPOR-min-generation.r $gen_hmm_scores)"
        		echo "GENERATION $s | generation min reached value of $gen_min" >> $generation_min
			#echo "GENERATION $s | generation min reached value of $gen_min"
		done
	fi
done

