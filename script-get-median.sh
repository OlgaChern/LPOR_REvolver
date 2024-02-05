#!/bin/bash


f_main="/project/olga-phylo/LPOR/2017.03.March-final-dataset/revolver-analysis"
f_out_end="revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod"
f_out="/project/olga-phylo/LPOR/2017.06.June/$f_out_end"
#f_scratch=$f_out
f_scratch="/scratch/olga/LPOR/2017.06.June/$f_out_end"

id=56          # root-sequence id
rev_test=6     # revolver model settings
f_res_revolver="results/revolver.$rev_test/root.${id}"


threshold_history="$f_out/$f_res_revolver/threshold.history"

for s in {1..502}
do
	gen_hmm_scores="$f_scratch/$f_res_revolver/iteration.$s/generation.hmm.scores"
        gen_median="$(~/Dropbox/R.codes/LPOR-selection-threshold-2.r $gen_hmm_scores)"
        echo "GENERATION $s | generation median reached value of $gen_median" >> $threshold_history
done

