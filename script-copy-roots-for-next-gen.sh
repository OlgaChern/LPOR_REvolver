#!/bin/bash
f_out_end="revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod"
f_out="/project/olga-phylo/LPOR/2017.06.June/$f_out_end"
f_scratch="/scratch/olga/LPOR/2017.06.June/$f_out_end"
#=====================================================
#	Parameters
#=====================================================
id=20		# root-sequence id
rev_test=6	# revolver model settings
f_res_revolver="results/revolver.$rev_test/root.${id}"

s=704 # iteration from which sample parents will be copied to iteration.s+1/roots
s_1="`echo $s+1 | bc`"
f_seq_next="$f_scratch/$f_res_revolver/iteration.$s_1/roots/"		
gen_hmm_scores="$f_scratch/$f_res_revolver/iteration.$s/generation.hmm.scores"

##################################################################
# Collect top k sequences for the next generation
##################################################################
		top_k_seq="$f_scratch/$f_res_revolver/iteration.$s/sample-next-parents"
                #~/Dropbox/R.codes/LPOR-selection-probabilistic-mod.r $gen_hmm_scores $top_k_seq	
                check_not_NA="`head -n 1 $top_k_seq | awk  -F " " '{print $3}'`"
                if [ "$check_not_NA" != "NA" ]
                then
			next_id="`ls $f_seq_next | wc -l | awk -F " " '{print $1}'`"
			cat $top_k_seq | while read -r line
			do
				root_new_id="`echo "$line" | awk -F " " '{print $1}'`"	# root from previous generation, which was chosen to be a parent
				i="`echo "$line" | awk -F " " '{print $3}'`"		# there is only one sequence for one branch tree, i is always 1
				scores_aux="`echo "$line" | awk -F " " '{print $4 " " $5}'`"
				next_r_file="$f_scratch/$f_res_revolver/iteration.$s/$root_new_id/revolver.run.1/out.fa"

				next_id="`echo $next_id+1 | bc`"	

				#new_r_name="seq_${i}_"
				#echo "`grep -A 1 "$new_r_name" $next_r_file`" > $f_seq_next/root_new.$next_id
				cp $next_r_file $f_seq_next/root_new.$next_id

                        	next_folder="$f_scratch/$f_res_revolver/iteration.$s_1/$next_id"
                        	mkdir -p $next_folder
			
				run_history_file="$f_scratch/$f_res_revolver/iteration.$s/$root_new_id/run.history"
                        	cp $run_history_file $next_folder

                        	echo "$f_scratch/$f_res_revolver/iteration.$s_1/roots/root_new.$next_id $scores_aux" >> $next_folder/run.history
			done

                else
                       	echo "- No sequences were selected for the next generation roots"
			exit
                fi

