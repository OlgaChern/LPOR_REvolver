#!/bin/bash
#start=`date +%s`
#=====================================================
#	Revolver binaries
#=====================================================
REvolver="/project/olga-phylo/LPOR/2016.12.December-REvolver/REvolver.jar"
hmm_binaries="/project/olga-phylo/LPOR/2016.12.December-REvolver/hmmer-3.1b2-linux-intel-x86_64/binaries"
hmmfetch_location="$hmm_binaries/hmmfetch"
hmmscan_location="$hmm_binaries/hmmscan"
#=====================================================
#	Working directories
#=====================================================
f_main="/project/olga-phylo/LPOR/2017.03.March-final-dataset/revolver-analysis"
f_out_end="revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod"
f_out="/project/olga-phylo/LPOR/2017.06.June/$f_out_end"
#f_scratch=$f_out
f_scratch="/scratch/olga/LPOR/2017.06.June/$f_out_end"
mkdir -p $f_scratch
#=====================================================
#	Parameters
#=====================================================
id=$1		# root-sequence id
sim_test=$2	# sim run settings
rev_test=$3	# revolver model settings
#----------------------------------
#	Check
#----------------------------------
if [ "$id" == "" ]
then
	echo "ERROR: You need to specify id of root-sequence!"
	echo "Exiting.."
	exit 1
fi
if [ "$sim_test" == "" ]
then
	echo "ERROR: You need to specify parameters for simulations!" 
	echo "Check possible sim_search_test values in sim-parameters."
	echo "Exiting.."
	exit 1
fi
if [ "$rev_test" == "" ]
then
	echo "ERROR: You need to specify parameters for revolver model!" 
	echo "Check possible rev_test values in sim-parameters."
	echo "Exiting.."
	exit 1
fi
#----------------------------------
# Reading file with parameters
#----------------------------------
source sim-parameters-1.sh $sim_test $rev_test
#----------------------------------
#       User datasets
#----------------------------------
d_hmm_folder="$f_main/databases/database.hmm"
d_hmm_db="aln-70_oxi-1SDRactino.fa.hmm"

root_file_orgn="$f_main/input/root-seq-${id}.fa"
query_orgn=$root_file_orgn
#----------------------------------
#	Results folders and files
#----------------------------------
f_main_hmmer="$f_out/results/hmmer.$rev_test"
f_main_revolver="$f_out/results/revolver.$rev_test"
[ -d $f_main_hmmer ] || mkdir -p $f_main_hmmer
[ -d $f_main_revolver ] || mkdir -p $f_main_revolver

f_res_hmmer="results/hmmer.$rev_test/root.${id}"
mkdir -p $f_scratch/$f_res_hmmer
f_res_revolver="results/revolver.$rev_test/root.${id}"
mkdir -p $f_scratch/$f_res_revolver

######################################################
# MAIN simulations
######################################################
if [ $run_sim == "yes" ]
then
	#----------------------------------------------------------
        # Get the score and the length of the root sequence
	#----------------------------------------------------------
        if [ $run_sim_search_hmmer_orgn == "yes" ]
        then

	hmm_score_comparison=()
        for d_hmm in $d_hmm_db # here d_hmm_db contains SDR hmm model [1] and LPOR hmm model [2]
        do
                #printf "\n-----------------------------------------------------------------\n"
                #printf "\nRunning HMMER search for original root sequence\n - ${query_orgn}\n - $d_hmm...\n"
                results_orgn="$f_scratch/$f_res_hmmer/original_sequence/$d_hmm/threshold.none/root.${id}.hmmer"
                mkdir -p "$f_scratch/$f_res_hmmer/original_sequence/$d_hmm/threshold.none"

                $hmm_binaries/hmmscan $d_hmm_folder/${d_hmm} ${query_orgn} > $results_orgn

                # Get sequence length   # wc -L - prints the length of the longest line
                len=`tail -n +2 $query_orgn | tr --delete '\n' | wc -L | awk -F " " '{print $1}'`
                #echo "the length is $len..."

                # Get the score of original sequence
                score_hmmer="$f_scratch/$f_res_hmmer/original_sequence/$d_hmm/threshold.none/results.root.${id}"
                n=`tail -n 5 $results_orgn | head -n 1 | awk -F " " '{print $5}'`
                if [ $n -gt 0 ]
                then
			# 16 because root sequences have one additional line "description of sequences", while evolved sequences don't have it
                	echo "${id} $len `head -n 16 $results_orgn | tail -n 1`" | cat >> $score_hmmer
                        #echo "the score is `head -n 16 $results_orgn | tail -n 1`"
			x="`head -n 16 $results_orgn | tail -n 1 | awk -F " " '{print $1}'`"
                else
                        echo "${id} $len  1 1 1 1 1 1 1 1 1" | cat >> $score_hmmer
                        #echo "no hit..."
			x="666"
                fi
		hmm_score_comparison=(${hmm_score_comparison[*]} $x)
	done
        fi
	
	###########################################################################
	# Start with original sequence
	###########################################################################
	root_file="$root_file_orgn"
	###########################################################################
	# Generations: one step is one evolutionary process on star-tree
	###########################################################################

	threshold_history="$f_scratch/$f_res_revolver/threshold.history"
	gen_final="`tail -n 1 $threshold_history | awk -F " " '{print $2}'`"
	all_it="$(seq 1 $gen_final)"
	for s in $all_it
	do
		#printf "\n===== GENERATION $s =============================================\n"
		f_seq_evolve="$f_scratch/$f_res_revolver/iteration.$s/roots/"
#####################################################################################################
	for root_file in $f_seq_evolve/*
	do
		root_new_id="`echo "$root_file" | awk -F "root_new." '{print $2}'`"
        	for i in $search_revolver
        	do
                        output_path="$f_scratch/$f_res_revolver/iteration.$s/$root_new_id/revolver.run.$i"
			############################################################################################
                        # (2) analysis of evolved sequence with HMMER
                        ############################################################################################
			query_fa=$output_path/out.fa
			hmm_score_comparison=()
			for d_hmm in $d_hmm_db
                        do
				results_h="$f_scratch/$f_res_hmmer/iteration.$s/$d_hmm/threshold.none/$root_new_id/root.${id}.$i.$s.hmmer"
				mkdir -p "$f_scratch/$f_res_hmmer/iteration.$s/$d_hmm/threshold.none/$root_new_id"
				$hmm_binaries/hmmscan $d_hmm_folder/${d_hmm} ${query_fa} > $results_h
				
				# Get sequence length   # wc -L - prints the length of the longest line
                                len=`tail -n +2 ${query_fa} | tr --delete '\n' | wc -L | awk -F " " '{print $1}'`
				#----------------------------------------------------------------------------------
				# Get the score of MANY evolved sequences - use this code for star tree
				#----------------------------------------------------------------------------------
                                #score_hmmer="$f_scratch/$f_res_hmmer/iteration.$s/$d_hmm/threshold.none/$root_new_id/all.root.${id}.$i"
				#score_hmmer_1="$f_scratch/$f_res_hmmer/iteration.$s/$d_hmm/threshold.none/$root_new_id/all.root.${id}.$i.1"
				#echo "`grep  -A 2 "Model" $results_h | sed "/Model/d" | sed "/------- ------ -----    ------- ------ -----   ---- --  --------   -----------/d" | sed "/--/d"`" >> $score_hmmer
				#awk -F " " '{print $1}' $score_hmmer >> $score_hmmer_1

				#----------------------------------------------------------------------------------
				# Get the score of ONE evolved sequence - use this code for one branch tree
				#----------------------------------------------------------------------------------
                		score_hmmer="$f_scratch/$f_res_hmmer/iteration.$s/$d_hmm/threshold.none/$root_new_id/all.root.${id}.$i"
                		n=`tail -n 5 $results_h | head -n 1 | awk -F " " '{print $5}'`
                		if [ $n -gt 0 ]
                		then
                        		echo "${id} $len `head -n 15 $results_h | tail -n 1`" | cat >> $score_hmmer
                        		x="`head -n 15 $results_h | tail -n 1 | awk -F " " '{print $1}'`"
                		else
                        		echo "${id} $len  1 1 1 1 1 1 1 1 1" | cat >> $score_hmmer
                        		x="666"
                		fi
                		hmm_score_comparison=(${hmm_score_comparison[*]} $x)
				gen_hmm_scores="$f_scratch/$f_res_hmmer/iteration.$s/$d_hmm/threshold.none/generation.hmm.scores.it_$s"
				echo "$root_new_id $s 1 ${hmm_score_comparison[*]}" >> $gen_hmm_scores
			done

		done	# lineages
		done	# roots_new
		done	# iterations

fi

#cp -r $f_scratch/$f_res_hmmer $f_main_hmmer
#cp -r $f_scratch/$f_res_revolver $f_main_revolver
#rm -r $f_scratch/$f_res_hmmer
#rm -r $f_scratch/$f_res_revolver
#end=`date +%s`

#echo "Time: $((end-start))"
