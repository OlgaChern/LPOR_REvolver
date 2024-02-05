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
#	Selection
#----------------------------------
#fitness_threshold="2.699404e-32"       # default initial fitness threshold - compute it as a mean score of all SDR sequences
#SDR_mean="$fitness_threshold"
#LPOR_mean="$fitness_threshold"
top_k=5
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
#d_hmm_db="aln-oxi-1SDRactino.fa.hmm hmm-LPOR-proteo-15.fasta"
d_hmm_db="aln-oxi-1SDRactino.fa.hmm hmm-LPOR-proteo-14-full"
#----------------------------------
#	Input files
#----------------------------------
hmmfile_sdr="$d_hmm_folder/aln-oxi-1SDRactino.fa.hmm"
#hmmfile_lpor="$d_hmm_folder/hmm-LPOR-proteo-15.fasta"


##############################################################################
#	For truncation    selection use a star tree with 50 branches
#	For probabilistic selection use a one branch tree
##############################################################################
tree="$f_main/input/tree-1-branch-0.005.tree"
#tree="/project/olga-phylo/LPOR/2017.03.March-final-dataset/revolver-all/evolution-domain-constraints/trees/tree-star-0.005"

root_file_orgn="$f_main/input/root-seq-${id}.fa"
#root_file_hmm_orgn_sdr="$f_main/input/root-seq-${id}.fa.sdr.hmm"
#root_file_hmm_orgn_lpor="$f_main/input/root-seq-${id}.fa.lpor.hmm"
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
        # Print revolver model details
        #----------------------------------------------------------
	model_info="$f_scratch/$f_res_revolver/model.info"
        echo "model_name $model_name"            > $model_info
        echo "insertion_rate $insertion_rate"   >> $model_info
        echo "deletion_rate $deletion_rate"     >> $model_info
        echo "len_distr_in $len_distr_in"       >> $model_info
        echo "len_distr_del $len_distr_del"     >> $model_info
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
        hmmfile=$hmmfile_sdr
        #hmmfile=$hmmfile_lpor

	# HMM annotations will be computed for root files later on
	# root_file_hmm="$root_file_hmm_orgn_sdr"
        # root_file_hmm="$root_file_hmm_orgn_lpor"

	f_seq_evolve="$f_scratch/$f_res_revolver/iteration.1/roots"
	mkdir -p $f_seq_evolve

	threshold_history="$f_scratch/$f_res_revolver/threshold.history"
	#next_threshold="$f_scratch/$f_res_revolver/iteration.1/compute-threshold"

	# You copy original sequence 5 times to make it consistent with the selection of 5 parents for the generation
	for id_k in {1..250}
        do
		cp $root_file $f_seq_evolve/root_new.$id_k
		run_history="$f_scratch/$f_res_revolver/iteration.1/$id_k/run.history"
		mkdir -p $f_scratch/$f_res_revolver/iteration.1/$id_k
		echo "$f_seq_evolve/root_new.$id_k ${hmm_score_comparison[*]}" >> $run_history
	done

	#for id_k in {1..5}
	#do
	#	echo "$f_scratch/$f_res_revolver/iteration.1/roots/root_new.$id_k ${hmm_score_comparison[*]}" >> $next_threshold
	#done
	###########################################################################
	# Generations: one step is one evolutionary process on star-tree
	###########################################################################
	for s in $iterations
	do
		#printf "\n===== GENERATION $s =============================================\n"
		f_seq_evolve="$f_scratch/$f_res_revolver/iteration.$s/roots/"
		s_1="`echo $s+1 | bc`"
		f_seq_next="$f_scratch/$f_res_revolver/iteration.$s_1/roots/"
		mkdir -p $f_seq_next
		gen_hmm_scores="$f_scratch/$f_res_revolver/iteration.$s/generation.hmm.scores"
		#printf "Number of sequences to evolve: `ls $f_seq_evolve | wc -l| awk -F " " '{print $1}'`\n"

		# Compute threshold------------------------------------------------------------------
		#next_threshold="$f_scratch/$f_res_revolver/iteration.$s/compute-threshold"
		#fitness_threshold="$(~/Dropbox/R.codes/LPOR-selection-threshold.r $next_threshold)"
		#echo "GENERATION $s | Fitness threshold is set to $fitness_threshold" >> $threshold_history


		####################################################################################
		# stop when fitness threshold reaches 1e-200 or after 1000 generations
		####################################################################################
		#stop_rule="1e-200"
		#check_fitness_threshold="$(~/Dropbox/R.codes/LPOR-compare-2-numbers.r $fitness_threshold $stop_rule)"
		#if [ $check_fitness_threshold -eq "1" ]
		#then
		#	echo "Reached threshold of 1e-200"
		#	exit
		#fi

#####################################################################################################
	for root_file in $f_seq_evolve/*
	do
		root_new_id="`echo "$root_file" | awk -F "root_new." '{print $2}'`"
		run_history="$f_scratch/$f_res_revolver/iteration.$s/$root_new_id/run.history"
		mkdir -p $f_scratch/$f_res_revolver/iteration.$s/$root_new_id

		# Used for domain constraints----------------------------------------------------------
		# root_file_hmm="$f_scratch/$f_res_revolver/iteration.$s/$root_new_id/$root_new_id.hmm"
		# $hmm_binaries/hmmscan --notextw $hmmfile $root_file > $root_file_hmm

		#printf "\n##### Evolving root sequence: iteration.$s/roots/root_new.$root_new_id ##########################################\n"
		#--------------------------------------------------------------------------
       		# Number of independent revolver "lineages" which will constitute "population" (number of individuals in each generation)
		#--------------------------------------------------------------------------
		# !!!!! Now I will run only one "lineage", but will do this not a star tre with many branches
		#--------------------------------------------------------------------------
        	for i in $search_revolver
        	do
			# MAIN REvolver simulation step: file+run
                        input_param_file="$f_scratch/$f_res_revolver/input_file.xml"
                        output_path="$f_scratch/$f_res_revolver/iteration.$s/$root_new_id/revolver.run.$i"
			results_hmm_comparison="$f_scratch/$f_res_revolver/iteration.$s/$root_new_id/summary.hmm.scores"
                        mkdir -p $output_path
                        ###########################################################################
                        #       (1) evolution with revolver
                        ###########################################################################
# SCRIPT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
printf "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>
<configdata  xsi:schemaLocation=\"http://www.cibiv.at/Revolver ./input_schema.xsd\" xmlns=\"http://www.cibiv.at/Revolver\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" >
        <config>
                <hmmdb path=\"$hmmfile\"/>
                <hmmfetch location=\"$hmmfetch_location\"/>
        </config>
        <model>
                <substitution name=\"$model_name\"/>
                <indel>
                        <insertion rate=\"$insertion_rate\">
                                <length distribution=\"geometric\" p=\"$len_distr_in\"/>
                        </insertion>
                        <deletion rate=\"$deletion_rate\">
                                <length distribution=\"geometric\" p=\"$len_distr_del\"/>
                        </deletion>
                </indel>
        </model>
        <tree path=\"$tree\" alpha=\"1.0\"  />
        <root>
                <inputSequence>
                        <fasta file=\"$root_file\"/>
                </inputSequence>
        </root>
        <output>
                <dir path=\"$output_path\" trueAlignment=\"false\" include=\"leaf\"/>
        </output>
</configdata>" > $input_param_file
# END of SCRIPT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			# For Domain constraint add: 
			# <hmmer file=\"$root_file_hmm\"/>
			# in <inputSequence> block
                        java -Xmx2G -Xms2G -cp "$REvolver" revolver $input_param_file
                        # delete white spaces in fasta file of new evolved sequences (introduced by revolver)
                        sed -i 's/ //g' $output_path/out.fa
                        rm $input_param_file

			############################################################################################
                        # (2) analysis of evolved sequence with HMMER
                        ############################################################################################
			query_fa=$output_path/out.fa
			#db_id=0
			hmm_score_comparison=()
			for d_hmm in $d_hmm_db
                        do
				#db_id="`echo $db_id+1 | bc`"
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
			done

			#score_SDR="$f_scratch/$f_res_hmmer/iteration.$s/aln-oxi-1SDRactino.fa.hmm/threshold.none/$root_new_id/all.root.${id}.$i.1"
			#score_LPOR="$f_scratch/$f_res_hmmer/iteration.$s/hmm-LPOR-proteo-15.fasta/threshold.none/$root_new_id/all.root.${id}.$i.1"


			#for y in {1..50}
			#do
			#	score_SDR_1="`head -n $y $score_SDR | tail -n 1`"
			#	score_LPOR_1="`head -n $y $score_LPOR | tail -n 1`"
			#	echo "$root_new_id $s $y $score_SDR_1 $score_LPOR_1" >> $results_hmm_comparison
			#	echo "$root_new_id $s $y $score_SDR_1 $score_LPOR_1" >> $gen_hmm_scores
			#done

			y=1
			echo "$root_new_id $s $y ${hmm_score_comparison[*]}" >> $results_hmm_comparison
			echo "$root_new_id $s $y ${hmm_score_comparison[*]}" >> $gen_hmm_scores

		done	# lineages
		done	# roots_new

		###################################################################
		# Check threshold: median reached proteo-LPOR median of 1e-186
		###################################################################
		gen_median="$(~/Dropbox/R.codes/LPOR-selection-threshold-2.r $gen_hmm_scores)"
		echo "GENERATION $s | generation median reached value of $gen_median" >> $threshold_history
		stop_rule="1e-186"
                check_fitness_threshold="$(~/Dropbox/R.codes/LPOR-compare-2-numbers.r $gen_median $stop_rule)"
                if [ $check_fitness_threshold -eq "1" ]
                then
                       echo "Reached threshold of median $stop_rule"
                       exit
                fi

		##################################################################
		# Collect top k sequences for the next generation
		##################################################################
		top_k_seq="$f_scratch/$f_res_revolver/iteration.$s/sample-next-parents"
                ~/Dropbox/R.codes/LPOR-selection-probabilistic-mod.r $gen_hmm_scores $top_k_seq	
		
		# This top 5 is used to compute the relative fitness of 5 best sequences in the generation
		# the simulation stops if this threshold reached <1e-200
		#top_5_seq="$f_scratch/$f_res_revolver/iteration.$s/top-k-sequences"
                #~/Dropbox/R.codes/LPOR-selection-of-top-k-evolved-seq-no-threshold.r $gen_hmm_scores $top_k $top_5_seq

		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		# - here call r code which will select 250 parents 
		#   for the next generation based on their fitness function:
		#   because fitness is based on E-scores, I choose probabilities as 1-E-score, 
		#   such that sequences with the lowest E-score have higher probability to be selected for the next generation
		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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

				# "Next threshold", i.e. relative fitness should be still computed from the best five sequences in the generation
				# - just copy stuff from the top_5_seq file
                        	#next_threshold="$f_scratch/$f_res_revolver/iteration.$s_1/compute-threshold"
                        	#echo "$f_scratch/$f_res_revolver/iteration.$s_1/roots/root_new.$next_id/out.fa $scores_aux" >> $next_threshold
                        
				#echo "- copying sequence - $next_r_file"
                        	#echo "- new root         - iteration.$s_1/roots/root_new.$next_id"
			done

			# next threshold
			#next_threshold="$f_scratch/$f_res_revolver/iteration.$s_1/compute-threshold"
			#cat $top_5_seq | while read -r line
                        #do
			#	root_new_id="`echo "$line" | awk -F " " '{print $1}'`"
			#	scores_aux="`echo "$line" | awk -F " " '{print $4 " " $5}'`"
			#	echo "$f_scratch/$f_res_revolver/iteration.$s/$root_new_id/revolver.run.1/out.fa $scores_aux" >> $next_threshold
			#done
				
                else
                       	echo "- No sequences were selected for the next generation roots"
			exit
                fi

		done	# iterations

fi

#cp -r $f_scratch/$f_res_hmmer $f_main_hmmer
#cp -r $f_scratch/$f_res_revolver $f_main_revolver
#rm -r $f_scratch/$f_res_hmmer
#rm -r $f_scratch/$f_res_revolver
#end=`date +%s`

#echo "Time: $((end-start))"
