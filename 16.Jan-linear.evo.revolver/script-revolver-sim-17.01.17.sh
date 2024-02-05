#!/bin/bash
#####################################################################################################
# PARAMETERS
#####################################################################################################
roots_ids="1 2 3 4 5"           # all roots         "$(seq 1 82)"
#search_type_trees="0 1 2"       # all types         "0 1 2"
#search_trees="$(seq 1 10)"      # all trees         "$(seq 1 50)"
search_revolver="$(seq 1 50)"  # all revolver runs "$(seq 1 100)"
iterations="$(seq 1 10)"
#----------------------------------------------------------------

# ---------------------------------------------------------------
sim_search_test="all"
if [ $sim_search_test == "1" ]
then
        roots_ids="1"
        search_revolver="$(seq 1 2)"
	iterations="$(seq 1 150)"
fi
if [ $sim_search_test == "2" ]
then
        roots_ids="6 7"
        search_revolver="$(seq 1 100)"
	iterations="$(seq 1 150)"
fi
if [ $sim_search_test == "all" ]
then
        roots_ids="$(seq 1 20)"
        search_revolver="$(seq 1 30)"
        iterations="$(seq 1 150)"
fi
#####################################################################################################
# Settings
#####################################################################################################
# for SDR aln
#--------------------------------------
	run_build_hmm="no"
	run_root_seq="no"
#-------------------------------------
# for trees
#-------------------------------------
	run_trees="no"
	run_trees_random="no"
#-------------------------------------
# for simulations
#-------------------------------------
	run_sim="yes"
	run_sim_rev="no"
#-------------------------------------
	run_sim_search="yes"
	run_sim_search_blast="no"
	run_sim_search_hmmer="yes"
	run_sim_search_hmmer_orgn="yes"
#-------------------------------------
	run_blast_refine="no"
	threshold_blast="1e-40"
#-------------------------------------
	run_hmmer_refine="no"
	threshold_hmm="1e-40"
	run_hmmer_refine_summary="no"
#-------------------------------------
# collecting seqs
#-------------------------------------
	run_collect="no"

#####################################################################################################
mkdir results
#mkdir trees
#===============================================================
REvolver="/project/olga-phylo/LPOR/2016.12.December/REvolver.jar"
#===============================================================
rootSEQs="input/sdr-oxi-alpha-proteo-ALN-3.fa"
input_aln="$rootSEQs"
#---------------------------------------------------------------
hmmfile="input/sdr-oxi-alpha-proteo-ALN-3.hmm"
hmm_binaries="/project/olga-phylo/LPOR/2016.12.December/hmmer-3.1b2-linux-intel-x86_64/binaries"
#---------------------------------------------------------------
hmmfetch_location="$hmm_binaries/hmmfetch"
hmmscan_location="$hmm_binaries/hmmscan"
#####################################################################################################
# MODEL details for REvolver
#####################################################################################################
rev_test="2"
if [ $rev_test == "1" ]
then
	model_name="LG"
	insertion_rate="0.001"
	deletion_rate="0.001"
	# Parameter for indel length distribution (geometric)
        len_distr_in="0.9"
        len_distr_del="0.9"	
fi
if [ $rev_test == "2" ]
then
        model_name="LG"
        insertion_rate="0.01"
        deletion_rate="0.01"
        # Parameter for indel length distribution (geometric)
        len_distr_in="0.9"
        len_distr_del="0.9"
fi
if [ $rev_test == "3" ]
then
        model_name="LG"
        insertion_rate="0.01"
        deletion_rate="0.01"
        # Parameter for indel length distribution (geometric)
        len_distr_in="0.5"
        len_distr_del="0.5"
fi



#len_distr_in="0`echo "scale=2; $(( ( RANDOM % 7 )  + 1 ))/10" | bc `"
#len_distr_del="0`echo "scale=2; $(( ( RANDOM % 7 )  + 1 ))/10" | bc `"
#for i in {1..10}
#do
#	echo "0`echo "scale=2; $(( ( RANDOM % 7 )  + 1 ))/10" | bc `"
#done
m="$model_name (in|del rates: $insertion_rate | $delition_rate)"
#####################################################################################################
# DETAILS for search with BLAST and HMMER against user datasets
#####################################################################################################
# HMMER --------------------------------------------------------
d_hmm_folder="/project/olga-phylo/LPOR/2017.01.January/16.Jan-linear.evo.revolver/databases/database.hmm"
d_hmm_db="LPOR-proteobacteria-15.hmm LPOR-bacteria-127.hmm"
#d_hmm_db="Aln-LPOR-bacteria-only-128-trimmed.fa.hmm"	#"LPOR-proteo-16.hmm" # LPOR-bacteria-128.hmm # LPOR-all-452.hmm"
#$hmm_binaries/hmmbuild $d_hmm $d_hmm_aln

# BLAST --------------------------------------------------------
d_blast_folder="/project/olga-phylo/LPOR/2017.01.January/analysis.LPOR/database.blast"
d_blast_db="LPOR-proteobacteria-sequences-16.fasta"
#####################################################################################################
# Building HMM model for SDR alignment
#####################################################################################################
if [ $run_build_hmm == "yes" ]
then
	rm $hmmfile
	rm $hmmfile.*

	printf "\n\nBuilding a profile HMM from an input alignment: $input_aln"
	$hmm_binaries/hmmbuild $hmmfile $input_aln

	printf "\n\nFormat an HMM database into a binary format for hmmscan.."
	$hmm_binaries/hmmpress $hmmfile
fi
#####################################################################################################
# Splitting root sequences from SDR alignment into separate files for further REvolver analysis
#####################################################################################################
if [ $run_root_seq == "yes" ]
then
	printf "\n\n"
	id=0
	cat $rootSEQs | while read line
	do
		q="$line"
		c="`echo $q | grep ">" | wc -l | awk -F " " '{print $1}'`"
		if [ $c -gt "0" ]
		then
			id="`echo "${id}+1" | bc `"
			root_seq="input/root-seq-${id}.fa"
			printf "Reading next root sequence from $rootSEQs: root-seq-${id}\n"
			echo "$q" > $root_seq
		else
			root_seq="input/root-seq-${id}.fa"
			echo "$q" >> $root_seq
		fi
	done
	printf "\n\n"
	for file in input/root-seq-*.fa
	do
		sed -i 's/-//g' $file # remove all gaps to obtain root sequences
		rm $file.hmm

		printf "Running hmmscan to profile $file\n"
                $hmm_binaries/hmmscan --notextw $hmmfile $file > $file.hmm
	done
fi
#####################################################################################################
# MAIN simulations: 
#	(1) evolution with revolver
#	(2) analysis of evolved sequences "Search":
#		- with HMMER
#		- with BLAST
#	(3) refinement of evolved sequences:
#               - with HMMER
#               - with BLAS
#####################################################################################################
if [ $run_sim == "yes" ]
then
	if [ $run_sim_rev == "yes" ]
	then
		#rm -r results/revolver
		mkdir -p results/revolver.$rev_test
		model_info="results/revolver.$rev_test/model.info"
        	echo "model_name $model_name"		 > $model_info
		echo "insertion_rate $insertion_rate"	>> $model_info
        	echo "deletion_rate $deletion_rate"	>> $model_info
		echo "len_distr_in $len_distr_in"	>> $model_info
        	echo "len_distr_del $len_distr_del"	>> $model_info
	fi

	#if [ $run_sim_search_hmmer == "yes" ]
        #then
        #        rm -r results/hmmer
        #fi
	for id in $roots_ids # Loop over different root sequences
	do
		#-----------------------------------------------------------
		printf "\n===============================================================\n"
		printf "\nStarting analysis for root-sequence-${id} ..\n"
		printf "\n===============================================================\n"
		#-----------------------------------------------------------
		root_file="input/root-seq-${id}.fa"
		root_file_hmm="input/root-seq-${id}.fa.hmm"
		
		# Get the score and the length of the root sequence
		if [ $run_sim_search_hmmer_orgn == "yes" ]
		then
			for d_hmm in $d_hmm_db
                        do
				query_orgn="input/root-seq-${id}.fa"
			
				printf "\n-----------------------------------------------------------------\n"
                        	printf "\nRunning HMMER search for original root sequence\n - ${query_orgn}\n - $d_hmm...\n"
                                results_orgn="results/hmmer.$rev_test/root.${id}/original_sequence/$d_hmm/threshold.none/root.${id}.hmmer"
                                mkdir -p "results/hmmer.$rev_test/root.${id}/original_sequence/$d_hmm/threshold.none"
                                $hmm_binaries/hmmsearch $d_hmm_folder/${d_hmm} ${query_orgn} > $results_orgn

				# Get sequence length   # wc -L - prints the length of the longest line
                                len=`tail -n +2 $query_orgn | tr --delete '\n' | wc -L | awk -F " " '{print $1}'`
				#echo "the length is $len..."
				
                                # Get the score of ONE evolved sequence
                                score_hmmer="results/hmmer.$rev_test/root.${id}/original_sequence/$d_hmm/threshold.none/results.root.${id}"

                                n=`tail -n 5 $results_orgn | head -n 1 | awk -F " " '{print $5}'`
                                if [ $n -gt 0 ]
                                then
                                	echo "${id} $len `head -n 15 $results_orgn | tail -n 1`" | cat >> $score_hmmer
					#echo "the score is `head -n 15 $results_orgn | tail -n 1`"
                                else
                                	echo "${id} $len  1 1 1 1 1 1 1 1 1" | cat >> $score_hmmer
					#echo "no hit..."
                                fi
			done
		fi	


		# Number of lineages revolver will simulate, i.e. number of revovler runs	
		for i in $search_revolver
		do
			printf "\n##### Revolver run $i ##########################################\n"
			root_file="input/root-seq-${id}.fa"
	                root_file_hmm="input/root-seq-${id}.fa.hmm"

		# Iterations: one step is one evolutionary process on one-branch tree
                for s in $iterations
                do	printf "\n===== Iteration $s =============================================\n"
			input_param_file="input/input_file_${id}.$i.$s.xml"
			output_path="results/revolver.$rev_test/root.${id}/revolver.run.$i/iteration.$s"
			mkdir -p $output_path

			tree="input/tree-1-branch.tree"
			###########################################################################
			#	(1) evolution with revolver
			###########################################################################
			if [ $run_sim_rev == "yes" ]
			then
		                printf "\n\nCreating the xml file for REvolver input parameters..\n - Iteration: rev.${i}.step.${s} \n"

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
                        <hmmer file=\"$root_file_hmm\"/>
                </inputSequence>
        </root>
        <output>
               	<dir path=\"$output_path\" trueAlignment=\"false\" include=\"leaf\"/>
        </output>
</configdata>" > $input_param_file
# END of SCRIPT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

				printf "\n\nRunning REvolver... \n - root sequence: \"$root_file\" \n - output: res.root.${id}/revolver.run.$i/iteration.$s/out.fa \n"
				#printf "\n\nRunning REvolver... \n - root sequence: \"$root_file\" \n - tree: \"$tree\" \n - model: \"$m\" \n - run $i (100 runs per tree) \n - output: res.root.${id}/revolver.run.$i/iteration.$s/out.fa \n"
				java -cp "$REvolver" revolver $input_param_file
				# delete white spaces in fasta file of new evolved sequences (introduced by revolver)
				sed -i 's/ //g' $output_path/out.fa

				root_file="$output_path/out.fa"
                		root_file_hmm="$output_path/out.fa.hmm"
				# Get hmm profile of evolved sequence
				#printf "\nRunning hmmscan for $root_file..\n"
                		$hmm_binaries/hmmscan --notextw $hmmfile $root_file > $root_file_hmm
				#mv $input_param_file $output_path/.
				rm $input_param_file

			fi
			############################################################################################
			# (2) analysis of evolved sequences "Search"
			############################################################################################
			if [ $run_sim_search == "yes" ]
			then

				query_fa=$output_path/out.fa

				#-----------------------------------#
				#	Search with HMMER	    #################################################################
				#-----------------------------------#
				if [ $run_sim_search_hmmer == "yes" ]
				then
				for d_hmm in $d_hmm_db
				do
					printf "\n-----------------------------------------------------------------\n"
					printf "\nRunning HMMER search for\n - ${query_fa}\n - $d_hmm...\n"
					results_h="results/hmmer.$rev_test/root.${id}/revolver.run.$i/$d_hmm/threshold.none/root.${id}.$i.$s.hmmer"
					mkdir -p "results/hmmer.$rev_test/root.${id}/revolver.run.$i/$d_hmm/threshold.none"
					$hmm_binaries/hmmsearch $d_hmm_folder/${d_hmm} ${query_fa} > $results_h

					# Get sequence length	# wc -L - prints the length of the longest line
                                        len=`tail -n +2 ${query_fa} | tr --delete '\n' | wc -L | awk -F " " '{print $1}'`
					#echo "length of sequence is $len"
					
					# Get the score of ONE evolved sequence
					score_hmmer="results/hmmer.$rev_test/root.${id}/revolver.run.$i/$d_hmm/threshold.none/all.root.${id}.$i"
					
					n=`tail -n 5 ${results_h} | head -n 1 | awk -F " " '{print $5}'`
                                        if [ $n -gt 0 ]
					then
						echo "${id} $i $s `head -n 15 $results_h | tail -n 1` $len" | cat >> $score_hmmer
						#echo "the hmmer score is `head -n 15 $results_h | tail -n 1`"
					else
						echo "${id} $i $s 1 1 1 1 1 1 1 1 1 $len" | cat >> $score_hmmer
						#echo "no hits.."
					fi
				done
				fi

				#-----------------------------------#
                                #       Refinement with HMMER       ################################################################
                                #-----------------------------------#
				if [ $run_hmmer_refine == "yes" ]
				then
				for d_hmm in $d_hmm_db
				do
				#==================================================================================
					file_hmm_result="results/hmmer.$rev_test/root.${id}/$d_hmm/threshold.none/root.${id}.t.${tr_id}.$i.hmmer"e
					# <file_hmm_result> is the same as <result_h>
					if [ -e $file_hmm_result ]
        				then
						file_hmm="results/hmmer.$rev_test/root.${id}/$d_hmm/threshold.${threshold_hmm}/root.${id}.t.${tr_id}.$i"
						file_all_info="results/hmmer.$rev_test/root.${id}/$d_hmm/threshold.${threshold_hmm}.all.info"
						mkdir -p results/hmmer.$rev_test/root.${id}/$d_hmm/threshold.${threshold_hmm}/
                				n=`tail -n 5 $file_hmm_result | head -n 1 | awk -F " " '{print $5}'`
						#printf "\n-----------------------------------\n"
                				#echo "(*) searching for hits in $query_fa"
                				if [ $n -gt 0 ]
                				then
                        				#echo "- $n hits found"
                        				head -n $((n+14)) $file_hmm_result | tail -n $((n)) > ${file_hmm}.hits.all
                  					sed -i '/inclusion/d' ${file_hmm}.hits.all
                        				#awk -F " " '{print $9}' > $file_hits.accNUM.all
                        				if [ -e ${file_hmm}.info ]
                        				then
                                				rm ${file_hmm}.info
                        				fi
                        				#echo "- collecting hits for threshold ${threshold_hmm}"
                        				cat ${file_hmm}.hits.all | while read line
                        				do
                                				s=`echo "$line" | awk -F " " '{print $1}' | awk -F "-" '{print $2}'`
                                				echo "$line" | awk -v t="$threshold_hmm" '$1 < t {print $9}' | cat >> ${file_hmm}.info
                        				done
							if [ -e ${file_hmm}.info ]
                        				then
                                				anCounts=`wc -l ${file_hmm}.info | awk -F " " '{print $1}'`
                                				#echo "$anCounts"
                                				if [ $anCounts -gt 0 ]
                                				then
                                        				file_seqs="${file_hmm}.sequences"
									echo "(*) searching for hits in $query_fa"
                                        				echo "- $anCounts hits matched chosen criteria"
                                        				if [ -e $file_seqs ]
                                        				then
                                                				rm $file_seqs
                                        				fi
                                        				cat ${file_hmm}.info | while read line
                                        				do
                                                				acc="$line"
                                                				# -m 1 -> grep only first match
                                                				grep -m 1 -A 1 "$acc" $query_fa | cat  >> $file_seqs
                                        				done
                                        				#echo "- sequences were collected in $file_seqs"
									echo "$anCounts sequences were collected in $file_seqs" >> $file_all_info
                                				else
                                        				#echo "- No hits found for the chosen threshold"
                                        				rm ${file_hmm}.info
                                				fi
                        				#else
                                				#echo "- No hits found for the chosen threshold"
                        				fi
                				#else
                        				#echo "- No hits found"
                				fi
					else
                				echo "ERROR: File with hmm search does not exist!"
                				exit 1
        				fi
				#==================================================================================
				done
				fi

				#-----------------------------------#
                                #       Search with BLAST           ###############################################
                                #-----------------------------------#
				if [ $run_sim_search_blast == "yes" ]
				then
				for d in $d_blast_db
				do
					printf "\n-----------------------------------------------------------------\n"
                        		printf "\nBlasting ${query_fa} against user dataset $d...\n"

                        		results_b="results/blast/root.${id}/$d/threshold.none/root.${id}.t.${tr_id}.$i.blast"
					mkdir -p "results/blast/root.${id}/$d/threshold.none"
                        		./blastp+ -db $d_blast_folder/$d -query ${query_fa} -outfmt 6 > $results_b
				done
				fi

				#-----------------------------------#
                                #       Refinement with BLAST       ###############################################
                                #-----------------------------------#
				if [ $run_blast_refine == "yes" ]
				then
					results_b="results/blast/root.${id}/threshold.none/root.${id}.t.${tr_id}.$i.blast"
                                        k=`wc -l $results_b | awk -F " " '{print $1}'`

                                        if [ $k -ne 0 ]
                                        then
                                                ##echo ">" >> $query
                                                s=`awk -F $'\t' '{print $1}' $results_b | sort | uniq | wc -l | awk -F " " '{print $1}'`
                                                echo "---> Found hits for $s sequences"
                                                echo "---> Results writen to root.${id}.t.${tr_id}.blast"
                                                seq=`awk -F $'\t' '{print $1}' $results_b | sort | uniq | xargs -n $s`
                                                echo "Sequences matching the threshold: $seq"
                                                ##for j in $seq
                                                ##do
                                                ##      seqAccNum="$j"

                                                        #echo "Found hits for sequence: $seqAccNum"
                                                        #echo "`grep -A 1 "${j}" ${i}`" | cat >> sequences.$threshold/${name}.fasta
                                                        #echo "`grep -A 1 "${j}" ${i}`" | cat >> sequences.$threshold/all.${folder}.fasta
                                                        #len=`grep -A 1 "${j}" ${i} | tail -n 1 | wc -L | awk -F " " '{print $1}'` # wc -L - the length of the longest line in the file

                                                        # -> copy line between /A/ and /B/, great for copying sequences in fasta and many lines
                                                        ##sed -n -e "/${seqAccNum}/,/>/ p" $i | sed -e '$d' >> sequences.$threshold/$folder/${name}.fasta
                                                        #echo "`sed -n -e "/${seqAccNum}/,/>/ p" $i | sed -e '$d'`"
                                                        ##sed -n -e "/${seqAccNum}/,/>/ p" $i | sed -e '$d' >> sequences.$threshold/$folder/all.fasta
                                                        # wc -L - prints the length of the longest line
                                                        ##len=`sed -n -e "/${seqAccNum}/,/>/ p" $i | sed -e '$d' | tail -n +2 | tr --delete '\n' | wc -L | awk -F " " '{print $1}'`
                                                ##done
                                        else
                                                echo "---> No hits for threshold e-value of ${threshold_blast} were found"
                                                ####rm $results_b
                                        fi
				fi # END of refinement BLAST
			#------------------------------------------------------------------------------------------------------
			fi # END of item (2) $run_sim_search
			#------------------------------------------------------------------------------------------------------
		done # END of iterations
		done # END of revolver runs
	done # END of root sequences ids
fi # END of simulations


############################################################################################
# (3) refinement of evolved sequences
############################################################################################
if [ $run_collect == "yes" ]
then
	printf "\n=========================================================================\n"
	printf '\nCollecting evolved sequences..\n'
	printf "\n=========================================================================\n"
        for id in $roots_ids
        do
		printf "for root.${id}..\n"
		for d_hmm in $d_hmm_db
                do
			query_folder="results/hmmer.$rev_test/root.${id}/$d_hmm/threshold.${threshold_hmm}"
                        file_seq_all="${query_folder}.seq_all"
			if [ -e $file_seq_all ]
			then
				echo "deleting existing file.."
				rm $file_seq_all
			fi
		done

		for z in $search_type_trees
                do
                for j in $search_trees
                do
                for i in $search_revolver
                do
			#-------------------------------------------------
			# HMMER - collect sequences for a given threshold
			if [ $run_hmmer_refine == "yes" ]
                        then
                                for d_hmm in $d_hmm_db
                                do
					file_hmm="results/hmmer.$rev_test/root.${id}/$d_hmm/threshold.${threshold_hmm}/root.${id}.t.$z.$j.$i"
					file_seqs="${file_hmm}.sequences"
					if [ -e $file_seqs ]
					then
						query_folder="results/hmmer.$rev_test/root.${id}/$d_hmm/threshold.${threshold_hmm}"
                                        	file_seq_all="${query_folder}.seq_all"
                                        	seq_id="root.${id}.t.$z.$j.$i"
						cp $file_seqs ${file_seqs}.auxiliary
						sed -i "s/>/>root.${id}.t.$z.$j.$i-/g" ${file_seqs}.auxiliary
						cat ${file_seqs}.auxiliary >> $file_seq_all
						rm ${file_seqs}.auxiliary
					fi
				done
			fi
			#------------------------------------------------
			# BLAST - collect sequences for a given threshold


		done
		done
		done

		if [ $run_hmmer_refine_summary == "yes" ]
		then
		for d_hmm in $d_hmm_db
                do
			printf "\n Making summary per tree type... root.${id}"
			file_summary="results/hmmer/$d_hmm/threshold.${threshold_hmm}.summary_per_tree_type"
			mkdir -p "results/hmmer/$d_hmm"
			file_all_info="results/hmmer/root.$id/$d_hmm/threshold.${threshold_hmm}.all.info"
			evolved_seqs=()
			for z in $search_type_trees
			do
				echo "tree type $z.."
                        	num_seq="`echo "$file_all_info" | grep "root.${id}.t.$z." | wc -l | awk -F " " '{print $1}'`"
				evolved_seqs=(${evolved_seqs[*]} $num_seq)
                        done
			echo "summary_for_root ${id} ${evolved_seqs[*]}" >> $file_summary
                done
		fi

	done

	if [ $run_hmmer_refine == "yes" ]
	then
	for d_hmm in $d_hmm_db
        do
		printf "\nSummary for HMMER search and threshold $threshold_hmm:\n"
		for id in $roots_ids 
		do
			query_folder="results/hmmer/root.${id}/$d_hmm/threshold.${threshold_hmm}"
                        file_seq_all="${query_folder}.seq_all"
			printf " - for root.${id}: `grep ">" $file_seq_all | wc -l | awk -F " " '{print $1}'` evolved sequences\n"
		done
	done
	fi
fi
#####################################################################################################
#	END of analysis
#####################################################################################################
printf "\n=========================================================================\n"
printf '\nAnalysis is completed.\n'
printf "\n=========================================================================\n"
























