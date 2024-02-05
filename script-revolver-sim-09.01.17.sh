#!/bin/bash

# PARAMETERS ====================================================
roots_ids="1 2 3 4 5"           # all roots         "$(seq 1 82)"
search_type_trees="0 1 2"       # all types         "0 1 2"
search_trees="$(seq 1 10)"      # all trees         "$(seq 1 50)"
search_revolver="$(seq 1 50)"  # all revolver runs "$(seq 1 100)"

sim_search_test="0"

if [ $sim_search_test == "1" ]
then
        roots_ids="1"
        search_type_trees="0"
        search_trees="1"
        search_revolver="$(seq 1 10)"
fi
if [ $sim_search_test == "2" ]
then
        roots_ids="1 2 3 4 5"
        search_type_trees="0 1 2"
        search_trees="$(seq 1 10)"
        search_revolver="$(seq 1 50)"
fi

# Settings =====================================================
run_build_hmm="no"
run_root_seq="no"
#--------------------
run_trees="no"
run_trees_random="yes"
#--------------------
run_sim="yes"
run_sim_rev="no"
run_sim_search="yes"
run_sim_search_blast="no"
run_sim_search_hmmer="no"

run_blast_refine="no"
threshold_blast="1e-40"
run_hmmer_refine="yes"
threshold_hmm="1e-40"
run_hmmer_refine_summary="no"
#--------------------
run_collect="yes"

#===============================================================
mkdir results
mkdir trees
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
#---------------------------------------------------------------
# MODEL details for REvolver
model_name="WAG"
insertion_rate="0.02"
delition_rate="0.01"

#len_distr_in="0`echo "scale=2; $(( ( RANDOM % 7 )  + 1 ))/10" | bc `"
#len_distr_del="0`echo "scale=2; $(( ( RANDOM % 7 )  + 1 ))/10" | bc `"
#for i in {1..10}
#do
#	echo "0`echo "scale=2; $(( ( RANDOM % 7 )  + 1 ))/10" | bc `"
#done

m="$model_name (in|del rates: $insertion_rate | $delition_rate)"
#---------------------------------------------------------------
# DETAILS for search with BLAST and HMMER against user datasets
# HMMER --------------------------------------------------------
d_hmm_folder="/project/olga-phylo/LPOR/2017.01.January/analysis.LPOR/database.hmm"
d_hmm_db="LPOR-proteo-16.hmm" # LPOR-all-452.hmm"

#d_hmm_aln="LPOR-seq-2016-10-09-sp-457.fasta"
#d_hmm="LPOR-proteo-16.hmm"
#$hmm_binaries/hmmbuild $d_hmm $d_hmm_aln

# BLAST --------------------------------------------------------
#d="LPOR-only+2newJuneSeqs.fasta"       # database for blast search !!!!!!! Maybe you should also include the new proteoseqs from ggkbase
#d="LPOR-seq-2016-10-09-sp-457.fasta"
#makeblastdb+ -in ${d} -dbtype prot
d_blast_folder="/project/olga-phylo/LPOR/2017.01.January/analysis.LPOR/database.blast"
#d_blast_db="LPOR-seq-2016-10-09-sp-452-ohne-SDR.fasta"
d_blast_db="LPOR-proteobacteria-sequences-16.fasta"
#===============================================================
if [ $run_build_hmm == "yes" ]
then
	rm $hmmfile
	rm $hmmfile.*

	printf "\n\nBuilding a profile HMM from an input alignment: $input_aln"
	$hmm_binaries/hmmbuild $hmmfile $input_aln

	printf "\n\nFormat an HMM database into a binary format for hmmscan.."
	$hmm_binaries/hmmpress $hmmfile
fi
#==============================================================
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
#==============================================================
if [ $run_trees == "yes" ]
then
	rm trees/*
	printf "==============================================================="
	printf "\n\nCreating a tree to run the evolution on...\n"
	printf "==============================================================="
        # A star tree with many different branch lengths: a star tree = "independent" lineages
	printf "(" > trees/tree-star-10
	i=0
	branches="0.001 0.005 0.01 0.05 0.1 0.5 1.0 1.5 2.0 5.0"
	for b in $branches
	do
		i="`echo "$i+1" | bc`"
		printf "T$i:$b," >> trees/tree-star-10	
	done
	sed -i '$ s/.$//' trees/tree-star-10
	printf ")root;" >> trees/tree-star-10
	
        # Random Yule-Harding trees --------------------------------
	min_len=(0.0001 0.01 0.01)
	mean_len=(0.01 0.3 0.5)
	max_len=(0.8 1.5 3)

	if [ $run_trees_random == "yes" ]
	then
		for z in 0 1 2            # tree type: short, mixed, long branches
        	do
			for j in {1..50} # 100 trees per type of branches 
			do
				iqtree -r 20 -rlen ${min_len[$z]} ${mean_len[$z]} ${max_len[$z]} trees/tree.$z.$j
			done
		done
	fi
fi
#=============================================================
if [ $run_sim == "yes" ]
then
	if [ $run_sim_rev == "yes" ]
	then
		rm -r results/revolver
	fi
	#for id in {1..82} # Loop over different root sequences
	for id in $roots_ids
	do
		#-----------------------------------------------------------
		printf "\n===============================================================\n"
		printf "\nStarting analysis for root-sequence-${id} ..\n"
		printf "\n===============================================================\n"
		#-----------------------------------------------------------
		root_file="input/root-seq-${id}.fa"
		root_file_hmm="input/root-seq-${id}.fa.hmm"
		
		for z in $search_type_trees #0 1 2
		do
		#for j in {1..50}
		#for j in {1..10}
		for j in $search_trees
		do
		#for i in {1..100} # for each tree run REvolver 20 times
		#for i in {11..20}
		for i in $search_revolver
		do
			tr_id="$z.$j"
			input_param_file="input/input_file_${id}_tree_${tr_id}.$i.xml"
			output_path="results/revolver/res.root.${id}.t.${tr_id}.$i"
			mkdir -p $output_path

			tree="trees/tree.$z.$j"

			if [ $run_sim_rev == "yes" ]
			then
		                printf "\n\nCreating the xml file for REvolver input parameters..\n"
				# Parameter for indel length distribution (geometric)
				#len_distr_in="0`echo "scale=2; $(( ( RANDOM % 7 )  + 1 ))/10" | bc `"
				#len_distr_del="0`echo "scale=2; $(( ( RANDOM % 7 )  + 1 ))/10" | bc `"
				len_distr_in="0.25"
				len_distr_del="0.25"

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
                        <deletion rate=\"$delition_rate\">
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
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

				printf "\n\nRunning REvolver... \n - root sequence: \"$root_file\" \n - tree: \"$tree\" \n - model: \"$m\" \n - run $i (100 runs per tree) \n - output: res.root.${id}.t.${tr_id}.$i \n"
				java -cp "$REvolver" revolver $input_param_file
				# delete white spaces in fasta file of new evolved sequences (introduced by revolver)
				sed -i 's/ //g' $output_path/out.fa

			fi
			# =================================================================================================
			if [ $run_sim_search == "yes" ]
			then

				query_fa=$output_path/out.fa

				# HMMER----------------------------------------------------------------------------
				if [ $run_sim_search_hmmer == "yes" ]
				then
				for d_hmm in $d_hmm_db
				do
					printf "\n-----------------------------------------------------------------\n"
					printf "\nRunning HMMER search for ${query_fa} against user data set $d_hmm...\n"
					results_h="results/hmmer/root.${id}/$d_hmm/threshold.none/root.${id}.t.${tr_id}.$i.hmmer"
					mkdir -p "results/hmmer/root.${id}/$d_hmm/threshold.none"
					$hmm_binaries/hmmsearch $d_hmm_folder/${d_hmm} ${query_fa} > $results_h
				done
				fi

				if [ $run_hmmer_refine == "yes" ]
				then
				for d_hmm in $d_hmm_db
				do
				#==================================================================================
					file_hmm_result="results/hmmer/root.${id}/$d_hmm/threshold.none/root.${id}.t.${tr_id}.$i.hmmer" 
					# <file_hmm_result> is the same as <result_h>
					if [ -e $file_hmm_result ]
        				then
						file_hmm="results/hmmer/root.${id}/$d_hmm/threshold.${threshold_hmm}/root.${id}.t.${tr_id}.$i"
						file_all_info="results/hmmer/root.${id}/$d_hmm/threshold.${threshold_hmm}.all.info"
						mkdir -p results/hmmer/root.${id}/$d_hmm/threshold.${threshold_hmm}/
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

				# BLAST----------------------------------------------------------------------------
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
				fi
			fi
		done
		done
		done
	done
fi
#===================================================================================
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
			query_folder="results/hmmer/root.${id}/$d_hmm/threshold.${threshold_hmm}"
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
					file_hmm="results/hmmer/root.${id}/$d_hmm/threshold.${threshold_hmm}/root.${id}.t.$z.$j.$i"
					file_seqs="${file_hmm}.sequences"
					if [ -e $file_seqs ]
					then
						query_folder="results/hmmer/root.${id}/$d_hmm/threshold.${threshold_hmm}"
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




printf "\n=========================================================================\n"
printf '\nAnalysis is completed.\n'
printf "\n=========================================================================\n"
























