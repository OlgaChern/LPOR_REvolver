#!/bin/bash

hmm_binaries="/project/olga-phylo/LPOR/2016.12.December/hmmer-3.1b2-linux-intel-x86_64/binaries"
d_hmm_folder="/project/olga-phylo/LPOR/2017.01.January/analysis.LPOR/database.hmm"
d_hmm_db="LPOR-all-452.hmm  LPOR-proteo-16.hmm"
query_folder="/project/olga-phylo/LPOR/2017.01.January/analysis.LPOR/database.blast"
queries="LPOR-proteobacteria-sequences-16.fasta LPOR-seq-2016-10-09-sp-452-ohne-SDR.fasta"

for d_hmm in $d_hmm_db
do
	for query_fa in $queries
	do
		results_h="./results.$query_fa-----$d_hmm"
		$hmm_binaries/hmmsearch $d_hmm_folder/${d_hmm} $query_folder/${query_fa} > $results_h
	done
done
