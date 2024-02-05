#!/bin/bash
f="/project/olga-phylo/LPOR/2017.06.June/revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod"
f_in="$f/results-collection"
f_new="$f/results-REvolver-sim-1/summary-mutant-sequences"

for folder in $f_in/*
do
	r_num="`ls $folder | grep "root" | wc -l | awk -F " " '{print $1}'`"
	f_name="`echo "$folder" |  awk -F "results-collection/" '{print  $2}'`"
	mkdir -p $f_new/$f_name-$r_num

	echo "================================================================"
	echo "Copying from $f_name"
	echo "================================================================"
	c="`ls $folder/ | grep "aln" | wc -l | awk -F " " '{print $1}'`"
	if [ $c -gt 0 ]
	then
		cp $folder/aln* $f_new/$f_name-$r_num
	fi
	for i in $f_in/$f_name/root.*
	do
		f_root="`echo "$i" | awk -F "$f_name/" '{print $2}'`"
#		mkdir $f_new/$f_name-$r_num/$f_root
#		echo "copying from $f_root..."
		h="$f_new/$f_name-$r_num/$f_root"
#		echo "/*plot*"
		q="`ls $i | grep "E_score-in-each-generation" | awk -F "scores.." '{print $2}'`"
#		cp $i/results.TIGRFAM.hmm.scores..*-plot-E_score-in-each-generation.pdf $h/$q
		q="`ls $i | grep "plot-TIGRFAM" | awk -F "scores.." '{print $2}'`"
#		cp $i/results.TIGRFAM.hmm.scores..*-plot-TIGRFAM.pdf $h/$q


#		echo "/initial-root-seq*"
#		cp $i/initial-root-seq* $h/

#		c="`ls $i/ | grep ".lvl.1" | wc -l | awk -F " " '{print $1}'`"
#		if [ $c -gt 0 ]
#		then
#			mkdir -p $h/lvl.1/
#			echo "/*.lvl.1*"
#			cp $i/seq*.lvl.1* $h/lvl.1/
#		fi
#
#		c="`ls $i/ | grep ".lvl.2" | wc -l | awk -F " " '{print $1}'`"
#		if [ $c -gt 0 ]
#		then
#			mkdir -p $h/lvl.2/
#			echo "/*.lvl.2*"
#			cp $i/seq*.lvl.2* $h/lvl.2/
#		fi
#
#		c="`ls $i/ | grep ".lvl.3" | wc -l | awk -F " " '{print $1}'`"
#		if [ $c -gt 0 ]
#		then
#			mkdir -p $h/lvl.3/
#			echo "/*.lvl.3*"
#			cp $i/seq*.lvl.3* $h/lvl.3
#		fi
#
#
#		echo "/sequences.*.TIGRFAM.all"
#		cp $i/sequences.*.TIGRFAM.all $h/
	done
done
