#!/bin/bash

t="$1"
f="$2"	# input file to be processed: names will be changed from t.COL1 to t.COL2

cat $t | while read line
do
	old_name="`echo $line | awk -F " " '{print $1}'`"
	new_name="`echo $line | awk -F " " '{print $2}'`"
	echo "$old_name ------> $new_name"
	#sed -i 's/$old_name/$new_name/g' $f
	perl -pi -e "s/$old_name/$new_name/g" $f
done






# Just a useful command line
# awk -F $"\t" '{print $2 "\t" $1 "_" $3}' Names-table-LPOR-sp-tag-1.txt >>Names-table-LPOR-sp-tag-1-mod.txt
# ./script-rename.sh Names-table-LPOR-sp-tag-1-mod.txt /project/olga-phylo/LPOR/January.2016/improved-tree-inference/LPOR/LPOR.trees/all.50.LPOR.trees
