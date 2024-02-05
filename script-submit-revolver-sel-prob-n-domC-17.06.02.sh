#!/bin/bash

f="/project/olga-phylo/LPOR/2017.06.June/revolver-all/simulations/sim-n-dom-con-w-sel-prob-mod"

sim_test="3"
#sim_test="1"


#root_ids="11 13 22 9 42 77 1 46 17 83 7 81 36 27 23 61 5 80 65 26 29 60 39 8 53 14 32 76 4 62 63 64 50 33 41 59 74 3 35 15 52 75 2 44 67 45 30 10 79 20 49 73 56 40 18 16 38 24 57 69 66 12 28 43 48 70 71 72 68 78 58 55 31 21 34 47 51 6 25 54 37 19 82"

# running on franklin - cleaned
#root_ids="20" #root_ids="68" #"73" #83 33 4

# running on zuse

# running on babbage - cleaned
#root_ids="63" # 1 "80"

# running on muller - cleaned
#root_ids="36" #"19"

# running on laplace - cleaned
#root_ids="10" #"51"

# running on haldane - cleaned
#root_ids="2" #"27" #"59"
# running on meselson - cleaned
#root_ids=#"25" #"53"
# running on stahl - cleaned
#root_ids="8" #"34"
# running on kimura - cleaned - rerun 36 on muller
# running on shannon - cleaned
#root_ids="56" #"67"
#running on dayhoff - cleaned
#root_ids="38" #"71"

# fitch
#root_ids="31 47 54 78 82"

# ritchie - only those from 71 dataset
#root_ids="14 15 16 21 24 28 29 30 37 39 40 41 44 48 49 50 6 62 65 66 72 74 75 79" # deleted: 58 76

# zuse - only those from 71 dataset
root_ids="18" # 32 35 43 5 52 61"
# 38 56 - partially on zuse, partially elsewhere
# removed from analysis: 69


#for id in {1..83}
#for id in 48 49 50 51 52 53 54 55 56 57 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 83

for id in $root_ids
do
for rev_mod in 6 #1 2 3
do
	#qsub -N "R$id-p-$rev_mod" -q compute -l compute -cwd -j yes -l h=!zuse.cibiv.univie.ac.at  $f/script-revolver-sel-prob-n-domC-17.06.1.sh $id $sim_test $rev_mod 
	#$f/script-revolver-sel-prob-n-domC-17.06.1.sh $id $sim_test $rev_mod &
	#$f/script-start-from-i-iteration.sh $id $sim_test $rev_mod &
	#qsub -N "R$id-p-$rev_mod" -q desktops $f/script-revolver-sel-prob-n-domC-17.06.1.sh $id $sim_test $rev_mod
	$f/script-hmmer-score-of-mutants-17.07.25.sh $id $sim_test $rev_mod &

done
done
