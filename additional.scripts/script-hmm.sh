#!/bin/bash

# Default values----------------
d="NA"
m="NA"
t="1" 		# get all hits
a="NA"
#-------------------------------
while getopts "hd:m:t:a" opt; do
  case $opt in
	h)
	echo "============================================================================"
	echo " To start use:" 
	echo "		./script-hmm.sh -d <DB_Name>"
	echo " and at least one of the following options: -m <HMM_Model> or -a"
	echo ""
	echo " Make sure that the your prefered database is present in the folder Databases/"
	echo "============================================================================"
	echo ""
	echo "	-d <DB_NAME> ----- the folder with corresponding name containing fasta files to be analysed (Databases/DB_NAME)"
	echo "	-m <HMM_MODEL> --- file with model to be used for hmmsearch"
	echo "	-t <VALUE> ------- threshold for hits inclusion (E-value), default all hits: t=1"
	echo "	-a --------------- get hits accession numbers and corresponding sequences from input files,  use if hmm search was already performed!"
	echo "============================================================================"
	exit 1
	;;
	d)
		d="$OPTARG"
	;;
    	m)
		m="$OPTARG"
	#	a="getAccNUM"
      	;;
        t)
                t="$OPTARG"
        ;;
	a)
		a="getAccNUM"
	;;
    	\?)
      		echo "Invalid option: -$OPTARG"
      		exit 1
      	;;
    	:)
      		echo "Option -$OPTARG requires an argument."
      		exit 1
      	;;
  esac
done

if [ $m = "NA" -a $a = "NA" ]
then
	echo "--------------------------------" 
	echo " So what are we planning to do?"
	echo " -h: help menu"
	echo "--------------------------------"
	exit 1
fi

echo ""
echo "-----------------------------------------------------------------------"
echo "Starting the analysis....."
if [ -d Databases/$d ]
then
	echo "DB name  : $d"
	echo "DB folder: Databases/$d"
else
	echo "ERROR: Such folder does not exit! -> Databases/$d"
	exit 1
fi

if [ $m != "NA" ]
then
	if [ -e $m ] 
	then
		echo "HMM model: $m"
	else
		echo "ERROR: Model file does not exist! -> $m "
		exit 1
	fi
fi

if [ $a != "NA" ]
then
	echo "Inclusion threshold: $t"
fi


ls Databases/$d > list.$d
mkdir -p results/$d

cat list.$d | while read line
do
echo "------------------------------------------------------------------------"
	q="$line"
	file_hmm_result="results/$d/results.$q"
	file_hits="results/$d/results.$q"
	file_accNUM="$file_hits.$t.accNUM"

	if [ $m != "NA" ]
	then
		echo "(*) searching for hmm hits in $q"
		programs/hmmsearch $m Databases/$d/$q > $file_hmm_result
	fi

	if [ $a = "getAccNUM" ]
	then
	if [ -e $file_hmm_result ]
	then
		n=`tail -n 5 $file_hmm_result | head -n 1 | awk -F " " '{print $5}'`
		echo "(*) searching for hits in $q"
		if [ $n -gt 0 ]
		then
			echo "- $n hits found"
			head -n $((n+15)) $file_hmm_result | tail -n $((n+1)) > $file_hits.hits.all
			sed -i '/inclusion/d' $file_hits.hits.all
			#awk -F " " '{print $9}' > $file_hits.accNUM.all
			if [ -e $file_accNUM ]
                        then
				rm $file_accNUM
			fi

			echo "- collecting hits for threshold $t"
			cat $file_hits.hits.all | while read line
			do
				s=`echo "$line" | awk -F " " '{print $1}' | awk -F "-" '{print $2}'`
				echo "$line" | awk -v t="$t" '$1 < t {print $9}' | cat >> $file_accNUM
			done

			if [ -e $file_accNUM ]
			then
				anCounts=`wc -l $file_accNUM | awk -F " " '{print $1}'`
				#echo "$anCounts"
				if [ $anCounts -gt 0 ]
				then
					file_seqs="$file_hits.$t.sequences"
					echo "- $anCounts hits matched chosen criteria"
					if [ -e $file_seqs ]
					then
						# The question did not work. Donno why..
						#echo -e "File $file_seqs already exists. Would you like to rewrite it? y/n"
						#read reply
						#if [ $reply -eq "y" ] || [ $reply -eq "yes" ]
						#then
						rm $file_seqs
						##else if [ $reply -eq "n" ] || [ $reply -eq "no" ]
						#else
						#	exit 1
						#	#fi
						#fi
					fi
                			cat $file_accNUM | while read line
                			do
                        			acc="$line"
						# -m 1 -> grep only first match
                        			grep -m 1 -A 1 "$acc" Databases/$d/$q | cat  >> $file_seqs
                			done
                			echo "- sequences were collected in $file_seqs"
				else
					echo "- No hits found for the chosen threshold"
					rm $file_accNUM
				fi
			else
				echo "- No hits found for the chosen threshold"
			fi	
 		else
			echo "- No hits found"
		fi
	else
		echo "ERROR: File with hmm search does not exist! Use -m option first. Check -h for help."
		exit 1
	fi
	fi
	
done
echo "------------------------------------------------------------------------"
echo "Output files written to: results/$d/ "
echo "	*.hits.all ----- hmmsearch results"
echo "	*.$t.accNUM ---- accession numbers of hits corresponding to inclusion threshold $t"
echo "	*.$t.sequences - sequences in fasta format for accession numbers from *.$t.accNUM"
echo "------------------------------------------------------------------------"
