#!/bin/bash

f="$1"          # file with names
out="$f-collected"

ff="$2"		# file with all sequences

rm $out

cat $f | while read line
do
	c="`echo "$line" | grep "bacteria" | wc -l | awk -F " " '{print $1}'`"
	#echo "$line" | grep "bacteria" | wc -l | awk -F " " '{print $1}'



	if [ "$c" == "1" ]
	then
		echo "$line"
        	echo "$line" >> $out
        	query=`echo $line | awk -F " " '{print $1}'`
        	sed -n -e "/$query/,/>/ p" $ff | sed '/>/d' | sed '/^$/d' >> $out
        else 
		c="`echo "$line" | grep "plant" | wc -l | awk -F " " '{print $1}'`"
		g="`echo "$line" | grep "algae" | wc -l | awk -F " " '{print $1}'`"
		

		#echo "$line" | grep "plant" | wc -l | awk -F " " '{print $1}'
		#echo "$line" | grep "algae" | wc -l | awk -F " " '{print $1}'

		if [ $c == "0" ]
		then
			if [ $g == "0" ]
			then
				echo "$line"
				echo "$line" >> $out
				query=`echo $line | awk -F " " '{print $1}'`
                		sed -n -e "/$query/,/>/ p" $ff | sed '/>/d' | sed '/^$/d' >> $out
			fi
		fi
	fi
done
