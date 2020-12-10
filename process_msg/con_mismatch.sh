#!/bin/bash

#### concatenate mismatch files and use the file name to keep track of origin

for i in *-matchMismatch.csv; do
	sed -e '/contig,start,end,par1_par2_diff,like_par1,like_par2/d' $i > tmp.txt && mv tmp.txt $i
	#sed '/^$/d' $i
	awk '{print FILENAME (NF?",":"") $0}' $i > tmp.txt && mv tmp.txt $i
	sed -e 's/-matchMismatch.csv//g' $i > tmp.txt && mv tmp.txt $i
	sed -e 's/indiv//g' $i > tmp.txt && mv tmp.txt $i
done

cat *-matchMismatch.csv > all-matchMismatch.csv
