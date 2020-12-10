#!/bin/bash

#### specify the plate for the multiplex samples

plate="Dys2"

#### clean up concatenated mismatch file for headers and name of sample

## remove headers

grep -E -v "indiv\w+-\w+.csv:ancestry,contig,start,end,par1_par2_diff,like_par1,like_par2" all-matchMismatch_example_preprocess.csv > all-matchMismatch_example_processInt1.csv 

## remove any samples homozygous for parent 2 (samples should be heterozygous or homozygous for parent 1 due to crossing scheme)

grep -v "homozygous_par2" all-matchMismatch_example_processInt1.csv  > all-matchMismatch_example_processInt2.csv

## replace name of fileand homozygous par1 part
sed 's/-matchMismatch.csv:homozygous_par1//g' all-matchMismatch_example_processInt2.csv > all-matchMismatch_example_processInt3.csv

## replace indiv with plate name (specify plate before running file)

sed "s/indiv/$plate,/g" all-matchMismatch_example_processInt3.csv > all-matchMismatch_example_processInt4.csv

#### add header

echo "plate,well,barcode,contig,start,end,markers1,markers2,likelihood" | cat - all-matchMismatch_example_processInt5.csv > all-matchMismatch_example_processed.csv

#### remove intermediate files

rm *processInt*.csv

#### add information about parent of origin and brood size manually or merging data files

