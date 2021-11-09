#!bin/bash

mkdir -p '03a_molecularClock_tree'

Rscript 03a_molecularClock_tree.R

cd 03a_molecularClock_tree


# Compare to basic UM tree
#python3 ../code/03a_make_ultrametric.py

## Tree string
treeString=$(cat Routput/filtered_tree.nwk)
#treeUMString=$(cat Routput/tree_Silva.nwk)

echo $treeString > treeString.txt
#echo $treeUMString > treeUMString.txt
echo -e "\n" > space.txt
### Make file


cat Routput/seqString.txt treeString.txt space.txt Routput/pathd8_nodes/* > input_file_pathd8.txt
#cat Routput/seqString.txt treeUMString.txt space.txt Routput/pathd8_nodes/* > input_file_pathd8_UM.txt
#cat qiime_output_deblur/molecularClockSeqs/seqLengthLine.txt qiime_output_deblur/molecularClockSeqs/treeString.txt qiime_output_deblur/molecularClockSeqs/space.txt qiime_output_deblur/molecularClockSeqs/pathd8_sets/* > qiime_output_deblur/molecularClockSeqs/input_file_pathd8.txt


cp ../code/PATHd8 .


./PATHd8 input_file_pathd8.txt calibrated_tree_99_ASV50.nwk
#./PATHd8 input_file_pathd8_UM.txt calibrated_tree_99_ASV50_UM.nwk
grep "d8 tree" calibrated_tree_99_ASV50.nwk | sed 's/d8 tree    : //g' > final_tree_99_ASV50.nwk
#grep "d8 tree" calibrated_tree_99_ASV50_UM.nwk | sed 's/d8 tree    : //g' > final_tree_99_ASV50_UM.nwk


cd ..



