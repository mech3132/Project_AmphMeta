#!bin/bash

. activate qiime2-2021.8

mkdir -p 04c_FAPROTAX
cd 04c_FAPROTAX

cp -r ../code/FAPROTAX_1.2.4/* .

sed 's/OTUID/\#OTUID/g' ../04b_data_filt_and_combine/downstream/otu_rare_withTaxa.txt > otu_rare_withTaxa.tsv

#./collapse_table.py \
#--input_table otu_rare_withTaxa.tsv \
#--row_names_are_in_column "#OTU ID" \
#--out_collapsed otu_rare_faprotax_collapsed.txt \
#--output_format_collapsed classical \
#-g FAPROTAX.txt

./collapse_table.py -i otu_rare_withTaxa.tsv \
-o func_table.tsv -g FAPROTAX.txt \
--column_names_are_in last_comment_line -c "#" -v \
-d "OTUID"
cd ..
