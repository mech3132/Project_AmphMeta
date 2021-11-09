#!bin/bash

### Picrust ####
. activate picrust2

mkdir -p 02b_picrust
mkdir -p 02b_picrust/input

picrusttree='02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/picrust_placed/tree.nwk'
otutable='02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt/feature-table.biom'

## Pathways for fierer
#cd output_data/deblur/
#mkdir -p 02b_picrust
#mkdir -p 02b_picrust/input

#picrusttree='vsearch_99/picrust_placed/tree.nwk'
#otutable='vsearch_99/clustered_table-ASV50-taxfilt/feature-table.biom'

# Get copy number 
hsp.py -i 16S -t $picrusttree -o 02b_picrust/marker_predicted_and_nsti.tsv.gz -p 1 -n
# Get EC
hsp.py -i EC -t $picrusttree -o 02b_picrust/EC_predicted.tsv.gz -p 1
# KEGG orthologs?
hsp.py -i KO -t $picrusttree -o 02b_picrust/KO_predicted.tsv.gz -p 1



## Make metagenome pipeline
metagenome_pipeline.py -i $otutable \
-m 02b_picrust/marker_predicted_and_nsti.tsv.gz \
-f 02b_picrust/EC_predicted.tsv.gz \
-o 02b_picrust/EC_metagenome_out #--strat_out ## NOTE: strat_out takes up too much memory on my laptop, ran on fierer server
                       
metagenome_pipeline.py -i $otutable \
-m 02b_picrust/marker_predicted_and_nsti.tsv.gz \
-f 02b_picrust/KO_predicted.tsv.gz \
-o 02b_picrust/KO_metagenome_out #--strat_out


# Pathway-level inference
pathway_pipeline.py -i 02b_picrust/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
-o 02b_picrust/pathways_out_EC -p 1
pathway_pipeline.py -i 02b_picrust/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
-o 02b_picrust/pathways_out_KO -p 1


# Add functional descriptions
add_descriptions.py -i 02b_picrust/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o 02b_picrust/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i 02b_picrust/pathways_out_EC/path_abun_unstrat.tsv.gz -m METACYC \
                    -o 02b_picrust/pathways_out_EC/path_abun_unstrat_descrip.tsv.gz
                    
                    
add_descriptions.py -i 02b_picrust/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o 02b_picrust/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz


                    
## unzip            
gunzip 02b_picrust/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
gunzip 02b_picrust/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
gunzip -c 02b_picrust/marker_predicted_and_nsti.tsv.gz > 02b_picrust/marker_predicted_and_nsti.tsv

   
                    
