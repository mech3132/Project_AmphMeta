#!bin/bash

. activate qiime2-2021.8

mkdir -p 04a_diversity

## First, filter out archaea 
mkdir -p 04a_diversity/filtered_tables
qiime taxa filter-table \
--i-table 02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt.qza \
--i-taxonomy 02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/taxonomy_vsearch_80.qza \
--p-exclude 'Archaea' \
--o-filtered-table 04a_diversity/filtered_tables/table_noarchaea.qza

# Now, cut off by 2500
qiime feature-table filter-samples \
--i-table 04a_diversity/filtered_tables/table_noarchaea.qza \
--p-min-frequency 2500 \
--p-filter-empty-features \
--o-filtered-table 04a_diversity/filtered_tables/table_cut2500.qza

qiime tools export \
--input-path 04a_diversity/filtered_tables/table_cut2500.qza \
--output-path 04a_diversity/filtered_tables/table_cut2500 


biom summarize-table -i 04a_diversity/filtered_tables/table_cut2500/feature-table.biom \
-o 04a_diversity/filtered_tables/table_cut2500/feature-table_summary.txt

biom convert -i 04a_diversity/filtered_tables/table_cut2500/feature-table.biom \
--to-tsv -o  04a_diversity/filtered_tables/table_cut2500/feature-table.txt

## Filter tree
mkdir -p 04a_diversity/tree/
qiime tools import \
--input-path 03a_molecularClock_tree/final_tree_99_ASV50.nwk \
--output-path 04a_diversity/tree/calibrated_tree.qza \
--type Phylogeny[Rooted]

qiime phylogeny filter-tree \
--i-tree 04a_diversity/tree/calibrated_tree.qza \
--i-table 04a_diversity/filtered_tables/table_cut2500.qza \
--o-filtered-tree 04a_diversity/tree/calibrated_tree_filt.qza

qiime phylogeny filter-tree \
--i-tree 02a_download_ENAEBI_data_amphibian/qiime_output_deblur/vsearch_99/tree_99_ASV50.qza \
--i-table 04a_diversity/filtered_tables/table_cut2500.qza \
--o-filtered-tree 04a_diversity/tree/tree_filt.qza


qiime tools export \
--input-path 04a_diversity/tree/tree_filt.qza \
--output-path 04a_diversity/tree/tree_filt


qiime tools export \
--input-path 04a_diversity/tree/calibrated_tree_filt.qza \
--output-path 04a_diversity/tree/calibrated_tree_filt


# Rarefy
qiime feature-table rarefy \
--i-table 04a_diversity/filtered_tables/table_cut2500.qza \
--p-sampling-depth 2500 \
--o-rarefied-table 04a_diversity/filtered_tables/table_rare2500.qza 

qiime tools export \
--input-path 04a_diversity/filtered_tables/table_rare2500.qza \
--output-path 04a_diversity/filtered_tables/table_rare2500


biom convert -i 04a_diversity/filtered_tables/table_rare2500/feature-table.biom \
--to-tsv -o  04a_diversity/filtered_tables/table_rare2500/feature-table.txt

qiime diversity alpha \
--i-table 04a_diversity/filtered_tables/table_rare2500.qza  \
--p-metric observed_features \
--o-alpha-diversity 04a_diversity/observed_otus_vector.qza


qiime diversity alpha-phylogenetic \
--i-table 04a_diversity/filtered_tables/table_rare2500.qza \
--i-phylogeny 04a_diversity/tree/tree_filt.qza \
--p-metric faith_pd \
--o-alpha-diversity 04a_diversity/faith_pd_vector.qza

qiime diversity beta \
--i-table 04a_diversity/filtered_tables/table_rare2500.qza  \
--p-metric braycurtis \
--o-distance-matrix 04a_diversity/braycurtis_matrix.qza

qiime diversity beta \
--i-table 04a_diversity/filtered_tables/table_rare2500.qza  \
--p-metric jaccard \
--o-distance-matrix 04a_diversity/jaccard_matrix_rare.qza

qiime diversity beta \
--i-table 04a_diversity/filtered_tables/table_cut2500.qza  \
--p-metric jaccard \
--o-distance-matrix 04a_diversity/jaccard_matrix_cut.qza

qiime diversity beta-phylogenetic \
--i-table 04a_diversity/filtered_tables/table_rare2500.qza  \
--i-phylogeny 04a_diversity/tree/tree_filt.qza \
--p-metric unweighted_unifrac \
--o-distance-matrix 04a_diversity/uwu_matrix.qza

qiime diversity beta-phylogenetic \
--i-table 04a_diversity/filtered_tables/table_rare2500.qza  \
--i-phylogeny 04a_diversity/tree/tree_filt.qza \
--p-metric weighted_unifrac \
--o-distance-matrix 04a_diversity/wu_matrix.qza

### EXPORT ####
qiime tools export \
--input-path 04a_diversity/observed_otus_vector.qza \
--output-path 04a_diversity/observed_otus_vector

qiime tools export \
--input-path 04a_diversity/faith_pd_vector.qza \
--output-path 04a_diversity/faith_pd_vector


qiime tools export \
--input-path 04a_diversity/braycurtis_matrix.qza \
--output-path 04a_diversity/braycurtis_matrix


qiime tools export \
--input-path 04a_diversity/uwu_matrix.qza \
--output-path 04a_diversity/uwu_matrix


qiime tools export \
--input-path 04a_diversity/wu_matrix.qza \
--output-path 04a_diversity/wu_matrix


qiime tools export \
--input-path 04a_diversity/jaccard_matrix_rare.qza \
--output-path 04a_diversity/jaccard_matrix_rare


qiime tools export \
--input-path 04a_diversity/jaccard_matrix_cut.qza \
--output-path 04a_diversity/jaccard_matrix_cut

### Move into neat directory
mkdir -p 04a_diversity/exported
mv 04a_diversity/observed_otus_vector/alpha-diversity.tsv 04a_diversity/exported/observed_otus.tsv
mv 04a_diversity/faith_pd_vector/alpha-diversity.tsv 04a_diversity/exported/faith_pd.tsv
mv 04a_diversity/braycurtis_matrix/distance-matrix.tsv 04a_diversity/exported/bray_curtis.tsv
mv 04a_diversity/uwu_matrix/distance-matrix.tsv 04a_diversity/exported/unweighted_unifrac.tsv
mv 04a_diversity/wu_matrix/distance-matrix.tsv 04a_diversity/exported/weighted_unifrac.tsv
mv 04a_diversity/jaccard_matrix_rare/distance-matrix.tsv 04a_diversity/exported/jaccard_matrix_rare.tsv
mv 04a_diversity/jaccard_matrix_cut/distance-matrix.tsv 04a_diversity/exported/jaccard_matrix_cut.tsv

rm -r 04a_diversity/observed_otus_vector
rm -r 04a_diversity/faith_pd_vector
rm -r 04a_diversity/braycurtis_matrix
rm -r 04a_diversity/uwu_matrix
rm -r 04a_diversity/wu_matrix
rm -r 04a_diversity/jaccard_matrix_rare
rm -r 04a_diversity/jaccard_matrix_cut

