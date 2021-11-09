#!bin/bash

#. activate qiime2-2020.11

. activate qiime2-2021.8
###### Progress bar ########
BAR='####################'

#ERRlistAcc='./01_parsed_meta/ERR_toDownload_sample_accessions.txt.txt'
SRRlistAcc='01_parsed_meta/SRR_toDownload_sample_accessions_COMPLETE.txt'
nonEBIlistAcc='01_parsed_meta/nonEBI_toDownload_sample_accessions_COMPLETE.txt'
kuenAdd='01_parsed_meta/Kuen_toAdd_COMPLETE.txt'

mkdir -p 02a_download_ENAEBI_data_amphibian
cd 02a_download_ENAEBI_data_amphibian
mkdir -p qiime_output_deblur
mkdir -p qiime_output_deblur/databases/
mkdir -p qiime_output_deblur/repsets
mkdir -p qiime_output_deblur/tables

mkdir -p allRawData

## Download classife rin advance

## Try assigning taxonomy temporarily
if [ -e qiime_output_deblur/databases/silva-138-99-nb-weighted-classifier.qza ]
then
echo "silva 138 classifier already downloaded"
else
	wget -O qiime_output_deblur/databases/silva-138-99-nb-weighted-classifier.qza https://data.qiime2.org/2021.8/common/silva-138-99-nb-weighted-classifier.qza
fi

if [ -e qiime_output_deblur/databases/sepp-refs-silva-128.qza ]
then
echo "silva sepp already downloaded"
else
	wget -O qiime_output_deblur/databases/sepp-refs-silva-128.qza https://data.qiime2.org/2021.8/common/sepp-refs-silva-128.qza
fi


### Download SILVA 138 for vsearch
#Michael S Robeson II, Devon R O’Rourke, Benjamin D Kaehler, Michal Ziemski, Matthew R Dillon, Jeffrey T Foster, Nicholas A Bokulich. RESCRIPt: Reproducible sequence taxonomy reference database management for the masses. bioRxiv 2020.10.05.326504; doi: https://doi.org/10.1101/2020.10.05.326504

if [ -e qiime_output_deblur/databases/silva-138-99-seqs-515-806.qza ]
then
echo "silva 138 repset"
else
	wget -O qiime_output_deblur/databases/silva-138-99-seqs-515-806.qza https://data.qiime2.org/2021.8/common/silva-138-99-seqs-515-806.qza
	wget -O qiime_output_deblur/databases/silva-138-99-tax-515-806.qza https://data.qiime2.org/2021.8/common/silva-138-99-tax-515-806.qza
	
	qiime tools export \
	--input-path qiime_output_deblur/databases/silva-138-99-seqs-515-806.qza \
	--output-path qiime_output_deblur/databases/silva-138-99-seqs-515-806
	
	qiime tools export \
	--input-path qiime_output_deblur/databases/silva-138-99-tax-515-806.qza \
	--output-path qiime_output_deblur/databases/silva-138-99-tax-515-806
	
#	bash 02a_get_molecClockSeqs.sh
fi


### Go through list of accessions and download from EBI/ENA


#while read sampacc
#do
#if [ -e $sampacc/dada2 ]
#then
#	rm -r $sampacc/dada2
#fi
#done < ../$SRRlistAcc

nSamps=$(cat ../$SRRlistAcc | wc -l) 
currSamp=1
while read sampacc
do
let "tempPerc=($currSamp*100/$nSamps)"
let "endBar=($tempPerc*20/100)"
echo -ne "\r $tempPerc % [ ${BAR:0:$endBar}"
# Get folders
fold1=${sampacc:0:6}
lastDig=${sampacc: 9}
if (( ${#lastDig}==1 ))
then
	nzeros=00
else
	nzeros=0
fi
#wget -q -P $sampacc/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$fold1/00$lastDig/$sampacc/*
# TEMP:SRR5835045 

if [ -e allRawData/$sampacc ]
then
	echo "$sampacc exists"
else
	#wget -q -P allRawData/$sampacc/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$fold1/$nzeros$lastDig/$sampacc/*
	echo "FTP NO WORKING RIGHT NOW; COMMENTING OUT"
# Store un-downloaded files here
: > MANUAL_DOWNLOAD.txt
if [ ! -e allRawData/$sampacc ]
then
echo "Saving to manual download"
echo $sampacc >> MANUAL_DOWNLOAD.txt
fi # manual downloade

fi # check if sampacc folder exists

### Run QIIME stuff

#### Check if imported; import if not
# Is single or F/R lane?
if [ -e allRawData/$sampacc/${sampacc}_1.fastq.gz ] 
then
	temppwd=$(echo "$(pwd)/allRawData/$sampacc/${sampacc}_1.fastq.gz")
else
	temppwd=$(echo "$(pwd)/allRawData/$sampacc/${sampacc}.fastq.gz")
fi

### For deblur, can't have underscores
sampacc2=$(echo $sampacc | sed 's/_/-/g')


if [ -e allRawData/$sampacc/$sampacc2-demux.qza ]
then
echo "Already ran deblur for $sampacc"

else
# Remove the underscores, if any
#gzip -dk $temppwd
#temppwdunzip=$(echo $temppwd |sed 's/.gz//g' )
#temppwd2unzip=$(echo $temppwd2 |sed 's/.gz//g' )
#cat $temppwdunzip | sed "s/$sampacc/$sampacc2/g" > $temppwd2unzip
#gzip $temppwd2unzip
#rm $temppwdunzip

echo -e "sample-id\tabsolute-filepath" > allRawData/$sampacc/${sampacc2}_manifest.txt
echo -e "${sampacc2}\t$temppwd" >> allRawData/$sampacc/${sampacc2}_manifest.txt


# Import R1 as single demultiplexed sample
qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path allRawData/$sampacc/${sampacc2}_manifest.txt \
--output-path allRawData/$sampacc/$sampacc2-demux.qza \
--input-format SingleEndFastqManifestPhred33V2

fi # deblur-demux exists?



if [ -e allRawData/$sampacc/deblur ] 
then
	echo "deblur already run for $sampacc"

else

# Re-name underscored sample names-- "deblur cannot operate on sample IDs that contain underscores"

# Run deblur -- default filter metrics, but included here for clarity
qiime deblur denoise-16S \
--i-demultiplexed-seqs allRawData/$sampacc/$sampacc2-demux.qza \
--p-trim-length 150 --p-min-reads 10 \
--p-min-size 2 \
--p-no-hashed-feature-ids \
--output-dir allRawData/$sampacc/deblur


mv allRawData/$sampacc/deblur/table.qza qiime_output_deblur/tables/$sampacc2-table.qza
mv allRawData/$sampacc/deblur/representative_sequences.qza qiime_output_deblur/repsets/$sampacc2-representative_sequences.qza

fi


if [ -e allRawData/$sampacc/${sampacc}_2.fastq.gz ]
then
rm -r allRawData/$sampacc/${sampacc}_2.fastq.gz ######### REHIGHLIGHT when doing full batch to save HD space
echo "Removed extra reverse fastq file $sampacc"
fi


let "currSamp++"
done < ../$SRRlistAcc










##### Run the dryad/prest/tim datasets ####

### MOVE PREST FILES THAT I WANT TO KEEP
while read nonEBI
do
echo $nonEBIt
#done < ../$nonEBIlistAcc 

tempSet=$(echo $nonEBI | sed 's/ .*$//g')
sampacc=$(echo $nonEBI | sed 's/^.* //g')


mkdir -p allRawData/$sampacc

if [ -e allRawData/$sampacc/$sampacc.fastq.gz ]
then
echo "already copied $sampacc"
else
cp ../00_nonEBI_data/$tempSet/split_fastq_by_sampleid/$sampacc.fastq.gz \
allRawData/$sampacc/$sampacc.fastq.gz
fi
temppwd=$(echo "$(pwd)/allRawData/$sampacc/${sampacc}.fastq.gz")

### RUN DEBLUR

### For deblur, can't have underscores
sampacc2=$(echo $sampacc | sed 's/_/-/g')

#if [ $sampacc2 == $sampacc ]
#then
#echo "No underscores to remove from sampleID for deblur"
#else

if [ -e allRawData/$sampacc/$sampacc2-demux.qza ]
then
echo "Already imported $sampacc"

else
# Remove the underscores, if any
#gzip -dk $temppwd
#temppwdunzip=$(echo $temppwd |sed 's/.gz//g' )
#temppwd2unzip=$(echo $temppwd2 |sed 's/.gz//g' )
#cat $temppwdunzip | sed "s/$sampacc/$sampacc2/g" > $temppwd2unzip
#gzip $temppwd2unzip
#rm $temppwdunzip

echo -e "sample-id\tabsolute-filepath" > allRawData/$sampacc/${sampacc2}_manifest.txt
echo -e "${sampacc2}\t$temppwd" >> allRawData/$sampacc/${sampacc2}_manifest.txt


# Import R1 as single demultiplexed sample
qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path allRawData/$sampacc/${sampacc2}_manifest.txt \
--output-path allRawData/$sampacc/$sampacc2-demux.qza \
--input-format SingleEndFastqManifestPhred33V2

fi # deblur-demux exists?
#fi # sampacc2==sampacc


if [ -e allRawData/$sampacc/deblur ] 
then
	echo "deblur already run for $sampacc"

else

# Re-name underscored sample names-- "deblur cannot operate on sample IDs that contain underscores"

# Run deblur -- default filter metrics, but included here for clarity
qiime deblur denoise-16S \
--i-demultiplexed-seqs allRawData/$sampacc/$sampacc2-demux.qza \
--p-trim-length 150 \
--p-min-reads 0 \
--p-min-size 2 \
--p-no-hashed-feature-ids \
--output-dir allRawData/$sampacc/deblur


mv allRawData/$sampacc/deblur/table.qza qiime_output_deblur/tables/$sampacc2-table.qza
mv allRawData/$sampacc/deblur/representative_sequences.qza qiime_output_deblur/repsets/$sampacc2-representative_sequences.qza

fi


done < ../$nonEBIlistAcc 





#### Now, finally do Jordan's sets

if [ ! -e qiime_output_deblur/tables/kuen-table-filt.qza ]
then
cd ../00_nonEBI_data/Kuen/

qiime tools import \
--input-path GlobalDeblurMiseqV4.final.biom \
--type FeatureTable[Frequency] \
--output-path kuen_table.qza

grep ">" GlobalDeblurMiseqV4.final.seqs.fa > features_only.txt

echo "Translating only sequences to upper case"
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' GlobalDeblurMiseqV4.final.seqs.fa > upper_case_seqs.fasta

#: > upper_case_seqs.fasta
#while read f 
#do
#grep $f -A1 GlobalDeblurMiseqV4.final.seqs.fa | grep -v $f | tr [:lower:] [:upper:] >> upper_case_seqs.fasta
#done <features_only.txt

#cat GlobalDeblurMiseqV4.final.seqs.fa | tr [:lower:] [:upper:] > upper_case_seqs.fasta
#This doesn't work so well because it changes the feature NAME.

qiime tools import \
--input-path upper_case_seqs.fasta \
--type FeatureData[Sequence] \
--output-path kuen_repset.qza


qiime feature-table filter-samples \
--i-table kuen_table.qza \
--m-metadata-file ../../01_parsed_meta/Kuen_toAdd_COMPLETE.txt \
--o-filtered-table kuen_table_filt.qza

qiime tools export \
--input-path kuen_table_filt.qza \
--output-path kuen_table_filt

biom summarize-table -i kuen_table_filt/feature-table.biom -o kuen_table_filt_summary.txt


qiime feature-table filter-seqs \
--i-data kuen_repset.qza \
--i-table kuen_table_filt.qza \
--o-filtered-data kuen-repset-filt.qza



cd ../../02a_download_ENAEBI_data_amphibian

cp ../00_nonEBI_data/Kuen/kuen-repset-filt.qza qiime_output_deblur/repsets/
cp ../00_nonEBI_data/Kuen/kuen_table_filt.qza qiime_output_deblur/tables/kuen-table-filt.qza

fi #Jordan import and move




########### Cross-checking files to keep
#ls qiime_output_deblur/repsets/ | sed 's/-representative.*$//g' > qiime_output_deblur/list_of_samples_downloaded.txt



#### Remove all samples that are NO LONGER relevant


#awk -F '\t' '{ print $2 }' ../$nonEBIlistAcc > allFoldersCheck.txt
#cat ../$SRRlistAcc >> allFoldersCheck.txt

#ls allRawData > allFoldersPresent.txt

#while read f
#do

#if grep -q $f allFoldersCheck.txt 
#then
#echo "keep $f"
#else
#echo "REMOVE $f"
#fi

#done <allFoldersPresent.txt

# COPY TO FIERER
#scp qiime_output_deblur/repsets/* mchen@microbe.colorado.edu:/data/mchen/Project_metaanalysis/input_data/deblur/repsets/

#scp qiime_output_deblur/tables/* mchen@microbe.colorado.edu:/data/mchen/Project_metaanalysis/input_data/deblur/tables/

#### MERGE ALL DATA


## Add path to molecular clock sequences
#qiime tools import \
#--type FeatureData[Sequence] \
#--input-path qiime_output_deblur/databases/silva-138-99-seqs-515-806/seqs_for_molecClock.fasta \
#--output-path qiime_output_deblur/databases/silva-138-99-seqs-515-806/seqs_for_molecClock.qza



### Too big to merge myself-- doing it on fierer comp

if [ ! -e qiime_output_deblur/merged-table.qza ]
then
echo "merging tables and repsets"
qiime feature-table merge \
--i-tables qiime_output_deblur/tables/* \
--o-merged-table qiime_output_deblur/merged-table.qza


qiime feature-table merge-seqs \
--i-data qiime_output_deblur/repsets/* \
--o-merged-data qiime_output_deblur/merged-repset.qza
else
echo "Already merged tables and repsets"
fi

#scp mchen@microbe.colorado.edu:/data/mchen/Project_metaanalysis/output_data/deblur/merged-table.qza qiime_output_deblur/merged-table.qza
#scp mchen@microbe.colorado.edu:/data/mchen/Project_metaanalysis/output_data/deblur/merged-repset.qza qiime_output_deblur/merged-repset.qza



## cluster at 99 de novo

if [ ! -e qiime_output_deblur/vsearch_99 ]
then
echo "Vsearch clustering 99 de novo"
qiime vsearch cluster-features-de-novo \
--i-sequences qiime_output_deblur/merged-repset.qza \
--i-table qiime_output_deblur/merged-table.qza \
--p-perc-identity 0.99 \
--output-dir qiime_output_deblur/vsearch_99

qiime tools export \
--input-path qiime_output_deblur/vsearch_99/clustered_table.qza \
--output-path qiime_output_deblur/vsearch_99/clustered_table

qiime tools export \
--input-path qiime_output_deblur/vsearch_99/clustered_sequences.qza \
--output-path qiime_output_deblur/vsearch_99/clustered_sequences
else
echo "Already clustered 99 de novo"
fi


## Remove 50
#if [ ! -e qiime_output_deblur/merged-table-ASV50.qza  ]
#then
#echo "filtering raw 50"
#qiime feature-table filter-features \
#--i-table qiime_output_deblur/merged-table.qza \
#--p-min-frequency 50 \
#--o-filtered-table qiime_output_deblur/merged-table-ASV50.qza 

#qiime featuer-table filter-seqs \
#--i-data qiime_output_deblur/merged-repset.qza \
#--i-table qiime_output_deblur/merged-table-ASV50.qza \
#--o-filtered-data qiime_output_deblur/merged-repset-ASV50.qza
#else 
#echo "already filtered ASV50 for raw seqs"
#fi

if [ ! -e qiime_output_deblur/vsearch_99/clustered_table-ASV50.qza  ]
then
echo "filtering clustered99 50"
qiime feature-table filter-features \
--i-table qiime_output_deblur/vsearch_99/clustered_table.qza \
--p-min-frequency 50 \
--o-filtered-table qiime_output_deblur/vsearch_99/clustered_table-ASV50.qza 

qiime feature-table filter-seqs \
--i-data qiime_output_deblur/vsearch_99/clustered_sequences.qza \
--i-table qiime_output_deblur/vsearch_99/clustered_table-ASV50.qza \
--o-filtered-data qiime_output_deblur/vsearch_99/clustered_repset-ASV50.qza

qiime tools export \
--input-path qiime_output_deblur/vsearch_99/clustered_table-ASV50.qza \
--output-path qiime_output_deblur/vsearch_99/clustered_table-ASV50
else 
echo "already filtered ASV50 for vsearch clustered 99"
fi



###### Do vsearch on dataset; try all at once #######
if [ ! -e qiime_output_deblur/vsearch_99/taxonomy_vsearch_80 ]
then
echo "Doing taxonomy vsearch 80"
qiime feature-classifier classify-consensus-vsearch \
--i-query qiime_output_deblur/vsearch_99/clustered_repset-ASV50.qza \
--i-reference-reads qiime_output_deblur/databases/silva-138-99-seqs-515-806.qza \
--i-reference-taxonomy qiime_output_deblur/databases/silva-138-99-tax-515-806.qza \
--p-perc-identity 0.8 \
--o-classification qiime_output_deblur/vsearch_99/taxonomy_vsearch_80

qiime tools export \
--input-path qiime_output_deblur/vsearch_99/taxonomy_vsearch_80.qza \
--output-path qiime_output_deblur/vsearch_99/taxonomy_vsearch_80
fi

# Now, filter ASVs by taxonomy

qiime taxa filter-table \
--i-table qiime_output_deblur/vsearch_99/clustered_table-ASV50.qza \
--i-taxonomy qiime_output_deblur/vsearch_99/taxonomy_vsearch_80.qza \
--p-exclude Unassigned,Chloroplast,Mitochondria,Eukaryote \
--o-filtered-table qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt.qza 

qiime feature-table filter-seqs \
--i-data qiime_output_deblur/vsearch_99/clustered_repset-ASV50.qza \
--i-table qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt.qza \
--o-filtered-data qiime_output_deblur/vsearch_99/clustered_repset-ASV50-taxfilt.qza 

qiime tools export \
--input-path qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt.qza \
--output-path qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt

biom convert --to-tsv \
-i qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt/feature-table.biom \
-o qiime_output_deblur/vsearch_99/clustered_table-ASV50-taxfilt/feature-table.txt

if [ ! -e qiime_output_deblur/vsearch_99/tree_99_ASV50.qza ]
then

qiime fragment-insertion sepp \
--i-representative-sequences qiime_output_deblur/vsearch_99/clustered_repset-ASV50-taxfilt.qza \
--i-reference-database qiime_output_deblur/databases/sepp-refs-silva-128.qza \
--o-tree qiime_output_deblur/vsearch_99/tree_99_ASV50.qza \
--o-placements qiime_output_deblur/vsearch_99/tree_placements_99_ASV50.qza

qiime tools export \
--input-path qiime_output_deblur/vsearch_99/tree_99_ASV50.qza \
--output-path qiime_output_deblur/vsearch_99/tree_99_ASV50

qiime tools export \
--input-path qiime_output_deblur/vsearch_99/tree_placements_99_ASV50.qza \
--output-path qiime_output_deblur/vsearch_99/tree_placements_99_ASV50
fi


#mkdir -p qiime_output_deblur/vsearch_99
#scp mchen@microbe.colorado.edu:/data/mchen/Project_metaanalysis/output_data/deblur/vsearch_99/* qiime_output_deblur/vsearch_99
#scp mchen@microbe.colorado.edu:/data/mchen/Project_metaanalysis/output_data/deblur/vsearch_99/* qiime_output_deblur/vsearch_99/

#bash ../02a_get_molecClockSeqs.sh

python ../02b_ultrametric_tree.py # NOTE: Can't be in qiime2 environment for this to run properly.


####### Create picrust tree ###########


wget http://kronos.pharmacology.dal.ca/public_files/picrust/picrust2_tutorial_files/picrust2_default_sepp_ref.qza \
-O qiime_output_deblur/databases/picrust2_default_sepp_ref.qza

if [ ! -e output_data/deblur/vsearch_99/picrust_placed ]
then
qiime fragment-insertion sepp \
--i-representative-sequences qiime_output_deblur/vsearch_99/clustered_repset-ASV50-taxfilt.qza \
--p-threads 1 \
--i-reference-database qiime_output_deblur/databases/picrust2_default_sepp_ref.qza \
--output-dir qiime_output_deblur/vsearch_99/picrust_placed

qiime tools export \
--input-path qiime_output_deblur/vsearch_99/picrust_placed/tree.qza \
--output-path qiime_output_deblur/vsearch_99/picrust_placed

fi

#scp -r mchen@microbe.colorado.edu:/data/mchen/Project_metaanalysis/output_data/deblur/vsearch_99/picrust_placed qiime_output_deblur/vsearch_99/


