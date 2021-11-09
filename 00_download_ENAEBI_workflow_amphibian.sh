#!bin/bash

###### Progress bar ########
BAR='####################'

#for i in {1..20}; do
#	let "tempPerc=($i*100/20)"
#	echo -ne "\r $tempPerc % [ ${BAR:0:$i}"
#	sleep 1
#done


#### Studies to download from ####

:>list_studies_to_include.txt
echo 'PRJEB35122'>> list_studies_to_include.txt 
echo 'PRJNA292303'>> list_studies_to_include.txt 
echo 'PRJNA292930'>> list_studies_to_include.txt 
echo 'PRJNA293463'>> list_studies_to_include.txt 
echo 'PRJNA299015'>> list_studies_to_include.txt 
echo 'PRJNA300286'>> list_studies_to_include.txt 
echo 'PRJNA300407'>> list_studies_to_include.txt 
echo 'PRJNA316224'>> list_studies_to_include.txt 
echo 'PRJNA320968'>> list_studies_to_include.txt 
echo 'PRJNA320969'>> list_studies_to_include.txt 
echo 'PRJNA320971'>> list_studies_to_include.txt 
echo 'PRJNA339746'>> list_studies_to_include.txt 
echo 'PRJNA368730'>> list_studies_to_include.txt 
echo 'PRJNA391810'>> list_studies_to_include.txt 
echo 'PRJNA393292'>> list_studies_to_include.txt 
echo 'PRJNA398139'>> list_studies_to_include.txt 
echo 'PRJNA430498'>> list_studies_to_include.txt 
echo 'PRJNA435631'>> list_studies_to_include.txt 
echo 'PRJNA437103'>> list_studies_to_include.txt 
echo 'PRJNA439189'>> list_studies_to_include.txt 
echo 'PRJNA477390'>> list_studies_to_include.txt 
echo 'PRJNA504463'>> list_studies_to_include.txt 
echo 'PRJNA504466'>> list_studies_to_include.txt 
echo 'PRJNA521543'>> list_studies_to_include.txt # F AND R
echo 'PRJNA549036'>> list_studies_to_include.txt 
echo 'PRJNA601697'>> list_studies_to_include.txt 
echo 'PRJNA603391'>> list_studies_to_include.txt 
# echo 'PRJNA613575'>> list_studies_to_include.txt This was hiseq
echo 'PRJNA641451'>> list_studies_to_include.txt 
echo 'PRJNA665482'>> list_studies_to_include.txt 
echo 'PRJNA720436'>> list_studies_to_include.txt 
echo 'PRJNA320970'>> list_studies_to_include.txt 
echo 'PRJNA326938'>> list_studies_to_include.txt 
echo 'PRJNA394790'>> list_studies_to_include.txt 
echo 'PRJNA474496'>> list_studies_to_include.txt 

#echo ''>> list_studies_to_include.txt  





## Also, got data from DRYAD

### #Get study accession numbers for samples #####
mkdir -p 00_study_acc

echo "LOADING STUDY ACCESSION LISTS"
while read study
do
echo $study

if [ -e 00_study_acc/$study.tsv ]
then
echo "Already downloaded $study"
else
wget -q "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${study}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,experiment_alias,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true" -O 00_study_acc/$study.tsv
fi

done < list_studies_to_include.txt

##### Remove any files that are NOT in list of studies ####
for f in 00_study_acc/*
do
studysearch=$(echo $f | sed 's/\.tsv//g' | sed 's/^.*\///g')

if grep -q $studysearch list_studies_to_include.txt
then
echo "keep $studysearch"
else
echo "remove $studysearch"
rm 00_study_acc/$studysearch.tsv
fi

done 

###### Get metadata xml files #######
mkdir -p 00_xml_meta

# sample accession is always column 2
echo 'GETTING ALL SAMPLE ACCESSIONS'
:>all_sample_accessions.txt
for filename in 00_study_acc/*.tsv; do
	echo ${filename}
	awk '{ print $2 }' $filename | grep -v 'sample_accession' >> all_sample_accessions.txt
done

#echo 'GETTING ALL RUN ACCESSIONS'
#:>all_run_accessions.txt
#for filename in 00_study_acc/*.tsv; do
#	echo ${filename}
#	awk '{ print $4 }' $filename | grep -v 'run_accession' >> all_run_accessions.txt
#done

# Look through all sample accessions and download xml file
nSamps=$(cat all_sample_accessions.txt | wc -l) 
currSamp=1
while read sampacc
do
	
let "tempPerc=($currSamp*100/$nSamps)"
let "endBar=($tempPerc*20/100)"
echo -ne "\r $tempPerc % [ ${BAR:0:$endBar}"

if [ -e "00_xml_meta/${sampacc}.xml" ]
then
	echo "$sampacc exists"
else
curl -s -C - -X POST "https://www.ebi.ac.uk/ena/browser/api/xml?accessions=${sampacc}&download=true" -H "accept: application/xml" -o 00_xml_meta/${sampacc}.xml
fi

let "currSamp++"

done < all_sample_accessions.txt



## Look through all run accessions and download xml file
#nSamps=$(cat all_run_accessions.txt | wc -l) 
#currSamp=1
#mkdir -p 00_xml_meta_run
#while read sampacc
#do
#	
#let "tempPerc=($currSamp*100/$nSamps)"
#let "endBar=($tempPerc*20/100)"
#echo -ne "\r $tempPerc % [ ${BAR:0:$endBar}"

#if [ -e "00_xml_meta_run/${sampacc}.xml" ]
#then
#	echo "$sampacc exists"
#else
#curl -s -C - -X POST "https://www.ebi.ac.uk/ena/browser/api/xml?accessions=${sampacc}&download=true" -H "accept: application/xml" -o 00_xml_meta_run/${sampacc}.xml
#fi

#let "currSamp++"

#done < all_run_accessions.txt


##### Remove any sample accessions that are NOT in list of studies ####
for f in 00_xml_meta/*
do
sampsearch=$(echo $f | sed 's/\.xml//g' | sed 's/^.*\///g')

if grep -q $sampsearch all_sample_accessions.txt
then
echo "keep $sampsearch"
else
rm 00_xml_meta/$sampsearch.xml
fi

done 



#curl -s -C - -X POST "https://www.ebi.ac.uk/ena/browser/api/xml?accessions=SAMEA5041711&download=true" -H "accept: application/xml" -o xml_meta/SAMEA5041711.xml

#### Convert XML to TSV ####

python3 code/xml_to_tsv.py -i 00_xml_meta -o all_sample_meta.txt

mkdir -p 00_all_sample_meta
mv all_sample_meta.txt ./00_all_sample_meta
mv all_sample_accessions.txt ./00_all_sample_meta
#### Filter desired samples #####



#### Get FASTA files ######





