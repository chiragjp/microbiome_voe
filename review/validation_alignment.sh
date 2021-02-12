#!/bin/bash

#run alignment for validation

num_threads=1

name=$(echo "${1}" | cut -f3 -d/ | sed 's/.screened.adapter.screened.hg19.pair.2.fq.gz//g')
index=${3} 
fastafile=${4}
outputname="${name}"_output
cat $1 $2 > ${name}

bam_filename="${name}".catalog.bam
echo "Starting Alignment." >&2
diamond blastp --db ${index}  --query ${name} --outfmt 101 -b 6 -p ${num_threads} -o ${outputname}

cat ${outputname} | samtools view -T ${fastafile} -@ ${num_threads} -b -h -o ${bam_filename} -

rm ${outputname}
rm ${name}
#
# Cleaning up fastq's to make space
#echo 'Removing FastQ files' >&2
#rm ${read1_filename} ${read2_filename}

#
# Sort the bam file
echo 'Sorting the bam file' >&2
samtools 'sort' \
    -l 9 \
    -o ${bam_filename%.*}'.sorted.bam' \
    -O bam \
    -@ ${num_threads} \
    ${bam_filename}

#
# Cleaning up unsorted bam
echo 'Removing unsorted bam' >&2
rm ${bam_filename}
bam_filename=${bam_filename%.*}'.sorted.bam'

#
# Index the bam
echo 'Indexing the bam file' >&2
bam_index_filename=${bam_filename%.*}'.bai'
samtools 'index' -@ ${num_threads} -b ${bam_filename} ${bam_index_filename}

echo 'Extracting aligned reads'
countsfile=${name}_alignment_data.tsv
samtools idxstats -@ ${num_threads} ${bam_filename} > ${countsfile}