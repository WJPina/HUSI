#!/bin/bash
dataset="E-MTAB-5403"
cp RNA_seq_PE.config.yaml RNA_seq_PE_$dataset.config.yaml
echo "DATASET: \"$dataset\"" >> RNA_seq_PE_$dataset.config.yaml
echo "SAMPLES: " >> RNA_seq_PE_$dataset.config.yaml
while read sample_id;do
	sample_id=$(echo $sample_id | tr -d '\r')
  	echo -e "    ${sample_id}: " >> RNA_seq_PE_$dataset.config.yaml
	echo -e "        R1: /mnt/data3/wangj2/AgingScore/Datasets/Train/${dataset}/Fastq/${sample_id}_1.fastq.gz" >> RNA_seq_PE_$dataset.config.yaml
	echo -e "        R2: /mnt/data3/wangj2/AgingScore/Datasets/Train/${dataset}/Fastq/${sample_id}_2.fastq.gz"	>> RNA_seq_PE_$dataset.config.yaml
done < sampleList.txt
