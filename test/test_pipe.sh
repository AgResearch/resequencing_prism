#!/bin/bash
export SEQ_PRISMS_BIN=/dataset/gseq_processing/active/bin/resequencing_prism/seq_prisms
export RESEQUENCING_PRISM_BIN=/dataset/gseq_processing/active/bin/resequencing_prism


sample_name=TEST 
#C3PF2ACXX	5	C3PF2ACXX-1143-09-5-1	932597	TGACCA	NZGL01143
#1769106617	M	932597
#1769106617	M	ANGUS	AAN	10	00	001769106617	12	AANNZLM001769106617

R1=/dataset/gseq_processing/active/bin/resequencing_prism/test/test_R1.fastq.gz 
R2=/dataset/gseq_processing/active/bin/resequencing_prism/test/test_R2.fastq.gz

mkdir -p /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/test/$sample_name
#$RESEQUENCING_PRISM_BIN/resequencing_prism.sh -f -n -r /dataset/datacache/scratch/misc/genomes/ARS-UCD1.2_Btau5.0.1Y.fa -b /dataset/datacache/scratch/misc/indexes/bwa/ARS-UCD1.2_Btau5.0.1Y.fa -v /dataset/datacache/scratch/misc/genomes/ARS1.2PlusY_BQSR.vcf.gz -s $sample_name -O /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/test/$sample_name $R1 $R2  
#$RESEQUENCING_PRISM_BIN/resequencing_prism.sh -f  -r /dataset/datacache/scratch/misc/genomes/ARS-UCD1.2_Btau5.0.1Y.fa -b /dataset/datacache/scratch/misc/indexes/bwa/ARS-UCD1.2_Btau5.0.1Y.fa -v /dataset/datacache/scratch/misc/genomes/ARS1.2PlusY_BQSR.vcf.gz -s $sample_name -O /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/test/$sample_name $R1 $R2  

$RESEQUENCING_PRISM_BIN/resequencing_prism.sh -f -a archive -A /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/archivetest/TEST  -r /dataset/datacache/scratch/misc/genomes/ARS-UCD1.2_Btau5.0.1Y.fa -b /dataset/datacache/scratch/misc/indexes/bwa/ARS-UCD1.2_Btau5.0.1Y.fa -v /dataset/datacache/scratch/misc/genomes/ARS1.2PlusY_BQSR.vcf.gz -s $sample_name -O /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/test/$sample_name $R1 $R2  
