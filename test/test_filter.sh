#!/bin/bash

REF_GENOME=/dataset/datacache/scratch/misc/indexes/bwa/ARS-UCD1.2_Btau5.0.1Y.fa
OUT_DIR=/dataset/gseq_processing/active/bin/resequencing_prism/test

# test real data 
set -x
RG=`../get_rg.py --subject AANNZLM020058096514 /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/*.fastq.gz`
prefix=`../get_rg.py --subject AANNZLM020058096514 -t common_prefix /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/*.fastq.gz`
OUT_DIR=/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514
#cp env.inc $OUT_DIR
#cp /dataset/gseq_processing/active/bin/resequencing_prism/etc/TruSeq3-PE.fa $OUT_DIR
#cp -s /dataset/AG_1000_bulls/archive/nzgl01143/NZGL01143_C3P58ACXX/Raw/C3P58ACXX-1143-13-5-1_NoIndex_L001_R1_001.fastq.gz $OUT_DIR
#cp -s /dataset/AG_1000_bulls/archive/nzgl01143/NZGL01143_C3P58ACXX/Raw/C3P58ACXX-1143-13-5-1_NoIndex_L001_R2_001.fastq.gz $OUT_DIR
R1=$OUT_DIR/C3P58ACXX-1143-13-5-1_NoIndex_L001_R1_001.fastq.gz
R2=$OUT_DIR/C3P58ACXX-1143-13-5-1_NoIndex_L001_R2_001.fastq.gz
summary=$OUT_DIR/${prefix}.trim_summary
R1_trimmed=$OUT_DIR/${prefix}1.trimmed.fastq
R2_trimmed=$OUT_DIR/${prefix}2.trimmed.fastq
R1_single=$OUT_DIR/${prefix}1.single.fastq
R2_single=$OUT_DIR/${prefix}2.single.fastq
paired_sam=$OUT_DIR/${prefix}_paired.sam
single1_sam=$OUT_DIR/${prefix}_1.sam
single2_sam=$OUT_DIR/${prefix}_2.sam

#time tardis -d $OUT_DIR -c 2000000 --shell-include-file env.inc trimmomatic PE -threads 8 -summary _condition_uncompressedtext_output_$summary _condition_fastq_input_$R1 _condition_fastq_input_$R2 _condition_throughput_$R1_trimmed  _condition_throughput_$R1_single  _condition_throughput_$R2_trimmed _condition_throughput_$R2_single ILLUMINACLIP:$OUT_DIR/TruSeq3-PE.fa:2:30:3:1:true LEADING:20 TRAILING:20 SLIDINGWINDOW:3:15 AVGQUAL:20 MINLEN:35 \; bwa mem -M -t 8 -R $RG $REF_GENOME _condition_throughput_$R1_trimmed  _condition_throughput_$R2_trimmed  \> _condition_uncompressedsam_output_$paired_sam \;  bwa mem -M -t 8 -R $RG $REF_GENOME _condition_throughput_$R1_single \> _condition_uncompressedsam_output_$single1_sam \; bwa mem -M -t 8 -R $RG $REF_GENOME _condition_throughput_$R2_single  \> _condition_uncompressedsam_output_$single2_sam
#real    386m2.120s
#user    0m39.970s
#sys     0m18.179s


#samtools sort -o ${OutputFile}-pe.sorted.sam -O BAM ${OutputFile}-pe.bam
#samtools index ${OutputFile}-pe.sorted.bam

#for samfile in $paired_sam $single1_sam $single2_sam; do
cd $OUT_DIR
#for samfile in $single2_sam; do
#for samfile in $paired_sam $single1_sam; do
#   moniker=`basename $samfile .sam`
#   sorted_bamfile=${moniker}.sorted.bam
#   time tardis -d $OUT_DIR --shell-include-file env.inc samtools sort -o $OUT_DIR/$sorted_bamfile -O BAM $samfile 
#   time tardis -d $OUT_DIR --shell-include-file env.inc samtools index $OUT_DIR/$sorted_bamfile
#done
#BAMlist="C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sorted.bam C3P58ACXX-1143-13-5-1_NoIndex_L001_R_2.sorted.bam C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sorted.bam"
INTERNATIONALID=AANNZLM020058096514
#time picard MergeSamFiles I=C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sorted.bam I=C3P58ACXX-1143-13-5-1_NoIndex_L001_R_2.sorted.bam I=C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sorted.bam  O=${INTERNATIONALID}.sorted.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true

#time picard MarkDuplicates I=${INTERNATIONALID}.sorted.bam O=${INTERNATIONALID}_dedup.bam M=${INTERNATIONALID}_dedup.metrics OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT > test_rmdup.log 2>&1


GATK=/dataset/gseq_processing/active/bin/resequencing_prism/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
VCF=/dataset/datacache/scratch/misc/genomes/ARS1.2PlusY_BQSR.vcf.gz
REF=/dataset/datacache/scratch/misc/genomes/ARS-UCD1.2_Btau5.0.1Y.fa

#picard CreateSequenceDictionary R=/dataset/datacache/scratch/misc/genomes/ARS-UCD1.2_Btau5.0.1Y.fa  O=/dataset/datacache/scratch/misc/genomes/ARS-UCD1.2_Btau5.0.1Y.dict 

#java -Xmx80G -jar $GATK -T BaseRecalibrator -nct 8 -R $REF -I ${INTERNATIONALID}_dedup.bam -knownSites:vcf ${VCF} -o ${INTERNATIONALID}.recal.table
#java -Xmx80G -jar $GATK -T PrintReads -nct 8 -R $REF -I ${INTERNATIONALID}_dedup.bam -BQSR ${INTERNATIONALID}.recal.table -o ${INTERNATIONALID}_dedup_recal.bam



#java -Xmx80G -jar $GATK -T HaplotypeCaller -nct 8 -R $REF -I ${INTERNATIONALID}_dedup_recal.bam -o ${INTERNATIONALID}_dedup_recal.g.vcf.gz -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 

java -Xmx80G -jar $GATK -T DepthOfCoverage -R $REF  -I ${INTERNATIONALID}_dedup_recal.bam --omitDepthOutputAtEachBase --logging_level ERROR --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 80 --summaryCoverageThreshold 90 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 --minBaseQuality 15 --minMappingQuality 30 --start 1 --stop 1000 --nBins 999 -dt NONE -o ${INTERNATIONALID}_dedup_recal.coverage

#java -Xmx80G -jar $GATK -T DepthOfCoverage -R ARS-UCD1.2_Btau5.0.1Y.fa -I ${INTERNATIONALID}_dedup_recal.bam --omitDepthOutputAtEachBase --logging_level ERROR --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 80 --summaryCoverageThreshold 90 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 --minBaseQuality 15 --minMappingQuality 30 --start 1 --stop 1000 --nBins 999 -dt NONE -o ${INTERNATIONALID}_dedup_recal.coverage



#samtools sort -o ${OutputFile}-pe.sorted.sam -O BAM ${OutputFile}-pe.bam
#samtools index ${OutputFile}-pe.sorted.bam


exit



----------------- samtools --------------------
+ R2_single=/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R2.single.fastq
+ paired_sam=/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sam
+ single1_sam=/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sam
+ single2_sam=/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_2.sam
+ cd /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514
+ for samfile in '$paired_sam' '$single1_sam'
++ basename /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sam .sam
+ moniker=C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired
+ sorted_bamfile=C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sorted.bam
+ tardis -d /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514 --shell-include-file env.inc samtools sort -o /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sorted.bam -O BAM /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sam
reading config from /etc/tardis/tardis.toml
tool args = ['samtools', 'sort', '-o', '/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sorted.bam', '-O', 'BAM', '/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sam']
tardis.py : logging this session to /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/tardis_nBJ3M8
[bam_sort_core] merging from 157 files and 1 in-memory blocks...
tardis.py : done logging this session to /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/tardis_nBJ3M8 , no errors detected

real    74m50.018s
user    0m0.775s
sys     0m0.166s
+ tardis -d /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514 --shell-include-file env.inc samtools index /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sorted.bam
reading config from /etc/tardis/tardis.toml
tool args = ['samtools', 'index', '/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_paired.sorted.bam']
tardis.py : logging this session to /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/tardis_CFo4wY
tardis.py : done logging this session to /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/tardis_CFo4wY , no errors detected

real    5m5.404s
user    0m0.080s
sys     0m0.028s
+ for samfile in '$paired_sam' '$single1_sam'
++ basename /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sam .sam
+ moniker=C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1
+ sorted_bamfile=C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sorted.bam
+ tardis -d /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514 --shell-include-file env.inc samtools sort -o /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sorted.bam -O BAM /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sam
reading config from /etc/tardis/tardis.toml
tool args = ['samtools', 'sort', '-o', '/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sorted.bam', '-O', 'BAM', '/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sam']
tardis.py : logging this session to /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/tardis_ewdCFW
[bam_sort_core] merging from 4 files and 1 in-memory blocks...
tardis.py : done logging this session to /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/tardis_ewdCFW , no errors detected

real    2m15.226s
user    0m0.048s
sys     0m0.024s
+ tardis -d /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514 --shell-include-file env.inc samtools index /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sorted.bam
reading config from /etc/tardis/tardis.toml
tool args = ['samtools', 'index', '/dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/C3P58ACXX-1143-13-5-1_NoIndex_L001_R_1.sorted.bam']
tardis.py : logging this session to /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/tardis_3DL_hD
tardis.py : done logging this session to /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/tardis_3DL_hD , no errors detected

real    0m30.114s
user    0m0.036s
sys     0m0.017s
(/dataset/gseq_processing/active/bin/resequencing_prism/conda/resequencing_prism) iramohio-01$




