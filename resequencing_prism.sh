#!/bin/bash

declare -a files_array

function get_opts() {

   DRY_RUN=no
   DEBUG=no
   HPC_TYPE=slurm
   OUT_DIR=
   MAX_TASKS=1
   FORCE=no
   ref_genome_index=
   ref_genome_sequence=
   variant_info=
   sample_name=
   gatk_jar=/dataset/gseq_processing/active/bin/resequencing_prism/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar

   OPTICAL_DUPLICATE_PIXEL_DISTANCE=100
   help_text="
\n
./resequencing_prism.sh  [-h (=help)] [-n (=dry run)] [-d (=debug - no clean up)] [-f (=overwrite any existing scritps)] [-C local|slurm] -s sample_namea -r ref_genome_index -b bwa_index -v variant_info -O outdir  input_R1 input_R2 [ input_R1 input_R2 . . . ] \n
\n
\n
examples:\n
# dry run - only generate scripts etc\n
./resequencing_prism.sh -n -r /dataset/datacache/scratch/misc/genomes/ARS-UCD1.2_Btau5.0.1Y.fa -b /dataset/datacache/scratch/misc/indexes/bwa/ARS-UCD1.2_Btau5.0.1Y.fa -v /dataset/datacache/scratch/misc/genomes/ARS1.2PlusY_BQSR.vcf.gz -s AANNZLM017427080180 -O /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM017427080180   /dataset/AG_1000_bulls/archive/nzgl01143/NZGL01143_C3P58ACXX/Raw/C3P58ACXX-1143-15-5-1_NoIndex_L003_R1_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01143/NZGL01143_C3P58ACXX/Raw/C3P58ACXX-1143-15-5-1_NoIndex_L003_R2_001.fastq.gz\n

# run everything on local machine (default is to schedule everything on slurm cluster) \n
./resequencing_prism.sh -C local -r /dataset/datacache/scratch/misc/genomes/ARS-UCD1.2_Btau5.0.1Y.fa -b /dataset/datacache/scratch/misc/indexes/bwa/ARS-UCD1.2_Btau5.0.1Y.fa -v /dataset/datacache/scratch/misc/genomes/ARS1.2PlusY_BQSR.vcf.gz -s AANNZLM017427080180 -O /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM017427080180   /dataset/AG_1000_bulls/archive/nzgl01143/NZGL01143_C3P58ACXX/Raw/C3P58ACXX-1143-15-5-1_NoIndex_L003_R1_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01143/NZGL01143_C3P58ACXX/Raw/C3P58ACXX-1143-15-5-1_NoIndex_L003_R2_001.fastq.gz\n

# process multiple pairs of input files 
./resequencing_prism.sh -r /dataset/datacache/scratch/misc/genomes/ARS-UCD1.2_Btau5.0.1Y.fa -b /dataset/datacache/scratch/misc/indexes/bwa/ARS-UCD1.2_Btau5.0.1Y.fa -v /dataset/datacache/scratch/misc/genomes/ARS1.2PlusY_BQSR.vcf.gz -s HERNZLM000281000146 -O /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/HERNZLM000281000146 /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9FHAADXX-1263-05-5-1_TAGCTT_L001_R1_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9FHAADXX-1263-05-5-1_TAGCTT_L001_R2_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9FHAADXX-1263-05-5-1_TAGCTT_L002_R1_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9FHAADXX-1263-05-5-1_TAGCTT_L002_R2_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9WBCADXX-1263-05-5-1_TAGCTT_L002_R1_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9WBCADXX-1263-05-5-1_TAGCTT_L002_R2_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9WUGADXX-1263-05-5-1_TAGCTT_L001_R1_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9WUGADXX-1263-05-5-1_TAGCTT_L001_R2_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9WUGADXX-1263-05-5-1_TAGCTT_L002_R1_001.fastq.gz /dataset/AG_1000_bulls/archive/nzgl01263_v2/Raw/H9WUGADXX-1263-05-5-1_TAGCTT_L002_R2_001.fastq.gz\n

\n
"

   # defaults:
   while getopts ":nhdfO:C:r:s:b:v:" opt; do
   case $opt in
       n)
         DRY_RUN=yes
         ;;
       d)
         DEBUG=yes
         ;;
       h)
         echo -e $help_text
         exit 0
         ;;
       f)
         FORCE=yes
         ;;
       O)
         OUT_DIR=$OPTARG
         ;;
       C)
         HPC_TYPE=$OPTARG
         ;;
       b)
         ref_genome_index=$OPTARG
         ;;
       r)
         ref_genome_sequence=$OPTARG
         ;;
       v)
         variant_info=$OPTARG
         ;;
       s)
         sample_name=$OPTARG
         ;;
       \?)
         echo "Invalid option: -$OPTARG" >&2
         exit 1
         ;;
       :)
         echo "Option -$OPTARG requires an argument." >&2
         exit 1
         ;;
     esac
   done

   shift $((OPTIND-1))

   FILE_STRING=$@

   # this is needed because of the way we process args a "$@" - which 
   # is needed in order to parse parameter sets to be passed to the 
   # aligner (which are space-separated)
   declare -a files="(${FILE_STRING})";
   NUM_FILES=${#files[*]}
   for ((i=0;$i<$NUM_FILES;i=$i+1)) do
      files_array[$i]=${files[$i]}     
   done
}


function check_opts() {
   if [  -z "$OUT_DIR" ]; then
      echo "must specify OUT_DIR ( -O )"
      exit 1
   fi
   if [ ! -d $OUT_DIR ]; then
      echo "OUT_DIR $OUT_DIR not found"
      exit 1
   fi
   if [[ $HPC_TYPE != "local" && $HPC_TYPE != "slurm" ]]; then
      echo "HPC_TYPE must be one of local, slurm"
      exit 1
   fi
   if [ ! -f ${ref_genome_index}.bwt  ]; then
      echo "bad index  (cant see ${ref_genome_index}.bwt ) (you might need to supply the full path ?)"
      exit 1
   fi
   if [ ! -f $ref_genome_sequence  ]; then
      echo "no such file $ref_genome_sequence"
      exit 1
   fi
   if [ ! -f $variant_info  ]; then
      echo "no such file $variant_info"
      exit 1
   fi

}

function echo_opts() {
  echo OUT_DIR=$OUT_DIR
  echo DRY_RUN=$DRY_RUN
  echo DEBUG=$DEBUG
  echo HPC_TYPE=$HPC_TYPE
  echo sample_name=$sample_name
  echo ref_genome_index=$ref_genome_index
}

#
# edit this method to set required environment (or set up
# before running this script)
#
function configure_env() {
   export CONDA_ENVS_PATH=$CONDA_ENVS_PATH:/dataset/bioinformatics_dev/active/conda-env

   cd $RESEQUENCING_PRISM_BIN
   cp ./resequencing_prism.sh $OUT_DIR
   cp ./resequencing_prism.mk $OUT_DIR
   cp ./get_rg.py  $OUT_DIR
   cp ./etc/tardis.toml.4g $OUT_DIR/tardis.toml    # this version points at slurm batch template with 4G spec - good for bwa/ trim step
   echo "
conda activate /dataset/gseq_processing/active/bin/resequencing_prism/conda/resequencing_prism
" > $OUT_DIR/resequencing_prism_env.src

   cd $OUT_DIR
}


function check_env() {
   if [ -z "$SEQ_PRISMS_BIN" ]; then
      echo "SEQ_PRISMS_BIN not set - exiting"
      exit 1
   fi
   if [ -z "$RESEQUENCING_PRISM_BIN" ]; then
      echo "RESEQUENCING_PRISM_BIN not set - exiting"
      exit 1
   fi
}

function get_targets() {

   rm -f $OUT_DIR/resequencing_targets.txt
   rm -f $OUT_DIR/input_file_list.txt
   picard_merge_input_string=""

   for ((j=0;$j<$NUM_FILES;j=$j+2)) do
      R1=${files_array[$j]}
      R2=${files_array[$j+1]}

      if [[ ( ! -f $R1 ) || ( ! -f $R2 ) ]]; then
         echo "could not find either file $R1 or file $R2"
         exit 1
      fi

      echo $R1 $R2 >> $OUT_DIR/input_file_list.txt

      file_base=`basename $R1`
      parameters_moniker=`basename $sample_name`
      parameters_moniker=${parameters_moniker}.`basename $ref_genome_index`

      moniker=${file_base}.${parameters_moniker}
      echo $OUT_DIR/${moniker}.resequencing_prism >> $OUT_DIR/resequencing_targets.txt

      prefix=`./get_rg.py --subject $sample_name  -t common_prefix $R1 $R2`
      picard_merge_string="$picard_merge_string I=${prefix}_paired.sorted.bam I=${prefix}_1.sorted.bam I=${prefix}_2.sorted.bam"

      # generate wrapper alignment script 
      script_filename=$OUT_DIR/${moniker}.sh

      if [ -f $script_filename ]; then
         if [ ! $FORCE == yes ]; then
            echo "found existing script $script_filename - will re-use (use -f to force rebuild ) "
            continue
         fi
      fi


      echo "#!/bin/bash
# set up shortcuts to paired files 
cp -s $R1  $OUT_DIR
cp -s $R2  $OUT_DIR
cd $OUT_DIR

# get read-goup info 
RG=\`./get_rg.py --subject $sample_name $R1 $R2 \`
cp $RESEQUENCING_PRISM_BIN/etc/env.inc .
cp $RESEQUENCING_PRISM_BIN/etc/TruSeq3-PE.fa .

R1link=$OUT_DIR/\`basename $R1\`
R2link=$OUT_DIR/\`basename $R2\`
summary=$OUT_DIR/$prefix.trim_summary
R1_trimmed=$OUT_DIR/${prefix}1.trimmed.fastq
R2_trimmed=$OUT_DIR/${prefix}2.trimmed.fastq
R1_single=$OUT_DIR/${prefix}1.single.fastq
R2_single=$OUT_DIR/${prefix}2.single.fastq
paired_sam=$OUT_DIR/${prefix}_paired.sam
single1_sam=$OUT_DIR/${prefix}_1.sam
single2_sam=$OUT_DIR/${prefix}_2.sam

time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc trimmomatic PE -threads 8 -summary _condition_uncompressedtext_output_\$summary _condition_fastq_input_\$R1link _condition_fastq_input_\$R2link _condition_throughput_\$R1_trimmed  _condition_throughput_\$R1_single  _condition_throughput_\$R2_trimmed _condition_throughput_\$R2_single ILLUMINACLIP:$OUT_DIR/TruSeq3-PE.fa:2:30:3:1:true LEADING:20 TRAILING:20 SLIDINGWINDOW:3:15 AVGQUAL:20 MINLEN:35 \; bwa mem -M -t 8 -R \$RG $ref_genome_index _condition_throughput_\$R1_trimmed  _condition_throughput_\$R2_trimmed  \> _condition_uncompressedsam_output_\$paired_sam \;  bwa mem -M -t 8 -R \$RG $ref_genome_index _condition_throughput_\$R1_single \> _condition_uncompressedsam_output_\$single1_sam \; bwa mem -M -t 8 -R \$RG $ref_genome_index _condition_throughput_\$R2_single  \> _condition_uncompressedsam_output_\$single2_sam


for samfile in \$paired_sam \$single1_sam \$single2_sam; do
   moniker=\`basename \$samfile .sam\`
   sorted_bamfile=\${moniker}.sorted.bam
   time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc samtools sort -o $OUT_DIR/\$sorted_bamfile -O BAM \$samfile
   time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc samtools index $OUT_DIR/\$sorted_bamfile
done
" > $script_filename
      chmod +x $script_filename 
   done    # end loop generating alignment script per pair


   ######## generate global sripts which use all output bams ##########
   prefix=`./get_rg.py --subject $sample_name  -t common_prefix \`cat $OUT_DIR/input_file_list.txt\``
   moniker="allfiles_$prefix"


   # generate merge and mark duplicates script
   merge_script_filename=$OUT_DIR/${moniker}.merge.sh

   if [ -f $merge_script_filename ]; then
      if [ ! $FORCE == yes ]; then
         echo "found existing merge script $merge_script_filename - will re-use (use -f to force rebuild ) "
      else 
         rm -f $merge_script_filename
      fi 
   fi

   if [ ! -f $merge_script_filename ]; then
      echo "#!/bin/bash
cd $OUT_DIR
cp $RESEQUENCING_PRISM_BIN/etc/tardis.toml $OUT_DIR/tardis.toml  # this version points at the larger memory slurm template, needed for merge steps
time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc picard MergeSamFiles $picard_merge_string O=${sample_name}.sorted.bam VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true

# either : 
#time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc picard MarkDuplicates TMP_DIR=$OUT_DIR I=${sample_name}.sorted.bam O=${sample_name}_dedup.bam M=${sample_name}_dedup.metrics OPTICAL_DUPLICATE_PIXEL_DISTANCE=$OPTICAL_DUPLICATE_PIXEL_DISTANCE CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

# or : (as above tends to memory fault)  
PICARD_JAR=/dataset/gseq_processing/active/bin/resequencing_prism/conda/resequencing_prism/share/picard-2.18.2-0/picard.jar
time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc java -Xmx80G -jar \$PICARD_JAR  MarkDuplicates TMP_DIR=$OUT_DIR I=${sample_name}.sorted.bam O=${sample_name}_dedup.bam M=${sample_name}_dedup.metrics OPTICAL_DUPLICATE_PIXEL_DISTANCE=$OPTICAL_DUPLICATE_PIXEL_DISTANCE CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

" > $merge_script_filename
      chmod +x $merge_script_filename
   fi



   # generate recalibrate script
   recal_script_filename=$OUT_DIR/${moniker}.recal.sh

   if [ -f $recal_script_filename ]; then
      if [ ! $FORCE == yes ]; then
         echo "found existing recal script $recal_script_filename - will re-use (use -f to force rebuild ) "
      else
         rm $recal_script_filename
      fi
   fi

   if [ ! -f $recal_script_filename ]; then
      echo "#!/bin/bash
cd $OUT_DIR
cp $RESEQUENCING_PRISM_BIN/etc/tardis.toml $OUT_DIR/tardis.toml  # this version points at the larger memory slurm template, needed for merge steps
time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc java -Xmx80G -jar $gatk_jar  -T BaseRecalibrator -nct 8 -R $ref_genome_sequence  -I ${sample_name}_dedup.bam -knownSites:vcf $variant_info -o ${sample_name}.recal.table

time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc java -Xmx80G -jar $gatk_jar -T PrintReads -nct 8 -R $ref_genome_sequence -I ${sample_name}_dedup.bam  -BQSR ${sample_name}.recal.table -o ${sample_name}_dedup_recal.bam

" > $recal_script_filename
      chmod +x $recal_script_filename
   fi


   # generate vcf script
   vcf_script_filename=$OUT_DIR/${moniker}.vcf.sh

   if [ -f $vcf_script_filename ]; then
      if [ ! $FORCE == yes ]; then
         echo "found existing vcf script $vcf_script_filename - will re-use (use -f to force rebuild ) "
      else 
         rm $vcf_script_filename
      fi
   fi

   if [ ! -f $vcf_script_filename ]; then
      echo "#!/bin/bash
cd $OUT_DIR
cp $RESEQUENCING_PRISM_BIN/etc/tardis.toml $OUT_DIR/tardis.toml  # this version points at the larger memory slurm template, needed for merge steps
time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc java -Xmx80G -jar $gatk_jar -T HaplotypeCaller -nct 8 -R $ref_genome_sequence -I ${sample_name}_dedup_recal.bam -o ${sample_name}_dedup_recal.g.vcf.gz -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000


" > $vcf_script_filename
      chmod +x $vcf_script_filename
   fi

   # generate post_vcf script
   post_vcf_script_filename=$OUT_DIR/${moniker}.post_vcf.sh

   if [ -f $post_vcf_script_filename ]; then
      if [ ! $FORCE == yes ]; then
         echo "found existing post_vcf script $post_vcf_script_filename - will re-use (use -f to force rebuild ) "
      else
         rm $post_vcf_script_filename
      fi
   fi

   if [ ! -f $post_vcf_script_filename ]; then
      echo "#!/bin/bash
cd $OUT_DIR
# coverage 
cp $RESEQUENCING_PRISM_BIN/etc/tardis.toml $OUT_DIR/tardis.toml  # this version points at the larger memory slurm template, needed for merge steps
time tardis --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file env.inc java -Xmx80G -jar $gatk_jar -T DepthOfCoverage -R $ref_genome_sequence -I ${sample_name}_dedup_recal.bam --omitDepthOutputAtEachBase --logging_level ERROR --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 80 --summaryCoverageThreshold 90 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 --minBaseQuality 15 --minMappingQuality 30 --start 1 --stop 1000 --nBins 999 -dt NONE -o ${sample_name}_dedup_recal.coverage

# checksums
for file in *${sample_name}* ; do
   md5sum -b \$file >> checksums.txt
done
   
" > $post_vcf_script_filename
      chmod +x $post_vcf_script_filename
   fi
}


function fake_prism() {
   echo "dry run ! "
   make -n -f resequencing_prism.mk -d -k  --no-builtin-rules -j 16 `cat $OUT_DIR/resequencing_targets.txt` > $OUT_DIR/resequencing_prism.log 2>&1
   echo "dry run : summary commands are 
   "
   exit 0
}

function run_prism() {
   # this prepares each file
   make -f resequencing_prism.mk -d -k  --no-builtin-rules -j 16 `cat $OUT_DIR/resequencing_targets.txt` > $OUT_DIR/resequencing_prism.log 2>&1
}

function clean() {
   if [ $DEBUG == "no" ]; then
      rm -rf $OUT_DIR/tardis_*
      rm $OUT_DIR/*.sam
   else 
      echo "debug mode, skipping clean"
   fi
}


function merge_and_gatk_prism() {
   cp $RESEQUENCING_PRISM_BIN/etc/tardis.toml $OUT_DIR/tardis.toml  # this version points at the larger memory slurm template, needed for merge steps
   $merge_script_filename > $OUT_DIR/merge.log 2>&1
   $recal_script_filename > $OUT_DIR/recal.log 2>&1
   $vcf_script_filename > $OUT_DIR/vcf.log 2>&1
}


function main() {
   get_opts "$@"
   check_opts
   echo_opts
   check_env
   configure_env
   get_targets
   if [ $DRY_RUN != "no" ]; then
      fake_prism
   else
      run_prism
      if [ $? == 0 ] ; then
         merge_and_gatk_prism
         if [ $? == 0 ] ; then
            clean
         fi
      else
         echo "error state from resequencing run - skipping clean and html page generation"
         exit 1
      fi
   fi
}


set -x
main "$@"
set +x

