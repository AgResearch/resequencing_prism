#!/bin/bash

declare -a files_array

function get_opts() {

   DRY_RUN=no
   DEBUG=no
   HPC_TYPE=slurm
   OUT_DIR=
   FORCE=no
   bed_file=
   help_text="
usage : 

extract_overlapping_reads.sh [-n for dry run] -b bed_file -O out_dir bam_file [bam_file . . . ]


"

   # defaults:
   while getopts ":nhdfO:C:b:" opt; do
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
         bed_file=$OPTARG
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

   if [ ! -f $bed_file  ]; then
      echo "no such file $bed_file"
      exit 1
   fi

}

function echo_opts() {
  echo OUT_DIR=$OUT_DIR
  echo DRY_RUN=$DRY_RUN
  echo DEBUG=$DEBUG
  echo HPC_TYPE=$HPC_TYPE
}

#
# edit this method to set required environment (or set up
# before running this script)
#
function configure_env() {

   # nothing here so far
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

   rm -f $OUT_DIR/command_file.txt

   for ((j=0;$j<$NUM_FILES;j=$j+2)) do
      bam_file=${files_array[$j]}

      file_base=`basename $bam_file`
      bed_base=`basename $bed_file`

      command="samtools view -b -L $bed_file --threads 4 -o $OUT_DIR/${file_base}_overlap_${bed_base}.bam $bam_file 1>$OUT_DIR/${file_base}_overlap_${bed_base}.stdout 2>${file_base}_overlap_${bed_base}.stderr"

      echo $command >> $OUT_DIR/command_file.txt

   done
}


function fake_prism() {
   echo "dry run ! "
   tardis -d $OUT_DIR --dry-run -c 1 source _condition_text_input_$OUT_DIR/command_file.txt > $OUT_DIR/command_file.log 2>&1
   echo "dry run !" 
   exit 0
}

function run_prism() {
   tardis -d $OUT_DIR -c 1 source _condition_text_input_$OUT_DIR/command_file.txt > $OUT_DIR/command_file.log 2>&1
}

function clean() {
   if [ $DEBUG == "no" ]; then
      rm -rf $OUT_DIR/tardis_*
   else 
      echo "debug mode, skipping clean"
   fi
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
         clean
      else
         echo "error state from run - skipping clean step"
         exit 1
      fi
   fi
}


set -x
main "$@"
set +x

