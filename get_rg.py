#!/usr/bin/env python

import os, re, argparse

def get_options():
    description = """
    """
    long_description = """
Get Read group string as per this spec : 

Map trimmed reads (pairs and singles that pass above QC) to reference using bwa mem (https://github.com/lh3/bwa) specifying read groups (see https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups for definitions) with the following options 
-R @RG\\tID:${RGID}\\tPL:${RGPL}\\tLB:${RGLB}\\tSM:${RGSM}
where RGID is usually the prefix for the fastq file, RGPL is the sequencing platform (ILLUMINA, SOLID or 454), RGLB is the library name and RGSM must be the international ID of the animal. Other read group tags can be populated but these are required. If your animal does not have an international ID, you should create one that conforms to Interbull standards, ie 3 character breed code + 3 character country code + sex code (M or F) + 12 character animal ID, eg HOLCANM000000352790. See http://www.interbull.org/ib/icarbreedcodes 
example : ./get_rg.py --subject AANNZLM020058096514 /dataset/gseq_processing/scratch/resequencing/1000_bulls_2018/AANNZLM020058096514/*.fastq.gz 

"""

    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file_names', type=str, nargs='+',metavar="filename", help='list of files to process')
    parser.add_argument('--subject', dest='subject', default=None, help="name of animal", required=True)
    parser.add_argument('--platform', dest='platform', default="ILLUMINA", help="name of platform")
    parser.add_argument('--library', dest='library', default=None, help="name of library")
    parser.add_argument('-t', dest='task', default="rg", choices=["rg","common_prefix"], help="task")
  
    args = vars(parser.parse_args())
    return args

        
    
def main():
    options = get_options()
    base_names = [ os.path.basename(file_name) for file_name in options["file_names"] ]
    options["moniker"] = os.path.commonprefix(base_names)
    if options["library"] is None:
        options["library"] = options["moniker"]

    
    rg="@RG\\\\tID:%(moniker)s\\\\tPL:%(platform)s\\\\tLB:%(library)s\\\\tSM:%(subject)s"%options

    if options["task"] == "rg":
       print(rg)
    else:
       print(options["moniker"])
        
                                
if __name__ == "__main__":
    main()
