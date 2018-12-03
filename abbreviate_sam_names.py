#!/bin/env pypy

# abbreviate names in a SAM file to be less than MAX_LEN
# but still unique , by including a hash of the entire name in the 
# abbreviated name 
# (this gets around e.g. 
#"
#[E::sam_parse1] query name too long
#[W::sam_read1] Parse error at line 7736052
#[main_samview] truncated file.
#"
# usage example (e.g. to avoid writing out the patched sam) : 
# cat long_names.sam | ./abbreviate_sam_names.py  | samtools view -Sb -@ 12  -o shorter_names.bam   1>abbreviation.stdout 2>abbreviation.sdterr


MAX_LEN=250   # found 252 did not work  - still get parsing error ( ref https://www.biostars.org/p/259891/ )

import sys
import hashlib
import re

for record in sys.stdin:
   if len(record.strip()) ==0:
      sys.stdout.write(record)		# output whitespace records unchanged
   else:
      fields=re.split("\s+", record.strip())
      if fields[0] == "@":
         sys.stdout.write(record)  	# output header records unchanged
      else:
         if len(fields[0]) <= MAX_LEN:
            sys.stdout.write(record) 	# output records that don't need abbreviating unchanged
         else:
            h=hashlib.md5()
            h.update(fields[0])
            hash=h.hexdigest()
            new_name="%s_%s"%(fields[0][0:MAX_LEN-len(hash)-1], hash)  	# i.e. make a digest of the entire name and make it part of the shortened
									# name so its unique
            sys.stdout.write("\t".join([new_name] + fields[1:])+"\n")
