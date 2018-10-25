# resequencing_prism main makefile
#***************************************************************************************
# references:
#***************************************************************************************
# make: 
#     http://www.gnu.org/software/make/manual/make.html
#
##############################################
# how to make resequencing 
##############################################
%.archive: %.post_vcf
	$@.sh
	date > $@

%.post_vcf: %.vcf
	$@.sh
	date > $@

%.vcf: %.recal
	$@.sh
	date >$@

%.recal: %.merge
	$@.sh
	date >$@

%.merge: %.allbwa
	$@.sh
	date >$@

%.allbwa:
	$@.sh

%.bwa: 
	$@.sh
	date > $@


##############################################
# specify the intermediate files to keep 
##############################################
.PRECIOUS: %.log %.resequencing_prism 

##############################################
# cleaning - not yet doing this using make  
##############################################
clean:
	echo "no clean for now" 
