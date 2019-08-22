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
	$@.sh > $@.mk.log 2>&1
	date > $@

%.post_vcf: %.vcf
	$@.sh > $@.mk.log 2>&1
	date > $@

%.vcf: %.recal
	$@.sh > $@.mk.log 2>&1
	date >$@

%.recal: %.merge
	$@.sh > $@.mk.log 2>&1
	date >$@

%.merge: %.allbwa
	$@.sh > $@.mk.log 2>&1
	date >$@

%.allbwa:
	$@.sh > $@.mk.log 2>&1

%.bwa: 
	$@.sh > $@.mk.log 2>&1
	date > $@


##############################################
# specify the intermediate files to keep 
##############################################
.PRECIOUS: %.log %.bwa %.allbwa %.merge %.recal %.vcf %.post_vcf %.archive %.archive %.resequencing_prism 

##############################################
# cleaning - not yet doing this using make  
##############################################
clean:
	echo "no clean for now" 
