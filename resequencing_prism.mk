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
%.resequencing_prism: %.gatk_gvcf
	$*.sh
	date > $@

%.gatk_gvcf: %.gatk_analyse 
	$*.sh
	date > $@

%.gatk_analyse: %.gatk_print 
	$*.sh
	date > $@

%.gatk_print: %.gatk_recal 
	$*.sh
	date > $@

%.gatk_recal: %.mark_duplicates  
	$*.sh
	date > $@

%.mark_duplicates: %.map   
	$*.sh
	date > $@

%.map: %.trim   
	$*.sh
	date > $@

%.fastqc:    
	$*.sh
	date > $@

##############################################
# specify the intermediate files to keep 
##############################################
.PRECIOUS: %.log %.resequencing_prism %.gatk_gvcf %.gatk_analyse %.gatk_print %.gatk_recal %.mark_duplicates  %.map %.trim

##############################################
# cleaning - not yet doing this using make  
##############################################
clean:
	echo "no clean for now" 

