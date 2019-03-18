# $1 = annovar input file name (eg. FILE.avinput)
# $2 = PED file name (eg. FAMILY.ped)

#module load R
PSAP_PATH="/share/terra/nz2274/Application/PSAP_cluster/" #/home/local/ARCS/nz2274/Application/ #INSERT PATH TO PSAP DIRECTORY HERE eg. "/scratch/dclab/"
echo "start"$1"\t"$2
FB=$(echo $1|sed 's/\.avinput//g')
echo $FB."\t".$1."\n";
if [ $# == 2 ]
then  
# annotate with popstat (requires ped file)
	echo "Annotating data with population model"
	echo "Rscript ${PSAP_PATH}psap/RScripts/generic_apply_popStat.R ${1}.avinput ${2} $PSAP_PATH"

	if [ ! -s ${FB}_popStat.txt ]; then
		Rscript ${PSAP_PATH}psap/RScripts/generic_apply_popStat.R ${1}.avinput ${2} $PSAP_PATH	
	fi 
	if [ ! -s ${FB}_popStat.txt ]; then
 		echo  ${FB}"_popStat.txt is not generated!\n"
		exit
	fi 
# perform family analysis and validate
	echo "Performing family analysis and candidate validation"
	echo "Rscript ${PSAP_PATH}psap/RScripts/generic_candidate_analysis.R ${1}.avinput ${2} $PSAP_PATH"
#	if [ ! -s $FB.report.txt ]; then
		Rscript ${PSAP_PATH}psap/RScripts/generic_candidate_analysis.R ${1}.avinput ${2} $PSAP_PATH
#	fi
	echo "Done.  All resutls can be found in the path ${PWD}/annotated/"
else
	echo $#
	echo "Incorrect number of arguments"
	echo "Please provide ANNOVAR input file name and family pedigree file"
fi
