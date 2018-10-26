#cd /home/local/ARCS/nz2274/Freeze8/
echo "input: vcf.gz pedigree"
if [ "$#" -lt 2 ];then
   echo "input:vcf.gz pedigree"
   exit
fi
fvcf=$1 #"/home/local/ARCS/nz2274/Freeze8/Columbia_Freeze_Eight.NF.pVCF.vcf.gz"
fped=$2 #"/home/local/ARCS/nz2274/Freeze8/RGN8.sampleHas.ped"
script_folder=$PSAP_PATH #/home/nz2274/Application/PSAP_cluster/  #"/home/local/ARCS/nz2274/Freeze8/script/"





#ANNOVAR_PATH=$ANNOVAR #/home/local/ARCS/nz2274/Application/annovar/; #/home/local/users/jw/software_packages/annovar/ #INSERT PATH TO ANNOVAR DIRECTORY HERE eg. /scratch/dclab/annovar/
#ANNOVAR_DB=$ANNHDB #/home/local/ARCS/nz2274/Application/annovar/
echo  "input example: file.vcf.gz file.ped [1] [$liftover/hg38ToHg19.nochr.over.chain.gz]"
if [ "$#" -lt 2 ]; then
	echo  "input example: file.vcf.gz file.ped [1] [$liftover/hg38ToHg19.nochr.over.chain.gz]"
	exit;
fi

if [ "$#" -eq 3 ]; then
	
	echo "please specify converttohg19  chain_file"
	echo "input examples: file.vcf.gz file.ped 1  $liftover/hg38ToHg19.nochr.over.chain.gz"
	exit;
fi
if [ "$#" -gt 3 ]; then
	
convert=$3 #1  ## 1: convert vcf to hg19 0: no convertion
fchain=$4 #"~/Application/liftover/hg38ToHg19.nochr.over.chain.gz"

fi

fname=$(basename $fvcf)
fb=${fname%.*}
cb=$(basename $chain|sed 's/.chain.gz//g')
mkdir annotated
fmerge="annotated/$fb.hg19.avinput.hg19_multianno.txt"
function lift_over ()
{

	
	vcf=$1
	chain=$2  ## "/home/local/ARCS/nz2274/Application/liftover/hg19ToHg38.over.nochr.chain.gz"
	i=$3

	
	#for i in {1..22} X Y
	#do  
	f=$fb.$cb.$i.hg19.vcf;
	tabix $vcf $i -h |python ~/Pipeline/scripts/lift_over_chain.py  --format vcf --chain $chain  > $f  2> nohup.$fbname.$i.log  
	annotation $f $fb.$cb.$i.hg19

}

function annotation ()
{

   f=$1
	OUTFILE=$2
	perl ${ANNOVAR}convert2annovar.pl -format vcf4old $f -outfile annotated/${OUTFILE}.avinput -includeinfo 2> annotated/$OUTFILE.convert.log
	echo "perl ${ANNOVAR}table_annovar.pl ${OUTFILE}.avinput -remove -outfile annotated/${OUTFILE}.avinput ${ANNHDB}/ -buildver hg19 -protocol wgEncodeGencodeBasicV19,mac63kFreq_ALL,esp6500si_all,1000g2014sep_all,snp137,cadd -operation g,f,f,f,f,f -nastring NA -otherinfo -argument -separate,,,,,-otherinfo  > nohup.$i.anno.log"
	#
	#exit
	perl ${ANNOVAR}table_annovar.pl ${OUTFILE}.avinput -remove -outfile annotated/${OUTFILE}.avinput ${ANNHDB}/ -buildver hg19 -protocol wgEncodeGencodeBasicV19,mac63kFreq_ALL,esp6500si_all,1000g2014sep_all,snp137,cadd -operation g,f,f,f,f,f -nastring NA -otherinfo -argument -separate,,,,,-otherinfo  > nohup.$i.anno.log


}

if [ $convert eq 1 ]; then
	for i in {1..22} X Y; do
		echo $i
		lift_over $fvcf $fchain $i &
		 ## lift over hg38 to hg19 by per chromosome
		 ## and annotated each converted hg19 chromosome
	done
 
	wait;
    

	egrep '^chr|^C' -m 1  annotated/$fb.$cb.1.avinput.hg19_multianno.txt  > $fmerge #annotated/Columbia_Freeze_Eight.avinput.hg19_multianno.txt
	tabix $fvcf -H |egrep '^#C' >>  $fmerge  #annotated/Columbia_Freeze_Eight.avinput.hg19_multianno.txt
	#sed ':a;N;$!ba;s/\n//g'  annotated/Columbia_Freeze_Eight.avinput.hg19_multianno.txt 
	sed ':a;N;$!ba;s/\n//g'  $fmerge |sed 's/Otherinfo//g' > temp.anno
	mv temp.anno $fmerge
	 for i in {1..22} X Y; do 
	 	awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1){next} if($6 ~/^exonic|^spli/ && $9 !~/^synony/ && ($11 <0.01||$11=="."||$11=="NA")){print}}'   annotated/$fb.$cb.$i.avinput.hg19_multianno.txt  >> $fmerge
	 done 
 
	 ## sort file
	 less $fmerge |vcf-sort -c > $fb.tem
	 egrep '^C'  $fb.tem -m 1 > $fmerge #annotated/Columbia_Freeze_Eight.hg19.avinput.hg19_multianno.txt 
	 egrep -v '^C'  $fb.tem  >> $fmerge
 
 else
	 
	annotation $fvcf $fb
	 
 	egrep '^chr|^C' -m 1  annotated/$fb.avinput.hg19_multianno.txt  > $fmerge #annotated/Columbia_Freeze_Eight.avinput.hg19_multianno.txt
 	tabix $fvcf -H |egrep '^#C' >>  $fmerge  #annotated/Columbia_Freeze_Eight.avinput.hg19_multianno.txt
 	#sed ':a;N;$!ba;s/\n//g'  annotated/Columbia_Freeze_Eight.avinput.hg19_multianno.txt 
 	sed ':a;N;$!ba;s/\n//g'  $fmerge |sed 's/Otherinfo//g' > $fb.temp.anno
 	mv $fb.temp.anno $fmerge
 	awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1){next} if($6 ~/^exonic|^spli/ && $9 !~/^synony/ && ($11 <0.01||$11=="."||$11=="NA")){print}}'   annotated/$fb.avinput.hg19_multianno.txt  >> $fmerge
 	 
 	 ## sort file
 	 less $fmerge |vcf-sort -c > $fb.tem
 	 egrep '^C'  $fb.tem -m 1 > $fmerge #annotated/Columbia_Freeze_Eight.hg19.avinput.hg19_multianno.txt 
 	 egrep -v '^C'  $fb.tem  >> $fmerge
fi

## done ***  prepared annotated files  **

## run psap per family 
 mkdir subfam
 cd subfam
 mkdir annotated
 
 ## separated to tri or others. run family-based psap for trio.
 
input=$(readlink -f $fmerge)

#tped="RGN.individual.ped"
#awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1){print} if(length($3)>1 && length($4)>1){next}print}' $fped > $tped
nrow=$(wc $fped|awk '{print $1}')

for n in $(eval echo {1..$nrow})
do
if  [ $(( n % 20 )) -ne 0 ] ; then

nohup bash $script_folder/script/call.sample.sh $n $input $fped &
else 
nohup bash $script_folder/script/call.sample.sh $n $input $fped  


fi

done 

wait
nohup bash $script_folder/script/PSAP_annotate.sh  $fvcf $fped  > nohup.anno.log &
