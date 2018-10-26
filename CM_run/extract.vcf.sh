
#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N VCF_Extractor
#$ -l h_rt=12:00:00
#$ -l h_vmem=6G
#$ -cwd

echo extract.vcf.sh  vcf_file.gz samplelist

nrow=$SGE_TASK_ID

if [ -z "$nrow" ]; then
	echo "please use qsub -t to set the task ID"
	exit
fi
InpFil==$1
echo  $InpFil
fvcf=$(head -n $nrow $InpFil|tail -n 1|awk '{print $1}')
flist=$(head -n $nrow $InpFil|tail -n 1|awk '{print $2}')


#flist=$2

#fvcf=$1

fp=$(basename $fvcf|sed 's/.vcf.gz//g')
echo $fp
bcftools view -Oz -S $flist $fvcf -c 1  -I --force-samples > $fp.CM.vcf.gz
echo $fp.CM.vcf.gz

f="$fp.CM.vcf.gz"
OUTFILE="$fp.CM"
echo  $f
ANNOVAR_PATH=$ANNOVAR
perl ${ANNOVAR_PATH}convert2annovar.pl -format vcf4old $f -outfile ${OUTFILE}.avinput -includeinfo 2> $fp.CM.convert.log

perl ${ANNOVAR_PATH}table_annovar.pl ${OUTFILE}.avinput -remove -outfile annotated/${OUTFILE}.avinput ${ANNHDB} -buildver hg19 -protocol wgEncodeGencodeBasicV19,mac63kFreq_ALL,esp6500si_all,1000g2014sep_all,snp137,cadd -operation g,f,f,f,f,f -nastring NA -otherinfo -argument -separate,,,,,-otherinfo  > nohup.$fp.anno.log

fmerge="$fp.CM.avinput.tmp.hg19_multianno.txt"
egrep '^chr|^C' -m 1  $annotated/$fp.CM.avinput.hg19_multianno.txt  > $fmerge #annotated/Columbia_Freeze_Eight.avinput.hg19_multianno.txt
tabix $f -H |egrep '^#C' >>  $fmerge  #annotated/Columbia_Freeze_Eight.avinput.hg19_multianno.txt
#sed ':a;N;$!ba;s/\n//g'  annotated/Columbia_Freeze_Eight.avinput.hg19_multianno.txt
sed ':a;N;$!ba;s/\n//g'  $fmerge |sed 's/Otherinfo//g' > $fp.CM.temp.anno
mv $fp.CM.temp.anno  $fmerge

awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==1){next} if($6 ~/^exonic|^spli/ && $9 !~/^synony/ && ($11 <0.01||$11=="."||$11=="NA")){print}}' annotated/$fp.CM.avinput.hg19_multianno.txt >> $fmerge

#mkdir annotated 
mv $fmerge annotated/$fp.CM.filtered.hg19_multianno.txt

echo "done"
