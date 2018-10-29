echo "input: vcf.gz pedigree"
if [ "$#" -lt 2 ];then
   echo "input:vcf.gz pedigree"
   exit
fi
fvcf=$1 #"/home/local/ARCS/nz2274/Freeze8/Columbia_Freeze_Eight.NF.pVCF.vcf.gz"
fped=$2 #"/home/local/ARCS/nz2274/Freeze8/RGN8.sampleHas.ped"
script_folder=$PSAP_PATH
fmerge=$3
#fname=$(basename $fvcf)
#fb=${fname%.*}
#cb=$(basename $chain|sed 's/.chain.gz//g')
#mkdir annotated
#fmerge="annotated/$fb.hg19.avinput.hg19_multianno.txt"


input=$(readlink -f $fmerge )

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
