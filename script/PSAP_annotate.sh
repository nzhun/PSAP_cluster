#fscript=$PSAP_PATH #/home/nz2274/Application/PSAP_cluster/ #/home/local/ARCS/nz2274/Freeze8/
#REGINDB="/home/nz2274/Resources/Region_Lib/region_db.0based.bed.gz" #/home/local/ARCS/nz2274/Resources/gnomad/2.0.2/coverage/region_db/region_db.0based.bed.gz"


echo "input: vcf.gz ped"
if [ "$#" -lt 2 ]; then
   echo "input: vcf.gz ped"
   exit
fi


vcf=$1 #"/home/local/ARCS/nz2274/Freeze8/Columbia_Freeze_Eight.hg19.vcf.gz" #"../Columbia";
ped=$2; # pedigree file




function anno_site {
  file=$1
  awk 'BEGIN{FS="\t";OFS="\t"}{
  a=$4;b=$5;
  if(length($4)>1 && length($4)<length($5)) {
  a=substr($4,0,1);
  b=substr($5,0,length($5)-length($4)+1);
  }
  if(length($5)>1 && length($5)<=length($4)) {
  a=substr($4,0,length($4)-length($5)+1);
  b=substr($5,0,1);};
  $4=a;
  $5=b;
  
  print $1,$2,$2+length($4)-1,$4,$5,$3,$_
  }'  $file|cut -f 1- > $file.format2.txt
  echo $file;
  perl $ANNOVAR/table_annovar.pl $file.format2.txt $ANNHDB --buildver hg19 --remove -protocol refGene,esp6500siv2_all,esp6500siv2_aa,esp6500siv2_ea,1000g2015aug_all,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_afr,1000g2015aug_sas,exac03,cadd13,genomicSuperDups,clinvar_20160302,mcap,revel,gnomad_exome,gnomad_genome -operation g,f,f,f,f,f,f,f,f,f,f,f,r,f,f,f,f,f -otherinfo  -nastring .
 # bash ~/Pipeline/NA_script/annobed.simple.sh  $file.format2.txt
  col=$(head -n 1 $file.format2.txt.hg19_multianno.txt|awk '{print NF}')
  
  head -n 1 $file.format2.txt.hg19_multianno.txt|cut -f 1-$((col-1)) > $file.anno.bed
  head -n 2 $file.format2.txt.hg19_multianno.txt|tail -n 1|cut -f ${col}- >> $file.anno.bed
  # 1258  less -S ../../run.flow.sh
  # 1259  sed ':a;N;$!ba;s/\n//g'  header.txt
  sed -i ':a;N;$!ba;s/\n/\t/g' $file.anno.bed;
  sed -i 's/^Chr|^chr/^#C/g'  $file.anno.bed
  
  egrep -v '^#|^C|^X.CHROM' -i $file.format2.txt.hg19_multianno.txt >> $file.anno.bed
  rm $file.format2.txt.hg19_multianno.txt
  rm $file.format2.txt
                
}

if [ "$#" -ne 2 ]
then
  echo "input format: vcf ped"
   exit
fi
for file in annotated/*report.txt
do
 #"COL-CHUNG_F14-1_ND-0234527700.report.txt";
 bn=$(basename $file|sed 's/.report.txt//g')
#if [ ! -f  "file.anno.bed" ]; then
# anno_site $file  &
#fi
# bn=$(basename $file|sed 's/.report.txt//g')
  if [ ! -e "$ANNOVAR" ]; then 
   echo "cannot find annovar, please set ANNOVAR=XXX (the annovar installed folder) and ANNHDB=XXX ( the annotation datasets for annovar) in your path"
   exit;
  fi
 if [ ! -f $file.anno.bed ]; then
 anno_site $file  &
 fi
wait 
 sz=$(wc -c $file.anno.bed|awk '{print $1}')

 if [ "$sz"  -lt 4000 ]; then
    anno_site $file &
 fi

done

wait

rm merge.bed

if [ ! -f $REGIONDB ]; then
 echo "cannot find REGIONDB in your path, please set REGIONDB=XXX ( the addtional region-based annotation file that you wantted to annotate the PSAP output, such as pLI score, Mappability...), otherwise, the annotation ends here!.\n"
 exit;
fi 


for f in annotated/*.anno.bed
do
cut -f 1-5 $f|egrep -v "^#|^Chr"|awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2-1,$3,$4,$5}' |vcf-sort -c >> merge.bed
done
less merge.bed|sort -k1,1n -k2,2n -k3,3n -u |vcf-sort -c > temp
mv temp   merge.bed
echo "merge bed generated"
tabix  $REGIONDB  -T merge.bed -h |bgzip -c > ask.db.bed.gz
tabix -f ask.db.bed.gz
echo $vcf
tabix  $vcf -T merge.bed -h |bgzip -c  > subset.vcf.gz
perl ${PSAP_PATH}/script/vcf2bed_full.pl subset.vcf.gz
#key=$(basename $f|cut -f 1 -d".")
#p=$(egrep $key $ped|cut -f 2)
#echo  $key."\t".$p."\n";
awk 'BEGIN{FS="\t";OFS="\t"}{fac=0;faf=0;fan=0;if(NR==1){for(i=7;i<10;i++){if($i=="AC"){fac=i;}if($i=="AF"){faf=i} if($i=="AN"){fan=i}} print $1,$2,$3,$4,$5,"AF\tAC\tAN",$(NF-2),$(NF-1),$NF;next};s=$1"\t"$2"\t"$3"\t"$4"\t"$5; if(faf >0){s=s"\t"$faf}else{s=s"\t"0;} if(fac>0){s=s"\t"$fac}else{s=s"\t0"}if(fan>0){s=s"\t"$fan}else{s=s"\t0"};  print s,$(NF-2),$(NF-1),$NF }'  subset.vcf.gz.bed|vcf-sort -c |bgzip -c > subset.vcf.bed.gz
tabix -f subset.vcf.bed.gz

tabix subset.vcf.bed.gz -H > ask.db.local.bed
tabix ask.db.bed.gz -H >> ask.db.local.bed
echo -e "overlap_local" >> ask.db.local.bed
sed  -i ':a;N;$!ba;s/\n/\t/g' ask.db.local.bed
bedtools intersect -a subset.vcf.bed.gz -b ask.db.bed.gz -wao  >> ask.db.local.bed
bgzip -f ask.db.local.bed
tabix -f ask.db.local.bed.gz

echo "sub query db generated"
tabix ask.db.local.bed.gz -H > db.header.txt

for f in annotated/*.anno.bed
do
head -n 1 $f > $f.addregion.bed;
cat db.header.txt >> $f.addregion.bed;
echo "overlap" >> $f.addregion.bed;
sed -i ':a;N;$!ba;s/\n/\t/g' $f.addregion.bed;
#set -i 's/^C/#C/g' $f.addregion.bed ;
echo "$f tested";
sed -i 's/^C/#C/g' $f;
sed -i 's/^c/#C/g' $f;
key=$(basename $f|cut -f 1 -d".");
p=$(egrep $key $ped|cut -f 2|egrep "$key"|uniq);
less ask.db.local.bed.gz|egrep "^#|^C|^c|^X.c|^X.C|$p" |bgzip -c > db.$p.local.bed.gz;
tabix -f db.$p.local.bed.gz;
s=$(awk -v id=$p '{for(i=1;i<NF+1;i++){if($i=="Dz.Model"||$i=="Dz.Model."id){print i;last}}}' $f.addregion.bed);
echo $s;
bedtools intersect -a $f -b db.$p.local.bed.gz -wao |sort -k1,1n -k2,2n -k3,3n -k4,4d -k5,5d -k$s,$s -t$'\t' -u|vcf-sort -c  >> $f.addregion.bed ;
rm db.$p.local.bed.gz*;
done 



wait

rm ask.db.local.bed.gz*
rm db.header.txt
rm subset.vcf*

freq=0.01
folder="OUTPUT_$freq"
mkdir $folder

for file in $(ls annotated/*.addregion.bed); do 
key=$(basename $file|cut -f 1 -d".")
p=$(cut -f 2 $ped|egrep "$key"|uniq)
perl ${PSAP_PATH}/script/extract_keys.pl $file ${PSAP_PATH}/src/format.header.txt $folder  $p $ped $freq
ftemp=$(ls $folder/*$key*.extract.txt)
 sort -k3,3g $ftemp > temp
 mv temp $ftemp

cp ${PSAP_PATH}/src/README.txt $folder/
done


freq=0.1
folder="OUTPUT_$freq"
mkdir $folder

for file in $(ls annotated/*.addregion.bed); do 
key=$(basename $file|cut -f 1 -d".")
p=$(egrep $key $ped|cut -f 2)
perl ${PSAP_PATH}/script/extract_keys.pl $file ${PSAP_PATH}/src/format.header.txt $folder  $p $ped $freq
ftemp=$(ls $folder/*$key*.extract.txt)
sort -k3,3g $ftemp > temp
mv temp $ftemp

cp ${PSAP_PATH}/src/README.txt $folder/
done



#mkdir OUTPUT_0.01

#for file in $(ls annotated/*.addregion.bed); do perl /home/local/ARCS/nz2274/Pipeline/NA_script/extract_keys.pl $file /home/local/ARCS/nz2274/Resources/PSAP/format.header.txt OUTPUT_0.01 0.01 & done

wait


echo "done"
