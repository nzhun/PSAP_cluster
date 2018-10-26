	#PSAP_PATH=/home/nz2274/Application/PSAP_cluster/ #/Users/nazhu/Documents/GitHub/PSAP_cluster/ #INSERT PATH TO PSAP DIRECTORY HERE eg. /scratch/dclab/
#	ANNOVAR_PATH=$ANNOVAR #/home/local/ARCS/nz2274/Application/annovar/;
	echo "input nrow input.hg19_multianno.txt input.ped"
	tped=$3; #"RGN.individual.ped"
	input=$2; #"../annotated/Columbia_Freeze_Eight.hg19.avinput.hg19_multianno.txt"
	 N=$1;
   	fam=0;
	bid=$(basename $tped|sed 's/\.ped//g')
	sped="annotated/$N.$fam.$bid.ped"
	head -n 1 $tped|cut -f 1-6  > $sped
	head -n $N $tped|tail -1 |cut -f 1-6 >> $sped
	fid=$(cut -f 1 $sped|tail -n1)
	proband=$(egrep -v '^#' $sped|cut -f 2|tail -n1)
	father=$(egrep -v '^#' $sped|cut -f 3|egrep -v '^0$|^X$|^-$' -i|tail -n1)
	rs=$(head -n 1 $input|egrep "$father")
	if [ -z "$rs" ]; then
 		father=""
		echo "$father does not include in the vcf"
	fi
	mother=$(egrep -v '^#' $sped|cut -f 4|egrep -v '^0$|^X$|^-$' -i|tail -n 1)
	rs=$(head -n 1 $input|egrep "$mother")
	if [ -z "$rs" ]; then
 		mother=""
		echo "$mother does not include in the vcf"
	fi

	if [ ! -z "$father"  ];
	then 
	echo -e "$fid\t$father\tX\tX\t1\t1" >> $sped
	echo $father
	fam=$(( fam + 1))
	fi
	if [ ! -z "$mother" ]; then
	echo -e "$fid\t$mother\tX\tX\t2\t1" >> $sped
	fam=$(( fam + 1))
	echo $mother
	fi
	mv $sped annotated/$N.$fam.$fid.ped
	sped="annotated/$N.$fam.$fid.ped"	
	echo $proband"\t"$father"\t"$mother
	#echo "head -n 1 $input| awk -v p=$proband -v f=$father -v m=$mother  'BEGIN{FS=\"\t\";OFS=\"\t\"}{pc=0;fc=0;mc=0;for(i=26;i<NF+1;i++){if(\$i==p){pc=i} if(\$i==f){fc=i} if(\$i==m){mc=i} } o=pc; if(fc!=0){o=o\",\"fc}if(mc!=0){o=o\",\"mc} print o}'"
	fcut=$(head -n 1 $input| awk -v p=$proband -v f=$father -v m=$mother  'BEGIN{FS="\t";OFS="\t"}{pc=0;fc=0;mc=0;for(i=26;i<NF+1;i++){if($i==p){pc=i} if(length(f)>1 && $i==f){fc=i} if(length(m)>1 && $i==m){mc=i} } o=pc; if(fc!=0){o=o","fc}if(mc!=0){o=o","mc} print o}')
	echo "$fcut found"
	head -n 1 $input|cut -f 17-25,$fcut > annotated/$proband.avinput.header 
	pcol=$(awk -v p=$proband   'BEGIN{FS="\t";OFS="\t"}{pc=0;fc=0;mc=0;for(i=10;i<NF+1;i++){if($i==p){pc=i-9} } print pc}' annotated/$proband.avinput.header)
	awk -v p=$pcol 'BEGIN{FS="\t";OFS="\t"}{if(p>1){t=$(9+p);$(9+p)=$10;$10=t;}print}' annotated/$proband.avinput.header > annotated/$proband.avinput.header.tmp
	mv annotated/$proband.avinput.header.tmp annotated/$proband.avinput.header
	sed -i 's/\t$//g' annotated/$proband.avinput.header
	echo "proband is the "$pcol" column"
	echo "header prepared"

	l1=$(head -n 1 $input|cut -f 1-16)

	echo -e ${l1}"\tOtherinfo"|sed 's/ /\t/g' > annotated/$proband.avinput.hg19_multianno.txt
	#cut -f 1-25,$fcut $input|egrep -v '^C'|awk 'BEGIN{FS="\t";OFS="\t"}{if($26 ~/^0\/0|^\./){next} alc=split($21,b,",");  n=split($26,a,":");split(a[1],al,"/");if(al[1]>0 && al[2] >0 &&al[1]!=al[2]){next;} if(al[1]>0){$21=b[al[1]];  }else{$21=b[al[2]]}f($21=="*"){next} split(a[2],ad,","); if(a[3]<10||a[3]=="."||a[4]<30){next} adp=ad[al[2]+1]; rdp=ad[1]; if(al[1]!=al[2]){af=ad[al[2]+1]/a[3]; if(af <0.2){next}}else{if(ad[al[2]+1]/a[3]!=1){next}}  if(alc>1){ for(col=26;col<NF+1;col++){  n=split($col,a,":");split(a[1],al,"/"); split(a[2],ad,",");adp=ad[al[2]+1]; rdp=ad[1]; if(al[1]>0){al[1]=1;} if(al[2]>0){al[2]=1;} $col=al[1]"/"al[2];$col=$col":"rdp","adp; for(i=3;i<n+1;i++){$col=$col":"a[i]}}} print  }'  >> annotated/$proband.avinput.hg19_multianno.txt

	echo "cut -f 1-25,$fcut $input   $pcol"

	 cut -f 1-25,$fcut $input|awk -v p=$pcol 'BEGIN{FS="\t";OFS="\t"}{if(p>1){t=$(25+p);$(25+p)=$26;$26=t;}print}'|egrep -v '^C'|awk 'BEGIN{FS="\t";OFS="\t"}{if($26 ~/^0\/0|^\./){next} alc=split($21,b,",");  
	 n=split($26,a,":");split(a[1],al,"/");if(al[1]>0 && al[2] >0 &&al[1]!=al[2]){next;}
	 if(alc>1){
		if(al[1]>0){
	                        $21=b[al[1]];
	                        nalt=al[1]
	                }else{
	                        $21=b[al[2]]
	                        nalt=al[2]
	                }
	                if($21=="*"){next}

	                split($24,infos,";AF=");
	                split(infos[2],fd,";");
	                split(fd[1],afs,",");
	                if(afs[nalt] >0.05){next}
	}
	  if(al[1]>0){$21=b[al[1]];  }else{$21=b[al[2]]} if($21=="*"){next}


	    split(a[2],ad,","); if(a[3]<10||a[3]=="."||a[4]<30){next} adp=ad[al[2]+1]; rdp=ad[1]; if(al[1]!=al[2]){af=ad[al[2]+1]/a[3]; if(af <0.2){next}}else{if(ad[al[2]+1]/a[3]!=1){next}}  
	 if(alc>1){ 
	  for(k=26;k<NF+1;k++){
	  n=split($k,a,":");
	  split(a[1],al,"/");
	  split(a[2],ad,",");
	  adp=ad[al[2]+1]; 
	  rdp=ad[1];
	 if(al[1]>0){al[1]=1;} 
	 if(al[2]>0){al[2]=1;} 
	 $k=al[1]"/"al[2];
	 $k=$k":"rdp","adp; 
	 for(i=3;i<n+1;i++){$k=$k":"a[i]}} 
	 }
 
	 print  }'   >> annotated/$proband.avinput.hg19_multianno.txt

#cut -f 1-25,$fcut $input|egrep -v '^C'|awk 'BEGIN{FS="\t";OFS="\t"}{if($26 ~/^0\/0|^\./){next} n=split($26,a,":");split(a[1],al,"/");split(a[2],ad,","); if(a[3]<10||a[3]=="."||a[4]<30){next} if(al[1]!=al[2]){af=ad[al[2]+1]/a[3]; if(af <0.2){next}}else{if(ad[al[2]+1]/a[3]!=1){next}} print  }'  >> annotated/$proband.avinput.hg19_multianno.txt
#if [] ; then
#nohup bash ${PSAP_PATH}psap/annotate_PSAP.sh $PWD/annotated/$proband $(readlink -f annotated/$N.ped) > nohup.$proband.psap.lo
#mv $sped annotated/$proband.ped
PED_FILE=$(readlink -f $sped)
OUTFILE=$(readlink -f annotated/$proband)

if [ $fam -lt 2 ]; then



	echo "PROGRESS: Extracting individual IDs"
	IDS=($(awk 'BEGIN{FS="\t"}{if($1 ~/^#/){next}print $2}' $PED_FILE))
	IDX=1

	# RUN apply_popStat_individual.R for each individual
	echo "PROGRESS: Starting PSAP annotation"
	for i in ${IDS[@]}
	do
		echo "Rscript ${PSAP_PATH}psap/RScripts/apply_popStat_individual.R ${OUTFILE}.avinput $i $PSAP_PATH $PED_FILE"
		Rscript ${PSAP_PATH}psap/RScripts/apply_popStat_individual.R ${OUTFILE}.avinput $i $PSAP_PATH $PED_FILE &
	   if [ `expr $IDX % 10` -eq 0 ]
	   then
	        echo "PROGRESS: Annotating individuals" $(($IDX - 10)) "-" $IDX "out of" ${#IDS[@]}
	        wait # Limit number of individuals annotated to no more than 20 at a time
	   fi
	   IDX=$(($IDX+1))
	done
	wait

	# Generate report file - will look for variants present in two or more affected with PSAP < 1e-3 and not observed in unaffected
	echo "PROGRESS: Generating report file for all individuals"
	echo " Rscript ${PSAP_PATH}psap/RScripts/unrelated_candidate_analysis.R ${OUTFILE}.avinput $PED_FILE $PSAP_PATH"
	Rscript ${PSAP_PATH}psap/RScripts/unrelated_candidate_analysis.R ${OUTFILE}.avinput $PED_FILE $PSAP_PATH
else
	
	echo "start psap"
	echo "bash ${PSAP_PATH}psap/annotate_PSAP.sh ${OUTFILE}.avinput $(readlink -f $PED_FILE)"

	#if [] ; then
	nohup bash ${PSAP_PATH}psap/annotate_PSAP.sh ${OUTFILE}.avinput $(readlink -f $PED_FILE) > nohup.$proband.psap.log
		
fi
