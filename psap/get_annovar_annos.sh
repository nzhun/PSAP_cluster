# Initialization script that downloads all necessary annotation tables prior to running the PSAP pipeline for the first time
ANNOVAR_PATH=$ANNOVAR #INSERT PATH TO ANNOVAR DIRECTORY HERE eg. /scratch/dclab/annovar/
PSAP_PATH=$PSAP #INSERT PATH TO PSAP DIRECTORY HERE eg. /scratch/dclab/

#cd $ANNOVAR_PATH
# Download Sep 2014 1000 Genomes allele frequencies from ANNOVAR
perl $ANNOVAR_PATH/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar 1000g2014sep $ANNOVAR/humandb/

# Download CADD V1.3 scores from University of Washington and formats the data into annovar format
perl $ANNOVAR_PATH/annotate_variation.pl -v -downdb cadd -buildver hg19 -webfrom annovar $ANNOVAR/humandb/

# Download ESP 6500 allele frequencies from ANNOVAR
perl $ANNOVAR_PATH/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar esp6500si_all $ANNOVAR/humandb/

# Download dbSNP137 annotations from ANNOVAR
perl $ANNOVAR_PAT/Hannotate_variation.pl -downdb -buildver hg19 -webfrom annovar snp137 $ANNOVAR/humandb/

# Download GencodeV19 gene map from UCSC
perl $ANNOVAR_PATH/annotate_variation.pl -downdb -buildver hg19 wgEncodeGencodeBasicV19 $ANNOVAR/humandb/
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y; do wget -P $ANNOVAR/humandb/ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${chr}.fa.gz; done
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y; do gunzip $ANNOVAR/humandb/chr${chr}.fa.gz; done
perl $ANNOVAR_PATH/retrieve_seq_from_fasta.pl -format genericGene -seqdir $ANNOVAR/humandb/ $ANNOVAR/humandb/hg19_wgEncodeGencodeBasicV19.txt --outfile hg19_wgEncodeGencodeBasicV19Mrna.fa

# Move ExAC allele frequencies provided with the PSAP pipeline to the ANNOVAR annotation folder and unzip
mv ${PSAP_PATH}psap/lookups/hg19_mac63kFreq_ALL.txt.gz $ANNOVAR/humandb/
gzip -d $ANNOVAR/humandb/hg19_mac63kFreq_ALL.txt.gz
