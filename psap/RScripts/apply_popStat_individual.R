#!/bin/Rscript

print(getwd())
## NAME OF FILE FOR ANALYSIS, PROVIDED AS AN ARGUMENT WHEN CALLING THE RSCRIPT
args<-commandArgs(trailingOnly=T) ## args[1] = family - must be annovar annotated (avinput.hg19_multianno.txt) and have a separate header file for the vcf colums(header); args[2] = individual ID
adir = args[3]

prefix<-args[1]
fam.id<-dirname(args[1]) #strsplit(args[1],".avinput",fixed=T)
ped = read.table(args[4],stringsAsFactors=F,skip = 1,header = F,sep="\t")
## Individual ID - ASSUMES only one individual is being analyzed/annotated

indv.id = args[2]

## at some point this may be changed to an argument based system but for now it's hard coded
score = "CADD_Phred"
scale = seq(0,70,0.05)
print(paste(adir,"psap/lookups/lookup_genes.txt",sep=""))
lookup.genes = scan(file=paste(adir,"psap/lookups/lookup_genes.txt",sep=""),"character")


#lookup.lof = scan(file=paste(adir,"psap/lookups/full.lof.pCADD.gencodeV19.allsites.txt.gz",sep=""),"character")
lookup.lof = read.table(file=paste(adir,"psap/lookups/full.lof.pCADD.gencodeV19.allsites.txt.gz",sep=""),stringsAsFactors=F)
## READ IN AND FORMAT DATA
print(paste(prefix,".hg19_multianno.txt",sep=""))
exome.raw<-read.table(file=paste(prefix,".hg19_multianno.txt",sep=""),sep="\t",comment.char="",quote="", stringsAsFactors=F,skip=1)
print ("data")
header<-read.table(file=paste(prefix,".hg19_multianno.txt",sep=""),sep="\t",stringsAsFactors=F,nrow=1)
vcf.header<-read.table(file=paste(prefix,".header",sep=""),sep="\t",stringsAsFactors=F,comment.char="@",strip.white = T)
stopifnot(indv.id %in% vcf.header) # CHECKS THAT SPECIFIED INDIVIDUAL IS IN THE DATA
n.annos=ncol(header)
vcf.header=vcf.header[which(vcf.header!="NA")]
names(exome.raw)=c(header[-n.annos],vcf.header)

if(any(grepl("cadd",names(exome.raw))) == T){
    exome.raw[,score] = as.numeric(sapply(exome.raw$cadd,function(x) unlist(strsplit(x,","))[2]))
}

stopifnot(any(grepl("CADD_Phred",names(exome.raw))))

# Extracts genotype info for specified individual
a1 = substr(exome.raw[,indv.id],1,1)
a2 = substr(exome.raw[,indv.id],3,3)
exome.raw[indv.id] = "NA"
if(length(which(a1 != a2)) > 0){
  exome.raw[which(a1 != a2),indv.id] = "het"
}
if(length(which(!a1 %in% c(0,".") & !a2 %in% c(0,".") & a1 == a2)) > 0){
  exome.raw[which(!a1 %in% c(0,".") & !a2 %in% c(0,".") & a1 == a2),indv.id] = "hom"
}
if(length(which(a1 == 0 & a2 == 0)) > 0){
  exome.raw[which(a1 == 0 & a2 == 0),indv.id] = "ref"
}
print( indv.id)
print(ped$V2)
print(ped$V5[which(ped$V2 == indv.id)])
if(ped$V5[which(ped$V2 == indv.id)] == 1){
  exome.raw[which(exome.raw[,indv.id] %in% c("het","hom") & exome.raw[,"Chr"] == "X"),indv.id] = "hom"
  exome.raw[which(exome.raw[,indv.id] %in% c("het","hom") & exome.raw[,"Chr"] == "Y"),indv.id] = "hom"
}

if(ped$V5[which(ped$V2 == indv.id)] == 2){
  exome.raw[which(exome.raw[,"Chr"] == "Y"),indv.id] = "ref"
}
print(dim(exome.raw))


maf=names(exome.raw)[grep("ExAC_ALL|mac63kFreq_ALL",names(exome.raw),ignore.case=T)[1]]

### CLEAN DATA: 1) REMOVE BLACKLIST GENES, 2) REMOVE VARIANTS WITH AF DISCREPANCIES, 3) REMOVE GENES NOT INCLUDED IN LOOKUP TABLES, 4) REMOVE REGIONS THAT ARE NOT COVERED IN ExAC, 5) REMOVE VARIANTS THAT DO NOT PASS TRANCHE FILTER 
# 1) REMOVE BLACKLISTED GENES-- I think it would be better to remove and output these genes to a separate file (like the missing data).  I also think this should include all low coverage genes because it's more generic.
bl<-scan(paste(adir,"psap/lookups/blacklist_122814.txt",sep=""),what="character")
bl.remove = unique(c(which(exome.raw$Gene.wgEncodeGencodeBasicV19 %in% bl),grep("^HLA", exome.raw$Gene.wgEncodeGencodeBasicV19), grep("^MUC", exome.raw$Gene.wgEncodeGencodeBasicV19), grep("^KRT", exome.raw$Gene.wgEncodeGencodeBasicV19),grep("^OR", exome.raw$Gene.wgEncodeGencodeBasicV19), grep("^TRBV", exome.raw$Gene.wgEncodeGencodeBasicV19)))

# 2) REMOVE AF DISCREPANCIES (ANYTHING THAT IS MISSING IN ExAC BUT PRESENT IN 1000GP OR ESP AT GREATER THAN 5% FREQUENCY
af.remove = which(is.na(exome.raw[,maf]) == T & exome.raw[,"1000g2014sep_all"] > 0.05 | is.na(exome.raw[,maf]) == T & exome.raw$esp6500si_all > 0.05)

# 3) REMOVE GENES NOT IN LOOKUP TABLES
multi_transcript_indel<-grep(";",exome.raw$Gene.wgEncodeGencodeBasicV19) #)intersect(which(exome.raw$Ref=="-"|exome.raw$Alt=="-"),
exome.raw$Gene.wgEncodeGencodeBasicV19[multi_transcript_indel]<-unlist(lapply(multi_transcript_indel,
      FUN <-function(i){
          annoGenes<-unlist(strsplit(exome.raw$Gene.wgEncodeGencodeBasicV19[i],split = ";"));
          ids<-which(annoGenes%in%lookup.genes); 
          if(length(ids)>1){
             lof_a<-lookup.lof[which(lookup.lof$V1%in%annoGenes[ids]),1:2];
             return(lof_a[order(lof_a$V2,decreasing = T)[1],1])
          } 
          if(length(ids)<1 || is.na(ids)){ids=1}; 
          return(annoGenes[ids])
      }
))
lookup.remove = intersect(which(! exome.raw$Gene.wgEncodeGencodeBasicV19 %in% lookup.genes),multi_transcript_indel)

# 4) REMOVE LINES WHERE ALL AFs ARE MISSING
af = ped$V2[which(ped$V6 == 2)]
af_index<-grep(af,names(exome.raw))
missing.remove = which(apply(exome.raw[af_index],1,function(row) return(sum(is.na(row)))) == length(af))
rms<-c(bl.remove,af.remove,lookup.remove,missing.remove)
if(length(rms) > 0){
  tmp.exome = exome.raw[-unique(rms),]
}else{
tmp.exome=exome.raw
}
# 4a) REMOVE REGIONS NOT COVERED IN ExAC - JUST RETAINING PROTEIN CODING SITES AND SPLICE SITES
keep<-unique(c(grep("splic",tmp.exome$Func.wgEncodeGencodeBasicV19),which(is.na(tmp.exome$ExonicFunc.wgEncodeGencodeBasicV19)==FALSE)))
#print(keep)
exome<-tmp.exome[keep,]

# 4b) SCORE INDELS

indels = grep("^frameshift",exome$ExonicFunc.wgEncodeGencodeBasicV19)
gene.index = as.integer(factor(exome$Gene.wgEncodeGencodeBasicV19[indels],levels=lookup.lof[,1]))
exome[,score][indels] = lookup.lof[gene.index,2]

# 5) REMOVE VARIANTS THAT DO NOT PASS QUALITY FILTER OR HAVE MISSING pCADD SCORES
info<-exome[which(exome$FILTER=="PASS" & is.na(exome[,score]) == F | exome$FILTER=="." & is.na(exome[,score]) == F),
            c(unlist(vcf.header[1:5]),"Chr","Start","Ref","Gene.wgEncodeGencodeBasicV19","Func.wgEncodeGencodeBasicV19",
              "ExonicFunc.wgEncodeGencodeBasicV19","AAChange.wgEncodeGencodeBasicV19",
             maf,"1000g2014sep_all","esp6500si_all","Alt",score,indv.id)]

## OUTPUT MISSING DATA/DATA NOT INCLUDED IN ANY OF THE ABOVE ANALYSES
id.raw = paste(exome.raw$Chr,exome.raw$Start,exome.raw$Ref,exome.raw$Alt,sep=":")
id.final = paste(info$Chr,as.numeric(info$Start),info$Ref,info$Alt,sep=":")
missing<-unique(exome.raw[which(! id.raw %in% id.final),
  c(unlist(vcf.header[1:5]),"Chr","Start","Ref","Gene.wgEncodeGencodeBasicV19",
    "Func.wgEncodeGencodeBasicV19","ExonicFunc.wgEncodeGencodeBasicV19","AAChange.wgEncodeGencodeBasicV19",
    maf,"1000g2014sep_all","esp6500si_all","Alt",score,indv.id)])

rm(list=c("keep","exome","tmp.exome","exome.raw","af.remove","lookup.remove","bl.remove","bl","lookup.genes"))
class(info[,score]) = "numeric"

print("data cleaned, beginning annotation")
tmp<-info
rms<-which(names(info) == indv.id)
if(length(rms)<1){
  tmp = info[-rms]
}
tmp["Geno"] = info[,which(names(info) == indv.id)]
print(dim(tmp))
out = data.frame()

## SOURCES CODE THAT WILL FORMAT AND ANALYZE THE DATA FOR EACH MODE OF INHERITANCE (AD, AR, CHET or X) MODEL
## AD MODEL CODE
dir<-adir
source(paste(dir,"psap/RScripts/apply_pop_stat_het.R",sep=""))
print("AD model complete")
  
## AR MODEL CODE
source(paste(adir,"psap/RScripts/apply_pop_stat_hom.R",sep=""))
print("AR-hom model complete")
  
## CHET MODEL CODE
source(paste(adir,"psap/RScripts/apply_pop_stat_chet_unphased.R",sep=""))
print("AR-het model complete")
print(dim(out))  
print("processing data")
out = out[which(!names(out) %in% c("i","j"))]
out[which(out$popScore == 0),"popScore"] = 1e-6

## WRITE OUTPUT FOR FAMILY
print("writing file")
write.out = out[which(is.na(out$popScore) == F),]
print(dim(write.out))
outfile=paste(fam.id,"/",indv.id,sep="")
print( paste(paste(outfile,"_popStat.txt",sep="")))
write.table(write.out,file=paste(outfile,"_popStat.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

missing.out = missing[-which(missing$Func.wgEncodeGencodeBasicV19 == "exonic" | missing$Func.wgEncodeGencodeBasicV19 == "splicing;intronic" | missing$Func.wgEncodeGencodeBasicV19 =="splicing;exonic" | missing$Func.wgEncodeGencodeBasicV19 =="exonic;splicing"),]
write.table(missing.out,file=paste(outfile,"_missing_data.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)

