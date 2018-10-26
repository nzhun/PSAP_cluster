
r<-read.table("~/Dropbox (CGC)/CM/CM_phenotype.txt",header=1,stringsAsFactors = F,sep="\t",quote = "")

r$Father="X"
r$Mother="X"

for(fam in unique(r$Family_ID)){
  fids<-which(r$Family_ID==fam)
  pid<-intersect(fids,grep("self|proband|sibling|twin",r$Relationship_to_Proband,ignore.case = T));
  fid<-intersect(fids,intersect(grep("^father|^parent",r$Relationship_to_Proband,ignore.case = T),grep("^m",r$Gender,ignore.case = T)));
  
  mid<-intersect(fids,intersect(grep("^mother|^parent",r$Relationship_to_Proband,ignore.case = T),grep("^f",r$Gender,ignore.case = T)));
  if(length(fid)>0){
    r$Father[pid]<-r$SampleID[fid]
  }
  if(length(mid)>0){
    r$Mother[pid]<-r$SampleID[mid]
  }
  
  cid<-intersect(fids,grep("child",r$Relationship_to_Proband,ignore.case = T));
  pfid<-intersect(fids,intersect(grep("^self|^proband",r$Relationship_to_Proband,ignore.case = T),grep("^m",r$Gender,ignore.case = T)));
  
  pmid<-intersect(fids,intersect(grep("^self|^proband",r$Relationship_to_Proband,ignore.case = T),grep("^f",r$Gender,ignore.case = T)));
  
  if(length(pfid)>0){
    r$Father[cid]<-r$SampleID[pfid]
  }
  if(length(pmid)>0){
    r$Mother[cid]<-r$SampleID[pmid]
  }
  
  #pgid<-intersect(fids,grep("child",r$Relationship_to_Proband,ignore.case = T));
  pgfid<-intersect(fids,intersect(grep("^PGF",r$Relationship_to_Proband,ignore.case = T),grep("^m",r$Gender,ignore.case = T)));
  
  pgmid<-intersect(fids,intersect(grep("^PGM",r$Relationship_to_Proband,ignore.case = T),grep("^f",r$Gender,ignore.case = T)));
  
  if(length(pgfid)>0){
    r$Father[fid]<-r$SampleID[pgfid]
  }
  if(length(pgmid)>0){
    r$Mother[fid]<-r$SampleID[pgmid]
  }
  
  
  mgfid<-intersect(fids,intersect(grep("^MGF",r$Relationship_to_Proband,ignore.case = T),grep("^m",r$Gender,ignore.case = T)));
  
  mgmid<-intersect(fids,intersect(grep("^MGM",r$Relationship_to_Proband,ignore.case = T),grep("^f",r$Gender,ignore.case = T)));
  
  if(length(mgfid)>0){
    r$Father[mid]<-r$SampleID[mgfid]
  }
  if(length(mgmid)>0){
    r$Mother[mid]<-r$SampleID[mgmid]
  }
  
}

r$Gender[grep("^m",r$Gender,ignore.case = T)]<-1

r$Gender[grep("^f",r$Gender,ignore.case = T)]<-2

r$Affected[grep("^y",r$Affected,ignore.case = T) ]<- 2
r$Affected[grep("^n",r$Affected,ignore.case = T)] <- 1
#ar<-r[which(r$Affected==2),c("Family_ID","SampleID","Father","Mother","Gender","Affected","Disease","Birth_Year","Relationship_to_Proband","Race")]
#write.table(ar[,1:6],row.names = F,file = "~/Dropbox (CGC)/CM/CM_phenotype.PSAP.ped",quote = F,sep = "\t")




ar=read.table("~/Dropbox (CGC)/CM/CM_phenotype.PSAP.ped",stringsAsFactors = F,header = 1)
f3=read.table("~/Dropbox (CGC)/CM/Freeze3.CM.txt",stringsAsFactors = F,header = F)
write.table(ar[which(ar$SampleID%in%f3$V1),1:6],row.names = F,file = "~/Dropbox (CGC)/CM/CM_phenotype.F3.PSAP.ped",quote = F,sep = "\t")


f4=read.table("~/Dropbox (CGC)/CM/Freeze4.CM.txt",stringsAsFactors = F,header = F)
write.table(ar[which(ar$SampleID%in%f4$V1),1:6],row.names = F,file = "~/Dropbox (CGC)/CM/CM_phenotype.F4.PSAP.ped",quote = F,sep = "\t")



f5=read.table("~/Dropbox (CGC)/CM/Freeze5.CM.txt",stringsAsFactors = F,header = F)
write.table(ar[which(ar$SampleID%in%f5$V1),1:6],row.names = F,file = "~/Dropbox (CGC)/CM/CM_phenotype.F5.PSAP.ped",quote = F,sep = "\t")



f6=read.table("~/Dropbox (CGC)/CM/Freeze6.CM.txt",stringsAsFactors = F,header = F)
write.table(ar[which(ar$SampleID%in%f6$V1),1:6],row.names = F,file = "~/Dropbox (CGC)/CM/CM_phenotype.F6.PSAP.ped",quote = F,sep = "\t")


done<-read.table("~/Dropbox (CGC)/CM/done.list",stringsAsFactors = F,header = F)

print(c(length(which(ar$SampleID%in%f3$V1)),length(which(done$V1%in%f3$V1))))
print(c(length(which(ar$SampleID%in%f4$V1)),length(which(done$V1%in%f4$V1))))
print(c(length(which(ar$SampleID%in%f5$V1)),length(which(done$V1%in%f5$V1))))
print(c(length(which(ar$SampleID%in%f6$V1)),length(which(done$V1%in%f6$V1))))


