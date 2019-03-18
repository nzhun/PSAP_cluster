args = commandArgs(trailingOnly=T)
adir = args[3]

## IMPORTANT FAMILY INFO
#cohort.id<-unlist(strsplit(args[1],".avinput",fixed=T))[1]
outf=args[1]
cohort.id<-dirname(outf)
emp=c("X","0","-")
ped<-read.table(file=args[2],header=F,stringsAsFactors=F,fill=T,sep="\t",blank.lines.skip= T )
uf = ped$V2[which(ped$V6==1)]
n.uf = length(uf)
if(n.uf==0){
uf<-unlist(lapply(which(ped$v6==2),FUN=function(x){if(nchar(ped$V3[x]) >1){return(ped$V3)}}))
uf<-c(uf,unlist(lapply(which(ped$v6==2),FUN=function(x){if(nchar(ped$V4[x]) >1){return(ped$V4)}})))
print(uf)
#if(nchar(ped$V3[which(ped$V6==2)])>1){uf=ped$V3[which(ped$V6==2)]}
#if(nchar(ped$V4[which(ped$V6==2)])>1){uf=c(uf,ped$V4[which(ped$V6==2)])}
#uf=uf[which(!uf %in% emp)]
n.uf=length(uf)
}
af<-ped$V2[which(ped$V6==2)]
n.af<-length(af)

id<-c(af,uf)
coverage.info = read.table(file=paste(adir,"psap/lookups/gene_coverage_stats_final_12172014.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)
low.coverage = coverage.info[which(coverage.info$Mean.Coverage < 10),"Gene"]
hgmd = read.table(file=paste(adir,"psap/lookups/hgmd_pro_2013_4.12202014.annotated.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)
hgmd.ad = unique(subset(hgmd,ModeInher == "AD")$Gene.wgEncodeGencodeBasicV19)
hgmd.ar = unique(subset(hgmd,ModeInher == "AR")$Gene.wgEncodeGencodeBasicV19)

rm(list = c("coverage.info","hgmd"))

models = c("DOM-het","REC-hom","REC-chet")
genos = c("het","hom")
print("databases loaded")
## ANALYSIS STEP 1: READ IN DATA FILES FOR EACH INDIVIDUAL IN THE ANALYSIS
##  data after applying the popStat
print("analysing affected individuals")

af.dat = list()
for(i in 1:length(af)){
  print(paste(cohort.id,"/",af[i],"_popStat.txt",sep=""))
  outfinal=paste(cohort.id,"/",af[i],".report.txt",sep="")
  dat<-read.table(file=paste(cohort.id,"/",af[i],"_popStat.txt",sep=""),sep="\t",header=T,stringsAsFactors=F,check.names=F,comment.char = "")
  dat$vid = paste(dat$Chr, dat$Start, dat$Alt,sep=":")
  dat$pid = af[i]
  names(dat)[which(names(dat)==af[i])]<-"pid.geno"
  print(af[i])
  
  for(m in 1:length(models)){
    tmp = dat[which(dat$Dz.Model == models[m]),]
    
    if(models[m] == "REC-chet" & nrow(tmp) > 0){
      a1 = dat[which(dat$Dz.Model == "DOM-het" & dat$Gene.wgEncodeGencodeBasicV19 %in% tmp$Gene.wgEncodeGencodeBasicV19),]
      a1 = merge(a1[-which(names(a1) %in% c("popScore","Dz.Model"))],tmp[c("Gene.wgEncodeGencodeBasicV19","popScore","Dz.Model")])
      tmp = rbind(tmp,a1)
      tmp$chet.id=paste(tmp$vid,tmp$popScore,sep=":")
    }
      
    if(i == 1){
      af.dat[[m]] = tmp
    }else{
      cnm<-intersect(names(af.dat[[m]]),names(tmp))
      af.dat[[m]] = rbind(af.dat[[m]][,cnm],tmp[,cnm])
    }
  }
}

candidates = data.frame()
for(m in 1:length(models)){
  if(models[m] == "REC-chet"){
      if(dim(af.dat[[m]])[1]>0){
      	tmp = do.call(rbind,by(af.dat[[m]],af.dat[[m]]$chet.id, function(dat) return(data.frame(unique(dat[which(!names(dat) %in% c("pid","chet.id"))]),pid=paste(dat$pid,collapse=","),n=nrow(dat)))))
      }
  }else{
    if(dim(af.dat[[m]])[1]>0){
      tmp = do.call(rbind,by(af.dat[[m]],af.dat[[m]]$vid, function(dat) return(data.frame(unique(dat[which(names(dat) != "pid")]),pid=paste(dat$pid,collapse=","),n=nrow(dat)))))
    }
  }
 candidates = rbind(candidates,tmp)
}
## ANALYSIS STEP 2: Validate inheritance models
if(n.uf > 0){
	print("validating against unrelated individuals")
	uf.dat = data.frame()
	for(i in uf){
		dat<-read.table(file=paste(cohort.id,"/",i,"_popStat.txt",sep=""),sep="\t",header=T,stringsAsFactors=F,check.names=F,comment.char = "")
		dat$vid = paste(dat$Chr, dat$Start, dat$Alt,sep=":")
		uf.dat = rbind(uf.dat,dat)
	}

	uf.dat = unique(uf.dat[c("Gene.wgEncodeGencodeBasicV19","Chr","Start","Ref","Alt","Dz.Model","Geno","vid","popScore")])
  
  candidates$validation = "ok"
  candidates$validation[which(candidates$Dz.Model == "DOM-het" & candidates$vid %in% uf.dat$vid)] = "violation"
	candidates$validation[which(candidates$Dz.Model != "DOM-het" & candidates$vid %in% uf.dat$vid[which(uf.dat$Dz.Model != "DOM-het")])] = "violation"
	
	# missing mate
	# 1) pull out chets tagged with violation
	violations = candidates[which(candidates$Dz.Model == "REC-chet" & candidates$validation == "violation"),]
	# 2) loop over violated chet genes
	for(g in violations$gene){
	#	3) look at subset of non-violated chets in gene
		ok.chets = which(candidates$Gene.wgEncodeGencodeBasicV19 == g & candidates$Dz.Model == "REC-chet" & candidates$validation == "ok")
	# 4) extract pid for violated chets in gene
		violated.pid = unique(unlist(strsplit(as.character(violations$pid[which(violations$Gene.wgEncodeGencodeBasicV19 == g)]),",",fixed=F)))
	# 5) look for extracted pid in non-violated chets - gsub pid for space
		for(p in violated.pid){
#			candidates$validation[ok.chets][grep(p,candidates$pid[ok.chets])] = "missing mate"
			if(length(ok.chets) > 0){
				for(row in ok.chets){
					if(length(grep(",",candidates$pid[row])) > 0){
						pids = unlist(strsplit(candidates$pid[row],",",fixed=T))
						print(pids)
						pids = pids[which(pids != p)]
						print(pids)					
						if(length(pids) > 0){ 
							candidates$pid[row] = paste(pids,collapse=",") 
						}else{
							candidates$pid[row] = NA
						}
					}else{
						pids = candidates$pid[row]
						pids = pids[-which(pids == p)]
						if(length(pids) > 0){ 
							candidates$pid[row] = pids 
						}else{
							candidates$pid[row] = NA
						}
					}
				}
			}
		}
	}
}else{
  candidates$validation = "ok"
}

print("validation complete")

validated = candidates[which(is.na(candidates$pid) == F & candidates$validation != "violation"),]
print(dim(validated))
rm(candidates)

## ANALYSIS STEP 3: FLAG SITES THAT HAVE LOW COVERAGE IN THE MAC61K DATA, MISSING ALLELE FREQUENCIES IN THE MAC61K DATA BUT NOT IN THE ESP AND 1000 GENOMES DATA, OR IS AN HGMD GENE WITH A MATCHING MODE OF INHERITANCE
validated$Flag = 0
validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% low.coverage & is.na(validated$ExAC_ALL)==F)] = validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% low.coverage & is.na(validated$ExAC_ALL)==F)] + 100
validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ad & validated$Dz.Model == "DOM-het")] = validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ad & validated$Dz.Model == "DOM-het")] + 20
validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ar & validated$Dz.Model %in% c("REC-hom","REC-chet"))] = validated$Flag[which(validated$Gene.wgEncodeGencodeBasicV19 %in% hgmd.ar & validated$Dz.Model %in% c("REC-hom","REC-chet"))] + 3

## WRITE VALIDATED CANDIDATES TO FILE
print("generating report file")

print (paste("output to",outfinal))
write.table(unique(validated[order(validated$popScore,validated$Gene.wgEncodeGencodeBasicV19,validated$Dz.Model),]),file=outfinal,sep="\t",col.names=T,row.names=F,quote=F)
