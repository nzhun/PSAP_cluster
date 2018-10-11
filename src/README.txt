## the method of PSAP: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5127779/
## Before applying PSAP, the following filter were applied:
## only variants in coding region in Genecode basic v19 are kept
## the synonymous variants are removed 
## ExACfreq < 0.01
## GQ <30 are removed
## DP < 10 are removed
## allelic balance (alt depth/total depth) <0.2 are removed for heterozygote
## cohort AF > 0.1 (depends on the number of cohorts) are removed
## gnomad Exome depth median <10 are removed 

## After applying PSAP, the following filter were applied:
## popscore >0.01 are removed

## columns:
## PROBAND, Proband ID in vcf
## GENE.WGENCODEGENCODEBASICV19, Gene symbol,Gene annotation used GENCODEBASICV19 https://www.gencodegenes.org/releases/19.html 
## POPSCORE, popscore from PSAP analysis.https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5127779/
## FUNC.WGENCODEGENCODEBASICV19, Function of the variants
## EXONICFUNC.WGENCODEGENCODEBASICV19, Exonic function of the variants
## AACHANGE.WGENCODEGENCODEBASICV19, the resulted amino change 
## DZ.MODEL, popscore based disease model. Three different models used for estimating popscore, DOM-Het: dominicant heterozygote, Rec-chet: recessive compound heterozygote, Rec-hom: recessive homozygote 
## FORMAT, AD:Allelic depths for the ref and alt alleles in the order listed, DP: Approximate read depth,GQ: Genotype Quality, GT:Genotype, MIN_DP: Minimum DP observed within the GVCF block, PGT: Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another, PID: Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group,PL: Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification,RGQ; Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong, SB;Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias
## GENOTYPE, the genotype of proband 
## #CHROM, chromosome
## POS: coordinates
## ID: dbsnp147 
## REF: reference allele
## ALT: alternated allele
## AC: altnated allele count in the cohort
## AF: alternated allele frequency in the cohort
## AN: allele count in the cohort
## EXAC_ALL: ExAC frequency,http://exac.broadinstitute.org/about  
## GNOMAD_EXOME_ALL: gnomad Exome frequency, http://gnomad.broadinstitute.org/about
## GNOMAD_GENOME_ALL: gnomad genome freuency
## 1000G2015AUG_ALL: 1000 genome frequency 
## CADD_PHRED, CADD score 1.0
## CADD13_PHRED, CADD score 1.3
## CLINSIG,clinical significant from clinvar
## MCAP, MCAP predicted missense deleterious score, the larger the more deleterious,http://bejerano.stanford.edu/mcap/
## REVEL: REVEL predicated missense deleterious score,similar performance as MCAP, the larger the more deleterious,https://www.ncbi.nlm.nih.gov/pubmed/27666373
## MPC: MPC predicated missense deleteriousness, higher sensitivity https://www.biorxiv.org/content/early/2017/06/12/148353, performance comparison among predictors,https://www.biorxiv.org/content/biorxiv/early/2018/02/02/259390.full.pdf
## GENOMICSUPERDUPS, Duplications of >1000 Bases of Non-RepeatMasked Sequence (>90 percent similar),http://varianttools.sourceforge.net/Annotation/GenomicSuperDups.
## GNOMAD_EXOME_MEAN, median sequence depth of coverg in gnomad Exome
## MAPPABILITY, it is prior estimated score to indicate how likely a position can be uniquely mapped to the genome on account of the sequenced reads length. ttp://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377
## MIS_Z, Z score of missense variants of gene, zscore is the deviation of ExAC dataset observation from mutation rate based expection .
## LOF_Z, Z score of loss of function variants in a gene 
## PLI, a posterior probability of a gene intolerant of a single loss of function variant,e.g. genes with pli>0.9 are roughly taken as  haploinsufficent genes. The estimation is based on the observed and expected number of loss-of-function variants per gene in ExAC dataset.
## PREC,the probability of a gene being intolerant of two loss-of-function variants, but not heterozygous loss-of-function variants
## VALIDATION, it indicates whether the genotype is valid in consider of the relatives' genotype, for example, whether compound heterozygotes have different inheritance, or whether the variant followed Mendelian law.
