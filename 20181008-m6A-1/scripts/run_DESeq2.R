## Parser for running DESeq2
## ZZ
## 7.25.2018

args <- commandArgs(TRUE)
samplelistA=as.character(args[1])     #e.g., Xist_knockdown//sleuth/MaleNoCre-FemaleNoCre/b1.txt
samplelistB=as.character(args[2])     #e.g., Xist_knockdown//sleuth/MaleNoCre-FemaleNoCre/b2.txt
output=as.character(args[3])          #e.g., Xist_knockdown//sleuth/MaleNoCre-FemaleNoCre/sleuth.txt
genome_type=as.character(args[4])   #e.g., mm10
input.IDconversion=as.character(args[5])    #e.g., Xist_knockdown//config/t2g_mm10.txt
#input.metadata=as.character(args[5])        #e.g., Xist_knockdown//config/metadata.txt
#input.IDconversion=as.character(args[6])    #e.g., Xist_knockdown//config/t2g_mm10.txt

print(samplelistA)
print(samplelistB)
print(output)
print(genome_type)
#print(input.metadata)
print(input.IDconversion)

rawlistA=strsplit(as.matrix(read.table(samplelistA)),split=",")[[1]]
rawlistB=strsplit(as.matrix(read.table(samplelistB)),split=",")[[1]]
getID=function(string){
  temp=strsplit(string,split="/")[[1]]
  return(temp[length(temp)-1])
}
listA=as.character(sapply(rawlistA,getID))
listB=as.character(sapply(rawlistB,getID))

############
#run DESeq#
############
library("DESeq2")

base_dir = paste(strsplit(rawlistA[1],split="/")[[1]][1:(length(strsplit(rawlistA[1],split="/")[[1]])-2)],collapse="/")

#Next get the list of sample IDs with
sample_id <- c(listA,listB)

#A list of paths to the kallisto results indexed by the sample IDs is collated with
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))

#The next step is to load an auxillary table that describes the experimental design and the relationship between the kallisto directories and the samples:

s2c = data.frame(
	sample=c(listA, listB),
	condition=c(rep(1, length(listA)), rep(0, length(listB))),
	path=kal_dirs,
	stringsAsFactors=F
	)
	
## COMMENT: DO NOT USE "dplyr" AS IT DOES NOT DO ANYTHING
## ZIJUN, JUNE 12 2018
#Now, we must enter the directories into a column in the table describing the experiment. 
#This column must be labeled path, otherwise sleuth will throw an error. 
#This is to ensure that the user can check which samples correspond to which kallisto runs.
#s2c <- dplyr::mutate(s2c, path = kal_dirs)

s2c

#library("biomaRt")
#if (genome_type=="mm10"){
  #mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl",host = 'www.ensembl.org')
#  mart <- try(suppressMessages(biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "mmusculus_gene_ensembl")), silent=T)
#}
#if (genome_type=="hg19"){
#  mart <- try(suppressMessages(biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = 'grch37.ensembl.org')), silent=T)
#}
#t2g=try(suppressMessages(biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name","transcript_version"), mart = mart)),silent=TRUE)   
#if (inherits(t2g,"try-error")){    #if we cannot connect to biomart
  print("cannot connect to BiomaRt")
  t2g=read.table(input.IDconversion,sep="\t",header=T)
#}
t2g$target_id <- t2g$ensembl_transcript_id
t2g[,c("ensembl_transcript_id","transcript_version")] <- list(NULL) # delete the ensembl transcript ID and transcript version columns
t2g <- dplyr::rename(t2g, target_id = target_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
t2g = t2g[,c('target_id', 'ens_gene', 'ext_gene')]


###
# Read in counts
###

library("tximport")
library("readr")

files = paste(s2c$path, 'abundance.tsv', sep='/')
names(files) = s2c$sample

txi = tximport(files, type="kallisto", tx2gene=t2g)

countTable = txi$counts
sampleTable = data.frame(s2c, type="paired-end")
rownames(sampleTable) = s2c$sample
sampleTable$condition = as.factor(sampleTable$condition)

dds <- DESeqDataSetFromMatrix(countData = round(countTable),
                                   colData = sampleTable,
                                   design = ~ condition)

#### run DESeq2 ####
dds=DESeq(dds)
res = results(dds)

idx = match(rownames(res), t2g$ens_gene)
res$ext_gene=t2g$ext_gene[idx]


write.table(res, output, row.names=T, sep="\t", quote=F)

print("finished!")

