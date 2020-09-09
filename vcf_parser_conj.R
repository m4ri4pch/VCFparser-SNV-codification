# Analisis variantes en genes 
# Parseo de VCFs, anotacion y priorizacion
# MPC
# Tue Jan 28 11:40:46 2020


fileList<-list.files(path="./data",pattern="zip",all.files=T, full.names = T, include.dirs=T,recursive=T)
r<-NULL
for(i in 1:length(fileList)){
  # if(i=1){
  #   print("0%")
  #   }else if(i=length(fileList)){
  #   print("=100%")
  #   }else{
      cat("=")
  #   }

  x<-fileList[i]  
  d<-gsub(".zip","",x)
  dir.create(d)
  unzip(x,exdir = d)
  v<-read.delim(list.files(list.dirs(file.path(d,"Variants")),full.names = T,pattern=".tsv"),header=T,fill=T,stringsAsFactors = F)
  if(nrow(v)==0) next()
  vars<-paste(v[,1],v[,3],sapply(strsplit(v[,5],"/"),"[[",2),sep=":")
  gen<-as.numeric(sapply(strsplit(v[,5],"/"),function(x) x[1]==x[2]))
  sam<-rep(gsub("./.*/","",d),times=length(gen))
  c<-data.frame("variant"=vars,"genotype"=gen,"sample"=sam,stringsAsFactors = F)
  r<-rbind(r,c)
}
# 0 = het
# 1 = homo
dim(r)
# 25428 variantes

sum(duplicated(r$variant))
# 24893

length(unique(r$variant))
# 535

variants<-unique(r$variant)

m<-data.frame("variants"=variants,stringsAsFactors = F)

l<-lapply(split(r,r$sample), function(x) x[,c(1,2)])
s<-sapply(l,function(x) as.numeric(m$variants%in%x[,1]))
rownames(s)<-m$variants
totalV<-colSums(s)

dim(s)
# 535 735

x<-l[[1]]
lHom<-lapply(l,function(x) x[x$genotype==1,])

sHom<-sapply(lHom,function(x) as.numeric(m$variants%in%x[,1]))
rownames(sHom)<-m$variants

dim(sHom)
# 535 735

sVar<-s+sHom
write.table(sVar,file="variants_matrix_total.txt",sep="\t",col.names = T,row.names = T)

# Anotacion de las variantes con CellBase

variants<-rownames(sVar)

urls<-paste("http://cellbase.clinbioinfosspa.es/cb/webservices/rest/v4/hsapiens/genomic/variant/",gsub("chr","",variants),
"/annotation?assembly=GRCh37",sep ="")

library(RCurl)
library(rjson)
p<-as.list(variants)
names(p)<-variants
for(i in 1:length(variants)){
  aux<-fromJSON(getURL(urls[i]))
  p[[i]]<-aux$response[[1]]$result
}

annotCB<-function(x){

  var=paste(x[[1]]$chromosome,x[[1]]$start,x[[1]]$reference,x[[1]]$alternate,sep=":")
  rs=x[[1]]$id
  if(is.null(rs)) rs<-NA
  hgvs=paste(x[[1]]$hgvs,collapse=", ")
  if(is.null(hgvs)) hgvs<-NA
  genes=paste(unique(unlist(sapply(x[[1]]$consequenceTypes,function(x) x$geneName))),collapse=", ")
  if(is.null(genes)) genes<-NA
  consequence=x[[1]]$displayConsequenceType
  if(is.null(consequence)) consequence<-NA
  scaledCADD=x[[1]]$functionalScore[[grep("cadd_scaled",x[[1]]$functionalScore)]]$score
  if(is.null(scaledCADD)) scaledCADD<-NA
  sift=try(min(unlist(lapply(x[[1]]$consequenceTypes, function(x) x$proteinVariantAnnotation$substitutionScores[[grep("sift",x$proteinVariantAnnotation$substitutionScores)]]$score))))
  if(sift=="Inf") sift<-NA
  polyphen=try(max(unlist(lapply(x[[1]]$consequenceTypes, function(x) x$proteinVariantAnnotation$substitutionScores[[grep("polyphen",x$proteinVariantAnnotation$substitutionScores)]]$score))))
  if(polyphen=="-Inf") polyphen<-NA
  EURf=try(x[[1]]$populationFrequencies[grep("EUR",x[[1]]$populationFrequencies)][[1]]$altAlleleFreq)
  if(is.null(EURf)) EURf<-NA
  EURhetf=try(x[[1]]$populationFrequencies[grep("EUR",x[[1]]$populationFrequencies)][[1]]$hetGenotypeFreq)
  if(is.null(EURhetf)) EURhetf<-NA
  EURhomf=try(x[[1]]$populationFrequencies[grep("EUR",x[[1]]$populationFrequencies)][[1]]$altHomGenotypeFreq)
  if(is.null(EURhomf)) EURhomf<-NA
  gerp=x[[1]]$conservation[grep("gerp",x[[1]]$conservation)][[1]]$score
  if(is.null(gerp)) gerp<-NA
  phylop=x[[1]]$conservation[grep("phylop",x[[1]]$conservation)][[1]]$score
  if(is.null(phylop)) phylop<-NA
  phastCons=x[[1]]$conservation[grep("phastCons",x[[1]]$conservation)][[1]]$score
  if(is.null(phastCons)) phastCons<-NA
  z<-cbind(var,rs,genes,hgvs,consequence,scaledCADD,sift,polyphen,EURf,EURhetf,EURhomf,gerp,phylop,phastCons)
  names(z)<-c("Variant","id","Gene","hgvs","Consequence","scCADD","sift","polyP","EUR_freq","EUR_hetFreq","EUR_homFreq","gerp","phyloP","PhastCons")
  
  return(z)
}
z<-sapply(p,annotCB)
z<-t(z)
z[grep("Error",z)]<-NA
all(rownames(z)==rownames(sVar))

sVarAnn<-cbind(sVar,z)
write.table(sVarAnn,"variants_matrix_total_Annotated.txt",col.names = T,row.names = T, sep="\t",quote = F)
