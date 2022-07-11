#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

if (length(args)==0){
  stop("missing arguments !\nUsage : ./getIntergenics.R <gff3 annotation>")
}

suppressMessages(library(GenomicFeatures))

suppressMessages(library(systemPipeR))

file<-args[1]

txdb<-suppressMessages(makeTxDbFromGFF(file,format = "gff"))

feat<-unlist(suppressMessages(genFeatures(txdb, featuretype="intergenic",reduce_ranges = T)))

intergenicGFF<-data.frame(chromosome=rep(as.character((feat@seqnames@values)),feat@seqnames@lengths),
                       source="annotation",
                       feature="intergenic",
                       start=unlist(start(feat)),
                       end=unlist(end(feat)),
                       empty1=".",
                       strand=".",
                       empty2=".",
                       ID=paste("ID=intergenic_",rep(as.character((feat@seqnames@values)),feat@seqnames@lengths),"_",unlist(start(feat)),"_",unlist(end(feat)),";gene_name=intergenic_",rep(as.character((feat@seqnames@values)),feat@seqnames@lengths),"_",unlist(start(feat)),"_",unlist(end(feat)),";gene_id=intergenic_",rep(as.character((feat@seqnames@values)),feat@seqnames@lengths),"_",unlist(start(feat)),"_",unlist(end(feat)),";transcript_id=intergenic_",rep(as.character((feat@seqnames@values)),feat@seqnames@lengths),"_",unlist(start(feat)),"_",unlist(end(feat)),sep=""))

write.table(intergenicGFF,stdout(),sep='\t',row.names=F,col.names=F,quote=F)
