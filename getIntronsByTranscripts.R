#!/usr/bin/env Rscript

#check this script for this issue : gene_id & transcrit finishing with "_1" (or globally "_X"), are truncated from the "_" -> gencode27lift37 for example

args <- commandArgs(TRUE)

if (length(args)==0){
  stop("missing arguments !\nUsage : ./getIntronsByTranscripts.R <gff3 annotation> <temp dir>")
}

suppressMessages(library(GenomicFeatures))
suppressMessages(library(data.table))
suppressMessages(library(rtracklayer))
library(data.table)
library(doParallel)
library(foreach)

file<-args[1]

#file<-"/home/marcgabriel/Documents/gencode27/gencode.v27.annotation_sorted.gff3"
#file<-"/media/marcgabriel/Transcend/Rocco_annotations/my_test/tmp_process/test_2.gff"
#file<-"/media/marcgabriel/Transcend/Rocco_annotations/my_test/test.gff"
#file<-"/media/marcgabriel/Transcend/Rocco_annotations/my_test/2_lvls.gff"
#file<-"/media/marcgabriel/saylar8/Rocco_cut_and_run_analysis/LNCaP_scallop_gencodeV32_all_lvl_02Ago2021_fixed_v2.gff3"

temp_dir=args[2]
#temp_dir<-"/media/marcgabriel/saylar8/Rocco_cut_and_run_analysis/"

temp_dir=paste(temp_dir,"/",sep="")
temp_dir<-gsub("//","/",temp_dir)
dir.create(temp_dir,showWarnings = F)



options(ucscChromosomeNames=F)

#be carefull, you need the parent having as feature "mRNA", "transcript", or "gene"
#The ID of the parent should be the same as the attribute of "Parent" of the exon
#file="/media/marcgabriel/homeborn/cillance_26meta_no_vs_31_meta_yes_dekupl_results/annotation_1st_version_final/diff_contigs_exon_mRNA.gff3"
#file="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/ving_genes/gene_exon_lvl.gff"
#issue with paraloguous genes (= same gene_id, same gene_name, but on different chromosomes...), so we will use a loop on the chromosomes

# ensure that transcript level comes before exon lvl (sort -k1,1 -k4,4n -k5,5nr -k3,3r)

#file="/home/marcgabriel/Documents/gencode32/gencode.v32.annotation_sorted.gff3"
all_chrs<-system(paste("cut -f1 ",file,"|grep -v ^# |sort -u",sep=""),intern = T)


files_to_delete<-c()

registerDoParallel(cores=detectCores())

all_sub_frames<-invisible(foreach(one_pos=1:length(all_chrs), .combine="c", .inorder=FALSE) %dopar% {
#for(one_chr in all_chrs){
  
      one_chr<-all_chrs[one_pos]
      chrom_gff<-paste(temp_dir,one_chr,"_sub_gff_for_intron.gff",sep="")
      sub_intron_gff<-paste(temp_dir,one_chr,"_sub_intron_for_intron.gff",sep="")
  
      system(paste("grep -P \"",one_chr,"\\t\" ",file, " >",chrom_gff,sep=""))

      txdb<-makeTxDbFromGFF(chrom_gff,format="gff3")
      
      IntronsByTranscripts<-unlist(intronsByTranscript(txdb,use.names=T))
      
      if(length(IntronsByTranscripts)!=0){
        
        
      
          #put as "ID" the Parent
          genecode<-readGFF(file, version=3, columns=NULL, tags=NULL, filter=NULL, nrows=-1, raw_data=FALSE)
          
          genecode<-genecode[which(genecode$type=="exon"),]
          
          #parent can be characterList, so put it in character
          genecode$Parent<-as.character(genecode$Parent)
          
          genecode$ID<-genecode$Parent
          
          
          
          if(!any(grepl("Parent",names(genecode))==T)){
            
            stop("no parent level in the gff !")
            
          }else{
            
            genecode$transcript_id<-genecode$Parent
            
            #makeTxDbFromGFF don't use ENSG00000181222.16_[0-9]+ as, it is the case for the backmapped annotations, we have to be consistent with it
            
            if(any(grepl("ENSG[0-9]+\\.[0-9]+_[0-9]+",genecode$transcript_id))){
              
              genecode$crossing_id<-gsub("_[0-9]+$","",genecode$transcript_id)
            
            }else{
              
              genecode$crossing_id<-genecode$transcript_id
            }
            
            #we concatenate the ID, the chr, and the strand, to take into account the paralogous genes
            #genecode$ID<-paste(genecode$ID,"_",genecode$seqid,"_",genecode$strand,sep="")
            genecode$ID<-genecode$Parent
            
            #vector that will store column names with the type "characterList" (we don't want them)
            col_to_remove=c()
            for(i in 9:(ncol(genecode)-1)){
              
              #if the column isn't a character, don't process it
              if(class(genecode[,i])!="CompressedCharacterList"){
              
                genecode[,i]<-paste(as.character(names(genecode[i])),"=",as.character(unlist(genecode[,i])),";",sep="")
              
              }else{
                
                col_to_remove<-c(col_to_remove,names(genecode[i]))
              }
            
              
            }
            
            
            #remove exon id and exon number, and columns that are "characterList"
            genecode<-genecode[,names(genecode)[!names(genecode)%in%c("exon_id","exon_number",col_to_remove)]]
            
            genecode<-genecode[,9:ncol(genecode)]
            
            genecode<-genecode[!duplicated(genecode), ]
            
            genecode<-genecode[,c(ncol(genecode),1:(ncol(genecode)-1))]
            
            
            genecode$attributes<- apply( genecode[,2:ncol(genecode)] , 1 , paste , collapse = "" )
            
            intronsGFF<-data.frame(transcript_id=names(IntronsByTranscripts),
                                   chromosome=seqnames(IntronsByTranscripts),
                                   source="annotation",
                                   feature="intron",
                                   start=start(IntronsByTranscripts),
                                   end=end(IntronsByTranscripts),
                                   empty1=".",
                                   strand=strand(IntronsByTranscripts),
                                   empty2=".")
            
            intronsGFF<-merge(intronsGFF,genecode[,c("crossing_id","attributes")],
                              by.x="transcript_id",
                              by.y="crossing_id",
                              all.x=TRUE,
                              all.y=FALSE)
            
            intronsGFF<-subset(intronsGFF,select=-c(transcript_id))
            
            #stdout()
            write.table(intronsGFF,sub_intron_gff, sep='\t',row.names=F,col.names=F,quote=F)
            
            file.remove(chrom_gff)
            
            return(sub_intron_gff)
            
          }
      
      }else{invisible(file.remove(chrom_gff))}

}
)


if(length(all_sub_frames)>0 & 
   all(is.logical(all_sub_frames)==F)){
  

final_frame<-data.frame()
for(i in 1:length(all_sub_frames)){
  
  one_file<-all_sub_frames[i]
  
  if(!is.logical(one_file) & file.exists(one_file)){
    
    final_frame <-rbind(final_frame,fread(one_file, sep = "\t",check.names = FALSE,header=F,data.table=F))
    
    invisible(file.remove(one_file))
    
  }
}

write.table(final_frame,stdout(), sep='\t',row.names=F,col.names=F,quote=F)

}else{
  
  cat("no results...")
  
}



