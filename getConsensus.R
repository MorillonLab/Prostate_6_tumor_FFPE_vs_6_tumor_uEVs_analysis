#!/usr/bin/env Rscript
#https://www.rdocumentation.org/packages/seqinr/versions/3.4-5/topics/consensus

args<-commandArgs(TRUE)

suppressPackageStartupMessages(library(seqinr))

suppressPackageStartupMessages(library(Biostrings))

consensus_name=args[1]

#aln_file="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/circant_results/circant_results_out/FFPE_T_1/circ/titi.fa"

#aln_file="/home/marcgabriel/Desktop/refine_dekupl_contigs/consensus_aln_refined.txt"
aln_file=args[2]

one_strand=args[3]

#ref_seq="/home/marcgabriel/Desktop/refine_dekupl_contigs/consensus_ref.fa"
ref_seq=args[4]

if(ref_seq!="no"){
  
  ref_aln<-read.alignment(file=ref_seq, format="fasta", forceToLower = FALSE)
  ref_aln<-unlist(strsplit(ref_aln$seq,""))

}else{
  
  ref_aln<-""
}


aligned_contigs<-read.alignment(file=aln_file, format="fasta", forceToLower = FALSE)

#profile : you have the score here for each nucleotide
#the matrix is like this :
#  1 2 3 4...
#- 0
#a 4
#t 0
#g 0
#c 0
#we have the number of sequences with the nucleotide at the tested position
res<-consensus(matali=aligned_contigs,method="profile")

final_consensus<-c()

for(i in 1:ncol(res)){
  
  #select the one(s) with the highest frequency (we don't count the gaps)
  my_max=max(res[,i][names(res[,i])!="-"])
  
  #if there's at least one nuc with a non-zero score, take it/them
  #
  if(my_max!=0){
    
    #take the corresponding nucleotides
    selected_nuc<-names(res[,i][res[,i]==my_max])
  
  #otherwise, give a gap  
  }else{
    
    selected_nuc<-"-"
    
  }
  

  
  #if many nucleotides are proposed, remove the one with the dash
  if(length(selected_nuc)!=1){
    
    selected_nogap<-selected_nuc[selected_nuc!="-"]
    
    #if we still have many propositions, take the reference at the given position (should we take a random one ?), if the ref file has been set
    if(length(selected_nogap)>1){
      
      if(ref_aln!=""){
      
        final_consensus<-c(final_consensus,ref_aln[i])
      
        #if the ref file isn't set, take the first one
      }else{
        
        final_consensus<-c(final_consensus,selected_nogap[1])
        
      }
    
    #otherwise, if we only have one proposition, take it  
    }else if(length(selected_nogap)==1){
      
      final_consensus<-c(final_consensus,selected_nogap)
    
    #otherwise, if we have no proposition after the removal of the gaps, that means the gaps were the only propositions, just just add it !  
    }else if(length(selected_nogap)==0){
      
      
      final_consensus<-c(final_consensus,"-")
      
    }
  
  #if only one proposition, just take it  
  }else{
    
    final_consensus<-c(final_consensus,selected_nuc)
  }
  
}

final_consensus<-toupper(final_consensus)

#remove the gaps
final_consensus<-final_consensus[final_consensus!="-"]

#remove the spaces
final_consensus<-final_consensus[final_consensus!="" ]



if(one_strand=="+"){
  
  cat(paste(">",consensus_name,sep=""),"\n",paste(final_consensus,collapse="",sep=""),"\n",sep="")

#if it's the minus strand, rev-comp the sequence
}else{
  
  cat(paste(">",consensus_name,sep=""),"\n",as.character(reverseComplement(DNAString(paste(final_consensus,collapse="",sep="")))),"\n",sep="")
  
}

