#!/usr/bin/env Rscript

library(ggplot2)
library(scales)

white_background<-ggplot2::theme(
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.background = element_blank())

args<-commandArgs(TRUE)

#"/media/marcgabriel/9d38aaa1-ce5d-4422-950f-0f083b797d13/dominika_hXrn1_sRNA/metamapping/all_samples_percent_reads_per_feature.tsv"

#"/media/marcgabriel/9d38aaa1-ce5d-4422-950f-0f083b797d13/dominika_hXrn1_miRNA/metamapping/all_samples_percent_reads_per_feature.tsv"

#data_frame<-read.delim(file="/media/marcgabriel/9d38aaa1-ce5d-4422-950f-0f083b797d13/dominika_hXrn1_miRNA/metamapping/all_samples_percent_reads_all_per_feature.tsv", sep = "\t",check.names = FALSE,header=T)

#data_frame<-read.delim(file="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/counts_distrib_per_features/AllFeat_AllConds_norm_AllReps_AllPerct.tsv", sep = "\t",check.names = FALSE,header=T)


data_frame<-read.delim(file=args[1], sep = "\t",check.names = FALSE,header=T)
home<-args[2]

cond_order<-args[3]

cond_order<-unlist(strsplit(cond_order,","))




#mimic ggplot colors
gg_color_hue <- function(n) {hues = seq(15, 375, length = n + 1);hcl(h = hues, l = 65, c = 100)[1:n]}

#take 15 colors
selected_colors<-gg_color_hue(15)

#select the ones that could be compatible
selected_colors<-selected_colors[c(6,8,11,12)]

#add the ones that could be useful
selected_colors<-c("grey",selected_colors,"#66cc66","orange","gold","red")

#remove the extra suffix
data_frame$feature<-gsub("_combined","",data_frame$feature)

wanted_conds<-sort(as.character(unique(data_frame$feature)))
selected_colors<-gg_color_hue(length(wanted_conds))
names(selected_colors)<-wanted_conds


gg_color_hue <- function(n) {hues = seq(15, 375, length = n + 1);hcl(h = hues, l = 65, c = 100)[1:n]}

#PCG<-gg_color_hue(4)[4]

#if we are in a particular config, adapt the colors
if(any(grepl("protein",data_frame$feature))){

  #"protein_coding"="#C77CFF"
    selected_colors<-c("protein_coding_gene"="purple","protein_coding"="purple","pseudogene"="#ffcc00","lncRNA"="#3eb489","snoRNA"="red","snRNA"="cornflowerblue",
      "sRNA"="cyan","miRNA"="#DF5286","rRNA"="grey40","tRNA"="brown","Murine_leukemia_virus"="darkgreen")

}else if(any(grepl("exon",data_frame$feature))){
  
  selected_colors<-c("exon"="#F8766D","intron"="#00BFC4","promoter"="#00BA38","five_prime_utr"="orange","three_prime_utr"="gold","intergenic"="grey")
  
}

unique_conds<-unique(data_frame$condition)

get_sem <- function(x) sqrt(var(x)/length(x))


final_frame<-data.frame()
for(i in 1:length(unique_conds)){
  
  tmp<-data_frame[which(data_frame$condition==unique_conds[i]),]
  tmp3<-data.frame()
  for(y in 1:length(unique(tmp$feature))){
    
    tmp2<-tmp[which(tmp$feature==unique(tmp$feature)[y]),]
    
    my_sd<-sd(tmp2$counts)
    my_mean<-mean(tmp2$counts)
    my_sem<-get_sem(tmp2$counts)
    
    tmp3<-rbind(tmp3,data.frame(condition=unique_conds[i],feature=unique(tmp$feature)[y],my_mean=my_mean,sd=my_sd,sem=my_sem))
    
  }
  
  final_frame<-rbind(final_frame,tmp3)
  
}

final_frame$label<-""

for(i in 1:nrow(final_frame)){
  
  if(final_frame$my_mean[i]<2){
    
    final_frame$label[i]<-""
    
    
  }else{
    
    final_frame$label[i]<-paste0(round(final_frame$my_mean[i],1),"%")
  }
}

feat_to_remove<-c()
for(i in 1:length(unique(final_frame$feature))){
  
  my_tmp<-data.frame()
  my_tmp<-final_frame[which(final_frame$feature==unique(final_frame$feature)[i]),]
  
  if(all(my_tmp$my_mean==0)==T){
    
    feat_to_remove<-c(feat_to_remove,as.character(unique(final_frame$feature)[i]))
    
  }
}

final_frame<-final_frame[which(!as.character(final_frame$feature)%in%feat_to_remove),]


#order by the features with the strongest mean
final_frame<-final_frame[order(-final_frame$my_mean),]

if(any(grepl("exon",data_frame$feature))){
  
  
  final_frame$feature<- factor(final_frame$feature, levels = c("exon","intron","three_prime_utr","five_prime_utr","intergenic","promoter"))
  
}else{
	
	final_frame$feature<- factor(final_frame$feature, levels = unique(as.character(final_frame$feature)))
	
	
}


final_frame$condition<- factor(final_frame$condition, levels = cond_order)

final_frame<-final_frame[which(final_frame$condition!="EV_NP40"),]


if(length(unique(final_frame$condition))<=4){

  #png(paste(home,"reads_genomic_distrib_model1.png",sep=""),width=500,height=800)
  
  png(paste(home,"reads_genomic_distrib_model1.png",sep=""),width=180,height=180,units="mm",res=1000)
  
}else{
  
  
  #units="mm",res=300,width=183,height=130,,res=300
  png(paste(home,"reads_genomic_distrib_model1.png",sep=""),width=1200,height=800)
  
  
}

  ggplot(data = final_frame, aes(x = condition, y = my_mean, fill = feature )) + 
    geom_bar(stat = "identity",position="stack") +
    #we put some extra values to the limit, if the sum of the percentages is slightly over 100
    scale_y_continuous(limits=c(0,100.2),expand=c(0,0))+
    scale_x_discrete(expand=c(0,0))+
    scale_fill_manual(values = selected_colors)+
    xlab("conditions")+
    guides(fill=guide_legend(title="feature_type"))+
    ylab("counts (%)")+
    geom_text(aes(label = label),position=position_stack(vjust =0.5),size=6)+
  theme(axis.text.x = element_text(size = 16,angle = 45,colour = "black",hjust=1,face="bold"),
        axis.text.y = element_text(size = 16,colour = "black",hjust=0.5,face="bold"),
        legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"),
        plot.title = element_text(hjust = 0.5,size=18),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

dev.off()



if(length(unique(final_frame$condition))<=4){
  
  png(paste(home,"reads_genomic_distrib_model1_NoLabel.png",sep=""),width=500,height=800)
  
}else{
  
  
  #units="mm",res=300,width=183,height=130,,res=300
  png(paste(home,"reads_genomic_distrib_model1_NoLabel.png",sep=""),width=1200,height=800)
  
  
}


ggplot(data = final_frame, aes(x = condition, y = my_mean, fill = feature )) + 
  geom_bar(stat = "identity",position="stack") +
  #we put some extra values to the limit, if the sum of the percentages is slightly over 100
  scale_y_continuous(limits=c(0,100.2),expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  scale_fill_manual(values = selected_colors)+
  xlab("conditions")+
  guides(fill=guide_legend(title="gene_type"))+
  ylab("counts (%)")+
  theme(axis.text.x = element_text(size = 16,angle = 45,colour = "transparent",hjust=1,face="bold"),
        axis.text.y = element_text(size = 16,colour = "transparent",hjust=0.5,face="bold"),
        legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"),
        plot.title = element_text(hjust = 0.5,size=18),
        axis.title=element_text(size=16,face="bold",colour="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

dev.off()

#per feature
# png(paste(home,"reads_genomic_distrib_model2.png",sep=""),units="mm",res=300,width=183,height=130)
# 
#   ggplot(data = final_frame, aes(x = condition, y = my_mean,fill=feature, color = feature )) + 
#     geom_bar(stat = "identity",position="dodge") +
#     geom_errorbar(aes(ymin=my_mean-my_sd, ymax=my_mean+my_sd), width=.5,
#                   position=position_dodge(.9),color="black")+
#     scale_y_continuous(limits=c(0,100),expand=c(0,0))+
#     scale_x_discrete(expand=c(0,0))+
#     ylab("counts percentage")+
#     xlab("")+
#     facet_wrap(~ feature,ncol =2,scales="free")+
#     theme(axis.text.x = element_text(size = 8,angle = 45,colour = "black",hjust=1,face="bold"),
#           axis.text.y = element_text(size = 8,colour = "black",hjust=0.5,face="bold"),
#           strip.text.x = element_text(size = 10))+
#     white_background+
#     theme(legend.position = "none")+
#     geom_vline(xintercept = 0.55)
# 
#   dev.off()

#per condition
# png(paste(home,"reads_genomic_distrib_model3.png",sep=""),units="mm",res=300,width=183,height=130)
#   
#   ggplot(data = final_frame, aes(x = feature, y = my_mean,color=feature, fill = feature )) + 
#     geom_bar(stat = "identity",position="dodge") +
#     geom_errorbar(aes(ymin=my_mean-my_sd, ymax=my_mean+my_sd), width=.5,
#                   position=position_dodge(.9),color="black")+
#     scale_y_continuous(limits=c(0,100),expand=c(0,0))+
#     white_background+
#     scale_fill_manual(values = selected_colors)+
#     scale_color_manual(values = selected_colors)+
#     scale_y_continuous(limits=c(0,100),expand=c(0,0))+
#     scale_x_discrete(expand=c(0,0))+
#     ylab("counts percentage")+
#     xlab("")+
#     facet_wrap(~ condition,ncol =2,scales="free")+
#     theme(axis.text.x = element_text(size = 8,angle = 45,colour = "black",hjust=1,face="bold"),
#           axis.text.y = element_text(size = 8,colour = "black",hjust=0.5,face="bold"),
#           plot.title = element_text(hjust = 0.5,size=18),
#           axis.title=element_text(size=8,face="bold"),
#           strip.text.x = element_text(size = 10))+
#     theme(legend.position = "none")+
#     geom_vline(xintercept = 0.55)
# 
# dev.off()  
  
# #pie chart via ggplot (multiplot)
# ggplot(data = final_frame, aes(x = "", y = my_mean, fill = feature )) + 
#   geom_bar(stat = "identity", position = position_fill()) +
#   geom_text(aes(label = my_mean), position = position_fill(vjust = 0.5)) +
#   coord_polar(theta = "y") +
#   facet_wrap(~ condition,nrow =3)  +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank()) + 
#   theme(legend.position='bottom') + 
#   guides(fill=guide_legend(nrow=2, byrow=TRUE))





  
