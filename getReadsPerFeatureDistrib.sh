#!/bin/bash


##### input variables #######################

export annotation="/home/marcgabriel/Documents/gencode32/gencode.v32.annotation_sorted_official_chromosomes.gff3"

export getIntronsByTranscripts="/home/marcgabriel/Desktop/scripts/getIntronsByTranscripts.R"

#create intergenic features
export getIntergenics="/home/marcgabriel/Desktop/scripts/getIntergenics.R"

export getMergedExonsPerGenes="/home/marcgabriel/Desktop/scripts/getMergedExonsPerGenes.sh"

export getTranscriptByExonsByStrand="/home/marcgabriel/Desktop/scripts/getTranscriptByExonsByStrand.sh"

export getClassifiedReads="/home/marcgabriel/Desktop/scripts/getClassifiedReads2.R"

export getClassifiedReadsAllConds="/home/marcgabriel/Desktop/scripts/getClassifiedReadsAllConds.R"

export genomic_distrib_stats="/home/marcgabriel/Desktop/scripts/genomic_distrib_stats_metagenome.R"

export bedtools="bedtools"

export genome_index="/home/marcgabriel/Documents/gencode32/GRCh38.primary_assembly_genome_official_chromosomes.fa.fai"

export bedops="/home/marcgabriel/Downloads/bedops-2.4.30/bin/bedops"

export samtools="samtools"

export featureCounts_prog="/home/marcgabriel/Downloads/subread-1.6.0-source_MAX_HIT_NUMBER_10e6/bin/featureCounts"

export threads=12

export lib_type="fr-firststrand"
#export lib_type="fr-secondstrand"

export paired_end="yes"
#export paired_end="no"

config_file="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/exon_intron_design_otherProstateCelllines.tsv"

output_dir="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/genes_expression2/counts_distrib_per_feature/"


################################################

#set -x

output_dir="${output_dir}/"
export output_dir=$(echo $output_dir | sed 's/\/\//\//g')

if [ ! -d $output_dir ];then mkdir $output_dir ;fi

export design_file="${output_dir}design.txt"
if [ -f $design_file ];then rm -rf $design_file;fi

cut -f1,2 $genome_index |LANG=en_EN sort -k1,1 -k2,2n >${output_dir}chrs_length.tsv

readarray my_config < $config_file

#all_conditions=($(echo "${my_config[*]}"|sed 's/^ //g'|grep -v "^$"|awk '{print $2}' |LANG=en_EN sort -u))


##look for all conds ; keep the order in the file
all_conditions=($(echo "${my_config[*]}"|sed 's/^ //g'|grep -v "^$"|awk '{print $2}' |awk '!seen[$0]++'))

echo -e "all conds in order are : ${all_conditions[*]}\n"


if [[ "$lib_type" == "fr-firststrand" ]];then

   orientation="-s 2"

elif [[ "$lib_type" == "fr-secondstrand" ]];then
  
  orientation="-s 1"

else

  orientation=""

fi

if [[ "$paired_end" == "yes" ]];then

	featurecount_pairs="-p"
  
else

	featurecount_pairs=""

fi
  
  
 


###############  dev #################


getCleanIntervals(){

  gff_to_remove=$1

  gff_to_extract=$2
  
  prefix=$3
  
  file_name=$4
  
  my_output=$5
  
  use_prefix_as_ID=$6
  
  #set -x
  
  
  if [ "$prefix" != "intergenic" ];then
  
  
     strandness="-s"
  
  
  else
  
  
     strandness=""
  
  
  fi
  
  if [ "$gff_to_remove" != "" ];then
  
        
  
           cat $gff_to_extract |LANG=en_EN sort -k1,1 -k4,4n| $bedtools merge -nonamecheck $strandness -c 7 -o "distinct" -i stdin |awk -v strandness="$strandness" 'OFS="\t"{if(strandness=="-s"){strand=$4}else{strand="."};print $1,$2,$3,".",".",strand}' >${my_output}gff_to_extract.bed
           
           cat $gff_to_remove |LANG=en_EN sort -k1,1 -k4,4n| $bedtools merge -nonamecheck $strandness -c 7 -o "distinct" -i stdin |awk -v strandness="$strandness" 'OFS="\t"{if(strandness=="-s"){strand=$4}else{strand="."};print $1,$2,$3,".",".",strand}' >${my_output}gff_to_remove.bed
   
		   if [ "$use_prefix_as_ID" == "T" ];then
		   
		   
			   #remove b from a
			   bedtools subtract $strandness -nonamecheck -a ${my_output}gff_to_extract.bed -b ${my_output}gff_to_remove.bed  |LANG=en_EN sort -k1,1 -k2,2n|$bedtools merge -nonamecheck $strandness -c 6 -o "distinct" -i stdin|awk -v prefix=$prefix 'OFS="\t"{print $1,"empty",prefix,$2+1,$3,".",$4,".","ID="prefix}' >${my_output}${file_name} && rm ${my_output}gff_to_extract.bed ${my_output}gff_to_remove.bed
			   
		   else
		   
			   #remove b from a
			   bedtools subtract $strandness -nonamecheck -a ${my_output}gff_to_extract.bed -b $gff_to_remove  |LANG=en_EN sort -k1,1 -k2,2n |$bedtools merge -nonamecheck $strandness -c 6 -o "distinct" -i stdin|awk -v prefix=$prefix 'OFS="\t"{print $1,"empty",prefix,$2+1,$3,".",$4,".","ID="NR}' >${my_output}${file_name} && rm ${my_output}gff_to_extract.bed ${my_output}gff_to_remove.bed
			   
		   fi
	  	  
  else
  
		   if [ "$use_prefix_as_ID" == "T" ];then
		 

			   #just merge the features (we create a tmp file, to avoid the risk input = output)
			   LANG=en_EN sort -k1,1 -k4,4n $gff_to_extract >${my_output}gff_to_extract.tmp && $bedtools merge -nonamecheck $strandness -c 7 -o "distinct" -i ${my_output}gff_to_extract.tmp |awk -v prefix=$prefix 'OFS="\t"{print $1,"empty",prefix,$2+1,$3,".",$4,".","ID="prefix}' >${my_output}${file_name} && rm ${my_output}gff_to_extract.tmp
			 
			   
		   else
		   
		   
		   
			   #just merge the features (we create a tmp file, to avoid the risk input = output)
			   LANG=en_EN sort -k1,1 -k4,4n $gff_to_extract >${my_output}gff_to_extract.tmp && $bedtools merge -nonamecheck $strandness -c 7 -o "distinct" -i ${my_output}gff_to_extract.tmp |awk -v prefix=$prefix 'OFS="\t"{print $1,"empty",prefix,$2+1,$3,".",$4,".","ID="NR}' >${my_output}${file_name} && rm ${my_output}gff_to_extract.tmp
			   
		   fi  
  
  fi
  

}

############### end of dev #############

#list of features to search
#exon_number=1
features_list=(exon promoter intron intergenic five_prime_utr three_prime_utr)


start=$(date)

if [ ! -f ${output_dir}ref_annotation_refined.gff ] || [[ $(wc -l ${output_dir}ref_annotation_refined.gff|awk '{print $1}') -lt 2 ]];then
### dev #####

	##construct the metatranscript (with exons) for ving
	$getMergedExonsPerGenes -a $annotation -o ${output_dir} >${output_dir}meta_transcript.gff

	##construct the gene lvl of the metatranscript
	$getTranscriptByExonsByStrand -a ${output_dir}meta_transcript.gff -f "gene" -u "yes" | sort -k1,1 -k4,4n >${output_dir}gene_lvl.gff

	##concatenate gene lvl and exon lvl in order to produce the introns
	##the highest lvl should be before the exons, otherwise, the exons will ot be connected...
	cat ${output_dir}gene_lvl.gff ${output_dir}meta_transcript.gff |sort -k1,1 -k4,4n -k5,5nr -k3,3r >${output_dir}gene_exon_lvl.gff
	
	##construct the introns
	$getIntronsByTranscripts ${output_dir}gene_exon_lvl.gff ${output_dir} >${output_dir}introns.gff

	##concatenate the metatranscript + the introns
	##remove artefact introns (those from snoRNAs, Y_RNA, etc)
	cat $annotation ${output_dir}introns.gff|sort -k1,1 -k4,4n|awk 'OFS="\t"{if($3!~/intron/){print}else{if($9!~/=SNOR[0-9]+;/ && $9!~/=SNORA[0-9]+;/ && $9!~/=snoRNA[0-9]+;/ && $9!~/=U[0-9]+;/ && $9!~/=Y_RNA;/){print}}}' |sort -k1,1 -k4,4n >${output_dir}ref_annotation_refined.gff
	
	$bedtools complement -i ${output_dir}ref_annotation_refined.gff -g ${output_dir}chrs_length.tsv|awk 'OFS="\t"{print $1,".","intergenic",$2+1,$3,".",".",".","ID=intergenic_"$2+1"_"$3}' |sort -k1,1 -k4,4n >>${output_dir}ref_annotation_refined.gff

### end of dev ###
fi

sort -k1,1 -k4,4n ${output_dir}ref_annotation_refined.gff >${output_dir}ref_annotation_refined.tmp && mv ${output_dir}ref_annotation_refined.tmp ${output_dir}ref_annotation_refined.gff


#create header from one bam in order to give it to the features, converted in bam files
one_bam_file=($(echo "${my_config[*]}"|sed 's/^ //g'|grep -v "^$"|head -n1|awk '{print $1}'))

$samtools view -H $one_bam_file |grep -E -v "^@PG|^@HD" | awk 'OFS="\t"{print $2,$3}' | sed 's/SN\://;s/LN\://g' | LANG=en_EN sort -k1,1 >${output_dir}genome_size.txt

$samtools view -H $one_bam_file >${output_dir}header.txt

genome_size=$(awk 'BEGIN{a=0}''{a=a+$2}''END{print a}' ${output_dir}genome_size.txt)

echo -e "\ngenome size is $genome_size\n"




declare -A link_feat_to_gff=(

[exon]="${output_dir}exon.gff"
[promoter]="${output_dir}promoter.gff"
[intron]="${output_dir}intron.gff"
[five_prime_utr]="${output_dir}five_prime_utr.gff"
[three_prime_utr]="${output_dir}three_prime_utr.gff"
[intergenic]="${output_dir}intergenic.gff"


)

all_features_gff=()

#construct the design file
for one_feat in ${features_list[*]};do

  echo -e "\n**** feature is :  $one_feat ****\n"

  if [ "$one_feat" != "promoter" ] ;then
  
  
         feat=$one_feat
  
         grep -P -i "${one_feat}\t" ${output_dir}ref_annotation_refined.gff|awk -v feat=$feat 'OFS="\t"{$3=feat;print}'|LANG=en_EN sort -k1,1 -k4,4n >${output_dir}${one_feat}.gff
         
         

         
     echo -e "${output_dir}${one_feat}.gff\t$one_feat\t${output_dir}${one_feat}.bam" >>$design_file

  elif [ $one_feat == "promoter" ];then
  
    
       #$5=a+500
       #$5-500
       grep -P -i "gene\t" ${output_dir}ref_annotation_refined.gff| awk 'OFS="\t"{if($7=="+"){a=$4-1;$4=$4-1001;if($4<=0){$4=1};$5=a}else if($7=="-"){$4=$5+1;if($4<=0){$4=1};$5=$5+1001};$2="promoter";print}' |LANG=en_EN sort -k1,1 -k4,4n >${output_dir}${one_feat}.gff
       
     
     echo -e "${output_dir}${one_feat}.gff\t$one_feat\t${output_dir}${one_feat}.bam" >>$design_file
     
     all_features_gff+=(${output_dir}${one_feat}.gff)


  fi
  
 
done


echo -e "\ntry to reduce features...\n"


#create the features that will be used to avoid overlapping
cat "${output_dir}five_prime_utr.gff" "${output_dir}three_prime_utr.gff" >${output_dir}remove_in_exon.gff

cat ${output_dir}exon.gff ${output_dir}five_prime_utr.gff ${output_dir}three_prime_utr.gff >${output_dir}remove_in_intron.gff

cat ${output_dir}five_prime_utr.gff ${output_dir}three_prime_utr.gff ${output_dir}exon.gff ${output_dir}intron.gff >${output_dir}remove_in_promoter.gff

#create non overlapping features

getCleanIntervals ${output_dir}remove_in_intron.gff ${output_dir}intron.gff "intron" "intron.gff" ${output_dir} "F"

getCleanIntervals ${output_dir}remove_in_exon.gff ${output_dir}exon.gff "exon" "exon.gff" ${output_dir} "F"

getCleanIntervals ${output_dir}remove_in_promoter.gff ${output_dir}promoter.gff "promoter" "promoter.gff" ${output_dir} "F"

getCleanIntervals ${output_dir}promoter.gff ${output_dir}intergenic.gff "intergenic" "intergenic.gff" ${output_dir} "F"

getCleanIntervals "" ${output_dir}five_prime_utr.gff "five_prime_utr" "five_prime_utr.gff" ${output_dir} "F"

getCleanIntervals "" ${output_dir}three_prime_utr.gff "three_prime_utr" "three_prime_utr.gff" ${output_dir} "F"


#parallelize this
for i in $(seq 0 $((${#all_conditions[*]}-1)));do

 one_condition=${all_conditions[$i]}

 #echo -e "\n====== working on condition : $one_condition  ======\n"

 bam_list=($(echo "${my_config[*]}"|sed 's/^ //g'|grep -v "^$"|awk -v condition="$one_condition" '{if($2==condition){print}}'|awk '{print $1}'))

 #echo -e "\tbam list : ${bam_list[*]}\n"

 rep_number=($(echo "${my_config[*]}"|sed 's/^ //g'|grep -v "^$"|awk -v condition="$one_condition" '{if($2==condition){print}}'|wc -l))

 echo -e "- condition : $one_condition\n"
 
      all_counts_all_feat_from_one_cond=()
      all_featLength_from_one_cond=()
 
     for l in $(seq 0 $((${#features_list[*]}-1)));do
     
         feat_name=${features_list[$l]}
         
         if [[ ! -f ${output_dir}${feat_name}_for_featureCounts.gff ]];then
         
         
                 #echo -e "\t- condition : no gff for counting for $feat_name, next\n"
				 paste <(cut -f1-8 ${link_feat_to_gff[$feat_name]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk -v feat="$feat_name" 'OFS="\t"{print "ID="feat"_combined;"}' ${link_feat_to_gff[$feat_name]}) >${output_dir}${feat_name}_for_featureCounts.gff
				 
			 
				 #allow multi-overlap for the exons (priority=1)
				 if [[ "$feat_name" == "exon" ]];then
				 
					   multi_overlap="-O"
				   
				   
				 #put the promoters with the UTRs (one same ID for all the promoters part, one same ID for all UTRs part, and so on for the exons)
				 #don't allow multi-overlap with the UTRs & exons
				 elif [[ "$feat_name" == "promoter" ]];then
				 
				 
				   
					   paste <(cut -f1-8 ${link_feat_to_gff["five_prime_utr"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=five_prime_utr""_combined;"}' ${link_feat_to_gff["five_prime_utr"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
					   
					   paste <(cut -f1-8 ${link_feat_to_gff["three_prime_utr"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=three_prime_utr""_combined;"}' ${link_feat_to_gff["three_prime_utr"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
					   
					   paste <(cut -f1-8 ${link_feat_to_gff["exon"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=exon""_combined;"}' ${link_feat_to_gff["exon"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
							   
					   multi_overlap=""
				   
				 #put the introns with the exons & UTRs (one same ID for all the introns parts, one same ID for all exons parts)
				 #don't allow multi-overlap		 
				 elif [[ "$feat_name" == "intron" ]];then
				 
				 
						paste <(cut -f1-8 ${link_feat_to_gff["exon"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=exon""_combined;"}' ${link_feat_to_gff["exon"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
						
					   paste <(cut -f1-8 ${link_feat_to_gff["five_prime_utr"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=five_prime_utr""_combined;"}' ${link_feat_to_gff["five_prime_utr"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
					   
					   paste <(cut -f1-8 ${link_feat_to_gff["three_prime_utr"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=three_prime_utr""_combined;"}' ${link_feat_to_gff["three_prime_utr"]}) >>${output_dir}${feat_name}_for_featureCounts.gff						
				 
						multi_overlap=""
				  
				  #for the UTRs, add the exons (one same ID for all the UTRs parts, one same ID for all exons parts)
				  #don't allow multi-overlap     
				 elif [[ "$feat_name" == "five_prime_utr" ]] || [[ "$feat_name" == "three_prime_utr" ]];then
				 
				 
					   paste <(cut -f1-8 ${link_feat_to_gff["exon"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=exon""_combined;"}' ${link_feat_to_gff["exon"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
				 
						multi_overlap=""
				 
					
				 #for the intergenic, add the exons, the UTRs, the promoters
				 #don't allow multi-overlap
				 elif [[ "$feat_name" == "intergenic" ]];then
				 
					   paste <(cut -f1-8 ${link_feat_to_gff["five_prime_utr"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=five_prime_utr""_combined;"}' ${link_feat_to_gff["five_prime_utr"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
					   
					   paste <(cut -f1-8 ${link_feat_to_gff["three_prime_utr"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=three_prime_utr""_combined;"}' ${link_feat_to_gff["three_prime_utr"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
					   
					   paste <(cut -f1-8 ${link_feat_to_gff["exon"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=exon""_combined;"}' ${link_feat_to_gff["exon"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
					   
					   paste <(cut -f1-8 ${link_feat_to_gff["promoter"]}|awk 'OFS="\t"{$3="meta_feature";print}') <(awk 'OFS="\t"{print "ID=promoter""_combined;"}' ${link_feat_to_gff["promoter"]}) >>${output_dir}${feat_name}_for_featureCounts.gff
					   
					   multi_overlap=""
				 
				 
				 fi
		 
		 fi
		 
		 
		 		 
		 my_annot="${output_dir}${feat_name}_for_featureCounts.gff"
		 
		 
		 echo -e "\t- feat : $feat_name, file : $my_annot\n"
		 
		 
		 #if [[ ! -f ${output_dir}${feat_name}_${one_condition}_counts_AllReps.txt ]];then
		 
		 
				 all_rep_one_cond=()
				 
				 #initialize rep number
				 k=1
			 
				 for j in ${bam_list[*]};do
				 
					echo -e "\t\t- rep : $k\n"
					
					if [[ ! -f ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.txt ]] || [[ ! -f ${output_dir}${feat_name}_${one_condition}_featLength_rep_${k}.txt ]];then
					
						#for intergenics, don't put the orientation
						if [[ "$feat_name" == "intergenic" ]];then
										
						#-L $orientation
						  $featureCounts_prog -F "GTF" -t "meta_feature" -g "ID" $featurecount_pairs -T $threads $multi_overlap -a $my_annot -o ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.txt $j
						  
						  cp ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.txt ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}_toGetLength.txt
						
						#take the length of the features with reads (only the parts with reads are taken into account)
						feat_length=$(grep -v "^#" ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}_toGetLength.txt|grep "${feat_name}_combined"|awk 'BEGIN{a=0}''{if($7>0){a=a+$6}}''END{print a}')
					
					   else
					   
						  $featureCounts_prog -F "GTF" -t "meta_feature" -g "ID" $orientation $featurecount_pairs -T $threads $multi_overlap -a $my_annot -o ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.txt $j
						  
						  cp ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.txt ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}_toGetLength.txt
						  
						  #take the length of the features with reads (only the parts with reads are taken into account)
						  feat_length=$(grep -v "^#" ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}_toGetLength.txt|grep "${feat_name}_combined"|awk 'BEGIN{a=0}''{if($7>0){a=a+$6}}''END{print a}')
					   
					   fi
					
						#echo -e "feature\tlength" >${output_dir}${feat_name}_${one_condition}_length_rep_${k}.txt
						echo -e "${feat_name}_combined\t${one_condition}\t$feat_length\t${k}" >${output_dir}${feat_name}_${one_condition}_featLength_rep_${k}.txt
						
						grep -v "^#" ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.txt|cut -f1,7|awk 'NR>1{print}' |grep "${feat_name}_combined" >${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.tmp && mv ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.tmp ${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.txt
						
						
					
					
					fi
					
					all_featLength_from_one_cond+=(${output_dir}${feat_name}_${one_condition}_featLength_rep_${k}.txt)
					
					all_rep_one_cond+=(${output_dir}${feat_name}_${one_condition}_counts_rep_${k}.txt)
					
					#increment rep number
					k=$((k+1))
				 
				 done
				 
				 echo -e "feature\tcounts\tcondition\trep" >${output_dir}${feat_name}_${one_condition}_counts_AllReps.txt
				 for one_rep in  $(seq 0 $((${#all_rep_one_cond[*]}-1)));do
				 
					  paste ${all_rep_one_cond[$one_rep]} <(awk -v feat=$feat_name -v condition=${one_condition} -v rep=$((one_rep+1)) 'OFS="\t"{print condition,rep}' ${all_rep_one_cond[$one_rep]}) >>${output_dir}${feat_name}_${one_condition}_counts_AllReps.txt
					
				 
				 done
		 
		 
		 #fi
		 
		 all_counts_all_feat_from_one_cond+=(${output_dir}${feat_name}_${one_condition}_counts_AllReps.txt)
		 
	 
	 
	 done
	 
	 #echo -e "all files for ${one_condition} are : ${all_counts_all_feat_from_one_cond[*]}\n"
	 
	 cat <(cat ${all_counts_all_feat_from_one_cond[*]}|sort -u |grep -P "feature\tcounts") <(cat ${all_counts_all_feat_from_one_cond[*]}|sort -u |grep -v -P "feature\tcounts") >${output_dir}AllFeat_${one_condition}_counts_AllReps.txt
	 
	 cat ${all_featLength_from_one_cond[*]}|sort -u >${output_dir}AllFeat_${one_condition}_featLength_AllReps.txt
	 
	 	 
	 echo -e "feature\tcounts\tcondition\trep" >${output_dir}AllFeat_${one_condition}_perct_AllReps.txt
	 
	 echo -e "feature\tcounts\tcondition\trep" >${output_dir}AllFeat_${one_condition}_norm_perct_AllReps.txt
	 
	
	 ##for each rep of one cond, compute the percentage of counts
	 for one_rep in $(cut -f4 ${output_dir}AllFeat_${one_condition}_counts_AllReps.txt|awk 'OFS="\t"{if(NR>1){print}}'|sort -n -u);do
	 
	 
	   echo -e "\t- rep is : $one_rep\n" 
	 
	    awk -v rep=$one_rep 'OFS="\t"{if($4==rep){print}}' ${output_dir}AllFeat_${one_condition}_counts_AllReps.txt >${output_dir}tmp.tsv
	    
	    #total counts for one rep
	    total_one_rep=$(awk 'BEGIN{a=0}''{a=a+$2}''END{print a}' ${output_dir}tmp.tsv)
	    #set -x
	    #length of all features for one rep
	    #featAllLengthOneRep=$(awk -v rep=$one_rep '{if($4==rep){print $3}}' ${output_dir}AllFeat_${one_condition}_featLength_AllReps.txt|awk 'BEGIN{a=0}''{a=a+$1}''END{print a}')
	      
	    #total norm counts for one rep  
	    #total_normCountsOneRep=$(awk -v total_one_rep=$total_one_rep -v featAllLengthOneRep=$featAllLengthOneRep 'BEGIN{print (total_one_rep*1000)/featAllLengthOneRep}')
	        
	    #compute the total norm counts across all features
	    total_normCountsOneRep=()
	    
	    for one_feat in $(cut -f1 ${output_dir}tmp.tsv|sort -u);do
	    
	      featLengthOneRep=$(awk -v rep=$one_rep -v one_feat=$one_feat '{if($1==one_feat){if($4==rep){print $3}}}' ${output_dir}AllFeat_${one_condition}_featLength_AllReps.txt)
	      
	      one_countNorm=$(grep -P "${one_feat}\t" ${output_dir}tmp.tsv|cut -f2|awk -v total_one_rep=$total_one_rep -v featLengthOneRep=$featLengthOneRep '{print ($1*1000)/featLengthOneRep}')
	      
	      total_normCountsOneRep+=($one_countNorm)
	    
	    
	    done 
	    
	   total_normCountsOneRep=$(echo -e ${total_normCountsOneRep[*]}|sed -e 's/ /\n/g'|awk 'BEGIN{a=0}''{a=a+$1}''END{print a}')
	   
	   echo -e "\t- total counts : $total_one_rep (norm : $total_normCountsOneRep)"
	    
	    
	    for one_feat in $(cut -f1 ${output_dir}tmp.tsv|sort -u);do
	    
	      echo -e "\t\t- $one_feat"
	      
	      #length of one feature for one rep
	      featLengthOneRep=$(awk -v rep=$one_rep -v one_feat=$one_feat '{if($1==one_feat){if($4==rep){print $3}}}' ${output_dir}AllFeat_${one_condition}_featLength_AllReps.txt)
	      	     
	    
	      #counts in perct (raw)
	      one_count_percent=$(grep -P "${one_feat}\t" ${output_dir}tmp.tsv|cut -f2|awk -v total_one_rep=$total_one_rep '{print ($1/total_one_rep)*100}')
	      
	      #counts per kb
	      one_countNorm=$(grep -P "${one_feat}\t" ${output_dir}tmp.tsv|cut -f2|awk -v total_one_rep=$total_one_rep -v featLengthOneRep=$featLengthOneRep '{print ($1*1000)/featLengthOneRep}')
	      
	      #counts per kb
	      one_countNorm_percent=$(awk -v one_countNorm=$one_countNorm -v total_normCountsOneRep=$total_normCountsOneRep 'BEGIN{print (one_countNorm/total_normCountsOneRep)*100}')	      
	      
	      echo -e "\t\t\t- % : $one_count_percent (norm : $one_countNorm_percent)"
	      
	      #exit
	      
	      
	      #exit
	      
	      echo -e "${one_feat}\t${one_count_percent}\t${one_condition}\t${one_rep}" >>${output_dir}AllFeat_${one_condition}_perct_AllReps.txt
	      
	      echo -e "${one_feat}\t${one_countNorm_percent}\t${one_condition}\t${one_rep}" >>${output_dir}AllFeat_${one_condition}_norm_perct_AllReps.txt
	    
	    done
	    
	 
	 
	 done
	 
 
     echo -e "\n******************\n"

done

cond_order=$(echo -e "${all_conditions[*]}"|sed 's/ /,/g')


cat <(find ${output_dir} -name "AllFeat_*_perct_AllReps.txt"|grep -v "_norm_"|xargs cat|grep -P "feature\t"|sort -u) <(find ${output_dir} -name "AllFeat_*_perct_AllReps.txt"|grep -v "_norm_"|xargs cat|grep -v -P "feature\t") >${output_dir}AllFeat_AllConds_AllReps_AllPerct.tsv

$genomic_distrib_stats ${output_dir}AllFeat_AllConds_AllReps_AllPerct.tsv ${output_dir} $cond_order

mv ${output_dir}reads_genomic_distrib_model1.png ${output_dir}reads_genomic_distrib_raw.png

mv ${output_dir}reads_genomic_distrib_model1_NoLabel.png ${output_dir}reads_genomic_distrib_raw_NoLabel.png


cat <(find ${output_dir} -name "AllFeat_*_norm_perct_AllReps.txt"|xargs cat|grep -P "feature\t"|sort -u) <(find ${output_dir} -name "AllFeat_*_norm_perct_AllReps.txt"|xargs cat|grep -v -P "feature\t") >${output_dir}AllFeat_AllConds_norm_AllReps_AllPerct.tsv

$genomic_distrib_stats ${output_dir}AllFeat_AllConds_norm_AllReps_AllPerct.tsv ${output_dir} $cond_order

mv ${output_dir}reads_genomic_distrib_model1.png ${output_dir}reads_genomic_distrib_normalized.png

mv ${output_dir}reads_genomic_distrib_model1_NoLabel.png ${output_dir}reads_genomic_distrib_normalized_NoLabel.png





