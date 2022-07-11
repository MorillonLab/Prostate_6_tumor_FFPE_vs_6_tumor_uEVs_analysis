#!/bin/bash

#this script thanks to a design file, will give for each sample of each condition, the counts with SEQ2HLA.
#the results will be in the output directory that you have set
#each sample has its own sub-directory, and in each, the counts are in *_abundance.tsv


############ inputs #############
#################################

#the design should be something like this (2 columns, tab-separated) :
# sample_name	condition
# A	cond1
# B	cond1
# O	cond2
# P	cond2
# ...


#design="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode_and_repeats/design_FFPE.tsv"
#design="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode_and_repeats/design_tumor.tsv"

design="/media/marcgabriel/saylar9/Dominika_HCT116_DIS3_HLA_type/seq2HLA_design.txt"


#For FFPE files
#reads_path="/media/marcgabriel/saylar4/urines_ffpe_files/fastq_files/"

#For urine files
#reads_path="/media/marcgabriel/saylar4/urines_files"

reads_path="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/"



#output files
#counter_outputs="/media/marcgabriel/saylar5/seq2HLA_results_FFPE_UrinaryEVs"
counter_outputs="/media/marcgabriel/saylar9/Dominika_HCT116_DIS3_HLA_type/"


#e.g : pre_R_1_, post_PD_1_
#separator_after_file_name="_"
separator_after_file_name="_"



#number of process (samples) to run in parallel
process_number_limit=4

#threads used by SEQ2HLA for one sample
SEQ2HLA_threads=2

#SEQ2HLA program
SEQ2HLA="/home/marcgabriel/Downloads/seq2HLA-master/seq2HLA.py"




###############  process the inputs ###################
#######################################################

#standardize the output dir name
counter_outputs="${counter_outputs}/"
counter_outputs=$(echo $counter_outputs |sed 's/\/\//\//g')
if [ ! -d $counter_outputs ]; then mkdir -p $counter_outputs; fi 

#directory of scripts for each sample (to be laucnhed in parallel)
sample_scripts_dir="${counter_outputs}/sample_scripts/"
if [ -d $sample_scripts_dir ]; then rm -rf $sample_scripts_dir; fi 
mkdir $sample_scripts_dir


#look for all conds
all_conds=($(grep -v "^#" $design|cut -f2 |sort -u))

echo -e "all conds are : ${all_conds[*]}\n"

#array to store all files
all_files=()

#array with all new names
all_new_names=()

#set -x
#for each sample of each cond, assign a name based on the condition (for conds normal & tumoral with 2 samples each : normal_1, normal_2...tumoral_1, tumoral_2)
for one_cond in $(seq 0 $((${#all_conds[*]}-1)));do

	#for a given cond, take all its files (samples)
	#files_one_cond=($(grep -P "\t${all_conds[$one_cond]}$" $design |cut -f1|sort -n))
	files_one_cond=($(grep -v "^#" $design |grep -P "\t${all_conds[$one_cond]}$" |cut -f1))
 
 
	#number of the 1st rep of the condition
	cond_start=1
 
	#number of the last rep of the condition
	cond_end=${#files_one_cond[*]}
 
	#store all the new names for this condition
	cond_new_names=()
 
    #loop across the samples of the cond, the for each of them concatenate the name of the condition & the number
	for ((i=$cond_start; i<=$cond_end; i++)); do cond_new_names+=(${all_conds[$one_cond]}_${i}) ; done
 
	all_files+=(${files_one_cond[*]})
 
	all_new_names+=(${cond_new_names[*]})

done

#echo -e "all files are : ${all_files[*]}"

#echo -e "new names are : ${all_new_names[*]}"

echo -e "new design : "

paste <(echo -e ${all_files[*]}|sed 's/ /\n/g') <(echo -e ${all_new_names[*]} |sed 's/ /\n/g')

paste <(echo -e ${all_files[*]}|sed 's/ /\n/g') <(echo -e ${all_new_names[*]} |sed 's/ /\n/g') >${counter_outputs}new_design.tsv




#summary of the count files
>${counter_outputs}summary.tsv


all_sig_classI=()

all_sig_classII=()


##############  run a loop across the samples in order to process them #################
########################################################################################

#create a script for each sample, that will run Kallisto & process the output file
for one_sample in $(seq 0 $((${#all_files[*]}-1)));do

          #for a given sample of a given cond, take its new name
		  one_sample_name=$(echo "${all_new_names[$one_sample]}")
		  
		  #path to the sample
		  counter_outputs_OneSample="${counter_outputs}${one_sample_name}/"
		  
		  if [ ! -d $counter_outputs_OneSample ]; then mkdir $counter_outputs_OneSample; fi 
		  
		   echo -e "#!/bin/bash\n\n" >${sample_scripts_dir}subscript_${one_sample_name}.sh
		   
		 
			  
		   echo -e "\nif [ ! -f ${counter_outputs_OneSample}${one_sample_name}-ClassI-class.HLAgenotype4digits ] || [ ! -f ${counter_outputs_OneSample}${one_sample_name}-ClassI-nonclass.HLAgenotype4digits ];then\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
		   
		   #set -x
		  		 
			  fastq_list_OneSample=($(find $reads_path -name "${all_files[$one_sample]}${separator_after_file_name}*" |grep ".*\.fastq"|sort -n)) || { echo "no fastq files 1 !!" 1>&2; exit; }
			  
			  #set +x
			  
			  echo -e "echo -e \"fastq files for condition $one_sample_name : \n ${fastq_list_OneSample[*]} \n--------------\n\"" >>${sample_scripts_dir}subscript_${one_sample_name}.sh

			  
			  if [[ ${#fastq_list_OneSample[*]} -eq 2 ]];then
			  
			 
			  
				  fastq_list_OneSample_R1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "1\.fastq")
				 
				  fastq_list_OneSample_R2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "2\.fastq")
				  
				 
			      
				  #run the quantif
			      echo -e "python ${SEQ2HLA} -1 ${fastq_list_OneSample_R1} -2 ${fastq_list_OneSample_R2} -p $SEQ2HLA_threads -r ${counter_outputs_OneSample}${one_sample_name}  || { echo \"seq2HLA failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			      
				  #check if sed \"s/'//g\" works perfectly well there 
				  echo -e "sed \"s/'//g\" ${counter_outputs_OneSample}${one_sample_name}-ClassI-class.HLAgenotype4digits |awk 'OFS=\"\\\t\"{if(NR>1){for(i=3;i<=NF;i++){if(i==3){if(\$i<=0.05){print \$2,\$3}};if(i==5){if(\$5<=0.05){print \$4,\$5}}}}}' |sort -k1,1 -k2,2g |sort -u -k1,1 -m >${counter_outputs_OneSample}${one_sample_name}_significant_expressed_ClassI_HLA.tsv\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
				  
				  echo -e "sed \"s/'//g\" ${counter_outputs_OneSample}${one_sample_name}-ClassI-nonclass.HLAgenotype4digits | awk 'OFS=\"\\\t\"{if(NR>1){for(i=3;i<=NF;i++){if(i==3){if(\$i<=0.05){print \$2,\$3}};if(i==5){if(\$5<=0.05){print \$4,\$5}}}}}' |sort -k1,1 -k2,2g |sort -u -k1,1 -m >>${counter_outputs_OneSample}${one_sample_name}_significant_expressed_ClassI_HLA.tsv\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
				  
				  
				  echo -e "sed \"s/'//g\" ${counter_outputs_OneSample}${one_sample_name}-ClassII.HLAgenotype4digits | awk 'OFS=\"\\\t\"{if(NR>1){for(i=3;i<=NF;i++){if(i==3){if(\$i<=0.05){print \$2,\$3}};if(i==5){if(\$5<=0.05){print \$4,\$5}}}}}' |sort -k1,1 -k2,2g |sort -u -k1,1 -m >${counter_outputs_OneSample}${one_sample_name}_significant_expressed_ClassII_HLA.tsv\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			      
			      
			      
			  
			  else
			  
			  
			    echo -e "issues with fastq files for ${one_sample_name} (you should have paired-end reads) ! \n"
			  
			  
			    exit
			  
			  
			  fi
			  
			  
			  echo -e "\nelse\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  #check if sed \"s/'//g\" works perfectly well there 
			  echo -e "sed \"s/'//g\" ${counter_outputs_OneSample}${one_sample_name}-ClassI-class.HLAgenotype4digits |awk 'OFS=\"\\\t\"{if(NR>1){for(i=3;i<=NF;i++){if(i==3){if(\$i<=0.05){print \$2,\$3}};if(i==5){if(\$5<=0.05){print \$4,\$5}}}}}' |sort -k1,1 -k2,2g |sort -u -k1,1 -m >${counter_outputs_OneSample}${one_sample_name}_significant_expressed_ClassI_HLA.tsv\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  echo -e "sed \"s/'//g\" ${counter_outputs_OneSample}${one_sample_name}-ClassI-nonclass.HLAgenotype4digits | awk 'OFS=\"\\\t\"{if(NR>1){for(i=3;i<=NF;i++){if(i==3){if(\$i<=0.05){print \$2,\$3}};if(i==5){if(\$5<=0.05){print \$4,\$5}}}}}' |sort -k1,1 -k2,2g |sort -u -k1,1 -m >>${counter_outputs_OneSample}${one_sample_name}_significant_expressed_ClassI_HLA.tsv\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  
			  echo -e "sed \"s/'//g\" ${counter_outputs_OneSample}${one_sample_name}-ClassII.HLAgenotype4digits | awk 'OFS=\"\\\t\"{if(NR>1){for(i=3;i<=NF;i++){if(i==3){if(\$i<=0.05){print \$2,\$3}};if(i==5){if(\$5<=0.05){print \$4,\$5}}}}}' |sort -k1,1 -k2,2g |sort -u -k1,1 -m >${counter_outputs_OneSample}${one_sample_name}_significant_expressed_ClassII_HLA.tsv\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  
			  all_sig_classI+=("${counter_outputs_OneSample}${one_sample_name}_significant_expressed_ClassI_HLA.tsv")

			  all_sig_classII+=("${counter_outputs_OneSample}${one_sample_name}_significant_expressed_ClassII_HLA.tsv")
			  			  
			  
			  echo -e "\nfi\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  

		
			  #give rights to the script	
			  chmod 755 ${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			
		
					   
done


#run the subscripts in parallel

#here, if you are on a cluster, and you want to use a qsub command, just use a loop on the content of ${sample_scripts_dir}, and give the result to qsub
find ${sample_scripts_dir} -name "*subscript_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "executing bash sorting and/or counting subsritpts failure !" 1>&2; exit; }

#take the alleles (no redundancy : best according to the confidence/p-value)
cat ${all_sig_classI[*]} |sort -k1,1 -k2,2g |sort -u -k1,1 -m >${counter_outputs}all_significant_expressed_ClassI_HLA.tsv
cat ${all_sig_classII[*]}  |sort -k1,1 -k2,2g |sort -u -k1,1 -m >${counter_outputs}all_significant_expressed_ClassII_HLA.tsv


echo -e "\n-> check ${counter_outputs}all_significant_expressed_ClassI_HLA.tsv & ${counter_outputs}all_significant_expressed_ClassII_HLA.tsv\n"









