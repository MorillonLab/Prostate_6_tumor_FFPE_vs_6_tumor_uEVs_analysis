#!/bin/bash

#rmk on the website :Your job 601D9FC5000039BB56455ADF have been killed
#This may be due to excessive time consumption of 36002 seconds. The time limit for this service is 36000 seconds.

netMHCpan="/home/marcgabriel/Downloads/netMHCpan/netMHCpan-4.1/netMHCpan"

#output_dir="/home/marcgabriel/Downloads/netMHCpan/netMHCpan-4.1/output_dir/"

#input_peptides="/home/marcgabriel/Downloads/netMHCpan/netMHCpan-4.1/test/B0702.fsa"
#input_peptides="/home/marcgabriel/Downloads/netMHCpan/netMHCpan-4.1/peptides_test.fa"

#input_peptides="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/Expressed_Transcripts_length/peptides_from_up_FFPE.fa"
#input_peptides="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/Expressed_Transcripts_length/peptides_from_up_urinaryEVs.fa"
#input_peptides="/media/marcgabriel/saylar5/neoantigens_test/supplementary_table_7_refined.fa"

#input_peptides="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/intersection_of_db/Up_transcripts_uEVs_NotIn_PCG_pseudogenes_highestExpressedTranscript_peptides_refined.fa"

#input_peptides="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/intersection_of_db/Up_transcripts_FFPE_NotIn_PCG_pseudogenes_highestExpressedTranscript_peptides_refined.fa"

#input_peptides="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/intersection_of_db/Up_transcripts_uEVs_NotIn_PCG_pseudogenes_highestExpressedTranscript_peptides_three_frames_peptides_refined.fa"

#input_peptides="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/intersection_of_db/Up_transcripts_FFPE_NotIn_PCG_pseudogenes_highestExpressedTranscript_peptides_three_frames_peptides_refined.fa"

#input_peptides="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_riboseq_DIS3_diagenode_protocol_16_12_2021/NETMHCpan_analysis/peptides_frame_1_from_1519_candidates_refined_for_netmhcpan.fa"
#input_peptides="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_riboseq_DIS3_diagenode_merged_runs_15_11_2021_and_16_12_2021/Dominika_ribotricer_analysis_gencode_scallop_eRNA_selected_length_from_unique_25_34nt/ORF_diff_final_in_DIS3_peptides_refined.fa"


#finally, we never used this result
#input_peptides="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_riboseq_DIS3_diagenode_merged_runs_15_11_2021_and_16_12_2021/Dominika_ribotricer_analysis_MM_enhancers/MM_enhancers_peptides_analysis/MM_enhancers_peptides_refined.fa"

input_peptides="/home/marcgabriel/Downloads/ORF_diff_translated_DIS3_A_peptides_refined.fa"



all_HLA="/home/marcgabriel/Downloads/netMHCpan/netMHCpan-4.1/Linux_x86_64/data/allelenames"

#HLA_from_RNAseq="/media/marcgabriel/saylar5/seq2HLA_results_FFPE_UrinaryEVs/all_significant_expressed_ClassI_HLA_FFPE.tsv /media/marcgabriel/saylar5/seq2HLA_results_FFPE_UrinaryEVs/all_significant_expressed_ClassI_HLA_EVs.tsv"
#HLA_from_RNAseq="/media/marcgabriel/saylar5/neoantigens_test/CD8_epitopes.tsv"

HLA_from_RNAseq="/home/marcgabriel/Downloads/HLA_data/combined_HLA_types.tsv"

#output_dir="/media/marcgabriel/saylar5/neoantigens_results_FFPE_UrinaryEVs/"
#output_dir="/media/marcgabriel/saylar5/neoantigens_test/output_results/"

#output_dir="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/intersection_of_db/neoantigens_1st_frame_uEVs/"

#output_dir="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/intersection_of_db/neoantigens_1st_frame_FFPE/"

#output_dir="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/intersection_of_db/neoantigens_three_frames_uEVs2/"

#output_dir="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/intersection_of_db/neoantigens_three_frames_FFPE/"

#output_dir="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_riboseq_DIS3_diagenode_protocol_16_12_2021/NETMHCpan_analysis/NETMHCpan_results/"

#output_dir="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_riboseq_DIS3_diagenode_merged_runs_15_11_2021_and_16_12_2021/Dominika_ribotricer_analysis_gencode_scallop_eRNA_selected_length_from_unique_25_34nt/NETMHCpan_classI_analysis/"

#output_dir="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/NETMHCpan_classI_analysis/"

#output_dir="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/NETMHCpan_classI_analysis_MM_enhancers/"

#this one gives tmpdir too long according to the tool
#output_dir="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_riboseq_DIS3_diagenode_merged_runs_15_11_2021_and_16_12_2021/Dominika_ribotricer_analysis_gencode_scallop_eRNA_selected_length_from_unique_25_34nt_extendedStartCodon/NETMHCpan_classI_analysis/"

output_dir="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_extendedStartCodon_NETMHCpan_classI_analysis/"

select_human_antigen="yes"

process_number_limit=30

#create output dir
if [ ! -d $output_dir ];then mkdir -p $output_dir;fi

#create a tmp dir for HLA files (split in many chunks to be parallelized)
if [ -d ${output_dir}HLA_list/ ];then rm -rf ${output_dir}HLA_list/;fi
mkdir -p ${output_dir}HLA_list/

#create a tmp dir for peptides files (split in many chunks to be parallelized)
tmp_for_peptides="${output_dir}peptides_tmp/"
if [[ -d $tmp_for_peptides ]];then rm -rf $tmp_for_peptides;fi
mkdir $tmp_for_peptides

sample_scripts_dir="${output_dir}/sample_scripts/"
if [ -d $sample_scripts_dir ]; then rm -rf $sample_scripts_dir; fi 
mkdir $sample_scripts_dir

>${output_dir}netMHCpan_global_summary.txt


#select only human antigens
if [[ $select_human_antigen == "yes" ]];then

	grep -i "^HLA" $all_HLA >${output_dir}human_antigenes.tsv
	
	selected_antigens=${output_dir}human_antigenes.tsv
	
	echo -e "human antigens are selected : $selected_antigens!\n"


else

	selected_antigens=$all_HLA
	
	echo -e "all antigens are selected (it will be very long) : $all_HLA !\n"

fi


selected_antigens_from_RNAseq=${output_dir}selected_antigens_from_RNAseq.tsv

cut -f1 $HLA_from_RNAseq | sed s/"'"//g |sort -u >${output_dir}HLA_from_RNAseq.tsv

grep -F -f ${output_dir}HLA_from_RNAseq.tsv $selected_antigens |cut -f1 >$selected_antigens_from_RNAseq

echo -e "selected antigens from RNAseq are there : $selected_antigens_from_RNAseq!\n"

#antigens list in chunks of 15 lines max
awk -v split_lines=15 -v HLA_list=${output_dir}HLA_list/ 'NR%split_lines==1{OFS="\t";x=++i"_subfile.txt"}{OFS="";print $1 > HLA_list x}' $selected_antigens_from_RNAseq


#sort numerically the chunks
split_HLA=($(cd ${output_dir}HLA_list/ ; ls |grep -E "[0-9]+_subfile.txt"|sort -n))

#check the pos of the last chunk in the list
last_chunk_pos=$((${#split_HLA[*]}-1))

#take the last chunk with its pos
last_chunk_file=${output_dir}HLA_list/${split_HLA[$last_chunk_pos]}



#give the number of lines of the last chunk
if [[ ${#split_HLA[*]} -gt 1 ]];then

	echo -e "-> ${#split_HLA[*]} chunks of HLA ; each of them has 50 lines, and the last one has $(wc -l $last_chunk_file|awk '{print $1}') line(s) (file $last_chunk_file)\n"
		
else

   echo -e "1 chunk of HLA with $(wc -l $last_chunk_file|awk '{print $1}') line(s)"

fi




#each line is a fasta seq of peptide
awk 'BEGIN{RS=">"}''{if(NR>1){sub("\n","\t"); gsub("\n",""); print RS$0}}' $input_peptides >${tmp_for_peptides}input_peptides.tsv

awk -v split_lines=200 -v peptide_dir=${tmp_for_peptides} 'NR%split_lines==1{OFS="\t";x=++i"_subfile.txt"}{OFS="\t";print $1,$2 > peptide_dir x}' ${tmp_for_peptides}input_peptides.tsv

#sort numerically the chunks
split_peptides=($(cd ${tmp_for_peptides} ; ls |grep -E "[0-9]+_subfile.txt"|sort -n))

#check the pos of the last chunk in the list
last_chunk_pos=$((${#split_peptides[*]}-1))

#take the last chunk with its pos
last_chunk_file=${tmp_for_peptides}${split_peptides[$last_chunk_pos]}


#give the number of lines of the last chunk
if [[ ${#split_peptides[*]} -gt 1 ]];then

	echo -e "-> ${#split_peptides[*]} chunks of peptides ; each of them has 200 lines, and the last one has $(wc -l $last_chunk_file|awk '{print $1}') line(s) (file $last_chunk_file)\n"
		
else

   echo -e "1 chunk of HLA with $(wc -l $last_chunk_file|awk '{print $1}') line(s)"

fi


all_chunks_HLA=()

#loop across the pos of each chunk of peptides
for one_chunk in $(seq 0 $((${#split_peptides[*]}-1)));do


   #select one chunk using its pos
   one_peptide=${split_peptides[$one_chunk]}
   
   
   	#create its dir
	#if [[ -d "${output_dir}${one_chunk}_chunk_peptide_dir/" ]];then rm -rf ${output_dir}${one_chunk}_chunk_peptide_dir/;fi
	
	if [[ ! -d "${output_dir}${one_chunk}_chunk_peptide_dir/" ]];then mkdir ${output_dir}${one_chunk}_chunk_peptide_dir/;fi
	
	
	awk -v OFS="\n" '{print $1,$2}' ${tmp_for_peptides}$one_peptide|grep -v "^$" >${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk}_chunk_peptide.tsv
	

  #loop across the pos of each chunk of HLA
  for one_chunk2 in $(seq 0 $((${#split_HLA[*]}-1)));do
  
            #remove existing files
           	#if [[ -d "${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/" ]];then rm -rf ${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/;fi
           	#mkdir ${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/
           	
           	
           	#create its dir
           	if [[ ! -d "${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/" ]];then mkdir ${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/;fi
           
	        

		   #select one chunk using its pos
		   one_HLA=${split_HLA[$one_chunk2]}


			#EL_Rank : <= 0.5 ->  strong binders
			#EL_Rank : 0.5 > X <= 2 ->  weak binders
			
			#if the xls file doesn't run the netmhcpan command
			if [[ ! -f ${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/Peptide_MHC_Class_I_Interaction_Predictions_chunk_peptide_${one_chunk}_chunk_HLA_${one_chunk2}.xls ]];then
			
					#put the commands common for all the chunks in a subscript
					echo -e "#!/bin/bash\n\n" >${sample_scripts_dir}subscript_peptide_${one_chunk}_HLA_${one_chunk2}.sh
					echo -e "start_date=\$(date)\n" >>${sample_scripts_dir}subscript_peptide_${one_chunk}_HLA_${one_chunk2}.sh
					
						#determine the binding of the supplied peptides using this chunk
						echo -e "$netMHCpan -tdir ${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/netMHCpanXXXXXX -l 8,9,10,11,12,13,14 -t 2.1 -a \$(cat ${output_dir}HLA_list/${one_HLA}|tr '\\\n' ','|sed 's/,\$//g'|awk '{print \$0}') -xls -xlsfile ${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/Peptide_MHC_Class_I_Interaction_Predictions_chunk_peptide_${one_chunk}_chunk_HLA_${one_chunk2}.xls ${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk}_chunk_peptide.tsv 1>${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/netMHCpan_summary_peptide_chunk_${one_chunk}_HLA_chunk_${one_chunk2}.txt\n" >>${sample_scripts_dir}subscript_peptide_${one_chunk}_HLA_${one_chunk2}.sh
					
					echo -e "end_date=\$(date)\n" >>${sample_scripts_dir}subscript_peptide_${one_chunk}_HLA_${one_chunk2}.sh
					
					echo -e "echo -e \"start : \$start_date ; end : \$end_date for ${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/netMHCpan_summary_peptide_chunk_${one_chunk}_HLA_chunk_${one_chunk2}.txt\\\n\"\n" >>${sample_scripts_dir}subscript_peptide_${one_chunk}_HLA_${one_chunk2}.sh
					
					chmod 755 ${sample_scripts_dir}subscript_peptide_${one_chunk}_HLA_${one_chunk2}.sh
			
			
			
			fi
			
			all_chunks_HLA+=(${output_dir}${one_chunk}_chunk_peptide_dir/${one_chunk2}_HLA_chunk/netMHCpan_summary_peptide_chunk_${one_chunk}_HLA_chunk_${one_chunk2}.txt)
			
			
			
			#command to check whether a sub-file is incrementing
			#this will give the number of HLA already processed
			#grep "using nearest neighbor HLA" /media/marcgabriel/saylar5/neoantigens_results_FFPE_UrinaryEVs/0_chunk_peptide_dir/0_HLA_chunk/netMHCpan_summary_peptide_chunk_0_HLA_chunk_0.txt
			
			#this will give the number of seq that is being processed
			#grep "Protein" /media/marcgabriel/saylar5/neoantigens_results_FFPE_UrinaryEVs/0_chunk_peptide_dir/0_HLA_chunk/netMHCpan_summary_peptide_chunk_0_HLA_chunk_0.txt|cut -d' ' -f2|sort -u
			
			#a way to visualize the results : https://dmnfarrell.github.io/bioinformatics/using-epitopepredict
			#plot where the x axis is the HLA proteine seq, the rectangles are the neantigens (epitopes), and the y-axis are the different allele names of the HLA
			
			
  done
	
	

done

#run the subscripts in parallel
find ${sample_scripts_dir} -name "*subscript_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "bash subscripts failure !" 1>&2; exit; }

#exit

#concatenate the results in a file
#find ${output_dir} -name "*netMHCpan_summary_part_" | xargs cat  >${output_dir}netMHCpan_global_summary.txt

>${output_dir}netMHCpan_global_summary.txt
for F in ${all_chunks_HLA[*]};do cat $F >>${output_dir}netMHCpan_global_summary.txt;done
#cat ${all_chunks_HLA[*]} >${output_dir}netMHCpan_global_summary.txt



#process the output result
#to know the meaning of the columns : http://www.cbs.dtu.dk/services/NetMHCpan-4.1/output.php
cat <(grep "\sIdentity\s" ${output_dir}netMHCpan_global_summary.txt|sort -u|sed -E 's/^ //g'|sed -E 's/\s+/\t/g') <(grep -v "#" ${output_dir}netMHCpan_global_summary.txt|grep -E "^\s*[0-9]+"|sed -E 's/^ //g'|sed -E 's/\s+/\t/g'|awk 'OFS="\t"{if($13<=0.5){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,"strong"}else if ($13<=2){print  $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,"weak"}}') >${output_dir}netMHCpan_final_results.tsv

echo -e "\n-> check final results : ${output_dir}netMHCpan_final_results.tsv\n"

#exit

#remove the chunks dir
#find ${output_dir} | grep -E "[0-9]+_chunk_peptide_dir.*" |while read line;do rm -rf $line;done




