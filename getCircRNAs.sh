#!/bin/bash

reads_path="/media/marcgabriel/saylar4/urines_files/ /media/marcgabriel/saylar4/urines_ffpe_files/fastq_files/"

config_file="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/circant_results/circant_config.yml"

output_dir="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/ciriquant_results/circant_results_out/"


my_gff="/home/marcgabriel/Documents/gencode32/gencode.v32.annotation_sorted.gff3"

CIRIquant="CIRIquant"

bedtools="bedtools"

gffread="/home/marcgabriel/Downloads/cufflinks-2.2.1.Linux_x86_64/gffread"

#1 : fr-secondstrand ; 2 : fr-firststrand ; 0 : unstranded
#lib_type=2

lib_type=2

process_number_limit=4

thread_num=4

design_list="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/circant_results/design_FFPE_urines.tsv"


#separator="_trimmed"
#separator="."
separator="_"


split_design=($(echo -e $design_list| sed 's/,/\n/g'))


output_dir_split=($(echo -e $output_dir_list| sed 's/,/\n/g'))

getConsensus_seq="/home/marcgabriel/Desktop/scripts/getConsensus.R"

#set -x
for one_design in $(seq 0 $((${#split_design[*]}-1))) ;do


		design=${split_design[$one_design]}
		
		echo -e "project : $design...\n"


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



		############ processing of variables #########
		#############################################
		
		
		
		output_dir=${output_dir_split[$one_design]}
		

		output_dir="${output_dir}/"
		output_dir=$(echo $output_dir |sed 's/\/\//\//g')
		if [ ! -d $output_dir ]; then mkdir -p $output_dir; fi 
		
	

		sample_scripts_dir="${output_dir}/sample_scripts/"
		if [ -d $sample_scripts_dir ]; then rm -rf $sample_scripts_dir; fi 
		mkdir $sample_scripts_dir
		
		
		echo -e "new design : "

		paste <(echo -e ${all_files[*]}|sed 's/ /\n/g') <(echo -e ${all_new_names[*]} |sed 's/ /\n/g')

		paste <(echo -e ${all_files[*]}|sed 's/ /\n/g') <(echo -e ${all_new_names[*]} |sed 's/ /\n/g') >${output_dir}new_design.tsv
		
	    #array to store the files containing junction seq
		all_consensus_junction_sequences=()


		########### process gff in order to have an annotation for circview : ######
		############################################################################

		#the construction of this annotation is based on the example of annotations supplied with circView

		#select transcript_IDs, keep chr, strand, gene_name, start (0-based), end
		#group exons by transcript_IDs, keep only the starts (0-based),ends, gene_type
		#rearrange the columns that are not (gene_name should be before transcript_id)
		LANG=en_EN join -t $'\t' -a1 -11 -21 <(grep -P "\ttranscript\t" $my_gff|awk 'OFS="\t"{print $9,$1,$7,$4-1,$5,$5,$5}'|awk 'OFS="\t"{split($1,a,";");for(i=1;i<=length(a);i++){if(a[i]~/transcript_id/){tp_id=a[i]};if(a[i]~/gene_name/){gene_name=a[i]}};print tp_id,gene_name,$2,$3,$4,$5,$6,$7}'|sed 's/transcript_id=//g'|sed 's/gene_name=//g'|LANG=en_EN sort -k1,1) <(grep -P "\texon\t" $my_gff|awk 'OFS="\t"{print $9,$4-1,$5}'|awk 'OFS="\t"{split($1,a,";");for(i=1;i<=length(a);i++){if(a[i]~/transcript_id/){tp_id=a[i]};if(a[i]~/gene_type/){gene_type=a[i]}};print tp_id,$2,$3,gene_type}'|sed 's/transcript_id=//g'|sed 's/gene_type=//g'|sort -k1,1 -k2,2n|$bedtools groupby -g 1 -c 2,2,3,4 -o count,collapse,collapse,distinct|awk 'OFS="\t"{print $1,$2,$3",",$4",",$5}'|LANG=en_EN sort -k1,1) |awk 'OFS="\t"{col_1=$2;col_2=$1;$1=col_1;$2=col_2;print}' >${output_dir}circview_annotation.txt


		############ analysis #######################
		#############################################

		ciriquant_design=${output_dir}circant_design.tsv

		#design of CIRIquant that will be used to summarize the results across all conditions
		>${ciriquant_design}


		#iteration across samples
		for one_sample in $(seq 0 $((${#all_files[*]}-1)));do



		  #sample
		  #one_sample_name=$(echo "${motif_list[$one_sample]}")
		  
		  #for a given sample of a given cond, take its new name
		  one_sample_name=$(echo "${all_new_names[$one_sample]}")  
		  
		  #retrieve the condition
		  one_cond=$(grep -P "^${all_files[$one_sample]}\t" $design |cut -f2)
		  
		  
		  #path to the sample
		  tophat_outputs_OneSample="${output_dir}${one_sample_name}/"
		  
		  if [ ! -d $tophat_outputs_OneSample ]; then mkdir $tophat_outputs_OneSample; fi 
		  
		  
			  fastq_list_OneSample=($(find ${reads_path} -name "${all_files[$one_sample]}${separator}*" |grep ".*\.fastq"|sort -u|head -n2)) || { echo "no fastq files 1 !!" 1>&2; exit; }
			  
			  
			  if [[ ${#fastq_list_OneSample[*]} -eq 0 ]];then
			  
			    echo -e "\nno files for ${all_files[$one_sample]}${separator} !\n\n"
			    
			    exit 1
			    
			  
			  
			  fi 
			  
			  echo -e "#!/bin/bash\n\n" >${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  #if the last file produced by ciriquant is not there, run ciriquant
			  #remove first all the temp files (could be very heavy, like 100+ G !!)
			  #[ ! -f ${tophat_outputs_OneSample}/circ/${one_sample_name}_index*.sa ]
			  echo -e "\nif [ ! -f ${tophat_outputs_OneSample}circ/${one_sample_name}_denovo.sorted.bam.bai ] ;then\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  			  
			  echo -e "\n\tcd ${tophat_outputs_OneSample} && rm -rf ${tophat_outputs_OneSample}* \n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh


			  echo -e "echo -e \"fastq files for condition $one_sample_name : \n ${fastq_list_OneSample[*]} \n--------------\n\"" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			 

			  
			  if [[ ${#fastq_list_OneSample[*]} -eq 2 ]];then
			  
				fastq_1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep -E "_1\.fastq")
				
				fastq_2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep -E "_2\.fastq")
				
				if [[ "$fastq_1" == "" ]];then
				
					fastq_1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep -i -E "\.R1\.fastq")
					
					fastq_2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep -i -E "\.R2\.fastq")				
				
				
				
				fi
					 
					 #if [[ "${just_counting}" == "yes" ]] && [[ -f "$known_circ_bed" ]];then
					 
							#echo -e "$CIRIquant -v -l $lib_type -t $thread_num -1 $fastq_1 -2 $fastq_2 --circ $known_circ_bed --config $config_file -o $tophat_outputs_OneSample -p $one_sample_name\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
							
							#echo -e "-> we're just going to do counting for $one_sample_name !!\n\n"
					 
					 
					 #else
					 
					 
					    echo -e "$CIRIquant -v -l $lib_type -t $thread_num -1 $fastq_1 -2 $fastq_2 --config $config_file -o $tophat_outputs_OneSample -p $one_sample_name\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
					 
					 
					 #fi
					 
			  
			  fi
			  
	  
							 
									
			  echo -e "\nfi\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  chmod 755 ${sample_scripts_dir}subscript_${one_sample_name}.sh
			  

			   
			   
		done
		
		#run in parallel the sub-scripts
		find ${sample_scripts_dir} -name "*subscript_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "executing bash sorting and/or counting subsritpts failure !" 1>&2; exit; }
		
		
		
		## dev ##

		#iteration across samples
		#pro
		for one_sample in $(seq 0 $((${#all_files[*]}-1)));do



		  #sample
		  #one_sample_name=$(echo "${motif_list[$one_sample]}")
		  
		  #for a given sample of a given cond, take its new name
		  one_sample_name=$(echo "${all_new_names[$one_sample]}")  
		  
		  #retrieve the condition
		  one_cond=$(grep -P "^${all_files[$one_sample]}\t" $design |cut -f2)
		  
		  
		  #path to the sample
		  tophat_outputs_OneSample="${output_dir}${one_sample_name}/"

			   if [ -f ${tophat_outputs_OneSample}/${one_sample_name}.gtf ];then
			   
			   
			   
					awk -F'\t' 'OFS="\t"{if($1!~/^#/){$9=gensub(":","_","g",$9);$9=gensub("\\|","_","g",$9);$9=gensub("circ_id \"","circ_id \"circRNA_","g",$9);print}else{print}}' ${tophat_outputs_OneSample}/${one_sample_name}.gtf >${tophat_outputs_OneSample}/${one_sample_name}_refined.gtf
					 
				   echo -e "${one_sample_name}\t${tophat_outputs_OneSample}/${one_sample_name}_refined.gtf\t${one_cond}" >>${ciriquant_design}
				  
				  
				  #process ciri output, apparently, there's a strange character at the end of the file that is not recognized by java
				  #nb_lines=$(wc -l ${tophat_outputs_OneSample}/circ/${one_sample_name}.ciri|awk '{print $1}')
				  
				  
				  #head -n $nb_lines ${tophat_outputs_OneSample}/circ/${one_sample_name}.ciri >${tophat_outputs_OneSample}/circ/${one_sample_name}.ciri.tmp && mv ${tophat_outputs_OneSample}/circ/${one_sample_name}.ciri.tmp ${tophat_outputs_OneSample}/circ/${one_sample_name}.ciri
				  
				  
				  tmp_dir="${tophat_outputs_OneSample}/circ/${one_sample_name}_tmp_dir/"
							  
			      if [ ! -d "${tmp_dir}" ];then mkdir "${tmp_dir}" ;fi
				  
				  
				  if [[ ! -f ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_5_3_orientation.fa ]] || [[ $(wc -l ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_5_3_orientation.fa|awk '{print $1}') -eq 0 ]];then
				  
							   if [ -d "${tophat_outputs_OneSample}/circ/${one_sample_name}_getConsensus_subscripts_dir/" ];then rm -rf "${tophat_outputs_OneSample}/circ/${one_sample_name}_getConsensus_subscripts_dir/";fi
							   
							   mkdir "${tophat_outputs_OneSample}/circ/${one_sample_name}_getConsensus_subscripts_dir/"         
							  

							   
							   ### to do : put all this par in a function, then export it in all subscripts
							  
							  
							   #select all junction reads
							   cut -f12 ${tophat_outputs_OneSample}/circ/${one_sample_name}.ciri|tr ',' '\n'|grep -v "^$"|sort -u|grep -v "junction_reads" >${tmp_dir}${one_sample_name}_junction_read_IDs.txt
							   
							   #convert them in bam
							   cat <(samtools view -H ${tophat_outputs_OneSample}/circ/${one_sample_name}_denovo.sorted.bam) <(samtools view -F4 ${tophat_outputs_OneSample}/circ/${one_sample_name}_denovo.sorted.bam|grep -F -f ${tmp_dir}${one_sample_name}_junction_read_IDs.txt|awk '{if($2<200){print}}')|samtools view -Sbh|samtools sort - >${tmp_dir}${one_sample_name}_junction_read_IDs.bam
							   
							   #index the bam
							   samtools index ${tmp_dir}${one_sample_name}_junction_read_IDs.bam
							   
							   #make a bed file with them (we use -splitD to treat "D" in the cigar as "N" )
							   bedtools bamtobed -bed12 -splitD -i ${tmp_dir}${one_sample_name}_junction_read_IDs.bam >${tmp_dir}${one_sample_name}_junction_read_IDs.bed 
							   
							   #convert them in bam, but for samtools tview (we have to convert all strange characters)
							   samtools view -h -F4 ${tmp_dir}${one_sample_name}_junction_read_IDs.bam|awk 'OFS="\t"{if($1!~/^@SQ/){$3=gensub(":","_","g",$3);$3=gensub("\\|","_","g",$3);$7=gensub(":","_","g",$7);$7=gensub("\\|","_","g",$7);print}else{$2=gensub(":","_","g",$2);$2=gensub("\\|","_","g",$2);$2=gensub("SN_","SN:","g",$2);print}}' |samtools view -Sbh|samtools sort - >${tmp_dir}${one_sample_name}_junction_read_IDs_tview.bam
							   
							   samtools index ${tmp_dir}${one_sample_name}_junction_read_IDs_tview.bam
							   
							   #for each circ detected, give the start, and the length of the left/right most end (for tview)
							   bedtools bamtobed -i ${tmp_dir}${one_sample_name}_junction_read_IDs_tview.bam|cut -f1,2,3|sort -k1,1|bedtools groupby -g 1 -c 2,3 -o min,max|awk 'OFS="\t"{print $1,$2+1,$3,($3-($2+1))+1}'|tee ${tmp_dir}${one_sample_name}_junction_read_IDs_starts_ends.txt|awk -v bam_file=${tmp_dir}${one_sample_name}_junction_read_IDs_tview.bam '{print "COLUMNS="$4" samtools tview -p "$1":"$2"-"$3" -d T "bam_file"|"}' >${tophat_outputs_OneSample}/circ/${one_sample_name}_getConsensus_subscripts_dir/${one_sample_name}_getConsensus_junctions.sh
							   
							   #COLUMNS=250 samtools tview -d T -p chr9_35657751_35658018:129-378 tata.bam|awk '{if(NR>3){print}}'|sed 's/ /-/g'|sed 's/</-/g'|sed 's/>/-/g'|sed 's/\./N/g'|sed 's/\*/-/g'|awk '{print toupper($0)}'|awk 'OFS="\t"{print ">seq_"NR-1"\n"$1"\n"}'
							   
							   nb_rows=$(wc -l ${tmp_dir}${one_sample_name}_junction_read_IDs_starts_ends.txt|awk '{print $1}')
							   
							   #print the samtools tview to awk in order to refine the alignment, used the store IDs to name the fasta file from the alignment ; give it to the script to construct the consensus sequence (from the reads of the junction)
							   paste -d' ' ${tophat_outputs_OneSample}/circ/${one_sample_name}_getConsensus_subscripts_dir/${one_sample_name}_getConsensus_junctions.sh <(yes $(echo -e "awk '{if(NR>3){print}}'|sed 's/ /-/g'|sed 's/</-/g'|sed 's/>/-/g'|sed 's/\./N/g'|sed 's/\*/-/g'|awk '{print toupper(\$0)}'|awk 'OFS=\"\\\t\"{print \">seq_\"NR-1\"\\\n\"\$1\"\\\n\"}' >")|head -n $nb_rows) <(cut -f1 ${tmp_dir}${one_sample_name}_junction_read_IDs_starts_ends.txt|awk -v getConsensus_seq=$getConsensus_seq -v my_path="${tmp_dir}${one_sample_name}_" '{print my_path$1"_getConsensus_fasta_seq.fa ; "getConsensus_seq" "$1" "my_path$1"_getConsensus_fasta_seq.fa"" ""+"" ""no"" >"my_path$1"_obtainedConsensus_fasta_seq.fa"}') >${tophat_outputs_OneSample}/circ/${one_sample_name}_getConsensus_subscripts_dir/${one_sample_name}_getConsensus_junctions.tmp && mv ${tophat_outputs_OneSample}/circ/${one_sample_name}_getConsensus_subscripts_dir/${one_sample_name}_getConsensus_junctions.tmp ${tophat_outputs_OneSample}/circ/${one_sample_name}_getConsensus_subscripts_dir/${one_sample_name}_getConsensus_junctions.sh
							   
							   #for each circ (in parallel), run the command lines in order to have the consensus seq
							   parallel -j $thread_num < ${tophat_outputs_OneSample}/circ/${one_sample_name}_getConsensus_subscripts_dir/${one_sample_name}_getConsensus_junctions.sh
							   
							   #add the same kind of ID as the gtf (circRNA_)
							   find ${tmp_dir} -name "*obtainedConsensus_fasta_seq.fa" | xargs cat|sed 's/>/>circRNA_/g' >${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences.fa
							   
							   #add the strand info (from the refined gtf) to the junction seq : col 1 -> id ; col 2 -> strand ; col 3 -> seq
							   join -t $'\t' -a1 -11 -21 <(grep -v "^#" ${tophat_outputs_OneSample}/${one_sample_name}_refined.gtf |awk 'OFS="\t"{print $10,$7}' |sed 's/\"//g'|sed 's/\;//g'|sed 's/:/_/g'|sed 's/|/_/g'|sort -k1,1) <(cat ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences.fa|awk 'BEGIN{RS=">"}''{if(NR>1){sub("\n","\t"); gsub("\n",""); print $0}}' |awk 'OFS="\t"{print $1,$2}'|sort -k1,1) >${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand.tsv
							   
							   #isolate junction seq on minus strand
							   awk 'OFS="\t"{if($2=="-"){print}}' ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand.tsv >${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand_minus.tsv
							   
							   #rev comp the junct seq on the minus strand
							   #1st manip : take IDs ; 2nd manip : rev comp the seq
							   paste <(cut -f1 ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand_minus.tsv) <(cut -f3 ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand_minus.tsv|tr 'ATGC' 'TACG' |rev) >${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand_minus.tmp && mv ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand_minus.tmp ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand_minus.tsv
							   
							   cat <(awk 'OFS="\t"{if($2=="+"){print ">"$1"\n"$3}}' ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand.tsv) <(awk 'OFS="\t"{print ">"$1"\n"$2}' ${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_with_strand_minus.tsv) >${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_5_3_orientation.fa
							   
							   
							   find ${tmp_dir} -name "*obtainedConsensus_fasta_seq.fa" | while read F;do rm $F;done
							   
							   find ${tmp_dir} -name "*getConsensus_fasta_seq.fa" | while read F;do rm $F;done
							   
							   all_consensus_junction_sequences+=(${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_5_3_orientation.fa)
							   
					else
					
						echo -e "${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_5_3_orientation.fa is already there, next now !\n"
						
						all_consensus_junction_sequences+=(${tophat_outputs_OneSample}/circ/${one_sample_name}_all_consensus_junction_sequences_5_3_orientation.fa)
					
					
					fi
				  
				  
			   fi
			   
		done



		#dir where to find the summary
		if [ ! -d ${output_dir}summarized_counts/ ]; then mkdir ${output_dir}summarized_counts/; fi

		#summarize the counts across samples
		prep_CIRIquant -i ${ciriquant_design} --lib ${output_dir}summarized_counts/library_info.csv --circ ${output_dir}summarized_counts/circRNA_info.csv --bsj ${output_dir}summarized_counts/circRNA_bsj.csv --ratio ${output_dir}summarized_counts/circRNA_ratio.csv

		#transform the infos into tabular data (be careful, in a same cell, you can have comma to delimitate many genes, and they are enclosed in double quotes)
		sed 's/","/\t/g' ${output_dir}summarized_counts/circRNA_info.csv|sed 's/\"//g'|awk 'OFS="\t"{$1=gensub("\\|","_","g",$1);print}'| awk 'OFS="\t"{$1=gensub("\\:","_","g",$1);print}'| awk 'OFS="\t"{if(NR>1){$1="circRNA_"$1};print}'|sed 's/circRNA_circRNA/circRNA/g'|awk '{if(NR==1){$0=gensub(",","\t","g",$0)};print}' >${output_dir}summarized_counts/circRNA_info.tsv

		head -n1 ${output_dir}summarized_counts/circRNA_info.tsv | awk 'OFS="\t"{print $1,$2,"circ_"$3,"circ_"$4,"circ_"$5}' >${output_dir}summarized_counts/circRNA_header.txt

		#concatenate the sequence of all circRNAs
		#we don't take the full sequences provided by CIRI, but just the junction from the reads
		#find ${output_dir} -name "*index.fa"
		cat ${all_consensus_junction_sequences[*]} >${output_dir}summarized_counts/all_circRNAs.fa.tmp && awk 'OFS="\t"{if($1~/^>/){$1=gensub("\\|","_","g",$1);$1=gensub("\\:","_","g",$1)};print}' ${output_dir}summarized_counts/all_circRNAs.fa.tmp|sed 's/circRNA_circRNA/circRNA/g' >${output_dir}summarized_counts/all_circRNAs.fa && rm ${output_dir}summarized_counts/all_circRNAs.fa.tmp

		#put them in tab (ID + seq)
		awk 'BEGIN{RS=">"}''{if(NR>1){sub("\n","\t"); gsub("\n",""); print $0}}' ${output_dir}summarized_counts/all_circRNAs.fa|sed 's/>//g' |sort -u -k1,1 -T ${output_dir}summarized_counts/ >${output_dir}summarized_counts/all_circRNAs.fa.tmp && mv ${output_dir}summarized_counts/all_circRNAs.fa.tmp ${output_dir}summarized_counts/all_circRNAs.fa


		LANG=en_EN join -a1 -t $'\t' -11 -21 <(awk 'OFS="\t"{if(NR>1){print}}' ${output_dir}summarized_counts/circRNA_info.tsv|LANG=en_EN sort -k1,1 -T ${output_dir}summarized_counts/) <(LANG=en_EN sort -k1,1 -T ${output_dir}summarized_counts/ ${output_dir}summarized_counts/all_circRNAs.fa) >${output_dir}summarized_counts/circRNA_info_with_sequences.tsv

		cat <(paste ${output_dir}summarized_counts/circRNA_header.txt <(echo -e "circ_junction_sequence")) ${output_dir}summarized_counts/circRNA_info_with_sequences.tsv >${output_dir}summarized_counts/circRNA_info_with_sequences.tmp && mv ${output_dir}summarized_counts/circRNA_info_with_sequences.tmp ${output_dir}summarized_counts/circRNA_info_with_sequences.tsv

		awk 'OFS="\t"{if(NR>1){print $1,length($NF)}}' ${output_dir}summarized_counts/circRNA_info_with_sequences.tsv|sort -u -k1,1 >${output_dir}summarized_counts/circRNA_length.tsv

		#create the combined GTF file
		#we keep only general informations (we remove infos about counts)
		#we make unique on chr, start, end, as it has be done in the program
		#|sed -E 's/junc_ratio [0-9]+\.[0-9]+;//g'|sed 's/gene_id/circ_gene_id/g'|sed 's/gene_name/circ_gene_name/g'
		cat $(cut -f2 ${ciriquant_design})|grep -v "^#"|sort -k1,1 -k4,4n -k5,5n|sort -u -k1,1 -k4,4n -k5,5n|sed -E 's/gene_type \".*\"/gene_type \"circRNA\"/g'|sed -E 's/bsj [0-9]+\.[0-9]+;//g'|sed -E 's/fsj [0-9]+\.[0-9]+;//g'|sed -E 's/junc_ratio [0-9]+\.[0-9]+;//g'|sed 's/gene_id/circ_gene_id/g'|sed 's/gene_name/circ_gene_name/g'|awk -F'\t' '{print $0" gene_id \"circRNA_"$1"_"$4"_"$5"\""";"}'|awk -F'\t' '{print $0" gene_name \"circRNA_"$1"_"$4"_"$5"\""";"}'|awk -F'\t' '{print $0" transcript_id \"circRNA_"$1"_"$4"_"$5"\""";"}'|awk -F'\t' '{print $0" transcript_name \"circRNA_"$1"_"$4"_"$5"\""";"}' |awk -F'\t' '{print $0" Parent \"circRNA_"$1"_"$4"_"$5"\""";"}' |awk -F'\t' 'OFS="\t"{$3="exon";print}' >${output_dir}summarized_counts/combined_circRNA.gtf


		#create the associated gff
		$gffread --force-exons -F -E ${output_dir}summarized_counts/combined_circRNA.gtf -o- |grep -v "^#" |sort -k1,1 -k4,4n| grep -v "^#"|awk 'OFS="\t"{if($3=="exon"){a=gensub("Parent","ID","g",$9);b=gensub("Parent","gene_name","g",$9);c=gensub("Parent","gene_id","g",$9);d=gensub("Parent","transcript_id","g",$9);print $0";"a";"b";"c";"d";gene_type=circRNA;"}else{$3="gene";print $0";"}}'|sed 's/geneID=/gene_id=/g' >${output_dir}summarized_counts/combined_circRNA.gff3


		#make tsv from the counts in csv
		cat ${output_dir}summarized_counts/circRNA_bsj.csv|tr ',' '\t'|awk 'OFS="\t"{if(NR>1){$1=gensub(":","_","g",$1);$1=gensub("\\|","_","g",$1);print $0}else{print}}' >${output_dir}summarized_counts/circRNA_bsj.tsv

		#retrieve the samples' name
		all_samples=($(awk '{if(NR==1){for(i=2;i<=NF;i++){print $i}}else{exit}}' ${output_dir}summarized_counts/circRNA_bsj.tsv))

		echo -e "all samples are : ${all_samples[*]}\n"

		#for each sample name, make a single table of counts (1st col : ID ; 2nd col : counts)
		 
		for one_pos in $(seq 0 $((${#all_samples[*]}-1)) );do

				 one_id=${all_samples[$one_pos]}
				 
				 #we add +1 to have 1 at the beginning (bash is 0-based position), then we add 1 to have the next column
				 col_pos=$((${one_pos}+1+1))
				 
				 
				 awk '{if(NR>1){print}}' ${output_dir}summarized_counts/circRNA_bsj.tsv|cut -f 1,$col_pos |awk 'OFS="\t"{print $1,$2}' >${output_dir}summarized_counts/${one_id}_CIRIquantCounts.tsv
		 
		done
		

done     

## end of dev ##


#in order to load circView

#https://github.com/GeneFeng/CircView
#~/Downloads/jre1.8.0_261/bin/java -jar /home/marcgabriel/Downloads/CircView/CircView.jar

#-> load annotation | load data (ciri/circexplorer)


#Ularcirc : shiny app for visualization

#https://www.bioconductor.org/packages/release/bioc/vignettes/Ularcirc/inst/doc/Ularcirc.html#splice-junction-files


