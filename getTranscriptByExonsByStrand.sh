

usage() { echo -e "Usage: $0 <Required arguments>\n\n
 \tRequired arguments :\n
                  -a < annotation >\n                 
                
 \tOptional arguments :\n
                  
                  -u < unique on chrs & 9th column (default=no ; others : yes) >\n
                  -e < feature to extract (default=\"exon\") >\n
                  -f < level to put (default is \"transcript\") >\n" 1>&2; exit 1;}

[[ $# -eq 0 ]] && usage

while getopts ":a:u:f:e:" opt; do
  case $opt in
  
    
      u)
      
	      export uniqness=$OPTARG
	      
	      
              ;;              
              
      a)
      
	      annotation=$OPTARG
	      
	      #echo -e "\n#annotation is : $annotation\n" >&2
	      
              ;;
        
      f)
      
	      feature=$OPTARG
	      
	      
              ;;
              
      e)
      
	      feat_to_extract=$OPTARG
	      
	      
              ;;              
              
              
      #invalid options (options not in the list)
      ######################
      
      
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    
  esac
done

bedtools="bedtools"



if [ "$annotation" == "" ]; then
      
echo -e "\none required argument is missing or is wrong !!\n"
      usage
      exit 1
fi


if [[ "$feature" == "" ]];then

  feature="transcript"

fi

if [[ "$feat_to_extract" == "" ]];then

  feat_to_extract="exon"

fi

if [[ "$uniqness" == "" ]];then

  uniqness="no"

fi

if [[ "$uniqness" == "yes" ]];then


	#enforce it by selecting gene_id, gene_name, parent or transcript id, etc
	#we group the exons by the chr, the strand, and the attribute column (here they all have the same, enforce it if necessary by selecting common attributes for the exons of a gene like gene_id & transcript_id)
	#we sort in the same way we want the grouping
	#we remove "Parent"
	#we give as "ID", the "Parent" used by the exons
	#we merge all exons, and take the min start and the min max as new start and end of the new lvl
	#we print the wanted attributes
	#make unique on attribute and chromosome
	grep -P "\t${feat_to_extract}\t" $annotation |sort -k1,1 -k7,7 -k9,9|awk 'OFS="\t"{gene_id="NA";gene_name="NA";Parent="NA";source="NA";gene_type="NA";split($9,a,";");for(i=1;i<=length(a);i++){if(a[i]~/gene_id/){gene_id=a[i]};if(a[i]~/gene_name/){gene_name=a[i]};if(a[i]~/gene_type/){gene_type=a[i]};if(a[i]~/source/){source=a[i]};if(a[i]~/Parent/){Parent=a[i]}};gene_id=gensub("gene_id=","","g",gene_id);gene_name=gensub("gene_name=","","g",gene_name);gene_type=gensub("gene_type=","","g",gene_type);source=gensub("source=","","g",source);Parent=gensub("Parent=","","g",Parent);print $1,$2,$3,$4,$5,$6,$7,$8,"ID="Parent";""gene_name="gene_name";""gene_id="gene_id";""source="source";""gene_type="gene_type}'|sort -k1,1 -k7,7 -k9,9|$bedtools groupby -i stdin -g 1,7,9 -c 4,5 -o min,max|awk -v feature=$feature 'OFS="\t"{print $1,".",feature,$4,$5,".",$2,".",$3}'|awk 'OFS="\t"{ID="NA";gene_id="NA";gene_name="NA";Parent="NA";source="NA";gene_type="NA";split($9,a,";");for(i=1;i<=length(a);i++){if(a[i]~/ID=/){ID=a[i]};if(a[i]~/gene_id/){gene_id=a[i]};if(a[i]~/gene_name/){gene_name=a[i]};if(a[i]~/gene_type/){gene_type=a[i]};if(a[i]~/source/){source=a[i]};if(a[i]~/Parent/){Parent=a[i]}};ID=gensub("ID=","","g",ID);gene_id=gensub("gene_id=","","g",gene_id);gene_name=gensub("gene_name=","","g",gene_name);gene_type=gensub("gene_type=","","g",gene_type);source=gensub("source=","","g",source);Parent=gensub("Parent=","","g",Parent);if(ID=="NA"){ID=gene_id};print $1,$2,$3,$4,$5,$6,$7,$8,"ID="ID";""gene_name="gene_name";""gene_id="gene_id";""source="source";""gene_type="gene_type}' |sort -u -k1,1 -k9,9|sort -k1,1 -k4,4n

else

	grep -P "\t${feat_to_extract}\t" $annotation |sort -k1,1 -k7,7 -k9,9|awk 'OFS="\t"{gene_id="NA";gene_name="NA";Parent="NA";source="NA";gene_type="NA";split($9,a,";");for(i=1;i<=length(a);i++){if(a[i]~/gene_id/){gene_id=a[i]};if(a[i]~/gene_name/){gene_name=a[i]};if(a[i]~/gene_type/){gene_type=a[i]};if(a[i]~/source/){source=a[i]};if(a[i]~/Parent/){Parent=a[i]}};gene_id=gensub("gene_id=","","g",gene_id);gene_name=gensub("gene_name=","","g",gene_name);gene_type=gensub("gene_type=","","g",gene_type);source=gensub("source=","","g",source);Parent=gensub("Parent=","","g",Parent);print $1,$2,$3,$4,$5,$6,$7,$8,"ID="Parent";""gene_name="gene_name";""gene_id="gene_id";""source="source";""gene_type="gene_type}'|sort -k1,1 -k7,7 -k9,9|$bedtools groupby -i stdin -g 1,7,9 -c 4,5 -o min,max|awk -v feature=$feature 'OFS="\t"{print $1,".",feature,$4,$5,".",$2,".",$3}'|awk 'OFS="\t"{ID="NA";gene_id="NA";gene_name="NA";Parent="NA";source="NA";gene_type="NA";split($9,a,";");for(i=1;i<=length(a);i++){if(a[i]~/ID=/){ID=a[i]};if(a[i]~/gene_id/){gene_id=a[i]};if(a[i]~/gene_name/){gene_name=a[i]};if(a[i]~/gene_type/){gene_type=a[i]};if(a[i]~/source/){source=a[i]};if(a[i]~/Parent/){Parent=a[i]}};ID=gensub("ID=","","g",ID);gene_id=gensub("gene_id=","","g",gene_id);gene_name=gensub("gene_name=","","g",gene_name);gene_type=gensub("gene_type=","","g",gene_type);source=gensub("source=","","g",source);Parent=gensub("Parent=","","g",Parent);if(ID=="NA"){ID=gene_id};print $1,$2,$3,$4,$5,$6,$7,$8,"ID="ID";""gene_name="gene_name";""gene_id="gene_id";""source="source";""gene_type="gene_type}' |sort -k1,1 -k4,4n



fi
