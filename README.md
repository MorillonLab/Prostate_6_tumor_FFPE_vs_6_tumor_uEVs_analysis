# Prostate_6_tumor_FFPE_vs_6_paired_tumor_uEVs_analysis
Description : Transcriptomic analysis of 6 FFPE tumor samples vs 6 paired uEVs samples (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9081490/)

link to have access to the processed data of the paper : https://filesender.renater.fr/?s=download&token=0886a70d-f868-4463-9475-9377b5a6a696


## Scripts

- Research of HLA alleles from RNAseq data (using seq2HLA) : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/getHLAAllotype2.sh

- Research of neoantigen ability from a supplied list of peptides (using netMHCpan-4.1) : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/getNeoantigen.sh

- Discovery of circRNAs + quantification + junction sequence (using CIRIquant) : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/getCircRNAs.sh

- Build metagenes from bigwig files : https://github.com/mgabriel01/meta_features

- Compute read counts distribution per genomic features : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/getReadsPerFeatureDistrib.sh ; modify the input variables to instert these scripts : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/getIntronsByTranscripts.R & https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/getTranscriptByExonsByStrand.sh


## Contacts 

marc.gabriel@curie.fr (bioinfo)

anna.almeida@curie.fr

antonin.morillon@curie.fr

