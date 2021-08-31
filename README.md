# Prostate_6_tumor_FFPE_vs_6_paired_tumor_uEVs_analysis
Description : Transcriptomic analysis of 6 FFPE tumor samples vs 6 paired uEVs samples

- Concatenated annotation of gencode32 + repeats + circRNAs (from 6 FFPE & 6 paired uEVs tumor samples) : http://xfer.curie.fr/get/fkwpKMDUSrS/concatenated_annot_gencode32_repeat_masker_circRNA.zip :

- Read counts from concatenated annotation of gencode32 + repeats + circRNAs (from 6 FFPE & 6 paired uEVs tumor samples) : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/gencode32_repeatMasker_circRNAs_FFPE_urines_counts.zip ; (**used in figure 1g & as Extended Data Table 2**)

- Differential expression table from DESeq2 between 6 tumor FFPE & 6 paired tumor uEVs : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/DESeq_output_FFPE_T_vs_urines_T.zip (features with full 0 counts are filtered out)

- Annotation of 54649 circRNAs from EVs of LNCaP, DU145 & PC3 cell lines : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/circRNAs_EVs_LNCap_DU_PC3_v2.zip

- Annotation of 31236 circRNAs from LNCaP, DU145 & PC3 cell lines : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/circRNAs_cellLines_LNCap_DU_PC3_v2.zip

- Annotated NetMHCpan results (neoantigens) from upregulated lncRNAs (highest expressed transcript per gene) for FFPE samples & uEVs samples : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/neoantigens_from_upregulated_lncRNAs_uEVs_FFPE.zip (**used as Extended Data Table 8**)

- Annotation of 309 common circRNAs in the overlap between 311 up uEVs circRNAs & 31236 circRNAs from LNCaP, DU145 & PC3 cell lines  : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/309_common_311_up_uEVs_26182_LNCaP_DU_PC3_cellLines.gff

- Annotation of 308 common circRNAs in the overlap between 311 up uEVs circRNAs & 54649 circRNAs from EVs of LNCaP, DU145 & PC3 cell lines : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/308_common_311_up_uEVs_26182_LNCaP_DU_PC3_EVs.gff


- 1248 expressed uEVs genes, specific to prostate : https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/1248_expressed_genes_uEVs_specific_to_prostate.xlsx (**used as Extended Data Table 3**)

- Samples information and quality control of RNA and RNAseq :
https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/Extended%20Data%20table%201_samples-QC%20RNA-QC%20RNAseq.xlsx (**used as Extended Data Table 1**)

- 4925 RNAs up in uEVs and 8336 RNAs up in Tumors (padj<=0.05; FC>=1.5): https://github.com/MorillonLab/Prostate_6_tumor_FFPE_vs_6_tumor_uEVs_analysis/blob/main/Extended%20Table%204_13261%20sig_diff_FFPE_T_vs_urines_T.xlsx (**used as Extended Data Table 4**)

