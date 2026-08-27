# Internship-Lode-Usja
My internship at the Sophie Janssens Lab in the VIB-IRC.

During this internship, we tried to optimize the integration of the samples and the treatment groups using totalVI, scVI, scANVI and Harmony. In the end, the Harmony integration was used for this step, as it represented the biology, while removing as much of the technical variability as possible. Next, we improved the celltype annotation based on previous attempts, and using the annotation of the individual datasets. This resulted in clusters representing different stages of the developmental trajectory of cDC1. 

Subsequently, we performed downstream processing to use these results in the RShiny tool, allowing others to interact with this dataset. The downstream processing included:
- DESeq2 analysis on the pseudobulk samples. Three general comparisons were made: Toxoplasma treatment vs WT, Toxoplasma treatment vs other (individual) immunogenic treatments (LNP) and individual LNP treatments vs WT
- TF inference based on DESeq2 data and raw counts data
- Conserved markers analysis based on Wilcoxon-rank test results
- GSEA and ORA for pathway inference
- Module analysis using (py)UCell

Finally, the RShiny app was created with the following tabs:
- Gene plots
- Metadata plots
- Marker filtering
- Forest plots for DESeq2 results
- TF inference results
- UCell module finding (does not work properly)
- Seurat based score module
- Seurat based scoring module based on target genes of given TFs
