# Skeletal_atlas
Scripts used to analyze data for the skeletal atlas. From scRNA-seq, RNA-seq and Proteomics data.

<b>RNA-seq</b>:<br>
Once raw fastq files have been downloaded, we used one of the following scripts:<br>
•	full_script_human_disease_paired.sh<br>
•	full_script_human_disease_single.sh<br>
•	full_script_mouse_disease_paired.sh<br>
•	full_script_mouse_disease_single.sh<br>
<br>
Depending of the organism and the layout of the sequencing. These scripts allow to perform quality check of raw reads, removing low quality reads with Trimmomatic, alignment with Hisat 2 and obtain counts with featureCounts.<br>
<br>
After obtaining the final results with featureCounts, the scripts do_counts.sh and do_lengths.sh allow us to get in the shape of matrix the counts and the lengths for each gene and for each sample.<br>
<br>
Finaly, the script change_matrix_gene_id.pl allow us to switch the ENSEMBL IDs for the gene symbols.<br>
