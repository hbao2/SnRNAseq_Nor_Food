#!/bin/bash
#!/bin/bash
#SBATCH --job-name=cellranger            # Job name
#SBATCH --ntasks=1                    # Number of tasks (1 task with multiple cores)
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --mem=256G                    # Memory per node
#SBATCH --time=48:00:00               # Time limit hrs:min:sec (adjust as needed)


cd /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60

# add cellranger
export PATH=/media/XLStorage/hbao2/SCRNA/cellranger-8.0.1:$PATH

FASTQS_DIR=/media/XLStorage/hbao2/SCRNA/Filterd_DATA/outs/fastq_path
TRANSCRIPTOME_DIR=/media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/FCA660_Cellranger801

# Run Cell Ranger count for each sample and echo "finish" upon completion

# cellranger count --id=Female_Ctrl_B \
# --sample=24252-01-01-01-01 \
# --fastqs=$FASTQS_DIR \
# --transcriptome=$TRANSCRIPTOME_DIR \
# --include-introns true \
# --create-bam true \
# --localcores 20 
# echo "Finished Female_Ctrl_B"

cellranger count --id=Female_Host_B \
--sample=24252-01-02-01-01 \
--fastqs=$FASTQS_DIR \
--transcriptome=$TRANSCRIPTOME_DIR \
--include-introns true \
--create-bam true \
--localcores 20 
echo "Finished Female_Host_B"

# cellranger count --id=Male_Ctrl_B \
# --sample=24252-01-03-01-01 \
# --fastqs=$FASTQS_DIR \
# --transcriptome=$TRANSCRIPTOME_DIR \
# --include-introns true \
# --create-bam true \
# --localcores 20 
# echo "Finished Male_Ctrl_B"

# cellranger count --id=Male_Host_B \
# --sample=24252-01-04-01-01 \
# --fastqs=$FASTQS_DIR \
# --transcriptome=$TRANSCRIPTOME_DIR \
# --include-introns true \
# --create-bam true \
# --localcores 20 
# echo "Finished Male_Host_B"

# file transfer
# rsync -av --exclude='*.bam' --exclude='*.bai' /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60/Female_Ctrl_B/ denglab@129.81.246.29:/mnt/DATA/RNAseq/SCRNA/cellranger_output/20241020_6.60_All/Female_Ctrl_B/
# rsync -av --exclude='*.bam' --exclude='*.bai' /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60/Female_Host_B/ denglab@129.81.246.29:/mnt/DATA/RNAseq/SCRNA/cellranger_output/20241020_6.60_All/Female_Host_B/
# rsync -av --exclude='*.bam' --exclude='*.bai' /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60/Male_Ctrl_B/ denglab@129.81.246.29:/mnt/DATA/RNAseq/SCRNA/cellranger_output/20241020_6.60_All/Male_Ctrl_B/
# rsync -av --exclude='*.bam' --exclude='*.bai' /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60/Male_Host_B/ denglab@129.81.246.29:/mnt/DATA/RNAseq/SCRNA/cellranger_output/20241020_6.60_All/Male_Host_B/


#----chemistry=SC3Pv1 
# error message:
#There were no reads to process. Please check whether the read lengths in the input fastqs satisfy the minumum read length requirements for the chemistry.

#summary
#https://www.10xgenomics.com/analysis-guides/quality-assessment-using-the-cell-ranger-web-summary

# # ###################################################################
# # # create gene table
# grep -E 'FlyBase.*gene' dmel-all-r6.60_EGFP_GAL4_mCD8GFP.gtf \
# | awk -F'\t' '{print $9}' | awk -F';' '{print $1"\t"$2}' \
# | sed 's/gene_id//g' |sed 's/gene_symbol//g' |sed 's/"//g'|sed 's/ //g' \
# | uniq > Gene_ID.txt
# # EGFP	EGFP
# # mCD8GFP	mCD8GFP
# # GAL4	GAL4
# # ###################################################################




# # mCD8GFP on positive 
# # GAL4 on negtive 
# # GAL80ts on positive
# ###################################################################
# # check the +/- for mCD8GFP
samtools view /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/Female_Host_B/outs/possorted_genome_bam.bam "mCD8GFP" \
| awk '{if ($2 == 16) {minus++} else if ($2 == 0) {plus++}} END {print "plus:", plus, "minus:", minus}'
samtools view /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60/Female_Host_B/outs/possorted_genome_bam.bam "mCD8GFP" \
| awk '{if ($2 == 16) {minus++} else if ($2 == 0) {plus++}} END {print "plus:", plus, "minus:", minus}'
# # 6.60 : plus: 308334 minus: 12171
# # 6.31  plus: 69975 minus: 32549
# # 6.60 : plus: 69530 minus: 32564


# # GAL80ts +
samtools view /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60/Female_Host_B/outs/possorted_genome_bam.bam "GAL80ts" \
| awk '{if ($2 == 16) {minus++} else if ($2 == 0) {plus++}} END {print "plus:", plus, "minus:", minus}'
# # none

# # GAL4 -

samtools view /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/Female_Host_B/outs/possorted_genome_bam.bam "GAL4" \
| awk '{if ($2 == 16) {minus++} else if ($2 == 0) {plus++}} END {print "plus:", plus, "minus:", minus}'
samtools view /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60/Female_Host_B/outs/possorted_genome_bam.bam "GAL4" \
| awk '{if ($2 == 16) {minus++} else if ($2 == 0) {plus++}} END {print "plus:", plus, "minus:", minus}'
# # plus: 69530 minus: 32564
# # 6.31+ : plus: 9444 minus: 21806
# # 6.60- : plus: 52122 minus: 6149


# # EGFP nothing
samtools view /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60/Female_Host_B/outs/possorted_genome_bam.bam "GAL80ts" \
| awk '{if ($2 == 16) {minus++} else if ($2 == 0) {plus++}} END {print "plus:", plus, "minus:", minus}'

# # gene LDH on negtive  3L:6,259,105..6,262,694 [-]
# samtools view -c -F 0x10 possorted_genome_bam.bam 3L:6259105-6262694  # positive 40563
# samtools view -c -f 0x10 possorted_genome_bam.bam 3L:6259105-6262694  # -f nagetive Strand # 132562

# # CG43168 on positive 3L:6,275,787..6,277,035 [+]
# samtools view -c -F 0x10 possorted_genome_bam.bam 3L:6275787-6277035 #1526
# samtools view -c -f 0x10 possorted_genome_bam.bam 3L:6275787-6277035 #943

# #read length 151
# samtools view possorted_genome_bam.bam | awk '{print length($10)}'

# ###################################################################

# ###################################################################
# # cannot tell for bulk RNA seq
# samtools index TM-1_S10Aligned.sortedByCoord.out.bam
# samtools view -c -F 0x10 TM-1_S10Aligned.sortedByCoord.out.bam "EGFP" 
# samtools view -c -f 0x10 TM-1_S10Aligned.sortedByCoord.out.bam "EGFP" 
# # bulk read length = 100 
# samtools view TM-1_S10Aligned.sortedByCoord.out.bam | awk '{print length($10)}'
# ###################################################################


# check in adata mCD8GFP counts
```python
# Check if mCD8GFP exists in .var (genes/features)
if "mCD8GFP" in adata.var_names:
    # Filter cells where 'sample' column in .obs is 'F_host'
    f_host_adata = adata[adata.obs['sample'] == 'F_Host']
    
    # Get the index of mCD8GFP in .var
    index = f_host_adata.var_names.get_loc("mCD8GFP")
    
    # Sum counts across all filtered cells for mCD8GFP
    total_counts = f_host_adata.X[:, index].sum()
    print(f"Total counts for mCD8GFP in F_Host: {total_counts}")
else:
    print("mCD8GFP not found in the .var names.")

```


```python
# Check if GAL4 exists in .var (genes/features)
if "GAL4" in adata.var_names:
    # Filter cells where 'sample' column in .obs is 'F_host'
    f_host_adata = adata[adata.obs['sample'] == 'F_Host']
    
    # Get the index of GAL4 in .var
    index = f_host_adata.var_names.get_loc("GAL4")
    
    # Sum counts across all filtered cells for GAL4
    total_counts = f_host_adata.X[:, index].sum()
    print(f"Total counts for GAL4 in F_Host: {total_counts}")
else:
    print("GAL4 not found in the .var names.")



# check sequence before gal4

cd /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.60/
samtools view Female_Host_B/outs/possorted_genome_bam.bam "GAL4" > GAL4_reads.sam
samtools fastq GAL4_reads.sam > GAL4_reads.fastq

samtools view -b Female_Host_B/outs/possorted_genome_bam.bam "GAL4" > GAL4_reads.bam
samtools fastq GAL4_reads.bam > GAL4_reads.fastq

samtools view -f 4 Female_Host_B/outs/possorted_genome_bam.bam > unmapped_reads.sam
samtools fastq unmapped_reads.sam > unmapped_reads.fastq
cat GAL4_reads.fastq unmapped_reads.fastq > combined_reads.fastq
fastqc combined_reads.fastq -o qc_output/
trimmomatic SE -phred33 combined_reads.fastq cleaned_combined_reads.fastq SLIDINGWINDOW:4:20 MINLEN:50





