#!/bin/bash

cd /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31


# add cellranger
export PATH=/media/XLStorage/hbao2/SCRNA/cellranger-8.0.1:$PATH

FASTQS_DIR=/media/XLStorage/hbao2/SCRNA/Filterd_DATA/outs/fastq_path
TRANSCRIPTOME_DIR=/media/XLStorage/hbao2/SCRNA/Ref/FCA_6.31/FCA631_Cellranger801

# Run Cell Ranger count for each sample and echo "finish" upon completion

cellranger count --id=Female_Ctrl_B \
--sample=24252-01-01-01-01 \
--fastqs=$FASTQS_DIR \
--transcriptome=$TRANSCRIPTOME_DIR \
--include-introns true \
--create-bam true \
--localcores 20 
echo "Finished Female_Ctrl_B"

cellranger count --id=Female_Host_B \
--sample=24252-01-02-01-01 \
--fastqs=$FASTQS_DIR \
--transcriptome=$TRANSCRIPTOME_DIR \
--include-introns true \
--create-bam true \
--localcores 20 
echo "Finished Female_Host_B"

cellranger count --id=Male_Ctrl_B \
--sample=24252-01-03-01-01 \
--fastqs=$FASTQS_DIR \
--transcriptome=$TRANSCRIPTOME_DIR \
--include-introns true \
--create-bam true \
--localcores 20 
echo "Finished Male_Ctrl_B"

cellranger count --id=Male_Host_B \
--sample=24252-01-04-01-01 \
--fastqs=$FASTQS_DIR \
--transcriptome=$TRANSCRIPTOME_DIR \
--include-introns true \
--create-bam true \
--localcores 20 
echo "Finished Male_Host_B"

# file transfer
# rsync -av --exclude='*.bam' --exclude='*.bai' /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/Female_Ctrl_B/ denglab@129.81.246.29:/mnt/DATA/RNAseq/SCRNA/cellranger_output/20241020_6.31_All/Female_Ctrl_B/
# rsync -av --exclude='*.bam' --exclude='*.bai' /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/Female_Host_B/ denglab@129.81.246.29:/mnt/DATA/RNAseq/SCRNA/cellranger_output/20241020_6.31_All/Female_Host_B/
# rsync -av --exclude='*.bam' --exclude='*.bai' /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/Male_Ctrl_B/ denglab@129.81.246.29:/mnt/DATA/RNAseq/SCRNA/cellranger_output/20241020_6.31_All/Male_Ctrl_B/
# rsync -av --exclude='*.bam' --exclude='*.bai' /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/Male_Host_B/ denglab@129.81.246.29:/mnt/DATA/RNAseq/SCRNA/cellranger_output/20241020_6.31_All/Male_Host_B/


#----chemistry=SC3Pv1 
# error message:
#There were no reads to process. Please check whether the read lengths in the input fastqs satisfy the minumum read length requirements for the chemistry.

#summary
#https://www.10xgenomics.com/analysis-guides/quality-assessment-using-the-cell-ranger-web-summary

# # ###################################################################
# # # create gene table
# grep -E 'FlyBase.*gene' dmel-all-r6.31_EGFP_GAL4_mCD8GFP.gtf \
# | awk -F'\t' '{print $9}' | awk -F';' '{print $1"\t"$2}' \
# | sed 's/gene_id//g' |sed 's/gene_symbol//g' |sed 's/"//g'|sed 's/ //g' \
# | uniq > Gene_ID.txt
# # EGFP	EGFP
# # mCD8GFP	mCD8GFP
# # GAL4	GAL4
# # ###################################################################




# # mCD8GFP on positive 
# # GAL4 on positive 
# ###################################################################
# # check the +/- for mCD8GFP
# samtools view /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/Female_Host_B/outs/possorted_genome_bam.bam "mCD8GFP" \
# | awk '{if ($2 == 16) {minus++} else if ($2 == 0) {plus++}} END {print "plus:", plus, "minus:", minus}'
# # plus: 308334 minus: 12171

# samtools view /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/Female_Host_B/outs/possorted_genome_bam.bam "GAL4" \
# | awk '{if ($2 == 16) {minus++} else if ($2 == 0) {plus++}} END {print "plus:", plus, "minus:", minus}'
# # plus: 52122 minus: 6149

# samtools view /media/XLStorage/hbao2/SCRNA/Cellranger_output_20241020_6.31/Female_Host_B/outs/possorted_genome_bam.bam "EGFP" \
# | awk '{if ($2 == 16) {minus++} else if ($2 == 0) {plus++}} END {print "plus:", plus, "minus:", minus}'

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


