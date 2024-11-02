# edit FASTA
# https://github.com/FlyCellAtlas/genome_references/tree/master/flybase/r6.31

cd /media/XLStorage/hbao2/SCRNA/Ref/
mkdir FCA_6.60/RAW
cd FCA_6.60
wget http://ftp.flybase.net/releases/current/dmel_r6.60/gtf/dmel-all-r6.60.gtf.gz
gunzip dmel-all-r6.60.gtf.gz
wget http://ftp.flybase.net/releases/current/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz
gunzip dmel-all-chromosome-r6.60.fasta.gz


cp /media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/RAW/dmel-all-chromosome-r6.31.fasta \
/media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/dmel-all-chromosome-r6.31_EGFP_GAL4_mCD8GFP.fasta
cat /media/XLStorage/hbao2/SCRNA/Ref/Other/EGFP.GAL4.mCD8GFP.fa >> /media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/dmel-all-chromosome-r6.31_EGFP_GAL4_mCD8GFP.fasta
grep ">" /media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/dmel-all-chromosome-r6.31_EGFP_GAL4_mCD8GFP.fasta

# edit GTF
cp /media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/RAW/dmel-all-r6.31.gtf /media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/dmel-all-r6.31_EGFP_GAL4_mCD8GFP.gtf
cat /media/XLStorage/hbao2/SCRNA/Ref/Other/EGFP_GAL4_mCD8GFP.gtf >> /media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/dmel-all-r6.31_EGFP_GAL4_mCD8GFP.gtf
tail -20 /media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/dmel-all-r6.31_EGFP_GAL4_mCD8GFP.gtf

# add cellranger
export PATH=/media/XLStorage/hbao2/SCRNA/cellranger-8.0.1:$PATH

GENOME_FILE=FCA660_Cellranger801
INPUT_GTF=/media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/dmel-all-r6.31_EGFP_GAL4_mCD8GFP.gtf
INPUT_FA=/media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/dmel-all-chromosome-r6.31_EGFP_GAL4_mCD8GFP.fasta


# Remove entries w/o strand information
# sed -i '/\tgene\t/d' $INPUT_GTF | awk '{ if($6 == "." && $7 == "." && $8 == ".") {} else {print}}'

awk -F'\t' '{ if($7 != ".") print }' $INPUT_GTF > fixed_file.gtf
mv fixed_file.gtf $INPUT_GTF

# fixed_file.gtf

# make ref
cd /media/XLStorage/hbao2/SCRNA/Ref/FCA_6.60/
cellranger mkref \
--nthreads 10 \
--genome=$GENOME_FILE --fasta=$INPUT_FA --genes=$INPUT_GTF

