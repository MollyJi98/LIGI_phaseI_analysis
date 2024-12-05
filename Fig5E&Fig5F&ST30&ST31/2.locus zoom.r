################手动绘制locus zoom#################
conda activate locusZoom_stand_alone


#提取要画图的位点上下游500kb的位点
#rs41268920 original_firstorder_90Percentile_CR_assoc
awk -v pos=32133344 -v flank=500000 'BEGIN{FS=","; OFS="\t"; print "MarkerName", "P-value"} $2 == "6" && $3 >= (pos - flank) && $3 <= (pos + flank) {print $18, $16}'./original_firstorder_90Percentile_CR_assoc.result > "90Percentile_CR_rs41268920.txt"
#rs41268920 COPD
awk -v pos=32133344 -v flank=500000 'BEGIN{FS="\t"; OFS="\t"; print "MarkerName", "P-value"} $1 == "6" && $2 >= (pos - flank) && $2 <= (pos + flank) {print $9, $6}' ./COPD_for_LDSC.txt > "COPD_rs41268920.txt"

#rs4693974 Maximum2DDiameterColumn_UR
awk -v pos=89810395 -v flank=500000 'BEGIN{FS=","; OFS="\t"; print "MarkerName", "P-value"} $2 == "4" && $3 >= (pos - flank) && $3 <= (pos + flank) {print $18, $16}'./original_shape_Maximum2DDiameterColumn_UR_assoc.result > "Maximum2DDiameterColumn_UR_rs4693974.txt"
#rs4693974 NSCLC
awk -v pos=89810395 -v flank=500000 'BEGIN{FS="\t"; OFS="\t"; print "MarkerName", "P-value"} $1 == "4" && $2 >= (pos - flank) && $2 <= (pos + flank) {print $9, $6}' ./NSCLC_for_LDSC.txt > "NSCLC_rs4693974.txt"


cd /data1/Imaging_Assoc_2023/calculation_result_JC_new/12.locuszoom/forFig4



#rs41268920 original_firstorder_90Percentile_CR_assoc
./locusZoom_stand_alone/locuszoom/bin/locuszoom \
    --metal 90Percentile_CR_rs41268920.txt \
    --plotonly \
    --build hg19 \
    --pop ASN \
    --source 1000G_March2012 \
    --flank 300kb \
    --gene-table gencode \
    --markercol MarkerName \
    --pvalcol Pvalue \
    --refsnp rs41268920 width=15 height=4 showPartialGenes=FALSE rfrows=1 legend='right'

#rs41268920 COPD
./locusZoom_stand_alone/locuszoom/bin/locuszoom \
    --metal COPD_rs41268920.txt \
    --plotonly \
    --build hg19 \
    --pop ASN \
    --source 1000G_March2012 \
    --flank 300kb \
    --gene-table gencode \
    --markercol MarkerName \
    --pvalcol P-value \
    --prefix  COPD\
    --refsnp rs41268920 width=15 height=4 showPartialGenes=FALSE rfrows=1 legend='right'


#rs4693974 Maximum2DDiameterColumn_UR
./locusZoom_stand_alone/locuszoom/bin/locuszoom \
    --metal Maximum2DDiameterColumn_UR_rs4693974.txt \
    --plotonly \
    --build hg19 \
    --pop ASN \
    --source 1000G_March2012 \
    --flank 300kb \
    --gene-table gencode \
    --markercol MarkerName \
    --pvalcol P-value \
    --prefix  Maximum2DDiameterColumn_UR\
    --refsnp rs4693974 width=15 height=4 showPartialGenes=FALSE rfrows=1 legend='right'

#rs4693974 NSCLC
./locusZoom_stand_alone/locuszoom/bin/locuszoom \
    --metal NSCLC_rs4693974.txt \
    --plotonly \
    --build hg19 \
    --pop ASN \
    --source 1000G_March2012 \
    --flank 300kb \
    --gene-table gencode \
    --markercol MarkerName \
    --pvalcol P-value \
    --prefix  NSCLC\
    --refsnp rs4693974 width=15 height=4 showPartialGenes=FALSE rfrows=1 legend='right'


#############导出后后在ai里拼图。区域在染色体上的定位图在ensemble上在线导出#############