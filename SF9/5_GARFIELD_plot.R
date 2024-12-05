

# firstorder
type=firstorder_minP
DATADIR=./garfield
PRUNETAGSDIR=$DATADIR/tags/r01
CLUMPTAGSDIR=$DATADIR/tags/r08
MAFTSSDDIR=$DATADIR/maftssd
PVALDIR=$DATADIR/pval/$type
ANNOTDIR=/data1/qtl/3QTL/enrichment/garfield-data/annotation
OUTDIR=$DATADIR/output/$type
mkdir -p $OUTDIR

ANNOTLINKFILE=$ANNOTDIR/link_file.txt
PTHRESH=5e-2,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,5e-8
BINNING=m5,n5,t5
CONDITION=0

# loop处理每个染色体，生成 garfield.prep.out 文件
F1=$OUTDIR/$type.garfield.prep.out
F0=$OUTDIR/$type.garfield.Meff.out

echo -n > $F1
for CHR in `seq 1 22` #X
do
echo 'CHR'$CHR
./garfield-prep-chr \
-ptags $PRUNETAGSDIR/chr$CHR \
-ctags $CLUMPTAGSDIR/chr$CHR \
-maftss $MAFTSSDDIR/chr$CHR \
-pval $PVALDIR/chr$CHR \
-ann $ANNOTDIR/chr$CHR \
-chr $CHR -o $F1 || { echo 'Failure!'; } 
done

# 计算Meff（多重比较中，调整后的有效检验数量）和 Padj
Rscript garfield-Meff-Padj.R -i $F1 -o $F0 
NEA=$(head -1 $F0 |awk '{print $2}')
Padj=$(tail -1 $F0 |awk '{print $2}')

# garfield-test.R 脚本执行测试分析，输入文件是 F1（预处理后的数据文件），输出文件为 F2（测试结果文件）
F2=$OUTDIR/$type.garfield.test.out 
Rscript garfield-test.R -i $F1 -o $F2 -l $ANNOTLINKFILE -pt $PTHRESH -b $BINNING -c $CONDITION

# plot
Rscript garfield-plot.R  -i $F2 -o $OUTDIR/$type -t "" -f 0 -padj 0 




# shape
type=shape_minP
DATADIR=./garfield
PRUNETAGSDIR=$DATADIR/tags/r01
CLUMPTAGSDIR=$DATADIR/tags/r08
MAFTSSDDIR=$DATADIR/maftssd
PVALDIR=$DATADIR/pval/$type
ANNOTDIR=/data1/qtl/3QTL/enrichment/garfield-data/annotation
OUTDIR=$DATADIR/output/$type
mkdir -p $OUTDIR

ANNOTLINKFILE=$ANNOTDIR/link_file.txt
PTHRESH=5e-2,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,5e-8
BINNING=m5,n5,t5
CONDITION=0

F1=$OUTDIR/$type.garfield.prep.out
F0=$OUTDIR/$type.garfield.Meff.out

echo 'Prune and Clump'
echo -n > $F1
for CHR in `seq 1 22` #X
do
echo 'CHR'$CHR
./garfield-prep-chr \
-ptags $PRUNETAGSDIR/chr$CHR \
-ctags $CLUMPTAGSDIR/chr$CHR \
-maftss $MAFTSSDDIR/chr$CHR \
-pval $PVALDIR/chr$CHR \
-ann $ANNOTDIR/chr$CHR \
-chr $CHR -o $F1 || { echo 'Failure!'; } 
done

Rscript garfield-Meff-Padj.R -i $F1 -o $F0
NEA=$(head -1 $F0 |awk '{print $2}')
Padj=$(tail -1 $F0 |awk '{print $2}')

F2=$OUTDIR/$type.garfield.test.out
Rscript garfield-test.R -i $F1 -o $F2 -l $ANNOTLINKFILE -pt $PTHRESH -b $BINNING -c $CONDITION

# plot
Rscript garfield-plot.R  -i $F2 -o $OUTDIR/$type -t "" -f 0 -padj 0 
