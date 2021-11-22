**Prepare work directory and QC**

```
cd /s3_d4/caorui/mouse_RNAseq/
mkdir second_stream_data/rawdata
cd second_stream_data/rawdata
ln -s ../../X201SC20123423-Z01-F004_01/raw_data/* .
ln -s ../../X201SC20123423-Z01-F004_02/raw_data/* .
ls > sampleList.txt ## prepare sampleList

for sample in `cat sampleList2.txt`
do
echo $sample
~/fqc.pl  adp_qual -p -f "$sample"/"$sample"_1.fq -r "$sample"/"$sample"_2.fq -o "$sample" 
done

```
**Build index**
```
nohup rsem-prepare-reference --gtf /s3_d4/caorui/mouse_RNAseq/genome/Mus_musculus.GRCm39.104.gtf /s3_d4/caorui/mouse_RNAseq/genome/Mus_musculus.GRCm39.dna.toplevel.fa --bowtie2 mouse_reference -p 40
```
**Quantify**
```
cd ..
mkdir quantify && cd quantify

for sample in `cat ../rawdata/sampleList.txt`
do      
echo $sample
rsem-calculate-expression --paired-end -no-bam-output --bowtie2 --append-names -p 20 \
../rawdata/"$sample"_1.fastq \
../rawdata/"$sample"_2.fastq \
../rawdata/mouse_reference \
$sample 
done
```
**Combine the matrix** 
use trinity_script to generate transcript-level experssion matrix
```
~/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix mouse_trans ../quantify/*.isoforms.results --gene_trans_map none
```

**Check the replicates**
```
 ~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/PtR -m mouse_trans.isoform.counts.matrix -s samples.txt --log2 --compare_replicates
 ~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/PtR -m mouse_trans.isoform.counts.matrix -s samples.txt --log2 --CPM --prin_comp 3
```
[B6E18.rep_compare.pdf](https://github.com/caorui12/mouse-RNAseq/files/7536361/B6E18.rep_compare.pdf)

[PCA.pdf](https://github.com/caorui12/mouse-RNAseq/files/7536369/mouse_trans.isoform.counts.matrix.CPM.log2.prcomp.principal_components.pdf)

Both figure show B6E18_5F is degarded, thus to be removed in the further analysis

**DE analysis(DEseq2)**
```
~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix mouse_trans.isoform.counts.matrix \
--method DESeq2 \
--samples_file samples.txt
```
these analysis obtain the pairwise comparsion of DE transcirpts

Here I keep two documents, one is DEs from pairwise comparison, one is only sequential comparison(two sequential comparison is to only retain the nearby stages, e.g. E10.5 vs E11.5, E11.5 vs E12.5 etc) 

**Sequential DE analysis** 
```
cd DESeq2.32273.dir(sequencial comparison)
~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../mouse_trans.isoform.TMM.EXPR.matrix --samples ../samples.txt -P 1e-3 -C 2 
```
total 1841 features identified 
**cut the tree**
```
~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData
```
**GO enrichment for DE analysis**
### GO enrichmet
```
library(clusterProfiler)
library(stringr)
library(org.Mm.eg.db)
library(DOSE)
library(topGO)
library(GSEABase)
library(GO.db)
library(ggplot2)
setwd('~/Desktop/RNAseq_practice/second_stream/mouse_genes_level/sample_without_W/DESeq(sequential)')
files=list.files(getwd(),pattern = "DE.subset") ## get the file_list
for (i in 1:length(files)){
target<-read.table(files[i], header=T, row.names=1)
ego <- enrichGO(
  gene  = rownames(target),
  keyType = "ENSEMBL", 
  ont ='ALL',
  OrgDb   = org.Mm.eg.db,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)
file_name=paste(files[i],'_','GO.pdf',sep='')
pdf(file_name)
p<-dotplot(ego,x ='Count',label_format=10)
print(p + scale_y_discrete(labels=function(ego) str_wrap(ego, width=40))) ##print is necessary, otherwise the figure will fail. 

dev.off()
write.table(ego, paste(file_name,'_','GO.matrix', sep=''), sep='\t', row.names = FALSE, quote = FALSE) 
}
```
