****(1) Prepare work directory and QC
```
cd /s3_d4/caorui/mouse_RNAseq/
mkdir second_stream_data/rawdata
cd second_stream_data/rawdata
ln -s ../../X201SC20123423-Z01-F004_01/raw_data/* .
ln -s ../../X201SC20123423-Z01-F004_02/raw_data/* .
ls > sampleList.txt ## prepare sampleList

for sample in `cat sampleList.txt`
do
echo $sample
~/fqc.pl  adp_qual -p -f "$sample"/"$sample"_1.fq -r "$sample"/"$sample"_2.fq -o "$sample" 
done

```
**(2) Build STAR index
```
nohup rsem-prepare-reference --gtf /s3_d4/caorui/mouse_RNAseq/genome/Mus_musculus.GRCm39.104.gtf /s3_d4/caorui/mouse_RNAseq/genome/Mus_musculus.GRCm39.dna.toplevel.fa --STAR mouse_reference -p 40
```
**(3) using RSEM to Quantify gene expresison
```
cd ..
mkdir quantify && cd quantify

for sample in `cat ../rawdata/sampleList.txt`
do      
echo $sample
rsem-calculate-expression --paired-end -no-bam-output --STAR --append-names -p 20 \
../rawdata/"$sample"_1.fastq \
../rawdata/"$sample"_2.fastq \
../rawdata/mouse_reference \
$sample 
done
```
**Combine the count matrix and TPM matrix** 
```
rsem-generate-data-matrix *gene* > mouse_trans.isoform.counts.matrix
python TPM_extract.py *genes*>TPM.matrix
```

**Revise the gene name and only keep the coding gene
```
join ../../genome/final_coding_gene.txt mouse_trans.isoform.counts.matrix -t $ '\t' > coding.counts.matrix
join ../../genome/final_coding_gene.txt TPM.matrix1 -t $'\t' > coding.TPM.matrix
```
**Check the replicates**
```
 ~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/PtR -m mouse_trans.isoform.counts.matrix -s samples.txt --log2 --compare_replicates
 ~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/PtR -m mouse_trans.isoform.counts.matrix -s samples.txt --log2 --CPM --prin_comp 3
```

**DE analysis(DEseq2)**
```
~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix mouse_trans.isoform.counts.matrix \
--method DESeq2 \
--samples_file samples.txt
```
these analysis obtain the pairwise comparsion of DE transcirpts


**Extract DE genes
```
cd DESeq2.32273.dir(sequencial comparison)
~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../mouse_trans.isoform.TMM.EXPR.matrix --samples ../samples.txt -P 1e-3 -C 1 
```


**clustering analysis**
```
library(cluster)
library(highcharter)
library(reshape2)
library(gplots)
source("~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/R/heatmap.3.R")
source("~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/R/misc_rnaseq_funcs.R")
source("~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/R/pairs3.R")
source("~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/R/vioplot2.R")
setwd('~/Dropbox/mouseRNAseq/145N1N2remove/DESeq2sequential/')
### prepare data, scale data
DF<-read.table('diffExpr.P1e-3_C1.matrix',header=T,row.names = 1) ### read the DE matrix
DF<-DF[,-c(1,2)] ### two sample need to be removed
data<-as.matrix(DF)
data<-t(scale(t(data)))

### prepare annotation 
myheatcol = colorpanel(75, 'blue','white','red')
samples_data = read.table("../sampleList.txt", header=F, check.names=F, fill=T)
samples_data = samples_data[samples_data[,2] != '',]
colnames(samples_data) = c('sample_name', 'replicate_name')
sample_types = as.character(unique(samples_data[,1]))
rep_names = as.character(samples_data[,2])
data = data[, colnames(data) %in% rep_names, drop=F ]
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
  samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
  sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
  sample_type = sample_types[i]
  replicates_want = sample_type_list[[sample_type]]
  sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)

### calculate distance 
gene_dist = dist(data, method='euclidean')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')
hc_genes<-hclust(gene_dist)


### split the gene cluster 
height=95 ### please specific the parameter of h (0-100)
gene_partition_assignments <- cutree(as.hclust(hc_genes),h=height/100*max(hc_genes$height)) ### please specific the parameter of h
write.table(gene_partition_assignments[hc_genes$order], file="clusters_fixed_P_90.heatmap.heatmap_gene_order.txt", quote=F, sep='	')
max_cluster_count = max(gene_partition_assignments)
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors_dframe = data.frame(clusters=gene_partition_assignments, colors=partition_colors[gene_partition_assignments])
gene_colors = as.matrix(partition_colors[gene_partition_assignments])
pdf(paste('diffExpr.P1e-3_C1.matrix',height,'pdf',sep = '_'))
heatmap.3(data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, 
          RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, cexCol=1, margins=c(10,10))
dev.off()
gene_names = rownames(data)
num_cols = length(data[1,])
for (i in 1:max_cluster_count) {
  partition_i = (gene_partition_assignments == i)
  partition_data = data[partition_i,,drop=F]
  outfile = paste("subcluster_", i, "_log2_medianCentered_fpkm.matrix", sep='')
  write.table(partition_data, file=outfile, quote=F, sep="\t")
}
files = list.files(path=getwd(),pattern='fpkm.matrix')
pdf("my_cluster_plots.pdf")
par(mfrow=c(2,1))
par(cex=0.6)
par(mar=c(7,4,4,2))
for (i in 1:length(files)) {
  data = read.table(files[i], header=T, row.names=1)
  plot_label = paste(files[i], ', ', length(data[,1]), " trans", sep='')
  boxplot(data,col = unique(partition_colors)[i], pch = 4,lwd = 0.5,varwidth = T,outline = F,main=plot_label,par(las="2"))
}
dev.off()
```
**GO enrichment for DE analysis**
### GO enrichmet [script can be found here](https://github.com/caorui12/mouse-RNAseq/blob/main/GO_KEGG.R)
```
R --slave --args file_list.txt < GO&KEGG.R
```

## WCGNA see code [here](https://github.com/caorui12/mouse-RNAseq/blob/main/WGCNA.R)


### local similarity 
we use ELSA (https://bitbucket.org/charade/elsa/src/master)

 ```
lsa_compute 166gene.txt 166gene.out -s 12 -d 3
   
lsa <- read.table('166gene.out', sep = '\t', stringsAsFactors = FALSE)
 
#select Q <=0.001 
lsa <- subset(lsa, Q <= 0.001)
 
#write table for cytoscape visulization
write.table(lsa, 'ARISA20.lsa.select.txt', row.names = FALSE, sep = '\t', quote = FALSE)
```
###  propotionality (validate the pearson correlation)
refer to [propr](https://cran.r-project.org/web/packages/propr/index.html)
```
library(propr)
source("~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/R/heatmap.3.R")
source("~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/R/misc_rnaseq_funcs.R")
source("~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/R/pairs3.R")
source("~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/R/vioplot2.R")
setwd('~/Dropbox/mouseRNAseq/145N1N2remove/DESeq2sequential/DE_result/')
myheatcol = colorpanel(75, 'blue','white','red')
count_matrix<-read.table('../../mouse.coding.counts.matrix',sep='\t',row.names = 1,header=T)
DE<-read.table('../diffExpr.P1e-3_C1.matrix',sep='\t',header=T,row.names = 1)
DF<-count_matrix[rownames(DE),]
TMM<-read.table('../diffExpr.P1e-3_C1.matrix',sep='\t',row.names = 1,header=T)
TMM<- TMM[,-c(1,2)]
TMM<- log2(TMM+1)
data<-t(DF)
rho <- propr(data, metric = "rho", symmetrize = TRUE)
matrix<-rho@matrix
best <- rho[">", .995]
best <- simplify(best)
gene_dist<- as.dist(1-abs(rho@matrix))
hc_genes<-hclust(gene_dist)
plot(hc_genes)

gene_partition_assignments <- cutree(as.hclust(hc_genes),k=5)
unique(gene_partition_assignments )
max_cluster_count = max(gene_partition_assignments)
gene_names = rownames(TMM)
num_cols = length(TMM[1,])
for (i in 1:max_cluster_count) {
  partition_i = (gene_partition_assignments == i)
  partition_data = TMM[partition_i,,drop=F]
  outfile = paste("subcluster_", i, "_log2_medianCentered_fpkm.matrix", sep='')
  write.table(partition_data, file=outfile, quote=F, sep="\t")
}
files = list.files(path=getwd(),pattern='fpkm.matrix')
pdf("my_cluster_plots.pdf")
par(mfrow=c(2,1))
par(cex=0.6)
par(mar=c(7,4,4,2))
for (i in 1:length(files)) {
  data = read.table(files[i], header=T, row.names=1)
  plot_label = paste(files[i], ', ', length(data[,1]), " trans", sep='')
  boxplot(data,col = colors()[i], pch = 4,lwd = 0.5,varwidth = T,outline = F,main=plot_label)
}
dev.off()
 ```
