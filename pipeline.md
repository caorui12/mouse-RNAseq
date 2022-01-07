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
use trinity_script to generate genes-level experssion matrix
```
~/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix mouse_trans ../quantify/*.isoforms.results --gene_trans_map none
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

Here I keep two documents, one is DEs from pairwise comparison, one is only sequential comparison(two sequential comparison is to only retain the nearby stages, e.g. E10.5 vs E11.5, E11.5 vs E12.5 etc) 

**Sequential DE analysis** 
```
cd DESeq2.32273.dir(sequencial comparison)
~/Desktop/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../mouse_trans.isoform.TMM.EXPR.matrix --samples ../samples.txt -P 1e-3 -C 1 
```
total 3323 features identified 
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

## WCGNA 

```
##step 1: data prepration
library(WGCNA)
library(reshape2)
library(stringr)
library(clusterProfiler)
options(stringsAsFactors = FALSE)
setwd('~/Dropbox/mouseRNAseq/sample_without_W/DESeq_seq/')
DF<-read.table('diffExpr.P1e-3_C1.matrix',header=T,row.names = 1)
datTraits<-data.frame(colnames(DF),c('E11.5','E11.5','E11.5','E12.5','E12.5','E12.5','E13.5','E13.5','E13.5','E14.5','E14.5','E14.5','E14.5','E14.5','E15.5','E15.5','E15.5','E16.5','E16.5','E16.5','E17.5','E17.5','E17.5','E18.5','E18.5','E18.5','P0','P0','P0','P1','P1','P1','P3','P3','P3','P7','P7','P7','P7'))
names(datTraits)[2]<-'sample_stage'
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
datExpr<-t(DF)

### step 2: determine soft thresholding power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
## step 3: one step to construct co-expression matrix
net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F, 
  verbose = 3
)
table(net$colors) 
### step 4: module visulization
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## assign all of the gene to their corresponding module 
## hclust for the genes.
### step 5: modules and phenotypes
table(datTraits$sample_stage)
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  design=model.matrix(~0+ datTraits$sample_stage)
  colnames(design)=levels(as.factor(datTraits$sample_stage))##务必转换成因子
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  sizeGrWindow(10,6)
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  png("step5-Module-trait-relationships.png",width = 800,height = 1200,res = 120)
  par(mar = c(6, 8.5, 3, 3));
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
## module analysis
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  stage = as.data.frame(design[,8]);
  names(stage) = "stageE18.5"
  geneTraitSignificance = as.data.frame(cor(datExpr, stage, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(stage), sep="")
  names(GSPvalue) = paste("p.GS.", names(stage), sep="")  
  module = "yellow"
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for stage E12.5",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
  
## network visulization  
  if(T){
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    geneTree = net$dendrograms[[1]]; 
    dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); 
    plotTOM = dissTOM^sft$powerEstimate; 
    diag(plotTOM) = NA; 
    #TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
    nSelect = 3378
    # For reproducibility, we set the random seed
    set.seed(10);
    select = sample(nGenes, size = nSelect);
    selectTOM = dissTOM[select, select];
    # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
    selectTree = hclust(as.dist(selectTOM), method = "average")
    selectColors = moduleColors[select];
    # Open a graphical window
    sizeGrWindow(9,9)
    # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
    # the color palette; setting the diagonal to NA also improves the clarity of the plot
    plotDiss = selectTOM^7;
    diag(plotDiss) = NA;
    
    png("step7-Network-heatmap.png",width = 800,height = 600)
    TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
    dev.off()
    
    # Recalculate module eigengenes
    MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

 
    P7 = as.data.frame(design[,12]);
    names(P7) = "P7"
    # Add the weight to existing module eigengenes
    MET = orderMEs(cbind(MEs, P7))
    # Plot the relationships among the eigengenes and the trait
    sizeGrWindow(5,7.5);
    
    par(cex = 0.9)
    png("step7-Eigengene-dendrogram.png",width = 800,height = 600)
    plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                          plotDendrograms = FALSE, xLabelsAngle = 90)
    dev.off()
    
    # Plot the dendrogram
    sizeGrWindow(6,6);
    par(cex = 1.0)
   
    png("step7-Eigengene-dendrogram-hclust.png",width = 800,height = 600)
    plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                          plotHeatmaps = FALSE)
    dev.off()
    # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
    par(cex = 1.0)
    
    
    png("step7-Eigengene-adjacency-heatmap.png",width = 800,height = 600)
    plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                          plotDendrograms = FALSE, xLabelsAngle = 90)
    dev.off()
    
  }
  
  datKME=signedKME(datExpr, MEs, outputColumnName="MM.") ###calculate connectivity
  FilterGenes= abs(geneTraitSignificance$GS.stageP7)> .2 & abs(datKME$MM.turquoise)>.8
  hub_list<-dimnames(data.frame(datExpr))[[2]][FilterGenes]
  hub_list ## identify hub gene
### step 8
  
    # Select module
    module = "blue"
    # Select module probes
    probes = colnames(datExpr) 
    inModule = (moduleColors==module);
    modProbes = probes[inModule]; 
    head(modProbes)
    
  
    library(pheatmap)
    dat=datExpr[,moduleColors==which.module ]
    pheatmap(dat ,show_colnames =F,show_rownames = F) 
    n=t(scale(t(log(dat+1)))) # 'scale'可以对log-ratio数值进行归一化
    n[n>2]=2 
    n[n< -2]= -2
    n[1:4,1:4]
    pheatmap(n,show_colnames =F,show_rownames = F)
    group_list=datTraits$sample_stage
    ac=data.frame(g=group_list)
    rownames(ac)=colnames(n) 
    pheatmap(n,show_colnames =F,show_rownames = F,
             annotation_col=ac )
    
### step 9 module export 
    # Recalculate topological overlap
    TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate); 
    # Select module
    module = "yellow";
    # Select module probes
    probes = colnames(datExpr) 
    inModule = (moduleColors==module);
    modProbes = probes[inModule]; 
   
    # Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule];
    dimnames(modTOM) = list(modProbes, modProbes)
   
    cyt = exportNetworkToCytoscape(
      modTOM,
      edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
      nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.02,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule]
    )
  ```
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
DF<-DF[,-c(1,2)]
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
