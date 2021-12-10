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
### GO enrichmet
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(topGO)
library(GSEABase)
library(GO.db)
library(ggplot2)
library(stringr)
input_file=commandArgs()[4]
files<-read.csv(input_file,header=F)


for (i in 1:nrow(files)){
print(paste('we are processing',files[i,],sep=' '))
target<-read.table(files[i,], header=T, row.names=1,sep='\t')
goCC <- enrichGO(rownames(target),OrgDb = org.Mm.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'ENSEMBL')
goBP <- enrichGO(rownames(target),OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'ENSEMBL')
goMF <- enrichGO(rownames(target),OrgDb = org.Mm.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'ENSEMBL')
ego_result_BP <- as.data.frame(goBP)
ego_result_CC <- as.data.frame(goCC)
ego_result_MF <- as.data.frame(goMF)
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)
write.table(ego, paste(files[i,],'_','GO.txt',sep=''), sep='\t', row.names = FALSE, quote = FALSE) 

display_number = 10 #这个数字分别代表选取的BP、CC、MF的通路条数
ego_result_BP <- as.data.frame(goBP)[1:display_number, ]
ego_result_CC <- as.data.frame(goCC)[1:display_number, ]
ego_result_MF <- as.data.frame(goMF)[1:display_number, ]

go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number), 
                rep("cellular component", display_number),
                rep("molecular function", display_number)), 
              levels=c("biological process", "cellular component","molecular function" )))

go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色
pdf(paste(files[i,],'_','GO.pdf',sep=''))
p<-ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") + 
  ylab("Gene_Number") +
  theme_bw()

print(p + scale_x_discrete(labels=function(x) str_wrap(x, width=50)))
dev.off()

### start KEGG 
kegg_list<-bitr(rownames(target),fromType = 'ENSEMBL', toType = 'ENTREZID',OrgDb = org.Mm.eg.db)
gene=kegg_list$ENTREZID
pdf(paste(files[i,],'_','KEGG.pdf',sep=''))
kk <- enrichKEGG(
  gene = gene,
  organism  = 'mmu',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05
)
k<-barplot(kk,showCategory=10)
print(k + scale_y_discrete(labels=function(x) str_wrap(x, width=50)))
dev.off()
write.table(kk, paste(files[i,],'_','KEGG.txt',sep=''), sep='\t', row.names = FALSE, quote = FALSE) 
print(paste(files[i,],'finished',sep=' '))
}

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
  ## 算出每个模块跟基因的皮尔森相关系数矩阵
  ## MEs是每个模块在每个样本里面的值
  ## datExpr是每个基因在每个样本的表达量
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
    ## 只有连续型性状才能只有计算
    ## 这里把是否属 P7 表型这个变量0,1进行数值化
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
    ## 模块的进化树
    png("step7-Eigengene-dendrogram-hclust.png",width = 800,height = 600)
    plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                          plotHeatmaps = FALSE)
    dev.off()
    # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
    par(cex = 1.0)
    ## 性状与模块热
    
    png("step7-Eigengene-adjacency-heatmap.png",width = 800,height = 600)
    plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                          plotDendrograms = FALSE, xLabelsAngle = 90)
    dev.off()
    
  }
  
  datKME=signedKME(datExpr, MEs, outputColumnName="MM.") ###计算connectivity
  FilterGenes= abs(geneTraitSignificance$GS.stageP7)> .2 & abs(datKME$MM.turquoise)>.8
  hub_list<-dimnames(data.frame(datExpr))[[2]][FilterGenes]
  hub_list
### step 8
  
    # Select module
    module = "blue"
    # Select module probes
    probes = colnames(datExpr) ## 我们例子里面的probe就是基因
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
    probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
    inModule = (moduleColors==module);
    modProbes = probes[inModule]; 
    ## 也是提取指定模块的基因名
    # Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule];
    dimnames(modTOM) = list(modProbes, modProbes)
    ## 模块对应的基因关系矩阵 
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


