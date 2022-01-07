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