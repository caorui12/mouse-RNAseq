##step 1: data prepration
library(WGCNA)
library(reshape2)
library(stringr)
library(clusterProfiler)
options(stringsAsFactors = FALSE)
setwd('~/Dropbox/mouseRNAseq/sample_without_W/DESeq_seq/')
DF<-read.table('diffExpr.P1e-3_C1.matrix',header=T,row.names = 1)
datTraits<-data.frame(colnames(DF),c('E11.5','E11.5','E11.5','E12.5','E12.5','E12.5','E13.5','E13.5','E13.5','E14.5','E14.5','E14.5','E15.5','E15.5','E15.5',
                                     'E16.5','E16.5','E16.5','E17.5','E17.5','E17.5','E18.5','E18.5','E18.5',
                                     'P0','P0','P0','P1','P1','P1','P3','P3','P3','P7','P7','P7','P7'))
names(datTraits)[2]<-'sample_stage'
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
datExpr<-t(DF)

###  determine soft thresholding power
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
##  one step to construct co-expression matrix
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
###  modules and phenotypes
table(datTraits$sample_stage)
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  design=model.matrix(~0+ datTraits$sample_stage)
  colnames(design)=levels(as.factor(datTraits$sample_stage))## transfer to factor
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); 
  moduleTraitCor = cor(MEs, design , use = "p"); ## gene vs phenotype
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  sizeGrWindow(10,6)
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  png("Module-trait-relationships.png",width = 800,height = 1200,res = 120)
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
  stage = as.data.frame(design[,8]); ##select interest module and phenotype to visulize 
  names(stage) = "stageE18.5" 
  geneTraitSignificance = as.data.frame(cor(datExpr, stage, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(stage), sep="")
  names(GSPvalue) = paste("p.GS.", names(stage), sep="")  
  module = "yellow" ##select corresponding 
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
  
  
### identify hub gene for each module  
  
  
  datKME=signedKME(datExpr, MEs, outputColumnName="MM.") ###calculate connectivity
  FilterGenes= abs(geneTraitSignificance$GS.stageE18.5)> .2 & abs(datKME$MM.turquoise)>.8
  hub_list<-dimnames(data.frame(datExpr))[[2]][FilterGenes]
  hub_list ## identify hub gene

    
###  module export 
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
      threshold = 0.05,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule]
    )
