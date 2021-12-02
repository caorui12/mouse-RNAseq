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
```
R --slave --args file_list.txt < GO&KEGG.R

### local similarity analysis
```
library(rlist)
rankNormalization <- function(x)
{
  xLength <- length(x)
  rankScore <- rep(0,xLength)
  rankScore <- scale(qnorm(rank(x)/(xLength+1)))
  return(rankScore)
}

##################################################################
# LocalSimilarity<- function(x, y, maxDelay=3, rankscale=FALSE)
#
#  This function computes the local similarity score for two sequences.
#
# INPUT:
# ======
#
#	x, y	: sequences to copute LS score
#	maxDelay	: maximum time shift allowed in computing LS score.
#	rankscale		: If TRUE, perform rankNormalization first; False, otherwise.
#
# RETURN:
# ======
#
#  A six element vector contains: c(scoreMax, startX, startY, delay, length, PosOrNeg)
#

LocalSimilarity <- function(x, y, maxDelay=3, rankScale = FALSE){
  
  if (rankScale == TRUE)
  {
    x <- rankNormalization(x)
    y <- rankNormalization(y)
  }
  
  scoreMax <- 0;
  PosOrNeg <- 0;
  startX <- 0;
  startY <- 0;
  numTimepoints <- length(x)
  posMatrix <- matrix(0, numTimepoints+1, numTimepoints+1)
  negMatrix <- matrix(0, numTimepoints+1, numTimepoints+1)
  for (i in 1:numTimepoints){
    for (j in max(1, i-maxDelay):min(i+maxDelay, numTimepoints)){
      posMatrix[i+1,j+1] <- max(0, posMatrix[i,j] + x[i] * y[j])
      negMatrix[i+1,j+1] <- max(0, negMatrix[i,j] - x[i] * y[j])
    }
  }
  posMatrix <- posMatrix[-1,]
  posMatrix <- posMatrix[,-1]
  negMatrix <- negMatrix[-1,]
  negMatrix <- negMatrix[,-1]
  scorePosmax <- max(posMatrix)
  scoreNegmax <- max(negMatrix)
  scoreMax <- max(scorePosmax, scoreNegmax)
  if( scorePosmax > scoreNegmax) PosOrNeg = 1 else PosOrNeg = -1
  if (scorePosmax > scoreNegmax) {
    Maxposition <- which(posMatrix == scorePosmax, arr.ind = TRUE)
  }else  Maxposition <- which(negMatrix == scoreNegmax, arr.ind = TRUE)
  delay <- Maxposition[1] - Maxposition[2]
  subinternal<-NULL
  for(i in max(1, delay):min(numTimepoints,numTimepoints + delay)){
    if (scorePosmax > scoreNegmax) {
      subinternal <- c(subinternal,posMatrix[i, i-delay])
    }else  subinternal <- c(subinternal,negMatrix[i,i-delay])
  }
  if (delay>0) {
    startX <- max(1,which(subinternal[1:Maxposition[1]] == 0) + 1) + delay
  }else{
    startX <- max(1,which(subinternal[1:Maxposition[1]] == 0) + 1)
  }
  startY <- startX - delay
  lengths <- Maxposition[1] - startX + 1
  value <- t(c(scoreMax, startX, startY, delay, lengths, PosOrNeg))
  colnames(value)<-c('Maxscore', 'startX', 'startY', 'delay', 'length', 'PosOrNeg')
  return(value)
}
R --slave --args GO_list.txt < GO&KEGG.R
######################
# permutationTest(x, y, numPermu=1000, maxDelay=3,scale=TRUE)
#
# This function computes the significance level of the LS score by permutation test
#
######################
permutationTest <- function(x, y, numPermu=1000, maxDelay=3,scale=TRUE){
  scoreArray <- rep(0.0, numPermu+1);
  
  if (scale == TRUE){
    x <- rankNormalization(x)
    y <- rankNormalization(y)
  }
  numTimePoints <- length(x)
  scoreMax1 <- LocalSimilarity(x, y, maxDelay)[1];
  scoreArray[1] <- scoreMax1;
  highScoreCnt <- 0;
  
  for(idx in 1:numPermu)
  {
    dtt1 <- x[sample(numTimePoints)]
    scoreTmp <- LocalSimilarity(dtt1, y, maxDelay)[1];
    
    scoreArray[idx+1] <- scoreTmp;
    highScoreCnt <- highScoreCnt + (scoreTmp >= scoreMax1);
  }
  
  pValue <- 1.0 * highScoreCnt / numPermu;
  
  return(pValue);
}

### main analysis
main<-function(DF,maxDelay){
  DF<-read.table(DF)
  for (i in 1:nrow(DF)){                      ## pairwise comparison 
    for (j in 1:nrow(DF)){
      if (i<j){
        X<-as.vector(DF[i,])
        Y<-as.vector(DF[j,])
        local_score<-LocalSimilarity(x=X, y=Y, maxDelay=maxDelay, rankScale = TRUE)
        P_value<-permutationTest(x=X, y=Y, numPermu=1000,maxDelay=maxDelay, scale = TRUE)
        result<-list.append (local_score, P_value) ## append to the result list
        LSA_score<-rbind(result)     ## output final score for each comparison
      } 
      else {
        next
      }
    }
  }
  return(LSA_score)
}
input_file=commandArgs()[4]
maxDelay=commandArgs()[5]
maxDelay=as.integer(maxDelay)
LSA_matrix<-main(DF = input_file,maxDelay = maxDelay)   
write.table(LSA_matrix,'LSA_matrix',sep='\t',quote=F)
```
R --slave --args 300meanTMM_DE_LSA 3 < LSA_compute.R
