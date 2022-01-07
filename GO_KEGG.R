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
goCC <- enrichGO(rownames(target),OrgDb = org.Mm.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'SYMBOL')
goBP <- enrichGO(rownames(target),OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'SYMBOL')
goMF <- enrichGO(rownames(target),OrgDb = org.Mm.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'SYMBOL')
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
kegg_list<-bitr(rownames(target),fromType = 'SYMBOL', toType = 'ENTREZID',OrgDb = org.Mm.eg.db)
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
