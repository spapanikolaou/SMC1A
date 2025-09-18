### Figures for Panels
library(dplyr)
library(ggplot2)
library(ggrepel)
#install.packages("ggvenn")
library(ggvenn)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)

# ---------- Paths (edit these) ----------
setwd("~/SofiaP/smc1a/")
data <- read.delim("final_df_new_all_rnaseq_comp_with_active_enh.txt", header=T, sep="\t")
colnames(data)[1]<-"Gene"
### save s6 table as xlsx - 10/09/2025 (log2FC values of lupus-responsive genes) #2,778 
lrg<-data[which(abs(data$log2FoldChange_ITLvsC_siC)>1&data$padjvalue_ITLvsC_siC<0.05),]
lrg_data<-lrg[-(which(duplicated(lrg$Gene))),c(1,4,5)]
#write.xlsx(lrg_data, file = "./final_fig/SupplTable6_lupus-responsive-genes.xlsx")


data$active_prom <- as.character(data$active_prom)
data$active_prom[is.na(data$active_prom)]<-"No"
data$active_prom <- as.factor(data$active_prom)
#
data$active_enh <- as.character(data$active_enh)
data$active_enh[is.na(data$active_enh)]<-"No"
data$active_enh <- as.factor(data$active_enh)
#
# Creating a grouping of Enhancer Status
## Using this method, if a gene has both an up-bound enhancer and a down-bound enhancer, it is reported as having a down-bound enhancer.
data$EnhancerStatus<-"None"
data$EnhancerStatus[which(data$active_enh=="Yes")]<-"Active"
data$EnhancerStatus[which(data$up_enhancers=="Yes")]<-"Up"
data$EnhancerStatus[which(data$up_active_enhancers=="Yes")]<-"UpActive"
data$EnhancerStatus[which(data$down_enhancers=="Yes")]<-"Down"
data$EnhancerStatus[which(data$down_active_enhancers=="Yes")]<-"DownActive"
data$EnhancerStatus<-factor(data$EnhancerStatus, levels=c("None", "Active", "Down", "Up", "DownActive", "UpActive"))
#
data$PromoterStatus<-"None"
data$PromoterStatus[which(data$active_prom=="Yes")]<-"Active"
data$PromoterStatus[which(data$up_promoters=="Yes")]<-"Up"
data$PromoterStatus[which(data$up_active_promoters=="Yes")]<-"UpActive"
data$PromoterStatus[which(data$down_promoters=="Yes")]<-"Down"
data$PromoterStatus[which(data$down_active_promoters=="Yes")]<-"DownActive"
data$PromoterStatus<-factor(data$PromoterStatus, levels=c("None", "Active", "Down", "Up", "DownActive", "UpActive"))

siUntrgenes<-arrange(data[which( data$pvalue_siS.vs.siC_untr<=0.05) , c(1,8,9)], log2FoldChange_siS.vs.siC_untr)
colnames(siUntrgenes)<-c("Gene", "log2FC", "pvalue")
siUntrgenes$expr[siUntrgenes$log2FC>0]<-"up"
siUntrgenes$expr[siUntrgenes$log2FC<0]<-"down"
siUntrgenes$expr<-as.factor(siUntrgenes$expr)

#
siUntrgenes$delabel <- NA
siUntrgenes$delabel[which(abs(siUntrgenes$log2FC)>=0.5)]<-as.character(siUntrgenes$Gene[which(abs(siUntrgenes$log2FC)>=0.5)])
#

#Figure 2B
g <- ggplot(siUntrgenes, aes(x = log2FC, y = -log10(pvalue), color = expr, label = delabel)) + 
     geom_point() + 
     scale_color_manual(values = c("steelblue4", "firebrick4")) +
     geom_text_repel(box.padding = 0.3, col = "black", size = 8 / .pt) + # Adjust size for publication
    # ggtitle(label = "si-SMC1A Untreated Monocytes", subtitle = "Gene Expression Changes (442 Genes at p<=0.05)") +
     xlab("log2FC") + 
     ylab("-log10(p-value)") + 
     theme_minimal() +
     theme(
        #plot.title = element_text(color = "black", size = 8, hjust = 0.5),
           #plot.subtitle = element_text(color = "black", size = 8, hjust = 0.5),
           axis.title.x = element_text(color = "black", size = 8),
         axis.title.y = element_text(color = "black", size = 8),
         axis.text.x = element_text(color = "black", size = 8),
         axis.text.y = element_text(color = "black", size = 8),
         panel.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "black")
       )
#ggsave("final_fig/Figure2B_Volcano - siSCMC1A untreated monocytes.png",g,width = 88, height = 88,units ="mm" ,dpi = 300)
#
#

sitrgenes<-arrange(data[which( data$pvalue_siS.vs.siC_ITL<=0.05) , c(1,6,7)], log2FoldChange_siS.vs.siC_ITL)
colnames(sitrgenes)<-c("Gene", "log2FC", "pvalue")
sitrgenes$expr[sitrgenes$log2FC>0]<-"up"
sitrgenes$expr[sitrgenes$log2FC<0]<-"down"
sitrgenes$expr<-as.factor(sitrgenes$expr)
#
sitrgenes$delabel <- NA
sitrgenes$delabel[which(abs(sitrgenes$log2FC)>=0.5)]<-as.character(sitrgenes$Gene[which(abs(sitrgenes$log2FC)>=0.5)])

##Figure 2C
g<-ggplot(sitrgenes, aes(x=log2FC, y=-log10(pvalue), color=expr, label=delabel)) + 
  geom_point() + 
  scale_color_manual(values=c("steelblue4", "firebrick4")) +
  geom_text_repel(box.padding = 0.3, col="black", size=8/.pt) +
#  ggtitle(label="si-SMC1A Treated Monocytes", subtitle ="Gene Expression Changes (418 Genes at p<=0.05)") +
  xlab("log2FC") + 
  ylab("-log10(p-value)") + 
  theme(
    #plot.title = element_text(color="black", size=19),
    #plot.subtitle = element_text(color="black", size=15),
    axis.title.x = element_text(color="black", size=8),
    axis.title.y = element_text(color="black", size=8),
    axis.text.x = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.background = element_rect(fill="gray97"))
#ggsave("final_fig/Figure2C_Volcano - siSCMC1A treated monocytes.png",g,width = 88, height = 88,units ="mm" ,dpi = 300)


##Figure 2D
treatedgenes <- arrange(data[which(data$pvalue_siS.vs.siC_ITL<=0.05), c(1,6,7)], log2FoldChange_siS.vs.siC_ITL)
untreatedgenes <- arrange(data[which(data$pvalue_siS.vs.siC_untr<=0.05), c(1,8,9)], log2FoldChange_siS.vs.siC_untr)
x <- list(Primed=untreatedgenes$Gene, Colonized=treatedgenes$Gene)
par(mar=c(5,5,10,5))
ggvenn(x, fill_color = c("#0073C2FF",  "#CD534CFF"), stroke_size = 2, text_size =6, set_name_size = 7)

allgenes<-union(treatedgenes$Gene, untreatedgenes$Gene)
data %>% filter(Gene %in% allgenes) %>% dplyr::select(c(1,6,8)) -> seldata
colnames(seldata)<-c("Gene", "logFCTr", "logFCUn")
#
commonsiSMC1A<-intersect(treatedgenes$Gene, untreatedgenes$Gene)
treatedSpecific<-setdiff(treatedgenes$Gene, untreatedgenes$Gene)
untreatedSpecific<-setdiff(untreatedgenes$Gene, treatedgenes$Gene)

seldata$type <- NA
seldata$type[which(seldata$Gene %in% commonsiSMC1A)] <- "Common"
seldata$type[which(seldata$Gene %in% setdiff(treatedgenes$Gene, untreatedgenes$Gene))] <- "TreatedSpecific"
seldata$type[which(seldata$Gene %in% setdiff(untreatedgenes$Gene,treatedgenes$Gene))] <- "UntreatedSpecific"
#
seldata$delabel <- NA
#seldata$delabel[which(abs(seldata$logFCTr)>=0.58 | abs(seldata$logFCUn)>=0.58)]<-as.character(seldata$Gene[which(abs(seldata$logFCTr)>=0.58 | abs(seldata$logFCUn)>=0.58)])
seldata$delabel[which(seldata$type=="Common")]<-as.character(seldata$Gene[which(seldata$type=="Common")])
#

library(RColorBrewer)
brewer_colors <- brewer.pal(n=3, name="Set2")

g<-ggplot(seldata, aes(x=logFCTr, y=logFCUn, label=delabel, col=type)) + 
  geom_point() + 
  scale_color_manual(values=c(c("#E69F00", "#56B4E9", "#009E73"))) +
 geom_text_repel(box.padding = 0.3, max.overlaps = Inf, col="black", size=8/.pt) +
  labs(x="Treated log2FC", y="Untreated log2FC")+
theme(
  axis.title.x = element_text(color="black", size=8),
  axis.title.y = element_text(color="black", size=8),
  axis.text.x = element_text(color="black", size=8),
  axis.text.y = element_text(color="black", size=8),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),  
  legend.background = element_rect(fill="gray97"))
#ggsave("final_fig/Figure2D_Gene Expression Changes in SMC1A Down-regulation.png",g,width = 130, height = 88,units ="mm" ,dpi = 300)

## for revisions
fig2e_data<-seldata
#write.table(fig2e_data,"~/SofiaP/smc1a/revisions/fig2e_data.txt",quote = F, row.names = F, col.names = T)

##
#

#Figure 2E
### compare cluster
primed_genes<-siUntrgenes$Gene
colonized_genes<-sitrgenes$Gene

enrich_gene_list_f<-list(Untreated=as.vector(na.omit(primed_genes)),Treated=as.vector(na.omit(colonized_genes)))
enrich_go_f <- compareCluster(geneCluster = enrich_gene_list_f,
                              ont           = "ALL",
                              OrgDb = org.Hs.eg.db,
                              keyType       = 'SYMBOL',
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              fun =  "enrichGO")
out_SIMPL<-simplify(enrich_go_f,cutoff = 0.7,semData=NULL)                     

d_a1 <- dotplot(out_SIMPL,font.size=8/.pt) +
  theme(
    #plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
ggsave("final_fig/Figure2E_GO enrichment analysis of SMC1a-Silenced Untreated and Treated DEGs.png",d_a1,width = 120, height = 100,units ="mm" ,dpi = 300)

# Figure 3
## diffbind output (all regions) smc1a
df_smc1a_volc<-read.table("~/Downloads/diffbind.deseq2.unique.all.no.untr.3 .bed",header = T)
df_h3k27ac_volc<-read.table("~/Downloads/diffbind.deseq2.h3k27ac.bed",header = T)

volcano.pval.ggplot <- function(data, pval_cut, logFC_cut, height_lim){
  # Define the data subsets
  data$regulation <- "no"  # Default to 'no' regulation
  data$regulation[data$p.value <= pval_cut & data$Fold >= logFC_cut] <- "up"
  data$regulation[data$p.value <= pval_cut & data$Fold <= -logFC_cut] <- "down"
  
  # Set the color mapping for the regulation status
  data$regulation <- factor(data$regulation, levels = c("no", "up", "down"))
  
  # Plot using ggplot2
  p <- ggplot(data, aes(x = Fold, y = -log10(p.value), color = regulation)) +
    geom_point(size = 3) +  # Points for the volcano plot
    scale_color_manual(values = c("no" = "gray40", "up" = "firebrick4", "down" = "steelblue4")) +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "black", size = 1) +  # p-value cutoff
    geom_vline(xintercept = c(-logFC_cut, logFC_cut), linetype = "dotted", color = "black", size = 1) +  # Fold change cutoff
    labs(
      title = "SMC1a binding",
      x = "Fold change",
      y = "-log10(p-value)"
    ) +
    xlim(min(data$Fold) - 0.5, max(data$Fold) + 0.5) +  # Set x limits
    ylim(0, height_lim) +  # Set y limits
    theme(
      plot.title = element_text(color="black", size=19),
      plot.subtitle = element_text(color="black", size=15),
      axis.title.x = element_text(color="black", size=15),
      axis.title.y = element_text(color="black", size=15),
      axis.text.x = element_text(color="black", size=12),
      axis.text.y = element_text(color="black", size=12),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12))  # Optional: adjust legend appearance
  
  # Print the plot
  print(p)
}

png("final_fig/Figure3A_Volcano - ChIP SMC1a.png",width = 88, height = 88,units ="mm" ,res = 300)
volcano.pval.ggplot(df_smc1a_volc,0.001,1,20)
# Close the device
dev.off()

######################
volcano.pval.ggplot <- function(data, pval_cut, logFC_cut, height_lim){
  # Define the data subsets
  data$regulation <- "no"  # Default to 'no' regulation
  data$regulation[data$p.value <= pval_cut & data$Fold >= logFC_cut] <- "up"
  data$regulation[data$p.value <= pval_cut & data$Fold <= -logFC_cut] <- "down"
  
  # Set the color mapping for the regulation status
  data$regulation <- factor(data$regulation, levels = c("no", "up", "down"))
  
  # Plot using ggplot2
  p <- ggplot(data, aes(x = Fold, y = -log10(p.value), color = regulation)) +
    geom_point(size = 3) +  # Points for the volcano plot
    scale_color_manual(values = c("no" = "gray40", "up" = "firebrick4", "down" = "steelblue4")) +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "black", size = 1) +  # p-value cutoff
    geom_vline(xintercept = c(-logFC_cut, logFC_cut), linetype = "dotted", color = "black", size = 1) +  # Fold change cutoff
    labs(
      title = "H3k27ac binding",
      x = "Fold change",
      y = "-log10(p-value)"
    ) +
    xlim(min(data$Fold) - 0.5, max(data$Fold) + 0.5) +  # Set x limits
    ylim(0, height_lim) +  # Set y limits
    theme(
      plot.title = element_text(color="black", size=19),
      plot.subtitle = element_text(color="black", size=15),
      axis.title.x = element_text(color="black", size=15),
      axis.title.y = element_text(color="black", size=15),
      axis.text.x = element_text(color="black", size=12),
      axis.text.y = element_text(color="black", size=12),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12))  # Optional: adjust legend appearance
  
  # Print the plot
  print(p)
}

png("final_fig/Figure3B_Volcano - ChIP H3k27ac.png", 
    width = 88, height = 88,units ="mm" ,res = 300)
volcano.pval.ggplot(df_h3k27ac_volc,0.001,1,20)
# Close the device
dev.off()


#Figure 3C
###smc1a
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
samplefiles <- list.files(path="~/SofiaP/smc1a/",pattern= "diffbind.deseq2.", full.names=T)
samplefiles <- as.list(samplefiles)
samplefiles<-samplefiles[c(2,3)]
names(samplefiles) <- c("smc1a_down_bound","smc1a_up_bound")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)
p<-plotAnnoBar(peakAnnoList)+theme(
  #plot.title = element_text(color = "black", size = 19),
  axis.title.x = element_text(color = "black", size = 8),
  axis.title.y = element_text(color = "black", size = 8),
  axis.text.x = element_text(color = "black", size = 8),
  axis.text.y = element_text(color = "black", size = 8),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black")
)
ggsave("final_fig/Figure3C_peak_annotation_smc1a.png",p,width = 140, height = 88,units ="mm", dpi = 300)


#Figure 3D-E
#Figure S5G
smc1a_up_bound_annot<-read.delim("peaks_annot_up_smc1a_from_diffbind_output_pval001_fc1.txt",header = T)
smc1a.up<-unique(smc1a_up_bound_annot$SYMBOL) #4839
smc1a_down_bound_annot<-read.delim("peaks_annot_down_smc1a_from_diffbind_output_pval001_fc1.txt",header = T)
smc1a.down<-unique(smc1a_down_bound_annot$SYMBOL) #1023

smc1a_down_bound_PROMOTERS<-smc1a_down_bound_annot[which(smc1a_down_bound_annot$annotation=="Promoter"),"SYMBOL"] #unique 823
ego_smc1a_down_bound_PROMOTERS <- enrichGO(gene = unique(smc1a_down_bound_PROMOTERS), #
                           keyType = "SYMBOL", 
                           OrgDb = org.Hs.eg.db, 
                           ont = "ALL", 
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.05, 
                           readable = TRUE)
out_SIMPL<-simplify(ego_smc1a_down_bound_PROMOTERS,cutoff = 0.7,semData=NULL)                     
d_a1 <- dotplot(out_SIMPL,showCategory=15,orderBy = "x")+
  theme(
    #plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
ggsave("final_fig/SupplFigure5G_Enriched terms of smc1a down bound PROMOTERS_from_diffbind_output_pval001_fc1.png",d_a1,width = 140, height = 130,units ="mm" ,dpi = 300)

### ENHANCERS
enhancers_up_smc1a<-read.table("./final_fig/genes_associated_with_up_enh_onlySMC1a_NEW.bed")
ego_smc1a_up_bound_ENHANCERS <- enrichGO(gene = unique(enhancers_up_smc1a$V1), #
                                           keyType = "SYMBOL", 
                                           OrgDb = org.Hs.eg.db, 
                                           ont = "ALL", 
                                           pAdjustMethod = "BH", 
                                           pvalueCutoff = 0.05, 
                                           readable = TRUE)
out_SIMPL<-simplify(ego_smc1a_up_bound_ENHANCERS,cutoff = 0.7,semData=NULL)                     
d_a1 <- dotplot(out_SIMPL,showCategory=15,orderBy = "x")+
  theme(
    #plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
ggsave("final_fig/SupplFigure5e_Enriched terms of smc1a up bound ENHANCERS_from_intersection_GENEHANCER_and_diffbind_output_pval001_fc1.png",d_a1,width = 140, height = 130,units ="mm" ,dpi = 300)



h3k27ac_up_bound_annot<-read.delim("peaks_annot_up_h3k27ac_from_diffbind_output_pval001_fc1.txt",header = T)
h3k27ac.up<-unique(h3k27ac_up_bound_annot$SYMBOL) #440
ego_h3k27ac_up <- enrichGO(gene = h3k27ac.up, #
                           keyType = "SYMBOL", 
                           OrgDb = org.Hs.eg.db, 
                           ont = "ALL", 
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.05, 
                           readable = TRUE)
out_SIMPL<-simplify(ego_h3k27ac_up,cutoff = 0.7,semData=NULL)                     
d_a1 <- dotplot(out_SIMPL,showCategory=15,orderBy = "x")+#(), title = "Enriched terms of \nH3k27ac up bound genes") +
  theme(
    #plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
#ggsave("./Enriched terms of h3k27ac up bound genes.png",d_a1,width = 6, height = 7, dpi = 300)
  ggsave("final_fig/SupplFigure5G_Enriched terms of h3k27ac up bound genes.png",d_a1,width = 140, height = 130,units ="mm" ,dpi = 300)
  
h3k27ac_down_bound_annot<-read.delim("peaks_annot_down_h3k27ac_from_diffbind_output_pval001_fc1.txt",header = T)
h3k27ac.down<-unique(h3k27ac_down_bound_annot$SYMBOL) #269
  
upup<-intersect(smc1a.up,h3k27ac.up) #203
downdown<-intersect(smc1a.down,h3k27ac.down) #72
upup<-setdiff(upup,downdown)
upup<-setdiff(upup,smc1a.down)
upup<-setdiff(upup,h3k27ac.down)
ego_upup <-enrichGO(gene = upup,
                    keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05, readable = FALSE)
ego_upup_SIMPL<-simplify(ego_upup,cutoff = 0.7,semData=NULL)
d<-dotplot(ego_upup_SIMPL, showCategory=15)+#, title = "Enriched terms among genes commonly \nup bound by SMC1A and H3K27ac")+
  theme(
    #plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))

ggsave("./final_fig/Figure3E_Enriched terms among genes commonly differentially upbound by SMC1a and H3K27ac.png",d,width = 110, height = 120,units ="mm" ,dpi = 300)

ego_downdown <- enrichGO(gene = downdown, #
                         keyType = "SYMBOL", 
                         OrgDb = org.Hs.eg.db, 
                         ont = "ALL", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)

#Figure 3D
############################ VENN DIAGRAM
# Calculate the values for the Venn diagram
# Load the dplyr package
#BiocManager::install("eulerr")
library(eulerr)
library(dplyr)
h3k27ac_bound_annot<-rbind(h3k27ac_up_bound_annot,h3k27ac_down_bound_annot)
colnames(h3k27ac_bound_annot)[23]<-"gene_name"
h3k27ac_bound_annot <- h3k27ac_bound_annot %>%
  group_by(gene_name) %>%                # Group by the 'SYMBOL' column
  slice_min(order_by = abs(distanceToTSS)) %>%  # Keep the row with the max 'distanceToTSS' in each group
  ungroup()  
smc1a_bound_annot<-rbind(smc1a_up_bound_annot,smc1a_down_bound_annot)
colnames(smc1a_bound_annot)[23]<-"gene_name"
smc1a_bound_annot <- smc1a_bound_annot %>%
  group_by(gene_name) %>%                # Group by the 'SYMBOL' column
  slice_min(order_by = abs(distanceToTSS)) %>%  # Keep the row with the max 'distanceToTSS' in each group
  ungroup()   
smc1a_up<-smc1a_bound_annot[which(smc1a_bound_annot$V9>0),]
smc1a_down<-smc1a_bound_annot[which(smc1a_bound_annot$V9<0),]
h3k27ac_up<-h3k27ac_bound_annot[which(h3k27ac_bound_annot$V9>0),]
h3k27ac_down<-h3k27ac_bound_annot[which(h3k27ac_bound_annot$V9<0),]
i_s<-intersect(smc1a_up$gene_name,smc1a_down$gene_name) #3
intersect(h3k27ac_up$gene_name,h3k27ac_down$gene_name)

# exclude common genes bound on promoters
smc1a_up<-smc1a_up[-(which(smc1a_up$gene_name%in%i_s)),]
smc1a_down<-smc1a_down[-(which(smc1a_down$gene_name%in%i_s)),]


library("eulerr")


# Calculate the values for the Euler diagram
# Print the results of each step
print("Step 1: SMC1a_up")
SMC1a_up_value <- length(smc1a_up$gene_name) - 
  length(intersect(smc1a_up$gene_name, h3k27ac_up$gene_name)) - 
  length(intersect(smc1a_up$gene_name, h3k27ac_down$gene_name)) + 
  length(Reduce(intersect, list(smc1a_up$gene_name, smc1a_down$gene_name, h3k27ac_up$gene_name, h3k27ac_down$gene_name)))
print(SMC1a_up_value)

print("Step 2: SMC1a_down")
SMC1a_down_value <- length(smc1a_down$gene_name) - 
  length(intersect(smc1a_down$gene_name, h3k27ac_up$gene_name)) - 
  length(intersect(smc1a_down$gene_name, h3k27ac_down$gene_name)) + 
  length(Reduce(intersect, list(smc1a_up$gene_name, smc1a_down$gene_name, h3k27ac_up$gene_name, h3k27ac_down$gene_name)))
print(SMC1a_down_value)

print("Step 3: H3k27ac_up")
H3k27ac_up_value <- length(h3k27ac_up$gene_name) - 
  length(intersect(h3k27ac_up$gene_name, smc1a_up$gene_name)) - 
  length(intersect(h3k27ac_up$gene_name, smc1a_down$gene_name)) - 
  length(Reduce(intersect, list(smc1a_up$gene_name, smc1a_down$gene_name, h3k27ac_up$gene_name, h3k27ac_down$gene_name)))
print(H3k27ac_up_value)

print("Step 4: H3k27ac_down")
H3k27ac_down_value <- length(h3k27ac_down$gene_name) - 
  length(intersect(h3k27ac_down$gene_name, smc1a_up$gene_name)) - 
  length(intersect(h3k27ac_down$gene_name, smc1a_down$gene_name)) - 
  length(Reduce(intersect, list(smc1a_up$gene_name, smc1a_down$gene_name, h3k27ac_up$gene_name, h3k27ac_down$gene_name)))
print(H3k27ac_down_value)

print("Step 5: SMC1a_up & H3k27ac_up")
SMC1a_up_H3k27ac_up_value <- length(intersect(smc1a_up$gene_name, h3k27ac_up$gene_name)) - 
  length(Reduce(intersect, list(smc1a_up$gene_name, smc1a_down$gene_name, h3k27ac_up$gene_name, h3k27ac_down$gene_name)))
print(SMC1a_up_H3k27ac_up_value)

# Step 6: SMC1a_up & H3k27ac_down
SMC1a_up_H3k27ac_down_value <- length(intersect(smc1a_up$gene_name, h3k27ac_down$gene_name)) - 
  length(Reduce(intersect, list(smc1a_up$gene_name, smc1a_down$gene_name, h3k27ac_up$gene_name, h3k27ac_down$gene_name)))
print(SMC1a_up_H3k27ac_down_value)

# Step 7: SMC1a_down & H3k27ac_up
SMC1a_down_H3k27ac_up_value <- length(intersect(smc1a_down$gene_name, h3k27ac_up$gene_name)) - 
  length(Reduce(intersect, list(smc1a_up$gene_name, smc1a_down$gene_name, h3k27ac_up$gene_name, h3k27ac_down$gene_name)))
print(SMC1a_down_H3k27ac_up_value)

# Step 8: SMC1a_down & H3k27ac_down
SMC1a_down_H3k27ac_down_value <- length(intersect(smc1a_down$gene_name, h3k27ac_down$gene_name)) - 
  length(Reduce(intersect, list(smc1a_up$gene_name, smc1a_down$gene_name, h3k27ac_up$gene_name, h3k27ac_down$gene_name)))
print(SMC1a_down_H3k27ac_down_value)


UU <- euler(c(
  SMC1a_up = SMC1a_up_value,
  SMC1a_down = SMC1a_down_value,
  H3k27ac_up = H3k27ac_up_value,
  H3k27ac_down = H3k27ac_down_value,
  "SMC1a_up&H3k27ac_up" = SMC1a_up_H3k27ac_up_value,
  "SMC1a_up&H3k27ac_down" = SMC1a_up_H3k27ac_down_value,
  "SMC1a_down&H3k27ac_up" = SMC1a_down_H3k27ac_up_value,
  "SMC1a_down&H3k27ac_down" = SMC1a_down_H3k27ac_down_value
))

# Check if any sets are empty and remove them if they are

#Figure 3D
eulerplot<-plot(UU, quantities = F, fills = c("#4F94CD","#36648B","#E3072A","#8B1A1A"),edges = F)
ggsave("final_fig/Figure3D_euler_plot.png", eulerplot, width = 120, height = 100,units ="mm" ,dpi = 300)
eulerplot<-plot(UU, quantities = T, fills = c("#4F94CD","#36648B","#E3072A","#8B1A1A"),edges = F)
ggsave("final_fig/Figure3D_euler_plot_quant.png", eulerplot, width = 120, height = 100,units ="mm" ,dpi = 300)
eulerplot<-plot(UU, quantities = T, fills = c("#4F94CD","#36648B","#E3072A","#8B1A1A"),edges = F,label=F)
ggsave("final_fig/Figure3D_euler_plot_nolabels.png", eulerplot, width = 120, height = 100,units ="mm" ,dpi = 300)



# Figure 4A-B
#######
mycolors<-c("firebrick4", "darkorange2", "gold2", "olivedrab", "steelblue4", "orchid")
dcolors<-c("azure2","burlywood4", "lightgreen", "plum1", "cornflowerblue", "brown2")

activatedData<-data[which(data$log2FoldChange_ITLvsC_siC>1&data$padjvalue_ITLvsC_siC<0.05),]
my_comparisons <- list( c("None", "Active"), c("Down", "Up"), c("DownActive", "UpActive"), c("Down", "DownActive"), c("Up", "UpActive"), c("Up", "None"), c("Down", "None"), c("UpActive", "Active"), c("DownActive", "Active"))

#max_y_value <- max(activatedData$log2FoldChange_ITLvsC_siC, na.rm = TRUE)
g <- ggboxplot(data = activatedData, x = "EnhancerStatus", y = "log2FoldChange_ITLvsC_siC", fill = "EnhancerStatus") +
  geom_point(alpha = 0.1, position = position_jitter(height = .5, width = .5)) +
  scale_fill_manual(values = dcolors) + ## (values = dcolors, guide=F) to remove the legend of the fill
  #facet_grid(.~type.y) +
 # stat_compare_means(comparisons = my_comparisons) +
  # ylim(c(NA, max_y_value)) + # Setting the y-axis limit
  #labs(title = "ITL Stimulation", subtitle = "Activated Genes", y = "log2(FC)", x = "SMC1A Increased in Enhancers", fill = "Increased SMC1A\nin Enhancers") +
  labs(y = "log2(FC)", x = "SMC1A Increased in Enhancers", fill = "Increased SMC1A\nin Enhancers") +
  theme(
    plot.title = element_text(color = "black", size = 8),
    plot.subtitle = element_text(color = "black", size = 8),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
#ggsave("final_fig/Figure4A_ITL stimulation in control - activated - SMC1A Increased in Enhancers.png",g,width = 130, height = 88,units ="mm", dpi = 300)
#
##for revisions
fig4a_data<-activatedData
write.table(fig4a_data,"~/SofiaP/smc1a/revisions/fig4a-b_data.txt",quote = F, row.names = F, col.names = T)
##

analogous_set1 <- c("aliceblue", "wheat4", "palegreen", "orchid1", "mediumslateblue", "indianred1")
analogous_set2 <- c("lightcyan", "tan4", "springgreen", "violet", "royalblue", "red3")
analogous_set3 <- c("lightblue1", "darkgoldenrod", "mediumseagreen", "palevioletred1", "slateblue", "firebrick1")

repressedData<-data[which(data$log2FoldChange_ITLvsC_siC<(-1)&data$padjvalue_ITLvsC_siC<0.05),]
g <- ggboxplot(data = repressedData, x = "EnhancerStatus", y = "log2FoldChange_ITLvsC_siC", fill = "EnhancerStatus") +
  geom_point(alpha = 0.1, position = position_jitter(height = .5, width = .5)) +
  scale_fill_manual(values = dcolors) + ## guide = F to remove the legend of the fill
  #facet_grid(.~type.y) +
  #stat_compare_means(comparisons = my_comparisons) +
  #labs(title = "ITL Stimulation", subtitle = "Repressed Genes", y = "log2(FC)", x = "SMC1A Increased in Enhancers", fill = "Increased SMC1A\nin Enhancers") +
  labs(y = "log2(FC)", x = "SMC1A Increased in Enhancers", fill = "Increased SMC1A\nin Enhancers") +
  theme(
    #plot.title = element_text(color = "black", size = 19),
    #plot.subtitle = element_text(color = "black", size = 15),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
#ggsave("final_fig/SupplFigure7A_ITL stimulation in control - repressed - SMC1A Increased in Enhancers.png",g,width = 130, height = 88,units ="mm", dpi = 300)

##for revisions
figS8ab_data<-repressedData
write.table(figS8ab_data,"~/SofiaP/smc1a/revisions/figS8a-b_data.txt",quote = F, row.names = F, col.names = T)
##

## promoters
#analogous_colors <- c("firebrick3", "orange", "goldenrod", "darkolivegreen", "steelblue3", "mediumorchid")
promoter_colors_set1 <- c("lightcyan", "sienna", "palegreen2", "orchid", "lightskyblue", "firebrick")
promoter_colors_set3 <- c("lightyellow", "tan", "mediumspringgreen", "thistle", "deepskyblue", "tomato")
complementary_colors <- c("dodgerblue", "blue", "purple", "red", "orange", "green")
g <- ggboxplot(data = activatedData, x = "PromoterStatus", y = "log2FoldChange_ITLvsC_siC", fill = "PromoterStatus") +
  geom_point(alpha = 0.1, position = position_jitter(height = .5, width = .5)) +
  scale_fill_manual(values = dcolors) + ## guide = F to remove the legend of the fill
  #facet_grid(.~type.y) +
  #stat_compare_means(comparisons = my_comparisons) +
#  labs(title = "ITL Stimulation", subtitle = "Activated Genes", y = "log2(FC)", x = "SMC1A Increased in Promoters", fill = "Increased SMC1A\nin Enhancers") +
  labs(y = "log2(FC)", x = "SMC1A Increased in Promoters", fill = "Increased SMC1A\nin Enhancers") +
  theme(
    #plot.title = element_text(color = "black", size = 8),
    #plot.subtitle = element_text(color = "black", size = 8),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
#ggsave("./ITL stimulation in control - activated - SMC1A Increased in Promoters.png",g,width = 8, height = 8, dpi = 300)
ggsave("final_fig/Figure4B_ITL stimulation in control - activated - SMC1A Increased in Promoters.png",g,width = 130, height = 88,units ="mm", dpi = 300)


complementary_colors <- c("dodgerblue", "blue", "purple", "red", "orange", "green")
g <- ggboxplot(data = repressedData, x = "PromoterStatus", y = "log2FoldChange_ITLvsC_siC", fill = "PromoterStatus") +
  geom_point(alpha = 0.1, position = position_jitter(height = .5, width = .5)) +
  #geom_point(shape = 16, alpha = 0.5, position = position_jitter(height = .5, width = .5)) + # Circle shape for promoters
  #geom_point(alpha = 0.5, size = 1.5, position = position_jitter(height = 0.3, width = 0.3)) + # Increase jitter and reduce point size
  scale_fill_manual(values = dcolors) + ## guide = F to remove the legend of the fill
  #facet_grid(.~type.y) +
  #stat_compare_means(comparisons = my_comparisons) +
 # labs(title = "ITL Stimulation", subtitle = "Repressed Genes", y = "log2(FC)", x = "SMC1A Increased in Promoters", fill = "Increased SMC1A\nin Enhancers") +
  labs(y = "log2(FC)", x = "SMC1A Increased in Promoters", fill = "Increased SMC1A\nin Enhancers") +
  theme(
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
#ggsave("./ITL stimulation in control - repressed - SMC1A Increased in Promoters.png",g,width = 8, height = 8, dpi = 300)
ggsave("final_fig/SupplFigure7B_ITL stimulation in control - repressed - SMC1A Increased in Promoters.png",g,width = 130, height = 88,units ="mm", dpi = 300)
#
###figure 4c
### apothhkeysh
## t test for up in both promoter and enhancer
both<-data[which(data$up_promoters=="Yes"|data$up_enhancers=="Yes"|data$up_active_enhancers=="Yes"|data$up_active_promoters=="Yes"),c(1,2,10,8,17)]
both<-data[which(data$up_active_enhancers=="Yes"|data$up_active_promoters=="Yes"),]
both$up_enh_prom<-"No"
both$up_enh_prom[which(both$up_active_promoters=="Yes"&both$up_active_enhancers=="Yes")]<-"Yes"
both<-both[,c("Gene","log2FoldChange_ITLvsC_siC","up_active_promoters","up_active_enhancers","up_enh_prom")]
library(reshape2)
both_melt_ITLvsC_siC <- melt(both, id=c("Gene","log2FoldChange_ITLvsC_siC"))
both_melt_ITLvsC_siC_yes<-both_melt_ITLvsC_siC[which(both_melt_ITLvsC_siC$value=="Yes"),]

x_labels <- c("up_active_promoters" = "Up Active Promoters",
              "up_active_enhancers" = "Up Active Enhancers",
              "up_enh_prom" = "Up Active Promoters & \nEnhancers")

gb <- ggboxplot(both_melt_ITLvsC_siC_yes, y = "log2FoldChange_ITLvsC_siC", x = "variable", fill = "variable") +
  geom_point(alpha = 0.1, position = position_jitter(height = .5, width = .5)) + stat_compare_means(comparisons = list(c("up_active_enhancers", "up_enh_prom"), c("up_active_promoters", "up_enh_prom"))) +
  #scale_fill_manual(values = c( "#57164e","#8f174e" ,"#ff1a4c"), guide=F) + ## guide = F to remove the legend of the fill
  scale_fill_manual(values = c( "goldenrod1", "lightsalmon", "firebrick1"), guide=F) + ## guide = F to remove the legend of the fill
  scale_x_discrete(labels = x_labels)+
 ## labs(title = "Expression upon ITL Stimulation", subtitle = "SMC1a bound Genes", y = "log2(FC)", x = "Enhancer/Promoter Status", fill = "Status") +
  labs(y = "log2(FC)", x = "Enhancer/Promoter Status", fill = "Status") +
  theme(
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))
#legend.position = "none" # Remove legend if not needed

#ggsave("./ITL stimulation up bound Promoters AND Enhancers.png",gb,width = 8, height = 8, dpi = 300)
ggsave("final_fig/Figure4C_ITL stimulation up bound Promoters AND Enhancers.png",gb,width = 120, height = 90,units ="mm" ,dpi = 300)

##for revisions
fig4c_data<-both_melt_ITLvsC_siC_yes
write.table(fig4c_data,"~/SofiaP/smc1a/revisions/fig4c_data.txt",quote = F, row.names = F, col.names = T)
##

## GROUP A B C
groupa<-read.table("~/SofiaP/smc1a/groupa.txt")
groupb<-read.table("~/SofiaP/smc1a/groupb.txt")
groupc<-read.table("~/SofiaP/smc1a/groupc.txt")
ego_groupa <- enrichGO(gene = unique(groupa$V1),
                       keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05, readable = FALSE)
d_groupa<-dotplot(ego_groupa, showCategory=15,orderBy = "x")+ #, title = "Enriched terms of group A genes"
  theme(
   # plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))
#ggsave("./Dotplot group a .png",d_groupa,width = 6, height = 5, dpi = 300)
ggsave("final_fig/Figure4F_Dotplot - Enriched terms of group A genes.png",d_groupa,width = 100, height = 88,units ="mm" ,dpi = 300)


ego_groupbc <- enrichGO(gene = c(unique(groupb$V1),unique(groupc$V1)),
                        keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05, readable = FALSE)
ego_groupbc_SIMPL<-simplify(ego_groupbc,cutoff = 0.7,semData=NULL)
d_groupbc_simpl<-dotplot(ego_groupbc_SIMPL, showCategory=7,orderBy = "x")+ #, title = "Enriched terms of group B & C genes
  theme(
    #plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))
#ggsave("./Dotplot group b and c _ simpl _ theme.png",d_groupbc_simpl,width = 6, height = 5, dpi = 300)
ggsave("final_fig/Figure4E_Dotplot - Enriched terms of group B&C genes.png",d_groupbc_simpl,width = 100, height = 88,units ="mm" ,dpi = 300)


#######
groupbc<-c(groupb$V1,groupc$V1) #345 -> unique 325
data_itl<-data[which( data$Gene%in%groupbc), c(1,4,5)]
colnames(data_itl)<-c("Gene", "log2FC", "pvalue")
data_itl$expr[data_itl$log2FC>0]<-"up"
data_itl$expr[data_itl$log2FC<0]<-"down"
data_itl$expr<-as.factor(data_itl$expr)
#
data_itl$delabel <- NA
data_itl$delabel[which(abs(data_itl$log2FC)>=1)]<-as.character(data_itl$Gene[which(abs(data_itl$log2FC)>=1)])

##Figure 4G
g<-ggplot(data_itl, aes(x=log2FC, y=-log10(pvalue), color=expr, label=delabel)) + 
  geom_point() + 
  scale_color_manual(values=c("steelblue4", "firebrick4")) +
  geom_text_repel(box.padding = 0.3, col="black", size=8/.pt) +
  #ggtitle(label="Group B & C genes", subtitle ="Gene Expression Changes upon ITL stimulation") +
  xlab("log2FC") + 
  ylab("-log10(p-value)") + 
  theme(
    #plot.title = element_text(color="black", size=19),
    #plot.subtitle = element_text(color="black", size=15),
    axis.title.x = element_text(color="black", size=8),
    axis.title.y = element_text(color="black", size=8),
    axis.text.x = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
#ggsave("./Volcano - group b and c - itl stimulation.png",g,width = 6, height = 6, dpi = 300)
ggsave("final_fig/Figure4G_Volcano - group b and c - itl stimulation.png",g,width = 88, height = 88,units ="mm" ,dpi = 300)



## Figure 6
######### SLE monocytes female vs male
mtsq_untr_a<-read.table("~/SofiaP/smc1a/metaseqr_all_out_untr_female_vs_untr_male.txt", header = T)

mtsq_untr_df_volc<-mtsq_untr_a
mtsq_untr_df_volc$expression<-"unregulated"
mtsq_untr_df_volc$expression[mtsq_untr_df_volc$p.value_deseq2<0.05&mtsq_untr_df_volc$log2_normalized_fold_change_untr_female_vs_untr_male>1]<-"up"
mtsq_untr_df_volc$expression[mtsq_untr_df_volc$p.value_deseq2<0.05&mtsq_untr_df_volc$log2_normalized_fold_change_untr_female_vs_untr_male<(-1)]<-"down"
mtsq_untr_df_volc$expression<-as.factor(mtsq_untr_df_volc$expression)
#
mtsq_untr_df_volc$delabel <- NA
mtsq_untr_df_volc$delabel[which(mtsq_untr_df_volc$expression!="unregulated")]<-as.character(mtsq_untr_df_volc$gene_name[which(mtsq_untr_df_volc$expression!="")])
#
g<-ggplot(mtsq_untr_df_volc, aes(x=log2_normalized_fold_change_untr_female_vs_untr_male, y=-log10(p.value_deseq2), color=expression, label=delabel)) + 
#g<-ggplot(mtsq_untr_df_volc, aes(x=log2_normalized_fold_change_untr_female_vs_untr_male, y=-log10(p.value_deseq2), color=expression)) + 
  geom_point() + 
  scale_color_manual(values = c("down" = "steelblue4", "up" = "firebrick4")) +
  geom_text_repel(box.padding = 0.3, col="black", size=8/.pt) +
 # ggtitle(label="SLE female vs. male monocytes") +
  xlab("log2FC") + 
  ylab("-log10(p-value)") + 
  theme(
   # plot.title = element_text(color="black", size=19),
    #plot.subtitle = element_text(color="black", size=15),
    axis.title.x = element_text(color="black", size=8),
    axis.title.y = element_text(color="black", size=8),
    axis.text.x = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

#ggsave("./Volcano - SLE female male monocytes -f.png",g,width = 6, height = 6, dpi = 300)
ggsave("final_fig/Figure6Α_Volcano - SLE female male monocytes.png",g,width = 88, height = 88,units ="mm" ,dpi = 300)



up_SLEfm<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male>1),] #596
down_SLEfm<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male<(-1)),] #149

ego_up_SLEfm <- enrichGO(gene = unique(up_SLEfm$gene_name),
                         keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05, readable = FALSE)
ego_up_SLEfm_SIMPL<-simplify(ego_up_SLEfm,cutoff = 0.7,semData=NULL)
d<-dotplot(ego_up_SLEfm_SIMPL, showCategory=10,orderBy = "x")+#, title = "Enriched terms SLE female vs. male \nup-regulated genes")+
  theme(
   # plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))
#ggsave("./enriched terms_SLE female vs male_up.png",d,width = 8, height = 9, dpi = 300)
ggsave("final_fig/Figure6B_Dotplot - Enriched terms of SLE female vs male_up degs.png",d,width = 120, height = 90,units ="mm" ,dpi = 300)
ggsave("final_fig/Figure6B_Dotplot - Enriched terms of SLE female vs male_up degs - SIMILAR TO C.png",d,width = 120, height = 110,units ="mm" ,dpi = 300)

ego_down_SLEfm <- enrichGO(gene = unique(down_SLEfm$gene_name),
                           keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05, readable = FALSE)
d<-dotplot(ego_down_SLEfm, showCategory=10,orderBy = "x")+#, title = "Enriched terms SLE female vs. male \ndown-regulated genes")+
  theme(
    #plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))
ggsave("final_fig/Figure6C_Dotplot - enriched terms_SLE female vs male_down.png",d,width = 120, height = 110,units ="mm" ,dpi = 300)


groupbc<-unique(rbind(groupb,groupc))
target_genes<-intersect(groupbc$V1,activatedData$Gene) #277
up_act_enh<-data[which(data$EnhancerStatus=="UpActive"),]

mtsq_untr_a_with_target_genes<-mtsq_untr_a
mtsq_untr_a_with_target_genes$SMC1Aregulated<-"No"
mtsq_untr_a_with_target_genes$SMC1Aregulated[which(mtsq_untr_a_with_target_genes$gene_name%in%target_genes)]<-"Yes"
#write.table(mtsq_untr_a_with_target_genes,"./final_fig/SupplTable8_metaseq_output_SLEfemalevs.male_with_SMC1Atarget_genes.txt", quote = F,col.names = T, row.names = F)


##  heatmap
mtsq_untr_a_with_target_genes

mtsq_untr_a_with_target_genes_heat<-mtsq_untr_a_with_target_genes
colnames(mtsq_untr_a_with_target_genes_heat)<-gsub("log2_normalized_counts_","",as.character(colnames(mtsq_untr_a_with_target_genes_heat)))
cold<-read.table("~/SofiaP/smc1a/coldata")
rownames(cold)<-cold[,1]
cold<-cold[,-1]
cold_no309<-cold[!(cold$V4 == "P309"),]
cold_no309<-cold_no309[which(cold_no309$V3=="untr"),]
mtsq_untr_a_with_target_genes_heat<-mtsq_untr_a_with_target_genes_heat[which(mtsq_untr_a_with_target_genes_heat$SMC1Aregulated=="Yes"),]
rownames(mtsq_untr_a_with_target_genes_heat)<-mtsq_untr_a_with_target_genes_heat$gene_name
mtsq_untr_a_with_target_genes_heat<-mtsq_untr_a_with_target_genes_heat[,rownames(cold_no309)]

colnames(cold_no309)<-c("sex","treatment","patient")
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
norm_zscore <- t(apply(mtsq_untr_a_with_target_genes_heat, 1, cal_z_score))
#annorow<-as.data.frame(t$event_type)
#rownames(annorow)<-rownames(t)
library(pheatmap)
ph<-pheatmap(norm_zscore, annotation_col=cold_no309[,c(1,2)],cluster_rows=T,cluster_cols = F,show_rownames=F,main= "heatmap of SMC1a target genes",breaks=seq(-2.5,2.5,length.out =100),col = colorRampPalette(c("navy", "white", "firebrick3"))(100))
#ggsave("~/SMC1A_project/quant_seq_SLE_females-males/deseq_no309/Heatmap_212_genes.png",ph,width = 7, height = 5, dpi = 120)

## or
ccols <- c("firebrick4", "darkorange2", "gold2", "olivedrab3", "forestgreen", "steelblue3", "royalblue4", "mediumorchid4")
cols <- colorRampPalette(c("steelblue4", "white", "firebrick4"))(100)

annotation_colors <- list(
  # Assume 'cold_no309[,1]' corresponds to "Response" and 'cold_no309[,2]' to "Colonized"
  sex = c("male" = "#00A9FF", "female" = "#F8766D"), # Custom colors for 'Response'
  treatment = c("untr" = NA)             # Custom colors for 'Colonized'
)

p<-pheatmap(
  as.matrix(norm_zscore),
  show_colnames = F,
  show_rownames = F,
  annotation_col=cold_no309[,c(1,2)],
  #labels_col = c("Response", "Colonized"), 
  color = cols,                   # Color palette
  breaks = seq(-2, 2, length.out = 101),
  key = TRUE,                     # Show color key
  key.title = "log2(FC)", 
  #main = "Gene Expression Changes in SMC1a Target Genes", 
  cluster_rows = TRUE,            # Cluster rows
  cluster_cols = FALSE,           # Do not cluster columns
  scale = "none",                 # Do not scale
  fontsize_row = 8,               # Row label font size
  fontsize_col = 8 ,              # Column label font size
  angle_col = 0,
  annotation_colors=annotation_colors
)
ggsave("final_fig/Figure6D_Heatmap quantseq SLE - SMC1a target genes.png",p,width = 110, height = 90,units ="mm" ,dpi = 300)

##for revisions
fig5e_annot_data<-cold_no309[,c(1,2)]
fig5e_normcounts_data<-as.data.frame(norm_zscore)
write.table(fig5e_annot_data,"~/SofiaP/smc1a/revisions/fig5e_annot_data.txt",quote = F, row.names = T, col.names = T)
write.table(fig5e_normcounts_data,"~/SofiaP/smc1a/revisions/fig5e_normcounts_data.txt",quote = F, row.names = T, col.names = T)
##

## Suppl Fig 4
lupus<-read.table("gene_statistics_SLEvsHEALTHY-comb-ex-rem.txt", header=T)
stimulation<-read.table("gene_statistics_IFN-TNF-LPS_vs_unstimulated_siControl.txt", header = T)

#lupus_sig<- lupus[which(lupus$padj<0.05),] #1827 
lupus_sig<-lupus[which(lupus$padj<0.05&abs(lupus$log2FoldChange)>1),] #1089

#stimulation_sig<-stimulation[which(stimulation$padj<0.05&abs(stimulation$log2FoldChange)>1),] #3632
#common_degs_lupus_stim<-intersect(lupus_sig$gene_name,stimulation_sig$gene_name) #414

stimulationDATA_sig<-data[which(data$padjvalue_ITLvsC_siC<0.05&abs(data$log2FoldChange_ITLvsC_siC)>1),] #3632
stimulationDATA_sig<-stimulationDATA_sig[!duplicated(stimulationDATA_sig$Gene),]
colnames(stimulationDATA_sig)[1]<-"gene_name"
#lupus_treat_sig<-merge(lupus_sig,stimulation_sig, by="gene_name") #
#lupus_treat_sig<-na.omit(lupus_treat_sig) #414

lupus_treat_sig<-merge(lupus_sig,stimulationDATA_sig, by="gene_name")
lupus_treat_sig<-na.omit(lupus_treat_sig) #386

#to y einai to df_c_treat.....
colnames(lupus_treat_sig)
#lupus_treat_sig<-lupus_treat_sig[,c(1,4,8,10,14)] # with the stimulation file, not data
lupus_treat_sig<-lupus_treat_sig[,c(1,4,8,11,12)]
colnames(lupus_treat_sig)<-c("gene_name","log2FoldChange-Lupus_Signature","padj-Lupus_Signature","log2FoldChange-ITL_Signature","padj-ITL_Signature")  


ggplot(lupus_treat_sig, aes(x=`log2FoldChange-Lupus_Signature`, y=`log2FoldChange-ITL_Signature`))+ geom_vline(xintercept=0, linetype="dashed", size=.1) +
  geom_hline(yintercept=0, linetype="dashed", size=.1) +
  geom_point()+labs(title="Common Differentially Expressed Genes \nin treatment comparison of si-Control samples and SLE vs. Healthy comparison", x="Log2FC in SLE vs. Healthy", y="Log2FC in IFN-TNF-LPS v. Untreated")+
  geom_text(aes(label=gene_name),hjust=.5, vjust=-.8, size=3)
#ggsave("/home/user/gender_bias/new_data/log2fc.png",u_plot,width = 10, height = 9, dpi = 120)

g<-ggplot(lupus_treat_sig, aes(x=`log2FoldChange-Lupus_Signature`, y=`log2FoldChange-ITL_Signature`))+ 
  geom_vline(xintercept=0, linetype="dashed", size=.1)+
  geom_hline(yintercept=0, linetype="dashed", size=.1) +
  geom_point() + 
  #labs(title="Common Differentially Expressed Genes \nin treatment comparison of si-Control samples and SLE vs. Healthy comparison", x="Log2FC in SLE vs. Healthy", y="Log2FC in IFN-TNF-LPS v. Untreated")+
 # geom_text(aes(label=gene_name),hjust=.5, vjust=-.8, size=8/.pt)+
  geom_text_repel(aes(label=gene_name),box.padding = 0.3, col="black", size=8/.pt) +
  #ggtitle(label="Group B & C genes", subtitle ="Gene Expression Changes upon ITL stimulation") +
  xlab("Log2FC in SLE vs. Healthy") + 
  ylab("Log2FC in IFN-TNF-LPS v. Untreated") + 
  theme(
    #plot.title = element_text(color="black", size=19),
    #plot.subtitle = element_text(color="black", size=15),
    axis.title.x = element_text(color="black", size=8),
    axis.title.y = element_text(color="black", size=8),
    axis.text.x = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
#ggsave("final_fig/SupplFigure4B_ommon Differentially Expressed Genes \nin treatment comparison of si-Control samples and SLE vs. Healthy comparison.png",g,width = 120, height = 120,units ="mm" ,dpi = 300)

gg <- ggplot(lupus_treat_sig, aes(x=`log2FoldChange-Lupus_Signature`, y=`log2FoldChange-ITL_Signature`)) + 
  geom_vline(xintercept=0, linetype="dashed", size=.1)+
  geom_hline(yintercept=0, linetype="dashed", size=.1) +
  geom_point() + 
  geom_text_repel(aes(label=gene_name),box.padding = 0.3, col="black", size=8/.pt) +
  # geom_smooth(method="loess", se=F) + 
  #xlim(c(0, 0.1)) + 
  #ylim(c(0, 500000)) + 
  geom_jitter() + 
  geom_smooth(method="lm", se=F)+
  stat_cor(method = "spearman")+
  xlab("Log2FC in SLE vs. Healthy") + 
  ylab("Log2FC in IFN-TNF-LPS v. Untreated") + 
  theme(
    #plot.title = element_text(color="black", size=19),
    #plot.subtitle = element_text(color="black", size=15),
    axis.title.x = element_text(color="black", size=8),
    axis.title.y = element_text(color="black", size=8),
    axis.text.x = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
#ggsave("final_fig/SupplFigure4B_Common Differentially Expressed Genes \nin treatment comparison of si-Control samples and SLE vs. Healthy comparison - spearman cor.png",gg,width = 120, height = 120,units ="mm" ,dpi = 300)

##for revisions
figS5b_data<-lupus_treat_sig
write.table(figS5b_data,"~/SofiaP/smc1a/revisions/figS5b_data.txt",quote = F,row.names = F, col.names = T)
##
### fischer test
up.itl.lupus<-matrix(c(length(com_up_lupus_itl),(length(up_itl_sig)-length(com_up_lupus_itl)),(length(up_lupus_sig)-length(com_up_lupus_itl)),
                       (length(lupus_treat_sig$gene_name)-length(up_lupus_sig)-length(up_itl_sig)+length(com_up_lupus_itl))), 2, 2, dimnames = list(itl = c("up.itl", "not.up.itl"),sle = c("up.sle", "not.up.sle")))
                                                                                                                                                                                                                                                                                                                            "not.high.rated")))
fisher.test(up.itl.lupus)                 

#enrich_gene_list<-list(Lupus=unique(lupus_sig$gene_name),Stimulation=unique(stimulation_sig$gene_name))
enrich_gene_list<-list(Lupus=unique(lupus_sig$gene_name),Stimulation=unique(stimulationDATA_sig$gene_name))

enrich_go_lupus.stim<- compareCluster(geneCluster = enrich_gene_list,
                                      keyType       = 'SYMBOL',
                                      ont           = "ALL",
                                      OrgDb = org.Hs.eg.db,
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.05,
                                      fun =  "enrichGO")
out_SIMPL<-simplify(enrich_go_lupus.stim,cutoff = 0.7,semData=NULL)
d_simpl_lupus.stim<-dotplot(out_SIMPL,title = "GO enrichment analysis \nLupus and ITL signature")
#ggsave("GO enrichment analysis Lupus and ITL signature.png",d_simpl_lupus.stim,width = 6, height = 6, dpi = 300)

d<-dotplot(out_SIMPL)+
  theme(
    #plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))
ggsave("final_fig/SupplFigure4A_GO enrichment analysis Lupus and ITL signature.png",d,width = 120, height = 110,units ="mm" ,dpi = 300)




















#Transcriptome analysis revealed a total of 745 genes with sexual differential expression in SLE monocytes (|log2-FC|>1, p-value<0.05) 
up_SLEfm<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male>1),] #596
down_SLEfm<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male<(-1)),] #149


uuu<-intersect(target_genes,up_SLEfm$gene_name) #62 534
#write.table(uuu, "intersection_of_target_genes_and_up_regulated_in_SLEfemales.txt", quote = F, row.names = F, col.names = F)
ddd<-intersect(target_genes,down_SLEfm$gene_name) #2 147
#The chi-square statistic is 12.4611. The p-value is .000416. The result is significant at p < .05.
du<-union(uuu,ddd)

up_SLEfm<-mtsq_untr_a[which(mtsq_untr_a$FDR_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male>1),] #392
down_SLEfm<-mtsq_untr_a[which(mtsq_untr_a$FDR_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male<(-1)),] #65

intersect(target_genes,up_SLEfm$gene_name) #51
intersect(target_genes,down_SLEfm$gene_name) #0





library(readxl)
cl_mon<-read_xlsx("class_monocytes_fvsm_sle.xlsx")
colnames(cl_mon) # [9] "Log2 fold change\nF vs. M"                      [10] "Adj. P value\nF vs. M"  
cl_mon_up<-cl_mon[which(cl_mon[,9]>1&cl_mon[,10]<0.05),] #11
intersect(cl_mon_up$`Gene ID`,target_genes) #il6
cl_mon_up<-cl_mon[which(cl_mon[,9]>0&cl_mon[,10]<0.05),] #223
intersect(cl_mon_up$`Gene ID`,target_genes) #27 196

cl_mon_down<-cl_mon[which(cl_mon[,9]<(-1)&cl_mon[,10]<0.05),] #17
cl_mon_down<-cl_mon[which(cl_mon[,9]<0&cl_mon[,10]<0.05),] #86
intersect(cl_mon_down$`Gene ID`,target_genes) #0 86

intersect(cl_mon_up$`Gene ID`,du)#9
intersect(cl_mon_down$`Gene ID`,du)#0

## FIsher test
#SMC1a-target genes and non-target genes vs. (female vs. male) Degs in SLE monocytes.
fisher.test(table(which(mtsq_untr_a$p.value_deseq2<0.05), data$up_enhancers))

##

## create heatmap
mtsq_untr_a_sig<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&abs(mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male)>1),]
targetgenes<-data[which(data$log2FoldChange_ITLvsC_siC>1&data$padjvalue_ITLvsC_siC<0.05&data$Gene%in%groupbc$V1),"Gene"]
# define targets

mtsq_untr_a_sig_smc1a<-mtsq_untr_a_sig[mtsq_untr_a_sig$gene_name%in%targetgenes,] #11 "GTPBP1"  "TRIM25"  "N4BP2L1" "PML"     "PMAIP1"  "GBP5"    "PARP14"  "TNFAIP2" "ADA"     "PELI1"   "SLC9A8" 

common_smc1a_untr_a<-intersect(mtsq_untr_a_sig$gene_name, targetgenes) #64
mtsq_untr_a_smc1a<-mtsq_untr_a_sig[which(mtsq_untr_a_sig$gene_name%in%common_smc1a_untr_a),] # only 5/102 down


mtsq_untr_a_smc1a_heat<-mtsq_untr_a_smc1a
colnames(mtsq_untr_a_smc1a_heat)<-gsub("log2_normalized_counts_","",as.character(colnames(mtsq_untr_a_smc1a_heat)))
rownames(mtsq_untr_a_smc1a_heat)<-mtsq_untr_a_smc1a_heat$gene_name
cold<-read.table("~/SofiaP/smc1a/coldata")
rownames(cold)<-cold[,1]
cold<-cold[,-1]
cold_no309<-cold[!(cold$V4 == "P309"),]
cold_no309<-cold_no309[which(cold_no309$V3=="untr"),]

mtsq_untr_a_smc1a_heat<-mtsq_untr_a_smc1a_heat[,rownames(cold_no309)]
mtsq_untr_a_smc1a_heat<-mtsq_untr_a_smc1a
colnames(mtsq_untr_a_smc1a_heat)<-gsub("log2_normalized_counts_","",as.character(colnames(mtsq_untr_a_smc1a_heat)))
rownames(mtsq_untr_a_smc1a_heat)<-mtsq_untr_a_smc1a_heat$gene_name
#cold<-read.table("coldata")
#rownames(cold)<-cold[,1]
#cold<-cold[,-1]
#cold_no309<-cold[!(cold$V4 == "P309"),]
#cold_no309<-cold_no309[which(cold_no309$V3=="untr"),]
colnames(cold_no309)<-c("sex","treatment","patient")
mtsq_untr_a_smc1a_heat<-mtsq_untr_a_smc1a_heat[,rownames(cold_no309)]


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
mtsq_untr_a
norm_zscore <- t(apply(mtsq_untr_a_smc1a_heat, 1, cal_z_score))
#annorow<-as.data.frame(t$event_type)
#rownames(annorow)<-rownames(t)
library(pheatmap)
ph<-pheatmap(norm_zscore, annotation_col=cold_no309[,c(1,2)],cluster_rows=T,cluster_cols = F,show_rownames=F,main= "heatmap of SMC1a target genes",breaks=seq(-2.5,2.5,length.out =100),col = colorRampPalette(c("navy", "white", "firebrick3"))(100))
#ggsave("~/SMC1A_project/quant_seq_SLE_females-males/deseq_no309/Heatmap_212_genes.png",ph,width = 7, height = 5, dpi = 120)



### for all 277 genes START
mtsq_untr_a_heat<-mtsq_untr_a
colnames(mtsq_untr_a_heat)<-gsub("log2_normalized_counts_","",as.character(colnames(mtsq_untr_a_heat)))
cold<-read.table("~/SofiaP/smc1a/coldata")
rownames(cold)<-cold[,1]
cold<-cold[,-1]
cold_no309<-cold[!(cold$V4 == "P309"),]
cold_no309<-cold_no309[which(cold_no309$V3=="untr"),]

mtsq_untr_a_heat<-mtsq_untr_a_heat[,rownames(cold_no309)]
colnames(mtsq_untr_a_heat)<-gsub("log2_normalized_counts_","",as.character(colnames(mtsq_untr_a_heat)))
rownames(mtsq_untr_a_heat)<-mtsq_untr_a_heat$gene_name
colnames(cold_no309)<-c("sex","treatment","patient")
#mtsq_untr_a_smc1a_heat<-mtsq_untr_a_smc1a_heat[,rownames(cold_no309)]

norm_zscore <- t(apply(mtsq_untr_a_heat, 1, cal_z_score))
#annorow<-as.data.frame(t$event_type)
#rownames(annorow)<-rownames(t)
library(pheatmap)
ph<-pheatmap(norm_zscore, annotation_col=cold_no309[,c(1,2)],cluster_rows=T,cluster_cols = F,show_rownames=F,main= "heatmap of SMC1a target genes",breaks=seq(-2.5,2.5,length.out =100),col = colorRampPalette(c("navy", "white", "firebrick3"))(100))
#ggsave("~/SMC1A_project/quant_seq_SLE_females-males/deseq_no309/Heatmap_212_genes.png",ph,width = 7, height = 5, dpi = 120)
### for all 277 genes END



ccols <- c("firebrick4", "darkorange2", "gold2", "olivedrab3", "forestgreen", "steelblue3", "royalblue4", "mediumorchid4")
cols <- colorRampPalette(c("steelblue4", "white", "firebrick4"))(100)


p<-pheatmap(
  as.matrix(norm_zscore),
  show_colnames = F,
  show_rownames = F,
  annotation_col=cold_no309[,c(1,2)],
  #labels_col = c("Response", "Colonized"), 
  color = cols,                   # Color palette
  breaks = seq(-2, 2, length.out = 101),
  key = TRUE,                     # Show color key
  key.title = "log2(FC)", 
  main = "Gene Expression Changes in SMC1a Target Genes", 
  cluster_rows = TRUE,            # Cluster rows
  cluster_cols = FALSE,           # Do not cluster columns
  scale = "none",                 # Do not scale
  fontsize_row = 8,               # Row label font size
  fontsize_col = 12 ,              # Column label font size
  angle_col = 0
)
ggsave("./heatmap quantseq sle - SMC1a target genes.png",p,width = 8, height = 8, dpi = 300)

## deseq2 no 309 untr female vs. male
core_set<-d_up.reg.ITL_Up.bound.enh_activeprom.and.enh$Gene
core_set<-df_coreset$Gene

#write.table(core_set,"./core_set.txt",quote = F)
core_set_si<-d_up.enh_down.si$Gene
core_set_si<-df_coreset[which(df_coreset$log2FoldChange_siS.vs.siC_ITL<0&df_coreset$pvalue_siS.vs.siC_ITL<0.05),"Gene"]



mtsq_untr_a_genes_df<-mtsq_untr_a[mtsq_untr_a$gene_name%in%core_set,]
mtsq_untr_a_core_set_sis_df<-mtsq_untr_a[mtsq_untr_a$gene_name%in%core_set_si,]

mtsq_untr_a_sig<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05),]
mtsq_untr_a_sig_smc1a<-mtsq_untr_a_sig[mtsq_untr_a_sig$gene_name%in%core_set,] #110
mtsq_untr_a_sig_up<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male>1),] #596
mtsq_untr_a_sig_down<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male<(-1)),] #149

### volcano plot
#Transcriptome analysis revealed a total of 745 genes with sexual differential expression in SLE monocytes 
#(|log2-FC|>1, p-value<0.05) (Figure 6A). Female-overexpressed genes (n=596; Suppl) showed overrepresentation 
#of inflammatory pathways including cellular activation, inflammatory response regulation, response to LPS and
#MAPK cascade activation (Figure 6B). On the other hand, genes with male-biased expression (n=149; Suppl)
#were enriched in ribosome biogenesis, cellular metabolism and mitochondrial function (Figure 6C).

#mtsq_untr_a_sig_up<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male>0.58),] #1238
#mtsq_untr_a_sig_down<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male<(-0.58)),] #727


intersect(mtsq_untr_a_sig_up$gene_name,core_set_si) # [1] "GRAMD1A" "ZC3HAV1" "IL1A"    "ZBTB17"  "TRIM25"  "N4BP2L1" "PMAIP1"  "GBP5"    "ENGASE"  "TNFAIP2" "ADA"     "PELI1"  [13] "SLC9A8"
intersect(mtsq_untr_a_sig_up$gene_name,core_set) #88 ----61

intersect(mtsq_untr_a_sig_down$gene_name,core_set_si) #0
intersect(mtsq_untr_a_sig_down$gene_name,core_set) #3 ----2

colnames(mtsq_untr_a)
mtsq_untr<-mtsq_untr_a[,c(7,9,12)]
data_with_SLErnaseq<-merge(data,mtsq_untr,by.x="Gene",by.y="gene_name")
data_with_SLErnaseq_active_enh<-data_with_SLErnaseq[which(data_with_SLErnaseq$active_enh=="Yes"),]
colnames(data_with_SLErnaseq)

gg <- ggplot(data_with_SLErnaseq_active_enh, aes(x=log2FoldChange_siS.vs.siC_ITL, y=log2_normalized_fold_change_untr_female_vs_untr_male)) + 
  geom_point(aes(col=up_active_enhancers)) + 
  # geom_smooth(method="loess", se=F) + 
  #xlim(c(0, 0.1)) + 
  #ylim(c(0, 500000)) + 
  geom_jitter(aes(col=up_active_enhancers)) + 
  geom_smooth(aes(col=up_active_enhancers), method="lm", se=F)+
  stat_cor(aes(color = up_active_enhancers),method = "spearman")+
  #ylab("Expression in siSMC1a vs. siControl in ITL condition")+
  #xlab("Expression in SLE females vs. males")+
  theme(
    #plot.title = element_text(color = "black", size = 19),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))
 # labs(subtitle="Enhancer Status of Sex biased genes", 
    #   y="Expression in siSMC1a vs. siControl in ITL condition", 
   #    x="Expression in SLE females vs. males", 
   #    title="Correlation of expression in siSMC1a vs. siControl in ITL condition and SLE females vs. males")

plot(gg)
#ggsave("./Correlation of expression in siSMC1a vs. siControl in ITL condition and SLE females vs. males_spearman cor.png",gg,width = 10, height = 7, dpi = 300)
ggsave("final_fig/SupplFigure9_Correlation of expression in siSMC1a vs. siControl in ITL condition and SLE females vs. males_spearman cor.png",gg,width = 140, height = 110,units ="mm" ,dpi = 300)


gg <- ggplot(data_with_SLErnaseq_active_enh, aes(x=log2FoldChange_siS.vs.siC_ITL, y=log2_normalized_fold_change_untr_female_vs_untr_male)) + 
  geom_point(aes(col=up_active_enhancers)) + 
  # geom_smooth(method="loess", se=F) + 
  #xlim(c(0, 0.1)) + 
  #ylim(c(0, 500000)) + 
  geom_jitter(aes(col=up_active_enhancers)) + 
  geom_smooth(aes(col=up_active_enhancers), method="lm", se=F)+
  stat_cor(aes(color = up_active_enhancers),method = "pearson")+
  labs(subtitle="Enhancer Status of Sex biased genes", 
       y="Expression in siSMC1a vs. siControl in ITL condition", 
       x="Expression in SLE females vs. males", 
       title="Correlation of expression in siSMC1a vs. siControl in ITL condition and SLE females vs. males")

plot(gg)
#ggsave("./Correlation of expression in siSMC1a vs. siControl in ITL condition and SLE females vs. males_pearson cor.png",gg,width = 10, height = 7, dpi = 300)

#####
df_siun<-data[which(data$pvalue_siS.vs.siC_untr<0.05),]
df_siun$expr<-"no"
df_siun$expr[which(df_siun$log2FoldChange_siS.vs.siC_untr>0&df_siun$pvalue_siS.vs.siC_untr<0.05)]<-"up"
table(df_siun$expr,df_siun$up_enhancers)
fisher.test(table(df_siun$expr,df_siun$up_enhancers))

df_siITL<-data[which(data$pvalue_siS.vs.siC_ITL<0.05),]
df_siITL$expr<-"up"
df_siITL$expr[which(df_siITL$log2FoldChange_siS.vs.siC_ITL<0)]<-"down"
table(df_siITL$expr,df_siITL$up_enhancers)
fisher.test(table(df_siITL$expr,df_siITL$up_enhancers))
#

colonized_genes<-colonizedData$Gene
primed_genes<-primedData$Gene
response_and_colonized_genes<-intersect(responseData$Gene, colonizedData$Gene)

activated_genes<-response_and_colonized_genes #156
length(data[which(data$padjvalue_ITLvsC_siC<0.05&data$log2FoldChange_ITLvsC_siC>1&data$up_active_enhancers=="Yes"),"Gene"]) #279
length(data[which(data$padjvalue_ITLvsC_siC<0.05&data$log2FoldChange_ITLvsC_siC>1&data$up_active_promoters=="Yes"),"Gene"]) #35


unique(up_bound_at_annot$SYMBOL) #564
unique(down_bound_at_annot$SYMBOL) #227
unique(up_bound_h3_annot$SYMBOL) #440
unique(down_bound_h3_annot$SYMBOL) #269

intersect(unique(up_bound_at_annot$SYMBOL) , primed_genes) #31
intersect(unique(up_bound_at_annot$SYMBOL) , colonized_genes) #36


## HEATMAP WITH LOG FC
m<-heatmap.2(as.matrix((response_and_colonized_genes_logFC[, 2:3])), hclustfun = function(x) hclust(x,method = 'ward.D2'), Colv=NA)

mt<-as.hclust(m$rowDendrogram);cutree(mt, k=4)->cluster;table(cluster)
ccols<-c("firebrick4","darkorange2","gold2","olivedrab3","forestgreen","steelblue3","royalblue4","mediumorchid4")
cols<-colorRampPalette(c("steelblue4","white", "firebrick4"))(100)
#png("~/Dropbox/Articles/collaborations/Bertsias_Kosmara_SMC1A/Figures/PrimedColonizedHeatmap.png", width=480, height=960)
heatmap.2(as.matrix((response_and_colonized_genes_logFC[,2:3])), hclustfun = function(x) hclust(x,method = 'ward.D2'), Colv=NA, labRow=T,
          labCol=c("Response","Colonized"), col=cols, RowSideColors=ccols[cluster], 
          breaks = seq(-1.5, 1.5, length.out = 101),
          key.title="log2(FC)", main="Gene Expression Changes\nin SMC1A Down-regulation", density.info="none", scale="none", mar=c(7,1), tracecol="black", vline=0, cexCol = 1.5)
#
m <- heatmap.2(as.matrix(response_and_colonized_genes_logFC[, 2:3]), 
               hclustfun = function(x) hclust(x, method = 'ward.D2'), 
               Colv = NA)

mt <- as.hclust(m$rowDendrogram)
cutree(mt, k = 4) -> cluster
table(cluster)

ccols <- c("firebrick4", "darkorange2", "gold2", "olivedrab3", "forestgreen", "steelblue3", "royalblue4", "mediumorchid4")
cols <- colorRampPalette(c("steelblue4", "white", "firebrick4"))(100)

heatmap.2(as.matrix(response_and_colonized_genes_logFC[, 2:3]), 
          labRow = TRUE,  # Show row labels
          labCol = c("Response", "Colonized"), 
          col = cols, 
          breaks = seq(-1.5, 1.5, length.out = 101),
          key.title = "log2(FC)", 
          main = "Gene Expression Changes\nin SMC1A Down-regulation", 
          density.info = "none", 
          scale = "none", 
          mar = c(7, 1), 
          tracecol = "black", 
          vline = 0, 
          cexCol = 1.5)

p<-pheatmap(
  as.matrix(response_and_colonized_genes_logFC[, 2:3]),
  labels_col = c("Response", "Colonized"), 
  color = cols,                   # Color palette
  breaks = seq(-1.5, 1.5, length.out = 101),
  key = TRUE,                     # Show color key
  key.title = "log2(FC)", 
  main = "Gene Expression Changes\nin Colonized and Response Genes", 
  cluster_rows = TRUE,            # Cluster rows
  cluster_cols = FALSE,           # Do not cluster columns
  scale = "none",                 # Do not scale
  fontsize_row = 8,               # Row label font size
  fontsize_col = 12 ,              # Column label font size
  angle_col = 0
)
#ggsave("./heatmap 156 genes.png",p,width = 8, height = 16, dpi = 120)


#wtITLgenes <- arrange(data[which( (abs(data$log2FoldChange_ITLvsC_siC)>=1) & (data$padjvalue_ITLvsC_siC<=0.05) ), c(1,4,5)], log2FoldChange_ITLvsC_siC)
#siITLgenes <- arrange(data[which( (abs(data$log2FoldChange_ITLvsC_siS)>=1) & (data$padjvalue_ITLvsC_siS<=0.05) ), c(1,2,3)], log2FoldChange_ITLvsC_siS)
#commonITL <- intersect(wtITLgenes$Gene, siITLgenes$Gene)
#responseData<-data[which(data$Gene %in% commonITL),]
responseData<-data[which( (abs(data$log2FoldChange_ITLvsC_siC)>=1) & (data$padjvalue_ITLvsC_siC<=0.05) ),]

ego_responseData <- enrichGO(gene = unique(responseData[,1]),
                             keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05, readable = FALSE)
d<-dotplot(ego_responseData, showCategory=15,orderBy = "x", title = "Enriched terms treated vs. untreated deregulated genes")

ego_responseData_SIMPL<-simplify(ego_responseData,cutoff = 0.7,semData=NULL)

d_s<-dotplot(ego_responseData_SIMPL, showCategory=15,orderBy = "x", title = "Enriched terms treated vs. untreated deregulated genes")

response_and_colonized_genes<-intersect(responseData$Gene, colonizedData$Gene) #156
response_and_colonized_genes_df<-data[which(data$Gene%in%response_and_colonized_genes),]
#write.table(response_and_colonized_genes_df,"response_colonized_genes.txt",sep="\t",quote = F,row.names = F,col.names = T)
response_and_colonized_genes_logFC<-data[which(data$Gene%in%response_and_colonized_genes),c("Gene","log2FoldChange_ITLvsC_siC","log2FoldChange_siS.vs.siC_ITL","log2FoldChange_siS.vs.siC_untr")]
library(pheatmap)
rownames(response_and_colonized_genes_logFC)<-response_and_colonized_genes_logFC$Gene
p<-pheatmap(as.matrix((response_and_colonized_genes_logFC[, 2:4])) , scale="row", color=colorRampPalette(c("#1874CD", "white", "#EE0000"))(100), show_rownames = T,
            clustering_method = "ward.D2", angle_col = 45, fontsize_col = 12, cluster_cols=F,     
            show_colnames = T, border_color = NA)
p<-pheatmap(as.matrix((response_and_colonized_genes_logFC[, 2:4])) , scale="none", color=colorRampPalette(c("#1874CD", "white", "#EE0000"))(100), show_rownames = T,
            clustering_method = "ward.D2", angle_col = 45, fontsize_col = 12, cluster_cols=F,     
            show_colnames = T, border_color = NA)

###
# Create a color palette with white at 0
color_palette <- colorRampPalette(c("#1874CD", "white", "#EE0000"), space = "rgb")(100)

# Set 0 as white in the color palette
zero_index <- round((0 - min(response_and_colonized_genes_logFC[, 2:4])) / (max(response_and_colonized_genes_logFC[, 2:4]) - min(response_and_colonized_genes_logFC[, 2:4])) * 100)
color_palette[zero_index] <- "#FFFFFF"  # Set 0 as white

# Create the heatmap
custom_color_palette <- function(n) {
  colors <- colorRampPalette(c("#1874CD", "#FFFFFF", "#EE0000"))(n)
  return(colors)
}
p <- pheatmap(
  as.matrix(response_and_colonized_genes_logFC[, 2:4]),
  scale = "none",
  color = custom_color_palette(100),
  show_rownames = TRUE,
  clustering_method = "ward.D2",
  angle_col = 45,
  fontsize_col = 12,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  border_color = NA)
##
m<-heatmap.2(as.matrix((response_and_colonized_genes_logFC[, 2:4])), hclustfun = function(x) hclust(x,method = 'ward.D2'), Colv=NA)

mt<-as.hclust(m$rowDendrogram);cutree(mt, k=4)->cluster;table(cluster)
ccols<-c("firebrick4","darkorange2","gold2","olivedrab3","forestgreen","steelblue3","royalblue4","mediumorchid4")
cols<-colorRampPalette(c("steelblue4","white", "firebrick4"))(256)
#png("~/Dropbox/Articles/collaborations/Bertsias_Kosmara_SMC1A/Figures/PrimedColonizedHeatmap.png", width=480, height=960)
heatmap.2(as.matrix((response_and_colonized_genes_logFC[,2:4])), hclustfun = function(x) hclust(x,method = 'ward.D2'), Colv=NA,
          labCol=c("Response","Colonized", "Primed"), col=cols, RowSideColors=ccols[cluster], 
          key.title="log2(FC)", main="Gene Expression Changes\nin SMC1A Down-regulation", density.info="none", scale="row", mar=c(7,1), tracecol="black", vline=0, labRow="", cexCol = 1.5)
#dev.off()
colors = c(seq(-1,-0.01,length=100),seq(-0.009,0.009,length=100),seq(0.01,1,length=100))

my_palette <- colorRampPalette(c("blue","red"))(n = 299)
heatmap.2(as.matrix((response_and_colonized_genes_logFC[,2:4])), hclustfun = function(x) hclust(x,method = 'ward.D2'), Colv=NA,
          labCol=c("Response","Colonized", "Primed"), breaks=colors,col=my_palette, RowSideColors=ccols[cluster], symm=F,symkey=F,symbreaks=T,
          key.title="log2(FC)", main="Gene Expression Changes\nin SMC1A Down-regulation", density.info="none", scale="none", mar=c(7,1), tracecol="black", vline=0, labRow="", cexCol = 1.5)

m<-heatmap.2(as.matrix((response_and_colonized_genes_logFC[, 2:3])), hclustfun = function(x) hclust(x,method = 'ward.D2'), Colv=NA)

mt<-as.hclust(m$rowDendrogram);cutree(mt, k=4)->cluster;table(cluster)

heatmap.2(as.matrix((response_and_colonized_genes_logFC[,2:3])), hclustfun = function(x) hclust(x,method = 'ward.D2'), Colv=NA,
          labCol=c("Response","Colonized", "Primed"), breaks=colors,col=my_palette, RowSideColors=ccols[cluster], symm=F,symkey=F,symbreaks=T,
          key.title="log2(FC)", main="Gene Expression Changes\nin SMC1A Down-regulation", density.info="none", scale="none", mar=c(7,1), tracecol="black", vline=0, labRow="", cexCol = 1.5)

response_and_colonized_genes_logFC_sign<-response_and_colonized_genes_logFC
response_and_colonized_genes_logFC_sign$sign_si<-response_and_colonized_genes_logFC_sign$log2FoldChange_siS.vs.siC_ITL*response_and_colonized_genes_logFC_sign$log2FoldChange_siS.vs.siC_untr



############ plot for figure s5
# itl stimulation for degs si treated and untreated
itlgenes<-arrange(data[which(data$pvalue_siS.vs.siC_ITL<=0.05&data$log2FoldChange_siS.vs.siC_ITL<0), c(1,4,5)], log2FoldChange_ITLvsC_siC)
colnames(itlgenes)<-c("Gene", "log2FC", "padj")
itlgenes$expr[itlgenes$log2FC>0]<-"up"
itlgenes$expr[itlgenes$log2FC<0]<-"down"
itlgenes$expr<-as.factor(itlgenes$expr)
#
itlgenes$delabel <- NA
itlgenes$delabel[which(abs(itlgenes$log2FC)>2&itlgenes$padj<0.05)]<-as.character(itlgenes$Gene[which(abs(itlgenes$log2FC)>2&itlgenes$padj<0.05)])
#
g<-ggplot(itlgenes, aes(x=log2FC, y=-log10(padj), color=expr, label=delabel)) + 
  geom_point() + 
  scale_color_manual(values=c("steelblue4", "firebrick4")) +
  geom_text_repel(box.padding = 0.3, max.overlaps = Inf, col="black", size=3) +
  labs(title="Down-regulated genes in si-SMC1A vs. si-Control treated Monocytes", subtitle="Gene Expression Changes in ITL stimulation (275 Genes at p<=0.05)", y="-log10(p-adjusted)", x="log2FC", size="|log2FC|", color="Expression")
#ggsave("./Volcano - Down-regulated genes in si-SMC1A vs. si-Control treated Monocytes.png",g,width = 10, height = 8, dpi = 300)


itlgenes<-arrange(data[which(data$pvalue_siS.vs.siC_untr<=0.05&data$log2FoldChange_siS.vs.siC_untr<0), c(1,4,5)], log2FoldChange_ITLvsC_siC)
colnames(itlgenes)<-c("Gene", "log2FC", "padj")
itlgenes$expr[itlgenes$log2FC>0]<-"up"
itlgenes$expr[itlgenes$log2FC<0]<-"down"
itlgenes$expr<-as.factor(itlgenes$expr)
#
itlgenes$delabel <- NA
itlgenes$delabel[which(abs(itlgenes$log2FC)>2&itlgenes$padj<0.05)]<-as.character(itlgenes$Gene[which(abs(itlgenes$log2FC)>2&itlgenes$padj<0.05)])
#
g<-ggplot(itlgenes, aes(x=log2FC, y=-log10(padj), color=expr, label=delabel)) + 
  geom_point() + 
  scale_color_manual(values=c("steelblue4", "firebrick4")) +
  geom_text_repel(box.padding = 0.3, max.overlaps = Inf, col="black", size=3) +
  labs(title="Down-regulated genes in si-SMC1A vs. si-Control untreated Monocytes", subtitle="Gene Expression Changes in ITL stimulation (223 Genes at p<=0.05)", y="-log10(p-adjusted)", x="log2FC", size="|log2FC|", color="Expression")
#ggsave("./Volcano - Down-regulated genes in si-SMC1A vs. si-Control untreated Monocytes.png",g,width = 10, height = 8, dpi = 300)



###
#Here we can perform a Fisher's test taking into account the following populations:
#a. Genes with active enhancers
#b. Genes without active enhancers (complementary)
#vs
#c. Genes with up-bound SMC1A in enhancers
#d. Genes without up-bound SMC1A in enhancers (complementary) 

table(data$active_enh)
fisher.test(table(data$active_enh, data$up_enhancers))
#p-value < 2.2e-16, odds ratio  54.03537
# No  Yes
#No  7249   29
#Yes 3908  845


fisher.test(table(data$active_prom, data$up_promoters))
#p-value = 0.1047, odds ratio = 1.51763 
# No  Yes
#No  2670   20
#Yes 9236  105




###
#dim(data)
#[1] 12031    21

### numbers for the paper
length(unique(data[which(data$active_prom=="Yes"&data$active_enh=="Yes"),1])) #4221
length(unique(data[which(data$active_prom=="Yes"|data$active_enh=="Yes"),1])) #9834
length(unique(data[which(data$active_enh=="Yes"),1])) #4737
length(unique(data[which((data$active_prom=="Yes"|data$active_enh=="Yes")&(data$EnhancerStatus=="Up"|data$EnhancerStatus=="UpActive")),1])) #851
fisher.test(table(data$active_prom=="Yes"|data$active_enh=="Yes", data$EnhancerStatus=="Up"|data$EnhancerStatus=="UpActive"))
length(unique(data[which((data$active_prom=="Yes"|data$active_enh=="Yes")&(data$PromoterStatus=="Down"|data$PromoterStatus=="DownActive")),1])) #559
fisher.test(table(data$active_prom=="Yes"|data$active_enh=="Yes", data$PromoterStatus=="Down"|data$PromoterStatus=="DownActive"))
fisher.test(table(data$active_prom=="Yes"|data$active_enh=="Yes", data$PromoterStatus=="Up"|data$PromoterStatus=="UpActive"))

#SMC1A binding at active genome regulatory elements induces the transcriptional activation 
#of inflammatory genes in SLE-like monocytes
library(tidyverse)
data_respon<-data[which(abs(data$log2FoldChange_ITLvsC_siC)>1&data$padjvalue_ITLvsC_siC<0.05),]
data_respon<-data_respon %>% distinct(Gene, .keep_all = TRUE) #2778
length(which((data_respon$EnhancerStatus=="UpActive"&data_respon$PromoterStatus=="UpActive"))) #24

data_activ<-data[which(data$log2FoldChange_ITLvsC_siC>1&data$padjvalue_ITLvsC_siC<0.05),]
data_activ<-data_activ %>% distinct(Gene, .keep_all = TRUE) #1521
length(which((data_activ$EnhancerStatus=="UpActive"&data_activ$PromoterStatus=="UpActive"))) #24

length(unique(data_activ[which((data_activ$active_prom=="Yes"|data_activ$active_enh=="Yes")&(data_activ$EnhancerStatus=="Up"|data_activ$EnhancerStatus=="UpActive")),1])) #301
length(unique(data_respon[which((data_respon$active_prom=="Yes"|data_respon$active_enh=="Yes")&(data_respon$EnhancerStatus=="Up"|data_respon$EnhancerStatus=="UpActive")),1])) #358

df_coreset<-data_activ[which((data_activ$active_prom=="Yes"|data_activ$active_enh=="Yes")&(data_activ$EnhancerStatus=="Up"|data_activ$EnhancerStatus=="UpActive")),] #301
df_coreset[which(df_coreset$log2FoldChange_siS.vs.siC_ITL<0&df_coreset$pvalue_siS.vs.siC_ITL<0.05),1]

length(unique(data_activ[which((data_activ$EnhancerStatus=="None"|data_activ$PromoterStatus=="None")&(data_activ$up_enhancers=="No"|data_activ$up_promoters=="No")),1])) #723
length(unique(data_respon[which((data_respon$EnhancerStatus=="None"|data_respon$PromoterStatus=="None")&(data_respon$up_enhancers=="No"|data_respon$up_promoters=="No")),1])) #1532


length(which((data$active_prom=="Yes"&data$active_enh=="Yes")&data$up_enhancers=="Yes")) #772

length(which((data$active_prom=="Yes"|data$active_enh=="Yes")&data$up_enhancers=="Yes")) #9857
d_up.enh<-data[which(data$up_enhancers=="Yes"),] #874
#up_bound_enhancers -> 874/12031

d_down.enh<-data[which(data$down_enhancers=="Yes"),] #285
intersect(d_up.enh$Gene,d_down.enh$Gene) #16

#1% of the genes (n=125) showed increased SMC1A binding and 4.9% (n=585 genes)
d_up.prom<-data[which(data$up_promoters=="Yes"),] #125
d_down.prom<-data[which(data$down_promoters=="Yes"),] #585

d_act.enh<-data[which(data$active_enh=="Yes"),] #4753
d_act.prom<-data[which(data$active_prom=="Yes"),] #9341
intersect(d_act.enh$Gene,d_act.prom$Gene)#4221
#4753+9341-4221 #9873

d_up.act.enh<-data[which(data$up_active_enhancer=="Yes"),] #747
#12031/747
length(d_up.act.enh[which(d_up.act.enh$log2FoldChange_ITLvsC_siC>1&d_up.act.enh$padjvalue_ITLvsC_siC<0.01),1])

d_up.enh_act.enh_or_act.prom<-data[which((data$active_enh=="Yes"|data$active_prom=="Yes")&data$up_active_enhancers=="Yes"),] #747
### group b
d_up.enh_act.enh_and_act.prom<-data[which((data$active_enh=="Yes"&data$active_prom=="Yes")&data$up_active_enhancers=="Yes"),] #685
d_down.act.prom<-data[which(data$down_active_promoters=="Yes"),] #545
#545/9873

d_no.bind_no.act<-data[which(data$EnhancerStatus=="None"&data$EnhancerStatus=="None"),] #7233

d_up.enh_act.enh_and_act.prom_up.itl<-d_up.enh_act.enh_and_act.prom[which(d_up.enh_act.enh_and_act.prom$padjvalue_ITLvsC_siC<0.05&d_up.enh_act.enh_and_act.prom$log2FoldChange_ITLvsC_siC>1),]
### group c
d_up_prom_up_enh<-data[which(data$up_active_enhancers=="Yes"&data$up_active_promoters=="Yes"),] #35
## setdiff with enh
setdiff(d_up.act.enh$Gene,d_up_prom_up_enh$Gene) #712
### group b + c (685, 35)
setdiff(d_up.enh_act.enh_and_act.prom$Gene,d_up_prom_up_enh$Gene) #650

d_up.prom_up.enh_down.si<-data[which(data$up_active_enhancers=="Yes"&data$up_active_promoters=="Yes"&data$pvalue_siS.vs.siC_ITL<0.05&data$log2FoldChange_siS.vs.siC_ITL<0),] #35
d_up.prom_up.enh_up.itl<-d_up_prom_up_enh[which(d_up_prom_up_enh$padjvalue_ITLvsC_siC<0.01&d_up_prom_up_enh$log2FoldChange_ITLvsC_siC>1),] #24

d_up_enh_activeprom.and.enh<-data[which(data$up_active_enhancers=="Yes"&data$active_prom=="Yes"),] #685

## CORE SET OF GENES
d_up.reg.ITL_Up.bound.enh_activeprom.and.enh<-data[which(data$up_active_enhancers=="Yes"&data$active_prom=="Yes"&data$padjvalue_ITLvsC_siC<0.01&data$log2FoldChange_ITLvsC_siC>=1),] #253
#d_up.reg.ITL_Up.bound.enh_activeprom.and.enh<-data[which(data$up_active_enhancers=="Yes"&data$padjvalue_ITLvsC_siC<0.01&data$log2FoldChange_ITLvsC_siC>1),] #
out <- enrichGO(gene = as.character(d_up.reg.ITL_Up.bound.enh_activeprom.and.enh$Gene), keyType = "SYMBOL", OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = FALSE); 
out_SIMPL<-simplify(out,cutoff = 0.7,semData=NULL)
d<-dotplot(out_SIMPL, showCategory=15) + ggtitle("Up bound by SMC1a genes in enhancer region \n Up regulated after IFN-TNF-LPS stimulation")
#ggsave("./DotplotUp bound by SMC1a genes in enhancer region \n Up regulated after IFN-TNF-LPS.png",d,width = 8, height = 8, dpi = 300)

d_up.enh_down.si<-d_up.reg.ITL_Up.bound.enh_activeprom.and.enh[which(d_up.reg.ITL_Up.bound.enh_activeprom.and.enh$pvalue_siS.vs.siC_ITL<0.05&d_up.reg.ITL_Up.bound.enh_activeprom.and.enh$log2FoldChange_siS.vs.siC_ITL<0),] #23




##### number of genes - 5/7/24
#Group B: should fulfill all the criteria below:
#1. Responsive genes  (up- or down-regulated, we don’t mind)
##2. Active enhancer: YES
#3. Active promoter: YES
#4. Up-bound SMC1A at the enhancer: YES
#5. Up-bound SMC1A at the promoter: NO

responseData<-data[which( (abs(data$log2FoldChange_ITLvsC_siC)>=1) & (data$padjvalue_ITLvsC_siC<=0.05) ),]
length(unique(responseData$Gene)) #2778
### group b
d_up.enh_act.enh_and_act.prom<-responseData[which((responseData$active_enh=="Yes"&responseData$active_prom=="Yes")&responseData$up_active_enhancers=="Yes"),] #685
d_up.enh_no.up.prom_act.enh_and_act.prom<-responseData[which(responseData$active_enh=="Yes"&responseData$active_prom=="Yes"&responseData$up_enhancers=="Yes"&responseData$up_promoters=="No"),] #301

length(unique(d_up.enh_act.enh_and_act.prom$Gene)) #293
length(unique(d_up.enh_no.up.prom_act.enh_and_act.prom$Gene)) #301
groupb<-unique(d_up.enh_no.up.prom_act.enh_and_act.prom$Gene)


#Group C: should fulfill all the criteria below:
#1. Responsive genes  (up- or down-regulated, we don’t mind)
#2. Active enhancer: YES
#3. Active promoter: YES
#4. Up-bound SMC1A at the enhancer: YES
#5. Up-bound SMC1A at the promoter: YES


### group c
d_up_prom_up_enh<-responseData[which(responseData$up_active_enhancers=="Yes"&responseData$up_active_promoters=="Yes"),] #35
length(unique(d_up_prom_up_enh$Gene)) #24
groupc<-unique(d_up_prom_up_enh$Gene)
#Group A: 
#1. Responsive genes  (up- or down-regulated, we don’t mind)
#2. Active enhancer: ΝΟ
#3. Active promoter: ΝΟ
#4. Up-bound SMC1A at the enhancer: ΝΟ
#5. Up-bound SMC1A at the promoter: ΝΟ

d_no.up.enh_no.up.prom_no.act.enh_and_no.act.prom<-responseData[which(responseData$active_enh=="No"&responseData$active_prom=="No"&responseData$up_enhancers=="No"&responseData$up_promoters=="No"),] #685
length(unique(d_no.up.enh_no.up.prom_no.act.enh_and_no.act.prom$Gene)) #530
groupa<-unique(d_no.up.enh_no.up.prom_no.act.enh_and_no.act.prom$Gene)
#SMC1A “target genes”: should fulfill all the criteria below
#1. Group B  OR  Group C genes
#2. log2FC ≥1
#3. Adj. p-value <0.05

u<-union(d_up.enh_no.up.prom_act.enh_and_act.prom$Gene,d_up_prom_up_enh$Gene) #325
r<-responseData[which(responseData$log2FoldChange_ITLvsC_siC>1),"Gene"]
i<-intersect(r,u) #277
targetgenes<-unique(i)
data_si<-data[which(data$pvalue_siS.vs.siC_ITL<0.05&data$log2FoldChange_siS.vs.siC_ITL<0),]
intersect(i,data_si$Gene)


#Group B: 301
#Group C: 24
#Group A: 530
#“Target”: 277
#write.table(groupb,"groupa.txt",quote = F,row.names = F,col.names = F)
#write.table(groupb,"groupb.txt",quote = F,row.names = F,col.names = F)
#write.table(groupc,"groupc.txt",quote = F,row.names = F,col.names = F)




## deseq2 no309 untr female vs. male
mtsq_untr_a<-read.table("metaseqr_all_out_untr_female_vs_untr_male.txt", header = T)
mtsq_untr_a_genes_df<-mtsq_untr_a[mtsq_untr_a$gene_name%in%targetgenes,]
mtsq_untr_a_sig<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05),]
mtsq_untr_a_sig_smc1a<-mtsq_untr_a_sig[mtsq_untr_a_sig$gene_name%in%targetgenes,] #11 "GTPBP1"  "TRIM25"  "N4BP2L1" "PML"     "PMAIP1"  "GBP5"    "PARP14"  "TNFAIP2" "ADA"     "PELI1"   "SLC9A8" 
mtsq_untr_a_sig_up<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male>0),] #1664
mtsq_untr_a_sig_down<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.05&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male<0),] #1315

mtsq_untr_a_sig_up<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.01&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male>0.58),] #836
mtsq_untr_a_sig_down<-mtsq_untr_a[which(mtsq_untr_a$p.value_deseq2<0.01&mtsq_untr_a$log2_normalized_fold_change_untr_female_vs_untr_male<(-0.58)),] #330
## deseq2 all LPS female vs. male
mtsq_lps_a<-read.table("metaseqr_all_out_LPS_female_vs_LPS_male.txt", header = T)
mtsq_lps_a_genes_df<-mtsq_lps_a[mtsq_lps_a$gene_name%in%targetgenes,]
mtsq_lps_a_sig<-mtsq_lps_a[which(mtsq_lps_a$p.value_deseq2<0.05),]
mtsq_lps_a_sig_smc1a<-mtsq_lps_a_sig[mtsq_lps_a_sig$gene_name%in%targetgenes,] #0
mtsq_lps_a_sig_up<-mtsq_lps_a[which(mtsq_lps_a$p.value_deseq2<0.05&mtsq_lps_a$log2_normalized_fold_change_LPS_female_vs_LPS_male>0),] #1671
mtsq_lps_a_sig_down<-mtsq_lps_a[which(mtsq_lps_a$p.value_deseq2<0.05&mtsq_lps_a$log2_normalized_fold_change_LPS_female_vs_LPS_male<0),] #591

## common degs untr - LPS
common_u_lps_up<-intersect(mtsq_untr_a_sig_up$gene_name,mtsq_lps_a_sig_up$gene_name) #659
common_u_lps_down<-intersect(mtsq_untr_a_sig_down$gene_name,mtsq_lps_a_sig_down$gene_name) #274

## induced by LPS
up_lps_uniq<-setdiff(mtsq_lps_a_sig_up$gene_name,mtsq_untr_a_sig_up$gene_name) #1012
down_lps_uniq<-setdiff(mtsq_lps_a_sig_down$gene_name,mtsq_untr_a_sig_down$gene_name) #316

#### common genes (212 and quant seq untr deseq2 with outliers)
common_smc1a_untr_a<-intersect(mtsq_untr_a_sig$gene_name, targetgenes) #102
mtsq_untr_a_smc1a<-mtsq_untr_a_sig[which(mtsq_untr_a_sig$gene_name%in%common_smc1a_untr_a),] # only 5/102 down

up_females_untr_smc1a_reg<-mtsq_untr_a_smc1a[which(mtsq_untr_a_smc1a$log2_normalized_fold_change_untr_female_vs_untr_male>0),"gene_name"]
#write.table(up_females_untr_smc1a_reg,"~/SMC1A_project/quant_seq_SLE_females-males/deseq_no309/up_females_untr_smc1a_reg.txt", row.names = F, col.names = F, quote = F)
ego_common_smc1a_untr_a <- enrichGO(gene = mtsq_untr_a_smc1a[which(mtsq_untr_a_smc1a$log2_normalized_fold_change_untr_female_vs_untr_male>0),"gene_name"],
                                    keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05, readable = FALSE)
#write.table(as.data.frame(ego_common_smc1a_untr_a),"~/SMC1A_project/quant_seq_SLE_females-males/deseq_no309/ego_up_females_untr_smc1a_reg.xlsx",row.names = F, col.names = T)
d<-dotplot(ego_common_smc1a_untr_a, showCategory=30,orderBy = "x", title = "Enriched terms untreated female vs. male up regulated genes and \nSMC1a target genes")
#ggsave("~/SMC1A_project/quant_seq_SLE_females-males/deseq_no309/Enriched terms untreated female vs. male up regulated genes and SMC1a up bound and active enhancer associated.png",d,width = 8, height = 12)

ego_common_smc1a_untr_a_SIMPL<-clusterProfiler::simplify(ego_common_smc1a_untr_a, cutoff = 0.7, semData=NULL)
d<-dotplot(ego_common_smc1a_untr_a_SIMPL, showCategory=15) + ggtitle("Gene Ontology Enrichment Analysis of \nfemale vs. male up regulated and SMC1a target genes")
#ggsave("Gene Ontology Enrichment Analysis of female vs. male up regulated and SMC1a target genes.png",d,width = 8, height = 7)


