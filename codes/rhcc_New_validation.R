rm(list = ls())
gc()

library(Seurat)
library(tidyverse)
library(scRepertoire)
library(openxlsx)
library(circlize)
library(scales)
library(Hmisc)
library(shazam)
library(data.table)
library(languageserver)
library(rstatix)
library(ggpubr)
library(ggforce)
library(ggalluvial)
library(RColorBrewer)
library(patchwork)
library(Hmisc)
library(ggsci)
library(pheatmap)

#####Visualization Parameters#####
## Set the theme for the plots
## ggthemes:
mytheme2 = theme(
  plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
  plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
  plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
  axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
  axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
  axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
  axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
  legend.title  = element_text(color = 'black', size  = 16),
  legend.text   = element_text(color = 'black', size   = 16),
  axis.line.y = element_line(color = 'black', linetype = 'solid'), 
  axis.line.x = element_line (color = 'black',linetype = 'solid'), 
  panel.background=element_rect(fill="white")
)

## Colors
color = brewer.pal(12,"Set3")
two.colors = c("#87E9C3","#FCCF7C")
three.colors = c()
four.colors = c()
five.colors = c()
color.s1 = c("#E64B35B2","#00A087B2","#3C5488B2","#DC0000B2","#F39B7FB2",
             "#8491B4B2","#4DBBD5B2","#91D1C2B2","#7E6148B2")


## linewidth
violin.box.lwd = .1
p.lwd = .1

source("tools/utils_plots.R")

#####WYQ colors#####
library(scales)
#color4=c("#F8C4AC",'#DD8385','#C4C1DE','#9A98C9')
color4=c("#fdc7cd",'#ee84a8','#c2d7f3','#778ccc')

my.colors<-c('#BF331D','#199198', '#EADCD3','#E4C755','#8FCBD9','#736294','#D65792','#76575D',
             '#FACFCE','#6D9891','#F69F83','#2E4E50', #T
             '#58A544','#EDAA25','#D5D3D4', #B
             '#40498F','#931D54','#CDDEAA','#3D6584','#D79B8B','#9393AC',
             '#94B4BD','#A3A095','#C62C21','#954232', #Myeloid
             '#AE9E8A','#2C292B','#C962A0','#8E66A1', #Endo
             '#5C96BB','#87693F','#4A62A3','#D97360', #mac
             '#67AE84','#9FA8AD','#EF4C6A',
             '#E95C59',"#20B2AA",'#E5D2DD',"#9370DB","#98FB98", '#53A85F','#00468BFF', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
             '#AD002AFF', '#E59CC4',"#808000", '#23452F', '#BD956A', '#8C549C', '#585658',
             '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#F7F398',
             '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
             '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
             '#968175', 'seashell','wheat',"#1F78B4", "#33A02C", "#fcd300","#FB9A99","#FF7F00", "#CAB2D6",'#0780cf',
             '#0780cf','#765005','blue', 'orange', 'green', 'yellow', 'gray',
             'darkred','palevioletred','darkmagenta','rosybrown',
             'darkkhaki','darkslategray','limegreen','greenyellow',
             'darkblue','dodgerblue','mediumturquoise','dimgray','darkkhaki','limegreen')
my.colors1<-c("#E64B35B2","#00A087B2","#3C5488B2","#4DBBD5B2"  ,"#F39B7FB2" ,
              "#8491B4B2","#91D1C2B2" ,"#DC0000B2" ,"#7E6148B2")

#####Other Tool Functions#####
source("tools/st_interaction_enrich_function.R")
source("tools/dimplot_inhouse_R.R")
source("tools/featureplot_inhouse_R.R")
source("tools/get_surround_bins.R")
source("tools/compare_tertiary_bcr_and_single_cell_bcr.R")
source("tools/compare_tertiary_bcr_and_single_cell_bcr_extract_cdrs.R")
source("tools/cjy_BGI_interaction_tools.R")
source("tools/utils_plots.R")
source("tools/cjy_single_cell_proportion_test.R")
source("tools/as_matrix.R")
source("cjy_MSA_getCirclized_from_filtered_contig.R")
source("tools/featureplot_inhouse_R.R")
source("tools/cjy_BGI_interaction_tools.R")
source("tools/cjy_MSA_getCirclized_from_filtered_contig.R")
Spatial_interaction_spatial_featureplot <- function(ligand_gene, receptor_gene, stRNA, max_distance = 1,
                                                    prefix, features) {
  # Extract the bins that the expression of ligand_gene larger than 0
  ligand_gene_bins = colnames(stRNA@assays$Spatial@counts)[which(stRNA@assays$Spatial@counts[ligand_gene,] > 0)]
  meta.data = stRNA@meta.data
  
  # Extract all coordinates of the bins
  all_coords = meta.data[,c("row","col")]
  all_coords$bin = rownames(all_coords)
  
  # Initialize interaction_product
  interaction_product = numeric(length = ncol(stRNA))
  names(interaction_product) = colnames(stRNA)
  
  # Acquire the expression value of ligand_gene in ligand_gene_bins
  ligand_gene_values = stRNA@assays$Spatial@counts[ligand_gene, ligand_gene_bins]
  
  for (i in seq_along(ligand_gene_bins)) {
    bin = ligand_gene_bins[i]
    bin_expr = ligand_gene_values[i]
    bin_coords = meta.data[bin, c("row","col"), drop = F]
    
    surrounding_coords = get_surrounding_coords(bin_coords, max_distance = max_distance)
    sur_coords_bins = merge(all_coords, surrounding_coords, by = c("row","col"))
    
    sur_receptor_expr = stRNA@assays$Spatial@counts[receptor_gene, sur_coords_bins$bin]
    
    interaction_product[bin] = sum(bin_expr * sur_receptor_expr)
  }
  # Add the interaction_product to the counts slot
  count_mat = stRNA@assays$Spatial@counts
  new_gene_name = paste0(ligand_gene, "_", receptor_gene)
  # Make sure that the new_gene_name is not in the count_mat
  if (new_gene_name %in% rownames(count_mat)) {
    warning(paste("Gene", new_gene_name, "already exists in the count matrix. Overwriting it..."))
    count_mat = count_mat[-which(rownames(count_mat) == new_gene_name),]
  }
  # Convert the interaction_product to a sparse matrix
  interaction_product_sparse = Matrix::Matrix(interaction_product, nrow = 1, sparse = TRUE)
  rownames(interaction_product_sparse) = new_gene_name
  colnames(interaction_product_sparse) = colnames(count_mat)
  
  count_mat_extended = rbind(count_mat, interaction_product_sparse)
  stRNA@assays$Spatial@counts = count_mat_extended
  
  prefix = prefix
  features = features
  featureplot_inhouse_R(obj = stRNA, feature = features, prefix = prefix, name = features)
}


#####Figure 2#####
######Fig. 2b######
Tcell_TCR = readRDS("TCR_group/Tcell_TCR.rds")
Expansion_group=Tcell_TCR@meta.data[,c("orig.ident",'cellbarcode',"Frequency","CTnt","uniqueCTnt","treatment")]
TCR_pair=read.csv("TCR_group/RHCC_tumor_treatmentpair.csv")
Expansion_group2=merge(Expansion_group,TCR_pair,by.x = "orig.ident")
needed=unique(Expansion_group2$pair)
needed=needed[-c(13)]
Expansion_group2=subset(Expansion_group2,pair %in% needed)
load("TCR_group/Expansion_group3_Post.RData")

expanded.ct2 = Expansion_group3_Post$CTnt[Expansion_group3_Post$Freq.2.group == "yes" & !is.na(Expansion_group3_Post$Freq.2.group)]
expanded.ct.df2 = Expansion_group[Expansion_group$CTnt %in% expanded.ct2,]
expanded.ct.df2 = expanded.ct.df2[grep("A",expanded.ct.df2$orig.ident),]
Expansion_Tcell_count2 = expanded.ct.df2 %>% group_by(orig.ident) %>% count()
names(Expansion_Tcell_count2)[2] = "expanded.2.clone.count"

expanded.ct5 = Expansion_group3_Post$CTnt[Expansion_group3_Post$Freq.5.group == "yes" & !is.na(Expansion_group3_Post$Freq.5.group)]
expanded.ct.df5 = Expansion_group[Expansion_group$CTnt %in% expanded.ct5,]
expanded.ct.df5 = expanded.ct.df5[grep("A",expanded.ct.df5$orig.ident),]
Expansion_Tcell_count5 = expanded.ct.df5 %>% group_by(orig.ident) %>% count()
names(Expansion_Tcell_count5)[2] = "expanded.5.clone.count"

# 
# Expansion_Tcell_count = merge(Expansion_Tcell_count2, Expansion_Tcell_count5, by = "orig.ident")
# save(Expansion_Tcell_count, file = "TCR_group/Expansion_Tcell_count.RData")
# 
# #Define expanded clone
# Expansion_group2$Freq.2.group=''
# Expansion_group2$Freq.2.group[Expansion_group2$PairCTnt %in% Expansion_group3_Post$PairCTnt[Expansion_group3_Post$Freq.2.group=='yes']]="Expanded"
# Expansion_group2$Freq.2.group[Expansion_group2$Freq.2.group=='']="Non-expanded"
# unique(Expansion_group2$Freq.2.group)
# table(Expansion_group2$Freq.2.group)
# Expansion_group2$Freq.5.group=''
# Expansion_group2$Freq.5.group[Expansion_group2$PairCTnt %in% Expansion_group3_Post$PairCTnt[Expansion_group3_Post$Freq.5.group=='yes']]="Expanded"
# Expansion_group2$Freq.5.group[Expansion_group2$Freq.5.group=='']="Non-expanded"
# 
# Expansion_sample1=as.data.frame(prop.table(table(Expansion_group2$Freq.2.group, Expansion_group2$orig.ident)))
# Expansion_sample2=as.data.frame(prop.table(table(Expansion_group2$Freq.5.group, Expansion_group2$orig.ident)))
# save(Expansion_sample1,Expansion_sample2,file="TCR_group/Expansion_sample.RData")
# save(Expansion_group, Expansion_group2, Expansion_group3, Expansion_group3_Post,file="TCR_group/Expansion_group.RData")
# Expansion_sample1=Expansion_group3_Post %>% group_by(Freq.2.group,orig.ident) %>% count()
# Expansion_sample2=Expansion_group3_Post %>% group_by(Freq.5.group,orig.ident) %>% count()
# save(Expansion_sample1,Expansion_sample2,file="TCR_group/Expansion_sample.RData")

## Add TCR grouping label based on the following criteria
# 1) proportion of Expansion_sample1 Freq.2.group yes count in each sample's total immune cell count greater than 0.15;
# 2) proportion of Expansion_sample2 Freq.5.group yes count in each sample's total immune cell count greater than 0.10
load("TCR_group/Expansion_sample.RData")
Expansion_sample1 = Expansion_sample1[Expansion_sample1$Freq.2.group == "yes" & !is.na(Expansion_sample1$Freq.2.group),]
Expansion_sample2 = Expansion_sample2[Expansion_sample2$Freq.5.group == "yes" & !is.na(Expansion_sample2$Freq.5.group),]

load("TCR_group/Expansion_Tcell_count.RData")
load("TCR_group/Expansion_group.RData")

# scRNA = readRDS("../data/scRNAseq_data/RHTumor_clean.rds")
# unique(scRNA$Majortype)
# sc.immune = subset(scRNA, Majortype %in% c("T/NK cell", "Myeloid cell", "Plasma cell",
#                                            "B cell", "Mast cell"))
# rm(scRNA)
# gc()
# # Calculate total immune cell count of each sample
# immune.meta = sc.immune@meta.data
# immune.count = immune.meta %>% group_by(orig.ident) %>% count()
# immune.count = as.data.frame(immune.count)
# names(immune.count)[2] = "immune.cell.count"
# save(immune.count, file = "TCR_group/immune_cell_count.RData")
load("TCR_group/immune_cell_count.RData")
Expansion_sample1 = merge(Expansion_sample1, immune.count, by.x = "orig.ident")
Expansion_sample1 = Expansion_sample1[order(factor(Expansion_sample1$orig.ident, levels = unique(Expansion_sample1$orig.ident))),]
Expansion_sample2 = Expansion_sample2[order(factor(Expansion_sample2$orig.ident, levels = unique(Expansion_sample2$orig.ident))),]
all(Expansion_sample1$orig.ident == Expansion_sample2$orig.ident)
names(Expansion_sample1)[3] = "Freq.2.clone.number"
names(Expansion_sample2)[3] = "Freq.5.clone.number"
Expansion_samples = merge(Expansion_sample1, Expansion_sample2, by = "orig.ident")
Expansion_samples = merge(Expansion_samples, Expansion_Tcell_count, by = "orig.ident")
Expansion_samples$proportion.2 = Expansion_samples$expanded.2.clone.count / Expansion_samples$immune.cell.count
Expansion_samples$proportion.5 = Expansion_samples$expanded.5.clone.count / Expansion_samples$immune.cell.count

Expansion_samples_longer = Expansion_samples[,c("orig.ident","proportion.2","proportion.5")] %>% 
  pivot_longer(cols = c("proportion.2","proportion.5"),
               names_to = "group", values_to = "proportion")
# Plot a barplot to show the proportion of expanded clones in each sample
sample_order = c("RH08A","RH17A_T2","RH17A_T1","RH07A_T1","RH06A","RH04A_T2","RH14A",
                 "RH05A","RH10A","RH13A","RH04A_T1","RH07A_T2",
                 "RH16A","RH01A","RH02A_T2","RH03A","RH15A","RH11A","RH02A_T1","RH09A")
Expansion_samples_longer$orig.ident = factor(Expansion_samples_longer$orig.ident, levels = sample_order)
p = ggplot(Expansion_samples_longer, aes(x = orig.ident, y = proportion, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#FF9200","#FFD379")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Proportion of expanded clones", fill = "Group") +
  mytheme2 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p

######Fig. 2c######
immune = readRDS('元老师/RHTumor_immunocyte.rds')
immune$cellchat_cluster[immune$Sub_cluster=='CD8_Tex_CXCL13']="PD1+CD8"
immune$cellchat_cluster[immune$Sub_cluster=='CD4_Tht_CXCL13']="PD1+CD4"
immune$cellchat_cluster[immune$cellchat_cluster=='Treg']="CD4+T"

immune$group=paste0(immune$Expansion_normalized,'_',immune$progress)
immune$group2='Non-res'
immune$group2[immune$group=='HighE_Non_progress']='T_res'
immune$group2[immune$group=='LowE_Non_progress']="B_res"
immune$group3=paste0(immune$group2,immune$treatment)

merge_meta = immune@meta.data
cell_type = levels(factor(merge_meta$cellchat_cluster))
a <- merge_meta %>% group_by(cellchat_cluster,group3) %>%summarise(a=n())
ac <- merge_meta %>% group_by(group3) %>% summarise(ac=n())
ab <- merge_meta %>% group_by(cellchat_cluster) %>% summarise(ab=n())
abcd = dim(merge_meta)[1]

dist.data <- merge(a, ac)
dist.data <- merge(dist.data, ab)

dist.data['b'] = dist.data$ab - dist.data$a       
dist.data['c'] = dist.data$ac - dist.data$a
dist.data['d'] = abcd - dist.data$ab - dist.data$c

dist.data[,c('p.value','Ea', 'Ec', 'Eb', 'Ed')] = t(apply(dist.data, 1, function(x){
  x=as.numeric(x[c('a', 'c', 'b', 'd')])
  k.test = chisq.test(matrix(x, nrow=2,ncol=2))
  return(c(k.test$p.value, as.vector(k.test$expected)))
}))
dist.data['Reo'] = dist.data$a / dist.data$Ea

dist.heatmap.data <- dist.data[,c("cellchat_cluster", 'group3', 'Reo','p.value')]

dist.heatmap<-data.frame(matrix(data=0, nrow = length(cell_type), ncol=7))
rownames(dist.heatmap) = cell_type
colnames(dist.heatmap) = c('B_resPost','B_resPre', 'Non-resPost',  'Non-resPre', 'T_resPost', 'T_resPre' ,
                           'p.value')
for(i in 1:dim(dist.heatmap.data)[1]){
  ct = as.character(dist.heatmap.data[i, 'cellchat_cluster'])
  ts = dist.heatmap.data[i, 'group3']
  reo = dist.heatmap.data[i, 'Reo']
  p=dist.heatmap.data[i, 'p.value']
  dist.heatmap[ct, ts] = reo
  dist.heatmap[ct, 7] = p}

dist.heatmap1<-dist.heatmap[,1:6]
annote.heatmap = dist.heatmap
annote.heatmap$annote[annote.heatmap$p.value>=0.01] = ''
annote.heatmap$annote[annote.heatmap$p.value<0.01] = '*'
annote.heatmap$annote[annote.heatmap$p.value<0.001] = '**'
annote.heatmap$annote[annote.heatmap$p.value<0.0001] = '***'

annote.heatmap<-annote.heatmap[,7:8]

rownames(dist.heatmap1)<-factor(rownames(dist.heatmap1),level=cell_type)
dist.heatmap1 <- dist.heatmap1 %>% arrange(B_resPost) %>% arrange(B_resPre)
dist.heatmap1 <- dist.heatmap1[c("CD4+T","PD1+CD8","CD8+T","OtherT","NK",
                                 "Monocyte","Neutrophil","Macrophage","ClassicalDC",
                                 "pDC","Mast cell","B cell","Plasma cell"),
                               c("B_resPre","B_resPost","T_resPre","T_resPost",
                                 "Non-resPre","Non-resPost")]
cell_cluster_dist=pheatmap(dist.heatmap1,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           treeheight_row = 0, 
                           color = colorRampPalette(c("#5CACEE","white","#F08080"))(50),
                           #display_numbers =as.matrix(annote.heatmap),
                           display_numbers =round(as.matrix(dist.heatmap1), 2),
                           cellwidth = 30,cellheight = 16, fontsize = 10,
                           border_color = '#ffffff',
                           angle_col='45',
                           main = 'cell type',
                           breaks = seq(0,2,by=2/50)
)
rm(immune, immunocyte)
gc()

######Fig. 2d######
# library(monocle)
# bcell = readRDS("元老师/bcell-all.rds")
# bcell@meta.data$new_group[bcell@meta.data$progress =='Progress'] <- 'Non-responder'
# bcell@meta.data$new_group[bcell@meta.data$progress =='Non_progress' & bcell@meta.data$Expansion_normalized == 'HighE'] <- 'T-responder'
# bcell@meta.data$new_group[bcell@meta.data$progress =='Non_progress' & bcell@meta.data$Expansion_normalized == 'LowE'] <- 'B-responder'
# bcell@meta.data$new_group2 <- paste0(bcell@meta.data$treatment,'_',bcell@meta.data$new_group)
# 
# bcell_matirx<-as(as.matrix(bcell@assays$RNA@counts), 'sparseMatrix')
# sample_ann<-bcell@meta.data
# identical(rownames(sample_ann),colnames(bcell_matirx))
# 
# bcell_pd<-new("AnnotatedDataFrame", data =sample_ann)
# bcell_feature<-data.frame(gene_id=rownames(bcell_matirx),gene_short_name=rownames(bcell_matirx))
# rownames(bcell_feature)<-rownames(bcell_matirx)
# bcell_fd<-new("AnnotatedDataFrame", data = bcell_feature)
# cds<-newCellDataSet(bcell_matirx,
#                     phenoData =bcell_pd,featureData =bcell_fd,
#                     lowerDetectionLimit = 0.5,
#                     expressionFamily = negbinomial.size())
# 
# cds <- estimateSizeFactors(cds)
# cds <- estimateDispersions(cds)
# 
# cds <- detectGenes(cds, min_expr = 0.1)
# expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
# 
# bcell_cds <- cds
# bcell_cds <- bcell_cds[expressed_genes,]
# Idents(bcell)<-"Sub_cluster"
# deg.cluster <- FindAllMarkers(bcell,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
# diff.genes <- subset(deg.cluster,p_val_adj<0.05)$gene
# 
# clustering_DEG <- differentialGeneTest(bcell_cds, fullModelFormulaStr = '~Sub_cluster')
# clustering_DEG_genes <- subset(clustering_DEG, qval < 0.01) 
# clustering_DEG_genes <- clustering_DEG_genes[order(clustering_DEG_genes$qval,decreasing=F),]
# 
# diff.genes <- intersect(diff.genes,row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:130])
# bcell_cds <- setOrderingFilter(bcell_cds, diff.genes)
# 
# # bcell_ordering_genes <- rownames(clustering_DEG_genes)
# # bcell_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:500]
# # 
# # bcell_cds <- setOrderingFilter(bcell_cds, ordering_genes = bcell_ordering_genes)
# 
# # # 或者 选择高度离散基因
# # disp_table <- dispersionTable(bcell_cds)
# # bcell_ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 *dispersion_fit)$gene_id)
# # bcell_cds <- setOrderingFilter(bcell_cds, bcell_ordering_genes)
# 
# bcell_cds <- reduceDimension(bcell_cds, max_components = 2, method = 'DDRTree')
# bcell_cds <- orderCells(bcell_cds)
# 
# library(ggridges)
# seurat<-pData(bcell_cds)
# seurat<-seurat[,c("Pseudotime","State")]
# bcell<-AddMetaData(bcell,seurat)

bcell_cds = readRDS("元老师/Bcell monocle241229.rds")
plotdf=pData(bcell_cds)
ggplot(plotdf, aes(x=Pseudotime,y=group,fill=group))+
  geom_density_ridges(scale=0.8) +
  #geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_classic()+
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour = 'black',size = 10),
    axis.text.y = element_text(colour = 'black',size = 10)
  )+
  scale_fill_igv()+
  NoLegend()+
  scale_x_continuous(position ="top")

## Set time group
plotdf$pseudotime_group = ifelse(plotdf$Pseudotime <= 5, "Early", 
                                 ifelse(plotdf$Pseudotime > 10, "Late", "Middle"))
time.sample.freq = table(plotdf$orig.ident, plotdf$pseudotime_group) %>% as.data.frame()
time.sample.freq = time.sample.freq %>%
  pivot_wider(names_from = Var2,
              values_from = Freq) %>%
  column_to_rownames("Var1")
time.sample.freq = time.sample.freq[,c("Early", "Middle", "Late")]
if (!dir.exists("Visualization/Figure2/Fig_2d")) {
  dir.create("Visualization/Figure2/Fig_2d", recursive = T)
}
write.xlsx(time.sample.freq, "Visualization/Figure2/Fig_2d/Bcell_monocle_Time_Group.xlsx",
           rowNames = T)

######Fig. 2e######
b_bcr = readRDS("元老师/补-单细胞BCR/contig out old/b_bcr-t p pbmc.rds")
b_bcr = subset(b_bcr,tissue=="Tumor")

plasma=b_bcr@meta.data[b_bcr@meta.data$Sub_cluster%in%c("pB_IGHG1",'pB_STMN1'),]
plasma$group[plasma$expansion_group=='LowE_Non_progress']='B_responder'
plasma$group[plasma$expansion_group=='HighE_Non_progress']='T_responder'
plasma$group[plasma$expansion_group%in% c('HighE_Progress','LowE_Progress')]='Non_responder'
plasma$group2=paste0(plasma$treatment,plasma$group)
plasma$group2=factor(plasma$group2,c('PreB_responder','PostB_responder','PreT_responder','PostT_responder','PreNon_responder','PostNon_responder'))
cell.prop <-as.data.frame(prop.table(table(plasma$c_gene,plasma$group2), margin = 2))
colnames(cell.prop)<-c("sample","celltype","proportion")
ggplot(cell.prop,aes(x=celltype,proportion,fill=sample))+
  geom_bar(stat = "identity",position="fill")+ggtitle("")+theme_classic()+ 
  theme(axis.ticks.length=unit(0.5,'cm'))+xlab(" ")+
  guides(fill=guide_legend(title=NULL))+RotatedAxis() + 
  scale_fill_manual(values = color.s1)

######Fig. 2g######
cell.prop <-as.data.frame(prop.table(table(plasma$cloneType,plasma$group2), margin = 2))
colnames(cell.prop)<-c("sample","celltype","proportion")
p = ggplot(cell.prop,aes(x=celltype,proportion,fill=sample))+
  geom_bar(stat = "identity",position="fill")+ggtitle("")+theme_classic()+
  theme(axis.ticks.length=unit(0.5,'cm'))+xlab(" ")+
  guides(fill=guide_legend(title=NULL))+RotatedAxis() + 
  scale_fill_manual(values = my.colors)
p

######Fig. 2h######
all = readRDS('元老师/RHTumor_immunocyte.rds')
Idents(all)="Majortype"
all=subset(all,idents=c("B cell","Myeloid cell","Plasma cell","T/NK cell"))
all_post=subset(all,treatment=="Post")
all_post=subset(all_post,Expansion_normalized=="LowE")
all_post=subset(all_post,progress=="Non_progress")

immune.prop <-as.data.frame(prop.table(table(all_post$Sub_cluster, all_post$orig.ident), margin = 2))
colnames(immune.prop)<-c("celltype","sample","proportion")
immune.prop2 <- reshape(immune.prop,direction = "wide",timevar = "celltype",idvar = "sample")

colnames(immune.prop2) <- c("sample",levels(immune.prop$celltype))
cell.prop <- as.data.frame(lapply(immune.prop2, as.numeric))
cell.prop=cell.prop[,-1]
cor_data <- cor(cell.prop,method='pearson')

temp_t<-data.frame(table(all_post$orig.ident,all_post$Sub_cluster))
colnames(temp_t) <- c("Samples",'Cell_type',"Percentage of Cells")
temp_t2<-data.frame(table(all_post$orig.ident))
colnames(temp_t2) <- c("Samples",'total')

temp_tall=merge.data.frame(temp_t,temp_t2,all.x = T,by = 'Samples')
temp_tall$`Percentage of Cells`=temp_tall$`Percentage of Cells`/temp_tall$total

temp_pc=temp_tall
temp_pc=temp_pc[temp_pc$Cell_type=='pB_IGHG1',]

temp_all_T=temp_tall
needed=unique(temp_t$Cell_type)[4:11]

temp_all_T=temp_all_T[temp_all_T$Cell_type %in% needed,]
tcell=readRDS("元老师/Tcell_TCR.rds")
tcell_needed=subset(tcell,Expansion_group=="LE_NP")
tcell_needed2=subset(tcell_needed,treatment=="Post")

needed=unique(temp_tall$Samples)
temp_all=rbind(temp_pc,temp_all_T)
temp_all=temp_all[temp_all$Samples%in%tcell_needed2$orig.ident,]
temp_all=temp_all[temp_all$Samples%in%needed,]

my_cell=as.vector(unique(temp_all$Cell_type))
j=1
my_cell[j]
my_cell_tmp=my_cell[-j]
my_rev=c('PCC','cell')
for (i in 1:length(my_cell_tmp)) {
  # i=1
  my_target1=my_cell[j]
  my_target2=my_cell_tmp[i]
  temp_all1=temp_all[temp_all$Cell_type==my_target1,]
  temp_all2=temp_all[temp_all$Cell_type==my_target2,]
  temp_all1=temp_all1$`Percentage of Cells`
  temp_all2=temp_all2$`Percentage of Cells`
  my_cor=cor.test(temp_all1,temp_all2)
  my_cor=my_cor$estimate
  my_cor=c(my_cor,my_target2)
  my_cor
  my_rev=rbind(my_rev,my_cor)
}
colnames(my_rev)=my_rev[1,]
my_rev=my_rev[-1,]
my_rev=as.data.frame(my_rev)
my_rev$PCC=as.numeric(my_rev$PCC)
ggdotchart(my_rev, x = "cell", y = "PCC",
           color = "cell",                            
           #palette = colorvec,
           sorting = "ascending",                        
           add = "segments",                           
           ggtheme = theme_pubr(),                     
           xlab="",
           ylab = "Correlation with Plasma_cell",
           add.params = list(color = "lightgray", size = 2),
           dot.size = 9,                                 
           label = round(my_rev$PCC,digits = 2),                      
           font.label = list(color = "black", size = 10, face="bold",
                             vjust = 0.5), 
)+ theme(axis.text.x = element_text(angle = 45, hjust = 0.9,vjust = 0.9),legend.position = 'none')+
  scale_color_manual(values = my.colors1)


######Fig. 2i######
st.rds = dir("../data/slice7.results.renew/reanno_cjy_T_B_mye_endo_fib/") %>% .[grep("^RH", .)]
patients = sapply(st.rds, function(x){
  patient = unlist(strsplit(x, "_"))[1]
})
patients = unname(patients)
fin.group = read.xlsx("sample.cor.group_3group_new.xlsx")

#####MKI67+Tfh inside TLS status#####
tls.selected = c("STRH05B_10","STRH05B_11",
                 "STRH05B_12","STRH05B_2",
                 "STRH05B_21","STRH05B_25",
                 "STRH05B_30","STRH05B_4",
                 "STRH05B_6","STRH05B_9")
tfh.density.all = data.frame()
for (patient in patients[-grep("P$",patients)]) {
  cat("Processing patient", patient, "...\n")
  rdsfile = dir("../data/slice7.results.renew/reanno_cjy_T_B_mye_endo_fib/") %>% .[grep(patient, .)]
  stRNA = readRDS(paste0("../data/slice7.results.renew/reanno_cjy_T_B_mye_endo_fib/", rdsfile))
  meta.data = stRNA@meta.data
  
  if (!"TLS_location" %in% names(meta.data)) {
    meta.data$TLS_location = "M-T_TLS"
  }
  if (!"Bin_Region" %in% names(meta.data)) {
    meta.data$Bin_Region = "Tumor_side_of_Margin_area"
  }

  tfh.bins = meta.data$bins[meta.data$celltype %in% c("Tfh_BCL6","Tpex_TCF7","Tm_GZMK")]
  mki67.bins = stRNA@assays$Spatial@counts["MKI67",] %>% .[. > 0]
  tfh.bins = intersect(tfh.bins, names(mki67.bins))

  meta.data$cell_id = rownames(meta.data)
  if (!"TLS_maturity" %in% names(meta.data) & "TLS_Final" %in% names(meta.data)) {
    meta.data$interaction = ifelse(grepl("_NA$", meta.data$TLS_Final), "Others", "TLS")
  } else if (!"TLS_Final" %in% names(meta.data)) {
    cat("No TLS info in", patient, "\n")
    next
  } else if (patient == "RH05B") {
    meta.data$interaction = ifelse(meta.data$TLS_raw %in% tls.selected, "TLS", "Others")
  } else {
    meta.data$interaction = ifelse(meta.data$TLS_maturity %in% c("Mature","Conforming","Deviating"), "TLS", "Others")
  }
  meta.data$interaction[rownames(meta.data) %in% tfh.bins] = "Tfh_MKI67"
  stRNA@meta.data = meta.data
  prefix = paste0("st_T_B_interaction_HZY/TFH_MKI67/", patient)
  plot_item = "interaction"
  if (plot_item %in% colnames(stRNA@meta.data)) {
    dimplot_inhouse_R(obj = stRNA, plot_item = plot_item, prefix = prefix,
                      tmp_color = c(ggsci::pal_nejm()(1), "#87CEEBB3"))
  } else {
    tryCatch({
      print(paste0(patient, " does not have ", plot_item, " information in the meta.data"))
    }, error = function(e) {
      print(paste0("Error occurred: ", e))
    })
  }
  
  # Calculate Tfh density inside TLS and outside TLS
  if (!"TLS_maturity" %in% names(meta.data)) {
    tfh.inside = nrow(meta.data[meta.data$interaction == "Tfh_MKI67" & !grepl("_NA$", meta.data$TLS_Final),])
    tfh.outside = nrow(meta.data[meta.data$interaction == "Tfh_MKI67" & grepl("_NA$", meta.data$TLS_Final),])
    tls.count = nrow(meta.data[!grepl("_NA$", meta.data$TLS_Final),])
    outside.count = nrow(meta.data[grep("_NA$", meta.data$TLS_Final),])
  } else if (patient == "RH05B") {
    tfh.inside = nrow(meta.data[meta.data$interaction == "Tfh_MKI67" & meta.data$TLS_raw %in% tls.selected,])
    tfh.outside = nrow(meta.data[meta.data$interaction == "Tfh_MKI67" & !meta.data$TLS_raw %in% tls.selected,])
    tls.count = nrow(meta.data[meta.data$TLS_raw %in% tls.selected,])
    outside.count = nrow(meta.data[!meta.data$TLS_raw %in% tls.selected,])
  } else {
    tfh.inside = nrow(meta.data[meta.data$interaction == "Tfh_MKI67" & meta.data$TLS_maturity %in% c("Mature","Conforming","Deviating"),])
    tfh.outside = nrow(meta.data[meta.data$interaction == "Tfh_MKI67" & grepl("_NA$", meta.data$TLS_Final),])
    tls.count = nrow(meta.data[meta.data$TLS_maturity %in% c("Mature","Conforming","Deviating"),])
    outside.count = nrow(meta.data[grep("_NA$", meta.data$TLS_Final),])
  }
  tfh.inside.density = tfh.inside / tls.count
  tfh.outside.density = tfh.outside / outside.count
  tfh.density.df = data.frame(orig.ident = patient, tfh.inside.density = tfh.inside.density, tfh.outside.density = tfh.outside.density)
  tfh.density.all = rbind(tfh.density.all, tfh.density.df)
}
if (!dir.exists("st_T_B_interaction_HZY/TFH_MKI67/")) {
  dir.create("st_T_B_interaction_HZY/TFH_MKI67/", recursive = T)
}
write.xlsx(tfh.density.all, "st_T_B_interaction_HZY/TFH_MKI67/tfh.density.all.xlsx", rowNames = F)

stRNA = readRDS("../data/slice7.results.renew/reanno_cjy_T_B_mye_endo_fib/RH08B_rhcc_reanno_T_B_mye_endo_fib.rds")
meta.data = stRNA@meta.data
meta.data$interaction = ifelse(meta.data$TLS_maturity %in% c("Mature","Conforming","Deviating"), "TLS", "Others")
tfh.bins = meta.data$bins[meta.data$celltype %in% c("Tfh_BCL6","Tpex_TCF7","Tm_GZMK")]
mki67.bins = stRNA@assays$Spatial@counts["MKI67",] %>% .[. > 0]
tfh.bins = intersect(tfh.bins, names(mki67.bins))
meta.data$interaction[rownames(meta.data) %in% tfh.bins] = "Tfh_MKI67"
meta.data$cell_id = rownames(meta.data)
meta.data = meta.data[,c("cell_id", "interaction")]

stRNA = readRDS("st_T_B_interaction_HZY/RH08B_Bin20_tissue_only.rds")
orig.meta = stRNA@meta.data
orig.meta$cell_id = rownames(orig.meta)
orig.meta = merge(orig.meta, meta.data, by = "cell_id", all.x = T)
orig.meta$interaction[is.na(orig.meta$interaction)] = "Others"
rownames(orig.meta) = orig.meta$cell_id
orig.meta = orig.meta[order(factor(orig.meta$cell_id, levels = colnames(stRNA))),]
all(rownames(orig.meta) == colnames(stRNA))
stRNA@meta.data = orig.meta
prefix = "st_T_B_interaction_HZY/TFH_MKI67/RH08B"
plot_item = "interaction"
dimplot_inhouse_R(obj = stRNA, plot_item = plot_item, prefix = prefix,
                  tmp_color = c(ggsci::pal_nejm()(1), "#87CEEBB3"))

######Fig. 2j######
tfh.density.all = read.xlsx("st_T_B_interaction_HZY/TFH_MKI67/tfh.density.all.xlsx", detectDates = F)
tfh.density.all = merge(tfh.density.all, fin.group, by = "orig.ident")
## make column "tfh.inside.density" and "tfh.outside.density" into long format
tfh.density.all = tfh.density.all %>%
  pivot_longer(cols = c("tfh.inside.density", "tfh.outside.density"),
               names_to = "density_type", values_to = "density")
tfh.density.all$density_type = ifelse(tfh.density.all$density_type == "tfh.inside.density", "Inside_TLS", "Outside_TLS")

tfh.density.bres = tfh.density.all[tfh.density.all$group == "B_cell_responder",]
stat.test = tfh.density.bres %>%
  t_test(density ~ density_type) %>%
  add_xy_position(x = "density_type")
ggplot(tfh.density.bres, aes(x = density_type, 
                            y = density, group = density_type)) +
  geom_boxplot(aes(color = density_type), width = .6, outlier.shape = NA) +
  geom_jitter(aes(color = density_type), size = 1, position = position_jitterdodge(.3)) +
  scale_color_manual(values = rev(two.colors)) +
  stat_pvalue_manual(stat.test, label = "p", bracket.size = p.lwd, label.size = 5) +
  labs(x = "",
       y = "Density of MKI67+Tfh") +
  mytheme2 + theme(legend.position = "none")
myggsave(plot = last_plot(), filename = "st_T_B_interaction_HZY/TFH_MKI67/tfh.density.b_res.pdf", width = 4, height = 6)

tfh.density.tres.nonres = tfh.density.all[tfh.density.all$group %in% c("T_cell_responder","Non-responder"),]
stat.test = tfh.density.tres.nonres %>%
  t_test(density ~ density_type) %>%
  add_xy_position(x = "density_type")
ggplot(tfh.density.tres.nonres, aes(x = density_type, 
                            y = density, group = density_type)) +
  geom_boxplot(aes(color = density_type), width = .6, outlier.shape = NA) +
  geom_jitter(aes(color = density_type), size = 1, position = position_jitterdodge(.3)) +
  scale_color_manual(values = rev(two.colors)) +
  stat_pvalue_manual(stat.test, label = "p", bracket.size = p.lwd, label.size = 5) +
  labs(x = "",
       y = "Density of MKI67+Tfh") +
  mytheme2 + theme(legend.position = "none")
myggsave(plot = last_plot(), filename = "st_T_B_interaction_HZY/TFH_MKI67/tfh.density.t_res.non_res.pdf", width = 4, height = 6)


#####Figure 3#####
######Fig. 3B######
raw.dir = "../data/bcr_data/new_data/"
raw.bcr = dir(raw.dir) %>% .[grep("^RH05",.)]
## load the bcr data that has been merged with single cell bcell set
load("scRepertoire/BCR/combined_83samples.RData")
new.bc.set = readRDS("scRepertoire/rhcc_bcr_mergebcr.rds")
meta = new.bc.set@meta.data
unique.rh = unique(new.bc.set$orig.ident) %>% 
  str_extract(., "RH05") %>% 
  unique()

pass.list = readRDS("BCR_sharing_cjy_chordDiagram/shm_file/rhcc_bcr_shm_all_IMGT_VDJ.rds")
load("../data/bcr_data/BCR_phenotype.RData")

all.clonotypes.expand.prop.df = data.frame()
clonetypes.quantity = data.frame()
rh = "RH05"
filter.contigs = raw.bcr[grep(rh, raw.bcr)]
## Check if the single RH sample contains pre , post and T, PBMC samples
filter.samples = gsub("_filtered_contig_annotations.csv","",filter.contigs)
groups.df = groups[groups$orig.ident %in% filter.samples,]
groups.df = table(groups.df$tissue, groups.df$treatment) %>% as.data.frame() %>%
  .[.$Var1 != "liver",]
if (all(groups.df$Freq != 0)) {
  cat(rh, "contains all the samples\n")
} else if (all(groups.df$Freq[groups.df$Var1 == "tumor"] != 0)) {
  cat(rh, "contains pre and post tumor samples\n")
} else {
  cat(rh, "does not contain", c(as.character(groups.df$Var1[groups.df$Freq == 0]),
                                as.character(groups.df$Var2[groups.df$Freq == 0])), "\n")
  next
}
## remove the P type contigs
filter.contigs = filter.contigs[!grepl("*P$|*P[1-9]$",filter.samples)]
## remove the PBMC type contigs
filter.contigs = filter.contigs[!grepl("PBMC",filter.samples)]
filter.contigs = filter.contigs[!is.na(filter.contigs)]
rh.contigs = data.frame()
rh.contigs.all = data.frame()
for (contig in filter.contigs){
  contig.df = read.csv(paste0(raw.dir, contig), header = T, check.names = F)
  sample = gsub("_filtered_contig_annotations.csv","",contig)
  contig.df$orig.ident = sample
  should.column = c("barcode","is_cell","contig_id","high_confidence","length","chain",
                    "v_gene", "d_gene", "j_gene", "c_gene", "full_length", "productive",
                    "fwr1","fwr1_nt","cdr1","cdr1_nt","fwr2","fwr2_nt","cdr2","cdr2_nt","fwr3",
                    "fwr3_nt","cdr3","cdr3_nt", "fwr4", "fwr4_nt","reads", "umis", "raw_clonotype_id", "raw_consensus_id",
                    "orig.ident")
  if (FALSE %in% (should.column %in% names(contig.df))){
    not.contain.columns = should.column[!should.column %in% names(contig.df)]
    contig.df[,not.contain.columns] = ""
    contig.df = contig.df[,should.column]
  } else {
    contig.df = contig.df[,should.column]
  }
  ## Add columns recording full-length cdr
  contig.df$cdrs = paste0(contig.df$cdr1, contig.df$cdr2, contig.df$cdr3)
  contig.df$cdrs_nt = paste0(contig.df$cdr1_nt, contig.df$cdr2_nt,
                             contig.df$cdr3_nt)
  ## Use all chain:
  rh.contigs.all = rbind(rh.contigs.all, contig.df)
  ## Use only Heavy chain:
  contig.df = contig.df[contig.df$chain == "IGH",]
  rh.contigs = rbind(rh.contigs, contig.df)
}
rh.msa.contigs = cjy_Identical_cdr3nt_from_filtered_contig(merged_contigs = rh.contigs)
expanded.not.de.novo.type = c()
expanded.not.de.novo.type.df = data.frame()
de.novo.type = c()
expanded.not.de.novo.type.prob = c()
for (type in unique(rh.msa.contigs$new_Clonotype_ID)){
  single.type.df = rh.msa.contigs[rh.msa.contigs$new_Clonotype_ID == type, ]
  pre.df = single.type.df[-grep("A", single.type.df$orig.ident),]
  post.df = single.type.df[grep("A", single.type.df$orig.ident),]
  pre.samples = unique(pre.df$orig.ident)
  post.samples = unique(post.df$orig.ident)
  for (post.s in post.samples) {
    post.df.single = post.df[post.df$orig.ident == post.s,]
    if (nrow(pre.df) == 0 & nrow(post.df.single) > 0){
      diff.value = nrow(post.df.single) - nrow(pre.df)
      names(diff.value) = type
      de.novo.type = c(de.novo.type, diff.value)
    } else if (nrow(pre.df) > 0 & (nrow(post.df.single) - nrow(pre.df) > 0)) {
      diff.value = nrow(post.df.single) - nrow(pre.df)
      diff.prob = diff.value/nrow(pre.df)
      names(diff.value) = type
      names(diff.prob) = type
      expanded.not.de.novo.type = c(expanded.not.de.novo.type, diff.value)
      expanded.not.de.novo.type.prob = c(expanded.not.de.novo.type.prob, diff.prob)
      expanded.not.de.novo.type.tmp = data.frame(from = paste0(pre.samples, "_", type),
                                                 to = paste0(post.samples, "_", type),
                                                 value1 = nrow(pre.df),
                                                 value2 = nrow(post.df.single),
                                                 new_Clonotype_ID = type)
      expanded.not.de.novo.type.df = rbind(expanded.not.de.novo.type.df,
                                           expanded.not.de.novo.type.tmp)
    }
  }
}
expanded.not.de.novo.type = sort(expanded.not.de.novo.type, decreasing = T)
expanded.not.de.novo.type.prob = sort(expanded.not.de.novo.type.prob, decreasing = T)
de.novo.type = sort(de.novo.type, decreasing = T)
## Export increased clone type proportion
# increased.proportion = length(c(expanded.not.de.novo.type)) / length(unique(rh.msa.contigs$new_Clonotype_ID))
exp.num = ifelse(is.null(expanded.not.de.novo.type), 0, length(expanded.not.de.novo.type))
exp.prop = as.numeric(exp.num) / as.numeric(length(unique(rh.msa.contigs$new_Clonotype_ID)))
de.novo.num = ifelse(is.null(de.novo.type[de.novo.type>=5]), 0, length(de.novo.type[de.novo.type>=5]))
de.novo.prop = as.numeric(de.novo.num) / as.numeric(length(unique(rh.msa.contigs$new_Clonotype_ID)))
not.increased.num = length(unique(rh.msa.contigs$new_Clonotype_ID)) - de.novo.num - exp.num
not.increased.prop = as.numeric(not.increased.num) / as.numeric(length(unique(rh.msa.contigs$new_Clonotype_ID)))
increased.df = data.frame(patient = rh,
                          Freq = c(as.numeric(exp.num), as.numeric(de.novo.num), 
                                   as.numeric(not.increased.num)),
                          type = c("Expanded","De_novo","Not_increased"),
                          proportion = c(exp.prop, de.novo.prop, not.increased.prop))
clonetypes.quantity = rbind(clonetypes.quantity, increased.df)
## Get the expand and not expand new_Clonotype_ID based on Post tumor
rh.a.tumor = rh.msa.contigs %>%
  .[grep("A", .$orig.ident),] #%>%
# .[-grep("PBMC|P$", .$orig.ident),]
bcr.freq = table(rh.a.tumor$new_Clonotype_ID) %>% as.data.frame()
exp.type = bcr.freq$Var1[bcr.freq$Freq > 10]
rh.expand = rh.msa.contigs[rh.msa.contigs$new_Clonotype_ID %in% exp.type,]
rh.not.expand = rh.msa.contigs[rh.msa.contigs$new_Clonotype_ID %nin% exp.type,]
if (length(expanded.not.de.novo.type.prob) >= 10){
  expanded.not.de.novo.type.prob.top = expanded.not.de.novo.type.prob[1:10]
  expanded.not.de.novo.type.top.df = expanded.not.de.novo.type.df[expanded.not.de.novo.type.df$new_Clonotype_ID %in% names(expanded.not.de.novo.type.prob.top),]
  rh.top.expand = rh.msa.contigs[rh.msa.contigs$new_Clonotype_ID %in% names(expanded.not.de.novo.type.prob)[1:10],]
} else {
  expanded.not.de.novo.type.prob.top = expanded.not.de.novo.type.prob[1:length(expanded.not.de.novo.type.prob)]
  expanded.not.de.novo.type.top.df = expanded.not.de.novo.type.df[expanded.not.de.novo.type.df$new_Clonotype_ID %in% names(expanded.not.de.novo.type.prob.top),]
  rh.top.expand = rh.msa.contigs[rh.msa.contigs$new_Clonotype_ID %in% names(expanded.not.de.novo.type.prob)[1:length(expanded.not.de.novo.type.prob)],]
}
if (length(de.novo.type) >= 10){
  rh.de.novo = rh.msa.contigs[rh.msa.contigs$new_Clonotype_ID %in% names(de.novo.type)[1:10],]
} else {
  rh.de.novo = rh.msa.contigs[rh.msa.contigs$new_Clonotype_ID %in% names(de.novo.type)[1:length(de.novo.type)],]
}
cjy_getCirclized_after_cjy_MSA_from_data_frame(expanded.not.de.novo.type.top.df,
                                               plot_dir = "Visualization/Figure3/Fig_3B/",
                                               sample = paste0(rh,"_top10_expanded"),
                                               chord_transparency = .3,
                                               link.not.display.threshold = 0,
                                               niceFacing = T, big.gap = 2,small.gap = 2)

#####Figure 4#####
######Fig. 4B######
st.bcr.dir = "st_ter_BCR_new/"
bc.set = readRDS("scRepertoire/rhcc_bcr_mergebcr.rds")
bc.meta = bc.set@meta.data
le.samples = c("RH10B", "RH10T", "RH03T", "RH05B", "RH06B", "RH14B","RH15B","RH16B","RH08B","RH17B")
le.group = data.frame(sample = le.samples,
                      group = c("B_cell_responder","B_cell_responder","Non_responder","B_cell_responder",
                                "B_cell_responder","B_cell_responder","Non-responder",
                                "Non-responder","T_cell_responder","T_cell_responder"))
hbv.core.id = readRDS("st_ter_BCR_new/talent_version2/talent_result.rds")
load("st_ter_BCR_new/rh_st_ter_BCR_heavy_chain.RData")
load("st_ter_BCR_new/meta_st_B.RData")
shm.df = bcr.df[,c("sequence_id","sequence","shm","clone_id")]
hbv.core.id = merge(hbv.core.id, shm.df, by = "sequence_id")
bcr.df$hbv_group = ifelse(bcr.df$sequence_id %in% hbv.core.id$sequence_id, 
                          "HBV_core", "Other_BCR")
bcr = bcr.df[,c("BinID","hbv_group")]
meta = merge(meta.st, bcr, by = "BinID")
cell.df = table(meta$Cellsubtype, meta$hbv_group) %>% as.data.frame()
cell.prop = cjy_single_cell_proportion_test(single_cell_freq_input = cell.df,
                                            column_one = "celltype",
                                            column_two = "clone")
p = ggplot(cell.df, aes(Var2, Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "", y = "Proportion", title = "") + 
  scale_fill_manual(values = rev(color.s1[1:4])) + 
  guides(fill = guide_legend(title = "")) + 
  mytheme2
p
myggsave("Visualization/Figure4/Fig_4B/HBV_core_celltype.pdf", width = 6, height = 6, plot = p)

######Fig. 4C & D######
raw.dir = "../data/bcr_data/new_data/"
raw.bcr = dir(raw.dir) %>% .[grep("^RH",.)]
hbv.ab = read.xlsx("antibody_similarity/HBcAb阳性抗体 20240910.xlsx")
hbv.ab$`全长重链核苷酸` = gsub("\n", "", hbv.ab$`全长重链核苷酸`)
hbv.ab$`全长轻链核苷酸` = gsub("\n", "", hbv.ab$`全长轻链核苷酸`)
hbv.ab$`全长重链氨基酸` = gsub("\n", "", hbv.ab$`全长重链氨基酸`)
hbv.ab$`全长轻链氨基酸` = gsub("\n", "", hbv.ab$`全长轻链氨基酸`)
hbv.ab$`全长重链核苷酸` = gsub(" ", "", hbv.ab$`全长重链核苷酸`)
hbv.ab$`全长轻链核苷酸` = gsub(" ", "", hbv.ab$`全长轻链核苷酸`)
hbv.ab$`全长重链氨基酸` = gsub(" ", "", hbv.ab$`全长重链氨基酸`)
hbv.ab$`全长轻链氨基酸` = gsub(" ", "", hbv.ab$`全长轻链氨基酸`)

non.hbv.ab = read.xlsx("antibody_similarity/全部合成抗体.xlsx")
non.hbv.ab = non.hbv.ab[is.na(non.hbv.ab$BLI.HBcAg),]

new.bc.set = readRDS("scRepertoire/rhcc_bcr_mergebcr.rds")
unique.rh = unique(new.bc.set$orig.ident) %>% 
  str_extract(., "RH[0-9][0-9]") %>% 
  unique()
all.contigs.to.check = data.frame()
bcr.contigs = data.frame()
for (rh in unique.rh){
  filter.contigs = raw.bcr[grep(rh, raw.bcr)]
  filter.samples = gsub("_filtered_contig_annotations.csv","",filter.contigs)
  rh.contigs = data.frame()
  rh.contigs.all = data.frame()
  for (contig in filter.contigs){
    contig.df = read.csv(paste0(raw.dir, contig), header = T, check.names = F)
    sample = gsub("_filtered_contig_annotations.csv","",contig)
    contig.df$orig.ident = sample
    should.column = c("barcode","is_cell","contig_id","high_confidence","length","chain",
                      "v_gene", "d_gene", "j_gene", "c_gene", "full_length", "productive",
                      "fwr1","fwr1_nt","cdr1","cdr1_nt","fwr2","fwr2_nt","cdr2","cdr2_nt","fwr3",
                      "fwr3_nt","cdr3","cdr3_nt", "fwr4", "fwr4_nt","reads", "umis", "raw_clonotype_id", "raw_consensus_id",
                      "orig.ident")
    if (FALSE %in% (should.column %in% names(contig.df))){
      not.contain.columns = should.column[!should.column %in% names(contig.df)]
      contig.df[,not.contain.columns] = ""
      contig.df = contig.df[,should.column]
    } else {
      contig.df = contig.df[,should.column]
    }
    ## Add columns recording full-length cdr
    contig.df$cdrs = paste0(contig.df$fwr1, contig.df$cdr1, 
                            contig.df$fwr2, contig.df$cdr2, 
                            contig.df$fwr3, contig.df$cdr3,
                            contig.df$fwr4)
    contig.df$cdrs_nt = paste0(contig.df$fwr1_nt, contig.df$cdr1_nt, 
                               contig.df$fwr2_nt, contig.df$cdr2_nt,
                               contig.df$fwr3_nt, contig.df$cdr3_nt,
                               contig.df$fwr4_nt)
    ## Use all chain:
    rh.contigs.all = rbind(rh.contigs.all, contig.df)
    ## Use only Heavy chain:
    contig.df = contig.df[contig.df$chain == "IGH",]
    rh.contigs = rbind(rh.contigs, contig.df)
  }
  all.contigs.to.check = rbind(all.contigs.to.check, rh.contigs.all)
  rh.msa.contigs = cjy_Identical_cdr3nt_from_filtered_contig(merged_contigs = rh.contigs)
  bcr.contigs = rbind(bcr.contigs, rh.msa.contigs)
}
bcr.contigs$antibody_type = "Non-APR"
for (i in 1:nrow(bcr.contigs)){
  if (bcr.contigs$cdrs_nt[i] %in% hbv.ab$`全长重链核苷酸`){
    bcr.contigs$antibody_type[i] = "HBc"
  } else if (bcr.contigs$cdrs_nt[i] %in% non.hbv.ab$`全长重链核苷酸`){
    bcr.contigs$antibody_type[i] = "Non-HBc"
  }
}
bcr.contigs$patient <- substr(bcr.contigs$orig.ident, 1, 4)

tissue_names <- c("Tumor", "Peri_tumor", "PBMC")
combinations <- data.frame(
  Tumor = c(1, 1, 1, 0, 1, 0, 0),
  Peri_tumor = c(1, 0, 1, 1, 0, 1, 0),
  PBMC = c(1, 1, 0, 1, 0, 0, 1)
)
combinations$id <- 1:nrow(combinations)
long_data <- reshape2::melt(combinations, id.vars = "id")
long_data <- long_data[long_data$value == 1,]
long_data$y <- as.numeric(factor(long_data$variable, levels = rev(tissue_names)))
long_data$x <- long_data$id
p1 = ggplot() + 
  geom_point(data = long_data, aes(x = x, y = y), size = 5) +
  geom_line(data = long_data, aes(x = x, y = y, group = x), linewidth = 1) +
  scale_y_continuous(breaks = 1:length(tissue_names), labels = rev(tissue_names)) + 
  theme_minimal() +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12)
  )

###### sharing conditions ######
bcr.contigs$tissue = ifelse(grepl("PBMC", bcr.contigs$orig.ident), "PBMC",
                            ifelse(grepl("P$|P[1-2]$", bcr.contigs$orig.ident), "Peri_tumor", "Tumor"))
# The sharing conditions are defined as follows:
# 1) Tumor-PBMC-Peri_tumor
# 2) Tumor-PBMC
# 3) Tumor-Peri_tumor
# 4) PBMC-Peri_tumor
# 5) Tumor
# 6) PBMC
# 7) Peri_tumor
bcr.contigs.hbc = bcr.contigs[bcr.contigs$antibody_type == "HBc",]
###### Calculate the sharing conditions in HBc group ######
bcr.contigs.hbc.list = split(bcr.contigs.hbc, bcr.contigs.hbc$patient)
total.bcr.hbc.freq = data.frame()
for (patient in names(bcr.contigs.hbc.list)) {
  patient.bcr = bcr.contigs.hbc.list[[patient]]
  patient.bcr.freq = table(patient.bcr$tissue, patient.bcr$new_Clonotype_ID) %>% as.data.frame()
  # Initialize a column to track the sharing conditions
  patient.bcr.freq$sharing_condition = NA
  # Determine sharing conditions for each clonotype
  for (clonotype in unique(patient.bcr.freq$Var2)) {
    clonotype_data = patient.bcr.freq[patient.bcr.freq$Var2 == clonotype, ]
    tissues_present = clonotype_data$Var1[clonotype_data$Freq > 0]
    if (all(c("Tumor", "PBMC", "Peri_tumor") %in% tissues_present)) {
      condition = "Tumor-PBMC-Peri_tumor"
    } else if (all(c("Tumor", "PBMC") %in% tissues_present)) {
      condition = "Tumor-PBMC"
    } else if (all(c("Tumor", "Peri_tumor") %in% tissues_present)) {
      condition = "Tumor-Peri_tumor"
    } else if (all(c("PBMC", "Peri_tumor") %in% tissues_present)) {
      condition = "PBMC-Peri_tumor"
    } else if ("Tumor" %in% tissues_present) {
      condition = "Tumor"
    } else if ("PBMC" %in% tissues_present) {
      condition = "PBMC"
    } else if ("Peri_tumor" %in% tissues_present) {
      condition = "Peri_tumor"
    }
    # Assign the condition to the data frame
    patient.bcr.freq$sharing_condition[patient.bcr.freq$Var2 == clonotype] <- condition
  }
  names(patient.bcr.freq) = c("tissue", "new_Clonotype_ID", "Freq", "sharing_condition")
  patient.bcr.freq$patient = patient
  total.bcr.hbc.freq = rbind(total.bcr.hbc.freq, patient.bcr.freq)
}

# Plot 1: Different sharing conditions freq barplot
sharing.unique.bcr = data.frame()
for (patient in unique(total.bcr.hbc.freq$patient)) {
  share.sub = total.bcr.hbc.freq[total.bcr.hbc.freq$patient == patient,]
  share.sub$new_Clonotype_ID = as.character(share.sub$new_Clonotype_ID)
  share.sub = share.sub[share.sub$Freq != 0,]
  unique.bcr.df = table(share.sub$sharing_condition) %>% as.data.frame()
  unique.bcr.df$patient = patient
  names(unique.bcr.df)[1:2] = c("sharing_condition", "unique_freq")
  sharing.unique.bcr = rbind(sharing.unique.bcr, unique.bcr.df)
}
sharing.unique.bcr = sharing.unique.bcr %>%
  group_by(sharing_condition) %>%
  summarise(unique_freq = sum(unique_freq))
share.df = data.frame(sharing_condition = c("Tumor-PBMC-Peri_tumor", "Tumor-PBMC", "Tumor-Peri_tumor", "PBMC-Peri_tumor", "Tumor", "PBMC", "Peri_tumor"),
                      freq = c(sum(total.bcr.hbc.freq$Freq[total.bcr.hbc.freq$sharing_condition == "Tumor-PBMC-Peri_tumor"]),
                               sum(total.bcr.hbc.freq$Freq[total.bcr.hbc.freq$sharing_condition == "Tumor-PBMC"]),
                               sum(total.bcr.hbc.freq$Freq[total.bcr.hbc.freq$sharing_condition == "Tumor-Peri_tumor"]),
                               sum(total.bcr.hbc.freq$Freq[total.bcr.hbc.freq$sharing_condition == "PBMC-Peri_tumor"]),
                               sum(total.bcr.hbc.freq$Freq[total.bcr.hbc.freq$sharing_condition == "Tumor"]),
                               sum(total.bcr.hbc.freq$Freq[total.bcr.hbc.freq$sharing_condition == "PBMC"]),
                               sum(total.bcr.hbc.freq$Freq[total.bcr.hbc.freq$sharing_condition == "Peri_tumor"])))
share.df = merge(share.df, sharing.unique.bcr, by = "sharing_condition", all.x = T)
share.df$unique_freq[is.na(share.df$unique_freq)] = 0
share.df.m = melt(share.df, id.vars = "sharing_condition")
share.df.m$sharing_condition = factor(share.df.m$sharing_condition, levels = c("Tumor-PBMC-Peri_tumor", "Tumor-PBMC", "Tumor-Peri_tumor", "PBMC-Peri_tumor", "Tumor", "PBMC", "Peri_tumor"))
p2 = ggplot(share.df.m, aes(x = sharing_condition, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = value), vjust = -.1, size = 4, position = position_dodge(.9)) +
  scale_fill_manual(values = two.colors) +
  mytheme2 +
  theme(axis.text.x = element_blank()) +
  labs(title = "",
       x = "",
       y = "",
       fill = "")
p2
# Plot 2: sharing pie chart
share.tumor = share.df[grep("Tumor",share.df$sharing_condition),]
share.tumor$label = paste0(round(share.tumor$unique_freq / sum(share.tumor$unique_freq) * 100, 2), "%")
p3 = ggplot(share.tumor, aes(x = "", y = unique_freq, fill = sharing_condition)) + 
  geom_bar(stat = "identity", width = 1) + 
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 5) + 
  coord_polar("y", start = 0) + 
  theme_void() + 
  scale_fill_manual(values = color.s1) + 
  labs(fill = "") + 
  theme(legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(-.75,-.75,-.75,-.75), "cm"))
myggsave(plot = p3, filename = "rhcc_BCR_similarity/HBc_sharing_conditions_pie.pdf", width = 8, height = 8)
fin.plot = (p2 / p1)
myggsave(plot = fin.plot, filename = "rhcc_BCR_similarity/HBc_sharing_conditions.pdf", width = 10, height = 8)

###### APR ######
orig.input = read.xlsx("../RHCC_all_clonotypes_expand_prop_meg.xlsx")
orig.input$patient = substr(orig.input$orig.ident, 1, 4)
apr.clonotype = orig.input
apr.clonotype = apr.clonotype[apr.clonotype$c_gene != "",]
apr.clonotype = apr.clonotype[apr.clonotype$Freq >= 2 & !is.na(apr.clonotype$Freq),]
apr.clonotype = apr.clonotype[apr.clonotype$expand_group != "Not_expanded" | is.na(apr.clonotype$expand_group),]
apr.clonotype = unique(apr.clonotype$cdr3_nt)

bcr.contigs = read.xlsx("../bcr.contigs.xlsx")
bcr.contigs$patient = substr(bcr.contigs$orig.ident, 1, 4)
bcr.contigs$new_Clonotype_ID = paste0(bcr.contigs$patient, "_", bcr.contigs$new_Clonotype_ID)
bcr.contigs.apr = bcr.contigs[bcr.contigs$cdr3_nt %in% apr.clonotype,]
bcr.contigs.apr$tissue = ifelse(grepl("PBMC", bcr.contigs.apr$orig.ident), "PBMC",
                            ifelse(grepl("P$|P[1-2]$", bcr.contigs.apr$orig.ident), "Peri_tumor", "Tumor"))
unique(bcr.contigs.apr$tissue)

bcr.contigs.apr.list = split(bcr.contigs.apr, bcr.contigs.apr$patient)
total.bcr.apr.freq = data.frame()
for (patient in names(bcr.contigs.apr.list)) {
  patient.bcr = bcr.contigs.apr.list[[patient]]
  patient.bcr.freq = table(patient.bcr$tissue, patient.bcr$new_Clonotype_ID) %>% as.data.frame()
  # Initialize a column to track the sharing conditions
  patient.bcr.freq$sharing_condition = NA
  # Determine sharing conditions for each clonotype
  for (clonotype in unique(patient.bcr.freq$Var2)) {
    clonotype_data = patient.bcr.freq[patient.bcr.freq$Var2 == clonotype, ]
    tissues_present = clonotype_data$Var1[clonotype_data$Freq > 0]
    if (all(c("Tumor", "PBMC", "Peri_tumor") %in% tissues_present)) {
      condition = "Tumor-PBMC-Peri_tumor"
    } else if (all(c("Tumor", "PBMC") %in% tissues_present)) {
      condition = "Tumor-PBMC"
    } else if (all(c("Tumor", "Peri_tumor") %in% tissues_present)) {
      condition = "Tumor-Peri_tumor"
    } else if (all(c("PBMC", "Peri_tumor") %in% tissues_present)) {
      condition = "PBMC-Peri_tumor"
    } else if ("Tumor" %in% tissues_present) {
      condition = "Tumor"
    } else if ("PBMC" %in% tissues_present) {
      condition = "PBMC"
    } else if ("Peri_tumor" %in% tissues_present) {
      condition = "Peri_tumor"
    }
    # Assign the condition to the data frame
    patient.bcr.freq$sharing_condition[patient.bcr.freq$Var2 == clonotype] <- condition
  }
  names(patient.bcr.freq) = c("tissue", "new_Clonotype_ID", "Freq", "sharing_condition")
  patient.bcr.freq$patient = patient
  total.bcr.apr.freq = rbind(total.bcr.apr.freq, patient.bcr.freq)
}

# Plot 1: Different sharing conditions freq barplot
sharing.unique.bcr = data.frame()
for (patient in unique(total.bcr.apr.freq$patient)) {
  share.sub = total.bcr.apr.freq[total.bcr.apr.freq$patient == patient,]
  share.sub$new_Clonotype_ID = as.character(share.sub$new_Clonotype_ID)
  share.sub = share.sub[share.sub$Freq != 0,]
  unique.bcr.df = table(share.sub$sharing_condition) %>% as.data.frame()
  unique.bcr.df$patient = patient
  names(unique.bcr.df)[1:2] = c("sharing_condition", "unique_freq")
  sharing.unique.bcr = rbind(sharing.unique.bcr, unique.bcr.df)
}
sharing.unique.bcr = sharing.unique.bcr %>%
  group_by(sharing_condition) %>%
  summarise(unique_freq = sum(unique_freq))
share.df = data.frame(sharing_condition = c("Tumor-PBMC-Peri_tumor", "Tumor-PBMC", "Tumor-Peri_tumor", "PBMC-Peri_tumor", "Tumor", "PBMC", "Peri_tumor"),
                      freq = c(sum(total.bcr.apr.freq$Freq[total.bcr.apr.freq$sharing_condition == "Tumor-PBMC-Peri_tumor"]),
                               sum(total.bcr.apr.freq$Freq[total.bcr.apr.freq$sharing_condition == "Tumor-PBMC"]),
                               sum(total.bcr.apr.freq$Freq[total.bcr.apr.freq$sharing_condition == "Tumor-Peri_tumor"]),
                               sum(total.bcr.apr.freq$Freq[total.bcr.apr.freq$sharing_condition == "PBMC-Peri_tumor"]),
                               sum(total.bcr.apr.freq$Freq[total.bcr.apr.freq$sharing_condition == "Tumor"]),
                               sum(total.bcr.apr.freq$Freq[total.bcr.apr.freq$sharing_condition == "PBMC"]),
                               sum(total.bcr.apr.freq$Freq[total.bcr.apr.freq$sharing_condition == "Peri_tumor"])))
share.df = merge(share.df, sharing.unique.bcr, by = "sharing_condition", all.x = T)
share.df$unique_freq[is.na(share.df$unique_freq)] = 0
share.df.m = melt(share.df, id.vars = "sharing_condition")
share.df.m$sharing_condition = factor(share.df.m$sharing_condition, levels = c("Tumor-PBMC-Peri_tumor", "Tumor-PBMC", "Tumor-Peri_tumor", "PBMC-Peri_tumor", "Tumor", "PBMC", "Peri_tumor"))
p2 = ggplot(share.df.m, aes(x = sharing_condition, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = value), vjust = -.1, size = 4, position = position_dodge(.9)) +
  scale_fill_manual(values = two.colors) +
  mytheme2 +
  theme(axis.text.x = element_blank()) +
  labs(title = "",
       x = "",
       y = "",
       fill = "") + 
  scale_y_break(c(650, 4100), scales = .2)
p2
# Plot 2: sharing pie chart
share.tumor = share.df[grep("Tumor",share.df$sharing_condition),]
share.tumor$label = paste0(round(share.tumor$unique_freq / sum(share.tumor$unique_freq) * 100, 2), "%")
p3 = ggplot(share.tumor, aes(x = "", y = unique_freq, fill = sharing_condition)) + 
  geom_bar(stat = "identity", width = 1) + 
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 5) + 
  coord_polar("y", start = 0) + 
  theme_void() + 
  scale_fill_manual(values = color.s1) + 
  labs(fill = "") + 
  theme(legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(-.75,-.75,-.75,-.75), "cm"))
myggsave(plot = p3, filename = "rhcc_BCR_similarity/APR_sharing_conditions_pie.pdf", width = 8, height = 8)

library(patchwork)
fin.plot = (p2 / p1)
myggsave(plot = fin.plot, filename = "rhcc_BCR_similarity/APR_sharing_conditions.pdf", width = 10, height = 8)

###### All antibody ######
bcr.contigs = read.xlsx("../bcr.contigs.xlsx")
bcr.contigs$patient = substr(bcr.contigs$orig.ident, 1, 4)
bcr.contigs$new_Clonotype_ID = paste0(bcr.contigs$patient, "_", bcr.contigs$new_Clonotype_ID)
bcr.contigs$tissue = ifelse(grepl("PBMC", bcr.contigs$orig.ident), "PBMC",
                            ifelse(grepl("P$|P[1-2]$", bcr.contigs$orig.ident), "Peri_tumor", "Tumor"))
bcr.contigs = bcr.contigs[bcr.contigs$chain == "IGH",]
bcr.contigs = bcr.contigs[bcr.contigs$c_gene != "",]

bcr.contigs.list = split(bcr.contigs, bcr.contigs$patient)
total.bcr.freq = data.frame()
for (patient in names(bcr.contigs.list)) {
  patient.bcr = bcr.contigs.list[[patient]]
  patient.bcr.freq = table(patient.bcr$tissue, patient.bcr$new_Clonotype_ID) %>% as.data.frame()
  # Initialize a column to track the sharing conditions
  patient.bcr.freq$sharing_condition = NA
  # Determine sharing conditions for each clonotype
  for (clonotype in unique(patient.bcr.freq$Var2)) {
    clonotype_data = patient.bcr.freq[patient.bcr.freq$Var2 == clonotype, ]
    tissues_present = clonotype_data$Var1[clonotype_data$Freq > 0]
    if (all(c("Tumor", "PBMC", "Peri_tumor") %in% tissues_present)) {
      condition = "Tumor-PBMC-Peri_tumor"
    } else if (all(c("Tumor", "PBMC") %in% tissues_present)) {
      condition = "Tumor-PBMC"
    } else if (all(c("Tumor", "Peri_tumor") %in% tissues_present)) {
      condition = "Tumor-Peri_tumor"
    } else if (all(c("PBMC", "Peri_tumor") %in% tissues_present)) {
      condition = "PBMC-Peri_tumor"
    } else if ("Tumor" %in% tissues_present) {
      condition = "Tumor"
    } else if ("PBMC" %in% tissues_present) {
      condition = "PBMC"
    } else if ("Peri_tumor" %in% tissues_present) {
      condition = "Peri_tumor"
    }
    # Assign the condition to the data frame
    patient.bcr.freq$sharing_condition[patient.bcr.freq$Var2 == clonotype] <- condition
  }
  names(patient.bcr.freq) = c("tissue", "new_Clonotype_ID", "Freq", "sharing_condition")
  patient.bcr.freq$patient = patient
  total.bcr.freq = rbind(total.bcr.freq, patient.bcr.freq)
}

# Plot 1: Different sharing conditions freq barplot
sharing.unique.bcr = data.frame()
for (patient in unique(total.bcr.freq$patient)) {
  share.sub = total.bcr.freq[total.bcr.freq$patient == patient,]
  share.sub$new_Clonotype_ID = as.character(share.sub$new_Clonotype_ID)
  share.sub = share.sub[share.sub$Freq != 0,]
  unique.bcr.df = table(share.sub$sharing_condition) %>% as.data.frame()
  unique.bcr.df$patient = patient
  names(unique.bcr.df)[1:2] = c("sharing_condition", "unique_freq")
  sharing.unique.bcr = rbind(sharing.unique.bcr, unique.bcr.df)
}
sharing.unique.bcr = sharing.unique.bcr %>%
  group_by(sharing_condition) %>%
  summarise(unique_freq = sum(unique_freq))
share.df = data.frame(sharing_condition = c("Tumor-PBMC-Peri_tumor", "Tumor-PBMC", "Tumor-Peri_tumor", "PBMC-Peri_tumor", "Tumor", "PBMC", "Peri_tumor"),
                      freq = c(sum(total.bcr.freq$Freq[total.bcr.freq$sharing_condition == "Tumor-PBMC-Peri_tumor"]),
                               sum(total.bcr.freq$Freq[total.bcr.freq$sharing_condition == "Tumor-PBMC"]),
                               sum(total.bcr.freq$Freq[total.bcr.freq$sharing_condition == "Tumor-Peri_tumor"]),
                               sum(total.bcr.freq$Freq[total.bcr.freq$sharing_condition == "PBMC-Peri_tumor"]),
                               sum(total.bcr.freq$Freq[total.bcr.freq$sharing_condition == "Tumor"]),
                               sum(total.bcr.freq$Freq[total.bcr.freq$sharing_condition == "PBMC"]),
                               sum(total.bcr.freq$Freq[total.bcr.freq$sharing_condition == "Peri_tumor"])))
share.df = merge(share.df, sharing.unique.bcr, by = "sharing_condition", all.x = T)
share.df$unique_freq[is.na(share.df$unique_freq)] = 0
share.df.m = melt(share.df, id.vars = "sharing_condition")
share.df.m$sharing_condition = factor(share.df.m$sharing_condition, levels = c("Tumor-PBMC-Peri_tumor", "Tumor-PBMC", "Tumor-Peri_tumor", "PBMC-Peri_tumor", "Tumor", "PBMC", "Peri_tumor"))
p2 = ggplot(share.df.m, aes(x = sharing_condition, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = value), vjust = -.1, size = 4, position = position_dodge(.9)) +
  scale_fill_manual(values = two.colors) +
  mytheme2 +
  theme(axis.text.x = element_blank()) +
  labs(title = "",
       x = "",
       y = "",
       fill = "") + 
  scale_y_break(c(600, 3000), scales = 1.5) + 
  scale_y_break(c(4200, 17500), scales = 2)
p2
# Plot 2: sharing pie chart
share.tumor = share.df[grep("Tumor",share.df$sharing_condition),]
share.tumor$label = paste0(round(share.tumor$unique_freq / sum(share.tumor$unique_freq) * 100, 2), "%")
p3 = ggplot(share.tumor, aes(x = "", y = unique_freq, fill = sharing_condition)) + 
  geom_bar(stat = "identity", width = 1) + 
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 5) + 
  coord_polar("y", start = 0) + 
  theme_void() + 
  scale_fill_manual(values = color.s1) + 
  labs(fill = "") + 
  theme(legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(-.75,-.75,-.75,-.75), "cm"))
myggsave(plot = p3, filename = "rhcc_BCR_similarity/All_sharing_conditions_pie.pdf", width = 8, height = 8)

fin.plot = (p2 / p1)
myggsave(plot = fin.plot, filename = "rhcc_BCR_similarity/All_sharing_conditions.pdf", width = 10, height = 8)


######Fig. 4E######
patient = "RH05B"
tls.selected = c("STRH05B_10","STRH05B_11",
                 "STRH05B_12","STRH05B_2",
                 "STRH05B_21","STRH05B_25",
                 "STRH05B_30","STRH05B_4",
                 "STRH05B_6","STRH05B_9")

hbc.bcr = readRDS("空转互作_三代测序/hbcab_map_clonefam_v3.rds")
hbc.bcr.05 = hbc.bcr[hbc.bcr$patient == "RH05", ]
hbc.bcr.05 = hbc.bcr.05[!is.na(hbc.bcr.05$Bin20),]
hbc.bcr.05 = hbc.bcr.05[grep("^m", hbc.bcr.05$clone_subid),]

stRNA = readRDS(paste0("../data/slice7.results.renew/reanno_cjy_T_B_mye_endo_fib/", 
                       dir("../data/slice7.results.renew/reanno_cjy_T_B_mye_endo_fib/") %>% .[grep(patient, .)]))

meta.data = stRNA@meta.data
meta.data$bins = gsub("20_",".",meta.data$bins)
meta.data$antibody_type = ifelse(meta.data$TLS_raw %in% tls.selected, 
                                 "TLS", "Others")
meta.data$antibody_type[meta.data$Bin_Region == "Tumor_capsule"] = "Boundary"
meta.data$antibody_type[meta.data$Bin_Region == "Paratumor_side_of_Margin_area"] = "Peri-tumor"

match.bins = match(meta.data$bins, hbc.bcr.05$Bin20)
meta.data$antibody_type[!is.na(match.bins)] = hbc.bcr.05$clone_subid[match.bins[!is.na(match.bins)]]
stRNA@meta.data = meta.data
plot_item = "antibody_type"
prefix = paste0("空转互作_三代测序/抗体Dimplot_v3/", patient,"_with_TLS")
if (plot_item %in% colnames(stRNA@meta.data)) {
  dimplot_inhouse_R(obj = stRNA, plot_item = plot_item, prefix = prefix,
                    tmp_color = c("#000000",ggsci::pal_igv()(9),"#F8F5EC","#87CEEBB3"))
} else {
  tryCatch({
    print(paste0(patient, " does not have ", plot_item, " information in the meta.data"))
  }, error = function(e) {
    print(paste0("Error occurred: ", e))
  })
}

######Fig. 4f######


######Fig. 4i######
hbv.core.id = readRDS("st_ter_BCR_new/talent_version2/talent_result.rds")
load("st_ter_BCR_new/st_bcr_h.RData")
st = readRDS("st_ter_BCR_new/st_T_B.rds")
meta = st@meta.data
shm.df = bcr.df[,c("sequence_id","sequence","shm","clone_id")]
hbv.core.id = merge(hbv.core.id, shm.df, by = "sequence_id")
bcr.df$hbv_group = ifelse(bcr.df$sequence_id %in% hbv.core.id$sequence_id, 
                          "HBV_core", "Other_BCR")
bcr = bcr.df[,c("BinID","hbv_group")]
names(bcr)[1] = "BinID"
bcr$rows = rownames(bcr)
rows.keep = c()
for (bin in unique(bcr$BinID)) {
  bin.df = bcr[bcr$BinID == bin,,drop = F]
  if (nrow(bin.df) > 1) {
    if (all(bin.df$hbv_group == "Other_BCR")) {
      row.index = bin.df$rows[1]
      rows.keep = c(rows.keep, row.index)
    } else {
      row.index = bin.df$rows[bin.df$hbv_group == "HBV_core"]
      rows.keep = c(rows.keep, row.index)
    }
  } else {
    rows.keep = c(rows.keep, bin.df$rows)
  }
}
bcr = bcr[as.numeric(rows.keep),]
st.sub = subset(st, Cellsubtype == "Plasma")
st.sub$raw_rownames = rownames(st.sub@meta.data)
tmp = st.sub@meta.data
plasma.meta = merge(tmp, bcr, by = "BinID", all.x = T)
plasma.meta = plasma.meta[!duplicated(plasma.meta$bins),]
rownames(plasma.meta) = plasma.meta$raw_rownames
plasma.meta = plasma.meta[order(factor(plasma.meta$raw_rownames,
                                                 levels = colnames(st.sub))),]
st.sub@meta.data = plasma.meta
plasma.degs = FindMarkers(st.sub, ident.1 = "HBV_core", ident.2 = "Other_BCR", 
                               group.by = "hbv_group")
degs = plasma.degs[order(plasma.degs$p_val),]
degs$gene = rownames(degs)
degs$logP = -log10(degs$p_val)
degs$Group = rep("not significantly changed", nrow(degs))
degs$Group[which(degs$avg_log2FC > .25 & degs$p_val < 0.05)] = "up-regulated"
degs$Group[which(degs$avg_log2FC < -.25 & degs$p_val < 0.05)] = "down-regulated"
degs$Label = ""
up.genes <- head(degs$gene[which(degs$Group == "up-regulated")],10)
down.genes <- head(degs$gene[which(degs$Group == "down-regulated")],10)
deg.top10.genes <- c(as.character(up.genes),as.character(down.genes))
degs$Label[match(deg.top10.genes,degs$gene)] <- deg.top10.genes
options(ggrepel.max.overlaps = Inf)
ggscatter(degs,x="avg_log2FC",y="logP",
          color = "Group",
          palette = c(color.s1[2],"#BBBBBB",color.s1[1]),
          size = 1,
          label = degs$Label,
          font.label = 12,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)",) + 
  guides(color = guide_legend(ncol = 3)) + 
  geom_hline(yintercept = 1.30,linetype="dashed")+
  geom_vline(xintercept = c(-.25,.25),linetype="dashed") + 
  mytheme2

######Fig. 4j######
require(clusterProfiler)
require(org.Hs.eg.db)
k = keys(org.Hs.eg.db, keytype = "SYMBOL")
glist = AnnotationDbi::select(org.Hs.eg.db, keys = k, columns = c("ENTREZID"),keytype = "SYMBOL")
names(glist)[1] = "gene"
degs = merge(degs, glist, by = "gene")
degs.eff = degs[degs$Group != "not significantly changed",]
go = enrichGO(degs.eff$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", 
              pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, readable = T)
barplot(go, split = "ONTOLOGY", showCategory = 5, color = "pvalue") + 
  facet_grid(ONTOLOGY~., scales = "free")


#####Ex Figure1#####
######Ex Fig.1a######
library(scRNAtoolVis)
library(ggsci)
scRNA = readRDS("../data/scRNAseq_data/RHTumor_clean.rds")
scRNA$treatment = factor(scRNA$treatment, levels = c("Pre", "Post"))
clusterCornerAxes(object = scRNA, reduction = "umap", noSplit = F,
                  groupFacet = "treatment",
                  clusterCol = "Majortype", aspect.ratio = 1, pSize = .1,
                  axes = "one", relLength = .3) + scale_color_npg()
myggsave(plot = last_plot(), "Visualization/Ex_Figure1/Ex_Fig_1a/RHCC_Tumor_Dimplot.pdf", width = 10, height = 6)

######Ex Fig.1b######
clusterCornerAxes(object = scRNA, reduction = "umap", 
                  clusterCol = "orig.ident")

#####Ex Figure3#####
######Ex Fig.3a######
tcell = readRDS("元老师/Tcell_TCR.rds")
DimPlot(tcell, reduction = "umap", label = TRUE,group.by = "Sub_cluster",pt.size = 0.0001,raster = F)+scale_color_d3("category20c")
######Ex Fig.3b######
tcell$TCR_existence = ifelse(is.na(tcell$CTgene), "No TCR detected", "With TCR")
DimPlot(tcell, reduction = "umap", label = F,group.by = "TCR_existence",
        pt.size = 0.0001,raster = F)+scale_color_manual(values = c("grey","#6CA6CD"))
tcell$TCR_expansion = ifelse(is.na(tcell$CTgene), "No TCR detected",
                             ifelse(tcell$Frequency >= 5, "Hypterexpanded (x ≥ 5)",
                             ifelse(tcell$Frequency < 2, "Single (x < 2)", "Medium (2 ≤ x < 5)")))
DimPlot(tcell, reduction = "umap", label = F,group.by = "TCR_expansion",
        pt.size = 0.0001,raster = F)

######Ex Fig.3c######
plot(density(tcr$Freq5.ratio,width=0.6),xlim=c(0,0.8),main="TCR.Freq5",xlab="TCR.Freq5",lwd=2,axes=F)
abline(v=0.1)
axis(side = 1,at=seq(0:0.8,by=0.1))
axis(side=2)
box()
plot(density(tcr$Freq2.ratio,width=0.15),xlim=c(0,0.6),main="TCR.Freq2",xlab="TCR.Freq2",lwd=2)
abline(v=0.15)
axis(side = 1,at=0.15)

######Ex Fig.3f######
tfh=subset(tcell,Sub_cluster%in%c("CD4_Tht_CXCL13","CD8_Tex_CXCL13"))
tfh$group='Non_responder'
tfh$group[tfh$Expansion_group=='LE_NP']='B_responder'
tfh$group[tfh$Expansion_group=='HE_NP']='T_responder'
tfh$group=paste0(tfh$group,"_",tfh$treatment)
tfh$group2=paste0(tfh$group,"_",tfh$Sub_cluster)
tfh$pd1=tfh@assays$RNA@data["PDCD1",]
tfh$pd1=tfh@assays$RNA@data["PDCD1",]
tfh$PDCD1="negative"
tfh$PDCD1[tfh$pd1>0]='positive'
tfh_pre=subset(tfh,treatment=="Pre")
tfh_post=subset(tfh,treatment=='Post')

tfh_pre$group2=factor(tfh_pre$group2,c("B_responder_Pre_CD4_Tht_CXCL13",'B_responder_Pre_CD8_Tex_CXCL13',
                                       'T_responder_Pre_CD4_Tht_CXCL13','T_responder_Pre_CD8_Tex_CXCL13',
                                       'Non_responder_Pre_CD4_Tht_CXCL13','Non_responder_Pre_CD8_Tex_CXCL13'))

tfh_post$group2=factor(tfh_post$group2,c("B_responder_Post_CD4_Tht_CXCL13",'B_responder_Post_CD8_Tex_CXCL13',
                                         'T_responder_Post_CD4_Tht_CXCL13','T_responder_Post_CD8_Tex_CXCL13',
                                         'Non_responder_Post_CD4_Tht_CXCL13','Non_responder_Post_CD8_Tex_CXCL13'))

porportion=as.data.frame(table(tfh_pre$PDCD1,tfh_pre$group2))
porportion=as.data.frame(table(tfh_post$PDCD1,tfh_post$group2))

porp=porportion %>% 
  pivot_wider(names_from = "Var1",
              values_from = "Freq")
porp$ratio=porp$positive/(porp$negative+porp$positive)
porp$Var2=unique(porportion$Var2)
ggdotchart(porp, x = "Var2", y = "ratio",
           color = "Var2",                            
           #palette = colorvec,
           #sorting = "ascending",                        
           add = "segments",                           
           ggtheme = theme_pubr(),                     
           xlab="",
           ylab = "PDCD1 positive ratio",
           add.params = list(color = "lightgray", size = 2),
           dot.size = 9
)+ theme(axis.text.x = element_text(angle = 45, hjust = 0.9,vjust = 0.9),legend.position = 'none')+
  scale_color_manual(values = c(my.colors[2],my.colors[5],
                                my.colors[6],my.colors[4],
                                my.colors[1],my.colors[3]))

######Ex Fig.3g######
tcell$TCR="nonexpanded"
tcell$TCR[tcell$Frequency>5]='expanded'

porportion = table(tcell$TCR, tcell$orig.ident, tcell$Sub_cluster) %>% as.data.frame()
porportion[,2] = as.character(porportion[,2])
porportion[,3] = as.character(porportion[,3])
porportion.split = split(porportion, f = porportion$Var3)
porportion.split = lapply(porportion.split, FUN = function(x){
  nrows = nrow(x)
  sample_number = length(unique(x[,2]))
  for (i in seq(1,nrows, by = nrows/sample_number)){
    sum = sum(x$Freq[i:(i+((nrows/sample_number)-1))])
    x$Proportion[i:(i+((nrows/sample_number)-1))] = x$Freq[i:(i+((nrows/sample_number)-1))] / sum
  }
  return(x)
})

porportion = Reduce(function(x,y) rbind(x,y), porportion.split)
porp2=porportion[porportion$Var1=='expanded',]
colnames(porp2)=c("TCR",'orig.ident','Sub_cluster','Freq','Ratio')
pheno = tcell@meta.data[,c('orig.ident','Expansion_group','treatment')]
pheno = pheno[!duplicated(pheno$orig.ident),]
porp3=merge(porp2,pheno,by="orig.ident")

porp4=porp3[porp3$Expansion_group%in%c("HE_NP"),]
porp4$group=paste0(porp4$treatment,porp4$Sub_cluster)
porp.henp = porp4[porp4$treatment == "Post",]
porp.henp.split = split(porp.henp, f = porp.henp$group)
median = c()
for (i in names(porp.henp.split)){
  df = porp.henp.split[[i]]
  m = median(df$Ratio)
  names(m) = i
  median = c(median, m)
}
median = sort(median)
cell.order = names(median)

combinations_list <- as.list(data.frame(combn(cell.order, 2)))

#pre
porp.henp = porp4[porp4$treatment == "Pre",]
porp.henp.split = split(porp.henp, f = porp.henp$group)
median = c()
for (i in names(porp.henp.split)){
  df = porp.henp.split[[i]]
  m = median(df$Ratio)
  names(m) = i
  median = c(median, m)
}
median = sort(median)
cell.order2 = names(median)
cell.order2 <- gsub("Post", "Pre", cell.order, ignore.case = TRUE)

combinations_list <- as.list(data.frame(combn(cell.order2, 2)))
final.order=c(cell.order,cell.order2)

ggplot(porp4,aes(x=factor(group,levels = final.order),y=Ratio,fill=factor(treatment,levels = c('Pre','Post'))))+
  geom_boxplot(outlier.colour = NA)+geom_jitter(alpha=0.3,position = position_jitterdodge())+
  theme_classic()+
  labs(x = "", y = "expanded clonotype ratio", title = "High-expanded Late/no-recurrence") +
  stat_compare_means(label  = "p.format",
                     method = "anova", size = 3, bracket.size = p.lwd)+
  scale_fill_manual(values = two.colors)+theme(axis.text.x = element_text(angle = 90))+NoLegend()

porp4=porp3[porp3$Expansion_group%in%c("HE_P"),]
porp4$group=paste0(porp4$treatment,porp4$Sub_cluster)
porp.hep = porp4[porp4$treatment == "Post",]
porp.hep.split = split(porp.hep, f = porp.hep$group)
median = c()
for (i in names(porp.hep.split)){
  df = porp.hep.split[[i]]
  m = median(df$Ratio)
  names(m) = i
  median = c(median, m)
}
median = sort(median)
cell.order = names(median)
combinations_list <- as.list(data.frame(combn(cell.order, 2)))

porp.hep = porp4[porp4$treatment == "Pre",]
porp.hep.split = split(porp.hep, f = porp.hep$group)
median = c()
for (i in names(porp.hep.split)){
  df = porp.hep.split[[i]]
  m = median(df$Ratio)
  names(m) = i
  median = c(median, m)
}
median = sort(median)
cell.order2 <- gsub("Post", "Pre", cell.order, ignore.case = TRUE)

combinations_list <- as.list(data.frame(combn(cell.order2, 2)))
final.order=c(cell.order,cell.order2)

ggplot(porp4,aes(x=factor(group,levels = final.order),y=Ratio,fill=factor(treatment,levels = c('Pre','Post'))))+
  geom_boxplot(outlier.colour = NA)+geom_jitter(alpha=0.3,position = position_jitterdodge())+
  theme_classic()+
  labs(x = "", y = "expanded clonotype ratio", title = "High-expanded Early-recurrence") +
  stat_compare_means(label  = "p.format",
                     method = "anova", size = 3, bracket.size = p.lwd)+
  scale_fill_manual(values = two.colors)+theme(axis.text.x = element_text(angle = 90))+NoLegend()


#####Ex Figure5#####
######Ex Fig.5a######
Tex=subset(tcell,Sub_cluster=="CD8_Tex_CXCL13")
Tex_he=subset(Tex,Expansion_normalized=='HighE')
Tex_he$pd1=Tex_he@assays$RNA@data["PDCD1",]
Tex_he$group = factor(Tex_he$group, levels = c("Pre_NP", "Post_NP",
                                               "Pre_P", "Post_P"))
Tex_he_meta = Tex_he@meta.data
ggplot(Tex_he_meta, aes(x = group, y = pd1, color = group)) + 
  geom_violin(aes(fill = group, color = group), scale = "width", trim = F) +
  geom_boxplot(aes(fill = group, color = group),
               width = .05, outlier.shape = NA, color = "black", 
               linewidth = violin.box.lwd)+
  labs(x = "", y = "Expression level", title = "PDCD1 (Tex in HighE)") +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c('Pre_NP','Post_NP'),c("Pre_P","Post_P")),
                     method = "wilcox.test", size = 5, bracket.size = p.lwd)+
  scale_fill_manual(values = my.colors1) +
  scale_color_manual(values = my.colors1) + 
  NoLegend() + RotatedAxis() + mytheme2

######Ex Fig.5d######
tc.set = readRDS("../data/scRNAseq_data/RH_Tumor_Tcells_3groups.rds")
tex.set = subset(tc.set, Sub_cluster == "CD8_Tex_CXCL13")
pheno = tex.set@meta.data
pheno$barcode = rownames(pheno) 
pheno = pheno[,c("orig.ident","treatment","progress","Expansion_normalized","barcode")]
pheno2 = pheno[,-which(names(pheno) == ("barcode"))]
pheno2 = pheno2[!duplicated(pheno2$orig.ident),]
# Get the expression data of PRF1, NKG7, GZMK, PDCD1, LAG3, HAVCR2, NR4A1, JUN, FOS
genes = c("PRF1","NKG7","GZMK","PDCD1","LAG3","HAVCR2","NR4A1","JUN","FOS")
expr = GetAssayData(tex.set, assay = "RNA", slot = "data") %>% as_matrix()
expr = expr[genes,]
expr = expr %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>%
  pivot_longer(cols = -genes, names_to = "barcode", values_to = "expr")
expr = merge(expr, pheno, by = "barcode")
# Keep only HighE Expansion_normalized
expr = expr[expr$Expansion_normalized == "HighE",]
# Keep only Post treatment
expr = expr[expr$treatment == "Post",]
expr$progress = ifelse(expr$progress == "Progress", "NR", "R")
# Calculate the mean expression of each gene in each orig.ident
expr.mean = expr %>% 
  group_by(orig.ident, genes) %>% 
  summarise(mean_expr = mean(expr))
expr.mean = merge(expr.mean, pheno2, by = "orig.ident")
expr.mean$progress = ifelse(expr.mean$progress == "Progress", "NR", "R")
# Plot each gene plot in progress
for (gene in unique(expr.mean$genes)) {
  gene.expr = expr.mean[expr.mean$genes == gene,]
  p1 = ggplot(gene.expr, aes(x = factor(progress, levels = c("R", "NR")), y = mean_expr, color = progress)) + 
    # geom_violin(aes(fill = progress, color = progress), scale = "width", trim = F) + 
    geom_boxplot(aes(fill = progress), position = position_dodge(.9), 
                 width = .6, outlier.shape = NA, color = "black", linewidth = violin.box.lwd) +
    # stat_summary(fun.y = mean, geom = "point", size = 2) + 
    # stat_summary(fun.data = mean_se, geom = "errorbar") + 
    geom_jitter(stroke = .3,
                position = position_jitterdodge(dodge.width = .9, jitter.width = 0), 
                size = 2, aes(fill = progress), color = "black",
                alpha = .6) +
    stat_compare_means(comparisons = list(c("R","NR")), label = "p.signif",
                       method = "t.test", size = 5, bracket.size = p.lwd) +
    labs(x = "", y = "Expression", title = gene) +
    scale_color_manual(values = rev(two.colors)) +
    scale_fill_manual(values = rev(two.colors)) +
    mytheme2 + theme(legend.position = "none")
  gene.expr = expr[expr$genes == gene,]
  p2 = ggplot(gene.expr, aes(x = factor(progress, levels = c("R", "NR")), y = expr, color = progress)) + 
    geom_violin(aes(fill = progress, color = progress), scale = "width", trim = F) +
    geom_boxplot(aes(fill = progress), position = position_dodge(.9), 
                 width = .05, outlier.shape = NA, color = "black", linewidth = violin.box.lwd) +
    # stat_summary(fun.y = mean, geom = "point", size = 2) + 
    # stat_summary(fun.data = mean_se, geom = "errorbar") + 
    stat_compare_means(comparisons = list(c("R","NR")), label = "p.signif",
                       method = "t.test", size = 5, bracket.size = p.lwd) +
    labs(x = "", y = "Expression", title = gene) +
    scale_color_manual(values = rev(two.colors)) +
    scale_fill_manual(values = rev(two.colors)) +
    mytheme2 + theme(legend.position = "none")
  myggsave(paste0("Visualization/Ex_Figure5/Ex_Fig_5d/T_genes/",gene,"_sample_level_expression_in_progress.pdf"), p1, width = 4, height = 6)
  myggsave(paste0("Visualization/Ex_Figure5/Ex_Fig_5d/T_genes/",gene,"_expression_in_progress.pdf"), p2, width = 4, height = 6)
}

######Ex Fig.5b######
genelist = read.csv("元老师/肿瘤细胞与TLS相关基因集.csv")
cd8 = subset(tcell,Sub_cluster=="CD8_Tex_CXCL13")
cd8 = subset(cd8,Expansion_normalized=='HighE')
cd8=subset(cd8,group%in%c("Post_NP",'Post_P'))
cd8 = AddModuleScore(cd8,features  = data.frame(genelist$Cytotoxicity),name = "cytotoxicity")
ggplot(cd8@meta.data,aes(x=group,y=cytotoxicity1,fill=group)) + 
  geom_violin(aes(color = group))+
  geom_boxplot(width = .1, outlier.shape = NA, color = "black", 
               linewidth = violin.box.lwd)+
  theme_classic()+
  labs(x = "", y = "Cell score", title = "Tex Cytotoxic score") +
  stat_compare_means(label = "p.signif",
                     method = "t.test", size = 5, bracket.size = p.lwd)+
  scale_fill_manual(values = two.colors)+
  scale_colour_manual(values = two.colors)+
  NoLegend()+RotatedAxis() + mytheme2


#####Ex Figure9#####
######Ex Fig.9b######
data=FetchData(tumor,vars=c('group','treatment','REACTOME_COMPLEMENT_CASCADE'))
median = data %>% group_by(group) %>% summarise(median = median(REACTOME_COMPLEMENT_CASCADE, na.rm = T)) %>% ungroup()
p=ggplot(data,aes(group,REACTOME_COMPLEMENT_CASCADE,fill=group))
p+geom_violin(aes(color = group, fill = group))+
  geom_boxplot(width = .1, outlier.shape = NA, color = "black", linewidth = violin.box.lwd)+
  theme_classic()+ 
  geom_line(data = median[1:2,], mapping = aes(x = group, y = median, group = 1), 
            linewidth = .1, color = "black", linetype = "dashed") +
  geom_line(data = median[3:4,], mapping = aes(x = group, y = median, group = 1), 
            linewidth = .1, color = "black", linetype = "dashed") +
  geom_line(data = median[5:6,], mapping = aes(x = group, y = median, group = 1), 
            linewidth = .1, color = "black", linetype = "dashed") +
  stat_compare_means(label = "p.signif",comparisons = list(c('T_responder_Pre','T_responder_Post'),c('Non_responder_Pre','Non_responder_Post'),
                                                           c('B_responder_Pre','B_responder_Post')),
                     method = "t.test",size = 5, bracket.size = p.lwd)+
  scale_fill_manual(values = c(my.colors1[1:2],
                               my.colors1[5:6],
                               my.colors1[3:4])) + 
  scale_colour_manual(values = c(my.colors1[1:2],
                                 my.colors1[5:6],
                                 my.colors1[3:4]))+ 
  RotatedAxis() + NoLegend() + mytheme2
