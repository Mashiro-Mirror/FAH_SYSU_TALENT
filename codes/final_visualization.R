#####Loading the Library#####
library(tidyverse)
library(ggpubr)
library(ggsci)
library(Seurat)
library(rstatix)
library(Hmisc)
library(circlize)
library(corrplot)
library(data.table)
library(RColorBrewer)
library(msigdbr)
library(patchwork)

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
two.colors = c("#87E9C3","#FCCF7C")
three.colors = c()
# four.colors = c("#FDC7CD","#EE84A8","#C2D7F3","#778CCC")
four.colors = c("#F8C4AC","#DD8385","#C4C1DE","#9A98C9")
five.colors = c()
color.s1 = c("#E64B35B2","#00A087B2","#3C5488B2","#4DBBD5B2","#F39B7FB2",
             "#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2")
color = brewer.pal(12,"Set3")

## linewidth
violin.box.lwd = .1
p.lwd = .1

#####Fig2 & S2#####
######Fig. 2B#####
#Match pre and post treatment samples
Tcell_TCR = readRDS("../data/single_cell_RNA/Tcell_TCR.rds")
Expansion_group=Tcell_TCR@meta.data[,c("orig.ident",'cellbarcode',"Frequency","CTnt","uniqueCTnt","treatment")]
TCR_pair=read.csv("../data/RHCC_tumor_treatmentpair.csv")
Expansion_group2=merge(Expansion_group2,TCR_pair,by.x = "orig.ident")
#Delete the one that doesn't have pre-treatment sample
needed=unique(Expansion_group2$pair)
needed=needed[-c(18)]
Expansion_group2=subset(Expansion_group2,pair %in% needed)
#Keep unique clonotype
Expansion_group3=Expansion_group2 %>% drop_na(CTnt)
index_uniqueCT=duplicated(Expansion_group3$uniqueCTnt)
Expansion_group3=Expansion_group3[!index_uniqueCT,]
Expansion_group3$PairCTnt=paste(Expansion_group3$pair,'_',Expansion_group3$CTnt)
#Determine which clones are de_novo
Expansion_group3_Pre=subset(Expansion_group3,treatment=="Pre")
Expansion_group3_Post=subset(Expansion_group3,treatment=="Post")
Expansion_group3_Post$change=''
Expansion_group3_Post$exist=''
for (i in (1:(nrow(Expansion_group3_Post)))){
  if (Expansion_group3_Post[i,"PairCTnt"] %in% Expansion_group3_Pre$PairCTnt){
    Expansion_group3_Post[i,"exist"]="Pre_existed"}
  else{Expansion_group3_Post[i,"exist"]="De_novo"}
}
#Calcaculate the change in frequency
for (i in (1:nrow(Expansion_group3_Post))){
  if (Expansion_group3_Post[i,"exist"]=="Pre_existed"){
    for (j in (1:nrow(Expansion_group3_Pre))){
      if (Expansion_group3_Post[i,"PairCTnt"]==Expansion_group3_Pre[j,"PairCTnt"]){
        Expansion_group3_Post[i,"change"]=Expansion_group3_Post[i,'Frequency']-Expansion_group3_Pre[j,'Frequency']
      }
    }
  }
  else{Expansion_group3_Post[i,"change"]=Expansion_group3_Post[i,"Frequency"]}
}
colnames(Expansion_group3_Post)[9]="Freq.change"
#Define expansion with two criteria
Expansion_group3_Post$Freq.2.group=''
Expansion_group3_Post$Freq.5.group=''
Expansion_group3_Post=Expansion_group3_Post %>% mutate(Freq.2.group=ifelse((Frequency>2 & Freq.change>0),'yes','no'))
Expansion_group3_Post=Expansion_group3_Post %>% mutate(Freq.5.group=ifelse((Frequency>5 & Freq.change>0),'yes','no'))
#Define expanded clone
Expansion_group2$Freq.2.group=''
Expansion_group2$Freq.2.group[Expansion_group2$PairCTnt %in% Expansion_group3_Post$PairCTnt[Expansion_group3_Post$Freq.2.group=='yes']]="Expanded"
Expansion_group2$Freq.2.group[Expansion_group2$Freq.2.group=='']="Non-expanded"
unique(Expansion_group2$Freq.2.group)
table(Expansion_group2$Freq.2.group)
Expansion_group2$Freq.5.group=''
Expansion_group2$Freq.5.group[Expansion_group2$PairCTnt %in% Expansion_group3_Post$PairCTnt[Expansion_group3_Post$Freq.5.group=='yes']]="Expanded"
Expansion_group2$Freq.5.group[Expansion_group2$Freq.5.group=='']="Non-expanded"
Expansion_sample1=as.data.frame(prop.table(table(Expansion_group2$Freq.2.group, Expansion_group2$orig.ident)))
Expansion_sample2=as.data.frame(prop.table(table(Expansion_group2$Freq.5.group, Expansion_group2$orig.ident)))

## Plot the proportion of expanded clones in graphPad

#####Fig3 & S3#####
######Fig. 3B#####
raw.dir = "../data/bcr_data/new_data/"
raw.bcr = dir(raw.dir) %>% .[grep("^RH05",.)]
## load the bcr data that has been merged with single cell bcell set
load("../data/scRepertoire/BCR/combined_83samples.RData")
new.bc.set = readRDS("../data/scRepertoire/rhcc_bcr_mergebcr.rds")
meta = new.bc.set@meta.data
unique.rh = unique(new.bc.set$orig.ident) %>% 
  str_extract(., "RH05") %>% 
  unique()

pass.list = readRDS("../data/bcr_data/shm_file/rhcc_bcr_shm_all.rds")
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
                                                 plot_dir = "../data/BCR_sharing_cjy_chordDiagram/cjy_style_6.0(100_Identical_0)/",
                                                 sample = paste0(rh,"_top10_expanded"),
                                                 chord_transparency = .3,
                                                 link.not.display.threshold = 0,
                                                 niceFacing = T, big.gap = 2,small.gap = 2)

#####Fig4 & S4#####
######Fig. 4B#####
st.bcr.dir = "../data/st_ter_BCR_new/"

bc.set = readRDS("../data/scRepertoire/rhcc_bcr_mergebcr.rds")
bc.meta = bc.set@meta.data
le.samples = c("RH10B", "RH10T", "RH03T", "RH05B", "RH06B", "RH14B","RH15B","RH16B","RH08B","RH17B")
le.group = data.frame(sample = le.samples,
                      group = c("B_cell_responder","B_cell_responder","Non_responder","B_cell_responder",
                                "B_cell_responder","B_cell_responder","Non-responder",
                                "Non-responder","T_cell_responder","T_cell_responder"))
hbv.core.id = readRDS("../data/st_ter_BCR_new/talent_result.rds")
load("../data/st_ter_BCR_new/rh_st_ter_BCR_heavy_chain.RData")
load("../data/st_ter_BCR_new/meta_st_B.RData")
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
myggsave("../data/st_ter_BCR_new/plots/HBV_core_celltype.pdf", width = 6, height = 6, plot = p)

######Fig. 4C#####
dist.df = readRDS("../data/st_ter_BCR_new/tumor_distance_density.rds")
ecdf = ggplot(dist.df, aes(x = distance, color = cell)) + 
  stat_ecdf(linewidth = .7) + labs(x = "Distance", y = "Cumulative fraction") + theme_classic() + 
  mytheme2 + scale_color_manual(values = color)
myggsave(plot = ecdf, filename = "../data/st_ter_BCR_new/plots/dist_ecdf_tumor.pdf",
         width = 7, height = 6)

#####Fig.5 & S5#####
######Fig. 5B######
tumor=readRDS('../data/single_cell_RNA/tumor cell.rds')
all_gene_sets <- msigdbr(species = "Homo sapiens")
all_gene_sets <- as.data.frame(all_gene_sets)
pathway_use2=c("GOMF_COMPLEMENT_COMPONENT_C1Q_COMPLEX_BINDING","GOMF_COMPLEMENT_COMPONENT_C3B_BINDING",
               "GOMF_COMPLEMENT_RECEPTOR_ACTIVITY","HALLMARK_COMPLEMENT","REACTOME_COMPLEMENT_CASCADE",
               "GOBP_COMPLEMENT_ACTIVATION","GOBP_COMPLEMENT_DEPENDENT_CYTOTOXICITY",
               "GOBP_REGULATION_OF_COMPLEMENT_DEPENDENT_CYTOTOXICITY","GOMF_COMPLEMENT_BINDING")        
pathway2 <- list()
for (genelist in pathway_use2) {
  pathway2[[genelist]] <- all_gene_sets$gene_symbol[all_gene_sets$gs_name ==genelist]
}
tumor=AddModuleScore(object = tumor,features = pathway2,name = 'feature')
colnames(tumor@meta.data)[(ncol(tumor@meta.data)-length(pathway2)+1):ncol(tumor@meta.data)]<-names(pathway2)
colnames(tumor@meta.data)
table(tumor$group)
table(tumor$expansion_group)
tumor$group='Non_responder'
tumor$group[tumor$expansion_group=='LowE_Non_progress']='B_responder'
tumor$group[tumor$expansion_group=='HighE_Non_progress']='T_responder'
tumor$group=paste0(tumor$group,"_",tumor$treatment)
table(tumor$group)
tumor$group=factor(tumor$group,c('B_responder_Pre','B_responder_Post','T_responder_Pre','T_responder_Post','Non_responder_Pre','Non_responder_Post'))
data=FetchData(tumor,vars=c('group','treatment','GOMF_COMPLEMENT_BINDING'))
mycolor.s1 = color.s1[c(1:2, 5:6, 3:4)]
p=ggplot(data,aes(group,GOMF_COMPLEMENT_BINDING,fill=group))
p+geom_violin()+geom_boxplot(width = .1, outlier.shape = NA, color = "black", linewidth = violin.box.lwd)+
  theme_classic()+ labs(y = "Complement binding")+ 
  stat_compare_means(label = "p.signif",comparisons = list(c('T_responder_Pre','T_responder_Post'),c('Non_responder_Pre','Non_responder_Post'),
                                                            c('B_responder_Pre','B_responder_Post')),
                     method = "t.test", size = 5, bracket.size = p.lwd)+
  scale_fill_manual(values = mycolor.s1)+RotatedAxis()+
  NoLegend()

######Fig. 5C######
obj=subset(tumor,group%in%c('B_responder_Post','B_responder_Pre'))
pathway_use=c("GOBP_AUTOPHAGIC_CELL_DEATH","GOBP_CELL_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS","GOBP_NECROTIC_CELL_DEATH",
              "GOBP_POSITIVE_REGULATION_OF_CELL_DEATH","GOBP_POSITIVE_REGULATION_OF_CELL_DEATH","GOBP_PROGRAMMED_CELL_DEATH","REACTOME_PROGRAMMED_CELL_DEATH",
              "GOCC_CASPASE_COMPLEX","REACTOME_ACTIVATION_OF_CASPASES_THROUGH_APOPTOSOME_MEDIATED_CLEAVAGE" )
pathway <- list()
for (genelist in pathway_use) {
  pathway[[genelist]] <- all_gene_sets$gene_symbol[all_gene_sets$gs_name ==genelist]
}
obj <- AddModuleScore(object = obj,features = pathway,ctrl = 100, name = 'feature')
colnames(obj@meta.data)[(ncol(obj@meta.data)-length(pathway)+1):ncol(obj@meta.data)]<-names(pathway)
plotdata <- obj@meta.data %>% 
  bind_cols(C3=obj@assays$RNA@data["C3",],PARP1=obj@assays$RNA@data["PARP1",],PARP2=obj@assays$RNA@data["PARP2",],
            PARP3=obj@assays$RNA@data["PARP3",],ANXA1=obj@assays$RNA@data["ANXA1",],HMGB1=obj@assays$RNA@data["HMGB1",],
            CASP3=obj@assays$RNA@data["CASP3",], C5 = obj@assays$RNA@data["C5",])
plotdata3=plotdata[,c(1,29:53)]
plotdata4=plotdata3%>%group_by(orig.ident)%>%summarise_all("mean")
group <- ifelse(str_detect(plotdata4$orig.ident,"A"),"Post","Pre")
plotdata4$group=group
plotdata4=subset(plotdata4,group=="Post")
p1 = ggplot(data=plotdata4,aes(x=C3, y=GOBP_POSITIVE_REGULATION_OF_CELL_DEATH))+geom_point()+
  stat_smooth(method = "lm",se=T)+stat_cor(data = plotdata4,method="spearman")+theme_classic()+ggtitle("Tumor_lenp group")
p2 = ggplot(data=plotdata4,aes(x=C5, y=GOBP_POSITIVE_REGULATION_OF_CELL_DEATH))+geom_point()+
  stat_smooth(method = "lm",se=T)+stat_cor(data = plotdata4,method="spearman")+theme_classic()+ggtitle("Tumor_lenp group")
p1 / p2

#####Fig. S4#####
neo.dir = "../data/pMTnet_output/neoantigen_output/"
neos = dir(neo.dir)
prediction.df = data.frame()
for (res in neos) {
  file.dir = dir(paste0(neo.dir, res))
  df = read.csv(paste0(neo.dir, res, "/", file.dir), header = T, check.names = F)
  df$sample = res
  prediction.df = rbind(prediction.df, df)
}
groups = openxlsx::read.xlsx("../data/sample_single_cell_new_3_groups.xlsx")
prediction.df = data.frame()
for (res in neos) {
  file.dir = dir(paste0(neo.dir, res))
  df = read.csv(paste0(neo.dir, res, "/", file.dir), header = T, check.names = F)
  df$sample = res
  prediction.df = rbind(prediction.df, df)
}
prediction.df = merge(prediction.df, groups, by = "sample")
prediction.df$Rank = as.numeric(prediction.df$Rank)
prediction.df = prediction.df[!is.na(prediction.df$Rank),]
neo.df = prediction.df
neo.df$source = "Neoantigen"
neo.he = neo.df[grep("T|Non", neo.df$group),]
stat.test = neo.he %>% 
  t_test(Rank ~ group, ref.group = "Non-responder") %>% 
  adjust_pvalue() %>% 
  add_xy_position(x = "group", dodge = .9)
stat.test$p.adj.signif = ifelse(stat.test$p.adj > .05, "ns",
                                ifelse(stat.test$p.adj > .01, "*",
                                       ifelse(stat.test$p.adj > .001, "**",
                                              ifelse(stat.test$p.adj > .0001, "***", "****"))))
stat.test$p.signif = ifelse(stat.test$p > .05, "ns",
                            ifelse(stat.test$p > .01, "*",
                                   ifelse(stat.test$p > .001, "**",
                                          ifelse(stat.test$p > .0001, "***", "****"))))
neo.he$group = factor(neo.he$group, levels = c("Non-responder","T_cell_responder"))
print(ggboxplot(data = neo.he, x = "group", y = "Rank", color = "group",
                xlab = "", outlier.shape = NA,ylab = "Rank", main = "Neoantigen pMTnet Prediction", 
                bxp.errorbar = T,size = 0.6, width = 0.6) +
        # guides(color = "none") + 
        stat_compare_means(comparisons = list(c("Non-responder","T_cell_responder")), 
                           label = "p.format",
                           method = "t.test", size = 5) + 
        scale_color_manual(values = two.colors) +
        mytheme2 + theme(legend.position = "none"))

## TAA
taa.dir = "../data/pMTnet_output/TAA_parallel_output/"
taas = dir(taa.dir)
prediction.df = data.frame()
for (res in taas) {
  file.dir = dir(paste0(taa.dir, res))
  df = read.csv(paste0(taa.dir, res, "/", file.dir), header = T, check.names = F)
  df$sample = res
  prediction.df = rbind(prediction.df, df)
}
prediction.df = merge(prediction.df, groups, by = "sample")
prediction.df$Rank = as.numeric(prediction.df$Rank)
prediction.df = prediction.df[!is.na(prediction.df$Rank),]
compare.list = list()
for (i in c(1:ncol(combn(unique(prediction.df$group),2)))){
  compare.list[[i]] = combn(unique(prediction.df$group),2)[,i]
}
taa.df = prediction.df
taa.he = taa.df[grep("T|Non", taa.df$group),]
stat.test = taa.he %>% 
  t_test(Rank ~ group, ref.group = "Non-responder") %>% 
  adjust_pvalue() %>% 
  add_xy_position(x = "group", dodge = .9)
stat.test$p.adj.signif = ifelse(stat.test$p.adj > .05, "ns",
                                ifelse(stat.test$p.adj > .01, "*",
                                       ifelse(stat.test$p.adj > .001, "**",
                                              ifelse(stat.test$p.adj > .0001, "***", "****"))))
stat.test$p.signif = ifelse(stat.test$p > .05, "ns",
                            ifelse(stat.test$p > .01, "*",
                                   ifelse(stat.test$p > .001, "**",
                                          ifelse(stat.test$p > .0001, "***", "****"))))
taa.he$group = factor(taa.he$group, levels = c("Non-responder","T_cell_responder"))
print(ggboxplot(data = taa.he, x = "group", y = "Rank", color = "group",
          xlab = "", outlier.shape = NA,ylab = "Rank", main = "TAA pMTnet Prediction", 
          bxp.errorbar = T,size = 0.6, width = 0.6) +
  # guides(color = "none") + 
  stat_compare_means(comparisons = list(c("Non-responder","T_cell_responder")), 
                     label = "p.format",
                     method = "t.test", size = 5) + 
  scale_color_manual(values = two.colors) +
  mytheme2 + theme(legend.position = "none"))

