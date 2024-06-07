syEnrich <- function(degs,
                     outpath="tmp"){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(xlsx)
  
  allcluster_go=data.frame()
  allcluster_kegg=data.frame()
  
  for (i in unique(degs$cluster)) {
    small_degs=filter(degs,degs$cluster==i)
    gene_name_df=bitr(small_degs$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
    
    #GO富集
    go <- enrichGO(gene         = unique(gene_name_df$ENTREZID),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)
    go=clusterProfiler::simplify(go) #给term去冗余，非必要步骤
    go_res=go@result
    if (dim(go_res)[1] != 0) {
      go_res$cluster=i
      allcluster_go=rbind(allcluster_go,go_res)
    }
    
    #KEGG富集
    kegg <- enrichKEGG(unique(gene_name_df$ENTREZID), organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                       minGSSize = 5,maxGSSize = 1000,qvalueCutoff = 0.2,use_internal_data = FALSE)
    kegg_res=kegg@result
    if (dim(kegg_res)[1] != 0) {
      symbol_v=c()
      for (k in 1:dim(kegg_res)[1]) {
        symbol_v=append(symbol_v,
                        paste(gene_name_df$SYMBOL[gene_name_df$ENTREZID %in% strsplit(kegg_res[k,"geneID"],"/")[[1]]],collapse = ", ")
        )
      }
      kegg_res$cluster=i
      kegg_res$symbol=symbol_v
      allcluster_kegg=rbind(allcluster_kegg,kegg_res)
    }
  }
  
  write.xlsx(allcluster_go,  file = paste(outpath,".GO.xls",sep = ""),row.names = F,col.names = T)
  write.xlsx(allcluster_kegg,file = paste(outpath,".KEGG.xls",sep = ""),   row.names = F,col.names = T)
}
