kegg_loop=function(
    enrich.res,
    filename = "kegg_circos"
){
  library(circlize)
  ### 继续整理数据 #############################################
  enrich.res$ID=paste0("ko",enrich.res$ID)
  enrich.res$left=0
  enrich.res$this_pathway_gene=enrich.res$BgRatio %>% sapply(function(x){str_split(x,"/")[[1]][1]}) %>% as.numeric()
  enrich.res$all_pathway_gene=enrich.res$BgRatio %>% sapply(function(x){str_split(x,"/")[[1]][2]}) %>% as.numeric()
  enrich.res$DEGnum=enrich.res$GeneRatio %>% sapply(function(x){str_split(x,"/")[[1]][2]}) %>% as.numeric()
  enrich.res$right=max(enrich.res$this_pathway_gene)
  enrich.res$RichFactor=enrich.res$Count / enrich.res$this_pathway_gene
  enrich.res$big.annotion=factor(enrich.res$big.annotion,levels = sort(unique(enrich.res$big.annotion)))
  enrich.res=enrich.res%>%arrange(big.annotion,desc(RichFactor))
  rownames(enrich.res)=enrich.res$ID
  
  ### 开始画图 #################################################
  pdf(paste0(filename,".pdf"),width = 10,height = 10)
  plotdata=enrich.res
  circos.par(
    clock.wise=TRUE, #The direction for adding sectors. Default is TRUE.
    start.degree = 90, #if it is set to 90, sectors start from the top center of the circle
    gap.degree = 0.8, #Gap between two neighbour sectors
    xaxis.clock.wise = TRUE #The direction in the x-axes for all sectors. Default is TRUE.
  )
  
  ### 第一圈：分类信息 #####################################################################
  library(RColorBrewer)
  library(scales)
  color7=c('#6b97e4','#f7cd0b','#fc837e',"#8dd3c7",
            "#9b79e2","#009c77","#f781bf")
  family.len=length(unique(plotdata$big.annotion))
  color.selected = color7[1:family.len]
  color.order=plyr::mapvalues(x = as.character(plotdata$big.annotion), from = as.character(unique(plotdata$big.annotion)), to = color.selected)
  
  circos.genomicInitialize(plotdata[,c("ID","left","right","big.annotion")], plotType = NULL)  #Initialize circular plot with any genomic data
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = color.order,
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter') #get.cell.meta.data()函数获取扇区的元数据，查看可以获取什么信息?get.cell.meta.data
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')
      circos.text(xlim, ylim, sector.name, cex = 1.2, niceFacing = T)  #cex: Font size; niceFacing: Should the facing of text be adjusted to fit human eyes?
    } )
  
  ### 第二圈：这个通路总共有多少基因 ##################################
  circos.genomicTrackPlotRegion(
    plotdata[,c("ID", "left", "this_pathway_gene")], track.height = 0.05, bg.border = NA, 
    stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "#9b79e2", border = NA, ...)
      ylim = get.cell.meta.data('ycenter')
      xlim = plotdata[,c("ID", "left", "this_pathway_gene")][get.cell.meta.data('sector.index'),3] / 2
      sector.name = plotdata[,c("ID", "left", "this_pathway_gene")][get.cell.meta.data('sector.index'),3]
      circos.text(xlim, ylim, sector.name, cex = 0.7, niceFacing = T)
    } )
  ### 第三圈：generatio相关 #########################################
  plotdata.3 <- plotdata[,c("ID","left","Count","DEGnum","right")]
  plotdata.3$ratio=plotdata.3$Count / plotdata.3$DEGnum
  plotdata.3$len=plotdata.3$ratio * plotdata.3$right
  plotdata.3$len2=plotdata.3$right - plotdata.3$len
  
  tmpdf1=plotdata.3[,c("ID","left","len")]
  colnames(tmpdf1)=c("ID","start","end")
  tmpdf1$type=1
  tmpdf2=plotdata.3[,c("ID","len","right")]
  colnames(tmpdf2)=c("ID","start","end")
  tmpdf2$type=2
  tmpdata=tmpdf1%>%rbind(tmpdf2)
  color_assign <- colorRamp2(breaks = c(1, 2), col = c('#009d78', '#9dbf57'))
  
  circos.genomicTrackPlotRegion(
    tmpdata, track.height = 0.05, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value), border = NA, ...)
      
      ylim = get.cell.meta.data('ycenter')
      
      xlim = plotdata.3[get.cell.meta.data('sector.index'),"len"] / 2
      sector.name = plotdata.3[get.cell.meta.data('sector.index'),"Count"]
      circos.text(xlim, ylim, sector.name, cex = 0.7, niceFacing = T)
      
      xlim = plotdata.3[get.cell.meta.data('sector.index'),"len"] + plotdata.3[get.cell.meta.data('sector.index'),"len2"] / 2
      sector.name = plotdata.3[get.cell.meta.data('sector.index'),"DEGnum"] -  plotdata.3[get.cell.meta.data('sector.index'),"Count"]
      circos.text(xlim, ylim, sector.name, cex = 0.7, niceFacing = T)
    } )
  ### 第四圈：p值 ##################################################
  plotdata.4 <- plotdata[,c("ID", "left","right", "p.adjust")]
  total.len=unique(plotdata.4$right)
  plotdata.4$p.adjust_log10_neg = -log10(plotdata.4$p.adjust)
  plotdata.4$relative_value = plotdata.4$p.adjust_log10_neg / max(plotdata.4$p.adjust_log10_neg) * total.len
  
  p_max <- max(plotdata.4$p.adjust_log10_neg) %>% ceiling()  #定义一个 p 值的极值，以方便后续作图
  colorsChoice <- colorRampPalette(c('#FF906F', '#861D30'))  #这两句用于定义 p 值的渐变颜色
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  
  plotdata.4=plotdata.4[,c("ID", "left", "relative_value", "p.adjust_log10_neg")]
  
  circos.genomicTrackPlotRegion(
    plotdata.4, track.height = 0.05, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value), border = NA, ...)
      
      ylim = get.cell.meta.data('ycenter')
      xlim = plotdata.4[get.cell.meta.data('sector.index'),"relative_value"] / 2
      sector.name = plotdata.4[get.cell.meta.data('sector.index'),"p.adjust_log10_neg"] %>% round(2)
      
      tmpx = -log10(0.01) / max(plotdata.4$p.adjust_log10_neg) * total.len
      circos.lines(c(tmpx, tmpx), c(get.cell.meta.data('ylim')[1], get.cell.meta.data('ylim')[2]), col = 'gray', lwd = 0.45)
      
      circos.text(xlim, ylim, sector.name, cex = 0.7, niceFacing = T)
    } )
  ### 第五圈：富集因子 ##############################################
  plotdata.5 <- plotdata[,c("ID", "left", "right", "RichFactor","big.annotion")]
  color_assign = color.selected
  names(color_assign)=as.character(unique(plotdata.5$big.annotion))
  
  circos.genomicTrack(
    plotdata.5, ylim = c(0, max(plotdata.5$RichFactor)), track.height = 0.4, bg.col = 'gray95', bg.border = NA,
    panel.fun = function(region, value, ...) {
      sector.name = get.cell.meta.data('sector.index')
      circos.genomicRect(region, value, col = color_assign[plotdata.5[sector.name,5]], border = NA, ytop.column = 1, ybottom = 0, ...)
    } )
  ###
  circos.clear()
  
  ### 绘制图例 ################################
  library(ComplexHeatmap)
  
  category_legend <- Legend(
    labels = as.character(unique(plotdata$big.annotion)),
    background = color.selected, 
    type = 'points', pch = NA, 
    labels_gp = gpar(fontsize = 8), 
    grid_height = unit(0.5, 'cm'), 
    grid_width = unit(0.5, 'cm')
  )
  
  thispathway_gene_legend <- Legend(
    labels = c('gene number of this pathway'), 
    background = c('#9b79e2'), 
    type = 'points', pch = NA, 
    labels_gp = gpar(fontsize = 8), 
    grid_height = unit(0.5, 'cm'), 
    grid_width = unit(0.5, 'cm')
  )
  
  generatio_legend <- Legend(
    labels = c('DEGs in this pathway',"other DEGs"), 
    background = c("#009d78","#9dbf57"), 
    type = 'points', pch = NA, 
    labels_gp = gpar(fontsize = 8), 
    grid_height = unit(0.5, 'cm'), 
    grid_width = unit(0.5, 'cm')
  )
  
  pvalue_legend <- Legend(
    col_fun = colorRamp2(
      breaks = 0:p_max, 
      col = colorRampPalette(c('#FF906F', '#861D30'))(p_max + 1)
    ),
    legend_height = unit(3, 'cm'), 
    labels_gp = gpar(fontsize = 8), 
    title = '-log10(p.adjust)',
    title_gp = gpar(fontsize = 9), 
    title_position = 'lefttop',
    direction = "horizontal"
  )
  
  pack_Legend <- packLegend(
    category_legend, 
    thispathway_gene_legend, 
    generatio_legend,
    pvalue_legend,
    row_gap = unit(0.2, "cm")
  )
  
  grid.draw(pack_Legend)
  dev.off()
  
  return("ok!")
}
