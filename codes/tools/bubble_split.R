bubble_split=function(
    seu.obj=NULL,
    celltype.column=NULL,
    celltype.order=NULL,
    group.column=NULL,
    deg=NULL,
    rect_color="black",
    rect_size=0.5,
    color1=brewer.pal(9, "Reds"),
    color2=brewer.pal(9, "Blues"),
    color3=brewer.pal(9, "Greens"),
    color4=brewer.pal(9, "Purples")
){
  library(tidyverse)
  library(RColorBrewer)
  library(scales)
  library(reshape2)
  library(ggnewscale)
  
  
  ct.index=which(colnames(seu.obj@meta.data) == celltype.column)
  colnames(seu.obj@meta.data)[ct.index] = "celltype"
  g.index=which(colnames(seu.obj@meta.data) == group.column)
  colnames(seu.obj@meta.data)[g.index] = "condition"
  seu.obj@meta.data$condition=as.character(seu.obj@meta.data$condition)
  group.level.n=seu.obj@meta.data$condition %>% unique() %>% length()
  group.level=seu.obj@meta.data$condition %>% unique() %>% sort()
  seu.obj@meta.data$celltype = factor(
    as.character(seu.obj@meta.data$celltype),
    levels = celltype.order
  )
  Idents(seu.obj) = "celltype"
  
  
  deg$cluster=factor(as.character(deg$cluster),levels = celltype.order)
  deg=deg%>%arrange(cluster,gene)
  
  
  seu.obj@meta.data$celltype_condition=paste(
    as.character(seu.obj@meta.data$celltype),
    seu.obj@meta.data$condition,
    sep = "_"
  )
  seu.obj@meta.data$celltype_condition=factor(
    seu.obj@meta.data$celltype_condition,
    levels = paste(rep(celltype.order,each=group.level.n),group.level,sep = "_")
  )
  
  
  Idents(seu.obj)="celltype_condition"
  bubble.exp=AverageExpression(seu.obj,assays = "RNA",features = deg$gene)
  bubble.exp=scale(t(bubble.exp$RNA))
  bubble.exp=as.data.frame(bubble.exp)
  bubble.exp$type=rownames(bubble.exp)
  bubble.exp=bubble.exp%>%reshape2::melt(id.vars=c("type"))
  colnames(bubble.exp)[2:3]=c("gene","average.exp.scaled")
  bubble.exp$index=paste(bubble.exp$type,bubble.exp$gene,sep = "_")
  
  
  percent.df=data.frame()
  for (typei in sort(unique(as.character(bubble.exp$type)))) {
    tmpCB=seu.obj@meta.data[seu.obj@meta.data$celltype_condition == typei,"CB"]
    gene.exp=seu.obj[["RNA"]]@counts[sort(unique(as.character(bubble.exp$gene))),tmpCB]
    gene.exp=as.matrix(gene.exp)
    percent=rowSums(gene.exp > 0) / ncol(gene.exp)
    percent=as.data.frame(percent)
    percent$index=paste(typei,rownames(percent),sep = "_")
    rownames(percent)=NULL
    
    percent.df=percent.df%>%rbind(percent)
  }
  
  
  bubble.exp=bubble.exp%>%inner_join(percent.df,by = "index")
  celltype.pattern=paste0("(",celltype.order%>%paste(collapse = "|"),")_")
  condition.pattern=paste0("_(",group.level %>% paste(collapse = "|"),")")
  bubble.exp$celltype=str_replace(bubble.exp$type,condition.pattern,"")
  bubble.exp$condition=str_replace(bubble.exp$type,celltype.pattern,"")
  
  ###补齐数据框##########################################
  expected.lines=length(celltype.order)*group.level.n*length(deg$gene)
  if(dim(bubble.exp)[1] < expected.lines){
    remain.df.index=c()
    remain.df.type=c()
    remain.df.gene=c()
    remain.df.celltype=c()
    remain.df.condition=c()
    for (i in sort(unique(as.character(bubble.exp$celltype)))) {
      for (j in sort(unique(as.character(bubble.exp$condition)))) {
        for (k in sort(unique(as.character(bubble.exp$gene)))) {
          one.index=paste(i,j,k,sep = "_")
          if (one.index %in% bubble.exp$index) {
          }else{
            remain.df.index=append(remain.df.index,one.index)
            remain.df.type=append(remain.df.type,paste(i,j,sep = "_"))
            remain.df.gene=append(remain.df.gene,k)
            remain.df.celltype=append(remain.df.celltype,i)
            remain.df.condition=append(remain.df.condition,j)
          }
        }
      }
    }
    remain.df=data.frame(type=remain.df.type,
                         gene=remain.df.gene,
                         average.exp.scaled=NA,
                         index=remain.df.index,
                         percent=NA,
                         celltype=remain.df.celltype,
                         condition=remain.df.condition)
    bubble.exp=bubble.exp%>%rbind(remain.df)
  }
  
  ###调整横纵轴########################################
  ### 本系列代码由公众号【TOP生物信息】原创
  bubble.exp$final.celltype=""
  bubble.exp$final.gene=""
  if(group.level.n == 2){
    for (linei in 1:dim(bubble.exp)[1]) {
      linegroup=bubble.exp[linei,"condition"]
      if(linegroup == group.level[1]){
        bubble.exp[linei,"final.celltype"]=paste0(bubble.exp[linei,"celltype"],"_1")
        bubble.exp[linei,"final.gene"]=paste0(bubble.exp[linei,"gene"],"_1")
      }
      if(linegroup == group.level[2]){
        bubble.exp[linei,"final.celltype"]=paste0(bubble.exp[linei,"celltype"],"_1")
        bubble.exp[linei,"final.gene"]=paste0(bubble.exp[linei,"gene"],"_2")
      }
    }
  }
  
  if(group.level.n == 3){
    for (linei in 1:dim(bubble.exp)[1]) {
      linegroup=bubble.exp[linei,"condition"]
      if(linegroup == group.level[1]){
        bubble.exp[linei,"final.celltype"]=paste0(bubble.exp[linei,"celltype"],"_1")
        bubble.exp[linei,"final.gene"]=paste0(bubble.exp[linei,"gene"],"_1")
      }
      if(linegroup == group.level[2]){
        bubble.exp[linei,"final.celltype"]=paste0(bubble.exp[linei,"celltype"],"_1")
        bubble.exp[linei,"final.gene"]=paste0(bubble.exp[linei,"gene"],"_2")
      }
      if(linegroup == group.level[3]){
        bubble.exp[linei,"final.celltype"]=paste0(bubble.exp[linei,"celltype"],"_1")
        bubble.exp[linei,"final.gene"]=paste0(bubble.exp[linei,"gene"],"_3")
      }
    }
  }
  
  if(group.level.n == 4){
    for (linei in 1:dim(bubble.exp)[1]) {
      linegroup=bubble.exp[linei,"condition"]
      if(linegroup == group.level[1]){
        bubble.exp[linei,"final.celltype"]=paste0(bubble.exp[linei,"celltype"],"_2")
        bubble.exp[linei,"final.gene"]=paste0(bubble.exp[linei,"gene"],"_1")
      }
      if(linegroup == group.level[2]){
        bubble.exp[linei,"final.celltype"]=paste0(bubble.exp[linei,"celltype"],"_2")
        bubble.exp[linei,"final.gene"]=paste0(bubble.exp[linei,"gene"],"_2")
      }
      if(linegroup == group.level[3]){
        bubble.exp[linei,"final.celltype"]=paste0(bubble.exp[linei,"celltype"],"_1")
        bubble.exp[linei,"final.gene"]=paste0(bubble.exp[linei,"gene"],"_1")
      }
      if(linegroup == group.level[4]){
        bubble.exp[linei,"final.celltype"]=paste0(bubble.exp[linei,"celltype"],"_1")
        bubble.exp[linei,"final.gene"]=paste0(bubble.exp[linei,"gene"],"_2")
      }
    }
  }
  
  ###
  if(group.level.n == 4){
    rep.n=2
    paste.v=1:2
  }else{
    rep.n=1
    paste.v=1
  }
  bubble.exp$final.celltype=factor(
    bubble.exp$final.celltype,
    levels = paste(rep(celltype.order,each=rep.n),paste.v,sep = "_")
  )
  
  if(group.level.n == 3){
    rep.n=3
    paste.v=1:3
  }else{
    rep.n=2
    paste.v=1:2
  }
  bubble.exp$final.gene=factor(
    bubble.exp$final.gene,
    levels = paste(rep(deg$gene,each=rep.n),paste.v,sep = "_")
  )
  
  
  ###
  bubble.exp=bubble.exp%>%arrange(final.gene,final.celltype)
  minmax=range(bubble.exp$average.exp.scaled[!is.na(bubble.exp$average.exp.scaled)])
  
  if(group.level.n==2 | group.level.n==3){
    cum.sum=cumsum(deg$cluster %>% table() * group.level.n) + 0.5
    celltype.n=length(celltype.order)
    
    rect.df=data.frame(
      xmin=c(0.5,cum.sum[-celltype.n]),
      xmax=cum.sum,
      ymin=seq(0.5,celltype.n-0.5,1),
      ymax=seq(1.5,celltype.n+0.5,1)
    )
  }
  if(group.level.n==4){
    cum.sum=cumsum(deg$cluster %>% table() * 2) + 0.5
    celltype.n=length(celltype.order)
    
    rect.df=data.frame(
      xmin=c(0.5,cum.sum[-celltype.n]),
      xmax=cum.sum,
      ymin=c(0.5,c(1:(celltype.n-1))*2+0.5),
      ymax=c(1:celltype.n)*2+0.5
    )
  }
  
  if(group.level.n == 2){
    breaks.x=c(1:length(deg$gene))*2-0.5
    breaks.y=1:length(celltype.order)
  } else if(group.level.n == 3){
    breaks.x=c(1:length(deg$gene))*3-1
    breaks.y=1:length(celltype.order)
  } else {
    breaks.x=c(1:length(deg$gene))*2-0.5
    breaks.y=c(1:length(celltype.order))*2 - 0.5
  }
  
  p1=ggplot(data = NULL,aes(x=as.numeric(final.gene),y=as.numeric(final.celltype),size=percent))+
    geom_point(data = bubble.exp[bubble.exp$condition == group.level[1],],mapping = aes(color=average.exp.scaled))+
    scale_color_gradientn(group.level[1],colors = color1,limits=minmax)+
    new_scale_color()+
    geom_point(data = bubble.exp[bubble.exp$condition == group.level[2],],mapping = aes(color=average.exp.scaled))+
    scale_color_gradientn(group.level[2],colors = color2,limits=minmax)
  if(group.level.n==3){
    p1=p1+
      new_scale_color()+
      geom_point(data = bubble.exp[bubble.exp$condition == group.level[3],],mapping = aes(color=average.exp.scaled))+
      scale_color_gradientn(group.level[3],colors = color3,limits=minmax)
  }
  if(group.level.n==4){
    p1=p1+
      new_scale_color()+
      geom_point(data = bubble.exp[bubble.exp$condition == group.level[3],],mapping = aes(color=average.exp.scaled))+
      scale_color_gradientn(group.level[3],colors = color3,limits=minmax)+
      new_scale_color()+
      geom_point(data = bubble.exp[bubble.exp$condition == group.level[4],],mapping = aes(color=average.exp.scaled))+
      scale_color_gradientn(group.level[4],colors = color4,limits=minmax)
  }
  p1+annotate("rect",xmin = rect.df$xmin,xmax = rect.df$xmax,ymin=rect.df$ymin,ymax=rect.df$ymax,alpha=0,color=rect_color,size=rect_size)+
    scale_x_continuous(breaks=breaks.x,labels=deg$gene,expand = c(0,0))+
    scale_y_continuous(breaks=breaks.y,labels=celltype.order,expand = c(0,0))+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_line(color = "grey"),
      axis.title = element_blank(),
      axis.ticks.length = unit(0.15,"cm"),
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,color="black",size=14),
      axis.text.y.left = element_text(color="black",size=14),
      legend.position = "right",
      legend.direction = "horizontal"
    )
}
