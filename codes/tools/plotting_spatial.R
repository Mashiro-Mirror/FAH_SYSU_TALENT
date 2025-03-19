# ---------------------------------------------------------------------------- #
#  
#  plotting_spatial.R
# 
# ---------------------------------------------------------------------------- #
#  
#  AUTHOR: WU wenrui (wuwenruiwwr@outlook.com)
#  
#  DATE: 2024-04-25 (v0.1)
#  DATE: 2024-05-14 (v0.2) [Feature] Make sr_add_feature() available when given
#                          either single or multiple features. 
#  DATE: 2024-05-04 (v0.3) Update *_point/circle() functions to automatically 
#                          determine colors or orders if not given. 
#  DATE: 2024-05-17 (v0.4) Update *_point/circle() functions to use a data.frame 
#                          as input parameter instead of a Seurat object. 
#  DATE: 2024-05-17 (v0.5) [Feature] Add functions to draw contour of regions. 
#  
#  FUNCTION
#  (1) get_celltype_color()
#      .get_celltype_color()
#  (2) st_plot_discrete_point()
#      st_plot_discrete_circle()
#      .st_plot_discrete_data()
#      .st_plot_discrete_point()
#      .st_plot_discrete_circle()
#  (3) st_plot_continuous_point()
#      st_plot_continuous_circle()
#      .st_plot_continuous_data()
#      .st_plot_continuous_point()
#      .st_plot_continuous_circle()
#  (4) sr_add_feature()
#  (5) st_plot_size()
#  (6) st_plot_contour_1_data()
#      st_plot_contour_1()
#      .st_plot_contour_data_arc()
#      .st_plot_contour_data_segment()
#      .st_plot_contour_data()
#  
#  NOTICE
#  (1) 使用前直接source()这个R文件，示例放在恒为FALSE的IF语句中，不影响函数导入。
#  (2) 标注了"@export"或不以"."开头的函数为推荐使用的主要函数；未标注"@export"或
#      以"."开头的其他函数为中间函数，当主要函数无法数据需求时可自行使用中间函数
#      组装新函数。
#  
# ---------------------------------------------------------------------------- #


# Some Colors #####

#' Colors for Celltype
#' 
#' 常用的细胞大类注释颜色。
#' 
COLOR_CELLTYPE <- c(
  "#ce1256" = "T",  
  "#ce1256" = "T/NK", 
  "#ce1256" = "T&NK", 
  "#238b45" = "B", 
  "#238b45" = "B/Plasma", 
  "#fdbf6f" = "Plasma", 
  "#6a3d9a" = "Myeloid", 
  "#F7CCC1" = "Neutrophil", 
  "#A1CFFA" = "pDC", 
  "#6F89A2" = "Mast", 
  "#df65b0" = "Fibroblast", 
  "#1f78b4" = "Endothelial", 
  "#FDD84C" = "Epithelial", 
  "#C9BEE8" = "Tumor", 
  "#C9BEE8" = "Hep/Tumor", 
  "#eaeaea" = "Hepatocyte", 
  "#eaeaea" = "Tumor&Hepatocyte", 
  "#eaeaea" = "Doublet"
)

#' Get Color for Celltype
#'
#' 根据提前设置好的常量`COLOR_CELLTYPE`，输出与输入细胞类型匹配的16进制颜色代码。
#' 
#' @param x 长度为1的细胞类型向量。
#'
#' @return 与输入细胞类型匹配的颜色代码。
#'
.get_celltype_color <- function(x) {
  color <- COLOR_CELLTYPE %>%
    purrr::map(~ stringr::str_subset(x, .x)) %>%
    unlist() %>%
    names()
  if (is.null(color)) {color <- NA}
  return(color[1])
}

#' Get Colors for Celltype
#'
#' 根据提前设置好的常量`COLOR_CELLTYPE`，输出与输入细胞类型匹配的16进制颜色代码。
#' 对没有与其匹配颜色代码的细胞类型给予提醒。
#'
#' @param x 细胞类型向量。
#'
#' @return 与输入细胞类型匹配的颜色代码。
#' @export
#'
get_celltype_color <- function(x) {
  color <- x %>%
    purrr::map(.get_celltype_color) %>%
    unlist()
  n_na <- sum(is.na(color))
  if (n_na) {
    x_na <- stringr::str_c(unique(x[is.na(color)]), collapse = ", ")
    warning(stringr::str_glue("Not matched: {x_na}. "), call. = FALSE)
  }
  return(color)
}


# Discrete Data #####

#' Spatial Plotting: Data Preprocess for Discrete Data
#'
#' @param dat_p 数据框，至少包含三列：绘制点的x轴坐标、y轴坐标、用于映射绘制颜色
#' 的分类变量。
#' @param x 字符串，`dat_p`中存放点x轴坐标的列名。
#' @param y 字符串，`dat_p`中存放点y轴坐标的列名。
#' @param color.by 字符串，`dat_p`中存放需要以不同颜色绘制的分类变量的列名。
#' @param order 向量，用于控制颜色legend中分类变量的展示顺序。
#' @param color 与`order`等长的向量，分类变量不同值各自对应的颜色。
#' @param y.rev 逻辑值，用于控制是否将y轴坐标顺序颠倒（在ggplot2中，较小的y值位
#' 于绘制图片的下方；而在空转数据meta.data中，较小的y值位于病理图片的上方）。
#'
#' @return 数据框，包含可用于绘图的规整数据。
#'
.st_plot_discrete_data <- function(dat_p, x = "col", y = "row", 
                                   color.by, order = NULL, color = NULL, 
                                   y.rev = TRUE) {
  dat_p <- dat_p %>%
    dplyr::select(dplyr::all_of(c("x" = x, "y" = y, "color.by" = color.by))) 
  if (is.null(order)) {order <- sort(unique(dat_p$color.by))}
  if (is.null(color)) {color <- color_Reds_Blues(length(order))}
  dat_p <- dat_p %>%
    dplyr::mutate(
      dplyr::across(c(x, y), as.numeric), 
      dplyr::across(color.by, ~ factor(.x, levels = order)), 
      color = setNames(color, order)[color.by]
    )
  if (y.rev) {dat_p$y <- -1 * dat_p$y}
  return(dat_p)
}

#' Spatial Plotting: Point Plot for Discrete Data
#' 
#' @param dat_p 数据框，经过`.st_plot_discrete_data()`函数预处理的规整数据。
#'
#' @return 使用不同颜色绘制不同分类变量的点图。
#'
.st_plot_discrete_point <- function(dat_p) {
  color_values <- dplyr::distinct(dat_p, color.by, color) %>% na.omit()
  color_values <- setNames(color_values$color, color_values$color.by)
  p <- dat_p %>% 
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = color.by)) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::labs(color = NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2)))
  return(p)
}

#' Spatial Plotting: Circle Plot for Discrete Data
#'
#' @param dat_p 数据框，经过`.st_plot_discrete_data()`函数预处理的规整数据。
#'
#' @return 使用不同颜色绘制不同分类变量的相切圆图。
#'
.st_plot_discrete_circle <- function(dat_p) {
  color_values <- dplyr::distinct(dat_p, color.by, color) %>% na.omit()
  color_values <- setNames(color_values$color, color_values$color.by)
  p <- dat_p %>% 
    ggplot2::ggplot() +
    ggforce::geom_circle(
      ggplot2::aes(x0 = x, y0 = y, fill = color.by, r = 0.5), 
      color = NA
    ) +
    ggplot2::scale_fill_manual(values = color_values) +
    ggplot2::labs(fill = NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 2)))
  return(p)
}

#' Spatial Plotting: Point Plot for Discrete Data
#'
#' 绘制速度较快，用于预览参数设置是否正确。
#'
#' @param dat_p 数据框，至少包含三列：绘制点的x轴坐标、y轴坐标、用于映射绘制颜色
#' 的分类变量。
#' @param x 字符串，`dat_p`中存放点x轴坐标的列名。
#' @param y 字符串，`dat_p`中存放点y轴坐标的列名。
#' @param color.by 字符串，`dat_p`中存放需要以不同颜色绘制的分类变量的列名。
#' @param order 向量，用于控制颜色legend中分类变量的展示顺序。
#' @param color 与`order`等长的向量，分类变量不同值各自对应的颜色。
#' @param y.rev 逻辑值，用于控制是否将y轴坐标顺序颠倒（在ggplot2中，较小的y值位
#' 于绘制图片的下方；而在空转数据meta.data中，较小的y值位于病理图片的上方）。
#'
#' @return 使用不同颜色绘制不同分类变量的点图。
#' @export
#'
st_plot_discrete_point <- function(dat_p, x = "col", y = "row", 
                                   color.by, order = NULL, color = NULL, 
                                   y.rev = TRUE) {
  dat_p %>%
    .st_plot_discrete_data(x, y, color.by, order, color, y.rev) %>%
    .st_plot_discrete_point()
}

#' Spatial Plotting: Circle Plot for Discrete Data
#'
#' 绘制速度较慢，使用点图预览确认参数设置无误后，用于正式图片的绘制。
#'
#' @param dat_p 数据框，至少包含三列：绘制点的x轴坐标、y轴坐标、用于映射绘制颜色
#' 的分类变量。
#' @param x 字符串，`dat_p`中存放点x轴坐标的列名。
#' @param y 字符串，`dat_p`中存放点y轴坐标的列名。
#' @param color.by 字符串，`dat_p`中存放需要以不同颜色绘制的分类变量的列名。
#' @param order 向量，用于控制颜色legend中分类变量的展示顺序。
#' @param color 与`order`等长的向量，分类变量不同值各自对应的颜色。
#' @param y.rev 逻辑值，用于控制是否将y轴坐标顺序颠倒（在ggplot2中，较小的y值位
#' 于绘制图片的下方；而在空转数据meta.data中，较小的y值位于病理图片的上方）。
#' 
#' @return 使用不同颜色绘制不同分类变量的相切圆图。
#' @export
#' 
st_plot_discrete_circle <- function(dat_p, x = "col", y = "row", 
                                    color.by, order = NULL, color = NULL, 
                                    y.rev = TRUE) {
  dat_p %>%
    .st_plot_discrete_data(x, y, color.by, order, color, y.rev) %>%
    .st_plot_discrete_circle()
}


# Demo

if (FALSE) {
  sr_2938 <- read_rds("/data01/home/wuwenrui/CRLM_20240320/data_input/ST2938B_spotlight_Spatial_TLS.rds")
  
  color.by <- "predict_CellType"
  order <- c(
    "T&NK cell", 
    "B/Plasma cell", 
    "Myeloid cell", 
    "Fibroblast", 
    "Endothelial cell", 
    "Epithelial cell", 
    "Hepatocyte"
  )
  color <- order %>% get_celltype_color()
  x <- "col"
  y <- "row"
  
  # To preview (fast)
  sr_2938@meta.data %>%
    st_plot_discrete_point(x, y, color.by, order, color)
  
  # To save (slow)
  p <- sr_2938@meta.data %>%
    st_plot_discrete_circle(x, y, color.by, order, color)
}


# Continuous Data #####

#' Spatial Plotting: Data Preprocess for Continuous Data
#'
#' @param dat_p 数据框，至少包含三列：绘制点的x轴坐标、y轴坐标、用于映射绘制颜色
#' 的连续变量。
#' @param x 字符串，`dat_p`中存放点x轴坐标的列名。
#' @param y 字符串，`dat_p`中存放点y轴坐标的列名。
#' @param color.by 字符串，`dat_p`中存放需要以不同颜色绘制的连续变量的列名。
#' @param y.rev 逻辑值，用于控制是否将y轴坐标顺序颠倒（在ggplot2中，较小的y值位
#' 于绘制图片的下方；而在空转数据meta.data中，较小的y值位于病理图片的上方）。
#'
#' @return 数据框，包含可用于绘图的规整数据。
#'
.st_plot_continuous_data <- function(dat_p, x = "col", y = "row", 
                                     color.by, y.rev = TRUE) {
  dat_p <- dat_p %>%
    dplyr::select(dplyr::all_of(c("x" = x, "y" = y, "color.by" = color.by))) %>%
    dplyr::mutate(
      dplyr::across(c(x, y, color.by), as.numeric)
    )
  if (y.rev) {dat_p$y <- -1 * dat_p$y}
  return(dat_p)
}

#' Spatial Plotting: Point Plot for Continuous Data
#' 
#' @param dat_p 数据框，经过`.st_plot_continuous_data()`函数预处理的规整数据。
#' @param color 向量，映射连续变量梯度值的若干颜色。
#'
#' @return 使用不同颜色绘制不同分类变量的点图。
#'
.st_plot_continuous_point <- function(dat_p, color = color_Spectral()) {
  p <- dat_p %>% 
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = color.by)) +
    ggplot2::scale_color_gradientn(colors = color) +
    ggplot2::labs(color = NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme_void()
  return(p)
}

#' Spatial Plotting: Circle Plot for Continuous Data
#'
#' @param dat_p 数据框，经过`.st_plot_continuous_data()`函数预处理的规整数据。
#' @param color 向量，映射连续变量梯度值的若干颜色。
#'
#' @return 使用不同颜色绘制不同分类变量的相切圆图。
#'
.st_plot_continuous_circle <- function(dat_p, color = color_Spectral()) {
  p <- dat_p %>% 
    ggplot2::ggplot() +
    ggforce::geom_circle(
      ggplot2::aes(x0 = x, y0 = y, fill = color.by, r = 0.5), 
      color = NA
    ) +
    ggplot2::scale_fill_gradientn(colors = color) +
    ggplot2::labs(fill = NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme_void()
  return(p)
}

#' Spatial Plotting: Point Plot for Continuous Data
#'
#' 绘制速度较快，用于预览参数设置是否正确。
#'
#' @param dat_p 数据框，至少包含三列：绘制点的x轴坐标、y轴坐标、用于映射绘制颜色
#' 的连续变量。
#' @param x 字符串，`dat_p`中存放点x轴坐标的列名。
#' @param y 字符串，`dat_p`中存放点y轴坐标的列名。
#' @param color.by 字符串，`dat_p`中存放需要以不同颜色绘制的连续变量的列名。
#' @param color 向量，映射连续变量梯度值的若干颜色。
#' @param y.rev 逻辑值，用于控制是否将y轴坐标顺序颠倒（在ggplot2中，较小的y值位
#' 于绘制图片的下方；而在空转数据meta.data中，较小的y值位于病理图片的上方）。
#'
#' @return 使用不同颜色绘制不同分类变量的点图。
#' @export
#'
st_plot_continuous_point <- function(dat_p, x = "col", y = "row", 
                                     color.by, color = color_Spectral(), 
                                     y.rev = TRUE) {
  dat_p %>%
    .st_plot_continuous_data(x, y, color.by, y.rev) %>%
    .st_plot_continuous_point(color)
}

#' Spatial Plotting: Circle Plot for Continuous Data
#'
#' 绘制速度较慢，使用点图预览确认参数设置无误后，用于正式图片的绘制。
#'
#' @param dat_p 数据框，至少包含三列：绘制点的x轴坐标、y轴坐标、用于映射绘制颜色
#' 的连续变量。
#' @param x 字符串，`dat_p`中存放点x轴坐标的列名。
#' @param y 字符串，`dat_p`中存放点y轴坐标的列名。
#' @param color.by 字符串，`dat_p`中存放需要以不同颜色绘制的连续变量的列名。
#' @param color 向量，映射连续变量梯度值的若干颜色。
#' @param y.rev 逻辑值，用于控制是否将y轴坐标顺序颠倒（在ggplot2中，较小的y值位
#' 于绘制图片的下方；而在空转数据meta.data中，较小的y值位于病理图片的上方）。
#'
#' @return 使用不同颜色绘制不同分类变量的相切圆图。
#' @export
#' 
st_plot_continuous_circle <- function(dat_p, x = "col", y = "row", 
                                      color.by, color = color_Spectral(), 
                                      y.rev = TRUE) {
  dat_p %>%
    .st_plot_continuous_data(x, y, color.by, y.rev) %>%
    .st_plot_continuous_circle(color)
}


#' Add Expression Information into MetaData
#'
#' @param sr_object Seurat对象。
#' @param feature 字符串或字符串向量，需要将表达量添加进meta.data的基因名或基因
#' 名向量。
#' @param slot 字符串，Specific assay data to get, one of "count", "data" 
#' (log-normalized counts), and "scale.data" (scaled data).
#' @param assay 字符串，Specific assay to get data from, defaults to the default 
#' assay.
#'
#' @return
#' @export
#'
sr_add_feature <- function(sr_object, feature, slot = "data", assay = NULL) {
  dat_feature <- SeuratObject::GetAssayData(
    object = sr_object, 
    slot = slot, 
    assay = assay
  )[feature, , drop = FALSE]
  df_feature <- as.data.frame(t(as.matrix(dat_feature)))
  sr_add <- SeuratObject::AddMetaData(
    object = sr_object, 
    metadata = df_feature, 
  )
  return(sr_add)
}

# Demo

if (FALSE) {
  sr_2938 <- read_rds("/data01/home/wuwenrui/CRLM_20240320/data_input/ST2938B_spotlight_Spatial_TLS.rds")
  
  # plot gene expression
  feature <- "ALB"
  color.by <- feature
  order <- c(
    "T&NK cell", 
    "B/Plasma cell", 
    "Myeloid cell", 
    "Fibroblast", 
    "Endothelial cell", 
    "Epithelial cell", 
    "Hepatocyte"
  )
  x <- "col"
  y <- "row"
  
  sr_2938 <- sr_2938 %>%
    NormalizeData() %>%
    sr_add_feature(feature = feature, slot = "data") 
  sr_2938@meta.data %>% 
    st_plot_continuous_point(color.by = color.by)
  
  # plot geneset score
  color.by <- "Hepatocyte"
  sr_2938@meta.data  %>%
    st_plot_continuous_point(color.by = color.by)
}


# Plot Saving #####

#' Size of Plot Saving
#' 
#' 由于在绘图时通过`ggplot2::coord_equal()`将图片的x轴和y轴设置成相等的尺度，如
#' 果在保存图片时不根据x轴和y轴的修改图片保存的尺寸，会出现不同图片中点的大小不
#' 一致、保存的图片周围出现大片空白等现象。
#' 
#' @param p 使用`st_plot_*()`系列函数绘制的ggplot2对象。
#' @param bin_unit 每长度单位（默认为"in"）对应的bin数量。默认为每有30个bin，保
#' 存图片时长和宽的尺寸增加1个长度单位。该参数越大，保存的图片中每个点越小。
#'
#' @return
#' @export
#'
st_plot_size <- function(p, bin_unit = 30) {
  dat_p <- p$data
  n_width <- max(dat_p$x) - min(dat_p$x)
  n_height <- max(dat_p$y) - min(dat_p$y)
  p_size <- list(width = n_width / bin_unit, height = n_height / bin_unit)
  return(p_size)
}

# Demo 

if (FALSE) {
  sr_2938 <- readr::read_rds("/data01/home/wuwenrui/CRLM_20240320/data_input/ST2938B_spotlight_Spatial_TLS.rds")
  
  color.by <- "Hepatocyte"
  p <- sr_2938@meta.data %>%
    st_plot_continuous_circle(color.by = color.by)
  
  ggplot2::ggsave(
    here::here("data_output/demo.png"), 
    plot = p, 
    dpi = 300, 
    bg = "white", 
    width = st_plot_size(p)$width, 
    height = st_plot_size(p)$height
  )
}


# Border of Region of Interest #####

#' Data for Arc Plotting 
#' 
#' 检测每个边缘同侧点和边缘对侧点的相邻关系（上/右/下/左），如果边缘同侧点和对侧
#' 点存在两个及以上连续的相邻点（相邻的对侧点不在同侧点的上-下/左-右位置），则在
#' 同侧点与对侧点相邻的位置绘制弧线。
#' 
#' @param df_ipsi 数据框，边缘同侧的一圈点，`x`和`y`列分别为点的x和y轴坐标。
#' @param df_cont 数据框，边缘对侧的一圈点，`x`和`y`列分别为点的x和y轴坐标。
#'
#' @return 包含绘制弧度相关参数的数据框。
#' 
.st_plot_contour_data_arc <- function(df_ipsi, df_cont) {
  df_arc <- vector("list", nrow(df_ipsi))
  cli::cli_progress_bar("Calculating Contour Arc", 
                        total = length(df_arc), clear = FALSE)
  for (i in seq_along(df_arc)) {
    df_point <- df_ipsi[i, ]
    df_adjac <- dplyr::bind_rows(
      df_point %>% dplyr::mutate(x = x, y = y + 1, tag = "top"), 
      df_point %>% dplyr::mutate(x = x, y = y - 1, tag = "bottom"), 
      df_point %>% dplyr::mutate(x = x + 1, y = y, tag = "right"), 
      df_point %>% dplyr::mutate(x = x - 1, y = y, tag = "left")
    )
    tag_adjac <- df_cont %>% 
      dplyr::left_join(df_adjac, by = c("x", "y")) %>% 
      dplyr::pull(tag) %>% 
      na.omit()
    status_adjac <- c("top", "right", "bottom", "left") %in% tag_adjac
    
    i_adjac <- which(c(status_adjac, status_adjac[1]))
    i_arc <- stringr::str_c(i_adjac, dplyr::lead(i_adjac)) %>% 
      purrr::map_chr(~ stringr::str_extract("12345", .x)) %>% 
      stringr::str_extract_all("[1-5]") 
    df_arc[[i]] <- i_arc %>% 
      purrr::map(~ {
        i_angle <- 0:3 * 1/2 * pi
        tibble::tibble(
          x0 = df_point$x, 
          y0 = df_point$y, 
          start = i_angle[min(as.numeric(.x))], 
          end = start + 1/2 * pi, 
          r = 0.5 
        )
      }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::filter(!is.na(start))
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  df_arc <- dplyr::bind_rows(df_arc)
  return(df_arc)
}


#' Data for Segment Plotting 
#' 
#' 检测每个边缘同侧点和边缘对侧点以及其他同侧点的相邻关系（上/右/下/左），同侧点
#' 和相邻对侧点的相邻位置为线段的潜在起始点，同侧点和相邻对侧点连线的垂直方向为
#' 线段的潜在延伸方向（仅在x轴和y轴方向进行延伸）；当同侧点在潜在延伸方向存在其
#' 他同侧点，则从潜在起始点向该方向绘制0.5个单位长度的线段。
#' 
#' @param df_ipsi 数据框，边缘同侧的一圈点，`x`和`y`列分别为点的x和y轴坐标。
#' @param df_cont 数据框，边缘对侧的一圈点，`x`和`y`列分别为点的x和y轴坐标。
#'
#' @return 包含绘制线段相关参数的数据框。
#' 
.st_plot_contour_data_segment <- function(df_ipsi, df_cont) {
  df_segment <- vector("list", nrow(df_ipsi))
  cli::cli_progress_bar("Calculating Contour Seg", 
                        total = length(df_segment), clear = FALSE)
  for (i in seq_along(df_segment)) {
    df_point <- df_ipsi[i, ]
    df_point_adjac <- .st_expand_adjacent_1(df_point, x = "x", y = "y")
    df_beg_ipsi <- df_point
    df_end_ipsi <- df_ipsi %>% 
      dplyr::semi_join(df_point_adjac, by = c("x", "y")) 
    
    df_diff_ipsi <- tidyr::expand_grid(
      df_beg_ipsi %>% dplyr::rename_with(~ str_c(.x, "_ipsi")), 
      df_end_ipsi %>% dplyr::rename_with(~ str_c(.x, "end_ipsi"))
    ) %>%
      dplyr::mutate(diff_x = xend_ipsi - x_ipsi, diff_y = yend_ipsi - y_ipsi) %>%
      dplyr::filter(abs(diff_x) + abs(diff_y) == 1)
    
    df_beg_cont <- df_cont %>% 
      dplyr::semi_join(df_point_adjac, by = c("x", "y"))
    list_beg_end <- vector("list", nrow(df_diff_ipsi))
    for (j in seq_along(list_beg_end)) {
      df_diff_cont <- df_beg_cont %>%
        dplyr::mutate(xend = x + df_diff_ipsi[j, ]$diff_x, 
                      yend = y + df_diff_ipsi[j, ]$diff_y) %>%
        dplyr::semi_join(df_cont, by = c("xend" = "x", "yend" = "y")) %>%
        dplyr::rename_with(~ str_c(.x, "_cont"))
      list_beg_end[[j]] <- dplyr::bind_cols(df_diff_ipsi[j, ], df_diff_cont)
    }
    df_beg_end <- dplyr::bind_rows(list_beg_end) 
    if (nrow(df_beg_end) > 0) {
      df_segment[[i]] <- df_beg_end %>%
        dplyr::mutate(
          x = (x_ipsi + x_cont) / 2, 
          y = (y_ipsi + y_cont) / 2, 
          xend = (xend_ipsi + xend_cont) / 2, 
          yend = (yend_ipsi + yend_cont) / 2, 
        ) %>%
        dplyr::distinct(x, y, xend, yend)
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  df_segment <- dplyr::bind_rows(df_segment)
  return(df_segment)
}


#' Data for Contour Plotting
#'
#' @param df_location 数据框，至少包含目标区域点的x轴和y轴坐标。
#' @param x 字符串，`df_location`中存放点x轴坐标的列名。
#' @param y 字符串，`df_location`中存放点y轴坐标的列名。
#'
#' @return 包含绘制弧度和线段相关参数的数据框。数据框中`tag`列为"arc"的为绘制弧
#' 线的相关参数，`tag`列为"segment"的为绘制线段的相关参数。
#'
.st_plot_contour_data <- function(df_location, x = "col", y = "row", y.rev = TRUE) {
  df <- df_location %>%
    dplyr::select(setNames(c(x, y), c("x", "y"))) %>%
    dplyr::mutate(dplyr::across(c(x, y), as.numeric))
  if (y.rev) {df$y <- -1 * df$y}
  df_expand_1p <- df %>% 
    st_expand_n(1, x = "x", y = "y") %>% 
    dplyr::filter(expand_n == 1) %>% 
    dplyr::select(x, y)
  df_outermost <- df %>% 
    dplyr::semi_join(
      .st_expand_adjacent_1(df_expand_1p, x = "x", y = "y"), 
      by = c("x", "y")
    )
  df_arc <- dplyr::bind_rows(
    .st_plot_contour_data_arc(df_outermost, df_expand_1p),
    .st_plot_contour_data_arc(df_expand_1p, df_outermost) ,
  ) %>%
    dplyr::mutate(tag = "arc")
  df_segment <- dplyr::bind_rows(
    .st_plot_contour_data_segment(df_outermost, df_expand_1p), 
  ) %>%
    dplyr::mutate(tag = "segment")
  df_plot <- dplyr::bind_rows(df_arc, df_segment) %>%
    dplyr::relocate(tag)
  return(df_plot)
}


#' Data for Contour Plotting (Method 1)
#'
#' @param df_location 数据框，至少包含目标区域点的x轴和y轴坐标。
#' @param x 字符串，`df_location`中存放点x轴坐标的列名。
#' @param y 字符串，`df_location`中存放点y轴坐标的列名。
#' @param group.by 字符串，`df_location`中定义不同目标区域的列名，当`group.by`为
#' NULL时，则将所有点视作属于同一个目标区域。
#' @param y.rev 逻辑值，用于控制是否将y轴坐标顺序颠倒（在ggplot2中，较小的y值位
#' 于绘制图片的下方；而在空转数据meta.data中，较小的y值位于病理图片的上方）。
#'
#' @return 数据框，包含可用于绘制目标区域边缘的规整数据。
#' @export
#'
st_plot_contour_1_data <- function(df_location, 
                                   x = "col", y = "row", 
                                   group.by = NULL, 
                                   y.rev = TRUE) {
  if (is.null(group.by)) {
    group.by <- "group.by"
    df_location <- df_location %>%
      dplyr::mutate(group.by = 1)
  }
  df_p <- df_location %>%
    dplyr::select(setNames(c(x, y, group.by), c("x", "y", "group.by"))) 
  groups <- unique(df_p$group.by)
  list_df_p <- vector("list", length(groups))
  for (i in seq_along(list_df_p)) {
    message(stringr::str_glue("Calculating Contour ({i}/{length(list_df_p)})"))
    list_df_p[[i]] <- df_p %>%
      dplyr::filter(group.by == groups[i]) %>%
      .st_plot_contour_data(x = "x", y = "y", y.rev = y.rev)
  }
  dat_p <- dplyr::bind_rows(list_df_p)
  return(dat_p)
}


#' Contour Plotting (Method 1)
#'
#' @param p 使用`st_plot_*()`系列函数绘制的ggplot2对象。
#' @param dat_p 数据框，经过`st_plot_contour_1_data()`函数预处理的规整数据。
#' @param linewidth 数值，用于控制绘制的边缘的宽度。
#' @param color 字符串，用于控制绘制的边缘的颜色。
#'
#' @return 根据边缘绘制参数添加边缘的ggplot2对象。
#'
st_plot_contour_1 <- function(p, dat_p, linewidth = 1, color = "black") {
  layer_contour <- list(
    ggforce::geom_arc(
      ggplot2::aes(x0 = x0, y0 = y0, start = start, end = end, r = r),
      data = dat_p %>% dplyr::filter(tag == "arc"), 
      linewidth = linewidth, color = color
    ), 
    ggplot2::geom_segment(
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend), 
      data = dat_p %>% dplyr::filter(tag == "segment"), 
      linewidth = linewidth, color = color
    )
  )
  p <- p + layer_contour
  return(p)
}

# Demo 

if (FALSE) {
  setwd("~/ST_std/")
  dat_contour <- readr::read_csv("000.Demo_data/ST2938B_location.csv")
  df_location <- dat_contour %>% 
    dplyr::filter(TLS %in% c("ST2938B_2", "ST2938B_8", "ST2938B_9"))
  dat_contour_p <- df_location %>% st_plot_contour_1_data(group.by = "TLS")
  dat_contour_p2 <- df_location %>% st_plot_contour_1_data()
  df_location %>%
    st_plot_discrete_point(color.by = "TLS") %>% 
    st_plot_contour_1(dat_contour_p, linewidth = 1, color = "red") %>% 
    st_plot_contour_1(dat_contour_p2, linewidth = 0.5, color = "blue") 
}
