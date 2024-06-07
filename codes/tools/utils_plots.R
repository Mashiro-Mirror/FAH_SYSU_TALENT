myggsave <- function (filename, width = 6, height = 6, ...) {
    args <- list(...)
    if ("path" %in% names(args)) {
        outpath <- args$path
    }
    else {
        outpath <- dirname(filename)
    }
    if (!dir.exists(outpath))
        dir.create(outpath, recursive = TRUE)
    outfiles <- gsub(".(png|pdf)$", "", filename)
    outfiles <- paste0(outfiles, c(".png", ".pdf"))
    ggsave(outfiles[1], width = width, height = height, ...)
    ggsave(outfiles[2], width = width, height = height, ...)
}

basesave <- function (filename, func, width = 6, height = 6, moreArgs = NULL) {
    outpath <- dirname(filename)
    if (!dir.exists(outpath))
        dir.create(outpath, recursive = TRUE)
    outfiles <- gsub(".(png|pdf)$", "", filename)
    outfiles <- paste0(outfiles, c(".png", ".pdf"))
    pdf(width = width, height = height, file = outfiles[2])
    do.call(func, moreArgs)
    dev.off()
    png(filename = outfiles[1], width = width, height = height,
        units = "in", res = 300)
    do.call(func, moreArgs)
    dev.off()
}


featcolor <- function() {
    scale_color_gradientn(colours = BuenColors::jdb_palette("brewer_spectra", type = "continuous"))
}

removeAxis <- function(element.x = NULL, element.y = NULL) {
    elements_all <- c("title", "text", "line", "ticks")
    element.x <- intersect(element.x, elements_all)
    element.y <- intersect(element.y, elements_all)
    
    args <- list(validate = TRUE)
    for (item in element.x) {
        args[[paste0("axis.", item, ".x")]] <- element_blank()
    }
    for (item in element.y) {
        args[[paste0("axis.", item, ".y")]] <- element_blank()
    }
    
    if (length(args) == 1) {
        stop("`element.x` or `element.y` must have at least one value in ['title', 'text', 'line', 'ticks'].")
    }
    no.axis.theme <- do.call(theme, args)
    return(no.axis.theme)
}


getLegend <- function(plt) {
    ggplotify::as.ggplot(cowplot::get_legend(plt))    
}


continuecolors <- function(n, name, prefix = 0, suffix = 5) {
    colors <- as.character(BuenColors::jdb_palette(name, n, type = "continuous"))
    idx <- seq(from = 1, to = length(colors), length.out = prefix + n + suffix)
    colors[idx[(prefix + 1):(prefix + n)]]
}

largerBolderTitle <- function(dtype, fontsize = 18, face = "bold", ...) {
    if (missing(dtype)) dtype <- "x"
    dtype <- intersect(dtype, c("x", "y", "title"))
    stopifnot(length(dtype) > 0)
    
    themeargs <- list()
    boldtext <- element_text(size = fontsize, face = face, ...)
    for (item in dtype) {
        if (item == "title") {
            themeargs[["plot.title"]] <- boldtext
        } else {
            themeargs[[paste0("axis.title.", item)]] <- boldtext
        }
    }
    returntheme <- do.call(theme, args = themeargs)
    returntheme
}


plot_sigbox <- function(cells_subtype, x = "subtype", y = "ratio_in_immune", color = "treatment", yoffset = 0) {
    data <- cells_subtype[, c(x, y, color)]
    colnames(data) <- c("x", "y", "group")
    if (!is.factor(data$x)) data$x <- as.factor(data$x)
    
    pvalue <- sapply(split(data, data$x), function(df) {
        res <- wilcox.test(y ~ group, data = df)
        res$p.value
    })
    
    pos_y <- sapply(split(data, data$x), function(df) {
        val <- sapply(split(df, df$group), function(dt) {
            boxplot(dt$y, plot = FALSE)$stats[5, ]
        })
        val * (1 + 0.2)
    })
    pos_y <- max(pos_y) + yoffset
    labels <- symnum(pvalue, corr = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "n.s."))
    
    #xval <- unique(data$x)
    xval <- 1:length(levels(data$x))
    
    gg <- ggplot(data, aes(x = x, y = y, color = group)) +
        geom_boxplot(outlier.alpha = NA) +
        geom_point(position = position_dodge(width = 0.75), aes(group = group)) +
        ggsignif::geom_signif(y_position = pos_y, xmin = xval - .2, xmax = xval + .2, annotations = labels) +
        cowplot::theme_cowplot()
    
    list(plot = gg, pval = pvalue)
}


rank_plot <- function (markers, topn = 10, pval.column = "p_val_adj", logfc.column = "avg_logFC") {
    markers$pval <- markers[, pval.column]
    markers$logfc <- markers[, logfc.column]
    markers <- markers[order(markers$logfc, decreasing = TRUE),
    ]
    markers$rank <- 1:nrow(markers)
    markers_label <- rbind(head(markers, topn), tail(markers,topn))
    ggplot(markers, aes(x = rank, y = logfc, color = -log10(pval + 1e-05))) +
        geom_point() + 
        scale_color_gradient2(
            low = "gray", mid = "yellow", high = "red",
            midpoint = -log10(0.05 + 1e-05), 
            name = sprintf("-log10(%s+1e-05)", pval.column)
        ) +
        ggrepel::geom_text_repel(
            data = markers_label, aes(label = gene),
            color = "blue", max.overlaps = 50
        ) + theme_bw() +
        theme(
            legend.position = c(0.05, 0.05),
            legend.justification = c("left", "bottom"),
            legend.direction = "horizontal",
            legend.background = element_blank()
        ) + labs(x = "Rank", y = logfc.column)
}

barPlot <- function(data, group, color, value, ptype = c("value", "percent"), label = FALSE, label.color = NULL) {
    ptype <- match.arg(ptype)
    df <- data[, c(group, color, value)]
    colnames(df) <- c("group", "color", "value")
    df <- df %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(percent = value / sum(value),
                      pos.value = cumsum(value) - 0.5 * value,
                      pos.percent = 1 - cumsum(percent) + 0.5 * percent,
                      label = paste0(value, "(", round(percent * 100), "%)"))
    
    if (ptype == "value") {
        df$pos <- df$pos.value
    } else {
        df$pos <- df$pos.percent
    }
    
    gg <- ggplot(df, aes_string(x = "group", y = ptype, fill = "color")) +
        geom_col(position = "stack") +
        labs(x = group, y = ptype) +
        theme_bw()
    
    if (label) {
        df_label <- df
        if (!is.null(label.color)) {
            df_label$label[!df_label$color %in% label.color] <- ""
        }
        gg <- gg + geom_text_repel(aes(x = group, y = pos, label = label), data = df_label, angle = 90)
    }
    
    gg
}

heatmap_genesets <- function(scrna, genesets, group = "seurat_clusters", 
                             center = TRUE, labelcolumn = NULL, splitcolumn = NULL) {
    
    
    hm <- DoHeatmap(scrna, features = unlist(genesets), group.by = group, combine = FALSE)
    hm_data <- reshape2::dcast(hm[[1]]$data %>% dplyr::filter(Cell %in% colnames(scrna)), Feature ~ Cell, value.var = "Expression")
    rownames(hm_data) <- hm_data$Feature
    
    anno_col <- scrna@meta.data[colnames(hm_data)[-1], unique(c(group, labelcolumn)), drop = FALSE]
    anno_row <- unlist(genesets) %>% {
        inner <- .
        data.frame(
            category = gsub("\\d+$", "", names(inner)),
            gene = inner,
            row.names = inner
        )
    }
    anno_row <- anno_row[hm_data$Feature, ]
    
    if (!is.null(splitcolumn)) {
        splitcolumn <- scrna@meta.data[colnames(hm_data)[-1], splitcolumn, drop = FALSE]
    }
    
    cols <- lapply(colnames(anno_col), function(col) {
        if (col == "seurat_clusters") values <- as.character(sort(unique(scrna@meta.data[, col])))
        else values <- as.character(unique(scrna@meta.data[, col]))
        setNames(scales::hue_pal()(length(values)), nm = values)
    })
    names(cols) <- colnames(anno_col)
    anno_col <- ComplexHeatmap::columnAnnotation(
        df = anno_col,
        col = cols
    )
    
    hm <- ComplexHeatmap::Heatmap(
        matrix = as.matrix(hm_data[anno_row$gene, -1]),
        name = "",
        col = circlize::colorRamp2(seq(-4, 4, length.out = 50), Seurat::PurpleAndYellow()),
        #left_annotation = anno_cate,
        row_names_side = "right",
        row_split = anno_row$category,
        row_title_rot = 00,
        cluster_row_slices = TRUE,
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        # column annotation
        column_names_side = "top",
        column_split = splitcolumn,
        show_column_names = FALSE,
        top_annotation = anno_col,
        #column_title = NULL,
        column_title_rot = 00,
        cluster_column_slices = identical(group, splitcolumn),
        cluster_columns = FALSE
    )
    #ComplexHeatmap::draw(hm)
    ggplotify::as.ggplot(hm)
}


plot_velo_emb <- function(scrna, group.by, arrows, reduction, cols = NULL, pt.size = NULL) {
    DimPlot(scrna, reduction = reduction, group.by = group.by, cols = cols, pt.size = pt.size) +
        geom_segment(
            data = arrows, 
            aes(x = x0, y = y0, xend = x1, yend = y1), 
            inherit.aes = F,
            arrow = arrow(length = unit(arrows$alen, "inches"))
        )
}


circleHighlight <- function(
    df, highlight, from = "from", to = "to", 
    ordvalue = NULL, interacting_pair = "interacting_pair", 
    colormaps = NULL, arr.length = 0.1, arr.alpha = 1, trackheights = NULL
) {
    require(circlize)
    data <- df[, c(from, to, interacting_pair)]
    colnames(data) <- c("from", "to", "pair")
    
    if (!is.null(ordvalue)) {
        data <- data %>%
            dplyr::mutate(from = factor(from, levels = ordvalue)) %>%
            dplyr::arrange(from) %>%
            dplyr::mutate(from = as.character(from))
    }
    
    if (is.null(trackheights)) {
        trackheights <- c(0.25, 0.015)
    }
    
    allvalues <- unique(c(data$from, data$to))
    
    from <- data
    to <- data
    from$sectors <- from$from
    to$sectors <- to$to
    from$dtype <- "ligand"
    to$dtype <- "receptor"
    
    data2 <- rbind(from, to)
    data2$x <- rnorm(n = nrow(data2), mean = 2)
    data2$y <- rnorm(n = nrow(data2), mean = 2)
    for (val in allvalues) {
        data2$x[data2$sectors == val] <- 1:sum(data2$sectors == val)
    }
    
    sector.width <- table(data2$sectors)
    
    circos.initialize(data2$sectors, x = data2$x, sector.width = sector.width)
    
    # from or to labels
    circos.track(
        data2$sectors,
        y = data2$y,
        track.height = trackheights[1],
        bg.col = NA,
        bg.border = NA,
        panel.fun = function(x, y) {
            circos.text(
                x = CELL_META$xcenter,
                y = CELL_META$ylim[1],
                labels = CELL_META$sector.index,
                niceFacing = TRUE,
                facing = "clockwise",
                adj = c(0, 0.5)
            )
        }
    )
    
    # color block for "from" or "to" 
    circos.track(
        data2$sectors,
        y = data2$y,
        track.height = trackheights[2],
        cell.padding = rep(0, 4),
        bg.col = setNames(nm = allvalues, BuenColors::jdb_palette("corona", n = length(allvalues))),
        bg.border = NA
    )
    
    # links for "from" to "to"
    dthigh <- subset(data2, subset = pair %in% highlight)
    if (is.null(colormaps)) {
        colormaps <- setNames(nm = highlight, BuenColors::jdb_palette("corona", length(highlight)))
    }
    dthigh$colors <- colormaps[dthigh$pair]
    dthigh$colors <- grDevices::adjustcolor(dthigh$colors, alpha.f = arr.alpha)
    
    fromhigh <- subset(dthigh, subset = dtype == "ligand")
    tohigh <- subset(dthigh, subset = dtype == "receptor")
    
    for (idx in 1:nrow(fromhigh)) {
        circos.link(
            fromhigh$sectors[idx], fromhigh$x[idx], tohigh$sectors[idx], tohigh$x[idx],
            # arr.width = 0.05, arr.length = 0.05,
            arr.length = arr.length,
            directional = 1, col = fromhigh$colors[idx]
        )
    }
    circos.clear()
}


myvolcano <- function(
    markers, title = "", pval = 0.05, logfc = 1, 
    label.size = 3, adjust.pval = 5e-5, topn = 10
) {
    addsig <- function(mk, logfc = 1, pval = 0.05) {
        mk$significance <- ifelse(
            mk$avg_logFC >= logfc & mk$p_val_adj <= pval, "UP",
            ifelse(
                mk$avg_logFC <= -logfc & mk$p_val_adj <= pval, "DOWN",
                "Not Sig."
            ))
        mk
    }
    
    df <- markers[, c("gene", "avg_logFC", "p_val_adj")]
    df <- addsig(df, pval = pval, logfc = logfc)
    cols <- c("UP" = "red", "DOWN" = "blue", "Not Sig." = "lightgray")
    
    posvalue <- function(x) {
        ifelse(x <= 0, -x, x)
    }
    
    df2 <- dplyr::filter(df, significance != "Not Sig.") %>%
        dplyr::group_by(significance) %>%
        dplyr::top_n(wt = abs(avg_logFC), n = topn)
    
    ggplot(df, aes(avg_logFC, -log10(p_val_adj + adjust.pval), color = significance)) +
        geom_point() +
        scale_x_continuous(labels = posvalue) +
        scale_color_manual(values = cols) +
        geom_hline(yintercept = -log10(pval + adjust.pval), lty = 4,col = "grey",lwd = 0.6) +
        geom_vline(xintercept = c(-logfc, logfc), lty = 4,col = "grey",lwd = 0.6) +
        theme_bw() +
        theme(
            #legend.position = "none",
            panel.grid = element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14)
        ) +
        ggrepel::geom_text_repel(
            aes(label = gene), data = df2,
            show.legend = FALSE, 
            color = "black", size = label.size,
            segment.alpha = 0.5, segment.color = "lightgray"
        ) + labs(
            x = "avg_logFC",y = "-log10 (p_val_adj)", title = title,
            subtitle = sprintf("logFC = %s, p.adj = %s", logfc, pval)
        )
}
