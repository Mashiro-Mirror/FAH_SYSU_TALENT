add_pptx_images <- function(doc, imgfiles, imgdir, is.pdf = T, layout = c("1x1", "2x1", "2x2", "3x2", "3x3", "4x3", "4x4")) {
  layout_infos <- officer::layout_summary(doc)
  pptx_wh <- c(width = 10, height = 7.5) #inch
  # image inches: width = 6.5, height = 5
  positions <- list(
    "1x1" = c(width = 19.5, height = 15),
    "2x1" = c(width = 15.6, height = 12),
    "2x2" = c(width = 11.7, height = 9),
    "3x2" = c(width = 10.4, height = 8),
    "3x3" = c(width = 7.8, height = 6),
    "4x3" = c(width = 5.85, height = 4.5),
    "4x4" = c(width = 5.85, height = 4.5)
  )
  cm2inch <- 0.3937
  layout <- match.arg(layout)
  wh <- positions[[layout]] * cm2inch
  layout <- as.numeric(unlist(strsplit(layout, "x")))
  
  gap_x <- (pptx_wh[1] - layout[1] * wh[1]) / (layout[1] + 1)
  gap_y <- (pptx_wh[2] - layout[2] * wh[2]) / (layout[2] + 1)
  
  pos_x <- gap_x * (1:layout[1]) + wh[1] * (1:layout[1] - 1)
  pos_y <- gap_y * (1:layout[2]) + wh[2] * (1:layout[2] - 1)
  
  total <- layout[1] * layout[2]
  mat <- matrix(1:total, byrow = TRUE, ncol = layout[1])
  
  if (is.pdf){
    for (img in imgfiles){
      pdftools::pdf_convert(img, format = "png",dpi = 150)
    }
    imgfiles = dir(imgdir) %>% .[grep("png",.)]
    for (idx in 1:length(imgfiles)) {
      idx2 <- (idx - 1) %% total
      if (idx2 == 0) {
        doc <- officer::add_slide(doc, layout = layout_infos[1, 1], layout_infos[1, 2])
      }
      
      idx2 <- which(mat == (idx2 + 1), arr.ind = TRUE)
      
      doc <- officer::ph_with(
        x = doc,
        value = officer::external_img(imgfiles[idx], width = wh[1], height = wh[2]),
        use_loc_size = TRUE,
        location = officer::ph_location(
          left = pos_x[idx2[1, "col"]],
          top = pos_y[idx2[1, "row"]],
          width = wh[1],
          height = wh[2]
        )
      )
    }
    system("rm *.png")
  } else {
    for (idx in 1:length(imgfiles)) {
      idx2 <- (idx - 1) %% total
      if (idx2 == 0) {
        doc <- officer::add_slide(doc, layout = layout_infos[1, 1], layout_infos[1, 2])
      }
      
      idx2 <- which(mat == (idx2 + 1), arr.ind = TRUE)
      
      doc <- officer::ph_with(
        x = doc,
        value = officer::external_img(imgfiles[idx], width = wh[1], height = wh[2]),
        use_loc_size = TRUE,
        location = officer::ph_location(
          left = pos_x[idx2[1, "col"]],
          top = pos_y[idx2[1, "row"]],
          width = wh[1],
          height = wh[2]
        )
      )
    }
  }
  doc
}