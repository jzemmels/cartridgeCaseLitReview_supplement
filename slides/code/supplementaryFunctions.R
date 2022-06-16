cmcPlot_color <- function(reference,
                          target,
                          cmcClassifs,
                          type = "faceted",
                          cmcCol = "originalMethod",
                          corrCol = "pairwiseCompCor"){

  #check that the necessary columns are in cmcClassifs

  stopifnot("Make sure that there is a column called 'cellHeightValues' that is the result of the comparison_alignedTargetCell() function." = any(stringr::str_detect(names(cmcClassifs),"cellHeightValues")))

  stopifnot("Make sure that there is a column called 'alignedTargetCell' that is the result of the comparison_alignedTargetCell() function." = any(stringr::str_detect(names(cmcClassifs),"alignedTargetCell")))

  stopifnot("Make sure that there is a column called 'cellIndex'" = any(stringr::str_detect(names(cmcClassifs),"cellIndex")))

  stopifnot("Make sure that there is a column called 'theta'" = any(stringr::str_detect(names(cmcClassifs),"theta")))

  stopifnot(any(stringr::str_detect(names(cmcClassifs),cmcCol)))

  stopifnot(any(stringr::str_detect(names(cmcClassifs),corrCol)))

  # get the indices for the necessary columns
  referenceCellCol <- which(stringr::str_detect(names(cmcClassifs),"cellHeightValues"))

  targetCellCol <- which(stringr::str_detect(names(cmcClassifs),"alignedTargetCell"))

  cellIndexCol <- which(stringr::str_detect(names(cmcClassifs),"cellIndex"))

  thetaCol <- which(stringr::str_detect(names(cmcClassifs),"theta"))

  cmcIndexCol <- which(stringr::str_detect(names(cmcClassifs),cmcCol))

  # cmcClassifs <- cmcClassifs %>%
  #   dplyr::group_by(cellIndex) %>%
  #   dplyr::filter(!!as.name(corrCol) == max(!!as.name(corrCol)))

  targetCellData <- cmcClassifs %>%
    dplyr::select(dplyr::all_of(c(targetCellCol,cellIndexCol,thetaCol,cmcIndexCol))) %>%
    purrr::pmap_dfr(~ cmcR:::targetCellCorners(alignedTargetCell = ..1,
                                               cellIndex = ..2,
                                               theta = ..3,
                                               cmcClassif = ..4,
                                               target = target))

  referenceCells <- cmcClassifs %>%
    dplyr::pull(referenceCellCol)

  cellData <- cmcClassifs %>%
    dplyr::select(dplyr::all_of(c(cellIndexCol,referenceCellCol,cmcIndexCol))) %>%
    purrr::pmap_dfr(~ {

      cellInds <- ..2$cmcR.info$cellRange %>%
        stringr::str_remove("rows: ") %>%
        stringr::str_remove("cols: ") %>%
        stringr::str_split(pattern = ", ")

      cellInds_rows <- stringr::str_split(cellInds[[1]][1]," - ")[[1]]
      cellInds_cols <- stringr::str_split(cellInds[[1]][2]," - ")[[1]]

      return(data.frame(rowStart = as.numeric(cellInds_rows[1]),
                        rowEnd = as.numeric(cellInds_rows[2]),
                        colStart = as.numeric(cellInds_cols[1]),
                        colEnd = as.numeric(cellInds_cols[2])) %>%
               dplyr::mutate(cellIndex = ..1,
                             cmcClassif = ..3))

    }) %>%
    dplyr::mutate(rowStart = max(.data$rowEnd) - .data$rowStart,
                  rowEnd = max(.data$rowEnd) - .data$rowEnd,
                  colMean = purrr::map2_dbl(.data$colStart,.data$colEnd,~ mean(c(.x,.y))),
                  rowMean = purrr::map2_dbl(.data$rowStart,.data$rowEnd,~ mean(c(.x,.y))))

  # ggplot2 complains about the guides
  suppressWarnings({

    refPlt <- x3pListPlot(list("reference" = reference),
                          height.colors =
                            rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b'))) +
      ggplot2::guides(fill = "none") +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_rect(data = cellData,
                         ggplot2::aes(xmin = .data$colStart,xmax = .data$colEnd,ymin = .data$rowStart,ymax = .data$rowEnd,fill = .data$cmcClassif),
                         alpha = .2,
                         inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = c("#313695","#a50026")) +
      ggplot2::geom_text(data = cellData,
                         ggplot2::aes(x = .data$colMean,y = .data$rowMean,label = .data$cellIndex),inherit.aes = FALSE) +
      ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
      ggplot2::theme(
        legend.direction = "horizontal"
      ) +
      ggplot2::labs(fill = "CMC Classif.")

    cmcLegend <- ggplotify::as.ggplot(cowplot::get_legend(refPlt)$grobs[[1]])

    refPlt <- refPlt +
      ggplot2::theme(legend.position = "none")

    plt <- x3pListPlot(list("target" = target),
                       height.colors =
                         rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b'))) +
      ggplot2::theme(legend.position = "none")

    plt <- plt +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_raster(data = targetCellData,
                           ggplot2::aes(x = .data$x,y = .data$y,fill = .data$cmcClassif),
                           alpha = .2) +
      ggplot2::scale_fill_manual(values = c("#313695","#a50026")) +
      ggplot2::geom_text(data = targetCellData %>%
                           dplyr::group_by(.data$cellIndex) %>%
                           dplyr::summarise(x = mean(.data$x),
                                            y = mean(.data$y),
                                            theta = unique(.data$theta)),
                         ggplot2::aes(x=.data$x,y=.data$y,label = .data$cellIndex,angle = -1*.data$theta))

  })

  # library(patchwork)
  # return((refPlt | plt))
  if(type == "list"){
    return(list("reference" = refPlt,
                "target" = plt,
                "legend" = cmcLegend))
  }

  return(patchwork::wrap_plots(refPlt,plt,cmcLegend,nrow = 2,heights = c(1,.1)))
}

fiveplot <- function(comparisonResults, referenceScan, targetScan, cell) {

  referenceCell <- comparisonResults %>%
    filter(comparisonName == sprintf("%s vs. %s", referenceScan, targetScan) &
             cellIndex == cell) %>%
    pull(cellHeightValues) %>%
    .[[1]]

  targetCell <- comparisonResults %>%
    filter(comparisonName == sprintf("%s vs. %s", referenceScan, targetScan) &
             cellIndex == cell) %>%
    pull(alignedTargetCell) %>%
    .[[1]]

  referenceCell$surface.matrix <- (referenceCell$surface.matrix*referenceCell$cmcR.info$scaleByVal + referenceCell$cmcR.info$centerByVal)#*1e6
  targetCell$surface.matrix <- (targetCell$surface.matrix*targetCell$cmcR.info$scaleByVal + targetCell$cmcR.info$centerByVal)#*1e6


  patchComparisonPlts <- impressions::x3pComparisonPlot(
    reference = referenceCell,
    target = targetCell,
    cutoffThresh = sd(c(c(referenceCell$surface.matrix),c(targetCell$surface.matrix)),na.rm = TRUE),
    plotNames = c(sprintf("%s Cell %s", referenceScan, cell),
                  sprintf("%s Aligned Cell", targetScan),
                  "Filtered Element-wise Average",
                  sprintf("%s Cell %s\nFiltered Differences", referenceScan, cell),
                  sprintf("%s Aligned Cell\nFiltered Differences", targetScan))
  )

  patchComparisonLegend_match <-
    cowplot::plot_grid(patchComparisonPlts$legend$grobs[[1]])

  combinedValues <-  referenceCell %>%
    impressions::x3pToDF() %>%
    rename(refValue = value) %>%
    left_join(targetCell %>%
                impressions::x3pToDF() %>%
                rename(targValue = value),
              by = c("x","y"))

  blobBoundaries <- comparisonResults %>%
    filter((comparisonName == sprintf("%s vs. %s", referenceScan, targetScan) &
              cellIndex == cell)) %>%
    as.data.frame() %>%
    dplyr::select(comparisonName,cellIndex,cellHeightValues,alignedTargetCell) %>%
    pmap(~ {

      reference <-  ..3

      target <- ..4

      reference$surface.matrix <- (reference$surface.matrix*reference$cmcR.info$scaleByVal + reference$cmcR.info$centerByVal)*1e6
      target$surface.matrix <- (target$surface.matrix*target$cmcR.info$scaleByVal + target$cmcR.info$centerByVal)*1e6

      averageBinarized <- bind_rows(reference %>%
                                      impressions::x3pToDF() %>%
                                      mutate(value = value),
                                    target %>%
                                      impressions::x3pToDF() %>%
                                      mutate(value = value)) %>%
        group_by(x,y) %>%
        summarise(difference = diff(value),
                  absDifference = abs(diff(value)),
                  average = mean(value),
                  .groups = "drop")  %>%
        mutate(comparisonName = ..1,
               cellIndex = ..2)%>%
        mutate(value = ifelse(absDifference > sd(c(c(reference$surface.matrix),c(target$surface.matrix)),na.rm = TRUE),TRUE,FALSE))

      suppressWarnings({

        averageMat <- averageBinarized %>%
          mutate(x = x+1,
                 y=y+1) %>%
          as.data.frame() %>%
          dplyr::select(x,y,value) %>%
          imager::as.cimg() %>%
          as.matrix()

      })

      averageMat[is.na(averageMat)] <- 0

      # we pad the matrix so that the contours one the edge blobs are properly
      # identified. the padding is removed in the last lines of the creation of
      # the outline object below
      averageMat  <- averageMat %>%
        imager::as.cimg() %>%
        imager::pad(nPix = 10,axes = "xy",val = 0)

      labels <- imager::label(averageMat)

      bounds <- map(unique(labels[labels > 0]),
                    function(lab){

                      imager::boundary(labels == lab)

                    })

      return(list(bounds,labels))

    })

  # combine all labeled blobs into one image
  boundaryPx <- Reduce("+",blobBoundaries[[1]][[1]] %>%
                         map(as.matrix)) %>%
    imager::as.cimg()

  # the mask used to dilate the blobs will grow them towards the bottom-right of
  # the matrix
  dilatedPx <- imager::dilate_rect(boundaryPx,sx = 2,sy = 2)
  dilatedPx_labels <- imager::dilate_rect(blobBoundaries[[1]][[2]],sx = 2,sy = 2)

  # flip the image and re-apply the dilation to grow the borders to the other
  # corners. flip back after dilation
  dilatedPx_mirrorx <- imager::mirror(imager::dilate_rect(imager::mirror(boundaryPx,axis="x"),sx = 2,sy = 2),axis="x")
  dilatedPx_mirrorx_labels <- imager::mirror(imager::dilate_rect(imager::mirror(blobBoundaries[[1]][[2]],axis="x"),sx = 2,sy = 2),axis="x")

  dilatedPx_mirrory <- imager::mirror(imager::dilate_rect(imager::mirror(boundaryPx,axis="y"),sx = 2,sy = 2),"y")
  dilatedPx_mirrory_labels <- imager::mirror(imager::dilate_rect(imager::mirror(blobBoundaries[[1]][[2]],axis="y"),sx = 2,sy = 2),"y")

  dilatedPx_mirrorxy <- imager::mirror(imager::dilate_rect(imager::mirror(boundaryPx,axis="xy"),sx = 3,sy = 3),"xy")
  dilatedPx_mirrorxy_labels <- imager::mirror(imager::dilate_rect(imager::mirror(blobBoundaries[[1]][[2]],axis="xy"),sx = 3,sy = 3),"xy")

  # combine all of the dilated images together into one image
  dilatedPx_comb <- dilatedPx + dilatedPx_mirrorx + dilatedPx_mirrory + dilatedPx_mirrorxy

  # we just want a binary labeling
  dilatedPx_comb[dilatedPx_comb > 0] <- 1

  # the dilated boundaries will have also grown into the blobs, so we take those
  # pixels out
  dilatedPx_comb[blobBoundaries[[1]][[2]] > 0] <- 0

  # from: https://stackoverflow.com/questions/34756755/plot-outline-around-raster-cells
  outline <- dilatedPx_comb %>%
    as.data.frame() %>%
    filter(value > 0) %>%
    mutate(x = x-1,
           y = y-1) %>%
    raster::rasterFromXYZ() %>%
    raster::rasterToPolygons(dissolve = TRUE) %>%
    fortify() %>%
    #the boundaries around the filtered blobs all share a common value in the
    #"hole" column of TRUE
    filter(hole) %>%
    # remove padding used previously
    mutate(lat = lat-5,
           long = long-5)


  `-.gg` <- function(plot, layer) {
    if (missing(layer)) {
      stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
    }
    if (!is.ggplot(plot)) {
      stop('Need a plot on the left side')
    }
    plot$layers = c(layer, plot$layers)
    plot
  }

  topLeft <- patchComparisonPlts[[1]] +
    cowplot::theme_nothing() +
    labs(subtitle = sprintf("%s Cell %s", referenceScan, cell)) +
    theme(plot.margin = margin(0,0,5,0),
          plot.subtitle = element_text(hjust = .5,size = 8,vjust = -1)) +
    geom_raster(data = combinedValues %>%
                  filter(is.na(refValue) & !is.na(targValue)),
                fill = "gray40")

  bottomLeft <-patchComparisonPlts[[2]] +
    cowplot::theme_nothing() +
    labs(subtitle = paste0(targetScan," Aligned Cell\nat ",unique(comparisonResults$theta),"°")) +
    theme(plot.margin = margin(-20,-100,30,-100),
          plot.subtitle = element_text(hjust = .5,vjust = -78,size = 8)) +
    geom_raster(data = combinedValues %>%
                  filter(!is.na(refValue) & is.na(targValue)),
                fill = "gray40")

  middle <- patchComparisonPlts[[3]] +
    cowplot::theme_nothing() +
    labs(subtitle = "Filtered Element-wise Average\nAbs. Differences at Most 1") +
    theme(plot.margin = margin(0,25,0,25),
          plot.subtitle = element_text(hjust = .5,size = 8,vjust = -5)) -
    geom_raster(fill = "gray80") +
    geom_path(data = outline,  color = "grey40",
              aes(x=long,y=lat,group=group),
              colour = "gray40",
              inherit.aes = FALSE,
              size = .2)

  topRight <- patchComparisonPlts[[4]] +
    cowplot::theme_nothing() +
    labs(subtitle = sprintf("Filtered %s Cell %s\nAbs. Differences Greater Than 1",referenceScan, cell)) +
    theme(plot.margin = margin(0,0,5,0),
          plot.subtitle = element_text(hjust = .5,size = 8)) -
    geom_raster(fill = "gray80") +
    geom_path(data = outline,  color = "grey40",
              aes(x=long,y=lat,group=group),
              colour = "gray40",
              inherit.aes = FALSE,
              size = .1)

  bottomRight <- patchComparisonPlts[[5]] +
    cowplot::theme_nothing() +
    labs(subtitle = sprintf("Filtered %s Aligned Cell\nAbs. Differences Greater Than 1",targetScan)) +
    theme(plot.margin = margin(-20,-100,30,-100),
          plot.subtitle = element_text(hjust = .5,vjust = -78,size = 8)) -
    geom_raster(fill = "gray80") +
    geom_path(data = outline, color = "grey40",
              aes(x=long,y=lat,group=group),
              colour = "gray40",
              inherit.aes = FALSE,
              size = .1)

  design <- "ACCD\nBCCE"

  patchwork::wrap_plots(topLeft,bottomLeft,middle,topRight,bottomRight,design = design) +
    inset_element(patchComparisonLegend_match,left = -2.15,bottom = 0,right = -2.15,top = 0,on_top = FALSE,align_to = 'full')
}
