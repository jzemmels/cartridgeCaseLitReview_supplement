x3pListPlot_vert <- function(x3pList,
                             type = "faceted",
                             legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                             height.quantiles = c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                             height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                             na.value = "gray65"){

  if(purrr::is_empty(names(x3pList))){
    x3pList <- setNames(x3pList,paste0("x3p",1:length(x3pList)))
  }

  if(type == "faceted"){
    surfaceMat_df <- purrr::pmap_dfr(.l = list(x3pList,
                                               names(x3pList)),
                                     function(x3p,name){

                                       x3p$header.info$incrementX <- 1
                                       x3p$header.info$incrementY <- 1

                                       x3p %>%
                                         x3ptools::x3p_to_df() %>%
                                         #perform some transformations on the
                                         #x,y values so that the plot is
                                         #representative of the actual surface
                                         #matrix (i.e., element [1,1] of the
                                         #surface matrix is in the top-left
                                         #corner)
                                         dplyr::mutate(xnew = max(.data$y) - .data$y,
                                                       ynew = max(.data$x) - .data$x,
                                                       value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                                         dplyr::select(-c(.data$x,.data$y)) %>%
                                         dplyr::rename(x=.data$xnew,
                                                       y=.data$ynew) %>%
                                         dplyr::mutate(x3p = rep(name,times = nrow(.)))
                                     }) %>%
      dplyr::mutate(x3p = factor(.data$x3p,levels = names(x3pList)))

    plts <- surfaceMat_df %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = .data$value))  +
      ggplot2::scale_fill_gradientn(colours = height.colors,
                                    values = scales::rescale(quantile(surfaceMat_df$value,height.quantiles,na.rm = TRUE)),
                                    breaks = function(lims){
                                      dat <- quantile(surfaceMat_df$value,legend.quantiles,na.rm = TRUE)

                                      dat <- dat %>%
                                        setNames(paste0(names(dat)," [",round(dat,3),"]"))

                                      return(dat)
                                    },
                                    na.value = na.value) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()) +
      ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(2.5,"in"),
                                                      label.theme = ggplot2::element_text(size = 8),
                                                      title.theme = ggplot2::element_text(size = 10),
                                                      frame.colour = "black",
                                                      ticks.colour = "black"),
                      colour = 'none') +
      ggplot2::labs(fill = expression("Rel. Height ["*mu*"m]")) +
      ggplot2::facet_wrap(~ x3p)
  }
  else if(type == "list"){
    plts <- purrr::pmap(.l = list(x3pList,
                                  names(x3pList)),
                        function(x3p,name){

                          surfaceMat_df <- x3p %>%
                            x3ptools::x3p_to_df() %>%
                            #perform some transformations on the
                            #x,y values so that the plot is
                            #representative of the actual surface
                            #matrix (i.e., element [1,1] of the
                            #surface matrix is in the top-left
                            #corner)
                            dplyr::mutate(xnew = max(.data$y) - .data$y,
                                          ynew = max(.data$x) - .data$x,
                                          value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                            dplyr::select(-c(.data$x,.data$y)) %>%
                            dplyr::rename(x=.data$xnew,
                                          y=.data$ynew) %>%
                            dplyr::mutate(value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                            dplyr::mutate(x3p = rep(name,times = nrow(.)))

                          ret <- surfaceMat_df %>%
                            ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
                            ggplot2::geom_raster(ggplot2::aes(fill = .data$value))  +
                            ggplot2::scale_fill_gradientn(colours = height.colors,
                                                          values = scales::rescale(quantile(surfaceMat_df$value,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                                                          breaks = function(lims){
                                                            dat <- quantile(surfaceMat_df$value,legend.quantiles,na.rm = TRUE)

                                                            dat <- dat %>%
                                                              setNames(paste0(names(dat),"\n[",round(dat,1),"]"))

                                                            return(dat)
                                                          },
                                                          na.value = na.value) +
                            ggplot2::theme_minimal() +
                            ggplot2::coord_fixed(expand = FALSE) +
                            ggplot2::theme(
                              axis.title.x = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_blank(),
                              axis.ticks.x = ggplot2::element_blank(),
                              axis.title.y = ggplot2::element_blank(),
                              axis.text.y = ggplot2::element_blank(),
                              axis.ticks.y = ggplot2::element_blank(),
                              panel.grid.major = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank(),
                              panel.background = ggplot2::element_blank(),
                              plot.title = ggplot2::element_text(hjust = .5,
                                                                 size = 9),
                              legend.position = "bottom") +
                            ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(.1,"in"),
                                                                            barwidth = grid::unit(.75,"in"),
                                                                            label.theme = ggplot2::element_text(size = 6),
                                                                            title.theme = ggplot2::element_text(size = 8),
                                                                            frame.colour = "black",
                                                                            ticks.colour = "black"),
                                            colour =  'none') +
                            ggplot2::labs(fill = expression("Rel. Height\n[micron]")) +
                            ggplot2::labs(title = name)

                          # browser()

                          retPlt <- ggplot2::ggplot_build(ret)

                          pltLegend <- cowplot::get_legend(retPlt$plot)

                          pltLegend <- cowplot::plot_grid(pltLegend$grobs[[1]])

                          ret <- ret + theme(legend.position = "none")

                          return(list(ret,pltLegend))
                        })
  }
  return(plts)
}