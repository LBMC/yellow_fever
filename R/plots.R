############################## palette function ################################
day_palette <- function(day){
  day_color <- list(D15 = "#ffcc03",
                    D90 = "#1565bd",
                    D136 = "#1565bd",
                    D593 = "#1a2634",
                    "In Vitro Restim" = "#1a2634",
                    "D100+" = "#1565bd")
  days_color <- lapply(as.list(day), FUN = function(x, day_color){
      day_color[[x]]
    }
    , day_color)
  days_color <- unlist(days_color)
  names(days_color) <- day
  return(days_color)
}

cell_type_palette <- function(cell_type){
  cell_type_color <- list(EM = "#c90000",
                          TEMRA = "#f4a582",
                          TSCM = "#92c5de",
                          CM = "#0571b0",
                          Naive = "gray",
                          NAIVE = "gray",
                          Effector = "#c90000",
                          EFF = "#c90000",
                          Memory = "#0571b0",
                          MEM = "#0571b0",
                          UNK = "gray")
  cell_types_color <- lapply(as.list(cell_type),
  FUN = function(x, cell_type_color){
      cell_type_color[[x]]
    }
    , cell_type_color)
  cell_types_color <- unlist(cell_types_color)
  names(cell_types_color) <- cell_type
  return(cell_types_color)
}

antigen_palette <- function(antigen){
  antigen_color <- list(A2 = "#d8b365",
                        B7 = "#5ab4ac")
  antigens_color <- lapply(as.list(antigen), FUN = function(x, antigen_color){
      antigen_color[[x]]
    }
    , antigen_color)
  antigens_color <- unlist(antigens_color)
  names(antigens_color) <- antigen
  return(antigens_color)
}

sex_palette <- function(sex){
  sex_color <- list(male = "#5ab4ac",
                    female = "#ee7593")
  sex_color <- lapply(as.list(sex), FUN = function(x, sex_color){
      sex_color[[x]]
    }
    , sex_color)
  sex_color <- unlist(sex_color)
  names(sex_color) <- sex
  return(sex_color)
}

LLC_palette <- function(LLC){
  LLC_color <- list(SLC = "#8dd3c7",
                    LLC = "#bebada")
  LLCs_color <- lapply(as.list(LLC), FUN = function(x, LLC_color){
      LLC_color[[x]]
    }
    , LLC_color)
  LLCs_color <- unlist(LLCs_color)
  names(LLCs_color) <- LLC
  return(LLCs_color)
}

cycling_palette <- function(cycling){
  cycling_color <- list(SLC = "white",
                    cycling = "black")
  cyclings_color <- lapply(as.list(cycling), FUN = function(x, cycling_color){
      cycling_color[[x]]
    }
    , cycling_color)
  cyclings_color <- unlist(cyclings_color)
  names(cyclings_color) <- cycling
  return(cyclings_color)
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ( (diff(h) %% 360) < 1) {
    h[2] <- h[2] - 360 / n
  }
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

clonality_palette <- function(clonality, other_set = TRUE){
  if (length(clonality) > 9 | other_set){
    clonality_color <- ggplotColours(length(clonality))
  } else {
    clonality_color <- brewer.pal(length(clonality), "Set1")
  }
  names(clonality_color) <- clonality
  clonality_color <- lapply(
    as.list(clonality),
    FUN = function(x, clonality_color){
      if (x %in% c("all", "other")){
        return("gray")
      }
      return(clonality_color[[x]])
    },
    clonality_color = clonality_color)
  clonality_color <- unlist(clonality_color)
  names(clonality_color) <- clonality
  return(clonality_color)
}

clonality_EFF_palette <- function(clonality, av_EM){
  clonality <- clonality[order(av_EM)]
  av_EM <- as.numeric(as.vector(av_EM[order(av_EM)]))
  r_select <- av_EM < 0.5
  clonality_EFF_color <- clonality
  names(clonality_EFF_color) <- clonality
  length_less <- length(which(r_select))
  length_more <- length(which(!r_select))
  if (length_less > 0){
    clonality_EFF_color[r_select] <- color_from_range(
      length_less,
      "#003298",
      "#79DCFF")
  }
  if (length_more > 0){
    clonality_EFF_color[!r_select] <- color_from_range(
      length_more,
      "#FFB8A5",
      "#A3080A")
  }
  return(clonality_EFF_color)
}

color_from_range <- function(size, start_col, stop_col){
  color_palette <- colorRampPalette(
    c(start_col, stop_col))(size)
  if (size == 1) {
    color_vector <- stop_col
  } else {
    color_vector <- colorRamp2(
      seq_len(size),
      color_palette)(seq_len(size))
  }
  return(color_vector)
}

#' @import circlize
count_scale_and_color <- function(counts,
  FUN = function(x){x<-ascb(x, to_zero = TRUE); x / max(x)}, quant=FALSE){
  counts <-apply(counts, 2, FUN)
  if (quant){
    color_ramp <- colorRampPalette(c("blue", "white", "red"),
                                        space = "Lab",
                                        interpolate = c("linear"))(50)
    colors <- colorRamp2(quantile(as.matrix(counts),
                                  seq(from=0,
                                      to=1,
                                      length.out=50),
                                  na.rm = TRUE),
                         color_ramp)
  }else{
    colors <- colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1),
                         c("white", "#fdd0bb", "#ff9472" ,"#f7553e", "#cc191d",
                         "#6b000a"))
    color_ramp <- colorRampPalette(c("white", "#780000"),
      space = "Lab",
      interpolate = c("spline"))(length(seq(from = 0, to = 1, by = 0.01)))
    colors <- colorRamp2(seq(from = 0, to = 1, length.out = 10),
                         c("white", RColorBrewer::brewer.pal(9, "Reds")))
  }
  return(list(counts = counts,
              colors = colors))
}


#' @import circlize
corr_scale_and_color <- function(counts, FUN=function(x){x}, quant=FALSE){
  counts <- apply(counts, 2, FUN)
  if (quant){
    colors <- colorRamp2(quantile(as.matrix(counts),
                                  seq(from=0,
                                      to=1,
                                      length.out=11),
                                  na.rm = TRUE),
                         rev(RColorBrewer::brewer.pal(11, "RdBu")))
    }else{
      counts <- counts * -1
      colors <- colorRamp2(c(1, 0.5, 0, -0.5, -1),
                           c("#09325d", "blue", "white", "red", "#6a000a"))
      if (min(counts[!is.na(counts)]) < 0){
        colors <- colorRamp2(seq(from = 1, to = -1, length.out = 11),
                             (RColorBrewer::brewer.pal(11, "RdBu")))
      }else{
        colors <- colorRamp2(seq(from = 0, to = 1, length.out = 11),
                             (RColorBrewer::brewer.pal(11, "RdBu")))
      }
    }
  return(list(counts = counts,
              colors = colors))
}


#' display single or multiple genes as violin plot
#'
#' @param scd is a scRNASeq data object
#' @param gene gene name to display
#' @param condition condition to compare (must be features of scd)
#' @param file output file where to save the graph (optional)
#' @param transform (bool) use anscomb transform or raw counts
#' @return a list() with the number of time each cells is classied as 'good' (
#' non is_blank looking) or 'bad' (blank looking)
#' @examples
#' \dontrun{
#' check_gene(scd, "genes_a", "sexe")
#' }
#' @import ggplot2
#' @export check_gene
check_gene <- function(scd, gene, condition, file, transform=TRUE){
  if (!is.null(ncol(scd$getfeature(condition)))){
    x <- data.frame(
      count = scd$getgene(gene),
      condition = scd$getfeature(condition)[, 1],
      condition2 = scd$getfeature(condition)[, 2]
    )
  } else {
    x <- data.frame(
      count = scd$getgene(gene),
      condition = scd$getfeature(condition)
    )
  }
  if (transform){
    x$count <- scRNAtools::ascb(x$count)
  }
  g <- ggplot2::ggplot(
    data = x,
    ggplot2::aes(x = condition, y = count, color = condition)
  )
  g <- g + ggplot2::geom_violin(
    alpha = 0.5,
    data = x[x$count != 0, ],
    ggplot2::aes(x = condition, y = ascb(count), color = condition)
  ) +
  ggplot2::geom_jitter(alpha = 0.5, height = 0) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = gene) +
  ggplot2::theme(legend.position = "none")
  if (!is.null(ncol(scd$getfeature(condition)))){
    g <- g + ggplot2::facet_wrap(~condition2)
  }
  print(g)
  if (!missing(file)) {
    ggsave(
      paste0(file, ".pdf"),
      width = 20, height = 15, units = "cm", dpi = 1200
    )
  }
}

#' display single or multiple genes as violin plot
#'
#' @param scd is a scRNASeq data object
#' @param gene gene name to display
#' @param condition condition to compare (must be features of scd)
#' @param file output file where to save the graph (optional)
#' @param transform (bool) use anscomb transform or raw counts
#' @return a list() with the number of time each cells is classied as 'good' (
#' non is_blank looking) or 'bad' (blank looking)
#' @examples
#' \dontrun{
#' check_gene(scd, "genes_a", "sexe")
#' }
#' @import ggplot2
#' @import ade4
#' @export pca_plot
pca_plot <- function(scd, color=NULL, shape=NULL, size=NULL, alpha=NULL,
  wrap=NULL, file, main = "", axes, res=FALSE, rainbow=FALSE, heatcolor=FALSE,
  is_contour, label=NULL, genes_list, arrow=FALSE, color_name="day",
  return_data = FALSE, color_scale = NULL, tmp_file) {
  if (!missing(tmp_file) & file.exists(tmp_file)) {
    print("tmp file found skipping pca...")
    load(tmp_file)
  } else {
    pca_out <- dudi.pca(scd$getcounts,
                        scan = F,
                        nf = scd$getncells - 1,
                        scale = TRUE,
                        center = TRUE)
    prop_of_var <- 100 * pca_out$eig[1:3] / sum(pca_out$eig)
    pca_data <- data.frame(x = pca_out$l1$RS1,
                           y = pca_out$l1$RS2,
                           var_x = round(prop_of_var[1], digit = 2),
                           var_y = round(prop_of_var[2], digit = 2),
                           axes = "axes 1 and 2")
    if (!is.null(color)){
      if (is.numeric(scd$getfeature(color))){
        if (mean(scd$getfeature(color)[pca_data$x <= mean(pca_data$x)]) > mean(scd$getfeature(color))){
          pca_data$x <- - pca_data$x
          pca_data$y <- - pca_data$y
        }
      }
    }
    if (!missing(tmp_file)){
      save(pca_out, prop_of_var, pca_data, file = tmp_file)
    }
  }
  x_lab <- paste0(
    "first axe: ", round(prop_of_var[1], digit = 2), "% of variance")
  y_lab <- paste0(
    "second axe: ", round(prop_of_var[2], digit = 2), "% of variance")
  if (missing(axes)){
    g <- plot_2_axes(
      scd,
      x = pca_data$x,
      y = pca_data$y,
      color = color,
      shape = shape,
      size = size,
      alpha = alpha,
      wrap = wrap,
      main = main,
      x_lab = x_lab,
      y_lab = y_lab,
      file = file,
      is_contour = is_contour,
      label = label,
      display = FALSE,
      color_name = color_name,
      color_scale = color_scale)
  }else{
    for (axe in axes){
      pca_data_sup <- data.frame(
        x = pca_out$l1[[paste0("RS", axe[1])]],
        y = pca_out$l1[[paste0("RS", axe[2])]],
        var_x = round(prop_of_var[axe[1]], digit = 2),
        var_y = round(prop_of_var[axe[2]], digit = 2),
        axes = paste0("axes ", axe[1], " and ", axe[2]))
      pca_data <- rbind(pca_data, pca_data_sup)
    }
    g <- plot_2_axes(
      scd,
      x = pca_data$x,
      y = pca_data$y,
      color = color,
      shape = shape,
      size = size,
      alpha = alpha,
      wrap = wrap,
      main = main,
      x_lab = x_lab,
      y_lab = y_lab,
      wrap = pca_data$axes,
      file = file,
      rainbow = rainbow,
      heatcolor = heatcolor,
      is_contour = is_contour,
      display = FALSE,
      color_name = color_name,
      color_scale = color_scale)
  }
  if (arrow | !missing(genes_list)){
    g <- scRNAtools::draw_arrow_dudi(g, pca_out, genes_list)
  }
  print(g)
  if (return_data){
    return(pca_out)
  }
}

draw_arrow_dudi <- function(g, dudi_obj, genes_list){
  genes_loading <- dudi_obj$c1[order(abs(dudi_obj$c1[, 1]), decreasing = T), ]
  genes_loading$genes <- rownames(genes_loading)
  genes_loading$ori <- 0
  if ("ls" %in% names(dudi_obj)){
    dudi_score <- sweep(
      dudi_obj$ls, 2, apply(dudi_obj$ls, 2, norm, type = "2"), "/")
  }else{
    dudi_score <- sweep(
      dudi_obj$li, 2, apply(dudi_obj$li, 2, norm, type = "2"), "/")
    names(dudi_score) <- paste0("CS", 1:length(dudi_score))
  }
  if (!missing(genes_list)){
    genes_loading <- genes_loading[genes_loading$genes %in% genes_list, ]
  }else{
    genes_loading <- rbind(genes_loading[1:min(nrow(genes_loading), 10), ],
      genes_loading[
        order(abs(dudi_obj$c1[, 2]),
          decreasing = T), ][1:min(nrow(genes_loading), 5), ])
  }
  genes_loading$CS1 <- (genes_loading$CS1 / max(abs(genes_loading$CS1))) *
    max(abs(dudi_obj$l1$RS1))
  genes_loading$CS2 <- (genes_loading$CS2 / max(abs(genes_loading$CS2))) *
    max(abs(dudi_obj$l1$RS2))
  genes_loading$shift <- genes_loading$CS2 / abs(genes_loading$CS2) * 0.02
  g_arrow <- g + geom_segment(
    data = genes_loading,
    aes(x = ori, y = ori, xend = CS1, yend = CS2, colour = NULL, alpha = NULL,
      shape = NULL, size = NULL),
      arrow = arrow(length = unit(0.3, "cm")))
  g_arrow <- g_arrow +
  geom_text(data = genes_loading,
    aes(x = CS1, y = CS2, label = genes, colour = NULL, alpha = NULL,
      shape = NULL, size = NULL),
    nudge_y = genes_loading$shift, check_overlap = TRUE)
  return(g_arrow)
}

#' @import grDevices
plot_2_axes <- function(scd, x, y, color = NULL, shape = NULL, size = NULL,
  alpha = NULL, wrap = NULL, main = "", x_lab = "", y_lab = "", file,
  rainbow = FALSE, heatcolor = FALSE, is_contour, label = NULL, display = TRUE,
  color_name = "day", color_scale = NULL) {
  data <- data.frame(x = x, y = y)
  aes_str <- "x=x, y=y"
  for (opt in c("color", "size", "alpha", "shape", "wrap", "label")){
    if (!is.null(get(opt))){
      data[[opt]] <- as.factor(scd$getfeature(get(opt)))
      if (opt != "wrap"){
        if (opt != "shape" | !is.numeric(shape)){
          aes_str <- paste0(aes_str, ", ", opt, "=", opt)
        }
      }
    }
  }
  if (!missing(is_contour)){
    g_cmd <- paste0("data = data[!is_contour,], aes(", aes_str, ")")
  }else{
    g_cmd <- paste0("data = data, aes(", aes_str, ")")
  }
  g <- eval(parse(text = paste0("ggplot(", g_cmd, ")")))
  if (!missing(is_contour)){
    g <- g +
    geom_density2d(data = data, aes(x = x, y = y, color = NULL, shape = NULL,
      size = NULL, alpha = NULL), colour = "gray")
  }
  if (!is.null(label)){
    g <- g +
    eval(parse(text = paste0(
      "geom_text(check_overlap = TRUE, fontface = \"bold\", ", g_cmd, ")")))
  }else{
    if (is.numeric(shape)){
      old_data <- data
      for (i in levels(factorize(shape))){
        r_select <- old_data$shape %in% i
        data <- old_data[r_select, ]
        g <- g + eval(parse(
          text = paste0("geom_point(", g_cmd, ", shape=", i, ")")))
      }
    }else{
      g <- g + eval(parse(text = paste0("geom_point(", g_cmd, ")")))
    }
  }
  if (!is.null(wrap)){
    g <- g + facet_wrap(~data$wrap)
  }
  if (color_name == "keep"){
    g <- g +
      scale_fill_manual(values = levels(as.factor(as.vector(data$color))))
    g <- g +
      scale_color_manual(values = levels(as.factor(as.vector(data$color))))
  }
  if (color_name == "keep_levels"){
    g <- g +
      scale_fill_manual(values = color_scale)
    g <- g +
      scale_color_manual(values = color_scale)
  }
  if (color_name == "day"){
    g <- g +
      scale_fill_manual(
        values = day_palette(levels(as.factor(as.vector(data$color))))
      )
    g <- g +
      scale_color_manual(
        values = day_palette(levels(as.factor(as.vector(data$color))))
      )
  }
  if (color_name == "p_EM"){
    g <- g +
      scale_fill_manual(
        values = cell_type_palette(levels(as.factor(as.vector(data$color))))
      )
    g <- g +
      scale_color_manual(
        values = cell_type_palette(levels(as.factor(as.vector(data$color))))
      )
  }
  if (color_name == "clonality"){
    g <- g +
      scale_fill_manual(
        values = clonality_palette(levels(as.factor(as.vector(data$color))))
      )
    g <- g +
    scale_color_manual(
      values = clonality_palette(levels(as.factor(as.vector(data$color))))
    )
  }
  if (color_name == "cell_type"){
    g <- g +
      scale_fill_manual(
        values = cell_type_palette(levels(as.factor(as.vector(data$color))))
      )
    g <- g +
    scale_color_manual(
      values = cell_type_palette(levels(as.factor(as.vector(data$color))))
    )
  }
  if (color_name == "sex"){
    g <- g +
      scale_fill_manual(
        values = sex_palette(levels(as.factor(as.vector(data$color))))
      )
    g <- g +
      scale_color_manual(
        values = sex_palette(levels(as.factor(as.vector(data$color))))
      )
  }
  if (is.numeric(color)){
    if (color_name == "cytotoxic_score"){
      g <- g + scale_colour_gradient(low = "#00ac00", high = "#ff0090")
      g <- g + scale_fill_gradient(low = "#00ac00", high = "#ff0090")
    }
    if (color_name == "p_EM"){
        g <- g +
          scale_colour_gradient2(low = "blue", mid = "gray", high = "red",
            midpoint = 0.5)
        g <- g +
          scale_fill_gradient2(low = "blue", mid = "gray", high = "red",
            midpoint = 0.5)
    }
  }
  if (rainbow){
    g <- g + scale_colour_gradientn(colours = rev(
        rainbow(20, s = 0.9, v = 0.9, start = .0, end = 0.65)))
  }
  if (heatcolor){
    colfunc <- colorRampPalette(c("red", "gray", "royalblue"))
    g <- g + scale_colour_gradientn(colours = colfunc(100))
  }
  g <- g +
       labs(title = main,
          x = x_lab,
          y = y_lab) +
       theme_bw()
  if (display) {
    print(g)
    if (!missing(file)) {
      ggsave(paste0(file, ".pdf"), width = 20, height = 15, units = "cm",
        dpi = 1200)
    }
  } else {
    return(g)
  }
}

#' plot the 2 first axis of a BCA analysis
#'
#' @param scd is a scRNASeq data object
#' @param gene gene name to display
#' @param condition condition to compare (must be features of scd)
#' @param file output file where to save the graph (optional)
#' @param transform (bool) use anscomb transform or raw counts
#' @return a list() with the number of time each cells is classied as 'good' (
#' non is_blank looking) or 'bad' (blank looking)
#' @examples
#' \dontrun{
#' check_gene(scd, "genes_a", "sexe")
#' }
#' @import ade4 ggplot2
#' @export bca_plot
bca_plot <- function(scd, by, color, ncomp=2, top, n_groups, n_cells, norm_by,
  file, main="", genes_list, xlimit, ylimit){
  by <- scRNAtools::factorize(by)
  data <- scd$copy()
  if (!missing(top)){
    genes_order <- order(
      abs(
        scRNAtools::bca_loading(ascb(scd$getcounts), by = by)$CS1
      ),
      decreasing = TRUE)
    data$order(genes = genes_order)
    data$order(genes = 1:min(top, ncol(data)))
  }
  if (!missing(n_groups)){
    by_size <- table(by)
    by_size <- by_size[order(by_size, decreasing = TRUE)]
    b_cells <- which(by %in% names(by_size)[1:min(n_groups, length(by_size))])
    data <- data$select(b_cells = b_cells)
    by <- scRNAtools::factorize(by[b_cells])
    if (!missing(norm_by)) {
      if(is.null(ncol(norm_by))){
        norm_by <- scRNAtools::factorize(norm_by[b_cells])
      } else {
        norm_by <- scRNAtools::factorize(norm_by[b_cells, ])
      }
    }
    if (!missing(color)){
      color <- color[b_cells]
    }
    print(dim(data))
  }
  if (missing(color)){
    color <- by
    circle_color <- by
  } else {
    color <- scRNAtools::factorize(color)
    circle_color <- rep(NA, length(by))
    for (i in levels(by)){
      circle_color[by %in% i] <- max.factor(
        scRNAtools::factorize(color[by %in% i]))
    }
  }
  if (!missing(norm_by)){
    data <- scRNAtools::wca_norm(data$getcounts, norm_by, ncomp)
  }
  pca_out <- dudi.pca(data$getcounts,
                      scan = F,
                      nf = data$getncells - 1)
  bca_out <- bca(pca_out,
                 by,
                 scan = F,
                 nf = ncomp)
  rbca_out <- randtest(bca_out)
  wca_out <- wca(pca_out,
                 by,
                 scan = F,
                 nf = ncomp)
  genes_loading <- bca_out$c1[order(abs(bca_out$c1[ ,1]), decreasing = T), ]
  bca_score <- sweep(bca_out$ls, 2, apply(bca_out$ls, 2, norm, type = "2"), "/")
  if (dim(bca_out$ls)[2] == 1) {
    dataToPlot <- data.frame(comp1 = bca_score$CS1, by = by)
    g <- ggplot(
      data = dataToPlot,
        aes(x = comp1, group = by, fill = by)) +
      geom_histogram() + ggtitle(main)
  } else {
    genes_loading$genes <- rownames(genes_loading)
    genes_loading$ori <- 0
    if (!missing(genes_list)){
      genes_loading <- genes_loading[genes_loading$genes %in% genes_list, ]
    } else {
      genes_loading <- rbind(genes_loading[1:min(nrow(genes_loading), 10), ],
        genes_loading[
          order(abs(bca_out$c1[, 2]), decreasing=T),
        ][1:min(nrow(genes_loading), 5), ])
    }
    genes_loading$CS1 <- genes_loading$CS1 *
      max(bca_score$CS1) / max(genes_loading$CS1)
    genes_loading$CS2 <- genes_loading$CS2 *
      max(bca_score$CS2) / max(genes_loading$CS2)
    genes_loading$shift <- genes_loading$CS2 / abs(genes_loading$CS2) * 0.02
    dataToPlot <- data.frame(
      comp1 = bca_score$CS1,
      comp2 = bca_score$CS2,
      by = by,
      color = color,
      circle_color = circle_color)
    dataToPlotby <- matrix(ncol = 5)
    for (i in levels(factorize(dataToPlot$by))) {
      r_select <- as.vector(dataToPlot$by) %in% i
      dataToPlotby <- rbind(
        dataToPlotby,
        c(
          mean(as.vector(dataToPlot$comp1)[r_select]),
          mean(as.vector(dataToPlot$comp2)[r_select]),
          i,
          as.vector(dataToPlot$color)[r_select][1],
          as.vector(dataToPlot$circle_color[r_select])[1]
        )
      )
    }
    colnames(dataToPlotby) <- colnames(dataToPlot)
    dataToPlotby <- as.data.frame(dataToPlotby[-1, ])
    dataToPlotby$comp1 <- vectorize(dataToPlotby$comp1)
    dataToPlotby$comp2 <- vectorize(dataToPlotby$comp2)
    g <- ggplot() +
      geom_segment(data = genes_loading,
        aes(x = ori, y = ori, xend = CS1, yend = CS2, group = NA, fill = NA,
          label = NA),
        arrow = arrow(length = unit(0.3, "cm")))
    g <- g + geom_text(data = genes_loading,
        aes(x = CS1, y = CS2, label = genes, group = NA, fill = NA),
        nudge_y = genes_loading$shift, check_overlap = TRUE)
    g <- g + stat_ellipse(data = dataToPlot,
        aes(x = comp1, y = comp2, group = by, fill = circle_color),
        geom = "polygon", alpha = 0.25, level = 0.8)
    g <- g + geom_point(data = dataToPlot,
      aes(x = comp1, y = comp2, color = color))
    g <- g + geom_text(data = dataToPlotby,
        aes(x = comp1, y = comp2, label = by),
        check_overlap = TRUE)
    g <- g + theme_bw()
    g <- g + labs(title = paste0(main,
      " between variability: ", round(bca_out$ratio, digits = 3),
      " within variability: ", round(wca_out$ratio, digits = 3),
      " for ", length(levels(by)), " groups"),
      x = "axe 1", y = "axe 2", fill = "cell-type")
  }
  if (!missing(xlimit)){
    if (length(xlimit) == 2){
      g <- g + xlim(xlimit[1], xlimit[2])
    }
  }
  if (!missing(ylimit)){
    if (length(ylimit) == 2){
      g <- g + ylim(ylimit[1], ylimit[2])
    }
  }

  print(g)
  if (!missing(file)) {
    print("saving file")
    ggsave(
      paste0(file, ".pdf"), width = 20, height = 15, units = "cm", dpi = 1200
    )
  }
}


#' plot the 2 first axis of a pCMF analysis
#'
#' @param scd is a scRNASeq data object
#' @param gene gene name to display
#' @param condition condition to compare (must be features of scd)
#' @param file output file where to save the graph (optional)
#' @param transform (bool) use anscomb transform or raw counts
#' @return a list() with the number of time each cells is classied as 'good' (
#' non is_blank looking) or 'bad' (blank looking)
#' @examples
#' \dontrun{
#' check_gene(scd, "genes_a", "sexe")
#' }
#' @import pCMF ggplot2
#' @export bca_plot
pCMF_plot <- function(scd, by, color, ncomp=2, top, n_groups, n_cells, norm_by,
  file, main="", genes_list, xlimit, ylimit, ncores = 4){
  if (!missing(tmp_file) & file.exists(tmp_file)) {
    print("tmp file found skipping pCMF...")
    load(tmp_file)
  } else {
    pCMF_out <- pCMF(
      X = scd$getcounts,
      ncomp = ncomp,
      iterMax = 500,
      iterMin = 100,
      epsilon = 1e-3,
      verbose = TRUE,
      sparse = TRUE,
      ZI = TRUE,
      ncores = ncores
    )
    U <- getU(pCMF_out)
    prop_of_var <- expDeviance(pCMF, scd$getcounts)
    pCMF_data <- data.frame(x = U[1, ],
                           y = U[2, ],
                           var_x = round(prop_of_var[1], digit = 2),
                           var_y = round(prop_of_var[2], digit = 2),
                           axes = "axes 1 and 2")
    if (!missing(tmp_file)){
      save(pCMF_out, pCMF_data, file = tmp_file)
    }
  }
  g <- plot_2_axes(
    scd,
    x = pCMF_data$x,
    y = pCMF_data$y,
    color = color,
    shape = shape,
    size = size,
    alpha = alpha,
    wrap = wrap,
    main = main,
    x_lab = paste0(
      "first axe: ", pCMF_data$var_x[1], "% of deviance"),
    y_lab = paste0(
      "first axe: ", pCMF_data$var_y[2], "% of deviance"),
    file = file,
    is_contour = is_contour,
    label = label,
    display = FALSE,
    color_name = color_name,
    color_scale = color_scale)
  print(g)
  if (return_data){
    return(pCMF_out)
  }
}
