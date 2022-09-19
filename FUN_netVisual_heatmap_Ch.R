netVisual_heatmap_Ch <- function (object, comparison = c(1, 2), measure = c("count",
                                                    "weight"), signaling = NULL, slot.name = c("netP", "net"),
          color.use = NULL, color.heatmap = c("#2166ac", "#b2182b"),
          title.name = NULL, width = NULL, height = NULL, font.size = 8,
          font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
          sources.use = NULL, targets.use = NULL, remove.isolate = FALSE,
          row.show = NULL, col.show = NULL,
          MinSet = -20, MaxSet = 20, ylimSet = c(0,2))
{
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- (object@net[[comparison[1]]][[measure]]+1) %>% log2()
    obj2 <- (object@net[[comparison[2]]][[measure]]+1) %>% log2()
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    }
    else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  }
  else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[, , signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    }
    else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      }
      else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(MinSet, 0, MaxSet), # color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)),
                                   c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
                                                                      "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat,
                                                                                                                      na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
                                                                                                                                                     max(mat, na.rm = T))) + 1))
  }
  else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, MinSet, MaxSet),
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(MinSet, MaxSet),
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9,
                                                                                 name = color.heatmap))))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
                                                                      "\\1", min(mat, na.rm = T))) + 1), round(max(mat,
                                                                                                                   na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
                                                                                                                                                  max(mat, na.rm = T))) + 1))
  }
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)),ylim = ylimSet,
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)),
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)),ylim = ylimSet,
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)),
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white",
                name = legend.name, bottom_annotation = col_annotation,
                left_annotation = row_annotation, top_annotation = ha2,
                right_annotation = ha1, cluster_rows = cluster.rows,
                cluster_columns = cluster.rows, row_names_side = "left",
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size),
                column_names_gp = gpar(fontsize = font.size), column_title = title.name,
                column_title_gp = gpar(fontsize = font.size.title),
                column_names_rot = 90, row_title = "Sources (Sender)",
                row_title_gp = gpar(fontsize = font.size.title), row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8,
                                                            fontface = "plain"), title_position = "leftcenter-rot",
                                            border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8),
                                            grid_width = unit(2, "mm")))
  return(ht1)
}
