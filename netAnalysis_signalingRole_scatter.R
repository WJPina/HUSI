netAnalysis_signalingRole_scatter_log <- function (object, signaling = NULL, color.use = NULL, slot.name = "netP", 
  group = NULL, weight.MinMax = NULL, dot.size = c(2, 6), 
  point.shape = c(21, 22, 24, 23, 25, 8, 3), label.size = 3, 
  dot.alpha = 0.6, xlabel = "log10(Outgoing interaction strength)", 
  ylabel = "log10(Incoming interaction strength)", title = NULL, 
  font.size = 14, font.size.title = 16, do.label = T, show.legend = T, 
  show.axes = T) 
{
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (is.null(signaling)) {
    message("Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways")
  }
  else {
    message("Signaling role analysis on the cell-cell communication network from user's input")
    signaling <- signaling[signaling %in% object@netP$pathways]
    if (length(signaling) == 0) {
      stop("There is no significant communication for the input signaling. All the significant signaling are shown in `object@netP$pathways`")
    }
    outgoing <- outgoing[, signaling, drop = FALSE]
    incoming <- incoming[, signaling, drop = FALSE]
  }
  outgoing.cells <- rowSums(outgoing)
  incoming.cells <- rowSums(incoming)
  num.link <- aggregateNet(object, signaling = signaling, 
    return.object = FALSE, remove.isolate = FALSE)$count
  num.link <- rowSums(num.link) + colSums(num.link) - diag(num.link)
  df <- data.frame(x = log10(outgoing.cells), y = log10(incoming.cells), 
    labels = names(incoming.cells), Count = num.link)
  if (!is.null(group)) {
    df$Group <- group
  }
  if (is.null(color.use)) {
    color.use <- scPalette(nlevels(object@idents))
  }
  if (!is.null(group)) {
    gg <- ggplot(data = df, aes(x, y)) + geom_point(aes(size = Count, 
      colour = labels, fill = labels, shape = Group))
  }
  else {
    gg <- ggplot(data = df, aes(x, y)) + geom_point(aes(size = Count, 
      colour = labels, fill = labels))
  }
  gg <- gg + CellChat_theme_opts() + theme(text = element_text(size = font.size), 
    legend.key.height = grid::unit(0.15, "in")) + labs(title = title, 
    x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, 
    face = "plain")) + theme(axis.line.x = element_line(size = 0.25), 
    axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, 
    alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE) + 
    guides(colour = FALSE)
  if (!is.null(group)) {
    gg <- gg + scale_shape_manual(values = point.shape[1:length(unique(df$Group))])
  }
  if (is.null(weight.MinMax)) {
    gg <- gg + scale_size_continuous(range = dot.size)
  }
  else {
    gg <- gg + scale_size_continuous(limits = weight.MinMax, 
      range = dot.size)
  }
  if (do.label) {
    gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, 
      colour = labels), size = label.size, show.legend = F, 
      segment.size = 0.2, segment.alpha = 0.5)
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}
