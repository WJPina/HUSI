library(R.utils)
library(data.table)
library(RColorBrewer)

ID2symbol <- function(counts = NULL,database = 'genecode',version='70',dir = './',GRCh ='37') {
  if(database == 'genecode'){
    refpath = paste(dir,'genecode.id2sym.human.v',version,'.csv',sep = '')
    filename = paste('gencode.v',version,'.annotation.gtf.gz',sep = '')
    link = paste('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_',version,'/',filename,sep = '')
  }
  else if (database == 'esemble') {
    refpath = paste(dir,'esemble.id2sym.human.v',version,'.csv',sep = '')
    filename = paste('Homo_sapiens.GRCh',GRCh,'.',version,'.gtf.gz',sep = '')
    link = paste('https://ftp.ensembl.org/pub/release-',version,'/gtf/homo_sapiens/Homo_sapiens.GRCh',GRCh,'.',version,'.gtf.gz',sep = '')
  }
  if(!file.exists(refpath)){
    filepath = paste(dir,filename,sep = '')
    if(!file.exists(filepath) & !file.exists(gsub('.gz','',filepath))){
      download.file(link, destfile =filepath, method = "auto")
    }
    GTF =  fread(filepath)
    if (database == 'genecode'){
      idx=grepl('gene',GTF$V3)
      annotation = GTF[idx,'V9'][[1]]    
      id2sym_df = data.frame(symbol = sub(".*gene_name\\s(\\S+);.*", "\\1", annotation) %>% gsub('\\"','',.))
      rownames(id2sym_df) = sub(".*gene_id\\s(\\S+);.*", "\\1", annotation) %>% gsub('\\"','',.)
    }
    else if (database == 'esemble') {
      annotation = GTF$V9
      genes = sub('.*gene_name\\s"(\\S+)";.*', "\\1", annotation)
      idx=!duplicated(genes)
      id2sym_df = data.frame(symbol=genes[idx])
      rownames(id2sym_df) = sub('.*gene_id\\s"(\\S+)";.*', "\\1", annotation[idx])
    }
    write.csv(id2sym_df,file = refpath)
  }
  id2sym_df = read.csv(refpath,header = T,row.names = 1)
  counts = data.frame(counts)
  counts = counts[rownames(counts) %in% rownames(id2sym_df),]
  symbol = id2sym_df[rownames(counts),'symbol']
  counts = counts[!duplicated(symbol),]
  counts$Symbol = symbol[!duplicated(symbol)]
  return(counts)
}

Counts2TPM <- function(counts, gfe){
  idx=intersect(rownames(counts),rownames(gfe))
  counts=counts[idx,]
  symbol = counts$Symbol
  effLen=gfe[idx,'length']
  counts=counts[,-ncol(counts)]
  tpm=apply(counts, 2, function(c){  
    rate <- log(c) - log(effLen);
    denom <- log(sum(exp(rate)));
    t=exp(rate - denom + log(1e6));
    return(t)})
  rownames(tpm) <- symbol
  return(tpm)
}

getPPI_String <- function (object = NULL, species = 9606, score_threshold = 600, 
    save = FALSE) 
{
	require("data.table")
	require("igraph")
	
    linkFiles <- paste("https://stringdb-static.org/download/protein.links.v11.0/", 
        species, ".protein.links.v11.0.txt.gz", sep = "")
    if (!file.exists(sub(pattern = ".gz", replacement = "", x = basename(linkFiles)))) {
        if (!file.exists(basename(linkFiles))) 
            download.file(linkFiles, destfile = basename(linkFiles))
        gf <- gzfile(basename(linkFiles), "rt")
    }
    PPI <- read.table(gf, header = T, sep = "")
	PPI[,1] <- as.factor(PPI[,1])
	PPI[,2] <- as.factor(PPI[,2])
    close(gf)
    infoFiles <- paste("https://stringdb-static.org/download/protein.info.v11.0/", 
        species, ".protein.info.v11.0.txt.gz", sep = "")
    if (!file.exists(sub(pattern = ".gz", replacement = "", x = basename(infoFiles)))) {
        if (!file.exists(basename(infoFiles))) 
            download.file(infoFiles, destfile = basename(infoFiles))
        gf <- gzfile(basename(infoFiles), "rt")
    }
    Pinfo <- read.table(gf, header = T, sep = "\t", colClasses = c("character", 
        "character", "NULL", "NULL"), quote = "", row.names = 1)
    close(gf)
    PPI <- subset(PPI, combined_score > score_threshold)
    ENSP1 <- levels(PPI[, 1])
    levels(PPI[, 1]) <- toupper(Pinfo[ENSP1, ])
    ENSP2 <- levels(PPI[, 2])
    levels(PPI[, 2]) <- toupper(Pinfo[ENSP2, ])
    if (!is.null(object)) {
        gene_data <- rownames(object)
        gene_data_upper <- toupper(gene_data)
        gene_data <- as.data.frame(unique(as.data.table(data.frame(gene_data, 
            gene_data_upper)), by = "gene_data_upper"))
        rownames(gene_data) <- gene_data[, 2]
        PPI <- PPI[which(is.element(PPI[, 1], gene_data[, 2])), 
            ]
        PPI <- PPI[which(is.element(PPI[, 2], gene_data[, 2])), 
            ]
        levels(PPI[, 1]) <- gene_data[levels(PPI[, 1]), 1]
        levels(PPI[, 2]) <- gene_data[levels(PPI[, 2]), 1]
    }
    nodes <- union(PPI[, 1], PPI[, 2])
    links <- PPI[, 1:2]
    net <- graph_from_data_frame(d = links, vertices = nodes, 
        directed = FALSE)
    net <- igraph::simplify(net)
    if (save) {
        saveRDS(as_adj(net), paste(species, "_ppi_matrix_STRING-11.0.Rda", 
            sep = ""))
    }
    # file.remove(paste(species, ".protein.links.v11.0.txt.gz", 
    #     sep = ""))
    # file.remove(paste(species, ".protein.info.v11.0.txt.gz", 
    #     sep = ""))
    return(as_adj(net))
}

mygseaplot2 <- function (x, geneSetID, title = "", color = "green", base_size = 11, 
          rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
          ES_geom = "line",leg_x =  0.45,leg_y = 0.4) 
{
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer + theme(legend.position = c(leg_x, leg_y),
                                legend.title = element_blank(), legend.background = element_rect(fill = alpha('white', 0.8)),
                                title = element_text(size = 16),legend.text = element_text(size = 16))
  p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                 unit = "cm"))
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                           ymax = ~ymax, color = ~Description)) + xlab(NULL) + 
    ylab(NULL) + theme_classic(base_size) + theme(legend.position = "none", 
                                                  plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
                                                  axis.ticks = element_blank(), axis.text = element_blank(), 
                                                  axis.line.x = element_blank()) + scale_x_continuous(expand = c(0, 
                                                                                                                 0)) + scale_y_continuous(expand = c(0, 0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, 
                                                     "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") + 
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
                               l = 0.2, unit = "cm"))
  if (!is.null(title) && !is.na(title) && title != "") 
    p.res <- p.res + ggtitle(title)
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      # p.res <- p.res + theme(legend.position = "none")
      p.res <- p.res
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    rownames(pd) <- pd$Description
    pd <- pd[, -1]
    pd <- round(pd, 4)
    tp <- tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
                                                                         xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
                                                                                                                             0.95), ymin = quantile(p.res$data$runningScore, 
                                                                                                                                                    0.75), ymax = quantile(p.res$data$runningScore, 
                                                                                                                                                                           0.9))
  }
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
}

gsInfo <- function (object, geneSetID) 
{
  geneList <- object@geneList
  if (is.numeric(geneSetID)) 
    geneSetID <- object@result[geneSetID, "ID"]
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

gseaScores <- function (geneList, geneSet, exponent = 1, fortify = FALSE) 
{
  geneSet <- intersect(geneSet, names(geneList))
  N <- length(geneList)
  Nh <- length(geneSet)
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  Pmiss[!hits] <- 1/(N - Nh)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  }
  else {
    ES <- min.ES
  }
  df <- data.frame(x = seq_along(runningES), runningScore = runningES, 
                   position = as.integer(hits))
  if (fortify == TRUE) {
    return(df)
  }
  df$gene = names(geneList)
  res <- list(ES = ES, runningES = df)
  return(res)
}

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


library(CellChat)
netVisual_bubble_my<-function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
    pairLR.use = NULL, color.heatmap = c("Reds"), 
    n.colors = 9, direction = 1, thresh = 0.05, comparison = NULL, 
    group = NULL, remove.isolate = FALSE, max.dataset = NULL, 
    min.dataset = NULL, min.quantile = 0, max.quantile = 1, line.on = TRUE, 
    line.size = 0.2, color.text.use = TRUE, color.text = NULL, 
    title.name = NULL, font.size = 10, font.size.title = 10, 
    show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
    angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE) 
{
    # color.heatmap <- match.arg(color.heatmap)
    if (is.list(object@net[[1]])) {
        message("Comparing communications on a merged object \n")
    }
    else {
        message("Comparing communications on a single object \n")
    }
    if (is.null(vjust.x) | is.null(hjust.x)) {
        angle = c(0, 45, 90)
        hjust = c(0, 1, 1)
        vjust = c(0, 1, 0.5)
        vjust.x = vjust[angle == angle.x]
        hjust.x = hjust[angle == angle.x]
    }
    if (length(color.heatmap) == 1) {
        color.use <- tryCatch({
            RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
        }, error = function(e) {
            (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
        })
    }
    else {
        color.use <- color.heatmap
    }
    if (direction == -1) {
        color.use <- rev(color.use)
    }
    if (is.null(comparison)) {
        cells.level <- levels(object@idents)
        if (is.numeric(sources.use)) {
            sources.use <- cells.level[sources.use]
        }
        if (is.numeric(targets.use)) {
            targets.use <- cells.level[targets.use]
        }
        df.net <- subsetCommunication(object, slot.name = "net", 
            sources.use = sources.use, targets.use = targets.use, 
            signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
        df.net$source.target <- paste(df.net$source, df.net$target, 
            sep = " -> ")
        source.target <- paste(rep(sources.use, each = length(targets.use)), 
            targets.use, sep = " -> ")
        source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
        if (length(source.target.isolate) > 0) {
            df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                ncol = ncol(df.net)))
            colnames(df.net.isolate) <- colnames(df.net)
            df.net.isolate$source.target <- source.target.isolate
            df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
            df.net.isolate$pval <- 1
            a <- stringr::str_split(df.net.isolate$source.target, 
                " -> ", simplify = T)
            df.net.isolate$source <- as.character(a[, 1])
            df.net.isolate$target <- as.character(a[, 2])
            df.net <- rbind(df.net, df.net.isolate)
        }
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
        idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 
            0)
        if (sum(idx1) > 0) {
            values.assign <- seq(max(df.net$prob, na.rm = T) * 
                1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
            position <- sort(prob.original[idx1], index.return = TRUE)$ix
            df.net$prob[idx1] <- values.assign[match(1:length(idx1), 
                position)]
        }
        df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
            unique(df.net$source)])
        df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
            unique(df.net$target)])
        group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
            levels(df.net$target), sep = " -> ")
        df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
        df.net <- with(df.net, df.net[order(interaction_name_2), 
            ])
        df.net$interaction_name_2 <- factor(df.net$interaction_name_2, 
            levels = unique(df.net$interaction_name_2))
        cells.order <- group.names
        df.net$source.target <- factor(df.net$source.target, 
            levels = cells.order)
        df <- df.net
    }
    else {
        dataset.name <- names(object@net)
        df.net.all <- subsetCommunication(object, slot.name = "net", 
            sources.use = sources.use, targets.use = targets.use, 
            signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
        df.all <- data.frame()
        for (ii in 1:length(comparison)) {
            cells.level <- levels(object@idents[[comparison[ii]]])
            if (is.numeric(sources.use)) {
                sources.use <- cells.level[sources.use]
            }
            if (is.numeric(targets.use)) {
                targets.use <- cells.level[targets.use]
            }
            df.net <- df.net.all[[comparison[ii]]]
            df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
            df.net$source.target <- paste(df.net$source, df.net$target, 
                sep = " -> ")
            source.target <- paste(rep(sources.use, each = length(targets.use)), 
                targets.use, sep = " -> ")
            source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
            if (length(source.target.isolate) > 0) {
                df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                  ncol = ncol(df.net)))
                colnames(df.net.isolate) <- colnames(df.net)
                df.net.isolate$source.target <- source.target.isolate
                df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
                df.net.isolate$pval <- 1
                a <- stringr::str_split(df.net.isolate$source.target, 
                  " -> ", simplify = T)
                df.net.isolate$source <- as.character(a[, 1])
                df.net.isolate$target <- as.character(a[, 2])
                df.net <- rbind(df.net, df.net.isolate)
            }
            df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                unique(df.net$source)])
            df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                unique(df.net$target)])
            group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                levels(df.net$target), sep = " -> ")
            group.names0 <- group.names
            group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], 
                ")")
            if (nrow(df.net) > 0) {
                df.net$pval[df.net$pval > 0.05] = 1
                df.net$pval[df.net$pval > 0.01 & df.net$pval <= 
                  0.05] = 2
                df.net$pval[df.net$pval <= 0.01] = 3
                df.net$prob[df.net$prob == 0] <- NA
                df.net$prob.original <- df.net$prob
                df.net$prob <- -1/log(df.net$prob)
            }
            else {
                df.net <- as.data.frame(matrix(NA, nrow = length(group.names), 
                  ncol = 5))
                colnames(df.net) <- c("interaction_name_2", "source.target", 
                  "prob", "pval", "prob.original")
                df.net$source.target <- group.names0
            }
            df.net$group.names <- as.character(df.net$source.target)
            df.net$source.target <- paste0(df.net$source.target, 
                " (", dataset.name[comparison[ii]], ")")
            df.net$dataset <- dataset.name[comparison[ii]]
            df.all <- rbind(df.all, df.net)
        }
        if (nrow(df.all) == 0) {
            stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
        }
        idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 
            0)
        if (sum(idx1) > 0) {
            values.assign <- seq(max(df.all$prob, na.rm = T) * 
                1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
            position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
            df.all$prob[idx1] <- values.assign[match(1:length(idx1), 
                position)]
        }
        df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
        df <- df.all
        df <- with(df, df[order(interaction_name_2), ])
        df$interaction_name_2 <- factor(df$interaction_name_2, 
            levels = unique(df$interaction_name_2))
        cells.order <- c()
        dataset.name.order <- c()
        for (i in 1:length(group.names0)) {
            for (j in 1:length(comparison)) {
                cells.order <- c(cells.order, paste0(group.names0[i], 
                  " (", dataset.name[comparison[j]], ")"))
                dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
            }
        }
        df$source.target <- factor(df$source.target, levels = cells.order)
    }
    min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
    max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
    df$prob[df$prob < min.cutoff] <- min.cutoff
    df$prob[df$prob > max.cutoff] <- max.cutoff
    if (remove.isolate) {
        df <- df[!is.na(df$prob), ]
        line.on <- FALSE
    }
    if (!is.null(max.dataset)) {
        signaling <- as.character(unique(df$interaction_name_2))
        for (i in signaling) {
            df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
            cell <- as.character(unique(df.i$group.names))
            for (j in cell) {
                df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
                values <- df.i.j$prob
                idx.max <- which(values == max(values, na.rm = T))
                idx.min <- which(values == min(values, na.rm = T))
                dataset.na <- c(df.i.j$dataset[is.na(values)], 
                  setdiff(dataset.name[comparison], df.i.j$dataset))
                if (length(idx.max) > 0) {
                  if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
                    df.i.j$prob <- NA
                  }
                  else if ((idx.max != idx.min) & !is.null(min.dataset)) {
                    if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
                      df.i.j$prob <- NA
                    }
                    else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% 
                      dataset.na)) > 0) {
                      df.i.j$prob <- NA
                    }
                  }
                }
                df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
            }
            df[df$interaction_name_2 == i, "prob"] <- df.i$prob
        }
    }
    if (remove.isolate) {
        df <- df[!is.na(df$prob), ]
        line.on <- FALSE
    }
    if (nrow(df) == 0) {
        stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
    df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target), 
        unique(df$source.target)))
    g <- ggplot(df, aes(x = source.target, y = interaction_name_2, 
        color = prob, size = pval)) + geom_point(pch = 19) + 
        theme_linedraw() + theme(panel.grid.major = element_blank()) + 
        theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, 
            vjust = vjust.x), axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
    values <- c(1, 2, 3)
    names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
    g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), 
        breaks = sort(unique(df$pval)), labels = names(values)[values %in% 
            sort(unique(df$pval))], name = "p-value")
    if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
        g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
            na.value = "white", limits = c(quantile(df$prob, 
                0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
            breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 
                1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5, 
            title = "Commun. Prob."))
    }
    else {
        g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
            na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5, 
            title = "Commun. Prob."))
    }
    g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
        theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
    if (grid.on) {
        if (length(unique(df$source.target)) > 1) {
            g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 
                0.5, 1), lwd = 0.1, colour = color.grid)
        }
        if (length(unique(df$interaction_name_2)) > 1) {
            g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                0.5, 1), lwd = 0.1, colour = color.grid)
        }
    }
    if (!is.null(title.name)) {
        g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
    }
    if (!is.null(comparison)) {
        if (line.on) {
            xintercept = seq(0.5 + length(dataset.name[comparison]), 
                length(group.names0) * length(dataset.name[comparison]), 
                by = length(dataset.name[comparison]))
            g <- g + geom_vline(xintercept = xintercept, linetype = "dashed", 
                color = "grey60", size = line.size)
        }
        if (color.text.use) {
            if (is.null(group)) {
                group <- 1:length(comparison)
                names(group) <- dataset.name[comparison]
            }
            if (is.null(color.text)) {
                color <- ggPalette(length(unique(group)))
            }
            else {
                color <- color.text
            }
            names(color) <- names(group[!duplicated(group)])
            color <- color[group]
            dataset.name.order <- levels(df$source.target)
            dataset.name.order <- stringr::str_match(dataset.name.order, 
                "\\(.*\\)")
            dataset.name.order <- stringr::str_sub(dataset.name.order, 
                2, stringr::str_length(dataset.name.order) - 
                  1)
            xtick.color <- color[dataset.name.order]
            g <- g + theme(axis.text.x = element_text(colour = xtick.color))
        }
    }
    if (!show.legend) {
        g <- g + theme(legend.position = "none")
    }
    if (return.data) {
        return(list(communication = df, gg.obj = g))
    }
    else {
        return(g)
    }
}