library(CellChat)
netVisual_bubble_my<-function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
    pairLR.use = NULL, color.heatmap = c("OrRd"), 
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
        color = prob, size = pval)) + geom_point(pch = 16) + 
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