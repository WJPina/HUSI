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
    file.remove(paste(species, ".protein.links.v11.0.txt.gz", 
        sep = ""))
    file.remove(paste(species, ".protein.info.v11.0.txt.gz", 
        sep = ""))
    return(as_adj(net))
}