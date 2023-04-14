.extractFirstGene <- function(x) {
    if (x == "-" || nchar(x) == 0) return(NA)
    else if (length(grep(",", x)) == 0) return(x)
    return(substring(x, 1, regexpr(",", x)-1))
}

removeVersionInfo <- function(gene_list) {
    return(sapply(gene_list, function(x) {unlist(strsplit(as.character(x), '\\.'))[1]}))
}

.convertMouse2Human <- function(x) {
    # borrow from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://asia.ensembl.org")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://asia.ensembl.org")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    return(humanx)
}

mapMouse2Human <- function(x, mouse_origin=TRUE) {
    if (mouse_origin) {
        m <- merge(data.frame(mouse=x), symbol_table, all.x=TRUE, by='mouse')
        gene_set <- m$human
    } else {
        m <- merge(data.frame(human=x), symbol_table, all.x=TRUE, by='human')
        gene_set <- m$mouse
    }
    return(gene_set[!is.na(gene_set)])
}

extractGeneSet <- function(dir, conversion=TRUE) {
    file_path <- list()
    file_path[["LYSOSOME"]] = LYSOSOME
    file_path[["AUTOPHAGY"]] = AUTOPHAGY
    file_path[["LYSOTF"]] = LYSOTF
    file_path[["REPROGRAMMING"]] = REPROGRAMMING
    pattern = 'geneset_.*.txt'
    gene_file_list <- c(list.files(file_path[[dir]], pattern=pattern), 'all')
    gene_list <- list()
    for (fname in gene_file_list) {
        if (fname == 'all') {
            gene_list[[fname]] <- unique(unlist(gene_list))
        } else {
            gene_list[[fname]] <- read.table(file.path(file_path[[dir]], fname), header=T, stringsAsFactors=F)[,1]
            gene_list[[fname]] <- unlist(sapply(gene_list[[fname]], function(x) {return(unlist(strsplit(toupper(x), ';')))}))
            if (conversion) {
                gene_list[[fname]] <- mapMouse2Human(gene_list[[fname]], conversion)
            }
        }
    }
    return(gene_list)
}

extractMappedFoldChange <- function(fold_change, gene_list) {
    selected_change <- as.numeric(fold_change)
    names(selected_change) <- as.character(gene_list)
    return(selected_change[!(names(selected_change) == "NULL" || is.na(selected_change))])
}

addGeneDescription <- function(gene_matrix, gene_style='Refseq', rowname_column='Gene', keep_gene_list_column=FALSE) {
    if (gene_style == 'Refseq') {
        return(.addRefseqDescription(gene_matrix, gene_style, rowname_column, keep_gene_list_column))
    } else if (gene_style %in% c('t2t', 'Ensembl')) {
        return(.addEnsemblDescription(gene_matrix, gene_style, rowname_column, keep_gene_list_column))
    }
    stopifnot(FALSE)
}

obtainEnhancerName <- function(enhancer_id, gene_id) {
    return(sapply(1:length(enhancer_id), function(x) {
        if (!is.na(gene_id[x]))
            return(paste0(enhancer_id[x], '_', gene_id[x]))
        else
            return(enhancer_id[x])
    }))
}