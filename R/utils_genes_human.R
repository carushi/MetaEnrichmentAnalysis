.storeDescriptionData <- function() {
    x <- read.table(file.path(HOMEDIR, 'ann/human_gene_description.txt'), header=F, stringsAsFactors=F, sep="\t", quote="\"")
    x <- x[,c(5, 10, 8, 9)]
    colnames(x) <- c("Ensembl", "t2t", "Symbol", "Description")
    human_gene_description <- x
    x <- read.table(file.path(HOMEDIR, 'ann/mouse_to_human_OMA_mapping.csv'), header=T, sep=",")
    symbol_table <- x[,c(1, 5)]
    colnames(symbol_table) <- c("mouse", "human")
    symbol_table[,1] <- toupper(symbol_table[,1])
    symbol_table[,2] <- toupper(symbol_table[,2])
    use_data(human_gene_description, symbol_table, internal=FALSE, overwrite=TRUE)
}

.addRefseqDescription <- function(gene_matrix, gene_style='Refseq', rowname_column='Gene', keep_gene_list_column=FALSE) {
    return(gene_matrix)
}

.addEnsemblDescription <- function(gene_matrix, gene_style='Ensembl', rowname_column='Gene', keep_gene_list_column=FALSE) {
    colnames(human_gene_description)[colnames(human_gene_description) == gene_style] = "Gene"
    print(head(gene_matrix))
    print(head(human_gene_description))
    merged <- merge(human_gene_description, gene_matrix, by='Gene', all.y=TRUE)
    rownames(merged) <- merged$rowname_column
    if (!keep_gene_list_column)
        merged <- merged[,colnames(merged) != rowname_column]
    return(merged)
}

ref2entr <- function(gene_list, gene_style = "") {
    gene_list <- as.character(gene_list)
    gene_list <- lapply(gene_list, .extractFirstGene)
    if (gene_style %in% c("t2t", "Ensembl")) {
        return(lapply(gene_list, function(x) {
            index = which(human_gene_description[,gene_style] == x)
            if (length(index) > 0) return(index[1])
            return(NA)
        }))
    } else if (gene_style == "Refseq" || gene_style == "") {
        stopifnot(FALSE)
    } else if (gene_style == "Entrez") {
        return(lapply(gene_list, function(x) {return(x)}))
    } else {
        stopifnot(FALSE)
    }
}


