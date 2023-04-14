.storeDescriptionData <- function() {
    x <- read.table(file.path(HOMEDIR, 'ann/refseq_description.txt'), header=F, colClasses=c("character", "character", "character"))
    colnames(x) <- c("Gene", "Symbol", "Description")
    refseq_description <- x
    x <- read.table(file.path(HOMEDIR, 'ann/mouse_to_human_OMA_mapping.csv'), header=T, sep=",")
    symbol_table <- x[,c(1, 5)]
    colnames(symbol_table) <- c("mouse", "human")
    symbol_table[,1] <- toupper(symbol_table[,1])
    symbol_table[,2] <- toupper(symbol_table[,2])
    use_data(refseq_description, symbol_table, internal=FALSE, overwrite=T)
}

.addRefseqDescription <- function(gene_matrix, gene_style='Refseq', rowname_column='Gene', keep_gene_list_column=FALSE) {
    stopifnot(gene_style == 'Refseq')
    colnames(gene_matrix)[1] = "Gene"
    merged <- merge(refseq_description, gene_matrix, by='Gene', all.y=TRUE)
    rownames(merged) <- merged$rowname_column
    if (!keep_gene_list_column)
        merged <- merged[,colnames(merged) != rowname_column]
    return(merged)
}

.addEnsemblDescription <- function(gene_matrix, gene_style='Refseq', rowname_column='Gene', keep_gene_list_column=FALSE) {
    return(gene_matrix)
}

#' @import org.Mm.eg.db

ref2entr <- function(gene_list, gene_style = "")
{
    gene_list <- as.character(gene_list)
    gene_list <- lapply(gene_list, .extractFirstGene)
    if (gene_style == "Ensembl") {
        x <- org.Mm.egENSEMBL2EG
    } else if (gene_style == "Refseq" || gene_style == "") {
        x <- org.Mm.egREFSEQ2EG
    } else if (gene_style == "Entrez") {
        return(lapply(gene_list, function(x) {return(x)}))
    } else {
        stopifnot(FALSE)
    }
    # Get the entrez gene identifiers that are mapped to any RefSeq ID
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- x[mapped_genes]
    xx <- as.list(xx)
    return(lapply(gene_list, function(x) {return(xx[[x]])}))
}