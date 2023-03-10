.extractFirstGene <- function(x) {
    if (x == "-" || nchar(x) == 0) return(NA)
    else if (length(grep(",", x)) == 0) return(x)
    return(substring(x, 1, regexpr(",", x)-1))
}

extractGeneSet <- function(dir) {
    file_path <- list()
    file_path[["LYSOSOME"]] = '~/project/GTP_RNA_seq/aws/Liver_data/ann/gene_set/lysosome/'
    file_path[["AUTOPHAGY"]] = '~/project/GTP_RNA_seq/aws/Liver_data/ann/gene_set/autophagy/'
    pattern = 'geneset_.*.txt'
    gene_file_list <- c(list.files(file_path[[dir]], pattern=pattern), 'all')
    gene_list <- list()
    for (fname in gene_file_list) {
        if (fname == 'all') {
            gene_list[[fname]] <- unique(unlist(gene_list))
        } else {
            gene_list[[fname]] <- read.table(file.path(file_path[[dir]], fname), header=T, stringsAsFactors=F)[,1]
        }
        gene_list[[fname]] <- unlist(sapply(gene_list[[fname]], function(x) {return(unlist(strsplit(toupper(x), ';')))}))
    }
    return(gene_list)
}


.storeRefseqData <- function() {
    x <- read.table('~/project/GTP_RNA_seq/aws/Liver_data/ann/refseq_description.txt', header=F, colClasses=c("character", "character", "character"))
    colnames(x) <- c("Gene", "Symbol", "Description")
    refseq_description = x
    use_data(refseq_description, internal=TRUE, overwrite=T)
}

addRefseqDescription <- function(gene_matrix, gene_style='Refseq', rowname_column='Gene', keep_gene_list_column=FALSE) {
    stopifnot(gene_style == 'Refseq')
    colnames(gene_matrix)[1] = "Gene"
    merged <- merge(refseq_description, gene_matrix, by='Gene', all.y=TRUE)
    rownames(merged) <- merged$rowname_column
    if (!keep_gene_list_column)
        merged <- merged[,colnames(merged) != rowname_column]
    return(merged)
}


extractMappedFoldChange <- function(fold_change, gene_list) {
    selected_change <- as.numeric(fold_change)
    names(selected_change) <- as.character(gene_list)
    return(selected_change[!(names(selected_change) == "NULL" || is.na(selected_change))])
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
