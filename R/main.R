
goToGeneList <- function(go, contents, file, ref)
{
    gfile <- paste0(file, "_", paste(go, contents, sep="_"), "_gene_list.txt")
    if (!file.exists(gfile)) {
        geneList <- lookUp(go)
        print(c("Get", head(geneList), head(go)))
        if (length(geneList) == 0) return(c())
        write.table(cbind(geneList, unlist(lapply(geneList, function(x) { return(ref[ref[,2] == x, 13][1]) }))), file=gfile)
    } else {
        geneList <- as.character(read.table(gfile)[,1])
    }
    for (sp in c("WT", "FL")) {
        plotFC(paste0(file, "_", go, "_", sp), paste0("fold_", sp, "_"), geneList)
    }
    return(geneList)
}


.filterFoldChange <- function(mat) {
    mat[mat == "inf" | mat == "Inf"] <- 0 #20
    mat[mat == "-inf" | mat == "-Inf"] <- 0 #-20
    mat[mat == "Nan" | mat == 'nan'] <- 0
    return(mat)
}

.filterFDR <- function(mat) {
    mat[mat == "inf" | mat == "Inf"] <- 1 #20
    mat[mat == "-inf" | mat == "-Inf"] <- 1 #-20
    mat[mat == "Nan" | mat == 'nan'] <- 1
    return(mat)
}

pathwayEnrichmentAnalysis <- function(gene_list,
                                entrez_gene,
                                fold_change,
                                fdr,
                                threshold,
                                output_prefix,
                                type="KEGG",
                                gene_style="Refseq")
{
    index <- (!is.na(fold_change))
    fold_change <- .filterFoldChange(as.numeric(fold_change[index]))
    fdr <- .filterFDR(as.numeric(fdr[index]))
    gene_list <- as.character(gene_list[index])
    entrez_gene <- as.character(entrez_gene[index])
    if (type == 'KEGG') {
        enrichmentAnalysisKEGG(gene_list,
                               entrez_gene,
                               fold_change,
                               fdr,
                               threshold,
                               paste0('kegg_', output_prefix),
                               gene_style='Refseq',
                               kegg_target_file='kegg_target_id.txt')
    } else if (type == 'GO') {
        return()
        # enrichmentAnalysisGO(gene_list, fold_change, fdr, paste0('go_', output_prefix))
    } else if (type == 'EnrichR') {
        enrichmentAnalysisEnrichR(gene_list,
                                  fold_change,
                                  fdr,
                                  paste0('enrichr_', output_prefix),
                                  threshold)
    } else {
        return()
    }
}

#' Compute fold change and p-value using edgeR
#'
#' @param case Case factor
#' @param control Control factor
#' @param group_total Total group factors
#' @param fit Fitted edgeR model
#' @param output_file Output file name
#' @return
#'
enrichmentAnalysis <- function(input_file='filelist.txt',
                               design_file='design.txt',
                               output_dir='./',
                               verbose=TRUE)
{
    output_prefix <- "fold_change"
    output_suffix <- ".txt"
    output_file <- "mat.txt"
    if (!file.exists(output_file)) {
        combineCountMatrix(input_file='filelist.txt',
                           output_file='mat.txt',
                           style='featureCounts',
                           order=c(),
                           verbose=TRUE)
    }
    pattern <- paste0('^', output_prefix, '_[0-9].*\\', output_suffix)
    if (length(list.files('./', pattern=pattern)) == 0) {
        normalizeExpression(input_file='mat.txt',
                           design_file='design.txt',
                           normalization_method='edgeR',
                           rep=6,
                           verbose=TRUE)
    }
    fc_file_list <- list.files('./', pattern=pattern)
    for (file_name in fc_file_list) {
        if (verbose) {
            message(paste("-- Reading...", file_name))
        }
        mat <- read.table(file_name, header=T, sep=" ", stringsAsFactors=F)
        refseq_gene <- rownames(mat)
        entrez_gene <- ref2entr(refseq_gene)
        fc <- .filterFoldChange(mat[,'logFC'])
        cpm <- mat[,'logCPM']
        fdr <- .filterFDR(mat[,'FDR'])
        output_header <- gsub(output_suffix, '', file_name)
        plotVolcano(refseq_gene, fc, fdr, output_header, threshold=0.05)
        plotVolcanoGeneSet(refseq_gene, fc, fdr, output_header, c('LYSOSOME', 'AUTOPHAGY'), threshold=0.05)
        plotMDplot(cpm, fc, fdr, output_header, threshold=0.05)

        pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='KEGG')
        pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='GO')
        pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='EnrichR')
    }
    pattern <- paste0('^', output_prefix, '_log.*')
    fc_file_list <- list.files('./', pattern=pattern)
    design = read.table(design_file, header=F, sep=" ")
    colnames(design) = c('Group', 'Genot', 'Cond')
    for (file_name in fc_file_list) {
        mat <- read.table(file_name, header=T, sep="\t", stringsAsFactors=F)
        if (verbose) {
            message(file_name)
            message(colnames(mat))
        }
    }
}


enrichmentAnalysis()

