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
                                gene_style="Refseq",
                                organism="mouse")
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
                               gene_style=gene_style,
                               organism=organism,
                               kegg_target_file='kegg_target_id.txt')
    } else if (type == 'GO') {
        return()
        # enrichmentAnalysisGO(gene_list, fold_change, fdr, paste0('go_', output_prefix))
    } else if (type == 'EnrichR') {
        enrichmentAnalysisEnrichR(gene_list,
                                  fold_change,
                                  fdr,
                                  paste0('enrichr_', output_prefix),
                                  threshold,
                                  gene_style=gene_style)
    } else {
        return()
    }
}

analyzeSingleComparisonData <- function(file_name, output_suffix, gene_style, gene_set, organism, output_annotation_file='', gene_annotation_flag=TRUE) {
    mat <- read.table(file_name, header=T, sep=" ", stringsAsFactors=F)
    symbol <- NULL
    if (output_annotation_file == '') {
        refseq_gene <- rownames(mat)
        entrez_gene <- ref2entr(refseq_gene, gene_style=gene_style)
    } else {
        amat <- read.table(output_annotation_file, header=T, sep=" ", stringsAsFactors=F, row.names=1)
        amat <- amat[rownames(mat),]
        if (any(colnames(amat) == 'Ensembl')) {
            refseq_gene <- amat[,'Ensembl']
            entrez_gene <- ref2entr(refseq_gene, gene_style=gene_style)
        } else {
            refseq_gene <- amat[,'name']
            entrez_gene <- NULL
        }
        symbol <- amat[,"group3"]
    }
    fc <- .filterFoldChange(mat[,'logFC'])
    cpm <- mat[,'logCPM']
    fdr <- .filterFDR(mat[,'FDR'])
    output_header <- gsub(output_suffix, '', file_name)
    plotMDplot(cpm, fc, fdr, output_header, threshold=0.05)
    plotVolcano(refseq_gene, fc, fdr, output_header, threshold=0.05, gene_style=gene_style, symbol=symbol)
    if (gene_annotation_flag) {
        plotVolcanoGeneSet(refseq_gene, fc, fdr, output_header, gene_set, threshold=0.05, gene_style=gene_style, organism=organism, symbol=symbol)
        pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='EnrichR', gene_style=gene_style)
        pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='KEGG', gene_style=gene_style, organism=organism)
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
                               style='featureCounts',
                               organism='mouse',
                               gene_style='Refseq',
                               gene_set=c('LYSOSOME', 'AUTOPHAGY', 'LYSOTF'),
                               enhancer_flag=FALSE,
                               gene_annotation=TRUE,
                               verbose=TRUE)
{
    output_prefix <- "fold_change"
    output_suffix <- ".txt"
    output_file <- "mat.txt"
    if (enhancer_flag) {
        basename <- unlist(strsplit(input_file, '\\.'))[1]
        output_middle <- gsub('filelist_', '', basename, fixed=TRUE)
        output_prefix <- paste(output_prefix, output_middle, sep="_")
        output_file <- paste0("mat_", output_middle, ".txt")
        output_annotation_file <- gsub('.txt', '_annotation.txt', output_file)
    }
    if (!file.exists(output_file) || TRUE) {
        combineCountMatrix(input_file=input_file,
                           output_file=output_file,
                           style=style,
                           order=c(),
                           enhancer_flag=enhancer_flag,
                           verbose=TRUE)
    }
    pattern <- paste0('^', output_prefix, '_[0-9].*\\', output_suffix)
    if (length(list.files('./', pattern=pattern)) == 0 || TRUE) {
        normalizeExpression(input_file=output_file,
                            output_prefix=output_prefix,
                            output_suffix=output_suffix,
                            design_file='design.txt',
                            normalization_method='edgeR',
                            rep=6,
                            cpm_filtering_row=1,
                            verbose=TRUE)
    }
    fc_file_list <- list.files('./', pattern=pattern)
    for (file_name in fc_file_list) {
        if (verbose) {
            message(paste("-- Reading...", file_name))
        }
        analyzeSingleComparisonData(file_name, output_suffix, gene_style, gene_set, organism, output_annotation_file, gene_annotation)
    }
    # pattern <- paste0('^', output_prefix, '_log.*')
    # fc_file_list <- list.files('./', pattern=pattern)
    # design = read.table(design_file, header=F, sep=" ")
    # colnames(design) = c('Group', 'Genot', 'Cond')
    # for (file_name in fc_file_list) {
    #     mat <- read.table(file_name, header=T, sep="\t", stringsAsFactors=F)
    #     if (verbose) {
    #         message(file_name)
    #         message(colnames(mat))
    #     }
    # }
}

metabolomeAnalysis <- function(input_file='filelist.txt',
                               design_file='design.txt',
                               output_dir='./',
                               verbose=TRUE)
{
    output_prefix <- "mfold_change"
    output_suffix <- ".txt"
    output_file <- "mat_metabolome.txt"
    design_file <- "design_blood_metab.txt"
    # if (!file.exists(output_file)) {
        normalizeMetabolome(input_file="filelist_metabo.txt",
                           design_file=design_file,
                           output_file=output_file,
                           output_prefix=output_prefix,
                           output_suffix=output_suffix,
                           normalization_method='pca',
                           annotation_columns=1:2,
                           used_columns=6:31,
                           order_levels=c('WT', 'FL', 'WTF', 'FLF', 'WTH', 'FLH'),
                           verbose=TRUE)
    # }

    # fc_file_list <- list.files('./', pattern=pattern)
    # for (file_name in fc_file_list) {
    #     if (verbose) {
    #         message(paste("-- Reading...", file_name))
    #     }
    #     mat <- read.table(file_name, header=T, sep=" ", stringsAsFactors=F)
    #     refseq_gene <- rownames(mat)
    #     entrez_gene <- ref2entr(refseq_gene)
    #     fc <- .filterFoldChange(mat[,'logFC'])
    #     cpm <- mat[,'logCPM']
    #     fdr <- .filterFDR(mat[,'FDR'])
    #     output_header <- gsub(output_suffix, '', file_name)
    #     # plotVolcano(refseq_gene, fc, fdr, output_header, threshold=0.05)
    #     # plotVolcanoGeneSet(refseq_gene, fc, fdr, output_header, c('LYSOSOME', 'AUTOPHAGY'), threshold=0.05)
    #     # plotMDplot(cpm, fc, fdr, output_header, threshold=0.05)

    #     pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='KEGG')
    #     # pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='GO')
    #     # pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='EnrichR')
    # }

    # output_prefix <- "ffold_change"
    # output_suffix <- ".txt"
    # output_file <- "mat_fattyacid.txt"
    # if (!file.exists(output_file)) {
    #     normalizeExpression(input_file='mat.txt',
    #                        design_file='design.txt',
    #                        normalization_method='edgeR',
    #                        rep=6,
    #                        verbose=TRUE)
    # }

    # fc_file_list <- list.files('./', pattern=pattern)
    # for (file_name in fc_file_list) {
    #     if (verbose) {
    #         message(paste("-- Reading...", file_name))
    #     }
    #     mat <- read.table(file_name, header=T, sep=" ", stringsAsFactors=F)
    #     refseq_gene <- rownames(mat)
    #     entrez_gene <- ref2entr(refseq_gene)
    #     fc <- .filterFoldChange(mat[,'logFC'])
    #     cpm <- mat[,'logCPM']
    #     fdr <- .filterFDR(mat[,'FDR'])
    #     output_header <- gsub(output_suffix, '', file_name)
    #     # plotVolcano(refseq_gene, fc, fdr, output_header, threshold=0.05)
    #     # plotVolcanoGeneSet(refseq_gene, fc, fdr, output_header, c('LYSOSOME', 'AUTOPHAGY'), threshold=0.05)
    #     # plotMDplot(cpm, fc, fdr, output_header, threshold=0.05)

    #     pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='KEGG')
    #     # pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='GO')
    #     # pathwayEnrichmentAnalysis(refseq_gene, entrez_gene, fc, fdr, 0.05, output_header, type='EnrichR')
    # }    
}


enrichmentAnalysis()

