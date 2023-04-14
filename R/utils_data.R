#' Compute fold change and p-value using edgeR
#'
#' @param case Case factor
#' @param control Control factor
#' @param group_total Total group factors
#' @param fit Fitted edgeR model
#' @param output_file Output file name
#'
.computeFoldChangeByedgeR <- function(case, control, group_total, fit, output_file, header='')
{
    cont <- rep(0, group_total)
    cont[control] <- -1
    cont[case] <- 1
    lrt <- glmLRT(fit, contrast=cont)
    p <- topTags(lrt, n=200000) # Display all genes
    if (nchar(header) > 0) {
        cat(paste0(header, "\n"),file=output_file)
    }
    write.table(file=output_file, p, row.names=T, append=TRUE)
}

.getFirstLine <- function(filename) {
    con <- file(filename, "r")
    first_line <- readLines(con, n=1)
    close(con)
    first_line <- gsub(' ', '_', gsub('# ', '', first_line))
    return(first_line)
}

#' Combine selected DEG matrix data
#'
#' @param design Design vector to annotate samples
#' @param output_prefix Character string for output
#' @param output_suffix Character string for output
#' @examples
#'
.combineDEGMatrix <- function(design,
                             output_prefix,
                             output_suffix) {
    group <- factor(design, levels=unique(design))
    file_list <- list.files('./', paste0('^', output_prefix, '_[0-9].*\\', output_suffix))
    file_list <- file_list[order(file_list)]
    for (i in 1:2) {
        value <- colnames(read.table(file_list[1], header=T, sep=" ", stringsAsFactors=F))[i]
        value2 <- "FDR"
        output_file <- paste0(output_prefix, '_', value, output_suffix)
        df <- NULL
        for (filename in file_list) {
            column_header <- .getFirstLine(filename)
            result <- read.table(filename, header=T, sep=" ", stringsAsFactors=F)
            temp_df <- cbind(cbind(rownames(result), result[,i]), result[,value2])
            colnames(temp_df) <- c('Refseq', paste0(column_header, '_', value), paste0(column_header, '_', value2))
            if (is.null(df)) df <- temp_df
            else df <- merge(df, temp_df, by='Refseq', all=TRUE)
        }
        df <- addGeneDescription(df)
        write.table(df, file=output_file, sep='\t', quote=F)
    }
}

#' Apply dispersion estimation and fold change computation
#'
#' @param mat Data matrix
#' @param design Design vector to annotate samples
#' @param output_prefix Charcter string for output (default=fold_change)
#' @param output_suffix Character string for output (default=.txt)
#' @param writeIntegrated Boolean flag to determine whether the integrated matrix is written
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCCommonDisp estimateGLMTrendedDisp estimateGLMTagwiseDisp
#'
.normalizeCountByedgeR <- function(mat,
                                   design,
                                   output_prefix='fold_change',
                                   output_suffix='.txt',
                                   cpm_filtering_row=2,
                                   writeIntegrated=FALSE)
{
    sizeFactors <- colSums(mat)
    group <- factor(design, levels=unique(design))
    y <- DGEList(counts=mat,group=group)
    keep <- rowSums(cpm(y) > 1) >= cpm_filtering_row #Cpm filtering
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    design_matrix <- model.matrix(~0+group, data=y$samples)
    colnames(design_matrix) <- levels(y$samples$group)
    y <- estimateGLMCommonDisp(y,design_matrix)
    y <- estimateGLMTrendedDisp(y,design_matrix)
    y <- estimateGLMTagwiseDisp(y,design_matrix)
    fit <- glmFit(y, design_matrix)
    group_total <- length(unique(group))
    for (control in 1:(group_total-1)) {
        for (case in (control+1):group_total) {
            .computeFoldChangeByedgeR(case, control, group_total, fit,
                                      paste0(output_prefix, '_', control, '_', case, output_suffix),
                                      header=paste("#", colnames(design_matrix)[case], 'vs', colnames(design_matrix)[control]))
        }
    }
    if (writeIntegrated) {
        .combineDEGMatrix(design, output_prefix, output_suffix)
    }
}

#' Read fold changes from DE results
#'
#' @param input_file Input DEG file
#' @param style DEG format (default=edgeR)
#' @return Fold change and pvalue matrix
#'
.readFoldChange <- function(input_file, style='edgeR')
{
    mat <- read.table(input_file, header=T, stringsAsFactors=FALSE)
    if (style == 'edgeR') {
        return(cbind(rownames(mat), mat[,c(1, 5)]))
    } else if (style == 'DESeq2') {
        return(cbind(rownames(mat), mat[,c(2, 6)]))
    } else if (style == 'Cufflinks') {
        return(mat[,c(1, 4, 5)])
    } else {
        stopifnot(FALSE)
    }
}


#' Read CPM data from DE results
#'
#' @param input_file Input DEG file
#' @param style DEG format (default=edgeR)
#' @return Gene ID and CPM matrix
#'
readCpm <- function(input_file, style='edgeR')
{
    mat <- read.table(input_file, header=T)
    if (style == 'edgeR') {
        return(cbind(rownames(mat), mat[,c(2)]))
    } else if (style == 'DESeq2') {
        return(cbind(rownames(mat), mat[,c(1)]))
    } else if (style == 'Cufflinks') {
        # No CPM data is available
        return(NULL)
    } else {
        stopifnot(FALSE)
    }
}

.getColumnIndex <- function(style, enhancer_flag=FALSE) {
    if (style == 'featureCounts') {
        return(list(aidx=1, cidx=7))
    } else if (style == 'STAR') {
        if (enhancer_flag) { # unstranded
            return(list(aidx=1, cidx=2))
        } else {
            return(list(aidx=1, cidx=4))
        }
    } else if (style == 'bed') {
        return(list(aidx=8, cidx=11))        
    }
}

.readCountMatrix <- function(file, style, idx_list=NULL) {
    if (style == 'featureCounts') {
        df <- read.table(file, header=T)
    } else if (style == 'STAR') {
        df <- read.table(file, header=F)
        df <- df[!(df[,idx_list$aidx] %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")),]
    } else if (style == 'bed') {
        df <- read.table(file, header=F, sep="\t")
        colnames(df) <- c("chr", "start", "end", "activity", "gchr", "gstart", "gend", "Ensembl", "Symbol", "distance", "count")
    }
    return(df)
}

#' Convert count data for differential gene expression analysis
#'
#' @param input_file The file name which contains the file list to read
#' @param output_file The file name to store the computed count matrix
#' @param style Input file format (default=featureCounts)
#' @param order Reorder column indices of output matrix
#' @param verbose Flag for debugging
#' @return Count matrix data
#' @examples
#' ReadCountMatrix('filelist.txt', 'mat.txt')
#' ReadCountMatrix('filelist.txt', 'mat.txt', order= c(1:3, 7:9, 4:6, 10:12))
#' @export
#'

combineCountMatrix <- function(input_file='filelist.txt',
                               output_file='mat.txt',
                               style='featureCounts',
                               order=c(),
                               enhancer_flag=FALSE,
                               verbose=TRUE)
{
    idx_list = .getColumnIndex(style, enhancer_flag)
    file_list = unlist(read.table(input_file, stringsAsFactors=F, header=F))
    mat <- NULL
    amat <- NULL
    for (file in file_list) {
        df <- .readCountMatrix(file, style, idx_list)
        if (is.null(mat)) {
            mat <- as.matrix(df[,idx_list$cidx])
        } else {
            mat <- cbind(mat, df[,idx_list$cidx])
        }
        if (enhancer_flag) {
            rownames(mat) <- paste0('region_', 1:dim(mat)[1])
            amat <- df[,-c(idx_list$cidx)]
            rownames(amat) <- rownames(mat)
        } else {
            rownames(mat) <- df[,idx_list$aidx]
        }
    }
    colnames(mat) <- gsub('\\..*$', '', sapply(file_list, basename))
    if (!is.null(amat)) {
        amat <- amat[rowSums(mat) > 0, ]
    }
    mat <- mat[rowSums(mat) > 0, ]
    if (length(order) > 0) {
        mat <- mat[,order]
    }
    write.table(mat, file=output_file)
    if (!is.null(amat)) {
        write.table(amat, file=gsub('.txt', '_annotation.txt', output_file))
    }
    if (verbose) {
        print(head(mat))
    }
}

#' Noramlize expression data for DEG computation
#'
#' DESCRIPTION
#' @param input_file_list The temperature in degrees Fahrenheit
#' @return Count matrix data
#' @examples
#' temp1 <- F_to_C(50);
#' temp2 <- F_to_C( c(50, 63, 23) );
#' @export
#'

normalizeExpression <- function(input_file='mat.txt',
                                design_file='design.txt',
                                normalization_method='edgeR',
                                output_prefix='fold_change',
                                output_suffix='.txt',
                                ordered=c(),
                                rep=3,
                                cpm_filtering_row=2,
                                verbose=TRUE)
{
    mat <- read.table(input_file, header=T, stringsAsFactors = F)
    if (is.null(design_file)) {
        design <- as.integer(0:((dim(mat)[2]-1)/replicates)+1)
    } else {
        design <- read.table(design_file, header=F, sep=" ", row.names=NULL)[,1]
    }
    rownames(mat) <- as.character(rownames(mat))
    mat <- mat[,!is.na(design)]
    design <- design[!is.na(design)]
    if (verbose) {
        print(head(mat[,1:min(dim(mat)[2], 5)]))
        print(head(design))
    }
    .normalizeCountByedgeR(mat, design, output_prefix=output_prefix, output_suffix=output_suffix, cpm_filtering_row=cpm_filtering_row, writeIntegrated=TRUE)
}

.concatenateMetaboMatrix <- function(input_file, annotation_columns, used_columns) {
    mat <- NULL
    amat <- NULL
    for (file in unlist(read.table(input_file, header=F))) {
        df <- read.table(file, header=T, sep=",", stringsAsFactors=F, quote="\"")
        adf <- df[,annotation_columns]
        df <- df[,used_columns]
        rownames(df) = adf[,2]
        df[df == "N.D."] <- 0
        if (is.null(mat)) {
            mat <- df
            amat <- adf
        } else {
            mat <- rbind(mat, df)
            amat <- rbind(amat, adf)
        }
    }
    rownames(mat) <- as.character(rownames(mat))
    return(list(data=mat, annotation=amat))
}

.normalizeCountByScaling <- function(mat,
                                     annotation,
                                     design,
                                     output_prefix='fold_change',
                                     output_suffix='.txt',
                                     writeIntegrated=FALSE)
{
    sizeFactors <- colSums(mat)
    print(head(mat))
    print(head(design))

    # group <- factor(design, levels=unique(design))
    # y <- DGEList(counts=mat,group=group)
    # keep <- rowSums(cpm(y) > 1) >= 2 #Cpm filtering
    # y <- y[keep, , keep.lib.sizes=FALSE]
    # y <- calcNormFactors(y)
    # design_matrix <- model.matrix(~0+group, data=y$samples)
    # colnames(design_matrix) <- levels(y$samples$group)
    # y <- estimateGLMCommonDisp(y,design_matrix)
    # y <- estimateGLMTrendedDisp(y,design_matrix)
    # y <- estimateGLMTagwiseDisp(y,design_matrix)
    # fit <- glmFit(y, design_matrix)
    # group_total <- length(unique(group))
    # for (control in 1:(group_total-1)) {
    #     for (case in (control+1):group_total) {
    #         .computeFoldChangeByedgeR(case, control, group_total, fit,
    #                                   paste(output_prefix, '_', control, '_', case, output_suffix, sep=""),
    #                                   header=paste("#", colnames(design_matrix)[case], 'vs', colnames(design_matrix)[control]))
    #     }
    # }
    if (writeIntegrated) {
        print('write')
        # .combineDEGMatrix(design, output_prefix, output_suffix)
    }
}

normalizeMetabolome <- function(input_file='filelist.txt',
                                 design_file='design_mat.txt',
                                 output_file='mat_metabo.txt',
                                 output_prefix='mfold_change',
                                 output_suffix='.txt',
                                 normalization_method='pca',
                                 annotation_columns=1:2,
                                 used_columns=6:31,
                                 order_levels=NULL,
                                 verbose=TRUE)

{
    result <- .concatenateMetaboMatrix(input_file, annotation_columns, used_columns)
    write.table(cbind(result$annotation, result$data), output_file)
    design <- read.table(design_file, header=F, sep=" ", row.names=NULL)[,1]
    print(head(design))
    mat <- result$data
    mat <- mat[,!is.na(design)]
    design <- design[!is.na(design)]
    if (!is.null(order_levels)) {
        mat <- mat[,order(factor(design, levels=order_levels))]
        print(order(factor(design, levels=order_levels)))
        design <- design[order(factor(design, levels=order_levels)),]
    }
    if (verbose) {
        print(head(mat[,1:min(dim(mat)[2], 5)]))
        print(head(design))
    }
    .normalizeCountByScaling(mat, result$annotation, design, output_prefix=output_prefix, output_suffix=output_suffix, writeIntegrated=TRUE)
}