require(ggplot2)
require(stringr)

groupAverage <- function(mat, design) {
    if (length(unique(design$Group)) == dim(mat)[1]) return(NULL)
    mat <- mat[rowSums(mat) > 0,]
    mat <- cbind(data.frame(mat), group=design[,"Group"])
    mat <- data.frame(mat) %>% group_by(group) %>% summarise(across(everything(), list(mean)))
    mat <- data.frame(mat)
    rownames(mat) = mat[,'group']
    mat <- apply(mat[,!(names(mat) %in% c('group'))], c(1), as.numeric)
    design <- design[!duplicated(design[,1]),]
    mat <- mat[,design$Group]
    return(list(mat=mat, design=design))
}

plotCorrelationMatrix <- function(mat, design, col, output_prefix)
{
    set.seed(10)
    ha <- getAnnotation(mat, design, col, pos='col')

    mat <- t(mat)
    # Slow
    # png(paste0("heatmap_", output_prefix, ".png"), width=1200, height=1200)
    # draw(Heatmap(data.frame(mat), name="log10 x+1 or fc", column_title="Samples", row_title="Genes", top_annotation=ha,
    #              show_row_names=(dim(mat)[1] < 100), col=col$score))
    # dev.off()

    # df <- melt(cbind(mat,group=design$Group), id.vars='group')
    # pdf(paste0("violinplot_", output_prefix, ".pdf"))
    # g <- ggplot(df, aes(x=Var1, y=value, fill=group))+geom_violin()+theme_half_open()
    # plot(g)
    # dev.off()

    for (method in c('spearman', 'pearson')) {
        cormat <- round(cor(as.matrix(mat), method=method),4)
        for (i in 1:dim(cormat)[1]) {
            cormat[i,i] = NA
        }
        melted_cormat <- melt(cormat)
        melted_cormat$Var1 <- factor(melted_cormat$Var1, levels=colnames(mat))
        melted_cormat$Var2 <- factor(melted_cormat$Var2, levels=colnames(mat))
        pdf(paste0('corr_', output_prefix, '_', method, '.pdf'))
        g <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
        geom_tile()+scale_fill_viridis()+theme(axis.text.x = element_text(angle = 90))
        plot(g)
        dev.off()
    }
}

plotClustering <- function(mat, output_prefix)
{
    rownames(mat) = paste0(rownames(mat), '_', 1:dim(mat)[1])
    d <- dist(mat, method = "euclidean")
    # Hierarchical clustering using Complete Linkage
    hc1 <- hclust(d, method = "ward.D2")
    dend <- as.dendrogram(hc1)
    # Plot the obtained dendrogram
    pdf(paste0("hclust_", output_prefix, ".pdf"))
    plot(dend, cex = 0.6, hang = -1)
    dev.off()
}

#' @import GGally
applyCorrelationPlot <- function(mat, design, col, output_prefix) {
    # GGally slow
    # set.seed(42)    
    # cols = sample(1:dim(mat)[2], 100, replace = FALSE, prob = NULL)    
    # print(data.frame(mat[,cols]))
    # pdf(paste0('pairs_', output_prefix, '.pdf'))
    # ggpairs(data.frame(t(mat[,cols])), title="correlogram with ggpairs()")
    # dev.off()

    for (i in 1:(dim(mat)[1]-1)) {
        for (j in (i+1):dim(mat)[1]) {
            if (design[i,1] == design[j,1]) {
                png(paste0('pairs_', output_prefix, '_', i, '_', j, '.png'))
                g <- ggplot(data.frame(x=mat[i,], y=mat[j,]), aes(x=x, y=y))+geom_point(alpha=0.5)+theme_classic()
                plot(g)
                dev.off()
            }
        }
    }
}

analyzeCountMatrix <- function(matrix_file, design_file, output_prefix, col, mat=NULL, design=NULL, log10_flag=1) {
    if (is.null(mat)) {
        mat <- read.table(matrix_file, header=T, sep=" ")
    }
    if (is.null(design)) {
        design <- read.table(design_file, header=F, sep=" ")
        rownames(design) <- colnames(mat)
    }
    colnames(design) = c("Group", "Genot", "Cond")
    design <- design[rownames(design) %in% colnames(mat),]
    mat[is.na(mat)] = 0
    mat <- mat[rowSums(mat) > 0,]
    if (log10_flag != 0) {
        if (log10_flag > 0) {
            mat <- log10(mat+1)
        } else {
            mat <- apply(mat, c(1, 2), function(x) { return(-log10(x)) })
            mat <- data.frame(mat)
        }
    }
    mat <- t(mat)
    for (method in c('', '_average')) {
        if (method != '') {
            result <- groupAverage(mat, design)
            if (is.null(result)) break
            mat <- t(result$mat)
            design <- result$design
            rownames(design) = colnames(mat)
        }
        plotClustering(mat, paste0(output_prefix, method))
        plotCorrelationMatrix(mat, design, col, paste0(output_prefix, method))
        applyCorrelationPlot(mat, design, col, paste0(output_prefix, method))
    }
}

getAnnotation <- function(mat, design, col, pos) {
    annot = matrix(nrow=dim(mat)[1], ncol=2)
    rownames(annot) = rownames(mat)
    colnames(annot) = names(col)[1:2]
    for (n in colnames(annot)) {
        annot[,n] = design[,n]
    }
    annot = data.frame(annot)
    ha <- HeatmapAnnotation(df=annot, col=col, which=pos)
    return(ha)
}

plotCorrelationBetweenStudies <- function(mat_a, design_a, mat_b, design_b, output_prefix, col)
{
    ha_a <- getAnnotation(mat_a, design_a, col, 'row')
    ha_b <- getAnnotation(mat_b, design_b, col, 'col')
    set.seed(10)
    for (method in c('spearman', 'pearson')) {
        cormat <- round(cor(as.matrix(mat_a), as.matrix(mat_b), method=method),4)
        pdf(paste0("heatmap_", output_prefix, '_', method, ".pdf"))
        draw(Heatmap(data.frame(mat), name="Correlation", column_title="Study A", row_title="Study B",
                     top_annotation=ha_b, left_annotation=ha_a, col=col$score))
        dev.off()
        melted_cormat <- melt(cormat)
        head(melted_cormat)
        pdf(paste0('corr_', output_prefix, '_', method, '.pdf'))
        g <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
        geom_tile()+scale_fill_viridis()+theme(axis.text.x = element_text(angle = 90))
        plot(g)
        dev.off()
    }
}

compareTwoStudies <- function(mat_a, design_a, mat_b, design_b, output_prefix, col) {
    common_genes <- intersect(colnames(mat_a), colnames(mat_b))
    mat_a <- mat_a[,common_genes]
    mat_b <- mat_B[,common_genes]
    plotCorrelationBetweenStudies(mat_a, design_a, mat_b, design_b, output_prefix, col)
}

plotPCA <- function(mat, design, output_prefix) {
    require(matrixStats)
    mat <- t(mat)
    mat <- mat[,colMeans(mat) != 0]
    mat[is.na(mat)] <- 0
    mat <- mat[,colVars(mat) > 1e-5]
    ir.pca <- prcomp(mat,
                    center = TRUE,
                    scale. = TRUE)
    aload <- abs(ir.pca$rotation)
    write.table(sweep(aload, 2, colSums(aload), "/"), file=paste(output_prefix, "_cont.txt", sep="_"), quote=F)
    write.table(ir.pca$rotation, file=paste(output_prefix, "_cont_dir.txt", sep="_"), quote=F)
    s <- summary(ir.pca)
    capture.output(s, file=paste(output_prefix, "_summary.txt", sep=""))
    aload <- NULL
    png(paste(output_prefix, "_variance.png", sep=""))
    plot(ir.pca, type = "l")
    dev.off()
    scores = as.data.frame(ir.pca$x)
    sample_condition = factor(design[,1], levels=unique(design[,1]))
    patient_condition = factor(design[,2], levels=unique(design[,2]))
    time_condition = factor(design[,3], levels=unique(design[,3]))
    for (i in 1:4) {
        print(c('PCA plot -- ', i))
        pdf(paste0(output_prefix, c("", "_mark", "_both", "_comp")[i], "_pca_manual_color.pdf"), width=7, height=6)
        g <- ggplot(data = scores, aes(x = PC1, y = PC2)) + theme_classic()+theme(plot.margin= unit(c(1, 2, 1, 1), "cm"))+
        # theme(panel.background = element_rect(fill = "white"))+
        geom_hline(yintercept = 0, colour = "gray65") +
        geom_vline(xintercept = 0, colour = "gray65") +
        ggtitle("PCA plot")
        if (i == 1) {
            g <- g+geom_point(aes(fill = sample_condition), size=3, pch=21)
            g <- g+geom_text(label=1:dim(design)[1], vjust="inward", hjust="inward")
        } else if (i == 2) {
            g <- g+geom_text(label = patient_condition, vjust="inward",hjust="inward")
            g <- g+geom_point(aes(fill = time_condition), size=3, pch=21)
        } else if (i == 3) {
            g <- g+geom_point(aes(fill = sample_condition), size=3, pch=21)
            g <- g+geom_text(label = patient_condition, vjust="inward",hjust="inward")
        } else {
            g <- g+geom_point(aes(fill = time_condition), size=3, pch=21)
            g <- g+geom_text(label = patient_condition, vjust="inward",hjust="inward")
        }
        plot(g)
        dev.off()
    }
    for (pc1 in 1:9) {
        for (pc2 in (pc1+1):10) {
            for (i in c(2, 3)) {
                if (pc2 > dim(scores)[2]) break
                pdf(paste0(output_prefix, c("", "_mark", "_both", "_comp")[i], "_", pc1, "_", pc2, ".pdf"), width=7, height=6)
                g <- ggplot(data = scores, mapping = aes(x = .data[[paste0('PC', pc1)]], y = .data[[paste0('PC', pc2)]])) + theme_classic()+theme(plot.margin= unit(c(1, 2, 1, 1), "cm"))+
                    geom_hline(yintercept = 0, colour = "gray65") +
                    geom_vline(xintercept = 0, colour = "gray65") +
                    ggtitle("PCA plot")
                if (i == 2) {
                    g <- g+geom_text(label = patient_condition, vjust="inward",hjust="inward")
                    g <- g+geom_point(aes(fill = time_condition), size=3, pch=21)
                } else if (i == 3) {
                    g <- g+geom_point(aes(fill = sample_condition), size=3, pch=21)
                    g <- g+geom_text(label = patient_condition, vjust="inward",hjust="inward")
                }
                plot(g)
                dev.off()            
            }
        }
    }
}

plotCountMatrix <- function(mat, design, output_prefix, log10_flag) {
    colnames(design) = c("Group", "Genot", "Cond")
    if (dim(mat)[2] > 1 && any(duplicated(design$Group))) {
        factor_index <- factor(design$Group, levels=unique(design$Group))
        ave_mat <- cbind(t(mat), index=factor_index) %>% group_by(index) %>% summarise_at(vars(col_to_aggregate), list(name = mean))
        df <- melt(t(as.matrix(ave_mat)))
    } else {
        colnames(mat) = design$Group
        df <- melt(as.matrix(mat))
    }
    colnames(df) = c('Gene', 'Sample', 'Count')
    df$Gene <- as.character(df$Gene)
    for (norm in c('_', '_log10')) {
        if (norm == '_log10') {
            if (log10_flag)
                df$Count <- log10(df$Count)
            else
                break
        }
        if (length(unique(df$Group)) > 1) {
            g <- ggplot(df, mapping=aes(x=Sample, y=Count, color=Gene, group=Gene))
            g <- g+geom_point(size=3)+geom_line(linetype='dotted')+scale_colour_viridis(discrete=TRUE)
        } else {
            g <- ggplot(df, mapping=aes(x=Gene, y=Count, fill=Gene))
            g <- g+geom_bar(stat='identity')+scale_fill_viridis(discrete=TRUE)
        }
        g <- g+theme_minimal()
        pdf(paste0('raw_exp_comp_', output_prefix, norm, '.pdf'))
        plot(g)
        dev.off()
    }
}

genesetComparison <- function(mat, design, output_prefix, dir, col, log10_flag, gene_style='Refseq', verbose=TRUE) {
    pattern = 'geneset_.*.txt'
    gene_file_list <- c(list.files(dir, pattern=pattern), 'all')
    for (fname in gene_file_list) {
        if (fname != 'all') {
            gene_vector <- read.table(file.path(dir, fname), header=T, stringsAsFactors=F)[,1]
            header <- gsub('\\.txt', '', gsub('geneset_', '', fname))
            tmat <- cbind(rownames(mat), mat)
            colnames(tmat)[1] = 'Gene'
            subset_mat <- addGeneDescription(tmat, gene_style=gene_style)
            subset_mat <- subset_mat[toupper(subset_mat$Symbol) %in% toupper(gene_vector),]
            subset_mat <- subset_mat[!duplicated(subset_mat$Symbol),]
            rownames(subset_mat) <- subset_mat$Symbol
            subset_mat <- subset_mat[,!(colnames(subset_mat) %in% c('Refseq', 'Gene', 't2t', 'Symbol', 'Description')), drop=F]
        }
        if (verbose) {
            print(paste0('Gene set - ', output_prefix))
            print(dim(subset_mat))
        }
        if (dim(subset_mat)[1] < 20) {
            if (verbose) print(subset_mat)
            plotCountMatrix(subset_mat, design, paste(output_prefix, header, sep="_"), log10_flag)
        }
        if (min(dim(subset_mat)) > 2)
            analyzeCountMatrix('', '', paste(output_prefix, header, sep="_"), col, mat=subset_mat, design=design, log10_flag=log10_flag)
    }
}