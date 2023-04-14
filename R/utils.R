
getFoldChangeLabels <- function(x) {
    if (x > 0) return("up")
    if (x < 0) return("down")
    return("neutral")
}

getSignificanceLabels <- function(x, threshold) {
    if (x <= threshold) return("Significant")
    else return("Non")
}

getSignificantChange <- function(fdr, fc) {
    return(unlist(sapply(1:length(fc), function(x) {
        if (fdr[x] == "Significant") return(fc[x])
        else return('neutral')
    })))
}

filterMostSignificant <- function(df, fdr_column, fc_column, gene_column, logFDR=TRUE) {
    df <- df[order(df[,fdr_column], decreasing=logFDR),]
    df[,gene_column] <- as.character(df[,gene_column])
    unique_gene_list <- unique(df[duplicated(df[,gene_column]),gene_column])
    for (x in unique_gene_list[!is.na(unique_gene_list)]) {
        target_index <- which(df[,gene_column] == x)
        target_index <- unlist(target_index)
        target_index <- target_index[2:length(target_index)]
        df[target_index, fc_column] <- 0
        df[target_index, gene_column] <- NA
    }
    return(df)
}



up <- function(x) { if(x > 0) { return(abs(x))} else { return(0) } }
down <- function(x) { if (x < 0) { return(abs(x))} else { return(0) } }

':=' <- function(lhs, rhs) {
    frame <- parent.frame()
    lhs <- as.list(substitute(lhs))
    if (length(lhs) > 1)
    lhs <- lhs[-1]
    if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL))
    }
    if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
    if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
    for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
    return(invisible(NULL))
}

writeVolcanoTopGenes <- function(output_header, df, gene_style)
{
    # Top fold change
    end <- min(50, dim(df)[1])
    merged <- df[!is.na(df$Gene),]
    merged <- merged[order(abs(merged$Fc), decreasing=T), ]
    output_df <- merged[1:min(end, dim(merged)[1]),]
    merged <- merged[order(merged$FDR, decreasing=T),]
    output_df <- rbind(output_df, merged[1:end,])
    write.table(file=paste0(output_header, "_selected_gene.txt"), x=output_df, quote=T, row.names=T)
}

#' @import ggsci
getColorPalette <- function() {
    mypal = pal_aaas("default")(9)
    col = list("red"=mypal[2], "blue"=mypal[1], "black"="black", "grey"="grey")
    return(col)
}



#' @import ggplot2 cowplot
plotMDplot <- function(cpm, fc, fdr, output_header, threshold=0.05) {
    palette <- getColorPalette()
    df <- data.frame(CPM=cpm, FDR=-log10(fdr), Fc=fc, 
                     Significance=sapply(fdr, getSignificanceLabels, threshold=threshold),
                     Change=sapply(fc, getFoldChangeLabels))
    df <- cbind(df, SignificantChange=getSignificantChange(df$Significance, df$Change))
    df$SignificantChange <- factor(df$SignificantChange, levels=c('neutral', 'down', 'up'))
    df <- df[order(df$SignificantChange),]
    g <- ggplot(df, aes(x=CPM, y=Fc, colour=SignificantChange))
    g <- g+geom_point()
    g <- g+scale_color_manual(values=c("up"=palette[["red"]], "down"=palette[["blue"]], "neutral"=palette[["black"]]))
    g <- g+theme_half_open()
    g <- g+xlab("Log2 fold change")+ylab("FDR")
    pdf(paste0('mdplot_', output_header, '.pdf'))
    plot(g)
    dev.off()
}


#' @import ggplot2 cowplot
plotVolcano <- function(gene_list, fc, fdr, output_header, threshold=0.05, xrange=NULL, yrange=NULL, gene_style='Refseq') {
    palette <- getColorPalette()
    df <- data.frame(Gene=gene_list, FDR=-log10(fdr), Fc=fc, 
                     Significance=sapply(fdr, getSignificanceLabels, threshold=threshold),
                     Change=sapply(fc, getFoldChangeLabels))
    df <- addGeneDescription(df, gene_style, keep_gene_list_column=TRUE)
    df <- cbind(df, SignificantChange=getSignificantChange(df$Significance, df$Change))
    g <- ggplot(df, aes(x=Fc, y=FDR, color=SignificantChange, alpha=Significance))
    g <- g+geom_point()
    g <- g+scale_color_manual(values=c("up"=palette[["red"]], "down"=palette[["blue"]], "neutral"=palette[["black"]]))
    g <- g+scale_alpha_manual(values=c("Significant"=1.0, "Non"=0.3))
    g <- g+theme_half_open()
    g <- g+xlab("Log2 fold change")+ylab("FDR")
    g <- g+geom_hline(yintercept=-log10(threshold), linetype="dashed")
    pdf(paste0('volcano_', output_header, '.pdf'))
    plot(g)
    dev.off()
    if (any(colnames(df) == 'Symbol')) {
        df$Symbol[df$'Symbol' == '<NA>'] <- '-'
        g <- ggplot()
        g <- g+geom_point(data=df, aes(x=Fc, y=FDR, color=SignificantChange), color='grey')
        g <- g+geom_text(data=subset(df, df$Significance == 'Significant'), aes(x=Fc, y=FDR, label=Symbol, color=SignificantChange))
        g <- g+scale_color_manual(values=c("up"=palette[["red"]], "down"=palette[["blue"]], "neutral"="grey"))
        g <- g+theme_half_open()
        g <- g+xlab("Log2 fold change")+ylab("FDR")
        g <- g+geom_hline(yintercept=-log10(threshold), linetype="dashed")
        if (!is.null(yrange)) g <- g+ylim(yrange)
        if (!is.null(xrange)) g <- g+ylim(xrange)
        pdf(paste0('volcano_wt_', output_header, ".pdf"))
        plot(g)
        dev.off()
    }
    writeVolcanoTopGenes(paste0('top_gene_', output_header), df, gene_style)
}

.countChangeType <- function(df, data_set) {
    count_data <- data.frame(table(df$SignificantChange))
    count_data <- cbind(count_data, data=rep(data_set, dim(count_data)[1]))
    return(count_data)
}

plotVolcanoGeneSet <- function(gene_list, fc, fdr, output_header, gene_set_list, gene_style='Refseq', threshold=0.05, xrange=NULL, yrange=NULL, organism='mouse', verbose=TRUE) {
    palette <- getColorPalette()
    df <- data.frame(Gene=gene_list, FDR=-log10(fdr), Fc=fc, 
                     Significance=sapply(fdr, getSignificanceLabels, threshold=threshold),
                     Change=sapply(fc, getFoldChangeLabels))
    df <- addGeneDescription(df, gene_style, keep_gene_list_column=TRUE)
    df$Symbol <- sapply(df$Symbol, toupper)
    df <- cbind(df, SignificantChange=getSignificantChange(df$Significance, df$Change))
    for (gene_set_series in gene_set_list) {
        symbol_vectors <- extractGeneSet(gene_set_series, (organism=='mouse'))
        if (verbose) {
            print(c('---', gene_set_series))
            print(symbol_vectors)
        }
        for (gene_set_file in names(symbol_vectors)) {
            header <- gsub('\\.txt', '', gsub('geneset_', '', gene_set_file))
            g <- ggplot(df, aes(x=Fc, y=FDR, alpha=Significance))
            g <- g+geom_point(color='grey')
            g <- g+geom_point(data=subset(df, df$Symbol %in% symbol_vectors[[gene_set_file]]), aes(x=Fc, y=FDR, color=Change))
            g <- g+scale_color_manual(values=c("up"=palette[["red"]], "down"=palette[["blue"]], "neutral"="grey"))
            g <- g+scale_alpha_manual(values=c("Significant"=1.0, "Non"=0.3))
            g <- g+theme_half_open()
            g <- g+xlab("Log2 fold change")+ylab("FDR")
            g <- g+geom_hline(yintercept=-log10(threshold), linetype="dashed")
            pdf(paste0('volcano_', output_header, "_", gene_set_series, "_", header, '.pdf'))
            plot(g)
            dev.off()            
            g <- ggplot()
            g <- g+geom_point(data=df, aes(x=Fc, y=FDR, color=SignificantChange), color='grey')
            g <- g+geom_text(data=subset(df, df$Symbol %in% symbol_vectors[[gene_set_file]]), aes(x=Fc, y=FDR, label=Symbol, color=Change))
            g <- g+scale_color_manual(values=c("up"=palette[["red"]], "down"=palette[["blue"]], "neutral"="grey"))
            g <- g+theme_half_open()
            g <- g+xlab("Log2 fold change")+ylab("FDR")
            g <- g+geom_hline(yintercept=-log10(threshold), linetype="dashed")
            if (!is.null(yrange)) g <- g+ylim(yrange)
            if (!is.null(xrange)) g <- g+ylim(xrange)
            pdf(paste0('volcano_wt_', output_header, "_", gene_set_series, "_", header, ".pdf"))
            plot(g)
            dev.off()
            stats_table <- .countChangeType(df, 'all')
            stats_table <- rbind(stats_table, .countChangeType(df[df$Symbol %in% symbol_vectors[[gene_set_file]],], gene_set_series))
            print(stats_table)
            write.table(stats_table, paste0('sigchange_', output_header, '_', gene_set_series, '_', header, '.txt'))
            stats_table <- stats_table[stats_table[,'Var1'] != 'neutral',]
            pdf(paste0('barplot_', output_header, "_", gene_set_series, "_", header, ".pdf"))
            g <- ggplot(stats_table, aes(x = data, y = Freq, fill = Var1)) +
                geom_col(colour = "black", position="fill")
            g <- g+scale_fill_manual(values=c("up"=palette[["red"]], "down"=palette[["blue"]], "neutral"="grey"))
            g <- g+theme_classic()
            plot(g)
            dev.off()
        }
    }
}
