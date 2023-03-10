
#' @import pathview

plotPathview <- function(map, fc_df, output_prefix, gene_style, cpd)
{
    for (plot_type in c('kegg_native', 'graphviz')) {
        if (plot_type == 'kegg_native') {
            kegg <- T
            layer <- F
            extension <- "png"
        } else {
            kegg <- F
            layer <- F
            extension <- "pdf"
        }
        if (gene_style == "Refseq" || nchar(gene_style) == 0) {
            pv.out <- pathview(gene.data = fc_df, gene.idtype = "REFSEQ", cpd.data = cpd, pathway.id = map, species = "mmu", out.suffix = map, keys.align = "y",
                kegg.native = kegg, match.data = T, key.pos = "topright", same.layer=layer)
        } else {
            pv.out <- pathview(gene.data = fc_df, gene.idtype = "ENSEMBL", cpd.data = cpd, pathway.id = map, species = "mmu", out.suffix = map, keys.align = "y",
                kegg.native = kegg, match.data = T, key.pos = "topright", same.layer=layer)
        }
        output_file <- paste(output_prefix, "_", map, ".", extension, sep="")
        system(paste("mv ", paste(paste0("mmu", map), map, extension, sep="."), output_file))
        system('sleep 0.01')
    }
}

plotKEGGPathway <- function(fc, output_prefix, gene_style, kegg_target_file, all=FALSE, cpd=NULL)
{
    id_list <- read.table(kegg_target_file, header=F, colClasses=c("character", "character"), sep=" ")[,2]
    fc_df <- data.frame(FC=fc)
    for (id in id_list) {
        plotPathview(id, fc_df, output_prefix, gene_style, cpd)
    }
}

#' @import clusterProfiler
computeKEGGEnrichment <- function(df, output_prefix, kegg_threshold=0.05)
{
    background <- unique(df$Entrez)
    for (change in c("up", "down", "change")) {
        if (change %in% c('up', 'down')) {
            index <- (df$SignificantChange == change)
        } else {
            index <- (df$SignificantChange != 'neutral')
        }
        file_name <- paste0(output_prefix, "_", change, ".txt")
        deg_fc <- df$Fc[index]
        names(deg_fc) <- df$Entrez[index]
        deg_fc <- deg_fc[order(abs(deg_fc), decreasing=TRUE)]
        append = F
        if (length(deg_fc) < 10) next
        for (percentage in seq(10, 100, 10)) {
            end <- as.integer(length(deg_fc)/100*percentage)
            if (end <= 1) next
            t <- enrichKEGG(as.character(as.integer(names(deg_fc)[1:end])), organism = "mmu", pvalueCutoff = kegg_threshold,
              pAdjustMethod = "BH", universe = names(deg_fc), minGSSize = 10, maxGSSize = 500, use_internal_data=TRUE, qvalueCutoff = 0.2)
            write(x=paste("#",percentage), file=file_name, append=append)
            append = T
            capture.output(as.data.frame(t), file=file_name, append=append)
            if (percentage == 100) {
                print(dim(t))
                if (!is.null(t) && dim(t)[1] >= 1) {
                    pdf(paste0(output_prefix, "_", change, "_barplot.pdf"))
                    plot(barplot(t, showCategory=20))
                    dev.off()
                }
            }
         }
    }
}


enrichmentAnalysisKEGG <- function(gene_list,
                                   entrez_gene,
                                   fc,
                                   fdr,
                                   threshold,
                                   output_prefix,
                                   cpd=NULL,
                                   gene_style='Refseq',
                                   kegg_target_file='kegg_target_id.txt') 
{
    df <- data.frame(Gene=gene_list, Entrez=entrez_gene, FDR=-log10(fdr), Fc=fc, 
                     Significance=sapply(fdr, getSignificanceLabels, threshold=threshold),
                     Change=sapply(fc, getFoldChangeLabels))    
    df <- cbind(df, SignificantChange=getSignificantChange(df$Significance, df$Change))
    df <- filterMostSignificant(df, 'FDR', 'Fc', 'Gene', logFDR=TRUE)

    # All genes
    fc <- df$Fc
    names(fc) =  df$Gene
    plotKEGGPathway(fc, paste0(output_prefix, '_all'), gene_style, kegg_target_file, all=TRUE, cpd=cpd)
    # DEG
    df$Fc[df$Significance == 'Non'] = 0
    fc <- df$Fc
    names(fc) =  df$Gene
    plotKEGGPathway(fc, output_prefix, gene_style, kegg_target_file, cpd=cpd)

    df <- filterMostSignificant(df, 'FDR', 'Fc', 'Entrez', logFDR=TRUE)
    df <- df[!is.na(df$Entrez),]
    computeKEGGEnrichment(df, output_prefix)
}


#' @import enrichR
enrichmentAnalysisEnrichR <- function(gene_list,
                                      fc,
                                      fdr,
                                      output_prefix,
                                      threshold) 
{
    palette <- getColorPalette()
    df <- data.frame(Gene=gene_list, FDR=-log10(fdr), Fc=fc, 
                     Significance=sapply(fdr, getSignificanceLabels, threshold=threshold),
                     Change=sapply(fc, getFoldChangeLabels))
    gene_style = 'Refseq'
    df <- addRefseqDescription(df, gene_style, keep_gene_list_column=TRUE)
    df <- cbind(df, SignificantChange=getSignificantChange(df$Significance, df$Change))
    setEnrichrSite("Enrichr") # human or mouse genes
    dbs_list <- c('KEGG_2021_Human', 'GTEx_Tissue_Expression_Up', 'GO_Biological_Process_2021')
    for (change in c('up', 'down', 'change')) {
        if (change %in% c('up', 'down')) {
            gene_set <- unique(df$Symbol[df$SignificantChange == change])
        } else {
            gene_set <- unique(df$Symbol[df$SignificantChange != 'neutral'])
        }
        if (length(gene_set) == 0) next
        enriched <- enrichr(gene_set, dbs_list)
        for (i in 1:length(dbs_list)) {
            db <- dbs_list[i]
            print(length(enriched[[db]]))
            print(head(enriched[[db]]))
            if (length(enriched[[db]]) == 0) next
            output_header <- paste0(output_prefix, '_', change, '_db', i)
            pdf(paste0(output_header, '_dotplot.pdf'), width=15, height=12)
            plot(plotEnrich(data.frame(enriched[[db]]), showTerms = 20, numChar = 75, y = "Count", orderBy = "P.value"))
            dev.off()
            write.table(enriched[[db]], paste0(output_header, '.tsv'))
        }
    }
}