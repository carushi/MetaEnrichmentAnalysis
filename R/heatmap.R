compareLiver <- function(threshold=0.05) {
    matrix_file = 'mat.txt'
    design_file = 'design.txt'
    mat = read.table(matrix_file, header=T, sep=" ")
    design = read.table(design_file, header=F, sep=" ")
    colnames(design) = c('Group', 'Genot', 'Cond')
    rownames(design) = colnames(mat[1:dim(mat)[2]])

    color_panel <- pal_aaas("default")(10)[c(10, 1, 2)]
    names(color_panel) = c('black', 'blue', 'red')
    col=list(Genot=c(WT=color_panel[['black']], FL=color_panel[['red']], BKO=color_panel[['blue']]), Cond=c(N=viridis(5)[5], F=viridis(5)[1], H="white"),
             score=viridis(256))
    output_prefix = 'liver'
    analyzeCountMatrix(matrix_file, design_file, output_prefix, col, log10_flag=1, mat=mat, design=design)
    genesetComparison(mat, design, paste0(output_prefix, '_lysosome'), col, dir=LYSOSOME, log10_flag=T)
    plotPCA(mat, design, output_prefix)

    # mat <- mat[,-c(28,29)]
    # design <- design[-c(28,29),]
    # output_prefix = 'liver_comp_drop'
    # analyzeCountMatrix(matrix_file, design_file, output_prefix, col, log10_flag=1, mat=mat, design=design)

    require("RColorBrewer")
    col=list(Genot=c(WT=color_panel[['black']], FL=color_panel[['red']], BKO=color_panel[['blue']]), Cond=c(N=viridis(5)[5], F=viridis(5)[1], H="white"),
             score=rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)))


    mat = read.table("fold_change_logFC.txt", header=T, sep="\t")
    mat = mat[,grep('_vs_WT_logFC', colnames(mat))]
    output_prefix = "liver_comp"
    design = read.table('design2.txt', header=F, sep=" ")
    design = design[1:5,]
    rownames(design) = colnames(mat)
    analyzeCountMatrix(matrix_file, 'design2.txt', output_prefix, col, mat=mat, log10_flag=F)
    plotPCA(mat, design, output_prefix)
    genesetComparison(mat, design, paste0(output_prefix, '_lysosome'), col, dir=LYSOSOME, log10_flag=0)

    for (tail in c('logFC', 'FDR')) {
        mat <- read.table("fold_change_logFC.txt", header=T, sep="\t")
        sample_list <- c('WTF_vs_WT', 'FLF_vs_FL', 'BKOF_vs_BKO')
        fc_mat <- mat[,paste0(sample_list, '_logFC')]
        if (tail == 'FDR') {
            pmat <- mat[,paste0(sample_list, '_FDR')]
            fc_mat[(pmat >= threshold)] = 0
            print(dim(fc_mat))
            fc_mat <- fc_mat[rowMins(as.matrix(fc_mat)) != rowMaxs(as.matrix(fc_mat)),]
            print(dim(fc_mat))
        }
        design = read.table('design2.txt', header=F, sep=" ")
        design = design[c(3, 8, 12),]
        rownames(design) = colnames(fc_mat)
        output_prefix = "liver_comp_cond"
        log10_flag = 0
        if (tail == 'FDR') {
            output_prefix = paste0(output_prefix, '_', tail)
        }
        analyzeCountMatrix(matrix_file, 'design2.txt', 'fasting_comp', col, mat=mat, log10_flag=log10_flag)
        genesetComparison(fc_mat, design, paste0(output_prefix, '_lysosome'), col, dir=LYSOSOME, log10_flag=log10_flag)
    }
}


compareBlood <- function(threshold=0.05) {
    matrix_file = 'mat.txt'
    design_file = 'design.txt'
    mat = read.table(matrix_file, header=T, sep=" ")
    design = read.table(design_file, header=F, sep=" ")
    colnames(design) = c('Group', 'Genot', 'Cond')
    rownames(design) = colnames(mat[1:dim(mat)[2]])
    color_panel <- pal_aaas("default")(10)[c(10, 1, 2)]
    names(color_panel) = c('black', 'blue', 'red')
    col=list(Genot=c(WT=color_panel[['black']], BKO=color_panel[['blue']]), Cond=c(cont=viridis(5)[5], link=viridis(5)[1]), score=viridis(256))
    output_prefix = 'pMEF_link'
    analyzeCountMatrix(matrix_file, design_file, output_prefix, col, log10_flag=T, mat=mat, design=design)
    genesetComparison(mat, design, paste0(output_prefix, '_lysosome'), col, dir=LYSOSOME, log10_flag=T)
    plotPCA(mat, design, output_prefix)

    mat = read.table("fold_change_logFC.txt", header=T, sep="\t")
    rownames(mat) = mat[,1]
    mat = mat[,-c(1:3)]
    mat = mat[,grep('_vs_WT_logFC', colnames(mat))]
    output_prefix = "pMEF_comp"
    require("RColorBrewer")
    col=list(Genot=c(WT=color_panel[['black']], BKO=color_panel[['blue']]), Cond=c(cont=viridis(5)[5], link=viridis(5)[1]),
             score=rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)))
    design = read.table('design2.txt', header=F, sep=" ")
    design = design[1:3,]
    rownames(design) = colnames(mat)
    analyzeCountMatrix(matrix_file, '', output_prefix, col, mat=mat, design=design, log10_flag=F)
    plotPCA(mat, design, output_prefix)

    for (tail in c('logFC', 'FDR')) {
        mat = read.table("fold_change_logFC.txt", header=T, sep="\t")
        rownames(mat) = mat[,1]
        mat = mat[,-c(1:3)]
        sample_list <- c('WT_link_vs_WT', 'BKO_link_vs_BKO')
        fc_mat = mat[,paste0(sample_list, '_logFC')]
        if (tail == 'FDR') {
            pmat <- mat[,paste0(sample_list, '_FDR')]
            fc_mat[pmat >= threshold] = 0
            fc_mat <- fc_mat[rowMins(as.matrix(fc_mat)) != rowMaxs(as.matrix(fc_mat)),]
            print(rowMeans(fc_mat))
        }
        output_prefix = 'pMEF_comp_cond'
        log10_flag = 0
        if (tail == 'FDR') {
            output_prefix = paste0(output_prefix, '_', tail)
        }
        design = read.table('design2.txt', header=F, sep=" ")
        design = design[c(2, 5),]
        rownames(design) = colnames(fc_mat)
        analyzeCountMatrix(matrix_file, '', 'pMEF_cond_comp', col, mat=mat, design=design, log10_flag=log10_flag)
        genesetComparison(fc_mat, design, paste0(output_prefix, '_lysosome'), col, dir=LYSOSOME, log10_flag=log10_flag)
    }
}

compareGuanosine <- function(threshold=0.05) {
    matrix_file = 'mat.txt'
    design_file = 'design.txt'
    mat = read.table(matrix_file, header=T, sep=" ")
    design = read.table(design_file, header=F, sep=" ")
    print(head(design))
    mat = mat[,!is.na(design[,1])]
    design = design[!is.na(design[,1]),]
    colnames(design) = c('Group', 'Genot', 'Cond')
    rownames(design) = colnames(mat[1:dim(mat)[2]])
    color_panel <- pal_aaas("default")(10)[c(10, 1, 2)]
    names(color_panel) = c('black', 'blue', 'red')
    col=list(Genot=c(WT=color_panel[['black']], FL=color_panel[['red']]), Cond=c(cont=viridis(5)[5], gua=viridis(5)[1]), score=viridis(256))
    output_prefix = 'pMEF_guanosine'
    analyzeCountMatrix(matrix_file, design_file, output_prefix, col, log10_flag=T, mat=mat, design=design)
    genesetComparison(mat, design, paste0(output_prefix, '_lysosome'), col, dir=LYSOSOME, log10_flag=T)
    plotPCA(mat, design, output_prefix)
}

compareMEF <- function(threshold=0.05) {
    matrix_file = 'mat.txt'
    design_file = 'design.txt'
    mat = read.table(matrix_file, header=T, sep=" ")
    design = read.table(design_file, header=F, sep=" ")
    colnames(design) = c('Group', 'Genot', 'Cond')
    rownames(design) = colnames(mat[1:dim(mat)[2]])
    color_panel <- pal_aaas("default")(10)[c(10, 1, 2)]
    names(color_panel) = c('black', 'blue', 'red')
    col=list(Genot=c(WT=color_panel[['black']], BKO=color_panel[['blue']]), Cond=c(cont=viridis(5)[5], link=viridis(5)[1]), score=viridis(256))
    output_prefix = 'MEF_link'
    analyzeCountMatrix(matrix_file, design_file, output_prefix, col, log10_flag=T, mat=mat, design=design)
    genesetComparison(mat, design, paste0(output_prefix, '_lysosome'), col, dir=LYSOSOME, log10_flag=T)
    plotPCA(mat, design, output_prefix)
}    


compareMultiome <- function() {
    return()
}

visualizeData <- function(mat, design, output_prefix) {
    return()
}

# compareLiver()
# compareBlood()
# compareMEF()
# compareGuanosine()


