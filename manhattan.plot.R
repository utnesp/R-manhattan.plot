plot.manhattan <- function(raw.count.matrix, design, name = "", file.loc = "", title.name = "", count.mean.cutoff = 10, logFC.cutoff = 0.8, Pvalue.cutoff = 0.05, write.all.genes = T, reduce.point.size.by = 6, 
                           dpi = 300, width = 297, height = 210, units = "mm", open.pdf = F, open.html = F, which.attribute = 1, reverseGroup = F,
                           highligt.top.genes = F, nr_top_genes = 10, nr_bottom_genes = 10, which.biotypes = "", inverse.biotypes = F, order.by.pvalue = F, write.cpm = F) {
    library(stringr) 
    print("Design data.frame consists with row.names = samples, and columns for each clinical attribute")
    print("The design can only have a maximum of two levels")
    print("Any NA attributes will be removed")
    
    # raw.count.matrix = all.gene.counts.txt.filtered.sorted2; design = t; name = ""; file.loc = "/Volumes/Peter/SRA/manhattan/"; title.name = ""; count.mean.cutoff = 10; logFC.cutoff = 0.8; Pvalue.cutoff = 0.05; write.all.genes = F; reduce.point.size.by = 6; 
    # dpi = 300; width = 297; height = 210; units = "mm"; which.attribute = 10; open.pdf = F; open.html = F
    
    if (substr(file.loc, nchar(file.loc), nchar(file.loc)) != "/") {
        file.loc <- paste(file.loc, "/", sep = "")
    }
    
    design_manhattan <- data.frame(attributes = factor(design[, which.attribute]), samples = row.names(design))
    
    if (name == "") {
        name <- colnames(design)[which.attribute]
        name <- gsub(as.character(levels(design_manhattan$attributes)[2]), "", name)
        name <- gsub(".", "", name, fixed = T)
    }
    
    count.matrix <- raw.count.matrix[colnames(raw.count.matrix) %in% design_manhattan$samples]
    
    if (reverseGroup == T) {
        print(paste("Testing condition: ", levels(design_manhattan$attributes)[2], " / ", levels(design_manhattan$attributes)[1], sep = ""))
        design_manhattan <- design_manhattan[design_manhattan$attributes == levels(design_manhattan$attributes)[1], "samples"]
        design_manhattan <- design_manhattan[!is.na(design_manhattan)]
        group <- ifelse(colnames(count.matrix) %in% design_manhattan, 1, 2)
    } else {
        print(paste("Testing condition: ", levels(design_manhattan$attributes)[1], " / ", levels(design_manhattan$attributes)[2], sep = ""))
        design_manhattan <- design_manhattan[design_manhattan$attributes == levels(design_manhattan$attributes)[2], "samples"]
        design_manhattan <- design_manhattan[!is.na(design_manhattan)]
        group <- ifelse(colnames(count.matrix) %in% design_manhattan, 1, 2)
    }
    
    y <- DGEList(count.matrix, group = group)
    y <- calcNormFactors(y)
    
    genes <- row.names(cpm(y)[rowMeans(cpm(y)) > count.mean.cutoff, ])
    print(paste("Keeping ", length(genes), " genes out of ", nrow(count.matrix), " from the raw count matrix", sep = ""))
         
        if ( exists("ensembl.attributes") == "FALSE") {
            print("Quering bioMart ... this may take a long time")
            ensembl.attributes <<- data.frame(ensembl_gene_id = row.names(count.matrix))
            ensembl.attributes <<- ensg2ext_name_biotype(ensembl.attributes$ensembl_gene_id, combine = T, all = T)
            ensembl.attributes <<- ensg2chr_start_end(ensembl.attributes$ensembl_gene_id, combine = T, all = T, as.IGV = T)
            ensembl.attributes <<- ensg2source_status(ensembl.attributes$ensembl_gene_id, combine = T, all = T)
        } else {
            print("Using existing object 'ensembl.attributes' and 'gene.position")
        }
        
        if (exists("gene.position") == "FALSE") {
            gene.position <<- data.frame(ensembl_gene_id = row.names(count.matrix))
            gene.position <<- ensg2chr_start_end(ensembl.attributes$ensembl_gene_id, combine = T, all = T)
        }
    
    if (which.biotypes != "") {
        count.matrix <- count.matrix[row.names(count.matrix) %in% genes, ]
        lncRNAs <- ensg2gene_biotype(row.names(count.matrix))
        print("Distribution of gene biotypes: ")
        print(table(lncRNAs$gene_biotype))
            if (inverse.biotypes == T) {
            lncRNAs <- lncRNAs[!lncRNAs$gene_biotype %in% which.biotypes, ]
            } else {
            lncRNAs <- lncRNAs[lncRNAs$gene_biotype %in% which.biotypes, ]
            }
        
        count.matrix <- count.matrix[row.names(count.matrix) %in% lncRNAs$ensembl_gene_id, ]
        } else {
        count.matrix <- count.matrix[row.names(count.matrix) %in% genes, ]
    }
    
    y <- DGEList(count.matrix, group = group)
    y <- calcNormFactors(y)
    
    design <- model.matrix(~group)
    y <- estimateDisp(y, design = design)
    
    y.fit <- glmFit(y, design)
    y.chr <- glmLRT(y.fit, coef=2)
    q <- topTags(y.chr, n = Inf)
    q <- q$table
    q$ensembl_gene_id <- row.names(q)
    
    q <- merge(ensembl.attributes, q, all.y = T, by = "ensembl_gene_id")
    
        if (order.by.pvalue == T) {
        q <- q[order(q$PValue), ]
        } else {
        q <- q[order(q$FDR), ]
        }
    
    q <- q[c('external_gene_name', 'gene_biotype', 'ensembl_gene_id', 'source','status','position', 'band', 'strand','logFC','logCPM','LR','PValue','FDR')]
    
    k <- merge(gene.position, q[c("ensembl_gene_id", "band", "strand", "logFC", "logCPM", "LR", "PValue", "FDR")], by = "ensembl_gene_id")
    
    if (write.all.genes == F) {
        q <- q[abs(q$logFC) > logFC.cutoff, ]
        q <- q[q$FDR < Pvalue.cutoff, ]
    }
    
    cat("Differentially expressed genes:\n")
    print(head(q))
    
    if (title.name == "") {
        if (reverseGroup == T) {
        title.name <- paste(name, " n=", length(design_manhattan), " (", round(length(design_manhattan) / ncol(count.matrix), 2) * 100, "%)", sep = "")
        } else {
        nr_conditions = ncol(count.matrix) - length(design_manhattan)
        title.name <- paste(name, " n=", nr_conditions, " (", round(nr_conditions / ncol(count.matrix), 2) * 100, "%)", sep = "")
        }
    }
    
    
        if (order.by.pvalue == T) {
        l <- data.frame(ensembl_gene_id = k$ensembl_gene_id, CHR = k$chromosome_name, BP = round(as.integer(k$end_position - ((k$end_position - k$start_position) / 2)), 0), P = k$PValue, biotype = k$gene_biotype, logFC = k$logFC, external_gene_name = k$external_gene_name, logCPM = k$logCPM, band = paste(k$chromosome_name, k$band, sep = ""))
        } else {
        l <- data.frame(ensembl_gene_id = k$ensembl_gene_id, CHR = k$chromosome_name, BP = round(as.integer(k$end_position - ((k$end_position - k$start_position) / 2)), 0), P = k$FDR, biotype = k$gene_biotype, logFC = k$logFC, external_gene_name = k$external_gene_name, logCPM = k$logCPM, band = paste(k$chromosome_name, k$band, sep = ""))
        }
    
    
    
    l <- l[!is.na(l$BP), ]
    l <- l[abs(l$logFC) > logFC.cutoff, ]
    l <- l[l$P < Pvalue.cutoff, ]
    l <- l[!is.na(l$ensembl_gene_id), ]
    l <- l[l$CHR != "MT", ]
    l <- l[grep("KI", l$CHR, invert = T), ]
    l <- l[grep("GL", l$CHR, invert = T), ]
    l <- l[grep("CHR_", l$CHR, invert = T), ]
    
    l$CHR <- as.character(l$CHR)
    l$CHR <- ifelse(l$CHR == "X", 23, l$CHR)
    l$CHR <- ifelse(l$CHR == "Y", 24, l$CHR)
    l$CHR <- as.integer(l$CHR)
    l <- l[l$CHR <= 22, ]

    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector <- col_vector[-4]
    
    l$u <- col_vector[l$CHR]
    
    l$pos <- as.numeric(paste(l$CHR, str_pad(l$BP, nchar(max(l$BP)), pad = "0"), sep = "."))
    l <- l[order(l$pos), ]
    l$pos_nr <- 1:length(l$pos)
    
    ## use z to position x axis tick marks and labels
    z <- as.data.frame(cumsum(table(l$CHR))) 
    z$z.half <- table(l$CHR) / 2
    colnames(z) <- c("cumsum", "z.half")
    z$z <- z$cumsum - z$z.half
    
    ## assign square to protein coding genes and triangles to the rest
    l$pch <- ifelse(l$biotype == "protein_coding", 16, 17)
    
    
    
    p <- ggplot(l, aes(pos_nr, logFC)) + 
        geom_point(shape = l$pch, size = -log10(l$P) / abs(reduce.point.size.by), color = l$u, aes(text = paste("GENE: ", external_gene_name, "P value: ", signif(P,2), "BIOTYPE: ", biotype, "GENE_ID:", ensembl_gene_id, "BAND: ", band, "logCPM:", signif(logCPM, 2), sep =" "))) + 
        scale_x_continuous(breaks = z$z) + labs(x = "", title = title.name) + theme_classic() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA, colour = "black", size = 1))
    
    if (file.loc == "") {
        file.loc <- paste(getwd(), "/", sep = "") 
    }
    
    pdf.out <- paste(file.loc, "pdf/", sep = "")
    print(paste("Sending pdf's to", pdf.out))
    system(paste("mkdir -p", file.loc))
    system(paste("mkdir -p", pdf.out))
    file.pdf <- paste(pdf.out, name, ".pdf", sep = "")
    system(paste("trash", file.pdf))
    ggsave(file.pdf, p, device = "pdf", dpi = dpi, width = width, height = height, units = units)
    if (open.pdf == T) {
        openPDF(file.pdf)
    }
    
    print("Write html's to separate folder")
    p <- ggplotly(p)
    html.out <- paste(file.loc, "html/", sep = "")
    print(paste("Sending html's to", html.out))
    system(paste("mkdir -p", html.out))
    file.html <- paste(html.out, name, ".html", sep = "")
    system(paste("trash", file.html))
    htmlwidgets::saveWidget(p, file.html, selfcontained = T, libdir = NULL)
    if (open.html == T) {
        system(paste("open", file.html))
    }
    
    if (highligt.top.genes == T) {
        top <- l[l$P < Pvalue.cutoff & l$logFC > 0 , c("CHR", "BP", "pos_nr","logFC", "P", "external_gene_name")]
        top <- top[order(top$P), ]
        top <- top[1:nr_top_genes, ] 
        top <- top[!is.na(top$pos_nr) | !is.na(top$logFC), ] 
        bottom <- l[l$P < Pvalue.cutoff & l$logFC < 0 , c("CHR", "BP", "pos_nr", "logFC", "P", "external_gene_name")]
        bottom <- bottom[order(bottom$P), ]
        bottom <- bottom[1:nr_bottom_genes, ] 
        bottom <- bottom[!is.na(bottom$pos_nr) | !is.na(bottom$logFC), ] 
        top.genes <- rbind(top, bottom)
        
        p <- ggplot(l, aes(pos_nr, logFC)) + 
            geom_point(shape = l$pch, size = -log10(l$P) / abs(reduce.point.size.by), color = l$u, aes(text = paste("GENE: ", external_gene_name, "P value: ", signif(P,2), "BIOTYPE: ", biotype, "GENE_ID:", ensembl_gene_id, "BAND: ", band, "logCPM:", signif(logCPM, 2), sep =" "))) + 
            scale_x_continuous(breaks = z$z) + labs(x = "", title = title.name) + theme_classic() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA, colour = "black", size = 1)) + 
            geom_text_repel(data = top.genes, aes(pos_nr, logFC, label = top.genes$external_gene_name))
        system(paste("mkdir -p ", pdf.out, "labeled/", sep =""))
        file.pdf <- paste(pdf.out, "labeled/", name, ".labeled.pdf", sep = "")
        system(paste("trash", file.pdf))
        ggsave(file.pdf, p, device = "pdf", dpi = dpi, width = width, height = height, units = units)
        
        if (open.pdf == T) {
            openPDF(file.pdf)
        }
    }
    
    ### Write differential tests to separate folder
    diff.out <- paste(file.loc, "edgeR.GLM.tests/", sep = "")
    system(paste("mkdir -p", diff.out))
    file.diff = paste(diff.out, "edgeR.glmLRT.", name, ".txt", sep = "")
    system(paste("trash", file.diff))
    write.delim(q, file.diff)

    if (write.cpm == T) {
        file.cpm = paste(diff.out, "normalized.counts.", name, ".txt", sep = "")
        write.delim(cpm(y), file.cpm, row.names = T)
    }
    
    print(paste(table(q$FDR < 0.05), "genes with FDR < 0.05"))
    print("First 20 genes:")
    print(head(q, 20))
}
