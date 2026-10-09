#!/usr/bin/env Rscript
# Secondary analysis of deposited public mouse counts; not a patient experiment.
suppressPackageStartupMessages({library(edgeR); library(limma)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2L)
inputs <- args[1]; out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep="\t",
                                         quote=FALSE, row.names=FALSE, na="NA")
raw <- read.delim(gzfile(file.path(inputs, "GSE134781_counts.txt.gz")),
                  row.names=1, check.names=FALSE)
counts <- as.matrix(raw)
stopifnot(ncol(counts) == 9L, nrow(counts) == 24062L,
          !anyDuplicated(rownames(counts)), all(is.finite(counts)),
          all(counts >= 0), all(counts == round(counts)))
samples <- colnames(counts)
group <- ifelse(grepl("hetX753", samples), "X753",
                ifelse(grepl("hetL1002P", samples), "L1002P",
                       ifelse(grepl("WT", samples), "WT", NA)))
stopifnot(!anyNA(group), all(table(group) == 3L))
group <- factor(group, levels=c("WT", "L1002P", "X753"))
gmt <- strsplit(readLines(file.path(inputs, "mh.all.v2025.1.Mm.symbols.gmt")), "\t")
sets <- setNames(lapply(gmt, function(x) unique(x[-c(1,2)])),
                 vapply(gmt, `[`, "", 1L))
stopifnot(length(sets) == 50L, !anyDuplicated(names(sets)))
design_for <- function(g) {d <- model.matrix(~0+g); colnames(d) <- levels(g); d}
contrasts_for <- function(d) makeContrasts(X753_vs_WT=X753-WT,
        L1002P_vs_WT=L1002P-WT, X753_vs_L1002P=X753-L1002P, levels=d)
design <- design_for(group)
primary_keep <- filterByExpr(DGEList(counts=counts, group=group), design=design)

analyse <- function(keep, columns=seq_along(group), normalization="TMM",
                    correlation=NA_real_, save_genes=FALSE, label="primary") {
    d <- design_for(group[columns]); contrasts <- contrasts_for(d)
    y <- DGEList(counts=counts[keep,columns,drop=FALSE], group=group[columns])
    y <- calcNormFactors(y, method=normalization)
    v <- voom(y, d, plot=FALSE)
    fit <- eBayes(contrasts.fit(lmFit(v, d), contrasts), robust=TRUE)
    idx <- ids2indices(sets, rownames(v))
    stopifnot(length(idx) == 50L, all(lengths(idx) >= 10L))
    gene_results <- list(); pathway_results <- list()
    for (j in seq_len(ncol(contrasts))) {
        gene <- topTable(fit, coef=j, number=Inf, sort.by="none", confint=TRUE)
        gene$gene <- rownames(gene); gene$contrast <- colnames(contrasts)[j]
        gene_results[[j]] <- gene
        cam <- camera(v, idx, d, contrast=contrasts[,j], inter.gene.cor=correlation,
                      allow.neg.cor=FALSE, use.ranks=FALSE, sort=FALSE)
        cam$set <- rownames(cam); cam$contrast <- colnames(contrasts)[j]
        names(cam)[names(cam)=="FDR"] <- "FDR_within_contrast"
        cam$median_log2FC <- vapply(idx, function(i) median(gene$logFC[i]), 0.0)
        cam$original_set_size <- lengths(sets)[cam$set]
        cam$analysis <- label
        # Fixed-correlation mode does not return a Correlation column.
        if (!"Correlation" %in% names(cam)) cam$Correlation <- correlation
        cam$effective_correlation <- pmax(0, cam$Correlation)
        pathway_results[[j]] <- cam
    }
    pathways <- do.call(rbind, pathway_results)
    pathways$FDR_global150 <- p.adjust(pathways$PValue, "BH")
    genes <- do.call(rbind, gene_results)
    genes$FDR_global_genes <- p.adjust(genes$P.Value, "BH")
    if(save_genes) {
        con <- gzfile(file.path(out, "gene-results.tsv.gz"), "wt")
        write.table(genes, con, sep="\t", quote=FALSE, row.names=FALSE); close(con)
    }
    list(pathways=pathways, genes=genes, v=v, y=y, idx=idx)
}
primary <- analyse(primary_keep, save_genes=TRUE)
pdf(file.path(out,"voom-quality-check.pdf"), width=8, height=6)
invisible(voom(primary$y, design, plot=TRUE)); dev.off()
png(file.path(out,"voom-quality-check.png"), width=1200, height=900, res=150)
invisible(voom(primary$y, design, plot=TRUE)); dev.off()
write_tsv(primary$pathways, "pathways.tsv")
write_tsv(data.frame(gene=rownames(counts), retained=primary_keep), "gene-universe.tsv")
write_tsv(data.frame(sample=samples, genotype=group, library_counts=colSums(counts),
    retained_counts=primary$y$samples$lib.size,
    TMM_factor=primary$y$samples$norm.factors), "samples.tsv")

# Prespecified sensitivity checks, each reruns the full 150-test family.
sens <- list()
for (i in seq_along(group)) {
    label <- paste0("omit_", samples[i])
    sens[[label]] <- analyse(primary_keep, columns=setdiff(seq_along(group),i),
                            label=label)$pathways
}
sens$paper_filter <- analyse(rowMeans(counts) >= 10, label="paper_mean_count_ge10")$pathways
sens$RLE <- analyse(primary_keep, normalization="RLE", label="RLE_normalization")$pathways
sens$fixed_correlation <- analyse(primary_keep, correlation=0.01,
                                  label="fixed_correlation_0.01")$pathways
write_tsv(do.call(rbind,sens), "sensitivity.tsv")

# Membership and overlap make connected signatures visible; no drug ranking.
membership <- do.call(rbind, lapply(names(primary$idx), function(n)
    data.frame(set=n, gene=rownames(primary$v)[primary$idx[[n]]])))
write_tsv(membership, "tested-membership.tsv")
ov <- expand.grid(set1=names(sets), set2=names(sets), stringsAsFactors=FALSE)
ov$intersection <- mapply(function(a,b) length(intersect(primary$idx[[a]], primary$idx[[b]])),
                           ov$set1, ov$set2)
ov$jaccard <- mapply(function(a,b) length(intersect(primary$idx[[a]], primary$idx[[b]])) /
                    length(union(primary$idx[[a]], primary$idx[[b]])), ov$set1, ov$set2)
write_tsv(ov, "set-overlap.tsv")
markers <- intersect(c("Bub1b","Mtor","Rptor","Rps6kb1","Eif4ebp1","Cdkn2a",
                        "Cdkn1a","Xist","Ddx3y","Eif2s3y","Kdm5d","Uty"), rownames(counts))
write_tsv(data.frame(gene=markers, counts[markers,,drop=FALSE], check.names=FALSE),
          "marker-counts.tsv")

# A focused display was chosen before results: all 50 sets are in the tables.
display <- paste0("HALLMARK_", c("MTORC1_SIGNALING","PI3K_AKT_MTOR_SIGNALING",
    "OXIDATIVE_PHOSPHORYLATION","REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "UNFOLDED_PROTEIN_RESPONSE","MYOGENESIS","P53_PATHWAY",
    "INFLAMMATORY_RESPONSE","IL6_JAK_STAT3_SIGNALING"))
draw <- function() {
    layout(matrix(c(1,2),1,2), widths=c(1,1.8))
    cols <- c(WT="#3A5668", L1002P="#D58A25", X753="#AB4370")
    par(mar=c(5,4,4,1))
    pca <- prcomp(t(primary$v$E), center=TRUE, scale.=FALSE)
    explained <- pca$sdev^2/sum(pca$sdev^2)*100
    yr <- range(pca$x[,2]); yr <- yr + c(-1,1)*diff(yr)*.12
    plot(pca$x[,1:2], pch=19, col=cols[group], cex=1.4, ylim=yr,
         xlab=sprintf("PC1 (%.1f%%)",explained[1]), ylab=sprintf("PC2 (%.1f%%)",explained[2]),
         main="Nine independent mouse samples")
    text(pca$x[,1:2], labels=sub(".*_", "", samples), pos=3, cex=.8)
    legend("bottomleft", names(cols), pch=19, col=cols, bty="n", cex=.8)
    par(mar=c(6,14,4,2))
    ct <- c("X753_vs_WT","L1002P_vs_WT","X753_vs_L1002P")
    plot(NA, xlim=c(.5,3.5), ylim=c(.5,length(display)+.5), axes=FALSE,
         xlab="", ylab="", main="Transcriptional programs; not kinase activity")
    axis(1, at=1:3, labels=gsub("_vs_", " vs ", ct), las=1, cex.axis=.75)
    axis(2, at=seq_along(display), labels=gsub("_"," ",sub("HALLMARK_","",display)),
         las=2, cex.axis=.75)
    abline(h=seq_along(display), col="#EEEEEE")
    for (x in 1:3) for (y in seq_along(display)) {
        z <- primary$pathways[primary$pathways$set==display[y] &
                               primary$pathways$contrast==ct[x],]
        color <- if(z$Direction=="Up") "#AD426D" else "#246D9B"
        points(x,y,pch=if(z$FDR_global150<.05)19 else 1, col=color, cex=1.5)
        text(x+.16,y,sprintf("%.2g",z$FDR_global150),adj=0,cex=.63)
    }
    mtext("Pink: up; blue: down. Filled: global FDR < 0.05; numbers: FDR.\n50 sets x 3 contrasts; CAMERA estimates within-set residual correlation.",
           side=1,line=3.8,cex=.67)
}
pdf(file.path(out,"mouse-pathway-analysis.pdf"),width=13,height=6); draw(); dev.off()
png(file.path(out,"mouse-pathway-analysis.png"),width=1950,height=900,res=150); draw(); dev.off()
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
           file.path(out,"session-info.txt"))
writeLines(c(paste("Input genes:", nrow(counts)), paste("Retained genes:", sum(primary_keep)),
    "Primary: TMM, filterByExpr, voom; robust eBayes for gene tables only.",
    "CAMERA uses its own moderation, estimates residual correlation and floors variance inflation at 1.",
    "Genotype-only model; unknown batch/litter/sex cannot be adjusted as known metadata.",
    "Hallmark sets: 50; contrasts: 3; BH family: 150.",
    "Sensitivity checks are robustness analyses, not extra discovery opportunities."),
    file.path(out,"analysis-summary.txt"))
