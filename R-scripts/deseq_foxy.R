library("tximport")
library("DESeq2")
library("pheatmap")
library("gplots")
library("RColorBrewer")
library("ggplot2")
library("calibrate")
library("tidyverse")
library("mixOmics")
# https://support.bioconductor.org/p/103885/


organism = "FOXY"
t2gmap = "GCF_000149955.1_ASM14995v2_genomic.gff.trans_gene_map"
sample_pattern = "^Fol_[0-9A]+"

projectpath = "/Users/cschu/Documents/amey/foxy"
datapath = file.path(projectpath, "setup", organism)
tx2gene = file.path(datapath, t2gmap)
tx2gene = data.frame(read.csv(tx2gene, header = F, sep="\t"))
resultspath = file.path(projectpath, "results", organism)

sample_ids = basename(list.files(datapath, sample_pattern))
files = list.files(file.path(datapath, sample_ids), pattern="quant.sf", full.names = T)
names(files) = sample_ids
files = files[order(names(files))]

coldata = read.csv(file.path(datapath, "metadata.tsv"), header=T, sep="\t")
coldata <- coldata[order(coldata$Sample),]

# check that order of metadata entries matches file order (otherwise: chaos)
if (any(coldata$Sample != names(files))) {
  print(names(files))
  print(coldata$Sample)
  print(coldata$Sample != names(files))
  stop(c("Sample Order not identical between files and metadata"))
}

experiments = list()
experiments[[1]] = "TreatmentTime_trt1_vs_ctrl"
experiments[[2]] = "TreatmentTime_trt2_vs_ctrl"
experiments[[3]] = "TreatmentTime_trt3_vs_ctrl"
experiments[[4]] = "TreatmentTime_trt7_vs_ctrl"

# we run gene-level analysis first, then transcript
# only difference is that gene-level requires a map of transcripts to genes for the import

# import expression data for genes/transcripts
colheader <- "Gene"
outpath <- file.path(resultspath, "genes")
txi.salmon <- tximport(files, type="salmon", txOut = F, tx2gene = tx2gene)
  
dir.create(outpath, showWarnings = F, recursive = T)

#exp_design = formula(~ Time + Treatment + Time:Treatment)
exp_design = formula(~TreatmentTime)
exp_intgroup = c("TreatmentTime") #c("Treatment", "Time")
dds = DESeqDataSetFromTximport(txi.salmon, coldata, design = exp_design)
  
# require more than one count per gene/transcript
dds = dds[ rowSums(counts(dds)) > 1, ]

rld = rlogTransformation(dds, blind = F)

  # now run pairwise comparisons
for (exp_contrast in experiments) {
    
    if (typeof(exp_contrast) == "list") {
      deseq_res <- results(deseq, contrast=exp_contrast) #, alpha=0.1, lfcThreshold=1)
      # filename prefix
      print(exp_contrast)
      prefix <- paste(substring(exp_contrast[[1]][1], 1), "_vs_", substring(exp_contrast[[2]][1], 1), sep="")
    } else {
      deseq_res <- results(deseq, name=exp_contrast)#, alpha=0.1, lfcThreshold=1)
      prefix <- exp_contrast
    }
    
    
    de = (deseq_res$padj < 0.00000000001 & deseq_res$log2FoldChange > 2)
    de[is.na(de)] = FALSE
    
    png(file.path(outpath, paste(prefix, "-cim.png", sep="")), 1500, 1500, pointsize = 20)
    cim(t(assay(rld)[de, ]), transpose = T, keysize = c(0.5, 0.5),
        xlab = "Genes", ylab = "Samples",
        color = colorRampPalette(brewer.pal(9, "YlGnBu"))(255), symkey = FALSE,
        lhei = c(1,3))
    dev.off()
  }    



sessionInfo()