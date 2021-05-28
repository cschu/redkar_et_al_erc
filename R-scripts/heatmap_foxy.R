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

transcripts_of_interest = read.csv("/Users/cschu/Documents/amey/foxy/DGE_SIGP_TRANSCRIPTS_UP_WITH_GOI.txt",
                                   header = F, sep="\t")
genes_of_interest = read.csv("/Users/cschu/Documents/amey/foxy/DGE_SIGP_GENES_UP_WITH_GOI.txt",
                             header = F, sep="\t")

organism = "FOXY"
t2gmap = "GCF_000149955.1_ASM14995v2_genomic.gff.trans_gene_map"
sample_pattern = "^Fol_[0-9A]+"


projectpath = "/Users/cschu/Documents/amey/foxy"
datapath = file.path(projectpath, "setup", organism)
tx2gene = file.path(datapath, t2gmap)
tx2gene = data.frame(read.csv(tx2gene, header = F, sep="\t"))
resultspath = file.path(projectpath, "results_december_xmas", organism)

sample_ids = basename(list.files(datapath, sample_pattern))
files = list.files(file.path(datapath, sample_ids), pattern="quant.sf", full.names = T)
names(files) = sample_ids

sample_order = c("Fol_A_1", "Fol_Ap_1", 
               "Fol1_1", "Fol1_2", "Fol1_3", 
               "Fol2_1", "Fol2_2", "Fol2_3",
               "Fol3_1", "Fol3_2", "Fol3_3",
               "Fol7_1", "Fol7_2", "Fol7_3")

files = files[order(match(names(files), sample_order))]


coldata = read.csv(file.path(datapath, "metadata.tsv"), header=T, sep="\t")
#coldata <- coldata[order(coldata$Sample),]
coldata = coldata[order(match(coldata$Sample, sample_order)), ]

# check that order of metadata entries matches file order (otherwise: chaos)
if (any(coldata$Sample != names(files))) {
  print(names(files))
  print(coldata$Sample)
  print(coldata$Sample != names(files))
  stop(c("Sample Order not identical between files and metadata"))
}

# transcript level analysis
colheader <- "Gene"
outpath <- file.path(resultspath, "genes")
#txi.salmon <- tximport(files, type="salmon", txOut = T)
txi.salmon <- tximport(files, type="salmon", txOut = F, tx2gene = tx2gene)

dir.create(outpath, showWarnings = F, recursive = T)

exp_design = formula(~TreatmentTime)
exp_intgroup = c("TreatmentTime") #c("Treatment", "Time")
dds = DESeqDataSetFromTximport(txi.salmon, coldata, design = exp_design)
  
# require more than one count per gene/transcript
dds = dds[ rowSums(counts(dds)) > 1, ]
  
rld = rlogTransformation(dds, blind = F)
keep = names(rld) %in% genes_of_interest$V1

#png(file.path(outpath, "heatmap_clustering.png"), 1500, 1500, pointsize = 20)
#cim(assay(rld)[keep, ], 
#xlab = "Genes", ylab = "Samples",
#color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), symkey = FALSE,
#cluster = "row", scale = 0.1,
#lhei = c(1,3))
#dev.off()

for (pal in c("RdBu", "RdGy", "PuOr", "PRGn", "PiYG", "BrBG")) {
  print(pal)
  pdf(file.path(resultspath, paste("heatmap_genes", pal, "pdf", sep=".")))
  hm=heatmap.2(assay(rld)[keep, ], 
               col=colorRampPalette(rev(brewer.pal(11, pal)))(100),
               trace = "none",  dendrogram = "row", cexRow=0.25, cexCol=0.5,
               lwid=c(0.05, 0.1), margins=c(5,5), Colv=F, key=F, Rowv=T, 
               offsetCol=0.1, offsetRow=0.05, labRow=F)
  
  dev.off()
}
pdf(file.path(resultspath, "heatmap_genes.pdf"))
hm=heatmap.2(assay(rld)[keep, ], 
             col=colorRampPalette(rev(brewer.pal(11, "RdYlGn")))(100),
             trace = "none",  dendrogram = "row", cexRow=0.25, cexCol=0.5,
             lwid=c(0.05, 0.1), margins=c(5,5), Colv=F, key=F, Rowv=T, 
             offsetCol=0.1, offsetRow=0.05, labRow=F)
dev.off()    
rowseps = c(2,5,8,9,19,25,35,39,51,52,183,186,191,192,207,208)
rowseps_visual = c(132, 155, 235, 300)
rowseps_visual = c(2,13,23,44,85,87,168,172,258,262)
rowseps = c(327, 337, 330, 296, 254, 252, 170, 168, 130, 129, 81, 78)
rowseps = 283 - c(2, 6, 9, 13, 32, 37, 84, 86, 160, 164, 222, 225)

rowseps = 275 - c(75,76,160,167,191,196,223,227,231,234)

pdf(file.path(resultspath, "heatmap_genes_with_marks.pdf"))
hm=heatmap.2(assay(rld)[keep, ], 
             col=colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), 
             trace = "none",  dendrogram = "row", cexRow=0.25, cexCol=0.5,
             lwid=c(0.05, 0.1), margins=c(5,5), Colv=T, key=F, Rowv=T, 
             offsetCol=0.1, offsetRow=0.05, labRow=F,
             rowsep = rowseps, sepwidth=c(0.01,0.01), sepcolor = "black")
dev.off()

genes = rownames(assay(rld)[keep,])[hm$rowInd]
write.table(
  data.frame(gene=genes),
  file.path(resultspath, "gene_order.txt"),
  row.names = F, quote= F, sep = "\t")

sessionInfo()


