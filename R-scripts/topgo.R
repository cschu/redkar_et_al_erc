source('http://bioconductor.org/biocLite.R')
biocLite('topGO')

BiocManager::install("Rgraphviz")

library(topGO)
library(Rgraphviz)

projectpath = "/Users/cschu/Documents/amey/foxy"
datapath = file.path(projectpath, "results_2021/FOXY/genes") #/1dpi.up.gids
outpath = file.path(datapath, "topgo")
files = list.files(datapath, pattern = "[1237]dpi.up.gids", full.names = T)

geneID2GO = readMappings(file.path(projectpath, "/setup/FOXY/GCF_000149955.1_ASM14995v2_genomic.gff.geneid2go.topgo"))
geneUniverse <- names(geneID2GO)

genes_of_interest = as.character(read.table(file.path(datapath, "genes_of_interest_for_topgo.txt"))$V1)
gene_list = factor(as.integer(geneUniverse %in% genes_of_interest))
names(gene_list) = geneUniverse

nnodes = 30

for (ont in c("BP", "MF", "CC")) {
  myGOdata = new("topGOdata", description=paste("foxy_up", ont, sep="_"), 
                 ontology=ont, allGenes=gene_list,
                 annot = annFUN.gene2GO, gene2GO=geneID2GO)
  sg = sigGenes(myGOdata)
  numSigGenes(myGOdata)
  resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher")
  allRes = GenTable(myGOdata, topgoFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = nnodes)
  write.table(allRes, file=file.path(outpath, paste("combined", ont, "topGO.results.tsv", sep=".")), sep="\t")
  
  showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = nnodes, useInfo ='all')
  printGraph(myGOdata, resultFisher, firstSigNodes = nnodes, fn.prefix = file.path(outpath, paste("combined", ont, "topGO", sep=".")), useInfo = "all", pdfSW = TRUE)
  
  # print out the genes that are annotated with the significantly enriched GO terms:
  myterms <- allRes$GO.ID
  mygenes <- genesInTerm(myGOdata, myterms)
  for (i in 1:length(myterms))
  {
    myterm <- myterms[i]
    mygenesforterm <- mygenes[myterm][[1]]
    myfactor <- mygenesforterm %in% genes_of_interest # find the genes that are in the list of genes of interest
    mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
    mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
    cat(paste("Term",myterm,"genes:",mygenesforterm2), 
        file=file.path(outpath, paste("combined", ont, "topGO.txt", sep=".") ),
        append=TRUE)
  }
  
  
}

for (gid_file in files) {
  for (ont in c("BP", "MF", "CC")) {
    # print(gid_file)
    genes_of_interest = as.character(read.table(gid_file)$V1)
    head(genes_of_interest)
    gene_list = factor(as.integer(geneUniverse %in% genes_of_interest))
    names(gene_list) = geneUniverse
  
    myGOdata = new("topGOdata", description=paste("foxy_up", ont, sep="_"), 
                   ontology=ont, allGenes=gene_list,
                   annot = annFUN.gene2GO, gene2GO=geneID2GO)
    sg = sigGenes(myGOdata)
    numSigGenes(myGOdata)
    resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher")
    allRes = GenTable(myGOdata, topgoFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = nnodes, rm.one = TRUE)
    write.table(allRes, file=file.path(outpath, paste(basename(gid_file), ont, "topGO.results.tsv", sep=".")), sep="\t")
    showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = nnodes, useInfo ='all')
    printGraph(myGOdata, resultFisher, firstSigNodes = nnodes, fn.prefix = file.path(outpath, paste(basename(gid_file), ont, "topGO", sep=".")), useInfo = "all", pdfSW = TRUE)
    
    
    
    # print out the genes that are annotated with the significantly enriched GO terms:
    myterms <- allRes$GO.ID
    mygenes <- genesInTerm(myGOdata, myterms)
    for (i in 1:length(myterms))
    {
      myterm <- myterms[i]
      mygenesforterm <- mygenes[myterm][[1]]
      myfactor <- mygenesforterm %in% genes_of_interest # find the genes that are in the list of genes of interest
      mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
      mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
      cat(paste("Term",myterm,"genes:",mygenesforterm2), 
          file=file.path(outpath, paste(basename(gid_file), ont, "topGO.txt", sep=".") ),
          append=TRUE)
    }
  }
    
}
