install.packages("ggVennDiagram")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

library(ggVennDiagram)
library(ggplot2)
library(ggvenn)

datapath = "/Users/cschu/Documents/amey/foxy/results_december/FOXY/genes"

set1 = read.table(file.path(datapath, "TreatmentTime_trt1_vs_ctrl.filtered.csv.sigp_hc.csv.l2fup.csv.tids"))
set2 = read.table(file.path(datapath, "TreatmentTime_trt2_vs_ctrl.filtered.csv.sigp_hc.csv.l2fup.csv.tids"))
set3 = read.table(file.path(datapath, "TreatmentTime_trt3_vs_ctrl.filtered.csv.sigp_hc.csv.l2fup.csv.tids"))
set7 = read.table(file.path(datapath, "TreatmentTime_trt7_vs_ctrl.filtered.csv.sigp_hc.csv.l2fup.csv.tids"))

x <- list(dpi1=set1$V1, dpi2=set2$V1, dpi3=set3$V1, dpi7=set7$V1)
          
  ggVennDiagram(x, 
                category.names=c("day1", "day2", "day3", "day7"),
                label="count", color="black") +
    scale_fill_gradient(low="blue", high="red")

  names(x) <- c("1dpi", "2dpi", "3dpi", "7dpi")
  pdf(file.path(resultspath, "venn.pdf"))
  ggvenn(x, show_percentage = F, stroke_size = 0, stroke_alpha=1, set_name_size = 4)
  dev.off()
  
  
  
  
datapath = "/Users/cschu/Documents/amey/foxy/results_december_end/SLYC/genes"
set1 = read.table(file.path(datapath, "TreatmentTime_trt1_ctrl1.filtered.up.csv.l2fup.csv.tids"))
set2 = read.table(file.path(datapath, "TreatmentTime_trt2_ctrl2.filtered.up.csv.l2fup.csv.tids"))
set3 = read.table(file.path(datapath, "TreatmentTime_trt3_ctrl3.filtered.up.csv.l2fup.csv.tids"))
set7 = read.table(file.path(datapath, "TreatmentTime_trt7_ctrl7.filtered.up.csv.l2fup.csv.tids"))

x <- list(dpi1=set1$V1, dpi2=set2$V1, dpi3=set3$V1, dpi7=set7$V1)
names(x) <- c("1dpi", "2dpi", "3dpi", "7dpi")
ggvenn(x, show_percentage = F, stroke_size = 0, stroke_alpha=1, set_name_size = 4)
pdf(file.path(resultspath, "tomato_genes_venn.pdf"))
ggvenn(x, show_percentage = F, stroke_size = 0, stroke_alpha=1, set_name_size = 4)
dev.off()