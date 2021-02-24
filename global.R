library(shiny)
library(shinythemes)
library(shinyWidgets)
library(SpatialDecon)
library(shinyjs)
library(imager)
library(DT)
library(gplots)
library(dendextend)
library(plotrix)
library(RColorBrewer)


expression <- read.csv("data/expression.csv")
ROI <- read.csv("data/ROI.csv")
color_celltypes <- data.frame(c("#7f5d77", "#73b306", "#f30bbd", "#13b367", "#5e0383", "#fdb62c",
                                "#0d47c0", "#0b5e03", "#fc7de1", "#6f7105", "#281f52", "#11c1cb",
                                "#db0b55", "#8c3902", "#d9b6fd"),
                              c("macrophages", "mast", "B", "plasma", "CD4.T.cells", "CD8.T.cells", 
                                "NK", "pDC", "mDCs", "monocytes", "neutrophils", "Treg", "endothelial.cells", 
                                "fibroblasts", "cancer"))
colnames(color_celltypes) <- c("color", "celltype")
color_cluster <- data.frame(c("#77b241", "#b759c0", "#52a875", "#d14150", "#4dadcf",
                              "#c86a41", "#7276cb", "#d0a245", "#c3628d", "#7f7b35"), 1:10)
colnames(color_cluster) <- c("color", "cluster_num")




