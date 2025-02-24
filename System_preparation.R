##############################################################
################## System preparation ########################
##############################################################

rm(list=ls())

if (!require("aricode", quietly = TRUE))
  install.packages("aricode")

if (!require("clevr", quietly = TRUE))
  install.packages("clevr")

if (!require("cluster", quietly = TRUE))
  install.packages("cluster")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("ggvis", quietly = TRUE))
  install.packages("ggvis")

if (!require("igraph", quietly = TRUE))
  install.packages("igraph")

if (!require("mcclust", quietly = TRUE))
  install.packages("mcclust")

if (!require("mclust", quietly = TRUE))
  install.packages("mclust")

if (!require("RColorBrewer", quietly = TRUE))
  install.packages("RColorBrewer")

if (!require("R6", quietly = TRUE))
  install.packages("R6")

if (!require("tibble", quietly = TRUE))
  install.packages("tibble")

if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")

if (!require("umap", quietly = TRUE))
  install.packages("umap")