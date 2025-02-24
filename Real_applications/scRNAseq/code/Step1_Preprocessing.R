##############################################################
####################### Pre-processing #######################
##############################################################

########################## Note ##############################
## Choose the working directory as the directory where 
## this code file locates.
##############################################################

########################## Note ##############################
## Please first install required R and Python packages.
## R: System_preparation.R
## Python: pip install -r requirements.txt
##############################################################



rm(list = ls())
filePath <- '../input_data/Raw_data_files/'

annotation = read.csv(file = paste0(filePath, 
                                    "GSE143791_cell.annotation.human.csv"), 
                      header = T)
cellNames = annotation$barcode
cellNames = gsub("-", "\\.", cellNames)
annotation$barcode = cellNames

fileName <- paste0(c('GSM4274678_BMET1-Tumor.count',
                     'GSM4274682_BMET10-Tumor.count',
                     'GSM4274683_BMET11-Tumor.count',
                     'GSM4274686_BMET2-Tumor.count',
                     'GSM4274689_BMET3-Tumor.count',
                     'GSM4274692_BMET5-Tumor.count',
                     'GSM4274697_BMET6-Tumor.count',
                     'GSM4274698_BMET7-Tumor.count',
                     'GSM4274701_BMET8-Tumor.count'), '.csv')

rawDataList <- Map(function(f){
  read.csv(file = paste0(filePath, f), header = T)
}, fileName)

rawDataAll = Reduce(cbind, Map(function(i){rawDataList[[i]][, -1]}, 2:9), rawDataList[[1]])
# dim: 32738 * 16994


## 1. Remove redundant genes.
geneNames = rawDataAll[, 1] 
uniGeneNames = unique(geneNames)  
redund_genes = Filter(function(x){sum(geneNames == x) > 1}, uniGeneNames)

rmInds = c()
for (i in seq_along(redund_genes)) {
  # red_gene = redund_genes[i]
  inds = which(geneNames == redund_genes[i])
  
  tmp = Filter(function(i){sum(rawDataAll[i, -1]) == 0}, inds)
  
  if (length(tmp) < length(inds) - 1) {
    rmInds = c(rmInds, inds)
  } else {
    rmInds = c(rmInds, tmp)
  }
}

rawDataAll_rmRedundGenes = rawDataAll[-rmInds, ]
# dim: 32602 * 16994


## 2. Remove external RNA controls consortium (ERCC) spike-in molecules.
geneNames_rmRedundGenes = rawDataAll_rmRedundGenes[, 1]
tmpERCC = grep("[Ee][Rr][Cc][Cc]", geneNames_rmRedundGenes)
rawDataAll_rmRedundGenesAndERCC = rawDataAll_rmRedundGenes[-tmpERCC, ]
# dim: 32591 * 16994


## 3. Remove genes with zero proportions greater than 95% across cells.
rawDataAll_2 = as.matrix(rawDataAll_rmRedundGenesAndERCC[, -1])
rownames(rawDataAll_2) <- rawDataAll_rmRedundGenesAndERCC[, 1]

p_genes = 0.05
numGenes = nrow(rawDataAll_2)
numCells = ncol(rawDataAll_2)

gene_proportion = unlist(Map(function(i){sum(rawDataAll_2[i, ] > 0) / numCells}, 
                             seq_len(numGenes)))

gene_rmInds = which(gene_proportion < p_genes)
rawDataAll_2_rmGenes = rawDataAll_2[-gene_rmInds, ]
# dim: 6626 * 16993


## 4. Normalization
rawDataAll_3 = rawDataAll_2_rmGenes
colSumsVal = colSums(rawDataAll_3)
s = 1e+6
rawDataAll_3_nm = log(t(apply(rawDataAll_3, 1, function(x){x / colSumsVal})) * s + 1)


## 5. Select the top 100 highly expressed genes.
numOfGenes = nrow(rawDataAll_3_nm)
gene_var = unlist(Map(function(i){var(rawDataAll_3_nm[i, ])}, seq_len(numOfGenes)))
numHEG = 100
rawDataAll_3_nm_HEG = rawDataAll_3_nm[order(gene_var, decreasing = T)[1:numHEG], ]
# dim: 100 * 16993


## 6. Utilize UMAP (McInnes et al., 2018) to reduce the data dimensionality to two.
library(umap)
settings = umap.defaults
settings$n_neighbors = 15
settings$min_dist = 0.1
settings$metric = "euclidean"
settings$random_state = 42

config_2D = settings
config_2D$n_components = 2

rawDataAll_4.umap2D = umap(t(rawDataAll_3_nm_HEG), config = config_2D)
rawDataAll_5.umap2D = rawDataAll_4.umap2D$layout
# dim: 16993    2


## 7. Add column of cell types 
dataFile = as.data.frame(rawDataAll_5.umap2D)
rowNames = rownames(dataFile)
dataFile$CellType = unlist(Map(function(x){annotation$cells[annotation$barcode == x]}, rowNames))
dataFile = dataFile[, c(ncol(dataFile), 1:(ncol(dataFile) - 1))]

# Input data
testData <- as.matrix(dataFile[, -1])
rownames(testData) <- NULL


## 9. Save the number of data for each patient
tmpNames = unlist(strsplit(rowNames, split = "\\."))[seq(1, 3*length(rowNames), by=3)]
uniNames = unique(tmpNames)

numOfSubjData = as.numeric(Map(function(x){sum(tmpNames == x)}, uniNames))


save(list = c("dataFile", "testData", "numOfSubjData", "annotation"), file = "Real_data.RData")
