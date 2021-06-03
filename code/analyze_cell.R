library(data.table)
library(ggplot2)
library(qs)
library(cba)
library(ggpubr)
library(patchwork)
theme_set(theme_bw(base_size = 20))

codeFolder = 'code'
dataFolder = 'data'
outputFolder = 'output'

source(file.path(codeFolder, 'cell_utils.R'))

glmnetCor = qread(file.path(dataFolder, 'glmnet_cor_dt.qs'))
glmnetGenes = unique(glmnetCor$gene1)

#monaco
cellDataMonaco = fread(file.path(dataFolder, 'rna_blood_cell_monaco.tsv.gz'))
cellDataMonaco = processCellData(cellDataMonaco, glmnetGenes)
cellDataMonaco[, source := 'Monaco']

#schmiedel
cellDataSchmiedel = fread(file.path(dataFolder, 'rna_blood_cell_schmiedel.tsv.gz'))
cellDataSchmiedel = processCellData(cellDataSchmiedel, glmnetGenes)
cellDataSchmiedel[, source := 'Schmiedel']

cellData = rbind(cellDataMonaco, cellDataSchmiedel)
pCellData = plotCellData(cellData, source, scales = 'free_y', ncol = 1)

ggexport(pCellData, filename = file.path(outputFolder, 'cell_heatmap.pdf')
  , width = 14, height = 18, units = 'in')
ggexport(pCellData, filename = file.path(outputFolder, 'fig4.png')
         , width = 1400, height = 1080)

