library(data.table)
library(ggplot2)
library(qs)
library(cba)
library(ggpubr)
library(patchwork)

codeFolder = 'code'
dataFolder = 'data'
outputFolder = 'output'

source(file.path(codeFolder, 'cell_utils.R'))

glmnetCor = qread(file.path(dataFolder, 'glmnet_cor_dt.qs'))
glmnetGenes = unique(glmnetCor$gene1)

#monaco
cellDataMonaco = fread(file.path(dataFolder, 'rna_blood_cell_monaco.tsv.gz'))
cellDataMonaco = processCellData(cellDataMonaco, glmnetGenes)

monacoPlt =  plotCellData(cellDataMonaco)

#schmiedel
cellDataSchmiedel = fread(file.path(dataFolder, 'rna_blood_cell_schmiedel.tsv.gz'))
cellDataSchmiedel = processCellData(cellDataSchmiedel, glmnetGenes)

schmiedelPlt =  plotCellData(cellDataSchmiedel)


cellFig = monacoPlt / schmiedelPlt
cellFig = cellFig + plot_annotation(tag_levels = 'A')
ggexport(cellFig, filename = file.path(outputFolder, 'cell_heatmap.pdf')
  , width = 14, height = 18, units = 'in')
ggexport(cellFig, filename = file.path(outputFolder, 'cell_heatmap.png')
         , width = 1080, height = 720)

