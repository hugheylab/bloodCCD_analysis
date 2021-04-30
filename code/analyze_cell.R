library(data.table)
library(ggplot2)
library(qs)
library(cba)
library(ggpubr)

dataFolder = 'data'
outputFolder = 'output'

glmnetCor = qread(file.path(dataFolder, 'glmnet_cor_dt.qs'))
glmnetGenes = unique(glmnetCor$gene1)

#monaco
cellDataMonaco = fread(file.path(dataFolder, 'rna_blood_cell_monaco.tsv.gz'))
cellDataMonaco = processCellData(cellDataMonaco, glmnetGenes)

monacoPlt =  plotCellData(cellDataMonaco) +
  ggtitle('Gene expression by cell type, Monaco data')

#schmiedel
cellDataSchmiedel = fread(file.path(dataFolder, 'rna_blood_cell_schmiedel.tsv.gz'))
cellDataSchmiedel = processCellData(cellDataSchmiedel, glmnetGenes)

schmiedelPlt =  plotCellData(cellDataSchmiedel) +
  ggtitle('Gene expression by cell type, Schmiedel data')


cellFig = ggarrange(monacoPlt, schmiedelPlt, nrow = 2)
ggexport(cellFig, filename = file.path(outputFolder, 'cell_heatmap.pdf')
  , width = 14, height = 18, units = 'in')
