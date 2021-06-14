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
pMonaco = plotCellData(cellDataMonaco) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

#schmiedel
cellDataSchmiedel = fread(file.path(dataFolder, 'rna_blood_cell_schmiedel.tsv.gz'))
cellDataSchmiedel = processCellData(cellDataSchmiedel, glmnetGenes)
pSchmiedel = plotCellData(cellDataSchmiedel)

# #y axis
# yax = ggplot(data.frame(l = 'Cell type', x = 1, y = 1)) +
#   geom_text(aes(x, y, label = l), angle = 90, size = 8) + 
#   theme_void() +
#   coord_cartesian(clip = "off")

pCellData =  pMonaco + pSchmiedel +
  plot_layout(heights = c(1, .5)) +
  plot_annotation(tag_levels = 'A')
ggexport(pCellData, filename = file.path(outputFolder, 'cell_heatmap.pdf')
  , width = 14, height = 18, units = 'in')
ggexport(pCellData, filename = file.path(outputFolder, 'fig4.pdf')
         , width = 18, height = 18)

