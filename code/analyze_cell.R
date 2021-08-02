source(file.path('code', 'utils.R'))

glmnetCor = qread(file.path(dataDir, 'glmnet_cor_dt.qs'))
glmnetGenes = unique(glmnetCor$gene1)

#monaco
cellDataMonaco = fread(file.path(dataDir, 'rna_blood_cell_monaco.tsv.gz'))
cellDataMonaco = processCellData(cellDataMonaco, glmnetGenes)
pMonaco = plotCellData(cellDataMonaco) +
  theme(axis.title.x = eb,
        axis.ticks.x = eb,
        axis.text.x = eb)

#schmiedel
cellDataSchmiedel = fread(file.path(dataDir, 'rna_blood_cell_schmiedel.tsv.gz'))
cellDataSchmiedel = processCellData(cellDataSchmiedel, glmnetGenes)
pSchmiedel = plotCellData(cellDataSchmiedel)

pCellData =  pMonaco + pSchmiedel +
  plot_layout(heights = c(1, 0.5)) +
  plot_annotation(tag_levels = 'A')
ggexport(pCellData, filename = file.path(outputDir, 'fig4.pdf'))

