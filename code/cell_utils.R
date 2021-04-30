processCellData = function(cellData, genes) {
  
  setnames(cellData, names(cellData)[1:3], c('probe', 'gene', 'cell'))

  cellDataScaled = cellData[
    gene %in% genes 
    , .(gene, cell, tpm = TPM)]
  cellDataScaled[
    , tpm := (tpm - mean(tpm))/sd(tpm)
    , by = gene]
  cellDataScaled[is.nan(tpm)
    , tpm := 0]
  cellDataScaled[, gene := factor(gene, levels = attributes(genes)$levels)]

  cellDataWide = dcast(cellDataScaled, cell ~ gene, value.var = 'tpm')
  d = dist(as.matrix(cellDataWide, rownames = 'cell'))
  hc = hclust(d)$merge
  opt = order.optimal(d, hc)$order

  cellDataScaled[, cell := factor(cell, levels = names(opt)[opt])]
  
  return(cellDataScaled)
  
}

plotCellData = function(cellData) {
  
  p = ggplot(cellData, aes(x = cell, y = gene, fill = tpm)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white' 
                         , midpoint = 0, space = 'Lab', name = 'scaled tpm') +
    xlab('Cell') +
    ylab('Gene') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)
          , axis.text.y = element_text(size = 6))
  
}
