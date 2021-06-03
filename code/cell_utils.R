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
  
  return(cellDataScaled)}

plotCellData = function(cellData, ..., scales = NULL, ncol = NULL,
                        nrow = NULL) {
  
  p = ggplot(cellData, aes(x = gene, y = cell, fill = tpm)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white' 
                         , midpoint = 0, space = 'Lab', name = 'scaled tpm') +
    xlab('Gene') +
    ylab('Cell type') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if (!missing(...)) {
    p = p + facet_wrap(vars(...), scales = scales, nrow = nrow, ncol = ncol)}
  
  return(p)}
