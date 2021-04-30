library(data.table)
library(Biobase)
library(annotate)

getResults()

getTimeCourseDt = function(emat, sm, genes) {
  
  timeCourseDt = as.data.table(t(emat[genes, sm$sample])
    , keep.rownames = 'sample')
  colnames(timeCourseDt)[-1] = as.character(
    lookUp(colnames(timeCourseDt)[-1], 'org.Hs.eg', 'SYMBOL', load = TRUE))
  
  timeCourseDt = merge(timeCourseDt, sm[, .(study, sample, ztFrac)]
    , by = 'sample')
  timeCourseDt = melt(timeCourseDt
    , id.vars = c('sample', 'study', 'ztFrac'), value.name = 'expression'
    , variable.name = 'gene')
  
  return(timeCourseDt)}
