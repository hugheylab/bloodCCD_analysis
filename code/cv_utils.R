library(data.table)
library(Biobase)
library(annotate)

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
  
  return(timeCourseDt)
  
  }

plotTimeCourse = function(timeCourseDt) {
  
  p = ggplot(timeCourseDt, aes(x = ztFrac*24, y = expression, color = study)) +
    geom_point(shape = 1) +
    ylab('Expression') +
    xlab('Hour of day') +
    facet_wrap(~ gene, scales = 'free_y') +
    scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24))
  
  return(p)
  
  }

plotCoefs = function(coefDt, ..., ncol = NULL, nrow = NULL) {
  
  p = ggplot(data = coefDt
             , aes(x = reorder(gene_sym, coef), y = coef)) +
    
    geom_point(size = 1) +
    geom_segment(aes(x = reorder(gene_sym, coef), xend = reorder(gene_sym, coef)
                     , y = 0, yend = coef)) +
    scale_y_continuous(expand = c(0.015, 0)) +
    coord_flip() +
    xlab('Gene') +
    ylab('Coefficient') +
    theme(axis.text.y = element_text(size = 8))
  
  if (!missing(...)) {
    p = p + 
     facet_wrap(vars(...), scales = 'free_y', ncol = ncol, nrow = nrow)
    
    }
  
  return(p)
  
  }
