plotCoefs = function(coefDt, ..., ncol = NULL, nrow = NULL) {
  
  p = ggplot(data = coefDt, aes(x = reorder(gene_sym, coef), y = coef)) +
    geom_point(size = 1) +
    geom_segment(aes(x = reorder(gene_sym, coef), xend = reorder(gene_sym, coef)
                     , y = 0, yend = coef)) +
    scale_y_continuous(expand = c(0.015, 0)) +
    coord_flip() +
    xlab('Gene') +
    ylab('Coefficient') +
  
  if (!missing(...)) {
    p = p + 
      facet_wrap(vars(...), scales = 'free_y', ncol = ncol, nrow = nrow)}
  
  return(p)}
