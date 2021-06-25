getTimeCourseDt = function(emat, sm, genes) {
  
  timeCourseDt = as.data.table(t(emat[genes, sm$sample]), 
                               keep.rownames = 'sample')
  colnames(timeCourseDt)[-1] = as.character(
    lookUp(colnames(timeCourseDt)[-1], 'org.Hs.eg', 'SYMBOL', load = TRUE))
  
  timeCourseDt = merge(timeCourseDt, sm[, .(study, sample, ztFrac)], 
                       by = 'sample')
  timeCourseDt = melt(timeCourseDt, 
                      id.vars = c('sample', 'study', 'ztFrac'), 
                      value.name = 'expression', 
                      variable.name = 'gene')
  
  return(timeCourseDt)}

plotTimeCourse = function(timeCourseDt, ncol = NULL, nrow = NULL, breaks = 4) {
  
  p = ggplot(timeCourseDt, aes(x = ztFrac*24, y = expression, color = study)) +
    geom_point(shape = 1) +
    ylab('Expression (norm.)') +
    xlab('Time of day (h)') +
    facet_wrap(vars(gene), scales = 'free_y', ncol = ncol, nrow = nrow) +
    scale_x_continuous(breaks = seq(0, 24, by = breaks), limits = c(0, 24))
  
  return(p)}

getCormat = function(e, genes, entrezID = FALSE) {
  
  genes = unique(genes)
  
  if(class(e) == 'ExpressionSet') {
    emat = exprs(e)
  } else { emat = e }
  
  cormat = cor(t(emat)[, genes], method = 'spearman')
  
  if (!isTRUE(entrezID)) {
    while(isTRUE(e)) {
      geneNames = try(as.character(
        lookUp(rownames(cormat), 'org.Hs.eg', 'SYMBOL', load = TRUE)))
      if (class(geneNames) == 'try-error') {
        e = TRUE
      } else {
        e = FALSE
        rownames(cormat) = geneNames}}
      
    colnames(cormat) = rownames(cormat)}
  
  return(cormat)}

sortCormat = function (cormat) {
  
  distmat = as.dist(1 - cormat)/2
  hc = hclust(distmat)$merge
  opt = order.optimal(distmat, hc)$order
  ord = unique(colnames(cormat[opt, opt]))
  
  cormatDt = as.data.table(cormat, 
                           keep.rownames = 'gene1')
  cormatDt = melt(cormatDt, 
                  variable.name = 'gene2', 
                  value.name = 'rho')
    
  cormatDt[, `:=`(gene1 = factor(gene1, levels = ord), 
                  gene2 = factor(gene2, levels = ord))]
  
  cormatDt[gene1 == gene2, 
           rho := NA]
  
  return(cormatDt)}

plotHeatmap = function (cormatDt, ..., ncol = NULL, nrow = NULL, scales = 'free') {
  
  hm = ggplot(cormatDt, aes(x = gene1, y = gene2,  fill = rho))+
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
      midpoint = 0, limit = c(-1,1), space = 'Lab', name = 'rho') +
    xlab('Gene') +
    ylab('Gene')
  
  if (!missing(...)) { 
    hm = hm + facet_wrap(vars(...), scales = scales, ncol = ncol, nrow = nrow)}
  
  return(hm)}

calcCCD = function(eset1, eset2, genes, scale = TRUE) {
  
  corMat1 = getCormat(eset1, genes)
  corMat2 = getCormat(eset2, genes)
  
  corVec1 = corMat1[upper.tri(corMat1)]
  corVec2 = corMat2[upper.tri(corMat2)]
  
  ccd = as.numeric(dist(rbind(corVec1, corVec2), method = 'Euclidean'))
  
  if (scale) {
    nPairs = choose(length(genes), 2)
    ccd = ccd/nPairs}
  
  return(ccd)}
