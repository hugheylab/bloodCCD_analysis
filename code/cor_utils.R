getCormat = function(e, genes) {
  
  genes = unique(genes)
  
  if(class(e) == 'ExpressionSet') {
    emat = exprs(e)
  } else { emat = e }
  
  rownames(emat) = as.character(
    lookUp(rownames(emat), 'org.Hs.eg', 'SYMBOL', load = TRUE)
    )
  cormat = cor(t(emat)[, genes], method = 'spearman')
  
  return(cormat)
  }

sortCormat = function (cormat) {
  
  distmat = as.dist(1 - cormat)/2
  hc = hclust(distmat)$merge
  opt = order.optimal(distmat, hc)$order
  ord = unique(colnames(cormat[opt, opt]))
  
  cormatDt = as.data.table(cormat, keep.rownames = 'gene1')
  cormatDt = melt(cormatDt, variable.name = 'gene2'
    , value.name = 'rho')
    
  cormatDt[
    , `:=`(gene1 = factor(gene1, levels = ord)
           , gene2 = factor(gene2, levels = ord)
           )
    ]
  
  cormatDt[gene1 == gene2
    , rho := NA
    ]
  
  return(cormatDt)
  }

plotHeatmap = function (cormatDt, ...) {
  
  hm = ggplot(cormatDt, aes(x = gene1, y = gene2,  fill = rho))+
    geom_tile(color = 'white') +
    scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
      midpoint = 0, limit = c(-1,1), space = 'Lab', name = 'rho') +
    xlab('Gene') +
    ylab('Gene')
  
  if (!missing(...)) { 
    hm = hm + facet_wrap(vars(...), scales = 'free')
    }
  
  return(hm)
  }

calcCCD = function(eset1, eset2, genes, scale = TRUE) {
  
  corMat1 = getCormat(eset1, genes)
  corMat2 = getCormat(eset2, genes)
  
  corVec1 = corMat1[upper.tri(corMat1)]
  corVec2 = corMat2[upper.tri(corMat2)]
  
  ccd = as.numeric(dist(rbind(corVec1, corVec2), method = 'Euclidean'))
  
  if (scale) {
    
    nPairs = choose(length(genes), 2)
    ccd = ccd/nPairs
    
    }
  
  return(ccd)
  }
