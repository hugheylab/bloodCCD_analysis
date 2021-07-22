library('annotate')
library('Biobase')
library('BiocParallel')
library('cba')
library('data.table')
library('deltaccd')
library('doParallel')
library('ggplot2')
library('ggpubr')
library('ggrepel')
library('glmnet')
library('glue')
library('lubridate')
library('limorhyde2')
library('metapredict')
library('patchwork')
library('qs')
library('RColorBrewer')
library('sva')
library('VennDiagram')
library('zeitzeiger')

theme_set(theme_bw())
registerDoParallel(cores = 2)

#setup vars
codeFolder = file.path('code')
outputFolder = file.path('output')
dataFolder = file.path('data')

studyMetadataPath = file.path(dataFolder, 'metadata', 'study_metadata.csv')
sampleMetadataPath = file.path(dataFolder, 'metadata', 'sample_metadata.csv')
ematPath = file.path(dataFolder, 'circadian_human_blood_emat.qs')
esetPath = file.path(dataFolder, 'circadian_human_blood.qs')
ematPerturbPath = file.path('data', 'perturb_emat.qs')

#ggplot convenience vars
eb = element_blank()

#cv utils
plotCoefs = function(coefDt, ncol = NULL, nrow = NULL, ...) {
  
  p = ggplot(coefDt) +
    geom_point(aes(x = reorder(gene_sym, coef), y = coef), size = 1) +
    geom_segment(aes(x = reorder(gene_sym, coef), xend = reorder(gene_sym, coef), 
                     y = 0, yend = coef)) +
    scale_y_continuous(expand = c(0.015, 0)) +
    coord_flip() +
    labs(x = 'Gene', y = 'Coefficient')
  
  if (!missing(...)) {
    p = p + facet_wrap(vars(...), scales = 'free_y', ncol = ncol, nrow = nrow)}
  
  return(p)}


convertZt = function(md) {
  md[, zt := as.duration(hm(clock_time) - hm(sunrise_time))/as.duration(hours())
     ][zt < 0, zt := zt + timeMax]
  md[, ztFrac := zt/timeMax
  ][, zt := NULL]
  
  return(md)}


#cor utils
getTimeCourseDt = function(emat, sm, genes) {
  
  timeCourseDt = as.data.table(t(emat[genes, sm$sample]), 
                               keep.rownames = 'sample')
  setnames(timeCourseDt, 2:length(colnames(timeCourseDt)), 
           lookup(colnames(timeCourseDt)[-1], 'org.Hs.eg', 'SYMBOL', 
                  load = TRUE))
  
  timeCourseDt = merge(timeCourseDt, sm[, .(study, sample, ztFrac)], 
                       by = 'sample')
  timeCourseDt = melt(timeCourseDt, 
                      id.vars = c('sample', 'study', 'ztFrac'), 
                      value.name = 'expression', 
                      variable.name = 'gene')
  
  return(timeCourseDt)}


plotTimeCourse = function(timeCourseDt, ncol = NULL, nrow = NULL, breaks = 4) {
  
  p = ggplot(timeCourseDt) +
    geom_point(aes(x = ztFrac*24, y = expression, color = study), shape = 1) +
    labs(x = 'Expression (norm.)', y = 'Time of day (h)') +
    facet_wrap(vars(gene), scales = 'free_y', ncol = ncol, nrow = nrow) +
    scale_x_continuous(breaks = seq(0, 24, by = breaks), limits = c(0, 24))
  
  return(p)}


getCormat = function(e, genes, entrezID = FALSE) {

  genes = unique(genes)
 
  if(inherits(e, 'ExpressionSet')) {
    emat = exprs(e)
  } else { emat = e }
  
  cormat = cor(t(emat)[, genes], method = 'spearman')
  
  geneNames = as.character(
        lookUp(rownames(cormat), 'org.Hs.eg', 'SYMBOL', load = TRUE))
  rownames(cormat) = geneNames
  colnames(cormat) = rownames(cormat)
  
  return(cormat)}


sortCormat = function (cormat) {
  
  distmat = as.dist(1 - cormat)/2
  hc = hclust(distmat)$merge
  opt = order.optimal(distmat, hc)$order
  ord = unique(colnames(cormat[opt, opt]))
  
  cormatDt = as.data.table(cormat, keep.rownames = 'gene1')
  cormatDt = melt(cormatDt, id.vars = 'gene1', 
                  variable.name = 'gene2', value.name = 'rho')
    
  cormatDt[, gene1 := factor(gene1, levels = ord)]
  cormatDt[, gene2 := factor(gene2, levels = rev(ord))]
  
  cormatDt[gene1 == gene2, rho := NA]
  
  return(cormatDt)}


plotHeatmap = function (cormatDt, ..., ncol = NULL, nrow = NULL, scales = 'free') {
  
  hm = ggplot(cormatDt) +
    geom_tile(aes(x = gene1, y = gene2,  fill = rho), color = 'white') +
    scale_fill_distiller(palette = 'PuOr', limits = c(-1,1)) +
    labs(x = 'Gene', y = 'Gene', fill = 'rho')
  
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


#cell_utils
processCellData = function(cellData, genes) {
  
  setnames(cellData, 1:3, c('probe', 'gene', 'cell'))

  cellDataScaled = cellData[gene %in% genes, .(gene, cell, tpm = TPM)]
  cellDataScaled[, tpm := (tpm - mean(tpm))/sd(tpm), by = gene]
  cellDataScaled[is.nan(tpm), tpm := 0]
  cellDataScaled[, gene := factor(gene, levels = levels(genes))]

  cellDataWide = dcast(cellDataScaled, cell ~ gene, value.var = 'tpm')
  d = dist(as.matrix(cellDataWide, rownames = 'cell'))
  hc = hclust(d)$merge
  opt = order.optimal(d, hc)$order

  cellDataScaled[, cell := factor(cell, levels = names(opt)[opt])]
  
  return(cellDataScaled)}


plotCellData = function(cellData, ..., scales = NULL, ncol = NULL,
                        nrow = NULL, drop = TRUE) {
  
  p = ggplot(cellData) +
    geom_tile(aes(x = gene, y = cell, fill = tpm), color = 'white') +
    scale_fill_distiller(palette = 'PuOr') +
    labs(x = 'Gene', y = 'Cell type', fill = 'scaled tpm') +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
  
  if (!missing(...)) {
    p = p + facet_wrap(vars(...), scales = scales, nrow = nrow, ncol = ncol,
                       drop = drop)}
  
  return(p)}
