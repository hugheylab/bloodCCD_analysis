library(zeitzeiger)
library(glmnet)
library(Biobase)
library(doParallel)
registerDoParallel(cores = 2)
library(data.table)
library(qs)
library(lubridate)
library(glue)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(annotate)
library(ggpubr)
library(cba)
library(VennDiagram)
library(RColorBrewer)


getCormat = function(eset, genes) {
  
  emat = exprs(eset)
  rownames(emat) = as.character(lookUp(rownames(emat)
                                       , 'org.Hs.eg', 'SYMBOL', load = TRUE))
  cormat = cor(t(emat)[, genes], method = 'spearman')
  
  return(cormat)}

calcCCD = function(eset1, eset2, genes, scale = TRUE) {
  
  corMat1 = getCormat(eset1, genes)
  corMat2 = getCormat(eset2, genes)
  
  corVec1 = corMat1[upper.tri(corMat1)]
  corVec2 = corMat2[upper.tri(corMat2)]
  
  ccd = as.numeric(dist(rbind(corVec1, corVec2), method = 'Euclidean'))
  
  if (scale) {ccd = ccd/length(genes)}
  
  return(ccd)}


theme_set(theme_bw())

outputFolder = file.path('output')
dataFolder = file.path('data')

studyMetadataPath = file.path('data', 'metadata', 'study_metadata.csv')
studyMetadata = fread(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = file.path('data', 'metadata', 'sample_metadata.csv')
sampleMetadata = fread(sampleMetadataPath, stringsAsFactors = FALSE)

ematPath = file.path('data', 'circadian_human_blood_emat.qs')
emat = qread(ematPath)

timeMax = 24
sampleMetadata[
  , zt := as.duration(hm(clock_time) - hm(sunrise_time))/as.duration(hours(1))
  ][zt < 0, zt := zt + timeMax]
sampleMetadata[
  , ztFrac := zt/timeMax
  ][, zt := NULL]
controlConds = c('Sleep Extension', 'In phase with respect to melatonin'
  , 'baseline')
controlMetadata = sampleMetadata[condition %in% controlConds]

perturbConds = c('Sleep Restriction', 'Out of phase with respect to melatonin'
  , 'sleep deprivation')
perturbMetadata = sampleMetadata[condition %in% perturbConds]


zzCoefs = qread(file.path(dataFolder, 'zeitzeiger_coefs.qs'))
setorderv(zzCoefs, c('sumabsv', 'spc', 'coeff'))

glmnetCoefs = qread(file.path(dataFolder, 'glmnet_coefs.qs'))
setorderv(glmnetCoefs, c('lambda', 'coef'))

genes2017Coefs = qread(file.path(dataFolder, 'genes2017.qs'))
setorderv(genes2017Coefs, c('coef'))

### merged data
#zeitzeiger
finalSumAbsv = 2:3
zzCormatMelt = foreach(absv = finalSumAbsv, .combine = rbind) %dopar% {
  
  vCo = zzCoefs[sumabsv == absv]
  
  zCormat = cor(t(emat)[, vCo$gene], method = 'spearman')
  colnames(zCormat) = as.character(lookUp(as.character(colnames(zCormat))
                                          , 'org.Hs.eg', 'SYMBOL', load=TRUE))
  rownames(zCormat) = as.character(lookUp(as.character(rownames(zCormat))
                                          , 'org.Hs.eg', 'SYMBOL', load=TRUE))
  
  zDist = as.dist(1 - zCormat)/2
  zHc = hclust(zDist)$merge
  zOpt = order.optimal(zDist, zHc)$order
  zCormatOrd = zCormat[zOpt, zOpt]
  
  zCormatDt = as.data.table(zCormat, keep.rownames = 'gene1')
  zCormatMelt = melt(zCormatDt, variable.name = 'gene2'
                     , value.name = 'rho')
  zCormatMelt[
    , `:=`(sumabsv = absv
           , gene1 = factor(gene1, levels = unique(colnames(zCormatOrd)))
           , gene2 = factor(gene2, levels = unique(colnames(zCormatOrd))))]
  
  return(zCormatMelt)}
zzCormatMelt[gene1 == gene2, rho := NA]

zzHeatmap =  ggplot(zzCormatMelt
                    , aes(x = gene1, y = gene2,  fill = rho))+
  geom_tile(color = 'white') +
  facet_wrap(~ sumabsv, scales = 'free') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
                       midpoint = 0, limit = c(-1,1), space = 'Lab', 
                       name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  ggtitle(glue('Overall correlation heatmaps, genes selected by zeitzeiger'
               , ', by sumabsv with nSpc = 2')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 16))


zzCormatCoefMelt = foreach(absv = finalSumAbsv, .combine = rbind) %dopar% {
  
  vCo = zzCoefs[sumabsv == absv 
    & coeff != 0]
  setorder(vCo, coeff)
  zOrd = 1:nrow(vCo)
  names(zOrd) = vCo$gene_sym 
  
  zCormat = cor(t(emat)[, vCo$gene], method = 'spearman')
  colnames(zCormat) = as.character(lookUp(as.character(colnames(zCormat))
                                          , 'org.Hs.eg', 'SYMBOL', load=TRUE))
  rownames(zCormat) = as.character(lookUp(as.character(rownames(zCormat))
                                          , 'org.Hs.eg', 'SYMBOL', load=TRUE))
  
  zCormatOrd = zCormat[zOrd, zOrd]
  
  zCormatDt = as.data.table(zCormat, keep.rownames = 'gene1')
  zCormatMelt = melt(zCormatDt, variable.name = 'gene2'
                     , value.name = 'rho')
  zCormatMelt[
    , `:=`(sumabsv = absv
           , gene1 = factor(gene1, levels = unique(colnames(zCormatOrd)))
           , gene2 = factor(gene2, levels = unique(colnames(zCormatOrd))))]
  
  return(zCormatMelt)}
zzCormatCoefMelt[gene1 == gene2, rho := NA]

zzCoefHeatmap =  ggplot(zzCormatCoefMelt
                    , aes(x = gene1, y = gene2,  fill = rho))+
  geom_tile(color = 'white') +
  facet_wrap(~ sumabsv, scales = 'free') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
                       midpoint = 0, limit = c(-1,1), space = 'Lab', 
                       name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  ggtitle(glue('Overall correlation heatmaps, genes selected by zeitzeiger'
               , ', by sumabsv with nSpc = 2')
          , subtitle = 'ordered by coefficients') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 16))

#glmnet
glmnetCormatMelt = foreach(lam = unique(glmnetCoefs$lambda)
  , .combine = rbind) %dopar% {
    
    geneSummGnet = glmnetCoefs[lambda == lam]
                             
    gnetCormat = cor(t(emat)[, geneSummGnet$gene], method = 'spearman')
    colnames(gnetCormat) = as.character(lookUp(as.character(colnames(gnetCormat))
      , 'org.Hs.eg', 'SYMBOL', load=TRUE))
    rownames(gnetCormat) = as.character(lookUp(as.character(rownames(gnetCormat))
      , 'org.Hs.eg', 'SYMBOL', load=TRUE))
                             
    gnetDist = as.dist(1 - gnetCormat)/2
    gnetHc = hclust(gnetDist)$merge
    gnetOpt = order.optimal(gnetDist, gnetHc)$order
    gnetCormatOrd = gnetCormat[gnetOpt, gnetOpt]
                             
    gnetCormatDt = as.data.table(gnetCormat, keep.rownames = 'gene1')
    gnetCormatMelt = melt(gnetCormatDt, variable.name = 'gene2'
      , value.name = 'rho')
    gnetCormatMelt[
      , `:=`(lambda = lam
            , gene1 = factor(gene1, levels = unique(colnames(gnetCormatOrd)))
            , gene2 = factor(gene2, levels = unique(colnames(gnetCormatOrd))))]
                             
    return(gnetCormatMelt)}
glmnetCormatMelt[gene1 == gene2, rho := NA]

pGlmnetHeatmap =  ggplot(glmnetCormatMelt
  , aes(x = gene1, y = gene2,  fill = rho)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
                       midpoint = 0, limit = c(-1,1), space = 'Lab', 
                       name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  facet_wrap(~ factor(lambda, levels = unique(glmnetCoefs$lambda))
    , scales = 'free') +
  ggtitle(bquote('Overall correlations, genes selected by glmnet, by '~lambda)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 16))


glmnetCormatCoefMelt = foreach(lam = unique(glmnetCoefs$lambda)
  , .combine = rbind) %dopar% {
  
  geneSummGnet = glmnetCoefs[lambda == lam]
  setorder(geneSummGnet, coef)
  gnetOrd = 1:nrow(geneSummGnet)
  names(gnetOrd) = geneSummGnet$gene_sym 
  
  gnetCormat = cor(t(emat)[, geneSummGnet$gene], method = 'spearman')
  colnames(gnetCormat) = as.character(lookUp(as.character(colnames(gnetCormat))
                                             , 'org.Hs.eg', 'SYMBOL', load=TRUE))
  rownames(gnetCormat) = as.character(lookUp(as.character(rownames(gnetCormat))
                                             , 'org.Hs.eg', 'SYMBOL', load=TRUE))
  
  gnetCormatOrd = gnetCormat[gnetOrd, gnetOrd]
  
  gnetCormatDt = as.data.table(gnetCormat, keep.rownames = 'gene1')
  gnetCormatMelt = melt(gnetCormatDt, variable.name = 'gene2'
                        , value.name = 'rho')
  gnetCormatMelt[
    , `:=`(lambda = lam
           , gene1 = factor(gene1, levels = unique(colnames(gnetCormatOrd)))
           , gene2 = factor(gene2, levels = unique(colnames(gnetCormatOrd))))]
  
  return(gnetCormatMelt)}
glmnetCormatCoefMelt[gene1 == gene2, rho := NA]

pGlmnetCoefHeatmap =  ggplot(glmnetCormatCoefMelt
  , aes(x = gene1, y = gene2,  fill = rho)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
                       midpoint = 0, limit = c(-1,1), space = 'Lab', 
                       name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  facet_wrap(~ factor(lambda, levels = unique(glmnetCoefs$lambda))
    , scales = 'free') +
  ggtitle(bquote('Overall correlations, genes selected by glmnet, by '~lambda)
    , subtitle = 'ordered by coefficients') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 16))

#2017
emat2017 = emat
rownames(emat2017) = as.character(lookUp(rownames(emat2017)
  , 'org.Hs.eg', 'SYMBOL', load = TRUE))

cormat2017 = cor(t(emat2017)[, unique(genes2017Coefs$gene_sym)]
  , method = 'spearman')
dist2017 = as.dist(1 - cormat2017)/2
hc2017 = hclust(dist2017)$merge
ord2017 = order.optimal(dist2017, hc2017)$order
cormatOrd2017 = cormat2017[ord2017, ord2017]

cormatDt2017 = as.data.table(cormat2017, keep.rownames = 'gene1')
cormatMelt2017 = melt(cormatDt2017, variable.name = 'gene2'
                      , value.name = 'rho')
cormatMelt2017[
  , `:=`(gene1 = factor(gene1, levels = unique(colnames(cormatOrd2017)))
         , gene2 = factor(gene2, levels = unique(colnames(cormatOrd2017))))]
cormatMelt2017[gene1 == gene2, rho := NA]

heatmap2017 =  ggplot(cormatMelt2017, aes(x = gene1, y = gene2,  fill = rho)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white' 
    , midpoint = 0, limit = c(-1,1), space = 'Lab', name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  ggtitle('Overall correlation matrix, genes selected by zeitzeiger, 2017') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 16))
heatmap2017 = plot_grid(heatmap2017, NULL, ncol = 2)

ord2017Coef = 1:nrow(genes2017Coefs[!is.na(coef)])
names(ord2017Coef) = genes2017Coefs[!is.na(coef)]$gene_sym
cormatOrd2017Coef = cormat2017[ord2017Coef, ord2017Coef]
cormatDt2017Coef = as.data.table(cormat2017, keep.rownames = 'gene1')
cormatMelt2017Coef = melt(cormatDt2017Coef, variable.name = 'gene2'
                      , value.name = 'rho')
cormatMelt2017Coef[
  , `:=`(gene1 = factor(gene1, levels = unique(colnames(cormatOrd2017Coef)))
         , gene2 = factor(gene2, levels = unique(colnames(cormatOrd2017Coef))))]
cormatMelt2017Coef[gene1 == gene2, rho := NA]

heatmap2017Coef =  ggplot(cormatMelt2017Coef, aes(x = gene1, y = gene2
  , fill = rho)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white' 
    , midpoint = 0, limit = c(-1,1), space = 'Lab', name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  ggtitle('Overall correlation matrix, genes selected by zeitzeiger, 2017'
    , subtitle = 'ordered by coefficients') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 16))
heatmap2017Coef = plot_grid(heatmap2017Coef, NULL, ncol = 2)

mergedFig = ggarrange(plotlist = list(zzHeatmap, pGlmnetHeatmap, heatmap2017)
  , nrow = 3)
ggexport(mergedFig, filename = file.path(outputFolder, 'overall_correlations.pdf')
  , width = 36, height = 36, unit = 'in')

mergedFigCoef = ggarrange(plotlist = list(zzCoefHeatmap, pGlmnetCoefHeatmap
  , heatmap2017Coef), nrow = 3)
ggexport(mergedFigCoef
         , filename = file.path(outputFolder, 'overall_correlations_coef.pdf')
         , width = 36, height = 36, unit = 'in')

###study data-sets
studyEsetList = qread(file.path(dataFolder, 'subj_norm_esetList.qs'))
studyEsetListPerturb = qread(file.path(dataFolder, 'subj_norm__pert_esetList.qs'))

#zeitzeiger
zzStudyCormat = foreach(i = 1:length(studyEsetList)
  , .combine = rbind) %dopar% {
              
  eset = studyEsetList[[i]]  
    
  ematTmp = exprs(eset)  
  
  foreach(absv = unique(zzCoefs$sumabsv), .combine = rbind) %dopar% {
                   
    vCo = zzCoefs[sumabsv == absv]
  
    zCormat = cor(t(ematTmp)[, vCo$gene], method = 'spearman')
    colnames(zCormat) = as.character(lookUp(as.character(colnames(zCormat))
      , 'org.Hs.eg', 'SYMBOL', load=TRUE))
    rownames(zCormat) = as.character(lookUp(as.character(rownames(zCormat))
      , 'org.Hs.eg', 'SYMBOL', load=TRUE))
  
    zDist = as.dist(1 - zCormat)/2
    zHc = hclust(zDist)$merge
    zOpt = order.optimal(zDist, zHc)$order
    zCormatOrd = zCormat[zOpt, zOpt]
  
    zCormatDt = as.data.table(zCormat, keep.rownames = 'gene1')
    zCormatMelt = melt(zCormatDt, variable.name = 'gene2'
      , value.name = 'rho')
    zCormatMelt[
      , `:=`(study = names(studyEsetList)[i]
        , sumabsv = absv
        , gene1 = factor(gene1, levels = unique(colnames(zCormatOrd)))
        , gene2 = factor(gene2, levels = unique(colnames(zCormatOrd))))]
  
    return(zCormatMelt)}}
zzStudyCormat[
  gene1 == gene2
  , rho := NA]

zzStudyHeatmap =  ggplot(zzStudyCormat
  , aes(x = gene1, y = gene2,  fill = rho)) +
  geom_tile(color = 'white') +
  facet_wrap(~ sumabsv + study, scales = 'free', nrow = 2) +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
                       midpoint = 0, limit = c(-1,1), space = 'Lab', 
                       name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  ggtitle(glue('Correlation heatmaps, genes selected by zeitzeiger'
               , ', by sumabsv with nSpc = 2')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 20))


#glmnet
glmnetStudyCormat = foreach(i = 1:length(studyEsetList)
  , .combine = rbind) %dopar% {
              
  eset = studyEsetList[[i]]
    
  ematTmp = exprs(eset)  
  
  foreach(lam = unique(glmnetCoefs$lambda), .combine = rbind) %dopar% {
                   
    geneSummGnet = glmnetCoefs[lambda == lam]
                             
    gnetCormat = cor(t(ematTmp)[, geneSummGnet$gene]
      , method = 'spearman')
    colnames(gnetCormat) = as.character(lookUp(as.character(colnames(gnetCormat))
      , 'org.Hs.eg', 'SYMBOL', load=TRUE))
    rownames(gnetCormat) = as.character(lookUp(as.character(rownames(gnetCormat))
      , 'org.Hs.eg', 'SYMBOL', load=TRUE))
                             
    gnetDist = as.dist(1 - gnetCormat)/2
    gnetHc = hclust(gnetDist)$merge
    gnetOpt = order.optimal(gnetDist, gnetHc)$order
    gnetCormatOrd = gnetCormat[gnetOpt, gnetOpt]
                             
    gnetCormatDt = as.data.table(gnetCormat, keep.rownames = 'gene1')
    gnetCormatMelt = melt(gnetCormatDt, variable.name = 'gene2'
    , value.name = 'rho')
    gnetCormatMelt[
      , `:=`(study = names(studyEsetList)[i]
        , lambda = lam
        , gene1 = factor(gene1, levels = unique(colnames(gnetCormatOrd)))
        , gene2 = factor(gene2, levels = unique(colnames(gnetCormatOrd))))]
                             
    return(gnetCormatMelt)}}
glmnetStudyCormat[
  gene1 == gene2
  , rho := NA]


glmnetStudyHeatmap =  ggplot(glmnetStudyCormat
  , aes(x = gene1, y = gene2,  fill = rho)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white'
    , midpoint = 0, limit = c(-1,1), space = 'Lab',  name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  facet_wrap(~ lambda + study
    , scales = 'free') +
  ggtitle(bquote('Correlations by study, genes selected by glmnet, by '~lambda)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 20))

#2017
study2017Cormat = foreach(i = 1:length(studyEsetList)
  , .combine = rbind) %dopar% {
                              
  eset = studyEsetList[[i]]  
                              
  ematTmp = exprs(eset)
  rownames(ematTmp) = as.character(lookUp(rownames(ematTmp)
    , 'org.Hs.eg', 'SYMBOL', load = TRUE))
                              
  cormatTmp = cor(t(ematTmp)[, unique(genes2017)], method = 'spearman')
  distTmp = as.dist(1 - cormatTmp)/2
  hcTmp = hclust(distTmp)$merge
  ordTmp = order.optimal(distTmp, hcTmp)$order
  cormatOrdTmp = cormatTmp[ordTmp, ordTmp]
  
  cormatDtTmp = as.data.table(cormatTmp, keep.rownames = 'gene1')
  cormatMeltTmp = melt(cormatDtTmp, variable.name = 'gene2'
                        , value.name = 'rho')
  cormatMeltTmp[
    , `:=`(study = names(studyEsetList)[i]
      , gene1 = factor(gene1, levels = unique(colnames(cormatOrdTmp)))
      , gene2 = factor(gene2, levels = unique(colnames(cormatOrdTmp))))]
  
  
  return(cormatMeltTmp)}
study2017Cormat[
  gene1 == gene2
  , rho := NA]

heatmapStudy2017 =  ggplot(study2017Cormat, aes(x = gene1, y = gene2
  ,  fill = rho)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white' 
    , midpoint = 0, limit = c(-1,1), space = 'Lab', name = 'rho') +
  facet_wrap(~ study, scales = 'free') +
  xlab('Gene') +
  ylab('Gene') +
  ggtitle('Correlations by study, genes selected by zeitzeiger, 2017') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 20))
# heatmapStudy2017 = plot_grid(heatmapStudy2017, NULL, NULL, ncol = 3)

studyFig = ggarrange(plotlist = list(zzStudyHeatmap, glmnetStudyHeatmap
  , heatmapStudy2017), nrow = 3, heights = c(2, 2, 1.1))
ggexport(studyFig, filename = file.path(outputFolder, 'study_correlations.pdf')
         , width = 48, height = 60, unit = 'in')

#overlap
combinedCoefs = unique(glmnetCoefs[, gene_sym, lambda])
setnames(combinedCoefs, old = 'lambda', new = 'params')
combinedCoefs[
  , `:=`(params = paste0('lambda_', round(params, 4))
         , model = 'glmnet')]
combinedCoefs = rbind(combinedCoefs
  , unique(zzCoefs[
    , .(gene_sym
        , params = paste0('sumabsv_', sumabsv)
        , model = 'zeitzeiger')]))
combinedCoefs = rbind(combinedCoefs
  , unique(genes2017Coefs[
    !is.na(coef)
    , .(gene_sym 
        , params = '2017'
        , model = 'zeitzeiger')]))
setorderv(combinedCoefs, cols = c('gene_sym', 'model', 'params'))

modelOrder = unique(combinedCoefs[, .(params, model)])
setorderv(modelOrder, 'model')
modelOrder[
  , ord := rev(seq(.N))
  , by = model]
combinedCoefs = merge(combinedCoefs, modelOrder, by = c('model', 'params'))

vennDt = dcast(unique(combinedCoefs[, .(model, params, gene_sym)])
  , gene_sym ~ model + params)
vennDt[, gene_sym := NULL]
vennDt = vennDt[, !duplicated(as.list(vennDt)), with = FALSE]

vennColors = brewer.pal(5, 'Set2')
vennObj = venn.diagram(
  list(vennDt[!is.na(glmnet_lambda_0.1101)]$glmnet_lambda_0.1101
  , vennDt[!is.na(glmnet_lambda_0.1923)]$glmnet_lambda_0.1923
  , vennDt[!is.na(zeitzeiger_2017)]$zeitzeiger_2017
  , vennDt[!is.na(zeitzeiger_sumabsv_2)]$zeitzeiger_sumabsv_2
  , vennDt[!is.na(zeitzeiger_sumabsv_3)]$zeitzeiger_sumabsv_3)
, category.names = c('lambda_0.1101', 'lambda_0.1923', 'zeitzeiger_2017'
, 'sumabsv_2', 'sumabsv_3')
, fill = vennColors, fontfamily = 'sans', cat.fontfamily = 'sans'
, cat.cex = 1, cat.default.pos = 'outer', main.fontfamily = 'sans'
, main = 'Gene overlap between models', filename = NULL)

pdf(file = file.path(outputFolder, 'gene_venn.pdf'), height = 14, width = 14)
grid.newpage()
grid.draw(vennObj)
dev.off()


#study CCDs
ccdDt = foreach(param = unique(combinedCoefs$params), .combine = rbind) %:%
  foreach(cond = c('control', 'perturbation'), .combine = rbind) %dopar% {
  
  genes = unique(combinedCoefs[params == param, gene_sym])
  
  if (cond == 'control') {
  
    esetGSE39445 = studyEsetList$GSE39445
    esetGSE48113 = studyEsetList$GSE48113
    esetGSE56931 = studyEsetList$GSE56931
    
  } else {
    
    esetGSE39445 = studyEsetListPerturb$GSE39445
    esetGSE48113 = studyEsetListPerturb$GSE48113
    esetGSE56931 = studyEsetListPerturb$GSE56931}
  
  GSE39445_GSE48113 = calcCCD(esetGSE39445, esetGSE48113, genes)
  GSE39445_GSE56931 = calcCCD(esetGSE39445, esetGSE56931, genes)
  GSE48113_GSE56931 = calcCCD(esetGSE48113, esetGSE56931, genes)

  distDt = data.table(params = param
    , GSE39445_GSE48113
    , GSE39445_GSE56931
    , GSE48113_GSE56931
    , cond = cond)
  
  return(distDt)}
ccdDtMelt = melt(ccdDt, id.vars = c('params', 'cond'), variable.name = 'studies'
  , value.name = 'ccd')

pCcd = ggplot(data = ccdDtMelt) +
  geom_point(aes(x = params, y = ccd, fill = cond)
    , position = position_dodge(width = 0.3), shape = 21) +
  scale_fill_brewer(type = 'div', direction = -1) +
  coord_flip() +
  ggtitle('Range of CCDs by parameters and conditions, scaled by number of gene pairs')

ggsave(filename = file.path(outputFolder, 'ccd_plot.pdf'), pCcd 
  , height = 12, width = 12, units = 'in')

# overlapGenes = dcast(combinedCoefs, gene_sym ~ model)
# overlapGenes = overlapGenes[, gene_sym := NULL]

# singleGenes = combinedCoefs[
#   , .SD[which.max(count)], by = gene_sym
#   ][count == 1]
# overlapGenes = combinedCoefs[!(gene_sym %in% singleGenes$gene_sym)]

# overlapPlot = ggplot(data = overlapGenes
#   , aes(x = gene_sym, y = model, fill = gene_sym)) +
#   geom_point(shape = 21, size = 2) +
#   # scale_fill_gradient(low = 'white', high = 'black') +
#   ggtitle('Model overlap by gene') +
#   ylab('Model') +
#   xlab('Gene')
# ggexport(overlapPlot, filename = file.path(outputFolder, 'overlap_plot.pdf')
#   , width = 12, height = 12, units = 'in')
