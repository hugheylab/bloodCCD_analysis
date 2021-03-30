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

theme_set(theme_bw())

outputFolder = file.path('output')

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

nFold = 10

set.seed(21056)
foldIds = unique(controlMetadata[, .(study, subject)])
foldIds[, foldId := sample(rep_len(1:nFold, .N))]  

sm = merge(controlMetadata, foldIds, by = c('study', 'subject'))
noClock = sm[is.na(clock_time), sample]
sm = sm[!is.na(clock_time)]

xClock = t(emat)[sm$sample, ]


#### zeitzeiger
sumabsv = seq(1, 3)
nTime = 12 
nSpc = 1:3

ztFracRange = seq(0, 1, 0.001)

fitResultList = zeitzeigerFitCv(xClock, sm$ztFrac, sm$foldId)

spcResultList = foreach(absv = sumabsv) %dopar% {
  
  spcResult = zeitzeigerSpcCv(fitResultList, nTime = nTime, sumabsv = absv)
  return(spcResult)}

predResultList = foreach(absv = sumabsv) %dopar%  {
  
  predResult = zeitzeigerPredictCv(x = xClock, time = sm$ztFrac
    , foldid = sm$foldId, spcResultList[[absv]], nSpc = nSpc
    , timeRange = ztFracRange)}

timePredList = foreach(pred = predResultList) %dopar% {
  
  return(pred$timePred)}

zzCv = data.table(do.call(rbind, timePredList))
zzCv[
  , `:=`(sample = rep.int(rownames(xClock), times = length(sumabsv))
         , foldId = rep.int(sm$foldId, times = length(sumabsv))
         , sumabsv = rep(sumabsv, each = nrow(xClock)))]
zzCv = melt(zzCv, id.vars = c('sample', 'foldId' ,'sumabsv')
  , measure.vars = glue('V{nSpc}'), variable.name = 'nSpc'
  , value.name = 'ztFracPred')
zzCv[, nSpc := as.integer(sub('.', '', nSpc))]
zzCv = merge(zzCv, controlMetadata[, .(study, subject, sample, ztFrac)]
  , by = 'sample')
zzCv[, diffFrac := getCircDiff(ztFrac, ztFracPred) * 24
  ][, `:=`(absDiffFrac = abs(diffFrac)
           , nSpc_fac = factor(nSpc)
           , sumabsv_fac = factor(sumabsv)
           , study_fac = factor(study))]

#### cross-validation summary

zzCvSumm = zzCv[
  , .(mse = mean(diffFrac^2)
      , mae = mean(abs(diffFrac)))
  , by = .(nSpc_fac, sumabsv_fac)]

# genes by hyperparams
# zzGenes = foreach(absv = 1, .combine = rbind) %dopar% {
#   
#   genesPerSpc = foreach(spcResult = spcResultList[[absv]]) %dopar% {
#     spcResult$v != 0}
#   
#   geneCumSum = foreach(gps = genesPerSpc) %dopar% {
#       t(sapply(1:nrow(gps), function(i) cumsum(gps[i, ])))} 
#   
#   # genesCumsumTmp = foreach(gps = genesPerSpc) %dopar% {
#   #   gpsTmp = 
#   #   
#   # }
#   
#   return(geneCumSum)}

zzGenes = foreach(ii=1:length(sumabsv), .combine=rbind) %dopar% {
  genesPerSpc = lapply(spcResultList[[ii]]
    , function(spcResult) spcResult$v!=0)
  genesCumsumTmp = lapply(genesPerSpc
    , function(a) t(sapply(1:nrow(a), function(kk) cumsum(a[kk,]))))
  return(do.call(rbind, lapply(genesCumsumTmp, function(a) colSums(a!=0))))}
zzGenes = data.table(zzGenes[, nSpc])
zzGenes[, sumabsv := rep(sumabsv, each = length(unique(foldIds$foldId)))]
zzGenes = melt(zzGenes, id.vars = 'sumabsv', variable.name = 'nSpc'
  , value.name = 'nGenes')  
zzGenes[, nSpc := as.integer(sub('.', '', nSpc))]
zzGenes[
  , `:=`(nSpc_fac = as.factor(nSpc)
         , sumabsv_fac = as.factor(sumabsv))]

zzGeneSumm = zzGenes[
  , .(meanGenes = mean(nGenes)
      , sdGenes = sd(nGenes))
  , by = .(sumabsv_fac, nSpc_fac)]

#### summary plots

# p1 = ggplot(zzCvSumm) +
#   geom_point(aes(x = nSpc_fac, y = sqrt(mse)
#     , fill = sumabsv_fac, shape = sumabsv_fac)
#     , size = 2.5) +
#   scale_fill_brewer(name = 'sumabsv', type = 'seq', palette = 'Greys') +
#   scale_shape_manual(name = 'sumabsv', values = c(21:23)) +
#   labs(x = 'Number of SPCs', y = 'RMSE (hours)') +
#   ggtitle('MSE by nSPCs and sumabsv')

p2 = ggplot(zzCvSumm) +
  geom_point(aes(x = nSpc_fac, y = mae
    , fill = sumabsv_fac, shape = sumabsv_fac)
    , size = 2.5) +
  scale_fill_brewer(name = 'sumabsv', type = 'seq', palette = 'Greys') +
  scale_shape_manual(name = 'sumabsv', values = c(21:23)) +
  labs(x = 'Number of SPCs', y = 'MAE (hours)') +
  ggtitle('MAE by nSPCs and sumabsv')

p3 = ggplot(zzGeneSumm) +
  geom_point(aes(x = nSpc_fac, y = meanGenes
    , fill = sumabsv_fac, shape = sumabsv_fac)
    , size = 2.5) +
  scale_fill_brewer(name = 'sumabsv', type = 'seq', palette = 'Greys') +
  scale_shape_manual(name = 'sumabsv', values = 21:23) +
  scale_y_continuous(trans = 'log2') +
  labs(x = 'Number of SPCs', y = 'Number of genes') +
  ggtitle('Number of genes selected by Zeitzeiger')

cvPlt = ggarrange(plotlist = list(p2, p3), nrow = 1, ncol = 2)
cvPlt = annotate_figure(cvPlt, top = text_grob('Zeitzeiger cross-validation'))
ggexport(cvPlt, filename = file.path(outputFolder, 'zeitzeiger_cv.pdf')
  , width = 12, height = 6, unit = 'in', dpi = 500)

#### fitting 'best' model from cv

fitResultFinal = zeitzeigerFit(xClock, sm$ztFrac) 
spcResultFinal = zeitzeigerSpc(fitResultFinal$xFitMean, fitResultFinal$xFitResid
  , sumabsv = 3)
finalNspc = 2

vDt = data.table(spc = 1:length(spcResultFinal$d)
  , propVar = spcResultFinal$d^2 / sum(spcResultFinal$d^2))

# plot to evaluate nspcs
p4 = ggplot(vDt) +
  geom_point(aes(x = spc, y = propVar), size = 2, shape = 1) +
  scale_x_continuous(breaks = seq(1, 10)) +
  labs(x = 'SPC', y = 'Proportion of\nvariance explained')

# getting nonzero coefs
vCoefs = data.table(spcResultFinal$v[, 1:finalNspc])
vCoefs[, gene := colnames(xClock)]
setnames(vCoefs, glue('V{1:finalNspc}'), glue('spc_{1:finalNspc}'))
spcCols = glue('spc_{1:finalNspc}')
vCoefs = vCoefs[Reduce('+', mget(spcCols)) > 0]
setorderv(vCoefs, order = -1)
vCoefs[
  , gene_sym := as.character(lookUp(as.character(gene)
      , 'org.Hs.eg', 'SYMBOL', load=TRUE))]
vCoefsMelt = melt(vCoefs, measure.vars = glue('spc_{1:finalNspc}')
  , variable.name = 'spc', value.name = 'coeff')
vCoefsMelt[
  , gene_fac := factor(gene_sym, levels = rev(vCoefs$gene_sym))]

#### plot of nonzero gene coefs
p5 = ggplot(vCoefsMelt) + facet_wrap(~ spc, nrow = 1) +
  geom_bar(aes(x = gene_fac, y = coeff), stat = 'identity') +
  # scale_x_discrete(labels = vCoefsMelt$gene_sym) +
  labs(x = 'Gene', y = 'Coefficient') + 
  coord_flip() +
  theme(panel.spacing = unit(1.2, 'lines')) +
  ggtitle('Gene coefficients for sumabsv = 3')
ggexport(p5, filename = file.path(outputFolder, 'gene_zeitzeiger_coefs.pdf'))

#### heatmap of gene correlations
zzCormat = cor(t(emat)[, vCoefs$gene], method = 'spearman')
colnames(zzCormat) = as.character(lookUp(as.character(colnames(zzCormat))
  , 'org.Hs.eg', 'SYMBOL', load=TRUE))
rownames(zzCormat) = as.character(lookUp(as.character(rownames(zzCormat))
  , 'org.Hs.eg', 'SYMBOL', load=TRUE))
zzDist = as.dist(1 - zzCormat)/2
zzHc = hclust(zzDist)$order
zzCormatOrd = zzCormat[zzHc, zzHc]

zzCormatDt = as.data.table(zzCormat, keep.rownames = 'gene1')
zzCormatMelt = melt(zzCormatDt, variable.name = 'gene2'
  , value.name = 'rho')

zzHeatmap =  ggplot(zzCormatMelt
  , aes(factor(gene1, levels = unique(colnames(zzCormatOrd)))
      , factor(gene2, levels = unique(colnames(zzCormatOrd)))
      ,  fill = rho))+
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
                       midpoint = 0, limit = c(-1,1), space = 'Lab', 
                       name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  ggtitle('Correlation matrix, genes selected by zeitzeiger') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1)
        , axis.text.y = element_text(size = 8)) +
  coord_fixed()

#zeitzeiger time courses
zzTimeCourseDt = as.data.table(t(emat[vCoefs$gene, sm$sample])
  , keep.rownames = 'sample')
colnames(zzTimeCourseDt)[-1] = as.character(lookUp(colnames(zzTimeCourseDt)[-1]
  , 'org.Hs.eg', 'SYMBOL', load = TRUE))
zzTimeCourseDt = merge(zzTimeCourseDt, sm[, .(study, sample, ztFrac)]
  , by = 'sample')
zzTimeCourseDt = melt(zzTimeCourseDt, id.vars = c('sample', 'study', 'ztFrac')
  , value.name = 'expression', variable.name = 'gene')

pZzTimeCourse = ggplot(zzTimeCourseDt, aes(x = ztFrac*24, y = expression
  , color = study)) +
  geom_point(shape = 1) +
  ylab('Expression') +
  xlab('Hour of day') +
  facet_wrap(~ gene, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  ggtitle('Gene time courses for zeitzeiger')


#compare 2017
genes2017 = readRDS(file.path('data', 'genes2017.rds'))
genes2017 = as.character(genes2017)

emat2017 = emat
rownames(emat2017) = as.character(lookUp(rownames(emat2017)
  , 'org.Hs.eg', 'SYMBOL', load = TRUE))

cormat2017 = cor(t(emat2017)[, unique(genes2017)], method = 'spearman')
dist2017 = as.dist(1 - cormat2017)/2
hc2017 = hclust(dist2017)$order
cormatOrd2017 = cormat2017[hc2017, hc2017]

cormatDt2017 = as.data.table(cormat2017, keep.rownames = 'gene1')
cormatMelt2017 = melt(cormatDt2017, variable.name = 'gene2'
  , value.name = 'rho')

heatmap2017 =  ggplot(cormatMelt2017
  , aes(factor(gene1, levels = unique(colnames(cormatOrd2017)))
      , factor(gene2, levels = unique(colnames(cormatOrd2017)))
      ,  fill = rho))+
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
                       midpoint = 0, limit = c(-1,1), space = 'Lab', 
                       name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  ggtitle('Correlation matrix, genes selected by zeitzeiger, 2017') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1)
        , axis.text.y = element_text(size = 8)) +
  coord_fixed()

timeCourseDt2017 = as.data.table(t(emat2017[unique(genes2017), sm$sample])
  , keep.rownames = 'sample')
timeCourseDt2017 = merge(timeCourseDt2017, sm[, .(study, sample, ztFrac)]
  , by = 'sample')
timeCourseDt2017 = melt(timeCourseDt2017
  , id.vars = c('sample', 'study', 'ztFrac'), value.name = 'expression'
  , variable.name = 'gene')

pTimeCourse2017 = ggplot(timeCourseDt2017, aes(x = ztFrac*24, y = expression
  , color = study)) +
  geom_point(shape = 1) +
  ylab('Expression') +
  xlab('Hour of day') +
  facet_wrap(~ gene, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  ggtitle('Gene time courses for zeitzeiger, 2017')


#### glmnet
alph = 0.9

y = cbind(sm[, cos(2*pi*ztFrac)], sm[, sin(2*pi*ztFrac)])

glmnetCvMse = cv.glmnet(x = xClock, y = y, alpha = alph
  , foldid = sm[, foldId], family = 'mgaussian')

glmnetCvMae = cv.glmnet(x = xClock, y = y, alpha = alph
  , foldid = sm[, foldId], family = 'mgaussian', type.measure = 'mae')


#gene coefficients by lambda
cvCoefs = foreach(lambda = glmnetCvMse$lambda
  , .combine = rbind) %do% {
    
    betas = coef(glmnetCvMse, s = lambda)
    
    cosBetas = data.table(lambda = lambda
      , gene = rownames(betas[[1]])[betas[[1]][, 1] != 0]
      , cos = betas[[1]][betas[[1]][, 1] != 0])
    sinBetas = data.table(lambda = lambda
      , gene = rownames(betas[[2]])[betas[[2]][, 1] != 0]
      , sin = betas[[2]][betas[[2]][, 1] != 0])
  
    betaDt = merge(cosBetas, sinBetas, by = c('lambda', 'gene'), all = TRUE)

    return(betaDt)}
  
cvCoefs = cvCoefs[gene != '(Intercept)']
cvCoefs = melt(cvCoefs, id.vars = c('lambda', 'gene')
  , measure.vars = c('cos', 'sin'), value.name = 'coef'
  , variable.name = 'param')
cvCoefs[
  , gene_sym := as.character(lookUp(gene, 'org.Hs.eg', 'SYMBOL', load = TRUE))]

cvPreds = foreach(lambda = glmnetCvMse$lambda, .combine = rbind) %dopar% {
    
    preds = predict(glmnetCvMse, newx = xClock, s = lambda)
    preds = as.data.table(preds)
    setnames(preds, c(glue('V{1:3}')), c('sample', 'var', 'int'))
    preds[, int := NULL]
    preds[, `lambda` := lambda]
    return(preds)}

cvPreds[, var := ifelse(var == 'y1' ,'cos_ztFrac', 'sin_ztFrac')]
cvPreds = dcast(cvPreds, lambda + sample ~ var, value.var = 'value')
cvPreds[
  , ztPred := atan2(sin_ztFrac, cos_ztFrac)/(2*pi)
  ][, ztPred := ifelse(ztPred < 0, ztPred + 1, ztPred)]
cvPreds = merge(cvPreds, sm[, .(sample, ztFrac)], by = 'sample')
cvPreds[, diffFrac := getCircDiff(ztPred, ztFrac)*24]

cvPredSumm = cvPreds[
  , .(mae = mean(abs(diffFrac))
      , mse = mean(diffFrac^2))
  , by = lambda]

lambdaSumm = cvCoefs[
  , .(nGenes = uniqueN(gene_sym))
  , by = lambda]

#plot dts
pltMse = data.table(lambda = glmnetCvMse$lambda, MSE = glmnetCvMse$cvm
  , upper = glmnetCvMse$cvup, lower = glmnetCvMse$cvlo)

pltMae = data.table(lambda = round(glmnetCvMae$lambda, 7)
  , MAE = glmnetCvMae$cvm, upper = glmnetCvMae$cvup, lower = glmnetCvMae$cvlo)

# summary plots
p6 = ggplot(data = pltMae) +
  geom_errorbar(aes(x = lambda, ymax = upper, ymin = lower)) +
  geom_point(aes(x = lambda, y = MAE),
    shape = 21, fill = 'white', color = 'black') +
  geom_vline(aes(xintercept = glmnetCvMae$lambda.1se), lty = 'dashed'
    , alpha = 0.8) +
  geom_vline(aes(xintercept = glmnetCvMae$lambda.min), lty = 'dashed'
    , alpha = 0.8) +  
  scale_x_log10(expression(lambda)) +
  # scale_y_continuous(breaks = seq(0, 6, by = 0.5)) +
  ggtitle('Cross-validation by MAE') +
  ylab('MAE (radians)') +
  annotate(geom = "text",
    label = c(bquote(lambda[.('1se')]), bquote(lambda[.('min')])),
    x = c(glmnetCvMae$lambda.1se, glmnetCvMae$lambda.min),
    y = c(1, 1),angle = 0, vjust = 1, parse = TRUE)

lambdaTest = merge(lambdaSumm, pltMae, by = 'lambda')
lambdaTest = lambdaTest[
  nGenes < 60
  ][which.min(MAE), lambda]
 

# p7 = ggplot(data = cvPredSumm) +
#   # geom_errorbar(aes(x = lambda, ymax = upper, ymin = lower)) +
#   geom_point(aes(x = lambda, y = mae),
#     shape = 21, fill = 'white', color = 'black') + 
#   scale_x_log10(expression(lambda)) +
#   scale_y_continuous(breaks = seq(0, 6, by = 0.5)) + 
#   ggtitle('Cross-validation by MAE') +
#   ylab('RMSE (hours)')

p8 = ggplot(data = lambdaSumm) +
  geom_point(aes(x = lambda, y = nGenes)) +
  geom_vline(aes(xintercept = glmnetCvMae$lambda.1se), lty = 'dashed'
    , alpha = 0.8) +
  geom_vline(aes(xintercept = glmnetCvMae$lambda.min), lty = 'dashed'
    , alpha = 0.8) + 
  scale_y_continuous(trans = 'log2', breaks = 2^c(0:9)) +
  ggtitle(expression('Number of genes by '~lambda)) +
  ylab('Number of genes') +
  xlab(expression(lambda)) +
  annotate(geom = "text",
    label = c(bquote(lambda[.('1se')]), bquote(lambda[.('min')])),
    x = c(glmnetCvMae$lambda.1se, glmnetCvMae$lambda.min),
    y = c(4, 4),angle = 0, vjust = 1, parse = TRUE)

cvPlt2 = ggarrange(plotlist = list(p6, p8)
  , nrow = 2, align = 'v')
cvPlt2 = annotate_figure(cvPlt2, top = text_grob('Glmnet cross-validation'))
ggexport(cvPlt2, filename = file.path(outputFolder, 'glmnet_cv.pdf')
  , width = 12, height = 12, units = 'in')

#gene selections
top50Coefs = cvCoefs[
  , .SD[abs(coef) %in% head(sort(abs(coef)), 50)]
  , by = .(lambda, param)]

geneSummGlmnet = cvCoefs[lambda == lambdaTest]
  # | lambda == round(glmnetCvMae$lambda.1se, 7)
  # | lambda == round(glmnetCvMse$lambda.min, 7)
  # | lambda == cvPredSumm[which.min(sqrt(mse)), lambda]
  # | lambda == cvPredSumm[which.min(mae), lambda]
  # ][
  # , lambda_fac := ifelse(lambda == round(glmnetCvMse$lambda.1se, 7), 'mse_1se'
      # , ifelse(lambda == round(glmnetCvMae$lambda.1se, 7), 'mae_1se'
      # , ifelse(lambda == round(glmnetCvMse$lambda.min, 7), 'min'
      # , ifelse(lambda == cvPredSumm[which.min(sqrt(mse)), lambda]
          # , 'pred_rmse', 'pred_mae'))))
  # , by = gene]


# pGeneList = foreach(fac = unique(geneSummGlmnet$lambda_fac)) %dopar% {
  
pGlmnetCoef = ggplot(data = geneSummGlmnet
    , aes(x = reorder(gene_sym, coef), y = coef)) +
  
  geom_point(size = 1) +
  geom_segment(aes(x = reorder(gene_sym, coef), xend = reorder(gene_sym, coef)
                   , y = 0, yend = coef)) +
  scale_y_continuous(expand = c(0.015, 0)) +
  facet_grid(~ param) +
  coord_flip() +
  ggtitle(glue('Gene coefficients for ', expression(lambda), ' = {lambdaTest}')) +
  xlab('Gene') +
  ylab('Coefficient') +
  theme(axis.text.y = element_text(size = 6))  
ggsave(filename = file.path(outputFolder, 'gene_glmnet_coefs.pdf')
  , plot = pGlmnetCoef, width = 18, height = 14, units = 'in', dpi = 500)
  
  # return(p)}


geneFig = ggarrange(plotlist = pGeneList, nrow = 1, ncol = 1)
geneFig$`1` = annotate_figure(geneFig$`1`
  , top = text_grob(bquote('50 largest parameter coefficients by '~lambda)))
geneFig$`2` = annotate_figure(geneFig$`2`
  , top = text_grob(bquote('50 largest parameter coefficients by '~lambda)))
geneFig$`3` = annotate_figure(geneFig$`3`
  , top = text_grob(bquote('50 largest parameter coefficients by '~lambda)))
geneFig$`4` = annotate_figure(geneFig$`4`
  , top = text_grob(bquote('50 largest parameter coefficients by '~lambda)))
geneFig$`5` = annotate_figure(geneFig$`5`
  , top = text_grob(bquote('50 largest parameter coefficients by '~lambda)))
ggexport(geneFig
  , filename = file.path(outputFolder, 'gene_glmnet_coefs.pdf'))


# genes for time courses
glmnetTimeCourseDt = as.data.table(t(emat[unique(geneSummGlmnet$gene)
  , sm$sample])
  , keep.rownames = 'sample')
colnames(glmnetTimeCourseDt)[-1] = as.character(
  lookUp(colnames(glmnetTimeCourseDt)[-1], 'org.Hs.eg', 'SYMBOL', load = TRUE))
glmnetTimeCourseDt = merge(glmnetTimeCourseDt, sm[, .(study, sample, ztFrac)]
  , by = 'sample')
glmnetTimeCourseDt = melt(glmnetTimeCourseDt
  , id.vars = c('sample', 'study', 'ztFrac'), value.name = 'expression'
  , variable.name = 'gene')

pGlmnetTimeCourse = ggplot(glmnetTimeCourseDt, aes(x = ztFrac*24, y = expression
  , color = study)) +
  geom_point(shape = 1) +
  ylab('Expression') +
  xlab('Hour of day') +
  facet_wrap(~ gene, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  ggtitle('Gene time courses for glmnet')

ggsave(filename = file.path(outputFolder, 'gene_time_courses.pdf')
  , plot = pTimeCourse, width = 24, height = 24, units = 'in', dpi = 500)
  

glmnetCormat = cor(t(emat)[, geneSummGlmnet$gene], method = 'spearman')
colnames(glmnetCormat) = as.character(lookUp(as.character(colnames(glmnetCormat))
  , 'org.Hs.eg', 'SYMBOL', load=TRUE))
rownames(glmnetCormat) = as.character(lookUp(as.character(rownames(glmnetCormat))
  , 'org.Hs.eg', 'SYMBOL', load=TRUE))
glmnetDist = as.dist(1 - glmnetCormat)/2
glmnetHc = hclust(glmnetDist)$merge
glmnetOpt = order.optimal(glmnetDist, glmnetHc)$order
glmnetCormatOrd = glmnetCormat[glmnetOpt, glmnetOpt]

glmnetCormatDt = as.data.table(glmnetCormat, keep.rownames = 'gene1')
glmnetCormatMelt = melt(glmnetCormatDt, variable.name = 'gene2'
  , value.name = 'rho')

pGlmnetHeatmap =  ggplot(glmnetCormatMelt
  , aes(factor(gene1, levels = unique(colnames(glmnetCormatOrd)))
      , factor(gene2, levels = unique(colnames(glmnetCormatOrd)))
      ,  fill = rho))+
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', 
                       midpoint = 0, limit = c(-1,1), space = 'Lab', 
                       name = 'rho') +
  xlab('Gene') +
  ylab('Gene') +
  ggtitle('Correlation matrix, genes selected by glmnet') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1)
        , axis.text.y = element_text(size = 8)) +
  coord_fixed()

glmnetCoefFig = ggarrange(pGlmnetCoef, pGlmnetHeatmap, pGlmnetTimeCourse
  , nrow = 1, ncol = 1)
ggexport(filename = file.path(outputFolder, 'glmnet_genes.pdf')
  , plot = glmnetCoefFig, width = 24, height = 18, units = 'in')
