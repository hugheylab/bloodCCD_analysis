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

theme_set(theme_bw())

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
sumabsv = seq(1, 3, by = 0.5)
nTime = 12 
nSpc = 1:5

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

p1 = ggplot(zzCvSumm) +
  geom_point(aes(x = nSpc_fac, y = sqrt(mse)
    , fill = sumabsv_fac, shape = sumabsv_fac)
    , size = 2.5) +
  scale_fill_brewer(name = 'sumabsv', type = 'seq', palette = 'Greys') +
  scale_shape_manual(name = 'sumabsv', values = c(3, 4, 21:24)) +
  labs(x = 'Number of SPCs', y = 'RMSE') +
  ggtitle('MSE by nSPCs and sumabsv')

p2 = ggplot(zzCvSumm) +
  geom_point(aes(x = nSpc_fac, y = mae
    , fill = sumabsv_fac, shape = sumabsv_fac)
    , size = 2.5) +
  scale_fill_brewer(name = 'sumabsv', type = 'seq', palette = 'Greys') +
  scale_shape_manual(name = 'sumabsv', values = c(3, 4, 21:24)) +
  labs(x = 'Number of SPCs', y = 'MAE') +
  ggtitle('MAE by nSPCs and sumabsv')

p3 = ggplot(zzGeneSumm) +
  geom_point(aes(x = nSpc_fac, y = meanGenes
    , fill = sumabsv_fac, shape = sumabsv_fac)
    , size = 2.5) +
  scale_fill_brewer(name = 'sumabsv', type = 'seq', palette = 'Greys') +
  scale_shape_manual(name = 'sumabsv', values = c(3, 4, 21:24)) +
  scale_y_continuous(trans = 'log2') +
  labs(x = 'Number of SPCs', y = 'log2-Number of genes') +
  ggtitle('Number of genes selected by Zeitzeiger')

cvPlt = plot_grid(p1, p2, p3, nrow = 2, labels = 'AUTO')

#### fitting 'best' model from cv

fitResultFinal = zeitzeigerFit(xClock, sm$ztFrac) 
spcResultFinal = zeitzeigerSpc(fitResultFinal$xFitMean, fitResultFinal$xFitResid
  , sumabsv = 3)

vDt = data.table(spc = 1:length(spcResultFinal$d)
  , propVar = spcResultFinal$d^2 / sum(spcResultFinal$d^2))

# plot to evaluate nspcs
p4 = ggplot(vDt) +
  geom_point(aes(x = spc, y = propVar), size = 2, shape = 1) +
  scale_x_continuous(breaks = seq(1, 10)) +
  labs(x = 'SPC', y = 'Proportion of\nvariance explained')

# getting nonzero coefs
vCoefs = data.table(spcResultFinal$v[, 1:3])
vCoefs[, gene := colnames(xClock)]
setnames(vCoefs, glue('V{1:3}'), glue('spc_{1:3}'))
spcCols = glue('spc_{1:3}')
vCoefs = vCoefs[Reduce('+', mget(spcCols)) > 0]
setorderv(vCoefs, order = -1)
vCoefs[
  , gene_sym := as.character(lookUp(as.character(gene)
      , 'org.Hs.eg', 'SYMBOL', load=TRUE))]
vCoefsMelt = melt(vCoefs, measure.vars = glue('spc_{1:3}')
  , variable.name = 'spc', value.name = 'coeff')
vCoefsMelt[
  , gene_fac := factor(gene_sym, levels = rev(vCoefs$gene_sym))]

#### plot of nonzero gene coefs
p5 = ggplot(vCoefsMelt) + facet_wrap(~ spc, nrow = 1) +
  geom_bar(aes(x = gene_fac, y = coeff), stat = 'identity') +
  # scale_x_discrete(labels = vCoefsMelt$gene_sym) +
  labs(x = 'Gene', y = 'Coefficient') + 
  coord_flip() +
  theme(panel.spacing = unit(1.2, 'lines'))


#### glmnet
alph = 0.9

y = cbind(sm[, cos(2*pi*ztFrac)], sm[, sin(2*pi*ztFrac)])

glmnetCvMse = cv.glmnet(x = xClock, y = y, alpha = alph
  , foldid = sm[, foldId], family = 'mgaussian')

glmnetCvMae = cv.glmnet(x = xClock, y = y, alpha = alph
  , foldid = sm[, foldId], family = 'mgaussian', type.measure = 'mae')


#gene coefficients by lambda
cvCoefs = foreach(lambda = round(glmnetCvMse$lambda, 7)
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

cvPreds = foreach(lambda = round(glmnetCvMse$lambda, 7)
  , .combine = rbind) %dopar% {
    
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
cvPreds[, diffFrac := ztPred - ztFrac]

cvPredSumm = cvPreds[
  , .(mae = mean(abs(diffFrac))
      , mse = mean(diffFrac^2))
  , by = lambda]

lambdaSumm = cvCoefs[
  , .(nGenes = uniqueN(gene_sym))
  , by = lambda]

#plot dts
# pltMse = data.table(lambda = glmnetCvMse$lambda, MSE = glmnetCvMse$cvm
#   , upper = glmnetCvMse$cvup, lower = glmnetCvMse$cvlo)
# 
# pltMae = data.table(lambda = glmnetCvMae$lambda, MAE = glmnetCvMae$cvm
#   , upper = glmnetCvMae$cvup, lower = glmnetCvMae$cvlo)

#summary plots
p6 = ggplot(data = cvPredSumm) +
 # geom_errorbar(aes(x = lambda, ymax = upper, ymin = lower)) +
  geom_point(aes(x = lambda, y = sqrt(mse)),
    shape = 21, fill = 'white', color = 'black') + 
  scale_x_log10('log(lambda)') +
  geom_vline(xintercept = c(glmnetCvMse$lambda.min, glmnetCvMse$lambda.1se)
    , linetype = 'dashed') +
  ggtitle('RMSE by lambda')

p7 = ggplot(data = cvPredSumm) +
  # geom_errorbar(aes(x = lambda, ymax = upper, ymin = lower)) +
  geom_point(aes(x = lambda, y = mae),
    shape = 21, fill = 'white', color = 'black') + 
  scale_x_log10('log(lambda)') +
  geom_vline(xintercept = c(glmnetCvMse$lambda.min, glmnetCvMse$lambda.1se)
    , linetype = 'dashed') +
  ggtitle('Cross-validation by MAE')

p8 = ggplot(data = lambdaSumm) +
  geom_point(aes(x = lambda, y = nGenes)) +
  scale_y_continuous(trans = 'log2') +
  ggtitle('Number of genes by lambda') +
  ylab('log2 - Number of genes')

cvPlt2 = plot_grid(p6, p7, p8, nrow = 2, labels = 'AUTO', align = 'v')

#gene selections
geneSummGlmnet = cvCoefs[
  lambda == round(glmnetCvMse$lambda.1se, 7)
  | lambda == round(glmnetCvMae$lambda.1se, 7)
  | lambda == round(glmnetCvMse$lambda.min, 7)
  | lambda == cvPredSumm[which.min(sqrt(mse)), lambda]
  | lambda == cvPredSumm[which.min(mae), lambda]
  ][
  , lambda_fac := ifelse(lambda == round(glmnetCvMse$lambda.1se, 7), 'mse_1se'
      , ifelse(lambda == round(glmnetCvMae$lambda.1se, 7), 'mae_1se'
      , ifelse(lambda == round(glmnetCvMse$lambda.min, 7), 'min'
      , ifelse(lambda == cvPredSumm[which.min(sqrt(mse)), lambda]
          , 'pred_rmse', 'pred_mae'))))
  , by = gene]

p9 = ggplot(data = geneSummGlmnet[
  lambda_fac == 'mse_1se' 
  & param == 'cos']
  , aes(x = reorder(gene_sym, coef), y = coef)) +
  geom_point(size = 1) +
  geom_segment(aes(x = reorder(gene_sym, coef), xend = reorder(gene_sym, coef)
    , y = 0, yend = coef)) +
  scale_y_continuous(expand = c(0.015, 0)) +
  # facet_grid(~ opt + param) +
  coord_flip() +
  ggtitle('Genes selected by glmnet with lambda.1se') +
  xlab('Gene') +
  ylab('Coefficient') +
  theme(axis.text = element_text(size = 7)
    , axis.text.x = element_text(angle = 90))

p10 = ggplot(data = geneSummGlmnet[
  lambda_fac == 'mae_1se']
  , aes(x = reorder(gene_sym, coef), y = coef)) +
  geom_point(size = 1) +
  geom_segment(aes(x = reorder(gene_sym, coef), xend = reorder(gene_sym, coef)
    , y = 0, yend = coef)) +
  scale_y_continuous(expand = c(0.015, 0)) +
  facet_grid(~ param) +
  coord_flip() +
  ggtitle('Genes selected by glmnet with lambda.1se') +
  xlab('Gene') +
  ylab('Coefficient') +
  theme(axis.text = element_text(size = 7)
    , axis.text.x = element_text(angle = 90))


p10 = ggplot(data = geneSummGlmnet[lambda_fac == 'min']
  , aes(x = reorder(gene_sym, amp), y = amp)) +
  geom_point(size = 1) +
  geom_segment(aes(x = reorder(gene_sym, amp), xend = reorder(gene_sym, amp)
    , y = 0, yend = amp)) +
  scale_y_continuous(expand = c(0.015, 0)) +
  ggtitle('Genes selected by glmnet with lambda.min') +
  xlab('Gene') +
  ylab('Amplitude') +
  theme(axis.text = element_text(size = 5)
    , axis.text.x = element_text(angle = 90))


#prediction
glmnetPredsMin = predict(glmnetCvMse, xClock, s = glmnetCvMse$lambda.min)
glmnetPredsMinDt = as.data.table(glmnetPredsMin)
setnames(glmnetPredsMinDt, c('sample', 'var', 'int', 'value'))
glmnetPredsMinDt = dcast(glmnetPredsMinDt, sample ~ var, value.var = 'value')

glmnetPreds1se = predict(glmnetCvMse, xClock, s = glmnetCvMse$lambda.min)
glmnetPreds1seDt = as.data.table(glmnetPreds1se)
setnames(glmnetPreds1seDt, c('sample', 'var', 'int', 'value'))
glmnetPreds1seDt = dcast(glmnetPreds1seDt, sample ~ var, value.var = 'value')


#comparison with zeitzeiger
geneCombn = dcast(geneSummGlmnet, gene_sym ~ lambda_fac, value.var = 'lambda')
setnames(geneCombn, c('1se', 'min'), c('lam_1se', 'lam_min'))
geneCombn = merge(geneCombn, vCoefs, by = 'gene_sym', all = TRUE)
geneCombn = geneCombn[
  , .(lam_1se = ifelse(!is.na(lam_1se), 1, 0)
      , lam_min = ifelse(!is.na(lam_min), 1, 0)
      , zz = ifelse(!is.na(spc_1) | !is.na(spc_2) | !is.na(spc_3), 1, 0)
      , gene_sym)]


sum(is.nan(glmnetPredsMinDt$y))

