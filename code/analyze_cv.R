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

codeFolder = file.path('code')
outputFolder = file.path('output')
dataFolder = file.path('data')

source(file.path(codeFolder, 'cv_utils.R'))

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
sumabsv = 2:3
nTime = 12 
nSpc = 1:2

ztFracRange = seq(0, 1, 0.001)

fitResultList = zeitzeigerFitCv(xClock, sm$ztFrac, sm$foldId)

spcResultList = foreach(absv = sumabsv) %dopar% {
  spcResult = zeitzeigerSpcCv(fitResultList, nTime = nTime, sumabsv = absv)
  return(spcResult)
  }

predResultList = foreach(spcRes = spcResultList) %dopar% { 
  predResult = zeitzeigerPredictCv(x = xClock
                                   , time = sm$ztFrac
                                   , foldid = sm$foldId
                                   , spcRes
                                   , nSpc = nSpc
                                   , timeRange = ztFracRange
                                   )
  }

timePredList = foreach(pred = predResultList) %dopar% {
  return(pred$timePred)
  }

zzCv = data.table(do.call(rbind, timePredList))
zzCv[
  , `:=`(sample = rep.int(rownames(xClock), times = length(sumabsv))
         , foldId = rep.int(sm$foldId, times = length(sumabsv))
         , sumabsv = rep(sumabsv, each = nrow(xClock))
         )
  ]
zzCv = melt(zzCv
            , id.vars = c('sample', 'foldId' ,'sumabsv')
            , measure.vars = glue('V{nSpc}')
            , variable.name = 'nSpc'
            , value.name = 'ztFracPred'
            )
zzCv[, nSpc := as.integer(sub('.', '', nSpc))]
zzCv = merge(zzCv
             , controlMetadata[, .(study, subject, sample, ztFrac)]
             , by = 'sample'
             )
zzCv[, diffFrac := getCircDiff(ztFrac, ztFracPred) * 24
  ][, `:=`(absDiffFrac = abs(diffFrac)
           , nSpc_fac = factor(nSpc)
           , sumabsv_fac = factor(sumabsv)
           , study_fac = factor(study)
           )
    ]

#### cross-validation summary

zzCvSumm = zzCv[
  , .(mse = mean(diffFrac^2)
      , mae = mean(abs(diffFrac))
      )
  , by = .(nSpc_fac, sumabsv_fac)
  ]

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


zzSumm = merge(zzGeneSumm, zzCvSumm, by = c('sumabsv_fac', 'nSpc_fac'))
zzSumm[, params := paste0('sumabsv = ', sumabsv_fac, ', nSpc = ', nSpc_fac)]
zzSumm = zzSumm[
  , .(params
      , nGenes = meanGenes
      , mae
      , model = 'zeitzeiger')]

#### summary plots
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
  scale_shape_manual(name = 'sumabsv', values = c(21:23)) +
  scale_y_continuous(trans = 'log2') +
  labs(x = 'Number of SPCs', y = 'Number of genes') +
  ggtitle('Number of genes selected by Zeitzeiger')

cvPlt = ggarrange(plotlist = list(p2, p3), nrow = 2)
cvPlt = annotate_figure(cvPlt, top = text_grob('Zeitzeiger cross-validation'))
ggexport(cvPlt, filename = file.path(outputFolder, 'zeitzeiger_cv.pdf')
  , width = 12, height = 6, unit = 'in', dpi = 500)

#### fitting 'best' model from cv
finalSumAbsv = 2:3
finalNspc = 2
fitResultFinal = zeitzeigerFit(xClock, sm$ztFrac) 
spcCols = glue('spc_{1:finalNspc}')

vCoefsMelt = foreach(absv = finalSumAbsv, .combine = rbind) %dopar% {

  spcResultFinal = zeitzeigerSpc(fitResultFinal$xFitMean, fitResultFinal$xFitResid
    , sumabsv = absv)

  vDt = data.table(spc = 1:length(spcResultFinal$d)
    , propVar = spcResultFinal$d^2 / sum(spcResultFinal$d^2))

  # getting nonzero coefs
  vCo = data.table(spcResultFinal$v[, 1:finalNspc])
  vCo[, gene := colnames(xClock)]
  setnames(vCo, glue('V{1:finalNspc}'), glue('spc_{1:finalNspc}'))
  vCo = vCo[Reduce('+', mget(spcCols)) != 0]
  setorderv(vCo, order = -1)
  vCo[
    , gene_sym := as.character(lookUp(as.character(gene)
      , 'org.Hs.eg', 'SYMBOL', load=TRUE))]
  vCoMelt = melt(vCo, measure.vars = glue('spc_{1:finalNspc}')
    , variable.name = 'spc', value.name = 'coeff')
  vCoMelt[
    , gene_fac := factor(gene_sym, levels = rev(vCo$gene_sym))]
  vCoMelt[
    , `:=`(nSpc = finalNspc
           , sumabsv = absv)]}
qsave(vCoefsMelt, file = file.path(dataFolder, 'zeitzeiger_coefs.qs')) 

#### plot of nonzero gene coefs
p5 = ggplot(data = vCoefsMelt 
    , aes(x = gene_fac, y = coeff)) +
  
  geom_point(size = 1) +
  geom_segment(aes(x = gene_fac, xend =gene_fac 
                   , y = 0, yend = coeff)) +
  scale_y_continuous(expand = c(0.015, 0)) +
  facet_wrap(~ sumabsv + spc, scales = 'free_y'
    , nrow = 2) +
  coord_flip() +
  xlab('Gene') +
  ylab('Coefficient') +
  theme(axis.text.y = element_text(size = 8)) +
  ggtitle(glue('Zeitzeiger coefficients for sumabsv = '
    , '{paste(finalSumAbsv, collapse = \', \')}'))
ggexport(p5, filename = file.path(outputFolder, 'gene_zeitzeiger_coefs.pdf'))

#zeitzeiger time courses
zzTimeCourseDt = getTimeCourseDt(emat, sm, unique(vCoefsMelt$gene))

pZzTimeCourse = ggplot(zzTimeCourseDt, aes(x = ztFrac*24, y = expression
  , color = study)) +
  geom_point(shape = 1) +
  ylab('Expression') +
  xlab('Hour of day') +
  facet_wrap(~ gene, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  ggtitle('Gene time courses for zeitzeiger')

ggsave(filename = file.path(outputFolder, 'gene_time_courses_zeitzeiger.pdf')
       , plot = pZzTimeCourse, width = 24, height = 24, units = 'in', dpi = 500)

#compare 2017
genes2017Dt = qread(file.path('data', 'genes2017.qs'))
genes2017 = unique(genes2017Dt$gene_sym)

timeCourseDt2017 = getTimeCourseDt(emat, sm, genes2017)

pTimeCourse2017 = ggplot(timeCourseDt2017, aes(x = ztFrac*24, y = expression
  , color = study)) +
  geom_point(shape = 1) +
  ylab('Expression') +
  xlab('Hour of day') +
  facet_wrap(~ gene, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  ggtitle('Gene time courses for zeitzeiger, 2017')

ggsave(filename = file.path(outputFolder, 'gene_time_courses_2017.pdf')
  , plot = pTimeCourse2017, width = 24, height = 24, units = 'in', dpi = 500)
 

#### glmnet
alph = 0.9

y = cbind(sm[, cos(2*pi*ztFrac)], sm[, sin(2*pi*ztFrac)])
rownames(y) = rownames(xClock)

glmnetCvMae = cv.glmnet(x = xClock, y = y, alpha = alph
  , foldid = sm[, foldId], family = 'mgaussian', type.measure = 'mae'
  , keep = TRUE)

glmnetCvDt = as.data.table(glmnetCvMae$fit.preval)
setnames(glmnetCvDt, glue('V{1:3}'), c('sample', 'var', 'lambda'))
glmnetCvDt = dcast(glmnetCvDt, lambda + sample ~ var, value.var = 'value')

lambdaLookupDt = data.table(s = glue('s{0:99}'), lambda = glmnetCvMae$lambda)
lambdaLookup = match(glmnetCvDt$lambda, lambdaLookupDt$s)

glmnetCvDt[, lambda := round(lambdaLookupDt$lambda[lambdaLookup], 7)]
glmnetCvDt = glmnetCvDt[
  , .(lambda
      , sample
      , ztPred = atan2(y2, y1)/(2*pi))]
glmnetCvDt[ztPred < 0
  , ztPred := ztPred + 1]

yDt = data.table(sample = rownames(y), y1 = y[, 1], y2 = y[, 2])
yDt = yDt[
  , .(sample
      , ztAct = atan2(y2, y1)/(2*pi))]
yDt[ztAct < 0
    , ztAct := ztAct + 1]

glmnetCvDt = merge(glmnetCvDt, yDt, by = 'sample')
glmnetCvDt[, ae := abs(getCircDiff(ztAct, ztPred))]
glmnetCvDt = glmnetCvDt[
  , .(mae = mean(ae)
      , ae_sd = sd(ae))
  , by = lambda]
glmnetCvDt[
  , `:=`(upper = mae + ae_sd
         , lower = mae - ae_sd)]
glmnetCvDt = glmnetCvDt[
  , 24*.SD[, .(mae, ae_sd, upper, lower)]
  , by = lambda]


#gene coefficients by lambda
cvCoefs = foreach(lambda = glmnetCvMae$lambda
  , .combine = rbind) %do% {
    
    betas = coef(glmnetCvMae, s = lambda)
    
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


lambdaSumm = cvCoefs[
  , .(nGenes = uniqueN(gene_sym)) 
  , by = lambda]
lambdaSumm[, lambda := round(lambda, 7)]

#plot dt
pltGlmnet = merge(glmnetCvDt, lambdaSumm, by = 'lambda')

# summary plots
# cross-validation by mae
p6 = ggplot(data = pltGlmnet) +
  geom_errorbar(aes(x = lambda, ymax = upper, ymin = lower)) +
  geom_point(aes(x = lambda, y = mae),
    shape = 21, fill = 'white', color = 'black') +
  geom_vline(aes(xintercept = glmnetCvMae$lambda.1se), lty = 'dashed'
    , alpha = 0.8) +
  geom_vline(aes(xintercept = glmnetCvMae$lambda.min), lty = 'dashed'
    , alpha = 0.8) +  
  scale_x_log10(expression(lambda)) +
  ggtitle('Cross-validation by MAE') +
  ylab('MAE (hours)') +
  annotate(geom = "text",
    label = c(bquote(lambda[.('1se')]), bquote(lambda[.('min')])),
    x = c(glmnetCvMae$lambda.1se, glmnetCvMae$lambda.min),
    y = c(5, 5),angle = 0, vjust = 1, parse = TRUE)

p7 = ggplot(data = pltGlmnet) +
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

cvPlt2 = ggarrange(plotlist = list(p6, p7)
  , nrow = 2, align = 'v')
cvPlt2 = annotate_figure(cvPlt2, top = text_grob('Glmnet cross-validation'))


glmnetSumm = pltGlmnet[
 , .(params = lambda
     , nGenes
     , mae
     , model = 'glmnet')]
pLambdas = glmnetSumm$params[c(65, 77)]

mdlSumm = rbind(zzSumm, glmnetSumm)

p8 = ggplot(data = mdlSumm[nGenes < 60]) +
  geom_point(aes(x = mae, y = nGenes, fill = model)
    , shape = 21, size = 2, alpha = 0.6) +
  geom_text_repel(data = mdlSumm[params %in% pLambdas
                                 | model == 'zeitzeiger']
    , aes(x = mae, y = nGenes, label = params), size = 3) +
  ylab('Number of genes') +
  xlab('MAE (hours)') +
  ggtitle('MAE vs number of genes')
  

cvPltFinal = ggarrange(plotlist = list(cvPlt, cvPlt2, p8)
  , ncol = 3, align = 'v')
ggexport(cvPltFinal, filename = file.path(outputFolder, 'bloodCCD_cv.pdf')
  , width = 24, height = 12, units = 'in')


geneSummGlmnet = cvCoefs[round(lambda, 7) %in% pLambdas]
qsave(geneSummGlmnet, file = file.path(dataFolder, 'glmnet_coefs.qs'))


pGlmnetCoef = ggplot(data = geneSummGlmnet
    , aes(x = reorder(gene_sym, coef), y = coef)) +
  
  geom_point(size = 1) +
  geom_segment(aes(x = reorder(gene_sym, coef), xend = reorder(gene_sym, coef)
                   , y = 0, yend = coef)) +
  scale_y_continuous(expand = c(0.015, 0)) +
  facet_wrap(~ as.factor(lambda) + param, scales = 'free_y'
    , ncol = 2) +
  coord_flip() +
  xlab('Gene') +
  ylab('Coefficient') +
  theme(axis.text.y = element_text(size = 8)) +
  ggtitle('Glmnet coefficients by lambda')
ggsave(filename = file.path(outputFolder, 'gene_glmnet_coefs.pdf')
  , plot = pGlmnetCoef, width = 18, height = 18, units = 'in', dpi = 500)

coefFig = ggarrange(p5, pGlmnetCoef, nrow = 2)
ggexport(coefFig, filename = file.path(outputFolder, 'bloodCCD_coefs.pdf')
  , width = 16, height = 20, units = 'in', dpi = 500)

# genes for time courses
getTimeCourseDt(emat, sm, unique(geneSummGlmnet$gene))

pGlmnetTimeCourse = ggplot(glmnetTimeCourseDt, aes(x = ztFrac*24, y = expression
  , color = study)) +
  geom_point(shape = 1) +
  ylab('Expression') +
  xlab('Hour of day') +
  facet_wrap(~ gene, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 24, by = 2), limits = c(0, 24)) +
  ggtitle('Gene time courses for glmnet')
ggsave(filename = file.path(outputFolder, 'gene_time_courses_glmnet.pdf')
  , plot = pGlmnetTimeCourse, width = 24, height = 24, units = 'in', dpi = 500)
