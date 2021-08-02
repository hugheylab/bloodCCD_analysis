source(file.path('code', 'utils.R'))

#loading data
studyMetadata = fread(studyMetadataPath, stringsAsFactors = FALSE)
sampleMetadata = fread(sampleMetadataPath, stringsAsFactors = FALSE)
emat = qread(ematPath)

timeMax = 24
sampleMetadata = convertZt(sampleMetadata)
controlConds = c('Sleep Extension', 'In phase with respect to melatonin', 
                 'baseline')
controlMetadata = sampleMetadata[condition %in% controlConds]

nFold = 10

set.seed(21056)
foldIds = unique(controlMetadata[, .(study, subject)])
foldIds[, foldId := sample(rep_len(1:nFold, .N))]  

sm = merge(controlMetadata, foldIds, by = c('study', 'subject'))
sm = sm[!is.na(clock_time)]

xClock = t(emat)[sm$sample, ]

#### zeitzeiger fit
sumabsv = 2:3
nTime = 12 
nSpc = 1:2

ztFracRange = seq(0, 1, 0.001)

fitResultList = zeitzeigerFitCv(xClock, sm$ztFrac, sm$foldId)

spcResultList = foreach(absv = sumabsv) %dopar% {
  spcResult = zeitzeigerSpcCv(fitResultList, nTime = nTime, sumabsv = absv)}

zzCv = foreach(spcRes = spcResultList, .combine = rbind) %dopar% { 
  predResult = zeitzeigerPredictCv(x = xClock, 
                                   time = sm$ztFrac, 
                                   foldid = sm$foldId, 
                                   spcRes, 
                                   nSpc = nSpc, 
                                   timeRange = ztFracRange)
  
  return(data.table(predResult$timePred))}
set(zzCv, j = 'sample', value = rep.int(rownames(xClock), times = length(sumabsv)))
set(zzCv, j = 'foldId', value = rep.int(sm$foldId, times = length(sumabsv)))
set(zzCv, j = 'sumabsv', value = rep(sumabsv, each = nrow(xClock)))
zzCv = melt(zzCv, 
            id.vars = c('sample', 'foldId' ,'sumabsv'), 
            measure.vars = glue('V{nSpc}'), 
            variable.name = 'nSpc', 
            value.name = 'ztFracPred')
zzCv[, nSpc := as.integer(sub('.', '', nSpc))]
zzCv = merge(zzCv, 
             controlMetadata[, .(study, subject, sample, ztFrac)], 
             by = 'sample')
zzCv[, diffFrac := getCircDiff(ztFrac, ztFracPred) * 24]
zzCv[, absDiffFrac := abs(diffFrac)]
zzCv[, nSpc_fac := factor(nSpc)]
zzCv[, sumabsv_fac := factor(sumabsv)]
zzCv[, study_fac := factor(study)]

#### cross-validation summary
zzCvSumm = zzCv[, .(mse = mean(diffFrac^2), 
                    mae = mean(abs(diffFrac))), 
                by = .(nSpc_fac, sumabsv_fac)]

zzGenes = foreach(ii = 1:length(sumabsv), .combine=rbind) %dopar% {
  genesPerSpc = lapply(spcResultList[[ii]], 
                       function(spcResult) spcResult$v!=0)
  genesCumsumTmp = lapply(genesPerSpc, 
                          function(a) {
                            t(sapply(1:nrow(a), 
                                     function(kk) { cumsum(a[kk,]) }))})
  return(do.call(rbind, lapply(genesCumsumTmp, function(a) colSums(a!=0))))}
zzGenes = data.table(zzGenes[, nSpc])
zzGenes[, sumabsv := rep(sumabsv, each = nFold)]
zzGenes = melt(zzGenes, id.vars = 'sumabsv', variable.name = 'nSpc', 
               value.name = 'nGenes')  
zzGenes[, nSpc := as.integer(sub('.', '', nSpc))]
zzGenes[, nSpc_fac := as.factor(nSpc)]
zzGenes[, sumabsv_fac := as.factor(sumabsv)]

zzGeneSumm = zzGenes[, .(meanGenes = mean(nGenes) , sdGenes = sd(nGenes)), 
                     by = .(sumabsv_fac, nSpc_fac)]

zzSumm = merge(zzGeneSumm, zzCvSumm, by = c('sumabsv_fac', 'nSpc_fac'))
zzSumm[, params := paste0('sumabsv = ', sumabsv_fac, ', nSpc = ', nSpc_fac)]
zzSumm = zzSumm[, .(params, nGenes = meanGenes, mae, model = 'zeitzeiger')]

#### zeitzeiger cv summary plots
pZzMae = ggplot(zzCvSumm) +
  geom_point(aes(x = nSpc_fac, y = mae, fill = sumabsv_fac, shape = sumabsv_fac), 
             size = 4) +
  scale_fill_brewer(name = 'sumabsv', type = 'seq', palette = 'Greys') +
  scale_shape_manual(name = 'sumabsv', values = c(21:23)) +
  labs(x = 'Number of SPCs', y = 'MAE (hours)')

pZzNGenes = ggplot(zzGeneSumm) +
  geom_point(aes(x = nSpc_fac, y = meanGenes, fill = sumabsv_fac, 
                 shape = sumabsv_fac), size = 4) +
  scale_fill_brewer(name = 'sumabsv', type = 'seq', palette = 'Greys') +
  scale_shape_manual(name = 'sumabsv', values = c(21:23)) +
  scale_y_continuous(trans = 'log2') +
  labs(x = 'Number of SPCs', y = 'Number of genes')

pZzCv =  pZzMae / pZzNGenes + 
  plot_annotation(title = 'Zeitzeiger cross-validation')
ggexport(pZzCv, filename = file.path(outputDir, 'zeitzeiger_cv.pdf'), 
         width = 12, height = 6, unit = 'in', dpi = 500)

#### fitting 'best' model from zeitzeiger cv
finalSumAbsv = 2:3
finalNspc = 2
fitResultFinal = zeitzeigerFit(xClock, sm$ztFrac) 
spcCols = glue('spc_{1:finalNspc}')

zzCoefsMelt = foreach(absv = finalSumAbsv, .combine = rbind) %dopar% {

  spcResultFinal = zeitzeigerSpc(fitResultFinal$xFitMean, 
                                 fitResultFinal$xFitResid, 
                                 sumabsv = absv)

  vDt = data.table(spc = 1:length(spcResultFinal$d), 
                   propVar = spcResultFinal$d^2 / sum(spcResultFinal$d^2))

  # getting nonzero coefs
  vCo = data.table(spcResultFinal$v[, 1:finalNspc])
  vCo[, gene := colnames(xClock)]
  setnames(vCo, glue('V{1:finalNspc}'), glue('spc_{1:finalNspc}'))
  vCo = vCo[Reduce('+', mget(spcCols)) != 0]
  setorderv(vCo, order = -1)
  vCo[, gene_sym := as.character(lookUp(as.character(gene), 
                                        'org.Hs.eg', 'SYMBOL', load=TRUE))]
  vCoMelt = melt(vCo, measure.vars = glue('spc_{1:finalNspc}'), 
                 variable.name = 'spc', value.name = 'coef')
  vCoMelt[, gene_fac := factor(gene_sym, levels = rev(vCo$gene_sym))]
  vCoMelt[, nSpc := finalNspc]
  vCoMelt[, sumabsv := absv]}
qsave(zzCoefsMelt, file = file.path(dataDir, 'zeitzeiger_coefs.qs')) 

#### plot of nonzero gene coefs for zeitzeiger
pZzCoefs = plotCoefs(zzCoefsMelt, nrow = 2, sumabsv, spc) +
  ggtitle(glue('Zeitzeiger coefficients for sumabsv = ', 
               '{paste(finalSumAbsv, collapse = \', \')}'))
ggexport(pZzCoefs, filename = file.path(outputDir, 'gene_zeitzeiger_coefs.pdf'))

#compare zeitzeiger 2017
genes2017Dt = qread(file.path('data', 'genes2017.qs'))
genes2017 = unique(genes2017Dt$gene_sym)
 
#### glmnet fits
alph = 0.9

y = cbind(sm[, cos(2*pi*ztFrac)], sm[, sin(2*pi*ztFrac)])
rownames(y) = rownames(xClock)

glmnetCvFit = cv.glmnet(x = xClock, y = y, alpha = alph, foldid = sm[, foldId], 
                        family = 'mgaussian', type.measure = 'mae', keep = TRUE)

glmnetCvDt = as.data.table(glmnetCvFit$fit.preval)
setnames(glmnetCvDt, glue('V{1:3}'), c('sample', 'var', 'lambda'))
glmnetCvDt = dcast(glmnetCvDt, lambda + sample ~ var, value.var = 'value')

lambdaDt = data.table(lambda = glue('s{0:99}'), lambda_val = glmnetCvFit$lambda)
glmnetCvDt = merge(glmnetCvDt, lambdaDt, by = 'lambda')
setnames(glmnetCvDt, c('lambda', 'lambda_val'), c('s', 'lambda'))

glmnetCvDt[, lambda := round(lambda, 7)]
glmnetCvDt = glmnetCvDt[, .(lambda, sample, 
                            ztPred = atan2(y2, y1)/(2*pi))]
glmnetCvDt[ztPred < 0, ztPred := ztPred + 1]

yDt = data.table(y, keep.rownames = 'sample')
yDt = yDt[, .(sample, 
              ztAct = atan2(V2, V1)/(2*pi))]
yDt[ztAct < 0, ztAct := ztAct + 1]

glmnetCvDt = merge(glmnetCvDt, yDt, by = 'sample')
glmnetCvDt[, ae := abs(getCircDiff(ztAct, ztPred))]
glmnetCvDt = glmnetCvDt[, .(mae = mean(ae),
                            ae_sd = sd(ae)), 
                        by = lambda]
glmnetCvDt[, upper := mae + ae_sd]
glmnetCvDt[, lower := mae - ae_sd]
glmnetCvDt = glmnetCvDt[, 24*.SD[, .(mae, ae_sd, upper, lower)], 
                        by = lambda]

#glmnet gene coefficients by lambda
glmnetCoefs = foreach(lambda = glmnetCvFit$lambda, .combine = rbind) %do% {
    
    betas = coef(glmnetCvFit, s = lambda)
    
    cosBetas = data.table(lambda = lambda, 
                          gene = rownames(betas[[1]])[betas[[1]][, 1] != 0], 
                          cos = betas[[1]][betas[[1]][, 1] != 0])
    sinBetas = data.table(lambda = lambda, 
                          gene = rownames(betas[[2]])[betas[[2]][, 1] != 0], 
                          sin = betas[[2]][betas[[2]][, 1] != 0])
  
    betaDt = merge(cosBetas, sinBetas, by = c('lambda', 'gene'), all = TRUE)}
  
glmnetCoefs = glmnetCoefs[gene != '(Intercept)']
glmnetCoefs = melt(glmnetCoefs, id.vars = c('lambda', 'gene'), 
                   measure.vars = c('cos', 'sin'), value.name = 'coef', 
                   variable.name = 'param')
glmnetCoefs[, gene_sym := as.character(lookUp(gene, 'org.Hs.eg', 'SYMBOL', 
                                              load = TRUE))]

lambdaSumm = glmnetCoefs[, .(nGenes = uniqueN(gene_sym)) , by = lambda]
lambdaSumm[, lambda := round(lambda, 7)]

#glmnet plot dt
glmnetPltDt = merge(glmnetCvDt, lambdaSumm, by = 'lambda')

# glmnet summary plots
# glmnet cv plot
pGlmnetMae = ggplot(glmnetPltDt) +
  geom_errorbar(aes(x = lambda, ymax = upper, ymin = lower)) +
  geom_point(aes(x = lambda, y = mae),shape = 21, fill = 'white', color = 'black') +
  scale_x_log10(expression(lambda)) +
  labs(y = 'MAE (h)')

#ngenes by lambda
pGlmnetNGenes = ggplot(glmnetPltDt) +
  geom_point(aes(x = lambda, y = nGenes)) +
  scale_y_continuous(trans = 'log2', breaks = 2^c(0:9)) +
  labs(x = expression(lambda), y = 'Number of genes')

pGlmnetCv = pGlmnetMae / pGlmnetNGenes +
  plot_annotation(title = 'Glmnet cross-validation')
ggexport(pGlmnetCv, filename = file.path(outputDir, 'glmnet_cv.pdf'), 
         width = 12, height = 6, unit = 'in', dpi = 500)

glmnetSumm = glmnetPltDt[, .(nGenes, mae, 
                             params = lambda)]
glmnetSumm[, model := 'glmnet']
pLambdas = glmnetSumm$params[c(65, 77)]

#comparing mae and ngenes between zeitzeiger and glmnet
mdlSumm = rbind(zzSumm, glmnetSumm)

pZzVsGlmnet = ggplot(mdlSumm[nGenes < 60]) +
  geom_point(aes(x = mae, y = nGenes, fill = model), 
             shape = 21, size = 4, alpha = 0.6) +
  geom_text_repel(
    data = mdlSumm[params %in% pLambdas, 
      .(mae, nGenes, model, 
        params = ifelse(params == min(pLambdas), '41 genes', '21 genes'))], 
    box.padding = 0.5, 
    aes(x = mae, y = nGenes, label = params)) +
  labs(x = 'MAE (h)', y = 'Number of genes')
  
# combining glmnet, zeitzeiger cv plots
pCvFinal = (pZzCv | pGlmnetCv | pZzVsGlmnet) +
  plot_annotation(title = 'Cross-validation for Zeitzeiger and elastic net')
ggexport(pCvFinal, filename = file.path(outputDir, 'bloodCCD_cv.pdf'), 
         width = 24, height = 12, units = 'in')

fig1 = pZzMae + pGlmnetMae + pZzVsGlmnet +
  plot_layout(ncol = 2, widths = c(1.25, 2)) +
  plot_annotation(tag_levels = 'A')
ggexport(fig1, filename = file.path(outputDir, 'fig1.pdf'), 
         width = 1080, height = 720)


geneSummGlmnet = glmnetCoefs[round(lambda, 7) %in% pLambdas]
qsave(geneSummGlmnet, file = file.path(dataDir, 'glmnet_coefs.qs'))

pGlmnetCoef = plotCoefs(geneSummGlmnet, ncol = 2, as.factor(lambda), param) +
  ggtitle('Glmnet coefficients by lambda')
ggsave(filename = file.path(outputDir, 'gene_glmnet_coefs.pdf'), 
       plot = pGlmnetCoef, width = 18, height = 18, units = 'in', dpi = 500)

suppFig1 = plotCoefs(geneSummGlmnet[, .SD[lambda == min(lambda)]], ncol = 2, 
                     param)
ggexport(filename = file.path(outputDir, 'suppFig1.pdf'), 
       plot = suppFig1, width = 900, height = 600)

coefsFinal = pZzCoefs | pGlmnetCoef
ggexport(coefsFinal, filename = file.path(outputDir, 'bloodCCD_coefs.pdf'), 
         width = 16, height = 20, units = 'in', dpi = 500)