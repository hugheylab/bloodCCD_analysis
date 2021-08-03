source(file.path('code', 'utils.R'))

id = ipcid()

studyMetadata = fread(studyMetadataPath, stringsAsFactors = FALSE)
sampleMetadata = fread(sampleMetadataPath, stringsAsFactors = FALSE)

#loading control data
emat = qread(ematPath)
esetList = qread(esetPath)

#loading perturbation data
ematPerturb = qread(ematPerturbPath)
esetPerturbList = qread(esetPerturbPath)

timeMax = 24
sampleMetadata = convertZt(sampleMetadata)
controlConds = c('Sleep Extension', 'In phase with respect to melatonin', 
                 'baseline')
controlMetadata = sampleMetadata[condition %in% controlConds]

perturbConds = c('Sleep Restriction', 'Out of phase with respect to melatonin', 
                 'sleep deprivation')
perturbMetadata = sampleMetadata[condition %in% perturbConds]

zzCoefs = qread(file.path(dataDir, 'zeitzeiger_coefs.qs'))
setorderv(zzCoefs, c('sumabsv', 'spc', 'coef'))

glmnetCoefs = qread(file.path(dataDir, 'glmnet_coefs.qs'))
setorderv(glmnetCoefs, c('lambda', 'coef'))

genes2017Coefs = qread(file.path(dataDir, 'genes2017.qs'))
setorderv(genes2017Coefs, 'coef')

### merged data
#zeitzeiger
finalSumAbsv = 2:3
zzCormat = foreach(absv = finalSumAbsv, .combine = rbind) %dopar% {
  
  vCo = zzCoefs[sumabsv == absv]
 
  #ipclock to prevent strange interaction with lookUp (Database error)
  ipclock(id)
  zCormat = getCormat(emat, vCo$gene)
  ipcunlock(id)
  
  zCormatSort = sortCormat(zCormat)
  zCormatSort[, sumabsv := absv]}

zzHeatmap =  plotHeatmap(zzCormat, sumabsv) +
  ggtitle(glue('Overall correlation heatmaps, genes selected by zeitzeiger', 
               ', by sumabsv with nSpc = 2')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        text = element_text(size = 16))
  
#glmnet
glmnetCormat = foreach(lam = unique(glmnetCoefs$lambda), 
                       .combine = rbind) %dopar% {
    
    vCo = glmnetCoefs[lambda == lam]
    
    #ipclock to prevent strange interaction with lookUp (Database error)
    ipclock(id)
    gnetCormat = getCormat(emat, vCo$gene)
    ipcunlock(id)
    
    gnetCormatSort = sortCormat(gnetCormat)
    gnetCormatSort[, lambda := lam]}
qsave(glmnetCormat, file = file.path(dataDir, 'glmnet_cor_dt.qs'))

glmnetHeatmap = plotHeatmap(glmnetCormat, lambda) +
  ggtitle(bquote('Overall correlations, genes selected by glmnet, by '~lambda)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        text = element_text(size = 16))

#2017
cormat2017 = getCormat(emat, genes2017Coefs$gene)
cormat2017 = sortCormat(cormat2017)

heatmap2017 =  plotHeatmap(cormat2017) +
  ggtitle('Overall correlation matrix, genes selected by zeitzeiger, 2017') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        text = element_text(size = 16))
heatmap2017 = plot_grid(heatmap2017, NULL, ncol = 2)

combinedCors = copy(glmnetCormat)
setnames(combinedCors, 'lambda', 'params')
combinedCors[, params := paste0('lambda_', round(params, 4))]
combinedCors[, model := 'glmnet']
combinedCors = rbind(combinedCors, 
                     zzCormat[, .(gene1, gene2, rho,
                                  params = paste0('sumabsv_', sumabsv), 
                                  model = 'zeitzeiger')], 
                     cormat2017[, .(gene1, gene2, rho, 
                                    params = '2017', 
                                    model = 'zeitzeiger')])
setorderv(combinedCors, cols = c('gene1', 'gene2', 'model', 'params'))


combinedHeatmap = plotHeatmap(combinedCors, model, params) +
  ggtitle('Comparison of correlations by model/parameter') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        text = element_text(size = 6))
ggexport(combinedHeatmap, 
         filename = file.path(outputDir, 'overall_correlations.pdf'),
         width = 10)

#zeitzeiger
zzStudyCormat = foreach(i = 1:length(esetList), .combine = rbind) %do% {
  
  esetTmp = esetList[[i]]
  esetPerturbTmp = esetPerturbList[[i]]
  
  ematTmp = exprs(esetTmp)  
  ematPerturbTmp = exprs(esetPerturbTmp)
  
  foreach(absv = unique(zzCoefs$sumabsv), .combine = rbind)  %dopar% {
                 
  vCo = zzCoefs[sumabsv == absv]
  
  #ipclock to prevent strange interaction with lookUp (Database error)
  ipclock(id)
  zCormat = getCormat(ematTmp, vCo$gene)
  zCormatPerturb = getCormat(ematPerturbTmp, vCo$gene)
  ipcunlock(id)
  
  zCormatSort = sortCormat(zCormat)
  zCormatSort[, `:=`(study = names(esetList)[i], 
                     sumabsv = absv,
                     condition = 'control')]
  
  zCormatPerturbSort = sortCormat(zCormatPerturb)
  zCormatPerturbSort[, `:=`(study = names(esetPerturbList)[i], 
                            sumabsv = absv,
                            condition = 'perturb')]
  
  zCormatSort = rbind(zCormatSort, zCormatPerturbSort)}}

zzStudyHeatmap = plotHeatmap(zzStudyCormat[condition == 'control'], 
                             sumabsv, study) +
  ggtitle(glue('Correlation heatmaps, genes selected by zeitzeiger', 
               ', by sumabsv with nSpc = 2')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        text = element_text(size = 20))

#glmnet
glmnetStudyCormat = foreach(i = 1:length(esetList), .combine = rbind) %do% { 
  
  esetTmp = esetList[[i]]
  esetPerturbTmp = esetPerturbList[[i]]
  
  ematTmp = exprs(esetTmp)  
  ematPerturbTmp = exprs(esetPerturbTmp)
    
  foreach(lam = unique(glmnetCoefs$lambda), .combine = rbind) %dopar% {
                   
    geneSummGnet = glmnetCoefs[lambda == lam]
                             
    #ipclock to prevent strange interaction with lookUp (Database error)
    ipclock(id)
    gnetCormat = getCormat(esetTmp, unique(geneSummGnet$gene))
    gnetCormatPerturb = getCormat(esetPerturbTmp, unique(geneSummGnet$gene))
    ipcunlock(id)
    
    gnetCormatSort = sortCormat(gnetCormat)
    gnetCormatSort[, `:=`(study = names(esetList)[i], 
                          lambda = lam, 
                          condition = 'control')]
    
    gnetCormatPerturbSort = sortCormat(gnetCormatPerturb)
    gnetCormatPerturbSort[, `:=`(study = names(esetPerturbList)[i], 
                                 lambda = lam, 
                                 condition = 'perturb')]
    
    gnetCormatSort = rbind(gnetCormatSort, gnetCormatPerturbSort)}}

glmnetStudyHeatmap =  plotHeatmap(glmnetStudyCormat[condition == 'control'], 
                                  lambda, study) +
  ggtitle('Correlations by study, genes selected by glmnet') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        text = element_text(size = 10))
ggexport(glmnetStudyHeatmap,
         filename = file.path(outputDir, 'study_cors_glmnet.pdf'), 
         width = 20, height = 10, units = 'in')
 
glmnetCormatFigDt = glmnetCormat
glmnetCormatFigDt[, study := 'Merged']
glmnetCormatFigDt[, condition := 'control']
glmnetCormatFigDt = rbind(glmnetCormatFigDt, glmnetStudyCormat)

fig3 =  plotHeatmap(
  glmnetCormatFigDt[condition == 'control', 
                    .SD[lambda == min(lambda)]], 
  study, ncol = 2, nrow = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
corLevels = levels(glmnetCormat$gene1)
ggexport(fig3, filename = file.path(outputDir, 'fig3.pdf'), 
         width = 12, height = 12, units = 'in')


#### perturb vs control
pGlmnetStudyCondList = foreach(lam = unique(glmnetStudyCormat$lambda)) %dopar% {
  p =  plotHeatmap(glmnetStudyCormat, condition, study) +
    ggtitle(bquote('Correlations by study and condition, genes selected by glmnet, '~lambda[.(lam)])) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
          text = element_text(size = 8))}
pGlmnetStudyCond = ggarrange(plotlist = pGlmnetStudyCondList, nrow = 2)
ggexport(pGlmnetStudyCond, filename = file.path(outputDir, 'glmnet_study_cond_corr.pdf'), 
         width = 15, height = 15, units = 'in')

#2017
study2017Cormat = foreach(i = 1:length(esetList), .combine = rbind) %dopar% {
    
  esetTmp = esetList[[i]]
  esetPerturbTmp = esetPerturbList[[i]]
    
  ematTmp = exprs(esetTmp)  
  ematPerturbTmp = exprs(esetPerturbTmp)
    
  #ipclock to prevent strange interaction with lookUp (Database error)
  ipclock(id)
  cormat = getCormat(esetTmp, unique(genes2017Coefs$gene_sym))
  cormatPerturb = getCormat(esetPerturbTmp, unique(genes2017Coefs$gene_sym))
  ipcunlock(id)
    
  cormatSort = sortCormat(cormat)
  cormatSort[, study := names(esetList)[i]]
  cormatSort[, condition := 'control']
  
  cormatPerturbSort = sortCormat(cormatPerturb)
  cormatPerturbSort[, study := names(esetPerturbList)[i]]
  cormatPerturbSort[, condition := 'perturb']
  
  cormatSort = rbind(cormatSort, cormatPerturbSort)}

heatmapStudy2017 =  plotHeatmap(study2017Cormat[condition == 'control'], study) + 
  ggtitle('Correlations by study, genes selected by zeitzeiger, 2017') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        text = element_text(size = 20))


combinedStudyCors = glmnetStudyCormat[condition == 'control', 
                                      .(gene1, gene2, rho, study, lambda)]
setnames(combinedStudyCors, old = 'lambda', new = 'params')
combinedStudyCors[, params := paste0('lambda_', round(params, 4))]
combinedStudyCors[, model := 'glmnet']
combinedStudyCors = rbind(
  combinedStudyCors, 
  zzStudyCormat[condition == 'control', 
                .(gene1, gene2, rho, study,
                  params = paste0('sumabsv_', sumabsv), model = 'zeitzeiger')])
combinedStudyCors = rbind(
  combinedStudyCors, 
  study2017Cormat[condition == 'control', 
                 .(gene1, gene2, rho, study, 
                   params = '2017', model = 'zeitzeiger')])
setorderv(combinedStudyCors, 
          cols = c('gene1', 'gene2', 'study', 'params', 'model'))


combinedStudyHeatmap = plotHeatmap(combinedStudyCors, model, params, study) +
  ggtitle('Comparison of correlations by model/parameter') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        text = element_text(size = 6))
ggexport(combinedStudyHeatmap, 
         filename = file.path(outputDir, 'study_correlations.pdf'),
         width = 12, height = 12)

#overlap
vennDt = unique(combinedCors[, .(model, params, gene_sym = gene1)])
vennDt[model == 'glmnet', model := 'enet']
vennDt[, params := gsub('\\.', '_', params)]
vennDt = dcast(vennDt, gene_sym ~ model + params)
set(vennDt, j = 'gene_sym', value = NULL)
vennDt = vennDt[, !duplicated(vennDt), with = FALSE]
vennList = lapply(1:5, function (j) 
  vennDt[!which(is.na(vennDt[, j, with = FALSE])), j, with = FALSE][[1]])

vennColors = brewer.pal(5, 'Set2')
suppFig2 = venn.diagram(
  vennList, 
  category.names = names(vennDt), 
  fill = vennColors, fontfamily = 'sans', cat.fontfamily = 'sans', 
  cat.cex = 0.75,  main.fontfamily = 'sans', cat.pos = c(0, 0, 300, 220, 170),
  cat.dist = c(0.175, 0.2, 0.275, 0.2, 0.25), cex = 1, filename = NULL)

pdf(file = file.path(outputDir, 'suppFig2.pdf'))
grid.newpage()
grid.draw(suppFig2)
dev.off()

#study CCDs
genes = unique(combinedCors[
  params == sort(unique((combinedCors[model == 'glmnet']$params))[1]), gene1])

ccdDt = foreach(param = unique(combinedCors$params), .combine = rbind) %:%
  foreach(cond = c('control', 'perturbation'), .combine = rbind) %dopar% {
  
    ids = c('GSE39445', 'GSE48113', 'GSE56931')
    
    if (cond == 'control') {
      
      pair1 = calcCCD(esetList[[ids[1L]]], esetList[[ids[2L]]], genes)
      pair2 = calcCCD(esetList[[ids[1L]]], esetList[[ids[3L]]], genes)
      pair3 = calcCCD(esetList[[ids[2L]]], esetList[[ids[3L]]], genes)

    } else {
    
      pair1 = calcCCD(esetPerturbList[[ids[1L]]], emat, genes)
      pair2 = calcCCD(esetPerturbList[[ids[2L]]], emat, genes)   
      pair3 = calcCCD(esetPerturbList[[ids[3L]]], emat, genes)}
    
    distDt = data.table(pair1, pair2, pair3, params = param, cond = cond)
    distDt = melt(distDt, id.vars = c('params', 'cond'), variable.name = 'pair', 
                  value.name = 'ccd')}

pCcd = ggplot(ccdDt) +
  geom_point(aes(x = params, y = ccd, fill = cond), 
             position = position_dodge(width = 0.3), shape = 21) +
  scale_fill_brewer(type = 'div', direction = -1) +
  coord_flip() +
  ggtitle('Range of CCDs by parameters and conditions, scaled by number of gene pairs')
ggexport(filename = file.path(outputDir, 'ccd_plot.pdf'), pCcd)


glmnetCcdDt = foreach(
  cond = c('control', 'perturbation'), .combine = rbind) %dopar% {
  
  ids = c('GSE39445', 'GSE48113', 'GSE56931')
  
  if (cond == 'control') {
  
    pair1 = calcCCD(esetList[[ids[1L]]], esetList[[ids[2L]]], genes, scale = FALSE)
    pair2 = calcCCD(esetList[[ids[1L]]], esetList[[ids[3L]]], genes, scale = FALSE)
    pair3 = calcCCD(esetList[[ids[2L]]], esetList[[ids[3L]]], genes, scale = FALSE)
    
  } else {
    
    pair1 = calcCCD(esetPerturbList[[ids[1L]]], emat, genes, scale = FALSE)
    pair2 = calcCCD(esetPerturbList[[ids[2L]]], emat, genes, scale = FALSE)   
    pair3 = calcCCD(esetPerturbList[[ids[3L]]], emat, genes, scale = FALSE)}
  
  distDt = data.table(pair1, pair2, pair3, cond = cond)
  distDt = melt(distDt, id.vars = 'cond', variable.name = 'pair', 
                value.name = 'ccd')}

pGlmnetCcd = ggplot(glmnetCcdDt) +
  geom_point(
    aes(x = cond, y = ccd, fill = cond), position = position_dodge(width = 0.3), 
    shape = 21, alpha = 0.5, size = 6) +
  coord_flip() +
  theme(legend.position = 'none') +
  labs(x = 'Condition', y = 'CCD', fill = 'Condition') +
  scale_fill_brewer(palette = 'Dark2', direction = 1)

pGlmnetStudyPerturb = plotHeatmap(
  glmnetCormatFigDt[condition == 'perturb'
                    & lambda == min(lambda)], 
  study, scales = 'fixed', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 5),
        axis.text.y = element_text(size = 5))

fig5 = pGlmnetCcd + pGlmnetStudyPerturb +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(ncol = 1, heights = c(1, 5))
ggexport(filename = file.path(outputDir, 'fig5.pdf'), plot = fig5)

#getting peak phase
controlConds = c('Sleep Extension', 'In phase with respect to melatonin', 
                 'baseline')
sm = sampleMetadata[condition %in% controlConds]

glmnetGenes = unique(glmnetCoefs[lambda == min(lambda), gene])

sm = sm[!is.na(clock_time)]
sm[, time := ztFrac * 24]

fit = getModelFit(emat[, sm$sample], sm)
mFit = getPosteriorFit(fit, covMethod = 'data-driven')
rhyStats = getRhythmStats(mFit, fitType = 'posterior_mean', features = glmnetGenes)

rhyStats[, gene := lookUp(feature, 'org.Hs.eg', 'SYMBOL', load = TRUE)]
rhyStats[, gene_fac := factor(gene, levels = .SD[order(peak_phase), gene])]

pPeakPhase = ggplot(rhyStats) +
  geom_point(aes(x = peak_phase, y = gene_fac), size = 4) +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  theme(
    panel.grid.major.x = eb, panel.grid.minor.x = eb, 
    panel.grid.major.y = element_line(color = 'lightgrey', linetype = 'dotted', 
                                      size = 1)) +
  labs(x = 'Peak phase (h)', y = 'Gene', shape = 'Gene set')

#getting time courses
glmnetTimeCourseDt = getTimeCourseDt(emat, sm, glmnetGenes)
glmnetTimeCourseDt[, gene := factor(gene, levels = rev(levels(rhyStats$gene_fac)))]
pTimeCourse = plotTimeCourse(glmnetTimeCourseDt, ncol = 5, breaks = 6, size = 0.5)
ggexport(filename = file.path(outputDir, 'suppFig3.pdf'), pTimeCourse)

sampGenes = c('HNRNPDL', 'NR1D2', 'PER2')
pGlmnetTimeCourseSamp = plotTimeCourse(glmnetTimeCourseDt[gene %in% sampGenes], 
                                       ncol = 1)

fig2 = pPeakPhase + pGlmnetTimeCourseSamp +
  plot_annotation(tag_levels = 'A') +
  plot_layout(ncol = 2, widths = c(2, 1))
ggexport(filename = file.path(outputDir, 'fig2.pdf'), plot = fig2)

#saving glmnet correlations
ref = getCormat(emat, glmnetGenes, entrezID = TRUE)
diag(ref) = NA
qsave(ref, file.path(dataDir, 'result_blood_ref.qs'))

refGenes = unique(glmnetCoefs[lambda == min(lambda),
                              .(entrez_hs = gene, symbol_hs = gene_sym)])
fwrite(refGenes, file.path(dataDir, 'genes_blood.csv'))
