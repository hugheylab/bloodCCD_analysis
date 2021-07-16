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
library(deltaccd)
library(patchwork)
library(limorhyde2)
library(BiocParallel)

theme_set(theme_bw(base_size = 25))
id = ipcid()

codeFolder = file.path('code')
outputFolder = file.path('output')
dataFolder = file.path('data')

source(file.path(codeFolder, 'utils.R'))

studyMetadataPath = file.path(dataFolder, 'metadata', 'study_metadata.csv')
studyMetadata = fread(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = file.path(dataFolder, 'metadata', 'sample_metadata.csv')
sampleMetadata = fread(sampleMetadataPath, stringsAsFactors = FALSE)

#loading control data
ematPath = file.path(dataFolder, 'circadian_human_blood_emat.qs')
emat = qread(ematPath)

esetPath = file.path(dataFolder, 'circadian_human_blood.qs')
esetList = qread(esetPath)

#loading perturbation data
esetPerturbPath = file.path(dataFolder, 'subj_norm_pert_esetList.qs')
esetPerturbList = qread(esetPerturbPath)

# ematPerturbPath = file.path('data', 'perturb_emat.qs')
# ematPerturb = qread(ematPerturbPath)

timeMax = 24
sampleMetadata[
  , zt := as.duration(hm(clock_time) - hm(sunrise_time))/as.duration(hours(1))
  ][zt < 0, zt := zt + timeMax]
sampleMetadata[
  , ztFrac := zt/timeMax
  ][, zt := NULL]
controlConds = c('Sleep Extension', 'In phase with respect to melatonin', 
                 'baseline')
controlMetadata = sampleMetadata[condition %in% controlConds]

perturbConds = c('Sleep Restriction', 'Out of phase with respect to melatonin', 
                 'sleep deprivation')
perturbMetadata = sampleMetadata[condition %in% perturbConds]

zzCoefs = qread(file.path(dataFolder, 'zeitzeiger_coefs.qs'))
setorderv(zzCoefs, c('sumabsv', 'spc', 'coef'))

glmnetCoefs = qread(file.path(dataFolder, 'glmnet_coefs.qs'))
setorderv(glmnetCoefs, c('lambda', 'coef'))

genes2017Coefs = qread(file.path(dataFolder, 'genes2017.qs'))
setorderv(genes2017Coefs, c('coef'))

### merged data
#zeitzeiger
finalSumAbsv = 2:3
zzCormat = foreach(absv = finalSumAbsv, .combine = rbind) %dopar% {
  
  vCo = zzCoefs[sumabsv == absv]
 
  ipclock(id)
  zCormat = getCormat(emat, vCo$gene)
  ipcunlock(id)
  
  zCormatSort = sortCormat(zCormat)
  zCormatSort[, sumabsv := absv]
  
  return(zCormatSort)}

zzHeatmap =  plotHeatmap(zzCormat, sumabsv) +
  ggtitle(glue('Overall correlation heatmaps, genes selected by zeitzeiger', 
               ', by sumabsv with nSpc = 2')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        text = element_text(size = 16))
  
#glmnet
glmnetCormat = foreach(lam = unique(glmnetCoefs$lambda), 
                       .combine = rbind) %dopar% {
    
    vCo = glmnetCoefs[lambda == lam]
    
    ipclock(id)
    gnetCormat = getCormat(emat, vCo$gene)
    ipcunlock(id)
    
    gnetCormatSort = sortCormat(gnetCormat)
    gnetCormatSort[, lambda := lam]
    
    return(gnetCormatSort)}
qsave(glmnetCormat, file = file.path(dataFolder, 'glmnet_cor_dt.qs'))

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
  ggtitle('Comparison of correlations by model/parameter')
ggexport(combinedHeatmap, 
         filename = file.path(outputFolder, 'overall_correlations.pdf'), 
         width = 36, height = 36, unit = 'in')

###study data-sets
# studyEsetList = qread(file.path(dataFolder, 'subj_norm_esetList.qs'))
# studyEsetListPerturb = qread(file.path(dataFolder, 'subj_norm__pert_esetList.qs'))

#zeitzeiger
zzStudyCormat = foreach(i = 1:length(esetList), .combine = rbind) %do% {
  
  esetTmp = esetList[[i]]
  esetPerturbTmp = esetPerturbList[[i]]
  
  ematTmp = exprs(esetTmp)  
  ematPerturbTmp = exprs(esetPerturbTmp)
  
  foreach(absv = unique(zzCoefs$sumabsv), .combine = rbind)  %dopar% {
                 
  vCo = zzCoefs[sumabsv == absv]
  
  ipclock(id)
  zCormat = getCormat(ematTmp, vCo$gene)
  zCormatPerturb = getCormat(ematPerturbTmp, vCo$gene)
  ipcunlock(id)
  
  zCormatSort = sortCormat(zCormat)
  zCormatSort[
    , `:=`(study = names(esetList)[i], 
           sumabsv = absv,
           condition = 'control')]
  
  zCormatPerturbSort = sortCormat(zCormatPerturb)
  zCormatPerturbSort[
    , `:=`(study = names(esetPerturbList)[i], 
           sumabsv = absv,
           condition = 'perturb')]
  
  zCormatSort = rbind(zCormatSort, zCormatPerturbSort)
  return(zCormatSort)}}

zzStudyHeatmap = plotHeatmap(zzStudyCormat[condition == 'control'], 
                             sumabsv, study) +
  ggtitle(glue('Correlation heatmaps, genes selected by zeitzeiger'
               , ', by sumabsv with nSpc = 2')) +
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
                             
    ipclock(id)
    gnetCormat = getCormat(esetTmp, unique(geneSummGnet$gene))
    gnetCormatPerturb = getCormat(esetPerturbTmp, unique(geneSummGnet$gene))
    ipcunlock(id)
    
    gnetCormatSort = sortCormat(gnetCormat)
    gnetCormatSort[
      , `:=`(study = names(esetList)[i], 
             lambda = lam, 
             condition = 'control')]
    
    gnetCormatPerturbSort = sortCormat(gnetCormatPerturb)
    gnetCormatPerturbSort[
      , `:=`(study = names(esetPerturbList)[i], 
             lambda = lam, 
             condition = 'perturb')]
    
    gnetCormatSort = rbind(gnetCormat, gnetCormatPerturb)
                             
    return(gnetCormatSort)}}

glmnetStudyHeatmap =  plotHeatmap(glmnetStudyCormat[condition == 'control'], 
                                  lambda, study) +
  ggtitle('Correlations by study, genes selected by glmnet') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
ggexport(glmnetStudyHeatmap,
         filename = file.path(outputFolder, 'study_cors_glmnet.pdf'), 
         width = 20, height = 10, units = 'in')
 

glmnetCormatFigDt = copy(glmnetCormat)
glmnetCormatFigDt[, study := 'Merged']
glmnetCormatFigDt[, condition := 'control']
glmnetCormatFigDt = rbind(glmnetCormatFigDt, glmnetStudyCormat)

fig3 =  plotHeatmap(
  glmnetCormatFigDt[condition == 'control', 
                    .SD[lambda == min(lambda)]], 
  study, ncol = 2, nrow = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
corLevels = levels(glmnetCormat$gene1)
ggexport(fig3, filename = file.path(outputFolder, 'fig3.pdf'), 
         width= 19, height = 19, units = 'in')


#### perturb vs control
pGlmnetStudyCondList = foreach(lam = unique(glmnetStudyCormat$lambda)) %dopar% {
  p =  plotHeatmap(glmnetStudyCormat, condition, study) +
    ggtitle(bquote('Correlations by study and condition, genes selected by glmnet, '~lambda[.(lam)])) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
          text = element_text(size = 8))
  
  return(p)}
pGlmnetStudyCond = ggarrange(plotlist = pGlmnetStudyCondList, nrow = 2)
ggexport(pGlmnetStudyCond, filename = file.path(outputFolder, 'glmnet_study_cond_corr.pdf')
  , width = 18, height = 18, units = 'in')

#2017
study2017Cormat = foreach(i = 1:length(esetList), .combine = rbind) %dopar% {
    
  esetTmp = esetList[[i]]
  esetPerturbTmp = esetPerturbList[[i]]
    
  ematTmp = exprs(esetTmp)  
  ematPerturbTmp = exprs(esetPerturbTmp)
    
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
  
  cormatSort = rbind(cormatSort, cormatPerturbSort)
  
  return(cormatSort)}

heatmapStudy2017 =  plotHeatmap(study2017Cormat[condition == 'control'], study) + 
  ggtitle('Correlations by study, genes selected by zeitzeiger, 2017') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        text = element_text(size = 20))


combinedStudyCors = glmnetStudyCormat[condition == 'control', 
                                      .(gene1, gene2, rho, study, lambda)]
setnames(combinedStudyCors, old = 'lambda', new = 'params')
combinedStudyCors[, params := paste0('lambda_', round(params, 4))]
combinedStudyCors[, model := 'glmnet']
combinedStudyCors = rbind(combinedStudyCors, 
                          zzStudyCormat[condition == 'control'
                            , .(gene1, gene2, rho, study,
                                params = paste0('sumabsv_', sumabsv), 
                                model = 'zeitzeiger')])
combinedStudyCors = rbind(combinedStudyCors, 
                          study2017Cormat[condition == 'control'
                            , .(gene1, gene2, rho, study, 
                                params = '2017', 
                                model = 'zeitzeiger')])
setorderv(combinedStudyCors, 
          cols = c('gene1', 'gene2', 'study', 'params', 'model'))


combinedStudyHeatmap = plotHeatmap(combinedStudyCors, model, params, study) +
  ggtitle('Comparison of correlations by model/parameter')
ggexport(combinedStudyHeatmap, 
         filename = file.path(outputFolder, 'study_correlations.pdf'), 
         width = 48, height = 60, unit = 'in')

#overlap
vennDt = unique(combinedCors[, .(model, params, gene_sym = gene1)])
vennDt[model == 'glmnet', model := 'enet']
vennDt[, params := gsub('\\.', '_', params)]
vennDt = dcast(vennDt, gene_sym ~ model + params)
set(vennDt, j = 'gene_sym', value = NULL)
vennDt = vennDt[, !duplicated(vennDt), with = FALSE]

vennColors = brewer.pal(5, 'Set2')
vennObj = venn.diagram(
  list(vennDt[which(!is.na(vennDt[, 1])), 1], 
       vennDt[which(!is.na(vennDt[, 2])), 2], 
       vennDt[which(!is.na(vennDt[, 3])), 3], 
       vennDt[which(!is.na(vennDt[, 4])), 4], 
       vennDt[which(!is.na(vennDt[, 5])), 5]), 
  category.names = names(vennDt), 
  fill = vennColors, fontfamily = 'sans', cat.fontfamily = 'sans', 
  cat.cex = 1, cat.default.pos = 'outer', main.fontfamily = 'sans', 
  main = 'Gene overlap between models', filename = NULL)

pdf(file = file.path(outputFolder, 'gene_venn.pdf'), height = 1400, width = 1400)
grid.newpage()
grid.draw(vennObj)
dev.off()

suppFig2 = venn.diagram(
  list(vennDt[which(!is.na(vennDt[, 1])), 1], 
       vennDt[which(!is.na(vennDt[, 2])), 2], 
       vennDt[which(!is.na(vennDt[, 3])), 3], 
       vennDt[which(!is.na(vennDt[, 4])), 4], 
       vennDt[which(!is.na(vennDt[, 5])), 5]), 
  category.names = names(vennDt), 
  fill = vennColors, fontfamily = 'sans', cat.fontfamily = 'sans', 
  cat.cex = 6,  main.fontfamily = 'sans', cat.pos = c(0, 0, 300, 220, 160),
  cat.dist = c(0.175, 0.2, 0.275, 0.2, 0.25), cex = 8, filename = NULL)

pdf(file = file.path(outputFolder, 'suppFig2.pdf'), width = 72, height = 72)
grid.newpage()
grid.draw(suppFig2)
dev.off()

#study CCDs
genes = unique(combinedCors[
  params == sort(unique((combinedCors[model == 'glmnet']$params))[1]), gene1])

ccdDt = foreach(param = unique(combinedCors$params), .combine = rbind) %:%
  foreach(cond = c('control', 'perturbation'), .combine = rbind) %dopar% {
  
  if (cond == 'control') {
  
    esetGSE39445 = esetList$GSE39445
    esetGSE48113 = esetList$GSE48113
    esetGSE56931 = esetList$GSE56931
    
    GSE39445_GSE48113 = calcCCD(esetGSE39445, esetGSE48113, genes)
    GSE39445_GSE56931 = calcCCD(esetGSE39445, esetGSE56931, genes)
    GSE48113_GSE56931 = calcCCD(esetGSE48113, esetGSE56931, genes)

    distDt = data.table(GSE39445_GSE48113, GSE39445_GSE56931, GSE48113_GSE56931,
                        params = param, cond = cond)
    
  } else {
    
    esetGSE39445 = esetPerturbList$GSE39445
    esetGSE48113 = esetPerturbList$GSE48113
    esetGSE56931 = esetPerturbList$GSE56931
    
    esetGSE39445_ref = calcCCD(esetGSE39445, emat, genes)
    esetGSE48113_ref = calcCCD(esetGSE48113, emat, genes)   
    esetGSE56931_ref = calcCCD(esetGSE56931, emat, genes)         
    
    distDt = data.table(esetGSE39445_ref, esetGSE48113_ref, esetGSE56931_ref,
                        params = param, cond = cond)}
  
  distDt = melt(distDt, id.vars = c('params', 'cond'), variable.name = 'studies', 
                value.name = 'ccd')
  
  return(distDt)}

pCcd = ggplot(data = ccdDt) +
  geom_point(aes(x = params, y = ccd, fill = cond)
    , position = position_dodge(width = 0.3), shape = 21) +
  scale_fill_brewer(type = 'div', direction = -1) +
  coord_flip() +
  ggtitle('Range of CCDs by parameters and conditions, scaled by number of gene pairs')
ggexport(filename = file.path(outputFolder, 'ccd_plot.pdf'), pCcd 
  , height = 18, width = 18, units = 'in')


glmnetCcdDt = foreach(cond = c('control', 'perturbation'), 
                      .combine = rbind) %dopar% {
  
  if (cond == 'control') {
  
    esetGSE39445 = esetList$GSE39445
    esetGSE48113 = esetList$GSE48113
    esetGSE56931 = esetList$GSE56931
    
    GSE39445_GSE48113 = calcCCD(esetGSE39445, esetGSE48113, genes, scale = FALSE)
    GSE39445_GSE56931 = calcCCD(esetGSE39445, esetGSE56931, genes, scale = FALSE)
    GSE48113_GSE56931 = calcCCD(esetGSE48113, esetGSE56931, genes, scale = FALSE)

    distDt = data.table(
      GSE39445_GSE48113
      , GSE39445_GSE56931
      , GSE48113_GSE56931
      , cond = cond)
    
  } else {
    
    esetGSE39445 = esetPerturbList$GSE39445
    esetGSE48113 = esetPerturbList$GSE48113
    esetGSE56931 = esetPerturbList$GSE56931
    
    esetGSE39445_ref = calcCCD(esetGSE39445, emat, genes, scale = FALSE)
    esetGSE48113_ref = calcCCD(esetGSE48113, emat, genes, scale = FALSE)   
    esetGSE56931_ref = calcCCD(esetGSE56931, emat, genes, scale = FALSE)         
    
    distDt = data.table(
      esetGSE39445_ref
      , esetGSE48113_ref
      , esetGSE56931_ref
      , cond = cond)}
  
  distDt = melt(distDt, id.vars = 'cond', variable.name = 'studies'
    , value.name = 'ccd')
  
  return(distDt)}

pGlmnetCcd = ggplot(data = glmnetCcdDt) +
  geom_point(aes(x = cond, y = ccd, fill = cond)
    , position = position_dodge(width = 0.3), shape = 21, alpha = 0.5
    , size = 6) +
  coord_flip() +
  theme(legend.position = 'none') +
  labs(fill = 'Condition') +
  scale_fill_brewer(palette = 'Dark2', direction = 1) +
  xlab('Condition') +
  ylab('CCD')

pGlmnetStudyPerturb = plotHeatmap(glmnetCormatFigDt[condition == 'perturb'
                                                    & lambda == min(lambda)], 
                                  study, scales = 'fixed', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14))

fig5 = pGlmnetCcd + pGlmnetStudyPerturb +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(ncol = 1, heights = c(1, 5))
ggexport(filename = file.path(outputFolder, 'fig5.pdf'), 
         plot = fig5, width = 22, height = 22, units = 'in')

#getting peak phase
controlConds = c('Sleep Extension', 'In phase with respect to melatonin', 
                 'baseline')
sm = sampleMetadata[condition %in% controlConds]

glmnetGenes = unique(glmnetCoefs[lambda == min(lambda), gene])

sm = sm[!is.na(clock_time)]
sm[, time := ztFrac * 24]

fit = getModelFit(emat[sm$sample], sm)
mFit = getPosteriorFit(fit, covMethod = 'data-driven')
rhyStats = getRhythmStats(mFit, fitType = 'posterior_mean', features = glmnetGenes)

rhyStats[, gene := lookUp(feature, 'org.Hs.eg', 'SYMBOL', load = TRUE)]
rhyStats[, gene_fac := factor(gene, levels = .SD[order(peak_phase), gene])]

pPeakPhase = ggplot(rhyStats) +
  geom_point(aes(x = peak_phase, y = gene_fac)
             , size = 4) +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  theme(panel.grid.major.x = element_blank()
        , panel.grid.minor.x = element_blank()
        , panel.grid.major.y = element_line(color = 'lightgrey', linetype = 'dotted'
                                            , size = 1)) +
  labs(shape = 'Gene set') +
  xlab('Peak phase (h)') + 
  ylab('Gene')

#getting time courses
glmnetTimeCourseDt = getTimeCourseDt(emat, sm, glmnetGenes)
glmnetTimeCourseDt[, gene := factor(gene, levels = rev(levels(rhyStats$gene_fac)))]
pTimeCourse = plotTimeCourse(glmnetTimeCourseDt, ncol = 5, breaks = 6)
ggexport(filename = file.path(outputFolder, 'suppFig3.pdf'), pTimeCourse
         , height = 24, width = 16, units = 'in')

sampGenes = c('HNRNPDL', 'NR1D2', 'PER2')
pGlmnetTimeCourseSamp = plotTimeCourse(glmnetTimeCourseDt[gene %in% sampGenes], 
                                       ncol = 1)

fig2 = pPeakPhase + pGlmnetTimeCourseSamp +
  plot_annotation(tag_levels = 'A') +
  plot_layout(ncol = 2, widths = c(2, 1))
ggexport(filename = file.path(outputFolder, 'fig2.pdf')
         , plot = fig2, width = 28, height = 16, units = 'in')

#saving glmnet correlations
ref = getCormat(emat, glmnetGenes, entrezID = TRUE)
diag(ref) = NA
qsave(ref, file.path(dataFolder, 'result_blood_ref.qs'))

refGenes = unique(glmnetCoefs[lambda == min(lambda),
                              .(entrez_hs = gene, symbol_hs = gene_sym)])
fwrite(refGenes, file.path(dataFolder, 'genes_blood.csv'))
