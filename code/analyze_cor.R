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

theme_set(theme_bw(base_size = 25))

codeFolder = file.path('code')
outputFolder = file.path('output')
dataFolder = file.path('data')

source(file.path(codeFolder, 'cor_utils.R'))

studyMetadataPath = file.path(dataFolder, 'metadata', 'study_metadata.csv')
studyMetadata = fread(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = file.path(dataFolder, 'metadata', 'sample_metadata.csv')
sampleMetadata = fread(sampleMetadataPath, stringsAsFactors = FALSE)

#loading control data
ematPath = file.path(dataFolder, 'circadian_human_blood_emat.qs')
emat = qread(ematPath)

esetPath = file.path(dataFolder, 'circadian_human_blood.qs')
eset = qread(esetPath)

#loading perturbation data
esetPerturbPath = file.path(dataFolder, 'subj_norm_pert_esetList.qs')
esetPerturb = qread(esetPerturbPath)

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
 
  zCormat = getCormat(emat, vCo$gene)
  
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
    
    gnetCormat = getCormat(emat, vCo$gene)
    
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
setnames(combinedCors, old = 'lambda', new = 'params')
combinedCors[, params := paste0('lambda_', round(params, 4))]
combinedCors[, model := 'glmnet']
combinedCors = rbind(combinedCors, 
                     zzCormat[, .(gene1, gene2, rho,
                                  params = paste0('sumabsv_', sumabsv), 
                                  model = 'zeitzeiger')])
combinedCors = rbind(combinedCors, 
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
zzStudyCormat = foreach(i = 1:length(eset)
  , .combine = rbind) %:% 
  foreach(cond = c('control', 'perturb'), .combine = rbind) %do% {
    
    if(cond == 'control') {
      esetTmp = eset[[i]]
      } else {esetTmp = esetPerturb[[i]]}
    
    ematTmp = exprs(esetTmp)  
    
    foreach(absv = unique(zzCoefs$sumabsv), .combine = rbind)  %dopar% {
                   
    vCo = zzCoefs[sumabsv == absv]
    
    zCormat = getCormat(ematTmp, vCo$gene)
    
    zCormatSort = sortCormat(zCormat)
    zCormatSort[
      , `:=`(study = names(eset)[i], 
             sumabsv = absv,
             condition = cond)]
  
    return(zCormatSort)}}

zzStudyHeatmap = plotHeatmap(zzStudyCormat[condition == 'control'], 
                             sumabsv, study) +
  ggtitle(glue('Correlation heatmaps, genes selected by zeitzeiger'
               , ', by sumabsv with nSpc = 2')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        text = element_text(size = 20))

#glmnet
glmnetStudyCormat = foreach(i = 1:length(eset)
  , .combine = rbind) %:% 
  foreach(cond = c('control', 'perturb'), .combine = rbind) %do% {
              
  if(cond == 'control') {
   
     esetTmp = eset[[i]]
    
  } else {esetTmp = esetPerturb[[i]]}
    
  foreach(lam = unique(glmnetCoefs$lambda), .combine = rbind) %dopar% {
                   
    geneSummGnet = glmnetCoefs[lambda == lam]
                             
    gnetCormat = getCormat(esetTmp, unique(geneSummGnet$gene))
    
    gnetCormatSort = sortCormat(gnetCormat)
    gnetCormatSort[
      , `:=`(study = names(eset)[i], 
             lambda = lam, 
             condition = cond)]
                             
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
study2017Cormat = foreach(i = 1:length(eset), .combine = rbind) %:% 
  foreach(cond = c('control', 'perturb'), .combine = rbind) %dopar% {
    
    if(cond == 'control') {
      esetTmp = eset[[i]]
      } else {esetTmp = esetPerturb[[i]]}
    
    esetTmp = eset[[i]]  
    
    cormat = getCormat(esetTmp, unique(genes2017Coefs$gene_sym))
    
    cormatSort = sortCormat(cormat)
    cormatSort[, study := names(eset)[i]]
    cormatSort[, condition := cond]
  
    return(cormatSort)}

heatmapStudy2017 =  plotHeatmap(study2017Cormat[condition == 'control'], study) + 
  ggtitle('Correlations by study, genes selected by zeitzeiger, 2017') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), 
        text = element_text(size = 20))


combinedStudyCors = copy(
  glmnetStudyCormat[condition == 'control'
    , .(gene1, gene2, rho, study, lambda)])
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
vennDt = dcast(unique(combinedCors[, .(model, params, gene_sym = gene1)]), 
               gene_sym ~ model + params)
vennDt[, gene_sym := NULL]
vennDt = vennDt[, !duplicated(as.list(vennDt)), with = FALSE]

vennColors = brewer.pal(5, 'Set2')
vennObj = venn.diagram(
  list(vennDt[!is.na(glmnet_lambda_0.1101)]$glmnet_lambda_0.1101, 
       vennDt[!is.na(glmnet_lambda_0.1923)]$glmnet_lambda_0.1923, 
       vennDt[!is.na(zeitzeiger_2017)]$zeitzeiger_2017, 
       vennDt[!is.na(zeitzeiger_sumabsv_2)]$zeitzeiger_sumabsv_2, 
       vennDt[!is.na(zeitzeiger_sumabsv_3)]$zeitzeiger_sumabsv_3), 
  category.names = c('enet_lambda_0.1101', 'enet_lambda_0.1923', 'zeitzeiger_2017', 
                     'zeitzeiger_sumabsv_2', 'zeitzeiger_sumabsv_3'), 
  fill = vennColors, fontfamily = 'sans', cat.fontfamily = 'sans', 
  cat.cex = 1, cat.default.pos = 'outer', main.fontfamily = 'sans', 
  main = 'Gene overlap between models', filename = NULL)

pdf(file = file.path(outputFolder, 'gene_venn.pdf'), height = 1400, width = 1400)
grid.newpage()
grid.draw(vennObj)
dev.off()

suppFig2 = venn.diagram(
  list(vennDt[!is.na(glmnet_lambda_0.1101)]$glmnet_lambda_0.1101, 
       vennDt[!is.na(glmnet_lambda_0.1923)]$glmnet_lambda_0.1923, 
       vennDt[!is.na(zeitzeiger_2017)]$zeitzeiger_2017, 
       vennDt[!is.na(zeitzeiger_sumabsv_2)]$zeitzeiger_sumabsv_2, 
       vennDt[!is.na(zeitzeiger_sumabsv_3)]$zeitzeiger_sumabsv_3), 
  category.names = c('enet_lambda_0.1101', 'enet_lambda_0.1923', 'zeitzeiger_2017', 
                     'zeitzeiger_sumabsv_2', 'zeitzeiger_sumabsv_3'), 
  fill = vennColors, fontfamily = 'sans', cat.fontfamily = 'sans', 
  cat.cex = 6,  main.fontfamily = 'sans', cat.pos = c(0, 0, 300, 220, 160),
  cat.dist = c(0.175, 0.2, 0.275, 0.2, 0.25), cex = 8, filename = NULL)

pdf(file = file.path(outputFolder, 'suppFig2.pdf'), width = 72, height = 72)
grid.newpage()
grid.draw(suppFig2)
dev.off()

#study CCDs
ccdDt = foreach(param = unique(combinedCors$params), .combine = rbind) %:%
  foreach(cond = c('control', 'perturbation'), .combine = rbind) %dopar% {
  
  genes = unique(combinedCors[params == param, gene1])
  
  ref = emat 
  
  if (cond == 'control') {
  
    esetGSE39445 = eset$GSE39445
    esetGSE48113 = eset$GSE48113
    esetGSE56931 = eset$GSE56931
    
    GSE39445_GSE48113 = calcCCD(esetGSE39445, esetGSE48113, genes)
    GSE39445_GSE56931 = calcCCD(esetGSE39445, esetGSE56931, genes)
    GSE48113_GSE56931 = calcCCD(esetGSE48113, esetGSE56931, genes)

    distDt = data.table(GSE39445_GSE48113, GSE39445_GSE56931, GSE48113_GSE56931,
                        params = param, cond = cond)
    
  } else {
    
    esetGSE39445 = esetPerturb$GSE39445
    esetGSE48113 = esetPerturb$GSE48113
    esetGSE56931 = esetPerturb$GSE56931
    
    esetGSE39445_ref = calcCCD(esetGSE39445, ref, genes)
    esetGSE48113_ref = calcCCD(esetGSE48113, ref, genes)   
    esetGSE56931_ref = calcCCD(esetGSE56931, ref, genes)         
    
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
  
  genes = unique(combinedCors[params == 'lambda_0.1101', gene1])
  
  ref = emat 
  
  if (cond == 'control') {
  
    esetGSE39445 = eset$GSE39445
    esetGSE48113 = eset$GSE48113
    esetGSE56931 = eset$GSE56931
    
    GSE39445_GSE48113 = calcCCD(esetGSE39445, esetGSE48113, genes, scale = FALSE)
    GSE39445_GSE56931 = calcCCD(esetGSE39445, esetGSE56931, genes, scale = FALSE)
    GSE48113_GSE56931 = calcCCD(esetGSE48113, esetGSE56931, genes, scale = FALSE)

    distDt = data.table(
      GSE39445_GSE48113
      , GSE39445_GSE56931
      , GSE48113_GSE56931
      , cond = cond)
    
  } else {
    
    esetGSE39445 = esetPerturb$GSE39445
    esetGSE48113 = esetPerturb$GSE48113
    esetGSE56931 = esetPerturb$GSE56931
    
    esetGSE39445_ref = calcCCD(esetGSE39445, ref, genes, scale = FALSE)
    esetGSE48113_ref = calcCCD(esetGSE48113, ref, genes, scale = FALSE)   
    esetGSE56931_ref = calcCCD(esetGSE56931, ref, genes, scale = FALSE)         
    
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

noClock = sm[is.na(clock_time), sample]
sm = sm[!is.na(clock_time)]
sm[, time := ztFrac * 24]

xClock = t(emat)[sm$sample, ]
fit = getModelFit(t(xClock), sm)
mFit = getPosteriorFit(fit, covMethod = 'data-driven')
mSamps = getPosteriorSamples(mFit)
rhyStats = getRhythmStats(mFit, fitType = 'posterior_mean', feature = glmnetGenes)

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
pGlmnetTimeCourseSamp = plotTimeCourse(glmnetTimeCourseDt[gene %in% sampGenes]
                                       , ncol = 1)

fig2 = pPeakPhase + pGlmnetTimeCourseSamp +
  plot_annotation(tag_levels = 'A') +
  plot_layout(ncol = 2, widths = c(2, 1))
ggexport(filename = file.path(outputFolder, 'fig2.pdf')
         , plot = fig2, width = 28, height = 16, units = 'in')

#saving glmnet correlations
ref = getCormat(emat, glmnetGenes, entrezID = TRUE)
diag(ref) = NA
qsave(ref, file.path(dataFolder, 'result_blood_ref.qs'))
saveRDS(ref, file.path(dataFolder, 'result_blood_ref.rds'))

refGenes = unique(glmnetCoefs[lambda == min(lambda),
                              .(entrez_hs = gene, symbol_hs = gene_sym)])
fwrite(refGenes, file.path(dataFolder, 'genes_blood.csv'))
