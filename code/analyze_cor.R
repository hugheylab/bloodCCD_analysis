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

theme_set(theme_bw())

codeFolder = file.path('code')
outputFolder = file.path('output')
dataFolder = file.path('data')

source(file.path(codeFolder, 'cor_utils.R'))

studyMetadataPath = file.path('data', 'metadata', 'study_metadata.csv')
studyMetadata = fread(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = file.path('data', 'metadata', 'sample_metadata.csv')
sampleMetadata = fread(sampleMetadataPath, stringsAsFactors = FALSE)

#loading control data
ematPath = file.path('data', 'circadian_human_blood_emat.qs')
emat = qread(ematPath)

esetPath = file.path('data', 'circadian_human_blood.qs')
eset = qread(esetPath)

#loading perturbation data
esetPerturbPath = file.path('data', 'perturb_esetList.qs')
esetPerturb = qread(esetPerturbPath)

ematPerturbPath = file.path('data', 'perturb_emat.qs')
ematPerturb = qread(ematPerturbPath)

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
zzCormat = foreach(absv = finalSumAbsv, .combine = rbind) %dopar% {
  
  vCo = zzCoefs[sumabsv == absv]
 
  zCormat = getCormat(emat, vCo$gene_sym)
  
  zCormatSort = sortCormat(zCormat)
  zCormatSort[, sumabsv := absv]
  
  return(zCormatSort)}

zzHeatmap =  plotHeatmap(zzCormat, sumabsv) +
  ggtitle(glue('Overall correlation heatmaps, genes selected by zeitzeiger'
               , ', by sumabsv with nSpc = 2')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 16))
  
#glmnet
glmnetCormat = foreach(lam = unique(glmnetCoefs$lambda)
  , .combine = rbind) %dopar% {
    
    vCo = glmnetCoefs[lambda == lam]
    
    gnetCormat = getCormat(emat, vCo$gene_sym)
    
    gnetCormatSort = sortCormat(gnetCormat)
    gnetCormatSort[, lambda := lam]
    
    return(gnetCormatSort)}
qsave(glmnetCormat, file.path(dataFolder, 'glmnet_cor_dt.qs'))

glmnetHeatmap = plotHeatmap(glmnetCormat, lambda) +
  ggtitle(bquote('Overall correlations, genes selected by glmnet, by '~lambda)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 16))

#2017
cormat2017 = getCormat(emat, genes2017Coefs$gene_sym)
cormat2017 = sortCormat(cormat2017)

heatmap2017 =  plotHeatmap(cormat2017) +
  ggtitle('Overall correlation matrix, genes selected by zeitzeiger, 2017') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 16))
heatmap2017 = plot_grid(heatmap2017, NULL, ncol = 2)

mergedFig = ggarrange(plotlist = list(zzHeatmap, glmnetHeatmap, heatmap2017)
  , nrow = 3)
ggexport(mergedFig, filename = file.path(outputFolder, 'overall_correlations.pdf')
  , width = 36, height = 36, unit = 'in')

###study data-sets
# studyEsetList = qread(file.path(dataFolder, 'subj_norm_esetList.qs'))
# studyEsetListPerturb = qread(file.path(dataFolder, 'subj_norm__pert_esetList.qs'))

#zeitzeiger
zzStudyCormat = foreach(i = 1:length(eset)
  , .combine = rbind) %dopar% {
              
  esetTmp = eset[[i]]  
    
  ematTmp = exprs(esetTmp)  
  
  foreach(absv = unique(zzCoefs$sumabsv), .combine = rbind) %dopar% {
                   
    vCo = zzCoefs[sumabsv == absv]
    
    zCormat = getCormat(ematTmp, vCo$gene_sym)
    
    zCormatSort = sortCormat(zCormat)
    zCormatSort[
      , `:=`(study = names(eset)[i]
        , sumabsv = absv)]
  
    return(zCormatSort)}}

zzStudyHeatmap = plotHeatmap(zzStudyCormat, sumabsv, study) +
  ggtitle(glue('Correlation heatmaps, genes selected by zeitzeiger'
               , ', by sumabsv with nSpc = 2')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 20))

#glmnet
glmnetStudyCormat = foreach(i = 1:length(eset)
  , .combine = rbind) %:% 
  foreach(cond = c('control', 'perturb'), .combine = rbind) %dopar% {
              
  if(cond == 'control') {
   
     esetTmp = eset[[i]]
    
  } else {esetTmp = esetPerturb[[i]]}
    
  foreach(lam = unique(glmnetCoefs$lambda), .combine = rbind) %dopar% {
                   
    geneSummGnet = glmnetCoefs[lambda == lam]
                             
    gnetCormat = getCormat(esetTmp, unique(glmnetCoefs$gene))
    
    gnetCormatSort = sortCormat(gnetCormat)
    gnetCormatSort[
      , `:=`(study = names(eset)[i]
        , lambda = lam
        , condition = cond)]
                             
    return(gnetCormatSort)}}

glmnetStudyHeatmap =  plotHeatmap(glmnetStudyCormat[condition == 'control']
                                  , lambda
                                  , study) +
  ggtitle('Correlations by study, genes selected by glmnet') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 8))
 
#### perturb vs control
pList = foreach(lam = unique(glmnetStudyCormat$lambda)) %dopar% {
  p =  plotHeatmap(glmnetStudyCormat, condition, study) +
    ggtitle(bquote('Correlations by study and condition, genes selected by glmnet, '~lambda[.(lam)])) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
          , text = element_text(size = 8))
  
  return(p)}
condFig = ggarrange(plotlist = pList, nrow = 2)
ggexport(condFig, filename = file.path(outputFolder, 'study_cond_corr.pdf')
  , width = 18, height = 18, units = 'in')

#2017
study2017Cormat = foreach(i = 1:length(eset)
  , .combine = rbind) %dopar% {
                              
  esetTmp = eset[[i]]  
  
  cormat = getCormat(esetTmp, unique(genes2017Coefs$gene_sym))
  
  cormatSort = sortCormat(cormat)
  cormatSort[, study := names(eset)[i]]
  
  return(cormatSort)}

heatmapStudy2017 =  plotHeatmap(study2017Cormat, study) + 
  ggtitle('Correlations by study, genes selected by zeitzeiger, 2017') +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
        , text = element_text(size = 20))

studyFig = ggarrange(plotlist = list(zzStudyHeatmap, glmnetStudyHeatmap
  , heatmapStudy2017), nrow = 3, heights = c(2, 2, 1.1))
ggexport(studyFig, filename = file.path(outputFolder, 'study_correlations.pdf')
         , width = 48, height = 60, unit = 'in')

#overlap
combinedCors = unique(glmnetCormatMelt[
  , .(gene1
      , gene2
      , lambda)])
setnames(combinedCors, old = 'lambda', new = 'params')
combinedCors[
  , `:=`(params = paste0('lambda_', round(params, 4))
         , model = 'glmnet')]
combinedCors = rbind(combinedCors
  , unique(zzCormatMelt[
    , .(gene1
        , gene2
        , params = paste0('sumabsv_', sumabsv)
        , model = 'zeitzeiger')]))
combinedCors = rbind(combinedCors
  , unique(cormatMelt2017[
    , .(gene1
        , gene2
        , params = '2017'
        , model = 'zeitzeiger')]))
setorderv(combinedCors, cols = c('gene1', 'gene2', 'model', 'params'))

vennDt = dcast(unique(combinedCors[, .(model, params, gene_sym = gene1)])
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

    distDt = data.table(params = param
      , GSE39445_GSE48113
      , GSE39445_GSE56931
      , GSE48113_GSE56931
      , cond = cond)
    
  } else {
    
    esetGSE39445 = esetPerturb$GSE39445
    esetGSE48113 = esetPerturb$GSE48113
    esetGSE56931 = esetPerturb$GSE56931
    
    esetGSE39445_ref = calcCCD(esetGSE39445, ref, genes)
    esetGSE48113_ref = calcCCD(esetGSE48113, ref, genes)   
    esetGSE56931_ref = calcCCD(esetGSE56931, ref, genes)         
    
    distDt = data.table(params = param
      , esetGSE39445_ref
      , esetGSE48113_ref
      , esetGSE56931_ref
      , cond = cond)}
  
  distDt = melt(distDt, id.vars = c('params', 'cond'), variable.name = 'studies'
    , value.name = 'ccd')
  
  return(distDt)}

pCcd = ggplot(data = ccdDt) +
  geom_point(aes(x = params, y = ccd, fill = cond)
    , position = position_dodge(width = 0.3), shape = 21) +
  geom_text_repel(aes(x = params, y = ccd, label = studies), size = 2) +
  scale_fill_brewer(type = 'div', direction = -1) +
  coord_flip() +
  ggtitle('Range of CCDs by parameters and conditions, scaled by number of gene pairs')

ggsave(filename = file.path(outputFolder, 'ccd_plot.pdf'), pCcd 
  , height = 12, width = 12, units = 'in')