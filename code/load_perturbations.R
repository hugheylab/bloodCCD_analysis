library(metapredict)
library(sva)
library(data.table)
library(doParallel)
library(Biobase)
library(qs)


dataFolder= 'data'

parentFolderPath = file.path(dataFolder, 'expression_data') 

discoveryStudyNames = c('GSE39445', 'GSE48113', 'GSE56931')

studyMetadataPath = file.path('data', 'metadata', 'study_metadata.csv')
studyMetadata = fread(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = file.path('data', 'metadata', 'sample_metadata.csv')
sampleMetadata = fread(sampleMetadataPath, stringsAsFactors = FALSE)
perturbConds = c('Sleep Restriction', 'Out of phase with respect to melatonin'
  , 'sleep deprivation')
perturbMetadata = sampleMetadata[condition %in% perturbConds]

esetList = getStudyDataList(parentFolderPath, studyMetadata)

#limiting to perturbation samples
perturbEsetList = foreach (eset = esetList) %do% {
  
  esetPerturb = eset[, sampleNames(eset) %in% perturbMetadata$sample]
  return(esetPerturb)}

#standardizing genewise expression values to mean zero, variance 1 
scaledEsetList = foreach(eset = perturbEsetList) %do% {
  
  esetScaled = eset
  exprs(esetScaled) = t(scale(t(exprs(eset))))
  
  return(esetScaled)}
names(scaledEsetList) = names(esetList)

#combat for individual-level batch correction 
cbEsetList = foreach(eset = scaledEsetList) %do% {
  
  cbEset = eset
  pheno = pData(eset)
  edata = exprs(eset)
  batch = sampleMetadata[sample %in% pheno[, 'geo_accession'], subject]
  
  modCombat = model.matrix(~1, data = pheno)
  combatEdata = ComBat(edata, batch = batch, mod = modCombat)
  
  exprs(cbEset) = combatEdata
  
  return(cbEset)}
names(cbEsetList) = names(esetList)
qsave(cbEsetList, file = file.path(dataFolder, 'subj_norm__pert_esetList.qs'))

#cross-study normalization
ematList = extractExpressionData(cbEsetList, sampleMetadata)
ematDiscovery = mergeStudyData(ematList[discoveryStudyNames], sampleMetadata)

qsave(cbEsetList, file = file.path('data', 'perturb_esetList.qs'))
qsave(ematDiscovery, file = file.path('data', 'perturb_emat.qs'))
