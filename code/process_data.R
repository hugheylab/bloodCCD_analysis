library(metapredict)
library(sva)
library(data.table)
library(doParallel)
library(Biobase)
library(qs)
library(ggplot2)


parentFolderPath = file.path('data', 'expression_data') 

discoveryStudyNames = c('GSE39445', 'GSE48113', 'GSE56931')

studyMetadataPath = file.path('data', 'metadata', 'study_metadata.csv')
studyMetadata = fread(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = file.path('data', 'metadata', 'sample_metadata.csv')
sampleMetadata = fread(sampleMetadataPath, stringsAsFactors = FALSE)

esetList = getStudyDataList(parentFolderPath, studyMetadata)


#standardizing genewise expression values to mean zero, variance 1 
scaledEsetList = foreach(eset = esetList) %do% {
    
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


#gene-level plots
ggplot() +
  geom_violin(aes(x = ))


#individual-level plots


#study-level plots


#cross-study normalization
ematList = extractExpressionData(cbEsetList, sampleMetadata)
ematDiscovery = mergeStudyData(ematList[discoveryStudyNames], sampleMetadata)

qsave(cbEsetList, file = file.path('data', 'circadian_human_blood.qs'))
qsave(ematDiscovery, file = file.path('data', 'circadian_human_blood_emat.qs'))
