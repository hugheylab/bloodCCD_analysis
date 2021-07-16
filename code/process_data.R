library(metapredict)
library(sva)
library(data.table)
library(doParallel)
library(Biobase)
library(qs)
library(ggplot2)
theme_set(theme_bw())

dataFolder= 'data'

parentFolderPath = file.path(dataFolder, 'expression_data') 

discoveryStudyNames = c('GSE39445', 'GSE48113', 'GSE56931')

studyMetadataPath = file.path(dataFolder, 'metadata', 'study_metadata.csv')
studyMetadata = fread(studyMetadataPath)

sampleMetadataPath = file.path(dataFolder, 'metadata', 'sample_metadata.csv')
sampleMetadata = fread(sampleMetadataPath)
controlConds = c('Sleep Extension', 'In phase with respect to melatonin', 
                 'baseline')
controlMetadata = sampleMetadata[condition %in% controlConds]

esetList = getStudyDataList(parentFolderPath, studyMetadata)



cbEsetList = foreach(eset = esetList) %do% {
  #limiting to control samples 
  eset = eset[, sampleNames(eset) %in% controlMetadata$sample]
  
  #standardizing genewise expression values to mean zero, variance 1  
  exprs(eset) = t(scale(t(exprs(eset))))
  
  #combat for individual-level batch correction  
  pheno = pData(eset)
  edata = exprs(eset)
  batch = controlMetadata[sample %in% pheno[, 'geo_accession'], subject]
  
  modCombat = model.matrix(~1, data = pheno)
  combatEdata = ComBat(edata, batch = batch, mod = modCombat)
  
  exprs(cbEset) = combatEdata
  
  return(cbEset)}
names(cbEsetList) = names(esetList)

#long data for plots 
eDt = foreach(eset = cbEsetList, .combine = rbind) %do% {
  
  emat = exprs(eset)
  names(dimnames(emat)) = c('gene', 'sample')
  
  ematDt = as.data.table(as.table(emat))
  setnames(ematDt, 'N', 'expression')
  return(ematDt)}

eDt = merge(controlMetadata[, .(study, subject, sample)], 
            eDt, by = 'sample')

#gene-level plots
p1 = ggplot(eDt, aes(x = gene, y = expression)) +
  geom_boxplot()

#gene expression faceted by study
p2 = ggplot(eDt, aes(x = gene, y = expression)) +
  geom_boxplot() +
  facet_grid(~study)
   
#gene expression by subject
p3 = ggplot(eDt, aes(x = subject, y = expression)) +
  geom_boxplot() +
  facet_wrap(~study, scales = 'free_x')

#study-level expression
p4 = ggplot(eDt, aes(x = study, y = expression)) +
  geom_boxplot()


#cross-study normalization
ematList = extractExpressionData(cbEsetList, sampleMetadata)
ematDiscovery = mergeStudyData(ematList[discoveryStudyNames], controlMetadata)


#post-cross-study normalization long data
names(dimnames(ematDiscovery)) = c('gene', 'sample')
  
ematDiscoveryDt = as.data.table(as.table(ematDiscovery))
setnames(ematDiscoveryDt, 'N', 'expression')

ematDiscoveryDt = merge(controlMetadata[, .(study, subject, sample)]
  , ematDiscoveryDt, by = 'sample')


#gene-level plots
p5 = ggplot(ematDiscoveryDt, aes(x = gene, y = expression)) +
  geom_boxplot()

#gene expression faceted by study
p6 = ggplot(ematDiscoveryDt, aes(x = gene, y = expression)) +
  geom_boxplot() +
  facet_grid(~study)
   
#gene expression by subject
p7 = ggplot(ematDiscoveryDt, aes(x = subject, y = expression)) +
  geom_boxplot() +
  facet_wrap(~study, scales = 'free_x')

#study-level expression
p8 = ggplot(ematDiscoveryDt, aes(x = study, y = expression)) +
  geom_boxplot()

qsave(cbEsetList, file = file.path(dataFolder, 'circadian_human_blood.qs'))
qsave(ematDiscovery, file = file.path(dataFolder, 'circadian_human_blood_emat.qs'))
