source(file.path('code', 'utils.R'))

parentFolderPath = file.path(dataFolder, 'expression_data') 

discoveryStudyNames = c('GSE39445', 'GSE48113', 'GSE56931')

studyMetadataPath = file.path('data', 'metadata', 'study_metadata.csv')
studyMetadata = fread(studyMetadataPath)

sampleMetadataPath = file.path('data', 'metadata', 'sample_metadata.csv')
sampleMetadata = fread(sampleMetadataPath)
perturbConds = c('Sleep Restriction', 'Out of phase with respect to melatonin', 
                 'sleep deprivation')
perturbMetadata = sampleMetadata[condition %in% perturbConds]

esetList = getStudyDataList(parentFolderPath, studyMetadata)

#limiting to perturbation samples
perturbEsetList = foreach (eset = esetList) %do% {
  
  eset = eset[, sampleNames(eset) %in% perturbMetadata$sample]

  #standardizing genewise expression values to mean zero, variance 1 
  exprs(eset) = t(scale(t(exprs(eset))))
  
  #combat for individual-level batch correction 
  cbEset = eset
  pheno = pData(eset)
  edata = exprs(eset)
  batch = perturbMetadata[sample %in% pheno[, 'geo_accession'], subject]
  
  modCombat = model.matrix(~1, data = pheno)
  combatEdata = ComBat(edata, batch = batch, mod = modCombat)
  
  exprs(cbEset) = combatEdata
  
  return(cbEset)}
names(cbEsetList) = names(esetList)

#cross-study normalization
ematList = extractExpressionData(cbEsetList, perturbMetadata)
ematDiscovery = mergeStudyData(ematList[discoveryStudyNames], perturbMetadata)

qsave(cbEsetList, file = file.path('data', 'perturb_esetList.qs'))
qsave(ematDiscovery, file = file.path('data', 'perturb_emat.qs'))
