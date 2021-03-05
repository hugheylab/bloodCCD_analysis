library(glmnet)
library(Biobase)
library(doParallel)
library(data.table)
library(lubridate)
library(qs)
library(ggplot2)


GstudyMetadataPath = file.path('data', 'metadata', 'study_metadata.csv')
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


nFold = 10
alph = 0.9
lamb = 1:100/100

set.seed(21056)
foldIds = unique(sampleMetadata[, .(study, subject)])
foldIds[, foldId := sample(rep_len(1:nFold, .N))]  

sm = merge(sampleMetadata, foldIds, by = c('study', 'subject'))
noClock = sm[is.na(clock_time), sample]
sm = sm[!is.na(clock_time)]


X = t(emat)[!(rownames(t(emat)) %in% noClock), ]
Y = cbind(sm[, cos(2*pi*ztFrac)], sm[, sin(2*pi*ztFrac)])


glmnetCv = cv.glmnet(x = X, y = Y, lambda = lamb, alpha = alph
  , foldid = sm[, foldId], family = 'mgaussian')
beta1se = coef(glmnetCv, s = 'lambda.1se')
betaMin = coef(glmnetCv, s = 'lambda.min')

all.equal(rownames(beta1se[[1]])[beta1se[[1]][, 1] != 0],
rownames(betaMin[[1]])[beta1se[[1]][, 1] != 0])

