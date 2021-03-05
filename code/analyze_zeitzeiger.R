library(zeitzeiger)
library(Biobase)
library(doParallel)
library(data.table)
library(qs)
library(lubridate)

studyMetadataPath = file.path('data', 'metadata', 'study_metadata.csv')
studyMetadata = fread(studyMetadataPath, stringsAsFactors = FALSE)

sampleMetadataPath = file.path('data', 'metadata', 'sample_metadata.csv')
sampleMetadata = fread(sampleMetadataPath, stringsAsFactors = FALSE)

ematPath = file.path('data', 'circadian_human_blood_emat.qs')
emat = qread(ematListPath)

timeMax = 24
sampleMetadata[
  , zt := as.duration(hm(clock_time) - hm(sunrise_time))/as.duration(hours(1))
  ][zt < 0, zt := zt + timeMax]
sampleMetadata[
  , ztFrac := zt/timeMax
  ][, zt := NULL]


nFold = 10
sumabsv = c(1, 2, 3)
nTime = 12 
nSpc = 1:3

set.seed(21056)
foldIds = unique(sampleMetadata[, .(study, subject)])
foldIds[, foldId := sample(rep_len(1:nFold, .N))]  

sm = merge(sampleMetadata, foldIds, by = c('study', 'subject'))
noClock = sm[is.na(clock_time), sample]
sm = sm[!is.na(clock_time)]

xClock = t(emat)[!(rownames(t(emat)) %in% noClock), ]

fitResultList = zeitzeigerFitCv(xClock, sm[, ztFrac], sm[, foldId])

spcResultList = foreach(absv = sumabsv) %do% {
  
  spcResult = zeitzeigerSpcCv(fitResultList, nTime = nTime, sumabsv = absv)
  return(spcResult)}

