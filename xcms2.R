args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(xcms))
suppressMessages(library(BiocParallel))
suppressMessages(library(gtools))

date()

# Input dir
input_dir <- paste0(args[1])

# Output dir
output_dir<- paste0(args[2])

# Ion Mode
ionMode<- paste0(args[3])

# Biocparallel setting
register(BPPARAM = MulticoreParam(workers=36))

# Metadata
meta <- read.delim(paste0("~/untarg_xcms/mzML_MSn2022/metadata_",ionMode,".txt"))

# Load the MS OnDisk object combined in previous script
load(file=paste0(input_dir,"/xset-",ionMode,".RData"))

# Add variable for subsetting
idx<-which(xset@phenoData$Sample.Name ==  paste0("BIOSSCOPE pool ",ionMode))
xset@phenoData$subset.name <- "sample"
xset@phenoData$subset.name[idx] <- "pool"

# RT correction
prm <- ObiwarpParam(subset= which(xset@phenoData$subset.name == "pool"), subsetAdjust="average", binSize = 0.1,distFun = "cor", gapInit = 0.3, gapExtend = 2.4)
xset_obi <- adjustRtime(xset, param = prm, msLevel = 1L)

save(list=c("xset_obi"), file = paste0(output_dir,"/xcms2_obi-",ionMode,".RData"))
rm(xset)
print("Completed xcms obiwarp")

# Add variable for grouping(subsetting)
idx<-which(xset_obi@phenoData$Sample.Name ==  paste0("BIOSSCOPE pool ",ionMode))
xset_obi@phenoData$subset.name <- "sample"
xset_obi@phenoData$subset.name[idx] <- "pool"

# Grouping
groupDiff <- 0.7
bsize <- 0.025
pdp<-PeakDensityParam(sampleGroups = xset_obi@phenoData$subset.name, minFraction = 0.1, minSamples = 1, bw = groupDiff, binSize = bsize)
xset_gc<-groupChromPeaks(xset_obi, param = pdp)
rm(xset_obi)
save(list=c("xset_gc"), file = paste0(output_dir,"/xcms2_gc-",ionMode,".RData"))
print("Completed xcms grouping")

# Fillpeaks
fillParam<-FillChromPeaksParam(expandMz = 0, expandRt = 0, ppm = 25)
xset_fp<-fillChromPeaks(xset_gc,fillParam)
#xset_fp <- xset_gc
rm(xset_gc)

# Save final product
processedData<-xset_fp
rm(xset_fp)
save(list=c("processedData"), file = paste0(output_dir,"/xcms2_final-",ionMode,".RData"))

# Output all peaks and save
allPeaks<-chromPeaks(processedData)
write.csv(allPeaks, file = paste0(output_dir,"/BATSuntarg_",ionMode,"_aligned.csv"))

# Output features and save
featuresDef<-featureDefinitions(processedData)
featuresIntensities<-featureValues(processedData, value = "into", method = "maxint")
dataTable<-merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
dataTable <-dataTable[, !(colnames(dataTable) %in% c("peakidx"))]
write.csv(dataTable, file = paste0(output_dir,"/BATSuntarg_",ionMode,"_picked.csv"))

# function to format for GNPS from Johannes Rainer
formatSpectraForGNPS <- function(x) {
    fids <- mcols(x)$feature_id
    if (!length(fids))
        stop("No column named 'feature_id' present in 'mcols(x)'")
    fids <- as.integer(sub("^FT", "", fids))
    mendoapply(x, fids, FUN = function(z, id) {
        z@acquisitionNum <- id
        z
    })
}

# export MS1 and MS2 features
filteredMs2Spec <- featureSpectra(processedData, return.type = "Spectra")
filteredMs2Spec <- clean(filteredMs2Spec, all = TRUE)
filteredMs2Spec <- formatSpectraForGNPS(filteredMs2Spec)
writeMgfData(filteredMs2Spec, paste0(output_dir,"/", ionMode, "_ms2spectra_all.mgf"))

# function to use maxTic with combineSpectra from Johannes Rainer
maxTic <- function(z) {
    z[[which.max(lapply(intensity(z), sum))]]
}

# export MS2 features only
filteredMs2Spec_maxTic <- combineSpectra(filteredMs2Spec, fcol = "feature_id", method = maxTic)
writeMgfData(filteredMs2Spec_maxTic, paste0(output_dir,"/",ionMode, "_ms2spectra_maxTic.mgf"))

# filter data table to contain only peaks with MSMS
filteredDT <-dataTable[which(dataTable$Row.names %in% filteredMs2Spec@elementMetadata$feature_id),]
write.table(filteredDT, file = paste0(output_dir,"/",ionMode, "_xcms_onlyMS2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# export MS2 features only - consensus
filteredMs2Spec_consensus <- combineSpectra(filteredMs2Spec, fcol = "feature_id", method = consensusSpectrum, mzd = 0, minProp = 0.8, ppm = 10)
writeMgfData(filteredMs2Spec_consensus, paste0(output_dir,"/",ionMode, "_ms2spectra_consensus.mgf"))

consensusDT <- dataTable[which(dataTable$Row.names %in% filteredMs2Spec@elementMetadata$feature_id),]
write.table(consensusDT, file = paste0(output_dir,"/",ionMode, "_xcms_consensusMS2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

