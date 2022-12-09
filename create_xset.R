suppressMessages(library(xcms))

setwd("~/untarg_xcms/mzML_MSn2022/output_dir/xcms2")

# polarity mode
mode <- "neg"

# Load the MS OnDisk object combined in previous script
load(file=paste0("xcms2_final-",mode,".RData"))

# This function retrieve a xset like object and fixes the error of sample class naming, @author Gildas Le Corguille lecorguille@sb-roscoff.fr
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    if (class(xobject) == "xcmsSet")
        return (xobject)
    # XCMS 3.x
    if (class(xobject) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, 'xcmsSet'))
        if (is.null(xset@phenoData$sample_group))
            sampclass(xset) = "."
        else
            sampclass(xset) <- xset@phenoData$sample_group
        return (xset)
    }
}

# Create the xcmsSet object
xset <- getxcmsSetObject(processedData)
# Use fill peaks 
xset <- fillPeaks(xset, method = "chrom")
dim(xset@groups)

# Save
saveRDS(xset, file=paste0(mode,"_xset.rds"))
