args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(dplyr))

# Read in list of all mz files and specifiy (in sbatch command) the ion Mode
usePath <- paste0(args[1])
ext <- ".mzML"
pre <- paste0(usePath,"/")

ionMode<-paste0(args[2])

mzdatafiles <- list.files(usePath,recursive = FALSE, full.names=TRUE, pattern=glob2rx(paste0("*",ionMode,"*",ext)))

# Read in list of all csv files from sequence methods and specificy ion Mode
csvfile <- list.files(usePath, recursive=FALSE, full.name=TRUE, pattern=glob2rx(paste0("*",ionMode,"*",".csv")))

# Combine all of your files into one
all_csv<-bind_rows(lapply(csvfile,read.csv,skip=1))

# This is just a blank row generated when exported that we don't need
all_csv$Sample.ID<-NULL

# Filter all to keep just the ionMode and goodData of interest
all_csv<-all_csv[which(all_csv$ionMode==ionMode & all_csv$goodData==1),]

# Add a column with the full file basename
all_csv$FileWithExtension<-paste0(all_csv$File.Name,ext)

# Check you didn't make any mistakes and then Write file
if(all(all_csv$FileWithExtension %in% basename(mzdatafiles))==TRUE){ write.table(all_csv,paste0("metadata_",ionMode,".txt"),append = FALSE, sep = "\t",row.names = FALSE,col.names=TRUE,quote=FALSE)}
