args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(xcms))
suppressMessages(library(gtools))

date()

input_dir <- paste0(args[1])

ionMode <- paste0(args[2])
ext <- ".rds"

# List all peak picked files and sort order by file number
input <- mixedsort(list.files(input_dir, pattern = glob2rx(paste0("xcms1-",ionMode,"*",ext)), full.names = T),decreasing=T)

# Combine into single object
input_l <- lapply(input, readRDS)
xset <- input_l[[1]]
for(i in 2:length(input_l)) {
  set <- input_l[[i]]
  xset <- c(xset, set)
  print(i)
  }
rm(input,input_l)

# Save as R object
save(list=c("xset"), file = paste0(input_dir,"/xset-",ionMode,".RData"))
