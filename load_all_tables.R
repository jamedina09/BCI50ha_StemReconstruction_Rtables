rm(list = ls())

library(data.table)

for (i in 1:9) {
    load(paste0("./BCI_stem_reconstruction/DATA/RTABLES/bci.stem", i, ".Rdata"))
    ## Combine all the data.tables into one
    if (i == 1) {
        bci.stem <- get(paste0("bci.stem", i))
    } else {
        bci.stem <- rbind(bci.stem, get(paste0("bci.stem", i)))
    }
    ## clean up the environment
    rm(list = paste0("bci.stem", i))
}

bci.stem <- as.data.table(bci.stem)

# 460211

length(unique(bci.stem$sp))

fwrite(data.frame(sp = sort(unique(bci.stem$sp))), "./sp.csv")
