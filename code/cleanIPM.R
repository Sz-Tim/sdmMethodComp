# file cleaning
library(doSNOW); library(foreach)
dirs <- paste0("out/sp1/", c("none/none/", "noise/none/", 
                           "nonEq/none/", "sampBias/none/",
                           "none/noSB/", "none/underDisp/",
                           "none/overDisp/", "none/wrongCov/"))

p.c <- makeCluster(8); registerDoSNOW(p.c)
foreach(i=1:8) %dopar% {
  d <- dirs[i]
  fits <- dir(d, "IPM_fit", full.names=T)
  if(length(fits)>0) {
    for(f in 1:length(fits)) {
      if("try-error" %in% class(try(readRDS(fits[f]), silent=TRUE))) {
        file.remove(fits[f])
      }
    }
  }
}
stopCluster(p.c)