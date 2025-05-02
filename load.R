## Use devtools to load most recent version
library(utils)
library("devtools")
devtools:::load_all(export_all = FALSE)
setDTthreads(1)
