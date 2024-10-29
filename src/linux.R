#--------------------------------------------------
# lỉbrary
library(unix)
library("BiocParallel")
library("doParallel")
# width of display in terminal
options(width=190)


# funtion --------------------------------------------
# convert bytes to gigabytes
bytes_to_gb <- function(bytes) {
  gb <- bytes*0.95 / 1e9    # only use 95% of limit
  return(gb)
}
# convert gigabytes to bytes
gb_to_bytes <- function(gb) {
  bytes <- gb * 1e9 * 0.95
  return(bytes)
}


# setting core --------------------------------------------------
# setting memory and core usage
mem=gb_to_bytes(as.numeric(gsub("[^0-9]", "", Sys.getenv("R_MAX_VSIZE"))))
core=as.numeric(gsub("[^0-9]", "", Sys.getenv("R_MEMORY_LIMIT")))
print(paste0('mem:', mem))
print(paste0('core:', core))


# set resource limits
rlimit_as(mem)
rlimit_core(core)
print(rlimit_all())

# register for MulticoreParam
if (is.na(core)) {
  print('Cannot detect core')
} else {
  register(MulticoreParam(core))
  registerDoParallel(core) 
  register(DoparParam(), default = TRUE)
}

