library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)