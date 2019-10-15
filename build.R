library("rstantools")
library("testthat")
options(mc.cores = parallel::detectCores())
## args(rstan_create_package)
## rstan_create_package(path = "rstanOcc", rstudio=FALSE)

## list.files(all.files = TRUE)
## file.show("Read-and-delete-me")
## file.remove("Read-and-delete-me")

## Use the next step to compile the code 
## pkgbuild::compile_dll() # see note below

## The next step creates documentation 

## file.remove("NAMESPACE")
devtools::document()
devtools::install(quick = TRUE)

## only run next line to recompile everything, which takes a while 
devtools::install()



