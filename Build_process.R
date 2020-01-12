library(usethis)
library(devtools)
usethis::use_build_ignore("Build_process.R")
.libPaths()
lapply(.libPaths(), dir)
# install.packages("styler")
Authors@mpraprofiler: person("Xiaoming", "Lu", email = "xiaoming-lu@hotmail.com", role = c("aut", "cre"))
usethis::use_mit_license("Xiaoming Lu")


devtools::document()
devtools::load_all()

usethis::use_package("DESeq2")
usethis::use_package("methods")
usethis::use_package("stats")
usethis::use_package("dplyr")
usethis::use_package("tibble")
usethis::use_package("tidyr")
usethis::use_package("base")

use_r("fold_enhancer_plot")
use_r("allgeneric")
use_r("allelic_analysis")
use_r("allelic_enhancer_dot_plot")

install()

 devtools::load_all()

devtools::document()

?allelic_compare
devtools::check()
session_info()

use_readme_rmd()
use_data_raw()


