# Prepare the dolphin dataset for inclusion in the mfmsbm package
#
# The original data comes from the supplementary material of Geng et al. (2019)
# and is stored as dolphindata.RData with objects A (adjacency matrix) and
# Z0 (known community labels).
#
# To regenerate data/dolphin.rda, place dolphindata.RData in data-raw/ and run:
#   source("data-raw/prepare-dolphin.R")

if (file.exists("data-raw/dolphindata.RData")) {
  load("data-raw/dolphindata.RData")
} else {
  stop("Place dolphindata.RData in data-raw/ before running this script.")
}

diag(A) <- 0

Z0 <- c(1, 2, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2,
         1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1,
         1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1,
         2, 1)

dolphin <- list(A = A, Z0 = Z0)

usethis::use_data(dolphin, overwrite = TRUE)
