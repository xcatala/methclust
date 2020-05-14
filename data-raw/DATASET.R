## code to prepare `DATASET` dataset goes here

library(minfi)

test_data <- minfiData::MsetEx %>%
  mapToGenome() %>%
  ratioConvert()

usethis::use_data(test_data, overwrite = TRUE, compress = 'xz')
