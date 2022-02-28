# Maslin2

tse_gen <- agglomerateByRank(tse, rank = "Genus")

# transpose assay
asv <- t(assay(tse_gen))
# assign OTU table with strain names
# it is strange that names weren't inherited
# upon TSE object generation (need to check why)
colnames(asv) <- rowData(tse_gen)$Genus
# store meta data for subset
meta_data <- data.frame(colData(tse_gen))

# run Maaslin2

fit_data <- Maaslin2(
  asv,
  meta_data,
  output = "./results/Maaslin2",
  transform = "AST",
  fixed_effects = c( "sleep_cat", "age", "bmi", "smoke"),  
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "sleep_cat,0",  
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)

# fit_data

# kable(head(filter(fit_data$results, qval <= 0.05)))