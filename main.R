# generate tse object
source("import_data.R")
# run alpha diversity analysis
rmarkdown::render("alpha.Rmd",
                  output_format = "pdf_document",
                  output_file = "./results/alpha.pdf")
# run beta diversity analysis
rmarkdown::render("beta.Rmd",
                  output_format = "pdf_document",
                  output_file = "./results/beta.pdf")
# run aldex2 DA analysis
rmarkdown::render("run_aldex2.Rmd",
                  output_format = "pdf_document",
                  output_file = "./results/aldex2.pdf")
# run maaslin2 DA analysis
# maaslin2 stores output in a directory, so Rmarkdown is probably not necessary
source("run_maaslin2.R")
# run ancombc DA analysis
rmarkdown::render("run_ancombc.Rmd",
                  output_format = "pdf_document",
                  output_file = "./results/ancombc.pdf")
