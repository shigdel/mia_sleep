# load libraries
library(RColorBrewer)
library(readr)
library(openxlsx)
library(tidyverse)
library(microbiome)
library(vegan) # adonis
library(RVAideMemoire) # Pairwise permanova
library(compositions)
library(magrittr)
library(qwraps2)
library(dplyr)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ANCOMBC)
library(DT)
library(readxl)
library(lubridate)
library(haven)
library(mia)
library(kableExtra)
library(DT)
library(openxlsx)
library(gtsummary)
library(ggplot2)
library(scater)
library(patchwork)
library(ALDEx2)
library(Maaslin2)
library(dplyr)
library(rmarkdown)

knitr::opts_chunk$set(message = FALSE, warning = FALSE)

############## colData #####################

sleep_metadata = read.xlsx("Rh_sleep.xlsx",
                           sep = ",", colNames = TRUE)

# calculate age by using the date of birth and date of participating 

sleep_metadata <- sleep_metadata %>% 
  mutate(date_parti = make_datetime(qyearo, qmontho,qdayo),
         date_birth = make_datetime(byearo, bmontho,bdayo))

sleep_metadata$N_age <- as.numeric(difftime(sleep_metadata$date_parti, sleep_metadata$date_birth, units = "weeks"))/52.25


# ls(sleep_metadata)


selected <- dplyr::select(sleep_metadata, SampleID, N_age, gender,
                         adult,  sexo, heighto, weighto, gero,
                         hoursleepo,emao, dmso, diso, "osaso.y","osaso.x" , educationo, antibiotics, "snoringo.x","snoringo.y","smoke")




# table (selected$adult)


# rename the variables 

selected <- selected%>%mutate(bmi = as.numeric(weighto) / (as.numeric(heighto) / 100)^2, age = as.numeric(N_age), gender =sexo , 
                             height =heighto, weight=weighto, snoring= as.numeric(snoringo.x),
                             diso =as.numeric(diso), gero = as.numeric(gero))





#Heart burn and belching
selected$gero[selected$gero == 99] <- NA
selected$gero[selected$gero==1] <- "0"
selected$gero[selected$gero ==2]<- "1"
selected$gero[selected$gero == 3]<- "1"
selected$gero[selected$gero == 4]<- "1" # merge lower two group
selected$gero[selected$gero == 5]<- "1"


#snoring
selected$snoring[selected$snoring == 99] <- NA
selected$snoring[selected$snoring ==1] <- "0"
selected$snoring[selected$snoring ==2]<- "1"
selected$snoring[selected$snoring == 3]<- "1"
selected$snoring[selected$snoring == 4]<- "1"
selected$snoring[selected$snoring == 5]<- "1"



selected$sleep_cat[selected$snoring ==0 & selected$gero==0]=0

selected$sleep_cat[selected$snoring ==1 & selected$gero==0]=1

selected$sleep_cat[selected$snoring ==0 & selected$gero==1]=2

selected$sleep_cat[selected$snoring ==1 & selected$gero ==1]=3

selected$sleep_cat <- as.factor(selected$sleep_cat)

############### assays ##############

otu_data <-  read_tsv("feature-table.txt", skip = 1)

##remove the column otu since it is now used as a row name
otu_mat <- otu_data %>% dplyr::select (-"#OTU ID")

############### rowData ##############

## Taxonomy table
tax <- read.table(file = 'taxonomy.tsv', sep = '\t', header = TRUE)
taxtable<-tax %>% as.tibble() %>% separate(Taxon, sep="; ", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxtable

##  Idem for the two other matrices

row.names(taxtable) <- taxtable$Feature.ID
tax_mat <- taxtable %>% dplyr::select(-Feature.ID)

##Transform into matrixes otu and tax tables (sample table can be left as data frame)

otu_mat <- as.matrix(otu_mat)

tax_mat <- as.matrix(tax_mat)

# summarised experiment 

tse <- SummarizedExperiment(assays = list(counts = otu_mat),
                            colData = selected,
                            rowData = tax_mat)
#    rowTree = sleep_tree,
#        rowNodeLab = sleep_tree$node.label)
# Error: rowNodeLab should be provided.

# dim(tse)


# missing information on antibiotic use 
sum(sapply(tse$antibiotics, is.na))

# Including only adult participants 

tse <- tse[ , which(colData(tse)$adult==1)]

# Removing those who has used the antibiotics 

tse <- tse[ , which(colData(tse)$antibiotics!= "last4weeks")]

tse <- tse[ , which(colData(tse)$sleep_cat!= "NA")]

# dim (tse)

tse <- as(tse, "TreeSummarizedExperiment")

# tse
