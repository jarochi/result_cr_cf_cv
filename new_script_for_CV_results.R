library(tidyverse)
library(readODS)
library(reshape2)

mydata <- read_ods("Example_CV.ods", col_names = FALSE)

# check number of rows, check number of separate tables
n_rows <- nrow(mydata)
n_tables <- n_rows/9
  
start_table <- seq(1, n_rows, 9)
end_table <- seq(9, n_rows, 9)

mediums <- c("LB10", "TSB", "BHI", "M63")

for (i in seq(1, n_tables)) {
  # read start and end of single table
  x <- start_table[i]
  y <- end_table[i]
  raw_data <- mydata[x:y,]
  # raw_data <- mydata[1:9,]     # TEST
  # grab conditions
  conditions <- raw_data[1,1]
  # transform single table
  strains_data <- raw_data[-c(2, 9), -c(1, 12)]
  strains <- unique(as.integer(strains_data[1,]))
  strains_results <- strains_data[-1,]
  strains_results <- strains_results %>% unlist %>% matrix(nrow = 6)

  mresults <- melt(strains_results, varnames = c("row", "col")) %>% 
    mutate(row = factor(row), 
           col = factor(col))
  
  # create plate scheme for single table
  strain_medium <- sapply(strains, function(x) (paste0(x, "-", mediums)))
  strain_medium <- rep(strain_medium, each = 3)
  plate_scheme <- matrix(strain_medium, nrow = 6)
  # transform plate scheme
  mplate_scheme <- melt(plate_scheme, varnames = c("row", "col"), value.name = "description") %>% 
    mutate(row = factor(row), 
           col = factor(col))
  
  # create single df from single table
  inner_join(mplate_scheme, mresults, by = c("row" = "row", "col" = "col")) %>% 
    filter(!is.na(description)) %>% 
    mutate(description = as.character(description), 
           strain = sapply(strsplit(description, split = "-"), first),
           medium = sapply(strsplit(description, split = "-"), last),
           temp = sapply(strsplit(conditions, split = " "), first), 
           replicate = rep(1:3, times = nrow(mresults)/3),
           surface = sapply(strsplit(conditions, split = " "), last))  %>%
    select(strain, medium, value, replicate, temp, surface) } %>% 
  bind_rows()  # bind rows doesnt work

  
