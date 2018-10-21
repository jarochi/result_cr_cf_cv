library(dplyr)
library(readODS)
library(reshape2)

# mydata <- read_ods("Example_CV.ods", col_names = FALSE)
mydata <- read_ods("measure1.ods", col_names = FALSE)

# check number of rows, check number of separate tables
n_rows <- nrow(mydata)
n_tables <- ceiling(n_rows/9)

start_table <- seq(1, n_rows, 9)
end_table <- start_table + 8
mediums <- c("LB10", "TSB", "BHI", "M63")

# przerobić
# res <- lapply(split(mydata, ceiling(1L:nrow(mydata)/9)), function(ith_table) {
#   raw_data <- as.matrix(ith_table[-1, ])
#   class(raw_data) <- "numeric"
#   colnames(raw_data) <- ith_table[1, ]
#   rownames(raw_data) <- NULL
# 
#   rbind(cbind(melt(raw_data[, -c(1, 12)]), mediums = unlist(lapply(mediums, rep, 4))),
#         cbind(melt(raw_data[, c(1, 12)]), mediums = unlist(lapply(mediums, rep, 2)))) %>%
#     select(-Var1) %>%
#     rename(strain = Var2)
# }) %>% bind_rows()


# empty list to override bind_rows
datalist = list()

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
  datalist[[i]] <- inner_join(mplate_scheme, mresults, by = c("row" = "row", "col" = "col")) %>% 
    filter(!is.na(description)) %>% 
    mutate(description = as.character(description), 
           strain = sapply(strsplit(description, split = "-"), first),
           medium = sapply(strsplit(description, split = "-"), last),
           temp = sapply(strsplit(conditions, split = " "), first), 
           replicate = rep(1:3, times = nrow(mresults)/3),
           surface = sapply(strsplit(conditions, split = " "), last))  %>%
    select(strain, medium, value, replicate, temp, surface) 
} 
# %>% 
#   bind_rows()  # bind rows doesnt work


all_results <- do.call(rbind, datalist)



library(ggplot2)
library(ggbeeswarm)

ggplot(all_results, aes(x = strain, y = value, color = temp, shape = surface)) +
  geom_boxplot() +
  facet_wrap(~ medium) +
  theme_bw()


ggplot(all_results, aes(x = strain, y = value, color = temp, shape = surface)) +
  geom_quasirandom() +
  facet_wrap(~ medium) +
  theme_bw()

group_by(all_results, strain, medium, temp, surface, replicate) %>%
  summarise(value = median(value)) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = strain, y = value, color = temp, shape = surface)) +
  geom_point(size = 3) +
  facet_wrap(~ medium) +
  theme_bw()

# group_by(all_results, strain, medium, temp, surface, replicate) %>%
#   mutate(value = value/max(value)) %>%
#   ggplot(aes(x = strain, y = value, color = temp, shape = surface)) +
#   geom_quasirandom() +
#   facet_wrap(~ medium) +
#   theme_bw()
