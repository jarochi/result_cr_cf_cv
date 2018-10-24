library(dplyr)
library(readODS)
library(reshape2)

# mydata <- read_ods("Example_CV.ods", col_names = FALSE)
mydata <- read_ods("measure1.ods", col_names = FALSE)
# plate_scheme <- read.csv("plate_scheme.csv")

start_table <- seq(1, nrow(mydata), 9)
end_table <- start_table + 8

mediums <- c("LB10", "BHI", "TSB", "M63")

res <- lapply(split(mydata, ceiling(1L:nrow(mydata)/9)), function(ith_table) {
  raw_data <- as.matrix(ith_table[-1, ])

  class(raw_data) <- "numeric"
  colnames(raw_data) <- ith_table[1, ]
  rownames(raw_data) <- NULL
  
  conditions <- strsplit(ith_table[1,1], split = " ") %>% unlist()
  # browser()
  # 
  # raw_data %>% matrix(nrow = 8)
  # 
  # raw_data <- raw_data %>%  mutate(temp = conditions[1],
  #        surface = conditions[2])

  not_NA_rows <- !apply(raw_data, 1, function(i) mean(is.na(i))) > 0.7
  not_NA_cols <- !apply(raw_data, 2, function(i) mean(is.na(i))) > 0.7
  if(sum(dim(raw_data[not_NA_rows, not_NA_cols])) == 0) {
    data.frame()
  } else {
    cbind(melt(raw_data[not_NA_rows, not_NA_cols]), medium = unlist(lapply(mediums, rep, 3))) %>%
      colnames(seq(1, nrow(raw_data))) %>% 
      inner_join(mplate_scheme, by = c("row" = "row", "col" = "col"))
      select(-Var1) %>%
      rename(strain = Var2)
    browser()
  }
}) %>%  bind_rows()



# empty list to override bind_rows
datalist = list()

for (i in seq(1, ceiling(nrow(mydata)/9))) {
  # read start and end of single table
  x <- start_table[i]
  y <- end_table[i]
  raw_data <- mydata[x:y,]
  # raw_data <- mydata[1:9,] #TEST
  # transform single table
  strains_data <- raw_data[-c(2, 9), -c(1, 12)]
  # strains_results <- strains_data[-1,] %>% unlist %>% matrix(nrow = 6)
  mresults <- melt(strains_data[-1,] %>% unlist %>% matrix(nrow = 6), varnames = c("row", "col")) %>% 
    mutate(row = factor(row), 
           col = factor(col))
  
  # create and transform plate scheme for single table
  plate_scheme <- sapply(unique(as.integer(strains_data[1,])), function(x) (paste0(x, "-", mediums))) %>% 
    rep(each = 3) %>% 
    matrix(nrow = 6)
  mplate_scheme <- melt(plate_scheme, varnames = c("row", "col"), value.name = "description") %>% 
    mutate(row = factor(row), 
           col = factor(col))
  
  # create single df from single table
  conditions <- raw_data[1,1] # grab conditions
  datalist[[i]] <- inner_join(mplate_scheme, mresults, by = c("row" = "row", "col" = "col")) %>% 
    filter(!is.na(description)) %>% 
    mutate(description = as.character(description), 
           strain = sapply(strsplit(description, split = "-"), first),
           medium = sapply(strsplit(description, split = "-"), last),
           temp = sapply(strsplit(conditions, split = " "), first), 
           replicate = rep(1:3, times = nrow(mresults)/3),
           surface = sapply(strsplit(conditions, split = " "), last))  %>%
    select(strain, medium, value, replicate, temp, surface) %>% 
    filter(!is.na(value))
} 

all_results <- do.call(rbind, datalist)


library(ggplot2)
library(ggbeeswarm)

all_results %>% mutate(conditions = paste(temp, surface)) %>% 
ggplot(aes(x = strain, y = value, color = conditions)) +
  geom_boxplot() +
  facet_wrap(~ medium) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_results %>% mutate(conditions = paste(temp, surface)) %>% 
ggplot(aes(x = strain, y = value, color = conditions)) +
  geom_quasirandom() +
  facet_wrap(~ medium) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_results %>% mutate(conditions = paste(temp, surface)) %>% 
group_by(strain, medium, temp, surface, replicate) %>%
  summarise(value = median(value)) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = strain, y = value, color = temp, shape = surface)) +
  geom_point(size = 3) +
  facet_wrap(~ medium) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# group_by(all_results, strain, medium, temp, surface, replicate) %>%
#   mutate(value = value/max(value)) %>%
#   ggplot(aes(x = strain, y = value, color = temp, shape = surface)) +
#   geom_quasirandom() +
#   facet_wrap(~ medium) +
#   theme_bw()


OD <- all_results %>% group_by(strain, medium, temp, surface) %>% summarise(avg = mean(value))
ODc <- 0.075 # temporary, need better csv preparation

formers <- as.character(sapply(OD$avg, function(x) {
  cut(x, breaks = c(0, ODc, 2*ODc, 4*ODc, 10), labels = c("no biofilm", "weak", "moderate", "strong"))
}))

OD <- bind_cols(OD, data_frame(formers))


OD %>% mutate(conditions = paste(temp, surface)) %>% 
ggplot(aes(x = strain, y = formers, fill = conditions)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(~ medium) +
  theme_bw() +
  scale_fill_discrete("Conditions") +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


