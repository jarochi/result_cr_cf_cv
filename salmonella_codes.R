library(dplyr)
library(readODS)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(tidyverse)

######################################################### Wczytanie danych

mydata <- read.csv("js.csv", header = FALSE, stringsAsFactors = FALSE)  # Asia

start_table <- seq(1, nrow(mydata), 9)
end_table <- start_table + 8

mediums <- c("LB10", "BHI", "TSB", "M63")

# empty list to override bind_rows
datalist = list()

for (i in seq(1, ceiling(nrow(mydata)/9))) {
  # read start and end of single table
  x <- start_table[i]
  y <- end_table[i]
  raw_data <- mydata[x:y,]
  # raw_data <- mydata[1:9,] #TEST
  # transform single table
  ODc <- 0.125 
  strains_data <- raw_data[-c(2, 9), -c(1, 12)]
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
  all_res <- inner_join(mplate_scheme, mresults, by = c("row" = "row", "col" = "col")) %>% 
    filter(!is.na(description)) %>% 
    mutate(description = as.character(description), 
           strain = sapply(strsplit(description, split = "-"), first),
           medium = sapply(strsplit(description, split = "-"), last),
           temp = sapply(strsplit(conditions, split = " "), first), 
           replicate = rep(1:3, times = nrow(mplate_scheme)/3),
           surface = sapply(strsplit(conditions, split = " "), last),
           plate_no = rep(i, times = nrow(mplate_scheme)))  %>%                       # DK liczpa powtórzeń
    # ) %>% 
    # rep(1:3, by = 12, each = 12, times = 4)
    # select(strain, medium, value, replicate, temp, surface) %>% 
    filter(!is.na(value))
  
  # add to table above strength of biofilm formation
  all_res$value <- as.numeric(as.character(all_res$value))
  
  OD <- all_res %>% group_by(strain, medium, description, temp, surface) %>% summarise(avg = mean(value)) %>% 
    mutate(strength = as.character(sapply(avg, function(x) {
      cut(x, breaks = c(0, ODc, 2*ODc, 4*ODc, 10), labels = c("absence", "weak", "moderate", "strong"))
    })))
  
  datalist[[i]] <- inner_join(all_res, OD, by = c("strain", "medium", "temp", "surface"))
} 

all_results <- do.call(rbind, datalist)

# Asia
all_results <- all_results %>% 
  mutate(rep_no = rep(1:3, each = 240, times = (nrow(all_results)/720)))

######################################################### koniec wczytywania danych

##### Wykresy

all_results %>% 
  group_by(strain, medium, rep_no, temp, surface) %>% 
  # filter(medium == "M63") %>% 
  # select(strain, value, temp, surface, medium, replicate) %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = value, color = conditions)) +
  geom_boxplot() +
  facet_wrap(c("medium", "rep_no"), ncol=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_results %>% 
  # select(strain, value, temp, surface, medium, replicate) %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = value, color = conditions)) +
  geom_boxplot() +
  facet_wrap(c("rep_no", "medium"), ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_results %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = value, color = rep_no)) +
  geom_quasirandom() +
  facet_wrap(~ medium) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_results %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  group_by(strain, medium, conditions, rep_no) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = strain, y = value, color = conditions, fill = rep_no)) +
  geom_point(size = 3) +
  facet_wrap(~ medium) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_results %>% 
  # filter(medium == "LB10") %>%
  mutate(conditions = paste(temp, surface)) %>% 
  group_by(strain, medium, conditions, rep_no) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = strain, y = value, color = conditions)) +
  geom_point(size = 3) +
  facet_wrap(c("medium", "rep_no"), ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


all_results %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  # filter(temp == "37") %>%
  ggplot(aes(x = strain, y = rep_no, fill = strength)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(c("medium", "conditions"), ncol = 4) +
  theme_bw() +
  scale_fill_discrete("Conditions") +
  scale_x_discrete("Strain") +
  # scale_y_discrete("Biofilm forming strength", limits=c("absence","weak", "moderate", "strong")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


all_results %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  filter(medium == "BHI") %>%
  ggplot(aes(x = strain, y = strength, fill = conditions)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(c("rep_no", "conditions"), ncol = 4) +
  theme_bw() +
  scale_fill_discrete("Conditions") +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength", limits=c("absence","weak", "moderate", "strong")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


all_results %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  filter(medium == "TSB") %>%
  ggplot(aes(x = strain, y = strength, fill = conditions)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(c("rep_no", "conditions"), ncol = 4) +
  theme_bw() +
  scale_fill_discrete("Conditions") +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength", limits=c("absence","weak", "moderate", "strong")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##### Koniec wykresów

  
  # number of replicates
  all_results %>%
    mutate(conditions = paste0(temp, "_", surface, "_", rep_no)) %>% 
    group_by(strain, medium, strength, conditions) %>%
    summarise(n = length(strength)) %>% 
    arrange(desc(n))
  
  frac_biofilm <- all_results %>%
    mutate(conditions = paste0(temp, "_", surface)) %>% 
    group_by(conditions, strain, medium, strength) %>% 
    summarise(n = length(strength)) %>% 
    mutate(frac = n/sum(n)) %>% 
    select(-n)
  
  # multiple types of biofilm formed
  filter(frac_biofilm, frac != 1)
  
  # only one type of biofilm formed
  filter(frac_biofilm, frac == 1)
  
  dcast(frac_biofilm, conditions + strain + medium ~ strength, fill = 0)
  

  
  
  
  
  
  
  filter(frac_biofilm, frac == 1) %>% 
    ggplot(aes(x = strain, y = strength, fill = conditions)) +
    geom_tile(color = "black", position="dodge") +
    # facet_wrap(c("rep_no", "conditions"), ncol = 4) +
    theme_bw() +
    scale_fill_discrete("Conditions") +
    scale_x_discrete("Strain") +
    scale_y_discrete("Biofilm forming strength") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  
  filter(frac_biofilm, frac == 1) %>% 
    ggplot(aes(x = strain, y = strength, fill = conditions)) +
    geom_tile(color = "black", position="dodge") +
    facet_wrap(c("conditions", "medium"), ncol = 4) +
    theme_bw() +
    scale_fill_discrete("Conditions") +
    scale_x_discrete("Strain") +
    scale_y_discrete("Biofilm forming strength") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  
  nr_different <- filter(frac_biofilm, frac == 1) %>% group_by(strain) %>% summarise(n=length(strain))
  
  
 filter(frac_biofilm, frac == 1) %>% group_by(strain) %>% mutate(cond = paste(conditions, medium)) %>% arrange(strain)
 
 filter(frac_biofilm, frac == 1) %>% group_by(strain) %>% arrange(strain) %>% select(strain, medium, conditions, strength)
 
 DT::datatable(filter(frac_biofilm, frac == 1) %>% group_by(strain) %>% arrange(strain) %>% select(strain, medium, conditions, strength))
 
 
 
 
 ## Weak biofilm formers
 
 Put here all strains that in every condition produce very weak biofilm
 
 ## Strong biofilm formers
 
 Put here all strains that in every condition produce very strong biofilm
 
 ## Largest change of biofilm production
 
 Put here strains that has very strong or very weak production of biofilm, depending on conditions.
 
 
 
 # only one type of biofilm formed | only weak
filter(frac_biofilm, frac == 1, strength == "absence") %>% 
   group_by(strain, conditions) %>% 
   summarise(n = length(medium)) %>% 
   filter(n == 4)

 
filter(frac_biofilm, frac == 1, strength == "strong") %>% 
   group_by(strain, conditions) %>% 
   summarise(n = length(medium)) %>% 
   filter(n == 4) %>% 
   summarise(n_c = length(conditions)) %>% 
   filter(n_c == 4)
 
  
 
aa <- filter(frac_biofilm, frac == 1, strength == "weak" | strength == "strong") 
  group_by(strain, conditions, medium)
 
  
  
  
  filter(frac_biofilm, frac == 1, strength == "medium") %>% 
    group_by(strain, conditions) %>% 
    summarise(n = length(medium)) %>% 
    filter(n == 4)
  
  
  
  
  filter(frac_biofilm, frac == 1) %>% 
    group_by(strain) %>% 
    summarise(n = length(strain)) %>% 
    filter(n == 4)
  
  
  %>% 
    group_by(strain, conditions) %>% 
    summarise(n = length(medium)) %>% 
    filter(n == 4)
  
  
  
  
  
  
  
  