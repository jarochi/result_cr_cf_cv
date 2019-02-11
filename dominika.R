library(dplyr)
library(readODS)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)

mydata <- read.csv("dk_0912.csv", header = FALSE, stringsAsFactors = FALSE)  # Dominika

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
  ODc <- mean(as.numeric(raw_data[7:8, 12]))        # neg control for biofilm strength assumption  
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
           plate_no = rep(i, times = nrow(mplate_scheme)),             # if we have strains on more than 1 plate, do przerobienia na zwykłe repeat zamiast numeru płytki
           rep_no = ifelse(i<51, "1", "2"))  %>%                       # DK liczpa powtórzeń
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

# sort by strain
all_results <- all_results[order(all_results$strain),] %>%
  mutate(ID = seq(1:nrow(all_results)))

# podzielenie na grupy do łatwiejszego wyświetlania na plotach
all_results <- all_results %>%
  mutate(group = as.numeric(cut(all_results$ID, breaks = c(1, 732, 1452, 2184, nrow(all_results)))))

all_results %>%
  ggplot(aes(x = strain, y = value)) +
  geom_boxplot() +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

all_results %>%
  ggplot(aes(x = strain, y = value, fill = rep_no)) +
  geom_boxplot() +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

all_results %>% 
  ggplot(aes(x = strain, y = value, color = rep_no)) +
  geom_quasirandom() +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

# all_results[1:2976,] %>% 
#   ggplot(aes(x = strain, y = value, color = rep_no)) +
#   geom_quasirandom() +
#   facet_wrap(~ medium, ncol = 1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# all_results[2977:nrow(all_results),] %>% 
#   ggplot(aes(x = strain, y = value, color = rep_no)) +
#   geom_quasirandom() +
#   facet_wrap(~ medium, ncol = 1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_results %>%
  group_by(strain, medium, rep_no) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = strain, y = value, color = rep_no)) +
  geom_point(size = 3) +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

# all_results[1:2976,] %>%
#   group_by(strain, medium, rep_no) %>%
#   summarise(value = median(value)) %>%
#   ggplot(aes(x = strain, y = value, color = rep_no)) +
#   geom_point(size = 3) +
#   facet_wrap(~ medium, ncol = 1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# all_results[2977:nrow(all_results),] %>%
#   group_by(strain, medium, rep_no) %>%
#   summarise(value = median(value)) %>%
#   ggplot(aes(x = strain, y = value, color = rep_no)) +
#   geom_point(size = 3) +
#   facet_wrap(~ medium, ncol = 1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


all_results %>% 
  # na.omit() %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = strength, fill = strength)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(c("medium", "rep_no"), ncol = 1) +
  theme_bw() +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength", limits=c("absence","weak", "moderate", "strong")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position="none")

all_results[1:2976,] %>% 
  na.omit() %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = strength, fill = strength)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(c("medium", "rep_no"), ncol = 1) +
  theme_bw() +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength", limits=c("absence","weak", "moderate", "strong")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none")



all_results %>% 
  na.omit() %>% 
  filter(medium == "BHI") %>%
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = rep_no, fill = strength)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(c("medium", "rep_no"), ncol = 1) +
  theme_bw() +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength", limits=c("absence","weak", "moderate", "strong")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position="none")



all_results %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = rep_no, fill = strength)) +
  geom_tile(color = "black") +
  facet_wrap(c("medium"), ncol = 1) +
  theme_bw() +
  scale_x_discrete("Strain") +
  # scale_y_discrete("Biofilm forming strength", limits=c("absence","weak", "moderate", "strong")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))


all_results[1:2976,] %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = rep_no, fill = strength)) +
  geom_tile(color = "black") +
  facet_wrap(c("medium"), ncol = 1) +
  theme_bw() +
  scale_x_discrete("Strain") +
  # scale_y_discrete("Biofilm forming strength", limits=c("absence","weak", "moderate", "strong")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

all_results[2977:nrow(all_results),] %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = rep_no, fill = strength)) +
  geom_tile(color = "black") +
  facet_wrap(c("medium"), ncol = 1) +
  theme_bw() +
  scale_x_discrete("Strain") +
  # scale_y_discrete("Biofilm forming strength", limits=c("absence","weak", "moderate", "strong")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))









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
  facet_wrap(c("conditions", "medium"), ncol = 1) +
  theme_bw() +
  scale_fill_discrete("Conditions") +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



nr_different <- filter(frac_biofilm, frac == 1) %>% group_by(strain) %>% summarise(n=length(strain))


filter(frac_biofilm, frac == 1) %>% group_by(strain) %>% mutate(cond = paste(conditions, medium)) %>% arrange(strain)

filter(frac_biofilm, frac == 1) %>% group_by(strain) %>% arrange(strain) %>% select(strain, medium, conditions, strength)

DT::datatable(filter(frac_biofilm, frac == 1) %>% group_by(strain) %>% arrange(strain) %>% select(strain, medium, conditions, strength))




# only one type of biofilm formed | only weak
filter(frac_biofilm, frac == 1, strength == "weak") %>% 
  group_by(strain, conditions) %>% 
  summarise(n = length(medium)) %>% 
  filter(n == 4)


filter(frac_biofilm, frac == 1, strength == "strong") %>% 
  group_by(strain, conditions) %>% 
  summarise(n = length(medium)) %>% 
  filter(n == 4) %>% 
  select(strain)

