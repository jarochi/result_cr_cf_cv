library(dplyr)
library(readODS)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)

# mydata <- read_ods("Example_CV.ods", col_names = FALSE)
mydata <- read_ods("measure1.ods", col_names = FALSE)       # Asia
mydata <- read_ods("CV_RESULTS_dk.ods", col_names = FALSE)  # Dominika
mydata <- read.csv("dk_0912.csv", header = FALSE, stringsAsFactors = FALSE)  # Dominika
mydata <- read.csv("js.csv", header = FALSE, stringsAsFactors = FALSE)  # Asia
# plate_scheme <- read.csv("plate_scheme.csv")

start_table <- seq(1, nrow(mydata), 9)
end_table <- start_table + 8

mediums <- c("LB10", "BHI", "TSB", "M63")

# res <- lapply(split(mydata, ceiling(1L:nrow(mydata)/9)), function(ith_table) {
#   raw_data <- as.matrix(ith_table[-1, ])
# 
#   class(raw_data) <- "numeric"
#   colnames(raw_data) <- ith_table[1, ]
#   rownames(raw_data) <- NULL
#   
#   conditions <- strsplit(ith_table[1,1], split = " ") %>% unlist()
#   # browser()
#   # 
#   # raw_data %>% matrix(nrow = 8)
#   # 
#   # raw_data <- raw_data %>%  mutate(temp = conditions[1],
#   #        surface = conditions[2])
# 
#   not_NA_rows <- !apply(raw_data, 1, function(i) mean(is.na(i))) > 0.7
#   not_NA_cols <- !apply(raw_data, 2, function(i) mean(is.na(i))) > 0.7
#   if(sum(dim(raw_data[not_NA_rows, not_NA_cols])) == 0) {
#     data.frame()
#   } else {
#     cbind(melt(raw_data[not_NA_rows, not_NA_cols]), medium = unlist(lapply(mediums, rep, 3))) %>%
#       colnames(seq(1, nrow(raw_data))) %>% 
#       inner_join(mplate_scheme, by = c("row" = "row", "col" = "col"))
#       select(-Var1) %>%
#       rename(strain = Var2)
#     browser()
#   }
# }) %>%  bind_rows()



# empty list to override bind_rows
datalist = list()

for (i in seq(1, ceiling(nrow(mydata)/9))) {
  # read start and end of single table
  x <- start_table[i]
  y <- end_table[i]
  raw_data <- mydata[x:y,]
  # raw_data <- mydata[1:9,] #TEST
  # transform single table
  ODc <- mean(as.numeric(raw_data[7:8, 12]))        # neg control for biofilm strength assumption       <!!!!!!!!!!!!!!!!!!!!!!!!
  # ODc <- 0.125                          # temp for Asia                                       <!!!!!!!!!!!!!!!!!!!!!!!!
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
      cut(x, breaks = c(0, ODc, 2*ODc, 4*ODc, 10), labels = c("no biofilm", "weak", "moderate", "strong"))
    })))
  
  datalist[[i]] <- inner_join(all_res, OD, by = c("strain", "medium", "temp", "surface"))
} 

all_results <- do.call(rbind, datalist)

# Asia
all_results <- all_results %>% 
  mutate(rep_no = rep(1:3, each = 240, times = (nrow(all_results)/720)))

#DK
all_results[["plate_no"]][1:2964]

# JS
all_results %>% 
  group_by(strain, medium, rep_no, temp, surface) %>% 
  filter(medium == "M63") %>% 
  # select(strain, value, temp, surface, medium, replicate) %>% 
  mutate(conditions = paste(temp, surface)) %>% 
ggplot(aes(x = strain, y = value, color = conditions)) +
  geom_boxplot() +
  facet_wrap(c("medium", "rep_no")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_results %>% 
  # select(strain, value, temp, surface, medium, replicate) %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = value, color = conditions)) +
  geom_boxplot() +
  facet_wrap(c("medium", "rep_no"), ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# JS
all_results %>% 
  mutate(conditions = paste(temp, surface)) %>% 
ggplot(aes(x = strain, y = value, color = rep_no)) +
  geom_quasirandom() +
  facet_wrap(~ medium) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# JS
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
  filter(medium == "LB10") %>%
  mutate(conditions = paste(temp, surface)) %>% 
  group_by(strain, medium, conditions, rep_no) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = strain, y = value, color = conditions)) +
  geom_point(size = 3) +
  facet_wrap(c("medium", "rep_no"), ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# DK

# all_results <- all_results[order(all_results$strain),] %>% 
#   mutate(ID = seq(1:nrow(all_results))) 
# 
# all_results <- all_results %>% 
#   mutate(group = as.numeric(cut(all_results$ID, breaks = c(1, 732, 1452, 2184, nrow(all_results)))))

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

#tylko 2 powtórzenie
all_results[2965:nrow(all_results),] %>%
  ggplot(aes(x = strain, y = value)) +
  geom_boxplot() +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# DK
all_results %>% 
  ggplot(aes(x = strain, y = value, color = rep_no)) +
  geom_quasirandom() +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))


all_results <- all_results %>% 
  arrange(strain) 

  all_results[1:2976,] %>% 
  ggplot(aes(x = strain, y = value, color = rep_no)) +
  geom_quasirandom() +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  all_results[2977:nrow(all_results),] %>% 
    ggplot(aes(x = strain, y = value, color = rep_no)) +
    geom_quasirandom() +
    facet_wrap(~ medium, ncol = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

# DK

all_results %>%
  group_by(strain, medium, rep_no) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = strain, y = value, color = rep_no)) +
  geom_point(size = 3) +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

all_results[1:2976,] %>%
  group_by(strain, medium, rep_no) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = strain, y = value, color = rep_no)) +
  geom_point(size = 3) +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_results[2977:nrow(all_results),] %>%
  group_by(strain, medium, rep_no) %>%
  summarise(value = median(value)) %>%
  ggplot(aes(x = strain, y = value, color = rep_no)) +
  geom_point(size = 3) +
  facet_wrap(~ medium, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# JS
all_results %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  filter(medium == "M63") %>%
ggplot(aes(x = strain, y = strength, fill = conditions)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(c("medium", "rep_no"), nrow = 3) +
  theme_bw() +
  scale_fill_discrete("Conditions") +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# DK
all_results %>% 
  na.omit() %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = strength, fill = strength)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(c("medium", "rep_no"), ncol = 1) +
  theme_bw() +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength") +
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
  scale_y_discrete("Biofilm forming strength") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none")

gogo[745:nrow(gogo),] %>% 
  mutate(conditions = paste(temp, surface)) %>% 
  ggplot(aes(x = strain, y = strength, fill = strength)) +
  geom_tile(color = "black", position="dodge") +
  facet_wrap(c("medium", "rep_no"), ncol = 1) +
  theme_bw() +
  scale_x_discrete("Strain") +
  scale_y_discrete("Biofilm forming strength") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none")




gogo <- all_results %>% filter(medium == "M63") %>% arrange(strain)

