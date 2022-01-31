library(gt)
library(tidyverse)

pigments <- c("19' hexanoyloxyfucoxanthin", "19' butanoyloxyfucoxanthin", "Alloxanthin", "Zeaxanthin", "Chlorophyll b+Divinyl-chlorophyll b", "Fucoxanthin", "Peridinin")
abbreviations <- c("19'-HF", "19'-BF", "Allo", "Zea", "Tchlb", "Fuco", "Peri")
taxonomic <- c("Prymnesiophytes", "Pelagophytes", "Cryptophytes", "Cyanobacteria", "Green Flagellates and Prochlorophytes", "Diatoms", "Dinoflagellates")

pig_info <- data.frame(pigments, abbreviations, taxonomic)


tab <- pig_info %>% 
  rename("taxonomic significance" = taxonomic) %>% 
  gt()

plot_table <- tab %>% 
  tab_options(table.font.size = px(25))

plot_table

gtsave(plot_table, "Output/paper_fig/table_pig.png")
