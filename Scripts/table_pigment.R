pigments <- c("19' hexanoyloxyfucoxanthin", "19' butanoyloxyfucoxanthin", "Alloxanthin", "Zeaxanthin", "Chlorophyll b+Divinyl-chlorophyll b")
abbreviations <- c("19'-HF", "19'-BF", "Allo", "Zea", "Tchlb")
taxonomic <- c("prymnesiophytes", "pelagophytes", "cryptophytes", "cyanobacteria", "green flagellates and prochlorophytes")

pig_info <- data.frame(pigments, abbreviations, taxonomic)


tab <- pig_info %>% 
  gt()

plot_table <- tab %>% 
  tab_options(table.font.size = px(25)) %>% 
  tab_header(title = "Taxonomic pigments used in this study",
             subtitle = "Pigment name, abbreviation and taxonomic sgnificance")
  gtsave("Output/paper_fig/table.png")

plot_table
