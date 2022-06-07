library(readxl)
library(tidyverse)

soclim <- read_delim("Data/soclim/soclim_btl_V9.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
soclim <- janitor::clean_names(soclim[, c(1:103)])
#open pigment absorption in solution data
spectre <- read_excel("Data/Spectres_annick.xlsx") %>% janitor::clean_names()

#create a vector with all pigments names
pignames <- c("chlorophyll_c1_c2_mg_dvp",
              "peridinin",
              "x19_butanoyloxyfucoxanthin",
              "fucoxanthin",
              "x19_hexanoyloxyfucoxanthin",
              "diadinoxanthin",
              "alloxanthin",
              "zeaxanthin",
              "divinyl_chlorophyll_b",
              "chlorophyll_b",
              "divinyl_chlorophyll_a",
              "chlorophyll_a")


soclim_pig <- select(soclim, pignames)
soclim_pig[soclim_pig == "LOD"] <- "0"
soclim_pig <- mutate_all(soclim_pig, function(x) as.numeric(x))

#select pigments
soclim_pig <- select(soclim_pig, chlc12 = "chlorophyll_c1_c2_mg_dvp",
                     peri = "peridinin",
                     x19bf = "x19_butanoyloxyfucoxanthin",
                     fuco = "fucoxanthin",
                     x19_bf = "x19_hexanoyloxyfucoxanthin",
                     diad ="diadinoxanthin",
                     allox = "alloxanthin",
                     zea = "zeaxanthin",
                     dv_chlb = "divinyl_chlorophyll_b",
                     chl_b = "chlorophyll_b",
                     dv_chla = "divinyl_chlorophyll_a",
                     chl_a = "chlorophyll_a")


#select only pigments that contribute to photosynthesis
soclim_pur <- select(soclim_pig, - zea, -diad)

#change df to matrix 
soclim_mat <- as.matrix(soclim_pig)
soclim_mat <- t(soclim_mat)

#idem for photosynthetical
soclim_pur_mat <- as.matrix(soclim_pur)
soclim_pur_mat <- t(soclim_pur_mat)

#clean spectre and create a df with photosynthetical spectra
spectre[is.na(spectre)] <- 0
spectre <- select(spectre, - ss_car, -a_car)
spectre_pur <- select(spectre, -zea, -diad)

#create matrix
spectre_mat <- as.matrix(select(spectre, - lambda))
spectre_pur_mat <- as.matrix(select(spectre_pur, -lambda))

#matrix multiplication to get a* * [pig1] + a* * [pig2] + ...
result <- spectre_mat %*%  soclim_mat

result <- t(result)

#idem for pur
result_pur <- spectre_pur_mat %*% soclim_pur_mat
result_pur <- t(result_pur)

#create two df with results
result_df <- as.data.frame(result)
pur_df <- as.data.frame(result_pur)

#name cols
colnames(result_df) = paste('a', spectre$lambda, sep = '')
colnames(pur_df) = paste('pur', spectre$lambda, sep = '')

#combine the two results df
soclim_tot <- bind_cols(soclim, result_df, pur_df)

#create a long one with wl
soclim_long <- soclim_tot %>% mutate(round_depth = round(depth_m), ratio = a440/a470) %>%
  pivot_longer(c(a400:pur700), names_to = "wavelength", values_to = "abs")

#separate the type of abs and the num wl by a _
soclim_long$wavelength <- gsub('^([a-z]+)([0-9]{3})$', '\\1_\\2', soclim_long$wavelength)

#split the column to get the wavelength as num and the type of abs in two different columns
soclim_long <- soclim_long %>% separate(wavelength, into = c('type', 'lambda'), '_')

#pivot wider
soclim_long <- soclim_long %>% pivot_wider(names_from = type, values_from = abs)

#compute the factor of correction and the corrected abs
soclim_long <- soclim_long %>% mutate(factor = pur/ a)

#create 4 different df for each type of abs
pur_df <- select(soclim_long, depth_m, lambda, pur)
a_df <- select(soclim_long, depth_m, lambda, a)

#create a fubction to get df with only absorption values, in wide, each column for each wl value
transform_wide <- function(data, value, name){
  df <- data %>% group_by(lambda) %>% 
    mutate(row = row_number()) %>% 
    pivot_wider(names_from = lambda, values_from = {{value}}) %>% 
    select(-row, - depth_m)
  colnames(df) <- gsub(' ', '',paste(name, colnames(df), sep = ''))
  return(df)
}

#apply it on each type of abs
pur_df <- transform_wide(pur_df, pur, 'pur')
a_df <- transform_wide(a_df, a, 'a')


#add them to lov dataset. Nwo you get 3 new abs to the dataset !
soclim_tot <- bind_cols(soclim, a_df, pur_df)

write_csv(soclim_tot, "Data/soclim_abs.csv")
