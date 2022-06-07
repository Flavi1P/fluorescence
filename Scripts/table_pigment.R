library(gt)
library(tidyverse)

pigments <- c("Zeaxanthin", "Chlorophyll b+Divinyl-chlorophyll b", "19' hexanoyloxyfucoxanthin", "19' butanoyloxyfucoxanthin", "Alloxanthin", "Fucoxanthin", "Peridinin")
abbreviations <- c( "Zea", "Tchlb","19'-HF", "19'-BF", "Allo", "Fuco", "Peri")
taxonomic <- c("Cyanobacteria", "Green Flagellates and Prochlorophytes", "Prymnesiophytes", "Pelagophytes", "Cryptophytes",  "Diatoms", "Dinoflagellates")
size <- c( "Pico", "Pico", "Nano", "Nano", "Nano", "Micro", "Micro")

pig_info <- data.frame(pigments, abbreviations, taxonomic, size)


tab <- pig_info %>% 
  rename("Pigments" = pigments, "Abbreviations" = abbreviations, "Taxonomic significance" = taxonomic, "Size class" = size) %>% 
  gt()

plot_table <- tab %>% 
  tab_options(table.font.size = px(25))

plot_table

gtsave(plot_table, "Output/paper_fig/table_pig.html")

name_dataset <- c("Glo-Argo", "Glo-aphy", "Med-Bouss")
param <- c("HPLC, Fluo", "HPLC, Absorption", "HPLC, Absorption, Fluo")
num_samples <- c(335, 3340, 843)
scale <- c("Global", "Global", "Regional")

data_info <- data.frame(name_dataset, scale, param, num_samples)

tab_data <- data_info %>% 
  rename("Name of the dataset" = name_dataset,
         "Spatial scale" = scale,
         "Available measurments" = param,
         "Number of samples" = num_samples) %>% 
  gt()

plot_table_data <- tab_data %>% 
  tab_options(table.font.size = px(25))

plot_table_data

gtsave(plot_table_data, "Output/paper_fig/table_dataset.png")

fitted_slope <- argo_new %>% group_by(code) %>% do(model = lm(fluo_biased ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "t_chla")

r_squarred <- argo_new %>% group_by(code) %>% do(model = lm(fluo_biased ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, glance)) %>% 
  unnest(coef) %>% 
  pull("adj.r.squared")

fitted_slope$r_sqr <- r_squarred

tab_slope <- fitted_slope %>% 
  select(-term, - statistic) %>% 
  mutate(estimate = round(estimate, 2),
         std.error = round(std.error, 3),
         p.value = "< 0.001",
         r_sqr = round(r_sqr, 2)) %>% 
  rename("Code" = code,
         "Slope Factor" = estimate,
         "Std Error" = std.error,
         "P.value" = p.value,
         "Adjusted R²" = r_sqr) %>% 
  gt() %>% 
  tab_header(
    title = "Slope Factor"
  )

tab_slope
gtsave(tab_slope, "Output/paper_fig/tab_slope.png")

fitted_yield <- argo_new %>% group_by(code) %>% do(model = lm(counts ~ a_spe_470 + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "a_spe_470")

r_squarred_yield <- argo_new %>% group_by(code) %>% do(model = lm(counts ~ a_spe_470 + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, glance)) %>% 
  unnest(coef) %>% 
  pull("adj.r.squared")

fitted_yield$r_sqr <- r_squarred_yield

tab_yield <- fitted_yield %>% 
  select(-term, - statistic) %>% 
  mutate(estimate = round(estimate),
         std.error = round(std.error, 3),
         p.value = "< 0.001",
         r_sqr = round(r_sqr, 2)) %>% 
  rename("Code" = code,
         "Quantum Yield" = estimate,
         "Std Error" = std.error,
         "P.value" = p.value,
         "Adjusted R²" = r_sqr) %>% 
  gt() %>% 
  tab_header(
    title = "Quantum yield of fluorescence"
  )

tab_yield
gtsave(tab_yield, "Output/paper_fig/tab_yield.png")

fitted_abs <- argo_new %>% group_by(code) %>% do(model = lm(a_470 ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "t_chla")

r_squarred_abs <- argo_new %>% group_by(code) %>% do(model = lm(a_470 ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, glance)) %>% 
  unnest(coef) %>% 
  pull("adj.r.squared")

fitted_abs$r_sqr <- r_squarred_abs

tab_abs <- fitted_abs %>% 
  select(-term, - statistic) %>% 
  mutate(estimate = round(estimate, 2),
         std.error = round(std.error, 4),
         p.value = "< 0.001",
         r_sqr = round(r_sqr, 2)) %>% 
  rename("Code" = code,
         "a*470" = estimate,
         "Std Error" = std.error,
         "P.value" = p.value,
         "Adjusted R²" = r_sqr) %>% 
  gt() %>% 
  tab_header(
    title = "Chla-specific absorption"
  )

tab_abs
gtsave(tab_abs, "Output/paper_fig/tab_abs.png")

fit_abs <- fitted_abs %>% 
  select(-term, - statistic) %>% 
  mutate(estimate = round(estimate, 2),
         std.error = round(std.error, 4),
         p.value = "<0.001",
         r_sqr = round(r_sqr, 2)) 

fit_slope <- fitted_slope %>%  select(-term, - statistic) %>% 
  mutate(estimate = round(estimate, 2),
         std.error = round(std.error, 3),
         p.value = "<0.001",
         r_sqr = round(r_sqr, 2))

fit_yield <- fitted_yield %>% 
  select(-term, - statistic) %>% 
  mutate(estimate = round(estimate),
         std.error = round(std.error, 3),
         p.value = "<0.001",
         r_sqr = round(r_sqr, 2)) 

full_table <- bind_cols(fit_slope, fit_abs) %>% bind_cols(fit_yield)

tab_full <-  full_table %>%
  select(-code...6, - code...11) %>% 
  janitor::clean_names() %>% 
  rename("Code" = code_1,
         "Slope Factor_ Slope Factor" = estimate_2,
         "Slope Factor_Std Error" = std_error_3,
         "Slope Factor_P.value" = p_value_4,
         "Slope Factor_Adjusted R²" = r_sqr_5,
         "a*470_a*470" = estimate_7,
         "a*470_Std Error" = std_error_8,
         "a*470_P.value" = p_value_9,
         "a*470_Adjusted R²" = r_sqr_10,
         "Quantum Yield_Quantum Yield" = estimate_12,
         "Quantum Yield_Std Error" = std_error_13,
         "Quantum Yield_P.value" = p_value_14,
         "Quantum Yield_Adjusted R²" = r_sqr_15) %>% 
  mutate(Province = case_when(
    Code == "ANTA" ~ "Antarctic Province",
    Code == "ARCH" ~ "Archipelagic Deep Bassins Province",
    Code == "ARCT" ~ "Atlantic Arctic Province",
    Code == "BPLR" ~ "Boreal Polar Province",
    Code == "EMED"~ "Eastern Mediterranean Sea",
    Code == "SANT" ~ "Subantarctic Province",
    Code == "SPSG" ~ "S. Pacifc Subtropical Gyre Province",
    Code == "WMED" ~ "Western Mediterranean Sea"
  )) %>% 
  select(Province, Code, everything()) %>% 
  gt() %>% 
  tab_spanner_delim(delim = "_")

  
  

tab_full
gtsave(tab_full, "Output/paper_fig/tab_full.png")

