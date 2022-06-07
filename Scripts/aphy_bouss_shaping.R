library(tidyverse)

aphy <- read_table2("Data/Boussole/Aphy_boussole.txt", 
                             col_names = FALSE, skip = 1)

aphy_long <- pivot_longer(aphy, c(11:1030), names_to = "col", values_to = "value")

lambda_row <- grep("^[0-9]{3}$", aphy_long$value)
a_row <- lambda_row + 1
a <- aphy_long$value[a_row]
lambda <- aphy_long$value[lambda_row]

abs_df <- data.frame("lamda" = lambda, "aphy" = a)

id_df <- dplyr::select(aphy_long, X1:X10)
id_df <- id_df[lambda_row,]

aphy_bouss <- bind_cols(id_df, abs_df)

names(aphy_bouss) <- c("file", "cruise", "day", "month", "year", "time", "ctd", "bottle", "depth", "tchla", "lambda", "aphy")

aphy_bouss <- mutate(aphy_bouss, sample = paste(ctd, depth))
ggplot(filter(aphy_bouss, cruise == "033"))+
  geom_path(aes(x = lambda, y = aphy, group = sample, colour = depth))

write_csv(aphy_bouss, "Data/Boussole/aphy_bouss.csv")
