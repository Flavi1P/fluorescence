library(tidyverse)
library(sf)
library(s2)
library(readxl)
library(ggrepel)
library(gridExtra)
library(grid)
library(janitor)
library(Metrics)
library(cowplot)
library(broom)
library(lmodel2)
library(patchwork)
library(RColorBrewer)

#open datas
path <- "Data/Longhurst"
argo <- read_csv("Data/Argo Floats output/argo_matchup_ze.csv")
ref <- read_csv("Data/Argo Floats output/ref.csv")

#correct longitude for takap floats
argo$lon[argo$lovbio %in% c("takapm014b", "takapm013b", "takapm009b", "takapm005b")] <- -argo$lon[argo$lovbio %in% c("takapm014b", "takapm013b", "takapm009b", "takapm005b")]

#remove bad match up on lovbio024c, lovibo026c and lovbio027b
argo <- filter(argo, lovbio != "lovbio024c" & lovbio != "lovbio026c" & lovbio != "lovbio027b" & lovbio != "lovbio040b")

#create bioregion dataframe
longhurst_sf <- read_sf(dsn = path.expand(path), quiet = TRUE)
names(longhurst_sf) <- c("code", "region", "geometry")

#define the bioregion associate to each argo sample
pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(argo),
                                     function(i) {st_point(as.numeric(argo[i,c("lon", "lat") ]))}), list("crs" = 4326))) 

pnts_trans <- st_transform(pnts_sf, 4326)
longhurst_trans <- st_transform(longhurst_sf, 4326)  

sf_use_s2(FALSE) #correction of a bug that appeared with a sf package update
argo$code <- apply(st_intersects(longhurst_trans, pnts_trans, sparse = FALSE), 2, 
                   function(col) { 
                     longhurst_trans[which(col), ]$code
                   })

argo$code <- as.character(argo$code)
argo$code[argo$code=="character(0)"] <- "BPLR"
argo <- argo %>% mutate(code = case_when(code == "MEDI" & lon > 18 ~ "EMED",
                                         code == "MEDI" & lon <= 18 ~ "WMED",
                                         code != "EMED" & code != "WMED" ~ code)) #split the med sea in two regions

#Create a custom color scale
polarcolor <- brewer.pal(4,"PuBu")
medcolor <- brewer.pal(4, "Greens")[c(3,4)]
equatcolor <- brewer.pal(3, "OrRd")[c(2,3)]

myColors <- c(polarcolor, medcolor, equatcolor)

names(myColors) <- c("SANT", "ARCT", "ANTA", "BPLR", "EMED", "WMED", "SPSG", "ARCH")

#create a df with the fluorescence value and raw ratio
argo_new <- argo %>% mutate(fluo_biased = fluo *2, ratio = fluo_biased/t_chla) %>% 
  filter(lovbio != "lovbio077b")

#select only positiv ratios (no negativ values in the fluo) and filter on the depth z/ze < 1.5
argo_new <- filter(argo_new, ratio > 0) %>% filter(code != "SARC") 
    # mutate(z_ze = depth/ze) %>% 
    # filter(z_ze < 1)

#compute basopriton and counts
argo_new <- argo_new %>% mutate(a_spe_470 = 0.0332 * (t_chla^-0.368),
                                counts = fluo_biased/0.007 + 55,
                                a_470 = a_spe_470 * t_chla)


#create the mean of the ratio by region
region_argo <- argo_new  %>% group_by(code) %>% summarise_at(vars(ratio), c(mean, sd, min, max), na.rm = TRUE) %>% ungroup()
names(region_argo) <- c("code", "mean", "sd", "min", "max")

#make sure sd does not exceed 0
region_argo$sd <- ifelse(region_argo$sd > region_argo$mean, region_argo$mean, region_argo$sd)

#add the code of the region
codref <- read_excel("Data/Longhurst/Longhurst_Province_Summary.xls", 
                     range = "A17:B70", col_names = FALSE)
names(codref) <- c("code", "region")

codref <- codref %>% separate(region, "-", into = c("biome", "region"))

region_argo <- left_join(region_argo, codref)
argo_new <- left_join(argo_new, codref)

#write_csv(argo_new, "Output/Data/argo_matchup.csv")

map_vec <- map_data("world")
g1 <- ggplot(filter(argo_new))+
  geom_point(aes(x = lon, y = lat, fill = code), size = 4, shape = 21)+
  geom_polygon(aes(x = long, y = lat, group = group), data = map_vec, fill = "gray28")+
  xlab("Longitude (°E)")+
  ylab("Latitude (°N)")+
  theme_bw(base_size = 18)+
  scale_fill_manual(name = "code",values = myColors)+
  coord_quickmap()

region_argo$code <- as.character(region_argo$code)
g2 <- ggplot(region_argo)+
  geom_col(aes(x = reorder(code, mean), y = mean, fill = code))+
  geom_errorbar(aes(x = code, ymin = mean - sd, ymax = mean + sd))+
  xlab("Oceanic province")+
  ylab("Fluo/Chla ratio")+
  geom_errorbar(aes(code, ymax = 1, ymin = 1),
                size=1, linetype = "longdash", inherit.aes = F, width = 1)+
  guides(fill = FALSE)+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 20)

g1+g2

g3 <- ggplot(argo_new)+
  geom_boxplot(aes(x = code, y = ratio, fill = code))+
  geom_jitter(aes(x = code, y = ratio, colour = code))+
  xlab("Longitude (°E)")+
  ylab("ratio")+
  scale_fill_brewer(palette = 'Set1')+
  coord_quickmap()
g2 + g3

#compute the ratio by using a linear regrssion
fitted_slope <- argo_new %>% group_by(code) %>% do(model = lm(fluo_biased ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "t_chla")

ggplot(filter(argo_new, code == "SANT"))+
  geom_point(aes(x = t_chla, y = fluo_biased, colour = lovbio))+
  scale_colour_brewer(palette = "Set1")


fitted_slope$code <- factor(fitted_slope$code, levels = c("ANTA", "SANT", "ARCT", "BPLR", "EMED", "WMED", "SPSG", "ARCH"))


gslope <-ggplot(fitted_slope)+
  geom_bar(aes(y = estimate, x = code, fill = code), stat = "identity")+
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error, x = code))+
  scale_fill_manual(name = "Code",values = myColors, guide = "none")+
  theme_bw(base_size = 16)+
  ylab("Slope Factor")+
  xlab("Bioregion")+
  geom_vline(xintercept = c(4.5, 6.5)) +
  annotate(geom = 'text', x = 2.5, y = 3, label = 'Polar') +
  annotate(geom = 'text', x = 5.5, y = 3, label = 'Temperate') +
  annotate(geom = 'text', x = 7.5, y = 3, label = 'Subtropical')


g1/ gslope +
  plot_layout(guide = "collect", height = c(1,2))
#ggsave("Output/Figures/slope_factor.jpg", dpi = 300, width = 20, height = 15, unit = c("cm"))
#ggsave("Output/paper_fig/slope_factor.png", height = 8, width = 10)

sample_size <- argo %>% 
  group_by(code) %>%
  summarise(min_chl = min(t_chla),
            max_chl = max(t_chla)) %>% 
  ungroup() %>% 
  left_join(count(argo, code)) %>% 
  mutate(province = case_when(
    code == "ANTA" ~ "Antarctic Province",
    code == "ARCH" ~ "Archipelagic Deep Bassins Province",
    code == "ARCT" ~ "Atlantic Arctic Province",
    code == "BPLR" ~ "Boreal Polar Province",
    code == "EMED"~ "Eastern Mediterranean Sea",
    code == "SANT" ~ "Subantarctic Province",
    code == "SPSG" ~ "S. Pacifc Subtropical Gyre Province",
    code == "WMED" ~ "Western Mediterranean Sea"
  ))

tab <- fitted_slope %>% 
  tibble() %>% 
  select(code, estimate, std.error, p.value) %>% 
  right_join(sample_size) %>% 
  filter(code != "SARC") %>% 
  select(province, "Code" = code,"Number of samples" = n, "[Chla] minimum" = min_chl, "[Chla] maximum" = max_chl, "slope factor" = estimate, "standard error" = std.error, "p-value" = p.value) %>% 
  gt(rowname_col = "province") %>% 
  fmt_number(
    columns = c("slope factor", "[Chla] minimum", "[Chla] maximum", "standard error", "p-value"),
    decimals = 3,
    use_seps = FALSE
  )

plot_table <- tab %>% 
  tab_options(table.font.size = px(25)) %>% 
  tab_header(title = "Argo Dataset",
             subtitle = "Dataset of concomittant measurement of [Chla] from BGC-Argo fluorometer and HPLC") %>% 
  tab_stubhead(label = "Longhurst Bioregion") %>% 
  gtsave("Output/paper_fig/table.png")

plot_table
#ggsave("Output/paper_fig/barplot.png", height = 8, width = 10)

argo <- argo %>% mutate( argo, fluo_biased = fluo * 2)
argo <- filter(argo, lovbio != "lovbio024c" & lovbio != "lovbio026c" & lovbio != "lovbio027b" & lovbio != "lovbio040b")
argo <- filter(argo, (code == "SPSG" & depth < 100) | (code != "SPSG"))

fitted_slope2 <- argo %>% group_by(code) %>% do(model = lmodel2(fluo_biased ~ t_chla, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "Slope" & method == "OLS")

ggplot(filter(fitted_slope2, code != "SPSG" & code != "SARC"))+
  geom_bar(aes(y = estimate, x = reorder(code, estimate), fill = code), stat = "identity")+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, x = reorder(code, estimate)))+
  scale_fill_brewer(palette = "Set1", guide = FALSE)+
  theme_bw()+
  ylab("Slope Factor")+
  xlab("Oceanic province")+
  ggtitle("Slope factor")

argo %>% 
  gt() %>% 
  tab_header

argo_new$biome[is.na(argo_new$biome)] <- "Westerlies "

fitted_yield <- argo_new %>% group_by(code) %>% do(model = lm(counts ~ a_470 + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "a_470")

fitted_abs <- argo_new %>% group_by(code) %>% do(model = lm(a_470 ~ t_chla + 0, data = .)) %>% ungroup() %>% 
  transmute(code, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "t_chla")

fitted_abs$code <- factor(fitted_abs$code, levels = c("ANTA", "SANT", "ARCT", "BPLR", "EMED", "WMED", "SPSG", "ARCH"))
fitted_yield$code <- factor(fitted_yield$code, levels = c("ANTA", "SANT", "ARCT", "BPLR", "EMED", "WMED", "SPSG", "ARCH"))



yname2 <- expression(atop("a*(470)",~(m^2%.%"(mg"%.%"chla)"^{"-1"})))
yname3 <- expression(atop(phi~("Counts"%.%m^-1)))


gslope <-ggplot(fitted_slope)+
  geom_bar(aes(y = estimate, x = code, fill = code), stat = "identity")+
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error, x = code))+
  scale_fill_manual(name = "code",values = myColors, guide = "none")+
  theme_bw(base_size = 16)+
  ylab("Slope factor")+
  xlab("")+
  geom_vline(xintercept = c(4.5, 6.5)) +
  annotate(geom = 'text', x = 2.5, y = 3, label = 'Polar') +
  annotate(geom = 'text', x = 5.5, y = 3, label = 'Temperate') +
  annotate(geom = 'text', x = 7.5, y = 3, label = 'Subtropical')+
  ylim(0,3.2)

gabs <- ggplot(fitted_abs)+
  geom_bar(aes(y = estimate, x = code, fill = code), stat = "identity")+
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error, x = code))+
  scale_fill_manual(name = "code",values = myColors, guide = "none")+
  theme_bw(base_size = 16)+
  ylab(yname2)+
  xlab("")+
  ylim(0,0.085)+
  geom_vline(xintercept = c(4.5, 6.5))+
  annotate(geom = 'text', x = 2.5, y = 0.08, label = 'Polar') +
  annotate(geom = 'text', x = 5.5, y = 0.08, label = 'Temperate') +
  annotate(geom = 'text', x = 7.5, y = 0.08, label = 'Subtropical')

gyield <- ggplot(fitted_yield)+
  geom_bar(aes(y = estimate, x = code, fill = code), stat = "identity")+
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error, x = code))+
  scale_fill_manual(name = "code",values = myColors, guide = "none")+
  theme_bw(base_size = 16)+
  ylab(yname3)+
  xlab("Bioregion")+
  geom_vline(xintercept = c(4.5, 6.5))+
  annotate(geom = 'text', x = 2.5, y = 14500, label = 'Polar') +
  annotate(geom = 'text', x = 5.5, y = 14500, label = 'Temperate') +
  annotate(geom = 'text', x = 7.5, y = 14500, label = 'Subtropical')+
  ylim(0,15100)

gslope/gabs/gyield

#ggplot2::ggsave("Output/Figures/full_barplot.jpg", width = 20, height = 15, unit = "cm", dpi = "print")
#ggsave("Output/paper_fig/full_barplot.png", height = 13, width = 13)


argo_new <- argo_new %>% mutate(dp_rate = (peri + but + hex + zea + fuco + allo + chlb) / chla)

summary_dp <- argo_new %>% group_by(code) %>% 
  summarise("mean_dp" = mean(dp_rate, na.rm = TRUE),
            "sd_dp" = sd(dp_rate, na.rm = TRUE)) %>% 
  ungroup() %>% 
  right_join(fitted_yield)

summary_dp$mean_dp[summary_dp$code == "BPLR"] <- 0

gyield/
ggplot(summary_dp)+
  geom_bar(aes(y = mean_dp, x = reorder(code, estimate), fill = code), stat = "identity")+
  geom_errorbar(aes(ymin = mean_dp - sd_dp, ymax = mean_dp + sd_dp, x = reorder(code, estimate)))+
  scale_fill_brewer(palette = "Set1", guide = FALSE)+
  ylab("AP/T-Chl A")+
  xlab("Oceanic province")+
  ggtitle("Ratio between accessories pigments and chlorophyll a")+
  theme_bw()

#ggsave("Output/paper_fig/yiled_ap.png", height = 7, width = 10)

ggplot(argo_new)+
  geom_boxplot(aes(x = code, y = dp_rate, fill = code))+
  geom_jitter(aes(x = code, y = dp_rate))+
  scale_fill_brewer(palette = "Set1")+
  theme_bw()


argo_new <- argo_new %>% mutate("yield" = counts/a_470)


        #write_csv(argo_new, "Data/argo_new.csv")

argo_neww <- argo_new %>% mutate(wdp = 1.56 * fuco + 0.92 * peri + 4.8 * allo + 1.02 * but + 1.12 * hex + 1.51 * zea + 0.69 * chlb,
                                 micro = (1.56 * fuco + 0.92 * peri)/wdp,
                                 nano = (4.8 * allo + 1.02 * but + 1.51 * hex)/wdp,
                                 pico = (1.51 * zea + 0.69 * chlb)/wdp)

ggplot(argo_neww)+
  coord_tern()+
  stat_interpolate_tern(geom = "polygon", formula = value~x+y,
                        method = lm, n = 50,
                        breaks = seq(0,7, by = 1),
                        aes(x = micro, z = nano, y = pico, value = ratio, fill =..level..), expand = 4)+
  geom_point(aes(x= micro, y = nano, z = pico, colour = code, size = ratio))+
  scale_color_brewer(palette = "Set1", name = "Oceanic region")+
  theme_bw(base_size = 20)+
  guides(colour = FALSE)+
  weight_percent()+
  scale_fill_gradient(name = "Fluo/Chla ratio", low = "lightgrey", high = "gray45")
