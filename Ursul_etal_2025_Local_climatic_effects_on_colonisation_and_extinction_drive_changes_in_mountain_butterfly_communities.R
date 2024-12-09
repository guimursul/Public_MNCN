
## Title: Local climatic effects on colonisation and extinction drive changes in mountain butterfly communities
## Authors: Guim Ursul, Mario Mingarro, Sara Castro-Cobo, Juan Pablo Cancela, Helena Romo and Robert J. Wilson
## Affiliation fist author: Museo Nacional de Ciencias Naturales (MNCN-CSIC), Madrid
## Corresponding author: Guim Ursul, guim.ursul@mncn.csic.es

## install packages if not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,readxl,writexl,ggplot2,tidyverse,MuMIn,
               sf,sp,spdep,ggpubr,lme4,emmeans,AICcmodavg,broom)

#### GET DIRECTORY AND CREATE FOLDERS TO SAVE OUTPUTS ####
directory <- rstudioapi::getActiveDocumentContext()$path
directory <- dirname(directory)
dir.create(paste0(directory,"/Figures")) # to save figures from script
dir.create(paste0(directory,"/Tables")) # to save tables from script

#### DATA PROCESSING ####
##### 1. Microclimate modelling example (NO NEED TO RUN) ####
# data obtained from microclima model, as shown below, is included in sheet "3. Main_analysis_data" in supp. data
devtools::install_github('mrke/NicheMapR') # v3.1.0.
devtools::install_github("ilyamaclean/microclima") # v0.1.0.
library(raster)
library(NicheMapR)
library(microclima)
library(foreach)

coords_list <-read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "1. coordinates_list") %>% 
  dplyr::select(Code, lon, lat)

coords_sf <- st_as_sf(x = coords_list, 
                        coords = c("lon", "lat"),
                        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

## create table to save output
tmax_table <- data.frame(matrix(ncol=5, nrow=0))

colnames(tmax_table) <- c("code","month","year","survey","tmax")

# From 1980 to 1989
fi <- data.frame(temp_date = seq(as.Date("1980-01-01"), length=120, by="month"))  
ff <- data.frame(temp_date = seq(as.Date("1980-02-01"), length=120, by="month"))-1  

f_start <- fi %>% 
  separate(temp_date, into = c("day", "month", "year")) %>%
  mutate(good_date = paste(year, month, day, sep = "/")) %>%
  dplyr::select(good_date)

f_end <- ff %>% 
  separate(temp_date, into = c("day", "month", "year")) %>%
  mutate(good_date = paste(year, month, day, sep = "/")) %>%
  dplyr::select(good_date)

dates <- f_end %>% 
  separate(good_date, into = c("day", "month", "year")) # separate dates to save them in each iteration

foreach (j=1:length(coords_list)) %dopar% {
  lat <- coords_list$lat[j]
  long <- coords_list$lon[j]
  mdt <- microclima::get_dem(lat = lat, long = long, resolution = 30) ## microclimate modelled across 5x5 km square around coord  
  for (i in 1:nrow(f_end))  {
    temp <- runauto(mdt,
                    f_start[i,], f_end[i,], 
                    hgt = 0.1, 
                    l = NA, x = NA,
                    coastal = FALSE,
                    habitat = 7, 
                    plot.progress = FALSE, save.memory = FALSE)
    tmax <- temp$tmax
    a <- coords_sf[j,]
    tmax_ext <- raster::extract(tmax, a, fun=mean) # extract tmax value from coordinates

    ## save tmax output to table
    temp_table <- data.frame(matrix(ncol = 5, nrow = 1))
    colnames(temp_table) <- colnames(tmax_table)
    temp_table[1,] <- c(coords_list$Code[j],dates$month[i],dates$year[i],"historical",tmax_ext)
    tmax_table  <- rbind(tmax_table, temp_table)
  }
}

# From 2013 to 2022
fi <- data.frame(temp_date = seq(as.Date("2013-01-01"), length=120, by="month"))  
ff <- data.frame(temp_date = seq(as.Date("2013-02-01"), length=120, by="month"))-1  

f_start <- fi %>% 
  separate(temp_date, into = c("day", "month", "year")) %>%
  mutate(good_date = paste(year, month, day, sep = "/")) %>%
  dplyr::select(good_date)

f_end <- ff %>% 
  separate(temp_date, into = c("day", "month", "year")) %>%
  mutate(good_date = paste(year, month, day, sep = "/")) %>%
  dplyr::select(good_date)

dates <- f_end %>% 
  separate(good_date, into = c("day", "month", "year")) # separate dates to save them in each iteration

foreach (j=1:length(coords_list)) %dopar% {
  lat <- coords_list$lat[j]
  long <- coords_list$lon[j]
  mdt <- microclima::get_dem(lat = lat, long = long, resolution = 30) ## microclimate modelled across 5x5 km square around coord  
  for (i in 1:nrow(f_end))  {
    temp <- runauto(mdt,
                    f_start[i,], f_end[i,], 
                    hgt = 0.1, 
                    l = NA, x = NA,
                    coastal = FALSE,
                    habitat = 7, 
                    plot.progress = FALSE, save.memory = FALSE)
    tmax <- temp$tmax
    a <- coords_sf[j,]
    tmax_ext <- raster::extract(tmax, a, fun=mean) # extract tmax value from coordinates    
    
    ## save tmax output to table
    temp_table <- data.frame(matrix(ncol = 5, nrow = 1))
    colnames(temp_table) <- colnames(tmax_table)
    temp_table[1,] <- c(coords_list$Code[j],dates$month[i],dates$year[i],"recent",tmax_ext)
    tmax_table  <- rbind(tmax_table, temp_table)
  }
}

## prepare tmax table to export adding information on season
tmax_table <- tmax_table %>% 
  mutate(across(month:year, as.numeric)) %>% 
  mutate(across(tmax, as.numeric)) %>% 
  mutate(across(tmax, round, 3))
  
rm(coords_list, coords_sf, dates, f_end, f_start, ff, fi, temp_table)

##### 2. Sample Coverage and species richness (iNEXT) ####
iNext_size_based <-read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "2. iNEXT_spp_richness")

# iNEXT_size_based dataframe comes from using iNEXT() function as shown below using abundances matrices for each study period (historical and recent)
# Abundances matrices are not provided as we don't have the rights to share them publicly
# q means which Hill number we are interested in (q=0)
# knots by default is 40, we are using 300 to increase the number of times the model is run
# Example of iNEXT function used: Output <- iNEXT(abundance_table, q=0, knots = 300, datatype="abundance") 

## Sample coverage filters calculation
### get maximum Sample Coverage (SC) for each assemblage in the historical survey
max_SC_hist <- iNext_size_based %>% 
  filter(Survey=="Historical") %>% 
  group_by(Code) %>% 
  summarise(Max=max(SC)) %>% 
  arrange(Max)

filter1 <- max_SC_hist$Max[1] # lowest SC filter which is the lowest max. SC for the historical survey (SC >= 0.66)
filter2 <- max_SC_hist$Max[2] # Intermediate SC filter which is the second lowest max. SC for the historical survey (SC >= 0.75)

## get maximum Sample Coverage (SC) for each assemblage in the recent survey
max_SC_rec <- iNext_size_based %>% 
  filter(Survey=="Recent") %>% 
  group_by(Code) %>% 
  summarise(Max=max(SC)) %>% 
  arrange(Max)

filter3 <- max_SC_rec$Max[1] # highest SC filter which is the lowest max. SC for the recent survey (SC >= 0.93)

## Calculate species richness based on SC filters
### SC filter 1 and filter 2
SC_filters <- c(filter1, filter2, filter3)

qD_table <- data.frame() # table to save estimated species richness

for (i in SC_filters) {
  
  Code_to_remove <- max_SC_hist %>% 
    filter(Max < i) %>% 
    dplyr::select(Code) %>% 
    unlist() ## list of locations with maximum SC above the SC filter, which will be removed from data
  
  filter_data <- iNext_size_based %>% 
    filter(SC>=i) %>% 
    arrange(SC) %>% # Arrange by SC to get the closest SC at the top of each site & survey
    group_by(Code, Survey) %>% 
    filter(row_number()==1) %>% # Get first row of data per site and survey, so the closest value above the SC filter
    filter(!Code %in% Code_to_remove) # remove locations with max. SC below the applied threshold
  
  filter_data$SC_filter <- i
  
  qD_table <- rbind(qD_table, filter_data)
  
}

qD_table # contains estimated species richness (qD) for each survey and for each Sample Coverage filter used

# change filters nomenclature
qD_table$SC_filter <- as.character(as.numeric(qD_table$SC_filter)) 

qD_table$SC_filter <- gsub("0.66", "qD1", qD_table$SC_filter)
qD_table$SC_filter <- gsub("0.75", "qD2", qD_table$SC_filter)
qD_table$SC_filter <- gsub("0.93", "qD3", qD_table$SC_filter)

# Number of sites where each method (extrapolation, observed or rarefaction) has been applied to calculate qD for each filter applied
qD_table %>% 
  group_by(SC_filter, Survey, Method) %>% 
  summarise(n=n())

## export table
write_xlsx(qD_table, paste0(directory,"/Tables/spp_richness_iNEXT_filtered_by_SC.xlsx"))

# Estimated species richness (qD) obtained in qD_table is included in sheet 3 (Main_analysis_data) in supp. material in several columns with names as:
## qD_SC_filter1 for the lowest filter on SC
## qD_SC_filter2 for the intermediate filter on SC
## qD_SC_filter3 for the highest filter on SC

##### 3. Paired test on observed Sample Coverage (SC) and species richness (qD) #####
###### 3.1 Observed SC ####
Obs_SC <- read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "2. iNEXT_spp_richness") %>% 
  filter(Method=="Observed") %>% 
  dplyr::select(Code, Method, SC, Survey, Region) %>% 
  pivot_wider(names_from = Survey, values_from = SC)

Regions <- c("Gredos", "Guadarrama", "Meridional", "Javalambre")

table1 <- data.frame(
  Region = character(),
  mean_hist = numeric(),
  SD_hist = numeric(),
  mean_rec = numeric(),
  SD_rec = numeric(),
  stringsAsFactors = FALSE # to save results
)

for (i in Regions) {
  print(i)
  
  filter_data <- Obs_SC %>% 
    filter(Region==i )
  
  res <- shapiro.test(filter_data$Recent - filter_data$Historical) 
  print(res)
  
  filter_data <- filter_data %>% 
    group_by(Region) %>% 
    summarise(mean_hist=mean(Historical), 
              SD_hist=sd(Historical),
              mean_rec=mean(Recent),
              SD_rec=sd(Recent)) %>% 
    mutate(across(where(is.numeric), ~ round(., 3)))
  
  table1 <- rbind(table1, filter_data)
  
} ## All normally distributed, except Guadarrama

paired_test_results <- Obs_SC %>%
  group_by(Region) %>%
  summarise(
    p_value = t.test(Recent, Historical, paired = TRUE)$p.value,
    statistic = t.test(Recent, Historical, paired = TRUE)$statistic) %>% 
  mutate(paired_test="t_test",
         variable="Obs_SC") %>% 
  filter(Region!="Guadarrama")

paired_wilcoxon_results <- Obs_SC %>%
  filter(Region=="Guadarrama") %>% 
  group_by(Region) %>%
  summarise(
    p_value = wilcox.test(Recent, Historical, paired = TRUE)$p.value,    
    statistic = wilcox.test(Recent, Historical, paired = TRUE)$statistic) %>% 
  mutate(paired_test="wilcoxon",
         variable="Obs_SC") 

paired_test_SC <- rbind(paired_test_results, paired_wilcoxon_results)

## overall
filter_data <- Obs_SC %>% 
  group_by(Method) %>% 
  summarise(mean_hist=mean(Historical), 
            SD_hist=sd(Historical),
            mean_rec=mean(Recent),
            SD_rec=sd(Recent)) %>% 
  mutate(Region="overall") %>% 
  dplyr::select(-Method) %>% 
  relocate(Region, .before = mean_hist) %>% 
  mutate(across(where(is.numeric), ~ round(., 3)))

table1 <- rbind(table1, filter_data)

shapiro.test(Obs_SC$Recent - Obs_SC$Historical) # not normally distributed

paired_wilcoxon_results <- Obs_SC %>%
  summarise(
    p_value = wilcox.test(Recent, Historical, paired = TRUE)$p.value,    
    statistic = wilcox.test(Recent, Historical, paired = TRUE)$statistic) %>% 
  mutate(paired_test="wilcoxon",
         variable="Obs_SC",
         Region="overall") %>% 
  relocate(Region, .before = "p_value")

paired_test_SC <- rbind(paired_test_SC, paired_wilcoxon_results)

table1 <- left_join(table1, paired_test_SC, by="Region")

###### 3.1 Observed qD ####
Obs_qD <- read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "2. iNEXT_spp_richness") %>% 
  filter(Method=="Observed") %>% 
  dplyr::select(Code, Method, qD, Survey, Region) %>% 
  pivot_wider(names_from = Survey, values_from = qD)

Regions <- c("Gredos", "Guadarrama", "Meridional", "Javalambre")

table_temp <- data.frame(
  Region = character(),
  mean_hist = numeric(),
  SD_hist = numeric(),
  mean_rec = numeric(),
  SD_rec = numeric(),
  stringsAsFactors = FALSE # to save results
)

for (i in Regions) {
  print(i)
  
  filter_data <- Obs_qD %>% 
    filter(Region==i)
  
  res <- shapiro.test(filter_data$Recent - filter_data$Historical) 
  print(res)
  
  filter_data <- filter_data %>% 
    group_by(Region) %>% 
    summarise(mean_hist=mean(Historical), 
              SD_hist=sd(Historical),
              mean_rec=mean(Recent),
              SD_rec=sd(Recent)) %>% 
    mutate(across(where(is.numeric), ~ round(., 3)))
  
  table_temp <- rbind(table_temp, filter_data)
  
} ## All normally distributed

paired_test_results <- Obs_qD %>%
  group_by(Region) %>%
  summarise(
    p_value = t.test(Recent, Historical, paired = TRUE)$p.value,
    statistic = t.test(Recent, Historical, paired = TRUE)$statistic) %>% 
  mutate(paired_test="t_test",
         variable="Obs_qD")

## overall
filter_data <- Obs_qD %>% 
  group_by(Method) %>% 
  summarise(mean_hist=mean(Historical), 
            SD_hist=sd(Historical),
            mean_rec=mean(Recent),
            SD_rec=sd(Recent)) %>% 
  mutate(Region="overall") %>% 
  dplyr::select(-Method) %>% 
  relocate(Region, .before = mean_hist) %>% 
  mutate(across(where(is.numeric), ~ round(., 3)))

table_temp <- rbind(table_temp, filter_data)

shapiro.test(Obs_qD$Recent - Obs_qD$Historical) # not normally distributed

paired_wilcoxon_results <- Obs_qD %>%
  summarise(
    p_value = wilcox.test(Recent, Historical, paired = TRUE)$p.value,    
    statistic = wilcox.test(Recent, Historical, paired = TRUE)$statistic) %>% 
  mutate(paired_test="wilcoxon",
         variable="Obs_qD",
         Region="overall") %>% 
  relocate(Region, .before = "p_value")

paired_test_results <- rbind(paired_test_results, paired_wilcoxon_results)

table_temp <- left_join(table_temp, paired_test_results, by="Region")

table1 <- rbind(table1, table_temp)

write_xlsx(table1, paste0(directory, "/Tables/Table1_main_paper_SC_qD.xlsx"))

rm(paired_test_results, paired_test_SC, paired_wilcoxon_results, filter_data, res, table_temp,
   filter1, filter2, filter3, SC_filters, Code_to_remove, max_SC_hist, max_SC_rec, Obs_qD, Obs_SC)

##### 4. Collinearity analysis #####
directory <- rstudioapi::getActiveDocumentContext()$path
directory <- dirname(directory)

data_def <-read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "3. Main_analysis_data")

perform_spearman_test <- function(data, x_var, y_var, test_name) {
  shapiro_x <- shapiro.test(data[[x_var]])
  shapiro_y <- shapiro.test(data[[y_var]])
  
  spearman_test <- cor.test(data[[x_var]], data[[y_var]], 
                            method = "spearman", exact = FALSE)
  
  result <- tidy(spearman_test) %>%
    mutate(
      test_name = test_name,
      shapiro_x_p = shapiro_x$p.value,
      shapiro_y_p = shapiro_y$p.value,
      n = nrow(data)
    ) %>%
    dplyr::select(test_name, estimate, p.value, n, shapiro_x_p, shapiro_y_p)
  
  return(result)
}

# Initialize results dataframe
spearman_results <- data.frame()

## Tmax tests
# Test 1: Change in Tmax vs historical Tmax (10y method)
filter_data <- data_def %>% 
  filter(variable == "CTI_IB" & Survey == "Historical") %>% 
  group_by(Code, Region, historical_10y, diff_10y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "historical_10y", "diff_10y", 
                                                "Change in Tmax vs historical Tmax (10y method)"))

# Test 2: Change in Tmax vs historical Tmax (survey years method)
filter_data <- data_def %>% 
  filter(variable == "CTI_IB" & Survey == "Historical") %>% 
  group_by(Code, Region, historical_surv_y, diff_surv_y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "historical_surv_y", "diff_surv_y", 
                                                "Change in Tmax vs historical Tmax (survey years method)"))

# Test 3: Change in Tmax for both methods (10 year and survey years)
filter_data <- data_def %>% 
  filter(variable == "CTI_sd_EU" & Survey == "diff") %>% 
  group_by(Code, Region, diff_surv_y, diff_10y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "diff_surv_y", "diff_10y", 
                                                "Change in Tmax for both methods (10y and survey years)"))

# Test 4: Change in Tmax for both methods (10 year and survey years)
filter_data <- data_def %>% 
  filter(variable == "CTI_sd_EU" & Survey == "Historical") %>% 
  group_by(Code, Region, historical_10y, historical_surv_y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "historical_10y", "historical_surv_y", 
                                                "Historical Tmax for both methods (10y and survey years)"))

## CTI tests
# Test 5: Change in Iberian CTI vs change in Tmax (10y method)
filter_data <- data_def %>% 
  filter(variable == "CTI_IB" & Survey == "diff") %>% 
  group_by(Code, Region, median, diff_10y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "median", "diff_10y", 
                                                "Change in Iberian CTI vs change in Tmax (10y method)"))

# Test 5: Change in Iberian CTI vs change in Tmax (survey years method)
filter_data <- data_def %>% 
  filter(variable == "CTI_IB" & Survey == "diff") %>% 
  group_by(Code, Region, median, diff_surv_y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "median", "diff_surv_y", 
                                                "Change in Iberian CTI vs change in Tmax (survey years method)"))

# Test 6: Change in European CTI vs change in Tmax (10y method)
filter_data <- data_def %>% 
  filter(variable == "CTI_EU" & Survey == "diff") %>% 
  group_by(Code, Region, median, diff_10y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "median", "diff_10y", 
                                                "Change in European CTI vs change in Tmax (10y method)"))

# Test 6: Change in European CTI vs change in Tmax (survey years method)
filter_data <- data_def %>% 
  filter(variable == "CTI_EU" & Survey == "diff") %>% 
  group_by(Code, Region, median, diff_surv_y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "median", "diff_surv_y", 
                                                "Change in European CTI vs change in Tmax (survey years method)"))

# Test 7: Change in Iberian CTIsd vs change in Tmax (10y method)
filter_data <- data_def %>% 
  filter(variable == "CTI_sd_IB" & Survey == "diff") %>% 
  group_by(Code, Region, median, diff_10y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "median", "diff_10y", 
                                                "Change in Iberian CTIsd vs change in Tmax (10y method)"))

# Test 7: Change in Iberian CTIsd vs change in Tmax (survey years method)
filter_data <- data_def %>% 
  filter(variable == "CTI_sd_IB" & Survey == "diff") %>% 
  group_by(Code, Region, median, diff_surv_y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "median", "diff_surv_y", 
                                                "Change in Iberian CTIsd vs change in Tmax (survey years method)"))

# Test 8: Change in European CTIsd vs change in Tmax (10y method)
filter_data <- data_def %>% 
  filter(variable == "CTI_sd_EU" & Survey == "diff") %>% 
  group_by(Code, Region, median, diff_10y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "median", "diff_10y", 
                                                "Change in European CTIsd vs change in Tmax (10y method)"))

# Test 8: Change in European CTIsd vs change in Tmax (survey years method)
filter_data <- data_def %>% 
  filter(variable == "CTI_sd_EU" & Survey == "diff") %>% 
  group_by(Code, Region, median, diff_surv_y) %>% 
  summarise(n = n(), .groups = "drop")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(filter_data, "median", "diff_surv_y", 
                                                "Change in European CTIsd vs change in Tmax (survey years method)"))

# Test for STIs
STIs_table <- read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "5. Occupancy_analysis")

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(STIs_table, "IBER_mean_temp", "EU_temp.mean", 
                                                "Iberian vs European STIs (mean)"))

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(STIs_table, "IBER_sd_temp", "EU_temp.sd", 
                                                "Iberian vs European STIs (sd)"))

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(STIs_table, "IBER_mean_temp", "IBER_sd_temp", 
                                                "Iberian STI_mean vs STI_sd"))

spearman_results <- rbind(spearman_results, 
                          perform_spearman_test(STIs_table, "EU_temp.mean", "EU_temp.sd", 
                                                "European STI_mean vs STI_sd"))

# Print and save results
print(spearman_results)
write_xlsx(spearman_results, paste0(directory,"/Tables/spearman_correlation_results.xlsx"))

##### 5. Temperature methods analysis (Figure S1) ####
directory <- rstudioapi::getActiveDocumentContext()$path
directory <- dirname(directory)

data_plot <-read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "7. Microclima_data")

mytheme <- theme(axis.line = element_line(size=.2),
                 axis.text.x = element_text(size=11),
                 axis.title.x = element_text(size=13, vjust=-2),
                 axis.text.y = element_text(size=11),
                 axis.title.y = element_text(size=13, vjust = 2.5),
                 axis.ticks.y = element_line(size=.2),
                 axis.ticks.x = element_line(size=.2),
                 legend.key.width = unit(1.2, "cm"),
                 legend.text = element_text(size=11),
                 legend.title = element_blank(),
                 legend.position = "none",
                 plot.margin = unit(c(.4,.4,.4,.4), "cm")) # theme to use for plotting

data_plot_red <- data_plot %>% 
  group_by(Code, Region, Survey, method) %>% 
  summarise(mean_tmax=mean(mean_tmax),
            tmax_sd=mean(tmax_sd)) %>% 
  mutate(year_per=method) %>% 
  relocate(year_per, .before=Survey) %>% 
  relocate(method, .after=tmax_sd) ## calculate totals per method and region

data_plot_red <- rbind(data_plot, data_plot_red)

sort(unique(data_plot_red$year_per))

## Gredos
filter_data <- data_plot_red %>% 
  filter(Region=="Gredos")

filter_data$year_per <- factor(filter_data$year_per, levels=c("1984_1985", "1985_1986", "1980_1989",
                                                              "2020_2021","2021_2022","2013_2022", "surveyed years","10-year period"))

filter_data$method <- factor(levels=c("surveyed years","10-year period"), filter_data$method)

gre <- filter_data %>% 
  ggplot(aes(x=year_per, y=mean_tmax)) +
  geom_boxplot(aes(color=Survey, fill=method), size=.4, width=.6) +
  scale_y_continuous(limits = c(15,29)) +
  ylab(bquote("Mean T"["max"])) +
  scale_color_manual(values = c("grey20","grey60")) +
  scale_fill_manual(values = c("white","grey")) +
  ggtitle("Sierra de Gredos") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, size=10, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.line = element_line(size=.2),
        plot.title = element_text(size=12, face="bold"),
        axis.ticks = element_line(size=.2),
        strip.background = element_blank())

## Guadarrama 
filter_data <- data_plot_red %>% 
  filter(Region=="Guadarrama")

filter_data$year_per <- factor(filter_data$year_per, levels=c("2003_2004", "2004_2005", "1980_1989",
                                                              "2016_2017","2020_2021","2013_2022","surveyed years","10-year period"))

filter_data$method <- factor(levels=c("surveyed years","10-year period"), filter_data$method)

gua <- filter_data %>% 
  ggplot(aes(x=year_per, y=mean_tmax)) +
  geom_boxplot(aes(color=Survey, fill=method), size=.4, width=.6) +
  scale_y_continuous(limits = c(15,29)) +
  ylab(bquote("Mean T"["max"])) +
  scale_color_manual(values = c("grey20","grey60")) +
  scale_fill_manual(values = c("white","grey")) +
  ggtitle("Sierra de Guadarrama") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, size=10, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.line = element_line(size=.2),
        plot.title = element_text(size=12, face="bold"),
        axis.ticks = element_line(size=.2),
        strip.background = element_blank())

## Meridional 
filter_data <- data_plot_red %>% 
  filter(Region=="Meridional")

filter_data$year_per <- factor(filter_data$year_per, levels=c("1985_1986", "1980_1989",
                                                              "2019_2020","2020_2021","2021_2022", "2013_2022", "surveyed years","10-year period"))

filter_data$method <- factor(levels=c("surveyed years","10-year period"), filter_data$method)

mer <- filter_data %>% 
  ggplot(aes(x=year_per, y=mean_tmax)) +
  geom_boxplot(aes(color=Survey, fill=method), size=.4, width=.6) +
  scale_y_continuous(limits = c(15,29)) +
  ylab(bquote("Mean T"["max"])) +
  scale_color_manual(values = c("grey20","grey60")) +
  scale_fill_manual(values = c("white","grey")) +
  ggtitle("Sistema Ibérico Meridional") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, size=10, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.line = element_line(size=.2),
        plot.title = element_text(size=12, face="bold"),
        axis.ticks = element_line(size=.2),
        strip.background = element_blank())

## Javalambre
filter_data <- data_plot_red %>% 
  filter(Region=="Javalambre")

filter_data$year_per <- factor(filter_data$year_per, levels=c("1990_1991", "1980_1989",
                                                              "2019_2020","2020_2021","2021_2022", "2013_2022","surveyed years","10-year period"))

filter_data$method <- factor(levels=c("surveyed years","10-year period"), filter_data$method)

jav <- filter_data %>% 
  ggplot(aes(x=year_per, y=mean_tmax)) +
  geom_boxplot(aes(color=Survey, fill=method), size=.4, width=.6) +
  scale_y_continuous(limits = c(15,29)) +
  ylab(bquote("Mean T"["max"])) +
  scale_color_manual(values = c("grey20","grey60")) +
  scale_fill_manual(values = c("white","grey")) +
  ggtitle("Javalambre") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, size=10, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.line = element_line(size=.2),
        plot.title = element_text(size=12, face="bold"),
        axis.ticks = element_line(size=.2),
        strip.background = element_blank())

## export plot
Fig_S1 <- ggarrange(gre,mer,gua,jav, ncol=2, nrow=2, common.legend = T, align = "hv")
plot(Fig_S1)

tiff(paste0(directory,"/Figures/FigS1_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)
plot(Fig_S1)
dev.off()

rm(data_plot, data_plot_red, filter_data, gre,mer,gua,jav)

#### STATISTICAL ANALYSIS ####
## common data for all statistical analysis
directory <- rstudioapi::getActiveDocumentContext()$path
directory <- dirname(directory)

data_def <-read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "3. Main_analysis_data")

## order regions
regions_order <- c("Gredos", "Guadarrama", "Meridional", "Javalambre")
data_def$Region <- factor(levels = regions_order, data_def$Region)

mytheme <- theme(axis.line = element_line(size=.2),
                 axis.text.x = element_text(size=11),
                 axis.title.x = element_text(size=13, vjust=-2),
                 axis.text.y = element_text(size=11),
                 axis.title.y = element_text(size=13, vjust = 2.5),
                 axis.ticks.y = element_line(size=.2),
                 axis.ticks.x = element_line(size=.2),
                 legend.key.width = unit(1.2, "cm"),
                 legend.text = element_text(size=11),
                 legend.title = element_blank(),
                 legend.position = "none",
                 plot.margin = unit(c(.4,.4,.4,.4), "cm")) # theme to use for plotting

##### 1. Paired test on temporal shifts on main variables ####
directory <- rstudioapi::getActiveDocumentContext()$path
directory <- dirname(directory)

diff_table <-read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "6. diff_table")

###### 1.1 Paired test (t-test or Wilcoxon test) ####
## create functions to run
perform_paired_test <- function(data, test_type) {
  if (test_type == "t_test") {
    res <- t.test(data$Recent, data$Historical, paired = TRUE)
  } else {
    res <- wilcox.test(data$Recent, data$Historical, paired = TRUE)
  }
  return(tidy(res))
}

paired_test_loop <- function(data, variable_name) {
  regions <- c("Gredos", "Guadarrama", "Meridional", "Javalambre", "Overall")
  
  results <- data.frame()
  
  for (region in regions) {
    if (region == "Overall") {
      subset_data <- data
    } else {
      subset_data <- data %>% filter(Region == region)
    }
    
    normality_test <- shapiro.test(subset_data$diff)
    test_type <- if (normality_test$p.value > 0.05) "t_test" else "wilcoxon"
    
    test_results <- perform_paired_test(subset_data, test_type)
    
    results <- rbind(results, data.frame(
      Region = region,
      p_value = test_results$p.value,
      statistic = test_results$statistic,
      paired_test = test_type,
      variable = variable_name
    ))
  }
  
  return(results)
}

# run functions for paired tests
variables <- sort(unique(diff_table$variable)) # 11 variables in total
paired_test_res <- data.frame() # to save results

for (var in variables) {
  filter_data <- diff_table %>% filter(variable == var)
  results <- paired_test_loop(filter_data, var)
  paired_test_res <- rbind(paired_test_res, results)
}

## prepare data to plot
kk <- diff_table %>% 
  group_by(Region, variable) %>% 
  summarise_at(vars(Historical:diff), list(mean = mean, sd = sd), na.rm = TRUE)

kk_ov <- diff_table %>% 
  group_by(variable) %>% 
  summarise_at(vars(Historical:diff), list(mean = mean, sd = sd), na.rm = TRUE) %>% 
  mutate(Region="Overall") %>% 
  relocate(Region, .before="variable")

plot_data <- rbind(kk, kk_ov)
rm(kk, kk_ov)

plot_data <- left_join(plot_data, paired_test_res, by=c("Region", "variable")) %>% 
  mutate(diff_sig =ifelse(p_value<0.05 & diff_mean<0,"sig_decrease",
                          ifelse(p_value<0.05 & diff_mean>0,"sig_increase","non_sig"))) %>% 
  relocate(Historical_sd, .after = Historical_mean) %>% 
  relocate(Recent_sd, .after = Recent_mean)

plot_data$diff_sig <- factor(levels=c("sig_decrease", "sig_increase","non_sig"), plot_data$diff_sig)
plot_data$Region <- factor(levels=c("Javalambre","Meridional","Guadarrama","Gredos","Overall"), plot_data$Region)

###### 1.2 Export table with paired test results (Table S6) ####
export_t <- plot_data %>% 
  mutate_if(is.numeric, round, 3) %>%  
  mutate(historical=paste0(Historical_mean," (±",Historical_sd,")"),
         recent=paste0(Recent_mean," (±",Recent_sd,")"),
         diff=paste0(diff_mean," (±",diff_sd,")")) %>% 
  dplyr::select(Region, variable, historical, recent, diff, paired_test, statistic, diff_sig)

p_values <- plot_data %>% 
  dplyr::select(p_value)
p_values <- p_values[,-1]
export_t <- cbind(export_t, p_values)

export_t %>% 
  group_by(variable) %>% 
  summarise(n=n()) # 5 rows per variable (one for each region and one overall)

export_t$variable <- factor(levels = c("CTI_IB", "CTI_EU","CTI_sd_IB", "CTI_sd_EU", "STI_sd_IB", "STI_sd_EU",
                                       "qD_filter1", "qD_filter2", "qD_filter3", "Annual_Tmax_10y", "Annual_Tmax_surv_y"),
                            export_t$variable)

export_t$Region <- factor(levels = c("Gredos", "Guadarrama", "Meridional", "Javalambre", "Overall"),
                            export_t$Region)

export_t <- export_t %>% 
  arrange(Region) %>% 
  arrange(variable)

write_xlsx(export_t, path = paste0(directory, "/Tables/TableS6_DDI-2024-0174.xlsx"))

rm(p_values, paired_test_res, filter_data)

###### 1.3 Paired test figure (Figure 3) #####
vars <- c("CTI_IB", "CTI_sd_IB", "qD_filter2", "STI_sd_IB")
Fig_3 <- plot_data %>% 
  filter(variable %in% vars) %>% 
  ggplot(aes(x=diff_mean, y=Region, color=diff_sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color="grey") +
  geom_point(size=3) +
  geom_segment(aes(x = diff_mean - diff_sd, xend = diff_mean + diff_sd, y = Region, yend = Region), size=.5) +
  scale_color_manual(values = c("blue","red","grey")) +
  mytheme +
  facet_wrap(vars(variable), scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=11, angle=0),
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = .4),
        axis.ticks = element_line(linewidth =.4),
        strip.background = element_blank(),
        strip.text = element_text(size=13, face="bold"),
        panel.spacing.x = unit(1, "lines"),
        legend.position = "none",
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))

plot(Fig_3)

## export plot 
tiff(paste0(directory,"/Figures/Fig3_DDI-2024-0174.tif"), 
     width = 30, height = 10, units = "cm", res = 300)
plot(Fig_3)
dev.off()

###### 1.4 Paired test European STI and other richness filters (Figure S9) #####
vars <- c("CTI_EU", "CTI_sd_EU", "STI_sd_EU", "qD_filter1", "qD_filter2", "qD_filter3")
Fig_S9 <- plot_data %>% 
  filter(variable %in% vars) %>% 
  ggplot(aes(x=diff_mean, y=Region, color=diff_sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color="grey") +
  geom_point(size=3) +
  geom_segment(aes(x = diff_mean - diff_sd, xend = diff_mean + diff_sd, y = Region, yend = Region), size=.5) +
  scale_color_manual(values = c("blue","red","grey")) +
  mytheme +
  facet_wrap(vars(variable), scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=11, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = .4),
        axis.ticks = element_line(linewidth =.4),
        strip.background = element_blank(),
        strip.text = element_text(size=13, face="bold"),
        panel.spacing.x = unit(1, "lines"),
        legend.position = "none",
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))

plot(Fig_S9)

## export plot 
tiff(paste0(directory,"/Figures/FigS9_DDI-2024-0174.tif"), 
     width = 30, height = 10, units = "cm", res = 300)
plot(Fig_S9)
dev.off()

###### 1.5 Paired test on Tmax (Figure S4) ####
vars <- c("Annual_Tmax_10y", "Annual_Tmax_surv_y")
Fig_S4 <- plot_data %>% 
  filter(variable %in% vars) %>% 
  ggplot(aes(x=diff_mean, y=Region, color=diff_sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color="grey") +
  geom_point(size=3) +
  geom_segment(aes(x = diff_mean - diff_sd, xend = diff_mean + diff_sd, y = Region, yend = Region), size=.5) +
  scale_color_manual(values = c("blue","grey")) +
  mytheme +
  facet_wrap(vars(variable), scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=11, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = .4),
        axis.ticks = element_line(linewidth =.4),
        strip.background = element_blank(),
        strip.text = element_text(size=13, face="bold"),
        panel.spacing.x = unit(1, "lines"),
        legend.position = "none",
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))

plot(Fig_S4)

## export plot 
tiff(paste0(directory,"/Figures/FigS4_DDI-2024-0174.tif"), 
     width = 20, height = 10, units = "cm", res = 300)
plot(Fig_S4)
dev.off()

rm(var, variables, vars, plot_data, results)

##### 2. Distance matrix for spatial autocorrelation test #####
coords <- read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "1. coordinates_list") %>% 
  dplyr::select(Code, lon, lat)

filter_data <- data_def %>% 
  filter(variable=="CTI_IB" & Survey=="Historical")

filter_data <- left_join(filter_data, coords, by="Code")
filter_data <- st_as_sf(filter_data, coords = c("lon","lat")) # convert dataframe to spatial data
st_crs(filter_data) <- st_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # define crs

## create neighbor list from points using k nearest as criteria:
knea <- knearneigh(filter_data,k=5) ## k means the number of neighbors to count 
neib <- knn2nb(knea)

## assign weights to the neighbors
lw <- nb2listw(neib, style="W", zero.policy=TRUE) # lw will be use for spatial autocorrelation test

## visualise output from neighbor list
plot(st_geometry(filter_data), border = "lightgray")
plot.nb(neib, st_geometry(filter_data), add = TRUE)

rm(neib, knea, coords)

##### FUNCTIONS AND TABLES for GLMs and LMERs ####
## empty dataframes to save models ouputs from dredge function and coefficients from model averaging
A_model_outputs_exp <- data.frame(
  variable = character(),
  formula = character(),
  R.2 = numeric(),
  adjR.2 = numeric(),
  df = numeric(),
  logLik = numeric(),
  AICc = numeric(),
  delta = numeric(),
  weight = numeric(),
  stringsAsFactors = FALSE)

A_coefs_export <-  data.frame(
  response_variable = character(),
  select_variables = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  t_value = numeric(),
  p_value = numeric(), 
  stringsAsFactors = FALSE)

A_spatial_autocorr_residuals <- data.frame(
  model = character(),
  statistic = numeric(),
  Obs_rank = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

# to save model outputs from dredge function (more than one model below delta AICc < 2)
get_filtered_models_GLM <- function(dredge.model, response_variable_name) {
  # Get models info
  filtered_models <- data.frame(dredge.model) %>% 
    filter(delta < 2) %>%
    dplyr::select('R.2':'weight') %>% 
    mutate(variable = response_variable_name) %>% 
    relocate(variable, .before = 'R.2')
  
  # Get model formulas
  subset_models <- get.models(dredge.model, subset = delta < 2)
  subset_formulas <- lapply(subset_models, formula)
  subset_formulas <- sapply(subset_formulas, function(x) Reduce(paste, deparse(x)))
  filtered_models$formula <- subset_formulas
  filtered_models <- filtered_models %>%
    relocate(formula, .after = variable)
  
  return(filtered_models)
}

get_filtered_models_LMER <- function(dredge.model, response_variable_name) {
  # Get models info
  filtered_models <- data.frame(dredge.model) %>% 
    filter(delta < 2) %>%
    dplyr::select('R.2':'weight') %>% 
    mutate(variable = response_variable_name) %>% 
    relocate(variable, .before = 'R.2')
  
  # Get model formula
  subset_models <- get.models(dredge.model, subset = delta < 2)
  m1 <- subset_models[[1]]
  filtered_models$formula <- as.character(formula(m1))[3]
  filtered_models <- filtered_models %>%
    relocate(formula, .after = variable)
  
  return(filtered_models)
}

get_coefficients_LMER <- function(model, response_variable_name) {
  coefs <- data.frame(coef(summary(model))) # Get model coefficients
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value))) # Calculate p-values using normal distribution approximation
  coefs <- coefs %>% 
    mutate(response_variable = response_variable_name,
           selected_variables = rownames(coefs)) # Add response variable and selected variables
  rownames(coefs) <- NULL # Remove rownames
  coefs <- coefs %>% 
    relocate(response_variable:selected_variables, .before = Estimate) # Reorder columns
  colnames(coefs) <- colnames(A_coefs_export)
  return(coefs)
}

# to save coefficients from averaged models (more than one model below delta AICc < 2)
get_avg_model_coefs  <- function(avg_model, response_variable_name) {
  
  # Get summary and convert to data frame
  avg_model_summ <- summary(avg_model)$coefmat.full
  avg_model_summ <- data.frame(avg_model_summ[, -3])
  avg_model_summ <- avg_model_summ %>% 
    mutate(response_variable = response_variable_name) # Add response variable name
  avg_model_summ$selected_variables <- rownames(avg_model_summ) # Add selected variables as a column
  rownames(avg_model_summ) <- NULL # Remove rownames
  avg_model_summ <- avg_model_summ %>% 
    relocate(response_variable:selected_variables, .before = Estimate) # Reorder columns
  colnames(avg_model_summ) <- colnames(A_coefs_export)
  return(avg_model_summ)
}

# When only one model in delta AICc < 2
get_best_model_info <- function(dredge.model, response_variable_name) {
  # Get model info
  filtered_models <- data.frame(dredge.model) %>% 
    filter(delta < 2) %>%
    dplyr::select('R.2':'weight') %>% 
    mutate(variable = response_variable_name) %>% 
    relocate(variable, .before = 'R.2')
  
  # Get model formula
  best_model <- get.models(dredge.model, subset = delta < 2)[[1]]
  filtered_models$formula <- as.character(formula(best_model))[3]
  
  filtered_models <- filtered_models %>%
    relocate(formula, .after = variable)
  
  return(filtered_models)
}

# When only one model in delta AICc < 2
get_best_model_coefs <- function(model, response_variable_name) {
  coefs <- data.frame(summary(model)$coefficients) # Get model coefficients
  coefs <- coefs %>% 
    mutate(response_variable = response_variable_name,
           selected_variables = rownames(coefs)) # Add response variable and selected variable
  rownames(coefs) <- NULL
  coefs <- coefs %>% 
    relocate(response_variable:selected_variables, .before = Estimate) # Reorder columns
  colnames(coefs) <- colnames(A_coefs_export)
  return(coefs)
}

##### 3. CTI and spp richness vs local temperature ####
###### 3.1 Figure 4 manuscript #####
filter_data <- data_def %>% 
  filter(variable=="CTI_IB" & Survey=="Recent") %>% 
  dplyr::select(Code, median) %>% 
  arrange(median)

min_val <- filter_data$median[[1]] ## minimum value in recent for plotting

## historical period
filter_data <- data_def %>% 
  filter(variable=="CTI_IB" & Survey=="Historical") %>% 
  dplyr::select(Code, median, q25,q75,historical_10y, Region)

m1 <- glm(median ~ scale(historical_10y)*Region,
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get models info
filtered_models <- get_filtered_models_GLM(dredge.model, "CTI_IB_hist")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "CTI_IB_hist")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

# predict values from averaged model
filter_data$predicted <- predict(avg_model, newdata=filter_data)
filter_data$residuals <- filter_data$median - filter_data$predicted

# spatial autocorrelation test
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("CTI_IB_hist",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

# plot results
plot_CTI_hist <- filter_data %>% 
  ggplot(aes(x=historical_10y, y=median, fill=Region)) +
  geom_hline(yintercept = min_val, linetype="dashed", color="grey", linewidth=.5) +
  geom_pointrange(aes(ymin=q25, ymax=q75, fill=Region, shape=Region), size=.5, stroke = .2, colour="grey50", alpha=.8) +
  geom_line(aes(y=predicted, linetype=Region, color=Region), linewidth=.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0)) +
  xlab(bquote("Historical T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(10.5, 13), "CTI", expand = c(0,0)) +
  mytheme +
  theme(axis.title.x = element_blank())

## Recent data
filter_data <- data_def %>% 
  filter(variable=="CTI_IB" & Survey=="Recent") %>% 
  dplyr::select(Code, median, q25, q75, recent_10y, Region)

m1 <- glm(median ~ scale(recent_10y)*Region,
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get models info
filtered_models <- get_filtered_models_GLM(dredge.model, "CTI_IB_rec")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "CTI_IB_rec")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

# predict values from averaged model
filter_data$predicted <- predict(avg_model, newdata=filter_data)
filter_data$residuals <- filter_data$median - filter_data$predicted

# spatial autocorrelation test
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("CTI_IB_rec",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

# plot results
plot_CTI_rec <- filter_data %>% 
  ggplot(aes(x=recent_10y, y=median)) +
  geom_hline(yintercept = min_val, linetype="dashed", color="grey", linewidth=.5) +
  geom_pointrange(aes(ymin=q25, ymax=q75, fill=Region, shape=Region), size=.5, stroke = .2, colour="grey50", alpha=.8) +
  geom_line(aes(y=predicted, linetype=Region, color=Region), linewidth =.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0)) +
  xlab(bquote("Recent T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(10.5, 13), "CTI", expand = c(0,0)) +
  mytheme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

## Species richness (qD index)
## Historical data
filter_data <- data_def %>% 
  filter(Survey=="Historical") %>% 
  group_by(Code, Region, historical_10y, qD_SC_filter2, qD.LCL_SC_filter2, qD.UCL_SC_filter2) %>% 
  summarise(n=n()) %>% 
  dplyr::select(-n) %>% 
  mutate(historical_10y_sq = historical_10y^2) %>% 
  filter(Code!="9PIN") %>% 
  rename(qD=qD_SC_filter2)

m1 <- glm(qD ~ scale(historical_10y)*Region + scale(historical_10y_sq),
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get models info
filtered_models <- get_filtered_models_GLM(dredge.model, "qD_filter2_hist")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "qD_filter2_hist")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

# predict values from averaged model
filter_data$predicted <- predict(avg_model, newdata=filter_data)
filter_data$residuals <- filter_data$qD - filter_data$predicted

# spatial autocorrelation test
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("qD_filter2_hist",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

# plot results
filter_data$Region <- factor(levels=c("Gredos","Guadarrama","Meridional","Javalambre"), filter_data$Region)
plot_qD_hist <- filter_data %>% 
  ggplot(aes(x=historical_10y, y=qD, fill=Region)) +
  geom_pointrange(aes(ymin=qD.LCL_SC_filter2, ymax=qD.UCL_SC_filter2, fill=Region, shape=Region),
                  size=.5, stroke = .2, colour="grey50", alpha=.8) +
  geom_line(aes(y=predicted, linetype = Region, color = Region), linewidth=.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0), bquote("Historical T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(0,42), expand = c(0,0), "qD") +
  mytheme

## Recent data 
filter_data <- data_def %>% 
  filter(Survey=="Recent") %>% 
  group_by(Code, Region, recent_10y, qD_SC_filter2, qD.LCL_SC_filter2, qD.UCL_SC_filter2) %>% 
  summarise(n=n()) %>% 
  dplyr::select(-n) %>% 
  mutate(recent_10y_sq = recent_10y^2) %>% 
  filter(Code!="9PIN") %>% 
  rename(qD=qD_SC_filter2)

m1 <- glm(qD ~ scale(recent_10y)*Region + scale(recent_10y_sq),
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get models info
filtered_models <- get_filtered_models_GLM(dredge.model, "qD_filter2_rec")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "qD_filter2_rec")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

# predict values from averaged model
filter_data$predicted <- predict(avg_model, newdata=filter_data)
filter_data$residuals <- filter_data$qD - filter_data$predicted

# spatial autocorrelation test
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("qD_filter2_rec",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

# plot results
filter_data$Region <- factor(levels=c("Gredos","Guadarrama","Meridional","Javalambre"), filter_data$Region)
plot_qD_rec <- filter_data %>% 
  ggplot(aes(x=recent_10y, y=qD, fill=Region)) +
  geom_pointrange(aes(ymin=qD.LCL_SC_filter2, ymax=qD.UCL_SC_filter2, fill=Region, shape=Region),
                  size=.5, stroke = .2, colour="grey50", alpha=.8) +
  geom_line(aes(y=predicted, linetype = Region, color = Region), linewidth=.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0), bquote("Recent T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(0,42), expand = c(0,0), "index qD") +
  mytheme +
  theme(axis.title.y = element_blank())

Fig_4 <- ggarrange(plot_CTI_hist, plot_CTI_rec, plot_qD_hist, plot_qD_rec,
                   ncol = 2, nrow = 2, align = "v", labels = c("A","","B",""))

plot(Fig_4)

## export plot 
tiff(paste0(directory,"/Figures/Fig4_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)
plot(Fig_4)
dev.off()

###### 3.2 Figure S10 European STI vs local temperature ####
## Historical data
filter_data <- data_def %>% 
  filter(variable=="CTI_EU" & Survey=="Historical") %>% 
  dplyr::select(Code, median, q25,q75,historical_10y, Region)

m1 <- glm(median ~ scale(historical_10y)*Region,
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2) # only one model within delta < 2

## get model info
filtered_models <- get_best_model_info(dredge.model, "CTI_EU_hist")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

## get best model coefficients
avg_model_summ <- get_best_model_coefs(subset_models[[1]], "CTI_EU_hist")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

# predict values from averaged model
filter_data$predicted <- predict(subset_models[[1]], newdata=filter_data)
filter_data$residuals <- filter_data$median - filter_data$predicted

# spatial autocorrelation test
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("CTI_EU_hist",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

# plot results
plot_CTI_hist <- filter_data %>% 
  ggplot(aes(x=historical_10y, y=median, fill=Region)) +
  geom_pointrange(aes(ymin=q25, ymax=q75, fill=Region, shape=Region), size=.5, stroke = .2, colour="grey50", alpha=.7) +
  geom_line(aes(y=predicted, linetype=Region, color=Region), linewidth=.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0)) +
  xlab(bquote("Historical T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(8.8, 12.5), "CTI", expand = c(0,0)) +
  mytheme

## Recent data
filter_data <- data_def %>% 
  filter(variable=="CTI_EU" & Survey=="Recent") %>% 
  dplyr::select(Code, median, q25, q75, recent_10y, Region)

m1 <- glm(median ~ scale(recent_10y)*Region,
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2) # only one model within delta < 2

## get model info
filtered_models <- get_best_model_info(dredge.model, "CTI_EU_rec")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

## get best model coefficients
avg_model_summ <- get_best_model_coefs(subset_models[[1]], "CTI_EU_rec")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

# predict values from averaged model
filter_data$predicted <- predict(subset_models[[1]], newdata=filter_data)
filter_data$residuals <- filter_data$median - filter_data$predicted

# spatial autocorrelation test
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("CTI_EU_rec",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

# plot results
plot_CTI_rec <- filter_data %>% 
  ggplot(aes(x=recent_10y, y=median)) +
  geom_pointrange(aes(ymin=q25, ymax=q75, fill=Region, shape=Region), size=.5, stroke = .2, colour="grey50", alpha=.7) +
  geom_line(aes(y=predicted, linetype=Region, color=Region), linewidth =.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0)) +
  xlab(bquote("Recent T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(8.8, 12.5), "CTI", expand = c(0,0)) +
  mytheme +
  theme(axis.title.y = element_blank())

Fig_S10 <- ggarrange(plot_CTI_hist, plot_CTI_rec)
plot(Fig_S10)

tiff(paste0(directory,"/Figures/FigS10_DDI-2024-0174.tif"), 
     width = 20, height = 12.5, units = "cm", res = 300)
plot(Fig_S10)
dev.off()

###### 3.3 Spp richness (filter1 and filter3) vs local temperature (Figure S5) ####
## filter 1 
## Historical data
filter_data <- data_def %>% 
  filter(Survey=="Historical") %>% 
  group_by(Code, Region, historical_10y, qD_SC_filter1, qD.LCL_SC_filter1, qD.UCL_SC_filter1) %>% 
  summarise(n=n()) %>% 
  dplyr::select(-n) %>% 
  mutate(historical_10y_sq = historical_10y^2) %>% 
  rename(qD=qD_SC_filter1)

m1 <- glm(qD ~ scale(historical_10y)*Region + scale(historical_10y_sq),
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get model info
filtered_models <- get_filtered_models_GLM(dredge.model, "qD_filter1_hist")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "qD_filter1_hist")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

filter_data$predicted <- predict(avg_model, newdata=filter_data) # predict values from averaged model

# plot results
filter_data$Region <- factor(levels=c("Gredos","Guadarrama","Meridional","Javalambre"), filter_data$Region)
plot_qD_hist <- filter_data %>% 
  ggplot(aes(x=historical_10y, y=qD, fill=Region)) +
  geom_pointrange(aes(ymin=qD.LCL_SC_filter1, ymax=qD.UCL_SC_filter1, fill=Region, shape=Region),
                  size=.5, stroke = .2, colour="grey50", alpha=.7) +
  geom_line(aes(y=predicted, linetype = Region, color = Region), linewidth=.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0), bquote("Historical T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(0,65), expand = c(0,0), bquote("qD "["SC ≥ 0.66"])) +
  mytheme +
  theme(axis.title.x = element_blank())

## Recent data 
filter_data <- data_def %>% 
  filter(Survey=="Recent") %>% 
  group_by(Code, Region, recent_10y, qD_SC_filter1, qD.LCL_SC_filter1, qD.UCL_SC_filter1) %>% 
  summarise(n=n()) %>% 
  dplyr::select(-n) %>% 
  mutate(recent_10y_sq = recent_10y^2) %>% 
  rename(qD=qD_SC_filter1)

m1 <- glm(qD ~ scale(recent_10y)*Region + scale(recent_10y_sq),
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get model info
filtered_models <- get_filtered_models_GLM(dredge.model, "qD_filter1_rec")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "qD_filter1_rec")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

filter_data$predicted <- predict(avg_model, newdata=filter_data) # predict values from averaged model

# plot results
filter_data$Region <- factor(levels=c("Gredos","Guadarrama","Meridional","Javalambre"), filter_data$Region)
plot_qD_rec <- filter_data %>% 
  ggplot(aes(x=recent_10y, y=qD, fill=Region)) +
  geom_pointrange(aes(ymin=qD.LCL_SC_filter1, ymax=qD.UCL_SC_filter1, fill=Region, shape=Region),
                  size=.5, stroke = .2, colour="grey50", alpha=.7) +
  geom_line(aes(y=predicted, linetype = Region, color = Region), linewidth=.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0), bquote("Recent T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(0,65), expand = c(0,0), bquote("qD "["SC ≥ 0.66"])) +
  mytheme +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())

Fig_S5A <- ggarrange(plot_qD_hist, plot_qD_rec)

## filter 3 
## Historical data
filter_data <- data_def %>% 
  filter(Survey=="Historical") %>% 
  group_by(Code, Region, historical_10y, qD_SC_filter3, qD.LCL_SC_filter3, qD.UCL_SC_filter3) %>% 
  summarise(n=n()) %>% 
  dplyr::select(-n) %>% 
  mutate(historical_10y_sq = historical_10y^2) %>% 
  rename(qD=qD_SC_filter3) %>% 
  na.omit()

m1 <- glm(qD ~ scale(historical_10y)*Region + scale(historical_10y_sq),
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get models info
filtered_models <- get_filtered_models_GLM(dredge.model, "qD_filter3_hist")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "qD_filter3_hist")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

filter_data$predicted <- predict(avg_model, newdata=filter_data) # predict values from averaged model

# plot results
filter_data$Region <- factor(levels=c("Gredos","Guadarrama","Meridional","Javalambre"), filter_data$Region)
plot_qD_hist <- filter_data %>% 
  ggplot(aes(x=historical_10y, y=qD, fill=Region)) +
  geom_pointrange(aes(ymin=qD.LCL_SC_filter3, ymax=qD.UCL_SC_filter3, fill=Region, shape=Region),
                  size=.5, stroke = .2, colour="grey50", alpha=.7) +
  geom_line(aes(y=predicted, linetype = Region, color = Region), linewidth=.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0), bquote("Historical T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(0,65), expand = c(0,0), bquote("qD "["SC ≥ 0.93"])) +
  mytheme

## Recent data 
filter_data <- data_def %>% 
  filter(Survey=="Recent") %>% 
  group_by(Code, Region, recent_10y, qD_SC_filter3, qD.LCL_SC_filter3, qD.UCL_SC_filter3) %>% 
  summarise(n=n()) %>% 
  dplyr::select(-n) %>% 
  mutate(recent_10y_sq = recent_10y^2) %>% 
  rename(qD=qD_SC_filter3) %>% 
  na.omit()

m1 <- glm(qD ~ scale(recent_10y)*Region + scale(recent_10y_sq),
          data = filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get models info
filtered_models <- get_filtered_models_GLM(dredge.model, "qD_filter3_rec")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "qD_filter3_rec")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

filter_data$predicted <- predict(avg_model, newdata=filter_data) # predict values from averaged model

# plot results
filter_data$Region <- factor(levels=c("Gredos","Guadarrama","Meridional","Javalambre"), filter_data$Region)
plot_qD_rec <- filter_data %>% 
  ggplot(aes(x=recent_10y, y=qD, fill=Region)) +
  geom_pointrange(aes(ymin=qD.LCL_SC_filter3, ymax=qD.UCL_SC_filter3, fill=Region, shape=Region),
                  size=.5, stroke = .2, colour="grey50", alpha=.7) +
  geom_line(aes(y=predicted, linetype = Region, color = Region), linewidth=.5) +
  theme_classic() +
  scale_fill_manual(values=c("white", "black", "#3182bd", "red")) +
  scale_color_manual(values=c("black", "black", "#3182bd", "red")) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "solid")) +
  scale_shape_manual(values=c(21,21,24, 24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0), bquote("Recent T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(0,65), expand = c(0,0), bquote("qD "["SC ≥ 0.93"])) +
  mytheme +
  theme(axis.title.y = element_blank())

Fig_S5B <- ggarrange(plot_qD_hist, plot_qD_rec)

Fig_S5 <- ggarrange(Fig_S5A, Fig_S5B, nrow = 2, align="hv", labels = c("A","B"))
plot(Fig_S5)

tiff(paste0(directory, "/Figures/FigS5_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)
plot(Fig_S5)
dev.off()

##### 4. Temporal shift CTI, CTIsd, qD vs Tmax ####
###### 4.1 Change in Iberian CTI vs Tmax (Figure 5) ####
filter_data <- data_def %>% 
  filter(Survey=="diff" & variable=="CTI_IB") %>% 
  dplyr::select(Code, Region, median, q25, q75, historical_10y, diff_10y)

m1 <- glm(median ~ scale(diff_10y)*scale(historical_10y)*Region,
          data=filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2) # only one model within delta < 2

## get best model info
filtered_models <- get_best_model_info(dredge.model, "ΔCTI_IB")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

# get best model coefficients
avg_model_summ <- get_best_model_coefs(subset_models[[1]], "ΔCTI_IB")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

# predict values from averaged model
filter_data$predicted <- predict(subset_models[[1]], newdata=filter_data)
filter_data$residuals <- filter_data$median - filter_data$predicted

# spatial autocorrelation test
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("ΔCTI_IB",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

## plot results
gre <- filter_data %>% 
  filter(Region=="Gredos") %>% 
  ggplot(aes(x=historical_10y, y=median)) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size=.3) +
  geom_pointrange(aes(ymin=q25, ymax=q75, fill=Region, shape=Region), size=.5, stroke = .2, colour="grey50") +
  geom_line(aes(y=predicted, linetype=Region, color=Region), linewidth=.5) +
  theme_classic() +
  ggtitle("Sierra de Gredos") +
  scale_fill_manual(values=c("white")) +
  scale_color_manual(values=c("black")) +
  scale_linetype_manual(values=c("dashed")) +
  scale_shape_manual(values=c(21)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0)) +
  xlab(bquote("Historical T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(-1, 1.5), "Δ CTI") +
  mytheme +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=14, face="bold"))

gua <- filter_data %>% 
  filter(Region=="Guadarrama") %>% 
  ggplot(aes(x=historical_10y, y=median)) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size=.3) +
  geom_pointrange(aes(ymin=q25, ymax=q75, fill=Region, shape=Region), size=.5, stroke = .2, colour="grey50") +
  geom_line(aes(y=predicted, linetype=Region, color=Region), linewidth=.5) +
  theme_classic() +
  ggtitle("Sierra de Guadarrama") +
  scale_fill_manual(values=c("black")) +
  scale_color_manual(values=c("black")) +
  scale_linetype_manual(values=c("solid")) +
  scale_shape_manual(values=c(21)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0)) +
  xlab(bquote("Historical T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(-1, 1.5), "Δ CTI") +
  mytheme +
  theme(plot.title = element_text(size=14, face="bold"))

mer <- filter_data %>% 
  filter(Region=="Meridional") %>% 
  ggplot(aes(x=historical_10y, y=median)) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size=.3) +
  geom_pointrange(aes(ymin=q25, ymax=q75, fill=Region, shape=Region), size=.5, stroke = .2, colour="grey50") +
  geom_line(aes(y=predicted, linetype=Region, color=Region), linewidth=.5) +
  theme_classic() +
  ggtitle("Sistema Ibérico Meridional") +
  scale_fill_manual(values=c("#3182bd")) +
  scale_color_manual(values=c("#3182bd")) +
  scale_linetype_manual(values=c("solid")) +
  scale_shape_manual(values=c(24)) +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0)) +
  xlab(bquote("Historical T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(-1, 1.5), "Δ CTI") +
  mytheme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=14, face="bold"))

jav <- filter_data %>% 
  filter(Region=="Javalambre") %>% 
  ggplot(aes(x=historical_10y, y=median)) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size=.3) +
  geom_pointrange(aes(ymin=q25, ymax=q75, fill=Region, shape=Region), size=.5, stroke = .2, colour="grey50") +
  geom_line(aes(y=predicted, linetype=Region, color=Region), linewidth=.5) +
  theme_classic() +
  ggtitle("Javalambre") +
  scale_x_continuous(limits= c(15.5,28), expand=c(0,0)) +
  scale_fill_manual(values=c("red")) +
  scale_color_manual(values=c("red")) +
  scale_linetype_manual(values=c("solid")) +
  scale_shape_manual(values=c(24)) +
  xlab(bquote("Historical T"["max"]*" (ºC)")) +
  scale_y_continuous(limits=c(-1, 1.5), "Δ CTI") +
  mytheme +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(size=14, face="bold"))

## export plot
Fig_5 <- ggarrange(gre, mer, gua, jav, ncol=2, nrow=2, align="v")
plot(Fig_5)

tiff(paste0(directory, "/Figures/Fig5_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)
plot(Fig_5)
dev.off()

rm(gre,gua,mer,jav)

###### 4.2 Change in European CTI vs Tmax (Figure S6) #####
filter_data <- data_def %>% 
  filter(Survey=="diff" & variable=="CTI_EU") %>% 
  dplyr::select(Code, Region, median, q25, q75, historical_10y, diff_10y)

m1 <- glm(median ~ scale(diff_10y)*scale(historical_10y)*Region,
          data=filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get models info
filtered_models <- get_filtered_models_GLM(dredge.model, "ΔCTI_EU")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "ΔCTI_EU")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

## get predictions for each region
x <-range(filter_data$diff_10y)
x <-seq(x[1],x[2], length.out=50)
y <-range(filter_data$historical_10y)
y <-seq(y[1],y[2], length.out=50)

regions <- c("Gredos","Guadarrama","Meridional","Javalambre")

# Set desired z-axis limits
z_min <- -1  # minimum limit for z-axis
z_max <- 1   # maximum limit for z-axis

# Generate predictions and clip z values
predictions <- lapply(regions, function(region) {
  z <- outer(
    x, y,
    function(diff_10y, historical_10y) {
      predict(avg_model, data.frame(
        diff_10y = diff_10y,
        historical_10y = historical_10y,
        Region = region
      ))
    }
  )
  
  # Clip the z values to be within the specified range
  z_clipped <- pmin(pmax(z, z_min), z_max)
  
  return(z_clipped)
})

names(predictions) <- regions

# Function to create the plot for a given region
create_plot3D <- function(region, x, y, z_min, z_max) {
  z_clipped <- predictions[[region]]
  
  col.pal <- colorRampPalette(c("#4575b4","#91bfdb", "#e0f3f8","#fee090", "#fc8d59", "#d73027"))
  colors <- col.pal(100)
  z.facet.center <- (z_clipped[-1, -1] + z_clipped[-1, -ncol(z_clipped)] + 
                       z_clipped[-nrow(z_clipped), -1] + z_clipped[-nrow(z_clipped), -ncol(z_clipped)])/4
  z.facet.range <- cut(z.facet.center, 100)
  
  persp(x, y, z_clipped, theta=223, phi=23,
        xlab="shift_Annual_Tmax", ylab="historical_Tmax", zlab="shift_CTI_EU",
        zlim = c(z_min, z_max),
        col=colors[z.facet.range], expand = 0.8, r=2, shade=0.4, border=NA,
        ticktype = "detailed",
        main = region)
}

# Open PDF file
pdf(paste0(directory,"/Figures/FigS6_DDI-2024-0174.pdf"), width = 8, height = 6)

# List of regions
regions <- c("Gredos", "Guadarrama", "Meridional", "Javalambre")

for (region in regions) {
  create_plot3D(region, x, y, z_min, z_max)
}

dev.off()

###### 4.3 Change in CTIsd vs Tmax ####
## Change in Iberian CTIsd
filter_data <- data_def %>% 
  filter(Survey=="diff" & variable=="CTI_sd_IB") %>% 
  dplyr::select(Code, Region, median, q25, q75, historical_10y, diff_10y)

m1 <- glm(median ~ scale(diff_10y)*scale(historical_10y)*Region,
          data=filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2) # only one model within delta < 2

## get best model info
filtered_models <- get_best_model_info(dredge.model, "ΔCTI_sd_IB")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

# get best model coefficients
avg_model_summ <- get_best_model_coefs(subset_models[[1]], "ΔCTI_sd_IB")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

# predict values from averaged model
filter_data$predicted <- predict(subset_models[[1]], newdata=filter_data)
filter_data$residuals <- filter_data$median - filter_data$predicted

# Spatial autocorrelation test on model residuals
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("ΔCTI_sd_IB",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

## Change in European CTIsd
filter_data <- data_def %>% 
  filter(Survey=="diff" & variable=="CTI_sd_EU") %>% 
  dplyr::select(Code, Region, median, q25, q75, historical_10y, diff_10y)

m1 <- glm(median ~ scale(diff_10y)*scale(historical_10y)*Region,
          data=filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2) 

## get best model info
filtered_models <- get_filtered_models_GLM(dredge.model, "ΔCTI_sd_EU")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

# get best model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "ΔCTI_sd_EU")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

# predict values from averaged model
filter_data$predicted <- predict(subset_models[[1]], newdata=filter_data)
filter_data$residuals <- filter_data$median - filter_data$predicted

# Spatial autocorrelation test on model residuals
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("ΔCTI_sd_EU",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

###### 4.4 Change in qD filter 2 vs Tmax ####
filter_data <- data_def %>% 
  filter(Survey=="diff" & variable=="CTI_IB") %>% 
  dplyr::select(Code, Region, qD_SC_filter2, historical_10y, diff_10y) %>% 
  na.omit()

aren <- data_def %>% 
  filter(Code=="6Aren") %>% 
  dplyr::select(Code, Survey, Region, qD_SC_filter2, historical_10y, diff_10y) %>% 
  pivot_wider(names_from = Survey, values_from = qD_SC_filter2) %>% 
  mutate(qD_SC_filter2 = Recent - Historical) %>% 
  dplyr::select(Code, Region, qD_SC_filter2, historical_10y, diff_10y)

filter_data <- rbind(filter_data, aren)
rm(aren)

m1 <- glm(qD_SC_filter2 ~ scale(diff_10y)*scale(historical_10y)*Region,
          data=filter_data,
          na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2"))
subset_models <- get.models(dredge.model, subset = delta < 2)

## get models info
filtered_models <- get_filtered_models_GLM(dredge.model, "ΔqD_filter2")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # append to main table

avg_model <- model.avg(dredge.model, subset = delta < 2, fit = TRUE)   # Calculate averaged model

## get averaged model coefficients
avg_model_summ <- get_avg_model_coefs(avg_model, "ΔqD_filter2")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # append to main table

filter_data$predicted <- predict(avg_model, newdata=filter_data) # predict values from averaged model
filter_data$residuals <- filter_data$qD_SC_filter2 - filter_data$predicted

# Spatial autocorrelation test on model residuals
MC <- moran.mc(filter_data$residuals,lw, alternative="greater", nsim = 999)
vec <- data.frame("ΔqD_filter2",MC$statistic, MC$parameter, MC$p.value)
colnames(vec) <- colnames(A_spatial_autocorr_residuals)
A_spatial_autocorr_residuals <- rbind(A_spatial_autocorr_residuals, vec)

## get predictions for each region
x <-range(filter_data$diff_10y)
x <-seq(x[1],x[2], length.out=50)
y <-range(filter_data$historical_10y)
y <-seq(y[1],y[2], length.out=50)

regions <- c("Gredos","Guadarrama","Meridional","Javalambre")

# Set desired z-axis limits
z_min <- -20  # minimum limit for z-axis
z_max <- 20   # maximum limit for z-axis

# Generate predictions and clip z values
predictions <- lapply(regions, function(region) {
  z <- outer(
    x, y,
    function(diff_10y, historical_10y) {
      predict(avg_model, data.frame(
        diff_10y = diff_10y,
        historical_10y = historical_10y,
        Region = region
      ))
    }
  )
  
  # Clip the z values to be within the specified range
  z_clipped <- pmin(pmax(z, z_min), z_max)
  
  return(z_clipped)
})

names(predictions) <- regions

# Function to create the plot for a given region
create_plot3D <- function(region, x, y, z_min, z_max) {
  z_clipped <- predictions[[region]]
  
  col.pal <- colorRampPalette(c("#4575b4","#91bfdb", "#e0f3f8","#fee090", "#fc8d59", "#d73027"))
  colors <- col.pal(100)
  z.facet.center <- (z_clipped[-1, -1] + z_clipped[-1, -ncol(z_clipped)] + 
                       z_clipped[-nrow(z_clipped), -1] + z_clipped[-nrow(z_clipped), -ncol(z_clipped)])/4
  z.facet.range <- cut(z.facet.center, 100)
  
  persp(x, y, z_clipped, theta=223, phi=23,
        xlab="shift_Annual_Tmax", ylab="historical_Tmax", zlab="shift_qD",
        zlim = c(z_min, z_max),
        col=colors[z.facet.range], expand = 0.8, r=2, shade=0.4, border=NA,
        ticktype = "detailed",
        main = region)
}

## export plots
pdf(paste0(directory,"/Figures/FigS7_DDI-2024-0174.pdf"), width = 8, height = 6) # Open PDF file

regions <- c("Gredos", "Guadarrama", "Meridional", "Javalambre") # List of regions

for (region in regions) {
  create_plot3D(region, x, y, z_min, z_max) # Loop through each region and create plot
}

dev.off()

##### 5. Persistence analysis ####
directory <- rstudioapi::getActiveDocumentContext()$path
directory <- dirname(directory)

pers_data <-read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "4. persistence_cat_analysis")

###### 5.1 Iberian STI (Figure 6A) ####
filter_data <- pers_data %>% 
  filter(variable=="CTI_IB") %>% 
  dplyr::select(Code, Region, persistence_cat, median, historical_10y)

m1 <- lmer(median ~ scale(historical_10y)*persistence_cat + (1|Region/Code),
           data = filter_data,
           REML = F,
           na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2")) 
dredge.model
subset_models <- get.models(dredge.model, subset = delta < 2) # only one model within delta < 2

# get model info from LMER
filtered_models <- get_filtered_models_LMER(dredge.model, "CTI_IB_persistence_cat")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # save to table to export

# get model coefficients from LMER
coefs <- get_coefficients_LMER(subset_models[[1]], "CTI_IB_persistence_cat")
A_coefs_export <- rbind(A_coefs_export, coefs) # save to table to export

## contrast analysis on estimated marginal means
em <- emmeans(m1, "persistence_cat")
Table_S9 <- data.frame(contrast(em, "pairwise", adjust = "Tukey")) %>% 
  mutate(variable="CTI_IB") %>% 
  relocate(variable, .before = contrast)

## predict values using best model
pred <- predictSE(m1, newdata=filter_data) # predict values adding SE
filter_data$fit = pred$fit
filter_data$se.fit = pred$se.fit
filter_data$residuals = filter_data$median - filter_data$fit

filter_data$persistence_cat <- gsub("colonize","colonise",filter_data$persistence_cat)
filter_data$persistence_cat <- factor(filter_data$persistence_cat, levels = c("colonise", "persist", "extinct"))

## plot results
Fig_6A <- filter_data %>% 
  ggplot(aes(x=historical_10y, y=fit)) +
  geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit, fill=persistence_cat), alpha=.7, linetype = 0) +
  geom_line(aes(linetype=persistence_cat), size=0.5) +
  scale_fill_manual(values = c("grey10", "grey50", "grey90"), name = "Spp persistence") +
  scale_linetype_manual(values = c("solid", "dotdash", "longdash"), name = "Spp persistence") +
  scale_y_continuous(limits = c(10.7,13), "CTI", expand = c(0,0)) +
  scale_x_continuous(limits = c(16.5, 27.5), breaks=c(17,20,23,26)) +
  xlab(bquote("Historical T"["max"]*" (ºC)")) +
  theme_classic() +
  mytheme +
  theme(legend.position="top",
        axis.title.x = element_text(vjust = -2),
        axis.title.y = element_text(vjust = 2.5))

plot(Fig_6A)

tiff(paste0(directory, "/Figures/Fig6A_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)
plot(Fig_6A)
dev.off()

## spatial autocorrelation test 
append_spatial_autocorr <- function(table, model_name, statistic, obs_rank, p_value) {
  new_row <- data.frame(
    model = model_name,
    statistic = statistic,
    Obs_rank = obs_rank,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
  return(rbind(table, new_row))
}

temp <- filter_data %>% 
  filter(persistence_cat=="colonise")
MC<- moran.mc(temp$residuals,lw, alternative="greater", nsim = 999) 
A_spatial_autocorr_residuals <- append_spatial_autocorr(A_spatial_autocorr_residuals, "CTI_IB_colonise",
                                                        MC$statistic, MC$parameter, MC$p.value)

temp <- filter_data %>% 
  filter(persistence_cat=="persist")
MC<- moran.mc(temp$residuals,lw, alternative="greater", nsim = 999) 
A_spatial_autocorr_residuals <- append_spatial_autocorr(A_spatial_autocorr_residuals, "CTI_IB_persist",
                                                        MC$statistic, MC$parameter, MC$p.value)

temp <- filter_data %>% 
  filter(persistence_cat=="extinct")
MC<- moran.mc(temp$residuals,lw, alternative="greater", nsim = 999) 
A_spatial_autocorr_residuals <- append_spatial_autocorr(A_spatial_autocorr_residuals, "CTI_IB_extinct",
                                                        MC$statistic, MC$parameter, MC$p.value)

###### 5.2 European STI (Figure S11A) ####
filter_data <- pers_data %>% 
  filter(variable=="CTI_EU") %>% 
  dplyr::select(Code, Region, persistence_cat, median, historical_10y)

m1 <- lmer(median ~ scale(historical_10y)*persistence_cat + (1|Region/Code),
           data = filter_data,
           REML = F,
           na.action = "na.fail")

dredge.model <- dredge(m1, trace=F, extra = c("R^2", "adjR^2")) 
dredge.model
subset_models <- get.models(dredge.model, subset = delta < 2) # only one model within delta < 2

# get model info from LMER
filtered_models <- get_filtered_models_LMER(dredge.model, "CTI_EU_persistence_cat")
A_model_outputs_exp <- rbind(A_model_outputs_exp, filtered_models) # save to table to export

# get model coefficients from LMER
coefs <- get_coefficients_LMER(subset_models[[1]], "CTI_EU_persistence_cat")
A_coefs_export <- rbind(A_coefs_export, coefs) # save to table to export

## contrast analysis on estimated marginal means
em <- emmeans(m1, "persistence_cat")
cont_em <- data.frame(contrast(em, "pairwise", adjust = "Tukey")) %>% 
  mutate(variable="CTI_EU") %>% 
  relocate(variable, .before = contrast)

Table_S9 <- rbind(Table_S9, cont_em)

## export contrast analysis table
write_xlsx(Table_S9, paste0(directory,"/Tables/TableS9_DDI-2024-0174.xlsx"))

## predict values from best model
pred <- predictSE(m1, newdata=filter_data) # predict values adding SE
filter_data$fit = pred$fit
filter_data$se.fit = pred$se.fit
filter_data$residuals = filter_data$median - filter_data$fit

filter_data$persistence_cat <- gsub("colonize","colonise",filter_data$persistence_cat)
filter_data$persistence_cat <- factor(filter_data$persistence_cat, levels = c("colonise", "persist", "extinct"))

## plot results
Fig_S11A <- filter_data %>% 
  ggplot(aes(x=historical_10y, y=fit)) +
  geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit, fill=persistence_cat), alpha=.7, linetype = 0) +
  geom_line(aes(linetype=persistence_cat), size=0.5) +
  scale_fill_manual(values = c("grey10", "grey50", "grey90"), name = "persistence_cat") +
  scale_linetype_manual(values = c("solid", "dotdash", "longdash"), name = "persistence_cat") +
  scale_y_continuous(limits = c(9,12), "CTI", expand = c(0,0)) +
  scale_x_continuous(limits = c(16.5, 27.5), breaks=c(17,20,23,26)) +
  xlab(bquote("Historical T"["max"]*" (ºC)")) +
  theme_classic() +
  mytheme +
  theme(legend.position="top",
        axis.title.x = element_text(vjust = -2),
        axis.title.y = element_text(vjust = 2.5))

plot(Fig_S11A)

tiff(paste0(directory, "/Figures/FigS11A_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)
plot(Fig_S11A)
dev.off()

## spatial autocorrelation test 
temp <- filter_data %>% 
  filter(persistence_cat=="colonise")
MC<- moran.mc(temp$residuals,lw,alternative="greater", nsim = 999) 
A_spatial_autocorr_residuals <- append_spatial_autocorr(A_spatial_autocorr_residuals, "CTI_EU_colonise",
                                                        MC$statistic, MC$parameter, MC$p.value)

temp <- filter_data %>% 
  filter(persistence_cat=="persist")
MC<- moran.mc(temp$residuals,lw, alternative="greater", nsim = 999) 
A_spatial_autocorr_residuals <- append_spatial_autocorr(A_spatial_autocorr_residuals, "CTI_EU_persist",
                                                        MC$statistic, MC$parameter, MC$p.value)

temp <- filter_data %>% 
  filter(persistence_cat=="extinct")
MC<- moran.mc(temp$residuals,lw, alternative="greater", nsim = 999) 
A_spatial_autocorr_residuals <- append_spatial_autocorr(A_spatial_autocorr_residuals, "CTI_EU_extinct",
                                                        MC$statistic, MC$parameter, MC$p.value)

##### export spatial autocorrelation table (Table S3) ####
print(A_spatial_autocorr_residuals)
write_xlsx(A_spatial_autocorr_residuals, paste0(directory,"/Tables/TableS3_DDI-2024-0174.xlsx"))

##### 6. Shift in occupancy vs spp thermal affinities (PGLS models) ####
directory <- rstudioapi::getActiveDocumentContext()$path
directory <- dirname(directory)

n_sites_diff <-read_xlsx(paste0(directory, "/supp_material_Ursul_etal_2025_DDI.xlsx"), "5. Occupancy_analysis")

## phylogenetic tree adapted from Wiemers et al (2019) - https://zenodo.org/records/3531555
# Wiemers, M., Chazot, N., Wheat, C. W., Schweiger, O., & Wahlberg, N. (2019). A complete time-calibrated multi-gene phylogeny of the European butterflies - data [Data set]

if (!require("pacman")) install.packages("pacman")
pacman::p_load(aod, ape, geiger, caper, phytools)

arbol.cort <- read.tree(paste0(directory,"/Phylogenetic_tree_butterflies_cropped.nwk"))

## export phylogenetic tree figure (Figure S3)
png(filename=paste0(directory, "/Figures/FigS3_DDI-2024-0174.png"),
    width     = 12,
    height    = 22,
    units     = "in",
    res       = 500)

plot(arbol.cort, cex=1, label.offset = 2, no.margin=T)
dev.off()

## prepare data for pgls model
filter_data <- n_sites_diff
filter_data$IBER_mean_temp_sq <- filter_data$IBER_mean_temp^2 # squared term for Iberian STI
filter_data$EU_temp.mean_sq <- filter_data$EU_temp.mean^2 # squared term for European STI

kk <- data.frame(arbol.cort$tip.label)
colnames(kk) <- "Name"

filter_data <- left_join(kk, filter_data, by="Name") ## ROWS IN OUR DATASET WITH SAME ORDER AS TREE!!

cdat <- comparative.data(arbol.cort, filter_data, Name, vcv=TRUE, vcv.dim=3) ## this is what we will use in the pgls function

## create function to run pgls models
pgls_model_function <- function(i, num_list, formulas, cdat, models_list, models_output) {
  print(i)
  model_formula <- as.formula(formulas[[i]])
  model_name <- paste0("model_", i)
  model <- pgls(model_formula, data = cdat, lambda = "ML")
  models_list[[model_name]] <- model

  # Extract model statistics
  r_squared <- summary(model)$r.squared
  n <- nrow(cdat$data)
  p <- length(model$model$coef) - 1
  adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))
  lambda_value <- model$param["lambda"]
  
  vec <- c(model_name, formulas[[i]], r_squared, adj_r_squared, summary(model)$df[1], model$model$log.lik, model$aicc, lambda_value)
  vec_t <- data.frame(matrix(ncol=8, nrow=1))
  colnames(vec_t) <- colnames(models_output)
  vec_t[1,] <- vec
  
  models_output <- rbind(models_output, vec_t)
  return(list(models_list = models_list, models_output = models_output))
}

## function to extract coefficients from averaged models after PGLS
pgls_avg_models_coefs <- function(model.avg, response_variable_name) {
  avg_model_summ <- summary(model_avg)
  avg_model_summ <- data.frame(avg_model_summ$coefmat.full)
  avg_model_summ <- avg_model_summ %>% 
    mutate(response_variable = response_variable_name) # Add response variable
  avg_model_summ$selected_variables <- rownames(avg_model_summ) # Add selected variables
  rownames(avg_model_summ) <- NULL # Remove rownames
  avg_model_summ <- avg_model_summ %>% 
    relocate(response_variable:selected_variables, .before = Estimate) # Reorder columns
  colnames(avg_model_summ) <- colnames(A_coefs_export)
  return(avg_model_summ)
}

###### 6.1 Iberian STI (Figure 6B) ####
# Define variables
vars <- c("IBER_mean_temp", "IBER_mean_temp_sq", "IBER_sd_temp")

# Get all possible combinations
combinations <- expand.grid(
  IBER_mean_temp = c(0, 1),
  IBER_mean_temp_sq = c(0, 1),
  IBER_mean_temp_IBER_sd_temp = c(0, 1),  # Add interaction term
  IBER_sd_temp = c(0, 1))

# Convert combinations to formulas
formulas <- apply(combinations, 1, function(x) {
  terms <- c()
  if (x["IBER_mean_temp"] == 1) terms <- c(terms, "IBER_mean_temp")
  if (x["IBER_mean_temp_sq"] == 1) terms <- c(terms, "IBER_mean_temp_sq")
  if (x["IBER_mean_temp_IBER_sd_temp"] == 1) terms <- c(terms, "IBER_mean_temp:IBER_sd_temp")  # Add interaction term
  if (x["IBER_sd_temp"] == 1) terms <- c(terms, "IBER_sd_temp")
  formula <- paste("median ~", paste(terms, collapse = " + "))
  return(formula)
})

formulas[[1]] <- "median ~ 1"
formulas <- formulas[-c(3,7,11,15)] # remove all formulas with IBER_mean_temp_sq and without IBER_mean_temp
formulas # Print the formulas

## Run pgls
models_list <- list()
models_output <- data.frame(matrix(ncol=8, nrow=0))
colnames(models_output) <- c("m_num","formula","R2","Adj_R2","df","LogLik","AICc","lambda")
num_list <- 1:12

for (i in num_list) {
  result <- pgls_model_function(i, num_list, formulas, cdat, models_list, models_output)
  models_list <- result$models_list
  models_output <- result$models_output
}

models_output$AICc <- as.numeric(as.character(models_output$AICc))
min_AICc <- min(models_output$AICc) # Find the minimum AICc
best_models <- models_output[models_output$AICc <= (min_AICc + 2), ] # Select models with AICc less than 2 from the best model

best_models <- best_models %>% 
  arrange(AICc)

best_models # print best models
best_model_names <- best_models$m_num # Extract model names from the best models

best_models$m_num <- "Change_Nsites_vs_STI_IB" 
A_model_outputs_Nsites <- best_models %>% 
  rename(response_variable=m_num)

best_model_list <- lapply(best_model_names, function(x) models_list[[x]]) # Extract the actual models from models_list
model_avg <- model.avg(best_model_list) # Perform model averaging

## get averaged model coefficients
avg_model_summ <- pgls_avg_models_coefs(best_model_list, "Change_Nsites_vs_STI_IB")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # save to table to export

## plot results
x <- range(filter_data$IBER_mean_temp)
x <- seq(x[1], x[2], length.out = 50)

y <- range(filter_data$IBER_sd_temp)
y <- seq(y[1], y[2], length.out = 50)

grid <- expand.grid(IBER_mean_temp = x, IBER_sd_temp = y) # Create a grid of x (STImean) and y (STIsd)
grid$IBER_mean_temp_sq <- grid$IBER_mean_temp^2 # Add the quadratic term for STImean to the grid

predictions <- predict(model_avg, newdata = grid) # Predict using the averaged model

z <- matrix(predictions, nrow = length(x), ncol = length(y)) # Reshape predictions into a matrix before persp plot

# Set desired z-axis limits
z_min <- -10  # minimum limit for z-axis
z_max <- 20   # maximum limit for z-axis

z_clipped <- pmin(pmax(z, z_min), z_max) # Clip the z values to be within the specified range

# Set palette for plot
col.pal <- colorRampPalette(c("#4575b4","#91bfdb", "#e0f3f8","#fee090", "#fc8d59", "#d73027"))
colors <- col.pal(100)

# Calculate facet center for coloring
z.facet.center <- (z_clipped[-1, -1] + z_clipped[-1, -ncol(z_clipped)] + z_clipped[-nrow(z_clipped), -1] + z_clipped[-nrow(z_clipped), -ncol(z_clipped)]) / 4
z.facet.range <- cut(z.facet.center, 100)  # Range of the facet center on a 100-scale (number of colors)

## export plot
tiff(paste0(directory, "/Figures/Fig6B_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)

persp(x,y,z_clipped, theta=140, phi=25,
      xlab="STI_mean", ylab="STI_sd", zlab="Nsites_diff",
      zlim = c(z_min, z_max),
      col=colors[z.facet.range], expand = .8, r=2, shade=.4, border=NA,
      ticktype = "detailed") 

dev.off()

###### 6.2 European STI (Figure S11B) ####
# Define variables
vars <- c("EU_temp.mean", "EU_temp.mean_sq", "EU_temp.sd")

# Get all possible combinations
combinations <- expand.grid(
  EU_temp.mean = c(0, 1),
  EU_temp.mean_sq = c(0, 1),
  EU_temp.mean_EU_temp.sd = c(0, 1),  # Add interaction term
  EU_temp.sd = c(0, 1))

# Convert combinations to formulas
formulas <- apply(combinations, 1, function(x) {
  terms <- c()
  if (x["EU_temp.mean"] == 1) terms <- c(terms, "EU_temp.mean")
  if (x["EU_temp.mean_sq"] == 1) terms <- c(terms, "EU_temp.mean_sq")
  if (x["EU_temp.mean_EU_temp.sd"] == 1) terms <- c(terms, "EU_temp.mean:EU_temp.sd")  # Add interaction term
  if (x["EU_temp.sd"] == 1) terms <- c(terms, "EU_temp.sd")
  formula <- paste("median ~", paste(terms, collapse = " + "))
  return(formula)
})

formulas[[1]] <- "median ~ 1"
formulas <- formulas[-c(3,7,11,15)] # remove all formulas with mean_temp_sq and without mean_temp
formulas # Print the formulas

## Run pgls
models_list <- list()
models_output <- data.frame(matrix(ncol=8, nrow=0))
colnames(models_output) <- c("m_num","formula","R2","Adj_R2","df","LogLik","AICc","lambda")
num_list <- 1:12

for (i in num_list) {
  result <- pgls_model_function(i, num_list, formulas, cdat, models_list, models_output)
  models_list <- result$models_list
  models_output <- result$models_output
}

models_output$AICc <- as.numeric(as.character(models_output$AICc))
min_AICc <- min(models_output$AICc) # Find the minimum AICc
best_models <- models_output[models_output$AICc <= (min_AICc + 2), ] # Select models with AICc less than 2 from the best model

best_models <- best_models %>% 
  arrange(AICc)

best_models # print best models
best_model_names <- best_models$m_num # Extract model names from the best models

best_models$m_num <- "Change_Nsites_vs_STI_EU" 
temp <- best_models %>% 
  rename(response_variable=m_num)
A_model_outputs_Nsites <- rbind(A_model_outputs_Nsites, temp) # append to main table

best_model_list <- lapply(best_model_names, function(x) models_list[[x]]) # Extract the actual models from models_list
model_avg <- model.avg(best_model_list) # Perform model averaging

## get averaged model coefficients
avg_model_summ <- pgls_avg_models_coefs(best_model_list, "Change_Nsites_vs_STI_EU")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # save to table to export

# Set ranges for 3D graph
x <- range(filter_data$EU_temp.mean)
x <- seq(x[1], x[2], length.out = 50)

y <- range(filter_data$EU_temp.sd)
y <- seq(y[1], y[2], length.out = 50)

# Create a grid of x (STImean) and y (STIsd)
grid <- expand.grid(EU_temp.mean = x, EU_temp.sd = y)

# Add the quadratic term for STImean to the grid
grid$EU_temp.mean_sq <- grid$EU_temp.mean^2

# Predict using the averaged model
predictions <- predict(model_avg, newdata = grid)

# Reshape predictions into a matrix before persp plot
z <- matrix(predictions, nrow = length(x), ncol = length(y))

# Set desired z-axis limits
z_min <- -10  # minimum limit for z-axis
z_max <- 20   # maximum limit for z-axis

# Clip the z values to be within the specified range
z_clipped <- pmin(pmax(z, z_min), z_max)

# Set palette for plot
col.pal <- colorRampPalette(c("#4575b4","#91bfdb", "#e0f3f8","#fee090", "#fc8d59", "#d73027"))
colors <- col.pal(100)

# Calculate facet center for coloring
z.facet.center <- (z_clipped[-1, -1] + z_clipped[-1, -ncol(z_clipped)] + z_clipped[-nrow(z_clipped), -1] + z_clipped[-nrow(z_clipped), -ncol(z_clipped)]) / 4
z.facet.range <- cut(z.facet.center, 100)  # Range of the facet center on a 100-scale (number of colors)

## export plot
tiff(paste0(directory, "/Figures/FigS11B_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)

persp(x,y,z_clipped, theta=140, phi=25,
      xlab="STI_mean", ylab="STI_sd", zlab="Nsites_diff",
      zlim = c(z_min, z_max),
      col=colors[z.facet.range], expand = .8, r=2, shade=.4, border=NA,
      ticktype = "detailed") 

dev.off()

##### 7. Historical occupancy vs spp thermal affinities (PGLS models) ####
###### 7.1 Iberian STI (Figure S8A) ####
# Define the variables
vars <- c("IBER_mean_temp", "IBER_mean_temp_sq", "IBER_sd_temp")

# Get all possible combinations
combinations <- expand.grid(
  IBER_mean_temp = c(0, 1),
  IBER_mean_temp_sq = c(0, 1),
  IBER_mean_temp_IBER_sd_temp = c(0, 1),  # Add interaction term
  IBER_sd_temp = c(0, 1))

# Convert combinations to formulas
formulas <- apply(combinations, 1, function(x) {
  terms <- c()
  if (x["IBER_mean_temp"] == 1) terms <- c(terms, "IBER_mean_temp")
  if (x["IBER_mean_temp_sq"] == 1) terms <- c(terms, "IBER_mean_temp_sq")
  if (x["IBER_mean_temp_IBER_sd_temp"] == 1) terms <- c(terms, "IBER_mean_temp:IBER_sd_temp")  # Add interaction term
  if (x["IBER_sd_temp"] == 1) terms <- c(terms, "IBER_sd_temp")
  formula <- paste("N_sites_hist ~", paste(terms, collapse = " + "))
  return(formula)
})

formulas[[1]] <- "N_sites_hist ~ 1"
formulas <- formulas[-c(3,7,11,15)] # remove all formulas with mean_temp_sq and without IBER_mean_temp
formulas # Print the formulas

## Run pgls
models_list <- list()
models_output <- data.frame(matrix(ncol=8, nrow=0))
colnames(models_output) <- c("m_num","formula","R2","Adj_R2","df","LogLik","AICc","lambda")
num_list <- 1:12

for (i in num_list) {
  result <- pgls_model_function(i, num_list, formulas, cdat, models_list, models_output)
  models_list <- result$models_list
  models_output <- result$models_output
}

models_output$AICc <- as.numeric(as.character(models_output$AICc))
min_AICc <- min(models_output$AICc) # Find the minimum AICc
best_models <- models_output[models_output$AICc <= (min_AICc + 2), ] # Select models with AICc less than 2 from the best model

best_models <- best_models %>% 
  arrange(AICc)

best_models # print best models
best_model_names <- best_models$m_num # Extract model names from the best models

best_models$m_num <- "Nsites_hist_vs_STI_IB" 
temp <- best_models %>% 
  rename(response_variable=m_num)
A_model_outputs_Nsites <- rbind(A_model_outputs_Nsites, temp) # append to main table

best_model_list <- lapply(best_model_names, function(x) models_list[[x]]) # Extract the actual models from models_list
model_avg <- model.avg(best_model_list) # Perform model averaging

## get averaged model coefficients
avg_model_summ <- pgls_avg_models_coefs(best_model_list, "Nsites_hist_vs_STI_IB")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # save to table to export

# Set ranges for 3D graph
x <- range(filter_data$IBER_mean_temp)
x <- seq(x[1], x[2], length.out = 50)

y <- range(filter_data$IBER_sd_temp)
y <- seq(y[1], y[2], length.out = 50)

grid <- expand.grid(IBER_mean_temp = x, IBER_sd_temp = y) # Create a grid of x (STImean) and y (STIsd)
grid$IBER_mean_temp_sq <- grid$IBER_mean_temp^2 # Add the quadratic term for STImean to the grid

predictions <- predict(model_avg, newdata = grid) # Predict using the averaged model

z <- matrix(predictions, nrow = length(x), ncol = length(y)) # Reshape predictions into a matrix before persp plot

# Set desired z-axis limits
z_min <- 0  # minimum limit for z-axis
z_max <- 40   # maximum limit for z-axis

# Clip the z values to be within the specified range
z_clipped <- pmin(pmax(z, z_min), z_max)

# Set palette for plot
col.pal <- colorRampPalette(c("#4575b4","#91bfdb", "#e0f3f8","#fee090", "#fc8d59", "#d73027"))
colors <- col.pal(100)

# Calculate facet center for coloring
z.facet.center <- (z_clipped[-1, -1] + z_clipped[-1, -ncol(z_clipped)] + z_clipped[-nrow(z_clipped), -1] + z_clipped[-nrow(z_clipped), -ncol(z_clipped)]) / 4
z.facet.range <- cut(z.facet.center, 100)  # Range of the facet center on a 100-scale (number of colors)

# export plot
tiff(paste0(directory, "/Figures/FigS8A_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)

persp(x,y,z_clipped, theta=220, phi=25,
      xlab="STI_mean", ylab="STI_sd", zlab="Nsites_hist",
      zlim = c(z_min, z_max),
      col=colors[z.facet.range], expand = .8, r=2, shade=.4, border=NA,
      ticktype = "detailed") 

dev.off()

###### 7.2 European STI (Figure S8B) ####
# define variables
vars <- c("EU_temp.mean", "EU_temp.mean_sq", "EU_temp.sd")

# Get all possible combinations
combinations <- expand.grid(
  EU_temp.mean = c(0, 1),
  EU_temp.mean_sq = c(0, 1),
  EU_temp.mean_EU_temp.sd = c(0, 1),  # Add interaction term
  EU_temp.sd = c(0, 1))

# Convert combinations to formulas
formulas <- apply(combinations, 1, function(x) {
  terms <- c()
  if (x["EU_temp.mean"] == 1) terms <- c(terms, "EU_temp.mean")
  if (x["EU_temp.mean_sq"] == 1) terms <- c(terms, "EU_temp.mean_sq")
  if (x["EU_temp.mean_EU_temp.sd"] == 1) terms <- c(terms, "EU_temp.mean:EU_temp.sd")  # Add interaction term
  if (x["EU_temp.sd"] == 1) terms <- c(terms, "EU_temp.sd")
  formula <- paste("N_sites_hist ~", paste(terms, collapse = " + "))
  return(formula)
})

formulas[[1]] <- "N_sites_hist ~ 1"
formulas <- formulas[-c(3,7,11,15)] # remove all formulas with mean_temp_sq and without mean_temp
formulas # Print the formulas

## Run pgls
models_list <- list()
models_output <- data.frame(matrix(ncol=8, nrow=0))
colnames(models_output) <- c("m_num","formula","R2","Adj_R2","df","LogLik","AICc","lambda")
num_list <- 1:12

for (i in num_list) {
  result <- pgls_model_function(i, num_list, formulas, cdat, models_list, models_output)
  models_list <- result$models_list
  models_output <- result$models_output
}

models_output$AICc <- as.numeric(as.character(models_output$AICc))
min_AICc <- min(models_output$AICc) # Find the minimum AICc
best_models <- models_output[models_output$AICc <= (min_AICc + 2), ] # Select models with AICc less than 2 from the best model

best_models <- best_models %>% 
  arrange(AICc)

best_models # print best models
best_model_names <- best_models$m_num # Extract model names from the best models

best_models$m_num <- "Nsites_hist_vs_STI_EU" 
temp <- best_models %>% 
  rename(response_variable=m_num)
A_model_outputs_Nsites <- rbind(A_model_outputs_Nsites, temp) # append to main table

best_model_list <- lapply(best_model_names, function(x) models_list[[x]]) # Extract the actual models from models_list
model_avg <- model.avg(best_model_list) # Perform model averaging

## get averaged model coefficients
avg_model_summ <- pgls_avg_models_coefs(best_model_list, "Nsites_hist_vs_STI_EU")
A_coefs_export <- rbind(A_coefs_export, avg_model_summ) # save to table to export

# Set ranges for 3D graph
x <- range(filter_data$EU_temp.mean)
x <- seq(x[1], x[2], length.out = 50)

y <- range(filter_data$EU_temp.sd)
y <- seq(y[1], y[2], length.out = 50)

grid <- expand.grid(EU_temp.mean = x, EU_temp.sd = y) # Create a grid of x (STImean) and y (STIsd)
grid$EU_temp.mean_sq <- grid$EU_temp.mean^2 # Add the quadratic term for STImean to the grid

predictions <- predict(model_avg, newdata = grid) # Predict using the averaged model

# Reshape predictions into a matrix before persp plot
z <- matrix(predictions, nrow = length(x), ncol = length(y))

# Set desired z-axis limits
z_min <- 0  # minimum limit for z-axis
z_max <- 40   # maximum limit for z-axis

# Clip the z values to be within the specified range
z_clipped <- pmin(pmax(z, z_min), z_max)

# Set palette for plot
col.pal <- colorRampPalette(c("#4575b4","#91bfdb", "#e0f3f8","#fee090", "#fc8d59", "#d73027"))
colors <- col.pal(100)

# Calculate facet center for coloring
z.facet.center <- (z_clipped[-1, -1] + z_clipped[-1, -ncol(z_clipped)] + z_clipped[-nrow(z_clipped), -1] + z_clipped[-nrow(z_clipped), -ncol(z_clipped)]) / 4
z.facet.range <- cut(z.facet.center, 100)  # Range of the facet center on a 100-scale (number of colors)

# Export plot
tiff(paste0(directory, "/Figures/FigS8B_DDI-2024-0174.tif"), 
     width = 20, height = 15, units = "cm", res = 300)

persp(x,y,z_clipped, theta=220, phi=25,
      xlab="STI_mean", ylab="STI_sd", zlab="Nsites_hist",
      zlim = c(z_min, z_max),
      col=colors[z.facet.range], expand = .8, r=2, shade=.4, border=NA,
      ticktype = "detailed") 

dev.off()

## remove temporary files (not necessary)
rm(avg_model, avg_model_summ, best_model_list, best_models, coefs, combinations,
   cont_em, dredge.model, em, filter_data, filtered_models, grid, kk, m1, MC, lw, model_avg,
   models_list, models_output, mytheme, plot_CTI_hist, plot_CTI_rec, plot_qD_hist, plot_qD_rec,
   pred,result, subset_models, temp, vec, z, z_clipped, z.facet.center, best_model_names, colors, formulas,
   i, min_AICc, min_val, num_list, predictions, region, regions, Regions, regions_order,vars,x,y,z_max,z_min, z.facet.range)

##### 8. Export tables ####
###### 8.1 Model outputs (Table 2 and Table S8) ####
print(A_model_outputs_exp) # contains model outputs from dredge function for all GLMs and LMERs

## change names in explanatory variables for readability
A_model_outputs_exp$formula <- gsub("scale\\(historical_10y_sq\\)", "Tmax_hist_sq", A_model_outputs_exp$formula)
A_model_outputs_exp$formula <- gsub("scale\\(recent_10y_sq\\)", "Tmax_rec_sq", A_model_outputs_exp$formula)
A_model_outputs_exp$formula <- gsub("scale\\(diff_10y\\)", "Tmax_diff", A_model_outputs_exp$formula)
A_model_outputs_exp$formula <- gsub("scale\\(historical_10y\\)", "Tmax_hist", A_model_outputs_exp$formula)
A_model_outputs_exp$formula <- gsub("scale\\(recent_10y\\)", "Tmax_rec", A_model_outputs_exp$formula)

## export 
A_model_outputs_exp <- distinct(A_model_outputs_exp) # to remove any duplicates if present

write_xlsx(A_model_outputs_exp, paste0(directory,"/Tables/Table_GLM_and_LMER_model_outputs.xlsx"))

print(A_model_outputs_Nsites) # contains PGLS model outputs from occupancy analysis

## change names in explanatory variables for readability
A_model_outputs_Nsites$formula <- gsub("IBER_mean_temp_sq","STImean_sq", A_model_outputs_Nsites$formula)
A_model_outputs_Nsites$formula <- gsub("IBER_mean_temp","STImean", A_model_outputs_Nsites$formula)
A_model_outputs_Nsites$formula <- gsub("IBER_sd_temp","STIsd", A_model_outputs_Nsites$formula)
A_model_outputs_Nsites$formula <- gsub("EU_temp.mean_sq","STImean_sq", A_model_outputs_Nsites$formula)
A_model_outputs_Nsites$formula <- gsub("EU_temp.mean","STImean", A_model_outputs_Nsites$formula)
A_model_outputs_Nsites$formula <- gsub("EU_temp.sd","STIsd", A_model_outputs_Nsites$formula)

A_model_outputs_Nsites <- A_model_outputs_Nsites %>% 
  mutate(across(c("R2","Adj_R2","df","LogLik","AICc","lambda"), as.numeric))

A_model_outputs_Nsites <- distinct(A_model_outputs_Nsites) # to remove any duplicates if present

## export
write_xlsx(A_model_outputs_Nsites, paste0(directory,"/Tables/Table_PGLS_model_outputs.xlsx"))

###### 8.2 Model coefficients (Table S7) ####
(A_coefs_export)

## change names in explanatory variables for readability
A_coefs_export$select_variables <- gsub("scale\\(historical_10y_sq\\)", "Tmax_hist_sq", A_coefs_export$select_variables)
A_coefs_export$select_variables <- gsub("scale\\(recent_10y_sq\\)", "Tmax_rec_sq", A_coefs_export$select_variables)
A_coefs_export$select_variables <- gsub("scale\\(diff_10y\\)", "Tmax_diff", A_coefs_export$select_variables)
A_coefs_export$select_variables <- gsub("scale\\(historical_10y\\)", "Tmax_hist", A_coefs_export$select_variables)
A_coefs_export$select_variables <- gsub("scale\\(recent_10y\\)", "Tmax_rec", A_coefs_export$select_variables)

A_coefs_export$select_variables <- gsub("IBER_mean_temp_sq","STImean_sq", A_coefs_export$select_variables)
A_coefs_export$select_variables <- gsub("IBER_mean_temp","STImean", A_coefs_export$select_variables)
A_coefs_export$select_variables <- gsub("IBER_sd_temp","STIsd", A_coefs_export$select_variables)
A_coefs_export$select_variables <- gsub("EU_temp.mean_sq","STImean_sq", A_coefs_export$select_variables)
A_coefs_export$select_variables <- gsub("EU_temp.mean","STImean", A_coefs_export$select_variables)
A_coefs_export$select_variables <- gsub("EU_temp.sd","STIsd", A_coefs_export$select_variables)

A_coefs_export <- distinct(A_coefs_export) # to remove any duplicates if present

write_xlsx(A_coefs_export, paste0(directory,"/Tables/TableS7_model_coefficients.xlsx"))

