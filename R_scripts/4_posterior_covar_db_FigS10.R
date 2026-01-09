#set working directory
setwd("/Users/abbeyporzucek/Library/Mobile Documents/com~apple~CloudDocs/Yale/Summer 2024/DENV Aim")

#load libraries 
library(tidyverse)
library(janitor)
library(rsoi)
library(jsonlite)
library(zoo)
library(rnaturalearth)
library(sf)
library(corrplot)
library(dplyr)
library(naniar)
library(DataExplorer)
library(ncdf4)

#Note: the final data file produced from this code is available and named "denv_covar_db_12.16.2025.csv"

#Generate posterior-------------------------------------------------------------
#Note: you can read in the output from the posterior below instead of generating it yourself
#Read in data set
gam_offset <- readRDS("poisson_updated_region_11.21.25.RDS")

#Get posterior sample for entire model
posterior.sample <- inla.posterior.sample(gam_offset, 
                                          n = 10000, 
                                          num.threads = 4, 
                                          add.names = T, 
                                          parallel.configs = T)

#Extract posterior for global random effect
#Return YearN
fun1 <- function() return (YearN)
#Another function to actually get the values from the posterior
posterior_YearN <- inla.posterior.sample.eval(fun = fun1, 
                                              samples = posterior.sample)
#Check posterior values with histogram
hist(posterior_YearN[1,], prob=TRUE, n=300)

#Extract posterior for regional random effect
#Return continent random effects
fun2 <- function() return (YearN_continent)
#Another function to actually get the values from the posterior
posterior_YearN_continent <- inla.posterior.sample.eval(fun = fun2, 
                                                        samples = posterior.sample)
#Make histogram of the posterior values
hist(posterior_YearN_continent[50,], prob=TRUE, n=300)

#Extract posterior for country random effect
#Return country random effects
fun3 <- function() return (YearN_country)
#Another function to get values
posterior_YearN_country <- inla.posterior.sample.eval(fun = fun3, 
                                                      samples = posterior.sample)
#Make histogram of values
hist(posterior_YearN_country[1900,], prob=TRUE, n=300)

#Get posterior mean and sd------------------------------------------------------  
#First pull out the IDs from the INLA output
YearN_IDs <- gam_offset$summary.random$YearN$ID
YearN_continent_IDs <- gam_offset$summary.random$YearN_continent$ID
YearN_country_IDs <- gam_offset$summary.random$YearN_country$ID

#Global level posterior - calculate observation level mean and sd
global_post<- as.data.frame(posterior_YearN) 
global_post <- global_post %>%
  rowwise() %>% #opperate across the data frame
  mutate(global_post_mean = mean(c_across(everything())), #calculate the mean
         global_post_sd = sd(c_across(everything()))) %>% #calculate the standard deviation
  ungroup() %>%
  mutate(YearN = YearN_IDs) %>%
  select(global_post_mean, global_post_sd, YearN) #keep on the needed columns

#Regional level posterior - calculate observation level mean and sd
region_post<- as.data.frame(posterior_YearN_continent) 
region_post <- region_post %>%
  rowwise() %>%
  mutate(region_post_mean = mean(c_across(everything())),
         region_post_sd = sd(c_across(everything()))) %>%
  ungroup() %>%
  mutate(region_join = YearN_continent_IDs) %>% #create unique id to join data later
  select(region_post_mean, region_post_sd, region_join)


#Country level posterior - calculate observation level mean and sd
country_post<- as.data.frame(posterior_YearN_country) 
country_post <- country_post %>%
  rowwise() %>%
  mutate(country_post_mean = mean(c_across(everything())),
         country_post_sd = sd(c_across(everything()))) %>%
  ungroup() %>%
  mutate(country_join = YearN_country_IDs) %>% #create unique id to join data later
  select(country_post_mean, country_post_sd, country_join)


#Join the mean posteriors back into the main data set 
#Read in data set with RISc values
offset_dis <- read.csv("final_risc_11.21.2025.csv")
#create join columns in offset_dis
offset_dis <- offset_dis %>%
  mutate(region_join = paste(continent, YearN, sep = "_"),
         country_join = paste(ISO_A0, YearN, sep = "_"))
#Join posterior data into offset_dis
output_with_post <-offset_dis %>%
  left_join(global_post, by = "YearN") %>%
  left_join(region_post, by = "region_join") %>%
  left_join(country_post, by = "country_join")


#Sum the posterior mean and sd
output_with_post <- output_with_post %>%
  rowwise() %>%
  mutate(total_post_mean = sum(c(global_post_mean, region_post_mean, country_post_mean)),
         total_post_sd = sum(c(global_post_sd, region_post_sd, country_post_sd)),
         post_variance = total_post_sd^2)

ggplot(data = output_with_post) + geom_point(aes(x = dis_variance, y = post_variance))

#save the output
# write.csv(output_with_post, "posterior_variance_12.16.2025.csv")

#Create database with covariate info--------------------------------------------
#---Loading in denv data--------------------------------------------------------
#Read in data
gam_output <- read.csv("posterior_variance_12.16.2025.csv") %>%
  mutate(intensity = case_when(
    dis_mean <= -0.5 ~ "low",
    dis_mean >= 0.5 ~ "high",
    dis_mean > -0.5 & dis_mean < 0.5 ~ "expected"
  ))

#create dataset with only relevant DENV output
denv_output <- subset(gam_output, 
                      select = c("adjusted_year", "continent", "ISO_A0",  "full_name", "cases", "dis_mean", "intensity", "post_variance", "total_post_mean")) %>%
  rename(Year = "adjusted_year")

#---Load and clean WB data------------------------------------------------------
#World Bank data 
#read in first half of world bank data and clean data
wb <- read.csv("WB_part_1.csv")
#get list of unique ISO_A0 codes from DENV data
unique_iso_codes <- unique(denv_output$ISO_A0)
#create list of WB categories of interest
#read in dataset with WB codes of interest
denv_covar_list <- read.csv("DENV Covariates.csv")
#create the list of WB categories of interest
wb_indicators <- unique(denv_covar_list$WB.Indicator.Code)
#filter to only include data on countries with DENV data and indicators of interest
wb_data <- wb %>%
  filter(Country.Code %in% unique_iso_codes) %>% #filter to selected countries
  filter(Indicator.Code %in% wb_indicators) %>%
  clean_names() #clean the column names using janitor package
#pivoting to get the data in a format compatible with DENV data
wb_data <- pivot_longer(wb_data, (5:68), names_to = "Year") 
#remove indicator_code column (if this becomes an issue later, remove the indicator_name column and use the code instead below)
wb_data = dplyr::select(wb_data, -4)
#remove the "x" from the years in the Year column
wb_data$Year = gsub("x", "", wb_data$Year)
#move the indicator codes to column names
wb_data <- pivot_wider(wb_data, names_from = indicator_name, values_from = value) 
#Note: so now the data is formatted so each country has one row per year with all the necessary data across the row


#Pull in second half of WB data 
#read in CSV
wb_v2 <- read.csv("WB data v2 12.10.2024.csv")
#filter to only include countries of interest (in the unique_iso_codes list)
wb_v2 <- wb_v2 %>%
  filter(Country.Code %in% unique_iso_codes)
#pivoting to get the data in a format compatible with DENV data
wb_data_v2 <- pivot_longer(wb_v2, (5:68), names_to = "Year") 
#Fixing the year - only keeping the characters for the year (character positions 2-5) 
wb_data_v2$Year = substr(wb_data_v2$Year, 2, 5)
#Remove the Series.Code variable 
wb_data_v2 <- wb_data_v2 %>%
  dplyr::select(-"Series.Code") %>%
  filter(Year >= 1980)
#move the indicator codes to column names
wb_data_v2 <- pivot_wider(wb_data_v2, names_from = Series.Name, values_from = value) %>%
  clean_names() #clean up the column names 

#Pull in population data from the WB
wb_pop <- read.csv("WB population data 2025.csv")
#pivoting to get the data in a format compatible with DENV data
wb_pop_data <- pivot_longer(wb_pop, (5:69), names_to = "Year") 
#Fixing the year - only keeping the characters for the year (character positions 2-5) 
wb_pop_data$Year = substr(wb_pop_data$Year, 2, 5)
#Remove the Indicator.Code and Name variables
wb_pop_data <- wb_pop_data %>%
  dplyr::select(-"Indicator.Code", -"Indicator.Name") %>%
  filter(Year >= 1990) %>%
  clean_names() #clean up column names 
#rename the "values" column to "population_total"
names(wb_pop_data)[names(wb_pop_data) == "value"] <- "population_total"

#Read in land area data from WB
wb_land_area <- read.csv("WB land area.csv")
#pivot to get the data in compatible format
wb_land <- pivot_longer(wb_land_area, (5:54), names_to = "year")
#Fix the year info
wb_land$year = substr(wb_land$year, 2,5)
#Remove the Series.Code variable 
wb_land <- wb_land %>%
  dplyr::select(-"Series.Code", -"Series.Name") %>%
  filter(year >= 1980) %>%
  clean_names()
#rename "values" variable to "land_area_sqkm"
names(wb_land)[names(wb_land) == "value"] <- "land_area_sqkm"


#---Combine dengue data with world bank data------------------------------------
#Combine denv_outputs with wb_data 
#create a column with unique identifiers in each data frame
denv_output$join <- paste(denv_output$ISO_A0, denv_output$Year)
wb_data$join <- paste(wb_data$country_code, wb_data$Year) 

#join the two data sets into one new one
denv_covar_db <- left_join(denv_output, wb_data, by = "join")
#clean up the new db
#rename the Year.x column
names(denv_covar_db)[names(denv_covar_db) == "Year.x"] <- "Year"
#remove redundant rows
denv_covar_db <- denv_covar_db %>%
  dplyr::select(-"country_code", -"Year.y")
#reorder the rows and clean column names 
denv_covar_db <- denv_covar_db[,c(9, 1:8, 10:17)] %>%
  clean_names() #put all column names into snake format


#Combine second half of world bank data 
#create a join column in wb_data_v2
wb_data_v2$join <- paste(wb_data_v2$country_code, wb_data_v2$year) 

#join the two data sets into one new one
denv_covar_db <- left_join(denv_covar_db, wb_data_v2, by = "join")

#clean up the new db
#rename the Year.x column
names(denv_covar_db)[names(denv_covar_db) == "year.x"] <- "year"
#rename the country_name.x column
names(denv_covar_db)[names(denv_covar_db) == "country_name.x"] <- "country_name"
#remove redundant rows
denv_covar_db <- denv_covar_db %>%
  dplyr::select(-"country_name.y", -"country_code", -"year.y")

#Join population data into the denv_covar_db
#create a join column in wb_pop_data
wb_pop_data$join <- paste(wb_pop_data$country_code, wb_pop_data$year)
#join the two datasets 
denv_covar_db <- left_join(denv_covar_db, wb_pop_data, by = "join")
#clean up the new db
#rename the Year.x column
names(denv_covar_db)[names(denv_covar_db) == "year.x"] <- "year"
#rename the country_name.x column
names(denv_covar_db)[names(denv_covar_db) == "country_name.x"] <- "country_name"
#remove redundant rows
denv_covar_db <- denv_covar_db %>%
  dplyr::select(-"country_name.y", -"country_code", -"year.y")

#join land area data into the denv_covar_db
#create a join column in wb_land 
wb_land$join = paste(wb_land$country_code, wb_land$year)
#join the two datasets 
denv_covar_db <- left_join(denv_covar_db, wb_land, by = "join")
#clean up the new db
#rename the Year.x column
names(denv_covar_db)[names(denv_covar_db) == "year.x"] <- "year"
#rename the country_name.x column
names(denv_covar_db)[names(denv_covar_db) == "country_name.x"] <- "country_name"
#remove redundant rows
denv_covar_db <- denv_covar_db %>%
  dplyr::select(-"country_name.y", -"country_code", -"year.y")



#--Add immunoprevalence data---------------------------------------------------
#Add in column with binary values for high intensity (0) or low/mild intensity (1) years
denv_covar_db <- denv_covar_db %>%
  arrange(iso_a0, year) %>%
  group_by(iso_a0) %>% #organize data by country 
  mutate(high_intensity = ifelse(intensity == "high", 0, 1)) %>%
  group_by(iso_a0, grp = cumsum(high_intensity == 0)) %>%
  mutate(time_from_hi = cumsum(high_intensity)) %>%
  ungroup()%>%
  dplyr::select(-grp, -high_intensity)
#reorder to get the time since last outbreak as the 10th column
denv_covar_db <- denv_covar_db[,c(1:9, 29, 10:28)] 

#---Add in climate data---------------------------------------------------------
#Add ENSO data using the rsoi package
enso = download_enso()
#Get the average SOI (air pressure anomaly) per year
enso_sum <- enso %>%
  filter(Year >= 1980) %>% #filtering to only include the needed years
  group_by(Year) %>%
  summarise(annual_soi = mean(SOI)) %>% 
  ungroup() %>%
  clean_names() #put column names in snake case
#Add the enso data to denv_covar_db
denv_covar_db <- left_join(denv_covar_db, enso_sum, by= "year")

#Add in ONI (sea surface temperature anomaly) data 
#Summarize ENSO data annually
enso_summary <- enso %>%
  filter(Year >= 1980) %>%
  group_by(Year) %>%
  summarise(
    cool_oni_avg = sum(ifelse(grepl("Cool Phase", phase), ONI, 0)) / 12, #get avg cool temp
    warm_oni_avg = sum(ifelse(grepl("Warm Phase", phase), ONI, 0)) / 12, #get avg warm temp
    cool_phase_months = sum(phase == "Cool Phase/La Nina"), #count cool phase months
    neutral_phase_months = sum(phase == "Neutral Phase"), #count neutral phase months
    warm_phase_months = sum(phase == "Warm Phase/El Nino") #count warm phase months 
  ) %>%
  clean_names()

#Add the enso phase data to denv_covar_db
denv_covar_db <- left_join(denv_covar_db, enso_summary, by = "year")


#Add in cleaned temperature anomaly data
temp_anomaly <- read.csv("final_temp_anomaly_data_10_17_2025.csv")  
#Add a join column into the temp_anomaly data
temp_anomaly <- temp_anomaly %>%
  mutate(join = paste(adm0_a3, adjusted_year, sep = " "))
#Add the temp anomaly data to denv_covar_db
denv_covar_db <- left_join(denv_covar_db, temp_anomaly %>% dplyr::select(temp_anomaly, join), by = "join")

#Add in cleaned precipitation anomaly data
precip_anomaly <- read.csv("final_precip_anomaly_data_10_21_2025.csv")  
#Add a join column into the temp_anomaly data
precip_anomaly <- precip_anomaly %>%
  mutate(join = paste(adm0_a3, adjusted_year, sep = " "))
#Add the temp anomaly data to denv_covar_db
denv_covar_db <- left_join(denv_covar_db, precip_anomaly %>% dplyr::select(precip_anomaly, join), by = "join")



#Add in conflict data-----------------------------------------------------------
#Read in conflict data
#This data is the UCDP/Prio Armed Conflict Dataset version 25.1
#Downloaded from https://ucdp.uu.se/downloads/ on 7.18.2025
conflict <- readRDS("UcdpPrioConflict_v25_1 2.rds")
#Keep only the variables we care about 
conflict_sub <- conflict %>%
  rename(conflict_intensity = intensity_level) %>% #rename so we don't confuse with denv intensity
  dplyr::select(location, year, conflict_intensity, type_of_conflict) %>%
  mutate(id = paste(location, year, sep = "_")) #create and id column so we can join the data with denv data
#Add this into denv_covar_db
#first need to make a matching id column for the denv_covar_db
denv_covar_db <- denv_covar_db %>%
  mutate(id = paste(full_name, year, sep="_"))
denv_covar_db <- denv_covar_db %>%
  left_join(dplyr::select(conflict_sub, id, conflict_intensity, type_of_conflict), by = "id")
#Replace NAs with 0s
denv_covar_db$conflict_intensity[is.na(denv_covar_db$conflict_intensity)] <- 0

#---Add in DIS from previous year------------------------------------------------
denv_covar_db <- denv_covar_db %>%
  group_by(iso_a0) %>%
  mutate(previous_dis = lag(dis_mean, n=1)) #n=1 because we want to lag by one unit of rr

#Check for data class
sapply(denv_covar_db, class)
#Annoyingly WB data shows NA as ".." and there for is read in as a character
#Fix this by replacing ".." with NA
denv_covar_db <- denv_covar_db %>%
  mutate(across(where(is.character), ~na_if(., ".."))) #if there is a character variable with ".." replace ".." with NA
#Read in columns 10:39 as numeric
denv_covar_db <- denv_covar_db %>%
  dplyr::select(-id, -country_name) %>% #this removes an unnecessary id and country_name columns
  mutate(across(10:38, ~as.numeric(.)))
#Now double check data classes 
sapply(denv_covar_db, class)

#save final covar dataset
# write.csv(denv_covar_db, "denv_covar_db_12.16.2025.csv") #this has the posterior variance

denv_covar_db <- read.csv("denv_covar_db_12.16.2025.csv")

#check for data completeness 
plot_missing(denv_covar_db, title = "Covariate Data Completeness")

#Select variables to include in glm---------------------------------------------
#To run a correlation matrix, I need to filter to only include numeric variables
#first I am going to only select covariates that will be in the final GLM
covar_sub <- denv_covar_db %>%
  dplyr::select(year, continent, iso_a0, full_name, cases, dis_mean, post_variance, intensity, 
                time_from_hi, access_to_electricity_percent_of_population,
                forest_area_percent_of_land_area, 
                gdp_per_capita_ppp_constant_2021_international,
                population_ages_0_14_percent_of_total_population,
                urban_population_percent_of_total_population,
                cool_phase_months, warm_phase_months, temp_anomaly, precip_anomaly,
                previous_dis, conflict_intensity)

#Now check to make sure all necessary variables are numeric
sapply(covar_sub, class)

#Run correlation matrix---------------------------------------------------------
#select only numeric covariates for correlation matrix 
#Rename columns in covar_numeric
covar_numeric <- covar_sub %>%
  ungroup() %>%
  dplyr::select(previous_dis, 
                access_to_electricity_percent_of_population,
                forest_area_percent_of_land_area, 
                gdp_per_capita_ppp_constant_2021_international,
                urban_population_percent_of_total_population,
                temp_anomaly,
                precip_anomaly,
                year,
                # previous_dis
  ) %>%
  rename("Previous year RISc" = previous_dis,
         "Access to electricity" = access_to_electricity_percent_of_population,
         "Forest area" = forest_area_percent_of_land_area,
         # "Percent young population" = population_ages_0_14_percent_of_total_population,
         "GDP" = gdp_per_capita_ppp_constant_2021_international,
         "Urbanicity" = urban_population_percent_of_total_population,
         "Temperature anomaly" = temp_anomaly,
         "Precipitation anomaly" = precip_anomaly,
         # "Previous year DIS" = previous_dis,
         "Year" = year)
# "Conflict intensity" = conflict_intensity)

plot_missing(covar_numeric, title = "Selected Covariate Data Completeness")

#Save output from final data for fig S10
# write.csv(covar_numeric, "covar_numeric.csv")


#Figure S10---------------------------------------------------------------------
# first run a matrix with everything 
denv_matrix <- cor(covar_numeric, use = "complete.obs")
# Plot figure S10
#Note: for ease, the code for this plot is also in main figure file
corrplot(
  denv_matrix,
  tl.col = "black"
)  



