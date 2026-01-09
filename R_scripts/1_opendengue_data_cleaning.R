#Set working directory
setwd("")

#Load in libraries 
library(tidyverse)
library(vroom)
library(rnaturalearth)
library(rnaturalearthdata)
library(janitor)

#Note: the final output from this file is available and named "denv_offset_final_updated_regions_11.21.2025.csv"
  
#Cleaning OpenDengue Data-------------------------------------------------------
#Import case data - downloaded from OpenDengue July 8, 2025
#Import case data - downloaded from OpenDengue July 8, 2025
denv_temp <- vroom::vroom("Temporal_extract_V1_3.csv")
denv_nat <- vroom::vroom("National_extract_V1_3.csv")
#date formatting 
denv_nat$calendar_start_date <- as.Date(denv_nat$calendar_start_date, "%Y-%m-%d")
denv_temp$calendar_start_date <- as.Date(denv_temp$calendar_start_date, "%Y-%m-%d")
#create new variables called week.date and month.date 
denv_nat$week.date <- floor_date(denv_nat$calendar_start_date, unit='week')
denv_temp$week.date <- floor_date(denv_temp$calendar_start_date, unit='week')
denv_nat$month.date <- floor_date(denv_nat$calendar_start_date, unit='month')
denv_temp$month.date <- floor_date(denv_temp$calendar_start_date, unit='month')

#Merge denv_temp 2024 data with denv_nat_data
  #Prep national data set by removing data from 2024
  denv_nat_sub <- denv_nat %>%
    filter(!Year == "2024")
  #Prep temporal data set by filtering to only include what is needed
  denv_temp_sub <- denv_temp%>%
    filter(str_detect(RNE_iso_code, "^[A-Za-z]{3}$") & #filter to include only national level data
             Year == "2024" & #only include data from 2024
             case_definition_standardised == "Total" & #only include total cases (not suspected, etc)
             !T_res == "Year") #remove any annual resolution 
  #Now merge the two data sets together 
  denv_merge <- rbind(denv_nat_sub, denv_temp_sub)

#Add info on hemispheres and geographic regions using rnaturalearth
  #get world map data
  world <- ne_countries(scale = "medium", returnclass = "sf")
  #prep for join by getting world name data in all caps to match outbreak data
  world$sovereignt = toupper(world$sovereignt)
  #create a new data frame with names of countries and their regions
  regions <- world[,c("adm0_iso", "continent", "label_y")]
  #now join the regions data frame to denv_nat
  denv_nat <- denv_merge %>%
    left_join(regions %>% select(adm0_iso, continent, label_y),
              by = c("ISO_A0" = "adm0_iso"))
  # Keep only the necessary variables
  denv_temp <- denv_nat %>%
    select(full_name, ISO_A0, dengue_total, month.date, continent, label_y) 

#Offset the date by 6 months in the southern hemisphere
  #Flag the hemisphere and calculate "adjusted year" based on location
  denv_temp_adjusted <- denv_temp %>%
    mutate(
      hemisphere = if_else(label_y >= 0, "north", "south"),
      adjusted_year = case_when(
        hemisphere == "north" ~ year(month.date), #if location is northern hemisphere, use the year from the month.date
        hemisphere == "south" & month(month.date) >= 7 ~ year(month.date), #if location is southern hemisphere, and month is 7-12 use year from month date
        hemisphere == "south" & month(month.date) < 7 ~ year(month.date) - 1))#if location is southern, and month is 1-6 subtract 1 from the year 
  #So now if you look at the data set, months 1-6 are assigned to the previous year for southern hemisphere only

#Add cases by adjusted_year and keep all other variables
denv_yearly <- denv_temp_adjusted %>%
  group_by(ISO_A0, adjusted_year) %>% #group by country and adjusted year
  summarise(
    cases = sum(dengue_total, na.rm = TRUE), #sum all dengue_total values in adjusted_year 
    full_name = sub(",.*", "", first(full_name)), # keep everything before the first comma so we only keep the country name
    continent = first(continent),
    label_y = first(label_y),
    hemisphere = first(hemisphere),
    .groups = "drop" 
  ) %>%
  complete(ISO_A0, adjusted_year = 1990:2023) %>% #make sure every location has a row for 1990-2023
  mutate(join = paste(ISO_A0, adjusted_year, sep = " "))

#Filter for years of interest and data completeness  
denv_data_temp <- denv_yearly %>%
  filter(adjusted_year >= 1990 & adjusted_year <=2023)%>% #filter for years of interest
  group_by(ISO_A0) %>%
  mutate(total_missing = sum(is.na(cases) | cases == 0)) %>% #Add total number of years with missing data or 0 cases
  filter(total_missing <= 15) #Filter for data completeness 


#Fill missing values introduced by the complete() command
  # Create a reference table mapping ISO_A0 codes to full country names
  unique_values <- denv_data_temp %>%
    filter(!is.na(full_name)) %>%
    distinct(ISO_A0, full_name, continent, hemisphere)
  # Fill missing full_name values using a left join
  denv_complete <- denv_data_temp %>%
    left_join(unique_values, by = "ISO_A0", suffix = c("", "_lookup")) %>% #adding "_lookup" to the end of unique_values variables
    mutate(full_name = coalesce(full_name, full_name_lookup),
           continent = coalesce(continent, continent_lookup),
           hemisphere = coalesce(hemisphere, hemisphere_lookup)) %>%
    select(-ends_with("_lookup")) #removing all columns that end with _lookup
  #Reformat country names
  denv_complete$full_name = gsub("\\.", " ", denv_complete$full_name)
  denv_complete$full_name <- str_to_title(denv_complete$full_name)
  
#Clean up regions
  #Replace the region for the Maldieves with "Asia"
  denv_complete$continent <- ifelse(denv_complete$continent == "Seven seas (open ocean)", "Asia", denv_complete$continent)
  #Split "North America" into Central America and Caribbean
  #get a list of only unique countries in North America  
  north_america <- denv_complete %>%
    filter(continent == "North America")
  #over write north_america to pull out only unique country names
  north_america <- unique(north_america$full_name)
  #manually make a list for "Central America"
  c_am <- c("Belize", "Costa Rica", "Guatemala", "Honduras", "Mexico", "Nicaragua",
            "Panama", "El Salvador")
  #make a list of caribbean that includes all north_america countries not in c_am
  caribbean <- setdiff(north_america, c_am)
  #Now adjust the denv_offset_final df to show the updated continents
  denv_complete <- denv_complete %>%
    mutate(
      continent = case_when(
        full_name %in% c_am ~ "Central America",
        full_name %in% caribbean ~ "Caribbean",
        TRUE ~ continent))  # keep existing value if no match
  

#Add in additional columns for years (can't reuse year columns in INLA)
denv_complete <- denv_complete %>%
  group_by(ISO_A0)%>%
  mutate(YearN = row_number(),
         YearN2 = YearN,
         YearN3 = YearN,
         YearN4 = YearN)

#Add in population data 
  #Note: World Bank population data was downloaded on July 8, 2025 
  #read in world bank data and clean data
  wb_pop <- read.csv("WB population data 2025.csv")
  #get unique list of ISO codes from denv_subset
  unique_iso <- unique(denv_complete$ISO_A0)
  #filter wb data to only include unique ISO codes
  wb_pop <- wb_pop %>%
    filter(Country.Code %in% unique_iso)
  #pivoting to get the data in a format compatible with DENV data
  wb_pop_data <- pivot_longer(wb_pop, (5:69), names_to = "Year") 
  #Fixing the year - only keeping the characters for the year (character positions 2-5) 
  wb_pop_data$Year = substr(wb_pop_data$Year, 2, 5)
  #Remove the Series.Code variable 
  wb_pop_data <- wb_pop_data %>%
    dplyr::select(-"Indicator.Code", -"Indicator.Name") %>%
    filter(Year >= 1990) %>%
    clean_names() #clean up column names 
  #rename the "values" column to "population_total"
  names(wb_pop_data)[names(wb_pop_data) == "value"] <- "population_total"
  #create join column in wb_pop_data and in denv_subset
  wb_pop_data$join = paste(wb_pop_data$country_code, wb_pop_data$year)
  denv_complete$join = paste(denv_complete$ISO_A0, denv_complete$adjusted_year)  
  #join the two datasets together
  denv_pop <- left_join(denv_complete, wb_pop_data, by = "join")
  #get rid of redundant rows
  denv_pop <- denv_pop %>%
    dplyr::select(-"country_name", -"country_code", -"year")
  #Remove places with missing population data
  denv_pop_clean <- denv_pop %>%
    filter(adjusted_year >= 1990 & adjusted_year <=2023)%>%
    group_by(ISO_A0) %>%
    filter(!is.na(population_total)) %>%
    ungroup()
  
#create final dataset for GAM
denv_offset_final <- denv_pop_clean %>%
  group_by(ISO_A0) %>%
  mutate(ave_cases = mean(cases, na.rm = T),
         avg_inc = mean((cases/as.numeric(population_total))*100000, na.rm = T)) %>% #calculate avg incidence per 100,000 pop
  filter(ave_cases >= 500 | avg_inc >= 10) %>% #filter based on avg cases OR avg incidence 
  ungroup() %>%
  mutate(continent = if_else(is.na(continent),'other',continent),
         YearN_continent = as.factor(paste(continent, YearN, sep='_')),
         YearN_country = as.factor(paste(ISO_A0, YearN, sep='_')),
         offset1=as.numeric(population_total)/100000, #rerun with out dividing
         country=ISO_A0,
         countryID= as.numeric(as.factor(country)),
         continentID= as.numeric(as.factor(continent)))

#save final data set
write.csv(denv_offset_final, "denv_offset_final_updated_regions_11.21.2025.csv")


