#set working directory
setwd("/Users/abbeyporzucek/Library/Mobile Documents/com~apple~CloudDocs/Yale/Summer 2024/DENV Aim")

#Load in libraries 
library(raster)
library(sf)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(exactextractr)

#Note: the final cleaned data files from this code are available and named:
  #"final_temp_anomaly_data_10_17_2025.csv"
  #"final_precip_anomaly_data_10_21_2025.csv"

#Initial setup data
#Read in raster stack  
anomaly_all <- stack("ERA5 Temp Anomaly Data.grib")
anomaly_all <- setMinMax(anomaly_all) #note: this takes a few minutes to run
#Read in country boundary data from rnatural earth
world_boundaries <- ne_countries(scale = "medium", returnclass = "sf")
#Re-project if CRS differ (example: re-projecting vector to raster CRS)
#This makes sure the temp data and the boundary have the same coordinates
if (sf::st_crs(world_boundaries) != crs(anomaly_all)) {
  world_boundaries <- st_transform(world_boundaries, crs(anomaly_all))
}

#Data exploration---------------------------------------------------------------
#Set up list to hold raster for each year
raster_year <- list()
temp_anom_year <- list()
temp_anom_final <-list()

#Iteration one - this is quick data exploration to figure out set up for function below
#Now want to extract just the first year (bands 1-12 - one band for each month)
raster_year$one<- anomaly_all[[c(1:12)]]
#Convert from Kelvin to Celsius 
raster_year$one <- raster_year$one - 273.15
#Extract mean temperature for each country boundary for each month
#So take the average of all raster_year$one values within the matching border of world_boundaries
mean_temperatures <- exact_extract(raster_year$one, world_boundaries, 'mean')
#Convert the list to a data frame
temp_anom_year$one <- as.data.frame(mean_temperatures)
colnames(temp_anom_year$one) <- paste0("mean_temp_", 1:12)
# Bind the calculated mean temperatures back to the corresponding countries in world_boundaries
temp_anom_final$one <- bind_cols(world_boundaries, temp_anom_year$one)
# Select relevant columns, add mean temperatures, calculate year mean temperature
temp_anom_final$one <- temp_anom_final$one %>%
  dplyr::select(adm0_a3, geometry, starts_with("mean_temp_")) %>%
  mutate(label_x = st_coordinates(st_centroid(geometry))[, 1],
         label_y = st_coordinates(st_centroid(geometry))[, 2]) 
temp_anom_final$one <- temp_anom_final$one %>%
  pivot_longer(2:13, names_to = "month.date", values_to = "month_mean") %>%
  mutate(month.date = gsub("mean_temp_", "1990_", month.date)) #replace the mean_temp value with year

#Iteration two
#Now want to extract just the first year (bands 1-12 - one band for each month)
raster_year$two<- anomaly_all[[c(13:24)]]
#Convert from Kelvin to Celsius 
raster_year$two <- raster_year$two - 273.15
#Extract mean temperature for each country boundary for each month
#So take the average of all raster_year$two values within the matching border of world_boundaries
mean_temperatures <- exact_extract(raster_year$two, world_boundaries, 'mean')
#Convert the list to a data frame
temp_anom_year$two <- as.data.frame(mean_temperatures)
colnames(temp_anom_year$two) <- paste0("mean_temp_", 1:12)
# Bind the calculated mean temperatures back to the corresponding countries in world_boundaries
temp_anom_final$two <- bind_cols(world_boundaries, temp_anom_year$two)
# Select relevant columns, add mean temperatures, calculate year mean temperature
temp_anom_final$two <- temp_anom_final$two %>%
  dplyr::select(adm0_a3, geometry, starts_with("mean_temp_")) %>%
  mutate(label_x = st_coordinates(st_centroid(geometry))[, 1],
         label_y = st_coordinates(st_centroid(geometry))[, 2]) 
temp_anom_final$two <- temp_anom_final$two %>%
  pivot_longer(2:13, names_to = "month.date", values_to = "month_mean") %>%
  mutate(month.date = gsub("mean_temp_", "1991_", month.date)) #replace the mean_temp value with year

#Extract temperature anomalies--------------------------------------------------
#Now turn the two iteration examples into a function 
temp_fun <- function(anomaly_stack, world_boundaries, start_year = 1990) {
  # Calculate number of years based on stack layers
  n_layers <- nlayers(anomaly_stack)
  n_years <- n_layers %/% 12
  
  # Initialize list to store final results
  temp_anom_final <- list()
  
  # Loop through each year
  for(i in 1:n_years) {
    cat("Processing year", start_year + i - 1, "...\n") #Extract rasters for each year
    
    # Calculate band indices for current year
    start_band <- (i - 1) * 12 + 1 #Pull out 12 bands at a time 
    end_band <- i * 12
    
    # Extract year's worth of data (12 months)
    year_raster <- anomaly_stack[[start_band:end_band]]
    
    # Convert from Kelvin to Celsius
    year_raster <- year_raster - 273.15
    
    # Extract mean temperature for each country boundary for each month
    mean_temperatures <- exact_extract(year_raster, world_boundaries, 'mean')
    
    # Convert to data frame
    temp_df <- as.data.frame(mean_temperatures) #change from raster to data frame
    colnames(temp_df) <- paste0("mean_temp_", 1:12) #Rename the columns with the month number 
    
    # Bind with country boundaries
    temp_with_boundaries <- bind_cols(world_boundaries, temp_df)
    
    # Process data: select columns, add centroids, pivot longer
    temp_processed <- temp_with_boundaries %>%
      dplyr::select(adm0_a3, geometry, starts_with("mean_temp_")) %>%
      mutate(label_x = st_coordinates(st_centroid(geometry))[, 1], #since I'm keeping the geometry the country centroid
             label_y = st_coordinates(st_centroid(geometry))[, 2]) %>% #get country centroid 
      pivot_longer(cols = starts_with("mean_temp_"), #go from wide format to long format 
                   names_to = "month.date", 
                   values_to = "month_mean") %>%
      mutate(month.date = gsub("mean_temp_", paste0(start_year + i - 1, "_"), month.date)) #put in a flag for the year
    
    # Store in list
    temp_anom_final[[i]] <- temp_processed
  }
  
  # Combine all years into one data frame
  final_result <- do.call(rbind, temp_anom_final) #bind each individual year into one df 
  
  return(final_result)
}

#Process full raster with temp_fun function
all_temp_data <- temp_fun(anomaly_all, world_boundaries, start_year = 1990)

#Save this output!
# saveRDS(all_temp_data, "global_country_level_temp_data_10.17.2025.RDS")


#Final clean up of data
  #Now make the dates read in as dates 
  all_temp_dates <- all_temp_data %>%
    mutate(month.date = ymd(paste0(substr(month.date, 1, 5), "-", substr(month.date, 6, 7), "-01")))
  #Get adjusted years for southern hemisphere
  temp_adjusted <- all_temp_dates %>%
    mutate(
      hemisphere = if_else(label_y >= 0, "north", "south"), #flag northern and southern hemispheres 
      adjusted_year = case_when(
        hemisphere == "north" ~ year(month.date), #if location is N, use the year from the month.date
        hemisphere == "south" & month(month.date) >= 7 ~ year(month.date), #if location is S, and month is 7-12 use year from month date
        hemisphere == "south" & month(month.date) < 7 ~ year(month.date) - 1))#if location is S, and month is 1-6 subtract 1 from the year 
  
  #Filter and process for what I need in the final analysis 
  #First read in final data set for analysis
  denv_offset_final <- read.csv("denv_offset_final_updated_regions_11.21.2025.csv")
  countries<- unique(denv_offset_final$ISO_A0)
  #Subset larger temp data set - select just the countries in our analysis
  temp_sub <- temp_adjusted %>%
    filter(adm0_a3 %in% countries)
  #Double check the number of filtered countries match number of analysis countries 
  unique(temp_sub$adm0_a3) #They match! 

  #Covert monthly data to annual data
  temp_annual <- temp_sub %>%
    group_by(adm0_a3, adjusted_year) %>% #group by country and then year 
    summarise(annual_temp_mean = mean(month_mean))
  #Save this data set!
  # saveRDS(temp_annual, "annual_temp_by_country_10.21.2025.RDS")
  
  #Now add in a column for the average temp between 1991-2020 (standard anomaly reference)
  final_temp_anomaly_data <- temp_annual %>%
    group_by(adm0_a3) %>%
    mutate(anom_ref = mean(annual_temp_mean[adjusted_year >= 1991 & adjusted_year <= 2020], na.rm = TRUE),
           temp_anomaly = annual_temp_mean - anom_ref)
  final_temp_anomaly_data <- as.data.frame(final_temp_anomaly_data) %>%
    dplyr::select(-geometry)
  #Save this data!
  # write.csv(final_temp_anomaly_data, "final_temp_anomaly_data_10_17_2025.csv")
  

#Extract precipitation anomalies------------------------------------------------
#Initial setup data
#Read in entire raster stack
precip_all <- stack("global precip anomaly ERA5 data.grib")
precip_all <- setMinMax(precip_all) #note: this takes a few mintues to run
#Read in country boundary data from rnatural earth (if you didn't above)
world_boundaries <- ne_countries(scale = "medium", returnclass = "sf")
# Re-project if CRS differ (example: re-projecting vector to raster CRS)
#this just makes sure the temp data and the boundary have the same coordinates
if (sf::st_crs(world_boundaries) != crs(precip_all)) {
  world_boundaries <- st_transform(world_boundaries, crs(precip_all))
}

#Set up list to hold rasters for each year
raster_year <- list()
precip_anom_year <- list()
precip_anom_final <-list()

#Again starting with quick data exploration
  #Iteration one
  #Now want to extract just the first year (bands 1-12 - one band for each month)
  raster_year$one<- precip_all[[c(1:12)]]
  #Extract mean temperature for each country boundary for each month
  #So take the average of all raster_year$one values within the matching border of world_boundaries
  mean_precip <- exact_extract(raster_year$one, world_boundaries, 'mean')
  #Convert the list to a data frame
  precip_anom_year$one <- as.data.frame(mean_precip)
  colnames(precip_anom_year$one) <- paste0("mean_precip_", 1:12)
  # Bind the calculated mean temperatures back to the corresponding countries in world_boundaries
  precip_anom_final$one <- bind_cols(world_boundaries, precip_anom_year$one)
  # Select relevant columns, add mean temperatures, calculate year mean temperature
  precip_anom_final$one <- precip_anom_final$one %>%
    dplyr::select(adm0_a3, geometry, starts_with("mean_precip_")) %>%
    mutate(label_x = st_coordinates(st_centroid(geometry))[, 1],
           label_y = st_coordinates(st_centroid(geometry))[, 2]) 
  precip_anom_final$one <- precip_anom_final$one %>%
    pivot_longer(2:13, names_to = "month.date", values_to = "month_mean") %>%
    mutate(month.date = gsub("mean_precip_", "1990_", month.date)) #replace the mean_temp value with year

#Now turn the above example into a function 
precip_fun <- function(anomaly_stack, world_boundaries, start_year = 1990) {
  # Calculate number of years based on stack layers
  n_layers <- nlayers(anomaly_stack)
  n_years <- n_layers %/% 12
  
  # Initialize list to store final results
  precip_anom_final <- list()
  
  # Loop through each year
  for(i in 1:n_years) {
    cat("Processing year", start_year + i - 1, "...\n") #Extract rasters for each year
    
    # Calculate band indices for current year
    start_band <- (i - 1) * 12 + 1 #Pull out 12 bands at a time 
    end_band <- i * 12
    
    # Extract year's worth of data (12 months)
    year_raster <- anomaly_stack[[start_band:end_band]]
    
    # Extract mean precipitation for each country boundary for each month
    mean_precip <- exact_extract(year_raster, world_boundaries, 'mean')
    
    # Convert to data frame
    precip_df <- as.data.frame(mean_precip) #change from raster to data frame
    colnames(precip_df) <- paste0("mean_precip_", 1:12) #Rename the columns with the month number 
    
    # Bind with country boundaries
    precip_with_boundaries <- bind_cols(world_boundaries, precip_df)
    
    # Process data: select columns, add centroids, pivot longer
    precip_processed <- precip_with_boundaries %>%
      dplyr::select(adm0_a3, geometry, starts_with("mean_precip_")) %>%
      mutate(label_x = st_coordinates(st_centroid(geometry))[, 1], #since I'm keeping the geometry the country centroid
             label_y = st_coordinates(st_centroid(geometry))[, 2]) %>% #get country centroid 
      pivot_longer(cols = starts_with("mean_precip_"), #go from wide formate to long formate 
                   names_to = "month.date", 
                   values_to = "month_mean") %>%
      mutate(month.date = gsub("mean_precip_", paste0(start_year + i - 1, "_"), month.date)) #put in a flag for the year
    
    # Store in list
    precip_anom_final[[i]] <- precip_processed
  }
  
  # Combine all years into one data frame
  final_result <- do.call(rbind, precip_anom_final) #bind each individual year into one df 
  
  return(final_result)
}

#Process full raster with temp_fun function
all_precip_data <- precip_fun(precip_all, world_boundaries, start_year = 1990)

#Save this output!
# saveRDS(all_precip_data, "global_country_level_precip_data_10.21.2025.RDS")

#Final clean up of data
  #Now make the dates read in as dates 
  all_precip_dates <- all_precip_data %>%
    mutate(month.date = ymd(paste0(substr(month.date, 1, 5), "-", substr(month.date, 6, 7), "-01")))
  #Get adjusted years for southern hemisphere
  precip_adjusted <- all_precip_dates %>%
    mutate(
      hemisphere = if_else(label_y >= 0, "north", "south"),
      adjusted_year = case_when(
        hemisphere == "north" ~ year(month.date), #if location is N, use the year from the month.date
        hemisphere == "south" & month(month.date) >= 7 ~ year(month.date), #if location is S, and month is 7-12 use year from month date
        hemisphere == "south" & month(month.date) < 7 ~ year(month.date) - 1))#if location is S, and month is 1-6 subtract 1 from the year 
  
  #From here just going to filter and process for what I need in the final analysis 
  #First read in final data set for analysis
  denv_offset_final <- read.csv("denv_offset_final_updated_regions_11.21.2025.csv")
  countries<- unique(denv_offset_final$ISO_A0)

  #Subset larger temp data set 
  precip_sub <- precip_adjusted %>%
    filter(adm0_a3 %in% countries)
  #Double check the number of filter countries match number of analysis countries 
  unique(precip_sub$adm0_a3) #They match! 


precip_annual <- precip_sub %>%
  group_by(adm0_a3, adjusted_year) %>% #group by country and then year 
  summarise(annual_precip_mean = mean(month_mean))
#Save this data set!
# saveRDS(precip_annual, "annual_precip_by_country_10.21.2025.RDS")

#Now add in a column for the average temp between 1991-2020 (standard anomaly reference)
final_precip_anomaly_data <- precip_annual %>%
  group_by(adm0_a3) %>%
  mutate(anom_ref = mean(annual_precip_mean[adjusted_year >= 1991 & adjusted_year <= 2020], na.rm = TRUE),
         precip_anomaly = annual_precip_mean - anom_ref)
final_precip_anomaly_data <- as.data.frame(final_precip_anomaly_data) %>%
  dplyr::select(-geometry)
#Save this data!
# write.csv(final_precip_anomaly_data, "final_precip_anomaly_data_10_21_2025.csv")



