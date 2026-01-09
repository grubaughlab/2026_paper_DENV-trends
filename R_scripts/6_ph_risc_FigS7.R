setwd("/Users/abbeyporzucek/Library/Mobile Documents/com~apple~CloudDocs/Yale/Summer 2024/DENV Aim")

library(tidyverse)
library(INLA)
library(vroom)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(patchwork)
library(kableExtra)
library(plotly)
library(janitor)
library(stats)
library(zoo)
library(colorspace)
library(htmlwidgets)
library(scales)

#Note: RISc is coded here as dis (due to a previous naming of the variable)

#Read in data-------------------------------------------------------------------
#Read in outbreaks data
offset_cutoff <- read.csv("risc_cutoff_11.21.2025.csv")
#Read in WHO data:
#This was downloaded from https://github.com/cghss/dons/blob/master/Data/DONdatabase.csv on 7.16.2025
#From Georgetown's cleaned WHO DON data (https://gida.ghscosting.org/)
#I then went through and manually verified/cleaned everything
don <- read.csv("denv_dons_cleaned_reviewed.csv") 

#Compare countries included in DON to those included in outbreak analysis
#Make list of countries in our analysis
dis_countries <- unique(offset_cutoff$ISO_A0)
#Flag countries in don that are not in dis
don <- don %>%
  mutate(in_dis = ifelse(ISO %in% dis_countries, "yes", "no"))
#Make df showing just don countries and if they are in dis_countries 
don_countries <- don %>%
  distinct(ISO, Country, in_dis)
#Reunion island is in there twice so get rid of one of the lines
don_countries <- don_countries %>%
  filter(!Country == "La R-union Island")

# #Clean up the don data
#   don$ReportDate <- as.Date(don$ReportDate, "%m/%d/%Y")
#   #add a year column
#   don$year <- as.numeric(format(don$ReportDate, "%Y"))
#   #now select only the columns of interest
#   don_clean <- don[, c(1:10, 53:54)]
#Add in a column that creates a unique id for each country/year with an outbreak
don_clean <- don %>%
  mutate(id = paste(ISO, year, sep = "_")) %>%
  filter(!id == "ECU_2002") #no record found for this observation so it is removed
#   
#   #write CSV to further clean in excel
#   write.csv(don_clean, "don_clean.csv")
# don_clean <- read.csv("don_clean.csv")

#Cross check don data against dis data------------------------------------------
#Make a data frame of outbreaks identified by DIS
dis_outbreaks <- offset_cutoff %>%
  filter(intensity == "high") %>% #filter to only include outbreaks
  mutate(id = paste(ISO_A0, adjusted_year, sep = "_"),
         id_full = paste(full_name, adjusted_year, sep = " ")) #create an id column that matches don_clean$id
dis_ids <- dis_outbreaks$id  

#Flag which DON outbreaks are also DIS outbreaks 
don_temp <- don_clean %>%
  mutate(dis_outbreak = case_when(
    in_dis == "no" ~ NA_character_,
    in_dis == "yes" & id %in% dis_ids ~ "yes",
    in_dis == "yes" & !(id %in% dis_ids) ~ "no"),
    id_full = paste(Country, year, sep = " "))

#Now I want to add in DIS scores for all years in DON
#First create ID column in offset_cutoff
offset_cutoff <- offset_cutoff %>%
  mutate(id = paste(ISO_A0, adjusted_year, sep = "_"))
#Add in DIS for each outbreak year
don_dis <- don_temp %>%
  left_join(offset_cutoff %>% select(population_total, id, dis_mean), by = "id")
#Now filter for only countries included in the dis analysis
don_dis_overlap <- don_dis %>%
  filter(in_dis == "yes") %>%
  arrange(year, ISO) %>%
  mutate(id = factor(id, levels = unique(id)))  # Order id by year and ISO

#Add in Gideon data-------------------------------------------------------------
#Read in the data from gideon
gid <- read.csv("PHoutbreaks.csv")

#create lookup codes
codes <- offset_cutoff %>%
  distinct(full_name, ISO_A0)

#Add in iso codes
gid <- gid %>%
  left_join(codes, by = c("Country" = "full_name")) %>%
  na.omit() #get rid of countries that aren't in our data set (only one ob for Solomon Islands)

#Clean up the data set
gid_clean <- gid %>%
  rename(gid_outbreak = Outbreak) %>% #rename the outbreak column to denote gideon outbreaks
  mutate(id = paste(ISO_A0, adjusted_year, sep = "_"),
         id_full = paste(Country, adjusted_year, sep = " "))  #Create column to join with larger data set

#Add in DIS for each outbreak year
gid_dis <- gid_clean %>%
  left_join(offset_cutoff %>% select(population_total, id, dis_mean), by = "id")



#Bar graphs showing RISc DON comparison with different cutoff values-------------
# Assuming don_dis_overlap is your data frame containing the necessary data
# Arrange the data frame based on year to maintain chronological order within country
don_dis_overlap <- don_dis_overlap %>%
  arrange(year, id_full)

# Convert id_full to a factor and order it alphabetically by country and then by year
don_dis_overlap$id_full <- factor(don_dis_overlap$id_full, levels = unique(don_dis_overlap$id_full))

#Polt Figure S7a 
outbreak_cal_.5 <- ggplot(don_dis_overlap %>% filter(year > 2004)) +
  geom_bar(aes(x = id_full, y= dis_mean), stat = "identity", fill = "lightgrey")+
  geom_hline(yintercept = 0.5, color = "red", show)+
  annotate("text", x = Inf, y = 0.5, label = "Outbreak cutoff: RISc >= 0.5",
           hjust = 1.1, vjust = -0.5, color = "red", size = 3.5) +
  labs(title = "WHO vs RISc Outbreaks", 
       tag = "A)",
       x = "WHO Outbreak location and year", y = "Relative Intensity Score (RISc)")+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
outbreak_cal_.5  

# Arrange the data frame based on year to maintain chronological order within country
gid_dis <- gid_dis %>%
  arrange(adjusted_year, id_full)

# Convert id_full to a factor and order it alphabetically by country and then by year
gid_dis$id_full <- factor(gid_dis$id_full, levels = unique(gid_dis$id_full))

#Plot S7b
outbreak_gid <- ggplot(don_dis_overlap %>% filter(year > 2004)) +
  geom_bar(data = gid_dis, aes(x = id_full, y = dis_mean), stat = "identity", fill = "lightgrey")+
  geom_hline(yintercept = 0.5, color = "red", show)+
  annotate("text", x = Inf, y = 0.5, label = "Outbreak cutoff: RISc >= 0.5",
           hjust = 1.1, vjust = -0.5, color = "red", size = 3.5) +
  labs(title = "Outbreaks Reported in Literature vs RISc Outbreaks",
       tag = "B)",
       x = "Literature Outbreak location and year", y = "Relative Intensity Score (RISc)")+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
outbreak_gid  

#Final Figure 7 combined
outbreak_cal_.5 / outbreak_gid




