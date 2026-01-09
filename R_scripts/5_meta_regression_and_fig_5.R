#set working directory
setwd("/Users/abbeyporzucek/Library/Mobile Documents/com~apple~CloudDocs/Yale/Summer 2024/DENV Aim")


#Global meta regression---------------------------------------------------------
library(metafor)

#Read in covariate data set with variance
denv_covar_var <- read.csv("denv_covar_db_12.16.2025.csv")
#Rename precip_anomaly to precipitation_anomaly
denv_covar_var <- denv_covar_var %>%
  rename(precipitation_anomaly = "precip_anomaly") 

#Scale all of the covariates
covar_scaled <- denv_covar_var %>%
  mutate(across(c("year", "access_to_electricity_percent_of_population", 
                  "forest_area_percent_of_land_area", 
                  "gdp_per_capita_ppp_constant_2021_international", 
                  "temp_anomaly",
                  "urban_population_percent_of_total_population", 
                  "precipitation_anomaly", 
                  "previous_dis"),
                ~ scale(.x)))


#Select covariates and remove rows with NAs
covar_scaled_clean <- covar_scaled %>% 
  dplyr::select(dis_mean, total_post_mean, post_variance, continent,
                year,
                access_to_electricity_percent_of_population,
                forest_area_percent_of_land_area,
                gdp_per_capita_ppp_constant_2021_international,
                temp_anomaly,
                urban_population_percent_of_total_population,
                precipitation_anomaly,
                previous_dis)%>%
  na.omit()


#Run the meta-regression with the clean dataset
global_inla_mr <- rma.uni(yi = total_post_mean,
                          vi = post_variance,
                          data = covar_scaled_clean,
                          method = "REML",
                          mods = ~ year +
                            access_to_electricity_percent_of_population +
                            forest_area_percent_of_land_area +
                            gdp_per_capita_ppp_constant_2021_international +
                            temp_anomaly +
                            urban_population_percent_of_total_population +
                            precipitation_anomaly +
                            previous_dis)

# Summarize model results
summary(global_inla_mr)

# Extract coefficients and confidence intervals
coef_summary <- summary(global_inla_mr)$beta
ci_lower <- summary(global_inla_mr)$ci.lb
ci_upper <- summary(global_inla_mr)$ci.ub

# Combine into a data frame
moderators <- c("Intercept", "Year", 
                "Access to Electricity",
                "Forest Area", 
                "GDP per Capita", "Temperature Anomaly", "Urbanicity", 
                "Precipitation Anomaly", "Previous Year RISc")

forest_data <- data.frame(
  Moderator = moderators,
  Estimate = coef_summary,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper
) %>% 
  slice(-1) %>%#this removed the intercept - if intercept is in final plot this is where to go to fix it
  arrange(Estimate) #Arrange data by effect size (Estimate)

# Ensure the Moderator is a factor with levels in the order of appearance
forest_data$Moderator <- factor(forest_data$Moderator, levels = forest_data$Moderator)


#Figure 5a - forest plot
global_mr_plot <- ggplot(forest_data, aes(x = Estimate, y = Moderator)) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), width = 0.5,
                linewidth = 0.5, color = "darkslategrey") +
  geom_point(color = "#172869FF", size = 3.5) +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
  labs(title = "Global association between covariates and RISc", 
       tag = "A)",
       x = "Effect Size", y = "") +
  theme_minimal()+
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        title = element_text(size = 14),
        plot.title = element_text(size = 18))
global_mr_plot


#Regional meta regresiion-------------------------------------------------------
#Create a vector with variable names to make plotting easier 
# Create named vector of labels
covar_labels <- c(
  temp_anomaly = "Temperature Anomaly",
  cool_oni_avg = "ENSO cool phase",
  conflict_intensity = "Conflict intensity",
  previous_dis = "Previous year RISc",
  urban_population_percent_of_total_population = "Urbanicity",
  precipitation_anomaly = "Precipitation anomaly",
  population_ages_0_14_percent_of_total_population = "Percent population under 14 years",
  gdp_per_capita_ppp_constant_2021_international = "GDP per Capita",
  forest_area_percent_of_land_area = "Forest Area",
  forest_change = "Deforestation",
  access_to_electricity_percent_of_population = "Access to Electricity",
  year = "Year",
  cool_phase_months = "ENSO cool phase months",
  time_from_hi = "Time since last outbreak",
  warm_oni_avg = "ENSO warm phase months"
)


#Run the meta-regression including the interaction with continent
region_inla_mr <- rma.uni(yi = total_post_mean,
                          vi = post_variance,
                          data = covar_scaled_clean,
                          method = "REML",
                          mods = ~ continent * (year +
                                                  access_to_electricity_percent_of_population +
                                                  forest_area_percent_of_land_area +
                                                  gdp_per_capita_ppp_constant_2021_international +
                                                  temp_anomaly +
                                                  urban_population_percent_of_total_population +
                                                  precipitation_anomaly +
                                                  previous_dis))

# Summarize model results
summary(region_inla_mr)

# Extract coefficients and confidence intervals
coef_r_summary <- summary(region_inla_mr)$beta
ci_lower_r <- summary(region_inla_mr)$ci.lb
ci_upper_r <- summary(region_inla_mr)$ci.ub

# Combine into a data frame
moderators <- rownames(coef_r_summary)

# Create a data frame including all coefficients and CIs
region_forest_data <- data.frame(
  Moderator = moderators,
  Estimate = coef_r_summary,
  CI_Lower = ci_lower_r,
  CI_Upper = ci_upper_r
)

#Now need to do some data cleaning
#Create a column for region
forest_region <- region_forest_data %>%
  #Make a region column by pulling region from the Moderator column, any place without a region is the refrence (Asia)
  mutate(region = ifelse(grepl("^continent.*:", Moderator), sub("continent(.*):.*", "\\1", Moderator), "Asia")) %>%
  #Remove the intercept data
  slice(-(1:5)) %>%
  #Remove the region info from the moderator since it has its own column now 
  mutate(Moderator = sub("^[^:]*:", "", Moderator))

#Add the estimates for Asia to each of the other regions
# Extract the estimates for Asia
asia_estimates <- forest_region %>%
  filter(region == "Asia") %>%
  #Flag the column names noting they are from Asia 
  dplyr::select(Moderator, Estimate_Asia = Estimate, CI_Lower_Asia = CI_Lower, CI_Upper_Asia = CI_Upper)
#Add the flagged asia data back to larger data set
region_temp <- forest_region %>%
  left_join(asia_estimates, by = "Moderator")
#Add the estimates from Asia to the other regions
region_forest_final <- region_temp %>%
  #Create estimate_final column, for all regions that are not Asia add the estimate to the estimate_asia
  #If the region is Asia, then leave the Estimate value in Estimate_final
  mutate(Estimate_final = ifelse(region != "Asia", Estimate + Estimate_Asia, Estimate),
         CI_Lower_final = ifelse(region !="Asia", CI_Lower+CI_Lower_Asia, CI_Lower),
         CI_Upper_final = ifelse(region !="Asia", CI_Upper+CI_Upper_Asia, CI_Upper)) 

#Arrange data
#Arrange Moderators in same order as global moderators 
region_forest_final$Moderator <- factor(region_forest_final$Moderator, 
                                        levels = rev(c("temp_anomaly", 
                                                       "previous_dis",
                                                       "gdp_per_capita_ppp_constant_2021_international",
                                                       "forest_area_percent_of_land_area",
                                                       "access_to_electricity_percent_of_population",
                                                       "urban_population_percent_of_total_population",
                                                       "precipitation_anomaly",
                                                       "year")))
#Arrange regions
region_forest_final$region <- factor(region_forest_final$region, 
                                     levels = c("Asia", "Oceania", "Caribbean", "Central America", "South America"))

#Make dummy variable to have a shaded region behind every other region facet
region_forest_shade <- region_forest_final %>%
  mutate(shade = case_when(region == "Asia" ~ 1,
                           region == "Oceania" ~ 0,
                           region == "Caribbean" ~ 1,
                           region == "Central America" ~ 0,
                           region == "South America" ~ 1))

#Figure 5b
region_mr_plot <- ggplot(region_forest_shade, aes(x = Estimate_final, y = Moderator, color = region)) +
  geom_rect( #add shading behind every other facet panel
    data = region_forest_shade %>% filter(shade == 1) %>% distinct(region, shade),
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, group = region),
    fill = "#C7C1B8", alpha = 0.25) +
  geom_errorbar(aes(xmin = CI_Lower_final, xmax = CI_Upper_final),
                width = 0.5, size = 0.5,
                position = position_dodge(width = 0.6),
                color = "darkslategray") +
  geom_point(position = position_dodge(width = 0.6), size = 3.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  # coord_flip() +
  labs(title = "Regional association between covariates and RISc",
       tag = "B)",
       x = "Effect Size", y = "") +
  scale_y_discrete(labels = covar_labels) +
  scale_color_paletteer_d("LaCroixColoR::Orange", direction = 1) +
  facet_grid(~region)+
  guides(color = guide_legend(reverse = T))+ #this reverses the order in the legend so it matches the order in the plot
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        title = element_text(size = 14),
        plot.title = element_text(size = 18))
region_mr_plot

#Final Figure 5 combined 
global_mr_plot/region_mr_plot
     
     






  
  
  
  
  
  