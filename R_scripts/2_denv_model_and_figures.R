setwd("/Users/abbeyporzucek/Library/Mobile Documents/com~apple~CloudDocs/Yale/Summer 2024/DENV Aim")

library(tidyverse)
library(INLA)
library(vroom)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(patchwork)
library(plotly)
library(janitor)
library(colorspace)
library(scales)
library(ggstatsplot)

#Poisson Regression-------------------------------------------------------------
#Note: you can read in the output from the regression below if you don't want to run the model yourself
#read in saved CSV file if necessary
denv_offset_final <- read.csv("denv_offset_final_updated_regions_11.21.2025.csv") 

#Set up the priors and the formula 
#Priors for RW2 model
hyper3.rw = list(prec = list(prior='pc.prec', param=c(1, 0.01))) # weaker (suggested INLA default)

#Define model formula
form1 <- as.formula('cases ~ 1 +
                          f(ISO_A0, model = "iid") + 
                          f(YearN, model = "iid")  +
                          f(YearN_continent, model = "iid")  +
                           f(YearN_country, model = "iid")  +
                          f(YearN2, model = "rw2",scale.model=TRUE, constr=T, hyper = hyper3.rw)+
                          f(YearN3, model = "rw2",  replicate=continentID,scale.model=TRUE, constr=T, hyper = hyper3.rw)+
                          f(YearN4, model = "rw2",  replicate=countryID,scale.model=TRUE, constr=T, hyper = hyper3.rw)
  ')

#Run the model
gam_offset <- inla(form1, data = denv_offset_final,  family = "poisson", E= offset1,
                   control.compute = list(dic = T,
                                          waic = T,
                                          config = T,
                                          return.marginals=F
                   ),
                   control.predictor=list(link=1),
                   
)    

summary(gam_offset)
#save output
# saveRDS(gam_offset, "poisson_updated_region_11.21.25.RDS")

#Combine GAM output with the original data set 
#get fitted values (these are the expected case counts )
gam_offset_output <- gam_offset$summary.fitted.values %>%
  cbind.data.frame(., denv_offset_final)
# write.csv(gam_offset_output, "poisson_model_output_11.21.2025.csv")

#Calculate RISc----------------------------------------------------------
#Goal is to add time-varying random effects together to get Relative Intensity Score (RISc)
#Note: RISc is coded as "dis" (due to naming of the variable that was updated later on)

#If necessary, read in data
gam_offset <- readRDS("poisson_updated_region_11.21.25.RDS")
gam_offset_output <- read.csv("poisson_model_output_11.21.2025.csv")

#rename columns in gam_output for clarity (columns getting renamed are fitted values)
offset_output_temp <- gam_offset_output %>% 
  rename(mean_fitted = mean,
         sd_fitted = sd,
         x0.025quant_fitted = 'X0.025quant', 
         x0.5quant_fitted = 'X0.5quant',
         x0.975quant_fitted = 'X0.975quant',
         mode_fitted = mode)

#Country/year random effects
#add in a column to join yearN_country data (so one value per observation in the data set)
offset_output_temp$join_yearN_country = paste(offset_output_temp$ISO_A0,offset_output_temp$YearN, sep = "_")

#create new data set with the yearN_country values from the model
yearN_country_random <- gam_offset$summary.random$YearN_country %>%
  rename(mean_ync = mean,
         sd_ync = sd,
         x0.025_ync = '0.025quant',
         x0.5_ync = '0.5quant',
         x0.975_ync = '0.975quant') %>%
  select(-mode, -kld) %>%
  mutate(x0.05_ync = mean_ync - (1.645 * sd_ync), #calculate 5th percentile
         x0.95_ync = mean_ync + (1.645 * sd_ync)) #calculate 95th percentile 



#Create a new dataset to hold all random effect values 
#do this by joining gam output with yearN_country_random
gam_rand_eff <- merge(offset_output_temp, yearN_country_random, 
                      by.x = "join_yearN_country", by.y = "ID")

#Yearly random effects 
#have one row for each YearN
year_random <- gam_offset$summary.random$YearN %>%
  rename(mean_year_random = mean,
         sd_year_random = sd,
         x0.025quant_year_random = '0.025quant',
         x0.5quant_year_random = '0.5quant',
         x0.975quant_year_random = '0.975quant') %>%
  select(-mode, -kld)%>%
  mutate(x0.05_year_random = mean_year_random - (1.645 * sd_year_random), #calculate 5th percentile
         x0.95_year_random = mean_year_random + (1.645 * sd_year_random)) #calculate 95th percentile 

#merge with inla_random_v2 data set
gam_rand_eff <- merge(gam_rand_eff, year_random , by.x="YearN", by.y="ID")

#add in continent level random effects 
cont_random <- gam_offset$summary.random$YearN_continent %>%
  rename(mean_cont_random = mean,
         sd_cont_random = sd,
         x0.025cont_random = '0.025quant',
         x0.5cont_random = '0.5quant',
         x0.975cont_random ='0.975quant') %>%
  select(-mode, -kld)
#merge with larger data set
gam_rand_eff <- merge(gam_rand_eff, cont_random, by.x = "YearN_continent", by.y = "ID")

#add in country level random effect
country_random <- gam_offset$summary.random$ISO_A0 %>%
  rename(mean_country_random = mean) %>%
  select(-sd, - '0.025quant', -'0.5quant', -'0.975quant', -mode, -kld)
#merge to larger data set
gam_rand_eff <- merge(gam_rand_eff, country_random, by.x = "ISO_A0", by.y = "ID")

#Add random effects together!
#create new data set for total random effects
offset_dis <- gam_rand_eff %>%
  mutate(dis_mean = (mean_ync +  mean_year_random + mean_cont_random),
         dis_variance = (sd_ync^2 +  sd_year_random^2 + sd_cont_random^2),
         dis_0.025 = (x0.025_ync + x0.025quant_year_random),
         dis_0.975 = (x0.975_ync +  x0.975quant_year_random),
         dis_0.05 = (x0.05_ync + x0.05_year_random),
         dis_0.95 = (x0.95_ync + x0.95_year_random),
         dis_all = (mean_ync + mean_year_random + mean_cont_random + mean_country_random),
         dis_time = (mean_ync + mean_year_random + mean_cont_random),
         continent = factor(continent, levels = c("Asia", "Oceania", "Caribbean","Central America", "South America"))
  ) %>%
  select("cases", "adjusted_year", "YearN", "ISO_A0", "full_name", "continent","population_total","mean_fitted", "dis_mean", "dis_variance", "dis_0.025", 
         "dis_0.975", "dis_0.05", "dis_0.95", "dis_time", "dis_all", "continentID", "countryID")
#Save final data set!
# write.csv(offset_dis, "final_risc_11.21.2025.csv")


#Figures 1 & 3 (heat maps)------------------------------------------------------
#Read in data set if needed
offset_dis <- read.csv("final_risc_11.21.2025.csv")
  #Factor continents for better visual organization 
  mutate(continent = factor(continent, levels = c("Asia", "Oceania","Caribbean" , "Central America", "South America")))

#Overwrite long names to make plots easier to read
offset_dis$full_name <- replace(offset_dis$full_name, offset_dis$full_name == "Lao People's Democratic Republic", "Laos")
offset_dis$full_name <- replace(offset_dis$full_name, offset_dis$full_name == "Saint Vincent And The Grenadines", "Saint Vincent")
offset_dis$full_name <- replace(offset_dis$full_name, offset_dis$full_name == "Micronesia (Federated States Of)", "Micronesia")


#make list to hold the heatmaps 
offset_heatmap <- list()

#Figure 3
  #Figure 3a
    offset_heatmap$heatmap <- ggplot() +
      geom_tile(data = offset_dis %>% filter(), 
                aes(x = adjusted_year, y = factor(full_name, levels = rev(unique(full_name))), fill = dis_mean)) +
      labs(title = "Dengue Intensity by Year and Location",
           tag = "A)",
           x = "Dengue Season (starting year)",
           y = "")+
      scale_fill_viridis_c(option = "viridis",
                           name = "Relative Intensity Score",
                           direction = 1,
                           guide = guide_colorbar(title.position = "top")) +
      theme_bw()+
      theme(strip.text = element_text(size = rel(0.85)),
            legend.position = "bottom",
            axis.text.y = element_text(size = 8.25),
            plot.title = element_text(size = 10))+
      facet_grid(continent ~ ., scales = "free_y", space = "free_y") 
    
    offset_heatmap$heatmap
  
  #Figure 3b - intensity categories
    #first make data frame with intensity cutoffs 
    offset_cutoff <- offset_dis %>%
      mutate(intensity = case_when(
        dis_mean <= -0.5 ~ "low",
        dis_mean >= 0.5 ~ "high",
        dis_mean > -0.5 & dis_mean < 0.5 ~ "expected"
      ))
    
    offset_heatmap$cutoff <- ggplot() +
      geom_tile(data = offset_cutoff, 
                aes(x = adjusted_year, y = factor(full_name, levels = rev(unique(full_name))), fill = intensity)) +
      labs(title = "Intensity Categories by Year and Location",
           tag = "B)",
           x = "Dengue Season (starting year)",
           y = "")+
      scale_fill_manual(values = c("expected" = "grey", "high" = "indianred", "low" = "cadetblue"),
                        name = "Intensity")+
      theme_bw()+
      theme(strip.text = element_text(size = rel(0.85)),
            legend.position = "bottom",
            axis.text.y = element_text(size = 8.25),
            plot.title = element_text(size = 10))+
      facet_grid(continent ~ ., scales = "free_y", space = "free_y")
    
    offset_heatmap$cutoff
  #Figure 3a and 3b
    offset_heatmap$heatmap + offset_heatmap$cutoff

#Figure 1
  #Make a heat map without the heat - just to show where data is available and where it isn't
  #Read in original cleaned data set, if necessary
    denv_offset_final <- read.csv("denv_offset_final_updated_regions_11.21.2025.csv")
  #Create case_status and border_color in the data set we used to make the poisson regression
  denv_offset_final$case_status <- ifelse(is.na(denv_offset_final$cases), "No", "Yes") #flag missing data
  denv_offset_final$border_color <- ifelse(is.na(denv_offset_final$cases), NA, "cadetblue4")
  #order the continents correctly
  denv_offset_final <- denv_offset_final %>%
    mutate(continent = factor(continent, levels = c("Asia", "Oceania","Caribbean" , "Central America", "South America")))
  
  #Overwrite long names
  denv_offset_final$full_name <- replace(denv_offset_final$full_name, denv_offset_final$full_name == "Lao People's Democratic Republic", "Laos")
  denv_offset_final$full_name <- replace(denv_offset_final$full_name, denv_offset_final$full_name == "Saint Vincent And The Grenadines", "Saint Vincent")
  denv_offset_final$full_name <- replace(denv_offset_final$full_name, denv_offset_final$full_name == "Micronesia (Federated States Of)", "Micronesia")
  
  #create a new column with status including 0s 
  denv_zero <- denv_offset_final %>%
    mutate(status_zero = case_when(
      cases == 0 ~ "Zero cases",
      cases > 0 ~ "Yes",
      is.na(cases) ~ "NA"
    ))
  
  #Figure 1a - heat map with data availability 
  offset_heatmap$no_heat_zero <- ggplot() +
    geom_tile(data = denv_zero, 
              aes(x = adjusted_year, y = factor(full_name, levels = rev(unique(full_name))), fill = status_zero, color = border_color), 
              size = 0.3) +
    labs(title = "Data Availability by Year and Location",
         tag = "A)",
         x= "Dengue Season (starting year)",
         # y = "Country")+
         y = "")+
    # ggtitle("Data Availability by Year and Location") +
    # xlab("Dengue Season (starting year)") + ylab("Country") +
    theme_bw() +
    theme(strip.text = element_text(size = rel(0.85)),
          legend.position = "bottom",
          axis.text.y = element_text(size = 8.25),
          plot.title = element_text(size = 10)) +
    scale_fill_manual(
      name = "Data available:",
      values = c("Yes" = "cadetblue", "NA" = "grey", "Zero cases" = "#C5EDD8"),
      labels = c("NA", "Yes", "Zero cases")
    ) +
    scale_color_identity() +
    facet_grid(continent ~ ., scales = "free_y", space = "free_y") 
  
  # theme(strip.text.y = element_blank())
  # Display plot
  offset_heatmap$no_heat_zero

  #Figure 1b - heat map of incidence rate
  #make new data set with incidence column
  incidence <- offset_cutoff %>%
    mutate(inc = (cases/population_total)*100000) %>%
    mutate(continent = factor(continent, levels = c("Asia", "Oceania","Caribbean" , "Central America", "South America")))
  
  # Plot with color outlines for cases
  offset_heatmap$inc_rate <- ggplot() +
    geom_tile(data = incidence,
              aes(x = adjusted_year, y = factor(full_name, levels = rev(unique(full_name))), fill = log1p(inc))) +
    labs(title = "Dengue Incidence by Year and Location",
         tag = "B)",
         x = "Dengue Season (starting year)",
         y = "")+
    scale_fill_viridis_c(option = "G",
                         name = "Dengue Incidence (log transformed)",
                         direction = -1,
                         guide = guide_colorbar(title.position = "top"),
                         na.value = "grey") +
    theme_bw()+
    theme(strip.text = element_text(size = rel(0.85)),
          legend.position = "bottom",
          axis.text.y = element_text(size = 8.25),
          plot.title = element_text(size = 10))+
    facet_grid(continent ~ ., scales = "free_y", space = "free_y")
  offset_heatmap$inc_rate

offset_heatmap$no_heat_zero + offset_heatmap$inc_rate

#save the offset_cutoff data frame
# write.csv(offset_cutoff, "risc_cutoff_11.21.2025.csv")

#Figure 2 (Nicaragua dot plots)-----------------------------------------------------------
  #Read in data sets, if necessary
    denv_offset_final <- read.csv("denv_offset_final_updated_regions_11.21.2025.csv")
    offset_dis <- read.csv("final_risc_11.21.2025.csv") %>%
      #Factor continents for better visual organization 
      mutate(continent = factor(continent, levels = c("Asia", "Oceania","Caribbean" , "Central America", "South America")))
    # offset_cutoff <- read.csv("risc_cutoff_11.21.2025.csv")
  #Make a list to hold the dot plots
    offset_dots <- list()
  #Recalculate the offset term
    offset_trend <- offset_dis %>%
      mutate(offset = population_total/100000)

  #Figure 2
   #Figure 2a
    nic_cases <- ggplot(data = denv_offset_final %>% filter(ISO_A0 == "NIC")) +
      geom_point(aes(x = adjusted_year, y = (cases)), size = 3) +
      theme_minimal()+
      theme(legend.position = "none")+
      labs(title = "Data input",
           subtitle = "Example data from Nicaragua",
           tag = "A)",
           x = "Dengue Season (starting year)",
           y = "Cases")
    nic_cases
    #Figure 2b
    offset_dots$nicaragua_dot_30 <- ggplot(data = offset_dis %>% filter(ISO_A0 == "NIC")) +
      geom_point(aes(x = adjusted_year, y = (cases), color = dis_mean), size = 3) +
      scale_color_viridis_c(option = "viridis",
                            name = "Relative Intensity Score",
                            guide = guide_colorbar(title.position = "top"))+
      geom_line(data = offset_trend %>% filter(ISO_A0 == "NIC"), 
                aes(x = adjusted_year, y = (mean_fitted * offset) / (exp(dis_mean)), 
                ), color = "black") +
      theme_minimal()+
      theme(legend.position = "bottom")+
      labs(title = "Quantifying deviations from baseline",
           tag = "B)",
           x = "Dengue Season (starting year)",
           y = "Cases")
    offset_dots$nicaragua_dot_30
    #Figure 2c
    offset_dots$nicaragua_cat <- ggplot() +
      geom_line(data = offset_trend %>% filter(ISO_A0 == "NIC"), 
                aes(x = adjusted_year, y = (mean_fitted * offset) / (exp(dis_mean)), 
                ), color = "black") +
      geom_point(data = offset_cutoff %>% filter(ISO_A0 == "NIC"),
                 aes(x = adjusted_year, y = cases, color = intensity), size = 3) +
      scale_color_brewer(palette = "Dark2")+
      theme_minimal()+
      theme(legend.position = "bottom")+
      labs(title = "Categorizing transmission intensity",
           tag = "C)",
           color = "Intensity",
           x = "Dengue Season (starting year)",
           y = "Cases")
    offset_dots$nicaragua_cat
    
    #Final Figure 2
    nic_cases/offset_dots$nicaragua_dot_30/offset_dots$nicaragua_cat

#Figure 4 (bar plots)-----------------------------------------------------------
#Figure 4a - global DIS over time
    #take mean dis of all countries 
    dis_stripe_data <- offset_dis %>%
      group_by(adjusted_year) %>%
      summarize(dis = mean(dis_mean))
    #Plot 4a - global
        stripe_mean2 <- ggplot(data = dis_stripe_data) +
          geom_col(aes(x = adjusted_year, y = dis, fill = dis)) +
          scale_fill_gradient2(
            low = "cadetblue", 
            mid = "#bfacb5", 
            high = "indianred", 
            midpoint = 0
          ) +
          labs(title = "Global Intensity Over Time",
               tag = "A)",
               x = "Dengue Season (starting year)",
               y = "Relative Intensity Score",
               fill = "RISc")+
          theme_minimal() +
          theme(axis.title.y = element_text(size = 11),
                axis.title.x = element_text(size = 11),
                plot.title = element_text(size = 15))
        stripe_mean2
    #Plot 4a - regional
      #stripe plot broken down by region 
      dis_stripe_cont_data <- offset_dis %>%
        group_by(continent, adjusted_year) %>%
        summarize(dis = mean(dis_mean))
      #plot the data
      stripe_mean_cont <- ggplot(data = dis_stripe_cont_data) +
        geom_col(aes(x = adjusted_year, y = dis, fill = dis)) +
        scale_fill_gradient2(
          low = "cadetblue", 
          mid = "#bfacb5", 
          high = "indianred", 
          midpoint = 0
        ) +
        ylim(-2.7, 2.7)+
        labs(title = "Regional Intensity Over Time",
             x = "Dengue Season (starting year)",
             y = "Relative Intensity Score",
             fill = "RISc")+
        facet_wrap(~ continent, scales = "free_y", ncol = 1) +
        theme_minimal()+
        theme(axis.title.y = element_text(size = 11),
              axis.title.x = element_text(size = 11),
              strip.text = element_text(size = 13),
              plot.title = element_text(size = 15))
      stripe_mean_cont
    #Figure 4a combined - plot the global and regional plots together
      stripe_plots <- stripe_mean2 / stripe_mean_cont
      stripe_plots_final <- stripe_plots + plot_layout(heights=c(1,5))
      stripe_plots_final

#Figure 4c - cases
  #Figure 4c part 1 - global cases
    #Make bar chart of outbreak categories 
    bar_prep <- offset_cutoff %>%
      select("adjusted_year", "cases", "intensity")
    bar_data <- bar_prep %>%
      group_by(adjusted_year, intensity) %>%
      summarise(cases = sum(cases, na.rm = TRUE)) %>%
      mutate(intensity = factor(intensity, levels = c("high", "expected", "low")))
    #Count how many rows in the original data fall into each Year + intensity combo
    intensity_counts <- bar_prep %>%
      group_by(adjusted_year, intensity) %>%
      summarise(intensity_count = n(), .groups = "drop")
    # Now join that count data to your bar_data summary
    bar_data <- bar_prep %>%
      filter(adjusted_year <= 2024) %>%
      group_by(adjusted_year, intensity) %>%
      summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
      left_join(intensity_counts, by = c("adjusted_year", "intensity"))
    #Make sure that intensity is being read in with levels
    bar_data <- bar_data %>%
      mutate(intensity = factor(intensity, levels = c("high", "expected", "low")))
    # Extract 5-year interval labels
    years_to_show <- unique(bar_data$adjusted_year)
    years_to_show <- years_to_show[years_to_show %% 10 == 0]
    #Plot figure 4c part 1 
    cases_intensity_scaled <- ggplot(bar_data, aes(x = factor(adjusted_year), y = cases / 1e6, fill = intensity)) +
      geom_bar(stat = "identity") +
      scale_x_discrete(breaks = as.character(years_to_show)) +
      scale_y_continuous(labels = label_number(suffix = "M", accuracy = 0.1)) +
      scale_fill_manual(values = c("expected" = "#bfacb5", "high" = "indianred", "low" = "cadetblue")) +
      labs(
        title = "Global Cases by Intensity",
        tag = "C)",
        x = "Dengue Season (starting year)",
        y = "Cases (Millions)",
        fill = "Intensity"
      ) +
      theme_minimal() +
      theme(axis.title.y = element_text(size = 11),
            axis.title.x = element_text(size = 11),
            plot.title = element_text(size = 15))
    cases_intensity_scaled
  #Figure 4c part 2 - regional cases
    #Make bar chart of outbreak categories 
    #Select relevant columns
    prep_cont <- offset_cutoff %>%
      select(continent, adjusted_year, cases, intensity)
    #Define the levels of intensity
    intensity_levels <- c("high", "expected", "low")
    #Summarize total cases
    bar_data_cont <- prep_cont %>%
      filter(adjusted_year <= 2023) %>%
      mutate(intensity = factor(intensity, levels = intensity_levels)) %>%
      group_by(continent, adjusted_year, intensity) %>%
      summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop") %>%
      # Step 4: Fill in missing combinations with 0 cases
      complete(continent, adjusted_year, intensity = intensity_levels, fill = list(cases = 0))
    #Count rows for each combo
    intensity_counts_cont <- prep_cont %>%
      filter(adjusted_year <= 2023) %>%
      mutate(intensity = factor(intensity, levels = intensity_levels)) %>%
      group_by(continent, adjusted_year, intensity) %>%
      summarise(intensity_count = n(), .groups = "drop") %>%
      complete(continent, adjusted_year, intensity = intensity_levels, fill = list(intensity_count = 0))
    #Join counts back to main data set
    bar_data_cont <- bar_data_cont %>%
      left_join(intensity_counts_cont, by = c("continent", "adjusted_year", "intensity"))
    #Make sure intensity is read in with levels 
    bar_data_cont <- bar_data_cont %>%
      mutate(intensity = factor(intensity, levels = c("high", "expected", "low")))
    #Plot 4c part 2
    cases_intensity_scaled_cont <- ggplot(bar_data_cont, aes(x = factor(adjusted_year), y = cases / 1e6, fill = intensity)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ continent, scales = "free_y", ncol = 1) +
      scale_y_continuous(labels = label_number(suffix= "M", accuracy = 0.01))+
      scale_x_discrete(breaks = as.character(years_to_show)) +
      scale_fill_manual(values = c("expected" = "#bfacb5", "high" = "indianred", "low" = "cadetblue")) +
      labs(
        title = "Regional Cases by Intensity",
        x = "Dengue Season (starting year)",
        y = "Cases (Millions)",
        fill = "Intensity"
      ) +
      theme_minimal() +
      theme(axis.title.y = element_text(size = 11),
            axis.title.x = element_text(size = 11),
            strip.text = element_text(size = 13),
            plot.title = element_text(size = 15))
    cases_intensity_scaled_cont
    #Figure 2c - combined 
    combined_cases_intensity <- cases_intensity_scaled / cases_intensity_scaled_cont +
      plot_layout(heights = c(1, 5))
    combined_cases_intensity
    
  #Figure 4b - proportions of country by intensity 
    # Plot intensity_count instead of cases
    #Plot 4b part 1 - global data
    count_intensity <- ggplot(bar_data, aes(x = factor(adjusted_year), y = intensity_count, fill = intensity)) +
      geom_bar(stat = "identity") +
      scale_x_discrete(breaks = as.character(years_to_show)) +
      scale_fill_manual(values = c("expected" = "#bfacb5", "high" = "indianred", "low" = "cadetblue"))+
      labs(
        title = "Global Intensity Proportion",
        tag = "B)",
        x = "Dengue Season (starting year)",
        y = "Number of Countries",
        fill = "Intensity"
      ) +
      theme_minimal() +
      theme(axis.title.y = element_text(size = 11),
            axis.title.x = element_text(size = 11),
            plot.title = element_text(size = 15))
    #Plot 4b part 2 - regional data
    count_intensity_cont <- ggplot(bar_data_cont, aes(x = factor(adjusted_year), y = intensity_count, fill = intensity)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ continent, scales = "free_y", ncol = 1) +
      scale_x_discrete(breaks = as.character(years_to_show)) +
      scale_fill_manual(values = c("expected" = "#bfacb5", "high" = "indianred", "low" = "cadetblue"))+
      labs(
        title = "Regional Intensity Proportion",
        x = "Dengue Season (starting year)",
        y = "Number of Countries",
        fill = "Intensity"
      ) +
      theme_minimal() +
      theme(axis.title.y = element_text(size = 11),
            axis.title.x = element_text(size = 11),
            strip.text = element_text(size = 13),
            plot.title = element_text(size = 15))
    count_intensity_cont
    #Plot 4b combined
    count_intensity_combined <- count_intensity / count_intensity_cont +  plot_layout(heights = c(1, 5))
    count_intensity_combined 
    
  #Final figure 4 - Plot all global/regional panels together!
    final_cc <- ((count_intensity_combined | combined_cases_intensity) + plot_layout(guides = 'collect'))
    final_fig_4 <- (stripe_plots_final  | final_cc) + plot_layout(widths = c(1, 2))
    final_fig_4
    
    #save plot
    # ggsave(final_fig_4, width = 14, height = 14, dpi = 300, filename = "final_fig_4_test.png")


#Figure S2-S6 (trend plots)-----------------------------------------------------
  #Make list to hold the plots 
    region_dots <- list()
  #Figure S2 - plots for all of Asia
    region_dots$asia <- ggplot(data = offset_dis %>% filter(continent == "Asia")) +
      geom_point(aes(x = adjusted_year, y = (cases), color = dis_mean), size = 2) +
      scale_color_viridis_c(option = "viridis",
                            name = "Relative Intensity Score",
                            guide = guide_colorbar(title.position = "top"))+
      geom_line(data = offset_trend %>% filter(continent == "Asia"), 
                aes(x = adjusted_year, y = (mean_fitted * offset) / (exp(dis_mean)), 
                ), color = "black") +
      facet_wrap(full_name ~ ., scales = "free_y", ncol = 3) + 
      theme_minimal()+
      theme(legend.position = "bottom")+
      ylab("Cases")+
      xlab("Dengue Season (starting year)")+
      ggtitle("Dengue in Asia")
    region_dots$asia
    
  #Figure S3 - plots for all of Oceania
    region_dots$oceania <- ggplot(data = offset_dis %>% filter(continent == "Oceania")) +
      geom_point(aes(x = adjusted_year, y = (cases), color = dis_mean), size = 2) +
      scale_color_viridis_c(option = "viridis",
                            name = "Relative Intensity Score",
                            guide = guide_colorbar(title.position = "top"))+
      geom_line(data = offset_trend %>% filter(continent == "Oceania"), 
                aes(x = adjusted_year, y = (mean_fitted * offset) / (exp(dis_mean)), 
                ), color = "black") +
      facet_wrap(full_name ~ ., scales = "free_y", ncol = 3) + #, scales = "free_y", space = "free_y") +
      theme_minimal()+
      theme(legend.position = "bottom")+
      ylab("Cases")+
      xlab("Dengue Season (starting year)")+
      ggtitle("Dengue in Oceania")
    region_dots$oceania
    
  #Figure S4 - plots for all of Caribbean
    region_dots$caribbean <- ggplot(data = offset_dis %>% filter(continent == "Caribbean")) +
      geom_point(aes(x = adjusted_year, y = (cases), color = dis_mean), size = 2) +
      scale_color_viridis_c(option = "viridis",
                            name = "Relative Intensity Score",
                            guide = guide_colorbar(title.position = "top"))+
      geom_line(data = offset_trend %>% filter(continent == "Caribbean"), 
                aes(x = adjusted_year, y = (mean_fitted * offset) / (exp(dis_mean)), 
                ), color = "black") +
      facet_wrap(full_name ~ ., scales = "free_y", ncol = 3) + #, scales = "free_y", space = "free_y") +
      theme_minimal()+
      theme(legend.position = "bottom")+
      ylab("Cases")+
      xlab("Dengue Season (starting year)")+
      ggtitle("Dengue in the Caribbean")
    region_dots$caribbean
    
    
  #Figure S5 - plots for all of Central America
    region_dots$ca <- ggplot(data = offset_dis %>% filter(continent == "Central America")) +
      geom_point(aes(x = adjusted_year, y = (cases), color = dis_mean), size = 2) +
      scale_color_viridis_c(option = "viridis",
                            name = "Relative Intensity Score",
                            guide = guide_colorbar(title.position = "top"))+
      geom_line(data = offset_trend %>% filter(continent == "Central America"), 
                aes(x = adjusted_year, y = (mean_fitted * offset) / (exp(dis_mean)), 
                ), color = "black") +
      facet_wrap(full_name ~ ., scales = "free_y", ncol = 3) + #, scales = "free_y", space = "free_y") +
      theme_minimal()+
      theme(legend.position = "bottom")+
      ylab("Cases")+
      xlab("Dengue Season (starting year)")+
      ggtitle("Dengue in Central America")
    region_dots$ca
    
  #Figure S6 - plots for all of South America
    region_dots$sa <- ggplot(data = offset_dis %>% filter(continent == "South America")) +
      geom_point(aes(x = adjusted_year, y = (cases), color = dis_mean), size = 2) +
      scale_color_viridis_c(option = "viridis",
                            name = "Relative Intensity Score",
                            guide = guide_colorbar(title.position = "top"))+
      geom_line(data = offset_trend %>% filter(continent == "South America"), 
                aes(x = adjusted_year, y = (mean_fitted * offset) / (exp(dis_mean)), 
                ), color = "black") +
      facet_wrap(full_name ~ ., scales = "free_y", ncol = 3) + #, scales = "free_y", space = "free_y") +
      theme_minimal()+
      theme(legend.position = "bottom")+
      ylab("Cases")+
      xlab("Dengue Season (starting year)")+
      ggtitle("Dengue in South America")
    region_dots$sa
    
    
    
    
#Figure S8 (imputation validation)----------------------------------------------
 #Read in data set, if necessary
  # denv_offset_final <- read.csv("denv_offset_final_updated_regions_11.21.2025.csv") 
  
  #Going to (1) find the pattern of missing data; (2) copy that pattern in places with out missing data; (3) validate after removing known data following pattern 
  #Figure out the average number of missing years per country
    total_missing_avg <- denv_offset_final %>%
      group_by(full_name, ISO_A0) %>% #group by country
      summarise(avg_na = mean(is.na(cases))) %>% #calculate average number of missing years
      filter(avg_na > 0) #select only places with missing data
    #Calculate average % of missing data
    mean(total_missing_avg$avg_na) #mean = 0.088 (0.88*34 ~3 years missing on average)
    
    #Figure out the average number of consecutive missing years
    #count years of consecutive missing data
    val <- denv_offset_final %>%
      group_by(ISO_A0) %>% #group by country
      arrange(adjusted_year, .by_group = TRUE) %>% #arrange the data by year 
      mutate(
        na_flag = ifelse(is.na(cases), 1, 0), #make a new column with a 1 value for each NA
        group = cumsum(lag(na_flag, default = 0) == 0 & na_flag==1), #group consecutive NAs
        missing = ifelse(na_flag == 1, sequence(rle(na_flag)$lengths), 0)) #create column with consecutive missing data
    #Summarize to find total consecutive missing years for each block of missing data
    missing_summary <- val %>%
      filter(na_flag == 1) %>% #this only includes observations with missing data
      group_by(full_name, ISO_A0, group) %>% #group by country and then block of consecutive missing data
      summarize(total_missing_years = max(missing), .groups = "drop") %>% #just get the maximum number of years for each consecutive block
      select(full_name, ISO_A0, total_missing_years) #select just the country and total missing years
    #Calculate the average # of consecutive missing years
    mean(missing_summary$total_missing_years) #mean = 2.15 years missing consecutively on average
    
    #Figure out the average number of countries with missing data per continent 
    # Summarize which places have missing data
    missing_data_table <- denv_offset_final %>%
      group_by(continent, full_name, ISO_A0) %>%
      mutate(missing = ifelse(is.na(cases), 1, 0)) %>% #flag observations with missing data
      summarise(total_missing = sum(missing)) #for each countries total the years of missing data 
    #Summarize the total missing by continent
    missing_cont <- missing_data_table %>%
      group_by(continent) %>%
      mutate(has_missing = ifelse(total_missing >= 1, 1, 0)) %>% #flag which countries have missing data
      summarise(total_countries = length(unique(ISO_A0)), #count the number of countries in each continent
                total_with_missing = sum(has_missing)) %>% #sum the number of countries with missing data
      filter(total_with_missing >= 1) #only include continents with missing data 
    #calculate the average number of countries with missing data per region 
    mean(missing_cont$total_with_missing) #on average each continent has 3.5 countries with missing data
    
    #Summary:
    #1. Countries that are missing data are missing ~3 years on average (9% of total data)
          #Going to round up to 4 to mesh better with #2
    #2. Missing data is missing in chunks of two years on average 
    #3. On average 4 countries per region are missing data
    #These will be the parameters to randomly drop data for the validation 
    
  #Now randomly remove values following the summary parameters above
    #Randomly choose 4 countries from central america and south america 
    ca_countries <- denv_offset_final %>%
      filter(continent == "Central America") %>%
      distinct(full_name) %>%
      sample_n(4) %>%
      pull(full_name)
    sa_countries <- denv_offset_final %>%
      filter(continent == "South America") %>%
      distinct(full_name) %>%
      sample_n(4) %>%
      pull(full_name)
    
    #Group the observations into pairs of two
    pairs <- denv_offset_final %>%
      group_by(ISO_A0) %>%
      arrange(adjusted_year)%>%
      mutate(year_pair = ceiling(row_number()/2))
    #Identify groups of years to remove data from
    years_to_remove <- pairs %>%
      filter(full_name %in% ca_countries | full_name %in% sa_countries) %>% #only using randomly selected countries
      distinct(ISO_A0, year_pair) %>% #get unique pairs of countries and years
      group_by(ISO_A0) %>%
      slice_sample(n=2) %>% #randomly select 2 of our distinct pairs in each country
      mutate(remove_cases = TRUE) #add in a flag to remove case data from selected pairs
    #Add the remove_cases flag back to full data set and replace cases with NA
    validation_data <- pairs %>%
      left_join(years_to_remove, by = c("ISO_A0", "year_pair")) %>%
      mutate(cases_val = ifelse(remove_cases == TRUE & !is.na(remove_cases), NA, cases))
    
    #Now rerun the model with validation_data and cases_val
    #Priors for RW2 model
    hyper3.rw = list(prec = list(prior='pc.prec', param=c(1, 0.01))) # weaker (suggested INLA default)
    
    #Formula
    #formula
    form1 <- as.formula('cases_val ~ 1 +
                          f(ISO_A0, model = "iid") + 
                          f(YearN, model = "iid")  +
                          f(YearN_continent, model = "iid")  +
                           f(YearN_country, model = "iid")  +
                          f(YearN2, model = "rw2",scale.model=TRUE, constr=T, hyper = hyper3.rw)+
                          f(YearN3, model = "rw2",  replicate=continentID,scale.model=TRUE, constr=T, hyper = hyper3.rw)+
                          f(YearN4, model = "rw2",  replicate=countryID,scale.model=TRUE, constr=T, hyper = hyper3.rw)
  ')
    
    #Run the GAM
    gam_validation <- inla(form1, data = validation_data,  family = "poisson", E= offset1,
                           control.compute = list(dic = T,
                                                  waic = T,
                                                  config = T,
                                                  return.marginals=F
                           ),
                           control.predictor=list(link=1),
                           
    )    
    
    summary(gam_validation)

    # Combine GAM output with the original data set
    #get fitted values (these are the expected case counts )
    gam_validation_output <- gam_validation$summary.fitted.values %>%
      cbind.data.frame(., validation_data)
    
    # write.csv(gam_validation_output, "validation_output_11.25.2025.csv")
    
    #If necessary, read in data set 
    gam_validation_output <- read.csv("validation_output_11.25.2025.csv")
    
    #Compare imputed values to original values 
    val_results <- gam_validation_output %>%
      filter(remove_cases == TRUE) %>%
      select(mean, cases, adjusted_year, full_name, ISO_A0, continent, offset1) %>%
      mutate(imputed = mean*offset1) #undo the offset
    
    cor.test(log1p(val_results$imputed), log1p(val_results$cases), method = "pearson")
  
  #Figure S8b - correlation plot
    val_plot <- ggscatterstats(data = val_results, x = cases, y=imputed, bf.message = F, marginal = F) +
      labs(x = "Reported Cases", y = "Imputed Values", tag = "B)") + 
      scale_x_continuous(transform = "log1p") +
      scale_y_continuous(transform = "log1p")+
      theme(axis.text.x = element_text(angle = 90))
    val_plot
    
  #Figure S8a - heat map of removed data
    #Factor continents so they show up in the right order
    denv_zero <- denv_zero %>%
      mutate(continent = factor(continent, levels = c("Asia", "Oceania","Caribbean" , "Central America", "South America")))
    
    #Plot S8a - heat map showing where data was removed from 
    validation_heatmap <- ggplot() +
      geom_tile(data = denv_zero, 
                aes(x = adjusted_year, y = factor(full_name, levels = rev(unique(full_name))), fill = status_zero, color = border_color), 
                size = 0.3) +
      geom_tile(data = val_results, 
                aes(x = adjusted_year, y = full_name, fill = "Removed data"), color = "#2d2b59") +
      labs(tag = "A)",
           x = "Dengue Season (starting year)",
           y = "") +
      theme(strip.text = element_text(size = rel(0.825))) +
      theme(legend.position = "bottom",
            axis.text.y = element_text(size = 8.25),
            plot.title = element_text(size = 10)) +
      scale_fill_manual(
        name = "Data available:",
        values = c("cadetblue", "grey", "#C5EDD8", "#2d2b59"),
        labels = c("Yes", "NA", "Zero cases", "Removed data"),
        breaks = c("Yes", "NA", "Zero cases", "Removed data")
      ) +
      scale_color_identity() +
      facet_grid(continent ~ ., scales = "free_y", space = "free_y")
    validation_heatmap    
    
    #Figure S8 combined
    (validation_heatmap + val_plot) + plot_layout(widths = c(1,1.5))
    
    
#Figure S9 (RISc histograms)----------------------------------------------------
  #Figure S9 part 1
    global_hist <- ggplot(data = offset_dis) +
      geom_histogram(aes(dis_mean), binwidth = 0.15) +
      labs(x = "Relative Intensity Score (RISc)",
           y = "Count",
           title = "Global RISc Histogram") +
      theme_bw()
    global_hist
  #Figure S9 part 2
    regional_hist <- ggplot(data = offset_dis) +
      geom_histogram(aes(dis_mean), binwidth = 0.15) +
      labs(x = "Relative Intensity Score (RISc)",
           y = "Count",
           title = "Regional RISc Histograms")+
      facet_wrap(~continent, nrow = 1) +
      theme_bw()
    regional_hist
  #Figure S9 combined
    global_hist/regional_hist
     
    
