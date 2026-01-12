# **Variations in annual dengue intensities are explained by temperature anomalies**

This paper is now avaliabe as a preprint at [medrxiv](https://www.medrxiv.org/content/10.64898/2025.12.19.25342670v1). 

Abigail J. Porzucek1,2,\*, Rafael Lopes1,2,\*, Yi Ting Chew1,2, Ke Li1,2, Oliver J. Brady3, Joshua L. Warren2,4, Colin J. Carlson1,2, Daniel M. Weinberger1,2,\*\*, Nathan D. Grubaugh1,2,5,\*\* 

1 Department of Epidemiology of Microbial Diseases, Yale School of Public Health, New Haven, Connecticut, USA
2 Public Health Modeling Unit, Yale School of Public Health, Yale University; New Haven, Connecticut, USA
3 Department of Infectious Disease Epidemiology and Dynamics, Faculty of Epidemiology and Population Health, London School of Hygiene and Tropical Medicine, London, UK
4 Department of Biostatistics, Yale School of Public Health, Yale University; New Haven, Connecticut, USA
5 Department of Ecology and Evolutionary Biology, Yale University; New Haven, Connecticut, USA

\* Co-first authors

\** Co-senior authors

Correspondence: abbey.porzucek@yale.edu (AJP); nathan.grubaugh@yale.edu (NDG)

## **Data Avaliablity**
All findings are based on publically avaliable from [OpenDengue](https://opendengue.org/), the [World Bank](https://data.worldbank.org/), and the [Copernicus Cimate Change Data Store](https://cds.climate.copernicus.eu/). All raw data files analyzed and used in this research are avaliable in the '/Input_files' folder. 

## **Pipeline running order**
All code to reproduce the analysis presented in the manuscript are in 'R_scripts/' folder. At 2026-01-12, the pipeline running order is:
- 1_opendengue_data_cleaning.R
- 2_denv_model_and_figures.R
- 3_climate_data_cleaning.R
- 4_posterior_covar_db_FigS10.R
- 5_meta_regression_and_fig_5.R
- 6_ph_risc_FigS7.R

## **Cleaned Data**
Cleaned data is avaliable in the '/cleaned_files' folder. The files contain:
- "DENV surveillance and RISc data.csv": reported cases, population, RISc, and RISc category by location and year
- "DENV covariate data": everything in the file above plus additional data on all covariates used in the meta regression analysis 
