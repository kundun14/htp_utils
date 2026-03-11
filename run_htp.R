




########## FUNCTION CALLING

# setwd("D:/DOCS/CIENCIA DE SUELOS/R_SUELOS/QUINOA/QUINOA_GITHUB/HTP_PROCESING/")
# source('htp_utils.R')




###########
###########  MODELING
###########

### 01. HTP RM ZONAL PRODUCTS


htp_zonal <- extract_htp_zonal(
  multi_path   = "PROCESING/MULTI/", 
  dsm_path     = "PROCESING/DSM/", 
  border_path  = "PROCESING/BORDERS/",
  dtm_file     = "PROCESING/DTM/dtm.tif",
  save_rasters = TRUE, # FALSE TRUE
  output_path  = "PROCESING/STACKS_07_03_26V2/",
  group_var       = "COD",
  daps = c("62", "86", "93", "121", "128")
)


#### PLOTS

plot026 <- plot_htp(
  target_cod = "CQC-026", 
  target_bloq = "B1",
  daps = c("62", "86", "93", "121", "128"),
  plot_path = 'PROCESING/PLOTS/GEOPRODUCTS'
)


### 02. HTP PHENO FEATURES

subGeno <- "CQC-026"

htp_pheno <- extract_htp_pheno(
  df         = htp_zonal, 
  group_var  = "COD", 
  bloq_var   = "BLOQ",
  time_var   = "DAP", 
  plot       = TRUE,  ## FALSE TRUE ///  FITTED PLOTS
  genotypes  = subGeno,
  plot_path  = "PROCESING/PLOTS/SUBSET/"
)


htp_features <- htp_pheno$features
htp_fitting_metrics  <- htp_pheno$quality

### BOXPLOT 

subGeno <- c("CQC-003", "CQC-026", "CQC-034", "CQC-051","CQC-067")

boxplot(
  df_features      = htp_features, 
  htp_feature     = "F19", 
  genotypes        = subGeno,
  save_plot        = TRUE,
  plot_path = "PROCESING/PLOTS/BOXPLOTS/" # # 800*1000
)

##########
########## STATISTICAL MODELING
##########


traits <- read.xlsx('MODELING_MATRIX_BLOQ.xlsx', 
                       sheet = 'MODELING_MATRIX_BLOQ')


### 03. CORRELATION

cor_table <- htp_correlations(df = traits, 
                              htp_features = htp_features, 
                              trait = "PMSEM")


### 04. ANOVA

anova_results <- run_glmer_anova(traits,
                                 htp_features,
                                 'PMSEM',
                                 'COD',
                                 'BLOQ')

# print(anova_results)

### 05. REGRESION

reg_results <- ht_regression(traits, 
                             htp_features, 
                             trait = "PMSEM")


plot_comparison_grid(reg_results, plot_path = "PROCESING/PLOTS/REGRES/")


#### TODO
#### 04. FACTOR ANALYSIS




