# *htp_utils*: High-Throughput Phenotyping Utilities for Agronomic Trials

## 1. What is htp_utils?

`htp_utils` is a collection of R scripts designed to perform common high-throughput phenotyping tasks using time series of multispectral imagery in the context of agronomic experimental trials. The toolkit enables the fast and simple derivation of vegetation indices and structural canopy variables from UAV imagery.

The package focuses on extracting multispectral indicators such as canopy cover, canopy volume, excessive greenness index, and normalized difference vegetation index. In addition, it provides tools for modeling crop phenology using spline-based growth curves and for conducting statistical analysis through regression and generalized linear mixed models. Several plotting utilities are also included to visualize spatial variables, phenological curves, and statistical relationships.

---

## 2. Problem Statement

Phenotyping is the process of measuring and analyzing observable plant traits in order to understand their characteristics, performance, and genetic makeup (Bian et al., 2022). In modern agricultural research, remote sensing technologies such as UAV platforms and satellite imagery have become important tools for collecting phenotypic information at high spatial and temporal resolution (Feng et al., 2021; Rallo et al., 2020).

Typical phenotyping workflows consist of repeated multispectral UAV surveys during the crop growing season. From these observations, vegetation indices and canopy structural metrics can be derived to characterize plant health and development. By analyzing time-series observations of these variables, it becomes possible to model crop growth dynamics and extract phenological parameters.

Because high-throughput phenotyping is often applied in agricultural experimental trials, the unit of analysis should consider not only the genotype or variety under study but also the experimental block or replication structure of the design.

Despite the growing availability of remote sensing data, there are still relatively few simple workflows that allow researchers to easily derive phenotypic parameters from multispectral imagery and integrate them with agronomic traits such as yield or productivity indicators. The `htp_utils` toolkit aims to provide a practical workflow for extracting phenotypic variables, modeling crop growth, and linking these features with traits of agronomic interest.

---

# Modeling Workflow

## 3. Zonal Extraction of HTP Variables

### `extract_htp_zonal()`

This function processes multispectral reflectance data and digital surface models (DSM) for each flight time and produces a set of geospatial products to sucesive analysis. 
Canopy cover (CC), canopy volume (CV), excessive greenness index (ExG), and normalized difference vegetation index (NDVI) can be generated from RGB and multispectral orthomosaic images together with digital surface models (DSM).

The following equations are used to derive canopy cover and vegetation indices from RGB and multispectral imagery.

#### Canopy detection rule

$$
\text{Canopy} =
\left(\frac{\text{Red}}{\text{Green}} < P_1 \right)
\;\land\;
\left(\frac{\text{Blue}}{\text{Green}} < P_2 \right)
\;\land\;
\left(2 \cdot \text{Green} - \text{Red} - \text{Blue} > P_3 \right)
$$

#### Excess Green Index (ExG)

$$
ExG = 2g - r - b
$$

where the normalized RGB components are defined as

$$
r = \frac{\text{Red}}{\text{Red} + \text{Green} + \text{Blue}}
$$

$$
g = \frac{\text{Green}}{\text{Red} + \text{Green} + \text{Blue}}
$$

$$
b = \frac{\text{Blue}}{\text{Red} + \text{Green} + \text{Blue}}
$$

#### Normalized Difference Vegetation Index (NDVI)

$$
NDVI = \frac{\text{NIR} - \text{Red}}{\text{NIR} + \text{Red}}
$$


This function generated the zonal statistics of the geospatial products, also writes the
correspoding files in the local folder.

```r
htp_zonal <- extract_htp_zonal(
  multi_path   = "PROCESING/MULTI/",
  dsm_path     = "PROCESING/DSM/",
  border_path  = "PROCESING/BORDERS/",
  dtm_file     = "PROCESING/DTM/dtm.tif",
  save_rasters = TRUE,
  output_path  = "PROCESING/STACKS_07_03_26V2/",
  group_var    = "COD",
  daps = c("62", "86", "93", "121", "128")
)
```

---

The geospatial products genereted can be ploted  across the full or partial crop growth period. The example below generates plots for genotype **CQC-026** and block **B1** for several observation dates expressed as days after planting (DAP).

```r
plot026 <- plot_htp(
  target_cod  = "CQC-026",
  target_bloq = "B1",
  daps        = c("62", "86", "93", "121", "128"),
  plot_path   = "PROCESING/PLOTS/GEOPRODUCTS"
)
```

Here we show the canopy indicators for all blocks of genotype **CQC-026**

<!-- ```
fig / geoProducts_B1.png
fig / geoProducts_B2.png
fig / geoProducts_B3.png
fig / geoProducts_B4.png
``` -->

![Geospatial Products Batch 1](figs/geoProducst_B1.png)
<!-- ![Geospatial Products Batch 4](figs/geoProducst_B2.png)
![Geospatial Products Batch 4](figs/geoProducst_B3.png)
![Geospatial Products Batch 4](figs/geoProducst_B4.png) -->

---

## 6. Phenological Feature Extraction

### `extract_htp_pheno()`

This function allows to derive phenotypic features extracted from crop growth models derived from time-series UAV measurements fitted with spline functions.

Maximum canopy cover and canopy volume values can be obtained. Growth-rate curves are then used to derive additional parameters such as maximum growth rate, the day after planting at which this maximum occurs, and the duration of the half-maximum growth period.

Additional phenological metrics are derived from ExG and NDVI fitted curves  curves. For the ExG curve, slopes representing canopy expansion and senescence are obtained by fitting linear models to the increasing and decreasing canopy phases. Metrics such as maximum values, duration, and area under each phase of the curve are then calculated. Similarly, NDVI progression curves provide features including maximum NDVI, the DAP at which this maximum occurs, and the slopes representing vegetation increase and decline throughout the crop cycle.

The full list of features are shox in the following table:

| Variable name | UAS product | Feature description |
| :--- | :---: | :--- |
| F1 | **CC** | Maximum value of CC |
| F2 | | Maximum growth rate of CC |
| F3 | | DAP at maximum growth rate of CC |
| F4 | | Duration over the half maximum period of CC |
| F5 | **CV** | Maximum value of CV |
| F6 | | Maximum growth rate of CV |
| F7 | | DAP at maximum growth rate of CV |
| F8 | | Duration over the half maximum period of CV |
| F9 | **ExG** | Maximum of ExG value |
| F10 | | DAP at maximum of ExG value |
| F11 | | Increasing slope of ExG |
| F12 | | Maximum value of increasing line at maximum ExG DAP |
| F13 | | Duration of increasing ExG period |
| F14 | | Area of increasing period of ExG |
| F15 | | Decreasing slope of ExG |
| F16 | | Maximum value of decreasing line at maximum ExG DAP |
| F17 | | Duration of decreasing ExG period |
| F18 | | Area of decreasing period of ExG |
| F19 | **NDVI** | Maximum NDVI value |
| F20 | | DAP at maximum NDVI value |
| F21 | | Increasing slope of NDVI |
| F22 | | Decreasing slope of NDVI |


```r
subGeno <- "CQC-026"

htp_pheno <- extract_htp_pheno(
  df         = htp_zonal,
  group_var  = "COD",
  bloq_var   = "BLOQ",
  time_var   = "DAP",
  plot       = TRUE,
  genotypes  = subGeno,
  plot_path  = "PROCESING/PLOTS/SUBSET/"
)
```

The function also returns fitting quality metrics such as **R²** and **RMSE** for each block of the experimental design.

```r
htp_features <- htp_pheno$features
htp_fitting_metrics <- htp_pheno$quality
```

Generated figures:

<!-- ```
fig / Grid_Fit_NDVI_subset.png
fig / Grid_Fit_CC_subset.png
``` -->

![Geospatial Products Batch 1](figs/Grid_Fit_NDVI_subset.png)
![Geospatial Products Batch 1](figs/GGrid_Fit_CC_subset.png)

---

## 7. Boxplots of Phenotypic Traits

Boxplots can be generated to visualize the distribution of specific phenotypic features across selected genotypes.

```r
subGeno <- c("CQC-003", "CQC-026", "CQC-034", "CQC-051", "CQC-067")

boxplot(
  df_features = htp_features,
  htp_feature = "F19",
  genotypes   = subGeno,
  save_plot   = TRUE,
  plot_path   = "PROCESING/PLOTS/BOXPLOTS/"
)
```

Generated figures:

<!-- ```
fig / Boxplot_F3.png
fig / Boxplot_F11.png
``` -->
![Geospatial Products Batch 1](figs/Boxplot_F3.png)
![Geospatial Products Batch 1](figs/Boxplot_F11.png)


---

# Statistical Modeling

## 8. Trait Dataset

Agronomic traits can be loaded from an external spreadsheet that contains experimental observations per genotype and block.

```r
traits <- read.xlsx(
  "MODELING_MATRIX_BLOQ.xlsx",
  sheet = "MODELING_MATRIX_BLOQ"
)
```

---

## 9. Correlation Analysis

### `htp_correlations()`

This function generates correlation tables between the extracted phenological features and a selected agronomic trait.

Example:

```r
cor_table <- htp_correlations(
  df = traits,
  htp_features = htp_features,
  trait = "PMSEM"
)
```

---

## 10. ANOVA

### `run_glmer_anova()`

Genetic differences among genotypes can be evaluated using analysis of variance based on generalized linear mixed models.

Agronomic trials typically follow a randomized complete block design (RCBD). In this framework, the blocks are treated as random effects within a generalized linear mixed model with Gaussian error distribution and identity link function.

Example:

```r
anova_results <- run_glmer_anova(
  traits,
  htp_features,
  "PMSEM",
  "COD",
  "BLOQ"
)
```

---

## 11. Regression Modeling

### `ht_regression()`

This function develops predictive models using the set of twenty-two phenotypic features. Two automatic feature-selection strategies are implemented.

The first strategy applies forward feature selection, starting from a null model and progressively adding predictors. The second strategy applies backward feature elimination, beginning with a full model and removing non-significant variables.

Example:

```r
reg_results <- ht_regression(
  traits,
  htp_features,
  trait = "PMSEM"
)
```

Model comparisons and diagnostic plots can be generated with:

```r
plot_comparison_grid(
  reg_results,
  plot_path = "PROCESING/PLOTS/REGRES/"
)
```

Generated figures:

<!-- ```
fig / regresion_plots.png
``` -->
![Geospatial Products Batch 1](figs/regresion_plots.png)