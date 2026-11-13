# PREEMPT

# PREEMPT Virome Analysis

This repository contains the data and scripts for the PREEMPT mosquito virome project. The primary analysis is conducted using R, focusing on viral diversity, geographical distribution, host-virus networks, and ecological drivers of mosquito-borne viruses across multiple global sites (Cambodia, Central African Republic, French Guiana, and Madagascar).

## 📋 Table of Contents
- [Overview](#overview)
- [Requirements & Dependencies](#requirements--dependencies)
- [Input Data](#input-data)
- [Code Guide: `stats.R` Breakdown](#code-guide-statsr-breakdown)
- [Usage](#usage)

## 🔬 Overview
The core script, `stats.R`, processes viral abundance and detection data (measured in Reads Per Million - RPM) to generate a comprehensive suite of ecological and virological analyses. It assesses how viral communities vary by mosquito host species, geographical location, year, and season.

## 📦 Requirements & Dependencies
To run the analysis, you will need **R** (version 4.0+ recommended) and the following packages:
*   **Data Wrangling:** `dplyr`, `tidyr`, `tibble`, `reshape2`, `readr`
*   **Visualization:** `ggplot2`, `ggrepel`, `ggtext`, `scatterpie`, `ggvenn`, `gridExtra`, `ComplexHeatmap`, `circlize`
*   **Spatial/Mapping:** `sf`, `ggmap`, `rnaturalearth`
*   **Ecology & Stats:** `vegan`, `brglm2`, `car`, `rstatix`, `multcompView`, `indicspecies`, `permute`
*   **Networks/Trees:** `igraph`, `ggraph`, `ggtree`

## 📂 Input Data
Ensure the following files are in your working directory before running the script:
*   `virus_data.csv`: Main dataset containing Sample Names, Mosquito Species, Sites, Countries, RPM, and Virus classifications.

---

## 🗺️ Code Guide: `stats.R` Breakdown

To make it easier to navigate the `stats.R` script, it is divided into distinct analytical modules. Here is a guide to what each section does and what it outputs:

### 1. Rarefaction (`#### RAREFACTION ####`)
*   **What it does:** Calculates and plots species accumulation curves using the `vegan` package to determine if the viral diversity has been adequately sampled across different mosquito species.
*   **Output:** `Virus_Rarefaction_Curves.pdf`

### 2. Geographical Mapping (`#### MAP MOSQ SP ####` & `#### MAP TOP 50...`)
*   **What it does:** Uses `sf` and `rnaturalearth` to map physical sampling locations. It generates pie charts overlaying the maps to show sampling effort (mosquito species collected) and the distribution of the top 50 most frequently detected viruses globally.
*   **Outputs:** `[Country]_sampling_map.pdf`, `[Country]_top50_virus_map.pdf`

### 3. Diversity Metrics (`#### ALPHA AND SHANNON ####`)
*   **What it does:** Filters viruses (n >= 3 detections) and calculates Alpha Richness (virus count) and Shannon Diversity (H'). Fits a robust linear model (with HC3 standard errors) to compare viral alpha richness and Shannon diversity across mosquito species, using Tukey HSD for pairwise post-hoc comparisons to account for heteroscedasticity.
*   **Outputs:** 
    *   `Mosquito_Alpha_FINAL.pdf` 
    *   `Mosquito_Shannon_FINAL.pdf`
    *   `Virus_Richness_Heatmap_Final.pdf`

### 4. Statistical Modeling (`#### GLM ####` & `#### MOSQ SPECIES PERMANOVA ####`)
*   **What it does:** 
    *   **GLM:** Fits a standard Poisson Generalized Linear Model to partition the variance in viral alpha richness across mosquito species, site, season, and collection year. The model is statistically validated for overdispersion, and the sequential Type I deviance explained by each factor is visualised as a proportional pie chart.
    *   **PERMANOVA:** Uses `adonis2` to test the significance of these factors on overall virome composition.
    *   **Indicator Species:** Uses `indicspecies` to find which specific viruses are strong indicators for specific mosquito hosts, sites, or seasons.
*   **Outputs:** `Viral_GLM_Pie_Chart.pdf`, `permanova.pdf`, `PERMANOVA_Results_Global.txt`, `Indicators_[Factor].pdf`, `Consolidated_Indicators_Optimized.pdf`

### 5. Community Structure (`#### PCOA by factor ####`)
*   **What it does:** Performs Principal Coordinate Analysis (PCoA) on a Jaccard distance matrix to visualize how similar viromes are across different metadata factors (Genus, Species, Country, Season).
*   **Outputs:** `Virome_PCoA_4Panel_Revision.pdf`

### 6. Network Analysis (`#### NETWORK ####` & `#### MADAGASCAR ####`)
*   **What it does:** Constructs a Jaccard-similarity network of viral communities. It also includes a subset analysis exclusively for Madagascar to create a bipartite virus-host interactome (showing vOTUs with broader and more specialized hosts).
*   **Outputs:** `virome_network.pdf`, `Madagascar_Bipartite.pdf`, and various Madagascar-specific charts (Seasonal Richness, Venn Diagram).

### 7. Taxonomic Summaries (`#### FIGURE STACKED BARPLOT (WHOLE DATASET) ####`)
*   **What it does:** Aggregates the data by viral kingdom (e.g., dsRNA, negRNA, posRNA) and generates stacked barplots broken down by Country, Mosquito Species, Site, and Year.
*   **Outputs:** `Global_Filtered_Viral_Kingdoms.pdf` (contains 4 plots).

### 8. Sampling Effort & Unclassifieds (`#### HEATMAP OF SAMPLING DATA ####` & `#### UNCLASSIFIEDS ####`)
*   **What it does:** 
    *   Creates a heatmap of sampling effort (mosquitoes collected) over time/seasons.
    *   Calculates the "Viral Discovery Ratio" (proportion of unclassified contigs relative to known hits) and visualizes amino acid identities of novel viruses.
*   **Outputs:** `Sampling_Effort_Heatmap.pdf`, `Barplot_unclassifieds.pdf`, `Unclassified_Identity_by_Species.pdf`

### 9. Top 50 RPM Heatmaps (`#### HEATMAP TOP 50 RPM BY SPECIES/SITE ####`)
*   **What it does:** Uses `ComplexHeatmap` to generate clustered heatmaps of the Log10 Mean RPM for the 50 most abundant viruses, clustered by Mosquito Host Species and by Sampling Site.
*   **Outputs:** `Virus_Host_Fixed_Heatmap.pdf`, `Virus_Site_Heatmap.pdf`

---
