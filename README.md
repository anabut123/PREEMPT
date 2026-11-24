# PREEMPT Virome Analysis

This repository contains the data and R analysis code for the PREEMPT mosquito virome project, supporting the manuscript *Spatiotemporal factors shape mosquito virome diversity and composition in a host- and virus-dependent manner* (NCOMMS-25-80983C). The analysis covers viral diversity, geographic distribution, host–virus networks, and ecological drivers of mosquito-associated viruses across Cambodia, the Central African Republic, French Guiana, and Madagascar.

## Repository contents

| Path | Description |
|---|---|
| `stats.R` | Main R analysis script, organised into modular sections |
| `virus_dataset/FINAL_DATASET.csv` | Complete dataset of 2,290 virus hallmark gene-encoding contigs with sample metadata and RPM values |
| `virus_dataset/Filtered_Viruses_Metadata.csv` | Filtered contig set (1,401 contigs) retained for comparative analyses |
| `virus_dataset/viral_contigs.fasta` | Nucleotide sequences of all identified viral contigs |
| `virus_dataset/filtered_viruses.fasta` | Nucleotide sequences of the filtered contig set |
| `COI/` | Mosquito *cox1* barcode sequences and host phylogeny |
| `alignments_trees/` | Virus superclade alignments and full-size phylogenetic trees |

## Requirements

- **R** (>= 4.2.0)
- **System libraries (Linux only):** `libgdal-dev`, `libgeos-dev`, `libproj-dev` (required for `sf` and `ggmap`)

## Usage

Find `virus_data.csv`, then run `stats.R` (or execute individual modules marked by `#### ... ####` headers):


## Module guide

1. **Rarefaction (`#### RAREFACTION ####`)**
Calculates species accumulation curves (`vegan::specaccum`, 50 random permutations) per mosquito species.
*Output:* `Virus_Rarefaction_Curves.pdf`
2. **Geographic mapping (`#### MAP MOSQ SP ####`, `#### MAP TOP 50 MOST FREQUENTLY DETECTED VIRUS SPECIES ####`)**
Maps sampling locations (`sf`, `rnaturalearth`) with pie charts showing sampling effort by mosquito species and the distribution of the 50 most frequently detected viruses.
*Outputs:* `[Country]_sampling_map.pdf`, `[Country]_top50_virus_map.pdf`
3. **Diversity metrics (`#### ALPHA AND SHANNON ####`)**
Computes viral alpha richness and RPM-based Shannon diversity on the filtered matrix (> 3 detections). Compares mosquito species using linear models with HC3 robust standard errors (`emmeans`, `vcovHC`) and Tukey-adjusted compact letter displays, plus a mean richness heatmap across sites.
*Outputs:* `Mosquito_Alpha_FINAL.pdf`, `Mosquito_Shannon_FINAL.pdf`, `Virus_Richness_Heatmap_Final.pdf`
4. **Poisson GLM (`#### GLM ####`)**
Fits a Poisson GLM (`alpha_richness ~ species_reorder + Site + Season + Collection_Y`), validates equidispersion (`AER::dispersiontest`), and partitions Type I deviance explained per predictor.
*Output:* `Viral_Richness_Pie_Chart.pdf`
5. **PERMANOVA & indicator species (`#### PERMANOVA ####`)**
Runs global hierarchical PERMANOVA (`vegan::adonis2`, Jaccard distance, 999 permutations) nesting species within genus and site within country, plus per-species PERMANOVA models. Identifies significant indicator viruses ($p \le 0.05$) across all factors (`indicspecies::multipatt`, `IndVal.g`) and groups results by viral family.
*Outputs:* `PERMANOVA_Results_Global.txt`, `PERMANOVA_Summary.csv`, `Indicator_Virus_Species_By_Family.csv`, `Consolidated_Indicators_By_Family.pdf`, `Indicators_Grouped_By_Family.pdf`
6. **Community ordination (`#### PCOA by factor ####`)**
Plots a four-panel Principal Coordinates Analysis (`cmdscale` on Jaccard distances) coloured and shaped by mosquito genus, species, country, and season.
*Output:* `Virome_PCoA_4Panel_Revision.pdf`
7. **Virome similarity network (`#### NETWORK ####`)**
Constructs an undirected Jaccard similarity network (`igraph`, edge threshold 0.25, Fruchterman–Reingold layout over 5,000 iterations).
*Output:* `virome_network_CLEAN_LABELS.pdf`
8. **Madagascar sub-analysis (`#### MADAGASCAR ####`)**
Builds a bipartite host–virus interactome (`ggraph`) alongside seasonal richness, site overlap (Venn diagram), richness boxplots, and a top-20 presence/absence heatmap.
*Outputs:* `Madagascar_Bipartite_No_Cutoff.pdf`, `Madagascar_Seasonal_Richness.pdf`, `Madagascar_Venn_Diagram.pdf`, `Madagascar_Richness_Boxplot.pdf`, `Madagascar_Virus_Heatmap.pdf`
9. **Taxonomic composition (`#### FIGURE STACKED BARPLOT (WHOLE DATASET) ####`)**
Maps filtered viral families to viral kingdoms (`updated_kingdom.csv`) and generates stacked proportion barplots by country, mosquito species, site, and year.
*Output:* `Global_Filtered_Viral_Kingdoms.pdf`
10. **Sampling effort & unclassified viruses (`#### HEATMAP OF SAMPLING DATA ####`, `#### UNCLASSIFIEDS ####`)**
Plots unique library counts across site, year, and season; quantifies the proportion of unclassified novel contigs per mosquito species and their amino acid identity distributions.
*Outputs:* `Sampling_Effort_Heatmap.pdf`, `Barplot_unclassifieds.pdf`, `Unclassified_Identity_by_Species.pdf`
11. **Top-50 abundance heatmaps (`#### HEATMAP TOP 50 RPM BY SPECIES ####`, `#### HEATMAP TOP 50 RPM BY SITE ####`)**
Selects the 50 most abundant vOTUs by cumulative RPM, applies log_{10}({mean RPM} + 1)$ transformation, and renders hierarchically clustered heatmaps (`ComplexHeatmap`) by mosquito species and by sampling site grouped by country.
*Outputs:* `Virus_Host_Fixed_Heatmap.pdf`, `Virus_Site_Heatmap.pdf`

## Data availability

Raw metatranscriptomic sequencing reads are deposited in the NCBI Sequence Read Archive under BioProject **PRJNA1501385**. Assembled hallmark gene-encoding viral contigs are available in GenBank under accession numbers **PZ761577–PZ763865**.

## Citation

Koh C, Frangeul L, Blanc H, Gibaud A, Henrion-Lacritick A, Chitale RA, Boyer S, Duchemin JB, Dussart P, Girod R, Grau N, Ngoagouni C, Butković A, Saleh MC. *Spatiotemporal factors shape mosquito virome diversity and composition in a host- and virus-dependent manner.* Manuscript NCOMMS-25-80983C.
