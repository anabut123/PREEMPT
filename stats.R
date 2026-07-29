library(dplyr)
library(tidyr)
library(ggplot2)

# Load data - converting empty strings and 'nan' to proper NAs
df <- read.csv("virus_data.csv", na.strings = c("", "NA", "nan", " "), stringsAsFactors = FALSE)


#### RAREFACTION ####
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(tibble)

# 1. Filter for species with enough data (e.g., at least 3 samples)
species_to_analyze <- df %>%
  group_by(Mosq_Species_Name) %>%
  summarise(n_samples = n_distinct(Sample_Name)) %>%
  filter(n_samples >= 3) %>%
  pull(Mosq_Species_Name)

# 2. Create a function to calculate the accumulation curve for a single mosquito species
get_rarefaction_data <- function(sp_name) {
  temp_mat <- df %>%
    filter(Mosq_Species_Name == sp_name) %>%
    group_by(Sample_Name, Virus_name_fin_cut) %>%
    summarise(abundance = 1, .groups = 'drop') %>%
    pivot_wider(names_from = Virus_name_fin_cut, values_from = abundance, values_fill = 0) %>%
    column_to_rownames("Sample_Name")
  
  # Calculate accumulation curve (random method)
  curve <- specaccum(temp_mat, method = "random", permutations = 50)
  
  # Convert to a data frame for ggplot
  data.frame(
    Samples = curve$sites,
    Richness = curve$richness,
    SD = curve$sd,
    Mosq_Species_Name = sp_name
  )
}

# 3. Run the function for all selected species and combine results
rarefaction_results <- lapply(species_to_analyze, get_rarefaction_data) %>%
  bind_rows()

# 4. Plot the curves
ggplot(rarefaction_results, aes(x = Samples, y = Richness, color = Mosq_Species_Name)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD, fill = Mosq_Species_Name), alpha = 0.1, color = NA) +
  theme_minimal() +
  labs(
    title = "Rarefaction Analysis of Viral Discovery",
    subtitle = "Unique Viruses found per number of mosquito samples",
    x = "Number of Mosquito Samples",
    y = "Cumulative Virus Species Richness",
    fill = "Mosquito Species",
    color = "Mosquito Species"
  ) +
  theme(legend.position = "bottom")

# 5. Save as PDF
ggsave("Virus_Rarefaction_Curves.pdf", width = 10, height = 7)





#### MAP MOSQ SP ####
library(ggplot2)
library(sf)
library(ggmap)
library(dplyr)
library(rnaturalearth)
library(tidyr)
library(scatterpie)
library(ggrepel)


# 1. Load Data
data <- read.csv("virus_data.csv", 
                 header = TRUE, na.strings = c("", "NA", "nan"))

# Clean names
data <- data %>%
  mutate(
    Site = trimws(Site),
    Country = trimws(Country)
  )

# 2. Setup Manual Coordinates
manual_coords <- tribble(
  ~Site,            ~lon,    ~lat,
  "Bayanga",         16.20,   2.92,
  "Bouar",           15.58,   5.95,
  "Batalimo",        18.48,   3.68,
  "Gbozo",           18.15,   4.50, 
  "Mbaiki",          17.98,   3.87,
  "Pissa",           18.50,   4.10, 
  "Remire Montjoly", -52.27,  4.91,
  "Rattanakiri",    107.00,  13.73,
  "Siem Reap",      103.86,  13.36,
  "Koh Kong",       102.98,  11.61,
  "Ambato Boeny",    46.72,  -16.11,
  "Amparafaravola",  48.22,  -17.58,
  "Mampikony",       47.63,  -16.09,
  "Prek Toal",       103.40,  13.15, 
  "Macouria",        -52.47,   5.15, 
  "Roura",           -52.32,   4.60  
)

# Create lookup
unique_sites <- data %>%
  distinct(Site, Country) %>%
  left_join(manual_coords, by = "Site")

# 3. Species Counts (Sampling Effort)
mosquito_wide <- data %>%
  group_by(Site, Country, Mosq_Species_Name) %>%
  summarise(count = n_distinct(Sample_Name), .groups = "drop") %>% 
  pivot_wider(names_from = Mosq_Species_Name, values_from = count, values_fill = 0) %>%
  left_join(unique_sites, by = c("Site", "Country")) %>% 
  filter(!is.na(lon))

mosquito_wide <- mosquito_wide %>%
  mutate(total_samples = rowSums(
    dplyr::select(., where(is.numeric), -lon, -lat), 
    na.rm = TRUE
  ))

set.seed(42) 
mosquito_wide <- mosquito_wide %>%
  mutate(
    lon = lon + runif(n(), -0.25, 0.25),
    lat = lat + runif(n(), -0.25, 0.25)
  )

# 4. Map Data Setup
world <- ne_countries(scale = "medium", returnclass = "sf")
french_guiana_map <- ne_states(country = "France", returnclass = "sf") %>% 
  filter(name %in% c("Guyane française", "French Guiana", "Guyane"))

# 5. Define Colors
species_colors <- c(
  "Aedes aegypti" = "#E41A1C", "Aedes albopictus" = "#FF7F00",
  "Anopheles coustani" = "#984EA3", "Anopheles funestus" = "#4DAF4A",
  "Anopheles gambiae" = "#A65628", "Coquillettidia venezuelensis" = "#377EB8",
  "Culex antennatus" = "#FFFF33", "Culex bitaeniorhynchus" = "#FFC300",
  "Culex portesi" = "#BCBD22", "Culex quinquefasciatus" = "#999999",
  "Culex tritaeniorhynchus" = "#666666", "Culex vishnui" = "#808000",
  "Mansonia indiana" = "#1F78B4", "Mansonia titillans" = "#A6CEE3",
  "Mansonia uniformis" = "#00CED1"
)

# 6. Loop and Plot
unique_countries <- unique(mosquito_wide$Country)
max_effort <- max(mosquito_wide$total_samples, na.rm = TRUE)

for (country_name in unique_countries) {
  
  country_pie_data <- mosquito_wide %>% 
    filter(Country == country_name)
  
  if(nrow(country_pie_data) == 0) next 
  
  if (country_name %in% c("French Guiana", "Guyane", "French Guyana")) {
    country_map <- french_guiana_map
  } else {
    country_map <- world %>% filter(name == country_name | admin == country_name)
  }
  
  if(nrow(country_map) == 0) next 
  
  bbox <- st_bbox(country_map)
  pie_cols <- intersect(names(species_colors), names(country_pie_data))
  
  country_pie_data <- country_pie_data %>%
    mutate(scaled_r = 0.2) 
  
  p <- ggplot() +
    geom_sf(data = world, fill = "grey95", color = "grey80", linetype = "dashed", linewidth = 0.1) +
    geom_sf(data = country_map, fill = "white", color = "black", linewidth = 0.5) +
    geom_scatterpie(data = country_pie_data, 
                    aes(x = lon, y = lat, r = scaled_r), 
                    cols = pie_cols, color = "black", alpha = 0.9) +
    geom_text_repel(data = country_pie_data, 
                    aes(x = lon, y = lat, label = Site),
                    size = 4, 
                    fontface = "bold",
                    box.padding = 1.5, 
                    point.padding = 0.5,
                    segment.color = "grey30",
                    min.segment.length = 0) +
    scale_fill_manual(values = species_colors) +
    coord_sf(xlim = c(bbox["xmin"]-1.5, bbox["xmax"]+1.5), 
             ylim = c(bbox["ymin"]-1.5, bbox["ymax"]+1.5)) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "aliceblue"),
          legend.position = "bottom",
          legend.text = element_text(size = 8)) +
    labs(title = paste("Sampling Effort:", country_name),
         fill = "Mosquito Species")
  
  # Save Output
  file_name <- paste0(gsub(" ", "_", country_name), "_sampling_map.pdf")
  ggsave(file_name, p, width = 11, height = 8.5)
  print(paste("Generated map for:", country_name))
}


#### MAP TOP 50 MOST FREQUENTLY DETECTED VIRUS SPECIES ####
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
library(tidyr)
library(scatterpie)
library(ggrepel)
library(grDevices) 

# 1. Load Data
data <- read.csv("virus_data.csv", 
                 header = TRUE, na.strings = c("", "NA", "nan"))

# Clean names
data <- data %>%
  mutate(
    Site = trimws(Site),
    Country = trimws(Country),
    Virus_name_fin_cut = trimws(Virus_name_fin_cut)
  )

# 2. Identify Top 50 Viruses (Globally)
top_50_viruses <- data %>%
  count(Virus_name_fin_cut) %>%
  slice_max(n, n = 50, with_ties = FALSE) %>%
  pull(Virus_name_fin_cut)

# 3. Setup Manual Coordinates
manual_coords <- tribble(
  ~Site,              ~lon,    ~lat,
  "Bayanga",          16.20,   2.92,
  "Bouar",            15.58,   5.95,
  "Batalimo",         18.48,   3.68,
  "Gbozo",            18.15,   4.50,
  "Mbaiki",           17.98,   3.87,
  "Pissa",            18.50,   4.10,
  "Remire Montjoly", -52.27,   4.91,
  "Rattanakiri",     107.00,  13.73,
  "Siem Reap",       103.86,  13.36,
  "Koh Kong",        102.98,  11.61,
  "Ambato Boeny",     46.72,  -16.11,
  "Amparafaravola",   48.22,  -17.58,
  "Mampikony",        47.63,  -16.09,
  "Prek Toal",       103.40,  13.15,
  "Macouria",        -52.47,   5.15,
  "Roura",           -52.32,   4.60
)

# 4. Generate Virus Counts for Top 50
virus_wide_50 <- data %>%
  filter(Virus_name_fin_cut %in% top_50_viruses) %>%
  group_by(Site, Country, Virus_name_fin_cut) %>%
  summarise(count = n(), .groups = "drop") %>% 
  pivot_wider(names_from = Virus_name_fin_cut, values_from = count, values_fill = 0) %>%
  left_join(manual_coords, by = "Site") %>%
  mutate(total_hits = rowSums(dplyr::select(., any_of(top_50_viruses)), na.rm = TRUE))


# 5. Define Dynamic Colors
virus_colors_50 <- colorRampPalette(c("royalblue4", "seagreen3", "gold", "darkorange", "firebrick3", "darkorchid4", "deeppink4"))(50)
names(virus_colors_50) <- top_50_viruses

# 6. Map Data Setup
world <- ne_countries(scale = "medium", returnclass = "sf")
french_guiana_map <- ne_states(country = "France", returnclass = "sf") %>% 
  filter(name %in% c("Guyane française", "French Guiana", "Guyane"))

# 7. Loop and Plot
unique_countries <- unique(virus_wide_50$Country)

for (country_name in unique_countries) {
  
  country_pie_data <- virus_wide_50 %>% 
    filter(Country == country_name, !is.na(lon))
  
  if(nrow(country_pie_data) == 0) next 
  
  if (country_name %in% c("French Guiana", "Guyane", "French Guyana")) {
    country_map <- french_guiana_map
  } else {
    country_map <- world %>% filter(name == country_name | admin == country_name)
  }
  
  if(nrow(country_map) == 0) next 
  
  bbox <- st_bbox(country_map)
  
  present_in_country <- country_pie_data %>% 
    dplyr::select(any_of(top_50_viruses)) %>% 
    dplyr::select(where(~sum(.) > 0)) %>% 
    names()
  
  country_pie_data <- country_pie_data %>%
    mutate(scaled_r = 0.2)  
  
  p <- ggplot() +
    geom_sf(data = world, fill = "grey95", color = "grey80", linetype = "dashed", linewidth = 0.1) +
    geom_sf(data = country_map, fill = "white", color = "black", linewidth = 0.5) +
    geom_scatterpie(data = country_pie_data, 
                    aes(x = lon, y = lat, r = scaled_r), 
                    cols = present_in_country, color = "black", alpha = 0.9) +
    geom_text_repel(data = country_pie_data, 
                    aes(x = lon, y = lat, label = Site),
                    size = 4, fontface = "bold", box.padding = 1.5, min.segment.length = 0) +
    scale_fill_manual(values = virus_colors_50) +
    coord_sf(xlim = c(bbox["xmin"]-2, bbox["xmax"]+2), 
             ylim = c(bbox["ymin"]-2, bbox["ymax"]+2)) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "aliceblue"),
          legend.position = "bottom",
          legend.text = element_text(size = 6)) + 
    guides(fill = guide_legend(ncol = 5)) +
    labs(title = paste("Virus Distribution (Top 50):", country_name),
         subtitle = "All sites shown with uniform marker size",
         fill = "Virus Species")
  
  ggsave(paste0(gsub(" ", "_", country_name), "_top50_virus_map.pdf"), p, width = 12, height = 9)
  print(paste("Saved map for", country_name))
}



#### ALPHA AND SHANNON ####
library(dplyr)
library(tidyr)
library(vegan)
library(car)

raw_data <- read.csv("virus_data.csv") 

# We use RPM to capture the "evenness" of the viral community
aggregated_data <- raw_data %>%
  group_by(Sample_Name, Mosq_Species_Name, Site, Country, Collection_Y, Season, Virus_name_fin_cut) %>%
  summarise(total_rpm = sum(RPM, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Virus_name_fin_cut, values_from = total_rpm, values_fill = 0)

# Identify metadata vs virus columns
metadata_cols <- c("Sample_Name", "Mosq_Species_Name", "Site", "Country", "Collection_Y", "Season")
virus_cols <- setdiff(names(aggregated_data), metadata_cols)
virus_matrix_all <- as.matrix(aggregated_data[, virus_cols])

# --- Step B: FILTER VIRUSES (Occurrence >= 3) ---
virus_counts <- colSums(virus_matrix_all > 0)

# Identify viruses that meet the threshold (e.g., detected in at least 3 mosquitoes)
keep_viruses <- names(virus_counts[virus_counts >= 3])

# Create the filtered abundance matrix
virus_matrix_filtered <- virus_matrix_all[, keep_viruses]

# --- Step C: Calculate Metrics based on Filtered Data ---

# 1. Alpha Richness: Counts how many 'filtered' viruses are present (RPM > 0)
aggregated_data$alpha_richness <- rowSums(virus_matrix_filtered > 0)

# 2. Shannon Diversity: Uses the RPM values of the 'filtered' viruses
aggregated_data$shannon <- vegan::diversity(virus_matrix_filtered, index = "shannon")

# --- Step D: Filter Mosquito Groups (n >= 3) ---
sampled_data <- aggregated_data %>%
  group_by(Mosq_Species_Name) %>%
  filter(n() >= 3) %>%
  ungroup()

# --- Step E: Reorder and Factorize ---
species_alphabetical <- rev(sort(unique(sampled_data$Mosq_Species_Name)))
sampled_data$species_reorder <- factor(sampled_data$Mosq_Species_Name, levels = species_alphabetical)
sampled_data <- sampled_data %>% mutate(across(c(Site, Season), as.factor))

# Print summary to see the impact
print(paste("Original virus count:", length(virus_cols)))
print(paste("Filtered virus count (n >= 3):", length(keep_viruses)))




library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)
library(rstatix)
library(multcompView)

# ---------------------------------------------------------
# ALPHA RICHNESS SECTION (Non-Parametric)
# ---------------------------------------------------------

# 1. Non-parametric Dunn's Test for Alpha Richness
dunn_alpha <- sampled_data %>% dunn_test(alpha_richness ~ species_reorder, p.adjust.method = "BH")
p_vals_alpha <- dunn_alpha$p.adj
names(p_vals_alpha) <- paste(dunn_alpha$group1, dunn_alpha$group2, sep = "-")
cld_alpha_vec <- multcompLetters(p_vals_alpha)$Letters

cld_alpha <- data.frame(
  species_reorder = names(cld_alpha_vec),
  .group = as.character(cld_alpha_vec),
  stringsAsFactors = FALSE
)

# 2. Clean and Sync Order
cld_alpha$.group <- trimws(as.character(cld_alpha$.group))
plot_order <- levels(sampled_data$species_reorder)
final_alpha_stats <- cld_alpha[match(plot_order, cld_alpha$species_reorder), ]

# 3. Colour Mapping (Shared Letters = Shared Colours)
unique_groups_alpha <- unique(final_alpha_stats$.group)
n_grps_alpha <- length(unique_groups_alpha)
alpha_pal <- rainbow(n_grps_alpha, s = 0.5, v = 0.9)
names(alpha_pal) <- unique_groups_alpha
final_alpha_colors <- alpha_pal[final_alpha_stats$.group]

# 4. Generate PDF
pdf("Mosquito_Alpha_NonParametric_FINAL.pdf", width = 13, height = 8)
par(mar = c(5, 15, 4, 6)) 

boxplot(alpha_richness ~ species_reorder, data = sampled_data,
        horizontal = TRUE, las = 1, 
        col = final_alpha_colors, 
        border = "black",
        xlab = "Alpha Richness (Virus Count)", ylab = "",
        main = "Virus Alpha Richness per Mosquito Species",
        sub = "Kruskal-Wallis test with Dunn's post-hoc pairwise comparisons (BH adjusted)",
        pch = 16, outcol = "red")

x_pos_alpha <- max(sampled_data$alpha_richness, na.rm = TRUE) * 1.05
text(x = x_pos_alpha, y = 1:nrow(final_alpha_stats), 
     labels = final_alpha_stats$.group, pos = 4, xpd = TRUE, cex = 1, font = 2)
dev.off()


# ---------------------------------------------------------
# SHANNON DIVERSITY SECTION (Non-Parametric)
# ---------------------------------------------------------

# 1. Non-parametric Dunn's Test for Shannon Diversity
dunn_shannon <- sampled_data %>% dunn_test(shannon ~ species_reorder, p.adjust.method = "BH")
p_vals_shannon <- dunn_shannon$p.adj
names(p_vals_shannon) <- paste(dunn_shannon$group1, dunn_shannon$group2, sep = "-")
cld_shannon_vec <- multcompLetters(p_vals_shannon)$Letters

cld_shannon <- data.frame(
  species_reorder = names(cld_shannon_vec),
  .group = as.character(cld_shannon_vec),
  stringsAsFactors = FALSE
)

# 2. Clean and Sync Order
cld_shannon$.group <- trimws(as.character(cld_shannon$.group))
final_shannon_stats <- cld_shannon[match(plot_order, cld_shannon$species_reorder), ]

# 3. Colour Mapping
unique_groups_shannon <- unique(final_shannon_stats$.group)
n_grps_shannon <- length(unique_groups_shannon)
shannon_pal <- rainbow(n_grps_shannon, s = 0.5, v = 0.9)
names(shannon_pal) <- unique_groups_shannon
final_shannon_colors <- shannon_pal[final_shannon_stats$.group]

# 4. Generate PDF
pdf("Mosquito_Shannon_NonParametric_FINAL.pdf", width = 13, height = 8)
par(mar = c(5, 15, 4, 6)) 

boxplot(shannon ~ species_reorder, data = sampled_data,
        horizontal = TRUE, las = 1, 
        col = final_shannon_colors, 
        border = "black",
        xlab = "Shannon Index (H')", ylab = "",
        main = "Virus Shannon Diversity per Mosquito Species",
        sub = "Kruskal-Wallis test with Dunn's post-hoc pairwise comparisons (BH adjusted)",
        pch = 16, outcol = "red")

x_pos_shannon <- max(sampled_data$shannon, na.rm = TRUE) * 1.05
text(x = x_pos_shannon, y = 1:nrow(final_shannon_stats), 
     labels = final_shannon_stats$.group, pos = 4, xpd = TRUE, cex = 1, font = 2)
dev.off()

print("Both Alpha and Shannon PDFs have been generated with Non-Parametric statistics.")


# ---------------------------------------------------------
# MEAN ALPHA RICHNESS HEATMAP
# ---------------------------------------------------------

# 1. Aggregate Data
heatmap_data <- sampled_data %>%
  filter(alpha_richness > 0) %>%
  group_by(Country, Site, Mosq_Species_Name) %>%
  summarise(
    mean_richness = mean(alpha_richness, na.rm = TRUE),
    n_samples = n_distinct(Sample_Name), 
    .groups = 'drop'
  )

# 2. Identify Species found in more than one Country (for Bolding)
multi_country_species <- heatmap_data %>%
  group_by(Mosq_Species_Name) %>%
  summarise(n_countries = n_distinct(Country)) %>%
  filter(n_countries > 1) %>%
  pull(Mosq_Species_Name)

# 3. Define Genus-specific Colours
genus_colors <- c(
  "Mansonia"   = "forestgreen", 
  "Culex"      = "firebrick3", 
  "Anopheles"  = "darkorchid4", 
  "Aedes"      = "dodgerblue3"
)

# 4. Create Styled Labels for the Y-Axis
heatmap_data <- heatmap_data %>%
  mutate(
    Genus = sub(" .*", "", Mosq_Species_Name),
    Genus_Col = ifelse(Genus %in% names(genus_colors), genus_colors[Genus], "black"),
    label_style = case_when(
      Mosq_Species_Name %in% multi_country_species ~ 
        paste0("<span style='color:", Genus_Col, "'>***I. ", Mosq_Species_Name, "***</span>"),
      TRUE ~ paste0("<span style='color:", Genus_Col, "'>*I. ", Mosq_Species_Name, "*</span>")
    )
  )

# 5. Set Geographical Order for Sites
site_order <- c(
  "Koh Kong", "Prek Toal", "Rattanakiri", "Siem Reap",      # Cambodia
  "Batalimo", "Bayanga", "Bouar", "Gbozo", "Mbaiki", "Pissa", # CAR
  "Macouria", "Remire Montjoly", "Roura",                  # French Guiana
  "Ambato Boeny", "Amparafaravola", "Mampikony"            # Madagascar
)
heatmap_data$Site <- factor(heatmap_data$Site, levels = site_order)

# 6. Generate the Plot
final_plot <- ggplot(heatmap_data, aes(x = Site, y = reorder(label_style, Mosq_Species_Name))) +
  geom_tile(aes(fill = mean_richness), color = "grey90", linewidth = 0.1) +
  geom_text(aes(label = n_samples), color = "white", size = 3, fontface = "bold") +
  scale_fill_viridis_c(option = "viridis", name = "Mean virus\nrichness") +
  facet_grid(. ~ Country, scales = "free_x", space = "free_x", switch = "x") +
  theme_minimal() +
  labs(x = NULL, y = "Mosquito species") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_markdown(size = 10),
    strip.placement = "outside",
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(0, "lines"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# 7. Save output
ggsave("Virus_Richness_Heatmap_Final.pdf", plot = final_plot, width = 12, height = 9)
print("Plot saved as Virus_Richness_Heatmap_Final.pdf")






#### GLM ####
library(brglm2)
library(dplyr)
library(ggplot2)
library(ggrepel)

# 1. Fit Penalised Poisson GLM
glm_richness_penalised <- glm(
  alpha_richness ~ species_reorder + Site + Season + Collection_Y, 
  data = sampled_data, 
  family = poisson(link = "log"),
  method = "brglmFit"
)

# 2. Extract Type I Sequential Deviance Contributions
dev_table <- anova(glm_richness_penalised, test = "Chisq")
null_dev <- glm_richness_penalised$null.deviance

final_pie_data <- data.frame(
  Factor = c(rownames(dev_table)[-1], "Unexplained"),
  Percentage = (c(dev_table$Deviance[-1], glm_richness_penalised$deviance) / null_dev) * 100
) %>%
  filter(!is.na(Percentage)) %>%
  mutate(Factor = gsub("species_reorder", "Species", Factor)) %>%
  arrange(desc(Factor)) %>%
  mutate(ypos = cumsum(Percentage) - 0.5 * Percentage)

print(paste("Total Sum:", sum(final_pie_data$Percentage)))

# 3. Generate and Save Pie Chart
p <- ggplot(final_pie_data, aes(x = "", y = Percentage, fill = Factor)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Contribution of Factors to Viral Richness",
    subtitle = "Proportion of Deviance Explained (Penalised GLM)"
  ) +
  geom_text_repel(
    aes(y = ypos, label = paste0(round(Percentage, 1), "%")),
    size = 4,
    nudge_x = 0.8, 
    show.legend = FALSE
  )

ggsave(
  "Viral_GLM_Pie_Chart.pdf", 
  plot = p, 
  width = 8, 
  height = 6
)





#### MOSQ SPECIES PERMANOVA ####
library(dplyr)
library(tidyr)
library(vegan)
library(indicspecies)

# --- Step 1: Load and Process Data ---
raw_data <- read.csv("virus_data.csv")

# Include MorphoID_Genus in the aggregation to preserve the hierarchy
aggregated_data <- raw_data %>%
  group_by(Sample_Name, MorphoID_Genus, Mosq_Species_Name, Site, Country, Collection_Y, Season) %>%
  distinct(Virus_name_fin_cut) %>% 
  mutate(presence = 1) %>%
  pivot_wider(names_from = Virus_name_fin_cut, values_from = presence, values_fill = 0) %>%
  ungroup()

metadata_cols <- c("Sample_Name", "MorphoID_Genus", "Mosq_Species_Name", "Site", "Country", "Collection_Y", "Season")
virus_cols <- setdiff(names(aggregated_data), metadata_cols)

# Create the binary matrix for PERMANOVA
virus_matrix_full <- as.matrix(aggregated_data[, virus_cols])

# --- Step 2: Filtering and Factorization ---
# Filter for groups with n >= 3
sampled_data <- aggregated_data %>%
  group_by(Mosq_Species_Name) %>%
  filter(n() >= 3) %>%
  ungroup()

# Reorder species factors alphabetically for consistency
species_alphabetical <- rev(sort(unique(sampled_data$Mosq_Species_Name)))
sampled_data$species_reorder <- factor(sampled_data$Mosq_Species_Name, levels = species_alphabetical)

# Ensure categorical variables are factors
sampled_data <- sampled_data %>% 
  mutate(across(c(MorphoID_Genus, Site, Season, Collection_Y, Country), as.factor))

# --- Step 3: Global Hierarchical PERMANOVA & Matrix Filtering ---

# 1. First, create a temporary matrix from the filtered samples
temp_matrix <- as.matrix(sampled_data[, virus_cols])

# 2. Identify viruses present in at least 3 mosquitoes within this subset
virus_counts <- colSums(temp_matrix > 0)
keep_viruses <- names(virus_counts[virus_counts >= 3])

# 3. Create the final virus_species_matrix using only the kept viruses
virus_species_matrix <- temp_matrix[, keep_viruses]

print(paste("Initial total viruses:", length(virus_cols)))
print(paste("Viruses remaining after <3 filter:", length(keep_viruses)))
print(paste("Mosquito samples in analysis:", nrow(sampled_data)))

row_sums <- rowSums(virus_species_matrix)
empty_samples <- which(row_sums == 0)

# 2. If there are empty samples, remove them
if(length(empty_samples) > 0) {
  print(paste("Removing", length(empty_samples), "samples that have no viruses after filtering."))
  
  # Remove from the matrix
  virus_species_matrix <- virus_species_matrix[-empty_samples, ]
  
  # Remove from the metadata
  sampled_data <- sampled_data[-empty_samples, ]
}

# 3. Now run the PERMANOVA
perm_result <- adonis2(virus_species_matrix ~ MorphoID_Genus / species_reorder + Site + Season + Collection_Y, 
                       data = sampled_data, 
                       method = "jaccard", 
                       permutations = 999)

print("--- Hierarchical PERMANOVA Results ---")
print(perm_result)

# --- Step 4: Indicator Species Analysis ---
# Identifies which specific viruses are indicators for specific mosquito species
library(indicspecies)
library(permute)
library(parallel)

inv_sp <- multipatt(virus_species_matrix, 
                    sampled_data$species_reorder, 
                    func = "IndVal.g", 
                    duleg = TRUE,            
                    control = how(nperm = 999))

# Once it finishes, look at the specialists
summary(inv_sp)


level_names <- levels(sampled_data$species_reorder)

# 2. Build the data frame carefully
indic_df <- data.frame(
  Virus = rownames(inv_sp$sign),
  Group_Index = inv_sp$sign$index,
  Stat = inv_sp$sign$stat,
  P_val = inv_sp$sign$p.value
) %>%
  filter(Group_Index > 0) %>%
  mutate(Mosquito_Species = level_names[Group_Index])

# 3. Filter for significant viruses and take top 5
top_indicators <- indic_df %>%
  filter(P_val <= 0.05) %>%
  group_by(Mosquito_Species) %>%
  slice_max(order_by = Stat, n = 20, with_ties = FALSE) %>%
  ungroup()

print(nrow(top_indicators))

p <- ggplot(top_indicators, aes(x = factor(Mosquito_Species), y = reorder(Virus, Stat))) +
  geom_point(aes(size = Stat, color = factor(Mosquito_Species))) +
  # Use a robust color scale
  scale_color_viridis_d() + 
  theme_bw() +
  labs(
    title = "Top Viral Indicators per Mosquito Species",
    x = "Mosquito Host",
    y = "Virus Species",
    size = "IndVal Strength"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  ) +
  guides(color = "none")
ggsave("permanova.pdf", p, width = 20, height = 8)







#### REST OF FACTORS PERMANOVA ####
library(dplyr)
library(tidyr)
library(vegan)
library(indicspecies)
library(ggplot2)

# --- Step 1 & 2: Process and Filter ---
raw_data <- read.csv("virus_data.csv")

aggregated_data <- raw_data %>%
  group_by(Sample_Name, MorphoID_Genus, Mosq_Species_Name, Site, Country, Collection_Y, Season) %>%
  distinct(Virus_name_fin_cut) %>% 
  mutate(presence = 1) %>%
  pivot_wider(names_from = Virus_name_fin_cut, values_from = presence, values_fill = 0) %>%
  ungroup()

metadata_cols <- c("Sample_Name", "MorphoID_Genus", "Mosq_Species_Name", "Site", "Country", "Collection_Y", "Season")
virus_cols <- setdiff(names(aggregated_data), metadata_cols)

# Filter for species with n >= 3
sampled_data <- aggregated_data %>%
  group_by(Mosq_Species_Name) %>%
  filter(n() >= 3) %>%
  ungroup()

virus_species_matrix <- as.matrix(sampled_data[, virus_cols])
rownames(virus_species_matrix) <- sampled_data$Sample_Name

# 2. Filter viruses with < 3 occurrences
virus_counts <- colSums(virus_species_matrix > 0)
keep_viruses <- names(virus_counts[virus_counts >= 3])
virus_matrix_aligned <- virus_species_matrix[, keep_viruses]

# 3. Remove samples that now have zero viruses
virus_matrix_aligned <- virus_matrix_aligned[rowSums(virus_matrix_aligned) > 0, ]

# 4. Final Metadata Alignment
metadata_aligned <- sampled_data %>%
  filter(Sample_Name %in% rownames(virus_matrix_aligned)) %>%
  arrange(match(Sample_Name, rownames(virus_matrix_aligned))) %>%
  mutate(across(all_of(metadata_cols), ~factor(as.character(.x))))

# --- PERMANOVA for Global Drivers ---
perm_global <- adonis2(virus_matrix_aligned ~ MorphoID_Genus + Country / Site + Season, 
                       data = metadata_aligned, 
                       method = "jaccard", 
                       permutations = 999)

print("--- Statistical Drivers of the Virome ---")
print(perm_global)

# Save the PERMANOVA table to a text file
write.table(as.data.frame(perm_global), "PERMANOVA_Results_Global.txt", sep="\t")

# --- Step 4: Loop for Indicators ---
factors_to_test <- c("MorphoID_Genus", "Site", "Country", "Season")

for (f in factors_to_test) {
  # Skip factor if it only has 1 level after filtering
  if(length(levels(metadata_aligned[[f]])) < 2) next
  
  # Run Analysis
  inv <- multipatt(virus_matrix_aligned, metadata_aligned[[f]], 
                   func = "IndVal.g", duleg = TRUE, control = how(nperm = 999))
  
  # Process and Plot
  lvl_names <- levels(metadata_aligned[[f]])
  res_df <- data.frame(
    Virus = rownames(inv$sign),
    Group_Index = inv$sign$index,
    Stat = inv$sign$stat,
    P_val = inv$sign$p.value
  ) %>%
    filter(Group_Index > 0 & P_val <= 0.05) %>%
    mutate(GroupName = lvl_names[Group_Index])
  
  if(nrow(res_df) > 0) {
    p <- ggplot(res_df %>% group_by(GroupName) %>% slice_max(Stat, n=10), 
                aes(x = GroupName, y = reorder(Virus, Stat))) +
      geom_point(aes(size = Stat, color = GroupName)) +
      scale_color_viridis_d() +
      theme_bw() +
      labs(title = paste("Indicators for", f), x = f, y = "Virus") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0("Indicators_", f, ".pdf"), p, width = 12, height = 8)
  }
}

# --- Step 3.5: PCoA Plots ---
# 1. Calculate the Jaccard distance matrix
dist_matrix <- vegdist(virus_matrix_aligned, method = "jaccard")

# 2. Run the PCoA
pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  PC1 = pcoa_res$points[,1],
  PC2 = pcoa_res$points[,2],
  Sample_Name = rownames(virus_matrix_aligned)
) %>%
  left_join(metadata_aligned, by = "Sample_Name")

# Calculate variance explained for axis labels
pc1_var <- round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1)
pc2_var <- round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1)

# 3. Plot for SITE
p_site <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Site)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) + # Draws the "similarity cloud"
  theme_bw() +
  labs(title = "Virome Similarity by Site (PCoA)",
       x = paste0("PCoA1 (", pc1_var, "%)"),
       y = paste0("PCoA2 (", pc2_var, "%)"))

ggsave("PCoA_Site_Similarity.pdf", p_site, width = 8, height = 6)

# 4. Plot for SEASON
p_season <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Season)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  theme_bw() +
  labs(title = "Virome Similarity by Season (PCoA)",
       x = paste0("PCoA1 (", pc1_var, "%)"),
       y = paste0("PCoA2 (", pc2_var, "%)"))

ggsave("PCoA_Season_Similarity.pdf", p_season, width = 8, height = 6)



#### ALL PLOTS COMBINED ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(indicspecies)

# 1. Prepare a unified data frame for all factors
all_factors_list <- list()

# --- Process Mosquito Species ---
all_factors_list[["Species"]] <- top_indicators %>%
  dplyr::select(Virus, GroupName = Mosquito_Species, Stat, P_val) %>%
  dplyr::mutate(Factor = "Mosquito Species")

# --- Process Other Factors ---
factors_to_loop <- c("MorphoID_Genus", "Site", "Country", "Season")

for (f in factors_to_loop) {
  if(exists("metadata_aligned") && length(levels(metadata_aligned[[f]])) >= 2) {
    inv <- multipatt(virus_matrix_aligned, metadata_aligned[[f]], 
                     func = "IndVal.g", duleg = TRUE, control = how(nperm = 999))
    
    lvl_names <- levels(metadata_aligned[[f]])
    res_df <- data.frame(
      Virus = rownames(inv$sign),
      Group_Index = inv$sign$index,
      Stat = inv$sign$stat,
      P_val = inv$sign$p.value
    ) %>%
      dplyr::filter(Group_Index > 0 & P_val <= 0.05) %>%
      dplyr::mutate(GroupName = lvl_names[Group_Index],
                    Factor = f) %>%
      dplyr::group_by(GroupName) %>%
      dplyr::slice_max(Stat, n = 5) %>% 
      dplyr::ungroup()
    
    all_factors_list[[f]] <- res_df %>% 
      dplyr::select(Virus, GroupName, Stat, P_val, Factor)
  }
}

# 2. Combine all into one master data frame
master_plot_df <- dplyr::bind_rows(all_factors_list)

# 3. Create the Consolidated Plot
p_combined_final <- ggplot(master_plot_df, aes(x = GroupName, y = reorder(Virus, Stat))) +
  geom_point(aes(size = Stat, color = Factor), alpha = 0.8) +
  scale_size_continuous(range = c(3, 9)) + 
  facet_grid(. ~ Factor, scales = "free_x", space = "free_x") +
  scale_color_brewer(palette = "Set1") + 
  theme_bw() +
  labs(
    title = "Consolidated Viral Indicators",
    x = "Levels per Factor",
    y = "Virus Species",
    size = "IndVal Strength",
    color = "Grouping Category"
  ) +
  theme(
    axis.text.y = element_text(size = 11, color = "black"), 
    
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11, face = "italic"),
    
    # BOLD and visible Factor Names
    strip.text = element_text(face = "bold", size = 10), 
    strip.background = element_rect(fill = "gray90", color = "black"),
    
    # SHRINK SPACE BETWEEN PANELS
    panel.spacing.x = unit(0.1, "lines"), 
    
    legend.position = "bottom"
  )

# 4. SAVE WITH NEW DIMENSIONS:
ggsave("Consolidated_Indicators_Optimized.pdf", 
       p_combined_final, 
       width = 14,  
       height = 18) 


#### PCOA by factor ####
# Re-extract coordinates and add ALL metadata columns
pcoa_df <- data.frame(
  Axis1 = pcoa_res$points[,1],
  Axis2 = pcoa_res$points[,2],
  Genus = metadata_final$MorphoID_Genus,
  Species = metadata_final$Mosq_Species_Name,
  Country = metadata_final$Country,    # Added
  Season = metadata_final$Season,      # Added
  Site = metadata_final$Site           # Added
)

# Calculate % variation
pc1_var <- round(100 * pcoa_res$eig[1] / sum(pcoa_res$eig), 1)
pc2_var <- round(100 * pcoa_res$eig[2] / sum(pcoa_res$eig), 1)

library(gridExtra)

# Common theme and axis labels
p_theme <- theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))

x_lab <- paste0("Axis 1 (", pc1_var, "%)")
y_lab <- paste0("Axis 2 (", pc2_var, "%)")

# Plot 1: Genus and Country
p1 <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Genus, shape = Country)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Mosquito genus and Country", x = x_lab, y = y_lab) +
  p_theme

# Plot 2: Genus and Season
p2 <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Genus, shape = Season)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Mosquito genus and Season", x = x_lab, y = y_lab) +
  p_theme

# Plot 3: Species and Country
p3 <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Species, shape = Country)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Mosquito species and Country", x = x_lab, y = y_lab) +
  p_theme + theme(legend.position = "none") # Hide giant legend

# Plot 4: Species and Season
p4 <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Species, shape = Season)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Mosquito species and Season", x = x_lab, y = y_lab) +
  p_theme + theme(legend.position = "none")

# Final Arrangement
grid.arrange(p1, p2, p3, p4, ncol = 2)

library(gridExtra)

pdf("Virome_PCoA_4Panel_Revision.pdf", width = 12, height = 8)

# 2. Print the grid.arrange object
grid.arrange(p1, p2, p3, p4, ncol = 2)

dev.off()







#### 

#### NETWORK ####
library(igraph)
library(ggrepel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

# --- Step 1: Use the Filtered PERMANOVA Data ---
virus_cols <- keep_viruses

nodes_data <- sampled_data %>%
  mutate(Country = trimws(as.character(Country))) %>%
  group_by(Mosq_Species_Name, Collection_Y, Season, Site, Country, MorphoID_Genus) %>%
  summarize(across(all_of(virus_cols), ~ max(.)), .groups = "drop")

# Create Unique IDs
nodes_data$node_id <- paste(nodes_data$Mosq_Species_Name, nodes_data$Collection_Y, 
                            nodes_data$Season, nodes_data$Country, sep="_")
nodes_data$node_id <- make.unique(as.character(nodes_data$node_id), sep = "_")

# --- Step 2: Matrix & Similarity Calculation ---
virus_matrix <- as.matrix(nodes_data[, virus_cols])
rownames(virus_matrix) <- nodes_data$node_id

n_nodes <- nrow(virus_matrix)
similarity_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
rownames(similarity_matrix) <- rownames(virus_matrix)
colnames(similarity_matrix) <- rownames(virus_matrix)

for (i in 1:(n_nodes-1)) {
  for (j in (i+1):n_nodes) {
    p1 <- virus_matrix[i, ]
    p2 <- virus_matrix[j, ]
    inter <- sum(p1 & p2)
    uni <- sum(p1 | p2)
    if (uni > 0) {
      similarity <- inter / uni
      similarity_matrix[i,j] <- similarity
      similarity_matrix[j,i] <- similarity 
    }
  }
}

edges_raw <- as.data.frame(as.table(similarity_matrix)) %>%
  rename(from = Var1, to = Var2, weight = Freq) %>%
  mutate(from = as.character(from), to = as.character(to)) %>%
  filter(from < to & weight > 0)

country_map <- setNames(nodes_data$Country, nodes_data$node_id)
edges_raw$country_from <- country_map[edges_raw$from]
edges_raw$country_to   <- country_map[edges_raw$to]

edge_threshold <- 0.05 
strong_edges <- edges_raw %>%
  filter(weight >= edge_threshold | 
           country_from == "Guyane" | 
           country_to == "Guyane")

# --- Step 4: Create Network and Layout ---
vertex_data <- data.frame(
  name = nodes_data$node_id,
  species = nodes_data$Mosq_Species_Name,
  country = nodes_data$Country,
  year = as.factor(nodes_data$Collection_Y),
  season = nodes_data$Season,
  site = nodes_data$Site,
  virus_richness = rowSums(virus_matrix)
)

# Create network
network <- graph_from_data_frame(d = strong_edges, vertices = vertex_data, directed = FALSE)
isolated_nodes <- which(degree(network) == 0)
connected_network <- delete_vertices(network, isolated_nodes)

set.seed(42)
# Using Fruchterman-Reingold layout
L <- layout_with_fr(connected_network, niter = 5000)

layout_df <- data.frame(
  x = L[,1] * 30, 
  y = L[,2] * 30,
  name = V(connected_network)$name
) %>% left_join(vertex_data, by = "name")

# Prepare segments for ggplot
edge_df <- igraph::as_data_frame(connected_network, what = "edges") %>%
  inner_join(layout_df %>% dplyr::select(name, x, y), by = c("from" = "name")) %>%
  rename(from_x = x, from_y = y) %>%
  inner_join(layout_df %>% dplyr::select(name, x, y), by = c("to" = "name")) %>%
  rename(to_x = x, to_y = y)


# --- Step 5: Visualization ---
country_colors <- c("Cambodia"="#DB70AE", "Madagascar"="#F8766D", 
                    "Central African Republic"="#0B55D6", "Guyane"="#00C0AF")
season_colors <- c("Dry" = "#E69F00", "Rainy" = "#56B4E9")
year_shapes <- c("2019"=16, "2020"=17, "2021"=15, "2022"=18)

final_plot <- ggplot() +
  # 1. Edges (Thickened)
  geom_segment(data = edge_df, aes(x=from_x, y=from_y, xend=to_x, yend=to_y, 
                                   linewidth=weight, alpha=weight), color="grey50") +
  
  # 2. Nodes
  geom_point(data = layout_df, aes(x=x, y=y, size=virus_richness, color=country, shape=year), stroke=1.5) +
  
  # 3. STACKED LABELS: Species on top (italic), Site on bottom (bold)
  geom_text_repel(data = layout_df, 
                  aes(x=x, y=y, 
                      label=paste0(species, "\n", site), 
                      color=country),                    
                  size=8, 
                  fontface="bold.italic",                
                  lineheight = 0.9,                      
                  max.overlaps=Inf, 
                  box.padding=1.5,                       
                  point.padding=1.0, 
                  direction = "both",                    
                  bg.color = "white",                    
                  bg.r = 0.15) +                         
  
  scale_color_manual(values = country_colors) +
  
  scale_linewidth_continuous(range = c(2.0, 8.0), name = "Jaccard Similarity") +
  
  scale_size_continuous(range = c(8, 22), name = "Virus Richness") +
  
  scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
  theme_void() +
  theme(legend.position = "right", plot.margin = margin(20, 20, 20, 20))

ggsave("virome_network_CLEAN_LABELS.pdf", final_plot, width=40, height=40)








#### MADAGASCAR ####

# bipartite

# Filter for Madagascar only
mad_data <- sampled_data %>% 
  filter(Country == "Madagascar")

mad_data$Season <- factor(mad_data$Season, levels = c("Dry", "Rainy"))

library(dplyr)
library(ggraph)
library(igraph)
library(Cairo)

# 1. Create the Edge List
current_virus_cols <- intersect(keep_viruses, names(mad_data))
edges <- mad_data %>%
  pivot_longer(cols = all_of(current_virus_cols), names_to = "Virus", values_to = "Pres") %>%
  filter(Pres == 1) %>%
  dplyr::select(from = Mosq_Species_Name, to = Virus) %>%
  distinct()

# 2. Build Graph
graph <- graph_from_data_frame(d = edges, directed = FALSE)

# 3. Set Node Types and Degree
V(graph)$type <- V(graph)$name %in% mad_data$Mosq_Species_Name
V(graph)$degree <- degree(graph)

# 4. Define Status (Generalist = linked to 2+ genera)
V(graph)$status <- case_when(
  V(graph)$type ~ "Mosquito Host",
  V(graph)$degree >= 3 ~ "Virus: Generalist",
  TRUE ~ "Virus: Specialist"
)

# 5. FILTER LABELS: Only keep names for Hosts and Generalist Viruses
V(graph)$label_display <- ifelse(
  V(graph)$status %in% c("Mosquito Host", "Virus: Generalist"), 
  V(graph)$name, 
  ""
)

# 6. Plot with Leader Lines
p_final_network_padded <- ggraph(graph, layout = "bipartite") +
  
  geom_edge_diagonal(alpha = 0.3, color = "grey70", width = 0.6) +
  
  geom_node_point(aes(color = status, size = degree), alpha = 0.8) +
  
  geom_node_text(
    aes(label = label_display),
    repel = TRUE,
    size = 3.5,
    fontface = "bold",
    max.overlaps = Inf,
    segment.color = "grey30",
    segment.size = 0.5,
    segment.alpha = 0.6,
    min.segment.length = 0,
    box.padding = 0.8,
    point.padding = 0.5
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  
  scale_color_manual(values = c(
    "Mosquito Host" = "grey50", 
    "Virus: Generalist" = "royalblue3", 
    "Virus: Specialist" = "gold"
  )) +
  scale_size_continuous(range = c(2, 10), name = "Host Range (Degree)") +
  
  coord_fixed(ratio = 15) + 
  theme_graph() +
  theme(
    legend.position = "right",
    plot.margin = margin(30, 30, 30, 30) 
  ) +
  labs(
    title = "Madagascar Virus-Host Interactome",
    subtitle = "Fixed vertical clipping with Y-axis expansion"
  )

# 7. Save with increased height
ggsave("Madagascar_Bipartite_No_Cutoff.pdf", 
       plot = p_final_network_padded, 
       device = cairo_pdf, 
       width = 14, 
       height = 20) 





# seasonal
# Filter for Madagascar
mad_data <- sampled_data %>% 
  filter(Country == "Madagascar")

mad_data$alpha_richness <- rowSums(mad_data[, intersect(keep_viruses, names(mad_data))] > 0)

total_seasonal_richness <- mad_data %>%
  group_by(Site, Season) %>%
  summarise(
    total_richness = sum(colSums(across(all_of(intersect(keep_viruses, names(mad_data))))) > 0),
    .groups = "drop"
  )

p1 <- ggplot(total_seasonal_richness, aes(x = Season, y = total_richness, color = Site, group = Site)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Ambato Boeny" = "#D73027", "Amparafaravola" = "#FC8D59", "Mampikony" = "#FEE090")) +
  theme_minimal() +
  labs(title = "Total Virus Family Richness by Season", y = "Distinct Virus Families Detected")

# Save to PDF
ggsave("Madagascar_Seasonal_Richness.pdf", plot = p1, width = 8, height = 6)



# Venn
library(ggvenn)
library(tidyr)
library(dplyr)

# 1. Prepare list of viruses present per site
venn_list <- mad_data %>%
  pivot_longer(
    cols = all_of(intersect(keep_viruses, names(mad_data))), 
    names_to = "Virus", 
    values_to = "Presence"
  ) %>%
  filter(Presence > 0) %>%
  group_by(Site) %>%
  summarise(virus_list = list(unique(Virus)), .groups = "drop") %>%
  tibble::deframe()

# 2. Plot the Venn Diagram
p <- ggvenn(
  venn_list, 
  fill_color = c("#D73027", "#FC8D59", "#FEE090"),
  stroke_size = 0.5, 
  set_name_size = 4
)



# boxplot
mad_data$alpha_richness <- rowSums(mad_data[, intersect(keep_viruses, names(mad_data))] > 0)

# 2. Create the Boxplot
p_boxplot <- ggplot(mad_data, aes(x = Site, y = alpha_richness, fill = Site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = c("Ambato Boeny" = "#D73027", 
                               "Amparafaravola" = "#FC8D59", 
                               "Mampikony" = "#FEE090")) +
  theme_minimal() +
  labs(title = "Virus Family Richness by Site (Madagascar)",
       subtitle = "Calculated from filtered virus list (n >= 3 detections)",
       x = "Sampling Site",
       y = "Number of Virus Families per Sample") +
  theme(legend.position = "none")

# 3. View the plot
print(p_boxplot)



# heatmap
mad_long <- mad_data %>%
  filter(Country == "Madagascar") %>%
  pivot_longer(cols = all_of(intersect(keep_viruses, names(mad_data))), 
               names_to = "Virus", 
               values_to = "Presence") %>%
  filter(Presence > 0)

# 2. Identify Top 20 Viruses
top_20_mad_viruses <- mad_long %>%
  count(Virus) %>%
  slice_max(n, n = 20, with_ties = FALSE) %>%
  pull(Virus)

# 3. Create the Matrix
heatmap_data <- mad_long %>%
  filter(Virus %in% top_20_mad_viruses) %>%
  group_by(Site, Virus) %>%
  summarise(Presence = 1, .groups = "drop") %>%
  mutate(Site = factor(Site)) %>% 
  complete(Site, Virus = top_20_mad_viruses, fill = list(Presence = 0)) %>%
  mutate(fill_value = ifelse(Presence == 1, as.character(Site), "Absent"))

# 4. Plot Heatmap with Site-specific colors
p_heatmap <- ggplot(heatmap_data, aes(x = Site, y = Virus, fill = fill_value)) +
  geom_tile(color = "white", lwd = 0.5) +
  scale_fill_manual(
    values = c(
      "Ambato Boeny" = "#D73027", 
      "Amparafaravola" = "#FC8D59", 
      "Mampikony" = "#FEE090",
      "Absent" = "grey95"
    ),
    name = "Detection Status"
  ) +
  theme_minimal() +
  labs(title = "Presence/Absence of Top 20 Viruses (Madagascar Only)",
       x = "Site", y = "Virus Family") +
  theme(axis.text.y = element_text(face = "italic", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

print(p_heatmap)



ggsave("Madagascar_Venn_Diagram.pdf", p, width = 8, height = 8)
ggsave("Madagascar_Richness_Boxplot.pdf", p_boxplot, width = 8, height = 6)
ggsave("Madagascar_Virus_Heatmap.pdf", p_heatmap, width = 8, height = 10)









#### FIGURE STACKED BARPLOT (WHOLE DATASET) ####
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Load Data
raw_hits <- read.csv("virus_data.csv", header = TRUE)
virus_kingdom <- read.csv("updated_kingdom.csv", header = TRUE)

# 2. APPLY GLOBAL FILTER (n >= 3 detections across whole dataset)
keep_viruses_global <- raw_hits %>%
  group_by(Sample_Name, Virus_name_fin_cut) %>%
  summarise(presence = 1, .groups = "drop") %>%
  group_by(Virus_name_fin_cut) %>%
  summarise(total_detections = sum(presence)) %>%
  filter(total_detections >= 3) %>%
  pull(Virus_name_fin_cut)

data_filtered_all <- raw_hits %>%
  filter(Virus_name_fin_cut %in% keep_viruses_global)

# 4. Clean up and merge using 'family.supergroup'
data_kingdom_all <- data_filtered_all %>%
  mutate(family_clean = trimws(as.character(family.supergroup))) %>%
  mutate(family_clean = ifelse(family_clean == "" | is.na(family_clean), "unclassified", family_clean)) %>%
  left_join(virus_kingdom, by = c("family_clean" = "fam")) %>%
  mutate(kingdom = ifelse(is.na(kingdom), "unassigned", kingdom))

# 5. Global Color Scheme
color_scheme <- c("dsRNA" = "dodgerblue1", 
                  "negRNA" = "burlywood3", 
                  "posRNA" = "olivedrab", 
                  "RT" = "orchid4", 
                  "ssDNA" = "lightblue2", 
                  "unassigned" = "grey")

create_global_barplot <- function(df, factor_name) {
  proportions <- df %>%
    group_by(!!sym(factor_name), kingdom) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(!!sym(factor_name)) %>%
    mutate(total = sum(count),
           proportion = count / total)
  
  ggplot(proportions, aes(x = !!sym(factor_name), y = proportion, fill = kingdom)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    scale_fill_manual(values = color_scheme) +
    theme_minimal() +
    labs(title = paste("Global Viral Kingdom Distribution:", factor_name),
         subtitle = paste("Viruses with n >= 3 detections (Total:", length(keep_viruses_global), ")"),
         x = factor_name, y = "Proportion of Contig Hits") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
}

# 7. Generate Plots for the Whole Dataset
country_plot  <- create_global_barplot(data_kingdom_all, "Country")
species_plot  <- create_global_barplot(data_kingdom_all, "Mosq_Species_Name")
site_all_plot <- create_global_barplot(data_kingdom_all, "Site")
year_all_plot <- create_global_barplot(data_kingdom_all, "Collection_Y")

# 8. Save all Global Plots
pdf("Global_Filtered_Viral_Kingdoms.pdf", width = 12, height = 8)
print(country_plot)
print(species_plot)
print(site_all_plot)
print(year_all_plot)
dev.off()

# Print summary to console
print(paste("Global analysis includes", length(keep_viruses_global), "validated virus species."))











#### HEATMAP OF SAMPLING DATA ####
library(ggplot2)
library(dplyr)
library(tidyr)

heatmap_data <- data %>%
  group_by(Country, Site, Collection_Y, Season) %>%
  summarise(n_samples = n_distinct(Sample_Name), .groups = "drop") %>%
  
  complete(nesting(Country, Site), 
           Collection_Y = unique(data$Collection_Y), 
           Season = c("Rainy", "Dry"), 
           fill = list(n_samples = 0)) %>%
  
  filter(!(Collection_Y == 2021 & Season == "Dry"))

heatmap_data$Season <- factor(heatmap_data$Season, levels = c("Rainy", "Dry"))
heatmap_data$Collection_Y <- factor(heatmap_data$Collection_Y)

heatmap_data <- heatmap_data %>%
  arrange(Country, desc(Site))
heatmap_data$Site <- factor(heatmap_data$Site, levels = unique(heatmap_data$Site))

# 5. Create the Heatmap
p_effort <- ggplot(heatmap_data, aes(x = Season, y = Site, fill = n_samples)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n_samples), size = 3.5, fontface = "bold") +
  
  scale_fill_gradient(low = "grey95", high = "#E41A1C", name = "Samples") +
  
  facet_grid(Country ~ Collection_Y, scales = "free_y", space = "free_y", switch = "y") +
  
  theme_minimal() +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "white", color = "grey80"),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
    strip.text.x = element_text(face = "bold", size = 12),
    
    panel.spacing = unit(0.2, "lines"),
    panel.grid = element_blank(),
    
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) +
  labs(title = "Sampling Effort: Number of Unique Mosquito Samples",
       subtitle = "Rows: Country & Site | Columns: Year & Season")

# 6. Save to PDF
ggsave("Sampling_Effort_Heatmap.pdf", p_effort, width = 11, height = 8.5)

print(p_effort)






#### UNCLASSIFIEDS ####
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

data <- read.csv("virus_data.csv", na.strings = c("", "NA", "nan", " "), stringsAsFactors = FALSE)

# This counts how many contigs per species was already curated
classified_counts <- data %>%
  group_by(Mosq_Species_Name) %>%
  summarise(n_classified = n(), .groups = "drop")

# 2. Identify and Map the Novel contigs
class_info <- read.table("unclassifieds_contig_list", header = TRUE, sep = "\t", fill = TRUE)
contig_map <- read.table("fasta_all_contigs_name_map.txt", header = FALSE, sep = "\t")
colnames(contig_map) <- c("Library_File", "NODE_ID")

raw_mapped_data <- class_info %>%
  inner_join(contig_map, by = "NODE_ID") %>%
  mutate(Sample_Name = sub("_scaffolds_1kb.fasta", "", Library_File))

# 3. Filter for "Novel" and "Unclassified" (The Numerator)
unclassified_novel_hits <- raw_mapped_data %>%
  filter(!(NODE_ID %in% data$NODE_ID)) %>%
  filter(family == "unclassified") %>%
  left_join(distinct(data, Sample_Name, Mosq_Species_Name), by = "Sample_Name") %>%
  group_by(Mosq_Species_Name) %>%
  summarise(n_unclassified = n(), .groups = "drop")


proportion_data <- classified_counts %>%
  left_join(unclassified_novel_hits, by = "Mosq_Species_Name") %>%
  mutate(n_unclassified = replace_na(n_unclassified, 0)) %>%
  mutate(proportion = (n_unclassified / (n_classified + n_unclassified)) * 100)

overall_avg <- mean(proportion_data$proportion, na.rm = TRUE)

p_novelty <- ggplot(proportion_data, aes(x = Mosq_Species_Name, y = proportion)) +
  geom_bar(stat = "identity", fill = "grey60", width = 0.7) +
  geom_hline(yintercept = overall_avg, linetype = "dashed", color = "red3", linewidth = 0.8) +
  annotate("text", x = nrow(proportion_data), y = overall_avg + 1.5, 
           label = paste0("Average: ", round(overall_avg, 1), "%"), 
           color = "red3", fontface = "bold", hjust = 1) +
  theme_minimal() +
  labs(
    title = "Viral Discovery Ratio by Mosquito Species",
    subtitle = "Proportion of Unclassified 'Leftover' contigs relative to Curated Hits in CSV",
    x = "Mosquito species",
    y = "Proportion (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major.x = element_blank()
  )

print(p_novelty)
ggsave("Barplot_unclassifieds.pdf", p_novelty, width = 11, height = 8.5)





# 1. LOAD DATA
viral_hits <- read.csv("unclass_identity.csv")

species_map <- read_delim("species_unclass", 
                          delim = "\t", 
                          trim_ws = TRUE)


# 2. JOIN AND CLEAN
plot_data <- viral_hits %>%
  left_join(species_map, by = "Sample_Name") %>%
  filter(!is.na(Mosq_Species_Name)) %>%
  mutate(pident = as.numeric(pident))

# 3. CREATE THE PLOT
identity_plot <- ggplot(plot_data, aes(x = Mosq_Species_Name, y = pident)) +
  geom_boxplot(fill = "#2196F3", color = "black", outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue", size = 1.5) +
  theme_minimal() +
  labs(
    title = "Amino Acid Identity of Unclassified Viruses",
    x = "Mosquito Species",
    y = "Amino Acid Identity (%)"
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.minor = element_blank()
  )

# 4. VIEW AND SAVE
print(identity_plot)
ggsave("Unclassified_Identity_by_Species.pdf", width = 10, height = 6)







#### HEATMAP TOP 50 RPM BY SPECIES ####
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# 1. LOAD DATA
analysis_data <- read.csv("virus_data.csv", stringsAsFactors = FALSE)

# 2. SELECT TOP 50 VIRUSES GLOBALLY
top_50_species <- analysis_data %>%
  group_by(Virus_name_fin_cut) %>%
  summarise(Global_Total = sum(RPM, na.rm = TRUE)) %>%
  slice_max(Global_Total, n = 50) %>%
  pull(Virus_name_fin_cut)

# 1. PREPARE THE DATA
heatmap_prep <- analysis_data %>%
  filter(Virus_name_fin_cut %in% top_50_species) %>%
  group_by(Virus_name_fin_cut, Mosq_Species_Name) %>%
  summarise(mean_rpm = mean(RPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(log_rpm = log10(mean_rpm + 1))

# 2. TRANSFORM TO WIDE MATRIX
matrix_data <- heatmap_prep %>%
  dplyr::select(Virus_name_fin_cut, Mosq_Species_Name, log_rpm) %>%
  pivot_wider(names_from = Mosq_Species_Name, values_from = log_rpm, values_fill = 0)

mat <- as.matrix(matrix_data[,-1])
rownames(mat) <- matrix_data$Virus_name_fin_cut

# 3. DEFINE THE COLOR SCALE
max_val <- max(mat, na.rm = TRUE)
col_fun = colorRamp2(c(0, 0.01, max_val), c("lightgrey", "#f7fbff", "#084594"))

# 4. DRAW AND SAVE
pdf("Virus_Host_Fixed_Heatmap.pdf", width = 10, height = 12)

ht <- Heatmap(mat, 
              name = "Log10 Mean RPM", 
              col = col_fun,
              cluster_rows = TRUE, 
              cluster_columns = TRUE,
              show_row_dend = TRUE, 
              show_column_dend = TRUE,
              
              # Sizing for rectangular cells
              width = unit(8, "cm"), 
              height = unit(16, "cm"),
              
              row_names_gp = gpar(fontsize = 7, fontface = "italic"),
              column_names_gp = gpar(fontsize = 8),
              
              # Thinner borders
              rect_gp = gpar(col = "white", lwd = 0.3), 
              
              column_title = "Mosquito Host Species",
              row_title = "Top 50 Viruses",
              
              heatmap_legend_param = list(title_gp = gpar(fontsize = 9))
)

draw(ht)
dev.off()





#### HEATMAP TOP 50 RPM BY SITE ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# 1. LOAD DATA
analysis_data <- read.csv("virus_data.csv", stringsAsFactors = FALSE)

# 2. SELECT TOP 50 VIRUSES GLOBALLY (Keep the same top viruses for consistency)
top_50_species <- analysis_data %>%
  group_by(Virus_name_fin_cut) %>%
  summarise(Global_Total = sum(RPM, na.rm = TRUE)) %>%
  slice_max(Global_Total, n = 50) %>%
  pull(Virus_name_fin_cut)

# 3. PREPARE THE DATA BY SITE
heatmap_prep <- analysis_data %>%
  filter(Virus_name_fin_cut %in% top_50_species) %>%
  # Group by Site instead of Species
  group_by(Virus_name_fin_cut, Site, Country) %>% 
  summarise(mean_rpm = mean(RPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(log_rpm = log10(mean_rpm + 1))

# 4. TRANSFORM TO WIDE MATRIX
matrix_data <- heatmap_prep %>%
  dplyr::select(Virus_name_fin_cut, Site, log_rpm) %>%
  pivot_wider(names_from = Site, values_from = log_rpm, values_fill = 0)

mat <- as.matrix(matrix_data[,-1])
rownames(mat) <- matrix_data$Virus_name_fin_cut

# 5. CREATE SITE METADATA (For grouping columns by Country)
site_info <- heatmap_prep %>% 
  dplyr::select(Site, Country) %>% 
  distinct() %>%
  arrange(Country) 

mat <- mat[, site_info$Site]

# 6. DEFINE THE COLOR SCALE
max_val <- max(mat, na.rm = TRUE)
col_fun = colorRamp2(c(0, 0.01, max_val), c("lightgrey", "#f7fbff", "#084594"))

# 7. DRAW AND SAVE
pdf("Virus_Site_Heatmap.pdf", width = 12, height = 12)

ht <- Heatmap(mat, 
              name = "Log10 Mean RPM", 
              col = col_fun,
              cluster_rows = TRUE, 
              cluster_columns = TRUE,
              
              column_split = site_info$Country,
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              
              show_row_dend = TRUE, 
              show_column_dend = TRUE,
              
              width = unit(12, "cm"), 
              height = unit(16, "cm"),
              
              row_names_gp = gpar(fontsize = 7, fontface = "italic"),
              column_names_gp = gpar(fontsize = 8),
              
              rect_gp = gpar(col = "white", lwd = 0.3), 
              
              row_title = "Top 50 Viruses",
              column_title = "Sampling Sites (Grouped by Country)",
              
              heatmap_legend_param = list(title_gp = gpar(fontsize = 9))
)

draw(ht)
dev.off()


#####################################