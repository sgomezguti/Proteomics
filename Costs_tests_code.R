# set a working directory - this is on my desktop PC
# This script is connected to run citT_relative_abundance.R script
rm(list=ls()) #clean the space

# load packages
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)
library(SummarizedExperiment)
library("DEP")
library(ggplot2)
library(ggrepel)

# load data with protein intensities
# Load my data
LTEE <- read.csv(file="Data_for_R.csv", header = TRUE) # do not use row.names=1 because the code below wont work well
LTEE

# Include localization information in the dataframe
# Load localization results from DeepLocPro1.0
localization <- read.csv(file="./Localization_prediction/results_20260619-202210.csv", header = TRUE)
localization

# Find corresponding proteinID in localization to include the prediction in prot_totals_localization
LTEE_localization <- merge(x = LTEE, y = localization[, c("ProteinID", "Localization")], by = "ProteinID", all.x=TRUE) 

# ── Identify replicate columns automatically ───────────────────────────────
# Matches any column ending in _R1, _R2, or _R3 (e.g., ZDB143_R1)
strain_cols <- grep("_R[123]$", colnames(LTEE_localization), value = TRUE)
strain_cols   # check this looks right - should list all your replicate columns


# ── Step 1: Filter to cytoplasmic membrane proteins only ──────────────────
cytmembrane <- LTEE_localization %>%
  filter(Localization == "Cytoplasmic Membrane")

# ── Step 2: Calculate total cytoplasmic membrane intensity per sample (replicate) ──
cytmembrane_totals <- cytmembrane %>%
  summarise(across(all_of(strain_cols), ~ sum(., na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "cytmembrane_total")

#################################################################################
#Additional here. Plot total proteome mass allocated to the cytoplasmic membrane
total_mass_cytproteome_long <-cytmembrane_totals %>%
  mutate(
    clone     = str_remove(sample, "_R[123]$"),
    replicate = str_extract(sample, "R[123]$")
  )
total_mass_cytproteome_long <- merge(total_mass_cytproteome_long, clone_meta, by = "clone")
#Plot
#Only for Cit+ clones
total_mass_cytproteome_long_citplus <- total_mass_cytproteome_long %>% filter(clone %in% citplus_clones)

ggplot(total_mass_cytproteome_long_citplus, aes(x = generation, y = cytmembrane_total, color = clone)) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05), show.legend = FALSE) +
  geom_text_repel(data = total_mass_cytproteome_long_citplus %>% filter(replicate == "R1"),
                  aes(label = clone), size = 7, show.legend = FALSE) +
  # scale_x_continuous(breaks = unique(citT_long_citplus$generation)) +
  scale_x_continuous(breaks = seq(0, max(total_mass_cytproteome_long_citplus$generation), by = 2000)) + # ticks even every 2000 generations
  scale_y_continuous(breaks = seq(0, max(total_mass_cytproteome_long_citplus$cytmembrane_total))) +
  scale_color_manual(values = citi_palette) +   # palette for more than 8 clones
  labs(x = "Generation", 
       y = "Total proteome mass of cytoplasmic membrane)",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 20, color = "black")
  )


################################################################################
# ── Step 3: Pull out the citT row (still within cytoplasmic membrane subset) ──
citT_row <- cytmembrane %>%
  filter(Genes == "citT")   # adjust to match your actual citT identifier

citT_abundance <- citT_row %>%
  pivot_longer(cols = all_of(strain_cols), names_to = "sample", values_to = "citT_intensity")

# 3.1. For dctA to see if htere is any correlation
dctA_row <- cytmembrane %>%
  filter(Genes == "dctA")

dctA_abundance <- dctA_row %>%
  pivot_longer(cols = all_of(strain_cols), names_to = "sample", values_to = "dctA_intensity")

# ── Step 4: Calculate citT relative abundance within the cytoplasmic membrane proteome ──
citT_long <- citT_abundance %>%
  left_join(cytmembrane_totals, by = "sample") %>%
  mutate(
    abundance = citT_intensity / cytmembrane_total,
    clone     = str_remove(sample, "_R[123]$"),
    replicate = str_extract(sample, "R[123]$")
  )

dctA_long <- dctA_abundance %>%
  left_join(cytmembrane_totals, by = "sample") %>%
  mutate(
    abundance = dctA_intensity / cytmembrane_total,
    clone = str_remove(sample, "_R[123]$"),
    replicate = str_extract(sample, "R[123]$")
  )

# ── Step 5: Bring in module counts per clone ───────────────────────────────
# Metadata table: one row per clone, with module count and generation
clone_meta <- data.frame(
  clone = c("ZDB564", "ZDB172", "ZDB143","CZB152", "CZB154", "ZDB96",
            "ZDB107", "REL10979"),
  modules    = c(1, 4, 3, 8, 4, 3, 3, 3),  # your actual counts
  generation = c(31500, 32000, 32500, 33000, 33000,
                 36000, 38000, 40000)
)

# ALternative clone_meta with all clones
clone_meta <- data.frame(
  clone = c("REL606", "ZDB199", "ZDB30", "CZB199", "ZDB99",
            "ZDB111", "REL10988", "CZB152", "CZB154", "REL10979",
            "ZDB107", "ZDB143", "ZDB172", "ZDB564", "ZDB96"),
  modules    = c(0, 0, 0, 0, 0, 0, 0, 8, 4, 3, 3, 3, 4, 1, 3),  # your actual counts
  generation = c(0, 31500, 32000, 33000, 36000, 38000, 40000,
                 33000, 33000, 40000, 38000, 32500, 32000, 31500,
                 36000)
)

citT_long <- merge(citT_long, clone_meta, by = "clone")

citT_long_efficiency <- citT_long %>%
  mutate(efficiency = abundance / modules)

dctA_long <- merge(dctA_long, clone_meta, by = "clone")

# ── Step 6: Plot with all three replicates shown ───────────────────────────
# ── Plot 1: Individual replicate points, colored by clone ─────────────────
# Modules on x axis
library(RColorBrewer)
citi_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(15)

# Filter citT_long to show only Cit+ clones
citplus_clones <- c("CZB152", "CZB154", "REL10979", "ZDB107", "ZDB143", "ZDB172", "ZDB564", "ZDB96")

citT_long_citplus <- citT_long %>% filter(clone %in% citplus_clones)

ggplot(citT_long_citplus, aes(x = modules, y = abundance, color = clone)) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05)) +
  geom_text_repel(data = citT_long_citplus %>% filter(replicate == "R1"),
                  aes(label = clone), size = 3.5, show.legend = FALSE) +
  scale_x_continuous(breaks = unique(citT_long_citplus$modules)) +
  scale_color_manual(values = citi_palette) +   # palette for more than 8 clones
  labs(x = "Number of rnk-citT modules", 
       y = "   Relative citT abundance \
   (cyt membrane proteome)",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black")
  )

# Generation on x axis
ggplot(citT_long_citplus, aes(x = generation, y = abundance, color = clone)) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05), show.legend = FALSE) +
  geom_text_repel(data = citT_long_citplus %>% filter(replicate == "R1"),
                  aes(label = clone), size = 7, show.legend = FALSE) +
 # scale_x_continuous(breaks = unique(citT_long_citplus$generation)) +
  scale_x_continuous(breaks = seq(0, max(citT_long_citplus$generation), by = 2000)) + # ticks even every 2000 generations
  scale_color_manual(values = citi_palette) +   # palette for more than 8 clones
  labs(x = "Generation", 
       y = "citT expression (proteome mass fraction \n of cytoplasmic membrane)",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 20, color = "black")
  )

# Cladude improvements to add refinement and actualization
ggplot(citT_long_citplus, aes(x = generation, y = abundance, color = clone)) +
  # Shaded band highlighting the actualization → refinement transition window
  annotate("rect", xmin = 31500, xmax = 33000, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.4) +
  
  # Vertical reference lines
  geom_vline(xintercept = 31500, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = 33000, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  
  # Data layers (unchanged)
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05), show.legend = FALSE) +
  geom_text_repel(data = citT_long_citplus %>% filter(replicate == "R1"),
                  aes(label = clone), size = 7, show.legend = FALSE) +
  
  # Annotation text, rotated and placed near the top of the plot
  annotate("text", x = 31500, y = max(citT_long_citplus$abundance) * 1.05,
           label = "Actualization", angle = 90, vjust = -0.5, hjust = 1,
           size = 5, fontface = "italic", color = "grey20") +
  annotate("text", x = 33000, y = max(citT_long_citplus$abundance) * 1.05,
           label = "Refinement", angle = 90, vjust = -0.5, hjust = 1,
           size = 5, fontface = "italic", color = "grey20") +
  
  scale_x_continuous(breaks = seq(0, max(citT_long_citplus$generation), by = 2000)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +  # extra headroom for the text labels
  scale_color_manual(values = citi_palette) +
  labs(x = "Generation", 
       y = "citT expression (proteome mass fraction \n of cytoplasmic membrane)",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 20, color = "black")
  )
# Second try
library(patchwork)
max_gen <- max(citT_long_citplus$generation)

# ── Main plot: vertical "Actualization" label next to the 31500 line ──────
main_plot <- ggplot(citT_long_citplus, aes(x = generation, y = abundance, color = clone)) +
  annotate("rect", xmin = 33000, xmax = 40000, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.4) +
  geom_vline(xintercept = 31500, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = 33000, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05), show.legend = FALSE) +
  geom_text_repel(data = citT_long_citplus %>% filter(replicate == "R1"),
                  aes(label = clone), size = 7, show.legend = FALSE) +
  annotate("text", x = 31500, y = max(citT_long_citplus$abundance) * 1.05,
           label = "Actualization", angle = 90, vjust = -0.5, hjust = 0.7,
           size = 7, fontface = "italic", color = "grey20") +
  scale_x_continuous(breaks = seq(32000, 40000, by = 2000),
                     limits = c(31000, 41000), expand = c(0.02, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_color_manual(values = citi_palette) +
  labs(x = NULL,
       y = "citT expression (proteome mass fraction \n of cytoplasmic membrane)",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 20, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 2)
  )

# ── Annotation strip: only "Refinement" now ────────────────────────────────
annotation_data <- tibble(
  xmin  = 33000,
  xmax  = max_gen,
  stage = "Refinement"
)

annotation_strip <- ggplot(annotation_data) +
  annotate("rect", xmin = 33000, xmax = 40000, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.4) +
  geom_vline(xintercept = 31500, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = 33000, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_text(data = annotation_data,
            aes(x = (xmin + xmax) / 2, y = 0, label = stage),
            size = 7, fontface = "italic", color = "grey20") +
  scale_x_continuous(breaks = seq(32000, 40000, by = 2000),
                     limits = c(31000, 41000), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Generation", y = NULL) +
  theme_void(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 8), vjust = -1),
    axis.text.x  = element_text(size = 20, color = "black", vjust = -1),
    axis.ticks.x = element_line(color = "black"),
    axis.line.x  = element_line(color = "black"),
    plot.margin  = margin(t = 2)
  )

# ── Stack them ───────────────────────────────────────────────────────────
main_plot / annotation_strip + plot_layout(heights = c(9, 1))


# ----------------------------
# Save to file
# ----------------------------
ggsave(
  filename = "./Costs_tests_poster/citTexpression_updated.svg",
  width    = 10,
  height   = 7,
  dpi      = 300,
  device   = "svg"
)
# Plot 1 but for dctA
# Step below maybe not good idea becaause we should leave generation as quantiative and not qualitative factor
#dctA_long <- dctA_long %>% mutate(generation_factor = factor(generation, levels = sort(unique(generation))))

library(RColorBrewer)
dctA_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(15)
# for all clones
ggplot(dctA_long, aes(x = generation, y = abundance, color = clone)) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05)) +
  geom_text_repel(data = dctA_long %>% filter(replicate == "R1"),
                  aes(label = clone), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = dctA_palette) +   # 8 distinguishable colors; swap palette if you have >8 clones
  labs(x = "Generation", 
       y = "Relative dctA abundance\n(cytoplasmic membrane proteome)",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black")
  )

# for Cit+ clones only
dctA_long_citplus <- dctA_long %>% filter(clone %in% citplus_clones)

ggplot(dctA_long_citplus, aes(x = generation, y = abundance, color = clone)) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05), show.legend = FALSE) +
  geom_text_repel(data = dctA_long_citplus %>% filter(replicate == "R1"),
                  aes(label = clone), size = 7, show.legend = FALSE) +
  scale_color_manual(values = dctA_palette) +   # 8 distinguishable colors; swap palette if you have >8 clones
  labs(x = "Generation", 
       y = "Relative dctA abundance\n(cytoplasmic membrane proteome)",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 20, color = "black")
  )

# Claude version, same improvements than for citT

# ── Main plot: vertical "Actualization" label next to the 31500 line ──────
main_plot_dctA <- ggplot(dctA_long_citplus, aes(x = generation, y = abundance, color = clone)) +
  annotate("rect", xmin = 33000, xmax = 40000, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.4) +
  geom_vline(xintercept = 31500, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = 33000, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05), show.legend = FALSE) +
  geom_text_repel(data = dctA_long_citplus %>% filter(replicate == "R1"),
                  aes(label = clone), size = 7, show.legend = FALSE) +
  annotate("text", x = 31500, y = max(dctA_long_citplus$abundance) * 1.05,
           label = "Actualization", angle = 90, vjust = -0.5, hjust = 0.7,
           size = 7, fontface = "italic", color = "grey20") +
  scale_x_continuous(breaks = seq(32000, 40000, by = 2000),
                     limits = c(31000, 41000), expand = c(0.02, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_color_manual(values = citi_palette) +
  labs(x = NULL,
       y = "dctA expression (proteome mass fraction \n of cytoplasmic membrane)",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 20, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 2)
  )
main_plot_dctA
# ── Annotation strip: only "Refinement" now ────────────────────────────────
annotation_data <- tibble(
  xmin  = 33000,
  xmax  = max_gen,
  stage = "Refinement"
)

annotation_strip <- ggplot(annotation_data) +
  annotate("rect", xmin = 33000, xmax = 40000, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.4) +
  geom_vline(xintercept = 31500, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = 33000, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  geom_text(data = annotation_data,
            aes(x = (xmin + xmax) / 2, y = 0, label = stage),
            size = 7, fontface = "italic", color = "grey20") +
  scale_x_continuous(breaks = seq(32000, 40000, by = 2000),
                     limits = c(31000, 41000), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Generation", y = NULL) +
  theme_void(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 8), vjust = -1),
    axis.text.x  = element_text(size = 20, color = "black", vjust = -1),
    axis.ticks.x = element_line(color = "black"),
    axis.line.x  = element_line(color = "black"),
    plot.margin  = margin(t = 2)
  )

# ── Stack them ───────────────────────────────────────────────────────────
main_plot_dctA / annotation_strip + plot_layout(heights = c(9, 1))


# ----------------------------
# Save to file
# ----------------------------
ggsave(
  filename = "./Costs_tests_poster/Succinatetransporter_generation.svg",
  width    = 10,
  height   = 7,
  dpi      = 300,
  device   = "svg"
)

#------------------------------------------------------------------------------#
# TRy sergio's idea of plotting numer of modules on the right side
# Get the ranges for rescaling
ab_range  <- range(citT_long_citplus$abundance, na.rm = TRUE)
mod_range <- range(citT_long_citplus$modules, na.rm = TRUE)

# Function to rescale modules onto the abundance scale (for plotting)
rescale_mod_to_ab <- function(mod) {
  (mod - mod_range[1]) / (mod_range[2] - mod_range[1]) * diff(ab_range) + ab_range[1]
}
# Inverse function (for the secondary axis labels — converts back to true module values)
rescale_ab_to_mod <- function(ab) {
  (ab - ab_range[1]) / diff(ab_range) * diff(mod_range) + mod_range[1]
}

# One row per clone/generation for the module markers (avoid overplotting across replicates)
module_points <- citT_long_citplus %>%
  distinct(clone, generation, modules) %>%
  mutate(modules_scaled = rescale_mod_to_ab(modules))

ggplot(citT_long_citplus, aes(x = generation, y = abundance, color = clone)) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05)) +
  geom_text_repel(data = citT_long_citplus %>% filter(replicate == "R1"),
                  aes(label = clone), size = 3.5, show.legend = FALSE) +
  geom_point(data = module_points, 
             aes(x = generation, y = modules_scaled), 
             inherit.aes = FALSE, shape = 8, size = 3, color = "black") +
  scale_x_continuous(breaks = seq(0, max(citT_long_citplus$generation), by = 2000)) +
  scale_y_continuous(
    name = "   Relative citT abundance \n   (cyt membrane proteome)",
    sec.axis = sec_axis(~ rescale_ab_to_mod(.), name = "Number of rnk-citT modules")
  ) +
  scale_color_manual(values = citi_palette) +
  labs(x = "Generation", color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text  = element_text(size = 14, color = "black")
  )
#------------------------------------------------------------------------------#

# Alternative plot, plot efficiency as citT rel abundance divided by modules
citT_long_efficiency <- citT_long_efficiency %>%
  mutate(generation_factor = factor(generation, levels = sort(unique(generation))))

citT_long_efficiency <- citT_long_efficiency %>%
  mutate(modules_factor = factor(modules, levels = sort(unique(modules))))

ggplot(citT_long_efficiency, aes(x = modules_factor, y = efficiency, color = clone)) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.05)) +
  geom_text_repel(data = citT_long_efficiency %>% filter(replicate == "R1"),
                  aes(label = clone), size = 3.5, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Number of Modules", 
       y = "citT abundance/N. modules",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black")
  )

library(lme4)
model <- lmer(efficiency ~ modules + (1|clone), data = citT_long_efficiency)
summary(model)

# ── Plot 2: Mean ± SD (or SEM) per clone, no individual points ────────────
citT_summary <- citT_long %>%
  group_by(clone, modules) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    sd_abundance   = sd(abundance, na.rm = TRUE),
    n              = n(),
    sem_abundance  = sd_abundance / sqrt(n),
    .groups = "drop"
  )

ggplot(citT_summary, aes(x = modules, y = mean_abundance, color = clone)) +
  geom_point(size = 3.5, position = position_jitter(width = 0.05)) +
  geom_errorbar(aes(ymin = mean_abundance - sd_abundance, 
                    ymax = mean_abundance + sd_abundance),
                width = 0.08, position = position_jitter(width = 0.05)) +
  geom_text_repel(aes(label = clone), size = 3.5, show.legend = FALSE) +
  scale_x_continuous(breaks = unique(citT_summary$modules)) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Number of rnk-citT modules", 
       y = "Relative CitT abundance (mean ± SD, within cytoplasmic membrane proteome)",
       color = "Clone") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black")
  )



# ----------------------------
# Save to file
# ----------------------------
ggsave(
  filename = "./Costs_tests_poster/efficiencypercopy_modules.svg",
  width    = 10,
  height   = 6,
  dpi      = 300,
  device   = "svg"
)
